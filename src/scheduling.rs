//! Scheduling: run one input through a per-record visitor, on one thread or on many.
//!
//! Never parses — it picks what each worker reads and hands those bytes to
//! [`crate::files::RandomAccess`] or [`crate::files::ForwardAccess`]. Three access tiers at
//! two levels of parallelism give six functions, and [`for_each_record_in_input`] is what a
//! tool's `main` calls to reach the right one. Deduplication schedules its own passes and
//! depends on nothing here.

use std::io::{self, Read};
use std::ops::Range;
use std::sync::atomic::{AtomicBool, Ordering};

use forkunion::{DynamicScheduler, IntoParallelIterator, ThreadPool, Topology};

#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
use crate::blocks::BlockAccess;
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
use crate::files::StreamStart;
use crate::files::{ForwardAccess, InputFile, RandomAccess};
use crate::records::{
    last_record_boundary, next_record_boundary, sequence_format, Record, SequenceFormat,
};

// region: Granularity

/// Decoded bytes one worker aims to take at a time.
///
/// A constant rather than a division of the input by the worker count: a 200 GB file split
/// per worker gives 25 GB extents, where one straggler stalls the join for a full worker's
/// time and the buffered output scales with the file rather than with a budget.
pub const TARGET_WORK_BYTES: usize = 4 << 20;

/// States a tool builds per thread when its output must mirror its input.
///
/// Several, so a slow extent is stolen rather than idling the pool, and few, because each one
/// buffers a whole extent's output until its turn to be written comes round.
const STATES_PER_THREAD: usize = 4;

/// Parse window for one run of blocks, which never has to hold more than one record.
///
/// The [`crate::files::DEFAULT_WINDOW_BYTES`] window would cost a megabyte per state, and
/// sixty-four states is then sixty-four megabytes of window for records that are kilobytes.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
const BLOCK_WINDOW_BYTES: usize = 256 << 10;

/// Whether a tool's output has to mirror the order of its input.
///
/// The pool indexes states differently for the two answers, which is worth a genuine amount
/// of memory, so a tool states which it needs rather than paying for the stricter one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordOrder {
    /// Output mirrors the input, so one state takes one extent and states retire in turn.
    Preserved,
    /// Nothing downstream depends on which state saw what, so states are per thread and the
    /// tool merges them once the run returns.
    Ignored,
}

/// How many worker threads a tool may use.
#[derive(clap::Args, Debug, Clone, Copy)]
pub struct Parallelism {
    /// Worker threads; 1 runs serially, 0 uses one per logical core
    #[arg(short = 'j', long = "jobs", default_value_t = 1)]
    pub jobs: usize,
}

impl Parallelism {
    /// Workers for a run whose output must mirror its input, record for record.
    ///
    /// Every tool that writes records wants this. Each unit of work fills its own buffer and
    /// the buffers are written out in order, so the bytes are the same however many ran.
    pub fn ordered(self) -> io::Result<Workers> {
        self.spawn(RecordOrder::Preserved)
    }

    /// Workers for a run where nothing downstream depends on which state saw what.
    ///
    /// Summaries merge, so one accumulator per thread suffices and none of them wait on
    /// another. `fastq-stats` is the whole category.
    pub fn unordered(self) -> io::Result<Workers> {
        self.spawn(RecordOrder::Ignored)
    }

    /// Serial by default, so behaviour is unchanged until someone asks, and a workflow engine
    /// already running one job per core is not oversubscribed by every tool it invokes.
    fn spawn(self, order: RecordOrder) -> io::Result<Workers> {
        let topology = Topology::new()
            .map_err(|error| io::Error::other(format!("cannot read the machine: {error:?}")))?;
        let threads = match self.jobs {
            0 => topology.logical_cores_count().max(1),
            asked => asked,
        };
        if threads <= 1 {
            return Ok(Workers {
                pool: None,
                topology,
                order,
            });
        }
        let pool = ThreadPool::try_spawn(&topology, threads).map_err(|error| {
            io::Error::other(format!("cannot spawn {threads} worker threads: {error:?}"))
        })?;
        Ok(Workers {
            pool: Some(pool),
            topology,
            order,
        })
    }
}

/// How a run is spread across threads, and how its accumulators are shaped.
///
/// A tool builds one of these from `-j` and hands it to the driver. It never names a thread
/// count, a pool, or an ordering again: `states` sizes the accumulators to match, and one
/// worker running serially is just the case where there is no pool.
pub struct Workers {
    /// Declared first, so the threads are joined before the view that placed them is dropped.
    pool: Option<ThreadPool>,
    /// Held for the pool's lifetime, because the pool reads it while it runs.
    #[allow(dead_code)]
    topology: Topology,
    order: RecordOrder,
}

impl Workers {
    /// One accumulator per unit of work in flight, each built by `make`.
    ///
    /// The count follows from the ordering and the thread count, which is why a tool asks for
    /// the accumulators rather than deriving how many it needs.
    pub fn states<T>(&self, make: impl FnMut() -> T) -> Vec<T> {
        let count = match (&self.pool, self.order) {
            (None, _) => 1,
            (Some(pool), RecordOrder::Preserved) => pool.threads_count() * STATES_PER_THREAD,
            (Some(pool), RecordOrder::Ignored) => pool.threads_count(),
        };
        std::iter::repeat_with(make).take(count).collect()
    }

    /// The pool and the ordering, for a driver about to dispatch.
    fn dispatch(&mut self) -> (Option<&mut ThreadPool>, RecordOrder) {
        (self.pool.as_mut(), self.order)
    }
}

// endregion: Granularity

// region: Serial

/// Visit every record of a buffer already in memory, folding into `state`.
pub fn for_each_record_in_bytes<T>(
    bytes: &[u8],
    state: &mut T,
    mut visit: impl FnMut(&mut T, Record<'_>) -> io::Result<()>,
) -> io::Result<()> {
    RandomAccess::new(bytes).for_each_record(|record| visit(state, record))
}

/// Visit every record of a block-addressable input, folding into `state`.
///
/// Reads forward from block zero, so nothing here asks for the block table and nothing walks
/// the block chain to build one.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
pub fn for_each_record_in_blocks<T>(
    access: &dyn BlockAccess,
    state: &mut T,
    mut visit: impl FnMut(&mut T, Record<'_>) -> io::Result<()>,
) -> io::Result<()> {
    ForwardAccess::new(access.reader(0)).for_each_record(|record| visit(state, record))
}

/// Visit every record of a forward-only input, folding into `state`, in bounded memory.
pub fn for_each_record_in_stream<R: Read, T>(
    stream: &mut ForwardAccess<R>,
    state: &mut T,
    mut visit: impl FnMut(&mut T, Record<'_>) -> io::Result<()>,
) -> io::Result<()> {
    stream.for_each_record(|record| visit(state, record))
}

// endregion: Serial

// region: Parallel

/// One state and whatever went wrong in the work it ran.
///
/// The failure rides in the slot the worker already owns exclusively, so reporting it needs no
/// synchronization, and it carries the index of the work that failed so the error the user
/// sees is the earliest one rather than whichever thread lost the race.
struct Slot<'a, T> {
    state: &'a mut T,
    failure: Option<(usize, io::Error)>,
}

impl<T> Slot<'_, T> {
    /// Run one piece of work unless the round has already failed, keeping the earliest error.
    fn attempt(
        &mut self,
        cancelled: &AtomicBool,
        task_index: usize,
        run: &(impl Fn(&mut T, usize) -> io::Result<()> + ?Sized),
    ) {
        // One relaxed load per piece of work, never per record. Work already running finishes
        // and work not yet started declines, so a failure ends the round without unwinding
        // across threads.
        if cancelled.load(Ordering::Relaxed) {
            return;
        }
        if let Err(error) = run(self.state, task_index) {
            if earlier_than(&self.failure, task_index) {
                self.failure = Some((task_index, error));
            }
            cancelled.store(true, Ordering::Relaxed);
        }
    }
}

/// Whether `task_index` names work ahead of whatever `held` already reports.
fn earlier_than(held: &Option<(usize, io::Error)>, task_index: usize) -> bool {
    match held {
        Some((reported, _)) => task_index < *reported,
        None => true,
    }
}

/// One slot per state, or an error when a tool brought none.
fn slots_of<T>(states: &mut [T]) -> io::Result<Vec<Slot<'_, T>>> {
    if states.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "a run needs at least one state to fold into",
        ));
    }
    Ok(states
        .iter_mut()
        .map(|state| Slot {
            state,
            failure: None,
        })
        .collect())
}

/// Dispatch `tasks` pieces of work across the pool and report the earliest failure.
///
/// Indices are unique within a dispatch, so the aliasing obligation is discharged by
/// construction and the raw pointer work stays inside the dependency.
fn run_round<T: Send + Sync>(
    pool: &mut ThreadPool,
    order: RecordOrder,
    slots: &mut [Slot<'_, T>],
    tasks: usize,
    run: impl Fn(&mut T, usize) -> io::Result<()> + Send + Sync,
) -> io::Result<()> {
    if tasks == 0 {
        return Ok(());
    }
    let cancelled = AtomicBool::new(false);
    // Folding indexes states by thread, so a tool that brought fewer states than the pool has
    // threads is served by the task-indexed dispatch, which is correct for either order.
    let indexing = match slots.len() < pool.threads_count() {
        true => RecordOrder::Preserved,
        false => order,
    };
    match indexing {
        RecordOrder::Preserved => {
            forkunion::for_each_prong_mut_dynamic(pool, &mut slots[..tasks], |slot, prong| {
                slot.attempt(&cancelled, prong.task_index, &run);
            });
        }
        RecordOrder::Ignored => {
            forkunion::fold_with_scratch(
                pool,
                (0..tasks).into_par_iter(),
                DynamicScheduler,
                slots,
                |slot, task_index, _prong| slot.attempt(&cancelled, task_index, &run),
            );
        }
    }

    // Ascending order, so the earliest failing piece of work reports, deterministically.
    let mut earliest: Option<(usize, io::Error)> = None;
    for slot in slots.iter_mut() {
        if let Some((task_index, error)) = slot.failure.take() {
            if earlier_than(&earliest, task_index) {
                earliest = Some((task_index, error));
            }
        }
    }
    match earliest {
        Some((_, error)) => Err(error),
        None => Ok(()),
    }
}

/// End of the byte range opening at `from`: the first provable record boundary at least
/// `target_bytes` in, or the end of the input when none can be proven.
///
/// Multi-line FASTQ cannot be resynchronized from an arbitrary offset, so there the whole
/// remainder becomes one range and the run is serial in everything but name.
fn byte_range_end(bytes: &[u8], from: usize, target_bytes: usize, format: SequenceFormat) -> usize {
    match from.checked_add(target_bytes) {
        Some(target) if target < bytes.len() => {
            next_record_boundary(bytes, target, format).unwrap_or(bytes.len())
        }
        _ => bytes.len(),
    }
}

/// Visit every record of a buffer already in memory, one byte range per worker.
///
/// `states` is one instance per range in flight and its length is the memory bound: ranges
/// dispatch in rounds of that many, and a round retires before the next dispatches. `retire`
/// runs on the calling thread, in input order, and must leave its state ready for reuse.
pub fn for_each_record_in_bytes_in_parallel<T: Send + Sync>(
    bytes: &[u8],
    pool: &mut ThreadPool,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    in_byte_ranges(bytes, TARGET_WORK_BYTES, pool, order, states, visit, retire)
}

/// [`for_each_record_in_bytes_in_parallel`] with the range size stated, which the tests use to
/// force many rounds over a small buffer.
fn in_byte_ranges<T: Send + Sync>(
    bytes: &[u8],
    target_bytes: usize,
    pool: &mut ThreadPool,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    mut retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    let mut slots = slots_of(states)?;
    // A buffer with no records is empty, not malformed, so there is nothing to schedule.
    let Some(format) = sequence_format(bytes)? else {
        return Ok(());
    };

    let mut ranges: Vec<Range<usize>> = Vec::with_capacity(slots.len());
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        ranges.clear();
        while ranges.len() < slots.len() && cursor < bytes.len() {
            let end = byte_range_end(bytes, cursor, target_bytes, format);
            ranges.push(cursor..end);
            cursor = end;
        }
        let tasks = ranges.len();
        run_round(pool, order, &mut slots, tasks, |state, task_index| {
            for_each_record_in_bytes(&bytes[ranges[task_index].clone()], state, &visit)
        })?;
        if order == RecordOrder::Preserved {
            for slot in slots[..tasks].iter_mut() {
                retire(slot.state)?;
            }
        }
    }
    if order == RecordOrder::Ignored {
        for slot in slots.iter_mut() {
            retire(slot.state)?;
        }
    }
    Ok(())
}

/// The run of blocks opening at `first`, and how many decoded bytes the run itself owns.
///
/// A greedy prefix sum over decoded rather than compressed sizes, since poly-A compresses ten
/// to one against ordinary sequence and one worker would otherwise take ten times the work.
///
/// The budget is `None` for a run reaching the last block, which has nothing after it to read
/// into twice. A container that only bounds a block's decoded size cannot say where such a run
/// ends either, so that run takes every block left and the input reads on one worker.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
fn next_block_run(
    access: &dyn BlockAccess,
    first: usize,
    block_count: usize,
    target_bytes: usize,
) -> (Range<usize>, Option<usize>) {
    let mut decoded = 0usize;
    let mut exact = true;
    let mut last = first;
    // At least one block, so a block larger than the target still makes progress.
    while last < block_count && decoded < target_bytes {
        let bytes = access.bytes_in_block(last);
        exact &= bytes.exact().is_some();
        decoded += bytes.upper_bound();
        last += 1;
    }
    match last == block_count || !exact {
        true => (first..block_count, None),
        false => (first..last, Some(decoded)),
    }
}

/// Visit every record beginning in a run of blocks of one input.
///
/// A block is not a record, so the tail of a record the run before owns opens this one and is
/// skipped, while the record straddling the far end began here and is finished by decoding on
/// into the blocks that follow. `decoded_bytes` is what the run itself holds, known before
/// anything is decoded, and `None` reads to the end of the input.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
fn for_each_record_in_block_run(
    access: &dyn BlockAccess,
    first_block: usize,
    decoded_bytes: Option<usize>,
    format: SequenceFormat,
    visit: impl FnMut(Record<'_>) -> io::Result<()>,
) -> io::Result<()> {
    let start = match first_block {
        0 => StreamStart::FirstByte,
        _ => StreamStart::FirstBoundary(format),
    };
    ForwardAccess::with_window(access.reader(first_block), BLOCK_WINDOW_BYTES)
        .for_each_record_within(start, decoded_bytes, visit)
}

/// Format of the first record of a block-addressable input, or `None` when it holds none.
///
/// Decodes only the leading blocks a prefix takes, and is asked once on the calling thread so
/// every run of blocks can be told rather than decoding block zero to decide for itself.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
pub fn sequence_format_in_blocks(access: &dyn BlockAccess) -> io::Result<Option<SequenceFormat>> {
    let mut reader = access.reader(0);
    let mut prefix = [0u8; 64];
    let mut filled = 0usize;
    while filled < prefix.len() {
        match reader.read(&mut prefix[filled..])? {
            0 => break,
            bytes => filled += bytes,
        }
    }
    sequence_format(&prefix[..filled])
}

/// Visit every record of a block-addressable input, one run of blocks per worker.
///
/// The blocks are already enumerated, so no boundary search happens here at all, and the
/// format is detected once on the calling thread and told to every run.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
pub fn for_each_record_in_blocks_in_parallel<T: Send + Sync>(
    access: &dyn BlockAccess,
    pool: &mut ThreadPool,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    in_block_runs(
        access,
        TARGET_WORK_BYTES,
        pool,
        order,
        states,
        visit,
        retire,
    )
}

/// [`for_each_record_in_blocks_in_parallel`] with the run size stated, which the tests use to
/// force many rounds over a small input.
#[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
fn in_block_runs<T: Send + Sync>(
    access: &dyn BlockAccess,
    target_bytes: usize,
    pool: &mut ThreadPool,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    mut retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    let mut slots = slots_of(states)?;
    let Some(format) = sequence_format_in_blocks(access)? else {
        return Ok(());
    };
    let block_count = access.block_count();

    // One entry per run in flight: which blocks it covers, and the decoded bytes it owns.
    let mut runs: Vec<(Range<usize>, Option<usize>)> = Vec::with_capacity(slots.len());
    let mut first = 0usize;
    while first < block_count {
        runs.clear();
        while runs.len() < slots.len() && first < block_count {
            let (blocks, decoded_bytes) = next_block_run(access, first, block_count, target_bytes);
            first = blocks.end;
            runs.push((blocks, decoded_bytes));
        }
        let tasks = runs.len();
        run_round(pool, order, &mut slots, tasks, |state, task_index| {
            let (blocks, decoded_bytes) = &runs[task_index];
            for_each_record_in_block_run(access, blocks.start, *decoded_bytes, format, |record| {
                visit(state, record)
            })
        })?;
        if order == RecordOrder::Preserved {
            for slot in slots[..tasks].iter_mut() {
                retire(slot.state)?;
            }
        }
    }
    if order == RecordOrder::Ignored {
        for slot in slots.iter_mut() {
            retire(slot.state)?;
        }
    }
    Ok(())
}

/// Append to `slab` until it holds `wanted` bytes or the stream ends, saying which happened.
fn fill_slab(source: &mut impl Read, slab: &mut Vec<u8>, wanted: usize) -> io::Result<bool> {
    if slab.len() >= wanted {
        return Ok(false);
    }
    let mut filled = slab.len();
    slab.resize(wanted, 0);
    let mut at_eof = false;
    while filled < wanted {
        match source.read(&mut slab[filled..])? {
            0 => {
                at_eof = true;
                break;
            }
            read_bytes => filled += read_bytes,
        }
    }
    slab.truncate(filled);
    Ok(at_eof)
}

/// Visit every record of a forward-only input, one byte range of one slab per worker.
///
/// A stream cannot be split, so decoding stays serial: the calling thread fills a slab, ends it
/// at the last provable record boundary, and carries the tail into the next one. The parse and
/// the tool's work go wide over what is left behind.
pub fn for_each_record_in_stream_in_parallel<R: Read, T: Send + Sync>(
    stream: &mut ForwardAccess<R>,
    pool: &mut ThreadPool,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    in_slabs(
        stream,
        TARGET_WORK_BYTES,
        pool,
        order,
        states,
        visit,
        retire,
    )
}

/// [`for_each_record_in_stream_in_parallel`] with the range size stated, which the tests use to
/// force many slabs over a small stream.
fn in_slabs<R: Read, T: Send + Sync>(
    stream: &mut ForwardAccess<R>,
    target_bytes: usize,
    pool: &mut ThreadPool,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    mut retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    // A slab holds one round's worth of ranges, so the parse has something to spread.
    let slab_bytes = target_bytes.saturating_mul(states.len().max(1));
    let mut format: Option<SequenceFormat> = None;
    let mut slab: Vec<u8> = Vec::new();
    let mut wanted = slab_bytes;
    let source = stream.source_mut();

    loop {
        let at_eof = fill_slab(source, &mut slab, wanted)?;
        if format.is_none() {
            format = sequence_format(&slab)?;
        }
        let Some(known) = format else {
            // Nothing but whitespace so far: an empty stream is not an error, and a longer run
            // of it needs more bytes before it can say which format it carries.
            if at_eof {
                return Ok(());
            }
            wanted = slab.len() + target_bytes;
            continue;
        };

        let split = match at_eof {
            true => slab.len(),
            false => match last_record_boundary(&slab, known) {
                Some(boundary) if boundary > 0 => boundary,
                // At most one record start in the slab, so there is no boundary to end it on
                // and nothing complete to dispatch. Read more rather than splitting blind.
                _ => {
                    wanted = slab.len() + target_bytes;
                    continue;
                }
            },
        };

        in_byte_ranges(
            &slab[..split],
            target_bytes,
            pool,
            order,
            states,
            &visit,
            &mut retire,
        )?;
        if at_eof {
            return Ok(());
        }
        slab.drain(..split);
        wanted = slab_bytes;
    }
}

// endregion: Parallel

// region: Dispatch

/// Run one input through `visit`, on the pool when there is one and on this thread when there
/// is not.
///
/// `states` is what the tool folds into and `retire` empties one, on the calling thread and in
/// input order. A serial run uses the first state alone and retires it whenever a target's
/// worth of records has gone through, so a tool that buffers its output stays bounded.
pub fn for_each_record_in_input<T: Send + Sync>(
    input: &mut InputFile,
    workers: &mut Workers,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    let (pool, order) = workers.dispatch();
    in_one_input(input, pool, order, states, visit, retire)
}

/// [`for_each_record_in_input`] with the pool and the ordering already unpacked, so the driver
/// over many inputs can reuse one borrow of the workers across all of them.
fn in_one_input<T: Send + Sync>(
    input: &mut InputFile,
    pool: Option<&mut ThreadPool>,
    order: RecordOrder,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    mut retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    let Some(pool) = pool else {
        return serially(input, states, visit, retire);
    };
    match input {
        InputFile::Bytes(bytes) => for_each_record_in_bytes_in_parallel(
            bytes.as_slice(),
            pool,
            order,
            states,
            visit,
            &mut retire,
        ),
        #[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
        InputFile::Blocks(blocks) => for_each_record_in_blocks_in_parallel(
            &**blocks,
            pool,
            order,
            states,
            visit,
            &mut retire,
        ),
        InputFile::Stream(stream) => {
            for_each_record_in_stream_in_parallel(stream, pool, order, states, visit, &mut retire)
        }
    }
}

/// One input through one state, retiring it at a bounded cadence.
///
/// Serial order is input order, so retiring between records is always legal, and it is what
/// keeps memory bounded for a tool whose state buffers what it writes.
fn serially<T>(
    input: &mut InputFile,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()>,
    mut retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    let Some(state) = states.first_mut() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "a run needs at least one state to fold into",
        ));
    };

    let mut pending = 0usize;
    let mut visit_and_retire = |state: &mut T, record: Record<'_>| {
        let bytes = record.header.len() + record.sequence.len() + record.quality.len();
        visit(state, record)?;
        pending += bytes;
        if pending >= TARGET_WORK_BYTES {
            pending = 0;
            retire(state)?;
        }
        Ok(())
    };

    match input {
        InputFile::Bytes(bytes) => {
            for_each_record_in_bytes(bytes.as_slice(), state, &mut visit_and_retire)?
        }
        #[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
        InputFile::Blocks(blocks) => {
            for_each_record_in_blocks(&**blocks, state, &mut visit_and_retire)?
        }
        InputFile::Stream(stream) => {
            for_each_record_in_stream(stream, state, &mut visit_and_retire)?
        }
    }
    retire(state)
}

/// Run every input through `visit`, in argument order.
///
/// Inputs holding no records contribute nothing and are not an error. A record carries its own
/// format, so nothing has to be detected here or handed to the tool ahead of time; a tool that
/// mirrors its input calls [`crate::files::RecordWriter::adopt`], which also catches an input
/// whose format disagrees with an earlier one.
///
/// A named file is memory-mapped and read in place; `-` is read through a bounded window and
/// may appear at most once, because there is only one standard input to consume.
pub fn for_each_record_in_inputs<T: Send + Sync>(
    paths: &[String],
    workers: &mut Workers,
    states: &mut [T],
    visit: impl Fn(&mut T, Record<'_>) -> io::Result<()> + Send + Sync,
    mut retire: impl FnMut(&mut T) -> io::Result<()>,
) -> io::Result<()> {
    if paths.iter().filter(|path| path.as_str() == "-").count() > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "standard input can only be read once, but '-' was given more than once",
        ));
    }

    let (mut pool, order) = workers.dispatch();
    for path in paths {
        let mut input = InputFile::open(Some(path))?;
        in_one_input(
            &mut input,
            pool.as_deref_mut(),
            order,
            states,
            &visit,
            &mut retire,
        )?;
    }
    Ok(())
}

// endregion: Dispatch

#[cfg(test)]
mod tests {
    use super::*;

    fn collect_headers(data: &[u8]) -> Vec<Vec<u8>> {
        let mut headers = Vec::new();
        for_each_record_in_bytes(data, &mut headers, |headers, record| {
            headers.push(record.header.to_vec());
            Ok(())
        })
        .unwrap();
        headers
    }

    #[test]
    fn every_record_is_visited_in_order() {
        let headers = collect_headers(b">a\nACGT\n>b\nTGCA\n>c\nGGGG\n");
        assert_eq!(
            headers,
            vec![b">a".to_vec(), b">b".to_vec(), b">c".to_vec()]
        );
    }

    #[test]
    fn an_input_with_no_records_visits_nothing() {
        assert!(collect_headers(b"").is_empty());
    }

    /// There is only one standard input, so a second `-` cannot be honoured.
    #[test]
    fn standard_input_named_twice_is_rejected() {
        let paths = vec!["-".to_string(), "-".to_string()];
        let mut workers = Parallelism { jobs: 1 }.ordered().unwrap();
        let mut states = [()];
        let error =
            for_each_record_in_inputs(&paths, &mut workers, &mut states, |_, _| Ok(()), |_| Ok(()))
                .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
    }

    fn sample_fastq(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            let width = 8 + (index % 5) * 4;
            data.extend_from_slice(format!("@read{index} sample\n").as_bytes());
            data.extend(std::iter::repeat_n(b"ACGT"[index % 4], width));
            data.extend_from_slice(b"\n+\n");
            // A quality line opening with '@' is legal at Phred 31, and is exactly what a
            // boundary search must not mistake for a header.
            data.push(b'@');
            data.extend(std::iter::repeat_n(b'I', width - 1));
            data.push(b'\n');
        }
        data
    }

    fn sample_fasta(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            let width = 9 + (index % 5) * 4;
            data.extend_from_slice(format!(">contig{index} sample\n").as_bytes());
            for chunk_start in (0..width).step_by(6) {
                data.extend(std::iter::repeat_n(
                    b"ACGT"[index % 4],
                    (width - chunk_start).min(6),
                ));
                data.push(b'\n');
            }
        }
        data
    }

    /// Header and sequence of one record, which is what an order comparison needs.
    fn render(state: &mut Vec<u8>, record: Record<'_>) -> io::Result<()> {
        state.extend_from_slice(record.header);
        state.push(b'\n');
        state.extend_from_slice(record.sequence);
        state.push(b'\n');
        Ok(())
    }

    fn rendered_serially(data: &[u8]) -> Vec<u8> {
        let mut state = Vec::new();
        for_each_record_in_bytes(data, &mut state, render).unwrap();
        state
    }

    fn pool_of(threads: usize) -> (Topology, ThreadPool) {
        let topology = Topology::new().unwrap();
        let pool = ThreadPool::try_spawn(&topology, threads).unwrap();
        (topology, pool)
    }

    /// Every state's output, concatenated as it was retired.
    fn rendered_in_parallel(
        data: &[u8],
        target_bytes: usize,
        threads: usize,
        states_count: usize,
    ) -> io::Result<Vec<u8>> {
        let (_topology, mut pool) = pool_of(threads);
        let mut states: Vec<Vec<u8>> = vec![Vec::new(); states_count];
        let mut written = Vec::new();
        in_byte_ranges(
            data,
            target_bytes,
            &mut pool,
            RecordOrder::Preserved,
            &mut states,
            render,
            |state| {
                written.extend_from_slice(state);
                state.clear();
                Ok(())
            },
        )?;
        Ok(written)
    }

    /// The acceptance property: the byte stream a parallel run produces is the byte stream a
    /// serial run produces, at every thread count and every number of states.
    #[test]
    fn parallel_output_matches_serial_output() {
        for data in [sample_fasta(400), sample_fastq(400)] {
            let expected = rendered_serially(&data);
            for threads in [1usize, 2, 3, 8] {
                for states_count in [1usize, 2, 5, 16] {
                    for target_bytes in [64usize, 512, 4096] {
                        let produced =
                            rendered_in_parallel(&data, target_bytes, threads, states_count)
                                .unwrap();
                        assert_eq!(
                            produced, expected,
                            "{threads} threads, {states_count} states, {target_bytes}-byte ranges"
                        );
                    }
                }
            }
        }
    }

    /// Folding indexes states by thread rather than by range, so a tool that merges at the end
    /// still sees every record exactly once.
    #[test]
    fn folding_by_thread_sees_every_record_once() {
        let data = sample_fastq(400);
        let expected = collect_headers(&data).len();
        for threads in [1usize, 2, 8] {
            let (_topology, mut pool) = pool_of(threads);
            let mut states: Vec<usize> = vec![0; threads];
            in_byte_ranges(
                &data,
                256,
                &mut pool,
                RecordOrder::Ignored,
                &mut states,
                |state, _record| {
                    *state += 1;
                    Ok(())
                },
                |_| Ok(()),
            )
            .unwrap();
            assert_eq!(states.iter().sum::<usize>(), expected, "{threads} threads");
        }
    }

    #[test]
    fn empty_input_is_not_an_error() {
        for threads in [1usize, 4] {
            assert!(rendered_in_parallel(b"", 64, threads, 4)
                .unwrap()
                .is_empty());
        }
    }

    /// `next_record_boundary` cannot resynchronize a record wrapped across several lines, so
    /// the whole input becomes one range. Benign, but it is a performance cliff worth naming:
    /// a wrapped FASTQ runs on one worker however many were asked for.
    #[test]
    fn wrapped_fastq_falls_back_to_a_single_range() {
        let mut data = Vec::new();
        for index in 0..200 {
            data.extend_from_slice(format!("@read{index}\nACGT\nACGT\n+\nIIII\nIIII\n").as_bytes());
        }
        let expected = rendered_serially(&data);
        assert_eq!(
            byte_range_end(&data, 0, 64, SequenceFormat::Fastq),
            data.len(),
            "a wrapped record must not yield a boundary to split on"
        );
        assert_eq!(rendered_in_parallel(&data, 64, 4, 8).unwrap(), expected);
    }

    /// Two defects, and the earlier one must win every time. A shared error cell would report
    /// whichever thread reached it first, which is a different answer run to run.
    #[test]
    fn the_earliest_failing_range_reports_every_time() {
        let mut data = sample_fastq(200);
        // A quality line shorter than its sequence is malformed wherever it lands.
        data.extend_from_slice(b"@early\nACGTACGT\n+\nII\n");
        data.extend_from_slice(&sample_fastq(200));
        data.extend_from_slice(b"@late\nACGTACGT\n+\nI\n");
        data.extend_from_slice(&sample_fastq(200));

        let mut serial_state = Vec::new();
        let expected = for_each_record_in_bytes(&data, &mut serial_state, render)
            .unwrap_err()
            .to_string();
        for _ in 0..16 {
            let error = rendered_in_parallel(&data, 128, 8, 8).unwrap_err();
            assert_eq!(error.to_string(), expected);
        }
    }

    /// One BGZF block over `payload`: gzip magic, a `BC` extra field, raw deflate, trailer.
    #[cfg(feature = "gzip")]
    fn bgzf_block(payload: &[u8]) -> Vec<u8> {
        let mut deflated = Vec::with_capacity(payload.len() + 64);
        let mut compress = flate2::Compress::new(flate2::Compression::fast(), false);
        compress
            .compress_vec(payload, &mut deflated, flate2::FlushCompress::Finish)
            .unwrap();
        let mut checksum = flate2::Crc::new();
        checksum.update(payload);

        let size = 12 + 6 + deflated.len() + 8;
        let mut block = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0];
        block.extend_from_slice(b"BC");
        block.extend_from_slice(&2u16.to_le_bytes());
        block.extend_from_slice(&((size - 1) as u16).to_le_bytes());
        block.extend_from_slice(&deflated);
        block.extend_from_slice(&checksum.sum().to_le_bytes());
        block.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        block
    }

    /// `plain` as BGZF in blocks of `payload_bytes`, small enough that records straddle them.
    #[cfg(feature = "gzip")]
    fn bgzf_bytes(plain: &[u8], payload_bytes: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for chunk in plain.chunks(payload_bytes) {
            out.extend_from_slice(&bgzf_block(chunk));
        }
        out.extend_from_slice(&bgzf_block(b""));
        out
    }

    /// A block-addressable view of `compressed`, through the seam an opened input uses, so no
    /// test here names a concrete container.
    #[cfg(feature = "gzip")]
    fn block_access_of(compressed: Vec<u8>) -> Box<dyn BlockAccess> {
        let Ok(access) = crate::blocks::block_access(compressed, crate::codecs::Codec::Bgzf) else {
            panic!("the fixture must be block-addressable");
        };
        access
    }

    /// Decoded bytes the run `first..last` owns, as the scheduler computes them.
    #[cfg(feature = "gzip")]
    fn run_budget(access: &dyn BlockAccess, first: usize, last: usize) -> Option<usize> {
        match last == access.block_count() {
            true => None,
            false => Some(
                (first..last)
                    .map(|index| access.bytes_in_block(index).upper_bound())
                    .sum(),
            ),
        }
    }

    /// Records of one fixed width, in blocks holding a whole number of them, so a run's far
    /// end lands exactly on a record start. That is the case the stop rule is exact about and
    /// the ragged fixtures never reach.
    #[cfg(feature = "gzip")]
    fn aligned_fastq(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            data.extend_from_slice(format!("@read{index:04}\nACGTACGT\n+\nIIIIIIII\n").as_bytes());
        }
        data
    }

    /// Consecutive runs of blocks tile the input: every record appears exactly once and in
    /// order, whatever the run length. A record straddling a run's far end is emitted by the
    /// run it began in, and the run after skips the tail it inherits — miss either rule and
    /// this drops a record or emits one twice.
    #[cfg(feature = "gzip")]
    #[test]
    fn runs_of_blocks_tile_the_input() {
        // The last pair is aligned: 30 bytes per record, 8 records per block.
        let fixtures = [
            (sample_fastq(300), 256usize),
            (sample_fasta(300), 256),
            (aligned_fastq(300), 240),
        ];
        for (plain, payload_bytes) in fixtures {
            let access = block_access_of(bgzf_bytes(&plain, payload_bytes));
            let format = sequence_format_in_blocks(&*access).unwrap().unwrap();
            let expected = collect_headers(&plain);
            assert!(access.block_count() > 16, "fixture must hold many blocks");

            for run in [1, 2, 3, 7, access.block_count()] {
                let mut headers = Vec::new();
                let mut first = 0usize;
                while first < access.block_count() {
                    let last = (first + run).min(access.block_count());
                    let budget = run_budget(&*access, first, last);
                    for_each_record_in_block_run(&*access, first, budget, format, |record| {
                        headers.push(record.header.to_vec());
                        Ok(())
                    })
                    .unwrap();
                    first = last;
                }
                assert_eq!(
                    headers, expected,
                    "{format:?} in {payload_bytes}-byte blocks, runs of {run}"
                );
            }
        }
    }

    /// A run holding no record start contributes nothing rather than running to the end of the
    /// file, which is what stops one long record being emitted once per run it spans.
    #[cfg(feature = "gzip")]
    #[test]
    fn a_run_inside_one_long_record_yields_nothing() {
        let mut plain = Vec::from(&b">long\n"[..]);
        plain.extend(std::iter::repeat_n(b'A', 8192));
        plain.push(b'\n');
        let access = block_access_of(bgzf_bytes(&plain, 256));

        let mut headers = Vec::new();
        let budget = run_budget(&*access, 4, 8);
        for_each_record_in_block_run(&*access, 4, budget, SequenceFormat::Fasta, |record| {
            headers.push(record.header.to_vec());
            Ok(())
        })
        .unwrap();
        assert!(headers.is_empty());
    }

    /// Reads long enough that a four-line proof outruns the parse window.
    ///
    /// A 20 kb nanopore read spans several 64 KB blocks and its cycle is tens of kilobytes, so
    /// the boundary search a run does at its first block cannot finish inside one window. This
    /// is the shape that dropped 41 of 2600 records while every short-read fixture passed.
    #[cfg(feature = "gzip")]
    #[test]
    fn long_reads_survive_a_boundary_search_that_outruns_the_window() {
        let mut plain = Vec::new();
        for index in 0..24 {
            let length = 4000 + index * 500;
            plain.extend_from_slice(format!("@read{index} length={length}\n").as_bytes());
            plain.extend(std::iter::repeat_n(b'A', length));
            plain.extend_from_slice(b"\n+\n");
            plain.extend(std::iter::repeat_n(b'I', length));
            plain.push(b'\n');
        }
        let expected = rendered_serially(&plain);
        let access = block_access_of(bgzf_bytes(&plain, 4096));

        for threads in [1usize, 2, 4, 8] {
            let (_topology, mut pool) = pool_of(threads);
            let mut states: Vec<Vec<u8>> = vec![Vec::new(); threads * 2];
            let mut written = Vec::new();
            in_block_runs(
                &*access,
                8192,
                &mut pool,
                RecordOrder::Preserved,
                &mut states,
                render,
                |state| {
                    written.extend_from_slice(state);
                    state.clear();
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(
                written, expected,
                "lost or duplicated records at {threads} threads"
            );
        }
    }

    /// The block tier reaches the same bytes as a plain buffer, whatever the run size and
    /// however many workers share the input.
    #[cfg(feature = "gzip")]
    #[test]
    fn blocks_in_parallel_match_the_serial_bytes() {
        for plain in [sample_fasta(300), sample_fastq(300)] {
            let expected = rendered_serially(&plain);
            let access = block_access_of(bgzf_bytes(&plain, 256));
            for threads in [1usize, 2, 3, 8] {
                for states_count in [1usize, 3, 8] {
                    let (_topology, mut pool) = pool_of(threads);
                    let mut states: Vec<Vec<u8>> = vec![Vec::new(); states_count];
                    let mut written = Vec::new();
                    in_block_runs(
                        &*access,
                        512,
                        &mut pool,
                        RecordOrder::Preserved,
                        &mut states,
                        render,
                        |state| {
                            written.extend_from_slice(state);
                            state.clear();
                            Ok(())
                        },
                    )
                    .unwrap();
                    assert_eq!(
                        written, expected,
                        "{threads} threads, {states_count} states"
                    );
                }
            }
        }
    }

    /// A stream is read once from front to back, so the only question is whether the slab
    /// boundaries lose or repeat a record.
    #[test]
    fn slabs_in_parallel_match_the_serial_bytes() {
        for data in [sample_fasta(400), sample_fastq(400)] {
            let expected = rendered_serially(&data);
            for threads in [1usize, 2, 8] {
                for states_count in [1usize, 3, 8] {
                    for target_bytes in [64usize, 700] {
                        let (_topology, mut pool) = pool_of(threads);
                        let mut states: Vec<Vec<u8>> = vec![Vec::new(); states_count];
                        let mut written = Vec::new();
                        let mut stream = ForwardAccess::new(&data[..]);
                        in_slabs(
                            &mut stream,
                            target_bytes,
                            &mut pool,
                            RecordOrder::Preserved,
                            &mut states,
                            render,
                            |state| {
                                written.extend_from_slice(state);
                                state.clear();
                                Ok(())
                            },
                        )
                        .unwrap();
                        assert_eq!(
                            written, expected,
                            "{threads} threads, {states_count} states, {target_bytes}-byte slabs"
                        );
                    }
                }
            }
        }
    }

    /// A stream holding nothing but whitespace ends without a format and without an error.
    #[test]
    fn a_stream_of_whitespace_is_not_an_error() {
        let data = vec![b'\n'; 4096];
        let (_topology, mut pool) = pool_of(2);
        let mut states: Vec<usize> = vec![0; 2];
        let mut stream = ForwardAccess::new(&data[..]);
        in_slabs(
            &mut stream,
            64,
            &mut pool,
            RecordOrder::Preserved,
            &mut states,
            |state, _record| {
                *state += 1;
                Ok(())
            },
            |_| Ok(()),
        )
        .unwrap();
        assert_eq!(states.iter().sum::<usize>(), 0);
    }

    /// A tool that brings no states has nothing to fold into, which is a mistake worth naming
    /// rather than a run that quietly visits nothing.
    #[test]
    fn a_run_without_states_is_rejected() {
        let (_topology, mut pool) = pool_of(2);
        let mut states: [Vec<u8>; 0] = [];
        let error = in_byte_ranges(
            b">a\nACGT\n",
            64,
            &mut pool,
            RecordOrder::Preserved,
            &mut states,
            render,
            |_| Ok(()),
        )
        .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
    }

    /// The serial path retires in pieces too, and must produce the same bytes as one pass.
    #[test]
    fn a_serial_run_retires_in_bounded_pieces() {
        let data = sample_fasta(400);
        let expected = rendered_serially(&data);
        let mut states = vec![Vec::new()];
        let mut written = Vec::new();
        let mut stream = ForwardAccess::new(&data[..]);
        let source: Box<dyn Read> = Box::new(io::Cursor::new(data.clone()));
        let mut input = InputFile::Stream(ForwardAccess::new(source));
        let mut workers = Parallelism { jobs: 1 }.ordered().unwrap();
        for_each_record_in_input(&mut input, &mut workers, &mut states, render, |state| {
            written.extend_from_slice(state);
            state.clear();
            Ok(())
        })
        .unwrap();
        assert_eq!(written, expected);
        let mut counted = 0usize;
        for_each_record_in_stream(&mut stream, &mut counted, |state, _record| {
            *state += 1;
            Ok(())
        })
        .unwrap();
        assert_eq!(counted, collect_headers(&data).len());
    }
}
