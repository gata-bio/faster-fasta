//! Collapse duplicate records, keeping the first occurrence of each.
//!
//! Deduplication is shaped differently from every other tool here, and the reason is where
//! the time goes. The others spend it in a per-record transform, so the shared driver puts
//! that on the pool and lets a serial fold drain a buffer. Here the transform is one hash
//! and the expense is probing an index that spans the whole input. So the worker does
//! everything that can be done without the index — parse, canonicalize, hash, render — and
//! hands over a flat batch, leaving the fold to decide and nothing else.
//!
//! Identity is settled by comparing bytes, never by a digest alone. A digest that matches is
//! a reason to look; only the bytes are a reason to drop. That costs one comparison per true
//! duplicate and roughly one per two billion distinct records otherwise, which is why a
//! 64-bit digest with verification beats a 128-bit digest without.
//!
//! __Memory__: the distinct keys themselves, plus 32 to 64 bytes of index per distinct
//! record depending on how full the table is, and 8 more per key for its end offset.
//! Nothing scales with total input size, only with the number of survivors.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fasta-dedup sequences.fasta --output unique.fasta
//! fasta-dedup reads.fastq --by name --output unique.fastq
//! fasta-dedup contigs.fa --by canonical --threads 16 --output unique.fa
//! cat sequences.fasta | fasta-dedup > unique.fasta
//! ```
//!
//! Exit: 0 retained a record, 1 ran and retained none, 2 could not run.

use std::io::{self, Write};

use clap::Parser;
use stringzilla::sz::{equal, hash_with_seed, lookup};

use fasterfasta::files::{
    finish_or_exit, Destination, InputFile, Presentation, RecordWriter, Rendering, RunOutcome,
};
use fasterfasta::records::{key_bytes, Complements, Record, RecordKey, SequenceFormat};
use fasterfasta::scheduling::{
    for_each_record_in_input, for_each_record_to_destination, Parallelism, Workers,
};

// region: Identity

/// Seed for every digest this tool takes.
///
/// Fixed rather than random so a run is reproducible across machines and across passes.
const SEED: u64 = 0;

/// Which part of a record decides whether two are the same record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
enum By {
    /// The sequence bytes, whitespace already removed.
    Sequence,
    /// The identifier, so two records with one name collapse whatever they hold.
    Name,
    /// The sequence or its reverse complement, whichever sorts first.
    Canonical,
}

impl From<By> for Identity {
    fn from(by: By) -> Self {
        match by {
            By::Sequence => Identity::Field(RecordKey::Sequence),
            By::Name => Identity::Field(RecordKey::Identifier),
            By::Canonical => Identity::CanonicalStrand,
        }
    }
}

/// What decides whether two records are the same record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Identity {
    /// One field of the record, taken verbatim.
    Field(RecordKey),
    /// The sequence or its reverse complement, whichever sorts first.
    ///
    /// A read and its reverse complement are the same molecule, so a corpus holding both is
    /// holding the same information twice. Taking the lesser of the two is what makes the
    /// choice independent of which strand happened to be sequenced.
    CanonicalStrand,
}

/// Reverse complement of `sequence`, written into `scratch`.
fn reverse_complement_into(sequence: &[u8], complements: &Complements, scratch: &mut Vec<u8>) {
    let table = complements.table_for(sequence);
    scratch.clear();
    scratch.resize(sequence.len(), 0);
    lookup(scratch.as_mut_slice(), sequence, *table);
    scratch.reverse();
}

// endregion: Identity

// region: Retention

/// Index of a retained key, and the value marking a slot nothing has claimed.
///
/// A tape index rather than a pointer, so growing the tape cannot invalidate the table.
type KeyIndex = u32;
const UNCLAIMED: KeyIndex = KeyIndex::MAX;

/// Retained keys, laid end to end in the order they were first seen.
///
/// One allocation for the bytes and one for the ends, so a million survivors cost two
/// allocations rather than a million.
#[derive(Debug, Default)]
struct KeyTape {
    bytes: Vec<u8>,
    ends: Vec<u64>,
}

impl KeyTape {
    /// An empty tape holding `entries` keys over `bytes` before either has to reallocate.
    fn with_capacity(entries: usize, bytes: usize) -> Self {
        Self {
            bytes: Vec::with_capacity(bytes),
            ends: Vec::with_capacity(entries),
        }
    }

    /// Append `key`, returning the index it can be read back at.
    fn push(&mut self, key: &[u8]) -> io::Result<KeyIndex> {
        let index = self.ends.len();
        if index >= UNCLAIMED as usize {
            return Err(io::Error::other(
                "more distinct records than this build can index; split the input and merge",
            ));
        }
        self.bytes.extend_from_slice(key);
        self.ends.push(self.bytes.len() as u64);
        Ok(index as KeyIndex)
    }

    /// The key stored at `index`.
    fn get(&self, index: KeyIndex) -> &[u8] {
        let end = self.ends[index as usize] as usize;
        let start = match index {
            0 => 0,
            _ => self.ends[index as usize - 1] as usize,
        };
        &self.bytes[start..end]
    }

    /// How many keys are retained.
    #[cfg(test)]
    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.ends.len()
    }
}

/// One table slot: a digest and where the key it stands for is stored.
#[derive(Debug, Clone, Copy)]
struct Slot {
    tag: u64,
    key: KeyIndex,
}

const VACANT: Slot = Slot {
    tag: 0,
    key: UNCLAIMED,
};

/// Digest to retained key, open-addressed with linear probing.
///
/// A slot whose tag matches but whose bytes differ is a collision like any other, so probing
/// continues and both records survive. A map keyed on the digest would overwrite one of them,
/// which at a billion distinct records is not a hypothetical: it is a silently dropped
/// sequence roughly once every two billion.
#[derive(Debug)]
struct Retention {
    slots: Vec<Slot>,
    live: usize,
}

impl Retention {
    /// An empty table that holds `expected` distinct keys before the first growth.
    ///
    /// Twice `expected` slots, because the table grows at half full and probing degrades
    /// sharply past that.
    fn with_capacity(expected: usize) -> io::Result<Self> {
        let slots = expected.saturating_mul(2).max(16).next_power_of_two();
        // Fallible, because the count comes from a flag: a mistyped `--expect-distinct` must
        // cost an error rather than an allocation the kernel answers by killing the process.
        let mut table = Vec::new();
        table.try_reserve_exact(slots).map_err(|_| {
            io::Error::other(format!(
                "cannot pre-size an index for {expected} distinct records, which needs {} GB",
                slots * size_of::<Slot>() / 1_000_000_000
            ))
        })?;
        table.resize(slots, VACANT);
        Ok(Self {
            slots: table,
            live: 0,
        })
    }

    /// Retain `key` unless a record with the same bytes is already retained, saying which.
    ///
    /// `tape` is written only when the key is new, so a duplicate costs a probe and a
    /// comparison and no allocation at all.
    fn retain(&mut self, tag: u64, key: &[u8], tape: &mut KeyTape) -> io::Result<bool> {
        // Grown before the probe rather than after, so the position found is a position in
        // the table that will still exist when it is written to.
        if (self.live + 1) * 2 > self.slots.len() {
            self.grow();
        }
        let mask = self.slots.len() - 1;
        let mut at = tag as usize & mask;
        loop {
            let slot = self.slots[at];
            if slot.key == UNCLAIMED {
                let index = tape.push(key)?;
                self.slots[at] = Slot { tag, key: index };
                self.live += 1;
                return Ok(true);
            }
            if slot.tag == tag && equal(tape.get(slot.key), key) {
                return Ok(false);
            }
            at = (at + 1) & mask;
        }
    }

    /// Double the table, rehoming every live slot.
    fn grow(&mut self) {
        let mut bigger = vec![VACANT; self.slots.len() * 2];
        let mask = bigger.len() - 1;
        for slot in self.slots.iter().filter(|slot| slot.key != UNCLAIMED) {
            let mut at = slot.tag as usize & mask;
            while bigger[at].key != UNCLAIMED {
                at = (at + 1) & mask;
            }
            bigger[at] = *slot;
        }
        self.slots = bigger;
    }
}

// endregion: Retention

// region: Scanning

/// The range entry `index` occupies, given the end offsets of every entry before it.
///
/// Ends rather than starts, so appending is one push and the last end is also the length.
fn span(ends: &[u32], index: usize) -> std::ops::Range<usize> {
    let start = if index == 0 {
        0
    } else {
        ends[index - 1] as usize
    };
    start..ends[index] as usize
}

/// One unit of work's scan: every record keyed, hashed and rendered once.
///
/// Everything here is a pure function of one record, so it runs on the pool with no
/// coordination. What it cannot do is decide, because deciding needs the index.
struct Scanned {
    writer: RecordWriter<io::Sink>,
    rendering: Rendering,
    identity: Identity,
    complements: Complements,
    /// Holds a reverse complement while it is being compared against the forward strand.
    strand: Vec<u8>,
    /// Format of the first record seen, so the fold can make the same cross-input check.
    format: Option<SequenceFormat>,
    keys: Vec<u8>,
    key_ends: Vec<u32>,
    tags: Vec<u64>,
    rendered: Vec<u8>,
    record_ends: Vec<u32>,
}

impl Scanned {
    fn new(identity: Identity, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(io::sink(), rendering, SequenceFormat::Fasta),
            rendering,
            identity,
            complements: Complements::new(),
            strand: Vec::new(),
            format: None,
            keys: Vec::new(),
            key_ends: Vec::new(),
            tags: Vec::new(),
            rendered: Vec::new(),
            record_ends: Vec::new(),
        }
    }

    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        self.writer.adopt(record.format())?;
        self.format = Some(record.format());

        // Destructured so the key may borrow from `strand` while the batch borrows itself.
        let Self {
            writer,
            identity,
            complements,
            strand,
            keys,
            key_ends,
            tags,
            rendered,
            record_ends,
            ..
        } = self;

        let key = match *identity {
            Identity::Field(field) => key_bytes(&record, field),
            Identity::CanonicalStrand => {
                reverse_complement_into(record.sequence, complements, strand);
                if strand.as_slice() < record.sequence {
                    strand.as_slice()
                } else {
                    record.sequence
                }
            }
        };

        tags.push(hash_with_seed(key, SEED));
        keys.extend_from_slice(key);
        writer.append_into(record, rendered);

        // A unit of work is a few megabytes, so 32-bit offsets are ample and halve what the
        // batch costs. The assertions catch the day someone raises the unit size instead of
        // letting an offset wrap into a slice of the wrong record.
        debug_assert!(
            keys.len() <= u32::MAX as usize,
            "keys outgrew a 32-bit offset"
        );
        debug_assert!(
            rendered.len() <= u32::MAX as usize,
            "rendered records outgrew a 32-bit offset"
        );
        key_ends.push(keys.len() as u32);
        record_ends.push(rendered.len() as u32);
        Ok(())
    }

    /// How many records this batch holds.
    fn len(&self) -> usize {
        self.tags.len()
    }

    /// The key and the rendered bytes of record `index` within this batch.
    fn record(&self, index: usize) -> (&[u8], &[u8]) {
        (
            &self.keys[span(&self.key_ends, index)],
            &self.rendered[span(&self.record_ends, index)],
        )
    }

    /// Empty the batch, keeping every allocation for the next unit of work.
    ///
    /// The format goes with it: a batch reports the format of the records it holds, and an
    /// emptied batch holds none. Carrying the last one forward would let a run over a FASTA
    /// input and a FASTQ input agree where it should have refused.
    fn clear(&mut self) {
        // The writer is rebuilt rather than reused, because `adopt` latches the first format
        // it sees and a state outlives the input it first saw one in.
        self.writer =
            RecordWriter::with_rendering(io::sink(), self.rendering, SequenceFormat::Fasta);
        self.format = None;
        self.keys.clear();
        self.key_ends.clear();
        self.tags.clear();
        self.rendered.clear();
        self.record_ends.clear();
    }
}

// endregion: Scanning

// region: Deduplication

/// Distinct records the index is sized for absent `--expect-distinct`.
///
/// Small on purpose: a streamed input has no size to read off, so a 10 KB file must not pay
/// for a guess aimed at a 10 GB one.
const EXPECTED_DISTINCT: u32 = 1 << 16;

/// Bytes reserved per expected key, near the length of a short read.
///
/// A seed for the first allocation and nothing else; what a run may hold is measured.
const EXPECTED_KEY_BYTES: usize = 64;

/// Most slices a bounded run will cut the digest space into.
///
/// A ceiling too small to hold one slice even here is refused rather than chased forever.
const MOST_PASSES: u32 = 1 << 20;

/// Distinct records to size the index for.
fn expected_distinct(asked: Option<u32>) -> usize {
    asked.unwrap_or(EXPECTED_DISTINCT) as usize
}

/// Bytes named by a size like `48G`, `512M`, `1K`, or a plain count.
fn memory_bytes(asked: &str) -> io::Result<usize> {
    let (digits, scale) = match asked.as_bytes().last() {
        Some(b'K' | b'k') => (&asked[..asked.len() - 1], 1 << 10),
        Some(b'M' | b'm') => (&asked[..asked.len() - 1], 1 << 20),
        Some(b'G' | b'g') => (&asked[..asked.len() - 1], 1 << 30),
        _ => (asked, 1),
    };
    digits
        .trim()
        .parse::<usize>()
        .ok()
        .and_then(|count| count.checked_mul(scale))
        .filter(|bytes| *bytes > 0)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("--memory {asked} is not a size; write 48G, 512M, or a byte count"),
            )
        })
}

/// How many passes a bounded run starts at, read off the pre-sizing hint alone.
///
/// A guess and only a guess: the run weighs what it holds and doubles this whenever the
/// ceiling is crossed, so a hint that is wrong costs restarts rather than the ceiling.
/// A power of two, because a pass is named by the high bits of a digest.
fn passes_for(ceiling: usize, expected: usize) -> u32 {
    // What one pass would hold: the keys themselves plus the index that finds them.
    let wanted = expected.saturating_mul(EXPECTED_KEY_BYTES + size_of::<Slot>() * 2 + 8);
    let passes = wanted.div_ceil(ceiling.max(1)).max(1);
    passes.next_power_of_two().min(MOST_PASSES as usize) as u32
}

/// A decide pass held more than `--memory` allows, so the run doubles its passes and retries.
///
/// A type rather than a message, so a restart recognizes its own signal and no other error.
#[derive(Debug)]
struct CeilingCrossed;

impl std::fmt::Display for CeilingCrossed {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("a decide pass held more than --memory allows")
    }
}

impl std::error::Error for CeilingCrossed {}

/// Whether `error` is a decide pass reporting that it outgrew the ceiling.
fn crossed_the_ceiling(error: &io::Error) -> bool {
    error
        .get_ref()
        .is_some_and(|inner| inner.is::<CeilingCrossed>())
}

/// Every input can be read again, which a bounded run needs and a pipe cannot promise.
///
/// Checked before any work rather than trusted to fail: a second read of a forward-only input
/// yields no records and no error, so an unchecked run would quietly emit everything.
fn every_input_is_rereadable(paths: &[String]) -> io::Result<()> {
    for path in paths {
        if matches!(InputFile::open(Some(path))?, InputFile::Stream(_)) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "--memory needs to read '{path}' more than once, but it is a stream; \
                     write it to a file first, or recompress it with bgzip"
                ),
            ));
        }
    }
    Ok(())
}

/// Which records of one input turned out to be duplicates.
///
/// A bit per record rather than a list, because the emit pass asks about every record in turn
/// and a billion of them is 125 MB either way it is answered.
#[derive(Debug, Default)]
struct DropBits {
    words: Vec<u64>,
}

impl DropBits {
    fn mark(&mut self, position: usize) {
        let word = position / 64;
        if word >= self.words.len() {
            self.words.resize(word + 1, 0);
        }
        self.words[word] |= 1 << (position % 64);
    }

    fn marked(&self, position: usize) -> bool {
        self.words
            .get(position / 64)
            .is_some_and(|word| word >> (position % 64) & 1 == 1)
    }

    /// Forget every mark, keeping the allocation for the attempt that follows.
    fn clear(&mut self) {
        self.words.clear();
    }
}

/// What a run over the inputs is doing with the records it reads.
///
/// A bounded run decides one slice of the digest space per pass, holding only that slice's
/// keys, and writes nothing until every slice has been decided.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    /// Deciding and writing in one read, which is what an unbounded run does.
    Immediate,
    /// Deciding the records whose digest falls in this pass, marking the duplicates.
    ///
    /// The ceiling rides along because only this phase can answer for it: a pass that crosses
    /// it has written nothing, so it can be abandoned and retried at twice the passes.
    Deciding { pass: u32, ceiling: usize },
    /// Writing the records no pass marked, in input order.
    Emitting,
}

/// The global answer: which keys have been retained.
///
/// Identity spans every input, so this outlives them all. Where a survivor is written does
/// not, which is why the writer is an argument rather than a field.
struct Deduplicator {
    index: Retention,
    tape: KeyTape,
    /// Distinct records the pre-sizing hint names, split across whatever passes a run takes.
    hinted: usize,
    /// The most the index and the keys have held during the attempt now running.
    held_peak: usize,
    /// How far a digest is shifted to name its pass; zero when one pass sees everything.
    pass_shift: u32,
    phase: Phase,
    /// Records seen so far in the input being read, which indexes its `DropBits`.
    position: usize,
    examined: u64,
    retained: u64,
}

impl Deduplicator {
    /// A fold deciding and writing in one read, pre-sized for `expected` distinct records.
    fn new(expected: usize) -> io::Result<Self> {
        Ok(Self {
            index: Retention::with_capacity(expected)?,
            tape: KeyTape::with_capacity(expected, expected.saturating_mul(EXPECTED_KEY_BYTES)),
            hinted: expected,
            held_peak: 0,
            pass_shift: 0,
            phase: Phase::Immediate,
            position: 0,
            examined: 0,
            retained: 0,
        })
    }

    /// Which pass owns `tag`.
    ///
    /// The high bits, because the table indexes on the low ones: sharing bits between the two
    /// would leave every pass using only a slice of its own slots.
    fn pass_of(&self, tag: u64) -> u32 {
        match self.pass_shift {
            0 => 0,
            shift => (tag >> shift) as u32,
        }
    }

    /// Bytes the index and the retained keys hold right now, which is what `--memory` bounds.
    ///
    /// Capacity rather than length, because reserved memory is memory the machine cannot give
    /// to anything else, so a pass trips the ceiling on what it took rather than what it filled.
    fn held_bytes(&self) -> usize {
        self.index.slots.capacity() * size_of::<Slot>()
            + self.tape.bytes.capacity()
            + self.tape.ends.capacity() * size_of::<u64>()
    }

    /// Begin deciding `pass` of `passes`, dropping the keys the last one held and rewinding.
    fn begin_deciding(&mut self, pass: u32, passes: u32, ceiling: usize) -> io::Result<()> {
        // Each pass holds one slice of the digest space, so it needs one slice of the keys.
        let expected = self.hinted / passes.max(1) as usize;
        self.phase = Phase::Deciding { pass, ceiling };
        self.pass_shift = (u64::BITS - passes.max(1).trailing_zeros()) % u64::BITS;
        self.position = 0;
        self.index = Retention::with_capacity(expected)?;
        self.tape = KeyTape::with_capacity(expected, expected.saturating_mul(EXPECTED_KEY_BYTES));
        Ok(())
    }

    /// Begin writing what no pass marked, releasing the keys: the bits carry the answer now.
    fn begin_emitting(&mut self) -> io::Result<()> {
        self.phase = Phase::Emitting;
        self.position = 0;
        self.index = Retention::with_capacity(0)?;
        self.tape = KeyTape::default();
        Ok(())
    }

    /// Decide every record of `inputs` while holding no more than `ceiling` bytes.
    ///
    /// An attempt that crosses the ceiling is abandoned for one at twice the passes, and
    /// doubling renames every pass, so nothing an abandoned attempt marked may survive it.
    fn decide_every_record(
        &mut self,
        inputs: &[String],
        workers: &mut Workers,
        states: &mut [Scanned],
        drops: &mut DropBits,
        ceiling: usize,
    ) -> io::Result<()> {
        let mut passes = passes_for(ceiling, self.hinted);
        loop {
            match self.decide_in_passes(inputs, workers, states, drops, ceiling, passes) {
                Ok(()) => return self.begin_emitting(),
                Err(error) if !crossed_the_ceiling(&error) => return Err(error),
                Err(_) if passes >= MOST_PASSES => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--memory of {ceiling} bytes cannot hold one slice of the index \
                             even cut {MOST_PASSES} ways; raise the ceiling"
                        ),
                    ))
                }
                // Every mark and every half-scanned batch names the partition being abandoned.
                Err(_) => {
                    passes *= 2;
                    drops.clear();
                    states.iter_mut().for_each(Scanned::clear);
                }
            }
        }
    }

    /// One attempt: `passes` slices of the digest space decided in turn, marking duplicates.
    fn decide_in_passes(
        &mut self,
        inputs: &[String],
        workers: &mut Workers,
        states: &mut [Scanned],
        drops: &mut DropBits,
        ceiling: usize,
        passes: u32,
    ) -> io::Result<()> {
        // Deciding writes nothing, so the rendering this carries never reaches a byte.
        let mut discard =
            RecordWriter::with_rendering(io::sink(), Rendering::PLAIN, SequenceFormat::Fasta);
        self.held_peak = 0;
        for pass in 0..passes {
            self.begin_deciding(pass, passes, ceiling)?;
            for path in inputs {
                let mut input = InputFile::open(Some(path))?;
                for_each_record_in_input(&mut input, workers, states, Scanned::push, |state| {
                    self.absorb(state, &mut discard, drops)
                })?;
            }
        }
        Ok(())
    }

    /// Fold one scanned batch into the answer, writing whatever survives.
    ///
    /// Batches arrive in input order and are decided in position order, so first-occurrence
    /// means the same thing however many workers scanned.
    fn absorb<W: Write>(
        &mut self,
        scanned: &mut Scanned,
        writer: &mut RecordWriter<W>,
        drops: &mut DropBits,
    ) -> io::Result<()> {
        if let Some(format) = scanned.format {
            writer.adopt(format)?;
        }
        for batched in 0..scanned.len() {
            let tag = scanned.tags[batched];
            let (key, rendered) = scanned.record(batched);
            match self.phase {
                Phase::Immediate => {
                    self.examined += 1;
                    if self.index.retain(tag, key, &mut self.tape)? {
                        writer.write_rendered(rendered)?;
                        self.retained += 1;
                    }
                }
                // One slice of the digest space per pass; the rest belongs to another one.
                Phase::Deciding { pass, ceiling } => {
                    if self.pass_of(tag) == pass {
                        if self.index.retain(tag, key, &mut self.tape)? {
                            self.held_peak = self.held_peak.max(self.held_bytes());
                            if self.held_peak > ceiling {
                                return Err(io::Error::other(CeilingCrossed));
                            }
                        } else {
                            drops.mark(self.position);
                        }
                    }
                }
                Phase::Emitting => {
                    self.examined += 1;
                    if !drops.marked(self.position) {
                        writer.write_rendered(rendered)?;
                        self.retained += 1;
                    }
                }
            }
            self.position += 1;
        }
        scanned.clear();
        Ok(())
    }

    /// How many distinct keys are retained.
    #[cfg(test)]
    fn distinct(&self) -> usize {
        self.index.live
    }
}

// endregion: Deduplication

/// Collapse duplicate records
#[derive(Parser)]
#[command(name = "fasta-dedup")]
#[command(
    version,
    about = "Collapse duplicate records, keeping the first occurrence"
)]
#[command(
    long_about = "Remove duplicate records from FASTA or FASTQ input, keeping the first occurrence.\nIdentity is decided by --by, and is always settled by comparing bytes rather\nthan by a digest alone.\nMemory scales with the number of survivors, not with input size."
)]
struct Args {
    /// Input files; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// What decides whether two records are the same record
    #[arg(long, value_enum, value_name = "FIELD", default_value_t = By::Sequence, help_heading = "Identity")]
    by: By,

    /// Size the index for this many distinct records up front; it still grows past it
    #[arg(long, value_name = "N", help_heading = "Memory")]
    expect_distinct: Option<u32>,

    /// Hold no more than this much of the index at once, reading the inputs once more per
    /// halving; accepts a plain byte count or a K, M, or G suffix
    #[arg(long, value_name = "SIZE", help_heading = "Memory")]
    memory: Option<String>,

    /// Output file; '-' or omitted writes standard output
    #[arg(
        long,
        value_name = "FILE",
        conflicts_with_all = ["output_dir", "in_place", "dry_run", "quiet"],
        help_heading = "Output"
    )]
    output: Option<String>,

    /// Write one output per input into this directory, keeping each input's file name
    #[arg(
        long,
        value_name = "DIR",
        conflicts_with_all = ["in_place", "dry_run", "quiet"],
        help_heading = "Output"
    )]
    output_dir: Option<String>,

    /// Rewrite each input, swapping the result in once it is whole on disk
    #[arg(long, conflicts_with_all = ["dry_run", "quiet"], help_heading = "Output")]
    in_place: bool,

    /// Accept that an in-place run over several inputs leaves each one without the records an
    /// earlier input already held
    #[arg(long, requires = "in_place", help_heading = "Output")]
    allow_cross_input_loss: bool,

    /// Report what would be written, on standard error, without writing it
    #[arg(long, conflicts_with = "quiet", help_heading = "Output")]
    dry_run: bool,

    /// Suppress all output; the exit code carries the answer
    #[arg(long, help_heading = "Output")]
    quiet: bool,

    #[command(flatten)]
    presentation: Presentation,

    #[command(flatten)]
    parallelism: Parallelism,
}

/// Refuse rewriting several inputs in place, where corpus-wide identity is destructive.
///
/// A later input is rewritten without the records an earlier one already held, so no file is
/// self-contained afterwards and `--in-place` leaves no original to recover them from.
fn refuse_destructive_in_place(args: &Args) -> io::Result<()> {
    if !args.in_place || args.inputs.len() < 2 || args.allow_cross_input_loss {
        return Ok(());
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "--in-place over {} inputs would rewrite '{}' without the records '{}' held first, \
             because identity spans the whole corpus, and would replace the original with it; \
             write the results elsewhere with --output-dir, or pass --allow-cross-input-loss",
            args.inputs.len(),
            args.inputs[1],
            args.inputs[0],
        ),
    ))
}

fn run(args: &Args) -> io::Result<RunOutcome> {
    refuse_destructive_in_place(args)?;
    let identity = Identity::from(args.by);
    let expected = expected_distinct(args.expect_distinct);

    let destination = if args.quiet || args.dry_run {
        Destination::Discard
    } else if args.in_place {
        Destination::InPlace
    } else if let Some(directory) = &args.output_dir {
        Destination::Directory(directory.clone())
    } else {
        Destination::Stream(args.output.clone())
    };
    let rendering = args.presentation.rendering(destination.terminal_path());
    // Ordering must be preserved: the fold decides in input order, and that is the whole
    // reason first-occurrence survives going wide.
    let mut workers = args.parallelism.ordered()?;
    let mut states = workers.states(|| Scanned::new(identity, rendering));
    let mut combined = Deduplicator::new(expected)?;

    // One bit per record of the whole run, not per input: every pass reads the inputs in the
    // same order, so a record has the same position in all of them.
    let mut drops = DropBits::default();

    if let Some(asked) = &args.memory {
        let ceiling = memory_bytes(asked)?;
        // A bounded run reads its inputs once per pass, which a forward-only stream cannot do.
        every_input_is_rereadable(&args.inputs)?;
        combined.decide_every_record(
            &args.inputs,
            &mut workers,
            &mut states,
            &mut drops,
            ceiling,
        )?;
    }

    for_each_record_to_destination(
        &args.inputs,
        &destination,
        rendering,
        &mut workers,
        &mut states,
        Scanned::push,
        |state, writer| combined.absorb(state, writer, &mut drops),
    )?;

    if args.presentation.summary {
        eprintln!(
            "examined {} records, retained {}, dropped {}",
            combined.examined,
            combined.retained,
            combined.examined - combined.retained
        );
    }
    Ok(RunOutcome::of(combined.retained as usize))
}

fn main() {
    finish_or_exit("Error deduplicating", run(&Args::parse()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use fasterfasta::files::RandomAccess;

    /// A writer over a buffer the test can read back.
    fn collecting() -> RecordWriter<Vec<u8>> {
        RecordWriter::with_rendering(Vec::new(), Rendering::PLAIN, SequenceFormat::Fasta)
    }

    /// Run one buffer through the same scan-and-fold the tool uses, on one thread.
    fn dedup(data: &[u8], identity: Identity) -> String {
        dedup_in_batches(data, identity, usize::MAX)
    }

    /// The same buffer, folded one batch at a time, which is what many workers produce.
    fn dedup_in_batches(data: &[u8], identity: Identity, batch: usize) -> String {
        let mut scanned = Scanned::new(identity, Rendering::PLAIN);
        let mut combined = Deduplicator::new(0).unwrap();
        let mut writer = collecting();
        RandomAccess::new(data)
            .for_each_record(|record| {
                scanned.push(record)?;
                if scanned.len() >= batch {
                    combined.absorb(&mut scanned, &mut writer, &mut DropBits::default())?;
                }
                Ok(())
            })
            .unwrap();
        combined
            .absorb(&mut scanned, &mut writer, &mut DropBits::default())
            .unwrap();
        writer.flush().unwrap();
        String::from_utf8(writer.into_inner()).unwrap()
    }

    /// A corpus of `distinct` sequences, each repeated once, so half of it has to go.
    fn repetitive_corpus(distinct: u32) -> Vec<u8> {
        let mut data = Vec::new();
        for number in 0..distinct {
            for copy in 0..2 {
                data.extend_from_slice(
                    format!(">read{number}_{copy}\nACGTACGTAC{number:012}TTTTGGGGCC\n").as_bytes(),
                );
            }
        }
        data
    }

    /// What one unbounded pass over `corpus` holds, which is what a ceiling is measured against.
    fn held_unbounded(corpus: &[u8]) -> usize {
        let mut scanned = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        let mut combined = Deduplicator::new(0).unwrap();
        RandomAccess::new(corpus)
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        combined
            .absorb(&mut scanned, &mut collecting(), &mut DropBits::default())
            .unwrap();
        combined.held_bytes()
    }

    /// The answer a bounded run gives for `corpus`, with the peak bytes its decide phase held.
    /// A ceiling smaller than an empty table gives up instead of doubling forever.
    ///
    /// The smallest index is sixteen slots, so no number of passes fits it under a byte. The
    /// path is only reachable through the restart loop, which is why it is driven directly.
    #[test]
    fn a_ceiling_no_number_of_passes_can_hold_is_refused() {
        let identity = Identity::Field(RecordKey::Sequence);
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("corpus.fa");
        std::fs::write(&path, repetitive_corpus(64)).unwrap();
        let inputs = vec![path.to_string_lossy().into_owned()];

        let mut workers = Parallelism { threads: 1 }.ordered().unwrap();
        let mut states = workers.states(|| Scanned::new(identity, Rendering::PLAIN));
        let mut combined = Deduplicator::new(0).unwrap();
        let error = combined
            .decide_every_record(
                &inputs,
                &mut workers,
                &mut states,
                &mut DropBits::default(),
                1,
            )
            .expect_err("one byte cannot hold a slice of any index");
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("raise the ceiling"), "{error}");
    }

    fn dedup_under_ceiling(corpus: &[u8], ceiling: usize) -> (String, usize) {
        let identity = Identity::Field(RecordKey::Sequence);
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("corpus.fa");
        std::fs::write(&path, corpus).unwrap();
        let inputs = vec![path.to_string_lossy().into_owned()];

        let mut workers = Parallelism { threads: 1 }.ordered().unwrap();
        let mut states = workers.states(|| Scanned::new(identity, Rendering::PLAIN));
        let mut combined = Deduplicator::new(0).unwrap();
        let mut drops = DropBits::default();
        combined
            .decide_every_record(&inputs, &mut workers, &mut states, &mut drops, ceiling)
            .unwrap();
        let peak = combined.held_peak;

        let mut scanned = Scanned::new(identity, Rendering::PLAIN);
        let mut writer = collecting();
        RandomAccess::new(corpus)
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        combined
            .absorb(&mut scanned, &mut writer, &mut drops)
            .unwrap();
        writer.flush().unwrap();
        (String::from_utf8(writer.into_inner()).unwrap(), peak)
    }

    #[test]
    fn drops_repeated_sequences_keeping_the_first() {
        let data = b">a\nACGT\n>b\nTGCA\n>c\nACGT\n";
        assert_eq!(
            dedup(data, Identity::Field(RecordKey::Sequence)),
            ">a\nACGT\n>b\nTGCA\n"
        );
    }

    #[test]
    fn distinct_sequences_all_survive() {
        let data = b">a\nACGT\n>b\nTGCA\n>c\nGGGG\n";
        assert_eq!(
            dedup(data, Identity::Field(RecordKey::Sequence))
                .matches('>')
                .count(),
            3
        );
    }

    /// Wrapping is normalized before hashing, so layout cannot change identity.
    #[test]
    fn wrapping_does_not_affect_identity() {
        let data = b">a\nACGTTGCA\n>b\nACGT\nTGCA\n";
        assert_eq!(
            dedup(data, Identity::Field(RecordKey::Sequence)),
            ">a\nACGTTGCA\n"
        );
    }

    #[test]
    fn by_name_ignores_sequence_differences() {
        let data = b">gene1 first\nACGT\n>gene1 second\nTTTT\n>gene2\nGGGG\n";
        let out = dedup(data, Identity::Field(RecordKey::Identifier));
        assert!(out.contains("gene1 first"), "{out}");
        assert!(!out.contains("gene1 second"), "{out}");
        assert!(out.contains("gene2"), "{out}");
    }

    #[test]
    fn by_sequence_ignores_name_differences() {
        let data = b">first\nACGT\n>second\nACGT\n";
        assert_eq!(
            dedup(data, Identity::Field(RecordKey::Sequence)),
            ">first\nACGT\n"
        );
    }

    #[test]
    fn fastq_keeps_quality_of_the_first_occurrence() {
        let data = b"@a\nACGT\n+\nIIII\n@b\nACGT\n+\n####\n";
        assert_eq!(
            dedup(data, Identity::Field(RecordKey::Sequence)),
            "@a\nACGT\n+\nIIII\n"
        );
    }

    #[test]
    fn empty_input_produces_nothing() {
        assert_eq!(dedup(b"", Identity::Field(RecordKey::Sequence)), "");
    }

    #[test]
    fn counts_are_tracked() {
        let data = b">a\nACGT\n>b\nTGCA\n>c\nACGT\n";
        let mut scanned = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        let mut combined = Deduplicator::new(0).unwrap();
        RandomAccess::new(&data[..])
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        combined
            .absorb(&mut scanned, &mut collecting(), &mut DropBits::default())
            .unwrap();
        assert_eq!(combined.examined, 3);
        assert_eq!(combined.retained, 2);
    }

    /// Memory must follow the number of survivors, not the size of the input.
    #[test]
    fn state_scales_with_survivors_not_input() {
        let mut data = Vec::new();
        for _ in 0..10_000 {
            data.extend_from_slice(b">x\nACGT\n");
        }
        let mut scanned = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        let mut combined = Deduplicator::new(0).unwrap();
        let mut writer = collecting();
        RandomAccess::new(&data[..])
            .for_each_record(|record| {
                scanned.push(record)?;
                if scanned.len() >= 128 {
                    combined.absorb(&mut scanned, &mut writer, &mut DropBits::default())?;
                }
                Ok(())
            })
            .unwrap();
        combined
            .absorb(&mut scanned, &mut writer, &mut DropBits::default())
            .unwrap();
        assert_eq!(combined.examined, 10_000);
        assert_eq!(combined.distinct(), 1);
    }

    /// However the input is cut into batches, the answer is the same bytes.
    ///
    /// This is the property that lets the scan go wide: a batch boundary is where one worker's
    /// unit of work ends, so if the output moved with the boundary it would move with `-j`.
    #[test]
    fn batching_does_not_change_the_answer() {
        let mut data = Vec::new();
        for number in 0..500u32 {
            let sequence = match number % 7 {
                0 => "ACGTACGT",
                1 => "TTTTGGGG",
                2 => "CCCCAAAA",
                _ => "GATTACAG",
            };
            data.extend_from_slice(format!(">read{number}\n{sequence}\n").as_bytes());
        }
        let whole = dedup(&data, Identity::Field(RecordKey::Sequence));
        for batch in [1usize, 2, 3, 17, 64, 499, 1000] {
            assert_eq!(
                dedup_in_batches(&data, Identity::Field(RecordKey::Sequence), batch),
                whole,
                "batch size {batch} changed the answer"
            );
        }
    }

    /// The same property under canonicalization, where the key is derived rather than borrowed.
    #[test]
    fn batching_does_not_change_a_canonical_answer() {
        let mut data = Vec::new();
        for number in 0..300u32 {
            let sequence = match number % 4 {
                0 => "AAAACGT",
                1 => "ACGTTTT",
                2 => "GGGGCGT",
                _ => "ACGCCCC",
            };
            data.extend_from_slice(format!(">read{number}\n{sequence}\n").as_bytes());
        }
        let whole = dedup(&data, Identity::CanonicalStrand);
        for batch in [1usize, 5, 99] {
            assert_eq!(
                dedup_in_batches(&data, Identity::CanonicalStrand, batch),
                whole,
                "batch size {batch} changed the answer"
            );
        }
    }

    /// A sequence and its reverse complement are one molecule, so one of them survives.
    #[test]
    fn canonical_collapses_the_two_strands() {
        let data = b">forward\nAAAACGT\n>reverse\nACGTTTT\n";
        assert_eq!(
            dedup(data, Identity::CanonicalStrand),
            ">forward\nAAAACGT\n"
        );
        assert_eq!(
            dedup(data, Identity::Field(RecordKey::Sequence))
                .matches('>')
                .count(),
            2
        );
    }

    /// Canonicalization must not merge sequences that are merely similar.
    #[test]
    fn canonical_keeps_unrelated_sequences_apart() {
        let data = b">a\nAAAACGT\n>b\nGGGGCGT\n>c\nTTTTAAA\n";
        assert_eq!(
            dedup(data, Identity::CanonicalStrand).matches('>').count(),
            3
        );
    }

    /// The choice cannot depend on which strand was seen first.
    #[test]
    fn canonical_is_independent_of_input_order() {
        let forward_first = dedup(b">f\nAAAACGT\n>r\nACGTTTT\n", Identity::CanonicalStrand);
        let reverse_first = dedup(b">r\nACGTTTT\n>f\nAAAACGT\n", Identity::CanonicalStrand);
        assert_eq!(forward_first.matches('>').count(), 1);
        assert_eq!(reverse_first.matches('>').count(), 1);
        assert!(forward_first.contains(">f"), "{forward_first}");
        assert!(reverse_first.contains(">r"), "{reverse_first}");
    }

    /// A palindromic sequence is its own reverse complement and must not be dropped twice.
    #[test]
    fn canonical_handles_a_self_complementary_sequence() {
        let data = b">a\nACGT\n>b\nACGT\n>c\nAACGTT\n";
        let out = dedup(data, Identity::CanonicalStrand);
        assert_eq!(out.matches('>').count(), 2, "{out}");
    }

    /// A state outlives the input it first saw a format in, so clearing must forget it.
    ///
    /// The batch reports its format through `format`, but the scan writer latches one of its
    /// own through `adopt`. Resetting only the first left a FASTA input poisoning every FASTQ
    /// input after it, even with a separate output writer per input.
    #[test]
    fn a_cleared_batch_accepts_the_other_format() {
        let mut scanned = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        RandomAccess::new(&b">a\nACGT\n"[..])
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        scanned.clear();
        RandomAccess::new(&b"@b\nTGCA\n+\nIIII\n"[..])
            .for_each_record(|record| scanned.push(record))
            .expect("a cleared batch carries no format to disagree with");
        assert_eq!(scanned.format, Some(SequenceFormat::Fastq));
    }

    /// An RNA sequence is not a duplicate of the DNA spelling of itself.
    ///
    /// One shared complement table folding uracil onto adenine made `ACGU` and `ACGT` produce
    /// the same key, so a mixed corpus silently lost whichever came second.
    #[test]
    fn canonical_keeps_rna_apart_from_dna() {
        let data = b">dna\nACGT\n>rna\nACGU\n";
        let out = dedup(data, Identity::CanonicalStrand);
        assert_eq!(out.matches('>').count(), 2, "{out}");
    }

    /// Reverse complementing twice returns the original, on both alphabets.
    ///
    /// This is what a single table cannot give: adenine pairs with thymine in DNA and with
    /// uracil in RNA, so folding both into one table stops the operation being an involution.
    #[test]
    fn reverse_complement_is_an_involution_on_both_alphabets() {
        let complements = Complements::new();
        let mut scratch = Vec::new();
        for sequence in [b"AAAACGT".as_slice(), b"AAAACGU".as_slice()] {
            reverse_complement_into(sequence, &complements, &mut scratch);
            let once = scratch.clone();
            reverse_complement_into(&once, &complements, &mut scratch);
            assert_eq!(scratch, sequence, "not an involution for {sequence:?}");
        }
    }

    /// IUPAC ambiguity codes complement to their own base sets, in either case.
    #[test]
    fn canonical_covers_the_ambiguity_codes() {
        let complements = Complements::new();
        let mut scratch = Vec::new();
        // `S` and `W` are self-complementary, so they only change places, never letters.
        reverse_complement_into(b"RYKMBVDHNswacgt", &complements, &mut scratch);
        assert_eq!(scratch, b"acgtwsNDHBVKMRY");
        // Reverse complementing twice returns the original, which is what makes the lesser
        // of the two a stable choice.
        let once = scratch.clone();
        reverse_complement_into(&once, &complements, &mut scratch);
        assert_eq!(scratch, b"RYKMBVDHNswacgt");
    }

    /// An emptied batch holds no records, so it reports no format.
    ///
    /// Carrying the last format forward would let a run over a FASTA input and a FASTQ input
    /// agree where it should have refused.
    #[test]
    fn clearing_a_batch_forgets_its_format() {
        let mut scanned = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        RandomAccess::new(&b">a\nACGT\n"[..])
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        assert_eq!(scanned.format, Some(SequenceFormat::Fasta));
        scanned.clear();
        assert_eq!(scanned.format, None);
    }

    /// One output stream cannot carry both formats, however the batches fall.
    #[test]
    fn mixing_formats_across_batches_is_refused() {
        let mut combined = Deduplicator::new(0).unwrap();
        let mut writer = collecting();

        let mut fasta = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        RandomAccess::new(&b">a\nACGT\n"[..])
            .for_each_record(|record| fasta.push(record))
            .unwrap();
        combined
            .absorb(&mut fasta, &mut writer, &mut DropBits::default())
            .unwrap();

        let mut fastq = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        RandomAccess::new(&b"@b\nTGCA\n+\nIIII\n"[..])
            .for_each_record(|record| fastq.push(record))
            .unwrap();
        assert!(
            combined
                .absorb(&mut fastq, &mut writer, &mut DropBits::default())
                .is_err(),
            "a FASTQ batch after a FASTA one must be refused"
        );
    }

    /// Two distinct records sharing a digest must both survive.
    ///
    /// The digest is a reason to compare, never a reason to drop. A map keyed on the digest
    /// would overwrite one of these, and at a billion distinct records that is a silently
    /// dropped sequence roughly once every two billion — which is why this drives the table
    /// directly rather than through a hash nobody can force a collision in.
    #[test]
    fn a_shared_digest_does_not_collapse_distinct_keys() {
        let mut index = Retention::with_capacity(0).unwrap();
        let mut tape = KeyTape::default();
        assert!(index.retain(7, b"ACGT", &mut tape).unwrap());
        assert!(index.retain(7, b"TGCA", &mut tape).unwrap());
        assert_eq!(index.live, 2);

        // ...and each is still recognized as itself afterwards.
        assert!(!index.retain(7, b"ACGT", &mut tape).unwrap());
        assert!(!index.retain(7, b"TGCA", &mut tape).unwrap());
        assert_eq!(index.live, 2);
    }

    /// Keys survive the table doubling that a long run of insertions forces.
    #[test]
    fn growth_preserves_every_retained_key() {
        let mut index = Retention::with_capacity(0).unwrap();
        let mut tape = KeyTape::default();
        let keys: Vec<Vec<u8>> = (0..5_000u32)
            .map(|number| format!("seq{number}").into_bytes())
            .collect();
        for key in &keys {
            assert!(index
                .retain(hash_with_seed(key, SEED), key, &mut tape)
                .unwrap());
        }
        assert_eq!(index.live, keys.len());
        for key in &keys {
            assert!(!index
                .retain(hash_with_seed(key, SEED), key, &mut tape)
                .unwrap());
        }
        assert_eq!(index.live, keys.len());
    }

    /// Reserving room up front is an optimization, so it must retain exactly what growing does.
    #[test]
    fn pre_sizing_the_index_changes_nothing() {
        let keys: Vec<Vec<u8>> = (0..2_000u32)
            .map(|number| format!("seq{number}").into_bytes())
            .collect();
        let grown = Retention::with_capacity(0).unwrap();
        let reserved = Retention::with_capacity(keys.len()).unwrap();
        for mut index in [grown, reserved] {
            let mut tape = KeyTape::with_capacity(keys.len(), keys.len() * 8);
            for key in &keys {
                assert!(index
                    .retain(hash_with_seed(key, SEED), key, &mut tape)
                    .unwrap());
            }
            for key in &keys {
                assert!(!index
                    .retain(hash_with_seed(key, SEED), key, &mut tape)
                    .unwrap());
            }
            assert_eq!(index.live, keys.len());
        }
    }

    /// Sizes are read as bytes, and anything that is not a size is refused.
    #[test]
    fn a_memory_ceiling_is_read_as_bytes() {
        assert_eq!(memory_bytes("1024").unwrap(), 1024);
        assert_eq!(memory_bytes("1K").unwrap(), 1 << 10);
        assert_eq!(memory_bytes("512M").unwrap(), 512 << 20);
        assert_eq!(memory_bytes("48G").unwrap(), 48 << 30);
        assert_eq!(memory_bytes("48g").unwrap(), 48 << 30);
        for refused in ["", "0", "banana", "-1", "1.5G", "G"] {
            assert!(memory_bytes(refused).is_err(), "accepted {refused:?}");
        }
    }

    /// Halving the ceiling doubles the passes a run starts at, and the guess is never zero.
    #[test]
    fn a_tighter_ceiling_costs_more_passes() {
        let distinct = 1_000_000;
        let wide = passes_for(1 << 30, distinct);
        let narrow = passes_for(512 << 20, distinct);
        let narrower = passes_for(256 << 20, distinct);
        assert!(
            wide <= narrow && narrow <= narrower,
            "{wide} {narrow} {narrower}"
        );
        // Always a power of two, because a pass is named by the high bits of a digest.
        for passes in [wide, narrow, narrower, passes_for(usize::MAX, distinct)] {
            assert!(passes.is_power_of_two(), "{passes} is not a power of two");
        }
    }

    /// A ceiling far under what the keys need still holds, and the answer is unchanged.
    ///
    /// The pass count is chosen before reading and the distinct count is known only after, so
    /// the ceiling is kept by weighing what is held rather than by trusting either estimate.
    #[test]
    fn a_ceiling_bounds_what_a_run_holds() {
        let corpus = repetitive_corpus(8_000);
        let ceiling = held_unbounded(&corpus) / 16;
        let (bounded, peak) = dedup_under_ceiling(&corpus, ceiling);
        assert!(peak > 0 && peak <= ceiling, "held {peak} of {ceiling}");
        assert_eq!(
            bounded,
            dedup(&corpus, Identity::Field(RecordKey::Sequence))
        );
        assert_eq!(bounded.matches('>').count(), 8_000);
    }

    /// An abandoned attempt leaves no marks behind, so a restart answers what one try would.
    ///
    /// Doubling the passes renames every pass, so a mark made under the old partition names a
    /// record that a later attempt never decided.
    #[test]
    fn a_restart_answers_what_one_attempt_would() {
        let corpus = repetitive_corpus(4_000);
        let needed = held_unbounded(&corpus);
        let (restarted, tight) = dedup_under_ceiling(&corpus, needed / 32);
        let (settled, wide) = dedup_under_ceiling(&corpus, needed * 4);
        assert!(
            tight < wide,
            "the tight ceiling held {tight}, the wide {wide}"
        );
        assert_eq!(restarted, settled);
        assert_eq!(
            restarted,
            dedup(&corpus, Identity::Field(RecordKey::Sequence))
        );
    }

    /// Identity spans the corpus, so rewriting several inputs in place drops records from each.
    #[test]
    fn in_place_over_several_inputs_is_refused() {
        let several = Args::try_parse_from(["fasta-dedup", "a.fa", "b.fa", "--in-place"]).unwrap();
        let refusal = refuse_destructive_in_place(&several)
            .expect_err("b.fa would lose every record a.fa held first")
            .to_string();
        assert!(refusal.contains("--output-dir"), "{refusal}");
        assert!(refusal.contains("--allow-cross-input-loss"), "{refusal}");
        assert!(
            refusal.contains("a.fa") && refusal.contains("b.fa"),
            "{refusal}"
        );

        // Per-file and corpus-wide identity coincide for one input, so there is no hazard.
        for allowed in [
            vec!["fasta-dedup", "a.fa", "--in-place"],
            vec![
                "fasta-dedup",
                "a.fa",
                "b.fa",
                "--in-place",
                "--allow-cross-input-loss",
            ],
            vec!["fasta-dedup", "a.fa", "b.fa", "--output-dir", "clean"],
        ] {
            let args = Args::try_parse_from(&allowed).unwrap();
            assert!(
                refuse_destructive_in_place(&args).is_ok(),
                "refused {allowed:?}"
            );
        }
    }

    /// Every record belongs to exactly one pass, or a bounded run would drop or double it.
    #[test]
    fn every_digest_belongs_to_exactly_one_pass() {
        for passes in [1u32, 2, 4, 16] {
            let mut combined = Deduplicator::new(0).unwrap();
            combined.begin_deciding(0, passes, usize::MAX).unwrap();
            let mut seen = vec![0usize; passes as usize];
            for step in 0..10_000u64 {
                let tag = hash_with_seed(step.to_le_bytes(), SEED);
                let pass = combined.pass_of(tag);
                assert!(pass < passes, "{pass} is not one of {passes} passes");
                seen[pass as usize] += 1;
            }
            assert!(
                seen.iter().all(|count| *count > 0),
                "a pass got nothing: {seen:?}"
            );
        }
    }

    /// A bit is set for one record and no other.
    #[test]
    fn drop_bits_answer_for_the_record_they_were_set_for() {
        let mut drops = DropBits::default();
        for position in [0usize, 1, 63, 64, 65, 4095] {
            drops.mark(position);
        }
        for position in [0usize, 1, 63, 64, 65, 4095] {
            assert!(drops.marked(position), "{position} was marked");
        }
        for position in [2usize, 62, 66, 4094, 4096, 100_000] {
            assert!(!drops.marked(position), "{position} was never marked");
        }
    }

    /// Keys of unequal length must never compare equal, whatever the tape holds around them.
    #[test]
    fn a_prefix_is_not_the_key_it_prefixes() {
        let mut index = Retention::with_capacity(0).unwrap();
        let mut tape = KeyTape::default();
        assert!(index.retain(1, b"ACGT", &mut tape).unwrap());
        assert!(index.retain(1, b"ACG", &mut tape).unwrap());
        assert!(index.retain(1, b"ACGTA", &mut tape).unwrap());
        assert_eq!(index.live, 3);
    }

    /// Every flag is spelled out, so a call site says what it does and nothing is remembered
    /// by letter. `-h` and `-V` are clap's own and stay.
    #[test]
    fn declares_no_short_flags() {
        assert!(Args::command()
            .get_arguments()
            .all(|argument| argument.get_short().is_none()
                || matches!(argument.get_short(), Some('h') | Some('V'))));
    }

    /// Pins the surface, so adding, renaming, or reordering a flag is a deliberate edit here
    /// rather than a drift nobody reviews.
    #[test]
    fn declares_the_expected_flags() {
        let mut command = Args::command();
        command.build();
        let longs: Vec<_> = command
            .get_arguments()
            .filter_map(|argument| argument.get_long())
            .collect();
        assert_eq!(
            longs,
            [
                "by",
                "expect-distinct",
                "memory",
                "output",
                "output-dir",
                "in-place",
                "allow-cross-input-loss",
                "dry-run",
                "quiet",
                "line-width",
                "color",
                "summary",
                "threads",
                "help",
                "version"
            ]
        );
    }
}
