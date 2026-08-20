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
//! fasta-dedup sequences.fasta -o unique.fasta
//! fasta-dedup reads.fastq --by-name -o unique.fastq
//! fasta-dedup contigs.fa --canonical -j 16 -o unique.fa
//! cat sequences.fasta | fasta-dedup > unique.fasta
//! ```

use std::io::{self, Write};

use clap::Parser;
use stringzilla::sz::{equal, hash_with_seed, lookup};

use fasterfasta::files::{finish_or_exit, open_output, RecordWriter, Rendering};
use fasterfasta::records::{key_bytes, Complements, Record, RecordKey, SequenceFormat};
use fasterfasta::scheduling::{for_each_record_in_inputs, Parallelism};

// region: Identity

/// Seed for every digest this tool takes.
///
/// Fixed rather than random so a run is reproducible across machines and across passes.
const SEED: u64 = 0;

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
const EXPECTED_KEY_BYTES: usize = 64;

/// Distinct records to size the index for.
fn expected_distinct(asked: Option<u32>) -> usize {
    asked.unwrap_or(EXPECTED_DISTINCT) as usize
}

/// The global answer: which keys have been retained, and the output they were written to.
struct Deduplicator<W: Write> {
    writer: RecordWriter<W>,
    index: Retention,
    tape: KeyTape,
    examined: u64,
    retained: u64,
}

impl<W: Write> Deduplicator<W> {
    /// A fold sized for `expected` distinct records, which it still grows past if there are more.
    fn new(writer: W, rendering: Rendering, expected: usize) -> io::Result<Self> {
        Ok(Self {
            writer: RecordWriter::with_rendering(writer, rendering, SequenceFormat::Fasta),
            index: Retention::with_capacity(expected)?,
            tape: KeyTape::with_capacity(expected, expected.saturating_mul(EXPECTED_KEY_BYTES)),
            examined: 0,
            retained: 0,
        })
    }

    /// Fold one scanned batch into the answer, writing whatever survives.
    ///
    /// Batches arrive in input order and are decided in position order, so first-occurrence
    /// means the same thing however many workers scanned.
    fn absorb(&mut self, scanned: &mut Scanned) -> io::Result<()> {
        if let Some(format) = scanned.format {
            self.writer.adopt(format)?;
        }
        self.examined += scanned.len() as u64;
        for position in 0..scanned.len() {
            let tag = scanned.tags[position];
            let (key, rendered) = scanned.record(position);
            if self.index.retain(tag, key, &mut self.tape)? {
                self.writer.write_rendered(rendered)?;
                self.retained += 1;
            }
        }
        scanned.clear();
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        self.writer.flush()
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
    long_about = "Remove duplicate records from FASTA or FASTQ input, keeping the first occurrence.\nIdentity is decided by sequence content, by identifier with --by-name, or by the\nlesser of a sequence and its reverse complement with --canonical, and is always\nsettled by comparing bytes rather than by a digest alone.\nMemory scales with the number of survivors, not with input size."
)]
struct Args {
    /// Input files; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// Output file; '-' or omitted writes standard output
    #[arg(short, long)]
    output: Option<String>,

    /// Decide identity by identifier rather than by sequence content
    #[arg(long, conflicts_with = "canonical")]
    by_name: bool,

    /// Treat a sequence and its reverse complement as the same record
    #[arg(long)]
    canonical: bool,

    /// Size the index for this many distinct records up front; it still grows past it
    #[arg(long, value_name = "N")]
    expect_distinct: Option<u32>,

    /// Report how many records were examined and retained, on standard error
    #[arg(long)]
    report: bool,

    #[command(flatten)]
    parallelism: Parallelism,
}

fn run(args: &Args) -> io::Result<()> {
    let identity = match (args.by_name, args.canonical) {
        (true, _) => Identity::Field(RecordKey::Identifier),
        (_, true) => Identity::CanonicalStrand,
        _ => Identity::Field(RecordKey::Sequence),
    };

    let rendering = Rendering::for_output(args.output.as_deref());
    let output = open_output(args.output.as_deref())?;
    // Ordering must be preserved: the fold decides in input order, and that is the whole
    // reason first-occurrence survives going wide.
    let mut workers = args.parallelism.ordered()?;
    let mut states = workers.states(|| Scanned::new(identity, rendering));
    let expected = expected_distinct(args.expect_distinct);
    let mut combined = Deduplicator::new(output, rendering, expected)?;

    for_each_record_in_inputs(
        &args.inputs,
        &mut workers,
        &mut states,
        Scanned::push,
        |state| combined.absorb(state),
    )?;
    combined.finish()?;

    if args.report {
        eprintln!(
            "examined {} records, retained {}, dropped {}",
            combined.examined,
            combined.retained,
            combined.examined - combined.retained
        );
    }
    Ok(())
}

fn main() {
    finish_or_exit("Error deduplicating", run(&Args::parse()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use fasterfasta::files::RandomAccess;

    /// Run one buffer through the same scan-and-fold the tool uses, on one thread.
    fn dedup(data: &[u8], identity: Identity) -> String {
        let mut scanned = Scanned::new(identity, Rendering::PLAIN);
        let mut combined = Deduplicator::new(Vec::new(), Rendering::PLAIN, 0).unwrap();
        RandomAccess::new(data)
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        combined.absorb(&mut scanned).unwrap();
        combined.finish().unwrap();
        String::from_utf8(combined.writer.into_inner()).unwrap()
    }

    /// The same buffer, folded one batch at a time, which is what many workers produce.
    fn dedup_in_batches(data: &[u8], identity: Identity, batch: usize) -> String {
        let mut scanned = Scanned::new(identity, Rendering::PLAIN);
        let mut combined = Deduplicator::new(Vec::new(), Rendering::PLAIN, 0).unwrap();
        RandomAccess::new(data)
            .for_each_record(|record| {
                scanned.push(record)?;
                if scanned.len() >= batch {
                    combined.absorb(&mut scanned)?;
                }
                Ok(())
            })
            .unwrap();
        combined.absorb(&mut scanned).unwrap();
        combined.finish().unwrap();
        String::from_utf8(combined.writer.into_inner()).unwrap()
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
        let mut combined = Deduplicator::new(Vec::new(), Rendering::PLAIN, 0).unwrap();
        RandomAccess::new(&data[..])
            .for_each_record(|record| scanned.push(record))
            .unwrap();
        combined.absorb(&mut scanned).unwrap();
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
        let mut combined = Deduplicator::new(Vec::new(), Rendering::PLAIN, 0).unwrap();
        RandomAccess::new(&data[..])
            .for_each_record(|record| {
                scanned.push(record)?;
                if scanned.len() >= 128 {
                    combined.absorb(&mut scanned)?;
                }
                Ok(())
            })
            .unwrap();
        combined.absorb(&mut scanned).unwrap();
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
        let mut combined = Deduplicator::new(Vec::new(), Rendering::PLAIN, 0).unwrap();

        let mut fasta = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        RandomAccess::new(&b">a\nACGT\n"[..])
            .for_each_record(|record| fasta.push(record))
            .unwrap();
        combined.absorb(&mut fasta).unwrap();

        let mut fastq = Scanned::new(Identity::Field(RecordKey::Sequence), Rendering::PLAIN);
        RandomAccess::new(&b"@b\nTGCA\n+\nIIII\n"[..])
            .for_each_record(|record| fastq.push(record))
            .unwrap();
        assert!(
            combined.absorb(&mut fastq).is_err(),
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
}
