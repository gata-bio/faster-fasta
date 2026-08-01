//! Random subsampling of sequence records.
//!
//! Two modes with different memory profiles. `--fraction` decides each record independently
//! and holds nothing; `--count` keeps a uniform reservoir of exactly N records.
//!
//! __Memory__: constant for `--fraction`, N records for `--count`.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fasta-sample sequences.fasta --count 1000 -o sample.fasta
//! fasta-sample reads.fastq --fraction 0.1 -o tenth.fastq
//! cat reads.fastq | fasta-sample --fraction 0.01 --seed 42 > sample.fastq
//! ```

use std::collections::BinaryHeap;
use std::io::{self, Write};

use clap::Parser;

use faster_fasta::{finish_or_exit, map_reduce, digest, open_output, Record, RecordKey, Rendering, RecordWriter, SequenceFormat};

/// Whether a record with this key survives a `--fraction` pass.
fn keeps(key: u64, probability: f64) -> bool {
    if probability <= 0.0 {
        return false;
    }
    if probability >= 1.0 {
        return true;
    }
    (key as f64) < probability * (u64::MAX as f64)
}

enum Mode {
    /// Keep each record independently with the given probability.
    Fraction(f64),
    /// Keep a uniform sample of exactly this many records.
    Count(usize),
}

/// One record kept by `Mode::Count`, ordered by digest so the largest is discarded first.
///
/// `arrival` is carried so the kept records can be emitted in input order rather than in
/// whatever order the heap holds them.
#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct Kept {
    /// Leads the struct so the derived ordering is by key, which is what the heap needs.
    key: u64,
    arrival: u64,
    rendered: Vec<u8>,
}

struct Sampler<W: Write> {
    writer: RecordWriter<W>,
    mode: Mode,
    seed: u64,
    /// Records seen so far, used only to stamp arrival order.
    seen: u64,
    /// Records held by `Mode::Count`; always empty under `Mode::Fraction`.
    kept: BinaryHeap<Kept>,
}

impl<W: Write> Sampler<W> {
    fn new(writer: W, mode: Mode, seed: u64, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(writer, rendering, SequenceFormat::Fasta),
            mode,
            seed,
            seen: 0,
            // Grown as records are kept, so `--count` cannot ask for an allocation up front.
            kept: BinaryHeap::new(),
        }
    }


    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        self.writer.adopt(record.format())?;
        // A record's fate is decided by its identity alone, so it does not depend on the
        // record's position, on how the inputs were ordered, or on how many workers ran.
        let key = digest(&record, RecordKey::Identifier, self.seed);

        match self.mode {
            Mode::Fraction(probability) => {
                if keeps(key, probability) {
                    self.writer.write_record(record)?;
                }
            }
            Mode::Count(count) => {
                // The `count` smallest keys form a uniform sample, and a max-heap discards
                // the current largest, so this merges across inputs by construction.
                let worst = self.kept.peek().map_or(u64::MAX, |top| top.key);
                if self.kept.len() < count || key < worst {
                    let mut rendered = Vec::new();
                    self.writer.render_into(record, &mut rendered);
                    self.kept.push(Kept {
                        key,
                        arrival: self.seen,
                        rendered,
                    });
                    if self.kept.len() > count {
                        self.kept.pop();
                    }
                }
            }
        }
        self.seen += 1;
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        // Selection is order-independent, but the output is not: records leave in the order
        // they arrived, so a coordinate-sorted input stays sorted.
        let mut kept: Vec<Kept> = std::mem::take(&mut self.kept).into_vec();
        kept.sort_unstable_by_key(|entry| entry.arrival);
        for entry in kept {
            self.writer.write_rendered(&entry.rendered)?;
        }
        self.writer.flush()
    }
}

/// Randomly subsample sequence records
#[derive(Parser)]
#[command(name = "fasta-sample")]
#[command(version, about = "Randomly subsample sequence records")]
#[command(
    long_about = "Subsample FASTA or FASTQ records.\n--fraction decides each record independently in constant memory.\n--count keeps a uniform reservoir of exactly N records.\nPass --seed to reproduce a sample exactly."
)]
struct Args {
    /// Input files; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// Output file; '-' or omitted writes standard output
    #[arg(short, long)]
    output: Option<String>,

    /// Keep exactly this many records, sampled uniformly
    #[arg(short = 'n', long, conflicts_with = "fraction")]
    count: Option<usize>,

    /// Keep each record with this probability, between 0.0 and 1.0
    #[arg(short = 'f', long, conflicts_with = "count")]
    fraction: Option<f64>,

    /// Seed for reproducible sampling
    #[arg(short = 's', long, default_value_t = 42)]
    seed: u64,
}

fn main() {
    let args = Args::parse();

    let mode = match (args.count, args.fraction) {
        (Some(0), _) => {
            eprintln!("Error: --count must be at least 1");
            std::process::exit(1);
        }
        (Some(count), _) => Mode::Count(count),
        (_, Some(fraction)) if (0.0..=1.0).contains(&fraction) => Mode::Fraction(fraction),
        (_, Some(fraction)) => {
            eprintln!("Error: --fraction must be between 0.0 and 1.0, got {fraction}");
            std::process::exit(1);
        }
        (None, None) => {
            eprintln!("Error: pass either --count or --fraction");
            std::process::exit(1);
        }
    };

    let result = (|| {
        let rendering = Rendering::for_output(args.output.as_deref());
        let output = open_output(args.output.as_deref())?;
        let mut sampler = Sampler::new(output, mode, args.seed, rendering);
        map_reduce::for_each_record_at_paths(
            &args.inputs,
            &mut sampler,
            Sampler::push,
        )?;
        sampler.finish()
    })();

    finish_or_exit("Error sampling records", result);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(data: &[u8], format: SequenceFormat, mode: Mode, seed: u64) -> String {
        let mut sampler = Sampler::new(Vec::new(), mode, seed, Rendering::PLAIN);
        map_reduce::for_each_record_in_bytes(data, &mut sampler, Sampler::push).unwrap();
        sampler.finish().unwrap();
        String::from_utf8(sampler.writer.into_inner()).unwrap()
    }

    fn fasta(records: usize) -> Vec<u8> {
        fasta_range(0, records)
    }

    fn fasta_range(first: usize, last: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in first..last {
            data.extend_from_slice(format!(">seq{index}\nACGT\n").as_bytes());
        }
        data
    }

    #[test]
    fn count_keeps_exactly_that_many() {
        let data = fasta(100);
        let out = sample(&data, SequenceFormat::Fasta, Mode::Count(10), 42);
        assert_eq!(out.matches('>').count(), 10);
    }

    #[test]
    fn count_larger_than_input_keeps_everything() {
        let data = fasta(5);
        let out = sample(&data, SequenceFormat::Fasta, Mode::Count(50), 42);
        assert_eq!(out.matches('>').count(), 5);
    }

    #[test]
    fn same_seed_reproduces_the_sample() {
        let data = fasta(100);
        let first = sample(&data, SequenceFormat::Fasta, Mode::Count(10), 7);
        let second = sample(&data, SequenceFormat::Fasta, Mode::Count(10), 7);
        assert_eq!(first, second);
    }

    #[test]
    fn different_seeds_give_different_samples() {
        let data = fasta(1000);
        let first = sample(&data, SequenceFormat::Fasta, Mode::Count(10), 1);
        let second = sample(&data, SequenceFormat::Fasta, Mode::Count(10), 2);
        assert_ne!(first, second);
    }

    #[test]
    fn fraction_zero_keeps_nothing_and_one_keeps_all() {
        let data = fasta(50);
        assert_eq!(sample(&data, SequenceFormat::Fasta, Mode::Fraction(0.0), 42), "");
        assert_eq!(
            sample(&data, SequenceFormat::Fasta, Mode::Fraction(1.0), 42)
                .matches('>')
                .count(),
            50
        );
    }

    #[test]
    fn fraction_lands_near_the_requested_share() {
        let data = fasta(10_000);
        let kept = sample(&data, SequenceFormat::Fasta, Mode::Fraction(0.1), 42)
            .matches('>')
            .count();
        assert!((900..1100).contains(&kept), "kept {kept} of 10000 at p=0.1");
    }

    /// Quality must survive the reservoir round trip alongside the sequence.
    #[test]
    fn fastq_records_keep_their_quality() {
        let data = b"@read0\nACGT\n+\nIIII\n@read1\nTGCA\n+\nHHHH\n";
        let out = sample(data, SequenceFormat::Fastq, Mode::Count(2), 42);
        assert!(out.contains("IIII"), "{out}");
        assert!(out.contains("HHHH"), "{out}");
    }

    /// The point of a seed: the same records selected however the input is divided.
    /// A position-driven RNG could not hold this, because every record's draw depended on
    /// how many records preceded it.
    #[test]
    fn selection_survives_splitting_the_input() {
        let head = fasta_range(0, 100);
        let tail = fasta_range(100, 200);
        let whole = [head.clone(), tail.clone()].concat();

        let from_whole = sample(&whole, SequenceFormat::Fasta, Mode::Fraction(0.2), 7);
        let mut from_halves = sample(&head, SequenceFormat::Fasta, Mode::Fraction(0.2), 7);
        from_halves.push_str(&sample(&tail, SequenceFormat::Fasta, Mode::Fraction(0.2), 7));

        assert_eq!(from_halves, from_whole);
    }

    /// Paired mates carry the same identifier, so both files select the same read pairs
    /// with no coordination beyond the shared seed.
    #[test]
    fn paired_mates_select_matching_pairs() {
        let mut first = Vec::new();
        let mut second = Vec::new();
        for index in 0..200 {
            first.extend_from_slice(format!("@read{index}/1\nACGT\n+\nIIII\n").as_bytes());
            second.extend_from_slice(format!("@read{index}/2\nTGCA\n+\nHHHH\n").as_bytes());
        }

        let strip = |text: String, suffix: &str| -> Vec<String> {
            text.lines()
                .filter(|line| line.starts_with('@'))
                .map(|line| line.trim_end_matches(suffix).to_string())
                .collect()
        };

        let kept_first = strip(sample(&first, SequenceFormat::Fastq, Mode::Fraction(0.25), 3), "/1");
        let kept_second = strip(sample(&second, SequenceFormat::Fastq, Mode::Fraction(0.25), 3), "/2");

        assert!(!kept_first.is_empty(), "sampled nothing");
        assert_eq!(kept_first, kept_second);
    }

    /// Selection is order-independent, but output is not: records leave as they arrived.
    #[test]
    fn output_preserves_input_order() {
        let data = fasta(100);
        let out = sample(&data, SequenceFormat::Fasta, Mode::Count(10), 42);
        let order: Vec<usize> = out
            .lines()
            .filter(|line| line.starts_with('>'))
            .map(|line| line.trim_start_matches(">seq").parse().unwrap())
            .collect();
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(order, sorted, "reservoir order leaked into the output");
    }

    /// A huge count must not ask for the allocation up front.
    #[test]
    fn absurd_count_does_not_preallocate() {
        let data = fasta(5);
        let out = sample(&data, SequenceFormat::Fasta, Mode::Count(1_000_000_000), 42);
        assert_eq!(out.matches('>').count(), 5);
    }

    #[test]
    fn empty_input_produces_nothing() {
        assert_eq!(sample(b"", SequenceFormat::Fasta, Mode::Count(10), 42), "");
        assert_eq!(sample(b"", SequenceFormat::Fasta, Mode::Fraction(0.5), 42), "");
    }
}
