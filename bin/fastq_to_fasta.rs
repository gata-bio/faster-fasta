//! FASTQ to FASTA conversion.
//!
//! Drops quality scores and rewrites `@` headers as `>`, optionally filtering on mean
//! quality first — which is the only reason to do this before dropping the scores.
//!
//! __Memory__: O(1) — no per-record allocation.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fastq-to-fasta reads.fastq --output reads.fasta
//! fastq-to-fasta reads.fastq --min-quality 30 --output high_quality.fasta
//! cat reads.fastq | fastq-to-fasta > reads.fasta
//! ```
//!
//! Exit: 0 converted a record, 1 ran and converted none, 2 could not run.

use std::io::{self, Write};

use clap::Parser;

use fasterfasta::files::{
    finish_or_exit, Destination, Presentation, RecordWriter, Rendering, RunOutcome,
};
use fasterfasta::records::{mean_quality, needs_fastq, Record, SequenceFormat};
use fasterfasta::scheduling::{for_each_record_to_destination, Parallelism};

struct FastaConverter<W: Write> {
    writer: RecordWriter<W>,
    minimum_quality: Option<f32>,
    converted: usize,
    skipped: usize,
}

impl<W: Write> FastaConverter<W> {
    fn new(writer: W, minimum_quality: Option<f32>, rendering: Rendering) -> Self {
        Self {
            // The output format alone performs the conversion: the sigil is rewritten
            // and quality is not emitted.
            writer: RecordWriter::with_rendering(writer, rendering, SequenceFormat::Fasta),
            minimum_quality,
            converted: 0,
            skipped: 0,
        }
    }

    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        if self.minimum_quality.is_some() {
            needs_fastq("--min-quality", record.format())?;
        }
        if self
            .minimum_quality
            .is_some_and(|threshold| mean_quality(record.quality) < threshold)
        {
            self.skipped += 1;
            return Ok(());
        }
        self.writer.write_record(record)?;
        self.converted += 1;
        Ok(())
    }
}

/// Convert FASTQ to FASTA
#[derive(Parser)]
#[command(name = "fastq-to-fasta")]
#[command(version, about = "Convert FASTQ to FASTA format")]
#[command(
    long_about = "Convert FASTQ input to FASTA by dropping quality scores.\nOptionally filters on mean quality before the scores are discarded.\nSingle pass, constant memory, works on files and pipes alike."
)]
struct Args {
    /// Input FASTQ files; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// Keep only records whose mean quality is at least this Phred score
    #[arg(long, value_name = "PHRED", help_heading = "Criteria")]
    min_quality: Option<f32>,

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

fn run(args: &Args) -> io::Result<RunOutcome> {
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
    let mut workers = args.parallelism.ordered()?;
    // Each state buffers the records of one unit of work, and `retire` writes that buffer
    // out in turn, so the output is the same bytes however many workers ran.
    let mut states =
        workers.states(|| FastaConverter::new(Vec::new(), args.min_quality, rendering));

    for_each_record_to_destination(
        &args.inputs,
        &destination,
        rendering,
        &mut workers,
        &mut states,
        FastaConverter::push,
        |state, writer| state.writer.drain_into(writer.inner_mut()),
    )?;

    let converted: usize = states.iter().map(|state| state.converted).sum();
    let skipped: usize = states.iter().map(|state| state.skipped).sum();
    if args.dry_run {
        eprintln!("would convert {converted} records and skip {skipped}; nothing was written");
    } else if args.presentation.summary {
        eprintln!("converted {converted} records, skipped {skipped}");
    }
    Ok(RunOutcome::of(converted))
}

fn main() {
    finish_or_exit("Error converting FASTQ", run(&Args::parse()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use fasterfasta::scheduling::for_each_record_in_bytes;

    /// A FASTA record has no quality, so a mean-quality threshold would silently reject
    /// every record and exit successfully. It is an error instead.
    #[test]
    fn quality_threshold_on_fasta_is_rejected() {
        let mut tool = FastaConverter::new(Vec::new(), Some(30.0), Rendering::PLAIN);
        let error =
            for_each_record_in_bytes(b">a\nACGT\n>b\nTTTT\n", &mut tool, FastaConverter::push)
                .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("--min-quality"), "{error}");
    }

    fn convert(data: &[u8], minimum_quality: Option<f32>) -> String {
        let mut converter = FastaConverter::new(Vec::new(), minimum_quality, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut converter, FastaConverter::push).unwrap();
        String::from_utf8(converter.writer.into_inner()).unwrap()
    }

    #[test]
    fn converts_headers_and_drops_quality() {
        let data = b"@seq1 description\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n";
        assert_eq!(
            convert(data, None),
            ">seq1 description\nACGT\n>seq2\nTGCA\n"
        );
    }

    #[test]
    fn filters_on_mean_quality() {
        // 'I' is Q40, '#' is Q2.
        let data = b"@high\nACGT\n+\nIIII\n@low\nTGCA\n+\n####\n";
        assert_eq!(convert(data, Some(30.0)), ">high\nACGT\n");
    }

    /// Multi-line input is joined before conversion, so a wrapped record converts to the same
    /// bytes an unwrapped one does and the quality filter sees one contiguous string.
    #[test]
    fn joins_wrapped_records() {
        let data = b"@seq1\nACGT\nTGCA\n+\nIIII\nHHHH\n";
        assert_eq!(convert(data, None), ">seq1\nACGTTGCA\n");
    }

    #[test]
    fn empty_input_produces_nothing() {
        assert_eq!(convert(b"", None), "");
    }

    #[test]
    fn counts_are_tracked() {
        let data = b"@high\nACGT\n+\nIIII\n@low\nTGCA\n+\n####\n";
        let mut converter = FastaConverter::new(Vec::new(), Some(30.0), Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut converter, FastaConverter::push).unwrap();
        assert_eq!(converter.converted, 1);
        assert_eq!(converter.skipped, 1);
    }

    /// Every flag is spelled out, so a call site says what it does and nothing is remembered
    /// by letter. `-h` and `-V` are clap's own and stay.
    #[test]
    fn declares_no_short_flags() {
        // Built first, because `-h` and `-V` are only added then and they are the exemption.
        let mut command = Args::command();
        command.build();
        assert!(command
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
                "min-quality",
                "output",
                "output-dir",
                "in-place",
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
