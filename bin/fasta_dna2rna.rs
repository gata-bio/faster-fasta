//! DNA to RNA transcription.
//!
//! Rewrites thymine as uracil through a 256-byte lookup table, preserving case. Headers
//! and quality strings pass through untouched, since transcription changes neither.
//!
//! __Memory__: O(longest record) — one scratch buffer reused across every record.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fasta-dna2rna sequences.fasta --output transcribed.fasta
//! cat reads.fastq | fasta-dna2rna > transcribed.fastq
//! ```
//!
//! Exit: 0 transcribed a record, 1 ran and transcribed none, 2 could not run.

use std::io::{self, Write};

use clap::Parser;
use stringzilla::sz::lookup;

use fasterfasta::files::{
    finish_or_exit, Destination, Presentation, RecordWriter, Rendering, RunOutcome,
};
use fasterfasta::records::{Record, SequenceFormat};
use fasterfasta::scheduling::{for_each_record_to_destination, Parallelism};

/// Thymine to uracil in both cases; every other byte maps to itself.
fn transcription_table() -> [u8; 256] {
    let mut table: [u8; 256] = core::array::from_fn(|index| index as u8);
    table[b'T' as usize] = b'U';
    table[b't' as usize] = b'u';
    table
}

struct Transcriber<W: Write> {
    writer: RecordWriter<W>,
    table: [u8; 256],
    scratch: Vec<u8>,
    converted: u64,
}

impl<W: Write> Transcriber<W> {
    fn new(writer: W, format: SequenceFormat, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(writer, rendering, format),
            table: transcription_table(),
            scratch: Vec::new(),
            converted: 0,
        }
    }

    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        self.writer.adopt(record.format())?;
        let length = record.sequence.len();
        if self.scratch.len() < length {
            self.scratch.resize(length, 0);
        }
        lookup(&mut self.scratch[..length], record.sequence, self.table);

        self.writer.write_parts(
            record.header_without_sigil(),
            &self.scratch[..length],
            record.quality,
        )?;
        self.converted += 1;
        Ok(())
    }
}

/// Transcribe DNA to RNA
#[derive(Parser)]
#[command(name = "fasta-dna2rna")]
#[command(version, about = "Transcribe DNA to RNA by rewriting T as U")]
#[command(
    long_about = "Transcribe FASTA or FASTQ sequences from DNA to RNA.\nThymine becomes uracil in both cases; headers and quality strings are unchanged.\nSingle pass, works on files and pipes alike."
)]
struct Args {
    /// Input files, FASTA or FASTQ; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

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
        workers.states(|| Transcriber::new(Vec::new(), SequenceFormat::Fasta, rendering));

    for_each_record_to_destination(
        &args.inputs,
        &destination,
        rendering,
        &mut workers,
        &mut states,
        Transcriber::push,
        |state, writer| state.writer.drain_into(writer.inner_mut()),
    )?;

    let converted: u64 = states.iter().map(|state| state.converted).sum();
    if args.dry_run {
        eprintln!("would convert {converted} records; nothing was written");
    } else if args.presentation.summary {
        eprintln!("converted {converted} records");
    }
    Ok(RunOutcome::of(converted as usize))
}

fn main() {
    finish_or_exit("Error transcribing", run(&Args::parse()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use fasterfasta::scheduling::for_each_record_in_bytes;

    fn transcribe(data: &[u8], format: SequenceFormat) -> String {
        let mut transcriber = Transcriber::new(Vec::new(), format, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut transcriber, Transcriber::push).unwrap();
        String::from_utf8(transcriber.writer.into_inner()).unwrap()
    }

    #[test]
    fn rewrites_thymine() {
        assert_eq!(
            transcribe(b">seq1\nACGT\n", SequenceFormat::Fasta),
            ">seq1\nACGU\n"
        );
    }

    #[test]
    fn leaves_sequences_without_thymine_alone() {
        assert_eq!(
            transcribe(b">seq1\nACGACG\n", SequenceFormat::Fasta),
            ">seq1\nACGACG\n"
        );
    }

    #[test]
    fn preserves_case() {
        assert_eq!(
            transcribe(b">seq1\nacgtACGT\n", SequenceFormat::Fasta),
            ">seq1\nacguACGU\n"
        );
    }

    #[test]
    fn handles_multiple_records() {
        assert_eq!(
            transcribe(b">a\nTTTT\n>b\nAAAA\n>c\nATAT\n", SequenceFormat::Fasta),
            ">a\nUUUU\n>b\nAAAA\n>c\nAUAU\n"
        );
    }

    #[test]
    fn preserves_quality_for_fastq() {
        assert_eq!(
            transcribe(b"@seq1\nACGT\n+\nIIII\n", SequenceFormat::Fastq),
            "@seq1\nACGT\n+\nIIII\n".replace("ACGT", "ACGU")
        );
    }

    #[test]
    fn joins_wrapped_sequences() {
        assert_eq!(
            transcribe(b">seq1\nAT\nGT\n", SequenceFormat::Fasta),
            ">seq1\nAUGU\n"
        );
    }

    #[test]
    fn counts_are_tracked() {
        let data = b">a\nTTTT\n>b\nAAAA\n>c\nATAT\n";
        let mut transcriber = Transcriber::new(Vec::new(), SequenceFormat::Fasta, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut transcriber, Transcriber::push).unwrap();
        assert_eq!(transcriber.converted, 3);
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
