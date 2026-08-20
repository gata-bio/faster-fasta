//! Quality-based and fixed-position read trimming.
//!
//! Low-quality bases are trimmed from both ends, then fixed offsets are removed, then the
//! result is capped to a maximum length. Records shorter than a minimum after trimming are
//! dropped.
//!
//! __Memory__: O(1) — trimming resolves to a byte range, so nothing is copied. Resolving to a
//! range rather than to owned buffers is what makes a record the minimum-length filter rejects
//! cost nothing but arithmetic.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fastq-trim reads.fastq --trim-below-quality 20 --trim-end 5 --min-length 50 --output trimmed.fastq
//! fastq-trim reads.fastq --trim-start 10 --truncate-to 100 --output trimmed.fastq
//! cat reads.fastq | fastq-trim --trim-below-quality 20 > trimmed.fastq
//! ```
//!
//! Exit: 0 kept a record, 1 ran and kept none, 2 could not run.

use std::io::{self, Write};
use std::ops::Range;

use clap::Parser;

use stringzilla::sz::{find_byteset, rfind_byteset, Byteset};

use fasterfasta::files::{
    finish_or_exit, Destination, Presentation, RecordWriter, Rendering, RunOutcome,
};
use fasterfasta::records::{needs_fastq, phred33_to_ascii, Record, SequenceFormat};
use fasterfasta::scheduling::{for_each_record_to_destination, Parallelism};

#[derive(Debug, Clone, Copy, Default)]
struct TrimSettings {
    /// Quality bytes at or above the cutoff, or `None` to leave quality alone.
    ///
    /// The set rather than the cutoff, so there is no derived field to fall out of step with
    /// its source, and so the scan does not rebuild it per record.
    acceptable_scores: Option<Byteset>,
    trim_front: usize,
    trim_tail: usize,
    maximum_length: Option<usize>,
    minimum_length: Option<usize>,
}

/// Every byte encoding a Phred score at or above `cutoff`.
fn acceptable_scores(cutoff: u8) -> Byteset {
    let mut set = Byteset::new();
    for score in cutoff..=93 {
        set.add_u8(phred33_to_ascii(score));
    }
    set
}

impl TrimSettings {
    /// The surviving byte range, as an offset into both sequence and quality.
    ///
    /// Returning a range rather than owned buffers is the whole optimization: a record
    /// that the minimum-length filter will reject costs nothing but arithmetic.
    fn range(&self, record: &Record<'_>) -> Range<usize> {
        let mut start = 0usize;
        let mut end = record.sequence.len();

        // Trimming to the first and last acceptable base is what a byteset scan computes.
        // The bound is the quality length, not the sequence length: a FASTA record carries
        // no quality, and indexing it by a sequence offset would run off the end.
        if let Some(acceptable) = self.acceptable_scores {
            let scores = &record.quality[..end.min(record.quality.len())];
            start = find_byteset(scores, acceptable).unwrap_or(scores.len());
            end = rfind_byteset(scores, acceptable).map_or(start, |last| last + 1);
        }

        start = start.saturating_add(self.trim_front);
        end = end.saturating_sub(self.trim_tail);

        start = start.min(record.sequence.len());
        end = end.clamp(start, record.sequence.len());

        if let Some(maximum) = self.maximum_length {
            end = end.min(start.saturating_add(maximum));
        }

        start..end
    }

    fn keeps(&self, length: usize) -> bool {
        self.minimum_length.is_none_or(|minimum| length >= minimum)
    }
}

struct Trimmer<W: Write> {
    writer: RecordWriter<W>,
    settings: TrimSettings,
    examined: usize,
    retained: usize,
}

impl<W: Write> Trimmer<W> {
    fn new(writer: W, settings: TrimSettings, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(writer, rendering, SequenceFormat::Fastq),
            settings,
            examined: 0,
            retained: 0,
        }
    }

    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        if self.settings.acceptable_scores.is_some() {
            needs_fastq("--quality-cutoff", record.format())?;
        }
        self.writer.adopt(record.format())?;
        self.examined += 1;
        let range = self.settings.range(&record);
        if !self.settings.keeps(range.len()) {
            return Ok(());
        }

        // A FASTA record carries no quality, and the range is an offset into the sequence, so
        // a range the quality cannot hold yields nothing rather than an out-of-bounds slice.
        let trimmed_quality = record.quality.get(range.clone()).unwrap_or_default();

        self.writer.write_parts(
            record.header_without_sigil(),
            &record.sequence[range],
            trimmed_quality,
        )?;
        self.retained += 1;
        Ok(())
    }
}

/// Trim reads by quality and position
#[derive(Parser)]
#[command(name = "fastq-trim")]
#[command(version, about = "Trim reads by quality and fixed position")]
#[command(
    long_about = "Trim FASTQ reads from both ends.\nLow-quality bases go first, then fixed offsets, then a maximum-length cap.\nReads shorter than the minimum after trimming are dropped."
)]
struct Args {
    /// Input files; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// Trim bases below this Phred score from both ends
    #[arg(long, value_name = "PHRED", help_heading = "Trimming")]
    trim_below_quality: Option<u8>,

    /// Bases to remove from the start of every read
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = "Trimming")]
    trim_start: usize,

    /// Bases to remove from the end of every read
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = "Trimming")]
    trim_end: usize,

    /// Truncate reads longer than this, keeping the leading bases
    #[arg(long, value_name = "N", help_heading = "Trimming")]
    truncate_to: Option<usize>,

    /// Drop reads shorter than this after trimming
    #[arg(long, value_name = "N", help_heading = "Trimming")]
    min_length: Option<usize>,

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
    let settings = TrimSettings {
        acceptable_scores: args.trim_below_quality.map(acceptable_scores),
        trim_front: args.trim_start,
        trim_tail: args.trim_end,
        maximum_length: args.truncate_to,
        minimum_length: args.min_length,
    };

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
    let mut states = workers.states(|| Trimmer::new(Vec::new(), settings, rendering));

    for_each_record_to_destination(
        &args.inputs,
        &destination,
        rendering,
        &mut workers,
        &mut states,
        Trimmer::push,
        |state, writer| state.writer.drain_into(writer.inner_mut()),
    )?;

    let examined: usize = states.iter().map(|state| state.examined).sum();
    let retained: usize = states.iter().map(|state| state.retained).sum();
    if args.dry_run {
        eprintln!(
            "would trim {retained} of {examined} records, dropping {} below the minimum length; \
             nothing was written",
            examined - retained
        );
    } else if args.presentation.summary {
        eprintln!(
            "examined {examined} records, trimmed {retained}, dropped {} below the minimum length",
            examined - retained
        );
    }
    Ok(RunOutcome::of(retained))
}

fn main() {
    finish_or_exit("Error trimming reads", run(&Args::parse()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use fasterfasta::scheduling::for_each_record_in_bytes;

    fn trim(data: &[u8], settings: TrimSettings) -> String {
        let mut trimmer = Trimmer::new(Vec::new(), settings, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut trimmer, Trimmer::push).unwrap();
        String::from_utf8(trimmer.writer.into_inner()).unwrap()
    }

    fn trim_fasta(data: &[u8], settings: TrimSettings) -> io::Result<String> {
        let mut trimmer = Trimmer::new(Vec::new(), settings, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut trimmer, Trimmer::push)?;
        Ok(String::from_utf8(trimmer.writer.into_inner()).unwrap())
    }

    /// A FASTA record carries no quality, so a quality cutoff has nothing to act on.
    /// Bounding the walk by the sequence length instead indexed an empty slice and aborted.
    #[test]
    fn quality_cutoff_on_fasta_is_rejected_not_a_panic() {
        let settings = TrimSettings {
            acceptable_scores: Some(acceptable_scores(20)),
            ..Default::default()
        };
        let error = trim_fasta(b">seq1\nACGTACGT\n", settings).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("FASTQ"), "{error}");
    }

    /// Position trimming needs no quality, so it must still work on FASTA.
    #[test]
    fn position_trimming_works_on_fasta() {
        let settings = TrimSettings {
            trim_front: 2,
            ..Default::default()
        };
        assert_eq!(
            trim_fasta(b">seq1\nAACGTT\n", settings).unwrap(),
            ">seq1\nCGTT\n"
        );
    }

    #[test]
    fn trims_fixed_offsets() {
        let data = b"@seq1\nAACGTT\n+\nIIIIII\n";
        let front = TrimSettings {
            trim_front: 2,
            ..Default::default()
        };
        assert_eq!(trim(data, front), "@seq1\nCGTT\n+\nIIII\n");

        let tail = TrimSettings {
            trim_tail: 2,
            ..Default::default()
        };
        assert_eq!(trim(data, tail), "@seq1\nAACG\n+\nIIII\n");
    }

    /// Quality trimming walks in from both ends. '#' is Q2, 'I' is Q40.
    #[test]
    fn trims_low_quality_ends() {
        let data = b"@seq1\nAACGTT\n+\n##IIII\n";
        let settings = TrimSettings {
            acceptable_scores: Some(acceptable_scores(20)),
            ..Default::default()
        };
        assert_eq!(trim(data, settings), "@seq1\nCGTT\n+\nIIII\n");

        let both = b"@seq1\nAACGTT\n+\n#IIII#\n";
        assert_eq!(trim(both, settings), "@seq1\nACGT\n+\nIIII\n");
    }

    /// Every base below the cutoff leaves the scan with nothing to find.
    #[test]
    fn a_read_entirely_below_the_cutoff_is_emptied() {
        let settings = TrimSettings {
            acceptable_scores: Some(acceptable_scores(40)),
            ..Default::default()
        };
        assert_eq!(trim(b"@seq1\nACGT\n+\n####\n", settings), "@seq1\n\n+\n\n");
    }

    #[test]
    fn drops_reads_below_minimum_length() {
        let data = b"@keep\nACGTACGT\n+\nIIIIIIII\n@drop\nACGT\n+\nIIII\n";
        let settings = TrimSettings {
            minimum_length: Some(5),
            ..Default::default()
        };
        assert_eq!(trim(data, settings), "@keep\nACGTACGT\n+\nIIIIIIII\n");
    }

    #[test]
    fn caps_at_maximum_length() {
        let data = b"@seq1\nACGTACGT\n+\nIIIIIIII\n";
        let settings = TrimSettings {
            maximum_length: Some(4),
            ..Default::default()
        };
        assert_eq!(trim(data, settings), "@seq1\nACGT\n+\nIIII\n");
    }

    /// Trimming more than the read holds must yield an empty read, not underflow.
    #[test]
    fn over_trimming_yields_an_empty_read() {
        let data = b"@seq1\nACGT\n+\nIIII\n";
        let settings = TrimSettings {
            trim_front: 3,
            trim_tail: 3,
            ..Default::default()
        };
        assert_eq!(trim(data, settings), "@seq1\n\n+\n\n");
    }

    /// Quality and sequence must stay the same length after every combination.
    #[test]
    fn sequence_and_quality_stay_aligned() {
        let data = b"@seq1\nAACGTTAA\n+\n##IIII##\n";
        let settings = TrimSettings {
            acceptable_scores: Some(acceptable_scores(20)),
            trim_front: 1,
            trim_tail: 1,
            maximum_length: Some(2),
            ..Default::default()
        };
        let output = trim(data, settings);
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[1].len(), lines[3].len(), "{output}");
    }

    #[test]
    fn counts_are_tracked() {
        let data = b"@keep\nACGTACGT\n+\nIIIIIIII\n@drop\nACGT\n+\nIIII\n";
        let settings = TrimSettings {
            minimum_length: Some(5),
            ..Default::default()
        };
        let mut trimmer = Trimmer::new(Vec::new(), settings, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut trimmer, Trimmer::push).unwrap();
        assert_eq!(trimmer.examined, 2);
        assert_eq!(trimmer.retained, 1);
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
                "trim-below-quality",
                "trim-start",
                "trim-end",
                "truncate-to",
                "min-length",
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
