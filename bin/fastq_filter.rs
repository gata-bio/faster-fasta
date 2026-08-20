//! Quality, length, and N-content filtering.
//!
//! Criteria combine with AND: a record is kept only if it satisfies every one supplied.
//!
//! __Memory__: O(1) — no per-record allocation.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fastq-filter reads.fastq --min-quality 25 --min-length 75 --output filtered.fastq
//! fastq-filter reads.fastq --max-n-fraction 0.05 --output low_ambiguity.fastq
//! cat reads.fastq | fastq-filter --min-quality 30 > filtered.fastq
//! ```
//!
//! Exit: 0 kept a record, 1 ran and kept none, 2 could not run.

use std::io::{self, Write};

use clap::Parser;

use fasterfasta::files::{
    finish_or_exit, Destination, Presentation, RecordWriter, Rendering, RunOutcome,
};
use fasterfasta::records::{mean_quality, needs_fastq, Record, SequenceFormat};
use fasterfasta::scheduling::{for_each_record_to_destination, Parallelism};

/// Every criterion is optional; a `None` never rejects anything.
#[derive(Debug, Clone, Copy, Default)]
struct Criteria {
    minimum_quality: Option<f32>,
    minimum_length: Option<usize>,
    maximum_length: Option<usize>,
    maximum_n_fraction: Option<f32>,
}

impl Criteria {
    /// Bounds that exclude each other would silently drop every record.
    fn validate(&self) -> io::Result<()> {
        let outside_unit_range = self
            .maximum_n_fraction
            .filter(|fraction| !(0.0..=1.0).contains(fraction));
        if let Some(fraction) = outside_unit_range {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("--max-n-fraction must be between 0.0 and 1.0, got {fraction}"),
            ));
        }
        let contradictory = self
            .minimum_length
            .zip(self.maximum_length)
            .filter(|(minimum, maximum)| minimum > maximum);
        if let Some((minimum, maximum)) = contradictory {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "--min-length {minimum} exceeds --max-length {maximum}, so no record can pass"
                ),
            ));
        }
        Ok(())
    }

    fn accepts(&self, record: &Record<'_>) -> bool {
        if self
            .minimum_quality
            .is_some_and(|threshold| mean_quality(record.quality) < threshold)
        {
            return false;
        }

        let length = record.sequence.len();
        if self.minimum_length.is_some_and(|minimum| length < minimum) {
            return false;
        }
        if self.maximum_length.is_some_and(|maximum| length > maximum) {
            return false;
        }

        // The divisor is floored at one, so an empty sequence divides no ambiguous bases by
        // one and clears every bound rather than dividing by zero.
        if self.maximum_n_fraction.is_some_and(|maximum| {
            let ambiguous = record
                .sequence
                .iter()
                .filter(|&&base| base == b'N' || base == b'n')
                .count();
            ambiguous as f32 / length.max(1) as f32 > maximum
        }) {
            return false;
        }

        true
    }
}

struct Filter<W: Write> {
    writer: RecordWriter<W>,
    criteria: Criteria,
    examined: usize,
    retained: usize,
}

impl<W: Write> Filter<W> {
    fn new(writer: W, criteria: Criteria, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(writer, rendering, SequenceFormat::Fastq),
            criteria,
            examined: 0,
            retained: 0,
        }
    }

    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        if self.criteria.minimum_quality.is_some() {
            needs_fastq("--min-quality", record.format())?;
        }
        self.writer.adopt(record.format())?;
        self.examined += 1;
        if !self.criteria.accepts(&record) {
            return Ok(());
        }
        self.writer.write_record(record)?;
        self.retained += 1;
        Ok(())
    }
}

/// Filter records by quality, length, and ambiguity
#[derive(Parser)]
#[command(name = "fastq-filter")]
#[command(version, about = "Filter records by quality, length, and N-content")]
#[command(
    long_about = "Filter FASTQ or FASTA records on mean quality, length bounds, and N-base fraction.\nCriteria combine with AND: a record is kept only if it satisfies every one given.\nSingle pass, constant memory, works on files and pipes alike."
)]
struct Args {
    /// Input files; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// Minimum mean Phred quality
    #[arg(long, value_name = "PHRED", help_heading = "Criteria")]
    min_quality: Option<f32>,

    /// Minimum sequence length
    #[arg(long, value_name = "N", help_heading = "Criteria")]
    min_length: Option<usize>,

    /// Maximum sequence length
    #[arg(long, value_name = "N", help_heading = "Criteria")]
    max_length: Option<usize>,

    /// Maximum fraction of N bases, between 0.0 and 1.0
    #[arg(long, value_name = "FRACTION", help_heading = "Criteria")]
    max_n_fraction: Option<f32>,

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
    let criteria = Criteria {
        minimum_quality: args.min_quality,
        minimum_length: args.min_length,
        maximum_length: args.max_length,
        maximum_n_fraction: args.max_n_fraction,
    };
    // Config is checked before any input is opened, so a contradiction fails at once.
    criteria.validate()?;

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
    let mut states = workers.states(|| Filter::new(Vec::new(), criteria, rendering));

    for_each_record_to_destination(
        &args.inputs,
        &destination,
        rendering,
        &mut workers,
        &mut states,
        Filter::push,
        |state, writer| state.writer.drain_into(writer.inner_mut()),
    )?;

    let examined: usize = states.iter().map(|state| state.examined).sum();
    let retained: usize = states.iter().map(|state| state.retained).sum();
    if args.dry_run {
        eprintln!(
            "would retain {retained} of {examined} records, dropping {}; nothing was written",
            examined - retained
        );
    } else if args.presentation.summary {
        eprintln!(
            "examined {examined} records, retained {retained}, dropped {}",
            examined - retained
        );
    }
    Ok(RunOutcome::of(retained))
}

fn main() {
    finish_or_exit("Error filtering records", run(&Args::parse()));
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
        let criteria = Criteria {
            minimum_quality: Some(30.0),
            ..Default::default()
        };
        let mut tool = Filter::new(Vec::new(), criteria, Rendering::PLAIN);
        let error =
            for_each_record_in_bytes(b">a\nACGT\n>b\nTTTT\n", &mut tool, Filter::push).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("--min-quality"), "{error}");
    }

    fn filter(data: &[u8], criteria: Criteria) -> String {
        let mut filter = Filter::new(Vec::new(), criteria, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut filter, Filter::push).unwrap();
        String::from_utf8(filter.writer.into_inner()).unwrap()
    }

    /// 'I' is Q40, '#' is Q2.
    const MIXED: &[u8] = b"@high\nACGT\n+\nIIII\n@low\nTGCA\n+\n####\n";

    #[test]
    fn filters_by_mean_quality() {
        let criteria = Criteria {
            minimum_quality: Some(30.0),
            ..Default::default()
        };
        assert_eq!(filter(MIXED, criteria), "@high\nACGT\n+\nIIII\n");
    }

    #[test]
    fn filters_by_length_bounds() {
        let data = b"@short\nAC\n+\nII\n@medium\nACGT\n+\nIIII\n@long\nACGTACGT\n+\nIIIIIIII\n";
        let minimum = Criteria {
            minimum_length: Some(4),
            ..Default::default()
        };
        assert!(!filter(data, minimum).contains("@short"));
        assert!(filter(data, minimum).contains("@medium"));

        let maximum = Criteria {
            maximum_length: Some(4),
            ..Default::default()
        };
        assert!(!filter(data, maximum).contains("@long"));
        assert!(filter(data, maximum).contains("@medium"));
    }

    #[test]
    fn filters_by_ambiguity_fraction() {
        // Half the bases are N in the second record.
        let data = b"@clean\nACGT\n+\nIIII\n@murky\nACNN\n+\nIIII\n";
        let criteria = Criteria {
            maximum_n_fraction: Some(0.25),
            ..Default::default()
        };
        let kept = filter(data, criteria);
        assert!(kept.contains("@clean"));
        assert!(!kept.contains("@murky"));
    }

    #[test]
    fn criteria_combine_with_and() {
        let data = b"@a\nACGT\n+\nIIII\n@b\nAC\n+\nII\n@c\nACGT\n+\n####\n";
        let criteria = Criteria {
            minimum_quality: Some(30.0),
            minimum_length: Some(4),
            ..Default::default()
        };
        // Only `a` satisfies both.
        assert_eq!(filter(data, criteria), "@a\nACGT\n+\nIIII\n");
    }

    #[test]
    fn no_criteria_keeps_everything() {
        assert_eq!(
            filter(MIXED, Criteria::default()),
            String::from_utf8(MIXED.to_vec()).unwrap()
        );
    }

    #[test]
    fn contradictory_length_bounds_are_rejected() {
        let criteria = Criteria {
            minimum_length: Some(100),
            maximum_length: Some(50),
            ..Default::default()
        };
        let error = criteria.validate().unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("no record can pass"), "{error}");
    }

    #[test]
    fn ambiguity_fraction_outside_the_unit_range_is_rejected() {
        let criteria = Criteria {
            maximum_n_fraction: Some(1.5),
            ..Default::default()
        };
        assert_eq!(
            criteria.validate().unwrap_err().kind(),
            io::ErrorKind::InvalidInput
        );
    }

    /// An empty sequence must not divide by zero when an ambiguity bound is set.
    #[test]
    fn empty_sequence_survives_ambiguity_check() {
        let criteria = Criteria {
            maximum_n_fraction: Some(0.1),
            ..Default::default()
        };
        assert_eq!(filter(b"@empty\n\n+\n\n", criteria), "@empty\n\n+\n\n");
    }

    #[test]
    fn counts_are_tracked() {
        let criteria = Criteria {
            minimum_quality: Some(30.0),
            ..Default::default()
        };
        let mut filter = Filter::new(Vec::new(), criteria, Rendering::PLAIN);
        for_each_record_in_bytes(MIXED, &mut filter, Filter::push).unwrap();
        assert_eq!(filter.examined, 2);
        assert_eq!(filter.retained, 1);
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
                "min-length",
                "max-length",
                "max-n-fraction",
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
