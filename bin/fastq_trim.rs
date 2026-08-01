//! Quality-based and fixed-position read trimming.
//!
//! Low-quality bases are trimmed from both ends, then fixed offsets are removed, then the
//! result is capped to a maximum length. Records shorter than a minimum after trimming are
//! dropped.
//!
//! __Memory__: O(1) — trimming resolves to a byte range, so nothing is copied. The previous
//! implementation allocated two vectors per record, before the minimum-length filter could
//! reject the record and discard them.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fastq-trim reads.fastq -q 20 -t 5 -l 50 -o trimmed.fastq
//! fastq-trim reads.fastq -f 10 -L 100 -o trimmed.fastq
//! cat reads.fastq | fastq-trim -q 20 > trimmed.fastq
//! ```

use std::io::{self, Write};
use std::ops::Range;

use clap::Parser;

use stringzilla::sz::{find_byteset, rfind_byteset, Byteset};

use faster_fasta::{
    finish_or_exit, map_reduce, open_output, phred33_to_ascii, Record, RecordWriter, Rendering,
    SequenceFormat,
};

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
        !self.minimum_length.is_some_and(|minimum| length < minimum)
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
        if record.format() == SequenceFormat::Fasta && self.settings.acceptable_scores.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--quality-cutoff needs FASTQ input, but this input is FASTA",
            ));
        }
        self.writer.adopt(record.format())?;
        self.examined += 1;
        let range = self.settings.range(&record);
        if !self.settings.keeps(range.len()) {
            return Ok(());
        }

        // Quality is empty for FASTA, so it cannot be sliced by the same range.
        let trimmed_quality = if record.quality.is_empty() {
            &[][..]
        } else {
            &record.quality[range.clone()]
        };

        self.writer.write_parts(
            record.header_without_sigil(),
            &record.sequence[range],
            trimmed_quality,
        )?;
        self.retained += 1;
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        self.writer.flush()
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

    /// Output file; '-' or omitted writes standard output
    #[arg(short, long)]
    output: Option<String>,

    /// Trim bases below this Phred score from both ends
    #[arg(short = 'q', long)]
    quality_cutoff: Option<u8>,

    /// Bases to remove from the start of every read
    #[arg(short = 'f', long, default_value_t = 0)]
    trim_front: usize,

    /// Bases to remove from the end of every read
    #[arg(short = 't', long, default_value_t = 0)]
    trim_tail: usize,

    /// Truncate reads longer than this, keeping the leading bases
    #[arg(short = 'T', long)]
    truncate_to: Option<usize>,

    /// Drop reads shorter than this after trimming
    #[arg(short = 'l', long)]
    min_length: Option<usize>,
}

fn main() {
    let args = Args::parse();

    let settings = TrimSettings {
        acceptable_scores: args.quality_cutoff.map(acceptable_scores),
        trim_front: args.trim_front,
        trim_tail: args.trim_tail,
        maximum_length: args.truncate_to,
        minimum_length: args.min_length,
    };

    let result = (|| {
        let rendering = Rendering::for_output(args.output.as_deref());
        let output = open_output(args.output.as_deref())?;
        let mut trimmer = Trimmer::new(output, settings, rendering);
        map_reduce::for_each_record_at_paths(
            &args.inputs,
            &mut trimmer,
            Trimmer::push,
        )?;
        trimmer.finish()
    })();

    finish_or_exit("Error trimming reads", result);
}

#[cfg(test)]
mod tests {
    use super::*;
    
    fn trim(data: &[u8], settings: TrimSettings) -> String {
        let mut trimmer = Trimmer::new(Vec::new(), settings, Rendering::PLAIN);
        map_reduce::for_each_record_in_bytes(
            data,
            &mut trimmer,
            Trimmer::push).unwrap();
        String::from_utf8(trimmer.writer.into_inner()).unwrap()
    }

    fn trim_fasta(data: &[u8], settings: TrimSettings) -> io::Result<String> {
        let mut trimmer = Trimmer::new(Vec::new(), settings, Rendering::PLAIN);
        map_reduce::for_each_record_in_bytes(
            data,
            &mut trimmer,
            Trimmer::push,
        )?;
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
        map_reduce::for_each_record_in_bytes(
            data,
            &mut trimmer,
            Trimmer::push).unwrap();
        assert_eq!(trimmer.examined, 2);
        assert_eq!(trimmer.retained, 1);
    }
}
