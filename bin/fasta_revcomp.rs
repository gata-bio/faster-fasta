//! Reverse complement of nucleotide sequences.
//!
//! Complements through a 256-byte lookup table and reverses in place. For FASTQ the
//! quality string is reversed alongside the sequence, since position _n_ from the start
//! becomes position _n_ from the end.
//!
//! __Memory__: O(longest record) — one scratch buffer reused across every record.
//! __Streaming__: yes, for files and pipes alike.
//!
//! # Examples
//!
//! ```bash
//! fasta-revcomp sequences.fasta -o revcomp.fasta
//! fasta-revcomp reads.fastq -o revcomp.fastq
//! cat sequences.fasta | fasta-revcomp | fasta-sort -l
//! ```

use std::io::{self, Write};

use clap::Parser;
use stringzilla::sz::lookup;

use fasterfasta::files::{finish_or_exit, open_output, RecordWriter, Rendering};
use fasterfasta::records::{Complements, Record, SequenceFormat};
use fasterfasta::scheduling::{for_each_record_in_inputs, Parallelism};

struct ReverseComplementer<W: Write> {
    writer: RecordWriter<W>,
    complements: Complements,
    sequence_scratch: Vec<u8>,
    quality_scratch: Vec<u8>,
    converted: u64,
}

impl<W: Write> ReverseComplementer<W> {
    fn new(writer: W, format: SequenceFormat, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(writer, rendering, format),
            complements: Complements::new(),
            sequence_scratch: Vec::new(),
            quality_scratch: Vec::new(),
            converted: 0,
        }
    }
}

/// Grows `buffer` only when it is too short, so the zero-fill happens once rather than
/// per record, and hands back exactly `length` bytes to overwrite.
fn scratch_of(buffer: &mut Vec<u8>, length: usize) -> &mut [u8] {
    if buffer.len() < length {
        buffer.resize(length, 0);
    }
    &mut buffer[..length]
}

impl<W: Write> ReverseComplementer<W> {
    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        self.writer.adopt(record.format())?;
        // Uracil standing in for thymine means RNA, where adenine must complement back to
        // uracil. Thymine anywhere settles it as DNA.
        let table = *self.complements.table_for(record.sequence);

        let complemented = scratch_of(&mut self.sequence_scratch, record.sequence.len());
        lookup(complemented, record.sequence, table);
        complemented.reverse();

        let reversed_quality = scratch_of(&mut self.quality_scratch, record.quality.len());
        reversed_quality.copy_from_slice(record.quality);
        reversed_quality.reverse();

        self.writer.write_parts(
            record.header_without_sigil(),
            &self.sequence_scratch[..record.sequence.len()],
            &self.quality_scratch[..record.quality.len()],
        )?;
        self.converted += 1;
        Ok(())
    }
}

/// Reverse complement nucleotide sequences
#[derive(Parser)]
#[command(name = "fasta-revcomp")]
#[command(version, about = "Reverse complement nucleotide sequences")]
#[command(
    long_about = "Reverse complement FASTA or FASTQ sequences over the full IUPAC alphabet.\nFor FASTQ the quality string is reversed to stay aligned with the sequence.\nSingle pass, works on files and pipes alike."
)]
struct Args {
    /// Input files, FASTA or FASTQ; '-' or omitted reads standard input
    #[arg(default_value = "-")]
    inputs: Vec<String>,

    /// Output file; '-' or omitted writes standard output
    #[arg(short, long)]
    output: Option<String>,

    /// Report how many records were converted, on standard error
    #[arg(long)]
    report: bool,

    #[command(flatten)]
    parallelism: Parallelism,
}

fn run(args: &Args) -> io::Result<()> {
    let rendering = Rendering::for_output(args.output.as_deref());
    let mut output = open_output(args.output.as_deref())?;
    let mut workers = args.parallelism.ordered()?;
    // Each state buffers the records of one unit of work, and `retire` writes that buffer
    // out in turn, so the output is the same bytes however many workers ran.
    let mut states =
        workers.states(|| ReverseComplementer::new(Vec::new(), SequenceFormat::Fasta, rendering));

    for_each_record_in_inputs(
        &args.inputs,
        &mut workers,
        &mut states,
        ReverseComplementer::push,
        |state| state.writer.drain_into(&mut output),
    )?;
    output.flush()?;

    if args.report {
        let converted: u64 = states.iter().map(|state| state.converted).sum();
        eprintln!("converted {converted} records");
    }
    Ok(())
}

fn main() {
    finish_or_exit("Error reverse complementing", run(&Args::parse()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use fasterfasta::scheduling::for_each_record_in_bytes;

    fn revcomp(data: &[u8], format: SequenceFormat) -> String {
        let mut complementer = ReverseComplementer::new(Vec::new(), format, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut complementer, ReverseComplementer::push).unwrap();
        String::from_utf8(complementer.writer.into_inner()).unwrap()
    }

    #[test]
    fn complements_and_reverses() {
        assert_eq!(
            revcomp(b">seq1\nACGT\n", SequenceFormat::Fasta),
            ">seq1\nACGT\n"
        );
        assert_eq!(
            revcomp(b">seq1\nAAAACGT\n", SequenceFormat::Fasta),
            ">seq1\nACGTTTT\n"
        );
    }

    #[test]
    fn preserves_case() {
        assert_eq!(
            revcomp(b">seq1\nacgtACGT\n", SequenceFormat::Fasta),
            ">seq1\nACGTacgt\n"
        );
    }

    #[test]
    fn handles_iupac_ambiguity_codes() {
        // R is A/G, so its complement is Y for T/C; S and W are self-complementary.
        assert_eq!(
            revcomp(b">seq1\nRYSWKM\n", SequenceFormat::Fasta),
            ">seq1\nKMWSRY\n"
        );
        assert_eq!(
            revcomp(b">seq1\nBDHVN\n", SequenceFormat::Fasta),
            ">seq1\nNBDHV\n"
        );
    }

    /// Quality runs backwards with the sequence, so a reversed read still pairs each base
    /// with its own score.
    #[test]
    fn reverses_quality_alongside_sequence() {
        assert_eq!(
            revcomp(b"@seq1\nACGT\n+\n!#$%\n", SequenceFormat::Fastq),
            "@seq1\nACGT\n+\n%$#!\n"
        );
    }

    /// Complementing twice must return the input for RNA as well as DNA. Mapping uracil to
    /// adenine while adenine mapped to thymine silently transcribed RNA into DNA.
    #[test]
    fn rna_round_trips_without_becoming_dna() {
        let original = ">seq1\nAUGCUUAA\n";
        let once = revcomp(original.as_bytes(), SequenceFormat::Fasta);
        let twice = revcomp(once.as_bytes(), SequenceFormat::Fasta);
        assert_eq!(twice, original, "one pass gave {once}");
        assert!(!once.contains('T'), "RNA gained a thymine: {once}");
    }

    #[test]
    fn round_trips_to_the_original() {
        let original = b">seq1\nACGTTGCANNRY\n";
        let once = revcomp(original, SequenceFormat::Fasta);
        let twice = revcomp(once.as_bytes(), SequenceFormat::Fasta);
        assert_eq!(twice, String::from_utf8(original.to_vec()).unwrap());
    }

    #[test]
    fn joins_wrapped_sequences() {
        assert_eq!(
            revcomp(b">seq1\nAAAA\nCGT\n", SequenceFormat::Fasta),
            ">seq1\nACGTTTT\n"
        );
    }

    #[test]
    fn scratch_is_reused_across_records() {
        let data = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n";
        let mut complementer =
            ReverseComplementer::new(Vec::new(), SequenceFormat::Fasta, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut complementer, ReverseComplementer::push).unwrap();
        // Sized to the longest record, never to the sum of them.
        assert_eq!(complementer.sequence_scratch.len(), 8);
    }

    #[test]
    fn counts_are_tracked() {
        let data = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n";
        let mut complementer =
            ReverseComplementer::new(Vec::new(), SequenceFormat::Fasta, Rendering::PLAIN);
        for_each_record_in_bytes(data, &mut complementer, ReverseComplementer::push).unwrap();
        assert_eq!(complementer.converted, 3);
    }
}
