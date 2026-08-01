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
//! fasta-dna2rna sequences.fasta -o transcribed.fasta
//! cat reads.fastq | fasta-dna2rna > transcribed.fastq
//! ```

use std::io::{self, Write};

use clap::Parser;
use stringzilla::sz::lookup;

use faster_fasta::{finish_or_exit, map_reduce, open_output, Record, Rendering, RecordWriter, SequenceFormat};

/// Thymine to uracil in both cases; every other byte maps to itself.
fn transcription_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    for (index, entry) in table.iter_mut().enumerate() {
        *entry = index as u8;
    }
    table[b'T' as usize] = b'U';
    table[b't' as usize] = b'u';
    table
}

struct Transcriber<W: Write> {
    writer: RecordWriter<W>,
    table: [u8; 256],
    scratch: Vec<u8>,
}

impl<W: Write> Transcriber<W> {
    fn new(writer: W, format: SequenceFormat, rendering: Rendering) -> Self {
        Self {
            writer: RecordWriter::with_rendering(writer, rendering, format),
            table: transcription_table(),
            scratch: Vec::new(),
        }
    }

    fn push(&mut self, record: Record<'_>) -> io::Result<()> {
        self.writer.adopt(record.format())?;
        let length = record.sequence.len();
        if self.scratch.len() < length {
            self.scratch.resize(length, 0);
        }
        lookup(&mut self.scratch[..length], record.sequence, self.table);

        self.writer
            .write_parts(record.header_without_sigil(), &self.scratch[..length], record.quality)
    }

    fn finish(&mut self) -> io::Result<()> {
        self.writer.flush()
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
    #[arg(short, long)]
    output: Option<String>,
}

fn main() {
    let args = Args::parse();

    let result = (|| {
        let rendering = Rendering::for_output(args.output.as_deref());
        let output = open_output(args.output.as_deref())?;
        // The placeholder format is replaced by `begin` before any record arrives.
        let mut transcriber = Transcriber::new(output, SequenceFormat::Fasta, rendering);
        map_reduce::for_each_record_at_paths(
            &args.inputs,
            &mut transcriber,
            Transcriber::push,
        )?;
        transcriber.finish()
    })();

    finish_or_exit("Error transcribing", result);
}

#[cfg(test)]
mod tests {
    use super::*;
    
    fn transcribe(data: &[u8], format: SequenceFormat) -> String {
        let mut transcriber = Transcriber::new(Vec::new(), format, Rendering::PLAIN);
        map_reduce::for_each_record_in_bytes(data, &mut transcriber, Transcriber::push).unwrap();
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
}
