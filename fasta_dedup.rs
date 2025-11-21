//! Sequence deduplication utility for FASTA and FASTQ files
//!
//! Remove duplicate sequences, keeping first occurrence.
//! Uses StringZilla's `hash()` for SIMD-accelerated hashing.
//! Auto-detects format and preserves it in output.
//!
//! **Memory**: O(n) - stores all unique sequences in HashSet
//! **Streaming**: No - materializes HashSet in memory
//!
//! # Examples
//!
//! ```bash
//! # FASTA files
//! fasta-dedup sequences.fasta -o unique.fasta
//!
//! # FASTQ files (preserves quality scores)
//! fasta-dedup reads.fastq -o unique.fastq
//!
//! # From stdin
//! cat sequences.fasta | fasta-dedup > unique.fasta
//! ```

use std::collections::HashSet;
use std::io::{self, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::BuildSzHasher;

mod shared;
use shared::*;

/// Deduplicate sequences (FASTA or FASTQ)
///
/// Removes duplicate sequences, keeping the first occurrence of each unique sequence.
/// Auto-detects format and preserves it in output.
/// Uses StringZilla for hashing and fast operations.
pub fn fasta_dedup(input: &[u8], mut output: impl Write) -> io::Result<()> {
    let format = detect_format(input)?;
    let mut seen_sequences: HashSet<SequenceData, BuildSzHasher> =
        HashSet::with_hasher(BuildSzHasher::default());

    match format {
        SeqFormat::Fasta => {
            let parser = FastaParser::new(input);
            for entry in parser {
                if seen_sequences.insert(entry.sequence.clone()) {
                    output.write_all(entry.header)?;
                    output.write_all(b"\n")?;
                    output.write_all(entry.sequence.as_bytes())?;
                    output.write_all(b"\n")?;
                }
            }
        }
        SeqFormat::Fastq => {
            let parser = FastqParser::new(input);
            for entry in parser {
                if seen_sequences.insert(entry.sequence.clone()) {
                    output.write_all(entry.header)?;
                    output.write_all(b"\n")?;
                    output.write_all(entry.sequence.as_bytes())?;
                    output.write_all(b"\n+\n")?;
                    output.write_all(entry.quality)?;
                    output.write_all(b"\n")?;
                }
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Remove duplicate sequences from FASTA or FASTQ files
#[derive(Parser)]
#[command(name = "fasta-dedup")]
#[command(
    version,
    about = "Remove duplicate sequences, keeping first occurrence"
)]
#[command(
    long_about = "Remove duplicate sequences from FASTA or FASTQ files.\nAuto-detects format and preserves it in output.\nUses SIMD-accelerated hashing for maximum performance."
)]
struct Args {
    /// Input file (FASTA or FASTQ, use '-' or omit for stdin)
    input: Option<String>,

    /// Output file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,
}

fn main() {
    let args = Args::parse();

    let input = match get_input(args.input.as_deref()) {
        Ok(input) => input,
        Err(e) => {
            eprintln!("Error reading input: {}", e);
            process::exit(1);
        }
    };

    let output = match get_output(args.output.as_deref()) {
        Ok(output) => output,
        Err(e) => {
            eprintln!("Error opening output: {}", e);
            process::exit(1);
        }
    };

    if let Err(e) = fasta_dedup(input.as_bytes(), output) {
        if e.kind() == io::ErrorKind::BrokenPipe {
            process::exit(0);
        }
        eprintln!("Error processing sequences: {}", e);
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_dedup() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nACGT\n";
        let mut output = Vec::new();
        fasta_dedup(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1"));
        assert!(result.contains(">seq2"));
        assert!(!result.contains(">seq3"));
    }

    #[test]
    fn test_dedup_with_whitespace_variations() {
        let data = b">seq1\nACGT\n>seq2\n  ACGT  \n>seq3\nA C G T\n";
        let mut output = Vec::new();
        fasta_dedup(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1"));
        assert!(!result.contains(">seq2"));
        assert!(!result.contains(">seq3"));
        assert_eq!(result.matches("ACGT").count(), 1);
    }

    #[test]
    fn test_fastq_dedup() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n@seq3\nACGT\n+\nJJJJ\n";
        let mut output = Vec::new();
        fasta_dedup(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("@seq1"));
        assert!(result.contains("@seq2"));
        assert!(!result.contains("@seq3")); // Duplicate sequence
        assert_eq!(result.matches('@').count(), 2); // Only 2 unique sequences
    }
}
