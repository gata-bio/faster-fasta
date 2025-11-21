//! DNA to RNA conversion utility for FASTA and FASTQ files
//!
//! Convert DNA sequences to RNA (T→U).
//! Uses StringZilla's `lookup()` for SIMD-accelerated character mapping.
//! Auto-detects format and preserves it (FASTQ quality scores are preserved).
//!
//! **Memory**: O(1) - processes one sequence at a time
//! **Streaming**: Yes - constant memory usage
//!
//! # Examples
//!
//! ```bash
//! # Convert FASTA
//! fasta-dna2rna sequences.fasta -o rna.fasta
//!
//! # Convert FASTQ (quality preserved)
//! fasta-dna2rna reads.fastq -o rna.fastq
//!
//! # From stdin
//! cat dna.fasta | fasta-dna2rna > rna.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::lookup;

mod shared;
use shared::*;

/// Convert DNA sequences to RNA (FASTA or FASTQ)
///
/// Converts T (thymine) to U (uracil) using StringZilla's lookup() function.
/// Handles both uppercase and lowercase bases.
/// Auto-detects format and preserves it (FASTQ quality scores are preserved).
pub fn fasta_dna2rna(input: &[u8], mut output: impl Write) -> io::Result<()> {
    let format = detect_format(input)?;

    // StringZilla translation table for DNA -> RNA
    let mut table = [0u8; 256];
    for i in 0..256 {
        table[i] = i as u8;
    }
    // Convert T to U
    table[b'T' as usize] = b'U';
    table[b't' as usize] = b'u';

    match format {
        SeqFormat::Fasta => {
            let parser = FastaParser::new(input);
            for entry in parser {
                let seq = entry.sequence.as_bytes();
                let mut rna = vec![0u8; seq.len()];
                lookup(&mut rna, seq, table);

                output.write_all(entry.header)?;
                output.write_all(b"\n")?;
                output.write_all(&rna)?;
                output.write_all(b"\n")?;
            }
        }
        SeqFormat::Fastq => {
            let parser = FastqParser::new(input);
            for entry in parser {
                let seq = entry.sequence.as_bytes();
                let mut rna = vec![0u8; seq.len()];
                lookup(&mut rna, seq, table);

                output.write_all(entry.header)?;
                output.write_all(b"\n")?;
                output.write_all(&rna)?;
                output.write_all(b"\n+\n")?;
                output.write_all(entry.quality)?;
                output.write_all(b"\n")?;
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Convert DNA sequences to RNA (T -> U)
#[derive(Parser)]
#[command(name = "fasta-dna2rna")]
#[command(version, about = "Convert DNA to RNA (T→U)")]
#[command(
    long_about = "Convert DNA sequences to RNA from FASTA or FASTQ files.\nAuto-detects format and preserves it (FASTQ quality scores are preserved).\nUses SIMD-accelerated character mapping."
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

    if let Err(e) = fasta_dna2rna(input.as_bytes(), output) {
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
    fn test_dna2rna_simple() {
        let data = b">seq1\nACGT\n";
        let mut output = Vec::new();
        fasta_dna2rna(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1"));
        assert!(result.contains("ACGU"));
    }

    #[test]
    fn test_dna2rna_no_thymine() {
        let data = b">seq1\nACGACG\n";
        let mut output = Vec::new();
        fasta_dna2rna(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("ACGACG"));
    }

    #[test]
    fn test_dna2rna_all_thymine() {
        let data = b">seq1\nTTTT\n";
        let mut output = Vec::new();
        fasta_dna2rna(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("UUUU"));
    }

    #[test]
    fn test_dna2rna_mixed_case() {
        let data = b">seq1\nAcGt\n";
        let mut output = Vec::new();
        fasta_dna2rna(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("AcGu"));
    }

    #[test]
    fn test_dna2rna_multiple() {
        let data = b">seq1\nATGC\n>seq2\nTTAA\n";
        let mut output = Vec::new();
        fasta_dna2rna(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[1], "AUGC");
        assert_eq!(lines[3], "UUAA");
    }

    #[test]
    fn test_fastq_dna2rna() {
        let data = b"@seq1\nATGC\n+\nIIII\n@seq2\nTTAA\n+\nHHHH\n";
        let mut output = Vec::new();
        fasta_dna2rna(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("AUGC"));
        assert!(result.contains("UUAA"));
        // Verify quality is preserved
        assert!(result.contains("IIII"));
        assert!(result.contains("HHHH"));
    }
}
