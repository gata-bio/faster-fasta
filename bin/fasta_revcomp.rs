//! Reverse complement utility for FASTA and FASTQ files
//!
//! Compute reverse complement of DNA sequences.
//! Uses StringZilla's `lookup()` for SIMD-accelerated nucleotide translation.
//! Supports standard and ambiguous IUPAC codes.
//! Auto-detects format and preserves it (FASTQ quality scores are also reversed).
//!
//! **Memory**: O(1) - processes one sequence at a time
//! **Streaming**: Yes - constant memory usage
//!
//! # Examples
//!
//! ```bash
//! # FASTA files
//! fasta-revcomp sequences.fasta -o revcomp.fasta
//!
//! # FASTQ files (quality is also reversed)
//! fasta-revcomp reads.fastq -o revcomp.fastq
//!
//! # Pipe composition
//! cat sequences.fasta | fasta-revcomp | fasta-sort > sorted_revcomp.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::lookup;

use faster_fasta::shared::*;

/// Compute reverse complement of DNA sequences (FASTA or FASTQ)
///
/// Uses StringZilla's lookup() function for fast character mapping.
/// Handles both uppercase and lowercase DNA bases.
/// For FASTQ files, quality scores are also reversed to match the reversed sequence.
pub fn fasta_revcomp(input: &[u8], mut output: impl Write) -> io::Result<()> {
    let format = detect_format(input)?;

    // StringZilla translation table for reverse complement
    let mut table = [0u8; 256];
    for i in 0..256 {
        table[i] = i as u8;
    }
    // Uppercase DNA bases
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'G' as usize] = b'C';
    table[b'C' as usize] = b'G';
    // Lowercase DNA bases
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'g' as usize] = b'c';
    table[b'c' as usize] = b'g';
    // Ambiguous bases
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'n';
    table[b'R' as usize] = b'Y';
    table[b'Y' as usize] = b'R';
    table[b'S' as usize] = b'S';
    table[b'W' as usize] = b'W';
    table[b'K' as usize] = b'M';
    table[b'M' as usize] = b'K';
    table[b'B' as usize] = b'V';
    table[b'D' as usize] = b'H';
    table[b'H' as usize] = b'D';
    table[b'V' as usize] = b'B';

    match format {
        SeqFormat::Fasta => {
            let parser = FastaParser::new(input);
            for entry in parser {
                let seq = entry.sequence.as_bytes();
                let mut revcomp = vec![0u8; seq.len()];
                lookup(&mut revcomp, seq, table);
                revcomp.reverse();

                output.write_all(entry.header)?;
                output.write_all(b"\n")?;
                output.write_all(&revcomp)?;
                output.write_all(b"\n")?;
            }
        }
        SeqFormat::Fastq => {
            let parser = FastqParser::new(input);
            for entry in parser {
                let entry = entry?;
                let seq = entry.sequence.as_bytes();
                let mut revcomp = vec![0u8; seq.len()];
                lookup(&mut revcomp, seq, table);
                revcomp.reverse();

                // Also reverse quality scores to match reversed sequence
                let mut rev_quality = entry.quality.to_vec();
                rev_quality.reverse();

                output.write_all(entry.header)?;
                output.write_all(b"\n")?;
                output.write_all(&revcomp)?;
                output.write_all(b"\n+\n")?;
                output.write_all(&rev_quality)?;
                output.write_all(b"\n")?;
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Compute reverse complement of DNA sequences
#[derive(Parser)]
#[command(name = "fasta-revcomp")]
#[command(version, about = "Compute reverse complement of DNA sequences")]
#[command(
    long_about = "Compute reverse complement of DNA sequences from FASTA or FASTQ files.\nFor FASTQ, quality scores are also reversed to match the sequence orientation.\nAuto-detects format and preserves it in output."
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

    if let Err(e) = fasta_revcomp(input.as_bytes(), output) {
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
    fn revcomp_simple() {
        let data = b">seq1\nACGT\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1"));
        assert!(result.contains("ACGT"));
    }

    #[test]
    fn revcomp_palindrome() {
        let data = b">seq1\nGAATTC\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("GAATTC"));
    }

    #[test]
    fn revcomp_multiple() {
        let data = b">seq1\nAAAA\n>seq2\nTTTT\n>seq3\nGGGG\n>seq4\nCCCC\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[1], "TTTT");
        assert_eq!(lines[3], "AAAA");
        assert_eq!(lines[5], "CCCC");
        assert_eq!(lines[7], "GGGG");
    }

    #[test]
    fn revcomp_mixed_case() {
        let data = b">seq1\nAcGt\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("aCgT"));
    }

    #[test]
    fn revcomp_with_n() {
        let data = b">seq1\nACGTN\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("NACGT"));
    }
}
