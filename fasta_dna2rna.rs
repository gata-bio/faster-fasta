//! FASTA DNA to RNA conversion utility
//!
//! Convert DNA sequences to RNA (T→U).
//! Uses StringZilla's `lookup()` for SIMD-accelerated character mapping.
//!
//! **Memory**: O(1) - processes one sequence at a time
//! **Streaming**: Yes - constant memory usage
//!
//! # Examples
//!
//! ```bash
//! # From file
//! fasta-dna2rna sequences.fasta -o rna.fasta
//!
//! # From stdin
//! cat dna.fasta | fasta-dna2rna > rna.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::lookup;

mod faster_fasta;
use faster_fasta::*;

/// Convert DNA sequences to RNA
///
/// Converts T (thymine) to U (uracil) using StringZilla's translate() function.
/// Handles both uppercase and lowercase bases.
pub fn fasta_dna2rna(input: &[u8], mut output: impl Write) -> io::Result<()> {
    // StringZilla translation table for DNA -> RNA
    let mut table = [0u8; 256];
    for i in 0..256 {
        table[i] = i as u8;
    }
    // Convert T to U
    table[b'T' as usize] = b'U';
    table[b't' as usize] = b'u';

    let parser = FastaParser::new(input);

    for entry in parser {
        let seq = entry.sequence.as_bytes();
        let mut rna = vec![0u8; seq.len()];

        // Use StringZilla lookup for DNA -> RNA conversion
        lookup(&mut rna, seq, table);

        output.write_all(entry.header)?;
        output.write_all(b"\n")?;
        output.write_all(&rna)?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Convert DNA sequences to RNA (T -> U)
#[derive(Parser)]
#[command(name = "fasta-dna2rna")]
#[command(version, about, long_about = None)]
struct Args {
    /// Input FASTA file (use '-' or omit for stdin)
    input: Option<String>,

    /// Output FASTA file (use '-' or omit for stdout)
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
        eprintln!("Error processing FASTA: {}", e);
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
}
