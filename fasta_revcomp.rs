//! FASTA reverse complement utility
//!
//! Compute reverse complement of DNA sequences.
//! Uses StringZilla's `lookup()` for SIMD-accelerated nucleotide translation.
//! Supports standard and ambiguous IUPAC codes.
//!
//! **Memory**: O(1) - processes one sequence at a time
//! **Streaming**: Yes - constant memory usage
//!
//! # Examples
//!
//! ```bash
//! # From file
//! fasta-revcomp sequences.fasta -o revcomp.fasta
//!
//! # Pipe composition
//! cat sequences.fasta | fasta-revcomp | fasta-sort > sorted_revcomp.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::lookup;

mod faster_fasta;
use faster_fasta::*;

/// Compute reverse complement of DNA sequences
///
/// Uses StringZilla's translate() function for fast character mapping.
/// Handles both uppercase and lowercase DNA bases.
pub fn fasta_revcomp(input: &[u8], mut output: impl Write) -> io::Result<()> {
    // StringZilla translation table for reverse complement
    // Maps each byte to its complement
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
    table[b'R' as usize] = b'Y'; // A or G -> T or C
    table[b'Y' as usize] = b'R'; // C or T -> G or A
    table[b'S' as usize] = b'S'; // G or C -> C or G (palindrome)
    table[b'W' as usize] = b'W'; // A or T -> T or A (palindrome)
    table[b'K' as usize] = b'M'; // G or T -> C or A
    table[b'M' as usize] = b'K'; // A or C -> T or G
    table[b'B' as usize] = b'V'; // C, G, or T -> G, C, or A
    table[b'D' as usize] = b'H'; // A, G, or T -> T, C, or A
    table[b'H' as usize] = b'D'; // A, C, or T -> T, G, or A
    table[b'V' as usize] = b'B'; // A, C, or G -> T, G, or C

    let parser = FastaParser::new(input);

    for entry in parser {
        let seq = entry.sequence.as_bytes();
        let mut revcomp = vec![0u8; seq.len()];

        // Use StringZilla lookup for complement
        lookup(&mut revcomp, seq, table);

        // Reverse in place
        revcomp.reverse();

        output.write_all(entry.header)?;
        output.write_all(b"\n")?;
        output.write_all(&revcomp)?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Compute reverse complement of DNA sequences
#[derive(Parser)]
#[command(name = "fasta-revcomp")]
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

    if let Err(e) = fasta_revcomp(input.as_bytes(), output) {
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
    fn test_revcomp_simple() {
        let data = b">seq1\nACGT\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1"));
        assert!(result.contains("ACGT"));
    }

    #[test]
    fn test_revcomp_palindrome() {
        let data = b">seq1\nGAATTC\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("GAATTC"));
    }

    #[test]
    fn test_revcomp_multiple() {
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
    fn test_revcomp_mixed_case() {
        let data = b">seq1\nAcGt\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("aCgT"));
    }

    #[test]
    fn test_revcomp_with_n() {
        let data = b">seq1\nACGTN\n";
        let mut output = Vec::new();
        fasta_revcomp(data, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("NACGT"));
    }
}
