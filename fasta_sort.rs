//! FASTA sequence sorting utility
//!
//! Sort sequences by name (header), sequence content, or length.
//! Uses StringZilla's `argsort_permutation_by()` for SIMD-accelerated string sorting.
//!
//! **Memory**: O(n) - loads all entries for sorting
//! **Streaming**: No - sorting requires all data in memory
//!
//! # Examples
//!
//! ```bash
//! # Sort by header name (default)
//! fasta-sort sequences.fasta -o sorted.fasta
//!
//! # Sort by sequence content
//! fasta-sort --sequence sequences.fasta -o sorted.fasta
//!
//! # Sort by length (shortest first)
//! fasta-sort --length sequences.fasta -o sorted.fasta
//!
//! # Reverse sort by length (longest first)
//! cat sequences.fasta | fasta-sort -l -r > sorted.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::{argsort_permutation_by, SortedIdx};

mod faster_fasta;
use faster_fasta::*;

/// Sort criterion for FASTA sequences
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SortBy {
    Name,
    Sequence,
    Length,
}

/// Sort FASTA sequences
///
/// Sorts sequences by name (header), sequence content, or length.
/// Uses StringZilla's argsort_permutation_by() for string sorting.
pub fn fasta_sort(
    input: &[u8],
    mut output: impl Write,
    sort_by: SortBy,
    reverse: bool,
) -> io::Result<()> {
    let entries: Vec<_> = FastaParser::new(input).collect();
    let count = entries.len();

    if count == 0 {
        return Ok(());
    }

    let mut order: Vec<SortedIdx> = vec![0; count];

    match sort_by {
        SortBy::Length => {
            // Use standard sort for numeric length comparison
            for i in 0..count {
                order[i] = i;
            }
            order.sort_by(|&a, &b| {
                let len_a = entries[a].sequence.as_bytes().len();
                let len_b = entries[b].sequence.as_bytes().len();
                if reverse {
                    len_b.cmp(&len_a)
                } else {
                    len_a.cmp(&len_b)
                }
            });
        }
        SortBy::Name | SortBy::Sequence => {
            // Use StringZilla for string sorting
            let mapper = |idx: usize| -> &[u8] {
                match sort_by {
                    SortBy::Name => entries[idx].header,
                    SortBy::Sequence => entries[idx].sequence.as_bytes(),
                    _ => unreachable!(),
                }
            };

            argsort_permutation_by(mapper, &mut order)
                .map_err(|_| io::Error::new(io::ErrorKind::Other, "Sort failed"))?;

            if reverse {
                order.reverse();
            }
        }
    }

    // Write sorted entries
    for &idx in &order {
        let entry = &entries[idx];
        output.write_all(entry.header)?;
        output.write_all(b"\n")?;
        output.write_all(entry.sequence.as_bytes())?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Sort FASTA sequences by various criteria
#[derive(Parser)]
#[command(name = "fasta-sort")]
#[command(version, about, long_about = None)]
struct Args {
    /// Input FASTA file (use '-' or omit for stdin)
    input: Option<String>,

    /// Output FASTA file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,

    /// Sort by name (header)
    #[arg(short = 'n', long, conflicts_with_all = ["sequence", "length"])]
    name: bool,

    /// Sort by sequence
    #[arg(short = 's', long, conflicts_with_all = ["name", "length"])]
    sequence: bool,

    /// Sort by sequence length
    #[arg(short = 'l', long, conflicts_with_all = ["name", "sequence"])]
    length: bool,

    /// Reverse sort order
    #[arg(short, long)]
    reverse: bool,
}

fn main() {
    let args = Args::parse();

    let sort_by = if args.sequence {
        SortBy::Sequence
    } else if args.length {
        SortBy::Length
    } else {
        SortBy::Name // Default
    };

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

    if let Err(e) = fasta_sort(input.as_bytes(), output, sort_by, args.reverse) {
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
    fn test_sort_by_name() {
        let data = b">seq3\nACGT\n>seq1\nTGCA\n>seq2\nAAAA\n";
        let mut output = Vec::new();
        fasta_sort(data, &mut output, SortBy::Name, false).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[0], ">seq1");
        assert_eq!(lines[2], ">seq2");
        assert_eq!(lines[4], ">seq3");
    }

    #[test]
    fn test_sort_by_sequence() {
        let data = b">seq1\nTGCA\n>seq2\nACGT\n>seq3\nGGGG\n";
        let mut output = Vec::new();
        fasta_sort(data, &mut output, SortBy::Sequence, false).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[1], "ACGT");
        assert_eq!(lines[3], "GGGG");
        assert_eq!(lines[5], "TGCA");
    }

    #[test]
    fn test_sort_by_length() {
        let data = b">seq1\nACGTACGT\n>seq2\nAA\n>seq3\nTGCA\n";
        let mut output = Vec::new();
        fasta_sort(data, &mut output, SortBy::Length, false).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[0], ">seq2");
        assert_eq!(lines[1], "AA");
        assert_eq!(lines[2], ">seq3");
        assert_eq!(lines[3], "TGCA");
        assert_eq!(lines[4], ">seq1");
        assert_eq!(lines[5], "ACGTACGT");
    }

    #[test]
    fn test_sort_reverse() {
        let data = b">seq1\nAA\n>seq2\nTTTT\n>seq3\nGG\n";
        let mut output = Vec::new();
        fasta_sort(data, &mut output, SortBy::Length, true).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[0], ">seq2");
        assert_eq!(lines[1], "TTTT");
    }
}
