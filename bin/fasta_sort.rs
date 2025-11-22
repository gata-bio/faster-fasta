//! Sequence sorting utility for FASTA and FASTQ files
//!
//! Sort sequences by name (header), sequence content, or length.
//! Uses StringZilla's `argsort_permutation_by()` for SIMD-accelerated string sorting.
//! Auto-detects format and preserves it in output.
//!
//! **Memory**: O(n) - loads all entries for sorting
//! **Streaming**: No - sorting requires all data in memory
//!
//! # Examples
//!
//! ```bash
//! # Sort FASTA by header name (default)
//! fasta-sort sequences.fasta -o sorted.fasta
//!
//! # Sort FASTQ by sequence content
//! fasta-sort --sequence reads.fastq -o sorted.fastq
//!
//! # Sort by length (shortest first)
//! fasta-sort --length sequences.fasta -o sorted.fasta
//!
//! # Reverse sort by length (longest first)
//! cat sequences.fasta | fasta-sort -l -r > sorted.fasta
//! ```

use std::io::{self, Read, Write};
use std::process;

use clap::Parser;
use stringzilla::sz::{argsort_permutation_by, SortedIdx};

use faster_fasta::shared::*;

/// Sort criterion for sequences
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SortBy {
    Name,
    Sequence,
    Length,
}

/// Stored entry for sorting (supports both FASTA and FASTQ)
enum StoredEntry {
    Fasta {
        header: Vec<u8>,
        sequence: Vec<u8>,
    },
    Fastq {
        header: Vec<u8>,
        sequence: Vec<u8>,
        quality: Vec<u8>,
    },
}

impl StoredEntry {
    fn header(&self) -> &[u8] {
        match self {
            StoredEntry::Fasta { header, .. } => header,
            StoredEntry::Fastq { header, .. } => header,
        }
    }

    fn sequence(&self) -> &[u8] {
        match self {
            StoredEntry::Fasta { sequence, .. } => sequence,
            StoredEntry::Fastq { sequence, .. } => sequence,
        }
    }
}

/// Sort sequences (FASTA or FASTQ)
///
/// Sorts sequences by name (header), sequence content, or length.
/// Uses StringZilla's argsort_permutation_by() for string sorting.
/// Auto-detects format and preserves it in output.
pub fn fasta_sort(
    input: &[u8],
    output: impl Write,
    sort_by: SortBy,
    reverse: bool,
) -> io::Result<()> {
    let format = detect_format(input)?;

    let entries: Vec<StoredEntry> = match format {
        SeqFormat::Fasta => FastaParser::new(input)
            .map(|e| StoredEntry::Fasta {
                header: e.header.to_vec(),
                sequence: e.sequence.as_bytes().to_vec(),
            })
            .collect(),
        SeqFormat::Fastq => FastqParser::new(input)
            .map(|e| {
                let e = e?;
                Ok(StoredEntry::Fastq {
                    header: e.header.to_vec(),
                    sequence: e.sequence.as_bytes().to_vec(),
                    quality: e.quality.to_vec(),
                })
            })
            .collect::<io::Result<_>>()?,
    };

    sort_entries(entries, output, sort_by, reverse)
}

fn sort_entries(
    entries: Vec<StoredEntry>,
    mut output: impl Write,
    sort_by: SortBy,
    reverse: bool,
) -> io::Result<()> {
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
                let len_a = entries[a].sequence().len();
                let len_b = entries[b].sequence().len();
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
                    SortBy::Name => entries[idx].header(),
                    SortBy::Sequence => entries[idx].sequence(),
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
        match &entries[idx] {
            StoredEntry::Fasta { header, sequence } => {
                output.write_all(header)?;
                output.write_all(b"\n")?;
                output.write_all(sequence)?;
                output.write_all(b"\n")?;
            }
            StoredEntry::Fastq {
                header,
                sequence,
                quality,
            } => {
                output.write_all(header)?;
                output.write_all(b"\n")?;
                output.write_all(sequence)?;
                output.write_all(b"\n+\n")?;
                output.write_all(quality)?;
                output.write_all(b"\n")?;
            }
        }
    }

    output.flush()?;
    Ok(())
}

fn collect_fastq_entries_stream(
    mut parser: FastqStreamParser<impl Read>,
) -> io::Result<Vec<StoredEntry>> {
    let mut entries = Vec::new();
    while let Some(entry) = parser.next_entry()? {
        entries.push(StoredEntry::Fastq {
            header: entry.header,
            sequence: entry.sequence,
            quality: entry.quality,
        });
    }
    Ok(entries)
}

/// Sort sequences by various criteria
#[derive(Parser)]
#[command(name = "fasta-sort")]
#[command(version, about = "Sort sequences by name, content, or length")]
#[command(
    long_about = "Sort sequences from FASTA or FASTQ files.\nAuto-detects format and preserves it in output.\nUses SIMD-accelerated string sorting for maximum performance."
)]
struct Args {
    /// Input file (FASTA or FASTQ, use '-' or omit for stdin)
    input: Option<String>,

    /// Output file (use '-' or omit for stdout)
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

    let use_stream = matches!(args.input.as_deref(), None | Some("-"));

    let result = if use_stream {
        let (format, mut reader) = match stdin_with_peek(8192) {
            Ok(v) => v,
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

        match format {
            SeqFormat::Fasta => {
                let mut data = Vec::new();
                if let Err(e) = reader.read_to_end(&mut data) {
                    Err(e)
                } else {
                    fasta_sort(&data, output, sort_by, args.reverse)
                }
            }
            SeqFormat::Fastq => {
                match collect_fastq_entries_stream(FastqStreamParser::new(reader)) {
                    Ok(entries) => sort_entries(entries, output, sort_by, args.reverse),
                    Err(e) => Err(e),
                }
            }
        }
    } else {
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

        fasta_sort(input.as_bytes(), output, sort_by, args.reverse)
    };

    if let Err(e) = result {
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
    fn sort_by_name() {
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
    fn sort_by_sequence() {
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
    fn sort_by_length() {
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
    fn sort_reverse() {
        let data = b">seq1\nAA\n>seq2\nTTTT\n>seq3\nGG\n";
        let mut output = Vec::new();
        fasta_sort(data, &mut output, SortBy::Length, true).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[0], ">seq2");
        assert_eq!(lines[1], "TTTT");
    }

    #[test]
    fn sort_fastq_by_length() {
        let data = b"@seq1\nACGTACGT\n+\nIIIIIIII\n@seq2\nAA\n+\nHH\n@seq3\nTGCA\n+\nJJJJ\n";
        let mut output = Vec::new();
        fasta_sort(data, &mut output, SortBy::Length, false).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("@seq2"));
        assert!(result.contains("@seq3"));
        assert!(result.contains("@seq1"));
        // Verify quality is preserved
        assert!(result.contains("HH"));
        assert!(result.contains("JJJJ"));
        assert!(result.contains("IIIIIIII"));
    }
}
