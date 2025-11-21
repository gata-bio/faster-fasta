//! FASTQ quality filtering utility
//!
//! Filter FASTQ reads based on quality scores, length, and N-content.
//! Streaming implementation with O(1) memory usage.
//!
//! **Memory**: O(1) - processes one read at a time
//! **Streaming**: Yes - no materialization required
//!
//! # Examples
//!
//! ```bash
//! # Filter by minimum mean quality
//! fastq-filter reads.fastq --min-quality 20 -o filtered.fastq
//!
//! # Filter by length range and N-content
//! fastq-filter reads.fastq --min-length 50 --max-n-fraction 0.1 -o filtered.fastq
//!
//! # Multiple filters (AND logic)
//! fastq-filter reads.fastq --min-quality 25 --min-length 75 --max-length 150 -o filtered.fastq
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;

mod shared;
use shared::*;

/// Filter FASTQ reads based on quality and content criteria
///
/// Applies multiple filter criteria with AND logic: a read passes only if it meets ALL criteria.
/// Uses streaming to process one read at a time with minimal memory usage.
pub fn fastq_filter(
    input: &[u8],
    mut output: impl Write,
    min_quality: Option<f32>,
    min_length: Option<usize>,
    max_length: Option<usize>,
    max_n_fraction: Option<f32>,
) -> io::Result<()> {
    let parser = FastqParser::new(input);

    for entry in parser {
        let seq = entry.sequence.as_bytes();
        let qual = entry.quality;

        // Filter by mean quality
        if let Some(min_q) = min_quality {
            if quality::mean_quality(qual) < min_q {
                continue;
            }
        }

        // Filter by length
        let seq_len = seq.len();
        if let Some(min_len) = min_length {
            if seq_len < min_len {
                continue;
            }
        }
        if let Some(max_len) = max_length {
            if seq_len > max_len {
                continue;
            }
        }

        // Filter by N-base content
        if let Some(max_n) = max_n_fraction {
            let n_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
            let n_frac = n_count as f32 / seq_len as f32;
            if n_frac > max_n {
                continue;
            }
        }

        // Read passed all filters - write it
        output.write_all(entry.header)?;
        output.write_all(b"\n")?;
        output.write_all(seq)?;
        output.write_all(b"\n+\n")?;
        output.write_all(qual)?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Filter FASTQ reads by quality and content
#[derive(Parser)]
#[command(name = "fastq-filter")]
#[command(
    version,
    about = "Filter FASTQ reads by quality, length, and N-content"
)]
#[command(
    long_about = "Filter FASTQ reads based on quality scores, length, and N-base content.\nAll filters use AND logic: reads must pass ALL criteria.\nStreaming implementation with minimal memory usage."
)]
struct Args {
    /// Input FASTQ file (use '-' or omit for stdin)
    input: Option<String>,

    /// Output FASTQ file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,

    /// Minimum mean quality score (Phred scale)
    #[arg(short = 'q', long)]
    min_quality: Option<f32>,

    /// Minimum read length
    #[arg(short = 'l', long)]
    min_length: Option<usize>,

    /// Maximum read length
    #[arg(short = 'L', long)]
    max_length: Option<usize>,

    /// Maximum fraction of N-bases (0.0-1.0)
    #[arg(short = 'n', long)]
    max_n_fraction: Option<f32>,
}

fn main() {
    let args = Args::parse();

    // Validate max_n_fraction if provided
    if let Some(f) = args.max_n_fraction {
        if f < 0.0 || f > 1.0 {
            eprintln!("Error: max-n-fraction must be between 0.0 and 1.0");
            process::exit(1);
        }
    }

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

    if let Err(e) = fastq_filter(
        input.as_bytes(),
        output,
        args.min_quality,
        args.min_length,
        args.max_length,
        args.max_n_fraction,
    ) {
        if e.kind() == io::ErrorKind::BrokenPipe {
            process::exit(0);
        }
        eprintln!("Error processing FASTQ: {}", e);
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_by_quality() {
        // Quality: IIII = Q40, HHHH = Q39, !!!! = Q0
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n@seq3\nAAAA\n+\n!!!!\n";
        let mut output = Vec::new();
        fastq_filter(&data[..], &mut output, Some(20.0), None, None, None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("@seq1"));
        assert!(result.contains("@seq2"));
        assert!(!result.contains("@seq3")); // Low quality filtered
    }

    #[test]
    fn test_filter_by_length() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nAA\n+\nII\n@seq3\nACGTACGT\n+\nIIIIIIII\n";
        let mut output = Vec::new();
        fastq_filter(&data[..], &mut output, None, Some(4), Some(6), None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("@seq1")); // Length 4: passes
        assert!(!result.contains("@seq2")); // Length 2: too short
        assert!(!result.contains("@seq3")); // Length 8: too long
    }

    #[test]
    fn test_filter_by_n_content() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nNNNN\n+\nIIII\n@seq3\nACNT\n+\nIIII\n";
        let mut output = Vec::new();
        fastq_filter(&data[..], &mut output, None, None, None, Some(0.5)).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("@seq1")); // 0% N: passes
        assert!(!result.contains("@seq2")); // 100% N: fails
        assert!(result.contains("@seq3")); // 25% N: passes (< 50%)
    }

    #[test]
    fn test_multiple_filters() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nAC\n+\nII\n@seq3\nACGT\n+\n!!!!\n";
        let mut output = Vec::new();
        fastq_filter(&data[..], &mut output, Some(20.0), Some(4), None, None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("@seq1")); // Passes both filters
        assert!(!result.contains("@seq2")); // Fails length
        assert!(!result.contains("@seq3")); // Fails quality
    }
}
