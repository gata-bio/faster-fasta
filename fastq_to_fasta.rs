//! FASTQ to FASTA conversion utility
//!
//! Convert FASTQ files to FASTA format by dropping quality scores.
//! Optional quality filtering before conversion.
//! Streaming implementation with O(1) memory usage.
//!
//! **Memory**: O(1) - processes one read at a time
//! **Streaming**: Yes - no materialization required
//!
//! # Examples
//!
//! ```bash
//! # Simple conversion
//! fastq-to-fasta reads.fastq -o sequences.fasta
//!
//! # Convert only high-quality reads
//! fastq-to-fasta reads.fastq --min-quality 25 -o high_quality.fasta
//!
//! # From stdin
//! cat reads.fastq | fastq-to-fasta > sequences.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;

mod shared;
use shared::*;

/// Convert FASTQ to FASTA format
///
/// Drops quality scores and converts @ headers to > headers.
/// Optional quality filtering: only convert reads meeting quality threshold.
pub fn fastq_to_fasta(
    input: &[u8],
    mut output: impl Write,
    min_quality: Option<f32>,
) -> io::Result<()> {
    let parser = FastqParser::new(input);

    for entry in parser {
        // Optional quality filter
        if let Some(min_q) = min_quality {
            if quality::mean_quality(entry.quality) < min_q {
                continue;
            }
        }

        // Convert @ to > for FASTA header
        output.write_all(b">")?;
        // Skip the '@' character from FASTQ header
        if entry.header.len() > 1 {
            output.write_all(&entry.header[1..])?;
        }
        output.write_all(b"\n")?;
        output.write_all(entry.sequence.as_bytes())?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Convert FASTQ to FASTA format
#[derive(Parser)]
#[command(name = "fastq-to-fasta")]
#[command(version, about = "Convert FASTQ to FASTA format")]
#[command(
    long_about = "Convert FASTQ files to FASTA format by dropping quality scores.\nOptional quality filtering before conversion.\nStreaming implementation with minimal memory usage."
)]
struct Args {
    /// Input FASTQ file (use '-' or omit for stdin)
    input: Option<String>,

    /// Output FASTA file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,

    /// Only convert reads with mean quality >= threshold
    #[arg(short = 'q', long)]
    min_quality: Option<f32>,
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

    if let Err(e) = fastq_to_fasta(input.as_bytes(), output, args.min_quality) {
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
    fn test_basic_conversion() {
        let data = b"@seq1 description\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n";
        let mut output = Vec::new();
        fastq_to_fasta(data, &mut output, None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1 description"));
        assert!(result.contains(">seq2"));
        assert!(result.contains("ACGT"));
        assert!(result.contains("TGCA"));
        // Ensure no FASTQ artifacts
        assert!(!result.contains("@"));
        assert!(!result.contains("+"));
        assert!(!result.contains("IIII"));
    }

    #[test]
    fn test_with_quality_filter() {
        // IIII = Q40, !!!! = Q0
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\n!!!!\n";
        let mut output = Vec::new();
        fastq_to_fasta(data, &mut output, Some(20.0)).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(">seq1"));
        assert!(result.contains("ACGT"));
        assert!(!result.contains(">seq2")); // Filtered out
        assert!(!result.contains("TGCA"));
    }

    #[test]
    fn test_empty_input() {
        let data = b"";
        let mut output = Vec::new();
        fastq_to_fasta(data, &mut output, None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert_eq!(result, "");
    }
}
