//! FASTQ paired-end deinterleaving utility
//!
//! Deinterleave a single FASTQ file into two files (R1 and R2).
//! Alternating reads are split: even positions → R1, odd positions → R2.
//! Streaming implementation with O(1) memory usage.
//!
//! **Memory**: O(1) - processes one read at a time
//! **Streaming**: Yes - no materialization required
//!
//! # Examples
//!
//! ```bash
//! # Deinterleave into R1 and R2 files
//! fastq-deinterleave interleaved.fastq -1 R1.fastq -2 R2.fastq
//!
//! # From stdin
//! cat interleaved.fastq | fastq-deinterleave -1 R1.fastq -2 R2.fastq
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;

mod shared;
use shared::*;

/// Deinterleave a FASTQ file into R1 and R2
pub fn fastq_deinterleave(
    input: &[u8],
    mut r1_out: impl Write,
    mut r2_out: impl Write,
) -> io::Result<()> {
    let parser = FastqParser::new(input);

    for (i, entry) in parser.enumerate() {
        if i % 2 == 0 {
            // R1 (even indices: 0, 2, 4, ...)
            r1_out.write_all(entry.header)?;
            r1_out.write_all(b"\n")?;
            r1_out.write_all(entry.sequence.as_bytes())?;
            r1_out.write_all(b"\n+\n")?;
            r1_out.write_all(entry.quality)?;
            r1_out.write_all(b"\n")?;
        } else {
            // R2 (odd indices: 1, 3, 5, ...)
            r2_out.write_all(entry.header)?;
            r2_out.write_all(b"\n")?;
            r2_out.write_all(entry.sequence.as_bytes())?;
            r2_out.write_all(b"\n+\n")?;
            r2_out.write_all(entry.quality)?;
            r2_out.write_all(b"\n")?;
        }
    }

    r1_out.flush()?;
    r2_out.flush()?;
    Ok(())
}

/// Deinterleave paired-end FASTQ files
#[derive(Parser)]
#[command(name = "fastq-deinterleave")]
#[command(version, about = "Deinterleave paired-end FASTQ files")]
#[command(
    long_about = "Deinterleave a single FASTQ file into two files (R1 and R2).\nAlternating reads are split: even positions → R1, odd positions → R2.\nStreaming implementation with minimal memory usage."
)]
struct Args {
    /// Input interleaved FASTQ file (use '-' or omit for stdin)
    input: Option<String>,

    /// R1 (forward) output file
    #[arg(short = '1', long = "r1", required = true)]
    r1_output: String,

    /// R2 (reverse) output file
    #[arg(short = '2', long = "r2", required = true)]
    r2_output: String,
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

    let r1_output = match get_output(Some(&args.r1_output)) {
        Ok(output) => output,
        Err(e) => {
            eprintln!("Error opening R1 output: {}", e);
            process::exit(1);
        }
    };

    let r2_output = match get_output(Some(&args.r2_output)) {
        Ok(output) => output,
        Err(e) => {
            eprintln!("Error opening R2 output: {}", e);
            process::exit(1);
        }
    };

    if let Err(e) = fastq_deinterleave(input.as_bytes(), r1_output, r2_output) {
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
    fn test_deinterleave_basic() {
        let data = b"@read1/1\nACGT\n+\nIIII\n@read1/2\nGGGG\n+\nJJJJ\n@read2/1\nTGCA\n+\nHHHH\n@read2/2\nCCCC\n+\nKKKK\n";
        let mut r1_out = Vec::new();
        let mut r2_out = Vec::new();

        fastq_deinterleave(data, &mut r1_out, &mut r2_out).unwrap();

        let r1 = String::from_utf8(r1_out).unwrap();
        let r2 = String::from_utf8(r2_out).unwrap();

        // Check R1 contains /1 reads
        assert!(r1.contains("@read1/1"));
        assert!(r1.contains("@read2/1"));
        assert!(r1.contains("ACGT"));
        assert!(r1.contains("TGCA"));

        // Check R2 contains /2 reads
        assert!(r2.contains("@read1/2"));
        assert!(r2.contains("@read2/2"));
        assert!(r2.contains("GGGG"));
        assert!(r2.contains("CCCC"));
    }
}
