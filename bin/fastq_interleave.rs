//! FASTQ paired-end interleaving utility
//!
//! Interleave two FASTQ files (R1 and R2) into a single file.
//! Reads are alternated: R1[0], R2[0], R1[1], R2[1], ...
//! Streaming implementation with O(1) memory usage.
//!
//! **Memory**: O(1) - processes one pair at a time
//! **Streaming**: Yes - no materialization required
//!
//! # Examples
//!
//! ```bash
//! # Interleave R1 and R2 files
//! fastq-interleave R1.fastq R2.fastq -o interleaved.fastq
//!
//! # To stdout
//! fastq-interleave R1.fastq R2.fastq > interleaved.fastq
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;

use faster_fasta::shared::*;

/// Interleave two FASTQ files
pub fn fastq_interleave(r1_data: &[u8], r2_data: &[u8], mut output: impl Write) -> io::Result<()> {
    let mut r1_parser = FastqParser::new(r1_data);
    let mut r2_parser = FastqParser::new(r2_data);

    loop {
        let r1_next = r1_parser.next();
        let r2_next = r2_parser.next();

        match (r1_next, r2_next) {
            (None, None) => break,
            (Some(Err(e)), _) | (_, Some(Err(e))) => return Err(e),
            (Some(Ok(entry1)), Some(Ok(entry2))) => {
                output.write_all(entry1.header)?;
                output.write_all(b"\n")?;
                output.write_all(entry1.sequence.as_bytes())?;
                output.write_all(b"\n+\n")?;
                output.write_all(entry1.quality)?;
                output.write_all(b"\n")?;

                output.write_all(entry2.header)?;
                output.write_all(b"\n")?;
                output.write_all(entry2.sequence.as_bytes())?;
                output.write_all(b"\n+\n")?;
                output.write_all(entry2.quality)?;
                output.write_all(b"\n")?;
            }
            (Some(Ok(_)), None) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FASTQ interleave mismatch: R1 has extra reads",
                ))
            }
            (None, Some(Ok(_))) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FASTQ interleave mismatch: R2 has extra reads",
                ))
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Interleave paired-end FASTQ files
#[derive(Parser)]
#[command(name = "fastq-interleave")]
#[command(version, about = "Interleave paired-end FASTQ files")]
#[command(
    long_about = "Interleave two FASTQ files (R1 and R2) into a single file.\nReads are alternated: R1[0], R2[0], R1[1], R2[1], ...\nStreaming implementation with minimal memory usage."
)]
struct Args {
    /// R1 (forward) FASTQ file
    r1: String,

    /// R2 (reverse) FASTQ file
    r2: String,

    /// Output interleaved FASTQ file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,
}

fn main() {
    let args = Args::parse();

    let r1_input = match get_input(Some(&args.r1)) {
        Ok(input) => input,
        Err(e) => {
            eprintln!("Error reading R1 file: {}", e);
            process::exit(1);
        }
    };

    let r2_input = match get_input(Some(&args.r2)) {
        Ok(input) => input,
        Err(e) => {
            eprintln!("Error reading R2 file: {}", e);
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

    if let Err(e) = fastq_interleave(r1_input.as_bytes(), r2_input.as_bytes(), output) {
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
    fn interleave_basic() {
        let r1 = b"@read1/1\nACGT\n+\nIIII\n@read2/1\nTGCA\n+\nHHHH\n";
        let r2 = b"@read1/2\nGGGG\n+\nJJJJ\n@read2/2\nCCCC\n+\nKKKK\n";
        let mut output = Vec::new();

        fastq_interleave(r1, r2, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Check alternation
        assert_eq!(lines[0], "@read1/1");
        assert_eq!(lines[4], "@read1/2");
        assert_eq!(lines[8], "@read2/1");
        assert_eq!(lines[12], "@read2/2");
    }

    #[test]
    fn interleave_mismatch_detected() {
        // R2 missing final mate
        let r1 = b"@r1/1\nAC\n+\nII\n@r2/1\nGG\n+\nII\n";
        let r2 = b"@r1/2\nTT\n+\nII\n";
        let mut output = Vec::new();
        let err = fastq_interleave(r1, r2, &mut output).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }
}
