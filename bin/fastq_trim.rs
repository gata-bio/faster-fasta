//! FASTQ trimming utility
//!
//! Trim reads based on quality scores or fixed positions.
//! Streaming implementation with O(1) memory usage.
//!
//! **Memory**: O(1) - processes one read at a time
//! **Streaming**: Yes - no materialization required
//!
//! # Examples
//!
//! ```bash
//! # Trim low-quality ends
//! fastq-trim reads.fastq --quality-cutoff 20 -o trimmed.fastq
//!
//! # Fixed-length trimming
//! fastq-trim reads.fastq --trim-front 5 --trim-tail 10 -o trimmed.fastq
//!
//! # Truncate to max length and filter short reads
//! fastq-trim reads.fastq --max-length 100 --min-length 50 -o trimmed.fastq
//! ```

use std::io::{self, Read, Write};
use std::process;

use clap::Parser;

use faster_fasta::shared::*;

/// Trim a read based on various criteria
fn trim_read(
    seq: &[u8],
    qual: &[u8],
    quality_cutoff: Option<u8>,
    trim_front: usize,
    trim_tail: usize,
    max_length: Option<usize>,
) -> (Vec<u8>, Vec<u8>) {
    let mut start = 0;
    let mut end = seq.len();

    // Quality-based trimming (trim from both ends until quality >= cutoff)
    if let Some(cutoff) = quality_cutoff {
        // Trim front
        while start < end && quality::ascii_to_phred33(qual[start]) < cutoff {
            start += 1;
        }

        // Trim tail
        while end > start && quality::ascii_to_phred33(qual[end - 1]) < cutoff {
            end -= 1;
        }
    }

    // Fixed-position trimming
    start += trim_front;
    if end > trim_tail {
        end -= trim_tail;
    } else {
        end = start; // Read becomes empty
    }

    // Ensure start doesn't exceed end
    start = start.min(end).min(seq.len());
    end = end.min(seq.len());

    // Truncate to max length
    if let Some(max_len) = max_length {
        end = end.min(start + max_len);
    }

    let trimmed_seq = seq[start..end].to_vec();
    let trimmed_qual = qual[start..end].to_vec();

    (trimmed_seq, trimmed_qual)
}

/// Trim FASTQ reads
pub fn fastq_trim(
    input: &[u8],
    mut output: impl Write,
    quality_cutoff: Option<u8>,
    trim_front: usize,
    trim_tail: usize,
    max_length: Option<usize>,
    min_length: Option<usize>,
) -> io::Result<()> {
    let parser = FastqParser::new(input);

    for entry in parser {
        let entry = entry?;
        let seq = entry.sequence.as_bytes();
        let qual = entry.quality;

        let (trimmed_seq, trimmed_qual) =
            trim_read(seq, qual, quality_cutoff, trim_front, trim_tail, max_length);

        // Filter by minimum length
        if let Some(min_len) = min_length {
            if trimmed_seq.len() < min_len {
                continue;
            }
        }

        // Write trimmed read
        output.write_all(entry.header)?;
        output.write_all(b"\n")?;
        output.write_all(&trimmed_seq)?;
        output.write_all(b"\n+\n")?;
        output.write_all(&trimmed_qual)?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

fn fastq_trim_stream(
    mut parser: FastqStreamParser<impl Read>,
    mut output: impl Write,
    quality_cutoff: Option<u8>,
    trim_front: usize,
    trim_tail: usize,
    max_length: Option<usize>,
    min_length: Option<usize>,
) -> io::Result<()> {
    while let Some(entry) = parser.next_entry()? {
        let (trimmed_seq, trimmed_qual) =
            trim_read(&entry.sequence, &entry.quality, quality_cutoff, trim_front, trim_tail, max_length);

        if let Some(min_len) = min_length {
            if trimmed_seq.len() < min_len {
                continue;
            }
        }

        output.write_all(&entry.header)?;
        output.write_all(b"\n")?;
        output.write_all(&trimmed_seq)?;
        output.write_all(b"\n+\n")?;
        output.write_all(&trimmed_qual)?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Trim FASTQ reads by quality or position
#[derive(Parser)]
#[command(name = "fastq-trim")]
#[command(version, about = "Trim FASTQ reads by quality or position")]
#[command(
    long_about = "Trim FASTQ reads based on quality scores or fixed positions.\nMultiple trimming operations applied in order: quality → front → tail → max-length.\nReads shorter than min-length after trimming are discarded."
)]
struct Args {
    /// Input FASTQ file (use '-' or omit for stdin)
    input: Option<String>,

    /// Output FASTQ file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,

    /// Quality cutoff (trim ends until quality >= cutoff)
    #[arg(short = 'q', long)]
    quality_cutoff: Option<u8>,

    /// Remove N bases from front
    #[arg(short = 'f', long, default_value = "0")]
    trim_front: usize,

    /// Remove N bases from tail
    #[arg(short = 't', long, default_value = "0")]
    trim_tail: usize,

    /// Truncate to maximum length
    #[arg(short = 'L', long)]
    max_length: Option<usize>,

    /// Discard reads shorter than minimum length (after trimming)
    #[arg(short = 'l', long)]
    min_length: Option<usize>,
}

fn main() {
    let args = Args::parse();

    let use_stream = matches!(args.input.as_deref(), None | Some("-"));

    if use_stream {
        let output = match get_output(args.output.as_deref()) {
            Ok(output) => output,
            Err(e) => {
                eprintln!("Error opening output: {}", e);
                process::exit(1);
            }
        };

        if let Err(e) = fastq_trim_stream(
            FastqStreamParser::new(io::stdin()),
            output,
            args.quality_cutoff,
            args.trim_front,
            args.trim_tail,
            args.max_length,
            args.min_length,
        ) {
            if e.kind() == io::ErrorKind::BrokenPipe {
                process::exit(0);
            }
            eprintln!("Error processing FASTQ: {}", e);
            process::exit(1);
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

        if let Err(e) = fastq_trim(
            input.as_bytes(),
            output,
            args.quality_cutoff,
            args.trim_front,
            args.trim_tail,
            args.max_length,
            args.min_length,
        ) {
            if e.kind() == io::ErrorKind::BrokenPipe {
                process::exit(0);
            }
            eprintln!("Error processing FASTQ: {}", e);
            process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trim_front() {
        let data = b"@seq1\nACGTACGT\n+\nIIIIIIII\n";
        let mut output = Vec::new();
        fastq_trim(&data[..], &mut output, None, 2, 0, None, None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("GTACGT"));
    }

    #[test]
    fn trim_tail() {
        let data = b"@seq1\nACGTACGT\n+\nIIIIIIII\n";
        let mut output = Vec::new();
        fastq_trim(&data[..], &mut output, None, 0, 2, None, None).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("ACGTAC"));
    }

    #[test]
    fn min_length_filter() {
        let data = b"@seq1\nACGT\n+\nIIII\n";
        let mut output = Vec::new();
        fastq_trim(&data[..], &mut output, None, 2, 0, None, Some(10)).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.is_empty()); // Read filtered out
    }
}
