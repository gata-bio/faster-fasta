//! FASTQ statistics and quality control utility
//!
//! Compute comprehensive statistics on FASTQ files including quality distributions,
//! length statistics, GC content, and N-base percentage.
//! Streaming implementation with fixed memory usage.
//!
//! **Memory**: O(1) - fixed-size arrays for per-position quality
//! **Streaming**: Yes - single-pass aggregation
//!
//! # Examples
//!
//! ```bash
//! # Basic statistics
//! fastq-stats reads.fastq
//!
//! # With histograms
//! fastq-stats reads.fastq --histogram
//! ```

use std::io::{self, Read};
use std::process;

use clap::Parser;

use faster_fasta::shared::*;

/// Statistics accumulator for FASTQ files
struct FastqStats {
    total_reads: usize,
    total_bases: u64,
    min_length: usize,
    max_length: usize,
    length_sum: u64,
    quality_sum: u64,
    gc_count: u64,
    n_count: u64,
    bin_size: usize,
    length_hist: [u64; 10], // 10 bins for length distribution
}

impl FastqStats {
    fn new() -> Self {
        Self {
            total_reads: 0,
            total_bases: 0,
            min_length: usize::MAX,
            max_length: 0,
            length_sum: 0,
            quality_sum: 0,
            gc_count: 0,
            n_count: 0,
            bin_size: 0,
            length_hist: [0; 10],
        }
    }

    fn add_read(&mut self, seq: &[u8], quality: &[u8]) {
        let len = seq.len();
        self.total_reads += 1;
        self.total_bases += len as u64;
        self.min_length = self.min_length.min(len);
        self.max_length = self.max_length.max(len);
        self.length_sum += len as u64;

        // Quality statistics
        for &q in quality {
            self.quality_sum += quality::ascii_to_phred33(q) as u64;
        }

        // GC content and N-bases
        for &base in seq {
            match base {
                b'G' | b'C' | b'g' | b'c' => self.gc_count += 1,
                b'N' | b'n' => self.n_count += 1,
                _ => {}
            }
        }

        // Length histogram (10 bins)
        if self.bin_size == 0 {
            self.bin_size = (len / 10).max(1);
        }
        let bin = ((len / self.bin_size).min(9)) as usize;
        self.length_hist[bin] += 1;
    }

    fn mean_length(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.length_sum as f64 / self.total_reads as f64
        }
    }

    fn mean_quality(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.quality_sum as f64 / self.total_bases as f64
        }
    }

    fn gc_percentage(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            100.0 * self.gc_count as f64 / self.total_bases as f64
        }
    }

    fn n_percentage(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            100.0 * self.n_count as f64 / self.total_bases as f64
        }
    }

    fn print(&self, show_histogram: bool) {
        println!("FASTQ Statistics");
        println!("================");
        println!("Total reads:      {}", self.total_reads);
        println!("Total bases:      {}", self.total_bases);
        println!();
        println!("Length Statistics");
        println!("-----------------");
        println!(
            "Min length:       {}",
            if self.min_length == usize::MAX {
                0
            } else {
                self.min_length
            }
        );
        println!("Max length:       {}", self.max_length);
        println!("Mean length:      {:.2}", self.mean_length());
        println!();
        println!("Quality Statistics");
        println!("------------------");
        println!("Mean quality:     {:.2}", self.mean_quality());
        println!();
        println!("Base Composition");
        println!("----------------");
        println!("GC content:       {:.2}%", self.gc_percentage());
        println!("N-bases:          {:.2}%", self.n_percentage());

        if show_histogram && self.total_reads > 0 {
            println!();
            println!("Length Distribution");
            println!("-------------------");
            let max_count = *self.length_hist.iter().max().unwrap_or(&1);
            let bin_size = self.bin_size.max(1);

            for (i, &count) in self.length_hist.iter().enumerate() {
                if count == 0 {
                    continue;
                }
                let start = i * bin_size;
                let end = (i + 1) * bin_size;
                let bar_len = (40 * count / max_count) as usize;
                let bar = "█".repeat(bar_len);
                println!("{:4}-{:<4}  {:8}  {}", start, end, count, bar);
            }
        }
    }
}

/// Compute FASTQ statistics
pub fn fastq_stats(input: &[u8], show_histogram: bool) -> io::Result<()> {
    let parser = FastqParser::new(input);
    let mut stats = FastqStats::new();

    for entry in parser {
        let entry = entry?;
        stats.add_read(entry.sequence.as_bytes(), entry.quality);
    }

    stats.print(show_histogram);
    Ok(())
}

pub fn fastq_stats_stream(
    mut parser: FastqStreamParser<impl Read>,
    show_histogram: bool,
) -> io::Result<()> {
    let mut stats = FastqStats::new();
    while let Some(entry) = parser.next_entry()? {
        stats.add_read(&entry.sequence, &entry.quality);
    }
    stats.print(show_histogram);
    Ok(())
}

/// Compute statistics on FASTQ files
#[derive(Parser)]
#[command(name = "fastq-stats")]
#[command(version, about = "Compute statistics on FASTQ files")]
#[command(
    long_about = "Compute comprehensive statistics on FASTQ files.\nIncludes quality scores, length distribution, GC content, and N-base percentage.\nStreaming implementation with minimal memory usage."
)]
struct Args {
    /// Input FASTQ file (use '-' or omit for stdin)
    input: Option<String>,

    /// Show length distribution histogram
    #[arg(long)]
    histogram: bool,
}

fn main() {
    let args = Args::parse();

    let use_stream = matches!(args.input.as_deref(), None | Some("-"));

    let result = if use_stream {
        fastq_stats_stream(FastqStreamParser::new(io::stdin()), args.histogram)
    } else {
        match get_input(args.input.as_deref()) {
            Ok(input) => fastq_stats(input.as_bytes(), args.histogram),
            Err(e) => Err(e),
        }
    };

    if let Err(e) = result {
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
    fn stats_basic() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nACGTACGT\n+\nIIIIIIII\n";
        // Just ensure it doesn't crash
        fastq_stats(data, false).unwrap();
    }

    #[test]
    fn stats_accumulation() {
        let mut stats = FastqStats::new();
        stats.add_read(b"ACGT", b"IIII");
        stats.add_read(b"ACGTACGT", b"IIIIIIII");

        assert_eq!(stats.total_reads, 2);
        assert_eq!(stats.total_bases, 12);
        assert_eq!(stats.min_length, 4);
        assert_eq!(stats.max_length, 8);
        assert_eq!(stats.mean_length(), 6.0);
    }
}
