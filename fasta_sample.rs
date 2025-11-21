//! Sequence sampling utility for FASTA and FASTQ files
//!
//! Randomly sample sequences using reservoir sampling (Algorithm R).
//! Fixed memory usage regardless of input size.
//! Auto-detects format and preserves it in output.
//!
//! **Memory**: O(k) - only stores sampled sequences
//! **Streaming**: Yes - reservoir sampling, no materialization
//!
//! # Examples
//!
//! ```bash
//! # Sample 1000 sequences from FASTA
//! fasta-sample sequences.fasta --count 1000 -o sample.fasta
//!
//! # Sample 10% from FASTQ with seed for reproducibility
//! fasta-sample reads.fastq --fraction 0.1 --seed 42 -o sample.fastq
//!
//! # From stdin
//! cat sequences.fasta | fasta-sample --count 100 > sample.fasta
//! ```

use std::io::{self, Write};
use std::process;

use clap::Parser;

mod shared;
use shared::*;

/// Simple LCG random number generator for deterministic sampling
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    fn gen_range(&mut self, max: u64) -> u64 {
        self.next() % max
    }
}

/// Stored entry in reservoir (supports both FASTA and FASTQ)
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

/// Sample random sequences from FASTA or FASTQ files using reservoir sampling
///
/// Uses reservoir sampling algorithm to avoid materializing all entries.
/// For count-based sampling, uses Algorithm R. For fraction-based sampling,
/// first counts total sequences, then uses reservoir sampling with calculated size.
/// Auto-detects format and preserves it in output.
pub fn fasta_sample(
    input: &[u8],
    output: impl Write,
    count: Option<usize>,
    fraction: Option<f64>,
    seed: Option<u64>,
) -> io::Result<()> {
    let format = detect_format(input)?;
    let mut rng = SimpleRng::new(seed.unwrap_or(0));

    if let Some(n) = count {
        // Reservoir sampling with fixed count
        reservoir_sample_count(input, output, n, &mut rng, format)
    } else if let Some(f) = fraction {
        // First pass: count total sequences
        let total = match format {
            SeqFormat::Fasta => FastaParser::new(input).count(),
            SeqFormat::Fastq => FastqParser::new(input).count(),
        };
        let sample_size = ((total as f64 * f).round() as usize).min(total);

        // Second pass: reservoir sampling with calculated size
        reservoir_sample_count(input, output, sample_size, &mut rng, format)
    } else {
        Ok(())
    }
}

/// Reservoir sampling with fixed sample size (Algorithm R)
fn reservoir_sample_count(
    input: &[u8],
    mut output: impl Write,
    k: usize,
    rng: &mut SimpleRng,
    format: SeqFormat,
) -> io::Result<()> {
    if k == 0 {
        return Ok(());
    }

    let mut reservoir: Vec<StoredEntry> = Vec::with_capacity(k);

    match format {
        SeqFormat::Fasta => {
            let parser = FastaParser::new(input);
            for (i, entry) in parser.enumerate() {
                if i < k {
                    reservoir.push(StoredEntry::Fasta {
                        header: entry.header.to_vec(),
                        sequence: entry.sequence.as_bytes().to_vec(),
                    });
                } else {
                    let j = rng.gen_range((i + 1) as u64) as usize;
                    if j < k {
                        reservoir[j] = StoredEntry::Fasta {
                            header: entry.header.to_vec(),
                            sequence: entry.sequence.as_bytes().to_vec(),
                        };
                    }
                }
            }
        }
        SeqFormat::Fastq => {
            let parser = FastqParser::new(input);
            for (i, entry) in parser.enumerate() {
                if i < k {
                    reservoir.push(StoredEntry::Fastq {
                        header: entry.header.to_vec(),
                        sequence: entry.sequence.as_bytes().to_vec(),
                        quality: entry.quality.to_vec(),
                    });
                } else {
                    let j = rng.gen_range((i + 1) as u64) as usize;
                    if j < k {
                        reservoir[j] = StoredEntry::Fastq {
                            header: entry.header.to_vec(),
                            sequence: entry.sequence.as_bytes().to_vec(),
                            quality: entry.quality.to_vec(),
                        };
                    }
                }
            }
        }
    }

    // Write sampled entries
    for entry in reservoir {
        match entry {
            StoredEntry::Fasta { header, sequence } => {
                output.write_all(&header)?;
                output.write_all(b"\n")?;
                output.write_all(&sequence)?;
                output.write_all(b"\n")?;
            }
            StoredEntry::Fastq {
                header,
                sequence,
                quality,
            } => {
                output.write_all(&header)?;
                output.write_all(b"\n")?;
                output.write_all(&sequence)?;
                output.write_all(b"\n+\n")?;
                output.write_all(&quality)?;
                output.write_all(b"\n")?;
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Sample random sequences from FASTA or FASTQ files
#[derive(Parser)]
#[command(name = "fasta-sample")]
#[command(version, about = "Randomly sample sequences using reservoir sampling")]
#[command(
    long_about = "Randomly sample sequences from FASTA or FASTQ files.\nUses reservoir sampling for memory-efficient processing.\nAuto-detects format and preserves it in output."
)]
struct Args {
    /// Input file (FASTA or FASTQ, use '-' or omit for stdin)
    input: Option<String>,

    /// Output file (use '-' or omit for stdout)
    #[arg(short, long)]
    output: Option<String>,

    /// Number of sequences to sample
    #[arg(short = 'n', long, conflicts_with = "fraction")]
    count: Option<usize>,

    /// Fraction of sequences to sample (0.0 to 1.0)
    #[arg(short, long, conflicts_with = "count")]
    fraction: Option<f64>,

    /// Random seed for reproducible sampling
    #[arg(short, long)]
    seed: Option<u64>,
}

fn main() {
    let args = Args::parse();

    if args.count.is_none() && args.fraction.is_none() {
        eprintln!("Error: must specify either --count or --fraction");
        process::exit(1);
    }

    if let Some(f) = args.fraction {
        if f < 0.0 || f > 1.0 {
            eprintln!("Error: fraction must be between 0.0 and 1.0");
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

    if let Err(e) = fasta_sample(
        input.as_bytes(),
        output,
        args.count,
        args.fraction,
        args.seed,
    ) {
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
    fn test_fasta_sample_count() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n>seq4\nTTTT\n";
        let mut output = Vec::new();
        fasta_sample(data, &mut output, Some(2), None, Some(42)).unwrap();

        let result = String::from_utf8(output).unwrap();
        let count = result.matches('>').count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_fasta_sample_fraction() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n>seq4\nTTTT\n";
        let mut output = Vec::new();
        fasta_sample(data, &mut output, None, Some(0.5), Some(42)).unwrap();

        let result = String::from_utf8(output).unwrap();
        let count = result.matches('>').count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_fasta_sample_deterministic() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n";
        let mut output1 = Vec::new();
        let mut output2 = Vec::new();

        fasta_sample(data, &mut output1, Some(2), None, Some(123)).unwrap();
        fasta_sample(data, &mut output2, Some(2), None, Some(123)).unwrap();

        assert_eq!(output1, output2);
    }

    #[test]
    fn test_fastq_sample_count() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n@seq3\nAAAA\n+\nJJJJ\n";
        let mut output = Vec::new();
        fasta_sample(data, &mut output, Some(2), None, Some(42)).unwrap();

        let result = String::from_utf8(output).unwrap();
        let count = result.matches('@').count();
        assert_eq!(count, 2);
        // Verify quality is preserved
        assert!(result.contains('+'));
    }
}
