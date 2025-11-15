//! FASTA sequence sampling utility

use std::io::{self, Write};
use std::process;

use clap::Parser;

mod faster_fasta;
use faster_fasta::*;

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

/// Sample random sequences from FASTA files using reservoir sampling
///
/// Uses reservoir sampling algorithm to avoid materializing all entries.
/// For count-based sampling, uses Algorithm R. For fraction-based sampling,
/// first counts total sequences, then uses reservoir sampling with calculated size.
pub fn fasta_sample(
    input: &[u8],
    output: impl Write,
    count: Option<usize>,
    fraction: Option<f64>,
    seed: Option<u64>,
) -> io::Result<()> {
    let mut rng = SimpleRng::new(seed.unwrap_or(0));

    if let Some(n) = count {
        // Reservoir sampling with fixed count
        reservoir_sample_count(input, output, n, &mut rng)
    } else if let Some(f) = fraction {
        // First pass: count total sequences
        let total = FastaParser::new(input).count();
        let sample_size = ((total as f64 * f).round() as usize).min(total);

        // Second pass: reservoir sampling with calculated size
        reservoir_sample_count(input, output, sample_size, &mut rng)
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
) -> io::Result<()> {
    if k == 0 {
        return Ok(());
    }

    let mut reservoir: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(k);
    let parser = FastaParser::new(input);

    for (i, entry) in parser.enumerate() {
        if i < k {
            // Fill reservoir
            reservoir.push((
                entry.header.to_vec(),
                entry.sequence.as_bytes().to_vec(),
            ));
        } else {
            // Randomly replace elements with decreasing probability
            let j = rng.gen_range((i + 1) as u64) as usize;
            if j < k {
                reservoir[j] = (
                    entry.header.to_vec(),
                    entry.sequence.as_bytes().to_vec(),
                );
            }
        }
    }

    // Write sampled entries
    for (header, sequence) in reservoir {
        output.write_all(&header)?;
        output.write_all(b"\n")?;
        output.write_all(&sequence)?;
        output.write_all(b"\n")?;
    }

    output.flush()?;
    Ok(())
}

/// Sample random sequences from FASTA files
#[derive(Parser)]
#[command(name = "fasta-sample")]
#[command(version, about, long_about = None)]
struct Args {
    /// Input FASTA file (use '-' or omit for stdin)
    input: Option<String>,

    /// Output FASTA file (use '-' or omit for stdout)
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
        eprintln!("Error processing FASTA: {}", e);
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
}
