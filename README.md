![Faster FASTA Thumbnail](https://github.com/ashvardanian/ashvardanian/blob/master/repositories/faster-fasta.jpg?raw=true)

__Faster FASTA__ is a collection of command-line utilities for processing FASTA and FASTQ files, memory-mapped or streamed from external storage or via `stdin`.
It's a faster SIMD-accelerated alternative to pure Go [`seqkit`](https://github.com/shenwei356/seqkit) and C++ [`fastp`](https://github.com/OpenGene/fastp) tools.
It's implemented in Rust with [StringZilla](https://github.com/ashvardanian/StringZilla) to provide high-performance functionality with __auto-format detection__ (@ vs > header inspection).

## Tools Overview

Multi-Format tools for FASTA & FASTQ:

- `fasta-dedup` - remove duplicate sequences
- `fasta-sample` - randomly sample sequences using reservoir sampling
- `fasta-sort` - sort by name, sequence, or length
- `fasta-revcomp` - reverse complement DNA sequences (quality also reversed for FASTQ)
- `fasta-dna2rna` - convert DNA to RNA (T → U)

FASTQ-specific tools:

- `fastq-filter` - filter by quality, length, and N-content
- `fastq-trim` - quality-based and fixed-position trimming
- `fastq-stats` - comprehensive statistics with histograms
- `fastq-to-fasta` - format conversion (drop quality scores)
- `fastq-interleave` - merge paired-end files (R1 + R2 → interleaved)
- `fastq-deinterleave` - split interleaved file (interleaved → R1 + R2)

## Installation

```bash
cargo install --git https://github.com/gata-bio/faster-fasta    # install from GitHub
cargo install --path . --force                                  # or install from local clone
```

## Usage

Remove duplicates:

```bash
fasta-dedup sequences.fasta -o unique.fasta
```

Sample 1000 sequences:

```bash
fasta-sample sequences.fasta --count 1000 -o sample.fasta
```

Sort by length:

```bash
fasta-sort --length sequences.fasta -o sorted.fasta
```

Reverse complement:

```bash
fasta-revcomp sequences.fasta -o revcomp.fasta
```

Convert to RNA:

```bash
fasta-dna2rna sequences.fasta -o rna.fasta
```

FASTQ quality filtering and trimming:

```bash
# Keep reads with mean Q≥25 and length ≥75
fastq-filter reads.fastq --min-quality 25 --min-length 75 -o filtered.fastq

# Trim low-quality tails and drop short reads
fastq-trim reads.fastq --quality-cutoff 20 --trim-tail 5 --min-length 50 -o trimmed.fastq
```

FASTQ stats and format conversions:

```bash
fastq-stats reads.fastq --histogram | head      # quick QC summary
fastq-to-fasta reads.fastq -o reads.fasta       # drop qualities
```

Paired-end juggling:

```bash
fastq-interleave R1.fastq R2.fastq -o interleaved.fastq
fastq-deinterleave interleaved.fastq -1 out_R1.fastq -2 out_R2.fastq
```

All utilities support `stdin` and `stdout` for composability:

```bash
cat sequences.fasta | fasta-dedup | fasta-sort -l -r > output.fasta
```

## Performance

Consider pulling some traditional dataset, like the UniProt Swiss-Prot database and the paired Escherichia coli  (E. coli) Illumina reads, to benchmark performance.

```bash
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    gunzip uniprot_sprot.fasta.gz && \
    grep -c '^>' uniprot_sprot.fasta    # contains 573'661 sequences

curl -L -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/013/SRR25083113/SRR25083113_1.fastq.gz SRR25083113_1.fastq.gz && \
    gunzip SRR25083113_1.fastq.gz && \
    grep -c '^@' SRR25083113_1.fastq    # contains 1'181'120 sequences
    
curl -L -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/013/SRR25083113/SRR25083113_2.fastq.gz SRR25083113_2.fastq.gz && \
    gunzip SRR25083113_2.fastq.gz && \
    grep -c '^@' SRR25083113_2.fastq    # contains 1'181'120 sequences
```

Run following commands to compare the performance of `fasta-dedup` against a traditional `awk` approach for removing duplicate sequences:

```bash
time fasta-dedup uniprot_sprot.fasta -o unique_faster.fasta
grep -c '^>' unique_faster.fasta    # prints 485'423 sequences after 0.4s

time awk '/^>/ {if (seq != "" && !seen[seq]++) {print header; print seq} header = $0; seq = ""; next} {seq = seq $0} END {if (seq != "" && !seen[seq]++) {print header;  print seq}}' uniprot_sprot.fasta > unique_awk.fasta
grep -c '^>' unique_awk.fasta       # prints 485'423 sequences after 11.3s
```

You can also compare against a popular toolkit like `seqkit`:

```bash
brew install seqkit hyperfine

# Deduplication: 0.4s vs 1.1s
hyperfine \
    'fasta-dedup uniprot_sprot.fasta -o /tmp/ff.fasta' \
    'seqkit rmdup -s uniprot_sprot.fasta -o /tmp/seqkit.fasta' --warmup 1

# Sorting by length: 0.8s vs 2.9s
hyperfine \
    'fasta-sort --length SRR25083113_1.fastq -o /tmp/ff_sorted.fastq' \
    'seqkit sort -l SRR25083113_1.fastq -o /tmp/seqkit_sorted.fastq' --warmup 1

# Sampling (10% fraction): 0.17s vs 0.20s
hyperfine \
    'fasta-sample SRR25083113_1.fastq --fraction 0.1 -o /tmp/ff_sample.fastq' \
    'seqkit sample -p 0.1 SRR25083113_1.fastq -o /tmp/seqkit_sample.fastq' --warmup 1

# FASTQ stats
hyperfine \
    'fastq-stats SRR25083113_1.fastq' \
    'seqkit stats SRR25083113_1.fastq' --warmup 1
```
