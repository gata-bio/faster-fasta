![Faster FASTA Thumbnail](https://github.com/ashvardanian/ashvardanian/blob/master/repositories/faster-fasta.jpg?raw=true)

__Faster FASTA__ is a collection of command-line utilities for processing memory-mapped plain-text FASTA files.
It's implemented in Rust with [StringZilla](https://github.com/ashvardanian/stringZilla/) to provide high-performance functionality:

- `fasta-dedup` - remove duplicate sequences
- `fasta-sample` - randomly sample sequences
- `fasta-sort` - sort by name, sequence, or length
- `fasta-revcomp` - reverse complement DNA sequences
- `fasta-dna2rna` - convert DNA to RNA (T to U)

## Installation

```bash
cargo install --git https://github.com/ashvardanian/faster-fasta
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

All utilities support `stdin` and `stdout` for composability:

```bash
cat sequences.fasta | fasta-dedup | fasta-sort -l -r > output.fasta
```
