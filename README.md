# Faster FASTA

Faster FASTA is a collection of command-line utilities for processing memory-mapped plain-text FASTA files.
It's implemented in Rust with StringZilla to provide high-performance functionality for:

- deduplication via `fasta-dedup`
- sample sequences randomly via `fasta-sample`
- sorting sequences via `fasta-sort`
- reversing and complementing sequences via `fasta-revcomp`
- transliterating from DNA to RNA via `fasta-dna2rna`
- merging FASTA files via `fasta-merge`
- format line wrapping via `fasta-wrap`
- filtering sequences by length via `fasta-filter-length`
- extracting sequences by ID via `fasta-extract-ids`

Deduplication is implemented using StringZilla's hash-function populating a flat hash-set, and checking for collisions of subsequent sequences with previously seen ones.
Transliteration and reverse-complementation are implemented using StringZilla's SIMD-accelerated byte-mapping functions.
Sorting is implemented using StringZilla's built-in sort for byte-strings.

## Quick Start

```bash
cargo install --git https://github.com/gata-bio/faster-fasta
```