# Faster FASTA

Faster FASTA is a collection of command-line utilities for processing memory-mapped plain-text FASTA files.
It leverages StringZilla to provide high-performance functionality for:

- deduplication via `fasta-dedup`
- sorting sequences via `fasta-sort`
- filtering sequences by length via `fasta-filter-length`
- extracting sequences by ID via `fasta-extract-ids`
- merging FASTA files via `fasta-merge`
- reversing and complementing sequences via `fasta-revcomp`
- transliterating from DNA to RNA via `fasta-dna2rna`
