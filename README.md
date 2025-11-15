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
cargo install --git https://github.com/ashvardanian/faster-fasta    # install from GitHub
cargo install --path .                                              # or install from local clone
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

## Performance

Consider pulling some traditional dataset, like the UniProt Swiss-Prot database:

```bash
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
grep -c '^>' unique_faster.fa   # contains 573'661 lines
```

Run following commands to compare the performance of `fasta-dedup` against a traditional `awk` approach for removing duplicate sequences:

```bash
time fasta-dedup uniprot_sprot.fasta -o unique_faster.fa
grep -c '^>' unique_faster.fa   # prints 485'423 lines after 0.4s

time awk '/^>/ {if (seq != "" && !seen[seq]++) {print header; print seq} header = $0; seq = ""; next} {seq = seq $0} END {if (seq != "" && !seen[seq]++) {print header;  print seq}}' uniprot_sprot.fasta > unique_awk.fa
grep -c '^>' unique_awk.fa      # prints 485'423 lines after 11.3s
```
