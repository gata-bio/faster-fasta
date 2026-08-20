# FasterFASTA

![Faster FASTA Thumbnail](https://github.com/ashvardanian/ashvardanian/blob/master/repositories/FasterFASTA.jpg?raw=true)

__Faster FASTA__ is a collection of command-line utilities for processing FASTA and FASTQ files, memory-mapped or streamed from external storage or via `stdin`.
It's a faster SIMD-accelerated alternative to pure Go [`seqkit`](https://github.com/shenwei356/seqkit) and C++ [`fastp`](https://github.com/OpenGene/fastp) tools.
It's implemented in Rust with [StringZilla](https://github.com/ashvardanian/StringZilla) to provide high-performance functionality with __auto-format detection__, inspecting the header sigil — `@` for FASTQ, `>` for FASTA.

## Tools

Grouped by how much memory each tool needs, since that is what decides whether it survives a file larger than RAM.
Every tool auto-detects FASTA or FASTQ and preserves the format on output, unless noted.

Named files are memory-mapped, so resident size tracks the file even for the constant-memory tools.
Those pages are page-cache backed and reclaimable rather than an allocation, so a file larger than RAM is fine.

### Constant Memory — Safe on Any Input Size

Nothing is held between records, so these run on a 10 KB file and a 10 TB file alike, and they spread across cores with no coordination.

- `fastq-filter` — keep records passing quality, length, and N-content thresholds
- `fastq-trim` — trim by quality score, by fixed offsets at either end, and to a maximum length
- `fasta-sample --fraction` — Bernoulli sampling, deciding each record independently with no reservoir

Each tool is a single pass that buffers nothing between records, so a pipeline of them costs one parse per stage and no memory beyond the pipe.

### Bounded Memory — You Choose the Ceiling

Memory follows the flag rather than the input, so these stay usable where a full pass would not fit.

- `fasta-sample --count N` — uniform reservoir sample of N records, holding N records

### Whole-Input Memory — Bounded by the Data

- `fasta-dedup` — collapse duplicate records, matched by `--by sequence`, `--by name`, or `--by canonical`, the last being the lesser of a sequence and its reverse complement

Identity is always settled by comparing bytes rather than by a digest alone, and there is no fuzzy or edit-distance mode.
Memory follows the survivors rather than the input: each retained key is kept in full, plus 8 bytes for its end offset and __32 to 64 bytes of index__ depending on how full the open-addressed table is.
A billion distinct 150 bp reads is therefore roughly 180 GB, nearly all of it the sequences themselves.

Scanning goes wide with `--threads` — parsing, canonicalizing, hashing and rendering are pure functions of one record — while the index itself is folded serially, in input order, which is what makes first-occurrence mean the same thing at any worker count.
Splitting that index across workers inside one process would bound nothing, since the parts together hold what one table would.
Survivors that do not fit in RAM need the input split into partitions folded at different times or on different machines, and the results merged.

### Also Included

Small, exact, and unlikely to be your bottleneck, but here so a pipeline need not leave the toolkit.

- `fastq-stats` — length, quality, GC, and N-content summary, with an optional histogram
- `fastq-to-fasta` — drop quality scores; `--min-quality` filters first and is the only part that requires FASTQ input
- `fasta-revcomp` — reverse complement over the full IUPAC alphabet, reversing quality alongside sequence
- `fasta-dna2rna` — transcribe DNA to RNA by rewriting T as U

## Installation

```bash
cargo install --git https://github.com/unum-bio/FasterFASTA    # install from GitHub
cargo install --path . --force                                 # or install from local clone
```

## Usage

Remove duplicates:

```bash
fasta-dedup sequences.fasta --output unique.fasta
fasta-dedup reads.fastq --by name --output unique.fastq                    # match on the identifier
fasta-dedup a.fa b.fa c.fa --by canonical --threads 16 --output unique.fa  # either strand, 16 workers
```

Sample 1000 sequences:

```bash
fasta-sample sequences.fasta --count 1000 --output sample.fasta
```

Reverse complement, and transcribe to RNA:

```bash
fasta-revcomp sequences.fasta --output revcomp.fasta
fasta-dna2rna sequences.fasta --output rna.fasta
```

FASTQ quality filtering and trimming:

```bash
# Keep reads with mean Q≥25 and length ≥75
fastq-filter reads.fastq --min-quality 25 --min-length 75 --output filtered.fastq

# Trim low-quality tails and drop short reads
fastq-trim reads.fastq --trim-below-quality 20 --trim-end 5 --min-length 50 --output trimmed.fastq
```

FASTQ stats and format conversions:

```bash
fastq-stats reads.fastq --histogram | head        # quick QC summary
fastq-stats corpus/*.fastq --format json          # one JSON object per input
fastq-to-fasta reads.fastq --output reads.fasta   # drop qualities
```

A whole corpus at once, one output per input:

```bash
fastq-trim corpus/*.fastq --trim-end 5 --output-dir trimmed/   # keeps each input's file name
fastq-trim corpus/*.fastq --trim-end 5 --in-place              # rewrites each file where it lies
```

Count what a run would do before committing to it:

```bash
$ fastq-trim reads.fastq --trim-below-quality 20 --min-length 50 --dry-run
would trim 9708 of 10000 records, dropping 292 below the minimum length; nothing was written
```

All utilities support `stdin` and `stdout` for composability:

```bash
cat sequences.fasta | fasta-dedup | fasta-sample --fraction 0.1 > output.fasta
```

### Composing With StringZilla-CLI

Every tool here understands records — a header, a sequence, and for FASTQ a quality string.
For byte-level work on unstructured text, reach for [StringZilla-CLI](https://github.com/ashvardanian/StringZilla-CLI), whose `sz-*` tools operate on lines and bytes and share the same spelled-out flag vocabulary.

```bash
# Count and measure without parsing anything
sz-find --show count '>' assembly.fa                # 3
sz-count --fields lines,bytes assembly.fa           # 7 lines, 55 bytes

# Headers are plain lines, so sorting and deduplicating them is line work
sz-find --show lines '>' assembly.fa | sz-sort --unique

# One file per record, cutting at the header rather than at a line count
sz-split assembly.fa contig- --chunk-pattern '>'    # contig-aa, contig-ab, contig-ac

# Checksum a pipeline's output so a rerun proves itself
fastq-trim reads.fastq --trim-end 5 --output trimmed.fastq && sz-sha256 trimmed.fastq
```

The division is by what the bytes mean, and FASTQ shows why it matters.
A quality line may legitimately begin with `@`, which is also the header sigil, so counting lines is not counting records:

```bash
$ sz-find --show count '@' reads.fastq    # 3 lines carry an '@'
3
$ fastq-stats reads.fastq --format plain | cut -f3    # 2 records are there
2
```

Anything that has to respect a record boundary belongs here; anything that treats the file as lines belongs there.

### Command-Line Conventions

One vocabulary across all eight tools, so a flag learned once is a flag learned everywhere.

- __No single-letter flags.__
  Every option is spelled out, `-h` and `-V` aside, because `-n` meaning one thing here and another there is a bug waiting on a tired afternoon.
  Two tests per binary enforce it — one refuses a short flag, one pins the exact list of long ones — so a rename is a deliberate edit rather than a drift.
- __`--output`, `--output-dir`, `--in-place`, `--dry-run` and `--quiet`__ name where bytes go, and mean the same thing in every tool that emits records.
  `--output-dir` and `--in-place` are absent from `fasta-sample` and `fastq-stats` on purpose: a reservoir and a summary row both span every input, so there is no per-input file to write.
- __`--line-width`, `--color` and `--summary`__ decide how output looks and what the run says about itself afterwards.
  Wrapping and colour are detected from whether the destination is a terminal; these override that detection rather than replacing it.
- __A closed choice is one flag with named values__, never a pile of booleans: `--by sequence|name|canonical`, `--color auto|always|never`, `--format table|plain|json`.
- __Exit codes carry the answer__: 0 produced something, 1 ran and produced nothing, 2 could not run.

That last distinction lets a pipeline tell "widen the filter" apart from "fix the command line", and `--quiet` turns any tool into a test:

```bash
$ fastq-filter reads.fastq --min-quality 30 --quiet; echo $?
0
$ fastq-filter reads.fastq --min-quality 99 --quiet; echo $?
1
$ fastq-filter absent.fastq --quiet; echo $?
Error filtering records: No such file or directory (os error 2)
2
```

## Scaling

### Input Shapes

The container decides how much parallelism is available, independently of which tool runs.

| Input                       | Access              | Parallel decode | Decode speed |
| :-------------------------- | :------------------ | --------------: | -----------: |
| `.fa`, `.fq`                | byte-addressable    |             yes |     416 MB/s |
| `.bgz`, `bgzip` output      | block-addressable   |             yes |     207 MB/s |
| `.xz` written with `-T0`    | block-addressable   |             yes |      56 MB/s |
| `.zst`                      | forward-only stream |              no |     264 MB/s |
| `.gz`, plain gzip           | forward-only stream |              no |     222 MB/s |
| `.xz` written single-thread | forward-only stream |              no |      56 MB/s |
| `.bz2`                      | forward-only stream |              no |      25 MB/s |

Speeds are end-to-end `fastq-stats` throughput over a 194 MB FASTQ of 150 bp reads, measured on one core.
The codec is detected from the bytes rather than the file name, so a `bgzip` file called `.gz` is still read as blocks and a `.fq` that turns out to be gzip still opens.

A plain `.gz` is one continuous deflate stream and decodes on a single core whatever `--threads` says.
Recompressing with `bgzip` yields independently decodable 64 KB blocks and costs nothing in compression ratio.
`xz` is the opposite case: slow enough that a single core is the bottleneck, but threaded `xz` already writes a block index, so a large file compressed the way large files usually are can be decoded in parallel without anybody opting in.

### What a Worker Owns

`--threads` sets the worker count, and a worker owns one unit of work whose state fits in memory.
What that unit is follows from the input shape.

- A __byte-addressable__ file has no natural grain, so it is cut into byte ranges resynchronized to the next record boundary — which makes splitting one large file and walking a directory of files the same code path.
- A __block-addressable__ container is already divided, so a worker takes a run of blocks and never has to search for a boundary.
- A __forward-only stream__ cannot be divided at all, so it is read in slabs and only the parsing goes wide.

How those units recombine follows from the memory group:

| Tool                  | Combining partial results                                            |
| :-------------------- | :------------------------------------------------------------------- |
| Constant-memory tools | Outputs concatenate; there is nothing to merge                       |
| `fastq-stats`         | Partial summaries merge, in any order                                |
| `fasta-dedup`         | Scans wide and folds serially, so the bytes match at any `--threads` |

### Choosing a Block Size

__A block-addressable input cannot use more workers than it has blocks__, and that is usually the binding constraint rather than the core count.
Threaded `xz` splits at three times the dictionary — 24 MB at the default preset — so a 194 MB file holds nine blocks and `--threads 16` is slower than `--threads 8`, the extra workers only contending.
Cutting the same file into 98 blocks changes the picture entirely:

| The same 194 MB file, compressed as | 1 thread | 8 threads |   16 threads | Best speedup |
| :---------------------------------- | -------: | --------: | -----------: | -----------: |
| 9 blocks, `xz -T0` at its default   |  63 MB/s |  134 MB/s |     111 MB/s |         2.1× |
| 98 blocks, `xz --block-size=2MiB`   |  46 MB/s |  184 MB/s | __233 MB/s__ |     __5.1×__ |

So if you control how the archive is written and intend to read it back in parallel, set a block size rather than accepting the default:

```bash
xz -T0 --block-size=2MiB reads.fastq
```

## Performance

Benchmark on the shape of file you actually have.
Every comparison below runs on decompressed input, where the parser is the whole cost — but real archives arrive compressed, and then __the codec sets the ceiling and the parser is nearly free__.
A tool reading `.xz` spends roughly seven bytes of its time budget on decompression for every one on parsing, which is why the block tiers under Scaling matter more than any parsing win.

### Getting the Datasets

The UniProt Swiss-Prot database and the paired _Escherichia coli_ Illumina reads are the traditional pair to measure against.

```bash
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    gunzip uniprot_sprot.fasta.gz && \
    grep -c '^>' uniprot_sprot.fasta    # contains 573'661 sequences

curl -L -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/013/SRR25083113/SRR25083113_1.fastq.gz && \
    gunzip SRR25083113_1.fastq.gz && \
    grep -c '^@' SRR25083113_1.fastq    # contains 1'181'120 sequences

curl -L -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/013/SRR25083113/SRR25083113_2.fastq.gz && \
    gunzip SRR25083113_2.fastq.gz && \
    grep -c '^@' SRR25083113_2.fastq    # contains 1'181'120 sequences
```

### Against awk and seqkit

Deduplication, against the `awk` one-liner everybody reaches for first:

```bash
time fasta-dedup uniprot_sprot.fasta --output unique_faster.fasta
grep -c '^>' unique_faster.fasta    # prints 485'423 sequences after 0.4s

time awk '/^>/ {if (seq != "" && !seen[seq]++) {print header; print seq} header = $0; seq = ""; next} {seq = seq $0} END {if (seq != "" && !seen[seq]++) {print header;  print seq}}' uniprot_sprot.fasta > unique_awk.fasta
grep -c '^>' unique_awk.fasta       # prints 485'423 sequences after 11.3s
```

And against a full toolkit:

```bash
brew install seqkit hyperfine

# Deduplication: 0.4s vs 1.1s
hyperfine \
    'fasta-dedup uniprot_sprot.fasta --output /tmp/ff.fasta' \
    'seqkit rmdup -s uniprot_sprot.fasta -o /tmp/seqkit.fasta' --warmup 1

# Sampling (10% fraction): 0.17s vs 0.20s
hyperfine \
    'fasta-sample SRR25083113_1.fastq --fraction 0.1 --output /tmp/ff_sample.fastq' \
    'seqkit sample -p 0.1 SRR25083113_1.fastq -o /tmp/seqkit_sample.fastq' --warmup 1

# FASTQ stats
hyperfine \
    'fastq-stats SRR25083113_1.fastq' \
    'seqkit stats SRR25083113_1.fastq' --warmup 1
```

## Out of Scope

- __Read mapping against a reference index.__
  That is an FM-index and minimizer problem rather than a string-scanning one; use `bwa`, `minimap2`, or `bowtie2`.
- __SAM, BAM, VCF, and GFF.__
  Sequence files only; use `samtools` and `bcftools`.
- __Assembly and variant calling.__
- __Workflow orchestration.__
  These are composable single-purpose tools, meant to be driven by Nextflow, Snakemake, or a shell pipeline.
- __Splitting a plain `.gz` across cores.__
  `rapidgzip` and `pugz` decode an unindexed deflate stream in parallel; recompressing with `bgzip` is the cheaper answer here.
- __Sorting.__
  Assemblies small enough to sort are already instant with any tool, and the large-file workloads that look like sorting want the longest or first N records instead.
- __Paired-end interleaving.__
  Most aligners take R1 and R2 as separate arguments, and a byte range of R1 does not hold the mates of the same byte range of R2, so it is the one operation that cannot be divided across workers at all.
