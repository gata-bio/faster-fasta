#!/usr/bin/env python3
"""
Faster FASTA is a collection of command-line utilities for processing memory-mapped plain-text FASTA files.
It leverages StringZilla to provide high-performance functionality for:

- deduplication via `fasta-dedup`
- sorting sequences via `fasta-sort`
- filtering sequences by length via `fasta-filter-length`
- extracting sequences by ID via `fasta-extract-ids`
- merging FASTA files via `fasta-merge`
- reversing and complementing sequences via `fasta-revcomp`
- transliterating from DNA to RNA via `fasta-dna2rna`
- format line wrapping via `fasta-wrap`
- sample sequences randomly via `fasta-sample`
"""

import stringzilla as sz

class FASTA:
    
    def __init__(self, filepath):
        self._filepath = filepath
        self._file = sz.File(filepath)
        self._str = sz.Str(self._file)
        
    