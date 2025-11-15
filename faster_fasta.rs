//! Faster FASTA - Core FASTA parsing utilities
//!
//! Provides high-performance FASTA file parsing with memory-mapped I/O
//! and StringZilla for maximum performance.

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};

use memmap2::Mmap;
use stringzilla::sz::{copy, find, find_byteset, hash, Byteset};

/// Represents the input source - either a memory-mapped file or buffered stdin
pub enum InputSource {
    /// Memory-mapped file for zero-copy access
    MappedFile(Mmap),
    /// Buffered stdin data
    Buffer(Vec<u8>),
}

impl InputSource {
    /// Get the input data as a byte slice
    pub fn as_bytes(&self) -> &[u8] {
        match self {
            InputSource::MappedFile(mmap) => &mmap[..],
            InputSource::Buffer(buf) => buf,
        }
    }
}

/// Create an InputSource from either a file path or stdin
pub fn get_input(path: Option<&str>) -> io::Result<InputSource> {
    match path {
        None | Some("-") => {
            let mut buffer = Vec::new();
            io::stdin().read_to_end(&mut buffer)?;
            Ok(InputSource::Buffer(buffer))
        }
        Some(path) => {
            let file = File::open(path)?;
            let mmap = unsafe { Mmap::map(&file)? };
            Ok(InputSource::MappedFile(mmap))
        }
    }
}

/// Create an output writer from either a file path or stdout
pub fn get_output(path: Option<&str>) -> io::Result<Box<dyn Write>> {
    match path {
        None | Some("-") => Ok(Box::new(BufWriter::new(io::stdout()))),
        Some(path) => {
            let file = File::create(path)?;
            Ok(Box::new(BufWriter::new(file)))
        }
    }
}

/// A FASTA entry containing a header and sequence
#[derive(Debug, Clone)]
pub struct FastaEntry<'a> {
    pub header: &'a [u8],
    pub sequence: SequenceData<'a>,
}

/// Represents sequence data that can be either borrowed or owned
#[derive(Debug, Clone)]
pub enum SequenceData<'a> {
    /// Borrowed slice from original data (clean, single-line sequence)
    Borrowed(&'a [u8]),
    /// Owned vector (sequence with whitespace removed)
    Owned(Vec<u8>),
}

impl<'a> SequenceData<'a> {
    pub fn as_bytes(&self) -> &[u8] {
        match self {
            SequenceData::Borrowed(s) => s,
            SequenceData::Owned(v) => v,
        }
    }
}

impl<'a> std::hash::Hash for SequenceData<'a> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let bytes = self.as_bytes();
        state.write_u64(hash(bytes));
    }
}

impl<'a> PartialEq for SequenceData<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl<'a> Eq for SequenceData<'a> {}

/// Parse FASTA entries from input data using StringZilla
pub struct FastaParser<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> FastaParser<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    #[inline]
    fn find_newline(&self, start: usize) -> Option<usize> {
        if start >= self.data.len() {
            return None;
        }
        let remaining = &self.data[start..];
        find(remaining, b"\n").map(|i| start + i)
    }

    fn parse_sequence(seq_region: &'a [u8]) -> SequenceData<'a> {
        let whitespace = Byteset::from(b" \t\r\n");

        if find_byteset(seq_region, whitespace).is_some() {
            // Has whitespace: strip using StringZilla byteset operations
            let mut seq = Vec::with_capacity(seq_region.len());

            let mut pos = 0;
            while pos < seq_region.len() {
                let next_ws = find_byteset(&seq_region[pos..], whitespace)
                    .map(|i| pos + i)
                    .unwrap_or(seq_region.len());

                if next_ws > pos {
                    let chunk = &seq_region[pos..next_ws];
                    let old_len = seq.len();
                    seq.resize(old_len + chunk.len(), 0);
                    copy(&mut seq[old_len..], chunk);
                }

                pos = next_ws;
                if pos < seq_region.len() {
                    pos += 1;
                }
            }

            SequenceData::Owned(seq)
        } else {
            SequenceData::Borrowed(seq_region)
        }
    }
}

impl<'a> Iterator for FastaParser<'a> {
    type Item = FastaEntry<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let raw_data = self.data;

        // Skip to next '>'
        while self.pos < raw_data.len() && raw_data[self.pos] != b'>' {
            self.pos += 1;
        }

        if self.pos >= raw_data.len() {
            return None;
        }

        // Parse header
        let header_start = self.pos;
        let header_end = self.find_newline(self.pos).unwrap_or(raw_data.len());
        let header = &raw_data[header_start..header_end];

        self.pos = if header_end < raw_data.len() {
            header_end + 1
        } else {
            header_end
        };

        // Find sequence region
        let seq_start = self.pos;
        let mut seq_end = seq_start;

        while self.pos < raw_data.len() && raw_data[self.pos] != b'>' {
            let line_end = self.find_newline(self.pos).unwrap_or(raw_data.len());
            if self.pos < line_end {
                seq_end = line_end;
            }
            self.pos = if line_end < raw_data.len() {
                line_end + 1
            } else {
                line_end
            };
        }

        // Parse sequence
        let sequence = if seq_start < seq_end {
            Self::parse_sequence(&raw_data[seq_start..seq_end])
        } else {
            SequenceData::Borrowed(&raw_data[seq_start..seq_start])
        };

        Some(FastaEntry { header, sequence })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_parser() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nACGT\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 3);
        assert_eq!(entries[0].header, b">seq1");
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGT");
        assert_eq!(entries[1].header, b">seq2");
        assert_eq!(entries[1].sequence.as_bytes(), b"TGCA");
    }

    #[test]
    fn test_multiline_sequences() {
        let data = b">seq1\nACGT\nTGCA\n>seq2\nAAAA\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCA");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAA");
    }

    #[test]
    fn test_multiline_with_whitespace() {
        let data = b">seq1\n  ACGT  \n\tTGCA\t\n  GGGG  \n>seq2\nAAAA \n TTTT\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCAGGGG");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAATTTT");
    }

    #[test]
    fn test_multiline_with_blank_lines() {
        let data = b">seq1\nACGT\n\nTGCA\n\n>seq2\nAAAA\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCA");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAA");
    }

    #[test]
    fn test_multiline_with_mixed_whitespace() {
        let data = b">seq1 description\nACGT TGCA\n  GGGG\n\tCCCC  \n>seq2\n  AAAA\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCAGGGGCCCC");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAA");
    }
}
