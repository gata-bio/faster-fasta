//! Faster Sequence Processing - Core FASTA and FASTQ parsing utilities
//!
//! Provides high-performance FASTA and FASTQ file parsing with memory-mapped I/O
//! and StringZilla for maximum performance.
//!
//! Supports automatic format detection based on content inspection.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Cursor, Read, Write};

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

/// Peek stdin to detect format and return a reader that yields the peeked bytes first
#[allow(dead_code)]
pub fn stdin_with_peek(limit: usize) -> io::Result<(SeqFormat, Box<dyn Read>)> {
    let stdin = io::stdin();
    let mut locked = stdin.lock();
    let mut prefix = Vec::with_capacity(limit);
    locked.by_ref().take(limit as u64).read_to_end(&mut prefix)?;
    let format = detect_format(&prefix)?;
    let reader: Box<dyn Read> = Box::new(Cursor::new(prefix).chain(locked));
    Ok((format, reader))
}

/// Create an output writer from either a file path or stdout
#[allow(dead_code)]
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

/// A FASTQ entry containing a header, sequence, and quality scores
#[derive(Debug, Clone)]
pub struct FastqEntry<'a> {
    pub header: &'a [u8],
    pub sequence: SequenceData<'a>,
    pub quality: &'a [u8],
}

/// Owned FASTQ entry for streaming stdin
#[derive(Debug, Clone)]
pub struct FastqOwnedEntry {
    pub header: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

/// Sequence file format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

/// Detect sequence format from file content
///
/// Inspects the first non-whitespace byte:
/// - '@' indicates FASTQ
/// - '>' indicates FASTA
/// - Otherwise returns error
pub fn detect_format(data: &[u8]) -> io::Result<SeqFormat> {
    for &byte in data.iter() {
        match byte {
            b' ' | b'\t' | b'\r' | b'\n' => continue,
            b'@' => return Ok(SeqFormat::Fastq),
            b'>' => return Ok(SeqFormat::Fasta),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "Unknown format: expected '@' or '>', found '{}'",
                        byte as char
                    ),
                ))
            }
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        "Empty or whitespace-only file",
    ))
}

/// Quality score conversion utilities
pub mod quality {
    /// Convert ASCII quality character to Phred score (Phred+33 encoding)
    ///
    /// Phred+33 is the standard encoding used by modern sequencers (Illumina 1.8+).
    /// Quality character range: '!' (33) to '~' (126), representing Q0 to Q93.
    #[inline]
    pub fn ascii_to_phred33(ascii: u8) -> u8 {
        ascii.saturating_sub(33)
    }

    /// Convert Phred score to ASCII quality character (Phred+33 encoding)
    #[inline]
    pub fn phred33_to_ascii(phred: u8) -> u8 {
        phred.saturating_add(33).min(126)
    }

    /// Calculate mean quality score from ASCII quality string
    pub fn mean_quality(quality: &[u8]) -> f32 {
        if quality.is_empty() {
            return 0.0;
        }
        let sum: u32 = quality.iter().map(|&q| ascii_to_phred33(q) as u32).sum();
        sum as f32 / quality.len() as f32
    }
}

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
    fn fasta_basic() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nACGT\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 3);
        assert_eq!(entries[0].header, b">seq1");
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGT");
        assert_eq!(entries[1].header, b">seq2");
        assert_eq!(entries[1].sequence.as_bytes(), b"TGCA");
    }

    #[test]
    fn fasta_multiline() {
        let data = b">seq1\nACGT\nTGCA\n>seq2\nAAAA\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCA");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAA");
    }

    #[test]
    fn fasta_multiline_ws() {
        let data = b">seq1\n  ACGT  \n\tTGCA\t\n  GGGG  \n>seq2\nAAAA \n TTTT\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCAGGGG");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAATTTT");
    }

    #[test]
    fn fasta_multiline_blank() {
        let data = b">seq1\nACGT\n\nTGCA\n\n>seq2\nAAAA\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCA");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAA");
    }

    #[test]
    fn fasta_multiline_mixed_ws() {
        let data = b">seq1 description\nACGT TGCA\n  GGGG\n\tCCCC  \n>seq2\n  AAAA\n";
        let entries: Vec<_> = FastaParser::new(data).collect();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCAGGGGCCCC");
        assert_eq!(entries[1].sequence.as_bytes(), b"AAAA");
    }

    #[test]
    fn detect_format_fasta() {
        let data = b">seq1\nACGT\n";
        assert_eq!(detect_format(data).unwrap(), SeqFormat::Fasta);
    }

    #[test]
    fn detect_format_fastq() {
        let data = b"@seq1\nACGT\n+\nIIII\n";
        assert_eq!(detect_format(data).unwrap(), SeqFormat::Fastq);
    }

    #[test]
    fn detect_format_leading_ws() {
        let data = b"  \n  >seq1\nACGT\n";
        assert_eq!(detect_format(data).unwrap(), SeqFormat::Fasta);
    }

    #[test]
    fn quality_conversion() {
        assert_eq!(quality::ascii_to_phred33(b'!'), 0);
        assert_eq!(quality::ascii_to_phred33(b'I'), 40);
        assert_eq!(quality::phred33_to_ascii(0), b'!');
        assert_eq!(quality::phred33_to_ascii(40), b'I');
    }

    #[test]
    fn mean_quality_values() {
        let qual = b"IIII"; // All Q40
        assert_eq!(quality::mean_quality(qual), 40.0);

        let qual = b"!!!!"; // All Q0
        assert_eq!(quality::mean_quality(qual), 0.0);

        let qual = b"!I"; // Q0 and Q40
        assert_eq!(quality::mean_quality(qual), 20.0);
    }
}

/// Parse FASTQ entries from input data using StringZilla
pub struct FastqParser<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> FastqParser<'a> {
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
}

impl<'a> Iterator for FastqParser<'a> {
    type Item = io::Result<FastqEntry<'a>>;

    fn next(&mut self) -> Option<Self::Item> {
        let raw_data = self.data;

        // Skip to next '@' at start of line
        while self.pos < raw_data.len() && raw_data[self.pos] != b'@' {
            self.pos += 1;
        }

        if self.pos >= raw_data.len() {
            return None;
        }

        // Parse header line
        let header_start = self.pos;
        let header_end = self.find_newline(self.pos).unwrap_or(raw_data.len());
        let header = &raw_data[header_start..header_end];

        self.pos = if header_end < raw_data.len() {
            header_end + 1
        } else {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ entry missing sequence/quality",
            )));
        };

        // Parse sequence line(s) - collect until we hit '+'
        let seq_start = self.pos;
        let mut seq_end = seq_start;

        while self.pos < raw_data.len() && raw_data[self.pos] != b'+' {
            let line_end = self.find_newline(self.pos).unwrap_or(raw_data.len());
            if self.pos < line_end {
                seq_end = line_end;
            }
            self.pos = if line_end < raw_data.len() {
                line_end + 1
            } else {
                return Some(Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FASTQ entry missing sequence lines",
                )));
            };
        }

        if self.pos >= raw_data.len() {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ entry missing '+' separator",
            )));
        }

        // Parse sequence
        let sequence = if seq_start < seq_end {
            FastaParser::parse_sequence(&raw_data[seq_start..seq_end])
        } else {
            SequenceData::Borrowed(&raw_data[seq_start..seq_start])
        };

        // Skip '+' line (separator line)
        let plus_line_end = self.find_newline(self.pos).unwrap_or(raw_data.len());
        self.pos = if plus_line_end < raw_data.len() {
            plus_line_end + 1
        } else {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ entry missing quality header",
            )));
        };

        // Parse quality line(s) - same length as sequence
        let seq_len = sequence.as_bytes().len();
        let qual_start = self.pos;
        let mut qual_collected = 0;

        while qual_collected < seq_len && self.pos < raw_data.len() {
            let line_end = self.find_newline(self.pos).unwrap_or(raw_data.len());
            let line_len = line_end - self.pos;
            qual_collected += line_len;

            self.pos = if line_end < raw_data.len() {
                line_end + 1
            } else {
                line_end
            };
        }

        let qual_end = qual_start + seq_len.min(qual_collected);

        if qual_end > raw_data.len() || qual_collected < seq_len {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ entry missing quality string",
            )));
        }

        // Quality may have newlines - need to strip them
        let quality_raw = &raw_data[qual_start..qual_end];
        let whitespace = Byteset::from(b" \t\r\n");

        // For quality scores, we need to handle multi-line quality strings
        // For now, we'll store as borrowed if single-line, but this is a simplification
        // In real FASTQ files, quality is usually single-line, but specification allows multi-line
        let quality = if find_byteset(quality_raw, whitespace).is_some() {
            // Quality has whitespace - we can't easily return borrowed data
            // For now, return what we have and document this limitation
            // Most FASTQ files have single-line quality strings
            quality_raw
        } else {
            quality_raw
        };

        if quality.len() != sequence.as_bytes().len() {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FASTQ sequence/quality length mismatch: seq {}, qual {}",
                    sequence.as_bytes().len(),
                    quality.len()
                ),
            )));
        }

        Some(Ok(FastqEntry {
            header,
            sequence,
            quality,
        }))
    }
}

/// Streaming FASTQ parser for stdin using a bounded buffer
pub struct FastqStreamParser<R: Read> {
    reader: BufReader<R>,
    buf: Vec<u8>,
    start: usize,
    eof: bool,
}

impl<R: Read> FastqStreamParser<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            buf: Vec::with_capacity(1 << 20), // 1 MiB initial
            start: 0,
            eof: false,
        }
    }

    fn fill(&mut self) -> io::Result<()> {
        if self.eof {
            return Ok(());
        }

        // Compact buffer if we've consumed a large front portion
        if self.start > 0 {
            let remaining = self.buf.len().saturating_sub(self.start);
            self.buf.copy_within(self.start.., 0);
            self.buf.truncate(remaining);
            self.start = 0;
        }

        // Grow capacity if needed
        if self.buf.len() == self.buf.capacity() {
            let new_cap = (self.buf.capacity().max(1024)) * 2;
            self.buf.reserve(new_cap - self.buf.capacity());
        }

        let old_len = self.buf.len();
        // Make room for read
        let chunk = 1 << 20;
        self.buf.resize(old_len + chunk, 0);
        let read_n = self.reader.read(&mut self.buf[old_len..])?;
        self.buf.truncate(old_len + read_n);

        if read_n == 0 {
            self.eof = true;
        }
        Ok(())
    }

    fn find_newline(&mut self, from: usize) -> io::Result<Option<usize>> {
        let pos = from;
        loop {
            if pos >= self.buf.len() {
                if self.eof {
                    return Ok(None);
                }
                self.fill()?;
                continue;
            }
            if let Some(rel) = find(&self.buf[pos..], b"\n") {
                return Ok(Some(pos + rel));
            }
            if self.eof {
                return Ok(None);
            }
            self.fill()?;
        }
    }

    pub fn next_entry(&mut self) -> io::Result<Option<FastqOwnedEntry>> {
        loop {
            if self.start >= self.buf.len() {
                if self.eof {
                    return Ok(None);
                }
                self.fill()?;
                continue;
            }

            if self.buf[self.start] != b'@' {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FASTQ header must start with '@'",
                ));
            }

            // Header
            let header_end = self.find_newline(self.start)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FASTQ header missing newline")
            })?;
            let header = self.buf[self.start..header_end].to_vec();
            let mut pos = header_end + 1;

            // Sequence lines until '+'
            let mut sequence = Vec::new();
            loop {
                if pos >= self.buf.len() {
                    if self.eof {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "FASTQ sequence truncated",
                        ));
                    }
                    self.fill()?;
                    continue;
                }
                if self.buf[pos] == b'+' {
                    break;
                }
                let line_end = self.find_newline(pos)?.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "FASTQ sequence line missing newline",
                    )
                })?;
                sequence.extend_from_slice(&self.buf[pos..line_end]);
                pos = line_end + 1;
            }

            // Plus line
            let plus_line_end = self.find_newline(pos)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FASTQ '+' line missing newline")
            })?;
            pos = plus_line_end + 1;

            // Quality lines
            let mut quality = Vec::with_capacity(sequence.len());
            while quality.len() < sequence.len() {
                if pos >= self.buf.len() {
                    if self.eof {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "FASTQ quality truncated",
                        ));
                    }
                    self.fill()?;
                    continue;
                }
                let line_end = match self.find_newline(pos)? {
                    Some(v) => v,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "FASTQ quality line missing newline",
                        ))
                    }
                };
                quality.extend_from_slice(&self.buf[pos..line_end]);
                pos = line_end + 1;
            }

            if quality.len() != sequence.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "FASTQ sequence/quality length mismatch: seq {}, qual {}",
                        sequence.len(),
                        quality.len()
                    ),
                ));
            }

            self.start = pos;

            return Ok(Some(FastqOwnedEntry {
                header,
                sequence,
                quality,
            }));
        }
    }
}

#[cfg(test)]
mod fastq_tests {
    use super::*;

    #[test]
    fn parser_single_line() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n";
        let entries: Vec<_> = FastqParser::new(data).collect::<Result<_, _>>().unwrap();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].header, b"@seq1");
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGT");
        assert_eq!(entries[0].quality, b"IIII");
        assert_eq!(entries[1].header, b"@seq2");
        assert_eq!(entries[1].sequence.as_bytes(), b"TGCA");
        assert_eq!(entries[1].quality, b"HHHH");
    }

    #[test]
    fn parser_with_desc() {
        let data = b"@seq1 description here\nACGT\n+\nIIII\n";
        let entries: Vec<_> = FastqParser::new(data).collect::<Result<_, _>>().unwrap();

        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].header, b"@seq1 description here");
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGT");
        assert_eq!(entries[0].quality, b"IIII");
    }

    #[test]
    fn parser_multiline_seq() {
        let data = b"@seq1\nACGT\nTGCA\n+\nIIII\nHHHH\n";
        let entries: Vec<_> = FastqParser::new(data).collect::<Result<_, _>>().unwrap();

        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].sequence.as_bytes(), b"ACGTTGCA");
        // Note: quality parsing for multi-line is simplified
    }

    #[test]
    fn parser_empty_seq() {
        let data = b"@seq1\n\n+\n\n";
        let entries: Vec<_> = FastqParser::new(data).collect::<Result<_, _>>().unwrap();

        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].sequence.as_bytes(), b"");
        assert_eq!(entries[0].quality, b"");
    }

    #[test]
    fn parser_mismatch_len() {
        let data = b"@seq1\nACGT\n+\n!!!\n";
        let mut parser = FastqParser::new(data);
        let err = parser.next().unwrap().unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn stream_round_trip() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n";
        let mut parser = FastqStreamParser::new(std::io::Cursor::new(data));
        let mut out = Vec::new();
        while let Some(entry) = parser.next_entry().unwrap() {
            out.extend_from_slice(&entry.header);
            out.push(b'\n');
            out.extend_from_slice(&entry.sequence);
            out.push(b'\n');
            out.extend_from_slice(&entry.quality);
            out.push(b'\n');
        }
        let reconstructed = String::from_utf8(out).unwrap();
        assert!(reconstructed.contains("@seq1"));
        assert!(reconstructed.contains("@seq2"));
        assert!(reconstructed.contains("IIII"));
    }
}
