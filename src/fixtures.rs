//! The inputs every test module reads, as one space rather than many constants.
//!
//! A hand-rolled fixture pins every axis at the one value its author had in mind, and an axis
//! pinned at one value is an axis the suite does not test.
//! [`Shape`] names the axes, [`ALL`] names the points the suite spans, and [`records`] turns a
//! point into bytes.
//!
//! Nothing here is production code: the module is declared under `cfg(test)` and compiled only
//! into the test binary.

// A fixture serves whichever tests are compiled, and turning a codec feature off leaves the
// builders that feed it with no caller. That is a feature combination, not an unused function.
#![allow(dead_code)]

use std::io::{self, Read};

use crate::records::SequenceFormat;

// region: The Space

/// One point in the space of inputs the parser has to survive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Shape {
    /// Which of the two record grammars the bytes follow.
    pub(crate) format: SequenceFormat,
    /// Bases in the narrowest record; widths cycle upward from here so records stay ragged.
    pub(crate) read_bases: usize,
    /// Columns one sequence line holds, or zero for one line per sequence.
    pub(crate) wrap_columns: usize,
    /// The bytes that end every line, `"\n"` or `"\r\n"`.
    pub(crate) line_ending: &'static str,
    /// Whether the last line of the last record carries its ending.
    pub(crate) final_newline: bool,
    /// Whether a FASTQ quality line opens with `@`, which is Phred 31 rather than a header.
    pub(crate) quality_opens_with_sigil: bool,
}

/// Short unwrapped reads, the everyday shape of a short-read FASTQ.
pub(crate) const SHORT_FASTQ: Shape = Shape {
    format: SequenceFormat::Fastq,
    read_bases: 8,
    wrap_columns: 0,
    line_ending: "\n",
    final_newline: true,
    quality_opens_with_sigil: false,
};

/// Short wrapped FASTA, so nearly every record is joined out of several lines and reaches the
/// scratch buffer rather than being borrowed whole.
pub(crate) const WRAPPED_FASTA: Shape = Shape {
    format: SequenceFormat::Fasta,
    read_bases: 9,
    wrap_columns: 6,
    line_ending: "\n",
    final_newline: true,
    quality_opens_with_sigil: false,
};

/// CRLF line endings, which is what a file that travelled through Windows carries.
pub(crate) const CARRIAGE_RETURN_FASTQ: Shape = Shape {
    line_ending: "\r\n",
    ..SHORT_FASTQ
};

/// No final newline, which is what a cut-short download and a hand-edited file both look like,
/// and the one shape whose last record has no terminator to be found by.
pub(crate) const UNTERMINATED_FASTQ: Shape = Shape {
    final_newline: false,
    ..SHORT_FASTQ
};

/// A quality line opening with `@` at Phred 31, which only the four-line proof tells apart from
/// a header.
pub(crate) const SIGIL_QUALITY_FASTQ: Shape = Shape {
    quality_opens_with_sigil: true,
    ..SHORT_FASTQ
};

/// A 4 kb read, wide enough that one record fills a small parse window on its own.
pub(crate) const FOUR_KILOBASE_FASTQ: Shape = Shape {
    read_bases: 4096,
    ..SHORT_FASTQ
};

/// A 20 kb nanopore read, whose four-line proof spans tens of kilobytes and so cannot finish
/// inside one window.
///
/// This is the shape that dropped 41 of 2600 records while every short-read fixture passed.
pub(crate) const TWENTY_KILOBASE_FASTQ: Shape = Shape {
    read_bases: 20_000,
    ..SHORT_FASTQ
};

/// A long FASTA contig wrapped at sixty columns, the layout a reference genome ships in.
pub(crate) const LONG_WRAPPED_FASTA: Shape = Shape {
    read_bases: 20_000,
    wrap_columns: 60,
    ..WRAPPED_FASTA
};

/// Every shape the suite spans, so a consumer iterates rather than picking a value.
pub(crate) const ALL: [Shape; 8] = [
    SHORT_FASTQ,
    WRAPPED_FASTA,
    CARRIAGE_RETURN_FASTQ,
    UNTERMINATED_FASTQ,
    SIGIL_QUALITY_FASTQ,
    FOUR_KILOBASE_FASTQ,
    TWENTY_KILOBASE_FASTQ,
    LONG_WRAPPED_FASTA,
];

// endregion: The Space

// region: Sequence Fixtures

/// Appends `line` and the line ending `shape` asks for.
fn push_line(data: &mut Vec<u8>, line: &[u8], shape: Shape) {
    data.extend_from_slice(line);
    data.extend_from_slice(shape.line_ending.as_bytes());
}

/// Appends `payload` as one line, or as several of `shape.wrap_columns` bytes each.
fn push_wrapped(data: &mut Vec<u8>, payload: &[u8], shape: Shape) {
    match shape.wrap_columns {
        0 => push_line(data, payload, shape),
        columns => {
            for line in payload.chunks(columns) {
                push_line(data, line, shape);
            }
        }
    }
}

/// `count` records in the shape `shape` names.
///
/// Widths cycle over five values from [`Shape::read_bases`] upward, so a record straddles a
/// block or byte-range boundary far more often than it lands on one. [`uniform_records`] is the
/// fixture for the other case.
pub(crate) fn records(shape: Shape, count: usize) -> Vec<u8> {
    let mut data = Vec::new();
    for index in 0..count {
        let bases = shape.read_bases + (index % 5) * 4;
        let sequence = vec![b"ACGT"[index % 4]; bases];
        match shape.format {
            SequenceFormat::Fasta => {
                push_line(
                    &mut data,
                    format!(">contig{index} sample").as_bytes(),
                    shape,
                );
                push_wrapped(&mut data, &sequence, shape);
            }
            SequenceFormat::Fastq => {
                push_line(&mut data, format!("@read{index} sample").as_bytes(), shape);
                push_wrapped(&mut data, &sequence, shape);
                push_line(&mut data, b"+", shape);
                let mut quality = vec![b'I'; bases];
                if shape.quality_opens_with_sigil {
                    if let Some(first) = quality.first_mut() {
                        *first = b'@';
                    }
                }
                push_wrapped(&mut data, &quality, shape);
            }
        }
    }
    if !shape.final_newline {
        let ending = shape.line_ending.as_bytes();
        if data.ends_with(ending) {
            data.truncate(data.len() - ending.len());
        }
    }
    data
}

/// Roughly how many bytes one record of `shape` occupies, averaged over the width cycle.
fn bytes_per_record(shape: Shape) -> usize {
    let bases = shape.read_bases + 8;
    let ending = shape.line_ending.len();
    let sequence_lines = match shape.wrap_columns {
        0 => 1,
        columns => bases.div_ceil(columns).max(1),
    };
    let header = 24 + ending;
    match shape.format {
        SequenceFormat::Fasta => header + bases + sequence_lines * ending,
        SequenceFormat::Fastq => header + 2 * (bases + sequence_lines * ending) + 1 + ending,
    }
}

/// How many records of `shape` fit in about `target_bytes`, never fewer than two.
///
/// A test that spans the whole table would otherwise run a 20 kb shape at the record count a
/// short one wants, and spend minutes proving what a handful of records already proves.
pub(crate) fn count_within(shape: Shape, target_bytes: usize) -> usize {
    (target_bytes / bytes_per_record(shape)).max(2)
}

/// The bytes one record of [`uniform_records`] occupies, so a caller can size a block to hold a
/// whole number of them.
pub(crate) const UNIFORM_RECORD_BYTES: usize = 32;

/// FASTQ records of one fixed byte width, so a block boundary can land exactly on a record
/// start.
///
/// [`records`] never does, and a budget off by one still passes there. Every record is
/// [`UNIFORM_RECORD_BYTES`] wide, which the zero-padded index keeps true up to a million.
pub(crate) fn uniform_records(count: usize) -> Vec<u8> {
    let mut data = Vec::new();
    for index in 0..count {
        data.extend_from_slice(format!("@read{index:06}\nACGTACGT\n+\nIIIIIIII\n").as_bytes());
    }
    data
}

// endregion: Sequence Fixtures

// region: Readers

/// Returns at most `limit` bytes per call, the way a pipe or socket does.
///
/// `Cursor` always satisfies the full request and so can never drive a window through a
/// mid-record refill.
pub(crate) struct ShortReader<R: Read> {
    source: R,
    limit: usize,
}

impl<R: Read> ShortReader<R> {
    /// Wraps `source` so no single read returns more than `limit` bytes.
    pub(crate) fn new(source: R, limit: usize) -> Self {
        ShortReader { source, limit }
    }
}

impl<R: Read> Read for ShortReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let capped = out.len().min(self.limit);
        self.source.read(&mut out[..capped])
    }
}

// endregion: Readers

// region: Container Fixtures

/// One BGZF block over `payload`: gzip magic, a `BC` extra field, raw deflate, trailer.
#[cfg(feature = "gzip")]
pub(crate) fn bgzf_block(payload: &[u8]) -> Vec<u8> {
    let mut deflated = Vec::with_capacity(payload.len() + 64);
    let mut compress = flate2::Compress::new(flate2::Compression::fast(), false);
    compress
        .compress_vec(payload, &mut deflated, flate2::FlushCompress::Finish)
        .unwrap();
    let mut checksum = flate2::Crc::new();
    checksum.update(payload);

    let size = 12 + 6 + deflated.len() + 8;
    let mut block = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0];
    block.extend_from_slice(b"BC");
    block.extend_from_slice(&2u16.to_le_bytes());
    block.extend_from_slice(&((size - 1) as u16).to_le_bytes());
    block.extend_from_slice(&deflated);
    block.extend_from_slice(&checksum.sum().to_le_bytes());
    block.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    block
}

/// `plain` as BGZF in blocks of `payload_bytes`, small enough that records straddle them.
#[cfg(feature = "gzip")]
pub(crate) fn bgzf_bytes(plain: &[u8], payload_bytes: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in plain.chunks(payload_bytes) {
        out.extend_from_slice(&bgzf_block(chunk));
    }
    out.extend_from_slice(&bgzf_block(b""));
    out
}

/// `plain` as one plain gzip member, which carries no block layout.
#[cfg(feature = "gzip")]
pub(crate) fn gzip_bytes(plain: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    io::Write::write_all(&mut encoder, plain).unwrap();
    encoder.finish().unwrap()
}

/// A block-addressable view of `compressed`, through the seam an opened input uses, so no test
/// naming this names a concrete container.
#[cfg(feature = "gzip")]
pub(crate) fn block_access_of(compressed: Vec<u8>) -> Box<dyn crate::blocks::BlockAccess> {
    let Ok(access) = crate::blocks::block_access(compressed, crate::codecs::Codec::Bgzf) else {
        panic!("the fixture must be block-addressable");
    };
    access
}

/// One xz stream over `payload`, in blocks of `block_bytes` when one is given.
#[cfg(feature = "xz")]
pub(crate) fn xz_stream(
    payload: &[u8],
    block_bytes: Option<u64>,
    check: lzma_rust2::CheckType,
) -> Vec<u8> {
    use io::Write;
    let mut options = lzma_rust2::XzOptions::with_preset(0);
    options.set_check_sum_type(check);
    options.set_block_size(block_bytes.and_then(std::num::NonZeroU64::new));
    let mut writer = lzma_rust2::XzWriter::new(Vec::new(), options).unwrap();
    writer.write_all(payload).unwrap();
    writer.finish().unwrap()
}

/// `plain` as one single-block stream per chunk, which is what `cat a.xz b.xz` writes and the
/// cheapest multi-block fixture there is.
#[cfg(feature = "xz")]
pub(crate) fn xz_streams(plain: &[u8], chunk_bytes: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in plain.chunks(chunk_bytes) {
        out.extend_from_slice(&xz_stream(chunk, None, lzma_rust2::CheckType::Crc64));
    }
    out
}

/// `plain` as one lz4 frame written under `frame_info`.
#[cfg(feature = "lz4")]
pub(crate) fn lz4_frame(plain: &[u8], frame_info: lz4_flex::frame::FrameInfo) -> Vec<u8> {
    use io::Write;
    let mut encoder = lz4_flex::frame::FrameEncoder::with_frame_info(frame_info, Vec::new());
    encoder.write_all(plain).unwrap();
    encoder.finish().unwrap()
}

// endregion: Container Fixtures

#[cfg(test)]
mod tests {
    use super::*;

    /// An axis fixed at one value is an axis the suite does not test.
    /// A four-line proof that only truncates on long reads is how that ends.
    #[test]
    fn every_axis_takes_more_than_one_value() {
        assert!(
            ALL.iter().any(|shape| shape.read_bases > 4096),
            "no long read"
        );
        assert!(
            ALL.iter().any(|shape| shape.read_bases <= 64),
            "no short read"
        );
        assert!(
            ALL.iter().any(|shape| shape.wrap_columns > 0),
            "no wrapped record"
        );
        assert!(
            ALL.iter().any(|shape| shape.wrap_columns == 0),
            "no unwrapped record"
        );
        assert!(
            ALL.iter().any(|shape| shape.line_ending == "\r\n"),
            "no CRLF"
        );
        assert!(ALL.iter().any(|shape| shape.line_ending == "\n"), "no LF");
        assert!(
            ALL.iter().any(|shape| !shape.final_newline),
            "no missing final newline"
        );
        assert!(
            ALL.iter().any(|shape| shape.final_newline),
            "no final newline"
        );
        assert!(
            ALL.iter().any(|shape| shape.quality_opens_with_sigil),
            "no quality sigil"
        );
        assert!(
            ALL.iter().any(|shape| !shape.quality_opens_with_sigil),
            "no plain quality"
        );
        assert!(
            ALL.iter()
                .any(|shape| shape.format == SequenceFormat::Fasta),
            "no FASTA"
        );
        assert!(
            ALL.iter()
                .any(|shape| shape.format == SequenceFormat::Fastq),
            "no FASTQ"
        );
    }
}
