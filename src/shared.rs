//! Faster FASTA — the sequence-file layer every tool shares.
//!
//! Four regions:
//!
//! - __Records__ — what a record is, plus format detection and Phred conversions.
//!   A record is always a borrow; nothing here owns per-record bytes.
//! - __Parsing__ — one parser core returning [`ParseOutcome`], whose
//!   [`ParseOutcome::Incomplete`] separates "the buffer ended mid-record" from "the input is
//!   malformed". That single
//!   distinction lets the same parser serve a whole memory map, a partially filled stream
//!   window, and a byte range handed to one worker, which is why there is no second parser
//!   anywhere in this crate.
//! - __Containers__ — the three input shapes, and how records are read out of each.
//! - __Output__ — record formatting and the boilerplate every `main` shares.
//!
//! Scheduling lives in [`map_reduce`]. Nothing here spawns a thread or picks a shard.

pub mod map_reduce;

use std::io::{self, Read};

use stringzilla::sz::{find, find_byteset, Byteset};

// region: Records

/// Which of the two sequence file formats a stream carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceFormat {
    Fasta,
    Fastq,
}

impl SequenceFormat {
    /// The byte that opens a record header in this format.
    #[inline]
    pub fn sigil(self) -> u8 {
        match self {
            SequenceFormat::Fasta => b'>',
            SequenceFormat::Fastq => b'@',
        }
    }
}

/// One parsed record.
///
/// `sequence` and `quality` are logical values with every line break and interior space
/// already removed, so consumers never repeat that work and never have to ask whether the
/// underlying file was wrapped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Record<'a> {
    /// Header line including its leading `>` or `@`, excluding the line terminator.
    pub header: &'a [u8],
    /// Sequence bytes, whitespace removed.
    pub sequence: &'a [u8],
    /// Quality bytes, whitespace removed, equal in length to `sequence`.
    /// Empty for FASTA, which also makes an empty FASTQ read indistinguishable — consult
    /// the format rather than this field when the difference matters.
    pub quality: &'a [u8],
}

impl<'a> Record<'a> {
    /// Header with the format sigil stripped, which is what identifier matching wants.
    #[inline]
    pub fn header_without_sigil(&self) -> &'a [u8] {
        match self.header.split_first() {
            Some((&(b'>' | b'@'), rest)) => rest,
            _ => self.header,
        }
    }

    /// Identifier: the header up to the first space or tab, sigil stripped.
    pub fn identifier(&self) -> &'a [u8] {
        let without_sigil = self.header_without_sigil();
        match find_byteset(without_sigil, FIELD_SEPARATORS) {
            Some(position) => &without_sigil[..position],
            None => without_sigil,
        }
    }
}

/// Reusable normalization buffers owned by a parser and refilled per record.
///
/// Sequence and quality are kept apart so a single record can borrow both at once; one
/// shared buffer would make that impossible. Capacity is retained across records, so a run
/// over a uniform file allocates only while warming up.
#[derive(Debug, Default)]
pub struct RecordScratch {
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

impl RecordScratch {
    /// Empty buffers that grow to the widest record seen.
    pub fn new() -> Self {
        Self::default()
    }

    /// Preallocates both buffers for records of up to `bytes` each.
    pub fn with_capacity(bytes: usize) -> Self {
        Self {
            sequence: Vec::with_capacity(bytes),
            quality: Vec::with_capacity(bytes),
        }
    }

    /// Retained sequence capacity, for tests asserting that reuse actually happens.
    pub fn sequence_capacity(&self) -> usize {
        self.sequence.capacity()
    }

    /// Retained quality capacity, for tests asserting that reuse actually happens.
    pub fn quality_capacity(&self) -> usize {
        self.quality.capacity()
    }
}

/// Space and tab, the header field separators.
const FIELD_SEPARATORS: Byteset = Byteset::from_bytes(b" \t");

/// Everything stripped from a sequence or quality region.
const WHITESPACE: Byteset = Byteset::from_bytes(b" \t\r\n");

/// Detect the format from the first non-whitespace byte.
pub fn detect_format(data: &[u8]) -> io::Result<SequenceFormat> {
    for &byte in data.iter() {
        match byte {
            b' ' | b'\t' | b'\r' | b'\n' => continue,
            b'@' => return Ok(SequenceFormat::Fastq),
            b'>' => return Ok(SequenceFormat::Fasta),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("unknown format: expected '@' or '>', found {:?}", byte as char),
                ))
            }
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        "empty or whitespace-only input",
    ))
}

/// ASCII quality character to Phred score, Phred+33 as used by Illumina 1.8 and later.
#[inline]
pub fn ascii_to_phred33(ascii: u8) -> u8 {
    ascii.saturating_sub(33)
}

/// Phred score back to its ASCII quality character.
#[inline]
pub fn phred33_to_ascii(phred: u8) -> u8 {
    phred.saturating_add(33).min(126)
}

/// Mean Phred score across a quality string.
pub fn mean_quality(quality: &[u8]) -> f32 {
    if quality.is_empty() {
        return 0.0;
    }
    let total: u32 = quality.iter().map(|&byte| ascii_to_phred33(byte) as u32).sum();
    total as f32 / quality.len() as f32
}

// endregion: Records

// region: Parsing

/// Outcome of attempting to parse one record from the front of a buffer.
#[derive(Debug)]
pub enum ParseOutcome<'a> {
    /// A complete record, and how many bytes of the buffer it consumed.
    Record(Record<'a>, usize),
    /// The buffer ends mid-record. Read more and retry; this is not an error.
    Incomplete,
    /// The input is malformed and no amount of additional data will help.
    Invalid(io::Error),
}

fn invalid<'a>(message: &str) -> ParseOutcome<'a> {
    ParseOutcome::Invalid(io::Error::new(io::ErrorKind::InvalidData, message.to_string()))
}

/// At end of input a partial record is a truncation, not a request for more bytes.
fn incomplete_or_truncated<'a>(at_eof: bool, what: &str) -> ParseOutcome<'a> {
    if at_eof {
        invalid(&format!("truncated input: {what} is incomplete"))
    } else {
        ParseOutcome::Incomplete
    }
}

/// End of the line starting at `from`, exclusive of the terminator.
#[inline]
fn find_line_end(buffer: &[u8], from: usize) -> Option<usize> {
    if from >= buffer.len() {
        return None;
    }
    find(&buffer[from..], b"\n").map(|offset| from + offset)
}

/// First byte that is not a line terminator or blank.
#[inline]
fn skip_leading_whitespace(buffer: &[u8]) -> usize {
    let mut cursor = 0;
    while cursor < buffer.len() && matches!(buffer[cursor], b' ' | b'\t' | b'\r' | b'\n') {
        cursor += 1;
    }
    cursor
}

/// Line length excluding a trailing carriage return, so CRLF files count correctly.
#[inline]
fn length_without_carriage_return(line: &[u8]) -> usize {
    match line.last() {
        Some(b'\r') => line.len() - 1,
        _ => line.len(),
    }
}

#[inline]
fn trim_carriage_return(line: &[u8]) -> &[u8] {
    match line.last() {
        Some(b'\r') => &line[..line.len() - 1],
        _ => line,
    }
}

/// `region` with all whitespace removed.
///
/// Returns `region` itself when it holds none, which is the common case for an unwrapped
/// FASTQ read and costs one `find_byteset` scan and no copy. Otherwise the bytes are
/// gathered into `scratch` and a view of that is returned. Sharing one lifetime across the
/// input, the scratch, and the result is what lets a single return type cover both.
fn strip_whitespace<'a>(region: &'a [u8], scratch: &'a mut Vec<u8>) -> &'a [u8] {
    if find_byteset(region, WHITESPACE).is_none() {
        return region;
    }

    scratch.clear();
    scratch.reserve(region.len());

    let mut cursor = 0;
    while cursor < region.len() {
        let next = find_byteset(&region[cursor..], WHITESPACE)
            .map(|offset| cursor + offset)
            .unwrap_or(region.len());
        if next > cursor {
            // Each run of payload bytes is gathered in a single write.
            scratch.extend_from_slice(&region[cursor..next]);
        }
        cursor = if next < region.len() { next + 1 } else { next };
    }
    scratch
}

/// Parse one record of either format from the front of `buffer`.
pub fn parse_leading_record<'a>(
    buffer: &'a [u8],
    format: SequenceFormat,
    at_eof: bool,
    scratch: &'a mut RecordScratch,
) -> ParseOutcome<'a> {
    match format {
        SequenceFormat::Fasta => parse_leading_fasta_record(buffer, at_eof, scratch),
        SequenceFormat::Fastq => parse_leading_fastq_record(buffer, at_eof, scratch),
    }
}

/// Parse one FASTA record.
///
/// A FASTA record ends only where the next one begins, so a record is complete either when
/// a `>` is seen at a line start or when the input is known to be exhausted. That is why
/// `at_eof` is a parameter rather than something the parser could infer.
pub fn parse_leading_fasta_record<'a>(
    buffer: &'a [u8],
    at_eof: bool,
    scratch: &'a mut RecordScratch,
) -> ParseOutcome<'a> {
    let start = skip_leading_whitespace(buffer);
    if start >= buffer.len() {
        return ParseOutcome::Incomplete;
    }
    if buffer[start] != b'>' {
        return invalid("FASTA record must start with '>'");
    }

    let header_end = match find_line_end(buffer, start) {
        Some(position) => position,
        None if at_eof => buffer.len(),
        None => return ParseOutcome::Incomplete,
    };
    let header = trim_carriage_return(&buffer[start..header_end]);

    let sequence_start = (header_end + 1).min(buffer.len());
    let mut cursor = sequence_start;
    let mut sequence_end = sequence_start;
    let consumed;

    loop {
        if cursor >= buffer.len() {
            if at_eof {
                sequence_end = buffer.len();
                consumed = buffer.len();
                break;
            }
            return ParseOutcome::Incomplete;
        }
        if buffer[cursor] == b'>' {
            // The next record starts here, so this one is complete.
            consumed = cursor;
            break;
        }
        match find_line_end(buffer, cursor) {
            Some(position) => {
                sequence_end = position;
                cursor = position + 1;
            }
            None if at_eof => {
                sequence_end = buffer.len();
                consumed = buffer.len();
                break;
            }
            None => return ParseOutcome::Incomplete,
        }
    }

    let region = &buffer[sequence_start..sequence_end.max(sequence_start)];
    let sequence = strip_whitespace(region, &mut scratch.sequence);

    ParseOutcome::Record(
        Record {
            header,
            sequence,
            quality: b"",
        },
        consumed,
    )
}

/// Parse one FASTQ record.
///
/// Sequence lines run until a line opening with `+`; quality lines are then collected until
/// they match the sequence length. Counting against the sequence length is what makes
/// multi-line quality work, and is why quality is normalized through the same path as
/// sequence rather than sliced contiguously.
pub fn parse_leading_fastq_record<'a>(
    buffer: &'a [u8],
    at_eof: bool,
    scratch: &'a mut RecordScratch,
) -> ParseOutcome<'a> {
    let start = skip_leading_whitespace(buffer);
    if start >= buffer.len() {
        return ParseOutcome::Incomplete;
    }
    if buffer[start] != b'@' {
        return invalid("FASTQ record must start with '@'");
    }

    let header_end = match find_line_end(buffer, start) {
        Some(position) => position,
        None => return incomplete_or_truncated(at_eof, "FASTQ header"),
    };
    let header = trim_carriage_return(&buffer[start..header_end]);

    // Sequence lines, up to the '+' separator.
    let sequence_start = header_end + 1;
    let mut cursor = sequence_start;
    let mut sequence_end = sequence_start;

    loop {
        if cursor >= buffer.len() {
            return incomplete_or_truncated(at_eof, "FASTQ sequence");
        }
        if buffer[cursor] == b'+' {
            break;
        }
        match find_line_end(buffer, cursor) {
            Some(position) => {
                sequence_end = position;
                cursor = position + 1;
            }
            None => return incomplete_or_truncated(at_eof, "FASTQ sequence"),
        }
    }

    // Destructured so sequence and quality borrow separate buffers simultaneously.
    let RecordScratch {
        sequence: sequence_scratch,
        quality: quality_scratch,
    } = scratch;

    let sequence_region = &buffer[sequence_start..sequence_end.max(sequence_start)];
    let sequence = strip_whitespace(sequence_region, sequence_scratch);
    let sequence_length = sequence.len();

    let separator_end = match find_line_end(buffer, cursor) {
        Some(position) => position,
        None => return incomplete_or_truncated(at_eof, "FASTQ '+' separator"),
    };

    // Quality lines, matched against the sequence length.
    let quality_start = separator_end + 1;
    let mut cursor = quality_start;
    let mut quality_end;
    let mut collected = 0usize;

    loop {
        if cursor > buffer.len() {
            return incomplete_or_truncated(at_eof, "FASTQ quality");
        }
        match find_line_end(buffer, cursor) {
            Some(position) => {
                quality_end = position;
                collected += length_without_carriage_return(&buffer[cursor..position]);
                cursor = position + 1;
            }
            None if at_eof => {
                quality_end = buffer.len();
                collected += length_without_carriage_return(&buffer[cursor.min(buffer.len())..]);
                cursor = buffer.len();
            }
            None => return ParseOutcome::Incomplete,
        }
        if collected >= sequence_length {
            break;
        }
        if cursor >= buffer.len() {
            return incomplete_or_truncated(at_eof, "FASTQ quality");
        }
    }

    let quality_region = &buffer[quality_start..quality_end.max(quality_start)];
    let quality = strip_whitespace(quality_region, quality_scratch);

    if sequence.len() != quality.len() {
        return ParseOutcome::Invalid(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "FASTQ sequence and quality length mismatch: sequence {}, quality {}",
                sequence.len(),
                quality.len()
            ),
        ));
    }

    ParseOutcome::Record(
        Record {
            header,
            sequence,
            quality,
        },
        cursor.min(buffer.len()),
    )
}

/// First record boundary at or after `from`, for a worker handed an arbitrary byte range.
///
/// FASTA is unambiguous: a `>` opening a line always starts a record. FASTQ is not, because
/// `@` is a legal Phred+33 quality byte at Q31, so a candidate is accepted only after a full
/// four-line cycle checks out — a `+` on the third line and a fourth line matching the
/// second in length. Records wrapped across several lines cannot be resynchronized this way
/// and are only safe to read from the start of the file.
pub fn next_record_boundary(bytes: &[u8], from: usize, format: SequenceFormat) -> Option<usize> {
    let sigil = format.sigil();
    if from == 0 && bytes.first() == Some(&sigil) {
        return Some(0);
    }

    // A record opens a line, so the boundary is a newline immediately followed by the sigil.
    // Searching for both bytes at once leaves the whole scan to one SIMD pass on FASTA,
    // rather than one pass per line.
    let needle = [b'\n', sigil];
    let mut cursor = from;
    loop {
        let candidate = find(&bytes[cursor..], &needle)? + cursor + 1;
        if candidate >= bytes.len() {
            return None;
        }
        if format == SequenceFormat::Fasta || fastq_cycle_holds(bytes, candidate) {
            return Some(candidate);
        }
        cursor = candidate;
    }
}

/// Whether a candidate `@` really opens a FASTQ record, judged by the shape of the next
/// four lines rather than by the byte alone.
fn fastq_cycle_holds(bytes: &[u8], candidate: usize) -> bool {
    let mut cursor = candidate;
    let mut lines = [0usize; 4];
    for slot in lines.iter_mut() {
        let Some(end) = find_line_end(bytes, cursor) else { return false };
        *slot = length_without_carriage_return(&bytes[cursor..end]);
        cursor = end + 1;
    }
    let separator = candidate + lines[0] + 1 + lines[1] + 1;
    bytes.get(separator) == Some(&b'+') && lines[1] == lines[3]
}

// endregion: Parsing

// region: Containers

/// A whole input held in memory or memory-mapped, addressable at any offset.
///
/// Byte-addressable inputs can be split into as many shards as there are workers, since
/// [`next_record_boundary`] recovers a record boundary from any offset.
pub struct RandomAccess<'a> {
    data: &'a [u8],
}

/// A forward-only input: a pipe, a socket, or a decompressor over a plain gzip stream.
///
/// Nothing can be revisited and nothing can be split, so a stream is always exactly one
/// shard, which the type states rather than leaving to a runtime check.
pub struct ForwardAccess<R: Read> {
    reader: R,
    window: Vec<u8>,
}

/// Initial parse window for a [`ForwardAccess`], grown only when a single record does not fit.
pub const DEFAULT_WINDOW_BYTES: usize = 1 << 20;

impl<'a> RandomAccess<'a> {
    /// Wraps a buffer that outlives every record read from it.
    pub fn new(data: &'a [u8]) -> Self {
        Self { data }
    }

    /// The underlying bytes, for a caller that needs offsets rather than records.
    pub fn as_slice(&self) -> &'a [u8] {
        self.data
    }

    /// Format of the first record, from the leading non-whitespace byte.
    pub fn detect_format(&self) -> io::Result<SequenceFormat> {
        detect_format(self.data)
    }

    /// Visit every record in `range`, which must start at a record boundary.
    ///
    /// Records borrow either the input or a scratch buffer reused across the whole call, so
    /// this allocates only while warming up.
    pub fn for_each_record_in(
        &self,
        range: core::ops::Range<usize>,
        format: SequenceFormat,
        mut visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut scratch = RecordScratch::new();
        let region = &self.data[range];
        let mut consumed = 0usize;

        while consumed < region.len() {
            match parse_leading_record(&region[consumed..], format, true, &mut scratch) {
                ParseOutcome::Record(record, used) => {
                    visit(record)?;
                    debug_assert!(used > 0, "parser must make progress");
                    consumed += used;
                }
                // At end of input this only happens on trailing whitespace.
                ParseOutcome::Incomplete => break,
                ParseOutcome::Invalid(error) => return Err(error),
            }
        }
        Ok(())
    }

    /// Visit every record in the whole input.
    pub fn for_each_record(
        &self,
        format: SequenceFormat,
        visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        self.for_each_record_in(0..self.data.len(), format, visit)
    }
}

impl<R: Read> ForwardAccess<R> {
    /// Wraps a reader with the default parse window.
    pub fn new(reader: R) -> Self {
        Self::with_window(reader, DEFAULT_WINDOW_BYTES)
    }

    /// Wraps a reader with an explicit window, which the tests use to force refills.
    pub fn with_window(reader: R, window_bytes: usize) -> Self {
        Self {
            reader,
            // Allocated once and reused for the life of the run.
            window: vec![0u8; window_bytes.max(64)],
        }
    }

    /// Visit every record, in bounded memory.
    ///
    /// The window is compacted at exactly one point, where the only live index is `consumed`
    /// and it is reset in the same block, so no stale offset can survive a compaction.
    pub fn for_each_record(
        &mut self,
        format: SequenceFormat,
        mut visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut scratch = RecordScratch::new();
        let mut filled = 0usize;
        let mut consumed = 0usize;
        let mut at_eof = false;

        /// What the parse decided, carrying no borrow of the window.
        enum Outcome {
            Consumed(usize),
            NeedMore,
        }

        loop {
            let outcome = match parse_leading_record(&self.window[consumed..filled], format, at_eof, &mut scratch)
            {
                ParseOutcome::Record(record, used) => {
                    visit(record)?;
                    debug_assert!(used > 0, "parser must make progress");
                    Outcome::Consumed(used)
                }
                ParseOutcome::Incomplete => Outcome::NeedMore,
                ParseOutcome::Invalid(error) => return Err(error),
            };
            // The immutable borrow of `window` ends here, before any mutation below.

            match outcome {
                Outcome::Consumed(used) => {
                    consumed += used;
                    continue;
                }
                Outcome::NeedMore => {}
            }

            if at_eof {
                // At end of input the parser reports `Incomplete` only for trailing
                // whitespace; a genuine truncation comes back as `Invalid`.
                break;
            }

            if consumed > 0 {
                self.window.copy_within(consumed..filled, 0);
                filled -= consumed;
                consumed = 0;
            }

            if filled == self.window.len() {
                // One record is larger than the whole window, so grow rather than stall.
                self.window.resize(self.window.len() * 2, 0);
            }

            let read_bytes = self.reader.read(&mut self.window[filled..])?;
            if read_bytes == 0 {
                at_eof = true;
            } else {
                filled += read_bytes;
            }
        }
        Ok(())
    }
}

/// Bytes sniffed from a stream to decide its format before parsing starts.
const PEEK_BYTES: usize = 8 << 10;

impl ForwardAccess<Box<dyn Read>> {
    /// Read a prefix of standard input, detect the format, and replay it into the stream.
        pub fn from_stdin() -> io::Result<(SequenceFormat, Self)> {
        let mut prefix = Vec::with_capacity(PEEK_BYTES);
        io::stdin()
            .lock()
            .by_ref()
            .take(PEEK_BYTES as u64)
            .read_to_end(&mut prefix)?;
        let format = detect_format(&prefix)?;
        // `Chain` stops at the end of its first source, so this short-reads exactly once.
        // Harmless: the window treats a short read like any other.
        let reader: Box<dyn Read> = Box::new(io::Cursor::new(prefix).chain(io::stdin()));
        Ok((format, ForwardAccess::new(reader)))
    }
}

/// Memory-map a file for byte-addressable access.
pub fn map_file(path: &str) -> io::Result<memmap2::Mmap> {
    let file = std::fs::File::open(path)?;
    unsafe { memmap2::Mmap::map(&file) }
}

// endregion: Containers

// region: Output

/// Append one formatted record to `target`. Quality is emitted only for FASTQ output.
fn format_into(
    target: &mut Vec<u8>,
    format: SequenceFormat,
    header_body: &[u8],
    sequence: &[u8],
    quality: &[u8],
) {
    target.push(format.sigil());
    target.extend_from_slice(header_body);
    target.push(b'\n');
    target.extend_from_slice(sequence);
    target.push(b'\n');
    if format == SequenceFormat::Fastq {
        target.extend_from_slice(b"+\n");
        target.extend_from_slice(quality);
        target.push(b'\n');
    }
}

/// bytes. Staging into one reusable buffer makes it one call, and costs one allocation
/// for the whole run.
pub struct RecordWriter<W: io::Write> {
    writer: W,
    staging: Vec<u8>,
    format: SequenceFormat,
}

impl<W: io::Write> RecordWriter<W> {
    /// Wraps a writer, emitting records in `format` until [`RecordWriter::set_format`] says otherwise.
    pub fn new(writer: W, format: SequenceFormat) -> Self {
        Self {
            writer,
            staging: Vec::with_capacity(4 << 10),
            format,
        }
    }

    /// Output format, which need not match the input's — this is all of `fastq-to-fasta`.
    pub fn format(&self) -> SequenceFormat {
        self.format
    }

    /// Choose the output format after construction, which is how a mirroring tool
    /// applies the format handed to it before records start arriving.
    pub fn set_format(&mut self, format: SequenceFormat) {
        self.format = format;
    }

    /// Write `record` verbatim, re-sigiling the header to the output format.
    pub fn write_record(&mut self, record: Record<'_>) -> io::Result<()> {
        self.write_parts(record.header_without_sigil(), record.sequence, record.quality)
    }

    /// Write a record from parts, where `header_body` excludes the sigil.
    ///
    /// Quality is emitted only for FASTQ output, and a FASTA source converted to FASTQ
    /// would have none — callers wanting that must supply it.
    pub fn write_parts(
        &mut self,
        header_body: &[u8],
        sequence: &[u8],
        quality: &[u8],
    ) -> io::Result<()> {
        self.staging.clear();
        format_into(
            &mut self.staging,
            self.format,
            header_body,
            sequence,
            quality,
        );
        self.writer.write_all(&self.staging)
    }

    /// Format a record into `target`, replacing its contents and keeping its capacity.
    ///
    /// For a tool that must retain records: a record borrows the input, so surviving past
    /// the next one means owning the bytes.
    pub fn render_into(&self, record: Record<'_>, target: &mut Vec<u8>) {
        target.clear();
        format_into(
            target,
            self.format,
            record.header_without_sigil(),
            record.sequence,
            record.quality,
        );
    }

    /// Write bytes produced by [`RecordWriter::render`], which are already formatted.
    pub fn write_rendered(&mut self, rendered: &[u8]) -> io::Result<()> {
        self.writer.write_all(rendered)
    }

    /// Flushes the underlying writer; the staging buffer never holds a partial record.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Consumes the writer, which the tests use to recover a `Vec<u8>` sink.
    pub fn into_inner(self) -> W {
        self.writer
    }
}

/// Writer for `path`, or standard output when it is `None` or `"-"`.
pub fn open_output(path: Option<&str>) -> io::Result<Box<dyn io::Write>> {
    match path {
        None | Some("-") => Ok(Box::new(io::BufWriter::new(io::stdout()))),
        Some(path) => Ok(Box::new(io::BufWriter::new(std::fs::File::create(path)?))),
    }
}

/// The epilogue every tool's `main` shares: a broken pipe is a normal end, because
/// `tool | head` is a routine way to use these; anything else is fatal.
pub fn finish_or_exit(context: &str, result: io::Result<()>) {
    if let Err(error) = result {
        if error.kind() == io::ErrorKind::BrokenPipe {
            std::process::exit(0);
        }
        eprintln!("{context}: {error}");
        std::process::exit(1);
    }
}


// endregion: Output

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_both_formats() {
        assert_eq!(detect_format(b">seq1\nACGT\n").unwrap(), SequenceFormat::Fasta);
        assert_eq!(
            detect_format(b"@seq1\nACGT\n+\nIIII\n").unwrap(),
            SequenceFormat::Fastq
        );
    }

    #[test]
    fn detects_past_leading_whitespace() {
        assert_eq!(
            detect_format(b"  \n  >seq1\nACGT\n").unwrap(),
            SequenceFormat::Fasta
        );
    }

    #[test]
    fn rejects_unknown_and_empty() {
        assert_eq!(
            detect_format(b"ACGT\n").unwrap_err().kind(),
            io::ErrorKind::InvalidData
        );
        assert_eq!(
            detect_format(b"   \n").unwrap_err().kind(),
            io::ErrorKind::InvalidData
        );
    }

    #[test]
    fn strips_sigil_and_splits_identifier() {
        let record = Record {
            header: b">sp|P12345| some description here",
            sequence: b"ACGT",
            quality: b"",
        };
        assert_eq!(record.header_without_sigil(), b"sp|P12345| some description here");
        assert_eq!(record.identifier(), b"sp|P12345|");
    }

    #[test]
    fn identifier_without_description() {
        let record = Record {
            header: b"@read1",
            sequence: b"ACGT",
            quality: b"IIII",
        };
        assert_eq!(record.identifier(), b"read1");
    }

    #[test]
    fn quality_conversion_round_trips() {
        assert_eq!(ascii_to_phred33(b'!'), 0);
        assert_eq!(ascii_to_phred33(b'I'), 40);
        assert_eq!(phred33_to_ascii(0), b'!');
        assert_eq!(phred33_to_ascii(40), b'I');
    }

    #[test]
    fn mean_quality_values() {
        assert_eq!(mean_quality(b"IIII"), 40.0);
        assert_eq!(mean_quality(b"!!!!"), 0.0);
        assert_eq!(mean_quality(b"!I"), 20.0);
        assert_eq!(mean_quality(b""), 0.0);
    }

    /// Parses to exhaustion over a complete slice.
    fn parse_all(
        data: &[u8],
        format: SequenceFormat,
    ) -> io::Result<Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>> {
        let mut scratch = RecordScratch::new();
        let mut out = Vec::new();
        let mut cursor = 0;
        loop {
            let remaining = &data[cursor..];
            if skip_leading_whitespace(remaining) >= remaining.len() {
                return Ok(out);
            }
            match parse_leading_record(remaining, format, true, &mut scratch) {
                ParseOutcome::Record(record, consumed) => {
                    out.push((
                        record.header.to_vec(),
                        record.sequence.to_vec(),
                        record.quality.to_vec(),
                    ));
                    assert!(consumed > 0, "parser made no progress");
                    cursor += consumed;
                }
                ParseOutcome::Incomplete => return Ok(out),
                ParseOutcome::Invalid(error) => return Err(error),
            }
        }
    }

    #[test]
    fn fasta_single_line_records_borrow() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n";
        let mut scratch = RecordScratch::new();
        match parse_leading_fasta_record(data, true, &mut scratch) {
            ParseOutcome::Record(record, consumed) => {
                assert_eq!(record.header, b">seq1");
                assert_eq!(record.sequence, b"ACGT");
                assert_eq!(consumed, 11);
            }
            other => panic!("expected a record, got {other:?}"),
        }
        // No whitespace inside the region, so nothing was copied.
        assert_eq!(scratch.sequence_capacity(), 0);
    }

    #[test]
    fn fasta_multiline_is_joined() {
        let records = parse_all(b">seq1\nACGT\nTGCA\n>seq2\nAAAA\n", SequenceFormat::Fasta).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].1, b"ACGTTGCA");
        assert_eq!(records[1].1, b"AAAA");
    }

    #[test]
    fn fasta_tolerates_blank_lines_and_indentation() {
        let data = b">seq1\n  ACGT  \n\n\tTGCA\t\n>seq2\nAAAA \n TTTT\n";
        let records = parse_all(data, SequenceFormat::Fasta).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].1, b"ACGTTGCA");
        assert_eq!(records[1].1, b"AAAATTTT");
    }

    #[test]
    fn fasta_incomplete_without_eof() {
        let mut scratch = RecordScratch::new();
        // Without the next '>' or end of input, the record could still continue.
        assert!(matches!(
            parse_leading_fasta_record(b">seq1\nACGT\n", false, &mut scratch),
            ParseOutcome::Incomplete
        ));
        assert!(matches!(
            parse_leading_fasta_record(b">seq1\nACGT\n", true, &mut scratch),
            ParseOutcome::Record(..)
        ));
    }

    #[test]
    fn fastq_single_line_records() {
        let records = parse_all(
            b"@seq1\nACGT\n+\nIIII\n@seq2 with description\nTGCA\n+\nHHHH\n",
            SequenceFormat::Fastq,
        )
        .unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, b"@seq1");
        assert_eq!(records[0].1, b"ACGT");
        assert_eq!(records[0].2, b"IIII");
        assert_eq!(records[1].0, b"@seq2 with description");
        assert_eq!(records[1].2, b"HHHH");
    }

    /// Multi-line quality is joined to match its sequence, not sliced across the newlines
    /// that separate its lines.
    #[test]
    fn fastq_multiline_quality_is_joined() {
        let records =
            parse_all(b"@seq1\nACGT\nTGCA\n+\nIIII\nHHHH\n", SequenceFormat::Fastq).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].1, b"ACGTTGCA");
        assert_eq!(records[0].2, b"IIIIHHHH");
    }

    #[test]
    fn fastq_empty_sequence() {
        let records = parse_all(b"@seq1\n\n+\n\n", SequenceFormat::Fastq).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].1, b"");
        assert_eq!(records[0].2, b"");
    }

    #[test]
    fn fastq_length_mismatch_is_invalid() {
        let error = parse_all(b"@seq1\nACGT\n+\n!!!\n", SequenceFormat::Fastq).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn fastq_crlf_is_handled() {
        let records = parse_all(b"@seq1\r\nACGT\r\n+\r\nIIII\r\n", SequenceFormat::Fastq).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].0, b"@seq1");
        assert_eq!(records[0].1, b"ACGT");
        assert_eq!(records[0].2, b"IIII");
    }

    /// Every prefix of a valid input must be `Incomplete`, never `Invalid`. This is the
    /// property the whole driver depends on.
    #[test]
    fn every_prefix_is_incomplete_not_invalid() {
        let data = b"@seq1\nACGT\n+\nIIII\n";
        let mut scratch = RecordScratch::new();
        for length in 1..data.len() {
            match parse_leading_fastq_record(&data[..length], false, &mut scratch) {
                ParseOutcome::Incomplete => {}
                ParseOutcome::Record(..) => panic!("prefix of {length} bytes parsed as a full record"),
                ParseOutcome::Invalid(error) => {
                    panic!("prefix of {length} bytes reported as invalid: {error}")
                }
            }
        }
    }

    /// Scratch reuse: capacity may grow while warming up to the widest record, but a second
    /// pass over the same input must not allocate at all.
    #[test]
    fn scratch_allocates_only_while_warming_up() {
        let data = b">seq1\nACGT\nTGCA\n>seq2\nGGGGG\nCCCCC\n>seq3\nAAAA\nTTTT\n";
        let mut scratch = RecordScratch::new();

        let count_pass = |scratch: &mut RecordScratch| {
            let mut cursor = 0;
            let mut records = 0;
            while cursor < data.len() {
                match parse_leading_fasta_record(&data[cursor..], true, scratch) {
                    ParseOutcome::Record(_, consumed) => {
                        cursor += consumed;
                        records += 1;
                    }
                    other => panic!("expected a record, got {other:?}"),
                }
            }
            records
        };

        assert_eq!(count_pass(&mut scratch), 3);
        let warmed = scratch.sequence_capacity();
        assert!(warmed > 0, "multi-line records should use scratch");

        assert_eq!(count_pass(&mut scratch), 3);
        assert_eq!(
            scratch.sequence_capacity(),
            warmed,
            "second pass reallocated scratch"
        );
    }

    /// Re-emits every record verbatim, so a driver run can be compared byte for byte
    /// against another driver run over the same input.
    #[derive(Debug, Default)]
    struct Collector {
        output: Vec<u8>,
        records: usize,
        finished: bool,
    }

    impl Collector {
        fn push(&mut self, record: Record<'_>) -> io::Result<()> {
            self.output.extend_from_slice(record.header);
            self.output.push(b'\n');
            self.output.extend_from_slice(record.sequence);
            self.output.push(b'\n');
            if !record.quality.is_empty() {
                self.output.extend_from_slice(b"+\n");
                self.output.extend_from_slice(record.quality);
                self.output.push(b'\n');
            }
            self.records += 1;
            Ok(())
        }

        fn finish(&mut self) -> io::Result<()> {
            self.finished = true;
            Ok(())
        }
    }

    /// Returns at most `limit` bytes per call, the way a pipe or socket does. `Cursor`
    /// always satisfies the full request and so can never drive the window through a
    /// mid-record refill.
    struct ShortReader<R: Read> {
        inner: R,
        limit: usize,
    }

    impl<R: Read> Read for ShortReader<R> {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            let capped = out.len().min(self.limit);
            self.inner.read(&mut out[..capped])
        }
    }

    fn sample_fastq(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            let width = 8 + (index % 7) * 3;
            let base = [b'A', b'C', b'G', b'T'][index % 4];
            data.extend_from_slice(format!("@read{index} sample record\n").as_bytes());
            data.extend(std::iter::repeat(base).take(width));
            data.extend_from_slice(b"\n+\n");
            data.extend(std::iter::repeat(b'I').take(width));
            data.push(b'\n');
        }
        data
    }

    fn sample_fasta(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            let width = 9 + (index % 5) * 4;
            let base = [b'A', b'C', b'G', b'T'][index % 4];
            data.extend_from_slice(format!(">contig{index} sample record\n").as_bytes());
            // Wrapped at 6 columns, so most records exercise the scratch path.
            for chunk_start in (0..width).step_by(6) {
                let chunk = (width - chunk_start).min(6);
                data.extend(std::iter::repeat(base).take(chunk));
                data.push(b'\n');
            }
        }
        data
    }

    fn drive_bytes(data: &[u8], format: SequenceFormat) -> io::Result<Collector> {
        let mut sink = Collector::default();
        RandomAccess::new(data).for_each_record(format, |record| sink.push(record))?;
        sink.finish()?;
        Ok(sink)
    }

    fn drive_stream(
        data: &[u8],
        format: SequenceFormat,
        limit: usize,
        window: usize,
    ) -> io::Result<Collector> {
        let mut sink = Collector::default();
        let reader = ShortReader {
            inner: io::Cursor::new(data),
            limit,
        };
        ForwardAccess::with_window(reader, window)
            .for_each_record(format, |record| sink.push(record))?;
        sink.finish()?;
        Ok(sink)
    }

    /// The property the whole design rests on: however the bytes arrive, the records come
    /// out identical.
    #[test]
    fn stream_matches_slice_for_fastq() {
        let data = sample_fastq(40);
        let expected = drive_bytes(&data, SequenceFormat::Fastq).unwrap();
        assert_eq!(expected.records, 40);

        for limit in [1, 3, 7, 13, 64, 4096] {
            for window in [64, 128, 1024, DEFAULT_WINDOW_BYTES] {
                let actual = drive_stream(&data, SequenceFormat::Fastq, limit, window)
                    .unwrap_or_else(|error| panic!("limit {limit}, window {window}: {error}"));
                assert_eq!(
                    actual.output, expected.output,
                    "limit {limit}, window {window}: output diverged"
                );
                assert_eq!(actual.records, expected.records);
            }
        }
    }

    #[test]
    fn stream_matches_slice_for_fasta() {
        let data = sample_fasta(30);
        let expected = drive_bytes(&data, SequenceFormat::Fasta).unwrap();
        assert_eq!(expected.records, 30);

        for limit in [1, 5, 11, 64, 4096] {
            for window in [64, 256, DEFAULT_WINDOW_BYTES] {
                let actual = drive_stream(&data, SequenceFormat::Fasta, limit, window)
                    .unwrap_or_else(|error| panic!("limit {limit}, window {window}: {error}"));
                assert_eq!(
                    actual.output, expected.output,
                    "limit {limit}, window {window}: output diverged"
                );
            }
        }
    }

    /// A record larger than the initial window must grow it, not stall or truncate.
    #[test]
    fn record_larger_than_window_grows_it() {
        let mut data = Vec::from(&b"@long\n"[..]);
        data.extend(std::iter::repeat(b'A').take(5000));
        data.extend_from_slice(b"\n+\n");
        data.extend(std::iter::repeat(b'I').take(5000));
        data.push(b'\n');

        let sink = drive_stream(&data, SequenceFormat::Fastq, 17, 64).unwrap();
        assert_eq!(sink.records, 1);
        assert!(sink.output.windows(5000).any(|w| w.iter().all(|&b| b == b'A')));
    }

    #[test]
    fn truncated_input_is_an_error() {
        // Quality is one byte short of the sequence and the stream ends.
        let data = b"@seq1\nACGT\n+\nIII";
        let result = drive_stream(data, SequenceFormat::Fastq, 4, 64);
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn empty_input_yields_no_records() {
        let sink = drive_stream(b"", SequenceFormat::Fastq, 8, 64).unwrap();
        assert_eq!(sink.records, 0);
        assert!(sink.finished);
    }

    #[test]
    fn trailing_whitespace_is_not_an_error() {
        let data = b"@seq1\nACGT\n+\nIIII\n\n\n  \n";
        let sink = drive_stream(data, SequenceFormat::Fastq, 3, 64).unwrap();
        assert_eq!(sink.records, 1);
    }

    #[test]
    fn finish_runs_even_with_no_records() {
        let mut sink = Collector::default();
        RandomAccess::new(b"")
            .for_each_record(SequenceFormat::Fasta, |record| sink.push(record))
            .unwrap();
        sink.finish().unwrap();
        assert!(sink.finished);
    }

    #[test]
    fn format_detection_feeds_the_driver() {
        let data = sample_fasta(3);
        let format = detect_format(&data).unwrap();
        assert_eq!(format, SequenceFormat::Fasta);
        assert_eq!(drive_bytes(&data, format).unwrap().records, 3);
    }
}
