//! Faster FASTA — the inputs an operating system hands a tool, and the bytes it hands back.
//!
//! Two regions:
//!
//! - __Containers__ — the three input shapes, and how records are read out of each.
//! - __Output__ — record formatting and the boilerplate every `main` shares.
//!
//! This layer sits above [`crate::records`] and reads every record through that parser.
//! The dependency runs one way only, so nothing in [`crate::records`] names anything here and
//! a record can be parsed with no file in sight.
//!
//! Scheduling lives in [`crate::scheduling`].
//! Nothing here spawns a thread or picks a range, and nothing here declares a sibling module.

use std::io::{self, Read};
use std::ops::Deref;

use crate::records::{
    next_record_start, parse_leading_record, sequence_format, ParseOutcome, Record, RecordScratch,
    RecordStart, SequenceFormat,
};

// region: Containers

/// A whole input held in memory or memory-mapped, addressable at any offset.
///
/// Any byte range can go to a worker, since [`crate::records::next_record_boundary`]
/// recovers a record
/// boundary from any offset. Generic over its container so one type serves the owner of a
/// mapping and the borrower of a view of it alike.
pub struct RandomAccess<B: Deref<Target = [u8]>> {
    data: B,
}

/// A forward-only input: a pipe, a socket, or a decompressor over a plain gzip stream.
///
/// Nothing can be revisited and nothing can be split, so a stream is read once from front to
/// back, which the type states rather than leaving to a runtime check.
pub struct ForwardAccess<R: Read> {
    reader: R,
    window: Vec<u8>,
}

/// Where a run over a stream begins, and what becomes of the bytes before that point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamStart {
    /// Byte zero, which is where a whole file or a whole pipe begins.
    FirstByte,
    /// The first record boundary of the given format, dropping the partial record before it.
    ///
    /// A run of blocks starts at a block, and a block boundary is not a record boundary, so
    /// those leading bytes are the tail of a record the run before owns and emits. The format
    /// comes from the whole input rather than from these bytes, because a window opening
    /// mid-record has no sigil to detect and finding the boundary needs the format anyway.
    FirstBoundary(SequenceFormat),
}

/// Initial parse window for a [`ForwardAccess`], grown only when a single record does not fit.
pub const DEFAULT_WINDOW_BYTES: usize = 1 << 20;

impl<B: Deref<Target = [u8]>> RandomAccess<B> {
    /// Wraps a container that outlives every record read from it.
    pub fn new(data: B) -> Self {
        Self { data }
    }

    /// A two-word view of the same bytes, for a worker that must not own the container.
    pub fn borrowed(&self) -> RandomAccess<&[u8]> {
        RandomAccess { data: &self.data }
    }

    /// The bytes themselves, for a scheduler that slices before anything is parsed.
    pub fn as_slice(&self) -> &[u8] {
        &self.data
    }

    /// Visit every record in `range`, which must start at a record boundary.
    ///
    /// Records borrow either the input or a scratch buffer reused across the whole call, so
    /// this allocates only while warming up.
    pub fn for_each_record_in(
        &self,
        range: core::ops::Range<usize>,
        mut visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut scratch = RecordScratch::new();
        let region = &self.data[range];
        // A region with no records is empty, not malformed, so there is nothing to visit.
        let Some(format) = sequence_format(region)? else {
            return Ok(());
        };
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
        visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        self.for_each_record_in(0..self.data.len(), visit)
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
    pub fn for_each_record(
        &mut self,
        visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        self.for_each_record_within(StreamStart::FirstByte, None, visit)
    }

    /// The source itself, for a driver that reads whole slabs rather than records.
    ///
    /// The parse window keeps whatever it holds, so a caller reads records or slabs from one
    /// stream, never both.
    pub fn source_mut(&mut self) -> &mut R {
        &mut self.reader
    }

    /// Visit the records beginning in this stretch of the stream, in bounded memory.
    ///
    /// `decoded_bytes` is how much of the stream belongs to this unit of work, or `None` to
    /// read to the end. It is consulted before each record rather than after: a record
    /// straddling the end began inside and is emitted whole, while the first record starting
    /// past it is left alone. Consulting afterwards would emit that one twice.
    ///
    /// The window is compacted at exactly one point, where the only live index is `consumed`
    /// and it is reset in the same block, so no stale offset can survive a compaction.
    pub fn for_each_record_within(
        &mut self,
        start: StreamStart,
        decoded_bytes: Option<usize>,
        mut visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut scratch = RecordScratch::new();
        let mut filled = 0usize;
        let mut consumed = 0usize;
        let mut at_eof = false;

        let mut parsed_total = 0usize;

        // A whole stream announces its format with its first record; a window opening
        // mid-record cannot, so it is told. Either way a stream with no records is empty
        // rather than malformed.
        let format = match start {
            StreamStart::FirstBoundary(format) => format,
            StreamStart::FirstByte => loop {
                match sequence_format(&self.window[..filled])? {
                    Some(format) => break format,
                    None if at_eof => return Ok(()),
                    None => {
                        if filled == self.window.len() {
                            self.window.resize(self.window.len() * 2, 0);
                        }
                        let read_bytes = self.reader.read(&mut self.window[filled..])?;
                        if read_bytes == 0 {
                            at_eof = true;
                        } else {
                            filled += read_bytes;
                        }
                    }
                }
            },
        };

        /// What the parse decided, carrying no borrow of the window.
        enum Outcome {
            Consumed(usize),
            NeedMore,
        }

        // A window opening mid-record inherits the tail of a record the unit of work before
        // it emits. Skipping it needs the boundary search rather than the parser, because the
        // bytes ahead of it are not a record.
        if matches!(start, StreamStart::FirstBoundary(_)) {
            let mut searched_through = 0usize;
            loop {
                match next_record_start(&self.window[..filled], searched_through, format) {
                    RecordStart::At(boundary) => {
                        consumed = boundary;
                        // The skipped prefix is decoded input too, and `decoded_bytes` counts
                        // from the front of this window, so it must be counted here.
                        parsed_total = boundary;
                        break;
                    }
                    // Nothing more will arrive, so an unproven candidate stays unproven.
                    RecordStart::Undecided(_) | RecordStart::NotFound if at_eof => return Ok(()),
                    outcome => {
                        // Resume where the proof ran out of bytes, not where the scan stopped.
                        // A FASTQ cycle spans four lines, which on a long read is tens of
                        // kilobytes, so a candidate near the end is undecided rather than
                        // rejected, and stepping over it drops every record up to the next.
                        searched_through = match outcome {
                            RecordStart::Undecided(resume_from) => resume_from,
                            // Every candidate was judged and rejected; only the newline a
                            // two-byte needle straddles has to be reconsidered.
                            _ => filled.saturating_sub(1),
                        };
                        if filled == self.window.len() {
                            self.window.resize(self.window.len() * 2, 0);
                        }
                        let read_bytes = self.reader.read(&mut self.window[filled..])?;
                        if read_bytes == 0 {
                            at_eof = true;
                        } else {
                            filled += read_bytes;
                        }
                    }
                }
            }
        }

        loop {
            if decoded_bytes.is_some_and(|budget| parsed_total >= budget) {
                return Ok(());
            }
            let outcome = match parse_leading_record(
                &self.window[consumed..filled],
                format,
                at_eof,
                &mut scratch,
            ) {
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
                    parsed_total += used;
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

/// One opened input, in whichever shape its bytes allow. Variant names are the tier names.
///
/// The shape follows from the pair of codec and addressability rather than from the path: a
/// plain regular file can be mapped and addressed at any offset, while anything forward-only
/// — standard input, a named pipe, a file being decompressed a window at a time — is read
/// once from front to back.
///
/// Each variant owns its container rather than borrowing one, so the enum carries no
/// lifetime and a worker takes a view through `borrowed` while the mapping stays here.
pub enum InputFile {
    /// BYTE tier: every offset is addressable, so any byte range can go to a worker.
    Bytes(RandomAccess<memmap2::Mmap>),
    /// BLOCK tier: addressable at block granularity, since every block decodes alone.
    ///
    /// Boxed because three containers qualify and they share only their trait. Dispatch is one
    /// virtual call per block against decoding tens of kilobytes, so it does not register.
    #[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
    Blocks(Box<dyn crate::blocks::BlockAccess>),
    /// STREAM tier: forward-only, so decoding stays serial.
    Stream(ForwardAccess<Box<dyn Read>>),
}

impl InputFile {
    /// Open `path`, or standard input when it is `None` or `"-"`, deciding from the bytes
    /// themselves whether they can be addressed or only streamed.
    ///
    /// Mirrors [`open_output`] and [`Rendering::for_output`], with one asymmetry worth
    /// naming: an output's shape follows from its destination, because the bytes are ours to
    /// choose; an input's shape follows from its bytes, because they are not.
    pub fn open(path: Option<&str>) -> io::Result<Self> {
        let (source, name): (Box<dyn Read>, &str) = match path {
            None | Some("-") => (Box::new(io::stdin()), "standard input"),
            Some(path) => {
                let file = std::fs::File::open(path)?;
                // One `fstat` in place of a mapping that fails on a fifo with an error
                // nobody can act on. A process substitution names a pipe, and reads as a
                // stream rather than failing.
                if file.metadata()?.file_type().is_file() {
                    return Self::from_mapping(unsafe { memmap2::Mmap::map(&file)? }, path);
                }
                (Box::new(file), path)
            }
        };
        // Only a forward-only source needs its prefix held aside and replayed.
        let sniffed = crate::codecs::PeekedReader::sniffing(source)?;
        let codec = crate::codecs::Codec::detect(sniffed.peeked());
        Ok(InputFile::Stream(ForwardAccess::new(
            codec
                .decode(sniffed)
                .map_err(|error| at_source(name, &error))?,
        )))
    }

    /// Decide a mapped regular file's tier from a prefix of the mapping itself.
    ///
    /// A mapping is random access, so the codec comes straight out of it with no prefix held
    /// aside, and a compressed one that is not block-addressable decodes from that same
    /// mapping through an [`io::Cursor`] rather than from a second pass over the file.
    fn from_mapping(mapping: memmap2::Mmap, name: &str) -> io::Result<Self> {
        let prefix_bytes = mapping.len().min(crate::codecs::Codec::PEEK_BYTES);
        let codec = crate::codecs::Codec::detect(&mapping[..prefix_bytes]);
        let mapping = match codec {
            crate::codecs::Codec::Plain => return Ok(InputFile::Bytes(RandomAccess::new(mapping))),
            // Whether a container is addressable is a second question from which container it
            // is, and only a whole mapping can answer it: a single-block xz file and a blocked
            // one share every byte of their magic. A container that declines hands the mapping
            // back and reads as a stream.
            #[cfg(any(feature = "gzip", feature = "xz", feature = "lz4"))]
            _ => match crate::blocks::block_access(mapping, codec) {
                Ok(blocks) => return Ok(InputFile::Blocks(blocks)),
                Err(mapping) => mapping,
            },
            #[cfg(not(any(feature = "gzip", feature = "xz", feature = "lz4")))]
            _ => mapping,
        };
        Ok(InputFile::Stream(ForwardAccess::new(
            codec
                .decode(io::Cursor::new(mapping))
                .map_err(|error| at_source(name, &error))?,
        )))
    }

    /// Visit every record, whichever shape the input has.
    ///
    /// The block tier reads forward from block zero, so a serial run never asks for a block
    /// table and never walks the chain to build one.
    pub fn for_each_record(
        &mut self,
        visit: impl FnMut(Record<'_>) -> io::Result<()>,
    ) -> io::Result<()> {
        match self {
            InputFile::Bytes(bytes) => bytes.for_each_record(visit),
            #[cfg(feature = "gzip")]
            InputFile::Blocks(blocks) => {
                ForwardAccess::new(blocks.reader(0)).for_each_record(visit)
            }
            InputFile::Stream(stream) => stream.for_each_record(visit),
        }
    }
}

/// Prefixes an error with the input it came from, so a run over many paths says which.
///
/// Borrows rather than consumes: the message is rebuilt from the kind and the text, and the
/// original is dropped by the caller either way.
fn at_source(name: &str, error: &io::Error) -> io::Error {
    io::Error::new(error.kind(), format!("{name}: {error}"))
}

// endregion: Containers

// region: Output

/// How records are rendered, decided by whether the destination is a terminal.
///
/// A human reading a chromosome wants it wrapped and coloured; a file or pipe wants neither,
/// so redirecting and piping produce identical bytes and nothing downstream meets an escape
/// sequence. This is the convention `ls` follows, and the reason the default is detected
/// rather than flagged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rendering {
    /// Columns per sequence line; `0` writes each sequence on one line.
    pub line_width: usize,
    /// Colour nucleotides by base.
    pub color: bool,
}

impl Rendering {
    /// Plain: one line per sequence, no colour. What a file or pipe receives.
    pub const PLAIN: Self = Self {
        line_width: 0,
        color: false,
    };

    /// Wrapped at 60 columns and coloured. What a terminal receives.
    pub const TERMINAL: Self = Self {
        line_width: 60,
        color: true,
    };

    /// [`Rendering::TERMINAL`] when `path` names standard output and that is a terminal,
    /// [`Rendering::PLAIN`] otherwise.
    pub fn for_output(path: Option<&str>) -> Self {
        let to_stdout = matches!(path, None | Some("-"));
        if to_stdout && std::io::IsTerminal::is_terminal(&io::stdout()) {
            Self::TERMINAL
        } else {
            Self::PLAIN
        }
    }
}

/// Foreground colour for one nucleotide, or `None` to leave it uncoloured.
fn base_color(base: u8) -> Option<&'static [u8]> {
    match base.to_ascii_uppercase() {
        b'A' => Some(b"\x1b[32m"),
        b'C' => Some(b"\x1b[34m"),
        b'G' => Some(b"\x1b[33m"),
        b'T' | b'U' => Some(b"\x1b[31m"),
        b'N' => Some(b"\x1b[90m"),
        _ => None,
    }
}

/// Append `sequence`, wrapped per `rendering`, followed by a line terminator.
///
/// Split by case so the plain form — what every file and pipe receives — is one copy per
/// line rather than a decision per byte.
fn sequence_into(target: &mut Vec<u8>, sequence: &[u8], rendering: Rendering) {
    match (rendering.line_width, rendering.color) {
        (0, false) => target.extend_from_slice(sequence),
        (width, false) => {
            for line in sequence.chunks(width) {
                target.extend_from_slice(line);
                target.push(b'\n');
            }
            return;
        }
        (width, true) => {
            colored_sequence_into(target, sequence, width);
            return;
        }
    }
    target.push(b'\n');
}

/// Append `sequence` wrapped at `width` and coloured by base.
///
/// Reached only when the destination is a terminal, so its cost never falls on a pipe.
/// Colour is emitted per run rather than per byte, so a homopolymer costs one escape.
fn colored_sequence_into(target: &mut Vec<u8>, sequence: &[u8], width: usize) {
    const RESET: &[u8] = b"\x1b[0m";
    let width = if width == 0 { usize::MAX } else { width };

    for line in sequence.chunks(width) {
        let mut active: Option<&'static [u8]> = None;
        for &base in line {
            let wanted = base_color(base);
            if wanted != active {
                if active.is_some() {
                    target.extend_from_slice(RESET);
                }
                if let Some(escape) = wanted {
                    target.extend_from_slice(escape);
                }
                active = wanted;
            }
            target.push(base);
        }
        if active.is_some() {
            target.extend_from_slice(RESET);
        }
        target.push(b'\n');
    }
}

/// Append one formatted record to `target`. Quality is emitted only for FASTQ output.
fn format_into(
    target: &mut Vec<u8>,
    format: SequenceFormat,
    rendering: Rendering,
    header_body: &[u8],
    sequence: &[u8],
    quality: &[u8],
) {
    if rendering.color {
        target.extend_from_slice(b"\x1b[1;36m");
    }
    target.push(format.sigil());
    target.extend_from_slice(header_body);
    if rendering.color {
        target.extend_from_slice(b"\x1b[0m");
    }
    target.push(b'\n');

    sequence_into(target, sequence, rendering);

    if format == SequenceFormat::Fastq {
        target.extend_from_slice(b"+\n");
        // Quality wraps with the sequence so the two stay line-aligned, and is never
        // coloured, since its bytes are scores rather than bases.
        sequence_into(
            target,
            quality,
            Rendering {
                color: false,
                ..rendering
            },
        );
    }
}

/// Renders records to a sink, in whichever format the output carries.
///
/// A record written field by field would be four or five small writes, and a pipe charges for
/// bytes. Staging into one reusable buffer makes it one call, and costs one allocation
/// for the whole run.
pub struct RecordWriter<W: io::Write> {
    writer: W,
    staging: Vec<u8>,
    format: SequenceFormat,
    /// Format settled by [`RecordWriter::adopt`], so a later disagreement is caught.
    adopted: Option<SequenceFormat>,
    rendering: Rendering,
}

impl<W: io::Write> RecordWriter<W> {
    /// Wraps a writer, emitting plain records in `format` until told otherwise.
    pub fn new(writer: W, format: SequenceFormat) -> Self {
        Self {
            writer,
            staging: Vec::with_capacity(4 << 10),
            format,
            adopted: None,
            rendering: Rendering::PLAIN,
        }
    }

    /// Wraps a writer that renders per `rendering`.
    pub fn with_rendering(writer: W, rendering: Rendering, format: SequenceFormat) -> Self {
        Self {
            rendering,
            ..Self::new(writer, format)
        }
    }

    /// Mirror `format` onto the output, settling on the first record's and refusing a later
    /// change, since one stream cannot carry both.
    ///
    /// Idempotent, so a tool can call it per record instead of needing a hook that fires
    /// before the records do.
    pub fn adopt(&mut self, format: SequenceFormat) -> io::Result<()> {
        match self.adopted {
            Some(first) if first != format => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("input is {format:?} but an earlier record was {first:?}; one output stream cannot carry both"),
            )),
            Some(_) => Ok(()),
            None => {
                self.adopted = Some(format);
                self.format = format;
                Ok(())
            }
        }
    }

    /// Write `record` verbatim, re-sigiling the header to the output format.
    pub fn write_record(&mut self, record: Record<'_>) -> io::Result<()> {
        self.write_parts(
            record.header_without_sigil(),
            record.sequence,
            record.quality,
        )
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
            self.rendering,
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
            self.rendering,
            record.header_without_sigil(),
            record.sequence,
            record.quality,
        );
    }

    /// Write bytes produced by [`RecordWriter::render_into`], which are already formatted.
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

    /// The sink itself, for a caller draining a buffer between units of work.
    ///
    /// A record never straddles two calls, so whatever the sink holds is a whole number of
    /// records and can be forwarded as it stands.
    pub fn inner_mut(&mut self) -> &mut W {
        &mut self.writer
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
            data.extend(std::iter::repeat_n(base, width));
            data.extend_from_slice(b"\n+\n");
            data.extend(std::iter::repeat_n(b'I', width));
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
                data.extend(std::iter::repeat_n(base, chunk));
                data.push(b'\n');
            }
        }
        data
    }

    fn drive_bytes(data: &[u8]) -> io::Result<Collector> {
        let mut sink = Collector::default();
        RandomAccess::new(data).for_each_record(|record| sink.push(record))?;
        sink.finish()?;
        Ok(sink)
    }

    fn drive_stream(data: &[u8], limit: usize, window: usize) -> io::Result<Collector> {
        let mut sink = Collector::default();
        let reader = ShortReader {
            inner: io::Cursor::new(data),
            limit,
        };
        ForwardAccess::with_window(reader, window).for_each_record(|record| sink.push(record))?;
        sink.finish()?;
        Ok(sink)
    }

    /// The property the whole design rests on: however the bytes arrive, the records come
    /// out identical.
    #[test]
    fn stream_matches_slice_for_fastq() {
        let data = sample_fastq(40);
        let expected = drive_bytes(&data).unwrap();
        assert_eq!(expected.records, 40);

        for limit in [1, 3, 7, 13, 64, 4096] {
            for window in [64, 128, 1024, DEFAULT_WINDOW_BYTES] {
                let actual = drive_stream(&data, limit, window)
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
        let expected = drive_bytes(&data).unwrap();
        assert_eq!(expected.records, 30);

        for limit in [1, 5, 11, 64, 4096] {
            for window in [64, 256, DEFAULT_WINDOW_BYTES] {
                let actual = drive_stream(&data, limit, window)
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
        data.extend(std::iter::repeat_n(b'A', 5000));
        data.extend_from_slice(b"\n+\n");
        data.extend(std::iter::repeat_n(b'I', 5000));
        data.push(b'\n');

        let sink = drive_stream(&data, 17, 64).unwrap();
        assert_eq!(sink.records, 1);
        assert!(sink
            .output
            .windows(5000)
            .any(|w| w.iter().all(|&b| b == b'A')));
    }

    #[test]
    fn truncated_input_is_an_error() {
        // Quality is one byte short of the sequence and the stream ends.
        let data = b"@seq1\nACGT\n+\nIII";
        let result = drive_stream(data, 4, 64);
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn empty_input_yields_no_records() {
        let sink = drive_stream(b"", 8, 64).unwrap();
        assert_eq!(sink.records, 0);
        assert!(sink.finished);
    }

    #[test]
    fn trailing_whitespace_is_not_an_error() {
        let data = b"@seq1\nACGT\n+\nIIII\n\n\n  \n";
        let sink = drive_stream(data, 3, 64).unwrap();
        assert_eq!(sink.records, 1);
    }

    #[test]
    fn finish_runs_even_with_no_records() {
        let mut sink = Collector::default();
        RandomAccess::new(&b""[..])
            .for_each_record(|record| sink.push(record))
            .unwrap();
        sink.finish().unwrap();
        assert!(sink.finished);
    }

    #[test]
    fn format_detection_feeds_the_driver() {
        let data = sample_fasta(3);
        let format = sequence_format(&data).unwrap().unwrap();
        assert_eq!(format, SequenceFormat::Fasta);
        assert_eq!(drive_bytes(&data).unwrap().records, 3);
    }

    /// A worker reads through a view while the container stays with its owner.
    #[test]
    fn a_borrowed_view_yields_the_same_records() {
        let data = sample_fasta(12);
        let owned = RandomAccess::new(data.clone());
        let mut through_view = Vec::new();
        owned
            .borrowed()
            .for_each_record(|record| {
                through_view.push(record.header.to_vec());
                Ok(())
            })
            .unwrap();
        assert_eq!(through_view.len(), 12);
        assert_eq!(through_view, headers_in(&data));
    }

    /// Headers of every record in a buffer, as the tier-independent answer to compare against.
    fn headers_in(data: &[u8]) -> Vec<Vec<u8>> {
        let mut headers = Vec::new();
        RandomAccess::new(data)
            .for_each_record(|record| {
                headers.push(record.header.to_vec());
                Ok(())
            })
            .unwrap();
        headers
    }

    /// Writes `bytes` under `directory` and names the file.
    fn written(directory: &tempfile::TempDir, name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = directory.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    /// Header and sequence of every record [`InputFile::open`] yields for `path`.
    fn records_at(path: &std::path::Path) -> io::Result<Vec<(Vec<u8>, Vec<u8>)>> {
        let mut records = Vec::new();
        InputFile::open(Some(path.to_str().unwrap()))?.for_each_record(|record| {
            records.push((record.header.to_vec(), record.sequence.to_vec()));
            Ok(())
        })?;
        Ok(records)
    }

    #[test]
    fn a_plain_file_opens_as_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let path = written(&directory, "plain.fasta", &sample_fasta(8));
        assert!(matches!(
            InputFile::open(Some(path.to_str().unwrap())).unwrap(),
            InputFile::Bytes(_)
        ));
        assert_eq!(records_at(&path).unwrap().len(), 8);
    }

    #[test]
    fn a_missing_path_is_an_error() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("absent.fasta");
        match InputFile::open(Some(path.to_str().unwrap())) {
            Ok(_) => panic!("a path that names nothing must not open"),
            Err(error) => assert_eq!(error.kind(), io::ErrorKind::NotFound),
        }
    }

    /// A named pipe is not a regular file, so the `fstat` must route it to the stream tier
    /// rather than let `mmap` fail on it. This is what a process substitution looks like
    /// from the inside.
    #[cfg(unix)]
    #[test]
    fn a_named_pipe_opens_as_a_stream() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("pipe.fasta");
        let created = std::process::Command::new("mkfifo")
            .arg(&path)
            .status()
            .unwrap();
        assert!(created.success(), "mkfifo failed");

        let plain = sample_fasta(4);
        let expected = headers_in(&plain);
        // Both ends block until the other arrives, so the writer runs alongside the open.
        let written_path = path.clone();
        let writer = std::thread::spawn(move || std::fs::write(&written_path, &plain).unwrap());

        let mut input = InputFile::open(Some(path.to_str().unwrap())).unwrap();
        assert!(matches!(input, InputFile::Stream(_)));
        let mut headers = Vec::new();
        input
            .for_each_record(|record| {
                headers.push(record.header.to_vec());
                Ok(())
            })
            .unwrap();
        writer.join().unwrap();
        assert_eq!(headers, expected);
    }

    /// One BGZF block over `payload`: gzip magic, a `BC` extra field, raw deflate, trailer.
    #[cfg(feature = "gzip")]
    fn bgzf_block(payload: &[u8]) -> Vec<u8> {
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

    /// `plain` as BGZF, in blocks small enough that records straddle them.
    #[cfg(feature = "gzip")]
    fn bgzf_bytes(plain: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        for chunk in plain.chunks(512) {
            out.extend_from_slice(&bgzf_block(chunk));
        }
        out.extend_from_slice(&bgzf_block(b""));
        out
    }

    /// `plain` as one plain gzip member, which carries no block layout.
    #[cfg(feature = "gzip")]
    fn gzip_bytes(plain: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, plain).unwrap();
        encoder.finish().unwrap()
    }

    /// BGZF on a regular file stays block-addressable, which is why the tier exists.
    #[cfg(feature = "gzip")]
    #[test]
    fn a_bgzip_file_opens_as_blocks() {
        let directory = tempfile::tempdir().unwrap();
        let path = written(
            &directory,
            "blocked.fasta.gz",
            &bgzf_bytes(&sample_fasta(8)),
        );
        assert!(matches!(
            InputFile::open(Some(path.to_str().unwrap())).unwrap(),
            InputFile::Blocks(_)
        ));
        assert_eq!(records_at(&path).unwrap().len(), 8);
    }

    /// Plain gzip has no block layout, so it reads forward only however it is stored.
    #[cfg(feature = "gzip")]
    #[test]
    fn a_plain_gzip_file_opens_as_stream() {
        let directory = tempfile::tempdir().unwrap();
        let path = written(
            &directory,
            "streamed.fasta.gz",
            &gzip_bytes(&sample_fasta(8)),
        );
        assert!(matches!(
            InputFile::open(Some(path.to_str().unwrap())).unwrap(),
            InputFile::Stream(_)
        ));
        assert_eq!(records_at(&path).unwrap().len(), 8);
    }

    /// The tier is an implementation detail, so all three shapes of the same file must yield
    /// exactly the same records.
    #[cfg(feature = "gzip")]
    #[test]
    fn every_input_shape_yields_the_same_records() {
        let plain = sample_fasta(64);
        let directory = tempfile::tempdir().unwrap();
        let expected = records_at(&written(&directory, "plain.fasta", &plain)).unwrap();
        assert_eq!(expected.len(), 64);

        let blocked = written(&directory, "blocked.fasta.gz", &bgzf_bytes(&plain));
        assert_eq!(records_at(&blocked).unwrap(), expected, "blocked diverged");

        let streamed = written(&directory, "streamed.fasta.gz", &gzip_bytes(&plain));
        assert_eq!(
            records_at(&streamed).unwrap(),
            expected,
            "streamed diverged"
        );
    }
}
