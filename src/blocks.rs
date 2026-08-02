//! Block-addressable containers: BGZF, xz, and lz4.
//!
//! A container is block-addressable when its compressed bytes fall into runs that each decode
//! alone, which is what lets several workers read one compressed input. The three formats
//! answer the same three questions in three different ways — BGZF walks a chain of block
//! headers, xz reads a footer index, lz4 walks a chain of length prefixes — so [`BlockAccess`]
//! is a trait with one implementation each.
//!
//! Nothing here parses records. [`BlockAccess::reader`] hands back a plain [`Read`], so a
//! record spanning two blocks is invisible to this module and needs no handling — to the
//! parser it is a buffer that ended mid-record, exactly as a record spanning two pipe reads is.
//!
//! __Integrity is not uniform across the three.__ A BGZF trailer and an xz index record each
//! state what their block decodes to, so a truncated payload is an error where it is decoded.
//! lz4 states no decoded size anywhere, so a payload cut in half decodes to `Ok` with half the
//! bytes; its only structural guard is the chain of length prefixes that locates the blocks.

use std::io::{self, Read};
use std::ops::Deref;
#[cfg(any(feature = "gzip", feature = "lz4"))]
use std::sync::OnceLock;

use crate::codecs::Codec;

// region: The block tier

/// A container addressable at block granularity, where every block decodes alone.
///
/// `Sync` because a worker pool queries one shared instance while each worker decodes into its
/// own buffer. Not implemented for streams: a forward-only source has no blocks.
pub trait BlockAccess: Sync {
    /// How many blocks the input holds.
    ///
    /// A walk that meets a malformed block stops there and counts what validated, so a
    /// truncated tail is missing from the count rather than an error at open.
    fn block_count(&self) -> usize;

    /// Decoded bytes in one block, as far as the container states them.
    ///
    /// Decoded rather than compressed because poly-A compresses ten to one against ordinary
    /// sequence, so budgeting on compressed bytes hands one worker ten times the work.
    ///
    /// # Panics
    ///
    /// When `index` is at or past [`BlockAccess::block_count`]. A block that is not in the
    /// table has no size to report, and a caller that asks for one has counted wrong.
    fn bytes_in_block(&self, index: usize) -> BlockBytes;

    /// A decode scratch for one worker, separate from `&self` because the container is shared
    /// and the decoder is not.
    fn decoder(&self) -> Box<dyn BlockDecoder + '_>;

    /// A forward reader over `from_block` and every block after it.
    ///
    /// Reads past a caller's last block deliberately: a record straddling that boundary began
    /// inside the caller's run and has to be finished. A `from_block` at or past the last
    /// block gives an empty reader, which is what makes an empty run read as empty rather than
    /// as the whole input.
    fn reader(&self, from_block: usize) -> Box<dyn Read + '_> {
        Box::new(BlockRunReader {
            blocks: from_block..self.block_count(),
            decoder: self.decoder(),
            decoded: Vec::new(),
            handed_out: 0,
        })
    }
}

/// What a container knows about a block's decoded size before decoding it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockBytes {
    /// The container declares it: a BGZF trailer, an xz index record, an lz4 stored block.
    Exact(usize),
    /// A ceiling only: an lz4 compressed block, whose frame states no decoded size at all.
    AtMost(usize),
}

impl BlockBytes {
    /// For a scheduler sizing a budget, where over-estimating costs memory and nothing else.
    pub fn upper_bound(self) -> usize {
        match self {
            BlockBytes::Exact(bytes) | BlockBytes::AtMost(bytes) => bytes,
        }
    }

    /// For a caller that needs a number it can report or slice on, which a ceiling is not.
    pub fn exact(self) -> Option<usize> {
        match self {
            BlockBytes::Exact(bytes) => Some(bytes),
            BlockBytes::AtMost(_) => None,
        }
    }
}

/// One worker's decoder over one container, holding whatever state its format needs.
pub trait BlockDecoder {
    /// Decode block `index` into `out`, replacing what was there.
    ///
    /// # Panics
    ///
    /// When `index` is at or past the container's block count, for the same reason
    /// [`BlockAccess::bytes_in_block`] does.
    fn decode_into(&mut self, index: usize, out: &mut Vec<u8>) -> io::Result<()>;
}

/// A forward reader over a run of blocks, decoding one at a time.
///
/// Serves every block from the one it opened on to the end of the input, so a caller owning
/// only some of them stops reading once it has taken the decoded bytes those blocks hold.
struct BlockRunReader<'a> {
    blocks: core::ops::Range<usize>,
    decoder: Box<dyn BlockDecoder + 'a>,
    decoded: Vec<u8>,
    handed_out: usize,
}

impl Read for BlockRunReader<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        // An empty block — a BGZF end-of-file marker, or a flush of an empty buffer — decodes
        // to nothing, and handing back zero would read as end of stream, so keep decoding
        // until there is something to give.
        while self.handed_out == self.decoded.len() {
            let Some(index) = self.blocks.next() else {
                return Ok(0);
            };
            self.decoder.decode_into(index, &mut self.decoded)?;
            self.handed_out = 0;
        }
        let take = out.len().min(self.decoded.len() - self.handed_out);
        out[..take].copy_from_slice(&self.decoded[self.handed_out..self.handed_out + take]);
        self.handed_out += take;
        Ok(take)
    }
}

/// A block-addressable view of `compressed`, or the bytes back when it has no block layout
/// this build can address.
///
/// The second question detection asks. [`Codec::detect`] names the container from a prefix and
/// that settles it for BGZF and lz4, while a plain xz stream and a blocked one share every byte
/// of their magic and only the footer index tells them apart. Three inputs come back as the
/// bytes they went in as and read as streams: a single-block xz file, an lz4 frame whose blocks
/// depend on the window before them, and any container with no block layout at all.
pub fn block_access<B: Deref<Target = [u8]> + Sync + 'static>(
    compressed: B,
    codec: Codec,
) -> Result<Box<dyn BlockAccess>, B> {
    match codec {
        #[cfg(feature = "gzip")]
        Codec::Bgzf => BgzfAccess::new(compressed).map(|access| Box::new(access) as Box<_>),
        #[cfg(feature = "xz")]
        Codec::Xz => XzAccess::new(compressed).map(|access| Box::new(access) as Box<_>),
        #[cfg(feature = "lz4")]
        Codec::Lz4 => Lz4Access::new(compressed).map(|access| Box::new(access) as Box<_>),
        _ => Err(compressed),
    }
}

// endregion: The block tier

// region: BGZF

/// The first four bytes of every BGZF block: gzip magic, deflate, and `FEXTRA` set.
///
/// Bytes four through nine are the modification time, extra flags, and operating system, which
/// `bgzip` zeroes but the gzip format lets vary, so they are not checked.
#[cfg(feature = "gzip")]
pub const BGZF_HEADER_PREFIX: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Total size of the BGZF block at the front of `block`, from the `BC` extra subfield.
///
/// Nothing is inflated; only the first eighteen bytes are touched. `None` means the bytes do
/// not open a BGZF block, which is not an error but how a caller tells the containers apart.
#[cfg(feature = "gzip")]
pub fn bgzf_block_size(block: &[u8]) -> Option<usize> {
    if block.len() < 18 || block[..4] != BGZF_HEADER_PREFIX {
        return None;
    }
    let extra_length = u16::from_le_bytes([block[10], block[11]]) as usize;
    let extra = block.get(12..12 + extra_length)?;
    let mut cursor = 0usize;
    // `BC` comes first in every writer seen in the wild, but the format allows other subfields
    // ahead of it, so the extras are walked rather than assumed.
    while cursor + 4 <= extra.len() {
        let payload_length = u16::from_le_bytes([extra[cursor + 2], extra[cursor + 3]]) as usize;
        let payload = extra.get(cursor + 4..cursor + 4 + payload_length)?;
        if extra[cursor] == b'B' && extra[cursor + 1] == b'C' && payload_length == 2 {
            let size = u16::from_le_bytes([payload[0], payload[1]]) as usize + 1;
            // A size not covering header, extras, and the eight-byte trailer would let a chain
            // walk stand still or run backwards.
            return (size >= 12 + extra_length + 8).then_some(size);
        }
        cursor += 4 + payload_length;
    }
    None
}

/// Compressed offset, compressed size, and decoded size of one block, read from its header and
/// its trailer. Nothing is inflated to fill this in.
#[cfg(feature = "gzip")]
struct BgzfBlockEntry {
    compressed_offset: u64,
    compressed_bytes: u32,
    decoded_bytes: u32,
}

/// A BGZF file: blocked gzip, where every block is an independent deflate stream.
///
/// Generic over its container like [`crate::files::RandomAccess`], so one type serves the
/// owner of a mapping and the borrower of a view of it alike.
#[cfg(feature = "gzip")]
pub struct BgzfAccess<B: Deref<Target = [u8]>> {
    compressed: B,
    /// Filled by the first query that needs it, never at open.
    table: OnceLock<Vec<BgzfBlockEntry>>,
}

#[cfg(feature = "gzip")]
impl<B: Deref<Target = [u8]>> BgzfAccess<B> {
    /// Wraps compressed bytes opening with a BGZF block, handing them back when they do not.
    ///
    /// Plain gzip is another container rather than an error, so the rejected bytes come back
    /// and the caller reaches for a different one without mapping the file twice.
    pub fn new(compressed: B) -> Result<Self, B> {
        match bgzf_block_size(&compressed) {
            Some(_) => Ok(Self {
                compressed,
                table: OnceLock::new(),
            }),
            None => Err(compressed),
        }
    }

    /// A two-word view of the same bytes, for a worker that must not own the container.
    ///
    /// The view starts with an empty table of its own, so whichever instance the block queries
    /// go through is the one that pays for the walk.
    pub fn borrowed(&self) -> BgzfAccess<&[u8]> {
        BgzfAccess {
            compressed: &self.compressed,
            table: OnceLock::new(),
        }
    }

    /// Every block's offset and sizes, walked from offset zero on the first query.
    ///
    /// A 200 GB file holds roughly 3.1 million blocks and the walk faults about 13 GB of
    /// headers and trailers, so a run that only reads forward from the front never asks and
    /// never pays.
    fn table(&self) -> &[BgzfBlockEntry] {
        self.table.get_or_init(|| {
            let mut entries = Vec::new();
            let mut offset = 0usize;
            while offset < self.compressed.len() {
                let Some(size) = bgzf_block_size(&self.compressed[offset..]) else {
                    break;
                };
                // A declared size running past the end is a truncated final block, which has
                // no trailer to read and so is not a block this table can name.
                let Some(block) = self.compressed.get(offset..offset + size) else {
                    break;
                };
                entries.push(BgzfBlockEntry {
                    compressed_offset: offset as u64,
                    compressed_bytes: size as u32,
                    // The trailer's last four bytes are `ISIZE`, the decoded length.
                    decoded_bytes: u32::from_le_bytes([
                        block[size - 4],
                        block[size - 3],
                        block[size - 2],
                        block[size - 1],
                    ]),
                });
                offset += size;
            }
            entries
        })
    }

    /// Compressed offset the reader for `from_block` opens on.
    ///
    /// Block zero is the one query that answers without a table, and it is the one a serial
    /// run over a 200 GB file asks.
    fn offset_of_block(&self, from_block: usize) -> usize {
        match from_block {
            0 => 0,
            index => {
                let table = self.table();
                match table.get(index) {
                    Some(entry) => entry.compressed_offset as usize,
                    // Past the last block, so the reader opens on the end of the chain and
                    // gives back nothing.
                    None => table.last().map_or(0, |last| {
                        (last.compressed_offset + u64::from(last.compressed_bytes)) as usize
                    }),
                }
            }
        }
    }
}

#[cfg(feature = "gzip")]
impl<B: Deref<Target = [u8]> + Sync> BlockAccess for BgzfAccess<B> {
    fn block_count(&self) -> usize {
        self.table().len()
    }

    /// Exact for every block, read from its trailer without inflating it.
    fn bytes_in_block(&self, index: usize) -> BlockBytes {
        BlockBytes::Exact(self.table()[index].decoded_bytes as usize)
    }

    fn decoder(&self) -> Box<dyn BlockDecoder + '_> {
        Box::new(BgzfDecoder {
            compressed: &self.compressed,
            table: self.table(),
            inflate: flate2::Decompress::new(false),
        })
    }

    /// Follows the chain rather than the table, because each header states where the next
    /// block begins. Opening at block zero therefore costs nothing, which is what keeps a
    /// serial run over a 200 GB file from walking three million headers first.
    fn reader(&self, from_block: usize) -> Box<dyn Read + '_> {
        Box::new(BgzfReader::new(
            &self.compressed,
            self.offset_of_block(from_block),
        ))
    }
}

/// Deflate state and the block table, so one worker decodes any block of one input.
///
/// The state is reset per block rather than rebuilt, so a run of blocks costs one allocation.
#[cfg(feature = "gzip")]
struct BgzfDecoder<'a> {
    compressed: &'a [u8],
    table: &'a [BgzfBlockEntry],
    inflate: flate2::Decompress,
}

#[cfg(feature = "gzip")]
impl BlockDecoder for BgzfDecoder<'_> {
    fn decode_into(&mut self, index: usize, out: &mut Vec<u8>) -> io::Result<()> {
        let offset = self.table[index].compressed_offset as usize;
        inflate_bgzf_block(&mut self.inflate, self.compressed, offset, out)?;
        Ok(())
    }
}

/// Inflate the block at compressed `offset` into `out`, replacing what was there, and give back
/// the block's compressed size so a walker can step to the next one.
#[cfg(feature = "gzip")]
fn inflate_bgzf_block(
    inflate: &mut flate2::Decompress,
    compressed: &[u8],
    offset: usize,
    out: &mut Vec<u8>,
) -> io::Result<usize> {
    let block = &compressed[offset..];
    let size = bgzf_block_size(block).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("not a BGZF block at compressed offset {offset}"),
        )
    })?;
    let extra_length = u16::from_le_bytes([block[10], block[11]]) as usize;
    let payload = &block[12 + extra_length..size - 8];
    // The trailer's last four bytes are the decoded size, so the buffer is sized exactly.
    let decoded_size = u32::from_le_bytes([
        block[size - 4],
        block[size - 3],
        block[size - 2],
        block[size - 1],
    ]) as usize;

    out.clear();
    out.reserve(decoded_size);
    inflate.reset(false);
    inflate
        .decompress_vec(payload, out, flate2::FlushDecompress::Finish)
        .map_err(|error| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("corrupt BGZF block at compressed offset {offset}: {error}"),
            )
        })?;
    // The block's trailer states what it should have produced. Checking turns a short inflate
    // into an error rather than a block that quietly decodes to nothing.
    if out.len() != decoded_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF block at compressed offset {offset} declares {decoded_size} bytes but \
                 inflated {}",
                out.len()
            ),
        ));
    }
    Ok(size)
}

/// A forward reader over BGZF blocks, inflating one at a time and following the chain.
///
/// Serves every block from the one it opened on to the end of the file. A caller owning only
/// some of them stops reading once it has taken the decoded bytes those blocks hold, which the
/// trailers state before anything is inflated.
#[cfg(feature = "gzip")]
struct BgzfReader<'a> {
    /// The whole file, so the reader can run past a run of blocks to finish a record.
    compressed: &'a [u8],
    /// Compressed offset of the block to inflate next.
    cursor: usize,
    decoded: Vec<u8>,
    handed_out: usize,
    inflate: flate2::Decompress,
}

#[cfg(feature = "gzip")]
impl<'a> BgzfReader<'a> {
    fn new(compressed: &'a [u8], from_offset: usize) -> Self {
        Self {
            compressed,
            cursor: from_offset.min(compressed.len()),
            decoded: Vec::new(),
            handed_out: 0,
            inflate: flate2::Decompress::new(false),
        }
    }
}

#[cfg(feature = "gzip")]
impl Read for BgzfReader<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        // An empty block — the end-of-file marker, or a flush of an empty buffer — decodes to
        // nothing, and handing back zero would read as end of stream, so keep inflating until
        // there is something to give.
        while self.handed_out == self.decoded.len() {
            if self.cursor >= self.compressed.len() {
                return Ok(0);
            }
            let size = inflate_bgzf_block(
                &mut self.inflate,
                self.compressed,
                self.cursor,
                &mut self.decoded,
            )?;
            self.handed_out = 0;
            self.cursor += size;
        }
        let take = out.len().min(self.decoded.len() - self.handed_out);
        out[..take].copy_from_slice(&self.decoded[self.handed_out..self.handed_out + take]);
        self.handed_out += take;
        Ok(take)
    }
}

// endregion: BGZF

// region: xz

/// The six bytes every xz stream opens with.
#[cfg(feature = "xz")]
const XZ_HEADER_MAGIC: [u8; 6] = [0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00];

/// The two bytes every xz stream ends with, ASCII `YZ`.
#[cfg(feature = "xz")]
const XZ_FOOTER_MAGIC: [u8; 2] = [0x59, 0x5A];

/// One xz block, located and sized from the stream index alone.
#[cfg(feature = "xz")]
#[derive(Clone, Copy)]
struct XzBlockEntry {
    /// Offset of the block header, which is where the padded run begins.
    compressed_offset: usize,
    /// Offset of the LZMA2 payload, past a header whose width varies per block.
    payload_offset: usize,
    /// End of the payload, before the check bytes.
    payload_end: usize,
    /// Exact, from the index record, which carries it for every block.
    uncompressed_size: usize,
    /// Window the payload was coded against, from the block header's LZMA2 property byte.
    dictionary_size: u32,
}

/// An xz file: one or more streams, each an index over blocks that decode alone.
///
/// The table is eager because deciding this file is block-addressable already requires it. A
/// plain stream and a blocked one share every byte of their magic, and only the count of index
/// records in the footer tells them apart.
#[cfg(feature = "xz")]
pub struct XzAccess<B: Deref<Target = [u8]>> {
    compressed: B,
    table: Vec<XzBlockEntry>,
}

#[cfg(feature = "xz")]
impl<B: Deref<Target = [u8]>> XzAccess<B> {
    /// Wraps compressed bytes holding more than one addressable block, handing them back when
    /// they do not.
    ///
    /// A one-block file comes back rather than being accepted, because a single-entry table
    /// buys nothing the stream tier does not already do and reads faster there.
    pub fn new(compressed: B) -> Result<Self, B> {
        match xz_block_table(&compressed) {
            Some(table) if table.len() > 1 => Ok(Self { compressed, table }),
            _ => Err(compressed),
        }
    }
}

#[cfg(feature = "xz")]
impl<B: Deref<Target = [u8]> + Sync> BlockAccess for XzAccess<B> {
    fn block_count(&self) -> usize {
        self.table.len()
    }

    /// Exact for every block, from the index record that names it.
    fn bytes_in_block(&self, index: usize) -> BlockBytes {
        BlockBytes::Exact(self.table[index].uncompressed_size)
    }

    fn decoder(&self) -> Box<dyn BlockDecoder + '_> {
        Box::new(XzDecoder {
            compressed: &self.compressed,
            table: &self.table,
        })
    }
}

/// The block table and the bytes it addresses, which is all an LZMA2 decode needs.
///
/// No decoder state rides along: the LZMA2 reader exposes neither reset nor reuse, and building
/// one costs 1.4 µs against roughly 11 ms to decode a one-mebibyte block.
#[cfg(feature = "xz")]
struct XzDecoder<'a> {
    compressed: &'a [u8],
    table: &'a [XzBlockEntry],
}

#[cfg(feature = "xz")]
impl BlockDecoder for XzDecoder<'_> {
    fn decode_into(&mut self, index: usize, out: &mut Vec<u8>) -> io::Result<()> {
        let entry = self.table[index];
        let payload = &self.compressed[entry.payload_offset..entry.payload_end];
        let mut reader = lzma_rust2::Lzma2Reader::new(payload, entry.dictionary_size, None);
        out.clear();
        out.reserve(entry.uncompressed_size);
        reader.read_to_end(out)?;
        // The index states what the block owes. Checking turns a truncated file into an error
        // rather than a run of blocks that quietly drops records.
        if out.len() != entry.uncompressed_size {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "xz block {index} at compressed offset {} declares {} bytes but decoded {}",
                    entry.compressed_offset,
                    entry.uncompressed_size,
                    out.len()
                ),
            ));
        }
        Ok(())
    }
}

/// Little-endian base-128, high bit meaning another byte follows, nine bytes at most.
#[cfg(feature = "xz")]
fn read_variable_integer(bytes: &[u8], cursor: &mut usize) -> Option<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    for _ in 0..9 {
        let byte = *bytes.get(*cursor)?;
        *cursor += 1;
        value |= u64::from(byte & 0x7F) << shift;
        if byte & 0x80 == 0 {
            return Some(value);
        }
        shift += 7;
    }
    None
}

/// Payload bounds and dictionary size for the block whose header starts at `header_offset`.
///
/// The property byte moves. A threaded writer states both optional sizes and a single-threaded
/// one states neither, so the offset differs between two blocks of one file, and a hardcoded
/// one misparses every file written with more than one thread. Anything but a lone LZMA2 filter
/// is rejected rather than guessed at.
#[cfg(feature = "xz")]
fn xz_payload_bounds(
    compressed: &[u8],
    header_offset: usize,
    unpadded_size: usize,
    check_bytes: usize,
) -> Option<(usize, usize, u32)> {
    let header_size = (*compressed.get(header_offset)? as usize + 1) * 4;
    let header = compressed.get(header_offset..header_offset + header_size)?;
    let flags = header[1];
    // Reserved bits set, or a filter chain this reader will not guess at.
    if flags & 0x3C != 0 || flags & 0x03 != 0 {
        return None;
    }
    let mut cursor = 2usize;
    if flags & 0x40 != 0 {
        read_variable_integer(header, &mut cursor)?;
    }
    if flags & 0x80 != 0 {
        read_variable_integer(header, &mut cursor)?;
    }
    if read_variable_integer(header, &mut cursor)? != 0x21 {
        return None; // a filter that is not LZMA2
    }
    if read_variable_integer(header, &mut cursor)? != 1 {
        return None; // LZMA2 carries exactly one property byte
    }
    let property = u32::from(*header.get(cursor)?);
    let dictionary_size = match property {
        40 => u32::MAX,
        41.. => return None,
        _ => (2 | (property & 1)) << (property / 2 + 11),
    };
    let payload_offset = header_offset + header_size;
    let payload_end = header_offset + unpadded_size.checked_sub(check_bytes)?;
    (payload_end >= payload_offset).then_some((payload_offset, payload_end, dictionary_size))
}

/// Every block of every stream, walking backwards from the last byte.
///
/// Backwards because only the footer states where the index begins, and only the index states
/// where the blocks are. `None` means the bytes are not an xz file this reader can address,
/// which is how the caller falls back to the stream tier rather than failing the run.
#[cfg(feature = "xz")]
fn xz_block_table(compressed: &[u8]) -> Option<Vec<XzBlockEntry>> {
    if compressed.len() < 32 || compressed[..6] != XZ_HEADER_MAGIC {
        return None;
    }
    let mut streams: Vec<Vec<XzBlockEntry>> = Vec::new();
    let mut tail = compressed.len();
    while tail > 0 {
        // Stream padding is any number of four-byte zero groups between streams.
        while tail >= 4 && compressed[tail - 4..tail] == [0, 0, 0, 0] {
            tail -= 4;
        }
        if tail == 0 {
            break;
        }
        if tail < 24 || compressed[tail - 2..tail] != XZ_FOOTER_MAGIC {
            return None;
        }
        let stream_flags = [compressed[tail - 4], compressed[tail - 3]];
        // The check trails every block's payload, and getting its width wrong shifts every
        // payload end in the stream.
        let check_bytes = match stream_flags[1] & 0x0F {
            0 => 0usize,
            1..=3 => 4,
            4..=6 => 8,
            7..=9 => 16,
            10..=12 => 32,
            _ => 64,
        };
        let backward_size = u32::from_le_bytes([
            compressed[tail - 8],
            compressed[tail - 7],
            compressed[tail - 6],
            compressed[tail - 5],
        ]) as usize;
        let index_size = (backward_size + 1) * 4;
        let index_start = tail.checked_sub(12 + index_size)?;
        if compressed[index_start] != 0x00 {
            return None;
        }

        let mut cursor = index_start + 1;
        let record_count = read_variable_integer(compressed, &mut cursor)? as usize;
        let mut records = Vec::new();
        let mut blocks_span = 0usize;
        for _ in 0..record_count {
            let unpadded_size = read_variable_integer(compressed, &mut cursor)? as usize;
            let uncompressed_size = read_variable_integer(compressed, &mut cursor)? as usize;
            if unpadded_size == 0 {
                return None;
            }
            records.push((unpadded_size, uncompressed_size));
            blocks_span += (unpadded_size + 3) & !3;
        }

        // The stream header sits exactly one span plus twelve bytes ahead of the index.
        // Asserting it is what proves this index describes these bytes rather than another
        // file's, and it is what catches a second stream concatenated ahead of this one.
        let stream_start = index_start.checked_sub(blocks_span + 12)?;
        if compressed.get(stream_start..stream_start + 6)? != XZ_HEADER_MAGIC
            || compressed[stream_start + 6..stream_start + 8] != stream_flags
        {
            return None;
        }

        let mut entries = Vec::with_capacity(records.len());
        let mut header_offset = stream_start + 12;
        for (unpadded_size, uncompressed_size) in records {
            let (payload_offset, payload_end, dictionary_size) =
                xz_payload_bounds(compressed, header_offset, unpadded_size, check_bytes)?;
            entries.push(XzBlockEntry {
                compressed_offset: header_offset,
                payload_offset,
                payload_end,
                uncompressed_size,
                dictionary_size,
            });
            header_offset += (unpadded_size + 3) & !3;
        }
        streams.push(entries);
        tail = stream_start;
    }
    streams.reverse();
    Some(streams.concat())
}

// endregion: xz

// region: lz4

/// Magic of a frame carrying blocks, little endian.
#[cfg(feature = "lz4")]
const LZ4_FRAME_MAGIC: u32 = 0x184D_2204;

/// Magic range of a skippable frame, which states its own length and carries nothing this
/// reader looks at.
#[cfg(feature = "lz4")]
const LZ4_SKIPPABLE_MAGIC: core::ops::RangeInclusive<u32> = 0x184D_2A50..=0x184D_2A5F;

/// One lz4 block, located by the forward walk of length prefixes.
#[cfg(feature = "lz4")]
#[derive(Clone, Copy)]
struct Lz4BlockEntry {
    /// First byte of the block payload, past its four-byte length prefix.
    payload_offset: usize,
    /// Bytes on disk, which is the decoded count when the payload is stored uncompressed.
    stored_size: usize,
    /// Set when the prefix's high bit says the payload was never compressed.
    stored_uncompressed: bool,
    /// Frame maximum of this block's own frame, since a file may hold several frames and each
    /// states its own.
    frame_maximum: usize,
}

/// An lz4 file: one or more frames, each a chain of length-prefixed blocks.
///
/// The table is deferred, because the walk touches four bytes per block across the whole file
/// while addressability is answerable from the eleven bytes of the first frame descriptor.
///
/// __No block carries an integrity check unless the writer asked for one.__ A payload cut in
/// half decodes to `Ok` with half the bytes rather than to an error, since a raw lz4 block has
/// no end sentinel and no declared decoded size. What guards the file is the chain of length
/// prefixes: a corrupt region breaks the chain, the walk stops there, and the blocks past it
/// are missing from the count.
#[cfg(feature = "lz4")]
pub struct Lz4Access<B: Deref<Target = [u8]>> {
    compressed: B,
    table: OnceLock<Vec<Lz4BlockEntry>>,
}

#[cfg(feature = "lz4")]
impl<B: Deref<Target = [u8]>> Lz4Access<B> {
    /// Wraps compressed bytes opening a frame whose blocks decode alone, handing them back when
    /// they do not.
    ///
    /// Only the frame descriptor is read here, so a serial run never pays the walk. The
    /// descriptor is enough: block independence and the frame maximum both live in it.
    pub fn new(compressed: B) -> Result<Self, B> {
        let addressable = lz4_first_frame_offset(&compressed)
            .is_some_and(|offset| lz4_frame_is_addressable(&compressed[offset + 4..]));
        match addressable {
            true => Ok(Self {
                compressed,
                table: OnceLock::new(),
            }),
            false => Err(compressed),
        }
    }

    /// The block table, walked once on first use and shared thereafter.
    fn table(&self) -> &[Lz4BlockEntry] {
        self.table.get_or_init(|| lz4_block_table(&self.compressed))
    }
}

#[cfg(feature = "lz4")]
impl<B: Deref<Target = [u8]> + Sync> BlockAccess for Lz4Access<B> {
    fn block_count(&self) -> usize {
        self.table().len()
    }

    /// A ceiling for a compressed block and exact for a stored one, because the frame states no
    /// decoded size and the format lets a writer emit a one-byte block under a four-mebibyte
    /// maximum.
    fn bytes_in_block(&self, index: usize) -> BlockBytes {
        let entry = self.table()[index];
        match entry.stored_uncompressed {
            true => BlockBytes::Exact(entry.stored_size),
            false => BlockBytes::AtMost(entry.frame_maximum),
        }
    }

    fn decoder(&self) -> Box<dyn BlockDecoder + '_> {
        Box::new(Lz4Decoder {
            compressed: &self.compressed,
            table: self.table(),
        })
    }
}

/// The block table and the bytes it addresses, which is all a raw lz4 decode needs.
///
/// No decoder state rides along, because the raw block decoder is a free function.
#[cfg(feature = "lz4")]
struct Lz4Decoder<'a> {
    compressed: &'a [u8],
    table: &'a [Lz4BlockEntry],
}

#[cfg(feature = "lz4")]
impl BlockDecoder for Lz4Decoder<'_> {
    /// `out` grows to the frame maximum and is truncated to what the block actually gave,
    /// because the frame states no decoded size for any block.
    fn decode_into(&mut self, index: usize, out: &mut Vec<u8>) -> io::Result<()> {
        let entry = self.table[index];
        let payload =
            &self.compressed[entry.payload_offset..entry.payload_offset + entry.stored_size];
        if entry.stored_uncompressed {
            out.clear();
            out.extend_from_slice(payload);
            return Ok(());
        }
        out.resize(entry.frame_maximum, 0);
        let produced = lz4_flex::block::decompress_into(payload, out).map_err(|error| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "corrupt lz4 block {index} at compressed offset {}: {error}",
                    entry.payload_offset
                ),
            )
        })?;
        out.truncate(produced);
        Ok(())
    }
}

/// Decoded ceiling of a frame whose block-descriptor byte is `sizes`, or `None` when the code is
/// one no writer this reader trusts emits.
#[cfg(feature = "lz4")]
fn lz4_frame_maximum(sizes: u8) -> Option<usize> {
    if sizes & 0x8F != 0 {
        return None; // reserved high bit or low nibble set
    }
    match (sizes >> 4) & 7 {
        4 => Some(64 << 10),
        5 => Some(256 << 10),
        6 => Some(1 << 20),
        7 => Some(4 << 20),
        _ => None,
    }
}

/// Whether a frame descriptor, starting at its flags byte, opens blocks that decode alone.
///
/// Block independence is the whole question. A writer told to link its blocks produces a
/// perfectly valid frame whose every block references the window before it, so none of them
/// decodes alone and the file belongs on the stream tier.
#[cfg(feature = "lz4")]
fn lz4_frame_is_addressable(descriptor: &[u8]) -> bool {
    let (Some(&flags), Some(&sizes)) = (descriptor.first(), descriptor.get(1)) else {
        return false;
    };
    flags >> 6 == 1 && flags & 0x02 == 0 && flags & 0x20 != 0 && lz4_frame_maximum(sizes).is_some()
}

/// Offset of the first frame carrying blocks, stepping over any skippable frames ahead of it.
#[cfg(feature = "lz4")]
fn lz4_first_frame_offset(compressed: &[u8]) -> Option<usize> {
    let mut cursor = 0usize;
    loop {
        let magic = u32::from_le_bytes(compressed.get(cursor..cursor + 4)?.try_into().ok()?);
        if !LZ4_SKIPPABLE_MAGIC.contains(&magic) {
            return (magic == LZ4_FRAME_MAGIC).then_some(cursor);
        }
        let skipped = u32::from_le_bytes(compressed.get(cursor + 4..cursor + 8)?.try_into().ok()?);
        cursor = cursor.checked_add(8)?.checked_add(skipped as usize)?;
    }
}

/// Every block of every frame, stepping over skippable frames.
///
/// The walk stops at the first thing it cannot address — a legacy frame, a frame of linked
/// blocks, a chain running off the end — and keeps the blocks it validated, so a truncated file
/// yields its intact prefix rather than nothing.
#[cfg(feature = "lz4")]
fn lz4_block_table(compressed: &[u8]) -> Vec<Lz4BlockEntry> {
    let mut table = Vec::new();
    let mut cursor = 0usize;
    while cursor < compressed.len() {
        match lz4_walk_frame(compressed, cursor, &mut table) {
            Some(next) => cursor = next,
            None => break,
        }
    }
    table
}

/// Append one frame's blocks to `table` and give back the offset just past the frame.
#[cfg(feature = "lz4")]
fn lz4_walk_frame(compressed: &[u8], from: usize, table: &mut Vec<Lz4BlockEntry>) -> Option<usize> {
    let mut cursor = from;
    let magic = u32::from_le_bytes(compressed.get(cursor..cursor + 4)?.try_into().ok()?);
    if LZ4_SKIPPABLE_MAGIC.contains(&magic) {
        // A skippable frame states its own length and carries nothing this reader wants.
        let skipped = u32::from_le_bytes(compressed.get(cursor + 4..cursor + 8)?.try_into().ok()?);
        return cursor.checked_add(8)?.checked_add(skipped as usize);
    }
    if magic != LZ4_FRAME_MAGIC {
        return None; // a legacy frame, trailing junk, or not lz4 at all
    }

    let descriptor = compressed.get(cursor + 4..)?;
    if !lz4_frame_is_addressable(descriptor) {
        return None;
    }
    let flags = descriptor[0];
    let frame_maximum = lz4_frame_maximum(descriptor[1])?;
    let block_checksum = flags & 0x10 != 0;

    // The optional fields move the first prefix: offset seven bare, fifteen with a content size,
    // and four further with a dictionary identifier.
    cursor += 7 + 8 * usize::from(flags & 0x08 != 0) + 4 * usize::from(flags & 0x01 != 0);

    loop {
        let prefix = u32::from_le_bytes(compressed.get(cursor..cursor + 4)?.try_into().ok()?);
        cursor += 4;
        if prefix == 0 {
            break; // the end mark, after which another frame may follow
        }
        let stored_size = (prefix & 0x7FFF_FFFF) as usize;
        if stored_size > frame_maximum {
            return None; // a block may never exceed its frame's maximum
        }
        let payload_end = cursor.checked_add(stored_size)?;
        if payload_end > compressed.len() {
            return None; // the chain ran off the end
        }
        table.push(Lz4BlockEntry {
            payload_offset: cursor,
            stored_size,
            stored_uncompressed: prefix & 0x8000_0000 != 0,
            frame_maximum,
        });
        cursor = payload_end.checked_add(4 * usize::from(block_checksum))?;
        if cursor > compressed.len() {
            return None;
        }
    }
    // A content checksum trails the end mark when the descriptor asked for one.
    cursor.checked_add(4 * usize::from(flags & 0x04 != 0))
}

// endregion: lz4

#[cfg(test)]
mod tests {
    #![allow(dead_code)]

    use super::*;

    // region: What every implementation owes

    /// FASTQ of ragged records, so a record straddles a block boundary far more often than not.
    fn ragged_fastq(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            let width = 8 + (index % 5) * 4;
            data.extend_from_slice(format!("@read{index} sample\n").as_bytes());
            data.extend(std::iter::repeat_n(b"ACGT"[index % 4], width));
            data.extend_from_slice(b"\n+\n");
            data.extend(std::iter::repeat_n(b'I', width));
            data.push(b'\n');
        }
        data
    }

    /// FASTQ of records exactly thirty-two bytes wide, so a block boundary can land exactly on a
    /// record start. A ragged fixture never does, and a budget off by one still passes there.
    fn aligned_fastq(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            data.extend_from_slice(format!("@read{index:06}\nACGTACGT\n+\nIIIIIIII\n").as_bytes());
        }
        data
    }

    /// Headers of every record a plain mapping yields, which is what every tier owes.
    fn headers_in_bytes(plain: &[u8]) -> Vec<Vec<u8>> {
        let mut headers = Vec::new();
        crate::files::RandomAccess::new(plain)
            .for_each_record(|record| {
                headers.push(record.header.to_vec());
                Ok(())
            })
            .unwrap();
        headers
    }

    /// Headers of every record the block tier yields, read forward from block zero.
    fn headers_in_blocks(access: &dyn BlockAccess) -> Vec<Vec<u8>> {
        let mut headers = Vec::new();
        crate::files::ForwardAccess::new(access.reader(0))
            .for_each_record(|record| {
                headers.push(record.header.to_vec());
                Ok(())
            })
            .unwrap();
        headers
    }

    /// The contract, checked the same way for all three containers.
    ///
    /// Three claims: the records are the ones a plain mapping yields, the blocks tile the
    /// decoded bytes with no gap and no repeat, and the budget a scheduler sums out of
    /// `bytes_in_block` covers what the input decodes to.
    fn conforms(access: &dyn BlockAccess, plain: &[u8], label: &str) {
        assert!(
            access.block_count() > 1,
            "{label}: the fixture must hold several blocks"
        );
        assert_eq!(
            headers_in_blocks(access),
            headers_in_bytes(plain),
            "{label}: records"
        );

        let mut decoder = access.decoder();
        let mut decoded = Vec::new();
        let mut rejoined = Vec::new();
        let mut budget = 0usize;
        for index in 0..access.block_count() {
            decoder.decode_into(index, &mut decoded).unwrap();
            let stated = access.bytes_in_block(index);
            assert!(
                decoded.len() <= stated.upper_bound(),
                "{label}: block {index} decoded past its ceiling"
            );
            if let Some(exact) = stated.exact() {
                assert_eq!(exact, decoded.len(), "{label}: block {index} declared size");
            }
            budget += stated.upper_bound();
            rejoined.extend_from_slice(&decoded);
        }
        assert_eq!(rejoined, plain, "{label}: blocks must tile the input");
        assert!(
            budget >= plain.len(),
            "{label}: the budget must cover the decoded length"
        );
    }

    /// A run opening at or past the last block has nothing to give, which is what makes an empty
    /// run read as empty rather than as the whole input.
    fn an_empty_run_reads_nothing(access: &dyn BlockAccess, label: &str) {
        for from_block in [access.block_count(), access.block_count() + 4] {
            let mut out = Vec::new();
            access.reader(from_block).read_to_end(&mut out).unwrap();
            assert!(out.is_empty(), "{label}: a run from block {from_block}");
        }
    }

    // endregion: What every implementation owes

    // region: BGZF

    /// Write `payload` as one BGZF block.
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
        let mut out = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0];
        out.extend_from_slice(b"BC");
        out.extend_from_slice(&2u16.to_le_bytes());
        out.extend_from_slice(&((size - 1) as u16).to_le_bytes());
        out.extend_from_slice(&deflated);
        out.extend_from_slice(&checksum.sum().to_le_bytes());
        out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        out
    }

    /// A file of deliberately tiny blocks, so a record straddling a block is the common case
    /// rather than a rare one — the same trick `ShortReader` plays for pipes.
    #[cfg(feature = "gzip")]
    fn bgzf_bytes(plain: &[u8], payload_bytes: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for chunk in plain.chunks(payload_bytes) {
            out.extend_from_slice(&bgzf_block(chunk));
        }
        out.extend_from_slice(&bgzf_block(b""));
        out
    }

    /// Everything a reader opening at `from_block` inflates.
    #[cfg(feature = "gzip")]
    fn bgzf_decoded_from(compressed: &[u8], from_block: usize) -> Vec<u8> {
        let access = BgzfAccess::new(compressed).unwrap();
        let mut out = Vec::new();
        access.reader(from_block).read_to_end(&mut out).unwrap();
        out
    }

    /// Whatever the block layout, the bytes must come back exactly as they went in.
    #[cfg(feature = "gzip")]
    #[test]
    fn bgzf_blocks_inflate_to_the_original_bytes() {
        let plain = ragged_fastq(200);
        for payload_bytes in [64, 512, 4096, 1 << 16] {
            let compressed = bgzf_bytes(&plain, payload_bytes);
            assert_eq!(
                bgzf_decoded_from(&compressed, 0),
                plain,
                "payload {payload_bytes}"
            );
        }
    }

    /// An empty block decodes to nothing, and must not read as end of stream.
    #[cfg(feature = "gzip")]
    #[test]
    fn an_empty_bgzf_block_mid_file_does_not_end_the_stream() {
        let plain = ragged_fastq(20);
        let half = plain.len() / 2;
        let boundary = plain[half..].iter().position(|&byte| byte == b'@').unwrap() + half;
        let mut compressed = bgzf_block(&plain[..boundary]);
        compressed.extend_from_slice(&bgzf_block(b""));
        compressed.extend_from_slice(&bgzf_block(&plain[boundary..]));

        assert_eq!(bgzf_decoded_from(&compressed, 0), plain);
    }

    /// A file without the end-of-file marker still reads to completion.
    #[cfg(feature = "gzip")]
    #[test]
    fn a_missing_bgzf_end_marker_still_reads() {
        let plain = ragged_fastq(10);
        assert_eq!(bgzf_decoded_from(&bgzf_block(&plain), 0), plain);
    }

    /// The table names every block a walk from offset zero finds, and names it the same way.
    #[cfg(feature = "gzip")]
    #[test]
    fn the_bgzf_table_matches_a_walk_from_offset_zero() {
        let plain = ragged_fastq(120);
        let compressed = bgzf_bytes(&plain, 256);
        let access = BgzfAccess::new(&compressed[..]).unwrap();

        let mut offset = 0usize;
        let mut walked = 0usize;
        while let Some(size) = bgzf_block_size(&compressed[offset..]) {
            let entry = &access.table()[walked];
            assert_eq!(
                entry.compressed_offset as usize, offset,
                "block {walked} offset"
            );
            assert_eq!(entry.compressed_bytes as usize, size, "block {walked} size");
            assert_eq!(
                entry.decoded_bytes as usize,
                plain.len().min((walked + 1) * 256) - plain.len().min(walked * 256),
                "block {walked} decoded"
            );
            offset += size;
            walked += 1;
        }
        assert_eq!(walked, access.block_count());
        assert_eq!(offset, compressed.len(), "the walk must reach the end");
    }

    /// The contract, over a ragged fixture and over one whose records tile its blocks exactly.
    #[cfg(feature = "gzip")]
    #[test]
    fn bgzf_conforms() {
        for (plain, payload_bytes, label) in [
            (ragged_fastq(300), 256usize, "ragged"),
            (aligned_fastq(300), 256, "aligned"),
        ] {
            let access = BgzfAccess::new(bgzf_bytes(&plain, payload_bytes)).unwrap();
            conforms(&access, &plain, label);
            an_empty_run_reads_nothing(&access, label);
        }
    }

    /// Plain gzip is not block-addressable, and saying so is how the caller picks a tier.
    #[cfg(feature = "gzip")]
    #[test]
    fn plain_gzip_is_not_block_addressable() {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, &ragged_fastq(10)).unwrap();
        let compressed = encoder.finish().unwrap();
        assert!(BgzfAccess::new(&compressed[..]).is_err());
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn sequence_data_is_not_bgzf() {
        assert!(BgzfAccess::new(&b">seq1\nACGT\n"[..]).is_err());
        assert!(BgzfAccess::new(&b""[..]).is_err());
    }

    /// Rejection hands the container back, so a caller holding a mapping can read it as a stream
    /// without mapping the file again.
    #[cfg(feature = "gzip")]
    #[test]
    fn rejected_bytes_come_back_to_the_caller() {
        let plain = ragged_fastq(4);
        match BgzfAccess::new(plain.clone()) {
            Ok(_) => panic!("plain FASTQ must not read as BGZF"),
            Err(returned) => assert_eq!(returned, plain),
        }
    }

    /// A borrowed view reads exactly what the owning container reads.
    #[cfg(feature = "gzip")]
    #[test]
    fn a_borrowed_bgzf_view_inflates_the_same_bytes() {
        let plain = ragged_fastq(50);
        let owned = BgzfAccess::new(bgzf_bytes(&plain, 256)).unwrap();
        let view = owned.borrowed();
        assert_eq!(view.block_count(), owned.block_count());

        let mut through_view = Vec::new();
        view.reader(0).read_to_end(&mut through_view).unwrap();
        assert_eq!(through_view, plain);
    }

    // endregion: BGZF

    // region: xz

    /// One xz stream over `payload`, in blocks of `block_bytes` when one is given.
    #[cfg(feature = "xz")]
    fn xz_stream(
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

    /// `plain` as one single-block stream per chunk, which is what `cat a.xz b.xz` writes and
    /// the cheapest multi-block fixture there is.
    #[cfg(feature = "xz")]
    fn xz_streams(plain: &[u8], chunk_bytes: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for chunk in plain.chunks(chunk_bytes) {
            out.extend_from_slice(&xz_stream(chunk, None, lzma_rust2::CheckType::Crc64));
        }
        out
    }

    /// The contract, over concatenated streams that are ragged and aligned in turn.
    #[cfg(feature = "xz")]
    #[test]
    fn xz_conforms() {
        for (plain, chunk_bytes, label) in [
            (ragged_fastq(300), 1024usize, "ragged"),
            (aligned_fastq(300), 1024, "aligned"),
        ] {
            let access = XzAccess::new(xz_streams(&plain, chunk_bytes)).unwrap();
            conforms(&access, &plain, label);
            an_empty_run_reads_nothing(&access, label);
        }
    }

    /// Two streams back to back are one file, and a reader that parses only the last footer
    /// silently loses everything the first one holds.
    #[cfg(feature = "xz")]
    #[test]
    fn xz_concatenated_streams_are_all_reachable() {
        let plain = ragged_fastq(120);
        let streams = 4usize;
        let chunk_bytes = plain.len().div_ceil(streams);
        let access = XzAccess::new(xz_streams(&plain, chunk_bytes)).unwrap();
        assert_eq!(access.block_count(), streams);
        conforms(&access, &plain, "four streams");
    }

    /// Stream padding between two streams is legal and must not be read as a footer.
    #[cfg(feature = "xz")]
    #[test]
    fn xz_stream_padding_is_stepped_over() {
        let plain = ragged_fastq(60);
        let half = plain.len() / 2;
        let mut compressed = xz_stream(&plain[..half], None, lzma_rust2::CheckType::Crc64);
        compressed.extend_from_slice(&[0, 0, 0, 0, 0, 0, 0, 0]);
        compressed.extend_from_slice(&xz_stream(
            &plain[half..],
            None,
            lzma_rust2::CheckType::Crc64,
        ));
        compressed.extend_from_slice(&[0, 0, 0, 0]);

        let access = XzAccess::new(compressed).unwrap();
        assert_eq!(access.block_count(), 2);
        conforms(&access, &plain, "padded");
    }

    /// The check trails every payload, so its width has to come from the stream flags rather
    /// than from a constant.
    #[cfg(feature = "xz")]
    #[test]
    fn every_xz_check_type_locates_its_payloads() {
        let plain = ragged_fastq(80);
        let chunk_bytes = plain.len().div_ceil(3);
        for check in [
            lzma_rust2::CheckType::None,
            lzma_rust2::CheckType::Crc32,
            lzma_rust2::CheckType::Crc64,
            lzma_rust2::CheckType::Sha256,
        ] {
            let mut compressed = Vec::new();
            for chunk in plain.chunks(chunk_bytes) {
                compressed.extend_from_slice(&xz_stream(chunk, None, check));
            }
            let access = XzAccess::new(compressed).unwrap();
            conforms(&access, &plain, &format!("{check:?}"));
        }
    }

    /// One stream cut into blocks by the writer, which is the layout an explicit block size
    /// gives and the one where index records rather than stream headers locate the blocks.
    #[cfg(feature = "xz")]
    #[test]
    fn xz_blocks_of_one_stream_conform() {
        // The writer clamps a block to at least the dictionary, which preset zero puts at
        // 256 KiB, so the fixture has to be larger than that to hold more than one block.
        let block_bytes = 256usize << 10;
        let plain = aligned_fastq(3 * block_bytes / 32);
        let compressed = xz_stream(
            &plain,
            Some(block_bytes as u64),
            lzma_rust2::CheckType::Crc64,
        );
        let access = XzAccess::new(compressed).unwrap();
        assert_eq!(access.block_count(), 3);
        conforms(&access, &plain, "one stream of three blocks");
    }

    /// A file of one block reads faster as a stream, so it is handed back rather than wrapped.
    #[cfg(feature = "xz")]
    #[test]
    fn a_single_block_xz_file_is_not_block_addressable() {
        let compressed = xz_stream(&ragged_fastq(40), None, lzma_rust2::CheckType::Crc64);
        assert_eq!(Codec::detect(&compressed), Codec::Xz);
        assert!(XzAccess::new(&compressed[..]).is_err());
    }

    #[cfg(feature = "xz")]
    #[test]
    fn sequence_data_is_not_xz() {
        assert!(XzAccess::new(&b">seq1\nACGT\n"[..]).is_err());
        assert!(XzAccess::new(&b""[..]).is_err());
    }

    /// A threaded writer states both optional sizes, so the property byte sits eight bytes
    /// further into the header than a single-threaded writer puts it. These are the twenty bytes
    /// of one such header, and a reader assuming a fixed offset reads `0x03` as the dictionary
    /// code and gets every payload bound after it wrong.
    #[cfg(feature = "xz")]
    #[test]
    fn a_threaded_xz_block_header_moves_the_property_byte() {
        let header = [
            0x04, 0xc0, 0xf0, 0xc1, 0x9b, 0x03, 0x80, 0x80, 0x80, 0x0c, 0x21, 0x01, 0x16, 0x00,
            0x00, 0x00, 0x19, 0xec, 0x1d, 0x01,
        ];
        let unpadded_size = header.len() + 64 + 8;
        let mut compressed = header.to_vec();
        compressed.resize(unpadded_size, 0);

        let (payload_offset, payload_end, dictionary_size) =
            xz_payload_bounds(&compressed, 0, unpadded_size, 8).unwrap();
        assert_eq!(payload_offset, 20, "the header is twenty bytes wide");
        assert_eq!(
            payload_end,
            unpadded_size - 8,
            "the check trails the payload"
        );
        assert_eq!(dictionary_size, 8 << 20, "property 0x16 is eight mebibytes");
    }

    /// A chain of two filters is a file this reader will not guess at, and rejecting it is how
    /// the caller falls back to the stream tier instead of decoding rubbish.
    #[cfg(feature = "xz")]
    #[test]
    fn an_xz_filter_chain_is_rejected() {
        use io::Write;
        let plain = ragged_fastq(60);
        let half = plain.len() / 2;
        let mut compressed = Vec::new();
        for chunk in [&plain[..half], &plain[half..]] {
            let mut options = lzma_rust2::XzOptions::with_preset(0);
            options.prepend_pre_filter(lzma_rust2::FilterType::Delta, 1);
            let mut writer = lzma_rust2::XzWriter::new(Vec::new(), options).unwrap();
            writer.write_all(chunk).unwrap();
            compressed.extend_from_slice(&writer.finish().unwrap());
        }
        assert!(XzAccess::new(&compressed[..]).is_err());
    }

    // endregion: xz

    // region: lz4

    /// `plain` as one frame written under `frame_info`.
    #[cfg(feature = "lz4")]
    fn lz4_frame(plain: &[u8], frame_info: lz4_flex::frame::FrameInfo) -> Vec<u8> {
        use io::Write;
        let mut encoder = lz4_flex::frame::FrameEncoder::with_frame_info(frame_info, Vec::new());
        encoder.write_all(plain).unwrap();
        encoder.finish().unwrap()
    }

    /// A frame descriptor with block independence set and the smallest frame maximum.
    #[cfg(feature = "lz4")]
    fn lz4_independent_blocks() -> lz4_flex::frame::FrameInfo {
        lz4_flex::frame::FrameInfo::new().block_size(lz4_flex::frame::BlockSize::Max64KB)
    }

    /// A skippable frame carrying `bytes` of payload nothing reads.
    #[cfg(feature = "lz4")]
    fn lz4_skippable_frame(bytes: usize) -> Vec<u8> {
        let mut frame = Vec::new();
        frame.extend_from_slice(&0x184D_2A50u32.to_le_bytes());
        frame.extend_from_slice(&(bytes as u32).to_le_bytes());
        frame.extend(std::iter::repeat_n(0xA5, bytes));
        frame
    }

    /// The contract, over a ragged fixture and over one whose records tile its blocks exactly.
    ///
    /// Each fixture is three frame maxima long, so the frame holds several blocks.
    #[cfg(feature = "lz4")]
    #[test]
    fn lz4_conforms() {
        let records = 3 * (64 << 10) / 32;
        for (plain, label) in [
            (ragged_fastq(records), "ragged"),
            (aligned_fastq(records), "aligned"),
        ] {
            let access = Lz4Access::new(lz4_frame(&plain, lz4_independent_blocks())).unwrap();
            conforms(&access, &plain, label);
            an_empty_run_reads_nothing(&access, label);
        }
    }

    /// A writer that flushes as it goes emits blocks far under the frame maximum, so the maximum
    /// is a ceiling and never a size. This is the pattern a streaming logger writes.
    #[cfg(feature = "lz4")]
    #[test]
    fn lz4_blocks_are_at_most_the_frame_maximum() {
        use io::Write;
        let plain = ragged_fastq(40);
        let mut encoder =
            lz4_flex::frame::FrameEncoder::with_frame_info(lz4_independent_blocks(), Vec::new());
        let mut boundary = 0usize;
        for chunk_bytes in [1usize, 6, 93, 400] {
            let chunk_end = (boundary + chunk_bytes).min(plain.len());
            encoder.write_all(&plain[boundary..chunk_end]).unwrap();
            encoder.flush().unwrap();
            boundary = chunk_end;
        }
        encoder.write_all(&plain[boundary..]).unwrap();
        let access = Lz4Access::new(encoder.finish().unwrap()).unwrap();

        assert!(access.block_count() >= 5, "one block per flush, at least");
        // A one-byte block is stored rather than compressed, since compressing it would grow
        // it, and a stored block is the one case where the count on disk is the decoded count.
        assert_eq!(access.bytes_in_block(0), BlockBytes::Exact(1));
        let ceilings: Vec<BlockBytes> = (0..access.block_count())
            .map(|index| access.bytes_in_block(index))
            .filter(|stated| stated.exact().is_none())
            .collect();
        assert!(
            ceilings.contains(&BlockBytes::AtMost(64 << 10)),
            "a compressed block states no size of its own, only its frame's ceiling"
        );
        conforms(&access, &plain, "one block per flush");
    }

    /// Frames back to back are one file, each with its own descriptor, so the frame maximum is a
    /// property of a block rather than of the file.
    #[cfg(feature = "lz4")]
    #[test]
    fn lz4_concatenated_frames_are_all_reachable() {
        let plain = ragged_fastq(200);
        let half = plain.len() / 2;
        let mut compressed = lz4_frame(&plain[..half], lz4_independent_blocks());
        let wider = lz4_flex::frame::FrameInfo::new()
            .block_size(lz4_flex::frame::BlockSize::Max256KB)
            .content_size(Some((plain.len() - half) as u64));
        compressed.extend_from_slice(&lz4_frame(&plain[half..], wider));

        let access = Lz4Access::new(compressed).unwrap();
        assert_eq!(access.block_count(), 2, "one block per frame here");
        assert_eq!(access.bytes_in_block(0), BlockBytes::AtMost(64 << 10));
        assert_eq!(access.bytes_in_block(1), BlockBytes::AtMost(256 << 10));
        conforms(&access, &plain, "two frames");
    }

    /// A skippable frame is stepped over wherever it sits, and a walk reading one as a
    /// descriptor loses every block after it.
    #[cfg(feature = "lz4")]
    #[test]
    fn an_lz4_skippable_frame_is_stepped_over() {
        let plain = ragged_fastq(200);
        let half = plain.len() / 2;
        let mut compressed = lz4_skippable_frame(33);
        compressed.extend_from_slice(&lz4_frame(&plain[..half], lz4_independent_blocks()));
        compressed.extend_from_slice(&lz4_skippable_frame(7));
        compressed.extend_from_slice(&lz4_frame(&plain[half..], lz4_independent_blocks()));

        let access = Lz4Access::new(compressed).unwrap();
        assert_eq!(access.block_count(), 2);
        conforms(&access, &plain, "skippable frames");
    }

    /// Linked blocks reference the window before them, so no block decodes alone and the file
    /// belongs on the stream tier however valid it is.
    #[cfg(feature = "lz4")]
    #[test]
    fn an_lz4_frame_of_linked_blocks_is_not_block_addressable() {
        let plain = ragged_fastq(200);
        let linked = lz4_flex::frame::FrameInfo::new()
            .block_size(lz4_flex::frame::BlockSize::Max64KB)
            .block_mode(lz4_flex::frame::BlockMode::Linked);
        let compressed = lz4_frame(&plain, linked);
        assert_eq!(Codec::detect(&compressed), Codec::Lz4);
        assert!(Lz4Access::new(&compressed[..]).is_err());
    }

    /// A legacy frame carries no descriptor and no end mark, so it is detected and refused.
    #[cfg(feature = "lz4")]
    #[test]
    fn a_legacy_lz4_frame_is_not_block_addressable() {
        let mut compressed = 0x184C_2102u32.to_le_bytes().to_vec();
        compressed.extend_from_slice(&64u32.to_le_bytes());
        compressed.extend(std::iter::repeat_n(0u8, 64));
        assert!(Lz4Access::new(&compressed[..]).is_err());
    }

    #[cfg(feature = "lz4")]
    #[test]
    fn sequence_data_is_not_lz4() {
        assert!(Lz4Access::new(&b">seq1\nACGT\n"[..]).is_err());
        assert!(Lz4Access::new(&b""[..]).is_err());
    }

    /// lz4 carries no per-block check, so a payload cut short decodes to fewer bytes rather
    /// than to an error whenever the cut lands between two tokens. A decoder that runs out of
    /// input mid-token does say so, which catches the other cuts and nothing more.
    ///
    /// The frame here is doctored to keep its chain intact — the prefix states the shortened
    /// length — because a chain left inconsistent is caught by the walk instead.
    #[cfg(feature = "lz4")]
    #[test]
    fn a_truncated_lz4_payload_can_decode_without_an_error() {
        let plain = aligned_fastq(2000);
        let whole = lz4_frame(
            &plain,
            lz4_flex::frame::FrameInfo::new()
                .block_size(lz4_flex::frame::BlockSize::Max64KB)
                .content_checksum(false),
        );
        let stored_size = u32::from_le_bytes(whole[7..11].try_into().unwrap()) as usize;
        assert_eq!(
            stored_size & 0x8000_0000,
            0,
            "the fixture must be one compressed block"
        );

        let cut_short = |cut: usize| {
            let mut frame = whole[..7].to_vec();
            frame.extend_from_slice(&(cut as u32).to_le_bytes());
            frame.extend_from_slice(&whole[11..11 + cut]);
            frame.extend_from_slice(&0u32.to_le_bytes());
            frame
        };

        let mut decoded = Vec::new();
        let quiet = (stored_size / 2..stored_size).find(|&cut| {
            let access = Lz4Access::new(cut_short(cut)).unwrap();
            let quiet = access.decoder().decode_into(0, &mut decoded).is_ok();
            quiet
        });
        let cut = quiet.expect("some cut lands between two tokens and decodes quietly");
        assert!(cut < stored_size, "the payload really is short");
        assert!(decoded.len() < plain.len(), "and it really did lose bytes");
        assert_eq!(decoded, plain[..decoded.len()], "what decoded is a prefix");
    }

    // endregion: lz4

    // region: Choosing a tier

    /// Detection names the container and this names the tier, so a container with no block
    /// layout hands its bytes straight back.
    #[test]
    fn a_container_without_blocks_is_not_block_addressable() {
        for codec in [Codec::Plain, Codec::Gzip, Codec::Zstd, Codec::Bzip2] {
            assert!(
                block_access(b">seq1\nACGT\n".to_vec(), codec).is_err(),
                "{codec:?}"
            );
        }
    }

    /// Every container that does have one is reachable through the same entry point.
    #[test]
    fn every_block_container_is_reachable_through_one_entry_point() {
        #[cfg(feature = "gzip")]
        {
            let plain = ragged_fastq(120);
            let access = block_access(bgzf_bytes(&plain, 256), Codec::Bgzf).unwrap();
            conforms(access.as_ref(), &plain, "BGZF through block_access");
        }
        #[cfg(feature = "xz")]
        {
            let plain = ragged_fastq(120);
            let access = block_access(xz_streams(&plain, 1024), Codec::Xz).unwrap();
            conforms(access.as_ref(), &plain, "xz through block_access");
        }
        #[cfg(feature = "lz4")]
        {
            let plain = ragged_fastq(3 * (64 << 10) / 32);
            let access =
                block_access(lz4_frame(&plain, lz4_independent_blocks()), Codec::Lz4).unwrap();
            conforms(access.as_ref(), &plain, "lz4 through block_access");
        }
    }

    // endregion: Choosing a tier
}
