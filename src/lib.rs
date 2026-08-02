//! Faster FASTA — command-line utilities over memory-mapped and streamed sequence files.
//!
//! Five modules, in dependency order:
//!
//! - [`codecs`] — which compression container a prefix of bytes opens, and how to decode it.
//! - [`blocks`] — the block-addressable containers, as block tables over compressed bytes.
//! - [`records`] — what a record is, plus format detection, Phred conversions, and the parser.
//! - [`files`] — the three access tiers, record rendering, and the epilogue every `main` shares.
//! - [`scheduling`] — choosing what each worker reads, and how many run at once.
//!
//! The root declares these and re-exports nothing, so every item has exactly one path and a
//! tool's `use` lines say which module owns what it imports.

// Every public item is documented today, so the cost of requiring it is zero now and rises
// every week it is deferred.
#![deny(missing_docs)]

pub mod codecs;

pub mod blocks;

pub mod files;
pub mod records;
pub mod scheduling;
