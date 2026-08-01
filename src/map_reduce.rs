//! Scheduling: run one input through a per-record visitor.
//!
//! Never parses — it picks byte ranges and hands them to [`RandomAccess`] or
//! [`ForwardAccess`]. Deduplication schedules its own passes and depends on nothing here.

use std::io;
use std::ops::Range;

use crate::{map_file, next_record_boundary, ForwardAccess, RandomAccess, Record, SequenceFormat};

/// Below twice this, splitting costs more than it saves.
pub const MINIMUM_SHARD_BYTES: usize = 1 << 20;

/// Fill `ranges` with shards that each begin at a record boundary and tile the input in order.
///
/// Yields fewer than `workers` when the input is small or no further boundary is found.
pub fn split_into_shard_ranges(
    bytes: &[u8],
    format: SequenceFormat,
    workers: usize,
    ranges: &mut Vec<Range<usize>>,
) {
    ranges.clear();
    if workers <= 1 || bytes.len() < MINIMUM_SHARD_BYTES.saturating_mul(2) {
        ranges.push(0..bytes.len());
        return;
    }

    let target = (bytes.len() / workers).max(MINIMUM_SHARD_BYTES);
    let mut start = 0usize;
    let mut cursor = target;
    while cursor < bytes.len() {
        match next_record_boundary(bytes, cursor, format) {
            // A boundary at or before `start` would make an empty or reversed shard.
            Some(next) if next > start => {
                ranges.push(start..next);
                start = next;
                cursor = next + target;
            }
            _ => break,
        }
    }
    ranges.push(start..bytes.len());
}

/// Run every input through `visit`, in argument order, folding into `state`.
///
/// Inputs holding no records contribute nothing and are not an error. A record carries its
/// own format, so nothing has to be detected here or handed to the tool ahead of time; a
/// tool that mirrors its input calls [`crate::RecordWriter::adopt`], which also catches an
/// input whose format disagrees with an earlier one.
///
/// A named file is memory-mapped and read in place; `-` is read through a bounded window and
/// may appear at most once, because there is only one standard input to consume.
pub fn for_each_record_at_paths<T>(
    paths: &[String],
    state: &mut T,
    mut visit: impl FnMut(&mut T, Record<'_>) -> io::Result<()>,
) -> io::Result<()> {
    if paths.iter().filter(|path| path.as_str() == "-").count() > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "standard input can only be read once, but '-' was given more than once",
        ));
    }

    for path in paths {
        if path == "-" {
            ForwardAccess::from_stdin().for_each_record(|record| visit(state, record))?;
        } else {
            let mapping = map_file(path)?;
            RandomAccess::new(&mapping).for_each_record(|record| visit(state, record))?;
        }
    }
    Ok(())
}

/// [`for_each_record_at_paths`] over a single buffer already in memory.
pub fn for_each_record_in_bytes<T>(
    bytes: &[u8],
    state: &mut T,
    mut visit: impl FnMut(&mut T, Record<'_>) -> io::Result<()>,
) -> io::Result<()> {
    RandomAccess::new(bytes).for_each_record(|record| visit(state, record))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fasta(records: usize, width: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            data.extend_from_slice(format!(">contig{index} sample\n").as_bytes());
            data.extend(std::iter::repeat(b"ACGT"[index % 4]).take(width));
            data.push(b'\n');
        }
        data
    }

    fn fastq(records: usize, width: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for index in 0..records {
            data.extend_from_slice(format!("@read{index} sample\n").as_bytes());
            data.extend(std::iter::repeat(b"ACGT"[index % 4]).take(width));
            data.extend_from_slice(b"\n+\n");
            data.extend(std::iter::repeat(b'I').take(width));
            data.push(b'\n');
        }
        data
    }

    fn shards(bytes: &[u8], format: SequenceFormat, workers: usize) -> Vec<Range<usize>> {
        let mut ranges = Vec::new();
        split_into_shard_ranges(bytes, format, workers, &mut ranges);
        ranges
    }

    fn headers_of(bytes: &[u8], range: Range<usize>) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        RandomAccess::new(bytes)
            .for_each_record_in(range, |record| {
                out.push(record.header.to_vec());
                Ok(())
            })
            .unwrap();
        out
    }

    #[test]
    fn small_inputs_stay_whole() {
        let data = fasta(4, 20);
        assert_eq!(shards(&data, SequenceFormat::Fasta, 8), vec![0..data.len()]);
    }

    #[test]
    fn one_worker_stays_whole() {
        let data = fasta(200_000, 60);
        assert_eq!(shards(&data, SequenceFormat::Fasta, 1), vec![0..data.len()]);
    }

    #[test]
    fn shards_are_contiguous_and_ordered() {
        let data = fasta(200_000, 60);
        let ranges = shards(&data, SequenceFormat::Fasta, 8);
        assert!(ranges.len() > 1, "expected a split, got {}", ranges.len());
        assert_eq!(ranges[0].start, 0);
        assert_eq!(ranges.last().unwrap().end, data.len());
        for pair in ranges.windows(2) {
            assert_eq!(pair[0].end, pair[1].start, "shards must tile the input");
        }
    }

    /// Every shard must begin exactly where a record begins, or a worker parses garbage.
    #[test]
    fn shards_begin_on_record_boundaries() {
        for (data, format, sigil) in [
            (fasta(200_000, 60), SequenceFormat::Fasta, b'>'),
            (fastq(200_000, 60), SequenceFormat::Fastq, b'@'),
        ] {
            for range in shards(&data, format, 8) {
                assert_eq!(data[range.start], sigil, "shard did not start at a record");
            }
        }
    }

    /// The union of the shards must be exactly the records of the whole input, in order.
    #[test]
    fn sharding_preserves_every_record_in_order() {
        for (data, format) in [
            (fasta(200_000, 60), SequenceFormat::Fasta),
            (fastq(200_000, 60), SequenceFormat::Fastq),
        ] {
            let whole = headers_of(&data, 0..data.len());
            let mut sharded = Vec::new();
            for range in shards(&data, format, 8) {
                sharded.extend(headers_of(&data, range));
            }
            assert_eq!(sharded, whole);
        }
    }

    #[test]
    fn reused_range_buffer_does_not_accumulate() {
        let data = fasta(200_000, 60);
        let mut ranges = Vec::new();
        split_into_shard_ranges(&data, SequenceFormat::Fasta, 8, &mut ranges);
        let first = ranges.len();
        split_into_shard_ranges(&data, SequenceFormat::Fasta, 8, &mut ranges);
        assert_eq!(ranges.len(), first);
    }

    /// A FASTQ quality line may open with '@', so a naive boundary search lands mid-record.
    #[test]
    fn next_record_boundary_rejects_fastq_quality_opening_with_at() {
        // Quality '@' is Phred 31, entirely legal, and here it opens the fourth line.
        let data = b"@read0\nACGT\n+\n@III\n@read1\nTGCA\n+\nIIII\n".to_vec();
        let boundary = next_record_boundary(&data, 1, SequenceFormat::Fastq).unwrap();
        assert_eq!(&data[boundary..boundary + 6], b"@read1");
    }

    #[test]
    fn next_record_boundary_finds_the_next_fasta_record() {
        let data = b">a\nACGT\n>b\nTGCA\n";
        assert_eq!(next_record_boundary(data, 0, SequenceFormat::Fasta), Some(0));
        assert_eq!(next_record_boundary(data, 1, SequenceFormat::Fasta), Some(8));
    }

    #[test]
    fn next_record_boundary_returns_none_past_the_end() {
        let data = b">a\nACGT\n";
        assert_eq!(next_record_boundary(data, 3, SequenceFormat::Fasta), None);
    }
}
