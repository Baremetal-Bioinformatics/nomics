//! HaplotypeCaller's read stream (spec §1): BED intervals, the padded
//! overlap query, the standard read filters in order, and the positional
//! reservoir downsampler drawing from GATK's shared `java.util.Random`.
//!
//! Gates: the per-filter counts GATK logs, and (later) the RNG draw count
//! implied by the QD-jitter stream.

use std::io;
use std::path::Path;

use nomics_io::bam::{RecordSet, RecordView};

use crate::jrandom::JavaRandom;

/// `--assembly-region-padding`.
pub const REGION_PADDING: i64 = 100;
/// `--max-reads-per-alignment-start`.
pub const MAX_READS_PER_START: usize = 50;
/// `--minimum-mapping-quality`.
pub const MIN_MAPQ: u8 = 20;

/// 1-based closed interval on one contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: i64,
    pub end: i64,
}

/// Merge sorted intervals that overlap or abut (`IntervalMergingRule.ALL`).
fn merge(mut v: Vec<Interval>) -> Vec<Interval> {
    v.sort_by_key(|i| (i.start, i.end));
    let mut out: Vec<Interval> = Vec::with_capacity(v.len());
    for i in v {
        match out.last_mut() {
            Some(l) if i.start - 1 <= l.end => l.end = l.end.max(i.end),
            _ => out.push(i),
        }
    }
    out
}

/// Per-contig (header order) unpadded and 100-padded intervals from a BED.
pub struct Intervals {
    pub unpadded: Vec<Vec<Interval>>,
    pub padded: Vec<Vec<Interval>>,
}

impl Intervals {
    /// Load a BED (`StartOffset.ONE`: start + 1, end), keyed by the BAM's
    /// contig order and clipped to contig lengths when padding.
    pub fn from_bed(path: &Path, contigs: &[(Vec<u8>, i64)]) -> io::Result<Self> {
        let index: std::collections::HashMap<&[u8], usize> = contigs
            .iter()
            .enumerate()
            .map(|(i, (n, _))| (n.as_slice(), i))
            .collect();
        let mut raw: Vec<Vec<Interval>> = vec![Vec::new(); contigs.len()];
        for line in std::fs::read_to_string(path)?.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let c = *index.get(f[0].as_bytes()).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("BED contig {} not in BAM header", f[0]),
                )
            })?;
            let parse = |s: &str| {
                s.trim()
                    .parse::<i64>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            };
            raw[c].push(Interval {
                start: parse(f[1])? + 1,
                end: parse(f[2])?,
            });
        }
        let unpadded: Vec<Vec<Interval>> = raw.into_iter().map(merge).collect();
        let padded = unpadded
            .iter()
            .zip(contigs)
            .map(|(v, (_, len))| {
                merge(
                    v.iter()
                        .map(|i| Interval {
                            start: (i.start - REGION_PADDING).max(1),
                            end: (i.end + REGION_PADDING).min(*len),
                        })
                        .collect(),
                )
            })
            .collect();
        Ok(Self { unpadded, padded })
    }
}

/// Filters in `makeStandardHCReadFilters` order; a read is counted by the
/// first filter it fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    MappingQuality,
    MappingQualityAvailable,
    Mapped,
    NotSecondary,
    NotDuplicate,
    PassesVendorQc,
    NonZeroReferenceLength,
    GoodCigar,
    Wellformed,
}

pub const FILTERS: [Filter; 9] = [
    Filter::MappingQuality,
    Filter::MappingQualityAvailable,
    Filter::Mapped,
    Filter::NotSecondary,
    Filter::NotDuplicate,
    Filter::PassesVendorQc,
    Filter::NonZeroReferenceLength,
    Filter::GoodCigar,
    Filter::Wellformed,
];

impl Filter {
    /// GATK's class name, as in the log's filter summary.
    pub fn name(self) -> &'static str {
        match self {
            Filter::MappingQuality => "MappingQualityReadFilter",
            Filter::MappingQualityAvailable => "MappingQualityAvailableReadFilter",
            Filter::Mapped => "MappedReadFilter",
            Filter::NotSecondary => "NotSecondaryAlignmentReadFilter",
            Filter::NotDuplicate => "NotDuplicateReadFilter",
            Filter::PassesVendorQc => "PassesVendorQualityCheckReadFilter",
            Filter::NonZeroReferenceLength => "NonZeroReferenceLengthAlignmentReadFilter",
            Filter::GoodCigar => "GoodCigarReadFilter",
            Filter::Wellformed => "WellformedReadFilter",
        }
    }
}

/// `CigarUtils.isGood`: htsjdk-valid (non-zero lengths, clipping only at the
/// ends), no adjacent indels, no deletion as the first or last non-clip op.
fn good_cigar(cigar: &[u32]) -> bool {
    let n = cigar.len();
    if n == 0 {
        return false;
    }
    let op = |i: usize| cigar[i] & 0xf;
    if cigar.iter().any(|&c| c >> 4 == 0) {
        return false;
    }
    // H only at the ends; S only at the ends or next to an end H.
    for i in 0..n {
        match op(i) {
            5 if i != 0 && i != n - 1 => return false,
            4 => {
                let left_ok = i == 0 || (i == 1 && op(0) == 5);
                let right_ok = i == n - 1 || (i == n - 2 && op(n - 1) == 5);
                if !left_ok && !right_ok {
                    return false;
                }
            }
            _ => {}
        }
    }
    let is_indel = |o: u32| o == 1 || o == 2;
    if (1..n).any(|i| is_indel(op(i)) && is_indel(op(i - 1))) {
        return false;
    }
    let core: Vec<u32> = (0..n).map(op).filter(|&o| o != 4 && o != 5).collect();
    !(core.first() == Some(&2) || core.last() == Some(&2))
}

/// `WellformedReadFilter` for mapped reads.
fn wellformed(r: &RecordView, contig_len: i64) -> bool {
    let cigar: Vec<u32> = r.cigar().collect();
    let start = r.pos() as i64 + 1;
    let end = start + r.reference_length() as i64 - 1;
    let read_len: u32 = cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 1 | 4 | 7 | 8))
        .map(|&c| c >> 4)
        .sum();
    start > 0
        && end <= contig_len
        && r.string_tag(b"RG").is_some()
        && r.quals().is_some_and(|q| q.len() == r.len())
        && r.len() as u32 == read_len
        && r.len() > 0
        && !cigar.iter().any(|&c| c & 0xf == 3)
}

/// Reads sharing one alignment start: one `ReservoirDownsampler` batch.
#[derive(Debug, Clone)]
pub struct Group {
    /// 1-based alignment start.
    pub pos: i64,
    /// Filtered reads at `pos`, in stream (file) order.
    pub reads: Vec<u32>,
}

/// The filtered read stream of a whole run, before downsampling.
pub struct ReadStream {
    /// Per contig (header order): start-position groups in stream order.
    pub groups: Vec<Vec<Group>>,
    /// Reads returned by the overlap query.
    pub processed: u64,
    /// Rejections per filter, `FILTERS` order.
    pub filtered: [u64; 9],
}

/// `ReservoirDownsampler` over one start position: the first 50 reads fill
/// the reservoir; read `t > 50` replaces slot `nextInt(t)` if it is < 50.
/// Returns the kept reads in slot order and the number of draws.
pub fn downsample(group: &Group, rng: &mut JavaRandom) -> (Vec<u32>, u64) {
    let mut kept: Vec<u32> = group
        .reads
        .iter()
        .take(MAX_READS_PER_START)
        .copied()
        .collect();
    let mut draws = 0;
    for (i, &r) in group.reads.iter().enumerate().skip(MAX_READS_PER_START) {
        let slot = rng.next_int(i as i32 + 1);
        draws += 1;
        if (slot as usize) < MAX_READS_PER_START {
            kept[slot as usize] = r;
        }
    }
    (kept, draws)
}

/// Collects, in one forward pass over a coordinate-sorted BAM on disk, the
/// records each interval group's overlap query would return (the same rule
/// as [`read_stream`]), so a group can be called on its records alone
/// instead of the whole file being held in memory. A group's records are
/// complete once the stream passes its last padded interval; groups that
/// have not been taken yet are filled along the way, so groups in genome
/// order hold about one group's records at a time.
pub struct GroupCollector<'a> {
    cursor: nomics_io::bamstream::RecordCursor,
    groups: &'a [Intervals],
    /// Each group's last padded interval end, as (contig, end).
    last: Vec<Option<(usize, i64)>>,
    /// Per group, per contig: first padded interval that can still overlap.
    at: Vec<Vec<usize>>,
    /// Per group: `block_size`-prefixed records; `None` once taken.
    buffers: Vec<Option<Vec<u8>>>,
    /// (contig, start) of the last positioned record read; contig
    /// `usize::MAX` once the unplaced records are reached.
    pos: (usize, i64),
    finished: bool,
    record: Vec<u8>,
}

impl<'a> GroupCollector<'a> {
    pub fn new(
        stream: nomics_io::bamstream::BamStream,
        groups: &'a [Intervals],
    ) -> GroupCollector<'a> {
        let last = groups
            .iter()
            .map(|g| {
                g.padded
                    .iter()
                    .enumerate()
                    .rev()
                    .find_map(|(c, v)| v.last().map(|i| (c, i.end)))
            })
            .collect();
        GroupCollector {
            cursor: nomics_io::bamstream::RecordCursor::new(stream),
            groups,
            last,
            at: groups.iter().map(|g| vec![0; g.padded.len()]).collect(),
            buffers: groups.iter().map(|_| Some(Vec::new())).collect(),
            pos: (0, 0),
            finished: false,
            record: Vec::new(),
        }
    }

    pub fn header(&self) -> &nomics_io::bam::Header {
        self.cursor.header()
    }

    fn complete(&self, g: usize) -> bool {
        match self.last[g] {
            None => true,
            Some(last) => self.finished || self.pos > last,
        }
    }

    /// The records of group `g`, in file order. Each group can be taken once.
    pub fn take(&mut self, g: usize) -> io::Result<RecordSet> {
        while !self.complete(g) {
            let Some(rec) = self.cursor.peek()? else {
                self.finished = true;
                break;
            };
            self.record.clear();
            self.record.extend_from_slice(rec);
            self.cursor.advance();
            self.add()?;
        }
        let records = self.buffers[g]
            .take()
            .ok_or_else(|| io::Error::other(format!("interval group {g} taken twice")))?;
        RecordSet::from_parts(self.header().clone(), &records)
    }

    fn add(&mut self) -> io::Result<()> {
        let r = RecordView::new(&self.record)?;
        if r.ref_id() < 0 {
            self.pos = (usize::MAX, 0);
            return Ok(());
        }
        if r.pos() < 0 {
            return Ok(());
        }
        let c = r.ref_id() as usize;
        let start = r.pos() as i64 + 1;
        let end = if r.flag() & 0x4 != 0 {
            start
        } else {
            start + r.reference_length() as i64 - 1
        };
        self.pos = (c, start);
        for (g, group) in self.groups.iter().enumerate() {
            let Some(buffer) = &mut self.buffers[g] else {
                continue;
            };
            let Some(ivs) = group.padded.get(c) else {
                continue;
            };
            let at = &mut self.at[g][c];
            while *at < ivs.len() && ivs[*at].end < start {
                *at += 1;
            }
            if *at < ivs.len() && ivs[*at].start <= end {
                buffer.extend_from_slice(&(self.record.len() as i32).to_le_bytes());
                buffer.extend_from_slice(&self.record);
            }
        }
        Ok(())
    }
}

/// Run the overlap query, filters and downsampler over `bam`, drawing from `rng`.
pub fn read_stream(
    bam: &RecordSet,
    intervals: &Intervals,
    contig_lens: &[i64],
) -> io::Result<ReadStream> {
    let n_contigs = contig_lens.len();
    let mut out = ReadStream {
        groups: vec![Vec::new(); n_contigs],
        processed: 0,
        filtered: [0; 9],
    };
    let mut iv_at = vec![0usize; n_contigs];
    for i in 0..bam.len() {
        let r = RecordView::new(bam.record(i))?;
        let flag = r.flag();
        if r.ref_id() < 0 || r.pos() < 0 {
            continue; // no assigned position: never returned by an interval query
        }
        let c = r.ref_id() as usize;
        let start = r.pos() as i64 + 1;
        // htsjdk's query treats an unmapped read placed at its mate's position
        // as spanning its start base; such reads reach (and fail) the filters.
        let end = if flag & 0x4 != 0 {
            start
        } else {
            start + r.reference_length() as i64 - 1
        };
        // Reads arrive by start, so intervals ending before this read starts
        // cannot overlap it or any later read.
        let ivs = &intervals.padded[c];
        let at = &mut iv_at[c];
        while *at < ivs.len() && ivs[*at].end < start {
            *at += 1;
        }
        if !(*at < ivs.len() && ivs[*at].start <= end) {
            continue;
        }
        out.processed += 1;
        let cigar: Vec<u32> = r.cigar().collect();
        let fail = if r.mapq() < MIN_MAPQ {
            Some(0)
        } else if r.mapq() == 255 {
            Some(1)
        } else if flag & 0x4 != 0 {
            Some(2)
        } else if flag & 0x100 != 0 {
            Some(3)
        } else if flag & 0x400 != 0 {
            Some(4)
        } else if flag & 0x200 != 0 {
            Some(5)
        } else if !cigar
            .iter()
            .any(|&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8) && c >> 4 > 0)
        {
            Some(6)
        } else if !good_cigar(&cigar) {
            Some(7)
        } else if !wellformed(&r, contig_lens[c]) {
            Some(8)
        } else {
            None
        };
        if let Some(f) = fail {
            out.filtered[f] += 1;
            continue;
        }
        let groups = &mut out.groups[c];
        match groups.last_mut() {
            Some(g) if g.pos == start => g.reads.push(i as u32),
            _ => groups.push(Group {
                pos: start,
                reads: vec![i as u32],
            }),
        }
    }
    Ok(out)
}
