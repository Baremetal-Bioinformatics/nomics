//! Coordinate sort and duplicate marking, equivalent to GATK 4.6.2.0
//! `MarkDuplicates` (Picard 3.4.0) as run by sarek. Section numbers (§) refer
//! to `docs/spec/markdup.md`.
//!
//! Gate: `preprocessing/markduplicates/*.md.bam` records and flags.

mod header;
mod metrics;
mod optical;

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::io;
use std::path::PathBuf;

use std::io::Write;

use nomics_io::bam::{BinRule, Header, IntType, Record, RecordSet, RecordView, Value};
use nomics_io::bamstream::{BamStream, RecordCursor};

pub use metrics::DuplicationMetrics;

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

pub struct Options {
    pub inputs: Vec<PathBuf>,
    pub output: PathBuf,
    pub metrics: PathBuf,
    /// `@PG VN` and `CL` of the MarkDuplicates line; header metadata only.
    pub pg_version: String,
    pub pg_command_line: String,
    /// Worker threads for re-encoding and BGZF compression.
    pub threads: usize,
}

// ---------------------------------------------------------------------------
// Input: k-way merge (§1.5)

struct Head {
    rec: Vec<u8>,
    input: usize,
    local: usize,
    original_flag: u16,
}

impl PartialEq for Head {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Head {}
impl PartialOrd for Head {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Head {
    // BinaryHeap is a max-heap: reverse so the smallest record pops first.
    // Full ties need the same read name in two inputs; fastp chunking rules
    // that out, so the input index is only a deterministic last resort.
    fn cmp(&self, other: &Self) -> Ordering {
        coordinate_cmp_with_flags(
            &other.rec,
            other.original_flag,
            &self.rec,
            self.original_flag,
        )
        .then(other.input.cmp(&self.input))
    }
}

fn coordinate_cmp_with_flags(a: &[u8], a_flag: u16, b: &[u8], b_flag: u16) -> Ordering {
    let (va, vb) = (RecordView::new(a).unwrap(), RecordView::new(b).unwrap());
    let ref_order = match (va.ref_id(), vb.ref_id()) {
        (-1, -1) => Ordering::Equal,
        (-1, _) => Ordering::Greater,
        (_, -1) => Ordering::Less,
        (x, y) => x.cmp(&y),
    };
    ref_order
        .then((va.pos() + 1).cmp(&(vb.pos() + 1)))
        .then_with(|| (a_flag & 0x10 != 0).cmp(&(b_flag & 0x10 != 0)))
        .then_with(|| va.read_name().cmp(vb.read_name()))
        .then_with(|| a_flag.cmp(&b_flag))
        .then_with(|| a[9].cmp(&b[9]))
        .then_with(|| {
            i32::from_le_bytes(a[20..24].try_into().unwrap())
                .cmp(&i32::from_le_bytes(b[20..24].try_into().unwrap()))
        })
        .then_with(|| {
            i32::from_le_bytes(a[24..28].try_into().unwrap())
                .cmp(&i32::from_le_bytes(b[24..28].try_into().unwrap()))
        })
        .then_with(|| {
            i32::from_le_bytes(a[28..32].try_into().unwrap())
                .cmp(&i32::from_le_bytes(b[28..32].try_into().unwrap()))
        })
}

/// All records of the merged stream, resident in one arena.
struct Merged {
    arena: Vec<u8>,
    offsets: Vec<usize>,
    global_indices: Vec<Vec<usize>>,
}

impl Merged {
    fn len(&self) -> usize {
        self.offsets.len() - 1
    }
    fn get(&self, i: usize) -> &[u8] {
        &self.arena[self.offsets[i]..self.offsets[i + 1]]
    }

    fn get_mut(&mut self, i: usize) -> &mut [u8] {
        &mut self.arena[self.offsets[i]..self.offsets[i + 1]]
    }
}

fn merge_sets(sets: Vec<RecordSet>) -> io::Result<Merged> {
    let flags: Vec<Vec<u16>> = sets
        .iter()
        .map(|set| {
            (0..set.len())
                .map(|i| RecordView::new(set.record(i)).map(|v| v.flag()))
                .collect()
        })
        .collect::<io::Result<_>>()?;
    merge_sets_with_flags(sets, flags)
}

fn merge_sets_with_flags(sets: Vec<RecordSet>, flags: Vec<Vec<u16>>) -> io::Result<Merged> {
    if sets.len() != flags.len()
        || sets
            .iter()
            .zip(&flags)
            .any(|(set, flags)| set.len() != flags.len())
    {
        return Err(invalid("record flags do not match input sets".into()));
    }
    let mut heap = BinaryHeap::new();
    let mut next = vec![0usize; sets.len()];
    for (i, set) in sets.iter().enumerate() {
        if !set.is_empty() {
            let rec = set.record(0).to_vec();
            RecordView::new(&rec)?;
            heap.push(Head {
                rec,
                input: i,
                local: 0,
                original_flag: flags[i][0],
            });
            next[i] = 1;
        }
    }
    let mut m = Merged {
        arena: Vec::new(),
        offsets: vec![0],
        global_indices: sets.iter().map(|set| vec![0; set.len()]).collect(),
    };
    while let Some(Head {
        rec, input, local, ..
    }) = heap.pop()
    {
        m.global_indices[input][local] = m.len();
        m.arena.extend_from_slice(&rec);
        m.offsets.push(m.arena.len());
        if next[input] < sets[input].len() {
            let local = next[input];
            let rec = sets[input].record(local).to_vec();
            RecordView::new(&rec)?;
            next[input] += 1;
            heap.push(Head {
                rec,
                input,
                local,
                original_flag: flags[input][local],
            });
        }
    }
    Ok(m)
}

// ---------------------------------------------------------------------------
// Pass 1: read ends (§2)

const F: u8 = 0;
const R: u8 = 1;
const FF: u8 = 2;
const FR: u8 = 3;
const RR: u8 = 4;
const RF: u8 = 5;

fn orientation_byte(r1_neg: bool, r2_neg: bool) -> u8 {
    match (r1_neg, r2_neg) {
        (false, false) => FF,
        (false, true) => FR,
        (true, true) => RR,
        (true, false) => RF,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ReadEnds {
    library: i16,
    r1_ref: i32,
    r1_coord: i32,
    orientation: u8,
    r2_ref: i32,
    r2_coord: i32,
    r1_index: i64,
    r2_index: i64,
    score: i16,
    /// Tile, x, y from the read name and the read group's header index
    /// (Picard `PhysicalLocation`), for optical duplicates.
    loc: optical::Location,
    first_of_pair: bool,
    /// `orientationForOpticalDuplicates` of a pair (unset: 255).
    orientation_opt: u8,
}

impl ReadEnds {
    fn is_paired(&self) -> bool {
        self.r2_ref != -1
    }

    /// `ReadEndsMDComparator` (§3.1). Tile, x and y are -1 unless the read
    /// name parses (Illumina names); they then decide score ties.
    fn cmp(&self, o: &ReadEnds) -> Ordering {
        self.library
            .cmp(&o.library)
            .then(self.r1_ref.cmp(&o.r1_ref))
            .then(self.r1_coord.cmp(&o.r1_coord))
            .then(self.orientation.cmp(&o.orientation))
            .then(self.r2_ref.cmp(&o.r2_ref))
            .then(self.r2_coord.cmp(&o.r2_coord))
            .then(self.loc.tile.cmp(&o.loc.tile))
            .then(self.loc.x.cmp(&o.loc.x))
            .then(self.loc.y.cmp(&o.loc.y))
            .then(self.r1_index.cmp(&o.r1_index))
            .then(self.r2_index.cmp(&o.r2_index))
    }

    fn same_pair_key(&self, o: &ReadEnds) -> bool {
        self.same_fragment_key(o) && self.r2_ref == o.r2_ref && self.r2_coord == o.r2_coord
    }

    fn same_fragment_key(&self, o: &ReadEnds) -> bool {
        self.library == o.library
            && self.r1_ref == o.r1_ref
            && self.r1_coord == o.r1_coord
            && self.orientation == o.orientation
    }
}

/// POS (1-based) minus leading soft/hard clips (`SAMUtils.getUnclippedStart`).
fn unclipped_start(pos1: i32, cigar: &[u32]) -> i32 {
    let clipped: u32 = cigar
        .iter()
        .take_while(|&&c| matches!(c & 0xf, 4 | 5))
        .map(|c| c >> 4)
        .sum();
    pos1 - clipped as i32
}

/// Alignment end (1-based inclusive) plus trailing soft/hard clips.
fn unclipped_end(pos1: i32, cigar: &[u32]) -> i32 {
    let ref_len: u32 = cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|c| c >> 4)
        .sum();
    let clipped: u32 = cigar
        .iter()
        .rev()
        .take_while(|&&c| matches!(c & 0xf, 4 | 5))
        .map(|c| c >> 4)
        .sum();
    pos1 + ref_len as i32 - 1 + clipped as i32
}

/// SUM_OF_BASE_QUALITIES (§2.3).
fn duplicate_score(v: &RecordView, flag: u16) -> io::Result<i16> {
    let quals = v
        .quals()
        .ok_or_else(|| invalid("record truncated in qualities".into()))?;
    let sum: i64 = if quals.first() == Some(&0xff) {
        0 // missing qualities decode to an empty array
    } else {
        quals
            .iter()
            .map(|&q| q as i8 as i64)
            .filter(|&q| q >= 15)
            .sum()
    };
    let mut score = sum.min(16383) as i16;
    if flag & 0x200 != 0 {
        score = score.wrapping_add(-16384);
    }
    Ok(score)
}

struct Ends {
    pairs: Vec<ReadEnds>,
    frags: Vec<ReadEnds>,
}

trait Records: Sync {
    fn len(&self) -> usize;
    fn get(&self, i: usize) -> &[u8];
}

impl Records for Merged {
    fn len(&self) -> usize {
        self.len()
    }

    fn get(&self, i: usize) -> &[u8] {
        self.get(i)
    }
}

impl Records for RecordSet {
    fn len(&self) -> usize {
        self.len()
    }

    fn get(&self, i: usize) -> &[u8] {
        self.record(i)
    }
}

/// Library id per read group: ids start at 1 in order of first library (§2.4).
/// Per read-group ID: library id (ids start at 1 in order of first library,
/// §2.4) and the group's index among the header's @RG lines (Picard counts
/// `header.getReadGroups()` up to the record's RG for optical duplicates).
struct ReadGroups {
    library: HashMap<Vec<u8>, i16>,
    index: HashMap<Vec<u8>, i16>,
    count: i16,
}

fn read_group_libraries(header: &Header) -> ReadGroups {
    let text = String::from_utf8_lossy(&header.text);
    let mut libraries: Vec<String> = Vec::new();
    let mut rg_library = HashMap::new();
    let mut rg_index = HashMap::new();
    let mut count: i16 = 0;
    for line in text.lines().filter(|l| l.starts_with("@RG")) {
        let field = |k: &str| {
            line.split('\t')
                .find_map(|f| f.strip_prefix(k))
                .map(str::to_owned)
        };
        let (Some(id), lb) = (
            field("ID:"),
            field("LB:").unwrap_or_else(|| "Unknown Library".into()),
        ) else {
            continue;
        };
        let lib_id = match libraries.iter().position(|l| *l == lb) {
            Some(p) => p as i16 + 1,
            None => {
                libraries.push(lb);
                libraries.len() as i16
            }
        };
        rg_index.entry(id.clone().into_bytes()).or_insert(count);
        count += 1;
        rg_library.insert(id.into_bytes(), lib_id);
    }
    ReadGroups {
        library: rg_library,
        index: rg_index,
        count,
    }
}

/// The fragment end of a primary mapped record (§2.2 step 1).
fn fragment_end(v: &RecordView, index: usize, rgs: &ReadGroups) -> io::Result<ReadEnds> {
    let flag = v.flag();
    let neg = flag & 0x10 != 0;
    let pos1 = v.pos() + 1;
    let cigar: Vec<u32> = v.cigar().collect();
    let paired_mate_mapped = flag & 0x1 != 0 && flag & 0x8 == 0;
    Ok(ReadEnds {
        library: v
            .string_tag(b"RG")
            .and_then(|g| rgs.library.get(g).copied())
            .unwrap_or(1),
        r1_ref: v.ref_id(),
        r1_coord: if neg {
            unclipped_end(pos1, &cigar)
        } else {
            unclipped_start(pos1, &cigar)
        },
        orientation: if neg { R } else { F },
        r2_ref: if paired_mate_mapped {
            v.mate_ref_id()
        } else {
            -1
        },
        r2_coord: -1,
        r1_index: index as i64,
        r2_index: -1,
        score: duplicate_score(v, flag)?,
        loc: match optical::parse_read_name(v.read_name()) {
            Some((tile, x, y)) => optical::Location {
                // No RG tag: 0; an RG absent from the header: past the last.
                read_group: v
                    .string_tag(b"RG")
                    .map_or(0, |g| rgs.index.get(g).copied().unwrap_or(rgs.count)),
                tile,
                x,
                y,
            },
            None => optical::Location::NONE,
        },
        first_of_pair: flag & 0x40 != 0,
        orientation_opt: 255,
    })
}

fn build_fragments(
    records: &impl Records,
    header: &Header,
    threads: usize,
) -> io::Result<Vec<ReadEnds>> {
    let rgs = read_group_libraries(header);
    // Pass 1 stops at the first unplaced unmapped record. Coordinate order puts
    // every refID -1 record last, so find that boundary by binary search.
    let ref_of = |i: usize| RecordView::new(records.get(i)).map(|v| v.ref_id());
    let (mut lo, mut hi) = (0usize, records.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if ref_of(mid)? == -1 {
            hi = mid
        } else {
            lo = mid + 1
        }
    }
    let stop = lo;

    // Fragment ends for every primary mapped record, computed in parallel and
    // concatenated in index order.
    let threads = threads.max(1);
    let per = stop.div_ceil(threads).max(1);
    let parts: Vec<io::Result<Vec<ReadEnds>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..stop)
            .step_by(per)
            .map(|start| {
                let rgs = &rgs;
                scope.spawn(move || -> io::Result<Vec<ReadEnds>> {
                    let mut out = Vec::new();
                    for index in start..(start + per).min(stop) {
                        let v = RecordView::new(records.get(index))?;
                        if v.flag() & 0x904 == 0 {
                            out.push(fragment_end(&v, index, rgs)?);
                        }
                    }
                    Ok(out)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("pass-1 thread panicked"))
            .collect()
    });
    let mut frags = Vec::new();
    for p in parts {
        frags.extend(p?);
    }
    Ok(frags)
}

fn pair_fragments(records: &impl Records, frags: &[ReadEnds]) -> io::Result<Vec<ReadEnds>> {
    let mut pairs = Vec::new();
    let mut pending: HashMap<(&[u8], &[u8]), ReadEnds> = HashMap::new();
    for frag in frags {
        if frag.r2_ref == -1 {
            continue;
        }
        let index = frag.r1_index as usize;
        let v = RecordView::new(records.get(index))?;
        let rg = v
            .string_tag(b"RG")
            .ok_or_else(|| invalid(format!("record {index} has no RG tag")))?;
        let key = (rg, v.read_name());
        let Some(pair) = pending.remove(&key) else {
            pending.insert(key, *frag);
            continue;
        };
        pairs.push(join_mates(pair, frag));
    }
    Ok(pairs)
}

/// The pair end of a fragment end whose mate (`frag`) arrives later in
/// coordinate order (`buildSortedReadEndLists`).
fn join_mates(mut pair: ReadEnds, frag: &ReadEnds) -> ReadEnds {
    {
        let (a_neg, neg) = (pair.orientation == R, frag.orientation == R);
        // orientationForOpticalDuplicates goes by first then second of pair.
        pair.orientation_opt = if frag.first_of_pair {
            orientation_byte(neg, a_neg)
        } else {
            orientation_byte(a_neg, neg)
        };
        // orientationForOpticalDuplicates is unused: no read has a location.
        let (b_ref, b_coord) = (frag.r1_ref, frag.r1_coord);
        if b_ref > pair.r1_ref || (b_ref == pair.r1_ref && b_coord >= pair.r1_coord) {
            pair.r2_ref = b_ref;
            pair.r2_coord = b_coord;
            pair.r2_index = frag.r1_index;
            pair.orientation = orientation_byte(a_neg, neg);
            if pair.r2_ref == pair.r1_ref
                && pair.r2_coord == pair.r1_coord
                && pair.orientation == RF
            {
                pair.orientation = FR;
            }
        } else {
            pair.r2_ref = pair.r1_ref;
            pair.r2_coord = pair.r1_coord;
            pair.r2_index = pair.r1_index;
            pair.r1_ref = b_ref;
            pair.r1_coord = b_coord;
            pair.r1_index = frag.r1_index;
            pair.orientation = orientation_byte(neg, a_neg);
        }
        pair.score = pair.score.wrapping_add(frag.score);
    }
    pair
}

fn build_read_ends(records: &impl Records, header: &Header, threads: usize) -> io::Result<Ends> {
    let frags = build_fragments(records, header, threads)?;
    let pairs = pair_fragments(records, &frags)?;
    Ok(Ends { pairs, frags })
}

// ---------------------------------------------------------------------------
// Duplicate selection (§4)

/// `markDuplicatePairs`/fragment best: first entry with score > max so far.
fn best_of(chunk: &[ReadEnds]) -> usize {
    let mut best: Option<usize> = None;
    let mut max = 0i16;
    for (i, e) in chunk.iter().enumerate() {
        if best.is_none() || e.score > max {
            max = e.score;
            best = Some(i);
        }
    }
    best.unwrap()
}

struct Selection {
    duplicate: Vec<bool>,
    hist: metrics::SetHistograms,
    /// Optical duplicate pairs (single library).
    optical: u64,
}

fn select_duplicates(mut ends: Ends, n_records: usize) -> Selection {
    let mut sel = Selection {
        duplicate: vec![false; n_records],
        hist: Default::default(),
        optical: 0,
    };
    // `ReadEnds::cmp` is a total order (it ends in record indices), so any
    // sort gives the same sequence.
    use rayon::slice::ParallelSliceMut;
    ends.pairs.par_sort_unstable_by(ReadEnds::cmp);
    ends.frags.par_sort_unstable_by(ReadEnds::cmp);
    for chunk in ends.pairs.chunk_by(|a, b| a.same_pair_key(b)) {
        if chunk.len() == 1 {
            // addSingletonToCount
            *sel.hist.all.entry(1).or_default() += 1;
            *sel.hist.non_optical.entry(1).or_default() += 1;
        }
        if chunk.len() > 1 {
            let best = best_of(chunk);
            // trackOpticalDuplicates, then trackDuplicateCounts.
            let locs: Vec<optical::Location> = chunk.iter().map(|e| e.loc).collect();
            let orients: Vec<u8> = chunk.iter().map(|e| e.orientation_opt).collect();
            let n_opt =
                optical::track_optical_duplicates(&locs, &orients, best).unwrap_or_else(|| {
                    panic!("Picard: unexpected orientation in an FR/RF duplicate set")
                }) as u64;
            let len = chunk.len() as u64;
            sel.optical += n_opt;
            *sel.hist.all.entry(len).or_default() += 1;
            if len - n_opt > 0 {
                *sel.hist.non_optical.entry(len - n_opt).or_default() += 1;
            }
            if n_opt > 0 {
                *sel.hist.optical.entry(n_opt + 1).or_default() += 1;
            }
            for (i, e) in chunk.iter().enumerate() {
                if i != best {
                    sel.duplicate[e.r1_index as usize] = true;
                    sel.duplicate[e.r2_index as usize] = true;
                }
            }
        }
    }
    for chunk in ends.frags.chunk_by(|a, b| a.same_fragment_key(b)) {
        let contains_pairs = chunk.iter().any(ReadEnds::is_paired);
        let contains_frags = chunk.iter().any(|e| !e.is_paired());
        if chunk.len() <= 1 || !contains_frags {
            continue;
        }
        if contains_pairs {
            for e in chunk.iter().filter(|e| !e.is_paired()) {
                sel.duplicate[e.r1_index as usize] = true;
            }
        } else {
            let best = best_of(chunk);
            for (i, e) in chunk.iter().enumerate() {
                if i != best {
                    sel.duplicate[e.r1_index as usize] = true;
                }
            }
        }
    }
    sel
}

// ---------------------------------------------------------------------------
// Pass 2: htsjdk re-encoding (§5)

/// htsjdk `BinaryTagCodec.getIntegerType`: smallest type for the value.
fn integer_type(v: i64) -> IntType {
    if v > i32::MAX as i64 {
        IntType::U32
    } else if v > 65535 {
        IntType::I32
    } else if v > 32767 {
        IntType::U16
    } else if v > 255 {
        IntType::I16
    } else if v > 127 {
        IntType::U8
    } else if v >= -128 {
        IntType::I8
    } else if v >= -32768 {
        IntType::I16
    } else {
        IntType::I32
    }
}

fn tag_key(t: &[u8; 2]) -> u16 {
    (t[1] as u16) << 8 | t[0] as u16
}

/// Apply MarkDuplicates' record changes and htsjdk's encoding rules.
fn rewrite_record(raw: &[u8], duplicate: bool, pg_id: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    let mut rec = Record::decode(raw)?;
    rec.flag = if duplicate {
        rec.flag | 0x400
    } else {
        rec.flag & !0x400
    };
    rec.tags.retain(|(t, _)| t != b"DT");
    rec.tags.retain(|(t, _)| t != b"PG");
    rec.tags.push((*b"PG", Value::String(pg_id.to_vec())));
    // Sorted by (second char, first char); stable, so a duplicate key keeps the first.
    rec.tags.sort_by_key(|(t, _)| tag_key(t));
    rec.tags.dedup_by_key(|(t, _)| tag_key(t));
    for (_, v) in rec.tags.iter_mut() {
        match v {
            Value::Int(it, x) => *it = integer_type(*x),
            Value::Hex(_) => {
                return Err(invalid(
                    "H tags are re-encoded as B:c by htsjdk; not implemented".into(),
                ));
            }
            _ => {}
        }
    }
    if rec.qual.first() == Some(&0xff) {
        rec.qual.fill(0xff);
    }
    rec.encode(BinRule::Htsjdk, out)
}

fn pg_id(header: &Header) -> Vec<u8> {
    String::from_utf8_lossy(&header.text)
        .lines()
        .rfind(|l| l.starts_with("@PG"))
        .and_then(|l| l.split('\t').find_map(|f| f.strip_prefix("ID:")))
        .unwrap()
        .as_bytes()
        .to_vec()
}

fn rewrite_set(set: &RecordSet, pg_id: &[u8], threads: usize) -> io::Result<RecordSet> {
    let threads = threads.max(1).min(set.len().max(1));
    let per = set.len().div_ceil(threads).max(1);
    let parts: Vec<io::Result<Vec<u8>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..set.len())
            .step_by(per)
            .map(|start| {
                scope.spawn(move || -> io::Result<Vec<u8>> {
                    let mut part = Vec::new();
                    let mut enc = Vec::new();
                    for i in start..(start + per).min(set.len()) {
                        enc.clear();
                        rewrite_record(set.record(i), false, pg_id, &mut enc)?;
                        part.extend_from_slice(&(enc.len() as i32).to_le_bytes());
                        part.extend_from_slice(&enc);
                    }
                    Ok(part)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("overlap encode thread panicked"))
            .collect()
    });
    let mut records = Vec::with_capacity(set.bytes());
    for part in parts {
        records.extend(part?);
    }
    RecordSet::from_parts(set.header.clone(), &records)
}

/// Duplicate-independent work prepared while later chunks are still in flight.
pub struct PreparedSet {
    input: usize,
    set: RecordSet,
    original_flags: Vec<u16>,
    ends: Ends,
}

/// A fixed-header plan for incremental MarkDuplicates preparation.
pub struct OverlapPlan {
    input_headers: Vec<Header>,
    header: Header,
    pg_id: Vec<u8>,
}

impl OverlapPlan {
    pub fn new(input_headers: Vec<Header>, opts: &Options) -> io::Result<Self> {
        let header =
            header::merged_header(&input_headers, &opts.pg_version, &opts.pg_command_line)?;
        let pg_id = pg_id(&header);
        Ok(Self {
            input_headers,
            header,
            pg_id,
        })
    }

    pub fn prepare_set(
        &self,
        input: usize,
        set: RecordSet,
        threads: usize,
    ) -> io::Result<PreparedSet> {
        let expected = self
            .input_headers
            .get(input)
            .ok_or_else(|| invalid(format!("input index {input} is out of range")))?;
        if &set.header != expected {
            return Err(invalid(format!(
                "input {input} header changed during alignment"
            )));
        }
        let original_flags = (0..set.len())
            .map(|i| RecordView::new(set.record(i)).map(|v| v.flag()))
            .collect::<io::Result<Vec<_>>>()?;
        let ends = Ends {
            pairs: Vec::new(),
            frags: build_fragments(&set, &self.header, threads)?,
        };
        let set = rewrite_set(&set, &self.pg_id, threads)?;
        Ok(PreparedSet {
            input,
            set,
            original_flags,
            ends,
        })
    }

    pub fn finish(
        self,
        mut prepared: Vec<PreparedSet>,
        opts: &Options,
    ) -> io::Result<(RecordSet, DuplicationMetrics)> {
        prepared.sort_by_key(|p| p.input);
        if prepared.len() != self.input_headers.len()
            || prepared.iter().enumerate().any(|(i, p)| p.input != i)
        {
            return Err(invalid(
                "prepared inputs are incomplete or duplicated".into(),
            ));
        }

        let clock = std::time::Instant::now();
        let mut sets = Vec::with_capacity(prepared.len());
        let mut flags = Vec::with_capacity(prepared.len());
        let mut local_ends = Vec::with_capacity(prepared.len());
        for p in prepared {
            sets.push(p.set);
            flags.push(p.original_flags);
            local_ends.push(p.ends);
        }
        let mut merged = merge_sets_with_flags(sets, flags)?;
        eprintln!(
            "markdup stage overlap merge: {:.1}s",
            clock.elapsed().as_secs_f64()
        );

        let mut ends = Ends {
            pairs: Vec::new(),
            frags: Vec::new(),
        };
        for (input, mut local) in local_ends.into_iter().enumerate() {
            let map = &merged.global_indices[input];
            for end in &mut local.frags {
                end.r1_index = map[end.r1_index as usize] as i64;
            }
            ends.frags.append(&mut local.frags);
        }
        ends.frags.sort_unstable_by_key(|end| end.r1_index);
        ends.pairs = pair_fragments(&merged, &ends.frags)?;
        let selected_at = std::time::Instant::now();
        let sel = select_duplicates(ends, merged.len());
        eprintln!(
            "markdup stage overlap select: {:.1}s",
            selected_at.elapsed().as_secs_f64()
        );

        let finish_at = std::time::Instant::now();
        let mut metrics = metrics::Accumulator::new(&self.header);
        for (i, &duplicate) in sel.duplicate.iter().enumerate() {
            let raw = merged.get_mut(i);
            let mut flag = u16::from_le_bytes([raw[14], raw[15]]);
            if duplicate {
                flag |= 0x400;
            } else {
                flag &= !0x400;
            }
            raw[14..16].copy_from_slice(&flag.to_le_bytes());
            metrics.add(raw, duplicate)?;
        }
        let mut stream = Vec::with_capacity(
            merged.arena.len() + 4 * merged.len() + self.header.text.len() + (1 << 20),
        );
        self.header.write(&mut stream)?;
        for i in 0..merged.len() {
            let raw = merged.get(i);
            stream.extend_from_slice(&(raw.len() as i32).to_le_bytes());
            stream.extend_from_slice(raw);
        }
        let set = RecordSet::from_stream(stream)?;
        let result = metrics.finish(sel.optical);
        metrics::write_file(&opts.metrics, &opts.pg_command_line, &result, &sel.hist)?;
        eprintln!(
            "markdup stage overlap flags + metrics: {:.1}s",
            finish_at.elapsed().as_secs_f64()
        );
        Ok((set, result))
    }
}

/// Records re-encoded per parallel batch in pass 2.
const REWRITE_BATCH: usize = 1 << 16;

/// Pass 2 for one batch: re-encode `record(i)` for `i < duplicate.len()` on
/// `bufs.len()` threads, then add them to the metrics and append them
/// (`block_size`-prefixed) to `out` in order.
fn rewrite_batch<'a>(
    record: &(dyn Fn(usize) -> &'a [u8] + Sync),
    duplicate: &[bool],
    pg_id: &[u8],
    bufs: &mut [Vec<u8>],
    metrics: &mut metrics::Accumulator,
    out: &mut impl io::Write,
) -> io::Result<()> {
    let n = duplicate.len();
    let per = n.div_ceil(bufs.len());
    let results: Vec<io::Result<()>> = std::thread::scope(|scope| {
        let handles: Vec<_> = bufs
            .iter_mut()
            .enumerate()
            .map(|(k, buf)| {
                let (lo, hi) = ((k * per).min(n), ((k + 1) * per).min(n));
                scope.spawn(move || -> io::Result<()> {
                    buf.clear();
                    let mut enc = Vec::new();
                    for (i, &d) in (lo..hi).zip(&duplicate[lo..hi]) {
                        enc.clear();
                        rewrite_record(record(i), d, pg_id, &mut enc)?;
                        buf.extend_from_slice(&(enc.len() as i32).to_le_bytes());
                        buf.extend_from_slice(&enc);
                    }
                    Ok(())
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("encode thread panicked"))
            .collect()
    });
    results.into_iter().collect::<io::Result<()>>()?;
    let mut i = 0;
    for buf in bufs.iter() {
        let mut at = 0;
        while at < buf.len() {
            let len = i32::from_le_bytes(buf[at..at + 4].try_into().unwrap()) as usize;
            metrics.add(&buf[at + 4..at + 4 + len], duplicate[i])?;
            at += 4 + len;
            i += 1;
        }
        out.write_all(buf)?;
    }
    Ok(())
}

/// A k-way merge head over streamed inputs; its buffer is reused for the
/// input's next record.
struct StreamHead {
    rec: Vec<u8>,
    flag: u16,
    input: usize,
}

impl PartialEq for StreamHead {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for StreamHead {}
impl PartialOrd for StreamHead {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for StreamHead {
    /// The same order as [`Head`].
    fn cmp(&self, other: &Self) -> Ordering {
        coordinate_cmp_with_flags(&other.rec, other.flag, &self.rec, self.flag)
            .then(other.input.cmp(&self.input))
    }
}

/// Visit the records of coordinate-sorted BAMs in merged order (§1.5),
/// streaming them from disk.
fn merge_files(
    inputs: &[PathBuf],
    mut visit: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let mut cursors = inputs
        .iter()
        .map(|p| BamStream::open(p, STREAM_BATCH, 2, 2).map(RecordCursor::new))
        .collect::<io::Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for (input, cursor) in cursors.iter_mut().enumerate() {
        if let Some(rec) = cursor.peek()? {
            let flag = RecordView::new(rec)?.flag();
            heap.push(StreamHead {
                rec: rec.to_vec(),
                flag,
                input,
            });
            cursor.advance();
        }
    }
    while let Some(mut head) = heap.pop() {
        visit(&head.rec)?;
        let cursor = &mut cursors[head.input];
        if let Some(rec) = cursor.peek()? {
            head.flag = RecordView::new(rec)?.flag();
            head.rec.clear();
            head.rec.extend_from_slice(rec);
            cursor.advance();
            heap.push(head);
        }
    }
    Ok(())
}

/// FxHash (rustc's hasher): fast for short keys; the pending-mate map is
/// exact either way.
#[derive(Default, Clone, Copy)]
struct FxHasher(u64);

impl std::hash::Hasher for FxHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        const K: u64 = 0x51_7c_c1_b7_27_22_0a_95;
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            let w = u64::from_le_bytes(c.try_into().unwrap());
            self.0 = (self.0.rotate_left(5) ^ w).wrapping_mul(K);
        }
        for &b in chunks.remainder() {
            self.0 = (self.0.rotate_left(5) ^ u64::from(b)).wrapping_mul(K);
        }
    }
}

/// Mates waiting for their partner, keyed (as Picard does) by read group and
/// read name. Read groups are interned; the map is keyed by a 64-bit hash of
/// (read group, name) with the name stored to confirm, so lookups allocate
/// nothing; the rare hash collision goes to an exact overflow map.
#[derive(Default)]
struct PendingMates {
    read_groups: Vec<Vec<u8>>,
    by_hash: HashMap<u64, (u32, Box<[u8]>, ReadEnds), std::hash::BuildHasherDefault<FxHasher>>,
    overflow: HashMap<(u32, Vec<u8>), ReadEnds>,
}

impl PendingMates {
    /// The waiting mate of (rg, name), removed; or `None` after storing `frag`.
    fn take_or_insert(&mut self, rg: &[u8], name: &[u8], frag: ReadEnds) -> Option<ReadEnds> {
        use std::hash::Hasher;
        let rg_id = match self.read_groups.iter().position(|r| r == rg) {
            Some(i) => i as u32,
            None => {
                self.read_groups.push(rg.to_vec());
                (self.read_groups.len() - 1) as u32
            }
        };
        if !self.overflow.is_empty()
            && let Some(mate) = self.overflow.remove(&(rg_id, name.to_vec()))
        {
            return Some(mate);
        }
        let mut h = FxHasher(u64::from(rg_id));
        h.write(name);
        let key = h.finish();
        match self.by_hash.entry(key) {
            std::collections::hash_map::Entry::Occupied(e) => {
                let (r, n, _) = e.get();
                if *r == rg_id && **n == *name {
                    Some(e.remove().2)
                } else {
                    self.overflow.insert((rg_id, name.to_vec()), frag);
                    None
                }
            }
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert((rg_id, name.into(), frag));
                None
            }
        }
    }
}

/// Records per merged batch in streamed pass 1 (classified in parallel).
const PASS1_BATCH: usize = 1 << 16;

/// A pass-1 record: where the serial pass stops, records it skips, and
/// fragment ends (paired ones with the byte ranges of their RG value and
/// read name within the record, for mate pairing).
enum Pass1 {
    Unplaced,
    Skip,
    Frag(ReadEnds),
    Paired(
        ReadEnds,
        std::ops::Range<usize>,
        std::ops::Range<usize>,
        u64,
    ),
}

/// Only routes work; PendingMates still confirms the complete key, including
/// collisions. Hash the bytes rather than a shard-local interned RG id.
fn mate_shard_hash(rg: &[u8], name: &[u8]) -> u64 {
    use std::hash::Hasher;
    let mut h = FxHasher::default();
    h.write(rg);
    h.write(name);
    h.finish()
}

#[derive(Default)]
struct MateShard {
    pending: PendingMates,
    entries: Vec<(usize, Pass1)>,
    frags: Vec<ReadEnds>,
    pairs: Vec<ReadEnds>,
}

struct ShardedPendingMates {
    shards: Vec<MateShard>,
}

impl ShardedPendingMates {
    fn new(workers: usize) -> Self {
        Self {
            shards: (0..workers.max(1).next_power_of_two())
                .map(|_| MateShard::default())
                .collect(),
        }
    }

    fn add_batch(
        &mut self,
        classified: Vec<io::Result<Pass1>>,
        batch: &[u8],
        offsets: &[usize],
        stopped: &mut bool,
        ends: &mut Ends,
    ) -> io::Result<()> {
        use rayon::prelude::*;
        let mask = self.shards.len() - 1;
        // Validate in merged order: errors beyond the first unplaced record
        // are ignored, as in the serial sweep. Queues retain that same order
        // across batches, so each key pairs occurrences 1+2, 3+4, etc.
        for (i, kind) in classified.into_iter().enumerate() {
            if *stopped {
                break;
            }
            let kind = kind?;
            let shard = match &kind {
                Pass1::Unplaced => {
                    *stopped = true;
                    break;
                }
                Pass1::Skip => continue,
                Pass1::Frag(_) => i & mask,
                Pass1::Paired(_, _, _, hash) => *hash as usize & mask,
            };
            self.shards[shard].entries.push((i, kind));
        }
        self.shards.par_iter_mut().for_each(|shard| {
            for (i, kind) in shard.entries.drain(..) {
                match kind {
                    Pass1::Frag(frag) => shard.frags.push(frag),
                    Pass1::Paired(frag, rg, name, _) => {
                        shard.frags.push(frag);
                        let rec = &batch[offsets[i]..offsets[i + 1]];
                        if let Some(first) =
                            shard.pending.take_or_insert(&rec[rg], &rec[name], frag)
                        {
                            shard.pairs.push(join_mates(first, &frag));
                        }
                    }
                    Pass1::Unplaced | Pass1::Skip => unreachable!("only fragment ends are queued"),
                }
            }
        });
        // select_duplicates sorts both vectors by the total ReadEnds order;
        // shard completion/concatenation order cannot change duplicate picks.
        // Keep the per-shard allocations for the next batch.
        for shard in &mut self.shards {
            ends.frags.append(&mut shard.frags);
            ends.pairs.append(&mut shard.pairs);
        }
        Ok(())
    }
}

fn pass1_classify(rec: &[u8], index: usize, rgs: &ReadGroups) -> io::Result<Pass1> {
    let v = RecordView::new(rec)?;
    if v.ref_id() == -1 {
        return Ok(Pass1::Unplaced);
    }
    if v.flag() & 0x904 != 0 {
        return Ok(Pass1::Skip);
    }
    let frag = fragment_end(&v, index, rgs)?;
    if frag.r2_ref == -1 {
        return Ok(Pass1::Frag(frag));
    }
    let rg = v
        .string_tag(b"RG")
        .ok_or_else(|| invalid(format!("record {index} has no RG tag")))?;
    let name = v.read_name();
    // Both are slices of `rec`; keep their offsets.
    let at = |s: &[u8]| {
        let start = s.as_ptr() as usize - rec.as_ptr() as usize;
        start..start + s.len()
    };
    Ok(Pass1::Paired(
        frag,
        at(rg),
        at(name),
        mate_shard_hash(rg, name),
    ))
}

/// Records per batch handed over by each input's reader thread.
const STREAM_BATCH: usize = 1 << 14;

/// [`merge_files`] on its own thread, handing `visit` batches of at most
/// `size` merged records (`batch[offsets[i]..offsets[i + 1]]` is record `i`)
/// in order, so the merge overlaps the per-record work.
fn merge_files_batched(
    inputs: &[PathBuf],
    size: usize,
    mut visit: impl FnMut(&[u8], &[usize]) -> io::Result<()>,
) -> io::Result<()> {
    type Batch = (Vec<u8>, Vec<usize>);
    let (tx, rx) = std::sync::mpsc::sync_channel::<Batch>(2);
    let (back_tx, back_rx) = std::sync::mpsc::channel::<Batch>();
    std::thread::scope(|scope| {
        let producer = scope.spawn(move || -> io::Result<()> {
            let stopped = || io::Error::other("merge consumer stopped");
            let (mut batch, mut offsets) = (Vec::new(), vec![0]);
            merge_files(inputs, |rec| {
                batch.extend_from_slice(rec);
                offsets.push(batch.len());
                if offsets.len() > size {
                    let (mut b, mut o) = back_rx.try_recv().unwrap_or_default();
                    b.clear();
                    o.clear();
                    o.push(0);
                    std::mem::swap(&mut batch, &mut b);
                    std::mem::swap(&mut offsets, &mut o);
                    tx.send((b, o)).map_err(|_| stopped())?;
                }
                Ok(())
            })?;
            if offsets.len() > 1 {
                tx.send((batch, offsets)).map_err(|_| stopped())?;
            }
            Ok(())
        });
        let mut result = Ok(());
        while let Ok((batch, offsets)) = rx.recv() {
            result = visit(&batch, &offsets);
            if result.is_err() {
                break;
            }
            let _ = back_tx.send((batch, offsets));
        }
        drop(rx);
        let produced = producer.join().expect("merge thread panicked");
        result.and(produced)
    })
}

/// MarkDuplicates over coordinate-sorted chunk BAMs on disk (`opts.inputs`,
/// in chunk order), writing `opts.output`: the same result as [`mark_sets`]
/// with only the read ends in memory, for inputs too large to hold. The
/// inputs are merged twice, once for the read ends and once to rewrite.
pub fn mark_files(opts: &Options) -> io::Result<DuplicationMetrics> {
    mark_files_with(opts, |_| Ok(()))
}

/// md.bam's header for [`mark_files`] over `opts.inputs`.
pub fn files_header(opts: &Options) -> io::Result<Header> {
    let headers = opts
        .inputs
        .iter()
        .map(|p| BamStream::open(p, 1, 1, 1).map(|s| s.header.clone()))
        .collect::<io::Result<Vec<_>>>()?;
    header::merged_header(&headers, &opts.pg_version, &opts.pg_command_line)
}

/// [`mark_files`], also handing `sink` every batch of md.bam's records as
/// they are written: `block_size`-prefixed, duplicate flags set, in file
/// order, together making up md.bam's uncompressed stream after the header
/// ([`files_header`]). `sink` may take the buffer's contents (it is cleared
/// afterwards either way).
pub fn mark_files_with(
    opts: &Options,
    mut sink: impl FnMut(&mut Vec<u8>) -> io::Result<()>,
) -> io::Result<DuplicationMetrics> {
    let header = files_header(opts)?;
    let rgs = read_group_libraries(&header);
    let mut t = std::time::Instant::now();

    // Pass 1: the same fragment ends and mate pairs as `build_read_ends`,
    // preserving merged order within each mate key.
    let mut ends = Ends {
        pairs: Vec::new(),
        frags: Vec::new(),
    };
    let mut pending = ShardedPendingMates::new(rayon::current_num_threads());
    let (mut n_records, mut stopped) = (0usize, false);
    let mut consumer = std::time::Duration::ZERO;
    // Each batch: records classified in parallel (fragment ends are
    // independent per record), then routed in merged order and paired in
    // parallel shards. Errors surface only for records the serial pass reaches.
    merge_files_batched(&opts.inputs, PASS1_BATCH, |batch, offsets| {
        use rayon::prelude::*;
        let started = std::time::Instant::now();
        let base = n_records;
        let n = offsets.len() - 1;
        let record = |i: usize| &batch[offsets[i]..offsets[i + 1]];
        let classified: Vec<io::Result<Pass1>> = (0..n)
            .into_par_iter()
            .with_min_len(1024)
            .map(|i| pass1_classify(record(i), base + i, &rgs))
            .collect();
        n_records += n;
        pending.add_batch(classified, batch, offsets, &mut stopped, &mut ends)?;
        consumer += started.elapsed();
        Ok(())
    })?;
    drop(pending);
    eprintln!(
        "markdup stage stream pass1 ({n_records} records, {} fragment ends, {} pairs): {:.1}s (consumer busy {:.1}s)",
        ends.frags.len(),
        ends.pairs.len(),
        t.elapsed().as_secs_f64(),
        consumer.as_secs_f64()
    );
    t = std::time::Instant::now();
    let sel = select_duplicates(ends, n_records);
    eprintln!(
        "markdup stage stream select: {:.1}s",
        t.elapsed().as_secs_f64()
    );
    t = std::time::Instant::now();

    // Pass 2: merge again, re-encode in parallel batches, compress as we go.
    let pg_id = pg_id(&header);
    let threads = opts.threads.max(1);
    // md.bam, or the same stream as a light spill file (`.nspill` output).
    let mut out = nomics_io::bamstream::StreamWriter::create(&opts.output, threads)?;
    // The header is its own write: BGZF blocks are cut by byte count, not
    // by write calls, so the file is the same.
    let mut stream = Vec::new();
    header.write(&mut stream)?;
    out.write_all(&stream)?;
    stream.clear();
    let mut m = metrics::Accumulator::new(&header);
    let mut bufs: Vec<Vec<u8>> = vec![Vec::new(); threads];
    let mut done = 0usize;
    merge_files_batched(&opts.inputs, REWRITE_BATCH, |batch, offsets| {
        let n = offsets.len() - 1;
        if done + n > sel.duplicate.len() {
            return Err(invalid(format!(
                "pass 2 saw more than pass 1's {n_records} records"
            )));
        }
        rewrite_batch(
            &|i| &batch[offsets[i]..offsets[i + 1]],
            &sel.duplicate[done..done + n],
            &pg_id,
            &mut bufs,
            &mut m,
            &mut stream,
        )?;
        out.write_all(&stream)?;
        sink(&mut stream)?;
        stream.clear();
        done += n;
        Ok(())
    })?;
    if done != n_records {
        return Err(invalid(format!(
            "pass 2 saw {done} records, pass 1 {n_records}"
        )));
    }
    out.finish()?;
    eprintln!(
        "markdup stage stream pass2 + write ({threads} threads): {:.1}s",
        t.elapsed().as_secs_f64()
    );
    let result = m.finish(sel.optical);
    metrics::write_file(&opts.metrics, &opts.pg_command_line, &result, &sel.hist)?;
    Ok(result)
}

/// MarkDuplicates, then write `opts.output`.
pub fn run(opts: &Options) -> io::Result<DuplicationMetrics> {
    let (set, metrics) = mark(opts)?;
    let t = std::time::Instant::now();
    set.write(&opts.output, opts.threads)?;
    eprintln!(
        "markdup stage write md.bam: {:.1}s",
        t.elapsed().as_secs_f64()
    );
    Ok(metrics)
}

/// MarkDuplicates in memory: returns the output records (as they would be
/// written to md.bam) and writes the metrics file.
pub fn mark(opts: &Options) -> io::Result<(RecordSet, DuplicationMetrics)> {
    let sets = opts
        .inputs
        .iter()
        .map(|p| RecordSet::load(p, opts.threads.max(1)))
        .collect::<io::Result<Vec<_>>>()?;
    mark_sets(sets, opts)
}

/// MarkDuplicates over coordinate-sorted record sets already in memory.
pub fn mark_sets(
    sets: Vec<RecordSet>,
    opts: &Options,
) -> io::Result<(RecordSet, DuplicationMetrics)> {
    let headers: Vec<Header> = sets.iter().map(|set| set.header.clone()).collect();
    let header = header::merged_header(&headers, &opts.pg_version, &opts.pg_command_line)?;
    let clock = std::time::Instant::now();
    let lap = |name: &str, t: &mut std::time::Instant| {
        eprintln!("markdup stage {name}: {:.1}s", t.elapsed().as_secs_f64());
        *t = std::time::Instant::now();
    };
    let mut t = clock;
    let merged = merge_sets(sets)?;
    lap("merge", &mut t);
    let ends = build_read_ends(&merged, &header, opts.threads)?;
    lap("pass1 read ends", &mut t);
    let sel = select_duplicates(ends, merged.len());
    lap("select", &mut t);

    let pg_id = String::from_utf8_lossy(&header.text)
        .lines()
        .rfind(|l| l.starts_with("@PG"))
        .and_then(|l| {
            l.split('\t')
                .find_map(|f| f.strip_prefix("ID:"))
                .map(str::to_owned)
        })
        .unwrap();
    let threads = opts.threads.max(1);
    // The output stream: header, then every re-encoded record in order.
    let mut stream =
        Vec::with_capacity(merged.arena.len() + 4 * merged.len() + header.text.len() + (1 << 20));
    header.write(&mut stream)?;
    let mut m = metrics::Accumulator::new(&header);
    let mut bufs: Vec<Vec<u8>> = vec![Vec::new(); threads];
    let mut t_encode = std::time::Duration::ZERO;
    for start in (0..merged.len()).step_by(REWRITE_BATCH) {
        let end = (start + REWRITE_BATCH).min(merged.len());
        let t0 = std::time::Instant::now();
        rewrite_batch(
            &|i| merged.get(start + i),
            &sel.duplicate[start..end],
            pg_id.as_bytes(),
            &mut bufs,
            &mut m,
            &mut stream,
        )?;
        t_encode += t0.elapsed();
    }
    drop(merged);
    let set = RecordSet::from_stream(stream)?;
    eprintln!(
        "markdup stage pass2 re-encode ({threads} threads): {:.1}s",
        t_encode.as_secs_f64()
    );
    lap("pass2 total", &mut t);
    let result = m.finish(sel.optical);
    metrics::write_file(&opts.metrics, &opts.pg_command_line, &result, &sel.hist)?;
    Ok((set, result))
}

#[cfg(test)]
mod tests {
    #[test]
    fn pending_mates_pair_by_read_group_and_name_through_collisions() {
        let end = |i: i64| ReadEnds {
            library: 1,
            r1_ref: 0,
            r1_coord: 1,
            orientation: FR,
            r2_ref: 0,
            r2_coord: 5,
            r1_index: i,
            r2_index: i + 1,
            score: 1,
            loc: optical::Location::NONE,
            first_of_pair: true,
            orientation_opt: 255,
        };
        let mut p = PendingMates::default();
        assert!(p.take_or_insert(b"rg1", b"a", end(1)).is_none());
        assert!(p.take_or_insert(b"rg2", b"a", end(2)).is_none());
        // Force a hash collision: occupy "b"'s slot with another name.
        let mut h = FxHasher(0);
        std::hash::Hasher::write(&mut h, b"b");
        let key = std::hash::Hasher::finish(&h);
        p.by_hash.insert(key, (0, b"zzz"[..].into(), end(9)));
        assert!(p.take_or_insert(b"rg1", b"b", end(3)).is_none());
        assert_eq!(p.overflow.len(), 1);
        assert_eq!(p.take_or_insert(b"rg1", b"b", end(4)).unwrap().r1_index, 3);
        assert_eq!(p.take_or_insert(b"rg2", b"a", end(5)).unwrap().r1_index, 2);
        assert_eq!(p.take_or_insert(b"rg1", b"a", end(6)).unwrap().r1_index, 1);
        assert!(p.overflow.is_empty());
    }

    #[test]
    fn sharded_pairing_matches_serial_across_batches_and_odd_counts() {
        let mut keys = Vec::new();
        for rg in 0..3 {
            for name in 0..19 {
                // Reuse names across RGs; include singletons and 3+ occurrences.
                for _ in 0..1 + (rg + name) % 7 {
                    keys.push((
                        format!("rg{rg}").into_bytes(),
                        format!("read{name}").into_bytes(),
                    ));
                }
            }
        }
        let mut seed = 123456789u64;
        for i in (1..keys.len()).rev() {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            keys.swap(i, seed as usize % (i + 1));
        }
        let frags: Vec<_> = (0..keys.len())
            .map(|i| ReadEnds {
                library: (i % 3) as i16,
                r1_ref: (i % 2) as i32,
                r1_coord: (i % 5) as i32,
                orientation: if i % 2 == 0 { F } else { R },
                r2_ref: if i % 11 == 0 { -1 } else { 0 },
                r2_coord: -1,
                r1_index: i as i64,
                r2_index: -1,
                score: (i * 127) as i16,
                loc: optical::Location {
                    read_group: (i % 3) as i16,
                    tile: (i % 7) as i16,
                    x: i as i32,
                    y: -(i as i32),
                },
                first_of_pair: i % 3 == 0,
                orientation_opt: 255,
            })
            .collect();
        let mut serial = PendingMates::default();
        let mut want_pairs = Vec::new();
        for ((rg, name), frag) in keys.iter().zip(&frags) {
            if !frag.is_paired() {
                continue;
            }
            if let Some(first) = serial.take_or_insert(rg, name, *frag) {
                want_pairs.push(join_mates(first, frag));
            }
        }
        let waiting = |p: &PendingMates| {
            p.by_hash
                .values()
                .map(|(rg, name, frag)| {
                    ((p.read_groups[*rg as usize].clone(), name.to_vec()), *frag)
                })
                .chain(p.overflow.iter().map(|((rg, name), frag)| {
                    ((p.read_groups[*rg as usize].clone(), name.clone()), *frag)
                }))
                .collect::<HashMap<_, _>>()
        };
        let want_pending = waiting(&serial);
        assert!(!want_pending.is_empty());
        want_pairs.sort_unstable_by(ReadEnds::cmp);
        let mut want_frags = frags.clone();
        want_frags.sort_unstable_by(ReadEnds::cmp);

        for workers in [1, 2, 3, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                for batch_size in [1, 7, 64, keys.len()] {
                    let mut sharded = ShardedPendingMates::new(workers);
                    let mut ends = Ends {
                        pairs: Vec::new(),
                        frags: Vec::new(),
                    };
                    let mut stopped = false;
                    for start in (0..keys.len()).step_by(batch_size) {
                        let mut batch = Vec::new();
                        let mut offsets = vec![0];
                        let mut classified = Vec::new();
                        for i in start..(start + batch_size).min(keys.len()) {
                            let (rg, name) = &keys[i];
                            batch.extend_from_slice(rg);
                            batch.extend_from_slice(name);
                            offsets.push(batch.len());
                            classified.push(Ok(if frags[i].is_paired() {
                                Pass1::Paired(
                                    frags[i],
                                    0..rg.len(),
                                    rg.len()..rg.len() + name.len(),
                                    mate_shard_hash(rg, name),
                                )
                            } else {
                                Pass1::Frag(frags[i])
                            }));
                        }
                        sharded
                            .add_batch(classified, &batch, &offsets, &mut stopped, &mut ends)
                            .unwrap();
                    }
                    assert!(!stopped);
                    ends.frags.sort_unstable_by(ReadEnds::cmp);
                    ends.pairs.sort_unstable_by(ReadEnds::cmp);
                    assert_eq!(ends.frags, want_frags);
                    assert_eq!(ends.pairs, want_pairs);
                    let got_pending: HashMap<_, _> = sharded
                        .shards
                        .iter()
                        .flat_map(|shard| waiting(&shard.pending))
                        .collect();
                    assert_eq!(got_pending, want_pending);
                }
            });
        }
    }

    #[test]
    fn sharded_pairing_preserves_stop_and_error_boundary() {
        let set = pair_set(0, b"before-stop", 30);
        let rgs = read_group_libraries(&set.header);
        let mut want = build_read_ends(&set, &set.header, 1).unwrap();
        let mut batch = Vec::new();
        let mut offsets = vec![0];
        let mut classified = Vec::new();
        for i in 0..set.len() {
            batch.extend_from_slice(set.record(i));
            offsets.push(batch.len());
            classified.push(pass1_classify(set.record(i), i, &rgs));
        }
        classified.extend([
            Ok(Pass1::Skip),
            Ok(Pass1::Unplaced),
            Ok(Pass1::Frag(want.frags[0])),
            Err(invalid("ignored".into())),
        ]);
        let mut sharded = ShardedPendingMates::new(2);
        let mut ends = Ends {
            pairs: Vec::new(),
            frags: Vec::new(),
        };
        let mut stopped = false;
        sharded
            .add_batch(classified, &batch, &offsets, &mut stopped, &mut ends)
            .unwrap();
        assert!(stopped);
        ends.frags.sort_unstable_by(ReadEnds::cmp);
        want.frags.sort_unstable_by(ReadEnds::cmp);
        assert_eq!(ends.frags, want.frags);
        assert_eq!(ends.pairs, want.pairs);
        sharded
            .add_batch(
                vec![Err(invalid("also ignored in later batch".into()))],
                &[],
                &[],
                &mut stopped,
                &mut ends,
            )
            .unwrap();
        stopped = false;
        let error = sharded
            .add_batch(
                vec![
                    Ok(Pass1::Skip),
                    Err(invalid("first error".into())),
                    Ok(Pass1::Unplaced),
                    Err(invalid("later error".into())),
                ],
                &[],
                &[],
                &mut stopped,
                &mut ends,
            )
            .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert_eq!(error.to_string(), "first error");
        assert!(!stopped);
    }

    use super::*;
    use nomics_io::bam::Reference;

    fn test_header(chunk: usize) -> Header {
        Header {
            text: format!(
                "@HD\tVN:1.5\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
                 @RG\tID:L1\tSM:sample\tLB:lib\n\
                 @PG\tID:nomics-align\tPN:nomics-align\tCL:chunk {chunk}\n"
            )
            .into_bytes(),
            refs: vec![Reference {
                name: b"chr1".to_vec(),
                len: 1000,
            }],
        }
    }

    fn pair_set(chunk: usize, name: &[u8], qual: u8) -> RecordSet {
        let make = |first: bool| Record {
            ref_id: 0,
            pos: if first { 99 } else { 199 },
            mapq: 60,
            flag: if first { 0x41 } else { 0x91 },
            next_ref_id: 0,
            next_pos: if first { 199 } else { 99 },
            tlen: if first { 110 } else { -110 },
            read_name: name.to_vec(),
            cigar: vec![10 << 4],
            seq: vec![b'A'; 10],
            qual: vec![qual; 10],
            tags: vec![(*b"RG", Value::String(b"L1".to_vec()))],
        };
        let mut records = Vec::new();
        let mut raw = Vec::new();
        for rec in [make(true), make(false)] {
            raw.clear();
            rec.encode(BinRule::Sam, &mut raw).unwrap();
            records.extend_from_slice(&(raw.len() as i32).to_le_bytes());
            records.extend_from_slice(&raw);
        }
        RecordSet::from_parts(test_header(chunk), &records).unwrap()
    }

    fn synthetic_sets() -> Vec<RecordSet> {
        vec![pair_set(0, b"first", 30), pair_set(1, b"second", 30)]
    }

    fn temp_metrics(label: &str) -> PathBuf {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "nomics-markdup-{label}-{}-{n}.metrics",
            std::process::id()
        ))
    }

    #[test]
    fn unclipped_coordinates() {
        // POS=1000 CIGAR=5S146M forward -> 995; POS=1000 3H140M8S reverse -> 1147 (spec §2.2).
        assert_eq!(unclipped_start(1000, &[(5 << 4) | 4, 146 << 4]), 995);
        assert_eq!(
            unclipped_end(1000, &[(3 << 4) | 5, 140 << 4, (8 << 4) | 4]),
            1147
        );
    }

    #[test]
    fn integer_types_follow_htsjdk() {
        let cases = [
            (0, 'c'),
            (127, 'c'),
            (128, 'C'),
            (255, 'C'),
            (256, 's'),
            (32767, 's'),
            (32768, 'S'),
            (65535, 'S'),
            (65536, 'i'),
            (-1, 'c'),
            (-128, 'c'),
            (-129, 's'),
            (-32768, 's'),
            (-32769, 'i'),
        ];
        for (v, want) in cases {
            let got = match integer_type(v) {
                IntType::I8 => 'c',
                IntType::U8 => 'C',
                IntType::I16 => 's',
                IntType::U16 => 'S',
                IntType::I32 => 'i',
                IntType::U32 => 'I',
            };
            assert_eq!(got, want, "value {v}");
        }
    }

    #[test]
    fn tag_keys_order_like_htsjdk() {
        let mut tags: Vec<[u8; 2]> = vec![
            *b"NM", *b"MD", *b"MC", *b"AS", *b"XS", *b"RG", *b"PG", *b"SA", *b"XA", *b"MQ",
        ];
        tags.sort_by_key(tag_key);
        assert_eq!(
            tags,
            [
                *b"SA", *b"XA", *b"MC", *b"MD", *b"PG", *b"RG", *b"NM", *b"MQ", *b"AS", *b"XS"
            ]
        );
    }

    #[test]
    fn best_is_first_strict_maximum() {
        let e = |score, idx| ReadEnds {
            library: 1,
            r1_ref: 0,
            r1_coord: 1,
            orientation: FR,
            r2_ref: 0,
            r2_coord: 5,
            r1_index: idx,
            r2_index: idx + 1,
            score,
            loc: optical::Location::NONE,
            first_of_pair: true,
            orientation_opt: 255,
        };
        assert_eq!(best_of(&[e(5, 0), e(9, 2), e(9, 4)]), 1);
        assert_eq!(best_of(&[e(-3, 0), e(-5, 2)]), 0);
    }

    #[test]
    fn streamed_files_match_in_memory() {
        let one_metrics = temp_metrics("one-mem");
        let file_metrics = temp_metrics("files");
        let dir = std::env::temp_dir().join(format!("nomics-markdup-files-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let sets = synthetic_sets();
        let inputs: Vec<PathBuf> = sets
            .iter()
            .enumerate()
            .map(|(i, set)| {
                let path = dir.join(format!("{i}.bam"));
                set.write(&path, 1).unwrap();
                path
            })
            .collect();
        let opts = |metrics, inputs, output| Options {
            inputs,
            output,
            metrics,
            pg_version: "Version:test".into(),
            pg_command_line: "nomics test".into(),
            threads: 2,
        };
        let (one, one_summary) =
            mark_sets(sets, &opts(one_metrics.clone(), Vec::new(), PathBuf::new())).unwrap();
        let md = dir.join("md.bam");
        let file_summary = mark_files(&opts(file_metrics.clone(), inputs, md.clone())).unwrap();
        let streamed = RecordSet::load(&md, 1).unwrap();
        assert_eq!(one_summary, file_summary);
        assert_eq!(one.header, streamed.header);
        assert_eq!(one.len(), streamed.len());
        for i in 0..one.len() {
            assert_eq!(one.record(i), streamed.record(i), "record {i}");
        }
        assert_eq!(
            std::fs::read(&one_metrics).unwrap(),
            std::fs::read(&file_metrics).unwrap()
        );
        std::fs::remove_dir_all(dir).unwrap();
        std::fs::remove_file(one_metrics).unwrap();
        std::fs::remove_file(file_metrics).unwrap();
    }

    #[test]
    fn overlap_matches_one_shot_across_chunks() {
        let one_metrics = temp_metrics("one");
        let overlap_metrics = temp_metrics("overlap");
        let opts = |metrics| Options {
            inputs: Vec::new(),
            output: PathBuf::new(),
            metrics,
            pg_version: "Version:test".into(),
            pg_command_line: "nomics test".into(),
            threads: 2,
        };
        let one_opts = opts(one_metrics.clone());
        let overlap_opts = opts(overlap_metrics.clone());
        let (one, one_summary) = mark_sets(synthetic_sets(), &one_opts).unwrap();

        let overlap_sets = synthetic_sets();
        let headers = overlap_sets.iter().map(|s| s.header.clone()).collect();
        let plan = OverlapPlan::new(headers, &overlap_opts).unwrap();
        let prepared = overlap_sets
            .into_iter()
            .enumerate()
            .map(|(i, set)| plan.prepare_set(i, set, 2).unwrap())
            .collect();
        let (overlap, overlap_summary) = plan.finish(prepared, &overlap_opts).unwrap();

        assert_eq!(one_summary, overlap_summary);
        assert_eq!(one.header, overlap.header);
        assert_eq!(one.len(), overlap.len());
        for i in 0..one.len() {
            assert_eq!(one.record(i), overlap.record(i), "record {i}");
        }
        assert_eq!(
            RecordView::new(overlap.record(0)).unwrap().flag() & 0x400,
            0
        );
        assert_ne!(
            RecordView::new(overlap.record(1)).unwrap().flag() & 0x400,
            0
        );
        assert_eq!(
            std::fs::read(&one_metrics).unwrap(),
            std::fs::read(&overlap_metrics).unwrap()
        );
        std::fs::remove_file(one_metrics).unwrap();
        std::fs::remove_file(overlap_metrics).unwrap();
    }
}
