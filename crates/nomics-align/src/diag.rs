//! `NOMICS_DIAG=1` region invariants (docs/DIAG.md).
//!
//! Every [`Region`] must have `0 <= rid < n_seqs`, `0 <= rb <= re <= 2 *
//! l_pac` and `0 <= qb <= qe <= l_query`. bwa never makes another (chains
//! keep `rid >= 0`, extension copies the chain's rid, rescue the anchor's),
//! so a violation means corrupted memory. The checks run at the stage
//! boundaries listed in docs/DIAG.md; the first failure prints the read, its
//! batch and chunk, the stage, the bad region and its neighbours, and the
//! addresses of the region buffer and the HostBufs nearest it, then aborts
//! the process. With `NOMICS_DIAG` unset none of this runs.
//!
//! Also here: the batch-range knobs for replays (`NOMICS_ALIGN_SKIP_BATCHES`,
//! `NOMICS_ALIGN_MAX_BATCHES`), which work with or without diagnostics.
use std::cell::RefCell;
use std::sync::Mutex;

use crate::extend::Region;
use crate::index::Index;
use crate::sam::Read;

/// Whether `NOMICS_DIAG` is set ([`nomics_mem::diag`]).
#[inline]
pub fn on() -> bool {
    nomics_mem::diag()
}

/// Why `r` is not a region bwa can produce for a read of `l_query` bases,
/// or `None`.
pub fn region_fault(idx: &Index, l_query: usize, r: &Region) -> Option<&'static str> {
    region_fault_in(idx.bns.anns.len(), idx.bns.l_pac, l_query, r)
}

/// [`region_fault`] for a reference of `n_seqs` contigs and `l_pac` bases.
pub fn region_fault_in(
    n_seqs: usize,
    l_pac: i64,
    l_query: usize,
    r: &Region,
) -> Option<&'static str> {
    if r.rid < 0 || r.rid as usize >= n_seqs {
        return Some("rid outside 0..n_seqs");
    }
    if r.rb < 0 || r.re > 2 * l_pac {
        return Some("rb/re outside 0..=2*l_pac");
    }
    if r.rb > r.re {
        return Some("rb > re");
    }
    if r.qb < 0 || r.qe < r.qb || r.qe as usize > l_query {
        return Some("qb/qe outside 0..=l_query or qb > qe");
    }
    None
}

/// The first region of `regs` that fails [`region_fault`].
pub fn first_fault(idx: &Index, l_query: usize, regs: &[Region]) -> Option<(usize, &'static str)> {
    first_fault_in(idx.bns.anns.len(), idx.bns.l_pac, l_query, regs)
}

fn first_fault_in(
    n_seqs: usize,
    l_pac: i64,
    l_query: usize,
    regs: &[Region],
) -> Option<(usize, &'static str)> {
    regs.iter()
        .enumerate()
        .find_map(|(i, r)| region_fault_in(n_seqs, l_pac, l_query, r).map(|why| (i, why)))
}

// ---- Context: which read, batch and chunk ----

thread_local! {
    /// The batch this thread aligns (set where batch threads start).
    static CHUNK: RefCell<Option<BatchId>> = const { RefCell::new(None) };
    /// The read this pool thread is extending: its address (to find its
    /// batch) and a copy of its name.
    static READ: RefCell<(usize, Vec<u8>)> = const { RefCell::new((0, Vec::new())) };
}

struct BatchEntry {
    reads: usize,
    n: usize,
    first_pair: u64,
    id: Option<BatchId>,
}

static BATCHES: Mutex<Vec<BatchEntry>> = Mutex::new(Vec::new());
/// Where a batch was cut: its chunk (stream index), its ordinal among the
/// chunk's batches as the reader cut them (skipped ones included, so it is
/// the `NOMICS_ALIGN_SKIP_BATCHES` value that replays it first), and the
/// chunk's read group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchId {
    pub chunk: usize,
    pub ordinal: usize,
    pub rg: String,
}

/// Note the batch this thread aligns next (diagnostics only).
pub fn set_chunk(chunk: usize, ordinal: usize, rg_id: &[u8]) {
    if on() {
        CHUNK.with(|c| {
            *c.borrow_mut() = Some(BatchId {
                chunk,
                ordinal,
                rg: String::from_utf8_lossy(rg_id).into_owned(),
            })
        });
    }
}

/// Note that this thread now extends `r` (diagnostics only; the name is
/// copied into a reused buffer, so this does not allocate per read).
#[inline]
pub fn at_read(r: &Read) {
    if on() {
        READ.with(|c| {
            let mut c = c.borrow_mut();
            c.0 = r as *const Read as usize;
            c.1.clear();
            c.1.extend_from_slice(&r.name);
        });
    }
}

/// A batch registered for diagnostic reports while it aligns.
pub struct BatchGuard {
    reads: usize,
}

/// Register the batch `reads` (first pair id `first_pair`) for reports, until
/// the guard drops. `None` without diagnostics.
pub fn batch(reads: &[Read], first_pair: u64) -> Option<BatchGuard> {
    if !on() {
        return None;
    }
    let entry = BatchEntry {
        reads: reads.as_ptr() as usize,
        n: reads.len(),
        first_pair,
        id: CHUNK.with(|c| c.borrow_mut().take()),
    };
    let guard = BatchGuard { reads: entry.reads };
    BATCHES
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .push(entry);
    Some(guard)
}

impl Drop for BatchGuard {
    fn drop(&mut self) {
        let mut batches = BATCHES.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(at) = batches.iter().position(|b| b.reads == self.reads) {
            batches.swap_remove(at);
        }
    }
}

fn describe_batch(b: &BatchEntry) -> String {
    let at = match &b.id {
        Some(id) => format!(
            "chunk {} batch {} (as for NOMICS_ALIGN_SKIP_BATCHES={}; chunk = stream index in this run, read group {})",
            id.chunk, id.ordinal, id.ordinal, id.rg
        ),
        None => "batch of the single-stream path (chunk and ordinal unknown)".to_owned(),
    };
    format!("{at}, first pair id {} ({} reads)", b.first_pair, b.n)
}

/// The batch holding the `Read` at address `addr`, and the read's index.
fn batch_of(addr: usize) -> String {
    let size = std::mem::size_of::<Read>();
    let batches = BATCHES.lock().unwrap_or_else(|e| e.into_inner());
    for b in batches.iter() {
        if addr >= b.reads && addr < b.reads + b.n * size {
            let i = (addr - b.reads) / size;
            return format!(
                "read index {i} in {} (pair id {}, mate {})",
                describe_batch(b),
                b.first_pair + (i / 2) as u64,
                i % 2 + 1
            );
        }
    }
    format!("read at {addr:#x} not in a registered batch")
}

/// Live batches whose pair ids include `id` (ids restart per chunk, so a
/// chunk boundary can give two).
fn batches_with_pair(id: u64) -> String {
    let batches = BATCHES.lock().unwrap_or_else(|e| e.into_inner());
    let found: Vec<String> = batches
        .iter()
        .filter(|b| id >= b.first_pair && id < b.first_pair + (b.n / 2) as u64)
        .map(|b| {
            format!(
                "read indices {}/{} in {}",
                2 * (id - b.first_pair),
                2 * (id - b.first_pair) + 1,
                describe_batch(b)
            )
        })
        .collect();
    if found.is_empty() {
        "no registered batch holds this pair id".to_owned()
    } else {
        found.join("; ")
    }
}

// ---- Reports ----

fn region_lines(regs: &[Region], at: usize) -> String {
    let mut out = String::new();
    let lo = at.saturating_sub(2);
    let hi = (at + 3).min(regs.len());
    for (i, r) in regs.iter().enumerate().take(hi).skip(lo) {
        let mark = if i == at { ">>" } else { "  " };
        out.push_str(&format!("  {mark} [{i}] {r:?}\n"));
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn fail(
    stage: &str,
    why: &str,
    idx: &Index,
    read: &str,
    l_query: usize,
    regs: &[Region],
    capacity: Option<usize>,
    at: usize,
) -> ! {
    let ptr = regs.as_ptr() as usize;
    let size = std::mem::size_of::<Region>();
    let capacity = capacity.map_or("unknown".to_owned(), |c| c.to_string());
    let report = format!(
        "region invariant failed {stage}: {why}\n\
         {read}\n\
         read length {l_query}; n_seqs {}, l_pac {}\n\
         region {at} of {} (Vec buffer {ptr:#x}..{:#x}, capacity {capacity}, {size} bytes per region; bad region at {:#x}):\n\
         {}{}",
        idx.bns.anns.len(),
        idx.bns.l_pac,
        regs.len(),
        ptr + std::mem::size_of_val(regs),
        ptr + at * size,
        region_lines(regs, at),
        nomics_mem::diag::describe_near(ptr, 8),
    );
    nomics_mem::diag::fatal(&report);
}

/// Check the regions of the read this thread extends ([`at_read`]); `seq`
/// is its sequence.
pub fn check_read(stage: &str, idx: &Index, seq: &[u8], regs: &Vec<Region>) {
    if !on() {
        return;
    }
    if let Some((at, why)) = first_fault(idx, seq.len(), regs) {
        let (addr, name) = READ.with(|c| {
            let c = c.borrow();
            (c.0, String::from_utf8_lossy(&c.1).into_owned())
        });
        let read = if addr == 0 {
            "read unknown (no at_read on this thread)".to_owned()
        } else {
            format!("read {name}: {}", batch_of(addr))
        };
        fail(
            stage,
            why,
            idx,
            &read,
            seq.len(),
            regs,
            Some(regs.capacity()),
            at,
        );
    }
}

/// Check every read's regions of a batch (`regs[i]` belongs to `reads[i]`).
pub fn check_batch(
    stage: &str,
    idx: &Index,
    reads: &[Read],
    regs: &[Vec<Region>],
    first_pair: u64,
) {
    if !on() {
        return;
    }
    for (i, (r, a)) in reads.iter().zip(regs).enumerate() {
        if let Some((at, why)) = first_fault(idx, r.seq.len(), a) {
            let read = format!(
                "read {} (index {i} of the batch, pair id {}): {}",
                String::from_utf8_lossy(&r.name),
                first_pair + (i / 2) as u64,
                batch_of(r as *const Read as usize)
            );
            fail(
                stage,
                why,
                idx,
                &read,
                r.seq.len(),
                a,
                Some(a.capacity()),
                at,
            );
        }
    }
}

/// Check both ends of pair `id` during pairing.
pub fn check_pair(stage: &str, idx: &Index, id: u64, s: &[Read; 2], a: &[Vec<Region>; 2]) {
    if !on() {
        return;
    }
    for m in 0..2 {
        if let Some((at, why)) = first_fault(idx, s[m].seq.len(), &a[m]) {
            fail(
                stage,
                why,
                idx,
                &pair_context(id, s, m),
                s[m].seq.len(),
                &a[m],
                Some(a[m].capacity()),
                at,
            );
        }
    }
}

fn pair_context(id: u64, s: &[Read; 2], m: usize) -> String {
    format!(
        "read {} (mate {} of pair id {id}): {}",
        String::from_utf8_lossy(&s[m].name),
        m + 1,
        batches_with_pair(id)
    )
}

/// `mem_pair` found a primary region whose `rid` indexes no contig: report
/// it as [`check_pair`] would (diagnostics on or off) and abort.
pub fn bad_rid_in_mem_pair(
    idx: &Index,
    id: u64,
    s: &[Read; 2],
    a: &[Vec<Region>; 2],
    m: usize,
    at: usize,
) -> ! {
    let why = region_fault(idx, s[m].seq.len(), &a[m][at]).unwrap_or("rid outside 0..n_seqs");
    fail(
        "in mem_pair (primary region's contig lookup)",
        why,
        idx,
        &pair_context(id, s, m),
        s[m].seq.len(),
        &a[m],
        Some(a[m].capacity()),
        at,
    );
}

// ---- Batch range for replays ----

/// The reader's view of one chunk's batches under [`batch_range`]: which
/// are skipped, which aligned, where the chunk stops, and each batch's
/// ordinal and first pair id (skipped batches count in both).
pub struct BatchCutter {
    skip: usize,
    max: usize,
    ordinal: usize,
    next_pair: u64,
}

/// What to do with a batch just read ([`BatchCutter::cut`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cut {
    Skip,
    Align { first_pair: u64, ordinal: usize },
}

impl BatchCutter {
    pub fn new((skip, max): (usize, usize)) -> Self {
        Self {
            skip,
            max,
            ordinal: 0,
            next_pair: 0,
        }
    }

    /// Whether to read another batch (`false`: `NOMICS_ALIGN_MAX_BATCHES`
    /// reached; the chunk ends as at EOF).
    pub fn wants_more(&self) -> bool {
        self.ordinal < self.skip.saturating_add(self.max)
    }

    /// Account for a nonempty batch of `pairs` pairs just read.
    pub fn cut(&mut self, pairs: u64) -> Cut {
        let (ordinal, first_pair) = (self.ordinal, self.next_pair);
        self.ordinal += 1;
        self.next_pair += pairs;
        if ordinal < self.skip {
            Cut::Skip
        } else {
            Cut::Align {
                first_pair,
                ordinal,
            }
        }
    }
}

/// The consumer's check that each chunk's batches arrive in input order:
/// the expected next pair id per chunk, tracked apart from the count of
/// pairs aligned (they differ by the pairs of skipped batches).
pub struct PairCursor {
    skipping: bool,
    next: Vec<Option<u64>>,
}

impl PairCursor {
    pub fn new(chunks: usize, skip: usize) -> Self {
        Self {
            skipping: skip > 0,
            next: vec![None; chunks],
        }
    }

    /// Record chunk `chunk`'s next finished batch; `false` if it is out of
    /// order. A chunk's first batch starts at pair 0, or anywhere when
    /// leading batches were skipped; every later one where the last ended.
    pub fn accept(&mut self, chunk: usize, first_pair: u64, pairs: u64) -> bool {
        let ok = match self.next[chunk] {
            None => first_pair == 0 || self.skipping,
            Some(next) => first_pair == next,
        };
        self.next[chunk] = Some(first_pair + pairs);
        ok
    }
}

/// Batches per chunk to skip (`NOMICS_ALIGN_SKIP_BATCHES`, read but not
/// aligned; pair ids still count them, so later batches align exactly as in
/// a full run) and to align after those (`NOMICS_ALIGN_MAX_BATCHES`; the
/// chunk then ends as if at EOF). Defaults: 0 and unlimited.
pub fn batch_range() -> (usize, usize) {
    static RANGE: std::sync::OnceLock<(usize, usize)> = std::sync::OnceLock::new();
    *RANGE.get_or_init(|| {
        let get = |k: &str| std::env::var(k).ok().and_then(|v| v.trim().parse::<usize>().ok());
        let skip = get("NOMICS_ALIGN_SKIP_BATCHES").unwrap_or(0);
        let max = get("NOMICS_ALIGN_MAX_BATCHES").filter(|&n| n > 0).unwrap_or(usize::MAX);
        if skip > 0 || max < usize::MAX {
            let end = if max == usize::MAX {
                "the end".to_owned()
            } else {
                (skip + max).to_string()
            };
            eprintln!(
                "WARNING: aligning only batches {skip}..{end} (0-based) of each chunk (NOMICS_ALIGN_SKIP_BATCHES / NOMICS_ALIGN_MAX_BATCHES): the output is partial by design"
            );
        }
        (skip, max)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cutting and consuming with skipped and capped batches (run in the
    /// default debug test build, where the consumer's order check is a
    /// debug assertion): ids and ordinals keep counting skipped batches.
    #[test]
    fn batch_range_skips_and_caps_in_order() {
        let sizes = [10u64, 12, 9, 11, 8];
        for (skip, max) in [(0, usize::MAX), (2, usize::MAX), (1, 2), (0, 1), (7, 3)] {
            let mut cutter = BatchCutter::new((skip, max));
            let mut cursor = PairCursor::new(1, skip);
            let mut aligned = Vec::new();
            let mut emitted = 0u64;
            for &n in &sizes {
                if !cutter.wants_more() {
                    break;
                }
                if let Cut::Align { first_pair, ordinal } = cutter.cut(n) {
                    assert!(cursor.accept(0, first_pair, n), "skip {skip} max {max}");
                    aligned.push((ordinal, first_pair));
                    emitted += n;
                }
            }
            let expect: Vec<(usize, u64)> = (0..sizes.len())
                .skip(skip)
                .take(max)
                .map(|i| (i, sizes[..i].iter().sum()))
                .collect();
            assert_eq!(aligned, expect, "skip {skip} max {max}");
            // With a skip, the first aligned id is past the emitted count
            // (what the old `first_pair == pairs[chunk]` assert compared).
            if skip > 0 && skip < sizes.len() {
                assert!(aligned[0].1 > 0 && emitted < aligned.last().unwrap().1 + 20);
            }
        }
        // Out of order is caught, and a chunk without a skip starts at 0.
        let mut cursor = PairCursor::new(2, 0);
        assert!(cursor.accept(0, 0, 10));
        assert!(!cursor.accept(0, 11, 5));
        assert!(!cursor.accept(1, 3, 5));
    }

    #[test]
    fn region_checker_flags_synthetic_bad_regions() {
        let good = Region {
            rid: 0,
            rb: 10,
            re: 60,
            qb: 0,
            qe: 50,
            ..Region::default()
        };
        assert_eq!(region_fault_in(1, 1000, 100, &good), None);
        // The reverse strand is in range, and qe == qb marks an excluded region.
        let rev = Region {
            rb: 1500,
            re: 1550,
            qb: 7,
            qe: 7,
            ..good
        };
        assert_eq!(region_fault_in(1, 1000, 100, &rev), None);
        let cases = [
            (Region { rid: -1, ..good }, "rid"),
            (Region { rid: 1, ..good }, "rid"),
            (Region { rb: -5, ..good }, "rb/re"),
            (Region { re: 2001, ..good }, "rb/re"),
            (Region { rb: 70, ..good }, "rb > re"),
            (Region { qb: -1, ..good }, "qb/qe"),
            (Region { qb: 60, ..good }, "qb/qe"),
            (Region { qe: 101, ..good }, "qb/qe"),
        ];
        for (bad, expect) in cases {
            let why =
                region_fault_in(1, 1000, 100, &bad).unwrap_or_else(|| panic!("{bad:?} passed"));
            assert!(why.contains(expect), "{bad:?}: {why}");
            assert_eq!(
                first_fault_in(1, 1000, 100, &[good, rev, bad, good]),
                Some((2, why))
            );
        }
        assert_eq!(first_fault_in(1, 1000, 100, &[good, rev]), None);
    }
}
