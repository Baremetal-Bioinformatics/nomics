//! Chain extension, redundancy removal and patching (`bwamem.c`
//! `mem_chain2aln`, `mem_sort_dedup_patch`, `mem_patch_reg`,
//! `mem_flt_chained_seeds`; `bwa.c` `bwa_gen_cigar2`; `bntseq.c`
//! `bns_get_seq`, `bns_fetch_seq`).
//!
//! Gate: the oracle's `R` dump (regions after `mem_sort_dedup_patch`).

use crate::chain::{Chain, Seed, mem_chain, mem_chain_flt};
use crate::index::Index;
use crate::ksort::introsort;
use crate::ksw;
use crate::seed::SmemAux;

/// bwa-mem scoring and band defaults (`mem_opt_init`).
pub const A: i32 = 1;
pub const B: i32 = 4;
pub const O_DEL: i32 = 6;
pub const E_DEL: i32 = 1;
pub const O_INS: i32 = 6;
pub const E_INS: i32 = 1;
pub const W: i32 = 100;
pub const ZDROP: i32 = 100;
pub const PEN_CLIP5: i32 = 5;
pub const PEN_CLIP3: i32 = 5;
const MAX_CHAIN_GAP: i64 = 10000;
const MASK_LEVEL_REDUN: f32 = 0.95;
pub const MAX_BAND_TRY: i32 = 2;
const PATCH_MAX_R_BW: f32 = 0.05;
const PATCH_MIN_SC_RATIO: f32 = 0.90;
const MEM_MINSC_COEF: f32 = 5.5;
const MEM_SEEDSW_COEF: f32 = 0.05;

/// `bwa_fill_scmat(1, 4, mat)`: match 1, mismatch -4, anything with N -1.
pub const MAT: [i8; 25] = {
    let mut m = [-1i8; 25];
    let mut i = 0;
    while i < 4 {
        let mut j = 0;
        while j < 4 {
            m[i * 5 + j] = if i == j { A as i8 } else { -(B as i8) };
            j += 1;
        }
        i += 1;
    }
    m
};

/// `mem_alnreg_t`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Region {
    pub rb: i64,
    pub re: i64,
    pub qb: i32,
    pub qe: i32,
    pub rid: i32,
    pub score: i32,
    pub truesc: i32,
    pub sub: i32,
    pub alt_sc: i32,
    pub csub: i32,
    pub sub_n: i32,
    pub w: i32,
    pub seedcov: i32,
    pub secondary: i32,
    pub secondary_all: i32,
    pub seedlen0: i32,
    pub n_comp: i32,
    pub is_alt: bool,
    pub frac_rep: f32,
    pub hash: u64,
}

impl Index {
    /// `bns_get_seq`: reference codes for `[beg, end)` in concatenated
    /// coordinates; empty if the range bridges the strand boundary.
    pub fn get_seq(&self, beg: i64, end: i64, out: &mut Vec<u8>) {
        get_seq_packed(self.pac.bytes(), self.bns.l_pac, beg, end, out);
    }

    /// `bns_fetch_seq`: clip `[beg, end)` to the contig containing `mid` and
    /// fetch it; returns the contig id.
    pub fn fetch_seq(&self, beg: &mut i64, mid: i64, end: &mut i64, out: &mut Vec<u8>) -> i32 {
        if *end < *beg {
            std::mem::swap(beg, end);
        }
        assert!(*beg <= mid && mid < *end);
        let (pos_f, is_rev) = self.bns.depos(mid);
        let rid = self.bns.pos2rid(pos_f);
        let ann = &self.bns.anns[rid as usize];
        let (mut far_beg, mut far_end) = (ann.offset, ann.offset + ann.len as i64);
        if is_rev {
            let tmp = far_beg;
            far_beg = (self.bns.l_pac << 1) - far_end;
            far_end = (self.bns.l_pac << 1) - tmp;
        }
        *beg = (*beg).max(far_beg);
        *end = (*end).min(far_end);
        self.get_seq(*beg, *end, out);
        assert_eq!(*end - *beg, out.len() as i64, "bns_fetch_seq length");
        rid
    }
}

/// `bns_get_seq` over the packed forward strand `pac` (`Index::get_seq`).
fn get_seq_packed(pac: &[u8], l_pac: i64, mut beg: i64, mut end: i64, out: &mut Vec<u8>) {
    out.clear();
    if end < beg {
        std::mem::swap(&mut beg, &mut end);
    }
    end = end.min(l_pac << 1);
    beg = beg.max(0);
    if beg >= l_pac || end <= l_pac {
        if beg >= l_pac {
            // bns_get_seq pushes `3 - pac[k]` for k = end_f down to
            // beg_f + 1: the forward bases (beg_f, end_f] reversed and
            // complemented (`3 - x == x ^ 3` on 2-bit codes).
            let beg_f = (l_pac << 1) - 1 - end;
            let end_f = (l_pac << 1) - 1 - beg;
            if end_f > beg_f {
                crate::index::unpack_2bit(pac, (beg_f + 1) as u64, (end_f + 1) as u64, out);
                out.reverse();
                for x in out.iter_mut() {
                    *x ^= 3;
                }
            }
        } else if end > beg {
            crate::index::unpack_2bit(pac, beg as u64, end as u64, out);
        }
    }
}

/// Global alignment between known end points (`bwa_gen_cigar2`).
pub struct GenCigar {
    pub score: i32,
    pub cigar: Vec<u32>,
    pub nm: i32,
    pub md: String,
}

/// `bwa_gen_cigar2` with bwa-mem's scoring. `None` where bwa returns without
/// setting the score (empty query/reference, strand-bridging, out of range).
pub fn gen_cigar2(
    idx: &Index,
    w_: i32,
    query: &[u8],
    rb: i64,
    re: i64,
    want_md: bool,
) -> Option<GenCigar> {
    GEN_SCRATCH.with_borrow_mut(|(q, rseq)| gen_cigar2_in(idx, w_, query, rb, re, want_md, q, rseq))
}

/// `gen_cigar2` with the query copy `q` and reference `rseq` in reused buffers.
#[allow(clippy::too_many_arguments)]
fn gen_cigar2_in(
    idx: &Index,
    w_: i32,
    query: &[u8],
    rb: i64,
    re: i64,
    want_md: bool,
    q: &mut Vec<u8>,
    rseq: &mut Vec<u8>,
) -> Option<GenCigar> {
    let l_pac = idx.bns.l_pac;
    let l_query = query.len() as i32;
    if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
        return None;
    }
    idx.get_seq(rb, re, rseq);
    let rlen = rseq.len() as i64;
    if re - rb != rlen {
        return None;
    }
    q.clear();
    q.extend_from_slice(query);
    if rb >= l_pac {
        // Reverse both so indels are placed leftmost on the forward strand.
        q.reverse();
        rseq.reverse();
    }
    // bwa's no-gap shortcut needs w_ == 0, which mem-mode never passes.
    let (score, cigar) = if l_query as i64 == re - rb && w_ == 0 {
        let score = (0..l_query as usize)
            .map(|i| MAT[rseq[i] as usize * 5 + q[i] as usize] as i32)
            .sum();
        (score, vec![(l_query as u32) << 4])
    } else {
        let mat0 = MAT[0] as i32;
        let half = (l_query + 1) >> 1;
        let max_ins = ((half * mat0 - O_INS) as f64 / E_INS as f64 + 1.0) as i32;
        let max_del = ((half * mat0 - O_DEL) as f64 / E_DEL as f64 + 1.0) as i32;
        let max_gap = max_ins.max(max_del).max(1);
        let diff = (rlen as i32 - l_query).abs();
        let w = ((max_gap + diff + 1) >> 1).min(w_).max(diff + 3);
        ksw::global2(
            l_query,
            q,
            rlen as i32,
            rseq,
            5,
            &MAT,
            O_DEL,
            E_DEL,
            O_INS,
            E_INS,
            w,
        )
    };
    let (mut nm, mut md) = (-1, String::new());
    if want_md {
        let int2base: &[u8; 5] = if rb < l_pac { b"ACGTN" } else { b"TGCAN" };
        let (mut x, mut y, mut u, mut n_mm, mut n_gap) = (0usize, 0usize, 0, 0, 0);
        let n = cigar.len();
        for (k, &c) in cigar.iter().enumerate() {
            let (op, len) = (c & 0xf, (c >> 4) as usize);
            match op {
                0 => {
                    for i in 0..len {
                        if q[x + i] != rseq[y + i] {
                            md.push_str(&u.to_string());
                            md.push(int2base[rseq[y + i] as usize] as char);
                            n_mm += 1;
                            u = 0;
                        } else {
                            u += 1;
                        }
                    }
                    x += len;
                    y += len;
                }
                2 => {
                    if k > 0 && k < n - 1 {
                        md.push_str(&u.to_string());
                        md.push('^');
                        md.extend((0..len).map(|i| int2base[rseq[y + i] as usize] as char));
                        u = 0;
                        n_gap += len as i32;
                    }
                    y += len;
                }
                1 => {
                    x += len;
                    n_gap += len as i32;
                }
                _ => {}
            }
        }
        md.push_str(&u.to_string());
        nm = n_mm + n_gap;
    }
    Some(GenCigar {
        score,
        cigar,
        nm,
        md,
    })
}

thread_local! {
    /// `gen_cigar2`'s and `gen_cigar2_score`'s query copy and reference.
    static GEN_SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<u8>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

/// `bwa_gen_cigar2(.., &score, 0, 0)` as `mem_patch_reg` calls it: the score
/// only, with BWA's no-CIGAR `ksw_global2` (`ksw::global2_score`) and
/// per-thread buffers. Same `None` cases as `gen_cigar2`.
fn gen_cigar2_score(idx: &Index, w_: i32, query: &[u8], rb: i64, re: i64) -> Option<i32> {
    let l_pac = idx.bns.l_pac;
    let l_query = query.len() as i32;
    if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
        return None;
    }
    GEN_SCRATCH.with_borrow_mut(|(q, rseq)| {
        idx.get_seq(rb, re, rseq);
        let rlen = rseq.len() as i64;
        if re - rb != rlen {
            return None;
        }
        q.clear();
        q.extend_from_slice(query);
        if rb >= l_pac {
            q.reverse();
            rseq.reverse();
        }
        Some(if l_query as i64 == re - rb && w_ == 0 {
            (0..l_query as usize)
                .map(|i| MAT[rseq[i] as usize * 5 + q[i] as usize] as i32)
                .sum()
        } else {
            let mat0 = MAT[0] as i32;
            let half = (l_query + 1) >> 1;
            let max_ins = ((half * mat0 - O_INS) as f64 / E_INS as f64 + 1.0) as i32;
            let max_del = ((half * mat0 - O_DEL) as f64 / E_DEL as f64 + 1.0) as i32;
            let max_gap = max_ins.max(max_del).max(1);
            let diff = (rlen as i32 - l_query).abs();
            let w = ((max_gap + diff + 1) >> 1).min(w_).max(diff + 3);
            ksw::global2_score(
                l_query,
                q,
                rlen as i32,
                rseq,
                5,
                &MAT,
                O_DEL,
                E_DEL,
                O_INS,
                E_INS,
                w,
            )
        })
    })
}

/// `cal_max_gap`.
#[inline]
fn cal_max_gap(qlen: i32) -> i32 {
    let l_del = ((qlen * A - O_DEL) as f64 / E_DEL as f64 + 1.0) as i32;
    let l_ins = ((qlen * A - O_INS) as f64 / E_INS as f64 + 1.0) as i32;
    l_del.max(l_ins).clamp(1, W << 1)
}

/// `mem_flt_chained_seeds`: a no-op unless `5.5 * ln(len) <= 0.05 * len`
/// (reads of roughly 700 bp and longer), where bwa re-scores seeds with
/// `ksw_align2`, which is not ported yet.
fn flt_chained_seeds(l_query: i32) {
    let min_l = (MEM_MINSC_COEF as f64) * crate::glibc_math::log(l_query as f64);
    assert!(
        min_l > (MEM_SEEDSW_COEF * l_query as f32) as f64,
        "reads of {l_query} bp need mem_seed_sw re-scoring, which is not implemented"
    );
}

/// Per-thread buffers.
#[derive(Default)]
pub struct AlignAux {
    pub smem: SmemAux,
    rseq: Vec<u8>,
    qs: Vec<u8>,
    rs: Vec<u8>,
}

/// Extension workload counters (NOMICS_TIMING), for sizing a GPU port:
/// chains, seeds, seeds extended, left/right `ksw_extend2` calls including
/// band retries, and query x target cells those calls cover.
pub mod ext_stats {
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
    pub static ON: AtomicBool = AtomicBool::new(false);
    pub static CHAINS: AtomicU64 = AtomicU64::new(0);
    pub static SEEDS: AtomicU64 = AtomicU64::new(0);
    pub static EXTENDED: AtomicU64 = AtomicU64::new(0);
    pub static LEFT: AtomicU64 = AtomicU64::new(0);
    pub static RIGHT: AtomicU64 = AtomicU64::new(0);
    pub static RETRIES: AtomicU64 = AtomicU64::new(0);
    pub static CELLS: AtomicU64 = AtomicU64::new(0);
    /// Extensions found precomputed, and ones recomputed on the CPU.
    pub static HITS: AtomicU64 = AtomicU64::new(0);
    /// Nanoseconds in chaining (SA lookups + chain building + filter), in
    /// chain2aln overall, and inside its ksw_extend2 calls.
    pub static CHAIN_NS: AtomicU64 = AtomicU64::new(0);
    pub static C2A_NS: AtomicU64 = AtomicU64::new(0);
    pub static KSW_NS: AtomicU64 = AtomicU64::new(0);
    pub fn on() -> bool {
        ON.load(Relaxed)
    }
    pub static MISSES: AtomicU64 = AtomicU64::new(0);
    /// Fused path: jobs launched, seeds the planner left to the CPU (a left
    /// or right query over `EXTEND_MAX_QUERY`), and band-1 calls the engine
    /// made (whether or not the replay consumed them).
    pub static FUSED_JOBS: AtomicU64 = AtomicU64::new(0);
    pub static FALLBACKS: AtomicU64 = AtomicU64::new(0);
    pub static FUSED_RETRIES: AtomicU64 = AtomicU64::new(0);
    pub fn add(c: &AtomicU64, n: u64) {
        if ON.load(Relaxed) {
            c.fetch_add(n, Relaxed);
        }
    }
    /// `sort_dedup_patch`: calls (n > 1), regions in, `patch_reg` calls,
    /// ones that reached the global alignment, its query x target cells, and
    /// patches made.
    pub static MERGE_CALLS: AtomicU64 = AtomicU64::new(0);
    pub static MERGE_REGS: AtomicU64 = AtomicU64::new(0);
    pub static PATCH_TRIES: AtomicU64 = AtomicU64::new(0);
    pub static PATCH_DP: AtomicU64 = AtomicU64::new(0);
    pub static PATCH_CELLS: AtomicU64 = AtomicU64::new(0);
    pub static PATCH_OK: AtomicU64 = AtomicU64::new(0);
    /// Region sorts whose keys repeat (sorted with the klib port).
    pub static MERGE_TIES: AtomicU64 = AtomicU64::new(0);
    pub fn merge_report() -> String {
        let g = |c: &AtomicU64| c.load(Relaxed);
        format!(
            "merge: {} calls, {} regions, {} sorts with tied keys, {} patch_reg calls, {} global alignments ({} cells), {} patches",
            g(&MERGE_CALLS),
            g(&MERGE_REGS),
            g(&MERGE_TIES),
            g(&PATCH_TRIES),
            g(&PATCH_DP),
            g(&PATCH_CELLS),
            g(&PATCH_OK)
        )
    }
    pub fn report(reads: u64) -> String {
        let g = |c: &AtomicU64| c.load(Relaxed) as f64;
        let r = reads.max(1) as f64;
        format!(
            "extension per record: {:.2} chains, {:.2} seeds, {:.2} extended ({:.0}% of seeds), {:.2} left + {:.2} right calls, {:.3} band retries, {:.0} cells; CPU per record: chaining {:.1} us, chain2aln {:.1} us (ksw_extend2 {:.1} us); precomputed lookups: {} hits, {} misses; fused: {} jobs, {} CPU fallback seeds, {} engine band retries",
            g(&CHAINS) / r,
            g(&SEEDS) / r,
            g(&EXTENDED) / r,
            100.0 * g(&EXTENDED) / g(&SEEDS).max(1.0),
            g(&LEFT) / r,
            g(&RIGHT) / r,
            g(&RETRIES) / r,
            g(&CELLS) / r,
            g(&CHAIN_NS) / r / 1e3,
            g(&C2A_NS) / r / 1e3,
            g(&KSW_NS) / r / 1e3,
            HITS.load(Relaxed),
            MISSES.load(Relaxed),
            FUSED_JOBS.load(Relaxed),
            FALLBACKS.load(Relaxed),
            FUSED_RETRIES.load(Relaxed)
        )
    }
}

/// `ksw_extend2` results computed ahead (on the GPU) for one read, keyed by
/// chain, seed, side, band and h0. A miss is computed on the CPU, so the set
/// only has to be likely, not complete.
///
/// Keys of each side come in ascending (chain, seed) order, as the planners
/// emit them; a lookup goes straight to its seed's keys and still compares
/// the whole [`ExtKey`].
#[derive(Clone, Default)]
pub struct Extended<'a> {
    /// Left extensions, then right extensions.
    keys: [&'a [ExtKey]; 2],
    results: [&'a [ksw::Extend]; 2],
    /// Per chain, the read-wide index of its first seed.
    seed0: Vec<u32>,
    /// Per read-wide seed, the first of its keys on each side (or `NONE`).
    first: Vec<[u32; 2]>,
    /// This read's lookups: hits, misses, and hits on band retries.
    lookups: std::cell::Cell<[u32; 3]>,
}

/// Identifies one `ksw_extend2` call of `chain2aln`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExtKey {
    pub chain: u32,
    pub seed: u32,
    pub right: bool,
    pub w: i32,
    pub h0: i32,
}

impl<'a> Extended<'a> {
    const NONE: u32 = u32::MAX;

    /// Index `keys` / `results` (left, then right) of the read whose chains
    /// are `chains`.
    pub fn new(chains: &[Chain], keys: [&'a [ExtKey]; 2], results: [&'a [ksw::Extend]; 2]) -> Self {
        let mut seed0 = Vec::with_capacity(chains.len());
        let mut n = 0u32;
        for c in chains {
            seed0.push(n);
            n += c.seeds.len() as u32;
        }
        let mut first = vec![[Self::NONE; 2]; n as usize];
        for side in 0..2 {
            assert_eq!(keys[side].len(), results[side].len());
            for (i, k) in keys[side].iter().enumerate() {
                debug_assert!(
                    i == 0
                        || (keys[side][i - 1].chain, keys[side][i - 1].seed) <= (k.chain, k.seed),
                    "extension keys out of (chain, seed) order"
                );
                debug_assert_eq!(k.right, side == 1);
                let c = k.chain as usize;
                assert!((k.seed as usize) < chains[c].seeds.len());
                let slot = &mut first[(seed0[c] + k.seed) as usize][side];
                if *slot == Self::NONE {
                    *slot = i as u32;
                }
            }
        }
        Self {
            keys,
            results,
            seed0,
            first,
            lookups: Default::default(),
        }
    }

    /// Lookups so far: hits, misses, and hits on band retries (`w != W`).
    pub fn lookups(&self) -> [u32; 3] {
        self.lookups.get()
    }

    fn get(&self, key: &ExtKey) -> Option<ksw::Extend> {
        let side = usize::from(key.right);
        let found = self
            .seed0
            .get(key.chain as usize)
            .and_then(|&s0| self.first.get((s0 + key.seed) as usize))
            .map(|f| f[side])
            .filter(|&at| at != Self::NONE)
            .and_then(|at| {
                let keys = &self.keys[side][at as usize..];
                keys.iter()
                    .take_while(|k| k.chain == key.chain && k.seed == key.seed)
                    .position(|k| k == key)
                    .map(|i| self.results[side][at as usize + i])
            });
        let mut n = self.lookups.get();
        n[usize::from(found.is_none())] += 1;
        n[2] += u32::from(found.is_some() && key.w != W);
        self.lookups.set(n);
        ext_stats::add(
            if found.is_some() {
                &ext_stats::HITS
            } else {
                &ext_stats::MISSES
            },
            1,
        );
        found
    }
}

/// One planned extension: its key and where its query then target bases
/// start in the planner's arena.
pub struct ExtPlan {
    pub key: ExtKey,
    pub at: usize,
    pub qlen: usize,
    pub tlen: usize,
    pub end_bonus: i32,
}

/// The reference window `chain2aln` extends chain `c` in; its bases go to
/// `rseq`.
fn chain_window(idx: &Index, query: &[u8], c: &Chain, rseq: &mut Vec<u8>) -> [i64; 2] {
    let seeds: &[Seed] = &c.seeds;
    let l_query = query.len() as i32;
    let l_pac = idx.bns.l_pac;
    let mut rmax = [l_pac << 1, 0i64];
    for t in seeds {
        let b = t.rbeg - (t.qbeg + cal_max_gap(t.qbeg)) as i64;
        let rest = l_query - t.qbeg - t.len;
        let e = t.rbeg + t.len as i64 + (rest + cal_max_gap(rest)) as i64;
        rmax[0] = rmax[0].min(b);
        rmax[1] = rmax[1].max(e);
    }
    rmax[0] = rmax[0].max(0);
    rmax[1] = rmax[1].min(l_pac << 1);
    if rmax[0] < l_pac && l_pac < rmax[1] {
        if seeds[0].rbeg < l_pac {
            rmax[1] = l_pac;
        } else {
            rmax[0] = l_pac;
        }
    }
    let [mut beg, mut end] = rmax;
    let rid = idx.fetch_seq(&mut beg, seeds[0].rbeg, &mut end, rseq);
    assert_eq!(c.rid, rid);
    [beg, end]
}

/// Plan the first-band left extension of every seed of every chain (the
/// ones `chain2aln` makes unless a seed is contained, or a band retry).
pub fn plan_left(
    idx: &Index,
    query: &[u8],
    chains: &[Chain],
    plans: &mut Vec<ExtPlan>,
    bases: &mut Vec<u8>,
    rseq: &mut Vec<u8>,
) {
    for (ci, c) in chains.iter().enumerate() {
        if c.seeds.is_empty() {
            continue;
        }
        let rmax = chain_window(idx, query, c, rseq);
        for (si, s) in c.seeds.iter().enumerate() {
            if s.qbeg <= 0 || s.qbeg as usize > nomics_gpu::EXTEND_MAX_QUERY {
                continue;
            }
            let tmp = (s.rbeg - rmax[0]) as usize;
            let at = bases.len();
            bases.extend(query[..s.qbeg as usize].iter().rev());
            bases.extend(rseq[..tmp].iter().rev());
            plans.push(ExtPlan {
                key: ExtKey {
                    chain: ci as u32,
                    seed: si as u32,
                    right: false,
                    w: W,
                    h0: s.len * A,
                },
                at,
                qlen: s.qbeg as usize,
                tlen: tmp,
                end_bonus: PEN_CLIP5,
            });
        }
    }
}

/// Plan the first-band right extension of every seed, with h0 from its left
/// extension in `left` (seeds whose left result is missing are left to the
/// CPU).
pub fn plan_right(
    idx: &Index,
    query: &[u8],
    chains: &[Chain],
    left: &Extended,
    plans: &mut Vec<ExtPlan>,
    bases: &mut Vec<u8>,
    rseq: &mut Vec<u8>,
) {
    let l_query = query.len() as i32;
    for (ci, c) in chains.iter().enumerate() {
        if c.seeds.is_empty() {
            continue;
        }
        let rmax = chain_window(idx, query, c, rseq);
        for (si, s) in c.seeds.iter().enumerate() {
            let s_end_q = s.qbeg + s.len;
            if s_end_q == l_query || (l_query - s_end_q) as usize > nomics_gpu::EXTEND_MAX_QUERY {
                continue;
            }
            let sc0 = if s.qbeg > 0 {
                let key = ExtKey {
                    chain: ci as u32,
                    seed: si as u32,
                    right: false,
                    w: W,
                    h0: s.len * A,
                };
                match left.get(&key) {
                    Some(x) => x.score,
                    None => continue,
                }
            } else {
                s.len * A
            };
            let re = (s.rbeg + s.len as i64 - rmax[0]) as usize;
            let tlen = (rmax[1] - rmax[0]) as usize - re;
            let at = bases.len();
            bases.extend_from_slice(&query[s_end_q as usize..]);
            bases.extend_from_slice(&rseq[re..re + tlen]);
            plans.push(ExtPlan {
                key: ExtKey {
                    chain: ci as u32,
                    seed: si as u32,
                    right: true,
                    w: W,
                    h0: sc0,
                },
                at,
                qlen: (l_query - s_end_q) as usize,
                tlen,
                end_bonus: PEN_CLIP3,
            });
        }
    }
}

/// Which sides of a seed (`qbeg`, `len` in a read of `l_query` bases) a
/// fused job extends, as `(left, right)`; `None` leaves the whole seed to the
/// CPU. A left query over `EXTEND_MAX_QUERY` sends the seed to the CPU,
/// because the right's h0 depends on it. A right query over it leaves a
/// left-only job, whose right then misses to the CPU. A seed with neither
/// side to extend gets no job.
pub fn fused_sides(qbeg: i32, len: i32, l_query: i32) -> Option<(bool, bool)> {
    let max = nomics_gpu::EXTEND_MAX_QUERY;
    if qbeg as usize > max {
        return None;
    }
    let rest = l_query - qbeg - len;
    let sides = (qbeg > 0, rest > 0 && rest as usize <= max);
    (sides.0 || sides.1).then_some(sides)
}

/// One fused job (both sides of seed `seed` of chain `chain`): its bases sit
/// in the planner's arena at `at` as left query, left target, right query,
/// right target. A zero query length means no call on that side.
pub struct FusedPlan {
    pub chain: u32,
    pub seed: u32,
    pub at: usize,
    pub lq: u32,
    pub lt: u32,
    pub rq: u32,
    pub rt: u32,
    pub h0: i32,
}

/// Plan the fused job of every seed of every chain (all the `ksw_extend2`
/// calls `chain2aln` makes for the seed, band retries included, unless it is
/// contained), in one pass with each chain's window computed once.
pub fn plan_fused(
    idx: &Index,
    query: &[u8],
    chains: &[Chain],
    plans: &mut Vec<FusedPlan>,
    bases: &mut Vec<u8>,
    rseq: &mut Vec<u8>,
) {
    let l_query = query.len() as i32;
    for (ci, c) in chains.iter().enumerate() {
        if c.seeds.is_empty() {
            continue;
        }
        let rmax = chain_window(idx, query, c, rseq);
        for (si, s) in c.seeds.iter().enumerate() {
            let s_end_q = s.qbeg + s.len;
            if s.qbeg as usize > nomics_gpu::EXTEND_MAX_QUERY
                || (l_query - s_end_q) as usize > nomics_gpu::EXTEND_MAX_QUERY
            {
                ext_stats::add(&ext_stats::FALLBACKS, 1);
            }
            let Some((left, right)) = fused_sides(s.qbeg, s.len, l_query) else {
                continue;
            };
            let at = bases.len();
            let (mut lt, mut rt) = (0, 0);
            if left {
                // As chain2aln: reversed prefix and reversed window before the seed.
                let tmp = (s.rbeg - rmax[0]) as usize;
                bases.extend(query[..s.qbeg as usize].iter().rev());
                bases.extend(rseq[..tmp].iter().rev());
                lt = tmp;
            }
            if right {
                let re = (s.rbeg + s.len as i64 - rmax[0]) as usize;
                let tlen = (rmax[1] - rmax[0]) as usize - re;
                bases.extend_from_slice(&query[s_end_q as usize..]);
                bases.extend_from_slice(&rseq[re..re + tlen]);
                rt = tlen;
            }
            plans.push(FusedPlan {
                chain: ci as u32,
                seed: si as u32,
                at,
                lq: if left { s.qbeg as u32 } else { 0 },
                lt: lt as u32,
                rq: if right { (l_query - s_end_q) as u32 } else { 0 },
                rt: rt as u32,
                h0: s.len * A,
            });
        }
    }
}

/// `mem_chain2aln`: extend the seeds of chain `c` (index `ci` of the read's
/// chains) into regions appended to `av`; extensions found in `pre` are used
/// as computed.
pub fn chain2aln(
    idx: &Index,
    query: &[u8],
    c: &Chain,
    av: &mut Vec<Region>,
    aux: &mut AlignAux,
    ci: u32,
    pre: &Extended,
) {
    let seeds: &[Seed] = &c.seeds;
    if seeds.is_empty() {
        return;
    }
    let l_query = query.len() as i32;
    ext_stats::add(&ext_stats::CHAINS, 1);
    ext_stats::add(&ext_stats::SEEDS, seeds.len() as u64);
    let rmax = chain_window(idx, query, c, &mut aux.rseq);
    let rseq = &aux.rseq;

    // Seeds in descending (score, index): ks_introsort_64 on distinct keys.
    let mut srt: Vec<u64> = seeds
        .iter()
        .enumerate()
        .map(|(i, s)| (s.score as u64) << 32 | i as u64)
        .collect();
    srt.sort_unstable();

    for k in (0..seeds.len()).rev() {
        let s = seeds[srt[k] as u32 as usize];
        let s_end_q = s.qbeg + s.len;
        let s_end_r = s.rbeg + s.len as i64;

        // Is the seed (almost) contained in an existing region?
        let mut contained = false;
        for p in av.iter() {
            if s.rbeg < p.rb || s_end_r > p.re || s.qbeg < p.qb || s_end_q > p.qe {
                continue;
            }
            if (s.len - p.seedlen0) as f64 > 0.1 * l_query as f64 {
                continue;
            }
            let qd = s.qbeg - p.qb;
            let rd = s.rbeg - p.rb;
            let max_gap = cal_max_gap(if (qd as i64) < rd { qd } else { rd as i32 });
            let w = max_gap.min(p.w);
            if (qd as i64 - rd) < w as i64 && (rd - qd as i64) < w as i64 {
                contained = true;
                break;
            }
            let qd = p.qe - s_end_q;
            let rd = p.re - s_end_r;
            let max_gap = cal_max_gap(if (qd as i64) < rd { qd } else { rd as i32 });
            let w = max_gap.min(p.w);
            if (qd as i64 - rd) < w as i64 && (rd - qd as i64) < w as i64 {
                contained = true;
                break;
            }
        }
        if contained {
            // Extend anyway if a long seed overlaps on a different diagonal.
            let mut overlapping = false;
            for &st in &srt[k + 1..] {
                if st == 0 {
                    continue;
                }
                let t = seeds[st as u32 as usize];
                if (t.len as f64) < s.len as f64 * 0.95 {
                    continue;
                }
                if s.qbeg <= t.qbeg
                    && s.qbeg + s.len - t.qbeg >= s.len >> 2
                    && (t.qbeg - s.qbeg) as i64 != t.rbeg - s.rbeg
                {
                    overlapping = true;
                    break;
                }
                if t.qbeg <= s.qbeg
                    && t.qbeg + t.len - s.qbeg >= s.len >> 2
                    && (s.qbeg - t.qbeg) as i64 != s.rbeg - t.rbeg
                {
                    overlapping = true;
                    break;
                }
            }
            if !overlapping {
                srt[k] = 0;
                continue;
            }
        }

        ext_stats::add(&ext_stats::EXTENDED, 1);
        let mut a = Region {
            w: W,
            score: -1,
            truesc: -1,
            rid: c.rid,
            ..Region::default()
        };
        let mut aw = [W, W];

        if s.qbeg > 0 {
            // Left extension on reversed query and reference.
            aux.qs.clear();
            aux.qs.extend(query[..s.qbeg as usize].iter().rev());
            let tmp = (s.rbeg - rmax[0]) as usize;
            aux.rs.clear();
            aux.rs.extend(rseq[..tmp].iter().rev());
            let mut r = None;
            for i in 0..MAX_BAND_TRY {
                let prev = a.score;
                aw[0] = W << i;
                ext_stats::add(
                    if i == 0 {
                        &ext_stats::LEFT
                    } else {
                        &ext_stats::RETRIES
                    },
                    1,
                );
                ext_stats::add(&ext_stats::CELLS, s.qbeg as u64 * tmp as u64);
                let key = ExtKey {
                    chain: ci,
                    seed: srt[k] as u32,
                    right: false,
                    w: aw[0],
                    h0: s.len * A,
                };
                let x = pre.get(&key).unwrap_or_else(|| {
                    let t = ext_stats::on().then(std::time::Instant::now);
                    let x = ksw::extend2(
                        s.qbeg,
                        &aux.qs,
                        tmp as i32,
                        &aux.rs,
                        5,
                        &MAT,
                        O_DEL,
                        E_DEL,
                        O_INS,
                        E_INS,
                        aw[0],
                        PEN_CLIP5,
                        ZDROP,
                        s.len * A,
                    );
                    if let Some(t) = t {
                        ext_stats::add(&ext_stats::KSW_NS, t.elapsed().as_nanos() as u64);
                    }
                    x
                });
                a.score = x.score;
                r = Some(x);
                if a.score == prev || x.max_off < (aw[0] >> 1) + (aw[0] >> 2) {
                    break;
                }
            }
            let x = r.unwrap();
            if x.gscore <= 0 || x.gscore <= a.score - PEN_CLIP5 {
                a.qb = s.qbeg - x.qle;
                a.rb = s.rbeg - x.tle as i64;
                a.truesc = a.score;
            } else {
                a.qb = 0;
                a.rb = s.rbeg - x.gtle as i64;
                a.truesc = x.gscore;
            }
        } else {
            a.score = s.len * A;
            a.truesc = a.score;
            a.qb = 0;
            a.rb = s.rbeg;
        }

        if s_end_q != l_query {
            // Right extension.
            let sc0 = a.score;
            let qe = s_end_q;
            let re = (s_end_r - rmax[0]) as usize;
            let tlen = (rmax[1] - rmax[0]) as usize - re;
            let mut r = None;
            for i in 0..MAX_BAND_TRY {
                let prev = a.score;
                aw[1] = W << i;
                ext_stats::add(
                    if i == 0 {
                        &ext_stats::RIGHT
                    } else {
                        &ext_stats::RETRIES
                    },
                    1,
                );
                ext_stats::add(&ext_stats::CELLS, (l_query - qe) as u64 * tlen as u64);
                let key = ExtKey {
                    chain: ci,
                    seed: srt[k] as u32,
                    right: true,
                    w: aw[1],
                    h0: sc0,
                };
                let x = pre.get(&key).unwrap_or_else(|| {
                    let t = ext_stats::on().then(std::time::Instant::now);
                    let x = ksw::extend2(
                        l_query - qe,
                        &query[qe as usize..],
                        tlen as i32,
                        &rseq[re..],
                        5,
                        &MAT,
                        O_DEL,
                        E_DEL,
                        O_INS,
                        E_INS,
                        aw[1],
                        PEN_CLIP3,
                        ZDROP,
                        sc0,
                    );
                    if let Some(t) = t {
                        ext_stats::add(&ext_stats::KSW_NS, t.elapsed().as_nanos() as u64);
                    }
                    x
                });
                a.score = x.score;
                r = Some(x);
                if a.score == prev || x.max_off < (aw[1] >> 1) + (aw[1] >> 2) {
                    break;
                }
            }
            let x = r.unwrap();
            if x.gscore <= 0 || x.gscore <= a.score - PEN_CLIP3 {
                a.qe = qe + x.qle;
                a.re = rmax[0] + re as i64 + x.tle as i64;
                a.truesc += a.score - sc0;
            } else {
                a.qe = l_query;
                a.re = rmax[0] + re as i64 + x.gtle as i64;
                a.truesc += x.gscore - sc0;
            }
        } else {
            a.qe = l_query;
            a.re = s_end_r;
        }

        a.seedcov = seeds
            .iter()
            .filter(|t| {
                t.qbeg >= a.qb
                    && t.qbeg + t.len <= a.qe
                    && t.rbeg >= a.rb
                    && t.rbeg + t.len as i64 <= a.re
            })
            .map(|t| t.len)
            .sum();
        a.w = aw[0].max(aw[1]);
        a.seedlen0 = s.len;
        a.frac_rep = c.frac_rep;
        av.push(a);
    }
}

/// `mem_patch_reg`: score of a global alignment joining `a` then `b` if it
/// is good enough to merge them, with the band used.
fn patch_reg(idx: &Index, query: &[u8], a: &Region, b: &Region) -> Option<(i32, i32)> {
    let l_pac = idx.bns.l_pac;
    ext_stats::add(&ext_stats::PATCH_TRIES, 1);
    debug_assert!(a.rid == b.rid && a.rb <= b.rb);
    if a.rb < l_pac && b.rb >= l_pac {
        return None;
    }
    if a.qb >= b.qb || a.qe >= b.qe || a.re >= b.re {
        return None;
    }
    let mut w = (((a.re - b.rb) - (a.qe - b.qb) as i64) as i32).abs();
    // bwamem.c:442-449 computes `r` first, but it is pure and only ever
    // or-ed after the `w` test, so it is computed only when `w` passes.
    let r = || {
        ((a.re - b.rb) as f64 / (b.re - a.rb) as f64 - (a.qe - b.qb) as f64 / (b.qe - a.qb) as f64)
            .abs()
    };
    if a.re < b.rb || a.qe < b.qb {
        if w > W << 1 || r() >= PATCH_MAX_R_BW as f64 {
            return None;
        }
    } else if w > W << 2 || r() >= (PATCH_MAX_R_BW * 2.0) as f64 {
        return None;
    }
    w += a.w + b.w;
    w = w.min(W << 2);
    ext_stats::add(&ext_stats::PATCH_DP, 1);
    ext_stats::add(
        &ext_stats::PATCH_CELLS,
        ((b.qe - a.qb) as i64 * (b.re - a.rb)) as u64,
    );
    let score = gen_cigar2_score(idx, w, &query[a.qb as usize..b.qe as usize], a.rb, b.re)
        .expect("mem_patch_reg: bwa_gen_cigar2 left the score unset");
    let sc = (b.score + a.score) as f64;
    // bwamem.c:455-456 compile to fmadd: (ratio * sc) + .499 in one rounding.
    let q_s =
        ((b.qe - a.qb) as f64 / ((b.qe - b.qb) + (a.qe - a.qb)) as f64).mul_add(sc, 0.499) as i32;
    let r_s =
        ((b.re - a.rb) as f64 / ((b.re - b.rb) + (a.re - a.rb)) as f64).mul_add(sc, 0.499) as i32;
    if (score as f64 / q_s.max(r_s) as f64) < PATCH_MIN_SC_RATIO as f64 {
        return None;
    }
    Some((score, w))
}

/// Per-thread buffers of `sort_dedup_patch`'s sorts.
#[derive(Default)]
struct MergeScratch {
    by_re: Vec<u64>,
    by_score: Vec<u128>,
    tmp: Vec<Region>,
}

thread_local! {
    static MERGE_SCRATCH: std::cell::Cell<MergeScratch> = const {
        std::cell::Cell::new(MergeScratch {
            by_re: Vec::new(),
            by_score: Vec::new(),
            tmp: Vec::new(),
        })
    };
}

/// `ks_introsort` of `a` under a comparator that compares one key, moving
/// packed `key << shift | index` tags instead of 96-byte regions. `tag`
/// packs region `i`; the caller guarantees the key fits above `shift` with
/// the same strict order as the comparator.
///
/// When all keys are distinct the sorted order is unique, so any correct
/// sort yields klib's permutation: the tags are sorted as integers with the
/// standard library's unstable sort. Only when some key repeats does klib's
/// order within the tie depend on its algorithm: then the tags are rebuilt
/// in input order and sorted with the literal `introsort` port, comparing
/// keys only. Its moves depend only on the comparator's outcomes, so the tag
/// permutation is the one `introsort(a, ..)` makes.
fn sort_regions<T>(
    a: &mut Vec<Region>,
    tags: &mut Vec<T>,
    tmp: &mut Vec<Region>,
    shift: u32,
    keep: impl Fn(&Region) -> bool,
    tag: impl Fn(&Region, u32) -> T,
    index: impl Fn(T) -> usize,
) -> bool
where
    T: Copy + Ord + std::ops::Shr<u32, Output = T>,
{
    let fill = |tags: &mut Vec<T>| {
        tags.clear();
        tags.extend(
            a.iter()
                .enumerate()
                .filter(|(_, r)| keep(r))
                .map(|(i, r)| tag(r, i as u32)),
        );
    };
    fill(tags);
    tags.sort_unstable();
    let ties = tags.windows(2).any(|t| t[0] >> shift == t[1] >> shift);
    if ties {
        ext_stats::add(&ext_stats::MERGE_TIES, 1);
        fill(tags);
        introsort(tags, |&x, &y| x >> shift < y >> shift);
    }
    tmp.clear();
    tmp.extend(tags.iter().map(|&t| a[index(t)]));
    std::mem::swap(a, tmp);
    ties
}

/// `ks_introsort(mem_ars2, ..)`: `re` ascending. Tags are `re << 31 | i`.
fn sort_by_re(a: &mut Vec<Region>, sc: &mut MergeScratch) {
    if a.len() < (1 << 31) && a.iter().all(|r| (0..1i64 << 33).contains(&r.re)) {
        sort_regions(
            a,
            &mut sc.by_re,
            &mut sc.tmp,
            31,
            |_| true,
            |r, i| (r.re as u64) << 31 | i as u64,
            |t| (t & ((1 << 31) - 1)) as usize,
        );
    } else {
        introsort(a, |x, y| x.re < y.re);
    }
}

/// bwamem.c:497-503: drop the excluded regions (`qe == qb`), then
/// `ks_introsort(mem_ars, ..)`: `score` descending, then `rb`, then `qb`
/// ascending. Tags (of the kept regions, in order) are
/// `!(score ^ 1 << 31) << 96 | rb << 62 | qb << 32 | i`, whose integer order
/// above bit 32 is that comparator's. Returns false only when no two kept
/// regions share (score, rb, qb), i.e. when bwa's "identical hits" pass
/// that follows can mark nothing.
fn filter_sort_by_score(a: &mut Vec<Region>, sc: &mut MergeScratch) -> bool {
    if a.iter()
        .all(|r| (0..1i64 << 34).contains(&r.rb) && (0..1 << 30).contains(&r.qb))
    {
        sort_regions(
            a,
            &mut sc.by_score,
            &mut sc.tmp,
            32,
            |r| r.qe > r.qb,
            |r, i| {
                let score = !(r.score as u32 ^ 0x8000_0000);
                (score as u128) << 96 | (r.rb as u128) << 62 | (r.qb as u128) << 32 | i as u128
            },
            |t| t as u32 as usize,
        )
    } else {
        a.retain(|r| r.qe > r.qb);
        introsort(a, |x, y| {
            x.score > y.score
                || (x.score == y.score && (x.rb < y.rb || (x.rb == y.rb && x.qb < y.qb)))
        });
        true
    }
}

/// `mem_sort_dedup_patch`. `query` is `None` after mate rescue, where bwa
/// passes no reference and so never patches.
pub fn sort_dedup_patch(idx: &Index, query: Option<&[u8]>, a: &mut Vec<Region>) {
    if a.len() <= 1 {
        return;
    }
    match merge_dump::sink() {
        None => sort_dedup_patch_core(idx, query, a),
        Some(sink) => {
            let b = merge_dump::begin(query, a);
            sort_dedup_patch_core(idx, query, a);
            merge_dump::finish(sink, b, a);
        }
    }
}

fn sort_dedup_patch_core(idx: &Index, query: Option<&[u8]>, a: &mut Vec<Region>) {
    let n = a.len();
    ext_stats::add(&ext_stats::MERGE_CALLS, 1);
    ext_stats::add(&ext_stats::MERGE_REGS, n as u64);
    let mut sc = MERGE_SCRATCH.take();
    sort_by_re(a, &mut sc);
    let (mut spans_ok, mut min_qe) = (true, i32::MAX);
    for r in a.iter_mut() {
        r.n_comp = 1;
        // Query coordinates far from overflow keep `w` below exact too.
        spans_ok &= r.re >= r.rb && (0..1 << 30).contains(&r.qb) && (0..1 << 30).contains(&r.qe);
        min_qe = min_qe.min(r.qe);
    }
    for i in 1..n {
        if a[i].rid != a[i - 1].rid || a[i].rb >= a[i - 1].re + MAX_CHAIN_GAP {
            continue;
        }
        // j counts down from i - 1 to 0, as `ju`.
        let mut j = i;
        while j > 0 {
            let ju = j - 1;
            let (p, q) = (&a[i], &a[ju]);
            if !(p.rid == q.rid && p.rb < q.re + MAX_CHAIN_GAP) {
                break;
            }
            // Stop early once no remaining j can change anything. `a` stays
            // sorted by `re` (nothing below changes `re`), so the remaining
            // q have `re <= q.re` and a reference gap `G = p.rb - re` at
            // least this one's. Then (bwamem.c:477-486, 447-449):
            // - not redundant: `or = q.re - p.rb <= 0`, and `0.95 * mr >= 0`
            //   since every region has `re >= rb` (`spans_ok`; a patch only
            //   lowers `p.rb` below `p.re`), so `or > 0.95 * mr` is false;
            // - not patched: with no query never; else `mem_patch_reg` needs
            //   `w = |G + (q.qe - p.qb)| <= opt->w << 2` on every path, and
            //   live regions keep their initial `qe` (only exclusions change
            //   it), so `G + min_qe - p.qb > w << 2` rules it out (all of it
            //   far inside i32, so bwa's int truncation of `w` is moot).
            // So p (and every q) is left as is and, by induction, so are the
            // later j: breaking here is what bwa's loop does, minus no-ops.
            if spans_ok
                && q.re <= p.rb
                && (query.is_none() || (p.rb - q.re) + (min_qe - p.qb) as i64 > (W << 2) as i64)
            {
                break;
            }
            j -= 1;
            if q.qe == q.qb {
                continue; // already excluded
            }
            let or = q.re - p.rb;
            let oq = if q.qb < p.qb {
                q.qe - p.qb
            } else {
                p.qe - q.qb
            } as i64;
            let mr = (q.re - q.rb).min(p.re - p.rb);
            let mq = ((q.qe - q.qb).min(p.qe - p.qb)) as i64;
            if or as f32 > MASK_LEVEL_REDUN * mr as f32 && oq as f32 > MASK_LEVEL_REDUN * mq as f32
            {
                if p.score < q.score {
                    a[i].qe = a[i].qb;
                    break;
                } else {
                    a[ju].qe = a[ju].qb;
                }
            } else if q.rb < p.rb
                && let Some(query) = query
                && let Some((score, w)) = patch_reg(idx, query, q, p)
                && score > 0
            {
                ext_stats::add(&ext_stats::PATCH_OK, 1);
                let q = *q;
                let pm = &mut a[i];
                pm.n_comp += q.n_comp + 1;
                pm.seedcov = pm.seedcov.max(q.seedcov);
                pm.sub = pm.sub.max(q.sub);
                pm.csub = pm.csub.max(q.csub);
                pm.qb = q.qb;
                pm.rb = q.rb;
                pm.truesc = score;
                pm.score = score;
                pm.w = w;
                a[ju].qb = a[ju].qe;
            }
        }
    }
    let ties = filter_sort_by_score(a, &mut sc);
    MERGE_SCRATCH.set(sc);
    if !ties {
        // Nothing to mark below, and every kept region has qe > qb.
        return;
    }
    for i in 1..a.len() {
        if a[i].score == a[i - 1].score && a[i].rb == a[i - 1].rb && a[i].qb == a[i - 1].qb {
            a[i].qe = a[i].qb;
        }
    }
    // bwa keeps a[0] unconditionally here (m starts at 1).
    let mut m = 1;
    for i in 1..a.len() {
        if a[i].qe > a[i].qb {
            a[m] = a[i];
            m += 1;
        }
    }
    a.truncate(m);
}

/// `mem_align1_core`: all regions of one read (codes 0..4).
pub fn align1_core(idx: &Index, seq: &[u8], aux: &mut AlignAux) -> Vec<Region> {
    let chains = filtered_chains(seq, mem_chain(idx, seq, &mut aux.smem));
    align1_extend(idx, seq, &chains, aux, &Extended::default())
}

/// `mem_align1_core` given the read's seed intervals (e.g. from the GPU).
pub fn align1_core_from(
    idx: &Index,
    seq: &[u8],
    mem: &[crate::seed::Intv],
    aux: &mut AlignAux,
) -> Vec<Region> {
    let chains = align1_chains(idx, seq, Some(mem), aux);
    align1_extend(idx, seq, &chains, aux, &Extended::default())
}

/// The first half of `mem_align1_core`: the read's filtered chains, from
/// seed intervals in `mem` or (`None`) seeded here.
pub fn align1_chains(
    idx: &Index,
    seq: &[u8],
    mem: Option<&[crate::seed::Intv]>,
    aux: &mut AlignAux,
) -> Vec<Chain> {
    let t = ext_stats::on().then(std::time::Instant::now);
    let chains = match mem {
        Some(mem) => crate::chain::mem_chain_from(idx, seq.len() as i32, mem),
        None => mem_chain(idx, seq, &mut aux.smem),
    };
    let chains = filtered_chains(seq, chains);
    if let Some(t) = t {
        ext_stats::add(&ext_stats::CHAIN_NS, t.elapsed().as_nanos() as u64);
    }
    chains
}

fn filtered_chains(seq: &[u8], mut chains: Vec<Chain>) -> Vec<Chain> {
    mem_chain_flt(&mut chains);
    if !chains.is_empty() {
        flt_chained_seeds(seq.len() as i32);
    }
    chains
}

/// The second half of `mem_align1_core`: extend the chains into regions,
/// using extensions computed ahead in `pre`.
pub fn align1_extend(
    idx: &Index,
    seq: &[u8],
    chains: &[Chain],
    aux: &mut AlignAux,
    pre: &Extended,
) -> Vec<Region> {
    let mut regs = Vec::new();
    let t = ext_stats::on().then(std::time::Instant::now);
    for (ci, c) in chains.iter().enumerate() {
        chain2aln(idx, seq, c, &mut regs, aux, ci as u32, pre);
    }
    if let Some(t) = t {
        ext_stats::add(&ext_stats::C2A_NS, t.elapsed().as_nanos() as u64);
    }
    crate::diag::check_read("after chain2aln (extension)", idx, seq, &regs);
    sort_dedup_patch(idx, Some(seq), &mut regs);
    crate::diag::check_read("after sort_dedup_patch (region merge)", idx, seq, &regs);
    for r in regs.iter_mut() {
        if r.rid >= 0 && idx.bns.anns[r.rid as usize].is_alt {
            r.is_alt = true;
        }
    }
    regs
}

/// `NOMICS_DUMP_MERGE=<dir>`: append every `sort_dedup_patch` call with two
/// or more regions (query, regions before, regions after; calls with fewer
/// return untouched and are not written) to `<dir>/merge.bin`, for the
/// offline `merge_bench` replay. Little-endian records:
/// `u8 has_query, u32 qlen, query, u32 n, n x region, u32 m, m x region`,
/// where a region is `i64 rb re, i32 qb qe rid score truesc sub alt_sc csub
/// sub_n w seedcov secondary secondary_all seedlen0 n_comp, u8 is_alt,
/// u32 frac_rep bits, u64 hash` (`REGION_BYTES`).
pub mod merge_dump {
    use super::Region;
    use std::fs::File;
    use std::io::Write as _;
    use std::sync::{Mutex, OnceLock};

    pub const REGION_BYTES: usize = 16 + 15 * 4 + 1 + 4 + 8;

    pub fn sink() -> Option<&'static Mutex<File>> {
        static SINK: OnceLock<Option<Mutex<File>>> = OnceLock::new();
        SINK.get_or_init(|| {
            let dir = std::path::PathBuf::from(std::env::var_os("NOMICS_DUMP_MERGE")?);
            Some(Mutex::new(
                File::create(dir.join("merge.bin")).expect("create merge dump"),
            ))
        })
        .as_ref()
    }

    pub fn put_region(b: &mut Vec<u8>, r: &Region) {
        b.extend_from_slice(&r.rb.to_le_bytes());
        b.extend_from_slice(&r.re.to_le_bytes());
        for x in [
            r.qb,
            r.qe,
            r.rid,
            r.score,
            r.truesc,
            r.sub,
            r.alt_sc,
            r.csub,
            r.sub_n,
            r.w,
            r.seedcov,
            r.secondary,
            r.secondary_all,
            r.seedlen0,
            r.n_comp,
        ] {
            b.extend_from_slice(&x.to_le_bytes());
        }
        b.push(r.is_alt as u8);
        b.extend_from_slice(&r.frac_rep.to_bits().to_le_bytes());
        b.extend_from_slice(&r.hash.to_le_bytes());
    }

    pub fn get_region(b: &[u8]) -> Region {
        assert!(b.len() >= REGION_BYTES);
        let i64_at = |k: usize| i64::from_le_bytes(b[k..k + 8].try_into().unwrap());
        let i32_at = |k: usize| i32::from_le_bytes(b[k..k + 4].try_into().unwrap());
        let f = |n: usize| i32_at(16 + 4 * n);
        Region {
            rb: i64_at(0),
            re: i64_at(8),
            qb: f(0),
            qe: f(1),
            rid: f(2),
            score: f(3),
            truesc: f(4),
            sub: f(5),
            alt_sc: f(6),
            csub: f(7),
            sub_n: f(8),
            w: f(9),
            seedcov: f(10),
            secondary: f(11),
            secondary_all: f(12),
            seedlen0: f(13),
            n_comp: f(14),
            is_alt: b[76] != 0,
            frac_rep: f32::from_bits(u32::from_le_bytes(b[77..81].try_into().unwrap())),
            hash: u64::from_le_bytes(b[81..89].try_into().unwrap()),
        }
    }

    /// Region equality with floats compared bitwise.
    pub fn same(a: &Region, b: &Region) -> bool {
        let (mut x, mut y) = (Vec::new(), Vec::new());
        put_region(&mut x, a);
        put_region(&mut y, b);
        x == y
    }

    pub fn regions(b: &mut Vec<u8>, a: &[Region]) {
        b.extend_from_slice(&(a.len() as u32).to_le_bytes());
        for r in a {
            put_region(b, r);
        }
    }

    pub fn begin(query: Option<&[u8]>, a: &[Region]) -> Vec<u8> {
        let mut b = vec![query.is_some() as u8];
        let q = query.unwrap_or(&[]);
        b.extend_from_slice(&(q.len() as u32).to_le_bytes());
        b.extend_from_slice(q);
        regions(&mut b, a);
        b
    }

    pub fn finish(sink: &Mutex<File>, mut b: Vec<u8>, a: &[Region]) {
        regions(&mut b, a);
        sink.lock()
            .unwrap()
            .write_all(&b)
            .expect("write merge dump");
    }

    /// One dumped call.
    pub struct Record {
        pub query: Option<Vec<u8>>,
        pub before: Vec<Region>,
        pub after: Vec<Region>,
    }

    /// The next record of a dump, or `None` at its end.
    pub fn read_record(r: &mut impl std::io::Read) -> std::io::Result<Option<Record>> {
        let mut has = [0u8; 1];
        if r.read(&mut has)? == 0 {
            return Ok(None);
        }
        let u32_of = |r: &mut dyn std::io::Read| -> std::io::Result<usize> {
            let mut b = [0u8; 4];
            r.read_exact(&mut b)?;
            Ok(u32::from_le_bytes(b) as usize)
        };
        let ql = u32_of(r)?;
        let mut q = vec![0u8; ql];
        r.read_exact(&mut q)?;
        let regs = |r: &mut dyn std::io::Read| -> std::io::Result<Vec<Region>> {
            let n = u32_of(r)?;
            let mut b = vec![0u8; n * REGION_BYTES];
            r.read_exact(&mut b)?;
            Ok(b.as_chunks::<REGION_BYTES>()
                .0
                .iter()
                .map(|c| get_region(c))
                .collect())
        };
        let before = regs(r)?;
        let after = regs(r)?;
        Ok(Some(Record {
            query: (has[0] != 0).then_some(q),
            before,
            after,
        }))
    }

    pub fn read_all(mut bytes: &[u8]) -> Vec<Record> {
        let mut out = Vec::new();
        while let Some(rec) = read_record(&mut bytes).expect("truncated merge dump") {
            out.push(rec);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain(seeds: usize) -> Chain {
        Chain {
            pos: 0,
            rid: 0,
            first: -1,
            w: 0,
            kept: 0,
            is_alt: false,
            frac_rep: 0.0,
            seeds: vec![
                Seed {
                    rbeg: 0,
                    qbeg: 0,
                    len: 19,
                    score: 19,
                };
                seeds
            ],
        }
    }

    fn key(chain: u32, seed: u32, right: bool, w: i32, h0: i32) -> ExtKey {
        ExtKey {
            chain,
            seed,
            right,
            w,
            h0,
        }
    }

    fn ext(score: i32) -> ksw::Extend {
        ksw::Extend {
            score,
            qle: 0,
            tle: 0,
            gtle: 0,
            gscore: 0,
            max_off: 0,
        }
    }

    /// The tag sorts give `introsort`'s exact permutation on whole regions,
    /// ties included (few distinct keys force the klib fallback and its
    /// combsort), and the score sort drops excluded regions first.
    #[test]
    fn region_sorts_match_introsort() {
        let mut x = 0x5851_f42d_4c95_7f2du64;
        let mut next = move |n: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x % n
        };
        let mut sc = MergeScratch::default();
        for case in 0..3000 {
            let n = 1 + next(if case % 10 == 0 { 600 } else { 40 }) as usize;
            let range = [3u64, 50, 1 << 33][case % 3];
            let a: Vec<Region> = (0..n)
                .map(|i| {
                    let rb = next(range) as i64;
                    let qb = next(if case % 2 == 0 { 3 } else { 150 }) as i32;
                    Region {
                        rb,
                        re: rb + next(range.min(300)) as i64,
                        qb,
                        qe: qb + next(3) as i32,
                        score: next(if case % 4 == 0 { 3 } else { 200 }) as i32 - 50,
                        hash: i as u64,
                        ..Region::default()
                    }
                })
                .collect();
            let mut want = a.clone();
            introsort(&mut want, |x, y| x.re < y.re);
            let mut got = a.clone();
            sort_by_re(&mut got, &mut sc);
            assert_eq!(got, want, "re case {case}");

            let mut want = a.clone();
            want.retain(|r| r.qe > r.qb);
            introsort(&mut want, |x, y| {
                x.score > y.score
                    || (x.score == y.score && (x.rb < y.rb || (x.rb == y.rb && x.qb < y.qb)))
            });
            let mut got = a.clone();
            let ties = filter_sort_by_score(&mut got, &mut sc);
            assert_eq!(got, want, "score case {case}");
            let dup = want
                .windows(2)
                .any(|w| (w[0].score, w[0].rb, w[0].qb) == (w[1].score, w[1].rb, w[1].qb));
            assert!(ties || !dup, "score case {case}: ties missed");
        }
    }

    /// `get_seq_packed` against a literal per-base port of `bns_get_seq`,
    /// on ranges inside, across and beyond both strands, swapped and empty.
    #[test]
    fn get_seq_matches_per_base() {
        fn reference(pac: &[u8], l_pac: i64, mut beg: i64, mut end: i64) -> Vec<u8> {
            let get = |k: i64| pac[(k >> 2) as usize] >> ((!k & 3) << 1) & 3;
            let mut out = Vec::new();
            if end < beg {
                std::mem::swap(&mut beg, &mut end);
            }
            end = end.min(l_pac << 1);
            beg = beg.max(0);
            if beg >= l_pac || end <= l_pac {
                if beg >= l_pac {
                    let beg_f = (l_pac << 1) - 1 - end;
                    let end_f = (l_pac << 1) - 1 - beg;
                    let mut k = end_f;
                    while k > beg_f {
                        out.push(3 - get(k));
                        k -= 1;
                    }
                } else {
                    out.extend((beg..end).map(get));
                }
            }
            out
        }
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |n: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x % n
        };
        let mut out = vec![7u8; 3];
        for l_pac in [1i64, 2, 3, 4, 5, 7, 8, 9, 31, 64, 133] {
            let pac: Vec<u8> = (0..l_pac / 4 + 1).map(|_| next(256) as u8).collect();
            let span = (l_pac << 1) + 8;
            for _ in 0..2000 {
                let beg = next(span as u64) as i64 - 4;
                let end = next(span as u64) as i64 - 4;
                get_seq_packed(&pac, l_pac, beg, end, &mut out);
                assert_eq!(
                    out,
                    reference(&pac, l_pac, beg, end),
                    "l_pac={l_pac} [{beg},{end})"
                );
            }
        }
    }

    /// The per-seed index finds exactly the keys a linear scan would, and
    /// still compares band and h0.
    #[test]
    fn extended_lookup_matches_linear_scan() {
        let chains = [chain(3), chain(0), chain(2)];
        let left = [
            key(0, 0, false, 100, 19),
            key(0, 0, false, 200, 19),
            key(0, 2, false, 100, 19),
            key(2, 1, false, 100, 30),
        ];
        let right = [
            key(0, 1, true, 100, 19),
            key(2, 0, true, 100, 44),
            key(2, 0, true, 200, 44),
            key(2, 1, true, 100, 50),
        ];
        let lr: Vec<_> = (0..4).map(ext).collect();
        let rr: Vec<_> = (10..14).map(ext).collect();
        let pre = Extended::new(&chains, [&left, &right], [&lr, &rr]);
        let mut probes = Vec::new();
        for c in 0..4 {
            for s in 0..4 {
                for right in [false, true] {
                    for w in [100, 200, 400] {
                        for h0 in [19, 30, 44, 50, 51] {
                            probes.push(key(c, s, right, w, h0));
                        }
                    }
                }
            }
        }
        let mut hits = 0;
        for k in &probes {
            let side = usize::from(k.right);
            let (keys, res) = if k.right { (&right, &rr) } else { (&left, &lr) };
            let want = keys.iter().position(|x| x == k).map(|i| res[i]);
            assert_eq!(pre.get(k), want, "{k:?}");
            hits += usize::from(want.is_some());
            let _ = side;
        }
        assert_eq!(hits, 8);
        let none = Extended::default();
        assert!(probes.iter().all(|k| none.get(k).is_none()));
    }
}
