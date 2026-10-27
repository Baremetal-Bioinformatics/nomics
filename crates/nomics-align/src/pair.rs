//! Paired-end resolution (`bwamem_pair.c` `mem_pestat`, `mem_matesw`,
//! `mem_pair`, `mem_sam_pe`; `bwamem.c` `mem_mark_primary_se`).
//!
//! Floating-point expressions follow the golden aarch64 object code, which
//! fuses several multiply-adds (see docs/spec/bwa.md section 10).

use crate::extend::{A, B, E_DEL, E_INS, MAT, O_DEL, O_INS, Region, sort_dedup_patch};
use crate::index::Index;
use crate::ksort::introsort;
use crate::ksw::{self, KSW_XBYTE, KSW_XSTART, KSW_XSUBO};
use crate::libm::erfc;
use crate::sam::{Aln, AlnSink, BamRecords, Read, T, approx_mapq_se, gen_alt, reg2aln, reg2out};

const MIN_SEED_LEN: i32 = 19;
const MASK_LEVEL: f32 = 0.50;
const MAX_INS: i64 = 10000;
const PEN_UNPAIRED: i32 = 17;
const MAX_MATESW: usize = 50;
const MIN_RATIO: f64 = 0.8;
const MIN_DIR_CNT: usize = 10;
const MIN_DIR_RATIO: f64 = 0.05;
const OUTLIER_BOUND: f64 = 2.0;
const MAPPING_BOUND: f64 = 3.0;
const MAX_STDDEV: f64 = 4.0;

/// `mem_pestat_t` for one orientation (0 FF, 1 FR, 2 RF, 3 RR).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PeStat {
    pub low: i32,
    pub high: i32,
    pub failed: bool,
    pub avg: f64,
    pub std: f64,
}

/// `hash_64` (utils.h).
#[inline]
pub fn hash_64(mut key: u64) -> u64 {
    key = key.wrapping_add(!(key << 32));
    key ^= key >> 22;
    key = key.wrapping_add(!(key << 13));
    key ^= key >> 8;
    key = key.wrapping_add(key << 3);
    key ^= key >> 15;
    key = key.wrapping_add(!(key << 27));
    key ^= key >> 31;
    key
}

/// `mem_infer_dir`: orientation and distance of two region starts.
#[inline]
fn infer_dir(l_pac: i64, b1: i64, b2: i64) -> (usize, i64) {
    let (r1, r2) = (b1 >= l_pac, b2 >= l_pac);
    let p2 = if r1 == r2 { b2 } else { (l_pac << 1) - 1 - b2 };
    let dist = if p2 > b1 { p2 - b1 } else { b1 - p2 };
    (((r1 != r2) as usize) ^ if p2 > b1 { 0 } else { 3 }, dist)
}

/// `cal_sub`.
fn cal_sub(r: &[Region]) -> i32 {
    for j in 1..r.len() {
        let b_max = r[j].qb.max(r[0].qb);
        let e_min = r[j].qe.min(r[0].qe);
        if e_min > b_max {
            let min_l = (r[j].qe - r[j].qb).min(r[0].qe - r[0].qb);
            if (e_min - b_max) as f32 >= min_l as f32 * MASK_LEVEL {
                return r[j].score;
            }
        }
    }
    MIN_SEED_LEN * A
}

/// `mem_pestat`: insert-size statistics of one batch from each pair's
/// regions (as returned by `align1_core`). Returns the stats and the log
/// lines bwa prints, for comparison with the golden's `.command.err`.
pub fn pestat(l_pac: i64, regs: &[Vec<Region>]) -> ([PeStat; 4], Vec<String>) {
    let mut pes = [PeStat::default(); 4];
    let mut log = Vec::new();
    let mut isize: [Vec<u64>; 4] = Default::default();
    for pair in regs.chunks_exact(2) {
        let (r0, r1) = (&pair[0], &pair[1]);
        if r0.is_empty() || r1.is_empty() {
            continue;
        }
        if cal_sub(r0) as f64 > MIN_RATIO * r0[0].score as f64 {
            continue;
        }
        if cal_sub(r1) as f64 > MIN_RATIO * r1[0].score as f64 {
            continue;
        }
        if r0[0].rid != r1[0].rid {
            continue;
        }
        let (dir, is) = infer_dir(l_pac, r0[0].rb, r1[0].rb);
        if is != 0 && is <= MAX_INS {
            isize[dir].push(is as u64);
        }
    }
    log.push(format!(
        "[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): ({}, {}, {}, {})",
        isize[0].len(),
        isize[1].len(),
        isize[2].len(),
        isize[3].len()
    ));
    let lab = |d: usize| format!("{}{}", ['F', 'R'][d >> 1 & 1], ['F', 'R'][d & 1]);
    for d in 0..4 {
        let r = &mut pes[d];
        let q = &mut isize[d];
        if q.len() < MIN_DIR_CNT {
            log.push(format!(
                "[M::mem_pestat] skip orientation {} as there are not enough pairs",
                lab(d)
            ));
            r.failed = true;
            continue;
        }
        log.push(format!(
            "[M::mem_pestat] analyzing insert size distribution for orientation {}...",
            lab(d)
        ));
        q.sort_unstable(); // ks_introsort_64 on plain integers
        let n = q.len() as f64;
        // bwamem_pair.c:103-105: (int)(f * n + .499), fused.
        let at = |f: f64| q[n.mul_add(f, 0.499) as i32 as usize] as i32;
        let (p25, p50, p75) = (at(0.25), at(0.50), at(0.75));
        let iqr = (p75 - p25) as f64;
        // :106/:108 fuse p25 - 2*iqr and p75 + 2*iqr, then add .499.
        r.low = ((-iqr).mul_add(OUTLIER_BOUND, p25 as f64) + 0.499) as i32;
        if r.low < 1 {
            r.low = 1;
        }
        r.high = (iqr.mul_add(OUTLIER_BOUND, p75 as f64) + 0.499) as i32;
        log.push(format!(
            "[M::mem_pestat] (25, 50, 75) percentile: ({p25}, {p50}, {p75})"
        ));
        log.push(format!(
            "[M::mem_pestat] low and high boundaries for computing mean and std.dev: ({}, {})",
            r.low, r.high
        ));
        let (lo, hi) = (r.low as i64 as u64, r.high as i64 as u64);
        let mut x = 0i32;
        r.avg = 0.0;
        for &v in q.iter() {
            if v >= lo && v <= hi {
                r.avg += v as f64;
                x += 1;
            }
        }
        r.avg /= x as f64;
        r.std = 0.0;
        for &v in q.iter() {
            if v >= lo && v <= hi {
                let t = v as f64 - r.avg;
                r.std = t.mul_add(t, r.std); // :117 fused
            }
        }
        r.std = (r.std / x as f64).sqrt();
        log.push(format!(
            "[M::mem_pestat] mean and std.dev: ({:.2}, {:.2})",
            r.avg, r.std
        ));
        r.low = ((-iqr).mul_add(MAPPING_BOUND, p25 as f64) + 0.499) as i32;
        r.high = (iqr.mul_add(MAPPING_BOUND, p75 as f64) + 0.499) as i32;
        // :122-123 are not fused: std * 4 is computed once and reused.
        let s4 = r.std * MAX_STDDEV;
        if r.low as f64 > r.avg - s4 {
            r.low = (r.avg - s4 + 0.499) as i32;
        }
        if (r.high as f64) < r.avg + s4 {
            r.high = (r.avg + s4 + 0.499) as i32;
        }
        if r.low < 1 {
            r.low = 1;
        }
        log.push(format!(
            "[M::mem_pestat] low and high boundaries for proper pairs: ({}, {})",
            r.low, r.high
        ));
    }
    let max = isize.iter().map(|q| q.len()).max().unwrap_or(0);
    for d in 0..4 {
        if !pes[d].failed && (isize[d].len() as f64) < max as f64 * MIN_DIR_RATIO {
            pes[d].failed = true;
            log.push(format!("[M::mem_pestat] skip orientation {}", lab(d)));
        }
    }
    (pes, log)
}

/// `mem_mark_primary_se_core`.
fn mark_primary_se_core(a: &mut [Region], z: &mut Vec<usize>) {
    let tmp = (A + B).max(O_DEL + E_DEL).max(O_INS + E_INS);
    z.clear();
    z.push(0);
    for i in 1..a.len() {
        let mut hit = None;
        for (k, &j) in z.iter().enumerate() {
            let b_max = a[j].qb.max(a[i].qb);
            let e_min = a[j].qe.min(a[i].qe);
            if e_min > b_max {
                let min_l = (a[i].qe - a[i].qb).min(a[j].qe - a[j].qb);
                if (e_min - b_max) as f32 >= min_l as f32 * MASK_LEVEL {
                    if a[j].sub == 0 {
                        a[j].sub = a[i].score;
                    }
                    if a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt) {
                        a[j].sub_n += 1;
                    }
                    hit = Some(k);
                    break;
                }
            }
        }
        match hit {
            None => z.push(i),
            Some(k) => a[i].secondary = z[k] as i32,
        }
    }
}

/// `mem_mark_primary_se`: hash, sort and mark secondaries; returns the
/// number of non-ALT regions.
pub fn mark_primary_se(a: &mut [Region], id: u64) -> usize {
    let n = a.len();
    if n == 0 {
        return 0;
    }
    let mut n_pri = 0;
    for (i, r) in a.iter_mut().enumerate() {
        r.sub = 0;
        r.alt_sc = 0;
        r.secondary = -1;
        r.secondary_all = -1;
        r.hash = hash_64(id.wrapping_add(i as u64));
        if !r.is_alt {
            n_pri += 1;
        }
    }
    introsort(a, |x, y| {
        x.score > y.score
            || (x.score == y.score
                && ((x.is_alt as u8) < y.is_alt as u8 || (x.is_alt == y.is_alt && x.hash < y.hash)))
    });
    let mut z = Vec::new();
    mark_primary_se_core(a, &mut z);
    for i in 0..n {
        a[i].secondary_all = i as i32;
        let s = a[i].secondary;
        if !a[i].is_alt && s >= 0 && a[s as usize].is_alt {
            a[i].alt_sc = a[s as usize].score;
        }
    }
    if n_pri < n {
        if n_pri > 0 {
            introsort(a, |x, y| {
                (x.is_alt as u8) < y.is_alt as u8
                    || (x.is_alt == y.is_alt
                        && (x.score > y.score || (x.score == y.score && x.hash < y.hash)))
            });
        }
        let mut zz = vec![0usize; n];
        for i in 0..n {
            zz[a[i].secondary_all as usize] = i;
        }
        for i in 0..n {
            if a[i].secondary >= 0 {
                a[i].secondary_all = zz[a[i].secondary as usize] as i32;
                if a[i].is_alt {
                    a[i].secondary = i32::MAX;
                }
            } else {
                a[i].secondary_all = -1;
            }
        }
        if n_pri > 0 {
            for r in a[..n_pri].iter_mut() {
                r.sub = 0;
                r.secondary = -1;
            }
            mark_primary_se_core(&mut a[..n_pri], &mut z);
        }
    } else {
        for r in a.iter_mut() {
            r.secondary_all = r.secondary;
        }
    }
    n_pri
}

/// The orientations `mem_matesw` skips: failed ones, and those where `ma`
/// already has a hit at a proper distance from `a`. `None`: all skipped.
fn rescue_skip(l_pac: i64, pes: &[PeStat; 4], a: &Region, ma: &[Region]) -> Option<[bool; 4]> {
    let mut skip = [false; 4];
    for r in 0..4 {
        skip[r] = pes[r].failed;
    }
    for m in ma {
        let (r, dist) = infer_dir(l_pac, a.rb, m.rb);
        if dist >= pes[r].low as i64 && dist <= pes[r].high as i64 {
            skip[r] = true;
        }
    }
    (!skip.iter().all(|&s| s)).then_some(skip)
}

/// `mem_matesw`'s window for orientation `r`: whether the mate is reverse-
/// complemented and the reference interval after `bns_fetch_seq` (its bases
/// in `refseq`), or `None` when no `ksw_align2` call is made.
fn rescue_window(
    idx: &Index,
    pes: &[PeStat; 4],
    a: &Region,
    l_ms: i32,
    r: usize,
    refseq: &mut Vec<u8>,
) -> Option<(bool, i64, i64)> {
    let l_pac = idx.bns.l_pac;
    let is_rev = (r >> 1) != (r & 1);
    let is_larger = (r >> 1) == 0;
    let p = &pes[r];
    let (mut rb, mut re);
    if !is_rev {
        rb = if is_larger {
            a.rb + p.low as i64
        } else {
            a.rb - p.high as i64
        };
        re = (if is_larger {
            a.rb + p.high as i64
        } else {
            a.rb - p.low as i64
        }) + l_ms as i64;
    } else {
        rb = (if is_larger {
            a.rb + p.low as i64
        } else {
            a.rb - p.high as i64
        }) - l_ms as i64;
        re = if is_larger {
            a.rb + p.high as i64
        } else {
            a.rb - p.low as i64
        };
    }
    rb = rb.max(0);
    re = re.min(l_pac << 1);
    refseq.clear();
    // bwa leaves `rid` unset when the window is empty; it is only read when
    // rb < re, which then also fails the length test below.
    let mut rid = -1;
    if rb < re {
        let mid = (rb + re) >> 1;
        rid = idx.fetch_seq(&mut rb, mid, &mut re, refseq);
    }
    (a.rid == rid && re - rb >= MIN_SEED_LEN as i64).then_some((is_rev, rb, re))
}

/// `ksw_align2`'s `xtra` for a mate of `l_ms` bases.
fn rescue_xtra(l_ms: i32) -> i32 {
    KSW_XSUBO | KSW_XSTART | if l_ms * A < 250 { KSW_XBYTE } else { 0 } | (MIN_SEED_LEN * A)
}

/// The mate sequence as aligned: reverse-complemented for `is_rev`.
fn rescue_query(ms: &[u8], is_rev: bool, out: &mut Vec<u8>) {
    out.clear();
    if is_rev {
        out.extend(ms.iter().rev().map(|&c| if c < 4 { 3 - c } else { 4 }));
    } else {
        out.extend_from_slice(ms);
    }
}

/// Mate-rescue alignments computed ahead (on the GPU) for one pair, keyed by
/// (mate, reverse, window start, window end). A lookup miss falls back to
/// the CPU, so the set only has to be likely, not complete.
#[derive(Clone, Copy, Default)]
pub struct Rescued<'a> {
    pub keys: &'a [RescueKey],
    pub results: &'a [ksw::Kswr],
}

/// Identifies one `ksw_align2` call of a pair's mate rescue.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RescueKey {
    pub mate: u8,
    pub is_rev: bool,
    pub rb: i64,
    pub re: i64,
}

impl Rescued<'_> {
    fn get(&self, key: &RescueKey) -> Option<ksw::Kswr> {
        self.keys
            .iter()
            .position(|k| k == key)
            .map(|i| self.results[i])
    }
}

/// One planned rescue alignment: its key and its query/target bases.
pub struct RescuePlan {
    pub key: RescueKey,
    pub xtra: i32,
}

/// The `ksw_align2` calls `sam_pe` will make for this pair, assuming mate
/// rescue adds no hit that changes a later call's skip set (a superset in
/// practice; misses are recomputed on the CPU). `s` and `a` are the pair's
/// two reads and their regions. Each plan's query then target is appended to
/// `bases`, with its (offset, query length, target length) in `spans`.
pub fn plan_rescue(
    idx: &Index,
    pes: &[PeStat; 4],
    s: &[Read],
    a: &[Vec<Region>],
    plans: &mut Vec<RescuePlan>,
    bases: &mut Vec<u8>,
    spans: &mut Vec<(usize, usize, usize)>,
) {
    let mut refseq = Vec::new();
    let mut query = Vec::new();
    let first = plans.len();
    for i in 0..2 {
        let ma = &a[1 - i];
        let ms = &s[1 - i].seq;
        let l_ms = ms.len() as i32;
        let top = a[i].first().map_or(0, |r| r.score);
        for reg in a[i]
            .iter()
            .filter(|r| r.score >= top - PEN_UNPAIRED)
            .take(MAX_MATESW)
        {
            let Some(skip) = rescue_skip(idx.bns.l_pac, pes, reg, ma) else {
                continue;
            };
            for r in 0..4 {
                if skip[r] {
                    continue;
                }
                let Some((is_rev, rb, re)) = rescue_window(idx, pes, reg, l_ms, r, &mut refseq)
                else {
                    continue;
                };
                let key = RescueKey {
                    mate: (1 - i) as u8,
                    is_rev,
                    rb,
                    re,
                };
                if plans[first..].iter().any(|p| p.key == key) {
                    continue;
                }
                rescue_query(ms, is_rev, &mut query);
                let at = bases.len();
                bases.extend_from_slice(&query);
                bases.extend_from_slice(&refseq);
                spans.push((at, query.len(), refseq.len()));
                plans.push(RescuePlan {
                    key,
                    xtra: rescue_xtra(l_ms),
                });
            }
        }
    }
}

/// `mem_matesw`: rescue mate `ms` (read `mate` of the pair) near region `a`;
/// adds hits to `ma`. Alignments found in `rescued` are used as computed.
fn matesw(
    idx: &Index,
    pes: &[PeStat; 4],
    a: &Region,
    ms: &[u8],
    mate: u8,
    ma: &mut Vec<Region>,
    rescued: &Rescued,
) -> i32 {
    let l_pac = idx.bns.l_pac;
    let l_ms = ms.len() as i32;
    let Some(skip) = rescue_skip(l_pac, pes, a, ma) else {
        return 0;
    };
    let mut n = 0;
    let mut seq = Vec::new();
    let mut refseq = Vec::new();
    for r in 0..4 {
        if skip[r] {
            continue;
        }
        if let Some((is_rev, rb, re)) = rescue_window(idx, pes, a, l_ms, r, &mut refseq) {
            let key = RescueKey {
                mate,
                is_rev,
                rb,
                re,
            };
            let aln = match rescued.get(&key) {
                Some(aln) => aln,
                None => {
                    rescue_query(ms, is_rev, &mut seq);
                    ksw::align2(
                        l_ms,
                        &seq,
                        (re - rb) as i32,
                        &refseq,
                        5,
                        &MAT,
                        O_DEL,
                        E_DEL,
                        O_INS,
                        E_INS,
                        rescue_xtra(l_ms),
                        None,
                    )
                }
            };
            if aln.score >= MIN_SEED_LEN && aln.qb >= 0 {
                let mut b = Region {
                    rid: a.rid,
                    is_alt: a.is_alt,
                    qb: if is_rev { l_ms - (aln.qe + 1) } else { aln.qb },
                    qe: if is_rev { l_ms - aln.qb } else { aln.qe + 1 },
                    rb: if is_rev {
                        (l_pac << 1) - (rb + aln.te as i64 + 1)
                    } else {
                        rb + aln.tb as i64
                    },
                    re: if is_rev {
                        (l_pac << 1) - (rb + aln.tb as i64)
                    } else {
                        rb + aln.te as i64 + 1
                    },
                    score: aln.score,
                    csub: aln.score2,
                    secondary: -1,
                    ..Region::default()
                };
                b.seedcov = ((b.re - b.rb).min((b.qe - b.qb) as i64) >> 1) as i32;
                // Insert keeping `ma` sorted by score (descending).
                let at = ma
                    .iter()
                    .position(|x| x.score < b.score)
                    .unwrap_or(ma.len());
                ma.insert(at, b);
            }
            n += 1;
        }
        if n > 0 {
            sort_dedup_patch(idx, None, ma);
        }
    }
    n
}

/// `pair64_t` sorted as `ks_introsort_128`: by `x`, then `y`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Pair64 {
    x: u64,
    y: u64,
}

/// `mem_pair`: best proper pair; returns (score, sub, n_sub, z).
fn mem_pair(
    idx: &Index,
    pes: &[PeStat; 4],
    a: &[Vec<Region>; 2],
    id: u64,
    n_pri: [usize; 2],
    s: &[Read; 2],
) -> (i32, i32, i32, [usize; 2]) {
    let l_pac = idx.bns.l_pac;
    let mut v: Vec<Pair64> = Vec::new();
    for r in 0..2 {
        for i in 0..n_pri[r] {
            let e = &a[r][i];
            let fwd = if e.rb < l_pac {
                e.rb
            } else {
                (l_pac << 1) - 1 - e.rb
            };
            // Checked: a corrupted rid reports the read and region, not a
            // bare index panic in one thread (docs/DIAG.md).
            let Some(ann) = usize::try_from(e.rid)
                .ok()
                .and_then(|rid| idx.bns.anns.get(rid))
            else {
                crate::diag::bad_rid_in_mem_pair(idx, id, s, a, r, i);
            };
            let x = (e.rid as u64) << 32 | (fwd - ann.offset) as u64;
            let y =
                (e.score as u64) << 32 | (i as u64) << 2 | ((e.rb >= l_pac) as u64) << 1 | r as u64;
            v.push(Pair64 { x, y });
        }
    }
    v.sort_unstable(); // keys are distinct
    // `int id` in bwa: (id << 8) wraps in 32 bits and is sign-extended.
    let id_term = ((id as i32).wrapping_shl(8)) as i64 as u64;
    let mut u: Vec<Pair64> = Vec::new();
    let mut y = [-1i64; 4];
    for i in 0..v.len() {
        for r in 0..2 {
            let dir = r << 1 | (v[i].y >> 1 & 1) as usize;
            if pes[dir].failed {
                continue;
            }
            let which = r << 1 | ((v[i].y & 1) ^ 1) as usize;
            if y[which] < 0 {
                continue;
            }
            let mut k = y[which];
            while k >= 0 {
                let ku = k as usize;
                k -= 1;
                if (v[ku].y & 3) as usize != which {
                    continue;
                }
                let dist = v[i].x as i64 - v[ku].x as i64;
                if dist > pes[dir].high as i64 {
                    break;
                }
                if dist < pes[dir].low as i64 {
                    continue;
                }
                let ns = (dist as f64 - pes[dir].avg) / pes[dir].std;
                let sc = ((v[i].y >> 32) + (v[ku].y >> 32)) as i64 as f64;
                let e = erfc(ns.abs() * std::f64::consts::FRAC_1_SQRT_2);
                // bwamem_pair.c:244: fma(log(2*erfc)*.721, a, sum) + .499.
                let l = crate::glibc_math::log(e + e) * 0.721;
                let mut q = (l.mul_add(A as f64, sc) + 0.499) as i32;
                if q < 0 {
                    q = 0;
                }
                let py = (ku as u64) << 32 | i as u64;
                let px = (q as u64) << 32 | (hash_64(py ^ id_term) & 0xffff_ffff);
                u.push(Pair64 { x: px, y: py });
            }
        }
        y[(v[i].y & 3) as usize] = i as i64;
    }
    if u.is_empty() {
        return (0, 0, 0, [0, 0]);
    }
    let tmp = (A + B).max(O_DEL + E_DEL).max(O_INS + E_INS);
    u.sort_unstable();
    let last = u[u.len() - 1];
    let (i, k) = ((last.y >> 32) as usize, (last.y & 0xffff_ffff) as usize);
    let mut z = [0usize; 2];
    z[(v[i].y & 1) as usize] = ((v[i].y & 0xffff_ffff) >> 2) as usize;
    z[(v[k].y & 1) as usize] = ((v[k].y & 0xffff_ffff) >> 2) as usize;
    let ret = (last.x >> 32) as i32;
    let sub = if u.len() > 1 {
        (u[u.len() - 2].x >> 32) as i32
    } else {
        0
    };
    let n_sub = u[..u.len() - 1]
        .iter()
        .filter(|p| sub - (p.x >> 32) as i32 <= tmp)
        .count() as i32;
    (ret, sub, n_sub, z)
}

/// `raw_mapq` (not fused in the golden).
#[inline]
fn raw_mapq(diff: i32) -> i32 {
    (6.02 * diff as f64 / A as f64 + 0.499) as i32
}

/// `mem_sam_pe`: resolve one pair and append its SAM lines (read 1's, then
/// read 2's) to `out`. `id` is the pair's global index.
/// Set to accumulate per-thread time spent in mate rescue (debugging aid).
pub static PROFILE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
thread_local! {
    /// Seconds this thread spent in mate rescue while `PROFILE` was set.
    pub static PROF_RESCUE: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
}

#[allow(clippy::too_many_arguments)]
pub fn sam_pe(
    idx: &Index,
    rg_id: &[u8],
    pes: &[PeStat; 4],
    id: u64,
    s: &[Read; 2],
    a: &mut [Vec<Region>; 2],
    rescued: &Rescued,
    out: &mut Vec<u8>,
) {
    sam_pe_to(idx, rg_id, pes, id, s, a, rescued, out);
}

#[allow(clippy::too_many_arguments)]
pub fn bam_pe(
    idx: &Index,
    rg_id: &[u8],
    pes: &[PeStat; 4],
    id: u64,
    s: &[Read; 2],
    a: &mut [Vec<Region>; 2],
    rescued: &Rescued,
    out: &mut BamRecords,
) {
    sam_pe_to(idx, rg_id, pes, id, s, a, rescued, out);
}

#[allow(clippy::too_many_arguments)]
fn sam_pe_to<S: AlnSink>(
    idx: &Index,
    rg_id: &[u8],
    pes: &[PeStat; 4],
    id: u64,
    s: &[Read; 2],
    a: &mut [Vec<Region>; 2],
    rescued: &Rescued,
    out: &mut S,
) {
    let prof = PROFILE.load(std::sync::atomic::Ordering::Relaxed);
    let t0 = prof.then(std::time::Instant::now);
    // Mate rescue from every good-enough hit of each end.
    let b: [Vec<Region>; 2] = std::array::from_fn(|i| {
        a[i].iter()
            .filter(|r| !a[i].is_empty() && r.score >= a[i][0].score - PEN_UNPAIRED)
            .copied()
            .collect()
    });
    for i in 0..2 {
        for reg in b[i].iter().take(MAX_MATESW) {
            let (left, right) = a.split_at_mut(1);
            let ma = if i == 0 { &mut right[0] } else { &mut left[0] };
            matesw(idx, pes, reg, &s[1 - i].seq, (1 - i) as u8, ma, rescued);
        }
    }
    if let Some(t0) = t0 {
        PROF_RESCUE.with(|c| c.set(c.get() + t0.elapsed().as_secs_f64()));
    }
    crate::diag::check_pair(
        "after mate rescue (matesw), before mark_primary_se",
        idx,
        id,
        s,
        a,
    );
    let n_pri = [
        mark_primary_se(&mut a[0], id << 1),
        mark_primary_se(&mut a[1], id << 1 | 1),
    ];
    crate::diag::check_pair("after mark_primary_se, before mem_pair", idx, id, s, a);

    if n_pri[0] > 0 && n_pri[1] > 0 {
        let (o, subo0, n_sub, mut z) = mem_pair(idx, pes, a, id, n_pri, s);
        if o > 0 {
            let is_multi: [bool; 2] = std::array::from_fn(|i| {
                (1..n_pri[i]).any(|j| a[i][j].secondary < 0 && a[i][j].score >= T)
            });
            if !is_multi[0] && !is_multi[1] {
                pe_paired(idx, rg_id, s, a, n_pri, o, subo0, n_sub, &mut z, out);
                return;
            }
        }
    }
    no_pairing(idx, rg_id, pes, s, a, n_pri, out);
}

#[allow(clippy::too_many_arguments)]
fn pe_paired<S: AlnSink>(
    idx: &Index,
    rg_id: &[u8],
    s: &[Read; 2],
    a: &mut [Vec<Region>; 2],
    n_pri: [usize; 2],
    o: i32,
    subo0: i32,
    n_sub: i32,
    z: &mut [usize; 2],
    out: &mut S,
) {
    let score_un = a[0][0].score + a[1][0].score - PEN_UNPAIRED;
    let subo = subo0.max(score_un);
    let mut q_pe = raw_mapq(o - subo);
    if n_sub > 0 {
        q_pe -= crate::glibc_math::log((n_sub + 1) as f64).mul_add(4.343, 0.499) as i32; // :326 fused
    }
    q_pe = q_pe.clamp(0, 60);
    // :329: float add, then 1 - .5*f fused, then fused with .499.
    let fr = (a[0][0].frac_rep + a[1][0].frac_rep) as f64;
    let scale = (-fr).mul_add(0.5, 1.0);
    q_pe = (q_pe as f64).mul_add(scale, 0.499) as i32;
    let mut q_se = [0i32; 2];
    let mut extra_flag = 1;
    if o > score_un {
        for i in 0..2 {
            let c = &mut a[i][z[i]];
            if c.secondary >= 0 {
                let sec = c.secondary as usize;
                let sub = a[i][sec].score;
                let c = &mut a[i][z[i]];
                c.sub = sub;
                c.secondary = -2;
            }
            q_se[i] = approx_mapq_se(&a[i][z[i]]);
        }
        for q in q_se.iter_mut() {
            *q = if *q > q_pe {
                *q
            } else if q_pe < *q + 40 {
                q_pe
            } else {
                *q + 40
            };
        }
        extra_flag |= 2;
        for i in 0..2 {
            let c = &a[i][z[i]];
            q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub));
        }
    } else {
        *z = [0, 0];
        q_se = [approx_mapq_se(&a[0][0]), approx_mapq_se(&a[1][0])];
    }
    for i in 0..2 {
        let k = a[i][z[i]].secondary_all;
        if k >= 0 && (k as usize) < n_pri[i] {
            let k = k as usize;
            assert!(a[i][k].secondary_all < 0);
            for j in 0..a[i].len() {
                if a[i][j].secondary_all == k as i32 || j == k {
                    a[i][j].secondary_all = z[i] as i32;
                }
            }
            a[i][z[i]].secondary_all = -1;
        }
    }
    let xa: [Option<Vec<Option<String>>>; 2] =
        std::array::from_fn(|i| gen_alt(idx, &a[i], &s[i].seq));
    let mut h: [Aln; 2] = Default::default();
    let mut aa: [Vec<Aln>; 2] = Default::default();
    for i in 0..2 {
        h[i] = reg2aln(idx, &s[i].seq, Some(&a[i][z[i]]));
        h[i].mapq = q_se[i] as u8;
        h[i].flag |= 0x40 << i | extra_flag;
        h[i].xa = xa[i].as_ref().and_then(|x| x[z[i]].clone());
        aa[i].push(h[i].clone());
        if n_pri[i] < a[i].len() {
            let p = &a[i][n_pri[i]];
            if p.score < T || p.secondary >= 0 || !p.is_alt {
                continue;
            }
            let mut g = reg2aln(idx, &s[i].seq, Some(p));
            g.flag |= 0x800 | 0x40 << i | extra_flag;
            g.xa = xa[i].as_ref().and_then(|x| x[n_pri[i]].clone());
            aa[i].push(g);
        }
    }
    for k in 0..aa[0].len() {
        out.emit(idx, rg_id, &s[0], &aa[0], k, Some(&h[1]));
    }
    for k in 0..aa[1].len() {
        out.emit(idx, rg_id, &s[1], &aa[1], k, Some(&h[0]));
    }
}

fn no_pairing<S: AlnSink>(
    idx: &Index,
    rg_id: &[u8],
    pes: &[PeStat; 4],
    s: &[Read; 2],
    a: &[Vec<Region>; 2],
    n_pri: [usize; 2],
    out: &mut S,
) {
    let mut extra_flag = 1;
    let h: [Aln; 2] = std::array::from_fn(|i| {
        let mut which = None;
        if !a[i].is_empty() {
            if a[i][0].score >= T {
                which = Some(0);
            } else if n_pri[i] < a[i].len() && a[i][n_pri[i]].score >= T {
                which = Some(n_pri[i]);
            }
        }
        reg2aln(idx, &s[i].seq, which.map(|w| &a[i][w]))
    });
    if h[0].rid == h[1].rid && h[0].rid >= 0 {
        let (d, dist) = infer_dir(idx.bns.l_pac, a[0][0].rb, a[1][0].rb);
        if !pes[d].failed && dist >= pes[d].low as i64 && dist <= pes[d].high as i64 {
            extra_flag |= 2;
        }
    }
    reg2out(
        idx,
        rg_id,
        out,
        &s[0],
        &a[0],
        0x41 | extra_flag,
        Some(&h[1]),
    );
    reg2out(
        idx,
        rg_id,
        out,
        &s[1],
        &a[1],
        0x81 | extra_flag,
        Some(&h[0]),
    );
}
