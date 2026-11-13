//! Chaining and chain filtering (`bwamem.c` `mem_chain`, `test_and_merge`,
//! `mem_chain_weight`, `mem_chain_flt`).
//!
//! Gates: the oracle's `c` (chains in B-tree order) and `C` (filtered) dumps.

use crate::index::{Bns, Index};
use crate::kbtree::KBTree;
use crate::ksort::introsort;
use crate::seed::{MIN_SEED_LEN, SmemAux, collect_intv};

/// bwa-mem defaults used here (`mem_opt_init`).
const MAX_OCC: u64 = 500;
const W: i64 = 100;
const MAX_CHAIN_GAP: i64 = 10000;
const MASK_LEVEL: f32 = 0.50;
const DROP_RATIO: f32 = 0.50;
const MIN_CHAIN_WEIGHT: u32 = 0;
const MAX_CHAIN_EXTEND: i32 = 1 << 30;

/// `KB_DEFAULT_SIZE` and `sizeof(mem_chain_t)` on LP64 (aarch64 and x86-64):
/// they fix the B-tree's node size, and with duplicate keys the node size
/// decides which chain a seed is offered to.
const KB_DEFAULT_SIZE: usize = 512;
const SIZEOF_MEM_CHAIN_T: usize = 40;

/// `mem_seed_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seed {
    pub rbeg: i64,
    pub qbeg: i32,
    pub len: i32,
    pub score: i32,
}

/// `mem_chain_t`.
#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub pos: i64,
    pub rid: i32,
    pub first: i32,
    pub w: u32,
    pub kept: u8,
    pub is_alt: bool,
    pub frac_rep: f32,
    pub seeds: Vec<Seed>,
}

impl Chain {
    fn beg(&self) -> i32 {
        self.seeds[0].qbeg
    }

    fn end(&self) -> i32 {
        let s = self.seeds.last().unwrap();
        s.qbeg + s.len
    }
}

/// `test_and_merge`: true if `p` is absorbed by (contained in, or appended
/// to) chain `c`.
fn test_and_merge(l_pac: i64, c: &mut Chain, p: &Seed, seed_rid: i32) -> bool {
    let last = *c.seeds.last().unwrap();
    let qend = (last.qbeg + last.len) as i64;
    let rend = last.rbeg + last.len as i64;
    if seed_rid != c.rid {
        return false;
    }
    let first = c.seeds[0];
    if p.qbeg >= first.qbeg
        && (p.qbeg + p.len) as i64 <= qend
        && p.rbeg >= first.rbeg
        && p.rbeg + p.len as i64 <= rend
    {
        return true; // contained seed
    }
    if (last.rbeg < l_pac || first.rbeg < l_pac) && p.rbeg >= l_pac {
        return false; // different strand
    }
    let x = (p.qbeg - last.qbeg) as i64;
    let y = p.rbeg - last.rbeg;
    if y >= 0
        && x - y <= W
        && y - x <= W
        && x - (last.len as i64) < MAX_CHAIN_GAP
        && y - (last.len as i64) < MAX_CHAIN_GAP
    {
        c.seeds.push(*p);
        return true;
    }
    false
}

/// `mem_chain_weight`: the smaller of the query and reference coverage.
pub fn chain_weight(c: &Chain) -> u32 {
    let cover = |b: &dyn Fn(&Seed) -> i64| {
        let (mut end, mut w) = (0i64, 0i64);
        for s in &c.seeds {
            let (sb, se) = (b(s), b(s) + s.len as i64);
            if sb >= end {
                w += s.len as i64;
            } else if se > end {
                w += se - end;
            }
            end = end.max(se);
        }
        w as i32
    };
    let w = cover(&|s| s.qbeg as i64).min(cover(&|s| s.rbeg));
    w.min((1 << 30) - 1) as u32
}

/// `mem_chain`: seeds from `collect_intv`, placed and chained in bwa's B-tree.
pub fn mem_chain(idx: &Index, seq: &[u8], aux: &mut SmemAux) -> Vec<Chain> {
    if (seq.len() as i32) < MIN_SEED_LEN {
        return Vec::new();
    }
    collect_intv(&idx.bwt, seq, aux);
    let mem = std::mem::take(&mut aux.mem);
    let chains = mem_chain_from(idx, seq.len() as i32, &mem);
    aux.mem = mem;
    chains
}

/// `mem_chain` given the read's seed intervals (from `collect_intv`, on the
/// CPU or the GPU).
pub fn mem_chain_from(idx: &Index, len: i32, mem: &[crate::seed::Intv]) -> Vec<Chain> {
    if len < MIN_SEED_LEN {
        return Vec::new();
    }
    let bns: &Bns = &idx.bns;

    // frac_rep: query coverage of intervals with more than max_occ hits.
    let (mut b, mut e, mut l_rep) = (0i32, 0i32, 0i32);
    for p in mem.iter().filter(|p| p.x[2] > MAX_OCC) {
        let (sb, se) = (p.qbeg(), p.qend());
        if sb > e {
            l_rep += e - b;
            b = sb;
            e = se;
        } else {
            e = e.max(se);
        }
    }
    l_rep += e - b;

    // The tree stores (pos, arena index) handles; `chain_cmp` compares pos only.
    let mut arena: Vec<Chain> = Vec::new();
    let mut tree = KBTree::new(
        KB_DEFAULT_SIZE,
        SIZEOF_MEM_CHAIN_T,
        |a: &(i64, u32), b: &(i64, u32)| (b.0 < a.0) as i32 - (a.0 < b.0) as i32,
    )
    .expect("B-tree degree >= 2");
    for p in mem {
        let slen = p.qend() - p.qbeg();
        let step = if p.x[2] > MAX_OCC {
            p.x[2] / MAX_OCC
        } else {
            1
        };
        let (mut k, mut count) = (0u64, 0u64);
        while k < p.x[2] && count < MAX_OCC {
            let rbeg = idx.sa_at(p.x[0] + k) as i64;
            let s = Seed {
                rbeg,
                qbeg: p.qbeg(),
                len: slen,
                score: slen,
            };
            let rid = bns.intv2rid(s.rbeg, s.rbeg + s.len as i64);
            if rid >= 0 {
                let merged = !tree.is_empty()
                    && match tree.intervalp(&(rbeg, 0)).0 {
                        Some(&(_, lower)) => {
                            test_and_merge(bns.l_pac, &mut arena[lower as usize], &s, rid)
                        }
                        None => false,
                    };
                if !merged {
                    let mut seeds = Vec::with_capacity(4);
                    seeds.push(s);
                    tree.putp((rbeg, arena.len() as u32));
                    arena.push(Chain {
                        pos: rbeg,
                        rid,
                        first: 0,
                        w: 0,
                        kept: 0,
                        is_alt: bns.anns[rid as usize].is_alt,
                        frac_rep: 0.0,
                        seeds,
                    });
                }
            }
            k += step;
            count += 1;
        }
    }
    let frac_rep = l_rep as f32 / len as f32;
    let mut order = Vec::with_capacity(arena.len());
    tree.traverse(|&(_, i)| order.push(i));
    let mut slots: Vec<Option<Chain>> = arena.into_iter().map(Some).collect();
    order
        .into_iter()
        .map(|i| {
            let mut c = slots[i as usize].take().unwrap();
            c.frac_rep = frac_rep;
            c
        })
        .collect()
}

/// `mem_chain_flt`: drop chains shadowed by heavier overlapping chains.
pub fn mem_chain_flt(a: &mut Vec<Chain>) {
    if a.is_empty() {
        return;
    }
    for c in a.iter_mut() {
        c.first = -1;
        c.kept = 0;
        c.w = chain_weight(c);
    }
    #[allow(clippy::absurd_extreme_comparisons)] // min_chain_weight is 0 by default
    a.retain(|c| c.w >= MIN_CHAIN_WEIGHT);
    if a.is_empty() {
        return;
    }
    // ks_introsort(mem_flt): moves depend only on the comparator, so sorting
    // (w, index) handles yields bwa's permutation.
    let mut h: Vec<(u32, u32)> = a.iter().enumerate().map(|(i, c)| (c.w, i as u32)).collect();
    introsort(&mut h, |x, y| x.0 > y.0);
    let mut slots: Vec<Option<Chain>> = std::mem::take(a).into_iter().map(Some).collect();
    a.extend(h.iter().map(|&(_, i)| slots[i as usize].take().unwrap()));
    a[0].kept = 3;
    let mut chains = vec![0usize];
    for i in 1..a.len() {
        let mut large_ovlp = false;
        let mut dropped = false;
        for &j in &chains {
            let b_max = a[j].beg().max(a[i].beg());
            let e_min = a[j].end().min(a[i].end());
            if e_min > b_max && (!a[j].is_alt || a[i].is_alt) {
                let li = a[i].end() - a[i].beg();
                let lj = a[j].end() - a[j].beg();
                let min_l = li.min(lj);
                if (e_min - b_max) as f32 >= min_l as f32 * MASK_LEVEL
                    && (min_l as i64) < MAX_CHAIN_GAP
                {
                    large_ovlp = true;
                    if a[j].first < 0 {
                        a[j].first = i as i32;
                    }
                    if (a[i].w as f32) < a[j].w as f32 * DROP_RATIO
                        && a[j].w as i64 - a[i].w as i64 >= (MIN_SEED_LEN << 1) as i64
                    {
                        dropped = true;
                        break;
                    }
                }
            }
        }
        if !dropped {
            chains.push(i);
            a[i].kept = if large_ovlp { 2 } else { 3 };
        }
    }
    for &j in &chains {
        let f = a[j].first;
        if f >= 0 {
            a[f as usize].kept = 1;
        }
    }
    // max_chain_extend cap (never reached with the default 1<<30).
    let mut i = 0;
    let mut k = 0;
    while i < a.len() {
        if a[i].kept != 0 && a[i].kept != 3 {
            k += 1;
            if k >= MAX_CHAIN_EXTEND {
                break;
            }
        }
        i += 1;
    }
    for c in a.iter_mut().skip(i) {
        if c.kept < 3 {
            c.kept = 0;
        }
    }
    a.retain(|c| c.kept != 0);
}
