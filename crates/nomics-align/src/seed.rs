//! Seeding: SMEMs, re-seeding and LAST-like seeds (`bwt.c` `bwt_extend`,
//! `bwt_smem1a`, `bwt_seed_strategy1`; `bwamem.c` `mem_collect_intv`).
//!
//! Gate: the oracle's `M` dump lines (every read, every interval, in order).

use crate::index::Bwt;
use crate::kmer_table::KmerTable;

/// `bwtintv_t`: a bi-interval `[x0, x0 + x2)` on the forward BWT and
/// `[x1, x1 + x2)` on the reverse; `info` holds `qbeg << 32 | qend`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Intv {
    pub x: [u64; 3],
    pub info: u64,
}

impl Intv {
    #[inline]
    pub fn qbeg(&self) -> i32 {
        (self.info >> 32) as i32
    }

    #[inline]
    pub fn qend(&self) -> i32 {
        self.info as u32 as i32
    }
}

/// bwa-mem defaults that seeding depends on (`mem_opt_init`).
pub const MIN_SEED_LEN: i32 = 19;
pub const SPLIT_WIDTH: u64 = 10;
pub const MAX_MEM_INTV: u64 = 20;

/// `(int)(opt->min_seed_len * opt->split_factor + .499)`: `int * float`
/// rounds to float, then the double addition truncates.
pub fn split_len() -> i32 {
    ((MIN_SEED_LEN as f32 * 1.5f32) as f64 + 0.499) as i32
}

/// `nst_nt4_table`: A/C/G/T (either case) to 0..3, everything else 4.
pub fn nt4(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// `bwt_set_intv`: the interval of a single base.
#[inline]
fn set_intv(bwt: &Bwt, c: u8) -> Intv {
    let c = c as usize;
    Intv {
        x: [bwt.l2[c] + 1, bwt.l2[3 - c] + 1, bwt.l2[c + 1] - bwt.l2[c]],
        info: 0,
    }
}

/// `bwt_extend`: extend `ik` by one base in all four ways; `is_back` selects
/// backward (prepend) or forward (append, via the reverse BWT) extension.
#[inline]
pub fn extend(bwt: &Bwt, ik: &Intv, is_back: bool) -> [Intv; 4] {
    // `x[!is_back]` is the side searched through occ; `x[is_back]` follows.
    let (b, nb) = if is_back { (1, 0) } else { (0, 1) };
    let k = ik.x[nb].wrapping_sub(1);
    let (tk, tl) = bwt.occ4_pair(k, k.wrapping_add(ik.x[2]));
    let mut ok = [Intv::default(); 4];
    for i in 0..4 {
        ok[i].x[nb] = bwt.l2[i] + 1 + tk[i];
        ok[i].x[2] = tl[i] - tk[i];
    }
    let spans_primary = ik.x[nb] <= bwt.primary && ik.x[nb] + ik.x[2] > bwt.primary;
    ok[3].x[b] = ik.x[b] + spans_primary as u64;
    ok[2].x[b] = ok[3].x[b] + ok[3].x[2];
    ok[1].x[b] = ok[2].x[b] + ok[2].x[2];
    ok[0].x[b] = ok[1].x[b] + ok[1].x[2];
    ok
}

/// Reusable buffers (`smem_aux_t`).
#[derive(Default)]
pub struct SmemAux {
    pub mem: Vec<Intv>,
    mem1: Vec<Intv>,
    prev: Vec<Intv>,
    curr: Vec<Intv>,
}

/// `bwt_smem1` (= `bwt_smem1a` with `max_intv = 0`): SMEMs covering `x`,
/// sorted by start; returns where the next search starts.
pub fn smem1(
    bwt: &Bwt,
    q: &[u8],
    x: usize,
    min_intv: u64,
    mem: &mut Vec<Intv>,
    prev: &mut Vec<Intv>,
    curr: &mut Vec<Intv>,
) -> usize {
    smem1_impl(bwt, q, x, min_intv, mem, prev, curr, None)
}

/// Exact table-accelerated [`smem1`]. The table supplies every prefix
/// endpoint, so the forward pass records the same interval-size changes as
/// BWA before running the unchanged backward pass.
pub fn smem1_with_table(
    bwt: &Bwt,
    q: &[u8],
    x: usize,
    min_intv: u64,
    mem: &mut Vec<Intv>,
    prev: &mut Vec<Intv>,
    curr: &mut Vec<Intv>,
    table: &KmerTable,
) -> usize {
    smem1_impl(bwt, q, x, min_intv, mem, prev, curr, Some(table))
}

#[allow(clippy::too_many_arguments)]
fn smem1_impl(
    bwt: &Bwt,
    q: &[u8],
    x: usize,
    min_intv: u64,
    mem: &mut Vec<Intv>,
    prev: &mut Vec<Intv>,
    curr: &mut Vec<Intv>,
    table: Option<&KmerTable>,
) -> usize {
    let len = q.len();
    mem.clear();
    if q[x] > 3 {
        return x + 1;
    }
    let min_intv = min_intv.max(1);
    let mut ik = set_intv(bwt, q[x]);
    ik.info = x as u64 + 1;

    // Forward search. With max_intv = 0 the "small enough" branch never fires.
    curr.clear();
    let mut i = x + 1;
    let mut stopped = false;
    if let Some(table) = table
        && x + table.k() <= len
        && q[x..x + table.k()].iter().all(|&c| c < 4)
    {
        let mut code = q[x] as u64;
        for prefix_len in 2..=table.k() {
            code = code << 2 | q[x + prefix_len - 1] as u64;
            let next = table.interval(prefix_len, code);
            if next.x[2] != ik.x[2] {
                curr.push(ik);
                if next.x[2] < min_intv {
                    i = x + prefix_len - 1;
                    stopped = true;
                    break;
                }
            }
            ik = next;
            ik.info = (x + prefix_len) as u64;
        }
        if !stopped {
            i = x + table.k();
        }
    }
    while !stopped && i < len {
        if q[i] < 4 {
            let c = (3 - q[i]) as usize;
            let ok = extend(bwt, &ik, false);
            if ok[c].x[2] != ik.x[2] {
                curr.push(ik);
                if ok[c].x[2] < min_intv {
                    break;
                }
            }
            ik = ok[c];
            ik.info = i as u64 + 1;
        } else {
            curr.push(ik);
            break;
        }
        i += 1;
    }
    if i == len {
        curr.push(ik);
    }
    curr.reverse(); // longer matches first
    let ret = curr[0].info as u32 as usize;
    std::mem::swap(prev, curr);

    // Backward search for MEMs.
    let mut i = x as isize - 1;
    while i >= -1 {
        let c: i32 = if i < 0 || q[i as usize] > 3 {
            -1
        } else {
            q[i as usize] as i32
        };
        curr.clear();
        for p in prev.iter() {
            let ok = if c >= 0 {
                Some(extend(bwt, p, true))
            } else {
                None
            };
            match ok {
                Some(ok) if ok[c as usize].x[2] >= min_intv => {
                    let o = ok[c as usize];
                    if curr.last().is_none_or(|l| o.x[2] != l.x[2]) {
                        curr.push(Intv { info: p.info, ..o });
                    }
                }
                _ => {
                    // Reached the start, an ambiguous base, or the interval
                    // vanished: `p` is a MEM unless a longer one survives.
                    if curr.is_empty() && mem.last().is_none_or(|l| ((i + 1) as u64) < l.info >> 32)
                    {
                        let mut m = *p;
                        m.info |= ((i + 1) as u64) << 32;
                        mem.push(m);
                    }
                }
            }
        }
        if curr.is_empty() {
            break;
        }
        std::mem::swap(prev, curr);
        i -= 1;
    }
    mem.reverse(); // sorted by start
    ret
}

#[allow(clippy::needless_range_loop)] // mirrors the C index loop
/// `bwt_seed_strategy1`: the shortest exact match from `x` of length at least
/// `min_len + 1` whose interval is smaller than `max_intv`.
pub fn seed_strategy1(
    bwt: &Bwt,
    q: &[u8],
    x: usize,
    min_len: usize,
    max_intv: u64,
) -> (usize, Intv) {
    seed_strategy1_impl(bwt, q, x, min_len, max_intv, None)
}

pub fn seed_strategy1_with_table(
    bwt: &Bwt,
    q: &[u8],
    x: usize,
    min_len: usize,
    max_intv: u64,
    table: &KmerTable,
) -> (usize, Intv) {
    seed_strategy1_impl(bwt, q, x, min_len, max_intv, Some(table))
}

fn seed_strategy1_impl(
    bwt: &Bwt,
    q: &[u8],
    x: usize,
    min_len: usize,
    max_intv: u64,
    table: Option<&KmerTable>,
) -> (usize, Intv) {
    let len = q.len();
    if q[x] > 3 {
        return (x + 1, Intv::default());
    }
    let mut ik = set_intv(bwt, q[x]);
    let mut start = x + 1;
    if let Some(table) = table
        // The strategy only observes an interval after `min_len` extensions.
        // Jumping through at most `min_len` therefore skips no observation.
        && table.k() <= min_len
        && x + table.k() <= len
        && q[x..x + table.k()].iter().all(|&c| c < 4)
    {
        let code = q[x..x + table.k()]
            .iter()
            .fold(0u64, |code, &c| code << 2 | c as u64);
        ik = table.interval(table.k(), code);
        start = x + table.k();
    }
    for i in start..len {
        if q[i] < 4 {
            let c = (3 - q[i]) as usize;
            let ok = extend(bwt, &ik, false);
            if ok[c].x[2] < max_intv && i - x >= min_len {
                let mut m = ok[c];
                m.info = (x as u64) << 32 | (i as u64 + 1);
                return (i + 1, m);
            }
            ik = ok[c];
        } else {
            return (i + 1, Intv::default());
        }
    }
    (len, Intv::default())
}

/// `mem_collect_intv`: all seed intervals of read `seq` (codes 0..4) into
/// `aux.mem`, sorted by `info`.
pub fn collect_intv(bwt: &Bwt, seq: &[u8], aux: &mut SmemAux) {
    collect_intv_impl(bwt, seq, aux, None);
}

/// Exact table-accelerated [`collect_intv`].
pub fn collect_intv_with_table(bwt: &Bwt, seq: &[u8], aux: &mut SmemAux, table: &KmerTable) {
    collect_intv_impl(bwt, seq, aux, Some(table));
}

fn collect_intv_impl(bwt: &Bwt, seq: &[u8], aux: &mut SmemAux, table: Option<&KmerTable>) {
    let len = seq.len();
    let min_seed_len = MIN_SEED_LEN as u32;
    let SmemAux {
        mem,
        mem1,
        prev,
        curr,
    } = aux;
    mem.clear();
    let slen = |p: &Intv| (p.info as u32).wrapping_sub((p.info >> 32) as u32);

    // Pass 1: SMEMs.
    let mut x = 0;
    while x < len {
        if seq[x] < 4 {
            x = smem1_impl(bwt, seq, x, 1, mem1, prev, curr, table);
            mem.extend(mem1.iter().filter(|p| slen(p) >= min_seed_len));
        } else {
            x += 1;
        }
    }

    // Pass 2: re-seed long SMEMs with few occurrences from their midpoint.
    let split_len = split_len();
    let old_n = mem.len();
    for k in 0..old_n {
        let p = mem[k];
        let (start, end) = (p.qbeg(), p.qend());
        if end - start < split_len || p.x[2] > SPLIT_WIDTH {
            continue;
        }
        smem1_impl(
            bwt,
            seq,
            ((start + end) >> 1) as usize,
            p.x[2] + 1,
            mem1,
            prev,
            curr,
            table,
        );
        mem.extend(mem1.iter().filter(|p| slen(p) >= min_seed_len));
    }

    // Pass 3: LAST-like seeds.
    let mut x = 0;
    while x < len {
        if seq[x] < 4 {
            let (nx, m) =
                seed_strategy1_impl(bwt, seq, x, MIN_SEED_LEN as usize, MAX_MEM_INTV, table);
            x = nx;
            if m.x[2] > 0 {
                mem.push(m);
            }
        } else {
            x += 1;
        }
    }

    // bwa uses ks_introsort on `info`; entries with equal `info` cover the
    // same query substring and so are identical intervals, so any sort
    // yields the same array.
    mem.sort_by_key(|p| p.info);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer_table::KmerTable;
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

    fn temp_path(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "nomics-{label}-{}-{}",
            std::process::id(),
            NEXT_FILE.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn tiny_bwt(text: &[u8]) -> (Bwt, PathBuf) {
        assert!(!text.is_empty() && text.len() < 16);
        let n = text.len() as u64;
        let mut suffixes: Vec<u64> = (0..=n).collect();
        suffixes.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
        let primary = suffixes.iter().position(|&p| p == 0).unwrap() as u64;
        let mut counts = [0u64; 4];
        for &base in text {
            counts[base as usize] += 1;
        }
        let l2 = [
            0,
            counts[0],
            counts[0] + counts[1],
            counts[0] + counts[1] + counts[2],
            n,
        ];
        let mut words = [0u32; 17];
        for (row, &pos) in suffixes.iter().enumerate() {
            if row as u64 == primary {
                continue;
            }
            let stored = row as u64 - (row as u64 > primary) as u64;
            let base = text[pos as usize - 1];
            words[8] |= (base as u32) << ((!stored & 15) << 1);
        }
        let path = temp_path("tiny.bwt");
        let mut file = File::create(&path).unwrap();
        file.write_all(&primary.to_le_bytes()).unwrap();
        for &value in &l2[1..] {
            file.write_all(&value.to_le_bytes()).unwrap();
        }
        for word in words {
            file.write_all(&word.to_le_bytes()).unwrap();
        }
        drop(file);
        (Bwt::load(&path).unwrap(), path)
    }

    #[test]
    fn prefix_table_matches_scalar_extension() {
        let (bwt, bwt_path) = tiny_bwt(&[0, 1, 0, 3, 2, 1, 0, 2, 3, 0, 1, 3, 2, 0, 1]);
        let cache = temp_path("kmer.cache");
        let table = KmerTable::load_or_build(&cache, &bwt, 4).unwrap();
        for len in 2..=4 {
            for code in 0..1u64 << (2 * len) {
                let first = (code >> (2 * (len - 1))) as u8;
                let mut interval = set_intv(&bwt, first);
                for at in 1..len {
                    let base = (code >> (2 * (len - at - 1)) & 3) as u8;
                    interval = extend(&bwt, &interval, false)[(3 - base) as usize];
                }
                assert_eq!(
                    table.interval(len, code).x,
                    interval.x,
                    "len={len} code={code}"
                );
            }
        }
        let loaded = KmerTable::load(&cache, &bwt, 4).unwrap();
        assert_eq!(loaded.bytes(), table.bytes());
        std::fs::remove_file(cache).unwrap();
        std::fs::remove_file(bwt_path).unwrap();
    }

    #[test]
    fn table_preserves_smem_events_and_collected_intervals() {
        let (bwt, bwt_path) = tiny_bwt(&[0, 1, 0, 3, 2, 1, 0, 2, 3, 0, 1, 3, 2, 0, 1]);
        let cache = temp_path("kmer.cache");
        let table = KmerTable::load_or_build(&cache, &bwt, 4).unwrap();
        let reads = [
            vec![
                0, 1, 0, 3, 2, 1, 0, 2, 3, 0, 1, 3, 2, 0, 1, 0, 1, 0, 3, 2, 1, 0,
            ],
            vec![
                3, 3, 3, 0, 1, 2, 4, 0, 1, 0, 3, 2, 1, 0, 2, 3, 0, 1, 3, 2, 0, 1, 2,
            ],
            vec![0; 35],
        ];
        for read in &reads {
            for x in 0..read.len() {
                for min_intv in [1, 2, 5, 11] {
                    let (mut a, mut ap, mut ac) = (Vec::new(), Vec::new(), Vec::new());
                    let (mut b, mut bp, mut bc) = (Vec::new(), Vec::new(), Vec::new());
                    let ar = smem1(&bwt, read, x, min_intv, &mut a, &mut ap, &mut ac);
                    let br =
                        smem1_with_table(&bwt, read, x, min_intv, &mut b, &mut bp, &mut bc, &table);
                    assert_eq!((ar, a), (br, b), "x={x} min_intv={min_intv}");
                }
            }
            let (mut baseline, mut accelerated) = (SmemAux::default(), SmemAux::default());
            collect_intv(&bwt, read, &mut baseline);
            collect_intv_with_table(&bwt, read, &mut accelerated, &table);
            assert_eq!(baseline.mem, accelerated.mem);
        }
        std::fs::remove_file(cache).unwrap();
        std::fs::remove_file(bwt_path).unwrap();
    }
}
