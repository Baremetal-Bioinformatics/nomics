//! Replay a `NOMICS_DUMP_MERGE` dump of `sort_dedup_patch` calls through the
//! current implementation and a literal copy of the one before (`reference`):
//! mismatches against the dump and against the reference (all region fields,
//! floats bitwise), and time. Since real reads rarely patch, each call with a
//! query is also replayed with one region split in two overlapping colinear
//! pieces, which the patch path merges, against the reference. Dumps over
//! 4 GiB are streamed in one pass. `MERGE_BENCH_PROFILE=1` only runs the
//! current implementation `reps` times, for a profiler.
//!
//! usage:
//!   merge_bench <index-hint> <merge.bin> [reps]
//!   NOMICS_DUMP_MERGE=<dir> merge_bench gen <index-hint> <reads> [seed] [threads]
//!   MERGE_BENCH_FASTQ=<prefix> merge_bench gen <index-hint> <reads> [seed]
//!
//! `gen` writes a WGS-like dump: `<reads>` synthetic 2x150 pairs drawn
//! uniformly from the reference (no ambiguous bases, ~1% substitutions, some
//! 1-3 bp indels, and in ~1 read in 8 a longer indel or a mismatch cluster)
//! through the paired `align_stream` pipeline (set
//! `NOMICS_GPU_EXTEND=fused` for the fused CPU stand-in). With
//! `MERGE_BENCH_FASTQ` it only writes the pairs as FASTQ.
use std::path::Path;
use std::time::Instant;

use nomics_align::extend::{Region, ext_stats, merge_dump, sort_dedup_patch};
use nomics_align::index::Index;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// Copy `src` into `out` with ~1% substitutions and ~0.4% 1-3 bp indels;
/// one read in ~8 also gets a split-making event (a 10-120 bp deletion, a
/// 10-40 bp insertion, or a dense mismatch cluster) so patching is exercised.
fn mutate(rng: &mut Rng, src: &[u8], len: usize, out: &mut Vec<u8>) -> bool {
    out.clear();
    let mut k = 0usize;
    let event_at = 20 + rng.below(len as u64 - 40) as usize;
    let mut event = rng.below(24);
    while out.len() < len && k < src.len() {
        if out.len() == event_at {
            let e = std::mem::replace(&mut event, 99);
            match e {
                0 => k += 10 + rng.below(111) as usize,
                1 => out.extend((0..10 + rng.below(31)).map(|_| rng.below(4) as u8)),
                2 => {
                    for _ in 0..12 {
                        if k < src.len() {
                            let mm = rng.below(2) == 0;
                            out.push(if mm {
                                (src[k] + 1 + rng.below(3) as u8) % 4
                            } else {
                                src[k]
                            });
                            k += 1;
                        }
                    }
                }
                _ => {}
            }
        }
        if out.len() >= len || k >= src.len() {
            break;
        }
        let x = rng.below(1000);
        if x < 10 {
            out.push(((src[k] as u64 + 1 + rng.below(3)) % 4) as u8);
            k += 1;
        } else if x < 12 {
            k += 1 + rng.below(3) as usize;
        } else if x < 14 {
            for _ in 0..1 + rng.below(3) {
                out.push(rng.below(4) as u8);
            }
        } else {
            out.push(src[k]);
            k += 1;
        }
    }
    out.truncate(len);
    out.len() == len
}

fn fastq(out: &mut Vec<u8>, name: usize, seq: &[u8]) {
    out.extend_from_slice(format!("@s{name}\n").as_bytes());
    out.extend(seq.iter().map(|&c| b"ACGTN"[c as usize]));
    out.extend_from_slice(b"\n+\n");
    out.extend(std::iter::repeat_n(b'I', seq.len()));
    out.push(b'\n');
}

/// Paired 2x150 reads from fragments of ~400 bp (sd ~50) drawn uniformly
/// over the forward strand, skipping ambiguous bases and contig ends,
/// through the default `align_stream` pipeline.
fn synth(idx: &Index, n: usize, seed: u64, threads: usize) {
    let fastq_out = std::env::var_os("MERGE_BENCH_FASTQ");
    assert!(
        fastq_out.is_some() || std::env::var_os("NOMICS_DUMP_MERGE").is_some(),
        "gen needs NOMICS_DUMP_MERGE=<dir> (or MERGE_BENCH_FASTQ=<prefix>)"
    );
    let l_pac = idx.bns.l_pac;
    let ambs = &idx.bns.ambs;
    let in_hole = |b: i64, e: i64| {
        ambs.iter()
            .any(|h| h.offset < e && b < h.offset + h.len as i64)
    };
    let mut rng = Rng(seed | 1);
    let (mut f1, mut f2) = (Vec::new(), Vec::new());
    let (mut frag, mut r, mut rc) = (Vec::new(), Vec::new(), Vec::new());
    const L: usize = 150;
    let mut made = 0usize;
    while made < n {
        // ~N(400, 50) from a sum of uniforms.
        let flen = 300 + (0..4).map(|_| rng.below(51) as i64).sum::<i64>() + 10;
        let p = rng.below((l_pac - 1000) as u64) as i64;
        if in_hole(p, p + flen) || idx.bns.pos2rid(p) != idx.bns.pos2rid(p + flen) {
            continue;
        }
        idx.get_seq(p, p + flen, &mut frag);
        if !mutate(&mut rng, &frag, L, &mut r) {
            continue;
        }
        rc.clear();
        rc.extend(frag.iter().rev().map(|&c| 3 - c));
        let mut r2 = Vec::new();
        if !mutate(&mut rng, &rc, L, &mut r2) {
            continue;
        }
        if rng.below(2) == 0 {
            std::mem::swap(&mut r, &mut r2);
        }
        fastq(&mut f1, made, &r);
        fastq(&mut f2, made, &r2);
        made += 1;
    }
    if let Some(prefix) = fastq_out {
        // Only write the pairs, as `<prefix>_R1.fastq` / `<prefix>_R2.fastq`.
        let p = prefix.to_str().unwrap();
        std::fs::write(format!("{p}_R1.fastq"), &f1).unwrap();
        std::fs::write(format!("{p}_R2.fastq"), &f2).unwrap();
        return;
    }
    let t = Instant::now();
    let (mut c1, mut c2) = (std::io::Cursor::new(f1), std::io::Cursor::new(f2));
    let mut sam = Vec::new();
    nomics_align::mem::align_stream(idx, b"synth", &mut c1, &mut c2, threads, false, &mut sam)
        .unwrap();
    eprintln!(
        "aligned {n} synthetic pairs in {:.1}s ({} SAM bytes)",
        t.elapsed().as_secs_f64(),
        sam.len()
    );
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a[1] == "gen" {
        let idx = Index::load(Path::new(&a[2])).unwrap();
        let n: usize = a[3].parse().unwrap();
        let seed = a.get(4).map_or(0x5eed, |s| s.parse().unwrap());
        let threads = a.get(5).map_or(8, |s| s.parse().unwrap());
        synth(&idx, n, seed, threads);
        return;
    }
    assert!(
        std::env::var_os("NOMICS_DUMP_MERGE").is_none(),
        "unset NOMICS_DUMP_MERGE for a replay"
    );
    let idx = Index::load(Path::new(&a[1])).unwrap();
    if std::fs::metadata(&a[2]).unwrap().len() > 4 << 30 {
        stream(&idx, Path::new(&a[2]));
        return;
    }
    let bytes = std::fs::read(&a[2]).unwrap();
    let reps: usize = a.get(3).map_or(3, |s| s.parse().unwrap());
    let recs = merge_dump::read_all(&bytes);
    let mut ns: Vec<usize> = recs.iter().map(|r| r.before.len()).collect();
    ns.sort_unstable();
    let total: usize = ns.iter().sum();
    let pct = |p: f64| ns[((ns.len() - 1) as f64 * p) as usize];
    println!(
        "{} calls ({} with query), {} regions; regions/call mean {:.2} p50 {} p99 {} p99.9 {} max {}",
        recs.len(),
        recs.iter().filter(|r| r.query.is_some()).count(),
        total,
        total as f64 / recs.len().max(1) as f64,
        pct(0.5),
        pct(0.99),
        pct(0.999),
        ns.last().copied().unwrap_or(0)
    );
    // Prefault the reference the patches read.
    let mut w = Vec::new();
    for r in &recs {
        if let (Some(_), Some(x)) = (&r.query, r.before.first()) {
            idx.get_seq(x.rb, x.re, &mut w);
        }
    }

    if std::env::var_os("MERGE_BENCH_PROFILE").is_some() {
        // Only the current implementation, `reps` times (for a profiler).
        let mut v = Vec::new();
        for _ in 0..reps {
            for r in &recs {
                v.clear();
                v.extend_from_slice(&r.before);
                sort_dedup_patch(&idx, r.query.as_deref(), &mut v);
                std::hint::black_box(&v);
            }
        }
        return;
    }
    // Mismatch pass, with the merge counters on: the current implementation
    // against the dump and against the literal reference copy below.
    ext_stats::ON.store(true, std::sync::atomic::Ordering::Relaxed);
    let (mut bad, mut bad_ref, mut ref_vs_dump) = (0usize, 0usize, 0usize);
    let (mut v, mut u) = (Vec::new(), Vec::new());
    let same = |x: &[Region], y: &[Region]| {
        x.len() == y.len() && x.iter().zip(y).all(|(p, q)| merge_dump::same(p, q))
    };
    for (k, r) in recs.iter().enumerate() {
        v.clear();
        v.extend_from_slice(&r.before);
        sort_dedup_patch(&idx, r.query.as_deref(), &mut v);
        u.clear();
        u.extend_from_slice(&r.before);
        reference::sort_dedup_patch(&idx, r.query.as_deref(), &mut u);
        if !same(&v, &r.after) {
            if bad < 5 {
                eprintln!(
                    "mismatch at call {k}: got {} regions, want {}\n  got  {:?}\n  want {:?}",
                    v.len(),
                    r.after.len(),
                    v,
                    r.after
                );
            }
            bad += 1;
        }
        bad_ref += usize::from(!same(&v, &u));
        ref_vs_dump += usize::from(!same(&u, &r.after));
    }
    ext_stats::ON.store(false, std::sync::atomic::Ordering::Relaxed);
    println!(
        "mismatches: {bad} vs dump, {bad_ref} vs reference / {} calls (reference vs dump: {ref_vs_dump})",
        recs.len()
    );
    println!("{}", ext_stats::merge_report());

    // Patches are rare in real reads, so also split one region of each call
    // with a query into two colinear pieces overlapping by a few bases (the
    // second shifted on the reference by -3..=3, an indel) and compare the
    // current implementation with the reference on those calls.
    let mut rng = Rng(0xfeed_beef);
    let (mut n_split, mut bad_split) = (0usize, 0usize);
    let mut splits: Vec<(&[u8], Vec<Region>)> = Vec::new();
    let ok0 = ext_stats::PATCH_OK.load(std::sync::atomic::Ordering::Relaxed);
    ext_stats::ON.store(true, std::sync::atomic::Ordering::Relaxed);
    for r in &recs {
        let Some(q) = r.query.as_deref() else {
            continue;
        };
        let Some(k) = r
            .before
            .iter()
            .position(|x| x.qe - x.qb >= 60 && x.re - x.rb >= 60)
        else {
            continue;
        };
        let x = r.before[k];
        let cut = x.qb + 20 + rng.below((x.qe - x.qb - 40) as u64) as i32;
        let o = 1 + rng.below(12) as i32;
        let d = rng.below(7) as i64 - 3;
        let part =
            |qb: i32, qe: i32| ((qe - qb) as f64 / (x.qe - x.qb) as f64 * x.score as f64) as i32;
        let mut a = x;
        a.qe = cut + o;
        a.re = x.rb + (a.qe - x.qb) as i64;
        a.score = part(a.qb, a.qe);
        a.truesc = a.score;
        let mut b = x;
        b.qb = cut - o;
        b.rb = x.rb + (b.qb - x.qb) as i64 + d;
        b.re = x.re + d;
        b.score = part(b.qb, b.qe);
        b.truesc = b.score;
        if !(a.rb < b.rb && a.re < b.re && b.rb < b.re) {
            continue;
        }
        let mut input = r.before.clone();
        input[k] = a;
        input.push(b);
        v.clear();
        v.extend_from_slice(&input);
        sort_dedup_patch(&idx, Some(q), &mut v);
        u.clear();
        u.extend_from_slice(&input);
        reference::sort_dedup_patch(&idx, Some(q), &mut u);
        n_split += 1;
        bad_split += usize::from(!same(&v, &u));
        splits.push((q, input));
    }
    ext_stats::ON.store(false, std::sync::atomic::Ordering::Relaxed);
    println!(
        "split calls: {n_split}, {} patches, {bad_split} mismatches vs reference",
        ext_stats::PATCH_OK.load(std::sync::atomic::Ordering::Relaxed) - ok0
    );
    bad_ref += bad_split;
    let time_splits = |f: &dyn Fn(&mut Vec<Region>, &[u8])| {
        let mut v = Vec::new();
        let mut best = f64::MAX;
        for _ in 0..reps {
            let t = Instant::now();
            for (q, input) in &splits {
                v.clear();
                v.extend_from_slice(input);
                f(&mut v, q);
                std::hint::black_box(&v);
            }
            best = best.min(t.elapsed().as_secs_f64());
        }
        best
    };
    let split_cur = time_splits(&|v, q| sort_dedup_patch(&idx, Some(q), v));
    let split_old = time_splits(&|v, q| reference::sort_dedup_patch(&idx, Some(q), v));
    println!(
        "split calls best of {reps}: current {:.1} ms, reference {:.1} ms",
        split_cur * 1e3,
        split_old * 1e3
    );

    type Merge<'a> = &'a dyn Fn(&mut Vec<Region>, Option<&[u8]>);
    let time = |f: Merge, patch: bool| {
        let mut v = Vec::new();
        let mut best = f64::MAX;
        for _ in 0..reps {
            let t = Instant::now();
            for r in &recs {
                v.clear();
                v.extend_from_slice(&r.before);
                f(&mut v, if patch { r.query.as_deref() } else { None });
                std::hint::black_box(&v);
            }
            best = best.min(t.elapsed().as_secs_f64());
        }
        best
    };
    let copy = time(&|_, _| {}, true);
    let cur = time(&|v, q| sort_dedup_patch(&idx, q, v), true);
    let cur_np = time(&|v, q| sort_dedup_patch(&idx, q, v), false);
    let old = time(&|v, q| reference::sort_dedup_patch(&idx, q, v), true);
    let ms = |x: f64| (x - copy) * 1e3;
    println!(
        "best of {reps}, input copy ({:.1} ms) subtracted: current {:.1} ms ({:.0} ns/call; without patching {:.1} ms), reference {:.1} ms; {:.2}x",
        copy * 1e3,
        ms(cur),
        (cur - copy) * 1e9 / recs.len().max(1) as f64,
        ms(cur_np),
        ms(old),
        (old - copy) / (cur - copy)
    );
    if bad + bad_ref > 0 {
        std::process::exit(1);
    }
}

/// A dump too large to hold: one pass, each call through the current
/// implementation and the reference, both checked against the dump and
/// timed per call.
fn stream(idx: &Index, path: &Path) {
    let mut r = std::io::BufReader::with_capacity(1 << 24, std::fs::File::open(path).unwrap());
    let same = |x: &[Region], y: &[Region]| {
        x.len() == y.len() && x.iter().zip(y).all(|(p, q)| merge_dump::same(p, q))
    };
    let (mut calls, mut regions, mut bad, mut bad_ref) = (0usize, 0usize, 0usize, 0usize);
    let (mut t_cur, mut t_old) = (0f64, 0f64);
    let (mut v, mut u) = (Vec::new(), Vec::new());
    while let Some(rec) = merge_dump::read_record(&mut r).unwrap() {
        calls += 1;
        regions += rec.before.len();
        v.clear();
        v.extend_from_slice(&rec.before);
        let t = Instant::now();
        sort_dedup_patch(idx, rec.query.as_deref(), &mut v);
        t_cur += t.elapsed().as_secs_f64();
        u.clear();
        u.extend_from_slice(&rec.before);
        let t = Instant::now();
        reference::sort_dedup_patch(idx, rec.query.as_deref(), &mut u);
        t_old += t.elapsed().as_secs_f64();
        bad += usize::from(!same(&v, &rec.after));
        bad_ref += usize::from(!same(&u, &rec.after));
    }
    println!(
        "streamed {calls} calls, {regions} regions: mismatches {bad} (current vs dump), {bad_ref} (reference vs dump); current {:.0} ms, reference {:.0} ms, {:.2}x",
        t_cur * 1e3,
        t_old * 1e3,
        t_old / t_cur
    );
    if bad + bad_ref > 0 {
        std::process::exit(1);
    }
}

/// The implementation before this change (commit 1cb8903), copied literally:
/// `mem_sort_dedup_patch` with whole-region `introsort`s, `mem_patch_reg`
/// through the CIGAR-building `gen_cigar2`, and a per-base `bns_get_seq`.
mod reference {
    use nomics_align::extend::{E_DEL, E_INS, MAT, O_DEL, O_INS, Region, W};
    use nomics_align::index::Index;
    use nomics_align::ksort::introsort;
    use nomics_align::ksw;

    const MAX_CHAIN_GAP: i64 = 10000;
    const MASK_LEVEL_REDUN: f32 = 0.95;
    const PATCH_MAX_R_BW: f32 = 0.05;
    const PATCH_MIN_SC_RATIO: f32 = 0.90;

    fn get_seq(idx: &Index, mut beg: i64, mut end: i64, out: &mut Vec<u8>) {
        out.clear();
        let l_pac = idx.bns.l_pac;
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
                    out.push(3 - idx.pac.get(k as u64));
                    k -= 1;
                }
            } else {
                out.extend((beg..end).map(|k| idx.pac.get(k as u64)));
            }
        }
    }

    /// `gen_cigar2(.., want_md = false)`'s score.
    fn gen_cigar2_score(idx: &Index, w_: i32, query: &[u8], rb: i64, re: i64) -> Option<i32> {
        let l_pac = idx.bns.l_pac;
        let l_query = query.len() as i32;
        if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
            return None;
        }
        let mut rseq = Vec::new();
        get_seq(idx, rb, re, &mut rseq);
        let rlen = rseq.len() as i64;
        if re - rb != rlen {
            return None;
        }
        let mut q = query.to_vec();
        if rb >= l_pac {
            q.reverse();
            rseq.reverse();
        }
        let (score, _cigar) = if l_query as i64 == re - rb && w_ == 0 {
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
                &q,
                rlen as i32,
                &rseq,
                5,
                &MAT,
                O_DEL,
                E_DEL,
                O_INS,
                E_INS,
                w,
            )
        };
        Some(score)
    }

    fn patch_reg(idx: &Index, query: &[u8], a: &Region, b: &Region) -> Option<(i32, i32)> {
        let l_pac = idx.bns.l_pac;
        if a.rb < l_pac && b.rb >= l_pac {
            return None;
        }
        if a.qb >= b.qb || a.qe >= b.qe || a.re >= b.re {
            return None;
        }
        let mut w = (((a.re - b.rb) - (a.qe - b.qb) as i64) as i32).abs();
        let r = ((a.re - b.rb) as f64 / (b.re - a.rb) as f64
            - (a.qe - b.qb) as f64 / (b.qe - a.qb) as f64)
            .abs();
        if a.re < b.rb || a.qe < b.qb {
            if w > W << 1 || r >= PATCH_MAX_R_BW as f64 {
                return None;
            }
        } else if w > W << 2 || r >= (PATCH_MAX_R_BW * 2.0) as f64 {
            return None;
        }
        w += a.w + b.w;
        w = w.min(W << 2);
        let score = gen_cigar2_score(idx, w, &query[a.qb as usize..b.qe as usize], a.rb, b.re)
            .expect("mem_patch_reg: bwa_gen_cigar2 left the score unset");
        let sc = (b.score + a.score) as f64;
        let q_s = ((b.qe - a.qb) as f64 / ((b.qe - b.qb) + (a.qe - a.qb)) as f64).mul_add(sc, 0.499)
            as i32;
        let r_s = ((b.re - a.rb) as f64 / ((b.re - b.rb) + (a.re - a.rb)) as f64).mul_add(sc, 0.499)
            as i32;
        if (score as f64 / q_s.max(r_s) as f64) < PATCH_MIN_SC_RATIO as f64 {
            return None;
        }
        Some((score, w))
    }

    pub fn sort_dedup_patch(idx: &Index, query: Option<&[u8]>, a: &mut Vec<Region>) {
        let n = a.len();
        if n <= 1 {
            return;
        }
        introsort(a, |x, y| x.re < y.re);
        for r in a.iter_mut() {
            r.n_comp = 1;
        }
        for i in 1..n {
            let p = a[i];
            if p.rid != a[i - 1].rid || p.rb >= a[i - 1].re + MAX_CHAIN_GAP {
                continue;
            }
            let mut j = i as isize - 1;
            while j >= 0 {
                let ju = j as usize;
                let (p, q) = (a[i], a[ju]);
                if !(p.rid == q.rid && p.rb < q.re + MAX_CHAIN_GAP) {
                    break;
                }
                j -= 1;
                if q.qe == q.qb {
                    continue;
                }
                let or = q.re - p.rb;
                let oq = if q.qb < p.qb {
                    q.qe - p.qb
                } else {
                    p.qe - q.qb
                } as i64;
                let mr = (q.re - q.rb).min(p.re - p.rb);
                let mq = ((q.qe - q.qb).min(p.qe - p.qb)) as i64;
                if or as f32 > MASK_LEVEL_REDUN * mr as f32
                    && oq as f32 > MASK_LEVEL_REDUN * mq as f32
                {
                    if p.score < q.score {
                        a[i].qe = a[i].qb;
                        break;
                    } else {
                        a[ju].qe = a[ju].qb;
                    }
                } else if q.rb < p.rb
                    && let Some(query) = query
                    && let Some((score, w)) = patch_reg(idx, query, &q, &p)
                {
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
        a.retain(|r| r.qe > r.qb);
        introsort(a, |x, y| {
            x.score > y.score
                || (x.score == y.score && (x.rb < y.rb || (x.rb == y.rb && x.qb < y.qb)))
        });
        for i in 1..a.len() {
            if a[i].score == a[i - 1].score && a[i].rb == a[i - 1].rb && a[i].qb == a[i - 1].qb {
                a[i].qe = a[i].qb;
            }
        }
        let mut m = 1;
        for i in 1..a.len() {
            if a[i].qe > a[i].qb {
                a[m] = a[i];
                m += 1;
            }
        }
        a.truncate(m);
    }
}
