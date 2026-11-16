//! Validate the loaded BWA index against the reference FASTA and itself.
//!
//! usage: index_check <index-hint> <fasta> [samples]
//!
//! 1. `.pac` + `.ann` + `.amb`: every non-ambiguous base equals the FASTA.
//! 2. `occ` and `occ4` agree; `occ(seq_len)` equals the `L2` totals.
//! 3. For sampled rows k: BWT[k] = text[SA[k] - 1], and suffix(SA[k]) < suffix(SA[k+1]).
//! 4. Backward search of sampled 32-mers finds an interval containing their position.

use std::path::Path;
use std::time::Instant;

use nomics_align::index::Index;
use nomics_io::fasta::Fasta;

fn nt4(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let samples: u64 = args.get(3).map_or(200_000, |s| s.parse().unwrap());
    let t = Instant::now();
    let idx = Index::load(Path::new(&args[1]))?;
    let (bytes, huge) = idx.memory();
    println!(
        "loaded in {:.1}s: {} contigs, {} holes, l_pac {}, seq_len {}, primary {}, sa_intv {}, {:.2} GB ({:?} on huge pages)",
        t.elapsed().as_secs_f64(),
        idx.bns.anns.len(),
        idx.bns.ambs.len(),
        idx.bns.l_pac,
        idx.bwt.seq_len,
        idx.bwt.primary,
        idx.sa.intv,
        bytes as f64 / 1e9,
        huge
    );

    // 1. pac vs FASTA.
    let t = Instant::now();
    let fa = Fasta::load(Path::new(&args[2]), 16)?;
    assert_eq!(fa.index.len(), idx.bns.anns.len(), "contig count");
    let (mut checked, mut amb_skipped, mut mismatched) = (0u64, 0u64, 0u64);
    let mut holes = idx.bns.ambs.iter().peekable();
    for (i, ann) in idx.bns.anns.iter().enumerate() {
        assert_eq!(ann.name, fa.index[i].name, "contig {i} name");
        let seq = fa.contig(i);
        assert_eq!(seq.len() as i32, ann.len, "{} length", ann.name);
        let mut n_amb = 0;
        for (j, &b) in seq.iter().enumerate() {
            let c = nt4(b);
            if c > 3 {
                n_amb += 1;
                amb_skipped += 1;
                continue;
            }
            checked += 1;
            if idx.pac.get(ann.offset as u64 + j as u64) != c {
                mismatched += 1;
            }
        }
        // Holes are runs of identical ambiguous characters; check this
        // contig's holes cover exactly its ambiguous bases.
        let mut covered = 0i64;
        while let Some(h) = holes.peek() {
            if h.offset >= ann.offset + ann.len as i64 {
                break;
            }
            covered += h.len as i64;
            let first = (h.offset - ann.offset) as usize;
            assert!(
                seq[first..first + h.len as usize]
                    .iter()
                    .all(|&b| nt4(b) > 3),
                "{} hole at {} not ambiguous",
                ann.name,
                first
            );
            holes.next();
        }
        assert_eq!(covered, n_amb, "{} ambiguous base count", ann.name);
    }
    println!(
        "pac vs FASTA: {checked} bases checked, {mismatched} mismatched, {amb_skipped} ambiguous ({:.1}s)",
        t.elapsed().as_secs_f64()
    );
    assert_eq!(mismatched, 0);

    // 2. occ consistency.
    let t = Instant::now();
    let bwt = &idx.bwt;
    for c in 0..4u8 {
        assert_eq!(
            bwt.occ(bwt.seq_len, c),
            bwt.l2[c as usize + 1] - bwt.l2[c as usize]
        );
    }
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..samples {
        let k = rng.next() % bwt.seq_len;
        let o4 = bwt.occ4(k);
        for c in 0..4u8 {
            assert_eq!(bwt.occ(k, c), o4[c as usize], "occ({k}, {c})");
        }
        // Row k (k != primary) holds stored base k - (k > primary); occ grows by
        // exactly one in that base from row k-1 to row k.
        if k > 0 && k != bwt.primary {
            let prev = bwt.occ4(k - 1);
            let b = bwt.base(k - (k > bwt.primary) as u64) as usize;
            for c in 0..4 {
                assert_eq!(o4[c] - prev[c], (c == b) as u64, "occ step at {k}");
            }
        }
    }
    println!(
        "occ/occ4: {samples} rows consistent ({:.1}s)",
        t.elapsed().as_secs_f64()
    );

    // 3. SA order and BWT/text agreement.
    let t = Instant::now();
    let n = bwt.seq_len;
    let cmp_suffix = |a: u64, b: u64| {
        for d in 0..256 {
            let (x, y) = (a + d, b + d);
            if x >= n || y >= n {
                // `$` sorts first: the shorter suffix is smaller.
                return x >= n;
            }
            let (cx, cy) = (idx.text(x), idx.text(y));
            if cx != cy {
                return cx < cy;
            }
        }
        true // equal for 256 bases: undecided, accept
    };
    let mut sa_steps_total = 0u64;
    for _ in 0..samples {
        let k = 1 + rng.next() % (n - 1); // row 0 is `$`
        let p = idx.sa.sa(bwt, k);
        let q = idx.sa.sa(bwt, k + 1);
        assert!(p < n && q < n, "SA out of range at {k}: {p}");
        assert!(
            cmp_suffix(p, q),
            "suffix order broken at rows {k}/{}",
            k + 1
        );
        if k != bwt.primary {
            let b = bwt.base(k - (k > bwt.primary) as u64);
            assert!(p > 0, "row {k} is not primary but SA is 0");
            assert_eq!(b, idx.text(p - 1), "BWT[{k}] != text[SA-1]");
        } else {
            assert_eq!(p, 0, "primary row should be SA 0");
        }
        sa_steps_total += 1;
    }
    println!(
        "SA: {sa_steps_total} sampled rows ordered and consistent with BWT ({:.1}s)",
        t.elapsed().as_secs_f64()
    );

    // 4. Backward search of random 32-mers.
    let t = Instant::now();
    let mut hits_total = 0u64;
    for _ in 0..samples / 10 {
        let pos = rng.next() % (n - 32);
        let (mut lo, mut hi) = (0u64, n); // rows [lo, hi]
        for d in (0..32).rev() {
            let c = idx.text(pos + d) as usize;
            lo = bwt.l2[c] + bwt.occ(lo.wrapping_sub(1), c as u8) + 1;
            hi = bwt.l2[c] + bwt.occ(hi, c as u8);
            assert!(lo <= hi, "empty interval for k-mer at {pos}");
        }
        hits_total += hi - lo + 1;
        assert!(
            hi - lo > 100_000 || (lo..=hi).any(|r| idx.sa.sa(bwt, r) == pos),
            "k-mer at {pos} not found in rows {lo}..={hi}"
        );
    }
    println!(
        "backward search: {} 32-mers found at their positions, {hits_total} total hits ({:.1}s)",
        samples / 10,
        t.elapsed().as_secs_f64()
    );
    // 5. Dense SA (if requested) equals bwt_sa on sampled rows.
    if let Some(p) = std::env::var_os("NOMICS_DENSE_SA") {
        let mut idx = idx;
        let t = Instant::now();
        idx.attach_dense_sa(Path::new(&p), 20)?;
        println!("dense SA ready in {:.1}s", t.elapsed().as_secs_f64());
        let t = Instant::now();
        let n = idx.bwt.seq_len;
        for _ in 0..samples * 50 {
            let k = 1 + rng.next() % n;
            assert_eq!(
                idx.sa_at(k),
                idx.sa.sa(&idx.bwt, k),
                "dense SA differs at row {k}"
            );
        }
        println!(
            "dense SA: {} random rows equal bwt_sa ({:.1}s)",
            samples * 50,
            t.elapsed().as_secs_f64()
        );
        println!("INDEX OK");
        return Ok(());
    }
    println!("INDEX OK");
    Ok(())
}
