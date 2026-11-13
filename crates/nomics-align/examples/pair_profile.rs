//! Single-threaded breakdown of the pairing phase on the first N pairs:
//! mate rescue vs everything else in sam_pe (primary marking, pairing, SAM).
//!
//! usage: pair_profile <index-hint> <R1.fastq.gz> <R2.fastq.gz> <pairs>
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use nomics_align::extend::{AlignAux, align1_core};
use nomics_align::index::Index;
use nomics_align::mem::read_batch_limited;
use nomics_align::pair::{PROF_RESCUE, PROFILE, pestat, sam_pe};
use nomics_io::gzip::GzReader;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut idx = Index::load(Path::new(&a[1])).unwrap();
    if let Some(p) = std::env::var_os("NOMICS_DENSE_SA") {
        idx.attach_dense_sa(Path::new(&p), 20).unwrap();
    }
    let n: usize = a[4].parse().unwrap();
    let open = |p: &str| BufReader::with_capacity(1 << 20, GzReader::open(Path::new(p)).unwrap());
    let reads = read_batch_limited(&mut open(&a[2]), &mut open(&a[3]), n).unwrap();
    let mut aux = AlignAux::default();
    let mut regs: Vec<_> = reads
        .iter()
        .map(|r| align1_core(&idx, &r.seq, &mut aux))
        .collect();
    let (pes, _) = pestat(idx.bns.l_pac, &regs);
    PROFILE.store(true, std::sync::atomic::Ordering::Relaxed);
    let mut out = Vec::new();
    let t = Instant::now();
    for (i, pair) in regs.chunks_exact_mut(2).enumerate() {
        let s = [reads[2 * i].clone(), reads[2 * i + 1].clone()];
        let mut ab = [std::mem::take(&mut pair[0]), std::mem::take(&mut pair[1])];
        sam_pe(&idx, b"x", &pes, i as u64, &s, &mut ab, &Default::default(), &mut out);
    }
    let total = t.elapsed().as_secs_f64();
    let rescue = PROF_RESCUE.with(|c| c.get());
    let us = |x: f64| x * 1e6 / reads.len() as f64;
    println!(
        "{} reads, per read: sam_pe {:.1}us = rescue {:.1}us + rest (mark/pair/reg2aln/SAM) {:.1}us",
        reads.len(),
        us(total),
        us(rescue),
        us(total - rescue)
    );
}
