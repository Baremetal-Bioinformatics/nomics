//! Single-threaded time breakdown of align1_core's stages on the first N pairs.
//!
//! usage: align_profile <index-hint> <R1.fastq.gz> <R2.fastq.gz> <pairs>
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use nomics_align::chain::{mem_chain, mem_chain_flt};
use nomics_align::extend::{AlignAux, chain2aln, sort_dedup_patch};
use nomics_align::index::Index;
use nomics_align::mem::read_fastq;
use nomics_align::seed::{SmemAux, collect_intv};
use nomics_io::gzip::GzReader;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut idx = Index::load(Path::new(&a[1])).unwrap();
    if let Some(p) = std::env::var_os("NOMICS_DENSE_SA") {
        let t = Instant::now();
        idx.attach_dense_sa(Path::new(&p), 20).unwrap();
        eprintln!("dense SA ready in {:.1}s", t.elapsed().as_secs_f64());
    }
    let n: usize = a[4].parse().unwrap();
    let mut reads = Vec::new();
    let mut line = Vec::new();
    for p in [&a[2], &a[3]] {
        let mut r = BufReader::new(GzReader::open(Path::new(p)).unwrap());
        for _ in 0..n {
            reads.push(read_fastq(&mut r, &mut line).unwrap().unwrap().seq);
        }
    }
    let mut smem = SmemAux::default();
    let t = Instant::now();
    let mut n_mem = 0;
    for s in &reads {
        collect_intv(&idx.bwt, s, &mut smem);
        n_mem += smem.mem.len();
    }
    let t_seed = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut chains: Vec<_> = reads
        .iter()
        .map(|s| mem_chain(&idx, s, &mut smem))
        .collect();
    let t_chain = t.elapsed().as_secs_f64() - t_seed;
    let n_c: usize = chains.iter().map(|c| c.len()).sum();
    let t = Instant::now();
    for c in chains.iter_mut() {
        mem_chain_flt(c);
    }
    let t_flt = t.elapsed().as_secs_f64();
    let n_cf: usize = chains.iter().map(|c| c.len()).sum();
    let mut aux = AlignAux::default();
    let t = Instant::now();
    let mut regs: Vec<_> = reads
        .iter()
        .zip(&chains)
        .map(|(s, cs)| {
            let mut r = Vec::new();
            for c in cs {
                chain2aln(&idx, s, c, &mut r, &mut aux, 0, &Default::default());
            }
            r
        })
        .collect();
    let t_ext = t.elapsed().as_secs_f64();
    let n_r: usize = regs.iter().map(|r| r.len()).sum();
    let t = Instant::now();
    for (s, r) in reads.iter().zip(regs.iter_mut()) {
        sort_dedup_patch(&idx, Some(s), r);
    }
    let t_dd = t.elapsed().as_secs_f64();
    let us = |x: f64| x * 1e6 / reads.len() as f64;
    println!(
        "{} reads, per read: seed {:.1}us ({} intervals), chain(SA+tree) {:.1}us ({} chains), flt {:.1}us ({}), extend {:.1}us ({} regions), dedup {:.1}us",
        reads.len(),
        us(t_seed),
        n_mem,
        us(t_chain),
        n_c,
        us(t_flt),
        n_cf,
        us(t_ext),
        n_r,
        us(t_dd)
    );
}
