//! Print `M` lines (oracle dump format) for a chunk's read pairs.
//!
//! usage: seed_dump <index-hint> <R1.fastq.gz> <R2.fastq.gz> <out-prefix>   (writes <out-prefix>.M.txt)
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use nomics_align::index::Index;
use nomics_align::seed::{MIN_SEED_LEN, SmemAux, collect_intv, nt4};
use nomics_io::gzip::GzReader;

fn seqs(p: &str) -> Vec<Vec<u8>> {
    let r = BufReader::with_capacity(1 << 20, GzReader::open(Path::new(p)).unwrap());
    r.lines()
        .skip(1)
        .step_by(4)
        .map(|l| l.unwrap().bytes().map(nt4).collect())
        .collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Index::load(Path::new(&a[1])).unwrap();
    let (r1, r2) = (seqs(&a[2]), seqs(&a[3]));
    let mut out = BufWriter::new(std::fs::File::create(format!("{}.M.txt", a[4])).unwrap());
    let mut aux = SmemAux::default();
    // NOMICS_GPU_SEED=1: intervals from the GPU kernel instead of the CPU.
    let gpu = std::env::var_os("NOMICS_GPU_SEED").map(|_| {
        let reads: Vec<nomics_align::sam::Read> = r1
            .iter()
            .zip(&r2)
            .flat_map(|(a, b)| [a, b])
            .map(|s| nomics_align::sam::Read {
                seq: s.clone(),
                ..Default::default()
            })
            .collect();
        nomics_align::mem::gpu_seeds(&idx, &reads)
    });
    for (pair, (s1, s2)) in r1.iter().zip(&r2).enumerate() {
        for (mate, s) in [s1, s2].into_iter().enumerate() {
            if (s.len() as i32) < MIN_SEED_LEN {
                continue;
            }
            match gpu.as_ref().and_then(|g| g.get(pair << 1 | mate)) {
                Some(m) => aux.mem = m.to_vec(),
                None => collect_intv(&idx.bwt, s, &mut aux),
            }
            write!(out, "M {}", pair << 1 | mate).unwrap();
            for p in &aux.mem {
                write!(
                    out,
                    " {},{},{},{},{}",
                    p.qbeg(),
                    p.qend(),
                    p.x[0],
                    p.x[1],
                    p.x[2]
                )
                .unwrap();
            }
            writeln!(out).unwrap();
        }
    }
}
