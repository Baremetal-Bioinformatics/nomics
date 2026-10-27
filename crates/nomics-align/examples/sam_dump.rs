//! Write SAM records (no header) for a chunk, as bwa would print them.
//!
//! usage: sam_dump <index-hint> <R1.fastq.gz> <R2.fastq.gz> <out-prefix>  (writes <out-prefix>.sam.txt)
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use nomics_align::index::Index;
use nomics_align::mem::{process_batch, read_batch};
use nomics_io::gzip::GzReader;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Index::load(Path::new(&a[1])).unwrap();
    let open = |p: &str| BufReader::with_capacity(1 << 20, GzReader::open(Path::new(p)).unwrap());
    let (mut r1, mut r2) = (open(&a[2]), open(&a[3]));
    let mut out = BufWriter::new(File::create(format!("{}.sam.txt", a[4])).unwrap());
    let mut first_pair = 0u64;
    loop {
        let reads = read_batch(&mut r1, &mut r2).unwrap();
        if reads.is_empty() {
            break;
        }
        let mut buf = Vec::new();
        let (_, log) = process_batch(
            &idx,
            b"HG002_WES.L1",
            &reads,
            first_pair,
            1,
            false,
            &mut buf,
        );
        for l in log {
            eprintln!("{l}");
        }
        out.write_all(&buf).unwrap();
        first_pair += (reads.len() / 2) as u64;
    }
}
