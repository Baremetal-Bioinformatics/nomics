//! Re-compress a BGZF file (BAM, VCF) with the htsjdk block profile.
//!
//! `cargo run --release -p nomics-io --example bgzf_recompress -- <in> <out> <htsjdk|fast> [threads]`
//!
//! The uncompressed stream is copied unchanged, so the output has the same
//! records and, for an htsjdk-written input, the same block payload sizes.
//! With `htsjdk` an htsjdk-written input is reproduced byte for byte.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use nomics_io::bgzf::{BamCompression, ParallelBgzfReader, ParallelBgzfWriter};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: bgzf_recompress <in> <out> <htsjdk|fast> [threads]");
        std::process::exit(2);
    }
    let mode = BamCompression::parse(&args[3]).expect("mode must be htsjdk or fast");
    let threads = args.get(4).map_or(8, |t| t.parse().unwrap());
    let t = Instant::now();
    let mut r = ParallelBgzfReader::new(File::open(&args[1]).unwrap(), threads);
    let out = BufWriter::with_capacity(1 << 22, File::create(&args[2]).unwrap());
    let mut w = ParallelBgzfWriter::with_compression(out, threads, mode);
    let n = std::io::copy(&mut r, &mut w).unwrap();
    w.finish().unwrap().flush().unwrap();
    let (a, b) = (
        std::fs::metadata(&args[1]).unwrap().len(),
        std::fs::metadata(&args[2]).unwrap().len(),
    );
    println!(
        "{n} bytes uncompressed; {a} -> {b} bytes ({:.4}x) in {:.1}s, {mode:?}, {threads} threads",
        b as f64 / a as f64,
        t.elapsed().as_secs_f64()
    );
}
