//! The light spill codec against BGZF on a real BAM: sizes, and one-core
//! compress and decode rates.
//!
//! `cargo run --release -p nomics-io --example spill_bench -- <file.bam> [MiB]`
//!
//! Decodes the BAM's stream (at most `MiB` of it, default all), then on one
//! thread: inflates the BAM's own BGZF members, re-compresses the stream
//! with the fast BAM deflater, and compresses and decodes it as spill
//! blocks, checking the round trip.

use std::time::Instant;

use nomics_io::bgzf::{BamCompression, ParallelBgzfWriter, decode_block, split_blocks};
use nomics_io::spill::{BLOCK, compress, decompress};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let cap = args
        .get(2)
        .map_or(usize::MAX, |m| m.parse::<usize>().unwrap() << 20);
    let file = std::fs::read(path).unwrap();
    let blocks = split_blocks(&file).unwrap();

    let t = Instant::now();
    let mut stream = Vec::new();
    let mut bgzf_bytes = 0;
    for r in &blocks {
        if stream.len() >= cap {
            break;
        }
        decode_block(&file[r.clone()], r.start as u64, &mut stream).unwrap();
        bgzf_bytes = r.end;
    }
    let inflate_s = t.elapsed().as_secs_f64();
    let raw = stream.len();
    let mib = raw as f64 / (1 << 20) as f64;

    let t = Instant::now();
    let mut w = ParallelBgzfWriter::with_compression(Vec::new(), 1, BamCompression::Fast);
    std::io::Write::write_all(&mut w, &stream).unwrap();
    let fast = w.finish().unwrap();
    let fast_s = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let mut table = Vec::new();
    let mut frames: Vec<(usize, Vec<u8>)> = Vec::new();
    for b in stream.chunks(BLOCK) {
        let mut c = Vec::new();
        compress(b, &mut table, &mut c);
        frames.push((b.len(), c));
    }
    let compress_s = t.elapsed().as_secs_f64();
    let spill: usize = frames.iter().map(|(_, c)| c.len() + 12).sum();

    let t = Instant::now();
    let mut back = Vec::with_capacity(raw);
    for (n, c) in &frames {
        decompress(c, *n, &mut back).unwrap();
    }
    let decode_s = t.elapsed().as_secs_f64();
    assert!(back == stream, "round trip");
    let crc_t = Instant::now();
    let crc: u32 = stream
        .chunks(BLOCK)
        .map(nomics_io::crc32::crc32)
        .fold(0, |a, c| a ^ c);
    let crc_s = crc_t.elapsed().as_secs_f64();

    println!(
        "{path}: {mib:.1} MiB uncompressed ({} BGZF members)",
        blocks.len()
    );
    println!(
        "  BGZF as written  {:>8.1} MiB  ratio {:.2}  inflate {:.0} MiB/s",
        bgzf_bytes as f64 / (1 << 20) as f64,
        raw as f64 / bgzf_bytes as f64,
        mib / inflate_s
    );
    println!(
        "  fast BGZF        {:>8.1} MiB  ratio {:.2}  deflate {:.0} MiB/s",
        fast.len() as f64 / (1 << 20) as f64,
        raw as f64 / fast.len() as f64,
        mib / fast_s
    );
    println!(
        "  spill            {:>8.1} MiB  ratio {:.2}  compress {:.0} MiB/s  decode {:.0} MiB/s  (+crc32 {:.0} MiB/s, {crc:08x})",
        spill as f64 / (1 << 20) as f64,
        raw as f64 / spill as f64,
        mib / compress_s,
        mib / decode_s,
        mib / crc_s
    );
}
