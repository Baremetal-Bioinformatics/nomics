//! BGZF block compression: throughput and size of htsjdk's level-2 replica
//! against the fast deflater, on real blocks of a BGZF file.
//!
//! `cargo run --release -p nomics-io --example deflate_bench -- <file.bam> [samples] [blocks] [threads]`
//!
//! Reads `blocks` consecutive members at each of `samples` evenly spaced
//! offsets (resynchronising on the BGZF header), then for each encoder:
//! compresses every payload on 1 and on `threads` threads, totals the member
//! sizes, inflates the result and checks it against the payload.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::time::Instant;

use nomics_io::bgzf::{Deflate, FOOTER_LEN, HEADER_LEN, decode_block, read_raw_block};
use nomics_io::deflate::fast::FastConfig;
use nomics_io::deflate::{FastDeflater, HtsjdkDeflater};
use nomics_io::inflate::inflate;

const MAGIC: [u8; 16] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0,
];

fn sample(path: &str, samples: u64, blocks: usize) -> (Vec<Vec<u8>>, u64) {
    let mut f = File::open(path).unwrap();
    let len = f.metadata().unwrap().len();
    let (mut payloads, mut original) = (Vec::new(), 0u64);
    for s in 0..samples {
        let at = len * s / samples;
        f.seek(SeekFrom::Start(at)).unwrap();
        let mut probe = vec![0u8; 1 << 18];
        let got = f.read(&mut probe).unwrap();
        let Some(skip) = probe[..got].windows(16).position(|w| w == MAGIC) else {
            continue;
        };
        let mut coffset = at + skip as u64;
        f.seek(SeekFrom::Start(coffset)).unwrap();
        let mut r = BufReader::with_capacity(1 << 20, &f);
        let mut raw = Vec::new();
        let mut taken = 0;
        while taken < blocks && read_raw_block(&mut r, &mut raw, coffset).unwrap() {
            let mut p = Vec::new();
            let info = decode_block(&raw, coffset, &mut p).unwrap();
            coffset += info.csize as u64;
            if p.is_empty() {
                continue;
            }
            original += info.csize as u64;
            payloads.push(p);
            taken += 1;
        }
    }
    (payloads, original)
}

trait Encoder: Deflate + Send {
    /// Fresh per-block state where the encoder has any (as the parallel writer does).
    fn prepare(&mut self) {}
}
impl Encoder for HtsjdkDeflater {
    fn prepare(&mut self) {
        self.reinit();
    }
}
impl Encoder for FastDeflater {}

fn compress_all<E: Encoder>(enc: &mut E, payloads: &[Vec<u8>]) -> Vec<Vec<u8>> {
    payloads
        .iter()
        .map(|p| {
            let mut out = Vec::new();
            enc.prepare();
            enc.deflate(p, &mut out);
            out
        })
        .collect()
}

fn run<E: Encoder>(
    name: &str,
    make: impl Fn() -> E + Sync,
    payloads: &[Vec<u8>],
    threads: usize,
    base: Option<u64>,
) -> u64 {
    let bytes: u64 = payloads.iter().map(|p| p.len() as u64).sum();
    let t = Instant::now();
    let outs = compress_all(&mut make(), payloads);
    let one = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let per = payloads.len().div_ceil(threads);
    std::thread::scope(|s| {
        for chunk in payloads.chunks(per) {
            let make = &make;
            s.spawn(move || compress_all(&mut make(), chunk));
        }
    });
    let many = t.elapsed().as_secs_f64();
    let size: u64 = outs
        .iter()
        .map(|o| (o.len() + HEADER_LEN + FOOTER_LEN) as u64)
        .sum();
    let t = Instant::now();
    let mut back = Vec::with_capacity(65536);
    for (o, p) in outs.iter().zip(payloads) {
        back.clear();
        inflate(o, &mut back, p.len()).unwrap();
        assert!(back == *p, "{name}: roundtrip mismatch");
    }
    let inf = t.elapsed().as_secs_f64();
    let mb = bytes as f64 / 1e6;
    let rel = base.map_or(String::from("1.000"), |b| {
        format!("{:.3}", size as f64 / b as f64)
    });
    println!(
        "{name:<34} {:>9.1} {:>10.1} {:>13} {:>6} {:>9.1} {:>8.3}",
        mb / one,
        mb / many,
        size,
        rel,
        mb / inf,
        bytes as f64 / size as f64
    );
    size
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: deflate_bench <file.bam> [samples] [blocks] [threads]");
    let samples = args.get(2).map_or(16, |s| s.parse().unwrap());
    let blocks = args.get(3).map_or(256, |s| s.parse().unwrap());
    let threads = args.get(4).map_or(20, |s| s.parse().unwrap());
    let (payloads, original) = sample(path, samples, blocks);
    let bytes: u64 = payloads.iter().map(|p| p.len() as u64).sum();
    println!(
        "{} blocks, {:.1} MB uncompressed, {:.1} MB as stored in the file; {threads} threads",
        payloads.len(),
        bytes as f64 / 1e6,
        original as f64 / 1e6
    );
    for (name, hw) in [
        ("crc32 slice-by-8", false),
        ("crc32 dispatch (hw if present)", true),
    ] {
        let t = Instant::now();
        let mut x = 0u32;
        for p in &payloads {
            let c = nomics_io::crc32::Crc32::new();
            x ^= if hw {
                c.update(p)
            } else {
                c.update_portable(p)
            }
            .finish();
        }
        let s = t.elapsed().as_secs_f64();
        println!(
            "{name}: {:.0} MB/s (1 thread) [{x:08x}]",
            bytes as f64 / 1e6 / s
        );
    }
    println!(
        "{:<34} {:>9} {:>10} {:>13} {:>6} {:>9} {:>8}",
        "encoder", "MB/s(1)", "MB/s(N)", "bytes", "size", "infl MB/s", "ratio"
    );
    let base = run(
        "htsjdk level 2 (exact)",
        HtsjdkDeflater::new,
        &payloads,
        threads,
        None,
    );
    let fast = |max_insert, block_symbols| {
        move || {
            FastDeflater::with_config(FastConfig {
                max_insert,
                block_symbols,
            })
        }
    };
    let variants: Vec<(usize, usize)> = match std::env::var("FAST_VARIANTS") {
        // "ins,blk;ins,blk;..."
        Ok(v) => v
            .split(';')
            .map(|t| {
                let f: Vec<usize> = t.split(',').map(|x| x.parse().unwrap()).collect();
                (f[0], f[1])
            })
            .collect(),
        Err(_) => vec![
            (258, 16383),
            (32, 16383),
            (16, 16383),
            (8, 16383),
            (0, 16383),
            (258, 4096),
            (258, 65000),
        ],
    };
    for (ins, bs) in variants {
        let name = format!("fast ins<={ins} blk={bs}");
        run(&name, fast(ins, bs), &payloads, threads, Some(base));
    }
}
