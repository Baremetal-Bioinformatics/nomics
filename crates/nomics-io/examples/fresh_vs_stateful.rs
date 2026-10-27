//! How often does a fresh deflater disagree with htsjdk's reused one?
//!
//! `cargo run --release -p nomics-io --example fresh_vs_stateful -- <htsjdk.bam>`
//!
//! htsjdk compresses every block with one reset-per-block deflater whose
//! window carries stale bytes across blocks. Parallel compression would need
//! independent (fresh) state per block. This counts the blocks where the two
//! disagree on the golden file.

use std::fs::File;
use std::io::BufReader;

use nomics_io::bgzf::{Deflate, FOOTER_LEN, HEADER_LEN, decode_block, read_raw_block};
use nomics_io::deflate::HtsjdkDeflater;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: fresh_vs_stateful <htsjdk.bam>");
    let mut r = BufReader::with_capacity(1 << 20, File::open(&path).unwrap());
    let (mut raw, mut payload, mut ours) = (Vec::new(), Vec::new(), Vec::new());
    let mut stateful = HtsjdkDeflater::new();
    let (mut n, mut coffset, mut stateful_bad, mut fresh_bad) = (0u64, 0u64, 0u64, 0u64);
    while read_raw_block(&mut r, &mut raw, coffset).unwrap() {
        payload.clear();
        let info = decode_block(&raw, coffset, &mut payload).unwrap();
        coffset += info.csize as u64;
        if payload.is_empty() {
            continue;
        }
        n += 1;
        let golden = &raw[HEADER_LEN..raw.len() - FOOTER_LEN];
        ours.clear();
        stateful.deflate(&payload, &mut ours);
        if ours != golden {
            stateful_bad += 1;
        }
        ours.clear();
        HtsjdkDeflater::new().deflate(&payload, &mut ours);
        if ours != golden {
            fresh_bad += 1;
            if fresh_bad <= 5 {
                eprintln!("block {n} at {}: fresh-state output differs", info.coffset);
            }
        }
    }
    println!("{n} blocks: stateful mismatches {stateful_bad}, fresh-state mismatches {fresh_bad}");
}
