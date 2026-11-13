//! Size and speed of the CRAM block codecs on one raw data-series file
//! (e.g. a slice's quality scores), each output checked by decoding it.
//!
//! `cargo run --release -p nomics-io --example cram_codec_bench -- <raw-file>...`

use std::time::Instant;

use nomics_io::cram::rans;
use nomics_io::deflate::{FastDeflater, ZlibDeflater};
use nomics_io::inflate::inflate;

fn main() {
    for path in std::env::args().skip(1) {
        let data = std::fs::read(&path).unwrap();
        println!("{path}: {} bytes", data.len());
        let t = Instant::now();
        let r0 = rans::encode_o0(&data);
        let e0 = t.elapsed();
        assert_eq!(rans::decode(&r0).unwrap(), data);
        let t = Instant::now();
        let r1 = rans::encode_o1(&data);
        let e1 = t.elapsed();
        assert_eq!(rans::decode(&r1).unwrap(), data);
        println!(
            "  rans o0    {:>10}  {:>8.1} MB/s",
            r0.len(),
            mbps(data.len(), e0)
        );
        println!(
            "  rans o1    {:>10}  {:>8.1} MB/s",
            r1.len(),
            mbps(data.len(), e1)
        );
        let mut fast = FastDeflater::new();
        let mut out = Vec::new();
        let t = Instant::now();
        fast.compress(&data, &mut out);
        let ef = t.elapsed();
        check(&out, &data);
        println!(
            "  fast       {:>10}  {:>8.1} MB/s",
            out.len(),
            mbps(data.len(), ef)
        );
        for level in 1..=9 {
            let mut z = ZlibDeflater::new(level);
            let mut out = Vec::new();
            let t = Instant::now();
            z.compress(&data, &mut out);
            let ez = t.elapsed();
            check(&out, &data);
            println!(
                "  zlib -{level}   {:>10}  {:>8.1} MB/s",
                out.len(),
                mbps(data.len(), ez)
            );
        }
    }
}

fn mbps(n: usize, d: std::time::Duration) -> f64 {
    n as f64 / d.as_secs_f64() / 1e6
}

fn check(deflated: &[u8], data: &[u8]) {
    let mut back = Vec::new();
    inflate(deflated, &mut back, data.len() + 1).unwrap();
    assert_eq!(back, data);
}
