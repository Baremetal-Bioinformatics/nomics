//! java.util.Random against a 200k-call stream from the reference JVM
//! (reference/oracle/hc/JavaRandomMathFixtures.java "random").
use std::io::{BufRead, BufReader};

use nomics_hc::jrandom::JavaRandom;
use nomics_io::gzip::GzReader;

#[test]
fn matches_reference_jvm_stream() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/java_random.tsv.gz"
    );
    let r = BufReader::new(GzReader::open(std::path::Path::new(path)).unwrap());
    let mut rng = JavaRandom::gatk();
    let mut n = 0;
    for line in r.lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "i" => {
                let bound: i32 = f[1].parse().unwrap();
                let want: i32 = f[2].parse().unwrap();
                assert_eq!(rng.next_int(bound), want, "line {n}: nextInt({bound})");
            }
            "g" => {
                let want = u64::from_str_radix(f[1], 16).unwrap();
                assert_eq!(
                    rng.next_gaussian().to_bits(),
                    want,
                    "line {n}: nextGaussian"
                );
            }
            "d" => {
                let want = u64::from_str_radix(f[1], 16).unwrap();
                assert_eq!(rng.next_double().to_bits(), want, "line {n}: nextDouble");
            }
            _ => panic!("bad line {line}"),
        }
        n += 1;
    }
    assert_eq!(n, 200_000);
}
