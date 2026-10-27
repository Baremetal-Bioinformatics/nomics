//! LoglessPairHMM port vs the reference jar (reference/oracle/hc/GatkPairHmm.java).
use std::io::{BufRead, BufReader};

use nomics_hc::pairhmm::PairHmm;
use nomics_io::gzip::GzReader;

fn quals(s: &str) -> Vec<u8> {
    s.split(',').map(|x| x.parse().unwrap()).collect()
}

#[test]
fn matches_gatk_jar() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/pairhmm.tsv.gz");
    let r = BufReader::new(GzReader::open(std::path::Path::new(path)).unwrap());
    let mut hmm = PairHmm::default();
    let (mut n, mut bad) = (0, 0);
    for line in r.lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split(' ').collect();
        let got = hmm.log10_likelihood(
            f[0].as_bytes(),
            f[1].as_bytes(),
            &quals(f[2]),
            &quals(f[3]),
            &quals(f[4]),
            &quals(f[5]),
        );
        let want = u64::from_str_radix(f[6], 16).unwrap();
        n += 1;
        if got.to_bits() != want {
            bad += 1;
            if bad <= 5 {
                eprintln!(
                    "case {n}: got {got:e} want {:e} (hap {} read {})",
                    f64::from_bits(want),
                    f[0].len(),
                    f[1].len()
                );
            }
        }
    }
    assert_eq!(bad, 0, "{bad} of {n} differ");
    assert_eq!(n, 20006);
}
