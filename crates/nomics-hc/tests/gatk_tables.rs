//! nomics-hc tables vs values dumped from the reference GATK jar
//! (reference/oracle/hc/GatkTables.java, run in the golden's container).
use std::io::{BufRead, BufReader};

use nomics_hc::tables::*;
use nomics_io::gzip::GzReader;

fn hex(s: &str) -> u64 {
    u64::from_str_radix(s, 16).unwrap()
}

#[test]
fn matches_gatk_jar() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/gatk_tables.tsv.gz"
    );
    let r = BufReader::new(GzReader::open(std::path::Path::new(path)).unwrap());
    let (mut n, mut bad) = (0, 0);
    for line in r.lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split(' ').collect();
        let (got, want) = match f[0] {
            "qualToProbLog10" => (qual_to_prob_log10(f[1].parse().unwrap()), hex(f[2])),
            "qualToErrorProbLog10" => (qual_to_error_prob_log10(f[1].parse().unwrap()), hex(f[2])),
            "qualToErrorProb" => (qual_to_error_prob(f[1].parse().unwrap()), hex(f[2])),
            "qualToProb" => (qual_to_prob(f[1].parse().unwrap()), hex(f[2])),
            "approxSum" => (
                approximate_log10_sum_log10(f64::from_bits(hex(f[1])), f64::from_bits(hex(f[2]))),
                hex(f[3]),
            ),
            "normalDistribution" => (
                normal_distribution(50.0, 17.0, f[3].parse().unwrap()),
                hex(f[4]),
            ),
            "LOG10_ONE_THIRD" => (log10_one_third(), hex(f[1])),
            "LOG10_E" => (log10_e(), hex(f[1])),
            _ => continue, // logGamma etc.: checked where they are used
        };
        n += 1;
        if got.to_bits() != want {
            bad += 1;
            if bad <= 10 {
                eprintln!(
                    "{line}: got {:x} ({got:e}) want {want:x} ({:e})",
                    got.to_bits(),
                    f64::from_bits(want)
                );
            }
        }
    }
    assert_eq!(bad, 0, "{bad} of {n} values differ");
    assert!(n > 280_000);
}
