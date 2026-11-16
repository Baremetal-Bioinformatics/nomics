//! PCR error model and tandem-repeat lengths vs the reference jar
//! (reference/oracle/hc/GatkReadPrep.java).
use std::io::{BufRead, BufReader};

use nomics_hc::readprep::{pcr_cache, tandem_repeat_length};
use nomics_io::gzip::GzReader;

#[test]
fn matches_gatk_jar() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/readprep.tsv.gz"
    );
    let r = BufReader::new(GzReader::open(std::path::Path::new(path)).unwrap());
    let (mut reads, mut bad) = (0, 0);
    for line in r.lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "C" => assert_eq!(
                pcr_cache()[f[1].parse::<usize>().unwrap()],
                f[2].parse::<u8>().unwrap(),
                "PCR cache {}",
                f[1]
            ),
            "T" => {
                reads += 1;
                let b = f[1].as_bytes();
                let want: Vec<usize> = f[2].split(',').map(|x| x.parse().unwrap()).collect();
                let got: Vec<usize> = (0..b.len()).map(|o| tandem_repeat_length(b, o)).collect();
                if got != want {
                    bad += 1;
                    if bad <= 3 {
                        eprintln!("{}\n got  {got:?}\n want {want:?}", f[1]);
                    }
                }
            }
            _ => panic!("bad line"),
        }
    }
    assert_eq!(bad, 0, "{bad} of {reads} reads differ");
    assert_eq!(reads, 30000);
}
