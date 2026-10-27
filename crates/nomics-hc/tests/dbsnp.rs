//! dbSNP ID/DB gate: for every record of a golden HaplotypeCaller VCF,
//! recompute the ID from dbSNP and compare ID and the DB flag.
//!   NOMICS_GOLDEN_VCF=<golden.vcf.gz> NOMICS_DBSNP=<dbsnp.vcf.gz> \
//!     cargo test -p nomics-hc --release --test dbsnp -- --nocapture
//! Skipped (passes) when the variables are unset: dbSNP lives on the Spark.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use nomics_hc::dbsnp::Dbsnp;
use nomics_io::bgzf::{decode_block, read_raw_block};

#[test]
fn golden_ids_match() {
    let (Some(golden), Some(dbsnp)) = (
        std::env::var_os("NOMICS_GOLDEN_VCF").map(PathBuf::from),
        std::env::var_os("NOMICS_DBSNP").map(PathBuf::from),
    ) else {
        eprintln!("NOMICS_GOLDEN_VCF / NOMICS_DBSNP unset: skipped");
        return;
    };
    let mut r = BufReader::new(File::open(golden).unwrap());
    let (mut raw, mut text, mut coff) = (Vec::new(), Vec::new(), 0u64);
    while read_raw_block(&mut r, &mut raw, coff).unwrap() {
        decode_block(&raw, coff, &mut text).unwrap();
        coff += raw.len() as u64;
    }
    let mut db = Dbsnp::open(&dbsnp).unwrap();
    let (mut n, mut with_id, mut bad) = (0, 0, 0);
    for line in text.split(|&b| b == b'\n') {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let f: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        let pos: i64 = std::str::from_utf8(f[1]).unwrap().parse().unwrap();
        let mut alleles: Vec<&[u8]> = vec![f[3]];
        alleles.extend(f[4].split(|&b| b == b','));
        let ours = db.id_for(f[0], pos, &alleles).unwrap();
        let golden_id = (f[2] != b".").then(|| String::from_utf8_lossy(f[2]).into_owned());
        let golden_db = f[7].split(|&b| b == b';').any(|k| k == b"DB");
        n += 1;
        with_id += golden_id.is_some() as usize;
        if ours != golden_id || golden_db != golden_id.is_some() {
            bad += 1;
            if bad <= 10 {
                eprintln!(
                    "{}:{} {} {}: ours {:?} golden {:?} DB {}",
                    String::from_utf8_lossy(f[0]),
                    pos,
                    String::from_utf8_lossy(f[3]),
                    String::from_utf8_lossy(f[4]),
                    ours,
                    golden_id,
                    golden_db
                );
            }
        }
    }
    eprintln!("dbSNP: {n} records, {with_id} with ID, {bad} differ");
    assert!(n > 0);
    assert_eq!(bad, 0);
}
