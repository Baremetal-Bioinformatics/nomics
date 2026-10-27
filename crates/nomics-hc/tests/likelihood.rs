//! Read likelihoods vs GATK's PairHMMLikelihoodCalculationEngine on tier0 regions
//! (reference/oracle/hc/GatkLikelihoods.java; fixture under `$NOMICS_FIXTURES`).
use std::io::{BufRead, BufReader};

use nomics_hc::likelihood::{LkRead, compute};
use nomics_io::gzip::GzReader;

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn fixture() -> String {
    format!("{}/hc_finalize/t0_likelihoods.tsv.gz", fixtures_dir())
}

#[derive(Default)]
struct Region {
    haps: Vec<Vec<u8>>,
    reads: Vec<(String, Vec<u8>, Vec<u8>, u8)>,
    lk: Vec<(String, Vec<u64>)>,
    disq: Vec<String>,
}

fn check(r: &Region, bad: &mut usize) {
    if r.reads.is_empty() {
        return;
    }
    let haps: Vec<&[u8]> = r.haps.iter().map(|h| h.as_slice()).collect();
    let reads: Vec<LkRead> = r
        .reads
        .iter()
        .map(|(_, b, q, m)| LkRead {
            bases: b,
            quals: q,
            mapq: *m,
        })
        .collect();
    let res = compute(&haps, &reads);
    let got_disq: Vec<String> = res
        .disqualified
        .iter()
        .map(|&i| r.reads[i].0.clone())
        .collect();
    let got_lk: Vec<(String, Vec<u64>)> = res
        .kept
        .iter()
        .enumerate()
        .map(|(k, &i)| {
            (
                r.reads[i].0.clone(),
                res.values.iter().map(|v| v[k].to_bits()).collect(),
            )
        })
        .collect();
    if got_disq != r.disq || got_lk != r.lk {
        *bad += 1;
        if *bad <= 3 {
            eprintln!("region differs: disq ours {:?} gatk {:?}", got_disq, r.disq);
            for (a, b) in got_lk.iter().zip(&r.lk) {
                if a != b {
                    eprintln!("  first lk diff {} ours {:?} gatk {:?}", a.0, a.1, b.1);
                    break;
                }
            }
        }
    }
}

#[test]
#[ignore = "needs the GatkLikelihoods fixture under `$NOMICS_FIXTURES`"]
fn matches_gatk_engine() {
    let rd = BufReader::new(GzReader::open(std::path::Path::new(&fixture())).unwrap());
    let (mut regions, mut bad) = (0, 0);
    let mut cur = Region::default();
    for line in rd.lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "R" => {
                check(&cur, &mut bad);
                if !cur.reads.is_empty() {
                    regions += 1;
                }
                cur = Region::default();
            }
            "H" => cur.haps.push(f[f.len() - 1].as_bytes().to_vec()),
            "IN" => cur.reads.push((
                format!("{}:{}", f[1], f[2]),
                f[4].as_bytes().to_vec(),
                f[5].split(',').map(|x| x.parse().unwrap()).collect(),
                f[3].parse().unwrap(),
            )),
            "LK" => cur.lk.push((
                format!("{}:{}", f[1], f[2]),
                f[3..]
                    .iter()
                    .map(|x| u64::from_str_radix(x, 16).unwrap())
                    .collect(),
            )),
            "DISQ" => cur.disq.push(format!("{}:{}", f[1], f[2])),
            _ => {}
        }
    }
    check(&cur, &mut bad);
    if !cur.reads.is_empty() {
        regions += 1;
    }
    eprintln!("{regions} regions with reads, {bad} differ");
    assert_eq!(bad, 0);
}
