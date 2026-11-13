//! Check ApplyBQSR's output against its input: same records in the same order
//! (a subset), identical except for base qualities; qualities < 6 unchanged.
//!
//! `cargo run --release -p nomics-io --example recal_vs_md -- <md.bam> <recal.bam>`

use std::fs::File;
use std::io::BufReader;

use nomics_io::bam::{Header, RecordView, read_raw_record};
use nomics_io::bgzf::BgzfReader;

fn qual_range(raw: &[u8]) -> std::ops::Range<usize> {
    let l_seq = u32::from_le_bytes(raw[16..20].try_into().unwrap()) as usize;
    let n_cigar = u16::from_le_bytes([raw[12], raw[13]]) as usize;
    let start = 32 + raw[8] as usize + 4 * n_cigar + l_seq.div_ceil(2);
    start..start + l_seq
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let open = |p: &str| BgzfReader::new(BufReader::with_capacity(1 << 20, File::open(p).unwrap()));
    let (mut md, mut rc) = (open(&a[0]), open(&a[1]));
    let (hm, hr) = (
        Header::read(&mut md).unwrap(),
        Header::read(&mut rc).unwrap(),
    );
    let (tm, tr) = (
        String::from_utf8_lossy(&hm.text).into_owned(),
        String::from_utf8_lossy(&hr.text).into_owned(),
    );
    let extra: Vec<&str> = tr
        .lines()
        .filter(|l| !tm.lines().any(|m| m == *l))
        .collect();
    println!("header lines only in recal: {}", extra.len());
    for l in &extra {
        let fields: Vec<&str> = l.split('\t').map(|f| &f[..f.len().min(60)]).collect();
        println!("  {}", fields.join(" | "));
    }
    let (mut rm, mut rr) = (Vec::new(), Vec::new());
    let (mut n, mut skipped, mut other_diff, mut qual_changed, mut low_changed, mut low_bases) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    while read_raw_record(&mut rc, &mut rr).unwrap() {
        n += 1;
        loop {
            assert!(
                read_raw_record(&mut md, &mut rm).unwrap(),
                "recal record {n} not found in md.bam order"
            );
            let (vm, vr) = (RecordView::new(&rm).unwrap(), RecordView::new(&rr).unwrap());
            if vm.read_name() == vr.read_name()
                && vm.flag() == vr.flag()
                && vm.pos() == vr.pos()
                && vm.ref_id() == vr.ref_id()
            {
                break;
            }
            skipped += 1;
        }
        let q = qual_range(&rr);
        if rm.len() != rr.len() || rm[..q.start] != rr[..q.start] || rm[q.end..] != rr[q.end..] {
            other_diff += 1;
        }
        if rm[q.clone()] != rr[q.clone()] {
            qual_changed += 1;
        }
        for (x, y) in rm[q.clone()].iter().zip(&rr[q]) {
            if *x < 6 {
                low_bases += 1;
                if x != y {
                    low_changed += 1;
                }
            }
        }
    }
    let mut tail = 0u64;
    while read_raw_record(&mut md, &mut rm).unwrap() {
        tail += 1;
    }
    println!("recal records {n}; md records skipped between them {skipped}, after last {tail}");
    println!(
        "records differing outside QUAL: {other_diff}; records with QUAL changed: {qual_changed}"
    );
    println!("bases with QUAL < 6: {low_bases}, changed: {low_changed}");
}
