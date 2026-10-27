//! Stage-1 gate: parse a HaplotypeCaller .vcf.gz into typed records,
//! re-encode header + records with nomics-hc's writer, and compare the text
//! and the compressed bytes with the original.
//!
//! usage: vcf_roundtrip <golden.vcf.gz> <out.vcf.gz>
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

use nomics_hc::vcf::{Header, parse_record, write_vcf_gz};
use nomics_io::gzip::GzReader;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut text = String::new();
    GzReader::open(Path::new(&a[1]))
        .unwrap()
        .read_to_string(&mut text)
        .unwrap();
    let (mut meta, mut contigs, mut sample) = (Vec::new(), Vec::new(), String::new());
    let mut records = Vec::new();
    let mut index = HashMap::new();
    for line in text.lines() {
        if let Some(m) = line.strip_prefix("##") {
            if let Some(c) = m.strip_prefix("contig=<ID=") {
                let (name, rest) = c.split_once(",length=").unwrap();
                let len: u64 = rest.trim_end_matches('>').parse().unwrap();
                index.insert(name.to_string(), contigs.len());
                contigs.push((name.to_string(), len));
            } else if m != "fileformat=VCFv4.2" {
                meta.push(m.to_string());
            }
        } else if let Some(h) = line.strip_prefix("#CHROM") {
            sample = h.rsplit('\t').next().unwrap().to_string();
        } else {
            records.push(parse_record(line, &index).unwrap());
        }
    }
    let header = Header {
        meta,
        contigs,
        sample,
    };
    let out = write_vcf_gz(Vec::new(), &header, &records, 20).unwrap();
    std::fs::write(&a[2], &out).unwrap();
    let mut ours = String::new();
    GzReader::open(Path::new(&a[2]))
        .unwrap()
        .read_to_string(&mut ours)
        .unwrap();
    let golden_bytes = std::fs::read(&a[1]).unwrap();
    let first_diff = ours.lines().zip(text.lines()).position(|(x, y)| x != y);
    println!(
        "{} records; text {}; compressed bytes {} ({} vs {} bytes)",
        records.len(),
        if ours == text {
            "IDENTICAL".to_string()
        } else {
            format!("DIFFER at line {first_diff:?}")
        },
        if out == golden_bytes {
            "IDENTICAL"
        } else {
            "DIFFER"
        },
        out.len(),
        golden_bytes.len()
    );
    if let Some(i) = first_diff {
        println!(
            "ours:   {}\ngolden: {}",
            ours.lines().nth(i).unwrap(),
            text.lines().nth(i).unwrap()
        );
    }
}
