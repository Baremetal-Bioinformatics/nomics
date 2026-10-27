use nomics_hc::finalize;

use finalize::{HcRead, finalize_region};
use nomics_io::bam::{CIGAR_OPS, RecordSet, RecordView};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn oracle_dir() -> String {
    format!("{}/hc_finalize", fixtures_dir())
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Key {
    name: Vec<u8>,
    flag: u16,
}

#[derive(Debug)]
struct Region {
    contig: Vec<u8>,
    start: i64,
    end: i64,
    members: Vec<Key>,
}

fn parse_key(value: &str) -> Key {
    let (name, flag) = value.rsplit_once(':').expect("region member has name:flag");
    Key {
        name: name.as_bytes().to_vec(),
        flag: flag.parse().expect("numeric SAM flag"),
    }
}

fn load_regions(path: &Path) -> Vec<Region> {
    BufReader::new(File::open(path).expect("open t0_regions.tsv"))
        .lines()
        .map(|line| {
            let line = line.expect("read region line");
            let fields: Vec<_> = line.split('\t').collect();
            Region {
                contig: fields[0].as_bytes().to_vec(),
                start: fields[1].parse().expect("region start"),
                end: fields[2].parse().expect("region end"),
                members: fields[3..].iter().map(|value| parse_key(value)).collect(),
            }
        })
        .collect()
}

fn cigar_text(cigar: &[u32]) -> String {
    cigar
        .iter()
        .map(|&c| {
            format!(
                "{}{}",
                c >> 4,
                CIGAR_OPS.get((c & 0xf) as usize).copied().unwrap_or(b'?') as char
            )
        })
        .collect()
}

fn output_line(read: &HcRead) -> String {
    let qualities = read
        .quals
        .iter()
        .map(u8::to_string)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{} {} {} {} {} {}",
        String::from_utf8_lossy(&read.name),
        read.flag,
        read.start,
        cigar_text(&read.cigar),
        String::from_utf8_lossy(&read.bases),
        qualities
    )
}

fn next_line(lines: &mut impl Iterator<Item = std::io::Result<String>>) -> String {
    lines
        .next()
        .expect("oracle ended early")
        .expect("read oracle line")
}

#[test]
#[ignore = "requires the local GATK finalizeRegion oracle fixture"]
fn all_tier0_regions_match_gatk_4_6_2_0() {
    let dir = PathBuf::from(oracle_dir());
    let bam_path = dir.join("HG002_WES.recal.bam");
    let regions_path = dir.join("t0_regions.tsv");
    let oracle_path = dir.join("t0_finalized.tsv");
    if !bam_path.exists() || !regions_path.exists() || !oracle_path.exists() {
        eprintln!(
            "skipping: finalizeRegion oracle files are absent under {}",
            dir.display()
        );
        return;
    }

    let regions = load_regions(&regions_path);
    assert_eq!(regions.len(), 1_110, "fixture region count changed");
    let wanted: HashSet<Key> = regions
        .iter()
        .flat_map(|region| region.members.iter().cloned())
        .collect();

    let bam = RecordSet::load(&bam_path, 4).expect("load oracle BAM");
    let mut reads_by_key: HashMap<Key, Vec<HcRead>> = HashMap::new();
    for record in 0..bam.len() {
        let view = RecordView::new(bam.record(record)).expect("valid BAM record");
        let key = Key {
            name: view.read_name().to_vec(),
            flag: view.flag(),
        };
        if wanted.contains(&key) {
            reads_by_key
                .entry(key)
                .or_default()
                .push(HcRead::from_record(&view, record as u32));
        }
    }

    let contigs: HashMap<Vec<u8>, (i32, i64)> = bam
        .header
        .refs
        .iter()
        .enumerate()
        .map(|(index, reference)| {
            (
                reference.name.clone(),
                (index as i32, i64::from(reference.len)),
            )
        })
        .collect();
    let mut oracle = BufReader::new(File::open(&oracle_path).expect("open oracle")).lines();

    for (region_index, region) in regions.iter().enumerate() {
        let &(contig, contig_len) = contigs.get(&region.contig).unwrap_or_else(|| {
            panic!("unknown contig {}", String::from_utf8_lossy(&region.contig))
        });
        let mut input = Vec::new();
        for key in &region.members {
            let matches = reads_by_key.get(key).unwrap_or_else(|| {
                panic!(
                    "region {} missing BAM member {}:{}",
                    region_index,
                    String::from_utf8_lossy(&key.name),
                    key.flag
                )
            });
            input.extend(matches.iter().filter(|read| read.contig == contig).cloned());
        }

        let actual = finalize_region(contig_len, region.start, region.end, &input);
        let header = next_line(&mut oracle);
        let header_fields: Vec<_> = header.split_whitespace().collect();
        assert_eq!(
            header_fields[0], "R",
            "oracle header at region {region_index}"
        );
        assert_eq!(header_fields[1].as_bytes(), region.contig);
        assert_eq!(header_fields[2].parse::<i64>().unwrap(), region.start);
        assert_eq!(header_fields[3].parse::<i64>().unwrap(), region.end);
        let expected_count: usize = header_fields[4].parse().expect("oracle read count");
        assert_eq!(
            actual.len(),
            expected_count,
            "read-count mismatch at region {region_index} {}:{}-{} ({} input reads)",
            String::from_utf8_lossy(&region.contig),
            region.start,
            region.end,
            input.len()
        );
        for (read_index, read) in actual.iter().enumerate() {
            let expected = next_line(&mut oracle);
            let got = output_line(read);
            assert_eq!(
                got,
                expected,
                "first mismatch at region {region_index} {}:{}-{}, read {read_index}, source BAM record {}",
                String::from_utf8_lossy(&region.contig),
                region.start,
                region.end,
                read.record
            );
        }
    }
    assert!(oracle.next().is_none(), "oracle has trailing lines");
}
