use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use nomics_hc::finalize::{HcRead, finalize_region};
use nomics_hc::trim::{TrimHap, TrimResult, trim_region};
use nomics_io::bam::{CIGAR_OPS, RecordSet, RecordView};

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

#[derive(Debug)]
struct AssemblyInput {
    contig: String,
    start: i64,
    end: i64,
    padded_ref_start: i64,
    padded_ref: Vec<u8>,
    haplotypes: Vec<TrimHap>,
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

fn parse_interval(value: &str) -> (String, i64, i64) {
    let (contig, coordinates) = value.rsplit_once(':').expect("interval contig:start-end");
    let (start, end) = coordinates.split_once('-').expect("interval start-end");
    (
        contig.to_owned(),
        start.parse().expect("interval start"),
        end.parse().expect("interval end"),
    )
}

fn parse_cigar(value: &str) -> Vec<u32> {
    let mut result = Vec::new();
    let mut length = 0u32;
    for byte in value.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
        } else {
            let operator = CIGAR_OPS
                .iter()
                .position(|&candidate| candidate == byte)
                .expect("known CIGAR operator") as u32;
            result.push(length << 4 | operator);
            length = 0;
        }
    }
    assert_eq!(length, 0, "CIGAR ends in a length");
    result
}

fn cigar_text(cigar: &[u32]) -> String {
    cigar
        .iter()
        .map(|&element| {
            format!(
                "{}{}",
                element >> 4,
                CIGAR_OPS[(element & 0xf) as usize] as char
            )
        })
        .collect()
}

fn load_assemblies(path: &Path) -> Vec<AssemblyInput> {
    let mut lines = BufReader::new(File::open(path).expect("open t0_assembled.tsv")).lines();
    let mut result = Vec::new();
    while let Some(header) = lines.next() {
        let header = header.expect("read assembly header");
        let fields: Vec<_> = header.split_whitespace().collect();
        assert_eq!(fields[0], "R");
        let contig = fields[1].to_owned();
        let start = fields[2].parse().expect("region start");
        let end = fields[3].parse().expect("region end");
        let haplotype_count: usize = fields[4].parse().expect("haplotype count");

        let reference = lines.next().expect("REF line").expect("read REF line");
        let reference_fields: Vec<_> = reference.split_whitespace().collect();
        assert_eq!(reference_fields[0], "REF");
        let (_, padded_ref_start, _) = parse_interval(reference_fields[1]);
        let padded_ref = reference_fields[2].as_bytes().to_vec();

        let mut haplotypes = Vec::with_capacity(haplotype_count);
        for _ in 0..haplotype_count {
            let line = lines.next().expect("H line").expect("read H line");
            let fields: Vec<_> = line.split_whitespace().collect();
            assert_eq!(fields[0], "H");
            let (hap_contig, hap_start, hap_end) = parse_interval(fields[4]);
            haplotypes.push(TrimHap {
                bases: fields[7].as_bytes().to_vec(),
                is_ref: fields[1].parse().expect("reference boolean"),
                alignment_start_hap_wrt_ref: fields[2].parse().expect("alignment start"),
                cigar: parse_cigar(fields[3]),
                contig: hap_contig,
                start: hap_start,
                end: hap_end,
            });
        }
        result.push(AssemblyInput {
            contig,
            start,
            end,
            padded_ref_start,
            padded_ref,
            haplotypes,
        });
    }
    result
}

fn load_expected(path: &Path) -> Vec<Vec<String>> {
    let mut regions = Vec::new();
    for line in BufReader::new(File::open(path).expect("open t0_trimmed.tsv")).lines() {
        let line = line.expect("read trim oracle line");
        if line.starts_with("R ") {
            regions.push(vec![line]);
        } else {
            regions.last_mut().expect("oracle starts with R").push(line);
        }
    }
    regions
}

fn output_lines(contig: &str, start: i64, end: i64, result: TrimResult) -> Vec<String> {
    let mut lines = vec![format!("R {contig} {start} {end}")];
    let TrimResult::Trimmed {
        span,
        padded_span,
        variation_present,
        haplotypes,
        reads,
        filtered,
    } = result
    else {
        lines.push("NOVAR".to_owned());
        return lines;
    };
    lines.push(format!(
        "SPAN {contig}:{}-{} {contig}:{}-{} {variation_present}",
        span.0, span.1, padded_span.0, padded_span.1
    ));
    for haplotype in haplotypes {
        lines.push(format!(
            "H {} {} {} {}:{}-{} {}",
            haplotype.is_ref,
            haplotype.alignment_start_hap_wrt_ref,
            cigar_text(&haplotype.cigar),
            haplotype.contig,
            haplotype.start,
            haplotype.end,
            String::from_utf8_lossy(&haplotype.bases)
        ));
    }
    for read in reads {
        lines.push(format!(
            "READ {} {} {} {}",
            String::from_utf8_lossy(&read.name),
            read.flag,
            read.start,
            cigar_text(&read.cigar)
        ));
    }
    for read in filtered {
        lines.push(format!(
            "FILTERED {} {}",
            String::from_utf8_lossy(&read.name),
            read.flag
        ));
    }
    lines
}

#[test]
#[ignore = "requires the local GATK trim oracle fixtures"]
fn all_tier0_regions_match_gatk_4_6_2_0() {
    let dir = PathBuf::from(oracle_dir());
    let bam_path = dir.join("HG002_WES.recal.bam");
    let regions_path = dir.join("t0_regions.tsv");
    let assembled_path = dir.join("t0_assembled.tsv");
    let oracle_path = dir.join("t0_trimmed.tsv");
    if !bam_path.exists()
        || !regions_path.exists()
        || !assembled_path.exists()
        || !oracle_path.exists()
    {
        eprintln!(
            "skipping: trim oracle files are absent under {}",
            dir.display()
        );
        return;
    }

    let regions = load_regions(&regions_path);
    let assemblies = load_assemblies(&assembled_path);
    let expected = load_expected(&oracle_path);
    assert_eq!(regions.len(), 1_110, "fixture region count changed");
    assert_eq!(assemblies.len(), regions.len());
    assert_eq!(expected.len(), regions.len());

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

    for (region_index, ((region, assembly), expected_lines)) in
        regions.iter().zip(&assemblies).zip(&expected).enumerate()
    {
        assert_eq!(assembly.contig.as_bytes(), region.contig);
        assert_eq!((assembly.start, assembly.end), (region.start, region.end));
        let &(contig_id, contig_len) = contigs.get(&region.contig).expect("known contig");
        let mut original_reads = Vec::new();
        for key in &region.members {
            let matches = reads_by_key.get(key).unwrap_or_else(|| {
                panic!(
                    "region {region_index} missing BAM member {}:{}",
                    String::from_utf8_lossy(&key.name),
                    key.flag
                )
            });
            original_reads.extend(
                matches
                    .iter()
                    .filter(|read| read.contig == contig_id)
                    .cloned(),
            );
        }
        let finalized = finalize_region(contig_len, region.start, region.end, &original_reads);
        let actual = trim_region(
            &assembly.contig,
            assembly.start,
            assembly.end,
            contig_len,
            assembly.padded_ref_start,
            &assembly.padded_ref,
            &assembly.haplotypes,
            &finalized,
        );
        let actual_lines = output_lines(&assembly.contig, assembly.start, assembly.end, actual);
        assert_eq!(
            actual_lines, *expected_lines,
            "first mismatch at region {region_index} {}:{}-{}",
            assembly.contig, assembly.start, assembly.end
        );
    }
}
