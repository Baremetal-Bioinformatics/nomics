//! Read realignment vs GATK 4.6.2.0 `GatkLikelihoods.java` tier0 output.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use nomics_hc::finalize::HcRead;
use nomics_hc::realign::realign_reads_to_best_haplotype;
use nomics_hc::trim::TrimHap;
use nomics_io::bam::CIGAR_OPS;
use nomics_io::gzip::GzReader;

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn default_fixtures() -> String {
    format!("{}/hc_finalize", fixtures_dir())
}
const EXPECTED_REALIGNMENTS: usize = 81_894;

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Key {
    name: Vec<u8>,
    flag: u16,
}

#[derive(Clone, Debug)]
struct ReadAlignment {
    key: Key,
    start: i64,
    cigar: Vec<u32>,
}

#[derive(Clone, Debug)]
struct InputRead {
    key: Key,
    mapq: u8,
    bases: Vec<u8>,
    quals: Vec<u8>,
}

#[derive(Clone, Debug)]
struct LkRow {
    key: Key,
    values: Vec<f64>,
}

#[derive(Default)]
struct Region {
    header: String,
    padded_reference_start: Option<i64>,
    haplotypes: Vec<TrimHap>,
    reads: Vec<ReadAlignment>,
    inputs: Vec<InputRead>,
    likelihoods: Vec<LkRow>,
    expected: Vec<ReadAlignment>,
}

fn fixture_dir() -> PathBuf {
    std::env::var_os("NOMICS_HC_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_fixtures()))
}

fn parse_key(fields: &[&str]) -> Key {
    Key {
        name: fields[1].as_bytes().to_vec(),
        flag: fields[2].parse().expect("numeric SAM flag"),
    }
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
    let mut length = 0_u32;
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
                CIGAR_OPS[(element & 15) as usize] as char
            )
        })
        .collect()
}

fn check_region(
    region: &Region,
    total: &mut usize,
    mismatch_count: &mut usize,
    mismatches: &mut Vec<String>,
) {
    if region.expected.is_empty() {
        assert!(region.likelihoods.is_empty());
        return;
    }
    assert_eq!(region.expected.len(), region.likelihoods.len());
    let input_by_key: HashMap<_, _> = region
        .inputs
        .iter()
        .enumerate()
        .map(|(index, input)| (input.key.clone(), index))
        .collect();
    let alignment_by_key: HashMap<_, _> = region
        .reads
        .iter()
        .enumerate()
        .map(|(index, read)| (read.key.clone(), index))
        .collect();
    assert_eq!(input_by_key.len(), region.inputs.len(), "duplicate IN key");
    assert_eq!(
        alignment_by_key.len(),
        region.reads.len(),
        "duplicate READ key"
    );

    let mut reads = Vec::with_capacity(region.likelihoods.len());
    for row in &region.likelihoods {
        let input = &region.inputs[*input_by_key.get(&row.key).expect("LK read has IN line")];
        let alignment = &region.reads[*alignment_by_key
            .get(&row.key)
            .expect("LK read has READ line")];
        reads.push(HcRead {
            name: row.key.name.clone(),
            flag: row.key.flag,
            contig: 0,
            start: alignment.start,
            cigar: alignment.cigar.clone(),
            bases: input.bases.clone(),
            quals: input.quals.clone(),
            mapq: input.mapq,
            mate_contig: -1,
            mate_start: 0,
            tlen: 0,
            record: 0,
        });
    }
    let mut likelihoods = vec![Vec::with_capacity(reads.len()); region.haplotypes.len()];
    for row in &region.likelihoods {
        assert_eq!(row.values.len(), region.haplotypes.len());
        for (haplotype, &value) in row.values.iter().enumerate() {
            likelihoods[haplotype].push(value);
        }
    }
    let actual = realign_reads_to_best_haplotype(
        &reads,
        &region.haplotypes,
        region
            .padded_reference_start
            .expect("REALIGN region has SPAN"),
        &likelihoods,
    );
    for ((row, actual), expected) in region.likelihoods.iter().zip(&actual).zip(&region.expected) {
        assert_eq!(row.key, expected.key);
        *total += 1;
        if actual.start != expected.start || actual.cigar != expected.cigar {
            *mismatch_count += 1;
            if mismatches.len() < 20 {
                mismatches.push(format!(
                    "{} {}:{} ours {} {} GATK {} {}",
                    region.header,
                    String::from_utf8_lossy(&row.key.name),
                    row.key.flag,
                    actual.start,
                    cigar_text(&actual.cigar),
                    expected.start,
                    cigar_text(&expected.cigar),
                ));
            }
        }
    }
}

#[test]
#[ignore = "needs the GatkLikelihoods tier0 fixture under `$NOMICS_FIXTURES`"]
fn matches_every_gatk_realign_line() {
    let fixture = fixture_dir().join("t0_likelihoods.tsv.gz");
    let reader = BufReader::with_capacity(
        1024 * 1024,
        GzReader::open(Path::new(&fixture)).expect("open t0_likelihoods.tsv.gz"),
    );
    let mut region = Region::default();
    let mut total = 0_usize;
    let mut mismatch_count = 0_usize;
    let mut mismatch_examples = Vec::new();
    for line in reader.lines() {
        let line = line.expect("read likelihood fixture line");
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        match fields[0] {
            "R" => {
                check_region(
                    &region,
                    &mut total,
                    &mut mismatch_count,
                    &mut mismatch_examples,
                );
                region = Region {
                    header: line,
                    ..Region::default()
                };
            }
            "SPAN" => {
                region.padded_reference_start = Some(parse_interval(fields[2]).1);
            }
            "H" => {
                let (contig, start, end) = parse_interval(fields[4]);
                region.haplotypes.push(TrimHap {
                    bases: fields[5].as_bytes().to_vec(),
                    is_ref: fields[1].parse().expect("reference boolean"),
                    alignment_start_hap_wrt_ref: fields[2]
                        .parse()
                        .expect("haplotype alignment start"),
                    cigar: parse_cigar(fields[3]),
                    contig,
                    start,
                    end,
                });
            }
            "READ" => region.reads.push(ReadAlignment {
                key: parse_key(&fields),
                start: fields[3].parse().expect("read start"),
                cigar: parse_cigar(fields[4]),
            }),
            "IN" => region.inputs.push(InputRead {
                key: parse_key(&fields),
                mapq: fields[3].parse().expect("mapping quality"),
                bases: fields[4].as_bytes().to_vec(),
                quals: fields[5]
                    .split(',')
                    .map(|quality| quality.parse().expect("base quality"))
                    .collect(),
            }),
            "LK" => region.likelihoods.push(LkRow {
                key: parse_key(&fields),
                values: fields[3..]
                    .iter()
                    .map(|value| {
                        f64::from_bits(u64::from_str_radix(value, 16).expect("hex double"))
                    })
                    .collect(),
            }),
            "REALIGN" => region.expected.push(ReadAlignment {
                key: parse_key(&fields),
                start: fields[3].parse().expect("realigned start"),
                cigar: parse_cigar(fields[4]),
            }),
            "NOVAR" | "FILTERED" | "DISQ" => {}
            kind => panic!("unexpected fixture line kind {kind}"),
        }
    }
    check_region(
        &region,
        &mut total,
        &mut mismatch_count,
        &mut mismatch_examples,
    );

    for mismatch in &mismatch_examples {
        eprintln!("{mismatch}");
    }
    assert_eq!(total, EXPECTED_REALIGNMENTS);
    assert_eq!(mismatch_count, 0, "REALIGN mismatches (first 20 shown)");
    eprintln!("{total}/{EXPECTED_REALIGNMENTS} REALIGN lines matched");
}
