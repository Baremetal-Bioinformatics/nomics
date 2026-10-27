//! Tier0 gate for the GATK 4.6.2.0 per-region calling loop.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use nomics_hc::call::{CallRegionInput, EvidenceIndex, call_region};
use nomics_hc::finalize::HcRead;
use nomics_hc::genotype::{GenotypingEngine, NO_CALL_ALLELE};
use nomics_hc::trim::{TrimHap, TrimResult, trim_region};
use nomics_hc::vcf::{format_qual, format_vcf_double};
use nomics_io::bam::CIGAR_OPS;
use nomics_io::gzip::GzReader;

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn default_fixtures() -> String {
    format!("{}/hc_finalize", fixtures_dir())
}
const REGION_COUNT: usize = 1_110;
const CALL_COUNT: usize = 1_047;
const EVIDENCE_COUNT: usize = 71_732;

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Key {
    name: Vec<u8>,
    flag: u16,
}

#[derive(Clone, Debug)]
struct ReadData {
    key: Key,
    mapq: u8,
    start: i64,
    cigar: Vec<u32>,
    bases: Vec<u8>,
    quals: Vec<u8>,
}

impl ReadData {
    fn hc_read(&self) -> HcRead {
        HcRead {
            name: self.key.name.clone(),
            flag: self.key.flag,
            contig: 0,
            start: self.start,
            cigar: self.cigar.clone(),
            bases: self.bases.clone(),
            quals: self.quals.clone(),
            mapq: self.mapq,
            mate_contig: -1,
            mate_start: 0,
            tlen: 0,
            record: 0,
        }
    }
}

#[derive(Default)]
struct InputRegion {
    header: String,
    contig: String,
    active_span: Option<(i64, i64)>,
    padded_start: Option<i64>,
    reference: Vec<u8>,
    haplotypes: Vec<TrimHap>,
    read_alignments: HashMap<Key, (i64, Vec<u32>)>,
    inputs: HashMap<Key, (u8, Vec<u8>, Vec<u8>)>,
    realignments: HashMap<Key, (i64, Vec<u32>)>,
    likelihood_rows: Vec<(Key, Vec<f64>)>,
    filtered: Vec<Key>,
}

struct AssemblyData {
    padded_start: i64,
    reference: Vec<u8>,
    haplotypes: Vec<TrimHap>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct VcFields {
    pos: i64,
    reference: Vec<u8>,
    alternates: Vec<Vec<u8>>,
    qual: String,
    gt: String,
    pl: Vec<i32>,
    gq: i32,
    mleac: String,
    mleaf: String,
}

#[derive(Default)]
struct ExpectedEvent {
    vc: Option<VcFields>,
    alleles: Vec<Vec<u8>>,
    evidence: Vec<ReadData>,
    likelihoods: Vec<Vec<f64>>,
}

#[derive(Default)]
struct ExpectedRegion {
    header: String,
    events: Vec<ExpectedEvent>,
    calls: Vec<VcFields>,
}

fn fixture_dir() -> PathBuf {
    std::env::var_os("NOMICS_HC_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_fixtures()))
}

fn parse_interval(value: &str) -> (String, i64, i64) {
    let (contig, coordinates) = value.rsplit_once(':').expect("contig:start-end");
    let (start, end) = coordinates.split_once('-').expect("start-end");
    (
        contig.to_owned(),
        start.parse().expect("interval start"),
        end.parse().expect("interval end"),
    )
}

fn parse_key(fields: &[&str]) -> Key {
    Key {
        name: fields[1].as_bytes().to_vec(),
        flag: fields[2].parse().expect("SAM flag"),
    }
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
        .map(|element| {
            format!(
                "{}{}",
                element >> 4,
                CIGAR_OPS[(element & 15) as usize] as char
            )
        })
        .collect()
}

fn parse_quals(value: &str) -> Vec<u8> {
    if value.is_empty() {
        Vec::new()
    } else {
        value
            .split(',')
            .map(|quality| quality.parse().expect("base quality"))
            .collect()
    }
}

fn load_assemblies(path: &Path) -> HashMap<String, AssemblyData> {
    let mut result = HashMap::new();
    let mut header = String::new();
    let mut current: Option<AssemblyData> = None;
    for line in BufReader::new(File::open(path).expect("open t0_assembled.tsv")).lines() {
        let line = line.expect("read assembly fixture");
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        match fields[0] {
            "R" => {
                if let Some(assembly) = current.take() {
                    result.insert(header.clone(), assembly);
                }
                header = format!("R {} {} {}", fields[1], fields[2], fields[3]);
                current = Some(AssemblyData {
                    padded_start: 0,
                    reference: Vec::new(),
                    haplotypes: Vec::new(),
                });
            }
            "REF" => {
                let (_, padded_start, _) = parse_interval(fields[1]);
                let assembly = current.as_mut().expect("REF after R");
                assembly.padded_start = padded_start;
                assembly.reference = fields[2].as_bytes().to_vec();
            }
            "H" => {
                let (contig, start, end) = parse_interval(fields[4]);
                current
                    .as_mut()
                    .expect("H after R")
                    .haplotypes
                    .push(TrimHap {
                        bases: fields[7].as_bytes().to_vec(),
                        is_ref: fields[1].parse().expect("reference boolean"),
                        alignment_start_hap_wrt_ref: fields[2].parse().expect("alignment start"),
                        cigar: parse_cigar(fields[3]),
                        contig,
                        start,
                        end,
                    });
            }
            kind => panic!("unexpected assembly line {kind}"),
        }
    }
    if let Some(assembly) = current {
        result.insert(header, assembly);
    }
    result
}

fn load_finalized(path: &Path) -> HashMap<String, Vec<ReadData>> {
    let mut result = HashMap::new();
    let mut header = String::new();
    for line in BufReader::new(File::open(path).expect("open t0_finalized.tsv")).lines() {
        let line = line.expect("read finalized fixture");
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        if fields[0] == "R" {
            header = format!("R {} {} {}", fields[1], fields[2], fields[3]);
            result.insert(header.clone(), Vec::new());
        } else {
            let key = Key {
                name: fields[0].as_bytes().to_vec(),
                flag: fields[1].parse().expect("SAM flag"),
            };
            result
                .get_mut(&header)
                .expect("read after R")
                .push(ReadData {
                    key,
                    mapq: 0,
                    start: fields[2].parse().expect("read start"),
                    cigar: parse_cigar(fields[3]),
                    bases: fields[4].as_bytes().to_vec(),
                    quals: parse_quals(fields[5]),
                });
        }
    }
    result
}

fn load_inputs(path: &Path, assemblies: &HashMap<String, AssemblyData>) -> Vec<InputRegion> {
    let reader = BufReader::with_capacity(
        1 << 20,
        GzReader::open(path).expect("open t0_likelihoods.tsv.gz"),
    );
    let mut regions: Vec<InputRegion> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("read likelihood fixture");
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        match fields[0] {
            "R" => {
                let header = line.clone();
                regions.push(InputRegion {
                    contig: fields[1].to_owned(),
                    reference: assemblies
                        .get(&header)
                        .unwrap_or_else(|| panic!("missing assembly for {header}"))
                        .reference
                        .clone(),
                    header,
                    ..InputRegion::default()
                });
            }
            "NOVAR" | "DISQ" => {}
            "SPAN" => {
                let region = regions.last_mut().expect("SPAN after R");
                let (_, active_start, active_end) = parse_interval(fields[1]);
                let (_, padded_start, _) = parse_interval(fields[2]);
                region.active_span = Some((active_start, active_end));
                region.padded_start = Some(padded_start);
            }
            "H" => {
                let region = regions.last_mut().expect("H after R");
                let (contig, start, end) = parse_interval(fields[4]);
                region.haplotypes.push(TrimHap {
                    bases: fields[5].as_bytes().to_vec(),
                    is_ref: fields[1].parse().expect("reference boolean"),
                    alignment_start_hap_wrt_ref: fields[2].parse().expect("alignment start"),
                    cigar: parse_cigar(fields[3]),
                    contig,
                    start,
                    end,
                });
            }
            "READ" => {
                let region = regions.last_mut().expect("READ after R");
                region.read_alignments.insert(
                    parse_key(&fields),
                    (
                        fields[3].parse().expect("read start"),
                        parse_cigar(fields[4]),
                    ),
                );
            }
            "FILTERED" => regions
                .last_mut()
                .expect("FILTERED after R")
                .filtered
                .push(parse_key(&fields)),
            "IN" => {
                let region = regions.last_mut().expect("IN after R");
                region.inputs.insert(
                    parse_key(&fields),
                    (
                        fields[3].parse().expect("mapping quality"),
                        fields[4].as_bytes().to_vec(),
                        parse_quals(fields[5]),
                    ),
                );
            }
            "LK" => regions
                .last_mut()
                .expect("LK after R")
                .likelihood_rows
                .push((
                    parse_key(&fields),
                    fields[3..]
                        .iter()
                        .map(|value| {
                            f64::from_bits(u64::from_str_radix(value, 16).expect("hex likelihood"))
                        })
                        .collect(),
                )),
            "REALIGN" => {
                let region = regions.last_mut().expect("REALIGN after R");
                region.realignments.insert(
                    parse_key(&fields),
                    (
                        fields[3].parse().expect("realigned start"),
                        parse_cigar(fields[4]),
                    ),
                );
            }
            kind => panic!("unexpected likelihood line {kind}"),
        }
    }
    regions
}

fn info_value(info: &str, key: &str) -> String {
    info.split(';')
        .find_map(|field| field.strip_prefix(&format!("{key}=")))
        .unwrap_or_else(|| panic!("missing {key} in {info}"))
        .to_owned()
}

fn parse_vc(value: &str) -> VcFields {
    let fields: Vec<_> = value.split('\t').collect();
    let sample: Vec<_> = fields[9].split(':').collect();
    VcFields {
        pos: fields[1].parse().expect("VC position"),
        reference: fields[3].as_bytes().to_vec(),
        alternates: fields[4]
            .split(',')
            .map(|allele| allele.as_bytes().to_vec())
            .collect(),
        qual: fields[5].to_owned(),
        gt: sample[0].to_owned(),
        pl: sample[4]
            .split(',')
            .map(|value| value.parse().expect("PL"))
            .collect(),
        gq: sample[3].parse().expect("GQ"),
        mleac: info_value(fields[7], "MLEAC"),
        mleaf: info_value(fields[7], "MLEAF"),
    }
}

fn load_expected(path: &Path) -> Vec<ExpectedRegion> {
    let reader =
        BufReader::with_capacity(1 << 20, GzReader::open(path).expect("open t0_calls.tsv.gz"));
    let mut regions: Vec<ExpectedRegion> = Vec::new();
    let mut event: Option<ExpectedEvent> = None;
    for line in reader.lines() {
        let line = line.expect("read call fixture");
        let (kind, value) = line.split_once(' ').expect("fixture line kind");
        match kind {
            "R" => {
                if let Some(event) = event.take() {
                    regions.last_mut().expect("event region").events.push(event);
                }
                regions.push(ExpectedRegion {
                    header: line,
                    ..ExpectedRegion::default()
                });
            }
            "VC" => {
                if let Some(event) = event.replace(ExpectedEvent {
                    vc: Some(parse_vc(value)),
                    ..ExpectedEvent::default()
                }) {
                    regions.last_mut().expect("VC region").events.push(event);
                }
            }
            "A" => {
                event.as_mut().expect("A after VC").alleles = value
                    .split_ascii_whitespace()
                    .map(|allele| allele.as_bytes().to_vec())
                    .collect();
            }
            "E" => {
                let fields: Vec<_> = line.split_ascii_whitespace().collect();
                event.as_mut().expect("E after VC").evidence.push(ReadData {
                    key: parse_key(&fields),
                    mapq: fields[3].parse().expect("mapping quality"),
                    start: fields[5].parse().expect("read start"),
                    cigar: parse_cigar(fields[6]),
                    bases: fields[7].as_bytes().to_vec(),
                    quals: parse_quals(fields[8]),
                });
                event.as_mut().expect("E after VC").likelihoods.push(
                    fields[9..]
                        .iter()
                        .map(|value| {
                            f64::from_bits(u64::from_str_radix(value, 16).expect("hex likelihood"))
                        })
                        .collect(),
                );
            }
            "CALL" => regions
                .last_mut()
                .expect("CALL region")
                .calls
                .push(parse_vc(value)),
            kind => panic!("unexpected call line {kind}"),
        }
    }
    if let Some(event) = event {
        regions
            .last_mut()
            .expect("last event region")
            .events
            .push(event);
    }
    regions
}

fn actual_vc(event: &nomics_hc::call::CalledEvent, final_record: bool) -> VcFields {
    let (pos, alleles) = if final_record {
        (event.final_start, &event.final_alleles)
    } else {
        (event.start, &event.call_alleles)
    };
    let gt = event
        .call
        .gt
        .iter()
        .map(|&allele| {
            if allele == NO_CALL_ALLELE {
                ".".to_owned()
            } else {
                allele.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("/");
    VcFields {
        pos,
        reference: alleles[0].clone(),
        alternates: alleles[1..].to_vec(),
        qual: format_qual(-10.0 * event.call.log10_p_error + 0.0),
        gt,
        pl: event.call.pl.clone(),
        gq: event.call.gq.min(99),
        mleac: event
            .call
            .mleac
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(","),
        mleaf: event
            .call
            .mleaf
            .iter()
            .map(|&value| format_vcf_double(value))
            .collect::<Vec<_>>()
            .join(","),
    }
}

fn filtered_data(
    region: &InputRegion,
    finalized: &HashMap<String, Vec<ReadData>>,
    assembly: &AssemblyData,
) -> Vec<ReadData> {
    let finalized_reads: Vec<_> = finalized
        .get(&region.header)
        .unwrap_or_else(|| panic!("missing finalized region {}", region.header))
        .iter()
        .map(ReadData::hc_read)
        .collect();
    let header: Vec<_> = region.header.split_ascii_whitespace().collect();
    let region_start = header[2].parse().expect("region start");
    let region_end = header[3].parse().expect("region end");
    let TrimResult::Trimmed { filtered, .. } = trim_region(
        &region.contig,
        region_start,
        region_end,
        i64::MAX / 4,
        assembly.padded_start,
        &assembly.reference,
        &assembly.haplotypes,
        &finalized_reads,
    ) else {
        panic!("{} unexpectedly has no variation", region.header);
    };
    let clipped: HashMap<_, _> = filtered
        .into_iter()
        .map(|read| {
            (
                Key {
                    name: read.name.clone(),
                    flag: read.flag,
                },
                ReadData {
                    key: Key {
                        name: read.name,
                        flag: read.flag,
                    },
                    mapq: read.mapq,
                    start: read.start,
                    cigar: read.cigar,
                    bases: read.bases,
                    quals: read.quals,
                },
            )
        })
        .collect();
    region
        .filtered
        .iter()
        .map(|key| {
            clipped
                .get(key)
                .unwrap_or_else(|| panic!("{} missing finalized {:?}", region.header, key))
                .clone()
        })
        .collect()
}

#[test]
#[ignore = "requires the local GATK tier0 calling fixtures"]
fn all_tier0_regions_match_gatk_calling_loop() {
    let dir = fixture_dir();
    let required = [
        dir.join("t0_assembled.tsv"),
        dir.join("t0_finalized.tsv"),
        dir.join("t0_likelihoods.tsv.gz"),
        dir.join("t0_calls.tsv.gz"),
    ];
    if required.iter().any(|path| !path.exists()) {
        eprintln!(
            "skipping: calling fixtures are absent under {}",
            dir.display()
        );
        return;
    }
    let assemblies = load_assemblies(&required[0]);
    let finalized = load_finalized(&required[1]);
    let inputs = load_inputs(&required[2], &assemblies);
    let expected = load_expected(&required[3]);
    assert_eq!(inputs.len(), REGION_COUNT);
    assert_eq!(expected.len(), REGION_COUNT);

    let mut contigs = HashMap::<String, usize>::new();
    let mut engine = GenotypingEngine::new();
    let mut total_calls = 0_usize;
    let mut total_evidence = 0_usize;
    for (region_index, (input, expected)) in inputs.iter().zip(&expected).enumerate() {
        assert_eq!(
            input.header, expected.header,
            "region {region_index} header"
        );
        if input.active_span.is_none() {
            assert!(
                expected.events.is_empty(),
                "{} unexpectedly called",
                input.header
            );
            continue;
        }
        let next_contig = contigs.len();
        let contig = *contigs.entry(input.contig.clone()).or_insert(next_contig);
        let mut reads = Vec::with_capacity(input.likelihood_rows.len());
        let mut matrix =
            vec![Vec::with_capacity(input.likelihood_rows.len()); input.haplotypes.len()];
        for (key, values) in &input.likelihood_rows {
            let &(mapq, ref bases, ref quals) = input
                .inputs
                .get(key)
                .unwrap_or_else(|| panic!("{} missing IN {:?}", input.header, key));
            let (start, cigar) = input
                .realignments
                .get(key)
                .unwrap_or_else(|| panic!("{} missing REALIGN {:?}", input.header, key));
            reads.push(
                ReadData {
                    key: key.clone(),
                    mapq,
                    start: *start,
                    cigar: cigar.clone(),
                    bases: bases.clone(),
                    quals: quals.clone(),
                }
                .hc_read(),
            );
            for (haplotype, &value) in values.iter().enumerate() {
                matrix[haplotype].push(value);
            }
        }
        let filtered_data = filtered_data(
            input,
            &finalized,
            assemblies
                .get(&input.header)
                .expect("assembly for input region"),
        );
        let filtered: Vec<_> = filtered_data.iter().map(ReadData::hc_read).collect();
        let scores = vec![0.0; input.haplotypes.len()];
        let reference_haplotype = input
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.is_ref)
            .expect("reference haplotype");
        let full_reference_start = input.padded_start.expect("SPAN padded interval")
            - reference_haplotype.alignment_start_hap_wrt_ref;
        let actual = call_region(
            &mut engine,
            CallRegionInput {
                contig,
                active_span: input.active_span.expect("checked"),
                padded_reference_start: full_reference_start,
                padded_reference: &input.reference,
                haplotypes: &input.haplotypes,
                haplotype_scores: &scores,
                reads: &reads,
                likelihoods: &matrix,
                filtered_reads: &filtered,
            },
        );
        assert_eq!(
            actual.len(),
            expected.events.len(),
            "first mismatch at region {region_index} {}: call count",
            input.header
        );
        assert_eq!(
            actual.len(),
            expected.calls.len(),
            "{} CALL count",
            input.header
        );
        for (event_index, ((actual, expected), expected_call)) in actual
            .iter()
            .zip(&expected.events)
            .zip(&expected.calls)
            .enumerate()
        {
            let context = format!("region {region_index} event {event_index} {}", input.header);
            assert_eq!(
                actual_vc(actual, false),
                *expected.vc.as_ref().expect("VC line"),
                "{context}: VC"
            );
            assert_eq!(actual.annotation_alleles, expected.alleles, "{context}: A");
            assert_eq!(
                actual.annotation_rows.len(),
                expected.evidence.len(),
                "{context}: E count"
            );
            total_evidence += actual.annotation_rows.len();
            for (row_index, (row, expected_read)) in actual
                .annotation_rows
                .iter()
                .zip(&expected.evidence)
                .enumerate()
            {
                let read = match row.evidence {
                    EvidenceIndex::Read(index) => &reads[index],
                    EvidenceIndex::Filtered(index) => &filtered[index],
                };
                assert_eq!(
                    read.name, expected_read.key.name,
                    "{context}: E {row_index} name"
                );
                assert_eq!(
                    read.flag, expected_read.key.flag,
                    "{context}: E {row_index} flag"
                );
                assert_eq!(
                    read.start, expected_read.start,
                    "{context}: E {row_index} start"
                );
                assert_eq!(
                    read.cigar,
                    expected_read.cigar,
                    "{context}: E {row_index} cigar; ours {} GATK {}",
                    cigar_text(&read.cigar),
                    cigar_text(&expected_read.cigar)
                );
                assert_eq!(
                    row.likelihoods.len(),
                    expected.likelihoods[row_index].len(),
                    "{context}: E {row_index} allele count"
                );
                for (allele, (&ours, &gatk)) in row
                    .likelihoods
                    .iter()
                    .zip(&expected.likelihoods[row_index])
                    .enumerate()
                {
                    assert_eq!(
                        ours.to_bits(),
                        gatk.to_bits(),
                        "{context}: E {row_index} allele {allele}; ours {:x} GATK {:x}",
                        ours.to_bits(),
                        gatk.to_bits()
                    );
                }
            }
            let ours_call = actual_vc(actual, true);
            assert_eq!(ours_call.pos, expected_call.pos, "{context}: CALL POS");
            assert_eq!(
                ours_call.reference, expected_call.reference,
                "{context}: CALL REF"
            );
            assert_eq!(
                ours_call.alternates, expected_call.alternates,
                "{context}: CALL ALT"
            );
        }
        total_calls += actual.len();
    }
    assert_eq!(total_calls, CALL_COUNT);
    assert_eq!(total_evidence, EVIDENCE_COUNT);
    eprintln!(
        "{REGION_COUNT}/{REGION_COUNT} regions, {total_calls}/{CALL_COUNT} calls, and \
         {total_evidence}/{EVIDENCE_COUNT} evidence rows matched"
    );
}
