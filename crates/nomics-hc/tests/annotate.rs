use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::process::Command;

use nomics_hc::annotate::{
    AnnotationCall, AnnotationGenotype, AnnotationRead, AnnotationResult, annotate,
};
use nomics_hc::genotype::{GenotypingEngine, genotype_likelihoods, pls};
use nomics_hc::jrandom::JavaRandom;
use nomics_hc::vcf::{InfoValue, format_vcf_double};

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn default_fixture() -> String {
    format!("{}/hc_finalize/t0_calls.tsv.gz", fixtures_dir())
}

#[derive(Clone)]
struct VcfLine {
    contig: String,
    pos: i64,
    alleles: Vec<Vec<u8>>,
    qual: f64,
    info: BTreeMap<String, String>,
    gt: [usize; 2],
    ad: Vec<i32>,
    dp: i32,
    gq: i32,
    pl: Vec<i32>,
}

fn parse_vcf(line: &str) -> VcfLine {
    let fields: Vec<&str> = line.split('\t').collect();
    assert_eq!(fields.len(), 10, "bad VCF line: {line}");
    let mut alleles = vec![fields[3].as_bytes().to_vec()];
    alleles.extend(fields[4].split(',').map(|x| x.as_bytes().to_vec()));
    let info = fields[7]
        .split(';')
        .map(|field| field.split_once('=').unwrap_or((field, "")))
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect();
    let keys: Vec<&str> = fields[8].split(':').collect();
    let vals: Vec<&str> = fields[9].split(':').collect();
    let format: HashMap<&str, &str> = keys.into_iter().zip(vals).collect();
    let gt: Vec<usize> = format["GT"]
        .split('/')
        .map(|x| x.parse().unwrap())
        .collect();
    VcfLine {
        contig: fields[0].to_string(),
        pos: fields[1].parse().unwrap(),
        alleles,
        qual: fields[5].parse().unwrap(),
        info,
        gt: [gt[0], gt[1]],
        ad: format["AD"]
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect(),
        dp: format["DP"].parse().unwrap(),
        gq: format["GQ"].parse().unwrap(),
        pl: format["PL"]
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect(),
    }
}

fn parse_cigar(text: &str) -> Vec<u32> {
    let mut result = Vec::new();
    let mut length = 0_u32;
    for b in text.bytes() {
        if b.is_ascii_digit() {
            length = length * 10 + u32::from(b - b'0');
            continue;
        }
        let op = match b {
            b'M' => 0,
            b'I' => 1,
            b'D' => 2,
            b'N' => 3,
            b'S' => 4,
            b'H' => 5,
            b'P' => 6,
            b'=' => 7,
            b'X' => 8,
            _ => panic!("bad CIGAR operator {b}"),
        };
        result.push((length << 4) | op);
        length = 0;
    }
    assert_eq!(length, 0);
    result
}

fn parse_read(line: &str) -> AnnotationRead {
    let f: Vec<&str> = line.split_ascii_whitespace().collect();
    assert_eq!(f[0], "E");
    AnnotationRead {
        flag: f[2].parse().unwrap(),
        mapq: f[3].parse().unwrap(),
        reverse: f[4].parse().unwrap(),
        start: f[5].parse().unwrap(),
        cigar: parse_cigar(f[6]),
        bases: f[7].as_bytes().to_vec(),
        quals: f[8].split(',').map(|x| x.parse().unwrap()).collect(),
        likelihoods: f[9..]
            .iter()
            .map(|x| f64::from_bits(u64::from_str_radix(x, 16).unwrap()))
            .collect(),
    }
}

fn rendered_info(result: &AnnotationResult) -> BTreeMap<String, String> {
    result
        .info
        .iter()
        .map(|(key, value)| {
            let value = match value {
                InfoValue::Flag => String::new(),
                InfoValue::Str(x) => x.clone(),
                InfoValue::Ints(x) => x
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                InfoValue::Doubles(x) => x
                    .iter()
                    .map(|&v| format_vcf_double(v))
                    .collect::<Vec<_>>()
                    .join(","),
            };
            (key.clone(), value)
        })
        .collect()
}

fn info_string(info: &BTreeMap<String, String>) -> String {
    info.iter()
        .map(|(key, value)| {
            if value.is_empty() {
                key.clone()
            } else {
                format!("{key}={value}")
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
#[ignore = "requires the local GATK tier0 annotation fixture (NOMICS_HC_FIXTURES)"]
fn tier0_gatk_default_annotations() {
    let fixture = std::env::var_os("NOMICS_HC_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_fixture()));
    let output = Command::new("gzip")
        .args(["-dc"])
        .arg(&fixture)
        .output()
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", fixture.display()));
    assert!(
        output.status.success(),
        "gzip failed for {}",
        fixture.display()
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = text.lines().collect();

    let mut computed = Vec::new();
    let mut expected = Vec::new();
    let mut genotyper = GenotypingEngine::new();
    let mut contigs: HashMap<String, usize> = HashMap::new();
    let mut i = 0;
    while i < lines.len() {
        if let Some(vc_text) = lines[i].strip_prefix("VC ") {
            let vc = parse_vcf(vc_text);
            i += 1;
            let matrix_alleles: Vec<Vec<u8>> = lines[i]
                .strip_prefix("A ")
                .unwrap()
                .split_ascii_whitespace()
                .map(|x| x.as_bytes().to_vec())
                .collect();
            i += 1;
            let mut reads = Vec::new();
            while i < lines.len() && lines[i].starts_with("E ") {
                reads.push(parse_read(lines[i]));
                i += 1;
            }
            let matrix: Vec<Vec<f64>> = (0..matrix_alleles.len())
                .map(|a| reads.iter().map(|r| r.likelihoods[a]).collect())
                .collect();
            let input_pls = pls(&genotype_likelihoods(&matrix));
            let next_contig = contigs.len();
            let contig_index = *contigs.entry(vc.contig.clone()).or_insert(next_contig);
            let exact_call = genotyper
                .calculate(
                    contig_index,
                    vc.pos,
                    vc.pos + vc.alleles[0].len() as i64 - 1,
                    &matrix_alleles,
                    &input_pls,
                )
                .expect("oracle VC was not reproduced by genotyper");
            let kept_alleles: Vec<&Vec<u8>> = exact_call
                .alleles_kept
                .iter()
                .map(|&a| &matrix_alleles[a])
                .collect();
            assert_eq!(kept_alleles, vc.alleles.iter().collect::<Vec<_>>());
            let mle_ac = vc.info["MLEAC"]
                .split(',')
                .map(|x| x.parse().unwrap())
                .collect();
            let mle_af = vc.info["MLEAF"]
                .split(',')
                .map(|x| x.parse().unwrap())
                .collect();
            let call = AnnotationCall {
                end: vc.pos + vc.alleles[0].len() as i64 - 1,
                alleles: vc.alleles.clone(),
                matrix_alleles,
                log10_p_error: exact_call.log10_p_error,
                genotype: AnnotationGenotype {
                    alleles: vc.gt,
                    pl: vc.pl,
                    gq: vc.gq,
                    mle_ac,
                    mle_af,
                },
                contig: vc.contig,
                pos: vc.pos,
                reads,
            };
            let mut rng = JavaRandom::gatk();
            computed.push((annotate(&call, &mut rng), vc.ad, vc.dp));
            continue;
        }
        if let Some(call_text) = lines[i].strip_prefix("CALL ") {
            expected.push(parse_vcf(call_text));
        }
        i += 1;
    }

    assert_eq!(computed.len(), 1_047);
    assert_eq!(expected.len(), computed.len());
    let mut mismatches: BTreeMap<String, usize> = BTreeMap::new();
    let mut info_mismatches = 0_usize;
    for ((actual, vc_ad, vc_dp), expected) in computed.iter().zip(&expected) {
        assert_eq!(
            &actual.ad, vc_ad,
            "VC AD at {}:{}",
            expected.contig, expected.pos
        );
        assert_eq!(
            actual.dp, *vc_dp,
            "VC DP at {}:{}",
            expected.contig, expected.pos
        );
        assert_eq!(
            actual.ad, expected.ad,
            "CALL AD at {}:{}",
            expected.contig, expected.pos
        );
        assert_eq!(
            actual.dp, expected.dp,
            "CALL DP at {}:{}",
            expected.contig, expected.pos
        );
        let actual_info = rendered_info(actual);
        let mut expected_info = expected.info.clone();
        expected_info.remove("DB");
        let keys: std::collections::BTreeSet<&String> =
            actual_info.keys().chain(expected_info.keys()).collect();
        for key in keys {
            mismatches.entry(key.clone()).or_default();
            if actual_info.get(key) != expected_info.get(key) {
                *mismatches.get_mut(key).unwrap() += 1;
                eprintln!(
                    "{}:{} {key}: actual={:?} expected={:?} QUAL={} DP={}",
                    expected.contig,
                    expected.pos,
                    actual_info.get(key),
                    expected_info.get(key),
                    expected.qual,
                    expected.dp
                );
            }
        }
        if info_string(&actual_info) != info_string(&expected_info) {
            info_mismatches += 1;
        }
    }
    eprintln!("per-key mismatches: {mismatches:?}; INFO strings: {info_mismatches}");
    assert!(
        mismatches.values().all(|&count| count == 0),
        "annotation mismatches: {mismatches:?}"
    );
    assert_eq!(info_mismatches, 0);
}
