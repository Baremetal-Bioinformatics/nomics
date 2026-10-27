use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use nomics_hc::genotype::{GenotypingEngine, genotype_likelihoods, pls};

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn default_fixture() -> String {
    format!("{}/hc_finalize/genotype.tsv.gz", fixtures_dir())
}

#[derive(Default)]
struct Mismatches {
    counts: BTreeMap<&'static str, usize>,
    examples: Vec<String>,
}

impl Mismatches {
    fn check(
        &mut self,
        field: &'static str,
        case: usize,
        matches: bool,
        detail: impl FnOnce() -> String,
    ) {
        if !matches {
            *self.counts.entry(field).or_default() += 1;
            if self.examples.len() < 20 {
                self.examples
                    .push(format!("case {case} {field}: {}", detail()));
            }
        }
    }

    fn count(&self, field: &'static str) -> usize {
        self.counts.get(field).copied().unwrap_or(0)
    }
}

fn fixture_path() -> PathBuf {
    std::env::var_os("NOMICS_HC_GENOTYPE_FIXTURE")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_fixture()))
}

fn parse_hex_f64(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).expect("invalid f64 bit pattern"))
}

fn parse_i32_list(text: &str) -> Vec<i32> {
    let inner = text
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .expect("invalid integer list");
    if inner.is_empty() {
        Vec::new()
    } else {
        inner.split(',').map(|x| x.parse().unwrap()).collect()
    }
}

fn parse_f64_list(text: &str) -> Vec<f64> {
    let inner = text
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .expect("invalid double list");
    if inner.is_empty() {
        Vec::new()
    } else {
        inner.split(',').map(|x| x.parse().unwrap()).collect()
    }
}

fn next_line(lines: &mut impl Iterator<Item = std::io::Result<String>>) -> String {
    lines
        .next()
        .expect("fixture ended early")
        .expect("fixture I/O error")
}

fn open_fixture(
    path: &Path,
) -> (
    std::process::Child,
    impl Iterator<Item = std::io::Result<String>>,
) {
    let mut child = Command::new("gzip")
        .arg("-dc")
        .arg(path)
        .stdout(Stdio::piped())
        .spawn()
        .expect("could not run gzip");
    let stdout = child.stdout.take().expect("gzip stdout was not piped");
    (child, BufReader::new(stdout).lines())
}

#[test]
fn spanning_deletion_state_is_kept_across_calls() {
    let deletion_alleles = vec![b"ATC".to_vec(), b"A".to_vec()];
    let downstream_alleles = vec![b"A".to_vec(), b"C".to_vec(), b"*".to_vec()];

    let mut engine = GenotypingEngine::new();
    assert!(
        engine
            .calculate(0, 100, 102, &deletion_alleles, &[100, 50, 0])
            .is_some()
    );
    let covered = engine
        .calculate(
            0,
            101,
            101,
            &downstream_alleles,
            &[100, 100, 100, 100, 0, 100],
        )
        .expect("the downstream SNP should emit");
    assert_eq!(covered.alleles_kept, vec![0, 1, 2]);

    let mut fresh_engine = GenotypingEngine::new();
    let uncovered = fresh_engine
        .calculate(
            0,
            101,
            101,
            &downstream_alleles,
            &[100, 100, 100, 100, 0, 100],
        )
        .expect("the downstream SNP should emit without the spurious star");
    assert_eq!(uncovered.alleles_kept, vec![0, 1]);
}

#[test]
#[ignore = "full 20,000-case GATK oracle"]
fn gatk_genotype_oracle() {
    let path = fixture_path();
    if !path.exists() {
        eprintln!("genotype oracle absent; skipped: {}", path.display());
        return;
    }

    let (mut child, mut lines) = open_fixture(&path);
    let mut mismatches = Mismatches::default();
    let mut cases = 0usize;
    let mut calls = 0usize;

    while let Some(line) = lines.next() {
        let line = line.expect("fixture I/O error");
        if line.is_empty() {
            continue;
        }
        cases += 1;
        let mut case_fields = line.split_whitespace();
        assert_eq!(case_fields.next(), Some("CASE"));
        let allele_text = case_fields.next().expect("CASE alleles missing");
        let alleles: Vec<Vec<u8>> = allele_text
            .split(',')
            .map(|allele| allele.as_bytes().to_vec())
            .collect();
        let read_count: usize = case_fields
            .next()
            .expect("CASE read count missing")
            .parse()
            .unwrap();

        let mut matrix = vec![vec![0.0; read_count]; alleles.len()];
        for _ in 0..alleles.len() * read_count {
            let matrix_line = next_line(&mut lines);
            let mut fields = matrix_line.split_whitespace();
            assert_eq!(fields.next(), Some("M"));
            let allele: usize = fields.next().unwrap().parse().unwrap();
            let read: usize = fields.next().unwrap().parse().unwrap();
            matrix[allele][read] = parse_hex_f64(fields.next().unwrap());
        }

        let oracle_gl_line = next_line(&mut lines);
        let oracle_gl: Vec<f64> = oracle_gl_line
            .split_whitespace()
            .skip(1)
            .map(parse_hex_f64)
            .collect();
        let oracle_pl_line = next_line(&mut lines);
        let oracle_pl: Vec<i32> = oracle_pl_line
            .split_whitespace()
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();

        let raw_gl = genotype_likelihoods(&matrix);
        let computed_pl = pls(&raw_gl);
        let computed_stored_gl: Vec<f64> =
            computed_pl.iter().map(|&pl| pl as f64 / -10.0).collect();
        mismatches.check(
            "GL",
            cases,
            bits_eq(&computed_stored_gl, &oracle_gl),
            || {
                format!(
                    "got {:?}, expected {:?}",
                    bits(&computed_stored_gl),
                    bits(&oracle_gl)
                )
            },
        );
        mismatches.check("PL", cases, computed_pl == oracle_pl, || {
            format!("got {computed_pl:?}, expected {oracle_pl:?}")
        });

        let oracle_call = next_line(&mut lines);
        let mut engine = GenotypingEngine::new();
        let call = engine.calculate(
            0,
            1000,
            1000 + alleles[0].len() as i64 - 1,
            &alleles,
            &computed_pl,
        );
        if oracle_call == "NULL" {
            mismatches.check("NULL/CALL", cases, call.is_none(), || {
                format!("got {call:?}, expected NULL")
            });
            continue;
        }

        calls += 1;
        let call = match call {
            Some(call) => call,
            None => {
                mismatches.check("NULL/CALL", cases, false, || {
                    "got NULL, expected CALL".into()
                });
                continue;
            }
        };
        let mut fields = oracle_call.split_whitespace();
        assert_eq!(fields.next(), Some("CALL"));
        let expected_error = parse_hex_f64(fields.next().unwrap());
        let expected_kept: Vec<usize> = fields
            .next()
            .unwrap()
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect();
        let expected_gt: Vec<usize> = fields
            .next()
            .unwrap()
            .split('/')
            .map(|x| x.parse().unwrap())
            .collect();
        let expected_pl = parse_i32_list(fields.next().unwrap());
        let expected_gq: i32 = fields.next().unwrap().parse().unwrap();
        let expected_mleac = parse_i32_list(fields.next().unwrap());
        let expected_mleaf = parse_f64_list(fields.next().unwrap());

        mismatches.check(
            "log10PError",
            cases,
            call.log10_p_error.to_bits() == expected_error.to_bits(),
            || {
                format!(
                    "got {:016x}, expected {:016x}",
                    call.log10_p_error.to_bits(),
                    expected_error.to_bits()
                )
            },
        );
        mismatches.check("kept", cases, call.alleles_kept == expected_kept, || {
            format!("got {:?}, expected {expected_kept:?}", call.alleles_kept)
        });
        mismatches.check("GT", cases, call.gt.as_slice() == expected_gt, || {
            format!("got {:?}, expected {expected_gt:?}", call.gt)
        });
        mismatches.check("call PL", cases, call.pl == expected_pl, || {
            format!("got {:?}, expected {expected_pl:?}", call.pl)
        });
        mismatches.check("GQ", cases, call.gq == expected_gq, || {
            format!("got {}, expected {expected_gq}", call.gq)
        });
        mismatches.check("MLEAC", cases, call.mleac == expected_mleac, || {
            format!("got {:?}, expected {expected_mleac:?}", call.mleac)
        });
        mismatches.check(
            "MLEAF",
            cases,
            bits_eq(&call.mleaf, &expected_mleaf),
            || format!("got {:?}, expected {expected_mleaf:?}", call.mleaf),
        );
    }

    assert!(
        child.wait().expect("gzip wait failed").success(),
        "gzip failed"
    );
    let fields = [
        "GL",
        "PL",
        "NULL/CALL",
        "log10PError",
        "kept",
        "GT",
        "call PL",
        "GQ",
        "MLEAC",
        "MLEAF",
    ];
    eprintln!("genotype oracle: {cases} cases, {calls} calls");
    for field in fields {
        eprintln!("  {field}: {} mismatches", mismatches.count(field));
    }
    for example in &mismatches.examples {
        eprintln!("  {example}");
    }
    assert!(
        mismatches.counts.is_empty(),
        "GATK genotype oracle mismatches"
    );
    assert_eq!(cases, 20_000);
    assert_eq!(calls, 10_213);
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn bits_eq(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

#[test]
fn gvcf_retains_non_ref_outside_the_six_alt_limit() {
    // Golden chr20:46128304: six concrete alts plus NON_REF, the two middle
    // alleles win. NON_REF must not consume a slot in the concrete-alt limit.
    let alleles = [
        "CTTT",
        "C",
        "CT",
        "CTT",
        "CTTTT",
        "CTTTTT",
        "CTTTTTT",
        "<NON_REF>",
    ]
    .map(|a| a.as_bytes().to_vec());
    let pl = [
        177, 97, 558, 114, 423, 432, 76, 206, 233, 230, 70, 75, 64, 0, 211, 105, 195, 158, 62, 195,
        389, 180, 371, 300, 174, 248, 431, 614, 246, 406, 379, 266, 249, 357, 484, 527,
    ];
    let call = GenotypingEngine::gvcf()
        .calculate(0, 46128304, 46128307, &alleles, &pl)
        .unwrap();
    assert_eq!(call.alleles_kept, (0..8).collect::<Vec<_>>());
    assert_eq!(call.gt, [3, 4]);
    assert_eq!(call.pl, pl);
    assert_eq!(call.mleac, [0, 0, 1, 1, 0, 0, 0]);
}
