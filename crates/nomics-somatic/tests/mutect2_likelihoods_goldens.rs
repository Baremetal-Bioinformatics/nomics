//! Exact likelihood bits and ordered biological BAM records. No tolerances.
use nomics_somatic::{
    mutect2_likelihoods::{self, records},
    mutect2_regions::Inputs,
};
use std::{
    fs,
    path::{Path, PathBuf},
};
fn compare(a: &Path, b: &Path) -> Result<(), String> {
    let (a, b) = (fs::read(a).unwrap(), fs::read(b).unwrap());
    if a == b {
        return Ok(());
    }
    let at = a
        .iter()
        .zip(&b)
        .position(|(x, y)| x != y)
        .unwrap_or(a.len().min(b.len()));
    let line = a[..at].iter().filter(|&&b| b == b'\n').count() + 1;
    Err(format!(
        "first difference byte {at}, line {line}; lengths {} / {}",
        a.len(),
        b.len()
    ))
}
#[test]
fn mutect2_stage_three_available_goldens() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap());
    let root = std::env::var_os("NOMICS_M2_LIKELIHOODS_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic"));
    let reference = std::env::var_os("NOMICS_M2_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    let output =
        std::env::temp_dir().join(format!("nomics-m2-likelihoods-test-{}", std::process::id()));
    fs::create_dir_all(&output).unwrap();
    let mut failures = Vec::new();
    for c in ["chr21", "chr22", "chr19", "chr20", "chr13"] {
        let d = root.join(format!("m2-{c}"));
        if !d.join("oracle/haplotypes.bam").exists() {
            eprintln!(
                "skipping Mutect2 likelihoods {c}: {} absent; override NOMICS_M2_LIKELIHOODS_GOLDENS",
                d.display()
            );
            continue;
        }
        let dump = output.join(format!("{c}.likelihoods.tsv"));
        let records = output.join(format!("{c}.realigned.tsv"));
        let time = std::time::Instant::now();
        let (_, stats) = mutect2_likelihoods::run(
            &Inputs {
                input: &[
                    d.join(format!("HCC1395BL_WES.{c}.bam")),
                    d.join(format!("HCC1395_WES.{c}.bam")),
                ],
                normal_sample: "HCC1395_HCC1395BL_WES",
                reference: &reference,
                intervals: &d.join(format!("{c}.bed")),
                panel_of_normals: &d.join(format!("1000g_pon.{c}.vcf.gz")),
                germline_resource: &d.join(format!("af-only-gnomad.{c}.vcf.gz")),
                regions: &output.join(format!("{c}.igv")),
                stats: &output.join(format!("{c}.stats")),
                dump_prefix: None,
            },
            &dump,
            &records,
        )
        .unwrap();
        let golden = output.join(format!("{c}.golden.tsv"));
        let count = records::project(&d.join("oracle/haplotypes.bam"), &golden).unwrap();
        eprintln!(
            "{c}: {stats:?}; golden {count}; {:.3}s",
            time.elapsed().as_secs_f64()
        );
        if let Err(e) = compare(&records, &golden) {
            failures.push(format!("{c} BAM: {e}"));
        }
        if c == "chr21" {
            if let Some(path) = std::env::var_os("NOMICS_M2_LIKELIHOODS_ORACLE") {
                if let Err(e) = compare(&dump, Path::new(&path)) {
                    failures.push(format!("chr21 internal: {e}"));
                }
            } else {
                eprintln!(
                    "skipping chr21 internal: set NOMICS_M2_LIKELIHOODS_ORACLE to qualified Java dump"
                );
            }
        }
        if failures.is_empty() {
            for entry in fs::read_dir(&output).unwrap() {
                let path = entry.unwrap().path();
                if path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with(&format!("{c}."))
                {
                    fs::remove_file(path).unwrap();
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{}; retained {}",
        failures.join("; "),
        output.display()
    );
    fs::remove_dir_all(output).unwrap();
}
