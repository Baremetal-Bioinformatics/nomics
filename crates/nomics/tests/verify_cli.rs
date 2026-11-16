//! `nomics verify` from the command line, without Docker or sarek: setup
//! faults give ERROR (exit 2) with a report that says why, never PASS, and
//! `verify check-report` catches an edited report.

use std::path::{Path, PathBuf};
use std::process::Command;

fn scratch(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nomics-verify-cli-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    for f in ["r1.fq.gz", "r2.fq.gz", "targets.bed", "manifest.tsv"] {
        std::fs::write(d.join(f), "x\n").unwrap();
    }
    std::fs::create_dir_all(d.join("ref")).unwrap();
    std::fs::create_dir_all(d.join("empty-path")).unwrap();
    d
}

fn nomics(d: &Path) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_nomics"));
    // No docker on PATH.
    c.env("PATH", d.join("empty-path"));
    c
}

fn verify_args(c: &mut Command, d: &Path, out: &str) {
    c.arg("verify")
        .arg("--fastq1")
        .arg(d.join("r1.fq.gz"))
        .arg("--fastq2")
        .arg(d.join("r2.fq.gz"))
        .arg("--reference")
        .arg(d.join("ref"))
        .arg("--targets")
        .arg(d.join("targets.bed"))
        .arg("--cache-dir")
        .arg(d.join("cache"))
        .arg("--outdir")
        .arg(d.join(out));
}

#[test]
fn missing_docker_is_an_error_with_a_sealed_report() {
    let d = scratch("docker");
    let mut c = nomics(&d);
    verify_args(&mut c, &d, "out");
    let status = c.status().unwrap();
    assert_eq!(status.code(), Some(2));
    let report = std::fs::read_to_string(d.join("out/verify.json")).unwrap();
    assert!(report.contains("\"verdict\": \"ERROR\""), "{report}");
    assert!(report.contains("Docker is required"), "{report}");
    let md = std::fs::read_to_string(d.join("out/verify.md")).unwrap();
    assert!(md.starts_with("# nomics verify: ERROR"));

    // The report is sealed; an edit shows.
    let check = |path: &Path| {
        Command::new(env!("CARGO_BIN_EXE_nomics"))
            .args(["verify", "check-report"])
            .arg(path)
            .status()
            .unwrap()
            .code()
    };
    assert_eq!(check(&d.join("out/verify.json")), Some(0));
    std::fs::write(d.join("edited.json"), report.replace("\"ERROR\"", "\"PASS\"")).unwrap();
    assert_eq!(check(&d.join("edited.json")), Some(1));

    // A non-empty output directory is refused before anything runs.
    let mut again = nomics(&d);
    verify_args(&mut again, &d, "out");
    assert_eq!(again.status().unwrap().code(), Some(2));
    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn an_unknown_manifest_is_refused() {
    let d = scratch("against");
    let mut c = nomics(&d);
    verify_args(&mut c, &d, "out");
    c.arg("--against").arg(d.join("manifest.tsv"));
    assert_eq!(c.status().unwrap().code(), Some(2));
    let report = std::fs::read_to_string(d.join("out/verify.json")).unwrap();
    assert!(report.contains("\"mode\": \"manifest\""));
    assert!(report.contains("is not a manifest verify knows"), "{report}");
    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn a_wrong_reference_is_refused_with_the_reason() {
    // `--against` checks the reference before running anything; a reference
    // directory with the wrong contigs is refused, not handled.
    let d = scratch("dict");
    let fasta = d.join("ref/Sequence/WholeGenomeFasta");
    std::fs::create_dir_all(&fasta).unwrap();
    std::fs::create_dir_all(d.join("ref/Annotation/GATKBundle")).unwrap();
    std::fs::write(fasta.join("Homo_sapiens_assembly38.fasta"), ">chr1\nACGT\n").unwrap();
    std::fs::write(fasta.join("Homo_sapiens_assembly38.dict"), "@SQ\tSN:chr1\tLN:4\tM5:x\n").unwrap();
    std::fs::write(d.join("ref/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"), "x").unwrap();
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../reference/tier0-wes-20261008-r1.digest.tsv");
    let mut c = nomics(&d);
    verify_args(&mut c, &d, "out");
    c.arg("--against").arg(&manifest).args(["--read-group-ds", "/x.fasta"]);
    assert_eq!(c.status().unwrap().code(), Some(2));
    let report = std::fs::read_to_string(d.join("out/verify.json")).unwrap();
    assert!(report.contains("refused: the reference .dict differs from GATK.GRCh38"), "{report}");
    std::fs::remove_dir_all(&d).unwrap();
}
