//! Build the release CLI first: cargo build --release -p nomics.
//! NOMICS_RUNS defaults to ~/omics/runs. NOMICS_QC_BIN can select a release
//! binary in another target directory. No tool or fixture is modified.
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn eight_release_cli_gates() {
    let runs = std::env::var_os("NOMICS_RUNS")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(std::env::var_os("HOME").unwrap()).join("omics/runs"));
    let binary = std::env::var_os("NOMICS_QC_BIN")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/release/nomics")
        });
    let mut failures = Vec::new();
    let mut compared = 0;
    for tier in ["tier0", "tier1"] {
        let results = runs.join(format!("{tier}-wes-20261008-r1/results"));
        let prefix = "HG002_WES.haplotypecaller";
        let vcf = results.join(format!(
            "variant_calling/haplotypecaller/HG002_WES/{prefix}.vcf.gz"
        ));
        let vc = results.join("reports/vcftools/haplotypecaller/HG002_WES");
        let bc = results.join("reports/bcftools/haplotypecaller/HG002_WES");
        let suffixes = [
            "TsTv.count",
            "TsTv.qual",
            "FILTER.summary",
            "bcftools_stats.txt",
        ];
        let goldens: Vec<_> = suffixes
            .iter()
            .enumerate()
            .map(|(i, s)| if i == 3 { &bc } else { &vc }.join(format!("{prefix}.{s}")))
            .collect();
        if !vcf.is_file() || goldens.iter().any(|p| !p.is_file()) {
            eprintln!(
                "{tier}: goldens absent at {}: skipped four QC gates",
                results.display()
            );
            continue;
        }
        assert!(
            binary.is_file(),
            "build the release CLI first: cargo build --release -p nomics (or set NOMICS_QC_BIN)"
        );
        let dir = std::env::temp_dir().join(format!(
            "nomics-qc-{}-{}-{tier}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let status = Command::new(&binary)
            .arg("qc")
            .arg("vcf-stats")
            .arg("--vcf")
            .arg(&vcf)
            .arg("--outdir")
            .arg(&dir)
            .status()
            .unwrap();
        assert!(status.success(), "{tier}: CLI failed");
        for (suffix, golden) in suffixes.iter().zip(goldens) {
            let ours = fs::read(dir.join(format!("{prefix}.{suffix}"))).unwrap();
            let theirs = fs::read(golden).unwrap();
            compared += 1;
            let result = if *suffix == "TsTv.qual" {
                nomics_qc::compare_tstv_qual(
                    std::str::from_utf8(&ours).unwrap(),
                    std::str::from_utf8(&theirs).unwrap(),
                )
                .map(|fields| {
                    eprintln!(
                        "{tier} {suffix}: identical under normalization, {} exempt fields",
                        fields.len()
                    );
                    for f in fields {
                        eprintln!(
                            "  line {} col {}: nomics {:?}, golden {:?}{}",
                            f.row,
                            f.column,
                            f.nomics,
                            f.golden,
                            if f.nomics == f.golden { " (equal)" } else { "" }
                        );
                    }
                })
            } else if ours == theirs {
                // Also run the requested literal cmp, rather than only equality.
                let status = Command::new("cmp")
                    .arg(dir.join(format!("{prefix}.{suffix}")))
                    .arg(
                        if *suffix == "bcftools_stats.txt" {
                            &bc
                        } else {
                            &vc
                        }
                        .join(format!("{prefix}.{suffix}")),
                    )
                    .status()
                    .unwrap();
                assert!(status.success());
                eprintln!("{tier} {suffix}: identical (cmp)");
                Ok(())
            } else {
                let a = String::from_utf8_lossy(&ours);
                let b = String::from_utf8_lossy(&theirs);
                let line = a
                    .split_inclusive('\n')
                    .zip(b.split_inclusive('\n'))
                    .position(|(x, y)| x != y)
                    .unwrap_or_else(|| a.lines().count().min(b.lines().count()));
                Err(format!(
                    "first differing line {}: nomics {:?}, golden {:?}",
                    line + 1,
                    a.lines().nth(line),
                    b.lines().nth(line)
                ))
            };
            if let Err(e) = result {
                eprintln!("{tier} {suffix}: {e}");
                failures.push(format!("{tier} {suffix}: {e}"));
            }
        }
        fs::remove_dir_all(dir).unwrap();
    }
    eprintln!(
        "QC gates: {}/{} passed",
        compared - failures.len(),
        compared
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
