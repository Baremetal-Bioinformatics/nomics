//! Release CLI byte gates. Inputs/goldens under NOMICS_RUNS are read-only.
use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{Instant, SystemTime, UNIX_EPOCH},
};

#[test]
fn six_release_bam_stats_gates() {
    let runs = std::env::var_os("NOMICS_RUNS")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(std::env::var_os("HOME").unwrap()).join("omics/runs"));
    let binary = std::env::var_os("NOMICS_QC_BIN")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/release/nomics")
        });
    let reference = PathBuf::from(std::env::var_os("HOME").unwrap())
        .join("omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta");
    let dir = std::env::temp_dir().join(format!(
        "nomics-bam-stats-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&dir).unwrap();
    let mut compared = 0;
    let mut failures = Vec::new();
    for (run, sample, sizes) in [
        (
            "tier0-wes-20261008-r1",
            "HG002_WES",
            [50_729_515, 80_694_756],
        ),
        (
            "wgs20-20261010-r1",
            "HG002_WGS20",
            [1_091_442_954, 1_477_781_384],
        ),
        (
            "tier1-wes-20261008-r1",
            "HG002_WES",
            [5_047_709_513, 3_516_306_871],
        ),
    ] {
        let results = runs.join(run).join("results");
        for (i, (stage, subdir)) in [("md", "markduplicates"), ("recal", "recalibrated")]
            .iter()
            .enumerate()
        {
            let name = format!("{sample}.{stage}.bam");
            let bam = results.join(format!("preprocessing/{subdir}/{sample}/{name}"));
            let golden = results.join(format!("reports/samtools/{sample}/{name}.stats"));
            if !bam.is_file() || !golden.is_file() {
                eprintln!("{run} {stage}: BAM/golden absent, skipped");
                continue;
            }
            assert_eq!(
                fs::metadata(&bam).unwrap().len(),
                sizes[i],
                "{run} {stage}: incomplete or changed fixture"
            );
            assert!(
                binary.is_file(),
                "build release CLI first or set NOMICS_QC_BIN"
            );
            let out = dir.join(format!("{run}-{stage}.stats"));
            let start = Instant::now();
            let status = Command::new(&binary)
                .args(["qc", "bam-stats", "--bam"])
                .arg(&bam)
                .arg("--reference")
                .arg(&reference)
                .arg("--out")
                .arg(&out)
                .status()
                .unwrap();
            let elapsed = start.elapsed();
            assert!(status.success(), "{run} {stage}: CLI failed");
            compared += 1;
            let ours = fs::read(&out).unwrap();
            let theirs = fs::read(&golden).unwrap();
            if ours == theirs {
                assert!(
                    Command::new("cmp")
                        .arg(&out)
                        .arg(&golden)
                        .status()
                        .unwrap()
                        .success()
                );
                eprintln!(
                    "{run} {stage}: identical (cmp), {:.3}s",
                    elapsed.as_secs_f64()
                );
            } else {
                let a = String::from_utf8_lossy(&ours);
                let b = String::from_utf8_lossy(&theirs);
                let line = a
                    .split_inclusive('\n')
                    .zip(b.split_inclusive('\n'))
                    .position(|(x, y)| x != y)
                    .unwrap_or_else(|| a.lines().count().min(b.lines().count()));
                let msg = format!(
                    "{run} {stage}: {:.3}s, first differing line {}: nomics {:?}, golden {:?}",
                    elapsed.as_secs_f64(),
                    line + 1,
                    a.lines().nth(line),
                    b.lines().nth(line)
                );
                eprintln!("{msg}");
                failures.push(msg);
            }
        }
    }
    fs::remove_dir_all(dir).unwrap();
    eprintln!(
        "BAM stats gates: {}/{} passed",
        compared - failures.len(),
        compared
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
