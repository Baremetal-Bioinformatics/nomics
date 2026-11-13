//! CLI byte gate, including exact argv provenance, against the optional Mac fixture.
#[cfg(unix)]
#[test]
fn featurecounts_cli_matches_all_five_samples() {
    use std::path::PathBuf;
    use std::process::Command;
    let configured = std::env::var_os("NOMICS_RNASEQ_FIXTURE");
    let fixture = configured.as_ref().map_or_else(
        || {
            PathBuf::from(
                std::env::var("NOMICS_FIXTURES")
                    .unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into()),
            )
            .join("rnaseq-test")
        },
        PathBuf::from,
    );
    if !fixture.exists() {
        assert!(
            configured.is_none(),
            "configured NOMICS_RNASEQ_FIXTURE does not exist"
        );
        eprintln!("SKIP: rnaseq fixture absent");
        return;
    }
    let fixture = std::fs::canonicalize(fixture).unwrap();
    let temp = std::env::temp_dir().join(format!("nomics-count-cli-{}", std::process::id()));
    std::fs::create_dir_all(&temp).unwrap();
    std::os::unix::fs::symlink(fixture.join("genome_gfp.gtf"), temp.join("genome_gfp.gtf"))
        .unwrap();
    for sample in [
        "RAP1_IAA_30M_REP1",
        "RAP1_UNINDUCED_REP1",
        "RAP1_UNINDUCED_REP2",
        "WT_REP1",
        "WT_REP2",
    ] {
        let bam = format!("{sample}.markdup.sorted.bam");
        let output = format!("{sample}.featureCounts.tsv");
        std::os::unix::fs::symlink(
            fixture.join("results/star_salmon").join(&bam),
            temp.join(&bam),
        )
        .unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_nomics"));
        cmd.current_dir(&temp).args([
            "featurecounts",
            "-B",
            "-C",
            "-g",
            "gene_biotype",
            "-t",
            "exon",
        ]);
        if !sample.starts_with("RAP1_UNINDUCED") {
            cmd.arg("-p");
        }
        let result = cmd
            .args([
                "-T",
                "6",
                "-a",
                "genome_gfp.gtf",
                "-s",
                "2",
                "-o",
                &output,
                &bam,
                "--mqc-prefix",
                sample,
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{sample}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        for name in [
            output.clone(),
            format!("{output}.summary"),
            format!("{sample}.biotype_counts_mqc.tsv"),
            format!("{sample}.biotype_counts_rrna_mqc.tsv"),
        ] {
            let ours = std::fs::read(temp.join(&name)).unwrap();
            let golden = std::fs::read(
                fixture
                    .join("results/star_salmon/featurecounts")
                    .join(&name),
            )
            .unwrap();
            if ours != golden {
                let left = String::from_utf8_lossy(&ours);
                let right = String::from_utf8_lossy(&golden);
                let diff = left
                    .lines()
                    .zip(right.lines())
                    .enumerate()
                    .find(|(_, (l, r))| l != r);
                panic!(
                    "{name}: first diff {diff:?}; lengths {} vs {}",
                    ours.len(),
                    golden.len()
                );
            }
        }
        eprintln!("{sample}: CLI 4/4 byte-identical (no normalization)");
    }
    std::fs::remove_dir_all(temp).unwrap();
}
