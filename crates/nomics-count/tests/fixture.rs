//! Optional external gate: runs automatically when NOMICS_RNASEQ_FIXTURE or the
//! documented Mac fixture exists. No reference executable is needed.
use nomics_count::{
    Annotation, CountOptions, count_bam, write_multiqc, write_summary, write_table,
};
use std::path::PathBuf;

fn identical(sample: &str, kind: &str, ours: &[u8], golden: &[u8]) {
    if ours == golden {
        return;
    }
    let ours = String::from_utf8_lossy(ours);
    let golden = String::from_utf8_lossy(golden);
    let left: Vec<_> = ours.split('\n').collect();
    let right: Vec<_> = golden.split('\n').collect();
    let at = (0..left.len().max(right.len()))
        .find(|&i| left.get(i) != right.get(i))
        .unwrap();
    panic!(
        "{sample} {kind}: first diff line {}\nours: {:?}\ngolden: {:?}",
        at + 1,
        left.get(at),
        right.get(at)
    );
}

#[test]
fn all_five_rnaseq_samples() {
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
        eprintln!("SKIP: rnaseq fixture absent; set NOMICS_RNASEQ_FIXTURE");
        return;
    }
    let annotation =
        Annotation::read(&fixture.join("genome_gfp.gtf"), "gene_biotype", "exon").unwrap();
    let temp = std::env::temp_dir().join(format!("nomics-count-fixture-{}", std::process::id()));
    std::fs::create_dir_all(&temp).unwrap();
    for sample in [
        "RAP1_IAA_30M_REP1",
        "RAP1_UNINDUCED_REP1",
        "RAP1_UNINDUCED_REP2",
        "WT_REP1",
        "WT_REP2",
    ] {
        let bam = format!("{sample}.markdup.sorted.bam");
        let output = format!("{sample}.featureCounts.tsv");
        let input = fixture.join("results/star_salmon").join(&bam);
        // The fixture includes two SE samples; their golden commands omit -p.
        let options = CountOptions {
            paired_end: !sample.starts_with("RAP1_UNINDUCED"),
            ..CountOptions::default()
        };
        let counts = count_bam(&annotation, &input, &options).unwrap();
        let single = count_bam(
            &annotation,
            &input,
            &CountOptions {
                threads: 1,
                ..options
            },
        )
        .unwrap();
        assert_eq!(single, counts, "{sample}: thread-dependent counts");
        let mut command: Vec<String> = [
            "featureCounts",
            "-B",
            "-C",
            "-g",
            "gene_biotype",
            "-t",
            "exon",
        ]
        .map(str::to_owned)
        .to_vec();
        if options.paired_end {
            command.push("-p".into());
        }
        command.extend(
            [
                "-T",
                "6",
                "-a",
                "genome_gfp.gtf",
                "-s",
                "2",
                "-o",
                &output,
                &bam,
            ]
            .map(str::to_owned),
        );
        let inputs = [bam];
        let counts = [counts];
        let gold = fixture.join("results/star_salmon/featurecounts");
        let mut table = Vec::new();
        write_table(&mut table, &annotation, &inputs, &counts, &command).unwrap();
        identical(
            sample,
            "table",
            &table,
            &std::fs::read(gold.join(&output)).unwrap(),
        );
        let mut summary = Vec::new();
        write_summary(&mut summary, &inputs, &counts).unwrap();
        identical(
            sample,
            "summary",
            &summary,
            &std::fs::read(gold.join(format!("{output}.summary"))).unwrap(),
        );
        write_multiqc(&annotation, &counts[0], sample, &temp.join(sample)).unwrap();
        for suffix in ["biotype_counts_mqc.tsv", "biotype_counts_rrna_mqc.tsv"] {
            let name = format!("{sample}.{suffix}");
            identical(
                sample,
                suffix,
                &std::fs::read(temp.join(&name)).unwrap(),
                &std::fs::read(gold.join(name)).unwrap(),
            );
        }
        eprintln!("{sample}: table + summary + 2 MultiQC tables IDENTICAL; 1/6 threads identical");
    }
    std::fs::remove_dir_all(temp).unwrap();
}
