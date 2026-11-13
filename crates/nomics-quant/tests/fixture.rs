//! The external nf-core fixture is optional; golden bytes are never read by
//! nomics-quant itself. NOMICS_RNASEQ_FIXTURE overrides the default location.
use std::path::{Path, PathBuf};
use std::time::Instant;

use nomics_quant::{Options, run};

const SAMPLES: [(&str, &str); 5] = [
    ("WT_REP1", "ISR"),
    ("WT_REP2", "ISR"),
    ("RAP1_UNINDUCED_REP1", "SR"),
    ("RAP1_UNINDUCED_REP2", "SR"),
    ("RAP1_IAA_30M_REP1", "ISR"),
];

fn fixture() -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("NOMICS_RNASEQ_FIXTURE") {
        let path = PathBuf::from(path);
        assert!(
            path.is_dir(),
            "explicit NOMICS_RNASEQ_FIXTURE is missing: {}",
            path.display()
        );
        return Some(path);
    }
    let path = PathBuf::from(std::env::var_os("HOME")?).join("omics-fixtures/rnaseq-test");
    path.is_dir().then_some(path)
}

fn first_difference(ours: &[u8], golden: &[u8]) -> String {
    let a = String::from_utf8_lossy(ours);
    let b = String::from_utf8_lossy(golden);
    let a: Vec<_> = a.lines().collect();
    let b: Vec<_> = b.lines().collect();
    for i in 0..a.len().max(b.len()) {
        if a.get(i) != b.get(i) {
            return format!(
                "line {}\nours:   {}\ngolden: {}",
                i + 1,
                a.get(i).unwrap_or(&"<EOF>"),
                b.get(i).unwrap_or(&"<EOF>")
            );
        }
    }
    "bytes differ (line endings or final newline)".to_owned()
}

fn uncompressed(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    nomics_io::gzip::gunzip_to(&std::fs::read(path).unwrap(), &mut |bytes| {
        out.extend_from_slice(bytes);
        Ok(())
    })
    .unwrap();
    out
}

#[test]
fn nf_core_five_samples_byte_gate() {
    let Some(fixture) = fixture() else {
        eprintln!("SKIP: optional rnaseq-test fixture absent (set NOMICS_RNASEQ_FIXTURE)");
        return;
    };
    let output =
        std::env::temp_dir().join(format!("nomics-quant-fixture-test-{}", std::process::id()));
    let mut failures = Vec::new();
    for (sample, lib_type) in SAMPLES {
        let golden = fixture.join("results/star_salmon").join(sample);
        let cmd: serde_json::Value =
            serde_json::from_slice(&std::fs::read(golden.join("cmd_info.json")).unwrap()).unwrap();
        assert_eq!(
            cmd["libType"], lib_type,
            "fixture command changed for {sample}"
        );
        let start = Instant::now();
        let q = run(&Options {
            transcripts: fixture.join("transcriptome.fasta"),
            alignments: fixture
                .join("results/star_salmon")
                .join(format!("{sample}.Aligned.toTranscriptome.out.bam")),
            gene_map: Some(fixture.join("genome_gfp.gtf")),
            output: output.join(sample),
            lib_type: lib_type.to_owned(),
            incompat_prior: 0.0,
            threads: 6,
        })
        .unwrap();
        eprintln!(
            "{sample}: {:.3}s, {} mapped, {} classes, {} iterations",
            start.elapsed().as_secs_f64(),
            q.num_mapped,
            q.num_eq_classes,
            q.iterations
        );
        for file in [
            "quant.sf",
            "quant.genes.sf",
            "lib_format_counts.json",
            "aux_info/ambig_info.tsv",
            "aux_info/expected_bias.gz",
            "aux_info/observed_bias.gz",
            "aux_info/observed_bias_3p.gz",
        ] {
            let ours = std::fs::read(output.join(sample).join(file)).unwrap();
            let reference = std::fs::read(golden.join(file)).unwrap();
            if ours != reference {
                failures.push(format!(
                    "{sample}/{file}: {}",
                    first_difference(&ours, &reference)
                ));
            }
        }
        // Documented GATES.md residual: zlib-fast level 3 vs zlib-rs level 6.
        let ours = uncompressed(&output.join(sample).join("aux_info/fld.gz"));
        let reference = uncompressed(&golden.join("aux_info/fld.gz"));
        if ours != reference {
            failures.push(format!("{sample}/fld.gz: decompressed histogram differs"));
        }
    }
    assert!(
        failures.is_empty(),
        "{}\noutputs retained at {}",
        failures.join("\n"),
        output.display()
    );
    std::fs::remove_dir_all(output).unwrap();
}

#[test]
fn human_iu_byte_gate() {
    let fixture = if let Some(path) = std::env::var_os("NOMICS_SALMON_IU_FIXTURE") {
        let path = PathBuf::from(path);
        assert!(
            path.is_dir(),
            "explicit IU fixture is missing: {}",
            path.display()
        );
        path
    } else if let Some(home) = std::env::var_os("HOME") {
        let path = PathBuf::from(home).join("omics-fixtures/salmon-iu");
        if !path.is_dir() {
            eprintln!("SKIP: optional salmon-iu fixture absent");
            return;
        }
        path
    } else {
        eprintln!("SKIP: no default IU fixture location");
        return;
    };
    let golden = fixture.join("golden-r1");
    let output =
        std::env::temp_dir().join(format!("nomics-quant-human-iu-test-{}", std::process::id()));
    let start = Instant::now();
    let q = run(&Options {
        transcripts: fixture.join("genome.transcripts.fa"),
        alignments: fixture.join("sub.bam"),
        gene_map: Some(fixture.join("annotation.gtf")),
        output: output.clone(),
        lib_type: "IU".to_owned(),
        incompat_prior: 0.0,
        threads: 6,
    })
    .unwrap();
    eprintln!(
        "human IU: {:.3}s, {} mapped, {} classes, {} iterations",
        start.elapsed().as_secs_f64(),
        q.num_mapped,
        q.num_eq_classes,
        q.iterations
    );
    let mut failures = Vec::new();
    for file in [
        "quant.sf",
        "quant.genes.sf",
        "lib_format_counts.json",
        "aux_info/ambig_info.tsv",
        "aux_info/expected_bias.gz",
        "aux_info/observed_bias.gz",
        "aux_info/observed_bias_3p.gz",
    ] {
        let ours = std::fs::read(output.join(file)).unwrap();
        let reference = std::fs::read(golden.join(file)).unwrap();
        if ours != reference {
            failures.push(format!(
                "IU/{file}: {}",
                first_difference(&ours, &reference)
            ));
        }
    }
    if uncompressed(&output.join("aux_info/fld.gz"))
        != uncompressed(&golden.join("aux_info/fld.gz"))
    {
        failures.push("IU/fld.gz: decompressed histogram differs".to_owned());
    }
    assert!(
        failures.is_empty(),
        "{}\noutputs retained at {}",
        failures.join("\n"),
        output.display()
    );
    std::fs::remove_dir_all(output).unwrap();
}
