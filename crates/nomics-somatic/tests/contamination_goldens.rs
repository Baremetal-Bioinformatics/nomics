//! Local S1-r1 gate. Override NOMICS_SOMATIC_GOLDENS to relocate the read-only
//! fixture, or NOMICS_SOMATIC_DICTIONARY for the GRCh38 .dict.
use nomics_somatic::{contamination, pileup_summary};
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};
fn temp(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "nomics-{name}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    p
}
fn cmp(a: &Path, b: &Path) {
    assert!(
        Command::new("cmp")
            .arg(a)
            .arg(b)
            .status()
            .unwrap()
            .success(),
        "cmp {} {}",
        a.display(),
        b.display()
    );
}
#[test]
fn contamination_goldens() {
    let home = std::env::var_os("HOME")
        .map(PathBuf::from)
        .unwrap_or_default();
    let golden = std::env::var_os("NOMICS_SOMATIC_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic/S1-r1/golden"));
    if !golden.is_dir() {
        eprintln!("SKIP contamination goldens: {} absent", golden.display());
        return;
    }
    let dictionary = std::env::var_os("NOMICS_SOMATIC_DICTIONARY")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict",
            )
        });
    assert!(
        dictionary.is_file(),
        "goldens present but dictionary missing: {}",
        dictionary.display()
    );
    let out = temp("contamination-goldens");
    let mut gathered = Vec::new();
    // Recorded .command.sh argument order: chr17 first, chr1 second.
    for (role, sample) in [("TUMOR", "HCC1395_WES"), ("NORMAL", "HCC1395BL_WES")] {
        let inputs=[golden.join(format!("GETPILEUPSUMMARIES_{role}-{sample}/{sample}.mutect2.chr17_38877931-38878296.pileups.table")),golden.join(format!("GETPILEUPSUMMARIES_{role}-{sample}.2/{sample}.mutect2.chr1_11981-12351.pileups.table"))];
        let output = out.join(format!("{sample}.pileups.table"));
        pileup_summary::gather(&inputs, &dictionary, &output).unwrap();
        cmp(
            &output,
            &golden.join(format!(
                "GATHERPILEUPSUMMARIES_{role}-{sample}/{sample}.mutect2.pileups.table"
            )),
        );
        gathered.push(output);
    }
    let c = out.join("contamination.table");
    let s = out.join("segmentation.table");
    contamination::calculate(&gathered[0], Some(&gathered[1]), &c, Some(&s), 0.5, 3.).unwrap();
    for (output, kind) in [(&c, "contamination"), (&s, "segmentation")] {
        cmp(output,&golden.join(format!("CALCULATECONTAMINATION-HCC1395_WES_vs_HCC1395BL_WES/HCC1395_WES_vs_HCC1395BL_WES.mutect2.{kind}.table")));
    }
    fs::remove_dir_all(out).unwrap();
}
#[test]
fn gather_semantics_and_tumor_only_smoke() {
    let out = temp("somatic-synthetic");
    let dictionary = out.join("ref.dict");
    fs::write(
        &dictionary,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n",
    )
    .unwrap();
    let row = |contig: &str, pos: i32, alt: i32| pileup_summary::PileupSummary {
        contig: contig.into(),
        position: pos,
        ref_count: 50 - alt,
        alt_count: alt,
        other_alt_count: 0,
        allele_frequency: 0.1,
    };
    let a = out.join("a");
    let b = out.join("b");
    let empty = out.join("empty");
    let gathered = out.join("gathered");
    pileup_summary::write(
        &a,
        Some("sample"),
        &[row("chr2", 10, 25), row("chr2", 20, 0)],
    )
    .unwrap();
    pileup_summary::write(
        &b,
        Some("sample"),
        &[row("chr1", 10, 25), row("chr1", 20, 0), row("chr1", 30, 25)],
    )
    .unwrap();
    pileup_summary::write(&empty, Some("different-but-empty"), &[]).unwrap();
    pileup_summary::gather(
        &[a.clone(), empty.clone(), b.clone()],
        &dictionary,
        &gathered,
    )
    .unwrap();
    assert_eq!(pileup_summary::read(&gathered).unwrap().1[0].contig, "chr1");
    let c = out.join("c");
    let seg = out.join("s");
    contamination::calculate(&gathered, None, &c, Some(&seg), 0.5, 3.).unwrap();
    assert_eq!(fs::read_to_string(&c).unwrap().lines().count(), 2);
    assert!(fs::read_to_string(&seg).unwrap().contains("chr1\t10\t30"));
    // All-empty gather closes the writer with just the columns, no metadata.
    pileup_summary::gather(std::slice::from_ref(&empty), &dictionary, &gathered).unwrap();
    assert_eq!(
        fs::read_to_string(&gathered).unwrap(),
        format!("{}\n", pileup_summary::HEADER)
    );
    contamination::calculate(&empty, None, &c, Some(&seg), 0.5, 3.).unwrap();
    assert!(fs::read_to_string(&c).unwrap().ends_with("\t0.0\t1.0\n"));
    // Equal first loci preserve file order, and overlapping/duplicate rows survive.
    pileup_summary::write(&a, Some("sample"), &[row("chr1", 10, 1)]).unwrap();
    pileup_summary::gather(&[a.clone(), b.clone()], &dictionary, &gathered).unwrap();
    assert_eq!(pileup_summary::read(&gathered).unwrap().1[0].alt_count, 1);
    pileup_summary::write(&a, Some("other-sample"), &[row("chr1", 10, 1)]).unwrap();
    assert!(pileup_summary::gather(&[a, b], &dictionary, &gathered).is_err());
    fs::remove_dir_all(out).unwrap();
}
