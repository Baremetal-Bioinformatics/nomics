//! FilterMutectCalls S1-r1 gate. Only GATKCommandLine is normalized.
use nomics_somatic::filtering::{self, Inputs};
use std::{
    fs,
    io::Read,
    path::{Path, PathBuf},
};
fn decompress(path: &Path) -> String {
    let mut r = nomics_io::gzip::GzReader::open(path).unwrap();
    let mut s = String::new();
    r.read_to_string(&mut s).unwrap();
    s
}
#[test]
fn filtering_goldens() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap_or_default());
    let golden = std::env::var_os("NOMICS_SOMATIC_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic/S1-r1/golden"));
    if !golden.is_dir() {
        eprintln!("SKIP FilterMutectCalls: {} absent", golden.display());
        return;
    }
    let reference = std::env::var_os("NOMICS_SOMATIC_FASTA")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    let x = "HCC1395_WES_vs_HCC1395BL_WES";
    let input = |dir: &str, suffix: &str| golden.join(format!("{dir}-{x}/{x}.mutect2{suffix}"));
    let output = std::env::temp_dir().join(format!(
        "nomics-filter-golden-{}.vcf.gz",
        std::process::id()
    ));
    filtering::run(&Inputs {
        variant: &input("MERGE_MUTECT2", ".vcf.gz"),
        stats: Some(&input("MERGEMUTECTSTATS", ".vcf.gz.stats")),
        reference: &reference,
        orientation_priors: &input("LEARNREADORIENTATIONMODEL", ".artifactprior.tar.gz"),
        segmentation: &input("CALCULATECONTAMINATION", ".segmentation.table"),
        contamination: &input("CALCULATECONTAMINATION", ".contamination.table"),
        output: &output,
    })
    .unwrap();
    let actual = decompress(&output);
    let expected = decompress(&input("FILTERMUTECTCALLS", ".filtered.vcf.gz"));
    let records = |s: &str| {
        s.split_inclusive('\n')
            .filter(|l| !l.starts_with('#'))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let (a, e) = (records(&actual), records(&expected));
    let identical = a.iter().zip(&e).filter(|(a, e)| a == e).count();
    let pass = a
        .iter()
        .filter(|l| l.split('\t').nth(6) == Some("PASS"))
        .count();
    eprintln!(
        "FilterMutectCalls: {identical}/{} records identical; {pass} PASS; output {}",
        e.len(),
        output.display()
    );
    let stats_path = PathBuf::from(format!("{}.filteringStats.tsv", output.display()));
    let a_stats = fs::read(&stats_path).unwrap();
    let e_stats = fs::read(input(
        "FILTERMUTECTCALLS",
        ".filtered.vcf.gz.filteringStats.tsv",
    ))
    .unwrap();
    eprintln!(
        "Stats: {}/{} bytes; identical={}",
        a_stats.len(),
        e_stats.len(),
        a_stats == e_stats
    );
    if let Some(i) = a.iter().zip(&e).position(|(a, e)| a != e) {
        eprintln!(
            "First mismatch record {}\nREF {}\nOUT {}",
            i + 1,
            e[i],
            a[i]
        );
    }
    assert_eq!(a_stats, e_stats, "filteringStats bytes");
    assert_eq!(e.len(), 6351);
    assert_eq!(identical, 6351);
    assert_eq!(pass, 2024);
    let header = |s: &str| {
        s.split_inclusive('\n')
            .filter(|l| l.starts_with('#') && !l.starts_with("##GATKCommandLine="))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        header(&actual),
        header(&expected),
        "headers except GATKCommandLine"
    );
    fs::remove_file(output).unwrap();
    fs::remove_file(stats_path).unwrap();
}

#[test]
fn filtering_perturbed_inputs() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap_or_default());
    let root = std::env::var_os("NOMICS_SOMATIC_FMC_PERTURB")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic/fmc-perturb"));
    if !root.is_dir() {
        eprintln!(
            "SKIP FilterMutectCalls perturbations: {} absent",
            root.display()
        );
        return;
    }
    let reference = std::env::var_os("NOMICS_SOMATIC_FASTA")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    for (case, variant, contamination, count) in [
        ("odd", "odd", "contam.table", 3176),
        ("chr1to7", "chr1to7", "contam.table", 2666),
        ("rest", "rest", "contam.table", 3685),
        ("contam10", "raw", "contam10.table", 6351),
        ("same", "raw", "contam.table", 6351),
    ] {
        let output = std::env::temp_dir().join(format!(
            "nomics-filter-perturb-{}-{case}.vcf.gz",
            std::process::id()
        ));
        filtering::run(&Inputs {
            variant: &root.join(format!("{variant}.vcf.gz")),
            stats: Some(&root.join(format!("{variant}.vcf.gz.stats"))),
            reference: &reference,
            orientation_priors: &root.join("priors.tar.gz"),
            segmentation: &root.join("seg.table"),
            contamination: &root.join(contamination),
            output: &output,
        })
        .unwrap();
        let stats_path = PathBuf::from(format!("{}.filteringStats.tsv", output.display()));
        let a_stats = fs::read(&stats_path).unwrap();
        let e_stats =
            fs::read(root.join(format!("{case}.filtered.vcf.gz.filteringStats.tsv"))).unwrap();
        // Keep the failure readable: exact first differing line, then all bytes.
        if a_stats != e_stats {
            for (i, (a, e)) in a_stats
                .split(|&b| b == b'\n')
                .zip(e_stats.split(|&b| b == b'\n'))
                .enumerate()
            {
                if a != e {
                    panic!(
                        "{case}: stats line {}\nREF {}\nOUT {}",
                        i + 1,
                        String::from_utf8_lossy(e),
                        String::from_utf8_lossy(a)
                    );
                }
            }
        }
        assert_eq!(a_stats, e_stats, "{case}: filteringStats bytes");
        let actual = decompress(&output);
        let expected = decompress(&root.join(format!("{case}.filtered.vcf.gz")));
        let records = |s: &str| {
            s.split_inclusive('\n')
                .filter(|l| !l.starts_with('#'))
                .map(str::to_owned)
                .collect::<Vec<_>>()
        };
        let (a, e) = (records(&actual), records(&expected));
        assert_eq!(e.len(), count, "{case}: fixture record count");
        assert_eq!(a.len(), e.len(), "{case}: record count");
        for (i, (a, e)) in a.iter().zip(&e).enumerate() {
            assert_eq!(a, e, "{case}: complete record {}", i + 1);
        }
        let header = |s: &str| {
            s.split_inclusive('\n')
                .filter(|l| l.starts_with('#') && !l.starts_with("##GATKCommandLine="))
                .map(str::to_owned)
                .collect::<Vec<_>>()
        };
        assert_eq!(
            header(&actual),
            header(&expected),
            "{case}: normalized header"
        );
        eprintln!("{case}: records {count}/{count} identical; stats identical; header identical");
        fs::remove_file(output).unwrap();
        fs::remove_file(stats_path).unwrap();
    }
}
