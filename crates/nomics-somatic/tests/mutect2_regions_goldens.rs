//! Optional local chr21 tumour-normal gate. No reference executable is called.
use nomics_somatic::mutect2_regions::{self, Inputs};
use std::{
    fs,
    io::{BufReader, Read},
    path::{Path, PathBuf},
};

fn compare(a: &Path, b: &Path) {
    let (mut a, mut b) = (
        BufReader::new(fs::File::open(a).unwrap()),
        BufReader::new(fs::File::open(b).unwrap()),
    );
    let (mut x, mut y) = ([0; 65536], [0; 65536]);
    let mut offset = 0;
    loop {
        let n = a.read(&mut x).unwrap();
        let m = b.read(&mut y).unwrap();
        assert_eq!(n, m, "length differs at byte {offset}");
        assert_eq!(x[..n], y[..m], "content differs in block at byte {offset}");
        if n == 0 {
            break;
        }
        offset += n;
    }
}
#[test]
fn mutect2_chr21_regions_stats_and_filters() {
    let root = std::env::var_os("NOMICS_M2_GOLDEN")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap()).join("omics/runs/somatic/m2-chr21")
        });
    if !root.is_dir() {
        eprintln!(
            "skipping Mutect2 regions golden: {} absent (override NOMICS_M2_GOLDEN)",
            root.display()
        );
        return;
    }
    let reference = std::env::var_os("NOMICS_M2_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap()).join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    let output = std::env::temp_dir().join(format!("nomics-m2-regions-{}", std::process::id()));
    fs::create_dir_all(&output).unwrap();
    let oracle = std::env::var_os("NOMICS_M2_ORACLE_PREFIX").map(PathBuf::from);
    let dump = oracle.as_ref().map(|_| output.join("dump"));
    let stats = mutect2_regions::run(&Inputs {
        input: &[
            root.join("HCC1395BL_WES.chr21.bam"),
            root.join("HCC1395_WES.chr21.bam"),
        ],
        normal_sample: "HCC1395_HCC1395BL_WES",
        reference: &reference,
        intervals: &root.join("chr21.bed"),
        panel_of_normals: &root.join("1000g_pon.chr21.vcf.gz"),
        germline_resource: &root.join("af-only-gnomad.chr21.vcf.gz"),
        regions: &output.join("regions.igv"),
        stats: &output.join("X.stats"),
        dump_prefix: dump.as_deref(),
    })
    .unwrap();
    compare(
        &output.join("regions.igv"),
        &root.join("oracle/regions.igv"),
    );
    let golden_stats = fs::read_dir(root.join("r1"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .find(|p| p.to_string_lossy().ends_with(".vcf.gz.stats"))
        .unwrap();
    compare(&output.join("X.stats"), &golden_stats);
    assert_eq!(
        (stats.regions, stats.active, stats.callable),
        (7167, 699, 1230643)
    );
    let log = fs::read_to_string(root.join("r1/gatk.log")).unwrap();
    let mut summary = String::new();
    for line in log.lines() {
        if let Some(end) = line.find(" read(s) filtered by:") {
            let start = line[..end].rfind(' ').map_or(0, |i| i + 1);
            summary.push_str(&line[start..]);
            summary.push('\n');
        } else if line.contains(" total reads filtered out of ") {
            summary.push_str(line);
            summary.push('\n');
        }
    }
    assert_eq!(stats.filter_summary(), summary);
    if let Some(oracle) = oracle {
        for suffix in ["activity.tsv", "reads.tsv", "prepared.bin"] {
            compare(
                &PathBuf::from(format!("{}.{}", dump.as_ref().unwrap().display(), suffix)),
                &PathBuf::from(format!("{}.{}", oracle.display(), suffix)),
            );
        }
    }
    eprintln!(
        "Mutect2: 7167/7167 regions; 699 active; stats 35/35 bytes; 12 filters and totals exact"
    );
    fs::remove_dir_all(output).unwrap();
}
