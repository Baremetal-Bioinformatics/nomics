//! MergeVcfs and MergeMutectStats against the S1-r1 golden: the two Mutect2 scatter groups of
//! the HCC1395 pair, in the order sarek gave them to GATK, and in the other order.
use nomics_somatic::merge;
use std::{
    io::Read,
    path::{Path, PathBuf},
};

fn decompress(path: &Path) -> Vec<u8> {
    let mut r = nomics_io::gzip::GzReader::open(path).unwrap();
    let mut s = Vec::new();
    r.read_to_end(&mut s).unwrap();
    s
}

fn without_command_line(text: &[u8]) -> Vec<&[u8]> {
    text.split_inclusive(|&b| b == b'\n')
        .filter(|l| !l.starts_with(b"##GATKCommandLine="))
        .collect()
}

#[test]
fn merge_goldens() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap_or_default());
    let golden = std::env::var_os("NOMICS_SOMATIC_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic/S1-r1/golden"));
    if !golden.is_dir() {
        eprintln!("SKIP MergeVcfs, MergeMutectStats: {} absent", golden.display());
        return;
    }
    let x = "HCC1395_WES_vs_HCC1395BL_WES";
    // sarek's order: the chr17 group first.
    let groups = [
        (format!("MUTECT2_PAIRED-{x}"), "chr17_38877931-38878296"),
        (format!("MUTECT2_PAIRED-{x}.2"), "chr1_11981-12351"),
    ];
    let vcfs: Vec<PathBuf> =
        groups.iter().map(|(d, g)| golden.join(format!("{d}/{x}.mutect2.{g}.vcf.gz"))).collect();
    let stats: Vec<PathBuf> =
        groups.iter().map(|(d, g)| golden.join(format!("{d}/{x}.mutect2.{g}.vcf.gz.stats"))).collect();
    let dict = home.join("omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict");
    let dict = dict.is_file().then_some(dict);
    let out = std::env::temp_dir().join(format!("nomics-merge-golden-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();

    let expected = decompress(&golden.join(format!("MERGE_MUTECT2-{x}/{x}.mutect2.vcf.gz")));
    let merged = out.join("merged.vcf.gz");
    let n = merge::merge_vcfs(&vcfs, dict.as_deref(), &merged).unwrap();
    let actual = decompress(&merged);
    // In sarek's order the output is GATK's in every byte, the command line included.
    assert!(actual == expected, "merged VCF differs from GATK's");
    assert_eq!(n, 6351);

    // In the other order only the kept ##GATKCommandLine changes.
    let reversed: Vec<PathBuf> = vcfs.iter().rev().cloned().collect();
    merge::merge_vcfs(&reversed, dict.as_deref(), &merged).unwrap();
    let other = decompress(&merged);
    assert!(other != expected);
    assert!(without_command_line(&other) == without_command_line(&expected));

    let expected = std::fs::read(golden.join(format!("MERGEMUTECTSTATS-{x}/{x}.mutect2.vcf.gz.stats"))).unwrap();
    let merged = out.join("merged.stats");
    merge::merge_mutect_stats(&stats, &merged).unwrap();
    assert_eq!(std::fs::read(&merged).unwrap(), expected);
    std::fs::remove_dir_all(out).unwrap();
    eprintln!("MergeVcfs: 6,351 records and the header identical; MergeMutectStats: bytes identical");
}
