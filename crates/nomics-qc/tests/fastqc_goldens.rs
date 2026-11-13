//! Byte gates against the two text members of the tier0 reference ZIPs.
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

// Central-directory sizes are authoritative (FastQC uses data descriptors).
// Stored and raw DEFLATE only; no compression dependency or external unzip.
fn zip_member(path: &Path, suffix: &str) -> Vec<u8> {
    let zip = fs::read(path).unwrap();
    let u16at = |p| u16::from_le_bytes(zip[p..p + 2].try_into().unwrap()) as usize;
    let u32at = |p| u32::from_le_bytes(zip[p..p + 4].try_into().unwrap()) as usize;
    let end = (0..zip.len() - 3)
        .rev()
        .find(|&i| zip[i..i + 4] == [80, 75, 5, 6])
        .unwrap();
    let mut at = u32at(end + 16);
    for _ in 0..u16at(end + 10) {
        assert_eq!(&zip[at..at + 4], b"PK\x01\x02");
        let name = std::str::from_utf8(&zip[at + 46..at + 46 + u16at(at + 28)]).unwrap();
        if name.ends_with(&format!("/{suffix}")) {
            let local = u32at(at + 42);
            assert_eq!(&zip[local..local + 4], b"PK\x03\x04");
            assert_eq!(u16at(at + 8) & 1, 0, "encrypted ZIP");
            let start = local + 30 + u16at(local + 26) + u16at(local + 28);
            let compressed = &zip[start..start + u32at(at + 20)];
            let mut data = Vec::new();
            match u16at(at + 10) {
                0 => data.extend_from_slice(compressed),
                8 => {
                    let used =
                        nomics_io::inflate::inflate(compressed, &mut data, u32at(at + 24)).unwrap();
                    assert_eq!(used, compressed.len());
                }
                method => panic!("unsupported ZIP method {method}"),
            }
            assert_eq!(data.len(), u32at(at + 24));
            assert_eq!(
                nomics_io::crc32::Crc32::new().update(&data).finish() as usize,
                u32at(at + 16)
            );
            return data;
        }
        at += 46 + u16at(at + 28) + u16at(at + 30) + u16at(at + 32);
    }
    panic!("missing {suffix} in {}", path.display());
}

#[test]
fn tier0_both_mates() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap());
    let runs = std::env::var_os("NOMICS_RUNS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs"));
    let inputs = std::env::var_os("NOMICS_FASTQC_INPUTS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/data/tier0"));
    let gold = runs.join("tier0-wes-20261008-r1/results/reports/fastqc/HG002_WES-L1");
    if !(1..=2).all(|m| {
        inputs.join(format!("tier0_{m}.fastq.gz")).exists()
            && gold.join(format!("HG002_WES-L1_{m}_fastqc.zip")).exists()
    }) {
        eprintln!(
            "skip FastQC goldens: both tier0 mates and ZIPs required (NOMICS_RUNS / NOMICS_FASTQC_INPUTS)"
        );
        return;
    }
    let tmp = std::env::temp_dir().join(format!("nomics-fastqc-gate-{}", std::process::id()));
    fs::create_dir_all(&tmp).unwrap();
    for m in 1..=2 {
        let start = Instant::now();
        let out = tmp.join(m.to_string());
        nomics_qc::fastqc::fastqc(
            &inputs.join(format!("tier0_{m}.fastq.gz")),
            &out,
            Some(&format!("HG002_WES-L1_{m}.gz")),
        )
        .unwrap();
        eprintln!("mate {m}: {:.3}s", start.elapsed().as_secs_f64());
        for member in ["fastqc_data.txt", "summary.txt"] {
            let expected = zip_member(&gold.join(format!("HG002_WES-L1_{m}_fastqc.zip")), member);
            let actual = fs::read(out.join(member)).unwrap();
            if actual != expected {
                let pos = actual
                    .iter()
                    .zip(&expected)
                    .position(|(a, b)| a != b)
                    .unwrap_or(actual.len().min(expected.len()));
                panic!(
                    "mate {m} {member}: first difference byte {pos}; files retained in {}",
                    out.display()
                );
            }
        }
    }
    fs::remove_dir_all(tmp).unwrap();
}
