//! Pre-annotation evidence and an explicitly limited annotation-field gate.
//! This is not a full VCF/header/archive/XA or full post-annotation gate.
use nomics_somatic::{mutect2_calls, mutect2_regions::Inputs};
use std::{
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

fn identical(a: &Path, b: &Path) {
    let (a, b) = (fs::read(a).unwrap(), fs::read(b).unwrap());
    if a != b {
        let at = a
            .iter()
            .zip(&b)
            .position(|(x, y)| x != y)
            .unwrap_or(a.len().min(b.len()));
        let line = a[..at].iter().filter(|&&c| c == b'\n').count() + 1;
        panic!(
            "first differing byte {at}, line {line}; lengths {}/{}",
            a.len(),
            b.len()
        );
    }
}

fn annotation_fields(golden: &Path) -> (String, String) {
    let mut text = Vec::new();
    nomics_io::bgzf::decode_parallel(&fs::read(golden).unwrap(), 1, &mut |b| {
        text.extend_from_slice(b);
        Ok(())
    })
    .unwrap();
    let text = String::from_utf8(text).unwrap();
    let samples: Vec<_> = text
        .lines()
        .find(|s| s.starts_with("#CHROM"))
        .unwrap()
        .split('\t')
        .skip(9)
        .collect();
    let mut out = String::new();
    let mut phase = String::new();
    for line in text.lines().filter(|s| !s.starts_with('#')) {
        let f: Vec<_> = line.split('\t').collect();
        out += &format!("V {}:{} {},{}\n", f[0], f[1], f[3], f[4]);
        phase += &format!("V {}:{} {},{}\n", f[0], f[1], f[3], f[4]);
        let info: std::collections::BTreeMap<_, _> =
            f[7].split(';').filter_map(|s| s.split_once('=')).collect();
        for k in ["AS_SB_TABLE", "DP", "MBQ", "MFRL", "MMQ", "MPOS"] {
            out += &format!("I {k} {}\n", info[k]);
        }
        let keys: Vec<_> = f[8].split(':').collect();
        for (name, data) in samples.iter().zip(&f[9..]) {
            let fields: std::collections::BTreeMap<_, _> =
                keys.iter().copied().zip(data.split(':')).collect();
            out += &format!("S {name}\nAD {}\nDP {}\n", fields["AD"], fields["DP"]);
            for k in ["AF", "F1R2", "F2R1", "FAD", "SB"] {
                out += &format!("G {k} {}\n", fields[k]);
            }
            phase += &format!("S {name}\nGT {}\n", fields["GT"]);
            for k in ["PGT", "PID", "PS"] {
                if let Some(v) = fields.get(k) {
                    phase += &format!("G {k} {v}\n");
                }
            }
        }
    }
    (out, phase)
}

fn raw_annotation_subset(path: &Path) -> String {
    let mut out = String::new();
    for line in fs::read_to_string(path).unwrap().lines() {
        let f: Vec<_> = line.split_whitespace().collect();
        if f[0] == "S" {
            out += &format!("S {}\n", f[1]);
        } else if f[0] == "V"
            || f[0] == "AD"
            || f[0] == "DP"
            || (f[0] == "I" && ["AS_SB_TABLE", "DP", "MBQ", "MFRL", "MMQ", "MPOS"].contains(&f[1]))
            || (f[0] == "G" && ["AF", "F1R2", "F2R1", "FAD", "SB"].contains(&f[1]))
        {
            out += line;
            out.push('\n');
        }
    }
    out
}
#[test]
fn mutect2_pre_annotation_available_goldens() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap());
    let root = std::env::var_os("NOMICS_M2_CALLS_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic"));
    let reference = std::env::var_os("NOMICS_M2_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    let oracle = std::env::var_os("NOMICS_M2_CALLS_ORACLE").map(PathBuf::from);
    let output = std::env::temp_dir().join(format!("nomics-m2-calls-test-{}", std::process::id()));
    fs::create_dir_all(&output).unwrap();
    for (c, expected) in [
        ("chr21", 38),
        ("chr22", 134),
        ("chr19", 335),
        ("chr20", 144),
        ("chr13", 97),
        ("chr16", 449),
        ("chr12", 337),
    ] {
        let d = root.join(format!("m2-{c}"));
        if !d.join("oracle/haplotypes.bam").exists() {
            eprintln!(
                "skipping Mutect2 pre-annotation {c}: {} absent; override NOMICS_M2_CALLS_GOLDENS",
                d.display()
            );
            continue;
        }
        let time = Instant::now();
        let stats = output.join(format!("{c}.stats"));
        let calls = output.join(format!("{c}.calls.tsv"));
        let (_, s) = mutect2_calls::run(
            &Inputs {
                input: &[
                    d.join(format!("HCC1395BL_WES.{c}.bam")),
                    d.join(format!("HCC1395_WES.{c}.bam")),
                ],
                normal_sample: "HCC1395_HCC1395BL_WES",
                reference: &reference,
                intervals: &d.join(format!("{c}.bed")),
                panel_of_normals: &d.join(format!("1000g_pon.{c}.vcf.gz")),
                germline_resource: &d.join(format!("af-only-gnomad.{c}.vcf.gz")),
                regions: &output.join(format!("{c}.igv")),
                stats: &stats,
                dump_prefix: None,
            },
            &calls,
            &output.join("likelihoods.tsv"),
            &output.join("reads.tsv"),
        )
        .unwrap();
        let golden = fs::read_dir(d.join("oracle"))
            .unwrap()
            .map(|e| e.unwrap().path())
            .find(|p| p.to_string_lossy().ends_with(".vcf.gz.stats"))
            .unwrap();
        identical(&stats, &golden);
        let vcf = golden
            .to_string_lossy()
            .strip_suffix(".stats")
            .unwrap()
            .to_owned();
        let (expected_fields, expected_phase) = annotation_fields(Path::new(&vcf));
        let expected_path = output.join(format!("{c}.expected-annotation-fields.tsv"));
        fs::write(&expected_path, expected_fields).unwrap();
        identical(
            &PathBuf::from(format!("{}.annotation-fields.tsv", calls.display())),
            &expected_path,
        );
        eprintln!(
            "PASS {c}: six INFO annotation fields and seven genotype fields on {expected} calls (subset, no GT/phasing)"
        );
        let expected_path = output.join(format!("{c}.expected-phasing.tsv"));
        fs::write(&expected_path, expected_phase).unwrap();
        identical(
            &PathBuf::from(format!("{}.phasing.tsv", calls.display())),
            &expected_path,
        );
        eprintln!("PASS {c}: GT/PGT/PID/PS on {expected} calls (all samples)");
        assert_eq!(
            s.passing_events, expected,
            "{c} passing count only; not record parity"
        );
        if let Some(oracle) = oracle
            .as_ref()
            .map(|p| p.join(format!("{c}.calls.tsv")))
            .filter(|p| p.is_file())
        {
            identical(&calls, &oracle);
            eprintln!("PASS {c}: internal pre-annotation dump bytes");
            let annotations = PathBuf::from(format!("{}.annotations.tsv", oracle.display()));
            if annotations.is_file() {
                let expected = output.join(format!("{c}.expected-raw-annotations.tsv"));
                fs::write(&expected, raw_annotation_subset(&annotations)).unwrap();
                identical(
                    &PathBuf::from(format!("{}.annotations.tsv", calls.display())),
                    &expected,
                );
                eprintln!(
                    "PASS {c}: raw annotation subset bytes, AF binary64 bits (not the full internal gate)"
                );
            } else {
                eprintln!("skipping {c}: raw annotation oracle absent");
            }
        } else {
            eprintln!(
                "skipping {c} internal dump: supply qualified dump directory in NOMICS_M2_CALLS_ORACLE"
            );
        }
        eprintln!(
            "{c}: {s:?}, stats {} bytes, {:.3}s (passing count, not VCF parity)",
            fs::metadata(&stats).unwrap().len(),
            time.elapsed().as_secs_f64()
        );
        // Each contig has already passed; avoid retaining gigabytes of matrices.
        for entry in fs::read_dir(&output).unwrap() {
            let path = entry.unwrap().path();
            if path
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with(&format!("{c}."))
            {
                fs::remove_file(path).unwrap();
            }
        }
    }
    fs::remove_dir_all(output).unwrap();
}
