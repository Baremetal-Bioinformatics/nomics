//! Optional collector boundary and downstream chr21 gate. No oracle executable.
use nomics_somatic::{
    mutect2_f1r2::{self, Inputs},
    orientation::{self, tar},
};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

fn members(path: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut result = BTreeMap::new();
    tar::read(path, |name, bytes| {
        assert!(result.insert(name.to_owned(), bytes.to_vec()).is_none());
        Ok(())
    })
    .unwrap();
    result
}
fn compare_alt(expected: &[u8], actual: &[u8]) {
    let a: Vec<_> = expected.split_inclusive(|&b| b == b'\n').collect();
    let b: Vec<_> = actual.split_inclusive(|&b| b == b'\n').collect();
    assert_eq!(a.len(), b.len());
    let h = a
        .iter()
        .position(|line| line.starts_with(b"depth\t"))
        .unwrap();
    let fields = |line: &[u8]| -> Vec<Vec<u8>> {
        assert!(line.ends_with(b"\n") && !line.ends_with(b"\r\n"));
        line[..line.len() - 1]
            .split(|&c| c == b'\t')
            .map(<[u8]>::to_vec)
            .collect()
    };
    let ha = fields(a[h]);
    let hb = fields(b[h]);
    assert_eq!(ha.len(), 385);
    assert_eq!(hb.len(), 385);
    let map: BTreeMap<_, _> = hb.iter().enumerate().map(|(i, c)| (c, i)).collect();
    assert_eq!(map.len(), 385);
    let order: Vec<_> = ha
        .iter()
        .map(|c| *map.get(c).expect("missing alt column"))
        .collect();
    assert_eq!(
        ha.iter().collect::<std::collections::BTreeSet<_>>().len(),
        385
    );
    let mut rows = 0;
    for i in 0..a.len() {
        if i == h {
            continue;
        }
        if i > h && a[i] != b"\n" {
            let (x, y) = (fields(a[i]), fields(b[i]));
            assert_eq!(x.len(), 385);
            assert_eq!(y.len(), 385);
            for (col, &j) in order.iter().enumerate() {
                assert_eq!(x[col], y[j], "line {}, column {:?}", i + 1, ha[col]);
            }
            rows += 1;
        } else {
            assert_eq!(a[i], b[i], "non-column line {}", i + 1);
        }
    }
    assert_eq!(rows, 200);
}
#[test]
fn chr21_collector_members_and_learned_priors() {
    let root = std::env::var_os("NOMICS_M2_GOLDEN")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap()).join("omics/runs/somatic/m2-chr21")
        });
    if !root.is_dir() {
        eprintln!(
            "skipping Mutect2 f1r2 golden: {} absent (override NOMICS_M2_GOLDEN)",
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
    let output = std::env::temp_dir().join(format!("nomics-f1r2-test-{}", std::process::id()));
    fs::create_dir_all(&output).unwrap();
    let archive = output.join("chr21.f1r2.tar.gz");
    mutect2_f1r2::run(&Inputs {
        input: &[
            root.join("HCC1395BL_WES.chr21.bam"),
            root.join("HCC1395_WES.chr21.bam"),
        ],
        normal_sample: "HCC1395_HCC1395BL_WES",
        reference: &reference,
        intervals: &root.join("chr21.bed"),
        panel_of_normals: &root.join("1000g_pon.chr21.vcf.gz"),
        germline_resource: &root.join("af-only-gnomad.chr21.vcf.gz"),
        output: &archive,
    })
    .unwrap();
    let golden = root.join("oracle/HCC1395_WES_vs_HCC1395BL_WES.mutect2.chr21.f1r2.tar.gz");
    let (a, b) = (members(&golden), members(&archive));
    assert_eq!(a.len(), 3);
    assert_eq!(a.keys().collect::<Vec<_>>(), b.keys().collect::<Vec<_>>());
    for (name, expected) in &a {
        let actual = &b[name];
        if name.ends_with(".alt_histogram") {
            compare_alt(expected, actual);
        } else {
            assert_eq!(expected, actual, "{name}");
        }
    }
    let mut order = Vec::new();
    tar::read(&archive, |name, _| {
        order.push(name.rsplit('.').next().unwrap().to_owned());
        Ok(())
    })
    .unwrap();
    assert_eq!(order, ["ref_histogram", "alt_histogram", "alt_table"]);
    let (port_prior, oracle_prior) = (
        output.join("port.priors.tar.gz"),
        output.join("oracle.priors.tar.gz"),
    );
    orientation::run(&[archive], &port_prior).unwrap();
    orientation::run(&[golden], &oracle_prior).unwrap();
    assert_eq!(members(&port_prior), members(&oracle_prior));
    // The deterministic transport also allows an exact whole-archive check.
    assert_eq!(
        fs::read(&port_prior).unwrap(),
        fs::read(&oracle_prior).unwrap()
    );
    fs::remove_dir_all(output).unwrap();
}
