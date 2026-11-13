//! Ordered member names and exact bytes: outer tar/gzip metadata is not gated.
use std::path::{Path, PathBuf};
use std::process::Command;
fn only(dir: &Path, suffix: &str) -> PathBuf {
    let paths: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.to_string_lossy().ends_with(suffix))
        .collect();
    assert_eq!(paths.len(), 1);
    paths[0].clone()
}
fn untar(path: &Path, dir: &Path) -> Vec<u8> {
    std::fs::create_dir_all(dir).unwrap();
    assert!(
        Command::new("tar")
            .arg("-xzf")
            .arg(path)
            .arg("-C")
            .arg(dir)
            .status()
            .unwrap()
            .success()
    );
    let output = Command::new("tar").arg("-tzf").arg(path).output().unwrap();
    assert!(output.status.success());
    output.stdout
}
#[test]
fn orientation_golden_and_reversed_inputs() {
    let root = std::env::var_os("NOMICS_SOMATIC_GOLDEN")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap()).join("omics/runs/somatic/S1-r1/golden")
        });
    if !root
        .join("LEARNREADORIENTATIONMODEL-HCC1395_WES_vs_HCC1395BL_WES")
        .is_dir()
    {
        eprintln!("SKIP orientation golden: {} is absent", root.display());
        return;
    }
    let inputs = [
        only(
            &root.join("MUTECT2_PAIRED-HCC1395_WES_vs_HCC1395BL_WES"),
            ".f1r2.tar.gz",
        ),
        only(
            &root.join("MUTECT2_PAIRED-HCC1395_WES_vs_HCC1395BL_WES.2"),
            ".f1r2.tar.gz",
        ),
    ];
    let golden = only(
        &root.join("LEARNREADORIENTATIONMODEL-HCC1395_WES_vs_HCC1395BL_WES"),
        ".artifactprior.tar.gz",
    );
    let dir = std::env::temp_dir().join(format!(
        "nomics-orientation-test-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("ordered.tar.gz");
    nomics_somatic::orientation::run(&inputs, &output).unwrap();
    let want = untar(&golden, &dir.join("golden"));
    let got = untar(&output, &dir.join("ordered"));
    assert_eq!(got, want);
    let member = "HCC1395_HCC1395_WES.orientation_priors";
    let expected = std::fs::read(dir.join("golden").join(member)).unwrap();
    let actual = std::fs::read(dir.join("ordered").join(member)).unwrap();
    assert_eq!(
        actual.len(),
        expected.len(),
        "payload length; kept {}",
        dir.display()
    );
    assert!(
        actual == expected,
        "payload differs; kept {}",
        dir.display()
    );
    let output = dir.join("reversed.tar.gz");
    nomics_somatic::orientation::run(&[inputs[1].clone(), inputs[0].clone()], &output).unwrap();
    assert_eq!(untar(&output, &dir.join("reversed")), want);
    let reversed = std::fs::read(dir.join("reversed").join(member)).unwrap();
    assert_ne!(
        actual, reversed,
        "ordered floating-point accumulation must remain observable"
    );
    // Counts, contexts and row order survive reversal; only floating low bits differ.
    for (a, b) in std::str::from_utf8(&actual)
        .unwrap()
        .lines()
        .zip(std::str::from_utf8(&reversed).unwrap().lines())
        .skip(2)
    {
        let a: Vec<_> = a.split('\t').collect();
        let b: Vec<_> = b.split('\t').collect();
        assert_eq!(&a[..2], &b[..2]);
        assert_eq!(&a[14..], &b[14..]);
    }
    std::fs::remove_dir_all(dir).unwrap();
}
