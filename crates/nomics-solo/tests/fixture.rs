//! Full PBMC1k gate, automatically active when the local oracle is ready.
//! `cargo test --release -p nomics-solo --test fixture -- --nocapture`.
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

const FILES: [&str; 10] = [
    "Barcodes.stats",
    "Gene/Features.stats",
    "Gene/Summary.csv",
    "Gene/UMIperCellSorted.txt",
    "Gene/raw/matrix.mtx",
    "Gene/raw/barcodes.tsv",
    "Gene/raw/features.tsv",
    "Gene/filtered/matrix.mtx",
    "Gene/filtered/barcodes.tsv",
    "Gene/filtered/features.tsv",
];

fn reader(path: &Path) -> Box<dyn BufRead> {
    if path.exists() {
        Box::new(BufReader::new(std::fs::File::open(path).unwrap()))
    } else {
        Box::new(BufReader::new(
            nomics_io::gzip::GzReader::open(path.with_extension(format!(
                "{}.gz",
                path.extension().unwrap().to_str().unwrap()
            )))
            .unwrap(),
        ))
    }
}

fn first_diff(expected: &Path, actual: &Path) -> Option<String> {
    let (mut a, mut b) = (reader(expected), reader(actual));
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for line in 1.. {
        x.clear();
        y.clear();
        a.read_until(b'\n', &mut x).unwrap();
        b.read_until(b'\n', &mut y).unwrap();
        if x != y {
            return Some(format!(
                "line {line}: expected {:?}; actual {:?}",
                String::from_utf8_lossy(&x),
                String::from_utf8_lossy(&y)
            ));
        }
        if x.is_empty() {
            return None;
        }
    }
    unreachable!()
}

#[test]
fn pbmc1k_all_files_byte_identical() {
    let fixture = std::env::var_os("NOMICS_SOLO_FIXTURE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(
                std::env::var("NOMICS_FIXTURES")
                    .unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into()),
            )
            .join("starsolo-pbmc1k")
        });
    if !fixture.join("ORACLE-READY").exists() {
        eprintln!("SKIP: oracle fixture not ready at {}", fixture.display());
        return;
    }
    let output = std::env::temp_dir().join(format!("nomics-solo-gate-{}", std::process::id()));
    let features = fixture.join("golden-Solo.out/Gene/raw/features.tsv");
    let features = if features.exists() {
        features
    } else {
        features.with_extension("tsv.gz")
    };
    let config = nomics_solo::Config {
        bam: fixture.join("oracle.bam"),
        whitelist: fixture.join("10x_V3_barcode_whitelist.txt.gz"),
        features,
        output: output.clone(),
        threads: 6,
    };
    let report = nomics_solo::count(&config).unwrap();
    eprintln!("{report:?}\nOutput: {}", output.display());
    let mut differences = Vec::new();
    for file in FILES {
        if let Some(diff) = first_diff(
            &fixture.join("golden-Solo.out").join(file),
            &output.join(file),
        ) {
            eprintln!("{file}: {diff}");
            differences.push(format!("{file}: {diff}"));
        } else {
            eprintln!("{file}: IDENTICAL");
        }
    }
    assert!(
        differences.is_empty(),
        "{}\nOutput retained at {}",
        differences.join("\n"),
        output.display()
    );
    // Retain the verified artifacts for inspection; caller owns the temp directory.
}
