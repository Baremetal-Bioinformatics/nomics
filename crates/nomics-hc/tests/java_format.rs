use std::io::{BufRead, BufReader};
use std::path::Path;

use nomics_hc::javafmt::{double_to_string, format_fixed, format_sci};
use nomics_io::gzip::GzReader;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/java_format.tsv.gz"
);
fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn heavy_fixture() -> String {
    format!("{}/java_format_400k.tsv.gz", fixtures_dir())
}

fn check_fixture(path: &Path) {
    let reader = GzReader::open(path).unwrap_or_else(|error| {
        panic!(
            "cannot open Java-format fixture {}: {error}",
            path.display()
        )
    });
    let mut rows = 0_usize;
    for (index, line) in BufReader::new(reader).lines().enumerate() {
        let line = line.unwrap_or_else(|error| {
            panic!("cannot read {} line {}: {error}", path.display(), index + 1)
        });
        let columns: Vec<_> = line.split('\t').collect();
        assert_eq!(
            columns.len(),
            6,
            "{} line {} has {} columns: {line:?}",
            path.display(),
            index + 1,
            columns.len()
        );
        let bits = u64::from_str_radix(columns[0], 16).unwrap_or_else(|error| {
            panic!(
                "{} line {} has invalid bits {:?}: {error}",
                path.display(),
                index + 1,
                columns[0]
            )
        });
        let value = f64::from_bits(bits);
        let actual = [
            format_fixed(value, 2),
            format_fixed(value, 3),
            format_fixed(value, 4),
            format_sci(value, 3),
            double_to_string(value),
        ];
        let names = ["%.2f", "%.3f", "%.4f", "%.3e", "Double.toString"];
        for column in 0..actual.len() {
            assert_eq!(
                actual[column],
                columns[column + 1],
                "{} line {}, bits=0x{bits:016x}, value={value:?}, column {} ({})",
                path.display(),
                index + 1,
                column + 2,
                names[column]
            );
        }
        rows += 1;
    }
    assert!(rows > 100_000, "fixture had only {rows} rows");
}

#[test]
fn documented_decimal_rounding_examples() {
    assert_eq!(format_fixed(1.005, 2), "1.01");
    assert_eq!(format_fixed(2.675, 2), "2.68");
    assert_eq!(format_fixed(9.995, 2), "10.00");
    assert_eq!(format_sci(9.9995, 3), "1.000e+01");
}

#[test]
fn signs_zero_and_exponents() {
    assert_eq!(format_fixed(-0.0, 2), "-0.00");
    assert_eq!(format_sci(-0.0, 3), "-0.000e+00");
    assert_eq!(double_to_string(-0.0), "-0.0");
    assert_eq!(format_sci(0.001, 3), "1.000e-03");
    assert_eq!(format_sci(10.0, 3), "1.000e+01");
}

#[test]
fn java_100k_oracle() {
    check_fixture(Path::new(FIXTURE));
}

#[test]
#[ignore = "400k-line extended Java oracle"]
fn java_400k_oracle() {
    let fixture = heavy_fixture();
    let path = Path::new(&fixture);
    if path.exists() {
        check_fixture(path);
    }
}
