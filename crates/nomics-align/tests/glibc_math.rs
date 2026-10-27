use std::io::{BufRead, BufReader};
use std::path::Path;

use nomics_io::gzip::GzReader;

use nomics_align::glibc_math;

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn big_fixture() -> String {
    format!("{}/glibc_logexp_big.tsv.gz", fixtures_dir())
}

fn check_fixture(path: &Path) {
    let reader = GzReader::open(path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let reader = BufReader::with_capacity(1 << 20, reader);
    let mut checked = 0_usize;

    for (line_index, line) in reader.lines().enumerate() {
        let line = line.unwrap_or_else(|error| {
            panic!(
                "failed to read {} line {}: {error}",
                path.display(),
                line_index + 1
            )
        });
        let mut fields = line.split_ascii_whitespace();
        let function = fields.next().expect("fixture function");
        let input_bits = u64::from_str_radix(fields.next().expect("fixture input"), 16)
            .expect("hex fixture input");
        let expected_bits = u64::from_str_radix(fields.next().expect("fixture result"), 16)
            .expect("hex fixture result");
        assert!(
            fields.next().is_none(),
            "extra field on fixture line {}",
            line_index + 1
        );

        let input = f64::from_bits(input_bits);
        let actual = match function {
            "log" => glibc_math::log(input),
            "exp" => glibc_math::exp(input),
            other => panic!(
                "unknown function {other:?} on fixture line {}",
                line_index + 1
            ),
        };
        let expected = f64::from_bits(expected_bits);
        assert!(
            actual.to_bits() == expected_bits || (actual.is_nan() && expected.is_nan()),
            "{} line {}: {function}(bits=0x{input_bits:016x}) = \
             bits 0x{:016x}, expected 0x{expected_bits:016x}",
            path.display(),
            line_index + 1,
            actual.to_bits(),
        );
        checked += 1;
    }

    assert!(checked > 0, "empty fixture: {}", path.display());
}

#[test]
fn glibc_logexp_fixture() {
    check_fixture(Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/glibc_logexp.tsv.gz"
    )));
}

#[test]
#[ignore = "six-million-case fixture"]
fn glibc_logexp_big_fixture() {
    let fixture = big_fixture();
    let path = Path::new(&fixture);
    if path.try_exists().expect("check heavy fixture existence") {
        check_fixture(path);
    }
}
