use nomics_hc::jmath;

use std::io::{BufRead, BufReader};

fn bits(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn assert_java(got: f64, want: f64, what: &str) {
    assert!(
        (got.is_nan() && want.is_nan()) || got.to_bits() == want.to_bits(),
        "{what}: got {:016x}, want {:016x}",
        got.to_bits(),
        want.to_bits()
    );
}

#[test]
fn java_strictmath_fixture_is_bit_exact() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/java_strictmath.tsv.gz"
    );
    let reader = BufReader::new(nomics_io::gzip::GzReader::open(path).unwrap());
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        let x = bits(fields[1]);
        let (got, want) = match fields[0] {
            "log" => (jmath::log(x), bits(fields[2])),
            "log10" => (jmath::log10(x), bits(fields[2])),
            "log1p" => (jmath::log1p(x), bits(fields[2])),
            "exp" => (jmath::exp(x), bits(fields[2])),
            "expm1" => (jmath::expm1(x), bits(fields[2])),
            "pow" => (jmath::pow(x, bits(fields[2])), bits(fields[3])),
            _ => panic!("unknown fixture function at line {line_no}"),
        };
        assert_java(got, want, &format!("line {}: {line}", line_no + 1));
    }
}

#[test]
fn special_values_and_exact_powers() {
    assert_java(jmath::log(0.0), f64::NEG_INFINITY, "log +0");
    assert_java(jmath::log(-0.0), f64::NEG_INFINITY, "log -0");
    assert!(jmath::log(-1.0).is_nan());
    assert_java(jmath::log(f64::INFINITY), f64::INFINITY, "log inf");
    assert!(jmath::log(f64::NAN).is_nan());
    assert!(jmath::log(f64::from_bits(1)).is_finite());
    assert_java(jmath::log10(-0.0), f64::NEG_INFINITY, "log10 -0");
    assert_java(jmath::log10(1.0), 0.0, "log10 one");
    assert_java(jmath::log1p(-1.0), f64::NEG_INFINITY, "log1p -1");
    assert!(jmath::log1p(-1.000_000_000_000_000_2).is_nan());
    assert_java(jmath::log1p(f64::INFINITY), f64::INFINITY, "log1p inf");
    assert_java(jmath::expm1(-0.0), -0.0, "expm1 -0");
    assert_java(jmath::expm1(f64::NEG_INFINITY), -1.0, "expm1 -inf");
    assert_java(jmath::exp(f64::NEG_INFINITY), 0.0, "exp -inf");
    assert_java(jmath::exp(f64::INFINITY), f64::INFINITY, "exp inf");
    assert_java(jmath::pow(-0.0, 3.0), -0.0, "negative zero odd power");
    assert_java(
        jmath::pow(-0.0, -3.0),
        f64::NEG_INFINITY,
        "negative zero reciprocal",
    );
    assert!(jmath::pow(-2.0, 0.5).is_nan());
    assert_java(jmath::pow(2.0, 10.0), 1024.0, "exact positive power");
    assert_java(jmath::pow(-2.0, 11.0), -2048.0, "exact negative odd power");
}
