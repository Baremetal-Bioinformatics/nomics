//! Run with --release for meaningful timings. Missing local fixtures skip cleanly.
use std::path::PathBuf;
use std::time::Instant;

#[test]
fn tier0_json_and_unchanged_chunks() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap_or_default());
    let root = std::env::var_os("NOMICS_FASTP_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics"));
    let a = root.join("data/tier0/tier0_1.fastq.gz");
    let b = root.join("data/tier0/tier0_2.fastq.gz");
    let golden = root
        .join("runs/tier0-wes-20261008-r1/results/reports/fastp/HG002_WES/HG002_WES-L1.fastp.json");
    if !a.is_file() || !b.is_file() || !golden.is_file() {
        eprintln!(
            "SKIP tier0 fastp report: inputs/golden absent under {}",
            root.display()
        );
        return;
    }
    let expected = std::fs::read_to_string(&golden).unwrap();
    let command = expected
        .lines()
        .find_map(|l| l.strip_prefix("\t\"command\": \""))
        .unwrap()
        .strip_suffix('"')
        .unwrap();
    let start = Instant::now();
    let (plain, plain_stats) = nomics_prep::run(&a, &b, 12).unwrap();
    let off = start.elapsed();
    let start = Instant::now();
    let mut report = nomics_prep::report::Report::new().unwrap();
    let (reported, stats) = nomics_prep::run_with_report(&a, &b, 12, &mut report).unwrap();
    let actual = report.json(&stats, command);
    let on = start.elapsed();
    assert_eq!(stats, plain_stats);
    assert_eq!(plain.r1, reported.r1);
    assert_eq!(plain.r2, reported.r2);
    eprintln!(
        "tier0 prep report off {off:?}; on {on:?}; added {:?}",
        on.saturating_sub(off)
    );
    let dir = std::env::temp_dir().join(format!("nomics-fastp-report-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let output = dir.join("tier0.fastp.json");
    std::fs::write(&output, &actual).unwrap();
    let cmp = std::process::Command::new("cmp")
        .arg(&golden)
        .arg(&output)
        .status()
        .unwrap();
    if !cmp.success() {
        // This golden's rate equals the serial one, so any difference fails here,
        // the rate included (docs/GATES.md exempts it only across other runs).
        // Say which part differs. Never parse/re-serialize the document.
        fn split(s: &str) -> (&str, &str, &str) {
            let (before, rest) = s
                .split_once("\t\"duplication\": {\n\t\t\"rate\": ")
                .unwrap();
            let (token, after) = rest.split_once('\n').unwrap();
            assert!(token.parse::<f64>().is_ok());
            assert!(!after.contains("\"duplication\""));
            (before, token, after)
        }
        let (ap, av, asuf) = split(&actual);
        let (ep, ev, esuf) = split(&expected);
        assert_eq!(
            ap.as_bytes(),
            ep.as_bytes(),
            "difference before duplication.rate; raw output {}",
            output.display()
        );
        assert_eq!(
            asuf.as_bytes(),
            esuf.as_bytes(),
            "difference after duplication.rate; raw output {}",
            output.display()
        );
        panic!(
            "duplication.rate serial={av} golden={ev}; every other byte identical; raw output {}. This golden matched exactly before, so the serial count has changed.",
            output.display()
        );
    }
    let html_golden = golden.with_extension("html");
    let log_golden = golden.with_extension("log");
    if html_golden.is_file() && log_golden.is_file() {
        let expected_html = std::fs::read_to_string(&html_golden).unwrap();
        let title_time = expected_html
            .split_once("<title>fastp report at ")
            .unwrap()
            .1
            .split_once(" </title>")
            .unwrap()
            .0;
        let footer_time = expected_html
            .split_once("fastp 1.1.0, at ")
            .unwrap()
            .1
            .split_once(" </div>")
            .unwrap()
            .0;
        let html = report.html(&stats, command, title_time, footer_time);
        let html_path = dir.join("tier0.fastp.html");
        std::fs::write(&html_path, &html).unwrap();
        let equal = std::process::Command::new("cmp")
            .arg(&html_golden)
            .arg(&html_path)
            .status()
            .unwrap()
            .success();
        assert!(equal, "HTML differs; raw output {}", html_path.display());
        let expected_log = std::fs::read_to_string(&log_golden).unwrap();
        let json_path = expected_log
            .split_once("JSON report: ")
            .unwrap()
            .1
            .lines()
            .next()
            .unwrap();
        let html_path = expected_log
            .split_once("HTML report: ")
            .unwrap()
            .1
            .lines()
            .next()
            .unwrap();
        let seconds = expected_log
            .split_once("time used: ")
            .unwrap()
            .1
            .split_once(' ')
            .unwrap()
            .0
            .parse()
            .unwrap();
        let log = report.log(&stats, command, json_path, html_path, seconds);
        let log_path = dir.join("tier0.fastp.log");
        std::fs::write(&log_path, &log).unwrap();
        let equal = std::process::Command::new("cmp")
            .arg(&log_golden)
            .arg(&log_path)
            .status()
            .unwrap()
            .success();
        assert!(equal, "log differs; raw output {}", log_path.display());
    }
}
