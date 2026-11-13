//! Tier0 end-to-end CLI and multi-call shim report gate. No external fastp.
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::Command;

fn replace_value(text: &str, start: &str, end: &str) -> String {
    let (before, rest) = text.split_once(start).unwrap();
    let (_, after) = rest.split_once(end).unwrap();
    format!("{before}{start}<wall-clock>{end}{after}")
}

fn compare(dir: &Path, golden: &Path) {
    let name = "HG002_WES-L1.fastp";
    assert!(
        Command::new("cmp")
            .arg(dir.join(format!("{name}.json")))
            .arg(golden.join(format!("{name}.json")))
            .status()
            .unwrap()
            .success()
    );
    let read = |d: &Path, ext| std::fs::read_to_string(d.join(format!("{name}.{ext}"))).unwrap();
    let html = |d: &Path| {
        let s = replace_value(&read(d, "html"), "<title>fastp report at ", " </title>");
        replace_value(&s, "fastp 1.1.0, at ", " </div>")
    };
    assert_eq!(html(dir), html(golden), "HTML outside named timestamps");
    let log = |d: &Path| replace_value(&read(d, "log"), "fastp v1.1.0, time used: ", " seconds");
    assert_eq!(log(dir), log(golden), "log outside elapsed integer");
}

#[test]
fn tier0_native_cli_and_shim_reports() {
    let root = std::env::var_os("NOMICS_FASTP_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap_or_default()).join("omics")
        });
    let inputs: Vec<_> = (1..=2)
        .map(|m| root.join(format!("data/tier0/tier0_{m}.fastq.gz")))
        .collect();
    let golden = root.join("runs/tier0-wes-20261008-r1/results/reports/fastp/HG002_WES");
    if inputs.iter().any(|p| !p.is_file())
        || ["json", "html", "log"]
            .iter()
            .any(|ext| !golden.join(format!("HG002_WES-L1.fastp.{ext}")).is_file())
    {
        eprintln!(
            "SKIP tier0 CLI fastp reports: fixtures absent under {}",
            root.display()
        );
        return;
    }
    let expected = std::fs::read_to_string(golden.join("HG002_WES-L1.fastp.json")).unwrap();
    let command = expected
        .lines()
        .find_map(|l| l.strip_prefix("\t\"command\": \""))
        .unwrap()
        .strip_suffix('"')
        .unwrap();
    let base = std::env::temp_dir().join(format!("nomics-fastp-cli-{}", std::process::id()));
    for mode in ["cli", "shim"] {
        let dir = base.join(mode);
        std::fs::create_dir_all(&dir).unwrap();
        for (i, input) in inputs.iter().enumerate() {
            std::os::unix::fs::symlink(
                input,
                dir.join(format!("HG002_WES-L1_R{}.fastq.gz", i + 1)),
            )
            .unwrap();
        }
        let mut cmd = if mode == "cli" {
            let mut cmd = Command::new(env!("CARGO_BIN_EXE_nomics"));
            cmd.args([
                "prep",
                "--fastq1",
                "HG002_WES-L1_R1.fastq.gz",
                "--fastq2",
                "HG002_WES-L1_R2.fastq.gz",
                "--outdir",
                ".",
                "--prefix",
                "HG002_WES-L1",
                "--json",
                "HG002_WES-L1.fastp.json",
                "--html",
                "HG002_WES-L1.fastp.html",
                "--log",
                "HG002_WES-L1.fastp.log",
                "--report-command",
                command,
            ]);
            cmd
        } else {
            let shim = dir.join("fastp");
            std::os::unix::fs::symlink(env!("CARGO_BIN_EXE_nomics"), &shim).unwrap();
            let mut cmd = Command::new(shim);
            cmd.arg0("fastp").args(command.split_whitespace().skip(1));
            cmd.env("NOMICS_SHIM_REAL_PATH", dir.join("no-real-tools"));
            cmd
        };
        let result = cmd.current_dir(&dir).output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        compare(&dir, &golden);
    }
}
