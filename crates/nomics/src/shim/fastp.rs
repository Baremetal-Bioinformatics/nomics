//! `fastp` 1.1.0 as sarek runs it:
//! `--in1 --in2 --out1 --out2 [--json] [--html] --thread N
//! [--detect_adapter_for_pe] --disable_adapter_trimming
//! --split_by_lines 200000000 [--length_required 15]`.
//! Options are parsed as fastp's `cmdline.h` (getopt_long) does: long names
//! exactly (no abbreviations), `--name=value`, short options with the value
//! attached or separate.
//!
//! The filtered chunks are ours (`nomics_prep`, gated on decompressed
//! content): with `--split_by_lines`, worker `t` writes
//! `<dir of out1>/<t+1 as 4 digits>.<basename of out1>` (and the same for
//! out2). If a chunk would exceed the line limit (fastp would roll over to a
//! further file, which is not implemented), the call is delegated.
//!
//! The JSON/HTML QC reports are not reproduced: when `--json`/`--html` are
//! given they are written by the real fastp, run alongside in report-only
//! mode (same inputs and filter options, no outputs). Without a real fastp
//! they are not written, and the shim says so.

use std::io::Write;
use std::path::{Path, PathBuf};

use super::args::{Kind, Spec, Unsupported, parse, spec};
use super::{Decision, find_real, threads};

const FASTP: &[Spec] = &[
    spec(&["--in1", "-i"], Kind::Value, "in1"),
    spec(&["--in2", "-I"], Kind::Value, "in2"),
    spec(&["--out1", "-o"], Kind::Value, "out1"),
    spec(&["--out2", "-O"], Kind::Value, "out2"),
    spec(&["--json", "-j"], Kind::Value, "json"),
    spec(&["--html", "-h"], Kind::Value, "html"),
    spec(&["--thread", "-w"], Kind::Value, "thread"),
    spec(&["--detect_adapter_for_pe"], Kind::Flag, "detect_adapter_for_pe"),
    spec(&["--disable_adapter_trimming", "-A"], Kind::Flag, "disable_adapter_trimming"),
    spec(&["--split_by_lines", "-S"], Kind::Value, "split_by_lines"),
    spec(&["--length_required", "-l"], Kind::Value, "length_required"),
];

/// The only `--split_by_lines` `nomics_prep` implements.
const SPLIT_LINES: &str = "200000000";

struct Plan {
    in1: PathBuf,
    in2: PathBuf,
    outs: [PathBuf; 2],
    fastp_threads: usize,
    /// Report-only argv for the real fastp, when JSON/HTML were requested.
    report_args: Option<Vec<String>>,
    all_args: Vec<String>,
}

pub fn decide(args: &[String]) -> Decision {
    match plan(args) {
        Ok(p) => Decision::Handle {
            label: "fastp".into(),
            run: Box::new(move || run(p)),
        },
        Err(why) => Decision::Delegate {
            label: "fastp".into(),
            why,
        },
    }
}

fn plan(args: &[String]) -> Result<Plan, Unsupported> {
    if args.is_empty() {
        return Err(Unsupported::new("usage query"));
    }
    if args.iter().any(|a| a == "--version" || a == "-v" || a == "--help" || a == "-?") {
        return Err(Unsupported::new("version/help query"));
    }
    let p = parse(args, FASTP, true)?;
    let in1 = PathBuf::from(p.req("in1")?);
    let in2 = PathBuf::from(
        p.one("in2")
            .ok_or_else(|| Unsupported::new("single-end input (no --in2)"))?,
    );
    let out1 = PathBuf::from(
        p.one("out1")
            .ok_or_else(|| Unsupported::new("no --out1 (report-only run)"))?,
    );
    let out2 = PathBuf::from(p.one("out2").ok_or_else(|| Unsupported::new("no --out2"))?);
    if !p.has("disable_adapter_trimming") {
        return Err(Unsupported::new("adapter trimming enabled"));
    }
    match p.one("split_by_lines") {
        Some(SPLIT_LINES) => {}
        Some(v) => return Err(Unsupported::new(format!("--split_by_lines {v}"))),
        None => return Err(Unsupported::new("no --split_by_lines (unsplit output)")),
    }
    // fastp's default is 15, the value nomics_prep implements.
    p.only("length_required", &["15"])?;
    let fastp_threads = match p.one("thread") {
        Some(v) => match v.parse::<usize>() {
            Ok(n) if (1..=16).contains(&n) => n,
            _ => return Err(Unsupported::new(format!("--thread {v}"))),
        },
        None => return Err(Unsupported::new("no --thread (default worker count)")),
    };
    for (o, what) in [(&out1, "--out1"), (&out2, "--out2")] {
        let name = o.to_string_lossy();
        if !name.ends_with(".gz") {
            return Err(Unsupported::new(format!("{what} {name} (uncompressed output)")));
        }
    }
    for i in [&in1, &in2] {
        if !i.to_string_lossy().ends_with(".gz") {
            return Err(Unsupported::new(format!("{} (not gzip input)", i.display())));
        }
    }
    let report_args = (p.has("json") || p.has("html")).then(|| {
        let mut a = vec![
            "--in1".to_owned(),
            p.one("in1").unwrap().to_owned(),
            "--in2".to_owned(),
            p.one("in2").unwrap().to_owned(),
        ];
        for key in ["json", "html", "thread", "length_required"] {
            if let Some(v) = p.one(key) {
                a.push(format!("--{key}"));
                a.push(v.to_owned());
            }
        }
        for key in ["detect_adapter_for_pe", "disable_adapter_trimming"] {
            if p.has(key) {
                a.push(format!("--{key}"));
            }
        }
        a
    });
    Ok(Plan {
        in1,
        in2,
        outs: [out1, out2],
        fastp_threads,
        report_args,
        all_args: args.to_vec(),
    })
}

/// `<dir>/<NNNN>.<basename>` for chunk `n` (1-based), as fastp names split output.
pub fn chunk_path(out: &Path, n: usize) -> PathBuf {
    let base = out.file_name().map(|b| b.to_string_lossy().into_owned()).unwrap_or_default();
    let name = format!("{n:04}.{base}");
    match out.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.join(name),
        _ => PathBuf::from(name),
    }
}

fn run(p: Plan) -> Result<(), String> {
    // The QC reports come from the real fastp, run alongside.
    let mut reporter = None;
    if let Some(args) = &p.report_args {
        match find_real("fastp") {
            Ok(real) => {
                eprintln!("nomics-shim: fastp --json/--html delegated to {}", real.display());
                let child = std::process::Command::new(&real)
                    .args(args)
                    .stdout(std::process::Stdio::null())
                    .stderr(std::process::Stdio::null())
                    .spawn()
                    .map_err(|e| format!("spawn {}: {e}", real.display()))?;
                reporter = Some(child);
            }
            Err(e) => eprintln!("nomics-shim: fastp --json/--html not written: {e}"),
        }
    }
    let t = std::time::Instant::now();
    let (chunks, stats) = match nomics_prep::run(&p.in1, &p.in2, p.fastp_threads) {
        Ok(x) => x,
        Err(e) if e.to_string().contains("--split_by_lines") => {
            if let Some(mut c) = reporter {
                let _ = c.kill();
                let _ = c.wait();
            }
            eprintln!("nomics-shim: delegated fastp (unsupported: {e})");
            let args: Vec<std::ffi::OsString> = p.all_args.iter().map(Into::into).collect();
            std::process::exit(match super::delegate("fastp", &args) {
                c if c == std::process::ExitCode::SUCCESS => 0,
                _ => 127,
            });
        }
        Err(e) => return Err(e.to_string()),
    };
    // fastp's own summary goes to stderr; ours is shorter.
    eprintln!(
        "nomics-shim fastp: {} pairs in, {} pairs out ({} reads low quality, {} too many N, {} too short) in {:.1}s",
        stats.pairs_in,
        stats.pairs_out,
        stats.low_quality_reads,
        stats.too_many_n_reads,
        stats.too_short_reads,
        t.elapsed().as_secs_f64()
    );
    let nthreads = threads(None);
    for (mate, texts) in [&chunks.r1, &chunks.r2].into_iter().enumerate() {
        for (c, text) in texts.iter().enumerate() {
            let path = chunk_path(&p.outs[mate], c + 1);
            let f = std::fs::File::create(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let mut w =
                nomics_io::bgzf::ParallelBgzfWriter::new(std::io::BufWriter::new(f), nthreads);
            w.write_all(text).map_err(|e| e.to_string())?;
            w.finish()
                .and_then(|mut f| f.flush())
                .map_err(|e| format!("{}: {e}", path.display()))?;
        }
    }
    if let Some(mut c) = reporter {
        let status = c.wait().map_err(|e| format!("real fastp: {e}"))?;
        if !status.success() {
            return Err(format!("real fastp (reports) failed: {status}"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_names_follow_fastp() {
        assert_eq!(
            chunk_path(Path::new("HG002_WES-L1_R1.fastp.fastq.gz"), 1),
            PathBuf::from("0001.HG002_WES-L1_R1.fastp.fastq.gz")
        );
        assert_eq!(
            chunk_path(Path::new("out/x_R2.fq.gz"), 12),
            PathBuf::from("out/0012.x_R2.fq.gz")
        );
    }

    #[test]
    fn sarek_command_is_supported_and_others_delegate() {
        let v = |s: &str| s.split_whitespace().map(str::to_owned).collect::<Vec<_>>();
        let sarek = "--in1 a_R1.fastq.gz --in2 a_R2.fastq.gz --out1 a_R1.fastp.fastq.gz \
             --out2 a_R2.fastp.fastq.gz --json a.json --html a.html --thread 12 \
             --detect_adapter_for_pe --disable_adapter_trimming --split_by_lines 200000000 \
             --length_required 15";
        let p = plan(&v(sarek)).unwrap();
        assert_eq!(p.fastp_threads, 12);
        assert!(p.report_args.is_some());
        assert_eq!(
            plan(&v(&sarek.replace("--length_required 15", "--length_required 30")))
                .err()
                .unwrap()
                .0,
            "--length_required 30"
        );
        assert!(plan(&v(&sarek.replace("--disable_adapter_trimming", "--trim_poly_x"))).is_err());
        assert!(plan(&v("--version")).is_err());
    }
}
