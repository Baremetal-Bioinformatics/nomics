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
//! Reports are collected in the existing filtering pass and written natively.
//! A JSON request also writes the sibling `.log`; stderr keeps the shim marker.
//! No real fastp is needed for the supported report path.

use std::io::Write;
use std::path::{Path, PathBuf};

use super::args::{Kind, Spec, Unsupported, parse, spec};
use super::{Decision, threads};

const FASTP: &[Spec] = &[
    spec(&["--in1", "-i"], Kind::Value, "in1"),
    spec(&["--in2", "-I"], Kind::Value, "in2"),
    spec(&["--out1", "-o"], Kind::Value, "out1"),
    spec(&["--out2", "-O"], Kind::Value, "out2"),
    spec(&["--json", "-j"], Kind::Value, "json"),
    spec(&["--html", "-h"], Kind::Value, "html"),
    spec(&["--thread", "-w"], Kind::Value, "thread"),
    spec(
        &["--detect_adapter_for_pe"],
        Kind::Flag,
        "detect_adapter_for_pe",
    ),
    spec(
        &["--disable_adapter_trimming", "-A"],
        Kind::Flag,
        "disable_adapter_trimming",
    ),
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
    json: Option<PathBuf>,
    html: Option<PathBuf>,
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
    if args
        .iter()
        .any(|a| a == "--version" || a == "-v" || a == "--help" || a == "-?")
    {
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
            return Err(Unsupported::new(format!(
                "{what} {name} (uncompressed output)"
            )));
        }
    }
    for i in [&in1, &in2] {
        if !i.to_string_lossy().ends_with(".gz") {
            return Err(Unsupported::new(format!(
                "{} (not gzip input)",
                i.display()
            )));
        }
    }
    Ok(Plan {
        in1,
        in2,
        outs: [out1, out2],
        fastp_threads,
        json: p.one("json").map(PathBuf::from),
        html: p.one("html").map(PathBuf::from),
        all_args: args.to_vec(),
    })
}

/// `<dir>/<NNNN>.<basename>` for chunk `n` (1-based), as fastp names split output.
pub fn chunk_path(out: &Path, n: usize) -> PathBuf {
    let base = out
        .file_name()
        .map(|b| b.to_string_lossy().into_owned())
        .unwrap_or_default();
    let name = format!("{n:04}.{base}");
    match out.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.join(name),
        _ => PathBuf::from(name),
    }
}

/// Scratch for the uncompressed chunks, beside the outputs (the task's work directory).
fn spill_dir(out: &Path) -> PathBuf {
    let name = format!(".nomics-fastp-{}", std::process::id());
    match out.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.join(name),
        _ => PathBuf::from(name),
    }
}

fn run(p: Plan) -> Result<(), String> {
    let t = std::time::Instant::now();
    let start_seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    let mut report = if p.json.is_some() || p.html.is_some() {
        Some(nomics_prep::report::Report::new().map_err(|e| e.to_string())?)
    } else {
        None
    };
    // Chunks go to disk as they are filtered, uncompressed, and are compressed from there:
    // the task has the memory sarek gives fastp (4 GB), not room for the whole input.
    let spill = spill_dir(&p.outs[0]);
    let result =
        nomics_prep::run_to_dir_with_report(&p.in1, &p.in2, p.fastp_threads, &spill, report.as_mut());
    let stats = match result {
        Ok(x) => x,
        Err(e) if e.to_string().contains("--split_by_lines") => {
            let _ = std::fs::remove_dir_all(&spill);
            super::say(format!("nomics-shim: delegated fastp (unsupported: {e})"));
            let args: Vec<std::ffi::OsString> = p.all_args.iter().map(Into::into).collect();
            std::process::exit(match super::delegate("fastp", &args) {
                c if c == std::process::ExitCode::SUCCESS => 0,
                _ => 127,
            });
        }
        Err(e) => {
            let _ = std::fs::remove_dir_all(&spill);
            return Err(e.to_string());
        }
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
    for c in 0..p.fastp_threads {
        for (mate, src) in nomics_prep::spill_paths(&spill, c).iter().enumerate() {
            let path = chunk_path(&p.outs[mate], c + 1);
            let f = std::fs::File::create(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let mut w =
                nomics_io::bgzf::ParallelBgzfWriter::new(std::io::BufWriter::new(f), nthreads);
            let mut text = std::fs::File::open(src).map_err(|e| format!("{}: {e}", src.display()))?;
            std::io::copy(&mut text, &mut w).map_err(|e| e.to_string())?;
            w.finish()
                .and_then(|mut f| f.flush())
                .map_err(|e| format!("{}: {e}", path.display()))?;
            let _ = std::fs::remove_file(src);
        }
    }
    let _ = std::fs::remove_dir_all(&spill);
    if let Some(report) = report {
        // argv[0] is significant provenance in fastp's hand-written reports.
        let argv0 = std::env::args().next().unwrap_or_else(|| "fastp".into());
        let command = format!("{argv0} {} ", p.all_args.join(" "));
        if let Some(path) = &p.json {
            std::fs::write(path, report.json(&stats, &command)).map_err(|e| e.to_string())?;
            let now = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs();
            let log = report.log(
                &stats,
                &command,
                &path.to_string_lossy(),
                &p.html
                    .as_deref()
                    .unwrap_or(Path::new("fastp.html"))
                    .to_string_lossy(),
                now.saturating_sub(start_seconds),
            );
            let log_path = if path.extension().is_some_and(|e| e == "json") {
                path.with_extension("log")
            } else {
                super::with_suffix(path, ".log")
            };
            std::fs::write(log_path, log).map_err(|e| e.to_string())?;
        }
        if let Some(path) = &p.html {
            let title = nomics_prep::report::timestamp_now();
            let footer = nomics_prep::report::timestamp_now();
            std::fs::write(path, report.html(&stats, &command, &title, &footer))
                .map_err(|e| e.to_string())?;
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
        assert_eq!(p.html, Some(PathBuf::from("a.html")));
        assert_eq!(p.json, Some(PathBuf::from("a.json")));
        assert_eq!(
            plan(&v(
                &sarek.replace("--length_required 15", "--length_required 30")
            ))
            .err()
            .unwrap()
            .0,
            "--length_required 30"
        );
        assert!(
            plan(&v(
                &sarek.replace("--disable_adapter_trimming", "--trim_poly_x")
            ))
            .is_err()
        );
        assert!(plan(&v("--version")).is_err());
    }
}
