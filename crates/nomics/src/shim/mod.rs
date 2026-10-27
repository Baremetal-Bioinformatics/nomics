//! Multi-call shims: `nomics` invoked as `fastp`, `bwa`, `gatk` or `samtools`
//! (a symlink in a shim `bin/` directory) parses argv exactly as that tool
//! would. **Exact or delegate**: if every subcommand, flag and value is in
//! the supported set, the nomics implementation runs (its outputs are gated
//! identical to the real tool's); otherwise the real tool is `exec`ed with
//! the untouched argv. Nothing is approximated.
//!
//! Every invocation logs one stderr line:
//! `nomics-shim: handled gatk HaplotypeCaller` or
//! `nomics-shim: delegated gatk VariantFiltration (unsupported: --filter-expression)`.
//!
//! The real tool is found on `NOMICS_SHIM_REAL_PATH` (a `PATH`-style list of
//! directories) or else on `PATH` with every directory that resolves the tool
//! name to this executable removed. See `docs/strategy/nfcore-dropin.md`.

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

mod args;
mod bwa;
mod fastp;
mod gatk;

pub use args::Unsupported;

/// Tools the multi-call binary answers to.
pub const TOOLS: [&str; 4] = ["fastp", "bwa", "gatk", "samtools"];

/// The shim tool name if `argv[0]`'s basename is one of [`TOOLS`].
pub fn tool_name(argv0: &OsStr) -> Option<&'static str> {
    let base = Path::new(argv0).file_name()?.to_str()?;
    TOOLS.iter().copied().find(|t| *t == base)
}

/// What a shim decided to do with one invocation.
pub enum Decision {
    /// Run our implementation; `label` names it in the log line
    /// (`gatk HaplotypeCaller`), `run` does the work.
    Handle {
        label: String,
        run: Box<dyn FnOnce() -> Result<(), String>>,
    },
    /// Exec the real tool; `label` names the call, `why` the first unsupported item.
    Delegate { label: String, why: Unsupported },
}

/// Entry point for a shim invocation: decide, log, then run or delegate.
pub fn main(tool: &'static str) -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().collect();
    let args = &argv[1..];
    // Arguments that are not UTF-8 are never in a supported set.
    let utf8: Option<Vec<String>> = args.iter().map(|a| a.to_str().map(str::to_owned)).collect();
    let decision = match utf8 {
        None => Decision::Delegate {
            label: tool.to_owned(),
            why: Unsupported::new("non-UTF-8 argument"),
        },
        Some(args) => match tool {
            "fastp" => fastp::decide(&args),
            "bwa" => bwa::decide(&args),
            "gatk" => gatk::decide(&args),
            "samtools" => Decision::Delegate {
                label: match args.first() {
                    Some(sub) => format!("samtools {sub}"),
                    None => "samtools".to_owned(),
                },
                why: Unsupported::new(match args.first() {
                    Some(_) => "no samtools subcommand is handled in v1",
                    None => "usage query",
                }),
            },
            _ => unreachable!("tool_name only returns TOOLS"),
        },
    };
    match decision {
        Decision::Handle { label, run } => {
            eprintln!("nomics-shim: handled {label}");
            match run() {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("nomics-shim: {label} failed: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Decision::Delegate { label, why } => {
            eprintln!("nomics-shim: delegated {label} (unsupported: {why})");
            delegate(tool, args)
        }
    }
}

/// Exec the real `tool` with `args` untouched. Returns only on failure.
pub(super) fn delegate(tool: &str, args: &[OsString]) -> ExitCode {
    let real = match find_real(tool) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("nomics-shim: cannot delegate {tool}: {e}");
            return ExitCode::from(127);
        }
    };
    use std::os::unix::process::CommandExt;
    let err = std::process::Command::new(&real)
        .arg0(tool)
        .args(args)
        .exec();
    eprintln!("nomics-shim: exec {}: {err}", real.display());
    ExitCode::from(126)
}

/// Locate the real tool: `NOMICS_SHIM_REAL_PATH`, else `PATH` minus every
/// directory whose `tool` resolves to this executable.
pub fn find_real(tool: &str) -> Result<PathBuf, String> {
    let me = std::env::current_exe()
        .and_then(std::fs::canonicalize)
        .map_err(|e| format!("current_exe: {e}"))?;
    let (var, dirs) = match std::env::var_os("NOMICS_SHIM_REAL_PATH") {
        Some(v) => ("NOMICS_SHIM_REAL_PATH", v),
        None => ("PATH", std::env::var_os("PATH").unwrap_or_default()),
    };
    for dir in std::env::split_paths(&dirs) {
        if dir.as_os_str().is_empty() {
            continue;
        }
        let cand = dir.join(tool);
        let Ok(meta) = std::fs::metadata(&cand) else {
            continue;
        };
        use std::os::unix::fs::PermissionsExt;
        if !meta.is_file() || meta.permissions().mode() & 0o111 == 0 {
            continue;
        }
        // A shim directory (or any other link to us): skip it, never exec ourselves.
        if std::fs::canonicalize(&cand).is_ok_and(|c| c == me) {
            continue;
        }
        return Ok(cand);
    }
    Err(format!("no real {tool} on {var} (shim directories excluded)"))
}

/// Threads to use for our implementation: the machine's parallelism. The
/// tools' own thread flags do not change their outputs, so they are accepted
/// and only bound the pool when smaller.
pub fn threads(requested: Option<usize>) -> usize {
    let n = std::thread::available_parallelism().map_or(4, |n| n.get());
    match requested {
        Some(r) if r > 0 => r.min(n).max(1),
        _ => n,
    }
}

/// Path helper: `p` with `ext` appended to its file name.
pub fn with_suffix(p: &Path, ext: &str) -> PathBuf {
    let mut s = p.as_os_str().to_owned();
    s.push(ext);
    PathBuf::from(s)
}
