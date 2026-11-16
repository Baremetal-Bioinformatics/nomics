//! Shared fixture decisions; included only by tests, with no dependency graph edges.
use std::fs::OpenOptions;
use std::io::Write;

/// Record a fixture decision before running its comparison. `wanted` must explain
/// the location or environment variable that supplies the fixture/tool.
/// Counts are decisions (a parameterized test can make several), not assertions.
pub fn check(present: bool, wanted: impl std::fmt::Display) -> bool {
    let wanted = wanted.to_string();
    let thread = std::thread::current();
    // Integration targets can contain identically named tests in one package.
    let test = format!(
        "{}::{}",
        env!("CARGO_CRATE_NAME"),
        thread.name().unwrap_or("unnamed-test")
    );
    let status = if present { "executed" } else { "skipped" };
    if let Some(path) = std::env::var_os("NOMICS_FIXTURE_LOG") {
        let escape = |s: &str| {
            s.replace('\\', "\\\\")
                .replace('\t', "\\t")
                .replace('\n', "\\n")
                .replace('\r', "\\r")
        };
        let line = format!(
            "{}\t{}\t{status}\t{}\n",
            escape(&test),
            env!("CARGO_PKG_NAME"),
            escape(&wanted)
        );
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .unwrap_or_else(|e| panic!("open NOMICS_FIXTURE_LOG {path:?}: {e}"));
        // One append write, also across independent libtest processes. Never
        // retry a short write as a second fragment that could interleave.
        let written = file
            .write(line.as_bytes())
            .expect("append NOMICS_FIXTURE_LOG");
        assert_eq!(written, line.len(), "short NOMICS_FIXTURE_LOG append");
    }
    if !present {
        assert!(
            std::env::var_os("NOMICS_REQUIRE_FIXTURES").as_deref()
                != Some(std::ffi::OsStr::new("1")),
            "required fixture/oracle/binary missing: {wanted}. Provide the named files/tool or configure the indicated environment variable."
        );
        eprintln!("SKIP: {wanted}");
    }
    present
}
