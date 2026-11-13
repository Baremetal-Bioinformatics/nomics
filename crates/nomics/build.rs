//! Record the commit the binary is built from, for `nomics verify` reports.
//! A build from a source archive (no `.git`) records nothing, and the report
//! says `unknown`.

fn main() {
    let commit = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned())
        .filter(|c| c.len() == 40);
    if let Some(c) = commit {
        let dirty = std::process::Command::new("git")
            .args(["status", "--porcelain", "--untracked-files=no"])
            .output()
            .ok()
            .is_some_and(|o| !o.stdout.is_empty());
        println!("cargo:rustc-env=NOMICS_COMMIT={c}{}", if dirty { "-dirty" } else { "" });
    }
    // Rerun on anything that can change the commit or the dirty flag: HEAD
    // and refs (commits, branch switches), the index (staging), and the
    // tracked sources themselves, since an unstaged edit touches neither.
    // Without the source paths a build after an edit kept the last clean
    // stamp. Directories are scanned recursively by mtime; the script itself
    // is two git calls.
    for p in [
        "../../.git/HEAD",
        "../../.git/index",
        "../../.git/refs/heads",
        "../../.git/packed-refs",
        "../../crates",
        "../../Cargo.toml",
        "../../Cargo.lock",
        "../../reference",
    ] {
        println!("cargo:rerun-if-changed={p}");
    }
}
