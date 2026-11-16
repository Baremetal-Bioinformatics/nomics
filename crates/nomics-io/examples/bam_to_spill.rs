//! Rewrite a run's chunk BAMs as light spill files (`nomics_io::spill`): the
//! same decoded stream, `<name>.bam` -> `<name>.nspill` next to it, then a
//! spill manifest with the BAM manifest's run identity and order. For timing
//! markdup on a spill with `nomics run --resume-from markdup`: convert a
//! run's `mapped/` in place, then move the BAMs (and `bam.manifest`) aside (a
//! mapped dir holding both is refused, and so is a spill set the manifest
//! does not describe exactly).
//!
//! `cargo run --release -p nomics-io --example bam_to_spill -- [--threads N] <mapped-dir>`

use std::path::{Path, PathBuf};
use std::time::Instant;

use nomics_io::spill::{BAM_MANIFEST, MANIFEST, check_set, read_manifest, write_manifest};

fn main() -> std::io::Result<()> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut threads = 8;
    if args.first().is_some_and(|a| a == "--threads") {
        threads = args[1].parse().expect("--threads N");
        args.drain(..2);
    }
    let [dir] = &args[..] else {
        panic!("give the run's mapped directory (with its {BAM_MANIFEST})");
    };
    let dir = Path::new(dir);
    // The BAM set as the run wrote it: every listed file, at its size.
    let manifest = read_manifest(dir, BAM_MANIFEST)?;
    let identity = manifest
        .identity
        .clone()
        .expect("a BAM manifest records its run");
    let listed: Vec<PathBuf> = manifest.files.iter().map(|(n, _)| dir.join(n)).collect();
    let bams = check_set(dir, BAM_MANIFEST, &manifest, &listed)?;
    let mut outs = Vec::new();
    for bam in &bams {
        let started = Instant::now();
        let out = bam.with_extension(nomics_io::spill::EXTENSION);
        let file = std::io::BufReader::with_capacity(1 << 22, std::fs::File::open(bam)?);
        let mut r = nomics_io::bgzf::ParallelBgzfReader::new(file, threads);
        let mut w = nomics_io::spill::SpillFile::create(&out, threads)?;
        let raw = std::io::copy(&mut r, &mut w)?;
        w.finish()?;
        let (a, b) = (
            std::fs::metadata(bam)?.len(),
            std::fs::metadata(&out)?.len(),
        );
        println!(
            "{} -> {}: {raw} bytes decoded, {a} -> {b} bytes ({:.2}x the BAM), {:.1}s",
            bam.display(),
            out.display(),
            b as f64 / a as f64,
            started.elapsed().as_secs_f64()
        );
        outs.push(out);
    }
    write_manifest(dir, MANIFEST, &identity, &outs)?;
    println!("{}: {} files", dir.join(MANIFEST).display(), outs.len());
    Ok(())
}
