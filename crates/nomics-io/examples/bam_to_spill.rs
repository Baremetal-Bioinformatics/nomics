//! Rewrite BAMs as light spill files (`nomics_io::spill`): the same decoded
//! stream, `<name>.bam` -> `<name>.nspill` next to it, then a spill manifest
//! naming them in their directory (all inputs must share one). For timing
//! markdup on a spill with `nomics run --resume-from markdup`: convert a
//! run's `mapped/*.bam` in place, then move the BAMs aside (a mapped dir
//! holding both is refused, and so is a spill set the manifest does not
//! describe exactly).
//!
//! `cargo run --release -p nomics-io --example bam_to_spill -- [--threads N] <file.bam>...`

use std::path::PathBuf;
use std::time::Instant;

fn main() -> std::io::Result<()> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut threads = 8;
    if args.first().is_some_and(|a| a == "--threads") {
        threads = args[1].parse().expect("--threads N");
        args.drain(..2);
    }
    let bams: Vec<PathBuf> = args.iter().map(PathBuf::from).collect();
    let dir = match bams.first().and_then(|b| b.parent()) {
        Some(d) if bams.iter().all(|b| b.parent() == Some(d)) => d.to_path_buf(),
        _ => panic!("give one or more BAMs, all in one directory"),
    };
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
    nomics_io::spill::write_manifest(&dir, &outs)?;
    println!(
        "{}: {} files",
        dir.join(nomics_io::spill::MANIFEST).display(),
        outs.len()
    );
    Ok(())
}
