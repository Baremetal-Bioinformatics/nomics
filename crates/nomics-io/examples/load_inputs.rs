//! Load a FASTA and stream VCFs, reporting counts and timings.
//!
//! `cargo run --release -p nomics-io --example load_inputs -- <ref.fasta> <a.vcf.gz>...`

use std::path::Path;
use std::time::Instant;

use nomics_io::fasta::Fasta;
use nomics_io::vcf::for_each_record;

fn main() {
    let mut args = std::env::args().skip(1);
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let fasta = args
        .next()
        .expect("usage: load_inputs <ref.fasta> <vcf>...");
    let t = Instant::now();
    let f = Fasta::load(Path::new(&fasta), threads).unwrap();
    let bases: usize = (0..f.len()).map(|i| f.contig(i).len()).sum();
    let chr20 = f
        .index_of("chr20")
        .map(|i| String::from_utf8_lossy(&f.contig(i)[30_000_000..30_000_060]).into_owned());
    println!(
        "fasta: {} contigs, {bases} bases in {:.1}s; chr20:30,000,001-30,000,060 = {chr20:?}",
        f.len(),
        t.elapsed().as_secs_f64()
    );
    for vcf in args {
        let t = Instant::now();
        let (mut n, mut snv) = (0u64, 0u64);
        let header = for_each_record(Path::new(&vcf), threads, &mut |r| {
            n += 1;
            if r.ref_allele.len() == 1 && r.alt.len() == 1 {
                snv += 1;
            }
            Ok(())
        })
        .unwrap();
        println!(
            "{vcf}: {n} records ({snv} single-base REF/ALT), {} header lines, {:.1}s",
            header.len(),
            t.elapsed().as_secs_f64()
        );
    }
}
