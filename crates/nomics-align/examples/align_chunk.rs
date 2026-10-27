//! Align one fastp chunk to a coordinate-sorted BAM (bwa mem | samtools sort).
//!
//! usage: align_chunk <index-hint> <R1.fastq.gz> <R2.fastq.gz> <out.bam> [threads] [dense-sa-cache]
//!
//! With a dense-SA cache path (built on first use), SA lookups use it.
//! NOMICS_GPU_SEED=1 runs seed search on the GPU (build with --features cuda).
//! NOMICS_GPU_EXTEND=fused runs the fused extension path (the CPU stand-in
//! without cuda). NOMICS_TIMING=1 also prints the extension counters.
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use nomics_align::index::Index;
use nomics_align::mem::align_stream;
use nomics_align::samsort::{header, sort_sam};
use nomics_io::gzip::GzReader;

/// The `DS:` value is part of the golden output bytes; override with `NOMICS_RG_DS`.
fn rg() -> String {
    let ds = std::env::var("NOMICS_RG_DS").unwrap_or_else(|_| {
        let reference = std::env::var("NOMICS_REF").unwrap_or_else(|_| {
            format!(
                "{}/omics/ref/GATK.GRCh38",
                std::env::var("HOME").unwrap_or_default()
            )
        });
        format!("{reference}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta")
    });
    format!("@RG\tID:HG002_WES.L1\tPU:L1\tSM:HG002_HG002_WES\tLB:HG002_WES\tDS:{ds}\tPL:ILLUMINA")
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let threads: usize = a.get(5).map_or(20, |s| s.parse().unwrap());
    let gpu = std::env::var_os("NOMICS_GPU_SEED").is_some();
    let timing = std::env::var_os("NOMICS_TIMING").is_some();
    nomics_align::extend::ext_stats::ON.store(timing, std::sync::atomic::Ordering::Relaxed);
    let t = Instant::now();
    let mut idx = Index::load(Path::new(&a[1])).unwrap();
    if let Some(p) = a.get(6).map(std::ffi::OsString::from) {
        let t = Instant::now();
        idx.attach_dense_sa(Path::new(&p), 20).unwrap();
        eprintln!("dense SA ready in {:.1}s", t.elapsed().as_secs_f64());
    }
    let t_load = t.elapsed();
    let open = |p: &str| BufReader::with_capacity(1 << 20, GzReader::open(Path::new(p)).unwrap());
    let (mut r1, mut r2) = (open(&a[2]), open(&a[3]));
    let mut sam = Vec::new();
    let (first_pair, _) = align_stream(
        &idx,
        b"HG002_WES.L1",
        &mut r1,
        &mut r2,
        threads,
        gpu,
        &mut sam,
    )
    .unwrap();
    let t_align = t.elapsed();
    let set = sort_sam(header(&idx, &rg(), &a[1..].join(" ")), &sam, threads).unwrap();
    let t_sort = t.elapsed();
    set.write(Path::new(&a[4]), threads).unwrap();
    eprintln!(
        "load {:.1}s, align {:.1}s, sort {:.1}s, write {:.1}s; {} pairs, {} records",
        t_load.as_secs_f64(),
        (t_align - t_load).as_secs_f64(),
        (t_sort - t_align).as_secs_f64(),
        (t.elapsed() - t_sort).as_secs_f64(),
        first_pair,
        set.len()
    );
    if timing {
        eprintln!(
            "{}",
            nomics_align::extend::ext_stats::report(set.len() as u64)
        );
    }
}
