//! Compare dense-SA intervals 4, 2, and 1 on random lookups and real reads.
//!
//! usage: sa_bench <index> <R1.fq.gz> <R2.fq.gz> <pairs>
//!                 [lookups=10000000] [verify=1000000] [threads=20]
//!                 [cache-dir=$HOME/omics/cache]

use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

use nomics_align::chain::{Chain, mem_chain};
use nomics_align::index::Index;
use nomics_align::mem::read_fastq;
use nomics_align::seed::{SmemAux, collect_intv};
use nomics_io::gzip::GzReader;
use nomics_mem::HostBuf;

const RANDOM_SEED: u64 = 0x9e37_79b9_7f4a_7c15;

#[inline]
fn next_random(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

fn cache_path(dir: &Path, interval: u64) -> PathBuf {
    dir.join(format!("sa{interval}.GRCh38.nsa"))
}

fn read_sequences(r1: &Path, r2: &Path, pairs: usize) -> io::Result<Vec<Vec<u8>>> {
    let mut reads = Vec::with_capacity(2 * pairs);
    let mut line = Vec::new();
    for path in [r1, r2] {
        let mut reader = BufReader::with_capacity(1 << 20, GzReader::open(path)?);
        for _ in 0..pairs {
            let read = read_fastq(&mut reader, &mut line)?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!("{} contains fewer than {pairs} reads", path.display()),
                )
            })?;
            reads.push(read.seq);
        }
    }
    Ok(reads)
}

fn fill_random_rows(rows: &mut HostBuf, last_row: u64) {
    let mut state = RANDOM_SEED;
    for row in rows.as_u64s_mut() {
        *row = next_random(&mut state) % (last_row + 1);
    }
}

fn time_chains(idx: &Index, reads: &[Vec<u8>]) -> (f64, Vec<Vec<Chain>>) {
    // Match align_profile's subtraction approach: measure collect_intv alone,
    // then measure mem_chain (which includes the same seed search) and subtract.
    let mut aux = SmemAux::default();
    let start = Instant::now();
    for read in reads {
        collect_intv(&idx.bwt, read, &mut aux);
    }
    let seed_seconds = start.elapsed().as_secs_f64();

    let start = Instant::now();
    let chains = reads
        .iter()
        .map(|read| mem_chain(idx, read, &mut aux))
        .collect();
    let mem_chain_seconds = start.elapsed().as_secs_f64();
    ((mem_chain_seconds - seed_seconds).max(0.0), chains)
}

fn main() -> io::Result<()> {
    let args: Vec<_> = std::env::args_os().collect();
    if !(5..=9).contains(&args.len()) {
        eprintln!(
            "usage: sa_bench <index> <R1.fq.gz> <R2.fq.gz> <pairs> [lookups=10000000] [verify=1000000] [threads=20] [cache-dir=$HOME/omics/cache]"
        );
        std::process::exit(2);
    }
    let parse = |i: usize, default: usize| -> io::Result<usize> {
        match args.get(i) {
            Some(value) => value
                .to_str()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid number")),
            None => Ok(default),
        }
    };
    let pairs = parse(4, 0)?;
    let lookups = parse(5, 10_000_000)?;
    let verify = parse(6, 1_000_000)?;
    let threads = parse(7, 20)?;
    if pairs == 0 || lookups == 0 || verify == 0 || threads == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "pairs, lookups, verify, and threads must all be positive",
        ));
    }
    let cache_dir = match args.get(8) {
        Some(path) => PathBuf::from(path),
        None => PathBuf::from(std::env::var_os("HOME").ok_or_else(|| {
            io::Error::new(io::ErrorKind::NotFound, "HOME is not set; pass cache-dir")
        })?)
        .join("omics/cache"),
    };

    let mut idx = Index::load(Path::new(&args[1]))?;
    let reads = read_sequences(Path::new(&args[2]), Path::new(&args[3]), pairs)?;
    let row_count = lookups.max(verify);
    let mut rows = HostBuf::new(row_count.checked_mul(8).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, "random row count is too large")
    })?)?;
    fill_random_rows(&mut rows, idx.bwt.seq_len);
    let rows = rows.as_u64s();
    let mut expected = HostBuf::new(verify.checked_mul(8).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "verification count is too large",
        )
    })?)?;
    let expected = expected.as_u64s_mut();
    let mut baseline_chains: Option<Vec<Vec<Chain>>> = None;

    println!(
        "{} reads, {lookups} timed lookups, {verify} verified rows, {threads} build threads",
        reads.len()
    );
    for interval in [4u64, 2, 1] {
        // Drop the preceding cache before allocating/building the next one.
        idx.dense = None;
        let path = cache_path(&cache_dir, interval);
        let start = Instant::now();
        idx.attach_dense_sa_with_interval(&path, interval, threads)?;
        let ready_seconds = start.elapsed().as_secs_f64();
        let dense = idx.dense.as_ref().unwrap();
        assert_eq!(dense.interval(), interval);
        let (memory, huge) = idx.memory();

        let start = Instant::now();
        for (i, &row) in rows[..verify].iter().enumerate() {
            let value = idx.sa_at(row);
            if interval == 4 {
                let bwa = idx.sa.sa(&idx.bwt, row);
                assert_eq!(value, bwa, "interval 4 differs from bwa SA at row {row}");
                expected[i] = value;
            } else {
                assert_eq!(
                    value, expected[i],
                    "interval {interval} differs from interval 4 at row {row}"
                );
            }
        }
        let verify_seconds = start.elapsed().as_secs_f64();

        let start = Instant::now();
        let mut checksum = 0u64;
        for &row in &rows[..lookups] {
            checksum = checksum.wrapping_add(std::hint::black_box(idx.sa_at(row)));
        }
        std::hint::black_box(checksum);
        let lookup_seconds = start.elapsed().as_secs_f64();

        let (chain_seconds, chains) = time_chains(&idx, &reads);
        if let Some(baseline) = &baseline_chains {
            assert_eq!(
                &chains, baseline,
                "interval {interval} produced different chains"
            );
        } else {
            baseline_chains = Some(chains);
        }

        println!(
            "interval {interval}: cache {} ({:.3} GiB), ready {ready_seconds:.2}s, index memory {:.3} GiB (huge {:?}); verify {verify_seconds:.2}s; random {:.3} us/lookup; chain(SA+tree) {:.3} us/read; checksum {checksum:016x}",
            path.display(),
            dense.bytes() as f64 / (1u64 << 30) as f64,
            memory as f64 / (1u64 << 30) as f64,
            huge.map(|n| format!("{:.3} GiB", n as f64 / (1u64 << 30) as f64)),
            lookup_seconds * 1e6 / lookups as f64,
            chain_seconds * 1e6 / reads.len() as f64,
        );
    }
    println!("all intervals returned identical SA values and identical chains");
    Ok(())
}
