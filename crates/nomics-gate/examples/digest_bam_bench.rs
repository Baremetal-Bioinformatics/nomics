//! Time `digest_bam` against `digest_bam_threads`: `digest_bam_bench BAM [THREADS...]`.
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("BAM path"));
    let threads: Vec<usize> = args.map(|a| a.parse().expect("thread count")).collect();
    for t in if threads.is_empty() { vec![1, 4] } else { threads } {
        let start = Instant::now();
        let d = nomics_gate::digest_bam_threads(&path, t).expect("digest");
        println!("threads={t} secs={:.1} count={} sha256={}", start.elapsed().as_secs_f64(), d.count, d.sha256);
    }
}
