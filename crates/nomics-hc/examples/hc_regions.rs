//! Stage 2-3 gate: read stream (filters, downsampler) and assembly regions.
//!
//! usage: hc_regions <recal.bam> <intervals.bed> <reference.fasta> <oracle regions.igv> <golden.vcf.gz>
//!
//! The golden VCF supplies the QD-jitter records: each active region draws one
//! nextGaussian per jittered record it contains, at GATK's call time, and the
//! resulting values are compared with the golden QDs.
//! Prints GATK-log-style filter counts and compares regions with the oracle's
//! `--assembly-region-out`.
use std::path::Path;
use std::time::Instant;

use nomics_hc::activity::traverse;
use nomics_hc::javafmt::format_fixed;
use nomics_hc::jrandom::JavaRandom;
use nomics_hc::reads::{FILTERS, Intervals, read_stream};
use nomics_io::bam::RecordSet;
use nomics_io::fasta::Fasta;
use nomics_io::gzip::GzReader;
use std::io::Read as _;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let t = Instant::now();
    let bam = RecordSet::load(Path::new(&a[1]), 20).unwrap();
    let contigs: Vec<(Vec<u8>, i64)> = bam
        .header
        .refs
        .iter()
        .map(|r| (r.name.clone(), r.len as i64))
        .collect();
    let lens: Vec<i64> = contigs.iter().map(|c| c.1).collect();
    let intervals = Intervals::from_bed(Path::new(&a[2]), &contigs).unwrap();
    let fasta = Fasta::load(Path::new(&a[3]), 20).unwrap();
    for (i, (n, _)) in contigs.iter().enumerate() {
        assert_eq!(
            fasta.index[i].name.as_bytes(),
            n.as_slice(),
            "FASTA and BAM contig order differ"
        );
    }
    eprintln!("loaded in {:.1}s", t.elapsed().as_secs_f64());
    let t = Instant::now();
    let stream = read_stream(&bam, &intervals, &lens).unwrap();
    let total: u64 = stream.filtered.iter().sum();
    for (f, n) in FILTERS.iter().zip(stream.filtered) {
        println!("{n} read(s) filtered by: {}", f.name());
    }
    println!(
        "{total} total reads filtered out of {} reads processed ({:.1}s)",
        stream.processed,
        t.elapsed().as_secs_f64()
    );
    // Golden QD-jitter records: (contig, pos, printed QD).
    let names: std::collections::HashMap<Vec<u8>, usize> = contigs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.clone(), i))
        .collect();
    let mut vcf = String::new();
    GzReader::open(Path::new(&a[5]))
        .unwrap()
        .read_to_string(&mut vcf)
        .unwrap();
    let mut jitter: Vec<(usize, i64, String)> = Vec::new();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let qual: f64 = f[5].parse().unwrap();
        let qd = f[7]
            .split(';')
            .find_map(|kv| kv.strip_prefix("QD="))
            .map(str::to_string);
        let g: Vec<&str> = f[9].split(':').collect();
        let ad: f64 = g[1].split(',').map(|x| x.parse::<f64>().unwrap()).sum();
        let Some(qd) = qd else { continue };
        if ad == 0.0 {
            continue;
        }
        let raw = qual / ad;
        let jittered = if (raw - 35.0).abs() < 0.01 {
            format_fixed(raw, 2) != qd
        } else {
            raw >= 35.0
        };
        if jittered {
            jitter.push((names[f[0].as_bytes()], f[1].parse().unwrap(), qd));
        }
    }
    println!("{} jittered QD records in the golden VCF", jitter.len());
    let mut j = 0usize;
    let (mut qd_ok, mut qd_bad, mut first_bad) = (0usize, 0usize, None);
    let mut rng = JavaRandom::gatk();
    let t = Instant::now();
    let reference = |c: usize, pos: i64| {
        let b = fasta.contig(c)[(pos - 1) as usize].to_ascii_uppercase();
        if matches!(b, b'A' | b'C' | b'G' | b'T') {
            b
        } else {
            b'N'
        }
    };
    let mut ours_state: Vec<String> = Vec::new();
    // NOMICS_DUMP_REGIONS=<path>: active regions with their reads (name:flag), for the GATK finalize oracle.
    let mut dump = std::env::var_os("NOMICS_DUMP_REGIONS")
        .map(|p| std::io::BufWriter::new(std::fs::File::create(p).unwrap()));
    let mut on_call = |r: &nomics_hc::activity::Region, reads: &[u32], rng: &mut JavaRandom| {
        if !r.active {
            return;
        }
        if let Some(w) = dump.as_mut()
            && !reads.is_empty()
        {
            use std::io::Write as _;
            write!(
                w,
                "{}\t{}\t{}",
                String::from_utf8_lossy(&contigs[r.contig].0),
                r.start,
                r.end
            )
            .unwrap();
            for &i in reads {
                let v = nomics_io::bam::RecordView::new(bam.record(i as usize)).unwrap();
                write!(
                    w,
                    "\t{}:{}",
                    String::from_utf8_lossy(v.read_name()),
                    v.flag()
                )
                .unwrap();
            }
            writeln!(w).unwrap();
        }
        if !reads.is_empty() {
            // `--debug-assembly-region-state` block: span, then name + flags per read.
            let mut b = format!(
                "{}:{}-{}\n",
                String::from_utf8_lossy(&contigs[r.contig].0),
                r.start,
                r.end
            );
            for &i in reads {
                let v = nomics_io::bam::RecordView::new(bam.record(i as usize)).unwrap();
                b.push_str(&format!(
                    "{}   {}\n",
                    String::from_utf8_lossy(v.read_name()),
                    v.flag()
                ));
            }
            ours_state.push(b);
        }
        // Skip golden records before this region (inactive regions hold none).
        while j < jitter.len() && (jitter[j].0, jitter[j].1) < (r.contig, r.start) {
            j += 1;
        }
        while j < jitter.len() && jitter[j].0 == r.contig && jitter[j].1 <= r.end {
            let g = rng.next_gaussian();
            let ours = format_fixed(30.0 + g * 3.0, 2);
            if ours == jitter[j].2 {
                qd_ok += 1;
            } else {
                qd_bad += 1;
                if first_bad.is_none() {
                    first_bad = Some((jitter[j].0, jitter[j].1, ours, jitter[j].2.clone()));
                }
            }
            j += 1;
        }
    };
    let tr = traverse(
        &bam,
        &stream,
        &intervals.unpadded,
        &reference,
        &lens,
        &mut rng,
        &mut on_call,
        1,
    )
    .unwrap();
    if let Some(path) = a.get(6) {
        // Oracle blocks: "<span>\nNumber of reads in region: N     they are:\n<name>   <flags>\n..."
        let text = std::fs::read_to_string(path).unwrap();
        let mut oracle_state = Vec::new();
        let lines: Vec<&str> = text.lines().collect();
        let mut i = 0;
        while i < lines.len() {
            if let Some(rest) = lines
                .get(i + 1)
                .and_then(|l| l.strip_prefix("Number of reads in region: "))
            {
                let n: usize = rest.split_whitespace().next().unwrap().parse().unwrap();
                let mut b = format!("{}\n", lines[i]);
                for l in &lines[i + 2..i + 2 + n] {
                    b.push_str(l);
                    b.push('\n');
                }
                oracle_state.push(b);
                i += 2 + n;
            } else {
                i += 1;
            }
        }
        let first = ours_state
            .iter()
            .zip(&oracle_state)
            .position(|(x, y)| x != y);
        println!(
            "region membership: ours {} oracle {} -> {}",
            ours_state.len(),
            oracle_state.len(),
            if ours_state == oracle_state {
                "IDENTICAL".to_string()
            } else {
                format!("DIFFER at {first:?}")
            }
        );
        if let Some(k) = first {
            println!("ours:\n{}\noracle:\n{}", ours_state[k], oracle_state[k]);
        }
    }
    let ours = tr.regions;
    println!("downsampler: {} nextInt draws", tr.draws);
    println!("QD jitter: {qd_ok} match, {qd_bad} differ; first difference {first_bad:?}");
    println!(
        "{} regions ({} active) in {:.1}s",
        ours.len(),
        ours.iter().filter(|r| r.active).count(),
        t.elapsed().as_secs_f64()
    );
    // Oracle: "chr  start0  end  size=N  ±1.00000" lines.
    let names: std::collections::HashMap<Vec<u8>, usize> = contigs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.clone(), i))
        .collect();
    let mut oracle = Vec::new();
    for line in std::fs::read_to_string(&a[4]).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() == 5 && f[3].starts_with("size=") {
            let start: i64 = f[1].parse::<i64>().unwrap() + 1;
            let end: i64 = f[2].parse().unwrap();
            oracle.push((names[f[0].as_bytes()], start, end, f[4].starts_with('1')));
        }
    }
    let ours_t: Vec<_> = ours
        .iter()
        .map(|r| (r.contig, r.start, r.end, r.active))
        .collect();
    let first = ours_t.iter().zip(&oracle).position(|(x, y)| x != y);
    println!(
        "regions: ours {} oracle {} -> {}",
        ours_t.len(),
        oracle.len(),
        if ours_t == oracle {
            "IDENTICAL".to_string()
        } else {
            format!("DIFFER at {first:?}")
        }
    );
    if let Some(i) = first {
        for j in i.saturating_sub(2)..(i + 3).min(ours_t.len()).min(oracle.len()) {
            println!("  ours {:?}  oracle {:?}", ours_t[j], oracle[j]);
        }
    }
}
