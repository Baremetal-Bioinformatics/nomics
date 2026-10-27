//! End-to-end HaplotypeCaller gate: BAM → VCF through every stage, compared
//! record-for-record with the golden VCF (header compared without
//! `##GATKCommandLine`, per `docs/GATES.md`).
//!
//! usage: hc_call <recal.bam> <intervals.bed> <reference.fasta> <dbsnp.vcf.gz> <out.vcf.gz> [golden.vcf.gz]
use std::io::Read as _;
use std::path::Path;
use std::time::Instant;

use nomics_hc::caller::{Context, header_meta, run};
use nomics_hc::dbsnp::Dbsnp;
use nomics_hc::reads::{Intervals, read_stream};
use nomics_hc::vcf::{Header, encode, write_vcf_gz};
use nomics_io::bam::RecordSet;
use nomics_io::fasta::Fasta;
use nomics_io::gzip::GzReader;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let t0 = Instant::now();
    // NOMICS_STREAM_GROUP=1: stream the BAM and keep only the records the
    // group's overlap query returns (whole-genome BAMs do not fit in memory).
    let bam = if std::env::var_os("NOMICS_STREAM_GROUP").is_some() {
        let stream =
            nomics_io::bamstream::BamStream::open(Path::new(&a[1]), 1 << 14, 4, 8).unwrap();
        let contigs: Vec<(Vec<u8>, i64)> = stream
            .header
            .refs
            .iter()
            .map(|r| (r.name.clone(), r.len as i64))
            .collect();
        let groups = [Intervals::from_bed(Path::new(&a[2]), &contigs).unwrap()];
        nomics_hc::reads::GroupCollector::new(stream, &groups)
            .take(0)
            .unwrap()
    } else {
        RecordSet::load(Path::new(&a[1]), 20).unwrap()
    };
    let contigs: Vec<(Vec<u8>, i64)> = bam
        .header
        .refs
        .iter()
        .map(|r| (r.name.clone(), r.len as i64))
        .collect();
    let lens: Vec<i64> = contigs.iter().map(|c| c.1).collect();
    let names: Vec<String> = contigs
        .iter()
        .map(|c| String::from_utf8_lossy(&c.0).into_owned())
        .collect();
    let sample = String::from_utf8_lossy(&bam.header.text)
        .lines()
        .find(|l| l.starts_with("@RG"))
        .and_then(|l| l.split('\t').find_map(|f| f.strip_prefix("SM:")))
        .expect("@RG SM")
        .to_string();
    let intervals = Intervals::from_bed(Path::new(&a[2]), &contigs).unwrap();
    let fasta = Fasta::load(Path::new(&a[3]), 20).unwrap();
    for (i, n) in names.iter().enumerate() {
        assert_eq!(&fasta.index[i].name, n, "FASTA and BAM contig order differ");
    }
    let dbsnp = Dbsnp::open(Path::new(&a[4])).unwrap();
    eprintln!("loaded in {:.1}s", t0.elapsed().as_secs_f64());

    let t = Instant::now();
    let stream = read_stream(&bam, &intervals, &lens).unwrap();
    let reference = |c: usize, pos: i64| {
        let b = fasta.contig(c)[(pos - 1) as usize].to_ascii_uppercase();
        if matches!(b, b'A' | b'C' | b'G' | b'T') {
            b
        } else {
            b'N'
        }
    };
    let contig_bases = |c: usize| fasta.contig(c);
    let ctx = Context {
        contig_names: names.clone(),
        contig_lens: lens.clone(),
        reference: &contig_bases,
        dbsnp: Some(&dbsnp),
    };
    let threads: usize = std::env::var("NOMICS_THREADS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(20);
    let gvcf = std::env::var_os("NOMICS_GVCF").is_some();
    let run = if gvcf {
        nomics_hc::caller::run_gvcf
    } else {
        run
    };
    let header_meta = if gvcf {
        nomics_hc::gvcf::header_meta
    } else {
        header_meta
    };
    let (tr, records, timing) = run(
        &ctx,
        &bam,
        &stream,
        &intervals.unpadded,
        &reference,
        threads,
    )
    .unwrap();
    let names_t = [
        "finalize",
        "assemble",
        "trim",
        "likelihoods",
        "realign",
        "call",
        "annotate",
    ];
    for (n, d) in names_t.iter().zip(timing) {
        eprintln!("  {n:<12} {:8.1}s (summed over threads)", d.as_secs_f64());
    }
    let big = stream
        .groups
        .iter()
        .flatten()
        .filter(|g| g.reads.len() > nomics_hc::reads::MAX_READS_PER_START)
        .count();
    eprintln!(
        "{} regions, {} records, {} downsampled start positions; {:.1}s on {} threads",
        tr.regions.len(),
        records.len(),
        big,
        t.elapsed().as_secs_f64(),
        threads
    );

    let header = Header {
        meta: header_meta(
            r#"GATKCommandLine=<ID=HaplotypeCaller,CommandLine="nomics hc",Version="4.6.2.0",Date="">"#,
        ),
        contigs: names
            .iter()
            .cloned()
            .zip(lens.iter().map(|&l| l as u64))
            .collect(),
        sample,
    };
    let out = std::fs::File::create(&a[5]).unwrap();
    write_vcf_gz(std::io::BufWriter::new(out), &header, &records, 8).unwrap();

    let Some(golden) = a.get(6) else { return };
    let mut text = String::new();
    GzReader::open(Path::new(golden))
        .unwrap()
        .read_to_string(&mut text)
        .unwrap();
    let ours_header = header.text();
    let strip = |h: &str| -> Vec<String> {
        h.lines()
            .filter(|l| l.starts_with('#') && !l.starts_with("##GATKCommandLine"))
            .map(str::to_string)
            .collect()
    };
    let header_ok = strip(&ours_header) == strip(&text);
    let golden_recs: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    let mut buf = Vec::new();
    for r in &records {
        encode(r, &header.contigs, &mut buf);
    }
    let ours = String::from_utf8(buf).unwrap();
    let ours_recs: Vec<&str> = ours.lines().collect();
    let (mut same, mut shown) = (0, 0);
    for (o, g) in ours_recs.iter().zip(&golden_recs) {
        if o == g {
            same += 1;
        } else if shown < 10 {
            shown += 1;
            eprintln!("ours:   {o}\ngolden: {g}");
        }
    }
    println!(
        "header (excl. GATKCommandLine): {}; records: {} ours / {} golden, {} identical in order",
        if header_ok { "IDENTICAL" } else { "DIFFERS" },
        ours_recs.len(),
        golden_recs.len(),
        same
    );
}
