//! Fault injection without sarek: a synthetic sarek tree and nomics run
//! directory in the real layouts, gated by the shipped descriptor. Each fault
//! must give FAIL or ERROR, name the artifact, and never PASS.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use nomics_io::bam::{BinRule, Header, IntType, Record, Reference, Value, write_raw_record};
use nomics_io::bgzf::{BgzfWriter, Profile};
use nomics_io::deflate::HtsjdkDeflater;

use super::descriptor::{Descriptor, SAREK_GERMLINE};
use super::gate::{self, Roots, Status};
use super::report::{self, Verdict};

const P: &str = "SAMPLE";

fn bgzf(path: &Path, bytes: &[u8]) {
    let mut w = BgzfWriter::new(File::create(path).unwrap(), HtsjdkDeflater::new(), Profile::HTSJDK);
    w.write_all(bytes).unwrap();
    w.finish().unwrap().flush().unwrap();
}

fn rec(i: i32, mapq: u8) -> Record {
    Record {
        ref_id: 0,
        pos: 100 + i * 10,
        mapq,
        flag: 99,
        next_ref_id: 0,
        next_pos: 300 + i * 10,
        tlen: 250,
        read_name: format!("r{i}").into_bytes(),
        cigar: vec![(4 << 4)],
        seq: b"ACGT".to_vec(),
        qual: vec![30; 4],
        tags: vec![(*b"NM", Value::Int(IntType::U8, 0))],
    }
}

fn bam(path: &Path, pg: &str, recs: &[Record]) {
    let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:{P}-L1\tSM:VERIFY_{P}\tDS:/ref/x.fasta\n@PG\tID:{pg}\n");
    let h = Header {
        text: text.into_bytes(),
        refs: vec![Reference { name: b"chr1".to_vec(), len: 100_000 }],
    };
    let mut w = BgzfWriter::new(File::create(path).unwrap(), HtsjdkDeflater::new(), Profile::HTSJDK);
    h.write(&mut w).unwrap();
    for r in recs {
        let mut raw = Vec::new();
        r.encode(BinRule::Htsjdk, &mut raw).unwrap();
        write_raw_record(&mut w, &raw).unwrap();
    }
    w.finish().unwrap().flush().unwrap();
}

fn vcf_text(command: &str, records: &[String]) -> String {
    let mut t = format!("##fileformat=VCFv4.2\n##GATKCommandLine={command}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVERIFY_{P}\n");
    for r in records {
        t.push_str(r);
        t.push('\n');
    }
    t
}

fn calls() -> Vec<String> {
    (0..5).map(|i| format!("chr1\t{}\t.\tA\tG\t{}.64\t.\tDP=9\tGT:PL\t0/1:10,0,20", 1000 + i, 800 + i)).collect()
}

fn metrics(started: &str) -> String {
    format!("## htsjdk.samtools.metrics.StringHeader\n# MarkDuplicates --INPUT x\n## htsjdk.samtools.metrics.StringHeader\n# Started on: {started}\nLIBRARY\tREADS\nlib\t10\n")
}

fn fastq(n: usize, tweak: Option<usize>) -> String {
    (0..n)
        .map(|i| format!("@r{i}\n{}\n+\nIIII\n", if tweak == Some(i) { "ACGA" } else { "ACGT" }))
        .collect()
}

struct Fixture {
    base: PathBuf,
    work: PathBuf,
    results: PathBuf,
    run: PathBuf,
}

impl Fixture {
    /// sarek's layout (published files plus fastp chunks and per-chunk BAMs in
    /// task directories) and nomics' run directory, with identical content.
    fn new(tag: &str) -> Fixture {
        let base = std::env::temp_dir().join(format!("nomics-verify-fault-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        let f = Fixture {
            work: base.join("sarek/work"),
            results: base.join("sarek/results"),
            run: base.join("nomics"),
            base,
        };
        let fastp = f.work.join("aa/bbbb");
        let bwa = f.work.join("cc/dddd");
        let r = &f.results;
        for d in [
            fastp.clone(),
            bwa.clone(),
            r.join(format!("preprocessing/markduplicates/{P}")),
            r.join(format!("preprocessing/recalibrated/{P}")),
            r.join(format!("preprocessing/recal_table/{P}")),
            r.join(format!("reports/markduplicates/{P}")),
            r.join(format!("variant_calling/haplotypecaller/{P}")),
            f.run.join("prep"),
            f.run.join("mapped"),
        ] {
            std::fs::create_dir_all(d).unwrap();
        }
        let recs: Vec<Record> = (0..6).map(|i| rec(i, 60)).collect();
        for (dir, pg) in [(&fastp, "sarek"), (&f.run.join("prep"), "nomics")] {
            let _ = pg;
            for mate in [1, 2] {
                bgzf(&dir.join(format!("0001.{P}-L1_R{mate}.fastp.fastq.gz")), fastq(3, None).as_bytes());
            }
        }
        bam(&bwa.join(format!("{P}.0001.bam")), "bwa", &recs);
        bam(&f.run.join(format!("mapped/{P}.0001.bam")), "nomics", &recs);
        bam(&r.join(format!("preprocessing/markduplicates/{P}/{P}.md.bam")), "MarkDuplicates", &recs);
        bam(&f.run.join(format!("{P}.md.bam")), "nomics", &recs);
        bam(&r.join(format!("preprocessing/recalibrated/{P}/{P}.recal.bam")), "ApplyBQSR", &recs);
        bam(&f.run.join(format!("{P}.recal.bam")), "nomics", &recs);
        std::fs::write(r.join(format!("reports/markduplicates/{P}/{P}.md.bam.metrics")), metrics("Mon")).unwrap();
        std::fs::write(f.run.join(format!("{P}.md.metrics")), metrics("Tue")).unwrap();
        let table = "#:GATKReport.v1.1:5\n#:GATKTable:2:1:%s:%s:;\nA B\n1 2\n";
        std::fs::write(r.join(format!("preprocessing/recal_table/{P}/{P}.recal.table")), table).unwrap();
        std::fs::write(f.run.join(format!("{P}.recal.table")), table).unwrap();
        bgzf(&f.reference_vcf(), vcf_text("sarek /work/a", &calls()).as_bytes());
        bgzf(&f.nomics_vcf(), vcf_text("nomics", &calls()).as_bytes());
        f
    }

    fn reference_vcf(&self) -> PathBuf {
        self.results.join(format!("variant_calling/haplotypecaller/{P}/{P}.haplotypecaller.vcf.gz"))
    }

    fn nomics_vcf(&self) -> PathBuf {
        self.run.join(format!("{P}.haplotypecaller.vcf.gz"))
    }

    fn gate(&self) -> (Vec<gate::ArtifactResult>, Verdict, String) {
        let d = Descriptor::parse(SAREK_GERMLINE).unwrap();
        let reference = Roots { work: Some(self.work.clone()), results: Some(self.results.clone()), run: None };
        let nomics = Roots { run: Some(self.run.clone()), ..Default::default() };
        let results = gate::gate_all(&d.artifacts, &reference, &nomics, P, "L1").unwrap();
        let (v, why) = report::verdict(&results, 1);
        (results, v, why)
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.base);
    }
}

fn status<'a>(results: &'a [gate::ArtifactResult], name: &str) -> &'a gate::ArtifactResult {
    results.iter().find(|r| r.name == name).unwrap_or_else(|| panic!("no {name}"))
}

#[test]
fn identical_outputs_pass_under_the_gates_normalizations() {
    // The fixture's sides differ only where GATES.md allows: @PG lines,
    // ##GATKCommandLine, and metrics lines 2 and 4.
    let f = Fixture::new("pass");
    let (results, v, why) = f.gate();
    assert_eq!(v, Verdict::Pass, "{why}");
    let names: Vec<&str> = results.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(
        names,
        [
            "chunk/0001.SAMPLE-L1_R1.fastp.fastq.gz",
            "chunk/0001.SAMPLE-L1_R2.fastp.fastq.gz",
            "mapped/SAMPLE.0001.bam",
            "md_bam",
            "md_metrics",
            "recal_table",
            "recal_bam",
            "vcf"
        ]
    );
    assert!(results.iter().all(|r| r.status == Status::Pass));
}

#[test]
fn a_one_record_vcf_change_fails_and_names_the_record() {
    let f = Fixture::new("vcf");
    let mut c = calls();
    c[2] = c[2].replace("802.64", "802.65");
    bgzf(&f.nomics_vcf(), vcf_text("nomics", &c).as_bytes());
    let (results, v, why) = f.gate();
    assert_eq!(v, Verdict::Fail);
    assert!(why.contains("vcf"), "{why}");
    let r = status(&results, "vcf");
    assert_eq!(r.status, Status::Fail);
    let d = r.first_difference.as_ref().unwrap();
    assert_eq!(d.index, Some(3));
    assert!(d.detail.contains("chr1:1002: QUAL 802.64 (reference) vs 802.65 (nomics)"), "{}", d.detail);
}

#[test]
fn a_one_record_bam_change_fails_and_names_the_field() {
    let f = Fixture::new("bam");
    let mut recs: Vec<Record> = (0..6).map(|i| rec(i, 60)).collect();
    recs[4] = rec(4, 59);
    bam(&f.run.join(format!("{P}.md.bam")), "nomics", &recs);
    let (results, v, _) = f.gate();
    assert_eq!(v, Verdict::Fail);
    let r = status(&results, "md_bam");
    let d = r.first_difference.as_ref().unwrap();
    assert_eq!(d.index, Some(5));
    assert!(d.detail.contains("mapq"), "{}", d.detail);
    // Everything else still passes; only md_bam is named.
    assert_eq!(results.iter().filter(|r| r.status != Status::Pass).count(), 1);
}

#[test]
fn a_changed_read_in_a_chunk_fails_at_that_read() {
    let f = Fixture::new("chunk");
    bgzf(&f.run.join(format!("prep/0001.{P}-L1_R2.fastp.fastq.gz")), fastq(3, Some(1)).as_bytes());
    let (results, v, _) = f.gate();
    assert_eq!(v, Verdict::Fail);
    let d = status(&results, "chunk/0001.SAMPLE-L1_R2.fastp.fastq.gz").first_difference.clone().unwrap();
    assert_eq!(d.index, Some(2));
}

#[test]
fn empty_outputs_on_both_sides_are_an_error_never_a_pass() {
    let f = Fixture::new("empty");
    bgzf(&f.reference_vcf(), vcf_text("sarek", &[]).as_bytes());
    bgzf(&f.nomics_vcf(), vcf_text("nomics", &[]).as_bytes());
    let (results, v, why) = f.gate();
    assert_eq!(v, Verdict::Error, "{why}");
    assert!(why.contains("empty output") && why.contains("vcf"), "{why}");
    assert_eq!(status(&results, "vcf").status, Status::Empty);
}

#[test]
fn a_truncated_bgzf_output_fails_and_names_the_artifact() {
    let f = Fixture::new("truncated");
    let bytes = std::fs::read(f.nomics_vcf()).unwrap();
    std::fs::write(f.nomics_vcf(), &bytes[..bytes.len() / 2]).unwrap();
    let (results, v, why) = f.gate();
    assert_eq!(v, Verdict::Fail, "{why}");
    let r = status(&results, "vcf");
    assert_eq!(r.status, Status::Unreadable);
    assert!(r.first_difference.as_ref().unwrap().detail.starts_with("nomics file unreadable"));
}

#[test]
fn a_missing_output_fails_and_a_missing_reference_is_an_error() {
    let f = Fixture::new("missing");
    std::fs::remove_file(f.run.join(format!("{P}.recal.table"))).unwrap();
    let (results, v, _) = f.gate();
    assert_eq!(v, Verdict::Fail);
    assert_eq!(status(&results, "recal_table").status, Status::Missing);

    let g = Fixture::new("noref");
    std::fs::remove_file(g.reference_vcf()).unwrap();
    let (_, v, why) = g.gate();
    assert_eq!(v, Verdict::Error, "{why}");
    assert!(why.contains("vcf"));
}

#[test]
fn metrics_differ_outside_lines_2_and_4_fails() {
    let f = Fixture::new("metrics");
    std::fs::write(f.run.join(format!("{P}.md.metrics")), metrics("Tue").replace("lib\t10", "lib\t11")).unwrap();
    let (results, v, _) = f.gate();
    assert_eq!(v, Verdict::Fail);
    let d = status(&results, "md_metrics").first_difference.clone().unwrap();
    assert_eq!(d.index, Some(6));
}
