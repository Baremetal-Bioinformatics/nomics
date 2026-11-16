//! sarek's `reports/` tree for `nomics run --reports`: the QC reports sarek 3.10.0 publishes
//! by default, at sarek's paths and under sarek's names, made from this run's own reads,
//! BAMs and VCF. Gate: `scripts/qc_run_gate.sh` (rules in `docs/GATES.md`).
//!
//! ```text
//! reports/fastp/<sample>/<sample>-<lane>.fastp.{json,html,log}
//! reports/fastqc/<sample>-<lane>/<sample>-<lane>_{1,2}_fastqc/{fastqc_data.txt,summary.txt}
//! reports/markduplicates/<sample>/<sample>.md.bam.metrics
//! reports/samtools/<sample>/<sample>.{md,recal}.bam.stats
//! reports/mosdepth/<sample>/<sample>.{md,recal}.*
//! reports/vcftools/haplotypecaller/<sample>/<sample>.haplotypecaller.{TsTv.count,TsTv.qual,FILTER.summary}
//! reports/bcftools/haplotypecaller/<sample>/<sample>.haplotypecaller.bcftools_stats.txt
//! ```
//!
//! FastQC's reports are the two text members of its zip, in the directory FastQC itself
//! writes with `--extract`; its images and HTML page are not reproduced. `<sample>` is the
//! run's `--prefix` (sarek's samplesheet `sample`).

use std::io;
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use nomics_prep::report::Report;

use crate::pipeline::Args;

fn other(message: impl Into<String>) -> io::Error {
    io::Error::other(message.into())
}

/// What `--reports` cannot do yet, said before any work is done.
pub fn check(args: &Args) -> io::Result<()> {
    if !args.reports {
        return Ok(());
    }
    if args.resume_from.is_some() {
        return Err(other(
            "--reports needs a whole run: the read reports come from the FASTQ pass, which --resume-from skips",
        ));
    }
    if !args.save_output_as_bam {
        return Err(other(
            "--reports needs --save-output-as-bam true: the samtools and mosdepth reports read BAM, not CRAM, for now",
        ));
    }
    Ok(())
}

/// The command line fastp records in its reports, as sarek's FASTP task runs it on a lane
/// staged as `<sample>-<lane>_R{1,2}.fastq.gz` (trailing space included). The thread count
/// is sarek's, whatever this run uses: fastp's insert-size sample depends on it.
pub fn fastp_command(p: &str, fastp_threads: usize) -> String {
    format!(
        "fastp --in1 {p}_R1.fastq.gz --in2 {p}_R2.fastq.gz --out1 {p}_R1.fastp.fastq.gz --out2 {p}_R2.fastp.fastq.gz \
         --json {p}.fastp.json --html {p}.fastp.html --thread {fastp_threads} --detect_adapter_for_pe \
         --disable_adapter_trimming --split_by_lines 200000000 --length_required 15 "
    )
}

/// One lane's fastp reports: the statistics are collected by prep's own pass.
pub struct Fastp {
    report: Option<Report>,
    started: u64,
}

fn epoch_seconds() -> u64 {
    SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default().as_secs()
}

impl Fastp {
    pub fn start(wanted: bool) -> io::Result<Fastp> {
        Ok(Fastp {
            report: if wanted { Some(Report::new()?) } else { None },
            started: epoch_seconds(),
        })
    }

    pub fn report(&mut self) -> Option<&mut Report> {
        self.report.as_mut()
    }

    /// Write `<p>.fastp.{json,html,log}` under `reports/fastp/<sample>/`.
    pub fn write(
        &self,
        outdir: &Path,
        sample: &str,
        p: &str,
        stats: &nomics_prep::Stats,
        fastp_threads: usize,
    ) -> io::Result<()> {
        let Some(report) = &self.report else { return Ok(()) };
        let dir = outdir.join("reports/fastp").join(sample);
        std::fs::create_dir_all(&dir)?;
        let command = fastp_command(p, fastp_threads);
        // fastp prints the report names as given on its command line: relative.
        let (json, html) = (format!("{p}.fastp.json"), format!("{p}.fastp.html"));
        std::fs::write(dir.join(&json), report.json(stats, &command))?;
        let now = nomics_prep::report::timestamp_now();
        std::fs::write(dir.join(&html), report.html(stats, &command, &now, &now))?;
        let seconds = epoch_seconds().saturating_sub(self.started);
        std::fs::write(
            dir.join(format!("{p}.fastp.log")),
            report.log(stats, &command, &json, &html, seconds),
        )
    }
}

/// FastQC running beside the pipeline, one thread per FASTQ.
pub struct Pending {
    fastqc: Vec<JoinHandle<io::Result<()>>>,
    started: Instant,
}

/// Start FastQC on every lane's two FASTQs. sarek stages them as `<sample>-<lane>_{1,2}.gz`,
/// the name FastQC prints in both files.
pub fn start_fastqc<'a>(outdir: &Path, lanes: impl Iterator<Item = (&'a str, &'a Path, &'a Path)>) -> Pending {
    let mut fastqc = Vec::new();
    for (p, fastq1, fastq2) in lanes {
        for (mate, fastq) in [(1, fastq1), (2, fastq2)] {
            let dir = outdir.join("reports/fastqc").join(p).join(format!("{p}_{mate}_fastqc"));
            let (fastq, name) = (fastq.to_owned(), format!("{p}_{mate}.gz"));
            fastqc.push(std::thread::spawn(move || {
                std::fs::create_dir_all(&dir)?;
                nomics_qc::fastqc::fastqc(&fastq, &dir, Some(&name))
            }));
        }
    }
    Pending { fastqc, started: Instant::now() }
}

fn timed(what: &str, run: impl FnOnce() -> io::Result<()>) -> io::Result<()> {
    let started = Instant::now();
    run().map_err(|e| other(format!("reports: {what}: {e}")))?;
    eprintln!("reports: {what} in {:.1}s", started.elapsed().as_secs_f64());
    Ok(())
}

/// The reports made from the run's published files, then FastQC's threads joined.
/// `vcf` is the HaplotypeCaller VCF, or `None` when sarek would publish no VCF reports.
pub fn finish(pending: Pending, args: &Args, fasta: &Path, vcf: Option<&Path>, threads: usize) -> io::Result<()> {
    let (outdir, sample) = (&args.outdir, args.prefix.as_str());
    let reports = outdir.join("reports");
    let dir = |tool: &str| -> io::Result<PathBuf> {
        let d = reports.join(tool).join(sample);
        std::fs::create_dir_all(&d)?;
        Ok(d)
    };
    // MarkDuplicates' metrics: sarek publishes the file here too, under the BAM's name.
    std::fs::copy(
        outdir.join(format!("{sample}.md.metrics")),
        dir("markduplicates")?.join(format!("{sample}.md.bam.metrics")),
    )?;

    let (samtools, mosdepth) = (dir("samtools")?, dir("mosdepth")?);
    // sarek: exomes `mosdepth --by <targets>` with a per-base file; genomes
    // `mosdepth -n --fast-mode --by 500`.
    let by = if args.wes { args.targets.to_string_lossy().into_owned() } else { "500".to_owned() };
    let mosdepth_threads = (threads / 4).clamp(1, 4);
    std::thread::scope(|scope| -> io::Result<()> {
        let mut jobs = Vec::new();
        for kind in ["md", "recal"] {
            let bam = outdir.join(format!("{sample}.{kind}.bam"));
            let stats = samtools.join(format!("{sample}.{kind}.bam.stats"));
            let (bam2, by, prefix) = (bam.clone(), by.as_str(), mosdepth.join(format!("{sample}.{kind}")));
            jobs.push(scope.spawn(move || {
                timed(&format!("samtools stats, {kind}.bam"), || {
                    nomics_qc::samtools_stats::bam_stats(&bam, fasta, &stats, None, None)
                })
            }));
            jobs.push(scope.spawn(move || {
                timed(&format!("mosdepth, {kind}.bam"), || {
                    nomics_qc::mosdepth::run(&bam2, &prefix, by, !args.wes, !args.wes, mosdepth_threads)
                })
            }));
        }
        if let Some(vcf) = vcf {
            let reports = &reports;
            jobs.push(scope.spawn(move || {
                timed("vcftools and bcftools statistics", || {
                    let vt = reports.join("vcftools/haplotypecaller").join(sample);
                    let bt = reports.join("bcftools/haplotypecaller").join(sample);
                    std::fs::create_dir_all(&bt)?;
                    let arch = args.reference_arch.unwrap_or_else(nomics_qc::ReferenceArch::host);
                    nomics_qc::vcf_stats(vcf, &vt, None, arch)?;
                    // One directory per tool in sarek's tree.
                    let name = format!("{sample}.haplotypecaller.bcftools_stats.txt");
                    std::fs::rename(vt.join(&name), bt.join(&name))
                })
            }));
        }
        jobs.into_iter()
            .try_for_each(|j| j.join().map_err(|_| other("a reports thread panicked"))?)
    })?;

    for handle in pending.fastqc {
        handle
            .join()
            .map_err(|_| other("a FastQC thread panicked"))?
            .map_err(|e| other(format!("reports: FastQC: {e}")))?;
    }
    eprintln!(
        "reports: FastQC finished within {:.1}s of the run's start (it runs beside the pipeline)",
        pending.started.elapsed().as_secs_f64()
    );
    Ok(())
}
