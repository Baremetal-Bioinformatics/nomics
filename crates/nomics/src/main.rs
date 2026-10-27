use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use clap::{Args as ClapArgs, Parser, Subcommand};

mod pipeline;
mod scatter;
mod shim;
mod tier0;
use nomics_mem::HostBuf;

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// GATK 4.6.2.0 GenotypeGVCFs on one sample's gVCF (sarek --joint_germline).
    JointGenotype {
        #[arg(long)]
        gvcf: PathBuf,
        #[arg(long)]
        dbsnp: Option<PathBuf>,
        /// Output VCF (`.gz` for BGZF).
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value_t = 6)]
        threads: usize,
    },
    /// STAR 2.7.10b counting from an oracle BAM with raw barcode/UMI and gx tags.
    SoloCount {
        #[arg(long)]
        bam: PathBuf,
        #[arg(long)]
        whitelist: PathBuf,
        /// STAR features.tsv in original index order (GTF derivation deferred).
        #[arg(long)]
        features: PathBuf,
        /// New output directory, usually Solo.out.
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value_t = 6)]
        threads: usize,
    },
    /// Quantify a transcriptome BAM (Salmon 2.7.0 deterministic alignment mode).
    SalmonQuant {
        #[arg(long = "geneMap")]
        gene_map: Option<PathBuf>,
        #[arg(long, default_value_t = 6)]
        threads: usize,
        #[arg(short = 't', long = "targets")]
        transcripts: PathBuf,
        #[arg(short = 'a', long = "alignments")]
        alignments: PathBuf,
        #[arg(short = 'l', long = "libType", default_value = "ISR")]
        lib_type: String,
        /// Prior weight for incompatible placements (<= 0 drops them).
        #[arg(
            long = "incompatPrior",
            default_value_t = 0.0,
            allow_hyphen_values = true
        )]
        incompat_prior: f64,
        #[arg(short = 'o', long = "output")]
        output: PathBuf,
    },
    /// Check the unified-memory premise on this machine: huge-page host buffers,
    /// and the GPU reading them in place.
    Selftest {
        /// Buffer size in MiB.
        #[arg(long, default_value_t = 1024)]
        mib: usize,
    },
    /// FASTQ pairs to VCF in one process.
    Run {
        /// R1 FASTQ; repeat once per lane, in samplesheet order.
        #[arg(long, required = true)]
        fastq1: Vec<PathBuf>,
        /// R2 FASTQ; repeat once per lane, in the same order as --fastq1.
        #[arg(long, required = true)]
        fastq2: Vec<PathBuf>,
        /// Lane names as in the samplesheet, one per FASTQ pair (default L1, L2, ...).
        #[arg(long = "lane")]
        lanes: Vec<String>,
        /// Reference directory laid out like iGenomes GATK.GRCh38.
        #[arg(long)]
        reference: PathBuf,
        /// Target intervals (BED), as passed to sarek with `--wes --intervals`.
        #[arg(long)]
        targets: PathBuf,
        #[arg(long)]
        outdir: PathBuf,
        /// dbSNP VCF; derived from the reference directory when omitted.
        #[arg(long)]
        dbsnp: Option<PathBuf>,
        /// BQSR known-sites VCFs; the iGenomes bundle defaults are used when omitted.
        #[arg(long = "known-sites")]
        known_sites: Vec<PathBuf>,
        #[arg(long)]
        known_sites_cache: Option<PathBuf>,
        #[arg(long)]
        dense_sa_cache: PathBuf,
        /// Dense-SA sampling interval (1, 2 or 4); the cache file must match.
        #[arg(long, default_value_t = 4)]
        dense_sa_interval: u64,
        /// The reference run's GatherBQSRReports and samtools-merge group
        /// orders: `gather <group>` / `merge <group>` lines (default: BED order).
        #[arg(long)]
        scatter_order: Option<PathBuf>,
        /// Read-group DS (the reference run's FASTA path); default: our FASTA path.
        #[arg(long)]
        read_group_ds: Option<String>,
        /// Sample name stored in the read group and VCF.
        #[arg(long, default_value = "HG002_HG002_WES")]
        sample: String,
        /// Published file prefix and read-group library name.
        #[arg(long, default_value = "HG002_WES")]
        prefix: String,
        #[arg(long, default_value_t = 20)]
        threads: usize,
        #[arg(long)]
        write_intermediates: bool,
        /// Write prep's chunks to `<outdir>/spill` and align from disk instead of
        /// memory (automatic when the gzipped FASTQs total over 16 GiB).
        #[arg(long)]
        spill_prep: bool,
        /// Publish the per-chunk sorted BAMs in <outdir>/mapped, as sarek's
        /// `--save_mapped` does. With `false`, whole genomes keep markdup's
        /// chunk spill in a lighter internal format, deleted once read.
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        save_mapped: bool,
        /// sarek's `--save_output_as_bam`. `false` (sarek's default) writes
        /// <outdir>/preprocessing/markduplicates/<prefix>/<prefix>.md.cram and
        /// <outdir>/preprocessing/recalibrated/<prefix>/<prefix>.recal.cram
        /// (+ .crai) decoding as sarek's do, instead of md.bam / recal.bam.
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        save_output_as_bam: bool,
        /// With --save-output-as-bam false on whole genomes (streamed
        /// stages), stream md.bam and recal.bam between the stages and keep
        /// them for --resume-from. Without it, the stages stream internal
        /// <prefix>.md.nspill / <prefix>.recal.nspill (a lighter format),
        /// each removed once the CRAM made from it is complete.
        #[arg(long)]
        keep_bams: bool,
        /// Only run HaplotypeCaller on `<outdir>/<prefix>.recal.bam` from an
        /// earlier run, streamed one interval group at a time.
        /// Restart an earlier run in --outdir at this stage from its own outputs
        /// (mapped chunk BAMs, md.bam or recal.bam, or their spill files;
        /// the newer if both are there), streaming them.
        #[arg(long, value_enum)]
        resume_from: Option<pipeline::Resume>,
        /// Same as `--resume-from hc`.
        #[arg(long)]
        hc_only: bool,
        /// Emit GATK reference-confidence blocks for subsequent joint genotyping.
        #[arg(long)]
        emit_gvcf: bool,
    },
    /// HaplotypeCaller on a recalibrated BAM, with one RNG reset per interval group.
    Hc {
        #[arg(long)]
        input: PathBuf,
        /// BED for one HC invocation; repeat for already scattered groups.
        #[arg(long, required = true)]
        intervals: Vec<PathBuf>,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        dbsnp: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 8)]
        threads: usize,
        #[arg(long)]
        emit_gvcf: bool,
    },
    /// Mark duplicates across coordinate-sorted chunk BAMs (GATK 4.6.2.0 MarkDuplicates semantics).
    Markdup {
        /// Input chunk BAMs, in any order.
        #[arg(long = "input", required = true)]
        inputs: Vec<PathBuf>,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        metrics: PathBuf,
    },
    /// Subread featureCounts 2.1.1 alignment counting (nf-core/rnaseq).
    Featurecounts {
        #[arg(short = 'a', long)]
        annotation: PathBuf,
        #[arg(short = 'o', long)]
        output: PathBuf,
        #[arg(short = 'g', long, default_value = "gene_id")]
        attribute: String,
        #[arg(short = 't', long, default_value = "exon")]
        feature_type: String,
        /// Declare paired input. Counts alignments, as in nf-core; not fragments.
        #[arg(short = 'p', long)]
        paired_end: bool,
        #[arg(short = 'B', long)]
        both_ends: bool,
        /// Inactive in alignment-counting mode, matching Subread.
        #[arg(short = 'C', long)]
        no_chimeric: bool,
        #[arg(short = 's', long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=2))]
        strand: u8,
        #[arg(short = 'T', long, default_value_t = 1)]
        threads: usize,
        /// Also emit nf-core biotype MultiQC tables (one input only).
        #[arg(long)]
        mqc_prefix: Option<PathBuf>,
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
    },
    /// Base quality recalibration (GATK 4.6.2.0 BaseRecalibrator + ApplyBQSR semantics).
    Bqsr {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long = "known-sites", required = true)]
        known_sites: Vec<PathBuf>,
        /// Target BED (the reference's `--intervals`).
        #[arg(long)]
        intervals: PathBuf,
        #[arg(long)]
        output_table: PathBuf,
        /// Recalibrated BAM (omit to only write the table).
        #[arg(long)]
        output_bam: Option<PathBuf>,
        /// Binary cache of parsed known sites (reused while the VCFs are unchanged).
        #[arg(long)]
        known_sites_cache: Option<PathBuf>,
    },
    /// Filter and chunk a FASTQ pair (fastp 1.1.0 semantics as run by sarek).
    Prep {
        #[arg(long)]
        fastq1: PathBuf,
        #[arg(long)]
        fastq2: PathBuf,
        #[arg(long)]
        outdir: PathBuf,
        /// Chunk file prefix, e.g. `HG002_WES-L1` -> `0001.HG002_WES-L1_R1.fastp.fastq.gz`.
        #[arg(long)]
        prefix: String,
        /// fastp's `--thread` value in the reference run (sets the chunk count).
        #[arg(long, default_value_t = 12)]
        fastp_threads: usize,
    },
    /// Duplicate marking then BQSR in one process: records stay in memory
    /// between the stages; md.bam is written alongside BQSR.
    Preprocess {
        /// Input chunk BAMs (coordinate sorted), in any order.
        #[arg(long = "input", required = true)]
        inputs: Vec<PathBuf>,
        #[arg(long)]
        md_output: PathBuf,
        #[arg(long)]
        metrics: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long = "known-sites", required = true)]
        known_sites: Vec<PathBuf>,
        #[arg(long)]
        intervals: PathBuf,
        #[arg(long)]
        known_sites_cache: Option<PathBuf>,
        #[arg(long)]
        output_table: PathBuf,
        #[arg(long)]
        output_bam: PathBuf,
    },
    /// Developer tools (fixtures, gates).
    #[command(hide = true, subcommand)]
    Dev(Dev),
}

#[derive(Subcommand)]
enum Dev {
    /// Split an intervals BED into sarek's scatter groups (CREATE_INTERVALS_BED).
    SplitIntervals {
        #[arg(long)]
        bed: PathBuf,
        #[arg(long)]
        outdir: PathBuf,
        #[arg(long, default_value_t = 200000.0)]
        nucleotides_per_second: f64,
    },
    /// Select paired FASTQ reads whose BWA-MEM alignments touch a contig.
    Slice {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        fastq1: PathBuf,
        #[arg(long)]
        fastq2: PathBuf,
        #[arg(long)]
        out1: PathBuf,
        #[arg(long)]
        out2: PathBuf,
        #[arg(long, default_value = "chr20")]
        contig: String,
        #[arg(long)]
        dense_sa_cache: PathBuf,
        #[arg(long, default_value_t = 4)]
        dense_sa_interval: u64,
        #[arg(long, default_value_t = 20)]
        threads: usize,
    },
    /// Build the tier0 fixture FASTQs from tier1 goldens (see src/tier0.rs).
    Tier0 {
        /// Golden htsjdk BAM from a tier1 reference run (md.bam).
        #[arg(long)]
        golden_bam: PathBuf,
        /// Target BED (the reference's `--intervals`).
        #[arg(long)]
        targets: PathBuf,
        #[arg(long, default_value = "chr20")]
        chrom: String,
        #[arg(long)]
        fastq1: PathBuf,
        #[arg(long)]
        fastq2: PathBuf,
        #[arg(long)]
        outdir: PathBuf,
        /// Keep 1 in N of the fastp-dropped and of the unmapped pairs.
        #[arg(long, default_value_t = 50)]
        sample_every: u64,
    },
    /// Compare fastp-style chunk FASTQs by decompressed content.
    GateChunks {
        #[arg(long)]
        ours: PathBuf,
        #[arg(long)]
        golden: PathBuf,
    },
    /// Compare a GATKReport (e.g. BQSR recal table) with a golden: bytes, then cells.
    GateReport {
        #[arg(long)]
        ours: PathBuf,
        #[arg(long)]
        golden: PathBuf,
        #[arg(long, default_value_t = 20)]
        examples: usize,
    },
    /// Compare a BAM with a golden under the GATES.md normalization.
    GateBam {
        #[arg(long)]
        ours: PathBuf,
        #[arg(long)]
        golden: PathBuf,
        #[arg(long, default_value_t = 10)]
        examples: usize,
    },
    /// Compare VCF or GVCF headers and records under the GATES.md normalization.
    #[command(alias = "gate-gvcf")]
    GateVcf {
        #[arg(long)]
        ours: PathBuf,
        #[arg(long)]
        golden: PathBuf,
        #[arg(long, default_value_t = 10)]
        examples: usize,
    },
    /// Write htsjdk's index for a BAM (`.bai`) or BGZF VCF (`.tbi`), as GATK
    /// writes it alongside the file.
    HtsjdkIndex {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// Convert a BAM to CRAM 3.0 (+ .crai) against a FASTA (records gate:
    /// scripts/cram_gate.sh).
    BamToCram {
        input: PathBuf,
        reference: PathBuf,
        output: PathBuf,
        /// Encoding threads (default: all cores).
        #[arg(long)]
        threads: Option<usize>,
        #[arg(long, default_value_t = 10_000)]
        records_per_slice: usize,
        /// Also time writing the same records as BAM (our default BGZF
        /// compression) to this path, for comparison.
        #[arg(long)]
        bam_bench: Option<PathBuf>,
        /// What the CRAM decodes to: lossless (the BAM), sarek-md or
        /// sarek-recal (sarek 3.10.0's md.cram / recal.cram).
        #[arg(long, default_value = "lossless")]
        profile: String,
        /// With sarek-recal: the input is the merge of this many interval
        /// groups (BAM mode, as recal.bam); write sarek's CRAM-mode merge
        /// (its ID suffixes, PG:Z after the other tags).
        #[arg(long, default_value_t = 1)]
        groups: usize,
    },
    /// Emit normalized SHA-256 digests for every tier1 stage artifact.
    Digest(DigestArgs),
    /// Compare two digest manifests, printing PASS/FAIL for every artifact.
    DigestCompare {
        #[arg(long)]
        ours: PathBuf,
        #[arg(long)]
        golden: PathBuf,
    },
}

#[derive(ClapArgs)]
struct DigestArgs {
    /// Standard `nomics run --write-intermediates` output directory.
    #[arg(long)]
    run_dir: Option<PathBuf>,
    /// Override the fastp chunk directory (required without --run-dir).
    #[arg(long)]
    chunks: Option<PathBuf>,
    /// Override the mapped per-chunk BAM directory.
    #[arg(long)]
    mapped: Option<PathBuf>,
    #[arg(long)]
    md_bam: Option<PathBuf>,
    #[arg(long)]
    md_metrics: Option<PathBuf>,
    #[arg(long)]
    recal_table: Option<PathBuf>,
    #[arg(long)]
    recal_bam: Option<PathBuf>,
    #[arg(long)]
    vcf: Option<PathBuf>,
    /// GVCF-mode runs: the HaplotypeCaller gVCF.
    #[arg(long)]
    gvcf: Option<PathBuf>,
    /// GVCF-mode runs: `joint_variant_calling/joint_germline.vcf.gz`.
    #[arg(long)]
    joint_vcf: Option<PathBuf>,
    #[arg(long, default_value = "HG002_WES")]
    prefix: String,
    /// Write the manifest here instead of stdout.
    #[arg(long)]
    output: Option<PathBuf>,
    /// Reuse digests of unchanged files (keyed by device, inode, size and
    /// mtime) from this cache (default: ~/omics/cache/digests.tsv).
    #[arg(long)]
    cache: Option<PathBuf>,
    /// Digest everything afresh, without reading or writing the cache.
    #[arg(long)]
    no_cache: bool,
}

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// `nomics dev bam-to-cram`: stream a BAM into a CRAM and its `.crai`.
#[allow(clippy::too_many_arguments)]
fn bam_to_cram(
    input: &Path,
    reference: &Path,
    output: &Path,
    threads: usize,
    records_per_slice: usize,
    bam_bench: Option<&Path>,
    profile: &str,
    groups: usize,
) -> std::io::Result<()> {
    use std::io::Write;
    let t = Instant::now();
    let fasta = nomics_io::fasta::Fasta::load(reference, threads)?;
    let t_ref = t.elapsed();
    let t = Instant::now();
    nomics_io::bamstream::warn_if_incomplete(input);
    let mut stream =
        nomics_io::bamstream::BamStream::open(input, 10_000, 4 * threads, threads.min(8))?;
    let header = stream.header.clone();
    let name = output
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let mut profile = nomics_io::cram::Profile::parse(profile)
        .ok_or_else(|| std::io::Error::other("--profile: lossless, sarek-md or sarek-recal"))?;
    let merge = (groups > 1)
        .then(|| {
            if profile != nomics_io::cram::Profile::SarekRecal {
                return Err(std::io::Error::other("--groups needs --profile sarek-recal"));
            }
            profile = nomics_io::cram::Profile::SarekRecalMerged;
            let fasta_name = reference.file_name().unwrap_or_default().to_string_lossy();
            let stem = name.strip_suffix(".recal.cram").unwrap_or("sample");
            let cl = format!("samtools view -Ch -T {fasta_name} -o {stem}.md.cram {stem}.md.bam");
            nomics_bqsr::CramMerge::from_recal_header(&String::from_utf8_lossy(&header.text), groups, &cl)
                .map_err(std::io::Error::other)
        })
        .transpose()?;
    let header = match &merge {
        Some(m) => nomics_io::bam::Header {
            text: m.header_text.clone().into_bytes(),
            refs: header.refs.clone(),
        },
        None => header,
    };
    let opts = nomics_io::cram::CramOptions {
        records_per_slice,
        threads,
        ..Default::default()
    }
    .with_file_name(&name);
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let (header, opts) = profile.prepare(
        &header,
        reference,
        Some((env!("CARGO_PKG_VERSION"), &command_line)),
        opts,
    )?;
    let file = std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(output)?);
    let (file, summary) = nomics_io::cram::write_cram(file, &header, &fasta, &opts, |sink| {
        let mut buf = Vec::new();
        while let Some(batch) = stream.next_batch()? {
            for i in 0..batch.len() {
                match &merge {
                    Some(m) => {
                        buf.clear();
                        m.translate(batch.record(i), &mut buf);
                        sink.push(&buf[4..])?;
                    }
                    None => sink.push(batch.record(i))?,
                }
            }
        }
        Ok(())
    })?;
    file.into_inner().map_err(|e| e.into_error())?.sync_all()?;
    let mut crai = output.as_os_str().to_owned();
    crai.push(".crai");
    std::fs::write(&crai, nomics_io::cram::crai_bytes(&summary.index))?;
    let t_cram = t.elapsed();
    eprintln!(
        "bam-to-cram: {} records, {} containers, {} bytes; reference load {:.2} s, encode {:.2} s ({threads} threads)",
        summary.records,
        summary.containers,
        summary.bytes,
        t_ref.as_secs_f64(),
        t_cram.as_secs_f64()
    );
    if let Some(bam_out) = bam_bench {
        let t = Instant::now();
        let mut stream =
            nomics_io::bamstream::BamStream::open(input, 10_000, 4 * threads, threads.min(8))?;
        let f = std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(bam_out)?);
        let mut w = nomics_io::bgzf::ParallelBgzfWriter::for_bam(f, threads);
        stream.header.write(&mut w)?;
        while let Some(batch) = stream.next_batch()? {
            for i in 0..batch.len() {
                nomics_io::bam::write_raw_record(&mut w, batch.record(i))?;
            }
        }
        let mut f = w.finish()?;
        f.flush()?;
        f.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        eprintln!(
            "bam-to-cram: BAM rewrite {:.2} s ({threads} threads), {} bytes",
            t.elapsed().as_secs_f64(),
            std::fs::metadata(bam_out)?.len()
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    // Multi-call: invoked as fastp/bwa/gatk/samtools through a shim symlink.
    if let Some(tool) = std::env::args_os()
        .next()
        .as_deref()
        .and_then(shim::tool_name)
    {
        return shim::main(tool);
    }
    match Cli::parse().command {
        Command::Dev(Dev::BamToCram {
            input,
            reference,
            output,
            threads,
            records_per_slice,
            bam_bench,
            profile,
            groups,
        }) => {
            let threads = threads
                .unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));
            match bam_to_cram(
                &input,
                &reference,
                &output,
                threads,
                records_per_slice,
                bam_bench.as_deref(),
                &profile,
                groups,
            ) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("bam-to-cram: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Dev(Dev::HtsjdkIndex { input, output }) => {
            let name = input.to_string_lossy();
            let result = if name.ends_with(".bam") {
                nomics_io::htsjdk_index::write_bai_htsjdk(&input, &output)
            } else if name.ends_with(".vcf.gz") {
                nomics_io::htsjdk_index::write_tbi_htsjdk_vcf(&input, &output)
            } else {
                Err(std::io::Error::other("input must be .bam or .vcf.gz"))
            };
            match result {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("htsjdk-index: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::JointGenotype {
            gvcf,
            dbsnp,
            output,
            threads,
        } => {
            let start = Instant::now();
            let result = std::fs::File::create(&output).and_then(|f| {
                if output.extension().is_some_and(|e| e == "gz") {
                    let mut w = nomics_io::bgzf::ParallelBgzfWriter::new(f, threads);
                    let n = nomics_joint::genotype_gvcf(&gvcf, dbsnp.as_deref(), threads, &mut w)?;
                    w.finish()?;
                    Ok(n)
                } else {
                    let mut w = std::io::BufWriter::new(f);
                    let n = nomics_joint::genotype_gvcf(&gvcf, dbsnp.as_deref(), threads, &mut w)?;
                    std::io::Write::flush(&mut w)?;
                    Ok(n)
                }
            });
            match result {
                Ok(n) => {
                    eprintln!(
                        "joint-genotype: {n} records; {:.3}s",
                        start.elapsed().as_secs_f64()
                    );
                    ExitCode::SUCCESS
                }
                Err(error) => {
                    eprintln!("joint-genotype: {error}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::SoloCount {
            bam,
            whitelist,
            features,
            output,
            threads,
        } => {
            match nomics_solo::count(&nomics_solo::Config {
                bam,
                whitelist,
                features,
                output,
                threads,
            }) {
                Ok(report) => {
                    eprintln!(
                        "solo-count: {} reads, {} counted, {} UMIs, {} cells; {:.3}s",
                        report.reads,
                        report.counted_reads,
                        report.umis,
                        report.cells,
                        report.elapsed_seconds
                    );
                    ExitCode::SUCCESS
                }
                Err(error) => {
                    eprintln!("solo-count: {error}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::SalmonQuant {
            gene_map,
            threads,
            transcripts,
            alignments,
            lib_type,
            incompat_prior,
            output,
        } => {
            let start = Instant::now();
            let options = nomics_quant::Options {
                gene_map,
                threads,
                transcripts,
                alignments,
                lib_type,
                incompat_prior,
                output,
            };
            match nomics_quant::run(&options) {
                Ok(q) => {
                    eprintln!(
                        "salmon-quant: {} / {} fragments, {} classes, {} VBEM iterations in {:.3}s",
                        q.num_mapped,
                        q.num_processed,
                        q.num_eq_classes,
                        q.iterations,
                        start.elapsed().as_secs_f64()
                    );
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("salmon-quant: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Dev(Dev::SplitIntervals {
            bed,
            outdir,
            nucleotides_per_second,
        }) => {
            let split = std::fs::read_to_string(&bed).and_then(|text| {
                let groups = scatter::split_bed(&text, nucleotides_per_second)?;
                scatter::write_groups(&groups, &outdir)
            });
            match split {
                Ok(paths) => {
                    for p in paths {
                        println!("{}", p.display());
                    }
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("split-intervals: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Dev(Dev::Slice {
            index,
            fastq1,
            fastq2,
            out1,
            out2,
            contig,
            dense_sa_cache,
            dense_sa_interval,
            threads,
        }) => match nomics_align::slice::run(&nomics_align::slice::Args {
            index,
            fastq1,
            fastq2,
            out1,
            out2,
            contig,
            dense_sa_cache,
            dense_sa_interval,
            threads,
        }) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("slice: {e}");
                ExitCode::FAILURE
            }
        },
        Command::Selftest { mib } => selftest(mib),
        Command::Dev(Dev::Tier0 {
            golden_bam,
            targets,
            chrom,
            fastq1,
            fastq2,
            outdir,
            sample_every,
        }) => {
            let args = tier0::Args {
                golden_bam,
                targets,
                chrom,
                fastq1,
                fastq2,
                outdir,
                sample_every,
            };
            match tier0::run(&args) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("tier0: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Dev(Dev::GateReport {
            ours,
            golden,
            examples,
        }) => {
            let read = |p: &PathBuf| {
                std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))
            };
            let result = read(&ours)
                .and_then(|a| read(&golden).map(|b| (a, b)))
                .and_then(|(a, b)| {
                    nomics_gate::compare_gatk_reports(&a, &b, examples).map_err(|e| e.to_string())
                });
            match result {
                Ok(c) if c.bytes_identical => {
                    println!("verdict: IDENTICAL (bytes)");
                    ExitCode::SUCCESS
                }
                Ok(c) => {
                    for d in &c.diffs {
                        println!("  {d}");
                    }
                    println!("verdict: DIFFERENT");
                    ExitCode::FAILURE
                }
                Err(e) => {
                    eprintln!("gate-report: {e}");
                    ExitCode::from(2)
                }
            }
        }
        Command::Dev(Dev::GateBam {
            ours,
            golden,
            examples,
        }) => match nomics_gate::compare_bams(&ours, &golden, examples) {
            Ok(report) => {
                println!("{report}");
                if report.identical() {
                    ExitCode::SUCCESS
                } else {
                    ExitCode::FAILURE
                }
            }
            Err(e) => {
                eprintln!("gate-bam: {e}");
                ExitCode::from(2)
            }
        },
        Command::Dev(Dev::GateVcf {
            ours,
            golden,
            examples,
        }) => match nomics_gate::compare_vcfs(&ours, &golden, examples) {
            Ok(report) => {
                println!("{report}");
                if report.identical() {
                    ExitCode::SUCCESS
                } else {
                    ExitCode::FAILURE
                }
            }
            Err(e) => {
                eprintln!("gate-vcf: {e}");
                ExitCode::from(2)
            }
        },
        Command::Dev(Dev::Digest(args)) => match write_digest_manifest(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("digest: {e}");
                ExitCode::from(2)
            }
        },
        Command::Dev(Dev::DigestCompare { ours, golden }) => {
            match compare_digest_manifests(&ours, &golden) {
                Ok(true) => ExitCode::SUCCESS,
                Ok(false) => ExitCode::FAILURE,
                Err(e) => {
                    eprintln!("digest-compare: {e}");
                    ExitCode::from(2)
                }
            }
        }
        Command::Featurecounts {
            annotation,
            output,
            attribute,
            feature_type,
            paired_end,
            both_ends,
            no_chimeric,
            strand,
            threads,
            mqc_prefix,
            inputs,
        } => {
            // Subread records argv in order, including a trailing space. Use its
            // executable label, and omit only our extra MultiQC option.
            let mut command = vec!["featureCounts".to_owned()];
            let mut args = std::env::args().skip(2);
            while let Some(arg) = args.next() {
                if arg == "--mqc-prefix" {
                    args.next();
                } else if !arg.starts_with("--mqc-prefix=") {
                    command.push(arg);
                }
            }
            let options = nomics_count::Options {
                annotation,
                inputs,
                output,
                attribute,
                feature_type,
                command,
                counting: nomics_count::CountOptions {
                    paired_end,
                    both_ends,
                    no_chimeric,
                    strand,
                    threads,
                },
            };
            let result = (|| -> std::io::Result<()> {
                if mqc_prefix.is_some() && options.inputs.len() != 1 {
                    return Err(std::io::Error::other(
                        "MultiQC requires exactly one input BAM",
                    ));
                }
                let counts = nomics_count::run(&options)?;
                if let Some(prefix) = mqc_prefix {
                    let annotation = nomics_count::Annotation::read(
                        &options.annotation,
                        &options.attribute,
                        &options.feature_type,
                    )?;
                    let sample = options.inputs[0]
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy();
                    let sample = sample
                        .strip_suffix(".markdup.sorted.bam")
                        .unwrap_or(&sample);
                    nomics_count::write_multiqc(&annotation, &counts[0], sample, &prefix)?;
                }
                Ok(())
            })();
            match result {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("featurecounts: {e}");
                    ExitCode::from(2)
                }
            }
        }
        Command::Markdup {
            inputs,
            output,
            metrics,
        } => {
            let command_line = format!(
                "nomics markdup {} --output {} --metrics {}",
                inputs
                    .iter()
                    .map(|p| format!("--input {}", p.display()))
                    .collect::<Vec<_>>()
                    .join(" "),
                output.display(),
                metrics.display()
            );
            let opts = nomics_markdup::Options {
                inputs,
                output,
                metrics,
                pg_version: format!("Version:{}", env!("CARGO_PKG_VERSION")),
                pg_command_line: command_line,
                threads: std::thread::available_parallelism().map_or(4, |n| n.get()),
            };
            let t = Instant::now();
            match nomics_markdup::run(&opts) {
                Ok(m) => {
                    eprintln!("markdup: {m:?} in {:.1}s", t.elapsed().as_secs_f64());
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("markdup: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Bqsr {
            input,
            reference,
            known_sites,
            intervals,
            output_table,
            output_bam,
            known_sites_cache,
        } => {
            let command_line = format!(
                "nomics bqsr --input {} --output-table {}",
                input.display(),
                output_table.display()
            );
            let opts = nomics_bqsr::Options {
                input,
                reference,
                known_sites,
                intervals,
                scatter: None,
                known_sites_cache,
                output_table,
                output_bam,
                pg_command_line: command_line,
                threads: std::thread::available_parallelism().map_or(4, |n| n.get()),
            };
            let t = Instant::now();
            match nomics_bqsr::run(&opts) {
                Ok(()) => {
                    eprintln!("bqsr: done in {:.1}s", t.elapsed().as_secs_f64());
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("bqsr: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Preprocess {
            inputs,
            md_output,
            metrics,
            reference,
            known_sites,
            intervals,
            known_sites_cache,
            output_table,
            output_bam,
        } => {
            let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
            let t = Instant::now();
            let md_cl = format!(
                "nomics preprocess {} --md-output {}",
                inputs
                    .iter()
                    .map(|p| format!("--input {}", p.display()))
                    .collect::<Vec<_>>()
                    .join(" "),
                md_output.display()
            );
            let md = nomics_markdup::Options {
                inputs,
                output: md_output.clone(),
                metrics,
                pg_version: format!("Version:{}", env!("CARGO_PKG_VERSION")),
                pg_command_line: md_cl,
                threads,
            };
            let bqsr = nomics_bqsr::Options {
                input: md_output.clone(),
                reference,
                known_sites,
                intervals,
                scatter: None,
                known_sites_cache,
                pg_command_line: format!(
                    "nomics preprocess --output-table {}",
                    output_table.display()
                ),
                output_table,
                output_bam: Some(output_bam),
                threads,
            };
            let result = nomics_markdup::mark(&md).and_then(|(set, m)| {
                eprintln!(
                    "preprocess: markdup done at {:.1}s ({m:?})",
                    t.elapsed().as_secs_f64()
                );
                // md.bam is a published output; write it while BQSR works on the same records.
                std::thread::scope(|scope| {
                    let writer = scope.spawn(|| {
                        let tw = Instant::now();
                        set.write(&md_output, threads).map(|()| tw.elapsed())
                    });
                    let recal = nomics_bqsr::run_on(&set, &bqsr);
                    let written = writer.join().expect("md.bam writer panicked");
                    if let Ok(d) = &written {
                        eprintln!(
                            "preprocess: md.bam written in {:.1}s (overlapped)",
                            d.as_secs_f64()
                        );
                    }
                    recal.and(written.map(|_| ()))
                })
            });
            match result {
                Ok(()) => {
                    eprintln!("preprocess: done in {:.1}s", t.elapsed().as_secs_f64());
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("preprocess: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Prep {
            fastq1,
            fastq2,
            outdir,
            prefix,
            fastp_threads,
        } => {
            let t = Instant::now();
            let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
            let result =
                nomics_prep::run(&fastq1, &fastq2, fastp_threads).and_then(|(chunks, stats)| {
                    eprintln!(
                        "prep: filtered in {:.1}s: {stats:?}",
                        t.elapsed().as_secs_f64()
                    );
                    nomics_prep::write_chunks(&chunks, &outdir, &prefix, threads)
                });
            match result {
                Ok(()) => {
                    eprintln!("prep: done in {:.1}s", t.elapsed().as_secs_f64());
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("prep: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Dev(Dev::GateChunks { ours, golden }) => {
            match nomics_gate::compare_fastq_chunks(&ours, &golden) {
                Ok((lines, all)) => {
                    for l in lines {
                        println!("{l}");
                    }
                    println!("verdict: {}", if all { "IDENTICAL" } else { "DIFFERENT" });
                    if all {
                        ExitCode::SUCCESS
                    } else {
                        ExitCode::FAILURE
                    }
                }
                Err(e) => {
                    eprintln!("gate-chunks: {e}");
                    ExitCode::from(2)
                }
            }
        }
        Command::Hc {
            input,
            intervals,
            reference,
            dbsnp,
            output,
            threads,
            emit_gvcf,
        } => {
            match pipeline::haplotypecaller(
                &input, &intervals, &reference, &dbsnp, &output, threads, emit_gvcf,
            ) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("hc: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Command::Run {
            fastq1,
            fastq2,
            lanes,
            reference,
            targets,
            outdir,
            dbsnp,
            known_sites,
            known_sites_cache,
            dense_sa_cache,
            dense_sa_interval,
            scatter_order,
            read_group_ds,
            sample,
            prefix,
            threads,
            write_intermediates,
            spill_prep,
            save_mapped,
            save_output_as_bam,
            keep_bams,
            resume_from,
            hc_only,
            emit_gvcf,
        } => match pipeline::run(pipeline::Args {
            emit_gvcf,
            fastq1,
            fastq2,
            lanes,
            reference,
            targets,
            outdir,
            dbsnp,
            known_sites,
            known_sites_cache,
            dense_sa_cache,
            dense_sa_interval,
            scatter_order,
            read_group_ds,
            sample,
            prefix,
            threads,
            write_intermediates,
            spill_prep,
            save_mapped,
            save_output_as_bam,
            keep_bams,
            resume_from: if hc_only {
                Some(pipeline::Resume::Hc)
            } else {
                resume_from
            },
        }) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("run: {e}");
                ExitCode::FAILURE
            }
        },
    }
}

fn digest_path(
    explicit: &Option<PathBuf>,
    run_dir: &Option<PathBuf>,
    relative: impl AsRef<Path>,
    flag: &str,
) -> Result<PathBuf, String> {
    explicit
        .clone()
        .or_else(|| run_dir.as_ref().map(|p| p.join(relative)))
        .ok_or_else(|| format!("{flag} is required when --run-dir is not supplied"))
}

fn matching_files(
    dir: &Path,
    mut matches: impl FnMut(&str) -> bool,
) -> Result<Vec<(String, PathBuf)>, String> {
    let entries = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let mut files = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|e| format!("{}: {e}", dir.display()))?;
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| format!("{} contains a non-UTF-8 file name", dir.display()))?;
        if matches(&name) {
            files.push((name, entry.path()));
        }
    }
    files.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(files)
}

fn write_digest_manifest(args: &DigestArgs) -> Result<(), String> {
    let chunks = digest_path(&args.chunks, &args.run_dir, "prep", "--chunks")?;
    let mapped = digest_path(&args.mapped, &args.run_dir, "mapped", "--mapped")?;
    let named = |explicit: &Option<PathBuf>, suffix: &str, flag: &str| {
        digest_path(
            explicit,
            &args.run_dir,
            format!("{}.{suffix}", args.prefix),
            flag,
        )
    };
    let md_bam = named(&args.md_bam, "md.bam", "--md-bam")?;
    let md_metrics = named(&args.md_metrics, "md.metrics", "--md-metrics")?;
    let recal_table = named(&args.recal_table, "recal.table", "--recal-table")?;
    let recal_bam = named(&args.recal_bam, "recal.bam", "--recal-bam")?;
    // A GVCF-mode run writes the gVCF and the joint VCF instead of the VCF.
    let gvcf_mode = args.gvcf.is_some()
        || args.vcf.is_none()
            && args.run_dir.as_ref().is_some_and(|d| {
                d.join(format!("{}.haplotypecaller.g.vcf.gz", args.prefix))
                    .exists()
            });
    let vcfs = if gvcf_mode {
        vec![
            (
                "gvcf",
                named(&args.gvcf, "haplotypecaller.g.vcf.gz", "--gvcf")?,
            ),
            (
                "joint_vcf",
                digest_path(
                    &args.joint_vcf,
                    &args.run_dir,
                    "joint_variant_calling/joint_germline.vcf.gz",
                    "--joint-vcf",
                )?,
            ),
        ]
    } else {
        vec![("vcf", named(&args.vcf, "haplotypecaller.vcf.gz", "--vcf")?)]
    };

    let chunk_files = matching_files(&chunks, |name| name.ends_with(".fastp.fastq.gz"))?;
    // No chunk FASTQs at all is a run that did not write them (no
    // `--write-intermediates`, or a resume without prep); digest-compare
    // then reports each golden chunk as missing.
    if chunk_files.len() % 2 != 0 {
        return Err(format!(
            "{}: expected R1/R2 chunk FASTQ pairs, found {}",
            chunks.display(),
            chunk_files.len()
        ));
    }
    let mapped_prefix = format!("{}.", args.prefix);
    let bam_files = matching_files(&mapped, |name| {
        name.strip_prefix(&mapped_prefix)
            .and_then(|s| s.strip_suffix(".bam"))
            .is_some_and(|number| number.len() == 4 && number.bytes().all(|b| b.is_ascii_digit()))
    })?;
    // No mapped BAMs at all is a run that did not publish them
    // (`--save-mapped false`); digest-compare then reports each golden
    // mapped BAM as missing.
    if !bam_files.is_empty() && !chunk_files.is_empty() && bam_files.len() * 2 != chunk_files.len() {
        return Err(format!(
            "{}: {} mapped BAMs for prefix {} but {} chunk FASTQs",
            mapped.display(),
            bam_files.len(),
            args.prefix,
            chunk_files.len()
        ));
    }
    // (artifact, file, digest kind); kinds name the normalization applied.
    type DigestFn = fn(&Path) -> std::io::Result<nomics_gate::ArtifactDigest>;
    let mut tasks: Vec<(String, PathBuf, &str, DigestFn)> = Vec::new();
    for (name, path) in chunk_files {
        tasks.push((
            format!("chunk/{name}"),
            path,
            "fastq_chunk",
            nomics_gate::digest_fastq_chunk,
        ));
    }
    for (name, path) in bam_files {
        tasks.push((
            format!("mapped/{name}"),
            path,
            "bam",
            nomics_gate::digest_bam,
        ));
    }
    // A `--save-output-as-bam false` run publishes md.cram/recal.cram
    // instead (gated by scripts/cram_gate.sh); digest-compare then reports
    // the golden BAMs as missing from ours.
    let cram_instead = |bam: &Path, kind: &str| {
        !bam.exists()
            && args.run_dir.as_ref().is_some_and(|d| {
                pipeline::sarek_cram_path(d, &args.prefix, kind).exists()
            })
    };
    if !cram_instead(&md_bam, "md") {
        tasks.push(("md_bam".to_owned(), md_bam, "bam", nomics_gate::digest_bam));
    }
    if !cram_instead(&recal_bam, "recal") {
        tasks.push((
            "recal_bam".to_owned(),
            recal_bam,
            "bam",
            nomics_gate::digest_bam,
        ));
    }
    tasks.push((
        "md_metrics".to_owned(),
        md_metrics,
        "md_metrics",
        nomics_gate::digest_md_metrics,
    ));
    tasks.push((
        "recal_table".to_owned(),
        recal_table,
        "raw",
        nomics_gate::digest_raw,
    ));
    for (artifact, path) in vcfs {
        tasks.push((artifact.to_owned(), path, "vcf", nomics_gate::digest_vcf));
    }

    let cache_path = (!args.no_cache).then(|| {
        args.cache.clone().unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap_or_default())
                .join("omics/cache/digests.tsv")
        })
    });
    let mut cache = match &cache_path {
        Some(p) => {
            nomics_gate::DigestCache::load(p).map_err(|e| format!("{}: {e}", p.display()))?
        }
        None => nomics_gate::DigestCache::default(),
    };
    // Cached digests first; the rest run on a few threads (each digest is
    // a single decompress-and-hash stream).
    let mut rows: Vec<Option<(String, nomics_gate::ArtifactDigest)>> = vec![None; tasks.len()];
    let mut todo = Vec::new();
    for (i, (artifact, path, kind, _)) in tasks.iter().enumerate() {
        match cache.get(kind, path) {
            Some(d) => rows[i] = Some((artifact.clone(), d)),
            None => todo.push(i),
        }
    }
    let cached = tasks.len() - todo.len();
    let queue = std::sync::Mutex::new(todo);
    let results = std::sync::Mutex::new(Vec::new());
    let workers = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(12);
    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| {
                loop {
                    let Some(i) = queue.lock().unwrap().pop() else {
                        break;
                    };
                    let (_, path, _, digest) = &tasks[i];
                    let r = digest(path).map_err(|e| format!("{}: {e}", path.display()));
                    results.lock().unwrap().push((i, r));
                }
            });
        }
    });
    for (i, r) in results.into_inner().unwrap() {
        let d = r?;
        let (artifact, path, kind, _) = &tasks[i];
        cache.insert(kind, path, &d);
        rows[i] = Some((artifact.clone(), d));
    }
    if let Some(p) = &cache_path {
        cache.save(p).map_err(|e| format!("{}: {e}", p.display()))?;
    }
    eprintln!(
        "digest: {} artifacts, {cached} from the cache, {} computed",
        tasks.len(),
        tasks.len() - cached
    );
    let mut rows: Vec<(String, nomics_gate::ArtifactDigest)> = rows
        .into_iter()
        .map(|r| r.expect("every task digested"))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));

    let mut manifest = String::from("artifact\tsha256\tcount\n");
    for (artifact, digest) in &rows {
        manifest.push_str(&format!(
            "{artifact}\t{}\t{}\n",
            digest.sha256, digest.count
        ));
    }
    if let Some(path) = &args.output {
        std::fs::write(path, manifest).map_err(|e| format!("{}: {e}", path.display()))?;
        eprintln!(
            "digest: wrote {} artifacts to {}",
            rows.len(),
            path.display()
        );
    } else {
        print!("{manifest}");
    }
    Ok(())
}

fn read_digest_manifest(path: &Path) -> Result<BTreeMap<String, (String, u64)>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut lines = text.lines();
    if lines.next() != Some("artifact\tsha256\tcount") {
        return Err(format!(
            "{}: invalid or missing manifest header",
            path.display()
        ));
    }
    let mut manifest = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 3
            || fields[0].is_empty()
            || fields[1].len() != 64
            || !fields[1]
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(format!(
                "{}:{}: invalid manifest row",
                path.display(),
                index + 2
            ));
        }
        let count = fields[2]
            .parse::<u64>()
            .map_err(|e| format!("{}:{}: invalid count: {e}", path.display(), index + 2))?;
        if manifest
            .insert(fields[0].to_owned(), (fields[1].to_owned(), count))
            .is_some()
        {
            return Err(format!(
                "{}:{}: duplicate artifact {}",
                path.display(),
                index + 2,
                fields[0]
            ));
        }
    }
    if manifest.is_empty() {
        return Err(format!("{}: manifest has no artifacts", path.display()));
    }
    Ok(manifest)
}

fn compare_digest_manifests(ours: &Path, golden: &Path) -> Result<bool, String> {
    let ours = read_digest_manifest(ours)?;
    let golden = read_digest_manifest(golden)?;
    let names: BTreeSet<&String> = ours.keys().chain(golden.keys()).collect();
    let mut all = true;
    for name in names {
        match (ours.get(name), golden.get(name)) {
            (Some(a), Some(b)) if a == b => println!("PASS\t{name}\tcount={}", a.1),
            (Some(a), Some(b)) => {
                all = false;
                println!(
                    "FAIL\t{name}\tours={} count={}\tgolden={} count={}",
                    a.0, a.1, b.0, b.1
                );
            }
            (Some(_), None) => {
                all = false;
                println!("FAIL\t{name}\tmissing from golden");
            }
            (None, Some(_)) => {
                all = false;
                println!("FAIL\t{name}\tmissing from ours");
            }
            (None, None) => unreachable!(),
        }
    }
    println!("verdict: {}", if all { "PASS" } else { "FAIL" });
    Ok(all)
}

fn selftest(mib: usize) -> ExitCode {
    let mut buf = match HostBuf::new(mib << 20) {
        Ok(buf) => buf,
        Err(e) => {
            eprintln!("allocation failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let t = Instant::now();
    for (i, w) in buf.as_u64s_mut().iter_mut().enumerate() {
        *w = i as u64;
    }
    let fill = t.elapsed();
    let n = buf.as_u64s().len() as u64;
    let gib = buf.len() as f64 / (1u64 << 30) as f64;
    println!(
        "host buffer: {mib} MiB, filled in {:.1} ms ({:.1} GiB/s)",
        fill.as_secs_f64() * 1e3,
        gib / fill.as_secs_f64()
    );
    match buf.huge_page_bytes() {
        Some(h) => println!("huge-page backed: {} of {} MiB", h >> 20, buf.len() >> 20),
        None => println!("huge-page backing: not observable on this OS"),
    }

    let t = Instant::now();
    let cpu: u64 = buf.as_u64s().iter().copied().fold(0u64, u64::wrapping_add);
    let cpu_s = t.elapsed().as_secs_f64();
    println!("cpu sum (1 thread): {:.1} GiB/s", gib / cpu_s);

    match nomics_gpu::device_info(0) {
        Ok(info) => {
            #[cfg(feature = "metal")]
            let arch = format!(
                "Apple{} family, {} GPU cores",
                info.compute_capability.0, info.multiprocessors
            );
            #[cfg(not(feature = "metal"))]
            let arch = format!(
                "sm_{}{} x{} SMs",
                info.compute_capability.0, info.compute_capability.1, info.multiprocessors
            );
            println!(
                "gpu: {} {arch}, pageable memory access: {}",
                info.name, info.pageable_memory_access
            );
            println!(
                "gpu memory: {:.1} GiB working set, largest buffer {}",
                info.total_global_mem as f64 / (1u64 << 30) as f64,
                if info.max_buffer_bytes == usize::MAX {
                    "unlimited".to_owned()
                } else {
                    format!(
                        "{:.1} GiB",
                        info.max_buffer_bytes as f64 / (1u64 << 30) as f64
                    )
                }
            );
        }
        Err(e) => {
            println!("gpu: {e}");
            return ExitCode::SUCCESS;
        }
    }
    // First kernel runs at cold clocks; report the best of a few runs.
    let mut best = f32::MAX;
    let mut best_wall = f64::MAX;
    for _ in 0..5 {
        let t = Instant::now();
        let result = nomics_gpu::sum_u64_in_place(buf.as_u64s());
        best_wall = best_wall.min(t.elapsed().as_secs_f64());
        match result {
            Ok((sum, ms)) => {
                if sum != cpu || sum != n.wrapping_mul(n.wrapping_sub(1)) / 2 {
                    eprintln!("gpu sum mismatch: gpu {sum}, cpu {cpu}");
                    return ExitCode::FAILURE;
                }
                best = best.min(ms);
            }
            Err(e) => {
                eprintln!("gpu sum failed: {e}");
                return ExitCode::FAILURE;
            }
        }
    }
    println!(
        "gpu in-place sum: {:.1} GiB/s (best of 5, matches cpu); call {:.2} ms for {:.2} ms of GPU time",
        gib / (best as f64 / 1e3),
        best_wall * 1e3,
        best
    );
    ExitCode::SUCCESS
}
