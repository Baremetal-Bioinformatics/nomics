use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender, SyncSender};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use nomics_align::index::Index;
use nomics_align::mem::{align_chunks_bam, align_chunks_bam_files, align_stream};
use nomics_align::samsort::{header, sort_sam};
use nomics_hc::caller::{Context, header_meta};
use nomics_hc::dbsnp::Dbsnp;
use nomics_hc::reads::{GroupCollector, Intervals, read_stream};
use nomics_hc::vcf::{Header, encode, write_vcf_gz};
use nomics_io::bam::RecordSet;
use nomics_io::bgzf::ParallelBgzfWriter;
use nomics_io::fasta::Fasta;

/// fastp `--thread 12`: one chunk per worker.
const PREP_CHUNKS: usize = 12;
/// Gzipped FASTQ bytes above which prep spills its chunks to disk (the
/// in-memory chunk text is ~5x the gzipped size).
const SPILL_PREP_ABOVE: u64 = 16 << 30;

pub struct Args {
    /// One FASTQ pair per lane, in samplesheet order.
    pub fastq1: Vec<PathBuf>,
    pub fastq2: Vec<PathBuf>,
    /// The samplesheet's lane names (default `L1`, `L2`, ...).
    pub lanes: Vec<String>,
    pub reference: PathBuf,
    pub targets: PathBuf,
    pub outdir: PathBuf,
    pub dbsnp: Option<PathBuf>,
    pub known_sites: Vec<PathBuf>,
    pub known_sites_cache: Option<PathBuf>,
    pub dense_sa_cache: PathBuf,
    pub dense_sa_interval: u64,
    /// The reference run's gather/merge orders (see `scatter::orders`).
    pub scatter_order: Option<PathBuf>,
    /// The read group's `DS:` (sarek writes the FASTA path as passed to it,
    /// so a reference run elsewhere has another path); default: our FASTA.
    pub read_group_ds: Option<String>,
    pub sample: String,
    pub prefix: String,
    pub threads: usize,
    pub write_intermediates: bool,
    pub spill_prep: bool,
    /// Publish the per-chunk sorted BAMs in `<outdir>/mapped` (sarek
    /// `--save_mapped`). Off, whole genomes keep markdup's chunk spill in
    /// the light format (`nomics_io::spill`) and delete it once markdup has
    /// read it, and in-memory runs write no chunk files at all.
    pub save_mapped: bool,
    /// sarek `--save_output_as_bam`. `false` publishes md.cram and recal.cram
    /// (see [`sarek_cram_path`]) instead of md.bam and recal.bam. Whole
    /// genomes encode md.cram from markdup's last pass and recal.cram from
    /// ApplyBQSR's output as those stages write them, and stream md and
    /// recal between the stages as internal spill files ([`stream_path`]),
    /// each removed once its CRAM is complete.
    pub save_output_as_bam: bool,
    /// With `save_output_as_bam` false on whole genomes: stream md.bam and
    /// recal.bam between the stages instead, and keep them (for
    /// `--resume-from`).
    pub keep_bams: bool,
    pub emit_gvcf: bool,
    /// Restart an earlier run in `outdir` at this stage, from its own outputs
    /// (mapped chunk BAMs, md.bam or recal.bam, or their spill files: see
    /// [`find_stream`]), streaming as in whole-genome
    /// mode; for iterating on later stages at whole-genome scale.
    pub resume_from: Option<Resume>,
}

/// Where `--resume-from` restarts a run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum Resume {
    /// Re-align chunks of the earlier run's fastp output (`<outdir>/prep`)
    /// into `<outdir>/mapped`, then stop: `NOMICS_ALIGN_CHUNKS` (comma
    /// list, 1-based) selects chunks, default all.
    Align,
    Markdup,
    Bqsr,
    Hc,
}

struct Stage {
    name: &'static str,
    wall: Duration,
    peak_rss: Option<u64>,
}

enum MarkdupInput {
    /// md.bam already on disk (resuming at BQSR or later).
    Done,
    /// Whole genomes: the chunk BAMs on disk, in chunk order.
    Files(Vec<PathBuf>),
    Raw(Vec<RecordSet>),
    Prepared(
        nomics_markdup::OverlapPlan,
        Vec<nomics_markdup::PreparedSet>,
    ),
}

fn other(message: impl Into<String>) -> io::Error {
    io::Error::other(message.into())
}

fn join<T>(handle: JoinHandle<io::Result<T>>, what: &str) -> io::Result<T> {
    handle
        .join()
        .map_err(|_| other(format!("{what} thread panicked")))?
}

fn join_scoped<T>(
    handle: std::thread::ScopedJoinHandle<'_, io::Result<T>>,
    what: &str,
) -> io::Result<T> {
    handle
        .join()
        .map_err(|_| other(format!("{what} thread panicked")))?
}

fn stage(stages: &mut Vec<Stage>, name: &'static str, started: Instant) {
    stages.push(Stage {
        name,
        wall: started.elapsed(),
        peak_rss: nomics_mem::peak_rss_bytes(),
    });
}

fn reference_paths(args: &Args) -> (PathBuf, PathBuf, PathBuf, Vec<PathBuf>) {
    let fasta = args
        .reference
        .join("Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta");
    let bwa = args
        .reference
        .join("Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64");
    let bundle = args.reference.join("Annotation/GATKBundle");
    let dbsnp = args
        .dbsnp
        .clone()
        .unwrap_or_else(|| bundle.join("dbsnp_146.hg38.vcf.gz"));
    let known_sites = if args.known_sites.is_empty() {
        vec![
            dbsnp.clone(),
            bundle.join("beta/Homo_sapiens_assembly38.known_indels.vcf.gz"),
            bundle.join("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
        ]
    } else {
        args.known_sites.clone()
    };
    (fasta, bwa, dbsnp, known_sites)
}

/// A file read while another thread is still appending to it: at the end of
/// the data it waits for more until `done` is set (after the writer's last
/// flush), then reports end of file. Opens the file once it exists.
struct FollowReader {
    path: PathBuf,
    file: Option<File>,
    done: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

impl FollowReader {
    fn new(path: PathBuf, done: std::sync::Arc<std::sync::atomic::AtomicBool>) -> Self {
        FollowReader {
            path,
            file: None,
            done,
        }
    }
}

impl std::io::Read for FollowReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        use std::sync::atomic::Ordering::Acquire;
        loop {
            // Check before reading: once done is seen, one more read drains
            // everything the writer flushed before setting it.
            let finished = self.done.load(Acquire);
            if self.file.is_none() {
                match File::open(&self.path) {
                    Ok(f) => self.file = Some(f),
                    Err(e) if e.kind() == io::ErrorKind::NotFound && !finished => {
                        std::thread::sleep(Duration::from_millis(5));
                        continue;
                    }
                    Err(e) => return Err(e),
                }
            }
            let n = self.file.as_mut().expect("opened above").read(buf)?;
            if n > 0 || finished {
                return Ok(n);
            }
            std::thread::sleep(Duration::from_millis(5));
        }
    }
}

/// A chunk's FASTQ text shared by the aligner and the intermediate writer;
/// it is freed when the last holder drops it.
#[derive(Clone)]
struct SharedText(std::sync::Arc<Vec<u8>>);

impl AsRef<[u8]> for SharedText {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

fn align_all<R: BufRead + Send>(
    idx: &Index,
    streams: Vec<(R, R)>,
    rg_ids: &[Vec<u8>],
    headers: &[nomics_io::bam::Header],
    threads: usize,
    gpu: bool,
    mut completed: impl FnMut(usize, RecordSet) -> io::Result<()> + Send,
) -> io::Result<()> {
    if std::env::var_os("NOMICS_SAM_TEXT").is_some() {
        for (chunk, (mut r1, mut r2)) in streams.into_iter().enumerate() {
            let mut sam = Vec::new();
            let (pairs, _) = align_stream(
                idx,
                &rg_ids[chunk],
                &mut r1,
                &mut r2,
                threads,
                gpu,
                &mut sam,
            )?;
            let set = sort_sam(headers[chunk].clone(), &sam, threads)?;
            eprintln!(
                "align chunk {:04}: {pairs} pairs, {} records",
                chunk + 1,
                set.len()
            );
            completed(chunk, set)?;
        }
        return Ok(());
    }
    align_chunks_bam(
        idx,
        rg_ids,
        streams,
        headers.to_vec(),
        threads,
        gpu,
        |chunk, pairs, set| {
            eprintln!(
                "align chunk {:04}: {pairs} pairs, {} records",
                chunk + 1,
                set.len()
            );
            completed(chunk, set)
        },
    )?;
    Ok(())
}

/// [`align_all`] for whole genomes: the sorter writes each chunk's mapped BAM
/// itself (`output(chunk)`: its path and BGZF threads) without ever holding
/// the sorted chunk in memory, then `written(chunk, records)` runs. The files
/// are byte-identical to `RecordSet::write` of [`align_all`]'s sets.
#[allow(clippy::too_many_arguments)]
fn align_all_files<R: BufRead + Send>(
    idx: &Index,
    streams: Vec<(R, R)>,
    rg_ids: &[Vec<u8>],
    headers: &[nomics_io::bam::Header],
    threads: usize,
    gpu: bool,
    output: impl Fn(usize) -> (PathBuf, usize) + Send + Sync,
    mut written: impl FnMut(usize, usize) -> io::Result<()> + Send,
) -> io::Result<()> {
    if std::env::var_os("NOMICS_SAM_TEXT").is_some() {
        return align_all(idx, streams, rg_ids, headers, threads, gpu, |chunk, set| {
            let (path, n) = output(chunk);
            set.write(&path, n)?;
            written(chunk, set.len())
        });
    }
    align_chunks_bam_files(
        idx,
        rg_ids,
        streams,
        headers.to_vec(),
        threads,
        gpu,
        &output,
        |chunk, pairs, records| {
            eprintln!("align chunk {:04}: {pairs} pairs, {records} records", chunk + 1);
            written(chunk, records)
        },
    )?;
    Ok(())
}

/// HaplotypeCaller once per interval group, each run with fresh state (its
/// own `java.util.Random`), then MergeVcfs. Groups are contiguous runs of the
/// BED in order, so concatenating their records keeps coordinate order.
/// HaplotypeCaller's input: the recalibrated set in memory, or recal.bam
/// streamed from disk one interval group at a time (whole genomes).
enum HcInput<'a> {
    Set(&'a RecordSet),
    File(&'a Path),
}

/// Call one explicitly supplied interval group, or a list of scatter groups.
pub fn haplotypecaller(
    input: &Path,
    groups: &[PathBuf],
    reference: &Path,
    dbsnp: &Path,
    output: &Path,
    threads: usize,
    emit_gvcf: bool,
) -> io::Result<()> {
    call_variants(
        HcInput::File(input),
        groups,
        reference,
        dbsnp,
        output,
        threads,
        emit_gvcf,
    )
    .map(drop)
}

fn call_variants(
    input: HcInput,
    groups: &[PathBuf],
    fasta_path: &Path,
    dbsnp_path: &Path,
    output: &Path,
    threads: usize,
    emit_gvcf: bool,
) -> io::Result<Option<String>> {
    let stream = match input {
        HcInput::File(path) => Some(nomics_io::bamstream::BamStream::open(path, 1 << 14, 4, 8)?),
        HcInput::Set(_) => None,
    };
    let bam_header = match (&input, &stream) {
        (HcInput::Set(bam), _) => bam.header.clone(),
        (_, Some(stream)) => stream.header.clone(),
        _ => unreachable!(),
    };
    let contigs: Vec<(Vec<u8>, i64)> = bam_header
        .refs
        .iter()
        .map(|r| (r.name.clone(), r.len as i64))
        .collect();
    let lens: Vec<i64> = contigs.iter().map(|c| c.1).collect();
    let names: Vec<String> = contigs
        .iter()
        .map(|c| String::from_utf8_lossy(&c.0).into_owned())
        .collect();
    let sample = String::from_utf8_lossy(&bam_header.text)
        .lines()
        .find(|line| line.starts_with("@RG"))
        .and_then(|line| line.split('\t').find_map(|field| field.strip_prefix("SM:")))
        .ok_or_else(|| other("recalibrated BAM header has no @RG SM"))?
        .to_owned();
    let fasta = Fasta::load(fasta_path, threads)?;
    for (i, name) in names.iter().enumerate() {
        if fasta.index.get(i).map(|entry| &entry.name) != Some(name) {
            return Err(other("FASTA and BAM contig order differ"));
        }
    }
    let dbsnp = Dbsnp::open(dbsnp_path)?;
    let reference = |contig: usize, pos: i64| {
        let base = fasta.contig(contig)[(pos - 1) as usize].to_ascii_uppercase();
        if matches!(base, b'A' | b'C' | b'G' | b'T') {
            base
        } else {
            b'N'
        }
    };
    let contig_bases = |contig: usize| fasta.contig(contig);
    let ctx = Context {
        contig_names: names.clone(),
        contig_lens: lens.clone(),
        reference: &contig_bases,
        dbsnp: Some(&dbsnp),
    };
    let group_intervals = groups
        .iter()
        .map(|group| Intervals::from_bed(group, &contigs))
        .collect::<io::Result<Vec<_>>>()?;
    let mut records = Vec::new();
    let run = if emit_gvcf {
        nomics_hc::caller::run_gvcf
    } else {
        nomics_hc::caller::run
    };
    let header_meta = if emit_gvcf {
        nomics_hc::gvcf::header_meta
    } else {
        header_meta
    };
    let call_group = |g: usize, bam: &RecordSet| {
        let (group, intervals) = (&groups[g], &group_intervals[g]);
        let stream = read_stream(bam, intervals, &lens)?;
        let (traversal, called, timing) =
            run(&ctx, bam, &stream, &intervals.unpadded, &reference, threads)?;
        if groups.len() == 1 {
            for (name, elapsed) in [
                "finalize",
                "assemble",
                "trim",
                "likelihoods",
                "realign",
                "call",
                "annotate",
            ]
            .iter()
            .zip(timing)
            {
                eprintln!(
                    "  {name:<12} {:8.1}s (summed over threads)",
                    elapsed.as_secs_f64()
                );
            }
        }
        eprintln!(
            "haplotypecaller {}: {} regions, {} records",
            group.file_name().unwrap_or_default().to_string_lossy(),
            traversal.regions.len(),
            called.len()
        );
        Ok::<_, io::Error>(called)
    };
    // Whole-genome groups are called two at a time (each with its own state,
    // as before), so one group's tail overlaps the next; NOMICS_HC_OVERLAP=0
    // calls them one at a time.
    let overlap = std::env::var("NOMICS_HC_OVERLAP").map_or(true, |v| v != "0");
    match (input, stream) {
        (HcInput::Set(bam), _) => {
            for g in 0..groups.len() {
                records.extend(call_group(g, bam)?);
            }
        }
        (HcInput::File(_), Some(stream)) => {
            // Groups are collected from disk on request, one ahead of the
            // calling: at most two in memory. Without overlap, group g + 1 is
            // collected while group g is called. With overlap, g and g + 1 are
            // called together and g + 2 is requested once g is done.
            let mut collector = GroupCollector::new(stream, &group_intervals);
            let n = groups.len();
            std::thread::scope(|scope| -> io::Result<()> {
                let (tx, rx) = std::sync::mpsc::sync_channel(1);
                let (want, wanted) = std::sync::mpsc::channel::<()>();
                let reader = scope.spawn(move || {
                    for g in 0..n {
                        if wanted.recv().is_err() {
                            return;
                        }
                        let set = collector.take(g);
                        let failed = set.is_err();
                        if tx.send(set).is_err() || failed {
                            return;
                        }
                    }
                });
                let mut result = Ok(());
                let call_group = &call_group;
                let mut running: Option<std::thread::ScopedJoinHandle<'_, io::Result<_>>> = None;
                let _ = want.send(());
                for g in 0..n {
                    let set = match rx.recv() {
                        Ok(set) => set,
                        Err(_) => Err(other("interval group reader stopped")),
                    };
                    let set = match set {
                        Ok(set) => set,
                        Err(e) => {
                            result = Err(e);
                            break;
                        }
                    };
                    if overlap {
                        let next = scope.spawn(move || call_group(g, &set));
                        if let Some(prev) = running.replace(next) {
                            result = join_scoped(prev, "HaplotypeCaller group")
                                .map(|called| records.extend(called));
                        }
                        if g + 1 < n {
                            let _ = want.send(());
                        }
                    } else {
                        if g + 1 < n {
                            let _ = want.send(());
                        }
                        result = call_group(g, &set).map(|called| records.extend(called));
                    }
                    if result.is_err() {
                        break;
                    }
                }
                if let Some(last) = running.take() {
                    let joined = join_scoped(last, "HaplotypeCaller group");
                    if result.is_ok() {
                        result = joined.map(|called| records.extend(called));
                    }
                }
                drop(want);
                drop(rx);
                reader
                    .join()
                    .map_err(|_| other("interval group reader panicked"))?;
                result
            })?;
        }
        (HcInput::File(_), None) => unreachable!("streamed input opens a stream"),
    }
    let vcf_header = Header {
        meta: header_meta(
            r#"GATKCommandLine=<ID=HaplotypeCaller,CommandLine="nomics hc",Version="4.6.2.0",Date="">"#,
        ),
        contigs: names
            .iter()
            .cloned()
            .zip(lens.iter().map(|&len| len as u64))
            .collect(),
        sample,
    };
    let writer = std::io::BufWriter::new(File::create(output)?);
    if !emit_gvcf {
        return write_vcf_gz(writer, &vcf_header, &records, 8)?
            .flush()
            .map(|()| None);
    }
    // The gVCF text is kept for joint genotyping. BGZF blocks have a fixed
    // payload, so writing it in one piece gives the same bytes.
    let mut text = vcf_header.text().into_bytes();
    for r in &records {
        encode(r, &vcf_header.contigs, &mut text);
    }
    drop(records);
    let mut w = ParallelBgzfWriter::new(writer, 8);
    w.write_all(&text)?;
    w.finish()?.flush()?;
    String::from_utf8(text)
        .map(Some)
        .map_err(|_| other("gVCF text is not UTF-8"))
}

/// Where sarek 3.10.0 publishes `<prefix>.<kind>.cram` (`kind` md or
/// recal): `preprocessing/markduplicates/<id>/` and
/// `preprocessing/recalibrated/<id>/` (conf/modules/markduplicates.config,
/// recalibrate.config), with the `.crai` beside it.
pub fn sarek_cram_path(outdir: &Path, prefix: &str, kind: &str) -> PathBuf {
    let dir = if kind == "md" {
        "markduplicates"
    } else {
        "recalibrated"
    };
    outdir
        .join("preprocessing")
        .join(dir)
        .join(prefix)
        .join(format!("{prefix}.{kind}.cram"))
}

/// Where [`write_sarek_cram`]'s records go: straight into the CRAM, or for a
/// multi-group recal.cram with RG:Z/PG:Z renamed as sarek's CRAM-mode merge
/// names them ([`nomics_bqsr::CramMerge`]).
struct CramFeed<'a, 's> {
    sink: &'a mut nomics_io::cram::RecordSink<'s>,
    merge: Option<&'a nomics_bqsr::CramMerge>,
    buf: Vec<u8>,
}

impl CramFeed<'_, '_> {
    /// Add one raw BAM record (without its `block_size`).
    fn push(&mut self, raw: &[u8]) -> io::Result<()> {
        match self.merge {
            None => self.sink.push(raw),
            Some(m) => {
                self.buf.clear();
                m.translate(raw, &mut self.buf);
                self.sink.push(&self.buf[4..])
            }
        }
    }
}

/// Write records as the CRAM sarek publishes for `kind` (md or recal), plus
/// its `.crai`. `feed` pushes the records in order. `source`: the
/// generation of the stream the records come from ([`StreamId`]; `-` for
/// records in memory), recorded in the completion marker. `groups`: the
/// interval groups ApplyBQSR ran over; for recal with more than one, the
/// records (and `header`) are BQSR's BAM-mode merge, and the CRAM is
/// written as sarek's CRAM-mode `samtools merge` of the groups' CRAMs
/// would be: its ID suffixes (`CramMerge`), and PG:Z after the other tags.
#[allow(clippy::too_many_arguments)]
fn write_sarek_cram(
    header: &nomics_io::bam::Header,
    fasta: &Fasta,
    fasta_path: &Path,
    outdir: &Path,
    prefix: &str,
    kind: &str,
    threads: usize,
    source: &str,
    groups: usize,
    feed: impl FnOnce(&mut CramFeed) -> io::Result<()>,
) -> io::Result<()> {
    let started = Instant::now();
    let path = sarek_cram_path(outdir, prefix, kind);
    let dir = path.parent().unwrap();
    create_dir_all_synced(dir)?;
    // The completion marker goes first (and durably): from here until it is
    // written again, the pair is not complete, whatever else is on disk.
    remove_if_present(&done_path(&path))?;
    File::open(dir)?.sync_all()?;
    let merge = (kind == "recal" && groups > 1)
        .then(|| {
            let fasta_name = fasta_path.file_name().unwrap_or_default().to_string_lossy();
            let cl = format!("samtools view -Ch -T {fasta_name} -o {prefix}.md.cram {prefix}.md.bam");
            let text = String::from_utf8_lossy(&header.text);
            nomics_bqsr::CramMerge::from_recal_header(&text, groups, &cl).map_err(other)
        })
        .transpose()?;
    let merged_header;
    let header = match &merge {
        Some(m) => {
            merged_header = nomics_io::bam::Header {
                text: m.header_text.clone().into_bytes(),
                refs: header.refs.clone(),
            };
            &merged_header
        }
        None => header,
    };
    let profile = match (kind, &merge) {
        ("md", _) => nomics_io::cram::Profile::SarekMd,
        (_, Some(_)) => nomics_io::cram::Profile::SarekRecalMerged,
        _ => nomics_io::cram::Profile::SarekRecal,
    };
    let name = path.file_name().unwrap().to_string_lossy().into_owned();
    let opts = nomics_io::cram::CramOptions {
        threads,
        ..Default::default()
    }
    .with_file_name(&name);
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let (header, opts) = profile.prepare(
        header,
        fasta_path,
        Some((env!("CARGO_PKG_VERSION"), &command_line)),
        opts,
    )?;
    // Both files are written under temporary names and synced, renamed into
    // place, the directory synced, and only then the marker (sizes, CRC32s,
    // a generation id) written, synced and renamed: a pair counts as
    // complete only if the marker matches both files (`cram_complete`).
    let crai = crai_path(&path);
    let tmp = with_suffix(&path, ".tmp");
    let crai_tmp = with_suffix(&crai, ".tmp");
    let file = CrcWriter::new(std::io::BufWriter::with_capacity(
        1 << 22,
        File::create(&tmp)?,
    ));
    let (file, summary) = nomics_io::cram::write_cram(file, &header, fasta, &opts, |sink| {
        feed(&mut CramFeed {
            sink,
            merge: merge.as_ref(),
            buf: Vec::new(),
        })
    })?;
    let (cram_len, cram_crc) = (file.len, file.crc.finish());
    file.inner.into_inner().map_err(|e| e.into_error())?.sync_all()?;
    let crai_bytes = nomics_io::cram::crai_bytes(&summary.index);
    let mut f = File::create(&crai_tmp)?;
    f.write_all(&crai_bytes)?;
    f.sync_all()?;
    drop(f);
    std::fs::rename(&tmp, &path)?;
    std::fs::rename(&crai_tmp, &crai)?;
    File::open(dir)?.sync_all()?;
    let generation = format!(
        "{:x}-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos()),
        std::process::id()
    );
    let marker = DoneMarker {
        generation,
        source: source.to_owned(),
        cram: (cram_len, cram_crc),
        crai: (crai_bytes.len() as u64, nomics_io::crc32::crc32(&crai_bytes)),
    };
    let done = done_path(&path);
    let done_tmp = with_suffix(&done, ".tmp");
    let mut f = File::create(&done_tmp)?;
    f.write_all(marker.text().as_bytes())?;
    f.sync_all()?;
    drop(f);
    std::fs::rename(&done_tmp, &done)?;
    sync_dirs(dir, outdir)?;
    eprintln!(
        "{kind}.cram: {} records, {} bytes in {:.1}s ({})",
        summary.records,
        summary.bytes,
        started.elapsed().as_secs_f64(),
        path.display()
    );
    Ok(())
}

fn crai_path(cram: &Path) -> PathBuf {
    with_suffix(cram, ".crai")
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(suffix);
    PathBuf::from(p)
}

/// fsync `dir` and each parent up to and including `top`, so new entries
/// (and the directories holding them) survive a crash.
fn sync_dirs(dir: &Path, top: &Path) -> io::Result<()> {
    let mut d = Some(dir);
    while let Some(cur) = d {
        File::open(cur)?.sync_all()?;
        if cur == top {
            break;
        }
        d = cur.parent();
    }
    Ok(())
}

/// `<cram>.done`: written last by [`write_sarek_cram`].
fn done_path(cram: &Path) -> PathBuf {
    with_suffix(cram, ".done")
}

fn remove_if_present(path: &Path) -> io::Result<()> {
    match std::fs::remove_file(path) {
        Err(e) if e.kind() != io::ErrorKind::NotFound => Err(e),
        _ => Ok(()),
    }
}

/// `create_dir_all`, then fsync every directory it created and the parent
/// of the topmost one, so the new hierarchy itself survives a crash.
fn create_dir_all_synced(dir: &Path) -> io::Result<()> {
    let mut missing = Vec::new();
    let mut d = Some(dir);
    while let Some(cur) = d {
        if cur.as_os_str().is_empty() || cur.is_dir() {
            break;
        }
        missing.push(cur.to_path_buf());
        d = cur.parent();
    }
    std::fs::create_dir_all(dir)?;
    for m in &missing {
        File::open(m)?.sync_all()?;
    }
    if let Some(top) = missing.last() {
        let parent = top.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
        File::open(parent)?.sync_all()?;
    }
    Ok(())
}

/// A `Write` that counts bytes and CRC32s them on the way through.
struct CrcWriter<W> {
    inner: W,
    crc: nomics_io::crc32::Crc32,
    len: u64,
}

impl<W> CrcWriter<W> {
    fn new(inner: W) -> Self {
        CrcWriter {
            inner,
            crc: nomics_io::crc32::Crc32::new(),
            len: 0,
        }
    }
}

impl<W: Write> Write for CrcWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.crc = self.crc.update(&buf[..n]);
        self.len += n as u64;
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// What `<cram>.done` records about a completed CRAM and `.crai` pair.
#[derive(Debug, PartialEq, Eq)]
struct DoneMarker {
    generation: String,
    /// The generation of the stream the CRAM was made from (`-`: records
    /// held in memory).
    source: String,
    /// (size, CRC32) of each file.
    cram: (u64, u32),
    crai: (u64, u32),
}

impl DoneMarker {
    fn text(&self) -> String {
        format!(
            "nomics-cram-done 2\ngeneration {}\nsource {}\ncram {} {:08x}\ncrai {} {:08x}\n",
            self.generation, self.source, self.cram.0, self.cram.1, self.crai.0, self.crai.1
        )
    }

    fn parse(text: &str) -> Option<DoneMarker> {
        let mut lines = text.lines();
        if lines.next()? != "nomics-cram-done 2" {
            return None;
        }
        let generation = lines.next()?.strip_prefix("generation ")?.to_owned();
        let source = lines.next()?.strip_prefix("source ")?.to_owned();
        let file = |line: Option<&str>, key: &str| -> Option<(u64, u32)> {
            let mut f = line?.strip_prefix(key)?.split(' ');
            let size = f.next()?.parse().ok()?;
            let crc = u32::from_str_radix(f.next()?, 16).ok()?;
            f.next().is_none().then_some((size, crc))
        };
        let cram = file(lines.next(), "cram ")?;
        let crai = file(lines.next(), "crai ")?;
        lines.next().is_none().then_some(DoneMarker {
            generation,
            source,
            cram,
            crai,
        })
    }
}

/// Size and CRC32 of a file's contents.
fn file_crc(path: &Path) -> io::Result<(u64, u32)> {
    let mut f = File::open(path)?;
    let mut buf = vec![0u8; 1 << 22];
    let (mut crc, mut len) = (nomics_io::crc32::Crc32::new(), 0u64);
    loop {
        let n = std::io::Read::read(&mut f, &mut buf)?;
        if n == 0 {
            return Ok((len, crc.finish()));
        }
        crc = crc.update(&buf[..n]);
        len += n as u64;
    }
}

/// The marker of a complete CRAM and `.crai`: `<cram>.done` exists and both
/// files' sizes and CRC32s match it (reads both files).
fn cram_done(outdir: &Path, prefix: &str, kind: &str) -> Option<DoneMarker> {
    let cram = sarek_cram_path(outdir, prefix, kind);
    let marker = std::fs::read_to_string(done_path(&cram))
        .ok()
        .and_then(|t| DoneMarker::parse(&t))?;
    (file_crc(&cram).is_ok_and(|c| c == marker.cram)
        && file_crc(&crai_path(&cram)).is_ok_and(|c| c == marker.crai))
    .then_some(marker)
}

/// Whether the `kind` CRAM is complete ([`cram_done`]) and, with `source`,
/// was made from the stream of that generation.
fn cram_complete(outdir: &Path, prefix: &str, kind: &str, source: Option<&str>) -> bool {
    cram_done(outdir, prefix, kind).is_some_and(|m| source.is_none_or(|s| m.source == s))
}

fn write_set_cram(
    set: &RecordSet,
    fasta: &Fasta,
    fasta_path: &Path,
    args: &Args,
    kind: &str,
    threads: usize,
    groups: usize,
) -> io::Result<()> {
    write_sarek_cram(
        &set.header,
        fasta,
        fasta_path,
        &args.outdir,
        &args.prefix,
        kind,
        threads,
        "-",
        groups,
        |sink| (0..set.len()).try_for_each(|i| sink.push(set.record(i))),
    )
}

/// Whole genomes: re-encode a stream an earlier run's stage wrote (md or
/// recal, as BAM or spill file) as sarek's CRAM, for resumed runs.
#[allow(clippy::too_many_arguments)]
fn write_stream_cram(
    bam: &Path,
    source: &str,
    fasta: &Fasta,
    fasta_path: &Path,
    args: &Args,
    kind: &str,
    threads: usize,
    groups: usize,
) -> io::Result<()> {
    // Our own stream: it must be complete.
    nomics_io::bamstream::check_complete(bam)?;
    let mut stream =
        nomics_io::bamstream::BamStream::open(bam, 10_000, 4 * threads, threads.min(8))?;
    let header = stream.header.clone();
    write_sarek_cram(
        &header,
        fasta,
        fasta_path,
        &args.outdir,
        &args.prefix,
        kind,
        threads,
        source,
        groups,
        |sink| {
            while let Some(batch) = stream.next_batch()? {
                (0..batch.len()).try_for_each(|i| sink.push(batch.record(i)))?;
            }
            Ok(())
        },
    )
}

/// A piece of a BAM's uncompressed stream (header bytes, then
/// `block_size`-prefixed records) handed from the stage writing it to a
/// consumer thread; `None` marks the end of the stream. A channel that
/// closes without `None` (the producing stage failed) is an error to the
/// consumer, so a stream cut short never makes a complete-looking CRAM.
type StreamChunk = Option<Arc<Vec<u8>>>;

/// A consumer of a stage's output stream, run on its own thread: its
/// chunks, and where to return each buffer once read.
type Consumer<'s> =
    Box<dyn FnOnce(Receiver<StreamChunk>, Sender<Arc<Vec<u8>>>) -> io::Result<()> + Send + 's>;

/// Hands each chunk of a stage's output stream to every consumer (the I-06
/// table builder, a CRAM encoder) over bounded channels, one shared buffer
/// for all. Buffers come back through `spare` once every consumer is done
/// with them, and are reused.
struct Tee {
    txs: Vec<SyncSender<StreamChunk>>,
    spare: Receiver<Arc<Vec<u8>>>,
    /// The consumer whose channel closed first: it stopped, and its error
    /// is the cause of the producer's.
    gone: Option<usize>,
}

impl Tee {
    /// Hand over `stream`'s contents (leaving it empty).
    fn send(&mut self, stream: &mut Vec<u8>) -> io::Result<()> {
        if self.txs.is_empty() {
            return Ok(());
        }
        let mut next = Vec::new();
        while let Ok(chunk) = self.spare.try_recv() {
            // The last holder of a chunk to return it gets it back whole.
            if let Ok(buf) = Arc::try_unwrap(chunk) {
                next = buf;
                break;
            }
        }
        next.clear();
        let chunk = Arc::new(std::mem::replace(stream, next));
        self.send_all(Some(chunk))
    }

    /// The end of the stream.
    fn finish(&mut self) -> io::Result<()> {
        self.send_all(None)
    }

    fn send_all(&mut self, chunk: StreamChunk) -> io::Result<()> {
        for (i, tx) in self.txs.iter().enumerate() {
            if tx.send(chunk.clone()).is_err() {
                self.gone = Some(i);
                return Err(other("a consumer of the stage's output stream stopped"));
            }
        }
        Ok(())
    }
}

/// `Read` over a [`StreamChunk`] channel, returning each chunk to `spare`
/// once read.
struct ChunkReader {
    rx: Receiver<StreamChunk>,
    spare: Sender<Arc<Vec<u8>>>,
    cur: Arc<Vec<u8>>,
    at: usize,
    ended: bool,
}

impl ChunkReader {
    fn new(rx: Receiver<StreamChunk>, spare: Sender<Arc<Vec<u8>>>) -> Self {
        ChunkReader {
            rx,
            spare,
            cur: Default::default(),
            at: 0,
            ended: false,
        }
    }
}

impl io::Read for ChunkReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.at == self.cur.len() {
            if self.ended || buf.is_empty() {
                return Ok(0);
            }
            match self.rx.recv() {
                Ok(Some(chunk)) => {
                    let done = std::mem::replace(&mut self.cur, chunk);
                    let _ = self.spare.send(done);
                    self.at = 0;
                }
                Ok(None) => self.ended = true,
                Err(_) => {
                    return Err(other(
                        "the stage writing this stream stopped before its end",
                    ));
                }
            }
        }
        let n = buf.len().min(self.cur.len() - self.at);
        buf[..n].copy_from_slice(&self.cur[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

/// Run `produce` (a stage writing a stream, handed a [`Tee`]) with
/// `consumers` reading that stream on their own threads, each named for
/// errors. The producer's own error wins unless it stopped because a
/// consumer had: that consumer's error is then the cause.
fn tee_stage<'s, T>(
    scope: &'s std::thread::Scope<'s, '_>,
    consumers: Vec<(&'static str, Consumer<'s>)>,
    produce: impl FnOnce(&mut Tee) -> io::Result<T>,
) -> io::Result<T> {
    let (spare_tx, spare_rx) = std::sync::mpsc::channel();
    let mut txs = Vec::new();
    let mut handles = Vec::new();
    for (what, consume) in consumers {
        let (tx, rx) = std::sync::mpsc::sync_channel::<StreamChunk>(2);
        let spare = spare_tx.clone();
        txs.push(tx);
        handles.push((what, scope.spawn(move || consume(rx, spare))));
    }
    drop(spare_tx);
    let mut tee = Tee {
        txs,
        spare: spare_rx,
        gone: None,
    };
    let produced = produce(&mut tee).and_then(|v| tee.finish().map(|()| v));
    let gone = tee.gone;
    // Closes every channel: a consumer not yet at the end now fails.
    drop(tee);
    let mut cause = None;
    let mut first = None;
    for (i, (what, handle)) in handles.into_iter().enumerate() {
        if let Err(e) = join_scoped(handle, what) {
            let e = other(format!("{what}: {e}"));
            if gone == Some(i) {
                cause = Some(e);
            } else if first.is_none() {
                first = Some(e);
            }
        }
    }
    match (produced, cause, first) {
        (Err(_), Some(e), _) => Err(e),
        (Err(e), None, _) => Err(e),
        (Ok(_), Some(e), _) | (Ok(_), None, Some(e)) => Err(e),
        (Ok(v), None, None) => Ok(v),
    }
}

/// The FASTA the CRAM encoders need, loaded on its own thread while the
/// stages start; the first encoder to need it waits for the load.
struct LazyFasta {
    loader: std::sync::Mutex<Option<JoinHandle<io::Result<Fasta>>>>,
    fasta: std::sync::OnceLock<Fasta>,
}

impl LazyFasta {
    fn load(path: &Path, threads: usize) -> LazyFasta {
        let path = path.to_path_buf();
        LazyFasta {
            loader: std::sync::Mutex::new(Some(std::thread::spawn(move || {
                Fasta::load(&path, threads)
            }))),
            fasta: std::sync::OnceLock::new(),
        }
    }

    fn get(&self) -> io::Result<&Fasta> {
        let mut loader = self
            .loader
            .lock()
            .map_err(|_| other("FASTA loader lock poisoned"))?;
        if let Some(fasta) = self.fasta.get() {
            return Ok(fasta);
        }
        let handle = loader
            .take()
            .ok_or_else(|| other("loading the FASTA for the CRAM encoder failed"))?;
        let fasta = join(handle, "FASTA loader")?;
        Ok(self.fasta.get_or_init(|| fasta))
    }
}

/// A [`Consumer`] encoding sarek's `kind` CRAM from the stream a stage is
/// writing, as it writes it: md.cram from markdup's last pass, recal.cram
/// from ApplyBQSR's output, instead of reading md.bam / recal.bam back.
/// `source`: the generation the stage's stream will be recorded with.
fn cram_consumer<'s>(
    fasta: &'s LazyFasta,
    fasta_path: &'s Path,
    args: &'s Args,
    kind: &'static str,
    threads: usize,
    source: &'s str,
    groups: usize,
) -> (&'static str, Consumer<'s>) {
    let what = if kind == "md" {
        "md.cram writer"
    } else {
        "recal.cram writer"
    };
    let consume: Consumer<'s> = Box::new(move |rx, spare| {
        let mut r = ChunkReader::new(rx, spare);
        let header = nomics_io::bam::Header::read(&mut r)?;
        let fasta = fasta.get()?;
        let mut rec = Vec::new();
        write_sarek_cram(
            &header,
            fasta,
            fasta_path,
            &args.outdir,
            &args.prefix,
            kind,
            threads,
            source,
            groups,
            |sink| {
                while nomics_io::bam::read_raw_record(&mut r, &mut rec)? {
                    sink.push(&rec)?;
                }
                Ok(())
            },
        )
    });
    (what, consume)
}

fn write_timings(outdir: &Path, stages: &[Stage], total: Duration) -> io::Result<()> {
    let mut text = String::from("stage\twall_s\tpeak_rss_gb\n");
    for row in stages {
        let rss = row
            .peak_rss
            .map(|bytes| format!("{:.3}", bytes as f64 / (1u64 << 30) as f64))
            .unwrap_or_default();
        text.push_str(&format!(
            "{}\t{:.3}\t{}\n",
            row.name,
            row.wall.as_secs_f64(),
            rss
        ));
    }
    let peak = nomics_mem::peak_rss_bytes()
        .map(|bytes| format!("{:.3}", bytes as f64 / (1u64 << 30) as f64))
        .unwrap_or_default();
    text.push_str(&format!("TOTAL\t{:.3}\t{peak}\n", total.as_secs_f64()));
    std::fs::write(outdir.join("stages.tsv"), &text)?;
    eprintln!("stage\twall_s\tpeak_rss_gb");
    for line in text.lines().skip(1) {
        eprintln!("{line}");
    }
    Ok(())
}

/// One sequencing lane (a samplesheet row) and sarek's names for it.
struct Lane {
    name: String,
    fastq1: PathBuf,
    fastq2: PathBuf,
    /// `<flowcell>.<sample>.<lane>`, or `<sample>.<lane>` when the first read
    /// name carries no flowcell (sarek `addReadgroupToMeta`).
    rg_id: String,
    /// fastp's output prefix, `<sample>-<lane>`.
    prep_prefix: String,
}

/// sarek's `flowcellLaneFromFastq`: the first line split on ':' (comment
/// included); five fields give the instrument, seven or more the flowcell.
fn flowcell(fastq: &Path) -> io::Result<Option<String>> {
    if std::fs::metadata(fastq)?.len() == 0 {
        return Ok(None);
    }
    let head = nomics_io::gzip::head(fastq, 1 << 16)?;
    let line = head.split(|&b| b == b'\n').next().unwrap_or_default();
    let line = String::from_utf8_lossy(line);
    let line = line.trim_end_matches('\r');
    // Java's split drops trailing empty strings.
    let mut fields: Vec<&str> = line.split(':').collect();
    while fields.last() == Some(&"") {
        fields.pop();
    }
    Ok(match fields.len() {
        5 => Some(fields[0].strip_prefix('@').unwrap_or(fields[0]).to_owned()),
        n if n >= 7 => Some(fields[2].to_owned()),
        _ => None,
    })
}

fn lanes(args: &Args) -> io::Result<Vec<Lane>> {
    if args.fastq1.len() != args.fastq2.len() || args.fastq1.is_empty() {
        return Err(other("give one --fastq1 and one --fastq2 per lane"));
    }
    if !args.lanes.is_empty() && args.lanes.len() != args.fastq1.len() {
        return Err(other("give one --lane per FASTQ pair"));
    }
    (0..args.fastq1.len())
        .map(|i| {
            let name = args
                .lanes
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("L{}", i + 1));
            let (fastq1, fastq2) = (args.fastq1[i].clone(), args.fastq2[i].clone());
            let fc = flowcell(&fastq1)?;
            if fc.is_some() && fc != flowcell(&fastq2)? {
                return Err(other(format!(
                    "flowcell differs between {} and {}",
                    fastq1.display(),
                    fastq2.display()
                )));
            }
            let rg_id = match fc {
                Some(fc) => format!("{fc}.{}.{name}", args.prefix),
                None => format!("{}.{name}", args.prefix),
            };
            let prep_prefix = format!("{}-{name}", args.prefix);
            Ok(Lane {
                name,
                fastq1,
                fastq2,
                rg_id,
                prep_prefix,
            })
        })
        .collect()
}

/// A mapped chunk BAM's file name: `<sample>.<chunk>.bam`, with the lane
/// (`<sample>-<lane>.<chunk>.bam`) when there are several; with `spill`,
/// the light spill file's extension instead of `bam`.
fn mapped_name(args: &Args, lanes: &[Lane], lane: usize, chunk: usize, spill: bool) -> String {
    let ext = if spill {
        nomics_io::spill::EXTENSION
    } else {
        "bam"
    };
    if lanes.len() == 1 {
        format!("{}.{:04}.{ext}", args.prefix, chunk + 1)
    } else {
        format!("{}.{:04}.{ext}", lanes[lane].prep_prefix, chunk + 1)
    }
}

/// A mapped chunk BAM's header: @SQ lines, the lane's read group, bwa's @PG.
fn chunk_header(
    args: &Args,
    lanes: &[Lane],
    lane: usize,
    idx: &Index,
    fasta_path: &Path,
    threads: usize,
    chunk: usize,
) -> nomics_io::bam::Header {
    let lane_flag = if lanes.len() > 1 {
        format!(" --lane {}", lanes[lane].name)
    } else {
        String::new()
    };
    let lane = &lanes[lane];
    let rg = format!(
        "@RG\tID:{}\tPU:{}\tSM:{}\tLB:{}\tDS:{}\tPL:ILLUMINA",
        lane.rg_id,
        lane.name,
        args.sample,
        args.prefix,
        args.read_group_ds
            .clone()
            .unwrap_or_else(|| fasta_path.display().to_string())
    );
    let pg_cl = format!(
        "nomics run --reference {}{lane_flag} --chunk {:04} --threads {} --dense-sa-cache {}",
        args.reference.display(),
        chunk + 1,
        threads,
        args.dense_sa_cache.display()
    );
    header(idx, &rg, &pg_cl)
}

/// With GPU seeding on, warm the GPU streams before anything large is
/// mapped ([`nomics_align::mem::gpu_warm_up`]). A failure is reported
/// loudly; alignment then falls back to the CPU per call (same output).
fn warm_gpu(threads: usize, whole_genome: bool) {
    if std::env::var_os("NOMICS_GPU_SEED").is_none() {
        return;
    }
    let started = Instant::now();
    match nomics_align::mem::gpu_warm_up(threads, whole_genome) {
        Ok(n) => eprintln!(
            "gpu: warmed {n} streams in {:.2}s",
            started.elapsed().as_secs_f64()
        ),
        Err(e) => {
            nomics_gpu::diag_gpu_error("GPU warm-up", &e);
            eprintln!(
                "WARNING: GPU warm-up failed ({e:?}); GPU calls will likely fall back to the CPU (same output, slower). Free or defragment memory and rerun for normal speed."
            )
        }
    }
}

/// The end-of-alignment line for GPU calls that fell back to the CPU.
fn report_gpu_fallbacks() {
    if let Some(line) = nomics_align::mem::gpu_fallback_summary() {
        eprintln!("WARNING: {line}");
    }
    if let Some(line) = nomics_mem::diag::summary() {
        eprintln!("{line}");
    }
}

/// `--resume-from align`: re-align chunks of `<outdir>/prep` into
/// `<outdir>/mapped`, then stop.
fn align_only(args: &Args, threads: usize, bwa_path: &Path, fasta_path: &Path) -> io::Result<()> {
    let total_started = Instant::now();
    let mut stages = Vec::new();
    let chunks: Vec<usize> = match std::env::var("NOMICS_ALIGN_CHUNKS") {
        Ok(list) => list
            .split(',')
            .map(|c| match c.trim().parse::<usize>() {
                Ok(k) if (1..=PREP_CHUNKS).contains(&k) => Ok(k - 1),
                _ => Err(other(format!("NOMICS_ALIGN_CHUNKS: bad chunk {c:?}"))),
            })
            .collect::<io::Result<_>>()?,
        Err(_) => (0..PREP_CHUNKS).collect(),
    };
    // Chunk BAMs are streamed to files as in whole-genome runs.
    warm_gpu(threads, true);
    let started = Instant::now();
    let mut idx = Index::load(bwa_path)?;
    idx.attach_dense_sa_with_interval(&args.dense_sa_cache, args.dense_sa_interval, threads)?;
    if let Err(e) = idx.lock_memory() {
        eprintln!("index: mlock failed ({e})");
    }
    let (index_bytes, huge) = idx.memory();
    eprintln!(
        "index: {:.1} GB, {} on 2 MiB pages",
        index_bytes as f64 / 1e9,
        huge.map_or("unknown".to_owned(), |h| format!(
            "{:.1} GB",
            h as f64 / 1e9
        ))
    );
    // This mode exists for timing; GPU seeding on 4 KiB pages is several
    // times slower, so a run there would only mislead.
    if let Some(h) = huge
        && (h as f64) < 0.9 * index_bytes as f64
        && std::env::var_os("NOMICS_ALLOW_SMALL_PAGES").is_none()
    {
        return Err(other(format!(
            "index only {:.0}% on 2 MiB pages (NOMICS_ALLOW_SMALL_PAGES=1 to run anyway)",
            100.0 * h as f64 / index_bytes as f64
        )));
    }
    stage(&mut stages, "index", started);
    let started = Instant::now();
    let lanes = lanes(args)?;
    let prep_dir = args.outdir.join("prep");
    let gpu = std::env::var_os("NOMICS_GPU_SEED").is_some();
    let mapped_dir = args.outdir.join("mapped");
    let timing = std::env::var_os("NOMICS_TIMING").is_some();
    nomics_align::extend::ext_stats::ON.store(timing, std::sync::atomic::Ordering::Relaxed);
    // Records written (reads plus a few supplementary alignments).
    let mut records = 0u64;
    // One chunk at a time: each chunk's gzip text is read whole.
    for (l, lane) in lanes.iter().enumerate() {
        for &c in &chunks {
            let prep_prefix = &lane.prep_prefix;
            let path =
                |r: u8| prep_dir.join(format!("{:04}.{prep_prefix}_R{r}.fastp.fastq.gz", c + 1));
            let stream = (
                BufReader::with_capacity(4 << 20, nomics_io::gzip::GzReader::open(path(1))?),
                BufReader::with_capacity(4 << 20, nomics_io::gzip::GzReader::open(path(2))?),
            );
            let header = chunk_header(args, &lanes, l, &idx, fasta_path, threads, c);
            let chunk_started = Instant::now();
            align_all_files(
                &idx,
                vec![stream],
                &[lane.rg_id.as_bytes().to_vec()],
                &[header],
                threads,
                gpu,
                |_| (mapped_dir.join(mapped_name(args, &lanes, l, c, false)), threads),
                |_, n| {
                    records += n as u64;
                    Ok(())
                },
            )?;
            eprintln!(
                "align-only {} chunk {:04}: {:.1}s",
                lane.name,
                c + 1,
                chunk_started.elapsed().as_secs_f64()
            );
        }
    }
    stage(&mut stages, "align", started);
    report_gpu_fallbacks();
    if timing {
        eprintln!("{}", nomics_align::extend::ext_stats::report(records));
    }
    write_timings(&args.outdir, &stages, total_started.elapsed())
}

pub fn run(args: Args) -> io::Result<()> {
    let threads = args.threads.max(1);
    let total_started = Instant::now();
    let mut stages = Vec::new();
    create_dir_all_synced(&args.outdir)?;
    // Resolve (and validate) the BAM compression mode before any work.
    eprintln!(
        "BAM compression: {:?} ({})",
        nomics_io::bgzf::BamCompression::from_env(),
        nomics_io::bgzf::BamCompression::ENV
    );
    eprintln!(
        "BQSR tables built during markdup (whole genomes): {} (NOMICS_BQSR_IN_MARKDUP)",
        bqsr_in_markdup()
    );
    eprintln!(
        "chunk BAMs: {}",
        if args.save_mapped {
            "published in mapped/ (--save-mapped true)"
        } else {
            "not published; whole genomes spill them in the light format (--save-mapped false)"
        }
    );
    let prep_dir = args.outdir.join("prep");
    let mapped_dir = args.outdir.join("mapped");
    std::fs::create_dir_all(&mapped_dir)?;
    let (fasta_path, bwa_path, dbsnp_path, known_sites) = reference_paths(&args);
    if args.resume_from == Some(Resume::Align) {
        return align_only(&args, threads, &bwa_path, &fasta_path);
    }
    if let Some(from) = args.resume_from {
        let input = if from == Resume::Markdup {
            // `<prefix>.<chunk>.bam`, or `<prefix>-<lane>.<chunk>.bam`; or
            // the same names as spill files, left by a run with
            // `--save-mapped false` that stopped before markdup finished.
            let prefix = &args.prefix;
            let spill_ext = format!(".{}", nomics_io::spill::EXTENSION);
            let mut inputs: Vec<PathBuf> = std::fs::read_dir(&mapped_dir)?
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| {
                    p.file_name().and_then(|n| n.to_str()).is_some_and(|n| {
                        n.strip_prefix(prefix.as_str())
                            .is_some_and(|r| r.starts_with('.') || r.starts_with('-'))
                            && (n.ends_with(".bam") || n.ends_with(&spill_ext))
                    })
                })
                .collect();
            inputs.sort();
            if inputs.is_empty() {
                return Err(other(format!(
                    "no {prefix}*.bam chunk BAMs (or {prefix}*{spill_ext} spill files) in {}",
                    mapped_dir.display()
                )));
            }
            let spilled = inputs
                .iter()
                .filter(|p| nomics_io::spill::is_spill_path(p))
                .count();
            if spilled != 0 && spilled != inputs.len() {
                return Err(other(format!(
                    "{} holds both chunk BAMs and spill files; remove one set",
                    mapped_dir.display()
                )));
            }
            // A spill set is only taken whole: its manifest, written when
            // alignment finished, must name exactly the files present (a
            // cleanup or run cut short leaves a partial set). Chunk BAMs
            // are taken as found.
            if spilled != 0 {
                nomics_io::spill::check_manifest(&mapped_dir, &inputs)?;
            }
            // Our own chunk files end as complete ones do (a chunk BAM cut
            // at a block boundary would otherwise read as a shorter one).
            for input in &inputs {
                nomics_io::bamstream::check_complete(input)?;
            }
            MarkdupInput::Files(inputs)
        } else {
            MarkdupInput::Done
        };
        return finish_run(
            &args,
            threads,
            stages,
            total_started,
            input,
            from,
            (&fasta_path, &dbsnp_path, known_sites),
        );
    }

    // Each lane is its own fastp run, split into PREP_CHUNKS chunks; chunk
    // `g` overall is lane `g / PREP_CHUNKS`, chunk `g % PREP_CHUNKS` of it.
    let lanes = lanes(&args)?;
    let gz_bytes = lanes
        .iter()
        .flat_map(|l| [&l.fastq1, &l.fastq2])
        .map(|p| std::fs::metadata(p).map(|m| m.len()))
        .sum::<io::Result<u64>>()?;
    let spill_dir =
        (args.spill_prep || gz_bytes > SPILL_PREP_ABOVE).then(|| args.outdir.join("spill"));
    // Spilled runs align as whole genomes (chunk BAMs streamed to files),
    // which sets how many batches are in flight and so the engines warmed.
    warm_gpu(threads, spill_dir.is_some());
    // The index (BWT, pac, SA, dense SA) does not depend on prep, so load it
    // while prep runs; loading is mostly page-cache copies on one thread.
    let index_loader = {
        let (bwa_path, cache, interval) = (
            bwa_path.clone(),
            args.dense_sa_cache.clone(),
            args.dense_sa_interval,
        );
        std::thread::spawn(move || -> io::Result<(Index, f64)> {
            let started = Instant::now();
            let mut idx = Index::load(&bwa_path)?;
            idx.attach_dense_sa_with_interval(&cache, interval, threads)?;
            // Keep the index resident and on 2 MiB pages while prep's spill
            // fills the page cache (the GPU's reads do not mark it in use).
            if let Err(e) = idx.lock_memory() {
                eprintln!("index: mlock failed ({e}); it may be swapped out under memory pressure");
            }
            Ok((idx, started.elapsed().as_secs_f64()))
        })
    };

    let mut index_loader = Some(index_loader);
    let mut early_index = None;
    let started = Instant::now();
    // Whole genomes: prep spills on its own thread while alignment follows
    // the chunk files as they grow (NOMICS_PREP_OVERLAP=0: prep first).
    let overlap_prep = std::env::var("NOMICS_PREP_OVERLAP").map_or(true, |v| v != "0");
    let prep_done = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
    let mut prep_thread = None;
    let chunks = match &spill_dir {
        Some(dir) => {
            let jobs: Vec<(String, PathBuf, PathBuf, PathBuf)> = lanes
                .iter()
                .map(|l| {
                    (
                        l.name.clone(),
                        l.fastq1.clone(),
                        l.fastq2.clone(),
                        dir.join(&l.name),
                    )
                })
                .collect();
            // A spill left by an aborted run would be followed before prep
            // truncates it: start clean.
            if dir.exists() {
                std::fs::remove_dir_all(dir)?;
            }
            for (_, _, _, lane_dir) in &jobs {
                std::fs::create_dir_all(lane_dir)?;
            }
            if overlap_prep {
                // Load the index before prep's spill churns the page cache, so
                // it lands on 2 MiB pages (GPU seeding is slow on 4 KiB pages).
                let loader = index_loader.take().expect("index loader");
                early_index = Some(join(loader, "index loader")?);
            }
            let done = prep_done.clone();
            let shown = dir.display().to_string();
            let spill = std::thread::spawn(move || -> io::Result<Duration> {
                let t = Instant::now();
                let result = jobs.iter().try_for_each(|(name, fq1, fq2, lane_dir)| {
                    let stats = nomics_prep::run_to_dir(fq1, fq2, PREP_CHUNKS, lane_dir)?;
                    eprintln!("prep {name} (spilled to {shown}): {stats:?}");
                    Ok(())
                });
                done.store(true, std::sync::atomic::Ordering::Release);
                result.map(|()| t.elapsed())
            });
            if overlap_prep {
                prep_thread = Some(spill);
            } else {
                join(spill, "prep")?;
            }
            None
        }
        None => {
            let mut all = nomics_prep::Chunks {
                r1: Vec::new(),
                r2: Vec::new(),
            };
            for lane in &lanes {
                let (chunks, stats) = nomics_prep::run(&lane.fastq1, &lane.fastq2, PREP_CHUNKS)?;
                eprintln!("prep {}: {stats:?}", lane.name);
                all.r1.extend(chunks.r1);
                all.r2.extend(chunks.r2);
            }
            Some(all)
        }
    };
    if prep_thread.is_none() {
        stage(&mut stages, "prep", started);
    }

    let started = Instant::now();
    let (idx, load_s) = match early_index {
        Some(loaded) => loaded,
        None => join(index_loader.take().expect("index loader"), "index loader")?,
    };
    // GPU seeding reads the index in place and is only fast on 2 MiB pages;
    // a fragmented page cache can leave it on 4 KiB pages.
    let (index_bytes, huge) = idx.memory();
    eprintln!(
        "index: {:.1} GB, {} on 2 MiB pages",
        index_bytes as f64 / 1e9,
        huge.map_or("unknown".to_owned(), |h| format!(
            "{:.1} GB",
            h as f64 / 1e9
        ))
    );
    if let Some(h) = huge
        && (h as f64) < 0.9 * index_bytes as f64
    {
        eprintln!(
            "WARNING: index only {:.0}% on 2 MiB pages (fragmented memory): GPU seeding will be slower and timings are not comparable",
            100.0 * h as f64 / index_bytes as f64
        );
    }
    eprintln!(
        "align setup: index loaded in {load_s:.1}s alongside prep; waited {:.1}s after prep",
        started.elapsed().as_secs_f64()
    );
    let gpu = std::env::var_os("NOMICS_GPU_SEED").is_some();
    let count = PREP_CHUNKS * lanes.len();
    let align_headers: Vec<_> = (0..count)
        .map(|g| {
            chunk_header(
                &args,
                &lanes,
                g / PREP_CHUNKS,
                &idx,
                &fasta_path,
                threads,
                g % PREP_CHUNKS,
            )
        })
        .collect();
    let rg_ids: Vec<Vec<u8>> = (0..count)
        .map(|g| lanes[g / PREP_CHUNKS].rg_id.as_bytes().to_vec())
        .collect();
    let md_output = args.outdir.join(format!("{}.md.bam", args.prefix));
    let metrics = args.outdir.join(format!("{}.md.metrics", args.prefix));
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let md_opts = nomics_markdup::Options {
        inputs: Vec::new(),
        output: md_output.clone(),
        metrics,
        pg_version: format!("Version:{}", env!("CARGO_PKG_VERSION")),
        pg_command_line: command_line.clone(),
        threads,
    };
    let overlap = std::env::var("NOMICS_MARKDUP_OVERLAP").map_or(true, |v| v != "0");
    // Spilled (whole-genome) runs stream markdup from the chunk BAMs instead
    // of holding every chunk in memory.
    let disk = spill_dir.is_some();
    // Unpublished whole-genome chunks are only markdup's spill.
    let spill_mapped = disk && !args.save_mapped;
    let save_mapped = args.save_mapped;
    if spill_mapped {
        // An earlier run's manifest describes an earlier set.
        match std::fs::remove_file(mapped_dir.join(nomics_io::spill::MANIFEST)) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => return Err(e),
            _ => {}
        }
    }
    // Whole genomes: each chunk's BAM is compressed while later chunks align
    // (the sorter streams it out, see `align_all_files`), on one thread per
    // batch in flight: 2 threads write a 30x chunk in ~190 s, within the
    // ~230 s a chunk aligns at 2 batches, but not the ~200 s at 3 (fused),
    // where the writer fell further behind each chunk and the chunks waiting
    // for it grew the process by ~1.5 GB per chunk; on 1 thread it took about
    // as long as the alignment at 2. Same CPU work, same bytes either way.
    let overlap_threads = std::env::var("NOMICS_MARKDUP_OVERLAP_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(if disk {
            nomics_align::mem::batches_in_flight(threads, true).max(2)
        } else {
            1
        })
        .clamp(1, threads);
    let overlap_plan = (overlap && !disk)
        .then(|| nomics_markdup::OverlapPlan::new(align_headers.clone(), &md_opts))
        .transpose()?;
    let (bam_tx, bam_rx) = std::sync::mpsc::channel::<(usize, Option<RecordSet>, PathBuf)>();
    let aligning = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
    let writer_aligning = aligning.clone();
    // Sorted chunks sent to the writer and not yet taken (each a whole
    // chunk in memory, ~21 GB at 30x).
    let queued = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let writer_queued = queued.clone();
    let bam_writer = std::thread::spawn(move || -> io::Result<MarkdupInput> {
        let mut raw: Vec<Option<RecordSet>> = (0..count).map(|_| None).collect();
        let mut files: Vec<Option<PathBuf>> = (0..count).map(|_| None).collect();
        let mut prepared: Vec<Option<nomics_markdup::PreparedSet>> =
            (0..count).map(|_| None).collect();
        for (chunk, set, output) in bam_rx {
            let behind = writer_queued
                .fetch_sub(1, std::sync::atomic::Ordering::AcqRel)
                .saturating_sub(1);
            // While alignment saturates the CPU budget, the writer and markdup
            // preparation stay within the overlap budget (one thread by
            // default), widened by one share per chunk already waiting behind
            // this one so a backlog of whole chunks cannot build up; the
            // chunks still queued when alignment ends get all.
            let budget = || {
                if writer_aligning.load(std::sync::atomic::Ordering::Relaxed) {
                    (overlap_threads * (1 + behind)).min(threads)
                } else {
                    threads
                }
            };
            // `None`: the sorter wrote this chunk's BAM itself (whole genomes).
            let Some(set) = set else {
                files[chunk] = Some(output);
                continue;
            };
            // In memory, an unpublished chunk is never written.
            if save_mapped {
                set.write(&output, budget())?;
            }
            if disk {
                files[chunk] = Some(output);
            } else if let Some(plan) = &overlap_plan {
                let started = Instant::now();
                let n = budget();
                prepared[chunk] = Some(plan.prepare_set(chunk, set, n)?);
                eprintln!(
                    "markdup overlap prepare chunk {:04} ({n} threads): {:.1}s",
                    chunk + 1,
                    started.elapsed().as_secs_f64()
                );
            } else {
                raw[chunk] = Some(set);
            }
        }
        if disk {
            Ok(MarkdupInput::Files(
                files
                    .into_iter()
                    .map(|p| p.ok_or_else(|| other("missing alignment chunk BAM")))
                    .collect::<io::Result<_>>()?,
            ))
        } else if let Some(plan) = overlap_plan {
            Ok(MarkdupInput::Prepared(
                plan,
                prepared
                    .into_iter()
                    .map(|p| p.ok_or_else(|| other("missing prepared alignment chunk")))
                    .collect::<io::Result<_>>()?,
            ))
        } else {
            Ok(MarkdupInput::Raw(
                raw.into_iter()
                    .map(|s| s.ok_or_else(|| other("missing alignment chunk")))
                    .collect::<io::Result<_>>()?,
            ))
        }
    });
    let mut completed = 0;
    let mut aligned = |chunk: usize, set: Option<RecordSet>| {
        // The last chunk reaches the writer before alignment returns; widen
        // the writer's budget before it starts on it.
        completed += 1;
        if completed == count {
            aligning.store(false, std::sync::atomic::Ordering::Relaxed);
        }
        let output = mapped_dir.join(mapped_name(
            &args,
            &lanes,
            chunk / PREP_CHUNKS,
            chunk % PREP_CHUNKS,
            spill_mapped,
        ));
        queued.fetch_add(1, std::sync::atomic::Ordering::AcqRel);
        bam_tx
            .send((chunk, set, output))
            .map_err(|_| other("mapped BAM writer stopped"))
    };
    // Chunk FASTQ text comes from memory, or for whole genomes from prep's
    // spill on disk; either way the intermediate writer reads it alongside
    // the aligner.
    std::thread::scope(|scope| -> io::Result<()> {
        let (prep_dir, lanes) = (&prep_dir, &lanes);
        match (chunks, &spill_dir) {
            (Some(chunks), _) => {
                let shared: Vec<(SharedText, SharedText)> = chunks
                    .r1
                    .into_iter()
                    .zip(chunks.r2)
                    .map(|(r1, r2)| {
                        (
                            SharedText(std::sync::Arc::new(r1)),
                            SharedText(std::sync::Arc::new(r2)),
                        )
                    })
                    .collect();
                let streams = shared
                    .iter()
                    .map(|(r1, r2)| (Cursor::new(r1.clone()), Cursor::new(r2.clone())))
                    .collect();
                let writer = args.write_intermediates.then(|| {
                    scope.spawn(move || -> io::Result<()> {
                        for (g, (r1, r2)) in shared.into_iter().enumerate() {
                            nomics_prep::write_chunk(
                                r1.as_ref(),
                                r2.as_ref(),
                                prep_dir,
                                &lanes[g / PREP_CHUNKS].prep_prefix,
                                g % PREP_CHUNKS + 1,
                                threads,
                            )?;
                        }
                        Ok(())
                    })
                });
                let aligned = align_all(
                    &idx,
                    streams,
                    &rg_ids,
                    &align_headers,
                    threads,
                    gpu,
                    |chunk, set| aligned(chunk, Some(set)),
                );
                if let Some(writer) = writer {
                    join_scoped(writer, "FASTQ writer")?;
                }
                aligned
            }
            (None, Some(dir)) => {
                let prep_done = &prep_done;
                let streams = (0..count)
                    .map(|g| {
                        let lane_dir = dir.join(&lanes[g / PREP_CHUNKS].name);
                        let [p1, p2] = nomics_prep::spill_paths(&lane_dir, g % PREP_CHUNKS);
                        (
                            BufReader::with_capacity(
                                4 << 20,
                                FollowReader::new(p1, prep_done.clone()),
                            ),
                            BufReader::with_capacity(
                                4 << 20,
                                FollowReader::new(p2, prep_done.clone()),
                            ),
                        )
                    })
                    .collect::<Vec<_>>();
                let writer = args.write_intermediates.then(|| {
                    scope.spawn(move || -> io::Result<()> {
                        // Compress the spill only once prep has written all of it.
                        while !prep_done.load(std::sync::atomic::Ordering::Acquire) {
                            std::thread::sleep(Duration::from_millis(200));
                        }
                        for g in 0..count {
                            let lane = &lanes[g / PREP_CHUNKS];
                            nomics_prep::compress_spill(
                                &dir.join(&lane.name),
                                g % PREP_CHUNKS,
                                prep_dir,
                                &lane.prep_prefix,
                                threads,
                            )?;
                        }
                        Ok(())
                    })
                });
                // The sorter writes each chunk's BAM (2 threads by default
                // while later chunks align, all of them for the last), so a
                // whole sorted chunk (~21 GB at 30x) is never held for the writer.
                let aligned = align_all_files(
                    &idx,
                    streams,
                    &rg_ids,
                    &align_headers,
                    threads,
                    gpu,
                    |chunk| {
                        let path = mapped_dir.join(mapped_name(
                            &args,
                            lanes,
                            chunk / PREP_CHUNKS,
                            chunk % PREP_CHUNKS,
                            spill_mapped,
                        ));
                        let n = if chunk + 1 == count {
                            threads
                        } else {
                            overlap_threads
                        };
                        (path, n)
                    },
                    |chunk, _| aligned(chunk, None),
                );
                if let Some(writer) = writer {
                    join_scoped(writer, "FASTQ writer")?;
                }
                aligned?;
                std::fs::remove_dir_all(dir)
            }
            (None, None) => unreachable!("prep keeps its chunks or spills them"),
        }
    })?;
    if let Some(spill) = prep_thread {
        let wall = join(spill, "prep")?;
        // Prep overlapped alignment: its row is its own wall time, inserted
        // before the align row.
        let at = stages
            .iter()
            .position(|s| s.name == "align")
            .unwrap_or(stages.len());
        stages.insert(
            at,
            Stage {
                name: "prep",
                wall,
                peak_rss: nomics_mem::peak_rss_bytes(),
            },
        );
        eprintln!(
            "prep overlapped with alignment ({:.1}s)",
            wall.as_secs_f64()
        );
    }
    aligning.store(false, std::sync::atomic::Ordering::Relaxed);
    report_gpu_fallbacks();
    nomics_align::mem::release_gpu_buffers();
    drop(bam_tx);
    let drain = Instant::now();
    let markdup_input = join(bam_writer, "mapped BAM writer")?;
    // Every chunk's spill file is complete: record the set.
    if spill_mapped && let MarkdupInput::Files(files) = &markdup_input {
        nomics_io::spill::write_manifest(&mapped_dir, files)?;
    }
    let (index_bytes, huge) = idx.memory();
    eprintln!(
        "index after alignment: {:.1} GB, {} on 2 MiB pages",
        index_bytes as f64 / 1e9,
        huge.map_or("unknown".to_owned(), |h| format!(
            "{:.1} GB",
            h as f64 / 1e9
        ))
    );
    eprintln!(
        "align tail: mapped BAM writer drained {:.1}s after alignment",
        drain.elapsed().as_secs_f64()
    );
    drop(idx);
    stage(&mut stages, "align", started);
    finish_run(
        &args,
        threads,
        stages,
        total_started,
        markdup_input,
        Resume::Markdup,
        (&fasta_path, &dbsnp_path, known_sites),
    )
}

/// Whole-genome runs build BaseRecalibrator's tables from md.bam's records
/// as markdup writes them, so BQSR reads md.bam once (to apply) instead of
/// twice; `NOMICS_BQSR_IN_MARKDUP=0` restores the two reads.
fn bqsr_in_markdup() -> bool {
    std::env::var("NOMICS_BQSR_IN_MARKDUP").map_or(true, |v| v != "0")
}

/// [`nomics_markdup::mark_files`], handing md.bam's records, as markdup
/// writes them, to consumers on their own threads: a
/// [`nomics_bqsr::TableBuilder`] (`bqsr_opts`, I-06) and an md.cram encoder
/// (`cram`). Each sees exactly the header and records (duplicate flags set)
/// written to the md stream. The builder loads the reference and known
/// sites while markdup's first pass runs; its table is written, and md.cram
/// completed, only once markdup has written its last record.
fn mark_files_teed<'s>(
    scope: &'s std::thread::Scope<'s, '_>,
    md_opts: &nomics_markdup::Options,
    bqsr_opts: Option<&'s nomics_bqsr::Options>,
    cram: Option<(&'static str, Consumer<'s>)>,
) -> io::Result<(
    nomics_markdup::DuplicationMetrics,
    Option<nomics_bqsr::BuiltTables>,
)> {
    let header = nomics_markdup::files_header(md_opts)?;
    let built = Arc::new(std::sync::Mutex::new(None));
    let mut consumers: Vec<(&'static str, Consumer<'s>)> = Vec::new();
    if let Some(bqsr_opts) = bqsr_opts {
        let header = header.clone();
        let out = built.clone();
        consumers.push((
            "BQSR table builder",
            Box::new(
                move |rx: Receiver<StreamChunk>, spare: Sender<Arc<Vec<u8>>>| {
                    let started = Instant::now();
                    let mut b = nomics_bqsr::TableBuilder::new(&header, bqsr_opts)?;
                    eprintln!(
                        "bqsr in markdup: reference and known sites loaded in {:.1}s",
                        started.elapsed().as_secs_f64()
                    );
                    let mut busy = Duration::ZERO;
                    // The md.cram encoder's copy carries the header first.
                    let mut first = true;
                    while let Ok(batch) = rx.recv() {
                        let Some(batch) = batch else {
                            eprintln!(
                                "bqsr in markdup: recalibrator over {} records, busy {:.1}s",
                                b.records(),
                                busy.as_secs_f64()
                            );
                            *out.lock().map_err(|_| other("table lock poisoned"))? =
                                Some(b.finish()?);
                            return Ok(());
                        };
                        if std::mem::take(&mut first) {
                            continue;
                        }
                        let started = Instant::now();
                        b.add_stream(&batch)?;
                        busy += started.elapsed();
                        let _ = spare.send(batch);
                    }
                    Err(other("markdup stopped before its last record"))
                },
            ),
        ));
    }
    consumers.extend(cram);
    let metrics = tee_stage(scope, consumers, |tee| {
        // The stream's header bytes, as markdup writes them, go first.
        let mut head = Vec::new();
        header.write(&mut head)?;
        tee.send(&mut head)?;
        nomics_markdup::mark_files_with(md_opts, |stream| tee.send(stream))
    })?;
    let built = built
        .lock()
        .map_err(|_| other("table lock poisoned"))?
        .take();
    if bqsr_opts.is_some() && built.is_none() {
        return Err(other("BQSR table builder finished without tables"));
    }
    Ok((metrics, built))
}

/// Where a whole-genome run streams `<prefix>.<kind>` (md or recal) between
/// its stages: as a BAM, or (`internal`: `--save-output-as-bam false`
/// without `--keep-bams`, where nothing publishes it) as a light spill file,
/// which is cheaper to write and to read back.
fn stream_path(outdir: &Path, prefix: &str, kind: &str, internal: bool) -> PathBuf {
    let ext = if internal {
        nomics_io::spill::EXTENSION
    } else {
        "bam"
    };
    outdir.join(format!("{prefix}.{kind}.{ext}"))
}

/// Which write of a stream a file is: recorded beside it in `<stream>.gen`
/// once it is complete, and in the `.done` marker of every CRAM made from
/// it, so a CRAM counts as complete only for the stream it was made from.
#[derive(Debug, Clone, PartialEq, Eq)]
struct StreamId {
    /// `None`: complete, but with no generation recorded for this file
    /// (written by an older binary, or the run stopped before recording
    /// it). Such a stream is never trusted as a CRAM's source, nor as
    /// HaplotypeCaller's input on resume.
    generation: Option<String>,
    /// Orders generations (when written).
    order: u128,
    /// For a recal stream: the generation of the md stream it was made from.
    parent: Option<String>,
}

/// A fresh generation: write time (ns, hex) and process.
fn new_generation() -> String {
    static N: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    format!("{ns:x}-{}-{n}", std::process::id())
}

fn generation_order(generation: &str) -> u128 {
    generation
        .split('-')
        .next()
        .and_then(|h| u128::from_str_radix(h, 16).ok())
        .unwrap_or(0)
}

/// `<stream>.gen`. The stream's writer removes it (durably) just before a
/// new stream takes the name ([`nomics_io::bamstream::RECORD_SUFFIX`]), so
/// it never describes another file; the run writes the new one once the
/// stream is complete.
fn gen_path(stream: &Path) -> PathBuf {
    with_suffix(stream, nomics_io::bamstream::RECORD_SUFFIX)
}

/// The identity of a complete stream at `path`, or `None` if it is missing
/// or incomplete (no BGZF EOF marker / spill footer). A complete stream
/// whose `.gen` is missing or does not match it is unrecorded
/// ([`StreamId::generation`] `None`), ordered by its mtime.
fn stream_id(path: &Path) -> Option<StreamId> {
    let meta = std::fs::metadata(path).ok()?;
    nomics_io::bamstream::check_complete(path).ok()?;
    let recorded = || -> Option<(String, Option<String>)> {
        let text = std::fs::read_to_string(gen_path(path)).ok()?;
        let mut lines = text.lines();
        (lines.next()? == "nomics-stream 1").then_some(())?;
        let generation = lines.next()?.strip_prefix("generation ")?.to_owned();
        let size: u64 = lines.next()?.strip_prefix("size ")?.parse().ok()?;
        let parent = lines.next()?.strip_prefix("parent ")?;
        (lines.next().is_none() && size == meta.len()).then_some(())?;
        Some((generation, (parent != "-").then(|| parent.to_owned())))
    };
    Some(match recorded() {
        Some((generation, parent)) => StreamId {
            order: generation_order(&generation),
            generation: Some(generation),
            parent,
        },
        None => StreamId {
            generation: None,
            order: meta
                .modified()
                .ok()?
                .duration_since(std::time::UNIX_EPOCH)
                .ok()?
                .as_nanos(),
            parent: None,
        },
    })
}

/// A resume at HaplotypeCaller takes only a recal stream with a recorded
/// generation and md parent, beside (if any) the md stream it was made
/// from: a run that rewrote md and stopped before recal leaves the old
/// recal beside the new md, which must not be published together.
fn check_hc_resume(
    md: &Path,
    md_id: Option<&StreamId>,
    recal: &Path,
    recal_id: Option<&StreamId>,
) -> io::Result<String> {
    let id = recal_id.ok_or_else(|| {
        other(format!("no complete recal stream to resume HaplotypeCaller from ({} or its spill file)", recal.display()))
    })?;
    let (Some(_), Some(parent)) = (&id.generation, &id.parent) else {
        return Err(other(format!(
            "{} has no recorded generation (an older nomics, or a run stopped before recording it); resume from bqsr",
            recal.display()
        )));
    };
    if let Some(m) = md_id
        && m.generation.as_ref() != Some(parent)
    {
        return Err(other(format!(
            "{} was not made from {} (a rewrite stopped part-way?); resume from bqsr",
            recal.display(),
            md.display()
        )));
    }
    Ok(parent.clone())
}

/// Record `generation` for the complete stream at `path` (durably: written,
/// synced, renamed, directory synced).
fn write_stream_id(path: &Path, generation: &str, parent: Option<&str>) -> io::Result<()> {
    nomics_io::bamstream::check_complete(path)?;
    let size = std::fs::metadata(path)?.len();
    let text = format!(
        "nomics-stream 1\ngeneration {generation}\nsize {size}\nparent {}\n",
        parent.unwrap_or("-")
    );
    let sidecar = gen_path(path);
    let tmp = with_suffix(&sidecar, ".tmp");
    let mut f = File::create(&tmp)?;
    f.write_all(text.as_bytes())?;
    f.sync_all()?;
    drop(f);
    std::fs::rename(&tmp, &sidecar)?;
    sync_parent(path)
}

fn sync_parent(path: &Path) -> io::Result<()> {
    let dir = path.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    File::open(dir)?.sync_all()
}

/// Before a run rewrites the `<prefix>.<kind>` stream: the `.done` of the
/// CRAM made from it goes. The old stream keeps its `.gen` until the new one
/// takes its name (its writer drops it then), and the other format's copy
/// stays until the new stream is recorded ([`retire_alternate`]).
fn begin_stream_rewrite(outdir: &Path, prefix: &str, kind: &str) -> io::Result<()> {
    let cram = sarek_cram_path(outdir, prefix, kind);
    remove_if_present(&done_path(&cram))?;
    if cram.parent().is_some_and(|d| d.is_dir()) {
        sync_parent(&cram)?;
    }
    Ok(())
}

/// Once the stream at `written` is complete and its generation recorded,
/// remove the other format's copy of it (and its `.gen`), if any.
fn retire_alternate(written: &Path, other: &Path) -> io::Result<()> {
    if written == other || !other.exists() {
        return Ok(());
    }
    remove_if_present(&gen_path(other))?;
    remove_if_present(other)?;
    eprintln!(
        "removed {} (replaced by {})",
        other.display(),
        written.display()
    );
    sync_parent(other)
}

/// The `<prefix>.<kind>` stream an earlier run left (resumed runs): its BAM
/// or its spill file, whichever is complete, the later generation if both
/// are; the BAM's path and `None` if neither is.
fn find_stream(outdir: &Path, prefix: &str, kind: &str) -> (PathBuf, Option<StreamId>) {
    let bam = stream_path(outdir, prefix, kind, false);
    let spill = stream_path(outdir, prefix, kind, true);
    match (stream_id(&bam), stream_id(&spill)) {
        (Some(b), Some(s)) => {
            let (path, id) = if s.order > b.order { (spill, s) } else { (bam, b) };
            eprintln!(
                "resume: {prefix}.{kind} is complete as BAM and as spill file; using the later, {}",
                path.display()
            );
            (path, Some(id))
        }
        (None, Some(s)) => (spill, Some(s)),
        (b, None) => {
            if b.is_none() && spill.exists() {
                eprintln!("resume: {} is incomplete; not used", spill.display());
            }
            (bam, b)
        }
    }
}

/// Markdup, BQSR and HaplotypeCaller from `from` on.
fn finish_run(
    args: &Args,
    threads: usize,
    mut stages: Vec<Stage>,
    total_started: Instant,
    markdup_input: MarkdupInput,
    from: Resume,
    (fasta_path, dbsnp_path, known_sites): (&Path, &Path, Vec<PathBuf>),
) -> io::Result<()> {
    // Whole genomes stream md and recal between the stages (on disk).
    let streamed = matches!(markdup_input, MarkdupInput::Files(_) | MarkdupInput::Done);
    // With CRAM output and no `--keep-bams`, those streams are internal:
    // light spill files, removed once the CRAM made from each is complete.
    let internal = streamed && !args.save_output_as_bam && !args.keep_bams;
    // The streams this run writes (each with a fresh generation), and the
    // complete ones a resumed run reads ([`find_stream`]).
    let (outdir, prefix) = (&args.outdir, &args.prefix);
    let (md_output, md_id) = if streamed && from > Resume::Markdup {
        find_stream(outdir, prefix, "md")
    } else {
        (stream_path(outdir, prefix, "md", internal), None)
    };
    let (recal_bam, recal_id) = if streamed && from > Resume::Bqsr {
        find_stream(outdir, prefix, "recal")
    } else {
        (stream_path(outdir, prefix, "recal", internal), None)
    };
    // A resume at BQSR may take an unrecorded md stream (BQSR rewrites recal
    // from it): it is recorded now, under a new generation no CRAM marker
    // names, so md.cram is regenerated from it.
    let md_id = match md_id {
        Some(id) if id.generation.is_none() && from == Resume::Bqsr => {
            let generation = new_generation();
            eprintln!(
                "resume: {} has no recorded generation; recording {generation} (md.cram is made again from it)",
                md_output.display()
            );
            write_stream_id(&md_output, &generation, None)?;
            Some(StreamId {
                generation: Some(generation),
                ..id
            })
        }
        id => id,
    };
    // A resume at HaplotypeCaller: recal's md parent.
    let recal_parent = if streamed && from == Resume::Hc {
        Some(check_hc_resume(&md_output, md_id.as_ref(), &recal_bam, recal_id.as_ref())?)
    } else {
        None
    };
    let md_gen = (from <= Resume::Markdup).then(new_generation);
    let recal_gen = (from <= Resume::Bqsr).then(new_generation);
    // The generation of the md stream BQSR reads (the parent of recal).
    let md_source = md_gen
        .clone()
        .or_else(|| md_id.as_ref().and_then(|i| i.generation.clone()));
    // A rewritten stream loses its recorded generation, and the CRAM made
    // from it its completion marker, before the rewrite starts.
    for (kind, rewritten) in [("md", md_gen.is_some()), ("recal", recal_gen.is_some())] {
        if rewritten {
            begin_stream_rewrite(outdir, prefix, kind)?;
        }
    }
    // Once a rewritten stream is durable: record its generation, then
    // retire the other format's copy (until then, older, it is ignored).
    let record_stream = |kind: &str, path: &Path, generation: &str, parent: Option<&str>| {
        write_stream_id(path, generation, parent)?;
        let other = stream_path(outdir, prefix, kind, !nomics_io::spill::is_spill_path(path));
        retire_alternate(path, &other)
    };
    let recal_table = args.outdir.join(format!("{}.recal.table", args.prefix));
    let vcf = args.outdir.join(format!(
        "{}.haplotypecaller.{}vcf.gz",
        args.prefix,
        if args.emit_gvcf { "g." } else { "" }
    ));
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let md_opts = nomics_markdup::Options {
        inputs: Vec::new(),
        output: md_output.clone(),
        metrics: args.outdir.join(format!("{}.md.metrics", args.prefix)),
        pg_version: format!("Version:{}", env!("CARGO_PKG_VERSION")),
        pg_command_line: command_line.clone(),
        threads,
    };
    let groups = crate::scatter::split_bed(&std::fs::read_to_string(&args.targets)?, 200000.0)?;
    let n_groups = groups.len();
    let group_paths = crate::scatter::write_groups(&groups, &args.outdir.join("intervals"))?;
    let order = crate::scatter::orders(&groups, args.scatter_order.as_deref())?;
    // samtools merge seeds its ID suffixes with the wall clock.
    let merge_seed = order.merge_seed.unwrap_or_else(|| {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs() as i64)
    });
    eprintln!(
        "scatter: {} interval groups, merge seed {merge_seed}",
        groups.len()
    );
    let bqsr_opts = nomics_bqsr::Options {
        input: md_output.clone(),
        reference: fasta_path.to_path_buf(),
        known_sites,
        intervals: args.targets.clone(),
        scatter: (groups.len() > 1).then(|| nomics_bqsr::Scatter {
            groups: group_paths.clone(),
            gather_order: order.gather,
            merge_order: order.merge,
            merge_seed,
        }),
        known_sites_cache: args.known_sites_cache.clone(),
        output_table: recal_table,
        output_bam: Some(recal_bam.clone()),
        pg_command_line: command_line,
        threads,
    };
    if streamed && !args.save_output_as_bam {
        // What this run writes (and, internal, removes) or was handed.
        let role = |writes: bool| match (writes, internal) {
            (false, _) => "read from an earlier run, kept",
            (true, true) => "internal, removed once its CRAM is complete",
            (true, false) => "kept (--keep-bams)",
        };
        eprintln!(
            "CRAM output: md.cram encoded from markdup's last pass, recal.cram from ApplyBQSR's output; md stream {} ({}), recal stream {} ({})",
            md_output.display(),
            role(from <= Resume::Markdup),
            recal_bam.display(),
            role(from <= Resume::Bqsr),
        );
    }
    // `--save-output-as-bam false`: the FASTA for the CRAM encoder, loaded
    // while markdup starts.
    let cram_fasta =
        (!args.save_output_as_bam).then(|| LazyFasta::load(fasta_path, threads.min(8)));
    let cram_fasta = cram_fasta.as_ref();
    let started = Instant::now();
    // Whole genomes: BaseRecalibrator's tables, built from md.bam's records
    // as markdup writes them (`NOMICS_BQSR_IN_MARKDUP`).
    let mut built_tables = None;
    // `None`: the files path wrote the md stream and later stages read it.
    let (md_set, metrics) = match markdup_input {
        MarkdupInput::Files(inputs) => {
            let opts = nomics_markdup::Options { inputs, ..md_opts };
            let metrics = std::thread::scope(|scope| -> io::Result<_> {
                let source = md_gen.as_deref().unwrap_or("-");
                let cram = cram_fasta
                    .map(|fasta| cram_consumer(fasta, fasta_path, args, "md", threads, source, 1));
                let tables = bqsr_in_markdup().then_some(&bqsr_opts);
                let (metrics, built) = mark_files_teed(scope, &opts, tables, cram)?;
                built_tables = built;
                Ok(metrics)
            })?;
            record_stream("md", &md_output, md_gen.as_deref().unwrap_or("-"), None)?;
            // A run's own spill (`--save-mapped false`) has been read for
            // the last time; spill files a resumed run was handed are kept.
            // The manifest goes first, so an interrupted cleanup leaves a
            // set no resume accepts.
            let spill: Vec<PathBuf> = opts
                .inputs
                .iter()
                .filter(|p| nomics_io::spill::is_spill_path(p))
                .cloned()
                .collect();
            if args.resume_from.is_none() && !spill.is_empty() {
                nomics_io::spill::remove_spill(&args.outdir.join("mapped"), &spill)?;
            }
            (None, metrics)
        }
        MarkdupInput::Raw(sets) => {
            let (set, metrics) = nomics_markdup::mark_sets(sets, &md_opts)?;
            (Some(set), metrics)
        }
        MarkdupInput::Prepared(plan, sets) => {
            let (set, metrics) = plan.finish(sets, &md_opts)?;
            (Some(set), metrics)
        }
        MarkdupInput::Done => (None, Default::default()),
    };
    if from <= Resume::Markdup {
        eprintln!("markdup: {metrics:?}");
        stage(&mut stages, "markdup", started);
    }

    let started = Instant::now();
    // `None`: streaming BQSR wrote the recal stream, and HaplotypeCaller reads it.
    let recal_set = match &md_set {
        None if from == Resume::Hc => {
            // A resumed run regenerates an md.cram an earlier run left
            // missing, incomplete or made from another md stream than recal's
            // parent ([`check_hc_resume`]: the md stream beside it, if any).
            if let Some(fasta) = cram_fasta {
                let parent = recal_parent.as_deref().unwrap_or("-");
                if !cram_complete(outdir, prefix, "md", Some(parent)) {
                    if md_id.is_none() {
                        return Err(other(format!(
                            "md.cram is not complete for this run's recal stream, and there is no complete md stream ({}) to make it from",
                            md_output.display()
                        )));
                    }
                    let fasta = fasta.get()?;
                    write_stream_cram(&md_output, parent, fasta, fasta_path, args, "md", threads, 1)?;
                }
            }
            None
        }
        None => std::thread::scope(|scope| -> io::Result<_> {
            // A run resumed at BQSR regenerates an incomplete md.cram from
            // the md stream, alongside BQSR (markdup's own pass wrote it
            // otherwise).
            if from == Resume::Bqsr && md_id.is_none() {
                return Err(other(format!(
                    "no complete md stream to resume BQSR from ({} or its spill file)",
                    md_output.display()
                )));
            }
            let md_source = md_source.as_deref().unwrap_or("-");
            let md_cram = cram_fasta
                .filter(|_| {
                    from == Resume::Bqsr && !cram_complete(outdir, prefix, "md", Some(md_source))
                })
                .map(|fasta| {
                    let md_output = &md_output;
                    scope.spawn(move || {
                        let fasta = fasta.get()?;
                        write_stream_cram(md_output, md_source, fasta, fasta_path, args, "md", threads, 1)
                    })
                });
            // recal.cram from ApplyBQSR's output as it is written.
            let recal_source = recal_gen.as_deref().unwrap_or("-");
            let consumers: Vec<_> = cram_fasta
                .map(|fasta| cram_consumer(fasta, fasta_path, args, "recal", threads, recal_source, n_groups))
                .into_iter()
                .collect();
            // BQSR's input is our own md stream: it must be complete.
            nomics_io::bamstream::check_complete(&md_output)?;
            let applied = tee_stage(scope, consumers, |tee| match built_tables {
                Some(built) => {
                    nomics_bqsr::apply_streaming_with(&bqsr_opts, built, |s| tee.send(s))
                }
                None => nomics_bqsr::run_streaming_with(&bqsr_opts, |s| tee.send(s)),
            });
            let md_cram = md_cram.map(|h| join_scoped(h, "md.cram writer"));
            applied?;
            md_cram.transpose()?;
            record_stream("recal", &recal_bam, recal_source, Some(md_source))?;
            Ok(None)
        })?,
        Some(md_set) => Some(std::thread::scope(|scope| -> io::Result<RecordSet> {
            let writer = scope.spawn(|| match cram_fasta {
                Some(fasta) => {
                    write_set_cram(md_set, fasta.get()?, fasta_path, args, "md", threads, 1)
                }
                None => nomics_io::bamstream::write_set(md_set, &md_output, threads),
            });
            let recal = nomics_bqsr::run_on_set(md_set, &bqsr_opts);
            let written = writer
                .join()
                .map_err(|_| other("markdup BAM writer thread panicked"))?;
            written?;
            if cram_fasta.is_none() {
                record_stream("md", &md_output, md_gen.as_deref().unwrap_or("-"), None)?;
            }
            recal
        })?),
    };
    drop(md_set);
    if from <= Resume::Bqsr {
        stage(&mut stages, "bqsr", started);
    }
    // The md stream has been read for the last time; an internal one this
    // run wrote goes now (md.cram is complete), halving the streams on disk
    // while HaplotypeCaller runs.
    if let Some(generation) = md_gen.as_deref().filter(|_| internal) {
        remove_internal_stream(outdir, prefix, "md", &md_output, generation)?;
    }
    // For timing experiments on resumed whole-genome runs.
    if std::env::var("NOMICS_STOP_AFTER").is_ok_and(|s| s == "bqsr") {
        if let Some(recal_set) = &recal_set {
            nomics_io::bamstream::write_set(recal_set, &recal_bam, threads)?;
            let generation = recal_gen.as_deref().unwrap_or("-");
            record_stream("recal", &recal_bam, generation, md_gen.as_deref())?;
        }
        return write_timings(&args.outdir, &stages, total_started.elapsed());
    }

    let started = Instant::now();
    let gvcf_text = match &recal_set {
        None => std::thread::scope(|scope| -> io::Result<Option<String>> {
            // A run resumed at HaplotypeCaller regenerates an incomplete
            // recal.cram from the recal stream, alongside it (ApplyBQSR's
            // own pass wrote it otherwise).
            // Recorded ([`check_hc_resume`]) on a resume at HaplotypeCaller.
            let source = recal_id.as_ref().and_then(|i| i.generation.as_deref()).unwrap_or("-");
            let need = from == Resume::Hc && !cram_complete(outdir, prefix, "recal", Some(source));
            let cram = cram_fasta.filter(|_| need).map(|fasta| {
                let recal_bam = &recal_bam;
                scope.spawn(move || {
                    let fasta = fasta.get()?;
                    write_stream_cram(recal_bam, source, fasta, fasta_path, args, "recal", threads, n_groups)
                })
            });
            // Our own recal stream: it must be complete.
            nomics_io::bamstream::check_complete(&recal_bam)?;
            let called = call_variants(
                HcInput::File(&recal_bam),
                &group_paths,
                fasta_path,
                dbsnp_path,
                &vcf,
                threads,
                args.emit_gvcf,
            );
            if let Some(cram) = cram {
                join_scoped(cram, "recal.cram writer")?;
            }
            called
        })?,
        Some(recal_set) => std::thread::scope(|scope| -> io::Result<Option<String>> {
            let writer = scope.spawn(|| match cram_fasta {
                Some(fasta) => {
                    write_set_cram(recal_set, fasta.get()?, fasta_path, args, "recal", threads, n_groups)
                }
                None => nomics_io::bamstream::write_set(recal_set, &recal_bam, threads),
            });
            let called = call_variants(
                HcInput::Set(recal_set),
                &group_paths,
                fasta_path,
                dbsnp_path,
                &vcf,
                threads,
                args.emit_gvcf,
            );
            let written = writer
                .join()
                .map_err(|_| other("recalibrated BAM writer thread panicked"))?;
            written?;
            if cram_fasta.is_none() {
                let generation = recal_gen.as_deref().unwrap_or("-");
                record_stream("recal", &recal_bam, generation, md_gen.as_deref())?;
            }
            called
        })?,
    };
    drop(recal_set);
    stage(&mut stages, "haplotypecaller", started);

    // sarek --joint_germline: GenotypeGVCFs on the one-sample gVCF.
    if let Some(text) = gvcf_text {
        let started = Instant::now();
        let dir = args.outdir.join("joint_variant_calling");
        std::fs::create_dir_all(&dir)?;
        let file = File::create(dir.join("joint_germline.vcf.gz"))?;
        let mut w = ParallelBgzfWriter::new(file, threads);
        let n = nomics_joint::genotype_gvcf_text(&text, Some(dbsnp_path), threads, &mut w)?;
        w.finish()?.flush()?;
        eprintln!("joint genotyping: {n} records");
        stage(&mut stages, "joint", started);
    }
    // The recal stream has been read for the last time. Streams a resumed
    // run was handed are kept (as with spill files).
    if let Some(generation) = recal_gen.as_deref().filter(|_| internal) {
        remove_internal_stream(outdir, prefix, "recal", &recal_bam, generation)?;
    }
    write_timings(&args.outdir, &stages, total_started.elapsed())
}

/// Remove an internal stream (a spill file nothing publishes, of
/// `generation`) and its `.gen` once the `kind` CRAM made from it and its
/// `.crai` are complete and durable; if the CRAM is not complete, or was
/// made from another stream, keep it (to regenerate the CRAM from).
fn remove_internal_stream(
    outdir: &Path,
    prefix: &str,
    kind: &str,
    stream: &Path,
    generation: &str,
) -> io::Result<()> {
    if !cram_complete(outdir, prefix, kind, Some(generation)) {
        eprintln!("keeping {}: {kind}.cram is not complete", stream.display());
        return Ok(());
    }
    // The CRAM, .crai and their directories are durable before the stream goes.
    let cram = sarek_cram_path(outdir, prefix, kind);
    File::open(&cram)?.sync_all()?;
    File::open(crai_path(&cram))?.sync_all()?;
    sync_dirs(cram.parent().unwrap(), outdir)?;
    remove_if_present(&gen_path(stream))?;
    match std::fs::remove_file(stream) {
        Ok(()) => eprintln!(
            "removed {} (internal stream, {kind}.cram complete; --keep-bams keeps BAMs instead)",
            stream.display()
        ),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {}
        Err(e) => return Err(e),
    }
    Ok(())
}

#[cfg(test)]
mod cram_output_tests {
    use super::*;

    fn tiny_fasta(dir: &Path) -> (PathBuf, Fasta, nomics_io::bam::Header) {
        let fa = dir.join("ref.fa");
        std::fs::write(&fa, ">c1\nACGTACGTAC\n").unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "c1\t10\t4\t10\t11\n").unwrap();
        let fasta = Fasta::load(&fa, 1).unwrap();
        let header = nomics_io::bam::Header {
            text: b"@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\n".to_vec(),
            refs: vec![nomics_io::bam::Reference {
                name: b"c1".to_vec(),
                len: 10,
            }],
        };
        (fa, fasta, header)
    }

    fn scratch(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nomics-cram-out-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn failed_cram_write_leaves_nothing_complete() {
        let dir = scratch("fail");
        let (fa, fasta, header) = tiny_fasta(&dir);
        // A complete pair from an earlier run, then a failing rewrite.
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, "-", 1, |_| Ok(())).unwrap();
        assert!(cram_complete(&dir, "s", "md", None));
        let err = write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, "-", 1, |_| {
            Err(io::Error::other("feed failed"))
        });
        assert!(err.is_err());
        assert!(!cram_complete(&dir, "s", "md", None), "a failed rewrite must not look complete");
        // A successful write leaves no temporary files behind.
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, "-", 1, |_| Ok(())).unwrap();
        assert!(cram_complete(&dir, "s", "md", None));
        let cram = sarek_cram_path(&dir, "s", "md");
        assert!(!with_suffix(&cram, ".tmp").exists());
        assert!(!with_suffix(&done_path(&cram), ".tmp").exists());
        assert!(!with_suffix(&crai_path(&cram), ".tmp").exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn completion_marker_must_match_both_files() {
        let dir = scratch("marker");
        let (fa, fasta, header) = tiny_fasta(&dir);
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, "-", 1, |_| Ok(())).unwrap();
        let cram = sarek_cram_path(&dir, "s", "md");
        let crai = crai_path(&cram);
        assert!(cram_complete(&dir, "s", "md", None));
        let text = std::fs::read_to_string(done_path(&cram)).unwrap();
        let m = DoneMarker::parse(&text).unwrap();
        assert_eq!(m.text(), text);
        assert_eq!(m.cram, file_crc(&cram).unwrap());
        // Truncated CRAM.
        let full = std::fs::read(&cram).unwrap();
        std::fs::write(&cram, &full[..full.len() - 1]).unwrap();
        assert!(!cram_complete(&dir, "s", "md", None), "truncated CRAM accepted");
        std::fs::write(&cram, &full).unwrap();
        assert!(cram_complete(&dir, "s", "md", None));
        // A stale .crai of the same size.
        let idx = std::fs::read(&crai).unwrap();
        let mut stale = idx.clone();
        *stale.last_mut().unwrap() ^= 1;
        std::fs::write(&crai, &stale).unwrap();
        assert!(!cram_complete(&dir, "s", "md", None), "stale .crai accepted");
        std::fs::write(&crai, &idx).unwrap();
        // No marker, or a garbled one.
        std::fs::write(done_path(&cram), "nomics-cram-done 1\n").unwrap();
        assert!(!cram_complete(&dir, "s", "md", None));
        std::fs::remove_file(done_path(&cram)).unwrap();
        assert!(!cram_complete(&dir, "s", "md", None));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn create_dir_all_synced_creates_nested_dirs() {
        let dir = scratch("mkdir");
        let deep = dir.join("a/b/c");
        create_dir_all_synced(&deep).unwrap();
        assert!(deep.is_dir());
        create_dir_all_synced(&deep).unwrap(); // already there
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn internal_stream_removed_only_after_its_cram() {
        let dir = scratch("rm");
        let (fa, fasta, header) = tiny_fasta(&dir);
        let md = stream_path(&dir, "s", "md", true);
        let recal = stream_path(&dir, "s", "recal", true);
        assert_eq!(md, dir.join("s.md.nspill"));
        for f in [&md, &recal] {
            write_stream(f, &header);
            write_stream_id(f, "g1", None).unwrap();
        }
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, "g1", 1, |_| Ok(())).unwrap();
        remove_internal_stream(&dir, "s", "recal", &recal, "g1").unwrap();
        assert!(recal.exists(), "recal.cram missing: its stream must stay");
        remove_internal_stream(&dir, "s", "md", &md, "g1").unwrap();
        assert!(!md.exists() && !gen_path(&md).exists());
        // An incomplete CRAM (no marker) keeps its stream.
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "recal", 1, "g1", 1, |_| Ok(())).unwrap();
        std::fs::remove_file(done_path(&sarek_cram_path(&dir, "s", "recal"))).unwrap();
        remove_internal_stream(&dir, "s", "recal", &recal, "g1").unwrap();
        assert!(recal.exists());
        // So does a CRAM made from another stream.
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "recal", 1, "g0", 1, |_| Ok(())).unwrap();
        remove_internal_stream(&dir, "s", "recal", &recal, "g1").unwrap();
        assert!(recal.exists());
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "recal", 1, "g1", 1, |_| Ok(())).unwrap();
        remove_internal_stream(&dir, "s", "recal", &recal, "g1").unwrap();
        assert!(!recal.exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// A complete stream (header only) at `path`, BAM or spill by extension.
    fn write_stream(path: &Path, header: &nomics_io::bam::Header) {
        let mut w = nomics_io::bamstream::StreamWriter::create(path, 1).unwrap();
        header.write(&mut w).unwrap();
        w.finish().unwrap();
    }

    /// Resume picks a complete stream, by recorded generation; a BAM cut
    /// at a block boundary (no EOF marker), or one whose `.gen` belongs to
    /// another file, is not complete.
    #[test]
    fn resume_takes_the_latest_complete_stream() {
        let dir = scratch("find");
        let (_, _, header) = tiny_fasta(&dir);
        let bam = stream_path(&dir, "s", "md", false);
        let spill = stream_path(&dir, "s", "md", true);
        let (path, id) = find_stream(&dir, "s", "md");
        assert_eq!((path, id), (bam.clone(), None), "neither: the BAM's path");
        write_stream(&spill, &header);
        write_stream_id(&spill, &new_generation(), None).unwrap();
        write_stream(&bam, &header);
        let older = new_generation();
        let newer = new_generation();
        write_stream_id(&bam, &older, None).unwrap();
        write_stream_id(&spill, &newer, None).unwrap();
        assert_eq!(find_stream(&dir, "s", "md").0, spill, "the later generation");
        write_stream_id(&bam, &new_generation(), None).unwrap();
        assert_eq!(find_stream(&dir, "s", "md").0, bam);
        // The BAM cut before its EOF marker: a complete-looking stream
        // that would read cleanly, but is not taken.
        let full = std::fs::read(&bam).unwrap();
        std::fs::write(&bam, &full[..full.len() - 28]).unwrap();
        assert!(stream_id(&bam).is_none());
        assert_eq!(find_stream(&dir, "s", "md").0, spill);
        // A `.gen` for another file (size differs) is not this file's.
        std::fs::write(&bam, &full).unwrap();
        let mut longer = full.clone();
        longer.extend_from_slice(&nomics_io::bgzf::EOF_BLOCK);
        std::fs::write(&bam, &longer).unwrap();
        let id = stream_id(&bam).unwrap();
        assert_eq!(id.generation, None, "a .gen for another file is not this file's");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// The reviewer's case: run A publishes md.cram from its md stream, run
    /// B rewrites the stream; A's CRAM must not count as complete for B's.
    #[test]
    fn cram_marker_is_bound_to_its_source_stream() {
        let dir = scratch("bound");
        let (fa, fasta, header) = tiny_fasta(&dir);
        let md = stream_path(&dir, "s", "md", false);
        write_stream(&md, &header);
        let a = new_generation();
        write_stream_id(&md, &a, None).unwrap();
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, &a, 1, |_| Ok(())).unwrap();
        assert!(cram_complete(&dir, "s", "md", Some(&a)));
        // Run B rewrites md (BAM mode: no CRAM written).
        begin_stream_rewrite(&dir, "s", "md").unwrap();
        assert!(!cram_complete(&dir, "s", "md", None), "rewrite drops the marker");
        write_stream(&md, &header);
        let b = new_generation();
        write_stream_id(&md, &b, None).unwrap();
        // Even with A's marker back, it names A's stream, not B's.
        write_sarek_cram(&header, &fasta, &fa, &dir, "s", "md", 1, &a, 1, |_| Ok(())).unwrap();
        let (_, id) = find_stream(&dir, "s", "md");
        assert_eq!(id.as_ref().unwrap().generation.as_ref(), Some(&b));
        assert!(!cram_complete(&dir, "s", "md", Some(&b)));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// The reviewer's case: run A leaves md and recal; run B rewrites md and
    /// stops before recal. Resuming at HaplotypeCaller must not publish A's
    /// recal beside B's md, whether B recorded its md or stopped before.
    #[test]
    fn hc_resume_refuses_recal_not_made_from_md() {
        let dir = scratch("parent");
        let (_, _, header) = tiny_fasta(&dir);
        let md = stream_path(&dir, "s", "md", false);
        let recal = stream_path(&dir, "s", "recal", false);
        let a = new_generation();
        write_stream(&md, &header);
        write_stream_id(&md, &a, None).unwrap();
        write_stream(&recal, &header);
        write_stream_id(&recal, &new_generation(), Some(&a)).unwrap();
        let ids = || (find_stream(&dir, "s", "md").1, find_stream(&dir, "s", "recal").1);
        let (m, r) = ids();
        assert_eq!(check_hc_resume(&md, m.as_ref(), &recal, r.as_ref()).unwrap(), a);
        // B's md takes the name: its writer drops A's record first, so
        // until B records it, md is unrecorded (not A's).
        write_stream(&md, &header);
        assert!(!gen_path(&md).exists());
        let (m, r) = ids();
        assert_eq!(m.as_ref().unwrap().generation, None);
        assert!(check_hc_resume(&md, m.as_ref(), &recal, r.as_ref()).is_err());
        // B recorded its md: still not recal's parent.
        write_stream_id(&md, &new_generation(), None).unwrap();
        let (m, r) = ids();
        assert!(check_hc_resume(&md, m.as_ref(), &recal, r.as_ref()).is_err());
        // An unrecorded recal is never an HC input, md or not.
        std::fs::remove_file(gen_path(&recal)).unwrap();
        std::fs::remove_file(&md).unwrap();
        let (m, r) = ids();
        assert!(m.is_none());
        assert!(check_hc_resume(&md, m.as_ref(), &recal, r.as_ref()).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// A rewrite keeps the other format's copy until the new stream is
    /// complete and recorded; until then resume still finds the old one.
    #[test]
    fn alternate_copy_retired_only_after_the_new_stream() {
        let dir = scratch("alt");
        let (_, _, header) = tiny_fasta(&dir);
        let bam = stream_path(&dir, "s", "md", false);
        let spill = stream_path(&dir, "s", "md", true);
        write_stream(&bam, &header);
        let old = new_generation();
        write_stream_id(&bam, &old, None).unwrap();
        // The rewrite (as spill) starts, and is cut short: nothing under
        // the spill name, and the BAM is still there and found.
        begin_stream_rewrite(&dir, "s", "md").unwrap();
        {
            let mut w = nomics_io::bamstream::StreamWriter::create(&spill, 1).unwrap();
            header.write(&mut w).unwrap();
            // dropped unfinished
        }
        assert!(!spill.exists());
        assert!(bam.exists());
        assert_eq!(find_stream(&dir, "s", "md").0, bam);
        // Finished and recorded: the BAM goes.
        write_stream(&spill, &header);
        write_stream_id(&spill, &new_generation(), None).unwrap();
        retire_alternate(&spill, &bam).unwrap();
        assert!(!bam.exists() && !gen_path(&bam).exists());
        assert_eq!(find_stream(&dir, "s", "md").0, spill);
        std::fs::remove_dir_all(&dir).unwrap();
    }
    /// A tee'd stream reaches every consumer whole, and the CRAM encoded
    /// from it equals one encoded from the same records directly.
    #[test]
    fn teed_stream_reaches_every_consumer() {
        let dir = scratch("tee");
        let (fa, fasta, header) = tiny_fasta(&dir);
        let mut records = Vec::new();
        // Unmapped, unplaced, no sequence: name `NNNN\0`.
        for i in 0..3000u32 {
            let mut raw = vec![0u8; 32];
            for at in [0, 4, 20, 24] {
                raw[at..at + 4].copy_from_slice(&(-1i32).to_le_bytes());
            }
            raw[8] = 5;
            raw[14..16].copy_from_slice(&4u16.to_le_bytes()); // flag: unmapped
            raw.extend_from_slice(format!("{i:04}\0").as_bytes());
            records.push(raw);
        }
        let rec_len = 4 + records[0].len();
        let mut head = Vec::new();
        header.write(&mut head).unwrap();
        let total: Vec<u8> = records
            .iter()
            .flat_map(|r| {
                (r.len() as i32)
                    .to_le_bytes()
                    .into_iter()
                    .chain(r.iter().copied())
            })
            .collect();
        let seen = std::sync::Mutex::new(Vec::new());
        std::thread::scope(|scope| {
            let echo: Consumer = Box::new(|rx, spare| {
                let mut r = ChunkReader::new(rx, spare);
                let mut all = Vec::new();
                io::Read::read_to_end(&mut r, &mut all)?;
                *seen.lock().unwrap() = all;
                Ok(())
            });
            tee_stage(scope, vec![("echo", echo)], |tee| {
                tee.send(&mut head.clone())?;
                for chunk in total.chunks(1000) {
                    tee.send(&mut chunk.to_vec())?;
                }
                Ok(())
            })
            .unwrap();
        });
        let mut want = head.clone();
        want.extend_from_slice(&total);
        assert_eq!(*seen.lock().unwrap(), want);

        // A producer that fails mid-stream leaves no complete CRAM.
        let r = std::thread::scope(|scope| {
            let consume: Consumer = Box::new(|rx, spare| {
                let mut r = ChunkReader::new(rx, spare);
                let h = nomics_io::bam::Header::read(&mut r)?;
                let mut rec = Vec::new();
                write_sarek_cram(&h, &fasta, &fa, &dir, "s", "md", 1, "-", 1, |sink| {
                    while nomics_io::bam::read_raw_record(&mut r, &mut rec)? {
                        sink.push(&rec)?;
                    }
                    Ok(())
                })
            });
            tee_stage(scope, vec![("md.cram writer", consume)], |tee| {
                tee.send(&mut head.clone())?;
                tee.send(&mut total[..total.len() / 2 / rec_len * rec_len].to_vec())?;
                Err::<(), _>(other("markdup failed"))
            })
        });
        assert_eq!(r.unwrap_err().to_string(), "markdup failed");
        assert!(!cram_complete(&dir, "s", "md", None));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
