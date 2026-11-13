//! Base quality score recalibration, equivalent to GATK 4.6.2.0
//! `BaseRecalibrator` + `ApplyBQSR` as run by sarek (one interval group, three
//! known-sites VCFs, default arguments). Section numbers (§) refer to
//! `docs/spec/bqsr.md`.
//!
//! Gates: `preprocessing/recal_table/*.recal.table` (bytes) and
//! `preprocessing/recalibrated/*.recal.bam` (records).

mod merge;
pub use merge::{CramMerge, add_samtools_pg};
mod reads;
mod stats;
mod tables;

use std::collections::HashMap;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use nomics_io::bam::{Header, RecordSet};
use nomics_io::bamstream::{BamStream, Batch};
use nomics_io::fasta::Fasta;
use nomics_io::vcf::for_each_record;

use reads::{
    Read, clip_for_recalibration, context_keys, cycle_key, is_regular_base, known_site_mask,
    snp_errors,
};
use tables::{Counts, N_QUAL, Recalibration, Tables, build_report, gathered_estimates};

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

const PRESERVE_QSCORES_LESS_THAN: u8 = 6;
const LOW_QUAL_TAIL: u8 = 2;

/// Sarek's interval scatter (§9): BaseRecalibrator and ApplyBQSR run once per
/// group, GatherBQSRReports combines the reports, and `samtools merge` joins
/// the BAMs. Both orders are Nextflow's task-completion order, recorded from
/// the reference run's command lines; they are indices into `groups`.
pub struct Scatter {
    pub groups: Vec<PathBuf>,
    pub gather_order: Vec<usize>,
    pub merge_order: Vec<usize>,
    /// samtools merge's `time(NULL)` seed for its @RG/@PG ID suffixes.
    pub merge_seed: i64,
}

pub struct Options {
    pub input: PathBuf,
    pub reference: PathBuf,
    pub known_sites: Vec<PathBuf>,
    /// The one interval group when `scatter` is `None`.
    pub intervals: PathBuf,
    pub scatter: Option<Scatter>,
    /// Binary cache of parsed known sites (created or refreshed as needed).
    pub known_sites_cache: Option<PathBuf>,
    pub output_table: PathBuf,
    /// ApplyBQSR output; `None` runs BaseRecalibrator only.
    pub output_bam: Option<PathBuf>,
    /// `@PG CL` of the `GATK ApplyBQSR` line; header metadata only.
    pub pg_command_line: String,
    pub threads: usize,
}

// ---------------------------------------------------------------------------
// Intervals and known sites

/// Merged 1-based closed intervals per contig, from a BED (§1.1).
struct Intervals(Vec<Vec<(i64, i64)>>);

impl Intervals {
    fn load(
        path: &Path,
        contig_index: &HashMap<&[u8], usize>,
        n_contigs: usize,
    ) -> io::Result<Intervals> {
        let mut per: Vec<Vec<(i64, i64)>> = vec![Vec::new(); n_contigs];
        for line in std::fs::read_to_string(path)?.lines() {
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let parse = |s: &str| {
                s.trim()
                    .parse::<i64>()
                    .map_err(|_| invalid(format!("bad BED line: {line}")))
            };
            let (Some(&chrom), Some(a), Some(b)) = (f.first(), f.get(1), f.get(2)) else {
                return Err(invalid(format!("bad BED line: {line}")));
            };
            let c = *contig_index
                .get(chrom.as_bytes())
                .ok_or_else(|| invalid(format!("BED contig {chrom} not in BAM header")))?;
            per[c].push((parse(a)? + 1, parse(b)?));
        }
        for iv in &mut per {
            iv.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::with_capacity(iv.len());
            for &(a, b) in iv.iter() {
                match merged.last_mut() {
                    Some(last) if a <= last.1 + 1 => last.1 = last.1.max(b),
                    _ => merged.push((a, b)),
                }
            }
            *iv = merged;
        }
        Ok(Intervals(per))
    }

    fn overlaps(&self, contig: usize, start: i64, end: i64) -> bool {
        let iv = &self.0[contig];
        let i = iv.partition_point(|&(_, e)| e < start);
        i < iv.len() && iv[i].0 <= end
    }
}

/// Known-site extents (1-based closed, VCF `POS..POS+len(REF)-1` or `INFO/END`;
/// FILTER ignored), sorted by start per contig, duplicates removed (§2.2).
struct KnownSites {
    sites: Vec<Vec<(u32, u32)>>,
    max_len: Vec<u32>,
}

/// Cache format (little-endian): magic, contig count, then per contig its name,
/// site count and `(start, end)` pairs; validated by the source files' sizes and
/// mtimes and by the contig names. Known sites are static, and parsing dbSNP
/// takes ~37 s while reading this takes ~1 s.
const CACHE_MAGIC: &[u8; 8] = b"NKSITES1";

fn cache_key(vcfs: &[PathBuf]) -> io::Result<String> {
    let mut key = String::new();
    for v in vcfs {
        let m = std::fs::metadata(v)?;
        let mtime = m
            .modified()?
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        key += &format!("{}\t{}\t{}\n", v.canonicalize()?.display(), m.len(), mtime);
    }
    Ok(key)
}

impl KnownSites {
    /// Load from `cache` if it matches the sources and contigs; otherwise parse
    /// the VCFs and (re)write the cache.
    fn load_cached(
        vcfs: &[PathBuf],
        cache: Option<&Path>,
        contigs: &[&[u8]],
        contig_index: &HashMap<&[u8], usize>,
        threads: usize,
    ) -> io::Result<KnownSites> {
        let key = cache_key(vcfs)?;
        if let Some(path) = cache
            && let Ok(bytes) = std::fs::read(path)
            && let Some(ks) = KnownSites::decode(&bytes, &key, contigs)
        {
            return Ok(ks);
        }
        let ks = KnownSites::load(vcfs, contig_index, contigs.len(), threads)?;
        if let Some(path) = cache {
            let tmp = path.with_extension("tmp");
            std::fs::write(&tmp, ks.encode(&key, contigs))?;
            std::fs::rename(&tmp, path)?;
        }
        Ok(ks)
    }

    fn encode(&self, key: &str, contigs: &[&[u8]]) -> Vec<u8> {
        let mut out = CACHE_MAGIC.to_vec();
        out.extend_from_slice(&(key.len() as u64).to_le_bytes());
        out.extend_from_slice(key.as_bytes());
        out.extend_from_slice(&(contigs.len() as u64).to_le_bytes());
        for (name, sites) in contigs.iter().zip(&self.sites) {
            out.extend_from_slice(&(name.len() as u64).to_le_bytes());
            out.extend_from_slice(name);
            out.extend_from_slice(&(sites.len() as u64).to_le_bytes());
            for &(a, b) in sites {
                out.extend_from_slice(&a.to_le_bytes());
                out.extend_from_slice(&b.to_le_bytes());
            }
        }
        out
    }

    fn decode(bytes: &[u8], key: &str, contigs: &[&[u8]]) -> Option<KnownSites> {
        let mut at = 0usize;
        let mut take = |n: usize| -> Option<&[u8]> {
            let s = bytes.get(at..at + n)?;
            at += n;
            Some(s)
        };
        let u64_at = |s: &[u8]| u64::from_le_bytes(s.try_into().unwrap()) as usize;
        if take(8)? != CACHE_MAGIC {
            return None;
        }
        let klen = u64_at(take(8)?);
        if take(klen)? != key.as_bytes() || u64_at(take(8)?) != contigs.len() {
            return None;
        }
        let mut sites = Vec::with_capacity(contigs.len());
        for name in contigs {
            let nlen = u64_at(take(8)?);
            if take(nlen)? != *name {
                return None;
            }
            let n = u64_at(take(8)?);
            let raw = take(n.checked_mul(8)?)?;
            sites.push(
                raw.as_chunks::<8>()
                    .0
                    .iter()
                    .map(|c| {
                        (
                            u32::from_le_bytes([c[0], c[1], c[2], c[3]]),
                            u32::from_le_bytes([c[4], c[5], c[6], c[7]]),
                        )
                    })
                    .collect::<Vec<_>>(),
            );
        }
        let max_len = sites
            .iter()
            .map(|s: &Vec<(u32, u32)>| s.iter().map(|(a, b)| b - a).max().unwrap_or(0))
            .collect();
        Some(KnownSites { sites, max_len })
    }

    fn load(
        vcfs: &[PathBuf],
        contig_index: &HashMap<&[u8], usize>,
        n_contigs: usize,
        threads: usize,
    ) -> io::Result<KnownSites> {
        let mut sites: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n_contigs];
        for vcf in vcfs {
            for_each_record(vcf, threads, &mut |r| {
                let Some(&c) = contig_index.get(r.chrom) else {
                    return Ok(());
                };
                let start = r.pos as u32;
                let end = if r.info.windows(4).any(|w| w == b"END=") {
                    r.info
                        .split(|&b| b == b';')
                        .find_map(|f| f.strip_prefix(b"END="))
                        .and_then(|v| std::str::from_utf8(v).ok()?.parse::<u32>().ok())
                        .unwrap_or(start + r.ref_allele.len() as u32 - 1)
                } else {
                    start + r.ref_allele.len() as u32 - 1
                };
                sites[c].push((start, end));
                Ok(())
            })?;
        }
        let mut max_len = vec![0u32; n_contigs];
        for (c, s) in sites.iter_mut().enumerate() {
            s.sort_unstable();
            s.dedup();
            max_len[c] = s.iter().map(|(a, b)| b - a).max().unwrap_or(0);
        }
        Ok(KnownSites { sites, max_len })
    }

    /// Sites overlapping `[start, end]` on `contig`.
    fn overlapping(
        &self,
        contig: usize,
        start: i64,
        end: i64,
    ) -> impl Iterator<Item = (i64, i64)> + '_ {
        let s = &self.sites[contig];
        let hi = s.partition_point(|&(a, _)| a as i64 <= end);
        let lo_start = start - self.max_len[contig] as i64;
        let lo = s.partition_point(|&(a, _)| (a as i64) < lo_start);
        s[lo..hi]
            .iter()
            .filter(move |&&(_, b)| b as i64 >= start)
            .map(|&(a, b)| (a as i64, b as i64))
    }
}

/// `@RG` identifiers (PU if present, else ID) in header order, and ID -> index (§3.1).
fn read_groups(header: &Header) -> (Vec<String>, HashMap<Vec<u8>, usize>) {
    let mut names = Vec::new();
    let mut by_id = HashMap::new();
    for line in String::from_utf8_lossy(&header.text)
        .lines()
        .filter(|l| l.starts_with("@RG"))
    {
        let field = |k: &str| {
            line.split('\t')
                .find_map(|f| f.strip_prefix(k))
                .map(str::to_owned)
        };
        let Some(id) = field("ID:") else { continue };
        let ident = field("PU:").unwrap_or_else(|| id.clone());
        let idx = match names.iter().position(|n| *n == ident) {
            Some(i) => i,
            None => {
                names.push(ident);
                names.len() - 1
            }
        };
        by_id.insert(id.into_bytes(), idx);
    }
    (names, by_id)
}

/// Run `f(k, range)` over record indices `[0, n)` split into contiguous
/// chunks, one per thread; `k < threads` is the chunk's index.
fn par_chunks<T: Send>(
    n: usize,
    threads: usize,
    f: impl Fn(usize, std::ops::Range<usize>) -> io::Result<T> + Sync,
) -> io::Result<Vec<T>> {
    let per = n.div_ceil(threads.max(1)).max(1);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n)
            .step_by(per)
            .enumerate()
            .map(|(k, s)| {
                scope.spawn({
                    let f = &f;
                    move || f(k, s..(s + per).min(n))
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("worker panicked"))
            .collect()
    })
}

/// Where the input records come from: a set in memory, or `opts.input`
/// streamed from disk (whole genomes).
enum Input<'a> {
    Set(&'a RecordSet),
    File(&'a Path),
}

/// A contiguous run of input records.
trait Records: Sync {
    fn len(&self) -> usize;
    fn record(&self, i: usize) -> &[u8];
}

struct SetSlice<'a>(&'a RecordSet, usize, usize);

impl Records for SetSlice<'_> {
    fn len(&self) -> usize {
        self.2 - self.1
    }
    fn record(&self, i: usize) -> &[u8] {
        self.0.record(self.1 + i)
    }
}

impl Records for Batch {
    fn len(&self) -> usize {
        Batch::len(self)
    }
    fn record(&self, i: usize) -> &[u8] {
        Batch::record(self, i)
    }
}

impl Input<'_> {
    fn header(&self) -> io::Result<Header> {
        match self {
            Input::Set(set) => Ok(set.header.clone()),
            Input::File(path) => {
                // Standalone BQSR may be fed any BAM; the pipeline checks its
                // own streams before it gets here.
                nomics_io::bamstream::warn_if_incomplete(path);
                Ok(BamStream::open(path, 1, 1, 1)?.header)
            }
        }
    }

    /// Call `f(first, records)` on consecutive batches of at most `size`
    /// records, in order; `first` is the index of the batch's first record.
    fn batches(
        &self,
        size: usize,
        mut f: impl FnMut(usize, &dyn Records) -> io::Result<()>,
    ) -> io::Result<()> {
        match self {
            Input::Set(set) => {
                for start in (0..set.len()).step_by(size) {
                    f(start, &SetSlice(set, start, (start + size).min(set.len())))?;
                }
            }
            Input::File(path) => {
                let mut stream = BamStream::open(path, size.min(STREAM_BATCH), 2, INFLATE_THREADS)?;
                let mut first = 0;
                while let Some(batch) = stream.next_batch()? {
                    f(first, &batch)?;
                    first += batch.len();
                }
            }
        }
        Ok(())
    }
}

/// Records per batch when streaming from disk.
const STREAM_BATCH: usize = 1 << 17;
/// Threads inflating the input when streaming from disk.
const INFLATE_THREADS: usize = 8;

/// BaseRecalibrator + ApplyBQSR streaming `opts.input` from disk and writing
/// `opts.output_bam` as it goes, for inputs too large to hold in memory.
pub fn run_streaming(opts: &Options) -> io::Result<()> {
    run_streaming_with(opts, |_| Ok(()))
}

/// [`run_streaming`], also handing `sink` the uncompressed stream written to
/// `opts.output_bam` as it goes: the header's bytes first, then
/// `block_size`-prefixed records, in file order. `sink` may take the
/// buffer's contents. `opts.output_bam` with the spill extension
/// (`nomics_io::spill`) writes the light format instead of BGZF.
pub fn run_streaming_with(
    opts: &Options,
    mut sink: impl FnMut(&mut Vec<u8>) -> io::Result<()>,
) -> io::Result<()> {
    run_on_impl(
        Input::File(&opts.input),
        opts,
        true,
        TableSource::Build,
        &mut sink,
    )
    .map(|_| ())
}

/// Whether a BaseRecalibrator report is one [`apply_from_report`] implements:
/// the five standard tables with exactly this port's `Arguments` values.
pub fn check_report(text: &str) -> Result<(), String> {
    tables::check_report_arguments(text)
}

/// ApplyBQSR alone (`gatk ApplyBQSR --bqsr-recal-file`): recalibrate
/// `opts.input` (streamed) with the report at `report`, parsed back as GATK's
/// `RecalibrationReport` does (§8.2), and write `opts.output_bam`.
/// `known_sites`, `reference` and `output_table` are not used.
pub fn apply_from_report(opts: &Options, report: &Path) -> io::Result<()> {
    if opts.scatter.is_some() || opts.output_bam.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "apply_from_report needs one interval group and output_bam",
        ));
    }
    let text = std::fs::read_to_string(report)?;
    run_on_impl(
        Input::File(&opts.input),
        opts,
        true,
        TableSource::Report(&text),
        &mut |_| Ok(()),
    )
    .map(|_| ())
}

/// Load `opts.input` and recalibrate it.
pub fn run(opts: &Options) -> io::Result<()> {
    let t = Instant::now();
    let bam = RecordSet::load(&opts.input, opts.threads.max(1))?;
    eprintln!(
        "bqsr stage load+decode BAM: {:.1}s",
        t.elapsed().as_secs_f64()
    );
    run_on(&bam, opts)
}

/// BaseRecalibrator + ApplyBQSR over records already in memory (`opts.input` unused).
pub fn run_on(bam: &RecordSet, opts: &Options) -> io::Result<()> {
    run_on_impl(Input::Set(bam), opts, true, TableSource::Build, &mut |_| {
        Ok(())
    })
    .map(|_| ())
}

/// BaseRecalibrator + ApplyBQSR in memory. The recalibration table is written,
/// but `output_bam` is left for the caller to write concurrently with its next stage.
pub fn run_on_set(bam: &RecordSet, opts: &Options) -> io::Result<RecordSet> {
    run_on_impl(
        Input::Set(bam),
        opts,
        false,
        TableSource::Build,
        &mut |_| Ok(()),
    )?
    .ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "run_on_set requires output_bam",
        )
    })
}

/// Interval groups, scatter orders and read groups for a run over a BAM
/// with header `in_header`.
struct Setup {
    contig_lens: Vec<u32>,
    groups: Vec<Intervals>,
    gather_order: Vec<usize>,
    merge_order: Vec<usize>,
    merge_seed: i64,
    rg_names: Vec<String>,
    rg_by_id: HashMap<Vec<u8>, usize>,
}

impl Setup {
    fn new(in_header: &Header, opts: &Options) -> io::Result<Setup> {
        let contig_lens: Vec<u32> = in_header.refs.iter().map(|r| r.len).collect();
        let contig_index = contig_index(in_header);
        let (groups, gather_order, merge_order, merge_seed) = match &opts.scatter {
            None => (
                vec![Intervals::load(
                    &opts.intervals,
                    &contig_index,
                    contig_lens.len(),
                )?],
                vec![0],
                vec![0],
                0,
            ),
            Some(sc) => (
                sc.groups
                    .iter()
                    .map(|p| Intervals::load(p, &contig_index, contig_lens.len()))
                    .collect::<io::Result<Vec<_>>>()?,
                sc.gather_order.clone(),
                sc.merge_order.clone(),
                sc.merge_seed,
            ),
        };
        let n_groups = groups.len();
        if n_groups > 64 {
            return Err(invalid(format!(
                "{n_groups} interval groups; at most 64 are supported"
            )));
        }
        for order in [&gather_order, &merge_order] {
            let mut seen = order.clone();
            seen.sort_unstable();
            if seen != (0..n_groups).collect::<Vec<_>>() {
                return Err(invalid(
                    "scatter order is not a permutation of the groups".into(),
                ));
            }
        }
        let (rg_names, rg_by_id) = read_groups(in_header);
        Ok(Setup {
            contig_lens,
            groups,
            gather_order,
            merge_order,
            merge_seed,
            rg_names,
            rg_by_id,
        })
    }

    /// The groups whose runs visit a record: it overlaps one of the group's
    /// intervals under htsjdk's rule and passes Wellformed (§1.2-1.3).
    fn visiting(&self, r: &Read) -> u64 {
        let Some((s, e)) = r.query_span() else {
            return 0;
        };
        if !r.is_wellformed(&self.contig_lens) {
            return 0;
        }
        self.groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.overlaps(r.ref_id as usize, s, e))
            .fold(0, |m, (g, _)| m | 1 << g)
    }
}

fn contig_index(header: &Header) -> HashMap<&[u8], usize> {
    header
        .refs
        .iter()
        .enumerate()
        .map(|(i, r)| (&r.name[..], i))
        .collect()
}

/// BaseRecalibrator's finished tables (already written to `output_table`),
/// for [`apply_streaming`].
pub struct BuiltTables {
    tables: Tables,
    quantized: Vec<u8>,
}

/// BaseRecalibrator (§2-§5) over records fed in input order, batch by batch:
/// the same tables whether the records are read back from the BAM or come
/// straight from the stage that writes it.
///
/// Every count is an integer (`Counts { obs: u64, err: u64 }`), with one
/// accumulator per worker kept across batches and summed at the end, so how
/// the records are split into batches and workers never changes the sums.
/// Floating point (`qual_summary`, `gathered_estimates`, `build_report`,
/// `Recalibration::new`) only ever sees the final merged integers.
pub struct TableBuilder {
    setup: Setup,
    known: KnownSites,
    fasta: Fasta,
    /// Reference contigs in BAM header order.
    reference: Vec<Option<usize>>,
    /// Per worker: the summed tables, and each group's quality rows
    /// (`[group][rg][q]`) for GatherBQSRReports' read-group estimates.
    accumulators: Vec<std::sync::Mutex<(Tables, Vec<Counts>)>>,
    /// Records added so far (for error messages).
    seen: usize,
    threads: usize,
    output_table: PathBuf,
}

/// The record bodies of a `block_size`-prefixed record stream.
struct StreamRecords<'a> {
    data: &'a [u8],
    /// Record `i` is `data[starts[i]..starts[i + 1] - 4]`; the last one ends
    /// at `data.len()`.
    starts: Vec<usize>,
}

impl Records for StreamRecords<'_> {
    fn len(&self) -> usize {
        self.starts.len()
    }
    fn record(&self, i: usize) -> &[u8] {
        let end = self.starts.get(i + 1).map_or(self.data.len(), |&s| s - 4);
        &self.data[self.starts[i]..end]
    }
}

impl TableBuilder {
    /// Load the known sites and the reference for a BAM with header `header`.
    pub fn new(header: &Header, opts: &Options) -> io::Result<TableBuilder> {
        let threads = opts.threads.max(1);
        let mut t = Instant::now();
        let setup = Setup::new(header, opts)?;
        let contig_index = contig_index(header);
        let contig_names: Vec<&[u8]> = header.refs.iter().map(|r| &r.name[..]).collect();
        let known = KnownSites::load_cached(
            &opts.known_sites,
            opts.known_sites_cache.as_deref(),
            &contig_names,
            &contig_index,
            threads,
        )?;
        lap("known sites", &mut t);
        let mut fasta = Fasta::load(&opts.reference, threads)?;
        fasta.normalize_upper_acgtn();
        let reference: Vec<Option<usize>> = header
            .refs
            .iter()
            .map(|r| fasta.index_of(&String::from_utf8_lossy(&r.name)))
            .collect();
        lap("reference", &mut t);
        let n = setup.groups.len() * setup.rg_names.len() * N_QUAL;
        let accumulators = (0..threads)
            .map(|_| {
                std::sync::Mutex::new((
                    Tables::new(setup.rg_names.clone()),
                    vec![Counts::default(); n],
                ))
            })
            .collect();
        Ok(TableBuilder {
            setup,
            known,
            fasta,
            reference,
            accumulators,
            seen: 0,
            threads,
            output_table: opts.output_table.clone(),
        })
    }

    /// Add the next records: `block_size`-prefixed BAM records, as they
    /// follow the header in the uncompressed BAM.
    pub fn add_stream(&mut self, data: &[u8]) -> io::Result<()> {
        let mut starts = Vec::new();
        let mut at = 0;
        while at < data.len() {
            let len = data
                .get(at..at + 4)
                .map(|b| i32::from_le_bytes(b.try_into().unwrap()))
                .and_then(|n| usize::try_from(n).ok())
                .filter(|&n| at + 4 + n <= data.len())
                .ok_or_else(|| {
                    invalid(format!("record {} truncated", self.seen + starts.len()))
                })?;
            starts.push(at + 4);
            at += 4 + len;
        }
        self.add(&StreamRecords { data, starts })
    }

    fn add(&mut self, bam: &dyn Records) -> io::Result<()> {
        let first = self.seen;
        let (setup, known, fasta, reference) =
            (&self.setup, &self.known, &self.fasta, &self.reference);
        let accumulators = &self.accumulators;
        let n_rg = setup.rg_names.len();
        par_chunks(bam.len(), self.threads, |k, range| -> io::Result<()> {
            let mut acc = accumulators[k].lock().expect("accumulator poisoned");
            let (tables, group_qual) = &mut *acc;
            for i in range {
                let Some(r) = Read::decode(bam.record(i)) else {
                    return Err(invalid(format!("record {} undecodable", first + i)));
                };
                let groups_visiting = setup.visiting(&r);
                if groups_visiting == 0 || !r.passes_recalibrator_filters(&setup.contig_lens) {
                    continue;
                }
                let Some(rg) = r.rg.and_then(|g| setup.rg_by_id.get(g)).copied() else {
                    return Err(invalid(format!("record {}: unknown read group", first + i)));
                };
                let Some(c) = clip_for_recalibration(&r) else {
                    continue;
                };
                let refseq = reference[r.ref_id as usize]
                    .map(|ci| fasta.contig(ci))
                    .ok_or_else(|| invalid(format!("contig {} missing from reference", r.ref_id)))?;
                let bases = &r.bases[c.lo..c.hi];
                let quals = &r.quals[c.lo..c.hi];
                let snp = snp_errors(bases, &c, refseq);
                let mask = known_site_mask(
                    &c,
                    known.overlapping(r.ref_id as usize, r.pos1, r.alignment_end()),
                );
                let ctx = context_keys(bases, quals, r.reverse(), LOW_QUAL_TAIL);
                let len = c.len();
                for k in 0..len {
                    if !is_regular_base(bases[k])
                        || quals[k] < PRESERVE_QSCORES_LESS_THAN
                        || mask[k]
                    {
                        continue;
                    }
                    let error = snp[k] == 1;
                    let q = quals[k] as usize;
                    tables.add_base(
                        rg,
                        q,
                        ctx[k],
                        cycle_key(k, len, r.flag),
                        error,
                        groups_visiting.count_ones() as u64,
                    );
                    let mut m = groups_visiting;
                    while m != 0 {
                        let g = m.trailing_zeros() as usize;
                        m &= m - 1;
                        let c = &mut group_qual[(g * n_rg + rg) * N_QUAL + q];
                        c.obs += 1;
                        c.err += error as u64;
                    }
                }
            }
            Ok(())
        })?;
        self.seen += bam.len();
        Ok(())
    }

    /// Records added so far.
    pub fn records(&self) -> usize {
        self.seen
    }

    /// Sum the workers' tables and write the report to `output_table`.
    pub fn finish(self) -> io::Result<BuiltTables> {
        let mut t = Instant::now();
        let n_rg = self.setup.rg_names.len();
        let mut parts = self
            .accumulators
            .into_iter()
            .map(|m| m.into_inner().expect("accumulator poisoned"));
        let (mut tables, mut group_qual) = parts.next().expect("at least one worker");
        for (p, gq) in parts {
            tables.merge(&p);
            for (a, b) in group_qual.iter_mut().zip(&gq) {
                a.obs += b.obs;
                a.err += b.err;
            }
        }
        if self.setup.groups.len() > 1 {
            tables.gathered_estimates = Some(gathered_estimates(
                &group_qual,
                n_rg,
                &self.setup.gather_order,
            ));
        }
        let report = build_report(&tables);
        std::fs::write(&self.output_table, &report.text)?;
        lap("report", &mut t);
        Ok(BuiltTables {
            tables,
            quantized: report.quantized,
        })
    }
}

fn lap(what: &str, t: &mut Instant) {
    eprintln!("bqsr stage {what}: {:.1}s", t.elapsed().as_secs_f64());
    *t = Instant::now();
}

/// Where ApplyBQSR's tables come from.
enum TableSource<'a> {
    /// BaseRecalibrator over the input first.
    Build,
    /// A report, parsed back as GATK's `RecalibrationReport` does.
    Report(&'a str),
    /// Built from the same records as they were written ([`TableBuilder`]).
    Built(BuiltTables),
}

/// ApplyBQSR alone over `opts.input` streamed from disk, with the tables a
/// [`TableBuilder`] built from the same records (it already wrote
/// `opts.output_table`); writes `opts.output_bam`.
pub fn apply_streaming(opts: &Options, built: BuiltTables) -> io::Result<()> {
    apply_streaming_with(opts, built, |_| Ok(()))
}

/// [`apply_streaming`], also handing `sink` the uncompressed output stream
/// as [`run_streaming_with`] does.
pub fn apply_streaming_with(
    opts: &Options,
    built: BuiltTables,
    mut sink: impl FnMut(&mut Vec<u8>) -> io::Result<()>,
) -> io::Result<()> {
    if opts.output_bam.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "apply_streaming needs output_bam",
        ));
    }
    run_on_impl(
        Input::File(&opts.input),
        opts,
        true,
        TableSource::Built(built),
        &mut sink,
    )
    .map(|_| ())
}

/// With `write_output`, `output_bam` is written and no set is returned for an
/// `Input::File`; for an `Input::Set` the set is returned either way.
/// `sink` sees what is streamed to `output_bam` (`Input::File` only).
fn run_on_impl(
    input: Input,
    opts: &Options,
    write_output: bool,
    source: TableSource,
    sink: &mut dyn FnMut(&mut Vec<u8>) -> io::Result<()>,
) -> io::Result<Option<RecordSet>> {
    let in_header = input.header()?;
    let threads = opts.threads.max(1);
    let mut t = Instant::now();
    let setup = Setup::new(&in_header, opts)?;
    let n_groups = setup.groups.len();
    let (merge_order, merge_seed) = (&setup.merge_order, setup.merge_seed);
    let rg_by_id = &setup.rg_by_id;
    let visiting = |r: &Read| setup.visiting(r);
    let (tables, quantized) = match source {
        TableSource::Report(text) => {
            let parsed = tables::parse_report(text, setup.rg_names.clone()).map_err(invalid)?;
            lap("parse recal table", &mut t);
            parsed
        }
        TableSource::Built(b) => (b.tables, b.quantized),
        TableSource::Build => {
            let mut builder = TableBuilder::new(&in_header, opts)?;
            t = Instant::now();
            // In memory the whole set is one batch.
            let pass1_batch = match input {
                Input::Set(_) => usize::MAX,
                Input::File(_) => STREAM_BATCH,
            };
            input.batches(pass1_batch, |_, bam| builder.add(bam))?;
            lap("recalibrator pass", &mut t);
            let built = builder.finish()?;
            t = Instant::now();
            (built.tables, built.quantized)
        }
    };

    // --- ApplyBQSR
    let Some(output_bam) = &opts.output_bam else {
        return Ok(None);
    };
    let recal = Recalibration::new(&tables, &quantized);
    let mut header = in_header.clone();
    let mut text = String::from_utf8_lossy(&header.text)
        .trim_end_matches('\0')
        .to_owned();
    if !text.ends_with('\n') && !text.is_empty() {
        text.push('\n');
    }
    text.push_str(&format!(
        "@PG\tID:GATK ApplyBQSR\tVN:4.6.2.0\tCL:{}\tPN:GATK ApplyBQSR\n",
        opts.pg_command_line
    ));
    // Every group's BAM has this header; merging translates their IDs.
    let merged =
        (n_groups > 1).then(|| merge::merge_headers(&vec![text.clone(); n_groups], merge_seed));
    header.text = match &merged {
        Some(m) => m.header_text.clone().into_bytes(),
        None => text.into_bytes(),
    };
    let mut stream = Vec::new();
    header.write(&mut stream)?;
    // Streaming input writes the output as it goes; in memory it is returned.
    let mut writer = match input {
        // recal.bam, or the same stream as a light spill file.
        Input::File(_) if write_output => Some(nomics_io::bamstream::StreamWriter::create(
            output_bam, threads,
        )?),
        Input::File(_) => {
            return Err(invalid("streaming BQSR must write output_bam".into()));
        }
        Input::Set(_) => None,
    };
    // Each group's ApplyBQSR writes the records its run visits, recalibrated
    // identically. samtools merge then interleaves the group BAMs by
    // (contig, position, strand), breaking ties by input position, so a run
    // of equal keys is emitted group by group in `merge_order`, each group's
    // records in input order. A record visited by two groups appears twice.
    // The workers write each copy already translated for its group, tagged
    // with the group's position in `merge_order`; the serial step only
    // orders a run's copies (a stable sort keeps input order within a group).
    // The open run is carried across batches, owned, since it may continue.
    let mut merge_pos = vec![0usize; n_groups];
    for (k, &g) in merge_order.iter().enumerate() {
        merge_pos[g] = k;
    }
    let mut carry: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut carry_key = u64::MAX;
    let flush = |run: &mut Vec<(usize, &[u8])>, stream: &mut Vec<u8>| {
        run.sort_by_key(|&(k, _)| k);
        for (_, rec) in run.iter() {
            stream.extend_from_slice(rec);
        }
        run.clear();
    };
    const BATCH: usize = 1 << 17;
    input.batches(BATCH, |start, bam| {
        let end = start + bam.len();
        type Part = (Vec<u8>, Vec<(u64, usize, usize, usize)>);
        let parts = par_chunks(end - start, threads, |_, range| -> io::Result<Part> {
            let mut out = Vec::new();
            let mut meta = Vec::new();
            let mut rec = Vec::new();
            for i in range {
                let raw = bam.record(i);
                let Some(r) = Read::decode(raw) else {
                    return Err(invalid(format!("record {} undecodable", start + i)));
                };
                let mask = visiting(&r);
                if mask == 0 {
                    continue;
                }
                let rg =
                    r.rg.and_then(|g| rg_by_id.get(g)).copied().ok_or_else(|| {
                        invalid(format!("record {}: unknown read group", start + i))
                    })?;
                let ctx = context_keys(&r.bases, r.quals, r.reverse(), LOW_QUAL_TAIL);
                let len = r.bases.len();
                rec.clear();
                rec.extend_from_slice(raw);
                for k in 0..len {
                    let q = r.quals[k];
                    if q < PRESERVE_QSCORES_LESS_THAN {
                        continue;
                    }
                    if let Some(nq) =
                        recal.recalibrate(rg, q as usize, ctx[k], cycle_key(k, len, r.flag))
                    {
                        rec[r.qual_offset + k] = nq;
                    }
                }
                // samtools' merge key: contig, 1-based position, strand.
                let key = (r.ref_id as u64) << 32 | (r.pos1 as u64) << 1 | r.reverse() as u64;
                match &merged {
                    None => {
                        let at = out.len();
                        out.extend_from_slice(&(rec.len() as i32).to_le_bytes());
                        out.extend_from_slice(&rec);
                        meta.push((key, 0, at, out.len()));
                    }
                    Some(m) => {
                        let mut bits = mask;
                        while bits != 0 {
                            let g = bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            let k = merge_pos[g];
                            let at = out.len();
                            merge::translate_record(&rec, &m.rg[k], &m.pg[k], &mut out);
                            meta.push((key, k, at, out.len()));
                        }
                    }
                }
            }
            Ok((out, meta))
        })?;
        let mut run: Vec<(usize, &[u8])> = carry.iter().map(|(k, b)| (*k, &b[..])).collect();
        let mut run_key = carry_key;
        for (buf, meta) in &parts {
            for &(key, k, a, b) in meta {
                if key != run_key {
                    flush(&mut run, &mut stream);
                    run_key = key;
                }
                run.push((k, &buf[a..b]));
            }
        }
        let open: Vec<(usize, Vec<u8>)> = run.iter().map(|(k, b)| (*k, b.to_vec())).collect();
        drop(run);
        carry = open;
        carry_key = run_key;
        if let Some(w) = &mut writer {
            w.write_all(&stream)?;
            sink(&mut stream)?;
            stream.clear();
        }
        Ok(())
    })?;
    let mut run: Vec<(usize, &[u8])> = carry.iter().map(|(k, b)| (*k, &b[..])).collect();
    flush(&mut run, &mut stream);
    if let Some(mut w) = writer {
        w.write_all(&stream)?;
        sink(&mut stream)?;
        w.finish()?;
        lap("apply pass", &mut t);
        return Ok(None);
    }
    let set = RecordSet::from_stream(stream)?;
    if write_output {
        set.write(output_bam, threads)?;
    }
    lap("apply pass", &mut t);
    Ok(Some(set))
}
