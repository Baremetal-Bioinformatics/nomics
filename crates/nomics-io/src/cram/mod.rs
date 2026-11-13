//! CRAM 3.0 writer (hts-specs `CRAMv3.tex`, codecs from `CRAMcodecs.tex`),
//! with its `.crai` index.
//!
//! **Gate: decoded content, not bytes.** htslib's CRAM bytes depend on
//! per-version codec heuristics and container layout, so byte identity is out
//! of scope. A [`Profile`] says what the CRAM must decode to under htslib,
//! htsjdk 4 and htsjdk 5 (`docs/GATES.md`):
//!
//! - `SarekMd` / `SarekRecal`: what sarek 3.10.0's md.cram (samtools 1.21
//!   `view -C`) and recal.cram (GATK ApplyBQSR, htsjdk 4.2.0) decode to. The
//!   tag modes [`TagMode::Htslib`] and [`TagMode::HtsjdkSorted`] store tags
//!   as those writers do (MD/NM dropped where htslib would regenerate them,
//!   RG as a data series, htsjdk's sorted tag order), and the header gets the
//!   same `M5`/`UR` additions plus our own `@PG` line.
//! - `Lossless`: the BAM's own records, stricter than htslib's writer:
//!   every tag stored verbatim in BAM order with its BAM type code, and
//!   nothing regenerated on decode:
//!
//! - `RG` stays a tag; the RG data series is constant -1 (no read group).
//! - `MD`/`NM` stay tags. A mapped read lacking one gets `cF:C` (htslib's CRAM
//!   flags tag; htslib and htsjdk 5 strip it and then do not invent the
//!   tag). htsjdk 4 would show `cF`; BWA output never needs it.
//! - Mate fields are attached only when htslib, htsjdk 4 and htsjdk 5 all
//!   rebuild them exactly (see `container::attachable`; the three disagree
//!   on TLEN ties); otherwise they are stored verbatim (detached).
//!
//! Layout: one slice per container, `records_per_slice` records (10,000, as
//! htslib), a new container at every reference change, single-reference
//! slices only (no multi-ref containers: small contigs each get their own),
//! reference-based sequence with slice MD5s, no embedded reference. Trial
//! containers compress each block with whichever of raw, rANS 4x8 order 0/1
//! and gzip is smallest (gzip probed with our fast deflater, then zlib level 5
//! when in contention); the containers after a trial reuse its choice per
//! block content id ([`TRIAL_EVERY`], [`Role`]). Containers are encoded in
//! parallel; the output does not depend on the thread count.

mod container;
pub mod itf8;
pub mod md5;
pub mod rans;

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::{Arc, Mutex};

use crate::bam::Header;
use crate::deflate::FastDeflater;
use crate::fasta::Fasta;
pub use container::TagMode;
use container::{ContainerInput, EncodeParams, Encoded, Encoder, SliceRef};

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// The largest raw BAM record (without `block_size`) the writer accepts.
pub const MAX_RECORD_BYTES: usize = 256 << 20;

#[derive(Debug, Clone)]
pub struct CramOptions {
    /// Records per slice (and container). htslib's default is 10,000.
    pub records_per_slice: usize,
    /// A slice also ends before its raw BAM records would pass this many
    /// bytes (64 MiB by default, at most [`MAX_RECORD_BYTES`]). A single
    /// record larger than this forms a slice of its own; records larger than
    /// [`MAX_RECORD_BYTES`] are rejected. A slice is thus under 512 MiB of
    /// raw records, which keeps every block and container size well inside
    /// CRAM's 32-bit fields (all of which are also checked when written).
    pub max_slice_bytes: usize,
    /// Encoding threads (containers in flight). At least 1.
    pub threads: usize,
    /// The 20-byte file id in the file definition (htslib writes the file name).
    pub file_id: [u8; 20],
    /// How tags are stored ([`TagMode::Lossless`] by default).
    pub tags: TagMode,
    /// `UR:` to add (with `M5:`) to `@SQ` lines lacking one, as htslib does
    /// with the reference path it was given. `None`: M5 only.
    pub ur: Option<String>,
}

impl Default for CramOptions {
    fn default() -> Self {
        CramOptions {
            records_per_slice: 10_000,
            max_slice_bytes: 64 << 20,
            threads: 1,
            file_id: [0; 20],
            tags: TagMode::Lossless,
            ur: None,
        }
    }
}

impl CramOptions {
    /// `file_id` from a file name, truncated or NUL-padded to 20 bytes.
    pub fn with_file_name(mut self, name: &str) -> Self {
        self.file_id = [0; 20];
        for (d, s) in self.file_id.iter_mut().zip(name.bytes()) {
            *d = s;
        }
        self
    }
}

/// One `.crai` line: a slice.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CraiEntry {
    pub ref_id: i32,
    pub start: i64,
    pub span: i64,
    pub container_offset: u64,
    pub slice_offset: u64,
    pub slice_size: u64,
}

/// What a finished CRAM looks like.
#[derive(Debug, Clone, Default)]
pub struct CramSummary {
    pub records: u64,
    pub containers: u64,
    pub bytes: u64,
    pub index: Vec<CraiEntry>,
    /// Most finished containers ever held waiting for an earlier one
    /// (at most [`reorder_window`]).
    pub reorder_peak: usize,
}

/// How far past the next container to write a worker may start one: a slow
/// container holds back at most this many finished ones in memory.
pub fn reorder_window(threads: usize) -> u64 {
    2 * threads.max(1) as u64 + 2
}

/// The containers written so far, for workers keeping within
/// [`reorder_window`] of them; `closed` once the writer has stopped.
#[derive(Default)]
struct Written {
    state: Mutex<(u64, bool)>,
    moved: std::sync::Condvar,
}

impl Written {
    /// Wait until container `idx` is within `window` of the next to write
    /// (or the writer has stopped).
    fn wait_for(&self, idx: u64, window: u64) {
        let mut st = self.state.lock().unwrap();
        while idx >= st.0 + window && !st.1 {
            st = self.moved.wait(st).unwrap();
        }
    }

    fn set(&self, written: u64) {
        self.state.lock().unwrap().0 = written;
        self.moved.notify_all();
    }
}

/// Marks [`Written`] closed when the writer stops, however it stops.
struct CloseOnDrop<'a>(&'a Written);

impl Drop for CloseOnDrop<'_> {
    fn drop(&mut self) {
        if let Ok(mut st) = self.0.state.lock() {
            st.1 = true;
        }
        self.0.moved.notify_all();
    }
}

/// `.crai` text (one tab-separated line per slice), gzip-compressed.
pub fn crai_bytes(index: &[CraiEntry]) -> Vec<u8> {
    let mut text = Vec::new();
    for e in index {
        // As `samtools index` writes: the slice header's start and span.
        let (start, span) = (e.start, e.span);
        text.extend_from_slice(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                e.ref_id, start, span, e.container_offset, e.slice_offset, e.slice_size
            )
            .as_bytes(),
        );
    }
    container::gzip(&mut FastDeflater::new(), &text)
}

/// Add `M5:` (MD5 of the upper-cased contig) to every `@SQ` line that has
/// none, checking `LN` against the FASTA. Contigs missing from the FASTA are
/// left without `M5` (records on them are rejected at encode time).
pub fn header_with_m5(
    text: &[u8],
    fasta: &Fasta,
    threads: usize,
    ur: Option<&str>,
) -> io::Result<Vec<u8>> {
    let text = match text.iter().position(|&b| b == 0) {
        Some(n) => &text[..n], // BAM l_text may carry NUL padding
        None => text,
    };
    let lines: Vec<&[u8]> = text.split_inclusive(|&b| b == b'\n').collect();
    // (line index, contig index) for lines needing an M5.
    let mut todo: Vec<(usize, usize)> = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        if !line.starts_with(b"@SQ\t") {
            continue;
        }
        let body = line.strip_suffix(b"\n").unwrap_or(line);
        let fields: Vec<&[u8]> = body.split(|&b| b == b'\t').collect();
        if fields.iter().any(|f| f.starts_with(b"M5:")) {
            continue;
        }
        let Some(sn) = fields.iter().find_map(|f| f.strip_prefix(b"SN:")) else {
            return Err(invalid("@SQ line without SN"));
        };
        let sn = String::from_utf8_lossy(sn);
        let Some(ci) = fasta.index_of(&sn) else {
            continue;
        };
        if let Some(ln) = fields.iter().find_map(|f| f.strip_prefix(b"LN:")) {
            let ln: u64 = std::str::from_utf8(ln)
                .ok()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| invalid(format!("@SQ {sn}: bad LN")))?;
            if ln != fasta.index[ci].len {
                return Err(invalid(format!(
                    "@SQ {sn}: LN {ln} but the FASTA contig has {} bases",
                    fasta.index[ci].len
                )));
            }
        }
        todo.push((i, ci));
    }
    // Hash contigs in parallel, largest first.
    let mut order: Vec<usize> = (0..todo.len()).collect();
    order.sort_by_key(|&k| std::cmp::Reverse(fasta.index[todo[k].1].len));
    let digests: Vec<Mutex<[u8; 16]>> = (0..todo.len()).map(|_| Mutex::new([0; 16])).collect();
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1) {
            scope.spawn(|| {
                loop {
                    let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(&j) = order.get(k) else { break };
                    let mut m = md5::Md5::new();
                    m.update_upper(fasta.contig(todo[j].1));
                    *digests[j].lock().unwrap() = m.finish();
                }
            });
        }
    });
    let mut m5: Vec<Option<[u8; 16]>> = vec![None; lines.len()];
    for (k, &(i, _)) in todo.iter().enumerate() {
        m5[i] = Some(*digests[k].lock().unwrap());
    }
    let mut out = Vec::with_capacity(text.len() + 40 * todo.len());
    for (i, line) in lines.iter().enumerate() {
        match m5[i] {
            Some(d) => {
                let (body, nl) = match line.strip_suffix(b"\n") {
                    Some(b) => (b, true),
                    None => (*line, false),
                };
                out.extend_from_slice(body);
                out.extend_from_slice(b"\tM5:");
                out.extend_from_slice(md5::hex(&d).as_bytes());
                if let Some(u) = ur.filter(|_| !body.windows(4).any(|w| w == b"\tUR:")) {
                    out.extend_from_slice(b"\tUR:");
                    out.extend_from_slice(u.as_bytes());
                }
                if nl {
                    out.push(b'\n');
                }
            }
            None => out.extend_from_slice(line),
        }
    }
    Ok(out)
}

/// `@RG` IDs in header order (htslib's and htsjdk's RG series index).
fn read_group_ids(text: &[u8]) -> Vec<Vec<u8>> {
    text.split(|&b| b == b'\n')
        .filter(|l| l.starts_with(b"@RG\t"))
        .filter_map(|l| {
            l.split(|&b| b == b'\t')
                .find_map(|f| f.strip_prefix(b"ID:"))
        })
        .map(|id| id.to_vec())
        .collect()
}

fn field<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    line.split(|&b| b == b'\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key))
}

/// Append `@PG ID:<name>[.N] PN:<name> VN:<version> CL:<command line>` after
/// the last `@PG` line (or at the end), the ID made unique with `.1`, `.2`,
/// ... : this program's own record of the CRAM step.
pub fn add_pg(text: &[u8], name: &str, version: &str, cl: &str) -> Vec<u8> {
    let text = text.split(|&b| b == 0).next().unwrap_or(&[]);
    let lines: Vec<&[u8]> = text.split_inclusive(|&b| b == b'\n').collect();
    let ids: std::collections::HashSet<&[u8]> = lines
        .iter()
        .filter(|l| l.starts_with(b"@PG\t"))
        .filter_map(|l| field(l.strip_suffix(b"\n").unwrap_or(l), b"ID:"))
        .collect();
    let mut id = name.to_owned();
    let mut n = 1;
    while ids.contains(id.as_bytes()) {
        id = format!("{name}.{n}");
        n += 1;
    }
    let line = format!("@PG\tID:{id}\tPN:{name}\tVN:{version}\tCL:{cl}\n");
    let last_pg = lines.iter().rposition(|l| l.starts_with(b"@PG\t"));
    let mut out = Vec::with_capacity(text.len() + line.len() + 1);
    for (i, l) in lines.iter().enumerate() {
        out.extend_from_slice(l);
        if Some(i) == last_pg {
            if !out.ends_with(b"\n") {
                out.push(b'\n');
            }
            out.extend_from_slice(line.as_bytes());
        }
    }
    if last_pg.is_none() {
        if !out.is_empty() && !out.ends_with(b"\n") {
            out.push(b'\n');
        }
        out.extend_from_slice(line.as_bytes());
    }
    out
}

/// What a CRAM should decode to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Profile {
    /// The BAM's records and header (plus `M5`).
    Lossless,
    /// sarek 3.10.0's `<prefix>.md.cram` (`samtools view -Ch`, htslib 1.21).
    SarekMd,
    /// sarek 3.10.0's `<prefix>.recal.cram` (GATK 4.6.2.0 ApplyBQSR writing
    /// CRAM from md.cram, htsjdk 4.2.0; one interval group, so no merge).
    SarekRecal,
    /// sarek's recal.cram with more than one interval group: the groups'
    /// ApplyBQSR CRAMs merged by `samtools merge --output-fmt
    /// cram,version=3.0` (samtools 1.24).
    SarekRecalMerged,
}

impl Profile {
    pub fn parse(s: &str) -> Option<Profile> {
        Some(match s {
            "lossless" => Profile::Lossless,
            "sarek-md" => Profile::SarekMd,
            "sarek-recal" => Profile::SarekRecal,
            "sarek-recal-merged" => Profile::SarekRecalMerged,
            _ => return None,
        })
    }

    /// The header to write and the options to write it with. `reference`
    /// is the FASTA (`UR:` gets its absolute path). `pg`: (version, command
    /// line) for this program's own `@PG` line on the sarek profiles; their
    /// reference CRAMs carry samtools'/ApplyBQSR's `@PG` lines, which the
    /// gate does not compare (`docs/GATES.md`), and we do not write lines
    /// for tools we did not run.
    pub fn prepare(
        self,
        header: &Header,
        reference: &std::path::Path,
        pg: Option<(&str, &str)>,
        mut opts: CramOptions,
    ) -> io::Result<(Header, CramOptions)> {
        let mut h = header.clone();
        match self {
            Profile::Lossless => opts.tags = TagMode::Lossless,
            Profile::SarekMd | Profile::SarekRecal | Profile::SarekRecalMerged => {
                opts.ur = Some(
                    std::path::absolute(reference)?
                        .to_string_lossy()
                        .into_owned(),
                );
                opts.tags = match self {
                    Profile::SarekMd => TagMode::Htslib,
                    Profile::SarekRecalMerged => TagMode::HtsjdkSortedMerged,
                    _ => TagMode::HtsjdkSorted,
                };
                if let Some((version, cl)) = pg {
                    h.text = add_pg(&header.text, "nomics-cram", version, cl);
                }
            }
        }
        Ok((h, opts))
    }
}

/// Check that the text `@SQ` lines list exactly the binary references, in
/// order: CRAM readers rebuild the reference list from the text alone.
fn check_sq_lines(header: &Header, text: &[u8]) -> io::Result<()> {
    let names: Vec<&[u8]> = text
        .split(|&b| b == b'\n')
        .filter(|l| l.starts_with(b"@SQ\t"))
        .filter_map(|l| {
            l.split(|&b| b == b'\t')
                .find_map(|f| f.strip_prefix(b"SN:"))
        })
        .collect();
    if names.len() != header.refs.len()
        || names
            .iter()
            .zip(&header.refs)
            .any(|(n, r)| **n != r.name[..])
    {
        return Err(invalid(
            "header @SQ lines do not match the BAM reference list; CRAM keeps only the text",
        ));
    }
    Ok(())
}

/// Codec trials (htslib's approach): a *trial* container tries every codec
/// on every block; the containers after it compress each block with the
/// codec a trial chose for that block's content id. Roles are assigned by
/// [`RecordSink`] from the container sequence alone (indices and record
/// counts), so the output does not depend on thread count or timing.
///
/// A container becomes a trial when it holds at least [`TRIAL_MIN_RECORDS`]
/// records (a few reads on a small contig are no guide) and no trial was
/// made in the last [`TRIAL_EVERY`] containers. Other containers use the
/// latest trial handed out at least [`TRIAL_LAG`] containers before them,
/// which is almost always finished by then; before there is one, they make
/// a full trial of their own.
const TRIAL_EVERY: u64 = 32;
const TRIAL_LAG: u64 = 16;
const TRIAL_MIN_RECORDS: usize = 1000;

/// How a container picks its codecs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Role {
    /// Try every codec, and publish the choices for later containers.
    Trial,
    /// Try every codec (no trial to follow yet).
    Full,
    /// Use the choices of the trial container with this index.
    Follow(u64),
}

/// The codec plans trial containers publish, for the containers after them.
#[derive(Default)]
struct Plans {
    done: Mutex<BTreeMap<u64, Option<Arc<container::CodecPlan>>>>,
    published: std::sync::Condvar,
}

impl Plans {
    /// The plan a container with `role` compresses with (`None`: a full
    /// trial), waiting for the trial container it follows.
    fn for_role(&self, role: Role) -> Option<Arc<container::CodecPlan>> {
        let Role::Follow(trial) = role else {
            return None;
        };
        let mut done = self.done.lock().unwrap();
        loop {
            if let Some(plan) = done.get(&trial) {
                // A failed trial (the file fails anyway) leaves full trials.
                return plan.clone();
            }
            done = self.published.wait(done).unwrap();
        }
    }

    /// Record what trial container `idx` chose (`None`: it failed). Blocks
    /// whose data was too small to judge by are left out of the plan, so
    /// their followers try every codec.
    fn publish(&self, idx: u64, codecs: Option<&[(i32, container::Codec, usize)]>) {
        let plan = codecs.map(|c| {
            Arc::new(
                c.iter()
                    .filter(|&&(_, _, len)| len >= 256)
                    .map(|&(id, codec, _)| (id, codec))
                    .collect(),
            )
        });
        // Kept for the whole file (a few KB per trial): a container's
        // worker may reach its lookup late, after later trials.
        self.done.lock().unwrap().insert(idx, plan);
        self.published.notify_all();
    }
}

/// Pushes records into a CRAM being written by [`write_cram`].
pub struct RecordSink<'s> {
    tx: Option<SyncSender<(u64, Role, ContainerInput)>>,
    /// Indices of the trial containers handed out so far.
    trials: Vec<u64>,
    cur: ContainerInput,
    per_slice: usize,
    max_bytes: usize,
    next_index: u64,
    records: u64,
    _scope: std::marker::PhantomData<&'s ()>,
}

impl RecordSink<'_> {
    /// Add one raw BAM record (without its `block_size`).
    pub fn push(&mut self, raw: &[u8]) -> io::Result<()> {
        if raw.len() < 32 {
            return Err(invalid("BAM record shorter than its fixed fields"));
        }
        let ref_id = i32::from_le_bytes(raw[0..4].try_into().unwrap());
        if raw.len() > MAX_RECORD_BYTES {
            return Err(invalid(format!(
                "BAM record of {} bytes is over the {} MiB CRAM record limit",
                raw.len(),
                MAX_RECORD_BYTES >> 20
            )));
        }
        if !self.cur.is_empty()
            && (ref_id != self.cur.ref_id
                || self.cur.len() >= self.per_slice
                || self.cur.data.len() + raw.len() > self.max_bytes)
        {
            self.flush()?;
        }
        if self.cur.is_empty() {
            self.cur.ref_id = ref_id;
            self.cur.record_counter = self.records;
        }
        self.cur.data.extend_from_slice(raw);
        self.cur.ends.push(self.cur.data.len());
        self.records += 1;
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.cur.is_empty() {
            return Ok(());
        }
        let cap = self.cur.data.len();
        let next = ContainerInput {
            data: Vec::with_capacity(cap),
            ends: Vec::with_capacity(self.per_slice),
            ref_id: 0,
            record_counter: 0,
        };
        let full = std::mem::replace(&mut self.cur, next);
        let tx = self.tx.as_ref().expect("sink used after finish");
        let idx = self.next_index;
        let since_trial = self.trials.last().map(|&t| idx - t);
        let role = if full.len() >= TRIAL_MIN_RECORDS && since_trial.is_none_or(|d| d >= TRIAL_EVERY) {
            self.trials.push(idx);
            Role::Trial
        } else {
            match self.trials.iter().rev().find(|&&t| idx - t >= TRIAL_LAG) {
                Some(&t) => Role::Follow(t),
                None => Role::Full,
            }
        };
        tx.send((idx, role, full))
            .map_err(|_| io::Error::other("CRAM encoder stopped early"))?;
        self.next_index += 1;
        Ok(())
    }
}

/// Write a coordinate- (or otherwise-) ordered stream of BAM records as CRAM
/// 3.0 to `out`. `feed` pushes the records, in order, into the sink; the
/// header is `header` with `M5:` added (see [`header_with_m5`]).
pub fn write_cram<W, F>(
    out: W,
    header: &Header,
    fasta: &Fasta,
    opts: &CramOptions,
    feed: F,
) -> io::Result<(W, CramSummary)>
where
    W: Write + Send,
    F: FnOnce(&mut RecordSink) -> io::Result<()>,
{
    let threads = opts.threads.max(1);
    let text = header_with_m5(&header.text, fasta, threads, opts.ur.as_deref())?;
    let params = EncodeParams {
        tags: opts.tags,
        read_groups: read_group_ids(&text),
    };
    let params = &params;
    check_sq_lines(header, &text)?;
    // BAM reference id -> FASTA contig.
    let contigs: Vec<Option<usize>> = header
        .refs
        .iter()
        .map(|r| fasta.index_of(&String::from_utf8_lossy(&r.name)))
        .collect();
    let contigs = &contigs;

    let mut prefix = Vec::with_capacity(text.len() / 4 + 128);
    prefix.extend_from_slice(b"CRAM\x03\x00");
    prefix.extend_from_slice(&opts.file_id);
    prefix.extend_from_slice(&container::header_container(
        &text,
        &mut FastDeflater::new(),
    )?);

    let plans = Plans::default();
    let written = Written::default();
    let window = reorder_window(threads);
    let (in_tx, in_rx) = sync_channel::<(u64, Role, ContainerInput)>(2 * threads);
    let in_rx: Arc<Mutex<Receiver<(u64, Role, ContainerInput)>>> = Arc::new(Mutex::new(in_rx));
    let (out_tx, out_rx) = sync_channel::<(u64, io::Result<Encoded>)>(2 * threads);

    std::thread::scope(|scope| {
        for _ in 0..threads {
            let rx = Arc::clone(&in_rx);
            let tx = out_tx.clone();
            let plans = &plans;
            let written = &written;
            scope.spawn(move || {
                let mut enc = Encoder::new();
                loop {
                    let job = rx.lock().unwrap().recv();
                    let Ok((idx, role, input)) = job else { break };
                    // Within the reorder window of the writer (the container it
                    // waits for is never held back, so this cannot deadlock).
                    written.wait_for(idx, window);
                    let plan = plans.for_role(role);
                    let bases = if input.ref_id >= 0 {
                        contigs
                            .get(input.ref_id as usize)
                            .copied()
                            .flatten()
                            .map(|c| fasta.contig(c))
                    } else {
                        None
                    };
                    // A panic in the encoder is a bug, but it must surface as
                    // this file's error, not abort the caller's process.
                    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        enc.encode(&input, &SliceRef { bases }, params, plan.as_deref())
                    }))
                    .unwrap_or_else(|p| {
                        let msg = p
                            .downcast_ref::<&str>()
                            .map(|s| s.to_string())
                            .or_else(|| p.downcast_ref::<String>().cloned())
                            .unwrap_or_default();
                        Err(io::Error::other(format!("CRAM encoder panicked: {msg}")))
                    });
                    if role == Role::Trial {
                        plans.publish(idx, r.as_ref().ok().map(|e| &e.codecs[..]));
                    }
                    let failed = r.is_err();
                    if tx.send((idx, r)).is_err() || failed {
                        break;
                    }
                }
            });
        }
        drop(out_tx);
        // Only the workers may hold the receiver, so that a sink blocked on a
        // full queue fails once they have all stopped.
        drop(in_rx);

        let written = &written;
        let writer = scope.spawn(move || -> io::Result<(W, CramSummary)> {
            let _close = CloseOnDrop(written);
            let mut out = out;
            let mut summary = CramSummary::default();
            out.write_all(&prefix)?;
            let mut offset = prefix.len() as u64;
            let mut pending: BTreeMap<u64, Encoded> = BTreeMap::new();
            let mut want = 0u64;
            for (idx, r) in out_rx.iter() {
                pending.insert(idx, r?);
                summary.reorder_peak = summary.reorder_peak.max(pending.len());
                while let Some(e) = pending.remove(&want) {
                    out.write_all(&e.bytes)?;
                    summary.index.push(CraiEntry {
                        ref_id: e.ref_id,
                        start: e.start,
                        span: e.span,
                        container_offset: offset,
                        slice_offset: e.landmark as u64,
                        slice_size: e.slice_size as u64,
                    });
                    offset += e.bytes.len() as u64;
                    summary.records += e.records as u64;
                    summary.containers += 1;
                    want += 1;
                    written.set(want);
                }
            }
            if !pending.is_empty() {
                return Err(io::Error::other("CRAM encoder lost a container"));
            }
            out.write_all(&container::EOF_CONTAINER)?;
            offset += container::EOF_CONTAINER.len() as u64;
            summary.bytes = offset;
            Ok((out, summary))
        });

        let mut sink = RecordSink {
            tx: Some(in_tx),
            cur: ContainerInput {
                data: Vec::new(),
                ends: Vec::new(),
                ref_id: 0,
                record_counter: 0,
            },
            per_slice: opts.records_per_slice.max(1),
            max_bytes: opts.max_slice_bytes.clamp(1, MAX_RECORD_BYTES),
            next_index: 0,
            trials: Vec::new(),
            records: 0,
            _scope: std::marker::PhantomData,
        };
        let fed = feed(&mut sink).and_then(|()| sink.flush());
        let pushed = sink.records;
        drop(sink.tx.take());
        drop(sink);
        let written = writer
            .join()
            .map_err(|_| io::Error::other("CRAM writer thread panicked"))?;
        // An encoder error is the root cause of a failed feed; report it first.
        let (w, summary) = written?;
        fed?;
        if summary.records != pushed {
            return Err(io::Error::other("CRAM record count mismatch"));
        }
        Ok((w, summary))
    })
}
