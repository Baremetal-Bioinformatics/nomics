//! BAM geometry, including the information lost by Salmon's intermediate RAD.
use std::collections::BTreeMap;
use std::io;
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, channel, sync_channel};
use std::time::{Duration, Instant};

use nomics_io::bam::{Header, RecordView, aux_value_len};
use nomics_io::bamstream::SliceBamReader;
use nomics_io::bgzf::ParallelBgzfReader;
use nomics_mem::HostArray;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    Single,
    Left,
    Right,
    Pair,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Hit {
    pub tid: usize,
    pub pos: usize,
    pub forward: bool,
    pub mate_forward: bool,
    pub length: usize,
    pub score: i32,
    pub status: Status,
}

// Grow in 2 MiB HostBuf mappings, without reallocating/copying an ever larger
// placement array. Each fragment is a hit count followed by two u64s per hit:
// tid/pos (u32 each), AS (i32), RAD length (u16), status and strand bits.
// A 20M-fragment, single-placement sample uses ~480 MB, independent of BAM
// record size. Small Vecs only index the mappings or hold one fragment.
const STORE_WORDS: usize = 1 << 18;

#[derive(Default)]
pub(crate) struct Fragments {
    chunks: Vec<HostArray<u64>>,
    words: usize,
    /// Word offsets where appended batches begin; each is a fragment boundary.
    marks: Vec<usize>,
}

/// Encode one fragment's RAD geometry, only after the caller has tallied the
/// original (unsaturated, per-hit status) geometry for FLD and detection.
fn encode(hits: &mut [Hit], words: &mut Vec<u64>) {
    rad_geometry(hits);
    words.push(hits.len() as u64);
    for hit in hits {
        words.push(hit.tid as u64 | ((hit.pos as u64) << 32));
        let status = match hit.status {
            Status::Single => 0,
            Status::Left => 1,
            Status::Right => 2,
            Status::Pair => 3,
        };
        words.push(
            u64::from(hit.score as u32)
                | ((hit.length as u64) << 32)
                | (status << 48)
                | (u64::from(hit.forward) << 50)
                | (u64::from(hit.mate_forward) << 51),
        );
    }
}

impl Fragments {
    /// Append whole encoded fragments, marking a fragment boundary.
    fn extend(&mut self, mut words: &[u64]) -> io::Result<()> {
        if words.is_empty() {
            return Ok(());
        }
        self.marks.push(self.words);
        while !words.is_empty() {
            let offset = self.words % STORE_WORDS;
            if offset == 0 {
                self.chunks.push(HostArray::zeroed(STORE_WORDS)?);
            }
            let n = words.len().min(STORE_WORDS - offset);
            self.chunks.last_mut().unwrap()[offset..offset + n].copy_from_slice(&words[..n]);
            self.words += n;
            words = &words[n..];
        }
        Ok(())
    }

    #[cfg(test)]
    pub fn push(&mut self, hits: &mut [Hit]) -> io::Result<()> {
        let mut words = Vec::new();
        encode(hits, &mut words);
        self.extend(&words)
    }

    /// Word ranges that each start and end on fragment boundaries, in order.
    pub fn ranges(&self) -> Vec<std::ops::Range<usize>> {
        let mut ends = self.marks.iter().skip(1).copied().chain([self.words]);
        self.marks
            .iter()
            .map(|&start| start..ends.next().unwrap())
            .collect()
    }

    #[cfg(test)]
    pub fn visit(&self, visit: impl FnMut(&[Hit])) {
        self.visit_range(0..self.words, visit);
    }

    pub fn visit_range(&self, range: std::ops::Range<usize>, mut visit: impl FnMut(&[Hit])) {
        let mut words = range.map(|i| self.chunks[i / STORE_WORDS][i % STORE_WORDS]);
        let mut hits = Vec::new();
        while let Some(count) = words.next() {
            hits.clear();
            for _ in 0..count {
                let a = words.next().unwrap();
                let b = words.next().unwrap();
                hits.push(Hit {
                    tid: a as u32 as usize,
                    pos: (a >> 32) as usize,
                    score: b as u32 as i32,
                    length: ((b >> 32) & 0xffff) as usize,
                    status: match (b >> 48) & 3 {
                        0 => Status::Single,
                        1 => Status::Left,
                        2 => Status::Right,
                        _ => Status::Pair,
                    },
                    forward: b & (1 << 50) != 0,
                    mate_forward: b & (1 << 51) != 0,
                });
            }
            visit(&hits);
        }
    }
}

impl Hit {
    pub fn format(&self) -> usize {
        if self.status != Status::Pair {
            return if self.forward { 10 } else { 11 };
        }
        let base = if self.forward != self.mate_forward {
            1
        } else {
            7
        };
        base + usize::from(!self.forward)
    }
}

#[derive(Debug)]
struct Alignment {
    tid: usize,
    pos: i32,
    mate_tid: i32,
    mate_pos: i32,
    flags: u16,
    length: usize,
    span: usize,
    score: i32,
    hi: Option<i32>,
}

// Read AS and HI in one traversal, retaining the first occurrence of each
// tag (including a first occurrence with a non-integer or out-of-range value).
fn alignment_tags(rec: RecordView<'_>) -> (Option<i32>, Option<i32>) {
    fn integer(ty: u8, value: &[u8]) -> Option<i32> {
        match ty {
            b'c' => Some(value[0] as i8 as i32),
            b'C' => Some(value[0] as i32),
            b's' => Some(i16::from_le_bytes(value.try_into().ok()?) as i32),
            b'S' => Some(u16::from_le_bytes(value.try_into().ok()?) as i32),
            b'i' => Some(i32::from_le_bytes(value.try_into().ok()?)),
            b'I' => i32::try_from(u32::from_le_bytes(value.try_into().ok()?)).ok(),
            _ => None,
        }
    }
    let mut values = [None; 2];
    let mut seen = [false; 2];
    let mut bytes = rec.aux().unwrap_or_default();
    while bytes.len() >= 3 {
        let ty = bytes[2];
        let body = &bytes[3..];
        let Some(size) = aux_value_len(ty, body) else {
            break;
        };
        let Some(value) = body.get(..size) else { break };
        let index = match &bytes[..2] {
            b"AS" => Some(0),
            b"HI" => Some(1),
            _ => None,
        };
        if let Some(i) = index {
            if !seen[i] {
                values[i] = integer(ty, value);
                seen[i] = true;
            }
            if seen == [true; 2] {
                break;
            }
        }
        bytes = &body[size..];
    }
    (values[0], values[1])
}

fn placements(records: &[Alignment], consumed: &mut Vec<bool>, hits: &mut Vec<Hit>) {
    consumed.clear();
    consumed.resize(records.len(), false);
    hits.clear();
    for (i, a) in records.iter().enumerate() {
        if consumed[i] || a.flags & 64 == 0 || a.mate_tid != a.tid as i32 {
            continue;
        }
        let mut partner = None;
        for (j, b) in records.iter().enumerate() {
            if consumed[j]
                || b.flags & 64 != 0
                || b.tid as i32 != a.mate_tid
                || b.pos != a.mate_pos
                || b.mate_tid != a.tid as i32
                || b.mate_pos != a.pos
            {
                continue;
            }
            if a.hi.is_none() || b.hi.is_none() || a.hi == b.hi {
                partner = Some(j);
                break;
            }
            partner.get_or_insert(j);
        }
        if let Some(j) = partner {
            consumed[i] = true;
            consumed[j] = true;
            let b = &records[j];
            hits.push(Hit {
                tid: a.tid,
                pos: a.pos.min(b.pos).max(0) as usize,
                forward: a.flags & 16 == 0,
                mate_forward: b.flags & 16 == 0,
                length: a.length.max(b.length),
                score: a.score + b.score,
                status: Status::Pair,
            });
        }
    }
    for (i, a) in records.iter().enumerate() {
        if !consumed[i] {
            hits.push(Hit {
                tid: a.tid,
                pos: a.pos.max(0) as usize,
                forward: a.flags & 16 == 0,
                mate_forward: false,
                length: a.span.min(u16::MAX as usize - 1),
                score: a.score,
                status: if a.flags & 1 == 0 {
                    Status::Single
                } else if a.flags & 64 != 0 {
                    Status::Left
                } else {
                    Status::Right
                },
            });
        }
    }
}

/// Per-fragment integer tallies taken from the original (pre-RAD) geometry:
/// fragment-length observations for both orientations and library detection.
/// All fields are integer sums, so per-batch tallies merge exactly.
pub(crate) struct Tally {
    pub num_processed: u64,
    pub opposite: Vec<u64>,
    pub same: Vec<u64>,
    pub detection: [u64; 12],
}

impl Default for Tally {
    fn default() -> Self {
        Self {
            num_processed: 0,
            opposite: vec![0; crate::fld::MAX + 1],
            same: vec![0; crate::fld::MAX + 1],
            detection: [0; 12],
        }
    }
}

impl Tally {
    fn add(&mut self, hits: &[Hit]) {
        self.num_processed += 1;
        if hits.len() == 1 && hits[0].status == Status::Pair && hits[0].length > 0 {
            let h = &hits[0];
            let bins = if h.forward != h.mate_forward {
                &mut self.opposite
            } else {
                &mut self.same
            };
            bins[h.length.min(crate::fld::MAX)] += 1;
            self.detection[h.format()] += 1;
        }
    }

    fn merge(&mut self, other: &Tally) {
        self.num_processed += other.num_processed;
        for (a, b) in self.opposite.iter_mut().zip(&other.opposite) {
            *a += b;
        }
        for (a, b) in self.same.iter_mut().zip(&other.same) {
            *a += b;
        }
        for (a, b) in self.detection.iter_mut().zip(&other.detection) {
            *a += b;
        }
    }
}

/// Raw records of whole read-name groups, copied out of the decoded stream.
/// `error` is the splitter's own (framing or record) error immediately after
/// the last record, reported only if no earlier record fails in a worker.
#[derive(Default)]
struct RawBatch {
    data: Vec<u8>,
    ends: Vec<usize>,
    error: Option<io::Error>,
}

impl RawBatch {
    fn len(&self) -> usize {
        self.ends.len()
    }

    fn record(&self, i: usize) -> &[u8] {
        let start = if i == 0 { 0 } else { self.ends[i - 1] };
        &self.data[start..self.ends[i]]
    }

    fn clear(&mut self) {
        self.data.clear();
        self.ends.clear();
        self.error = None;
    }
}

const BATCH_RECORDS: usize = 1 << 15;
const BATCH_BYTES: usize = 8 << 20;

fn group_name(rec: RecordView<'_>) -> &[u8] {
    let cname = rec.read_name();
    if cname.len() > 2 && cname[cname.len() - 2] == b'/' {
        &cname[..cname.len() - 2]
    } else {
        cname
    }
}

struct Encoded {
    words: Vec<u64>,
    tally: Tally,
}

/// Exactly the serial per-record consumer, over a batch that starts and ends
/// on read-name group boundaries.
fn process(batch: &mut RawBatch, refs: usize) -> io::Result<Encoded> {
    let mut out = Encoded {
        words: Vec::new(),
        tally: Tally::default(),
    };
    let mut name: &[u8] = &[];
    let mut records = Vec::new();
    let mut consumed = Vec::new();
    let mut hits = Vec::new();
    let mut have_group = false;
    for i in 0..batch.len() {
        let rec = RecordView::new(batch.record(i))?;
        let cname = group_name(rec);
        if have_group && name != cname {
            placements(&records, &mut consumed, &mut hits);
            out.tally.add(&hits);
            encode(&mut hits, &mut out.words);
            records.clear();
        }
        name = cname;
        have_group = true;
        if rec.flag() & 4 == 0 && rec.ref_id() >= 0 {
            if rec.ref_id() as usize >= refs {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "BAM reference id out of range",
                ));
            }
            let (score, hi) = alignment_tags(rec);
            records.push(Alignment {
                tid: rec.ref_id() as usize,
                pos: rec.pos().max(0),
                mate_tid: if rec.flag() & 8 == 0 {
                    rec.mate_ref_id()
                } else {
                    -1
                },
                mate_pos: rec.mate_pos(),
                flags: rec.flag(),
                length: rec.tlen().unsigned_abs() as usize,
                span: rec.reference_length() as usize,
                score: score.unwrap_or(0),
                hi,
            });
        }
    }
    if let Some(e) = batch.error.take() {
        return Err(e);
    }
    if have_group {
        placements(&records, &mut consumed, &mut hits);
        out.tally.add(&hits);
        encode(&mut hits, &mut out.words);
    }
    Ok(out)
}

#[derive(Default)]
struct SplitTiming {
    wait: Duration,
    split: Duration,
    blocked: Duration,
    records: u64,
}

fn lap(slot: &mut Duration, mark: &mut Option<Instant>) {
    if let Some(m) = mark {
        let now = Instant::now();
        *slot += now - *m;
        *m = now;
    }
}

/// Serial: cut the record stream into batches at read-name changes. Only
/// framing, the read name, and a copy happen here.
fn split(
    mut input: SliceBamReader<ParallelBgzfReader>,
    work: SyncSender<(usize, RawBatch)>,
    tokens: Receiver<()>,
    free: Receiver<RawBatch>,
    abort: &AtomicBool,
    timing: bool,
    batch_records: usize,
) -> SplitTiming {
    let mut t = SplitTiming::default();
    let mut mark = timing.then(Instant::now);
    let mut seq = 0;
    let mut batch = RawBatch::default();
    let mut name = Vec::new();
    let mut have_group = false;
    loop {
        if abort.load(Ordering::Relaxed) {
            return t;
        }
        let raw = input.next_record();
        lap(&mut t.wait, &mut mark);
        let raw = match raw {
            Ok(Some(raw)) => raw,
            Ok(None) => break,
            Err(e) => {
                batch.error = Some(e);
                break;
            }
        };
        let rec = match RecordView::new(raw) {
            Ok(rec) => rec,
            Err(e) => {
                batch.error = Some(e);
                break;
            }
        };
        let cname = group_name(rec);
        let new_group = !have_group || name != cname;
        if new_group
            && have_group
            && (batch.len() >= batch_records || batch.data.len() >= BATCH_BYTES)
        {
            let next = free.try_recv().unwrap_or_default();
            let full = std::mem::replace(&mut batch, next);
            lap(&mut t.split, &mut mark);
            if tokens.recv().is_err() || work.send((seq, full)).is_err() {
                return t;
            }
            lap(&mut t.blocked, &mut mark);
            seq += 1;
        }
        if new_group {
            name.clear();
            name.extend_from_slice(cname);
        }
        have_group = true;
        batch.data.extend_from_slice(raw);
        batch.ends.push(batch.data.len());
        t.records += 1;
        lap(&mut t.split, &mut mark);
    }
    if (batch.len() > 0 || batch.error.is_some()) && tokens.recv().is_ok() {
        let _ = work.send((seq, batch));
    }
    lap(&mut t.blocked, &mut mark);
    t
}

/// One BAM pass: a serial splitter cuts read-name groups into batches,
/// workers run tags/placements/tally/RAD encoding in parallel, and this
/// thread appends the encoded batches to `Fragments` in input order. The
/// result (store contents, tallies, first error) equals the serial scan.
pub(crate) fn scan(path: &Path) -> io::Result<(Header, Fragments, Tally)> {
    scan_batched(path, BATCH_RECORDS)
}

fn scan_batched(path: &Path, batch_records: usize) -> io::Result<(Header, Fragments, Tally)> {
    let threads = std::thread::available_parallelism().map_or(1, usize::from);
    let input = SliceBamReader::open(path, threads)?;
    let header = input.header.clone();
    if header.refs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BAM has no references",
        ));
    }
    let header_text = String::from_utf8_lossy(&header.text);
    if header_text.lines().any(|s| {
        s.starts_with("@HD")
            && s.split('\t').any(|f| f == "SO:coordinate")
            && !s.split('\t').any(|f| f == "GO:query")
    }) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "alignment quantification requires read-name-grouped BAM",
        ));
    }
    // NOMICS_QUANT_INFERENCE_TIMING: record wait (BGZF decode handoff), the
    // serial splitter, splitter time blocked on workers, summed worker time
    // and in-order append time. Record wait and splitter time are per-record
    // laps (two clock reads per record, as the serial scan's split was), so
    // they include that overhead; the scan phase itself grows accordingly.
    let timing = std::env::var_os("NOMICS_QUANT_INFERENCE_TIMING").is_some();
    let refs = header.refs.len();
    let inflight = 4 * threads;
    let mut fragments = Fragments::default();
    let mut tally = Tally::default();
    let mut first_error = None;
    let abort = AtomicBool::new(false);
    let work_time = Mutex::new(Duration::ZERO);
    let mut append = Duration::ZERO;
    let (work_tx, work_rx) = sync_channel::<(usize, RawBatch)>(threads);
    let work_rx = Mutex::new(work_rx);
    let split_timing = std::thread::scope(|s| {
        let (done_tx, done_rx) = channel::<(usize, io::Result<Encoded>)>();
        let (token_tx, token_rx) = sync_channel::<()>(inflight);
        for _ in 0..inflight {
            token_tx.send(()).expect("token channel has capacity");
        }
        let (free_tx, free_rx) = channel::<RawBatch>();
        let abort = &abort;
        let splitter = s.spawn(move || {
            split(
                input,
                work_tx,
                token_rx,
                free_rx,
                abort,
                timing,
                batch_records,
            )
        });
        for _ in 0..threads {
            let (done_tx, free_tx) = (done_tx.clone(), free_tx.clone());
            let (work_rx, work_time) = (&work_rx, &work_time);
            s.spawn(move || {
                let mut busy = Duration::ZERO;
                loop {
                    let job = work_rx.lock().unwrap().recv();
                    let Ok((seq, mut batch)) = job else { break };
                    let start = timing.then(Instant::now);
                    let result = process(&mut batch, refs);
                    if let Some(start) = start {
                        busy += start.elapsed();
                    }
                    batch.clear();
                    let _ = free_tx.send(batch);
                    if done_tx.send((seq, result)).is_err() {
                        break;
                    }
                }
                *work_time.lock().unwrap() += busy;
            });
        }
        drop((done_tx, free_tx));
        let mut pending = BTreeMap::new();
        let mut next = 0;
        for (seq, result) in done_rx {
            pending.insert(seq, result);
            while let Some(result) = pending.remove(&next) {
                next += 1;
                let start = timing.then(Instant::now);
                if first_error.is_none() {
                    let appended = result.and_then(|e| {
                        tally.merge(&e.tally);
                        fragments.extend(&e.words)
                    });
                    if let Err(e) = appended {
                        first_error = Some(e);
                        abort.store(true, Ordering::Relaxed);
                    }
                }
                if let Some(start) = start {
                    append += start.elapsed();
                }
                let _ = token_tx.send(());
            }
        }
        splitter.join().expect("BAM splitter panicked")
    });
    if timing {
        eprintln!(
            "salmon-quant timing: scan record wait {:.3}s, serial splitter {:.3}s, \
             splitter blocked {:.3}s, workers {:.3}s (summed, {threads} threads), \
             in-order append {:.3}s; {} records (wait and splitter include \
             2 clock reads per record)",
            split_timing.wait.as_secs_f64(),
            split_timing.split.as_secs_f64(),
            split_timing.blocked.as_secs_f64(),
            work_time.into_inner().unwrap().as_secs_f64(),
            append.as_secs_f64(),
            split_timing.records
        );
    }
    match first_error {
        Some(e) => Err(e),
        None => Ok((header, fragments, tally)),
    }
}

/// Reproduce the RAD writer/reader round trip without materializing a RAD file.
/// The fragment-level status replaces each hit's own status (even mixed hits).
pub(crate) fn rad_geometry(hits: &mut [Hit]) {
    let mut status = Status::Single;
    for hit in hits.iter() {
        // Left and right orphans have equal rank; the last one wins.
        if hit.status == Status::Pair || (status != Status::Pair && hit.status != Status::Single) {
            status = hit.status;
        }
    }
    for hit in hits {
        hit.length = hit.length.min(u16::MAX as usize);
        hit.status = status;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn placements(records: &[Alignment]) -> Vec<Hit> {
        let mut hits = Vec::new();
        super::placements(records, &mut Vec::new(), &mut hits);
        hits
    }

    #[test]
    fn compact_fragments_preserve_rad_fields_and_order_across_mappings() {
        let mut store = Fragments::default();
        store.push(&mut []).unwrap();
        let make_hit = |i: usize| Hit {
            tid: u32::MAX as usize - i,
            pos: i as u32 as usize,
            score: if i.is_multiple_of(2) {
                i32::MIN
            } else {
                i32::MAX
            },
            length: usize::MAX,
            forward: i & 1 != 0,
            mate_forward: i & 2 != 0,
            status: match i % 4 {
                0 => Status::Single,
                1 => Status::Left,
                2 => Status::Right,
                _ => Status::Pair,
            },
        };
        let count = STORE_WORDS / 3 + 3;
        for i in 0..count {
            store.push(&mut [make_hit(i)]).unwrap();
        }
        // Mixed per-hit statuses become one fragment-level status, after FLD
        // has used the original geometry; the final orphan wins equal ranks.
        let mut mixed = [make_hit(0), make_hit(1), make_hit(2)];
        store.push(&mut mixed).unwrap();
        assert!(mixed.iter().all(|h| h.status == Status::Right));
        let mut visited = 0;
        store.visit(|hits| {
            match visited {
                0 => assert!(hits.is_empty()),
                n if n <= count => {
                    let mut expected = make_hit(n - 1);
                    expected.length = u16::MAX as usize;
                    assert_eq!(hits, &[expected]);
                }
                _ => assert_eq!(hits, &mixed),
            }
            visited += 1;
        });
        assert_eq!(visited, count + 2);
        assert_eq!(store.chunks.len(), 2);
    }

    #[test]
    fn as_and_hi_keep_integer_types_first_occurrences_and_malformed_tag_behavior() {
        let tags = |aux: &[u8]| {
            let mut raw = vec![0; 32];
            raw.extend_from_slice(aux);
            alignment_tags(RecordView::new(&raw).unwrap())
        };
        for (ty, bytes, expected) in [
            (b'c', vec![255], -1),
            (b'C', vec![255], 255),
            (b's', (-123i16).to_le_bytes().to_vec(), -123),
            (b'S', u16::MAX.to_le_bytes().to_vec(), 65535),
            (b'i', i32::MIN.to_le_bytes().to_vec(), i32::MIN),
            (b'I', (123u32).to_le_bytes().to_vec(), 123),
        ] {
            let mut aux = b"ZZZskip\0AS".to_vec();
            aux.push(ty);
            aux.extend(bytes);
            aux.extend_from_slice(b"HIc\x07ASc\x09HIc\x08");
            assert_eq!(tags(&aux), (Some(expected), Some(7)));
        }
        assert_eq!(tags(b"HIZx\0HIc\x07ASc\x05"), (Some(5), None));
        assert_eq!(tags(b"ASI\xff\xff\xff\xffHIc\x07"), (None, Some(7)));
        assert_eq!(tags(b"ASc\x05ZZi\x00HIc\x07"), (Some(5), None));
        assert_eq!(tags(b"ZZqASc\x05HIc\x07"), (None, None));
        assert_eq!(tags(b"HIc\x07ASi\x00"), (None, Some(7)));
        assert_eq!(tags(&[]), (None, None));
    }

    #[test]
    fn batch_boundaries_do_not_change_the_store_or_tallies() {
        let Some(home) = std::env::var_os("HOME") else {
            return;
        };
        let bam = std::path::PathBuf::from(home).join(
            "omics-fixtures/rnaseq-test/results/star_salmon/WT_REP1.Aligned.toTranscriptome.out.bam",
        );
        if !bam.is_file() {
            eprintln!("SKIP: optional rnaseq-test fixture absent");
            return;
        }
        let snapshot = |batch_records| {
            let (_, store, tally) = scan_batched(&bam, batch_records).unwrap();
            let mut fragments = Vec::new();
            store.visit(|hits| fragments.push(hits.to_vec()));
            let mut split = Vec::new();
            for range in store.ranges() {
                store.visit_range(range, |hits| split.push(hits.to_vec()));
            }
            assert_eq!(split, fragments);
            let t = (
                tally.num_processed,
                tally.opposite,
                tally.same,
                tally.detection,
            );
            (fragments, t, store.marks.len())
        };
        let (whole, tally, marks) = snapshot(usize::MAX);
        assert_eq!(whole.len() as u64, tally.0);
        // Byte-capped batches only: `marks` already exceeds one on this BAM.
        assert!(marks > 1);
        for batch_records in [1, 2, 7, 1000] {
            let (fragments, t, marks) = snapshot(batch_records);
            // One read-name group per batch when every group change cuts.
            if batch_records == 1 {
                assert_eq!(marks, fragments.len());
            }
            assert_eq!(fragments, whole);
            assert_eq!(t, tally);
        }
    }

    fn alignment(read1: bool, pos: i32, mate_pos: i32, hi: i32) -> Alignment {
        Alignment {
            tid: 0,
            pos,
            mate_tid: 0,
            mate_pos,
            flags: 1 | if read1 { 64 | 16 } else { 128 },
            length: 200,
            span: 75,
            score: 50,
            hi: Some(hi),
        }
    }

    #[test]
    fn reciprocal_pairs_do_not_form_a_cross_product() {
        let records = [
            alignment(true, 100, 200, 1),
            alignment(true, 120, 220, 2),
            alignment(false, 200, 100, 1),
            alignment(false, 220, 120, 2),
        ];
        let hits = placements(&records);
        assert_eq!(hits.len(), 2);
        assert!(
            hits.iter()
                .all(|h| h.status == Status::Pair && crate::LibraryType::Isr.compatible(h))
        );
        assert_eq!(hits[0].score, 100);
        assert_eq!(hits[1].pos, 120);
        // An unmatched alignment must survive as an orphan, never as an invented pair.
        let hits = placements(&records[..3]);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[1].status, Status::Left);
        assert_eq!(hits[1].length, 75);
    }

    #[test]
    fn rad_pair_status_overrides_orphans_without_restoring_lost_geometry() {
        let records = [
            alignment(true, 200, 100, 1),
            alignment(false, 100, 200, 1),
            alignment(true, 300, 400, 2),
        ];
        let mut hits = placements(&records);
        // This pair's coordinate orientation is outward; RAD retains strands only.
        assert_eq!(hits[0].format(), 2);
        assert!(crate::LibraryType::Isr.compatible(&hits[0]));
        assert_eq!(hits[1].status, Status::Left);
        rad_geometry(&mut hits);
        assert!(hits.iter().all(|h| h.status == Status::Pair));
        assert_eq!(hits[1].length, 75);
        assert!(!crate::LibraryType::Isr.compatible(&hits[1]));
    }
}
