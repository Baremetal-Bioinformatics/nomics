//! Streaming BAM input: borrow slices of parallel-decoded buffers with
//! [`SliceBamReader`], or receive owned record batches with [`BamStream`].

use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread::JoinHandle;

use crate::bam::{Header, MAX_ALLOC, read_raw_record};
use crate::bgzf::ParallelBgzfReader;

/// A BAM cursor borrowing complete records directly from decoded BGZF buffers.
/// Only records (including length prefixes) crossing buffers use the reusable
/// seam buffer. The header is read once through the same decoder.
pub struct SliceBamReader<R> {
    pub header: Header,
    input: R,
    consumed: usize,
    seam: Vec<u8>,
}

impl<R: BufRead> SliceBamReader<R> {
    pub fn new(mut input: R) -> io::Result<Self> {
        let header = Header::read(&mut input)?;
        Ok(Self {
            header,
            input,
            consumed: 0,
            seam: Vec::new(),
        })
    }

    pub fn next_record(&mut self) -> io::Result<Option<&[u8]>> {
        self.input.consume(self.consumed);
        self.consumed = 0;
        let bytes = self.input.fill_buf()?;
        if bytes.is_empty() {
            return Ok(None);
        }
        if bytes.len() >= 4 {
            let size = i32::from_le_bytes(bytes[..4].try_into().unwrap());
            let size = usize::try_from(size)
                .ok()
                .filter(|&n| (32..=MAX_ALLOC).contains(&n))
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, format!("bad block_size {size}"))
                })?;
            if bytes.len() >= size + 4 {
                self.consumed = size + 4;
                return Ok(Some(&self.input.fill_buf()?[4..size + 4]));
            }
        }
        // read_raw_record also checks partial prefixes and truncated payloads.
        read_raw_record(&mut self.input, &mut self.seam)?;
        Ok(Some(&self.seam))
    }
}

impl SliceBamReader<ParallelBgzfReader> {
    pub fn open(path: &Path, inflate_threads: usize) -> io::Result<Self> {
        let file = BufReader::with_capacity(1 << 22, std::fs::File::open(path)?);
        Self::new(ParallelBgzfReader::new(file, inflate_threads))
    }
}

/// Raw records (without `block_size`), in file order.
#[derive(Default)]
pub struct Batch {
    data: Vec<u8>,
    /// Record `i` is `data[ends[i - 1]..ends[i]]` (`ends[-1]` = 0).
    ends: Vec<usize>,
}

impl Batch {
    pub fn len(&self) -> usize {
        self.ends.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ends.is_empty()
    }

    pub fn record(&self, i: usize) -> &[u8] {
        let start = if i == 0 { 0 } else { self.ends[i - 1] };
        &self.data[start..self.ends[i]]
    }
}

/// A BAM read on a background thread, `batch_records` records at a time,
/// with up to `depth` batches decoded ahead of the consumer and blocks
/// inflated on `inflate_threads` threads.
pub struct BamStream {
    pub header: Header,
    rx: Receiver<io::Result<Batch>>,
    reader: Option<JoinHandle<()>>,
}

impl BamStream {
    pub fn open(
        path: &Path,
        batch_records: usize,
        depth: usize,
        inflate_threads: usize,
    ) -> io::Result<BamStream> {
        // A BAM, or the same stream as a light spill file (`crate::spill`).
        // The BGZF EOF marker is only recommended (SAMv1 4.1.2), so any BAM
        // opens; the pipeline checks its own streams ([`check_complete`]).
        let mut r = BufReader::with_capacity(
            1 << 20,
            crate::spill::open_decoded(path, inflate_threads)?,
        );
        let header = Header::read(&mut r)?;
        let (tx, rx) = sync_channel(depth.max(1));
        let batch_records = batch_records.max(1);
        let reader = std::thread::spawn(move || {
            let mut rec = Vec::new();
            loop {
                let mut batch = Batch::default();
                let mut result = Ok(());
                while batch.len() < batch_records {
                    match read_raw_record(&mut r, &mut rec) {
                        Ok(true) => {
                            batch.data.extend_from_slice(&rec);
                            batch.ends.push(batch.data.len());
                        }
                        Ok(false) => break,
                        Err(e) => {
                            result = Err(e);
                            break;
                        }
                    }
                }
                let done = batch.len() < batch_records || result.is_err();
                let message = result.map(|()| batch);
                if (message.as_ref().is_ok_and(|b| !b.is_empty()) || message.is_err())
                    && tx.send(message).is_err()
                {
                    return;
                }
                if done {
                    return;
                }
            }
        });
        Ok(BamStream {
            header,
            rx,
            reader: Some(reader),
        })
    }

    /// The next non-empty batch, or `None` at the end of the file.
    pub fn next_batch(&mut self) -> io::Result<Option<Batch>> {
        match self.rx.recv() {
            Ok(batch) => batch.map(Some),
            Err(_) => {
                if let Some(reader) = self.reader.take() {
                    reader
                        .join()
                        .map_err(|_| io::Error::other("BAM reader thread panicked"))?;
                }
                Ok(None)
            }
        }
    }
}

/// One record at a time over a [`BamStream`].
pub struct RecordCursor {
    stream: BamStream,
    batch: Batch,
    at: usize,
}

impl RecordCursor {
    pub fn new(stream: BamStream) -> RecordCursor {
        RecordCursor {
            stream,
            batch: Batch::default(),
            at: 0,
        }
    }

    pub fn header(&self) -> &Header {
        &self.stream.header
    }

    /// The current record, refilling from the stream as needed; `None` at the end.
    pub fn peek(&mut self) -> io::Result<Option<&[u8]>> {
        if self.at == self.batch.len() {
            match self.stream.next_batch()? {
                Some(batch) => {
                    self.batch = batch;
                    self.at = 0;
                }
                None => return Ok(None),
            }
        }
        Ok(Some(self.batch.record(self.at)))
    }

    pub fn advance(&mut self) {
        self.at += 1;
    }
}

/// Suffix of a BAM [`StreamWriter`] (or [`write_set`]) is still writing.
const PARTIAL_SUFFIX: &str = ".partial";

fn partial_path(path: &Path) -> std::path::PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(PARTIAL_SUFFIX);
    std::path::PathBuf::from(p)
}

/// Suffix of the record the pipeline keeps beside a stream it wrote (its
/// generation). A stream's writer removes it, durably, just before the new
/// stream takes the name, so the record never describes another file.
pub const RECORD_SUFFIX: &str = ".gen";

/// Remove `path`'s [`RECORD_SUFFIX`] record (if any) and sync the directory.
fn drop_record(path: &Path) -> io::Result<()> {
    let mut p = path.as_os_str().to_owned();
    p.push(RECORD_SUFFIX);
    match std::fs::remove_file(std::path::PathBuf::from(p)) {
        Ok(()) => sync_parent(path),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e),
    }
}

/// For tools reading BAMs from elsewhere: warn loudly (naming the file) when
/// `path` lacks the BGZF EOF marker, which SAMv1 4.1.2 only recommends; a
/// file cut at a block boundary would read as a shorter, valid one.
pub fn warn_if_incomplete(path: &Path) {
    if let Err(e) = check_complete(path) {
        eprintln!("WARNING: {e}; reading it anyway, but its records may be cut short");
    }
}

/// fsync the directory holding `path`, so a rename into it is durable.
fn sync_parent(path: &Path) -> io::Result<()> {
    let dir = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::File::open(dir)?.sync_all()
}

/// Whether `path` ends as a complete stream of ours does: a BAM with the
/// BGZF EOF marker (our writers always add it), or a spill file with its
/// footer magic. A BAM cut at a block (and record) boundary otherwise reads
/// back cleanly as a shorter file. For selecting and resuming the
/// pipeline's own streams; [`BamStream::open`] does not require it.
pub fn check_complete(path: &Path) -> io::Result<()> {
    use std::io::{Read, Seek, SeekFrom};
    let mut f = std::fs::File::open(path)?;
    let len = f.metadata()?.len();
    let (tail, what): (&[u8], &str) = if crate::spill::is_spill_path(path) {
        (crate::spill::FOOTER_MAGIC, "spill footer")
    } else {
        (&crate::bgzf::EOF_BLOCK, "BGZF EOF marker")
    };
    let mut got = vec![0u8; tail.len()];
    let ok = len >= tail.len() as u64 && {
        f.seek(SeekFrom::Start(len - tail.len() as u64))?;
        f.read_exact(&mut got)?;
        got == tail
    };
    if ok {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: no {what} at the end (truncated, or still being written)",
                path.display()
            ),
        ))
    }
}

/// A BAM's uncompressed stream (header, then `block_size`-prefixed records)
/// going to disk as a BAM, or as a light spill file ([`crate::spill`]) when
/// the path has the spill extension. [`BamStream`] reads either back.
///
/// Either way the stream is written under a partial name and given its own
/// name only by [`StreamWriter::finish`], once complete and synced (with
/// the directory synced after the rename); dropped unfinished, the partial
/// file is removed. A file with the final name is thus always complete.
pub enum StreamWriter {
    /// The writer every streamed BAM has always had (same block cuts).
    Bam(BamFile),
    /// Written under a partial name, synced and renamed by `finish`.
    Spill(crate::spill::SpillFile),
}

/// A BAM [`StreamWriter`] in progress.
pub struct BamFile {
    writer: Option<crate::bgzf::ParallelBgzfWriter<std::io::BufWriter<std::fs::File>>>,
    partial: std::path::PathBuf,
    path: std::path::PathBuf,
}

impl Drop for BamFile {
    fn drop(&mut self) {
        if self.writer.is_some() {
            drop(self.writer.take());
            let _ = std::fs::remove_file(&self.partial);
        }
    }
}

impl StreamWriter {
    /// Create `path` (under its partial name), compressing on `threads` threads.
    pub fn create(path: &Path, threads: usize) -> io::Result<StreamWriter> {
        Ok(if crate::spill::is_spill_path(path) {
            StreamWriter::Spill(crate::spill::SpillFile::create(path, threads)?)
        } else {
            let partial = partial_path(path);
            StreamWriter::Bam(BamFile {
                writer: Some(crate::bgzf::ParallelBgzfWriter::for_bam(
                    std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(&partial)?),
                    threads,
                )),
                partial,
                path: path.to_path_buf(),
            })
        })
    }

    /// Write the end of the stream (BGZF EOF block, or the spill footer),
    /// sync it, give it its name and sync the directory.
    pub fn finish(self) -> io::Result<()> {
        use std::io::Write;
        match self {
            StreamWriter::Bam(mut b) => {
                let w = b.writer.take().expect("BAM stream writer");
                let mut f = w.finish()?;
                f.flush()?;
                let file = f.into_inner().map_err(|e| e.into_error())?;
                let renamed = file.sync_all().and_then(|()| {
                    drop(file);
                    drop_record(&b.path)?;
                    std::fs::rename(&b.partial, &b.path)
                });
                if let Err(e) = renamed {
                    let _ = std::fs::remove_file(&b.partial);
                    return Err(e);
                }
                sync_parent(&b.path)
            }
            StreamWriter::Spill(s) => {
                let path = s.path().to_path_buf();
                drop_record(&path)?;
                s.finish()?;
                sync_parent(&path)
            }
        }
    }

    fn inner(&mut self) -> &mut dyn std::io::Write {
        match self {
            StreamWriter::Bam(b) => b.writer.as_mut().expect("BAM stream writer"),
            StreamWriter::Spill(s) => s,
        }
    }
}

impl std::io::Write for StreamWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner().write(buf)
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        self.inner().write_all(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner().flush()
    }
}

/// Write `set` to `path` as [`StreamWriter`] would (BAM or spill by
/// extension, complete under its name or not there at all). The BAM bytes
/// are [`RecordSet::write`]'s.
pub fn write_set(set: &crate::bam::RecordSet, path: &Path, threads: usize) -> io::Result<()> {
    if crate::spill::is_spill_path(path) {
        drop_record(path)?;
        set.write(path, threads)?;
        return sync_parent(path);
    }
    let partial = partial_path(path);
    let written = set.write(&partial, threads).and_then(|()| {
        std::fs::File::open(&partial)?.sync_all()?;
        drop_record(path)?;
        std::fs::rename(&partial, path)
    });
    if let Err(e) = written {
        let _ = std::fs::remove_file(&partial);
        return Err(e);
    }
    sync_parent(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::RecordSet;

    fn bam_bytes() -> (Header, Vec<Vec<u8>>, Vec<u8>) {
        let header = Header {
            text: b"@HD\tVN:1.6\tGO:query\n".to_vec(),
            refs: Vec::new(),
        };
        let records: Vec<_> = (0..40).map(|i| vec![i; 32 + i as usize]).collect();
        let mut bytes = Vec::new();
        header.write(&mut bytes).unwrap();
        for record in &records {
            crate::bam::write_raw_record(&mut bytes, record).unwrap();
        }
        (header, records, bytes)
    }

    #[test]
    fn slice_records_borrow_input_and_copy_only_at_buffer_seams() {
        let (header, records, bytes) = bam_bytes();
        let mut direct = SliceBamReader::new(bytes.as_slice()).unwrap();
        assert_eq!(direct.header, header);
        for record in &records {
            let raw = direct.next_record().unwrap().unwrap();
            assert_eq!(raw, record);
            assert!(bytes.as_ptr_range().contains(&raw.as_ptr()));
        }
        assert!(direct.seam.is_empty());
        assert!(direct.next_record().unwrap().is_none());
        assert!(direct.next_record().unwrap().is_none());

        // All prefix and payload offsets, including records spanning many
        // buffers, and the header itself crossing buffer boundaries.
        for capacity in 1..100 {
            let mut reader =
                SliceBamReader::new(BufReader::with_capacity(capacity, bytes.as_slice())).unwrap();
            assert_eq!(reader.header, header);
            for record in &records {
                assert_eq!(reader.next_record().unwrap().unwrap(), record);
            }
            assert!(reader.next_record().unwrap().is_none());
        }
    }

    #[test]
    fn slice_records_reject_truncated_and_invalid_framing() {
        let (header, _, _) = bam_bytes();
        let mut prefix = Vec::new();
        header.write(&mut prefix).unwrap();
        for tail in [
            vec![32],
            vec![32, 0, 0],
            31i32.to_le_bytes().to_vec(),
            (-1i32).to_le_bytes().to_vec(),
            (MAX_ALLOC as i32 + 1).to_le_bytes().to_vec(),
            [32i32.to_le_bytes().as_slice(), &[0; 31]].concat(),
        ] {
            let mut bytes = prefix.clone();
            bytes.extend(tail);
            for capacity in [1, 3, 7, 1000] {
                let mut reader =
                    SliceBamReader::new(BufReader::with_capacity(capacity, bytes.as_slice()))
                        .unwrap();
                assert!(reader.next_record().is_err());
            }
        }
    }

    #[test]
    fn stream_returns_every_record_in_order() {
        let header = Header {
            text: b"@HD\tVN:1.6\tSO:coordinate\n".to_vec(),
            refs: Vec::new(),
        };
        let mut records = Vec::new();
        for i in 0..1000u32 {
            let mut raw = vec![0u8; 36];
            raw[0..4].copy_from_slice(&(-1i32).to_le_bytes());
            raw[4..8].copy_from_slice(&(-1i32).to_le_bytes());
            raw[8] = 4; // l_read_name: "rNN\0" padded below
            raw[32..36].copy_from_slice(&i.to_le_bytes());
            records.extend_from_slice(&(raw.len() as i32).to_le_bytes());
            records.extend_from_slice(&raw);
        }
        let set = RecordSet::from_parts(header, &records).unwrap();
        let path =
            std::env::temp_dir().join(format!("nomics-bamstream-{}.bam", std::process::id()));
        set.write(&path, 2).unwrap();
        let mut cursor = RecordCursor::new(BamStream::open(&path, 7, 2, 3).unwrap());
        let mut n = 0;
        while let Some(rec) = cursor.peek().unwrap() {
            assert_eq!(rec, set.record(n));
            cursor.advance();
            n += 1;
        }
        assert_eq!(n, set.len());
        std::fs::remove_file(&path).unwrap();
    }

    /// A stream is under its own name only once finished: dropped
    /// unfinished it leaves nothing, and readers refuse a BAM without its
    /// BGZF EOF marker (cut at a block boundary, it would read cleanly).
    #[test]
    fn stream_writer_names_only_complete_streams() {
        let header = Header {
            text: b"@HD\tVN:1.6\n".to_vec(),
            refs: Vec::new(),
        };
        for ext in ["bam", "nspill"] {
            let path = std::env::temp_dir()
                .join(format!("nomics-streamwriter-{}.{ext}", std::process::id()));
            {
                let mut w = StreamWriter::create(&path, 2).unwrap();
                header.write(&mut w).unwrap();
            }
            assert!(!path.exists(), "{ext}: unfinished stream under its name");
            assert!(!partial_path(&path).exists(), "{ext}: partial left behind");
            let mut w = StreamWriter::create(&path, 2).unwrap();
            header.write(&mut w).unwrap();
            w.finish().unwrap();
            check_complete(&path).unwrap();
            assert!(BamStream::open(&path, 10, 1, 1).is_ok());
            let full = std::fs::read(&path).unwrap();
            std::fs::write(&path, &full[..full.len() - 8]).unwrap();
            assert!(check_complete(&path).is_err(), "{ext}: truncated accepted");
            if ext == "bam" {
                // Without only its EOF marker (optional in SAMv1 4.1.2), an
                // external BAM still opens and reads; only the pipeline's
                // own streams must have it.
                std::fs::write(&path, &full[..full.len() - 28]).unwrap();
                assert!(check_complete(&path).is_err());
                let mut s = BamStream::open(&path, 10, 1, 1).unwrap();
                assert_eq!(s.header.text, header.text);
                assert!(s.next_batch().unwrap().is_none());
            }
            std::fs::remove_file(&path).unwrap();
        }
    }

    /// A rewrite that fails (writer dropped before `finish`) leaves the
    /// finished stream of that name whole, BAM or spill: only `finish`'s
    /// rename replaces it.
    #[test]
    fn failed_rewrite_keeps_the_old_stream() {
        let header = Header {
            text: b"@HD\tVN:1.6\n".to_vec(),
            refs: Vec::new(),
        };
        for ext in ["bam", "nspill"] {
            let path = std::env::temp_dir()
                .join(format!("nomics-rewrite-{}.{ext}", std::process::id()));
            let mut w = StreamWriter::create(&path, 1).unwrap();
            header.write(&mut w).unwrap();
            w.finish().unwrap();
            let old = std::fs::read(&path).unwrap();
            {
                let mut w = StreamWriter::create(&path, 1).unwrap();
                header.write(&mut w).unwrap();
                std::io::Write::write_all(&mut w, &[0u8; 1000]).unwrap();
            }
            assert_eq!(std::fs::read(&path).unwrap(), old, "{ext}: old stream lost");
            check_complete(&path).unwrap();
            assert!(!partial_path(&path).exists());
            std::fs::remove_file(&path).unwrap();
        }
    }
}
