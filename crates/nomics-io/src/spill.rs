//! A light, internal spill format for BAM streams that are written once and
//! read back by this process only (whole-genome chunk BAMs when they are not
//! published, `--save-mapped false`).
//!
//! The decoded bytes are exactly a BAM's uncompressed stream (header, then
//! `block_size`-prefixed records), so readers above the byte layer (header
//! parsing, [`crate::bamstream::BamStream`]) are unchanged. What changes is
//! the compression: instead of DEFLATE (whose Huffman decoding is the cost
//! of reading a BAM back), blocks are compressed with a byte-aligned LZ77
//! code in the style of LZ4, which decodes at close to memory-copy speed.
//!
//! File layout (little-endian):
//! - [`MAGIC`];
//! - frames: `raw_len: u32`, `stored_len: u32` (bit 31 set: the payload is
//!   the raw bytes), `crc32(raw): u32`, then the payload;
//! - an end frame of three zero words;
//! - a footer: total decoded bytes (`u64`), frame count (`u64`),
//!   [`FOOTER_MAGIC`]; then the end of the file. A reader checks all three,
//!   so a truncated file, a zeroed frame header or trailing bytes are
//!   errors, never a short stream.
//!
//! Files are written under a `.partial` name and renamed once complete
//! ([`SpillFile`]); a directory's set of chunk spill files is described by a
//! [`MANIFEST`] written once the set is complete.
//!
//! A compressed payload is a sequence of LZ77 sequences: a token (high
//! nibble: literal count, low nibble: match length - 4; 15 means more length
//! bytes follow, each added, until one is not 255), the literals, then a
//! 2-byte match offset and any match length bytes. The last sequence has
//! literals only and ends the payload.

use std::io::{self, BufRead, Read, Write};
use std::path::{Path, PathBuf};

/// First bytes of a spill file.
pub const MAGIC: &[u8; 8] = b"NOMSPL02";
/// The magic's prefix shared by every version (older ones are refused).
const MAGIC_FAMILY: &[u8] = b"NOMSPL";
/// Last bytes of a complete spill file.
pub const FOOTER_MAGIC: &[u8; 8] = b"NOMSPEND";
const FOOTER_LEN: usize = 24;
/// Suffix of a spill file (or manifest) still being written.
const PARTIAL_SUFFIX: &str = ".partial";
/// File extension of spill files: [`is_spill_path`].
pub const EXTENSION: &str = "nspill";
/// Raw bytes per block (the last one may be shorter).
pub const BLOCK: usize = 1 << 20;

const STORED: u32 = 1 << 31;
const MIN_MATCH: usize = 4;
const MAX_OFFSET: usize = 65535;
const HASH_BITS: u32 = 16;
/// No match starts in the last bytes of a block: 8-byte reads stay in bounds.
const TAIL: usize = 8;

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// Whether `path` names a spill file (extension [`EXTENSION`]).
pub fn is_spill_path(path: &Path) -> bool {
    path.extension().is_some_and(|e| e == EXTENSION)
}

fn read_u32(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

fn read_u64(b: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(b[at..at + 8].try_into().unwrap())
}

fn push_len(mut n: usize, out: &mut Vec<u8>) {
    while n >= 255 {
        out.push(255);
        n -= 255;
    }
    out.push(n as u8);
}

/// One sequence: `literals`, then (unless `None`) a match.
fn emit(literals: &[u8], m: Option<(usize, usize)>, out: &mut Vec<u8>) {
    let lit = literals.len();
    let ml = m.map_or(0, |(_, len)| len - MIN_MATCH);
    out.push(((lit.min(15) as u8) << 4) | ml.min(15) as u8);
    if lit >= 15 {
        push_len(lit - 15, out);
    }
    out.extend_from_slice(literals);
    if let Some((offset, _)) = m {
        out.extend_from_slice(&(offset as u16).to_le_bytes());
        if ml >= 15 {
            push_len(ml - 15, out);
        }
    }
}

/// Compress `src` as one payload, appended to `out`. `table` is scratch
/// space reused between calls.
pub fn compress(src: &[u8], table: &mut Vec<u32>, out: &mut Vec<u8>) {
    table.clear();
    table.resize(1 << HASH_BITS, u32::MAX);
    let n = src.len();
    let limit = n.saturating_sub(TAIL);
    let (mut anchor, mut ip, mut misses) = (0usize, 0usize, 0usize);
    let hash = |v: u32| (v.wrapping_mul(2_654_435_761) >> (32 - HASH_BITS)) as usize;
    while ip < limit {
        let v = read_u32(src, ip);
        let h = hash(v);
        let cand = table[h] as usize;
        table[h] = ip as u32;
        if cand == u32::MAX as usize || ip - cand > MAX_OFFSET || read_u32(src, cand) != v {
            // Skip faster through data that does not compress.
            misses += 1;
            ip += 1 + (misses >> 6);
            continue;
        }
        misses = 0;
        let (mut start, mut from) = (ip, cand);
        while start > anchor && from > 0 && src[start - 1] == src[from - 1] {
            start -= 1;
            from -= 1;
        }
        let mut len = ip + MIN_MATCH - start;
        loop {
            let (a, b) = (from + len, start + len);
            if b + 8 > n {
                while start + len < n && src[from + len] == src[start + len] {
                    len += 1;
                }
                break;
            }
            let x = read_u64(src, a) ^ read_u64(src, b);
            if x != 0 {
                len += (x.trailing_zeros() / 8) as usize;
                break;
            }
            len += 8;
        }
        emit(&src[anchor..start], Some((start - from, len)), out);
        ip = start + len;
        anchor = ip;
        if ip >= 2 && ip - 2 < limit {
            table[hash(read_u32(src, ip - 2))] = (ip - 2) as u32;
        }
    }
    emit(&src[anchor..], None, out);
}

fn read_len(src: &[u8], at: &mut usize, mut n: usize) -> io::Result<usize> {
    loop {
        let b = *src
            .get(*at)
            .ok_or_else(|| invalid("spill block truncated"))?;
        *at += 1;
        n += b as usize;
        if b != 255 {
            return Ok(n);
        }
    }
}

/// Slack past the end of a decoded block for 16-byte copies.
const WILD: usize = 32;

/// Decompress one payload of `raw_len` bytes, appended to `out`.
///
/// Short literal runs and matches are copied 16 bytes at a time into the
/// slack past what has been decoded (overwritten by the next sequence and
/// cut at the end), which is what makes this faster than inflate.
pub fn decompress(src: &[u8], raw_len: usize, out: &mut Vec<u8>) -> io::Result<()> {
    let base = out.len();
    out.resize(base + raw_len + WILD, 0);
    let result = decode_into(src, &mut out[base..], raw_len);
    out.truncate(base + raw_len);
    result
}

fn decode_into(src: &[u8], dst: &mut [u8], raw_len: usize) -> io::Result<()> {
    let truncated = || invalid("spill block truncated");
    let corrupt = || invalid("spill block corrupt");
    let (mut s, mut d) = (0usize, 0usize);
    loop {
        let token = *src.get(s).ok_or_else(truncated)?;
        s += 1;
        let mut lit = (token >> 4) as usize;
        if lit == 15 {
            lit = read_len(src, &mut s, lit)?;
        }
        if d + lit > raw_len {
            return Err(corrupt());
        }
        if lit <= 16 && s + 16 <= src.len() {
            dst[d..d + 16].copy_from_slice(&src[s..s + 16]);
        } else {
            dst[d..d + lit].copy_from_slice(src.get(s..s + lit).ok_or_else(truncated)?);
        }
        s += lit;
        d += lit;
        if s >= src.len() {
            if s > src.len() {
                return Err(truncated());
            }
            break;
        }
        let offset = src
            .get(s..s + 2)
            .map(|b| u16::from_le_bytes([b[0], b[1]]) as usize)
            .ok_or_else(truncated)?;
        s += 2;
        let mut len = (token & 15) as usize;
        if len == 15 {
            len = read_len(src, &mut s, len)?;
        }
        len += MIN_MATCH;
        if offset == 0 || offset > d || d + len > raw_len {
            return Err(corrupt());
        }
        let from = d - offset;
        if offset >= 16 {
            // Each 16-byte step reads only bytes already final.
            let mut k = 0;
            while k < len {
                dst.copy_within(from + k..from + k + 16, d + k);
                k += 16;
            }
        } else {
            // A repeat of the last `offset` bytes: copy runs that only read
            // bytes already written, doubling each time.
            let mut done = 0;
            while done < len {
                let n = (len - done).min(d + done - from);
                dst.copy_within(from..from + n, d + done);
                done += n;
            }
        }
        d += len;
    }
    if d != raw_len {
        return Err(invalid("spill block has the wrong length"));
    }
    Ok(())
}

/// Encode one block as a frame (header and payload), appended to `out`.
fn encode_frame(raw: &[u8], table: &mut Vec<u32>, out: &mut Vec<u8>) {
    let at = out.len();
    out.extend_from_slice(&[0; 12]);
    compress(raw, table, out);
    let mut stored = (out.len() - at - 12) as u32;
    if stored as usize >= raw.len() {
        out.truncate(at + 12);
        out.extend_from_slice(raw);
        stored = raw.len() as u32 | STORED;
    }
    out[at..at + 4].copy_from_slice(&(raw.len() as u32).to_le_bytes());
    out[at + 4..at + 8].copy_from_slice(&stored.to_le_bytes());
    out[at + 8..at + 12].copy_from_slice(&crate::crc32::crc32(raw).to_le_bytes());
}

/// A message to the compressing thread.
enum Msg {
    Blocks(Vec<Vec<u8>>),
    /// Every block has been sent: write the end frame and the footer.
    Finish,
}

/// Writes a spill stream: blocks of [`BLOCK`] bytes compressed on `threads`
/// threads by a background thread, while the caller fills the next batch.
/// The end frame and footer are written only by [`SpillWriter::finish`]; a
/// writer dropped without it (an error, a panic) leaves a stream readers
/// reject. To a file, use [`SpillFile`], which also names the file only
/// once it is complete.
pub struct SpillWriter<W: Write + Send + 'static> {
    threads: usize,
    buf: Vec<u8>,
    batch: Vec<Vec<u8>>,
    tx: Option<std::sync::mpsc::SyncSender<Msg>>,
    recycled: std::sync::mpsc::Receiver<Vec<u8>>,
    worker: Option<std::thread::JoinHandle<io::Result<W>>>,
}

impl<W: Write + Send + 'static> SpillWriter<W> {
    pub fn new(mut inner: W, threads: usize) -> Self {
        let threads = threads.max(1);
        let (tx, rx) = std::sync::mpsc::sync_channel::<Msg>(1);
        let (recycle, recycled) = std::sync::mpsc::channel::<Vec<u8>>();
        let worker = std::thread::spawn(move || -> io::Result<W> {
            inner.write_all(MAGIC)?;
            let mut tables: Vec<Vec<u32>> = vec![Vec::new(); threads];
            let mut frames: Vec<Vec<u8>> = vec![Vec::new(); threads];
            let (mut n_frames, mut raw_total) = (0u64, 0u64);
            loop {
                let blocks = match rx.recv() {
                    Ok(Msg::Blocks(blocks)) => blocks,
                    Ok(Msg::Finish) => break,
                    Err(_) => {
                        return Err(io::Error::other(
                            "spill writer dropped before finish: stream left incomplete",
                        ));
                    }
                };
                let per = blocks.len().div_ceil(threads).max(1);
                std::thread::scope(|scope| {
                    for ((frame, table), chunk) in frames
                        .iter_mut()
                        .zip(tables.iter_mut())
                        .zip(blocks.chunks(per))
                    {
                        scope.spawn(move || {
                            frame.clear();
                            for raw in chunk {
                                encode_frame(raw, table, frame);
                            }
                        });
                    }
                });
                for frame in &frames[..blocks.len().div_ceil(per)] {
                    inner.write_all(frame)?;
                }
                n_frames += blocks.len() as u64;
                for b in blocks {
                    raw_total += b.len() as u64;
                    let _ = recycle.send(b);
                }
            }
            inner.write_all(&[0; 12])?;
            inner.write_all(&footer(raw_total, n_frames))?;
            inner.flush()?;
            Ok(inner)
        });
        SpillWriter {
            threads,
            buf: Vec::with_capacity(BLOCK),
            batch: Vec::new(),
            tx: Some(tx),
            recycled,
            worker: Some(worker),
        }
    }

    fn send(&mut self, msg: Msg) -> io::Result<()> {
        if self.tx.as_ref().is_some_and(|tx| tx.send(msg).is_ok()) {
            return Ok(());
        }
        // The worker stopped on a write error; surface it.
        Err(self
            .join()
            .err()
            .unwrap_or_else(|| io::Error::other("spill compression thread stopped")))
    }

    fn send_batch(&mut self) -> io::Result<()> {
        if self.batch.is_empty() {
            return Ok(());
        }
        let blocks = std::mem::take(&mut self.batch);
        self.send(Msg::Blocks(blocks))
    }

    fn join(&mut self) -> io::Result<W> {
        drop(self.tx.take());
        match self.worker.take() {
            Some(w) => w
                .join()
                .map_err(|_| io::Error::other("spill compression thread panicked"))?,
            None => Err(io::Error::other("spill compression thread already joined")),
        }
    }

    /// Write the last block, the end frame and the footer; returns the
    /// inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        if !self.buf.is_empty() {
            let last = std::mem::take(&mut self.buf);
            self.batch.push(last);
        }
        self.send_batch()?;
        self.send(Msg::Finish)?;
        self.join()
    }
}

impl<W: Write + Send + 'static> Write for SpillWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let n = data.len().min(BLOCK - self.buf.len());
        self.buf.extend_from_slice(&data[..n]);
        if self.buf.len() == BLOCK {
            let mut next = self.recycled.try_recv().unwrap_or_default();
            next.clear();
            next.reserve(BLOCK);
            self.batch.push(std::mem::replace(&mut self.buf, next));
            if self.batch.len() == self.threads * 2 {
                self.send_batch()?;
            }
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The footer after the end frame: decoded length, frame count, marker.
fn footer(raw_total: u64, frames: u64) -> [u8; FOOTER_LEN] {
    let mut f = [0u8; FOOTER_LEN];
    f[..8].copy_from_slice(&raw_total.to_le_bytes());
    f[8..16].copy_from_slice(&frames.to_le_bytes());
    f[16..].copy_from_slice(FOOTER_MAGIC);
    f
}

/// A spill file being written: the stream goes to `<path>.partial`, renamed
/// to `path` by [`SpillFile::finish`] once the footer is on disk. Dropped
/// unfinished (an error, a panic), it removes the partial file, so a file
/// named `path` always holds a complete stream.
pub struct SpillFile {
    writer: Option<SpillWriter<io::BufWriter<std::fs::File>>>,
    partial: PathBuf,
    path: PathBuf,
    done: bool,
}

impl SpillFile {
    pub fn create(path: &Path, threads: usize) -> io::Result<SpillFile> {
        let mut partial = path.as_os_str().to_owned();
        partial.push(PARTIAL_SUFFIX);
        let partial = PathBuf::from(partial);
        // A complete file of this name from an earlier run stays until
        // `finish` renames the new one over it: a rewrite that fails leaves
        // the old stream whole (callers drop what depends on it, such as
        // the pipeline's `.gen` and CRAM markers, themselves).
        let f = io::BufWriter::with_capacity(1 << 22, std::fs::File::create(&partial)?);
        Ok(SpillFile {
            writer: Some(SpillWriter::new(f, threads)),
            partial,
            path: path.to_path_buf(),
            done: false,
        })
    }

    /// Finish the stream, flush it to disk and give it its name.
    pub fn finish(mut self) -> io::Result<()> {
        let writer = self.writer.take().expect("spill file writer");
        let file = writer.finish()?.into_inner().map_err(|e| e.into_error())?;
        file.sync_all()?;
        std::fs::rename(&self.partial, &self.path)?;
        self.done = true;
        Ok(())
    }

    /// The name the file gets once finished.
    pub fn path(&self) -> &Path {
        &self.path
    }

    fn writer(&mut self) -> &mut SpillWriter<io::BufWriter<std::fs::File>> {
        self.writer.as_mut().expect("spill file writer")
    }
}

impl Write for SpillFile {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.writer().write(data)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Drop for SpillFile {
    fn drop(&mut self) {
        if self.done {
            return;
        }
        // Stop the compressing thread (it leaves the stream incomplete),
        // then remove what it wrote.
        if let Some(mut w) = self.writer.take() {
            let _ = w.join();
        }
        let _ = std::fs::remove_file(&self.partial);
    }
}

/// A frame's header, then its payload, read from `r`; `None` at the end frame.
fn read_frame<R: Read>(r: &mut R, payload: &mut Vec<u8>) -> io::Result<Option<(usize, u32, u32)>> {
    let mut head = [0u8; 12];
    r.read_exact(&mut head).map_err(|e| match e.kind() {
        io::ErrorKind::UnexpectedEof => invalid("spill file truncated (no end frame)"),
        _ => e,
    })?;
    let (raw_len, stored, crc) = (read_u32(&head, 0), read_u32(&head, 4), read_u32(&head, 8));
    if raw_len == 0 && stored == 0 && crc == 0 {
        return Ok(None);
    }
    let len = (stored & !STORED) as usize;
    if raw_len as usize > BLOCK || len > BLOCK + BLOCK / 2 {
        return Err(invalid("spill frame too large"));
    }
    payload.clear();
    payload.resize(len, 0);
    r.read_exact(payload)
        .map_err(|_| invalid("spill file truncated"))?;
    Ok(Some((raw_len as usize, stored, crc)))
}

/// After the end frame: the footer must match what was read, and nothing
/// may follow it.
fn check_end<R: Read>(r: &mut R, raw_total: u64, frames: u64) -> io::Result<()> {
    let mut f = [0u8; FOOTER_LEN];
    r.read_exact(&mut f)
        .map_err(|_| invalid("spill file truncated (no footer)"))?;
    if f != footer(raw_total, frames) {
        return Err(invalid(format!(
            "spill footer does not match the stream ({frames} frames, {raw_total} bytes read)"
        )));
    }
    let mut extra = [0u8; 1];
    if read_up_to(r, &mut extra)? != 0 {
        return Err(invalid("bytes after the spill footer"));
    }
    Ok(())
}

fn decode_frame(
    payload: &[u8],
    raw_len: usize,
    stored: u32,
    crc: u32,
    out: &mut Vec<u8>,
) -> io::Result<()> {
    out.clear();
    if stored & STORED != 0 {
        if payload.len() != raw_len {
            return Err(invalid("stored spill frame has the wrong length"));
        }
        out.extend_from_slice(payload);
    } else {
        decompress(payload, raw_len, out)?;
    }
    if crate::crc32::crc32(out) != crc {
        return Err(invalid("spill block CRC mismatch"));
    }
    Ok(())
}

/// Reads a spill file's decoded stream: a background thread reads frames
/// and decodes them on `threads` threads, ahead of the consumer. The end
/// frame, footer and physical end of file are checked before the last
/// frames are handed over, so a damaged file is an error, never a short
/// stream.
pub struct SpillReader {
    rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    recycle: std::sync::mpsc::Sender<Vec<u8>>,
    out: Vec<u8>,
    pos: usize,
    eof: bool,
    worker: Option<std::thread::JoinHandle<()>>,
}

impl SpillReader {
    pub fn new<R: Read + Send + 'static>(mut inner: R, threads: usize) -> Self {
        let threads = threads.max(1);
        let (tx, rx) = std::sync::mpsc::sync_channel::<io::Result<Vec<u8>>>(threads * 2);
        let (recycle, recycled) = std::sync::mpsc::channel::<Vec<u8>>();
        let worker = std::thread::spawn(move || {
            let mut magic = [0u8; 8];
            if let Err(e) = inner.read_exact(&mut magic) {
                let _ = tx.send(Err(e));
                return;
            }
            if &magic != MAGIC {
                let _ = tx.send(Err(invalid("not a spill file (or an older format)")));
                return;
            }
            let mut payloads: Vec<Vec<u8>> = vec![Vec::new(); threads * 2];
            let mut frames: Vec<(usize, u32, u32)> = Vec::new();
            let (mut n_frames, mut raw_total) = (0u64, 0u64);
            loop {
                frames.clear();
                let mut end = false;
                while frames.len() < payloads.len() {
                    let read = read_frame(&mut inner, &mut payloads[frames.len()]).and_then(|f| {
                        if f.is_none() {
                            check_end(&mut inner, raw_total, n_frames)?;
                        }
                        Ok(f)
                    });
                    match read {
                        Ok(Some(f)) => {
                            n_frames += 1;
                            raw_total += f.0 as u64;
                            frames.push(f);
                        }
                        Ok(None) => {
                            end = true;
                            break;
                        }
                        Err(e) => {
                            let _ = tx.send(Err(e));
                            return;
                        }
                    }
                }
                let mut outs: Vec<Vec<u8>> = (0..frames.len())
                    .map(|_| recycled.try_recv().unwrap_or_default())
                    .collect();
                let per = frames.len().div_ceil(threads).max(1);
                let decoded: Vec<io::Result<()>> = std::thread::scope(|scope| {
                    let handles: Vec<_> = outs
                        .chunks_mut(per)
                        .zip(frames.chunks(per))
                        .zip(payloads.chunks(per))
                        .map(|((outs, frames), payloads)| {
                            scope.spawn(move || -> io::Result<()> {
                                for ((out, &(raw_len, stored, crc)), payload) in
                                    outs.iter_mut().zip(frames).zip(payloads)
                                {
                                    decode_frame(payload, raw_len, stored, crc, out)?;
                                }
                                Ok(())
                            })
                        })
                        .collect();
                    handles
                        .into_iter()
                        .map(|h| h.join().expect("spill decode thread panicked"))
                        .collect()
                });
                if let Some(e) = decoded.into_iter().find_map(Result::err) {
                    let _ = tx.send(Err(e));
                    return;
                }
                for out in outs {
                    if !out.is_empty() && tx.send(Ok(out)).is_err() {
                        return;
                    }
                }
                if end {
                    return;
                }
            }
        });
        SpillReader {
            rx,
            recycle,
            out: Vec::new(),
            pos: 0,
            eof: false,
            worker: Some(worker),
        }
    }
}

impl BufRead for SpillReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.pos == self.out.len() {
            if self.eof {
                return Ok(&[]);
            }
            match self.rx.recv() {
                Ok(Ok(next)) => {
                    let used = std::mem::replace(&mut self.out, next);
                    let _ = self.recycle.send(used);
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => {
                    self.eof = true;
                    if let Some(w) = self.worker.take() {
                        w.join()
                            .map_err(|_| io::Error::other("spill reader thread panicked"))?;
                    }
                }
            }
        }
        Ok(&self.out[self.pos..])
    }

    fn consume(&mut self, amount: usize) {
        self.pos += amount.min(self.out.len() - self.pos);
    }
}

impl Read for SpillReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        let available = self.fill_buf()?;
        let n = buf.len().min(available.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

/// `path`'s decoded stream: a spill file ([`MAGIC`]) or a BGZF file,
/// inflated on `threads` threads either way.
pub fn open_decoded(path: &Path, threads: usize) -> io::Result<Box<dyn Read + Send>> {
    let mut file = std::fs::File::open(path)?;
    let mut magic = [0u8; 8];
    let n = read_up_to(&mut file, &mut magic)?;
    use std::io::Seek;
    file.seek(io::SeekFrom::Start(0))?;
    let file = io::BufReader::with_capacity(1 << 22, file);
    Ok(if n == magic.len() && magic.starts_with(MAGIC_FAMILY) {
        Box::new(SpillReader::new(file, threads))
    } else {
        Box::new(crate::bgzf::ParallelBgzfReader::new(file, threads))
    })
}

fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut got = 0;
    while got < buf.len() {
        match r.read(&mut buf[got..])? {
            0 => break,
            n => got += n,
        }
    }
    Ok(got)
}

/// Write `data` (a BAM's uncompressed stream) to `path` as a spill file.
pub fn write_file(path: &Path, data: &[u8], threads: usize) -> io::Result<()> {
    let mut w = SpillFile::create(path, threads)?;
    w.write_all(data)?;
    w.finish()
}

// ---------------------------------------------------------------------------
// Manifest

/// A directory's spill manifest: every chunk spill file of one run, written
/// once they are all complete. Spill files are only read back (`--resume-from
/// markdup`) with a manifest naming exactly the files present, and cleanup
/// removes the manifest before any file, so a partial set is never taken
/// for a whole one.
pub const MANIFEST: &str = "spill.manifest";
const MANIFEST_HEAD: &str = "nomics spill manifest 1";

fn file_name(p: &Path) -> io::Result<String> {
    p.file_name()
        .and_then(|n| n.to_str())
        .map(str::to_owned)
        .ok_or_else(|| invalid(format!("{}: no file name", p.display())))
}

/// Write `dir/`[`MANIFEST`] naming `files` (all in `dir`) with their sizes.
pub fn write_manifest(dir: &Path, files: &[PathBuf]) -> io::Result<()> {
    let mut text = format!("{MANIFEST_HEAD}\n");
    for f in files {
        if f.parent() != Some(dir) {
            return Err(invalid(format!(
                "{} is not in {}",
                f.display(),
                dir.display()
            )));
        }
        text += &format!("{}\t{}\n", file_name(f)?, std::fs::metadata(f)?.len());
    }
    let tmp = dir.join(format!("{MANIFEST}{PARTIAL_SUFFIX}"));
    std::fs::write(&tmp, text)?;
    std::fs::File::open(&tmp)?.sync_all()?;
    std::fs::rename(&tmp, dir.join(MANIFEST))
}

/// Check that `dir/`[`MANIFEST`] names exactly `files`, at their sizes.
pub fn check_manifest(dir: &Path, files: &[PathBuf]) -> io::Result<()> {
    let path = dir.join(MANIFEST);
    let text = std::fs::read_to_string(&path).map_err(|e| {
        invalid(format!(
            "{}: {e}; spill files without a manifest may be an incomplete set",
            path.display()
        ))
    })?;
    let mut lines = text.lines();
    if lines.next() != Some(MANIFEST_HEAD) {
        return Err(invalid(format!("{}: not a spill manifest", path.display())));
    }
    let mut listed = Vec::new();
    for line in lines {
        let (name, size) = line
            .split_once('\t')
            .and_then(|(n, s)| Some((n.to_owned(), s.parse::<u64>().ok()?)))
            .ok_or_else(|| invalid(format!("{}: bad line {line:?}", path.display())))?;
        listed.push((name, size));
    }
    let mut present = files
        .iter()
        .map(|f| Ok((file_name(f)?, std::fs::metadata(f)?.len())))
        .collect::<io::Result<Vec<_>>>()?;
    listed.sort();
    present.sort();
    if listed != present {
        let names = |v: &[(String, u64)]| {
            v.iter()
                .map(|(n, s)| format!("{n} ({s} B)"))
                .collect::<Vec<_>>()
                .join(", ")
        };
        return Err(invalid(format!(
            "{} lists {} spill files [{}] but {} are present [{}]: an incomplete set",
            path.display(),
            listed.len(),
            names(&listed),
            present.len(),
            names(&present)
        )));
    }
    Ok(())
}

/// Remove a run's spill: the manifest first (so an interrupted cleanup
/// leaves files no resume accepts), then `files`.
pub fn remove_spill(dir: &Path, files: &[PathBuf]) -> io::Result<()> {
    match std::fs::remove_file(dir.join(MANIFEST)) {
        Err(e) if e.kind() != io::ErrorKind::NotFound => return Err(e),
        _ => {}
    }
    for f in files {
        std::fs::remove_file(f)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(data: &[u8]) -> usize {
        let mut table = Vec::new();
        let mut c = Vec::new();
        compress(data, &mut table, &mut c);
        let mut d = vec![7u8; 3];
        decompress(&c, data.len(), &mut d).unwrap();
        assert_eq!(&d[3..], data, "round trip of {} bytes", data.len());
        c.len()
    }

    /// xorshift64*.
    fn rng(seed: u64) -> impl FnMut() -> u64 {
        let mut s = seed.max(1);
        move || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            s.wrapping_mul(0x2545_f491_4f6c_dd1d)
        }
    }

    /// Records shaped like BAM: fixed fields, a name, packed bases, quals.
    fn bam_like(n: usize, seed: u64) -> Vec<u8> {
        let mut r = rng(seed);
        let mut out = Vec::new();
        for i in 0..n {
            let mut rec = Vec::new();
            rec.extend_from_slice(&(i as u32 / 3).to_le_bytes());
            rec.extend_from_slice(&(1000 + 7 * i as u32).to_le_bytes());
            rec.extend_from_slice(&[60, 0x49, 0x12, 0x00, 1, 0, 0x63, 0]);
            rec.extend_from_slice(format!("A00{}:8:H{:x}:1:{}", i % 4, r() % 4096, i).as_bytes());
            for _ in 0..75 {
                rec.push([0x11, 0x12, 0x14, 0x18, 0x21, 0x42, 0x84, 0x88][(r() % 8) as usize]);
            }
            for _ in 0..150 {
                rec.push([37, 37, 37, 25, 11, 37][(r() % 6) as usize]);
            }
            rec.extend_from_slice(b"RGZHG002.L1\0NMC\x01ASC\x96XSC\x00");
            out.extend_from_slice(&(rec.len() as i32).to_le_bytes());
            out.extend_from_slice(&rec);
        }
        out
    }

    #[test]
    fn round_trips() {
        roundtrip(b"");
        roundtrip(b"a");
        roundtrip(b"abcdefghijklmnop");
        // Long runs: matches with offset 1 and lengths needing extension bytes.
        roundtrip(&[0u8; 100_000]);
        roundtrip(&b"ab".repeat(40_000));
        roundtrip(&b"abcdefg".repeat(10_000));
        // Literal runs needing extension bytes, then a match.
        let mut r = rng(3);
        let noise: Vec<u8> = (0..70_000).map(|_| r() as u8).collect();
        let mut mixed = noise.clone();
        mixed.extend_from_slice(&noise[..5000]);
        mixed.extend_from_slice(&noise[100..300]);
        roundtrip(&mixed);
        for len in [4, 5, 7, 8, 9, 12, 13, 15, 16, 19, 20, 271, 4096] {
            let v: Vec<u8> = (0..len).map(|i| (i % 5) as u8).collect();
            roundtrip(&v);
            roundtrip(&noise[..len]);
        }
        // Matches reaching back almost the whole window, and beyond it.
        let mut far = noise[..66_000].to_vec();
        far.extend_from_slice(&noise[1000..2000]);
        far.extend_from_slice(&noise[..500]);
        roundtrip(&far);
        let bam = bam_like(5000, 9);
        let c = roundtrip(&bam);
        assert!(
            c * 4 < bam.len() * 3,
            "BAM-like data compresses ({c} of {})",
            bam.len()
        );
    }

    #[test]
    fn corrupt_payloads_are_errors() {
        let bam = bam_like(500, 1);
        let mut table = Vec::new();
        let mut c = Vec::new();
        compress(&bam, &mut table, &mut c);
        let mut out = Vec::new();
        assert!(decompress(&c[..c.len() - 1], bam.len(), &mut out).is_err());
        assert!(decompress(&c, bam.len() + 1, &mut out.clone()).is_err());
        assert!(decompress(&c, bam.len() - 1, &mut Vec::new()).is_err());
        // Offset reaching before the block.
        assert!(decompress(&[0x00, 0x05, 0x00], 4, &mut Vec::new()).is_err());
    }

    #[test]
    fn files_round_trip_across_blocks_and_threads() {
        let dir = std::env::temp_dir().join(format!("nomics-spill-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut data = bam_like(12_000, 5);
        // Incompressible blocks are stored.
        let mut r = rng(11);
        data.extend((0..BLOCK + 17).map(|_| r() as u8));
        data.extend(bam_like(3000, 6));
        for (threads, len) in [
            (1, 0),
            (1, data.len()),
            (3, data.len()),
            (8, BLOCK),
            (4, BLOCK * 2 + 1),
        ] {
            let len = len.min(data.len());
            let path = dir.join(format!("t{threads}-{len}.{EXTENSION}"));
            // Odd write sizes cross block boundaries.
            let f = io::BufWriter::new(std::fs::File::create(&path).unwrap());
            let mut w = SpillWriter::new(f, threads);
            for piece in data[..len].chunks(77_777) {
                w.write_all(piece).unwrap();
            }
            w.finish().unwrap().flush().unwrap();
            let mut back = Vec::new();
            open_decoded(&path, threads)
                .unwrap()
                .read_to_end(&mut back)
                .unwrap();
            assert_eq!(back.len(), len);
            assert!(back == data[..len], "{threads} threads, {len} bytes");
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    fn read_all(path: &Path) -> io::Result<Vec<u8>> {
        let mut back = Vec::new();
        open_decoded(path, 2)?.read_to_end(&mut back)?;
        Ok(back)
    }

    #[test]
    fn damaged_files_are_errors() {
        let dir = std::env::temp_dir().join(format!("nomics-spill-bad-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let data = bam_like(12_000, 8);
        assert!(data.len() > 2 * BLOCK);
        let path = dir.join(format!("whole.{EXTENSION}"));
        write_file(&path, &data, 4).unwrap();
        assert!(read_all(&path).unwrap() == data);
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[bytes.len() - 8..], FOOTER_MAGIC);
        let cut = dir.join(format!("cut.{EXTENSION}"));
        let check = |bytes: &[u8], what: &str| {
            std::fs::write(&cut, bytes).unwrap();
            assert!(read_all(&cut).is_err(), "{what} must be an error");
        };
        // Truncated: inside the footer, at the end frame, inside a frame.
        check(&bytes[..bytes.len() - 1], "a cut footer");
        check(&bytes[..bytes.len() - FOOTER_LEN], "no footer");
        check(&bytes[..bytes.len() - FOOTER_LEN - 12], "no end frame");
        check(&bytes[..bytes.len() / 2], "a cut frame");
        // Trailing bytes after the footer.
        let mut extra = bytes.clone();
        extra.push(0);
        check(&extra, "a trailing byte");
        let mut twice = bytes.clone();
        twice.extend_from_slice(&bytes);
        check(&twice, "a second stream");
        // A zeroed frame header reads as an early end frame: the footer
        // check and the trailing frames catch it.
        let first = MAGIC.len();
        let second = first + 12 + (read_u32(&bytes, first + 4) & !STORED) as usize;
        let mut zeroed = bytes.clone();
        zeroed[second..second + 12].fill(0);
        check(&zeroed, "a zeroed frame header");
        // A footer that disagrees with the frames.
        let mut footer_off = bytes.clone();
        let at = footer_off.len() - FOOTER_LEN;
        footer_off[at] ^= 1;
        check(&footer_off, "a wrong decoded length");
        // A flipped payload byte fails the CRC or the decoder.
        let mut bad = bytes.clone();
        bad[MAGIC.len() + 12 + 1000] ^= 0x40;
        check(&bad, "a flipped payload byte");
        // An older format version is refused, not read as BGZF.
        let mut old = bytes.clone();
        old[..8].copy_from_slice(b"NOMSPL01");
        check(&old, "a version-1 file");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn unfinished_writers_leave_no_valid_file() {
        let dir = std::env::temp_dir().join(format!("nomics-spill-drop-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let data = bam_like(12_000, 4);
        // SpillFile dropped before finish (an error or a panic upstream):
        // no partial file left, and an older file of that name untouched
        // (only `finish`'s rename replaces it).
        let path = dir.join(format!("dropped.{EXTENSION}"));
        std::fs::write(&path, b"an older complete file").unwrap();
        {
            let mut w = SpillFile::create(&path, 3).unwrap();
            w.write_all(&data).unwrap();
        }
        let names: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(names, vec![std::ffi::OsString::from(format!("dropped.{EXTENSION}"))]);
        assert_eq!(std::fs::read(&path).unwrap(), b"an older complete file");
        // While being written, the new stream exists only under the
        // partial name; the old one keeps its name until the rename.
        let mut w = SpillFile::create(&path, 3).unwrap();
        w.write_all(&data).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"an older complete file");
        w.finish().unwrap();
        assert!(read_all(&path).unwrap() == data);
        // A bare SpillWriter dropped mid-stream writes no end frame or
        // footer, and buffered blocks are not passed off as the whole.
        let bare = dir.join(format!("bare.{EXTENSION}"));
        {
            let f = io::BufWriter::new(std::fs::File::create(&bare).unwrap());
            let mut w = SpillWriter::new(f, 2);
            w.write_all(&data).unwrap();
        }
        // The worker is detached on drop; wait for it to stop writing.
        let mut last = u64::MAX;
        for _ in 0..200 {
            std::thread::sleep(std::time::Duration::from_millis(10));
            let n = std::fs::metadata(&bare).unwrap().len();
            if n == last {
                break;
            }
            last = n;
        }
        assert!(read_all(&bare).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn manifests_refuse_partial_sets() {
        let dir = std::env::temp_dir().join(format!("nomics-spill-set-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let files: Vec<PathBuf> = (1..=3)
            .map(|i| {
                let p = dir.join(format!("S.{i:04}.{EXTENSION}"));
                write_file(&p, &bam_like(100 * i, i as u64), 2).unwrap();
                p
            })
            .collect();
        assert!(check_manifest(&dir, &files).is_err(), "no manifest yet");
        write_manifest(&dir, &files).unwrap();
        check_manifest(&dir, &files).unwrap();
        // A missing file (an interrupted cleanup), an extra one, a changed one.
        assert!(check_manifest(&dir, &files[..2]).is_err());
        let extra = dir.join(format!("S.0004.{EXTENSION}"));
        write_file(&extra, &bam_like(10, 1), 1).unwrap();
        let mut more = files.clone();
        more.push(extra);
        assert!(check_manifest(&dir, &more).is_err());
        write_file(&files[1], &bam_like(150, 9), 2).unwrap();
        assert!(check_manifest(&dir, &files).is_err());
        write_manifest(&dir, &files).unwrap();
        check_manifest(&dir, &files).unwrap();
        // Cleanup removes the manifest first: what an interrupted cleanup
        // leaves is refused.
        std::fs::remove_file(dir.join(MANIFEST)).unwrap();
        std::fs::remove_file(&files[0]).unwrap();
        assert!(check_manifest(&dir, &files[1..]).is_err());
        write_manifest(&dir, &files[1..]).unwrap();
        remove_spill(&dir, &files[1..]).unwrap();
        assert!(!dir.join(MANIFEST).exists() && !files[1].exists() && !files[2].exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn bam_stream_reads_a_spill_as_the_bam() {
        use crate::bam::{Header, RecordSet};
        use crate::bamstream::BamStream;
        let dir = std::env::temp_dir().join(format!("nomics-spill-bam-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let header = Header {
            text: b"@HD\tVN:1.6\tSO:coordinate\n".to_vec(),
            refs: Vec::new(),
        };
        let mut stream = Vec::new();
        header.write(&mut stream).unwrap();
        stream.extend(bam_like(20_000, 2));
        let set = RecordSet::from_stream(stream).unwrap();
        let (bam, spill) = (dir.join("c.bam"), dir.join(format!("c.{EXTENSION}")));
        set.write(&bam, 3).unwrap();
        set.write(&spill, 3).unwrap();
        assert_eq!(&std::fs::read(&spill).unwrap()[..8], MAGIC);
        let read = |p: &Path| {
            let mut s = BamStream::open(p, 777, 2, 2).unwrap();
            let mut records = Vec::new();
            while let Some(b) = s.next_batch().unwrap() {
                records.extend((0..b.len()).map(|i| b.record(i).to_vec()));
            }
            (s.header.text.clone(), records)
        };
        let (from_bam, from_spill) = (read(&bam), read(&spill));
        assert_eq!(from_bam.1.len(), 20_000);
        assert!(from_bam == from_spill);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
