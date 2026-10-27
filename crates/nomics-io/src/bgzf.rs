//! BGZF: concatenated gzip members of at most 64 KiB, each carrying its
//! compressed size in a `BC` extra subfield (SAMv1 §4.1).
//!
//! Block layout is part of the byte-identity contract, so the writer follows a
//! [`Profile`] taken from the reference writer:
//!
//! - [`Profile::HTSJDK`]: GATK/Picard (`BlockCompressedOutputStream`). Fills
//!   blocks to exactly 65,498 uncompressed bytes regardless of record
//!   boundaries; zlib level 2 (measured on the tier1 goldens).
//!
//! The deflate engine is pluggable ([`Deflate`]) so that a bit-exact zlib
//! replica can be swapped in and verified block by block against goldens.
//! BAM outputs choose between that replica and a faster deflater with the
//! same block cuts ([`BamCompression`]).

use std::io::{self, BufRead, Read, Write};

use crate::crc32::crc32;
use crate::deflate::{FastDeflater, HtsjdkDeflater};
use crate::inflate::{InflateError, inflate};

/// Fixed 18-byte BGZF member header up to and including the BSIZE field.
pub const HEADER_LEN: usize = 18;
/// CRC32 + ISIZE.
pub const FOOTER_LEN: usize = 8;
/// Largest legal member (BSIZE + 1).
pub const MAX_BLOCK_LEN: usize = 65536;
/// The empty member every BGZF file ends with.
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Profile {
    /// Uncompressed bytes per full block.
    pub block_payload: usize,
}

impl Profile {
    /// htsjdk `BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE`.
    pub const HTSJDK: Profile = Profile {
        block_payload: 65498,
    };
}

/// How a BAM output compresses its BGZF blocks. Both modes cut blocks at the
/// same 65,498-byte payloads ([`Profile::HTSJDK`]), so the uncompressed
/// stream, the records and the payload layout that `.bai` virtual offsets
/// depend on are identical; only the compressed bytes differ.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BamCompression {
    /// htsjdk's zlib level 2, byte for byte: with the reference's header, the
    /// file equals the reference's. For byte-comparison workflows.
    Htsjdk,
    /// [`crate::deflate::FastDeflater`]: greedy single-probe matching.
    Fast,
}

impl BamCompression {
    /// Environment variable selecting the mode for BAM outputs.
    pub const ENV: &'static str = "NOMICS_BAM_COMPRESSION";
    /// Mode when [`Self::ENV`] is unset. Fast: at 30x it compresses BAM
    /// blocks 2.3x faster than the htsjdk replica and writes files ~8%
    /// smaller; the gates are on records and the uncompressed stream.
    pub const DEFAULT: BamCompression = BamCompression::Fast;

    pub fn parse(s: &str) -> Option<BamCompression> {
        match s {
            "htsjdk" => Some(BamCompression::Htsjdk),
            "fast" => Some(BamCompression::Fast),
            _ => None,
        }
    }

    /// The mode for BAM outputs (md.bam, recal.bam, chunk BAMs):
    /// `NOMICS_BAM_COMPRESSION=htsjdk|fast`, else [`Self::DEFAULT`]. Read
    /// once per process. Any other value panics rather than silently
    /// writing something the caller did not ask for.
    pub fn from_env() -> BamCompression {
        static MODE: std::sync::OnceLock<BamCompression> = std::sync::OnceLock::new();
        *MODE.get_or_init(|| match std::env::var(Self::ENV) {
            Err(std::env::VarError::NotPresent) => Self::DEFAULT,
            Ok(v) => Self::parse(&v)
                .unwrap_or_else(|| panic!("{}={v:?}: expected \"htsjdk\" or \"fast\"", Self::ENV)),
            Err(e) => panic!("{}: {e}", Self::ENV),
        })
    }
}

/// A raw-DEFLATE compressor for one block's payload.
pub trait Deflate {
    /// Compress `input` as one complete raw DEFLATE stream (final block set),
    /// appending to `out`.
    fn deflate(&mut self, input: &[u8], out: &mut Vec<u8>);
}

/// Stored (uncompressed) DEFLATE blocks. Valid output, but not yet verified
/// against zlib level 0 or htsjdk's no-compression fallback.
#[derive(Debug, Default, Clone, Copy)]
pub struct Stored;

impl Deflate for Stored {
    fn deflate(&mut self, input: &[u8], out: &mut Vec<u8>) {
        let mut chunks = input.chunks(0xffff).peekable();
        if chunks.peek().is_none() {
            out.extend_from_slice(&[0x01, 0x00, 0x00, 0xff, 0xff]);
            return;
        }
        while let Some(chunk) = chunks.next() {
            let last = chunks.peek().is_none();
            let len = chunk.len() as u16;
            out.push(last as u8);
            out.extend_from_slice(&len.to_le_bytes());
            out.extend_from_slice(&(!len).to_le_bytes());
            out.extend_from_slice(chunk);
        }
    }
}

#[derive(Debug)]
pub enum BgzfError {
    Io(io::Error),
    /// Not a BGZF member header (bad magic, flags, or missing `BC` subfield).
    BadHeader {
        offset: u64,
    },
    Inflate {
        offset: u64,
        error: InflateError,
    },
    /// Decompressed size or CRC32 does not match the member footer.
    BadFooter {
        offset: u64,
    },
    /// The stream ended inside a member.
    Truncated {
        offset: u64,
    },
}

impl std::fmt::Display for BgzfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BgzfError::Io(e) => write!(f, "{e}"),
            BgzfError::BadHeader { offset } => write!(f, "bad BGZF header at byte {offset}"),
            BgzfError::Inflate { offset, error } => {
                write!(f, "inflate failed in block at byte {offset}: {error}")
            }
            BgzfError::BadFooter { offset } => {
                write!(f, "CRC32/ISIZE mismatch in block at byte {offset}")
            }
            BgzfError::Truncated { offset } => write!(f, "truncated BGZF block at byte {offset}"),
        }
    }
}

impl std::error::Error for BgzfError {}

impl From<io::Error> for BgzfError {
    fn from(e: io::Error) -> Self {
        BgzfError::Io(e)
    }
}

impl From<BgzfError> for io::Error {
    fn from(e: BgzfError) -> Self {
        match e {
            BgzfError::Io(e) => e,
            other => io::Error::new(io::ErrorKind::InvalidData, other),
        }
    }
}

/// Where a block sits in the compressed file and how big it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockInfo {
    /// Byte offset of the member in the compressed file.
    pub coffset: u64,
    /// Total member length including header and footer.
    pub csize: u32,
    /// Uncompressed payload length.
    pub isize: u32,
}

/// Read one whole compressed member into `raw` (replacing its contents).
/// Returns `Ok(false)` at a clean end of input.
pub fn read_raw_block<R: Read>(
    r: &mut R,
    raw: &mut Vec<u8>,
    coffset: u64,
) -> Result<bool, BgzfError> {
    raw.clear();
    raw.resize(HEADER_LEN, 0);
    let mut got = 0;
    while got < HEADER_LEN {
        match r.read(&mut raw[got..])? {
            0 if got == 0 => return Ok(false),
            0 => return Err(BgzfError::Truncated { offset: coffset }),
            n => got += n,
        }
    }
    let bsize = parse_header(raw, coffset)?;
    raw.resize(bsize, 0);
    r.read_exact(&mut raw[HEADER_LEN..])
        .map_err(|e| match e.kind() {
            io::ErrorKind::UnexpectedEof => BgzfError::Truncated { offset: coffset },
            _ => BgzfError::Io(e),
        })?;
    Ok(true)
}

/// Validate the fixed BGZF header and return the member length (BSIZE + 1).
fn parse_header(h: &[u8], offset: u64) -> Result<usize, BgzfError> {
    // ID1 ID2 CM FLG(FEXTRA) ... XLEN=6, SI1='B' SI2='C' SLEN=2.
    if h[0] != 0x1f
        || h[1] != 0x8b
        || h[2] != 8
        || h[3] & 4 == 0
        || h[10..16] != [6, 0, b'B', b'C', 2, 0]
    {
        return Err(BgzfError::BadHeader { offset });
    }
    let bsize = u16::from_le_bytes([h[16], h[17]]) as usize + 1;
    if bsize < HEADER_LEN + FOOTER_LEN {
        return Err(BgzfError::BadHeader { offset });
    }
    Ok(bsize)
}

/// Decode a whole member (as read by [`read_raw_block`]), appending the
/// payload to `out` after checking ISIZE and CRC32.
pub fn decode_block(raw: &[u8], coffset: u64, out: &mut Vec<u8>) -> Result<BlockInfo, BgzfError> {
    let bsize = parse_header(raw, coffset)?;
    if raw.len() != bsize {
        return Err(BgzfError::Truncated { offset: coffset });
    }
    let cdata = &raw[HEADER_LEN..bsize - FOOTER_LEN];
    let footer = &raw[bsize - FOOTER_LEN..];
    let crc = u32::from_le_bytes(footer[0..4].try_into().unwrap());
    let isize = u32::from_le_bytes(footer[4..8].try_into().unwrap());
    // A BGZF member never holds more than 64 KiB; don't trust ISIZE beyond that.
    if isize as usize > MAX_BLOCK_LEN {
        return Err(BgzfError::BadFooter { offset: coffset });
    }
    let start = out.len();
    out.reserve(isize as usize);
    inflate(cdata, out, isize as usize).map_err(|error| BgzfError::Inflate {
        offset: coffset,
        error,
    })?;
    let payload = &out[start..];
    if payload.len() != isize as usize || crc32(payload) != crc {
        return Err(BgzfError::BadFooter { offset: coffset });
    }
    Ok(BlockInfo {
        coffset,
        csize: bsize as u32,
        isize,
    })
}

/// Byte ranges of every BGZF member in `data` (headers validated, payloads
/// not decoded).
pub fn split_blocks(data: &[u8]) -> Result<Vec<std::ops::Range<usize>>, BgzfError> {
    let mut blocks = Vec::new();
    let mut at = 0;
    while at < data.len() {
        let h = data
            .get(at..at + HEADER_LEN)
            .ok_or(BgzfError::Truncated { offset: at as u64 })?;
        let bsize = parse_header(h, at as u64)?;
        if at + bsize > data.len() {
            return Err(BgzfError::Truncated { offset: at as u64 });
        }
        blocks.push(at..at + bsize);
        at += bsize;
    }
    Ok(blocks)
}

/// Decode a whole in-memory BGZF file on `threads` threads, handing the
/// decompressed payload to `on_batch` in file order, one batch of
/// consecutive blocks at a time. Memory is bounded by the batch size.
pub fn decode_parallel(
    data: &[u8],
    threads: usize,
    on_batch: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let threads = threads.max(1);
    let blocks = split_blocks(data)?;
    // ~64 blocks (~4 MiB decoded) per thread per batch.
    let per_thread = 64;
    let mut outs: Vec<Vec<u8>> = vec![Vec::new(); threads];
    for batch in blocks.chunks(threads * per_thread) {
        decode_chunks(data, batch, &mut outs, per_thread)?;
        for out in outs.iter().take(batch.len().div_ceil(per_thread)) {
            on_batch(out)?;
        }
    }
    Ok(())
}

// Each output owns consecutive decoded blocks. Keep these buffers separate so
// slice consumers do not pay for concatenating the entire decoded segment.
fn decode_chunks(
    data: &[u8],
    blocks: &[std::ops::Range<usize>],
    outs: &mut [Vec<u8>],
    per_thread: usize,
) -> io::Result<()> {
    let results: Vec<Result<(), BgzfError>> = std::thread::scope(|scope| {
        let handles: Vec<_> = outs
            .iter_mut()
            .zip(blocks.chunks(per_thread))
            .map(|(out, chunk)| {
                scope.spawn(move || -> Result<(), BgzfError> {
                    out.clear();
                    for r in chunk {
                        decode_block(&data[r.clone()], r.start as u64, out)?;
                    }
                    Ok(())
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("decode thread panicked"))
            .collect()
    });
    for r in results {
        r?;
    }
    Ok(())
}

/// BGZF decoder implementing [`Read`] and [`BufRead`] for streaming large files: a
/// background thread reads segments of `threads * 64` blocks and inflates
/// each on `threads` threads while the consumer reads the previous one. Decoded
/// buffers are handed off without concatenation and recycled after consumption.
pub struct ParallelBgzfReader {
    rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    recycle: std::sync::mpsc::Sender<Vec<u8>>,
    out: Vec<u8>,
    pos: usize,
    eof: bool,
    worker: Option<std::thread::JoinHandle<()>>,
}

impl ParallelBgzfReader {
    pub fn new<R: Read + Send + 'static>(mut inner: R, threads: usize) -> Self {
        let threads = threads.max(1);
        let (tx, rx) = std::sync::mpsc::sync_channel::<io::Result<Vec<u8>>>(threads * 2);
        let (recycle, recycled) = std::sync::mpsc::channel::<Vec<u8>>();
        let worker = std::thread::spawn(move || {
            let (mut raw, mut segment, mut coffset) = (Vec::new(), Vec::new(), 0u64);
            let mut outs = vec![Vec::new(); threads];
            loop {
                segment.clear();
                let (mut blocks, mut eof) = (0, false);
                while blocks < threads * 64 {
                    match read_raw_block(&mut inner, &mut raw, coffset) {
                        Ok(true) => {}
                        Ok(false) => {
                            eof = true;
                            break;
                        }
                        Err(e) => {
                            let _ = tx.send(Err(e.into()));
                            return;
                        }
                    }
                    coffset += raw.len() as u64;
                    segment.extend_from_slice(&raw);
                    blocks += 1;
                }
                let decoded = split_blocks(&segment)
                    .map_err(io::Error::from)
                    .and_then(|ranges| {
                        for out in &mut outs {
                            *out = recycled.try_recv().unwrap_or_default();
                            out.clear();
                        }
                        decode_chunks(&segment, &ranges, &mut outs, 64)
                    });
                if let Err(e) = decoded {
                    let _ = tx.send(Err(e));
                    return;
                }
                for out in &mut outs {
                    if !out.is_empty() && tx.send(Ok(std::mem::take(out))).is_err() {
                        return;
                    }
                }
                if eof {
                    return;
                }
            }
        });
        ParallelBgzfReader {
            rx,
            recycle,
            out: Vec::new(),
            pos: 0,
            eof: false,
            worker: Some(worker),
        }
    }
}

impl BufRead for ParallelBgzfReader {
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
                            .map_err(|_| io::Error::other("BGZF reader thread panicked"))?;
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

impl Read for ParallelBgzfReader {
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

/// Sequential BGZF decoder implementing [`Read`].
pub struct BgzfReader<R> {
    inner: R,
    raw: Vec<u8>,
    block: Vec<u8>,
    pos: usize,
    next_coffset: u64,
    current: Option<BlockInfo>,
}

impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R) -> Self {
        BgzfReader {
            inner,
            raw: Vec::new(),
            block: Vec::new(),
            pos: 0,
            next_coffset: 0,
            current: None,
        }
    }

    /// Decode the next member into the internal buffer. `Ok(None)` at EOF.
    pub fn next_block(&mut self) -> Result<Option<BlockInfo>, BgzfError> {
        let coffset = self.next_coffset;
        if !read_raw_block(&mut self.inner, &mut self.raw, coffset)? {
            self.current = None;
            return Ok(None);
        }
        self.block.clear();
        self.pos = 0;
        let info = decode_block(&self.raw, coffset, &mut self.block)?;
        self.next_coffset += info.csize as u64;
        self.current = Some(info);
        Ok(Some(info))
    }

    /// Payload of the most recently decoded block.
    pub fn block(&self) -> &[u8] {
        &self.block
    }

    /// BAM virtual offset of the next byte `read` will return.
    pub fn virtual_offset(&self) -> u64 {
        match self.current {
            Some(info) if self.pos < self.block.len() => (info.coffset << 16) | self.pos as u64,
            _ => self.next_coffset << 16,
        }
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // Empty members (including the EOF marker) are skipped.
        while self.pos == self.block.len() {
            if self.next_block()?.is_none() {
                return Ok(0);
            }
        }
        let n = buf.len().min(self.block.len() - self.pos);
        buf[..n].copy_from_slice(&self.block[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Append one complete BGZF member (header, deflated `payload`, CRC32, ISIZE).
fn encode_member<D: Deflate>(deflater: &mut D, payload: &[u8], out: &mut Vec<u8>) {
    let start = out.len();
    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0, 0, 0,
    ]);
    deflater.deflate(payload, out);
    if out.len() - start + FOOTER_LEN > MAX_BLOCK_LEN {
        // Only reachable for compressed sizes 65,511-65,517 (HtsjdkDeflater
        // already falls back to stored at >= 65,518). htsjdk would emit a
        // block whose BSIZE overflows u16, i.e. a corrupt file; we store instead.
        // FastDeflater never gets here: for payloads up to 65,535 bytes its
        // output is at most the payload + 5 (one stored block).
        out.truncate(start + HEADER_LEN);
        Stored.deflate(payload, out);
    }
    let bsize = (out.len() - start + FOOTER_LEN - 1) as u16;
    out[start + 16..start + 18].copy_from_slice(&bsize.to_le_bytes());
    out.extend_from_slice(&crc32(payload).to_le_bytes());
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
}

/// BGZF encoder. Call [`BgzfWriter::finish`] to flush and append the EOF member;
/// dropping without it leaves the file without the EOF marker.
pub struct BgzfWriter<W: Write, D: Deflate> {
    inner: W,
    deflater: D,
    profile: Profile,
    buf: Vec<u8>,
    cdata: Vec<u8>,
    coffset: u64,
}

impl<W: Write, D: Deflate> BgzfWriter<W, D> {
    pub fn new(inner: W, deflater: D, profile: Profile) -> Self {
        BgzfWriter {
            inner,
            deflater,
            profile,
            buf: Vec::with_capacity(profile.block_payload),
            cdata: Vec::with_capacity(MAX_BLOCK_LEN),
            coffset: 0,
        }
    }

    /// BAM virtual offset of the next byte to be written.
    pub fn virtual_offset(&self) -> u64 {
        (self.coffset << 16) | self.buf.len() as u64
    }

    fn emit_block(&mut self) -> io::Result<()> {
        self.cdata.clear();
        encode_member(&mut self.deflater, &self.buf, &mut self.cdata);
        self.inner.write_all(&self.cdata)?;
        self.coffset += self.cdata.len() as u64;
        self.buf.clear();
        Ok(())
    }

    /// Write any partial block, then the EOF member. Returns the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        if !self.buf.is_empty() {
            self.emit_block()?;
        }
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()?;
        Ok(self.inner)
    }
}

impl<W: Write, D: Deflate> Write for BgzfWriter<W, D> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let room = self.profile.block_payload - self.buf.len();
        let n = data.len().min(room);
        self.buf.extend_from_slice(&data[..n]);
        if self.buf.len() == self.profile.block_payload {
            self.emit_block()?;
        }
        Ok(n)
    }

    /// Pushes completed blocks to the inner writer. Does not cut a partial
    /// block: block boundaries are part of the byte-identity contract.
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Multi-threaded BGZF writer, byte-identical to
/// `BgzfWriter<_, HtsjdkDeflater>` with `Profile::HTSJDK`.
///
/// htsjdk reuses one deflater per file, and `deflateReset` keeps the window
/// and `high_water`, so in general a block's output can depend on bytes left
/// over from earlier blocks. Here it cannot, except for the last block:
///
/// - Every full block's 65,498 input bytes are copied into window bytes
///   `[0, 65498)` by the first `fill_window`. The tail `[65498, 65536)` is
///   zeroed by `high_water` during the first block and never written by any
///   full block afterwards (post-slide reads land below 65,498). So each full
///   block sees exactly the window a fresh deflater would.
/// - Bytes past the end of the data after a mid-block slide are this block's
///   own pre-slide data, identical in both cases.
/// - The final partial block of `n` bytes can read stale bytes in `[n, 65536)`
///   left by the previous full block. It is compressed by a deflater primed
///   with that previous block, which by the argument above reproduces the
///   reused deflater's state exactly.
///
/// Measured: fresh-state compression matched all 668,293 full blocks of the
/// tier1 goldens; full-file rewrites through this writer are compared with
/// `examples/bam_roundtrip --rewrite --parallel`.
///
/// With [`BamCompression::Fast`] ([`ParallelBgzfWriter::with_compression`])
/// the block cuts are the same but each block is compressed independently by
/// a [`FastDeflater`], which has no state to carry.
pub struct ParallelBgzfWriter<W: Write + Send + 'static> {
    threads: usize,
    mode: BamCompression,
    batch: Vec<Vec<u8>>,
    buf: Vec<u8>,
    /// Full-block batches to the compressing thread, which owns the inner
    /// writer: batch N is compressed and written while the caller fills N+1.
    tx: Option<std::sync::mpsc::SyncSender<Vec<Vec<u8>>>>,
    worker: Option<std::thread::JoinHandle<io::Result<(W, Option<Vec<u8>>)>>>,
}

/// One compressing thread's deflater.
enum BlockDeflater {
    Htsjdk(HtsjdkDeflater),
    Fast(FastDeflater),
}

impl BlockDeflater {
    fn new(mode: BamCompression) -> Self {
        match mode {
            BamCompression::Htsjdk => BlockDeflater::Htsjdk(HtsjdkDeflater::new()),
            BamCompression::Fast => BlockDeflater::Fast(FastDeflater::new()),
        }
    }

    /// The fresh state a full block is compressed from (see the type docs).
    fn prepare(&mut self) {
        if let BlockDeflater::Htsjdk(d) = self {
            d.reinit();
        }
    }
}

impl Deflate for BlockDeflater {
    fn deflate(&mut self, input: &[u8], out: &mut Vec<u8>) {
        match self {
            BlockDeflater::Htsjdk(d) => d.deflate(input, out),
            BlockDeflater::Fast(d) => d.deflate(input, out),
        }
    }
}

impl<W: Write + Send + 'static> ParallelBgzfWriter<W> {
    /// htsjdk-exact output ([`BamCompression::Htsjdk`]).
    pub fn new(inner: W, threads: usize) -> Self {
        Self::with_compression(inner, threads, BamCompression::Htsjdk)
    }

    /// For BAM outputs: the mode from [`BamCompression::from_env`].
    pub fn for_bam(inner: W, threads: usize) -> Self {
        Self::with_compression(inner, threads, BamCompression::from_env())
    }

    pub fn with_compression(inner: W, threads: usize, mode: BamCompression) -> Self {
        let threads = threads.max(1);
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Vec<u8>>>(1);
        let worker = std::thread::spawn(move || -> io::Result<(W, Option<Vec<u8>>)> {
            let mut inner = inner;
            let mut deflaters: Vec<BlockDeflater> =
                (0..threads).map(|_| BlockDeflater::new(mode)).collect();
            let mut members: Vec<Vec<u8>> = (0..threads * 8).map(|_| Vec::new()).collect();
            let mut last_full = None;
            for mut blocks in rx {
                let n = blocks.len();
                let per = n.div_ceil(threads);
                std::thread::scope(|scope| {
                    for ((members, deflater), chunk) in members[..n]
                        .chunks_mut(per)
                        .zip(deflaters.iter_mut())
                        .zip(blocks.chunks(per))
                    {
                        scope.spawn(move || {
                            for (out, payload) in members.iter_mut().zip(chunk) {
                                out.clear();
                                deflater.prepare();
                                encode_member(deflater, payload, out);
                            }
                        });
                    }
                });
                for m in &members[..n] {
                    inner.write_all(m)?;
                }
                last_full = blocks.pop();
            }
            Ok((inner, last_full))
        });
        ParallelBgzfWriter {
            threads,
            mode,
            batch: Vec::new(),
            buf: Vec::with_capacity(Profile::HTSJDK.block_payload),
            tx: Some(tx),
            worker: Some(worker),
        }
    }

    fn batch_len(&self) -> usize {
        self.threads * 8
    }

    /// Hand the queued full blocks to the compressing thread.
    fn flush_batch(&mut self) -> io::Result<()> {
        if self.batch.is_empty() {
            return Ok(());
        }
        let blocks = std::mem::take(&mut self.batch);
        let sent = self.tx.as_ref().map(|tx| tx.send(blocks));
        if !matches!(sent, Some(Ok(()))) {
            // The worker stopped on a write error; surface it.
            return Err(self
                .join()
                .err()
                .unwrap_or_else(|| io::Error::other("BGZF compression thread stopped")));
        }
        Ok(())
    }

    fn join(&mut self) -> io::Result<(W, Option<Vec<u8>>)> {
        drop(self.tx.take());
        match self.worker.take() {
            Some(w) => w
                .join()
                .map_err(|_| io::Error::other("BGZF compression thread panicked"))?,
            None => Err(io::Error::other("BGZF compression thread already joined")),
        }
    }

    /// Write the final partial block and the EOF member; returns the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.flush_batch()?;
        let (mut inner, last_full) = self.join()?;
        if !self.buf.is_empty() {
            let mut d = BlockDeflater::new(self.mode);
            // htsjdk's last partial block can see the previous block's
            // window bytes: prime a fresh deflater with that block.
            if let (BlockDeflater::Htsjdk(h), Some(prev)) = (&mut d, &last_full) {
                let mut scratch = Vec::new();
                h.deflate(prev, &mut scratch);
            }
            let mut member = Vec::new();
            encode_member(&mut d, &self.buf, &mut member);
            inner.write_all(&member)?;
        }
        inner.write_all(&EOF_BLOCK)?;
        inner.flush()?;
        Ok(inner)
    }
}

impl<W: Write + Send + 'static> Write for ParallelBgzfWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let room = Profile::HTSJDK.block_payload - self.buf.len();
        let n = data.len().min(room);
        self.buf.extend_from_slice(&data[..n]);
        if self.buf.len() == Profile::HTSJDK.block_payload {
            let full = std::mem::replace(
                &mut self.buf,
                Vec::with_capacity(Profile::HTSJDK.block_payload),
            );
            self.batch.push(full);
            if self.batch.len() == self.batch_len() {
                self.flush_batch()?;
            }
        }
        Ok(n)
    }

    /// Does not cut a partial block (block boundaries are part of the
    /// contract) and does not wait for queued blocks; `finish` does.
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(data: &[u8], profile: Profile) -> (Vec<u8>, Vec<BlockInfo>) {
        let mut w = BgzfWriter::new(Vec::new(), Stored, profile);
        w.write_all(data).unwrap();
        let file = w.finish().unwrap();
        let mut r = BgzfReader::new(&file[..]);
        let mut blocks = Vec::new();
        while let Some(info) = r.next_block().unwrap() {
            blocks.push(info);
        }
        let mut out = Vec::new();
        BgzfReader::new(&file[..]).read_to_end(&mut out).unwrap();
        assert_eq!(out, data);
        (file, blocks)
    }

    #[test]
    fn htsjdk_profile_fills_blocks_exactly() {
        let data: Vec<u8> = (0..200_000u32).map(|i| (i % 253) as u8).collect();
        let (file, blocks) = roundtrip(&data, Profile::HTSJDK);
        let sizes: Vec<u32> = blocks.iter().map(|b| b.isize).collect();
        assert_eq!(sizes, [65498, 65498, 65498, 3506, 0]);
        assert!(file.ends_with(&EOF_BLOCK));
        let total: u64 = blocks.iter().map(|b| b.csize as u64).sum();
        assert_eq!(total, file.len() as u64);
    }

    #[test]
    fn parallel_writer_matches_sequential_htsjdk_writer() {
        // Mixed compressible/incompressible content; sizes cover empty, < 1
        // block, exact multiples, and partial last blocks after many batches.
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut data = Vec::new();
        while data.len() < 3_000_000 {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            let r = x.wrapping_mul(0x2545_F491_4F6C_DD1D);
            if r & 3 == 0 && data.len() > 1000 {
                let d = 1 + (r >> 8) as usize % 900;
                for _ in 0..(r >> 40) as usize % 200 {
                    data.push(data[data.len() - d]);
                }
            } else {
                data.push(b"ACGTN"[(r >> 16) as usize % 5]);
                if r & 0x30 == 0 {
                    data.push((r >> 24) as u8);
                }
            }
        }
        let b = Profile::HTSJDK.block_payload;
        for len in [0, 100, b, b + 1, 2 * b, 3 * b - 7, 40 * b + 123, data.len()] {
            let input = &data[..len];
            let mut seq = BgzfWriter::new(Vec::new(), HtsjdkDeflater::new(), Profile::HTSJDK);
            seq.write_all(input).unwrap();
            let seq = seq.finish().unwrap();
            for threads in [1, 3, 8] {
                let mut par = ParallelBgzfWriter::new(Vec::new(), threads);
                // Uneven write sizes exercise the block-filling path.
                for piece in input.chunks(7919) {
                    par.write_all(piece).unwrap();
                }
                assert!(par.finish().unwrap() == seq, "len {len} threads {threads}");
            }
        }
    }

    #[test]
    fn fast_mode_keeps_block_cuts_and_roundtrips() {
        // BAM-like mix: repeated record skeletons, random bases and qualities.
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut data = Vec::new();
        while data.len() < 2_500_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            data.extend_from_slice(b"\x00\x01\x02\x03readname:1234:");
            data.extend_from_slice(&x.to_le_bytes());
            for k in 0..150 {
                data.push(b"ACGT"[((x >> (k % 60)) & 3) as usize]);
            }
            for k in 0..150u64 {
                data.push(b'!' + ((x.rotate_left(k as u32) & 31) as u8));
            }
        }
        let b = Profile::HTSJDK.block_payload;
        for len in [0, 1, b - 1, b, b + 1, 7 * b + 5, data.len()] {
            let input = &data[..len];
            let mut seq = BgzfWriter::new(Vec::new(), FastDeflater::new(), Profile::HTSJDK);
            seq.write_all(input).unwrap();
            let seq = seq.finish().unwrap();
            for threads in [1, 3, 8] {
                let mut par =
                    ParallelBgzfWriter::with_compression(Vec::new(), threads, BamCompression::Fast);
                for piece in input.chunks(7919) {
                    par.write_all(piece).unwrap();
                }
                assert!(par.finish().unwrap() == seq, "len {len} threads {threads}");
            }
            assert!(seq.ends_with(&EOF_BLOCK));
            let mut r = BgzfReader::new(&seq[..]);
            let mut sizes = Vec::new();
            while let Some(info) = r.next_block().unwrap() {
                assert!(info.csize as usize <= MAX_BLOCK_LEN);
                sizes.push(info.isize as usize);
            }
            let mut want: Vec<usize> = input.chunks(b).map(|c| c.len()).collect();
            want.push(0);
            assert_eq!(sizes, want, "len {len}");
            let mut out = Vec::new();
            BgzfReader::new(&seq[..]).read_to_end(&mut out).unwrap();
            assert!(out == input);
        }
        // Incompressible payloads still fit one member each.
        let noise: Vec<u8> = (0..3 * b as u64)
            .map(|i| (i.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 56) as u8)
            .collect();
        let mut w = ParallelBgzfWriter::with_compression(Vec::new(), 2, BamCompression::Fast);
        w.write_all(&noise).unwrap();
        let file = w.finish().unwrap();
        let mut got = Vec::new();
        BgzfReader::new(&file[..]).read_to_end(&mut got).unwrap();
        assert!(got == noise);
    }

    #[test]
    fn compression_mode_names() {
        assert_eq!(BamCompression::parse("fast"), Some(BamCompression::Fast));
        assert_eq!(
            BamCompression::parse("htsjdk"),
            Some(BamCompression::Htsjdk)
        );
        assert_eq!(BamCompression::parse("Fast"), None);
    }

    #[test]
    fn parallel_decode_matches_input() {
        let payload: Vec<u8> = (0..5_000_000u32)
            .map(|i| (i.wrapping_mul(2654435761) >> 11) as u8 % 7 + b'0')
            .collect();
        let mut w = ParallelBgzfWriter::new(Vec::new(), 4);
        w.write_all(&payload).unwrap();
        let file = w.finish().unwrap();
        for threads in [1, 3, 16] {
            let mut got = Vec::new();
            decode_parallel(&file, threads, &mut |b| {
                got.extend_from_slice(b);
                Ok(())
            })
            .unwrap();
            assert!(got == payload, "threads {threads}");
        }
        let mut truncated = file.clone();
        truncated.truncate(file.len() - 40);
        assert!(decode_parallel(&truncated, 2, &mut |_| Ok(())).is_err());
    }

    #[test]
    fn parallel_reader_preserves_slices_partial_batches_and_errors() {
        // Many small blocks exercise several complete segments and a partial
        // final segment without a large test allocation. Empty members do not
        // terminate the stream.
        let payload: Vec<u8> = (0..(64 * 5 * 3 + 7) * 257).map(|i| i as u8).collect();
        let mut writer = BgzfWriter::new(Vec::new(), Stored, Profile { block_payload: 257 });
        writer.write_all(&payload).unwrap();
        let mut file = writer.finish().unwrap();
        file.splice(0..0, EOF_BLOCK);
        for threads in [1, 3, 5] {
            let mut reader = ParallelBgzfReader::new(std::io::Cursor::new(file.clone()), threads);
            // An empty Read must not consume or wait for a decoded buffer.
            assert_eq!(reader.read(&mut []).unwrap(), 0);
            let mut at = 0;
            loop {
                let bytes = reader.fill_buf().unwrap();
                if bytes.is_empty() {
                    break;
                }
                let n = bytes.len().min(7919);
                assert_eq!(&bytes[..n], &payload[at..at + n]);
                at += n;
                reader.consume(n);
            }
            assert_eq!(at, payload.len());
            assert!(reader.fill_buf().unwrap().is_empty());
            let mut decoded = Vec::new();
            ParallelBgzfReader::new(std::io::Cursor::new(file.clone()), threads)
                .read_to_end(&mut decoded)
                .unwrap();
            assert_eq!(decoded, payload);
        }
        let mut corrupt = file.clone();
        let crc_at = corrupt.len() - EOF_BLOCK.len() - 8;
        corrupt[crc_at] ^= 1;
        for bad in [corrupt, file[..file.len() - 5].to_vec()] {
            let mut decoded = Vec::new();
            assert!(
                ParallelBgzfReader::new(std::io::Cursor::new(bad), 3)
                    .read_to_end(&mut decoded)
                    .is_err()
            );
        }
    }

    #[test]
    fn empty_input_is_just_eof() {
        let (file, blocks) = roundtrip(&[], Profile::HTSJDK);
        assert_eq!(file, EOF_BLOCK);
        assert_eq!(blocks.len(), 1);
    }

    #[test]
    fn corrupt_crc_is_detected() {
        let (mut file, _) = roundtrip(b"some payload", Profile::HTSJDK);
        let crc_at = file.len() - EOF_BLOCK.len() - 8;
        file[crc_at] ^= 1;
        let err = BgzfReader::new(&file[..]).next_block().unwrap_err();
        assert!(matches!(err, BgzfError::BadFooter { offset: 0 }));
    }

    #[test]
    fn virtual_offsets_track_blocks() {
        let mut w = BgzfWriter::new(Vec::new(), Stored, Profile::HTSJDK);
        assert_eq!(w.virtual_offset(), 0);
        w.write_all(&[7u8; 65498 + 10]).unwrap();
        let first_block_len = 18 + 5 + 65498 + 8;
        assert_eq!(w.virtual_offset(), ((first_block_len as u64) << 16) | 10);
    }
}
