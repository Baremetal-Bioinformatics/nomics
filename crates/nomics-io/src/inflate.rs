//! Raw DEFLATE decoder (RFC 1951).
//!
//! Packed primary/subtable Huffman lookup, word-sized bit refill, and safe
//! bulk match copies. The input and output tails use the checked scalar path.

use std::fmt;
use std::io;
use std::sync::OnceLock;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InflateError {
    /// Input ended before the final block was complete.
    Truncated,
    /// Reserved block type 3.
    BadBlockType,
    /// Stored block LEN/NLEN mismatch.
    BadStoredLength,
    /// Code lengths describe an over-subscribed or empty Huffman code.
    BadCode,
    /// A decoded symbol is invalid in context (e.g. length 286/287, distance 30/31).
    BadSymbol,
    /// A back-reference points before the start of this stream's output.
    DistanceTooFar,
    /// The stream would decode to more than the caller's limit.
    OutputTooLarge,
    /// The streaming sink returned an error (see [`inflate_to`]).
    Sink,
}

impl fmt::Display for InflateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl std::error::Error for InflateError {}

const FAST_BITS: u32 = 10;
const MAX_BITS: usize = 15;

const LEN_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LEN_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
/// Order in which code-length code lengths are transmitted.
const CLEN_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
    buf: u64,
    nbits: u32,
}

impl<'a> Bits<'a> {
    fn new(data: &'a [u8]) -> Self {
        Bits {
            data,
            pos: 0,
            buf: 0,
            nbits: 0,
        }
    }

    /// Refill with one little-endian word when there is input margin. Only
    /// whole bytes that fit are committed; the overlapping load is intentional.
    /// At the tail, zero padding permits table lookahead, but `consumed` must
    /// validate every decoded symbol before it can produce output.
    #[inline]
    fn refill(&mut self) {
        if let Some(word) = self.data.get(self.pos..).and_then(|s| s.get(..8)) {
            let word = u64::from_le_bytes(word.try_into().unwrap());
            let bytes = (63 - self.nbits) / 8;
            self.buf |= word << self.nbits;
            self.pos += bytes as usize;
            self.nbits += bytes * 8;
        } else {
            while self.nbits <= 56 {
                let b = self.data.get(self.pos).copied().unwrap_or(0);
                self.buf |= (b as u64) << self.nbits;
                self.nbits += 8;
                self.pos += 1;
            }
        }
    }

    #[inline]
    fn peek(&mut self, n: u32) -> u32 {
        if self.nbits < n {
            self.refill();
        }
        (self.buf & ((1u64 << n) - 1)) as u32
    }

    #[inline]
    fn consume(&mut self, n: u32) {
        self.buf >>= n;
        self.nbits -= n;
    }

    #[inline]
    fn bits(&mut self, n: u32) -> u32 {
        if n == 0 {
            return 0;
        }
        self.peek(n);
        self.buffered_bits(n)
    }

    #[inline]
    fn buffered_bits(&mut self, n: u32) -> u32 {
        let v = (self.buf & ((1u64 << n) - 1)) as u32;
        self.consume(n);
        v
    }

    fn align_to_byte(&mut self) {
        let r = self.nbits % 8;
        self.consume(r);
    }

    /// Bytes of real input consumed so far, or `Truncated` if decoding read
    /// past the end.
    fn consumed(&self) -> Result<usize, InflateError> {
        let used = self.pos - (self.nbits / 8) as usize;
        if used > self.data.len() {
            Err(InflateError::Truncated)
        } else {
            Ok(used)
        }
    }
}

// Entry layout: value:16, reserved:6, kind:2, extra:4, code length:4.
// Subtable pointers use kind 3 and store their width in `extra`.
const LENGTH: u32 = 1 << 8;
const END: u32 = 2 << 8;
const LINK: u32 = 3 << 8;
#[derive(Clone, Copy)]
enum Alphabet {
    Symbols,
    Litlen,
    Distance,
}

struct Huffman {
    table: Vec<u32>,
}

impl Huffman {
    fn new(lengths: &[u8], alphabet: Alphabet) -> Result<Self, InflateError> {
        let mut counts = [0u16; MAX_BITS + 1];
        for &len in lengths {
            let count = counts.get_mut(len as usize).ok_or(InflateError::BadCode)?;
            *count += 1;
        }
        counts[0] = 0;
        let mut left = 1i32;
        for &count in &counts[1..] {
            left = (left << 1) - count as i32;
            if left < 0 {
                return Err(InflateError::BadCode);
            }
        }
        // An empty distance alphabet is permitted in literal-only blocks.
        // Unassigned entries remain zero and fail if actually decoded.
        let mut next = [0u32; MAX_BITS + 1];
        let mut code = 0;
        for len in 1..=MAX_BITS {
            code = (code + u32::from(counts[len - 1])) << 1;
            next[len] = code;
        }
        let mut codes = Vec::with_capacity(lengths.len());
        let mut widths = [0u8; 1 << FAST_BITS];
        for (symbol, &len) in lengths.iter().enumerate() {
            if len == 0 {
                continue;
            }
            let code = next[len as usize].reverse_bits() >> (32 - len);
            next[len as usize] += 1;
            if u32::from(len) > FAST_BITS {
                let prefix = code as usize & ((1 << FAST_BITS) - 1);
                widths[prefix] = widths[prefix].max(len - FAST_BITS as u8);
            }
            let (value, extra, kind) = match alphabet {
                Alphabet::Litlen if (257..=285).contains(&symbol) => {
                    (LEN_BASE[symbol - 257], LEN_EXTRA[symbol - 257], LENGTH)
                }
                Alphabet::Litlen if symbol == 256 => (0, 0, END),
                Alphabet::Litlen if symbol > 285 => (0, 0, LINK),
                Alphabet::Distance if symbol < 30 => {
                    (DIST_BASE[symbol], DIST_EXTRA[symbol], LENGTH)
                }
                Alphabet::Distance => (0, 0, LINK),
                _ => (symbol as u16, 0, 0),
            };
            let entry = (u32::from(value) << 16) | kind | (u32::from(extra) << 4) | u32::from(len);
            codes.push((code as usize, len as u32, entry));
        }
        let mut table = vec![0; 1 << FAST_BITS];
        for (prefix, &width) in widths.iter().enumerate() {
            if width != 0 {
                table[prefix] = ((table.len() as u32) << 16) | LINK | (u32::from(width) << 4);
                table.resize(table.len() + (1 << width), 0);
            }
        }
        for (code, len, entry) in codes {
            let (base, start, stride, size) = if len <= FAST_BITS {
                (0, code, 1 << len, 1 << FAST_BITS)
            } else {
                let prefix = code & ((1 << FAST_BITS) - 1);
                (
                    (table[prefix] >> 16) as usize,
                    code >> FAST_BITS,
                    1 << (len - FAST_BITS),
                    1 << widths[prefix],
                )
            };
            for i in (start..size).step_by(stride) {
                table[base + i] = entry;
            }
        }
        Ok(Self { table })
    }

    /// Return a packed literal, length, distance, or end entry. Lookahead may
    /// include padding; callers check actual consumption at the input tail.
    #[inline]
    fn decode(&self, br: &mut Bits) -> Result<u32, InflateError> {
        br.peek(MAX_BITS as u32);
        self.decode_buffered(br)
    }

    /// The fast loop has at least 48 bits before starting a symbol, enough
    /// for the entire literal/length + extras + distance + extras sequence.
    #[inline]
    fn decode_buffered(&self, br: &mut Bits) -> Result<u32, InflateError> {
        let mut e = self.table[br.buf as usize & ((1 << FAST_BITS) - 1)];
        if e & 15 == 0 && e & (3 << 8) == LINK {
            let width = (e >> 4) & 15;
            let suffix = (br.buf >> FAST_BITS) as usize & ((1 << width) - 1);
            e = self.table[(e >> 16) as usize + suffix];
        }
        let len = e & 15;
        if len == 0 {
            return Err(InflateError::BadCode);
        }
        if e & (3 << 8) == LINK {
            return Err(InflateError::BadSymbol);
        }
        br.consume(len);
        Ok(e)
    }
}

fn fixed_tables() -> &'static (Huffman, Huffman) {
    static FIXED: OnceLock<(Huffman, Huffman)> = OnceLock::new();
    FIXED.get_or_init(|| {
        let mut lit = [0u8; 288];
        lit[..144].fill(8);
        lit[144..256].fill(9);
        lit[256..280].fill(7);
        lit[280..].fill(8);
        (
            Huffman::new(&lit, Alphabet::Litlen).unwrap(),
            Huffman::new(&[5u8; 32], Alphabet::Distance).unwrap(),
        )
    })
}

fn dynamic_tables(br: &mut Bits) -> Result<(Huffman, Huffman), InflateError> {
    let hlit = br.bits(5) as usize + 257;
    let hdist = br.bits(5) as usize + 1;
    // RFC 1951 §3.2.7: at most 286 literal/length and 30 distance codes.
    if hlit > 286 || hdist > 30 {
        return Err(InflateError::BadCode);
    }
    let hclen = br.bits(4) as usize + 4;
    let mut clen = [0u8; 19];
    for &i in &CLEN_ORDER[..hclen] {
        clen[i] = br.bits(3) as u8;
    }
    let clen_code = Huffman::new(&clen, Alphabet::Symbols)?;
    let mut lengths = [0u8; 286 + 30];
    let mut i = 0;
    while i < hlit + hdist {
        let sym = clen_code.decode(br)? >> 16;
        br.consumed()?;
        let (value, repeat) = match sym {
            0..=15 => (sym as u8, 1),
            16 => {
                let prev = *lengths
                    .get(i.wrapping_sub(1))
                    .filter(|_| i > 0)
                    .ok_or(InflateError::BadSymbol)?;
                (prev, 3 + br.bits(2) as usize)
            }
            17 => (0, 3 + br.bits(3) as usize),
            18 => (0, 11 + br.bits(7) as usize),
            _ => return Err(InflateError::BadSymbol),
        };
        if i + repeat > hlit + hdist {
            return Err(InflateError::BadSymbol);
        }
        lengths[i..i + repeat].fill(value);
        i += repeat;
    }
    if lengths[256] == 0 {
        return Err(InflateError::BadCode);
    }
    Ok((
        Huffman::new(&lengths[..hlit], Alphabet::Litlen)?,
        Huffman::new(&lengths[hlit..hlit + hdist], Alphabet::Distance)?,
    ))
}

/// Output bounds for one stream: back-references may not reach before `base`,
/// and the stream may not grow `out` past `end`. When `out` reaches
/// `flush_at`, everything but the last [`WSIZE`] bytes goes to the sink.
struct Window<'s> {
    base: usize,
    end: usize,
    flush_at: usize,
    sink: Option<ByteSink<'s>>,
    sink_error: Option<io::Error>,
}

/// Receives decoded bytes in order; an error stops decoding.
pub type ByteSink<'a> = &'a mut dyn FnMut(&[u8]) -> io::Result<()>;

/// DEFLATE's maximum back-reference distance.
const WSIZE: usize = 32768;

impl Window<'_> {
    #[inline]
    fn maybe_flush(&mut self, out: &mut Vec<u8>) -> Result<(), InflateError> {
        if out.len() >= self.flush_at {
            self.flush(out, WSIZE)
        } else {
            Ok(())
        }
    }

    /// Hand all but the last `keep` bytes of `out` to the sink.
    fn flush(&mut self, out: &mut Vec<u8>, keep: usize) -> Result<(), InflateError> {
        let Some(sink) = self.sink.as_mut() else {
            return Ok(());
        };
        let n = out.len().saturating_sub(keep);
        if n == 0 {
            return Ok(());
        }
        if let Err(e) = sink(&out[..n]) {
            self.sink_error = Some(e);
            return Err(InflateError::Sink);
        }
        out.drain(..n);
        self.base = self.base.saturating_sub(n);
        self.end = self.end.saturating_sub(n);
        Ok(())
    }
}

/// Copy a match in logarithmically many bulk copies, including distances 1..7.
/// Each copy reads only bytes already produced, so memcpy overlap is impossible.
#[inline]
fn copy_match(out: &mut [u8], at: usize, distance: usize, len: usize) {
    if distance == 1 {
        let byte = out[at - 1];
        out[at..at + len].fill(byte);
        return;
    }
    let mut copied = 0;
    let mut available = distance;
    while copied < len {
        let n = available.min(len - copied);
        out.copy_within(at - distance..at - distance + n, at + copied);
        copied += n;
        available += n;
    }
}

#[inline]
fn match_distance(br: &mut Bits, dist: &Huffman) -> Result<usize, InflateError> {
    let e = dist.decode(br)?;
    Ok((e >> 16) as usize + br.bits((e >> 4) & 15) as usize)
}

fn inflate_block(
    br: &mut Bits,
    lit: &Huffman,
    dist: &Huffman,
    out: &mut Vec<u8>,
    w: &mut Window,
) -> Result<(), InflateError> {
    loop {
        w.maybe_flush(out)?;
        // Initialize one output slab, then fill by index: no per-byte Vec
        // capacity checks. Stop before the largest possible match (258).
        // Eight real input bytes suffice for a complete length/distance pair
        // (15+5+15+13 bits), even with an empty bit buffer.
        let start = out.len();
        let room = (w.end - start).min(w.flush_at - start).min(16 * 1024);
        if room >= 258 && br.data.len().saturating_sub(br.pos) >= 8 {
            out.resize(start + room, 0);
            let mut at = start;
            let result = (|| {
                while start + room - at >= 258 && br.data.len().saturating_sub(br.pos) >= 8 {
                    if br.nbits < 48 {
                        br.refill();
                    }
                    let e = lit.decode_buffered(br)?;
                    match e & (3 << 8) {
                        0 => {
                            out[at] = (e >> 16) as u8;
                            at += 1;
                        }
                        END => return Ok(true),
                        LENGTH => {
                            let len = (e >> 16) as usize + br.buffered_bits((e >> 4) & 15) as usize;
                            let de = dist.decode_buffered(br)?;
                            let d = (de >> 16) as usize + br.buffered_bits((de >> 4) & 15) as usize;
                            if d > at - w.base {
                                return Err(InflateError::DistanceTooFar);
                            }
                            copy_match(out, at, d, len);
                            at += len;
                        }
                        _ => return Err(InflateError::BadSymbol),
                    }
                }
                Ok(false)
            })();
            out.truncate(at);
            if result? {
                return Ok(());
            }
        }
        // Checked tail also handles tiny output limits and exact flush boundaries.
        w.maybe_flush(out)?;
        let e = lit.decode(br)?;
        br.consumed()?;
        match e & (3 << 8) {
            0 => {
                if out.len() >= w.end {
                    return Err(InflateError::OutputTooLarge);
                }
                out.push((e >> 16) as u8);
            }
            END => return Ok(()),
            LENGTH => {
                let len = (e >> 16) as usize + br.bits((e >> 4) & 15) as usize;
                let d = match_distance(br, dist)?;
                br.consumed()?;
                if d > out.len() - w.base {
                    return Err(InflateError::DistanceTooFar);
                }
                if len > w.end - out.len() {
                    return Err(InflateError::OutputTooLarge);
                }
                let at = out.len();
                out.resize(at + len, 0);
                copy_match(out, at, d, len);
            }
            _ => return Err(InflateError::BadSymbol),
        }
    }
}

/// Decode one raw DEFLATE stream from `input`, appending at most `limit`
/// bytes to `out`. Returns the number of input bytes consumed.
pub fn inflate(input: &[u8], out: &mut Vec<u8>, limit: usize) -> Result<usize, InflateError> {
    let mut w = Window {
        base: out.len(),
        end: out.len().saturating_add(limit),
        flush_at: usize::MAX,
        sink: None,
        sink_error: None,
    };
    inflate_stream(input, out, &mut w)
}

/// Decode one raw DEFLATE stream from `input`, passing the output to `sink`
/// in chunks (with bounded memory). Returns the number of input bytes consumed.
pub fn inflate_to(
    input: &[u8],
    sink: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<usize> {
    const CHUNK: usize = 1 << 20;
    let mut out = Vec::with_capacity(CHUNK + WSIZE + 65536);
    let mut w = Window {
        base: 0,
        end: usize::MAX,
        flush_at: CHUNK + WSIZE,
        sink: Some(sink),
        sink_error: None,
    };
    let result =
        inflate_stream(input, &mut out, &mut w).and_then(|n| w.flush(&mut out, 0).map(|_| n));
    match result {
        Ok(n) => Ok(n),
        Err(InflateError::Sink) => Err(w.sink_error.take().expect("sink error recorded")),
        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn inflate_stream(input: &[u8], out: &mut Vec<u8>, w: &mut Window) -> Result<usize, InflateError> {
    let mut br = Bits::new(input);
    loop {
        let last = br.bits(1) == 1;
        let kind = br.bits(2);
        br.consumed()?;
        match kind {
            0 => {
                br.align_to_byte();
                let len = br.bits(16) as usize;
                let nlen = br.bits(16) as usize;
                br.consumed()?;
                if len != !nlen & 0xffff {
                    return Err(InflateError::BadStoredLength);
                }
                if len > w.end - out.len() {
                    return Err(InflateError::OutputTooLarge);
                }
                // The bit buffer is byte-aligned; drain it, then copy directly.
                let mut remaining = len;
                while remaining > 0 && br.nbits >= 8 {
                    out.push(br.bits(8) as u8);
                    remaining -= 1;
                }
                let start = br.consumed()?;
                let end = start
                    .checked_add(remaining)
                    .filter(|&e| e <= input.len())
                    .ok_or(InflateError::Truncated)?;
                out.extend_from_slice(&input[start..end]);
                w.maybe_flush(out)?;
                br = Bits {
                    data: input,
                    pos: end,
                    buf: 0,
                    nbits: 0,
                };
            }
            1 => {
                let (lit, dist) = fixed_tables();
                inflate_block(&mut br, lit, dist, out, w)?;
            }
            2 => {
                let (lit, dist) = dynamic_tables(&mut br)?;
                inflate_block(&mut br, &lit, &dist, out, w)?;
            }
            _ => return Err(InflateError::BadBlockType),
        }
        if last {
            return br.consumed();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Raw deflate streams produced by zlib (Python `zlib.compressobj(L, 8, -15)`).
    #[test]
    fn stored_block() {
        // zlib level 0 of b"hello"
        let s = [0x01, 0x05, 0x00, 0xfa, 0xff, b'h', b'e', b'l', b'l', b'o'];
        let mut out = Vec::new();
        assert_eq!(inflate(&s, &mut out, usize::MAX), Ok(s.len()));
        assert_eq!(out, b"hello");
    }

    #[test]
    fn fixed_huffman_with_overlapping_match() {
        // zlib level 6 of b"abcabcabcabcabcabc"
        let s = [0x4b, 0x4c, 0x4a, 0x4e, 0x44, 0x45, 0x00];
        let mut out = Vec::new();
        assert_eq!(inflate(&s, &mut out, usize::MAX), Ok(s.len()));
        assert_eq!(out, b"abcabcabcabcabcabc");
    }

    fn fixture_input() -> Vec<u8> {
        let mut d: Vec<u8> = (0..3000u32)
            .map(|i| ((i * 7 + i / 13) % 251) as u8)
            .collect();
        for _ in 0..500 {
            d.extend_from_slice(b"ACGT");
        }
        d
    }

    #[test]
    fn dynamic_huffman_fixtures() {
        for f in ["dynamic_l6.deflate", "dynamic_l2.deflate"] {
            let path = format!("{}/tests/fixtures/{f}", env!("CARGO_MANIFEST_DIR"));
            let s = std::fs::read(path).unwrap();
            let mut out = Vec::new();
            assert_eq!(inflate(&s, &mut out, usize::MAX), Ok(s.len()), "{f}");
            assert_eq!(out, fixture_input(), "{f}");
        }
    }

    #[test]
    fn empty_stream() {
        // zlib level 6 of b""
        let s = [0x03, 0x00];
        let mut out = Vec::new();
        assert_eq!(inflate(&s, &mut out, usize::MAX), Ok(2));
        assert!(out.is_empty());
    }

    #[test]
    fn output_limit_is_enforced() {
        let s = [0x4b, 0x4c, 0x4a, 0x4e, 0x44, 0x45, 0x00];
        assert_eq!(
            inflate(&s, &mut Vec::new(), 17),
            Err(InflateError::OutputTooLarge)
        );
        assert_eq!(inflate(&s, &mut Vec::new(), 18), Ok(s.len()));
        let stored = [0x01, 0x05, 0x00, 0xfa, 0xff, b'h', b'e', b'l', b'l', b'o'];
        assert_eq!(
            inflate(&stored, &mut Vec::new(), 4),
            Err(InflateError::OutputTooLarge)
        );
    }

    #[test]
    fn back_reference_cannot_reach_previous_output() {
        // Same stream, but `out` already holds bytes from an earlier stream:
        // the first match (distance 3) is still valid only within this stream.
        let s = [0x4b, 0x4c, 0x4a, 0x4e, 0x44, 0x45, 0x00];
        let mut out = b"previous block".to_vec();
        assert_eq!(inflate(&s, &mut out, usize::MAX), Ok(s.len()));
        assert_eq!(&out[14..], b"abcabcabcabcabcabc");
    }

    #[test]
    fn oversized_code_counts_are_rejected() {
        // BFINAL=1, BTYPE=2 (dynamic), HLIT=31 (288 codes), HDIST=31 (32 codes).
        let s = [0b1111_1101, 0b1111_1111, 0xff, 0xff];
        assert_eq!(
            inflate(&s, &mut Vec::new(), usize::MAX),
            Err(InflateError::BadCode)
        );
    }

    #[test]
    fn streaming_matches_one_shot_across_many_flushes() {
        // ~5 MiB of output with long-distance matches crosses several 1 MiB flushes.
        let mut data = Vec::new();
        for i in 0u32..2_500_000 {
            data.extend_from_slice(&(i.wrapping_mul(2654435761) >> 7).to_le_bytes()[..2]);
            if i % 997 == 0 {
                let from = data.len().saturating_sub(30_000);
                let chunk: Vec<u8> = data[from..from + 250.min(data.len() - from)].to_vec();
                data.extend_from_slice(&chunk);
            }
        }
        let mut z = crate::deflate::ZlibFast::new(2);
        let mut compressed = Vec::new();
        z.deflate_finish(&data, &mut compressed, usize::MAX);
        let mut streamed = Vec::new();
        let mut chunks = 0;
        let used = inflate_to(&compressed, &mut |c: &[u8]| {
            chunks += 1;
            streamed.extend_from_slice(c);
            Ok(())
        })
        .unwrap();
        assert_eq!(used, compressed.len());
        assert!(chunks > 2, "expected several flushes, got {chunks}");
        assert!(streamed == data);
    }

    #[test]
    fn sink_errors_propagate() {
        let s = [0x4b, 0x4c, 0x4a, 0x4e, 0x44, 0x45, 0x00];
        let err = inflate_to(&s, &mut |_| Err(io::Error::other("disk full"))).unwrap_err();
        assert_eq!(err.to_string(), "disk full");
    }

    #[test]
    fn truncated_input_is_an_error() {
        let s = [0x4b, 0x4c, 0x4a];
        assert_eq!(
            inflate(&s, &mut Vec::new(), usize::MAX),
            Err(InflateError::Truncated)
        );
    }
}

#[cfg(test)]
mod oracle;

#[cfg(test)]
mod checks;
