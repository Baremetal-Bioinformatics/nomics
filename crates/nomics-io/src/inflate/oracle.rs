//! Raw DEFLATE decoder (RFC 1951).
//!
//! Huffman decoding uses a 2^10-entry lookup table for short codes and a
//! canonical (counts/symbols) walk for longer ones.

// Frozen pre-optimization implementation, used only for differential tests and
// throughput comparisons. Keep its decoding behavior independent of the new code.
#![allow(dead_code)]

use std::fmt;
use std::io;

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

    /// Top up the bit buffer to at least 57 bits. Past the end of input, zero
    /// bytes are fed; `consumed()` detects reads beyond the real input.
    #[inline]
    fn refill(&mut self) {
        while self.nbits <= 56 {
            let b = self.data.get(self.pos).copied().unwrap_or(0);
            self.buf |= (b as u64) << self.nbits;
            self.nbits += 8;
            self.pos += 1;
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
        let v = self.peek(n);
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

struct Huffman {
    /// `(symbol << 4) | length` for codes of at most FAST_BITS bits; 0 = longer code.
    fast: Vec<u16>,
    counts: [u16; MAX_BITS + 1],
    /// Symbols ordered by (code length, symbol value).
    symbols: Vec<u16>,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Result<Self, InflateError> {
        let mut counts = [0u16; MAX_BITS + 1];
        for &l in lengths {
            counts[l as usize] += 1;
        }
        counts[0] = 0;
        // Reject over-subscribed codes; incomplete codes are legal (e.g. a
        // single distance code) and simply never match the unused prefixes.
        let mut left: i32 = 1;
        for &c in &counts[1..] {
            left = (left << 1) - c as i32;
            if left < 0 {
                return Err(InflateError::BadCode);
            }
        }
        let mut offs = [0u16; MAX_BITS + 2];
        for l in 1..=MAX_BITS {
            offs[l + 1] = offs[l] + counts[l];
        }
        let mut symbols = vec![0u16; offs[MAX_BITS + 1] as usize];
        for (sym, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbols[offs[l as usize] as usize] = sym as u16;
                offs[l as usize] += 1;
            }
        }
        let mut fast = vec![0u16; 1 << FAST_BITS];
        let mut code: u32 = 0;
        let mut idx = 0usize;
        for len in 1..=MAX_BITS as u32 {
            for _ in 0..counts[len as usize] {
                if len <= FAST_BITS {
                    let rev = code.reverse_bits() >> (32 - len);
                    let entry = (symbols[idx] << 4) | len as u16;
                    let mut i = rev as usize;
                    while i < fast.len() {
                        fast[i] = entry;
                        i += 1 << len;
                    }
                }
                code += 1;
                idx += 1;
            }
            code <<= 1;
        }
        Ok(Huffman {
            fast,
            counts,
            symbols,
        })
    }

    #[inline]
    fn decode(&self, br: &mut Bits) -> Result<u16, InflateError> {
        let e = self.fast[br.peek(FAST_BITS) as usize];
        if e != 0 {
            br.consume((e & 0xf) as u32);
            return Ok(e >> 4);
        }
        // Canonical walk, one bit at a time (codes are MSB-first).
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for len in 1..=MAX_BITS {
            code |= br.bits(1) as i32;
            let count = self.counts[len] as i32;
            if code - first < count {
                return Ok(self.symbols[(index + code - first) as usize]);
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        Err(InflateError::BadCode)
    }
}

fn fixed_tables() -> (Huffman, Huffman) {
    let mut lit = [0u8; 288];
    lit[..144].fill(8);
    lit[144..256].fill(9);
    lit[256..280].fill(7);
    lit[280..].fill(8);
    (
        Huffman::new(&lit).unwrap(),
        Huffman::new(&[5u8; 30]).unwrap(),
    )
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
    let clen_code = Huffman::new(&clen)?;
    let mut lengths = [0u8; 286 + 30];
    let mut i = 0;
    while i < hlit + hdist {
        let sym = clen_code.decode(br)?;
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
        Huffman::new(&lengths[..hlit])?,
        Huffman::new(&lengths[hlit..hlit + hdist])?,
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

fn inflate_block(
    br: &mut Bits,
    lit: &Huffman,
    dist: &Huffman,
    out: &mut Vec<u8>,
    w: &mut Window,
) -> Result<(), InflateError> {
    loop {
        w.maybe_flush(out)?;
        let sym = lit.decode(br)?;
        match sym {
            0..=255 => {
                if out.len() >= w.end {
                    return Err(InflateError::OutputTooLarge);
                }
                out.push(sym as u8)
            }
            256 => return Ok(()),
            257..=285 => {
                let i = (sym - 257) as usize;
                let len = LEN_BASE[i] as usize + br.bits(LEN_EXTRA[i] as u32) as usize;
                let dsym = dist.decode(br)? as usize;
                if dsym >= 30 {
                    return Err(InflateError::BadSymbol);
                }
                let d = DIST_BASE[dsym] as usize + br.bits(DIST_EXTRA[dsym] as u32) as usize;
                if d > out.len() - w.base {
                    return Err(InflateError::DistanceTooFar);
                }
                if len > w.end - out.len() {
                    return Err(InflateError::OutputTooLarge);
                }
                let start = out.len() - d;
                if d >= len {
                    out.extend_from_within(start..start + len);
                } else {
                    // Overlapping copy: the run repeats with period `d`.
                    for k in 0..len {
                        let b = out[start + k];
                        out.push(b);
                    }
                }
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
    let mut fixed: Option<(Huffman, Huffman)> = None;
    loop {
        let last = br.bits(1) == 1;
        match br.bits(2) {
            0 => {
                br.align_to_byte();
                let len = br.bits(16) as usize;
                let nlen = br.bits(16) as usize;
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
                let (lit, dist) = fixed.get_or_insert_with(fixed_tables);
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
