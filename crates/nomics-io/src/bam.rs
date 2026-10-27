//! BAM header and record codec (SAMv1 §4.2).
//!
//! Records have two forms:
//! - raw bytes (`block_size` excluded), as read from or written to the stream;
//! - [`Record`], fully decoded and owned. [`Record::encode`] recomputes the
//!   derived fields (`bin`, `l_read_name`, `n_cigar_op`, `l_seq`, packed
//!   sequence) rather than copying them, so a decode→encode round trip against
//!   a golden checks our understanding of the format, not just a memcpy.
//!
//! Tag values keep their on-disk type code, so re-encoding never changes a
//! type. Choosing types for *new* tags is writer policy and belongs to the
//! stage that creates them (htsjdk picks the smallest fitting integer type).

use std::io::{self, Read, Write};

pub const MAGIC: &[u8; 4] = b"BAM\x01";

/// Upper bound on any single length field we will allocate for (header text,
/// one record). Far above real data; stops a corrupt length from allocating GiBs.
pub const MAX_ALLOC: usize = 256 << 20;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    pub name: Vec<u8>,
    pub len: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// SAM header text exactly as stored (`l_text` bytes, including any NUL padding).
    pub text: Vec<u8>,
    pub refs: Vec<Reference>,
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn read_i32<R: Read>(r: &mut R) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

fn read_len<R: Read>(r: &mut R, what: &str) -> io::Result<usize> {
    let v = read_i32(r)?;
    usize::try_from(v)
        .ok()
        .filter(|&n| n <= MAX_ALLOC)
        .ok_or_else(|| invalid(format!("implausible {what}: {v}")))
}

impl Header {
    pub fn read<R: Read>(r: &mut R) -> io::Result<Header> {
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(invalid("not a BAM stream (bad magic)"));
        }
        let mut text = vec![0u8; read_len(r, "l_text")?];
        r.read_exact(&mut text)?;
        let n_ref = read_len(r, "n_ref")?;
        let mut refs = Vec::with_capacity(n_ref.min(1 << 16));
        for _ in 0..n_ref {
            let mut name = vec![0u8; read_len(r, "l_name")?];
            r.read_exact(&mut name)?;
            if name.pop() != Some(0) {
                return Err(invalid("reference name not NUL-terminated"));
            }
            let len = read_i32(r)?;
            let len = u32::try_from(len).map_err(|_| invalid(format!("negative l_ref: {len}")))?;
            refs.push(Reference { name, len });
        }
        Ok(Header { text, refs })
    }

    pub fn write<W: Write>(&self, w: &mut W) -> io::Result<()> {
        w.write_all(MAGIC)?;
        w.write_all(&(self.text.len() as i32).to_le_bytes())?;
        w.write_all(&self.text)?;
        w.write_all(&(self.refs.len() as i32).to_le_bytes())?;
        for r in &self.refs {
            w.write_all(&(r.name.len() as i32 + 1).to_le_bytes())?;
            w.write_all(&r.name)?;
            w.write_all(&[0])?;
            w.write_all(&(r.len as i32).to_le_bytes())?;
        }
        Ok(())
    }
}

/// Read the next raw record into `buf` (without its `block_size` prefix).
/// Returns `Ok(false)` at a clean end of stream.
pub fn read_raw_record<R: Read>(r: &mut R, buf: &mut Vec<u8>) -> io::Result<bool> {
    let mut len = [0u8; 4];
    let mut got = 0;
    while got < 4 {
        match r.read(&mut len[got..])? {
            0 if got == 0 => return Ok(false),
            0 => return Err(io::ErrorKind::UnexpectedEof.into()),
            n => got += n,
        }
    }
    let n = i32::from_le_bytes(len);
    let n = usize::try_from(n)
        .ok()
        .filter(|&n| (32..=MAX_ALLOC).contains(&n))
        .ok_or_else(|| invalid(format!("bad block_size {n}")))?;
    buf.clear();
    buf.resize(n, 0);
    r.read_exact(buf)?;
    Ok(true)
}

pub fn write_raw_record<W: Write>(w: &mut W, raw: &[u8]) -> io::Result<()> {
    w.write_all(&(raw.len() as i32).to_le_bytes())?;
    w.write_all(raw)
}

/// Zero-copy view of the fixed fields, name and CIGAR of a raw record (as
/// returned by [`read_raw_record`]), for scans that do not need full decoding.
#[derive(Clone, Copy)]
pub struct RecordView<'a>(&'a [u8]);

impl<'a> RecordView<'a> {
    /// Checks only that the fixed fields, name and CIGAR are in bounds.
    pub fn new(raw: &'a [u8]) -> io::Result<Self> {
        if raw.len() < 32 {
            return Err(invalid("record shorter than its fixed fields"));
        }
        let v = RecordView(raw);
        if raw.len() < 32 + v.l_read_name() + 4 * v.n_cigar_op() {
            return Err(invalid("record truncated before end of CIGAR"));
        }
        Ok(v)
    }

    fn i32_at(&self, at: usize) -> i32 {
        i32::from_le_bytes(self.0[at..at + 4].try_into().unwrap())
    }

    fn l_read_name(&self) -> usize {
        self.0[8] as usize
    }

    fn n_cigar_op(&self) -> usize {
        u16::from_le_bytes([self.0[12], self.0[13]]) as usize
    }

    pub fn ref_id(&self) -> i32 {
        self.i32_at(0)
    }

    pub fn pos(&self) -> i32 {
        self.i32_at(4)
    }

    pub fn flag(&self) -> u16 {
        u16::from_le_bytes([self.0[14], self.0[15]])
    }

    /// Read name without the trailing NUL.
    pub fn read_name(&self) -> &'a [u8] {
        let n = self.l_read_name();
        &self.0[32..32 + n.saturating_sub(1)]
    }

    pub fn cigar(&self) -> impl Iterator<Item = u32> + 'a {
        let start = 32 + self.l_read_name();
        let raw = self.0;
        (0..self.n_cigar_op()).map(move |i| {
            u32::from_le_bytes(raw[start + 4 * i..start + 4 * i + 4].try_into().unwrap())
        })
    }

    /// Reference bases spanned (ops M, D, N, =, X).
    pub fn reference_length(&self) -> u32 {
        self.cigar()
            .filter(|&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
            .map(|c| c >> 4)
            .sum()
    }

    pub fn mate_ref_id(&self) -> i32 {
        self.i32_at(20)
    }

    pub fn mapq(&self) -> u8 {
        self.0[9]
    }

    /// 0-based mate position (-1 if none).
    pub fn mate_pos(&self) -> i32 {
        self.i32_at(24)
    }

    /// Signed observed template length.
    pub fn tlen(&self) -> i32 {
        self.i32_at(28)
    }

    /// Read length (`l_seq`).
    pub fn len(&self) -> usize {
        self.l_seq()
    }

    pub fn is_empty(&self) -> bool {
        self.l_seq() == 0
    }

    /// Base `i` as an uppercase letter from `=ACMGRSVTWYHKDBN`.
    #[inline]
    pub fn base(&self, i: usize) -> u8 {
        let b = self.0[self.seq_start() + i / 2];
        let code = if i % 2 == 0 { b >> 4 } else { b & 0xf };
        SEQ_CODES[code as usize]
    }

    fn l_seq(&self) -> usize {
        u32::from_le_bytes(self.0[16..20].try_into().unwrap()) as usize
    }

    fn seq_start(&self) -> usize {
        32 + self.l_read_name() + 4 * self.n_cigar_op()
    }

    /// Phred qualities (no +33), `None` if the record is truncated.
    pub fn quals(&self) -> Option<&'a [u8]> {
        let start = self.seq_start() + self.l_seq().div_ceil(2);
        self.0.get(start..start + self.l_seq())
    }

    /// Bytes of the auxiliary tag block, `None` if the record is truncated.
    pub fn aux(&self) -> Option<&'a [u8]> {
        self.0.get(self.aux_offset()..)
    }

    /// Offset of the auxiliary tag block within the record payload.
    pub fn aux_offset(&self) -> usize {
        self.seq_start() + self.l_seq().div_ceil(2) + self.l_seq()
    }

    /// Value of a `Z` (string) tag without its NUL, scanning the raw tags.
    /// `None` if absent, of another type, or if the tag block is malformed.
    pub fn string_tag(&self, tag: &[u8; 2]) -> Option<&'a [u8]> {
        let mut aux = self.aux()?;
        while aux.len() >= 3 {
            let (t, ty) = ([aux[0], aux[1]], aux[2]);
            let body = &aux[3..];
            let len = aux_value_len(ty, body)?;
            if t == *tag {
                return if ty == b'Z' {
                    Some(&body[..len - 1])
                } else {
                    None
                };
            }
            aux = body.get(len..)?;
        }
        None
    }
}

/// Length of a tag value of type `ty` at the start of `body` (the bytes after
/// the tag and type), `None` for an unknown type or a truncated value.
pub fn aux_value_len(ty: u8, body: &[u8]) -> Option<usize> {
    Some(match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => body.iter().position(|&c| c == 0)? + 1,
        b'B' => {
            let width = match *body.first()? {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            let n = u32::from_le_bytes(body.get(1..5)?.try_into().unwrap()) as usize;
            5 + n * width
        }
        _ => return None,
    })
}

/// A whole BAM held in memory as its uncompressed stream (header followed by
/// `block_size`-prefixed records) plus the offset of every record's payload.
/// Stages hand these to each other instead of writing and re-reading BAMs;
/// writing one is just compressing `data`.
pub struct RecordSet {
    pub header: Header,
    data: Vec<u8>,
    /// Payload start of record `i` is `offsets[i]`; `offsets[len]` is one past the end + 4.
    offsets: Vec<usize>,
}

impl RecordSet {
    /// Index an uncompressed BAM stream.
    pub fn from_stream(data: Vec<u8>) -> io::Result<RecordSet> {
        let mut cursor = &data[..];
        let header = Header::read(&mut cursor)?;
        let mut at = data.len() - cursor.len();
        let mut offsets = Vec::new();
        while at + 4 <= data.len() {
            let n = i32::from_le_bytes(data[at..at + 4].try_into().unwrap());
            let n = usize::try_from(n)
                .ok()
                .filter(|&n| n >= 32)
                .ok_or_else(|| invalid(format!("bad block_size at {at}")))?;
            if at + 4 + n > data.len() {
                return Err(invalid("truncated BAM record"));
            }
            offsets.push(at + 4);
            at += 4 + n;
        }
        if at != data.len() {
            return Err(invalid("trailing bytes after last BAM record"));
        }
        offsets.push(at + 4);
        Ok(RecordSet {
            header,
            data,
            offsets,
        })
    }

    /// Build from a header and already-encoded, `block_size`-prefixed records.
    pub fn from_parts(header: Header, records: &[u8]) -> io::Result<RecordSet> {
        let mut data = Vec::new();
        header.write(&mut data)?;
        data.extend_from_slice(records);
        RecordSet::from_stream(data)
    }

    /// Read and decode a BGZF-compressed BAM using `threads` threads.
    pub fn load(path: &std::path::Path, threads: usize) -> io::Result<RecordSet> {
        let compressed = std::fs::read(path)?;
        let mut data = Vec::new();
        crate::bgzf::decode_parallel(&compressed, threads, &mut |b| {
            data.extend_from_slice(b);
            Ok(())
        })?;
        drop(compressed);
        RecordSet::from_stream(data)
    }

    /// Write as BGZF with the htsjdk block profile on `threads` threads,
    /// compressed as [`crate::bgzf::BamCompression::from_env`] selects.
    /// Write the set as a BAM, or as a light spill file when `path` has the
    /// spill extension ([`crate::spill::is_spill_path`]).
    pub fn write(&self, path: &std::path::Path, threads: usize) -> io::Result<()> {
        if crate::spill::is_spill_path(path) {
            return crate::spill::write_file(path, &self.data, threads);
        }
        let f = std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(path)?);
        let mut w = crate::bgzf::ParallelBgzfWriter::for_bam(f, threads);
        w.write_all(&self.data)?;
        w.finish()?.flush()
    }

    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Raw record `i` (without `block_size`).
    pub fn record(&self, i: usize) -> &[u8] {
        &self.data[self.offsets[i]..self.offsets[i + 1] - 4]
    }

    /// Bytes held (the uncompressed stream).
    pub fn bytes(&self) -> usize {
        self.data.len()
    }
}

/// CIGAR operation codes, in BAM numeric order.
pub const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";

/// Reference bases consumed by a CIGAR (ops M, D, N, =, X).
pub fn reference_length(cigar: &[u32]) -> u32 {
    cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|&c| c >> 4)
        .sum()
}

/// UCSC binning scheme for the half-open interval [beg, end) (SAMv1 §5.3).
pub fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    for (shift, offset) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return (offset + (beg >> shift)) as u16;
        }
    }
    0
}

/// How the `bin` field is filled. Writers disagree for unplaced reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinRule {
    /// SAMv1 §5.3 / htslib: always reg2bin; unplaced reads get 4680.
    Sam,
    /// htsjdk `BAMRecordCodec.encode`: 0 when the read has no alignment start
    /// (pos -1). Verified on tier1 `md.bam`: all 8,830 unplaced records have bin 0.
    /// htsjdk also writes 0 when the reference is too long for the bin field
    /// (> 2^29), which cannot occur on GRCh38 and is not implemented.
    Htsjdk,
}

const SEQ_CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

fn base_code(b: u8) -> Option<u8> {
    SEQ_CODES
        .iter()
        .position(|&c| c == b.to_ascii_uppercase())
        .map(|i| i as u8)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntType {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
}

impl IntType {
    fn code(self) -> u8 {
        match self {
            IntType::I8 => b'c',
            IntType::U8 => b'C',
            IntType::I16 => b's',
            IntType::U16 => b'S',
            IntType::I32 => b'i',
            IntType::U32 => b'I',
        }
    }

    fn from_code(c: u8) -> Option<IntType> {
        Some(match c {
            b'c' => IntType::I8,
            b'C' => IntType::U8,
            b's' => IntType::I16,
            b'S' => IntType::U16,
            b'i' => IntType::I32,
            b'I' => IntType::U32,
            _ => return None,
        })
    }

    fn width(self) -> usize {
        match self {
            IntType::I8 | IntType::U8 => 1,
            IntType::I16 | IntType::U16 => 2,
            IntType::I32 | IntType::U32 => 4,
        }
    }

    fn decode(self, b: &[u8]) -> i64 {
        match self {
            IntType::I8 => b[0] as i8 as i64,
            IntType::U8 => b[0] as i64,
            IntType::I16 => i16::from_le_bytes([b[0], b[1]]) as i64,
            IntType::U16 => u16::from_le_bytes([b[0], b[1]]) as i64,
            IntType::I32 => i32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i64,
            IntType::U32 => u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i64,
        }
    }

    fn encode(self, v: i64, out: &mut Vec<u8>) {
        match self {
            IntType::I8 | IntType::U8 => out.push(v as u8),
            IntType::I16 | IntType::U16 => out.extend_from_slice(&(v as u16).to_le_bytes()),
            IntType::I32 | IntType::U32 => out.extend_from_slice(&(v as u32).to_le_bytes()),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ArrayValue {
    Int(IntType, Vec<i64>),
    /// IEEE-754 bit patterns, so NaN payloads round-trip exactly.
    Float(Vec<u32>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Char(u8),
    Int(IntType, i64),
    /// IEEE-754 bit pattern.
    Float(u32),
    String(Vec<u8>),
    Hex(Vec<u8>),
    Array(ArrayValue),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub ref_id: i32,
    /// 0-based leftmost position, -1 if none.
    pub pos: i32,
    pub mapq: u8,
    pub flag: u16,
    pub next_ref_id: i32,
    pub next_pos: i32,
    pub tlen: i32,
    /// Without the trailing NUL.
    pub read_name: Vec<u8>,
    /// `len << 4 | op`, as in BAM.
    pub cigar: Vec<u32>,
    /// Bases as uppercase IUPAC letters from `=ACMGRSVTWYHKDBN`.
    pub seq: Vec<u8>,
    /// Phred scores (no +33). All 0xff when absent.
    pub qual: Vec<u8>,
    pub tags: Vec<([u8; 2], Value)>,
}

struct Cursor<'a> {
    b: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> io::Result<&'a [u8]> {
        let s = self
            .b
            .get(self.at..self.at + n)
            .ok_or_else(|| invalid("record truncated"))?;
        self.at += n;
        Ok(s)
    }
    fn u8(&mut self) -> io::Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> io::Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> io::Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn cstr(&mut self) -> io::Result<&'a [u8]> {
        let rest = &self.b[self.at..];
        let n = rest
            .iter()
            .position(|&c| c == 0)
            .ok_or_else(|| invalid("unterminated string tag"))?;
        self.at += n + 1;
        Ok(&rest[..n])
    }
}

impl Record {
    /// Fully decode a raw record (as returned by [`read_raw_record`]).
    pub fn decode(raw: &[u8]) -> io::Result<Record> {
        let mut c = Cursor { b: raw, at: 0 };
        let ref_id = c.i32()?;
        let pos = c.i32()?;
        let l_read_name = c.u8()? as usize;
        let mapq = c.u8()?;
        let _bin = c.u16()?;
        let n_cigar = c.u16()? as usize;
        let flag = c.u16()?;
        let l_seq = c.u32()? as usize;
        let next_ref_id = c.i32()?;
        let next_pos = c.i32()?;
        let tlen = c.i32()?;
        let name = c.take(l_read_name)?;
        let read_name = name
            .strip_suffix(&[0])
            .ok_or_else(|| invalid("read name not NUL-terminated"))?
            .to_vec();
        let cigar = (0..n_cigar)
            .map(|_| c.u32())
            .collect::<io::Result<Vec<_>>>()?;
        let packed = c.take(l_seq.div_ceil(2))?;
        let seq = (0..l_seq)
            .map(|i| SEQ_CODES[((packed[i / 2] >> (4 * (1 - i % 2))) & 0xf) as usize])
            .collect();
        let qual = c.take(l_seq)?.to_vec();
        let mut tags = Vec::new();
        while c.at < raw.len() {
            let tag: [u8; 2] = c.take(2)?.try_into().unwrap();
            let ty = c.u8()?;
            let value = match ty {
                b'A' => Value::Char(c.u8()?),
                b'f' => Value::Float(c.u32()?),
                b'Z' => Value::String(c.cstr()?.to_vec()),
                b'H' => Value::Hex(c.cstr()?.to_vec()),
                b'B' => {
                    let sub = c.u8()?;
                    let n = c.u32()? as usize;
                    if sub == b'f' {
                        Value::Array(ArrayValue::Float(
                            (0..n).map(|_| c.u32()).collect::<io::Result<_>>()?,
                        ))
                    } else {
                        let it = IntType::from_code(sub)
                            .ok_or_else(|| invalid(format!("bad B subtype {sub}")))?;
                        let bytes = c.take(n * it.width())?;
                        Value::Array(ArrayValue::Int(
                            it,
                            bytes
                                .chunks_exact(it.width())
                                .map(|b| it.decode(b))
                                .collect(),
                        ))
                    }
                }
                _ => {
                    let it = IntType::from_code(ty)
                        .ok_or_else(|| invalid(format!("bad tag type {ty}")))?;
                    Value::Int(it, it.decode(c.take(it.width())?))
                }
            };
            tags.push((tag, value));
        }
        Ok(Record {
            ref_id,
            pos,
            mapq,
            flag,
            next_ref_id,
            next_pos,
            tlen,
            read_name,
            cigar,
            seq,
            qual,
            tags,
        })
    }

    /// The indexing bin, from position and CIGAR, under `rule`. Under
    /// [`BinRule::Sam`] unplaced reads (pos -1) get reg2bin(-1, 0) = 4680;
    /// reads with no reference span use [pos, pos+1).
    pub fn bin(&self, rule: BinRule) -> u16 {
        if rule == BinRule::Htsjdk && self.pos < 0 {
            return 0;
        }
        let beg = self.pos as i64;
        let span = reference_length(&self.cigar) as i64;
        let end = if self.pos < 0 { 0 } else { beg + span.max(1) };
        reg2bin(beg, end)
    }

    /// Append the raw encoding (without `block_size`) to `out`.
    pub fn encode(&self, rule: BinRule, out: &mut Vec<u8>) -> io::Result<()> {
        if self.read_name.len() > 254 {
            return Err(invalid("read name longer than 254 bytes"));
        }
        if self.qual.len() != self.seq.len() {
            return Err(invalid("qual and seq lengths differ"));
        }
        out.extend_from_slice(&self.ref_id.to_le_bytes());
        out.extend_from_slice(&self.pos.to_le_bytes());
        out.push(self.read_name.len() as u8 + 1);
        out.push(self.mapq);
        out.extend_from_slice(&self.bin(rule).to_le_bytes());
        out.extend_from_slice(&(self.cigar.len() as u16).to_le_bytes());
        out.extend_from_slice(&self.flag.to_le_bytes());
        out.extend_from_slice(&(self.seq.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.next_ref_id.to_le_bytes());
        out.extend_from_slice(&self.next_pos.to_le_bytes());
        out.extend_from_slice(&self.tlen.to_le_bytes());
        out.extend_from_slice(&self.read_name);
        out.push(0);
        for op in &self.cigar {
            out.extend_from_slice(&op.to_le_bytes());
        }
        for pair in self.seq.chunks(2) {
            let hi = base_code(pair[0])
                .ok_or_else(|| invalid(format!("bad base {:?}", pair[0] as char)))?;
            let lo = match pair.get(1) {
                Some(&b) => {
                    base_code(b).ok_or_else(|| invalid(format!("bad base {:?}", b as char)))?
                }
                None => 0,
            };
            out.push(hi << 4 | lo);
        }
        out.extend_from_slice(&self.qual);
        for (tag, value) in &self.tags {
            out.extend_from_slice(tag);
            match value {
                Value::Char(ch) => out.extend_from_slice(&[b'A', *ch]),
                Value::Int(it, v) => {
                    out.push(it.code());
                    it.encode(*v, out);
                }
                Value::Float(bits) => {
                    out.push(b'f');
                    out.extend_from_slice(&bits.to_le_bytes());
                }
                Value::String(s) | Value::Hex(s) => {
                    out.push(if matches!(value, Value::String(_)) {
                        b'Z'
                    } else {
                        b'H'
                    });
                    out.extend_from_slice(s);
                    out.push(0);
                }
                Value::Array(ArrayValue::Int(it, vs)) => {
                    out.extend_from_slice(&[b'B', it.code()]);
                    out.extend_from_slice(&(vs.len() as u32).to_le_bytes());
                    for &v in vs {
                        it.encode(v, out);
                    }
                }
                Value::Array(ArrayValue::Float(vs)) => {
                    out.extend_from_slice(b"Bf");
                    out.extend_from_slice(&(vs.len() as u32).to_le_bytes());
                    for v in vs {
                        out.extend_from_slice(&v.to_le_bytes());
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Record {
        Record {
            ref_id: 0,
            pos: 99,
            mapq: 60,
            flag: 99,
            next_ref_id: 0,
            next_pos: 299,
            tlen: 351,
            read_name: b"r001".to_vec(),
            cigar: vec![(5 << 4) | 4, (90 << 4), (2 << 4) | 1, (4 << 4)],
            seq: b"NACGT".iter().copied().cycle().take(101).collect(),
            qual: (0..101).map(|i| (i % 41) as u8).collect(),
            tags: vec![
                (*b"RG", Value::String(b"HG002_WES.L1".to_vec())),
                (*b"NM", Value::Int(IntType::U8, 3)),
                (*b"AS", Value::Int(IntType::I16, -300)),
                (*b"XF", Value::Float(1.5f32.to_bits())),
                (*b"XA", Value::Char(b'z')),
                (
                    *b"ZB",
                    Value::Array(ArrayValue::Int(IntType::I32, vec![-1, 0, 70000])),
                ),
                (
                    *b"ZF",
                    Value::Array(ArrayValue::Float(vec![f32::NAN.to_bits()])),
                ),
                (*b"ZH", Value::Hex(b"1AE301".to_vec())),
            ],
        }
    }

    #[test]
    fn record_roundtrip_with_odd_length_and_all_tag_types() {
        let rec = sample();
        let mut raw = Vec::new();
        rec.encode(BinRule::Sam, &mut raw).unwrap();
        let back = Record::decode(&raw).unwrap();
        assert_eq!(back, rec);
        let mut again = Vec::new();
        back.encode(BinRule::Sam, &mut again).unwrap();
        assert_eq!(again, raw);
    }

    #[test]
    fn bin_matches_spec_examples() {
        // 0-based [99, 196): 5S 90M 2I 4M spans 94 ref bases -> bin of level 5.
        assert_eq!(sample().bin(BinRule::Sam), reg2bin(99, 193));
        assert_eq!(sample().bin(BinRule::Htsjdk), reg2bin(99, 193));
        let unplaced = Record {
            ref_id: -1,
            pos: -1,
            cigar: vec![],
            ..sample()
        };
        assert_eq!(unplaced.bin(BinRule::Sam), 4680);
        assert_eq!(unplaced.bin(BinRule::Htsjdk), 0);
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 16384), 4681);
        assert_eq!(reg2bin(0, 16385), 585);
    }

    #[test]
    fn record_view_matches_decode() {
        let rec = sample();
        let mut raw = Vec::new();
        rec.encode(BinRule::Sam, &mut raw).unwrap();
        let v = RecordView::new(&raw).unwrap();
        assert_eq!(
            (v.ref_id(), v.pos(), v.flag()),
            (rec.ref_id, rec.pos, rec.flag)
        );
        assert_eq!(v.read_name(), &rec.read_name[..]);
        assert_eq!(v.cigar().collect::<Vec<_>>(), rec.cigar);
        assert_eq!(v.reference_length(), reference_length(&rec.cigar));
        assert_eq!(v.mate_ref_id(), rec.next_ref_id);
        assert_eq!(v.quals(), Some(&rec.qual[..]));
        assert_eq!(v.string_tag(b"RG"), Some(&b"HG002_WES.L1"[..]));
        assert_eq!(v.string_tag(b"ZH"), None); // H, not Z
        assert_eq!(v.string_tag(b"QQ"), None); // after skipping every tag type
        assert!(RecordView::new(&raw[..40]).is_err());
    }

    #[test]
    fn record_set_indexes_and_writes_round_trip() {
        let h = Header {
            text: b"@HD\tVN:1.6\n".to_vec(),
            refs: vec![Reference {
                name: b"chr1".to_vec(),
                len: 1000,
            }],
        };
        let mut recs = Vec::new();
        for i in 0..3 {
            let mut r = sample();
            r.pos = 100 + i;
            let mut raw = Vec::new();
            r.encode(BinRule::Htsjdk, &mut raw).unwrap();
            write_raw_record(&mut recs, &raw).unwrap();
        }
        let set = RecordSet::from_parts(h.clone(), &recs).unwrap();
        assert_eq!(set.len(), 3);
        assert_eq!(RecordView::new(set.record(2)).unwrap().pos(), 102);
        let dir = std::env::temp_dir().join(format!("nomics-rs-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.bam");
        set.write(&path, 3).unwrap();
        let back = RecordSet::load(&path, 2).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(back.header, h);
        assert_eq!(
            (0..3).map(|i| back.record(i).to_vec()).collect::<Vec<_>>(),
            (0..3).map(|i| set.record(i).to_vec()).collect::<Vec<_>>()
        );
        assert!(RecordSet::from_stream(set.data[..set.data.len() - 3].to_vec()).is_err());
    }

    #[test]
    fn header_roundtrip() {
        let h = Header {
            text: b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n".to_vec(),
            refs: vec![Reference {
                name: b"chr1".to_vec(),
                len: 248_956_422,
            }],
        };
        let mut buf = Vec::new();
        h.write(&mut buf).unwrap();
        assert_eq!(Header::read(&mut &buf[..]).unwrap(), h);
    }

    #[test]
    fn raw_record_framing() {
        let mut raw = Vec::new();
        sample().encode(BinRule::Sam, &mut raw).unwrap();
        let mut stream = Vec::new();
        write_raw_record(&mut stream, &raw).unwrap();
        write_raw_record(&mut stream, &raw).unwrap();
        let mut r = &stream[..];
        let mut buf = Vec::new();
        assert!(read_raw_record(&mut r, &mut buf).unwrap());
        assert_eq!(buf, raw);
        assert!(read_raw_record(&mut r, &mut buf).unwrap());
        assert!(!read_raw_record(&mut r, &mut buf).unwrap());
    }
}
