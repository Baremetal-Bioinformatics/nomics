//! bwa SAM text -> coordinate-sorted BAM, as `samtools sort` 1.22.1 does it
//! (htslib `sam_parse1`/`aux_parse`, `bam_sort.c`). See docs/spec/bwa.md §9.
//!
//! Gate: per-chunk `preprocessing/mapped/*/*.NNNN.bam` records.

use std::collections::HashMap;
use std::io;

use nomics_io::bam::{BinRule, Header, IntType, Record, RecordSet, Reference, Value};

use crate::index::Index;
use crate::sam::{BamRecords, ChunkRecords};

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// The sorted-BAM header: bwa's lines with `@HD` rewritten by samtools
/// (`SO:coordinate`, no `GO`), then our own `@PG`.
pub fn header(idx: &Index, rg_line: &str, pg_cl: &str) -> Header {
    let mut text = String::from("@HD\tVN:1.5\tSO:coordinate\n");
    let mut refs = Vec::with_capacity(idx.bns.anns.len());
    for a in &idx.bns.anns {
        text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", a.name, a.len));
        refs.push(Reference {
            name: a.name.as_bytes().to_vec(),
            len: a.len as u32,
        });
    }
    text.push_str(rg_line);
    text.push('\n');
    // Two program records per chunk, aligner then sorter, as in the
    // reference's `bwa mem | samtools sort`. The header's @PG count matters
    // downstream: samtools merge draws one random ID suffix per @PG line.
    text.push_str(&format!(
        "@PG\tID:nomics-align\tPN:nomics-align\tVN:{v}\tCL:{pg_cl}\n\
         @PG\tID:nomics-sort\tPN:nomics-sort\tPP:nomics-align\tVN:{v}\tCL:coordinate sort\n",
        v = env!("CARGO_PKG_VERSION"),
    ));
    Header {
        text: text.into_bytes(),
        refs,
    }
}

/// htslib's choice of integer type for a SAM `i` tag value.
fn int_type(v: i64) -> IntType {
    if v < 0 {
        if v >= -128 {
            IntType::I8
        } else if v >= -32768 {
            IntType::I16
        } else {
            IntType::I32
        }
    } else if v <= 255 {
        IntType::U8
    } else if v <= 65535 {
        IntType::U16
    } else {
        IntType::U32
    }
}

fn num<T: std::str::FromStr>(f: &[u8], what: &str) -> io::Result<T> {
    std::str::from_utf8(f)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| invalid(format!("bad {what}: {:?}", String::from_utf8_lossy(f))))
}

fn parse_cigar(f: &[u8]) -> io::Result<Vec<u32>> {
    let mut out = Vec::new();
    if f == b"*" {
        return Ok(out);
    }
    let mut n = 0u32;
    for &c in f {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else {
            let op = nomics_io::bam::CIGAR_OPS
                .iter()
                .position(|&o| o == c)
                .ok_or_else(|| invalid(format!("bad CIGAR op {:?}", c as char)))?;
            out.push(n << 4 | op as u32);
            n = 0;
        }
    }
    Ok(out)
}

/// `sam_parse1` for the fields and tags bwa writes.
pub fn parse_line(line: &[u8], tids: &HashMap<&[u8], i32>) -> io::Result<Record> {
    let f: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
    if f.len() < 11 {
        return Err(invalid("SAM line with fewer than 11 fields"));
    }
    let tid = |name: &[u8]| -> io::Result<i32> {
        if name == b"*" {
            Ok(-1)
        } else {
            tids.get(name).copied().ok_or_else(|| {
                invalid(format!(
                    "unknown reference {:?}",
                    String::from_utf8_lossy(name)
                ))
            })
        }
    };
    let ref_id = tid(f[2])?;
    let mut flag: u16 = num(f[1], "FLAG")?;
    let cigar = parse_cigar(f[5])?;
    if ref_id < 0 || cigar.is_empty() {
        flag |= 0x4;
    }
    let next_ref_id = if f[6] == b"=" { ref_id } else { tid(f[6])? };
    let seq = f[9].to_vec();
    let qual = if f[10] == b"*" {
        vec![0xff; seq.len()]
    } else {
        f[10].iter().map(|&q| q - 33).collect()
    };
    let mut tags = Vec::new();
    for t in &f[11..] {
        if t.len() < 5 || t[2] != b':' || t[4] != b':' {
            return Err(invalid("bad aux field"));
        }
        let key = [t[0], t[1]];
        let v = &t[5..];
        let value = match t[3] {
            b'i' => {
                let x: i64 = num(v, "integer tag")?;
                Value::Int(int_type(x), x)
            }
            b'Z' => Value::String(v.to_vec()),
            b'f' => Value::Float(num::<f32>(v, "float tag")?.to_bits()),
            c => return Err(invalid(format!("unsupported aux type {:?}", c as char))),
        };
        tags.push((key, value));
    }
    Ok(Record {
        ref_id,
        pos: num::<i32>(f[3], "POS")? - 1,
        mapq: num(f[4], "MAPQ")?,
        flag,
        next_ref_id,
        next_pos: num::<i32>(f[7], "PNEXT")? - 1,
        tlen: num(f[8], "TLEN")?,
        read_name: f[0].to_vec(),
        cigar,
        seq,
        qual,
        tags,
    })
}

/// Parse bwa's SAM records and sort them as `samtools sort` does: stably by
/// (tid with -1 last, POS, reverse strand). Parsing and encoding run on
/// `threads` threads over line-aligned slices; the sort key carries each
/// record's input index, so the order is the stable one.
pub fn sort_sam(header: Header, sam: &[u8], threads: usize) -> io::Result<RecordSet> {
    let tids: HashMap<&[u8], i32> = header
        .refs
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.as_slice(), i as i32))
        .collect();
    let n_ref = header.refs.len() as u64;

    // Line-aligned slices.
    let threads = threads.max(1);
    let target = sam.len().div_ceil(threads * 8).max(1);
    let mut slices = Vec::new();
    let mut at = 0;
    while at < sam.len() {
        let mut end = (at + target).min(sam.len());
        while end < sam.len() && sam[end - 1] != b'\n' {
            end += 1;
        }
        slices.push(&sam[at..end]);
        at = end;
    }

    // Encode each slice: block_size-prefixed records plus (key, offset).
    type Part = io::Result<(Vec<u8>, Vec<(u64, u64)>)>;
    let mut parts: Vec<Option<Part>> = (0..slices.len()).map(|_| None).collect();
    {
        let work: Vec<_> = slices.iter().zip(parts.iter_mut()).collect();
        let queue = std::sync::Mutex::new(work);
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(|| {
                    loop {
                        let Some((slice, slot)) = queue.lock().unwrap().pop() else {
                            break;
                        };
                        let encode = || -> Part {
                            let mut bytes = Vec::with_capacity(slice.len() / 2);
                            let mut keys = Vec::new();
                            let mut buf = Vec::new();
                            for line in slice.split(|&c| c == b'\n').filter(|l| !l.is_empty()) {
                                let r = parse_line(line, &tids)?;
                                let tid = if r.ref_id < 0 { n_ref } else { r.ref_id as u64 };
                                let key = tid << 33
                                    | ((r.pos as i64 + 1) as u64) << 1
                                    | (r.flag >> 4 & 1) as u64;
                                keys.push((key, bytes.len() as u64));
                                buf.clear();
                                r.encode(BinRule::Sam, &mut buf)?;
                                bytes.extend_from_slice(&(buf.len() as u32).to_le_bytes());
                                bytes.extend_from_slice(&buf);
                            }
                            Ok((bytes, keys))
                        };
                        *slot = Some(encode());
                    }
                });
            }
        });
    }
    let parts: Vec<(Vec<u8>, Vec<(u64, u64)>)> = parts
        .into_iter()
        .map(|p| p.unwrap())
        .collect::<Result<_, _>>()?;

    // Stable order: (key, slice, index within slice).
    let mut order: Vec<(u64, u32, u64)> = Vec::new();
    for (pi, (_, keys)) in parts.iter().enumerate() {
        order.extend(keys.iter().map(|&(k, off)| (k, pi as u32, off)));
    }
    order.sort_unstable();
    let total: usize = parts.iter().map(|p| p.0.len()).sum();
    let mut out = Vec::with_capacity(total);
    for &(_, pi, off) in &order {
        let b = &parts[pi as usize].0;
        let off = off as usize;
        let n = u32::from_le_bytes(b[off..off + 4].try_into().unwrap()) as usize;
        out.extend_from_slice(&b[off..off + 4 + n]);
    }
    RecordSet::from_parts(header, &out)
}

/// Sort records emitted directly by the aligner with the same stable key as
/// [`sort_sam`]. Their encodings already use htslib's SAM conversion rules.
pub fn sort_bam(header: Header, mut records: BamRecords) -> io::Result<RecordSet> {
    if let Some(error) = records.take_error() {
        return Err(error);
    }
    let mut order: Vec<(u64, u32, u64)> = records
        .keys
        .iter()
        .enumerate()
        .map(|(ordinal, &(key, off))| (key, ordinal as u32, off))
        .collect();
    order.sort_unstable();
    // Build the set's stream in place (header, then records in order) and
    // free the unsorted copy before indexing: at most two copies of a chunk
    // are resident, which matters for whole-genome chunks.
    let mut out = Vec::new();
    header.write(&mut out)?;
    out.reserve_exact(records.bytes.len());
    for &(_, _, off) in &order {
        let off = off as usize;
        let n = u32::from_le_bytes(records.bytes[off..off + 4].try_into().unwrap()) as usize;
        out.extend_from_slice(&records.bytes[off..off + 4 + n]);
    }
    drop(order);
    drop(records);
    RecordSet::from_stream(out)
}

/// [`sort_bam`] over a chunk accumulated by reference ([`ChunkRecords`]):
/// the same order, since a key's top bits are its reference id, and within
/// one reference `(key, batch, record)` orders exactly as `(key, record
/// ordinal)` does over the whole chunk. Each reference's pages are returned
/// to the system as soon as its records are in the sorted stream, so the
/// chunk is held about once rather than twice.
pub fn sort_chunk(header: Header, chunk: ChunkRecords) -> io::Result<RecordSet> {
    // Header plus records, so the last record never doubles the allocation.
    let mut head = Vec::new();
    header.write(&mut head)?;
    let mut out = Vec::new();
    out.reserve_exact(head.len() + chunk_bytes(&chunk));
    sort_chunk_into(header, chunk, &mut out)?;
    RecordSet::from_stream(out)
}

/// [`sort_chunk`] straight into a BGZF BAM file on `threads` threads, the
/// same bytes [`RecordSet::write`] gives: BGZF blocks are cut by byte count,
/// not by write calls. The sorted chunk is never held in memory as a whole:
/// records stream from the chunk's pages, which are freed reference by
/// reference. Returns the number of records. A path with the spill
/// extension ([`nomics_io::spill::is_spill_path`]) gets the same stream in
/// the light spill format instead.
pub fn write_sorted_chunk(
    header: Header,
    chunk: ChunkRecords,
    path: &std::path::Path,
    threads: usize,
) -> io::Result<usize> {
    // An unpublished chunk is only markdup's spill: the same stream in the
    // light spill format, which is cheaper to write and to read back. The
    // file gets its name only once complete.
    if nomics_io::spill::is_spill_path(path) {
        let mut w = nomics_io::spill::SpillFile::create(path, threads)?;
        let records = sort_chunk_into(header, chunk, &mut w)?;
        w.finish()?;
        return Ok(records);
    }
    let f = std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(path)?);
    let mut w = nomics_io::bgzf::ParallelBgzfWriter::for_bam(f, threads);
    let records = sort_chunk_into(header, chunk, &mut w)?;
    std::io::Write::flush(&mut w.finish()?)?;
    Ok(records)
}

/// Record bytes in a chunk (without the sort keys).
fn chunk_bytes(chunk: &ChunkRecords) -> usize {
    chunk
        .refs
        .iter()
        .flatten()
        .map(|(_, s)| s.len - 8 * s.count)
        .sum()
}

/// The sorted stream (header, then records) into `out`; returns the
/// number of records.
fn sort_chunk_into(
    header: Header,
    mut chunk: ChunkRecords,
    out: &mut impl std::io::Write,
) -> io::Result<usize> {
    if let Some(error) = chunk.error.take() {
        return Err(error);
    }
    // Tie order equals sort_bam's (key, u32 ordinal) only below 2^32
    // records per chunk; a 30x chunk holds ~58M. Refuse rather than diverge.
    let total: usize = chunk.refs.iter().flatten().map(|(_, s)| s.count).sum();
    if total > u32::MAX as usize {
        return Err(io::Error::other(format!(
            "chunk of {total} records exceeds the 2^32 supported by the sort order"
        )));
    }
    let mut head = Vec::new();
    header.write(&mut head)?;
    out.write_all(&head)?;
    let page = nomics_mem::page_size();
    let mut records = 0;
    // (key, bin, record, batch, offset of the record in its batch).
    let mut order: Vec<(u64, u32, u32, u32, usize)> = Vec::new();
    for bins in &chunk.refs {
        order.clear();
        for (i, &(b, span)) in bins.iter().enumerate() {
            let buf = chunk.batches[b].as_slice();
            let keys = &buf[span.start..span.start + 8 * span.count];
            let mut at = span.start + 8 * span.count;
            for (r, key) in keys.as_chunks::<8>().0.iter().enumerate() {
                let key = u64::from_le_bytes(*key);
                order.push((key, i as u32, r as u32, b as u32, at));
                at += 4 + u32::from_le_bytes(buf[at..at + 4].try_into().unwrap()) as usize;
            }
            debug_assert_eq!(at, span.start + span.len);
        }
        order.sort_unstable();
        records += order.len();
        for &(_, _, _, b, at) in &order {
            let buf = chunk.batches[b as usize].as_slice();
            let n = 4 + u32::from_le_bytes(buf[at..at + 4].try_into().unwrap()) as usize;
            out.write_all(&buf[at..at + n])?;
        }
        for &(b, span) in bins {
            // Spans start on a page and the next one starts on the next page.
            chunk.batches[b].release(span.start, span.len.next_multiple_of(page));
        }
    }
    Ok(records)
}
