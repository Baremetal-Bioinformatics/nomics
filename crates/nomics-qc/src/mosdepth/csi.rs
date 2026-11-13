//! hts-nim's on-the-fly CSI layout, including htslib's bin table order.
use nomics_io::bgzf::{BgzfReader, BgzfWriter, Profile};
use nomics_io::deflate::FastDeflater;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

const PAYLOAD: Profile = Profile {
    block_payload: 65_280,
};
fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}
fn first(level: u32) -> u32 {
    ((1u32 << (3 * level)) - 1) / 7
}
fn bin(beg: u32, end: u32, levels: u32) -> u32 {
    let end = i64::from(end) - 1;
    for level in (1..=levels).rev() {
        let shift = 14 + 3 * (levels - level);
        if i64::from(beg) >> shift == end >> shift {
            return first(level) + (beg >> shift);
        }
    }
    0
}
#[derive(Default, Clone)]
struct Bin {
    loff: u64,
    chunks: Vec<(u64, u64)>,
}

// Integer identity hash, triangular probing, 0.77 load, and in-place rehash
// kick-outs determine the serialized bin order. No upstream code is linked.
#[derive(Default)]
struct Bins {
    slots: Vec<Option<(u32, Bin)>>,
    count: usize,
}
impl Bins {
    fn find(&self, key: u32) -> Option<usize> {
        if self.slots.is_empty() {
            return None;
        }
        let mask = self.slots.len() - 1;
        let mut at = key as usize & mask;
        let mut step = 0;
        loop {
            match &self.slots[at] {
                None => return None,
                Some((k, _)) if *k == key => return Some(at),
                _ => {}
            }
            step += 1;
            at = (at + step) & mask;
        }
    }
    fn grow(&mut self) {
        let old_len = self.slots.len();
        let new_len = (old_len * 2).max(4);
        self.slots.resize_with(new_len, || None);
        let mut occupied = vec![false; new_len];
        for j in 0..old_len {
            if occupied[j] {
                continue;
            }
            if let Some(mut item) = self.slots[j].take() {
                loop {
                    let mut at = item.0 as usize & (new_len - 1);
                    let mut step = 0;
                    while occupied[at] {
                        step += 1;
                        at = (at + step) & (new_len - 1);
                    }
                    occupied[at] = true;
                    let displaced = self.slots[at].replace(item);
                    match displaced {
                        Some(next) => item = next,
                        None => break,
                    }
                }
            }
        }
    }
    fn append(&mut self, key: u32, chunk: (u64, u64)) {
        // kh_put grows even for an existing key at the load threshold.
        if self.count >= (self.slots.len() as f64 * 0.77 + 0.5) as usize {
            self.grow();
        }
        let at = if let Some(at) = self.find(key) {
            at
        } else {
            let mut at = key as usize & (self.slots.len() - 1);
            let mut step = 0;
            while self.slots[at].is_some() {
                step += 1;
                at = (at + step) & (self.slots.len() - 1);
            }
            self.slots[at] = Some((key, Bin::default()));
            self.count += 1;
            at
        };
        self.slots[at].as_mut().unwrap().1.chunks.push(chunk);
    }
    fn finish(&mut self, linear: &mut [u64], levels: u32) {
        for i in (0..linear.len().saturating_sub(1)).rev() {
            if linear[i] == u64::MAX {
                linear[i] = linear[i + 1];
            }
        }
        let n_bins = first(levels + 1);
        for (key, b) in self.slots.iter_mut().flatten() {
            if *key >= n_bins {
                continue;
            }
            let mut level = 0;
            let mut parent = *key;
            while parent > 0 {
                parent = (parent - 1) >> 3;
                level += 1;
            }
            let bottom = (*key - first(level)) << (3 * (levels - level));
            b.loff = linear.get(bottom as usize).copied().unwrap_or(0);
        }
        // No new bins are inserted during compaction, so deleted slots need
        // no tombstone lookup: locate parents by their keys directly.
        for level in (1..=levels).rev() {
            for at in 0..self.slots.len() {
                let Some((key, b)) = &mut self.slots[at] else {
                    continue;
                };
                if *key >= n_bins || *key < first(level) {
                    continue;
                }
                if level < levels && b.chunks.len() > 1 {
                    b.chunks.sort_unstable();
                }
                if (b.chunks.last().unwrap().1 >> 16) - (b.chunks[0].0 >> 16) >= 0x10000 {
                    continue;
                }
                let parent = (*key - 1) >> 3;
                let Some(p) = self
                    .slots
                    .iter()
                    .position(|x| x.as_ref().is_some_and(|x| x.0 == parent))
                else {
                    continue;
                };
                let child = self.slots[at].take().unwrap().1;
                self.slots[p]
                    .as_mut()
                    .unwrap()
                    .1
                    .chunks
                    .extend(child.chunks);
                self.count -= 1;
            }
        }
        for (key, b) in self.slots.iter_mut().flatten() {
            if *key >= n_bins {
                continue;
            }
            if *key == 0 {
                b.chunks.sort_unstable();
            }
            let mut merged: Vec<(u64, u64)> = Vec::new();
            for &(s, e) in &b.chunks {
                if let Some(last) = merged.last_mut()
                    && last.1 >> 16 >= s >> 16
                {
                    last.1 = last.1.max(e);
                    continue;
                }
                merged.push((s, e));
            }
            b.chunks = merged;
        }
    }
}
#[derive(Default)]
struct Contig {
    name: String,
    bins: Bins,
    linear: Vec<u64>,
    last_bin: Option<u32>,
    save_off: u64,
    begin: u64,
    count: u64,
    last_start: u32,
}
impl Contig {
    fn add(&mut self, start: u32, stop: u32, before: u64, levels: u32) -> io::Result<()> {
        if self.count > 0 && start < self.last_start {
            return Err(invalid("CSI input starts out of order"));
        }
        if stop < start || u64::from(stop) > 1u64 << (14 + levels * 3) {
            return Err(invalid("CSI interval outside index range"));
        }
        let stop = stop.max(1);
        let lo = (start >> 14) as usize;
        let hi = ((stop - 1) >> 14) as usize;
        if hi >= self.linear.len() {
            self.linear.resize(hi + 1, u64::MAX);
        }
        if lo <= hi {
            for v in &mut self.linear[lo..=hi] {
                if *v == u64::MAX {
                    *v = before;
                }
            }
        }
        let b = bin(start, stop, levels);
        if self.last_bin != Some(b) {
            if let Some(old) = self.last_bin {
                self.bins.append(old, (self.save_off, before));
            }
            self.save_off = before;
            self.last_bin = Some(b);
        }
        self.last_start = start;
        self.count += 1;
        Ok(())
    }
    fn finish(&mut self, end: u64, levels: u32) {
        if let Some(b) = self.last_bin {
            self.bins.append(b, (self.save_off, end));
        }
        let meta = first(levels + 1) + 1;
        self.bins.append(meta, (self.begin, end));
        self.bins.append(meta, (self.count, 0));
        self.bins.finish(&mut self.linear, levels);
    }
}

/// Rebuild the index from the file's own BGZF virtual offsets. `levels` is
/// computed from the BAM header by mosdepth; standalone callers must supply
/// that same value (5 for a GRCh38 header), which BED text cannot always infer.
pub fn build(bed_gz: &Path, out: &Path, levels: u32) -> io::Result<()> {
    if levels > 9 {
        return Err(invalid("CSI levels exceeds 32-bit bin range"));
    }
    let mut reader = BgzfReader::new(File::open(bed_gz)?);
    let mut contigs: Vec<Contig> = Vec::new();
    let mut line = Vec::new();
    let mut before = 0;
    let mut final_offset = 0;
    while let Some(info) = reader.next_block()? {
        let data = reader.block();
        if data.is_empty() {
            continue;
        }
        final_offset = (info.coffset + u64::from(info.csize)) << 16;
        let mut at = 0;
        while at < data.len() {
            let n = data[at..].iter().position(|&b| b == b'\n');
            let end = n.map_or(data.len(), |n| at + n);
            line.extend_from_slice(&data[at..end]);
            at = end;
            if n.is_none() {
                break;
            }
            at += 1;
            let after = if at == data.len() {
                final_offset
            } else {
                info.coffset << 16 | at as u64
            };
            let text = std::str::from_utf8(&line).map_err(|_| invalid("non-UTF8 BED"))?;
            let mut fields = text.split('\t');
            let name = fields.next().unwrap_or("");
            let start: u32 = fields
                .next()
                .ok_or_else(|| invalid("missing BED start"))?
                .parse()
                .map_err(|_| invalid("invalid BED start"))?;
            let stop: u32 = fields
                .next()
                .ok_or_else(|| invalid("missing BED end"))?
                .parse()
                .map_err(|_| invalid("invalid BED end"))?;
            if contigs.last().is_none_or(|c| c.name != name) {
                if contigs.iter().any(|c| c.name == name) {
                    return Err(invalid("CSI chromosomes not continuous"));
                }
                if let Some(c) = contigs.last_mut() {
                    c.finish(before, levels);
                }
                contigs.push(Contig {
                    name: name.into(),
                    begin: before,
                    ..Contig::default()
                });
            }
            contigs
                .last_mut()
                .unwrap()
                .add(start, stop, before, levels)?;
            before = after;
            line.clear();
        }
    }
    if !line.is_empty() {
        return Err(invalid("BED missing final newline"));
    }
    if let Some(c) = contigs.last_mut() {
        c.finish(final_offset, levels);
    }
    let mut w = BgzfWriter::new(
        BufWriter::new(File::create(out)?),
        FastDeflater::default(),
        PAYLOAD,
    );
    w.write_all(b"CSI\x01")?;
    let names_len: usize = contigs.iter().map(|c| c.name.len() + 1).sum();
    for v in [
        14,
        levels,
        28 + names_len as u32,
        0x10000,
        1,
        2,
        3,
        35,
        0,
        names_len as u32,
    ] {
        w.write_all(&v.to_le_bytes())?;
    }
    for c in &contigs {
        w.write_all(c.name.as_bytes())?;
        w.write_all(&[0])?;
    }
    w.write_all(&(contigs.len() as u32).to_le_bytes())?;
    for c in &contigs {
        w.write_all(&(c.bins.count as u32).to_le_bytes())?;
        for (key, b) in c.bins.slots.iter().flatten() {
            w.write_all(&key.to_le_bytes())?;
            w.write_all(&b.loff.to_le_bytes())?;
            w.write_all(&(b.chunks.len() as u32).to_le_bytes())?;
            for (s, e) in &b.chunks {
                w.write_all(&s.to_le_bytes())?;
                w.write_all(&e.to_le_bytes())?;
            }
        }
    }
    w.write_all(&0u64.to_le_bytes())?;
    w.finish()?;
    Ok(())
}
