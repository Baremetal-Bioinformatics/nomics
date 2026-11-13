//! BAI-selected BAM stream. Chunks from all requested intervals are unioned
//! before decoding, so a record is returned at most once, in file order.
//! The caller must still test exact record/interval overlap (BAI bins are coarse).
use crate::{
    bam::{self, Header},
    bgzf::BgzfReader,
};
use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
};

/// A query is `(reference index, 0-based inclusive start, exclusive end)`.
pub struct BamQuery {
    pub header: Header,
    reader: BgzfReader<BufReader<File>>,
    chunks: Vec<(u64, u64)>,
    at: usize,
    positioned: bool,
}

fn u32le(r: &mut impl Read) -> io::Result<u32> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
fn u64le(r: &mut impl Read) -> io::Result<u64> {
    let mut b = [0; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
impl BamQuery {
    /// Open the BAM and its `.bam.bai` (or sibling `.bai`). Empty queries
    /// return no records. Coordinates must be nonempty and fit BAI's 2^29 range.
    pub fn open(path: &Path, intervals: &[(usize, u32, u32)]) -> io::Result<Self> {
        let mut reader = BgzfReader::new(BufReader::new(File::open(path)?));
        let header = Header::read(&mut reader)?;
        if intervals
            .iter()
            .any(|&(tid, beg, end)| tid >= header.refs.len() || beg >= end || end > 1 << 29)
        {
            return Err(io::Error::other("invalid BAI query interval"));
        }
        let mut index_path = path.as_os_str().to_owned();
        index_path.push(".bai");
        let index = File::open(Path::new(&index_path)).or_else(|e| {
            if e.kind() == io::ErrorKind::NotFound {
                File::open(path.with_extension("bai"))
            } else {
                Err(e)
            }
        })?;
        let mut index = BufReader::new(index);
        if u32le(&mut index)? != u32::from_le_bytes(*b"BAI\x01") {
            return Err(io::Error::other("not a BAI index"));
        }
        let nref = u32le(&mut index)? as usize;
        if nref != header.refs.len() {
            return Err(io::Error::other("BAI reference count mismatch"));
        }
        let mut chunks = Vec::new();
        for tid in 0..nref {
            // No linear-index cutoff: retaining coarse-bin chunks is conservative
            // for long alignments and sparse intervals, with exact filtering later.
            let mut wanted_bins = std::collections::BTreeSet::new();
            for &(_, beg, end) in intervals.iter().filter(|i| i.0 == tid) {
                wanted_bins.insert(0);
                for (shift, offset) in [(26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)] {
                    wanted_bins.extend(offset + (beg >> shift)..=offset + ((end - 1) >> shift));
                }
            }
            for _ in 0..u32le(&mut index)? {
                let bin = u32le(&mut index)?;
                let wanted = wanted_bins.contains(&bin);
                for _ in 0..u32le(&mut index)? {
                    let chunk = (u64le(&mut index)?, u64le(&mut index)?);
                    if wanted && chunk.0 < chunk.1 {
                        chunks.push(chunk);
                    }
                }
            }
            for _ in 0..u32le(&mut index)? {
                u64le(&mut index)?;
            }
        }
        chunks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (beg, end) in chunks {
            if let Some(last) = merged.last_mut().filter(|last| beg <= last.1) {
                last.1 = last.1.max(end);
            } else {
                merged.push((beg, end));
            }
        }
        Ok(Self {
            header,
            reader,
            chunks: merged,
            at: 0,
            positioned: false,
        })
    }
    /// Virtual offset immediately after the most recently returned record.
    /// Together with the source file, this identifies an input record exactly.
    pub fn virtual_offset(&self) -> u64 {
        self.reader.virtual_offset()
    }
    /// Raw BAM records without block_size, borrowing the caller's reusable buffer.
    pub fn next_record(&mut self, raw: &mut Vec<u8>) -> io::Result<bool> {
        while let Some(&(beg, end)) = self.chunks.get(self.at) {
            if !self.positioned {
                self.reader.seek_virtual(beg)?;
                self.positioned = true;
            }
            if self.reader.virtual_offset() < end {
                return bam::read_raw_record(&mut self.reader, raw);
            }
            self.at += 1;
            self.positioned = false;
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::{BinRule, Record, RecordView};
    #[test]
    fn indexed_union_matches_sequential_and_crosses_blocks() {
        let path = std::env::temp_dir().join(format!("nomics-query-{}.bam", std::process::id()));
        let header = Header {
            text: b"@HD\tSO:coordinate\n".to_vec(),
            refs: vec![bam::Reference {
                name: b"chr1".to_vec(),
                len: 100000,
            }],
        };
        let mut w = crate::bgzf::BgzfWriter::new(
            File::create(&path).unwrap(),
            crate::bgzf::Stored,
            crate::bgzf::Profile::HTSJDK,
        );
        header.write(&mut w).unwrap();
        let mut record = Record {
            ref_id: 0,
            pos: 0,
            mapq: 60,
            flag: 0,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
            read_name: b"r".to_vec(),
            cigar: vec![100 << 4],
            seq: vec![b'A'; 100],
            qual: vec![30; 100],
            tags: vec![],
        };
        let mut raw = Vec::new();
        for pos in (0..100000).step_by(10) {
            record.pos = pos;
            raw.clear();
            record.encode(BinRule::Htsjdk, &mut raw).unwrap();
            bam::write_raw_record(&mut w, &raw).unwrap();
        }
        w.finish().unwrap();
        let mut bai = path.as_os_str().to_owned();
        bai.push(".bai");
        crate::htsjdk_index::write_bai_htsjdk(&path, Path::new(&bai)).unwrap();
        let intervals = [(0, 16380, 16400), (0, 16390, 16450), (0, 99999, 100000)];
        let mut q = BamQuery::open(&path, &intervals).unwrap();
        let mut found = Vec::new();
        let mut prev = -1;
        while q.next_record(&mut raw).unwrap() {
            let r = RecordView::new(&raw).unwrap();
            assert!(r.pos() > prev);
            prev = r.pos();
            if intervals.iter().any(|&(_, b, e)| {
                r.pos() < e as i32 && r.pos() + r.reference_length() as i32 > b as i32
            }) {
                found.push(r.pos());
            }
        }
        let expected: Vec<_> = (0..100000)
            .step_by(10)
            .filter(|&p| {
                intervals
                    .iter()
                    .any(|&(_, b, e)| p < e as i32 && p + 100 > b as i32)
            })
            .collect();
        assert_eq!(found, expected);
        assert!(
            !BamQuery::open(&path, &[])
                .unwrap()
                .next_record(&mut raw)
                .unwrap()
        );
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(bai).unwrap();
    }
}
