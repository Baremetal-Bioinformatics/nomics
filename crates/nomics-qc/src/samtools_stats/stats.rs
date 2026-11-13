use super::{invalid, ring::Ring};
use nomics_io::crc32::crc32;
use std::io;

/// A decoded record boundary independent of BAM framing, for a future CRAM
/// adapter. `packed` is BAM's sequence encoding, including its odd padding
/// nibble: CHK hashes those bytes, not decoded letters.
pub struct Record<'a> {
    pub name: &'a [u8],
    pub packed: &'a [u8],
    pub quals: &'a [u8],
    pub cigar: &'a [u32],
    pub tid: i32,
    pub pos: i64,
    pub mate_tid: i32,
    pub mate_pos: i64,
    pub tlen: i64,
    pub flag: u16,
    pub mapq: u8,
    pub nm: Option<i64>,
    pub barcode: [Option<&'a [u8]>; 4],
    pub barqual: [Option<&'a [u8]>; 4],
}
impl Record<'_> {
    fn base(&self, i: usize) -> usize {
        ((self.packed[i / 2] >> (4 * (1 - i % 2))) & 15) as usize
    }
}

#[derive(Default)]
pub(super) struct Barcode {
    pub bases: Vec<[u64; 6]>,
    pub quals: Vec<[u64; 257]>,
    pub sep: Option<usize>,
    pub maxq: Option<usize>,
}
impl Barcode {
    fn add(&mut self, seq: &[u8], qual: Option<&[u8]>) {
        if seq.is_empty() {
            return;
        }
        if self.bases.is_empty() {
            self.bases.resize(seq.len(), [0; 6]);
            self.quals.resize(seq.len(), [0; 257]);
        }
        if seq.len() > self.bases.len() {
            return;
        }
        for (i, &c) in seq.iter().enumerate() {
            let b = match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'N' => 4,
                _ => 5,
            };
            if b == 5 {
                if let Some(s) = self.sep {
                    if s != i {
                        return;
                    }
                } else {
                    self.sep = Some(i);
                }
            } else {
                self.bases[i][b] += 1;
            }
        }
        if let Some(q) = qual.filter(|q| q.len() == seq.len()) {
            for (i, &x) in q.iter().enumerate() {
                if let Some(v) = x.checked_sub(33) {
                    self.quals[i][v as usize] += 1;
                    self.maxq = Some(self.maxq.unwrap_or(0).max(v as usize));
                }
            }
        }
    }
}

/// Streaming default-mode accumulator. No region, flag, read-group or length
/// filters are applied by the pinned command.
pub struct Stats {
    pub(super) chk: [u32; 3],
    pub(super) nbases: usize,
    pub(super) maxlen: [usize; 3],
    pub(super) maxq: usize,
    pub(super) lens: Vec<[u64; 3]>,
    pub(super) quals: [Vec<[u64; 256]>; 2],
    pub(super) bases: [Vec<[u64; 6]>; 2],
    pub(super) oriented: Vec<[u64; 4]>,
    pub(super) gc: [[u64; 200]; 2],
    pub(super) barcodes: [Barcode; 4],
    pub(super) reads: [u64; 3],
    pub(super) length: [u64; 3],
    pub(super) sumq: f64,
    pub(super) dup: u64,
    pub(super) duplen: u64,
    pub(super) unmapped: u64,
    pub(super) mapped: u64,
    pub(super) paired_mapped: u64,
    pub(super) proper: u64,
    pub(super) paired: u64,
    pub(super) anomalous: u64,
    pub(super) mq0: u64,
    pub(super) qcfail: u64,
    pub(super) secondary: u64,
    pub(super) supplementary: u64,
    pub(super) mapped_bases: u64,
    pub(super) cigar_bases: u64,
    pub(super) nm: u64,
    pub(super) mapq: [u64; 256],
    pub(super) isize: Vec<[u64; 3]>,
    pub(super) indels: [[u64; 2]; 300],
    pub(super) cycles: Vec<[u64; 4]>,
    pub(super) sorted: bool,
    tid: i32,
    pos: i64,
    gcd_pos: i64,
    pub(super) gcd: Vec<(f32, u32)>,
    pub(super) ring: Ring,
}
impl Default for Stats {
    fn default() -> Self {
        Self::new()
    }
}
impl Stats {
    pub fn new() -> Self {
        Self {
            chk: [0; 3],
            nbases: 300,
            maxlen: [0; 3],
            maxq: 0,
            lens: vec![[0; 3]; 301],
            quals: std::array::from_fn(|_| vec![[0; 256]; 300]),
            bases: std::array::from_fn(|_| vec![[0; 6]; 300]),
            oriented: vec![[0; 4]; 300],
            gc: [[0; 200]; 2],
            barcodes: std::array::from_fn(|_| Barcode::default()),
            reads: [0; 3],
            length: [0; 3],
            sumq: 0.,
            dup: 0,
            duplen: 0,
            unmapped: 0,
            mapped: 0,
            paired_mapped: 0,
            proper: 0,
            paired: 0,
            anomalous: 0,
            mq0: 0,
            qcfail: 0,
            secondary: 0,
            supplementary: 0,
            mapped_bases: 0,
            cigar_bases: 0,
            nm: 0,
            mapq: [0; 256],
            isize: vec![[0; 3]; 8001],
            indels: [[0; 2]; 300],
            cycles: vec![[0; 4]; 301],
            sorted: true,
            tid: -1,
            pos: 0,
            gcd_pos: -1,
            gcd: vec![(0., 0)],
            ring: Ring::new(),
        }
    }
    /// Account for one record, preserving upstream's counter/filter ordering.
    pub fn add(&mut self, r: &Record<'_>) -> io::Result<()> {
        let n = r.quals.len();
        if r.packed.len() != n.div_ceil(2) {
            return Err(invalid("packed sequence length disagrees with qualities"));
        }
        for (sum, bytes) in self.chk.iter_mut().zip([r.name, r.packed, r.quals]) {
            *sum = sum.wrapping_add(crc32(bytes));
        }
        let f = r.flag;
        if f & 0x100 != 0 {
            self.secondary += 1;
            return Ok(());
        }
        if f & 0x800 != 0 {
            self.supplementary += 1;
        }
        if n == 0 {
            return Ok(());
        }
        let original = f & 0x900 == 0;
        if f & 0x400 != 0 {
            self.dup += 1;
            self.duplen += n as u64;
        }
        let order = if f & 1 == 0 {
            1
        } else {
            usize::from(f & 64 != 0) + 2 * usize::from(f & 128 != 0)
        };
        let group = match order {
            1 => Some(0),
            2 => Some(1),
            _ => None,
        };
        let reverse = f & 16 != 0;
        let readlen = n + r
            .cigar
            .iter()
            .filter(|c| **c & 15 == 5)
            .map(|c| (c >> 4) as usize)
            .sum::<usize>();
        if readlen >= self.nbases {
            self.nbases += 2 * (1 + readlen - self.nbases);
            self.lens.resize(self.nbases + 1, [0; 3]);
            for q in &mut self.quals {
                q.resize(self.nbases, [0; 256]);
            }
            for b in &mut self.bases {
                b.resize(self.nbases, [0; 6]);
            }
            self.oriented.resize(self.nbases, [0; 4]);
            self.cycles.resize(self.nbases + 1, [0; 4]);
            self.ring.grow(readlen);
        }
        self.maxlen[0] = self.maxlen[0].max(readlen);
        if let Some(g) = group {
            self.maxlen[g + 1] = self.maxlen[g + 1].max(readlen);
        }
        if f & (4 | 0x100 | 0x800 | 0x200 | 0x400) == 0 {
            self.mapq[r.mapq as usize] += 1;
        }
        let paired_mapped = f & 1 != 0 && f & (4 | 8) == 0;
        let mut gc = 0;
        if original {
            self.lens[readlen][0] += 1;
            if let Some(g) = group {
                self.lens[readlen][g + 1] += 1;
            }
            self.length[0] += n as u64;
            self.qcfail += u64::from(f & 0x200 != 0);
            self.paired += u64::from(f & 1 != 0);
            if let Some(g) = group {
                self.reads[g] += 1;
                self.length[g + 1] += n as u64;
                for i in 0..n {
                    let cycle = if reverse { n - i - 1 } else { i };
                    let b = match r.base(i) {
                        1 => 0,
                        2 => 1,
                        4 => 2,
                        8 => 3,
                        15 => 4,
                        _ => 5,
                    };
                    self.bases[g][cycle][b] += 1;
                    if b < 4 {
                        self.oriented[cycle][if reverse { 3 - b } else { b }] += 1;
                    }
                    if b == 1 || b == 2 {
                        gc += 1;
                    }
                    let q = r.quals[i] as usize;
                    self.maxq = self.maxq.max(q);
                    self.quals[g][cycle][q] += 1;
                    self.sumq += q as f64;
                }
                let lo = gc * 199 / n;
                let hi = ((gc + 1) * 199 / n).min(199);
                for bin in lo..hi {
                    self.gc[g][bin] += 1;
                }
                if g == 0 {
                    for i in 0..4 {
                        if let Some(b) = r.barcode[i] {
                            self.barcodes[i].add(b, r.barqual[i]);
                        }
                    }
                }
            } else {
                self.reads[2] += 1;
            }
            if f & 4 != 0 {
                self.unmapped += 1;
            } else {
                self.mapped += 1;
                self.mapped_bases += n as u64;
                self.mq0 += u64::from(r.mapq == 0);
                if paired_mapped {
                    self.paired_mapped += 1;
                    self.proper += u64::from(f & 2 != 0);
                    self.anomalous += u64::from(r.tid != r.mate_tid);
                }
            }
        }
        if f & 4 != 0 {
            return Ok(());
        }
        let mut cycle = 0i64;
        for &c in r.cigar {
            let op = c & 15;
            let len = (c >> 4) as i64;
            if len == 0 {
                continue;
            }
            match op {
                1 => {
                    let idx = if reverse {
                        n as i64 - cycle - len
                    } else {
                        cycle
                    };
                    if idx < 0 || idx as usize >= self.nbases {
                        return Err(invalid("insertion cycle outside read"));
                    }
                    if let Some(g) = group {
                        self.cycles[idx as usize][g] += 1;
                    }
                    if len <= 300 {
                        self.indels[len as usize - 1][0] += 1;
                    }
                    cycle += len;
                }
                2 => {
                    let idx = if reverse {
                        n as i64 - cycle - 1
                    } else {
                        cycle - 1
                    };
                    if idx >= 0 {
                        if idx as usize >= self.nbases {
                            return Err(invalid("deletion cycle outside read"));
                        }
                        if let Some(g) = group {
                            self.cycles[idx as usize][g + 2] += 1;
                        }
                        if len <= 300 {
                            self.indels[len as usize - 1][1] += 1;
                        }
                    }
                }
                3 | 5 | 6 => {}
                _ => cycle += len,
            }
            if matches!(op, 0 | 1 | 7 | 8) {
                self.cigar_bases += len as u64;
            }
        }
        if paired_mapped && original {
            let insert = r.tlen.abs().min(8000) as usize;
            if insert > 0 || r.tid == r.mate_tid {
                let dir = if reverse { -1 } else { 1 };
                let mate_dir = if f & 32 != 0 { -1 } else { 1 };
                let first = if f & 64 != 0 { 1 } else { -1 };
                let delta = r.mate_pos - r.pos;
                let orientation = if dir * mate_dir > 0 {
                    2
                } else if first * delta > 0 {
                    usize::from(first * dir <= 0)
                } else if first * delta < 0 {
                    usize::from(first * dir > 0)
                } else {
                    0
                };
                self.isize[insert][orientation] += 1;
            }
        }
        self.nm = self.nm.wrapping_add(r.nm.unwrap_or(0) as u64);
        if r.cigar.is_empty() {
            return Err(invalid("mapped read with no CIGAR"));
        }
        if self.tid == r.tid && r.pos < self.pos {
            self.sorted = false;
        }
        self.pos = r.pos;
        if self.sorted {
            if self.tid == -1 || self.tid != r.tid {
                self.ring.flush(-1)?;
            }
            if self.gcd_pos == -1 || self.tid != r.tid || r.pos - self.gcd_pos > 20000 {
                self.tid = r.tid;
                self.gcd_pos = r.pos;
                self.gcd.push((0., 0));
            }
            let bin = self.gcd.last_mut().unwrap();
            bin.1 = bin.1.wrapping_add(1);
            bin.0 += gc as f32 / n as f32;
            self.ring.flush(r.pos)?;
            let mut p = r.pos;
            for &c in r.cigar {
                let op = c & 15;
                let len = (c >> 4) as i64;
                if matches!(op, 0 | 7 | 8) {
                    self.ring.insert(p, p + len)?;
                }
                if matches!(op, 0 | 2 | 3 | 7 | 8) {
                    p += len;
                }
            }
        }
        Ok(())
    }
}

pub(super) fn percentile(bins: &[(f32, u32)], p: u32) -> f32 {
    let n = p as f32 * (bins.len() + 1) as f32 / 100.;
    let k = n as usize;
    if k == 0 {
        return bins[0].1 as f32;
    }
    if k >= bins.len() {
        return bins[bins.len() - 1].1 as f32;
    }
    bins[k - 1].1 as f32 + (n - k as f32) * bins[k].1.wrapping_sub(bins[k - 1].1) as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    fn record(flag: u16) -> Record<'static> {
        Record {
            name: b"r",
            packed: &[0x12, 0x48],
            quals: &[10, 20, 30, 40],
            cigar: &[64],
            tid: 0,
            pos: 10,
            mate_tid: 0,
            mate_pos: 20,
            tlen: 14,
            flag,
            mapq: 60,
            nm: Some(2),
            barcode: [None; 4],
            barqual: [None; 4],
        }
    }
    #[test]
    fn read_categories_follow_counter_boundaries() {
        let mut s = Stats::new();
        for flag in [0, 0x100, 0x800, 0x200, 0x400, 4] {
            s.add(&record(flag)).unwrap();
        }
        assert_eq!(s.reads, [4, 0, 0]);
        assert_eq!(s.secondary, 1);
        assert_eq!(s.supplementary, 1);
        assert_eq!(s.length[0], 16);
        assert_eq!(s.mapped_bases, 12);
        assert_eq!(s.cigar_bases, 16);
        assert_eq!(s.nm, 8);
        assert_eq!(s.mapq[60], 1);
        assert_eq!(s.qcfail, 1);
        assert_eq!(s.dup, 1);
        assert_eq!(s.unmapped, 1);
        assert_eq!(s.chk[0], crc32(b"r").wrapping_mul(6));
        assert_eq!(s.bases[0][0], [4, 0, 0, 0, 0, 0]);
    }
    #[test]
    fn single_end_reverse_missing_qualities_and_ambiguity() {
        let mut r = record(16 | 4);
        r.packed = &[0x1f, 0x30];
        r.quals = &[255; 4];
        let mut s = Stats::new();
        s.add(&r).unwrap();
        assert_eq!(s.quals[0][0][255], 1);
        assert_eq!(s.maxq, 255);
        assert_eq!(s.bases[0][2][4], 1);
        assert_eq!(s.bases[0][1][5], 1);
        assert_eq!(s.oriented[3][3], 1);
        assert!(
            s.render("x.bam", "r.fa")
                .unwrap()
                .contains("SN\taverage quality:\t255.0\n")
        );
    }
    #[test]
    fn hard_and_soft_clips_indels_and_long_read_growth() {
        let packed = vec![0x11; 200];
        let qual = vec![30; 400];
        let cigar = [
            (10 << 4) | 5,
            (5 << 4) | 4,
            100 << 4,
            (2 << 4) | 1,
            100 << 4,
            (3 << 4) | 2,
            193 << 4,
        ];
        let mut r = record(1 | 64);
        r.packed = &packed;
        r.quals = &qual;
        r.cigar = &cigar;
        let mut s = Stats::new();
        s.add(&r).unwrap();
        assert_eq!(s.maxlen, [410, 410, 0]);
        assert_eq!(s.lens[410], [1, 1, 0]);
        assert_eq!(s.indels[1], [1, 0]);
        assert_eq!(s.indels[2], [0, 1]);
        assert_eq!(s.cycles[105], [1, 0, 0, 0]);
        assert_eq!(s.cycles[206], [0, 0, 1, 0]);
        assert_eq!(s.cigar_bases, 395);
        assert!(s.render("x.bam", "r.fa").unwrap().contains("FRL\t410\t1\n"));
    }
    #[test]
    fn gcd_keeps_supplementary_zero_gc_and_contig_gap_rules() {
        let mut s = Stats::new();
        let mut r = record(0);
        s.add(&r).unwrap();
        r.pos = 20010;
        s.add(&r).unwrap();
        assert_eq!(s.gcd.len(), 2); // strict >20000, not >=
        r.pos += 1;
        r.flag = 0x800;
        s.add(&r).unwrap();
        assert_eq!(s.gcd.len(), 3);
        assert_eq!(s.gcd[2], (0., 1));
        r.tid = 2;
        r.pos = 0;
        s.add(&r).unwrap();
        assert_eq!(s.gcd.len(), 4);
        assert!(s.sorted);
        r.pos = -1;
        s.add(&r).unwrap();
        assert!(!s.sorted);
        let text = s.render("x", "y").unwrap();
        assert!(!text.contains("# GC-depth."));
    }
    #[test]
    fn gcd_final_bin_is_unscaled_and_output_omits_the_largest_sorted_bin() {
        let mut s = Stats::new();
        let mut r = record(0);
        s.add(&r).unwrap();
        r.pos += 20001;
        s.add(&r).unwrap();
        let text = s.render("x", "y").unwrap();
        assert!(text.contains("GCD\t0.0\t66.667\t0.000\t0.000\t0.000\t0.000\t0.000\n"));
        assert!(text.contains("GCD\t0.5\t100.000\t0.000\t0.000\t0.000\t0.000\t0.000\n"));
        assert!(!text.contains("GCD\t50.0"));
    }
    #[test]
    fn barcode_separator_and_quality_rules() {
        let mut s = Stats::new();
        let mut r = record(0);
        r.barcode[0] = Some(b"AC-GT");
        r.barqual[0] = Some(b"II!II");
        s.add(&r).unwrap();
        let text = s.render("x", "y").unwrap();
        assert!(text.contains("BCC2\t1\t0.00\t0.00\t100.00\t0.00\t0.00\n"));
        assert!(!text.contains("BCC1\t3\t"));
        assert!(text.contains("# Barcode Qualities."));
    }
}
