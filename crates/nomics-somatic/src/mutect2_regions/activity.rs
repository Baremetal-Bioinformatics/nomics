//! Mutect2Engine.PileupQualBuffer and LocusIteratorByState evidence.
use super::{Read, numerics::log_likelihood_ratio};

/// (current CIGAR index, read offset, reference offset within element).
pub(crate) fn element(r: &Read, pos: i64) -> Option<(usize, usize, i64)> {
    let a = &r.alignment;
    let (mut p, mut off) = (a.start, 0usize);
    for (i, &c) in a.cigar.iter().enumerate() {
        let (op, n) = (c & 15, i64::from(c >> 4));
        if matches!(op, 0 | 2 | 3 | 7 | 8) {
            if pos < p + n {
                return (op != 3).then_some((
                    i,
                    if op == 2 {
                        off.saturating_sub(1)
                    } else {
                        off + (pos - p) as usize
                    },
                    pos - p,
                ));
            }
            p += n;
        }
        if matches!(op, 0 | 1 | 4 | 7 | 8) {
            off += n as usize;
        }
    }
    None
}
pub fn covered(r: &Read, pos: i64) -> bool {
    if r.alignment.tlen <= 100
        && let Some((rev, boundary)) = super::reads::adaptor(&r.alignment)
        && if rev {
            pos <= boundary
        } else {
            pos >= boundary
        }
    {
        return false;
    }
    element(r, pos).is_some()
}

/// Reusable buffers; preserve pileup order within each sample/category.
#[derive(Default)]
pub struct Evidence {
    quals: [Vec<u8>; 6],
    pub depth: usize,
}
impl Evidence {
    pub fn clear(&mut self) {
        self.depth = 0;
        for q in &mut self.quals {
            q.clear();
        }
    }
    pub fn add(&mut self, r: &Read, pos: i64, ref_base: u8) {
        if !covered(r, pos) {
            return;
        }
        let Some((i, off, within)) = element(r, pos) else {
            return;
        };
        self.depth += 1;
        let a = &r.alignment;
        let c = a.cigar[i];
        let op = c & 15;
        let at_end = within == i64::from(c >> 4) - 1;
        let next = a.cigar.get(i + 1).copied();
        let indel = if op == 2 {
            c >> 4
        } else if at_end {
            let next_on_genome = a.cigar[i + 1..]
                .iter()
                .copied()
                .find(|c| matches!(c & 15, 0 | 2 | 7 | 8));
            if let Some(deletion) = next_on_genome.filter(|c| c & 15 == 2) {
                deletion >> 4
            } else {
                next.filter(|c| c & 15 == 1).map_or(0, |c| c >> 4)
            }
        } else {
            0
        };
        if indel > 0 {
            self.quals[5].push((30 + 10 * (indel - 1)).min(127) as u8);
            return;
        }
        let q = a.quals[off];
        let before_soft = at_end && next.is_some_and(|c| c & 15 == 4);
        let after_soft = within == 0 && i > 0 && a.cigar[i - 1] & 15 == 4;
        if (q as i8) > 6
            && ((before_soft && (a.quals[off + 1] as i8) > 6)
                || (after_soft && (a.quals[off - 1] as i8) > 6))
        {
            self.quals[5].push(30);
        } else if a.bases[off] != ref_base && q > 6 {
            let mate = if a.flag & 3 != 3 || a.flag & 8 != 0 {
                i64::from(i32::MAX)
            } else {
                a.mate_start
            };
            // Java int overflow in mateStart + read.length is intentional.
            let mate_end = (mate as i32).wrapping_add(a.bases.len() as i32) as i64;
            let q = if mate <= pos && pos < mate_end {
                q.min(20)
            } else {
                q
            };
            let index = match a.bases[off] {
                b'A' | b'a' | b'*' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            };
            self.quals[index].push(if index == 4 {
                q
            } else {
                (u16::from(q) + 5).min(254) as u8
            });
        }
    }
    fn best(&self) -> (usize, &[u8], i64) {
        let (mut index, mut sum) = (0, 0);
        for (i, qs) in self.quals.iter().enumerate() {
            // GATK sums signed bytes, even after capping at MAX_QUAL=254.
            let s = qs.iter().map(|&q| i64::from(q as i8)).sum();
            if s > sum {
                index = i;
                sum = s;
            }
        }
        (index, &self.quals[index], sum)
    }
}
pub fn is_active(tumour: &Evidence, normal: &Evidence, pon: bool, callable: &mut u64) -> bool {
    let depth = tumour.depth + normal.depth;
    if depth == 0 {
        return false;
    }
    if depth >= 10 {
        *callable += 1;
    }
    let (t, quals, _) = tumour.best();
    if log_likelihood_ratio(tumour.depth - quals.len(), quals)
        < 2.0 * crate::orientation::jmath::log(10.0)
    {
        return false;
    }
    let (n, nq, sum) = normal.best();
    if n == t && nq.len() as f64 > normal.depth as f64 * 0.3 && sum > 100 {
        return false;
    }
    if pon {
        return false;
    }
    if depth < 10 {
        *callable += 1;
    }
    true
}
