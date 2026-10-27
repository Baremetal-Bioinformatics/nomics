//! dbSNP `ID`/`DB` (spec §8.12, GATK `VariantOverlapAnnotator`): dbSNP
//! records starting at the call's start, split into biallelic events
//! (`splitVariantContextToEvents(vc, trimLeft = true, …)`), matched against
//! the call's events. Gate: golden `ID`/`DB` of every call (`tests/dbsnp.rs`).

use std::io;
use std::path::Path;

use nomics_io::tabix::Tabix;

const SPAN_DEL: &[u8] = b"*";

/// `Event(contig, start, ref, alt)`: the constructor's common-suffix trim
/// when neither allele is a single base.
fn event(start: i64, reference: &[u8], alternate: &[u8]) -> (i64, Vec<u8>, Vec<u8>) {
    let (mut r, mut a) = (reference.to_vec(), alternate.to_vec());
    if r.len() != 1 && a.len() != 1 {
        while r.len() > 0 && a.len() > 0 && r.last() == a.last() {
            r.pop();
            a.pop();
        }
    }
    (start, r, a)
}

fn symbolic(a: &[u8]) -> bool {
    a.first() == Some(&b'<') || a.contains(&b'[') || a.contains(&b']') || a == b"."
}

/// `GATKVariantContextUtils.trimAlleles(vc, true, true)`.
fn trim_alleles(start: i64, alleles: &[&[u8]]) -> (i64, Vec<Vec<u8>>) {
    let unchanged = || (start, alleles.iter().map(|a| a.to_vec()).collect());
    if alleles.len() <= 1 || alleles.iter().any(|a| a.len() == 1 && *a != SPAN_DEL) {
        return unchanged();
    }
    let seqs: Vec<&[u8]> = alleles
        .iter()
        .copied()
        .filter(|a| !symbolic(a) && *a != SPAN_DEL)
        .collect();
    // AlignmentUtils.normalizeAlleles(seqs, [0, len), maxShift 0, trim true).
    let mut bounds: Vec<(usize, usize)> = seqs.iter().map(|s| (0, s.len())).collect();
    let mut min_size = bounds.iter().map(|b| b.1 - b.0).min().unwrap();
    let (mut end_trim, mut start_trim) = (0usize, 0usize);
    while min_size > 0
        && seqs
            .iter()
            .zip(&bounds)
            .all(|(s, b)| s[b.1 - 1] == seqs[0][bounds[0].1 - 1])
    {
        bounds.iter_mut().for_each(|b| b.1 -= 1);
        min_size -= 1;
        end_trim += 1;
    }
    while min_size > 0
        && seqs
            .iter()
            .zip(&bounds)
            .all(|(s, b)| s[b.0] == seqs[0][bounds[0].0])
    {
        bounds.iter_mut().for_each(|b| b.0 += 1);
        min_size -= 1;
        start_trim += 1;
    }
    let empty = bounds.iter().any(|b| b.0 == b.1);
    let end_clip = if empty && start_trim == 0 {
        end_trim - 1
    } else {
        end_trim
    };
    let start_clip = if empty && start_trim > 0 {
        start_trim - 1
    } else {
        start_trim
    };
    if start_clip == 0 && end_clip == 0 {
        return unchanged();
    }
    let trimmed = alleles
        .iter()
        .map(|a| {
            if symbolic(a) || *a == SPAN_DEL {
                a.to_vec()
            } else {
                a[start_clip..a.len() - end_clip].to_vec()
            }
        })
        .collect();
    (start + start_clip as i64, trimmed)
}

/// Biallelic events of a site: a biallelic site as is (no trimming), a
/// multiallelic one split per alt and trimmed both ways.
fn events(start: i64, alleles: &[&[u8]]) -> Vec<(i64, Vec<u8>, Vec<u8>)> {
    let alts: Vec<&[u8]> = alleles[1..]
        .iter()
        .copied()
        .filter(|a| *a != b".")
        .collect();
    match alts.len() {
        0 => Vec::new(),
        1 => vec![event(start, alleles[0], alts[0])],
        _ => alts
            .iter()
            .map(|alt| {
                let (s, t) = trim_alleles(start, &[alleles[0], alt]);
                event(s, &t[0], &t[1])
            })
            .collect(),
    }
}

pub struct Dbsnp {
    index: Tabix,
}

/// dbSNP records starting inside one span, prefetched off the ordered path.
#[derive(Default)]
pub struct DbsnpSlice {
    lines: Vec<(i64, Vec<u8>)>,
}

impl Dbsnp {
    pub fn open(vcf_gz: &Path) -> io::Result<Self> {
        Ok(Self {
            index: Tabix::open(vcf_gz)?,
        })
    }

    /// Records starting in `[beg, end]` (1-based) of `contig`.
    pub fn slice(&self, contig: &[u8], beg: i64, end: i64) -> io::Result<DbsnpSlice> {
        let mut lines = Vec::new();
        for line in self
            .index
            .lines_starting_in(contig, beg as u64, end as u64)?
        {
            let pos = line
                .split(|&b| b == b'\t')
                .nth(1)
                .and_then(|p| std::str::from_utf8(p).ok()?.parse().ok())
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "bad dbSNP POS"))?;
            lines.push((pos, line));
        }
        Ok(DbsnpSlice { lines })
    }

    /// The `ID` to write for a call (`None` = `.`, and no `DB` flag).
    pub fn id_for(
        &self,
        contig: &[u8],
        start: i64,
        alleles: &[&[u8]],
    ) -> io::Result<Option<String>> {
        Ok(self.slice(contig, start, start)?.id_for(start, alleles))
    }
}

impl DbsnpSlice {
    /// The `ID` for a call at `start` (which must lie in the slice's span).
    pub fn id_for(&self, start: i64, alleles: &[&[u8]]) -> Option<String> {
        let call_events = events(start, alleles);
        let mut ids: Vec<String> = Vec::new();
        for (_, line) in self.lines.iter().filter(|(p, _)| *p == start) {
            let f: Vec<&[u8]> = line.split(|&b| b == b'\t').take(7).collect();
            if f.len() < 7 {
                continue;
            }
            let filtered = !(f[6] == b"." || f[6] == b"PASS");
            if filtered {
                continue;
            }
            let mut site: Vec<&[u8]> = vec![f[3]];
            site.extend(f[4].split(|&b| b == b','));
            if events(start, &site).iter().any(|e| call_events.contains(e)) {
                let id = String::from_utf8_lossy(f[2]).into_owned();
                if id != "." && !ids.contains(&id) {
                    ids.push(id);
                }
            }
        }
        (!ids.is_empty()).then(|| ids.join(";"))
    }
}
