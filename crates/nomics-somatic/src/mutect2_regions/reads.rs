//! Query/filter/merge and post-filter palindrome clipping at GATK 4.6.2.0.
use super::{Read, Stats};
use nomics_hc::{
    finalize::{HcRead, hard_clip_bases},
    reads::Intervals,
};
use nomics_io::{
    bam::{Header, RecordView},
    bam_query::BamQuery,
};
use std::{cmp::Ordering, collections::HashMap, io, path::PathBuf};

pub const FILTERS: [&str; 12] = [
    "MappingQualityReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappingQualityNotZeroReadFilter",
    "MappedReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "NotDuplicateReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "NonChimericOriginalAlignmentReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "ReadLengthReadFilter",
    "GoodCigarReadFilter",
    "WellformedReadFilter",
];

// Local source-confirmed filters: do not change the HC path.
fn good_cigar(cigar: &[u32]) -> bool {
    let n = cigar.len();
    if n == 0 {
        return false;
    }
    let op = |i: usize| cigar[i] & 0xf;
    if cigar.iter().any(|&c| c >> 4 == 0) {
        return false;
    }
    let real = |o: u32| matches!(o, 0 | 1 | 2 | 3 | 7 | 8);
    if !cigar.iter().any(|&c| real(c & 15)) {
        return false;
    }
    // htsjdk Cigar.isValid padding and repeated-indel checks.
    for i in 0..n {
        if op(i) == 6 && i != 0 && (i == n - 1 || !real(op(i - 1)) || !real(op(i + 1))) {
            return false;
        }
        if matches!(op(i), 1 | 2) {
            for j in i + 1..n {
                if (real(op(j)) && !matches!(op(j), 1 | 2)) || op(j) == 6 {
                    break;
                }
                if op(j) == op(i) {
                    return false;
                }
            }
        }
    }
    // H only at the ends; S only at the ends or next to an end H.
    for i in 0..n {
        match op(i) {
            5 if i != 0 && i != n - 1 => return false,
            4 => {
                let left_ok = i == 0 || (i == 1 && op(0) == 5);
                let right_ok = i == n - 1 || (i == n - 2 && op(n - 1) == 5);
                if !left_ok && !right_ok {
                    return false;
                }
            }
            _ => {}
        }
    }
    let is_indel = |o: u32| o == 1 || o == 2;
    if (1..n).any(|i| is_indel(op(i)) && is_indel(op(i - 1))) {
        return false;
    }
    let core: Vec<u32> = (0..n).map(op).filter(|&o| o != 4 && o != 5).collect();
    !(core.first() == Some(&2) || core.last() == Some(&2))
}

/// `WellformedReadFilter` for mapped reads.
fn wellformed(r: &RecordView, contig_len: i64) -> bool {
    let cigar: Vec<u32> = r.cigar().collect();
    let start = r.pos() as i64 + 1;
    let end = start + r.reference_length() as i64 - 1;
    let read_len: u32 = cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 1 | 4 | 7 | 8))
        .map(|&c| c >> 4)
        .sum();
    start > 0
        && end >= start - 1
        && start <= contig_len
        && r.string_tag(b"RG").is_some()
        && r.quals().is_some_and(|q| q.len() == r.len())
        && r.len() as u32 == read_len
        && !r.is_empty()
        && !cigar.iter().any(|&c| c & 0xf == 3)
}

fn filter(r: &RecordView<'_>, len: i64) -> Option<usize> {
    let c: Vec<_> = r.cigar().collect();
    let f = r.flag();
    if r.mapq() < 20 {
        Some(0)
    } else if r.mapq() == 255 {
        Some(1)
    } else if r.mapq() == 0 {
        Some(2)
    } else if f & 4 != 0 {
        Some(3)
    } else if f & 0x100 != 0 {
        Some(4)
    } else if f & 0x400 != 0 {
        Some(5)
    } else if f & 0x200 != 0 {
        Some(6)
    } else if r
        .string_tag(b"OA")
        .zip(r.string_tag(b"XM"))
        .is_some_and(|(oa, xm)| oa.split(|&b| b == b',').next() != Some(xm))
    {
        Some(7)
    } else if r.reference_length() == 0 {
        Some(8)
    } else if r.len() < 30 {
        Some(9)
    } else if !good_cigar(&c) {
        Some(10)
    } else if !wellformed(r, len) {
        Some(11)
    } else {
        None
    }
}
fn compare(a: &RecordView<'_>, b: &RecordView<'_>) -> Ordering {
    (a.ref_id(), a.pos(), a.flag() & 16)
        .cmp(&(b.ref_id(), b.pos(), b.flag() & 16))
        .then_with(|| a.read_name().cmp(b.read_name()))
        .then_with(|| {
            (a.flag(), a.mapq(), a.mate_ref_id(), a.mate_pos(), a.tlen()).cmp(&(
                b.flag(),
                b.mapq(),
                b.mate_ref_id(),
                b.mate_pos(),
                b.tlen(),
            ))
        })
}
pub(super) fn samples(header: &Header) -> io::Result<HashMap<Vec<u8>, String>> {
    let mut out = HashMap::new();
    for line in header
        .text
        .split(|&b| b == b'\n')
        .filter(|l| l.starts_with(b"@RG\t"))
    {
        let fields: Vec<_> = line.split(|&b| b == b'\t').collect();
        let id = fields.iter().find_map(|f| f.strip_prefix(b"ID:"));
        let sm = fields.iter().find_map(|f| f.strip_prefix(b"SM:"));
        let (id, sm) = id
            .zip(sm)
            .ok_or_else(|| io::Error::other("RG requires ID and SM"))?;
        out.insert(
            id.to_vec(),
            String::from_utf8(sm.to_vec()).map_err(io::Error::other)?,
        );
    }
    Ok(out)
}
/// GATK ReadUtils.getAdaptorBoundary (no pileup adaptor-size restriction).
pub fn adaptor(r: &HcRead) -> Option<(bool, i64)> {
    let rev = r.flag & 16 != 0;
    let well = r.tlen != 0
        && r.flag & 1 != 0
        && r.flag & 12 == 0
        && rev != (r.flag & 32 != 0)
        && if rev {
            r.end() > r.mate_start
        } else {
            r.start <= r.mate_start + i64::from(r.tlen)
        };
    well.then(|| {
        (
            rev,
            if rev {
                r.mate_start - 1
            } else {
                r.start + i64::from(r.tlen).abs()
            },
        )
    })
}
fn palindrome(mut r: HcRead, reference: &[u8]) -> HcRead {
    let Some((_, boundary)) = adaptor(&r).filter(|_| r.flag & 3 == 3) else {
        return r;
    };
    let upstream = r.tlen > 0;
    let Some(&ce) = (if upstream {
        r.cigar.first()
    } else {
        r.cigar.last()
    }) else {
        return r;
    };
    if !matches!(ce & 15, 1 | 4) {
        return r;
    }
    let artifact = (ce >> 4) as usize;
    let n = (artifact + 5).min(r.bases.len());
    let (rs, re) = if upstream {
        (boundary - n as i64, boundary - 1)
    } else {
        (boundary + 1, boundary + n as i64)
    };
    if rs < 1
        || re > reference.len() as i64
        || (upstream && rs < r.start)
        || (!upstream && r.end() < re)
    {
        return r;
    }
    let last = if upstream { n - 1 } else { r.bases.len() - 1 };
    let matches = reference[rs as usize - 1..re as usize]
        .iter()
        .enumerate()
        .filter(|&(i, &b)| {
            let comp = match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            };
            comp == r.bases[last - i]
        })
        .count();
    if matches as f64 / n as f64 >= 0.9 {
        let len = r.bases.len();
        r = if upstream {
            hard_clip_bases(r, 0, artifact - 1)
        } else {
            hard_clip_bases(r, len - artifact, len - 1)
        };
    }
    r
}

pub fn load(
    inputs: &[PathBuf],
    intervals: &Intervals,
    reference: &[Option<nomics_mem::HostBuf>],
    normal: &str,
    stats: &mut Stats,
) -> io::Result<Vec<Read>> {
    let query: Vec<_> = intervals
        .padded
        .iter()
        .enumerate()
        .flat_map(|(c, ivs)| {
            ivs.iter()
                .map(move |iv| (c, (iv.start - 1) as u32, iv.end as u32))
        })
        .collect();
    let mut sources: Vec<_> = inputs
        .iter()
        .map(|p| BamQuery::open(p, &query))
        .collect::<io::Result<_>>()?;
    let maps: Vec<_> = sources
        .iter()
        .map(|s| samples(&s.header))
        .collect::<io::Result<_>>()?;
    let all_samples: std::collections::BTreeSet<_> = maps.iter().flat_map(|m| m.values()).collect();
    if all_samples.len() != 2 || !all_samples.contains(&normal.to_owned()) {
        return Err(io::Error::other(
            "only one tumour and one named normal are supported",
        ));
    }
    for s in &sources[1..] {
        if s.header.refs != sources[0].header.refs {
            return Err(io::Error::other("BAM dictionaries differ: unsupported"));
        }
    }
    let mut buffers = vec![Vec::new(); sources.len()];
    let mut present: Vec<_> = sources
        .iter_mut()
        .zip(&mut buffers)
        .map(|(s, b)| s.next_record(b))
        .collect::<io::Result<_>>()?;
    let mut reads = Vec::new();
    while let Some(k) = (0..sources.len()).filter(|&i| present[i]).min_by(|&a, &b| {
        compare(
            &RecordView::new(&buffers[a]).unwrap(),
            &RecordView::new(&buffers[b]).unwrap(),
        )
    }) {
        let r = RecordView::new(&buffers[k])?;
        if r.ref_id() >= 0 && r.pos() >= 0 {
            let c = r.ref_id() as usize;
            let start = i64::from(r.pos()) + 1;
            let end = if r.flag() & 4 != 0 {
                start
            } else {
                start + i64::from(r.reference_length()) - 1
            };
            let ivs = &intervals.padded[c];
            let at = ivs.partition_point(|iv| iv.end < start);
            if ivs.get(at).is_some_and(|iv| iv.start <= end) {
                stats.processed += 1;
                if let Some(f) = filter(&r, i64::from(sources[k].header.refs[c].len)) {
                    stats.filtered[f] += 1;
                } else {
                    let sample = maps[k]
                        .get(r.string_tag(b"RG").unwrap())
                        .ok_or_else(|| io::Error::other("unknown read group"))?;
                    let hc = HcRead::from_record(&r, reads.len() as u32);
                    let old_len = hc.bases.len();
                    let hc = palindrome(hc, reference[c].as_ref().unwrap().as_slice());
                    let clipped = old_len - hc.bases.len();
                    let left = if r.tlen() > 0 { clipped } else { 0 };
                    let indel = |tag| -> io::Result<Vec<u8>> {
                        match r.string_tag(tag) {
                            Some(q) if q.len() == old_len => Ok(q[left..left + hc.bases.len()]
                                .iter()
                                .map(|b| b - 33)
                                .collect()),
                            Some(_) => Err(io::Error::other("invalid BI/BD length")),
                            None => Ok(vec![45; hc.bases.len()]),
                        }
                    };
                    let insertion_quals = indel(b"BI")?;
                    let deletion_quals = indel(b"BD")?;
                    reads.push(Read {
                        source: (k, sources[k].virtual_offset()),
                        overlap_indels: false,
                        alignment: hc,
                        sample: sample.clone(),
                        is_normal: sample == normal,
                        insertion_quals,
                        deletion_quals,
                    });
                }
            }
        }
        present[k] = sources[k].next_record(&mut buffers[k])?;
    }
    Ok(reads)
}
