//! Diploid GATK 4.6.2.0 ReferenceConfidenceModel, including the cached indel
//! informativeness calculation and its read-length boundary quirks.
use crate::finalize::{HcRead, coordinate_cmp};
use crate::genotype::pls;
use crate::tables::{
    approximate_log10_sum_log10 as logsum, log10_one_third, qual_to_error_prob_log10,
    qual_to_prob_log10,
};
use crate::vcf::{Genotype, GvcfGenotype, Record};

#[derive(Clone, Copy)]
struct Element {
    base: u8,
    qual: u8,
    deletion: bool,
    before_indel: bool,
    offset: usize,
}

struct AlignedRead<'a> {
    read: &'a HcRead,
    pileup: Vec<Option<Element>>,
    bases: Vec<u8>,
    quals: Vec<u8>,
    informative: Option<Vec<bool>>,
    adaptor: Option<(bool, i64)>,
}

impl<'a> AlignedRead<'a> {
    fn new(read: &'a HcRead) -> Self {
        let mut pileup = Vec::new();
        let (mut bases, mut quals) = (Vec::new(), Vec::new());
        let mut r = 0;
        for (i, &c) in read.cigar.iter().enumerate() {
            let (n, op) = ((c >> 4) as usize, c & 15);
            let next = read.cigar.get(i + 1).map(|c| c & 15);
            match op {
                0 | 7 | 8 => {
                    for j in 0..n {
                        pileup.push(Some(Element {
                            base: read.bases[r],
                            qual: read.quals[r],
                            deletion: false,
                            before_indel: j + 1 == n && matches!(next, Some(1 | 2)),
                            offset: bases.len(),
                        }));
                        bases.push(read.bases[r]);
                        quals.push(read.quals[r]);
                        r += 1;
                    }
                }
                2 => {
                    for _ in 0..n {
                        pileup.push(Some(Element {
                            base: b'-',
                            qual: 30,
                            deletion: true,
                            before_indel: false,
                            offset: bases.len(),
                        }));
                        bases.push(b'-');
                        quals.push(0);
                    }
                }
                3 => {
                    for _ in 0..n {
                        pileup.push(None);
                        bases.push(b'-');
                        quals.push(0);
                    }
                }
                4 => {
                    bases.extend_from_slice(&read.bases[r..r + n]);
                    quals.extend_from_slice(&read.quals[r..r + n]);
                    r += n;
                }
                1 => r += n,
                _ => {}
            }
        }
        Self {
            read,
            pileup,
            bases,
            quals,
            informative: None,
            adaptor: adaptor_boundary(read),
        }
    }
}

// LocusIteratorByState applies ReadUtils.isBaseInsideAdaptor again after
// realignment. Its fragment-length cutoff is signed, deliberately not abs().
fn adaptor_boundary(read: &HcRead) -> Option<(bool, i64)> {
    let reverse = read.flag & 0x10 != 0;
    let well_defined = read.tlen != 0
        && read.flag & 0x1 != 0
        && read.flag & 0xc == 0
        && reverse != (read.flag & 0x20 != 0)
        && if reverse {
            read.end() > read.mate_start
        } else {
            read.start <= read.mate_start + i64::from(read.tlen)
        };
    (well_defined && read.tlen <= 100).then(|| {
        (
            reverse,
            if reverse {
                read.mate_start - 1
            } else {
                read.start + i64::from(read.tlen).abs()
            },
        )
    })
}

fn mismatch(a: u8, b: u8) -> bool {
    fn mask(b: u8) -> u8 {
        match b.to_ascii_uppercase() {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            b'R' => 5,
            b'Y' => 10,
            b'S' => 6,
            b'W' => 9,
            b'K' => 12,
            b'M' => 3,
            b'B' => 14,
            b'D' => 13,
            b'H' => 11,
            b'V' => 7,
            b'N' => 15,
            _ => 0,
        }
    }
    a != b'-' && mask(a) & mask(b) == 0
}

/// Source's bitset is initialized on the first eligible pileup element, then
/// reused at subsequent loci (including across CIGAR indels).
fn informative(
    read: &AlignedRead<'_>,
    start: usize,
    reference: &[u8],
    ref_start: usize,
) -> Vec<bool> {
    let b = &read.bases;
    let q = &read.quals;
    let mut bits = vec![false; b.len()];
    let original_len = read.read.bases.len();
    if original_len < start + 10 || reference.len() < ref_start + 10 || b.len() <= start + 10 {
        return bits;
    }
    let secondary = original_len - 10;
    let ref_shorter = b.len() > reference.len() - ref_start + start;
    let last = if !ref_shorter {
        b.len() - 10
    } else {
        reference.len() - ref_start + start - 10 + 1
    };
    let n = (b.len() - start).min(reference.len() - ref_start);
    let mut baseline = vec![0_i32; n];
    let mut sum = 0;
    for i in (0..n).rev() {
        if mismatch(b[start + i], reference[ref_start + i]) {
            sum += i32::from(q[start + i]);
        }
        baseline[i] = sum;
    }
    for size in 1..=10 {
        for insertion in [false, true] {
            let ins = if insertion { size } else { 0 };
            let del = if insertion { 0 } else { size };
            let count = (b.len() - start - ins).min(reference.len() - ref_start - del);
            let mut ro = (count + ins) as isize - 1;
            let mut fo = (count + del) as isize - 1;
            let mut cost = 0;
            while ro >= 0 && fo >= 0 {
                let (r, f) = (ro as usize, fo as usize);
                if mismatch(b[start + r], reference[ref_start + f]) {
                    cost += i32::from(q[start + r]);
                    if cost > baseline[0] {
                        break;
                    }
                }
                let site = r.min(f);
                if b[start + site] != b'-'
                    && start + site < last
                    && start + site <= secondary
                    && baseline[site] >= cost
                {
                    bits[start + site] = true;
                }
                ro -= 1;
                fo -= 1;
            }
        }
    }
    let flip = if last <= secondary {
        last
    } else {
        secondary + 1
    };
    for bit in &mut bits[..flip] {
        *bit = !*bit;
    }
    if last <= secondary && ref_shorter {
        bits[last - 1] = false;
    }
    bits
}

fn indel_gls(n: usize) -> [f64; 3] {
    let ref_lk = qual_to_prob_log10(45);
    let alt_lk = qual_to_error_prob_log10(45);
    let log2 = crate::jmath::log10(2.0);
    let n = n.min(40) as f64;
    [
        n * ref_lk,
        n * (logsum(ref_lk, alt_lk) - log2),
        n * (logsum(f64::NEG_INFINITY, alt_lk + log2) - log2),
    ]
}

/// Emit every reference-confidence locus; the consumer combines GQ bands.
/// `reference` is the padded region's reference haplotype, not the full contig.
/// One reference-confidence locus: everything a hom-ref `<NON_REF>` record
/// carries, without its allocations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RefSite {
    pub pos: i64,
    pub base: u8,
    pub ad: [i32; 2],
    pub dp: i32,
    pub gq: i32,
    pub pl: [i32; 3],
}

impl RefSite {
    pub fn to_record(&self, contig: usize) -> Record {
        Record {
            contig,
            pos: self.pos,
            id: None,
            ref_allele: (self.base as char).to_string(),
            alts: vec!["<NON_REF>".into()],
            qual: f64::NAN,
            info: Vec::new(),
            genotype: Genotype {
                alleles: vec![0, 0],
                ad: self.ad.to_vec(),
                dp: self.dp,
                gq: self.gq,
                pl: self.pl.to_vec(),
                gvcf: Some(GvcfGenotype::default()),
            },
        }
    }
}

/// What [`emit_sites`] produces, in position order.
pub enum Site {
    Ref(RefSite),
    Variant(Record),
}

/// [`emit_sites`] with every reference locus as a full [`Record`].
pub fn emit(
    contig: usize,
    span: (i64, i64),
    padded_start: i64,
    reference: &[u8],
    reads: &[HcRead],
    variants: Vec<Record>,
    mut submit: impl FnMut(Record),
) {
    emit_sites(span, padded_start, reference, reads, variants, |s| match s {
        Site::Ref(r) => submit(r.to_record(contig)),
        Site::Variant(v) => submit(v),
    });
}

/// Reference confidence over `span`, with `variants` passed through at their
/// positions.
pub fn emit_sites(
    span: (i64, i64),
    padded_start: i64,
    reference: &[u8],
    reads: &[HcRead],
    variants: Vec<Record>,
    mut submit: impl FnMut(Site),
) {
    let mut ordered: Vec<&HcRead> = reads.iter().collect();
    ordered.sort_by(|a, b| coordinate_cmp(a, b));
    let mut aligned: Vec<_> = ordered.into_iter().map(AlignedRead::new).collect();
    let mut variants = variants.into_iter().peekable();
    let log2 = crate::jmath::log10(2.0);
    let third = log10_one_third();
    for pos in span.0..=span.1 {
        if variants.peek().is_some_and(|v| v.pos == pos) {
            submit(Site::Variant(variants.next().unwrap()));
            continue;
        }
        let offset = (pos - padded_start) as usize;
        let base = reference[offset];
        let mut gl = [0.0_f64; 3];
        let mut ad = [0; 2];
        let mut n = 0;
        for read in &mut aligned {
            let r = read.read;
            if pos < r.start {
                break;
            }
            if read.adaptor.is_some_and(|(reverse, boundary)| {
                if reverse {
                    pos <= boundary
                } else {
                    pos >= boundary
                }
            }) {
                continue;
            }
            let Some(Some(Element {
                base: b,
                qual: q,
                deletion: del,
                before_indel: before,
                offset: igv,
            })) = read.pileup.get((pos - r.start) as usize).copied()
            else {
                continue;
            };
            if del || q > 6 {
                let alt = del || b != base;
                let yes = qual_to_prob_log10(q);
                let no = qual_to_error_prob_log10(q) + third;
                let (rl, al) = if alt { (no, yes) } else { (yes, no) };
                gl[0] += rl + log2;
                gl[1] += logsum(rl, al);
                gl[2] += al + log2;
                ad[usize::from(alt)] += 1;
            }
            if !del && !before && n <= 40 {
                if read.informative.is_none() {
                    read.informative = Some(informative(read, igv, reference, offset));
                }
                if read.informative.as_ref().unwrap().get(igv) == Some(&true) {
                    n += 1;
                }
            }
        }
        let dp = ad.iter().sum::<i32>();
        let denom = f64::from(dp) * log2;
        for x in &mut gl {
            *x -= denom;
        }
        gl[1] = gl[1].min(gl[0]);
        gl[2] = gl[2].min(gl[0]);
        let indel = indel_gls(n);
        // Both vectors are capped by hom-ref. GQLog10 is second-best minus ref.
        let snp_gq = gl[1].max(gl[2]) - gl[0];
        let indel_gq = indel[1].max(indel[2]) - indel[0];
        let pl = pls(if indel_gq > snp_gq { &indel } else { &gl });
        let gq = pl[1].min(pl[2]) - pl[0];
        submit(Site::Ref(RefSite {
            pos,
            base,
            ad,
            dp,
            gq,
            pl: [pl[0], pl[1], pl[2]],
        }));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pileups_exclude_adaptors_at_realigned_boundaries() {
        let reference = b"ACGTCAGTGCATCGATGCTAGCATCGATGCAT";
        let mut read = HcRead {
            name: b"r".to_vec(),
            flag: 0x21,
            contig: 0,
            start: 1,
            cigar: vec![30 << 4],
            bases: reference.to_vec(),
            quals: vec![30; 30],
            mapq: 60,
            mate_contig: 0,
            mate_start: 10,
            tlen: 10,
            record: 0,
        };
        let depths = |r: &HcRead| {
            let mut depths = Vec::new();
            emit(
                0,
                (1, 30),
                1,
                reference,
                std::slice::from_ref(r),
                vec![],
                |v| {
                    depths.push(v.genotype.dp);
                },
            );
            depths
        };
        let forward = depths(&read);
        assert_eq!(&forward[..10], &[1; 10]);
        assert_eq!(&forward[10..], &[0; 20]);
        // A reverse read's negative TLEN passes the signed >100 cutoff.
        read.flag = 0x11;
        read.tlen = -1000;
        read.mate_start = 15;
        let reverse = depths(&read);
        assert_eq!(&reverse[..14], &[0; 14]);
        assert_eq!(&reverse[14..], &[1; 16]);
        read.flag |= 0x8;
        assert_eq!(depths(&read), vec![1; 30]);
    }

    #[test]
    fn empty_pileup_and_read_end_indel_uncertainty() {
        let reference = b"ACGTCAGTGCATCGATGCTAGCATCGATGCAT";
        let read = HcRead {
            name: b"r".to_vec(),
            flag: 0,
            contig: 0,
            start: 1,
            cigar: vec![30 << 4],
            bases: reference.to_vec(),
            quals: vec![30; 30],
            mapq: 60,
            mate_contig: -1,
            mate_start: 0,
            tlen: 0,
            record: 0,
        };
        let mut records = Vec::new();
        emit(0, (1, 30), 1, reference, &[read], vec![], |r| {
            records.push(r)
        });
        // Q30 SNP hom-alt PL rounds 30 + 10*log10(3) to 35; SNP GQ
        // is slightly worse than indel GQ, though both round to 3.
        assert_eq!(records[0].genotype.pl, [0, 3, 35]);
        assert_eq!(records[0].genotype.dp, 1);
        assert_eq!(records[29].genotype.pl, [0, 0, 0]);
        assert_eq!(records[29].genotype.dp, 1);
        let mut empty = Vec::new();
        emit(0, (1, 2), 1, reference, &[], vec![], |r| empty.push(r));
        assert_eq!(empty[0].genotype.pl, [0, 0, 0]);
        assert_eq!(empty[0].genotype.dp, 0);
    }
}
