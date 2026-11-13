//! Read/fragment annotations in GATK 4.6.2.0's StandardMutectAnnotation group.
//! No rounding of AF is performed at this boundary. Matrix rows are explicitly
//! labelled with their emitted-allele indices: Java marginalization can reorder
//! rows through HashMap, even when the VCF allele order stays unchanged.
use super::{evidence, fragments::Fragment};
use nomics_hc::{finalize::HcRead, jmath::log};
use std::io;

/// First maximum, reference-priority tie breaking, and the *literal* .2
/// informative threshold (not the ln-scale tie window).
pub fn best(matrix: &[Vec<f64>], column: usize, order: &[usize]) -> (usize, f64) {
    best_with_reference(matrix, column, order, 0)
}

pub(super) fn best_with_reference(
    matrix: &[Vec<f64>],
    column: usize,
    order: &[usize],
    reference: usize,
) -> (usize, f64) {
    assert!(!order.is_empty());
    let mut b = 0;
    let mut second = 0;
    let mut bv = matrix[order[0]][column];
    let mut sv = f64::NEG_INFINITY;
    for a in 1..order.len() {
        let v = matrix[order[a]][column];
        if v > bv {
            second = b;
            b = a;
            sv = bv;
            bv = v;
        } else if v > sv {
            second = a;
            sv = v;
        }
    }
    let window = log(10.0) * 0.2;
    if bv - sv < window {
        let mut priority = i32::from(order[b] == reference);
        let mut second_priority = i32::from(order[second] == reference);
        for a in 0..order.len() {
            if a == b || bv - matrix[order[a]][column] > window {
                continue;
            }
            let candidate = i32::from(order[a] == reference);
            if candidate > priority {
                second = b;
                b = a;
                second_priority = priority;
                priority = candidate;
            } else if candidate > second_priority {
                second = a;
                second_priority = candidate;
            }
        }
    }
    let bv = matrix[order[b]][column];
    let sv = if b == second {
        f64::NEG_INFINITY
    } else {
        matrix[order[second]][column]
    };
    let confidence = if bv == sv { 0.0 } else { bv - sv };
    (order[b], confidence)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_gatk_annotation_bits() {
        let mut matrix = Vec::new();
        let mut order = Vec::new();
        for line in include_str!("annotations.tsv").lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            match f[0] {
                "E" => {
                    let a: usize = f[1].parse().unwrap();
                    let n: usize = f[2].parse().unwrap();
                    let bits: Vec<_> = f[3..]
                        .iter()
                        .map(|s| u64::from_str_radix(s, 16).unwrap())
                        .collect();
                    matrix = (0..a)
                        .map(|i| {
                            bits[i * n..(i + 1) * n]
                                .iter()
                                .map(|&b| f64::from_bits(b))
                                .collect()
                        })
                        .collect();
                    let mut actual = effective_counts(&matrix);
                    actual.extend_from_slice(&fractions(&matrix)[1..]);
                    assert_eq!(
                        actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        bits[a * n..],
                        "{a}x{n}"
                    );
                }
                "O" => {
                    let a: usize = f[1].parse().unwrap();
                    let alleles: Vec<_> = (0..a)
                        .map(|i| format!("A{}", "C".repeat(i)).into_bytes())
                        .collect();
                    order = allele_map_order(&alleles).unwrap();
                    assert_eq!(
                        order,
                        f[2..]
                            .iter()
                            .map(|s| s.parse::<usize>().unwrap())
                            .collect::<Vec<_>>()
                    );
                }
                "B" => {
                    for (column, pair) in f[2..].as_chunks::<2>().0.iter().enumerate() {
                        let (a, confidence) = best(&matrix, column, &order);
                        assert_eq!(a, pair[0].parse::<usize>().unwrap());
                        assert_eq!(
                            confidence.to_bits(),
                            u64::from_str_radix(pair[1], 16).unwrap(),
                            "best column {column}"
                        );
                    }
                }
                "M" => {
                    let n: usize = f[1].parse().unwrap();
                    let mut values: Vec<_> =
                        f[2..2 + n].iter().map(|s| s.parse().unwrap()).collect();
                    assert_eq!(median(&mut values, 17), f[2 + n].parse::<i32>().unwrap());
                }
                "R" => {
                    let mut cigar = Vec::new();
                    let mut n = 0;
                    for c in f[1].bytes() {
                        if c.is_ascii_digit() {
                            n = n * 10 + u32::from(c - b'0');
                        } else {
                            let op = b"MIDNSHP=X".iter().position(|&x| x == c).unwrap() as u32;
                            cigar.push(n << 4 | op);
                            n = 0;
                        }
                    }
                    let length = cigar
                        .iter()
                        .filter(|c| matches!(**c & 15, 0 | 1 | 4 | 7 | 8))
                        .map(|c| (c >> 4) as usize)
                        .sum();
                    let r = HcRead {
                        name: b"r".to_vec(),
                        flag: 0,
                        contig: 0,
                        start: 100,
                        cigar,
                        bases: vec![b'A'; length],
                        quals: (0..length).map(|i| (i + 15) as u8).collect(),
                        mapq: 60,
                        mate_contig: 0,
                        mate_start: 0,
                        tlen: 0,
                        record: 0,
                    };
                    let pos = f[2].parse().unwrap();
                    let expected = |s: &str| {
                        if s == "null" {
                            None
                        } else {
                            Some(s.parse::<i32>().unwrap())
                        }
                    };
                    assert_eq!(base_quality(&r, pos), expected(f[3]), "BQ {line}");
                    assert_eq!(read_position(&r, pos), expected(f[4]), "RP {line}");
                }
                _ => panic!("unknown fixture row {line}"),
            }
        }
    }
}

/// Collectors.toMap (putIfAbsent) on SimpleAllele keys. No identity-hash
/// fallback is guessed: a tree bin is an explicit unsupported condition.
pub fn allele_map_order(alleles: &[Vec<u8>]) -> io::Result<Vec<usize>> {
    let hashes: Vec<_> = alleles
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let array = a.iter().fold(1u32, |h, &b| {
                h.wrapping_mul(31).wrapping_add((b as i8 as i32) as u32)
            });
            let h = array
                .wrapping_mul(31)
                .wrapping_add(if i == 0 { 1231 } else { 1237 });
            h ^ (h >> 16)
        })
        .collect();
    let mut capacity = 16;
    while alleles.len() > capacity * 3 / 4 {
        capacity *= 2;
    }
    loop {
        let mut bins = vec![Vec::new(); capacity];
        for (i, &h) in hashes.iter().enumerate() {
            bins[h as usize & (capacity - 1)].push(i);
        }
        // putVal requests treeification on the ninth node. Below 64 it resizes.
        if bins.iter().any(|b| b.len() >= 9) {
            if capacity < 64 {
                capacity *= 2;
                continue;
            }
            return Err(io::Error::other(
                "Mutect2 allele HashMap tree-bin order is not supported",
            ));
        }
        return Ok(bins.into_iter().flatten().collect());
    }
}

/// Normalized posterior means in emitted order, including reference.
pub fn fractions(matrix: &[Vec<f64>]) -> Vec<f64> {
    let prior = vec![1.0; matrix.len()];
    let alpha = if matrix[0].is_empty() {
        prior
    } else {
        evidence::posterior(matrix, &prior)
    };
    let total = alpha.iter().fold(0.0, |a, &b| a + b);
    alpha.iter().map(|&a| a / total).collect()
}

/// Flat per-column soft counts used for initial AD (later overwritten).
pub fn effective_counts(matrix: &[Vec<f64>]) -> Vec<f64> {
    if matrix[0].is_empty() {
        return vec![0.0; matrix.len()];
    }
    let column = |i| evidence::normalize(&matrix.iter().map(|r| r[i]).collect::<Vec<_>>());
    let mut total = column(0);
    for i in 1..matrix[0].len() {
        for (sum, value) in total.iter_mut().zip(column(i)) {
            *sum += value;
        }
    }
    total
}

/// Commons Median followed by FastMath.round; widening happens before addition.
pub fn median(values: &mut [i32], empty: i32) -> i32 {
    if values.is_empty() {
        return empty;
    }
    values.sort_unstable();
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        ((f64::from(values[n / 2 - 1]) + f64::from(values[n / 2])) / 2.0 + 0.5).floor() as i32
    }
}

fn index_at(read: &HcRead, pos: i64) -> Option<(usize, u32)> {
    let soft: i64 = read
        .cigar
        .iter()
        .take_while(|c| matches!(**c & 15, 4 | 5))
        .filter(|c| **c & 15 == 4)
        .map(|c| i64::from(c >> 4))
        .sum();
    let mut reference = read.start - soft;
    let mut index = 0usize;
    for &c in &read.cigar {
        let op = c & 15;
        let n = i64::from(c >> 4);
        let cr = matches!(op, 0 | 1 | 4 | 7 | 8);
        let cf = matches!(op, 0 | 2 | 3 | 4 | 7 | 8);
        if cf && reference <= pos && pos < reference + n {
            return Some((index + if cr { (pos - reference) as usize } else { 0 }, op));
        }
        if cr {
            index += n as usize;
        }
        if cf {
            reference += n;
        }
    }
    None
}

pub fn base_quality(read: &HcRead, pos: i64) -> Option<i32> {
    if pos < read.start || pos > read.end() {
        return None;
    }
    let (i, op) = index_at(read, pos)?;
    matches!(op, 0 | 1 | 4 | 7 | 8).then(|| i32::from(read.quals[i] as i8))
}

pub fn read_position(read: &HcRead, pos: i64) -> Option<i32> {
    // ReadPosition.getPosition rejects out-of-alignment positions *before*
    // ReadPosRankSumTest's special leading-insertion case can run.
    if pos < read.start || pos > read.end() {
        return None;
    }
    let (i, _) = index_at(read, pos)?;
    let hard = |c: Option<&u32>| c.filter(|c| **c & 15 == 5).map_or(0, |c| (c >> 4) as i32);
    let left = hard(read.cigar.first()) + i as i32;
    let right = read.bases.len() as i32 - 1 - i as i32 + hard(read.cigar.last());
    Some(left.min(right))
}

#[derive(Debug)]
pub struct SampleAnnotation {
    pub af: Vec<f64>,
    pub ad: Vec<i32>,
    pub dp: i32,
    pub fad: Vec<i32>,
    pub f1r2: Vec<i32>,
    pub f2r1: Vec<i32>,
    pub sb: [i32; 4],
    pub strands: Vec<[i32; 2]>,
}

#[derive(Debug)]
pub struct ReadAnnotations {
    pub dp: usize,
    pub mbq: Vec<i32>,
    pub mmq: Vec<i32>,
    pub mfrl: Vec<i32>,
    pub mpos: Vec<i32>,
    pub strands: Vec<[i32; 2]>,
}

pub struct SampleInput<'a> {
    pub reads: &'a [HcRead],
    /// Full sample reads, indexed by Fragment::reads before overlap retention.
    pub fragment_reads: &'a [HcRead],
    pub read_matrix: &'a [Vec<f64>],
    pub fragments: &'a [&'a Fragment],
    pub fragment_matrix: &'a [Vec<f64>],
}

/// Matrices are in emitted order; `order` is trimmedToUntrimmedAlleleMap's
/// HashMap iteration order. These are different domains and never conflated.
pub fn annotate(
    pos: i64,
    order: &[usize],
    inputs: &[SampleInput<'_>],
) -> (ReadAnnotations, Vec<SampleAnnotation>) {
    let a = order.len();
    let mut bq = vec![Vec::new(); a];
    let mut mq = vec![Vec::new(); a];
    let mut fl = vec![Vec::new(); a];
    let mut rp = vec![Vec::new(); a];
    let mut total = ReadAnnotations {
        dp: 0,
        mbq: Vec::new(),
        mmq: Vec::new(),
        mfrl: Vec::new(),
        mpos: Vec::new(),
        strands: vec![[0; 2]; a],
    };
    let mut samples = Vec::new();
    for input in inputs {
        total.dp += input.reads.len();
        let mut sample = SampleAnnotation {
            af: fractions(input.fragment_matrix)[1..].to_vec(),
            ad: vec![0; a],
            dp: 0,
            fad: vec![0; a],
            f1r2: vec![0; a],
            f2r1: vec![0; a],
            sb: [0; 4],
            strands: vec![[0; 2]; a],
        };
        for (i, read) in input.reads.iter().enumerate() {
            let (allele, confidence) = best(input.read_matrix, i, order);
            if confidence <= 0.2 || confidence.is_nan() {
                continue;
            }
            sample.ad[allele] += 1;
            sample.dp += 1;
            let strand = usize::from(read.flag & 16 != 0);
            sample.sb[usize::from(allele != 0) * 2 + strand] += 1;
            sample.strands[allele][strand] += 1;
            if matches!(read.mapq, 0 | 255) {
                continue;
            }
            if let Some(q) = base_quality(read, pos) {
                bq[allele].push(q);
            }
            mq[allele].push(i32::from(read.mapq));
            fl[allele].push(read.tlen.wrapping_abs());
            if let Some(p) = read_position(read, pos) {
                rp[allele].push(p);
            }
        }
        // AS strand counts have a per-sample strict minimum of two.
        if sample.dp > 2 {
            for (sum, counts) in total.strands.iter_mut().zip(&sample.strands) {
                for (s, &c) in sum.iter_mut().zip(counts) {
                    *s += c;
                }
            }
        }
        for (i, fragment) in input.fragments.iter().enumerate() {
            let (allele, confidence) = best(input.fragment_matrix, i, order);
            if confidence <= 0.2 || confidence.is_nan() {
                continue;
            }
            sample.fad[allele] += 1;
            let read = &input.fragment_reads[fragment.reads[0]];
            if matches!(read.mapq, 0 | 255) || base_quality(read, pos).unwrap_or(0) < 20 {
                continue;
            }
            if (read.flag & 16 != 0) != (read.flag & 64 != 0) {
                sample.f1r2[allele] += 1;
            } else {
                sample.f2r1[allele] += 1;
            }
        }
        samples.push(sample);
    }
    total.mbq = bq.iter_mut().map(|v| median(v, 0)).collect();
    total.mmq = mq.iter_mut().map(|v| median(v, 60)).collect();
    total.mfrl = fl.iter_mut().map(|v| median(v, 0)).collect();
    total.mpos = rp.iter_mut().skip(1).map(|v| median(v, 50)).collect();
    (total, samples)
}
