//! GATK 4.6.2.0 `AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype`.
//! Gate: `tests/realign.rs` against `GatkLikelihoods.java` tier0 output.

use crate::assembly::cigar::{
    D, Elem, H, I, M, S, Strategy, cigar_builder, consumes_read, consumes_ref, left_align, pack,
    smith_waterman, trim_by_bases, unpack,
};
use crate::finalize::HcRead;
use crate::trim::TrimHap;

const INFORMATIVE_THRESHOLD: f64 = 0.2;

#[derive(Clone, Copy, Debug, PartialEq)]
struct BestAllele {
    haplotype: usize,
    second_haplotype: usize,
    likelihood: f64,
    second_likelihood: f64,
}

impl BestAllele {
    fn is_informative(self) -> bool {
        let confidence = if self.likelihood == self.second_likelihood {
            0.0
        } else {
            self.likelihood - self.second_likelihood
        };
        confidence > INFORMATIVE_THRESHOLD
    }
}

fn priority(haplotype: &TrimHap) -> f64 {
    let reference_term = i32::from(haplotype.is_ref);
    let cigar_term = 1_i32.wrapping_sub(
        i32::try_from(haplotype.cigar.len()).expect("haplotype CIGAR has too many elements"),
    );
    f64::from(reference_term.wrapping_add(cigar_term))
}

fn best_allele(haplotypes: &[TrimHap], likelihoods: &[Vec<f64>], read: usize) -> BestAllele {
    assert!(!haplotypes.is_empty(), "cannot choose from no haplotypes");
    assert_eq!(likelihoods.len(), haplotypes.len());

    let mut best_haplotype = 0_usize;
    let mut second_haplotype = 0_usize;
    let mut best_likelihood = likelihoods[0][read];
    let mut second_likelihood = f64::NEG_INFINITY;
    for (haplotype, values) in likelihoods.iter().enumerate().skip(1) {
        let candidate = values[read];
        if candidate > best_likelihood {
            second_haplotype = best_haplotype;
            best_haplotype = haplotype;
            second_likelihood = best_likelihood;
            best_likelihood = candidate;
        } else if candidate > second_likelihood {
            second_haplotype = haplotype;
            second_likelihood = candidate;
        }
    }

    if best_likelihood - second_likelihood < INFORMATIVE_THRESHOLD {
        let priorities: Vec<f64> = haplotypes.iter().map(priority).collect();
        let mut best_priority = priorities[best_haplotype];
        let mut second_priority = priorities[second_haplotype];
        for (haplotype, values) in likelihoods.iter().enumerate() {
            let candidate = values[read];
            if haplotype == best_haplotype || best_likelihood - candidate > INFORMATIVE_THRESHOLD {
                continue;
            }
            let candidate_priority = priorities[haplotype];
            if candidate_priority > best_priority {
                second_haplotype = best_haplotype;
                best_haplotype = haplotype;
                second_priority = best_priority;
                best_priority = candidate_priority;
            } else if candidate_priority > second_priority {
                second_haplotype = haplotype;
                second_priority = candidate_priority;
            }
        }
    }

    best_likelihood = likelihoods[best_haplotype][read];
    second_likelihood = if second_haplotype != best_haplotype {
        likelihoods[second_haplotype][read]
    } else {
        f64::NEG_INFINITY
    };
    BestAllele {
        haplotype: best_haplotype,
        second_haplotype,
        likelihood: best_likelihood,
        second_likelihood,
    }
}

fn read_start_on_reference(cigar: &[Elem], read_start_on_haplotype: usize) -> usize {
    if read_start_on_haplotype == 0 {
        return 0;
    }
    let mut reference_bases = 0_usize;
    let mut haplotype_bases = 0_usize;
    for &element in cigar {
        if consumes_ref(element.op) {
            reference_bases += element.len;
        }
        if consumes_read(element.op) {
            haplotype_bases += element.len;
        }
        if haplotype_bases >= read_start_on_haplotype {
            let excess = if consumes_ref(element.op) {
                haplotype_bases - read_start_on_haplotype
            } else {
                0
            };
            return reference_bases - excess;
        }
    }
    panic!("haplotype CIGAR does not reach the read start");
}

fn cigar_class(op: u8) -> u8 {
    match op {
        M | 7 | 8 => M,
        I | S => I,
        D => D,
        _ => panic!("unsupported CIGAR operator {op} in composition"),
    }
}

fn apply_cigar_to_cigar(first_to_second: &[Elem], second_to_third: &[Elem]) -> Vec<Elem> {
    let mut result = Vec::new();
    let (mut first_index, mut second_index) = (0_usize, 0_usize);
    let (mut first_offset, mut second_offset) = (0_usize, 0_usize);
    while first_index < first_to_second.len() && second_index < second_to_third.len() {
        let first = cigar_class(first_to_second[first_index].op);
        let second = cigar_class(second_to_third[second_index].op);
        let (output, advance_first, advance_second) = match (first, second) {
            (M, M) => (Some(M), 1, 1),
            (M, I) => (Some(I), 1, 1),
            (M, D) => (Some(D), 0, 1),
            (D, M) => (Some(D), 1, 1),
            (D, D) => (Some(D), 0, 1),
            (D, I) => (None, 1, 1),
            (I, M) | (I, D) | (I, I) => (Some(I), 1, 0),
            _ => unreachable!(),
        };
        if let Some(op) = output {
            result.push(Elem { len: 1, op });
        }
        first_offset += advance_first;
        second_offset += advance_second;
        if first_offset == first_to_second[first_index].len {
            first_index += 1;
            first_offset = 0;
        }
        if second_offset == second_to_third[second_index].len {
            second_index += 1;
            second_offset = 0;
        }
    }
    cigar_builder(result).cigar
}

fn append_original_clips(aligned: &[Elem], original: &[u32]) -> Vec<u32> {
    let original = unpack(original);
    assert!(!original.is_empty(), "read CIGAR is empty");
    let mut first = 0_usize;
    let mut last = original.len() - 1;
    let mut result = Vec::new();
    while matches!(original[first].op, S | H) && first != last {
        result.push(original[first]);
        first += 1;
    }
    result.extend_from_slice(aligned);
    let mut end = Vec::new();
    while matches!(original[last].op, S | H) && first != last {
        end.push(original[last]);
        last -= 1;
    }
    result.extend(end.into_iter().rev());
    pack(&result)
}

fn create_read_aligned_to_ref(
    original: &HcRead,
    haplotype: &TrimHap,
    reference_haplotype: &TrimHap,
    padded_reference_start: i64,
) -> HcRead {
    let (read_bases, _) = crate::assembly::cigar::hard_clip_soft_clips(
        &original.bases,
        &original.quals,
        &original.cigar,
    );
    let (raw_sw_cigar, offset) = smith_waterman(
        &haplotype.bases,
        &read_bases,
        10,
        -15,
        -30,
        -5,
        Strategy::Softclip,
    );
    let sw_cigar = cigar_builder(raw_sw_cigar).cigar;
    let offset = usize::try_from(offset).expect("Smith-Waterman returned a negative offset");

    let mut padded_haplotype_cigar = unpack(&haplotype.cigar);
    padded_haplotype_cigar.push(Elem { len: 1000, op: M });
    let padded_haplotype_cigar = cigar_builder(padded_haplotype_cigar).cigar;
    let read_start_on_reference = read_start_on_reference(&padded_haplotype_cigar, offset);
    let haplotype_read_length: usize = padded_haplotype_cigar
        .iter()
        .filter(|element| consumes_read(element.op))
        .map(|element| element.len)
        .sum();
    let haplotype_to_reference =
        trim_by_bases(&padded_haplotype_cigar, offset, haplotype_read_length - 1).cigar;
    let read_to_reference = apply_cigar_to_cigar(&sw_cigar, &haplotype_to_reference);
    let left_aligned = left_align(
        &read_to_reference,
        &reference_haplotype.bases,
        &read_bases,
        read_start_on_reference,
    );

    let mut result = original.clone();
    let reference_start = i32::try_from(padded_reference_start)
        .expect("padded reference start does not fit Java int");
    let alignment_start = i32::try_from(haplotype.alignment_start_hap_wrt_ref)
        .expect("haplotype alignment start does not fit Java int");
    let read_start = i32::try_from(read_start_on_reference)
        .expect("read start within haplotype does not fit Java int");
    let leading_deletions = i32::try_from(left_aligned.leading_deletions)
        .expect("leading deletion count does not fit Java int");
    result.start = i64::from(
        reference_start
            .wrapping_add(alignment_start)
            .wrapping_add(read_start)
            .wrapping_add(leading_deletions),
    );
    result.cigar = append_original_clips(&left_aligned.cigar, &original.cigar);
    result
}

/// Realign each retained likelihood-matrix read to its priority-adjusted best haplotype.
/// `padded_reference_start` is the trimmed genotyping region's padded-span start.
pub fn realign_reads_to_best_haplotype(
    reads: &[HcRead],
    haplotypes: &[TrimHap],
    padded_reference_start: i64,
    likelihoods: &[Vec<f64>],
) -> Vec<HcRead> {
    assert_eq!(likelihoods.len(), haplotypes.len());
    assert!(
        likelihoods.iter().all(|values| values.len() == reads.len()),
        "likelihood matrix dimensions do not match reads"
    );
    let reference_haplotype = haplotypes
        .iter()
        .find(|haplotype| haplotype.is_ref)
        .expect("missing reference haplotype");
    // AssemblyResultSet.trimTo keeps the original paddedReferenceLoc while the
    // emitted SPAN and haplotype location are trimmed.  Recover that retained
    // coordinate frame from the reference haplotype's alignment offset.
    let assembly_reference_start = i64::from(
        i32::try_from(padded_reference_start)
            .expect("padded reference start does not fit Java int")
            .wrapping_sub(
                i32::try_from(reference_haplotype.alignment_start_hap_wrt_ref)
                    .expect("reference haplotype alignment start does not fit Java int"),
            ),
    );
    reads
        .iter()
        .enumerate()
        .map(|(read_index, read)| {
            let best = best_allele(haplotypes, likelihoods, read_index);
            let _is_informative = best.is_informative();
            create_read_aligned_to_ref(
                read,
                &haplotypes[best.haplotype],
                reference_haplotype,
                assembly_reference_start,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hap(is_ref: bool, cigar_elements: usize) -> TrimHap {
        TrimHap {
            bases: vec![b'A'],
            is_ref,
            alignment_start_hap_wrt_ref: 0,
            cigar: vec![16; cigar_elements],
            contig: String::new(),
            start: 1,
            end: 1,
        }
    }

    #[test]
    fn gatk_priority_tie_break_and_informative_boundary() {
        let haplotypes = vec![hap(false, 1), hap(false, 1), hap(true, 1)];
        let likelihoods = vec![vec![0.0], vec![-0.19], vec![-0.2]];
        let best = best_allele(&haplotypes, &likelihoods, 0);
        assert_eq!(best.haplotype, 2);
        assert_eq!(best.second_haplotype, 0);
        assert!(!best.is_informative());

        let haplotypes = vec![hap(false, 1), hap(true, 1)];
        let likelihoods = vec![vec![0.0], vec![-0.2]];
        let best = best_allele(&haplotypes, &likelihoods, 0);
        assert_eq!(best.haplotype, 0);
        assert_eq!(best.second_haplotype, 1);
        assert!(!best.is_informative());
    }
}
