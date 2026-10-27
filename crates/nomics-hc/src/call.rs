//! GATK 4.6.2.0 `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`.
//! Gate: `tests/call.rs` over every tier0 genotyping region and call.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::finalize::HcRead;
use crate::genotype::{Call, GenotypingEngine, genotype_likelihoods, pls};
use crate::trim::{Event, TrimHap, events_for_haplotype};

const NON_REF: &[u8] = b"<NON_REF>";
const SPAN_DEL: &[u8] = b"*";
const INFORMATIVE_READ_OVERLAP_MARGIN: i64 = 2;
const MAX_PRACTICAL_ALLELE_COUNT: usize = 44;

/// Identifies a row without copying the region's reads into each event.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EvidenceIndex {
    Read(usize),
    Filtered(usize),
}

/// One row of the allele-likelihood matrix handed to annotations.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnotationRow {
    pub evidence: EvidenceIndex,
    pub likelihoods: Vec<f64>,
}

/// A called event at the boundary immediately before annotations.
#[derive(Clone, Debug, PartialEq)]
pub struct CalledEvent {
    pub start: i64,
    pub end: i64,
    /// Alleles of the genotyped call, before annotation and reverse trimming.
    pub call_alleles: Vec<Vec<u8>>,
    pub call: Call,
    /// Haplotype indices mapped to any emitted allele, for physical phasing.
    pub called_haplotypes: Vec<usize>,
    /// Allele order of `annotation_rows` (the merged, not subsetted, alleles).
    pub annotation_alleles: Vec<Vec<u8>>,
    pub annotation_rows: Vec<AnnotationRow>,
    /// Position and alleles after `makeAnnotatedCall`'s reverse trimming.
    pub final_start: i64,
    pub final_end: i64,
    pub final_alleles: Vec<Vec<u8>>,
}

/// Inputs at the `assignGenotypeLikelihoods` boundary for one region.
pub struct CallRegionInput<'a> {
    pub contig: usize,
    pub active_span: (i64, i64),
    /// Start of `fullReferenceWithPadding`, not the trimmed padded span.
    pub padded_reference_start: i64,
    pub padded_reference: &'a [u8],
    pub haplotypes: &'a [TrimHap],
    /// GATK haplotype scores, in haplotype order. They matter only if more
    /// than 44 merged alleles must be reduced.
    pub haplotype_scores: &'a [f64],
    /// Realigned likelihood evidence, in matrix-column order.
    pub reads: &'a [HcRead],
    /// Haplotype-major likelihood matrix (`[haplotype][read]`).
    pub likelihoods: &'a [Vec<f64>],
    /// Reads removed after region trimming, in their original list order.
    pub filtered_reads: &'a [HcRead],
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Variant {
    event: Event,
}

fn overlapping_events(events: &[Event], loc: i64) -> Vec<&Event> {
    let mut result: Vec<_> = events
        .iter()
        .filter(|event| event.start <= loc && loc <= event.end())
        .collect();
    if result.iter().any(|event| event.is_insertion()) {
        result.retain(|event| !(event.is_deletion() && event.end() == loc));
    }
    result
}

fn variants_at(events: &[Vec<Event>], loc: i64) -> Vec<Variant> {
    let mut result = Vec::new();
    for hap_events in events {
        for event in overlapping_events(hap_events, loc) {
            if !result
                .iter()
                .any(|variant: &Variant| variant.event == *event)
            {
                result.push(Variant {
                    event: event.clone(),
                });
            }
        }
    }
    result
}

fn merge_variants(variants: &[Variant], loc: i64, reference_base: u8) -> Vec<Vec<u8>> {
    if variants.is_empty() {
        return Vec::new();
    }
    let single_reference = [reference_base];
    let longest_reference = variants
        .iter()
        .filter(|variant| variant.event.start == loc)
        .map(|variant| variant.event.reference.as_slice())
        .max_by_key(|reference| reference.len())
        .unwrap_or(&single_reference);
    let mut alleles = vec![longest_reference.to_vec()];
    for variant in variants {
        let alternate = if variant.event.start < loc {
            SPAN_DEL.to_vec()
        } else if variant.event.reference.len() < longest_reference.len() {
            let mut padded = variant.event.alternate.clone();
            padded.extend_from_slice(&longest_reference[variant.event.reference.len()..]);
            padded
        } else {
            variant.event.alternate.clone()
        };
        if !alleles.contains(&alternate) {
            alleles.push(alternate);
        }
    }
    alleles
}

fn allele_mapper(alleles: &[Vec<u8>], loc: i64, events: &[Vec<Event>]) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); alleles.len()];
    for (haplotype, hap_events) in events.iter().enumerate() {
        let overlapping = overlapping_events(hap_events, loc);
        if overlapping.is_empty() {
            result[0].push(haplotype);
            continue;
        }
        for event in overlapping {
            let allele = if event.start < loc {
                SPAN_DEL.to_vec()
            } else if event.reference.len() == alleles[0].len() {
                event.alternate.clone()
            } else if event.reference.len() < alleles[0].len() {
                let mut padded = event.alternate.clone();
                padded.extend_from_slice(&alleles[0][event.reference.len()..]);
                padded
            } else {
                continue;
            };
            if let Some(index) = alleles.iter().position(|candidate| *candidate == allele) {
                result[index].push(haplotype);
            }
            if event.start < loc {
                break;
            }
        }
    }
    result
}

fn score_cmp(
    first: &(usize, f64, f64),
    second: &(usize, f64, f64),
    alleles: &[Vec<u8>],
) -> Ordering {
    second
        .1
        .total_cmp(&first.1)
        .then_with(|| second.2.total_cmp(&first.2))
        .then_with(|| alleles[first.0].cmp(&alleles[second.0]))
}

fn reduce_alleles(
    alleles: &mut Vec<Vec<u8>>,
    mapper: &mut Vec<Vec<usize>>,
    haplotype_scores: &[f64],
) {
    if alleles.len() <= MAX_PRACTICAL_ALLELE_COUNT {
        return;
    }
    let mut ranked = Vec::with_capacity(alleles.len() - 1);
    for allele in 1..alleles.len() {
        let mut scores: Vec<_> = mapper[allele]
            .iter()
            .map(|&haplotype| haplotype_scores[haplotype])
            .collect();
        scores.sort_by(f64::total_cmp);
        ranked.push((
            allele,
            scores.last().copied().unwrap_or(f64::NEG_INFINITY),
            scores
                .get(scores.len().wrapping_sub(2))
                .copied()
                .unwrap_or(f64::NEG_INFINITY),
        ));
    }
    ranked.sort_by(|first, second| score_cmp(first, second, alleles));
    let keep: BTreeSet<_> = std::iter::once(0)
        .chain(
            ranked
                .into_iter()
                .take(MAX_PRACTICAL_ALLELE_COUNT - 1)
                .map(|entry| entry.0),
        )
        .collect();
    let mut old = 0_usize;
    alleles.retain(|_| {
        let retained = keep.contains(&old);
        old += 1;
        retained
    });
    old = 0;
    mapper.retain(|_| {
        let retained = keep.contains(&old);
        old += 1;
        retained
    });
}

fn marginalize(
    mapper: &[Vec<usize>],
    likelihoods: &[Vec<f64>],
    read_count: usize,
) -> Vec<Vec<f64>> {
    mapper
        .iter()
        .map(|haplotypes| {
            (0..read_count)
                .map(|read| {
                    haplotypes
                        .iter()
                        .fold(f64::NEG_INFINITY, |best, &haplotype| {
                            best.max(likelihoods[haplotype][read])
                        })
                })
                .collect()
        })
        .collect()
}

fn overlaps(read: &HcRead, start: i64, end: i64) -> bool {
    read.start <= end && read.end() >= start
}

fn reverse_trim(alleles: &[Vec<u8>]) -> Vec<Vec<u8>> {
    if alleles.len() <= 1
        || alleles.iter().any(|allele| {
            allele.as_slice() != SPAN_DEL && allele.as_slice() != NON_REF && allele.len() == 1
        })
    {
        return alleles.to_vec();
    }
    let concrete: Vec<_> = alleles
        .iter()
        .filter(|allele| allele.as_slice() != SPAN_DEL && allele.as_slice() != NON_REF)
        .collect();
    let minimum = concrete
        .iter()
        .map(|allele| allele.len())
        .min()
        .unwrap_or(1);
    let mut common = 0_usize;
    while common < minimum
        && concrete.iter().all(|allele| {
            allele[allele.len() - common - 1] == concrete[0][concrete[0].len() - common - 1]
        })
    {
        common += 1;
    }
    let trim = common.min(minimum - 1);
    alleles
        .iter()
        .map(|allele| {
            if allele.as_slice() == SPAN_DEL || allele.as_slice() == NON_REF {
                allele.clone()
            } else {
                allele[..allele.len() - trim].to_vec()
            }
        })
        .collect()
}

/// `AlleleLikelihoods.updateNonRefAlleleLikelihoods`: median of likelihoods
/// strictly below the best, or the best itself when every concrete allele ties.
pub fn add_non_ref_likelihoods(matrix: &mut Vec<Vec<f64>>) {
    let concrete = matrix.len();
    matrix.push(vec![f64::NEG_INFINITY; matrix[0].len()]);
    update_non_ref_likelihoods(matrix, &(0..concrete).collect::<Vec<_>>());
}

fn update_non_ref_likelihoods(matrix: &mut [Vec<f64>], considered: &[usize]) {
    let concrete = matrix.len() - 1;
    let row = (0..matrix[0].len())
        .map(|r| {
            // searchBestAllele considers all concrete alleles, even those
            // excluded from the median by allelesToConsider.
            let best = matrix[..concrete]
                .iter()
                .map(|a| a[r])
                .fold(f64::NEG_INFINITY, f64::max);
            let mut worse: Vec<f64> = considered
                .iter()
                .filter(|&&a| a < concrete)
                .map(|&a| matrix[a][r])
                .filter(|&v| v < best)
                .collect();
            worse.sort_by(f64::total_cmp);
            match worse.len() {
                0 if concrete <= 1 => f64::NAN,
                0 => best,
                n if n % 2 == 1 => worse[n / 2],
                n => worse[n / 2 - 1] + 0.5 * (worse[n / 2] - worse[n / 2 - 1]),
            }
        })
        .collect();
    matrix[concrete] = row;
}

/// Call all assembled events whose starts lie in the region's active span.
/// Reuse `engine` across regions in contig/traversal order so spanning-deletion
/// state has the same lifetime as GATK's `GenotypingEngine`.
pub fn call_region(engine: &mut GenotypingEngine, input: CallRegionInput<'_>) -> Vec<CalledEvent> {
    assert_eq!(input.haplotypes.len(), input.likelihoods.len());
    assert_eq!(input.haplotypes.len(), input.haplotype_scores.len());
    assert!(
        input
            .likelihoods
            .iter()
            .all(|row| row.len() == input.reads.len())
    );

    let events: Vec<_> = input
        .haplotypes
        .iter()
        .map(|haplotype| {
            events_for_haplotype(
                haplotype,
                input.padded_reference,
                input.padded_reference_start,
            )
        })
        .collect();
    let event_starts: BTreeSet<_> = events
        .iter()
        .flat_map(|haplotype| haplotype.iter().map(|event| event.start))
        .collect();
    let mut result = Vec::new();

    for loc in event_starts {
        if loc < input.active_span.0 || loc > input.active_span.1 {
            continue;
        }
        let variants = variants_at(&events, loc);
        let reference_offset = usize::try_from(loc - input.padded_reference_start)
            .expect("event starts before padded reference");
        let mut alleles = merge_variants(&variants, loc, input.padded_reference[reference_offset]);
        if alleles.len() < 2 {
            continue;
        }
        let original_allele_count = alleles.len() + usize::from(engine.emit_gvcf);
        let mut mapper = allele_mapper(&alleles, loc, &events);
        reduce_alleles(&mut alleles, &mut mapper, input.haplotype_scores);
        let end = loc + alleles[0].len() as i64 - 1;
        let relevant_start = loc.saturating_sub(INFORMATIVE_READ_OVERLAP_MARGIN);
        let relevant_end = end.saturating_add(INFORMATIVE_READ_OVERLAP_MARGIN);
        let mut allele_likelihoods = marginalize(&mapper, input.likelihoods, input.reads.len());
        if engine.emit_gvcf {
            alleles.push(NON_REF.to_vec());
            add_non_ref_likelihoods(&mut allele_likelihoods);
        }
        let retained: Vec<_> = input
            .reads
            .iter()
            .enumerate()
            .filter(|(_, read)| overlaps(read, relevant_start, relevant_end))
            .map(|(index, _)| index)
            .collect();
        let genotyping_matrix: Vec<Vec<f64>> = allele_likelihoods
            .iter()
            .map(|row| retained.iter().map(|&read| row[read]).collect())
            .collect();
        let input_pls = pls(&genotype_likelihoods(&genotyping_matrix));
        let Some(call) = engine.calculate(input.contig, loc, end, &alleles, &input_pls) else {
            continue;
        };

        let call_alleles: Vec<_> = call
            .alleles_kept
            .iter()
            .map(|&allele| alleles[allele].clone())
            .collect();
        // prepareReadAlleleLikelihoodsForAnnotation updates NON_REF after
        // subsetting, before adding the zero-likelihood filtered evidence.
        if engine.emit_gvcf && call_alleles.len() != alleles.len() {
            update_non_ref_likelihoods(&mut allele_likelihoods, &call.alleles_kept);
        }
        let mut annotation_rows = Vec::with_capacity(retained.len() + input.filtered_reads.len());
        for &read in &retained {
            annotation_rows.push(AnnotationRow {
                evidence: EvidenceIndex::Read(read),
                likelihoods: allele_likelihoods.iter().map(|row| row[read]).collect(),
            });
        }
        for (filtered, read) in input.filtered_reads.iter().enumerate() {
            if overlaps(read, relevant_start, relevant_end) {
                annotation_rows.push(AnnotationRow {
                    evidence: EvidenceIndex::Filtered(filtered),
                    likelihoods: vec![0.0; alleles.len()],
                });
            }
        }
        let final_alleles = if call_alleles.len() == original_allele_count {
            call_alleles.clone()
        } else {
            reverse_trim(&call_alleles)
        };
        let final_end = loc + final_alleles[0].len() as i64 - 1;
        result.push(CalledEvent {
            start: loc,
            end,
            call_alleles,
            called_haplotypes: call
                .alleles_kept
                .iter()
                .filter_map(|&a| mapper.get(a))
                .flatten()
                .copied()
                .collect(),
            call,
            annotation_alleles: alleles,
            annotation_rows,
            final_start: loc,
            final_end,
            final_alleles,
        });
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_non_ref_median_uses_retained_alleles_and_global_best() {
        let mut matrix = vec![vec![0.0], vec![-1.0], vec![-3.0]];
        add_non_ref_likelihoods(&mut matrix);
        assert_eq!(matrix[3], [-2.0]);
        update_non_ref_likelihoods(&mut matrix, &[0, 1, 3]);
        assert_eq!(matrix[3], [-1.0]);
        // Even when the best concrete allele is dropped, it still defines
        // the strict upper bound for the qualified likelihoods.
        update_non_ref_likelihoods(&mut matrix, &[1, 2, 3]);
        assert_eq!(matrix[3], [-2.0]);
        assert_eq!(&matrix[..3], &[vec![0.0], vec![-1.0], vec![-3.0]]);
    }

    #[test]
    fn reverse_trim_leaves_one_concrete_base_and_star_unchanged() {
        assert_eq!(
            reverse_trim(&[b"CAA".to_vec(), b"TA".to_vec(), b"*".to_vec()]),
            vec![b"CA".to_vec(), b"T".to_vec(), b"*".to_vec()]
        );
    }
}
