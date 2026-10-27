//! GATK 4.6.2.0 post-assembly region trimming.
//!
//! This is the default (non-legacy) path through `HaplotypeCallerEngine.callRegion`:
//! event extraction, `AssemblyRegionTrimmer`, `AssemblyResultSet.trimTo`, read
//! clipping, and the final genotyping read filters. CIGAR elements use BAM's
//! numeric encoding (`length << 4 | operator`). Coordinates are 1-based closed.

use std::collections::BTreeMap;

use crate::finalize::{HcRead, coordinate_cmp, hard_clip_to_region};

const M: u32 = 0;
const I: u32 = 1;
const D: u32 = 2;
const N: u32 = 3;
const S: u32 = 4;
const H: u32 = 5;
const P: u32 = 6;
const EQ: u32 = 7;
const X: u32 = 8;

#[inline]
fn op(element: u32) -> u32 {
    element & 0xf
}

#[inline]
fn len(element: u32) -> usize {
    (element >> 4) as usize
}

#[inline]
fn ce(length: usize, operator: u32) -> u32 {
    (length as u32) << 4 | operator
}

#[inline]
#[cfg(test)]
fn consumes_read(operator: u32) -> bool {
    matches!(operator, M | I | S | EQ | X)
}

#[inline]
fn consumes_reference(operator: u32) -> bool {
    matches!(operator, M | D | N | EQ | X)
}

/// An assembled haplotype before or after region trimming.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TrimHap {
    pub bases: Vec<u8>,
    pub is_ref: bool,
    pub alignment_start_hap_wrt_ref: i64,
    pub cigar: Vec<u32>,
    pub contig: String,
    pub start: i64,
    pub end: i64,
}

/// The post-assembly outcome at the `callRegion` trimming boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TrimResult {
    NoVariation,
    Trimmed {
        span: (i64, i64),
        padded_span: (i64, i64),
        variation_present: bool,
        haplotypes: Vec<TrimHap>,
        reads: Vec<HcRead>,
        filtered: Vec<HcRead>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Event {
    pub(crate) start: i64,
    pub(crate) reference: Vec<u8>,
    pub(crate) alternate: Vec<u8>,
}

impl Event {
    fn new(start: i64, mut reference: Vec<u8>, mut alternate: Vec<u8>) -> Self {
        if reference.len() != 1 && alternate.len() != 1 {
            let mut common = 0;
            let limit = reference.len().min(alternate.len());
            while common < limit
                && reference[reference.len() - common - 1]
                    == alternate[alternate.len() - common - 1]
            {
                common += 1;
            }
            if common != 0 {
                reference.truncate(reference.len() - common);
                alternate.truncate(alternate.len() - common);
            }
        }
        Self {
            start,
            reference,
            alternate,
        }
    }

    pub(crate) fn end(&self) -> i64 {
        self.start + self.reference.len() as i64 - 1
    }

    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternate.len() == 1
    }

    fn is_indel(&self) -> bool {
        self.reference.len() != self.alternate.len()
    }

    pub(crate) fn is_insertion(&self) -> bool {
        self.reference.len() == 1 && self.alternate.len() > 1
    }

    pub(crate) fn is_deletion(&self) -> bool {
        self.reference.len() > 1 && self.alternate.len() == 1
    }
}

fn regular(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

fn all_regular(bases: &[u8]) -> bool {
    bases.iter().copied().all(regular)
}

fn compound(first: &Event, second: &Event) -> Event {
    assert_eq!(first.start, second.start);
    if first.is_snp() || second.is_snp() {
        assert!(!(first.is_snp() && second.is_snp()));
        let (snp, indel) = if first.is_snp() {
            (first, second)
        } else {
            (second, first)
        };
        if snp.reference == indel.reference {
            let mut alt = snp.alternate.clone();
            alt.extend_from_slice(&indel.alternate[1..]);
            Event::new(snp.start, snp.reference.clone(), alt)
        } else {
            Event::new(snp.start, indel.reference.clone(), snp.alternate.clone())
        }
    } else {
        assert!(
            (first.is_deletion() && second.is_insertion())
                || (first.is_insertion() && second.is_deletion())
        );
        let insertion = if first.is_insertion() { first } else { second };
        let deletion = if first.is_deletion() { first } else { second };
        Event::new(
            first.start,
            deletion.reference.clone(),
            insertion.alternate.clone(),
        )
    }
}

fn add_event(events: &mut BTreeMap<i64, Event>, event: Event) {
    if let Some(previous) = events.remove(&event.start) {
        events.insert(event.start, compound(&previous, &event));
    } else {
        events.insert(event.start, event);
    }
}

pub(crate) fn events_for_haplotype(
    haplotype: &TrimHap,
    reference: &[u8],
    ref_start: i64,
) -> Vec<Event> {
    if haplotype.alignment_start_hap_wrt_ref < 0 {
        return Vec::new();
    }
    let mut events = BTreeMap::new();
    let mut ref_pos = haplotype.alignment_start_hap_wrt_ref as usize;
    let mut alignment_pos = 0usize;

    for (cigar_index, &element) in haplotype.cigar.iter().enumerate() {
        let length = len(element);
        match op(element) {
            I => {
                if ref_pos > 0 && cigar_index > 0 && cigar_index + 1 < haplotype.cigar.len() {
                    let pad = reference[ref_pos - 1];
                    let mut alternate = Vec::with_capacity(length + 1);
                    alternate.push(pad);
                    alternate
                        .extend_from_slice(&haplotype.bases[alignment_pos..alignment_pos + length]);
                    if all_regular(&alternate) {
                        add_event(
                            &mut events,
                            Event::new(ref_start + ref_pos as i64 - 1, vec![pad], alternate),
                        );
                    }
                }
                alignment_pos += length;
            }
            S => alignment_pos += length,
            D => {
                if ref_pos > 0 {
                    let reference_allele = reference[ref_pos - 1..ref_pos + length].to_vec();
                    let pad = reference[ref_pos - 1];
                    if regular(pad) && all_regular(&reference_allele) {
                        add_event(
                            &mut events,
                            Event::new(ref_start + ref_pos as i64 - 1, reference_allele, vec![pad]),
                        );
                    }
                }
                ref_pos += length;
            }
            M | EQ | X => {
                for offset in 0..length {
                    let ref_base = reference[ref_pos + offset];
                    let alt_base = haplotype.bases[alignment_pos + offset];
                    if ref_base != alt_base && regular(ref_base) && regular(alt_base) {
                        add_event(
                            &mut events,
                            Event::new(
                                ref_start + (ref_pos + offset) as i64,
                                vec![ref_base],
                                vec![alt_base],
                            ),
                        );
                    }
                }
                ref_pos += length;
                alignment_pos += length;
            }
            N | H | P => panic!("unsupported haplotype CIGAR operator {}", op(element)),
            operator => panic!("invalid haplotype CIGAR operator {operator}"),
        }
    }
    events.into_values().collect()
}

fn find_repeated_substring(bases: &[u8]) -> usize {
    for period in 1..=bases.len() {
        let candidate = &bases[..period];
        let mut matches = true;
        for start in (period..bases.len()).step_by(period) {
            let end = start + period;
            if end > bases.len() || &bases[start..end] != candidate {
                matches = false;
                break;
            }
        }
        if matches {
            return period;
        }
    }
    bases.len() + 1
}

fn leading_repetitions(unit: &[u8], bases: &[u8]) -> usize {
    let mut count = 0;
    let mut start = 0;
    while start + unit.len() <= bases.len() && &bases[start..start + unit.len()] == unit {
        count += 1;
        start += unit.len();
    }
    count
}

/// Return `(reference repeats, alternate repeats, repeat-unit length)` when
/// GATK recognizes the event as a tandem repeat.
fn tandem_repeat(event: &Event, context: &[u8]) -> Option<(usize, usize, usize)> {
    if !event.is_indel() {
        return None;
    }
    let reference = &event.reference[1..];
    let alternate = &event.alternate[1..];
    let longer = if alternate.len() > reference.len() {
        alternate
    } else {
        reference
    };
    let unit_len = find_repeated_substring(longer);
    let unit = &longer[..unit_len];
    let allele_ref_repeats = leading_repetitions(unit, reference);

    let mut ref_with_context = Vec::with_capacity(reference.len() + context.len());
    ref_with_context.extend_from_slice(reference);
    ref_with_context.extend_from_slice(context);
    let mut alt_with_context = Vec::with_capacity(alternate.len() + context.len());
    alt_with_context.extend_from_slice(alternate);
    alt_with_context.extend_from_slice(context);

    let ref_count = leading_repetitions(unit, &ref_with_context) - allele_ref_repeats;
    let alt_count = leading_repetitions(unit, &alt_with_context) - allele_ref_repeats;
    (ref_count != 0 && alt_count != 0).then_some((ref_count, alt_count, unit_len))
}

#[derive(Default)]
struct CigarBuilder {
    elements: Vec<u32>,
    last: Option<u32>,
    remove_end_deletions: bool,
}

impl CigarBuilder {
    fn new(remove_end_deletions: bool) -> Self {
        Self {
            remove_end_deletions,
            ..Self::default()
        }
    }

    fn last_two_deletion_insertion(&self) -> bool {
        self.last == Some(I)
            && self.elements.len() > 1
            && op(self.elements[self.elements.len() - 2]) == D
    }

    fn add(&mut self, element: u32) {
        let length = len(element);
        if length == 0 {
            return;
        }
        let operator = op(element);
        if self.remove_end_deletions
            && operator == D
            && (self.last.is_none()
                || self.last.is_some_and(|last| matches!(last, S | H))
                || (self.last == Some(I)
                    && (self.elements.len() == 1
                        || matches!(op(self.elements[self.elements.len() - 2]), S | H))))
        {
            return;
        }
        if self.last == Some(operator) {
            let last = self.elements.len() - 1;
            self.elements[last] = ce(len(self.elements[last]) + length, operator);
        } else if operator == D && self.last == Some(I) {
            let size = self.elements.len();
            if size > 1 && op(self.elements[size - 2]) == D {
                self.elements[size - 2] = ce(len(self.elements[size - 2]) + length, D);
            } else {
                self.elements.insert(size - 1, element);
            }
        } else {
            self.elements.push(element);
            self.last = Some(operator);
        }
    }

    fn finish(mut self) -> Vec<u32> {
        if self.remove_end_deletions {
            if self.last == Some(D) {
                self.elements.pop();
            } else if self.last_two_deletion_insertion() {
                let deletion = self.elements.len() - 2;
                self.elements.remove(deletion);
            }
        }
        self.elements
    }
}

fn bases_covering_reference_interval(
    start: usize,
    end: usize,
    bases: &[u8],
    cigar: &[u32],
) -> Option<Vec<u8>> {
    let mut ref_pos = 0usize;
    let mut bases_pos = 0usize;
    let mut bases_start = None;
    let mut bases_stop = None;

    'cigar: for &element in cigar {
        match op(element) {
            I => bases_pos += len(element),
            M | X | EQ => {
                for _ in 0..len(element) {
                    if ref_pos == start {
                        bases_start = Some(bases_pos);
                    }
                    if ref_pos == end {
                        bases_stop = Some(bases_pos);
                        break 'cigar;
                    }
                    ref_pos += 1;
                    bases_pos += 1;
                }
            }
            D => {
                for _ in 0..len(element) {
                    if ref_pos == start || ref_pos == end {
                        return None;
                    }
                    ref_pos += 1;
                }
            }
            operator => panic!("unsupported trim CIGAR operator {operator}"),
        }
    }
    let from = bases_start.expect("trim start not reached by haplotype CIGAR");
    let through = bases_stop.expect("trim end not reached by haplotype CIGAR");
    Some(bases[from..=through].to_vec())
}

fn trim_cigar_by_reference(cigar: &[u32], start: usize, end: usize) -> Vec<u32> {
    let mut builder = CigarBuilder::new(true);
    let mut element_end = 0usize;
    for &element in cigar {
        let element_start = element_end;
        element_end = element_start + usize::from(consumes_reference(op(element))) * len(element);
        if element_end < start || (element_end == start && element_start < start) {
            continue;
        }
        if element_start > end && element_end > end + 1 {
            break;
        }
        let overlap = if element_end == element_start {
            len(element)
        } else {
            (end + 1).min(element_end) - start.max(element_start)
        };
        builder.add(ce(overlap, op(element)));
    }
    assert!(element_end > end, "haplotype CIGAR does not reach trim end");
    builder.finish()
}

fn trim_haplotype(haplotype: &TrimHap, start: i64, end: i64) -> Option<TrimHap> {
    let new_start = usize::try_from(start - haplotype.start).expect("trim starts before haplotype");
    let new_end = new_start + usize::try_from(end - start).expect("invalid trim interval");
    let bases =
        bases_covering_reference_interval(new_start, new_end, &haplotype.bases, &haplotype.cigar)?;
    if bases.is_empty() {
        return None;
    }

    let cigar = trim_cigar_by_reference(&haplotype.cigar, new_start, new_end);
    let leading_insertion = cigar.first().is_some_and(|&e| !consumes_reference(op(e)));
    let trailing_insertion = cigar.last().is_some_and(|&e| !consumes_reference(op(e)));
    let first = usize::from(leading_insertion);
    let last = cigar.len() - usize::from(trailing_insertion);
    if last <= first {
        return None;
    }
    let cigar = if leading_insertion || trailing_insertion {
        let mut builder = CigarBuilder::new(false);
        for &element in &cigar[first..last] {
            builder.add(element);
        }
        builder.finish()
    } else {
        cigar
    };

    Some(TrimHap {
        bases,
        is_ref: false,
        alignment_start_hap_wrt_ref: haplotype.alignment_start_hap_wrt_ref + new_start as i64,
        cigar,
        contig: haplotype.contig.clone(),
        start,
        end,
    })
}

fn trim_haplotypes(haplotypes: &[TrimHap], start: i64, end: i64) -> Vec<TrimHap> {
    let mut collapsed: Vec<(TrimHap, bool)> = Vec::new();
    for original in haplotypes {
        let Some(trimmed) = trim_haplotype(original, start, end) else {
            assert!(!original.is_ref, "trimming eliminated reference haplotype");
            continue;
        };
        if let Some(existing) = collapsed
            .iter_mut()
            .find(|(candidate, _)| candidate.bases == trimmed.bases)
        {
            if original.is_ref {
                *existing = (trimmed, true);
            }
        } else {
            collapsed.push((trimmed, original.is_ref));
        }
    }
    let mut result: Vec<_> = collapsed
        .into_iter()
        .map(|(mut haplotype, is_ref)| {
            haplotype.is_ref = is_ref;
            haplotype
        })
        .collect();
    result.sort_by(|left, right| {
        left.bases
            .len()
            .cmp(&right.bases.len())
            .then_with(|| left.bases.cmp(&right.bases))
    });
    assert!(result.iter().any(|haplotype| haplotype.is_ref));
    result
}

fn unclipped_read_length(read: &HcRead) -> usize {
    let soft_clipped: usize = read
        .cigar
        .iter()
        .filter(|&&element| op(element) == S)
        .map(|&element| len(element))
        .sum();
    read.bases.len() - soft_clipped
}

fn mate_on_same_contig_or_no_mapped_mate(read: &HcRead) -> bool {
    read.flag & 0x1 == 0
        || read.flag & 0x8 != 0
        || (read.flag & 0x4 == 0 && read.contig == read.mate_contig)
}

pub(crate) fn trim_reads(reads: &[HcRead], start: i64, end: i64) -> (Vec<HcRead>, Vec<HcRead>) {
    let mut clipped: Vec<_> = reads
        .iter()
        .cloned()
        .map(|read| hard_clip_to_region(read, start, end))
        .filter(|read| {
            read.flag & 0x4 == 0
                && !read.bases.is_empty()
                && read.start <= end
                && read.end() >= start
        })
        .collect();
    clipped.sort_by(coordinate_cmp);

    // `callRegion` removes these stubs before `filterNonPassingReads`, so they
    // deliberately do not appear in the returned filtered collection.
    clipped.retain(|read| unclipped_read_length(read) >= 10);
    let mut passing = Vec::with_capacity(clipped.len());
    let mut filtered = Vec::new();
    for read in clipped {
        if unclipped_read_length(&read) < 10
            || read.mapq < 20
            || !mate_on_same_contig_or_no_mapped_mate(&read)
        {
            filtered.push(read);
        } else {
            passing.push(read);
        }
    }
    (passing, filtered)
}

/// Reproduce HaplotypeCaller's default post-assembly trim for one active region.
#[allow(clippy::too_many_arguments)]
pub fn trim_region(
    contig: &str,
    start: i64,
    end: i64,
    contig_len: i64,
    padded_ref_start: i64,
    padded_ref: &[u8],
    haplotypes: &[TrimHap],
    reads: &[HcRead],
) -> TrimResult {
    assert!(start <= end && end <= contig_len);
    assert!(
        haplotypes
            .iter()
            .all(|haplotype| haplotype.contig == contig)
    );
    let mut events: Vec<_> = haplotypes
        .iter()
        .flat_map(|haplotype| events_for_haplotype(haplotype, padded_ref, padded_ref_start))
        .filter(|event| event.start <= end && event.end() >= start)
        .collect();
    if events.is_empty() {
        return TrimResult::NoVariation;
    }
    // This is the ordering used by AssemblyResultSet's TreeSet. It is not
    // needed for min/max, but de-duplicating by the same comparator preserves
    // the exact event collection semantics before STR padding.
    events.sort_by(|left, right| {
        left.start
            .cmp(&right.start)
            .then_with(|| left.reference.len().cmp(&right.reference.len()))
            .then_with(|| left.alternate.cmp(&right.alternate))
    });
    events.dedup_by(|left, right| {
        left.start == right.start
            && left.reference.len() == right.reference.len()
            && left.alternate == right.alternate
    });

    let variant_start = events
        .iter()
        .map(|event| event.start)
        .min()
        .unwrap()
        .max(start);
    let variant_end = events.iter().map(Event::end).max().unwrap().min(end);
    let old_padded_start = (start - 100).max(1);
    let old_padded_end = (end + 100).min(contig_len);
    let mut padded_start = events.iter().map(|event| event.start).min().unwrap();
    let mut padded_end = events.iter().map(Event::end).max().unwrap();

    for event in &events {
        let mut padding = if event.is_indel() { 75 } else { 20 };
        if event.is_indel() {
            let context_start = usize::try_from(event.start + 1 - padded_ref_start)
                .expect("event before padded reference context");
            if let Some((ref_repeats, alt_repeats, unit_len)) =
                tandem_repeat(event, &padded_ref[context_start..])
            {
                padding = 75 + ref_repeats.max(alt_repeats) as i64 * unit_len as i64;
            }
        }
        padded_start = padded_start.min((event.start - padding).max(1));
        padded_end = padded_end.max(event.end() + padding);
    }
    padded_start = padded_start.max(old_padded_start);
    padded_end = padded_end.min(old_padded_end);

    let trimmed_haplotypes = trim_haplotypes(haplotypes, padded_start, padded_end);
    let variation_present =
        haplotypes.iter().any(|haplotype| !haplotype.is_ref) && trimmed_haplotypes.len() > 1;
    let (trimmed_reads, filtered) = trim_reads(reads, padded_start, padded_end);

    TrimResult::Trimmed {
        span: (variant_start, variant_end),
        padded_span: (padded_start, padded_end),
        variation_present,
        haplotypes: trimmed_haplotypes,
        reads: trimmed_reads,
        filtered,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn event_minimal_representation_removes_common_suffix() {
        let event = Event::new(10, b"ACG".to_vec(), b"ATG".to_vec());
        assert_eq!(event.reference, b"AC");
        assert_eq!(event.alternate, b"AT");
        assert_eq!(event.end(), 11);
    }

    #[test]
    fn tandem_repeat_counts_context_after_alleles() {
        let deletion = Event::new(10, b"ACACA".to_vec(), b"A".to_vec());
        assert_eq!(tandem_repeat(&deletion, b"CACACACAG"), Some((4, 2, 2)));
    }

    #[test]
    fn equal_trimmed_bases_keep_reference_original() {
        let common = TrimHap {
            bases: b"AACAA".to_vec(),
            is_ref: false,
            alignment_start_hap_wrt_ref: 0,
            cigar: vec![ce(5, M)],
            contig: "chr1".into(),
            start: 1,
            end: 5,
        };
        let mut reference = common.clone();
        reference.is_ref = true;
        let result = trim_haplotypes(&[common, reference], 2, 4);
        assert_eq!(result.len(), 1);
        assert!(result[0].is_ref);
        assert_eq!(result[0].bases, b"ACA");
    }

    #[test]
    fn cigar_trim_keeps_internal_boundary_insertion() {
        let cigar = vec![ce(2, M), ce(1, I), ce(3, M)];
        assert_eq!(
            trim_cigar_by_reference(&cigar, 2, 4),
            vec![ce(1, I), ce(3, M)]
        );
    }

    #[test]
    fn operator_helpers_agree_for_haplotype_ops() {
        assert!(consumes_read(M));
        assert!(consumes_read(I));
        assert!(!consumes_read(D));
    }
}
