//! GATK 4.6.2.0 `AssemblyBasedCallerUtils.finalizeRegion` read finalization.
//!
//! This module intentionally keeps GATK's clipping and overlapping-fragment
//! behaviour, including a few surprising edge cases. CIGAR operations use
//! BAM's numeric encoding (`length << 4 | operator`).

use std::cmp::Ordering;
use std::collections::HashMap;

/// The part of a BAM record used by HaplotypeCaller's assembly preparation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HcRead {
    pub name: Vec<u8>,
    pub flag: u16,
    pub contig: i32,
    pub start: i64,
    pub cigar: Vec<u32>,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    pub mapq: u8,
    pub mate_contig: i32,
    pub mate_start: i64,
    pub tlen: i32,
    pub record: u32,
}

impl HcRead {
    pub fn from_record(view: &nomics_io::bam::RecordView<'_>, record: u32) -> HcRead {
        HcRead {
            name: view.read_name().to_vec(),
            flag: view.flag(),
            contig: view.ref_id(),
            start: i64::from(view.pos()) + 1,
            cigar: view.cigar().collect(),
            bases: (0..view.len())
                .map(|i| view.base(i).to_ascii_uppercase())
                .collect(),
            quals: view.quals().unwrap_or(&[]).to_vec(),
            mapq: view.mapq(),
            mate_contig: view.mate_ref_id(),
            mate_start: i64::from(view.mate_pos()) + 1,
            tlen: view.tlen(),
            record,
        }
    }

    /// Alignment end as htsjdk 4.2.0's `SAMRecord.getAlignmentEnd`:
    /// `start + referenceLength - 1`, so `start - 1` for a read that aligns
    /// no reference bases (an all-insertion read after realignment), which
    /// then overlaps nothing; 0 for an unplaced read.
    pub fn end(&self) -> i64 {
        let reference_length: i64 = self
            .cigar
            .iter()
            .filter(|&&c| consumes_reference(op(c)))
            .map(|&c| i64::from(len(c)))
            .sum();
        if self.start == 0 {
            0
        } else {
            self.start + reference_length - 1
        }
    }

    fn cigar_read_length(&self) -> usize {
        self.cigar
            .iter()
            .filter(|&&c| consumes_read(op(c)))
            .map(|&c| len(c) as usize)
            .sum()
    }
}

const M: u32 = 0;
const I: u32 = 1;
const D: u32 = 2;
const N: u32 = 3;
const S: u32 = 4;
const H: u32 = 5;
const EQ: u32 = 7;
const X: u32 = 8;
const READ_INDEX_NOT_FOUND: i64 = -1;

#[inline]
fn op(c: u32) -> u32 {
    c & 0xf
}

#[inline]
fn len(c: u32) -> u32 {
    c >> 4
}

#[inline]
fn ce(length: u32, operator: u32) -> u32 {
    length << 4 | operator
}

#[inline]
fn consumes_read(operator: u32) -> bool {
    matches!(operator, M | I | S | EQ | X)
}

#[inline]
fn consumes_reference(operator: u32) -> bool {
    matches!(operator, M | D | N | EQ | X)
}

#[inline]
fn is_clipping(operator: u32) -> bool {
    matches!(operator, S | H)
}

/// The output-affecting normalization performed by GATK's CigarBuilder.
#[derive(Default)]
struct CigarBuilder {
    elements: Vec<u32>,
    last: Option<u32>,
}

impl CigarBuilder {
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

        if operator == D {
            if self.last.is_none() || self.last.is_some_and(is_clipping) {
                return;
            }
            if self.last == Some(I)
                && (self.elements.len() == 1
                    || is_clipping(op(self.elements[self.elements.len() - 2])))
            {
                return;
            }
        }

        if self.last == Some(operator) {
            let n = self.elements.len() - 1;
            self.elements[n] = ce(len(self.elements[n]) + length, operator);
        } else if self.last.is_none() {
            self.elements.push(element);
            self.last = Some(operator);
        } else if is_clipping(operator) {
            let last = self.last.expect("checked above");
            if !consumes_read(last) && !is_clipping(last) {
                self.elements.pop();
                self.elements.push(element);
                self.last = Some(operator);
            } else if self.last_two_deletion_insertion() {
                let insert = self.elements.pop().expect("insertion");
                self.elements.pop();
                self.elements.push(insert);
                self.elements.push(element);
                // Deliberately retain `last == I`, matching CigarBuilder.
            } else {
                self.elements.push(element);
                self.last = Some(operator);
            }
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
        if self.last == Some(D) {
            self.elements.pop();
        } else if self.last_two_deletion_insertion() {
            let n = self.elements.len() - 2;
            self.elements.remove(n);
        }
        self.elements
    }
}

fn clip_cigar(cigar: &[u32], start: usize, stop: usize, clipping_operator: u32) -> Vec<u32> {
    let clip_left = start == 0;
    let mut builder = CigarBuilder::default();
    let mut element_start = 0usize;

    for &element in cigar {
        let operator = op(element);
        let length = len(element) as usize;
        if operator == H {
            builder.add(element);
            continue;
        }
        let element_end = element_start + usize::from(consumes_read(operator)) * length;
        if element_end <= start || element_start >= stop {
            if consumes_read(operator) || (element_start != start && element_start != stop) {
                builder.add(element);
            }
        } else {
            let unclipped_length = if clip_left {
                element_end as isize - stop as isize
            } else {
                start as isize - element_start as isize
            };
            let clipped_length = length as isize - unclipped_length;
            if unclipped_length <= 0 {
                if consumes_read(operator) {
                    builder.add(ce(length as u32, clipping_operator));
                }
            } else if clip_left {
                builder.add(ce(clipped_length as u32, clipping_operator));
                builder.add(ce(unclipped_length as u32, operator));
            } else {
                builder.add(ce(unclipped_length as u32, operator));
                builder.add(ce(clipped_length as u32, clipping_operator));
            }
        }
        element_start = element_end;
    }
    builder.finish()
}

fn alignment_start_shift(cigar: &[u32], num_clipped: usize) -> i64 {
    let mut reference_bases_clipped = 0i64;
    let mut element_start = 0usize;
    for &element in cigar {
        let operator = op(element);
        if operator == H {
            continue;
        }
        let length = len(element) as usize;
        let element_end = element_start + usize::from(consumes_read(operator)) * length;
        if element_end <= num_clipped {
            if consumes_reference(operator) {
                reference_bases_clipped += length as i64;
            }
        } else if element_start < num_clipped {
            if consumes_reference(operator) {
                reference_bases_clipped += (num_clipped - element_start) as i64;
            }
            break;
        }
        element_start = element_end;
    }
    reference_bases_clipped
}

fn empty_read(mut read: HcRead) -> HcRead {
    read.flag |= 0x4;
    read.contig = -1;
    read.start = 0;
    read.cigar.clear();
    read.bases.clear();
    read.quals.clear();
    read
}

fn hard_clip_bases(mut read: HcRead, start: usize, stop: usize) -> HcRead {
    if start > stop || stop >= read.bases.len() {
        return read;
    }
    let new_length = read.bases.len() - (stop - start + 1);
    if new_length == 0 {
        return empty_read(read);
    }

    let old_cigar = std::mem::take(&mut read.cigar);
    read.cigar = clip_cigar(&old_cigar, start, stop + 1, H);
    if start == 0 {
        read.bases.drain(..=stop);
        read.quals.drain(..=stop);
        read.start += alignment_start_shift(&old_cigar, stop + 1);
    } else {
        read.bases.truncate(start);
        read.quals.truncate(start);
    }
    read
}

fn soft_start(read: &HcRead) -> i64 {
    let mut result = read.start;
    for &element in &read.cigar {
        match op(element) {
            S => result -= i64::from(len(element)),
            H => {}
            _ => break,
        }
    }
    result
}

fn revert_soft_clipped_bases(mut read: HcRead) -> HcRead {
    if read.bases.is_empty() {
        return read;
    }
    let has_terminal_clip = read
        .cigar
        .first()
        .into_iter()
        .chain(read.cigar.last())
        .any(|&c| is_clipping(op(c)));
    if !has_terminal_clip {
        return read;
    }

    let new_start = soft_start(&read);
    let mut builder = CigarBuilder::default();
    for &element in &read.cigar {
        builder.add(if op(element) == S {
            ce(len(element), M)
        } else {
            element
        });
    }
    read.cigar = builder.finish();

    if new_start <= 0 {
        read.start = 1;
        read = hard_clip_bases(read, 0, (-new_start) as usize);
        if read.flag & 0x4 == 0 {
            read.start = 1;
        }
    } else {
        read.start = new_start;
    }
    read
}

fn hard_clip_soft_clipped_bases(mut read: HcRead) -> HcRead {
    if read.bases.is_empty() {
        return read;
    }
    let mut read_index = 0usize;
    let mut cut_left = None;
    let mut cut_right = None;
    let mut right_tail = false;
    for &element in &read.cigar {
        let operator = op(element);
        let length = len(element) as usize;
        if operator == S {
            if right_tail {
                cut_right = Some(read_index);
            } else {
                cut_left = Some(read_index + length - 1);
            }
        } else if operator != H {
            right_tail = true;
        }
        if consumes_read(operator) {
            read_index += length;
        }
    }
    if let Some(right) = cut_right {
        let last = read.bases.len() - 1;
        read = hard_clip_bases(read, right, last);
    }
    if let Some(left) = cut_left {
        read = hard_clip_bases(read, 0, left);
    }
    read
}

fn hard_clip_low_quality_ends(mut read: HcRead, low_quality: u8) -> HcRead {
    if read.bases.is_empty() {
        return read;
    }
    let read_length = read.bases.len();
    let mut right = read_length;
    while right > 0 && read.quals[right - 1] <= low_quality {
        right -= 1;
    }
    let mut left = 0usize;
    while left < read_length && read.quals[left] <= low_quality {
        left += 1;
    }
    if left >= right {
        return empty_read(read);
    }
    if right < read_length {
        read = hard_clip_bases(read, right, read_length - 1);
    }
    if left > 0 {
        read = hard_clip_bases(read, 0, left - 1);
    }
    read
}

fn has_well_defined_fragment_size(read: &HcRead) -> bool {
    if read.tlen == 0
        || read.flag & 0x1 == 0
        || read.flag & 0x4 != 0
        || read.flag & 0x8 != 0
        || (read.flag & 0x10 != 0) == (read.flag & 0x20 != 0)
    {
        return false;
    }
    if read.flag & 0x10 != 0 {
        read.end() > read.mate_start
    } else {
        read.start <= read.mate_start + i64::from(read.tlen)
    }
}

fn adaptor_boundary(read: &HcRead) -> Option<i64> {
    if !has_well_defined_fragment_size(read) {
        None
    } else if read.flag & 0x10 != 0 {
        Some(read.mate_start - 1)
    } else {
        Some(read.start + i64::from(read.tlen).abs())
    }
}

fn read_index_for_reference_coordinate(read: &HcRead, coordinate: i64) -> (i64, Option<u32>) {
    let alignment_start = soft_start(read);
    if coordinate < alignment_start {
        return (READ_INDEX_NOT_FOUND, None);
    }
    let mut last_read = 0i64;
    let mut last_reference = alignment_start;
    for &element in &read.cigar {
        let operator = op(element);
        let first_read = last_read;
        let first_reference = last_reference;
        if consumes_read(operator) {
            last_read += i64::from(len(element));
        }
        if consumes_reference(operator) || operator == S {
            last_reference += i64::from(len(element));
        }
        if first_reference <= coordinate && coordinate < last_reference {
            let offset = first_read
                + if consumes_read(operator) {
                    coordinate - first_reference
                } else {
                    0
                };
            return (offset, Some(operator));
        }
    }
    (READ_INDEX_NOT_FOUND, None)
}

fn hard_clip_by_reference_left(read: HcRead, reference_stop: i64) -> HcRead {
    let (offset, operator) = read_index_for_reference_coordinate(&read, reference_stop);
    if offset == READ_INDEX_NOT_FOUND {
        return read;
    }
    let stop = offset - i64::from(!operator.is_some_and(consumes_read));
    if stop < 0 {
        read
    } else {
        hard_clip_bases(read, 0, stop as usize)
    }
}

fn hard_clip_by_reference_right(read: HcRead, reference_start: i64) -> HcRead {
    let (offset, _) = read_index_for_reference_coordinate(&read, reference_start);
    if offset == READ_INDEX_NOT_FOUND {
        return read;
    }
    let last = read.bases.len() - 1;
    hard_clip_bases(read, offset as usize, last)
}

fn hard_clip_both_reference(read: HcRead, left: i64, right: i64) -> HcRead {
    if read.bases.is_empty() || left == right {
        return empty_read(read);
    }
    let original = read.clone();
    let right_clipped = hard_clip_by_reference_right(read, right);
    if left > right_clipped.end() {
        empty_read(original)
    } else {
        hard_clip_by_reference_left(right_clipped, left)
    }
}

fn hard_clip_adaptor(read: HcRead) -> HcRead {
    let Some(boundary) = adaptor_boundary(&read) else {
        return read;
    };
    if boundary < read.start || boundary > read.end() {
        return read;
    }
    if read.flag & 0x10 != 0 {
        hard_clip_by_reference_left(read, boundary)
    } else {
        hard_clip_by_reference_right(read, boundary)
    }
}

pub(crate) fn hard_clip_to_region(read: HcRead, region_start: i64, region_end: i64) -> HcRead {
    let alignment_start = read.start;
    let alignment_end = read.end();
    if alignment_start <= region_end && alignment_end >= region_start {
        if alignment_start < region_start && alignment_end > region_end {
            hard_clip_both_reference(read, region_start - 1, region_end + 1)
        } else if alignment_start < region_start {
            hard_clip_by_reference_left(read, region_start - 1)
        } else if alignment_end > region_end {
            hard_clip_by_reference_right(read, region_end + 1)
        } else {
            read
        }
    } else {
        empty_read(read)
    }
}

pub(crate) fn coordinate_cmp(first: &HcRead, second: &HcRead) -> Ordering {
    first
        .contig
        .cmp(&second.contig)
        .then_with(|| first.start.cmp(&second.start))
        .then_with(|| (first.flag & 0x10 != 0).cmp(&(second.flag & 0x10 != 0)))
        .then_with(|| first.name.cmp(&second.name))
        .then_with(|| first.flag.cmp(&second.flag))
        .then_with(|| first.mapq.cmp(&second.mapq))
        .then_with(|| {
            if first.flag & 0x1 != 0 && second.flag & 0x1 != 0 {
                first
                    .mate_contig
                    .cmp(&second.mate_contig)
                    .then_with(|| first.mate_start.cmp(&second.mate_start))
            } else {
                Ordering::Equal
            }
        })
        .then_with(|| first.tlen.cmp(&second.tlen))
}

fn adjust_overlapping_pair(reads: &mut [HcRead], left: usize, right: usize) {
    let left_first = soft_start(&reads[left]) < soft_start(&reads[right]);
    let (first_index, second_index) = if left_first {
        (left, right)
    } else {
        (right, left)
    };
    let (first, second) = if first_index < second_index {
        let (before, after) = reads.split_at_mut(second_index);
        (&mut before[first_index], &mut after[0])
    } else {
        let (before, after) = reads.split_at_mut(first_index);
        (&mut after[0], &mut before[second_index])
    };

    if first.end() < second.start || first.contig != second.contig {
        return;
    }
    let (offset, operator) = read_index_for_reference_coordinate(first, second.start);
    if offset == READ_INDEX_NOT_FOUND || operator.is_some_and(is_clipping) {
        return;
    }
    let (first_end, _) = read_index_for_reference_coordinate(first, first.end());
    let (second_end, _) = read_index_for_reference_coordinate(second, second.end());
    let (second_offset, _) = read_index_for_reference_coordinate(second, second.start);
    let overlap = (first_end - offset).min(second_end - second_offset) + 1;
    if overlap <= 0 {
        return;
    }
    for i in 0..overlap as usize {
        let first_i = offset as usize + i;
        let second_i = second_offset as usize + i;
        if first.bases[first_i] == second.bases[second_i] {
            first.quals[first_i] = first.quals[first_i].min(20);
            second.quals[second_i] = second.quals[second_i].min(20);
        } else {
            first.quals[first_i] = 0;
            second.quals[second_i] = 0;
        }
    }
}

fn clean_overlapping_read_pairs(reads: &mut [HcRead]) {
    let mut waiting: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut pairs = Vec::new();
    for (index, read) in reads.iter().enumerate() {
        if read.flag & 0x1 == 0
            || read.flag & 0x8 != 0
            || read.mate_start == 0
            || read.mate_start > read.end()
        {
            continue;
        }
        if let Some(first) = waiting.remove(&read.name) {
            pairs.push((first, index));
        } else {
            waiting.insert(read.name.clone(), index);
        }
    }
    for (first, second) in pairs {
        adjust_overlapping_pair(reads, first, second);
    }
}

/// Reproduce the `readsToUse` output of GATK 4.6.2.0's `finalizeRegion` for
/// the HaplotypeCaller options used by this project.
pub fn finalize_region(contig_len: i64, start: i64, end: i64, reads: &[HcRead]) -> Vec<HcRead> {
    let padded_start = (start - 100).max(1);
    let padded_end = (end + 100).min(contig_len);
    let mut result = Vec::with_capacity(reads.len());

    for original in reads {
        let read = if has_well_defined_fragment_size(original) {
            revert_soft_clipped_bases(original.clone())
        } else {
            hard_clip_soft_clipped_bases(original.clone())
        };
        let read = hard_clip_low_quality_ends(read, 9);
        if read.flag & 0x4 != 0 || read.start > read.end() {
            continue;
        }
        let adaptor_clipped = hard_clip_adaptor(read);
        if adaptor_clipped.bases.is_empty() || adaptor_clipped.cigar_read_length() == 0 {
            continue;
        }
        let overlaps = adaptor_clipped.start <= padded_end && adaptor_clipped.end() >= padded_start;
        let clipped = hard_clip_to_region(adaptor_clipped, padded_start, padded_end);
        if clipped.flag & 0x4 == 0
            && clipped.start <= clipped.end()
            && !clipped.bases.is_empty()
            && overlaps
        {
            result.push(clipped);
        }
    }

    result.sort_by(coordinate_cmp);
    clean_overlapping_read_pairs(&mut result);
    result
}
