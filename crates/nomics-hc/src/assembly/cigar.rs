pub(crate) const M: u8 = 0;
pub(crate) const I: u8 = 1;
pub(crate) const D: u8 = 2;
pub(crate) const S: u8 = 4;
pub(crate) const H: u8 = 5;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Elem {
    pub(crate) len: usize,
    pub(crate) op: u8,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Strategy {
    Softclip,
    LeadingIndel,
}

pub(crate) fn unpack(cigar: &[u32]) -> Vec<Elem> {
    cigar
        .iter()
        .map(|value| Elem {
            len: (value >> 4) as usize,
            op: (value & 15) as u8,
        })
        .collect()
}

pub(crate) fn pack(cigar: &[Elem]) -> Vec<u32> {
    cigar
        .iter()
        .filter(|elem| elem.len != 0)
        .map(|elem| {
            (u32::try_from(elem.len).expect("CIGAR element too long") << 4) | u32::from(elem.op)
        })
        .collect()
}

pub(crate) fn consumes_read(op: u8) -> bool {
    matches!(op, M | I | S | 7 | 8)
}

pub(crate) fn consumes_ref(op: u8) -> bool {
    matches!(op, M | D | 3 | 7 | 8)
}

pub(super) fn cigar_reference_len(cigar: &[u32]) -> usize {
    cigar
        .iter()
        .filter(|value| consumes_ref((**value & 15) as u8))
        .map(|value| (*value >> 4) as usize)
        .sum()
}

pub(crate) fn hard_clip_soft_clips(
    bases: &[u8],
    quals: &[u8],
    cigar: &[u32],
) -> (Vec<u8>, Vec<u8>) {
    let left = cigar
        .iter()
        .take_while(|value| matches!((**value & 15) as u8, H | S))
        .filter(|value| (**value & 15) as u8 == S)
        .map(|value| (*value >> 4) as usize)
        .sum::<usize>();
    let right = cigar
        .iter()
        .rev()
        .take_while(|value| matches!((**value & 15) as u8, H | S))
        .filter(|value| (**value & 15) as u8 == S)
        .map(|value| (*value >> 4) as usize)
        .sum::<usize>();
    let end = bases.len().saturating_sub(right).max(left);
    (bases[left..end].to_vec(), quals[left..end].to_vec())
}

fn push(cigar: &mut Vec<Elem>, elem: Elem) {
    if elem.len == 0 {
        return;
    }
    if let Some(last) = cigar.last_mut()
        && last.op == elem.op
    {
        last.len += elem.len;
    } else {
        cigar.push(elem);
    }
}

pub(super) fn java_last_index_of(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(haystack.len());
    }
    haystack
        .windows(needle.len())
        .rposition(|window| window == needle)
}

/// GATK `SmithWatermanJavaAligner.align`. Dispatches to the wavefront
/// kernel in [`super::sw`], which is exact by construction against
/// [`smith_waterman_reference`]; `NOMICS_DUMP_SW=<path>` records every call.
pub(crate) fn smith_waterman(
    reference: &[u8],
    alternate: &[u8],
    match_value: i32,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
    strategy: Strategy,
) -> (Vec<Elem>, i32) {
    let weights = [match_value, mismatch, gap_open, gap_extend];
    let result = super::sw::smith_waterman_fast(reference, alternate, weights, strategy);
    if let Some(sink) = super::sw::dump_sink() {
        // Record the oracle's answer, not the fast path's.
        let oracle = smith_waterman_reference(
            reference,
            alternate,
            match_value,
            mismatch,
            gap_open,
            gap_extend,
            strategy,
        );
        super::sw::dump_call(sink, reference, alternate, weights, strategy, &oracle);
    }
    result
}

/// The straight port of `SmithWatermanJavaAligner` (full `i32` back-pointer
/// matrix, row by row). Kept as the oracle for the fast paths.
pub(crate) fn smith_waterman_reference(
    reference: &[u8],
    alternate: &[u8],
    match_value: i32,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
    strategy: Strategy,
) -> (Vec<Elem>, i32) {
    if matches!(strategy, Strategy::Softclip)
        && let Some(offset) = java_last_index_of(reference, alternate)
    {
        return (
            vec![Elem {
                len: alternate.len(),
                op: M,
            }],
            i32::try_from(offset).expect("alignment offset too large"),
        );
    }

    let rows = reference.len() + 1;
    let cols = alternate.len() + 1;
    let mut back = vec![0_i32; rows * cols];
    let mut previous = vec![0_i32; cols];
    let mut current = vec![0_i32; cols];
    if matches!(strategy, Strategy::LeadingIndel) {
        if cols > 1 {
            previous[1] = gap_open;
            for j in 2..cols {
                previous[j] = previous[j - 1] + gap_extend;
            }
        }
    }
    let low = i32::MIN / 2;
    let mut best_v = vec![low; cols + 1];
    let mut size_v = vec![0_i32; cols + 1];
    let mut p1 = 0_usize;
    let mut max_score = i32::MIN;
    for i in 1..rows {
        current[0] = if matches!(strategy, Strategy::LeadingIndel) {
            if i == 1 {
                gap_open
            } else {
                previous[0] + gap_extend
            }
        } else {
            0
        };
        let mut best_h = low;
        let mut size_h = 0_i32;
        for j in 1..cols {
            let diagonal = previous[j - 1]
                + if reference[i - 1] == alternate[j - 1] {
                    match_value
                } else {
                    mismatch
                };
            let opened_v = previous[j] + gap_open;
            best_v[j] += gap_extend;
            if opened_v > best_v[j] {
                best_v[j] = opened_v;
                size_v[j] = 1;
            } else {
                size_v[j] += 1;
            }
            let opened_h = current[j - 1] + gap_open;
            best_h += gap_extend;
            if opened_h > best_h {
                best_h = opened_h;
                size_h = 1;
            } else {
                size_h += 1;
            }
            let index = i * cols + j;
            if diagonal >= best_v[j] && diagonal >= best_h {
                current[j] = diagonal.max(-100_000_000);
                back[index] = 0;
            } else if best_h >= best_v[j] {
                current[j] = best_h.max(-100_000_000);
                back[index] = -size_h;
            } else {
                current[j] = best_v[j].max(-100_000_000);
                back[index] = size_v[j];
            }
        }
        let score = current[alternate.len()];
        if score >= max_score {
            p1 = i;
            max_score = score;
        }
        std::mem::swap(&mut previous, &mut current);
    }

    let mut p2 = alternate.len();
    let mut trailing = 0_usize;
    if matches!(strategy, Strategy::Softclip) {
        for j in 1..cols {
            let score = previous[j];
            if score > max_score
                || (score == max_score && reference.len().abs_diff(j) < p1.abs_diff(p2))
            {
                p1 = reference.len();
                p2 = j;
                max_score = score;
                trailing = alternate.len() - j;
            }
        }
    }

    // The Java routine builds reverse-order elements and reverses the list at
    // the end.  Preserve that behavior rather than coalescing across states.
    let mut reverse = Vec::<Elem>::new();
    if trailing > 0 && matches!(strategy, Strategy::Softclip) {
        reverse.push(Elem {
            len: trailing,
            op: S,
        });
    }
    let mut state = M;
    let mut segment = 0_usize;
    while p1 > 0 && p2 > 0 {
        let btr = back[p1 * cols + p2];
        let (new_state, step) = if btr > 0 {
            (D, btr as usize)
        } else if btr < 0 {
            (I, (-btr) as usize)
        } else {
            (M, 1)
        };
        match new_state {
            M => {
                p1 -= 1;
                p2 -= 1;
            }
            I => p2 -= step,
            D => p1 -= step,
            _ => unreachable!(),
        }
        if new_state == state {
            segment += step;
        } else {
            if segment != 0 {
                reverse.push(Elem {
                    len: segment,
                    op: state,
                });
            }
            segment = step;
            state = new_state;
        }
    }
    let offset = if matches!(strategy, Strategy::Softclip) {
        reverse.push(Elem {
            len: segment,
            op: state,
        });
        if p2 > 0 {
            reverse.push(Elem { len: p2, op: S });
        }
        i32::try_from(p1).expect("offset too large")
    } else {
        reverse.push(Elem {
            len: segment,
            op: state,
        });
        if p1 > 0 {
            reverse.push(Elem { len: p1, op: D });
        } else if p2 > 0 {
            reverse.push(Elem { len: p2, op: I });
        }
        0
    };
    reverse.reverse();
    (reverse, offset)
}

#[derive(Default)]
pub(crate) struct BuildResult {
    pub(crate) cigar: Vec<Elem>,
    pub(crate) leading_deletions: usize,
    pub(crate) trailing_deletions: usize,
}

pub(crate) fn cigar_builder(elements: impl IntoIterator<Item = Elem>) -> BuildResult {
    let mut result = BuildResult::default();
    for elem in elements {
        if elem.len == 0 {
            continue;
        }
        if elem.op == D
            && (result.cigar.is_empty()
                || result.cigar.last().is_some_and(|e| e.op == S || e.op == H))
        {
            result.leading_deletions += elem.len;
            continue;
        }
        if elem.op == D
            && result.cigar.last().is_some_and(|e| e.op == I)
            && (result.cigar.len() == 1 || matches!(result.cigar[result.cigar.len() - 2].op, S | H))
        {
            result.leading_deletions += elem.len;
            continue;
        }
        if elem.op == D && result.cigar.last().is_some_and(|e| e.op == I) {
            let at = result.cigar.len() - 1;
            if at > 0 && result.cigar[at - 1].op == D {
                result.cigar[at - 1].len += elem.len;
            } else {
                result.cigar.insert(at, elem);
            }
        } else {
            push(&mut result.cigar, elem);
        }
    }
    if result.cigar.last().is_some_and(|elem| elem.op == D) {
        result.trailing_deletions += result.cigar.pop().expect("last existed").len;
    } else if result.cigar.last().is_some_and(|elem| elem.op == I)
        && result.cigar.len() > 1
        && result.cigar[result.cigar.len() - 2].op == D
    {
        let at = result.cigar.len() - 2;
        result.trailing_deletions += result.cigar.remove(at).len;
    }
    result
}

pub(crate) fn trim_by_bases(cigar: &[Elem], start: usize, end: usize) -> BuildResult {
    let mut kept = Vec::new();
    let mut element_end = 0_usize;
    for elem in cigar {
        let element_start = element_end;
        element_end += if consumes_read(elem.op) { elem.len } else { 0 };
        if element_end < start || (element_end == start && element_start < start) {
            continue;
        }
        if element_start > end && element_end > end + 1 {
            break;
        }
        let overlap = if element_end == element_start {
            elem.len
        } else {
            (end + 1).min(element_end) - start.max(element_start)
        };
        kept.push(Elem {
            len: overlap,
            op: elem.op,
        });
    }
    cigar_builder(kept)
}

#[derive(Clone, Copy)]
struct Range {
    start: usize,
    end: usize,
}

impl Range {
    fn len(self) -> usize {
        self.end - self.start
    }
    fn shift_left(&mut self, amount: usize) {
        self.start -= amount;
        self.end -= amount;
    }
}

fn normalize_alleles(
    reference: &[u8],
    alternate: &[u8],
    ref_range: &mut Range,
    alt_range: &mut Range,
    max_shift: usize,
) -> (isize, usize) {
    let mut start_shift = 0_isize;
    let mut end_shift = 0_usize;
    let mut min_size = ref_range.len().min(alt_range.len());
    while min_size > 0 && reference[ref_range.end - 1] == alternate[alt_range.end - 1] {
        ref_range.end -= 1;
        alt_range.end -= 1;
        min_size -= 1;
        end_shift += 1;
    }
    while min_size > 0 && reference[ref_range.start] == alternate[alt_range.start] {
        ref_range.start += 1;
        alt_range.start += 1;
        min_size -= 1;
        start_shift -= 1;
    }
    while start_shift < max_shift as isize
        && reference[ref_range.start - 1] == alternate[alt_range.start - 1]
        && reference[ref_range.end - 1] == alternate[alt_range.end - 1]
    {
        ref_range.shift_left(1);
        alt_range.shift_left(1);
        start_shift += 1;
        end_shift += 1;
    }
    (start_shift, end_shift)
}

pub(crate) fn left_align(
    cigar: &[Elem],
    reference: &[u8],
    alternate: &[u8],
    read_start: usize,
) -> BuildResult {
    if cigar.iter().all(|elem| !matches!(elem.op, I | D)) {
        return BuildResult {
            cigar: cigar.to_vec(),
            ..BuildResult::default()
        };
    }
    let ref_len: usize = cigar
        .iter()
        .filter(|e| consumes_ref(e.op))
        .map(|e| e.len)
        .sum();
    let mut ref_range = Range {
        start: read_start + ref_len,
        end: read_start + ref_len,
    };
    let mut alt_range = Range {
        start: alternate.len(),
        end: alternate.len(),
    };
    let mut reverse = Vec::new();
    for (index, elem) in cigar.iter().enumerate().rev() {
        if matches!(elem.op, I | D) {
            if consumes_ref(elem.op) {
                ref_range.start -= elem.len;
            }
            if consumes_read(elem.op) {
                alt_range.start -= elem.len;
            }
        } else if ref_range.len() == 0 && alt_range.len() == 0 {
            reverse.push(*elem);
            if consumes_ref(elem.op) {
                ref_range.shift_left(elem.len);
            }
            if consumes_read(elem.op) {
                alt_range.shift_left(elem.len);
            }
        } else {
            let max_shift = if matches!(elem.op, M | 7 | 8) {
                elem.len
            } else {
                0
            };
            let (start_shift, end_shift) = normalize_alleles(
                reference,
                alternate,
                &mut ref_range,
                &mut alt_range,
                max_shift,
            );
            reverse.push(Elem {
                len: end_shift,
                op: M,
            });
            let emit =
                index == 0 || start_shift < max_shift as isize || !matches!(elem.op, M | 7 | 8);
            let new_left_match = if start_shift < 0 {
                (-start_shift) as usize
            } else {
                0
            };
            let remaining = if start_shift < 0 {
                elem.len
            } else {
                elem.len - start_shift as usize
            };
            if emit {
                reverse.push(Elem {
                    len: ref_range.len(),
                    op: D,
                });
                reverse.push(Elem {
                    len: alt_range.len(),
                    op: I,
                });
                ref_range.end = ref_range.start;
                alt_range.end = alt_range.start;
                let ref_shift = new_left_match + if consumes_ref(elem.op) { remaining } else { 0 };
                let alt_shift = new_left_match + if consumes_read(elem.op) { remaining } else { 0 };
                ref_range.shift_left(ref_shift);
                alt_range.shift_left(alt_shift);
            }
            reverse.push(Elem {
                len: new_left_match,
                op: M,
            });
            reverse.push(Elem {
                len: remaining,
                op: elem.op,
            });
        }
    }
    reverse.push(Elem {
        len: ref_range.len(),
        op: D,
    });
    reverse.push(Elem {
        len: alt_range.len(),
        op: I,
    });
    reverse.reverse();
    cigar_builder(reverse)
}

pub(super) fn calculate_cigar(reference: &[u8], alternate: &[u8]) -> Option<Vec<u32>> {
    if alternate.is_empty() {
        return Some(vec![
            (u32::try_from(reference.len()).ok()? << 4) | u32::from(D),
        ]);
    }
    if reference == alternate {
        return Some(vec![u32::try_from(reference.len()).ok()? << 4]);
    }
    let mut padded_ref = Vec::with_capacity(reference.len() + 20);
    let mut padded_alt = Vec::with_capacity(alternate.len() + 20);
    padded_ref.extend_from_slice(b"NNNNNNNNNN");
    padded_ref.extend_from_slice(reference);
    padded_ref.extend_from_slice(b"NNNNNNNNNN");
    padded_alt.extend_from_slice(b"NNNNNNNNNN");
    padded_alt.extend_from_slice(alternate);
    padded_alt.extend_from_slice(b"NNNNNNNNNN");
    let (raw, offset) = smith_waterman(
        &padded_ref,
        &padded_alt,
        200,
        -150,
        -260,
        -11,
        Strategy::Softclip,
    );
    if offset > 0 || raw.iter().any(|elem| elem.op == S) {
        return None;
    }
    let mut trimmed = trim_by_bases(&raw, 10, padded_alt.len() - 11);
    let leading = trimmed.leading_deletions;
    if trimmed.trailing_deletions != 0 {
        trimmed.cigar.push(Elem {
            len: trimmed.trailing_deletions,
            op: D,
        });
    }
    let aligned = left_align(&trimmed.cigar, reference, alternate, leading);
    let mut final_cigar = Vec::new();
    push(
        &mut final_cigar,
        Elem {
            len: leading + aligned.leading_deletions,
            op: D,
        },
    );
    for elem in aligned.cigar {
        push(&mut final_cigar, elem);
    }
    push(
        &mut final_cigar,
        Elem {
            len: aligned.trailing_deletions,
            op: D,
        },
    );
    Some(pack(&final_cigar))
}

pub(super) fn dangling_cigar(reference: &[u8], alternate: &[u8]) -> Vec<u32> {
    let (mut cigar, _) = smith_waterman(
        reference,
        alternate,
        25,
        -50,
        -110,
        -6,
        Strategy::LeadingIndel,
    );
    if cigar.last().is_some_and(|elem| elem.op == D) {
        cigar.pop();
    }
    pack(&cigar)
}
