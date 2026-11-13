//! AssemblyBasedCallerUtils physical phasing. Sets are used only for membership,
//! size and containment: no hash-table iteration establishes a phase or tie.
use nomics_hc::trim::Event;
use std::{
    collections::BTreeSet,
    io::{self, Write},
};

pub(super) struct Call {
    pub start: i64,
    pub alleles: Vec<Vec<u8>>,
    pub haplotypes: Vec<usize>,
}

/// (phase group, alt in first slot). A contradictory new connection clears all
/// mappings, exactly as the source; existing partner mappings are not rewritten.
pub fn phase_sets(sets: &[BTreeSet<usize>]) -> Vec<Option<(usize, bool)>> {
    let total = sets
        .iter()
        .flatten()
        .copied()
        .collect::<BTreeSet<_>>()
        .len();
    let mut result = vec![None; sets.len()];
    let mut next = 0;
    for i in 0..sets.len().saturating_sub(1) {
        let call = &sets[i];
        if call.is_empty() {
            continue;
        }
        let all = call.len() == total;
        let mut available = call.clone();
        for j in i + 1..sets.len() {
            let comp = &sets[j];
            if comp.is_empty() {
                continue;
            }
            let same = call == comp || (all && comp.is_subset(&available)) || comp.len() == total;
            let opposite = !same && call.len() + comp.len() == total && call.is_disjoint(comp);
            if !same && !opposite {
                continue;
            }
            if let Some((group, phase)) = result[i] {
                if result[j].is_none() {
                    result[j] = Some((group, phase ^ opposite));
                }
            } else {
                if result[j].is_some() {
                    return vec![None; sets.len()];
                }
                result[i] = Some((next, false));
                result[j] = Some((next, opposite));
                if same {
                    available.retain(|h| comp.contains(h));
                }
                next += 1;
            }
        }
    }
    result
}

pub(super) fn write(
    calls: &[Call],
    events: &[Vec<Event>],
    chrom: &str,
    samples: &[String],
    normal: &str,
    w: &mut impl Write,
) -> io::Result<()> {
    let called: BTreeSet<_> = calls
        .iter()
        .flat_map(|c| c.haplotypes.iter().copied())
        .collect();
    let sets: Vec<BTreeSet<_>> = calls
        .iter()
        .map(|c| {
            if c.alleles.len() != 2 {
                return BTreeSet::new();
            }
            called
                .iter()
                .copied()
                .filter(|&h| {
                    events[h]
                        .iter()
                        .any(|e| e.start == c.start && e.alternate == c.alleles[1])
                })
                .collect()
        })
        .collect();
    let phases = phase_sets(&sets);
    for (i, c) in calls.iter().enumerate() {
        writeln!(
            w,
            "V {chrom}:{} {}",
            c.start,
            c.alleles
                .iter()
                .map(|a| String::from_utf8_lossy(a))
                .collect::<Vec<_>>()
                .join(",")
        )?;
        for sample in samples {
            writeln!(w, "S {sample}")?;
            let gt = if let Some((_, alt_first)) = phases[i] {
                if sample == normal {
                    "0|0".to_owned()
                } else if alt_first {
                    "1|0".to_owned()
                } else {
                    "0|1".to_owned()
                }
            } else if sample == normal {
                "0/0".to_owned()
            } else {
                (0..c.alleles.len())
                    .map(|a| a.to_string())
                    .collect::<Vec<_>>()
                    .join("/")
            };
            writeln!(w, "GT {gt}")?;
            if let Some((group, alt_first)) = phases[i] {
                let first = phases
                    .iter()
                    .position(|p| p.is_some_and(|(g, _)| g == group))
                    .unwrap();
                let anchor = &calls[first];
                writeln!(w, "G PGT {}", if alt_first { "1|0" } else { "0|1" })?;
                writeln!(
                    w,
                    "G PID {}_{}_{}",
                    anchor.start,
                    String::from_utf8_lossy(&anchor.alleles[0]),
                    String::from_utf8_lossy(&anchor.alleles[1])
                )?;
                writeln!(w, "G PS {}", anchor.start)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_jar_phase_groups() {
        for line in include_str!("phasing.tsv").lines() {
            let (masks, expected) = line.split_once(' ').unwrap();
            let sets: Vec<BTreeSet<_>> = masks
                .split(',')
                .map(|s| {
                    let mask: u32 = s.parse().unwrap();
                    (0..8).filter(|h| mask & (1 << h) != 0).collect()
                })
                .collect();
            let actual = phase_sets(&sets)
                .iter()
                .map(|p| {
                    p.map_or_else(
                        || ".".to_string(),
                        |(g, b)| format!("{g}:{}", usize::from(b)),
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(actual, expected, "{masks}");
        }
    }
}
