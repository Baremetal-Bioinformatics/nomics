//! AlleleLikelihoods.groupEvidence and Fragment.createAndAvoidFailure.
use nomics_hc::{finalize::HcRead, jmath::log};

#[derive(Debug)]
pub struct Fragment {
    /// All evidence contributes to likelihoods, even when not represented below.
    pub evidence: Vec<usize>,
    /// The one/two reads exposed to downstream fragment annotations.
    pub reads: Vec<usize>,
    pub start: i64,
    pub end: i64,
}

/// Apply the second, natural-log-scale Q45 cap in place. This literal operation
/// order is qual * (-Math.log(10)/10), not -4.5 * Math.log(10).
pub fn cap(matrix: &mut [Vec<f64>]) {
    let n = matrix.first().map_or(0, Vec::len);
    assert!(matrix.iter().all(|row| row.len() == n));
    let difference = 45.0 * (-log(10.0) / 10.0);
    for column in 0..n {
        let best = matrix
            .iter()
            .map(|row| row[column])
            .fold(f64::NEG_INFINITY, f64::max);
        let floor = best + difference;
        for row in matrix.iter_mut() {
            if row[column] < floor {
                row[column] = floor;
            }
        }
    }
}

/// Groups are in OpenJDK groupingBy HashMap.values order; members retain the
/// original evidence order. Neither sorted names nor insertion order is valid.
pub fn group(reads: &[HcRead], matrix: &[Vec<f64>]) -> (Vec<Fragment>, Vec<Vec<f64>>) {
    assert!(matrix.iter().all(|row| row.len() == reads.len()));
    let names: Vec<_> = reads
        .iter()
        .map(|r| String::from_utf8_lossy(&r.name))
        .collect();
    let groups = crate::java_hashmap::group_indices(names.iter().map(|s| s.as_ref()));
    let mut values = vec![vec![0.0; groups.len()]; matrix.len()];
    let fragments = groups
        .into_iter()
        .enumerate()
        .map(|(column, group)| {
            for (a, row) in values.iter_mut().enumerate() {
                for &r in &group {
                    row[column] += matrix[a][r];
                }
            }
            let chosen = if group.len() <= 2 {
                group.clone()
            } else {
                let usable: Vec<_> = group
                    .iter()
                    .copied()
                    .filter(|&i| reads[i].flag & (0x100 | 0x400 | 0x800) == 0)
                    .take(2)
                    .collect();
                if usable.is_empty() {
                    vec![group[0]]
                } else {
                    usable
                }
            };
            let (start, end) = if chosen.len() == 1 {
                let r = &reads[chosen[0]];
                (r.start.min(r.end()), r.start.max(r.end()))
            } else {
                let start = reads[chosen[0]].start.min(reads[chosen[1]].start);
                let end = reads[chosen[0]].end().max(reads[chosen[1]].end());
                (start.min(end), start.max(end))
            };
            Fragment {
                evidence: group,
                reads: chosen,
                start,
                end,
            }
        })
        .collect();
    (fragments, values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn independent_gatk_fragment_selection_and_hash_bins() {
        let fixture = include_str!("fragments.tsv");
        let mut reads = Vec::new();
        let mut m = vec![Vec::new(), Vec::new()];
        for line in fixture.lines().take_while(|l| l.starts_with("R ")) {
            let f: Vec<_> = line.split_whitespace().collect();
            let length: usize = f[4].parse().unwrap();
            reads.push(HcRead {
                name: f[1].as_bytes().to_vec(),
                flag: f[2].parse().unwrap(),
                contig: 0,
                start: f[3].parse().unwrap(),
                cigar: vec![(length as u32) << 4],
                bases: vec![b'A'; length],
                quals: vec![0; length],
                mapq: 0,
                mate_contig: -1,
                mate_start: 0,
                tlen: 0,
                record: reads.len() as u32,
            });
            for a in 0..2 {
                m[a].push(f64::from_bits(u64::from_str_radix(f[5 + a], 16).unwrap()));
            }
        }
        cap(&mut m);
        let (fragments, grouped) = group(&reads, &m);
        let mut actual = Vec::new();
        super::super::matrix(&mut actual, "C", &m).unwrap();
        for f in fragments {
            writeln!(
                actual,
                "F {} {} {} {}",
                f.start,
                f.end,
                f.evidence
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                f.reads
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
            .unwrap();
        }
        super::super::matrix(&mut actual, "G", &grouped).unwrap();
        let expected = &fixture[fixture.find("C 2 ").unwrap()..];
        assert_eq!(String::from_utf8(actual).unwrap(), expected);
    }
}
