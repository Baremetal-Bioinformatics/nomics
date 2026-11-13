//! GATK 4.6.2.0 AdaptiveChainPruner; graph-owned chain order is retained.
use crate::mutect2_regions::numerics::constant_error_log_likelihood_ratio;
use nomics_hc::assembly::PruningChain;
use std::{cmp::Ordering, collections::BinaryHeap};

/// Indices of chains to remove, in their original order.
pub fn chains_to_remove(chains: &[PruningChain]) -> Vec<usize> {
    if chains.is_empty() {
        return Vec::new();
    }
    let errors = likely_errors(chains, 0.001);
    let error_count: i32 = errors
        .iter()
        .map(|&i| *chains[i].multiplicities.last().unwrap())
        .sum();
    let total: i32 = chains.iter().flat_map(|c| &c.multiplicities).sum();
    likely_errors(chains, f64::from(error_count) / f64::from(total))
        .into_iter()
        .filter(|&i| !chains[i].any_ref)
        .collect()
}

struct Candidate<'a> {
    index: usize,
    score: f64,
    chain: &'a PruningChain,
}
impl PartialEq for Candidate<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Candidate<'_> {}
impl PartialOrd for Candidate<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Candidate<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .total_cmp(&other.score)
            .then_with(|| other.chain.first_bases.cmp(&self.chain.first_bases))
            .then_with(|| other.chain.bases.cmp(&self.chain.bases))
    }
}
fn likely_errors(chains: &[PruningChain], error: f64) -> Vec<usize> {
    let odds: Vec<_> = chains
        .iter()
        .map(|c| {
            let left = c.multiplicities[0];
            let right = *c.multiplicities.last().unwrap();
            (
                if c.is_source {
                    0.0
                } else {
                    constant_error_log_likelihood_ratio(c.left_total - left, left, error)
                },
                if c.is_sink {
                    0.0
                } else {
                    constant_error_log_likelihood_ratio(c.right_total - right, right, error)
                },
            )
        })
        .collect();
    let vertices = chains.iter().map(|c| c.source.max(c.target)).max().unwrap() + 1;
    let mut incoming = vec![Vec::new(); vertices];
    let mut outgoing = vec![Vec::new(); vertices];
    let mut seedable = vec![0; vertices];
    let threshold = nomics_hc::jmath::log(10.0);
    for (i, c) in chains.iter().enumerate() {
        let (l, r) = odds[i];
        if l >= threshold || c.first_is_ref {
            outgoing[c.source].push(i);
        }
        if r >= threshold || c.first_is_ref {
            incoming[c.target].push(i);
        }
        if l >= threshold * 4.0 && r >= threshold * 4.0 {
            seedable[c.source] += 1;
            seedable[c.target] += 1;
        }
    }
    // Stream.max retains the first chain when its comparator returns equality.
    let mut best = 0;
    for i in 1..chains.len() {
        let key = |i: usize| {
            (
                *chains[i].multiplicities.iter().max().unwrap(),
                chains[i].multiplicities.len(),
                &chains[i].first_bases,
            )
        };
        if key(i) > key(best) {
            best = i;
        }
    }
    let mut queue = BinaryHeap::new();
    queue.push(Candidate {
        index: best,
        score: f64::INFINITY,
        chain: &chains[best],
    });
    let mut processed = vec![false; vertices];
    let enqueue = |v: usize, queue: &mut BinaryHeap<_>| {
        for &i in &outgoing[v] {
            queue.push(Candidate {
                index: i,
                score: odds[i].0,
                chain: &chains[i],
            });
        }
        for &i in &incoming[v] {
            queue.push(Candidate {
                index: i,
                score: odds[i].1,
                chain: &chains[i],
            });
        }
    };
    for v in 0..vertices {
        if seedable[v] > 2 {
            enqueue(v, &mut queue);
            processed[v] = true;
        }
    }
    let mut good = vec![false; chains.len()];
    let mut used_source = vec![false; vertices];
    let mut variants = 0;
    while variants <= 100 {
        let Some(candidate) = queue.pop() else {
            break;
        };
        let i = candidate.index;
        if good[i] {
            continue;
        }
        good[i] = true;
        let c = &chains[i];
        let new_variant = used_source[c.source];
        used_source[c.source] = true;
        if new_variant {
            variants += 1;
        }
        if new_variant && variants > 100 {
            continue;
        }
        for v in [c.source, c.target] {
            if !processed[v] {
                enqueue(v, &mut queue);
                processed[v] = true;
            }
        }
    }
    (0..chains.len()).filter(|&i| !good[i]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unmodified_jar_small_graphs_and_likelihood_bits() {
        let mut chains = Vec::new();
        let mut expected = Vec::new();
        let check = |chains: &[PruningChain], expected: &[usize]| {
            assert_eq!(chains_to_remove(chains), expected, "chains={chains:?}");
        };
        for line in include_str!("pruner_oracle.tsv").lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "CASE" => {
                    if !chains.is_empty() {
                        check(&chains, &expected);
                    }
                    chains.clear();
                    expected.clear();
                }
                "CHAIN" => {
                    if f[12] == "true" {
                        expected.push(chains.len());
                    }
                    chains.push(PruningChain {
                        source: f[1].parse().unwrap(),
                        target: f[2].parse().unwrap(),
                        first_bases: f[3].as_bytes().to_vec(),
                        bases: f[4].as_bytes().to_vec(),
                        multiplicities: f[5].split(',').map(|n| n.parse().unwrap()).collect(),
                        first_is_ref: f[6].parse().unwrap(),
                        any_ref: f[7].parse().unwrap(),
                        left_total: f[8].parse().unwrap(),
                        right_total: f[9].parse().unwrap(),
                        is_source: f[10].parse().unwrap(),
                        is_sink: f[11].parse().unwrap(),
                    });
                }
                "LLR" => {
                    let actual = constant_error_log_likelihood_ratio(
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                        f[3].parse().unwrap(),
                    );
                    assert_eq!(
                        actual.to_bits(),
                        u64::from_str_radix(f[4], 16).unwrap(),
                        "{line}: {actual}"
                    );
                }
                _ => panic!("invalid oracle row"),
            }
        }
        check(&chains, &expected);
    }
}
