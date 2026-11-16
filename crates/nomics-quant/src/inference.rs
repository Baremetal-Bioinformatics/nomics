//! Range-factorized equivalence classes and data-fixed VBEM reductions.
//! Set NOMICS_QUANT_INFERENCE_TIMING to print inference wall time (including
//! packing and the final M-step) and accumulated per-iteration phase times.
use std::collections::BTreeMap;

use nomics_mem::HostArray;
use rayon::prelude::*;
use statrs::function::gamma::digamma;
use std::io;

const SCALE: f64 = (1u64 << 40) as f64;

#[derive(Default)]
pub(crate) struct Builder(BTreeMap<(Vec<usize>, Vec<u32>), Accumulator>);

struct Accumulator {
    count: u64,
    weights: Vec<u128>,
}

pub(crate) struct Class {
    pub tids: Vec<usize>,
    pub count: u64,
    pub weights: Vec<f64>,
}

impl Builder {
    pub fn add(&mut self, tids: Vec<usize>, weights: Vec<f64>) {
        let ranges = (tids.len() as f64).sqrt() as u32 + 4;
        let bins = weights.iter().map(|w| (w * ranges as f64) as u32).collect();
        let acc = self.0.entry((tids, bins)).or_insert_with(|| Accumulator {
            count: 0,
            weights: vec![0; weights.len()],
        });
        acc.count += 1;
        for (sum, w) in acc.weights.iter_mut().zip(weights) {
            *sum += (w.max(0.0) * SCALE).round() as u128;
        }
    }

    /// Classes are ordered by key and accumulate only integers (count and
    /// fixed-point weights), so merging per-thread builders by key is exact.
    pub fn merge(mut self, other: Builder) -> Builder {
        if self.0.len() < other.0.len() {
            return other.merge(self);
        }
        for (key, acc) in other.0 {
            match self.0.entry(key) {
                std::collections::btree_map::Entry::Vacant(e) => {
                    e.insert(acc);
                }
                std::collections::btree_map::Entry::Occupied(mut e) => {
                    let mine = e.get_mut();
                    mine.count += acc.count;
                    for (sum, w) in mine.weights.iter_mut().zip(acc.weights) {
                        *sum += w;
                    }
                }
            }
        }
        self
    }

    pub fn finish(self, effective: &[f64]) -> Vec<Class> {
        self.0
            .into_iter()
            .map(|((tids, _), acc)| {
                let mut weights: Vec<f64> = acc
                    .weights
                    .iter()
                    .zip(&tids)
                    .map(|(&w, &t)| (w as f64 / SCALE) / effective[t].max(1.0))
                    .collect();
                let total: f64 = weights.iter().sum();
                if total > 0.0 {
                    for w in &mut weights {
                        *w /= total;
                    }
                }
                Class {
                    tids,
                    weights,
                    count: acc.count,
                }
            })
            .collect()
    }
}

/// The class order, incidence order, and logical shard boundaries are frozen.
/// Only independent transcripts may be scheduled in a different order.
struct Packed {
    offsets: HostArray<usize>,
    counts: HostArray<f64>,
    tids: HostArray<usize>,
    weights: HostArray<f64>,
    /// Dense inference ID -> original BAM transcript ID, in ascending order.
    active: HostArray<usize>,
}

impl Packed {
    fn new(classes: &[Class], n: usize) -> io::Result<Self> {
        let mut used = HostArray::<usize>::zeroed(n)?;
        for c in classes {
            for &t in &c.tids {
                used[t] = 1;
            }
        }
        let mut active = HostArray::zeroed(used.iter().sum())?;
        let mut rank = 0;
        for (t, used) in used.iter_mut().enumerate() {
            if *used != 0 {
                active[rank] = t;
                *used = rank;
                rank += 1;
            }
        }
        let mut offsets = HostArray::<usize>::zeroed(classes.len() + 1)?;
        for (i, c) in classes.iter().enumerate() {
            offsets[i + 1] = offsets[i] + c.tids.len();
        }
        let entries = offsets[classes.len()];
        let mut tids = HostArray::zeroed(entries)?;
        let mut weights = HostArray::zeroed(entries)?;
        for (i, c) in classes.iter().enumerate() {
            for (dst, &t) in tids[offsets[i]..offsets[i + 1]].iter_mut().zip(&c.tids) {
                *dst = used[t];
            }
            weights[offsets[i]..offsets[i + 1]].copy_from_slice(&c.weights);
        }
        Ok(Self {
            offsets,
            counts: HostArray::collect_exact(classes.iter().map(|c| c.count as f64))?,
            tids,
            weights,
            active,
        })
    }
}

struct Shard {
    start: usize,
    end: usize,
    touched: HostArray<usize>,
    local_indices: HostArray<usize>,
    sums: HostArray<f64>,
}

fn shard_plan(packed: &Packed) -> io::Result<Vec<Shard>> {
    let offsets = &packed.offsets;
    let classes = packed.counts.len();
    let entries = packed.tids.len();
    let count = entries.div_ceil(4096).clamp(1, 128.min(classes.max(1)));
    let mut boundaries = vec![0];
    for i in 1..count {
        let target = (entries as u128 * i as u128).div_ceil(count as u128) as usize;
        let low = boundaries[i - 1] + 1;
        let high = classes - (count - i);
        let upper = (low + offsets[low..=high].partition_point(|&v| v < target)).min(high);
        let lower = upper.saturating_sub(1).max(low);
        boundaries.push(
            if offsets[lower].abs_diff(target) <= offsets[upper].abs_diff(target) {
                lower
            } else {
                upper
            },
        );
    }
    boundaries.push(classes);
    boundaries
        .windows(2)
        .map(|range| {
            let (start, end) = (range[0], range[1]);
            let tids = &packed.tids[offsets[start]..offsets[end]];
            let mut sorted = HostArray::collect_exact(tids.iter().copied())?;
            sorted.sort_unstable();
            let unique = sorted.chunk_by(|a, b| a == b).count();
            let mut touched = HostArray::zeroed(unique)?;
            for (dst, group) in touched.iter_mut().zip(sorted.chunk_by(|a, b| a == b)) {
                *dst = group[0];
            }
            let local_indices =
                HostArray::collect_exact(tids.iter().map(|t| touched.binary_search(t).unwrap()))?;
            Ok(Shard {
                start,
                end,
                sums: HostArray::zeroed(touched.len())?,
                touched,
                local_indices,
            })
        })
        .collect()
}

/// Transpose the shard merge, preserving ascending shard order separately for
/// every transcript. Each output has one writer and starts with the same +0.0.
struct Merge {
    offsets: HostArray<usize>,
    shards: HostArray<usize>,
    locals: HostArray<usize>,
}

impl Merge {
    fn new(shards: &[Shard], active: &[usize], n: usize) -> io::Result<Self> {
        let mut offsets = HostArray::<usize>::zeroed(n + 1)?;
        for shard in shards {
            for &t in &shard.touched {
                offsets[active[t] + 1] += 1;
            }
        }
        for t in 0..n {
            offsets[t + 1] += offsets[t];
        }
        let entries = offsets[n];
        let mut shard_ids = HostArray::zeroed(entries)?;
        let mut locals = HostArray::zeroed(entries)?;
        let mut cursor = HostArray::collect_exact(offsets[..n].iter().copied())?;
        for (s, shard) in shards.iter().enumerate() {
            for (local, &t) in shard.touched.iter().enumerate() {
                let slot = &mut cursor[active[t]];
                shard_ids[*slot] = s;
                locals[*slot] = local;
                *slot += 1;
            }
        }
        Ok(Self {
            offsets,
            shards: shard_ids,
            locals,
        })
    }

    fn apply(&self, shards: &[Shard], next: &mut [f64]) {
        next.par_iter_mut()
            .enumerate()
            .with_min_len(1024)
            .for_each(|(t, dst)| {
                let mut sum = 0.0;
                for i in self.offsets[t]..self.offsets[t + 1] {
                    sum += shards[self.shards[i]].sums[self.locals[i]];
                }
                *dst = sum;
            });
    }
}

pub(crate) fn optimize(classes: &[Class], n: usize) -> io::Result<(HostArray<f64>, u32, bool)> {
    optimize_observed(classes, n, |_, _| {})
}

fn optimize_observed(
    classes: &[Class],
    n: usize,
    mut observe: impl FnMut(f64, &[f64]),
) -> io::Result<(HostArray<f64>, u32, bool)> {
    let timer =
        std::env::var_os("NOMICS_QUANT_INFERENCE_TIMING").map(|_| std::time::Instant::now());
    let total: u64 = classes.iter().map(|c| c.count).sum();
    let mut alpha = HostArray::zeroed(n)?;
    alpha.fill(total as f64 / n as f64);
    let mut next = HostArray::<f64>::zeroed(n)?;
    let packed = Packed::new(classes, n)?;
    let mut theta = HostArray::<f64>::zeroed(packed.active.len())?;
    // digamma has no iteration-dependent inputs other than alpha. Reuse its
    // bits only when the argument's input bits are unchanged. exp still sees
    // this iteration's norm; no sums, products, or divisions are reassociated.
    let mut previous = HostArray::<f64>::zeroed(theta.len())?;
    previous.fill(f64::NAN);
    let mut psi = HostArray::<f64>::zeroed(theta.len())?;
    let mut shards = shard_plan(&packed)?;
    let merge = Merge::new(&shards, &packed.active, n)?;
    let mut iterations = 0;
    let mut converged = false;
    let mut phases = [std::time::Duration::ZERO; 5];
    while iterations < 10_000 {
        let phase = timer.map(|_| std::time::Instant::now());
        // Includes unused transcripts and the prior in precisely the original
        // transcript order. Even a run of identical 0.01 terms cannot be grouped.
        let norm = digamma(alpha.iter().map(|a| a + 0.01).sum());
        if let Some(p) = phase {
            phases[0] += p.elapsed();
        }
        let phase = timer.map(|_| std::time::Instant::now());
        theta
            .par_iter_mut()
            .zip(previous.par_iter_mut())
            .zip(psi.par_iter_mut())
            .zip(packed.active.par_iter())
            .with_min_len(256)
            .for_each(|(((dst, previous), psi), &t)| {
                let a = alpha[t];
                *dst = if a + 0.01 > 1e-10 {
                    if a.to_bits() != previous.to_bits() {
                        *psi = digamma(a + 0.01);
                        *previous = a;
                    }
                    (*psi - norm).exp()
                } else {
                    0.0
                };
            });
        if let Some(p) = phase {
            phases[1] += p.elapsed();
        }
        let phase = timer.map(|_| std::time::Instant::now());
        shards.par_iter_mut().for_each(|shard| {
            shard.sums.fill(0.0);
            let base = packed.offsets[shard.start];
            for i in shard.start..shard.end {
                let start = packed.offsets[i];
                let end = packed.offsets[i + 1];
                let tids = &packed.tids[start..end];
                let weights = &packed.weights[start..end];
                let locals = &shard.local_indices[start - base..end - base];
                if tids.len() == 1 {
                    shard.sums[locals[0]] += packed.counts[i];
                    continue;
                }
                let denominator: f64 = tids
                    .iter()
                    .zip(weights)
                    .filter(|(t, _)| theta[**t] > 0.0)
                    .map(|(&t, &w)| theta[t] * w)
                    .sum();
                if denominator > f64::MIN_POSITIVE {
                    let factor = packed.counts[i] / denominator;
                    for ((&t, &w), &local) in tids.iter().zip(weights).zip(locals) {
                        if theta[t] > 0.0 {
                            shard.sums[local] += theta[t] * w * factor;
                        }
                    }
                }
            }
        });
        if let Some(p) = phase {
            phases[2] += p.elapsed();
        }
        let phase = timer.map(|_| std::time::Instant::now());
        merge.apply(&shards, &mut next);
        if let Some(p) = phase {
            phases[3] += p.elapsed();
        }
        observe(norm, &next);
        let phase = timer.map(|_| std::time::Instant::now());
        iterations += 1;
        if iterations >= 50 {
            let difference = alpha
                .iter()
                .zip(&next)
                .filter(|(a, b)| **a > 0.01 && **b > 0.0)
                .map(|(a, b)| (b - a).abs() / b)
                .fold(f64::NEG_INFINITY, f64::max);
            converged = difference.is_finite() && difference < 0.01;
        }
        if let Some(p) = phase {
            phases[4] += p.elapsed();
        }
        std::mem::swap(&mut alpha, &mut next);
        if converged {
            break;
        }
    }
    // The last M-step retains the original sequential class/incidence order.
    // Transcripts absent from all classes never contribute to this calculation.
    theta
        .par_iter_mut()
        .zip(packed.active.par_iter())
        .with_min_len(256)
        .for_each(|(dst, &t)| {
            *dst = if alpha[t] < 1e-8 {
                0.0
            } else {
                digamma(alpha[t] + 0.01).exp()
            };
        });
    next.fill(0.0);
    for i in 0..packed.counts.len() {
        let start = packed.offsets[i];
        let end = packed.offsets[i + 1];
        let tids = &packed.tids[start..end];
        let weights = &packed.weights[start..end];
        if tids.len() == 1 {
            let t = packed.active[tids[0]];
            if alpha[t] >= 1e-8 {
                next[t] += packed.counts[i];
            }
        } else {
            let denominator: f64 = tids.iter().zip(weights).map(|(&t, &w)| theta[t] * w).sum();
            if denominator > f64::MIN_POSITIVE {
                let factor = packed.counts[i] / denominator;
                for (&t, &w) in tids.iter().zip(weights) {
                    let value = theta[t] * w;
                    if value > 0.0 {
                        next[packed.active[t]] += value * factor;
                    }
                }
            }
        }
    }
    if let Some(start) = timer {
        eprintln!(
            "nomics-quant inference: {:.6}s, {iterations} iterations, {} shards, {n} transcripts, {} active",
            start.elapsed().as_secs_f64(),
            shards.len(),
            packed.active.len()
        );
        eprintln!(
            "nomics-quant inference phases: norm={:.6}s theta={:.6}s shards={:.6}s merge={:.6}s convergence={:.6}s",
            phases[0].as_secs_f64(),
            phases[1].as_secs_f64(),
            phases[2].as_secs_f64(),
            phases[3].as_secs_f64(),
            phases[4].as_secs_f64()
        );
    }
    Ok((next, iterations, converged))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ambiguous_fragments_follow_unique_evidence() {
        let classes = vec![
            Class {
                tids: vec![0],
                weights: vec![1.0],
                count: 10,
            },
            Class {
                tids: vec![0, 1],
                weights: vec![0.5, 0.5],
                count: 100,
            },
            Class {
                tids: vec![1],
                weights: vec![1.0],
                count: 90,
            },
        ];
        let (counts, iterations, converged) = optimize(&classes, 3).unwrap();
        assert!(converged && iterations >= 50);
        assert_eq!(counts[2], 0.0);
        assert!((counts.iter().sum::<f64>() - 200.0).abs() < 1e-10);
        assert!(counts[0] < 21.0 && counts[1] > 179.0);
    }
}

#[cfg(test)]
mod determinism_tests {
    use super::*;

    #[test]
    fn multiple_shards_are_identical_at_one_and_six_threads() {
        let mut classes = Vec::new();
        for i in 0..5000 {
            classes.push(Class {
                tids: vec![i % 257, 257 + i % 251],
                weights: vec![0.3, 0.7],
                count: (1 + i % 13) as u64,
            });
        }
        classes.sort_by(|a, b| a.tids.cmp(&b.tids));
        assert!(
            shard_plan(&Packed::new(&classes, 508).unwrap())
                .unwrap()
                .len()
                > 1
        );
        let estimate = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| optimize(&classes, 508).unwrap())
        };
        let (one, it1, conv1) = estimate(1);
        let (six, it6, conv6) = estimate(6);
        assert_eq!(it1, it6);
        assert_eq!(conv1, conv6);
        assert_eq!(
            one.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            six.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
        let total: u64 = classes.iter().map(|c| c.count).sum();
        assert!((one.iter().sum::<f64>() - total as f64).abs() < 1e-8);
    }

    #[test]
    fn fixed_point_classes_do_not_depend_on_fragment_order() {
        let weights: Vec<Vec<f64>> = (1..1000)
            .map(|i| {
                let a = 0.3 + i as f64 * 0.00001;
                vec![a, 1.0 - a]
            })
            .collect();
        let collapse = |reverse| {
            let mut b = Builder::default();
            let mut order: Vec<_> = weights.iter().collect();
            if reverse {
                order.reverse();
            }
            for w in order {
                b.add(vec![0, 1], w.clone());
            }
            b.finish(&[100.0, 200.0])
        };
        let a = collapse(false);
        let b = collapse(true);
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(&b) {
            assert_eq!(a.count, b.count);
            assert_eq!(
                a.weights.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                b.weights.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
        }
    }
}

#[cfg(test)]
mod reference;

#[cfg(test)]
mod equivalence_tests {
    use super::*;

    fn check(classes: &[Class], n: usize) {
        let packed = Packed::new(classes, n).unwrap();
        let shards = shard_plan(&packed).unwrap();
        assert_eq!(
            shards.iter().map(|s| (s.start, s.end)).collect::<Vec<_>>(),
            reference::boundaries(classes),
            "logical shard boundaries must never depend on worker count or layout"
        );
        let mut trace = Vec::new();
        let (expected, iterations, converged) = reference::optimize(classes, n, |norm, next| {
            trace.push((
                norm.to_bits(),
                next.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            ));
        })
        .unwrap();
        for threads in [1, 6, 20] {
            let mut step = 0;
            let (actual, actual_iterations, actual_converged) = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    optimize_observed(classes, n, |norm, next| {
                        assert_eq!(norm.to_bits(), trace[step].0, "norm at iteration {step}");
                        for (t, (value, &expected)) in next.iter().zip(&trace[step].1).enumerate() {
                            assert_eq!(
                                value.to_bits(),
                                expected,
                                "iteration {step}, transcript {t}, {threads} threads"
                            );
                        }
                        step += 1;
                    })
                    .unwrap()
                });
            assert_eq!(step, trace.len());
            assert_eq!(
                (actual_iterations, actual_converged),
                (iterations, converged)
            );
            assert_eq!(
                actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                "final sequential masked M-step"
            );
        }
    }

    #[test]
    fn every_iteration_matches_original_with_sparse_and_tiny_transcripts() {
        let mut classes = Vec::new();
        let mut state = 0x428a_2f98_d728_ae22u64;
        for i in 0..1800 {
            let len = 1 + i % 23;
            let mut tids = Vec::new();
            let mut weights = Vec::new();
            for j in 0..len {
                // Nonnegative values spanning zero, subnormals, and unequal
                // weights exercise filters and the final transcript mask.
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                tids.push(((i * 17 + j * 31) % 1021) * 5);
                weights.push(match state % 13 {
                    0 => 0.0,
                    1 => 1e-320,
                    2 => 1e-20,
                    _ => (1 + state % 997) as f64 / 997.0,
                });
            }
            classes.push(Class {
                tids,
                weights,
                count: 1 + state % 10_000,
            });
        }
        // Identical transcript lists in separate classes must not be combined.
        classes.push(Class {
            tids: vec![5999],
            weights: vec![1.0],
            count: 9007199254740993,
        });
        classes.push(Class {
            tids: vec![5999],
            weights: vec![1.0],
            count: 3,
        });
        // An ambiguous transcript that falls below the last M-step's mask.
        classes.push(Class {
            tids: vec![0, 5998],
            weights: vec![1.0, 1e-100],
            count: 23,
        });
        check(&classes, 6000);
    }

    #[test]
    fn capped_shards_and_one_large_class_keep_the_original_addition_order() {
        let classes: Vec<_> = (0..520)
            .map(|i| Class {
                tids: (0..1024).map(|j| (i + j) % 1024).collect(),
                weights: (0..1024)
                    .map(|j| (1 + (i + j) % 7) as f64 / 4096.0)
                    .collect(),
                count: (1 + i % 13) as u64,
            })
            .collect();
        assert_eq!(
            shard_plan(&Packed::new(&classes, 1024).unwrap())
                .unwrap()
                .len(),
            128
        );
        check(&classes, 1024);
        check(&classes[..1], 1024);
    }

    #[test]
    fn empty_plan_and_nearest_boundary_ties_match_original() {
        for lengths in [
            vec![],
            vec![2048, 4096, 2048],
            vec![0, 8192, 0, 0],
            vec![1, 1, 1],
        ] {
            let classes: Vec<_> = lengths
                .into_iter()
                .map(|len| Class {
                    tids: vec![0; len],
                    weights: vec![1.0; len],
                    count: 1,
                })
                .collect();
            let packed = Packed::new(&classes, 1).unwrap();
            let shards = shard_plan(&packed).unwrap();
            assert_eq!(
                shards.iter().map(|s| (s.start, s.end)).collect::<Vec<_>>(),
                reference::boundaries(&classes)
            );
        }
    }
}
