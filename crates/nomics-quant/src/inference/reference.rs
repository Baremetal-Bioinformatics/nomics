//! Frozen pre-optimization arithmetic oracle (5e1d349), test-only.
use super::*;

struct Shard {
    start: usize,
    end: usize,
    touched: Vec<usize>,
    local_indices: Vec<Vec<usize>>,
    sums: HostArray<f64>,
}

fn shard_plan(classes: &[Class]) -> io::Result<Vec<Shard>> {
    let mut offsets = vec![0usize];
    for c in classes {
        offsets.push(offsets.last().unwrap() + c.tids.len());
    }
    let entries = *offsets.last().unwrap();
    let count = entries
        .div_ceil(4096)
        .clamp(1, 128.min(classes.len().max(1)));
    let mut boundaries = vec![0];
    for i in 1..count {
        let target = (entries as u128 * i as u128).div_ceil(count as u128) as usize;
        let low = boundaries[i - 1] + 1;
        let high = classes.len() - (count - i);
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
    boundaries.push(classes.len());
    boundaries
        .windows(2)
        .map(|range| {
            let (start, end) = (range[0], range[1]);
            let mut touched: Vec<usize> = classes[start..end]
                .iter()
                .flat_map(|c| c.tids.iter().copied())
                .collect();
            touched.sort_unstable();
            touched.dedup();
            let local_indices = classes[start..end]
                .iter()
                .map(|c| {
                    c.tids
                        .iter()
                        .map(|t| touched.binary_search(t).unwrap())
                        .collect()
                })
                .collect();
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

pub(super) fn optimize(
    classes: &[Class],
    n: usize,
    mut observe: impl FnMut(f64, &[f64]),
) -> io::Result<(HostArray<f64>, u32, bool)> {
    let total: u64 = classes.iter().map(|c| c.count).sum();
    let mut alpha = HostArray::zeroed(n)?;
    alpha.fill(total as f64 / n as f64);
    let mut next = HostArray::<f64>::zeroed(n)?;
    let mut theta = HostArray::<f64>::zeroed(n)?;
    let mut shards = shard_plan(classes)?;
    let mut iterations = 0;
    let mut converged = false;
    while iterations < 10_000 {
        let norm = digamma(alpha.iter().map(|a| a + 0.01).sum());
        for (dst, &a) in theta.iter_mut().zip(alpha.iter()) {
            *dst = if a + 0.01 > 1e-10 {
                (digamma(a + 0.01) - norm).exp()
            } else {
                0.0
            };
        }
        shards.par_iter_mut().for_each(|shard| {
            shard.sums.fill(0.0);
            for (c, locals) in classes[shard.start..shard.end]
                .iter()
                .zip(&shard.local_indices)
            {
                if c.tids.len() == 1 {
                    shard.sums[locals[0]] += c.count as f64;
                    continue;
                }
                let denominator: f64 = c
                    .tids
                    .iter()
                    .zip(&c.weights)
                    .filter(|(t, _)| theta[**t] > 0.0)
                    .map(|(&t, &w)| theta[t] * w)
                    .sum();
                if denominator > f64::MIN_POSITIVE {
                    let factor = c.count as f64 / denominator;
                    for ((&t, &w), &local) in c.tids.iter().zip(&c.weights).zip(locals) {
                        if theta[t] > 0.0 {
                            shard.sums[local] += theta[t] * w * factor;
                        }
                    }
                }
            }
        });
        next.fill(0.0);
        for shard in &shards {
            for (&t, &value) in shard.touched.iter().zip(&shard.sums) {
                next[t] += value;
            }
        }
        observe(norm, &next);
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
        std::mem::swap(&mut alpha, &mut next);
        if converged {
            break;
        }
    }
    // The last M-step masks tiny transcripts and has no global normalization.
    let basis = HostArray::collect_exact(alpha.iter().map(|a| {
        if *a < 1e-8 {
            0.0
        } else {
            digamma(a + 0.01).exp()
        }
    }))?;
    next.fill(0.0);
    for c in classes {
        if c.tids.len() == 1 {
            let t = c.tids[0];
            if alpha[t] >= 1e-8 {
                next[t] += c.count as f64;
            }
        } else {
            let denominator: f64 = c
                .tids
                .iter()
                .zip(&c.weights)
                .map(|(&t, &w)| basis[t] * w)
                .sum();
            if denominator > f64::MIN_POSITIVE {
                let factor = c.count as f64 / denominator;
                for (&t, &w) in c.tids.iter().zip(&c.weights) {
                    let value = basis[t] * w;
                    if value > 0.0 {
                        next[t] += value * factor;
                    }
                }
            }
        }
    }
    Ok((next, iterations, converged))
}

pub(super) fn boundaries(classes: &[Class]) -> Vec<(usize, usize)> {
    shard_plan(classes)
        .unwrap()
        .iter()
        .map(|s| (s.start, s.end))
        .collect()
}
