//! GATK 4.6.2.0 ContaminationSegmenter, KernelSegmenter and PersistenceOptimizer.
use crate::{
    math::{fastmath, stream_sum, svd::singular_value_decomposition},
    pileup_summary::PileupSummary,
};
#[derive(Clone)]
struct Cost {
    d: f64,
    w: Vec<f64>,
    v: f64,
    c: f64,
}
fn cost(start: usize, end: usize, z: &[Vec<f64>], diag: &[f64]) -> Cost {
    let mut d = diag[start];
    let mut w = z[start].clone();
    let mut v = stream_sum(w.iter().map(|x| x * x));
    let count = (end + z.len() - start) % z.len() + 1;
    for offset in 1..count {
        let i = (start + offset) % z.len();
        d += diag[i];
        let mut dot = 0.0;
        for (wj, zj) in w.iter_mut().zip(&z[i]) {
            dot += zj * *wj;
            *wj += zj;
        }
        v += 2.0 * dot + diag[i];
    }
    Cost {
        d,
        w,
        v,
        c: d - v / count as f64,
    }
}
fn shift(c: &mut Cost, remove: usize, add: usize, z: &[Vec<f64>], diag: &[f64]) {
    c.d -= diag[remove];
    let mut dot = 0.0;
    for (wj, zj) in c.w.iter_mut().zip(&z[remove]) {
        dot += zj * *wj;
        *wj -= zj;
    }
    c.v += -2.0 * dot + diag[remove];
    c.d += diag[add];
    let mut dot = 0.0;
    for (wj, zj) in c.w.iter_mut().zip(&z[add]) {
        dot += zj * *wj;
        *wj += zj;
    }
    c.v += 2.0 * dot + diag[add];
}
fn window_costs(z: &[Vec<f64>], diag: &[f64]) -> Vec<f64> {
    let n = z.len();
    let mut start = (n - 49) % n;
    let mut end = 50 % n;
    let mut left = cost(start, 0, z, diag);
    let mut right = cost(1, end, z, diag);
    let mut total = cost(start, end, z, diag);
    let mut out = vec![0.0; n];
    out[0] = left.c + right.c - total.c;
    let reciprocal = 1.0 / 50.0;
    for center in 0..n {
        let next = (center + 1) % n;
        let en = (end + 1) % n;
        shift(&mut left, start, next, z, diag);
        left.c = left.d - left.v * reciprocal;
        shift(&mut right, next, en, z, diag);
        right.c = right.d - right.v * reciprocal;
        shift(&mut total, start, en, z, diag);
        total.c = total.d - 0.5 * total.v * reciprocal;
        out[next] = left.c + right.c - total.c;
        start = (start + 1) % n;
        end = en;
    }
    out
}
#[derive(Clone)]
struct Component {
    left: usize,
    right: usize,
    min: usize,
}
pub(crate) fn persistence_optimizer(data: &[f64]) -> Vec<usize> {
    let mut indices: Vec<_> = (0..data.len()).collect();
    indices.sort_by(|&a, &b| data[a].total_cmp(&data[b]));
    let mut colors = vec![usize::MAX; data.len()];
    let mut comp: Vec<Component> = Vec::new();
    let mut pairs = Vec::new();
    for &i in &indices {
        let l = if i == 0 { usize::MAX } else { colors[i - 1] };
        let r = if i + 1 == data.len() {
            usize::MAX
        } else {
            colors[i + 1]
        };
        if l == usize::MAX && r == usize::MAX {
            colors[i] = comp.len();
            comp.push(Component {
                left: i,
                right: i,
                min: i,
            });
        } else if l == usize::MAX || r == usize::MAX {
            let c = l.min(r);
            if i + 1 == comp[c].left {
                comp[c].left = i;
            } else if i == comp[c].right + 1 {
                comp[c].right = i;
            }
            colors[i] = c;
        } else {
            let paired = if data[comp[r].min] < data[comp[l].min] {
                comp[l].min
            } else {
                comp[r].min
            };
            let (min, max) = if data[paired] > data[i] || (data[paired] == data[i] && paired > i) {
                (i, paired)
            } else {
                (paired, i)
            };
            pairs.push((min, data[max] - data[min]));
            let (keep, merge) = if data[comp[l].min] < data[comp[r].min]
                || (data[comp[l].min] == data[comp[r].min] && l < r)
            {
                (l, r)
            } else {
                (r, l)
            };
            colors[comp[merge].left] = keep;
            colors[comp[merge].right] = keep;
            if comp[keep].min > comp[merge].min {
                comp[keep].left = comp[merge].left;
            } else {
                comp[keep].right = comp[merge].right;
            }
            colors[i] = colors[i - 1];
        }
    }
    pairs.sort_by(|a, b| b.1.total_cmp(&a.1));
    std::iter::once(indices[0])
        .chain(pairs.into_iter().map(|p| p.0))
        .collect()
}
/// java.util.Random.nextInt(bound), rejection test uses wrapping signed arithmetic.
pub(crate) struct JavaRandom(u64);
impl JavaRandom {
    pub(crate) fn new(seed: u64) -> Self {
        Self((seed ^ 0x5deece66d) & ((1 << 48) - 1))
    }
    fn next(&mut self) -> i32 {
        self.0 = (self.0.wrapping_mul(0x5deece66d).wrapping_add(11)) & ((1 << 48) - 1);
        (self.0 >> 17) as i32
    }
    pub(crate) fn next_int(&mut self, bound: usize) -> usize {
        let b = bound as i32;
        if b.isolate_lowest_one() == b {
            return ((i64::from(b) * i64::from(self.next())) >> 31) as usize;
        }
        loop {
            let bits = self.next();
            let val = bits % b;
            if bits.wrapping_sub(val).wrapping_add(b - 1) >= 0 {
                return val as usize;
            }
        }
    }
}
fn observations(data: &[f64]) -> Vec<Vec<f64>> {
    let p = data.len().min(100);
    let mut rng = JavaRandom::new(1216);
    let sample: Vec<_> = if p == data.len() {
        data.to_vec()
    } else {
        (0..p).map(|_| data[rng.next_int(data.len())]).collect()
    };
    let kernel = |a: f64, b: f64| {
        let d = a - b;
        fastmath::exp(-(d * d) / (2.0 * 0.025))
    };
    let mut matrix = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..i {
            let x = kernel(sample[i], sample[j]);
            matrix[i][j] = x;
            matrix[j][i] = x;
        }
        matrix[i][i] = kernel(sample[i], sample[i]);
    }
    let (s, mut u) = singular_value_decomposition(matrix);
    for row in &mut u {
        for (v, s) in row.iter_mut().zip(&s) {
            *v *= 1.0 / (s.sqrt() + 1e-10);
        }
    }
    data.iter()
        .map(|&d| {
            let k: Vec<_> = sample.iter().map(|&s| kernel(d, s)).collect();
            (0..p)
                .map(|j| {
                    let mut sum = 0.0;
                    for i in 0..p {
                        sum += k[i] * u[i][j];
                    }
                    sum
                })
                .collect()
        })
        .collect()
}
pub(crate) fn kernel_segmenter(data: &[f64]) -> Vec<usize> {
    if data.len() < 100 {
        return Vec::new();
    }
    let z = observations(data);
    // ArrayRealVector.getNorm is a plain sequential sum, followed by sqrt;
    // upstream then squares that rounded sqrt. Do not replace with sum of squares.
    let diag: Vec<_> = z
        .iter()
        .map(|r| {
            let norm = r.iter().map(|x| x * x).sum::<f64>().sqrt();
            norm * norm
        })
        .collect();
    let costs = window_costs(&z, &diag);
    let mut candidates: Vec<_> = persistence_optimizer(&costs)
        .into_iter()
        .filter(|&i| i != 0 && i != data.len() - 1)
        .take(10)
        .collect();
    candidates.sort_unstable();
    candidates.dedup();
    let mut start = 0;
    let mut segments = Vec::new();
    for end in candidates
        .into_iter()
        .chain(std::iter::once(data.len() - 1))
    {
        segments.push((start, end, cost(start, end, &z, &diag).c));
        start = end + 1;
    }
    let mut totals = vec![stream_sum(segments.iter().map(|s| s.2))];
    let mut changes = Vec::new();
    // Recomputing the adjacent pair's sequential cost is arithmetic-identical to
    // the upstream cached value (its inputs never change).
    while segments.len() > 1 {
        let merged: Vec<_> = segments
            .windows(2)
            .map(|s| cost(s[0].0, s[1].1, &z, &diag).c)
            .collect();
        let mut best = 0;
        for i in 1..merged.len() {
            if segments[i].2 + segments[i + 1].2 - merged[i]
                > segments[best].2 + segments[best + 1].2 - merged[best]
            {
                best = i;
            }
        }
        let left = segments[best];
        let right = segments.remove(best + 1);
        segments[best] = (left.0, right.1, merged[best]);
        changes.insert(0, left.1);
        totals.insert(0, stream_sum(segments.iter().map(|s| s.2)));
    }
    let mut best = 0;
    let mut value = totals[0];
    for (i, &total) in totals.iter().enumerate().skip(1) {
        let k = i as f64;
        let v = total + k + k * crate::math::java_log::log(data.len() as f64 / (k + 1e-10));
        if v < value {
            best = i;
            value = v;
        }
    }
    changes.truncate(best);
    changes.sort_unstable();
    changes
}
pub(crate) fn find_segments(sites: &[PileupSummary]) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    for group in crate::java_hashmap::group_indices(sites.iter().map(|s| s.contig.as_str())) {
        let het: Vec<_> = group
            .iter()
            .copied()
            .filter(|&i| (0.1..=0.9).contains(&sites[i].alt_fraction()))
            .collect();
        if het.is_empty() {
            continue;
        }
        let mafs: Vec<_> = het
            .iter()
            .map(|&i| {
                let f = sites[i].alt_fraction();
                f.min(1.0 - f)
            })
            .collect();
        let mut start = 0;
        for end in kernel_segmenter(&mafs)
            .into_iter()
            .chain(std::iter::once(het.len() - 1))
        {
            let mut segment: Vec<_> = group
                .iter()
                .copied()
                .filter(|&i| {
                    sites[i].position >= sites[het[start]].position
                        && sites[i].position <= sites[het[end]].position
                })
                .collect();
            segment.sort_by_key(|&i| sites[i].position);
            result.push(segment);
            start = end + 1;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn java_random_and_persistence_oracles() {
        // java.util.Random(1216), continuous nextInt calls. Last bound exercises rejection.
        let expected = [
            [1, 0, 0, 1, 0, 0, 1, 0],
            [85, 11, 97, 55, 53, 76, 37, 68],
            [304, 597, 940, 326, 870, 426, 472, 629],
            [
                10696943, 324157366, 623700311, 196473308, 53717894, 258552828, 943213532,
                100750171,
            ],
        ];
        let mut random = JavaRandom::new(1216);
        for (n, e) in [2, 100, 1000, 1073741825].into_iter().zip(expected) {
            assert_eq!(std::array::from_fn::<_, 8, _>(|_| random.next_int(n)), e);
        }
        // Unmodified GATK PersistenceOptimizer Java source, including equal minima.
        for (data, e) in [
            (vec![1., 0., 1., 0., 1.], vec![1, 3]),
            (vec![3., 1., 2., 1., 3., 0., 4.], vec![5, 1, 3]),
            (vec![0., 0., 0.], vec![0]),
            (vec![1.], vec![0]),
        ] {
            assert_eq!(persistence_optimizer(&data), e);
        }
    }
    #[test]
    fn recurrence_and_no_candidates() {
        let z = vec![vec![1., 2.], vec![2., 0.], vec![0., 3.]];
        let diag = vec![5., 4., 9.];
        let c = cost(0, 2, &z, &diag);
        assert_eq!(
            (c.d, c.w, c.v, c.c),
            (18., vec![3., 5.], 34., 18. - 34. / 3.)
        );
        assert!(kernel_segmenter(&vec![0.4; 99]).is_empty());
        assert!(kernel_segmenter(&vec![0.4; 101]).is_empty());
    }
}
