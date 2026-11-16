use super::{above, float::java};
use std::collections::HashMap;
use std::fmt::Write;

pub(super) struct Gc {
    counts: [f64; 101],
    models: HashMap<usize, Vec<Vec<(usize, f64)>>>,
}
impl Default for Gc {
    fn default() -> Self {
        Self {
            counts: [0.0; 101],
            models: HashMap::new(),
        }
    }
}
impl Gc {
    pub fn add(&mut self, seq: &[u8]) {
        // PerSequenceGCContent.truncateSequence rounds DOWN, not to nearest.
        let n = seq.len();
        let n = if n > 1000 {
            n / 1000 * 1000
        } else if n > 100 {
            n / 100 * 100
        } else {
            n
        };
        if n == 0 {
            return;
        }
        let gc = seq[..n].iter().filter(|&&b| b == b'G' || b == b'C').count();
        let model = self.models.entry(n).or_insert_with(|| {
            let mut claims = [0u32; 101];
            let bounds = |pos: usize| {
                (
                    (((pos as f64 - 0.5).max(0.0) * 100.0 / n as f64) + 0.5).floor() as usize,
                    (((pos as f64 + 0.5).min(n as f64) * 100.0 / n as f64) + 0.5).floor() as usize,
                )
            };
            for pos in 0..=n {
                let (a, b) = bounds(pos);
                for c in &mut claims[a..=b] {
                    *c += 1;
                }
            }
            (0..=n)
                .map(|pos| {
                    let (a, b) = bounds(pos);
                    (a..=b).map(|p| (p, 1.0 / claims[p] as f64)).collect()
                })
                .collect()
        });
        for &(p, v) in &model[gc] {
            self.counts[p] += v;
        }
    }
    pub fn render(&self) -> (String, u8) {
        let mut total = 0.0;
        let mut mode_count = 0.0;
        let mut first = 0;
        for (i, &v) in self.counts.iter().enumerate() {
            total += v;
            if v > mode_count {
                mode_count = v;
                first = i;
            }
        }
        let cutoff = self.counts[first] - self.counts[first] / 10.0; // comment says 95%, code uses 90%.
        let mut mode = 0.0;
        let mut n = 0;
        let mut top = true;
        let mut bottom = true;
        for i in first..=100 {
            if self.counts[i] > cutoff {
                mode += i as f64;
                n += 1;
            } else {
                top = false;
                break;
            }
        }
        for i in (0..first).rev() {
            if self.counts[i] > cutoff {
                mode += i as f64;
                n += 1;
            } else {
                bottom = false;
                break;
            }
        }
        if top || bottom {
            mode = first as f64;
        } else {
            mode /= n as f64;
        }
        let mut stdev = 0.0;
        for (i, &v) in self.counts.iter().enumerate() {
            stdev += (i as f64 - mode).powi(2) * v;
        }
        stdev = (stdev / (total - 1.0)).sqrt();
        let lhs = 1.0 / (2.0 * std::f64::consts::PI * stdev * stdev).sqrt();
        let mut deviation = 0.0;
        let mut body = String::from("#GC Content\tCount\n");
        for (i, &v) in self.counts.iter().enumerate() {
            let rhs =
                std::f64::consts::E.powf(0.0 - ((i as f64 - mode).powi(2) / (2.0 * stdev * stdev)));
            deviation += (lhs * rhs * total - v).abs();
            writeln!(body, "{i}\t{}", java(v)).unwrap();
        }
        (body, above(deviation / total * 100.0, 15.0, 30.0))
    }
}
pub(super) fn lengths(counts: &[u64]) -> String {
    let min = counts
        .iter()
        .position(|&n| n > 0)
        .unwrap_or(0)
        .saturating_sub(1);
    let max = counts.iter().rposition(|&n| n > 0).unwrap_or(0) + 1;
    let mut base = 1;
    let interval = 'find: loop {
        for d in [1, 2, 5] {
            if (max - min) / (base * d) <= 50 {
                break 'find base * d;
            }
        }
        base *= 10;
    };
    let start = min / interval * interval;
    let mut bins = Vec::new();
    for a in (start..=max).step_by(interval) {
        let b = (a + interval - 1).min(max);
        let mut sum = 0.0;
        for i in a..=b {
            sum += counts.get(i).copied().unwrap_or(0) as f64;
        }
        bins.push((a, b, sum));
    }
    let mut body = String::from("#Length\tCount\n");
    for (i, &(a, b, n)) in bins.iter().enumerate() {
        if (i == 0 || i == bins.len() - 1) && n == 0.0 {
            continue;
        }
        let label = if interval == 1 {
            a.to_string()
        } else {
            format!("{a}-{b}")
        };
        writeln!(body, "{label}\t{}", java(n)).unwrap();
    }
    body
}
