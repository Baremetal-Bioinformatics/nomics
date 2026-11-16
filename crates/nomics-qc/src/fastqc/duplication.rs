use super::{contaminants, float::java, java_map};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write;
#[derive(Default)]
pub(super) struct Dup {
    seqs: HashMap<Key, (i32, usize)>,
    limit: u64,
    count: u64,
}
// All default duplication keys are <=50 ASCII bases. Store them inline so
// admitting a sequence does not allocate a per-read heap object.
#[derive(Hash, PartialEq, Eq)]
struct Key {
    bases: [u8; 50],
    len: u8,
}
impl Key {
    fn new(seq: &[u8]) -> Self {
        let len = seq.len().min(50);
        let mut bases = [0; 50];
        bases[..len].copy_from_slice(&seq[..len]);
        Self {
            bases,
            len: len as u8,
        }
    }
    fn seq(&self) -> &[u8] {
        &self.bases[..self.len as usize]
    }
}
fn spread(seq: &[u8]) -> u32 {
    let h = seq
        .iter()
        .fold(0u32, |h, &c| h.wrapping_mul(31).wrapping_add(c as u32));
    h ^ (h >> 16)
}
impl Dup {
    pub fn add(&mut self, seq: &[u8]) {
        self.count += 1;
        let key = Key::new(seq);
        let n = self.seqs.len();
        if let Some((count, _)) = self.seqs.get_mut(&key) {
            *count = count.wrapping_add(1);
        } else if n < 100000 {
            self.seqs.insert(key, (1, n));
        }
        if n < 100000 {
            self.limit = self.count;
        }
    }
    pub fn render(&self, total: u64) -> (String, u8, String, u8) {
        let mut insertion: Vec<_> = self.seqs.iter().collect();
        insertion.sort_by_key(|(_, (_, order))| *order);
        let entries: Vec<_> = insertion
            .iter()
            .map(|(s, _)| (spread(s.seq()), s.seq()))
            .collect();
        let mut seqs: Vec<_> = java_map::order(&entries)
            .into_iter()
            .map(|i| insertion[i])
            .collect();
        let mut collated = BTreeMap::<i32, (usize, usize)>::new();
        for &(_, &(count, _)) in &seqs {
            let n = collated.len();
            let v = collated.entry(count).or_insert((0, n));
            v.0 += 1;
        }
        let mut insertion: Vec<_> = collated.iter().collect();
        insertion.sort_by_key(|&(_, &(_, order))| order);
        let entries: Vec<_> = insertion
            .iter()
            .map(|&(&d, _)| ((d as u32) ^ (d as u32 >> 16), d))
            .collect();
        let collated_order: Vec<_> = java_map::order(&entries)
            .into_iter()
            .map(|i| insertion[i])
            .collect();
        // correctedCounts is another default HashMap populated in collated order.
        let entries: Vec<_> = collated_order
            .iter()
            .map(|&(&d, _)| ((d as u32) ^ (d as u32 >> 16), d))
            .collect();
        let levels: Vec<_> = java_map::order(&entries)
            .into_iter()
            .map(|i| collated_order[i])
            .collect();
        let mut dedup = 0.0;
        let mut raw = 0.0;
        let mut slots = [0.0; 16];
        for (&d, &(n, _)) in levels {
            let count = correct(self.limit, total, d, n as f64);
            dedup += count;
            raw += count * d as f64;
            let slot = match d {
                1..=9 => (d - 1) as usize,
                10..=50 => 9,
                51..=100 => 10,
                101..=500 => 11,
                501..=1000 => 12,
                1001..=5000 => 13,
                5001..=10000 => 14,
                _ => 15,
            };
            slots[slot] += count * d as f64;
        }
        let pct = if raw == 0.0 {
            100.0
        } else {
            (dedup / raw) * 100.0
        };
        let mut body = format!(
            "#Total Deduplicated Percentage\t{}\n#Duplication Level\tPercentage of total\n",
            java(pct)
        );
        let labels = [
            "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k",
            ">5k", ">10k+",
        ];
        for (label, n) in labels.iter().zip(slots) {
            writeln!(body, "{label}\t{}", java((n / raw) * 100.0)).unwrap();
        }
        // Arrays.sort(Object[]) is stable: count ties preserve HashMap iteration.
        seqs.sort_by_key(|(_, (count, _))| std::cmp::Reverse(*count));
        let mut over = String::new();
        let mut status = 0;
        for (seq, &(count, _)) in seqs {
            let pct = (count as f64 / total as f64) * 100.0;
            if pct <= 0.1 {
                continue;
            }
            if over.is_empty() {
                over.push_str("#Sequence\tCount\tPercentage\tPossible Source\n");
            }
            status = status.max(if pct > 1.0 { 2 } else { 1 });
            writeln!(
                over,
                "{}\t{count}\t{}\t{}",
                String::from_utf8_lossy(seq.seq()),
                java(pct),
                contaminants::find(seq.seq())
            )
            .unwrap();
        }
        (
            body,
            if pct < 50.0 {
                2
            } else if pct < 70.0 {
                1
            } else {
                0
            },
            over,
            status,
        )
    }
}
fn correct(limit: u64, total: u64, dup: i32, n: f64) -> f64 {
    if limit == total || total - (n as u64) < limit {
        return n;
    }
    let cutoff = 1.0 - n / (n + 0.01);
    let mut p = 1.0;
    for i in 0..limit {
        p *= ((total - i) as i64 - dup as i64) as f64 / (total - i) as f64;
        if p < cutoff {
            p = 0.0;
            break;
        }
    }
    n / (1.0 - p)
}

#[cfg(test)]
mod tests {
    #[test]
    fn keys_truncate_at_50_and_freeze_after_100000_distinct() {
        let mut d = super::Dup::default();
        let a = vec![b'A'; 50];
        let mut b = a.clone();
        b.push(b'C');
        d.add(&a);
        d.add(&b);
        assert_eq!(d.seqs.len(), 1);
        assert_eq!(d.seqs[&super::Key::new(&a)].0, 2);
        let mut key = [b'A'; 9];
        for n in 0..100000 {
            for (pos, v) in key.iter_mut().enumerate() {
                *v = b"ACGT"[(n >> (pos * 2)) & 3];
            }
            d.add(&key);
        }
        assert_eq!(d.seqs.len(), 100000);
        assert_eq!(d.limit, 100001);
        let limit = d.limit;
        let count = d.seqs[&super::Key::new(&a)].0;
        d.add(&a);
        d.add(b"NNNNNNNNN");
        assert_eq!(d.limit, limit);
        assert_eq!(d.seqs[&super::Key::new(&a)].0, count + 1);
        assert!(!d.seqs.contains_key(&super::Key::new(b"NNNNNNNNN")));
    }
}
