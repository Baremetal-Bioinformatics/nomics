use crate::{invalid, open_text};
use nomics_mem::HostArray;
use std::io::{self, BufRead};
use std::path::Path;

pub(crate) fn encode(seq: &[u8]) -> (u32, Vec<usize>) {
    let mut value = 0;
    let mut ns = Vec::new();
    for (i, &b) in seq.iter().enumerate() {
        value = (value << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    ns.push(i);
                    0
                }
            };
    }
    (value, ns)
}

pub(crate) fn decode(mut value: u32, len: usize) -> String {
    let mut seq = vec![b'A'; len];
    for b in seq.iter_mut().rev() {
        *b = b"ACGT"[(value & 3) as usize];
        value >>= 2;
    }
    String::from_utf8(seq).unwrap()
}

pub(crate) struct Whitelist {
    pub values: HostArray<u32>,
    pub len: usize,
    prefix: Vec<usize>,
    pub exact: HostArray<u32>,
}

#[derive(Clone, Debug)]
pub(crate) struct Match {
    // STAR's cbMatch: negative rejection, zero exact, positive candidate count.
    pub code: i32,
    pub candidates: Vec<(u32, u8)>,
}

impl Match {
    pub fn stat(&self) -> usize {
        match self.code {
            -2 => 3,
            -23 => 4,
            -24 => 5,
            -1 => 6,
            -3 => 8,
            0 => 9,
            1 => 10,
            _ => 11,
        }
    }
    pub fn resolve(&self, wl: &Whitelist) -> Result<u32, usize> {
        if self.code <= 1 {
            let cb = self.candidates[0].0;
            if self.code == 1 && wl.exact[cb as usize] == 0 {
                return Err(5);
            }
            return Ok(cb);
        }
        // SoloReadFeature_inputRecords: each product and each accumulation rounds to float.
        let (mut total, mut maximum, mut best) = (0f32, 0f32, 0);
        for &(cb, qual) in &self.candidates {
            let prior = wl.exact[cb as usize];
            if prior == 0 {
                continue;
            }
            let q = (i32::from(qual) - 33).min(33);
            let p = (f64::from(prior) * 10f64.powf(-f64::from(q) / 10.0)) as f32;
            total += p;
            if p > maximum {
                maximum = p;
                best = cb;
            }
        }
        if total > 0.0 && maximum >= 0.975f32 * total {
            Ok(best)
        } else {
            Err(4)
        }
    }
}

impl Whitelist {
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut n = 0;
        for line in open_text(path)?.lines() {
            line?;
            n += 1;
        }
        let mut values = HostArray::zeroed(n)?;
        for (i, line) in open_text(path)?.lines().enumerate() {
            let line = line?;
            if line.len() != 16 || !encode(line.as_bytes()).1.is_empty() {
                return Err(invalid("whitelist requires 16 A/C/G/T bases per line"));
            }
            values[i] = encode(line.as_bytes()).0;
        }
        values.sort_unstable();
        let mut len = 0;
        for i in 0..n {
            if len == 0 || values[i] != values[len - 1] {
                values[len] = values[i];
                len += 1;
            }
        }
        if len == 0 || len > (1 << 23) {
            return Err(invalid("whitelist size outside 1..=2^23"));
        }
        let mut prefix = vec![0; 65537];
        let mut at = 0;
        for (p, entry) in prefix.iter_mut().enumerate() {
            while at < len && (values[at] >> 16) < p as u32 {
                at += 1;
            }
            *entry = at;
        }
        Ok(Self {
            values,
            len,
            prefix,
            exact: HostArray::zeroed(len)?,
        })
    }
    fn find(&self, value: u32) -> Option<u32> {
        let p = (value >> 16) as usize;
        let start = self.prefix[p];
        self.values[start..self.prefix[p + 1]]
            .binary_search(&value)
            .ok()
            .map(|i| (start + i) as u32)
    }
    pub fn match_cb(&self, seq: &[u8], qual: &[u8]) -> Match {
        let (value, ns) = encode(seq);
        let mut m = Match {
            code: -1,
            candidates: Vec::new(),
        };
        if ns.len() > 1 {
            m.code = -2;
            return m;
        }
        if ns.is_empty() {
            if let Some(cb) = self.find(value) {
                m.code = 0;
                m.candidates.push((cb, 0));
                return m;
            }
            // Reverse base order, then XOR 1,2,3: posterior accumulation order matters.
            for i in 0..16 {
                for delta in 1..4 {
                    if let Some(cb) = self.find(value ^ (delta << (2 * i))) {
                        m.candidates.push((cb, qual[15 - i]));
                    }
                }
            }
        } else {
            let i = ns[0];
            for base in 0..4 {
                if let Some(cb) = self.find(value ^ (base << (2 * (15 - i)))) {
                    m.candidates.push((cb, qual[i]));
                }
            }
            if m.candidates.len() > 1 {
                m.code = -3;
                m.candidates.clear();
                return m;
            }
        }
        if !m.candidates.is_empty() {
            m.code = m.candidates.len() as i32;
        }
        m
    }
    pub fn match_read(&self, cb: &[u8], cq: &[u8], umi: &[u8]) -> Match {
        let (_, ns) = encode(umi);
        // UMI rejection takes precedence even when the CB itself is invalid.
        if !ns.is_empty() {
            return Match {
                code: -23,
                candidates: vec![],
            };
        }
        if umi.iter().all(|b| *b == umi[0]) {
            return Match {
                code: -24,
                candidates: vec![],
            };
        }
        self.match_cb(cb, cq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn wl(seqs: &[&str]) -> Whitelist {
        let p = std::env::temp_dir().join(format!(
            "solo-wl-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::write(&p, seqs.join("\n") + "\n").unwrap();
        let w = Whitelist::load(&p).unwrap();
        std::fs::remove_file(p).unwrap();
        w
    }
    #[test]
    fn precedence_and_n_ambiguity() {
        let w = wl(&["AAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAC"]);
        assert_eq!(w.match_cb(b"AAAAAAAAAAAAAAAN", &[73; 16]).code, -3);
        assert_eq!(
            w.match_read(b"NNNNNNNNNNNNNNNN", &[73; 16], b"AAAAAAAAAAAN")
                .code,
            -23
        );
        assert_eq!(
            w.match_read(b"NNNNNNNNNNNNNNNN", &[73; 16], b"CCCCCCCCCCCC")
                .code,
            -24
        );
    }
    #[test]
    fn posterior_and_exact_prior() {
        let mut w = wl(&["AAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAC"]);
        let m = w.match_cb(b"AAAAAAAAAAAAAAAG", &[73; 16]);
        assert_eq!(m.code, 2);
        assert_eq!(m.resolve(&w), Err(4));
        w.exact[0] = 39;
        w.exact[1] = 1;
        // Test against the float arithmetic, including the 0.975 boundary.
        let p0 = (39.0 * 10f64.powf(-3.3)) as f32;
        let p1 = 10f64.powf(-3.3) as f32;
        assert_eq!(m.resolve(&w).is_ok(), p0 >= 0.975f32 * (p0 + p1));
        w.exact[0] = 100;
        assert_eq!(m.resolve(&w), Ok(0));
        let one = w.match_cb(b"CAAAAAAAAAAAAAAA", &[73; 16]);
        assert_eq!(one.resolve(&w), Ok(0));
        w.exact[0] = 0;
        assert_eq!(one.resolve(&w), Err(5));
    }
}
