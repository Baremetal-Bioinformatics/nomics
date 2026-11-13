use std::io;

/// Preserve samtools' actual circular pileup semantics, including the full
/// flush's size-1 endpoint and realloc_buffers' byte-count memcpy.
pub(super) struct Ring {
    pub pos: i64,
    start: usize,
    data: Vec<i32>,
    pub cov: [u64; 1002],
}
impl Ring {
    pub fn new() -> Self {
        Self {
            pos: 0,
            start: 0,
            data: vec![0; 1500],
            cov: [0; 1002],
        }
    }
    fn index(&self, p: i64) -> usize {
        (self.start + (p - self.pos) as usize % self.data.len()) % self.data.len()
    }
    pub fn flush(&mut self, new: i64) -> io::Result<()> {
        if new == self.pos {
            return Ok(());
        }
        let p = if new == -1 || new - self.pos >= self.data.len() as i64 {
            self.pos + self.data.len() as i64 - 1
        } else {
            new
        };
        if p < self.pos {
            return Err(super::invalid("coverage coordinates out of order"));
        }
        let n = (p - self.pos) as usize;
        for i in 0..n {
            let idx = (self.start + i) % self.data.len();
            let depth = self.data[idx];
            if depth != 0 {
                self.cov[(depth as usize).min(1001)] += 1;
                self.data[idx] = 0;
            }
        }
        self.start = if new == -1 { 0 } else { self.index(p) };
        self.pos = new;
        Ok(())
    }
    pub fn insert(&mut self, from: i64, to: i64) -> io::Result<()> {
        if to - from > self.data.len() as i64 || from < self.pos {
            return Err(super::invalid("coverage interval outside ring"));
        }
        let mut a = self.index(from);
        let b = self.index(to);
        if a > b {
            for x in &mut self.data[a..] {
                *x = x.wrapping_add(1);
            }
            a = 0;
        }
        for x in &mut self.data[a..b] {
            *x = x.wrapping_add(1);
        }
        Ok(())
    }
    pub fn grow(&mut self, read_len: usize) {
        let n = self.data.len() - self.start;
        let mut bytes = vec![0u8; read_len * 5 * 4];
        let old: Vec<u8> = self.data.iter().flat_map(|v| v.to_ne_bytes()).collect();
        bytes[..n].copy_from_slice(&old[self.start * 4..self.start * 4 + n]);
        if self.start > 1 {
            bytes[n * 4..n * 4 + self.start].copy_from_slice(&old[..self.start]);
        }
        self.data = bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|s| i32::from_ne_bytes(*s))
            .collect();
        self.start = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn resize_reproduces_upstream_byte_counts() {
        let mut r = Ring::new();
        r.data.fill(0x01020304);
        r.grow(300);
        assert!(r.data[..375].iter().all(|&x| x == 0x01020304));
        assert!(r.data[375..].iter().all(|&x| x == 0));
        r.data.fill(0x01020304);
        r.start = 2;
        let old: Vec<_> = r.data.iter().flat_map(|x| x.to_ne_bytes()).collect();
        r.grow(600);
        let bytes: Vec<_> = r.data.iter().flat_map(|x| x.to_ne_bytes()).collect();
        assert_eq!(&bytes[..1498], &old[8..1506]);
        assert!(bytes[1498..5992].iter().all(|&x| x == 0));
        assert_eq!(&bytes[5992..5994], &old[..2]);
    }
    #[test]
    fn overlaps_wrap_and_overflow_bin() {
        let mut r = Ring::new();
        r.flush(1498).unwrap();
        r.insert(1498, 1502).unwrap();
        r.insert(1500, 1503).unwrap();
        r.flush(1503).unwrap();
        assert_eq!(r.cov[1], 3);
        assert_eq!(r.cov[2], 2);
        for _ in 0..1002 {
            r.insert(1503, 1504).unwrap();
        }
        r.flush(-1).unwrap();
        assert_eq!(r.cov[1001], 1);
    }
    #[test]
    fn upstream_full_flush_leaves_last_slot() {
        let mut r = Ring::new();
        r.insert(1499, 1500).unwrap();
        r.flush(-1).unwrap();
        assert_eq!(r.cov[1], 0);
    }
}
