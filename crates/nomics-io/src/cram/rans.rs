//! rANS 4x8 (CRAM 3.0 block method 4), orders 0 and 1.
//!
//! Spec: hts-specs `CRAMcodecs.tex`, "rANS 4x8". The encoder follows htscodecs
//! `rANS_static.c` (`rans_compress_O0`/`_O1`): the same frequency
//! normalisation (tables summing to 4095), the same table RLE and the same
//! interleaving of the four states, so every table it writes is one htslib's
//! decoder accepts. Byte identity with htscodecs is not a goal; the decoder
//! here exists to test the encoder and accepts any valid stream.

use std::io;

const TF_SHIFT: u32 = 12;
const TOTFREQ: u32 = 1 << TF_SHIFT;
/// Lower bound of the state interval.
const RANS_L: u32 = 1 << 23;

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("rANS: {msg}"))
}

#[derive(Clone, Copy, Default)]
struct EncSym {
    start: u32,
    freq: u32,
}

/// Bytes are produced back to front; `out` holds them reversed until
/// [`finish`] flips them into stream order.
#[inline]
fn put(x: &mut u32, out: &mut Vec<u8>, s: EncSym) {
    let x_max = ((RANS_L >> TF_SHIFT) << 8) * s.freq;
    while *x >= x_max {
        out.push(*x as u8);
        *x >>= 8;
    }
    *x = ((*x / s.freq) << TF_SHIFT) + (*x % s.freq) + s.start;
}

/// Flush the states (`states[0]` is read first by the decoder) and reverse.
fn finish(states: [u32; 4], mut rev: Vec<u8>) -> Vec<u8> {
    for &x in states.iter().rev() {
        rev.extend_from_slice(&x.to_be_bytes());
    }
    rev.reverse();
    rev
}

/// `n` as a frequency-table value: 1 byte below 128, else 2 (ITF8 form).
fn put_freq(t: &mut Vec<u8>, f: u32) {
    if f < 128 {
        t.push(f as u8);
    } else {
        t.push(0x80 | (f >> 8) as u8);
        t.push(f as u8);
    }
}

/// Normalise order-0 counts as `rans_compress_O0` does (sum 4095).
#[allow(clippy::needless_range_loop)] // indexed as in htscodecs
///
/// When the first pass overshoots by more than half the largest frequency,
/// htscodecs scales the already-normalised frequencies again by 0.98
/// (`goto normalise_harder`); each pass shrinks them, so this terminates.
fn normalise_o0(f: &mut [u32; 256], n: usize) {
    let mut tr: u64 = ((TOTFREQ as u64) << 31) / n as u64 + (1u64 << 30) / n as u64;
    loop {
        let (mut fsum, mut m, mut big) = (0u32, 0u32, 0usize);
        for j in 0..256 {
            if f[j] == 0 {
                continue;
            }
            if m < f[j] {
                m = f[j];
                big = j;
            }
            f[j] = ((f[j] as u64 * tr) >> 31) as u32;
            if f[j] == 0 {
                f[j] = 1;
            }
            fsum += f[j];
        }
        fsum += 1;
        if fsum < TOTFREQ {
            f[big] += TOTFREQ - fsum;
        } else if fsum - TOTFREQ > f[big] / 2 {
            tr = 2104533975; // ×0.98, as htscodecs
            continue;
        } else {
            f[big] -= fsum - TOTFREQ;
        }
        return;
    }
}

/// Normalise one order-1 context row as `rans_compress_O1` does.
#[allow(clippy::needless_range_loop)] // indexed as in htscodecs
/// Retries rescale the current frequencies by 0.98, as htscodecs does.
fn normalise_o1(f: &mut [u32; 256], total: u32) {
    let mut p = TOTFREQ as f64 / total as f64;
    loop {
        let (mut t2, mut m, mut big) = (0u32, 0u32, 0usize);
        for j in 0..256 {
            if f[j] == 0 {
                continue;
            }
            if m < f[j] {
                m = f[j];
                big = j;
            }
            // C: `F[i][j] *= p` on an int truncates the double product.
            f[j] = (f[j] as f64 * p) as u32;
            if f[j] == 0 {
                f[j] = 1;
            }
            t2 += f[j];
        }
        t2 += 1;
        if t2 < TOTFREQ {
            f[big] += TOTFREQ - t2;
        } else if t2 - TOTFREQ >= f[big] / 2 {
            p = 0.98;
            continue;
        } else {
            f[big] -= t2 - TOTFREQ;
        }
        return;
    }
}

/// Write an order-0 table with the alphabet RLE; returns the symbols.
fn write_table0(t: &mut Vec<u8>, f: &[u32; 256], syms: &mut [EncSym; 256]) {
    let mut rle = 0usize;
    let mut x = 0u32;
    for j in 0..256 {
        if f[j] == 0 {
            continue;
        }
        if rle > 0 {
            rle -= 1;
        } else {
            t.push(j as u8);
            if j > 0 && f[j - 1] != 0 {
                let mut k = j + 1;
                while k < 256 && f[k] != 0 {
                    k += 1;
                }
                rle = k - (j + 1);
                t.push(rle as u8);
            }
        }
        put_freq(t, f[j]);
        syms[j] = EncSym {
            start: x,
            freq: f[j],
        };
        x += f[j];
    }
    t.push(0);
}

fn header(order: u8, out_body_len: usize, in_len: usize) -> [u8; 9] {
    let mut h = [0u8; 9];
    h[0] = order;
    h[1..5].copy_from_slice(&(out_body_len as u32).to_le_bytes());
    h[5..9].copy_from_slice(&(in_len as u32).to_le_bytes());
    h
}

/// rANS 4x8 order 0.
pub fn encode_o0(input: &[u8]) -> Vec<u8> {
    let n = input.len();
    let mut table = Vec::with_capacity(1024);
    let mut syms = [EncSym::default(); 256];
    if n > 0 {
        let mut f = [0u32; 256];
        for &b in input {
            f[b as usize] += 1;
        }
        normalise_o0(&mut f, n);
        write_table0(&mut table, &f, &mut syms);
    } else {
        // Nothing to code; a valid table still needs a full alphabet.
        table.extend_from_slice(&[0x00, 0x8f, 0xff, 0x00]);
    }
    let mut rev = Vec::with_capacity(n + n / 4 + 16);
    let mut x = [RANS_L; 4];
    // Symbol i is decoded by state i % 4, in increasing i.
    for i in (0..n).rev() {
        put(&mut x[i & 3], &mut rev, syms[input[i] as usize]);
    }
    let body = finish(x, rev);
    let mut out = Vec::with_capacity(9 + table.len() + body.len());
    out.extend_from_slice(&header(0, table.len() + body.len(), n));
    out.extend_from_slice(&table);
    out.extend_from_slice(&body);
    out
}

/// rANS 4x8 order 1. Inputs shorter than 4 bytes are written as order 0, as
/// the spec requires.
pub fn encode_o1(input: &[u8]) -> Vec<u8> {
    let n = input.len();
    if n < 4 {
        return encode_o0(input);
    }
    let q = n >> 2;
    // Counts per context; segment starts (0, q, 2q, 3q) have context 0.
    let mut f = vec![[0u32; 256]; 256];
    let mut t = [0u32; 256];
    for seg in 0..4 {
        let start = seg * q;
        let end = if seg == 3 { n } else { start + q };
        let mut ctx = 0usize;
        for &b in &input[start..end] {
            f[ctx][b as usize] += 1;
            t[ctx] += 1;
            ctx = b as usize;
        }
    }
    let mut table = Vec::with_capacity(4096);
    let mut syms = vec![[EncSym::default(); 256]; 256];
    let mut rle_i = 0usize;
    for i in 0..256 {
        if t[i] == 0 {
            continue;
        }
        normalise_o1(&mut f[i], t[i]);
        if rle_i > 0 {
            rle_i -= 1;
        } else {
            table.push(i as u8);
            if i > 0 && t[i - 1] != 0 {
                let mut k = i + 1;
                while k < 256 && t[k] != 0 {
                    k += 1;
                }
                rle_i = k - (i + 1);
                table.push(rle_i as u8);
            }
        }
        write_table0(&mut table, &f[i], &mut syms[i]);
    }
    table.push(0);

    let mut rev = Vec::with_capacity(n + n / 4 + 16);
    let mut x = [RANS_L; 4];
    let ctx_of = |p: usize| -> usize {
        if p == 0 || p == q || p == 2 * q || p == 3 * q {
            0
        } else {
            input[p - 1] as usize
        }
    };
    // The decoder reads the q interleaved rounds (states 0..3), then the
    // remainder on state 3; encode in exactly the reverse order.
    for p in (4 * q..n).rev() {
        put(&mut x[3], &mut rev, syms[ctx_of(p)][input[p] as usize]);
    }
    for i in (0..q).rev() {
        for j in (0..4).rev() {
            let p = i + j * q;
            put(&mut x[j], &mut rev, syms[ctx_of(p)][input[p] as usize]);
        }
    }
    let body = finish(x, rev);
    let mut out = Vec::with_capacity(9 + table.len() + body.len());
    out.extend_from_slice(&header(1, table.len() + body.len(), n));
    out.extend_from_slice(&table);
    out.extend_from_slice(&body);
    out
}

struct Reader<'a> {
    b: &'a [u8],
    at: usize,
}

impl Reader<'_> {
    fn u8(&mut self) -> io::Result<u8> {
        let v = *self.b.get(self.at).ok_or_else(|| invalid("truncated"))?;
        self.at += 1;
        Ok(v)
    }

    fn u32(&mut self) -> io::Result<u32> {
        let s = self
            .b
            .get(self.at..self.at + 4)
            .ok_or_else(|| invalid("truncated"))?;
        self.at += 4;
        Ok(u32::from_le_bytes(s.try_into().unwrap()))
    }

    fn freq(&mut self) -> io::Result<u32> {
        let a = self.u8()? as u32;
        Ok(if a < 128 {
            a
        } else {
            ((a & 0x7f) << 8) | self.u8()? as u32
        })
    }
}

/// Decoding table for one context: frequency, cumulative start, slot→symbol.
struct Dec {
    freq: [u32; 256],
    start: [u32; 256],
    sym: Vec<u8>,
}

fn read_table0(r: &mut Reader) -> io::Result<Dec> {
    let mut freq = [0u32; 256];
    let mut s = r.u8()? as usize;
    let mut last = s;
    let mut rle = 0u32;
    loop {
        freq[s] = r.freq()?;
        if rle > 0 {
            rle -= 1;
            s += 1;
            if s > 255 {
                return Err(invalid("alphabet RLE past 255"));
            }
        } else {
            s = r.u8()? as usize;
            if s == last + 1 {
                rle = r.u8()? as u32;
            }
        }
        last = s;
        if s == 0 {
            break;
        }
    }
    let mut start = [0u32; 256];
    let mut sym = vec![0u8; TOTFREQ as usize];
    let mut x = 0u32;
    for j in 0..256 {
        start[j] = x;
        if x + freq[j] > TOTFREQ {
            return Err(invalid("frequencies sum past 4096"));
        }
        for slot in x..x + freq[j] {
            sym[slot as usize] = j as u8;
        }
        x += freq[j];
    }
    if x < TOTFREQ - 1 {
        return Err(invalid("frequencies sum below 4095"));
    }
    Ok(Dec { freq, start, sym })
}

#[inline]
fn step(x: &mut u32, d: &Dec, r: &mut Reader) -> io::Result<u8> {
    let slot = *x & (TOTFREQ - 1);
    let s = d.sym[slot as usize];
    let f = d.freq[s as usize];
    if f == 0 {
        return Err(invalid("slot maps to a zero-frequency symbol"));
    }
    *x = f * (*x >> TF_SHIFT) + slot - d.start[s as usize];
    while *x < RANS_L {
        *x = (*x << 8) | r.u8()? as u32;
    }
    Ok(s)
}

/// Decode a whole rANS 4x8 stream (either order).
pub fn decode(data: &[u8]) -> io::Result<Vec<u8>> {
    let mut r = Reader { b: data, at: 0 };
    let order = r.u8()?;
    let body_len = r.u32()? as usize;
    let n = r.u32()? as usize;
    if 9 + body_len != data.len() {
        return Err(invalid("compressed size does not match the block"));
    }
    let mut out = vec![0u8; n];
    match order {
        0 => {
            let d = read_table0(&mut r)?;
            let mut x = [r.u32()?, r.u32()?, r.u32()?, r.u32()?];
            for (i, o) in out.iter_mut().enumerate() {
                *o = step(&mut x[i & 3], &d, &mut r)?;
            }
        }
        1 => {
            let mut tables: Vec<Option<Dec>> = (0..256).map(|_| None).collect();
            let mut c = r.u8()? as usize;
            let mut last = c;
            let mut rle = 0u32;
            loop {
                tables[c] = Some(read_table0(&mut r)?);
                if rle > 0 {
                    rle -= 1;
                    c += 1;
                    if c > 255 {
                        return Err(invalid("context RLE past 255"));
                    }
                } else {
                    c = r.u8()? as usize;
                    if c == last + 1 {
                        rle = r.u8()? as u32;
                    }
                }
                last = c;
                if c == 0 {
                    break;
                }
            }
            let mut x = [r.u32()?, r.u32()?, r.u32()?, r.u32()?];
            let mut ctx = [0usize; 4];
            let q = n >> 2;
            for i in 0..q {
                for j in 0..4 {
                    let d = tables[ctx[j]]
                        .as_ref()
                        .ok_or_else(|| invalid("missing context table"))?;
                    let s = step(&mut x[j], d, &mut r)?;
                    out[i + j * q] = s;
                    ctx[j] = s as usize;
                }
            }
            for o in out.iter_mut().skip(4 * q) {
                let d = tables[ctx[3]]
                    .as_ref()
                    .ok_or_else(|| invalid("missing context table"))?;
                let s = step(&mut x[3], d, &mut r)?;
                *o = s;
                ctx[3] = s as usize;
            }
        }
        _ => return Err(invalid("order is neither 0 nor 1")),
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *seed >> 33
    }

    fn inputs() -> Vec<Vec<u8>> {
        let mut v: Vec<Vec<u8>> = vec![
            vec![],
            vec![7],
            vec![1, 2],
            vec![1, 2, 3],
            b"abra".to_vec(),
            b"abracadabra".to_vec(),
            b"abracadabraabracadabraabracadabraabracadabrad".to_vec(),
            vec![0; 1000],
            vec![255; 4097],
            (0..=255u8).collect(),
            (0..=255u8).cycle().take(100_003).collect(),
        ];
        let mut seed = 42u64;
        for &n in &[5usize, 6, 7, 64, 1023, 65_537, 300_001] {
            // Skewed, quality-like data.
            v.push(
                (0..n)
                    .map(|_| {
                        let r = lcg(&mut seed) % 100;
                        if r < 80 { 37 } else { (r % 41) as u8 + 2 }
                    })
                    .collect(),
            );
            // Uniform bytes.
            v.push((0..n).map(|_| lcg(&mut seed) as u8).collect());
        }
        // One very frequent symbol plus many singletons (the corner case
        // that forces renormalisation in htscodecs).
        let mut x = vec![b'A'; 10_000];
        x.extend(0..=255u8);
        v.push(x);
        v
    }

    #[test]
    fn order0_round_trips() {
        for input in inputs() {
            let enc = encode_o0(&input);
            assert_eq!(decode(&enc).unwrap(), input, "len {}", input.len());
        }
    }

    #[test]
    fn order1_round_trips() {
        for input in inputs() {
            let enc = encode_o1(&input);
            assert_eq!(decode(&enc).unwrap(), input, "len {}", input.len());
        }
    }

    #[test]
    fn spec_order0_table_example() {
        // CRAMcodecs.tex, "Order-0 encoding": the table for "abracadabra",
        // byte for byte (htscodecs' normalisation gives the spec's numbers).
        let enc = encode_o0(b"abracadabra");
        assert_eq!(enc[0], 0);
        assert_eq!(&enc[5..9], &11u32.to_le_bytes());
        let table = [
            0x61, 0x87, 0x47, 0x62, 0x02, 0x82, 0xe8, 0x81, 0x74, 0x81, 0x74, 0x72, 0x82, 0xe8,
            0x00,
        ];
        assert_eq!(&enc[9..9 + table.len()], &table);
    }

    /// Run `f` on a thread; fail if it does not finish in time.
    fn within<T: Send + 'static>(secs: u64, f: impl FnOnce() -> T + Send + 'static) -> T {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(f());
        });
        rx.recv_timeout(std::time::Duration::from_secs(secs))
            .expect("rANS encode did not finish (normalisation loop?)")
    }

    #[test]
    fn renormalisation_terminates() {
        // 200 singleton symbols plus 56 symbols ×1000: the first pass
        // overshoots 4095 by more than half the largest frequency.
        let mut v = Vec::new();
        for s in 0..200u8 {
            v.push(s);
        }
        for s in 200..=255u8 {
            v.extend(std::iter::repeat_n(s, 1000));
        }
        let input = v.clone();
        let (e0, e1) = within(10, move || (encode_o0(&input), encode_o1(&input)));
        assert_eq!(decode(&e0).unwrap(), v);
        assert_eq!(decode(&e1).unwrap(), v);
    }

    #[test]
    fn skewed_distributions_round_trip() {
        let mut seed = 7u64;
        for case in 0..300 {
            let symbols = 1 + (lcg(&mut seed) % 256) as usize;
            let mut v = Vec::new();
            for s in 0..symbols {
                // Mix of singletons and heavy symbols.
                let reps = match lcg(&mut seed) % 4 {
                    0 => 1,
                    1 => 1 + lcg(&mut seed) % 5,
                    2 => 1 + lcg(&mut seed) % 200,
                    _ => 1 + lcg(&mut seed) % 5000,
                };
                v.extend(std::iter::repeat_n(s as u8, reps as usize));
            }
            // Shuffle a little so order-1 contexts vary.
            for i in (1..v.len()).rev().step_by(3) {
                let j = (lcg(&mut seed) as usize) % (i + 1);
                v.swap(i, j);
            }
            let input = v.clone();
            let (e0, e1) = within(20, move || (encode_o0(&input), encode_o1(&input)));
            assert_eq!(decode(&e0).unwrap(), v, "case {case} o0");
            assert_eq!(decode(&e1).unwrap(), v, "case {case} o1");
        }
    }

    #[test]
    fn rejects_garbage() {
        assert!(decode(&[2, 0, 0, 0, 0, 0, 0, 0, 0]).is_err());
        assert!(decode(&[0, 1, 0]).is_err());
    }
}
