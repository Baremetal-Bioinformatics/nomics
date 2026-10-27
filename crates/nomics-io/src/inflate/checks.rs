//! Differential, boundary, and throughput tests. No external codec dependency.
use super::*;
use crate::deflate::{Finish, ZlibFast};
use std::io::Write;
use std::process::{Command, Stdio};

fn next(seed: &mut u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed
}

fn compress(data: &[u8], level: u32) -> Vec<u8> {
    let mut out = Vec::new();
    assert_eq!(
        ZlibFast::new(level).deflate_finish(data, &mut out, usize::MAX),
        Finish::Done
    );
    out
}

fn compare(raw: &[u8], expected: &[u8]) {
    let mut old = b"prefix".to_vec();
    let mut new = old.clone();
    let old_used = oracle::inflate(raw, &mut old, expected.len()).unwrap();
    let new_used = inflate(raw, &mut new, expected.len()).unwrap();
    assert_eq!(new_used, old_used);
    assert_eq!(new, old);
    assert_eq!(&new[6..], expected);
    let mut streamed = Vec::new();
    assert_eq!(
        inflate_to(raw, &mut |s| {
            streamed.extend_from_slice(s);
            Ok(())
        })
        .unwrap(),
        old_used
    );
    assert_eq!(streamed, expected);
}

#[test]
fn randomized_encoder_differential_all_levels() {
    let mut seed = 0x1d31_eabb_9876_5432;
    for case in 0..160 {
        let n = match case {
            0 => 0,
            1 => 1,
            2 => 32768,
            3 => 65536,
            _ => next(&mut seed) as usize % 25000,
        };
        let mut data = Vec::new();
        while data.len() < n {
            let r = next(&mut seed);
            if r & 7 == 0 && !data.is_empty() {
                let d = 1 + next(&mut seed) as usize % data.len().min(32768);
                for _ in 0..(r as usize % 258).min(n - data.len()) {
                    data.push(data[data.len() - d]);
                }
            } else {
                data.push((r >> 32) as u8 % [1, 5, 64, 255][case % 4]);
            }
        }
        for level in 1..=3 {
            compare(&compress(&data, level), &data);
        }
        // The fast BGZF deflater, decoded by both inflaters.
        let mut fast = Vec::new();
        crate::deflate::FastDeflater::new().compress(&data, &mut fast);
        compare(&fast, &data);
        let mut stored = Vec::new();
        crate::bgzf::Deflate::deflate(&mut crate::bgzf::Stored, &data, &mut stored);
        compare(&stored, &data);
    }
}

struct Writer {
    data: Vec<u8>,
    bit: usize,
}
impl Writer {
    fn new() -> Self {
        Self {
            data: Vec::new(),
            bit: 0,
        }
    }
    fn bits(&mut self, value: u32, n: u32) {
        for i in 0..n {
            if self.bit / 8 == self.data.len() {
                self.data.push(0);
            }
            self.data[self.bit / 8] |= (((value >> i) & 1) as u8) << (self.bit % 8);
            self.bit += 1;
        }
    }
    fn code(&mut self, symbol: usize, lengths: &[u8]) {
        let len = lengths[symbol];
        let mut code = 0u32;
        for l in 1..=len {
            code = (code + lengths.iter().filter(|&&n| n == l - 1 && n != 0).count() as u32) << 1;
        }
        code += lengths[..symbol].iter().filter(|&&n| n == len).count() as u32;
        self.bits(code.reverse_bits() >> (32 - len), len as u32);
    }
}
fn fixed_lengths() -> [u8; 288] {
    let mut lit = [8; 288];
    lit[144..256].fill(9);
    lit[256..280].fill(7);
    lit
}

#[test]
fn every_distance_and_match_length_across_blocks() {
    let mut w = Writer::new();
    // A stored history followed by a fixed block exercises distance 32768,
    // which our encoder deliberately does not emit.
    w.bits(0, 3);
    w.bits(0, 5);
    w.bits(32768, 16);
    w.bits(32767, 16);
    let mut expected: Vec<u8> = (0..32768).map(|i| (i * 73 + i / 7) as u8).collect();
    for &b in &expected {
        w.bits(b as u32, 8);
    }
    w.bits(3, 3);
    let lit = fixed_lengths();
    for (d, len) in (1..=32768)
        .map(|d| (d, 258))
        .chain((3..=258).map(|len| (1, len)))
    {
        let l = LEN_BASE.iter().rposition(|&b| b as usize <= len).unwrap();
        w.code(257 + l, &lit);
        w.bits((len - LEN_BASE[l] as usize) as u32, LEN_EXTRA[l] as u32);
        let ds = DIST_BASE.iter().rposition(|&b| b as usize <= d).unwrap();
        w.code(ds, &[5; 32]);
        w.bits((d - DIST_BASE[ds] as usize) as u32, DIST_EXTRA[ds] as u32);
        for _ in 0..len {
            expected.push(expected[expected.len() - d]);
        }
    }
    w.code(256, &lit);
    compare(&w.data, &expected);
    // Slow path for every overlap distance and every length, at an exact limit.
    for d in 1..=258 {
        for len in 3..=258 {
            let at = d + 7;
            let mut out: Vec<u8> = (0..at).map(|i| (i * 29) as u8).collect();
            let mut expected = out.clone();
            for _ in 0..len {
                expected.push(expected[expected.len() - d]);
            }
            out.resize(at + len, 0);
            copy_match(&mut out, at, d, len);
            assert_eq!(out, expected);
        }
    }
}

fn dynamic_header(w: &mut Writer, lengths: &[u8; 258], last: bool) {
    w.bits(if last { 5 } else { 4 }, 3);
    w.bits(0, 5);
    w.bits(0, 5);
    w.bits(15, 4);
    let mut clen = [0; 19];
    clen[..16].fill(4);
    for &s in &CLEN_ORDER {
        w.bits(clen[s] as u32, 3);
    }
    for &len in lengths {
        w.code(len as usize, &clen);
    }
}

#[test]
fn randomized_huffman_trees_and_bit_alignments() {
    let mut seed = 0x2197_654a_bcad;
    let mut w = Writer::new();
    let mut expected = Vec::new();
    for block in 0..40 {
        let count = 2 + next(&mut seed) as usize % 255;
        let mut depths = vec![0u8];
        while depths.len() < count {
            let i = next(&mut seed) as usize % depths.len();
            if depths[i] < 15 {
                depths[i] += 1;
                depths.push(depths[i]);
            }
        }
        let mut symbols: Vec<usize> = (0..256).collect();
        for i in (1..symbols.len()).rev() {
            symbols.swap(i, next(&mut seed) as usize % (i + 1));
        }
        symbols.truncate(count - 1);
        let mut lengths = [0; 258];
        lengths[256] = depths[0];
        // Deliberately empty distance alphabets in some literal-only blocks.
        lengths[257] = if block % 2 == 0 { 0 } else { 1 };
        for (&symbol, &depth) in symbols.iter().zip(&depths[1..]) {
            lengths[symbol] = depth;
        }
        dynamic_header(&mut w, &lengths, block == 39);
        for _ in 0..512 {
            let symbol = symbols[next(&mut seed) as usize % symbols.len()];
            w.code(symbol, &lengths[..257]);
            expected.push(symbol as u8);
        }
        w.code(256, &lengths[..257]);
    }
    compare(&w.data, &expected);
    // Table lookahead must not consume a following gzip trailer/member.
    w.data.extend_from_slice(&[0xff; 16]);
    compare(&w.data, &expected);
}

#[test]
fn full_48_bit_length_distance_tokens() {
    let mut w = Writer::new();
    w.bits(0, 8);
    w.bits(32768, 16);
    w.bits(32767, 16);
    let mut expected: Vec<u8> = (0..32768).map(|i| (i * 37 + i / 31) as u8).collect();
    for &b in &expected {
        w.bits(b as u32, 8);
    }
    w.bits(5, 3);
    w.bits(29, 5);
    w.bits(29, 5);
    w.bits(15, 4);
    let mut clen = [0; 19];
    clen[..16].fill(4);
    for &s in &CLEN_ORDER {
        w.bits(clen[s] as u32, 3);
    }
    let mut lit = [0; 286];
    let mut dist = [0; 30];
    for i in 0..14 {
        lit[i] = i as u8 + 1;
        dist[i] = i as u8 + 1;
    }
    lit[256] = 15;
    lit[284] = 15;
    dist[14] = 15;
    dist[29] = 15;
    for &len in lit.iter().chain(&dist) {
        w.code(len as usize, &clen);
    }
    for _ in 0..512 {
        w.code(284, &lit);
        w.bits(30, 5); // 257 bytes
        w.code(29, &dist);
        w.bits(8191, 13); // distance 32768
        for _ in 0..257 {
            expected.push(expected[expected.len() - 32768]);
        }
    }
    w.code(256, &lit);
    compare(&w.data, &expected);
}

#[test]
fn malformed_dynamic_headers() {
    for lengths in [[0; 258], [1; 258]] {
        let mut w = Writer::new();
        dynamic_header(&mut w, &lengths, true);
        assert!(inflate(&w.data, &mut Vec::new(), 1000).is_err());
    }
    for repeat in [16, 18] {
        let mut w = Writer::new();
        w.bits(5, 3);
        w.bits(0, 5);
        w.bits(0, 5);
        w.bits(0, 4);
        // Code-length symbols 16 and 18, both one bit. 16 without a
        // preceding length is invalid; two maximum 18 runs overrun HLIT+HDIST.
        for len in [1, 0, 1, 0] {
            w.bits(len, 3);
        }
        for _ in 0..2 {
            w.bits(if repeat == 16 { 0 } else { 1 }, 1);
            w.bits(127, if repeat == 16 { 2 } else { 7 });
        }
        assert_eq!(
            inflate(&w.data, &mut Vec::new(), 1000),
            Err(InflateError::BadSymbol)
        );
    }
    assert_eq!(
        inflate(&[1, 1, 0, 0, 0], &mut Vec::new(), 1000),
        Err(InflateError::BadStoredLength)
    );
}

#[test]
fn fifteen_bit_codes_and_unused_prefixes() {
    let mut w = Writer::new();
    w.bits(5, 3); // final dynamic
    w.bits(0, 5);
    w.bits(0, 5);
    w.bits(15, 4);
    let mut clen = [0; 19];
    clen[..16].fill(4);
    for &s in &CLEN_ORDER {
        w.bits(clen[s] as u32, 3);
    }
    let mut lengths = [0; 258];
    for (i, l) in lengths[..15].iter_mut().enumerate() {
        *l = i as u8 + 1;
    }
    lengths[256] = 15;
    lengths[257] = 1;
    for &l in &lengths {
        w.code(l as usize, &clen);
    }
    let expected: Vec<u8> = (0..4096).map(|i| (i % 15) as u8).collect();
    for &s in &expected {
        w.code(s as usize, &lengths[..257]);
    }
    w.code(256, &lengths[..257]);
    compare(&w.data, &expected);
    let h = Huffman::new(&[1], Alphabet::Symbols).unwrap();
    assert_eq!(
        h.decode(&mut Bits::new(&[1, 0])),
        Err(InflateError::BadCode)
    );
    assert!(Huffman::new(&[16], Alphabet::Symbols).is_err());
    assert!(Huffman::new(&[1, 1, 1], Alphabet::Symbols).is_err());
}

#[test]
fn malformed_and_truncated_streams_are_bounded_errors() {
    let data = b"abcabcabcabcabcabcabcabcabcabcabcabc".repeat(100);
    for level in 1..=3 {
        let raw = compress(&data, level);
        for end in 0..raw.len() {
            assert!(inflate(&raw[..end], &mut Vec::new(), data.len()).is_err());
        }
        for limit in [0, 1, 257, data.len() - 1] {
            assert_eq!(
                inflate(&raw, &mut Vec::new(), limit),
                Err(InflateError::OutputTooLarge)
            );
        }
    }
    let lit = fixed_lengths();
    for reserved in [286, 287] {
        let mut w = Writer::new();
        w.bits(3, 3);
        w.code(reserved, &lit);
        assert_eq!(
            inflate(&w.data, &mut Vec::new(), 1000),
            Err(InflateError::BadSymbol)
        );
    }
    for ds in [0, 30, 31] {
        let mut w = Writer::new();
        w.bits(3, 3);
        w.code(257, &lit);
        w.code(ds, &[5; 32]);
        assert!(inflate(&w.data, &mut Vec::new(), 1000).is_err());
    }
    // Deterministic malformed-input smoke fuzz: a bounded output limit prevents
    // adversarial expansion; any successful prefix must agree with the oracle.
    let mut seed = 874934;
    for _ in 0..20000 {
        let n = next(&mut seed) as usize % 128;
        let raw: Vec<u8> = (0..n).map(|_| next(&mut seed) as u8).collect();
        let mut out = Vec::new();
        if let Ok(used) = inflate(&raw, &mut out, 4096) {
            let mut old = Vec::new();
            assert_eq!(oracle::inflate(&raw, &mut old, 4096).unwrap(), used);
            assert_eq!(out, old);
        }
    }
}

fn gzip(data: &[u8], level: u32) -> Option<Vec<u8>> {
    let mut child = match Command::new("gzip")
        .args(["-n", "-c", &format!("-{level}")])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Ok(c) => c,
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            eprintln!("gzip unavailable; CLI fixtures skipped");
            return None;
        }
        Err(e) => panic!("gzip: {e}"),
    };
    let mut stdin = child.stdin.take().unwrap();
    let result = std::thread::scope(|scope| {
        let writer = scope.spawn(move || stdin.write_all(data).unwrap());
        let result = child.wait_with_output().unwrap();
        writer.join().unwrap();
        result
    });
    assert!(result.status.success());
    Some(result.stdout)
}

#[test]
fn bgzf_members_match_oracle() {
    use crate::bgzf::{BamCompression, ParallelBgzfWriter};
    let data = corpus(true, 1000);
    for mode in [BamCompression::Htsjdk, BamCompression::Fast] {
        let mut writer = ParallelBgzfWriter::with_compression(Vec::new(), 3, mode);
        writer.write_all(&data).unwrap();
        bgzf_members_decode(&writer.finish().unwrap(), &data);
    }
}

fn bgzf_members_decode(members: &[u8], data: &[u8]) {
    let mut at = 0;
    for range in crate::bgzf::split_blocks(members).unwrap() {
        let member = &members[range];
        let size = u32::from_le_bytes(member[member.len() - 4..].try_into().unwrap()) as usize;
        compare(&member[18..member.len() - 8], &data[at..at + size]);
        let mut decoded = Vec::new();
        crate::bgzf::decode_block(member, 0, &mut decoded).unwrap();
        assert_eq!(decoded, data[at..at + size]);
        at += size;
    }
    assert_eq!(at, data.len());
}

// Deterministic reads drawn from a shared genome; varying qualities and record
// fields prevent this from degenerating into a benchmark of a single long run.
fn corpus(binary: bool, records: usize) -> Vec<u8> {
    let mut seed = 0x6512_efaa_7891;
    let genome: Vec<u8> = (0..65536)
        .map(|_| b"ACGT"[next(&mut seed) as usize % 4])
        .collect();
    let mut out = Vec::new();
    for id in 0..records {
        let pos = next(&mut seed) as usize % (genome.len() - 150);
        if binary {
            out.extend_from_slice(&[0x20, 0x01, 0, 0]);
            for n in [
                id as u32 % 24,
                pos as u32,
                0x1249_3c10,
                0x0096_0001,
                150,
                u32::MAX,
                0,
                0,
            ] {
                out.extend_from_slice(&n.to_le_bytes());
            }
            out.extend_from_slice(format!("read{id:010}\0").as_bytes());
            out.extend_from_slice(&2400u32.to_le_bytes()); // 150M
            for pair in genome[pos..pos + 150].chunks(2) {
                let encode = |b| match b {
                    b'A' => 1,
                    b'C' => 2,
                    b'G' => 4,
                    _ => 8,
                };
                out.push(encode(pair[0]) << 4 | encode(pair[1]));
            }
            for _ in 0..150 {
                out.push(25 + (next(&mut seed) % 16) as u8);
            }
            out.extend_from_slice(b"NMi\0\0\0\0RGZsample\0");
        } else {
            out.extend_from_slice(
                format!(
                    "@instrument:1:flowcell:1:{}:{}:{} 1:N:0:ACGT\n",
                    id / 1000,
                    id % 1000,
                    pos
                )
                .as_bytes(),
            );
            out.extend_from_slice(&genome[pos..pos + 150]);
            out.extend_from_slice(b"\n+\n");
            for _ in 0..150 {
                out.push(58 + (next(&mut seed) % 16) as u8);
            }
            out.push(b'\n');
        }
    }
    out
}

#[test]
fn gzip_cli_all_levels_and_concatenated_members() {
    let data = corpus(false, 200);
    let mut members = Vec::new();
    for level in 1..=9 {
        let Some(member) = gzip(&data, level) else {
            return;
        };
        compare(&member[10..member.len() - 8], &data);
        members.extend_from_slice(&member);
    }
    let mut got = Vec::new();
    crate::gzip::gunzip_to(&members, &mut |s| {
        got.extend_from_slice(s);
        Ok(())
    })
    .unwrap();
    assert_eq!(got, data.repeat(9));
}

#[test]
#[ignore = "throughput benchmark: cargo test -p nomics-io --release inflate_benchmark -- --ignored --nocapture"]
fn inflate_benchmark() {
    use std::hint::black_box;
    use std::time::{Duration, Instant};
    for (name, binary) in [("FASTQ-like", false), ("BAM-like", true)] {
        let data = corpus(binary, 12000);
        let mut inputs = vec![("nomics-2", compress(&data, 2))];
        if let Some(gz) = gzip(&data, 6) {
            inputs.push(("gzip-6", gz[10..gz.len() - 8].to_vec()));
        }
        for (encoder, raw) in inputs {
            compare(&raw, &data);
            type Decoder = fn(&[u8], &mut Vec<u8>, usize) -> Result<usize, String>;
            let old: Decoder = |s, o, n| oracle::inflate(s, o, n).map_err(|e| e.to_string());
            let new: Decoder = |s, o, n| inflate(s, o, n).map_err(|e| e.to_string());
            let mut rates = [0.; 2];
            // Alternate order across five rounds and report median MB/s, decimal.
            let mut samples = [Vec::new(), Vec::new()];
            for round in 0..5 {
                for index in if round % 2 == 0 { [0, 1] } else { [1, 0] } {
                    let decoder = [old, new][index];
                    let mut out = Vec::with_capacity(data.len());
                    let start = Instant::now();
                    let mut iterations = 0;
                    while start.elapsed() < Duration::from_millis(200) {
                        out.clear();
                        assert_eq!(
                            decoder(black_box(&raw), &mut out, data.len()).unwrap(),
                            raw.len()
                        );
                        black_box(&out);
                        iterations += 1;
                    }
                    samples[index].push(
                        data.len() as f64 * iterations as f64 / start.elapsed().as_secs_f64() / 1e6,
                    );
                    assert_eq!(out, data);
                }
            }
            for i in 0..2 {
                samples[i].sort_by(f64::total_cmp);
                rates[i] = samples[i][2];
            }
            eprintln!(
                "{name} {encoder}: {} bytes -> {} bytes; old {:.1} MB/s, new {:.1} MB/s, {:.2}x",
                raw.len(),
                data.len(),
                rates[0],
                rates[1],
                rates[1] / rates[0]
            );
        }
    }
}
