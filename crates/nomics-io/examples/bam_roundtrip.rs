//! Golden round trip for the BAM codec and BGZF block layout.
//!
//! `cargo run --release -p nomics-io --example bam_roundtrip -- <golden.bam> [--rewrite]`
//!
//! For every record: decode to `Record`, re-encode, compare with the golden
//! bytes. The header is checked the same way. The BGZF block-size sequence is
//! compared with what `Profile::HTSJDK` would produce for the same stream.
//! With `--rewrite`, the re-encoded stream is also compressed with our
//! htsjdk-exact BGZF writer (zlib-exact level 2) and every output byte is
//! compared with the golden file.
//! Exits non-zero on any difference.

use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::process::ExitCode;
use std::time::Instant;

use nomics_io::bam::{BinRule, Header, Record, read_raw_record, write_raw_record};
use nomics_io::bgzf::{BgzfWriter, ParallelBgzfWriter, Profile, decode_block, read_raw_block};
use nomics_io::deflate::HtsjdkDeflater;

enum Writer {
    Sequential(Box<BgzfWriter<CompareSink, HtsjdkDeflater>>),
    Parallel(Box<ParallelBgzfWriter<CompareSink>>),
}

impl Write for Writer {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        match self {
            Writer::Sequential(w) => w.write(data),
            Writer::Parallel(w) => w.write(data),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Writer {
    fn finish(self) -> io::Result<CompareSink> {
        match self {
            Writer::Sequential(w) => w.finish(),
            Writer::Parallel(w) => w.finish(),
        }
    }
}

/// A writer that compares everything written to it with a file's bytes.
struct CompareSink {
    golden: BufReader<File>,
    buf: Vec<u8>,
    offset: u64,
    first_diff: Option<u64>,
}

impl Write for CompareSink {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.buf.resize(data.len(), 0);
        let mut got = 0;
        while got < data.len() {
            match self.golden.read(&mut self.buf[got..])? {
                0 => break,
                n => got += n,
            }
        }
        if self.first_diff.is_none() {
            let pos = self.buf[..got].iter().zip(data).position(|(a, b)| a != b);
            if let Some(p) = pos {
                self.first_diff = Some(self.offset + p as u64);
            } else if got < data.len() {
                self.first_diff = Some(self.offset + got as u64);
            }
        }
        self.offset += data.len() as u64;
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Decompressed BGZF stream that logs each block's payload size and can
/// capture the bytes handed out (to compare the header).
struct Stream<R> {
    inner: R,
    raw: Vec<u8>,
    block: Vec<u8>,
    pos: usize,
    coffset: u64,
    isizes: Vec<u32>,
    capture: Option<Vec<u8>>,
}

impl<R: Read> Read for Stream<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.pos == self.block.len() {
            if !read_raw_block(&mut self.inner, &mut self.raw, self.coffset)? {
                return Ok(0);
            }
            self.block.clear();
            self.pos = 0;
            let info = decode_block(&self.raw, self.coffset, &mut self.block)?;
            self.coffset += info.csize as u64;
            self.isizes.push(info.isize);
        }
        let n = buf.len().min(self.block.len() - self.pos);
        buf[..n].copy_from_slice(&self.block[self.pos..self.pos + n]);
        if let Some(c) = &mut self.capture {
            c.extend_from_slice(&buf[..n]);
        }
        self.pos += n;
        Ok(n)
    }
}

fn describe_diff(golden: &Record, ours: &[u8], raw: &[u8]) -> String {
    let first = ours
        .iter()
        .zip(raw)
        .position(|(a, b)| a != b)
        .unwrap_or(ours.len().min(raw.len()));
    let field = match first {
        0..=3 => "refID",
        4..=7 => "pos",
        8 => "l_read_name",
        9 => "mapq",
        10..=11 => "bin",
        12..=13 => "n_cigar_op",
        14..=15 => "flag",
        16..=19 => "l_seq",
        20..=31 => "next/tlen",
        _ => "variable part (name/cigar/seq/qual/tags)",
    };
    format!(
        "{} first differing byte {first} ({field}); len ours {} golden {}",
        String::from_utf8_lossy(&golden.read_name),
        ours.len(),
        raw.len()
    )
}

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: bam_roundtrip <golden.bam> [--rewrite]");
        return ExitCode::from(2);
    };
    let rewrite = std::env::args().any(|a| a == "--rewrite");
    let parallel = std::env::args().any(|a| a == "--parallel");
    let mut writer = rewrite.then(|| {
        let golden = File::open(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let sink = CompareSink {
            golden: BufReader::with_capacity(1 << 20, golden),
            buf: Vec::new(),
            offset: 0,
            first_diff: None,
        };
        if parallel {
            let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
            Writer::Parallel(Box::new(ParallelBgzfWriter::new(sink, threads)))
        } else {
            Writer::Sequential(Box::new(BgzfWriter::new(
                sink,
                HtsjdkDeflater::new(),
                Profile::HTSJDK,
            )))
        }
    });
    let file = File::open(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut s = Stream {
        inner: BufReader::with_capacity(1 << 20, file),
        raw: Vec::new(),
        block: Vec::new(),
        pos: 0,
        coffset: 0,
        isizes: Vec::new(),
        capture: Some(Vec::new()),
    };
    let t = Instant::now();
    let mut ok = true;

    let header = Header::read(&mut s).expect("header");
    let golden_header = s.capture.take().unwrap();
    let mut ours = Vec::new();
    header.write(&mut ours).unwrap();
    if let Some(w) = &mut writer {
        w.write_all(&ours).unwrap();
    }
    let header_match = ours == golden_header;
    ok &= header_match;
    println!(
        "header: {} bytes, {} refs, re-encode {}",
        golden_header.len(),
        header.refs.len(),
        if header_match { "IDENTICAL" } else { "DIFFERS" }
    );

    let (mut n, mut mismatches, mut bin_mismatches, mut bytes) =
        (0u64, 0u64, 0u64, golden_header.len() as u64);
    let mut raw = Vec::new();
    let mut enc = Vec::new();
    while read_raw_record(&mut s, &mut raw).expect("record") {
        n += 1;
        bytes += 4 + raw.len() as u64;
        let rec = Record::decode(&raw).unwrap_or_else(|e| panic!("record {n}: {e}"));
        let golden_bin = u16::from_le_bytes([raw[10], raw[11]]);
        if rec.bin(BinRule::Htsjdk) != golden_bin {
            bin_mismatches += 1;
            if bin_mismatches <= 5 {
                eprintln!(
                    "bin mismatch: {} pos {} flag {} ours {} golden {golden_bin}",
                    String::from_utf8_lossy(&rec.read_name),
                    rec.pos,
                    rec.flag,
                    rec.bin(BinRule::Htsjdk)
                );
            }
        }
        enc.clear();
        rec.encode(BinRule::Htsjdk, &mut enc)
            .unwrap_or_else(|e| panic!("record {n}: {e}"));
        if let Some(w) = &mut writer {
            write_raw_record(w, &enc).unwrap();
        }
        if enc != raw {
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!("record {n}: {}", describe_diff(&rec, &enc, &raw));
            }
        }
    }
    ok &= mismatches == 0 && bin_mismatches == 0;
    println!("records: {n}, re-encode mismatches {mismatches}, bin mismatches {bin_mismatches}");

    // Expected layout under Profile::HTSJDK: full blocks, one partial, EOF.
    let payload = Profile::HTSJDK.block_payload as u64;
    let mut expected: Vec<u32> = vec![payload as u32; (bytes / payload) as usize];
    if bytes % payload != 0 {
        expected.push((bytes % payload) as u32);
    }
    expected.push(0);
    let layout_match = s.isizes == expected;
    ok &= layout_match;
    println!(
        "bgzf: {} blocks, {bytes} uncompressed bytes, htsjdk layout {}",
        s.isizes.len(),
        if layout_match { "IDENTICAL" } else { "DIFFERS" }
    );
    if !layout_match {
        let first = s.isizes.iter().zip(&expected).position(|(a, b)| a != b);
        eprintln!(
            "first layout difference at block {first:?}: golden {:?}",
            first.map(|i| &s.isizes[i..(i + 4).min(s.isizes.len())])
        );
    }
    if let Some(w) = writer {
        let mut sink = w.finish().unwrap();
        let mut rest = Vec::new();
        sink.golden.read_to_end(&mut rest).unwrap();
        let identical = sink.first_diff.is_none() && rest.is_empty();
        ok &= identical;
        match sink.first_diff {
            None if rest.is_empty() => println!(
                "rewrite: {} compressed bytes, IDENTICAL to golden",
                sink.offset
            ),
            None => println!(
                "rewrite: ours ends at byte {}; golden has {} more bytes",
                sink.offset,
                rest.len()
            ),
            Some(at) => println!(
                "rewrite: DIFFERS from golden at compressed byte {at} (ours {} bytes)",
                sink.offset
            ),
        }
    }
    println!("elapsed {:.1}s", t.elapsed().as_secs_f64());
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
