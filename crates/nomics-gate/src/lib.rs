//! Comparison of our outputs (or a second reference run) against goldens at
//! each stage boundary, applying only the normalizations listed in
//! `docs/GATES.md`. Every difference is reported; nothing is silently masked.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread::JoinHandle;

use nomics_io::bam::{Header, Record, read_raw_record};
use nomics_io::bgzf::BgzfReader;
use nomics_io::gzip::GzReader;

/// A normalized artifact digest and the number of records or lines it contains.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactDigest {
    pub sha256: String,
    pub count: u64,
}

/// Streaming SHA-256 (FIPS 180-4), with runtime-dispatched ARMv8 SHA2 compression.
#[derive(Clone)]
pub struct Sha256 {
    state: [u32; 8],
    block: [u8; 64],
    block_len: usize,
    bytes: u64,
}

impl Default for Sha256 {
    fn default() -> Self {
        Self::new()
    }
}

impl Sha256 {
    pub fn new() -> Self {
        Self {
            state: [
                0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
                0x5be0cd19,
            ],
            block: [0; 64],
            block_len: 0,
            bytes: 0,
        }
    }

    pub fn update(&mut self, input: &[u8]) {
        self.update_with(input, nomics_simd::sha256_compress);
    }

    // A static kernel parameter lets tests force the portable path while using
    // exactly the same buffering and padding as production.
    #[inline]
    fn update_with(&mut self, mut input: &[u8], mut compress: impl FnMut(&mut [u32; 8], &[u8])) {
        self.bytes = self
            .bytes
            .checked_add(input.len() as u64)
            .expect("SHA-256 input is too long");
        if self.block_len != 0 {
            let take = (64 - self.block_len).min(input.len());
            self.block[self.block_len..self.block_len + take].copy_from_slice(&input[..take]);
            self.block_len += take;
            input = &input[take..];
            if self.block_len < 64 {
                return;
            }
            let block = self.block;
            compress(&mut self.state, &block);
            self.block_len = 0;
        }
        let full_len = input.len() / 64 * 64;
        if full_len != 0 {
            compress(&mut self.state, &input[..full_len]);
            input = &input[full_len..];
        }
        self.block[..input.len()].copy_from_slice(input);
        self.block_len = input.len();
    }

    pub fn finalize(self) -> [u8; 32] {
        self.finalize_with(nomics_simd::sha256_compress)
    }

    #[inline]
    fn finalize_with(mut self, mut compress: impl FnMut(&mut [u32; 8], &[u8])) -> [u8; 32] {
        let bits = self.bytes.wrapping_mul(8);
        self.block[self.block_len] = 0x80;
        self.block_len += 1;
        if self.block_len > 56 {
            self.block[self.block_len..].fill(0);
            let block = self.block;
            compress(&mut self.state, &block);
            self.block_len = 0;
        }
        self.block[self.block_len..56].fill(0);
        self.block[56..].copy_from_slice(&bits.to_be_bytes());
        let block = self.block;
        compress(&mut self.state, &block);
        let mut out = [0; 32];
        for (chunk, word) in out.as_chunks_mut::<4>().0.iter_mut().zip(self.state) {
            chunk.copy_from_slice(&word.to_be_bytes());
        }
        out
    }
}

fn hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for &byte in bytes {
        out.push(DIGITS[(byte >> 4) as usize] as char);
        out.push(DIGITS[(byte & 0xf) as usize] as char);
    }
    out
}

fn finish_digest(hash: Sha256, count: u64) -> ArtifactDigest {
    ArtifactDigest {
        sha256: hex(&hash.finalize()),
        count,
    }
}

fn hash_part(hash: &mut Sha256, bytes: &[u8]) {
    hash.update(&(bytes.len() as u64).to_le_bytes());
    hash.update(bytes);
}

/// Hash the exact decompressed bytes of one fastp chunk.
pub fn digest_fastq_chunk(path: &Path) -> io::Result<ArtifactDigest> {
    let mut reader = GzReader::open(path)?;
    let mut hash = Sha256::new();
    let mut newlines = 0u64;
    let mut buf = [0u8; 1 << 20];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        newlines += buf[..n].iter().filter(|&&b| b == b'\n').count() as u64;
        hash.update(&buf[..n]);
    }
    Ok(finish_digest(hash, newlines / 4))
}

/// Hash a BAM under the same normalization as [`compare_bams`].
pub fn digest_bam(path: &Path) -> io::Result<ArtifactDigest> {
    let mut reader = BgzfReader::new(BufReader::with_capacity(1 << 20, File::open(path)?));
    let header = Header::read(&mut reader)?;
    let lines = normalized_header_lines(&header);
    let mut hash = Sha256::new();
    hash.update(b"nomics-bam-digest-v1\0");
    hash.update(&(lines.len() as u64).to_le_bytes());
    for line in lines {
        hash_part(&mut hash, line.as_bytes());
    }
    hash.update(&(header.refs.len() as u64).to_le_bytes());
    for reference in header.refs {
        hash_part(&mut hash, &reference.name);
        hash.update(&reference.len.to_le_bytes());
    }
    let mut count = 0u64;
    let mut raw = Vec::new();
    while read_raw_record(&mut reader, &mut raw)? {
        hash_part(&mut hash, &raw);
        count += 1;
    }
    hash.update(b"records\0");
    hash.update(&count.to_le_bytes());
    Ok(finish_digest(hash, count))
}

fn digest_reader(mut reader: impl Read) -> io::Result<ArtifactDigest> {
    let mut hash = Sha256::new();
    let mut lines = 0u64;
    let mut any = false;
    let mut ended_in_newline = false;
    let mut buf = [0u8; 1 << 20];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        any = true;
        ended_in_newline = buf[n - 1] == b'\n';
        lines += buf[..n].iter().filter(|&&b| b == b'\n').count() as u64;
        hash.update(&buf[..n]);
    }
    if any && !ended_in_newline {
        lines += 1;
    }
    Ok(finish_digest(hash, lines))
}

/// Hash raw file bytes, reporting its physical line count.
pub fn digest_raw(path: &Path) -> io::Result<ArtifactDigest> {
    digest_reader(BufReader::new(File::open(path)?))
}

/// Hash mark-duplicate metrics after dropping physical lines 2 and 4, exactly
/// like the `sed '2d;4d'` gate in `scripts/full_exome.sh`.
pub fn digest_md_metrics(path: &Path) -> io::Result<ArtifactDigest> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut hash = Sha256::new();
    let mut physical = 0u64;
    let mut kept = 0u64;
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        physical += 1;
        if physical != 2 && physical != 4 {
            hash.update(&line);
            kept += 1;
        }
    }
    Ok(finish_digest(hash, kept))
}

/// Hash VCF header lines except `##GATKCommandLine`, followed by records,
/// using the same UTF-8 and line normalization as [`compare_vcfs`].
pub fn digest_vcf(path: &Path) -> io::Result<ArtifactDigest> {
    let mut text = String::new();
    GzReader::open(path)?.read_to_string(&mut text)?;
    let headers: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with('#') && !line.starts_with("##GATKCommandLine"))
        .collect();
    let records: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    let mut hash = Sha256::new();
    hash.update(b"nomics-vcf-digest-v1\0");
    hash.update(&(headers.len() as u64).to_le_bytes());
    for line in headers {
        hash_part(&mut hash, line.as_bytes());
    }
    hash.update(&(records.len() as u64).to_le_bytes());
    for line in &records {
        hash_part(&mut hash, line.as_bytes());
    }
    Ok(finish_digest(hash, records.len() as u64))
}

#[derive(Debug, Default)]
pub struct VcfReport {
    pub header_identical: bool,
    pub records_ours: usize,
    pub records_golden: usize,
    pub mismatches: usize,
    pub examples: Vec<(usize, String, String)>,
}

impl VcfReport {
    pub fn identical(&self) -> bool {
        self.header_identical && self.records_ours == self.records_golden && self.mismatches == 0
    }
}

impl std::fmt::Display for VcfReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "header (excl. GATKCommandLine): {}",
            if self.header_identical {
                "IDENTICAL"
            } else {
                "DIFFERS"
            }
        )?;
        writeln!(
            f,
            "records: {} ours / {} golden, {} mismatches",
            self.records_ours, self.records_golden, self.mismatches
        )?;
        for (i, ours, golden) in &self.examples {
            writeln!(f, "  record {i}\n    ours:   {ours}\n    golden: {golden}")?;
        }
        write!(
            f,
            "verdict: {}",
            if self.identical() {
                "IDENTICAL"
            } else {
                "DIFFERENT"
            }
        )
    }
}

/// Compare VCF headers without `##GATKCommandLine`, then records in order.
pub fn compare_vcfs(ours: &Path, golden: &Path, max_examples: usize) -> io::Result<VcfReport> {
    use std::io::Read;

    fn header(text: &str) -> Vec<&str> {
        text.lines()
            .filter(|line| line.starts_with('#') && !line.starts_with("##GATKCommandLine"))
            .collect()
    }

    fn records(text: &str) -> Vec<&str> {
        text.lines().filter(|line| !line.starts_with('#')).collect()
    }

    let read = |path: &Path| -> io::Result<String> {
        let mut text = String::new();
        GzReader::open(path)?.read_to_string(&mut text)?;
        Ok(text)
    };
    let (ours_text, golden_text) = (read(ours)?, read(golden)?);
    let (a, b) = (records(&ours_text), records(&golden_text));
    let mut report = VcfReport {
        header_identical: header(&ours_text) == header(&golden_text),
        records_ours: a.len(),
        records_golden: b.len(),
        ..Default::default()
    };
    report.mismatches = a.len().abs_diff(b.len());
    for (i, (ours, golden)) in a.iter().zip(&b).enumerate() {
        if ours != golden {
            report.mismatches += 1;
            if report.examples.len() < max_examples {
                report
                    .examples
                    .push((i, (*ours).to_owned(), (*golden).to_owned()));
            }
        }
    }
    Ok(report)
}

/// Header lines removed before comparing (GATES.md: `@PG` order and content
/// vary with task completion order in the reference itself).
fn normalized_header_lines(h: &Header) -> Vec<String> {
    String::from_utf8_lossy(&h.text)
        .trim_end_matches('\0')
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(str::to_owned)
        .collect()
}

/// Records decoded on a background thread, in batches.
struct RecordStream {
    rx: Receiver<Vec<Vec<u8>>>,
    batch: std::vec::IntoIter<Vec<u8>>,
    worker: Option<JoinHandle<io::Result<()>>>,
}

impl RecordStream {
    fn open(path: &Path) -> io::Result<(Header, RecordStream)> {
        let mut r = BgzfReader::new(BufReader::with_capacity(1 << 20, File::open(path)?));
        let header = Header::read(&mut r)?;
        let (tx, rx) = sync_channel(16);
        let worker = std::thread::spawn(move || {
            loop {
                let mut batch = Vec::with_capacity(4096);
                let mut raw = Vec::new();
                while batch.len() < 4096 && read_raw_record(&mut r, &mut raw)? {
                    batch.push(std::mem::take(&mut raw));
                }
                let done = batch.len() < 4096;
                if tx.send(batch).is_err() || done {
                    return Ok(());
                }
            }
        });
        Ok((
            header,
            RecordStream {
                rx,
                batch: Vec::new().into_iter(),
                worker: Some(worker),
            },
        ))
    }

    fn next(&mut self) -> io::Result<Option<Vec<u8>>> {
        loop {
            if let Some(r) = self.batch.next() {
                return Ok(Some(r));
            }
            match self.rx.recv() {
                Ok(b) if !b.is_empty() => self.batch = b.into_iter(),
                _ => {
                    if let Some(w) = self.worker.take() {
                        w.join()
                            .map_err(|_| io::Error::other("reader thread panicked"))??;
                    }
                    return Ok(None);
                }
            }
        }
    }
}

fn cigar_string(cigar: &[u32]) -> String {
    use nomics_io::bam::CIGAR_OPS;
    cigar
        .iter()
        .map(|&c| {
            format!(
                "{}{}",
                c >> 4,
                CIGAR_OPS
                    .get((c & 0xf) as usize)
                    .map_or('?', |&o| o as char)
            )
        })
        .collect()
}

/// Fields that differ between two decoded records.
fn field_diff(a: &Record, b: &Record) -> Vec<String> {
    let mut d = Vec::new();
    macro_rules! cmp {
        ($f:ident) => {
            if a.$f != b.$f {
                d.push(format!(
                    "{}: ours {:?} golden {:?}",
                    stringify!($f),
                    a.$f,
                    b.$f
                ));
            }
        };
    }
    if a.read_name != b.read_name {
        d.push(format!(
            "read_name: ours {} golden {}",
            String::from_utf8_lossy(&a.read_name),
            String::from_utf8_lossy(&b.read_name)
        ));
    }
    cmp!(flag);
    cmp!(ref_id);
    cmp!(pos);
    cmp!(mapq);
    if a.cigar != b.cigar {
        d.push(format!(
            "cigar: ours {} golden {}",
            cigar_string(&a.cigar),
            cigar_string(&b.cigar)
        ));
    }
    cmp!(next_ref_id);
    cmp!(next_pos);
    cmp!(tlen);
    if a.seq != b.seq {
        d.push("seq differs".into());
    }
    if a.qual != b.qual {
        d.push("qual differs".into());
    }
    if a.tags != b.tags {
        let names = |r: &Record| {
            r.tags
                .iter()
                .map(|(t, _)| String::from_utf8_lossy(t).into_owned())
                .collect::<Vec<_>>()
        };
        if names(a) != names(b) {
            d.push(format!(
                "tag order/set: ours {:?} golden {:?}",
                names(a),
                names(b)
            ));
        } else {
            for ((t, va), (_, vb)) in a.tags.iter().zip(&b.tags) {
                if va != vb {
                    d.push(format!(
                        "tag {}: ours {va:?} golden {vb:?}",
                        String::from_utf8_lossy(t)
                    ));
                }
            }
        }
    }
    if d.is_empty() {
        d.push("raw bytes differ but decoded fields are equal (bin or encoding)".into());
    }
    d
}

#[derive(Debug, Default)]
pub struct BamReport {
    pub ours: PathBuf,
    pub golden: PathBuf,
    /// `(only in ours, only in golden)` normalized header lines.
    pub header_diff: (Vec<String>, Vec<String>),
    pub refs_equal: bool,
    pub records_ours: u64,
    pub records_golden: u64,
    pub mismatches: u64,
    /// Up to `max_examples` mismatches: record index and field differences.
    pub examples: Vec<(u64, Vec<String>)>,
}

impl BamReport {
    pub fn identical(&self) -> bool {
        self.header_diff.0.is_empty()
            && self.header_diff.1.is_empty()
            && self.refs_equal
            && self.records_ours == self.records_golden
            && self.mismatches == 0
    }
}

impl std::fmt::Display for BamReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "ours:   {}\ngolden: {}",
            self.ours.display(),
            self.golden.display()
        )?;
        writeln!(
            f,
            "header (excl. @PG): {} | refs: {}",
            if self.header_diff.0.is_empty() && self.header_diff.1.is_empty() {
                "identical"
            } else {
                "DIFFERS"
            },
            if self.refs_equal {
                "identical"
            } else {
                "DIFFER"
            }
        )?;
        for l in &self.header_diff.0 {
            writeln!(f, "  + {l}")?;
        }
        for l in &self.header_diff.1 {
            writeln!(f, "  - {l}")?;
        }
        writeln!(
            f,
            "records: ours {} golden {}, mismatches {}",
            self.records_ours, self.records_golden, self.mismatches
        )?;
        for (i, d) in &self.examples {
            writeln!(f, "  record {i}: {}", d.join("; "))?;
        }
        write!(
            f,
            "verdict: {}",
            if self.identical() {
                "IDENTICAL"
            } else {
                "DIFFERENT"
            }
        )
    }
}

/// Compare two coordinate-sorted BAMs record by record (row-paired).
pub fn compare_bams(ours: &Path, golden: &Path, max_examples: usize) -> io::Result<BamReport> {
    let (ha, mut a) = RecordStream::open(ours)?;
    let (hb, mut b) = RecordStream::open(golden)?;
    let (la, lb) = (normalized_header_lines(&ha), normalized_header_lines(&hb));
    let mut report = BamReport {
        ours: ours.to_owned(),
        golden: golden.to_owned(),
        header_diff: (
            la.iter().filter(|l| !lb.contains(l)).cloned().collect(),
            lb.iter().filter(|l| !la.contains(l)).cloned().collect(),
        ),
        refs_equal: ha.refs == hb.refs,
        ..Default::default()
    };
    if report.header_diff == (vec![], vec![]) && la != lb {
        report
            .header_diff
            .0
            .push("(same lines, different order)".into());
    }
    loop {
        match (a.next()?, b.next()?) {
            (None, None) => break,
            (Some(_), None) => report.records_ours += 1,
            (None, Some(_)) => report.records_golden += 1,
            (Some(ra), Some(rb)) => {
                report.records_ours += 1;
                report.records_golden += 1;
                if ra != rb {
                    report.mismatches += 1;
                    if report.examples.len() < max_examples {
                        let diff = match (Record::decode(&ra), Record::decode(&rb)) {
                            (Ok(x), Ok(y)) => field_diff(&x, &y),
                            (Err(e), _) | (_, Err(e)) => vec![format!("undecodable record: {e}")],
                        };
                        report.examples.push((report.records_golden - 1, diff));
                    }
                }
            }
        }
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    #[test]
    fn digest_cache_hits_unchanged_files_only() {
        let dir = std::env::temp_dir().join(format!("nomics-digest-cache-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("a.txt");
        let cache_file = dir.join("cache.tsv");
        std::fs::write(&file, "one\n").unwrap();
        let d = digest_raw(&file).unwrap();
        let mut cache = DigestCache::load(&cache_file).unwrap();
        assert!(cache.get("raw", &file).is_none());
        cache.insert("raw", &file, &d);
        cache.save(&cache_file).unwrap();
        let cache = DigestCache::load(&cache_file).unwrap();
        assert_eq!(cache.get("raw", &file), Some(d));
        // Another kind of digest of the same file is a different entry.
        assert!(cache.get("vcf", &file).is_none());
        // A rewrite (new mtime, here also a new size) misses.
        std::thread::sleep(std::time::Duration::from_millis(20));
        std::fs::write(&file, "two!\n").unwrap();
        assert!(cache.get("raw", &file).is_none());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    use super::*;
    use nomics_io::bam::{BinRule, IntType, Reference, Value, write_raw_record};
    use nomics_io::bgzf::{BgzfWriter, Profile, Stored};
    use std::io::Write;

    fn rec(name: &str, pos: i32, flag: u16) -> Record {
        Record {
            ref_id: 0,
            pos,
            mapq: 60,
            flag,
            next_ref_id: 0,
            next_pos: pos + 100,
            tlen: 250,
            read_name: name.as_bytes().to_vec(),
            cigar: vec![(4 << 4)],
            seq: b"ACGT".to_vec(),
            qual: vec![30; 4],
            tags: vec![(*b"NM", Value::Int(IntType::U8, 0))],
        }
    }

    fn write_bam(path: &Path, header_text: &str, recs: &[Record]) {
        let h = Header {
            text: header_text.as_bytes().to_vec(),
            refs: vec![Reference {
                name: b"chr1".to_vec(),
                len: 1000,
            }],
        };
        let mut w = BgzfWriter::new(File::create(path).unwrap(), Stored, Profile::HTSJDK);
        h.write(&mut w).unwrap();
        for r in recs {
            let mut raw = Vec::new();
            r.encode(BinRule::Htsjdk, &mut raw).unwrap();
            write_raw_record(&mut w, &raw).unwrap();
        }
        w.finish().unwrap().flush().unwrap();
    }

    fn write_vcf(path: &Path, command: &str, record: &str) {
        let mut w = BgzfWriter::new(File::create(path).unwrap(), Stored, Profile::HTSJDK);
        writeln!(w, "##fileformat=VCFv4.2").unwrap();
        writeln!(w, "##GATKCommandLine={command}").unwrap();
        writeln!(w, "#CHROM\tPOS").unwrap();
        writeln!(w, "{record}").unwrap();
        w.finish().unwrap().flush().unwrap();
    }

    #[test]
    fn sha256_standard_vectors() {
        let cases = [
            (
                "",
                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            ),
            (
                "abc",
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            ),
            (
                "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
                "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
            ),
            (
                concat!(
                    "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn",
                    "hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu"
                ),
                "cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1",
            ),
        ];
        for (input, expected) in cases {
            for chunk_len in [1, 3, 64, 65, input.len().max(1)] {
                let mut hash = Sha256::new();
                let mut portable = Sha256::new();
                for chunk in input.as_bytes().chunks(chunk_len) {
                    hash.update(chunk);
                    portable.update_with(chunk, nomics_simd::sha256_compress_portable);
                }
                assert_eq!(hex(&hash.finalize()), expected);
                assert_eq!(
                    hex(&portable.finalize_with(nomics_simd::sha256_compress_portable)),
                    expected
                );
            }
        }
        let mut million_a = Sha256::new();
        let mut portable = Sha256::new();
        for _ in 0..1000 {
            million_a.update(&[b'a'; 1000]);
            portable.update_with(&[b'a'; 1000], nomics_simd::sha256_compress_portable);
        }
        let expected = "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0";
        assert_eq!(hex(&million_a.finalize()), expected);
        assert_eq!(
            hex(&portable.finalize_with(nomics_simd::sha256_compress_portable)),
            expected
        );
    }

    #[test]
    fn sha256_randomized_split_updates_match_portable() {
        let mut seed = 0xcafe_0123_4567_89abu64;
        let mut random = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        let lengths = (0..=300).chain([
            511,
            512,
            513,
            1023,
            1024,
            1025,
            4095,
            4096,
            4097,
            65535,
            65536,
            65537,
            1048576,
            1048576 + 55,
            1048576 + 56,
            1048576 + 63,
        ]);
        for len in lengths {
            // Offset the slice to exercise unaligned vector loads as well.
            let data: Vec<u8> = (0..=len).map(|_| random() as u8).collect();
            let input = &data[1..];
            let mut portable = Sha256::new();
            portable.update_with(input, nomics_simd::sha256_compress_portable);
            let expected = portable.finalize_with(nomics_simd::sha256_compress_portable);
            let mut whole = Sha256::new();
            whole.update(input);
            assert_eq!(whole.finalize(), expected, "whole length {len}");

            // Every possible two-part split for short messages, including
            // empty calls and the padding boundaries at 55/56 and 63/64.
            if len <= 300 {
                for split in 0..=len {
                    let mut hash = Sha256::new();
                    hash.update(&input[..split]);
                    hash.update(&[]);
                    hash.update(&input[split..]);
                    assert_eq!(hash.finalize(), expected, "length {len}, split {split}");
                }
            }
            for _ in 0..3 {
                let mut hash = Sha256::new();
                let mut portable = Sha256::new();
                let mut offset = 0;
                while offset < len {
                    let end = (offset + 1 + random() % 4096).min(len);
                    hash.update(&[]);
                    hash.update(&input[offset..end]);
                    portable
                        .update_with(&input[offset..end], nomics_simd::sha256_compress_portable);
                    offset = end;
                }
                assert_eq!(hash.finalize(), expected, "random split length {len}");
                assert_eq!(
                    portable.finalize_with(nomics_simd::sha256_compress_portable),
                    expected,
                    "portable random split length {len}"
                );
            }
        }
    }

    // cargo test --release -p nomics-gate sha256_throughput -- --ignored --nocapture
    #[test]
    #[ignore = "hashes 1 GiB per path; run in release mode for throughput"]
    fn sha256_throughput() {
        use std::hint::black_box;
        use std::time::Instant;

        fn hash_repeated<const PORTABLE: bool>(input: &[u8], repeats: usize) -> [u8; 32] {
            let mut hash = Sha256::new();
            for _ in 0..repeats {
                let input = black_box(input);
                if PORTABLE {
                    hash.update_with(input, nomics_simd::sha256_compress_portable);
                } else {
                    hash.update(input);
                }
            }
            black_box(if PORTABLE {
                hash.finalize_with(nomics_simd::sha256_compress_portable)
            } else {
                hash.finalize()
            })
        }

        // Reuse a modest buffer so the measurement covers hashing, with no
        // allocation or file I/O in the timed region. GB/s uses decimal bytes.
        let input = std::array::from_fn::<_, 65536, _>(|i| (i % 251) as u8);
        let repeats = (1 << 30) / input.len();
        #[cfg(target_arch = "aarch64")]
        eprintln!(
            "SHA2 detected: {}",
            std::arch::is_aarch64_feature_detected!("sha2")
        );
        hash_repeated::<true>(&input, 256);
        hash_repeated::<false>(&input, 256);
        let start = Instant::now();
        let portable = hash_repeated::<true>(&input, repeats);
        let portable_secs = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let dispatched = hash_repeated::<false>(&input, repeats);
        let dispatched_secs = start.elapsed().as_secs_f64();
        assert_eq!(portable, dispatched);
        let gb = (input.len() * repeats) as f64 / 1e9;
        eprintln!(
            "portable: {:.3} GB/s ({portable_secs:.3} s)",
            gb / portable_secs
        );
        eprintln!(
            "dispatched: {:.3} GB/s ({dispatched_secs:.3} s)",
            gb / dispatched_secs
        );
        eprintln!(
            "speedup: {:.2}x; digest: {}",
            portable_secs / dispatched_secs,
            hex(&dispatched)
        );
    }

    #[test]
    fn pg_lines_are_ignored_and_field_diffs_reported() {
        let dir = std::env::temp_dir().join(format!("nomics-gate-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (a, b) = (dir.join("a.bam"), dir.join("b.bam"));
        let recs = vec![rec("r1", 10, 99), rec("r2", 20, 147)];
        write_bam(&a, "@HD\tVN:1.6\n@PG\tID:x\n", &recs);
        write_bam(&b, "@HD\tVN:1.6\n@PG\tID:y\n", &recs);
        assert!(compare_bams(&a, &b, 5).unwrap().identical());
        assert_eq!(digest_bam(&a).unwrap(), digest_bam(&b).unwrap());

        let mut dup = recs.clone();
        dup[1].flag |= 0x400;
        write_bam(&b, "@HD\tVN:1.6\n", &dup);
        let r = compare_bams(&a, &b, 5).unwrap();
        assert!(!r.identical());
        assert_ne!(digest_bam(&a).unwrap(), digest_bam(&b).unwrap());
        assert_eq!(r.mismatches, 1);
        assert!(
            r.examples[0].1[0].starts_with("flag: ours 147 golden 1171"),
            "{:?}",
            r.examples
        );

        write_bam(&b, "@HD\tVN:1.6\n", &recs[..1]);
        let r = compare_bams(&a, &b, 5).unwrap();
        assert_eq!((r.records_ours, r.records_golden), (2, 1));
        assert!(!r.identical());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn vcf_gate_ignores_only_gatk_command_line() {
        let dir = std::env::temp_dir().join(format!("nomics-vcf-gate-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (a, b) = (dir.join("a.vcf.gz"), dir.join("b.vcf.gz"));
        write_vcf(&a, "ours", "chr1\t10");
        write_vcf(&b, "golden", "chr1\t10");
        assert!(compare_vcfs(&a, &b, 2).unwrap().identical());
        assert_eq!(digest_vcf(&a).unwrap(), digest_vcf(&b).unwrap());

        write_vcf(&b, "golden", "chr1\t11");
        let report = compare_vcfs(&a, &b, 2).unwrap();
        assert!(!report.identical());
        assert_ne!(digest_vcf(&a).unwrap(), digest_vcf(&b).unwrap());
        assert_eq!(report.mismatches, 1);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

/// One table of a GATKReport v1.1 file: name, column names, rows of
/// whitespace-separated cells (values, not padding), and its raw lines.
#[derive(Debug, Clone, PartialEq)]
pub struct ReportTable {
    pub name: String,
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
    pub raw: Vec<String>,
}

/// Parse GATKReport v1.1 text into tables (`#:GATKTable:` blocks).
pub fn parse_gatk_report(text: &str) -> io::Result<Vec<ReportTable>> {
    let mut tables = Vec::new();
    let mut lines = text.lines().peekable();
    if !lines
        .next()
        .is_some_and(|l| l.starts_with("#:GATKReport.v1.1:"))
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a GATKReport v1.1 file",
        ));
    }
    while let Some(line) = lines.next() {
        if !line.starts_with("#:GATKTable:") {
            continue;
        }
        // First header: #:GATKTable:<ncols>:<nrows>:<formats>:; then #:GATKTable:<name>:<desc>
        let name_line = lines.next().unwrap_or("");
        let name = name_line
            .trim_start_matches("#:GATKTable:")
            .split(':')
            .next()
            .unwrap_or("")
            .to_owned();
        let mut raw = vec![line.to_owned(), name_line.to_owned()];
        let header_line = lines.next().unwrap_or("");
        raw.push(header_line.to_owned());
        let header = header_line.split_whitespace().map(str::to_owned).collect();
        let mut rows = Vec::new();
        while let Some(l) = lines.peek() {
            if l.is_empty() || l.starts_with("#:") {
                break;
            }
            raw.push((*l).to_owned());
            rows.push(l.split_whitespace().map(str::to_owned).collect());
            lines.next();
        }
        tables.push(ReportTable {
            name,
            header,
            rows,
            raw,
        });
    }
    Ok(tables)
}

#[derive(Debug, Default)]
pub struct ReportComparison {
    pub bytes_identical: bool,
    /// Human-readable differences, most specific first, capped by the caller.
    pub diffs: Vec<String>,
}

/// Compare two GATKReport files: bytes first, then table by table, row by
/// row, cell by cell, reporting padding-only differences separately.
pub fn compare_gatk_reports(ours: &str, golden: &str, max: usize) -> io::Result<ReportComparison> {
    let mut c = ReportComparison {
        bytes_identical: ours == golden,
        ..Default::default()
    };
    if c.bytes_identical {
        return Ok(c);
    }
    let (a, b) = (parse_gatk_report(ours)?, parse_gatk_report(golden)?);
    let names = |t: &[ReportTable]| t.iter().map(|x| x.name.clone()).collect::<Vec<_>>();
    if names(&a) != names(&b) {
        c.diffs.push(format!(
            "tables: ours {:?} golden {:?}",
            names(&a),
            names(&b)
        ));
    }
    for tb in &b {
        let Some(ta) = a.iter().find(|t| t.name == tb.name) else {
            continue;
        };
        if ta.raw[0] != tb.raw[0] {
            c.diffs.push(format!(
                "{}: table header ours {:?} golden {:?}",
                tb.name, ta.raw[0], tb.raw[0]
            ));
        }
        if ta.header != tb.header {
            c.diffs.push(format!(
                "{}: columns ours {:?} golden {:?}",
                tb.name, ta.header, tb.header
            ));
        }
        if ta.rows.len() != tb.rows.len() {
            c.diffs.push(format!(
                "{}: {} rows, golden {}",
                tb.name,
                ta.rows.len(),
                tb.rows.len()
            ));
        }
        for (i, (ra, rb)) in ta.rows.iter().zip(&tb.rows).enumerate() {
            if ra != rb {
                let cells: Vec<String> = tb
                    .header
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| ra.get(*k) != rb.get(*k))
                    .map(|(k, col)| format!("{col}: ours {:?} golden {:?}", ra.get(k), rb.get(k)))
                    .collect();
                c.diffs.push(format!(
                    "{} row {i} ({}): {}",
                    tb.name,
                    rb.join(" "),
                    cells.join("; ")
                ));
            }
        }
        if ta.rows == tb.rows && ta.header == tb.header && ta.raw != tb.raw {
            let first = ta.raw.iter().zip(&tb.raw).position(|(x, y)| x != y);
            c.diffs.push(format!(
                "{}: values equal but padding/format differs, first at line {first:?}: ours {:?} golden {:?}",
                tb.name,
                first.map(|i| &ta.raw[i]),
                first.map(|i| &tb.raw[i])
            ));
        }
    }
    if c.diffs.is_empty() {
        c.diffs
            .push("bytes differ outside tables (blank lines or file ending)".into());
    }
    c.diffs.truncate(max);
    Ok(c)
}

#[cfg(test)]
mod report_tests {
    use super::*;

    const GOLDEN: &str = "#:GATKReport.v1.1:2\n#:GATKTable:2:2:%s:%s:;\n#:GATKTable:Arguments:Args\nArgument  Value\nfoo       1    \nbar       x    \n\n#:GATKTable:3:2:%d:%d:%.2f:;\n#:GATKTable:RecalTable1:\nQ  N   E\n11  5  1.00\n37  9  2.50\n\n";

    #[test]
    fn identical_reports() {
        let c = compare_gatk_reports(GOLDEN, GOLDEN, 10).unwrap();
        assert!(c.bytes_identical && c.diffs.is_empty());
    }

    #[test]
    fn cell_and_padding_differences_are_located() {
        let value = GOLDEN.replace("37  9  2.50", "37  9  2.49");
        let c = compare_gatk_reports(&value, GOLDEN, 10).unwrap();
        assert!(c.diffs[0].starts_with("RecalTable1 row 1"), "{:?}", c.diffs);
        assert!(c.diffs[0].contains("E: ours Some(\"2.49\") golden Some(\"2.50\")"));

        let padding = GOLDEN.replace("11  5  1.00", "11 5 1.00");
        let c = compare_gatk_reports(&padding, GOLDEN, 10).unwrap();
        assert!(
            c.diffs[0].contains("padding/format differs"),
            "{:?}",
            c.diffs
        );
    }
}

/// Compare fastp-style chunk files (`NNNN.<prefix>_R{1,2}.fastp.fastq.gz`) by
/// decompressed content. Returns one line per file, and whether all matched.
pub fn compare_fastq_chunks(ours: &Path, golden: &Path) -> io::Result<(Vec<String>, bool)> {
    use std::io::Read;
    let mut names: Vec<String> = std::fs::read_dir(golden)?
        .filter_map(|e| e.ok()?.file_name().into_string().ok())
        .filter(|n| n.ends_with(".fastp.fastq.gz"))
        .collect();
    names.sort();
    let mut lines = Vec::new();
    let mut all = !names.is_empty();
    for name in names {
        let read = |p: PathBuf| -> io::Result<Vec<u8>> {
            let mut v = Vec::new();
            GzReader::open(&p)?.read_to_end(&mut v)?;
            Ok(v)
        };
        let g = read(golden.join(&name))?;
        let line = match read(ours.join(&name)) {
            Err(e) => {
                all = false;
                format!("{name}: missing in ours ({e})")
            }
            Ok(o) if o == g => format!(
                "{name}: IDENTICAL ({} records)",
                g.iter().filter(|&&c| c == b'\n').count() / 4
            ),
            Ok(o) => {
                all = false;
                let at = o
                    .iter()
                    .zip(&g)
                    .position(|(a, b)| a != b)
                    .unwrap_or(o.len().min(g.len()));
                let line_no = g[..at.min(g.len())].iter().filter(|&&c| c == b'\n').count() + 1;
                format!(
                    "{name}: DIFFERS at byte {at} (golden line {line_no}); sizes ours {} golden {}",
                    o.len(),
                    g.len()
                )
            }
        };
        lines.push(line);
    }
    Ok((lines, all))
}

/// Digests of files that have not changed since they were last digested.
/// An entry is keyed by the digest kind and the resolved file's device,
/// inode, size and modification time (nanoseconds), so a rewrite or a new
/// file misses; nothing is reused across different content. Stored as TSV.
#[derive(Default)]
pub struct DigestCache {
    entries: std::collections::HashMap<(String, u64, u64, u64, i64), ArtifactDigest>,
    dirty: bool,
}

/// The cache key of `path` for digest `kind`: (kind, dev, inode, size, mtime ns).
fn digest_cache_key(kind: &str, path: &Path) -> io::Result<(String, u64, u64, u64, i64)> {
    use std::os::unix::fs::MetadataExt;
    let m = std::fs::metadata(path)?;
    let mtime = m.mtime() * 1_000_000_000 + m.mtime_nsec();
    Ok((kind.to_owned(), m.dev(), m.ino(), m.size(), mtime))
}

impl DigestCache {
    /// Load `path`; a missing file is an empty cache.
    pub fn load(path: &Path) -> io::Result<DigestCache> {
        let text = match std::fs::read_to_string(path) {
            Ok(t) => t,
            Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(DigestCache::default()),
            Err(e) => return Err(e),
        };
        let mut entries = std::collections::HashMap::new();
        for line in text.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 7 {
                continue;
            }
            let (Ok(dev), Ok(ino), Ok(size), Ok(mtime), Ok(count)) = (
                f[1].parse(),
                f[2].parse(),
                f[3].parse(),
                f[4].parse(),
                f[6].parse(),
            ) else {
                continue;
            };
            entries.insert(
                (f[0].to_owned(), dev, ino, size, mtime),
                ArtifactDigest {
                    sha256: f[5].to_owned(),
                    count,
                },
            );
        }
        Ok(DigestCache {
            entries,
            dirty: false,
        })
    }

    pub fn get(&self, kind: &str, path: &Path) -> Option<ArtifactDigest> {
        let key = digest_cache_key(kind, path).ok()?;
        self.entries.get(&key).cloned()
    }

    pub fn insert(&mut self, kind: &str, path: &Path, digest: &ArtifactDigest) {
        if let Ok(key) = digest_cache_key(kind, path) {
            self.entries.insert(key, digest.clone());
            self.dirty = true;
        }
    }

    /// Write the cache to `path` (atomically) if anything was added.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        if !self.dirty {
            return Ok(());
        }
        let mut text = String::from("kind\tdev\tinode\tsize\tmtime_ns\tsha256\tcount\n");
        let mut rows: Vec<_> = self.entries.iter().collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        for ((kind, dev, ino, size, mtime), d) in rows {
            text.push_str(&format!(
                "{kind}\t{dev}\t{ino}\t{size}\t{mtime}\t{}\t{}\n",
                d.sha256, d.count
            ));
        }
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir)?;
        }
        let tmp = path.with_extension("tsv.tmp");
        std::fs::write(&tmp, text)?;
        std::fs::rename(tmp, path)
    }
}
