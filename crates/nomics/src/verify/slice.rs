//! The read slice both pipelines run on: the first N pairs of the user's
//! FASTQs, records copied byte for byte. It needs no alignment, so it is the
//! same on every machine and every nomics version; the 1/12-genome tier was
//! made the same way (REPRODUCE.md). BWA's batching and insert-size
//! estimation then see a prefix of a real run.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use nomics_gate::Sha256;
use nomics_io::bgzf::{BgzfWriter, Profile};
use nomics_io::deflate::HtsjdkDeflater;
use nomics_io::gzip::GzReader;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Slice {
    pub pairs_requested: u64,
    pub pairs_taken: u64,
    /// SHA-256 over the decompressed slice (section 2 of docs/VERIFY.md).
    pub sha256: String,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// One FASTQ record (four lines, each with its newline), or `None` at a clean end.
fn next_record<R: BufRead>(r: &mut R, buf: &mut Vec<u8>, what: &str, n: u64) -> io::Result<bool> {
    buf.clear();
    for line in 0..4 {
        let start = buf.len();
        let got = r.read_until(b'\n', buf)?;
        if got == 0 {
            if line == 0 {
                return Ok(false);
            }
            return Err(invalid(format!("{what}: file ends inside record {}", n + 1)));
        }
        if buf.last() != Some(&b'\n') {
            return Err(invalid(format!(
                "{what}: record {} has no final newline",
                n + 1
            )));
        }
        let first = buf[start];
        if (line == 0 && first != b'@') || (line == 2 && first != b'+') {
            return Err(invalid(format!("{what}: record {} is not FASTQ", n + 1)));
        }
    }
    Ok(true)
}

/// The read name as BWA pairs it: up to the first whitespace, without `/1` or `/2`.
fn pair_name(record: &[u8]) -> &[u8] {
    let line = &record[1..];
    let end = line
        .iter()
        .position(|&b| b == b' ' || b == b'\t' || b == b'\n' || b == b'\r')
        .unwrap_or(line.len());
    let name = &line[..end];
    match name {
        [rest @ .., b'/', b'1' | b'2'] => rest,
        _ => name,
    }
}

/// Write the first `pairs` pairs of `fastq1`/`fastq2` (gzip or BGZF) as
/// BGZF to `out1`/`out2`. Fewer pairs in the input is not an error: all of
/// them are taken and the count says so.
pub fn head(fastq1: &Path, fastq2: &Path, pairs: u64, out1: &Path, out2: &Path) -> io::Result<Slice> {
    let open = |p: &Path| -> io::Result<BufReader<GzReader>> {
        Ok(BufReader::with_capacity(1 << 20, GzReader::open(p)?))
    };
    let create = |p: &Path| -> io::Result<_> {
        Ok(BgzfWriter::new(
            BufWriter::new(File::create(p)?),
            HtsjdkDeflater::new(),
            Profile::HTSJDK,
        ))
    };
    let (mut r1, mut r2) = (open(fastq1)?, open(fastq2)?);
    let (mut w1, mut w2) = (create(out1)?, create(out2)?);
    let (mut h1, mut h2) = (Sha256::new(), Sha256::new());
    let (mut a, mut b) = (Vec::new(), Vec::new());
    let mut taken = 0u64;
    while taken < pairs {
        let more1 = next_record(&mut r1, &mut a, "--fastq1", taken)?;
        let more2 = next_record(&mut r2, &mut b, "--fastq2", taken)?;
        match (more1, more2) {
            (false, false) => break,
            (true, true) => {}
            _ => {
                return Err(invalid(format!(
                    "--fastq1 and --fastq2 have different numbers of reads (one ends after {taken} pairs)"
                )));
            }
        }
        if pair_name(&a) != pair_name(&b) {
            return Err(invalid(format!(
                "pair {} has different read names in --fastq1 and --fastq2",
                taken + 1
            )));
        }
        w1.write_all(&a)?;
        w2.write_all(&b)?;
        h1.update(&a);
        h2.update(&b);
        taken += 1;
    }
    w1.finish()?.flush()?;
    w2.finish()?.flush()?;
    if taken == 0 {
        return Err(invalid("the FASTQs hold no reads".to_owned()));
    }
    Ok(Slice {
        pairs_requested: pairs,
        pairs_taken: taken,
        sha256: combine(&h1.finalize(), &h2.finalize(), taken),
    })
}

/// `sha256("nomics-slice/1" || sha256(R1) || sha256(R2) || pairs as u64 LE)`.
fn combine(r1: &[u8; 32], r2: &[u8; 32], pairs: u64) -> String {
    let mut h = Sha256::new();
    h.update(b"nomics-slice/1");
    h.update(r1);
    h.update(r2);
    h.update(&pairs.to_le_bytes());
    super::hex(&h.finalize())
}

/// MD5 and size of a file's bytes as stored (what `md5sum` prints).
pub fn file_md5(path: &Path) -> io::Result<(String, u64)> {
    let mut f = File::open(path)?;
    let mut md5 = nomics_io::cram::md5::Md5::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut size = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        md5.update(&buf[..n]);
        size += n as u64;
    }
    Ok((nomics_io::cram::md5::hex(&md5.finish()), size))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gz(path: &Path, text: &str) {
        let mut w = BgzfWriter::new(File::create(path).unwrap(), HtsjdkDeflater::new(), Profile::HTSJDK);
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap().flush().unwrap();
    }

    fn reads(n: usize, mate: u8) -> String {
        (0..n)
            .map(|i| format!("@r{i}/{mate} extra\nACGT\n+\nIIII\n"))
            .collect()
    }

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("nomics-slice-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn takes_the_first_pairs_byte_for_byte() {
        let d = dir("head");
        gz(&d.join("1.fq.gz"), &reads(10, 1));
        gz(&d.join("2.fq.gz"), &reads(10, 2));
        let s = head(&d.join("1.fq.gz"), &d.join("2.fq.gz"), 4, &d.join("o1.gz"), &d.join("o2.gz")).unwrap();
        assert_eq!(s.pairs_taken, 4);
        let mut out = String::new();
        GzReader::open(d.join("o1.gz")).unwrap().read_to_string(&mut out).unwrap();
        assert_eq!(out, reads(4, 1));
        // Same input, same digest; more pairs, another digest.
        let again = head(&d.join("1.fq.gz"), &d.join("2.fq.gz"), 4, &d.join("p1.gz"), &d.join("p2.gz")).unwrap();
        assert_eq!(again, s);
        let all = head(&d.join("1.fq.gz"), &d.join("2.fq.gz"), 100, &d.join("q1.gz"), &d.join("q2.gz")).unwrap();
        assert_eq!(all.pairs_taken, 10);
        assert_ne!(all.sha256, s.sha256);
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn refuses_bad_pairs() {
        let d = dir("bad");
        gz(&d.join("1.fq.gz"), &reads(3, 1));
        gz(&d.join("short.fq.gz"), &reads(2, 2));
        gz(&d.join("other.fq.gz"), &reads(3, 2).replace("@r1/2", "@x1/2"));
        gz(&d.join("cut.fq.gz"), "@r0/2\nACGT\n+\n");
        gz(&d.join("empty.fq.gz"), "");
        let o = (d.join("o1.gz"), d.join("o2.gz"));
        let try_with = |r2: &str| head(&d.join("1.fq.gz"), &d.join(r2), 10, &o.0, &o.1);
        assert!(try_with("short.fq.gz").unwrap_err().to_string().contains("different numbers"));
        assert!(try_with("other.fq.gz").unwrap_err().to_string().contains("different read names"));
        assert!(try_with("cut.fq.gz").unwrap_err().to_string().contains("ends inside"));
        let none = head(&d.join("empty.fq.gz"), &d.join("empty.fq.gz"), 10, &o.0, &o.1);
        assert!(none.unwrap_err().to_string().contains("no reads"));
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn pair_names() {
        assert_eq!(pair_name(b"@SRR1.5/1\n"), b"SRR1.5");
        assert_eq!(pair_name(b"@SRR1.5 1:N:0\n"), b"SRR1.5");
        assert_eq!(pair_name(b"@A:1:B 2:N\n"), b"A:1:B");
    }

    #[test]
    fn md5_of_file() {
        let d = dir("md5");
        std::fs::write(d.join("a"), b"abc").unwrap();
        assert_eq!(file_md5(&d.join("a")).unwrap(), ("900150983cd24fb0d6963f7d28e17f72".to_owned(), 3));
        std::fs::remove_dir_all(&d).unwrap();
    }
}
