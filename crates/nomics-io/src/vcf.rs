//! Streaming VCF reader (plain or BGZF-compressed), yielding borrowed fields
//! of the eight fixed columns. Interpretation of the alleles (e.g. which
//! reference span a known site masks) is left to the consumer.

use std::fs;
use std::io;
use std::path::Path;

use crate::bgzf::decode_parallel;

/// The fixed VCF columns of one data line; `pos` is 1-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VcfRecord<'a> {
    pub chrom: &'a [u8],
    pub pos: u64,
    pub id: &'a [u8],
    pub ref_allele: &'a [u8],
    pub alt: &'a [u8],
    pub qual: &'a [u8],
    pub filter: &'a [u8],
    pub info: &'a [u8],
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

fn parse_line(line: &[u8]) -> io::Result<VcfRecord<'_>> {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    let mut f = line.splitn(9, |&c| c == b'\t');
    let mut next = || {
        f.next()
            .ok_or_else(|| invalid(format!("short VCF line: {}", String::from_utf8_lossy(line))))
    };
    let chrom = next()?;
    let pos_s = next()?;
    let pos = std::str::from_utf8(pos_s)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| invalid(format!("bad POS {:?}", String::from_utf8_lossy(pos_s))))?;
    Ok(VcfRecord {
        chrom,
        pos,
        id: next()?,
        ref_allele: next()?,
        alt: next()?,
        qual: next()?,
        filter: next()?,
        info: next()?,
    })
}

/// Call `f` for every data line of the VCF at `path`, in file order. BGZF
/// input is decoded on `threads` threads. Returns the header lines (`#...`).
pub fn for_each_record(
    path: &Path,
    threads: usize,
    f: &mut dyn FnMut(&VcfRecord) -> io::Result<()>,
) -> io::Result<Vec<String>> {
    let data = fs::read(path)?;
    let mut header = Vec::new();
    let mut carry: Vec<u8> = Vec::new();
    let mut handle =
        |text: &[u8], header: &mut Vec<String>, carry: &mut Vec<u8>| -> io::Result<()> {
            let mut rest = text;
            if !carry.is_empty() {
                match rest.iter().position(|&c| c == b'\n') {
                    Some(nl) => {
                        carry.extend_from_slice(&rest[..nl]);
                        let line = std::mem::take(carry);
                        process(&line, header, f)?;
                        rest = &rest[nl + 1..];
                    }
                    None => {
                        carry.extend_from_slice(rest);
                        return Ok(());
                    }
                }
            }
            let last_nl = rest.iter().rposition(|&c| c == b'\n');
            let (complete, tail) = match last_nl {
                Some(i) => (&rest[..i], &rest[i + 1..]),
                None => (&rest[..0], rest),
            };
            if last_nl.is_some() {
                for line in complete.split(|&c| c == b'\n') {
                    process(line, header, f)?;
                }
            }
            carry.extend_from_slice(tail);
            Ok(())
        };
    if data.starts_with(&[0x1f, 0x8b]) {
        decode_parallel(&data, threads, &mut |b| handle(b, &mut header, &mut carry))?;
    } else {
        handle(&data, &mut header, &mut carry)?;
    }
    if !carry.is_empty() {
        let line = std::mem::take(&mut carry);
        process(&line, &mut header, f)?;
    }
    Ok(header)
}

fn process(
    line: &[u8],
    header: &mut Vec<String>,
    f: &mut dyn FnMut(&VcfRecord) -> io::Result<()>,
) -> io::Result<()> {
    if line.is_empty() {
        Ok(())
    } else if line[0] == b'#' {
        header.push(String::from_utf8_lossy(line).into_owned());
        Ok(())
    } else {
        f(&parse_line(line)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf::ParallelBgzfWriter;
    use std::io::Write;

    fn sample(n: usize) -> String {
        let mut s =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for i in 0..n {
            s += &format!(
                "chr{}\t{}\trs{i}\tAC\tA,ACT\t.\tPASS\tRS={i};VC=DIV\n",
                1 + i % 3,
                1000 + i * 7
            );
        }
        s
    }

    fn collect(path: &Path) -> (Vec<String>, Vec<(String, u64, String)>) {
        let mut recs = Vec::new();
        let header = for_each_record(path, 4, &mut |r| {
            recs.push((
                String::from_utf8_lossy(r.chrom).into_owned(),
                r.pos,
                String::from_utf8_lossy(r.alt).into_owned(),
            ));
            assert_eq!(r.ref_allele, b"AC");
            assert_eq!(r.filter, b"PASS");
            Ok(())
        })
        .unwrap();
        (header, recs)
    }

    #[test]
    fn plain_and_bgzf_give_the_same_records() {
        let dir = std::env::temp_dir().join(format!("nomics-vcf-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // Enough lines to span many BGZF blocks and decode batches, so lines
        // straddle block and batch boundaries.
        let text = sample(200_000);
        let plain = dir.join("x.vcf");
        std::fs::write(&plain, &text).unwrap();
        let gz = dir.join("x.vcf.gz");
        let mut w = ParallelBgzfWriter::new(Vec::new(), 3);
        w.write_all(text.as_bytes()).unwrap();
        std::fs::write(&gz, w.finish().unwrap()).unwrap();

        let (h1, r1) = collect(&plain);
        let (h2, r2) = collect(&gz);
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(h1.len(), 2);
        assert_eq!(h1, h2);
        assert_eq!(r1.len(), 200_000);
        assert_eq!(r1, r2);
        assert_eq!(r1[5], ("chr3".to_owned(), 1035, "A,ACT".to_owned()));
    }

    #[test]
    fn malformed_pos_is_an_error() {
        assert!(parse_line(b"chr1\tx\t.\tA\tC\t.\t.\t.").is_err());
        assert!(parse_line(b"chr1\t5").is_err());
    }
}
