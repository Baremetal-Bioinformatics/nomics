//! The two merge steps sarek 3.10.0 runs after a scattered Mutect2: GATK 4.6.2.0 `MergeVcfs`
//! (Picard) on the per-group VCFs and `MergeMutectStats` on their stats files.
//!
//! Both are exact for what sarek gives them, the outputs of one Mutect2 configuration run on
//! disjoint interval groups, and refuse the rest instead of guessing:
//!
//! * `MergeVcfs` merges headers with htsjdk's `VCFUtils.smartMergeHeaders` and re-encodes
//!   every record. For Mutect2's scatter outputs the headers differ only in
//!   `##GATKCommandLine`, where the first input's line is kept, and re-encoding is the
//!   identity (the S1-r1 golden: 6,351 records). So the inputs' headers must be equal but
//!   for that line, and records are copied as text.
//! * Records come out in the header's contig order, then by position. Two inputs with a
//!   record at the same position are refused: Picard's order between them is not confirmed,
//!   and disjoint interval groups cannot produce it.
//! * `MergeMutectStats` sums each statistic over the inputs. Mutect2 writes one, `callable`;
//!   any other name is refused (GATK holds them in a hash map).

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use nomics_io::bgzf::ParallelBgzfWriter;
use nomics_io::gzip::GzReader;

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// A VCF as lines, gzip (BGZF) or plain.
fn open_lines(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let n = File::open(path)?.read(&mut magic)?;
    Ok(if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(GzReader::open(path)?))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    })
}

struct Input {
    name: String,
    lines: Box<dyn BufRead>,
    /// The next record, with its line ending, and its (contig index, position).
    next: Option<(Vec<u8>, (usize, u64))>,
    last: Option<(usize, u64)>,
}

impl Input {
    /// Read the header (returned, with line endings) and stop at the first record.
    fn open(path: &Path, first_record: &mut Vec<u8>) -> io::Result<(Input, Vec<Vec<u8>>)> {
        let mut lines = open_lines(path)?;
        let mut header = Vec::new();
        first_record.clear();
        loop {
            let mut line = Vec::new();
            if lines.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line.starts_with(b"#") {
                header.push(line);
            } else {
                *first_record = line;
                break;
            }
        }
        if !header.last().is_some_and(|l| l.starts_with(b"#CHROM\t")) {
            return Err(invalid(format!("{}: no #CHROM line before the records", path.display())));
        }
        let input = Input { name: path.display().to_string(), lines, next: None, last: None };
        Ok((input, header))
    }

    fn set(&mut self, line: Vec<u8>, contigs: &[Vec<u8>]) -> io::Result<()> {
        if line.is_empty() {
            self.next = None;
            return Ok(());
        }
        if !line.ends_with(b"\n") {
            return Err(invalid(format!("{}: last record has no line ending", self.name)));
        }
        let mut fields = line.split(|&b| b == b'\t');
        let chrom = fields.next().unwrap_or_default();
        let pos = fields
            .next()
            .and_then(|p| std::str::from_utf8(p).ok())
            .and_then(|p| p.parse::<u64>().ok())
            .ok_or_else(|| invalid(format!("{}: record without a position", self.name)))?;
        let contig = contigs.iter().position(|c| c == chrom).ok_or_else(|| {
            invalid(format!(
                "{}: contig {} is not in the header",
                self.name,
                String::from_utf8_lossy(chrom)
            ))
        })?;
        let key = (contig, pos);
        if self.last.is_some_and(|last| key < last) {
            return Err(invalid(format!("{}: records are not sorted", self.name)));
        }
        self.last = Some(key);
        self.next = Some((line, key));
        Ok(())
    }

    fn advance(&mut self, contigs: &[Vec<u8>]) -> io::Result<()> {
        let mut line = Vec::new();
        self.lines.read_until(b'\n', &mut line)?;
        self.set(line, contigs)
    }
}

fn is_command_line(line: &[u8]) -> bool {
    line.starts_with(b"##GATKCommandLine=")
}

/// Contig names of the header's `##contig=<ID=...>` lines, in order.
fn header_contigs(header: &[Vec<u8>]) -> Vec<Vec<u8>> {
    header
        .iter()
        .filter_map(|l| l.strip_prefix(b"##contig=<ID="))
        .map(|rest| rest.split(|&b| b == b',' || b == b'>').next().unwrap_or_default().to_vec())
        .collect()
}

/// `MergeVcfs` for the VCFs of a scattered Mutect2. Returns the number of records written.
pub fn merge_vcfs(inputs: &[PathBuf], sequence_dictionary: Option<&Path>, output: &Path) -> io::Result<u64> {
    if inputs.is_empty() {
        return Err(invalid("no input VCF"));
    }
    let mut opened = Vec::new();
    let mut headers = Vec::new();
    let mut firsts = Vec::new();
    for path in inputs {
        let mut first = Vec::new();
        let (input, header) = Input::open(path, &mut first)?;
        opened.push(input);
        headers.push(header);
        firsts.push(first);
    }
    let header = &headers[0];
    if !header.iter().any(|l| l.starts_with(b"##source=Mutect2")) {
        return Err(invalid(format!(
            "{}: not a Mutect2 VCF; only Mutect2's scatter outputs are supported",
            inputs[0].display()
        )));
    }
    let plain = |h: &[Vec<u8>]| h.iter().filter(|l| !is_command_line(l)).cloned().collect::<Vec<_>>();
    let reference = plain(header);
    for (path, other) in inputs.iter().zip(&headers).skip(1) {
        if plain(other) != reference {
            return Err(invalid(format!(
                "{}: its header differs from the first input's in more than ##GATKCommandLine; \
                 merging different headers is not supported",
                path.display()
            )));
        }
    }
    let contigs = header_contigs(header);
    if contigs.is_empty() {
        return Err(invalid("the header has no ##contig lines"));
    }
    if let Some(dictionary) = sequence_dictionary {
        let text = std::fs::read_to_string(dictionary)?;
        let names: Vec<Vec<u8>> = text
            .lines()
            .filter(|l| l.starts_with("@SQ\t"))
            .filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("SN:")))
            .map(|s| s.as_bytes().to_vec())
            .collect();
        if names != contigs {
            return Err(invalid(format!(
                "{}: its sequences are not the header's contigs; a dictionary that changes the \
                 header is not supported",
                dictionary.display()
            )));
        }
    }
    for (input, first) in opened.iter_mut().zip(firsts) {
        input.set(first, &contigs)?;
    }

    let mut out = ParallelBgzfWriter::new(File::create(output)?, 1);
    for line in header {
        out.write_all(line)?;
    }
    let mut written = 0u64;
    loop {
        let mut best: Option<(usize, (usize, u64))> = None;
        for (i, input) in opened.iter().enumerate() {
            let Some((_, key)) = &input.next else { continue };
            match best {
                Some((b, k)) if *key == k => {
                    return Err(invalid(format!(
                        "{} and {} both have a record at {}:{}; their order in Picard's merge \
                         is not confirmed",
                        opened[b].name,
                        input.name,
                        String::from_utf8_lossy(&contigs[key.0]),
                        key.1
                    )));
                }
                Some((_, k)) if *key > k => {}
                _ => best = Some((i, *key)),
            }
        }
        let Some((i, _)) = best else { break };
        let (line, _) = opened[i].next.take().expect("chosen input has a record");
        out.write_all(&line)?;
        written += 1;
        opened[i].advance(&contigs)?;
    }
    out.finish()?;
    Ok(written)
}

/// `MergeMutectStats`: the sum of `callable` over the inputs, printed as Java prints a double.
pub fn merge_mutect_stats(inputs: &[PathBuf], output: &Path) -> io::Result<f64> {
    if inputs.is_empty() {
        return Err(invalid("no stats file"));
    }
    let mut callable = 0.0f64;
    for path in inputs {
        let text = std::fs::read_to_string(path)?;
        let mut lines = text.lines();
        if lines.next() != Some("statistic\tvalue") {
            return Err(invalid(format!("{}: not a Mutect2 stats file", path.display())));
        }
        let mut seen = false;
        for line in lines {
            let Some(value) = line.strip_prefix("callable\t") else {
                return Err(invalid(format!(
                    "{}: statistic {:?} is not supported, only callable",
                    path.display(),
                    line.split('\t').next().unwrap_or_default()
                )));
            };
            // Decimal forms only: the forms Java's and Rust's parsers both read, the same way.
            let decimal = !value.is_empty()
                && value.bytes().all(|b| b.is_ascii_digit() || matches!(b, b'.' | b'E' | b'-'));
            let value: f64 = match value.parse() {
                Ok(v) if decimal && !seen => v,
                _ => return Err(invalid(format!("{}: bad callable line {line:?}", path.display()))),
            };
            seen = true;
            callable += value;
        }
        if !seen {
            return Err(invalid(format!("{}: no callable statistic", path.display())));
        }
    }
    std::fs::write(
        output,
        format!("statistic\tvalue\ncallable\t{}\n", crate::math::javafmt::double_to_string(callable)),
    )?;
    Ok(callable)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n##GATKCommandLine=<ID=Mutect2,CommandLine=\"X\">\n##contig=<ID=chr1,length=100>\n##contig=<ID=chr2,length=100>\n##source=Mutect2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nomics-merge-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn vcf(d: &Path, name: &str, header: &str, records: &str) -> PathBuf {
        let p = d.join(name);
        std::fs::write(&p, format!("{header}{records}")).unwrap();
        p
    }

    fn text(path: &Path) -> String {
        let mut s = String::new();
        GzReader::open(path).unwrap().read_to_string(&mut s).unwrap();
        s
    }

    #[test]
    fn records_follow_contig_order_whatever_the_input_order() {
        let d = dir("order");
        let a = vcf(&d, "a.vcf", HEADER, "chr2\t5\t.\tA\tC\t.\t.\tX\n");
        let other = HEADER.replace("\"X\"", "\"Y\"");
        let b = vcf(&d, "b.vcf", &other, "chr1\t9\t.\tA\tC\t.\t.\tX\nchr2\t4\t.\tA\tG\t.\t.\tX\n");
        let out = d.join("out.vcf.gz");
        assert_eq!(merge_vcfs(&[a, b], None, &out).unwrap(), 3);
        // The first input's header, command line included; records by contig then position.
        assert_eq!(
            text(&out),
            format!("{HEADER}chr1\t9\t.\tA\tC\t.\t.\tX\nchr2\t4\t.\tA\tG\t.\t.\tX\nchr2\t5\t.\tA\tC\t.\t.\tX\n")
        );
        std::fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn what_is_not_confirmed_is_refused() {
        let d = dir("refuse");
        let a = vcf(&d, "a.vcf", HEADER, "chr1\t9\t.\tA\tC\t.\t.\tX\n");
        let out = d.join("out.vcf.gz");
        let same_position = vcf(&d, "b.vcf", HEADER, "chr1\t9\t.\tA\tG\t.\t.\tX\n");
        assert!(merge_vcfs(&[a.clone(), same_position], None, &out).is_err());
        let other_header = vcf(&d, "c.vcf", &HEADER.replace("length=100>\n##source", "length=101>\n##source"), "");
        assert!(merge_vcfs(&[a.clone(), other_header], None, &out).is_err());
        let unsorted = vcf(&d, "d.vcf", HEADER, "chr2\t5\t.\tA\tC\t.\t.\tX\nchr1\t1\t.\tA\tC\t.\t.\tX\n");
        assert!(merge_vcfs(&[a.clone(), unsorted], None, &out).is_err());
        let not_mutect = vcf(&d, "e.vcf", &HEADER.replace("##source=Mutect2\n", ""), "");
        assert!(merge_vcfs(&[not_mutect], None, &out).is_err());
        let dict = d.join("x.dict");
        std::fs::write(&dict, "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:100\n").unwrap();
        assert!(merge_vcfs(std::slice::from_ref(&a), Some(&dict), &out).is_err());
        std::fs::write(&dict, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:100\n").unwrap();
        assert_eq!(merge_vcfs(&[a], Some(&dict), &out).unwrap(), 1);
        std::fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn stats_sum_prints_as_java_does() {
        let d = dir("stats");
        let (a, b, out) = (d.join("a"), d.join("b"), d.join("out"));
        std::fs::write(&a, "statistic\tvalue\ncallable\t2.3995465E7\n").unwrap();
        std::fs::write(&b, "statistic\tvalue\ncallable\t9.640902E7\n").unwrap();
        merge_mutect_stats(&[a.clone(), b], &out).unwrap();
        assert_eq!(std::fs::read_to_string(&out).unwrap(), "statistic\tvalue\ncallable\t1.20404485E8\n");
        let other = d.join("c");
        std::fs::write(&other, "statistic\tvalue\ncallable\t1.0\nelse\t2.0\n").unwrap();
        assert!(merge_mutect_stats(&[a.clone(), other.clone()], &out).is_err());
        std::fs::write(&other, "statistic\tvalue\ncallable\tNaN\n").unwrap();
        assert!(merge_mutect_stats(&[a, other], &out).is_err());
        std::fs::remove_dir_all(d).unwrap();
    }
}
