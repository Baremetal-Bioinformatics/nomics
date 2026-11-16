//! FastQC 0.12.1 default FASTQ text reports. See docs/spec/qc-fastqc.md.
mod contaminants;
mod distribution;
mod duplication;
mod float;
mod java_map;

use float::java;
use nomics_io::gzip::GzReader;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

// FastQC v0.12.1 Configuration/adapter_list.txt, in upstream order.
const ADAPTERS: [(&str, &[u8]); 6] = [
    ("Illumina Universal Adapter", b"AGATCGGAAGAG"),
    ("Illumina Small RNA 3' Adapter", b"TGGAATTCTCGG"),
    ("Illumina Small RNA 5' Adapter", b"GATCGTCGGACT"),
    ("Nextera Transposase Sequence", b"CTGTCTCTTATA"),
    ("PolyA", b"AAAAAAAAAAAA"),
    ("PolyG", b"GGGGGGGGGGGG"),
];
fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}

/// Write `fastqc_data.txt` and `summary.txt` into `outdir`, one input/mate per call.
pub fn fastqc(input: &Path, outdir: &Path, input_name: Option<&str>) -> io::Result<()> {
    let name = input_name
        .or_else(|| input.file_name()?.to_str())
        .ok_or_else(|| invalid("input filename must be UTF-8"))?;
    let mut file = File::open(input)?;
    let mut magic = [0; 2];
    let n = file.read(&mut magic)?;
    let reader: Box<dyn Read> = if n == 2 && magic == [31, 139] {
        Box::new(GzReader::open(input)?)
    } else {
        Box::new(File::open(input)?)
    };
    let [data, summary] = reports(BufReader::with_capacity(1 << 20, reader), name)?;
    fs::create_dir_all(outdir)?;
    fs::write(outdir.join("fastqc_data.txt"), data)?;
    fs::write(outdir.join("summary.txt"), summary)
}

/// Reusable four-line buffers, as in nomics-prep's private `Fastq` reader;
/// decoding is the existing nomics-io reader. FastQC permits blank inter-record lines.
pub fn reports(mut reader: impl BufRead, name: &str) -> io::Result<[String; 2]> {
    let mut lines: [Vec<u8>; 4] = std::array::from_fn(|_| Vec::new());
    let mut stats = Stats::default();
    let mut colorspace = None;
    loop {
        for (i, line) in lines.iter_mut().enumerate() {
            loop {
                line.clear();
                if reader.read_until(b'\n', line)? == 0 {
                    if i == 0 {
                        return stats.render(name);
                    }
                    return Err(invalid("truncated FASTQ record"));
                }
                while matches!(line.last(), Some(b'\n' | b'\r')) {
                    line.pop();
                }
                if i != 0 || !line.is_empty() {
                    break;
                }
            }
        }
        if !lines[0].starts_with(b"@") || !lines[2].starts_with(b"+") {
            return Err(invalid("malformed FASTQ record"));
        }
        lines[1].make_ascii_uppercase();
        let cs = *colorspace.get_or_insert_with(|| {
            lines[1].len() > 1
                && b"GATCN".contains(&lines[1][0])
                && lines[1][1..].iter().all(|b| b".0123456".contains(b))
        });
        if cs && !lines[1].is_empty() {
            // FastQFile.convertColorspaceToBases: missing colour poisons the rest.
            let seq = &mut lines[1];
            let mut base = seq[0];
            if !b"ACGT".contains(&base) {
                return Err(invalid("invalid colorspace primer"));
            }
            for i in 1..seq.len() {
                base = if base == b'N' || b".456".contains(&seq[i]) {
                    b'N'
                } else if let Some(color) = seq[i].checked_sub(b'0').filter(|&v| v < 4) {
                    let b = b"ACGT"
                        .iter()
                        .position(|&v| v == base)
                        .ok_or_else(|| invalid("invalid colorspace primer"))?;
                    b"ACGT"[b ^ color as usize]
                } else {
                    return Err(invalid("invalid colorspace base"));
                };
                seq[i - 1] = base;
            }
            seq.pop();
        }
        stats.colorspace = cs;
        stats.add(&lines[0], &lines[1], &lines[3])?;
    }
}

#[derive(Clone)]
struct Quality {
    counts: [u64; 150],
    total: u64,
}
impl Default for Quality {
    fn default() -> Self {
        Self {
            counts: [0; 150],
            total: 0,
        }
    }
}
impl Quality {
    fn add(&mut self, q: u8) {
        self.counts[q as usize] += 1;
        self.total += 1;
    }
    fn min(&self) -> usize {
        self.counts.iter().position(|&v| v > 0).unwrap_or(1000)
    }
    fn mean(&self, offset: usize) -> f64 {
        let mut total = 0;
        let mut n = 0;
        for (q, &v) in self.counts.iter().enumerate().skip(offset) {
            total += v * (q - offset) as u64;
            n += v;
        }
        total as f64 / n as f64
    }
    fn percentile(&self, offset: usize, p: u64) -> f64 {
        let threshold = self.total * p / 100;
        let mut n = 0;
        for (q, &v) in self.counts.iter().enumerate().skip(offset) {
            n += v;
            if n >= threshold {
                return (q - offset) as f64;
            }
        }
        -1.0
    }
}
fn encoding(low: usize) -> io::Result<(usize, &'static str)> {
    match low {
        33..64 => Ok((33, "Sanger / Illumina 1.9")),
        65 => Ok((64, "Illumina 1.3")),
        64..=126 => Ok((64, "Illumina 1.5")),
        _ => Err(invalid("no known FASTQ quality encoding")),
    }
}

#[derive(Default)]
struct Stats {
    count: u64,
    total_bases: u64,
    colorspace: bool,
    bases: Vec<[u64; 6]>, // G A T C N other
    quals: Vec<Quality>,
    lengths: Vec<u64>,
    averages: BTreeMap<usize, u64>,
    gc: distribution::Gc,
    dup: duplication::Dup,
    adapter: [Vec<u64>; 6],
    longest: usize,
    tiles: BTreeMap<i32, Vec<Quality>>,
    tile_ignore: bool,
    tile_count: u64,
    tile_split: Option<usize>,
    tile_length: usize,
}
impl Stats {
    fn add(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) -> io::Result<()> {
        if qual.iter().any(|&q| !(33..=126).contains(&q)) {
            return Err(invalid("invalid quality character"));
        }
        self.count += 1;
        self.total_bases += seq.len() as u64;
        self.bases.resize(self.bases.len().max(seq.len()), [0; 6]);
        if self.lengths.len() <= seq.len() {
            self.lengths.resize(seq.len() + 1, 0);
        }
        self.lengths[seq.len()] += 1;
        for (i, &b) in seq.iter().enumerate() {
            self.bases[i][match b {
                b'G' => 0,
                b'A' => 1,
                b'T' => 2,
                b'C' => 3,
                b'N' => 4,
                _ => 5,
            }] += 1;
        }
        self.quals
            .resize_with(self.quals.len().max(qual.len()), Quality::default);
        let mut sum = 0;
        for (i, &q) in qual.iter().enumerate() {
            self.quals[i].add(q);
            sum += q as usize;
        }
        if !qual.is_empty() {
            *self.averages.entry(sum / qual.len()).or_default() += 1;
        }
        self.gc.add(seq);
        self.dup.add(seq);
        // Store only first-hit increments; expansion inherits the previous tail.
        // This is equivalent to Adapter.expandLengthTo's cumulative counts.
        if seq.len() > self.longest && seq.len() > 12 {
            self.longest = seq.len();
            for v in &mut self.adapter {
                let tail = v.last().copied().unwrap_or(0);
                v.resize(self.longest - 11, tail);
            }
        }
        for (a, (_, adapter)) in ADAPTERS.iter().enumerate() {
            if let Some(i) = seq.windows(adapter.len()).position(|v| v == *adapter) {
                for v in self.adapter[a].iter_mut().skip(i) {
                    *v += 1;
                }
            }
        }
        self.add_tile(id, qual);
        Ok(())
    }
    fn add_tile(&mut self, id: &[u8], qual: &[u8]) {
        if self.tile_ignore || qual.is_empty() {
            return;
        }
        self.tile_count += 1;
        if self.tile_count > 10000 && !self.tile_count.is_multiple_of(10) {
            return;
        }
        // String.split(":") discards empty trailing fields.
        let end = id.iter().rposition(|&b| b != b':').map_or(0, |i| i + 1);
        let fields = id[..end].split(|&b| b == b':');
        let count = fields.clone().count();
        let index = *self.tile_split.get_or_insert(if count >= 7 {
            4
        } else if count >= 5 {
            2
        } else {
            usize::MAX
        });
        let tile = fields
            .clone()
            .nth(index)
            .and_then(|v| std::str::from_utf8(v).ok()?.parse::<i32>().ok());
        let Some(tile) = tile else {
            self.tile_ignore = true;
            return;
        };
        if qual.len() > self.tile_length {
            self.tile_length = qual.len();
            for v in self.tiles.values_mut() {
                v.resize_with(qual.len(), Quality::default);
            }
        }
        if !self.tiles.contains_key(&tile) && self.tiles.len() > 2500 {
            self.tile_ignore = true;
            self.tiles.clear();
            return;
        }
        let counts = self
            .tiles
            .entry(tile)
            .or_insert_with(|| vec![Quality::default(); self.tile_length]);
        for (q, &v) in counts.iter_mut().zip(qual) {
            q.add(v);
        }
    }
    fn render(self, name: &str) -> io::Result<[String; 2]> {
        let low = self.quals.iter().map(Quality::min).min().unwrap_or(126);
        let (offset, enc) = encoding(low)?;
        let mut report = Report {
            data: "##FastQC\t0.12.1\n".into(),
            summary: String::new(),
            name,
        };
        let minlen = self.lengths.iter().position(|&v| v > 0).unwrap_or(0);
        let maxlen = self.lengths.len().saturating_sub(1);
        let length = label(minlen, maxlen);
        let mut totals = [0u64; 6];
        for b in &self.bases {
            for (a, b) in totals.iter_mut().zip(b) {
                *a += b;
            }
        }
        let called: u64 = totals[..4].iter().sum();
        let gc = ((totals[0] + totals[3]) * 100)
            .checked_div(called)
            .unwrap_or(0);
        let basic = format!(
            "#Measure\tValue\nFilename\t{}\nFile type\t{}\nEncoding\t{enc}\nTotal Sequences\t{}\nTotal Bases\t{}\nSequences flagged as poor quality\t0\nSequence length\t{length}\n%GC\t{gc}\n",
            name.replacen("stdin:", "", 1),
            if self.count == 0 {
                "null"
            } else if self.colorspace {
                "Colorspace converted to bases"
            } else {
                "Conventional base calls"
            },
            self.count,
            format_length(self.total_bases)
        );
        report.module("Basic Statistics", 0, &basic);
        let quality_groups = groups(self.quals.len());
        if !self.quals.is_empty() {
            let mut body = "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n".to_owned();
            let mut status = 0;
            for &(a, b) in &quality_groups {
                let mean = group_mean(&self.quals, a, b, offset);
                let pct = |p| {
                    let mut sum = 0.0;
                    let mut count = 0;
                    for q in &self.quals[a - 1..b] {
                        if q.total > 100 {
                            sum += q.percentile(offset, p);
                            count += 1;
                        }
                    }
                    if count > 0 {
                        sum / count as f64
                    } else {
                        f64::NAN
                    }
                };
                let median = pct(50);
                let lower = pct(25);
                // Configuration/limits.txt: strict lower/median comparisons.
                if !lower.is_nan() {
                    status = status.max(if lower < 5.0 || median < 20.0 {
                        2
                    } else if lower < 10.0 || median < 25.0 {
                        1
                    } else {
                        0
                    });
                }
                write!(body, "{}", label(a, b)).unwrap();
                for v in [mean, median, lower, pct(75), pct(10), pct(90)] {
                    write!(body, "\t{}", java(v)).unwrap();
                }
                body.push('\n');
            }
            report.module("Per base sequence quality", status, &body);
        }
        if !self.tile_ignore && self.tile_length > 0 {
            let tile_offset = encoding(
                self.tiles
                    .values()
                    .flatten()
                    .map(Quality::min)
                    .min()
                    .unwrap_or(126),
            )?
            .0;
            let g = groups(self.tile_length);
            let mut means: Vec<_> = self
                .tiles
                .values()
                .map(|q| {
                    g.iter()
                        .map(|&(a, b)| group_mean(q, a, b, tile_offset))
                        .collect::<Vec<_>>()
                })
                .collect();
            let mut body = "#Tile\tBase\tMean\n".to_owned();
            let mut max = 0.0f64;
            for i in 0..g.len() {
                let mut sum = 0.0;
                for m in &means {
                    sum += m[i];
                }
                sum /= means.len() as f64;
                for m in &mut means {
                    m[i] -= sum;
                    if m[i].abs() > max {
                        max = m[i].abs();
                    }
                }
            }
            for ((tile, _), m) in self.tiles.iter().zip(means) {
                for ((a, b), v) in g.iter().zip(m) {
                    writeln!(body, "{tile}\t{}\t{}", label(*a, *b), java(v)).unwrap();
                }
            }
            report.module("Per tile sequence quality", above(max, 5.0, 10.0), &body);
        }
        if let (Some((&lo, _)), Some((&hi, _))) = (
            self.averages.first_key_value(),
            self.averages.last_key_value(),
        ) {
            let mut body = "#Quality\tCount\n".to_owned();
            let mut best = 0;
            let mut mode = 0;
            for q in lo..=hi {
                let n = *self.averages.get(&q).unwrap_or(&0);
                if n > best {
                    best = n;
                    mode = q - offset;
                }
                writeln!(body, "{}\t{}", q - offset, java(n as f64)).unwrap();
            }
            report.module(
                "Per sequence quality scores",
                if mode <= 20 {
                    2
                } else if mode <= 27 {
                    1
                } else {
                    0
                },
                &body,
            );
        }
        let mut body = "#Base\tG\tA\tT\tC\n".to_owned();
        let mut ns = "#Base\tN-Count\n".to_owned();
        let mut status = 0;
        let mut nstatus = 0;
        for (a, b) in groups(self.bases.len()) {
            let mut counts = [0u64; 6];
            for c in &self.bases[a - 1..b] {
                for (d, v) in counts.iter_mut().zip(c) {
                    *d += v;
                }
            }
            let called: u64 = counts[..4].iter().sum();
            let pct = counts[..4]
                .iter()
                .map(|&v| (v as f64 / called as f64) * 100.0)
                .collect::<Vec<_>>();
            for diff in [(pct[0] - pct[3]).abs(), (pct[1] - pct[2]).abs()] {
                status = status.max(above(diff, 10.0, 20.0));
            }
            write!(body, "{}", label(a, b)).unwrap();
            for p in pct {
                write!(body, "\t{}", java(p)).unwrap();
            }
            body.push('\n');
            let n = 100.0 * (counts[4] as f64 / counts.iter().sum::<u64>() as f64);
            nstatus = nstatus.max(above(n, 5.0, 20.0));
            writeln!(ns, "{}\t{}", label(a, b), java(n)).unwrap();
        }
        report.module("Per base sequence content", status, &body);
        let (body, status) = self.gc.render();
        report.module("Per sequence GC content", status, &body);
        report.module("Per base N content", nstatus, &ns);
        report.module(
            "Sequence Length Distribution",
            if self.lengths.first().copied().unwrap_or(0) > 0 {
                2
            } else if self.lengths.iter().filter(|&&n| n > 0).count() > 1 {
                1
            } else {
                0
            },
            &distribution::lengths(&self.lengths),
        );
        let (dup, status, overrep, overstatus) = self.dup.render(self.count);
        report.module("Sequence Duplication Levels", status, &dup);
        report.module("Overrepresented sequences", overstatus, &overrep);
        let mut body = String::new();
        let mut status = 0;
        if self.longest >= 12 {
            body.push_str("#Position");
            for (name, _) in ADAPTERS {
                write!(body, "\t{name}").unwrap();
            }
            body.push('\n');
            for (a, b) in groups(self.longest - 11) {
                write!(body, "{}", label(a, b)).unwrap();
                for positions in &self.adapter {
                    let mut total = 0.0;
                    for &n in &positions[a - 1..b] {
                        total += (n as f64 * 100.0) / self.count as f64;
                    }
                    total /= (b - a + 1) as f64;
                    status = status.max(above(total, 5.0, 10.0));
                    write!(body, "\t{}", java(total)).unwrap();
                }
                body.push('\n');
            }
        } else {
            status = 1;
        }
        report.module("Adapter Content", status, &body);
        Ok([report.data, report.summary])
    }
}
fn group_mean(q: &[Quality], a: usize, b: usize, offset: usize) -> f64 {
    let mut sum = 0.0;
    let mut n = 0;
    for q in &q[a - 1..b] {
        if q.total > 0 {
            sum += q.mean(offset);
            n += 1;
        }
    }
    if n > 0 { sum / n as f64 } else { 0.0 }
}
fn above(v: f64, warn: f64, fail: f64) -> u8 {
    if v > fail {
        2
    } else if v > warn {
        1
    } else {
        0
    }
}
fn label(a: usize, b: usize) -> String {
    if a == b {
        a.to_string()
    } else {
        format!("{a}-{b}")
    }
}
fn groups(n: usize) -> Vec<(usize, usize)> {
    let mut interval = 1;
    if n > 75 {
        let mut multiplier = 1;
        'outer: loop {
            for base in [2, 5, 10] {
                interval = base * multiplier;
                if 9 + (n - 9).div_ceil(interval) < 75 {
                    break 'outer;
                }
            }
            multiplier *= 10;
        }
    }
    let mut out = Vec::new();
    let mut a = 1;
    while a <= n {
        let b = if a < 10 {
            a
        } else if a == 10 && interval > 10 {
            interval - 1
        } else {
            a + interval - 1
        };
        let b = b.min(n);
        out.push((a, b));
        a = b + 1;
    }
    out
}
fn format_length(n: u64) -> String {
    let (scale, unit) = if n >= 1_000_000_000 {
        (1e9, "Gbp")
    } else if n >= 1_000_000 {
        (1e6, "Mbp")
    } else if n >= 1000 {
        (1e3, "kbp")
    } else {
        (1.0, "bp")
    };
    let raw = java(n as f64 / scale);
    let dot = raw.find('.').unwrap_or(raw.len() - 1);
    let end = if raw.as_bytes().get(dot + 1).is_some_and(|&c| c != b'0') {
        dot + 2
    } else {
        dot
    };
    format!("{} {unit}", &raw[..end])
}
struct Report<'a> {
    data: String,
    summary: String,
    name: &'a str,
}
impl Report<'_> {
    fn module(&mut self, name: &str, status: u8, body: &str) {
        let status = ["pass", "warn", "fail"][status as usize];
        writeln!(self.data, ">>{name}\t{status}").unwrap();
        self.data.push_str(body);
        self.data.push_str(">>END_MODULE\n");
        writeln!(
            self.summary,
            "{}\t{name}\t{}",
            status.to_uppercase(),
            self.name
        )
        .unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_qualities_and_encoding_are_inferred_from_characters() {
        // PerBaseQualityScores excludes positions with <=100 observations
        // from percentiles, but still reports their means.
        let record = b"@r\nACGT\n+\n!!!!\n";
        let [data, summary] = reports(record.as_slice(), "staged.gz").unwrap();
        assert!(data.contains("1\t0.0\tNaN\tNaN\tNaN\tNaN\tNaN\n"));
        assert!(summary.starts_with("PASS\tBasic Statistics\tstaged.gz\n"));
        assert!(!data.contains(">>Per tile"));
        assert!(!data.contains(">>Kmer"));
        let mut reads = Vec::new();
        for _ in 0..101 {
            reads.extend_from_slice(record);
        }
        let [data, _] = reports(reads.as_slice(), "x").unwrap();
        assert!(data.contains(">>Per base sequence quality\tfail\n"));
        assert!(data.contains("1\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n"));
        assert_eq!(encoding(64).unwrap(), (64, "Illumina 1.5"));
        assert_eq!(encoding(65).unwrap(), (64, "Illumina 1.3"));
        assert!(reports(b"@r\nA\n+\n".as_slice(), "x").is_err());
    }
    #[test]
    fn linear_bins_have_first_nine_and_fewer_than_75_groups() {
        assert_eq!(groups(75).len(), 75);
        assert_eq!(&groups(151)[8..11], [(9, 9), (10, 14), (15, 19)]);
        assert_eq!(groups(151).last(), Some(&(150, 151)));
        assert_eq!(&groups(1001)[8..11], [(9, 9), (10, 19), (20, 39)]);
    }
    #[test]
    fn tile_sampling_uses_first_10000_then_each_tenth() {
        let mut s = Stats::default();
        for _ in 0..10000 {
            s.add_tile(b"@a:1:f:1:1101:2:3", b"!!");
        }
        for _ in 0..9 {
            s.add_tile(b"@a:1:f:1:1102:2:3", b"II");
        }
        assert!(!s.tiles.contains_key(&1102));
        s.add_tile(b"@a:1:f:1:1102:2:3", b"II");
        assert_eq!(s.tiles[&1102][0].total, 1);
        for _ in 0..10 {
            s.add_tile(b"@invalid", b"II");
        }
        assert!(s.tile_ignore);
    }
}
