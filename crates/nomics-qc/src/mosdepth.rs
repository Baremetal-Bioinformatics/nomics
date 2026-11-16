//! mosdepth 0.3.14's BED/window reports. Source ledger: docs/spec/qc-mosdepth.md.
pub mod csi;
use nomics_io::bam::{self, RecordView};
use nomics_io::bgzf::{BgzfWriter, ParallelBgzfReader, Profile};
use nomics_io::deflate::FastDeflater;
use nomics_mem::HostArray;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}
fn path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}
type BedWriter = BgzfWriter<BufWriter<File>, FastDeflater>;
fn writer(prefix: &Path, suffix: &str) -> io::Result<BedWriter> {
    Ok(BgzfWriter::new(
        BufWriter::new(File::create(path(prefix, suffix))?),
        FastDeflater::default(),
        Profile {
            block_payload: 65_280,
        },
    ))
}
#[derive(Debug)]
struct Region {
    start: usize,
    stop: usize,
    name: String,
}
fn bed_integer(s: &str) -> io::Result<i64> {
    let digits = s.strip_prefix(['+', '-']).unwrap_or(s);
    if !digits.as_bytes().first().is_some_and(u8::is_ascii_digit)
        || !digits.bytes().all(|c| c.is_ascii_digit() || c == b'_')
    {
        return Err(invalid("invalid BED integer"));
    }
    s.replace('_', "")
        .parse()
        .map_err(|_| invalid("BED integer outside signed 64-bit range"))
}
fn bed(file: &Path) -> io::Result<HashMap<String, Vec<Region>>> {
    let mut magic = [0; 2];
    let n = File::open(file)?.read(&mut magic)?;
    let r: Box<dyn BufRead> = if n == 2 && magic == [31, 139] {
        Box::new(BufReader::new(nomics_io::gzip::GzReader::open(file)?))
    } else {
        Box::new(BufReader::new(File::open(file)?))
    };
    let mut regions: HashMap<String, Vec<Region>> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        // hts_getline's >0 loop terminates at an empty line.
        if line.is_empty() {
            break;
        }
        if line.starts_with('#') || line.starts_with("track ") {
            continue;
        }
        let f: Vec<_> = line.splitn(6, '\t').collect();
        if f.len() < 3 {
            eprintln!("mosdepth: skipping bad BED line: {line}");
            continue;
        }
        let start = bed_integer(f[1])?;
        let stop = bed_integer(f[2])?;
        if start > stop {
            return Err(invalid("BED start > end"));
        }
        regions.entry(f[0].into()).or_default().push(Region {
            start: start as u32 as usize,
            stop: stop as u32 as usize,
            name: f
                .get(3)
                .unwrap_or(&"")
                .trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{b}' | '\u{c}'))
                .into(),
        });
    }
    for regions in regions.values_mut() {
        regions.sort_by_key(|r| r.start);
    }
    Ok(regions)
}
#[derive(Clone, Copy)]
struct Stats {
    length: u64,
    bases: u64,
    min: u32,
    max: u32,
}
impl Default for Stats {
    fn default() -> Self {
        Self {
            length: 0,
            bases: 0,
            min: u32::MAX,
            max: 0,
        }
    }
}
impl Stats {
    fn add(&mut self, values: &[u32]) {
        for &d in values {
            self.length += 1;
            self.bases = self.bases.wrapping_add((d as i32) as u64);
            self.min = self.min.min(d);
            self.max = self.max.max(d);
        }
    }
    fn merge(&mut self, other: Self) {
        self.length += other.length;
        self.bases = self.bases.wrapping_add(other.bases);
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
    }
    fn write(self, name: &str, w: &mut impl Write) -> io::Result<()> {
        let mean = if self.length > 0 {
            self.bases as f64 / self.length as f64
        } else {
            0.0
        };
        writeln!(
            w,
            "{name}\t{}\t{}\t{mean:.2}\t{}\t{}",
            self.length,
            self.bases,
            if self.min == u32::MAX { 0 } else { self.min },
            self.max
        )
    }
}
fn increment(hist: &mut Vec<u64>, values: &[u32]) {
    for &d in values {
        // Nim first clamps >400000, then ignores signed-negative coverage.
        let d = d as i32;
        let d = if d > 400_000 { 399_990 } else { d };
        if d < 0 {
            continue;
        }
        let d = d as usize;
        if d >= hist.len() - 1 {
            hist.resize(d + 10, 0);
        }
        hist[d] += 1;
    }
}
fn distribution(name: &str, hist: &[u64], w: &mut impl Write) -> io::Result<()> {
    let sum: u64 = hist.iter().sum();
    if sum == 0 {
        return Ok(());
    }
    let mut cum = 0.0;
    for i in (0..hist.len()).rev() {
        let v = hist[i];
        if i > 300 && v == 0 {
            continue;
        }
        cum += v as f64 / sum as f64;
        if cum < 8e-5 {
            continue;
        }
        writeln!(w, "{name}\t{i}\t{cum:.2}")?;
    }
    Ok(())
}
fn sum_into(total: &mut Vec<u64>, hist: &[u64]) {
    if total.len() < hist.len() {
        total.resize(hist.len(), 0);
    }
    for (a, b) in total.iter_mut().zip(hist) {
        *a += b;
    }
}

#[derive(Clone)]
struct Mate {
    pos: usize,
    stop: usize,
    cigar: Vec<u32>,
}
fn segments(cigar: impl IntoIterator<Item = u32>, start: usize) -> Vec<(usize, usize)> {
    let mut intervals: Vec<(usize, usize)> = Vec::new();
    let mut pos = start;
    for c in cigar {
        let len = (c >> 4) as usize;
        match c & 15 {
            0 | 7 | 8 => {
                if let Some(last) = intervals.last_mut()
                    && last.1 == pos
                {
                    last.1 += len;
                } else {
                    intervals.push((pos, pos + len));
                }
                pos += len;
            }
            2 | 3 => pos += len,
            _ => {}
        }
    }
    intervals
}
fn change(arr: &mut [u32], start: usize, stop: usize, subtract: bool) -> io::Result<()> {
    if start >= arr.len() || stop >= arr.len() {
        return Err(invalid("alignment extends outside chromosome depth array"));
    }
    let v = if subtract { u32::MAX } else { 1 };
    arr[start] = arr[start].wrapping_add(v);
    arr[stop] = arr[stop].wrapping_sub(v);
    Ok(())
}
// htslib expands a soft-clip placeholder from CG:B:I/i before mosdepth
// sees it. Keep this local to QC; the raw BAM codec stays lossless.
fn decoded_cigar(rec: RecordView<'_>) -> io::Result<Vec<u32>> {
    let cigar: Vec<_> = rec.cigar().collect();
    if cigar.first().copied() != Some(((rec.len() as u32) << 4) | 4) {
        return Ok(cigar);
    }
    let mut aux = rec.aux().ok_or_else(|| invalid("truncated BAM sequence"))?;
    while aux.len() >= 3 {
        let ty = aux[2];
        let body = &aux[3..];
        let len = bam::aux_value_len(ty, body).ok_or_else(|| invalid("invalid BAM aux tag"))?;
        if len > body.len() {
            return Err(invalid("truncated BAM aux tag"));
        }
        if aux[..2] == *b"CG" && ty == b'B' && matches!(body[0], b'I' | b'i') {
            let n = u32::from_le_bytes(body[1..5].try_into().unwrap()) as usize;
            if n >= cigar.len() && n < 1 << 29 {
                return Ok(body[5..len]
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|b| u32::from_le_bytes(*b))
                    .collect());
            }
            return Ok(cigar);
        }
        aux = &body[len..];
    }
    Ok(cigar)
}
fn add_read(
    arr: &mut [u32],
    rec: RecordView<'_>,
    seen: &mut HashMap<Vec<u8>, Mate>,
    fast: bool,
) -> io::Result<()> {
    let start = rec.pos() as usize;
    let cigar = decoded_cigar(rec)?;
    let stop = start + bam::reference_length(&cigar).max(1) as usize;
    if !fast && rec.flag() & 2 != 0 && rec.flag() & 2048 == 0 {
        let mpos = rec.mate_pos();
        if rec.ref_id() == rec.mate_ref_id()
            && stop as i64 > i64::from(mpos)
            && (rec.pos() < mpos || (rec.pos() == mpos && !seen.contains_key(rec.read_name())))
        {
            seen.insert(
                rec.read_name().to_vec(),
                Mate {
                    pos: start,
                    stop,
                    cigar: cigar.clone(),
                },
            );
        } else if let Some(mate) = seen.remove(rec.read_name()) {
            // Upstream's single-CIGAR shortcut does not test the opcode or
            // use min(mate.stop, rec.stop), including for contained mates.
            if cigar.len() == 1 && mate.cigar.len() == 1 {
                change(arr, start, mate.stop, true)?;
            } else {
                let a = segments(cigar.iter().copied(), start);
                let b = segments(mate.cigar, mate.pos);
                for &(s, e) in &a {
                    for &(ms, me) in &b {
                        let lo = s.max(ms);
                        let hi = e.min(me);
                        if lo < hi {
                            change(arr, lo, hi, true)?;
                        }
                    }
                }
            }
        }
    }
    if fast {
        change(arr, start, stop, false)?;
    } else {
        for (s, e) in segments(cigar, start) {
            change(arr, s, e, false)?;
        }
    }
    Ok(())
}

/// Emit the two sarek forms. `by` is an integer window or a BED filename.
/// BAM decoding uses nomics-io; no reference sequence is consulted.
pub fn run(
    bam_path: &Path,
    prefix: &Path,
    by: &str,
    no_per_base: bool,
    fast_mode: bool,
    threads: usize,
) -> io::Result<()> {
    let window = if !by.is_empty() && by.bytes().all(|c| c.is_ascii_digit()) {
        Some(by.parse::<i64>().map_err(|_| invalid("invalid window"))? as u32 as usize)
    } else {
        None
    };
    if window == Some(0) {
        return Err(invalid("window must be positive"));
    }
    let mut regions = if window.is_none() {
        Some(bed(Path::new(by))?)
    } else {
        None
    };
    if let Some(parent) = prefix.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let mut reader = ParallelBgzfReader::new(File::open(bam_path)?, threads);
    let header = bam::Header::read(&mut reader)?;
    let max_len = header.refs.iter().map(|r| r.len).max().unwrap_or(0);
    let mut levels = 0;
    while u64::from(max_len) > 1u64 << (14 + 3 * levels) {
        levels += 1;
    }
    let mut base = if no_per_base {
        None
    } else {
        Some(writer(prefix, ".per-base.bed.gz")?)
    };
    let mut region = writer(prefix, ".regions.bed.gz")?;
    let mut summary = BufWriter::new(File::create(path(prefix, ".mosdepth.summary.txt"))?);
    let mut global_dist = BufWriter::new(File::create(path(prefix, ".mosdepth.global.dist.txt"))?);
    let mut region_dist = BufWriter::new(File::create(path(prefix, ".mosdepth.region.dist.txt"))?);
    writeln!(summary, "chrom\tlength\tbases\tmean\tmin\tmax")?;
    let mut global = Stats::default();
    let mut global_region = Stats::default();
    let mut ghist = vec![0; 512];
    let mut rhist = vec![0; 512];
    let mut chist = vec![0; 512];
    let mut crhist = vec![0; 512];
    let mut raw = Vec::new();
    let mut have = bam::read_raw_record(&mut reader, &mut raw)?;
    for (tid, target) in header.refs.iter().enumerate() {
        let name =
            std::str::from_utf8(&target.name).map_err(|_| invalid("non-UTF8 chromosome name"))?;
        let len = target.len as usize;
        let skipped = no_per_base && regions.as_ref().is_some_and(|r| !r.contains_key(name));
        // Unplaced records are outside every indexed chromosome query.
        let mut found = false;
        let mut arr: Option<HostArray<u32>> = None;
        let mut seen = HashMap::new();
        while have {
            let rec = RecordView::new(&raw)?;
            if rec.ref_id() < 0 || rec.ref_id() as usize > tid {
                break;
            }
            if (rec.ref_id() as usize) < tid {
                return Err(invalid("BAM chromosomes out of order"));
            }
            if !skipped && rec.pos() >= 0 && (rec.pos() as usize) < len {
                if !found {
                    arr = Some(HostArray::zeroed(len + 1)?);
                    found = true;
                }
                if rec.flag() & 1796 == 0 {
                    add_read(arr.as_mut().unwrap(), rec, &mut seen, fast_mode)?;
                }
            }
            have = bam::read_raw_record(&mut reader, &mut raw)?;
        }
        if skipped {
            continue;
        }
        if let Some(arr) = &mut arr {
            let mut depth = 0u32;
            for d in arr.iter_mut() {
                depth = depth.wrapping_add(*d);
                *d = depth;
            }
        }
        chist.fill(0);
        crhist.fill(0);
        let values = arr.as_deref().unwrap_or(&[]);
        let rs: Box<dyn Iterator<Item = Region>> = if let Some(size) = window {
            Box::new((0..len).step_by(size).map(move |start| Region {
                start,
                stop: (start + size).min(len),
                name: String::new(),
            }))
        } else {
            Box::new(
                regions
                    .as_mut()
                    .unwrap()
                    .remove(name)
                    .unwrap_or_default()
                    .into_iter(),
            )
        };
        let mut stat = Stats::default();
        for r in rs {
            let start = r.start.min(values.len());
            let stop = r.stop.min(values.len());
            let vals = if start <= stop {
                &values[start..stop]
            } else {
                &[]
            };
            let mut mean = 0.0;
            if found {
                let length = (r.stop - r.start) as f64;
                for &v in vals {
                    mean += (v as i32) as f64 / length;
                }
                stat.add(vals);
                if window.is_some() {
                    let bin = (mean.round() as usize).min(crhist.len() - 1);
                    crhist[bin] += 1;
                } else {
                    increment(&mut crhist, vals);
                }
            }
            if r.name.is_empty() {
                writeln!(region, "{name}\t{}\t{}\t{mean:.2}", r.start, r.stop)?;
            } else {
                writeln!(
                    region,
                    "{name}\t{}\t{}\t{}\t{mean:.2}",
                    r.start, r.stop, r.name
                )?;
            }
        }
        if found {
            increment(&mut chist, &values[..len]);
            let mut cstat = Stats::default();
            cstat.add(&values[..len]);
            cstat.write(name, &mut summary)?;
            stat.write(&format!("{name}_region"), &mut summary)?;
            global.merge(cstat);
            global_region.merge(stat);
        }
        distribution(name, &chist, &mut global_dist)?;
        distribution(name, &crhist, &mut region_dist)?;
        if found {
            sum_into(&mut ghist, &chist);
            sum_into(&mut rhist, &crhist);
        }
        if let Some(base) = &mut base {
            if !found {
                writeln!(base, "{name}\t0\t{len}\t0")?;
            } else {
                let mut start = 0;
                for i in 1..=len {
                    if i == len || values[i] != values[start] {
                        writeln!(base, "{name}\t{start}\t{i}\t{}", values[start])?;
                        start = i;
                    }
                }
            }
        }
    }
    for (name, rs) in regions.unwrap_or_default() {
        eprintln!(
            "mosdepth: BED chromosome {name} ({} regions) not found",
            rs.len()
        );
    }
    global.write("total", &mut summary)?;
    global_region.write("total_region", &mut summary)?;
    distribution("total", &ghist, &mut global_dist)?;
    distribution("total", &rhist, &mut region_dist)?;
    summary.flush()?;
    global_dist.flush()?;
    region_dist.flush()?;
    region.finish()?;
    csi::build(
        &path(prefix, ".regions.bed.gz"),
        &path(prefix, ".regions.bed.gz.csi"),
        levels,
    )?;
    if let Some(base) = base {
        base.finish()?;
        csi::build(
            &path(prefix, ".per-base.bed.gz"),
            &path(prefix, ".per-base.bed.gz.csi"),
            levels,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn record(pos: i32, mpos: i32, cigar: Vec<u32>, flag: u16) -> Vec<u8> {
        let r = bam::Record {
            ref_id: 0,
            pos,
            mapq: 0,
            flag,
            next_ref_id: 0,
            next_pos: mpos,
            tlen: 0,
            read_name: b"pair".to_vec(),
            cigar,
            seq: vec![],
            qual: vec![],
            tags: vec![],
        };
        let mut raw = Vec::new();
        r.encode(bam::BinRule::Sam, &mut raw).unwrap();
        raw
    }
    fn depth(mut d: Vec<u32>) -> Vec<u32> {
        let mut n = 0u32;
        for x in &mut d {
            n = n.wrapping_add(*x);
            *x = n;
        }
        d
    }
    #[test]
    fn cigar_and_overlap_rules() {
        // Insertion and clipping leave a continuous match; D/N create gaps.
        assert_eq!(
            segments(
                [
                    2 << 4,
                    (1 << 4) | 1,
                    3 << 4,
                    (2 << 4) | 2,
                    1 << 4,
                    (3 << 4) | 3,
                    (1 << 4) | 7,
                    (1 << 4) | 8
                ],
                0
            ),
            [(0, 5), (7, 8), (11, 13)]
        );
        let first = record(2, 5, vec![8 << 4], 2);
        let second = record(5, 2, vec![2 << 4], 2);
        let mut arr = vec![0; 13];
        let mut seen = HashMap::new();
        for r in [&first, &second] {
            add_read(&mut arr, RecordView::new(r).unwrap(), &mut seen, false).unwrap();
        }
        // The single-op shortcut subtracts through first.stop=10 even though
        // the contained second ends at 7. This upstream quirk must survive.
        assert_eq!(depth(arr), [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]);
        let mut arr = vec![0; 13];
        for r in [&first, &second] {
            add_read(
                &mut arr,
                RecordView::new(r).unwrap(),
                &mut HashMap::new(),
                true,
            )
            .unwrap();
        }
        assert_eq!(depth(arr), [0, 0, 1, 1, 1, 2, 2, 1, 1, 1, 0, 0, 0]);
        let first = record(1, 3, vec![2 << 4, (1 << 4) | 2, 3 << 4], 2);
        let second = record(3, 1, vec![4 << 4], 2);
        let mut arr = vec![0; 10];
        let mut seen = HashMap::new();
        for r in [&first, &second] {
            add_read(&mut arr, RecordView::new(r).unwrap(), &mut seen, false).unwrap();
        }
        assert_eq!(depth(arr), [0, 1, 1, 1, 1, 1, 1, 0, 0, 0]);
    }
    #[test]
    fn histogram_clamp_and_empty_stats() {
        let mut h = vec![0; 512];
        increment(&mut h, &[511, 512, 400_000, 400_001, u32::MAX]);
        assert_eq!(h[511], 1);
        assert_eq!(h[512], 1);
        assert_eq!(h[400_000], 1);
        assert_eq!(h[399_990], 1);
        assert_eq!(h.iter().sum::<u64>(), 4);
        let mut w = Vec::new();
        Stats::default().write("empty", &mut w).unwrap();
        assert_eq!(w, b"empty\t0\t0\t0.00\t0\t0\n");
        // No reads => no distribution row, including at zero depth.
        let mut w = Vec::new();
        distribution("empty", &[0; 512], &mut w).unwrap();
        assert!(w.is_empty());
    }
    #[test]
    fn bed_stable_sort_names_and_empty_line_termination() {
        let dir = std::env::temp_dir().join(format!("nomics-mosdepth-bed-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let p = dir.join("regions.bed");
        fs::write(&p,"track name=test\n#comment\nchrA\t7\t9\t last \nchrA\t2\t5\tone\nchrA\t2\t4\ttwo\textra\nchrB\t1\t1\n\nchrA\t0\t2\n").unwrap();
        let regions = bed(&p).unwrap();
        let a = &regions["chrA"];
        assert_eq!(
            a.iter()
                .map(|r| (r.start, r.stop, r.name.as_str()))
                .collect::<Vec<_>>(),
            [(2, 5, "one"), (2, 4, "two"), (7, 9, "last")]
        );
        assert_eq!(regions["chrB"][0].start, 1);
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn absent_and_filtered_contigs_and_overlapping_targets() {
        use nomics_io::bgzf::{BgzfReader, Stored};
        let dir =
            std::env::temp_dir().join(format!("nomics-mosdepth-small-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("small.bam");
        let refs = ["A", "absent", "filtered"]
            .map(|s| bam::Reference {
                name: s.as_bytes().to_vec(),
                len: 10,
            })
            .to_vec();
        let mut w = BgzfWriter::new(
            File::create(&bam).unwrap(),
            Stored,
            Profile {
                block_payload: 65280,
            },
        );
        bam::Header { text: vec![], refs }.write(&mut w).unwrap();
        let a = record(2, -1, vec![4 << 4], 0);
        bam::write_raw_record(&mut w, &a).unwrap();
        let mut f = record(2, -1, vec![4 << 4], 1024);
        f[..4].copy_from_slice(&2i32.to_le_bytes());
        bam::write_raw_record(&mut w, &f).unwrap();
        w.finish().unwrap();
        let targets = dir.join("targets.bed");
        fs::write(&targets,"A\t4\t8\tlater\nA\t1\t5\tearlier\nabsent\t1\t8\tmissing\nfiltered\t0\t11\tfiltered\nunknown\t0\t3\n").unwrap();
        let prefix = dir.join("out");
        run(&bam, &prefix, targets.to_str().unwrap(), false, false, 1).unwrap();
        let summary = fs::read_to_string(path(&prefix, ".mosdepth.summary.txt")).unwrap();
        assert_eq!(
            summary,
            "chrom\tlength\tbases\tmean\tmin\tmax\nA\t10\t4\t0.40\t0\t1\nA_region\t8\t5\t0.62\t0\t1\nfiltered\t10\t0\t0.00\t0\t0\nfiltered_region\t11\t0\t0.00\t0\t0\ntotal\t20\t4\t0.20\t0\t1\ntotal_region\t19\t5\t0.26\t0\t1\n"
        );
        let mut text = String::new();
        BgzfReader::new(File::open(path(&prefix, ".regions.bed.gz")).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(
            text,
            "A\t1\t5\tearlier\t0.75\nA\t4\t8\tlater\t0.50\nabsent\t1\t8\tmissing\t0.00\nfiltered\t0\t11\tfiltered\t0.00\n"
        );
        let dist = fs::read_to_string(path(&prefix, ".mosdepth.global.dist.txt")).unwrap();
        assert!(!dist.contains("absent"));
        assert!(dist.contains("filtered\t0\t1.00\n"));
        let mut text = String::new();
        BgzfReader::new(File::open(path(&prefix, ".per-base.bed.gz")).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(
            text,
            "A\t0\t2\t0\nA\t2\t6\t1\nA\t6\t10\t0\nabsent\t0\t10\t0\nfiltered\t0\t10\t0\n"
        );
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn long_cigar_placeholder_is_expanded_from_cg() {
        let r = bam::Record {
            ref_id: 0,
            pos: 1,
            mapq: 0,
            flag: 0,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
            read_name: b"long".to_vec(),
            cigar: vec![(8 << 4) | 4, (6 << 4) | 3],
            seq: vec![b'A'; 8],
            qual: vec![30; 8],
            tags: vec![(
                *b"CG",
                bam::Value::Array(bam::ArrayValue::Int(
                    bam::IntType::U32,
                    vec![2 << 4, (2 << 4) | 1, 4 << 4],
                )),
            )],
        };
        let mut raw = Vec::new();
        r.encode(bam::BinRule::Sam, &mut raw).unwrap();
        let cigar = decoded_cigar(RecordView::new(&raw).unwrap()).unwrap();
        assert_eq!(cigar, [2 << 4, (2 << 4) | 1, 4 << 4]);
        assert_eq!(segments(cigar, 1), [(1, 7)]);
    }
}
