//! First Mutect2 stage, default sarek tumour-normal mode, GATK 4.6.2.0.
//! No assembly or variant calling. `run_with_sink` is the next-stage boundary:
//! reads are shared across regions and must be consumed at the callback, in order.
pub(crate) mod activity;
pub mod numerics;
mod reads;

use nomics_hc::{
    activity::{BinaryActivityProfile, Region},
    finalize::HcRead,
    jrandom::JavaRandom,
    reads::Intervals,
};
use nomics_io::{bam_query::BamQuery, fasta::read_fai, tabix::Tabix, vcf::parse_line};
use nomics_mem::HostBuf;
use std::{
    collections::{HashMap, VecDeque},
    fs::File,
    io::{self, BufWriter, Read as _, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

/// Read after palindrome clipping and downsampling. Indel qualities are carried
/// separately because Mutect2's overlap adjustment mutates BI and BD as well.
#[derive(Clone, Debug)]
pub struct Read {
    pub alignment: HcRead,
    /// Input file ordinal and BGZF virtual offset just after this record.
    pub source: (usize, u64),
    pub sample: String,
    pub is_normal: bool,
    /// Whether overlap correction materialized BI/BD tags (even if later clipped).
    pub overlap_indels: bool,
    pub insertion_quals: Vec<u8>,
    pub deletion_quals: Vec<u8>,
}
#[derive(Default, Debug)]
pub struct Stats {
    pub processed: u64,
    pub filtered: [u64; 12],
    pub callable: u64,
    pub regions: usize,
    pub active: usize,
}
impl Stats {
    pub fn filter_summary(&self) -> String {
        let mut out = String::new();
        for (name, count) in reads::FILTERS.iter().zip(self.filtered) {
            out.push_str(&format!("{count} read(s) filtered by: {name} \n"));
        }
        out.push_str(&format!(
            "{} total reads filtered out of {} reads processed\n",
            self.filtered.iter().sum::<u64>(),
            self.processed
        ));
        out
    }
}
pub struct Inputs<'a> {
    pub input: &'a [PathBuf],
    pub normal_sample: &'a str,
    pub reference: &'a Path,
    pub intervals: &'a Path,
    pub panel_of_normals: &'a Path,
    /// Accepted for the sarek contract; tumour-normal activity never queries it.
    pub germline_resource: &'a Path,
    pub regions: &'a Path,
    pub stats: &'a Path,
    /// Optional exact internal-state dumps, matching the Java test harness.
    pub dump_prefix: Option<&'a Path>,
}
/// Called once per region, active or inactive. References include callRegion's
/// overlap edits, before any downstream preparation copies. Do not retain or
/// independently replay regions: the arena carries edits into subsequent calls.
pub trait RegionSink {
    fn on_region(&mut self, region: &Region, reads: &[&Read]) -> io::Result<()>;
}
impl<F> RegionSink for F
where
    F: FnMut(&Region, &[&Read]) -> io::Result<()>,
{
    fn on_region(&mut self, region: &Region, reads: &[&Read]) -> io::Result<()> {
        self(region, reads)
    }
}
/// Activity-pileup hook, called in traversal order before LOD/normal/PoN rejection.
/// The arena includes prior callRegion overlap edits. `active` is the downsampled
/// locus candidate list; use `activity::covered` to apply LIBS adapter/CIGAR rules.
/// Only nonempty joint pileups are delivered. Do not retain read references.
pub trait LocusSink {
    fn on_samples(&mut self, _tumour_samples: &[String]) -> io::Result<()> {
        Ok(())
    }
    fn on_locus(
        &mut self,
        position: i64,
        reference: &[u8],
        reads: &[Read],
        active: &[usize],
    ) -> io::Result<()>;
}
/// The locus sink of a run that collects nothing per locus.
pub(crate) struct NoLoci;
impl LocusSink for NoLoci {
    fn on_locus(&mut self, _: i64, _: &[u8], _: &[Read], _: &[usize]) -> io::Result<()> {
        Ok(())
    }
}
pub fn run(inputs: &Inputs<'_>) -> io::Result<Stats> {
    run_with_sink(inputs, &mut |_: &Region, _: &[&Read]| Ok(()))
}

fn reference(
    path: &Path,
    names: &[(Vec<u8>, i64)],
    ivs: &Intervals,
) -> io::Result<Vec<Option<HostBuf>>> {
    let index = read_fai(Path::new(&format!("{}.fai", path.display())))?;
    let mut file = File::open(path)?;
    let mut out: Vec<Option<HostBuf>> = (0..names.len()).map(|_| None).collect();
    for (i, (name, len)) in names
        .iter()
        .enumerate()
        .filter(|(i, _)| !ivs.unpadded[*i].is_empty())
    {
        let e = index
            .iter()
            .find(|e| e.name.as_bytes() == name)
            .ok_or_else(|| io::Error::other("missing reference contig"))?;
        if e.len != *len as u64 {
            return Err(io::Error::other("reference length mismatch"));
        }
        file.seek(SeekFrom::Start(e.offset))?;
        let n = e.len.div_ceil(e.line_bases) * e.line_width;
        let mut raw = HostBuf::new(n as usize)?;
        let mut count = 0;
        while count < raw.as_slice().len() {
            let n = file.read(&mut raw.as_mut_slice()[count..])?;
            if n == 0 {
                break;
            }
            count += n;
        }
        let mut seq = HostBuf::new(e.len as usize)?;
        let mut at = 0;
        for b in raw.as_slice()[..count]
            .iter()
            .copied()
            .filter(|b| !b.is_ascii_whitespace())
            .take(e.len as usize)
        {
            seq.as_mut_slice()[at] = match b.to_ascii_uppercase() {
                x @ (b'A' | b'C' | b'G' | b'T') => x,
                _ => b'N',
            };
            at += 1;
        }
        if at != e.len as usize {
            return Err(io::Error::other("short FASTA"));
        }
        out[i] = Some(seq);
    }
    Ok(out)
}
fn downsample(all: &[Read]) -> Vec<usize> {
    let mut rng = JavaRandom::gatk();
    let mut kept = Vec::new();
    let mut at = 0;
    while at < all.len() {
        let a = &all[at].alignment;
        let end = at
            + all[at..].partition_point(|r| {
                (r.alignment.contig, r.alignment.start) == (a.contig, a.start)
            });
        if end - at <= 50 {
            kept.extend(at..end);
        } else {
            let mut reservoir = Vec::new();
            let mut seen = 0;
            for (id, _) in all
                .iter()
                .enumerate()
                .take(end)
                .skip(at)
                .filter(|(_, r)| r.alignment.mapq > 50)
            {
                seen += 1;
                if reservoir.len() < 50 {
                    reservoir.push(id);
                } else {
                    let slot = rng.next_int(seen) as usize;
                    if slot < 50 {
                        reservoir[slot] = id;
                    }
                }
            }
            kept.extend(reservoir);
        }
        at = end;
    }
    kept
}
fn soft_start(r: &HcRead) -> i64 {
    let soft = r
        .cigar
        .iter()
        .take_while(|c| matches!(**c & 15, 4 | 5))
        .filter(|c| **c & 15 == 4)
        .map(|c| i64::from(c >> 4))
        .sum::<i64>();
    r.start - soft
}
/// ReadUtils.getReadIndexForReferenceCoordinate: soft clips consume reference
/// here, deletions return the *following* read offset (preserve upstream bug).
fn offset(r: &HcRead, pos: i64) -> Option<(usize, u32)> {
    let (mut p, mut off) = (soft_start(r), 0usize);
    for &c in &r.cigar {
        let (op, n) = (c & 15, (c >> 4) as usize);
        let consumes_read = matches!(op, 0 | 1 | 4 | 7 | 8);
        if matches!(op, 0 | 2 | 3 | 4 | 7 | 8) {
            if p <= pos && pos < p + n as i64 {
                return Some((off + if consumes_read { (pos - p) as usize } else { 0 }, op));
            }
            p += n as i64;
        }
        if consumes_read {
            off += n;
        }
    }
    None
}
fn clean_pairs(all: &mut [Read], ids: &[usize]) {
    let mut pending = HashMap::new();
    for &id in ids {
        let r = &all[id];
        let a = &r.alignment;
        if a.flag & 1 == 0 || a.flag & 8 != 0 || a.mate_start == 0 || a.mate_start > a.end() {
            continue;
        }
        let key = (r.is_normal, a.name.clone());
        if let Some(first) = pending.remove(&key) {
            let (left, right) = if first < id {
                let (l, r) = all.split_at_mut(id);
                (&mut l[first], &mut r[0])
            } else {
                let (l, r) = all.split_at_mut(first);
                (&mut r[0], &mut l[id])
            };
            // FragmentUtils orders by *soft* start and reverses on ties.
            let (left, right) = if soft_start(&left.alignment) < soft_start(&right.alignment) {
                (left, right)
            } else {
                (right, left)
            };
            let (a, b) = (&left.alignment, &right.alignment);
            if a.contig != b.contig || a.end() < b.start {
                continue;
            }
            let Some((ao, op)) = offset(a, b.start) else {
                continue;
            };
            if matches!(op, 4 | 5) {
                continue;
            }
            let Some((ae, _)) = offset(a, a.end()) else {
                continue;
            };
            let Some((bo, _)) = offset(b, b.start) else {
                continue;
            };
            let Some((be, _)) = offset(b, b.end()) else {
                continue;
            };
            let n = (ae as i64 - ao as i64).min(be as i64 - bo as i64) + 1;
            left.overlap_indels = true;
            right.overlap_indels = true;
            for k in 0..n.max(0) as usize {
                let (i, j) = (ao + k, bo + k);
                if left.alignment.bases[i] == right.alignment.bases[j] {
                    left.alignment.quals[i] = (left.alignment.quals[i] as i8).min(20) as u8;
                    right.alignment.quals[j] = (right.alignment.quals[j] as i8).min(20) as u8;
                }
                left.insertion_quals[i] = (left.insertion_quals[i] as i8).min(20) as u8;
                left.deletion_quals[i] = (left.deletion_quals[i] as i8).min(20) as u8;
                right.insertion_quals[j] = (right.insertion_quals[j] as i8).min(20) as u8;
                right.deletion_quals[j] = (right.deletion_quals[j] as i8).min(20) as u8;
            }
        } else {
            pending.insert(key, id);
        }
    }
}
fn call(
    region: Region,
    ids: &[usize],
    all: &mut [Read],
    sink: &mut dyn RegionSink,
    prepared: &mut Option<BufWriter<File>>,
) -> io::Result<()> {
    clean_pairs(all, ids);
    if let Some(w) = prepared {
        for &id in ids {
            let r = &all[id];
            w.write_all(&(r.alignment.bases.len() as u32).to_be_bytes())?;
            for bytes in [
                &r.alignment.bases,
                &r.alignment.quals,
                &r.insertion_quals,
                &r.deletion_quals,
            ] {
                w.write_all(bytes)?;
            }
        }
    }
    sink.on_region(&region, &ids.iter().map(|&i| &all[i]).collect::<Vec<_>>())
}
fn cigar_string(r: &HcRead) -> String {
    r.cigar
        .iter()
        .map(|c| {
            format!(
                "{}{}",
                c >> 4,
                nomics_io::bam::CIGAR_OPS[(c & 15) as usize] as char
            )
        })
        .collect()
}

pub fn run_with_sink(inputs: &Inputs<'_>, sink: &mut dyn RegionSink) -> io::Result<Stats> {
    run_with_sinks(inputs, sink, &mut NoLoci)
}
/// Run region calling and activity collection in a single BAM traversal.
pub fn run_with_sinks(
    inputs: &Inputs<'_>,
    sink: &mut dyn RegionSink,
    loci: &mut dyn LocusSink,
) -> io::Result<Stats> {
    if inputs.input.len() != 2 {
        return Err(io::Error::other(
            "only two BAM tumour-normal mode is supported",
        ));
    }
    File::open(inputs.germline_resource)?;
    let header = BamQuery::open(&inputs.input[0], &[])?.header;
    let names: Vec<_> = header
        .refs
        .iter()
        .map(|r| (r.name.clone(), i64::from(r.len)))
        .collect();
    let ivs = Intervals::from_bed(inputs.intervals, &names)?;
    let reference = reference(inputs.reference, &names, &ivs)?;
    let mut stats = Stats::default();
    let mut all = reads::load(
        inputs.input,
        &ivs,
        &reference,
        inputs.normal_sample,
        &mut stats,
    )?;
    let mut tumour_samples = std::collections::BTreeSet::new();
    for path in inputs.input {
        let header = BamQuery::open(path, &[])?.header;
        for sample in reads::samples(&header)?.into_values() {
            if sample != inputs.normal_sample {
                tumour_samples.insert(sample);
            }
        }
    }
    loci.on_samples(&tumour_samples.into_iter().collect::<Vec<_>>())?;
    let kept = downsample(&all);
    let pon = Tabix::open(inputs.panel_of_normals)?;
    let mut out = BufWriter::new(File::create(inputs.regions)?);
    writeln!(
        out,
        "#track graphType=line\nChromosome\tStart\tEnd\tFeature\tAssemblyRegions"
    )?;
    let mut activities = inputs
        .dump_prefix
        .map(|p| File::create(format!("{}.activity.tsv", p.display())).map(BufWriter::new))
        .transpose()?;
    let mut members = inputs
        .dump_prefix
        .map(|p| File::create(format!("{}.reads.tsv", p.display())).map(BufWriter::new))
        .transpose()?;
    let mut prepared = inputs
        .dump_prefix
        .map(|p| File::create(format!("{}.prepared.bin", p.display())).map(BufWriter::new))
        .transpose()?;
    for (c, intervals) in ivs
        .unpadded
        .iter()
        .enumerate()
        .filter(|(_, i)| !i.is_empty())
    {
        let contig = String::from_utf8_lossy(&names[c].0);
        let ids: Vec<_> = kept
            .iter()
            .copied()
            .filter(|&i| all[i].alignment.contig == c as i32)
            .collect();
        let mut pon_spans = Vec::new();
        for iv in intervals {
            for line in pon.vcf_lines_overlapping(&names[c].0, iv.start as u64, iv.end as u64)? {
                let v = parse_line(&line)?;
                let mut end = v.pos as i64 + v.ref_allele.len() as i64 - 1;
                for field in v.info.split(|&b| b == b';') {
                    if let Some(e) = field.strip_prefix(b"END=") {
                        end = std::str::from_utf8(e)
                            .map_err(io::Error::other)?
                            .parse()
                            .map_err(io::Error::other)?;
                    }
                }
                pon_spans.push((v.pos as i64, end));
            }
        }
        pon_spans.sort_unstable();
        pon_spans.dedup();
        let mut profile = BinaryActivityProfile::new(names[c].1);
        let mut pending = VecDeque::new();
        let mut previous = Vec::new();
        let mut cached_at = 0;
        let mut ready: Option<(Region, Vec<usize>)> = None;
        let mut admitted = 0;
        let mut active: Vec<usize> = Vec::new();
        let (mut tumour, mut normal) =
            (activity::Evidence::default(), activity::Evidence::default());
        let mut emit = |r: Region, all: &mut [Read]| -> io::Result<()> {
            let (ps, pe) = ((r.start - 100).max(1), (r.end + 100).min(names[c].1));
            let mut m: Vec<usize> = previous
                .iter()
                .copied()
                .filter(|&id: &usize| {
                    all[id].alignment.start <= pe && all[id].alignment.end() >= ps
                })
                .collect();
            while cached_at < ids.len() && all[ids[cached_at]].alignment.start <= pe {
                let id = ids[cached_at];
                if all[id].alignment.end() >= ps {
                    m.push(id);
                }
                cached_at += 1;
            }
            previous = m.clone();
            stats.regions += 1;
            stats.active += usize::from(r.active);
            writeln!(
                out,
                "{contig}\t{}\t{}\tend-marker\t0.00000",
                r.start - 1,
                r.start
            )?;
            writeln!(
                out,
                "{contig}\t{}\t{}\tsize={}\t{}",
                r.start - 1,
                r.end,
                r.end - r.start + 1,
                if r.active { "1.00000" } else { "-1.00000" }
            )?;
            if let Some(w) = &mut members {
                for &id in &m {
                    let a = &all[id].alignment;
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        stats.regions,
                        all[id].sample,
                        String::from_utf8_lossy(&a.name),
                        a.flag,
                        a.start,
                        cigar_string(a)
                    )?;
                }
            }
            if let Some((pr, pm)) = ready.replace((r, m)) {
                call(pr, &pm, all, sink, &mut prepared)?;
            }
            Ok(())
        };
        for iv in intervals {
            for pos in iv.start..=iv.end {
                while admitted < ids.len() && all[ids[admitted]].alignment.start <= pos {
                    active.push(ids[admitted]);
                    admitted += 1;
                }
                active.retain(|&id| all[id].alignment.end() >= pos);
                tumour.clear();
                normal.clear();
                for &id in &active {
                    let r = &all[id];
                    if r.is_normal {
                        normal.add(
                            r,
                            pos,
                            reference[c].as_ref().unwrap().as_slice()[pos as usize - 1],
                        );
                    } else {
                        tumour.add(
                            r,
                            pos,
                            reference[c].as_ref().unwrap().as_slice()[pos as usize - 1],
                        );
                    }
                }
                if tumour.depth + normal.depth > 0 {
                    loci.on_locus(
                        pos,
                        reference[c].as_ref().unwrap().as_slice(),
                        &all,
                        &active,
                    )?;
                }
                let in_pon = pon_spans
                    .iter()
                    .take_while(|s| s.0 <= pos)
                    .any(|s| s.1 >= pos);
                let probability = f64::from(activity::is_active(
                    &tumour,
                    &normal,
                    in_pon,
                    &mut stats.callable,
                ));
                if let Some(w) = &mut activities {
                    writeln!(w, "{contig}\t{pos}\t{}", probability as u8)?;
                }
                pending.extend(profile.push(c, pos, probability));
                if let Some(r) =
                    pending.pop_front_if(|r: &mut Region| pos > (r.end + 100).min(names[c].1))
                {
                    emit(r, &mut all)?;
                }
            }
        }
        pending.extend(profile.finish(c));
        while let Some(r) = pending.pop_front() {
            emit(r, &mut all)?;
        }
        if let Some((r, m)) = ready.take() {
            call(r, &m, &mut all, sink, &mut prepared)?;
        }
    }
    if let Some(w) = &mut prepared {
        w.flush()?;
    }
    out.flush()?;
    if let Some(w) = &mut activities {
        w.flush()?;
    }
    if let Some(w) = &mut members {
        w.flush()?;
    }
    std::fs::write(
        inputs.stats,
        format!(
            "statistic\tvalue\ncallable\t{}\n",
            crate::math::javafmt::double_to_string(stats.callable as f64)
        ),
    )?;
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn joint_sample_downsampling_matches_unmodified_jar() {
        let all: Vec<_> = (0..143)
            .map(|id| {
                let start = if id < 80 {
                    10
                } else if id < 140 {
                    20
                } else {
                    30
                };
                let mapq = if id < 80 {
                    if id % 7 == 0 { 50 } else { 60 }
                } else if id < 140 {
                    if id < 90 { 60 } else { 50 }
                } else {
                    20
                };
                Read {
                    source: (0, id as u64),
                    alignment: HcRead {
                        name: id.to_string().into_bytes(),
                        flag: 0,
                        contig: 0,
                        start,
                        cigar: vec![16],
                        bases: vec![b'A'],
                        quals: vec![40],
                        mapq,
                        mate_contig: -1,
                        mate_start: 0,
                        tlen: 0,
                        record: id as u32,
                    },
                    sample: if id % 2 == 0 { "normal" } else { "tumour" }.into(),
                    is_normal: id % 2 == 0,
                    overlap_indels: false,
                    insertion_quals: vec![45],
                    deletion_quals: vec![45],
                }
            })
            .collect();
        let expected: Vec<usize> = include_str!("numerics_oracle.tsv")
            .lines()
            .find(|l| l.starts_with("DS\t"))
            .unwrap()
            .split('\t')
            .skip(1)
            .map(|s| s.parse().unwrap())
            .collect();
        assert_eq!(downsample(&all), expected);
        assert_eq!(expected.len(), 63); // 50 from joint pool, ten retained with no low-MQ refill, then all three.
    }
}
