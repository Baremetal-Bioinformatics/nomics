//! GATK 4.6.2.0 walkers/contamination/{PileupSummary,GatherPileupSummaries}.java.
use crate::math::javafmt::double_to_string;
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};
pub const HEADER: &str =
    "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency";
#[derive(Clone, Debug)]
pub struct PileupSummary {
    pub contig: String,
    pub position: i32,
    pub ref_count: i32,
    pub alt_count: i32,
    pub other_alt_count: i32,
    pub allele_frequency: f64,
}
impl PileupSummary {
    pub fn total(&self) -> i32 {
        self.ref_count
            .wrapping_add(self.alt_count)
            .wrapping_add(self.other_alt_count)
    }
    pub fn alt_fraction(&self) -> f64 {
        if self.total() == 0 {
            0.0
        } else {
            f64::from(self.alt_count) / f64::from(self.total())
        }
    }
}
pub fn read(path: &Path) -> io::Result<(String, Vec<PileupSummary>)> {
    let mut sample = None;
    let mut columns = None;
    let mut rows = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if let Some(s) = line.strip_prefix("#<METADATA>SAMPLE=") {
            sample = Some(s.to_owned());
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let values: Vec<_> = line.split('\t').collect();
        if columns.is_none() {
            let c: io::Result<Vec<_>> = HEADER
                .split('\t')
                .map(|name| {
                    values
                        .iter()
                        .position(|v| *v == name)
                        .ok_or_else(|| io::Error::other(format!("missing column {name}")))
                })
                .collect();
            columns = Some(c?);
            continue;
        }
        let cols = columns.as_ref().unwrap();
        let field = |i: usize| {
            values
                .get(cols[i])
                .copied()
                .ok_or_else(|| io::Error::other("short pileup row"))
        };
        let integer = |i| field(i)?.parse::<i32>().map_err(io::Error::other);
        rows.push(PileupSummary {
            contig: field(0)?.to_owned(),
            position: integer(1)?,
            ref_count: integer(2)?,
            alt_count: integer(3)?,
            other_alt_count: integer(4)?,
            allele_frequency: field(5)?.parse().map_err(io::Error::other)?,
        });
    }
    if columns.is_none() {
        return Err(io::Error::other("missing pileup header"));
    }
    Ok((
        sample.ok_or_else(|| io::Error::other("missing SAMPLE metadata"))?,
        rows,
    ))
}
pub fn write(path: &Path, sample: Option<&str>, rows: &[PileupSummary]) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write_header(&mut w, sample)?;
    for r in rows {
        write_record(&mut w, r)?;
    }
    w.flush()
}
pub(crate) fn write_header(w: &mut impl Write, sample: Option<&str>) -> io::Result<()> {
    if let Some(s) = sample {
        writeln!(w, "#<METADATA>SAMPLE={s}")?;
    }
    writeln!(w, "{HEADER}")
}
pub(crate) fn write_record(w: &mut impl Write, r: &PileupSummary) -> io::Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}",
        r.contig,
        r.position,
        r.ref_count,
        r.alt_count,
        r.other_alt_count,
        double_to_string(r.allele_frequency)
    )
}
pub fn gather(inputs: &[PathBuf], dictionary: &Path, output: &Path) -> io::Result<()> {
    let dict = std::fs::read_to_string(dictionary)?;
    let contigs: Vec<_> = dict
        .lines()
        .filter(|s| s.starts_with("@SQ\t"))
        .filter_map(|s| s.split('\t').find_map(|v| v.strip_prefix("SN:")))
        .collect();
    if contigs.is_empty() {
        return Err(io::Error::other(
            "sequence dictionary has no @SQ SN entries",
        ));
    }
    let mut tables = Vec::new();
    for input in inputs {
        let t = read(input)?;
        if !t.1.is_empty() {
            tables.push(t);
        }
    }
    let index = |name: &str| {
        contigs
            .iter()
            .position(|s| *s == name)
            .map_or(-1, |i| i as i32)
    };
    tables.sort_by_key(|(_, r)| (index(&r[0].contig), r[0].position)); // stable; unknown contigs have Java indexOf=-1
    let sample = tables.first().map(|t| t.0.clone());
    let mut rows = Vec::new();
    for (s, mut r) in tables {
        if Some(&s) != sample.as_ref() {
            return Err(io::Error::other(
                "Combining PileupSummaryTables from different samples is not supported",
            ));
        }
        rows.append(&mut r);
    }
    write(output, sample.as_deref(), &rows)
}
