//! Find each descriptor artifact on both sides, digest it under its
//! normalization, and on a mismatch run the detailed comparator for the first
//! difference. Fails closed: an artifact missing from either side, or empty
//! on either side, never passes.

use std::collections::BTreeMap;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use nomics_gate::ArtifactDigest;
use nomics_io::gzip::GzReader;

use super::descriptor::{ArtifactSpec, Comparator, Location, glob_match};

/// The directories a location's roots name.
#[derive(Debug, Default, Clone)]
pub struct Roots {
    /// sarek's `work/` (task directories two levels down).
    pub work: Option<PathBuf>,
    /// sarek's `--outdir`.
    pub results: Option<PathBuf>,
    /// The nomics run directory.
    pub run: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Pass,
    Fail,
    /// On one side only.
    Missing,
    /// Zero records on a side: an early harness once passed two empty outputs.
    Empty,
    /// A side's file could not be read (truncated, corrupt); the detail says which.
    Unreadable,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Fail => "FAIL",
            Status::Missing => "MISSING",
            Status::Empty => "EMPTY",
            Status::Unreadable => "UNREADABLE",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Difference {
    /// 1-based record (or line) number, when the comparator gives one.
    pub index: Option<u64>,
    pub detail: String,
}

#[derive(Debug, Clone)]
pub struct ArtifactResult {
    pub name: String,
    pub comparator: Comparator,
    pub normalization: String,
    pub reference: Option<ArtifactDigest>,
    pub nomics: Option<ArtifactDigest>,
    pub status: Status,
    pub first_difference: Option<Difference>,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Files under `dir` (recursively) whose name matches `pattern`.
fn walk(dir: &Path, pattern: &str, out: &mut Vec<PathBuf>) -> io::Result<()> {
    let mut entries: Vec<_> = std::fs::read_dir(dir)?.collect::<Result<_, _>>()?;
    entries.sort_by_key(|e| e.file_name());
    for e in entries {
        let path = e.path();
        let ft = e.file_type()?;
        if ft.is_dir() {
            walk(&path, pattern, out)?;
        } else if e.file_name().to_str().is_some_and(|n| glob_match(pattern, n)) {
            out.push(path);
        }
    }
    Ok(())
}

/// Files in `dir` itself whose name matches `pattern`.
fn list(dir: &Path, pattern: &str) -> io::Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    if !dir.is_dir() {
        return Ok(out);
    }
    for e in std::fs::read_dir(dir)? {
        let e = e?;
        if e.file_name().to_str().is_some_and(|n| glob_match(pattern, n)) && e.path().is_file() {
            out.push(e.path());
        }
    }
    out.sort();
    Ok(out)
}

fn search(roots: &Roots, root: &str, pattern: &str) -> io::Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    match root {
        "results" => {
            if let Some(d) = roots.results.as_deref().filter(|d| d.is_dir()) {
                walk(d, pattern, &mut out)?;
            }
        }
        // Nextflow task directories: work/ab/cdef.../
        "work" => {
            if let Some(d) = roots.work.as_deref().filter(|d| d.is_dir()) {
                let mut tops: Vec<_> = std::fs::read_dir(d)?.collect::<Result<_, _>>()?;
                tops.sort_by_key(|e| e.file_name());
                for top in tops.into_iter().filter(|e| e.path().is_dir()) {
                    let mut tasks: Vec<_> = std::fs::read_dir(top.path())?.collect::<Result<_, _>>()?;
                    tasks.sort_by_key(|e| e.file_name());
                    for task in tasks.into_iter().filter(|e| e.path().is_dir()) {
                        out.extend(list(&task.path(), pattern)?);
                    }
                }
            }
        }
        "run" => {
            if let Some(d) = roots.run.as_deref() {
                let (sub, file) = pattern.rsplit_once('/').unwrap_or(("", pattern));
                out.extend(list(&d.join(sub), file)?);
            }
        }
        other => return Err(invalid(format!("unknown location root {other}"))),
    }
    Ok(out)
}

/// The files a location names, by base name. The same file reached twice (a
/// task's staged input is a link to an earlier task's output) counts once;
/// two different files with one name is an error, never a silent pick.
pub fn locate(
    location: &Location,
    roots: &Roots,
    prefix: &str,
    lane: &str,
) -> io::Result<BTreeMap<String, PathBuf>> {
    for root in &location.roots {
        for pattern in location.patterns_for(prefix, lane) {
            let found = search(roots, root, &pattern)?;
            if found.is_empty() {
                continue;
            }
            let mut by_name: BTreeMap<String, (PathBuf, PathBuf)> = BTreeMap::new();
            for path in found {
                let name = path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .ok_or_else(|| invalid(format!("{}: unreadable name", path.display())))?
                    .to_owned();
                let real = std::fs::canonicalize(&path)?;
                match by_name.get(&name) {
                    Some((_, seen)) if *seen != real => {
                        return Err(invalid(format!(
                            "two different files named {name} ({} and {}); refusing to pick one",
                            seen.display(),
                            real.display()
                        )));
                    }
                    Some(_) => {}
                    None => {
                        by_name.insert(name, (path, real));
                    }
                }
            }
            return Ok(by_name.into_iter().map(|(n, (_, real))| (n, real)).collect());
        }
    }
    Ok(BTreeMap::new())
}

/// The files for one spec on one side, keyed by artifact name.
pub fn artifact_files(
    spec: &ArtifactSpec,
    location: &Location,
    roots: &Roots,
    prefix: &str,
    lane: &str,
) -> io::Result<BTreeMap<String, PathBuf>> {
    let found = locate(location, roots, prefix, lane)?;
    match spec.group() {
        Some(group) => Ok(found
            .into_iter()
            .map(|(name, path)| (format!("{group}/{name}"), path))
            .collect()),
        None => {
            if found.len() > 1 {
                return Err(invalid(format!(
                    "{}: {} files match ({}); expected one",
                    spec.name,
                    found.len(),
                    found.keys().cloned().collect::<Vec<_>>().join(", ")
                )));
            }
            Ok(found
                .into_values()
                .map(|path| (spec.name.clone(), path))
                .collect())
        }
    }
}

/// One prepass entry: the comparator, the file, and its digest once taken.
type DigestSlot = (Comparator, PathBuf, Option<io::Result<ArtifactDigest>>);

/// Digests computed ahead, keyed by comparator and path; [`Digested::take`] falls back to
/// computing one it does not hold.
pub struct Digested(std::sync::Mutex<Vec<DigestSlot>>);

impl Digested {
    /// Digest `jobs` on up to 12 threads (as `nomics dev digest`): each is a
    /// single decompress-and-hash stream, so a whole genome's BAMs, chunks and
    /// VCF otherwise run one after another (30 min on the 30x WGS).
    pub fn all(jobs: Vec<(Comparator, PathBuf)>) -> Digested {
        let mut jobs = jobs;
        jobs.sort_by(|a, b| a.1.cmp(&b.1));
        jobs.dedup_by(|a, b| a.1 == b.1 && a.0 == b.0);
        // Largest first, so the long BAM digests start at once.
        jobs.sort_by_key(|(_, p)| std::cmp::Reverse(std::fs::metadata(p).map_or(0, |m| m.len())));
        let workers = std::thread::available_parallelism().map_or(4, |n| n.get()).min(12);
        let queue = std::sync::Mutex::new(jobs.into_iter());
        let done = std::sync::Mutex::new(Vec::new());
        std::thread::scope(|scope| {
            for _ in 0..workers {
                scope.spawn(|| {
                    loop {
                        let Some((c, p)) = queue.lock().unwrap_or_else(|e| e.into_inner()).next() else {
                            break;
                        };
                        let d = digest(c, &p);
                        done.lock().unwrap_or_else(|e| e.into_inner()).push((c, p, Some(d)));
                    }
                });
            }
        });
        Digested(done)
    }

    pub fn take(&self, comparator: Comparator, path: &Path) -> io::Result<ArtifactDigest> {
        let held = self
            .0
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .iter_mut()
            .find(|(c, p, _)| *c == comparator && p == path)
            .and_then(|(_, _, d)| d.take());
        held.unwrap_or_else(|| digest(comparator, path))
    }
}

pub fn digest(comparator: Comparator, path: &Path) -> io::Result<ArtifactDigest> {
    let r = match comparator {
        Comparator::FastqChunk => nomics_gate::digest_fastq_chunk(path),
        // A whole-genome BAM (40-70 GB) is the gate's critical path; inflating
        // on 4 threads cut a 41 GB md.bam from 392 s to 159 s on the Spark
        // (8 gained nothing: the sequential hash is then the limit). The small
        // per-chunk BAMs already run side by side, one per worker.
        Comparator::Bam => {
            let big = std::fs::metadata(path).is_ok_and(|m| m.len() > 4 << 30);
            nomics_gate::digest_bam_threads(path, if big { 4 } else { 1 })
        }
        Comparator::MdMetrics => nomics_gate::digest_md_metrics(path),
        Comparator::Raw => nomics_gate::digest_raw(path),
        Comparator::Vcf => nomics_gate::digest_vcf(path),
    };
    r.map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
}

/// The first difference between two artifacts that digest differently.
pub fn first_difference(comparator: Comparator, nomics: &Path, reference: &Path) -> io::Result<Difference> {
    match comparator {
        Comparator::Bam => {
            let r = nomics_gate::compare_bams(nomics, reference, 1)?;
            if !r.header_diff.0.is_empty() || !r.header_diff.1.is_empty() || !r.refs_equal {
                let line = r
                    .header_diff
                    .1
                    .first()
                    .map(|l| format!("reference has {}", clip(l)))
                    .or_else(|| r.header_diff.0.first().map(|l| format!("nomics has {}", clip(l))))
                    .unwrap_or_else(|| "reference sequences differ".to_owned());
                return Ok(Difference {
                    index: None,
                    detail: format!("header differs: {line}"),
                });
            }
            if let Some((i, fields)) = r.examples.first() {
                return Ok(Difference {
                    index: Some(i + 1),
                    detail: fields.join("; "),
                });
            }
            Ok(Difference {
                index: None,
                detail: format!(
                    "record counts differ: reference {}, nomics {}",
                    r.records_golden, r.records_ours
                ),
            })
        }
        Comparator::Vcf => {
            let r = nomics_gate::compare_vcfs(nomics, reference, 1)?;
            if let Some((i, ours, golden)) = r.examples.first() {
                return Ok(Difference {
                    index: Some(*i as u64 + 1),
                    detail: vcf_difference(golden, ours),
                });
            }
            if !r.header_identical {
                return Ok(Difference {
                    index: None,
                    detail: "header differs (outside ##GATKCommandLine)".to_owned(),
                });
            }
            Ok(Difference {
                index: None,
                detail: format!(
                    "record counts differ: reference {}, nomics {}",
                    r.records_golden, r.records_ours
                ),
            })
        }
        Comparator::FastqChunk => {
            let read = |p: &Path| -> io::Result<Vec<u8>> {
                let mut v = Vec::new();
                GzReader::open(p)?.read_to_end(&mut v)?;
                Ok(v)
            };
            let (o, g) = (read(nomics)?, read(reference)?);
            let at = o
                .iter()
                .zip(&g)
                .position(|(a, b)| a != b)
                .unwrap_or(o.len().min(g.len()));
            let line = g[..at.min(g.len())].iter().filter(|&&c| c == b'\n').count() as u64;
            Ok(Difference {
                index: Some(line / 4 + 1),
                detail: format!(
                    "read {} differs (records: reference {}, nomics {})",
                    line / 4 + 1,
                    g.iter().filter(|&&c| c == b'\n').count() / 4,
                    o.iter().filter(|&&c| c == b'\n').count() / 4
                ),
            })
        }
        Comparator::MdMetrics | Comparator::Raw => {
            // md_metrics: physical lines 2 and 4 are not compared (GATES.md).
            let skip = |n: usize| comparator == Comparator::MdMetrics && (n == 2 || n == 4);
            let lines = |p: &Path| -> io::Result<Vec<String>> {
                BufReader::new(std::fs::File::open(p)?).lines().collect()
            };
            let (o, g) = (lines(nomics)?, lines(reference)?);
            for n in 1..=o.len().max(g.len()) {
                if skip(n) {
                    continue;
                }
                let (a, b) = (o.get(n - 1), g.get(n - 1));
                if a != b {
                    let show = |l: Option<&String>| l.map_or_else(|| "(no line)".to_owned(), |l| clip(l));
                    return Ok(Difference {
                        index: Some(n as u64),
                        detail: format!("line {n}: reference {} | nomics {}", show(b), show(a)),
                    });
                }
            }
            // Same lines, different digest: a final newline or line ending.
            Ok(Difference {
                index: None,
                detail: "same lines; bytes differ (line endings or final newline)".to_owned(),
            })
        }
    }
}

fn clip(s: &str) -> String {
    const MAX: usize = 120;
    if s.len() <= MAX {
        return s.to_owned();
    }
    let mut end = MAX;
    while !s.is_char_boundary(end) {
        end -= 1;
    }
    format!("{}…", &s[..end])
}

/// `chr20:31022441: QUAL 812.64 (reference) vs 812.65 (nomics)`: the site and
/// the columns that differ, not the whole record.
fn vcf_difference(reference: &str, nomics: &str) -> String {
    const COLUMNS: [&str; 9] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"];
    let (a, b): (Vec<&str>, Vec<&str>) = (reference.split('\t').collect(), nomics.split('\t').collect());
    let site = |c: &[&str]| format!("{}:{}", c.first().unwrap_or(&"?"), c.get(1).unwrap_or(&"?"));
    let mut parts = Vec::new();
    for i in 0..a.len().max(b.len()) {
        let (x, y) = (a.get(i).copied(), b.get(i).copied());
        if x != y {
            let name = COLUMNS.get(i).map_or_else(|| format!("sample {}", i - 8), |c| (*c).to_owned());
            parts.push(format!(
                "{name} {} (reference) vs {} (nomics)",
                clip(x.unwrap_or("(none)")),
                clip(y.unwrap_or("(none)"))
            ));
        }
    }
    if site(&a) != site(&b) {
        return format!("site {} (reference) vs {} (nomics)", site(&a), site(&b));
    }
    format!("{}: {}", site(&a), parts.join("; "))
}

/// Gate every descriptor artifact. `Err` only for a harness problem (an
/// unreadable file, an ambiguous match); differences are results.
pub fn gate_all(
    specs: &[ArtifactSpec],
    reference: &Roots,
    nomics: &Roots,
    prefix: &str,
    lane: &str,
) -> io::Result<Vec<ArtifactResult>> {
    let mut jobs = Vec::new();
    for spec in specs {
        for (loc, roots) in [(&spec.reference, reference), (&spec.nomics, nomics)] {
            for p in artifact_files(spec, loc, roots, prefix, lane)?.into_values() {
                jobs.push((spec.comparator, p));
            }
        }
    }
    let digested = Digested::all(jobs);
    let mut out = Vec::new();
    for spec in specs {
        let refs = artifact_files(spec, &spec.reference, reference, prefix, lane)?;
        let ours = artifact_files(spec, &spec.nomics, nomics, prefix, lane)?;
        let mut names: Vec<&String> = refs.keys().chain(ours.keys()).collect();
        names.sort();
        names.dedup();
        if names.is_empty() {
            // Nothing on either side: still a result, so the report shows it.
            out.push(ArtifactResult {
                name: spec.name.clone(),
                comparator: spec.comparator,
                normalization: spec.normalization.clone(),
                reference: None,
                nomics: None,
                status: Status::Missing,
                first_difference: None,
            });
            continue;
        }
        for name in names {
            let (r, o) = (refs.get(name), ours.get(name));
            // A file that cannot be digested is a result for that artifact,
            // not an abort: a truncated nomics output must show up as such.
            let rd = r.map(|p| digested.take(spec.comparator, p));
            let od = o.map(|p| digested.take(spec.comparator, p));
            let unreadable = [("reference", &rd), ("nomics", &od)]
                .into_iter()
                .find_map(|(side, d)| match d {
                    Some(Err(e)) => Some(format!("{side} file unreadable: {e}")),
                    _ => None,
                });
            let (rd, od) = (rd.and_then(Result::ok), od.and_then(Result::ok));
            let (status, first_difference) = if let Some(why) = unreadable {
                (Status::Unreadable, Some(Difference { index: None, detail: why }))
            } else {
                match (&rd, &od) {
                    (Some(a), Some(b)) if a.count == 0 || b.count == 0 => (Status::Empty, None),
                    (Some(a), Some(b)) if a == b => (Status::Pass, None),
                    (Some(_), Some(_)) => (
                        Status::Fail,
                        Some(first_difference(spec.comparator, o.unwrap(), r.unwrap())?),
                    ),
                    _ => (Status::Missing, None),
                }
            };
            out.push(ArtifactResult {
                name: name.clone(),
                comparator: spec.comparator,
                normalization: spec.normalization.clone(),
                reference: rd,
                nomics: od,
                status,
                first_difference,
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::descriptor::{Descriptor, SAREK_GERMLINE};

    #[test]
    fn vcf_differences_name_the_site_and_columns() {
        let r = "chr20\t100\t.\tA\tG\t812.64\t.\tDP=5\tGT:PL\t0/1:1,2,3";
        let n = "chr20\t100\t.\tA\tG\t812.65\t.\tDP=5\tGT:PL\t0/1:1,2,4";
        assert_eq!(
            vcf_difference(r, n),
            "chr20:100: QUAL 812.64 (reference) vs 812.65 (nomics); sample 1 0/1:1,2,3 (reference) vs 0/1:1,2,4 (nomics)"
        );
        let moved = n.replace("\t100\t", "\t101\t");
        assert_eq!(vcf_difference(r, &moved), "site chr20:100 (reference) vs chr20:101 (nomics)");
    }

    #[test]
    fn locates_sarek_and_nomics_layouts() {
        let base = std::env::temp_dir().join(format!("nomics-locate-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        let work = base.join("work");
        let results = base.join("results");
        let run = base.join("run");
        let task_a = work.join("ab/cdef");
        let task_b = work.join("12/3456");
        for d in [&task_a, &task_b, &results.join("preprocessing/markduplicates/S"), &run.join("prep")] {
            std::fs::create_dir_all(d).unwrap();
        }
        std::fs::write(task_a.join("0001.S-L1_R1.fastp.fastq.gz"), b"x").unwrap();
        std::fs::write(task_a.join(".command.sh"), b"x").unwrap();
        // A later task's staged input links back to the first task's output.
        std::os::unix::fs::symlink(task_a.join("0001.S-L1_R1.fastp.fastq.gz"), task_b.join("0001.S-L1_R1.fastp.fastq.gz")).unwrap();
        std::fs::write(results.join("preprocessing/markduplicates/S/S.md.bam"), b"x").unwrap();
        std::fs::write(run.join("prep/0001.S-L1_R1.fastp.fastq.gz"), b"x").unwrap();
        std::fs::write(run.join("S.md.bam"), b"x").unwrap();

        let d = Descriptor::parse(SAREK_GERMLINE).unwrap();
        let roots = Roots { work: Some(work.clone()), results: Some(results.clone()), run: Some(run.clone()) };
        let spec = |n: &str| d.artifacts.iter().find(|a| a.name == n).unwrap().clone();
        let chunk = spec("chunk/*");
        let refs = artifact_files(&chunk, &chunk.reference, &roots, "S", "L1").unwrap();
        assert_eq!(refs.keys().collect::<Vec<_>>(), ["chunk/0001.S-L1_R1.fastp.fastq.gz"]);
        let ours = artifact_files(&chunk, &chunk.nomics, &roots, "S", "L1").unwrap();
        assert_eq!(ours.keys().collect::<Vec<_>>(), ["chunk/0001.S-L1_R1.fastp.fastq.gz"]);
        let md = spec("md_bam");
        assert_eq!(artifact_files(&md, &md.reference, &roots, "S", "L1").unwrap().len(), 1);
        assert_eq!(artifact_files(&md, &md.nomics, &roots, "S", "L1").unwrap().len(), 1);

        // Two different files with one name is an error, not a pick.
        std::fs::remove_file(task_b.join("0001.S-L1_R1.fastp.fastq.gz")).unwrap();
        std::fs::write(task_b.join("0001.S-L1_R1.fastp.fastq.gz"), b"y").unwrap();
        assert!(artifact_files(&chunk, &chunk.reference, &roots, "S", "L1").is_err());
        std::fs::remove_dir_all(&base).unwrap();
    }
}
