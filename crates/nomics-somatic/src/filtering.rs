//! Fresh GATK 4.6.2.0 FilterMutectCalls tumour-normal default traversal.
#![allow(clippy::needless_range_loop)]
//! Reference source: mutect/filtering and mutect/clustering, tag 4.6.2.0.
mod cluster;
mod numerics;
mod record;
mod strand;
use crate::{
    math::{
        binomial::{binomial_probability, log_binomial_probability},
        javafmt::double_to_string,
    },
    orientation::{artifact_probability, jmath, tar},
};
use cluster::Model;
use numerics::{log_sum, log1mexp, normalize, posterior, regularized_beta};
use record::{Record, floats, ints};
use std::{
    collections::{BTreeMap, HashMap},
    fs::{self, File},
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};
use strand::{Step, Strand};
fn invalid(s: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s.into())
}
const NAMES: [&str; 19] = [
    "weak_evidence",
    "base_qual",
    "map_qual",
    "duplicate",
    "strand_bias",
    "contamination",
    "strict_strand",
    "position",
    "low_allele_frac",
    "normal_artifact",
    "n_ratio",
    "panel_of_normals",
    "orientation",
    "clustered_events",
    "multiallelic",
    "fragment",
    "slippage",
    "haplotype",
    "germline",
];
fn allele_filter(i: usize) -> bool {
    i < 9
}
fn non_somatic(i: usize) -> bool {
    i == 5 || i == 18
}
struct Probabilities {
    p: Vec<Vec<f64>>,
    artifact: Vec<f64>,
    non_somatic: Vec<f64>,
    combined: Vec<f64>,
    steps: Vec<Step>,
}
struct Engine {
    model: Model,
    strand: Strand,
    contamination: f64,
    tumor: usize,
    phased: HashMap<String, Vec<(u64, f64)>>,
    error_order: [usize; 3],
}
impl Engine {
    fn probabilities(&self, r: &Record) -> Probabilities {
        let n = r.alts.len();
        let mut p = vec![vec![0.0; n]; 19];
        let mut steps = Vec::new();
        for i in 0..n {
            p[0][i] = self.model.sequencing(&r.datum(i, 0.0, 0.0));
        }
        for (filter, key, offset, limit) in [
            (1, "MBQ", 1, 20.0),
            (2, "MMQ", 1, 30.0),
            (3, "AS_UNIQ_ALT_READ_COUNT", 0, 0.0),
            (7, "MPOS", 0, 1.0),
        ] {
            let v = r.values(key);
            if v.is_empty() {
                p[filter].clear();
                continue;
            }
            for i in 0..n {
                let j = if filter == 2 && r.len(i) >= 5 {
                    0
                } else {
                    i + offset
                };
                let x = v[j];
                p[filter][i] = f64::from(if filter == 3 {
                    x <= limit
                } else if filter == 7 {
                    x > -1.0 && x < limit
                } else {
                    x < limit
                });
            }
        }
        if let Some(sb) = r.info.get("AS_SB_TABLE") {
            let counts: Vec<_> = sb.split('|').map(ints).collect();
            let nf = counts.iter().map(|v| v[0]).sum();
            let nr = counts.iter().map(|v| v[1]).sum();
            for i in 0..n {
                let s = self.strand.step(
                    nf,
                    nr,
                    counts[i + 1][0],
                    counts[i + 1][1],
                    r.len(i).unsigned_abs(),
                );
                p[4][i] = s.p();
                steps.push(s);
            }
        } else {
            p[4].clear();
        }
        p[6].clear();
        if let Some(af) = r.gt[self.tumor].get("AF") {
            for (i, x) in floats(af).iter().enumerate() {
                p[8][i] = f64::from(*x < 0.0);
            }
        }

        let pop = r.values("POPAF");
        if pop.is_empty() {
            p[5].clear();
        } else {
            for i in 0..n {
                let f = jmath::pow(10.0, -pop[i]);
                let c = self.contamination.clamp(0.0, 1.0 - 1e-10);
                let k = r.ad[i + 1];
                let d = r.depth();
                let single =
                    2.0 * f * (1.0 - f) * binomial_probability(d as i32, k as i32, c / 2.0)
                        + f * f * binomial_probability(d as i32, k as i32, c);
                let many = binomial_probability(d as i32, k as i32, c * f);
                let ll = self.model.likelihood(d, k);
                p[5][i] = posterior(
                    ll - jmath::log(single.max(many)),
                    self.model.prior(r.len(i)),
                );
            }
        }
        let i = r.max_lod();
        let depth = r.depth();
        let k = r.ad[i + 1];
        if r.info.contains_key("NALOD") {
            let nd: u32 = r.normal_ad.iter().sum();
            let nk = r.normal_ad[i + 1];
            let nf = if nd == 0 { 0.0 } else { nk as f64 / nd as f64 };
            // Java's early `<` return must fall through when tumour AD is
            // zero (0/0 is NaN); `>=` is not its logical complement.
            let negligible = nf < 0.1 * (k as f64 / depth as f64);
            if !negligible {
                let q = r.values("MBQ").first().copied().unwrap_or(30.0);
                let error = jmath::pow(10.0, -q / 10.0);
                let cdf = if nk == 0 {
                    0.0
                } else if nk > nd {
                    1.0
                } else {
                    1.0 - regularized_beta(error, nk as f64, (nd - nk + 1) as f64)
                };
                let normal = if 1.0 - cdf < 0.001 {
                    1.0
                } else {
                    posterior(
                        -r.values("NALOD")[i] * std::f64::consts::LN_10,
                        self.model.artifact_prior,
                    )
                };
                p[9].fill(normal);
            }
        }
        if let Some(nc) = r.values("NCount").first() {
            let alt: u32 = r.ad[1..].iter().chain(&r.normal_ad[1..]).sum();
            if alt > 0 {
                p[10].fill(f64::from(*nc / alt as f64 >= f64::INFINITY));
            }
        }
        p[11].fill(f64::from(r.info.contains_key("PON")));
        p[12].fill(r.orientation);
        if r.info.contains_key("ECNT") && r.info.contains_key("ECNTH") {
            p[13].fill(f64::from(
                r.values("ECNT")[0] > 3.0 || r.values("ECNTH").iter().any(|&x| x > 2.0),
            ));
        }
        p[14].fill(f64::from(r.lod.iter().filter(|&&x| x > 5.0).count() > 1));
        let mfrl = r.values("MFRL");
        if mfrl.len() >= 2 {
            p[15].fill(f64::from((mfrl[0] - mfrl[1]).abs() > 10000.0));
        }
        if let (Some(rpa), Some(ru)) = (r.info.get("RPA"), r.info.get("RU")) {
            let rpa = ints(rpa);
            if rpa.len() >= 2 && ru.len() * rpa[0] as usize >= 8 && rpa[0].abs_diff(rpa[1]) == 1 {
                let alt = depth - r.ad[0];
                let artifact = regularized_beta(0.1, r.ad[1] as f64 + 1.0, r.ad[0] as f64 + 1.0);
                p[16].fill(posterior(
                    self.model.likelihood(depth, alt) - jmath::log(artifact),
                    self.model.prior(r.len(0)),
                ));
            }
        }
        if let Some(values) = r.phase.as_ref().and_then(|s| self.phased.get(s)) {
            p[17].fill(
                values
                    .iter()
                    .filter(|(pos, _)| pos.abs_diff(r.pos) <= 100)
                    .map(|x| x.1)
                    .fold(0.0, f64::max),
            );
        }
        if !pop.is_empty() && depth > 0 {
            let f = jmath::pow(10.0, -pop[i]);
            let germ = if f < 1e-10 {
                0.0
            } else if f > 1.0 - 1e-10 {
                1.0
            } else {
                let maf = (r.maf * depth as f64) / depth as f64;
                let het = -std::f64::consts::LN_2
                    + log_sum(&[
                        log_binomial_probability(depth as i32, k as i32, maf),
                        log_binomial_probability(depth as i32, k as i32, 1.0 - maf),
                    ])
                    - self.model.likelihood(depth, k);
                let af = r.gt[self.tumor].get("AF").map_or(0.0, |s| floats(s)[i]);
                let af = af * depth as f64 * (1.0 / depth as f64);
                let hom = if af < 0.9 { f64::NEG_INFINITY } else { 0.0 };
                let nl = r.values("NLOD").get(i).copied().unwrap_or(0.0) * -std::f64::consts::LN_10;
                let prior = self.model.prior(r.len(i));
                let not = log1mexp(prior);
                let g = log_sum(&[
                    jmath::log(2.0 * f * (1.0 - f)) + het + nl + not,
                    jmath::log(f * f) + hom + nl + not,
                ]);
                normalize(&[g, jmath::log((1.0 - f) * (1.0 - f)) + prior])[0]
            };
            p[18].fill(germ);
        }
        for v in &mut p {
            for x in v {
                *x = x.clamp(0.0, 1.0);
            }
        }
        let mut artifact = vec![0.0; n];
        let mut non = vec![0.0; n];
        for (j, v) in p.iter().enumerate() {
            if j == 0 {
                continue;
            }
            for (i, &x) in v.iter().enumerate() {
                if non_somatic(j) {
                    non[i] = f64::max(non[i], x);
                } else {
                    artifact[i] = f64::max(artifact[i], x);
                }
            }
        }
        // GATK iterates a HashMap keyed by ErrorType (JVM identity hashes).
        // The supported golden's complete output is gated for all six orders;
        // there is deliberately no claim that intermediate binary64 bits agree.
        let combined = (0..n)
            .map(|i| {
                let by_type = [p[0][i], artifact[i], non[i]];
                numerics::combined_error(by_type, self.error_order)
            })
            .collect();
        Probabilities {
            p,
            artifact,
            non_somatic: non,
            combined,
            steps,
        }
    }
    fn accumulate(
        &mut self,
        r: &Record,
        p: &Probabilities,
        phased: &mut HashMap<String, Vec<(u64, f64)>>,
    ) {
        self.strand.steps.extend_from_slice(&p.steps);
        if let Some(key) = &r.phase {
            let ignore = p.p[18].iter().any(|&x| x > 0.25);
            let prob =
                p.p.iter()
                    .enumerate()
                    .filter(|(i, _)| !non_somatic(*i) && *i != 17 && !(ignore && *i == 9))
                    .flat_map(|(_, v)| v)
                    .copied()
                    .fold(0.0, f64::max);
            phased.entry(key.clone()).or_default().push((r.pos, prob));
        }
        for i in 0..r.alts.len() {
            if p.artifact[i] > 0.9 {
                self.model.obvious += 1;
            } else if p.non_somatic[i] <= 0.9 {
                self.model
                    .data
                    .push(r.datum(i, p.artifact[i], p.non_somatic[i]));
            }
        }
    }
}
/// Paths for the supported two-sample default command. The merged Mutect stats
/// sidecar defaults to `<variant>.stats` and can be supplied explicitly.
pub struct Inputs<'a> {
    pub variant: &'a Path,
    pub reference: &'a Path,
    pub orientation_priors: &'a Path,
    pub segmentation: &'a Path,
    pub contamination: &'a Path,
    pub stats: Option<&'a Path>,
    pub output: &'a Path,
}
fn read_vcf(path: &Path) -> io::Result<Vec<String>> {
    // Keep individual lines, rather than a large contiguous String buffer.
    BufReader::new(nomics_io::gzip::GzReader::open(path)?)
        .lines()
        .collect()
}
fn prepare(
    records: &mut [Record],
    inputs: &Inputs<'_>,
    tumor: usize,
    sample: &str,
) -> io::Result<()> {
    let mut priors: BTreeMap<String, [f64; 12]> = BTreeMap::new();
    let mut members = 0;
    tar::read(inputs.orientation_priors, |_, bytes| {
        members += 1;
        if members > 1 {
            return Err(invalid(
                "only a single-sample orientation archive is supported",
            ));
        }
        let text = std::str::from_utf8(bytes).map_err(io::Error::other)?;
        if !text
            .lines()
            .any(|l| l == format!("#<METADATA>SAMPLE={sample}"))
        {
            return Err(invalid("orientation prior sample mismatch"));
        }
        for l in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.starts_with("context\t"))
        {
            let f: Vec<_> = l.split('\t').collect();
            if f.len() != 16 {
                return Err(invalid("bad prior row"));
            }
            let mut p = [0.0; 12];
            for (p, f) in p.iter_mut().zip(&f[2..14]) {
                *p = f.parse().map_err(|_| invalid("bad prior probability"))?;
            }
            priors.insert(f[0].to_owned(), p);
        }
        Ok(())
    })?;
    let seg = fs::read_to_string(inputs.segmentation)?;
    if !seg
        .lines()
        .any(|l| l == format!("#<METADATA>SAMPLE={sample}"))
    {
        return Err(invalid("segmentation sample mismatch"));
    }
    let segments: Vec<(String, u64, u64, f64)> = seg
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("contig\t"))
        .map(|l| {
            let f: Vec<_> = l.split('\t').collect();
            if f.len() != 4 {
                return Err(invalid("bad segmentation row"));
            }
            Ok((
                f[0].to_owned(),
                f[1].parse().map_err(io::Error::other)?,
                f[2].parse().map_err(io::Error::other)?,
                f[3].parse().map_err(io::Error::other)?,
            ))
        })
        .collect::<io::Result<_>>()?;
    let fai = nomics_io::fasta::read_fai(&PathBuf::from(format!(
        "{}.fai",
        inputs.reference.display()
    )))?;
    let mut fasta = File::open(inputs.reference)?;
    let base = |b: u8| b"ACGT".iter().position(|&x| x == b.to_ascii_uppercase());
    for r in records {
        let end = r.pos + r.fields[3].len() as u64 - 1;
        let overlapping: Vec<_> = segments
            .iter()
            .filter(|s| s.0 == r.fields[0] && s.1 <= end && s.2 >= r.pos)
            .collect();
        if overlapping.len() > 1 {
            return Err(invalid("overlapping tumour segments are not supported"));
        }
        if let Some(s) = overlapping.first() {
            r.maf = s.3;
        }
        if r.alts.iter().any(|a| a.len() != r.fields[3].len())
            || r.gt[tumor]
                .get("GT")
                .is_some_and(|s| s == "0/0" || s == "0|0")
        {
            continue;
        }
        let e = fai
            .iter()
            .find(|e| e.name == r.fields[0])
            .ok_or_else(|| invalid("contig absent from FASTA"))?;
        let i = r.max_lod();
        let f1 = r.gt[tumor].get("F1R2").map(|s| ints(s));
        let f2 = r.gt[tumor].get("F2R1").map(|s| ints(s));
        let (Some(f1), Some(f2)) = (f1, f2) else {
            continue;
        };
        let k = f1[i + 1] + f2[i + 1];
        let d = f1[0] + f2[0] + k;
        for (j, &alt) in r.alts[i].as_bytes().iter().enumerate() {
            let pos = r.pos + j as u64;
            if pos < 2 || pos >= e.len {
                continue;
            }
            let mut ctx = [0; 3];
            for (off, b) in ctx.iter_mut().enumerate() {
                let idx = pos - 2 + off as u64;
                fasta.seek(SeekFrom::Start(
                    e.offset + idx / e.line_bases * e.line_width + idx % e.line_bases,
                ))?;
                fasta.read_exact(std::slice::from_mut(b))?;
                *b = b.to_ascii_uppercase();
            }
            let (Some(reference), Some(alt)) = (base(ctx[1]), base(alt)) else {
                continue;
            };
            let Some(prior) = std::str::from_utf8(&ctx).ok().and_then(|s| priors.get(s)) else {
                continue;
            };
            r.orientation =
                r.orientation
                    .max(artifact_probability(reference, alt, d, k, f1[i + 1], prior));
        }
    }
    Ok(())
}
#[derive(Default)]
struct Stats {
    pass: u32,
    tp: f64,
    fp: f64,
    fn_: f64,
    filter_fp: [f64; 19],
    filter_fn: [f64; 19],
}
impl Stats {
    fn record(&mut self, p: &Probabilities, t: f64) {
        for (i, &combined) in p.combined.iter().enumerate() {
            let filtered = combined > t;
            if filtered {
                self.fn_ += 1.0 - combined;
            } else {
                self.pass += 1;
                self.fp += combined;
                self.tp += 1.0 - combined;
            }
            for (j, v) in p.p.iter().enumerate() {
                let Some(&v) = v.get(i) else { continue };
                if v > 1e-10 && v > t - 1e-10 {
                    self.filter_fn[j] += 1.0 - combined;
                } else if !filtered {
                    self.filter_fp[j] += v;
                }
            }
        }
    }
    fn text(&self, model: &Model, t: f64) -> String {
        let mut s = model.metadata();
        s += &format!(
            "#<METADATA>threshold={}\n#<METADATA>fdr={}\n#<METADATA>sensitivity={}\nfilter\tFP\tFDR\tFN\tFNR\n",
            numerics::round(t, 3),
            numerics::round(self.fp / self.pass as f64, 3),
            numerics::round(self.tp / (self.tp + self.fn_), 3)
        );
        for (i, name) in NAMES.iter().enumerate() {
            if self.filter_fp[i] > 0.0 || self.filter_fn[i] > 0.0 {
                s += name;
                for x in [
                    self.filter_fp[i],
                    self.filter_fp[i] / self.pass as f64,
                    self.filter_fn[i],
                    self.filter_fn[i] / (self.tp + self.fn_),
                ] {
                    s.push('\t');
                    s += &numerics::round(x, 2);
                }
                s.push('\n');
            }
        }
        s
    }
}
fn output_record(r: &Record, p: &Probabilities, t: f64) -> String {
    let t = t.clamp(1e-10, 1.0 - 1e-10);
    let mut filters = BTreeMap::new();
    let mut statuses = Vec::new();
    for i in 0..r.alts.len() {
        let names: Vec<_> =
            p.p.iter()
                .enumerate()
                .filter(|(j, v)| allele_filter(*j) && v.get(i).is_some_and(|&x| x > t))
                .map(|(j, _)| NAMES[j])
                .collect();
        statuses.push(if names.is_empty() {
            "SITE".to_owned()
        } else {
            names.join(",")
        });
    }
    for (j, v) in p.p.iter().enumerate() {
        if v.is_empty() {
            continue;
        }
        if allele_filter(j) {
            if v.iter().all(|&x| x > t) {
                filters.insert(NAMES[j], 1.0);
            }
        } else if v[0] > t {
            filters.insert(NAMES[j], v[0]);
        }
    }
    if filters.is_empty() && statuses.iter().all(|s| s != "SITE") {
        filters.insert("FAIL", 1.0);
    }
    let max = filters.values().copied().fold(0.0, f64::max).min(0.1);
    let names: Vec<_> = filters
        .iter()
        .filter(|(_, p)| **p >= max)
        .map(|(s, _)| *s)
        .collect();
    let mut fields = r.fields.clone();
    fields[6] = if names.is_empty() {
        "PASS".to_owned()
    } else {
        names.join(";")
    };
    let mut info = r.info.clone();
    info.insert("AS_FilterStatus".to_owned(), statuses.join("|"));
    for (j, key, required) in [
        (12, "ROQ", true),
        (18, "GERMQ", r.info.contains_key("POPAF")),
        (
            16,
            "STRQ",
            r.info.contains_key("RPA") && r.info.contains_key("RU"),
        ),
    ] {
        if required {
            info.insert(key.to_owned(), numerics::qual(p.p[j][0]).to_string());
        }
    }
    fields[7] = info
        .iter()
        .map(|(k, v)| {
            if v.is_empty() {
                k.clone()
            } else {
                format!("{k}={v}")
            }
        })
        .collect::<Vec<_>>()
        .join(";");
    fields.join("\t") + "\n"
}
/// Run the four traversals and write valid BGZF. No tabix index is emitted.
pub fn run(inputs: &Inputs<'_>) -> io::Result<()> {
    run_with_error_order(inputs, [0, 1, 2])
}
fn run_with_error_order(inputs: &Inputs<'_>, error_order: [usize; 3]) -> io::Result<()> {
    let text = read_vcf(inputs.variant)?;
    let header: Vec<_> = text
        .iter()
        .map(String::as_str)
        .filter(|s| s.starts_with('#'))
        .collect();
    let samples: Vec<_> = header
        .iter()
        .find(|s| s.starts_with("#CHROM\t"))
        .ok_or_else(|| invalid("missing VCF header"))?
        .split('\t')
        .skip(9)
        .collect();
    let normals: Vec<_> = header
        .iter()
        .filter_map(|s| s.strip_prefix("##normal_sample="))
        .collect();
    if samples.len() != 2 || normals.len() != 1 {
        return Err(invalid("only one tumour and one normal sample supported"));
    }
    let normal = samples
        .iter()
        .position(|s| *s == normals[0])
        .ok_or_else(|| invalid("normal sample absent"))?;
    let tumor = 1 - normal;
    let inferred = PathBuf::from(format!("{}.stats", inputs.variant.display()));
    let stats_text = fs::read_to_string(inputs.stats.unwrap_or(&inferred))?;
    let callable = stats_text
        .lines()
        .find_map(|s| s.strip_prefix("callable\t"))
        .ok_or_else(|| invalid("callable statistic required"))?
        .parse::<f64>()
        .map_err(io::Error::other)?;
    if callable < 1.0 {
        return Err(invalid("no callable sites is not supported"));
    }
    let contamination_text = fs::read_to_string(inputs.contamination)?;
    let first = contamination_text
        .lines()
        .find(|s| !s.starts_with('#') && !s.starts_with("sample\t"))
        .ok_or_else(|| invalid("empty contamination table"))?;
    let f: Vec<_> = first.split('\t').collect();
    if f.len() != 3 || f[0] != samples[tumor] {
        return Err(invalid("contamination sample mismatch"));
    }
    let contamination = f[1].parse().map_err(io::Error::other)?;
    let mut records: Vec<_> = text
        .iter()
        .map(String::as_str)
        .filter(|s| !s.starts_with('#'))
        .map(|s| Record::parse(s, tumor, normal))
        .collect::<io::Result<_>>()?;
    prepare(&mut records, inputs, tumor, samples[tumor])?;
    let mut engine = Engine {
        model: Model::new(callable),
        strand: Strand::new(),
        contamination,
        tumor,
        phased: HashMap::new(),
        error_order,
    };
    for r in &records {
        for i in 0..r.alts.len() {
            engine.model.register_length(r.len(i));
        }
    }
    let mut threshold = 0.1;
    for pass in 0..3 {
        let mut probabilities = Vec::new();
        let mut phased = HashMap::new();
        for r in &records {
            let p = engine.probabilities(r);
            engine.accumulate(r, &p, &mut phased);
            probabilities.extend(p.combined);
        }
        if pass < 2 {
            engine.strand.learn()?;
            engine.phased = phased;
            engine.model.learn();
        }
        threshold = numerics::threshold(probabilities);
        eprintln!(
            "FilterMutectCalls pass {pass}: threshold {}",
            double_to_string(threshold)
        );
    }
    let mut metadata: BTreeMap<String, String> = header
        .iter()
        .filter(|s| {
            s.starts_with("##")
                && !s.starts_with("##filtering_status=")
                && !s.starts_with("##fileformat=")
        })
        .map(|s| (s.to_string(), s.to_string()))
        .collect();
    for s in include_str!("filtering/header.txt").lines() {
        metadata.insert(s.to_owned(), s.to_owned());
    }
    let mut out = nomics_io::bgzf::ParallelBgzfWriter::new(File::create(inputs.output)?, 2);
    out.write_all(b"##fileformat=VCFv4.2\n")?;
    let mut lines: Vec<_> = metadata.values().collect();
    let contig_order: HashMap<_, _> = header
        .iter()
        .filter(|s| s.starts_with("##contig="))
        .enumerate()
        .map(|(i, s)| (*s, i))
        .collect();
    lines.sort_by(|a, b| {
        if a.starts_with("##contig=") && b.starts_with("##contig=") {
            contig_order[a.as_str()].cmp(&contig_order[b.as_str()])
        } else {
            a.cmp(b)
        }
    });
    for s in lines {
        writeln!(out, "{s}")?;
    }
    writeln!(
        out,
        "{}",
        header.iter().find(|s| s.starts_with("#CHROM\t")).unwrap()
    )?;
    let mut stats = Stats::default();
    let mut pass = 0;
    for r in &records {
        let p = engine.probabilities(r);
        stats.record(&p, threshold - 1e-10);
        let line = output_record(r, &p, threshold);
        if line.split('\t').nth(6) == Some("PASS") {
            pass += 1;
        }
        out.write_all(line.as_bytes())?;
    }
    out.finish()?;
    fs::write(
        format!("{}.filteringStats.tsv", inputs.output.display()),
        stats.text(&engine.model, threshold),
    )?;
    eprintln!("FilterMutectCalls: {} records, {pass} PASS", records.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn error_type_hash_iteration_orders_have_identical_golden_boundaries() {
        let home = PathBuf::from(std::env::var_os("HOME").unwrap_or_default());
        let golden = std::env::var_os("NOMICS_SOMATIC_GOLDENS")
            .map(PathBuf::from)
            .unwrap_or_else(|| home.join("omics/runs/somatic/S1-r1/golden"));
        if !golden.is_dir() {
            eprintln!(
                "SKIP FilterMutectCalls error-order gate: {} absent",
                golden.display()
            );
            return;
        }
        let reference = std::env::var_os("NOMICS_SOMATIC_FASTA")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                home.join(
                    "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
                )
            });
        let x = "HCC1395_WES_vs_HCC1395BL_WES";
        let input = |dir: &str, suffix: &str| golden.join(format!("{dir}-{x}/{x}.mutect2{suffix}"));
        let expected = read_vcf(&input("FILTERMUTECTCALLS", ".filtered.vcf.gz")).unwrap();
        let expected: Vec<_> = expected
            .into_iter()
            .filter(|s| !s.starts_with("##GATKCommandLine="))
            .collect();
        let stats = fs::read(input(
            "FILTERMUTECTCALLS",
            ".filtered.vcf.gz.filteringStats.tsv",
        ))
        .unwrap();
        let output = std::env::temp_dir().join(format!(
            "nomics-filter-orders-{}.vcf.gz",
            std::process::id()
        ));
        for order in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            run_with_error_order(
                &Inputs {
                    variant: &input("MERGE_MUTECT2", ".vcf.gz"),
                    stats: Some(&input("MERGEMUTECTSTATS", ".vcf.gz.stats")),
                    reference: &reference,
                    orientation_priors: &input(
                        "LEARNREADORIENTATIONMODEL",
                        ".artifactprior.tar.gz",
                    ),
                    segmentation: &input("CALCULATECONTAMINATION", ".segmentation.table"),
                    contamination: &input("CALCULATECONTAMINATION", ".contamination.table"),
                    output: &output,
                },
                order,
            )
            .unwrap();
            let actual: Vec<_> = read_vcf(&output)
                .unwrap()
                .into_iter()
                .filter(|s| !s.starts_with("##GATKCommandLine="))
                .collect();
            assert!(actual == expected, "VCF error-type order {order:?}");
            assert_eq!(
                fs::read(format!("{}.filteringStats.tsv", output.display())).unwrap(),
                stats,
                "stats error-type order {order:?}"
            );
        }
        fs::remove_file(format!("{}.filteringStats.tsv", output.display())).unwrap();
        fs::remove_file(output).unwrap();
    }
    #[test]
    fn long_indel_priors_stay_at_first_access_value() {
        let mut model = Model::new(100.0);
        model.register_length(60);
        let first = model.prior(60);
        model.data.push(cluster::Datum {
            lod: 20.0,
            artifact: 0.0,
            non_seq: 0.0,
            k: 30,
            n: 50,
            len: 0,
        });
        model.learn();
        assert_eq!(model.prior(60), first);
    }
}
