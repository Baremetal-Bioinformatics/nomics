//! GATK 4.6.2.0 `GenotypeGVCFs` as sarek 3.10.0 runs it for
//! `--joint_germline` on one sample: GenomicsDBImport → GenotypeGVCFs
//! (`--dbsnp`, by-variant traversal) → `bcftools sort` → `MergeVcfs`.
//!
//! Each gVCF variant record (a record with an alt besides `<NON_REF>`) goes
//! through `ReferenceConfidenceVariantContextMerger.merge` (drop
//! `<NON_REF>`, remap PL/AD, median the INFO values), then
//! `GenotypeGVCFsEngine.regenotypeVC` (re-genotype from PLs, finalize
//! `RAW_MQandDP`, reverse-trim, annotate, dbSNP), then the htslib/htsjdk
//! text round trip of [`htsfmt`]. Reference blocks never become variants.
//! Scope: one diploid sample, the default `StandardAnnotation` group. See
//! `docs/spec/joint.md`.

pub mod htsfmt;

use std::collections::{BTreeMap, HashMap};
use std::io::{self, Write};
use std::path::Path;

use nomics_hc::annotate::{fisher_strand, strand_odds_ratio};
use nomics_hc::dbsnp::Dbsnp;
use nomics_hc::genotype::{GenotypingEngine, NO_CALL_ALLELE};
use nomics_hc::javafmt::format_fixed;
use nomics_hc::jrandom::{GATK_SEED, JavaRandom};
use nomics_hc::vcf::{format_qual, format_vcf_double};

use htsfmt::{bcf_float, merged_qual};

const NON_REF: &[u8] = b"<NON_REF>";
const SPAN_DEL: &[u8] = b"*";

/// INFO keys `GenotypeGVCFs` writes as VCF Float (re-encoded by bcftools).
const FLOAT_INFO: &[&str] = &[
    "AF",
    "BaseQRankSum",
    "ExcessHet",
    "FS",
    "InbreedingCoeff",
    "MLEAF",
    "MQ",
    "MQRankSum",
    "QD",
    "ReadPosRankSum",
    "SOR",
];

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// A merged INFO value, typed as the merger and annotators produce it.
#[derive(Debug, Clone, PartialEq)]
enum Value {
    /// `Integer.parseInt` from the gVCF.
    Int(i64),
    /// `Double.parseDouble` from the gVCF; written with `formatVCFDouble`.
    Double(f64),
    Ints(Vec<i64>),
    Doubles(Vec<f64>),
    /// An annotation that formats itself.
    Str(String),
    Flag,
}

impl Value {
    fn render(&self) -> Option<String> {
        let join = |v: Vec<String>| v.join(",");
        Some(match self {
            Value::Int(x) => x.to_string(),
            Value::Double(x) => format_vcf_double(*x),
            Value::Ints(v) => join(v.iter().map(|x| x.to_string()).collect()),
            Value::Doubles(v) => join(v.iter().map(|&x| format_vcf_double(x)).collect()),
            Value::Str(s) => s.clone(),
            Value::Flag => return None,
        })
    }
}

/// The single sample of one gVCF record.
#[derive(Debug, Clone)]
struct Sample {
    /// Allele indices into the record's alleles; `None` is `.`.
    gt: Vec<Option<usize>>,
    phased: bool,
    ad: Option<Vec<i32>>,
    dp: Option<i32>,
    gq: Option<i32>,
    pl: Option<Vec<i32>>,
    /// Other FORMAT keys, verbatim (PGT, PID, PS, SB, MIN_DP, ...).
    extended: BTreeMap<String, String>,
}

struct GvcfRecord<'a> {
    contig: &'a str,
    pos: i64,
    alleles: Vec<Vec<u8>>,
    info: Vec<(&'a str, Option<&'a str>)>,
    sample: Sample,
}

fn parse_record(line: &str) -> io::Result<GvcfRecord<'_>> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() != 10 {
        return Err(invalid(format!("expected one sample: {line}")));
    }
    let pos = f[1]
        .parse()
        .map_err(|_| invalid(format!("bad POS: {line}")))?;
    let mut alleles = vec![f[3].as_bytes().to_vec()];
    alleles.extend(f[4].split(',').map(|a| a.as_bytes().to_vec()));
    let info = if f[7] == "." {
        Vec::new()
    } else {
        f[7].split(';')
            .map(|kv| match kv.split_once('=') {
                Some((k, v)) => (k, Some(v)),
                None => (kv, None),
            })
            .collect()
    };
    let ints = |s: &str| -> io::Result<Vec<i32>> {
        s.split(',')
            .map(|x| x.parse().map_err(|_| invalid(format!("bad integer list: {line}"))))
            .collect()
    };
    let mut sample = Sample {
        gt: Vec::new(),
        phased: false,
        ad: None,
        dp: None,
        gq: None,
        pl: None,
        extended: BTreeMap::new(),
    };
    for (key, value) in f[8].split(':').zip(f[9].split(':')) {
        if value == "." && key != "GT" {
            continue;
        }
        match key {
            "GT" => {
                sample.phased = value.contains('|');
                sample.gt = value
                    .split(['/', '|'])
                    .map(|a| if a == "." { Ok(None) } else { a.parse().map(Some) })
                    .collect::<Result<_, _>>()
                    .map_err(|_| invalid(format!("bad GT: {line}")))?;
            }
            "AD" => sample.ad = Some(ints(value)?),
            "DP" => sample.dp = Some(ints(value)?[0]),
            "GQ" => sample.gq = Some(ints(value)?[0]),
            "PL" => sample.pl = Some(ints(value)?),
            _ => {
                sample.extended.insert(key.to_string(), value.to_string());
            }
        }
    }
    Ok(GvcfRecord {
        contig: f[0],
        pos,
        alleles,
        info,
        sample,
    })
}

#[inline]
fn genotype_index(a: usize, b: usize) -> usize {
    let (low, high) = if a <= b { (a, b) } else { (b, a) };
    high * (high + 1) / 2 + low
}

/// Header INFO types, for `parseNumericInfoAttributeValue`.
#[derive(Clone, Copy, PartialEq)]
enum InfoType {
    Integer,
    Float,
    Other,
}

fn info_types(header: &[String]) -> HashMap<String, InfoType> {
    let mut types = HashMap::new();
    for line in header {
        let Some(body) = line.strip_prefix("##INFO=<") else {
            continue;
        };
        let field = |name: &str| {
            body.split(',')
                .find_map(|kv| kv.strip_prefix(name).and_then(|v| v.strip_prefix('=')))
        };
        if let Some(id) = field("ID") {
            let t = match field("Type") {
                Some("Integer") => InfoType::Integer,
                Some("Float") => InfoType::Float,
                _ => InfoType::Other,
            };
            types.insert(id.to_string(), t);
        }
    }
    types
}

/// One output record before text rendering.
struct Called {
    contig: String,
    pos: i64,
    id: String,
    alleles: Vec<Vec<u8>>,
    qual: f64,
    info: BTreeMap<String, Value>,
    gt: [usize; 2],
    phased: bool,
    ad: Option<Vec<i32>>,
    dp: Option<i32>,
    gq: Option<i32>,
    pl: Vec<i32>,
    extended: BTreeMap<String, String>,
}

pub struct Joint {
    engine: GenotypingEngine,
    rng: JavaRandom,
    contigs: HashMap<String, usize>,
    types: HashMap<String, InfoType>,
}

impl Joint {
    pub fn new(header: &[String]) -> Self {
        let mut contigs = HashMap::new();
        for line in header {
            if let Some(body) = line.strip_prefix("##contig=<ID=") {
                let id = body.split([',', '>']).next().unwrap_or_default();
                let next = contigs.len();
                contigs.insert(id.to_string(), next);
            }
        }
        Self {
            engine: GenotypingEngine::new(),
            rng: JavaRandom::new(GATK_SEED),
            contigs,
            types: info_types(header),
        }
    }

    /// Genotype one gVCF data line; `None` when nothing is emitted.
    fn line(&mut self, line: &str) -> io::Result<Option<Called>> {
        let rec = parse_record(line)?;
        if !rec.alleles[1..].iter().any(|a| a != NON_REF) {
            return Ok(None);
        }
        self.genotype(&rec)
    }

    fn genotype(&mut self, rec: &GvcfRecord<'_>) -> io::Result<Option<Called>> {
        // ---- ReferenceConfidenceVariantContextMerger.merge
        let non_ref = rec.alleles.iter().position(|a| a == NON_REF);
        let target: Vec<usize> = (0..rec.alleles.len())
            .filter(|&i| Some(i) != non_ref)
            .collect();
        // getIndexesOfRelevantAllelesForGVCF: each target allele's index in
        // the record (all present here, so NON_REF is never substituted).
        let g = &rec.sample;
        let pl = match &g.pl {
            Some(pl) => {
                let mut out = Vec::with_capacity(target.len() * (target.len() + 1) / 2);
                for high in 0..target.len() {
                    for low in 0..=high {
                        out.push(pl[genotype_index(target[low], target[high])]);
                    }
                }
                out
            }
            None => return Err(invalid(format!("variant without PL at {}:{}", rec.contig, rec.pos))),
        };
        let ad = g
            .ad
            .as_ref()
            .map(|ad| target.iter().map(|&i| ad[i]).collect::<Vec<i32>>());

        let mut info: BTreeMap<String, Value> = BTreeMap::new();
        let mut raw_mq: Option<(i64, i64)> = None;
        for &(key, value) in &rec.info {
            if key == "RAW_MQandDP" {
                let v = value.unwrap_or_default();
                let (a, b) = v
                    .split_once(',')
                    .ok_or_else(|| invalid(format!("bad RAW_MQandDP {v}")))?;
                raw_mq = Some((
                    a.parse().map_err(|_| invalid(format!("bad RAW_MQandDP {v}")))?,
                    b.parse().map_err(|_| invalid(format!("bad RAW_MQandDP {v}")))?,
                ));
                continue;
            }
            let Some(value) = value else {
                continue; // a Flag is not numeric: the merger drops it
            };
            if value.contains(',') {
                continue; // multi-valued keys here are only the stale MLEAC/MLEAF
            }
            let parsed = match self.types.get(key) {
                Some(InfoType::Integer) => value.parse().ok().map(Value::Int),
                Some(InfoType::Float) => value.parse().ok().map(Value::Double),
                _ => None,
            };
            if let Some(v) = parsed {
                info.insert(key.to_string(), v);
            }
        }
        let depth = match rec.info.iter().find(|(k, _)| *k == "DP") {
            Some((_, Some(v))) => v.parse::<i64>().unwrap_or(0),
            _ => g
                .extended
                .get("MIN_DP")
                .and_then(|v| v.parse().ok())
                .or(g.dp.map(i64::from))
                .unwrap_or(0),
        };
        if depth > 0 {
            info.insert("DP".into(), Value::Str(depth.to_string()));
        }
        for k in ["AC", "AF", "AN", "MLEAC", "MLEAF", "END"] {
            info.remove(k);
        }

        // ---- GenotypeGVCFsEngine.regenotypeVC
        if depth <= 0 || target.len() < 2 {
            return Ok(None);
        }
        let alleles: Vec<Vec<u8>> = target.iter().map(|&i| rec.alleles[i].clone()).collect();
        let contig = *self
            .contigs
            .get(rec.contig)
            .ok_or_else(|| invalid(format!("contig {} not in header", rec.contig)))?;
        let Some(call) = self
            .engine
            .calculate(contig, rec.pos, rec.pos, &alleles, &pl)
        else {
            return Ok(None);
        };
        let kept = &call.alleles_kept;
        info.insert(
            "MLEAC".into(),
            Value::Ints(call.mleac.iter().map(|&x| i64::from(x)).collect()),
        );
        info.insert("MLEAF".into(), Value::Doubles(call.mleaf.clone()));
        // finalizeAnnotations: RMSMappingQuality from RAW_MQandDP.
        if let Some((sum_sq, n)) = raw_mq {
            let mq = (sum_sq as f64 / n as f64).sqrt();
            info.insert("MQ".into(), Value::Str(format_fixed(mq, 2)));
        }
        let mut out_alleles: Vec<Vec<u8>> = kept.iter().map(|&i| alleles[i].clone()).collect();
        reverse_trim(&mut out_alleles);
        let out_ad = ad.map(|ad| kept.iter().map(|&i| ad[i]).collect::<Vec<i32>>());
        let qual = -10.0 * call.log10_p_error;

        // ---- annotateContext (no reads: genotype- and SB-based annotations)
        let gt = call.gt;
        let called: Vec<usize> = gt.iter().copied().filter(|&a| a != NO_CALL_ALLELE).collect();
        let n_alt = out_alleles.len() - 1;
        let is_het = called.len() == 2 && called[0] != called[1];
        let is_hom_var = called.len() == 2 && called[0] == called[1] && called[0] != 0;
        if !called.is_empty() {
            let mut ac = vec![0_i64; n_alt];
            for &a in &called {
                if a != 0 {
                    ac[a - 1] += 1;
                }
            }
            let an = called.len() as i64;
            info.insert("AC".into(), Value::Ints(ac.clone()));
            info.insert(
                "AF".into(),
                Value::Doubles(ac.iter().map(|&x| x as f64 / an as f64).collect()),
            );
            info.insert("AN".into(), Value::Ints(vec![an]));
            // ExcessHet for a single diploid sample: the exact test's p-value is 1.
            info.insert("ExcessHet".into(), Value::Str("0.0000".into()));
        }
        if let Some(sb) = g.extended.get("SB") {
            let v: Vec<i32> = sb
                .split(',')
                .map(|x| x.parse().map_err(|_| invalid(format!("bad SB {sb}"))))
                .collect::<io::Result<_>>()?;
            let table = |min: i32| {
                if v.iter().sum::<i32>() > min {
                    [[v[0], v[1]], [v[2], v[3]]]
                } else {
                    [[0, 0], [0, 0]]
                }
            };
            info.insert(
                "FS".into(),
                Value::Str(format_fixed(fisher_strand(table(2)), 3)),
            );
            info.insert(
                "SOR".into(),
                Value::Str(format_fixed(strand_odds_ratio(table(0)), 3)),
            );
        }
        // QualByDepth.getDepth.
        let mut qd_depth = 0;
        let mut ad_restricted = 0;
        if is_het || is_hom_var {
            match out_ad.as_ref().map(|ad| (ad.iter().sum::<i32>(), ad[0])) {
                Some((total, r)) if total != 0 => {
                    if total - r > 1 {
                        ad_restricted += total;
                    }
                    qd_depth += total;
                }
                _ => qd_depth += g.dp.unwrap_or(0),
            }
        }
        if ad_restricted > 0 {
            qd_depth = ad_restricted;
        }
        if qd_depth != 0 {
            let mut qd = qual / f64::from(qd_depth);
            if qd >= 35.0 {
                qd = 30.0 + self.rng.next_gaussian() * 3.0;
            }
            info.insert("QD".into(), Value::Str(format_fixed(qd, 2)));
        }

        // GenomicsDB does not carry the gVCF ID column through the import;
        // the ID comes from dbSNP alone ([`annotate_dbsnp`]).
        let id = ".".to_string();

        // ---- cleanupGenotypeAnnotations(keepSB = false)
        let mut extended = g.extended.clone();
        extended.remove("SB");
        let mut dp = g.dp;
        if let Some(min_dp) = extended.remove("MIN_DP") {
            dp = min_dp.parse().ok();
        }
        if is_hom_var && extended.contains_key("PGT") {
            extended.insert("PGT".into(), "1|1".into());
        }

        Ok(Some(Called {
            contig: rec.contig.to_string(),
            pos: rec.pos,
            id,
            alleles: out_alleles,
            qual,
            info,
            gt,
            phased: g.phased,
            ad: out_ad,
            dp,
            gq: Some(call.gq),
            pl: call.pl,
            extended,
        }))
    }
}

/// `GATKVariantContextUtils.reverseTrimAlleles`: drop trailing bases shared
/// by every allele, keeping at least one base each. `*` and symbolic
/// alleles do not take part.
fn reverse_trim(alleles: &mut [Vec<u8>]) {
    let real = |a: &Vec<u8>| a.as_slice() != SPAN_DEL && !a.starts_with(b"<");
    let mut clip = 0;
    'outer: loop {
        let mut base = None;
        for a in alleles.iter().filter(|a| real(a)) {
            if a.len() <= clip + 1 {
                break 'outer;
            }
            let b = a[a.len() - 1 - clip];
            match base {
                None => base = Some(b),
                Some(x) if x != b => break 'outer,
                _ => {}
            }
        }
        clip += 1;
    }
    if clip > 0 {
        for a in alleles.iter_mut().filter(|a| real(a)) {
            a.truncate(a.len() - clip);
        }
    }
}

fn render(c: &Called) -> String {
    let mut s = String::with_capacity(256);
    s.push_str(&c.contig);
    s.push('\t');
    s.push_str(&c.pos.to_string());
    s.push('\t');
    s.push_str(&c.id);
    s.push('\t');
    s.push_str(std::str::from_utf8(&c.alleles[0]).expect("ASCII"));
    s.push('\t');
    let alts: Vec<&str> = c.alleles[1..]
        .iter()
        .map(|a| std::str::from_utf8(a).expect("ASCII"))
        .collect();
    s.push_str(&alts.join(","));
    s.push('\t');
    s.push_str(&merged_qual(&format_qual(c.qual)));
    s.push_str("\t.\t");
    let mut fields = Vec::with_capacity(c.info.len());
    for (k, v) in &c.info {
        match v.render() {
            None => fields.push(k.clone()),
            Some(text) => {
                let text = if FLOAT_INFO.contains(&k.as_str()) {
                    text.split(',').map(bcf_float).collect::<Vec<_>>().join(",")
                } else {
                    text
                };
                fields.push(format!("{k}={text}"));
            }
        }
    }
    s.push_str(&fields.join(";"));
    // htsjdk: GT first, then every other key sorted.
    let mut keys: Vec<(String, String)> = Vec::new();
    if let Some(ad) = &c.ad {
        keys.push(("AD".into(), join_ints(ad)));
    }
    if let Some(dp) = c.dp {
        keys.push(("DP".into(), dp.to_string()));
    }
    if let Some(gq) = c.gq {
        keys.push(("GQ".into(), gq.min(99).to_string()));
    }
    keys.push(("PL".into(), join_ints(&c.pl)));
    for (k, v) in &c.extended {
        keys.push((k.clone(), v.clone()));
    }
    keys.sort_by(|a, b| a.0.cmp(&b.0));
    let sep = if c.phased { "|" } else { "/" };
    let gt: Vec<String> = c
        .gt
        .iter()
        .map(|&a| if a == NO_CALL_ALLELE { ".".into() } else { a.to_string() })
        .collect();
    s.push('\t');
    s.push_str("GT");
    for (k, _) in &keys {
        s.push(':');
        s.push_str(k);
    }
    s.push('\t');
    s.push_str(&gt.join(sep));
    for (_, v) in &keys {
        s.push(':');
        s.push_str(v);
    }
    s
}

fn join_ints(v: &[i32]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

/// Header lines GenotypeGVCFs adds (`setupVCFWriter`) and the `PASS` line
/// bcftools adds, verbatim.
const ADDED_HEADER: &[&str] = &[
    r#"FILTER=<ID=PASS,Description="All filters passed">"#,
    r#"FILTER=<ID=LowQual,Description="Low quality">"#,
    r#"FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">"#,
    r#"INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">"#,
    r#"INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">"#,
    r#"INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">"#,
    r#"INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">"#,
    r#"INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">"#,
    r#"INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">"#,
    r#"INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">"#,
    r#"INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">"#,
    r#"INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">"#,
    // The run-specific command line is outside the gate (docs/GATES.md).
    r#"GATKCommandLine=<ID=GenomicsDBImport,CommandLine="nomics joint-genotype",Version="4.6.2.0",Date="">"#,
    "source=GenomicsDBImport",
];

/// The published header: the gVCF's lines without GVCF blocks, HaplotypeCaller's
/// command line and source, plus [`ADDED_HEADER`], in htsjdk `VCFHeader`
/// order (fileformat first, lines sorted by text, contigs as one block in
/// dictionary order where `contig=` sorts).
fn joint_header(gvcf_header: &[String]) -> io::Result<String> {
    let mut meta: Vec<String> = Vec::new();
    let mut contigs: Vec<&str> = Vec::new();
    let mut columns = None;
    for line in gvcf_header {
        if let Some(c) = line.strip_prefix("##contig=") {
            contigs.push(c);
        } else if let Some(m) = line.strip_prefix("##") {
            if !(m.starts_with("fileformat=")
                || m.starts_with("GVCFBlock")
                || m.starts_with("GATKCommandLine=")
                || m.starts_with("source="))
            {
                meta.push(m.to_string());
            }
        } else if line.starts_with("#CHROM") {
            columns = Some(line.as_str());
        }
    }
    for line in ADDED_HEADER {
        let id = |l: &str| l.split(',').next().map(str::to_string);
        if !meta.iter().any(|m| id(m) == id(line)) {
            meta.push((*line).to_string());
        }
    }
    meta.sort();
    let mut s = String::from("##fileformat=VCFv4.2\n");
    let mut contigs_done = false;
    for m in &meta {
        if !contigs_done && m.as_str() > "contig=" {
            for c in &contigs {
                s.push_str("##contig=");
                s.push_str(c);
                s.push('\n');
            }
            contigs_done = true;
        }
        s.push_str("##");
        s.push_str(m);
        s.push('\n');
    }
    s.push_str(columns.ok_or_else(|| invalid("gVCF has no #CHROM line".into()))?);
    s.push('\n');
    Ok(s)
}

/// Genotype a whole gVCF; writes the joint VCF (header and records) to `out`.
pub fn genotype_gvcf(
    gvcf: &Path,
    dbsnp: Option<&Path>,
    threads: usize,
    out: &mut dyn Write,
) -> io::Result<u64> {
    genotype_gvcf_text(&read_text(gvcf, threads)?, dbsnp, threads, out)
}

/// [`genotype_gvcf`] on gVCF text already in memory (`nomics run`).
pub fn genotype_gvcf_text(
    text: &str,
    dbsnp: Option<&Path>,
    threads: usize,
    out: &mut dyn Write,
) -> io::Result<u64> {
    let db = dbsnp.map(Dbsnp::open).transpose()?;
    let header: Vec<String> = text
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(str::to_string)
        .collect();
    let mut joint = Joint::new(&header);

    // The gVCF's per-contig runs, in file (dictionary) order.
    let mut runs: Vec<(&str, Vec<&str>)> = Vec::new();
    for line in text.lines().skip(header.len()) {
        if line.is_empty() {
            continue;
        }
        let contig = line.split('\t').next().unwrap_or_default();
        match runs.last_mut() {
            Some((c, lines)) if *c == contig => lines.push(line),
            _ => runs.push((contig, vec![line])),
        }
    }
    // GenotypeGVCFs reads the GenomicsDB workspace array by array, in sorted
    // array-name order (`<contig>$<start>$<end>`, one array per contig under
    // --merge-input-intervals). That order decides which QD jitter draw each
    // record takes; `bcftools sort` then restores dictionary order.
    let mut traversal: Vec<usize> = (0..runs.len()).collect();
    traversal.sort_by(|&a, &b| format!("{}$", runs[a].0).cmp(&format!("{}$", runs[b].0)));
    let mut called: Vec<Vec<Called>> = (0..runs.len()).map(|_| Vec::new()).collect();
    for &r in &traversal {
        for line in &runs[r].1 {
            if let Some(c) = joint.line(line)? {
                called[r].push(c);
            }
        }
    }
    if let Some(db) = &db {
        annotate_dbsnp(db, &mut called, threads)?;
    }

    out.write_all(joint_header(&header)?.as_bytes())?;
    let mut n = 0;
    let mut buf = String::with_capacity(1 << 20);
    for c in called.iter().flatten() {
        buf.push_str(&render(c));
        buf.push('\n');
        n += 1;
        if buf.len() >= 1 << 20 {
            out.write_all(buf.as_bytes())?;
            buf.clear();
        }
    }
    out.write_all(buf.as_bytes())?;
    Ok(n)
}

/// `VariantOverlapAnnotator.annotateRsID` + `annotateOverlaps`: the dbSNP
/// IDs whose events match the final alleles, and the `DB` flag. Independent
/// per record, so it runs after genotyping on `threads` threads, one tabix
/// query per window of nearby records.
fn annotate_dbsnp(db: &Dbsnp, called: &mut [Vec<Called>], threads: usize) -> io::Result<()> {
    const WINDOW: i64 = 1 << 16;
    let mut windows: Vec<&mut [Called]> = Vec::new();
    for run in called.iter_mut() {
        let mut rest: &mut [Called] = run;
        while !rest.is_empty() {
            let first = rest[0].pos;
            let len = rest.iter().take_while(|c| c.pos - first < WINDOW).count();
            let (head, tail) = rest.split_at_mut(len);
            windows.push(head);
            rest = tail;
        }
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let windows: Vec<std::sync::Mutex<&mut [Called]>> =
        windows.into_iter().map(std::sync::Mutex::new).collect();
    std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads.max(1))
            .map(|_| {
                scope.spawn(|| -> io::Result<()> {
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(window) = windows.get(i) else {
                            return Ok(());
                        };
                        let mut window = window.lock().expect("window lock");
                        let (first, last) = (window[0].pos, window[window.len() - 1].pos);
                        let slice = db.slice(window[0].contig.as_bytes(), first, last)?;
                        for c in window.iter_mut() {
                            let alleles: Vec<&[u8]> = c.alleles.iter().map(Vec::as_slice).collect();
                            if let Some(rs) = slice.id_for(c.pos, &alleles) {
                                if c.id == "." {
                                    c.id = rs;
                                } else if !c.id.contains(&rs) {
                                    c.id = format!("{};{rs}", c.id);
                                }
                                c.info.insert("DB".into(), Value::Flag);
                            }
                        }
                    }
                })
            })
            .collect();
        workers
            .into_iter()
            .try_for_each(|w| w.join().expect("dbSNP worker panicked"))
    })
}

fn read_text(path: &Path, threads: usize) -> io::Result<String> {
    let data = std::fs::read(path)?;
    let bytes = if data.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        nomics_io::bgzf::decode_parallel(&data, threads, &mut |b| {
            out.extend_from_slice(b);
            Ok(())
        })?;
        out
    } else {
        data
    };
    String::from_utf8(bytes).map_err(|_| invalid("gVCF is not UTF-8".into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trimmed(alleles: &[&str]) -> Vec<String> {
        let mut a: Vec<Vec<u8>> = alleles.iter().map(|s| s.as_bytes().to_vec()).collect();
        reverse_trim(&mut a);
        a.into_iter().map(|x| String::from_utf8(x).unwrap()).collect()
    }

    #[test]
    fn reverse_trim_keeps_one_base_and_skips_star() {
        assert_eq!(trimmed(&["GAA", "GA"]), ["GA", "G"]);
        assert_eq!(trimmed(&["CA", "C"]), ["CA", "C"]);
        assert_eq!(trimmed(&["AT", "*", "GT"]), ["A", "*", "G"]);
        assert_eq!(trimmed(&["A", "G"]), ["A", "G"]);
    }

    #[test]
    fn header_follows_htsjdk_order() {
        let gvcf: Vec<String> = [
            "##fileformat=VCFv4.2",
            "##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)",
            "##FILTER=<ID=LowQual,Description=\"Low quality\">",
            "##GATKCommandLine=<ID=HaplotypeCaller,CommandLine=\"x\">",
            "##contig=<ID=chr2,length=2>",
            "##contig=<ID=chr10,length=10>",
            "##source=HaplotypeCaller",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let h = joint_header(&gvcf).unwrap();
        let keys: Vec<&str> = h.lines().map(|l| l.split(['=', '\t']).next().unwrap()).collect();
        assert_eq!(keys.first(), Some(&"##fileformat"));
        assert!(!h.contains("GVCFBlock") && !h.contains("HaplotypeCaller"));
        assert!(h.find("ID=chr2,").unwrap() < h.find("ID=chr10,").unwrap());
        assert!(h.find("##INFO=<ID=SOR").unwrap() < h.find("##contig").unwrap());
        assert!(h.find("##contig").unwrap() < h.find("##source=GenomicsDBImport").unwrap());
        assert_eq!(h.matches("ID=LowQual").count(), 1);
    }
}
