//! Recalibration tables: accumulation (§4), finalisation (§5), the GATKReport
//! v1.1 text (§7) and the ApplyBQSR quality model (§8.3).

use crate::reads::{context_from_key, cycle_from_key, key_from_context};
use crate::stats::{empirical_quality, qual_to_error_prob, quantize, round4};

pub const N_QUAL: usize = 94;
pub const N_CONTEXT: usize = 1012;
pub const N_CYCLE: usize = 1002;

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Counts {
    pub obs: u64,
    pub err: u64,
}

impl Counts {
    fn add(&mut self, o: Counts) {
        self.obs += o.obs;
        self.err += o.err;
    }
}

/// `M`-event tables for every read group, indexed `[rg][q]`, `[rg][q][key]`.
#[derive(Clone)]
pub struct Tables {
    pub read_groups: Vec<String>,
    pub qual: Vec<Counts>,
    pub context: Vec<Counts>,
    pub cycle: Vec<Counts>,
    /// Per read group `EstimatedQReported` as GatherBQSRReports combines it
    /// across scatter groups; `None` for a single BaseRecalibrator run.
    pub gathered_estimates: Option<Vec<Option<f64>>>,
}

impl Tables {
    pub fn new(read_groups: Vec<String>) -> Tables {
        let n = read_groups.len();
        Tables {
            read_groups,
            qual: vec![Counts::default(); n * N_QUAL],
            context: vec![Counts::default(); n * N_QUAL * N_CONTEXT],
            cycle: vec![Counts::default(); n * N_QUAL * N_CYCLE],
            gathered_estimates: None,
        }
    }

    /// Count one base `times` times: once per scatter group whose
    /// BaseRecalibrator run visits the read.
    #[inline]
    pub fn add_base(
        &mut self,
        rg: usize,
        q: usize,
        context: i32,
        cycle: i32,
        error: bool,
        times: u64,
    ) {
        let e = error as u64 * times;
        let qi = rg * N_QUAL + q;
        self.qual[qi].obs += times;
        self.qual[qi].err += e;
        if context >= 0 {
            let c = &mut self.context[qi * N_CONTEXT + context as usize];
            c.obs += times;
            c.err += e;
        }
        let c = &mut self.cycle[qi * N_CYCLE + cycle as usize];
        c.obs += times;
        c.err += e;
    }

    /// Counts are exact integers, so merging in any order is exact (§4.2).
    pub fn merge(&mut self, o: &Tables) {
        for (a, b) in self.qual.iter_mut().zip(&o.qual) {
            a.add(*b);
        }
        for (a, b) in self.context.iter_mut().zip(&o.context) {
            a.add(*b);
        }
        for (a, b) in self.cycle.iter_mut().zip(&o.cycle) {
            a.add(*b);
        }
    }

    /// Per read group: observations, errors and the rounded `EstimatedQReported`
    /// from the iterative combine over qualities ascending (§5.1-5.2).
    pub fn read_group_summary(&self, rg: usize) -> Option<(Counts, f64)> {
        let (c, estimated) = qual_summary(&self.qual[rg * N_QUAL..(rg + 1) * N_QUAL])?;
        match &self.gathered_estimates {
            Some(g) => g[rg].map(|e| (c, e)),
            None => Some((c, estimated)),
        }
    }
}

/// One read group's quality rows collapsed into its read-group row:
/// observations, errors and the rounded `EstimatedQReported` (§5.1-5.2).
pub fn qual_summary(qual: &[Counts]) -> Option<(Counts, f64)> {
    let mut acc: Option<(Counts, f64)> = None;
    for (q, &c) in qual.iter().enumerate() {
        if c.obs == 0 {
            continue;
        }
        acc = Some(match acc {
            None => (c, q as f64),
            Some((mut n, r)) => {
                let expected = n.obs as f64 * qual_to_error_prob(r)
                    + c.obs as f64 * qual_to_error_prob(q as f64);
                n.add(c);
                (n, -10.0 * (expected / n.obs as f64).log10())
            }
        });
    }
    acc.map(|(c, r)| (c, round4(r)))
}

/// GatherBQSRReports' read-group estimates. Each scatter group's report is
/// parsed back from text (so its estimate is the rounded one), then the
/// reports are combined in `order` with `RecalDatum.combine`: expected errors
/// add, and the estimate is recomputed unrounded until it is printed.
/// `group_qual` is `[group][rg][q]`.
pub fn gathered_estimates(
    group_qual: &[Counts],
    n_rg: usize,
    order: &[usize],
) -> Vec<Option<f64>> {
    (0..n_rg)
        .map(|rg| {
            let mut acc: Option<(u64, f64)> = None;
            for &g in order {
                let at = (g * n_rg + rg) * N_QUAL;
                let Some((c, e)) = qual_summary(&group_qual[at..at + N_QUAL]) else {
                    continue;
                };
                acc = Some(match acc {
                    None => (c.obs, e),
                    Some((n, r)) => {
                        let expected =
                            n as f64 * qual_to_error_prob(r) + c.obs as f64 * qual_to_error_prob(e);
                        let n = n + c.obs;
                        (n, -10.0 * (expected / n as f64).log10())
                    }
                });
            }
            acc.map(|(_, r)| round4(r))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// GATKReport writing (§7)

#[derive(Clone, Debug, PartialEq)]
enum Cell {
    Str(String),
    Int(i64),
    Dbl(f64, usize),
}

impl Cell {
    fn text(&self) -> String {
        match self {
            Cell::Str(s) => s.clone(),
            Cell::Int(v) => v.to_string(),
            Cell::Dbl(v, p) => format!("{v:.p$}"),
        }
    }

    fn cmp(&self, o: &Cell) -> std::cmp::Ordering {
        match (self, o) {
            (Cell::Int(a), Cell::Int(b)) => a.cmp(b),
            (Cell::Dbl(a, _), Cell::Dbl(b, _)) => a.total_cmp(b),
            _ => self.text().cmp(&o.text()),
        }
    }
}

/// commons-lang3 `NumberUtils.isCreatable` for the strings these tables hold.
fn looks_numeric(s: &str) -> bool {
    let s = s.strip_prefix(['-', '+']).unwrap_or(s);
    let (mantissa, exp) = match s.find(['e', 'E']) {
        Some(i) => (&s[..i], Some(&s[i + 1..])),
        None => (s, None),
    };
    let digits_ok = !mantissa.is_empty()
        && mantissa.chars().filter(|&c| c == '.').count() <= 1
        && mantissa.chars().all(|c| c.is_ascii_digit() || c == '.')
        && mantissa.chars().any(|c| c.is_ascii_digit());
    let exp_ok = exp.is_none_or(|e| {
        let e = e.strip_prefix(['-', '+']).unwrap_or(e);
        !e.is_empty() && e.chars().all(|c| c.is_ascii_digit())
    });
    digits_ok && exp_ok
}

struct Table {
    name: &'static str,
    description: &'static str,
    columns: Vec<(&'static str, &'static str)>,
    rows: Vec<Vec<Cell>>,
}

impl Table {
    fn write(mut self, out: &mut String) {
        self.rows.sort_by(|a, b| {
            a.iter()
                .zip(b)
                .map(|(x, y)| x.cmp(y))
                .find(|o| o.is_ne())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let ncol = self.columns.len();
        let texts: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|r| r.iter().map(Cell::text).collect())
            .collect();
        let mut widths: Vec<usize> = self.columns.iter().map(|(n, _)| n.len()).collect();
        let mut right: Vec<bool> = vec![true; ncol];
        for row in &texts {
            for (k, t) in row.iter().enumerate() {
                widths[k] = widths[k].max(t.len());
                if !matches!(t.as_str(), "null" | "NA" | "Infinity" | "-Infinity" | "NaN")
                    && !looks_numeric(t)
                {
                    right[k] = false;
                }
            }
        }
        out.push_str(&format!("#:GATKTable:{}:{}", ncol, self.rows.len()));
        for (_, f) in &self.columns {
            out.push_str(&format!(":{f}"));
        }
        out.push_str(":;\n");
        out.push_str(&format!("#:GATKTable:{}:{}\n", self.name, self.description));
        let header: Vec<String> = self
            .columns
            .iter()
            .zip(&widths)
            .map(|((n, _), &w)| format!("{n:<w$}"))
            .collect();
        out.push_str(&header.join("  "));
        out.push('\n');
        for row in &texts {
            let cells: Vec<String> = row
                .iter()
                .enumerate()
                .map(|(k, t)| {
                    if right[k] {
                        format!("{t:>w$}", w = widths[k])
                    } else {
                        format!("{t:<w$}", w = widths[k])
                    }
                })
                .collect();
            out.push_str(&cells.join("  "));
            out.push('\n');
        }
        out.push('\n');
    }
}

/// The finished recalibration report: the five tables plus what ApplyBQSR needs.
pub struct Report {
    pub text: String,
    pub quantized: Vec<u8>,
}

/// `Arguments` table values (fixed for this pipeline's command line).
const REPORT_ARGUMENTS: [(&str, &str); 17] = [
    ("binary_tag_name", "null"),
    (
        "covariate",
        "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate",
    ),
    ("default_platform", "null"),
    ("deletions_default_quality", "45"),
    ("force_platform", "null"),
    ("indels_context_size", "3"),
    ("insertions_default_quality", "45"),
    ("low_quality_tail", "2"),
    ("maximum_cycle_value", "500"),
    ("mismatches_context_size", "2"),
    ("mismatches_default_quality", "-1"),
    ("no_standard_covs", "false"),
    ("quantizing_levels", "16"),
    ("recalibration_report", "null"),
    ("run_without_dbsnp", "false"),
    ("solid_nocall_strategy", "THROW_EXCEPTION"),
    ("solid_recal_mode", "SET_Q_ZERO"),
];

/// Build the GATKReport text and quantization map from accumulated tables.
pub fn build_report(t: &Tables) -> Report {
    let n_rg = t.read_groups.len();
    let mut out = String::from("#:GATKReport.v1.1:5\n");

    // Arguments (fixed for this pipeline's command line).
    let args = REPORT_ARGUMENTS;
    Table {
        name: "Arguments",
        description: "Recalibration argument collection values used in this run",
        columns: vec![("Argument", "%s"), ("Value", "%s")],
        rows: args
            .iter()
            .map(|(k, v)| vec![Cell::Str((*k).into()), Cell::Str((*v).into())])
            .collect(),
    }
    .write(&mut out);

    // Quantized: histogram of quality-table observations by empirical quality (§6).
    let mut hist = vec![0u64; N_QUAL];
    for rg in 0..n_rg {
        for q in 0..N_QUAL {
            let c = t.qual[rg * N_QUAL + q];
            if c.obs > 0 {
                hist[empirical_quality(c.obs, c.err as f64, q as f64) as usize] += c.obs;
            }
        }
    }
    let quantized = quantize(&hist, 16, 6);
    Table {
        name: "Quantized",
        description: "Quality quantization map",
        columns: vec![
            ("QualityScore", "%d"),
            ("Count", "%d"),
            ("QuantizedScore", "%d"),
        ],
        rows: (0..N_QUAL)
            .map(|q| {
                vec![
                    Cell::Int(q as i64),
                    Cell::Int(hist[q] as i64),
                    Cell::Int(quantized[q] as i64),
                ]
            })
            .collect(),
    }
    .write(&mut out);

    let num = |c: Counts| (Cell::Int(c.obs as i64), Cell::Dbl(c.err as f64, 2));
    let mut rt0 = Vec::new();
    let mut rt1 = Vec::new();
    let mut rt2 = Vec::new();
    for rg in 0..n_rg {
        let name = Cell::Str(t.read_groups[rg].clone());
        if let Some((c, estimated)) = t.read_group_summary(rg) {
            let emp = empirical_quality(c.obs, c.err as f64, estimated) as f64;
            let (o, e) = num(c);
            rt0.push(vec![
                name.clone(),
                Cell::Str("M".into()),
                Cell::Dbl(emp, 4),
                Cell::Dbl(estimated, 4),
                o,
                e,
            ]);
        }
        for q in 0..N_QUAL {
            let qi = rg * N_QUAL + q;
            let c = t.qual[qi];
            if c.obs == 0 {
                continue;
            }
            let emp =
                |c: Counts| Cell::Dbl(empirical_quality(c.obs, c.err as f64, q as f64) as f64, 4);
            let (o, e) = num(c);
            rt1.push(vec![
                name.clone(),
                Cell::Int(q as i64),
                Cell::Str("M".into()),
                emp(c),
                o,
                e,
            ]);
            for key in 0..N_CONTEXT {
                let c = t.context[qi * N_CONTEXT + key];
                if c.obs > 0 {
                    let (o, e) = num(c);
                    rt2.push(vec![
                        name.clone(),
                        Cell::Int(q as i64),
                        Cell::Str(context_from_key(key)),
                        Cell::Str("Context".into()),
                        Cell::Str("M".into()),
                        emp(c),
                        o,
                        e,
                    ]);
                }
            }
            for key in 0..N_CYCLE {
                let c = t.cycle[qi * N_CYCLE + key];
                if c.obs > 0 {
                    let (o, e) = num(c);
                    rt2.push(vec![
                        name.clone(),
                        Cell::Int(q as i64),
                        Cell::Str(cycle_from_key(key).to_string()),
                        Cell::Str("Cycle".into()),
                        Cell::Str("M".into()),
                        emp(c),
                        o,
                        e,
                    ]);
                }
            }
        }
    }
    let common = [("ReadGroup", "%s")];
    Table {
        name: "RecalTable0",
        description: "",
        columns: [
            &common[..],
            &[
                ("EventType", "%s"),
                ("EmpiricalQuality", "%.4f"),
                ("EstimatedQReported", "%.4f"),
                ("Observations", "%d"),
                ("Errors", "%.2f"),
            ],
        ]
        .concat(),
        rows: rt0,
    }
    .write(&mut out);
    Table {
        name: "RecalTable1",
        description: "",
        columns: [
            &common[..],
            &[
                ("QualityScore", "%d"),
                ("EventType", "%s"),
                ("EmpiricalQuality", "%.4f"),
                ("Observations", "%d"),
                ("Errors", "%.2f"),
            ],
        ]
        .concat(),
        rows: rt1,
    }
    .write(&mut out);
    Table {
        name: "RecalTable2",
        description: "",
        columns: [
            &common[..],
            &[
                ("QualityScore", "%d"),
                ("CovariateValue", "%s"),
                ("CovariateName", "%s"),
                ("EventType", "%s"),
                ("EmpiricalQuality", "%.4f"),
                ("Observations", "%d"),
                ("Errors", "%.2f"),
            ],
        ]
        .concat(),
        rows: rt2,
    }
    .write(&mut out);
    Report {
        text: out,
        quantized,
    }
}

// ---------------------------------------------------------------------------
// ApplyBQSR model (§8.3)

/// Precomputed recalibrated qualities: for each read group and quality,
/// `Qq` plus the per-context and per-cycle deltas, with hierarchical priors.
pub struct Recalibration {
    /// `[rg][q]` -> Qq (or Qrg if the quality has no row); None if the read
    /// group has no data.
    qq: Vec<Option<i32>>,
    ctx_delta: Vec<i32>,
    cyc_delta: Vec<i32>,
    has_ctx: Vec<bool>,
    has_cyc: Vec<bool>,
    q93: u8,
}

impl Recalibration {
    pub fn new(t: &Tables, quantized: &[u8]) -> Recalibration {
        let n_rg = t.read_groups.len();
        let mut r = Recalibration {
            qq: vec![None; n_rg * N_QUAL],
            ctx_delta: vec![0; n_rg * N_QUAL * N_CONTEXT],
            cyc_delta: vec![0; n_rg * N_QUAL * N_CYCLE],
            has_ctx: vec![false; n_rg * N_QUAL * N_CONTEXT],
            has_cyc: vec![false; n_rg * N_QUAL * N_CYCLE],
            q93: quantized[93],
        };
        for rg in 0..n_rg {
            let Some((rg_counts, estimated)) = t.read_group_summary(rg) else {
                continue;
            };
            let q_rg = empirical_quality(rg_counts.obs, rg_counts.err as f64, estimated);
            for q in 0..N_QUAL {
                let qi = rg * N_QUAL + q;
                let c = t.qual[qi];
                let q_q = if c.obs > 0 {
                    empirical_quality(c.obs, c.err as f64, q_rg as f64)
                } else {
                    q_rg
                };
                r.qq[qi] = Some(q_q);
                for key in 0..N_CONTEXT {
                    let c = t.context[qi * N_CONTEXT + key];
                    if c.obs > 0 {
                        r.has_ctx[qi * N_CONTEXT + key] = true;
                        r.ctx_delta[qi * N_CONTEXT + key] =
                            empirical_quality(c.obs, c.err as f64, q_q as f64) - q_q;
                    }
                }
                for key in 0..N_CYCLE {
                    let c = t.cycle[qi * N_CYCLE + key];
                    if c.obs > 0 {
                        r.has_cyc[qi * N_CYCLE + key] = true;
                        r.cyc_delta[qi * N_CYCLE + key] =
                            empirical_quality(c.obs, c.err as f64, q_q as f64) - q_q;
                    }
                }
            }
        }
        r
    }

    /// New quality for one base (QUAL >= 6), `None` if the read group has no data.
    #[inline]
    pub fn recalibrate(&self, rg: usize, q: usize, context: i32, cycle: i32) -> Option<u8> {
        let qi = rg * N_QUAL + q;
        let q_q = self.qq[qi]?;
        let mut raw = q_q;
        if context >= 0 && self.has_ctx[qi * N_CONTEXT + context as usize] {
            raw += self.ctx_delta[qi * N_CONTEXT + context as usize];
        }
        if self.has_cyc[qi * N_CYCLE + cycle as usize] {
            raw += self.cyc_delta[qi * N_CYCLE + cycle as usize];
        }
        let bounded = raw.clamp(1, 93) as u8;
        Some(if bounded == 93 { self.q93 } else { bounded })
    }
}

// ---------------------------------------------------------------------------
// Reading a report back (ApplyBQSR, §8.2)

/// One parsed GATKReport table: name, column names, rows of cells.
struct ParsedTable {
    name: String,
    columns: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn parse_tables(text: &str) -> Result<Vec<ParsedTable>, String> {
    let mut lines = text.lines();
    match lines.next() {
        Some(l) if l.starts_with("#:GATKReport.v1.1:") => {}
        other => return Err(format!("not a GATKReport v1.1: {other:?}")),
    }
    let mut out = Vec::new();
    while let Some(line) = lines.next() {
        if line.is_empty() {
            continue;
        }
        let Some(format) = line.strip_prefix("#:GATKTable:") else {
            return Err(format!("unexpected line {line:?}"));
        };
        let nrows: usize = format
            .split(':')
            .nth(1)
            .and_then(|n| n.parse().ok())
            .ok_or_else(|| format!("bad table format line {line:?}"))?;
        let name_line = lines.next().ok_or("truncated report")?;
        let name = name_line
            .strip_prefix("#:GATKTable:")
            .and_then(|r| r.split(':').next())
            .ok_or_else(|| format!("bad table name line {name_line:?}"))?
            .to_owned();
        let columns: Vec<String> = lines
            .next()
            .ok_or("truncated report")?
            .split_whitespace()
            .map(str::to_owned)
            .collect();
        let mut rows = Vec::with_capacity(nrows);
        for _ in 0..nrows {
            let row: Vec<String> = lines
                .next()
                .ok_or("truncated table")?
                .split_whitespace()
                .map(str::to_owned)
                .collect();
            if row.len() != columns.len() {
                return Err(format!(
                    "{name}: row with {} of {} cells",
                    row.len(),
                    columns.len()
                ));
            }
            rows.push(row);
        }
        out.push(ParsedTable {
            name,
            columns,
            rows,
        });
    }
    Ok(out)
}

impl ParsedTable {
    fn col(&self, name: &str) -> Result<usize, String> {
        self.columns
            .iter()
            .position(|c| c == name)
            .ok_or_else(|| format!("{}: no column {name}", self.name))
    }
}

/// Check that a report has the five standard tables and that its
/// `Arguments` are exactly the ones this port writes (and so implements).
/// ApplyBQSR re-reads them from the table, so any other value is outside
/// the supported set.
pub fn check_report_arguments(text: &str) -> Result<(), String> {
    let tables = parse_tables(text)?;
    let names: Vec<&str> = tables.iter().map(|t| t.name.as_str()).collect();
    if names != ["Arguments", "Quantized", "RecalTable0", "RecalTable1", "RecalTable2"] {
        return Err(format!("report tables {names:?}"));
    }
    let got: Vec<(&str, &str)> = tables[0]
        .rows
        .iter()
        .map(|r| (r[0].as_str(), r[1].as_str()))
        .collect();
    if got != REPORT_ARGUMENTS {
        let diff = got
            .iter()
            .zip(REPORT_ARGUMENTS.iter())
            .find(|(a, b)| a != b)
            .map_or_else(
                || format!("{} argument rows", got.len()),
                |(a, _)| format!("{} {}", a.0, a.1),
            );
        return Err(format!("report argument {diff}"));
    }
    Ok(())
}

fn int(s: &str, what: &str) -> Result<u64, String> {
    s.parse().map_err(|_| format!("bad {what} {s:?}"))
}

/// `Errors` is printed `%.2f` but counts whole mismatches for `M` events.
fn errors(s: &str) -> Result<u64, String> {
    let v: f64 = s.parse().map_err(|_| format!("bad Errors {s:?}"))?;
    if v < 0.0 || v.fract() != 0.0 || v > (1u64 << 53) as f64 {
        return Err(format!("non-integral Errors {s:?}"));
    }
    Ok(v as u64)
}

/// Parse a BaseRecalibrator report back into tables for ApplyBQSR, as
/// `RecalibrationReport` does (§8.2): counts from every `M` row, the read
/// group's parsed `EstimatedQReported` as its prior, and the `Quantized`
/// column. `read_groups` are the BAM header's identifiers (PU, else ID).
pub fn parse_report(text: &str, read_groups: Vec<String>) -> Result<(Tables, Vec<u8>), String> {
    check_report_arguments(text)?;
    let tables = parse_tables(text)?;
    let n_rg = read_groups.len();
    let rg_of = |name: &str| -> Result<usize, String> {
        read_groups
            .iter()
            .position(|g| g == name)
            .ok_or_else(|| format!("report read group {name:?} is not in the BAM header"))
    };
    let mut t = Tables::new(read_groups.clone());
    // Quantized
    let qt = &tables[1];
    let (qc, qq) = (qt.col("QualityScore")?, qt.col("QuantizedScore")?);
    if qt.rows.len() != N_QUAL {
        return Err(format!("Quantized has {} rows", qt.rows.len()));
    }
    let mut quantized = vec![0u8; N_QUAL];
    for r in &qt.rows {
        let q = int(&r[qc], "QualityScore")? as usize;
        let v = int(&r[qq], "QuantizedScore")?;
        if q >= N_QUAL || v > 93 {
            return Err(format!("Quantized row {r:?}"));
        }
        quantized[q] = v as u8;
    }
    // RecalTable1: quality rows.
    let t1 = &tables[3];
    let (rg1, q1, e1, o1, x1) = (
        t1.col("ReadGroup")?,
        t1.col("QualityScore")?,
        t1.col("EventType")?,
        t1.col("Observations")?,
        t1.col("Errors")?,
    );
    for r in &t1.rows {
        if r[e1] != "M" {
            continue;
        }
        let rg = rg_of(&r[rg1])?;
        let q = int(&r[q1], "QualityScore")? as usize;
        if q >= N_QUAL {
            return Err(format!("RecalTable1 quality {q}"));
        }
        let c = &mut t.qual[rg * N_QUAL + q];
        c.obs = int(&r[o1], "Observations")?;
        c.err = errors(&r[x1])?;
    }
    // RecalTable2: context and cycle rows.
    let t2 = &tables[4];
    let (rg2, q2, v2, n2, e2, o2, x2) = (
        t2.col("ReadGroup")?,
        t2.col("QualityScore")?,
        t2.col("CovariateValue")?,
        t2.col("CovariateName")?,
        t2.col("EventType")?,
        t2.col("Observations")?,
        t2.col("Errors")?,
    );
    for r in &t2.rows {
        if r[e2] != "M" {
            continue;
        }
        let rg = rg_of(&r[rg2])?;
        let q = int(&r[q2], "QualityScore")? as usize;
        if q >= N_QUAL {
            return Err(format!("RecalTable2 quality {q}"));
        }
        let qi = rg * N_QUAL + q;
        let counts = Counts {
            obs: int(&r[o2], "Observations")?,
            err: errors(&r[x2])?,
        };
        match r[n2].as_str() {
            "Context" => {
                let key = key_from_context(r[v2].as_bytes());
                if key < 0 || key as usize >= N_CONTEXT || context_from_key(key as usize) != r[v2]
                {
                    return Err(format!("context {:?}", r[v2]));
                }
                t.context[qi * N_CONTEXT + key as usize] = counts;
            }
            "Cycle" => {
                let c: i64 = r[v2].parse().map_err(|_| format!("cycle {:?}", r[v2]))?;
                let key = ((c.unsigned_abs() << 1) | (c < 0) as u64) as usize;
                if key >= N_CYCLE || cycle_from_key(key) != c {
                    return Err(format!("cycle {c}"));
                }
                t.cycle[qi * N_CYCLE + key] = counts;
            }
            other => return Err(format!("covariate {other:?}")),
        }
    }
    // RecalTable0: the read-group priors, checked against the quality rows.
    let t0 = &tables[2];
    let (rg0, e0, est0, o0, x0) = (
        t0.col("ReadGroup")?,
        t0.col("EventType")?,
        t0.col("EstimatedQReported")?,
        t0.col("Observations")?,
        t0.col("Errors")?,
    );
    let mut estimates = vec![None; n_rg];
    for r in &t0.rows {
        if r[e0] != "M" {
            continue;
        }
        let rg = rg_of(&r[rg0])?;
        let est: f64 = r[est0]
            .parse()
            .map_err(|_| format!("EstimatedQReported {:?}", r[est0]))?;
        let summed = qual_summary(&t.qual[rg * N_QUAL..(rg + 1) * N_QUAL]).map(|(c, _)| c);
        let row = Counts {
            obs: int(&r[o0], "Observations")?,
            err: errors(&r[x0])?,
        };
        if summed != Some(row) {
            return Err(format!(
                "RecalTable0 {} counts {row:?} differ from its quality rows {summed:?}",
                r[rg0]
            ));
        }
        estimates[rg] = Some(est);
    }
    t.gathered_estimates = Some(estimates);
    Ok((t, quantized))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_detection() {
        for s in ["37", "-151", "38.0000", "1e5", "0.5"] {
            assert!(looks_numeric(s), "{s}");
        }
        for s in ["L1", "M", "AA", "Context", "", "-", "1.2.3", "null"] {
            assert!(!looks_numeric(s), "{s}");
        }
    }

    #[test]
    fn read_group_summary_matches_tier0_golden() {
        let mut t = Tables::new(vec!["L1".into()]);
        t.qual[11] = Counts {
            obs: 1_320_847,
            err: 100_924,
        };
        t.qual[25] = Counts {
            obs: 2_851_349,
            err: 5_427,
        };
        t.qual[37] = Counts {
            obs: 87_437_550,
            err: 14_884,
        };
        let (c, est) = t.read_group_summary(0).unwrap();
        assert_eq!((c.obs, c.err), (91_609_746, 121_235));
        assert_eq!(format!("{est:.4}"), "28.4341");
    }

    /// ApplyBQSR's view of a report (§8.2) equals BaseRecalibrator's tables.
    #[test]
    fn parsed_report_round_trips_tables() {
        let mut t = Tables::new(vec!["L1".into(), "L2".into()]);
        let ctx = key_from_context(b"AC");
        for (rg, q, cyc, err, times) in [
            (0, 11, 2, true, 3),
            (0, 11, 2, false, 40),
            (0, 37, 3, false, 1000),
            (0, 37, 3, true, 7),
            (1, 25, 1001, false, 12),
        ] {
            t.add_base(rg, q, ctx, cyc, err, times);
        }
        t.add_base(1, 25, -1, 5, true, 2);
        let report = build_report(&t);
        let (p, quantized) = parse_report(&report.text, t.read_groups.clone()).unwrap();
        assert_eq!(p.qual, t.qual);
        assert_eq!(p.context, t.context);
        assert_eq!(p.cycle, t.cycle);
        assert_eq!(quantized, report.quantized);
        for rg in 0..2 {
            assert_eq!(p.read_group_summary(rg), t.read_group_summary(rg));
        }
        // A read group missing from the header, or other arguments, are refused.
        assert!(parse_report(&report.text, vec!["L1".into()]).is_err());
        let other = report.text.replace(
            "maximum_cycle_value         500",
            "maximum_cycle_value         501",
        );
        assert!(check_report_arguments(&other).is_err());
    }

    #[test]
    fn table_layout_pads_and_aligns() {
        let mut s = String::new();
        Table {
            name: "T",
            description: "",
            columns: vec![("ReadGroup", "%s"), ("Errors", "%.2f")],
            rows: vec![vec![Cell::Str("L1".into()), Cell::Dbl(123456.0, 2)]],
        }
        .write(&mut s);
        assert_eq!(
            s,
            "#:GATKTable:2:1:%s:%.2f:;\n#:GATKTable:T:\nReadGroup  Errors   \nL1         123456.00\n\n"
        );
    }
}
