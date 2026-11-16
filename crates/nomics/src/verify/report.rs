//! `verify.json` (schema `nomics-verify/1`, docs/VERIFY.md) and `verify.md`.
//! JSON is written by hand: the workspace takes no serialization crates, and
//! a fixed key order is what makes `report_sha256` reproducible.

use std::fmt::Write as _;

use super::gate::{ArtifactResult, Status};

pub const SCHEMA: &str = "nomics-verify/1";

/// A JSON value whose objects keep insertion order.
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Int(i64),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn str(s: impl Into<String>) -> Json {
        Json::Str(s.into())
    }

    pub fn obj(pairs: Vec<(&str, Json)>) -> Json {
        Json::Obj(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
    }

    pub fn opt_str(s: Option<&str>) -> Json {
        s.map_or(Json::Null, Json::str)
    }

    /// Pretty-printed with two-space indents and a final newline.
    pub fn render(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out.push('\n');
        out
    }

    fn write(&self, out: &mut String, depth: usize) {
        let pad = |out: &mut String, d: usize| {
            for _ in 0..d {
                out.push_str("  ");
            }
        };
        match self {
            Json::Null => out.push_str("null"),
            Json::Int(i) => {
                let _ = write!(out, "{i}");
            }
            // Timings: three decimals is a millisecond, and keeps the text stable.
            Json::Num(x) if x.is_finite() => {
                let _ = write!(out, "{x:.3}");
            }
            Json::Num(_) => out.push_str("null"),
            Json::Str(s) => escape(s, out),
            Json::Arr(items) if items.is_empty() => out.push_str("[]"),
            Json::Arr(items) => {
                out.push_str("[\n");
                for (i, item) in items.iter().enumerate() {
                    pad(out, depth + 1);
                    item.write(out, depth + 1);
                    out.push_str(if i + 1 < items.len() { ",\n" } else { "\n" });
                }
                pad(out, depth);
                out.push(']');
            }
            Json::Obj(pairs) if pairs.is_empty() => out.push_str("{}"),
            Json::Obj(pairs) => {
                out.push_str("{\n");
                for (i, (k, v)) in pairs.iter().enumerate() {
                    pad(out, depth + 1);
                    escape(k, out);
                    out.push_str(": ");
                    v.write(out, depth + 1);
                    out.push_str(if i + 1 < pairs.len() { ",\n" } else { "\n" });
                }
                pad(out, depth);
                out.push('}');
            }
        }
    }
}

fn escape(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    Pass,
    Fail,
    Error,
}

impl Verdict {
    pub fn name(self) -> &'static str {
        match self {
            Verdict::Pass => "PASS",
            Verdict::Fail => "FAIL",
            Verdict::Error => "ERROR",
        }
    }

    pub fn exit_code(self) -> u8 {
        match self {
            Verdict::Pass => 0,
            Verdict::Fail => 1,
            Verdict::Error => 2,
        }
    }
}

/// The verdict from gate results. FAIL needs both sides to have run; anything
/// that prevents a real comparison is ERROR, and PASS needs every artifact
/// present, non-empty and identical, with enough VCF records to mean something.
pub fn verdict(results: &[ArtifactResult], min_vcf_records: u64) -> (Verdict, String) {
    let count = |s: Status| results.iter().filter(|r| r.status == s).count();
    let total = results.len();
    if total == 0 {
        return (Verdict::Error, "no artifacts were gated".to_owned());
    }
    let empty: Vec<&str> = results
        .iter()
        .filter(|r| r.status == Status::Empty)
        .map(|r| r.name.as_str())
        .collect();
    if !empty.is_empty() {
        return (
            Verdict::Error,
            format!("empty output, nothing to compare: {}", empty.join(", ")),
        );
    }
    let missing_ref: Vec<&str> = results
        .iter()
        .filter(|r| {
            (r.status == Status::Missing && r.reference.is_none())
                || (r.status == Status::Unreadable
                    && r.first_difference.as_ref().is_some_and(|d| d.detail.starts_with("reference")))
        })
        .map(|r| r.name.as_str())
        .collect();
    if !missing_ref.is_empty() {
        return (
            Verdict::Error,
            format!("the reference run did not produce a readable: {}", missing_ref.join(", ")),
        );
    }
    let (fail, missing) = (count(Status::Fail), count(Status::Missing) + count(Status::Unreadable));
    if fail + missing > 0 {
        let names: Vec<&str> = results
            .iter()
            .filter(|r| matches!(r.status, Status::Fail | Status::Missing | Status::Unreadable))
            .map(|r| r.name.as_str())
            .collect();
        return (
            Verdict::Fail,
            format!("{} of {total} artifacts differ or are missing: {}", fail + missing, names.join(", ")),
        );
    }
    let vcf: u64 = results
        .iter()
        .filter(|r| r.comparator == super::descriptor::Comparator::Vcf)
        .filter_map(|r| r.reference.as_ref().map(|d| d.count))
        .sum();
    if vcf < min_vcf_records {
        return (
            Verdict::Error,
            format!(
                "all {total} artifacts identical, but the reference VCF has only {vcf} records (fewer than {min_vcf_records}); take a larger slice"
            ),
        );
    }
    (
        Verdict::Pass,
        format!("{total}/{total} artifacts identical; {} VCF records", thousands(vcf)),
    )
}

pub fn thousands(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

pub fn artifact_json(r: &ArtifactResult) -> Json {
    let side = |d: &Option<nomics_gate::ArtifactDigest>| match d {
        Some(d) => Json::obj(vec![("sha256", Json::str(&d.sha256)), ("count", Json::Int(d.count as i64))]),
        None => Json::Null,
    };
    Json::obj(vec![
        ("name", Json::str(&r.name)),
        ("comparator", Json::str(r.comparator.name())),
        ("normalization", Json::str(&r.normalization)),
        ("reference", side(&r.reference)),
        ("nomics", side(&r.nomics)),
        ("status", Json::str(r.status.name())),
        (
            "first_difference",
            match &r.first_difference {
                Some(d) => Json::obj(vec![
                    ("index", d.index.map_or(Json::Null, |i| Json::Int(i as i64))),
                    ("detail", Json::str(&d.detail)),
                ]),
                None => Json::Null,
            },
        ),
    ])
}

/// Set `report_sha256` to the SHA-256 of the report rendered with that
/// field empty, and return the final text. Anyone can recompute it by
/// blanking the field and hashing the result (docs/VERIFY.md).
pub fn seal(mut report: Json) -> (String, String) {
    let Json::Obj(pairs) = &mut report else {
        panic!("report is an object");
    };
    let at = match pairs.iter().position(|(k, _)| k == "report_sha256") {
        Some(i) => i,
        None => {
            pairs.push(("report_sha256".to_owned(), Json::str("")));
            pairs.len() - 1
        }
    };
    pairs[at].1 = Json::str("");
    let hash = super::sha256_hex(report.render().as_bytes());
    if let Json::Obj(pairs) = &mut report {
        pairs[at].1 = Json::str(&hash);
    }
    (report.render(), hash)
}

/// Check a sealed report's text: blank `report_sha256` and rehash.
pub fn check_seal(text: &str) -> Option<bool> {
    let key = "\"report_sha256\": \"";
    let start = text.rfind(key)? + key.len();
    let end = start + text[start..].find('"')?;
    let claimed = &text[start..end];
    let blank = format!("{}{}", &text[..start], &text[end..]);
    Some(super::sha256_hex(blank.as_bytes()) == claimed)
}

/// Seconds since the epoch as `YYYY-MM-DDTHH:MM:SSZ`.
pub fn utc(secs: u64) -> String {
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    // Howard Hinnant's civil_from_days.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

/// Readable lines for the artifact table: per-file groups (`chunk/*`) that
/// all pass collapse to one row; any other group lists each file.
pub fn markdown_rows(results: &[ArtifactResult]) -> Vec<[String; 4]> {
    let mut rows = Vec::new();
    let mut i = 0;
    let count = |d: &Option<nomics_gate::ArtifactDigest>| d.as_ref().map_or_else(|| "none".to_owned(), |d| thousands(d.count));
    while i < results.len() {
        let r = &results[i];
        let group = r.name.split_once('/').map(|(g, _)| g);
        let mut j = i + 1;
        if let Some(g) = group {
            while j < results.len() && results[j].name.split_once('/').map(|(h, _)| h) == Some(g) {
                j += 1;
            }
        }
        let members = &results[i..j];
        if let Some(g) = group.filter(|_| members.len() > 1 && members.iter().all(|m| m.status == Status::Pass)) {
            let records: u64 = members.iter().filter_map(|m| m.reference.as_ref().map(|d| d.count)).sum();
            rows.push([
                format!("{g}/ ({} files)", members.len()),
                r.normalization.clone(),
                format!("{} / {} in total", thousands(records), thousands(records)),
                "PASS".to_owned(),
            ]);
        } else {
            for m in members {
                let status = if m.status == Status::Pass {
                    "PASS".to_owned()
                } else {
                    format!("**{}**", m.status.name())
                };
                rows.push([
                    m.name.clone(),
                    m.normalization.clone(),
                    format!("{} / {}", count(&m.reference), count(&m.nomics)),
                    status,
                ]);
            }
        }
        i = j;
    }
    rows
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::descriptor::Comparator;
    use crate::verify::gate::Difference;
    use nomics_gate::ArtifactDigest;

    fn result(name: &str, comparator: Comparator, r: Option<u64>, n: Option<u64>, status: Status) -> ArtifactResult {
        let d = |c: Option<u64>, h: char| c.map(|count| ArtifactDigest { sha256: h.to_string().repeat(64), count });
        ArtifactResult {
            name: name.to_owned(),
            comparator,
            normalization: "n".to_owned(),
            reference: d(r, 'a'),
            nomics: d(n, if status == Status::Fail { 'b' } else { 'a' }),
            status,
            first_difference: (status == Status::Fail).then(|| Difference { index: Some(3), detail: "x".to_owned() }),
        }
    }

    #[test]
    fn verdicts_fail_closed() {
        let pass = vec![
            result("md_bam", Comparator::Bam, Some(10), Some(10), Status::Pass),
            result("vcf", Comparator::Vcf, Some(500), Some(500), Status::Pass),
        ];
        assert_eq!(verdict(&pass, 200).0, Verdict::Pass);
        // Too few calls to mean anything.
        assert_eq!(verdict(&pass, 1000).0, Verdict::Error);
        // Empty on both sides is an ERROR, never a PASS (LEARNINGS.md).
        let empty = vec![result("vcf", Comparator::Vcf, Some(0), Some(0), Status::Empty)];
        assert_eq!(verdict(&empty, 0).0, Verdict::Error);
        let mut fail = pass.clone();
        fail[1] = result("vcf", Comparator::Vcf, Some(500), Some(500), Status::Fail);
        let (v, why) = verdict(&fail, 200);
        assert_eq!(v, Verdict::Fail);
        assert!(why.contains("vcf"));
        // Missing from nomics is a FAIL; missing from the reference is an ERROR.
        let mut gone = pass.clone();
        gone[0] = result("md_bam", Comparator::Bam, Some(10), None, Status::Missing);
        assert_eq!(verdict(&gone, 200).0, Verdict::Fail);
        gone[0] = result("md_bam", Comparator::Bam, None, Some(10), Status::Missing);
        assert_eq!(verdict(&gone, 200).0, Verdict::Error);
        assert_eq!(verdict(&[], 0).0, Verdict::Error);
    }

    #[test]
    fn seal_is_checkable_and_edits_show() {
        let report = Json::obj(vec![
            ("schema", Json::str(SCHEMA)),
            ("verdict", Json::str("PASS")),
            ("n", Json::Int(3)),
            ("t", Json::Num(1.5)),
            ("s", Json::str("a\"b\\c\nd")),
            ("report_sha256", Json::str("")),
        ]);
        let (text, hash) = seal(report.clone());
        assert_eq!(hash.len(), 64);
        assert_eq!(check_seal(&text), Some(true));
        assert_eq!(seal(report).0, text, "sealing is deterministic");
        let edited = text.replace("\"PASS\"", "\"FAIL\"");
        assert_eq!(check_seal(&edited), Some(false));
        assert!(text.contains("\"a\\\"b\\\\c\\nd\""));
        assert!(text.contains("\"t\": 1.500"));
    }

    #[test]
    fn utc_dates() {
        assert_eq!(utc(0), "1970-01-01T00:00:00Z");
        assert_eq!(utc(951_782_400), "2000-02-29T00:00:00Z");
        assert_eq!(utc(1_793_318_400 + 3_723), "2026-10-30T01:02:03Z");
    }

    #[test]
    fn groups_collapse_only_when_all_pass() {
        let rs = vec![
            result("chunk/1", Comparator::FastqChunk, Some(4), Some(4), Status::Pass),
            result("chunk/2", Comparator::FastqChunk, Some(6), Some(6), Status::Pass),
            result("mapped/1", Comparator::Bam, Some(4), Some(4), Status::Pass),
            result("mapped/2", Comparator::Bam, Some(6), Some(6), Status::Fail),
            result("vcf", Comparator::Vcf, Some(1234), Some(1234), Status::Pass),
        ];
        let rows = markdown_rows(&rs);
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0][0], "chunk/ (2 files)");
        assert_eq!(rows[2][3], "**FAIL**");
        assert_eq!(rows[3][2], "1,234 / 1,234");
    }
}
