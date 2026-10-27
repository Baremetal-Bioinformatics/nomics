//! `gatk` (GATK 4.6.2.0 wrapper + Barclay argument parser).
//!
//! The wrapper removes the first `--java-options <value>` wherever it
//! appears; any other wrapper option (`--list`, `--dry-run`, `--version`,
//! `--spark-runner`, ...) is delegated. The tool name follows, matched
//! exactly. Its arguments are parsed as Barclay does: `--long value`,
//! `--long=value`, `-short value`; booleans take an optional `true`/`false`.
//! Picard tools (MarkDuplicates) use their `UPPER_CASE` long names; the
//! legacy `KEY=value` syntax is not accepted by GATK and is delegated.
//!
//! Supported (the sarek 3.10.0 germline command lines):
//! - `MarkDuplicates`: `--INPUT/-I`+, `--OUTPUT/-O`, `--METRICS_FILE/-M`,
//!   `--TMP_DIR`, `--REFERENCE_SEQUENCE/-R` (BAM inputs only),
//!   `--REMOVE_DUPLICATES false`, `--VALIDATION_STRINGENCY LENIENT|SILENT`,
//!   `--CREATE_INDEX`.
//! - `BaseRecalibrator`: `--input/-I`, `--output/-O`, `--reference/-R`,
//!   `--intervals/-L` (one BED), `--known-sites`+, `--tmp-dir`.
//! - `ApplyBQSR`: `--input/-I`, `--output/-O`, `--reference/-R`,
//!   `--bqsr-recal-file`, `--intervals/-L` (one BED), `--tmp-dir`,
//!   `--create-output-bam-index`.
//! - `HaplotypeCaller`: `--input/-I`, `--output/-O` (`.vcf.gz`),
//!   `--reference/-R`, `--dbsnp/-D`, `--intervals/-L` (one BED),
//!   `--native-pair-hmm-threads`, `--tmp-dir`, `--pcr-indel-model
//!   CONSERVATIVE`, `--emit-ref-confidence/-ERC NONE|GVCF`,
//!   `--create-output-variant-index`.
//!
//! Indexes (`.bai` named like Picard/GATK: `x.bam` -> `x.bai`; `.vcf.gz.tbi`)
//! are written with htsjdk's layout (`nomics_io::htsjdk_index`).

use std::path::{Path, PathBuf};

use super::args::{Kind, Spec, Unsupported, parse, spec};
use super::{Decision, threads};

const MARKDUPLICATES: &[Spec] = &[
    spec(&["--INPUT", "-I"], Kind::Multi, "INPUT"),
    spec(&["--OUTPUT", "-O"], Kind::Value, "OUTPUT"),
    spec(&["--METRICS_FILE", "-M"], Kind::Value, "METRICS_FILE"),
    spec(&["--TMP_DIR"], Kind::Multi, "TMP_DIR"),
    spec(&["--REFERENCE_SEQUENCE", "-R"], Kind::Value, "REFERENCE_SEQUENCE"),
    spec(&["--REMOVE_DUPLICATES"], Kind::Bool, "REMOVE_DUPLICATES"),
    spec(&["--VALIDATION_STRINGENCY"], Kind::Value, "VALIDATION_STRINGENCY"),
    spec(&["--CREATE_INDEX"], Kind::Bool, "CREATE_INDEX"),
];

const BASERECALIBRATOR: &[Spec] = &[
    spec(&["--input", "-I"], Kind::Multi, "input"),
    spec(&["--output", "-O"], Kind::Value, "output"),
    spec(&["--reference", "-R"], Kind::Value, "reference"),
    spec(&["--intervals", "-L"], Kind::Multi, "intervals"),
    spec(&["--known-sites"], Kind::Multi, "known-sites"),
    spec(&["--tmp-dir"], Kind::Value, "tmp-dir"),
];

const APPLYBQSR: &[Spec] = &[
    spec(&["--input", "-I"], Kind::Multi, "input"),
    spec(&["--output", "-O"], Kind::Value, "output"),
    spec(&["--reference", "-R"], Kind::Value, "reference"),
    spec(&["--bqsr-recal-file", "-bqsr"], Kind::Value, "bqsr-recal-file"),
    spec(&["--intervals", "-L"], Kind::Multi, "intervals"),
    spec(&["--tmp-dir"], Kind::Value, "tmp-dir"),
    spec(
        &["--create-output-bam-index", "-OBI"],
        Kind::Bool,
        "create-output-bam-index",
    ),
];

const HAPLOTYPECALLER: &[Spec] = &[
    spec(&["--input", "-I"], Kind::Multi, "input"),
    spec(&["--output", "-O"], Kind::Value, "output"),
    spec(&["--reference", "-R"], Kind::Value, "reference"),
    spec(&["--dbsnp", "-D"], Kind::Value, "dbsnp"),
    spec(&["--intervals", "-L"], Kind::Multi, "intervals"),
    spec(
        &["--native-pair-hmm-threads"],
        Kind::Value,
        "native-pair-hmm-threads",
    ),
    spec(&["--tmp-dir"], Kind::Value, "tmp-dir"),
    spec(&["--pcr-indel-model"], Kind::Value, "pcr-indel-model"),
    spec(
        &["--emit-ref-confidence", "-ERC"],
        Kind::Value,
        "emit-ref-confidence",
    ),
    spec(
        &["--create-output-variant-index", "-OVI"],
        Kind::Bool,
        "create-output-variant-index",
    ),
];

/// Split the wrapper's own arguments from the tool's: drop the first
/// `--java-options <value>`, then the tool name must come first.
fn split_wrapper(args: &[String]) -> Result<(String, Vec<String>), Unsupported> {
    let mut rest = args.to_vec();
    if let Some(i) = rest.iter().position(|a| a == "--java-options") {
        if i + 1 >= rest.len() {
            return Err(Unsupported::new("--java-options without a value"));
        }
        rest.drain(i..i + 2);
    }
    match rest.first() {
        None => Err(Unsupported::new("no tool (usage)")),
        Some(t) if t.starts_with('-') => Err(Unsupported::new(format!("wrapper option {t}"))),
        Some(t) => {
            let tool = t.clone();
            rest.remove(0);
            Ok((tool, rest))
        }
    }
}

fn label(tool: &str) -> String {
    format!("gatk {tool}")
}

pub fn decide(args: &[String]) -> Decision {
    let (tool, rest) = match split_wrapper(args) {
        Ok(x) => x,
        Err(why) => {
            return Decision::Delegate {
                label: "gatk".into(),
                why,
            };
        }
    };
    let planned = match tool.as_str() {
        "MarkDuplicates" => mark_duplicates(&rest),
        "BaseRecalibrator" => base_recalibrator(&rest),
        "ApplyBQSR" => apply_bqsr(&rest),
        "HaplotypeCaller" => haplotype_caller(&rest),
        _ => Err(Unsupported::new(format!("tool {tool}"))),
    };
    match planned {
        Ok(run) => Decision::Handle {
            label: label(&tool),
            run,
        },
        Err(why) => Decision::Delegate {
            label: label(&tool),
            why,
        },
    }
}

type Run = Box<dyn FnOnce() -> Result<(), String>>;

fn ends_with_ci(p: &str, ext: &str) -> bool {
    p.len() >= ext.len() && p[p.len() - ext.len()..].eq_ignore_ascii_case(ext)
}

/// One `--intervals` value naming an existing BED file.
fn one_bed(values: &[String], what: &str) -> Result<PathBuf, Unsupported> {
    match values {
        [] => Err(Unsupported::new(format!("no --{what} (whole-genome traversal)"))),
        [one] if ends_with_ci(one, ".bed") && Path::new(one).is_file() => Ok(PathBuf::from(one)),
        [one] => Err(Unsupported::new(format!("--{what} {one} (only a BED file)"))),
        _ => Err(Unsupported::new(format!("more than one --{what}"))),
    }
}

fn one_input(values: &[String], what: &str) -> Result<PathBuf, Unsupported> {
    match values {
        [one] if ends_with_ci(one, ".bam") => Ok(PathBuf::from(one)),
        [one] => Err(Unsupported::new(format!("--{what} {one} (only BAM)"))),
        [] => Err(Unsupported::new(format!("missing --{what}"))),
        _ => Err(Unsupported::new(format!("more than one --{what}"))),
    }
}

/// Picard/GATK index name for a BAM: `x.bam` -> `x.bai`.
fn bai_path(bam: &Path) -> PathBuf {
    bam.with_extension("bai")
}

fn io_err(what: &str) -> impl Fn(std::io::Error) -> String + '_ {
    move |e| format!("{what}: {e}")
}

fn mark_duplicates(args: &[String]) -> Result<Run, Unsupported> {
    let p = parse(args, MARKDUPLICATES, false)?;
    let inputs: Vec<PathBuf> = p.all("INPUT").iter().map(PathBuf::from).collect();
    if inputs.is_empty() {
        return Err(Unsupported::new("missing --INPUT"));
    }
    if let Some(i) = p.all("INPUT").iter().find(|i| !ends_with_ci(i, ".bam")) {
        return Err(Unsupported::new(format!("--INPUT {i} (only BAM)")));
    }
    let output = PathBuf::from(p.req("OUTPUT")?);
    if !ends_with_ci(&output.to_string_lossy(), ".bam") {
        return Err(Unsupported::new(format!(
            "--OUTPUT {} (only BAM)",
            output.display()
        )));
    }
    let metrics = PathBuf::from(p.req("METRICS_FILE")?);
    if p.flag("REMOVE_DUPLICATES", false) {
        return Err(Unsupported::new("--REMOVE_DUPLICATES true"));
    }
    match p.one("VALIDATION_STRINGENCY") {
        Some("LENIENT" | "SILENT") => {}
        Some(v) => return Err(Unsupported::new(format!("--VALIDATION_STRINGENCY {v}"))),
        None => return Err(Unsupported::new("default VALIDATION_STRINGENCY STRICT")),
    }
    let create_index = p.flag("CREATE_INDEX", false);
    let cl = format!("MarkDuplicates {}", args.join(" "));
    Ok(Box::new(move || {
        let opts = nomics_markdup::Options {
            inputs,
            output: output.clone(),
            metrics,
            pg_version: format!("Version:nomics-{}", env!("CARGO_PKG_VERSION")),
            pg_command_line: cl,
            threads: threads(None),
        };
        nomics_markdup::run(&opts).map_err(io_err("MarkDuplicates"))?;
        if create_index {
            nomics_io::htsjdk_index::write_bai_htsjdk(&output, &bai_path(&output))
                .map_err(io_err("index"))?;
        }
        Ok(())
    }))
}

fn base_recalibrator(args: &[String]) -> Result<Run, Unsupported> {
    let p = parse(args, BASERECALIBRATOR, false)?;
    let input = one_input(p.all("input"), "input")?;
    let output = PathBuf::from(p.req("output")?);
    let reference = PathBuf::from(p.req("reference")?);
    let intervals = one_bed(p.all("intervals"), "intervals")?;
    let known: Vec<PathBuf> = p.all("known-sites").iter().map(PathBuf::from).collect();
    if known.is_empty() {
        return Err(Unsupported::new("no --known-sites"));
    }
    if let Some(k) = p.all("known-sites").iter().find(|k| !ends_with_ci(k, ".vcf.gz")) {
        return Err(Unsupported::new(format!("--known-sites {k} (only .vcf.gz)")));
    }
    let cl = format!("BaseRecalibrator {}", args.join(" "));
    Ok(Box::new(move || {
        let opts = nomics_bqsr::Options {
            input,
            reference,
            known_sites: known,
            intervals,
            scatter: None,
            known_sites_cache: None,
            output_table: output,
            output_bam: None,
            pg_command_line: cl,
            threads: threads(None),
        };
        nomics_bqsr::run_streaming(&opts).map_err(io_err("BaseRecalibrator"))
    }))
}

fn apply_bqsr(args: &[String]) -> Result<Run, Unsupported> {
    let p = parse(args, APPLYBQSR, false)?;
    let input = one_input(p.all("input"), "input")?;
    let output = PathBuf::from(p.req("output")?);
    if !ends_with_ci(&output.to_string_lossy(), ".bam") {
        return Err(Unsupported::new(format!(
            "--output {} (only BAM)",
            output.display()
        )));
    }
    let reference = PathBuf::from(p.req("reference")?);
    let table = PathBuf::from(p.req("bqsr-recal-file")?);
    let intervals = one_bed(p.all("intervals"), "intervals")?;
    // ApplyBQSR re-reads its arguments from the table: only ours are implemented.
    let text = std::fs::read_to_string(&table)
        .map_err(|e| Unsupported::new(format!("--bqsr-recal-file unreadable: {e}")))?;
    nomics_bqsr::check_report(&text).map_err(Unsupported::new)?;
    let create_index = p.flag("create-output-bam-index", true);
    let cl = format!("ApplyBQSR {}", args.join(" "));
    Ok(Box::new(move || {
        let opts = nomics_bqsr::Options {
            input,
            reference,
            known_sites: Vec::new(),
            intervals,
            scatter: None,
            known_sites_cache: None,
            output_table: PathBuf::new(),
            output_bam: Some(output.clone()),
            pg_command_line: cl,
            threads: threads(None),
        };
        nomics_bqsr::apply_from_report(&opts, &table).map_err(io_err("ApplyBQSR"))?;
        if create_index {
            nomics_io::htsjdk_index::write_bai_htsjdk(&output, &bai_path(&output))
                .map_err(io_err("index"))?;
        }
        Ok(())
    }))
}

fn haplotype_caller(args: &[String]) -> Result<Run, Unsupported> {
    let p = parse(args, HAPLOTYPECALLER, false)?;
    let input = one_input(p.all("input"), "input")?;
    let output = PathBuf::from(p.req("output")?);
    if !ends_with_ci(&output.to_string_lossy(), ".vcf.gz") {
        return Err(Unsupported::new(format!(
            "--output {} (only .vcf.gz)",
            output.display()
        )));
    }
    let reference = PathBuf::from(p.req("reference")?);
    let dbsnp = PathBuf::from(
        p.one("dbsnp")
            .ok_or_else(|| Unsupported::new("no --dbsnp"))?,
    );
    let intervals = one_bed(p.all("intervals"), "intervals")?;
    if let Some(n) = p.one("native-pair-hmm-threads")
        && n.parse::<u32>().is_err()
    {
        return Err(Unsupported::new(format!("--native-pair-hmm-threads {n}")));
    }
    // GATK's default is CONSERVATIVE.
    p.only("pcr-indel-model", &["CONSERVATIVE"])?;
    let gvcf = match p.one("emit-ref-confidence") {
        None | Some("NONE") => false,
        Some("GVCF") => true,
        Some(v) => return Err(Unsupported::new(format!("--emit-ref-confidence {v}"))),
    };
    if gvcf && !ends_with_ci(&output.to_string_lossy(), ".g.vcf.gz") {
        return Err(Unsupported::new(format!(
            "-ERC GVCF with --output {} (not .g.vcf.gz)",
            output.display()
        )));
    }
    let create_index = p.flag("create-output-variant-index", true);
    Ok(Box::new(move || {
        crate::pipeline::haplotypecaller(
            &input,
            &[intervals],
            &reference,
            &dbsnp,
            &output,
            threads(None),
            gvcf,
        )
        .map_err(io_err("HaplotypeCaller"))?;
        if create_index {
            nomics_io::htsjdk_index::write_tbi_htsjdk_vcf(
                &output,
                &super::with_suffix(&output, ".tbi"),
            )
            .map_err(io_err("index"))?;
        }
        Ok(())
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_owned).collect()
    }

    #[test]
    fn wrapper_java_options_anywhere() {
        let (t, rest) =
            split_wrapper(&v("--java-options -Xmx3g ApplyBQSR --input a.bam")).unwrap();
        assert_eq!(t, "ApplyBQSR");
        assert_eq!(rest, v("--input a.bam"));
        let (t, _) = split_wrapper(&v("HaplotypeCaller --java-options -Xmx1g -I a.bam")).unwrap();
        assert_eq!(t, "HaplotypeCaller");
        assert!(split_wrapper(&v("--list")).is_err());
        assert!(split_wrapper(&v("--version")).is_err());
    }

    #[test]
    fn unsupported_tools_and_values_delegate() {
        assert!(matches!(
            decide(&v("VariantFiltration -V x.vcf.gz")),
            Decision::Delegate { .. }
        ));
        assert!(
            haplotype_caller(&v("--input a.bam --output o.vcf.gz --pcr-indel-model NONE"))
                .is_err()
        );
        let e = mark_duplicates(&v("--INPUT a.bam --OUTPUT o.bam --METRICS_FILE m"))
            .err()
            .unwrap();
        assert_eq!(e.0, "default VALIDATION_STRINGENCY STRICT");
        let e = mark_duplicates(&v(
            "--INPUT a.bam --OUTPUT o.bam --METRICS_FILE m --VALIDATION_STRINGENCY LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500",
        ))
        .err()
        .unwrap();
        assert_eq!(e.0, "--OPTICAL_DUPLICATE_PIXEL_DISTANCE");
        // Legacy Picard syntax is not GATK syntax.
        assert!(mark_duplicates(&v("I=a.bam O=o.bam M=m")).is_err());
    }

    #[test]
    fn index_names() {
        assert_eq!(bai_path(Path::new("x/HG002.md.bam")), PathBuf::from("x/HG002.md.bai"));
    }
}
