//! FASTQ pair orchestration. No caller algorithms live here.
use crate::{pipeline, scatter, somatic_chain::Chain};
use std::{
    collections::BTreeSet,
    fs, io,
    path::{Path, PathBuf},
    time::Instant,
};

fn bad(row: usize, message: impl std::fmt::Display) -> io::Error {
    io::Error::other(format!("samplesheet row {row}: {message}"))
}
#[derive(Debug)]
struct Sample {
    patient: String,
    name: String,
    status: String,
    sex: String,
    row: usize,
    lanes: Vec<String>,
    fq1: Vec<PathBuf>,
    fq2: Vec<PathBuf>,
}
// Single-line CSV records, including quoted commas and escaped quotes. Multiline
// fields cannot be identifiers/FASTQ paths here and are explicitly refused.
fn csv(line: &str, row: usize) -> io::Result<Vec<String>> {
    let mut fields = Vec::new();
    let mut chars = line.trim_end_matches('\r').chars().peekable();
    loop {
        let mut field = String::new();
        if chars.peek() == Some(&'"') {
            chars.next();
            loop {
                match chars.next() {
                    Some('"') if chars.peek() == Some(&'"') => {
                        chars.next();
                        field.push('"');
                    }
                    Some('"') => break,
                    Some(c) => field.push(c),
                    None => {
                        return Err(bad(
                            row,
                            "unterminated quote / multiline CSV is unsupported",
                        ));
                    }
                }
            }
            if chars.peek().is_some_and(|c| *c != ',') {
                return Err(bad(row, "characters after closing quote"));
            }
        } else {
            while chars.peek().is_some_and(|c| *c != ',') {
                let c = chars.next().unwrap();
                if c == '"' {
                    return Err(bad(row, "quote in unquoted field"));
                }
                field.push(c);
            }
        }
        fields.push(field);
        if chars.next().is_none() {
            break;
        }
    }
    Ok(fields)
}
fn parse(text: &str) -> io::Result<[Sample; 2]> {
    let mut lines = text.lines().enumerate();
    let (_, header) = lines.next().ok_or_else(|| bad(1, "empty samplesheet"))?;
    let header = csv(header, 1)?;
    let columns = [
        "patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2",
    ];
    if header.len() != columns.len()
        || header.iter().collect::<BTreeSet<_>>().len() != columns.len()
    {
        return Err(bad(
            1,
            "expected exactly patient,sex,status,sample,lane,fastq_1,fastq_2",
        ));
    }
    let positions = columns
        .iter()
        .map(|c| {
            header
                .iter()
                .position(|h| h == c)
                .ok_or_else(|| bad(1, format!("missing column {c}")))
        })
        .collect::<io::Result<Vec<_>>>()?;
    let mut samples: Vec<Sample> = Vec::new();
    for (index, line) in lines {
        let row = index + 1;
        let fields = csv(line, row)?;
        if fields.len() != header.len() {
            return Err(bad(row, "wrong number of columns"));
        }
        let f: Vec<&str> = positions.iter().map(|&i| fields[i].as_str()).collect();
        for i in [0, 3, 4] {
            pipeline::check_identifier(columns[i], f[i]).map_err(|e| bad(row, e))?;
            // Chain publication labels currently permit this smaller alphabet.
            if f[i].contains('+') {
                return Err(bad(row, "'+' is unsupported in pair identifiers"));
            }
        }
        if !["0", "1"].contains(&f[2]) {
            return Err(bad(row, "status must be 0 (normal) or 1 (tumour)"));
        }
        if !["XX", "XY", "NA"].contains(&f[1]) {
            return Err(bad(row, "sex must be XX, XY or NA"));
        }
        for path in &f[5..] {
            if path.is_empty() || path.chars().any(char::is_control) {
                return Err(bad(row, "empty or control-containing FASTQ path"));
            }
        }
        if f[5] == f[6] {
            return Err(bad(row, "mates must name different FASTQs"));
        }
        if let Some(first) = samples.first() {
            if first.patient != f[0] {
                return Err(bad(row, "exactly one patient is supported"));
            }
            if first.sex != f[1] {
                return Err(bad(row, "patient sex differs between rows"));
            }
        }
        let at = if let Some(i) = samples.iter().position(|s| s.name == f[3]) {
            if samples[i].status != f[2] {
                return Err(bad(row, "sample appears with two statuses"));
            }
            i
        } else {
            if samples.iter().any(|s| s.status == f[2]) {
                return Err(bad(
                    row,
                    "exactly one normal and one tumour sample are supported",
                ));
            }
            samples.push(Sample {
                patient: f[0].into(),
                name: f[3].into(),
                sex: f[1].into(),
                status: f[2].into(),
                row,
                lanes: Vec::new(),
                fq1: Vec::new(),
                fq2: Vec::new(),
            });
            samples.len() - 1
        };
        let s = &mut samples[at];
        if s.lanes.iter().any(|l| l == f[4]) {
            return Err(bad(row, "duplicate lane for sample"));
        }
        s.lanes.push(f[4].into());
        s.fq1.push(f[5].into());
        s.fq2.push(f[6].into());
    }
    if samples.len() != 2 {
        return Err(bad(
            samples.first().map_or(1, |s| s.row),
            "exactly one normal and one tumour are required",
        ));
    }
    samples.sort_by(|a, b| a.status.cmp(&b.status));
    Ok(samples.try_into().unwrap())
}
fn sample_args(base: &pipeline::Args, s: &Sample, outdir: PathBuf) -> pipeline::Args {
    let mut args = base.clone();
    args.bwa_mismatch_penalty =
        nomics_align::options::Options::new(if s.status == "1" { 3 } else { 4 }).unwrap();
    eprintln!(
        "sample {} status {}: BWA mismatch penalty {}",
        s.name,
        s.status,
        args.bwa_mismatch_penalty.mismatch_penalty()
    );
    args.sample = format!("{}_{}", s.patient, s.name);
    args.prefix = s.name.clone();
    args.fastq1 = s.fq1.clone();
    args.fastq2 = s.fq2.clone();
    args.lanes = s.lanes.clone();
    args.outdir = outdir;
    args.stop_after_bqsr = true;
    args
}
fn publish(root: &Path, work: &Path, sample: &str, mapped: bool) -> io::Result<PathBuf> {
    let pre = root.join("results/preprocessing");
    for (stage, suffix) in [
        ("markduplicates", "md.bam"),
        ("recalibrated", "recal.bam"),
        ("recal_table", "recal.table"),
    ] {
        let dir = pre.join(stage).join(sample);
        fs::create_dir_all(&dir)?;
        let name = format!("{sample}.{suffix}");
        fs::rename(work.join(&name), dir.join(&name))?;
    }
    if mapped {
        fs::create_dir_all(pre.join("mapped"))?;
        fs::rename(work.join("mapped"), pre.join("mapped").join(sample))?;
    }
    let reports = root.join("results/reports/markduplicates").join(sample);
    fs::create_dir_all(&reports)?;
    fs::rename(
        work.join(format!("{sample}.md.metrics")),
        reports.join(format!("{sample}.md.bam.metrics")),
    )?;
    let bam = pre
        .join("recalibrated")
        .join(sample)
        .join(format!("{sample}.recal.bam"));
    nomics_io::htsjdk_index::write_bai_htsjdk(&bam, &bam.with_extension("bai"))?;
    Ok(bam)
}
#[allow(clippy::too_many_arguments)]
pub fn run(
    base: pipeline::Args,
    sheet: &Path,
    germline_resource: PathBuf,
    panel_of_normals: PathBuf,
    pileup_resource: PathBuf,
    mut orientation_order: Vec<String>,
    scatter_orders: [Option<PathBuf>; 2],
) -> io::Result<()> {
    let started = Instant::now();
    let samples = parse(&fs::read_to_string(sheet)?)?;
    if base.resume_from.is_some()
        || base.emit_gvcf
        || base.reports
        || base.keep_bams
        || !base.save_output_as_bam
    {
        return Err(io::Error::other(
            "pair mode requires BAM output and does not support --resume-from, --hc-only, --emit-gvcf, --reports or --keep-bams",
        ));
    }
    if std::env::var_os("NOMICS_STOP_AFTER").is_some() {
        return Err(io::Error::other(
            "NOMICS_STOP_AFTER is unsupported in pair mode",
        ));
    }
    if base.outdir.symlink_metadata().is_ok() {
        return Err(io::Error::other(
            "pair outdir must not exist (no resume/overwrite)",
        ));
    }
    let groups = scatter::split_bed(&fs::read_to_string(&base.targets)?, 200000.0)?;
    let names: Vec<String> = groups
        .iter()
        .map(|g| g.name.trim_end_matches(".bed").into())
        .collect();
    if names.is_empty() {
        return Err(io::Error::other("empty target BED"));
    }
    scatter::orders(&groups, base.scatter_order.as_deref())?;
    if orientation_order.is_empty() {
        orientation_order = names.clone();
        eprintln!(
            "pair: --orientation-order absent; using scatter-group order: {} (sarek completion order may differ)",
            names.join(",")
        );
    }
    if orientation_order.len() != names.len()
        || orientation_order.iter().collect::<BTreeSet<_>>()
            != names.iter().collect::<BTreeSet<_>>()
    {
        return Err(io::Error::other(
            "orientation-order must contain every BED stem exactly once",
        ));
    }
    for order in scatter_orders.iter().flatten() {
        scatter::orders(&groups, Some(order))?;
    }
    // Validate both samples before creating staging or starting the normal.
    for sample in &samples {
        let args = sample_args(&base, sample, base.outdir.clone());
        pipeline::validate_preprocessing(&args).map_err(|e| bad(sample.row, e))?;
    }
    for resource in [&germline_resource, &panel_of_normals, &pileup_resource] {
        fs::File::open(resource)?;
    }
    let parent = base
        .outdir
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let staging = parent.join(format!(".nomics-pair-{}.incomplete", std::process::id()));
    fs::create_dir(&staging)?;
    let result = (|| {
        let mut bams = Vec::new();
        let mut timing = String::from("stage\twall_s\tpeak_rss_bytes\n");
        for (i, sample) in samples.iter().enumerate() {
            let t = Instant::now();
            let work = staging.join("work/preprocessing").join(&sample.name);
            eprintln!(
                "pair: preprocessing {} (SM={}_{})",
                sample.name, sample.patient, sample.name
            );
            let mut args = sample_args(&base, sample, work.clone());
            if scatter_orders[i].is_some() {
                args.scatter_order = scatter_orders[i].clone();
            }
            pipeline::run(args)?;
            // run returned: all per-sample Rust buffers and worker threads are gone.
            bams.push(publish(&staging, &work, &sample.name, base.save_mapped)?);
            let stats = staging.join("timings");
            fs::create_dir_all(&stats)?;
            fs::copy(
                work.join("stages.tsv"),
                stats.join(format!("{}.tsv", sample.name)),
            )?;
            timing.push_str(&format!(
                "{}\t{:.3}\t{}\n",
                sample.name,
                t.elapsed().as_secs_f64(),
                nomics_mem::peak_rss_bytes().unwrap_or(0)
            ));
            eprintln!(
                "pair: {} buffers released; process peak RSS {} bytes",
                sample.name,
                nomics_mem::peak_rss_bytes().unwrap_or(0)
            );
        }
        let beds = scatter::write_groups(&groups, &staging.join("intervals"))?;
        let t = Instant::now();
        Chain {
            normal: bams[0].clone(),
            tumour: bams[1].clone(),
            reference: pipeline::reference_paths(&base).0,
            intervals: beds,
            intervals_dir: None,
            germline_resource,
            panel_of_normals,
            pileup_resource,
            outdir: staging.join("chain"),
            normal_label: Some(samples[0].name.clone()),
            tumour_label: Some(samples[1].name.clone()),
            orientation_order,
            keep_intermediates: base.write_intermediates,
            threads: base.threads,
        }
        .run()?;
        fs::rename(
            staging.join("chain/variant_calling"),
            staging.join("results/variant_calling"),
        )?;
        if base.write_intermediates {
            fs::rename(staging.join("chain/work"), staging.join("work/chain"))?;
        }
        fs::rename(
            staging.join("chain/chain-timings.tsv"),
            staging.join("chain-timings.tsv"),
        )?;
        fs::remove_dir(staging.join("chain"))?;
        if !base.write_intermediates {
            fs::remove_dir_all(staging.join("work"))?;
        }
        timing.push_str(&format!(
            "somatic_chain\t{:.3}\t{}\nTOTAL\t{:.3}\t{}\n",
            t.elapsed().as_secs_f64(),
            nomics_mem::peak_rss_bytes().unwrap_or(0),
            started.elapsed().as_secs_f64(),
            nomics_mem::peak_rss_bytes().unwrap_or(0)
        ));
        eprintln!("{timing}");
        fs::write(staging.join("stages.tsv"), timing)?;
        fs::rename(&staging, &base.outdir)
    })();
    if result.is_err() {
        eprintln!(
            "pair: incomplete staging retained at {} (not resumable)",
            staging.display()
        );
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    const HEADER: &str = "patient,sex,status,sample,lane,fastq_1,fastq_2\n";
    const NORMAL: &str = "P,XX,0,N,L1,n1,n2\n";
    const TUMOUR: &str = "P,XX,1,T,L1,t1,t2\n";
    #[test]
    fn sample_scoring_follows_status() {
        let dir = std::env::temp_dir().join(format!("nomics-pair-scoring-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let base = crate::pipeline::identity_tests::args(&dir, &["L1"], 1);
        let samples = parse(&format!("{HEADER}{TUMOUR}{NORMAL}")).unwrap();
        for (sample, expected) in samples.iter().zip([4, 3]) {
            let args = sample_args(&base, sample, dir.clone());
            assert_eq!(args.bwa_mismatch_penalty.mismatch_penalty(), expected);
            assert!(args.stop_after_bqsr);
        }
        assert_eq!(base.bwa_mismatch_penalty.mismatch_penalty(), 4);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn lanes_and_csv() {
        let p = parse(&format!(
            "{HEADER}{TUMOUR}{NORMAL}P,XX,0,N,L2,\"n,3\",\"n\"\"4\"\n"
        ))
        .unwrap();
        assert_eq!(p[0].name, "N");
        assert_eq!(p[0].lanes, ["L1", "L2"]);
        assert_eq!(p[0].fq1[1], Path::new("n,3"));
        assert_eq!(p[0].fq2[1], Path::new("n\"4"));
    }
    #[test]
    fn refusals_name_the_row() {
        for (body, row, message) in [
            (TUMOUR.to_owned(), 2, "exactly one"),
            (NORMAL.to_owned(), 2, "exactly one"),
            (
                format!("{NORMAL}{TUMOUR}P,XX,1,T2,L1,a,b\n"),
                4,
                "exactly one",
            ),
            (format!("{NORMAL}Q,XX,1,T,L1,a,b\n"), 3, "one patient"),
            (format!("{NORMAL}P,XX,2,T,L1,a,b\n"), 3, "status must"),
            (format!("{NORMAL}P,XX,1,N,L2,a,b\n"), 3, "two statuses"),
            (format!("{NORMAL}{NORMAL}{TUMOUR}"), 3, "duplicate lane"),
            (format!("{NORMAL}P,XX,1,../T,L1,a,b\n"), 3, "sample"),
            (format!("{NORMAL}P,XY,1,T,L1,a,b\n"), 3, "sex differs"),
            (format!("{NORMAL}P,XX,1,T,L1,a,a\n"), 3, "different FASTQs"),
            (format!("{NORMAL}P,XX,1,T,L1,,a\n"), 3, "empty"),
        ] {
            let e = parse(&format!("{HEADER}{body}")).unwrap_err().to_string();
            assert!(
                e.contains(&format!("row {row}:")) && e.contains(message),
                "{e}"
            );
        }
        for line in [
            "P,XX,1,T,L1,a",
            "P,XX,1,T,L1,\"a,b",
            "P,XX,1,T,L1,\"a\"x,b",
            "P,XX,1,T,L1,a\"x,b",
            "P,bad,1,T,L1,a,b",
            "../P,XX,1,T,L1,a,b",
            "P,XX,1,T,../L,a,b",
            "P,XX,1,T+,L,a,b",
        ] {
            assert!(
                parse(&format!("{HEADER}{NORMAL}{line}\n"))
                    .unwrap_err()
                    .to_string()
                    .contains("row 3:")
            );
        }
        for header in [
            "",
            "patient,patient,status,sample,lane,fastq_1,fastq_2\n",
            "patient,bad,status,sample,lane,fastq_1,fastq_2\n",
        ] {
            assert!(parse(header).is_err());
        }
    }
}
