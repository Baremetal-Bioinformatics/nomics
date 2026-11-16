//! One shared queue across groups and contigs; ordered, single-writer reduction.
use super::*;
use s::mutect2_regions::{Inputs, plan};

struct Part {
    group: usize,
    contig: usize,
    prefix: PathBuf,
}
fn combine(paths: &[PathBuf], output: &Path, vcf: bool) -> io::Result<()> {
    if vcf {
        return combine_vcf(paths, output);
    }
    let mut writer = io::BufWriter::new(fs::File::create(output)?);
    for (i, path) in paths.iter().enumerate() {
        for (line_no, line) in BufReader::new(fs::File::open(path)?).lines().enumerate() {
            let line = line?;
            if i == 0 || line_no >= 2 {
                writeln!(writer, "{line}")?;
            }
        }
    }
    writer.flush()
}
fn combine_vcf(paths: &[PathBuf], output: &Path) -> io::Result<()> {
    let mut writer = nomics_io::bgzf::ParallelBgzfWriter::new(fs::File::create(output)?, 1);
    for (i, path) in paths.iter().enumerate() {
        let mut reader = BufReader::new(nomics_io::gzip::GzReader::open(path)?);
        let mut line = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if i == 0 || !line.starts_with(b"#") {
                writer.write_all(&line)?;
            }
        }
    }
    writer.finish()?;
    Ok(())
}
fn identical(a: &Path, b: &Path) -> io::Result<()> {
    let mut a = BufReader::new(fs::File::open(a)?);
    let mut b = BufReader::new(fs::File::open(b)?);
    loop {
        let x = a.fill_buf()?;
        let y = b.fill_buf()?;
        let n = x.len().min(y.len());
        if x[..n] != y[..n] || (n == 0 && x.len() != y.len()) {
            return Err(bad("contig self-check: serial output differs"));
        }
        if n == 0 {
            return Ok(());
        }
        a.consume(n);
        b.consume(n);
    }
}
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    chain: &Chain,
    work: &Path,
    beds: &[PathBuf],
    names: &[String],
    normal: &str,
    inputs: &[PathBuf],
    vcfs: &[PathBuf],
    stats: &[PathBuf],
    archives: &[PathBuf],
    tp: &[PathBuf],
    np: &[PathBuf],
    labels: [&str; 2],
    timing: &Timings,
) -> io::Result<()> {
    let self_check = std::env::var_os("NOMICS_CONTIG_SELF_CHECK").is_some();
    // These HC helpers own process-wide diagnostic sinks (including a shared
    // SW record limit). Preserve their group/contig order as well as outputs.
    let ordered_debug = [
        "NOMICS_DUMP_SW",
        "NOMICS_DEBUG_ASSEMBLY",
        "NOMICS_DEBUG_GRAPH_DUMP",
    ]
    .iter()
    .any(|key| std::env::var_os(key).is_some());
    if ordered_debug && self_check {
        return Err(bad(
            "NOMICS_CONTIG_SELF_CHECK cannot be combined with SW/assembly/graph dumps: the verification traversal would duplicate global diagnostics",
        ));
    }
    let workers = if ordered_debug {
        eprintln!("contig diagnostics: using one call worker to preserve global dump order");
        1
    } else {
        chain.threads
    };
    let scratch = work.join(".contigs");
    fs::create_dir(&scratch)?;
    let header = bam(&chain.tumour)?;
    let mut parts = Vec::new();
    let mut groups = Vec::new();
    for (group, bed_path) in beds.iter().enumerate() {
        let contigs: BTreeSet<_> = bed(bed_path, &header)?.into_iter().map(|r| r.0).collect();
        let start = parts.len();
        for contig in contigs {
            parts.push(Part {
                group,
                contig,
                prefix: scratch.join(format!("{group}-{contig}")),
            });
        }
        groups.push(start..parts.len());
    }
    let input = |g: usize| Inputs {
        input: inputs,
        normal_sample: normal,
        reference: &chain.reference,
        intervals: beds[g].as_path(),
        panel_of_normals: &chain.panel_of_normals,
        germline_resource: &chain.germline_resource,
        regions: Path::new("/dev/null"),
        stats: Path::new("/dev/null"),
        dump_prefix: None,
    };
    let tapes: Vec<_> = parts.iter().map(|_| Mutex::new(None)).collect();
    timing.measure("all", "contig_plan", || {
        tasks(parts.len(), chain.threads, |i| {
            let p = &parts[i];
            let tape = timing.measure(&names[p.group], &format!("plan.{}", p.contig), || {
                plan::record(&input(p.group), p.contig, &suffix(&p.prefix, ".rng"))
            })?;
            eprintln!(
                "contig_rng\t{}\t{}\t{}",
                names[p.group],
                String::from_utf8_lossy(&header.refs[p.contig].name),
                tape.calls
            );
            *tapes[i].lock().unwrap() = Some(tape);
            Ok(())
        })
    })?;
    let mut starts = Vec::new();
    timing.measure("all", "contig_replay", || {
        for group in &groups {
            let tapes: Vec<_> = group
                .clone()
                .map(|i| tapes[i].lock().unwrap().take().unwrap())
                .collect();
            starts.extend(plan::replay(&tapes)?);
        }
        Ok(())
    })?;
    let collectors: Vec<_> = parts.iter().map(|_| Mutex::new(None)).collect();
    let af: Vec<_> = (0..parts.len() * 2).map(|_| Mutex::new(None)).collect();
    // Per-contig traversal and pileup counters, summed per group after the pool has
    // drained and printed there under the labels the gates read.
    let counts: Vec<_> = parts
        .iter()
        .map(|_| Mutex::new(s::mutect2_regions::Stats::default()))
        .collect();
    let pileup_counts: Vec<_> = (0..parts.len() * 2)
        .map(|_| Mutex::new(s::get_pileup_summaries::Stats::default()))
        .collect();
    let phase_totals = Mutex::new(vec![[std::time::Duration::ZERO; 8]; groups.len()]);
    let spans = Mutex::new(vec![[None::<(Instant, Instant)>; 3]; groups.len()]);
    tasks(parts.len() * 3, workers, |task| {
        // Start the long caller tasks first. Pileups fill slots as callers
        // finish; at one thread each group's step stays contiguous.
        let (i, kind) = (task % parts.len(), task / parts.len());
        let p = &parts[i];
        let start = Instant::now();
        let step = ["caller", "pileup_tumour", "pileup_normal"][kind];
        timing.measure(&names[p.group], &format!("{step}.{}", p.contig), || {
            if kind == 0 {
                let (phases, collector, traversal) = s::mutect2_calls::vcf::run_part(
                    &input(p.group),
                    &suffix(&p.prefix, ".vcf.gz"),
                    true,
                    Some(&starts[i]),
                    None,
                    None,
                )?;
                *collectors[i].lock().unwrap() = collector;
                *counts[i].lock().unwrap() = traversal;
                for (total, elapsed) in phase_totals.lock().unwrap()[p.group].iter_mut().zip(phases)
                {
                    *total += elapsed;
                }
                for (phase, time) in s::timing::PHASES.iter().zip(phases) {
                    timing.add(
                        &names[p.group],
                        &format!("caller.{}.{phase}", p.contig),
                        time.as_secs_f64(),
                    );
                }
            } else {
                let bam = if kind == 1 {
                    &chain.tumour
                } else {
                    &chain.normal
                };
                let (stats, state) = s::get_pileup_summaries::run_part(
                    bam,
                    &chain.pileup_resource,
                    &chain.reference,
                    &beds[p.group],
                    &suffix(&p.prefix, &format!(".{kind}.table")),
                    Some(p.contig),
                )?;
                *pileup_counts[i * 2 + kind - 1].lock().unwrap() = stats;
                *af[i * 2 + kind - 1].lock().unwrap() = Some(state);
            }
            Ok(())
        })?;
        let end = Instant::now();
        let mut spans = spans.lock().unwrap();
        let span = &mut spans[p.group][kind];
        *span = Some(span.map_or((start, end), |(a, b)| (a.min(start), b.max(end))));
        Ok(())
    })?;
    for (g, group) in groups.iter().enumerate() {
        let start_merge = Instant::now();
        combine(
            &group
                .clone()
                .map(|i| suffix(&parts[i].prefix, ".vcf.gz"))
                .collect::<Vec<_>>(),
            &vcfs[g],
            true,
        )?;
        let mut traversal = s::mutect2_regions::Stats::default();
        for i in group.clone() {
            let part = counts[i].lock().unwrap();
            traversal.processed += part.processed;
            traversal.callable += part.callable;
            traversal.regions += part.regions;
            traversal.active += part.active;
            for (total, n) in traversal.filtered.iter_mut().zip(part.filtered) {
                *total += n;
            }
        }
        eprintln!("chain: calling and collecting pileups for {}", names[g]);
        eprint!("{}", traversal.filter_summary());
        traversal.write_callable(&stats[g])?;
        let mut collector = s::mutect2_f1r2::Collector::default();
        for i in group.clone() {
            collector.append(collectors[i].lock().unwrap().take().unwrap())?;
        }
        collector.write_archive(&archives[g])?;
        for (kind, output) in [(1, &tp[g]), (2, &np[g])] {
            let states: Vec<_> = group
                .clone()
                .map(|i| af[i * 2 + kind - 1].lock().unwrap().take().unwrap())
                .collect();
            s::get_pileup_summaries::finish_af(&states)?;
            let mut total = s::get_pileup_summaries::Stats::default();
            for i in group.clone() {
                let part = pileup_counts[i * 2 + kind - 1].lock().unwrap();
                total.processed += part.processed;
                total.rows += part.rows;
                for (sum, n) in total.filtered.iter_mut().zip(part.filtered) {
                    *sum += n;
                }
            }
            eprintln!("chain: pileup {}.{}", labels[kind - 1], names[g]);
            total.log();
            combine(
                &group
                    .clone()
                    .map(|i| suffix(&parts[i].prefix, &format!(".{kind}.table")))
                    .collect::<Vec<_>>(),
                output,
                false,
            )?;
        }
        timing.add(
            &names[g],
            "contig_merge",
            start_merge.elapsed().as_secs_f64(),
        );
        let spans = spans.lock().unwrap()[g];
        for (kind, step) in ["caller", "pileup_tumour", "pileup_normal"]
            .iter()
            .enumerate()
        {
            let (a, b) = spans[kind].unwrap();
            timing.add(&names[g], step, b.duration_since(a).as_secs_f64());
        }
        let a = spans.iter().flatten().map(|x| x.0).min().unwrap();
        let b = spans.iter().flatten().map(|x| x.1).max().unwrap();
        timing.add(&names[g], "group", b.duration_since(a).as_secs_f64());
        for (phase, elapsed) in s::timing::PHASES
            .iter()
            .zip(phase_totals.lock().unwrap()[g])
        {
            timing.add(&names[g], &format!("caller.{phase}"), elapsed.as_secs_f64());
        }
        if self_check {
            // The serial caller prints its own read-filter summary; give it a label so
            // it does not run on under the last pileup block.
            eprintln!("chain: self-check, serial run of {}", names[g]);
            let vcf = scratch.join(format!("serial-{g}.vcf.gz"));
            let stat = scratch.join(format!("serial-{g}.stats"));
            let archive = scratch.join(format!("serial-{g}.tar.gz"));
            let mut serial = input(g);
            serial.stats = &stat;
            s::mutect2_calls::vcf::run_part(
                &serial,
                &vcf,
                true,
                None,
                Some(&starts[group.clone()]),
                Some(&archive),
            )?;
            for (a, b) in [
                (&vcf, &vcfs[g]),
                (&stat, &stats[g]),
                (&archive, &archives[g]),
            ] {
                identical(a, b)?;
            }
            for (bam, output) in [(&chain.tumour, &tp[g]), (&chain.normal, &np[g])] {
                let table = scratch.join("serial.table");
                s::get_pileup_summaries::run(
                    bam,
                    &chain.pileup_resource,
                    &chain.reference,
                    &beds[g],
                    &table,
                )?;
                identical(&table, output)?;
            }
            eprintln!(
                "contig_self_check\t{}\tPASS\tcontigs={}\tcalls={}",
                names[g],
                group.len(),
                starts[group.clone()].iter().map(|s| s.calls).sum::<u64>()
            );
        }
    }
    fs::remove_dir_all(scratch)?;
    Ok(())
}
