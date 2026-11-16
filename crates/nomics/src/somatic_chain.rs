//! Transactional orchestration of the existing gated somatic tools.
mod contigs;
use clap::Args;
use nomics_io::{bam::Header, bam_query::BamQuery};
use nomics_somatic as s;
use std::{
    collections::BTreeSet,
    fs,
    io::{self, BufRead, BufReader, Read, Write},
    path::{Path, PathBuf},
    sync::Mutex,
    time::Instant,
};

#[derive(Args, Debug)]
pub struct Chain {
    #[arg(long)]
    pub tumour: PathBuf,
    #[arg(long)]
    pub normal: PathBuf,
    #[arg(long)]
    pub reference: PathBuf,
    #[arg(
        long,
        required_unless_present = "intervals_dir",
        conflicts_with = "intervals_dir"
    )]
    pub intervals: Vec<PathBuf>,
    #[arg(long)]
    pub intervals_dir: Option<PathBuf>,
    #[arg(long)]
    pub panel_of_normals: PathBuf,
    #[arg(long)]
    pub germline_resource: PathBuf,
    /// GetPileupSummaries resource; sarek uses af-only-gnomad.hg38.vcf.gz.
    #[arg(long)]
    pub pileup_resource: PathBuf,
    #[arg(long)]
    pub outdir: PathBuf,
    /// Publication label only; never overrides BAM SM.
    #[arg(long)]
    pub tumour_label: Option<String>,
    #[arg(long)]
    pub normal_label: Option<String>,
    /// Order in which LearnReadOrientationModel reads the groups' f1r2 archives: a complete
    /// permutation of BED stems. Default: the interval-group order. sarek passes them in the
    /// order its Mutect2 tasks finished; give that order to reproduce its priors table byte
    /// for byte. The order moves the priors in their last digits only.
    #[arg(long, value_delimiter = ',')]
    pub orientation_order: Vec<String>,
    #[arg(long)]
    pub keep_intermediates: bool,
    /// CPU task budget (defaults to the available cores).
    #[arg(long, default_value_t = available_threads())]
    pub threads: usize,
}
fn available_threads() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}
fn tasks(
    count: usize,
    threads: usize,
    task: impl Fn(usize) -> io::Result<()> + Sync,
) -> io::Result<()> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        let mut workers = Vec::new();
        for _ in 0..threads.min(count) {
            let (next, task) = (&next, &task);
            workers.push(
                std::thread::Builder::new()
                    .stack_size(8 * 1024 * 1024)
                    .spawn_scoped(scope, move || {
                        let mut results = Vec::new();
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            if i >= count {
                                break;
                            }
                            results.push((i, task(i)));
                        }
                        results
                    })?,
            );
        }
        // Join every worker before propagating errors or cleaning staging.
        let mut results = Vec::new();
        let mut panicked = false;
        for worker in workers {
            match worker.join() {
                Ok(rows) => results.extend(rows),
                Err(_) => panicked = true,
            }
        }
        if panicked {
            return Err(bad("chain task panicked"));
        }
        results.sort_by_key(|r| r.0);
        for (_, result) in results {
            result?;
        }
        Ok(())
    })
}
fn bad(s: impl Into<String>) -> io::Error {
    io::Error::other(s.into())
}
fn suffix(p: &Path, s: &str) -> PathBuf {
    let mut p = p.as_os_str().to_owned();
    p.push(s);
    p.into()
}
fn safe(s: &str) -> io::Result<()> {
    if s.is_empty()
        || s == "."
        || s == ".."
        || !s
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
    {
        return Err(bad("unsafe or empty sample/group label"));
    }
    Ok(())
}
fn bam(p: &Path) -> io::Result<Header> {
    let cram = || {
        bad(
            "BAM is required; a sarek run must use --save_output_as_bam (CRAM input is unsupported)",
        )
    };
    if p.extension()
        .is_some_and(|s| s.eq_ignore_ascii_case("cram"))
    {
        return Err(cram());
    }
    let mut magic = [0; 4];
    fs::File::open(p)?.read_exact(&mut magic)?;
    if &magic == b"CRAM" {
        return Err(cram());
    }
    Ok(BamQuery::open(p, &[])?.header)
}
fn sample(h: &Header) -> io::Result<String> {
    let text = std::str::from_utf8(&h.text).map_err(|_| bad("non-UTF8 BAM header"))?;
    let mut names = BTreeSet::new();
    for line in text.lines().filter(|l| l.starts_with("@RG\t")) {
        names.insert(
            line.split('\t')
                .find_map(|f| f.strip_prefix("SM:"))
                .ok_or_else(|| bad("read group without SM"))?,
        );
    }
    if names.len() != 1 || names.contains("") {
        return Err(bad("each BAM must have exactly one distinct read-group SM"));
    }
    Ok(names.into_iter().next().unwrap().into())
}
fn permutation(names: &[String], order: &[String]) -> io::Result<Vec<usize>> {
    if order.len() != names.len()
        || order.iter().collect::<BTreeSet<_>>() != names.iter().collect::<BTreeSet<_>>()
    {
        return Err(bad(
            "orientation-order must contain every BED stem exactly once",
        ));
    }
    Ok(order
        .iter()
        .map(|s| names.iter().position(|n| n == s).unwrap())
        .collect())
}
fn bed(p: &Path, h: &Header) -> io::Result<Vec<(usize, u32, u32)>> {
    let mut rows = Vec::new();
    for line in fs::read_to_string(p)?
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
    {
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() < 3 {
            return Err(bad(format!("{}: invalid BED record", p.display())));
        }
        let tid = h
            .refs
            .iter()
            .position(|r| r.name == f[0].as_bytes())
            .ok_or_else(|| bad("BED contig absent from BAM"))?;
        let start = f[1].parse::<u32>().map_err(|_| bad("invalid BED start"))?;
        let end = f[2].parse::<u32>().map_err(|_| bad("invalid BED end"))?;
        if end <= start || end as u64 > h.refs[tid].len as u64 || end > 1 << 29 {
            return Err(bad(
                "BED end must be after start and within the contig and BAI range",
            ));
        }
        rows.push((tid, start, end));
    }
    if rows.is_empty() {
        return Err(bad("empty interval group"));
    }
    Ok(rows)
}
fn resource(p: &Path, h: &Header, needs_af: bool) -> io::Result<()> {
    nomics_io::tabix::Tabix::open(p)?;
    let mut reader = BufReader::new(nomics_io::gzip::GzReader::open(p)?);
    let mut line = String::new();
    let mut af = false;
    let mut contigs = 0;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(bad("resource lacks #CHROM header"));
        }
        if line.starts_with("##INFO=<ID=AF,") {
            af = true;
        }
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let name = rest.split([',', '>']).next().unwrap_or("");
            let r = h
                .refs
                .iter()
                .find(|r| r.name == name.as_bytes())
                .ok_or_else(|| bad(format!("resource contig {name} absent from BAM")))?;
            if let Some(len) = rest
                .split([',', '>'])
                .find_map(|f| f.strip_prefix("length="))
                && len.parse::<u64>().ok() != Some(r.len as u64)
            {
                return Err(bad("resource contig length differs from BAM"));
            }
            contigs += 1;
        }
        if line.starts_with("#CHROM\t") {
            break;
        }
        if !line.starts_with('#') {
            return Err(bad("invalid resource VCF header"));
        }
    }
    if contigs == 0 || (needs_af && !af) {
        return Err(bad(
            "resource requires contig headers and germline/pileup resource requires INFO/AF",
        ));
    }
    Ok(())
}
fn layout(root: &Path, tumour: &str, normal: &str) -> (PathBuf, String) {
    let pair = format!("{tumour}_vs_{normal}");
    (
        root.join("variant_calling/mutect2").join(&pair),
        format!("{pair}.mutect2"),
    )
}
#[derive(Default)]
struct Timings(Mutex<Vec<(String, String, f64)>>);
impl Timings {
    fn add(&self, group: &str, step: &str, seconds: f64) {
        eprintln!("chain_timing\t{group}\t{step}\t{seconds:.6}");
        self.0
            .lock()
            .unwrap()
            .push((group.into(), step.into(), seconds));
    }
    fn measure<T>(
        &self,
        group: &str,
        step: &str,
        f: impl FnOnce() -> io::Result<T>,
    ) -> io::Result<T> {
        let start = Instant::now();
        let result = f();
        self.add(group, step, start.elapsed().as_secs_f64());
        result
    }
    fn write(&self, path: &Path) -> io::Result<()> {
        let mut rows = self.0.lock().unwrap();
        rows.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
        let mut out = fs::File::create(path)?;
        writeln!(out, "group\tstep\twall_seconds")?;
        for (group, step, seconds) in rows.iter() {
            writeln!(out, "{group}\t{step}\t{seconds:.6}")?;
        }
        Ok(())
    }
}
impl Chain {
    pub fn run(&self) -> io::Result<()> {
        if self.threads == 0 {
            return Err(bad("threads must be positive"));
        }
        // All configuration validation precedes creation of any output directory.
        let nh = bam(&self.normal)?;
        let th = bam(&self.tumour)?;
        let normal = sample(&nh)?;
        let tumour = sample(&th)?;
        if normal == tumour
            || nh
                .refs
                .iter()
                .map(|r| (&r.name, r.len))
                .ne(th.refs.iter().map(|r| (&r.name, r.len)))
        {
            return Err(bad(
                "BAM samples must differ and sequence dictionaries must agree",
            ));
        }
        let nl = self.normal_label.as_deref().unwrap_or(&normal);
        let tl = self.tumour_label.as_deref().unwrap_or(&tumour);
        safe(nl)?;
        safe(tl)?;
        if nl == tl {
            return Err(bad("publication labels must differ"));
        }
        let reference_file = fs::File::open(&self.reference)?;
        let reference_metadata = reference_file.metadata()?;
        if !reference_metadata.is_file() {
            return Err(bad("reference must be a regular FASTA file"));
        }
        let fai = nomics_io::fasta::read_fai(&suffix(&self.reference, ".fai"))?;
        if fai
            .iter()
            .map(|r| (r.name.as_bytes(), r.len))
            .ne(nh.refs.iter().map(|r| (r.name.as_slice(), r.len as u64)))
        {
            return Err(bad("FASTA/BAM dictionaries differ"));
        }
        for entry in &fai {
            let last = entry.len.saturating_sub(1);
            let end = (last / entry.line_bases)
                .checked_mul(entry.line_width)
                .and_then(|x| entry.offset.checked_add(x))
                .and_then(|x| x.checked_add(last % entry.line_bases + 1));
            if entry.len == 0 || end.is_none_or(|end| end > reference_metadata.len()) {
                return Err(bad("FASTA index extends beyond reference file"));
            }
        }
        let dict = self.reference.with_extension("dict");
        let dt = fs::read_to_string(&dict)?;
        let dn: Vec<_> = dt
            .lines()
            .filter(|l| l.starts_with("@SQ\t"))
            .map(|l| {
                let name = l
                    .split('\t')
                    .find_map(|f| f.strip_prefix("SN:"))
                    .ok_or_else(|| bad("dictionary missing SN"))?;
                let len = l
                    .split('\t')
                    .find_map(|f| f.strip_prefix("LN:"))
                    .and_then(|v| v.parse::<u64>().ok())
                    .ok_or_else(|| bad("dictionary missing LN"))?;
                Ok((name.as_bytes(), len))
            })
            .collect::<io::Result<_>>()?;
        if dn
            .into_iter()
            .ne(nh.refs.iter().map(|r| (r.name.as_slice(), r.len as u64)))
        {
            return Err(bad("sequence dictionary differs from BAM"));
        }
        for p in [
            &self.panel_of_normals,
            &self.germline_resource,
            &self.pileup_resource,
        ] {
            resource(p, &nh, p != &self.panel_of_normals)?;
        }
        let mut beds = self.intervals.clone();
        if let Some(dir) = &self.intervals_dir {
            beds = fs::read_dir(dir)?
                .map(|e| e.map(|e| e.path()))
                .collect::<io::Result<Vec<_>>>()?;
            beds.retain(|p| p.extension().is_some_and(|x| x == "bed"));
            beds.sort();
        }
        if beds.is_empty() {
            return Err(bad("no interval groups"));
        }
        let names: Vec<String> = beds
            .iter()
            .map(|p| {
                p.file_stem()
                    .and_then(|x| x.to_str())
                    .unwrap_or("")
                    .to_owned()
            })
            .collect();
        for n in &names {
            safe(n)?;
        }
        if names.iter().collect::<BTreeSet<_>>().len() != names.len() {
            return Err(bad("duplicate BED stems"));
        }
        let order = if self.orientation_order.is_empty() {
            eprintln!(
                "chain: no --orientation-order; LearnReadOrientationModel reads the f1r2 archives in interval-group order"
            );
            (0..names.len()).collect()
        } else {
            permutation(&names, &self.orientation_order)?
        };
        let mut spans = Vec::new();
        for (g, p) in beds.iter().enumerate() {
            for (t, b, e) in bed(p, &nh)? {
                spans.push((t, b, e, g));
            }
        }
        spans.sort();
        for (i, a) in spans.iter().enumerate() {
            for b in &spans[i + 1..] {
                if b.0 != a.0 || b.1 >= a.2 {
                    break;
                }
                if a.3 != b.3 {
                    return Err(bad("overlapping scatter groups are unsupported"));
                }
            }
        }
        if self.outdir.symlink_metadata().is_ok() {
            return Err(bad("outdir must not exist (refusing overwrite)"));
        }
        let parent = self
            .outdir
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        if !parent.is_dir() {
            return Err(bad("outdir parent must exist"));
        }
        let staging = parent.join(format!(".nomics-chain-{}.incomplete", std::process::id()));
        fs::create_dir(&staging)?;
        let result = self
            .execute(&staging, &beds, &names, &order, &normal, tl, nl)
            .and_then(|()| fs::rename(&staging, &self.outdir));
        if result.is_err() {
            let _ = fs::remove_dir_all(&staging);
        }
        result
    }
    #[allow(clippy::too_many_arguments)]
    fn execute(
        &self,
        root: &Path,
        beds: &[PathBuf],
        names: &[String],
        order: &[usize],
        normal: &str,
        tl: &str,
        nl: &str,
    ) -> io::Result<()> {
        let timing = Timings::default();
        let total = Instant::now();
        let base = root.join("variant_calling/mutect2");
        let (out, p) = layout(root, tl, nl);
        let work = root.join("work");
        for d in [&out, &work, &base.join(tl), &base.join(nl)] {
            fs::create_dir_all(d)?;
        }
        let inputs = vec![self.normal.clone(), self.tumour.clone()];
        // Construct reduction inputs in interval order, never completion order.
        let vcfs: Vec<_> = names
            .iter()
            .map(|g| work.join(format!("{p}.{g}.vcf.gz")))
            .collect();
        let stats: Vec<_> = vcfs.iter().map(|v| suffix(v, ".stats")).collect();
        let archives: Vec<_> = names
            .iter()
            .map(|g| work.join(format!("{p}.{g}.f1r2.tar.gz")))
            .collect();
        let tp: Vec<_> = names
            .iter()
            .map(|g| work.join(format!("{tl}.mutect2.{g}.pileups.table")))
            .collect();
        let np: Vec<_> = names
            .iter()
            .map(|g| work.join(format!("{nl}.mutect2.{g}.pileups.table")))
            .collect();
        contigs::run(
            self,
            &work,
            beds,
            names,
            normal,
            &inputs,
            &vcfs,
            &stats,
            &archives,
            &tp,
            &np,
            [tl, nl],
            &timing,
        )?;
        let raw = out.join(format!("{p}.vcf.gz"));
        let stat = suffix(&raw, ".stats");
        let priors = out.join(format!("{p}.artifactprior.tar.gz"));
        let dict = self.reference.with_extension("dict");
        eprintln!("chain: merging calls, learning orientation and estimating contamination");
        timing.measure("all", "merge_vcfs", || {
            s::merge::merge_vcfs(&vcfs, Some(&dict), &raw)
        })?;
        timing.measure("all", "merge_stats", || {
            s::merge::merge_mutect_stats(&stats, &stat)
        })?;
        timing.measure("all", "orientation", || {
            s::orientation::run(
                &order
                    .iter()
                    .map(|&i| archives[i].clone())
                    .collect::<Vec<_>>(),
                &priors,
            )
        })?;
        let tumor_table = base.join(tl).join(format!("{tl}.mutect2.pileups.table"));
        let normal_table = base.join(nl).join(format!("{nl}.mutect2.pileups.table"));
        timing.measure("all", "gather_tumour", || {
            s::pileup_summary::gather(&tp, &dict, &tumor_table)
        })?;
        timing.measure("all", "gather_normal", || {
            s::pileup_summary::gather(&np, &dict, &normal_table)
        })?;
        let contamination = out.join(format!("{p}.contamination.table"));
        let segmentation = out.join(format!("{p}.segmentation.table"));
        timing.measure("all", "contamination", || {
            s::contamination::calculate(
                &tumor_table,
                Some(&normal_table),
                &contamination,
                Some(&segmentation),
                0.5,
                3.0,
            )
        })?;
        eprintln!("chain: filtering calls");
        timing.measure("all", "filter", || {
            s::filtering::run(&s::filtering::Inputs {
                variant: &raw,
                output: &out.join(format!("{p}.filtered.vcf.gz")),
                reference: &self.reference,
                orientation_priors: &priors,
                segmentation: &segmentation,
                contamination: &contamination,
                stats: Some(&stat),
            })
        })?;
        timing.add("all", "total", total.elapsed().as_secs_f64());
        timing.write(&root.join("chain-timings.tsv"))?;
        if !self.keep_intermediates {
            fs::remove_dir_all(work)?;
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[derive(clap::Parser)]
    struct TestCli {
        #[command(flatten)]
        chain: Chain,
    }
    #[test]
    fn arguments_and_no_writes() {
        use clap::Parser;
        let root = std::env::temp_dir().join(format!("chain-args-{}", std::process::id()));
        let out = root.join("output");
        let args = [
            "chain",
            "--tumour",
            "t.bam",
            "--normal",
            "absent.cram",
            "--reference",
            "ref.fa",
            "--panel-of-normals",
            "pon.gz",
            "--germline-resource",
            "germ.gz",
            "--pileup-resource",
            "sites.gz",
            "--outdir",
            out.to_str().unwrap(),
            "--intervals",
            "g.bed",
            "--orientation-order",
            "g",
        ];
        let cli = TestCli::try_parse_from(args).unwrap();
        assert!(
            cli.chain
                .run()
                .unwrap_err()
                .to_string()
                .contains("--save_output_as_bam")
        );
        assert!(!root.exists());
        let mut conflict = args.to_vec();
        conflict.extend(["--intervals-dir", "beds"]);
        assert!(TestCli::try_parse_from(conflict).is_err());
        // The orientation order is optional; without it the group order is used.
        let cli = TestCli::try_parse_from(&args[..args.len() - 2]).unwrap();
        assert!(cli.chain.orientation_order.is_empty());
    }
    #[test]
    fn invalid_files() {
        let dir = std::env::temp_dir().join(format!("chain-validation-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let p = dir.join("disguised.bam");
        fs::write(&p, b"CRAMjunk").unwrap();
        assert!(bam(&p).unwrap_err().to_string().contains("BAM is required"));
        let h = Header {
            text: Vec::new(),
            refs: vec![nomics_io::bam::Reference {
                name: b"chr1".to_vec(),
                len: 100,
            }],
        };
        let p = dir.join("bad.bed");
        fs::write(&p, "chr1\t20\t10\n").unwrap();
        assert!(
            bed(&p, &h)
                .unwrap_err()
                .to_string()
                .contains("end must be after start")
        );
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn task_budget_and_errors_join() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let active = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let done = AtomicUsize::new(0);
        let result = tasks(17, 2, |i| {
            let n = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(n, Ordering::SeqCst);
            std::thread::yield_now();
            active.fetch_sub(1, Ordering::SeqCst);
            done.fetch_add(1, Ordering::SeqCst);
            if i == 3 {
                Err(bad("intentional"))
            } else {
                Ok(())
            }
        });
        assert!(result.is_err());
        assert_eq!(done.load(Ordering::SeqCst), 17);
        assert_eq!(active.load(Ordering::SeqCst), 0);
        assert!(peak.load(Ordering::SeqCst) <= 2);
    }
    #[test]
    fn orders() {
        let n = vec!["a".into(), "b".into()];
        assert_eq!(
            permutation(&n, &["b".into(), "a".into()]).unwrap(),
            vec![1, 0]
        );
        assert!(permutation(&n, &["a".into(), "a".into()]).is_err());
    }
    #[test]
    fn read_group_samples() {
        let mut h = Header {
            text: b"@RG\tID:a\tSM:patient_sample\n@RG\tID:b\tSM:patient_sample\n".to_vec(),
            refs: vec![],
        };
        assert_eq!(sample(&h).unwrap(), "patient_sample");
        h.text.extend_from_slice(b"@RG\tID:c\tSM:other\n");
        assert!(sample(&h).is_err());
        h.text = b"@RG\tID:a\n".to_vec();
        assert!(sample(&h).is_err());
        h.text = b"@RG\tID:a\tSM:\n".to_vec();
        assert!(sample(&h).is_err());
    }
    #[test]
    fn published_layout() {
        let (out, p) = layout(Path::new("results"), "T", "N");
        assert_eq!(out, Path::new("results/variant_calling/mutect2/T_vs_N"));
        assert_eq!(p, "T_vs_N.mutect2");
        assert_eq!(
            out.join(format!("{p}.filtered.vcf.gz")),
            Path::new("results/variant_calling/mutect2/T_vs_N/T_vs_N.mutect2.filtered.vcf.gz")
        );
    }
    #[test]
    fn labels() {
        for s in ["", "..", "a/b"] {
            assert!(safe(s).is_err());
        }
        safe("HCC1395_WES").unwrap();
    }
    #[test]
    fn cram_extension() {
        assert!(
            bam(Path::new("absent.CRAM"))
                .unwrap_err()
                .to_string()
                .contains("--save_output_as_bam")
        );
    }
}
