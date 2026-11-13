//! `nomics verify`: is nomics identical to sarek on my data, on my machine?
//!
//! Takes a deterministic slice of the user's reads, runs pinned sarek 3.10.0
//! (from containers) and nomics on it, gates every artifact under the
//! `docs/GATES.md` normalizations, and writes `verify.json` and `verify.md`.
//! `--against` instead checks a nomics run on a public fixture against a
//! shipped digest manifest, with no sarek. See `docs/VERIFY.md`.

mod descriptor;
#[cfg(test)]
mod fault_tests;
mod gate;
mod host;
mod report;
mod runner;
mod slice;

use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::time::Instant;

use clap::{Args, Subcommand, ValueEnum};

use descriptor::Descriptor;
use gate::{ArtifactResult, Roots, Status};
use report::{Json, Verdict};

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

pub fn sha256_hex(bytes: &[u8]) -> String {
    let mut h = nomics_gate::Sha256::new();
    h.update(bytes);
    hex(&h.finalize())
}

/// SHA-256 over each `@SQ` line's `SN`, `LN` and `M5` of the iGenomes
/// GATK.GRCh38 `Homo_sapiens_assembly38.dict` (3,366 contigs, ALT included).
pub const GRCH38_CONTIGS_SHA256: &str = "d9e37f99db847b3f5aa7a985195bc2faeca705e3195c00def061223ea2649935";

/// The commit this binary was built from (`build.rs`), or `unknown`.
pub const COMMIT: &str = match option_env!("NOMICS_COMMIT") {
    Some(c) => c,
    None => "unknown",
};

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Device {
    Auto,
    Cpu,
    Cuda,
}

#[derive(Args)]
#[command(args_conflicts_with_subcommands = true)]
pub struct VerifyCli {
    #[command(subcommand)]
    action: Option<Action>,
    #[command(flatten)]
    args: VerifyArgs,
}

#[derive(Subcommand)]
enum Action {
    /// Pull the pinned Nextflow image and sarek source, Nextflow's plugins and
    /// (on arm64) build the GATK image, so later runs need no network.
    Fetch {
        #[arg(long)]
        cache_dir: Option<PathBuf>,
    },
    /// Check that a verify.json has not been edited since it was written.
    CheckReport { report: PathBuf },
}

#[derive(Args)]
pub struct VerifyArgs {
    /// R1 FASTQ (gzip or BGZF). One lane.
    #[arg(long)]
    fastq1: Option<PathBuf>,
    /// R2 FASTQ, in the same order as --fastq1.
    #[arg(long)]
    fastq2: Option<PathBuf>,
    /// Reference directory laid out like iGenomes GATK.GRCh38.
    #[arg(long)]
    reference: Option<PathBuf>,
    /// Exome targets (BED), as sarek's `--wes --intervals`. Omit with --wgs, and with an
    /// `--against` manifest whose golden run called on the whole genome.
    #[arg(long)]
    targets: Option<PathBuf>,
    /// Whole-genome input: no targets; calling is limited to --region.
    #[arg(long)]
    wgs: bool,
    /// With --wgs: the contig HaplotypeCaller calls on (sarek `--intervals`).
    #[arg(long, default_value = "chr20")]
    region: String,
    /// `head:N` (the first N pairs) or `contig:NAME` (pairs whose nomics
    /// alignment touches NAME). Default contig:chr20 for exome, head:4000000
    /// with --wgs.
    #[arg(long)]
    slice: Option<String>,
    /// New or empty output directory.
    #[arg(long)]
    outdir: Option<PathBuf>,
    /// nomics threads (default: all cores). sarek keeps the goldens' task caps.
    #[arg(long)]
    threads: Option<usize>,
    #[arg(long, value_enum, default_value_t = Device::Auto)]
    device: Device,
    /// nomics index caches, sarek source and Nextflow home (default ~/.cache/nomics).
    #[arg(long)]
    cache_dir: Option<PathBuf>,
    /// Check nomics against a shipped digest manifest instead of running sarek.
    #[arg(long)]
    against: Option<PathBuf>,
    /// With --against: the golden run's read-group DS when the manifest's
    /// record does not carry it.
    #[arg(long)]
    read_group_ds: Option<String>,
    /// Never use the network; fail if anything pinned is missing.
    #[arg(long, conflicts_with = "online")]
    offline: bool,
    /// Allow network use beyond pinned pulls (arm64: Seqera Wave image builds).
    #[arg(long)]
    online: bool,
    /// x86: also run sarek's stock configuration (Intel GKL) and report its
    /// VCF difference as EXPECTED.
    #[arg(long)]
    also_sarek_default: bool,
    /// Keep the slice, sarek's work/ and nomics' intermediates after a PASS.
    #[arg(long)]
    keep_work: bool,
    /// Reference VCF records needed for a PASS to mean anything.
    #[arg(long, default_value_t = 200)]
    min_records: u64,
}

/// A shipped digest manifest `--against` accepts, with what its golden run used.
struct KnownManifest {
    sha256: &'static str,
    name: &'static str,
    fastq1_md5: &'static str,
    fastq2_md5: &'static str,
    patient_sample: &'static str,
    prefix: &'static str,
    /// The golden run's FASTA path (bwa's `@RG DS`); `None`: pass --read-group-ds.
    read_group_ds: Option<&'static str>,
    /// What the golden run called on.
    intervals: GoldenIntervals,
    /// The golden run's scatter order (`nomics run --scatter-order`), when its
    /// recal BAM or table depend on it.
    scatter_order: Option<&'static str>,
}

/// The intervals a golden run called on.
#[derive(Clone, Copy, PartialEq, Eq)]
enum GoldenIntervals {
    /// The exome targets, given as `--targets`.
    Targets,
    /// The reference's whole-genome calling regions
    /// (`Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed`); no `--targets`.
    WgsCallingRegions,
}

const KNOWN_MANIFESTS: &[KnownManifest] = &[
    KnownManifest {
        sha256: "ab37b9432e3021801bfb83b7972bf7403224fac0cc1807b2fcfe64e36a80e5c0",
        name: "reference/tier0-wes-20261008-r1.digest.tsv",
        fastq1_md5: "ecefd51acf7932c7c1e7ec3e059e4d66",
        fastq2_md5: "0d28435ce16d800d60b0b48ca8d61076",
        patient_sample: "HG002_HG002_WES",
        prefix: "HG002_WES",
        read_group_ds: None,
        intervals: GoldenIntervals::Targets,
        scatter_order: None,
    },
    // 30x HG002 WGS (SRR14724532) on x86 sarek; its VCF is the Java-mode
    // HaplotypeCaller merge (docs/REFERENCE-VARIABILITY.md).
    KnownManifest {
        sha256: "d280eec89680ae77bea4b00e1d8439d658cdb11ce23a65f2e1c6216e3f0d974a",
        name: "reference/wgs30-x86-r1.digest.tsv",
        fastq1_md5: "32334668cf250778a209e0b0777e0d05",
        fastq2_md5: "6bc959a5d2cb2cf5ae07e826b621e523",
        patient_sample: "HG002_HG002_WGS30",
        prefix: "HG002_WGS30",
        read_group_ds: Some("/data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"),
        intervals: GoldenIntervals::WgsCallingRegions,
        scatter_order: Some(include_str!("../../../../reference/wgs30-x86-r1.scatter-order.txt")),
    },
];

const PATIENT: &str = "VERIFY";
const SAMPLE: &str = "SAMPLE";
const LANE: &str = "L1";

fn home_cache() -> PathBuf {
    std::env::var_os("HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp"))
        .join(".cache/nomics")
}

fn absolute(p: &Path) -> Result<PathBuf, String> {
    std::fs::canonicalize(p).map_err(|e| format!("{}: {e}", p.display()))
}

/// What the report knows; filled in as the run goes, so an ERROR report still
/// says how far it got.
#[derive(Default)]
struct Ctx {
    mode: &'static str,
    started: u64,
    host: Option<host::Host>,
    device: String,
    inputs: Vec<(&'static str, Json)>,
    reference_dict_sha256: Option<String>,
    slice: Option<(String, slice::Slice, Option<String>)>,
    reference: Option<Json>,
    manifest: Option<(String, String)>,
    artifacts: Vec<ArtifactResult>,
    sarek_default: Option<Json>,
    timing: Vec<(&'static str, f64)>,
    notes: Vec<String>,
}

impl Ctx {
    fn time<T>(&mut self, step: &'static str, f: impl FnOnce(&mut Ctx) -> Result<T, String>) -> Result<T, String> {
        let t = Instant::now();
        let r = f(self);
        self.timing.push((step, t.elapsed().as_secs_f64()));
        r
    }
}

pub fn main(cli: VerifyCli) -> ExitCode {
    if cfg!(feature = "verify-fault-tiebreak") {
        eprintln!("verify: this is a test-only build with a wrong duplicate tie-break; it must FAIL");
    }
    let mut d = Descriptor::parse(descriptor::SAREK_GERMLINE).expect("the shipped descriptor parses");
    // Nextflow runs from a local build of the pinned image: on arm64 because the
    // image is amd64-only, on x86 to replace its Docker CLI (too old for Engine 29).
    let upstream_nextflow = d.nextflow_image.clone();
    d.nextflow_image = runner::local_nextflow_image(std::env::consts::ARCH, &d.nextflow_version).0;
    if let Some(Action::CheckReport { report }) = &cli.action {
        let text = match std::fs::read_to_string(report) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("verify check-report: {}: {e}", report.display());
                return ExitCode::from(2);
            }
        };
        return match report::check_seal(&text) {
            Some(true) => {
                println!("report_sha256 matches: unedited since it was written");
                ExitCode::SUCCESS
            }
            Some(false) => {
                println!("report_sha256 does NOT match: the report was edited");
                ExitCode::FAILURE
            }
            None => {
                eprintln!("verify check-report: {} has no report_sha256", report.display());
                ExitCode::from(2)
            }
        };
    }
    if let Some(Action::Fetch { cache_dir }) = cli.action {
        let cache = cache_dir.unwrap_or_else(home_cache);
        return match fetch(&d, &upstream_nextflow, &cache) {
            Ok(()) => {
                eprintln!("verify fetch: done; later runs can use --offline");
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("verify fetch: {e}");
                ExitCode::from(2)
            }
        };
    }
    let a = cli.args;
    let outdir = match prepare_outdir(a.outdir.as_deref()) {
        Ok(o) => o,
        Err(e) => {
            eprintln!("verify: {e}");
            return ExitCode::from(2);
        }
    };
    let mut ctx = Ctx {
        mode: if a.against.is_some() { "manifest" } else { "reference-run" },
        started: now(),
        ..Default::default()
    };
    let result = if a.against.is_some() {
        run_against(&d, &a, &outdir, &mut ctx)
    } else {
        run_reference(&d, &a, &outdir, &mut ctx)
    };
    let (verdict, reason) = match result {
        Ok(()) => report::verdict(&ctx.artifacts, if ctx.mode == "manifest" { 1 } else { a.min_records }),
        Err(e) => (Verdict::Error, e),
    };
    match write_report(&d, &ctx, &outdir, verdict, &reason) {
        Ok(hash) => {
            eprintln!("verify: {} ({reason})", verdict.name());
            eprintln!("verify: report {} (sha256 {hash})", outdir.join("verify.md").display());
        }
        Err(e) => {
            eprintln!("verify: could not write the report: {e}");
            return ExitCode::from(2);
        }
    }
    if verdict == Verdict::Pass && !a.keep_work {
        for sub in ["slice", "sarek/work", "nomics/prep", "nomics/mapped"] {
            let _ = std::fs::remove_dir_all(outdir.join(sub));
        }
    }
    ExitCode::from(verdict.exit_code())
}

fn now() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
}

fn prepare_outdir(outdir: Option<&Path>) -> Result<PathBuf, String> {
    let outdir = outdir.ok_or("--outdir is required")?;
    if outdir.exists() && std::fs::read_dir(outdir).map_err(|e| e.to_string())?.next().is_some() {
        return Err(format!("{} is not empty; give a new directory", outdir.display()));
    }
    for sub in ["slice", "sarek", "nomics", "digests", "logs", "inputs"] {
        std::fs::create_dir_all(outdir.join(sub)).map_err(|e| format!("{}: {e}", outdir.display()))?;
    }
    absolute(outdir)
}

/// Contigs of a `.dict` as `(SN, LN, M5)` lines, and their SHA-256.
fn dict_contigs(dict: &Path) -> Result<(Vec<(String, u64)>, String), String> {
    let text = std::fs::read_to_string(dict).map_err(|e| format!("{}: {e}", dict.display()))?;
    let mut h = nomics_gate::Sha256::new();
    let mut contigs = Vec::new();
    for line in text.lines().filter(|l| l.starts_with("@SQ\t")) {
        let field = |k: &str| line.split('\t').find_map(|f| f.strip_prefix(k)).unwrap_or("");
        let (sn, ln, m5) = (field("SN:"), field("LN:"), field("M5:"));
        h.update(format!("{sn}\t{ln}\t{m5}\n").as_bytes());
        contigs.push((sn.to_owned(), ln.parse::<u64>().unwrap_or(0)));
    }
    Ok((contigs, hex(&h.finalize())))
}

/// Reference checks shared by both modes: the iGenomes layout and the pinned
/// GRCh38 contigs.
fn check_reference(d: &Descriptor, reference: &Path, ctx: &mut Ctx) -> Result<Vec<(String, u64)>, String> {
    let dict = reference.join("Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict");
    for p in [
        dict.clone(),
        reference.join("Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"),
        reference.join("Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"),
    ] {
        if !p.exists() {
            return Err(format!("{} is missing; --reference must be an iGenomes GATK.GRCh38 directory", p.display()));
        }
    }
    let (contigs, sha) = dict_contigs(&dict)?;
    ctx.reference_dict_sha256 = Some(sha.clone());
    if sha != GRCH38_CONTIGS_SHA256 {
        return Err(format!("refused: {}", d.refusal("dict")));
    }
    Ok(contigs)
}

/// BED intervals and their contigs, checked against the reference.
fn check_bed(bed: &Path, contigs: &[(String, u64)]) -> Result<u64, String> {
    let text = std::fs::read_to_string(bed).map_err(|e| format!("{}: {e}", bed.display()))?;
    let known: std::collections::HashMap<&str, u64> = contigs.iter().map(|(n, l)| (n.as_str(), *l)).collect();
    let mut n = 0;
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let bad = || format!("{}:{}: not a BED interval", bed.display(), i + 1);
        if f.len() < 3 {
            return Err(bad());
        }
        let (start, end) = (f[1].parse::<u64>().map_err(|_| bad())?, f[2].parse::<u64>().map_err(|_| bad())?);
        let len = known
            .get(f[0])
            .ok_or_else(|| format!("{}:{}: contig {} is not in the reference", bed.display(), i + 1, f[0]))?;
        if start >= end || end > *len {
            return Err(bad());
        }
        n += 1;
    }
    if n == 0 {
        return Err(format!("{}: no intervals", bed.display()));
    }
    Ok(n)
}

fn file_json(path: &Path) -> Result<Json, String> {
    let (md5, bytes) = slice::file_md5(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Json::obj(vec![
        ("basename", Json::str(path.file_name().map_or_else(String::new, |n| n.to_string_lossy().into_owned()))),
        ("bytes", Json::Int(bytes as i64)),
        ("md5", Json::str(md5)),
    ]))
}

fn require<'a>(v: &'a Option<PathBuf>, flag: &str) -> Result<&'a PathBuf, String> {
    v.as_ref().ok_or_else(|| format!("{flag} is required"))
}

fn docker_version() -> Result<String, String> {
    let out = Command::new("docker")
        .args(["version", "--format", "{{.Server.Version}}"])
        .output()
        .map_err(|e| format!("Docker is required and was not found ({e})"))?;
    if !out.status.success() {
        return Err(format!(
            "Docker is installed but not usable: {}",
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    Ok(String::from_utf8_lossy(&out.stdout).trim().to_owned())
}

/// Run a command with its output in a log file; an error names the log.
fn run_logged(mut c: Command, log: &Path, what: &str) -> Result<(), String> {
    let f = std::fs::File::create(log).map_err(|e| format!("{}: {e}", log.display()))?;
    let f2 = f.try_clone().map_err(|e| e.to_string())?;
    let status = c
        .stdin(Stdio::null())
        .stdout(f)
        .stderr(f2)
        .status()
        .map_err(|e| format!("{what}: could not start: {e}"))?;
    if !status.success() {
        return Err(format!("{what} failed ({status}); see {}", log.display()));
    }
    Ok(())
}

fn sarek_checkout(d: &Descriptor, cache: &Path) -> PathBuf {
    cache.join("sarek").join(&d.commit)
}

/// Clone the pinned sarek with Nextflow's own git (no git on the host), and
/// check it is at the pinned commit.
/// Copy the Nextflow framework the pinned image ships into the cached
/// Nextflow home. verify points `NXF_HOME` at the cache so plugins persist
/// between runs; without the framework there, Nextflow would download it
/// again, which fails offline.
fn seed_nextflow_home(d: &Descriptor, cache: &Path) -> Result<(), String> {
    let nxf = cache.join("nextflow");
    // Seeded means the framework jar is there: a failed copy can leave the
    // version directory behind, empty.
    let seeded = || {
        let dir = nxf.join("framework").join(&d.nextflow_version);
        std::fs::read_dir(&dir).is_ok_and(|mut it| {
            it.any(|e| e.is_ok_and(|e| e.file_name().to_string_lossy().ends_with(".jar")))
        })
    };
    if seeded() {
        return Ok(());
    }
    std::fs::create_dir_all(&nxf).map_err(|e| e.to_string())?;
    let (uid, gid) = runner::uid_gid();
    // The image ships the framework jar as root-only (0600), so copy as root
    // and hand the copy to the caller.
    let mut c = Command::new("docker");
    c.args(["run", "--rm", "--entrypoint", "sh"])
        .args(["-v", &format!("{}:/seed", nxf.display())])
        .arg(&d.nextflow_image)
        .args([
            "-c",
            &format!("cp -R \"${{NXF_HOME:-/.nextflow}}/framework\" /seed/ && chown -R {uid}:{gid} /seed/framework"),
        ]);
    run_logged(c, &cache.join("nextflow-seed.log"), "copying the Nextflow framework from its image")?;
    if seeded() {
        Ok(())
    } else {
        Err(format!(
            "{} does not carry Nextflow {}; the descriptor's image and version disagree",
            d.nextflow_image, d.nextflow_version
        ))
    }
}

fn ensure_sarek(d: &Descriptor, cache: &Path, offline: bool) -> Result<PathBuf, String> {
    let dir = sarek_checkout(d, cache);
    if runner::git_head(&dir).as_deref() == Some(d.commit.as_str()) {
        return Ok(dir);
    }
    if offline {
        return Err(format!(
            "the pinned sarek source is not in {}; run `nomics verify fetch` first",
            dir.display()
        ));
    }
    let _ = std::fs::remove_dir_all(&dir);
    let parent = dir.parent().expect("has a parent");
    std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    let nxf = cache.join("nextflow");
    std::fs::create_dir_all(&nxf).map_err(|e| e.to_string())?;
    let mut c = runner::nextflow_container(d, &nxf, &nxf, &[cache.to_path_buf()], true);
    // Nextflow writes .nextflow/ under its working directory; the image's is /.
    c.args(["-w", &nxf.display().to_string()]);
    c.arg(&d.nextflow_image)
        .args(["nextflow", "clone", "-r", &d.source_tag, &d.source_url])
        .arg(&dir);
    run_logged(c, &cache.join("sarek-clone.log"), "cloning sarek")?;
    match runner::git_head(&dir) {
        Some(h) if h == d.commit => Ok(dir),
        h => Err(format!(
            "sarek {} checked out at {}, not the pinned {}",
            d.source_tag,
            h.unwrap_or_else(|| "an unknown commit".into()),
            d.commit
        )),
    }
}

/// Install the Nextflow plugins sarek's `nextflow.config` names, into the
/// cached Nextflow home, so `-offline` runs find them.
fn ensure_plugins(d: &Descriptor, cache: &Path, sarek: &Path) -> Result<(), String> {
    let config = std::fs::read_to_string(sarek.join("nextflow.config")).map_err(|e| e.to_string())?;
    let plugins: Vec<String> = config
        .lines()
        .filter_map(|l| l.trim().strip_prefix("id '"))
        .filter_map(|l| l.split('\'').next())
        .map(str::to_owned)
        .collect();
    if plugins.is_empty() {
        return Ok(());
    }
    let nxf = cache.join("nextflow");
    let mut c = runner::nextflow_container(d, &nxf, &nxf, &[cache.to_path_buf()], true);
    c.args(["-w", &nxf.display().to_string()]);
    c.arg(&d.nextflow_image).args(["nextflow", "plugin", "install", &plugins.join(",")]);
    run_logged(c, &cache.join("nextflow-plugins.log"), "installing Nextflow plugins")
}

/// Pull the task containers the descriptor pins for this architecture, by
/// digest; offline, require them to be present already.
fn ensure_pins(d: &Descriptor, arch: &str, offline: bool, cache: &Path) -> Result<(), String> {
    for pin in d.pins(arch) {
        if runner::image_digest(&pin.image).is_some() {
            continue;
        }
        if offline {
            return Err(format!(
                "{} ({}) is not in Docker; run `nomics verify fetch` first",
                pin.image, pin.process
            ));
        }
        let mut c = Command::new("docker");
        c.args(["pull", &pin.image]);
        run_logged(c, &cache.join("pin-pull.log"), &format!("pulling {}", pin.image))?;
    }
    Ok(())
}

/// Build this architecture's Nextflow image from the pinned upstream image.
fn ensure_local_nextflow(d: &Descriptor, upstream: &str, cache: &Path) -> Result<(), String> {
    if runner::image_digest(&d.nextflow_image).is_some() {
        return Ok(());
    }
    let dockerfile = runner::local_nextflow_image(std::env::consts::ARCH, &d.nextflow_version).1;
    let ctx = cache.join("nextflow-image");
    std::fs::create_dir_all(&ctx).map_err(|e| e.to_string())?;
    std::fs::write(ctx.join("Dockerfile"), dockerfile).map_err(|e| e.to_string())?;
    let mut c = Command::new("docker");
    c.args(["build", "--build-arg", &format!("NXF_IMAGE={upstream}"), "-t", &d.nextflow_image]).arg(&ctx);
    run_logged(c, &cache.join("nextflow-image-build.log"), "building the Nextflow image")
}

/// A local image verify builds in `fetch` must exist before a run.
fn require_local_image(image: &str) -> Result<(), String> {
    match runner::image_digest(image) {
        Some(_) => Ok(()),
        None => Err(format!("{image} is not built; run `nomics verify fetch` first")),
    }
}

fn fetch(d: &Descriptor, upstream_nextflow: &str, cache: &Path) -> Result<(), String> {
    docker_version()?;
    std::fs::create_dir_all(cache).map_err(|e| format!("{}: {e}", cache.display()))?;
    let cache = absolute(cache)?;
    ensure_local_nextflow(d, upstream_nextflow, &cache)?;
    seed_nextflow_home(d, &cache)?;
    let sarek = ensure_sarek(d, &cache, false)?;
    ensure_plugins(d, &cache, &sarek)?;
    ensure_pins(d, std::env::consts::ARCH, false, &cache)?;
    if d.pins(std::env::consts::ARCH).is_empty() {
        eprintln!(
            "verify fetch: task containers are pulled by the first verify run; a second run is then fully offline"
        );
    }
    Ok(())
}

/// The WGS calling regions of the bundle, restricted to one contig.
fn wgs_intervals(reference: &Path, region: &str, out: &Path) -> Result<(), String> {
    let bed = reference.join("Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed");
    let text = std::fs::read_to_string(&bed).map_err(|e| format!("{}: {e}", bed.display()))?;
    let lines: String = text
        .lines()
        .filter(|l| l.split('\t').next() == Some(region))
        .map(|l| format!("{l}\n"))
        .collect();
    if lines.is_empty() {
        return Err(format!("--region {region}: no calling regions in {}", bed.display()));
    }
    std::fs::write(out, lines).map_err(|e| e.to_string())
}

fn nomics_env(c: &mut Command, device: Device) -> Result<String, String> {
    let gpu = cfg!(feature = "cuda") && nomics_gpu::device_info(0).is_ok();
    let use_cuda = match device {
        Device::Cuda if !gpu => {
            return Err("--device cuda: this build or machine has no usable CUDA GPU".to_owned());
        }
        Device::Cuda => true,
        Device::Auto => gpu,
        Device::Cpu => false,
    };
    if use_cuda {
        c.env("NOMICS_GPU_SEED", "1");
        Ok("cuda".to_owned())
    } else {
        c.env_remove("NOMICS_GPU_SEED").env("NOMICS_GPU", "0").env("NOMICS_GPU_EXTEND", "0");
        Ok("cpu".to_owned())
    }
}

struct NomicsRun<'a> {
    fastq1: &'a Path,
    fastq2: &'a Path,
    reference: &'a Path,
    intervals: &'a Path,
    outdir: &'a Path,
    cache: &'a Path,
    scatter_order: Option<&'a Path>,
    read_group_ds: &'a str,
    sample: &'a str,
    prefix: &'a str,
    threads: usize,
    device: Device,
}

fn run_nomics(r: &NomicsRun, log: &Path) -> Result<String, String> {
    let exe = std::env::current_exe().map_err(|e| e.to_string())?;
    let bundle = r.reference.join("Annotation/GATKBundle");
    let dbsnp = bundle.join("dbsnp_146.hg38.vcf.gz");
    // The known-sites cache key includes canonical paths: one per reference.
    let tag = &sha256_hex(r.reference.display().to_string().as_bytes())[..12];
    let mut c = Command::new(exe);
    c.arg("run")
        .arg("--fastq1")
        .arg(r.fastq1)
        .arg("--fastq2")
        .arg(r.fastq2)
        .args(["--lane", LANE])
        .arg("--reference")
        .arg(r.reference)
        .arg("--targets")
        .arg(r.intervals)
        .arg("--outdir")
        .arg(r.outdir)
        .arg("--dbsnp")
        .arg(&dbsnp)
        .arg("--known-sites")
        .arg(&dbsnp)
        .arg("--known-sites")
        .arg(bundle.join("beta/Homo_sapiens_assembly38.known_indels.vcf.gz"))
        .arg("--known-sites")
        .arg(bundle.join("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"))
        .arg("--known-sites-cache")
        .arg(r.cache.join(format!("known_sites.{tag}.nks")))
        .arg("--dense-sa-cache")
        .arg(r.cache.join(format!("sa4.{tag}.nsa")))
        .args(["--dense-sa-interval", "4"])
        .args(["--read-group-ds", r.read_group_ds])
        .args(["--sample", r.sample])
        .args(["--prefix", r.prefix])
        .args(["--threads", &r.threads.to_string()])
        .arg("--write-intermediates");
    if let Some(s) = r.scatter_order {
        c.arg("--scatter-order").arg(s);
    }
    let device = nomics_env(&mut c, r.device)?;
    run_logged(c, log, "nomics run")?;
    Ok(device)
}

fn parse_slice(spec: Option<&str>, wgs: bool) -> Result<(String, SliceKind), String> {
    // Exome: chr20 reads give a full-depth exome on one chromosome (Spark,
    // 30x exome: 58x mean target depth, 1,060 calls), where head:500000 gave
    // 0.8x. Whole genome: a chr20 slice takes ~50 min of HaplotypeCaller, so
    // the default stays a head slice.
    let spec = spec.map_or_else(|| if wgs { "head:4000000".to_owned() } else { "contig:chr20".to_owned() }, str::to_owned);
    match spec.split_once(':') {
        Some(("head", n)) => {
            let n: u64 = n.parse().map_err(|_| format!("--slice {spec}: N must be a positive number"))?;
            if n == 0 {
                return Err(format!("--slice {spec}: N must be positive"));
            }
            Ok((spec.clone(), SliceKind::Head(n)))
        }
        Some(("contig", name)) if !name.is_empty() => Ok((spec.clone(), SliceKind::Contig(name.to_owned()))),
        _ => Err(format!("--slice {spec}: use head:N or contig:NAME")),
    }
}

enum SliceKind {
    Head(u64),
    Contig(String),
}

/// `contig:NAME`: pairs whose nomics alignment touches NAME (`nomics dev
/// slice`); the digest is taken over what it wrote.
#[allow(clippy::too_many_arguments)]
fn contig_slice(fq1: &Path, fq2: &Path, reference: &Path, contig: &str, cache: &Path, threads: usize, device: Device, out: (&Path, &Path), log: &Path) -> Result<slice::Slice, String> {
    let exe = std::env::current_exe().map_err(|e| e.to_string())?;
    let tag = &sha256_hex(reference.display().to_string().as_bytes())[..12];
    let mut c = Command::new(exe);
    c.args(["dev", "slice", "--index"])
        .arg(reference.join("Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64"))
        .arg("--fastq1")
        .arg(fq1)
        .arg("--fastq2")
        .arg(fq2)
        .arg("--out1")
        .arg(out.0)
        .arg("--out2")
        .arg(out.1)
        .args(["--contig", contig])
        .arg("--dense-sa-cache")
        .arg(cache.join(format!("sa4.{tag}.nsa")))
        .args(["--dense-sa-interval", "4", "--threads", &threads.to_string()]);
    // Same device as the nomics run: GPU seeding selects the same pairs
    // (Spark, 30x exome: identical slice, 220 s vs 500 s on the CPU).
    nomics_env(&mut c, device)?;
    run_logged(c, log, "nomics dev slice")?;
    // Digest what it wrote, with the same definition as a head slice.
    let tmp = (out.0.with_extension("check.gz"), out.1.with_extension("check.gz"));
    let s = slice::head(out.0, out.1, u64::MAX, &tmp.0, &tmp.1).map_err(|e| e.to_string())?;
    let _ = std::fs::remove_file(&tmp.0);
    let _ = std::fs::remove_file(&tmp.1);
    Ok(slice::Slice { pairs_requested: s.pairs_taken, ..s })
}

fn run_reference(d: &Descriptor, a: &VerifyArgs, outdir: &Path, ctx: &mut Ctx) -> Result<(), String> {
    let fastq1 = absolute(require(&a.fastq1, "--fastq1")?)?;
    let fastq2 = absolute(require(&a.fastq2, "--fastq2")?)?;
    let reference = absolute(require(&a.reference, "--reference")?)?;
    if a.wgs == a.targets.is_some() {
        return Err("give --targets for an exome, or --wgs for a genome".to_owned());
    }
    let cache = a.cache_dir.clone().unwrap_or_else(home_cache);
    std::fs::create_dir_all(&cache).map_err(|e| format!("{}: {e}", cache.display()))?;
    let cache = absolute(&cache)?;
    let threads = a.threads.unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    let (slice_spec, slice_kind) = parse_slice(a.slice.as_deref(), a.wgs)?;

    // Preflight.
    let host = host::probe();
    ctx.host = Some(host.clone());
    let arch = host.arch.clone();
    let intervals = outdir.join("inputs/intervals.bed");
    ctx.time("preflight", |ctx| {
        let docker = docker_version()?;
        ctx.notes.push(format!("Docker {docker}"));
        if d.profile(&arch).is_none() {
            return Err(format!("no reference configuration for {arch}; verify runs on x86_64 and aarch64 Linux"));
        }
        if a.also_sarek_default && arch != "x86_64" {
            return Err("--also-sarek-default is for x86 (Intel GKL); arm64 has no stock configuration that differs".to_owned());
        }
        let contigs = check_reference(d, &reference, ctx)?;
        if a.wgs {
            if !contigs.iter().any(|(n, _)| *n == a.region) {
                return Err(format!("--region {}: not a reference contig", a.region));
            }
            wgs_intervals(&reference, &a.region, &intervals)?;
        } else {
            let targets = absolute(a.targets.as_ref().unwrap())?;
            std::fs::copy(&targets, &intervals).map_err(|e| format!("{}: {e}", targets.display()))?;
        }
        let n = check_bed(&intervals, &contigs)?;
        ctx.inputs.push(("fastq1", file_json(&fastq1)?));
        ctx.inputs.push(("fastq2", file_json(&fastq2)?));
        let mut bed = match file_json(&intervals)? {
            Json::Obj(p) => p,
            _ => unreachable!(),
        };
        if let Some(t) = &a.targets {
            bed[0].1 = Json::str(t.file_name().map_or_else(String::new, |n| n.to_string_lossy().into_owned()));
        } else {
            bed[0].1 = Json::str(format!("wgs_calling_regions_noseconds.hg38.bed ({})", a.region));
        }
        bed.push(("intervals".to_owned(), Json::Int(n as i64)));
        ctx.inputs.push(("intervals", Json::Obj(bed)));
        Ok(())
    })?;

    // Slice.
    let (s1, s2) = (outdir.join("slice/slice_1.fastq.gz"), outdir.join("slice/slice_2.fastq.gz"));
    let sliced = ctx.time("slice", |_| match &slice_kind {
        SliceKind::Head(n) => slice::head(&fastq1, &fastq2, *n, &s1, &s2).map_err(|e| format!("slice: {e}")),
        SliceKind::Contig(c) => contig_slice(&fastq1, &fastq2, &reference, c, &cache, threads, a.device, (&s1, &s2), &outdir.join("logs/slice.log")),
    })?;
    if sliced.pairs_taken < sliced.pairs_requested {
        ctx.notes.push(format!(
            "the input has {} pairs, fewer than the {} asked for; all were used",
            sliced.pairs_taken, sliced.pairs_requested
        ));
    }
    let region = a.wgs.then(|| a.region.clone());
    ctx.slice = Some((slice_spec, sliced, region));

    // Reference run.
    let sarek = ctx.time("setup", |_| {
        require_local_image(&d.nextflow_image)?;
        ensure_pins(d, &arch, a.offline, &cache)?;
        seed_nextflow_home(d, &cache)?;
        let sarek = ensure_sarek(d, &cache, a.offline)?;
        if !a.offline {
            ensure_plugins(d, &cache, &sarek)?;
        }
        Ok(sarek)
    })?;
    let plan = runner::Plan {
        descriptor: d,
        arch: &arch,
        dir: outdir.join("sarek"),
        sarek: sarek.clone(),
        nxf_home: cache.join("nextflow"),
        reference: reference.clone(),
        fastq1: s1.clone(),
        fastq2: s2.clone(),
        intervals: intervals.clone(),
        wes: !a.wgs,
        patient: PATIENT.to_owned(),
        sample: SAMPLE.to_owned(),
        lane: LANE.to_owned(),
        host_cpus: host.cores,
        host_memory_gib: host.memory_gib,
        online: a.online,
        offline: a.offline,
        extra_mounts: vec![cache.clone()],
    };
    plan.write_inputs().map_err(|e| e.to_string())?;
    let config = plan.config_name().map_err(|e| e.to_string())?;
    let mut reference_json = vec![
        ("pipeline", Json::str(&d.pipeline)),
        ("version", Json::str(&d.version)),
        ("commit", Json::str(&d.commit)),
        ("nextflow", Json::str(&d.nextflow_version)),
        ("nextflow_image", Json::str(&d.nextflow_image)),
        ("nextflow_image_digest", Json::opt_str(runner::image_digest(&d.nextflow_image).as_deref())),
        ("config", Json::str(&config)),
    ];
    let sarek_result = ctx.time("reference", |_| {
        run_logged(plan.command().map_err(|e| e.to_string())?, &outdir.join("logs/sarek.log"), "the sarek reference run")
    });
    let tasks = std::fs::read_to_string(plan.trace())
        .ok()
        .and_then(|t| runner::parse_trace(&t).ok())
        .unwrap_or_default();
    let containers: Vec<Json> = runner::containers(&tasks)
        .into_iter()
        .map(|(process, image)| {
            Json::obj(vec![
                ("process", Json::str(process)),
                ("digest", Json::opt_str(runner::image_digest(&image).as_deref())),
                ("image", Json::str(image)),
            ])
        })
        .collect();
    reference_json.push(("containers", Json::Arr(containers)));
    ctx.reference = Some(Json::obj(reference_json));
    sarek_result?;
    runner::check_tasks(&tasks).map_err(|e| e.to_string())?;
    runner::check_pins(&tasks, &d.pins(&arch)).map_err(|e| e.to_string())?;

    // What nomics takes from the reference run.
    let results = plan.results();
    let md_bam = results.join(format!("preprocessing/markduplicates/{SAMPLE}/{SAMPLE}.md.bam"));
    let rg_ds = runner::read_group_ds(&md_bam)
        .map_err(|e| format!("{}: {e}", md_bam.display()))?
        .ok_or_else(|| format!("{}: no @RG DS", md_bam.display()))?;
    let recal_bam = results.join(format!("preprocessing/recalibrated/{SAMPLE}/{SAMPLE}.recal.bam"));
    let order = runner::scatter_order(&tasks, &plan.work(), recal_bam.exists().then_some(recal_bam.as_path()), SAMPLE)
        .map_err(|e| format!("scatter order: {e}"))?;
    let order_path = outdir.join("inputs/scatter-order.txt");
    if let Some(o) = &order {
        std::fs::write(&order_path, o).map_err(|e| e.to_string())?;
    }

    // nomics on the same slice.
    let nomics_dir = outdir.join("nomics");
    let run = NomicsRun {
        fastq1: &s1,
        fastq2: &s2,
        reference: &reference,
        intervals: &intervals,
        outdir: &nomics_dir,
        cache: &cache,
        scatter_order: order.as_ref().map(|_| order_path.as_path()),
        read_group_ds: &rg_ds,
        sample: &format!("{PATIENT}_{SAMPLE}"),
        prefix: SAMPLE,
        threads,
        device: a.device,
    };
    ctx.device = ctx.time("nomics", |_| run_nomics(&run, &outdir.join("logs/nomics.log")))?;

    // Gate.
    let reference_roots = Roots {
        work: Some(plan.work()),
        results: Some(results),
        run: None,
    };
    let nomics_roots = Roots {
        run: Some(nomics_dir),
        ..Default::default()
    };
    ctx.artifacts = ctx.time("gate", |_| {
        gate::gate_all(&d.artifacts, &reference_roots, &nomics_roots, SAMPLE, LANE).map_err(|e| format!("gate: {e}"))
    })?;
    write_manifests(&ctx.artifacts, outdir)?;

    if a.also_sarek_default {
        ctx.sarek_default = Some(sarek_default(d, &plan, outdir, ctx)?);
    }
    Ok(())
}

/// Write each side's digests in the `nomics dev digest` manifest format.
fn write_manifests(results: &[ArtifactResult], outdir: &Path) -> Result<(), String> {
    for (file, side) in [("reference.tsv", true), ("nomics.tsv", false)] {
        let mut text = String::from("artifact\tsha256\tcount\n");
        for r in results {
            if let Some(d) = if side { &r.reference } else { &r.nomics } {
                text.push_str(&format!("{}\t{}\t{}\n", r.name, d.sha256, d.count));
            }
        }
        std::fs::write(outdir.join("digests").join(file), text).map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// `--also-sarek-default`: sarek's stock x86 configuration on the same
/// slice. Its VCF is expected to differ from the pinned run's (Intel GKL);
/// the count is reported, never gated.
fn sarek_default(d: &Descriptor, pinned: &runner::Plan, outdir: &Path, ctx: &mut Ctx) -> Result<Json, String> {
    let dir = outdir.join("sarek-default");
    std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let plan = runner::Plan {
        dir: dir.clone(),
        descriptor: d,
        arch: pinned.arch,
        sarek: pinned.sarek.clone(),
        nxf_home: pinned.nxf_home.clone(),
        reference: pinned.reference.clone(),
        fastq1: pinned.fastq1.clone(),
        fastq2: pinned.fastq2.clone(),
        intervals: pinned.intervals.clone(),
        wes: pinned.wes,
        patient: pinned.patient.clone(),
        sample: pinned.sample.clone(),
        lane: pinned.lane.clone(),
        host_cpus: pinned.host_cpus,
        host_memory_gib: pinned.host_memory_gib,
        online: pinned.online,
        offline: pinned.offline,
        extra_mounts: pinned.extra_mounts.clone(),
    };
    plan.write_inputs().map_err(|e| e.to_string())?;
    // The stock configuration: drop the Java PairHMM overlay.
    let overlay = plan.overlay().map_err(|e| e.to_string())?;
    let stock = overlay.split(runner::HCJAVA_CONFIG).next().unwrap_or(&overlay).to_owned();
    std::fs::write(dir.join("verify.config"), stock).map_err(|e| e.to_string())?;
    ctx.time("reference_default", |_| {
        run_logged(plan.command().map_err(|e| e.to_string())?, &outdir.join("logs/sarek-default.log"), "the stock sarek run")
    })?;
    let vcf = |p: &runner::Plan| p.results().join(format!("variant_calling/haplotypecaller/{SAMPLE}/{SAMPLE}.haplotypecaller.vcf.gz"));
    let r = nomics_gate::compare_vcfs(&vcf(&plan), &vcf(pinned), 0).map_err(|e| format!("stock VCF: {e}"))?;
    Ok(Json::obj(vec![
        ("status", Json::str("EXPECTED")),
        ("config", Json::str("x86_64 stock: -profile docker, Intel GKL PairHMM")),
        ("vcf_records_pinned", Json::Int(r.records_golden as i64)),
        ("vcf_records_stock", Json::Int(r.records_ours as i64)),
        ("vcf_records_differing", Json::Int(r.mismatches as i64)),
        (
            "note",
            Json::str("sarek's default x86 HaplotypeCaller (Intel GKL) differs from the Java implementation by design; this is not a nomics difference"),
        ),
    ]))
}

fn run_against(d: &Descriptor, a: &VerifyArgs, outdir: &Path, ctx: &mut Ctx) -> Result<(), String> {
    let manifest = absolute(a.against.as_ref().unwrap())?;
    let fastq1 = absolute(require(&a.fastq1, "--fastq1")?)?;
    let fastq2 = absolute(require(&a.fastq2, "--fastq2")?)?;
    let reference = absolute(require(&a.reference, "--reference")?)?;
    let cache = a.cache_dir.clone().unwrap_or_else(home_cache);
    std::fs::create_dir_all(&cache).map_err(|e| format!("{}: {e}", cache.display()))?;
    let cache = absolute(&cache)?;
    let threads = a.threads.unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    ctx.host = Some(host::probe());

    let text = std::fs::read(&manifest).map_err(|e| format!("{}: {e}", manifest.display()))?;
    let sha = sha256_hex(&text);
    let known = KNOWN_MANIFESTS.iter().find(|k| k.sha256 == sha).ok_or_else(|| {
        format!(
            "refused: {} is not a manifest verify knows (known: {})",
            manifest.display(),
            KNOWN_MANIFESTS.iter().map(|k| k.name).collect::<Vec<_>>().join(", ")
        )
    })?;
    ctx.manifest = Some((known.name.to_owned(), sha));
    let targets = match known.intervals {
        GoldenIntervals::Targets => absolute(require(&a.targets, "--targets")?)?,
        GoldenIntervals::WgsCallingRegions => {
            if a.targets.is_some() {
                return Err(format!(
                    "--targets: {} called on the reference's whole-genome calling regions; omit --targets",
                    known.name
                ));
            }
            absolute(&reference.join("Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed"))?
        }
    };
    let scatter_order = match known.scatter_order {
        Some(text) => {
            let p = outdir.join("inputs/scatter-order.txt");
            std::fs::write(&p, text).map_err(|e| format!("{}: {e}", p.display()))?;
            Some(p)
        }
        None => None,
    };
    let golden = nomics_gate::read_digest_manifest(&manifest)?;
    let rg_ds = match (known.read_group_ds, &a.read_group_ds) {
        (_, Some(ds)) => ds.clone(),
        (Some(ds), None) => ds.to_owned(),
        (None, None) => {
            return Err(format!(
                "{}: pass --read-group-ds with the golden run's FASTA path (bwa's @RG DS)",
                known.name
            ));
        }
    };
    ctx.time("preflight", |ctx| {
        let contigs = check_reference(d, &reference, ctx)?;
        let n = check_bed(&targets, &contigs)?;
        for (key, path, want) in [("fastq1", &fastq1, known.fastq1_md5), ("fastq2", &fastq2, known.fastq2_md5)] {
            let f = file_json(path)?;
            let got = match &f {
                Json::Obj(p) => p.iter().find(|(k, _)| k == "md5").map(|(_, v)| v.clone()),
                _ => None,
            };
            if got != Some(Json::str(want)) {
                return Err(format!("--{key}: md5 differs from the FASTQ {} was made from (want {want})", known.name));
            }
            ctx.inputs.push((key, f));
        }
        let mut bed = match file_json(&targets)? {
            Json::Obj(p) => p,
            _ => unreachable!(),
        };
        bed.push(("intervals".to_owned(), Json::Int(n as i64)));
        ctx.inputs.push(("intervals", Json::Obj(bed)));
        Ok(())
    })?;

    let nomics_dir = outdir.join("nomics");
    let run = NomicsRun {
        fastq1: &fastq1,
        fastq2: &fastq2,
        reference: &reference,
        intervals: &targets,
        outdir: &nomics_dir,
        cache: &cache,
        scatter_order: scatter_order.as_deref(),
        read_group_ds: &rg_ds,
        sample: known.patient_sample,
        prefix: known.prefix,
        threads,
        device: a.device,
    };
    ctx.device = ctx.time("nomics", |_| run_nomics(&run, &outdir.join("logs/nomics.log")))?;

    ctx.artifacts = ctx.time("gate", |_| {
        gate_against(d, &golden, &nomics_dir, known.prefix).map_err(|e| format!("gate: {e}"))
    })?;
    write_manifests(&ctx.artifacts, outdir)?;
    Ok(())
}

/// Digest the nomics run per descriptor and compare with the manifest rows.
fn gate_against(
    d: &Descriptor,
    golden: &std::collections::BTreeMap<String, nomics_gate::ArtifactDigest>,
    run: &Path,
    prefix: &str,
) -> std::io::Result<Vec<ArtifactResult>> {
    let roots = Roots {
        run: Some(run.to_path_buf()),
        ..Default::default()
    };
    let mut jobs = Vec::new();
    for spec in &d.artifacts {
        for p in gate::artifact_files(spec, &spec.nomics, &roots, prefix, LANE)?.into_values() {
            jobs.push((spec.comparator, p));
        }
    }
    let digested = gate::Digested::all(jobs);
    let mut out = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for spec in &d.artifacts {
        let ours = gate::artifact_files(spec, &spec.nomics, &roots, prefix, LANE)?;
        let mut names: Vec<String> = ours.keys().cloned().collect();
        let in_group = |n: &str| match spec.group() {
            Some(g) => n.strip_prefix(g).is_some_and(|r| r.starts_with('/')),
            None => n == spec.name,
        };
        names.extend(golden.keys().filter(|n| in_group(n)).cloned());
        names.sort();
        names.dedup();
        for name in names {
            seen.insert(name.clone());
            let nd_result = ours.get(&name).map(|p| digested.take(spec.comparator, p));
            let unreadable = match &nd_result {
                Some(Err(e)) => Some(format!("nomics file unreadable: {e}")),
                _ => None,
            };
            let nd = nd_result.and_then(Result::ok);
            let gd = golden.get(&name).cloned();
            let status = match (&gd, &nd) {
                (Some(_), None) if ours.contains_key(&name) => Status::Unreadable,
                (Some(a), Some(b)) if a.count == 0 || b.count == 0 => Status::Empty,
                (Some(a), Some(b)) if a == b => Status::Pass,
                (Some(_), Some(_)) => Status::Fail,
                _ => Status::Missing,
            };
            out.push(ArtifactResult {
                first_difference: match (status, unreadable) {
                    (_, Some(why)) => Some(gate::Difference { index: None, detail: why }),
                    (Status::Fail, None) => Some(gate::Difference {
                        index: None,
                        detail: "digest differs from the manifest; with no reference files, verify cannot show the record (run without --against to see it)".to_owned(),
                    }),
                    _ => None,
                },
                name,
                comparator: spec.comparator,
                normalization: spec.normalization.clone(),
                reference: gd,
                nomics: nd,
                status,
            });
        }
    }
    // A manifest row no descriptor artifact covers is a harness mismatch.
    if let Some(extra) = golden.keys().find(|k| !seen.contains(*k)) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("manifest artifact {extra} is not in the descriptor"),
        ));
    }
    Ok(out)
}

fn write_report(d: &Descriptor, ctx: &Ctx, outdir: &Path, verdict: Verdict, reason: &str) -> Result<String, String> {
    let host = ctx.host.clone().unwrap_or_else(host::probe);
    let features: Vec<Json> = [
        cfg!(feature = "cuda").then_some("cuda"),
        cfg!(feature = "metal").then_some("metal"),
        cfg!(feature = "mimalloc").then_some("mimalloc"),
        cfg!(feature = "verify-fault-tiebreak").then_some("verify-fault-tiebreak"),
    ]
    .into_iter()
    .flatten()
    .map(Json::str)
    .collect();
    let summary = {
        let n = |s: Status| ctx.artifacts.iter().filter(|r| r.status == s).count() as i64;
        Json::obj(vec![
            ("artifacts", Json::Int(ctx.artifacts.len() as i64)),
            ("pass", Json::Int(n(Status::Pass))),
            ("fail", Json::Int(n(Status::Fail))),
            ("missing", Json::Int(n(Status::Missing))),
            ("empty", Json::Int(n(Status::Empty))),
            ("unreadable", Json::Int(n(Status::Unreadable))),
        ])
    };
    let slice = match &ctx.slice {
        Some((spec, s, region)) => {
            let (method, _) = spec.split_once(':').unwrap_or((spec, ""));
            Json::obj(vec![
                ("method", Json::str(method)),
                ("spec", Json::str(spec)),
                ("pairs_requested", Json::Int(s.pairs_requested.min(i64::MAX as u64) as i64)),
                ("pairs_taken", Json::Int(s.pairs_taken as i64)),
                ("sha256", Json::str(&s.sha256)),
                ("region", Json::opt_str(region.as_deref())),
            ])
        }
        None => Json::Null,
    };
    let mut inputs: Vec<(&str, Json)> = ctx.inputs.clone();
    inputs.push(("reference_contigs_sha256", Json::opt_str(ctx.reference_dict_sha256.as_deref())));
    let mut pairs = vec![
        ("schema", Json::str(report::SCHEMA)),
        ("verdict", Json::str(verdict.name())),
        ("verdict_reason", Json::str(reason)),
        ("created_at", Json::str(report::utc(ctx.started))),
        ("mode", Json::str(ctx.mode)),
        ("descriptor", Json::obj(vec![("id", Json::str(&d.id)), ("sha256", Json::str(&d.sha256))])),
        (
            "tools",
            Json::obj(vec![
                (
                    "nomics",
                    Json::obj(vec![
                        ("version", Json::str(env!("CARGO_PKG_VERSION"))),
                        ("commit", Json::str(COMMIT)),
                        ("features", Json::Arr(features)),
                        ("device", Json::str(if ctx.device.is_empty() { "not run" } else { &ctx.device })),
                    ]),
                ),
                ("reference", ctx.reference.clone().unwrap_or(Json::Null)),
                (
                    "manifest",
                    match &ctx.manifest {
                        Some((name, sha)) => Json::obj(vec![("name", Json::str(name)), ("sha256", Json::str(sha))]),
                        None => Json::Null,
                    },
                ),
            ]),
        ),
        (
            "host",
            Json::obj(vec![
                ("arch", Json::str(&host.arch)),
                ("cpu", Json::str(&host.cpu)),
                ("cores", Json::Int(host.cores as i64)),
                ("memory_gib", Json::Int(host.memory_gib as i64)),
                ("gpu", Json::opt_str(host.gpu.as_deref())),
                ("os", Json::str(&host.os)),
                ("kernel", Json::str(&host.kernel)),
            ]),
        ),
        ("inputs", Json::obj(inputs)),
        ("slice", slice),
        ("artifacts", Json::Arr(ctx.artifacts.iter().map(report::artifact_json).collect())),
        ("summary", summary),
        ("sarek_default", ctx.sarek_default.clone().unwrap_or(Json::Null)),
        (
            "timing_s",
            Json::Obj(ctx.timing.iter().map(|(k, v)| ((*k).to_owned(), Json::Num(*v))).collect()),
        ),
        ("notes", Json::Arr(ctx.notes.iter().map(Json::str).collect())),
    ];
    pairs.push(("report_sha256", Json::str("")));
    let (text, hash) = report::seal(Json::obj(pairs));
    std::fs::write(outdir.join("verify.json"), &text).map_err(|e| e.to_string())?;
    std::fs::write(outdir.join("verify.md"), markdown(d, ctx, &host, verdict, reason, &hash)).map_err(|e| e.to_string())?;
    Ok(hash)
}

fn markdown(d: &Descriptor, ctx: &Ctx, host: &host::Host, verdict: Verdict, reason: &str, hash: &str) -> String {
    let mut m = format!("# nomics verify: {}\n\n{reason}.\n\n| | |\n|---|---|\n", verdict.name());
    let device = if ctx.device.is_empty() { "not run".to_owned() } else { ctx.device.to_uppercase() };
    m.push_str(&format!(
        "| nomics | {} ({}), {device}{} |\n",
        env!("CARGO_PKG_VERSION"),
        &COMMIT[..COMMIT.len().min(12)],
        host.gpu.as_deref().map_or_else(String::new, |g| format!(" on {g}"))
    ));
    if ctx.mode == "manifest" {
        if let Some((name, sha)) = &ctx.manifest {
            m.push_str(&format!("| Compared with | {name} (sha256 {}…) |\n", &sha[..12]));
        }
    } else {
        let config = match &ctx.reference {
            Some(Json::Obj(p)) => p
                .iter()
                .find(|(k, _)| k == "config")
                .and_then(|(_, v)| if let Json::Str(s) = v { Some(s.clone()) } else { None }),
            _ => None,
        };
        m.push_str(&format!(
            "| Reference | {} {} ({}), Nextflow {}{} |\n",
            d.pipeline,
            d.version,
            &d.commit[..8],
            d.nextflow_version,
            config.map_or_else(String::new, |c| format!(", {c}"))
        ));
    }
    m.push_str(&format!(
        "| Host | {}, {}, {} cores, {} GiB, {}, kernel {} |\n",
        host.arch, host.cpu, host.cores, host.memory_gib, host.os, host.kernel
    ));
    let input = |k: &str| {
        ctx.inputs.iter().find(|(n, _)| *n == k).and_then(|(_, v)| match v {
            Json::Obj(p) => {
                let get = |f: &str| p.iter().find(|(n, _)| n == f).map(|(_, v)| v.clone());
                let name = match get("basename") {
                    Some(Json::Str(s)) => s,
                    _ => String::new(),
                };
                let md5 = match get("md5") {
                    Some(Json::Str(s)) => format!(" (md5 {}…)", &s[..8.min(s.len())]),
                    _ => String::new(),
                };
                let n = match get("intervals") {
                    Some(Json::Int(n)) => format!(" ({} intervals)", report::thousands(n as u64)),
                    _ => String::new(),
                };
                Some(format!("{name}{md5}{n}"))
            }
            _ => None,
        })
    };
    let inputs: Vec<String> = ["fastq1", "fastq2", "intervals"].iter().filter_map(|k| input(k)).collect();
    if !inputs.is_empty() {
        m.push_str(&format!("| Input | {} |\n", inputs.join(", ")));
    }
    if let Some((spec, s, region)) = &ctx.slice {
        m.push_str(&format!(
            "| Slice | {spec}: {} pairs{}, sha256 {}… |\n",
            report::thousands(s.pairs_taken),
            region.as_deref().map_or_else(String::new, |r| format!(", calling on {r}")),
            &s.sha256[..12]
        ));
    }
    let t = |k: &str| ctx.timing.iter().find(|(n, _)| *n == k).map(|(_, v)| *v);
    let mins = |s: f64| format!("{} min {:02} s", (s as u64) / 60, (s as u64) % 60);
    let mut times = Vec::new();
    if let Some(s) = t("reference") {
        times.push(format!("sarek {}", mins(s)));
    }
    if let Some(s) = t("nomics") {
        times.push(format!("nomics {}", mins(s)));
    }
    if let Some(s) = t("setup").filter(|s| *s > 5.0) {
        times.push(format!("first-use setup {}", mins(s)));
    }
    if !times.is_empty() {
        m.push_str(&format!("| Time | {} |\n", times.join(", ")));
    }
    if !ctx.artifacts.is_empty() {
        m.push_str("\n| Artifact | Normalization | Records (reference / nomics) | Result |\n|---|---|---|---|\n");
        for row in report::markdown_rows(&ctx.artifacts) {
            m.push_str(&format!("| {} | {} | {} | {} |\n", row[0], row[1], row[2], row[3]));
        }
        for r in ctx.artifacts.iter().filter(|r| r.first_difference.is_some()) {
            let diff = r.first_difference.as_ref().unwrap();
            m.push_str(&format!("\n## {}: first difference\n\n", r.name));
            if let Some(i) = diff.index {
                m.push_str(&format!("Record {i}: "));
            }
            m.push_str(&format!("{}\n", diff.detail));
        }
    }
    if let Some(Json::Obj(p)) = &ctx.sarek_default {
        let get = |k: &str| p.iter().find(|(n, _)| n == k).map_or(0, |(_, v)| if let Json::Int(i) = v { *i } else { 0 });
        m.push_str(&format!(
            "\n## Stock sarek (Intel GKL): EXPECTED\n\n{} of {} VCF records differ from the pinned run; this is GKL against the Java PairHMM, not nomics.\n",
            report::thousands(get("vcf_records_differing") as u64),
            report::thousands(get("vcf_records_pinned") as u64)
        ));
    }
    if !ctx.notes.is_empty() {
        m.push_str("\n## Notes\n\n");
        for n in &ctx.notes {
            m.push_str(&format!("- {n}\n"));
        }
    }
    m.push_str(&format!("\nReport sha256: {hash}\n"));
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_specs() {
        assert!(matches!(parse_slice(None, false).unwrap().1, SliceKind::Contig(ref c) if c == "chr20"));
        assert!(matches!(parse_slice(None, true).unwrap().1, SliceKind::Head(4_000_000)));
        assert!(matches!(parse_slice(Some("contig:chr20"), false).unwrap().1, SliceKind::Contig(ref c) if c == "chr20"));
        for bad in ["head:0", "head:x", "tail:5", "contig:", "500"] {
            assert!(parse_slice(Some(bad), false).is_err(), "{bad}");
        }
    }

    #[test]
    fn dict_digest_and_bed_checks() {
        let dir = std::env::temp_dir().join(format!("nomics-verify-ref-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let dict = dir.join("x.dict");
        std::fs::write(&dict, "@HD\tVN:1.5\n@SQ\tSN:chr1\tLN:100\tM5:abc\tUR:/x\n").unwrap();
        let (contigs, sha) = dict_contigs(&dict).unwrap();
        assert_eq!(contigs, [("chr1".to_owned(), 100)]);
        assert_eq!(sha, sha256_hex(b"chr1\t100\tabc\n"));
        let bed = dir.join("t.bed");
        std::fs::write(&bed, "chr1\t0\t10\nchr1\t20\t30\textra\n").unwrap();
        assert_eq!(check_bed(&bed, &contigs).unwrap(), 2);
        std::fs::write(&bed, "chr2\t0\t10\n").unwrap();
        assert!(check_bed(&bed, &contigs).unwrap_err().contains("not in the reference"));
        std::fs::write(&bed, "chr1\t0\t200\n").unwrap();
        assert!(check_bed(&bed, &contigs).is_err());
        std::fs::write(&bed, "").unwrap();
        assert!(check_bed(&bed, &contigs).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn known_manifests_match_the_shipped_files() {
        for k in KNOWN_MANIFESTS {
            let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..").join(k.name);
            let sha = sha256_hex(&std::fs::read(&path).unwrap());
            assert_eq!(sha, k.sha256, "{}", k.name);
            assert_eq!(k.fastq1_md5.len(), 32, "{}", k.name);
            assert_eq!(k.fastq2_md5.len(), 32, "{}", k.name);
            if let Some(order) = k.scatter_order {
                assert!(order.lines().any(|l| l.starts_with("merge-seed ")), "{}", k.name);
            }
        }
    }

    #[test]
    fn descriptor_covers_the_shipped_manifests() {
        // Every row of every shipped manifest belongs to a descriptor artifact.
        let d = Descriptor::parse(descriptor::SAREK_GERMLINE).unwrap();
        for k in KNOWN_MANIFESTS {
            let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..").join(k.name);
            let m = nomics_gate::read_digest_manifest(&path).unwrap();
            for name in m.keys() {
                assert!(
                    d.artifacts.iter().any(|a| match a.group() {
                        Some(g) => name.starts_with(&format!("{g}/")),
                        None => *name == a.name,
                    }),
                    "{}: {name}",
                    k.name
                );
            }
        }
    }
}
