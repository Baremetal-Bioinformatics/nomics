//! The reference run: pinned sarek under pinned Nextflow, both from
//! containers, so Docker is the only prerequisite. Nextflow runs in its own
//! container with the Docker socket mounted and every path mounted at the
//! same location, so the task containers it starts see the same paths.

use std::collections::BTreeMap;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;

use super::descriptor::Descriptor;

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// Everything the reference run needs to know.
pub struct Plan<'a> {
    pub descriptor: &'a Descriptor,
    /// `x86_64` or `aarch64`.
    pub arch: &'a str,
    /// `<outdir>/sarek`: samplesheet, params, overlay, results, work.
    pub dir: PathBuf,
    /// The checked-out sarek source.
    pub sarek: PathBuf,
    /// Nextflow's home (plugins, assets), kept in the cache between runs.
    pub nxf_home: PathBuf,
    pub reference: PathBuf,
    pub fastq1: PathBuf,
    pub fastq2: PathBuf,
    /// sarek `--intervals`; with `wes` true, also `--wes`.
    pub intervals: PathBuf,
    pub wes: bool,
    pub patient: String,
    pub sample: String,
    pub lane: String,
    pub host_cpus: usize,
    pub host_memory_gib: u64,
    /// Allow Seqera Wave (arm64 image builds) and other network use.
    pub online: bool,
    /// --offline: the Nextflow container and every task container run with
    /// no network, so "offline" is enforced rather than requested.
    pub offline: bool,
    /// Directories the Nextflow container must see, beyond the above.
    pub extra_mounts: Vec<PathBuf>,
}

/// The x86 overlay the 30x golden VCF used (reference/cloud-hcjava.config).
pub const HCJAVA_CONFIG: &str = include_str!("../../../../reference/cloud-hcjava.config");
/// nextflow/nextflow is published for linux/amd64 only; on arm64 verify builds
/// this from the pinned image's launcher (reference/containers/nextflow).
pub const ARM64_NEXTFLOW_DOCKERFILE: &str = include_str!("../../../../reference/containers/nextflow/Dockerfile");

/// On x86 the pinned image runs as published except for its Docker CLI, which
/// is too old for Docker Engine 29 (reference/containers/nextflow/Dockerfile.amd64).
pub const AMD64_NEXTFLOW_DOCKERFILE: &str = include_str!("../../../../reference/containers/nextflow/Dockerfile.amd64");

/// The local Nextflow image for a descriptor's Nextflow version on this
/// architecture (`aarch64` or `x86_64`), and the Dockerfile that builds it.
pub fn local_nextflow_image(arch: &str, version: &str) -> (String, &'static str) {
    if arch == "aarch64" {
        (format!("nomics-verify/nextflow-arm64:{version}"), ARM64_NEXTFLOW_DOCKERFILE)
    } else {
        (format!("nomics-verify/nextflow-amd64:{version}"), AMD64_NEXTFLOW_DOCKERFILE)
    }
}

pub fn uid_gid() -> (u32, u32) {
    // Our own ids, from /proc (the crate takes no libc dependency).
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let field = |key: &str| {
        status
            .lines()
            .find_map(|l| l.strip_prefix(key))
            .and_then(|v| v.split_whitespace().next()?.parse().ok())
            .unwrap_or(0)
    };
    (field("Uid:"), field("Gid:"))
}

/// The group that owns the Docker socket, so the Nextflow container's user
/// may use it.
fn docker_socket_gid() -> Option<u32> {
    use std::os::unix::fs::MetadataExt;
    std::fs::metadata("/var/run/docker.sock").ok().map(|m| m.gid())
}

impl Plan<'_> {
    pub fn results(&self) -> PathBuf {
        self.dir.join("results")
    }

    pub fn work(&self) -> PathBuf {
        self.dir.join("work")
    }

    pub fn trace(&self) -> PathBuf {
        self.dir.join("trace.txt")
    }

    fn profile(&self) -> io::Result<&super::descriptor::Profile> {
        self.descriptor
            .profile(self.arch)
            .ok_or_else(|| invalid(format!("{} has no reference configuration for {}", self.descriptor.id, self.arch)))
    }

    /// The configuration name the report shows.
    pub fn config_name(&self) -> io::Result<String> {
        let p = self.profile()?;
        Ok(match p.overlay.as_str() {
            "hcjava" => format!("{} pinned: -profile {}, Java PairHMM and Smith-Waterman (cloud-hcjava.config)", self.arch, p.profiles),
            "arm64-pinned" => format!(
                "{} pinned: -profile {}, Wave off, {} task containers pinned by digest (the arm64 goldens' images)",
                self.arch,
                p.profiles,
                self.descriptor.pins(self.arch).len()
            ),
            other => format!("{} -profile {} + {other}", self.arch, p.profiles),
        })
    }

    pub fn samplesheet(&self) -> String {
        format!(
            "patient,sex,status,sample,lane,fastq_1,fastq_2\n{},NA,0,{},{},{},{}\n",
            self.patient,
            self.sample,
            self.lane,
            self.fastq1.display(),
            self.fastq2.display()
        )
    }

    /// sarek parameters, as `reference/params.tier1.yaml`.
    pub fn params(&self) -> String {
        let r = self.reference.display();
        let bundle = format!("{r}/Annotation/GATKBundle");
        let mut p = format!(
            "input: {dir}/samplesheet.csv\noutdir: {dir}/results\nigenomes_ignore: true\n\
             fasta: {r}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta\n\
             fasta_fai: {r}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai\n\
             dict: {r}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict\n\
             bwa: {r}/Sequence/BWAIndex/\n\
             dbsnp: {bundle}/dbsnp_146.hg38.vcf.gz\n\
             dbsnp_tbi: {bundle}/dbsnp_146.hg38.vcf.gz.tbi\n\
             known_indels: {bundle}/{{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}}.vcf.gz\n\
             known_indels_tbi: {bundle}/{{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}}.vcf.gz.tbi\n\
             wes: {wes}\nintervals: {intervals}\n\
             custom_config_base: {dir}/no-institutional-config\n",
            dir = self.dir.display(),
            wes = self.wes,
            intervals = self.intervals.display(),
        );
        for (k, v) in &self.descriptor.params {
            p.push_str(&format!("{k}: {v}\n"));
        }
        p
    }

    /// verify's Nextflow overlay: the goldens' per-task caps (fastp's thread
    /// count decides its chunks, so FASTP must get the goldens' 12 CPUs on
    /// any machine), the platform's pinned tool configuration, and a trace
    /// that names each task's container.
    pub fn overlay(&self) -> io::Result<String> {
        let (uid, gid) = uid_gid();
        // The goldens cap tasks at 20 CPUs (reference/spark.config); the
        // executor may oversubscribe a smaller machine rather than refuse.
        let executor_cpus = self.host_cpus.max(20);
        let memory = (self.host_memory_gib * 9 / 10).max(8);
        let mut c = format!(
            "// Written by nomics verify ({}).\n\
             process.resourceLimits = [cpus: 20, memory: {memory}.GB, time: 48.h]\n\
             executor {{ name = 'local'; cpus = {executor_cpus}; memory = {memory}.GB }}\n\
             docker.runOptions = '-u {uid}:{gid}{net}'\n\
             trace.fields = 'task_id,hash,native_id,name,status,exit,submit,duration,realtime,container'\n",
            self.descriptor.id,
            net = if self.offline { " --network none" } else { "" },
        );
        match self.profile()?.overlay.as_str() {
            "hcjava" => {
                c.push_str(HCJAVA_CONFIG);
                c.push('\n');
            }
            "arm64-pinned" => {
                // sarek's arm64 profile has Wave build images from each
                // module's conda spec; the GATK one (gcnvkernel, MKL-only)
                // cannot build for arm64. Run the goldens' images instead.
                c.push_str("wave.enabled = false\nprocess {\n");
                for pin in self.descriptor.pins(self.arch) {
                    c.push_str(&format!("    withName: '{}' {{ container = '{}' }}\n", pin.process, pin.image));
                }
                c.push_str("}\n");
            }
            other => return Err(invalid(format!("unknown overlay {other}"))),
        }
        Ok(c)
    }

    /// Write the samplesheet, params, overlay and an empty institutional
    /// config (sarek would otherwise fetch nf-core/configs).
    pub fn write_inputs(&self) -> io::Result<()> {
        std::fs::create_dir_all(self.dir.join("no-institutional-config/pipeline"))?;
        std::fs::write(self.dir.join("no-institutional-config/nfcore_custom.config"), "")?;
        std::fs::write(self.dir.join("no-institutional-config/pipeline/sarek.config"), "")?;
        std::fs::write(self.dir.join("samplesheet.csv"), self.samplesheet())?;
        std::fs::write(self.dir.join("params.yaml"), self.params())?;
        std::fs::write(self.dir.join("verify.config"), self.overlay()?)?;
        std::fs::create_dir_all(self.dir.join("home"))?;
        Ok(())
    }

    fn mounts(&self) -> Vec<PathBuf> {
        let mut m = vec![
            self.dir.clone(),
            self.sarek.clone(),
            self.nxf_home.clone(),
            self.reference.clone(),
        ];
        for p in [&self.fastq1, &self.fastq2, &self.intervals] {
            if let Some(d) = p.parent() {
                m.push(d.to_path_buf());
            }
        }
        m.extend(self.extra_mounts.iter().cloned());
        m.sort();
        m.dedup();
        // A path inside another mount needs no mount of its own.
        let all = m.clone();
        m.retain(|p| !all.iter().any(|q| q != p && p.starts_with(q)));
        m
    }

    /// `docker run … nextflow run <sarek> …`.
    pub fn command(&self) -> io::Result<Command> {
        let profile = self.profile()?;
        let mut c = nextflow_container(self.descriptor, &self.nxf_home, &self.dir.join("home"), &self.mounts(), self.online);
        if self.offline {
            // The Docker socket is a unix socket: task launches still work.
            c.args(["--network", "none"]);
        }
        c.current_dir(&self.dir);
        c.args(["-w", &self.dir.display().to_string()]);
        c.arg(&self.descriptor.nextflow_image);
        c.args(["nextflow", "-log"]).arg(self.dir.join("nextflow.log"));
        c.arg("run").arg(&self.sarek);
        c.args(["-profile", &profile.profiles]);
        c.arg("-params-file").arg(self.dir.join("params.yaml"));
        c.arg("-c").arg(self.dir.join("verify.config"));
        c.arg("-work-dir").arg(self.work());
        c.arg("-with-trace").arg(self.trace());
        if !self.online {
            c.arg("-offline");
        }
        Ok(c)
    }
}

/// The `docker run` prefix for the pinned Nextflow image (image and command
/// follow).
pub fn nextflow_container(d: &Descriptor, nxf_home: &Path, home: &Path, mounts: &[PathBuf], online: bool) -> Command {
    let (uid, gid) = uid_gid();
    let mut c = Command::new("docker");
    c.args(["run", "--rm", "-i"]);
    c.args(["-u", &format!("{uid}:{gid}")]);
    if let Some(g) = docker_socket_gid() {
        c.args(["--group-add", &g.to_string()]);
    }
    c.args(["-v", "/var/run/docker.sock:/var/run/docker.sock"]);
    for m in mounts {
        let s = m.display().to_string();
        c.args(["-v", &format!("{s}:{s}")]);
    }
    c.args(["-e", &format!("HOME={}", home.display())]);
    c.args(["-e", &format!("NXF_HOME={}", nxf_home.display())]);
    c.args(["-e", &format!("NXF_VER={}", d.nextflow_version)]);
    c.args(["-e", "NXF_DISABLE_CHECK_LATEST=true"]);
    c.args(["-e", "NXF_ANSI_LOG=false"]);
    if online {
        // A proxy the host needs, the containers need too.
        for v in ["HTTP_PROXY", "HTTPS_PROXY", "NO_PROXY", "http_proxy", "https_proxy", "no_proxy"] {
            if std::env::var_os(v).is_some() {
                c.args(["-e", v]);
            }
        }
    } else {
        c.args(["-e", "NXF_OFFLINE=true"]);
    }
    c
}

/// The commit a git checkout is at, read from `.git` without git.
pub fn git_head(dir: &Path) -> Option<String> {
    let git = dir.join(".git");
    let head = std::fs::read_to_string(git.join("HEAD")).ok()?;
    let head = head.trim();
    let Some(r) = head.strip_prefix("ref: ") else {
        return Some(head.to_owned());
    };
    if let Ok(s) = std::fs::read_to_string(git.join(r)) {
        return Some(s.trim().to_owned());
    }
    let packed = std::fs::read_to_string(git.join("packed-refs")).ok()?;
    packed
        .lines()
        .find_map(|l| l.strip_suffix(r).map(|h| h.trim().to_owned()))
}

/// One row of a Nextflow trace.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Task {
    pub hash: String,
    pub name: String,
    pub status: String,
    pub exit: String,
    pub submit: String,
    pub container: String,
}

pub fn parse_trace(text: &str) -> io::Result<Vec<Task>> {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().ok_or_else(|| invalid("empty trace"))?.split('\t').collect();
    let col = |name: &str| head.iter().position(|h| *h == name);
    let (hash, name, status) = match (col("hash"), col("name"), col("status")) {
        (Some(a), Some(b), Some(c)) => (a, b, c),
        _ => return Err(invalid("trace lacks hash, name or status columns")),
    };
    let (exit, submit, container) = (col("exit"), col("submit"), col("container"));
    let get = |f: &[&str], i: Option<usize>| i.and_then(|i| f.get(i)).map_or_else(String::new, |s| (*s).to_owned());
    Ok(lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Task {
                hash: get(&f, Some(hash)),
                name: get(&f, Some(name)),
                status: get(&f, Some(status)),
                exit: get(&f, exit),
                submit: get(&f, submit),
                container: get(&f, container),
            }
        })
        .collect())
}

/// Every task must have completed with exit 0. A failed task that Nextflow
/// retried successfully still left the run correct; one that never
/// succeeded did not.
pub fn check_tasks(tasks: &[Task]) -> io::Result<()> {
    if tasks.is_empty() {
        return Err(invalid("the reference run ran no tasks"));
    }
    let mut last: BTreeMap<&str, &Task> = BTreeMap::new();
    for t in tasks {
        last.insert(t.name.as_str(), t);
    }
    let bad: Vec<String> = last
        .values()
        .filter(|t| !matches!(t.status.as_str(), "COMPLETED" | "CACHED") || (t.exit != "0" && !t.exit.is_empty() && t.exit != "-"))
        .map(|t| format!("{} ({} exit {})", t.name, t.status, t.exit))
        .collect();
    if !bad.is_empty() {
        return Err(invalid(format!("reference tasks did not complete: {}", bad.join(", "))));
    }
    Ok(())
}

/// With pins, every task must have run in its process's pinned image: an
/// unpinned process would fall back to the module's own (amd64) image.
pub fn check_pins(tasks: &[Task], pins: &[&super::descriptor::ContainerPin]) -> io::Result<()> {
    if pins.is_empty() {
        return Ok(());
    }
    let mut bad = Vec::new();
    for (process, container) in containers(tasks) {
        match pins.iter().find(|p| p.process == process) {
            Some(p) if p.image == container => {}
            Some(p) => bad.push(format!("{process} ran {container}, pinned {}", p.image)),
            None => bad.push(format!("{process} ran {container}, which no pin covers")),
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(invalid(format!("reference tasks left their pinned containers: {}", bad.join("; "))))
    }
}

/// The containers the run used, by process name (the part of the task name
/// before ` (`), sorted.
pub fn containers(tasks: &[Task]) -> Vec<(String, String)> {
    let mut m: BTreeMap<String, String> = BTreeMap::new();
    for t in tasks {
        if t.container.is_empty() || t.container == "-" {
            continue;
        }
        let process = t.name.split(" (").next().unwrap_or(&t.name);
        let process = process.rsplit(':').next().unwrap_or(process);
        m.insert(process.to_owned(), t.container.clone());
    }
    m.into_iter().collect()
}

/// `docker image inspect` digest of an image, when Docker knows it.
pub fn image_digest(image: &str) -> Option<String> {
    let out = Command::new("docker")
        .args(["image", "inspect", "--format", "{{index .RepoDigests 0}}{{if not .RepoDigests}}{{.Id}}{{end}}", image])
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8_lossy(&out.stdout).trim().to_owned();
    let digest = s.rsplit('@').next().unwrap_or(&s).to_owned();
    (!digest.is_empty()).then_some(digest)
}

/// `@RG` `DS:` of a BAM's header: the FASTA path as bwa saw it, which
/// `nomics run --read-group-ds` must reproduce.
pub fn read_group_ds(bam: &Path) -> io::Result<Option<String>> {
    let header = read_bam_header(bam)?;
    let text = String::from_utf8_lossy(&header.text);
    Ok(text
        .lines()
        .filter(|l| l.starts_with("@RG\t"))
        .find_map(|l| l.split('\t').find_map(|f| f.strip_prefix("DS:").map(str::to_owned))))
}

fn read_bam_header(bam: &Path) -> io::Result<nomics_io::bam::Header> {
    let mut r = io::BufReader::new(nomics_io::gzip::GzReader::open(bam)?);
    nomics_io::bam::Header::read(&mut r)
}

/// The run's scatter as `nomics run --scatter-order` input, ported from
/// `scripts/scatter_order.py`: GatherBQSRReports' `--input` order, the
/// recal `samtools merge` input order, and merge's `time(NULL)` seed found
/// from the merged BAM's second `@RG` suffix. `None` when the run had one
/// interval group (no gather or merge task), where BED order is right.
pub fn scatter_order(tasks: &[Task], work: &Path, recal_bam: Option<&Path>, sample: &str) -> io::Result<Option<String>> {
    let command = |t: &Task| -> io::Result<String> {
        let (top, rest) = t.hash.split_once('/').ok_or_else(|| invalid(format!("bad task hash {}", t.hash)))?;
        let dir = std::fs::read_dir(work.join(top))?
            .filter_map(Result::ok)
            .find(|e| e.file_name().to_string_lossy().starts_with(rest))
            .ok_or_else(|| invalid(format!("no work directory for task {}", t.hash)))?;
        std::fs::read_to_string(dir.path().join(".command.sh"))
    };
    let gather = tasks.iter().find(|t| t.name.contains(":GATK4_GATHERBQSRREPORTS ("));
    let Some(gather) = gather else {
        return Ok(None);
    };
    let input = format!("--input {sample}_");
    let gather_order = groups_after(&command(gather)?, &input, ".recal.table");
    let mut out = String::from("# from the verify reference run\n");
    let merge = tasks
        .iter()
        .filter(|t| t.name.contains(":MERGE_BAM ("))
        .map(|t| command(t).map(|c| (t, c)))
        .collect::<io::Result<Vec<_>>>()?
        .into_iter()
        .find(|(_, c)| c.contains(".recal.bam"));
    let mut merge_order = Vec::new();
    if let Some((task, cmd)) = merge {
        merge_order = groups_after(&cmd, &format!("/{sample}_"), ".recal.bam");
        if let Some(bam) = recal_bam {
            let header = read_bam_header(bam)?;
            let text = String::from_utf8_lossy(&header.text);
            let suffix = text
                .lines()
                .filter_map(|l| l.strip_prefix("@RG\tID:"))
                .nth(1)
                .and_then(|id| id.split('\t').next()?.rsplit_once('-'))
                .and_then(|(_, hex)| u64::from_str_radix(hex, 16).ok())
                .ok_or_else(|| invalid("merged recal.bam has no suffixed second @RG"))?;
            let approx = submit_seconds(&task.submit).ok_or_else(|| invalid(format!("bad submit time {}", task.submit)))?;
            let seed = merge_seed(approx, suffix).ok_or_else(|| invalid("samtools merge seed not found"))?;
            out.push_str(&format!("merge-seed {seed}\n"));
        }
    }
    if merge_order.is_empty() {
        merge_order = gather_order.clone();
    }
    for g in &gather_order {
        out.push_str(&format!("gather {g}\n"));
    }
    for g in &merge_order {
        out.push_str(&format!("merge {g}\n"));
    }
    Ok(Some(out))
}

/// Group names in `text` between `start` and `end`, in order.
fn groups_after(text: &str, start: &str, end: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(i) = rest.find(start) {
        rest = &rest[i + start.len()..];
        let token_end = rest.find(char::is_whitespace).unwrap_or(rest.len());
        if let Some(name) = rest[..token_end].strip_suffix(end) {
            out.push(name.to_owned());
        }
    }
    out
}

/// `2026-10-10 04:41:58.022` as seconds, read as UTC (the host's zone is
/// unknown, so the seed search spans ±26 h).
fn submit_seconds(s: &str) -> Option<i64> {
    let (date, time) = s.split_once(' ')?;
    let mut d = date.split('-').map(|p| p.parse::<i64>());
    let (y, m, day) = (d.next()?.ok()?, d.next()?.ok()?, d.next()?.ok()?);
    let mut t = time.split('.').next()?.split(':').map(|p| p.parse::<i64>());
    let (hh, mm, ss) = (t.next()?.ok()?, t.next()?.ok()?, t.next()?.ok()?);
    // days_from_civil (Howard Hinnant).
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = (m + 9) % 12;
    let doy = (153 * mp + 2) / 5 + day - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    let days = era * 146_097 + doe - 719_468;
    Some(days * 86_400 + hh * 3600 + mm * 60 + ss)
}

/// The `hts_srand48(time(NULL))` seed whose first `lrand48()` draw is `suffix`.
fn merge_seed(approx: i64, suffix: u64) -> Option<i64> {
    const A: u64 = 0x5DEECE66D;
    const C: u64 = 0xB;
    const M: u64 = 1 << 48;
    (approx - 26 * 3600..approx + 26 * 3600).find(|&s| {
        let x = (((s as u64 & 0xFFFF_FFFF) << 16) | 0x330E) % M;
        let x = (A.wrapping_mul(x).wrapping_add(C)) % M;
        x >> 17 == suffix
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TRACE: &str = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\tcontainer\n\
        1\tab/cdef12\t1\tNFCORE_SAREK:SAREK:FASTP (S-L1)\tCOMPLETED\t0\t2026-10-10 04:41:58.022\t1m\t1m\tquay.io/biocontainers/fastp:1.1.0\n\
        2\t12/345678\t2\tNFCORE_SAREK:SAREK:BWAMEM1_MEM (S-L1)\tFAILED\t137\t2026-10-10 04:42:58.022\t1m\t1m\tquay.io/x/bwa:1\n\
        3\t13/345678\t3\tNFCORE_SAREK:SAREK:BWAMEM1_MEM (S-L1)\tCOMPLETED\t0\t2026-10-10 04:43:58.022\t1m\t1m\tquay.io/x/bwa:1\n";

    #[test]
    fn traces_and_task_checks() {
        let tasks = parse_trace(TRACE).unwrap();
        assert_eq!(tasks.len(), 3);
        // The retried BWA task finally completed.
        check_tasks(&tasks).unwrap();
        // Without the retry, the run did not complete.
        assert!(check_tasks(&tasks[..2]).unwrap_err().to_string().contains("BWAMEM1_MEM"));
        assert!(check_tasks(&[]).is_err());
        assert_eq!(
            containers(&tasks),
            [
                ("BWAMEM1_MEM".to_owned(), "quay.io/x/bwa:1".to_owned()),
                ("FASTP".to_owned(), "quay.io/biocontainers/fastp:1.1.0".to_owned())
            ]
        );
    }

    #[test]
    fn pinned_containers_are_checked() {
        use super::super::descriptor::ContainerPin;
        let tasks = parse_trace(TRACE).unwrap();
        let pin = |p: &str, i: &str| ContainerPin { arch: "aarch64".into(), process: p.into(), image: i.into() };
        let fastp = pin("FASTP", "quay.io/biocontainers/fastp:1.1.0");
        let bwa = pin("BWAMEM1_MEM", "quay.io/x/bwa:1");
        check_pins(&tasks, &[]).unwrap();
        check_pins(&tasks, &[&fastp, &bwa]).unwrap();
        // A process with no pin, and one that ran another image, both fail.
        assert!(check_pins(&tasks, &[&fastp]).unwrap_err().to_string().contains("BWAMEM1_MEM ran quay.io/x/bwa:1, which no pin covers"));
        let other = pin("BWAMEM1_MEM", "quay.io/x/bwa:2");
        assert!(check_pins(&tasks, &[&fastp, &other]).unwrap_err().to_string().contains("pinned quay.io/x/bwa:2"));
    }

    #[test]
    fn merge_seed_matches_the_recorded_golden() {
        // reference/wgs30-x86-r1.scatter-order.txt: merge-seed 1791676187. Its
        // first lrand48 draw is the suffix; search from an hour off.
        let seed: i64 = 1_791_676_187;
        let x = (((seed as u64 & 0xFFFF_FFFF) << 16) | 0x330E) % (1 << 48);
        let suffix = (0x5DEECE66D_u64.wrapping_mul(x).wrapping_add(0xB) % (1 << 48)) >> 17;
        assert_eq!(merge_seed(seed + 3600, suffix), Some(seed));
        assert_eq!(submit_seconds("1970-01-02 00:00:01.5"), Some(86_401));
        assert_eq!(submit_seconds("2026-10-30 00:00:00.000"), Some(1_793_318_400));
    }

    #[test]
    fn group_names_from_commands() {
        let cmd = "gatk GatherBQSRReports --input S_chr1_1-100.recal.table --input S_chr2_5-9.recal.table --output S.recal.table";
        assert_eq!(groups_after(cmd, "--input S_", ".recal.table"), ["chr1_1-100", "chr2_5-9"]);
        let merge = "samtools merge out.bam 1/S_chr2_5-9.recal.bam 2/S_chr1_1-100.recal.bam";
        assert_eq!(groups_after(merge, "/S_", ".recal.bam"), ["chr2_5-9", "chr1_1-100"]);
    }

    #[test]
    fn git_heads() {
        let d = std::env::temp_dir().join(format!("nomics-git-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join(".git/refs/heads")).unwrap();
        std::fs::write(d.join(".git/HEAD"), "8ccac7ad37b05dd792447763bf9671b719824587\n").unwrap();
        assert_eq!(git_head(&d).as_deref(), Some("8ccac7ad37b05dd792447763bf9671b719824587"));
        std::fs::write(d.join(".git/HEAD"), "ref: refs/heads/main\n").unwrap();
        std::fs::write(d.join(".git/packed-refs"), "# pack-refs\nabc123 refs/heads/main\n").unwrap();
        assert_eq!(git_head(&d).as_deref(), Some("abc123"));
        std::fs::remove_dir_all(&d).unwrap();
    }
}
