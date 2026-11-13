//! A pipeline descriptor: the reference pipeline, its pins, the values taken
//! from its run, and for each published artifact the comparator and where
//! each side writes it. Data, not code paths, so a later harness can add
//! pipelines without changing verify (`docs/VERIFY.md`).

use std::io;

/// The descriptor shipped in the binary.
pub const SAREK_GERMLINE: &str = include_str!("../../descriptors/sarek-3.10.0-germline.tsv");

/// How an artifact is compared; each names a digest and a detailed comparator
/// in `gate.rs`, and a normalization row in `docs/GATES.md`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comparator {
    FastqChunk,
    Bam,
    MdMetrics,
    Raw,
    Vcf,
}

impl Comparator {
    fn parse(s: &str) -> Option<Comparator> {
        Some(match s {
            "fastq_chunk" => Comparator::FastqChunk,
            "bam" => Comparator::Bam,
            "md_metrics" => Comparator::MdMetrics,
            "raw" => Comparator::Raw,
            "vcf" => Comparator::Vcf,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            Comparator::FastqChunk => "fastq_chunk",
            Comparator::Bam => "bam",
            Comparator::MdMetrics => "md_metrics",
            Comparator::Raw => "raw",
            Comparator::Vcf => "vcf",
        }
    }
}

/// Where one side writes an artifact: roots tried in order, then file name
/// patterns (`*` and `?`) tried in order within each root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Location {
    pub roots: Vec<String>,
    pub patterns: Vec<String>,
}

impl Location {
    fn parse(s: &str) -> Option<Location> {
        let (roots, patterns) = s.split_once(':')?;
        let split = |s: &str| -> Vec<String> {
            s.split(',')
                .map(str::trim)
                .filter(|p| !p.is_empty())
                .map(str::to_owned)
                .collect()
        };
        let location = Location {
            roots: split(roots),
            patterns: split(patterns),
        };
        let known = |r: &String| matches!(r.as_str(), "work" | "results" | "run");
        (!location.roots.is_empty()
            && !location.patterns.is_empty()
            && location.roots.iter().all(known))
        .then_some(location)
    }

    /// The patterns with `{prefix}` and `{lane}` filled in.
    pub fn patterns_for(&self, prefix: &str, lane: &str) -> Vec<String> {
        self.patterns
            .iter()
            .map(|p| p.replace("{prefix}", prefix).replace("{lane}", lane))
            .collect()
    }
}

#[derive(Debug, Clone)]
pub struct ArtifactSpec {
    /// `md_bam`, or `chunk/*` for one artifact per matching file.
    pub name: String,
    pub comparator: Comparator,
    pub reference: Location,
    pub nomics: Location,
    pub normalization: String,
}

impl ArtifactSpec {
    /// The `chunk` of `chunk/*`, for per-file artifacts.
    pub fn group(&self) -> Option<&str> {
        self.name.strip_suffix("/*")
    }
}

#[derive(Debug, Clone)]
pub struct Profile {
    pub arch: String,
    pub profiles: String,
    pub overlay: String,
}

#[derive(Debug, Clone)]
pub struct ContainerPin {
    pub arch: String,
    /// Process simple name, e.g. `GATK4_MARKDUPLICATES`.
    pub process: String,
    /// `repository:tag@sha256:<digest>`.
    pub image: String,
}

#[derive(Debug, Clone)]
pub struct Descriptor {
    pub id: String,
    pub pipeline: String,
    pub version: String,
    pub commit: String,
    pub source_url: String,
    pub source_tag: String,
    pub nextflow_version: String,
    pub nextflow_image: String,
    pub profiles: Vec<Profile>,
    pub params: Vec<(String, String)>,
    pub captures: Vec<(String, String)>,
    pub artifacts: Vec<ArtifactSpec>,
    pub refusals: Vec<(String, String)>,
    /// Task containers pinned per architecture, by process name, as
    /// `image@sha256:<digest>`.
    pub containers: Vec<ContainerPin>,
    /// SHA-256 of the descriptor text, recorded in the report.
    pub sha256: String,
}

fn invalid(line: usize, msg: impl std::fmt::Display) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("descriptor line {line}: {msg}"),
    )
}

impl Descriptor {
    pub fn parse(text: &str) -> io::Result<Descriptor> {
        let mut fields: std::collections::BTreeMap<&str, Vec<&str>> = Default::default();
        let mut d = Descriptor {
            id: String::new(),
            pipeline: String::new(),
            version: String::new(),
            commit: String::new(),
            source_url: String::new(),
            source_tag: String::new(),
            nextflow_version: String::new(),
            nextflow_image: String::new(),
            profiles: Vec::new(),
            params: Vec::new(),
            captures: Vec::new(),
            artifacts: Vec::new(),
            refusals: Vec::new(),
            containers: Vec::new(),
            sha256: super::sha256_hex(text.as_bytes()),
        };
        for (i, line) in text.lines().enumerate() {
            let n = i + 1;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            let (key, rest) = (cols[0], &cols[1..]);
            let want = |k: usize| -> io::Result<()> {
                if rest.len() == k {
                    Ok(())
                } else {
                    Err(invalid(n, format!("{key} takes {k} values, found {}", rest.len())))
                }
            };
            match key {
                "format" | "id" | "pipeline" | "version" | "commit" => {
                    want(1)?;
                    if fields.insert(key, rest.to_vec()).is_some() {
                        return Err(invalid(n, format!("duplicate {key}")));
                    }
                }
                "source" | "nextflow" => {
                    want(2)?;
                    if fields.insert(key, rest.to_vec()).is_some() {
                        return Err(invalid(n, format!("duplicate {key}")));
                    }
                }
                "profile" => {
                    want(3)?;
                    d.profiles.push(Profile {
                        arch: rest[0].to_owned(),
                        profiles: rest[1].to_owned(),
                        overlay: rest[2].to_owned(),
                    });
                }
                "param" => {
                    want(2)?;
                    d.params.push((rest[0].to_owned(), rest[1].to_owned()));
                }
                "capture" => {
                    want(2)?;
                    d.captures.push((rest[0].to_owned(), rest[1].to_owned()));
                }
                "container" => {
                    want(3)?;
                    let image = rest[2];
                    let digest = image.rsplit_once("@sha256:").map(|(_, h)| h).unwrap_or("");
                    if digest.len() != 64 || !digest.bytes().all(|b| b.is_ascii_hexdigit()) {
                        return Err(invalid(n, format!("container {image} is not pinned by @sha256 digest")));
                    }
                    if d.containers.iter().any(|c| c.arch == rest[0] && c.process == rest[1]) {
                        return Err(invalid(n, format!("duplicate container for {} on {}", rest[1], rest[0])));
                    }
                    d.containers.push(ContainerPin {
                        arch: rest[0].to_owned(),
                        process: rest[1].to_owned(),
                        image: image.to_owned(),
                    });
                }
                "refuse" => {
                    want(2)?;
                    d.refusals.push((rest[0].to_owned(), rest[1].to_owned()));
                }
                "artifact" => {
                    want(5)?;
                    let comparator = Comparator::parse(rest[1])
                        .ok_or_else(|| invalid(n, format!("unknown comparator {}", rest[1])))?;
                    let location = |s: &str| {
                        Location::parse(s).ok_or_else(|| invalid(n, format!("bad location {s}")))
                    };
                    let spec = ArtifactSpec {
                        name: rest[0].to_owned(),
                        comparator,
                        reference: location(rest[2])?,
                        nomics: location(rest[3])?,
                        normalization: rest[4].to_owned(),
                    };
                    if d.artifacts.iter().any(|a| a.name == spec.name) {
                        return Err(invalid(n, format!("duplicate artifact {}", spec.name)));
                    }
                    d.artifacts.push(spec);
                }
                other => return Err(invalid(n, format!("unknown key {other}"))),
            }
        }
        let one = |k: &str| -> io::Result<String> {
            fields
                .get(k)
                .map(|v| v[0].to_owned())
                .ok_or_else(|| invalid(0, format!("missing {k}")))
        };
        if one("format")? != "nomics-descriptor/1" {
            return Err(invalid(0, "unsupported format"));
        }
        d.id = one("id")?;
        d.pipeline = one("pipeline")?;
        d.version = one("version")?;
        d.commit = one("commit")?;
        let source = fields.get("source").ok_or_else(|| invalid(0, "missing source"))?;
        d.source_url = source[0].to_owned();
        d.source_tag = source[1].to_owned();
        let nextflow = fields
            .get("nextflow")
            .ok_or_else(|| invalid(0, "missing nextflow"))?;
        d.nextflow_version = nextflow[0].to_owned();
        d.nextflow_image = nextflow[1].to_owned();
        if d.commit.len() != 40 || !d.commit.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(invalid(0, "commit must be a full 40-character hash"));
        }
        if d.artifacts.is_empty() {
            return Err(invalid(0, "no artifacts"));
        }
        Ok(d)
    }

    pub fn profile(&self, arch: &str) -> Option<&Profile> {
        self.profiles.iter().find(|p| p.arch == arch)
    }

    /// The pinned task containers for an architecture.
    pub fn pins(&self, arch: &str) -> Vec<&ContainerPin> {
        self.containers.iter().filter(|c| c.arch == arch).collect()
    }

    pub fn refusal(&self, key: &str) -> String {
        self.refusals
            .iter()
            .find(|(k, _)| k == key)
            .map_or_else(|| key.to_owned(), |(_, why)| why.clone())
    }
}

/// `*` matches any run of characters and `?` exactly one; nothing else is special.
pub fn glob_match(pattern: &str, name: &str) -> bool {
    let (p, n) = (pattern.as_bytes(), name.as_bytes());
    let (mut pi, mut ni) = (0, 0);
    let (mut star, mut mark) = (None, 0);
    while ni < n.len() {
        if pi < p.len() && (p[pi] == b'?' || p[pi] == n[ni]) {
            pi += 1;
            ni += 1;
        } else if pi < p.len() && p[pi] == b'*' {
            star = Some(pi);
            mark = ni;
            pi += 1;
        } else if let Some(s) = star {
            pi = s + 1;
            mark += 1;
            ni = mark;
        } else {
            return false;
        }
    }
    while pi < p.len() && p[pi] == b'*' {
        pi += 1;
    }
    pi == p.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shipped_descriptor_parses() {
        let d = Descriptor::parse(SAREK_GERMLINE).unwrap();
        assert_eq!(d.id, "sarek-3.10.0-germline");
        assert_eq!(d.commit, "8ccac7ad37b05dd792447763bf9671b719824587");
        let names: Vec<&str> = d.artifacts.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(
            names,
            [
                "chunk/*",
                "mapped/*",
                "md_bam",
                "md_metrics",
                "recal_table",
                "recal_bam",
                "vcf"
            ]
        );
        assert!(d.profile("x86_64").is_some() && d.profile("aarch64").is_some());
        assert_eq!(d.sha256.len(), 64);
        // Every arm64 process the goldens ran is pinned; x86 pins by tag in sarek.
        assert_eq!(d.pins("aarch64").len(), 19);
        assert!(d.pins("x86_64").is_empty());
        assert!(d.pins("aarch64").iter().any(|p| p.process == "GATK4_HAPLOTYPECALLER"
            && p.image.contains("gatk4_htslib_samtools:6fd46722be032bdd@sha256:43aa6d13")));
    }

    #[test]
    fn rejects_bad_descriptors() {
        let base = "format\tnomics-descriptor/1\nid\tx\npipeline\tp\nversion\t1\n\
                    commit\t8ccac7ad37b05dd792447763bf9671b719824587\nsource\tu\tt\nnextflow\t1\timg\n";
        let ok = format!("{base}artifact\tvcf\tvcf\tresults:a.vcf.gz\trun:a.vcf.gz\tn\n");
        assert!(Descriptor::parse(&ok).is_ok());
        // No artifacts, an unknown comparator, an unknown root, a duplicate.
        assert!(Descriptor::parse(base).is_err());
        assert!(Descriptor::parse(&ok.replace("\tvcf\tresults", "\tcram\tresults")).is_err());
        assert!(Descriptor::parse(&ok.replace("results:", "elsewhere:")).is_err());
        let dup = format!("{ok}artifact\tvcf\tvcf\tresults:b\trun:b\tn\n");
        assert!(Descriptor::parse(&dup).is_err());
        assert!(Descriptor::parse(&ok.replace("8ccac7ad37b0", "8ccac7ad")).is_err());
        // Containers must be pinned by digest, once per process and architecture.
        let pin = format!("container\taarch64\tFASTP\tr/fastp:1@sha256:{}\n", "a".repeat(64));
        assert_eq!(Descriptor::parse(&format!("{ok}{pin}")).unwrap().pins("aarch64").len(), 1);
        assert!(Descriptor::parse(&format!("{ok}container\taarch64\tFASTP\tr/fastp:1\n")).is_err());
        assert!(Descriptor::parse(&format!("{ok}{}", pin.replace(&"a".repeat(64), "abc"))).is_err());
        assert!(Descriptor::parse(&format!("{ok}{pin}{pin}")).is_err());
    }

    #[test]
    fn globs() {
        assert!(glob_match("S.????.bam", "S.0001.bam"));
        assert!(!glob_match("S.????.bam", "S.md.bam"));
        assert!(glob_match("*.fastp.fastq.gz", "0001.S-L1_R1.fastp.fastq.gz"));
        assert!(!glob_match("*.fastp.fastq.gz", "0001.S-L1_R1.fastp.json"));
        assert!(glob_match("*", ""));
        assert!(glob_match("a*b*c", "aXXbYc"));
        assert!(!glob_match("a*b*c", "aXXbY"));
    }
}
