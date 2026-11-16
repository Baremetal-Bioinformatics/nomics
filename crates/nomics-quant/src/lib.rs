//! Fresh alignment-mode quantifier for the nf-core/rnaseq 3.27.0 Salmon 2.7.0
//! command. The spec is upstream tag v2.7.0, commit d3be8b31d31242683c01fc5510c48fbd8bc917fb.
//! See `docs/research/pipelines/nomics-quant.md` for the source ledger and gate.
//!
//! Supported scope: name-grouped transcriptome BAM, all Salmon library types, deterministic
//! AS-score weights (missing AS = 0), VBEM, default fragment model, no bias flags.
//! Unsupported options are rejected by the CLI rather than silently ignored.
use std::io;
use std::path::PathBuf;

mod fld;
mod fragment;
mod inference;
mod input;
mod library;
pub use library::LibraryType;
mod output;

use fragment::{Hit, Status};
use nomics_mem::HostArray;

/// The nf-core alignment quantification inputs. `threads` controls inference
/// scheduling, never reduction boundaries. BGZF inflate, the BAM scan workers,
/// the placement-weight pass and the gene-map parse use all available cores;
/// their integer/keyed reductions are exact in any grouping.
#[derive(Clone, Debug)]
pub struct Options {
    pub transcripts: PathBuf,
    pub alignments: PathBuf,
    pub gene_map: Option<PathBuf>,
    pub output: PathBuf,
    pub lib_type: String,
    /// Incompatible placements are dropped when <= 0; otherwise weighted by this prior.
    pub incompat_prior: f64,
    pub threads: usize,
}

/// Unrounded transcript estimates, in BAM `@SQ` order.
#[derive(Debug)]
pub struct Quantification {
    pub names: Vec<String>,
    pub lengths: HostArray<u32>,
    pub effective_lengths: HostArray<f64>,
    pub counts: HostArray<f64>,
    pub tpm: HostArray<f64>,
    pub num_processed: u64,
    pub num_mapped: u64,
    pub num_compatible: u64,
    pub num_incompatible: u64,
    pub format_counts: [u64; 12],
    /// Canonical format used by filtering; A (or a) if detection was unresolved.
    pub expected_format: String,
    pub num_eq_classes: usize,
    pub iterations: u32,
    pub converged: bool,
    pub unique: HostArray<u32>,
    pub ambiguous: HostArray<u32>,
    /// Normalized fragment probabilities indexed by lengths 0..=1000.
    pub fragment_pmf: Vec<f64>,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Quantify without writing files. One BAM pass tallies the fragment model and
/// retains compact placements in HostBuf-backed chunks; no BAM records, FASTA
/// sequences, or intermediate RAD need to stay in memory.
pub fn quantify(options: &Options) -> io::Result<Quantification> {
    let requested: LibraryType = options.lib_type.parse()?;
    if options.threads == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "threads must be positive",
        ));
    }
    // NOMICS_QUANT_INFERENCE_TIMING also prints input and output phase times.
    let mut clock =
        std::env::var_os("NOMICS_QUANT_INFERENCE_TIMING").map(|_| std::time::Instant::now());
    let phase = |label: &str, start: &mut Option<std::time::Instant>| {
        if let Some(t) = start {
            eprintln!(
                "salmon-quant timing: {label} {:.3}s",
                t.elapsed().as_secs_f64()
            );
            *t = std::time::Instant::now();
        }
    };
    let (header, fragments, tally) = fragment::scan(&options.alignments)?;
    let fragment::Tally {
        num_processed,
        opposite,
        same,
        detection: detection_counts,
    } = tally;
    phase("input scan (decode+placements+store)", &mut clock);
    input::validate_fasta(&options.transcripts, &header)?;
    phase("validate fasta", &mut clock);
    let distribution =
        fld::Distribution::new(if opposite.iter().sum::<u64>() >= same.iter().sum() {
            &opposite
        } else {
            &same
        });
    let lengths = HostArray::collect_exact(header.refs.iter().map(|r| r.len))?;
    let effective_lengths = HostArray::collect_exact(
        lengths
            .iter()
            .map(|&l| distribution.effective_length(l as usize)),
    )?;
    // Keep the upstream raw-string mode decisions (including lowercase A's
    // unresolved behavior), separately from its case-insensitive format parser.
    let paired = !matches!(options.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    let detected = library::detect(&detection_counts);
    let expected = if options.lib_type == "A" {
        detected
    } else if requested == LibraryType::Auto {
        None
    } else {
        Some(requested)
    };
    let expected_format = if requested == LibraryType::Auto {
        expected
            .map(|f| f.as_str().to_owned())
            .unwrap_or_else(|| options.lib_type.clone())
    } else {
        requested.as_str().to_owned()
    };
    // Each fragment's weights are a pure function of the frozen model; every
    // accumulator is an integer (or a key-ordered integer class map), so
    // fragment ranges are processed in parallel and merged exactly.
    let second = Second {
        expected,
        incompat_prior: options.incompat_prior,
        lengths: &lengths,
        distribution: &distribution,
        paired,
    };
    let Visited {
        builder,
        format_counts,
        num_mapped,
        num_compatible,
        num_incompatible,
    } = {
        use rayon::prelude::*;
        fragments
            .ranges()
            .into_par_iter()
            .fold(Visited::default, |mut acc, range| {
                fragments.visit_range(range, |hits| second.visit(hits, &mut acc));
                acc
            })
            .reduce(Visited::default, Visited::merge)
    };
    phase("placement weights + eq-class build", &mut clock);
    drop(fragments);
    let classes = builder.finish(&effective_lengths);
    phase("eq-class finish", &mut clock);
    let n = lengths.len();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()
        .map_err(io::Error::other)?;
    let (counts, iterations, converged) = pool.install(|| inference::optimize(&classes, n))?;
    phase("inference", &mut clock);
    let mut unique = HostArray::<u32>::zeroed(n)?;
    let mut ambiguous = HostArray::<u32>::zeroed(n)?;
    for c in &classes {
        let bins = if c.tids.len() == 1 {
            &mut unique
        } else {
            &mut ambiguous
        };
        for &t in &c.tids {
            bins[t] = bins[t].wrapping_add(c.count as u32);
        }
    }
    let rates = HostArray::collect_exact(
        counts
            .iter()
            .zip(effective_lengths.iter())
            .map(|(c, l)| if *l > 0.0 { c / l } else { 0.0 }),
    )?;
    let total: f64 = rates.iter().sum();
    let tpm = HostArray::collect_exact(
        rates
            .iter()
            .map(|r| if total > 0.0 { r / total * 1e6 } else { 0.0 }),
    )?;
    Ok(Quantification {
        names: header
            .refs
            .iter()
            .map(|r| String::from_utf8_lossy(&r.name).into_owned())
            .collect(),
        lengths,
        effective_lengths,
        counts,
        tpm,
        num_processed,
        num_mapped,
        num_compatible,
        num_incompatible,
        format_counts,
        expected_format,
        num_eq_classes: classes.len(),
        iterations,
        converged,
        unique,
        ambiguous,
        fragment_pmf: distribution.pmf.iter().map(|p| p.exp()).collect(),
    })
}

/// The frozen inputs of the second (placement-weight) pass.
struct Second<'a> {
    expected: Option<LibraryType>,
    incompat_prior: f64,
    lengths: &'a [u32],
    distribution: &'a fld::Distribution,
    paired: bool,
}

/// Second-pass accumulators: all integer, merged exactly in any grouping.
#[derive(Default)]
struct Visited {
    builder: inference::Builder,
    format_counts: [u64; 12],
    num_mapped: u64,
    num_compatible: u64,
    num_incompatible: u64,
}

impl Visited {
    fn merge(mut self, other: Visited) -> Visited {
        self.builder = self.builder.merge(other.builder);
        for (a, b) in self.format_counts.iter_mut().zip(other.format_counts) {
            *a += b;
        }
        self.num_mapped += other.num_mapped;
        self.num_compatible += other.num_compatible;
        self.num_incompatible += other.num_incompatible;
        self
    }
}

impl Second<'_> {
    fn visit(&self, hits: &[Hit], acc: &mut Visited) {
        let expected = self.expected;
        if hits.is_empty() {
            return;
        }
        let mut mask = 0u16;
        let mut compatible = false;
        for hit in hits {
            mask |= 1 << hit.format();
            compatible |= expected.is_some_and(|f| f.compatible(hit));
        }
        for (i, count) in acc.format_counts.iter_mut().enumerate() {
            *count += u64::from(mask & (1 << i) != 0);
        }
        if expected.is_some() {
            if compatible {
                acc.num_compatible += 1;
            } else {
                acc.num_incompatible += 1;
            }
        }
        let best = hits.iter().map(|h| h.score).max().unwrap_or(0);
        let mut weights: Vec<(usize, f64)> = hits
            .iter()
            .filter_map(|hit| {
                let incompatible = expected.is_some_and(|f| !f.compatible(hit));
                if incompatible && self.incompat_prior <= 0.0 {
                    return None;
                }
                let mut weight =
                    placement_log_weight(hit, best, self.lengths, self.distribution, self.paired);
                if incompatible {
                    weight += self.incompat_prior.ln();
                }
                Some((hit.tid, weight))
            })
            .collect();
        if weights.is_empty() {
            return;
        }
        acc.num_mapped += 1;
        weights.sort_by_key(|w| w.0);
        let mut tids = Vec::new();
        let mut logs = Vec::new();
        for run in weights.chunk_by(|a, b| a.0 == b.0) {
            tids.push(run[0].0);
            let max = run.iter().map(|w| w.1).fold(f64::NEG_INFINITY, f64::max);
            logs.push(if max == f64::NEG_INFINITY {
                max
            } else {
                max + run.iter().map(|w| (w.1 - max).exp()).sum::<f64>().ln()
            });
        }
        let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut weights: Vec<f64> = logs.iter().map(|l| (l - max).exp()).collect();
        let total: f64 = weights.iter().sum();
        if total > 0.0 {
            for w in &mut weights {
                *w /= total;
            }
        }
        acc.builder.add(tids, weights);
    }
}

fn placement_log_weight(
    h: &Hit,
    best: i32,
    lengths: &[u32],
    d: &fld::Distribution,
    paired: bool,
) -> f64 {
    let len = lengths[h.tid] as usize;
    let score = -(best - h.score) as f64;
    let reference_mass = d.cmf[len.clamp(1, fld::MAX)];
    let fragment_mass = if h.status == Status::Pair && h.length > 0 && h.length != u16::MAX as usize
    {
        d.pmf[h.length.min(len.max(1)).min(fld::MAX)] - reference_mass
    } else if !paired {
        0.0
    } else {
        let read_length = if h.length == u16::MAX as usize {
            0
        } else {
            h.length
        };
        let bound = if h.forward {
            len.saturating_sub(h.pos)
        } else {
            (h.pos + read_length).min(len)
        };
        if reference_mass == f64::NEG_INFINITY {
            -1e10
        } else {
            d.cmf[bound.min(fld::MAX)] - reference_mass
        }
    };
    score + fragment_mass
}

/// Quantify and write quant.sf, optional quant.genes.sf, and auxiliary files.
pub fn run(options: &Options) -> io::Result<Quantification> {
    let quantification = quantify(options)?;
    let start =
        std::env::var_os("NOMICS_QUANT_INFERENCE_TIMING").map(|_| std::time::Instant::now());
    output::write(options, &quantification)?;
    if let Some(t) = start {
        eprintln!(
            "salmon-quant timing: output {:.3}s",
            t.elapsed().as_secs_f64()
        );
    }
    Ok(quantification)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn as_penalty_and_orphan_bounds_affect_likelihoods() {
        let d = fld::Distribution::new(&vec![0; fld::MAX + 1]);
        let mut h = Hit {
            tid: 0,
            pos: 100,
            forward: false,
            mate_forward: true,
            length: 250,
            score: 100,
            status: Status::Pair,
        };
        let best = placement_log_weight(&h, 100, &[2000], &d, true);
        h.score = 97;
        let worse = placement_log_weight(&h, 100, &[2000], &d, true);
        assert_eq!(worse, best - 3.0);
        h.status = Status::Left;
        h.length = 75;
        h.pos = 1000;
        let ample = placement_log_weight(&h, 100, &[2000], &d, true);
        h.pos = 0;
        let bounded = placement_log_weight(&h, 100, &[2000], &d, true);
        assert!(bounded < ample - 1.0);
        // A saturated proper-pair length becomes RAD's unpaired sentinel;
        // the ambiguous length bound then has no read-length information.
        h.status = Status::Pair;
        h.length = u16::MAX as usize;
        assert!(placement_log_weight(&h, 100, &[2000], &d, true) < ample - 1.0);
        // Single-end libraries do not use an orphan-length penalty.
        h.status = Status::Single;
        assert_eq!(placement_log_weight(&h, 100, &[2000], &d, false), -3.0);
    }
}
