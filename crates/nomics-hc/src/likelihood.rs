//! Read × haplotype likelihoods (spec §6): `PairHMMLikelihoodCalculationEngine.
//! computeReadLikelihoods` with HaplotypeCaller's defaults — read
//! preparation, the Java PairHMM, symmetric normalization to within 4.5 of
//! each read's best haplotype, and poorly-modeled-read disqualification.
//! Gate: `tests/likelihood.rs` against the reference jar (GatkLikelihoods.java).

use crate::pairhmm_trie::{self, HmmRead};
use crate::readprep::hmm_qualities;

/// `--pair-hmm-gap-continuation-penalty`.
const GCP: u8 = 10;
/// `log10GlobalReadMismappingRate = qualToErrorProbLog10(45)`.
const MISMAPPING_CAP: f64 = -4.5;
/// `--expected-mismatch-rate-for-read-disqualification`.
const EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;

/// A read as the likelihood engine sees it.
pub struct LkRead<'a> {
    pub bases: &'a [u8],
    pub quals: &'a [u8],
    pub mapq: u8,
}

/// Result: `values[h][r]` for the kept reads (input order), and the input
/// indices of kept and disqualified reads.
pub struct Likelihoods {
    pub values: Vec<Vec<f64>>,
    pub kept: Vec<usize>,
    pub disqualified: Vec<usize>,
}

/// `log10MinTrueLikelihood(0.02, capLikelihoods = true)` for a read of `len` bases.
fn min_true_likelihood(len: usize) -> f64 {
    let max_errors = 2.0f64.min((len as f64 * EXPECTED_ERROR_RATE_PER_BASE).ceil());
    max_errors * -4.0
}

pub fn compute(haplotypes: &[&[u8]], reads: &[LkRead]) -> Likelihoods {
    if let Some(dump) = dump_sink() {
        dump_batch(dump, haplotypes, reads);
    }
    let prepared: Vec<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)> = reads
        .iter()
        .map(|r| {
            let (q, ins, del) = hmm_qualities(r.bases, r.quals, r.mapq);
            (q, ins, del, vec![GCP; r.bases.len()])
        })
        .collect();
    let hmm_reads: Vec<HmmRead> = reads
        .iter()
        .zip(&prepared)
        .map(|(r, (q, ins, del, gcp))| HmmRead {
            bases: r.bases,
            quals: q,
            ins,
            del,
            gcp,
        })
        .collect();
    // NOMICS_PAIRHMM_LANES=1: the per-read haplotype lanes instead of the
    // prefix-sharing traversal (a diagnostic; both are bit-identical).
    let values = if std::env::var_os("NOMICS_PAIRHMM_LANES").is_some() {
        let mut hmm = crate::pairhmm::PairHmmLanes::default();
        let lanes = crate::pairhmm::LANES;
        let mut values = vec![Vec::with_capacity(reads.len()); haplotypes.len()];
        for r in &hmm_reads {
            let mut out = vec![0.0; lanes];
            for (k, haps) in haplotypes.chunks(lanes).enumerate() {
                let out = &mut out[..haps.len()];
                hmm.log10_likelihoods(haps, r.bases, r.quals, r.ins, r.del, r.gcp, out);
                for (h, &v) in out.iter().enumerate() {
                    values[k * lanes + h].push(v);
                }
            }
        }
        values
    } else {
        pairhmm_trie::likelihoods(haplotypes, &hmm_reads)
    };
    let lens: Vec<usize> = reads.iter().map(|r| r.bases.len()).collect();
    finish(values, &lens)
}

/// Normalization and disqualification of raw PairHMM log10 likelihoods
/// (`values[h][r]`, reads of lengths `read_lens`), shared by the CPU and GPU
/// paths.
pub fn finish(mut values: Vec<Vec<f64>>, read_lens: &[usize]) -> Likelihoods {
    let nh = values.len();
    // normalizeLikelihoods(-4.5, symmetric = true): cap below the first-max allele.
    if nh > 1 {
        for ri in 0..read_lens.len() {
            let mut best = values[0][ri];
            for v in values.iter().skip(1) {
                if v[ri] > best {
                    best = v[ri];
                }
            }
            let cap = best + MISMAPPING_CAP;
            for v in values.iter_mut() {
                if v[ri] < cap {
                    v[ri] = cap;
                }
            }
        }
    }
    // filterPoorlyModeledEvidence
    let (mut kept, mut disqualified) = (Vec::new(), Vec::new());
    for (ri, &len) in read_lens.iter().enumerate() {
        let mut max = f64::NEG_INFINITY;
        for v in &values {
            if v[ri] > max {
                max = v[ri];
            }
        }
        if max < min_true_likelihood(len) {
            disqualified.push(ri);
        } else {
            kept.push(ri);
        }
    }
    let values = values
        .into_iter()
        .map(|v| kept.iter().map(|&i| v[i]).collect())
        .collect();
    Likelihoods {
        values,
        kept,
        disqualified,
    }
}

/// Concatenate `compute` results of consecutive read chunks (`chunk` reads
/// each, the last possibly shorter). Every per-read step is independent, so
/// this is bit-identical to one `compute` over all reads.
pub fn merge(parts: Vec<Likelihoods>, chunk: usize, haplotypes: usize) -> Likelihoods {
    let mut out = Likelihoods {
        values: vec![Vec::new(); haplotypes],
        kept: Vec::new(),
        disqualified: Vec::new(),
    };
    for (k, p) in parts.into_iter().enumerate() {
        let offset = k * chunk;
        for (dst, src) in out.values.iter_mut().zip(p.values) {
            dst.extend(src);
        }
        out.kept.extend(p.kept.iter().map(|&i| i + offset));
        out.disqualified
            .extend(p.disqualified.iter().map(|&i| i + offset));
    }
    out
}

/// `NOMICS_DUMP_PAIRHMM=<path>`: append every likelihood batch (haplotypes
/// and prepared reads) for offline PairHMM benchmarks. Little-endian
/// records: `u32 n_haps, (u32 len, bases)*, u32 n_reads,
/// (u32 len, bases, quals, ins, del)*`; the gap-continuation penalty is
/// always 10.
fn dump_sink() -> Option<&'static std::sync::Mutex<std::fs::File>> {
    static SINK: std::sync::OnceLock<Option<std::sync::Mutex<std::fs::File>>> =
        std::sync::OnceLock::new();
    SINK.get_or_init(|| {
        let path = std::env::var_os("NOMICS_DUMP_PAIRHMM")?;
        let f = std::fs::File::create(path).expect("create PairHMM dump");
        Some(std::sync::Mutex::new(f))
    })
    .as_ref()
}

fn dump_batch(sink: &std::sync::Mutex<std::fs::File>, haplotypes: &[&[u8]], reads: &[LkRead]) {
    use std::io::Write as _;
    let mut b = Vec::new();
    let u32le = |b: &mut Vec<u8>, x: usize| b.extend_from_slice(&(x as u32).to_le_bytes());
    u32le(&mut b, haplotypes.len());
    for h in haplotypes {
        u32le(&mut b, h.len());
        b.extend_from_slice(h);
    }
    u32le(&mut b, reads.len());
    for r in reads {
        let (q, ins, del) = hmm_qualities(r.bases, r.quals, r.mapq);
        u32le(&mut b, r.bases.len());
        b.extend_from_slice(r.bases);
        b.extend_from_slice(&q);
        b.extend_from_slice(&ins);
        b.extend_from_slice(&del);
    }
    sink.lock()
        .unwrap()
        .write_all(&b)
        .expect("write PairHMM dump");
}
