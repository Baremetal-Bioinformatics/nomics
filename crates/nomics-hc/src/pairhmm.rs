//! GATK 4.6.2.0's Java `LoglessPairHMM` (spec §6.3), bit-exact: all double
//! arithmetic in Java's evaluation order, no FMA, transcendental functions
//! through `jmath`. Gate: `tests/pairhmm.rs` against the reference jar
//! (`reference/oracle/hc/GatkPairHmm.java`).
//!
//! Because HaplotypeCaller always passes `recacheReadValues = true`, every
//! (read, haplotype) likelihood is computed from scratch: a pure function of
//! the read arrays and the haplotype bases.

use std::sync::OnceLock;

use crate::jmath;
use crate::tables::{MAX_QUAL, approximate_log10_sum_log10, qual_to_error_prob, qual_to_prob};

const TRISTATE_CORRECTION: f64 = 3.0;

/// `PairHMMModel.matchToMatchProb` cache, indexed `max*(max+1)/2 + min`.
fn match_to_match_table() -> &'static [f64] {
    static T: OnceLock<Vec<f64>> = OnceLock::new();
    T.get_or_init(|| {
        let inv_ln10 = 1.0 / jmath::log(10.0);
        let n = ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1;
        let mut t = vec![0.0; n];
        let mut offset = 0usize;
        for i in 0..=MAX_QUAL {
            for j in 0..=i {
                let log10_sum = approximate_log10_sum_log10(-0.1 * i as f64, -0.1 * j as f64);
                let m2m_log10 = jmath::log1p(-(jmath::pow(10.0, log10_sum)).min(1.0)) * inv_ln10;
                t[offset + j] = jmath::pow(10.0, m2m_log10);
            }
            offset += i + 1;
        }
        t
    })
}

/// `PairHMMModel.matchToMatchProb(insQual, delQual)` (qualities <= 254).
fn match_to_match(ins: u8, del: u8) -> f64 {
    let (min, max) = if ins <= del {
        (ins as usize, del as usize)
    } else {
        (del as usize, ins as usize)
    };
    match_to_match_table()[((max * (max + 1)) >> 1) + min]
}

pub fn initial_condition() -> (f64, f64) {
    static C: OnceLock<(f64, f64)> = OnceLock::new();
    *C.get_or_init(|| {
        let c = jmath::pow(2.0, 1020.0);
        (c, jmath::log10(c))
    })
}

/// The seven per-read-row values consumed by the PairHMM recurrence, in the
/// order expected by the CUDA experiment: `m2m`, `m2i`, `m2d`, `i2m`, `i2i`,
/// `p_match`, and `p_mismatch`.
///
/// Keeping this precomputation here makes the GPU path use the same lookup
/// tables and `jmath` operations as the gated CPU implementation.
pub fn row_params(ins_q: u8, del_q: u8, q: u8, gcp: u8) -> [f64; 7] {
    [
        match_to_match(ins_q, del_q),
        qual_to_error_prob(ins_q),
        qual_to_error_prob(del_q),
        qual_to_prob(gcp),
        qual_to_error_prob(gcp),
        qual_to_prob(q),
        qual_to_error_prob(q) / TRISTATE_CORRECTION,
    ]
}

/// Reusable rows: the recurrence only reads the previous row (and the
/// current row's previous column), so two rows per state suffice. Every
/// cell keeps Java's operand order.
#[derive(Default)]
pub struct PairHmm {
    m: [Vec<f64>; 2],
    ins: [Vec<f64>; 2],
    del: [Vec<f64>; 2],
}

impl PairHmm {
    /// `computeReadLikelihoodGivenHaplotypeLog10(hap, read, quals, ins, del,
    /// gcp, recacheReadValues = true, …)`.
    #[allow(clippy::too_many_arguments)]
    pub fn log10_likelihood(
        &mut self,
        hap: &[u8],
        read: &[u8],
        quals: &[u8],
        ins_q: &[u8],
        del_q: &[u8],
        gcp: &[u8],
    ) -> f64 {
        let hl = hap.len() + 1;
        let (init, init_log10) = initial_condition();
        let initial_value = init / hap.len() as f64;
        for k in 0..2 {
            for v in [&mut self.m[k], &mut self.ins[k], &mut self.del[k]] {
                v.clear();
                v.resize(hl, 0.0);
            }
        }
        // Row 0: deletions start anywhere along the haplotype.
        self.del[0].fill(initial_value);
        let mut cur = 0usize;
        for (i, (&x, &q)) in read.iter().zip(quals).enumerate() {
            let prev = cur;
            cur ^= 1;
            // PairHMMModel.qualToTransProbs for read base i.
            let m2m = match_to_match(ins_q[i], del_q[i]);
            let m2i = qual_to_error_prob(ins_q[i]);
            let m2d = qual_to_error_prob(del_q[i]);
            let i2m = qual_to_prob(gcp[i]);
            let i2i = qual_to_error_prob(gcp[i]);
            let d2d = i2i;
            // initializePriors
            let p_match = qual_to_prob(q);
            let p_mismatch = qual_to_error_prob(q) / TRISTATE_CORRECTION;
            let [m0, m1] = &mut self.m;
            let [i0, i1] = &mut self.ins;
            let [d0, d1] = &mut self.del;
            let (mp, mc) = if prev == 0 { (&*m0, m1) } else { (&*m1, m0) };
            let (ip, ic) = if prev == 0 { (&*i0, i1) } else { (&*i1, i0) };
            let (dp, dc) = if prev == 0 { (&*d0, d1) } else { (&*d1, d0) };
            mc[0] = 0.0;
            ic[0] = 0.0;
            dc[0] = 0.0;
            for j in 1..hl {
                let y = hap[j - 1];
                let prior = if x == y || x == b'N' || y == b'N' {
                    p_match
                } else {
                    p_mismatch
                };
                mc[j] = prior * (mp[j - 1] * m2m + ip[j - 1] * i2m + dp[j - 1] * i2m);
                ic[j] = mp[j] * m2i + ip[j] * i2i;
                dc[j] = mc[j - 1] * m2d + dc[j - 1] * d2d;
            }
        }
        let mut sum = 0.0;
        for j in 1..hl {
            sum += self.m[cur][j] + self.ins[cur][j];
        }
        jmath::log10(sum) - init_log10
    }
}

/// Haplotypes per vector pass.
pub const LANES: usize = 8;
type V = [f64; LANES];

/// Rows for `log10_likelihood_lanes`.
#[derive(Default)]
pub struct PairHmmLanes {
    m: [Vec<V>; 2],
    ins: [Vec<V>; 2],
    del: [Vec<V>; 2],
    hap: Vec<[u8; LANES]>,
    /// `hap` widened for the NEON row.
    hap64: Vec<[u64; LANES]>,
}

impl PairHmmLanes {
    /// `PairHmm::log10_likelihood` for up to `LANES` haplotypes against one
    /// read at once. Each lane performs exactly the scalar recurrence, so
    /// the results are bit-identical; shorter haplotypes are padded, and
    /// columns past a lane's own length never feed its earlier columns.
    #[allow(clippy::too_many_arguments)]
    pub fn log10_likelihoods(
        &mut self,
        haps: &[&[u8]],
        read: &[u8],
        quals: &[u8],
        ins_q: &[u8],
        del_q: &[u8],
        gcp: &[u8],
        out: &mut [f64],
    ) {
        self.likelihoods(haps, read, quals, ins_q, del_q, gcp, out, false);
    }

    /// [`Self::log10_likelihoods`] on the portable array lanes, also on
    /// aarch64 (where the NEON row is used otherwise); for equivalence tests.
    #[allow(clippy::too_many_arguments)]
    pub fn log10_likelihoods_portable(
        &mut self,
        haps: &[&[u8]],
        read: &[u8],
        quals: &[u8],
        ins_q: &[u8],
        del_q: &[u8],
        gcp: &[u8],
        out: &mut [f64],
    ) {
        self.likelihoods(haps, read, quals, ins_q, del_q, gcp, out, true);
    }

    #[allow(clippy::too_many_arguments)]
    fn likelihoods(
        &mut self,
        haps: &[&[u8]],
        read: &[u8],
        quals: &[u8],
        ins_q: &[u8],
        del_q: &[u8],
        gcp: &[u8],
        out: &mut [f64],
        portable: bool,
    ) {
        let neon = cfg!(target_arch = "aarch64") && !portable;
        let n = haps.len();
        assert!((1..=LANES).contains(&n) && out.len() == n);
        let max_len = haps.iter().map(|h| h.len()).max().unwrap();
        let hl = max_len + 1;
        let (init, init_log10) = initial_condition();
        self.hap.clear();
        self.hap.resize(max_len, [0; LANES]);
        let mut initial: V = [0.0; LANES];
        for l in 0..LANES {
            let h = haps[l.min(n - 1)];
            for (j, &b) in h.iter().enumerate() {
                self.hap[j][l] = b;
            }
            initial[l] = init / h.len() as f64;
        }
        for k in 0..2 {
            for v in [&mut self.m[k], &mut self.ins[k], &mut self.del[k]] {
                v.clear();
                v.resize(hl, [0.0; LANES]);
            }
        }
        self.del[0].fill(initial);
        if neon {
            self.hap64.clear();
            self.hap64
                .extend(self.hap.iter().map(|y| y.map(|b| b as u64)));
        }
        let mut cur = 0usize;
        for (i, (&x, &q)) in read.iter().zip(quals).enumerate() {
            let prev = cur;
            cur ^= 1;
            let m2m = match_to_match(ins_q[i], del_q[i]);
            let m2i = qual_to_error_prob(ins_q[i]);
            let m2d = qual_to_error_prob(del_q[i]);
            let i2m = qual_to_prob(gcp[i]);
            let i2i = qual_to_error_prob(gcp[i]);
            let d2d = i2i;
            let p_match = qual_to_prob(q);
            let p_mismatch = qual_to_error_prob(q) / TRISTATE_CORRECTION;
            let [m0, m1] = &mut self.m;
            let [i0, i1] = &mut self.ins;
            let [d0, d1] = &mut self.del;
            let (mp, mc) = if prev == 0 { (&*m0, m1) } else { (&*m1, m0) };
            let (ip, ic) = if prev == 0 { (&*i0, i1) } else { (&*i1, i0) };
            let (dp, dc) = if prev == 0 { (&*d0, d1) } else { (&*d1, d0) };
            mc[0] = [0.0; LANES];
            ic[0] = [0.0; LANES];
            dc[0] = [0.0; LANES];
            #[cfg(target_arch = "aarch64")]
            if neon {
                nomics_simd::pairhmm::row(
                    x,
                    nomics_simd::pairhmm::RowParams {
                        m2m,
                        m2i,
                        m2d,
                        i2m,
                        i2i,
                        p_match,
                        p_mismatch,
                    },
                    &self.hap64,
                    (mp, ip, dp),
                    (mc, ic, dc),
                    n,
                );
                continue;
            }
            for j in 1..hl {
                let y = self.hap[j - 1];
                let (mpj, ipj, dpj) = (mp[j - 1], ip[j - 1], dp[j - 1]);
                let (mpu, ipu) = (mp[j], ip[j]);
                let (mcl, dcl) = (mc[j - 1], dc[j - 1]);
                let mut mv = [0.0; LANES];
                let mut iv = [0.0; LANES];
                let mut dv = [0.0; LANES];
                for l in 0..LANES {
                    let prior = if x == y[l] || x == b'N' || y[l] == b'N' {
                        p_match
                    } else {
                        p_mismatch
                    };
                    mv[l] = prior * (mpj[l] * m2m + ipj[l] * i2m + dpj[l] * i2m);
                    iv[l] = mpu[l] * m2i + ipu[l] * i2i;
                    dv[l] = mcl[l] * m2d + dcl[l] * d2d;
                }
                // `dv` uses this row's previous column, computed last iteration.
                mc[j] = mv;
                ic[j] = iv;
                dc[j] = dv;
            }
        }
        for (l, o) in out.iter_mut().enumerate() {
            let mut sum = 0.0;
            for j in 1..=haps[l].len() {
                sum += self.m[cur][j][l] + self.ins[cur][j][l];
            }
            *o = jmath::log10(sum) - init_log10;
        }
    }
}
