//! Read preparation for the PairHMM (spec §6.2): GATK 4.6.2.0's PCR indel
//! error model (`--pcr-indel-model CONSERVATIVE`) and quality capping, as
//! `PairHMMLikelihoodCalculationEngine.modifyReadQualities` applies them.
//! Gate: `tests/readprep.rs` against the reference jar.

use std::sync::OnceLock;

use crate::jmath;
use crate::tables::fast_round;

/// `ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH`.
pub const MAX_REPEAT_LENGTH: usize = 20;
const MAX_STR_UNIT_LENGTH: usize = 8;
const INITIAL_QSCORE: f64 = 40.0;
const MIN_ADJUSTED_QSCORE: i32 = 10;
/// `PCRErrorModel.CONSERVATIVE`.
const RATE_FACTOR: f64 = 3.0;
/// Default insertion/deletion quality when the read has no BI/BD tags.
pub const DEFAULT_INDEL_QUAL: u8 = 45;
/// `PairHMM.BASE_QUALITY_SCORE_THRESHOLD`.
const BASE_QUALITY_SCORE_THRESHOLD: u8 = 18;
/// `QualityUtils.MIN_USABLE_Q_SCORE`.
const MIN_USABLE_Q_SCORE: u8 = 6;

/// `pcrIndelErrorModelCache[r] = max(10, fastRound(40 - exp(r / (3π)) + 1))`.
pub fn pcr_cache() -> &'static [u8; MAX_REPEAT_LENGTH + 1] {
    static C: OnceLock<[u8; MAX_REPEAT_LENGTH + 1]> = OnceLock::new();
    C.get_or_init(|| {
        let mut c = [0u8; MAX_REPEAT_LENGTH + 1];
        for (i, v) in c.iter_mut().enumerate() {
            let q = fast_round(
                INITIAL_QSCORE - jmath::exp(i as f64 / (RATE_FACTOR * std::f64::consts::PI)) + 1.0,
            );
            *v = q.max(MIN_ADJUSTED_QSCORE) as u8;
        }
        c
    })
}

/// `GATKVariantContextUtils.findNumberOfRepetitions(unit, test, leadingRepeats)`
/// on slices.
fn repetitions(unit: &[u8], test: &[u8], leading: bool) -> usize {
    if test.is_empty() {
        return 0;
    }
    let ul = unit.len();
    if test.len() < ul {
        return 0;
    }
    let diff = test.len() - ul;
    let mut n = 0;
    if leading {
        let mut start = 0;
        while start <= diff {
            if &test[start..start + ul] != unit {
                return n;
            }
            n += 1;
            start += ul;
        }
    } else {
        let mut start = diff as isize;
        while start >= 0 {
            let s = start as usize;
            if &test[s..s + ul] != unit {
                return n;
            }
            n += 1;
            start -= ul as isize;
        }
    }
    n
}

/// `ReadLikelihoodCalculationEngine.findTandemRepeatUnits(bases, offset).right`.
pub fn tandem_repeat_length(bases: &[u8], offset: usize) -> usize {
    let mut max_bw = 0usize;
    let mut best_bw: &[u8] = &bases[offset..offset + 1];
    for s in 1..=MAX_STR_UNIT_LENGTH {
        if offset + 1 < s {
            break;
        }
        // backward unit ending at offset, counted against bases[0..=offset]
        max_bw = repetitions(
            &bases[offset + 1 - s..offset + 1],
            &bases[..offset + 1],
            false,
        );
        if max_bw > 1 {
            best_bw = &bases[offset + 1 - s..offset + 1];
            break;
        }
    }
    let mut max_rl = max_bw;
    if offset < bases.len() - 1 {
        let mut best_fw: &[u8] = &bases[offset + 1..offset + 2];
        let mut max_fw = 0usize;
        for s in 1..=MAX_STR_UNIT_LENGTH {
            if offset + s + 1 > bases.len() {
                break;
            }
            max_fw = repetitions(
                &bases[offset + 1..offset + 1 + s],
                &bases[offset + 1..],
                true,
            );
            if max_fw > 1 {
                best_fw = &bases[offset + 1..offset + s + 1];
                break;
            }
        }
        if best_fw == best_bw {
            max_rl = max_bw + max_fw;
        } else {
            let bw = repetitions(best_fw, &bases[..offset + 1], false);
            max_rl = max_fw + bw;
        }
    }
    max_rl.min(MAX_REPEAT_LENGTH)
}

/// Qualities the PairHMM sees for one read: (base quals, insertion quals,
/// deletion quals), from `modifyReadQualities` with soft clips kept.
pub fn hmm_qualities(bases: &[u8], quals: &[u8], mapq: u8) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let n = bases.len();
    let mut q = quals.to_vec();
    let mut ins = vec![DEFAULT_INDEL_QUAL; n];
    let mut del = vec![DEFAULT_INDEL_QUAL; n];
    let cache = pcr_cache();
    for i in 1..n {
        let rl = tandem_repeat_length(bases, i - 1);
        ins[i - 1] = ins[i - 1].min(cache[rl]);
        del[i - 1] = del[i - 1].min(cache[rl]);
    }
    for i in 0..n {
        q[i] = q[i].min(mapq);
        if q[i] < BASE_QUALITY_SCORE_THRESHOLD {
            q[i] = MIN_USABLE_Q_SCORE;
        }
        if ins[i] < MIN_USABLE_Q_SCORE {
            ins[i] = MIN_USABLE_Q_SCORE;
        }
        if del[i] < MIN_USABLE_Q_SCORE {
            del[i] = MIN_USABLE_Q_SCORE;
        }
    }
    (q, ins, del)
}
