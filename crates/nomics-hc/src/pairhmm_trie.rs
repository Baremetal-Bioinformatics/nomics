//! The Java PairHMM over all (read, haplotype) pairs of a region, computing
//! each distinct haplotype prefix once.
//!
//! GATK 4.6.2.0 computes every pair from scratch (its `hapStartIndex` reuse
//! is unreachable: `recacheReadValues` is always true). For a fixed read,
//! cell (i, j) depends only on the read, haplotype bases `0..j` and, through
//! row 0 of D (`2^1020 / hapLen`), the haplotype length. So for haplotypes of
//! equal length sharing a prefix of length `p`, columns `0..=p` hold the same
//! values, and so does the running final sum up to column `p`: computing them
//! once is exact. Every computed cell keeps Java's operand order and the final
//! sum stays sequential over columns, so results are bit-identical to
//! [`crate::pairhmm::PairHmm::log10_likelihood`].
//!
//! Haplotypes are grouped by length and sorted; haplotype k starts at column
//! `s[k]`, its common prefix with the previous one, from the column state of
//! the latest earlier haplotype that computed column `s[k]`. Lanes hold
//! [`LANES`] reads (each with its own per-row parameters); the haplotype base
//! of a column is shared by all lanes.

use crate::jmath;
use crate::pairhmm::{initial_condition, row_params};

/// Reads per pass.
pub const LANES: usize = 8;
type V = [f64; LANES];

/// One read's PairHMM inputs.
pub struct HmmRead<'a> {
    pub bases: &'a [u8],
    pub quals: &'a [u8],
    pub ins: &'a [u8],
    pub del: &'a [u8],
    pub gcp: &'a [u8],
}

/// Per-row parameters of a group of reads, lane-major per row.
struct Rows {
    /// Read base per row and lane (0 past a read's end).
    x: Vec<[u8; LANES]>,
    /// `row_params` per row: m2m, m2i, m2d, i2m, i2i, p_match, p_mismatch.
    p: Vec<[V; 7]>,
}

/// Column state of all rows at one column, saved for a later haplotype.
#[derive(Clone)]
struct Checkpoint {
    m: Vec<V>,
    i: Vec<V>,
    d: Vec<V>,
    /// Each lane's partial final sum up to this column.
    sum: V,
}

/// `values[h][r]`: log10 likelihood of read `r` given haplotype `h`.
pub fn likelihoods(haps: &[&[u8]], reads: &[HmmRead]) -> Vec<Vec<f64>> {
    let mut values = vec![vec![0.0; reads.len()]; haps.len()];
    if reads.is_empty() || haps.is_empty() {
        return values;
    }
    // Reads by length (longest first) so a pass pads few rows.
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&r| std::cmp::Reverse(reads[r].bases.len()));
    // Haplotypes grouped by length, sorted within a group.
    let mut by_len: Vec<usize> = (0..haps.len()).collect();
    by_len.sort_by(|&a, &b| haps[a].len().cmp(&haps[b].len()).then(haps[a].cmp(haps[b])));
    let groups: Vec<&[usize]> = by_len
        .chunk_by(|&a, &b| haps[a].len() == haps[b].len())
        .collect();

    for lanes in order.chunks(LANES) {
        let rows = lane_rows(reads, lanes);
        let lens: Vec<usize> = lanes.iter().map(|&r| reads[r].bases.len()).collect();
        for group in &groups {
            let sums = group_pass(haps, group, &rows, &lens);
            for (g, &h) in group.iter().enumerate() {
                for (l, &r) in lanes.iter().enumerate() {
                    values[h][r] = sums[g][l];
                }
            }
        }
    }
    values
}

fn lane_rows(reads: &[HmmRead], lanes: &[usize]) -> Rows {
    let n = lanes.iter().map(|&r| reads[r].bases.len()).max().unwrap();
    let mut rows = Rows {
        x: vec![[0; LANES]; n],
        p: vec![[[0.0; LANES]; 7]; n],
    };
    for (l, &r) in lanes.iter().enumerate() {
        let read = &reads[r];
        for i in 0..read.bases.len() {
            rows.x[i][l] = read.bases[i];
            let p = row_params(read.ins[i], read.del[i], read.quals[i], read.gcp[i]);
            for k in 0..7 {
                rows.p[i][k][l] = p[k];
            }
        }
    }
    rows
}

/// Final log10 likelihoods of every haplotype of one equal-length group
/// (`group`, sorted) for one pass of reads: `[haplotype in group][lane]`.
fn group_pass(haps: &[&[u8]], group: &[usize], rows: &Rows, lens: &[usize]) -> Vec<V> {
    let hl = haps[group[0]].len();
    let n_rows = rows.x.len();
    let (init, init_log10) = initial_condition();
    let initial = init / hl as f64;

    // Start column of each haplotype, and which earlier haplotype supplies it.
    let start: Vec<usize> = (0..group.len())
        .map(|k| {
            if k == 0 {
                0
            } else {
                common_prefix(haps[group[k - 1]], haps[group[k]])
            }
        })
        .collect();
    let mut source = vec![usize::MAX; group.len()];
    let mut stack: Vec<usize> = Vec::new();
    for k in 0..group.len() {
        while let Some(&j) = stack.last() {
            if start[j] > start[k] {
                stack.pop();
            } else {
                break;
            }
        }
        if let Some(&j) = stack.last() {
            source[k] = j;
        }
        stack.push(k);
    }
    // Checkpoints each haplotype must save: (column, consumer).
    let mut saves: Vec<Vec<(usize, usize)>> = vec![Vec::new(); group.len()];
    for k in 0..group.len() {
        if start[k] > 0 {
            saves[source[k]].push((start[k], k));
        }
    }
    let mut checkpoints: Vec<Option<Checkpoint>> = vec![None; group.len()];

    let mut out = vec![[0.0; LANES]; group.len()];
    let zero: V = [0.0; LANES];
    // Row buffers: previous and current row of M, I, D over columns 0..=hl.
    let mut m = [vec![zero; hl + 1], vec![zero; hl + 1]];
    let mut ins = [vec![zero; hl + 1], vec![zero; hl + 1]];
    let mut del = [vec![zero; hl + 1], vec![zero; hl + 1]];
    for k in 0..group.len() {
        let hap = haps[group[k]];
        let s = start[k];
        let from = if s > 0 {
            checkpoints[k]
                .take()
                .expect("checkpoint saved by an earlier haplotype")
        } else {
            Checkpoint {
                m: vec![zero; n_rows + 1],
                i: vec![zero; n_rows + 1],
                d: vec![zero; n_rows + 1],
                sum: zero,
            }
        };
        let mut saving: Vec<(usize, usize, Checkpoint)> = saves[k]
            .iter()
            .map(|&(c, consumer)| {
                (
                    c,
                    consumer,
                    Checkpoint {
                        m: vec![zero; n_rows + 1],
                        i: vec![zero; n_rows + 1],
                        d: vec![zero; n_rows + 1],
                        sum: zero,
                    },
                )
            })
            .collect();
        // Row 0: deletions start anywhere along the haplotype.
        m[0][s..].fill(zero);
        ins[0][s..].fill(zero);
        del[0][s..].fill([initial; LANES]);
        if s > 0 {
            m[0][s] = from.m[0];
            ins[0][s] = from.i[0];
            del[0][s] = from.d[0];
        }
        for (c, _, cp) in saving.iter_mut() {
            cp.m[0] = m[0][*c];
            cp.i[0] = ins[0][*c];
            cp.d[0] = del[0][*c];
        }
        // Each lane's running sum over its last row, from column s.
        let mut sum = from.sum;
        let mut cur = 0usize;
        for i in 0..n_rows {
            let prev = cur;
            cur ^= 1;
            let [m0, m1] = &mut m;
            let [i0, i1] = &mut ins;
            let [d0, d1] = &mut del;
            let (mp, mc) = if prev == 0 { (&*m0, m1) } else { (&*m1, m0) };
            let (ip, ic) = if prev == 0 { (&*i0, i1) } else { (&*i1, i0) };
            let (dp, dc) = if prev == 0 { (&*d0, d1) } else { (&*d1, d0) };
            // Column s of this row: column 0 (zeros) or the checkpoint.
            mc[s] = from.m[i + 1];
            ic[s] = from.i[i + 1];
            dc[s] = from.d[i + 1];
            #[cfg(target_arch = "aarch64")]
            nomics_simd::pairhmm_reads::row(
                &rows.x[i],
                &rows.p[i],
                hap,
                s,
                (mp, ip, dp),
                (mc, ic, dc),
            );
            #[cfg(not(target_arch = "aarch64"))]
            row(&rows.x[i], &rows.p[i], hap, s, (mp, ip, dp), (mc, ic, dc));
            for (c, _, cp) in saving.iter_mut() {
                cp.m[i + 1] = mc[*c];
                cp.i[i + 1] = ic[*c];
                cp.d[i + 1] = dc[*c];
            }
            // Lanes whose read ends on this row: their final sum, sequential
            // over columns, recording it at every checkpoint column.
            for (l, &len) in lens.iter().enumerate() {
                if len == i + 1 {
                    let mut acc = sum[l];
                    let mut next_save = 0;
                    let mut saves_by_col: Vec<usize> = (0..saving.len()).collect();
                    saves_by_col.sort_by_key(|&t| saving[t].0);
                    while next_save < saves_by_col.len() && saving[saves_by_col[next_save]].0 == s {
                        saving[saves_by_col[next_save]].2.sum[l] = acc;
                        next_save += 1;
                    }
                    for j in s + 1..=hl {
                        acc += mc[j][l] + ic[j][l];
                        while next_save < saves_by_col.len()
                            && saving[saves_by_col[next_save]].0 == j
                        {
                            saving[saves_by_col[next_save]].2.sum[l] = acc;
                            next_save += 1;
                        }
                    }
                    sum[l] = acc;
                }
            }
        }
        for (l, o) in out[k].iter_mut().enumerate() {
            if l < lens.len() {
                *o = jmath::log10(sum[l]) - init_log10;
            }
        }
        for (_, consumer, cp) in saving {
            checkpoints[consumer] = Some(cp);
        }
    }
    out
}

fn common_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
/// One read row over columns `s + 1..=hap.len()`, every lane exactly the
/// scalar recurrence. `cur[s]` must already hold column `s`.
fn row(
    x: &[u8; LANES],
    p: &[V; 7],
    hap: &[u8],
    s: usize,
    prev: (&[V], &[V], &[V]),
    cur: (&mut [V], &mut [V], &mut [V]),
) {
    let (mp, ip, dp) = prev;
    let (mc, ic, dc) = cur;
    let [m2m, m2i, m2d, i2m, i2i, p_match, p_mismatch] = p;
    for j in s + 1..=hap.len() {
        let y = hap[j - 1];
        let (mpj, ipj, dpj) = (mp[j - 1], ip[j - 1], dp[j - 1]);
        let (mpu, ipu) = (mp[j], ip[j]);
        let (mcl, dcl) = (mc[j - 1], dc[j - 1]);
        let mut mv = [0.0; LANES];
        let mut iv = [0.0; LANES];
        let mut dv = [0.0; LANES];
        for l in 0..LANES {
            let prior = if x[l] == y || x[l] == b'N' || y == b'N' {
                p_match[l]
            } else {
                p_mismatch[l]
            };
            mv[l] = prior * (mpj[l] * m2m[l] + ipj[l] * i2m[l] + dpj[l] * i2m[l]);
            iv[l] = mpu[l] * m2i[l] + ipu[l] * i2i[l];
            dv[l] = mcl[l] * m2d[l] + dcl[l] * i2i[l];
        }
        mc[j] = mv;
        ic[j] = iv;
        dc[j] = dv;
    }
}
