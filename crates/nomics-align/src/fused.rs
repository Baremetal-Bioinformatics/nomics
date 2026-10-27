//! Fused left+right seed extension (`docs/strategy/gpu-extend-fused-design.md`):
//! one job per seed runs both of `chain2aln`'s band-retry loops and returns
//! every `ksw_extend2` call it made, which the replay then finds by
//! [`ExtKey`].
//!
//! Engines: the CUDA kernel behind [`nomics_gpu::GpuFusedExtend`], and a CPU
//! stand-in ([`extend_cpu`]) that runs the same loops with [`ksw::extend2`].
//! The stand-in is the Mac gate and the reference for the kernel.

use std::time::Instant;

use nomics_gpu::{
    GpuError, GpuExtendResult, GpuFusedExtend, GpuFusedInput, GpuFusedParams, GpuFusedResult,
    KswTimings,
};

use crate::extend::{
    E_DEL, E_INS, ExtKey, MAT, MAX_BAND_TRY, O_DEL, O_INS, PEN_CLIP3, PEN_CLIP5, W, ZDROP,
};
use crate::ksw;

/// bwa-mem's scoring for the fused launch.
pub fn params() -> GpuFusedParams {
    let mut mat = [0i8; 64];
    mat[..MAT.len()].copy_from_slice(&MAT);
    GpuFusedParams {
        m: 5,
        o_del: O_DEL,
        e_del: E_DEL,
        o_ins: O_INS,
        e_ins: E_INS,
        w: W,
        max_band_try: MAX_BAND_TRY,
        end_bonus_left: PEN_CLIP5,
        end_bonus_right: PEN_CLIP3,
        zdrop: ZDROP,
        mat,
    }
}

fn result(x: ksw::Extend) -> GpuExtendResult {
    GpuExtendResult {
        score: x.score,
        qle: x.qle,
        tle: x.tle,
        gtle: x.gtle,
        gscore: x.gscore,
        max_off: x.max_off,
    }
}

fn extend(x: &GpuExtendResult) -> ksw::Extend {
    ksw::Extend {
        score: x.score,
        qle: x.qle,
        tle: x.tle,
        gtle: x.gtle,
        gscore: x.gscore,
        max_off: x.max_off,
    }
}

/// The CPU stand-in for one job: `chain2aln`'s left loop (`prev` starts at
/// `a.score = -1`), then its right loop from `sc0` (`prev` starts at `sc0`).
/// Unused result slots stay zero.
pub fn extend_cpu(p: &GpuFusedParams, job: &GpuFusedInput) -> GpuFusedResult {
    assert!((1..=2).contains(&p.max_band_try));
    let mat = &p.mat[..(p.m * p.m) as usize];
    let call = |q: &[u8], t: &[u8], w: i32, end_bonus: i32, h0: i32| {
        ksw::extend2(
            q.len() as i32,
            q,
            t.len() as i32,
            t,
            p.m,
            mat,
            p.o_del,
            p.e_del,
            p.o_ins,
            p.e_ins,
            w,
            end_bonus,
            p.zdrop,
            h0,
        )
    };
    let mut r = GpuFusedResult::default();
    let mut score = job.h0;
    if !job.lq.is_empty() {
        score = -1;
        for i in 0..p.max_band_try {
            let prev = score;
            let w = p.w << i;
            let x = call(job.lq, job.lt, w, p.end_bonus_left, job.h0);
            r.left[i as usize] = result(x);
            r.n_left = i + 1;
            score = x.score;
            if score == prev || x.max_off < (w >> 1) + (w >> 2) {
                break;
            }
        }
    }
    let sc0 = score;
    r.right_h0 = sc0;
    if !job.rq.is_empty() {
        for i in 0..p.max_band_try {
            let prev = score;
            let w = p.w << i;
            let x = call(job.rq, job.rt, w, p.end_bonus_right, sc0);
            r.right[i as usize] = result(x);
            r.n_right = i + 1;
            score = x.score;
            if score == prev || x.max_off < (w >> 1) + (w >> 2) {
                break;
            }
        }
    }
    r
}

/// Append the calls of seed `seed` of chain `chain`'s fused result as
/// [`ExtKey`]s (call `i` used band `w << i`; left h0 is `h0`, right h0
/// `right_h0`), left to `keys[0]` / `results[0]`, right to index 1.
pub fn push_keys(
    p: &GpuFusedParams,
    chain: u32,
    seed: u32,
    h0: i32,
    r: &GpuFusedResult,
    keys: &mut [Vec<ExtKey>; 2],
    results: &mut [Vec<ksw::Extend>; 2],
) {
    assert!((0..=2).contains(&r.n_left) && (0..=2).contains(&r.n_right));
    for (side, n, calls, h0) in [
        (0, r.n_left, &r.left, h0),
        (1, r.n_right, &r.right, r.right_h0),
    ] {
        for (i, x) in calls[..n as usize].iter().enumerate() {
            keys[side].push(ExtKey {
                chain,
                seed,
                right: side == 1,
                w: p.w << i,
                h0,
            });
            results[side].push(extend(x));
        }
    }
}

/// A fused extension engine.
pub enum FusedEngine {
    /// The CPU stand-in, with its result buffer.
    Cpu(Vec<GpuFusedResult>),
    Gpu(GpuFusedExtend),
}

impl FusedEngine {
    /// The GPU kernel when built with `cuda`, unless `NOMICS_EXTEND_ENGINE=cpu`;
    /// otherwise the stand-in.
    pub fn from_env() -> Self {
        let cpu = std::env::var("NOMICS_EXTEND_ENGINE").is_ok_and(|v| v == "cpu");
        if cfg!(any(feature = "cuda", feature = "metal")) && !cpu {
            Self::Gpu(GpuFusedExtend::new())
        } else {
            Self::Cpu(Vec::new())
        }
    }

    /// Create the GPU engine's stream now and run an empty kernel on it (a
    /// no-op for the stand-in).
    pub fn warm(&mut self) -> Result<(), nomics_gpu::GpuError> {
        match self {
            Self::Cpu(_) => Ok(()),
            Self::Gpu(engine) => engine.warm(),
        }
    }

    /// Bytes of GPU host buffers held (diagnostics; 0 for the stand-in).
    pub fn bytes(&self) -> usize {
        match self {
            Self::Cpu(_) => 0,
            Self::Gpu(engine) => engine.bytes(),
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::Cpu(_) => "cpu",
            Self::Gpu(_) => "gpu",
        }
    }

    /// Every job of `groups`, results in group then iterator order. The
    /// stand-in runs groups in parallel on `threads` threads.
    pub fn compute<'a, I>(
        &mut self,
        p: &GpuFusedParams,
        groups: &[I],
        threads: usize,
    ) -> Result<(&[GpuFusedResult], KswTimings), GpuError>
    where
        I: Iterator<Item = GpuFusedInput<'a>> + Clone + Send + Sync,
    {
        match self {
            Self::Gpu(gpu) => gpu.compute_groups(p, groups, threads),
            Self::Cpu(out) => {
                let start = Instant::now();
                let sizes: Vec<usize> = groups.iter().map(|g| g.clone().count()).collect();
                out.clear();
                out.resize(sizes.iter().sum(), GpuFusedResult::default());
                let mut rest = &mut out[..];
                let mut work = Vec::with_capacity(groups.len());
                for (g, &n) in groups.iter().zip(&sizes) {
                    let (chunk, tail) = rest.split_at_mut(n);
                    rest = tail;
                    work.push((g.clone(), chunk));
                }
                crate::mem::for_each_parallel(threads, work, |(g, chunk)| {
                    for (slot, job) in chunk.iter_mut().zip(g) {
                        *slot = extend_cpu(p, &job);
                    }
                });
                Ok((
                    out,
                    KswTimings {
                        kernel: start.elapsed(),
                        ..KswTimings::default()
                    },
                ))
            }
        }
    }
}

/// `NOMICS_DUMP_FUSED=<dir>`: write every fused launch (jobs with their bases,
/// and the engine's results) to `<dir>/fused.bin`, for `fused_bench` and the
/// kernel gate. Little-endian:
/// - header: magic `NOMFUSD1`, then the params as `nomics_fused_params` lays
///   them out: `i32 m o_del e_del o_ins e_ins w max_band_try end_bonus_left
///   end_bonus_right zdrop`, `mat` (64 i8, `m*m` used);
/// - then one record per job until EOF: `i32 h0, u32 lq_len, lq,
///   u32 lt_len, lt, u32 rq_len, rq, u32 rt_len, rt, i32 n_left n_right
///   right_h0`, then `n_left` left and `n_right` right calls as
///   `i32 score qle tle gtle gscore max_off` each.
///
/// Launches from concurrent lanes are written whole, in completion order.
pub mod dump {
    use std::fs::File;
    use std::io::Write as _;
    use std::ops::Range;
    use std::sync::{Mutex, OnceLock};

    use nomics_gpu::{GpuExtendResult, GpuFusedInput, GpuFusedParams, GpuFusedResult};

    pub const MAGIC: &[u8; 8] = b"NOMFUSD1";

    fn ints(b: &mut Vec<u8>, v: &[i32]) {
        for x in v {
            b.extend_from_slice(&x.to_le_bytes());
        }
    }

    fn seq(b: &mut Vec<u8>, s: &[u8]) {
        b.extend_from_slice(&(s.len() as u32).to_le_bytes());
        b.extend_from_slice(s);
    }

    pub fn encode_header(b: &mut Vec<u8>, p: &GpuFusedParams) {
        b.extend_from_slice(MAGIC);
        ints(
            b,
            &[
                p.m,
                p.o_del,
                p.e_del,
                p.o_ins,
                p.e_ins,
                p.w,
                p.max_band_try,
                p.end_bonus_left,
                p.end_bonus_right,
                p.zdrop,
            ],
        );
        b.extend(p.mat.iter().map(|&x| x as u8));
    }

    pub fn encode_job(b: &mut Vec<u8>, job: &GpuFusedInput, r: &GpuFusedResult) {
        ints(b, &[job.h0]);
        for s in [job.lq, job.lt, job.rq, job.rt] {
            seq(b, s);
        }
        ints(b, &[r.n_left, r.n_right, r.right_h0]);
        let calls = r.left[..r.n_left as usize]
            .iter()
            .chain(&r.right[..r.n_right as usize]);
        for x in calls {
            ints(b, &[x.score, x.qle, x.tle, x.gtle, x.gscore, x.max_off]);
        }
    }

    /// One dumped job: ranges into [`Dump::bytes`] and the recorded result
    /// (unused slots zero).
    pub struct Job {
        pub h0: i32,
        pub lq: Range<usize>,
        pub lt: Range<usize>,
        pub rq: Range<usize>,
        pub rt: Range<usize>,
        pub result: GpuFusedResult,
    }

    pub struct Dump {
        pub params: GpuFusedParams,
        pub bytes: Vec<u8>,
        pub jobs: Vec<Job>,
    }

    impl Dump {
        pub fn parse(bytes: Vec<u8>) -> Self {
            assert!(bytes.len() >= 8 && &bytes[..8] == MAGIC, "not a fused.bin");
            let mut at = 8usize;
            let take = |at: &mut usize, n: usize| {
                let r = *at..*at + n;
                assert!(r.end <= bytes.len(), "truncated fused.bin");
                *at += n;
                r
            };
            let int = |at: &mut usize| {
                let r = take(at, 4);
                i32::from_le_bytes(bytes[r].try_into().unwrap())
            };
            let h: [i32; 10] = std::array::from_fn(|_| int(&mut at));
            let mut mat = [0i8; 64];
            for (m, &x) in mat.iter_mut().zip(&bytes[take(&mut at, 64)]) {
                *m = x as i8;
            }
            let params = GpuFusedParams {
                m: h[0],
                o_del: h[1],
                e_del: h[2],
                o_ins: h[3],
                e_ins: h[4],
                w: h[5],
                max_band_try: h[6],
                end_bonus_left: h[7],
                end_bonus_right: h[8],
                zdrop: h[9],
                mat,
            };
            let mut jobs = Vec::new();
            while at < bytes.len() {
                let h0 = int(&mut at);
                let mut s = || {
                    let n = int(&mut at) as u32 as usize;
                    take(&mut at, n)
                };
                let (lq, lt, rq, rt) = (s(), s(), s(), s());
                let mut result = GpuFusedResult {
                    n_left: int(&mut at),
                    n_right: int(&mut at),
                    right_h0: int(&mut at),
                    ..GpuFusedResult::default()
                };
                assert!((0..=2).contains(&result.n_left) && (0..=2).contains(&result.n_right));
                let mut call = || {
                    let v: [i32; 6] = std::array::from_fn(|_| int(&mut at));
                    GpuExtendResult {
                        score: v[0],
                        qle: v[1],
                        tle: v[2],
                        gtle: v[3],
                        gscore: v[4],
                        max_off: v[5],
                    }
                };
                for i in 0..result.n_left as usize {
                    result.left[i] = call();
                }
                for i in 0..result.n_right as usize {
                    result.right[i] = call();
                }
                jobs.push(Job {
                    h0,
                    lq,
                    lt,
                    rq,
                    rt,
                    result,
                });
            }
            Self {
                params,
                bytes,
                jobs,
            }
        }

        pub fn input(&self, j: &Job) -> GpuFusedInput<'_> {
            GpuFusedInput {
                lq: &self.bytes[j.lq.clone()],
                lt: &self.bytes[j.lt.clone()],
                rq: &self.bytes[j.rq.clone()],
                rt: &self.bytes[j.rt.clone()],
                h0: j.h0,
            }
        }
    }

    fn sink(p: &GpuFusedParams) -> Option<&'static Mutex<File>> {
        static SINK: OnceLock<Option<(GpuFusedParams, Mutex<File>)>> = OnceLock::new();
        let (params, file) = SINK
            .get_or_init(|| {
                let dir = std::path::PathBuf::from(std::env::var_os("NOMICS_DUMP_FUSED")?);
                let mut file = File::create(dir.join("fused.bin")).expect("create fused dump");
                let mut b = Vec::new();
                encode_header(&mut b, p);
                file.write_all(&b).expect("write fused dump");
                Some((*p, Mutex::new(file)))
            })
            .as_ref()?;
        assert_eq!(params, p, "fused dump: launches differ in params");
        Some(file)
    }

    fn env_u64(name: &str) -> Option<u64> {
        let v = std::env::var(name).ok()?;
        Some(v.parse().unwrap_or_else(|_| panic!("{name}={v}: not an integer")))
    }

    /// Append one launch (`groups` of jobs, `results` in the same order) if
    /// `NOMICS_DUMP_FUSED` is set. `NOMICS_DUMP_FUSED_EVERY=k` writes only
    /// every k-th launch (0, k, 2k, ... in completion order), to sample a
    /// whole-genome chunk; `NOMICS_DUMP_FUSED_MAX_JOBS=n` stops writing once
    /// n jobs are in the file. Launches are always written whole.
    pub fn launch<'a, I>(p: &GpuFusedParams, groups: &[I], results: &[GpuFusedResult])
    where
        I: Iterator<Item = GpuFusedInput<'a>> + Clone,
    {
        static LIMITS: OnceLock<(u64, u64)> = OnceLock::new();
        // (launches seen, jobs written), updated under the file lock.
        static COUNTS: Mutex<(u64, u64)> = Mutex::new((0, 0));
        let Some(file) = sink(p) else {
            return;
        };
        let (every, max_jobs) = *LIMITS.get_or_init(|| {
            (
                env_u64("NOMICS_DUMP_FUSED_EVERY").unwrap_or(1).max(1),
                env_u64("NOMICS_DUMP_FUSED_MAX_JOBS").unwrap_or(u64::MAX),
            )
        });
        let mut file = file.lock().unwrap();
        {
            let mut counts = COUNTS.lock().unwrap();
            let launch = counts.0;
            counts.0 += 1;
            if !launch.is_multiple_of(every) || counts.1 >= max_jobs {
                return;
            }
            counts.1 += results.len() as u64;
        }
        let mut results = results.iter();
        let mut b = Vec::new();
        for g in groups {
            b.clear();
            for job in g.clone() {
                encode_job(&mut b, &job, results.next().expect("fused dump: results"));
            }
            file.write_all(&b).expect("write fused dump");
        }
        assert!(results.next().is_none(), "fused dump: results");
    }
}
