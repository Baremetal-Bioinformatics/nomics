//! Host side of the fused left+right `ksw_extend2` launch (`shim/fused.cu`,
//! design in `docs/strategy/gpu-extend-fused-design.md`). One job per seed:
//! the kernel runs BWA's two band-retry loops and returns every call it made.
//! No alignment arithmetic lives here; jobs keep their input order.
use crate::backend::{Pending, Stream};
use crate::extend::{ensure, parallel, typed, typed_mut};
use crate::{EXTEND_MAX_QUERY, GpuError, GpuExtendResult, KswTimings};
use nomics_mem::HostBuf;
use std::time::{Duration, Instant};

/// `nomics_fused_params`: scoring shared by every job of a launch.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GpuFusedParams {
    /// Alphabet size, `m * m <= 64`.
    pub m: i32,
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    /// First band; band `i` is `w << i`.
    pub w: i32,
    /// At most 2 (the result holds two calls per side).
    pub max_band_try: i32,
    pub end_bonus_left: i32,
    pub end_bonus_right: i32,
    pub zdrop: i32,
    /// `m * m` used, the rest zero.
    pub mat: [i8; 64],
}

/// `nomics_fused_job`. Offsets index the launch's query (`lq`, `rq`) and
/// target (`lt`, `rt`) arenas.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GpuFusedJob {
    pub lq_off: u64,
    pub lt_off: u64,
    pub rq_off: u64,
    pub rt_off: u64,
    /// 0: no left (`qbeg == 0`).
    pub lq_len: u32,
    pub lt_len: u32,
    /// 0: no right on the GPU.
    pub rq_len: u32,
    pub rt_len: u32,
    /// Seed `len * A`.
    pub h0: i32,
    pub pad: i32,
}

/// `nomics_fused_result`: call `i` on a side used band `w << i`; left calls
/// used `h0 = job.h0`, right calls `right_h0`. Slots at or past `n_left` /
/// `n_right` are unspecified.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GpuFusedResult {
    pub left: [GpuExtendResult; 2],
    pub right: [GpuExtendResult; 2],
    pub n_left: i32,
    pub n_right: i32,
    /// `sc0`: the last left call's score, or `h0` without a left.
    pub right_h0: i32,
    pub pad: i32,
}

const _: () = assert!(std::mem::size_of::<GpuFusedParams>() == 104);
const _: () = assert!(std::mem::size_of::<GpuFusedJob>() == 56);
const _: () = assert!(std::mem::size_of::<GpuFusedResult>() == 112);

/// One seed's fused job: the left query/target (already reversed) and the
/// right query/target. An empty query means no call on that side.
#[derive(Debug, Clone, Copy)]
pub struct GpuFusedInput<'a> {
    pub lq: &'a [u8],
    pub lt: &'a [u8],
    pub rq: &'a [u8],
    pub rt: &'a [u8],
    pub h0: i32,
}

#[derive(Clone, Copy, Default)]
struct Size {
    q: usize,
    t: usize,
    n: usize,
}

/// Reusable huge-page-backed buffers and a stream for
/// [`GpuFusedExtend::compute_groups`] / [`GpuFusedExtend::launch`]. Keep one
/// per lane for the life of the alignment: its buffers stay on pages the GPU
/// has already translated.
pub struct GpuFusedExtend {
    queries: Option<HostBuf>,
    targets: Option<HostBuf>,
    jobs: Option<HostBuf>,
    /// Workspace: the launcher's length-binned job order.
    order: Option<HostBuf>,
    out: Option<HostBuf>,
    /// Created on first launch and reused: nonblocking, the greatest priority
    /// unless NOMICS_FUSED_PRIORITY=0 (a lane waits on its own tile).
    stream: Option<Stream>,
    /// Jobs of the last waited launch.
    last: usize,
}

impl GpuFusedExtend {
    fn ensure_stream(&mut self) -> Result<(), GpuError> {
        if self.stream.is_none() {
            let high = std::env::var("NOMICS_FUSED_PRIORITY").map_or(true, |v| v != "0");
            self.stream = Some(Stream::new(high)?);
        }
        Ok(())
    }

    /// Create this engine's stream now and run an empty kernel on it, so its
    /// driver allocations happen before memory fragments (see
    /// `nomics_stream_warm`).
    pub fn warm(&mut self) -> Result<(), GpuError> {
        if !cfg!(any(feature = "cuda", feature = "metal")) {
            return Err(GpuError::Unavailable);
        }
        self.ensure_stream()?;
        self.stream.as_ref().unwrap().warm()
    }
}

/// A [`GpuFusedExtend`] launch in flight. It owns the engine, whose buffers
/// the GPU reads and writes, until [`Self::wait`] hands it back; dropping it
/// unwaited waits for the launch before freeing them.
pub struct GpuFusedTicket {
    // Declared first: dropped (waited) before the engine's buffers.
    pending: Pending,
    engine: GpuFusedExtend,
    n: usize,
    prep: Duration,
}

impl GpuFusedTicket {
    /// Wait for the launch; the engine's [`GpuFusedExtend::results`] then
    /// hold its results (on success).
    pub fn wait(self) -> (GpuFusedExtend, Result<KswTimings, GpuError>) {
        let GpuFusedTicket {
            pending,
            mut engine,
            n,
            prep,
        } = self;
        match pending.wait() {
            Ok(kernel) => {
                engine.last = n;
                let timings = KswTimings {
                    prep,
                    kernel,
                    post: Duration::ZERO,
                };
                (engine, Ok(timings))
            }
            Err(e) => {
                crate::diag_gpu_fatal("fused extension wait", &e, &engine.diag_buffers());
                (engine, Err(e))
            }
        }
    }
}

impl GpuFusedExtend {
    pub fn new() -> Self {
        Self {
            queries: None,
            targets: None,
            jobs: None,
            order: None,
            out: None,
            stream: None,
            last: 0,
        }
    }

    /// The results of the last waited launch, in group then iterator order.
    pub fn results(&self) -> &[GpuFusedResult] {
        typed(self.out.as_ref(), self.last)
    }

    /// This engine's buffers by field, for diagnostic reports.
    fn diag_buffers(&self) -> Vec<(&'static str, &HostBuf)> {
        crate::extend::named([
            ("fused queries", &self.queries),
            ("fused targets", &self.targets),
            ("fused jobs", &self.jobs),
            ("fused order", &self.order),
            ("fused results", &self.out),
        ])
    }

    /// Bytes of host buffers held (diagnostics).
    pub fn bytes(&self) -> usize {
        [&self.queries, &self.targets, &self.jobs, &self.order, &self.out]
            .iter()
            .filter_map(|b| b.as_ref().map(HostBuf::len))
            .sum()
    }

    fn pack<'a, I>(&mut self, groups: &[I], threads: usize) -> Result<usize, GpuError>
    where
        I: Iterator<Item = GpuFusedInput<'a>> + Clone + Send + Sync,
    {
        let mut sizes = vec![Size::default(); groups.len()];
        parallel(
            groups.iter().zip(&mut sizes).collect(),
            threads,
            |(group, size)| {
                for input in group.clone() {
                    assert!(input.lq.len() <= EXTEND_MAX_QUERY);
                    assert!(input.rq.len() <= EXTEND_MAX_QUERY);
                    assert!(!input.lq.is_empty() || !input.rq.is_empty());
                    assert!(input.lt.len() <= u32::MAX as usize);
                    assert!(input.rt.len() <= u32::MAX as usize);
                    assert!(input.h0 > 0, "BWA ksw_extend2 asserts h0 > 0");
                    size.q += input.lq.len() + input.rq.len();
                    size.t += input.lt.len() + input.rt.len();
                    size.n += 1;
                }
            },
        );
        let mut total = Size::default();
        for size in &sizes {
            total.q += size.q;
            total.t += size.t;
            total.n += size.n;
        }
        assert!(total.n <= u32::MAX as usize);
        ensure(&mut self.queries, total.q, "GPU fused queries")?;
        ensure(&mut self.targets, total.t, "GPU fused targets")?;
        ensure(
            &mut self.jobs,
            total.n * std::mem::size_of::<GpuFusedJob>(),
            "GPU fused jobs",
        )?;
        ensure(
            &mut self.order,
            total.n * std::mem::size_of::<u32>(),
            "GPU fused order",
        )?;
        ensure(
            &mut self.out,
            total.n * std::mem::size_of::<GpuFusedResult>(),
            "GPU fused results",
        )?;
        let mut queries = self.queries.as_mut().unwrap().as_mut_slice();
        let mut targets = self.targets.as_mut().unwrap().as_mut_slice();
        let mut jobs = typed_mut::<GpuFusedJob>(self.jobs.as_mut().unwrap(), total.n);
        let mut work = Vec::with_capacity(groups.len());
        let (mut qo, mut to) = (0, 0);
        for (group, size) in groups.iter().zip(sizes) {
            let (q, rest) = queries.split_at_mut(size.q);
            queries = rest;
            let (t, rest) = targets.split_at_mut(size.t);
            targets = rest;
            let (j, rest) = jobs.split_at_mut(size.n);
            jobs = rest;
            work.push((group.clone(), q, t, j, qo, to));
            qo += size.q;
            to += size.t;
        }
        parallel(work, threads, |(group, q, t, jobs, qo, to)| {
            let (mut qi, mut ti) = (0, 0);
            for (raw, input) in jobs.iter_mut().zip(group) {
                let put = |arena: &mut [u8], at: &mut usize, base: usize, s: &[u8]| {
                    arena[*at..*at + s.len()].copy_from_slice(s);
                    *at += s.len();
                    (base + *at - s.len()) as u64
                };
                let lq_off = put(q, &mut qi, qo, input.lq);
                let rq_off = put(q, &mut qi, qo, input.rq);
                let lt_off = put(t, &mut ti, to, input.lt);
                let rt_off = put(t, &mut ti, to, input.rt);
                *raw = GpuFusedJob {
                    lq_off,
                    lt_off,
                    rq_off,
                    rt_off,
                    lq_len: input.lq.len() as u32,
                    lt_len: input.lt.len() as u32,
                    rq_len: input.rq.len() as u32,
                    rt_len: input.rt.len() as u32,
                    h0: input.h0,
                    pad: 0,
                };
            }
        });
        Ok(total.n)
    }

    /// Pack independent groups in parallel and run one launch over all of
    /// them, then wait. The returned results borrow this engine's huge-page
    /// buffer, in group then iterator order; consume them before reusing the
    /// engine. Separate engines may execute concurrently.
    pub fn compute_groups<'a, I>(
        &mut self,
        params: &GpuFusedParams,
        groups: &[I],
        threads: usize,
    ) -> Result<(&[GpuFusedResult], KswTimings), GpuError>
    where
        I: Iterator<Item = GpuFusedInput<'a>> + Clone + Send + Sync,
    {
        let ticket = match std::mem::take(self).launch(params, groups, threads) {
            Ok(ticket) => ticket,
            Err((engine, e)) => {
                *self = engine;
                return Err(e);
            }
        };
        let (engine, timings) = ticket.wait();
        *self = engine;
        Ok((self.results(), timings?))
    }

    /// [`Self::compute_groups`] without the wait: pack, then enqueue one
    /// launch on this engine's stream. The ticket owns the engine until it is
    /// waited; on error the engine comes back with nothing in flight.
    // The engine itself comes back on failure.
    #[allow(clippy::result_large_err)]
    pub fn launch<'a, I>(
        mut self,
        params: &GpuFusedParams,
        groups: &[I],
        threads: usize,
    ) -> Result<GpuFusedTicket, (Self, GpuError)>
    where
        I: Iterator<Item = GpuFusedInput<'a>> + Clone + Send + Sync,
    {
        match self.launch_in_place(params, groups, threads) {
            Ok((pending, n, prep)) => Ok(GpuFusedTicket {
                pending,
                engine: self,
                n,
                prep,
            }),
            Err(e) => Err((self, e)),
        }
    }

    /// The launch; the caller moves `self` into the ticket with `Pending`.
    fn launch_in_place<'a, I>(
        &mut self,
        params: &GpuFusedParams,
        groups: &[I],
        threads: usize,
    ) -> Result<(Pending, usize, Duration), GpuError>
    where
        I: Iterator<Item = GpuFusedInput<'a>> + Clone + Send + Sync,
    {
        // Preserve the stub's cheap failure, without allocating buffers.
        if !cfg!(any(feature = "cuda", feature = "metal")) {
            return Err(GpuError::Unavailable);
        }
        assert!((1..=8).contains(&params.m));
        assert!((1..=2).contains(&params.max_band_try));
        assert!(params.w > 0);
        self.last = 0;
        self.ensure_stream()?;
        let start = Instant::now();
        let n = self.pack(groups, threads)?;
        let prep = start.elapsed();
        let pending = crate::backend::fused_extend_launch(
            self.stream.as_ref().unwrap(),
            self.queries.as_ref().unwrap().as_slice(),
            self.targets.as_ref().unwrap().as_slice(),
            params,
            typed_mut::<GpuFusedJob>(self.jobs.as_mut().unwrap(), n),
            typed_mut::<u32>(self.order.as_mut().unwrap(), n),
            typed_mut::<GpuFusedResult>(self.out.as_mut().unwrap(), n),
        )?;
        Ok((pending, n, prep))
    }
}

impl Default for GpuFusedExtend {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BASES: [u8; 1024] = {
        let mut bases = [0; 1024];
        let mut i = 0;
        while i < bases.len() {
            bases[i] = ((i ^ (i >> 3) ^ (i >> 5)) & 3) as u8;
            i += 1;
        }
        bases
    };

    // Covers empty left, empty right, empty targets and 256-base queries.
    fn input(i: usize) -> GpuFusedInput<'static> {
        let lq = if i.is_multiple_of(5) { 0 } else { 1 + i % 256 };
        let rq = if i.is_multiple_of(7) && lq > 0 {
            0
        } else {
            1 + (3 * i) % 256
        };
        let lt = if lq == 0 { 0 } else { (5 * i) % 300 };
        let rt = if rq == 0 { 0 } else { (7 * i) % 300 };
        GpuFusedInput {
            lq: &BASES[i % 200..i % 200 + lq],
            lt: &BASES[i % 300..i % 300 + lt],
            rq: &BASES[(2 * i) % 200..(2 * i) % 200 + rq],
            rt: &BASES[(3 * i) % 300..(3 * i) % 300 + rt],
            h0: 1 + i as i32 % 90,
        }
    }

    #[test]
    fn grouped_pack_preserves_all_inputs_and_reuses_buffers() {
        let mut engine = GpuFusedExtend::new();
        for (lengths, threads) in [
            (vec![0, 1, 513, 0, 43, 256, 0], 4),
            (vec![2, 0, 1], 0),
            (vec![0, 0], 3),
            (vec![33, 3000, 0, 19], 7),
            (vec![7, 18], 1),
        ] {
            let mut at = 0;
            let groups: Vec<_> = lengths
                .into_iter()
                .map(|n| {
                    let start = at;
                    at += n;
                    (start..at).map(input)
                })
                .collect();
            let n = engine.pack(&groups, threads).unwrap();
            assert_eq!(n, at);
            let queries = engine.queries.as_ref().unwrap().as_slice();
            let targets = engine.targets.as_ref().unwrap().as_slice();
            let jobs = typed_mut::<GpuFusedJob>(engine.jobs.as_mut().unwrap(), n);
            let (mut qo, mut to) = (0, 0);
            for (i, job) in jobs.iter().enumerate() {
                let x = input(i);
                let want = GpuFusedJob {
                    lq_off: qo as u64,
                    lt_off: to as u64,
                    rq_off: (qo + x.lq.len()) as u64,
                    rt_off: (to + x.lt.len()) as u64,
                    lq_len: x.lq.len() as u32,
                    lt_len: x.lt.len() as u32,
                    rq_len: x.rq.len() as u32,
                    rt_len: x.rt.len() as u32,
                    h0: x.h0,
                    pad: 0,
                };
                assert_eq!(*job, want, "job {i}");
                let q = |o: u64, l: u32| &queries[o as usize..o as usize + l as usize];
                let t = |o: u64, l: u32| &targets[o as usize..o as usize + l as usize];
                assert_eq!(q(job.lq_off, job.lq_len), x.lq);
                assert_eq!(q(job.rq_off, job.rq_len), x.rq);
                assert_eq!(t(job.lt_off, job.lt_len), x.lt);
                assert_eq!(t(job.rt_off, job.rt_len), x.rt);
                qo += x.lq.len() + x.rq.len();
                to += x.lt.len() + x.rt.len();
            }
        }
    }

    #[test]
    fn abi_layout_matches_header() {
        use std::mem::offset_of;
        assert_eq!(offset_of!(GpuFusedParams, mat), 40);
        assert_eq!(offset_of!(GpuFusedJob, lq_len), 32);
        assert_eq!(offset_of!(GpuFusedJob, rq_len), 40);
        assert_eq!(offset_of!(GpuFusedJob, h0), 48);
        assert_eq!(offset_of!(GpuFusedResult, right), 48);
        assert_eq!(offset_of!(GpuFusedResult, n_left), 96);
        assert_eq!(offset_of!(GpuFusedResult, right_h0), 104);
    }

    #[test]
    #[cfg(not(any(feature = "cuda", feature = "metal")))]
    fn stub_fails_without_retaining_buffers() {
        let params = GpuFusedParams {
            m: 5,
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            w: 100,
            max_band_try: 2,
            end_bonus_left: 5,
            end_bonus_right: 5,
            zdrop: 100,
            mat: [0; 64],
        };
        let mut engine = GpuFusedExtend::new();
        assert!(matches!(
            engine.compute_groups(&params, &[(1..5).map(input)], 2),
            Err(GpuError::Unavailable)
        ));
        assert!(engine.queries.is_none() && engine.out.is_none());
        // A refused launch hands the engine back.
        let Err((engine, e)) = engine.launch(&params, &[(1..5).map(input)], 2) else {
            panic!("the stub launched");
        };
        assert_eq!(e, GpuError::Unavailable);
        assert!(engine.queries.is_none() && engine.results().is_empty());
    }

    #[test]
    #[cfg(any(feature = "cuda", feature = "metal"))]
    fn tickets_match_synchronous_results() {
        let mut mat = [0i8; 64];
        for (k, m) in mat[..25].iter_mut().enumerate() {
            *m = if k / 5 == 4 || k % 5 == 4 { -1 } else if k / 5 == k % 5 { 1 } else { -4 };
        }
        let params = GpuFusedParams {
            m: 5,
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            w: 100,
            max_band_try: 2,
            end_bonus_left: 5,
            end_bonus_right: 5,
            zdrop: 100,
            mat,
        };
        let groups = [(1..3000).map(input), (3000..3001).map(input)];
        let mut serial = GpuFusedExtend::new();
        let want = serial.compute_groups(&params, &groups, 4).unwrap().0.to_vec();
        assert_eq!(want.len(), 3000);
        let launch = |engine: GpuFusedExtend| {
            engine
                .launch(&params, &groups, 4)
                .map_err(|(_, e)| e)
                .unwrap()
        };
        // Two engines in flight at once, waited out of order.
        let a = launch(GpuFusedExtend::new());
        let b = launch(GpuFusedExtend::new());
        let (b, rb) = b.wait();
        rb.unwrap();
        assert_eq!(b.results(), &want[..]);
        let (a, ra) = a.wait();
        ra.unwrap();
        assert_eq!(a.results(), &want[..]);
        // A dropped ticket waits for its launch; the engine is then free.
        drop(launch(a));
        let (c, rc) = launch(b).wait();
        rc.unwrap();
        assert_eq!(c.results(), &want[..]);
    }
}
