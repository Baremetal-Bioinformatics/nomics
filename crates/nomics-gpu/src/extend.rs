//! Host packing shared by the CUDA path and Mac tests. No alignment arithmetic
//! lives here: groups keep their input order and all call parameters verbatim.
use crate::backend::{Pending, Stream};
use crate::{
    EXTEND_MAX_QUERY, GpuError, GpuExtendInput, GpuExtendResult, GpuKswExtend2, KswTimings,
};
use nomics_mem::HostBuf;
use std::time::{Duration, Instant};

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Job {
    pub query_offset: u64,
    pub target_offset: u64,
    pub mat_offset: u64,
    pub query_len: u32,
    pub target_len: u32,
    pub m: i32,
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub w: i32,
    pub end_bonus: i32,
    pub zdrop: i32,
    pub h0: i32,
}

#[derive(Clone, Copy, Default)]
struct Size {
    q: usize,
    t: usize,
    m: usize,
    n: usize,
}

/// Packing work on the process-wide host pool ([`crate::pool`]).
pub(crate) fn parallel<T: Send>(work: Vec<T>, threads: usize, f: impl Fn(T) + Sync) {
    crate::pool::for_each(threads, work, f);
}

pub(crate) fn ensure(
    buffer: &mut Option<HostBuf>,
    bytes: usize,
    label: &str,
) -> Result<(), GpuError> {
    if buffer.as_ref().is_none_or(|b| b.len() < bytes) {
        // The need plus 1/4, never above the next power of two (the old
        // sizing): pooled engines keep their largest buffers for the run.
        let need = bytes.max(1);
        let capacity = need.saturating_add(need / 4).min(need.next_power_of_two());
        let mut buf = HostBuf::new(capacity)
            .map_err(|e| GpuError::HostMemory(e.to_string()))?;
        // Engines are pooled, so this runs only while buffers grow.
        buf.touch_huge(label);
        *buffer = Some(buf);
    }
    Ok(())
}

// Used only for Job, GpuExtendResult and the fused job/result types: repr(C), all bit patterns valid, no
// destructors. HostBuf is zero-filled and 2 MiB aligned. The mutable borrow
// excludes every other CPU access while it lives; launched GPU work is
// excluded by ownership (the engine owning the buffer sits in the ticket).
pub(crate) fn typed_mut<T: Copy>(buffer: &mut HostBuf, n: usize) -> &mut [T] {
    assert!(std::mem::align_of::<T>() <= nomics_mem::HUGE_PAGE);
    assert!(std::mem::size_of::<T>().checked_mul(n).unwrap() <= buffer.len());
    unsafe { std::slice::from_raw_parts_mut(buffer.as_mut_ptr().cast(), n) }
}

/// [`typed_mut`] for reading results after a launch was waited.
pub(crate) fn typed<T: Copy>(buffer: Option<&HostBuf>, n: usize) -> &[T] {
    let Some(buffer) = buffer else {
        assert_eq!(n, 0);
        return &[];
    };
    assert!(std::mem::align_of::<T>() <= nomics_mem::HUGE_PAGE);
    assert!(std::mem::size_of::<T>().checked_mul(n).unwrap() <= buffer.len());
    unsafe { std::slice::from_raw_parts(buffer.as_ptr().cast(), n) }
}

/// The present buffers of `fields`, by name (diagnostic reports).
pub(crate) fn named<'a, const N: usize>(
    fields: [(&'static str, &'a Option<HostBuf>); N],
) -> Vec<(&'static str, &'a HostBuf)> {
    fields
        .into_iter()
        .filter_map(|(name, b)| b.as_ref().map(|b| (name, b)))
        .collect()
}

pub(super) struct KswExtend2 {
    queries: Option<HostBuf>,
    targets: Option<HostBuf>,
    matrices: Option<HostBuf>,
    jobs: Option<HostBuf>,
    out: Option<HostBuf>,
    /// Created on first launch and reused: nonblocking, least priority.
    stream: Option<Stream>,
    /// Jobs of the last waited launch.
    pub(crate) last: usize,
}

impl KswExtend2 {
    pub fn new() -> Self {
        Self {
            queries: None,
            targets: None,
            matrices: None,
            jobs: None,
            out: None,
            stream: None,
            last: 0,
        }
    }

    fn pack<'a, I>(&mut self, groups: &[I], threads: usize) -> Result<usize, GpuError>
    where
        I: Iterator<Item = GpuExtendInput<'a>> + Clone + Send + Sync,
    {
        let mut sizes = vec![Size::default(); groups.len()];
        parallel(
            groups.iter().zip(&mut sizes).collect(),
            threads,
            |(group, size)| {
                for input in group.clone() {
                    assert!(!input.query.is_empty() && input.query.len() <= EXTEND_MAX_QUERY);
                    assert!(input.target.len() <= i32::MAX as usize);
                    assert!((1..=8).contains(&input.m));
                    assert_eq!(input.mat.len(), input.m as usize * input.m as usize);
                    assert!(input.h0 > 0, "BWA ksw_extend2 asserts h0 > 0");
                    size.q += input.query.len();
                    size.t += input.target.len();
                    size.m += input.mat.len();
                    size.n += 1;
                }
            },
        );
        let mut total = Size::default();
        for size in &sizes {
            total.q += size.q;
            total.t += size.t;
            total.m += size.m;
            total.n += size.n;
        }
        assert!(total.n <= u32::MAX as usize);
        ensure(&mut self.queries, total.q, "GPU extend queries")?;
        ensure(&mut self.targets, total.t, "GPU extend targets")?;
        ensure(&mut self.matrices, total.m, "GPU extend matrices")?;
        ensure(
            &mut self.jobs,
            total.n * std::mem::size_of::<Job>(),
            "GPU extend jobs",
        )?;
        ensure(
            &mut self.out,
            total.n * std::mem::size_of::<GpuExtendResult>(),
            "GPU extend results",
        )?;
        let mut queries = self.queries.as_mut().unwrap().as_mut_slice();
        let mut targets = self.targets.as_mut().unwrap().as_mut_slice();
        let mut matrices = self.matrices.as_mut().unwrap().as_mut_slice();
        let mut jobs = typed_mut::<Job>(self.jobs.as_mut().unwrap(), total.n);
        let mut work = Vec::with_capacity(groups.len());
        let (mut qo, mut to, mut mo) = (0, 0, 0);
        for (group, size) in groups.iter().zip(sizes) {
            let (q, rest) = queries.split_at_mut(size.q);
            queries = rest;
            let (t, rest) = targets.split_at_mut(size.t);
            targets = rest;
            let (m, rest) = matrices.split_at_mut(size.m);
            matrices = rest;
            let (j, rest) = jobs.split_at_mut(size.n);
            jobs = rest;
            work.push((group.clone(), q, t, m, j, qo, to, mo));
            qo += size.q;
            to += size.t;
            mo += size.m;
        }
        parallel(work, threads, |(group, q, t, m, jobs, qo, to, mo)| {
            let (mut qi, mut ti, mut mi) = (0, 0, 0);
            for (raw, input) in jobs.iter_mut().zip(group) {
                q[qi..qi + input.query.len()].copy_from_slice(input.query);
                t[ti..ti + input.target.len()].copy_from_slice(input.target);
                for (slot, &v) in m[mi..mi + input.mat.len()].iter_mut().zip(input.mat) {
                    *slot = v as u8;
                }
                *raw = Job {
                    query_offset: (qo + qi) as u64,
                    target_offset: (to + ti) as u64,
                    mat_offset: (mo + mi) as u64,
                    query_len: input.query.len() as u32,
                    target_len: input.target.len() as u32,
                    m: input.m,
                    o_del: input.o_del,
                    e_del: input.e_del,
                    o_ins: input.o_ins,
                    e_ins: input.e_ins,
                    w: input.w,
                    end_bonus: input.end_bonus,
                    zdrop: input.zdrop,
                    h0: input.h0,
                };
                qi += input.query.len();
                ti += input.target.len();
                mi += input.mat.len();
            }
        });
        Ok(total.n)
    }

    /// Pack `groups` and enqueue one launch on this engine's stream. The
    /// ticket covers this engine's buffers: keep the engine, untouched, with
    /// it until it is waited. Returns the ticket, the job count and the
    /// packing time.
    pub(crate) fn launch<'a, I>(
        &mut self,
        groups: &[I],
        threads: usize,
    ) -> Result<(Pending, usize, Duration), GpuError>
    where
        I: Iterator<Item = GpuExtendInput<'a>> + Clone + Send + Sync,
    {
        // Preserve the stub's cheap failure, without allocating GPU buffers.
        if !cfg!(feature = "cuda") {
            return Err(GpuError::Unavailable);
        }
        self.last = 0;
        if self.stream.is_none() {
            self.stream = Some(Stream::new(false)?);
        }
        let start = Instant::now();
        let n = self.pack(groups, threads)?;
        let prep = start.elapsed();
        let pending = crate::backend::ksw_extend2_launch(
            self.stream.as_ref().unwrap(),
            self.queries.as_ref().unwrap().as_slice(),
            self.targets.as_ref().unwrap().as_slice(),
            self.matrices.as_ref().unwrap().as_slice(),
            typed_mut::<Job>(self.jobs.as_mut().unwrap(), n),
            typed_mut::<GpuExtendResult>(self.out.as_mut().unwrap(), n),
        )?;
        Ok((pending, n, prep))
    }

    /// This engine's buffers by field, for diagnostic reports.
    pub(crate) fn diag_buffers(&self) -> Vec<(&'static str, &HostBuf)> {
        named([
            ("extend queries", &self.queries),
            ("extend targets", &self.targets),
            ("extend matrices", &self.matrices),
            ("extend jobs", &self.jobs),
            ("extend results", &self.out),
        ])
    }

    /// The results of the last waited launch, in group then iterator order.
    pub(crate) fn results(&self) -> &[GpuExtendResult] {
        typed(self.out.as_ref(), self.last)
    }

    /// Synchronous [`Self::launch`] and wait.
    pub fn compute<'a, I>(
        &mut self,
        groups: &[I],
        threads: usize,
    ) -> Result<(&[GpuExtendResult], KswTimings), GpuError>
    where
        I: Iterator<Item = GpuExtendInput<'a>> + Clone + Send + Sync,
    {
        let (pending, n, prep) = self.launch(groups, threads)?;
        let kernel = pending.wait().inspect_err(|e| {
            crate::diag_gpu_fatal("two-pass extension wait", e, &self.diag_buffers())
        })?;
        self.last = n;
        Ok((
            self.results(),
            KswTimings {
                prep,
                kernel,
                post: Duration::ZERO,
            },
        ))
    }
}

/// A [`GpuKswExtend2`] launch in flight. It owns the engine, whose buffers
/// the GPU reads and writes, until [`Self::wait`] hands it back; dropping it
/// unwaited waits for the launch before freeing them.
pub struct GpuExtendTicket {
    // Declared first: dropped (waited) before the engine's buffers.
    pending: Pending,
    engine: GpuKswExtend2,
    n: usize,
    prep: Duration,
}

impl GpuExtendTicket {
    pub(crate) fn launched(
        engine: GpuKswExtend2,
        (pending, n, prep): (Pending, usize, Duration),
    ) -> Self {
        Self {
            pending,
            engine,
            n,
            prep,
        }
    }

    /// Wait for the launch; the engine's [`GpuKswExtend2::results`] then
    /// hold its results (on success).
    pub fn wait(self) -> (GpuKswExtend2, Result<KswTimings, GpuError>) {
        let GpuExtendTicket {
            pending,
            mut engine,
            n,
            prep,
        } = self;
        match pending.wait() {
            Ok(kernel) => {
                engine.inner.last = n;
                let timings = KswTimings {
                    prep,
                    kernel,
                    post: Duration::ZERO,
                };
                (engine, Ok(timings))
            }
            Err(e) => {
                crate::diag_gpu_fatal("two-pass extension wait", &e, &engine.inner.diag_buffers());
                (engine, Err(e))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(i: usize) -> GpuExtendInput<'static> {
        const BASES: [u8; 512] = {
            let mut bases = [0; 512];
            let mut i = 0;
            while i < bases.len() {
                bases[i] = ((i ^ (i >> 3)) & 1) as u8;
                i += 1;
            }
            bases
        };
        const MAT: [i8; 4] = [2, -3, -5, 7];
        GpuExtendInput {
            query: &BASES[i % 257..i % 257 + 1 + i % 256],
            target: &BASES[(3 * i) % 257..(3 * i) % 257 + i % 257], // Includes zero-length targets.
            m: 2,
            mat: &MAT,
            o_del: 6,
            e_del: 2,
            o_ins: 8,
            e_ins: 1,
            w: 5 + i as i32 % 100,
            end_bonus: 3,
            zdrop: 20,
            h0: 1 + i as i32,
        }
    }

    #[test]
    fn grouped_pack_preserves_all_inputs_and_reuses_buffers() {
        assert_eq!(std::mem::size_of::<Job>(), 72);
        assert_eq!(std::mem::size_of::<GpuExtendResult>(), 24);
        let mut engine = KswExtend2::new();
        for (lengths, threads) in [
            (vec![0, 1, 513, 0, 43, 256, 0], 4),
            (vec![2, 0, 1], 0),
            (vec![0, 0], 3),
            (vec![33, 5000, 0, 19], 7),
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
            let matrices = engine.matrices.as_ref().unwrap().as_slice();
            let jobs = typed_mut::<Job>(engine.jobs.as_mut().unwrap(), n);
            let (mut qo, mut to, mut mo) = (0, 0, 0);
            for (i, job) in jobs.iter().enumerate() {
                let expected = input(i);
                let want = Job {
                    query_offset: qo as u64,
                    target_offset: to as u64,
                    mat_offset: mo as u64,
                    query_len: expected.query.len() as u32,
                    target_len: expected.target.len() as u32,
                    m: expected.m,
                    o_del: expected.o_del,
                    e_del: expected.e_del,
                    o_ins: expected.o_ins,
                    e_ins: expected.e_ins,
                    w: expected.w,
                    end_bonus: expected.end_bonus,
                    zdrop: expected.zdrop,
                    h0: expected.h0,
                };
                assert_eq!(*job, want, "job {i}");
                assert_eq!(&queries[qo..qo + expected.query.len()], expected.query);
                assert_eq!(&targets[to..to + expected.target.len()], expected.target);
                assert_eq!(matrices[mo..mo + 4], [2, 253, 251, 7]);
                qo += expected.query.len();
                to += expected.target.len();
                mo += expected.mat.len();
            }
        }
    }

    #[test]
    #[cfg(not(feature = "cuda"))]
    fn stub_fails_without_retaining_buffers() {
        let mut engine = KswExtend2::new();
        assert!(matches!(
            engine.compute(&[(0..4).map(input)], 2),
            Err(GpuError::Unavailable)
        ));
        assert!(engine.queries.is_none() && engine.out.is_none());
    }

    #[test]
    #[cfg(feature = "cuda")]
    fn concurrent_engines_match_serial_results() {
        let groups = [
            (0..4000).map(input),
            (0..0).map(input),
            (4000..5000).map(input),
        ];
        let mut serial = KswExtend2::new();
        let expected = serial.compute(&groups, 1).unwrap().0.to_vec();
        std::thread::scope(|s| {
            for _ in 0..4 {
                s.spawn(|| {
                    let mut engine = KswExtend2::new();
                    for _ in 0..3 {
                        let (out, _) = engine.compute(&groups, 3).unwrap();
                        assert_eq!(out, expected);
                    }
                });
            }
        });
    }
}
