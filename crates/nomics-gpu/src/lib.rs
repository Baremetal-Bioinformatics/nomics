//! GPU kernels behind a C ABI (`shim/nomics_gpu.h`).
//!
//! Built with `--features cuda` on the Spark. `--features metal` (Apple
//! Silicon, `shim/metal/`) implements the same ABI behind the same Rust
//! backend; kernels it has not ported (or cannot: PairHMM is FP64) return
//! [`GpuError::Unavailable`], so their callers keep the CPU path. So far:
//! the device query, the in-place smoke sum and mate rescue (`ksw_align2`).
//! Without either, every entry point returns [`GpuError::Unavailable`], so the
//! workspace builds and tests anywhere and every stage keeps a CPU path (the
//! control arm that GPU speedups are measured against).
#![allow(unsafe_code)]

#[cfg(all(feature = "cuda", feature = "metal"))]
compile_error!("features `cuda` and `metal` are exclusive");

mod extend;
mod fused;
pub mod pool;

pub use extend::GpuExtendTicket;
pub use fused::{
    GpuFusedExtend, GpuFusedInput, GpuFusedJob, GpuFusedParams, GpuFusedResult, GpuFusedTicket,
};

use std::fmt;
use std::time::Duration;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GpuError {
    /// Built without a GPU backend, or the backend lacks this kernel.
    Unavailable,
    /// A CUDA runtime call failed with this `cudaError_t`.
    Cuda(i32),
    /// The Metal shim failed with this status (`shim/metal/nomics_metal.m`).
    Metal(i32),
    /// A huge-page-friendly host buffer could not be allocated.
    HostMemory(String),
}

impl fmt::Display for GpuError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GpuError::Unavailable => f.write_str("not available in this build's GPU backend"),
            GpuError::Cuda(code) => write!(f, "CUDA error {code}"),
            GpuError::Metal(code) => write!(f, "Metal shim error {code}"),
            GpuError::HostMemory(message) => write!(f, "host memory allocation failed: {message}"),
        }
    }
}

impl std::error::Error for GpuError {}

/// With `NOMICS_DIAG` (docs/DIAG.md): a failed GPU call is fatal, never a
/// silent CPU fallback. A fault on a guard page leaves a sticky CUDA error,
/// so falling back would run the rest of the job on the CPU and still pass
/// every gate. Prints `what`, the error, `buffers` (the failed launch's own,
/// by field) with their canary state and every live HostBuf, then aborts.
/// A no-op without `NOMICS_DIAG` and for [`GpuError::Unavailable`] (a kernel
/// the backend does not have, which is not a failure).
pub fn diag_gpu_error(what: &str, err: &GpuError) {
    diag_gpu_fatal(what, err, &[]);
}

pub(crate) fn diag_gpu_fatal(what: &str, err: &GpuError, buffers: &[(&str, &nomics_mem::HostBuf)]) {
    if !nomics_mem::diag() || *err == GpuError::Unavailable {
        return;
    }
    let hint = match err {
        GpuError::Cuda(700) => {
            " (cudaErrorIllegalAddress: a kernel touched an unmapped or PROT_NONE address)"
        }
        GpuError::Cuda(716) => " (cudaErrorMisalignedAddress)",
        GpuError::Cuda(719) => " (cudaErrorLaunchFailure)",
        GpuError::Cuda(_) => " (cudaError_t; see cuda_runtime_api.h)",
        _ => "",
    };
    let mut report = format!("GPU error in {what}: {err}{hint}\n");
    if !buffers.is_empty() {
        report.push_str("buffers of the failed launch:\n");
        for (name, buf) in buffers {
            let canaries = match buf.check_canaries(name) {
                Ok(()) => "canaries intact".to_owned(),
                Err(d) => d.to_string(),
            };
            report.push_str(&format!(
                "  {name}: data {:#x}..{:#x} (len {}), accessible to {:#x}; {canaries}\n",
                buf.as_ptr() as usize,
                buf.as_ptr() as usize + buf.len(),
                buf.len(),
                buf.as_ptr() as usize + buf.accessible_len()
            ));
        }
    }
    report.push_str(&nomics_mem::diag::describe_all());
    report.push_str(
        "A GPU fault on a guard page should also appear in the kernel log (dmesg: an SMMU/ATS or NVRM Xid fault with the faulting address); match that address against the ranges above.",
    );
    nomics_mem::diag::fatal(&report);
}

#[derive(Debug, Clone)]
pub struct DeviceInfo {
    pub name: String,
    pub compute_capability: (i32, i32),
    pub multiprocessors: i32,
    /// The GPU can dereference ordinary host pointers (required by our design).
    pub pageable_memory_access: bool,
    pub total_global_mem: usize,
    /// Largest single buffer the GPU can map (`usize::MAX`: no limit).
    pub max_buffer_bytes: usize,
}

pub fn device_info(device: i32) -> Result<DeviceInfo, GpuError> {
    backend::device_info(device)
}

/// A [`HostBuf`](nomics_mem::HostBuf) registered with the GPU backend for its
/// whole life (`nomics_host_register`): Metal wraps it once instead of on every
/// launch that reads it (~85 ms per launch for the 3.2 GB BWT). CUDA reads host
/// memory in place and builds without a GPU do nothing, so this is free there.
/// Owning the buffer makes the registration impossible to outlive; dropping it
/// unregisters before the memory is unmapped. Its owner keeps it until every
/// launch that touched it was waited, as for any buffer a launch reads.
pub struct GpuHostBuf {
    buf: nomics_mem::HostBuf,
    registered: bool,
}

impl GpuHostBuf {
    /// Register `buf`. If the backend refuses, the buffer stays usable and
    /// unregistered: launches then wrap it per call (same results).
    pub fn register(buf: nomics_mem::HostBuf) -> Self {
        let registered = host_register(&buf);
        Self { buf, registered }
    }

    /// Whether the backend holds a registration for this buffer.
    pub fn is_registered(&self) -> bool {
        self.registered
    }
}

impl std::ops::Deref for GpuHostBuf {
    type Target = nomics_mem::HostBuf;
    fn deref(&self) -> &nomics_mem::HostBuf {
        &self.buf
    }
}

impl std::ops::DerefMut for GpuHostBuf {
    fn deref_mut(&mut self) -> &mut nomics_mem::HostBuf {
        &mut self.buf
    }
}

impl Drop for GpuHostBuf {
    fn drop(&mut self) {
        if self.registered {
            host_unregister(&self.buf);
        }
    }
}

#[cfg(any(feature = "cuda", feature = "metal"))]
fn host_register(buf: &nomics_mem::HostBuf) -> bool {
    // The mapping covers whole 2 MiB pages past `len` (only `len`'s own pages
    // with NOMICS_DIAG's guard pages).
    let mapped = buf.accessible_len();
    unsafe { ffi_common::nomics_host_register(buf.as_ptr().cast(), mapped) == 0 }
}

#[cfg(any(feature = "cuda", feature = "metal"))]
fn host_unregister(buf: &nomics_mem::HostBuf) {
    unsafe { ffi_common::nomics_host_unregister(buf.as_ptr().cast()) };
}

#[cfg(not(any(feature = "cuda", feature = "metal")))]
fn host_register(_buf: &nomics_mem::HostBuf) -> bool {
    false
}

#[cfg(not(any(feature = "cuda", feature = "metal")))]
fn host_unregister(_buf: &nomics_mem::HostBuf) {}

/// Sum `data` on the GPU, reading it in place from host memory. Returns the sum
/// and the kernel time in milliseconds. Smoke test for the unified-memory path.
pub fn sum_u64_in_place(data: &[u64]) -> Result<(u64, f32), GpuError> {
    backend::sum_u64_in_place(data)
}

/// bwa's BWT as the seed kernel reads it (see `shim/nomics_gpu.h`).
#[derive(Debug, Clone, Copy)]
pub struct GpuBwt<'a> {
    pub primary: u64,
    pub l2: [u64; 5],
    pub seq_len: u64,
    pub words: &'a [u32],
}

/// `bwtintv_t` (`nomics_intv`).
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GpuIntv {
    pub x: [u64; 3],
    pub info: u64,
}

/// Reads longer than this are left to the CPU.
pub const SEED_MAX_LEN: usize = 256;
/// `out_n` value for reads the GPU did not finish (too many intervals or too long).
pub const SEED_OVERFLOW: u32 = u32::MAX;
pub const KMER_MAX_K: usize = 13;

/// Exact prefix-endpoint table as consumed by the seed kernel. The intervals
/// are three `u64`s per entry and remain in caller-owned host memory.
#[derive(Debug, Clone, Copy)]
pub struct GpuKmerTable<'a> {
    pub k: u32,
    pub level_offsets: [u64; KMER_MAX_K + 1],
    pub intervals: &'a [u64],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedWork {
    Cap64,
    Cap96,
    Cap128,
    /// Run the 64-entry kernel first, then retry only its overflow reads with
    /// the 128-entry kernel.
    TwoPass64Then128,
}

impl SeedWork {
    #[cfg(any(feature = "cuda", feature = "metal"))]
    fn abi(self) -> u32 {
        match self {
            Self::Cap64 => 64,
            Self::Cap96 => 96,
            Self::Cap128 => 128,
            Self::TwoPass64Then128 => 192,
        }
    }
}

/// Doubles per read row in the PairHMM ABI, matching `nomics_pairhmm_pair`.
pub const PAIRHMM_ROW_PARAMS: usize = 7;
/// Longest read handled by the experimental CUDA PairHMM kernel.
pub const PAIRHMM_MAX_READ_LEN: usize = 160;

/// One `(read, haplotype)` PairHMM job. Base offsets count bytes/bases, while
/// `params_offset` counts seven-double rows.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GpuPairHmmPair {
    pub read_offset: u64,
    pub hap_offset: u64,
    pub params_offset: u64,
    pub read_len: u32,
    pub hap_len: u32,
    /// `2^1020 / hap_len`, precomputed by the exact CPU implementation.
    pub initial: f64,
}

/// Host-side preparation/copy, CUDA-event kernel time, and output copy.
#[derive(Debug, Clone, Copy, Default)]
pub struct PairHmmTimings {
    pub prep: Duration,
    pub kernel: Duration,
    pub post: Duration,
}

/// Reusable huge-page-backed buffers for the experimental PairHMM kernel.
pub struct GpuPairHmm {
    inner: backend::PairHmm,
}

impl GpuPairHmm {
    pub fn new() -> Self {
        Self {
            inner: backend::PairHmm::new(),
        }
    }

    /// Compute raw PairHMM final-row sums. Row parameters must be laid out as
    /// `[m2m, m2i, m2d, i2m, i2i, p_match, p_mismatch]` per read row. The
    /// caller performs `jmath::log10(sum) - log10(initial_condition)` on CPU.
    pub fn compute(
        &mut self,
        read_bases: &[u8],
        hap_bases: &[u8],
        row_params: &[[f64; PAIRHMM_ROW_PARAMS]],
        pairs: &[GpuPairHmmPair],
        out_sum: &mut [f64],
    ) -> Result<PairHmmTimings, GpuError> {
        assert!(out_sum.len() >= pairs.len());
        assert!(pairs.len() <= u32::MAX as usize);
        for pair in pairs {
            let read_end = pair.read_offset as usize + pair.read_len as usize;
            let hap_end = pair.hap_offset as usize + pair.hap_len as usize;
            let params_end = pair.params_offset as usize + pair.read_len as usize;
            assert!(
                pair.read_len > 0 && pair.read_len as usize <= PAIRHMM_MAX_READ_LEN,
                "PairHMM read length {} is outside 1..={PAIRHMM_MAX_READ_LEN}",
                pair.read_len
            );
            assert!(pair.hap_len > 0, "PairHMM haplotype must not be empty");
            assert!(read_end <= read_bases.len());
            assert!(hap_end <= hap_bases.len());
            assert!(params_end <= row_params.len());
        }
        self.inner
            .compute(read_bases, hap_bases, row_params, pairs, out_sum)
    }
}

impl Default for GpuPairHmm {
    fn default() -> Self {
        Self::new()
    }
}

/// One experimental BWA `ksw_align2` call. `size` is the effective query
/// profile width from the dump (1 for saturating bytes, 2 for signed i16).
#[derive(Debug, Clone, Copy)]
pub struct GpuKswAlign2Input<'a> {
    pub query: &'a [u8],
    pub target: &'a [u8],
    pub m: i32,
    pub mat: &'a [i8],
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub xtra: i32,
    pub size: u8,
}

/// All fields returned by BWA 0.7.19's `ksw_align2`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GpuKswResult {
    pub score: i32,
    pub te: i32,
    pub qe: i32,
    pub score2: i32,
    pub te2: i32,
    pub tb: i32,
    pub qb: i32,
}

/// One BWA `ksw_extend2` call for [`GpuKswExtend2`].
#[derive(Debug, Clone, Copy)]
pub struct GpuExtendInput<'a> {
    pub query: &'a [u8],
    pub target: &'a [u8],
    pub m: i32,
    pub mat: &'a [i8],
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub w: i32,
    pub end_bonus: i32,
    pub zdrop: i32,
    pub h0: i32,
}

/// All fields returned by BWA 0.7.19's `ksw_extend2`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GpuExtendResult {
    pub score: i32,
    pub qle: i32,
    pub tle: i32,
    pub gtle: i32,
    pub gscore: i32,
    pub max_off: i32,
}

/// Longest query [`GpuKswExtend2`] accepts.
pub const EXTEND_MAX_QUERY: usize = 256;

/// Reusable huge-page-backed buffers for `ksw_extend2` on the GPU (one warp
/// per call, row by row; exact). One launch handles every input.
pub struct GpuKswExtend2 {
    inner: extend::KswExtend2,
}

impl GpuKswExtend2 {
    pub fn new() -> Self {
        Self {
            inner: extend::KswExtend2::new(),
        }
    }

    pub fn compute(
        &mut self,
        inputs: &[GpuExtendInput<'_>],
        out: &mut [GpuExtendResult],
    ) -> Result<KswTimings, GpuError> {
        assert!(out.len() >= inputs.len());
        let (results, mut timings) = self.compute_groups(&[inputs.iter().copied()], 1)?;
        let post = std::time::Instant::now();
        out[..results.len()].copy_from_slice(results);
        timings.post += post.elapsed();
        Ok(timings)
    }

    /// Pack independent groups in parallel without materializing an input
    /// vector, launch and wait. The returned results borrow this engine's
    /// huge-page buffer; consume them before reusing the engine. Group and
    /// iterator order are preserved. Separate engines may execute
    /// concurrently (each has its own stream).
    pub fn compute_groups<'a, I>(
        &mut self,
        groups: &[I],
        threads: usize,
    ) -> Result<(&[GpuExtendResult], KswTimings), GpuError>
    where
        I: Iterator<Item = GpuExtendInput<'a>> + Clone + Send + Sync,
    {
        self.inner.compute(groups, threads)
    }

    /// [`Self::compute_groups`] without the wait: the ticket owns the engine
    /// until it is waited. On error the engine comes back with nothing in
    /// flight.
    // The engine itself comes back on failure.
    #[allow(clippy::result_large_err)]
    pub fn launch<'a, I>(
        mut self,
        groups: &[I],
        threads: usize,
    ) -> Result<GpuExtendTicket, (Self, GpuError)>
    where
        I: Iterator<Item = GpuExtendInput<'a>> + Clone + Send + Sync,
    {
        match self.inner.launch(groups, threads) {
            Ok(launched) => Ok(GpuExtendTicket::launched(self, launched)),
            Err(e) => Err((self, e)),
        }
    }

    /// The results of the last waited launch.
    pub fn results(&self) -> &[GpuExtendResult] {
        self.inner.results()
    }
}

impl Default for GpuKswExtend2 {
    fn default() -> Self {
        Self::new()
    }
}

/// Host packing, CUDA-event kernel time, and result copy for one KSW batch.
#[derive(Debug, Clone, Copy, Default)]
pub struct KswTimings {
    pub prep: Duration,
    pub kernel: Duration,
    pub post: Duration,
}

/// Reusable huge-page-backed buffers for the experimental `ksw_align2`
/// kernel. One launch handles every input in [`Self::compute`].
///
/// Each engine has its own stream (created on first launch) and buffers, so
/// separate engines run concurrently without a lock.
pub struct GpuKswAlign2 {
    inner: backend::KswAlign2,
    /// Inputs of the last waited launch (what [`Self::results`] returns).
    last: usize,
}

/// A [`GpuKswAlign2`] launch in flight. It owns the engine, whose buffers
/// the GPU reads and writes, until [`Self::wait`] hands it back; dropping it
/// unwaited waits for the launch before freeing them.
pub struct GpuKswTicket {
    // Declared first: dropped (waited) before the engine's buffers.
    pending: backend::Pending,
    engine: GpuKswAlign2,
    n: usize,
    prep: Duration,
}

impl GpuKswTicket {
    /// Wait for the launch; the engine's [`GpuKswAlign2::results`] then hold
    /// its results (on success).
    pub fn wait(self) -> (GpuKswAlign2, Result<KswTimings, GpuError>) {
        let GpuKswTicket {
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
                diag_gpu_fatal("mate-rescue wait", &e, &engine.inner.diag_buffers());
                (engine, Err(e))
            }
        }
    }
}

/// Create the GPU seeding stream and run an empty kernel on it now: call
/// before loading the index, so the driver's per-stream context buffers are
/// allocated while memory is still contiguous (2026-10-22: a 30x run whose
/// context allocation failed ran every batch on the CPU).
pub fn warm_seeding() -> Result<(), GpuError> {
    backend::warm_seed_stream()
}

/// Query plus target bases per mate-rescue launch in
/// [`GpuKswAlign2::compute`] (`NOMICS_RESCUE_LAUNCH_MB`, default 2 GiB): a
/// bound on each pooled engine's buffers for outsized batches. A 30x batch
/// rescues ~1 GB of bases (~750k calls) and stays one launch: rescue is on
/// each batch's serial tail, and 256 MiB launches (4 per batch, each waited
/// before the next is packed) made the 30x chunk ~10% slower (204 -> 225 s,
/// 2026-10-23) for ~3 GB less memory at 3 batches in flight.
pub fn rescue_launch_bytes() -> usize {
    static BYTES: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *BYTES.get_or_init(|| {
        std::env::var("NOMICS_RESCUE_LAUNCH_MB")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&mb| mb > 0)
            .map_or(2 << 30, |mb| mb << 20)
    })
}

/// Most calls per mate-rescue launch (bounds the job and result buffers).
pub const RESCUE_LAUNCH_JOBS: usize = 1 << 20;

impl GpuKswAlign2 {
    /// Create this engine's stream now and run an empty kernel on it.
    pub fn warm(&mut self) -> Result<(), GpuError> {
        self.inner.warm()
    }

    pub fn new() -> Self {
        Self {
            inner: backend::KswAlign2::new(),
            last: 0,
        }
    }

    /// The results of the last waited launch, in input order. After
    /// [`Self::compute`], which may split a batch into several launches,
    /// that is only the last slice: use the `out` it fills instead.
    pub fn results(&self) -> &[GpuKswResult] {
        self.inner.results(self.last)
    }

    /// Bytes of host buffers this engine holds (diagnostics).
    pub fn bytes(&self) -> usize {
        self.inner.bytes()
    }

    /// [`Self::compute_bounded`] with [`rescue_launch_bytes`] per launch.
    pub fn compute(
        &mut self,
        inputs: &[GpuKswAlign2Input<'_>],
        out: &mut [GpuKswResult],
    ) -> Result<KswTimings, GpuError> {
        self.compute_bounded(inputs, out, rescue_launch_bytes())
    }

    /// Synchronous [`Self::launch`] + [`GpuKswTicket::wait`] over consecutive
    /// slices of `inputs`, each at most `max_bytes` of query plus target
    /// bases (and [`RESCUE_LAUNCH_JOBS`] jobs; always at least one job), so
    /// the engine's buffers stay bounded whatever the batch size; results
    /// copied to `out` in input order. Every call is independent of the
    /// others in its launch, so the slicing does not change any result.
    pub fn compute_bounded(
        &mut self,
        inputs: &[GpuKswAlign2Input<'_>],
        out: &mut [GpuKswResult],
        max_bytes: usize,
    ) -> Result<KswTimings, GpuError> {
        assert!(out.len() >= inputs.len());
        let mut total = KswTimings::default();
        // [`Self::results`] then covers only the last slice; an empty call
        // leaves none. Callers read `out`.
        self.last = 0;
        let mut start = 0;
        while start < inputs.len() {
            let mut end = start;
            let mut bytes = 0usize;
            while end < inputs.len() && end - start < RESCUE_LAUNCH_JOBS {
                let call = inputs[end].query.len() + inputs[end].target.len();
                if end > start && bytes + call > max_bytes {
                    break;
                }
                bytes += call;
                end += 1;
            }
            let ticket = match std::mem::take(self).launch(&inputs[start..end]) {
                Ok(ticket) => ticket,
                Err((engine, e)) => {
                    *self = engine;
                    return Err(e);
                }
            };
            let (engine, timings) = ticket.wait();
            *self = engine;
            let timings = timings?;
            let post = std::time::Instant::now();
            out[start..end].copy_from_slice(self.results());
            total.prep += timings.prep;
            total.kernel += timings.kernel;
            total.post += post.elapsed();
            start = end;
        }
        Ok(total)
    }

    /// Pack `inputs` and enqueue one launch on this engine's stream. The
    /// ticket owns the engine until it is waited; on error the engine comes
    /// back with nothing in flight.
    // The engine itself comes back on failure.
    #[allow(clippy::result_large_err)]
    pub fn launch(
        mut self,
        inputs: &[GpuKswAlign2Input<'_>],
    ) -> Result<GpuKswTicket, (Self, GpuError)> {
        assert!(inputs.len() <= u32::MAX as usize);
        for input in inputs {
            assert!(!input.query.is_empty(), "BWA ksw_align2 requires a query");
            assert!(input.query.len() <= u32::MAX as usize);
            assert!(input.target.len() <= u32::MAX as usize);
            assert!(input.m > 0);
            assert_eq!(input.mat.len(), input.m as usize * input.m as usize);
            assert!(input.size == 1 || input.size == 2);
            assert!(
                input
                    .query
                    .iter()
                    .chain(input.target)
                    .all(|&base| base < input.m as u8)
            );
        }
        self.last = 0;
        match self.inner.launch(inputs) {
            Ok((pending, prep)) => Ok(GpuKswTicket {
                pending,
                engine: self,
                n: inputs.len(),
                prep,
            }),
            Err(e) => Err((self, e)),
        }
    }
}

impl Default for GpuKswAlign2 {
    fn default() -> Self {
        Self::new()
    }
}

/// `mem_collect_intv` for a batch of reads on the GPU, all buffers in host
/// memory. Read `r` is `seqs[offs[r] as usize..][..lens[r] as usize]` (codes
/// 0..4); its sorted intervals land in `out[r * cap..][..out_n[r]]`, or
/// `out_n[r] == SEED_OVERFLOW`. Returns the kernel time in milliseconds.
pub fn seed_collect(
    bwt: &GpuBwt,
    seqs: &[u8],
    offs: &[u64],
    lens: &[u32],
    cap: u32,
    out: &mut [GpuIntv],
    out_n: &mut [u32],
) -> Result<f32, GpuError> {
    seed_collect_with_options(
        bwt,
        None,
        seqs,
        offs,
        lens,
        cap,
        SeedWork::Cap64,
        out,
        out_n,
    )
}

/// Configurable seed kernel used by the k-mer/cap experiment.
#[allow(clippy::too_many_arguments)]
pub fn seed_collect_with_options(
    bwt: &GpuBwt,
    table: Option<&GpuKmerTable<'_>>,
    seqs: &[u8],
    offs: &[u64],
    lens: &[u32],
    cap: u32,
    work: SeedWork,
    out: &mut [GpuIntv],
    out_n: &mut [u32],
) -> Result<f32, GpuError> {
    let n = lens.len();
    assert_eq!(offs.len(), n);
    assert!(cap > 0);
    assert!(out_n.len() >= n && out.len() >= n * cap as usize);
    assert!(
        offs.iter()
            .zip(lens)
            .all(|(&o, &l)| o as usize + l as usize <= seqs.len())
    );
    assert_eq!(
        bwt.words.as_ptr() as usize % 64,
        0,
        "BWT words must be 64-byte aligned"
    );
    if let Some(table) = table {
        assert!((2..=KMER_MAX_K as u32).contains(&table.k));
        let entries = table.level_offsets[table.k as usize] + (1u64 << (2 * table.k));
        assert!(table.intervals.len() >= entries as usize * 3);
        assert_eq!(
            table.intervals.as_ptr() as usize % 64,
            0,
            "k-mer intervals must be 64-byte aligned"
        );
    }
    backend::seed_collect(bwt, table, seqs, offs, lens, cap, work, out, out_n)
}

/// Huge-page-backed result storage for large seed batches, for up to
/// `capacity()` reads per call. Reuse one across batches: the kernel writes
/// it in place, and fresh pages cost the GPU first-touch translation work
/// (measured on a 30x batch: ~1 s, more than the kernel itself).
pub struct GpuSeedBuffers {
    intervals: GpuHostBuf,
    counts: GpuHostBuf,
    reads: usize,
    cap: u32,
    /// Reads in the last `collect` call.
    used: usize,
}

impl GpuSeedBuffers {
    pub fn new(reads: usize, cap: u32) -> Result<Self, GpuError> {
        if cap == 0 {
            return Err(GpuError::HostMemory("seed result cap is zero".into()));
        }
        let interval_bytes = reads
            .checked_mul(cap as usize)
            .and_then(|n| n.checked_mul(std::mem::size_of::<GpuIntv>()))
            .ok_or_else(|| GpuError::HostMemory("seed result size overflow".into()))?;
        let count_bytes = reads
            .checked_mul(std::mem::size_of::<u32>())
            .ok_or_else(|| GpuError::HostMemory("seed count size overflow".into()))?;
        Ok(Self {
            intervals: GpuHostBuf::register(
                nomics_mem::HostBuf::new(interval_bytes)
                    .map_err(|e| GpuError::HostMemory(e.to_string()))?,
            ),
            counts: GpuHostBuf::register(
                nomics_mem::HostBuf::new(count_bytes)
                    .map_err(|e| GpuError::HostMemory(e.to_string()))?,
            ),
            reads,
            cap,
            used: 0,
        })
    }

    /// Most reads one `collect` call can take.
    pub fn capacity(&self) -> usize {
        self.reads
    }

    #[allow(clippy::too_many_arguments)]
    pub fn collect(
        &mut self,
        bwt: &GpuBwt,
        table: Option<&GpuKmerTable<'_>>,
        seqs: &[u8],
        offs: &[u64],
        lens: &[u32],
        work: SeedWork,
    ) -> Result<f32, GpuError> {
        assert!(lens.len() <= self.reads);
        self.used = 0;
        let intervals = unsafe {
            std::slice::from_raw_parts_mut(
                self.intervals.as_mut_ptr().cast::<GpuIntv>(),
                self.reads * self.cap as usize,
            )
        };
        let ms = seed_collect_with_options(
            bwt,
            table,
            seqs,
            offs,
            lens,
            self.cap,
            work,
            intervals,
            self.counts.as_u32s_mut(),
        )
        .inspect_err(|e| {
            let buffers = [
                ("seed intervals", &*self.intervals),
                ("seed counts", &*self.counts),
            ];
            diag_gpu_fatal("GPU seed collect", e, &buffers);
        })?;
        self.used = lens.len();
        Ok(ms)
    }

    /// Touch both buffers on the CPU, then promote them to 2 MiB pages
    /// (advisory; see [`nomics_mem::HostBuf::collapse_huge`]). The kernel
    /// writes them in place, so 4 KiB pages cost a GPU translation miss per
    /// read.
    pub fn touch_and_collapse(&mut self) {
        self.intervals.touch_huge("GPU seed intervals");
        self.counts.touch_huge("GPU seed counts");
    }

    /// Read `read`'s interval count from the last `collect` call, or
    /// [`SEED_OVERFLOW`]. Panics past that call's reads.
    pub fn count(&self, read: usize) -> u32 {
        self.counts.as_u32s()[..self.used][read]
    }

    pub fn get(&self, read: usize) -> Option<&[GpuIntv]> {
        let count = self.count(read);
        if count == SEED_OVERFLOW {
            return None;
        }
        let at = read * self.cap as usize;
        let intervals = unsafe {
            std::slice::from_raw_parts(
                self.intervals.as_ptr().cast::<GpuIntv>(),
                self.reads * self.cap as usize,
            )
        };
        Some(&intervals[at..at + count as usize])
    }

    pub fn overflow_count(&self) -> usize {
        self.counts.as_u32s()[..self.used]
            .iter()
            .filter(|&&count| count == SEED_OVERFLOW)
            .count()
    }
}

/// The entry points every GPU backend implements (`nomics_gpu.h`).
#[cfg(any(feature = "cuda", feature = "metal"))]
mod ffi_common {
    use std::ffi::CStr;
    use std::mem::MaybeUninit;
    use std::os::raw::{c_char, c_int, c_void};

    use super::{DeviceInfo, GpuError};

    #[repr(C)]
    pub struct NomicsDeviceInfo {
        pub name: [c_char; 256],
        pub major: c_int,
        pub minor: c_int,
        pub multiprocessors: c_int,
        pub pageable_memory_access: c_int,
        pub total_global_mem: usize,
        pub max_buffer_bytes: usize,
    }

    unsafe extern "C" {
        pub fn nomics_device_info_get(device: c_int, out: *mut NomicsDeviceInfo) -> c_int;
        pub fn nomics_sum_u64_in_place(
            data: *const u64,
            n: usize,
            out_sum: *mut u64,
            out_ms: *mut f32,
        ) -> c_int;
        pub fn nomics_host_register(ptr: *const c_void, len: usize) -> c_int;
        pub fn nomics_host_unregister(ptr: *const c_void) -> c_int;
    }

    pub fn check(code: i32) -> Result<(), GpuError> {
        match code {
            0 => Ok(()),
            #[cfg(feature = "cuda")]
            code => Err(GpuError::Cuda(code)),
            // NM_UNSUPPORTED: a kernel the Metal shim has not ported.
            #[cfg(feature = "metal")]
            1000 => Err(GpuError::Unavailable),
            #[cfg(feature = "metal")]
            code => Err(GpuError::Metal(code)),
        }
    }

    pub fn device_info(device: i32) -> Result<DeviceInfo, GpuError> {
        let mut raw = MaybeUninit::<NomicsDeviceInfo>::uninit();
        check(unsafe { nomics_device_info_get(device, raw.as_mut_ptr()) })?;
        let raw = unsafe { raw.assume_init() };
        let name = unsafe { CStr::from_ptr(raw.name.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        Ok(DeviceInfo {
            name,
            compute_capability: (raw.major, raw.minor),
            multiprocessors: raw.multiprocessors,
            pageable_memory_access: raw.pageable_memory_access != 0,
            total_global_mem: raw.total_global_mem,
            max_buffer_bytes: raw.max_buffer_bytes,
        })
    }

    pub fn sum_u64_in_place(data: &[u64]) -> Result<(u64, f32), GpuError> {
        let (mut sum, mut ms) = (0u64, 0f32);
        check(unsafe { nomics_sum_u64_in_place(data.as_ptr(), data.len(), &mut sum, &mut ms) })?;
        Ok((sum, ms))
    }
}

#[cfg(any(feature = "cuda", feature = "metal"))]
mod ffi {
    use std::os::raw::{c_int, c_void};

    #[repr(C)]
    pub struct NomicsBwt {
        pub primary: u64,
        pub l2: [u64; 5],
        pub seq_len: u64,
        pub words: *const u32,
    }

    #[repr(C)]
    pub struct NomicsKmerTable {
        pub k: u32,
        pub level_offsets: [u64; super::KMER_MAX_K + 1],
        pub intervals: *const u64,
    }

    /// `nomics_seed_scratch`: per-list-size (64, 96, 128) scratch buffers.
    #[repr(C)]
    pub struct NomicsSeedScratch {
        pub lists: [*mut super::GpuIntv; 3],
        pub bytes: [usize; 3],
    }

    #[repr(C)]
    #[derive(Clone, Copy)]
    pub struct NomicsKswJob {
        pub query_offset: u64,
        pub target_offset: u64,
        pub mat_offset: u64,
        pub scratch_offset: u64,
        pub row_offset: u64,
        pub query_len: u32,
        pub target_len: u32,
        pub m: i32,
        pub o_del: i32,
        pub e_del: i32,
        pub o_ins: i32,
        pub e_ins: i32,
        pub xtra: i32,
        pub size: i32,
    }

    unsafe extern "C" {
        pub fn nomics_stream_create(priority: c_int, out: *mut *mut c_void) -> c_int;
        pub fn nomics_stream_destroy(stream: *mut c_void) -> c_int;
        pub fn nomics_stream_warm(stream: *mut c_void) -> c_int;
        pub fn nomics_ticket_wait(ticket: *mut c_void, out_ms: *mut f32) -> c_int;
        pub fn nomics_fused_extend_launch(
            stream: *mut c_void,
            queries: *const u8,
            targets: *const u8,
            params: *const super::GpuFusedParams,
            jobs: *const super::GpuFusedJob,
            n: u32,
            order: *mut u32,
            out: *mut super::GpuFusedResult,
            ticket: *mut *mut c_void,
        ) -> c_int;
        pub fn nomics_ksw_extend2_launch(
            stream: *mut c_void,
            queries: *const u8,
            targets: *const u8,
            matrices: *const i8,
            jobs: *const super::extend::Job,
            n: u32,
            out: *mut super::GpuExtendResult,
            ticket: *mut *mut c_void,
        ) -> c_int;
        pub fn nomics_pairhmm(
            read_bases: *const u8,
            hap_bases: *const u8,
            row_params: *const f64,
            pairs: *const super::GpuPairHmmPair,
            n: u32,
            out_sum: *mut f64,
            out_ms: *mut f32,
        ) -> c_int;
        pub fn nomics_ksw_align2_launch(
            stream: *mut c_void,
            queries: *const u8,
            targets: *const u8,
            matrices: *const i8,
            jobs: *const NomicsKswJob,
            n: u32,
            scratch: *mut u16,
            row_maxima: *mut u16,
            out: *mut super::GpuKswResult,
            ticket: *mut *mut c_void,
        ) -> c_int;
        pub fn nomics_ksw_align2_wave_launch(
            stream: *mut c_void,
            queries: *const u8,
            targets: *const u8,
            matrices: *const i8,
            jobs: *const NomicsKswJob,
            n: u32,
            row_maxima: *mut u16,
            out: *mut super::GpuKswResult,
            ticket: *mut *mut c_void,
        ) -> c_int;
        #[allow(clippy::too_many_arguments)]
        pub fn nomics_seed_collect_launch(
            stream: *mut c_void,
            bwt: *const NomicsBwt,
            table: *const NomicsKmerTable,
            seqs: *const u8,
            offs: *const u64,
            lens: *const u32,
            n: u32,
            cap: u32,
            work_mode: u32,
            out: *mut super::GpuIntv,
            out_n: *mut u32,
            scratch: *const NomicsSeedScratch,
            ticket: *mut *mut c_void,
        ) -> c_int;
        pub fn nomics_seed_scratch_bytes(list_size: u32, bytes: *mut usize) -> c_int;
    }
}

#[cfg(any(feature = "cuda", feature = "metal"))]
mod backend {
    use super::{
        GpuError, GpuExtendResult, GpuKswAlign2Input, GpuKswResult, GpuPairHmmPair,
        PAIRHMM_ROW_PARAMS, PairHmmTimings, ffi,
    };
    use std::os::raw::{c_int, c_void};
    use std::time::{Duration, Instant};

    use nomics_mem::HostBuf;

    use super::ffi_common::check;
    pub use super::ffi_common::{device_info, sum_u64_in_place};

    /// The backend's engine buffers: touched, then registered for life.
    pub type RegBuf = super::GpuHostBuf;

    fn reg_buf(bytes: usize, what: &str) -> Result<RegBuf, GpuError> {
        let mut buf = HostBuf::new(bytes.max(1)).map_err(|e| GpuError::HostMemory(e.to_string()))?;
        buf.touch_huge(what);
        Ok(RegBuf::register(buf))
    }

    fn ensure_registered(
        buffer: &mut Option<RegBuf>,
        bytes: usize,
        what: &str,
    ) -> Result<(), GpuError> {
        if buffer.as_ref().is_some_and(|b| b.len() >= bytes) {
            return Ok(());
        }
        // Rescue buffers are sized by a whole batch (several GB at 30x,
        // touched in full): grow to the need plus 1/4 rather than to the next
        // power of two, which held up to twice the memory per engine. Never
        // above that power of two, so no device buffer limit (Metal's
        // maxBufferLength) is exceeded where the old sizing fit.
        let need = bytes.max(1);
        let capacity = need.saturating_add(need / 4).min(need.next_power_of_two());
        *buffer = None; // unregistered before its replacement is mapped
        *buffer = Some(reg_buf(capacity, what)?);
        Ok(())
    }

    fn ensure_buffer(buffer: &mut Option<HostBuf>, bytes: usize) -> Result<(), GpuError> {
        if buffer.as_ref().is_some_and(|b| b.len() >= bytes) {
            return Ok(());
        }
        let capacity = bytes.max(1).next_power_of_two();
        let mut buf = HostBuf::new(capacity).map_err(|e| GpuError::HostMemory(e.to_string()))?;
        buf.touch_huge("GPU PairHMM buffer");
        *buffer = Some(buf);
        Ok(())
    }

    fn copy_typed<T: Copy>(buffer: &mut HostBuf, values: &[T]) {
        let bytes = std::mem::size_of_val(values);
        assert!(bytes <= buffer.len());
        // HostBuf is 2 MiB aligned and owns `buffer.len()` writable bytes.
        let destination = unsafe {
            std::slice::from_raw_parts_mut(buffer.as_mut_ptr().cast::<T>(), values.len())
        };
        destination.copy_from_slice(values);
    }

    fn typed_mut<T>(buffer: &mut HostBuf, len: usize) -> &mut [T] {
        assert!(len * std::mem::size_of::<T>() <= buffer.len());
        // HostBuf is 2 MiB aligned and owns `buffer.len()` writable bytes.
        unsafe { std::slice::from_raw_parts_mut(buffer.as_mut_ptr().cast::<T>(), len) }
    }

    pub struct PairHmm {
        read_bases: Option<HostBuf>,
        hap_bases: Option<HostBuf>,
        row_params: Option<HostBuf>,
        pairs: Option<HostBuf>,
        out: Option<HostBuf>,
    }

    impl PairHmm {
        pub fn new() -> Self {
            Self {
                read_bases: None,
                hap_bases: None,
                row_params: None,
                pairs: None,
                out: None,
            }
        }

        pub fn compute(
            &mut self,
            read_bases: &[u8],
            hap_bases: &[u8],
            row_params: &[[f64; PAIRHMM_ROW_PARAMS]],
            pairs: &[GpuPairHmmPair],
            out_sum: &mut [f64],
        ) -> Result<PairHmmTimings, GpuError> {
            if pairs.is_empty() {
                return Ok(PairHmmTimings::default());
            }

            let prep_start = Instant::now();
            ensure_buffer(&mut self.read_bases, read_bases.len())?;
            ensure_buffer(&mut self.hap_bases, hap_bases.len())?;
            ensure_buffer(&mut self.row_params, std::mem::size_of_val(row_params))?;
            ensure_buffer(&mut self.pairs, std::mem::size_of_val(pairs))?;
            ensure_buffer(&mut self.out, pairs.len() * std::mem::size_of::<f64>())?;
            self.read_bases.as_mut().unwrap().as_mut_slice()[..read_bases.len()]
                .copy_from_slice(read_bases);
            self.hap_bases.as_mut().unwrap().as_mut_slice()[..hap_bases.len()]
                .copy_from_slice(hap_bases);
            copy_typed(self.row_params.as_mut().unwrap(), row_params);
            copy_typed(self.pairs.as_mut().unwrap(), pairs);
            let prep = prep_start.elapsed();

            let mut kernel_ms = 0.0f32;
            check(unsafe {
                ffi::nomics_pairhmm(
                    self.read_bases.as_ref().unwrap().as_ptr(),
                    self.hap_bases.as_ref().unwrap().as_ptr(),
                    self.row_params.as_ref().unwrap().as_ptr().cast(),
                    self.pairs.as_ref().unwrap().as_ptr().cast(),
                    pairs.len() as u32,
                    self.out.as_mut().unwrap().as_mut_ptr().cast(),
                    &mut kernel_ms,
                )
            })?;

            let post_start = Instant::now();
            let result = unsafe {
                std::slice::from_raw_parts(
                    self.out.as_ref().unwrap().as_ptr().cast::<f64>(),
                    pairs.len(),
                )
            };
            out_sum[..pairs.len()].copy_from_slice(result);
            let post = post_start.elapsed();
            Ok(PairHmmTimings {
                prep,
                kernel: Duration::from_secs_f64(f64::from(kernel_ms) * 1e-3),
                post,
            })
        }
    }

    /// A nonblocking CUDA stream, created once by its owner and reused for
    /// every launch (`nomics_stream_create`). Dropping it waits for its work.
    pub struct Stream(*mut c_void);

    // CUDA streams may be used and destroyed from any thread.
    unsafe impl Send for Stream {}
    unsafe impl Sync for Stream {}

    impl Stream {
        /// `high`: the greatest stream priority, else the least (the default).
        pub fn new(high: bool) -> Result<Self, GpuError> {
            let mut raw = std::ptr::null_mut();
            check(unsafe { ffi::nomics_stream_create(c_int::from(high), &mut raw) })?;
            Ok(Self(raw))
        }

        /// Run an empty kernel on this stream and wait (`nomics_stream_warm`).
        pub fn warm(&self) -> Result<(), GpuError> {
            check(unsafe { ffi::nomics_stream_warm(self.0) })
        }
    }

    impl Drop for Stream {
        fn drop(&mut self) {
            unsafe { ffi::nomics_stream_destroy(self.0) };
        }
    }

    /// One launch in flight (`nomics_ticket`). Its owner keeps every buffer
    /// the launch touches alive and untouched until [`Pending::wait`]
    /// returns; dropping it unwaited waits (and drains on error) too, so
    /// owners declare it before the buffers it covers.
    pub struct Pending(*mut c_void);

    // A ticket's events may be waited on from any thread.
    unsafe impl Send for Pending {}

    impl Pending {
        /// Nothing launched: waits as an immediate success.
        pub fn none() -> Self {
            Self(std::ptr::null_mut())
        }

        fn launched(launch: impl FnOnce(*mut *mut c_void) -> c_int) -> Result<Self, GpuError> {
            let mut raw = std::ptr::null_mut();
            // On error the shim has drained its stream and returns no ticket.
            check(launch(&mut raw))?;
            Ok(Self(raw))
        }

        /// The launch's GPU time once it has finished.
        pub fn wait(mut self) -> Result<Duration, GpuError> {
            let raw = std::mem::replace(&mut self.0, std::ptr::null_mut());
            let mut ms = 0f32;
            let code = unsafe { ffi::nomics_ticket_wait(raw, &mut ms) };
            // NOMICS_DIAG: every live HostBuf's canaries after every launch.
            nomics_mem::check_all_canaries("a GPU ticket wait");
            check(code)?;
            Ok(Duration::from_secs_f64(f64::from(ms) * 1e-3))
        }
    }

    impl Drop for Pending {
        fn drop(&mut self) {
            if !self.0.is_null() {
                let code = unsafe { ffi::nomics_ticket_wait(self.0, std::ptr::null_mut()) };
                nomics_mem::check_all_canaries("a dropped GPU ticket's wait");
                if let Err(e) = check(code) {
                    super::diag_gpu_error("a dropped (unwaited) GPU ticket", &e);
                }
            }
        }
    }

    /// Enqueue `ksw_extend2` for `jobs`. The caller keeps every slice alive
    /// and untouched until the returned ticket is waited.
    pub fn ksw_extend2_launch(
        stream: &Stream,
        queries: &[u8],
        targets: &[u8],
        matrices: &[u8],
        jobs: &[super::extend::Job],
        out: &mut [GpuExtendResult],
    ) -> Result<Pending, GpuError> {
        assert!(out.len() >= jobs.len());
        Pending::launched(|ticket| unsafe {
            ffi::nomics_ksw_extend2_launch(
                stream.0,
                queries.as_ptr(),
                targets.as_ptr(),
                matrices.as_ptr().cast(),
                jobs.as_ptr(),
                jobs.len() as u32,
                out.as_mut_ptr(),
                ticket,
            )
        })
    }

    /// Enqueue the fused extension for `jobs`. The caller keeps every slice
    /// alive and untouched until the returned ticket is waited; `params` is
    /// copied.
    pub fn fused_extend_launch(
        stream: &Stream,
        queries: &[u8],
        targets: &[u8],
        params: &super::GpuFusedParams,
        jobs: &[super::GpuFusedJob],
        order: &mut [u32],
        out: &mut [super::GpuFusedResult],
    ) -> Result<Pending, GpuError> {
        assert_eq!(jobs.len(), out.len());
        assert_eq!(jobs.len(), order.len());
        Pending::launched(|ticket| unsafe {
            ffi::nomics_fused_extend_launch(
                stream.0,
                queries.as_ptr(),
                targets.as_ptr(),
                params,
                jobs.as_ptr(),
                jobs.len() as u32,
                order.as_mut_ptr(),
                out.as_mut_ptr(),
                ticket,
            )
        })
    }

    pub struct KswAlign2 {
        queries: Option<RegBuf>,
        targets: Option<RegBuf>,
        matrices: Option<RegBuf>,
        jobs: Option<RegBuf>,
        scratch: Option<RegBuf>,
        rows: Option<RegBuf>,
        out: Option<RegBuf>,
        /// Created on first launch: the least priority (as the legacy
        /// default stream had), or the greatest with NOMICS_RESCUE_PRIORITY=1.
        stream: Option<Stream>,
    }

    impl KswAlign2 {
        fn ensure_stream(&mut self) -> Result<(), GpuError> {
            if self.stream.is_none() {
                let high = std::env::var("NOMICS_RESCUE_PRIORITY").is_ok_and(|v| v == "1");
                self.stream = Some(Stream::new(high)?);
            }
            Ok(())
        }

        /// Create the stream now and run an empty kernel on it.
        pub fn warm(&mut self) -> Result<(), GpuError> {
            self.ensure_stream()?;
            self.stream.as_ref().unwrap().warm()
        }

        pub fn new() -> Self {
            Self {
                queries: None,
                targets: None,
                matrices: None,
                jobs: None,
                scratch: None,
                rows: None,
                out: None,
                stream: None,
            }
        }

        /// This engine's buffers by field, for diagnostic reports.
        pub fn diag_buffers(&self) -> Vec<(&'static str, &HostBuf)> {
            [
                ("rescue queries", &self.queries),
                ("rescue targets", &self.targets),
                ("rescue matrices", &self.matrices),
                ("rescue jobs", &self.jobs),
                ("rescue scratch", &self.scratch),
                ("rescue row maxima", &self.rows),
                ("rescue results", &self.out),
            ]
            .into_iter()
            .filter_map(|(name, b)| b.as_ref().map(|b| (name, &**b)))
            .collect()
        }

        /// Bytes of host buffers held (diagnostics).
        pub fn bytes(&self) -> usize {
            [
                &self.queries,
                &self.targets,
                &self.matrices,
                &self.jobs,
                &self.scratch,
                &self.rows,
                &self.out,
            ]
            .iter()
            .filter_map(|b| b.as_ref().map(|b| b.len()))
            .sum()
        }

        /// The results of the last waited launch of `n` inputs.
        pub fn results(&self, n: usize) -> &[GpuKswResult] {
            if n == 0 {
                return &[];
            }
            let out = self.out.as_ref().unwrap();
            assert!(n * std::mem::size_of::<GpuKswResult>() <= out.len());
            // repr(C) i32 fields, any bit pattern valid; HostBuf is 2 MiB aligned.
            unsafe { std::slice::from_raw_parts(out.as_ptr().cast::<GpuKswResult>(), n) }
        }

        /// Pack `inputs` into this engine's buffers and enqueue them on its
        /// stream. The ticket covers the buffers: the caller keeps this
        /// engine with it, untouched, until it is waited. Returns the ticket
        /// and the packing time.
        pub fn launch(
            &mut self,
            inputs: &[GpuKswAlign2Input<'_>],
        ) -> Result<(Pending, Duration), GpuError> {
            if inputs.is_empty() {
                return Ok((Pending::none(), Duration::ZERO));
            }
            self.ensure_stream()?;
            let prep_start = Instant::now();
            let mut query_bytes = 0usize;
            let mut target_bytes = 0usize;
            let mut matrix_bytes = 0usize;
            let mut scratch_words = 0usize;
            let mut row_words = 0usize;
            for input in inputs {
                query_bytes = query_bytes
                    .checked_add(input.query.len())
                    .expect("query bytes");
                target_bytes = target_bytes
                    .checked_add(input.target.len())
                    .expect("target bytes");
                matrix_bytes = matrix_bytes
                    .checked_add(input.mat.len())
                    .expect("matrix bytes");
                let lanes = if input.size == 1 { 16 } else { 8 };
                let padded = input.query.len().div_ceil(lanes) * lanes;
                scratch_words = scratch_words
                    .checked_add(4usize.checked_mul(padded).expect("KSW scratch"))
                    .expect("KSW scratch");
                row_words = row_words
                    .checked_add(input.target.len())
                    .expect("KSW row maxima");
            }
            // The wavefront kernel; NOMICS_KSW_STRIPED=1 selects the striped one.
            let wave = std::env::var_os("NOMICS_KSW_STRIPED").is_none();
            const WHAT: &str = "GPU mate-rescue buffer";
            ensure_registered(&mut self.queries, query_bytes, WHAT)?;
            ensure_registered(&mut self.targets, target_bytes, WHAT)?;
            ensure_registered(&mut self.matrices, matrix_bytes, WHAT)?;
            let job_bytes = std::mem::size_of::<ffi::NomicsKswJob>() * inputs.len();
            ensure_registered(&mut self.jobs, job_bytes, WHAT)?;
            // Only the striped kernel uses the scratch and row-maxima buffers
            // (the wavefront kernel keeps rows in shared memory and ignores
            // its `row_maxima` argument), so the default path never maps
            // them: at 30x they were ~3/4 of each engine's memory.
            if !wave {
                let scratch_bytes = std::mem::size_of::<u16>() * scratch_words;
                ensure_registered(&mut self.scratch, scratch_bytes, WHAT)?;
                ensure_registered(&mut self.rows, std::mem::size_of::<u16>() * row_words, WHAT)?;
            }
            let out_bytes = std::mem::size_of::<GpuKswResult>() * inputs.len();
            ensure_registered(&mut self.out, out_bytes, WHAT)?;

            let query_buffer = self.queries.as_mut().unwrap().as_mut_slice();
            let target_buffer = self.targets.as_mut().unwrap().as_mut_slice();
            let matrix_buffer = self.matrices.as_mut().unwrap().as_mut_slice();
            let jobs = typed_mut::<ffi::NomicsKswJob>(self.jobs.as_mut().unwrap(), inputs.len());
            let (mut query_offset, mut target_offset, mut mat_offset) = (0usize, 0usize, 0usize);
            let (mut scratch_offset, mut row_offset) = (0usize, 0usize);
            for (raw, input) in jobs.iter_mut().zip(inputs) {
                query_buffer[query_offset..query_offset + input.query.len()]
                    .copy_from_slice(input.query);
                target_buffer[target_offset..target_offset + input.target.len()]
                    .copy_from_slice(input.target);
                for (slot, &score) in matrix_buffer[mat_offset..mat_offset + input.mat.len()]
                    .iter_mut()
                    .zip(input.mat)
                {
                    *slot = score as u8;
                }
                let lanes = if input.size == 1 { 16 } else { 8 };
                let padded = input.query.len().div_ceil(lanes) * lanes;
                *raw = ffi::NomicsKswJob {
                    query_offset: query_offset as u64,
                    target_offset: target_offset as u64,
                    mat_offset: mat_offset as u64,
                    scratch_offset: scratch_offset as u64,
                    row_offset: row_offset as u64,
                    query_len: input.query.len() as u32,
                    target_len: input.target.len() as u32,
                    m: input.m,
                    o_del: input.o_del,
                    e_del: input.e_del,
                    o_ins: input.o_ins,
                    e_ins: input.e_ins,
                    xtra: input.xtra,
                    size: i32::from(input.size),
                };
                query_offset += input.query.len();
                target_offset += input.target.len();
                mat_offset += input.mat.len();
                scratch_offset += 4 * padded;
                row_offset += input.target.len();
            }
            let prep = prep_start.elapsed();

            let stream = self.stream.as_ref().unwrap().0;
            let pending = Pending::launched(|ticket| unsafe {
                if wave {
                    ffi::nomics_ksw_align2_wave_launch(
                        stream,
                        self.queries.as_ref().unwrap().as_ptr(),
                        self.targets.as_ref().unwrap().as_ptr(),
                        self.matrices.as_ref().unwrap().as_ptr().cast(),
                        self.jobs.as_ref().unwrap().as_ptr().cast(),
                        inputs.len() as u32,
                        // Unused by both shims' wavefront kernels (nomics_gpu.h).
                        std::ptr::null_mut(),
                        self.out.as_mut().unwrap().as_mut_ptr().cast(),
                        ticket,
                    )
                } else {
                    ffi::nomics_ksw_align2_launch(
                        stream,
                        self.queries.as_ref().unwrap().as_ptr(),
                        self.targets.as_ref().unwrap().as_ptr(),
                        self.matrices.as_ref().unwrap().as_ptr().cast(),
                        self.jobs.as_ref().unwrap().as_ptr().cast(),
                        inputs.len() as u32,
                        self.scratch.as_mut().unwrap().as_mut_ptr().cast(),
                        self.rows.as_mut().unwrap().as_mut_ptr().cast(),
                        self.out.as_mut().unwrap().as_mut_ptr().cast(),
                        ticket,
                    )
                }
            })?;
            Ok((pending, prep))
        }
    }

    pub fn seed_collect(
        bwt: &super::GpuBwt,
        table: Option<&super::GpuKmerTable<'_>>,
        seqs: &[u8],
        offs: &[u64],
        lens: &[u32],
        cap: u32,
        work: super::SeedWork,
        out: &mut [super::GpuIntv],
        out_n: &mut [u32],
    ) -> Result<f32, GpuError> {
        let raw = ffi::NomicsBwt {
            primary: bwt.primary,
            l2: bwt.l2,
            seq_len: bwt.seq_len,
            words: bwt.words.as_ptr(),
        };
        let raw_table = table.map(|table| ffi::NomicsKmerTable {
            k: table.k,
            level_offsets: table.level_offsets,
            intervals: table.intervals.as_ptr(),
        });
        let mut scratch = ffi::NomicsSeedScratch {
            lists: [std::ptr::null_mut(); 3],
            bytes: [0; 3],
        };
        let sizes: &[usize] = match work {
            super::SeedWork::Cap64 => &[0],
            super::SeedWork::Cap96 => &[1],
            super::SeedWork::Cap128 => &[2],
            super::SeedWork::TwoPass64Then128 => &[0, 2],
        };
        for &at in sizes {
            if let Some((ptr, len)) = seed_scratch(at) {
                // The kernel writes this memory; Rust never reads or writes it.
                scratch.lists[at] = ptr as *mut super::GpuIntv;
                scratch.bytes[at] = len;
            }
        }
        // Lengths and bounds were checked by the caller (`super::seed_collect`).
        // Every slice is borrowed until the ticket is waited below; `raw` and
        // `raw_table` are copied at launch.
        let stream = seed_stream()?;
        let pending = Pending::launched(|ticket| unsafe {
            ffi::nomics_seed_collect_launch(
                stream.0,
                &raw,
                raw_table
                    .as_ref()
                    .map_or(std::ptr::null(), |table| table as *const _),
                seqs.as_ptr(),
                offs.as_ptr(),
                lens.as_ptr(),
                lens.len() as u32,
                cap,
                work.abi(),
                out.as_mut_ptr(),
                out_n.as_mut_ptr(),
                &scratch,
                ticket,
            )
        })?;
        Ok(pending.wait()?.as_secs_f32() * 1e3)
    }

    /// The process's seeding stream, created once: nonblocking (seeding no
    /// longer serializes with rescue on the legacy default stream) at the
    /// least priority, as before. One stream keeps seed kernels in order, one
    /// resident at a time, which the scratch slot pool is sized for.
    /// Create the seeding stream now and run an empty kernel on it.
    pub fn warm_seed_stream() -> Result<(), GpuError> {
        seed_stream()?.warm()
    }

    fn seed_stream() -> Result<&'static Stream, GpuError> {
        static STREAM: std::sync::OnceLock<Result<Stream, GpuError>> = std::sync::OnceLock::new();
        STREAM
            .get_or_init(|| Stream::new(false))
            .as_ref()
            .map_err(Clone::clone)
    }

    /// Process-wide scratch for the seed kernel's interval lists, one buffer
    /// per list size (64, 96, 128), allocated on first use: sized by how many
    /// blocks can be resident (302 MB for 64 entries on the GB10: 384
    /// slots), touched and collapsed to 2 MiB pages. Concurrent calls share it (the kernel hands
    /// out disjoint slots). `None` (local-memory kernel, same output) when
    /// the allocation fails or `NOMICS_GPU_SEED_SCRATCH=0`.
    fn seed_scratch(at: usize) -> Option<(usize, usize)> {
        static POOLS: [std::sync::OnceLock<Option<(usize, usize)>>; 3] =
            [const { std::sync::OnceLock::new() }; 3];
        *POOLS[at].get_or_init(|| {
            if std::env::var_os("NOMICS_GPU_SEED_SCRATCH").is_some_and(|v| v == "0") {
                return None;
            }
            let mut bytes = 0usize;
            let size = [64u32, 96, 128][at];
            if unsafe { ffi::nomics_seed_scratch_bytes(size, &mut bytes) } != 0 {
                return None;
            }
            // Registered (Metal binds it without re-wrapping) and never freed.
            let buf: &'static mut RegBuf =
                Box::leak(Box::new(reg_buf(bytes, "GPU seed scratch").ok()?));
            buf.as_mut_slice().fill(0);
            let _ = buf.collapse_huge();
            // On 4 KiB pages the steady kernel time is unchanged (measured),
            // but the GPU's first pass over them costs ~20 s once per
            // process. Still far better than the local-memory kernel.
            if let Some(huge) = buf.huge_page_bytes()
                && huge + (2 << 20) < buf.len()
            {
                eprintln!(
                    "WARNING: GPU seed scratch {:.0} MB, only {:.0}% on 2 MiB pages \
                     (the first seed batch will be slow; output is unaffected)",
                    buf.len() as f64 / 1e6,
                    100.0 * huge as f64 / buf.len() as f64
                );
            }
            Some((buf.as_ptr() as usize, buf.len()))
        })
    }
}

#[cfg(not(any(feature = "cuda", feature = "metal")))]
mod backend {
    use super::{
        DeviceInfo, GpuError, GpuExtendResult, GpuKswAlign2Input, GpuKswResult, GpuPairHmmPair,
        PAIRHMM_ROW_PARAMS, PairHmmTimings,
    };
    use std::time::Duration;

    pub struct Stream;

    impl Stream {
        pub fn new(_high: bool) -> Result<Self, GpuError> {
            Err(GpuError::Unavailable)
        }

        pub fn warm(&self) -> Result<(), GpuError> {
            Err(GpuError::Unavailable)
        }
    }

    pub struct Pending;

    impl Pending {
        pub fn wait(self) -> Result<Duration, GpuError> {
            Ok(Duration::ZERO)
        }
    }

    pub fn ksw_extend2_launch(
        _stream: &Stream,
        _queries: &[u8],
        _targets: &[u8],
        _matrices: &[u8],
        _jobs: &[super::extend::Job],
        _out: &mut [GpuExtendResult],
    ) -> Result<Pending, GpuError> {
        Err(GpuError::Unavailable)
    }

    pub fn fused_extend_launch(
        _stream: &Stream,
        _queries: &[u8],
        _targets: &[u8],
        _params: &super::GpuFusedParams,
        _jobs: &[super::GpuFusedJob],
        _order: &mut [u32],
        _out: &mut [super::GpuFusedResult],
    ) -> Result<Pending, GpuError> {
        Err(GpuError::Unavailable)
    }

    pub struct PairHmm;

    impl PairHmm {
        pub fn new() -> Self {
            Self
        }

        pub fn compute(
            &mut self,
            _read_bases: &[u8],
            _hap_bases: &[u8],
            _row_params: &[[f64; PAIRHMM_ROW_PARAMS]],
            _pairs: &[GpuPairHmmPair],
            _out_sum: &mut [f64],
        ) -> Result<PairHmmTimings, GpuError> {
            Err(GpuError::Unavailable)
        }
    }

    pub struct KswAlign2;

    pub fn warm_seed_stream() -> Result<(), GpuError> {
        Err(GpuError::Unavailable)
    }

    impl KswAlign2 {
        pub fn new() -> Self {
            Self
        }

        pub fn warm(&mut self) -> Result<(), GpuError> {
            Err(GpuError::Unavailable)
        }

        pub fn results(&self, _n: usize) -> &[GpuKswResult] {
            &[]
        }

        pub fn bytes(&self) -> usize {
            0
        }

        pub fn diag_buffers(&self) -> Vec<(&'static str, &nomics_mem::HostBuf)> {
            Vec::new()
        }

        pub fn launch(
            &mut self,
            _inputs: &[GpuKswAlign2Input<'_>],
        ) -> Result<(Pending, Duration), GpuError> {
            Err(GpuError::Unavailable)
        }
    }

    pub fn device_info(_device: i32) -> Result<DeviceInfo, GpuError> {
        Err(GpuError::Unavailable)
    }

    pub fn sum_u64_in_place(_data: &[u64]) -> Result<(u64, f32), GpuError> {
        Err(GpuError::Unavailable)
    }

    pub fn seed_collect(
        _bwt: &super::GpuBwt,
        _table: Option<&super::GpuKmerTable<'_>>,
        _seqs: &[u8],
        _offs: &[u64],
        _lens: &[u32],
        _cap: u32,
        _work: super::SeedWork,
        _out: &mut [super::GpuIntv],
        _out_n: &mut [u32],
    ) -> Result<f32, GpuError> {
        Err(GpuError::Unavailable)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(not(any(feature = "cuda", feature = "metal")))]
    #[test]
    fn stub_reports_unavailable() {
        assert_eq!(sum_u64_in_place(&[1, 2, 3]), Err(GpuError::Unavailable));
        let input = GpuKswAlign2Input {
            query: &[0],
            target: &[0],
            m: 1,
            mat: &[1],
            o_del: 1,
            e_del: 1,
            o_ins: 1,
            e_ins: 1,
            xtra: 0,
            size: 2,
        };
        let mut engine = GpuKswAlign2::new();
        let mut out = [GpuKswResult::default()];
        assert!(matches!(
            engine.compute(&[input], &mut out),
            Err(GpuError::Unavailable)
        ));
    }

    #[cfg(any(feature = "cuda", feature = "metal"))]
    #[test]
    fn gpu_sums_host_buffer_in_place() {
        let info = device_info(0).unwrap();
        assert!(
            info.pageable_memory_access,
            "{} cannot read host pointers",
            info.name
        );
        let mut buf = nomics_mem::HostBuf::new(64 << 20).unwrap();
        for (i, w) in buf.as_u64s_mut().iter_mut().enumerate() {
            *w = i as u64;
        }
        let n = buf.as_u64s().len() as u64;
        let (sum, _ms) = sum_u64_in_place(buf.as_u64s()).unwrap();
        assert_eq!(sum, n * (n - 1) / 2);
    }
}
