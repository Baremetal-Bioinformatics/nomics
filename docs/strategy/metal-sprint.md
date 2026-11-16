# Metal sprint: Apple Silicon port, 24 h rented window (plan, 2026-10-21)

*The cloud and kit scripts this document names (`scripts/cloud/`, `scripts/kit/`) drive one owner's AWS account and are not part of the public repository; the measurements and code changes they produced are.*

Goal: an exact `nomics run` (FASTQ → VCF, every gate identical to sarek) on Apple Silicon.
It needs a measured CPU control arm, and a Metal path that ships only if it beats that arm
(docs/METHODOLOGY.md). The rented box is for full-scale runs. Everything that fits in 24 GB is done
on the local M4 Pro first, before the 24 h clock starts.

## Where we start

- **CPU path is mostly portable already.** NEON kernels (PairHMM rows, ksw u8/i16, HC
  wavefront SW) have portable twins. The aligner's libm is ported (glibc `log/exp/erfc`,
  bit-exact on macOS) and HC uses its own fdlibm. The workspace builds with the GPU stub on
  the Mac. A full run has never been done on macOS.
- **GPU code is CUDA-only.** About 2.5 k lines of `.cu` sit behind a C ABI
  (`crates/nomics-gpu/shim/nomics_gpu.h`), and Rust only sees that header. A Metal backend
  can implement the same ABI, which leaves the callers in `nomics-align` and `nomics-hc`
  untouched.
- **Memory.** The tier1 exome in memory peaks at 57.6 GB (`bench/full.md`). The FM-index
  is 5.6 GB and the dense SA 8.0 GB. Whole-genome mode peaked at 22.6 GB on wgs20 and
  74.5–80 GB at 30×.
- **Exactness gate without shipping goldens.** `scripts/kit/` + `nomics dev digest` checks
  41 normalized artifact digests. `run.sh` is Linux/CUDA-specific: it uses `stat -c`,
  GNU time, numactl and nvidia-smi.
- **Kernel oracles.** `crates/nomics-align/examples/{ksw_bench,extend_bench}.rs` replay
  dumped calls (`NOMICS_DUMP_KSW`, `extend2.bin`) and compare the GPU against the dump with
  0 mismatches. The dumps live on the Spark.

## Metal vs CUDA, what matters for our kernels

| Concern | CUDA (today) | Metal | Consequence |
|---|---|---|---|
| FP64 | used in `pairhmm.cu` (30 lines), `extend_*.cuh` `max_ins/max_del` | **none** | PairHMM stays on CPU (NEON). `max_ins/max_del` depend only on qlen + scoring: compute host-side in Rust, pass in the plan struct. |
| Warp ops | `__shfl*` (ksw, extend, pairhmm), `__syncwarp`, `__popc` | `simd_shuffle*`, `simdgroup_barrier`, `popcount` | Direct map. Check `threadExecutionWidth == 32` at runtime instead of assuming it. |
| Atomics + spin | `atomicCAS` slot claim + `__nanosleep` spin (`seed.cu:429`) | 32-bit atomics; no sleep; no forward-progress guarantee for spin-waits | **Redesign, don't translate.** Fetch-add slot claiming or a two-pass count/scan. |
| Streams/priority | `cudaStreamCreateWithPriority`, 3 streams | `MTLCommandQueue` per stream; no priorities | Seed/rescue/extend each get their own queue. Seeding's priority over rescue becomes an ordering question. |
| Async tickets | `cudaEvent` (`ticket.cuh`, `stream.cu`) | `MTLSharedEvent` or `addCompletedHandler` | Same ticket ABI; `nomics_ticket_wait` blocks on the event value. |
| Zero-copy host memory | GPU reads HostBuf in place (2 MiB THP) | `newBufferWithBytesNoCopy` (page-aligned ptr and length) | HostBuf is already 2 MiB aligned and rounded. Watch `maxBufferLength` (index + SA = 13.6 GB; may need splitting). |
| Page size | 4 K + THP 2 MiB, `MADV_HUGEPAGE` | 16 K pages, no THP advice | HostBuf degrades gracefully. The huge-page warnings must not fire as false alarms on macOS. |
| Shared mem | `__shared__` (fused, ksw, seed) | `threadgroup` memory, 32 KB | Check ksw/fused tile sizes fit. |
| Build | nvcc in `build.rs`, `--features cuda` | `xcrun metal` → `.metallib`, ObjC shim via clang, `--features metal` | Metal toolchain is a separate component on macOS 26 (present locally). |
| Profiling | nsys | `xctrace` Metal System Trace (headless) | Learn it locally; it is slow to set up on the clock. |

What does **not** port: FP64 stages (PairHMM, bwa MAPQ/pairing, HC genotyping). They stay
on CPU, as the original plan decided.

GPU targets, in value order: seeding (`seed_collect`, 150 GPU-s per 30× chunk on the Spark)
and mate rescue (`ksw_align2`). The fused extension does not pay yet even on the Spark
(`gpu-alignment-scope.md`), so it is stretch only.

## Hardware choice

| Type | CPU | GPU | Memory | $/h | 24 h |
|---|---|---|---|---:|---:|
| local M4 Pro | 14 | 20-core | 24 GB | 0 | 0 |
| `mac-m4pro.metal` | 14 | 20-core | 48 GB | 1.97 | $47 |
| **`mac-m4max.metal`** | 16 | 40-core | 128 GB | 6.25 | **$150** |
| `mac-m3ultra.metal` | 28 | 80-core | 256 GB | 12.50 | $300 |

**M4 Max.** 128 GB holds the in-memory tier1 config the Spark uses. It is the newest GPU
family and the machine people actually own. The M4 Pro is the same chip as the local Mac,
and 48 GB cannot hold tier1 in memory. The M3 Ultra only earns its price for a "Mac Studio
vs Spark" throughput headline.

The 24 h minimum is real. Billing starts at **host allocation**, not instance launch.

**Quota.** All Mac host quotas were 0. M4 Max quota increases (1 host) were requested
2026-10-21 in us-east-1 and us-west-2. Both are PENDING; they
usually take days. M4 Max is offered in only two AZs per region, so expect capacity
errors. Check with:
`aws service-quotas list-requested-service-quota-change-history-by-quota --service-code ec2 --quota-code L-D82CB68A --region us-east-1`.

## Phase A: local M4 Pro, while the quota is pending (no clock)

### Progress

- **A1 done (2026-10-21).** `cargo test --workspace --release` passes on macOS (M4 Pro). The
  only failure was `nomics-hc/tests/annotate.rs`, whose fixture lived in an old session's
  scratch dir; it is now `#[ignore]` like the other local-oracle tests. HostBuf's huge-page
  warnings already stay quiet off Linux.
- **A3 done (2026-10-21).** `--features metal` (on `nomics-gpu` and `nomics`): ObjC shim
  `shim/metal/nomics_metal.m` implements `nomics_device_info_get` and
  `nomics_sum_u64_in_place` from `nomics_gpu.h`. Kernels are compiled with
  `-fno-fast-math` into one metallib, embedded in the static lib. `nomics selftest` on
  the local M4 Pro (Apple9 family, 20 GPU cores):
  - in-place GPU sum of a HostBuf: **231–236 GiB/s** (1–12 GiB; Spark 154 GiB/s), exact;
  - `recommendedMaxWorkingSetSize` 17.8 GiB, **`maxBufferLength` 13.3 GiB**: the FM-index
    (5.6 GB) and dense SA (8.0 GB) fit as separate buffers;
  - **wrapping is not free: ~17 ms per GiB per call** (12 GiB: 260 ms call for 52 ms of GPU
    time). Resident buffers (index, SA, scratch pools) must be wrapped once. That needs a
    register/unregister ABI tied to HostBuf's lifetime, not a transparent cache: a freed
    and remapped address would leave a stale wrapper. Design it before A5/A6.
  - `nomics_gpu.h` gained `max_buffer_bytes` (CUDA sets `SIZE_MAX`; not yet compiled on
    the Spark).
- **Register-once ABI done (2026-10-21).** `nomics_host_register/unregister` in
  `nomics_gpu.h` (CUDA: no-ops in `stream.cu`). The Metal shim wraps a registered range
  once and binds it for any launch pointer inside it. Rust owns registration through
  `RegBuf` (a HostBuf registered for its whole life; dropped only after its launches are
  waited). The same Rust backend module now builds for `cuda` or `metal`. Kernels the
  Metal shim has not ported return status 1000, which maps to `GpuError::Unavailable`,
  so callers keep their CPU path exactly as with the stub. Streams are `MTLCommandQueue`s
  (no priorities). Tickets are committed command buffers.
- **A4 done (2026-10-21): seeding is a go on Metal.** `bench/experiments/occ_microbench_metal.mm`
  is the twin of `occ_microbench.cu`: dependent random `occ4` lookups (LF-walk chains)
  over the real 3.2 GB `.bwt`, read in place.

  | | lookups/s |
  |---|---:|
  | M4 Pro CPU, 14 threads | 97 M |
  | **M4 Pro GPU (20 cores), 1M–4M chains** | **2.72 G (~28×)** |
  | Spark CPU, 20 threads | 149 M |
  | Spark GB10 GPU | 2.36 G |

  All 4,096 sampled chains equal the CPU walk. The no-copy wrap of 3.2 GB itself took
  ~0 ms; the first GPU run is ~15% slower, as its pages get mapped. So the A3 "~17 ms per
  GiB per call" is better read as the cost of a launch's first GPU touch of freshly
  wrapped pages. Registering once (above) still removes it from every later launch. On the
  Spark the real kernel ran ~15× below this microbench, because of per-read control-flow
  divergence (`bench/align.md`). Expect the same shape on Metal.
- **A5 done (2026-10-21): mate rescue on Metal, exact.** `shim/metal/ksw.metal` is a
  line-for-line port of `ksw_wave_kernel` (simdgroup shuffles for warp shuffles; the
  launcher checks `threadExecutionWidth == 32`). Threadgroup memory is 32 KiB vs CUDA's
  48 KiB, so the launcher fits up to four alignments per threadgroup and drops to fewer
  for long windows.
  - `ksw_bench` exome `align2.bin`: **0 mismatches / 62,000 calls**. Kernel 0.157 s vs
    0.176 s on 14 CPU threads (M4 Pro).
  - `ksw_bench` wgs20 `align2.bin`: **0 / 70,248**. Kernel 0.274 s vs CPU 0.305 s.
  - `nomics-align/tests/ksw_align2_gpu.rs` matches the CPU port on synthetic calls:
    - windows up to the 4,000-base GPU cap (two alignments per threadgroup);
    - reads of 100–512 bp, so **all five kernel widths** run. The dumps' 151 bp reads use
      only the 4- and 5-column kernels;
    - tandem repeats and windows with no mate copy (weak hits), so tie-breaks matter.
  - Mutation check: the test fails on three injected tie-break changes (column argmax,
    row max, second best). An `ksw_bench` dump also catches the first. The one surviving
    mutant (`rev.gmax == score` → `>= score - 1`) is equivalent in practice, because the
    reverse pass re-scores the same alignment.
  - Every width's pipeline allows 1,024 threads per threadgroup on the M4 Pro
    (`NOMICS_METAL_INFO=1` prints this). The launcher still clamps simdgroups to that
    limit, because dispatching past it would fail a whole batch back to the CPU.
  - In `nomics run`, rescue goes to the GPU only when GPU seeding succeeded; with A6 below,
    `NOMICS_GPU_SEED=1` turns on both under Metal.
- **Cross-model review (gpt-6.1-sol, 2026-10-21) of the shim: three findings, all fixed.**
  - A nil compute encoder is now an error. Before, an empty command buffer would "succeed"
    with stale results.
  - Rescue threadgroup sizing is computed in `size_t` and bounded before narrowing. A
    huge target could wrap it. `ksw.cu` has the same pattern; it was left alone because
    targets above 4,000 never reach it.
  - A failed ticket drains its queue, as `nomics_gpu.h` requires, and
    `nomics_stream_destroy` reports a failed drain.
  - The review found no issue in the registry or ARC lifetimes, the Rust ticket ownership,
    or the kernel-vs-CUDA semantics.
- **A6 done (2026-10-21): seed search on Metal, exact.** `shim/metal/seed.metal` ports
  `collect_slot_kernel`. CUDA's slot pool (an `atomicCAS` spin with `__nanosleep`) is
  replaced by serial dispatches of ≤131,072 reads inside one command buffer, with scratch
  indexed by thread position. That needs no atomics and no forward-progress assumption.
  Seed launches share the scratch, so each waits on the CPU for the previous one before
  committing (deadlock-free on error). The split length is the constant 28: Metal has no
  double, and the CUDA kernel notes it is 28 either way. The scratch pools are
  `RegBuf`s, and there is no local-memory kernel: no scratch means CPU fallback.
  - `seed_gate` on **all of tier0 (991,118 reads)**: every read the GPU finishes is
    identical to the CPU. This holds for cap64 (3,421 overflow to the CPU), two-pass (473),
    two-pass with the k = 12 table (473), and plain lists (`NOMICS_SEED_PACK=0`, same
    counts).
  - Speed on the M4 Pro: two-pass, 662k-read batch in 658 ms (~1.0 µs/read). The Spark
    did tier0 at ~0.73 µs/read.
  - Still to do: register the 3.2 GB BWT words once. Today each launch wraps them, which
    costs only first-touch mapping per batch.
- **A8 done (2026-10-21): tier0 FASTQ → VCF on the local M4 Pro, both arms exact.**
  `scripts/mac_tier0.sh cpu|metal` runs `nomics run` and digests every artifact. It
  compares them with `reference/tier0-wes-20261008-r1.digest.tsv` (41 artifacts, made on
  the Spark from the sarek tier0 golden with `nomics dev digest`), so no goldens ship.

  | tier0, 14 threads | total | align | peak RSS | digests |
  |---|---:|---:|---:|---|
  | CPU, run 1 / 2 | 17.1 / 15.6 s | 9.9 / 8.4 s | 14.9 / 15.0 GB | 41/41, 41/41 |
  | Metal seeds + rescue, run 1 / 2 | 13.4 / 12.9 s | 6.8 / 6.3 s | 15.7 / 16.3 GB | 41/41, 41/41 |

  - ~4 s of align is loading the index (alongside prep), so alignment proper is ~4.5 s
    (CPU) vs ~2.5 s (Metal). Per 82k-read batch, CPU busy time falls from ~10 s to
    ~1.3 s with GPU seeds; rescue runs ~3,000 calls in 0.02 s.
  - The first run failed 14/41 on the BAMs only: bwa's `@RG DS` holds the reference
    *path*, and the golden's is the Spark's. The run passes `--read-group-ds` with the
    golden's path, as the x86 goldens do. **The rented box needs the same flag**
    (tier1 golden: the Spark path).
  - The first run also built the 8 GB dense-SA cache. That took ~23 min on 24 GB, under
    swap (Spark: ~4 min). For the rented box, put `sa4.GATK.GRCh38.nsa` in the S3 kit or
    budget the build in hour 0–1.

Each item ends in something checkable. The kernel ports (A5–A7) are independent files with
their own oracle gates, so they can run as parallel workers once A3 lands. They are
exactness-critical ports: `gpt-6.1-sol` at high/max, lead reviews diffs and gate output.

1. **macOS CPU build and tests.** `cargo test --workspace --release` on the Mac with no GPU
   features. Silence the THP and huge-page warnings that are false alarms at 16 K pages.
2. **Kit on macOS.** Port `scripts/kit/run.sh` (BSD `stat`, `/usr/bin/time -l`,
   `sysctl hw.ncpu`, no numactl, a Metal device check instead of nvidia-smi). Done when
   `--digest-only` works on macOS against a Spark run dir.
3. **Metal backend skeleton.** Add the `metal` feature in `nomics-gpu`: `shim/*.metal`,
   an ObjC shim implementing `nomics_gpu.h`, and `build.rs` compiling the metallib. Queues
   act as streams and `MTLSharedEvent` as tickets. Smoke test: `sum_u64_in_place` on a
   HostBuf via no-copy, with bandwidth printed (Spark: 154 GiB/s; M4 Pro peak ~273 GB/s).
   **Record `maxBufferLength` and `recommendedMaxWorkingSetSize`.** If 13.6 GB cannot be one
   no-copy buffer, split the index before designing seed_collect.
4. **Dependent-lookup microbench** (the original plan asked for it before porting seed search):
   FM-index occ chasing, Metal vs NEON CPU, random-access latency at index scale.
5. **ksw mate rescue → Metal.** Integer and simd shuffles. Gate: `ksw_bench` on the Spark
   `align2.bin` dumps, 0 mismatches.
6. **seed_collect → Metal.** Redesign the slot claim (see the table). Gate: Metal seeds equal
   CPU seeds, interval for interval, on tier0 reads. The full GRCh38 index (13.6 GB) on
   24 GB is tight; if it doesn't fit, this gate moves to the rented box.
7. **Stretch: fused extension → Metal**, with `max_ins/max_del` moved host-side. Gate:
   `extend_bench` / `fused_bench` dumps, 0 mismatches.
8. **tier0 end-to-end on the local Mac** (`--spill-prep`, CPU arm, then Metal arm). This
   flushes macOS bugs (16 K pages, mmap, mimalloc, file APIs) cheaply.
9. **Stage data.** Pull the oracle dumps and tier0 inputs to the Mac. Push the tier1 kit
   (35–50 GiB) to an S3 bucket in the target region, so the box pulls at EBS speed. First
   check whether the Spark has `aws` now and what the uplink does.
10. **Launch runbook.** Write a script for host allocate → instance launch → bootstrap
    (APFS resize, Xcode CLT + Metal toolchain, rustup pin, S3 pull) and dry-run everything
    that does not allocate.

- **A9 (2026-10-21): kit staged in S3.** A private bucket
  (us-east-1, AES256, objects expire after 45 days). `scripts/kit/stage_s3.sh`
  writes `kit/` in the layout `run.sh` expects:
  - `git archive HEAD` plus vendored crates, the tier1 golden manifest and its `@RG DS`;
  - prebuilt Metal kernels;
  - the reference, index, known sites and tier1 FASTQs (the FASTQs stream from the
    Spark over ssh; no AWS credentials go to the Spark);
  - the path-independent 8 GB dense-SA cache.

  Staged: **24.2 GiB at source `90cdb59`**, every object size-checked after upload. Uplink
  from the Mac is ~11 MB/s, so ~40 min. Two kit bugs found and fixed
  on the way:
  - `run.sh` never passed `--read-group-ds`. Every BAM digest would have failed off the
    Spark; the kit's earlier 41/41 was `--digest-only` against a Spark-path run.
  - `--gpu metal` died without the Metal compiler. The compiler is an Xcode-only
    download, so the kit now carries the kernels prebuilt from the same source.
- **A10 (2026-10-21): runbook `scripts/cloud/mac_ec2.sh`.** `setup` has run: a key pair, an ssh-only security group
  for one client IP, and an instance role that can read the bucket and write
  `results/` only. The region had no default VPC, so one was created. `preflight` checks the quota, offering zones, the
  newest macOS 27 AMI, resources and kit. It fails only on the quota and the
  in-progress upload. As of the end of A9, only the quota is not ready.
- **Pre-rent check left:** re-digest the tier1 golden on the Spark when its lock is free.
  The manifest is from 2026-10-09; digest code changes since then are performance-only.

## Phase B log (2026-10-25, `mac-m4max.metal`, us-east-1c)

**Boot (07:29–10:20 PDT, ~3 h lost).** The host was allocated at 07:29
(releasable Mon 07:29; reminder task at 07:15). The macOS 27.0 AMI never reached the
network: the instance reachability check failed and port 22 stayed closed, and a reboot did
not help. In-place root-volume replacement (no terminate, so no host scrub) to 26.7 did
the same. The swap to **15.8 came up 25 min after the swap**. The network config was audited
and is standard. Mac instances are bare metal, so CloudWatch shows no CPU and its EBS
counters are not a boot signal. The hour budget counts from 10:20.

**Box:** Apple M4 Max, 16 cores, 128 GiB, macOS **15.8**, Command Line Tools only, **no Metal
compiler** (the kit's prebuilt kernels are used). `nomics selftest`: Apple9, 40 GPU
cores, 96 GiB GPU working set, 72 GiB largest buffer. **In-place GPU sum: 433 GiB/s**
(Spark 154 GiB/s), exact.

**Fixes found on the live box** (`metal/phase-b`):
- the prebuilt metallib needs `-mmacosx-version-min=15.0`, or macOS 15 refuses it and the
  Metal arm silently runs on the CPU;
- `nohup` cannot detach over a tty-less ssh on macOS;
- Homebrew's `aws` is not on the ssh PATH;
- the APFS container may already fill the disk;
- `cargo vendor -q` leaves `.cargo/config.toml` empty.

**macOS 15.8 vs 27, caveats for every number below:**
- Exactness is unaffected: integer GPU kernels, our own libm and fdlibm ports, and every
  run is digest-gated.
- Timings are partly an OS measurement, because the GPU backend compiler ships with macOS.
  M4 Max (15.8) vs local M4 Pro (27) comparisons mix hardware and OS. Same-box CPU-vs-Metal
  comparisons do not.

| tier1 FASTQ → VCF, 16 threads | total | prep | align | markdup | BQSR | HC | peak RSS | digests |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| CPU arm, run 1 | 706.6 s | 19.5 | 581.8 | 32.8 | 52.2* | 20.3 | 110 GB | PASS (41/41) |
| CPU arm, run 2 | 693.2 s | 36.2 | 582.1 | 33.0 | 20.4 | 21.5 | 110 GB | PASS |
| Metal seeds + rescue (cap64), runs 1 / 2 / 3 | 379.2 / 380.0 / 380.6 s | 19.7–19.9 | 285.0–285.8 | 32.5–32.8 | 20.2–20.4 | 21.4–22.6 | 111–115 GB | PASS ×3 |
| + BWT and seed buffers registered once, two-pass seeding (`832c03b`) | **372.0 s** | 19.8 | 276.8 | 32.6 | 20.1 | 22.6 | 115 GB | PASS |
| + fused GPU extension (`4fa728d`) | 409.9 s | 20.0 | 314.4 | 32.6 | 20.2 | 22.6 | 115 GB | PASS |

\* includes building the known-sites cache. The CPU arm's 132 "GPU seed search failed"
fallbacks came from `run.sh` setting `NOMICS_GPU_SEED` for both arms; fixed (`7a7daa2`).

**Go/no-go (12:00): GO.** The Metal arm is 1.83–1.86× the CPU arm end to end (align
2.04×), exact on every run, with zero GPU fallbacks.

**Profile of the Metal arm** (`NOMICS_TIMING=1`, 132 batches of 662k reads):
- seeding: ~0.55 s of GPU per batch (~70–82 GPU-s), fully overlapped;
- rescue: 46 s;
- CPU align (chaining, `ksw_extend2`, CPU seeding of overflow reads): **7,550 CPU-s**,
  the bottleneck. The GPU is ~40% busy.
- Two-pass seeding (now affordable with the BWT registered) cuts CPU-seeded reads from
  ~19.9k to ~0.8k per batch and CPU align work by 4%.

**Fused GPU extension on Metal (`shim/metal/fused.metal`), exact everywhere:**
- `fused_bench` replays with 0 mismatches: the tier0 dump (200,888 jobs), the band-retry
  fixture (390 jobs, with second bands) and the 30× dump (5,434,925 jobs, 7.84 M calls);
- tier1 run: 482 M lookups, 0 misses, 41/41;
- in isolation on the M4 Max: ~5 M jobs/s, 3.5× the 16-thread CPU; threadgroup size
  (32–256) makes no difference.

In the pipeline it does not pay yet, the same finding as on the Spark. CPU align work
falls 30% (5,049 CPU-s), but the GPU becomes the bottleneck: seeding 110 s + fused
kernels 219 s + rescue 52 s of GPU time in a 309 s batch timeline. The fused kernel time
is ~2× its isolated rate, because seeding, rescue and extension contend for one GPU, and
Metal has no queue priorities (CUDA runs extension at the highest priority).

**More tiles in flight fixes it** (`NOMICS_GPU_EXTEND_LANES`; two-pass seeding, fused
extension, `NOMICS_TIMING=1`; every run PASS):

| GPU extension lanes | 2 | 4 | 5 | 6 | 7 | 8 | 4 + 2 CPU lanes | 2 + 2 CPU (steal) | 2 + 2 CPU (50%) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| align (s) | 314.4 | 250.1 | 239.2 | 237.8 / 240.0 | 245.1 | 260.7 | 249.0 | 262.1 | 272.1 |

The Metal default is now 6 lanes (`4150fb3`; CUDA keeps 2), and two-pass seeding is the
default on every backend again (`3f3cd6d`). End-to-end totals in this sweep vary in
**prep** (20–50 s): runs peak at 122–127 GB of 128, which evicts the 5.2 GB of FASTQs from
the page cache, so prep re-reads them from EBS. There is no swapping. Align is the
like-for-like column.

**Final configuration** (branch head, defaults + `NOMICS_GPU_EXTEND=fused`, no
`NOMICS_TIMING`), two runs:

| tier1, M4 Max (macOS 15.8) | total | prep | align | markdup | BQSR | HC | peak RSS | digests |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| final Metal, run 1 / 2 | **320.0 / 333.0 s** | 29.1 / 43.2 | **216.4 / 215.2** | 32.8 | 20.3 | 21.3–21.4 | 125 GB | PASS / PASS |
| CPU arm (above) | 693.2–706.6 s | 19.5–36.2 | 582 | 33.0 | 20.4 | 21.5 | 110 GB | PASS |

**2.1–2.2× the CPU arm end to end, 2.7× in alignment, exact.** For scale (same pipeline,
other hardware and OS): the Spark runs tier1 in 298 s with fused extension. Memory is the
limit to watch: at this configuration a WES uses ~125 of 128 GB.

**30× whole genome on the Mac, exact** (`scripts/cloud/mac_wgs30.sh`, HG002 SRR14724532,
350 M pairs; whole-genome mode turns on by itself; x86 golden's scatter order and `@RG DS`;
gated by the committed `reference/wgs30-x86-r1.digest.tsv`):

| HG002 30×, FASTQ → VCF | wall | prep | align | markdup | BQSR | HC | peak RSS | digests |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **M4 Max, Metal (seeds + rescue + fused extension), clean** | **4,973 s (1 h 23 m)** | 502 | 3,465 | 284 | 214 | 997 | 105 GB | **41/41** |
| M4 Max, Metal, same, with `NOMICS_TIMING=1` | 5,592 s | 444 | 4,072 | 284 | 215 | 1,009 | 105 GB | 41/41 |
| M4 Max, **CPU arm** (same box, same source) | 10,421 s (2 h 54 m) | 487 | 8,903 | 283 | 212 | 1,010 | 90 GB | 41/41 |
| Spark GB10, fused, for scale (other hardware and OS) | 4,003 s | | | | | | 80 GB | 41/41 |
| sarek on the Spark | 41,044 s | | | | | | | — |

**At 30×, Metal is 2.10× the CPU arm end to end and 2.57× in alignment** on the same box,
exact on every run. Alignment is 70% of the Metal run. HaplotypeCaller (997 s, CPU; PairHMM
is FP64, so no Metal) and prep (502 s) are the next largest pieces.

Inputs (52 GB of FASTQs) are staged at `s3://…/wgs30/`. The VCF, digests and logs are
archived at `s3://…/results/wgs30-{metal-20261025-203620,cpu-20261025-221933,metal-20261026-012431}/`. A 30× run with
`--write-intermediates` takes 212 GB of disk.


The CPU arm's PASS also validates the 2026-10-09 tier1 manifest and the `@RG DS` fix off
the Spark.

### Phase B outcome (2026-10-26)

Final numbers, configurations and the profile are in **`bench/metal.md`**. In short, on an
M4 Max (macOS 15.8), every run is 41/41 against sarek:
- tier1 exome: Metal 320–333 s vs CPU arm 693–722 s (2.1–2.2×);
- 30× whole genome at matched settings: **Metal 4,592 s (1 h 17 m) vs CPU arm 10,161 s
  (2.21×, align 2.81×)**.

At 30× the GPU is saturated (94–100% active). The two levers that pay both move work back
to idle cores: CPU extension lanes and `NOMICS_GPU_RESCUE_SHARE`. Dense SA 2 and `global2`
do not matter there. Cost: one host for 24 h (~$150). About 3 h of it was lost to an AMI
that never booted on this host (macOS 27.0 and 26.7; 15.8 booted).

### Phase B commands

```
scripts/cloud/mac_ec2.sh preflight              # must print READY
scripts/cloud/mac_ec2.sh allocate --yes         # BILLING STARTS; note the release time
scripts/cloud/mac_ec2.sh launch                 # 10–20 min to ssh
scripts/cloud/mac_ec2.sh bootstrap              # APFS, toolchains, ~25 GB kit from S3
scripts/cloud/mac_ec2.sh run cpu                # control arm (twice); then:
scripts/cloud/mac_ec2.sh run metal              # Metal seeds + rescue (twice)
scripts/cloud/mac_ec2.sh status | fetch         # progress; digests and timings back here
scripts/cloud/mac_ec2.sh release --yes          # after 24 h; terminates, waits out scrubbing
```

**Ready to rent when:** A1–A3, A5 and A9–A10 are done, and A6 is gated or blocked only on
memory.

## Phase B: the 24 h on `mac-m4max.metal`

Times are hours from host allocation.

| Hours | Work | Output |
|---|---|---|
| 0–1 | Launch instance (10–20 min boot). Bootstrap: APFS resize, toolchains, S3 pull, `cargo build --release --features metal`. Start a timer for host release at 23:30. | box ready |
| 1–3 | **CPU control arm**: tier1 FASTQ → VCF, in memory, twice; kit digests. The first run builds the dense-SA and known-sites caches. | "sarek-identical exome on a Mac in N min" |
| 3–4 | On M4 Max: Metal smoke, bandwidth (peak ~546 GB/s), dependent-lookup bench, buffer limits; seed gate if A6 was memory-blocked. | numbers for the doc |
| 4–12 | **Metal arm**: seeds + rescue in `nomics run`. Gate every tier1 run by digest. Profile with xctrace; tune threadgroup sizes, queue layout and the CPU thread count. macOS has no core pinning; P/E-core balance goes through QoS. | exact Metal runs, timings |
| **12** | **Go/no-go.** Does the Metal arm beat the Mac CPU arm on tier1 wall time? | decision |
| 12–20 | Go: tune further, stretch to fused extension, run the 1/12 genome (wgs20) in whole-genome mode. No-go: profile why, and spend the time on CPU-arm wins and the wgs20 run. | |
| 20–23 | Final runs: each arm twice, identical digests. Write `bench/metal.md`. Update the status page (Apple Silicon row). | results |
| 23–24 | Copy results, logs and traces to S3, then pull them to the NAS. Release the host. | host released |

### Gotchas that cost hours on the clock

- **Never stop or terminate the instance to "reset" it.** That triggers host scrubbing
  (hours, and still billed). Reboot only.
- M4 Max has **no instance store**. Size the EBS root at 400–500 GB and raise gp3
  throughput (default 125 MB/s is unusable). Run `diskutil apfs resizeContainer` on first
  boot; macOS AMIs do not grow APFS.
- Check that the AMI's Xcode or CLT includes the Metal toolchain
  (`xcodebuild -downloadComponent MetalToolchain` on macOS 26).
- Long jobs use `nohup` + `timeout`; foreground ssh sessions die.

## Success criteria

1. Tier1 on macOS, CPU arm: 41/41 digests identical, wall time and peak RSS recorded.
2. At least seeds and mate rescue on Metal, each exact on its oracle gate and in a full
   tier1 run.
3. The Metal arm timed against the Mac CPU arm on the same box. It ships as the default
   only if it is faster.
4. Every residual or limit (buffer size, FP64, no priorities) is written down in
   `bench/metal.md`.
