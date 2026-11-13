# Design: fused left+right GPU extension with band retries (option 1, 2026-10-18)

Scope doc: `docs/strategy/gpu-alignment-scope.md`, option 1. This is the first half: one
GPU job per seed instead of two whole-batch launches. The three-deep pipeline (seed GPU →
chain+plan CPU → extend GPU → replay+pair CPU) waits for the Spark measurement of this half.

## What BWA does per seed (`chain2aln`, `crates/nomics-align/src/extend.rs`)

```
if qbeg > 0:                                   # left, on reversed prefix and window
    prev = -1                                  # a.score starts at -1
    for i in 0..MAX_BAND_TRY (2):
        w = 100 << i
        x = ksw_extend2(qL, tL, w, h0 = len*A, end_bonus = PEN_CLIP5)
        if x.score == prev or x.max_off < (w>>1)+(w>>2): break
        prev = x.score
    sc0 = x.score
else:
    sc0 = len*A
if qbeg+len != l_query:                        # right
    prev = sc0
    for i in 0..2:
        w = 100 << i
        x = ksw_extend2(qR, tR, w, h0 = sc0, end_bonus = PEN_CLIP3)
        if x.score == prev or x.max_off < (w>>1)+(w>>2): break
        prev = x.score
```

Each call is looked up in `Extended` by `ExtKey {chain, seed, right, w, h0}` and computed on
the CPU on a miss. Everything except containment (which seeds get extended at all) is a pure
function of the seed and the chain's reference window. So a GPU job can run a seed's whole
sequence speculatively. Contained seeds waste a job, but outputs are unchanged, because the
replay looks up only the calls BWA makes.

Today's path launches left (w=100 only), waits, plans right with h0 from left's first-band
score, launches again, then replays. A seed whose left retried gets a right key whose h0 is
wrong for BWA's call: the lookup misses and the CPU redoes it. That is exact but costs a launch
and a barrier per batch.

## Fused job

The GPU runs the loop above verbatim, up to four `ksw_extend2` calls per seed, and returns
every call it made. The replay is unchanged and still matches the full `ExtKey`, including `w`
and `h0`. What the key guards is control flow: a call BWA did not make (a band-1 retry it
skipped, a right extension computed from the wrong h0) can never be consumed. It does not
guard arithmetic. A miscomputed result under the right key is a hit and gets used.
Computational exactness comes only from the oracle gates (3 and 4 below), as for every GPU
kernel so far.

Seeds go to the CPU (no job) when the left query is over `EXTEND_MAX_QUERY` (256), because the
right's h0 depends on it. When only the right query is too long, the job runs the left only and
the right misses to the CPU.

### C ABI (`shim/nomics_gpu.h`)

```c
// Scoring shared by every job of a launch (BWA mem defaults today).
typedef struct {
    int32_t m;                 // alphabet size, m*m <= 64
    int32_t o_del, e_del, o_ins, e_ins;
    int32_t w;                 // first band, 100; band i is w << i
    int32_t max_band_try;      // 2
    int32_t end_bonus_left;    // PEN_CLIP5
    int32_t end_bonus_right;   // PEN_CLIP3
    int32_t zdrop;
    int8_t  mat[64];           // m*m used
} nomics_fused_params;         // 104 bytes

typedef struct {
    uint64_t lq_off, lt_off;   // left query (reversed prefix) and target (reversed window)
    uint64_t rq_off, rt_off;   // right query (suffix) and target
    uint32_t lq_len, lt_len;   // lq_len == 0: no left (qbeg == 0)
    uint32_t rq_len, rt_len;   // rq_len == 0: no right on the GPU
    int32_t  h0;               // seed len * A
    int32_t  pad;
} nomics_fused_job;            // 56 bytes

typedef struct {
    nomics_extend_result left[2], right[2];
    int32_t n_left, n_right;   // calls made per side, 0..2
    int32_t right_h0;          // sc0
    int32_t pad;
} nomics_fused_result;         // 112 bytes

// Queries up to 256 bases per side; lq_len + rq_len > 0. Same stream and
// synchronisation contract as nomics_ksw_extend2. `order`: n entries of
// caller-owned workspace (kernel v3's length-binned job order).
int nomics_fused_extend(const uint8_t* queries, const uint8_t* targets,
                        const nomics_fused_params* params,
                        const nomics_fused_job* jobs, uint32_t n, uint32_t* order,
                        nomics_fused_result* out, float* out_ms);
```

Call i on a side used band `w << i`; left calls used h0 = `job.h0`, and right calls used
`right_h0`. The host turns each result into `ExtKey`s, so the replay needs no new logic.

### Kernel (`shim/fused.cu`)

- Reuse `extend_group<W, C>` from `extend.cu`, refactored into a device function that takes
  pointers, lengths, w, h0 and end_bonus, and returns the result in registers. Every lane
  already holds the reduced scalars. The refactor must leave `nomics_ksw_extend2`'s behaviour
  bit-identical; the existing `extend_bench` GPU gate rechecks it.
- One group per job. The group width (8/16/32) and the column count C are chosen from
  `max(lq_len, rq_len)`, with the same launch filters and per-warp uniform C as today.
- A group runs left band 0, at most band 1, then right band 0 and at most band 1, with the
  break test above evaluated identically on every lane.
- Retries are rare: 0 of 1.96M dumped calls (exome and wgs20) used w=200. Divergence from them
  is not a concern; correctness on them is, so they get synthetic tests.

## Host (`nomics-align`, `nomics-gpu`)

- `plan_fused` in `extend.rs` replaces `plan_left` + `plan_right`. It makes one pass per read
  over chains and seeds, writes bases into the arena, and computes the chain window once per
  chain.
- `FusedExtend` engine, with two implementations:
  - `nomics_gpu::GpuFusedExtend`: HostBuf arenas, the FFI call, and a stub without `cuda`;
  - a CPU stand-in in `nomics-align` that runs the same loop with `ksw::extend2`, in parallel.
  The stand-in is the Mac gate and the reference for the kernel.
- `align1_gpu_extend_tile` does one fused launch per tile: chains → plan_fused → engine →
  replay through `align1_extend`, with `Extended` built from the fused results.
- Lookup by direct per-seed index instead of a linear scan of the read's keys, still comparing
  the full `ExtKey`.
- Selection: `NOMICS_GPU_EXTEND=fused`. `NOMICS_GPU_EXTEND=1` keeps the current two-pass path
  until the Spark measurement decides. `NOMICS_EXTEND_ENGINE=cpu` forces the stand-in, which is
  the only engine without `cuda`.
- `NOMICS_DUMP_FUSED=<dir>` writes `fused.bin`: params, every job with its bases, and the
  results. A `fused_bench` example replays a dump through the stand-in and, with `cuda` and
  `EXTEND_GPU=1`, through the kernel, reporting mismatches field by field.

## Gates

**Mac**

1. Unit tests: stand-in vs a scalar transcription of `chain2aln`'s two loops on synthetic
   seeds, including:
   - forced retries on both sides (indels of 80+ bp, so max_off ≥ 75);
   - left only, right only, both;
   - a right query over 256;
   - tlen 0;
   - qlen 1 and 256.
2. End to end: tier0 chunks (`~/omics-fixtures/sarek-tier0/r1/work/*/0001.*fastp.fastq.gz`,
   index `~/omics-fixtures/gvcf-tier0/ref/BWAIndex`) through `align_chunk` with the fused
   stand-in. The output must be byte-identical to the CPU path's on the same chunk, and match
   sarek's `mapped/HG002_WES/HG002_WES.0001.bam` under the existing chunk gate. Report hits,
   misses, CPU fallbacks and retries.
3. `fused.bin` from gate 2, plus the synthetic retry cases, replayed through the stand-in: 0
   mismatches.

**Spark (after sarek wgs30 r1 finishes)**

4. `fused_bench` on gate 3's dump with the kernel: 0 mismatches. `extend_bench` with
   `EXTEND_GPU=1` on both dumps: still 0 (refactor guard).
5. Tier1 18/18, the wgs20 chunk gate, and the 30× chunk (`scripts/wgs30_align_gate.sh`), all
   with `NOMICS_GPU_EXTEND=fused`.
6. Measurement as in the scope doc: chunk wall time vs CPU and vs the two-pass path, GPU
   timeline, huge pages ≥ 90%, a timing-off control, and peak memory (≤ 90 GB at 30×).

## Memory

Per job: 56 B job + 112 B result + bases (~2 × (qlen + window)). A 331k-pair batch has ~6M
jobs, so ~1 GB of results and jobs, on the order of today's two-pass arenas, which are freed
per tile. Gate 6 checks the 30× peak.

## Results on the Spark (2026-10-18, commits c7cacdf + 363cc97 + fixture fix)

Run dir: `~/omics/runs/fused-gate-20261018-2321` (`summary.txt`, per-mode `run.log`,
`smi-*.csv`).

### Exactness: every gate passes

| Gate | Result |
|---|---|
| `extend_bench EXTEND_GPU=1`, exome dump (refactor guard) | 0 mismatches of 1,202,012 |
| `extend_bench EXTEND_GPU=1`, wgs20 dump | 0 mismatches of 754,265 |
| `fused_bench`, retry fixture (both sides retrying) | 0 mismatching jobs |
| `fused_bench`, tier0 chunk 0001 (200,888 jobs) | 0 mismatching jobs |
| tier1 end to end, fused | 18/18, 506 s |
| wgs20 end to end, fused | 18/18, 265 s |
| 30× chunk 1, all four modes | byte-identical mapped BAM |

30× chunk 1 with timing on:
- 622,676,225 lookups, 0 misses, 0 CPU fallbacks;
- 595,802,750 fused jobs;
- band retries made by BWA's control flow: 3; made by the engine: 6. They are real at
  whole-genome scale, and the GPU handled them exactly.

### Speed: exact but slower

| 30× chunk 1 | Align phase | Wall | Mean GPU |
|---|---:|---:|---:|
| CPU extension (default) | **379 s** | 426 s | 47% |
| Two-pass GPU extension | 392 s | 429 s | 63% |
| Fused GPU extension | 423 s | 461 s | 69% |
| Fused, timing on | 441 s | 478 s | 69% |

From the timing run, per chunk:
- 970 tiles;
- extension "kernel" time 759 s, from the host-side timer. This is not GPU time: it includes
  waiting in the queue behind other kernels. The nsys trace below measures 143 s;
- seed kernel 191 s;
- rescue 140 s;
- CPU busy 3,090 core-s, about 7 of 20 cores over the align phase. The CPU path keeps ~18.7
  busy.

What it means:
- The offload worked as designed: the CPU's per-read work fell from ~130 µs to ~45 µs
  (chaining 34 µs, chain2aln replay 9 µs).
- The GPU is now the bottleneck. *Corrected 2026-10-19 from the nsys trace:* the extension
  kernel runs close to its isolated rate (143 GPU-s in the pipeline vs ~120 s at the tier0
  bench's 5M jobs/s). It is slow end to end because it adds 143 s to a GPU that already spends
  ~200 s per chunk on seeding and rescue, and the kernels run one at a time.
- Each tile still waits on its own launch: chains + plan + replay is ~0.3 s of CPU per tile,
  against 0.6–2.1 s of kernel.

### Decision

Keep `NOMICS_GPU_EXTEND=fused` opt-in. CPU extension stays the default.

### Next, in order

1. **Split extension between CPU and GPU.** About 13 CPU cores sit idle on the fused path and
   about half the GPU sits idle on the CPU path. Sending a fraction of tiles (or reads) to the
   CPU path, fixed or work-stealing, is exact by construction and cheap to try. It should beat
   both 379 s and 423 s.
2. **nsys timeline** of the fused chunk: split kernel time into execution vs queueing, and see
   how much extension slows seeding.
3. **Kernel v3** (several jobs per warp for short bands, fewer row syncs). This is the option-2
   work the scope ordered after option 1; it pays once the split shows how much GPU time
   extension can take.

## CPU/GPU split (2026-10-19, commit 19f4086)

CPU lanes run next to the GPU lanes and align whole tiles as the default path does, so they
are exact by construction. There are two modes:
- `NOMICS_GPU_EXTEND_CPU_LANES=N`: work stealing from one queue;
- adding `NOMICS_GPU_EXTEND_CPU_SHARE=f`: a fixed fraction of the tiles goes to the CPU lanes.

Sweep: `scripts/extend_split_sweep.sh` on 30× chunk 1, timing off. Run dir:
`~/omics/runs/extend-split-20261019`. Every config gave a byte-identical BAM.

| Config | Align phase | Mean GPU |
|---|---:|---:|
| CPU only (default) | 370.1 s | 46% |
| fused | 431.2 s | 68% |
| fused + 1 CPU lane | 397.0 s | 63% |
| **fused + 2 CPU lanes** | **369.9 s** | 65% |
| fused + 4 CPU lanes | 415.6 s | 63% |
| fused + 2 CPU lanes, 16k-read tiles | 392.1 s | 62% |
| fused + 4 CPU lanes, 16k-read tiles | 379.9 s | 58% |
| fused + 2 CPU lanes, share 0.3 | 398.4 s | 66% |
| fused + 2 CPU lanes, share 0.5 | 406.6 s | 67% |

The best split ties the CPU path (369.9 vs 370.1 s, within run-to-run noise of ~10 s). It
does not beat it. In the fused-only run, ~13 cores were idle and the GPU sat at ~68%. Even so,
giving work to both sides doesn't add their throughputs. My first guess was shared memory
bandwidth; the next section measured it and found otherwise.

### Decision

Keep extension on the CPU by default. Fused GPU extension and the split stay opt-in. Both are
exact and gated, and they may pay on a discrete GPU with its own memory.

### Before more GPU-extension work

Measure the memory contention instead of assuming it. One nsys/ncu pass on the fused chunk
should give:
- DRAM throughput of the extension and seed kernels;
- the CPU's memory bandwidth from perf uncore counters, with and without the kernels running.

If the box is near its bandwidth ceiling, kernel v3 (fewer row syncs) won't help either. Work
that cuts bytes moved would: register-resident target windows, compact 2-bit bases, and fewer
speculative jobs. Measure how many fused results the replay never uses (jobs for contained seeds) first.

## Is the GB10 bandwidth-bound? (2026-10-19)

No. The fused path is bound by **GPU time**, with kernels running one at a time. What
interferes is the CPU's **memory latency** under heavy GPU traffic, not exhausted bandwidth.
Run dir: `~/omics/runs/bw-profile-20261019`. GPU hardware counters are admin-only here
(`RmProfilingAdminOnly=1`), so ncu and nsys GPU metrics were unavailable. Instead:
- an nsys CUDA trace for GPU occupancy;
- per-core PMU counters, calibrated against a triad, for CPU traffic;
- direct interference tests.

**Ceilings** (`~/omics/bench/bw`: an OpenMP triad and a CUDA copy):

| | Alone | Both at once |
|---|---:|---:|
| CPU triad, 20 threads | 115 GB/s | 33 GB/s |
| GPU copy | 218 GB/s | 196 GB/s |

The combined ceiling is ~230 GB/s, and the GPU wins the arbitration.

**GPU timeline of the fused 30× chunk** (nsys, 396 s span):

| Kernel | GPU time | Launches | Mean |
|---|---:|---:|---:|
| `collect_kernel` (SMEM seeding) | 150 s | 89 | 1.68 s |
| `fused_kernel` (extension) | 143 s | 2,910 | 49 ms |
| `ksw_wave_kernel` (mate rescue) | 52 s | 89 | 0.58 s |

- The GPU is busy 331 s, 84% of the span.
- Only one kernel runs for 318 of those seconds; two overlap for 13 s.
- Seeding plus rescue (~200 s) plus extension (143 s) is about the whole span, so the GPU is
  the critical path.

**CPU traffic** (`l3d_cache_refill` / `bus_access` over the whole gate process). These are
estimates: no DRAM PMU is exposed, so counts were converted with the triad's calibration
(~158 B per refill, ~15 B per bus access):
- fused chunk: ~5.7–6.3 TB in 450 s, **~13–14 GB/s**;
- CPU-path chunk: ~3.7–4.1 TB in 408 s, **~9–10 GB/s**.

Both are about a tenth of the CPU's 115 GB/s ceiling. The aligner's memory access is random
(FM index, SA, reference windows), so it is latency-bound, not bandwidth-bound.

**Interference:**
- **The CPU path is latency-sensitive to GPU traffic.** The CPU-path chunk takes 370 s alone
  and 593 s (+60%) while a GPU copy loop saturates memory. That hog is far heavier than our
  kernels; it bounds the effect and agrees with the 2026-10-10 observation that CPU work slows
  by about a third during seeding.
- **The extension kernel is compute-bound.** On the tier0 dump it takes 40 ms alone and 41 ms
  under a 118 GB/s CPU triad, and the triad didn't slow either.

**What this changes:**
- Kernel v3 is back on the table. Extension is compute-bound and runs near its isolated rate,
  so fewer row syncs and several jobs per warp translate straight into GPU time.
- Seeding is the largest GPU consumer (150 s per chunk, 44% of GPU time). It is the one kernel
  likely to put heavy random traffic on memory. A faster seed kernel frees GPU time and
  probably lowers the CPU's latency penalty.
- The split tied because GPU and CPU time were both near full: each tile moved to the CPU
  costs CPU time, which latency interference inflates. It should start paying once kernel v3
  and seeding free GPU time.
- The host "kernel" timers measure enqueue to completion, so queue wait is included. Event
  timing on a per-call stream cannot separate the two; use nsys for GPU time.

Unverified without GPU counters: the seed kernel's actual DRAM throughput. Enabling profiling
for non-admin users (`NVreg_RestrictProfilingToAdminUsers=0`, then a reboot), or running ncu
under sudo, would settle it. That is a system setting, for the owner to change.

## Kernel v3 (2026-10-19)

Run dir: `~/omics/runs/kernel-v3-20261019` (sweep logs, `nsys-v3.*`, `split/`,
`pipe-*/`). Code: `shim/fused.cu` (launcher and bins), `shim/extend_v3.cuh` (the call body).
`extend_group` and `nomics_ksw_extend2` are untouched.

### What changed

- **Length bins.** The launcher counting-sorts job indices by `max(lq_len, rq_len)` into bins,
  each one launch of `fused_bin_kernel<W, C>` with a fixed group width and column count
  (8×1..4 up to 32, 16×3..5 up to 80, 32×3..8 up to 256). Within a bin, jobs run in
  descending target rows. Results are written to `out[index]`, so they stay in input order.
  The order lives in a new caller-owned `order` workspace (`HostBuf`, 4 bytes per job): an
  ABI change.
- **A new body, `extend_group_v3`.** It does the same arithmetic as `extend_group`, with
  fewer instructions per cell:
  - one profile word per column (five 6-bit fields);
  - branch-free columns;
  - F carried as `f = running - j*e_ins`, which removes the per-column multiplies.
  The header comment argues the equivalence. It needs `m <= 5`, matrix entries in
  [-32, 31] and gap penalties in range. The launcher refuses other scoring, and the caller
  then extends on the CPU.
- **Registers.** Scoring is copied into scalars. v2 put the 104-byte params struct on the
  stack of every thread. The job is read at each call site instead of held in registers.
  `__maxnreg__(80)` applies for C ≤ 5. Blocks are 32 threads.
- **Highest-priority stream** for the extension launches (`NOMICS_FUSED_PRIORITY=0` turns
  it off). This changes scheduling only. It is the change that moved the align phase; see
  below.

### Isolated kernel time (fused_bench, steady reps, ms)

| Step | tier0 (200,888 jobs) | WGS 30× dump (5.43M jobs) |
|---|---:|---:|
| v2 | 40.0 | 785 |
| bins + sort + per-(W,C) instantiation + params off the stack | 28.7 | 584 |
| + `extend_group_v3` body | 22.4 | 479 |
| + 32-thread blocks, job by reference | 21.7 | 452 |
| + `__maxnreg__(80)` for C ≤ 5 | 21.5 | 437 |
| + 16×5 bin for 65–80 | 21.4 | 428 |
| final (packed reduce removed) | **21.4** | **426** |

The final kernel is 47% faster on tier0 and 46% faster on the WGS dump.

Measured and rejected (WGS dump, ms, against the step it modified):

| Variant | Before | After |
|---|---:|---:|
| Narrower groups for 65–128 (16×5..8) | 584 | 660 |
| Narrower groups for 65–160 (16×5..10) | 584 | 681 |
| Narrower groups for 65–160 (8×10..16) | 584 | 846 |
| W=4 for ≤32 | 584 | 589 |
| W=8 for 33–64 | 584 | 588 |
| W=32 for ≤64 | 437 | 449 |
| 256-thread blocks | 584 | 693 |
| No sort within bins | 437 | 446 |
| Band bit masks instead of compares | 452 | 513 |
| One packed score+column reduce (`row_max << 9 \| mj + 1`, bounded) | 424 | 428 |

Per-cell instruction count dominates, not the collectives: halving the group width doubles
each lane's columns and loses.

### Exactness

All gates show 0 mismatches:
- `fused_bench` on the retry fixture (390 jobs, re-blessed);
- tier0;
- the new WGS dump.

`extend_bench` on both ksw dumps also shows 0 mismatches.

The retry fixture gained 135 synthetic jobs (255 → 390):
- every multiple of 8 and the length after it, up to 256 (bin edges);
- the long side left or right, with the other side 1 base or random;
- large h0 (2^22 − qlen) on each side.

New WGS dump: `~/omics/bench/fused/wgs30-0001.fused.bin`, 1.93 GB, 5,434,925 jobs. It was
generated from 30× chunk 1 through the CPU stand-in:
- `NOMICS_DUMP_FUSED_EVERY=110` writes every 110th launch, so the sample spans the whole
  chunk;
- `NOMICS_DUMP_FUSED_MAX_JOBS=5200000` caps the total.

That run passed the chunk gate the same way (records identical, `@PG` differs).

### In the pipeline (30× chunk 1, align-only, timing off)

Run-to-run noise is about ±10 s, so the table lists every run.

| Config | Align phase |
|---|---:|
| CPU extension | 378.6, 365.0 s |
| fused v2 (earlier session) | 431.2 s |
| fused v3, default-priority stream | 425.1 s |
| fused v3, default priority, 128-thread blocks | 430.3 s |
| fused v3, default priority, 3 GPU lanes | 437.2 s |
| fused v3, default priority, 4 GPU lanes | 437.4 s |
| fused v3, default priority, + 1 CPU lane | 401.8 s |
| fused v3, default priority, + 2 CPU lanes | 366.9 s |
| fused v2 kernel + high-priority stream | 412.2, 399.9 s |
| **fused v3 + high-priority stream (final)** | **357.9, 369.1, 375.9 s** (mean 367.6) |
| fused v3 + high priority, + 2 CPU lanes | 362.5, 392.1 s |

Every run passed the chunk gate as "normalized": all 58,621,557 records and refs identical, and only the `@PG` header line differs. The CPU control shows the same.

On its own, the faster kernel did not move the align phase (425 vs 431 s). nsys showed why:

| nsys, 30× chunk 1 | v2 | v3, default priority | v3, high priority |
|---|---:|---:|---:|
| Align phase under the profiler | 419.6 s | 449.0 s | 378.5 s |
| Extension kernels, union | 138.9 s | 99.1 s | 151.7 s* |
| `collect_kernel` (seeding), union | 149.9 s | 174.1 s | 194.3 s* |
| Seeding × extension overlap | 7.9 s | 36.0 s | 95.0 s |
| GPU busy (union) / span | 330.9 / 396 s | 286.4 / 418.6 s | 301.9 / 352.4 s (86%) |
| Extension launch → start, p50 / p90 | 27 ms / 2.16 s | 39 ms / 2.22 s | 53 ms / 0.48 s |

\* Co-running kernels stretch each other's durations, so per-kernel GPU-seconds stop being
additive once they overlap. Read the kernel's own cost from the isolated table (−46%).

- With the default priority, a GPU lane waits on its own tile's launch. That launch queued
  behind the pending blocks of a ~1.8 s seeding kernel on the legacy stream (p90 2.2 s), so
  the GPU sat 30% idle and the time the kernel saved went to waiting.
- The highest-priority stream lets extension blocks dispatch as soon as seeding blocks
  retire. The p90 wait falls to 0.48 s and the span shrinks by 66 s.
- The two changes need each other:

  | Kernel | Default priority | High priority |
  |---|---:|---:|
  | v2 | ~431 s | ~406 s |
  | v3 | ~425 s | ~368 s |

### Review and hardening

A GPT-6 Astra adversarial review of `extend_v3.cuh` and `fused.cu` found the following
correct on BWA's inputs:
- the arithmetic, including profile sign extension, the F rewrite, ties, z-drop, trimming and
  `max_off`;
- input-order results through the sort;
- warp collectives;
- the CPU fallback on a refused launch.

It found four gaps for inputs BWA never produces, now closed in the launcher:

| Gap | Fix |
|---|---|
| Score headroom: with h0 near 2^31, the F rewrite's translated maxima overflow differently from the reference | h0 must be in [1, 2^24]; gap penalties and end bonuses ≤ 2^16; w ≤ 2^20. Every score, F term and band bound stays far below 2^31 |
| A query symbol ≥ m reads past the matrix (and the 6-bit packing differs from the reference there) | the profile is all zeros for such a symbol: memory-safe, but outside the exactness claim, as in BWA itself (its contract is 0..m-1) |
| An end bonus near INT_MAX makes the band bound's double-to-int conversion undefined (shared with `extend_group`) | end bonuses limited to [0, 2^16] |
| `NOMICS_FUSED_BINS` checked capacity on longs, then narrowed them | values are range-checked before narrowing |

### Gates on master (integrated build, 2026-10-19)

`scripts/fused_spark_gate.sh`, run dir `~/omics/runs/fused-gate-v3-20261019`:
- kernel gates, 0 mismatches each: the extension dumps (1,202,012 and 754,265 calls), the
  retry fixture (390 jobs), tier0 (200,888) and the whole-genome dump (5,434,925);
- tier1 18/18 and wgs20 18/18, both fused;
- 30× chunk 1, byte-identical in all four modes: CPU 374.9 s, two-pass 389.9 s, fused 372.8 s,
  fused with timing on 376.8 s;
- 622,676,225 lookups, 0 misses.

The hardened build was re-gated separately (`~/omics/runs/fused-harden-20261019`). Every
kernel gate again had 0 mismatches, and kernel times were unchanged (tier0 21.25 ms, WGS dump
427.7 ms). Its 30× chunk passed but took 395.8 s, with the index fully on 2 MiB pages.

All v3 + priority chunk runs so far: 357.9, 369.1, 375.9, 372.8 and 395.8 s, a mean of 374 s.
CPU path: 378.6, 365.0 and 374.9 s, a mean of 373 s. Still a tie, with more run-to-run spread
than the ±10 s assumed so far.

### Decision

The fused path now ties the CPU path (mean 367.6 vs 371.8 s, within noise), against 50–60 s
slower before. CPU extension stays the default for now. The fused path (v3 + priority) is
GPU-bound again (86% busy) and leaves ~13 CPU cores free for the next pipeline stage.

### Next

- Seeding is the largest GPU consumer, and it now co-runs with extension. A faster or
  shorter-launch seeding kernel pays directly.
- Double-buffer each GPU lane, so it plans tile k+1 while tile k's kernel runs.
- Rerun the CPU/GPU split once seeding changes. With priority on, adding 2 CPU lanes did not
  help (362.5 and 392.1 s).

## Faster seed kernel (2026-10-19)

This isn't extension work, but it was the prerequisite for extension to pay.
`collect_kernel` (BWA's `mem_collect_intv`, one thread per read) was not bound by FM-index
traffic. It touched ~916M host lines per batch at ~38 GB/s, against a measured random-line
ceiling of ~160 GB/s. The bound was the per-thread interval lists: 6.4 KB per thread in local
memory, interleaved word by word across the warp, so one 32-byte interval cost 8 sectors.

The new `collect_slot_kernel` keeps each thread's three lists contiguous in a huge-page host
scratch (302 MB). Each block claims a slot from a device-global pool. Intervals are packed to
16 bytes, and `put()` refuses any value out of range, sending that read to the CPU. The L1
carveout is 0, and `GpuSeeds` buffers are pooled.

| | Before | After |
|---|---:|---:|
| Seed kernel, 662k-read 30× batch | ~1,540 ms | ~590 ms (2.6×) |
| nsys seeding GPU-s, fused 30× chunk | 194 s | 86 s |
| 30× chunk 1 align, fused | ~374 s | 293–297 s |
| 30× chunk 1 align, CPU extension (default) | ~373 s | 321–333 s |
| tier1 FASTQ → VCF, default | 460 s (bake-off) | 379 s |
| tier1 FASTQ → VCF, fused | — | 392 s |
| wgs20 FASTQ → VCF, default / fused | — | 265 / 278 s |

A GPT-6 Astra review found no issue on the single-device path. It flagged four latent hazards,
now fixed:
- slot regions now use the plain 32-byte stride for packed and plain launches alike, so a
  slot owns the same bytes either way;
- a failed second launch in two-pass mode drains the device before returning;
- the slot-count cache is atomic for C ABI callers;
- read lengths are compared before narrowing to int, which fixes a hole inherited from the
  old kernel.

Gates on the final code (`~/omics/runs/seed-master2-20261019`):
- `seed_gate`: 10M reads (30× chunk 1, tier1, wgs20) identical to the CPU, in cap-64 and
  two-pass modes. The overflow sets are identical to the old kernel's, and packed vs plain
  lists agree.
- Every kernel dump: 0 mismatches.
- tier1 and wgs20 18/18, on both the default and fused paths.
- 30× chunk 1 byte-identical in all four modes.

### Decision

The CPU extension path stays the default for now. Fused is ahead on the 30× chunk (~295 vs
~327 s), but behind on tier1 (392 vs 379 s) and wgs20 (278 vs 265 s). Next:
- a full 30× run, to decide the default per input size;
- rescue (`ksw_wave_kernel`, 81 GPU-s fused) is now the next GPU consumer;
- seeding still queues behind other work on the legacy stream (p90 0.9 s).

## Full 30× end to end (2026-10-20, commit 5137125)

Two full runs of `scripts/full_exome.sh wgs30`, using the bake-off's procedure: page caches
evicted and a 60 s settle before each run, one run at a time under the lock. Run dirs:
`~/omics/runs/nomics-full30-{default,fused}-20261020`, summary
`~/omics/runs/full30-20261020.txt`. Both gated **41/41 identical**.

| Stage | Bake-off (2026-10-17) | Default (CPU extension) | Fused GPU extension |
|---|---:|---:|---:|
| prep (overlaps alignment) | 521 s | 460 s | 330 s |
| align | 4,681 s | 4,117 s | **3,578 s** |
| markdup | 294 s | 297 s | 294 s |
| BQSR | 392 s | 393 s | 392 s |
| HaplotypeCaller | 966 s | 956 s | 955 s |
| **FASTQ → VCF** | **6,349 s** (1 h 46 min) | **5,781 s** (1 h 36 min), −9% | **5,240 s** (1 h 27 min), −17% |
| Peak RSS | 74 GB | 78.8 GB | 80.1 GB |

What it means:
- Alignment (3,578 s) is still most of the run.
- Against sarek's 41,044 s, nomics is now 7.8× faster.
- Peak memory stays inside the 90 GB budget.
