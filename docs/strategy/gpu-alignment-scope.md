# Scope: moving alignment work to the GPU, exactly (2026-10-18)

> **Status 2026-10-19.** Option 1 is built and exact on every gate
> (`docs/strategy/gpu-extend-fused-design.md`). It does not yet pay: on the 30× chunk, align
> takes 413–431 s with the fused GPU extension vs 370–379 s on the CPU, and a CPU/GPU split only
> ties. The cause is GPU time. In nsys, kernels run one at a time, and the extension kernel
> (143 GPU-s per chunk) adds to seeding (150 s) and rescue (52 s). The box is not
> bandwidth-bound.
>
> The −1,400 s estimate below assumed GPU headroom that seeding and rescue already use. That
> makes option 2 (faster extension and seeding kernels) the prerequisite rather than the
> follow-up. Extension kernel v3 is in progress.

## Where the gap is (Spark bake-off, HG002 30×, two runs each, identical outputs every time)

| Stage | nomics |
|---|---:|
| Alignment + sort | 4,681 s (prep 521 s overlapped inside it) |
| Duplicate marking + BQSR | 686 s |
| HaplotypeCaller | 966 s |
| **Total** | **6,349 s** |

Resource use during alignment (system sampler, 5 s, first 4,700 s):

| | nomics |
|---|---:|
| Busy CPU cores (of 20) | 18.7 |
| GPU utilization | 53% |

nomics alignment is **CPU-bound with the GPU half idle**. The lever is moving CPU work to the GPU without changing a
single output.

## Today's split, per read at 30×

From the 2026-10-13/15 profiles; the S-08 loop change since then makes extension ~8% cheaper.

- **GPU:** SMEM seeding (`seed_collect`, overlapped one batch ahead) and mate rescue
  (`ksw_align2`, warp-wavefront).
- **CPU, per record:**
  - `ksw_extend2` ~85 µs, about 65% of align1;
  - chaining ~28 µs (SA lookups at dense-SA interval 2, chain building);
  - the rest of `chain2aln` ~15 µs;
  - pairing and output.
- Model check: ~130 µs × 700M reads ÷ 20 threads ≈ 4,550 s, against 4,681 s measured.

To bring alignment to ~2,700 s, the CPU work per read has to drop by about 42%
(~55 µs per read), while the GPU absorbs it inside its ~2,200 GPU-seconds of headroom.

## Options, all exact by construction

GPU results must equal the CPU port's call for call. This is the discipline already used for
seeding, rescue and the opt-in GPU extension:
- GPU kernels compute keyed results;
- the CPU replays BWA's control flow, looking results up;
- any miss falls back to the CPU;
- oracle dumps gate every kernel (`extend_bench`, `ksw_bench`).

### 1. Extension offload, done right (recommended first)

The existing opt-in path (`NOMICS_GPU_EXTEND=1`) is exact (403M lookups, 0 misses) but
neutral: the chunk takes 417 s on GPU, 417 s on CPU. It doesn't pay because:
- whole-batch barriers: left launch, then right launch (right h0 = left score), then CPU replay;
- launches queue behind the other batch's seeding and rescue;
- CPU planning and packing (~0.3 s per side per batch).

Changes:
- Do each seed's left and right extensions in one kernel, the same warp running left and then
  right with h0 from the left result. This is exact, and removes a barrier and a launch.
- Restructure into a three-deep pipeline: seed (GPU) → chain + plan (CPU) → extend (GPU) →
  replay + pair (CPU). The GPU then always has queued work and the CPU never waits on it.
- Kernel v3: several calls per warp for short bands, fewer per-row sync points.
- Plan straight into pinned HostBufs (most of this already landed 2026-10-15).

Expected: the CPU saves ~70 of 85 µs per read (replay and planning remain). The GPU then
needs about 2,490 GPU-s (today's seeding and rescue) plus the extension kernel. At kernel
v2 rates that is roughly 820 GPU-s for the 7.4G calls of a 30× run, an estimate. Alignment
becomes GPU-bound at roughly 3,300 s, down from 4,681 s: **about −1,400 s**.
Effort: 4–8 days with the agent workflow.

Exactness is not the open risk: the existing path already showed 0 mismatches over 403M
lookups (I-19, P-04). The risks are performance and resources:
- **Band retries.** Right h0 is left's *final* score after BWA's band-retry loop (`w <<= 1`
  when `max_off` is large). A fused kernel must run that loop itself, or send those seeds to
  the CPU; either way the lookup key changes from today's path.
- **Memory.** The opt-in path added +8.5 GB peak, and a third in-flight batch adds more
  (I-20). 30× peaks at 74 GB of 121 today; budget ≤ 90 GB.
- **Contention.** Unified memory is shared: while the seed kernel runs, the CPU's own
  alignment work slows by about a third (recorded 2026-10-10). Concurrent kernels may
  overlap usefully or fight for bandwidth.

### After option 1, the bottleneck flips

With extension on the GPU, CPU work falls to ~60 µs per read: ~42,000 core-s, about 2,100 s on
20 cores. The GPU then needs ~3,300 s, so it becomes the bottleneck and the CPU sits ~36% idle.
From there, only work that **frees GPU time** helps. GPU chaining (option 3) moves CPU work
onto the new bottleneck and would make alignment slower, not faster. So the order to parity is:

1. option 1 (extension offload);
2. option 2 (faster seeding and a faster extension kernel), which buys back GPU time;
3. only then option 3 or 4, if the measurements show GPU headroom again.

### 2. Faster seeding and extension kernels (frees GPU capacity after option 1)

Seeding is the largest GPU consumer today, ~1.5 s per 331k-pair batch: FM-index walks bound by
memory latency. A k-mer table was tried and made the GPU 4–15% slower, so it is ruled out.
Untried:
- order reads by leading k-mer so walks share cache lines;
- warp-cooperative SMEM extension;
- overlap the backward and forward passes.

Every 10% off seeding is ~160 GPU-s. Effort is uncertain, 3–10 days. Do it once option 1 has
made the GPU the bound. The extension kernel v3 (several calls per warp, fewer row syncs)
belongs here too: the original target was ~3× the v2 kernel.

### 3. GPU chaining (~28 µs per read of CPU)

SA lookups into the GPU-resident dense SA, plus BWA's chain building. That is ordered b-tree
insertion, so exactness means reproducing insertion order. It needs a chain oracle dump first
(board V-05).
- Saves ~25 µs per read of CPU, ~875 s of 20-thread time, but only if the GPU has headroom.
  After option 1 it does not (see above).
- Effort: 5–10 days.

### 4. All of BWA-MEM on the GPU (Parabricks' architecture)

Seeding, chaining, extension, rescue and pairing on the GPU; the CPU only formats output.
This is the end state if options 1–3 leave the GPU underused. Effort: 3–6 weeks. Each stage
gets its own oracle dump and gate.

### 5. Outside alignment: the ~320 s in duplicate marking and BQSR

Part of this is a sarek requirement, not inefficiency. nomics' 392 s of BQSR includes
ApplyBQSR writing recal.bam, which sarek publishes; Parabricks never writes one, because its
HaplotypeCaller applies the table on the fly. What is fixable: whole-genome mode streams
md.bam and recal.bam through disk several times, and building the BQSR tables during markdup
pass 2 (board I-06) removes one full read of md.bam. CPU-side, 3–6 days.

## Measure before building (1 day, after sarek's 30× bake-off run frees the Spark)

`scripts/wgs30_align_gate.sh` (one 30× chunk, ~7 min, byte-identical gate) with
`NOMICS_TIMING` and a GPU timeline (`nsys`, or `nvidia-smi` at 100 ms):
- GPU-seconds per batch for seeding, rescue and extension (with `NOMICS_GPU_EXTEND=1`);
- CPU split after S-08 (`perf`);
- how long each batch waits on the other's kernels.

These numbers turn the estimates above into a plan. Option 1's −1,400 s rests on the
kernel-v2 extension rate and today's seeding share. Both are measurable in that one run.

- **Measure concurrent throughput, not isolated kernels.** `nvidia-smi` utilization is
  kernel-active time, not SM saturation, and latency-bound FM-index walks read as busy while
  SMs idle. So measure extension-kernel throughput *while* seeding runs, and CPU cost per read
  with and without concurrent kernels. Summing isolated kernel times could be pessimistic
  (streams overlap) or optimistic (bandwidth contention).
- **Split today's 2,490 GPU-s into seeding vs rescue with `nsys`.** The estimates assume
  seeding dominates; that is unmeasured.
- **Hygiene:**
  - check that index huge-page coverage is ≥ 90% (the harness refuses below it);
  - pair the `NOMICS_TIMING` run with a timing-off control, since that flag cost ~25% in
    HaplotypeCaller;
  - run under the resource lock with nothing else on the box.

## What "closing the gap" would and would not mean

- With options 1 and 2, alignment ~2,700–3,300 s gives a 30× total of ~4,350–4,950 s, with
  every output identical to sarek's.
- nomics runs fastp semantics because sarek does, and prep's CPU competes with alignment.

## Profile on the fused default (2026-10-26, `scripts/profile_chunk.sh`, 30× chunk 1)

The align stage runs 199.9 s plain and 202.3 s under nsys. All three runs are byte-identical to
the reference chunk BAM. **The GPU is now the limit:** kernels are busy 181.3 s of the
202.5 s span, 90% as the union over streams.

| GPU work | Kernel-seconds | Launches |
|---|---:|---:|
| Mate rescue (`ksw_wave_kernel<5>`) | 85.5 | 89 |
| Seeding, first pass (`collect_slot_kernel<64>`) | 82.7 | 89 |
| Seeding, two-pass retry (`collect_slot_kernel<128>`) | 16.1 | 89 |
| Fused extension (`fused_bin_kernel`, all length bins) | ~60 | 970 per bin |

CPU self time, perf:

| Function | Share |
|---|---:|
| `mem_chain_from` | 22% |
| `unpack_2bit` (reference fetch) | 9% |
| `mem_chain_flt` | 9% |
| `global2`, `chain2aln`, `occ4`, `sort_dedup_patch_core`, `chain_window` | ~4% each |

What this changes:
- **Options that add GPU work now lose:** GPU chaining (option 3) and GPU region merge, until
  GPU time comes down.
- **The levers are GPU efficiency, and handing work back to the CPU:**
  - the seed kernel (99 s with the retry pass) and the rescue wavefront kernel (86 s) are the
    two largest kernels;
  - the cap-128 retry costs 16 GPU-s, about what the CPU spends re-seeding the same overflow
    reads, so it can move back to the CPU at no loss;
  - mate rescue can be split between the CPU and GPU with work stealing, as the extension
    split did (`NOMICS_GPU_EXTEND_CPU_LANES`).
- **A stale binary first produced a misleading profile:** the CPU-extension path, align 266 s
  and 43% GPU busy. The repo sync to the Spark kept the Mac's mtimes and cargo skipped the rebuild; it now
  syncs by checksum.

## Mate rescue: the kernel and a CPU/GPU split (2026-10-26, branch `perf/rescue`)

Run dirs: `~/omics/runs/rescuek-*`. The real-shape replay is the first rescue batch of
30× chunk 1, written with `NOMICS_DUMP_RESCUE` (`rescuek-dump/batch0.bin`, 744,199 calls, 1.2 GB,
`align2.bin` format): `ksw_bench batch0.bin 20 1000000`. Every call there has a 151-bp query and
99.3% have a 1,334-base window, so length binning is not a lever.

### The limiter

ncu cannot read counters here (`RmProfilingAdminOnly=1`, a system setting for the owner), so
the kernel was taken apart with ptxas, SASS, occupancy padding and clock64 phase timers:
- **Instruction issue in the step loop.** `ksw_wave_kernel<5>` spent ~172 instructions per
  wavefront step for 5 cells (~34 per cell): a 64-bit shift per score lookup, the row maximum
  and its column as two compares and selects, 5 shuffles plus a stop broadcast every step, the
  stripe-reset test on every column, and per-step rebuilds of lane constants and the warp's
  shared-memory pointer. Registers 76 and 16 KB of shared memory per block cap it at 24 warps
  per SM.
- **Occupancy matters but is not the limiter.** Padding the pair kernel's shared memory: 12
  warps/SM 0.246 s, 8 warps 0.274 s (+11%), 6 warps 0.339 s, 4 warps 0.451 s.
- **The sequential second-best replay** (lane 0 walks every row of the window after the
  forward pass) was **31% of the kernel** once the DP was cheaper (phases: load 2%, forward
  56%, second-best scan 31%, reverse 11%).

### Changes (all bit-exact with `ksw_align2`: score, te, qe, score2, te2, tb, qb)

1. One sign-extending PRMT per score; the row maximum and its column as one key
   (`vp * 1024 + 1023 - k`, max = largest vp at the smallest k); H and F in one shuffle; the stop
   flag checked every 8 steps with the last lane frozen at its stopping row; the stripe reset
   only on a lane's first column when P = 32·C.
2. **Two calls per warp in 16-bit halves** (`VIADD.16x2`, `VIMNMX.S16x2[.RELU]`,
   `VIMNMX.U16x2` are single instructions on sm_121): byte-sized calls with P = 32·C and
   windows of bases 0-3 pair up (the batch: all 372,099 pairs). One PRMT gives both halves'
   scores. The row key is `min(vp, 255) * 256 + 255 - k` per half; saturating at 255 is exact
   because a byte-sized call stops at the first row reaching `min(255 - shift, endsc) <= 255`.
   H <= qlen · max score (padding columns score 0, gap penalties are non-negative), and
   `pair_ok` admits only calls with qlen · max score < 30,000 and gap terms below 8,192, so no
   half wraps. Other
   calls (N in a window, 16-bit calls, other widths, odd counts) run one per warp.
3. **Second-best by the whole warp.** Row i sets BWA's (S, E) exactly when
   `q(i) && (!active(i-1) || s(i) > s(i-1))`, a two-state map per row; lanes take contiguous
   chunks, compose the maps with a shuffle scan, find their run ends, and two reductions pick
   the largest eligible S at its first end.
4. F from vp: `max(f - e_ins, h - oe_ins) = max(f - e_ins, vp - oe_ins)` for o_ins >= 0
   (negative gaps are now refused at launch and fall back to the CPU).

| Replay (30× batch, uncontended) | Kernel | |
|---|---:|---:|
| `ksw_wave_kernel<5>` (base) | 0.527 s | 1.0× |
| + PRMT, keys, 3 shuffles | 0.461 s | 1.14× |
| + two calls per warp | 0.313 s | 1.68× |
| + warp second-best | 0.239 s | 2.2× |
| + F from vp (shipped) | 0.232-0.240 s | 2.2× |
| shipped, one call per warp (`NOMICS_KSW_PAIR=0`) | 0.383 s | 1.38× |

Per step the pair loop runs ~117 instructions for 10 cells (~12 per cell). Phases after:
load 3%, forward 77%, scan 5%, reverse 16% (a third of the pairs run their reverse passes one by
one, when a qe puts a reverse query outside P = 160).

### On the chunk (`scripts/profile_chunk.sh`, chunk 1, both byte-identical)

| | Base (1ac7dd4) | `perf/rescue` (4f02c9f) |
|---|---:|---:|
| Align, plain | 202.0 s | 191.6 s |
| Align, under nsys | 205.6 s | 195.1 s |
| Rescue kernel GPU-s (89 launches) | 62.9 | 32.0 |
| All kernels GPU-s | 210.2 | 180.3 |
| GPU busy (union) | 167.7 of 205.7 s (82%) | 139.2 of 197.6 s (70%) |
| Rescue on each batch's tail, mean | 1.50 s | 0.96 s |

(This morning's profile billed rescue 85.5 GPU-s; the base rerun bills 62.9. nsys charges each
kernel its wall span, which grows with whatever shares the SMs.) Under perf both run ~200 s:
with the GPU at 70% the CPU is now the limit, and perf's sampling costs CPU.

### The CPU/GPU split

`NOMICS_RESCUE_CPU_THREADS=N` lets N CPU workers (from the batch's own) take calls from the
long end of the length-sorted batch while GPU lanes launch from the short end
(`NOMICS_RESCUE_GPU_GUIDED`, default 3/4 of what is left per launch, at least
`NOMICS_RESCUE_GPU_TILE`, default 16,384; `NOMICS_RESCUE_GPU_LANES`, default 2).
`NOMICS_RESCUE_CPU_SHARE=f` fixes the longest f of the calls to the CPU instead (1 = all on the
CPU) and runs the GPU part as one launch. Results are stored by call index, so who computes what
cannot change them (`rescue_split_matches_cpu_whatever_the_split`; on the chunk, the
share-0.25 batch dump recomputed on the CPU).

Sweep: `scripts/rescue_split_sweep.sh` on chunk 1 with the shipped kernel (fb1dce2), timing on,
each config under the lock for its own run. Run dir `rescuek-sweep2`. Every config gave a
byte-identical BAM, and `ksw_bench` on the share-0.25 config's batch dump (its longest quarter
computed by CPU workers) found 0 mismatches in 744,199 calls.

| Config | Align | Mean GPU (smi) | Calls on the GPU | Rescue on the batch tail |
|---|---:|---:|---:|---:|
| **all on the GPU, one launch (default)** | **193.4 s, 192.0 s** | 53%, 52% | 100% | 1.12 s, 1.24 s |
| 8 CPU workers stealing | 193.9 s | 53% | 95.2% | 1.20 s |
| 20 CPU workers stealing | 198.1 s | 49% | 93.8% | 1.25 s |
| fixed CPU share 0.25 (+ batch dump) | 224.2 s | 46% | 75.0% | 1.62 s |
| fixed CPU share 1 (all planned calls on CPU workers) | 351.1 s | 26% | 0% | 5.92 s |
| `NOMICS_GPU_RESCUE=0` (rescue inside pairing, CPU) | 337.9 s | 25% | — | — |

A first sweep with the first kernel change (sweep1, partly overlapped by builds and tests, so
not used above) agreed: 4 and 8 CPU workers with fixed 64k tiles lost 11-22 s.

The split does not pay: one batch's calls take 1.95 s on 20 threads (39 core-s) against 0.24 s
of the whole GPU, so every call moved to the CPU costs far more than it frees, and with the GPU
at 70% busy the CPU is the tighter side. The default stays all on the GPU in one launch.

### Gates (shipped tree, fb1dce2 kernel)

- Replays: the 30× batch dump and `bench/ksw/align2.bin`, 0 mismatches, paired and unpaired;
  `ksw_align2_gpu` (new: paired calls, a match score that reaches 255, unequal gaps) passes on
  CUDA with and without `NOMICS_KSW_PAIR=0`, and on Metal (`--features metal`, Mac).
- 30× chunk 1 byte gate: every sweep config above, identical.
- `NOMICS_DIAG=1` on the chunk: identical, 197.2 s, 1,148 canary sweeps with no damage.
- `scripts/full_exome.sh tier1`: 18/18 PASS (total 300.8 s, align 209.1 s).
- CPU-extension path (`NOMICS_GPU_EXTEND=0`), chunk 1, three runs each: 275.2, 262.0, 272.4 s
  against base 270.7, 281.4, 271.1 s (means 269.9 vs 274.4 s), all identical. No regression;
  its rescue tail drops from 1.30 to 0.97 s per batch.

### Tried and dropped

- Fixed 64k-call GPU tiles with two lanes per batch: +11 s on the chunk (sweep1, 214 vs 203 s
  with the first kernel change); each launch's pack and wait sits on the batch's tail.
- Length binning: the batch's windows are 99.3% one length.
- Templating the pair loop on equal gaps and on "no cell can reach 255": 3% on the replay,
  not worth four more instantiations.

## Seed kernel limiter (2026-10-26, branch `perf/seed-kernel`)

Run dirs: `~/omics/runs/seedk-*` (`seedk-bench` holds the variant logs). ncu and nsys GPU
metrics are unavailable (`RmProfilingAdminOnly=1`, `ERR_NVGPUCTRPERM`), so the limiter was
pinned down with in-kernel counters, microbenchmarks and controlled variants of
`collect_slot_kernel<64, PackedList>` on one 662,000-read 30× batch (`seed_bench`, GPU only,
digest-checked: every variant gave the same 26,329 overflows and 7,617,869 intervals).

**What the kernel does per read:** 856 `extend` calls (1,712 `occ4`), split p1 forward 148,
p1 backward 143, reseeding forward 33, reseeding backward 393 (46%), LAST seeds 138. Only 38%
of extends read one BWT block for both `occ4`s, so a batch touches ~916M distinct 64-byte
blocks: ~1.5G per second at 595 ms. SIMT efficiency at the extend call sites is 8–11% in
`smem1` and 48% in the LAST-seed loop. There are 1,800 interval-list operations per read.

**Ceilings** (`seedk-tmp/rnd*.cu`, random dependent reads over a 3.2 GB huge-page host buffer):
unloaded latency 0.73 µs; the fabric saturates at ~2.6G random requests per second whether a
request is 16 or 64 bytes, from about 3,000 requests in flight.

**Levers tried, isolated batch (baseline 585–600 ms between rounds):**

| Variant | ms | Verdict |
|---|---:|---|
| Occupancy: 4 / 6 / 8 (default) / 10 / 12 blocks per SM | 588 / 605 / 597 / 633 / 758 | flat, then spills |
| 2 blocks per SM (256 threads) | 610 | still flat |
| Both `occ4` blocks loaded as `uint4`s before counting | 613–1,390 | registers 128–167, spills |
| `prefetch.global.L1` of the second block, and of every backward entry's blocks | 600 | none |
| Same, `L2` | 640 | worse |
| +16 / +64 dummy instructions per `occ4` | 607 / 609 | not issue-bound |
| Reads ordered by leading 12-mer | 597 | none |
| BWT and/or list scratch in device memory (probe only) | 584–620 | none |
| Shared-memory carveout 25% instead of 0 | 687 | L1 size matters |
| Per-thread work pulling (persistent, 5 ms … ∞ budget) | 586–602 | none |
| k-mer table k=12 (exact, existing `seed_gate` path) | 563–576 | −3 to −5% |
| Fewer list reads in `smem1` (kept) | 582–595 | −0 to −2% |

**Reading:** the kernel is insensitive to threads in flight (2× fewer or more), memory-level
parallelism, instruction count, placement and block residency. It responds only to L1
capacity and to the number of BWT lookups, so it sits on the memory system's random-miss
throughput. Its requests are ~60% of the synthetic ceiling, but the synthetic kernel has no
divergence and no list traffic. The 2.6× of 2026-10-19 removed a pathology (word-interleaved
local memory). No second one turned up.

**A curiosity that is not a lever:** in steady state only ~73% of the resident-block slots are
busy (uniform across SMs; replacement blocks arrive 2–9 ms after a block ends). A synthetic
kernel with the same block shape, slot pool, carveout, local memory, stream priority and host
traffic stays at 100%. 32-thread blocks also stay at 100%. Work pulling, which keeps every
lane busy, did not make the batch faster.

**Kept** (exact; `seed_gate` tier0 and a 30× batch, all modes plus the table path, identical to
the CPU): `smem1` serves the last `curr` size and last `mem` start from registers, skips the
final `mem` reversal (unobservable: the caller sorts by `info`) and reads the forward list from
its end instead of reversing it.

**30× chunk 1, fused default** (`profile_chunk.sh` and interleaved runs; every BAM
byte-identical to the base build's):

| | Base (`1ac7dd4`) | This branch |
|---|---:|---:|
| Align, plain runs | 199.0, 201.9, 202.4 s | 193.5, 195.8, 190.2 s |
| Seed kernel, nsys | 85.2 GPU-s | 82.6 GPU-s |
| GPU busy | 81% | 85% |

`NOMICS_SEED_TABLE=<cache>` (opt-in) adds the k-mer table: 81.9 GPU-s, align 192.9 s, the
same BAM. That is within noise of the branch, so it stays off.

**What would move seeding further:** fewer BWT requests per read, not more parallelism. That
means bwa's algorithm with less work. One example is reusing pass 1's forward walk for the
first LAST seed. GPU counters would settle the miss breakdown; they need
`NVreg_RestrictProfilingToAdminUsers=0`, which is the owner's call.

## Combined branch: rescue kernel, seed trims, seeding default (2026-10-26, `perf/gpu-kernels`)

Every row below is gated: the chunk BAM is byte-identical to the base build's, tier1 is 18/18,
and both 30× runs are 41/41 (42 PASS lines).

**30× chunk 1, fused default** (plain runs of `wgs30_align_gate.sh`):

| Build | Seeding | Align |
|---|---|---:|
| Base `1ac7dd4` | cap64 | 199.0, 201.9, 202.4 s |
| Rescue kernel only | cap64 | 191.6 s |
| `f244f4e`, rescue + seed trims | cap64 | 184.3 s (GPU busy 71%) |
| Same build, interleaved A/B | cap64 | 180.4, 181.2 s |
| Same build, interleaved A/B | two-pass | **168.4, 168.6 s** |

The earlier cap64 default (`3a26b41`) was chosen when the GPU was the bottleneck: cap64 trades
CPU seeding work for GPU work. Once the rescue kernel had halved its GPU time, the GPU had room
again. Two-pass then won by 12 s per chunk (−6.8%), so `7e8cda7` makes two-pass the default
again. Metal keeps cap64.

**Full 30×** (both runs on a Spark with 29–34 h of uptime; the GPU buffers fell back to small
pages, and the second run's index did too):

| Build | Seeding | Prep | Align | Markdup | BQSR | HC | Total |
|---|---|---:|---:|---:|---:|---:|---:|
| `5ddb62f` main, clean after reboot | two-pass | — | 2,383 | — | — | — | 3,818 s |
| `f244f4e` | cap64 | 342 | 2,406 | 304 | 204 | 938 | 3,871 s |
| `7e8cda7` | two-pass | 361 | 2,451 | 334 | 203 | 948 | 3,978 s |

Fragmented memory hides the chunk-level gain. Main's own fused runs (two-pass) on fragmented
memory ranged 3,974–4,159 s; this branch's two runs fall at or below the bottom of that range. The second run was slower in every stage that the seeding change does not touch,
too (prep +19 s, markdup +30 s). These two runs therefore do not separate cap64 from two-pass at
30×.

**Interleaved full-30× A/B after a reboot** (2026-10-27, the merged branch build, runs back to
back starting 81 s after boot; `~/omics/runs/seedab30-20261027`, archived under
`runs-archive/session-20261027`). All four are 41/41 (42 PASS):

| Run | Seeding | Prep | Align | Markdup | BQSR | HC | Total | Memory |
|---|---|---:|---:|---:|---:|---:|---:|---|
| 1 | two-pass | 323 | **2,172** | 286 | 201 | 938 | **3,619 s** | clean |
| 2 | cap64 | 333 | 2,471 | 290 | 203 | 939 | 3,930 s | index on small pages (78%) |
| 3 | cap64 | 351 | 2,474 | 294 | 205 | 939 | 3,934 s | index on small pages (93–95%) |
| 4 | two-pass | 443 | 2,570 | 294 | 205 | 934 | 4,024 s | index on small pages (91–94%) |

Only run 1 had the index on huge pages; from run 2 on, the 2 MiB collapse failed with ENOMEM,
and the GPU buffers were on small pages in all four. Memory state therefore dominates: the
fragmented runs cost 300–400 s of alignment, while the seeding mode moves it by less than the
run-to-run drift (run 4 also lost 120 s in prep, which seeding does not touch). The A/B does not
overturn the chunk result, so two-pass stays the default.

Run 1 is the best full 30× so far: **3,619 s (1 h 00 min)**, 199 s under main's clean 3,818 s,
almost all of it in alignment (2,172 s vs 2,383 s, −8.9%). That is 11.3× sarek's 41,044 s.
A cleaner seeding comparison would need a reboot (or a successful compaction) before every run.
