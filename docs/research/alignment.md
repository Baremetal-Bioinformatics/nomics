# Alignment: where the time goes and what could cut it (2026-10-09)

Alignment is 575 s of the 704 s tier1 `nomics run` (82 %). These notes are research only; nothing here has been tried yet.

## Measurements (tier1 chunk 0001, 3.37M pairs, Spark)

Per 662k-read `-K` batch, 20 threads, GPU seeds overlapped:

| phase | wall | CPU busy | per read |
|---|---:|---:|---:|
| align (chain + extend; CPU seeding for ~3 % GPU overflow) | ~2.5 s | ~50 s | ~75 µs |
| pestat | 0.03 s | | |
| pair (mate rescue + pairing + SAM text) | ~1.2 s | ~24 s | ~36 µs |
| batch wall incl. gaps | ~4.5 s | | |

- **Whole chunk:** 52.6 s wall, 768 user-s. Align 46.3 s, sort 2.2 s, write 1.5 s, load 2.1 s.
- **Gaps:** about 0.8 s per batch is outside the two parallel phases: reading, waiting on the GPU, phase tails.
- **Single thread, uncontended, 20k pairs (`align_profile`, `pair_profile`), per read:**
  - seed 57.5 µs (now on the GPU, except for overflow)
  - chain (SA + tree) 15.8 µs
  - chain filter 1.6 µs
  - **extend 22.9 µs**
  - **mate rescue 26.4 µs**
  - rest of pairing / SAM 2.7 µs
- **Contention.** The CPU parts add up to about 71 µs per read single-threaded, against about 111 µs per read in the 20-thread pipeline: a 1.55× inflation. About 1.35× of that is the platform (10 of the 20 cores are A725s, plus all-core clocks; HaplotypeCaller shows the same). The rest is consistent with GPU seeding competing for memory bandwidth (bench/align.md).

## Constraints

- **Batch structure is fixed.** Output depends on the fastp chunks and on bwa's `-K` batches, because pestat (insert-size statistics) is estimated per batch. The ~120 batches (12 chunks × ~10) are therefore independent units. They can be reordered or overlapped freely, as long as each one is processed as a whole.
- **Integer DP is GPU-safe.** ksw extend, global and align2 are integer dynamic programming, so a GPU port has no floating-point exactness risk. Only the tie-breaking semantics, which follow evaluation order, have to be matched; that is the same discipline as PairHMM.

## Ideas, ranked by expected payoff / risk

1. **Batch-level pipelining across all 120 batches.**
   - Keep two batches in flight, so phase barriers, the pestat barrier and batch-boundary gaps fill with the other batch's work.
   - Start aligning while prep is still running. A chunk's first batch exists as soon as its first fastp packs arrive (pack `p` goes to chunk `p % 12`).
   - Exact by construction, because batch content is unchanged.
   - Expected: most of the 0.8 s per batch (~90 s total), plus prep's 30 s hidden.
2. **Mate rescue on the GPU, batched like PairHMM.** Rescue is 26 µs per read, the largest single CPU item.
   - The jobs are small, independent integer SW problems (read vs. a reference window). They are compute-bound with little memory traffic, so they contend less with the CPU than GPU seeding does.
   - The GPU has to reproduce the striped algorithm's results: the score, the end positions, and the second-best rows (`second_best`).
   - Gate: dump rescue calls and compare, as for PairHMM.
   - Could take about a third of the CPU work off. A hybrid split, as in HaplotypeCaller, keeps the CPU busy too.
3. **Extension (`ksw_extend2`, 23 µs per read).**
   - The per-seed control flow (seed coverage checks between extensions) is sequential within a chain, but extensions across reads are independent.
   - Two options:
     - Batch all the left extensions of a pass on the GPU. This changes control flow and needs care.
     - Vectorize bwa's scalar banded `ksw_extend2` across *reads* (several reads per SIMD lane, like the PairHMM lanes). This is exact per lane.
4. **SAM text → BAM records directly.** Today the pair phase formats SAM text, and `sort_sam` parses it back.
   - Emitting BAM records directly saves the formatting plus 2.2 s of sort/parse per chunk, about 30–40 s in total.
   - Gate: chunk BAM records, unchanged.
5. **Overlap MarkDuplicates with alignment.** Pass 1 (read ends) can be built per chunk as each one finishes; only the final selection needs all chunks. That hides most of markdup's 42 s.
6. **Seeding traffic.**
   - A k-mer → SA-interval lookup table (e.g. 12-mers: 16.8M entries, about 0.4 GB) would skip the first backward-extension steps of every SMEM. The intervals are identical, so this is exact.
   - It cuts GPU seed-kernel time and, more importantly, the memory traffic that slows the CPU.
   - Separately, raising the 64-interval GPU cap would bring the 3 % CPU-seeding fallback (about 1.7 µs per read) towards 0.
7. **Chain SA lookups (15.8 µs per read).** A full SA (sampling 1 instead of every 4th, about 16 GB) would cut the SA walk. It fits in memory; it is a memory-for-time trade.

## Suggested order

1. Batch pipelining (1) and BAM-direct output (4): cheap and exact by construction.
2. Measure again.
3. GPU mate rescue (2), with a dump-based oracle.
4. Extension vectorized across reads (3), then the seeding table (6).

Rough ceiling if 1–5 land: alignment from ~575 s to ~300 s, and FASTQ → VCF from ~11.7 min to ~7 min.

## Outcomes (2026-10-09, all on tier1 unless noted, every pipeline run gated 18/18)

| # | Idea | Result | Verdict |
|---|---|---|---|
| 1 | Two `-K` batches in flight | align 610.6 → 594.7 s on its own; 556.5 s with idea 2 | kept |
| 2 | BAM-direct output | align 610.6 → 573.8 s on its own. Its first version wrote bwa's `MIDSH` op codes as BAM codes, turning soft clips into `N`; the equivalence test and the gates caught it, and it is fixed. | kept |
| 3 | GPU mate rescue (`ksw_align2`) | first try (one thread per alignment emulating the striped vectors in global memory): exact, but 2.99e9 cells/s against 2.53e10 for the 20-thread CPU. **Superseded 2026-10-13, see below: a warp-wavefront kernel is 4.4× the CPU and is kept.** | ruled out, then kept in a new form |
| 4 | Extension, reads in SIMD lanes | exact (1.2M calls), 0.26–0.56× of scalar; lane utilization only 31–57 % | ruled out |
| 5 | Markdup overlapped with alignment | chunk-local preparation (6–12 s per chunk) now hides behind alignment. The global merge (13.0 s), selection (7.4 s) and flags/metrics (8.3 s) remain, and the stage is still ~42 s. | kept (exact); next: parallelize the global steps |
| 6 | k-mer table + GPU cap | the table is exact, but CPU seeding only gains 6–9 % on the ~3 % fallback, and the GPU gets 4–15 % slower. Caps 96/128 make the kernel 7.6–9.0 s per batch instead of 2.0 s. CPU contention from seeding is only 1.07–1.10×. | ruled out |
| 7 | Dense SA interval 1 | chaining 15.5 → 6.9 µs per read, exact. Tier0 passes, but on tier1 GPU seeding fails with CUDA error 400 when the 30 GB SA is resident. **2026-10-13:** interval 2 kept (30× chunk 430 → 414.5 s); interval 1 ran without the CUDA error but was slower overall (475.6 s, 79% huge pages). | interval 2 kept |
| 8 | NEON rows for the striped `ksw_u8`/`ksw_i16` (`nomics-simd`), replacing array lanes | exact on the 62k-call dump and on a 70k-call WGS dump; 1.5× single-threaded (Spark 4.60e9 cells/s); a shorter F dependency chain (exact rewrite) gains nothing further in the full kernel | kept: tier1 align 521 → 468 s, wgs20 398 → 254 s |
| 9 | Extension rows in NEON (F as a prefix max of `t[m] + (m+1)·e_ins`, H and E per column; argmax per lane) | exact on 1.96M dumped calls (exome + WGS), but no faster: 623 vs 653 Mcell/s (exome), 661 vs 675 (WGS) on Apple M; rows average ~40 cells and the wide scalar cores already run bwa's loop near the vector rate | ruled out |
| 10 | Two rescue alignments interleaved per thread (resumable u8 pass; across a pair's two independent rescue streams, and speculatively across consecutive calls of one stream) | exact (dumps and 18/18 gates), 1.45× per alignment on Apple M, but only 1.08× on the Spark's cores, and rescue per WGS read stayed at ~90 µs | ruled out; reverted |

Best fully gated tier1 end to end: 679.3 s (ideas 1, 2 and 5, SA interval 4), down from 704.4 s.

Best fully gated tier1 end to end is now 660.6 s (2026-10-10, below).

## Regressions (resolved 2026-10-10)

- **Memory:** peak RSS rose to ~73–79 GB because each chunk's FASTQ text was held for the whole run. Fixed in `d0b3ece`; it is now 57–61 GB.
- **Align 575 → 610 s: noise, not code.** In an interleaved tier1 A/B, each build run twice, the same binary varied by up to 39 s. The commit that produced 575 s ran 587 and 626 s; later builds ran 556–572 s. Align claims under ~30 s need three or more runs, or a direct timer.
- **A real tail regression, found on the way.** Since markdup overlap, the mapped-BAM writer and markdup preparation ran on one thread for the whole run. The last chunk's write and prep (~21 s) therefore ran after alignment, on one core. The writer now takes the full thread budget once the last chunk completes, and the tail is 3.6 s.
- **Index load during prep.** The index and dense SA (11–14 s) now load in a thread while prep runs, and alignment never waits for them. Prep's peak RSS rises from 30.7 to 43.4 GB, under the run's 60.6 GB peak.

| tier1, 18/18 gates | prep | align | total |
|---|---:|---:|---:|
| before (2 runs) | 30.8 / 30.5 | 570.1 / 556.2 | 700.4 / 685.7 |
| writer widened + index during prep (2 runs) | 32.4 / 35.2 | 528.5 / 524.1 | 660.9 / 660.6 |

## Whole-genome reads (2026-10-11, wgs20 chunk 0001, 20k pairs, single thread, NEON rescue)

| per read | tier1 exome | wgs20 |
|---|---:|---:|
| seed (GPU in the pipeline) | 57.3 µs | 60.5 µs |
| chain (SA + tree) | 15.3 µs | 21.6 µs |
| chain filter | 1.6 µs | 6.9 µs |
| extend | 23.7 µs | 80.9 µs |
| mate rescue | 9.8 µs | 91.2 µs |
| rest of pairing / SAM | 2.7 µs | 4.4 µs |

- WGS rescue: 1.76 `ksw_align2` calls per read over 1,339 bp windows (wider insert spread than the exome's 785 bp); 91 % of calls find the mate, so none can be skipped. ~355k cells per read; the forward pass is ~80 % of rescue time and runs near the single-problem limit of the kernel.
- Next levers: interleave two independent rescue problems per thread (the kernel is latency-bound per problem), and extension (3.4× the exome's cost per read on WGS).

## GPU mate rescue, second design (2026-10-13; kept, `7eed225`, `db1b166`)

- **Kernel.** One warp per alignment over BWA's padded query; lane `l` owns C columns and handles target row `s - l` at step `s`, taking its left boundary from lane `l-1` by shuffle. Only the striped kernel's observable recurrences are reproduced, not its lane order: E from the pre-lazy-F value (insertions within a stripe only, reset at stripe starts), H including all insertions, row maxima and the first column of the maximum (BWA's qe tie-break picks the smallest index, so order does not matter), score-0 padding columns, the u8 255 and XSTOP breaks, `record_row`/`second_best`, and the XSTART reverse pass.
- **Exactness.** 0 field-complete mismatches on 132,248 dumped calls (tier1 and wgs20), first try.
- **Speed.** v1 ran at 2.2e10 cells/s. Register score tables and exact-fit columns per lane gave 2.4–3.1e10. Staging each target and its row maxima in shared memory gave **2.4–2.7e11** (4.4× the 20-thread NEON CPU). Per-step reads of host-mapped memory were the bottleneck, and probably the old kernel's too.
- **Integration.** After pestat, a parallel planner lists every pair's rescue windows from the pre-rescue regions, using `mem_matesw`'s own window code. One launch per batch runs them, and pairing looks each call up by (mate, strand, window). A miss is computed on the CPU, so speculation costs time, never output. Per 30× batch: plan 0.11 s, pack 0.13 s, kernel 0.57 s.
- **Result.** wgs20 align 252 → 207 s; 30× chunk 1 594.5 → 491.8 s. With 20 workers per in-flight batch (oversubscribed, so the batch waiting on the GPU yields its cores) it is 439.4 s. All byte-identical.
- **Caveat.** GPU seeding is 7× slower when the index lands on 4 KiB pages; rescue launches then queue behind it (one 30× chunk took 1,688 s). Timing runs now refuse to start below 90% huge-page coverage.
- **Next lever.** Seed extension (~6 s per 30× batch, CPU). A row-at-a-time warp kernel (prefix-max F, ballot band trim, per-row z-drop) could be exact. First measure how much speculative work `chain2aln`'s containment skips and band retries would waste.

## GPU seed extension (2026-10-13; exact, neutral, opt-in `NOMICS_GPU_EXTEND=1`)

- **Workload (30× chunk 1, counters under `NOMICS_TIMING`).** Per record: 10.0 chains, 13.1 seeds, 10.1 extended (77% of seeds), 4.9 left + 5.7 right `ksw_extend2` calls, ~0 band retries. Speculating every seed's first-band extension therefore costs ~1.3× the executed work.
- **CPU profile per record** (CPU path, 40 oversubscribed threads, so wall per thread): chaining 42 µs; `chain2aln` 100 µs, of which `ksw_extend2` is 85 µs. Extension is ~60% of align1. The older single-thread profile above understated chaining.
- **Kernel.** One warp per call, one target row at a time. F is an exclusive warp prefix-max of `max(M − oe_ins, 0) + (k+1)·e_ins`. E and M are per lane. H/E persist in shared memory exactly like BWA's `eh[]`, including stale cells beyond the band. The band trim comes from registers; the row argmax (ties to the last column) is a packed `redux.sync`; z-drop is a per-row break. Exact on its first run: 0 mismatches in all six fields over 1,956,277 dumped calls.
- **Integration.** Per batch: chains for every read, then all left extensions in one launch, then all right extensions (h0 = left score) in a second, then `chain2aln` with lookups by (chain, seed, side, band, h0) and CPU fallback. tier1 and wgs20 18/18; 403M lookups, 0 misses.
- **Why neutral.** The batch's CPU phases halve, but each batch now waits on two launches. Per side on wgs20 (7–9M calls): plan ~0.3 s, pack 0.2–0.3 s, kernel 0.7–1.1 s, 1.5–2.9 s wall including queueing behind the other batch's launch. Align: tier1 499 vs 432 s, wgs20 203 vs 200 s.
- **To pay off.** A kernel ~3× faster (several calls per warp, fewer per-row sync points), planning straight into the GPU's host buffers, and launches overlapped with the other phases.

## End of the 2026-10-12/14 iteration window

The 30× genome went from 9,994 s to 7,135 s end to end (final run; two earlier near-final builds 7,186 and 7,281 s), 41/41 identical. Alignment went from 6,926 s to ~4,970 s, with prep (~550 s) now overlapped inside it. What remains of alignment per record, on the CPU: ksw_extend2 ~85 µs, chaining ~28 µs (SA interval 2), the rest of chain2aln ~15 µs, plus pairing.

Open items:
1. **GPU extension needs a ~3× faster kernel** before it pays: several calls per warp, fewer per-row sync points, planning straight into the GPU's buffers, launches overlapped with CPU phases.
2. **The root cause of index paging is unfixed.** Coverage ranged 19–100%; the GPU also failed to allocate once when the page cache was full (board V-09).
3. **Profile with perf first.** It now works on the Spark; this window's hand timers misled three times.

## GPU extension host scheduling (2026-10-15; Spark measurement pending)

The reported 30× chunk baseline is CPU align **415 s**, GPU extension **458 s**
with **+8.5 GB** peak RSS (20 CPU threads). These are pre-change measurements;
there is no measured speedup for the changes below.

Source inspection confirms the global extension mutex covered input packing,
CUDA execution and the result copy across all in-flight batches. Unreported work
also included serial input/key/read-range assembly; the old `total` stopped
*before* result conversion and planner buffer destruction. The planner used
repeated `Vec::remove` for oversized queries. Within a batch all chains, both
extension sides, and final extension replay ran as whole-batch barriers. Right
h0 really depends on the left score; those calls cannot simply run together.

The opt-in path now schedules independent tiles of **65,536 reads**, with **two
lanes per in-flight batch** by default. Lanes divide `NOMICS_ALIGN_WORKERS` (or
its existing default) rather than each taking another full CPU thread budget.
Each lane owns one reusable engine, with no global extension mutex. The CUDA
shim uses a private nonblocking stream and synchronizes it even after errors
before the Rust buffers can be reused. The kernel arithmetic and C ABI layout
are unchanged. Original `-K` batches, per-read output slots, pestat and pairing
are unchanged.

Plans remain in 512-read blocks. Iterator groups feed parallel validation and
packing directly into the engine's `HostBuf`s, without a batch-sized
`GpuExtendInput` vector. Keys/read ranges stay in blocks, and results convert in
parallel directly from the engine's output buffer (no intermediate result
copy). Oversized queries are rejected before copying their bases. Planner base
arenas and packed bases **still coexist during a tile**; tiling bounds their
lifetime and the number of retained chains/results, but does not remove that
remaining copy. Engine capacities grow to a lane's largest tile and are freed
at batch completion. A read-count tile limit is not a hard byte limit on
repeat-heavy input.

`NOMICS_TIMING` retains the `gpu extend` lines and adds execution wall time and
post time (conversion plus planner cleanup). `lock-wait none (private engine)`
means the old GPU mutex no longer exists, not a sampled zero. `execute` includes
GPU scheduling/event/stream overhead; `kernel` is the CUDA-event interval and
is a subset of `execute`. `total` now includes post work. Concurrent tile times
must not be added and called align wall time.

Spark follow-up:

- Build and test with `--features cuda`. The new packing tests also run on the
  Mac; `tiled_extension_matches_cpu_in_input_order` compares full regions to the
  CPU on a generated tiny index, including long-query fallbacks, partial tiles,
  and several worker/lane counts. On the Mac it exercises unavailable-GPU
  fallback; on Spark it also exercises GPU caching. The GPU-only
  `concurrent_engines_match_serial_results` test checks four simultaneous
  engines and repeated buffer reuse.
- Replay the dumped extension calls (`extend_bench`), then run tier1/wgs20 and
  the same 30× chunk through the existing exact-output gates. Check all six
  extension fields, final output bytes/digests, and fallback warnings/misses.
- Under the resource lock, compare repeated saturated CPU control runs against
  the new opt-in path with identical input, seed/rescue settings, in-flight
  batch count, thread count, and huge-page coverage. Record align wall time,
  CPU-seconds, peak RSS, summed calls, and distributions of plan/pack/kernel/
  execute/post/total. Concurrent streams can contend with seeding and rescue;
  inspect those timings too.
- Sweep `NOMICS_GPU_EXTEND_READS=32768,65536,131072` and
  `NOMICS_GPU_EXTEND_LANES=1,2`. Read limits round up to 512-read blocks. Lanes
  cap at the worker count and available tiles. A large read limit with one
  lane gives whole-batch barriers for isolating the tiling effect. More tiles
  trade retained memory and CPU/GPU overlap for more launches/thread creation;
  the defaults are candidates to measure, not demonstrated winners.

Do not enable GPU extension by default unless it beats the saturated CPU path
end to end while passing the exact-output gates. No stage or gate status has
been promoted by this host-scheduling change.
