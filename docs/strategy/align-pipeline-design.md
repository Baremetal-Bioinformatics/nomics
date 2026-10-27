# Design: an alignment pipeline without barriers (2026-10-20)

Goal: shorten 30× alignment by keeping the GPU and the
CPU busy together. Every output stays byte-identical. Evidence and the gap breakdown are in
`docs/strategy/structural-gap.md` and `docs/strategy/gpu-extend-fused-design.md`.

## Where the time goes today

On 30× chunk 1, fused path, both sides are about three-quarters busy:

| Measure | Value | Source |
|---|---|---|
| GPU busy | 206 s of a 271 s span (76%): extension 101 s, rescue 84 s, seeding 80 s | clean nsys trace |
| CPU busy, align phase | 14.4 of 20 cores on average (72%); per-second median 14.7, p90 18.5; ≥ 19 cores for only 19 of 276 s | system sampler at 1 s, 2026-10-20 |
| CPU memory traffic | ~16 GB/s | calibrated counters |

That puts the work per chunk at ~4,100 core-seconds of CPU (~205 s on 20 cores) and ~206 s
of GPU. **The two sides are balanced, and neither is saturated.** Alignment waits on
synchronization. Removing the barriers can at best bring a chunk from ~292 s to the ~210 s
floor, max(CPU work, GPU work). That would be about 2,500–2,600 s per genome. Realistically, with ~90% utilization on both sides, ~230–250 s per
chunk, or ~2,800–3,000 s.

Because the sides are balanced, both kinds of follow-up lever keep their value: CPU work (the
region merge, G-09, ~30–35% of CPU) and GPU work (rescue, extension). Two-pass seeding (G-08)
only moves work from CPU to GPU, so it pays only if one side ends up clearly heavier. Mapped from the code, the waits are:

| # | Barrier or overhead | Where |
|---|---|---|
| B1 | Only chunk 0 aligns while prep runs: `FollowReader` reads chunks in order, so alignment runs at 1/12 of prep's output rate for the first ~330–460 s | `nomics/src/pipeline.rs:139-179, 1077-1092` |
| B2 | GPU seeding runs on the single producer thread, which blocks on `cudaEventSynchronize` before reading the next batch | `mem.rs:1439-1464`, `seed.cu:594` |
| B3 | Seeding and rescue both use the legacy default stream, so they serialize on the GPU; rescue also holds a global mutex across in-flight batches | `mem.rs:400-405`, `ksw.cu:623` |
| B4 | Each tile runs five barriers in sequence (chains → plan → GPU, blocking → post → replay). A lane's CPU idles during its own kernel; only two lanes overlap | `mem.rs:904-1030`, `fused.cu:343-370` |
| B5 | Per-batch setup: the fused engine (HostBuf arenas) is rebuilt per lane per batch and the GPU first-touches new pages; a CUDA stream is created and destroyed per launch; `for_each_parallel` spawns OS threads about 5× per tile and 3× per batch | `mem.rs:743-746, 107`, `fused.cu` |
| B6 | Batch tail: pestat (serial), rescue (one launch per batch, mutex) and pairing run strictly in sequence, and all CPU waits on the rescue kernel | `mem.rs:232-263` |
| B7 | Output: `sort_bam` is single-threaded, and BGZF compression during alignment uses 1 thread by default | `samsort.rs:246`, `bgzf.rs:593-630` |

## What must not change (exactness)

- Batch composition: `read_batch` cuts at `CHUNK_SIZE` = 1e8 bases per prep chunk (BWA
  `-K`).
- pestat is computed over exactly one batch's final regions, and rescue and pairing use that
  batch's `pes`.
- Per-read results are pure functions of the read, its seeds and the index. Scheduling may
  change; values may not.
- Within each output chunk, records go to the sorter in batch order, and within a batch in
  read order.
- **Pair ids are per chunk.** BWA breaks score ties in `mem_sam_pe` / `mem_mark_primary_se`
  with `hash_64(id)`, where `id` is the read counter of one bwa process. sarek runs one bwa
  per chunk, so `first_pair` counts from 0 within each chunk. Each chunk keeps its own counter,
  assigned when `read_batch` cuts the batch, never a global or completion-order counter.
  Single-chunk gates (tier0/tier1 chunks) cannot catch a mistake here; the 30× byte gates can.

With these held, every change below is exact by construction. The existing gates check it:
- tier0 and tier1 chunks on the Mac (CPU paths, and the fused stand-in);
- on the Spark, tier1 and wgs20 18/18, the 30× chunk byte gate, and a full 30× run at 41/41.

## The design

### 1. Persistent workers and engines (removes B5)

- One process-wide CPU pool replaces `for_each_parallel`'s per-call spawning. Use a
  work-stealing scope, rayon-style; the workspace already uses rayon in places.
- One fused engine per GPU lane, kept for the life of the alignment, so its HostBufs and
  stream are reused.
- Streams are created once per lane, per seeding context and per rescue context.

### 2. Asynchronous GPU calls with completion events (B2, B4)

Split each GPU entry point into `launch(...) -> Ticket` (packing, enqueue, record an event)
and `wait(Ticket)`. Synchronous wrappers stay for benches and tests. The C ABI gains
`*_launch` / `*_wait` (or a stream and event out-parameter). Every buffer a ticket touches
stays owned by the ticket until `wait` returns: the Rust type system enforces it, with the
buffers moved into the ticket.

### 3. A three-stage batch pipeline (B2, B3, B6)

```
reader ──► S1 seed (GPU, own stream) ──► S2 align1 (tiles: CPU ⇄ GPU, double-buffered)
       ──► S3 batch tail: pestat → rescue (GPU, own stream) → pairing ──► sorter/writer
```

- **S1** takes batch N+2 while S2 aligns N+1 and S3 finishes N. The reader never blocks on
  the GPU. Seeding gets its own stream, so it no longer serializes with rescue. S1 keeps one
  seeding launch in flight at a time: the seed slot pool is sized for one resident kernel.
- **S2 runs double-buffered lanes.** Each GPU lane holds two tiles. While tile k's kernel
  runs, the lane's workers chain and plan tile k+1, and replay tile k−1. A lane only waits
  when both of its buffers are in flight. Lanes and CPU lanes (the existing split) pull tiles
  from the batch's queue.
- **S3 rescue gets its own stream and no global mutex.** Pairing for the batch overlaps the
  next batch's S2 on the shared pool.
- Queues between stages are bounded. Depth is a knob (default: 1 batch queued per stage, so
  about three batches are live), sized against memory.
- **Order is kept by sequence numbers.** S3 completes batches in order per chunk, using a
  small reorder buffer if S2 finishes out of order.

### 4. Prep feeding without the chunk-0 throttle (B1)

Align every prep chunk whose next batch is ready, not chunk 0 first. Each chunk's batches are
still cut exactly as now (the same `read_batch` over that chunk's stream), and its outputs go
to its own sorter slot, so composition and order per chunk are unchanged. The scheduler takes
the next ready batch from any chunk, keeping round-robin fairness so no chunk starves.

Value is bounded: the throttle costs at most ~60–100 s at 30× (12 × chunk-1 align ≈ the
full-run align phase). Memory is the risk, since each open chunk holds ~1/12 of the genome's
records until it completes. **So the design is a cap of 2–3 open chunks, with peak RSS
measured first**, not a free-for-all.

### 5. Output (B7, small)

`sort_bam` gets a parallel key sort. BGZF deflate threads during alignment follow free CPU
(`NOMICS_MARKDUP_OVERLAP_THREADS`); only measure, since deflate is ~1.4% of CPU here.

### 6. Knobs for the follow-up levers

- **Two-pass GPU seeding (G-08)** becomes an S1 option (`NOMICS_GPU_SEED_TWO_PASS`), so
  "CPU-bound vs GPU-bound" is decided by measurement.
- **The CPU-lane share in S2** stays.

## Memory budget

The 30× peak is 80 GB today; the budget is ≤ 90 GB. Per live batch:
- seeds 1.7 GB (intervals at cap 64) + 0.13 GB of sequence;
- reads 0.3 GB;
- regions ~0.7 GB;
- fused arenas per tile in flight, ~0.2–0.4 GB.

Three live batches plus four in-flight tiles is about 9 GB, close to today's set (2 in
flight + 2 queued + 1 reading). The open question is the chunk-0 fix (section 4). Measure
peak RSS on every full run, and add a sum of live-batch bytes under `NOMICS_TIMING`.

## Delivery, in gated steps (each kept only if exact and not slower)

1. Persistent pool, engines and streams (section 1). Expected small gain; it also removes
   noise.
2. Async GPU tickets + seeding on its own stream + rescue without the mutex (sections 2, 3
   S1/S3).
3. Double-buffered S2 lanes (section 3 S2).
4. Prep feeding across chunks (section 4), measured for memory first.
5. Output sort (section 5) only if the profile shows it.

Each step is measured on 30× chunk 1, fused and default, with **at least 2 runs per config**
(noise is ±15 s on ~292 s). Report chunk align time with both CPU busy (sampler) and GPU busy
(nsys); align time is the target, not either utilization. A full 30× run comes at the end.
Target: alignment at 30× from 3,578 s to ~2,800–3,000 s (floor ~2,550 s).

Delivery is in reviewed increments:
- steps 1–2 are one task: gated, reviewed (GPT-6 Astra on the ticket / ABI change: buffer
  ownership until `wait`, every error path drains) and committed before step 3 starts;
- then step 3, then step 4.

## Risks

- **Out-of-order completion bugs** could change output order. They are caught by the byte
  gates, and sequence numbers plus a reorder buffer make them structural.
- **Buffer lifetime with async GPU work.** Tickets own the buffers, and the C ABI keeps every
  error path draining, as hardened in the kernel reviews.
- **CPU oversubscription** from a shared pool plus GPU-host threads. Use measured worker
  counts per stage, not one guess.
- **More batches in flight raise peak memory.** Bounded queues, and RSS is checked at every
  step.

## Increments 1–2: landed (2026-10-20)

What changed:
- **Pool:** one persistent rayon pool replaces per-call thread spawning.
- **Engines:** fused, two-pass and rescue engines are pooled and reused across batches,
  keeping their HostBufs and streams.
- **Streams:** stream creates per 30× chunk fell from 970 to 7.
- **Async GPU calls:** a `*_launch` → ticket → `nomics_ticket_wait` ABI. A ticket owns its
  engine until it is waited; every error path after a launch drains.
- **Seeding:** a seeder thread between reader and aligners, with seeding on its own stream.
- **Rescue:** its own stream, with no global mutex.

Results:

| | Before | After |
|---|---:|---:|
| 30× chunk 1 align, fused | ~297 s | ~265 s |
| 30× chunk 1 align, default | ~327 s | ~327 s |
| tier1 FASTQ → VCF, fused | 392 s | **344–345 s** |
| tier1 FASTQ → VCF, default | 379 s | 389 s (noise) |
| Host time in GPU event waits (chunk) | 434 s | 242 s |

**Wrong-seed event and fix.**
- **What happened:** one wgs20 fused run on an intermediate build gave 64 contiguous reads
  the seeds of the read 2,048 earlier, which is exactly 16 KiB back in `offs`.
- **Suspected cause:** `offs` and `lens` were the only GPU-read inputs on ordinary 4 KiB heap
  pages. The committed code had the same pattern.
- **Fix:** both now live in pooled HostBufs. The seed gate and bench tools were changed the
  same way.
- **Audit (GPT-6 Astra):** every production GPU buffer is now a HostBuf, with no in-flight
  lifetime, ordering or pool-deadlock defect. PairHMM's post-launch error path now drains.
- **Stress gate (`NOMICS_SEED_VERIFY=1`):**
  - full 30× run: 41/41 identical, 1,068 batches with 0 seed differences, peak 84.8 GB;
  - 4 wgs20 runs: 18/18 each, 144 batches with 0 differences.
- **Unproven root cause.** A HostBuf only *requests* 2 MiB pages. Follow-up: warn when
  per-batch GPU buffers are not actually huge-page backed.

**Huge-page fragility.** A tier1 fused run with the index only 84% on 2 MiB pages took 500 s
instead of 345 s. Every timing run must check the run.log index line, and harnesses should
evict caches first.

## Increment 3: shelved (2026-10-21)

Double-buffered fused lanes were built on branch `shelved/align-double-buffered-lanes`
(0f202cf). Each lane chains and plans tile k+1 and replays tile k−1 while tile k's kernel
runs, with `NOMICS_GPU_LANE_DEPTH` (1 = serial). Exactness held on every gate that ran: Mac
tier0 with ~120 interleaved tiles, the kernel dumps, tier1 and wgs20 fused 18/18, and the 30×
chunk in all four modes. The default-path and seed-verify runs did not happen because the
disk was full.

| 30× chunk 1, fused | Before | After |
|---|---:|---:|
| Align phase | 257.5 / 252.7 s | 257.5 / 254.0 s |
| Lane waits on own kernel | 1,148 waits | 0.4 s in total |
| Host `cudaEventSynchronize` | 233 s | 145 s |
| GPU busy (nsys) | 70.5% | 74.4% |
| System memory peak | 88.4 GB | 90.5–91.3 GB (> 90 GB budget) |

**Why it doesn't help:** with 2 batches in flight, chunk time is set by batch latency. Per
batch: align1 3.0 s, then rescue 1.27 s, then pairing 0.76 s, in series. The removed waits
were not on that path; the batch tail is.

**Next:** overlap the batch tail with the next batch (a third batch in flight). That needs
memory headroom first. Candidates:
- `SEED_POOL_KEEP` 6 → 3: −3.3 GB;
- compact the cap-64 seed interval buffers (1.7 GB per batch);
- drop per-tile arenas earlier.

The CPU-side region merge (G-09) is independent and still pays.

## Region merge, exact CPU speedup (G-09, 2026-10-21)

On real reads `sort_dedup_patch` almost never patches. The cost was the two sorts moving whole
96-byte regions, plus an O(n²) dedup loop. Mate rescue makes ~46% of the calls, at ~90
regions per call.

Changes:
- sort packed integer tags, falling back to the literal klib introsort when any key ties;
- a provably safe early break in the dedup loop;
- skip the identical-hit pass when nothing ties;
- a byte-wise 2-bit `get_seq`;
- a score-only `global2_score` for the patch path;
- reused buffers in `global2` and `gen_cigar2`.

The added `score > 0` merge condition matches bwamem.c (docs/spec/bwa.md:338); the old port
lacked it, and it never fired in any test.

Evidence:
- **Real-genome replay:** `merge_bench` on 1,982,423 calls (162M regions) from 500k real 30×
  pairs, fused path: 0 mismatches against the dump and against a frozen copy of the old code,
  2.84× faster.
- **CIGARs:** `global2` CIGARs checked against BWA's C (`ksw_vs_c`, including the 50k-case
  campaign).
- **Review:** GPT-6 Astra found one regression, a `global2` panic on an empty-query
  narrow-band row. It is fixed and pinned by a test.
- **Gates:** kernel dumps 0 mismatches; tier1 and wgs20 18/18 on both paths; the 30× chunk
  byte-identical in all modes.

| | Before | After |
|---|---:|---:|
| 30× chunk 1 align, fused | ~255 s | 231.2 / 231.1 s |
| 30× chunk 1 align, default | ~321 s | 291.3 / 292.9 s |
| tier1 FASTQ → VCF, default / fused | 389 / ~345 s | 372 / 343 s |
| wgs20 FASTQ → VCF, default / fused | 265 / 211 s | 255 / 203 s |

## Memory headroom and a third batch in flight (2026-10-21)

Measured on 30× chunk 1 (align-only) with `[mem]` lines under `NOMICS_TIMING` (RssAnon,
live HostBuf bytes, seed sets, chunk record bytes) and the system sampler.

Where the memory went (fused, 2 batches in flight, before):

| Item | Size |
|---|---:|
| Index | 21.7 GB |
| GPU seed buffer sets: 5 live × 1.75 GB (2 aligning, 2 queued, 1 seeding) | 8.8 GB |
| Rescue engines: per engine 4.4 GB (rows 2.1, targets 1.1, scratch 1.0), 2–3 engines | 9–13 GB |
| Fused engines (4 lanes) | ~2.5 GB |
| Chunk records awaiting the sorter (59M records) | 21.2 GB used, 30.6 GB allocated; the 15 → 30 GB growth copied 15 GB |
| `sort_bam` | a second full copy of the chunk |
| Per live batch: reads 0.26 GB, regions 0.97 GB | |

The process peak was at the chunk sort, not in the batch loop: RssAnon 87.4 GB, system 93.3 GB
(inflight 3: 100.4 GB).

Changes, all exact by construction (same records in the same order):
- **Compact GPU seeds.** After the seed kernel, the cap-64 result buffer is copied into a
  CPU-only `GpuSeeds` (intervals + per-read spans, 0.25 GB per batch, ~25 ms on the seeder
  thread) and the HostBuf set returns to the pool at once. One set is live instead of five;
  `SEED_POOL_KEEP` 6 → 2. `GpuSeeds::get` returns a slice, so the per-read clone is gone.
  A batch's seeds are dropped right after align1 (rescue and pairing only need "GPU-seeded").
- **Chunk records by reference.** A batch's records are regrouped by reference id on its own
  thread (`BamRecords::into_bins`) into one mapping of their own (`nomics_mem::PageBuf`), each
  reference page-aligned with its sort keys in front; the chunk keeps every reference's bins
  in batch order (`ChunkRecords`): no growing 20 GB buffer. `sort_chunk` sorts one reference
  at a time and returns its pages at once (`MADV_DONTNEED`), so the chunk is held about once.
  Ordering: a key's top bits are the reference id, and within a reference `(key, batch,
  record)` orders exactly as `(key, record ordinal)` did. Record memory bypasses the
  allocator: batches are built on their own threads and freed on the sorter thread, and with
  heap bins mimalloc kept such frees resident (+16 GB during one chunk sort).
- **Rescue buffers** grow to the need + 1/8 instead of the next power of two (avoids a 4 GB
  rows buffer when a batch crosses 2 GB).
- **After alignment** the seed, fused, two-pass and rescue pools are unmapped
  (`release_gpu_buffers`), so markdup/BQSR/HC no longer carry ~10–15 GB of idle engines.
- **Whole-genome chunk BAMs are streamed.** Full 30× runs showed what align-only runs do
  not: each sorted chunk (~21 GB) was held as a `RecordSet` while the mapped BAM writer
  compressed it on 1 thread. The first full run with the changes above aligned faster (3,942
  s vs 4,117 s), the writer fell behind, sorted chunks queued (drained 206 s after alignment
  instead of ~10 s) and the system peaked at 103.6 GB (process 95.4 GB). With 2 writer
  threads it drained in 9.9 s and peaked at 80.1 GB (process 73.2 GB), but the fused path at
  3 batches in flight still reached 104.0 GB (process 95.8 GB). Now, for spilled (whole-genome)
  runs, the sorter writes each chunk's BAM itself (`align_chunks_bam_files`,
  `samsort::write_sorted_chunk`): records stream from the chunk's pages into the BGZF writer,
  each reference freed once written, so no sorted copy exists. BGZF blocks are cut by byte
  count, so the file is the one `RecordSet::write` gave (unit-tested byte for byte). The
  writer uses 2 threads while later chunks align (`NOMICS_MARKDUP_OVERLAP_THREADS`, the same
  CPU work done in half the time) and all threads for the last chunk; exomes keep the
  in-memory path, whose writer also widens by one share per chunk waiting behind.

30× chunk 1 align-only, timing off. "Before" (fa3e154) rows are controls run in the same
session, interleaved with the "after" runs; they came out 5–10% slower than the region-merge
numbers above for the same commit (fused 231 s, default 291 s), so compare within this table.
System peak is MemTotal − MemAvailable; cores and GPU are medians/means over the whole process
(index load, sort and write included), not the align window alone.

| Build, mode, batches in flight | Align (s) | System peak | RssAnon peak | Cores (median) | GPU |
|---|---:|---:|---:|---:|---:|
| before, default, 2 | 324.8 / 320.9 | 90.9 / 90.2 GB | 86.9 GB | 19.9 | 28% |
| after, default, 2 | 294.9 / 292.9 (296.2 / 308.0) | 64.5 / 63.9 GB | 59.2 GB | 19.9 | 30% |
| after, default, 3 | 306.2 / 303.7 | 73.3 / 73.0 GB | 68.2 GB | 20.0 | 30% |
| after, default, 4 | 292.7 | 79.4 GB | 73.9 GB | 20.0 | 30% |
| before, fused, 2 | 258.1 / 244.1 | 93.6 / 91.7 GB | 89.3 GB | 13.0 | 52% |
| before, fused, 3 | 229.0 | 100.4 GB | 97.2 GB | 16.4 | 57% |
| after, fused, 2 | 233.6 / 235.9 | 70.1 / 70.1 GB | 64.3 GB | 13.7 | 54% |
| after, fused, 3 | 210.5 / 214.3 | 79.5 / 78.7 GB | 73.6 GB | 17.1 | 60% |
| after, fused, 4 | 209.4 | 87.4 GB | 81.1 GB | 18.6 | 60% |

The "after" rows were measured before records moved to one mapping per batch and before the
writer change; with them, a 3-chunk default align-only run peaked at 62.8 GB system and the
gate's fused chunk took 209.5 s, in line with the table.

Full 30× runs with the final code (streamed chunk BAMs, rescue headroom 1/4), sampler at 1 s:

| Run | Gates | FASTQ → VCF | Align | Process peak | System peak |
|---|---|---:|---:|---:|---:|
| default, 2 in flight (2026-10-20, before) | 41/41 | 5,781 s | 4,117 s | 78.8 GB | — |
| default, 2 in flight | **41/41** | **5,379 s** | 3,706 s | 72.3 GB (align phase 60.8) | 79.2 GB (markdup start; 57–67 GB while aligning) |
| fused, 2 in flight (2026-10-20, before) | 41/41 | 5,240 s | 3,578 s | 80.1 GB | — |
| fused, 3 in flight | **41/41** | **4,374 s** | 2,672 s | 84.3 GB | **92.3 GB** |

Both final runs flagged small GPU buffers partly on 4 KiB pages (a 22 MiB rescue buffer;
64–256 MiB fused engine buffers); timing-suspect by the harness rule, so read the times as
indicative.

**Decision: `NOMICS_BATCH_INFLIGHT` stays 2.** A third batch makes the fused path much faster
(30× chunk −10%, whole genome 4,374 s), but a fused whole genome then peaks at 92.3 GB system
(84.3 GB process), over the 90 GB budget; its memory also creeps up across chunks (80 → 90 GB)
rather than cycling as the default path does. On the CPU-extension path a third batch gains
nothing (CPU-saturated), and exomes do not need it (tier1 fused 306.0 s at 2 vs 313–319 s at
3). `NOMICS_BATCH_INFLIGHT=3` is the switch to flip for fused whole genomes once ~3 GB more is
freed. With the final default the fused gate gave tier1/wgs20 fused 18/18 (306.0 / 203.5 s)
and the 30× chunk identical in all four modes (cpu 297.7, two-pass 254.0, fused 230.2 s).

Gates with the final code: Mac tier0 0001/0002/0012 default and fused stand-in identical to
sarek; `chunk_sort_matches_whole_sort` (ties across batches; streamed file byte-identical to
`RecordSet::write`); Spark cuda tests; tier1 and wgs20 default 18/18; the two full 30× runs
above 41/41; `fused_spark_gate.sh` (kernel dumps, tier1/wgs20 fused 18/18, the 30× chunk in
all four modes) at 3 in flight and again with the final default of 2.

**Next memory levers** (for 3 in flight on whole genomes): the rescue engines, now the largest
item after the index (one per concurrent batch tail, ~5 GB each at 1/4 headroom: rows 2.1,
targets 1.1, scratch 1.0 GB), which could run in bounded sub-launches; the seeded-batch queue
(`sync_channel(inflight)`, 0.5 GB per queued batch; depth 1 suffices); fused engine buffers
still grow to powers of two; and the fused path's slow growth across chunks.

## Full 30× after this week's work (2026-10-22, d6095cc)

| FASTQ → VCF, 30×, Spark | Bake-off (10-17) | 10-20 | **10-22** |
|---|---:|---:|---:|
| Fused GPU extension | — | 5,240 s | **4,572 s** (align 2,875 s, peak 72.5 GB) |
| Default (CPU extension) | 6,349 s | 5,781 s | **5,379 s** (peak 72.3 GB) |

Both runs are 41/41 identical. The fused run warmed 7 GPU streams and had 0 fallbacks. It was
flagged only because a few small GPU buffers were partly on 4 KiB pages.

## Re-baseline on the Linux 7.0 kernel (2026-10-22, b80201d)

The Spark moved to kernel 7.0.0-1019-nvidia; driver 580.178.04 and CUDA 13.0 are unchanged.
There is no regression:

| | Before | After |
|---|---:|---:|
| Kernel dumps | 0 mismatches | 0 mismatches |
| tier1 fused / default | 308 / 359 s | 307 / 358 s |
| wgs20 fused / default | 207 / 252 s | 207 / 254 s |
| 30× chunk align, fused / default | 232 / 287–294 s | 230–232 / 292–310 s |
| Full 30× fused | 4,572 s | **4,530 s**, 41/41 |
| Full 30× default | 5,379 s | **5,358 s**, 41/41 |

The default full run was TIMING SUSPECT. It started right after a 220 GB archive to the NAS,
and the index reached only 85% on 2 MiB pages despite eviction plus `compact_memory`.
Compaction does not undo that much churn; keep archives away from timed runs.

MPAM: the kernel has MPAM and resctrl built in, but the GB10's firmware lists no MSC devices
("No MSC devices found in firmware"), so there are still no memory-bandwidth monitors.

## A third batch in flight for fused whole genomes (2026-10-23)

`NOMICS_BATCH_INFLIGHT` now defaults to 3 for the fused extension on whole genomes (chunk BAMs
streamed to files: spilled `nomics run`, and `--resume-from align`), 2 everywhere else; the
variable still overrides. `gpu_warm_up` takes the same flag, so a fused whole genome warms 10
streams (seeding, 3 rescue, 6 fused lanes). Output is byte-identical: the changes below only
move or drop buffers, and every gate passed.

Where the memory went (30× chunks 1-3, fused, 3 in flight, align-only, `[mem]` lines):
HostBuf 42.1 GB before, 27.4 GB after (the index is 21.7 of it). Changes:
- **Rescue buffers the kernel never read.** The default (wavefront) mate-rescue kernel keeps
  rows in shared memory and ignores `row_maxima`; the scratch buffer is only for the striped
  kernel (`NOMICS_KSW_STRIPED=1`). Each rescue engine still mapped both, ~3.9 GB of its ~5 GB
  at 1/4 headroom, one engine per batch in flight. They are now mapped only for the striped
  kernel, and the wavefront launch gets NULL (both shims ignore it; noted in `nomics_gpu.h`).
- **Bounded rescue launches** (`GpuKswAlign2::compute_bounded`, `NOMICS_RESCUE_LAUNCH_MB`):
  consecutive slices of the calls, results in input order (unit-tested against one launch on
  CUDA and Metal; tier0 Metal BAMs identical at 1 MiB). At 256 MiB per launch (4 per 30×
  batch) the chunk was ~10% slower (224.8 / 222.0 s vs 203.9 s one launch), because rescue is
  on each batch's serial tail and each slice is packed only after the last one is waited. The
  default is 2 GiB, so a 30× batch (~1 GB of bases, ~750k calls) stays one launch and only an
  outsized batch is split.
- **Seeded-batch queue depth 1** (was the batch count): ~0.5 GB per queued batch; depth 3 was
  no faster (229.5 s).
- **Fused (and two-pass) engine buffers** grow to the need + 1/4, capped at the old power of
  two; no time change (power of two: 224.6 s).

30× chunk 1, fused, timing off (`wgs30_align_gate.sh`, PASS each time):

| Build, batches in flight | Align (s) | System peak | Process peak |
|---|---:|---:|---:|
| master, 2 (time: re-baseline; memory: 10-21 table) | 230–232 | 70 GB | 64 GB |
| master, 3 (10-21 table) | 210.5 / 214.3 | 79 GB | 74 GB |
| new, 2 (with 256 MiB rescue launches) | 241.2 / 241.4 | 59 GB | 56 GB |
| new, 3, 256 MiB rescue launches | 221.4 / 225.4 / 224.8 / 222.0 | 63 GB | 60 GB |
| new, 3, one rescue launch (final) | 203.9; gate 202.0 | 66 GB | 63 GB |

**The cross-chunk creep was the mapped-BAM writer, not a leak.** The 3-chunk align-only runs
cycle (process back to 28.6 GB after each chunk, idle engines and HostBuf flat), because
align-only writes each chunk on all threads. A full run writes chunk BAMs on 2 threads while
later chunks align, ~190 s per 30× chunk (~270 s while prep still runs). At 2 batches a chunk
aligns in ~230 s and the writer keeps up; at 3 (~200 s) it never caught up: in the full fused
run each chunk BAM was created the moment the previous one was closed, sorted chunks waited
longer and longer for it, and the process floor rose ~1.5 GB per chunk (60 → 70 GB over chunks
3-11). The writer now gets one thread per batch in flight on whole genomes (3 for fused, still
2 for the default path; `NOMICS_MARKDUP_OVERLAP_THREADS` overrides). BGZF blocks are cut by
byte count, so the files do not depend on the thread count. With it, chunk BAMs take
~146–217 s with idle gaps between them, and the process floor stays flat (52–55 GB) through
all 12 chunks; the system peak (78.2 GB) moved to the end of alignment and markdup's start,
where the default path peaks too.

Full 30× (sampler at 1 s; system = MemTotal − MemAvailable):

| Run | Gates | FASTQ → VCF | Align | Process peak | System peak |
|---|---|---:|---:|---:|---:|
| fused, 2 in flight (master, 10-22) | 41/41 | 4,530 s | — | 72.6 GB | — |
| fused, 3 in flight (10-21, before) | 41/41 | 4,374 s | 2,672 s | 84.3 GB | 92.3 GB |
| fused, 3, writer on 2 threads | 41/41 | 4,331 s | 2,671 s | 74.4 GB | 81.9 GB |
| fused, 3, writer on 3 threads (final) | **41/41** | **4,262 s** | 2,595 s | 72.6 GB (align phase 60.8) | **78.2 GB** (67.5 while aligning) |
| default, 2 in flight (final) | 41/41 | 5,343 s | 3,676 s | 72.4 GB | 78.5 GB |

Both fused runs were flagged only for a few small fused engine buffers partly on 4 KiB pages,
as the 4,530 s baseline was; the default run was not flagged. Other gates with the final code:
Mac tier0 0001/0002/0012 default and fused stand-in PASS against sarek (Metal GPU seeding and
rescue too); Spark cuda tests 70/70; `fused_spark_gate.sh` (kernel dumps 0 mismatches; tier1 /
wgs20 fused 18/18 in 306.5 / 204.7 s; the 30× chunk cpu 294.7, two-pass 254.9, fused 202.0 s);
tier1 / wgs20 default 18/18 (358.5 / 254.8 s).

## Clean full 30× after a reboot (2026-10-23, a9c55a1)

Freshly rebooted Spark (uptime 45 s at start), kernel 7.0.0-1019-nvidia, with fast BAM
deflate and 3 batches in flight for fused whole genomes. Both runs are 41/41 identical, with
0 GPU fallbacks.

| Stage | Fused | Default |
|---|---:|---:|
| prep (overlaps alignment) | 329 s | 394 s |
| align | 2,495 s | 3,598 s |
| markdup | 222 s | 220 s |
| BQSR | 319 s | 316 s |
| HaplotypeCaller | 944 s | 935 s |
| **FASTQ → VCF** | **4,003 s (1 h 07 min)** | **5,086 s (1 h 25 min)** |
| Peak RSS | 72.5 GB | 72.3 GB |

The default run is clean. The fused run carries only the "GPU buffers on small pages" flag,
for a few small engine buffers; the index was fully on 2 MiB pages.

**Comparison on the same box:**
- sarek: 41,044 s, so nomics is 10.3× faster.
- The bake-off's nomics time was 6,349 s, so today's run is 37% faster than that.
