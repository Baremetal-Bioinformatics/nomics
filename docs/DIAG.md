# `NOMICS_DIAG`: memory-corruption diagnostics

The mode exists to catch the 2026-10-24 panic at its source. One full 30x run with
`NOMICS_GPU_EXTEND=0` hit `index out of bounds: the len is 3366 but the index is
18446744073709551615` in `mem_pair`. That is a primary region with `rid == -1`, which no
code path produces.

The leading hypothesis is a stray GPU or DMA write into live host memory. HostBufs had no
guard pages, so such a write would land silently in a neighbouring mapping.

The mode is cheap enough to leave on for full runs. It has two settings, and the first line
of the run log says which is active:

| `NOMICS_DIAG` | What runs |
|---|---|
| unset, empty or `0` | Nothing below. The HostBuf layout is the original one, there is no registry and no check, and output is byte-identical. |
| `regions` | Region invariants, the checked `mem_pair` report and fatal GPU errors only, on the original HostBuf layout: no guard pages, canaries, registry or fault handler. |
| `1` (anything else) | Everything: the `regions` checks plus guard pages, canaries, canary sweeps, the registry and the fault handler. |

**Why `regions` exists.** Guard pages move buffers apart, so a bug that depends on which
mappings are adjacent can vanish in full mode. A clean `NOMICS_DIAG=1` run cannot rule
that out, but a `regions` run keeps the normal layout. In `regions` mode the region
reports list no HostBufs, since there is no registry.

tier0 digests are the same in all three settings.

## What is checked, and where

**Guard pages and canaries around every HostBuf** (`crates/nomics-mem/src/lib.rs`,
`HostBuf::with_guards`). The mapping layout is:

```
[PROT_NONE head guard ...][head canary page][data: 2 MiB aligned, len][slack canary to page end][PROT_NONE tail guard ...]
```

- **Data range.** The start keeps its 2 MiB alignment, and `MADV_HUGEPAGE` is still applied
  before first touch.
- **Guards.** The guards are the over-mapped slack, at least one page on each side
  (4 KiB on the Spark, 16 KiB on the Mac).
- **Huge-page cost.** Only the last, partial 2 MiB of each buffer loses its huge-page
  backing, because the tail guard splits it. Buffers whose length is a multiple of 2 MiB
  lose nothing. `touch_huge` leaves that tail out of its "<90% on 2 MiB pages" warning, so
  the gate's TIMING SUSPECT grep does not fire on it.
- **Registry.** Every guarded buffer is in a lock-free registry (`crates/nomics-mem/src/diag.rs`)
  with its label (the `touch_huge` label, or the engine field name).

**CPU faults.** A CPU access to a guard page delivers SIGSEGV (SIGBUS on macOS). The
handler (`diag.rs`, `on_fault`) acts only on kernel-reported faults inside a registered
buffer's guard or canary pages. It prints one line, then the default action runs and the
process dumps core:

```
==== NOMICS_DIAG: SIGSEGV at 0x... in the tail guard page of HostBuf 'GPU seed intervals' (data 0x..., len N, offset +M) ====
```

Every other SIGSEGV/SIGBUS goes to the previous handler, called directly, and ours stays
installed. That includes faults outside our mappings (for example a recoverable one in the
CUDA driver) and software-sent signals. If the previous disposition is `SIG_DFL`, the signal
is re-raised with it; `SIG_IGN` returns.

**Canary sweeps after every GPU wait.** `nomics_mem::check_all_canaries` checks every live
HostBuf's head canary page and tail slack. It runs:
- after every ticket wait, in `Pending::wait` (`crates/nomics-gpu/src/lib.rs`, ~line 1208);
- after an unwaited ticket is dropped, in `Pending::drop` (~1218).

That covers the seed collect (both passes), mate rescue, fused extension and two-pass
extension. On damage the process aborts. The report gives the buffer name, the offset from
the data start (negative means the head canary), the damaged byte count, up to 32
found/expected bytes, and every live HostBuf's address ranges.

**GPU errors are fatal.** With the mode on, these all print the error, the failed launch's
buffers by field (address ranges and canary state) and every live HostBuf, then abort. No
CPU fallback runs.
- **Ticket waits:**
  - mate rescue: `GpuKswTicket::wait`, lib.rs ~439;
  - fused extension: `GpuFusedTicket::wait`, fused.rs ~160;
  - two-pass extension: extend.rs ~273 and ~331;
  - seed collect: `GpuSeedBuffers::collect`, lib.rs ~743.
- **Fallback sites** (backstops, which also catch launch and allocation errors) in
  `crates/nomics-align/src/mem.rs`:
  - rescue: ~591;
  - two-pass: ~742;
  - fused: ~1180;
  - seed: ~1567;
  - seed retry: ~1636.
- **GPU warm-up:** `crates/nomics/src/pipeline.rs` ~1339.
- **Exception:** `GpuError::Unavailable` (a kernel the backend lacks) is not an error.

Without the mode, these paths keep their old behaviour: they fall back to the CPU and print
`GPU FALLBACK` at the end.

**Region invariants** (`crates/nomics-align/src/diag.rs`). Every region must satisfy:
- `0 <= rid < n_seqs`;
- `0 <= rb <= re <= 2*l_pac`;
- `0 <= qb <= qe <= l_query`.

`qe == qb` is allowed: it is merge's "excluded" marker. The checks run at these points:

| Check | Where |
|---|---|
| after `chain2aln` (extension; CPU, two-pass and fused paths all end in `align1_extend`) | `extend.rs` ~1370 |
| after `sort_dedup_patch` (region merge) | `extend.rs` ~1372 |
| whole batch, after align1, before pestat and GPU mate rescue | `mem.rs` ~304 |
| whole batch, after the GPU mate-rescue kernel, before pairing | `mem.rs` ~316 |
| per pair, after `matesw` (rescue merge), before `mark_primary_se` | `pair.rs` ~738 |
| per pair, after `mark_primary_se`, before `mem_pair` | `pair.rs` ~749 |
| `mem_pair`'s contig lookup is now a checked access (always, not only in diag mode) | `pair.rs` ~589 |

**Region failure report.** On the first failure the process prints a report, then aborts
(a core dump; no panic in one thread). The report covers:
- the stage;
- the read name;
- its index in the batch, and its pair id and mate;
- the batch, as `chunk C batch B (as for NOMICS_ALIGN_SKIP_BATCHES=B)`. `B` is the reader's
  ordinal within the chunk, counting skipped batches, so a replay with
  `NOMICS_ALIGN_SKIP_BATCHES=B` aligns that batch first. `C` is the stream index in this
  run: with `NOMICS_ALIGN_CHUNKS=8` alone, stream 0 is chunk 8. The report also gives the
  read group and the batch's first pair id;
- the bad region and its neighbours, with every field;
- the region Vec's buffer address and capacity;
- the eight HostBufs nearest that buffer.

Per-read context comes from `diag::at_read` in the align1 loops, a name copy into a
thread-local. Batch context comes from `diag::batch` in `process_batch_seeded_to` and from
`diag::set_chunk`, which records the reader's chunk and ordinal, where batch threads start
(`mem.rs` ~2209).

At the end of the alignment the run log carries a line like:

```
NOMICS_DIAG: <n> canary sweeps after GPU waits found no damage; <k> guarded HostBufs live
```

## How a GPU fault on a guard page surfaces (expected, not yet confirmed on GB10)

On GB10 the GPU reaches host memory through ATS, which walks the CPU page tables. A kernel
access to a `PROT_NONE` page should therefore fault instead of landing in a neighbouring
mapping:
1. **The fault is sticky.** The context gets a sticky error, most likely
   `cudaErrorIllegalAddress` (700) or `cudaErrorLaunchFailure` (719).
2. **The ticket wait reports it.** The error comes back from `nomics_ticket_wait`
   (`cudaEventSynchronize`) of that launch.
3. **The canary sweep runs first.** It reports any short overrun that stopped before the
   guard.
4. **The ticket wait aborts.** It prints `GPU error in <launch> wait: CUDA error 700 (...)`,
   that launch's buffers and every live HostBuf.
5. **The kernel log has the address.** It should hold an SMMU/ATS or NVRM Xid line with the
   faulting address (`sudo dmesg | tail`). Match it against the printed ranges.

If GB10's ATS turns out not to honour `PROT_NONE` for GPU writes, the guard pages catch
only CPU accesses. In that case the canaries and region checks remain, and memcheck (below)
is the GPU-side tool.

## Limits

- **Writes past `len` only.** Guards and canaries see only writes past a buffer's `len`.
  Engine buffers are sized at the need plus 1/4 and kept across batches, and kernels write
  only the current launch's `n` elements. An overrun of the used part therefore lands in
  unused capacity, which this mode cannot see. For that case the region invariants are the
  only catch.
- **Allocator memory is not guarded.** Region Vecs live in allocator (mimalloc) memory. A
  stray write that skips over a guard into allocator memory is caught only by its effect:
  the region checks.
- **PairHMM is not swept.** PairHMM waits do not go through `Pending`, so no canary sweep
  follows them. That does not matter for alignment.
- **Some degradations now abort.** In diag mode, a GPU warm-up failure (memory too
  fragmented for the driver's stream buffers) and a seed retry buffer that cannot be
  allocated abort instead of degrading. An abort right after start is that, not the bug.
  Evict caches, compact memory and rerun.

## Commands (Spark)

Build: `cargo build --release --features cuda -p nomics`.

**(a) Full 30x run on the default (fused) path with diagnostics.** This is the full-run
script with the env variable set:

```
cd ~/omics/new_omics && NOMICS_DIAG=1 nohup flock "$NOMICS_LOCK" \
  timeout 4h scripts/full_exome.sh wgs30 ~/omics/runs/wgs30-diag-$(date +%Y%m%d-%H%M) \
  > ~/omics/logs/wgs30-diag.log 2>&1 &
```

To reproduce the panic's configuration, add `NOMICS_GPU_EXTEND=0` (CPU extension, GPU
seeding and rescue). Afterwards, check the following:
- `grep -n "NOMICS_DIAG" <out>/run.log`: the first line says the mode was on, and the
  summary line should report no damage;
- any `==== NOMICS_DIAG FATAL` block;
- `dmesg` on an abort.

A core dump needs `ulimit -c unlimited` in the shell that starts the job.

**(b) Chunk-8 replay under compute-sanitizer memcheck.** Run it with diagnostics, on the
panic's CPU-extension path:

```
cd ~/omics/new_omics && CHUNKS=8 NOMICS_DIAG=1 NOMICS_GPU_EXTEND=0 NOMICS_ALIGN_MAX_BATCHES=4 \
  NOMICS_WRAP="compute-sanitizer --tool memcheck --error-exitcode 99 --print-limit 100 --log-file $HOME/omics/logs/memcheck-chunk8.%p.log" \
  nohup flock "$NOMICS_LOCK" timeout 2h \
  scripts/wgs30_align_gate.sh ~/omics/runs/memcheck-chunk8-$(date +%Y%m%d-%H%M) > ~/omics/logs/memcheck-chunk8.log 2>&1 &
```

- **Batch knobs.** They apply to the chunked align path (2 or more batches in flight, the
  default) and work without `NOMICS_DIAG` too.
  - `NOMICS_ALIGN_MAX_BATCHES=N` aligns N batches per chunk, then ends the chunk as at EOF.
  - `NOMICS_ALIGN_SKIP_BATCHES=K` reads and numbers the first K batches but does not align
    them, so later batches get the same pair ids (and so the same output) as in a full run.
    A 30x chunk has ~87 batches of ~333k pairs. To scan the chunk in under 2 h per job, step
    `K` (e.g. 0, 4, 8, ...).
- **Calibrate first.** Start with `N=4` and read the per-batch `[timing]` lines in
  `run.log`, then size `N` to the ~2 h budget. memcheck's slowdown on these kernels has not
  been measured.
- **Gates.** The mapped-BAM gates FAIL by design on a partial replay. The result is
  memcheck's log and exit code (99 on any error) plus any `NOMICS_DIAG` block in `run.log`.
- **Smoke test first.** Before the long job, check that memcheck accepts in-place accesses
  to plain mmap'd host memory, which our kernels use. Under ATS/HMM it should, but this is
  unverified here:
  ```
  cargo test --release --features cuda -p nomics-gpu --no-run   # prints the test binary path
  compute-sanitizer --tool memcheck --error-exitcode 99 <that binary> gpu_sums_host_buffer_in_place
  ```
  If that reports invalid accesses to host addresses, memcheck cannot judge these kernels;
  use (a) instead. If init prints benign API-error reports (e.g. from attribute or occupancy
  queries), add `--report-api-errors no`.

## Status: not reproduced on 082fd26 (2026-11-03, parked)

The one `rid -1` panic (`default1`, 2026-10-25, 9c64a30, CPU extension, after chunk 8) has
not recurred. Replaying its chunk passed (2026-10-26, memcheck 0 errors), and so did every
full 30x run since, including these back-to-back runs on 082fd26 with
`NOMICS_DIAG=regions NOMICS_GPU_EXTEND=0`, 42/42 gates each, no kernel lines:

| Runs | Condition | Wall per run | Pages migrated/s | Compaction stalls/run |
|---|---|---|---|---|
| 4 | page cache full (114 of 121 GB), box up 4.4 days | 89-90 min | ~19.6k | ~435k |
| 4 + 2 | the same, plus a compaction stressor (`MADV_HUGEPAGE` faults over a fragmented 4 KiB ballast, at most 12 GB) | 95.6 min | ~47.9k | ~930k |

What is known:

- Every path that makes a `Region` drops `rid < 0` (chaining, as BWA's `mem_chain` does)
  or copies an id already checked (extension from the chain, mate rescue gated on
  `a.rid == fetch_seq`'s rid). A `-1` at `mem_pair` therefore needs a write from outside
  that logic, which is consistent with the stray-write hypothesis above.
- `vm.compact_unevictable_allowed` is 1 on the Spark, so compaction migrates even
  `mlock`ed pages. Buffers the GPU writes into are not locked. A GPU write through a
  stale translation during migration would land outside every guard. Raising the
  migration rate 2.4x did not reproduce the panic, which weighs against migration alone.
- At a rate of about 1 in 12 runs, 6 clean runs happen about 59% of the time by chance,
  so this is not evidence of a fix.

If it recurs: keep the run directory and the kernel log from `journalctl -k`, note the
commit and `NOMICS_*` settings, and rerun its chunk with `NOMICS_ALIGN_SKIP_BATCHES`.

The migration test needs two things together: `sudo sysctl vm.compact_unevictable_allowed=0`
(lost at reboot) and `NOMICS_HOSTBUF_LOCK=1`, which locks every `HostBuf` in RAM as its pages
are touched (`mlock2` with `MLOCK_ONFAULT`), so the kernel never moves a page the GPU reads or
writes in place. The soak runs it with a `lock` word in its tag (`[soak: 6 extend0=6 lock]`);
a run is unclean if the sysctl is not 0 or a lock fails. The process also needs an unlimited
locked-memory limit (`ulimit -l`): the GitHub runner's systemd service has 8 MiB on the Spark
unless its unit sets `LimitMEMLOCK=infinity`, and every lock then fails with ENOMEM (the first
locked 30x soak run, 2026-11-11, came back unclean for that reason after 79 minutes; the
script now refuses to start instead). What the outcomes mean: a `rid -1`
panic with the buffers locked rules migration out. Clean runs do not rule it in, because the
panic was already rare without the lock (not seen again in the 30 or so full runs since).

Result, 2026-11-12: six locked 30x runs with CPU extension, started from an ssh shell
(`ulimit -l unlimited`, sysctl 0, commit 2f42f07), all clean: 41/41 gates each, no panic, no
failed lock, 25.4 GiB locked at peak, 5,346 to 5,394 s per run including the digest. So the
panic is parked as not reproduced, cause unknown: about 36 full runs since it, locked or not,
without a recurrence. `NOMICS_DIAG=regions` and the lock switch stay for the next one.
