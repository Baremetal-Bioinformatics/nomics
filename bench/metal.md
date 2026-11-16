# Apple Silicon / Metal: `nomics run` on an M4 Max (2026-10-25/26)

24 h on an EC2 `mac-m4max.metal` (Apple M4 Max, 16 cores, 40-core GPU, 128 GiB,
macOS **15.8**), branch `metal/phase-b`. Plan, log and caveats:
`docs/strategy/metal-sprint.md`. Every run below is gated against the sarek 3.10.0
reference run by `nomics dev digest` and the committed manifests: 41 normalized artifact
digests (24 chunk FASTQs, 12 chunk BAMs, md.bam, metrics, recal table, recal.bam, VCF).
No run failed a gate.

## What runs on the GPU (Metal, `--features metal`)

| Kernel | Port of | Exactness evidence |
|---|---|---|
| seed search (`mem_collect_intv`, two-pass 64/128 lists) | `seed.cu` | `seed_gate`: all 991,118 tier0 reads identical to the CPU (cap64, two-pass, k = 12 table, plain lists) |
| mate rescue (`ksw_align2` wavefront) | `ksw.cu` | `ksw_bench`: 0 / 62,000 (exome) and 0 / 70,248 (wgs20) mismatches; synthetic windows up to 4 kb and reads up to 512 bp, all five kernel widths |
| fused seed extension (left + right `ksw_extend2`, band retries) | `fused.cu` + `extend_v3.cuh` | `fused_bench`: 0 mismatches on the tier0 dump (200,888 jobs), the band-retry fixture and the 30× dump (5,434,925 jobs, 7.84 M calls) |
| PairHMM | — | stays on the CPU: FP64, and Metal has no double |

Host memory is read in place (`newBufferWithBytesNoCopy`). Long-lived buffers (the
3.2 GB BWT, seed and engine buffers) are registered once (`nomics_gpu::GpuHostBuf`). The
kit carries the kernels prebuilt (minimum macOS 15.0), because EC2's macOS images have no
Metal compiler.

## Tier1 exome (HG002 WES, 40.9 M pairs), FASTQ → VCF, 16 threads

| Configuration | total | prep | align | markdup | BQSR | HC | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| CPU arm | 693.2–706.6 s | 19.5–36.2 | 582 | 33 | 20.4 | 21 | 110 GB |
| Metal seeds + rescue | 379.2–380.6 s | 19.7 | 285 | 33 | 20.3 | 21 | 111–115 GB |
| + BWT registered, two-pass seeding | 372.0 s | 19.8 | 277 | 33 | 20.1 | 23 | 115 GB |
| **+ fused extension, 6 GPU lanes (final)** | **320.0 / 333.0 s** | 29 / 43 | **216 / 215** | 33 | 20.3 | 21 | 125 GB |

The final configuration is **2.1–2.2× the CPU arm (align 2.7×)**. Totals vary in prep,
because at 122–127 GB peak the FASTQs fall out of the page cache. There is no swapping.

GPU extension lanes (tiles in flight; fused, two-pass, align in s): 2: 314.4, 4: 250.1,
5: 239.2, 6: 237.8 / 240.0, 7: 245.1, 8: 260.7. Metal has no queue priorities, so with
2 lanes extension tiles queue behind seeding and fused extension costs more than it saves.
The Metal default is 6.

## 30× whole genome (HG002 SRR14724532, 350 M pairs), FASTQ → VCF, 16 threads

Whole-genome mode turns on by itself. The x86 golden's scatter order and `@RG DS` are used,
and the gate is `reference/wgs30-x86-r1.digest.tsv`.

| Configuration (dense SA interval 4) | total | prep | align | markdup | BQSR | HC | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| CPU arm | 10,421 s (2 h 54 m) | 487 | 8,903 | 283 | 212 | 1,010 | 90 GB |
| **Metal, fused extension** | **4,973 s (1 h 23 m)** | 502 | 3,465 | 284 | 214 | 997 | 105 GB |
| Metal, same, `NOMICS_TIMING=1` | 5,592 s | 444 | 4,072 | 284 | 215 | 1,009 | 105 GB |

**Metal is 2.10× the CPU arm end to end (align 2.57×).** `NOMICS_TIMING=1` itself costs
~15% of alignment at 30×, because its per-lookup counters contend across 16 threads, so
profiles taken with it overstate CPU time.

| Metal, fused, **dense SA interval 2** (15 GB; cache build ~1 min included) | 5,015 s | 502 | 3,444 | 284 | 216 | 997 | 108 GB |

**Dense SA interval 2 does not shorten the Metal run** (align 3,444 vs 3,465 s), although
it cuts chaining CPU (32% → 11.5% of samples). The Metal arm at 30× is therefore **not
CPU-bound**: `sample` shows ~57% of thread time outside wait states. The timing log's
"align busy" counts lanes waiting on GPU tiles as busy, so do not read it as CPU use.
The limit is GPU throughput or how the batch pipeline overlaps CPU and GPU work.

**The GPU is saturated:** `powermetrics` during a 30× Metal run shows GPU HW active
residency 94–100% at its top clock (1,578 MHz), while a P-core cluster is 76–99% active.
Summed kernel time agrees. Mate rescue is 799 M calls and 1,936 s, seeding 1,300 s, fused
extension ~1,000–1,400 s, against a ~4,000 s timeline.

Matched configuration, dense SA interval 2 (every run 41/41):

| HG002 30×, dense SA 2 | total | prep | align | markdup | BQSR | HC | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| CPU arm | 10,161 s (2 h 49 m) | 489 | 8,647 | 283 | 213 | 998 | 96 GB |
| Metal, fused, 6 GPU lanes | 5,015 s | 502 | 3,444 | 284 | 216 | 997 | 108 GB |
| **Metal, fused, 6 GPU + 2 CPU extension lanes** (`NOMICS_GPU_EXTEND_CPU_LANES=2`) | **4,707 s (1 h 18 m)** | 502 | **3,191** | 283 | 216 | 997 | 108 GB |
| Metal, fused, 6 GPU lanes, a quarter of mate rescue on the CPU (`NOMICS_GPU_RESCUE_SHARE=0.75`) | 4,707 s | 500 | 3,190 | 283 | 216 | 997 | 108 GB |
| **Metal, fused, 6 GPU + 2 CPU lanes and rescue share 0.75 (best)** | **4,592 s (1 h 17 m)** | 388 | **3,080** | 283 | 215 | 997 | 99 GB |

**At matched configuration, Metal is 2.21× the CPU arm end to end (align 2.81×)** with
both levers: `NOMICS_GPU_EXTEND=fused NOMICS_GPU_EXTEND_CPU_LANES=2
NOMICS_GPU_RESCUE_SHARE=0.75` (4,592 vs 10,161 s). With
the GPU saturated, CPU lanes that steal extension tiles add throughput at 30×, while on the
GPU-light exome they did not. `NOMICS_GPU_RESCUE_SHARE` (`16dc7c7`) moves a fraction of mate
rescue back to idle cores. A quarter of rescue on the CPU gives the same −7% align as the
CPU extension lanes: both relieve the saturated GPU, and together they give another −3.5%.

## Where the CPU goes at 30× (Metal arm, `sample`, interval 2)

Share of non-waiting CPU samples during alignment: `chain2aln` / region patching 18%;
BAM and chunk deflate 12%; chaining 11.5% (32% on tier1 at interval 4); allocator /
memmove 9%; reference fetch 8%; `global2` (final CIGAR) 6%; pairing / rescue / MAPQ 6%.
The deflate is mostly `write_chunk`'s zlib-exact chunk FASTQs, which exist only for the
chunk gates (`--write-intermediates`), so production runs skip it.

## Comparison points (other hardware and OS; not a speedup claim)

- Spark (GB10, 20 Arm cores, Linux): tier1 298 s and 30× 4,003 s with fused extension.
- sarek on the Spark: 30× 41,044 s.

## Recommended Metal invocation

`nomics run` built with `--features metal`, plus `NOMICS_GPU_SEED=1 NOMICS_GPU_EXTEND=fused`.
For whole genomes, also `NOMICS_GPU_EXTEND_CPU_LANES=2 NOMICS_GPU_RESCUE_SHARE=0.75`
(30×: 4,592 s, 41/41). Exomes keep the defaults: on the exome the GPU is not saturated, and
the in-memory run peaks at ~125 GB on a 128 GB machine. Defaults stay conservative; the
whole-genome knobs were measured on one machine and one sample.

## Next levers (not started)

- The GPU is the 30× bound. Mate rescue is its largest consumer (799 M calls), so a faster
  rescue kernel (e.g. two problems per simdgroup) or an adaptive rescue share would pay
  first.
- `global2` (final CIGARs, 6% of 30× CPU) on the GPU needs a two-phase output path
  (`reg2aln` → `gen_cigar2` runs per read at emission). Not worth it while the GPU is the
  bound.
- Exome memory: whole-genome (spill) mode would let the exome use dense SA 2 and drop the
  ~125 GB peak. Not measured.
