# Results

Every number on this page comes from a run whose outputs were gated identical to sarek
3.10.0 at every stage boundary (41 artifacts for a 30× genome, 18 for an exome) unless it
says otherwise. Hardware: one NVIDIA DGX Spark (GB10 GPU, 10 Cortex-X925 + 10 Cortex-A725
cores, 128 GB unified memory). Input: HG002, GRCh38 (iGenomes GATK bundle).

## Headline

| FASTQ → VCF | sarek 3.10.0 | nomics | nomics vs sarek |
|---|---:|---:|---:|
| 30× whole genome (350M pairs) | 11 h 24 min (41,044 s) | **1 h 00 min (3,619 s)** | **11.3×** |
| Whole exome, Agilent v7 (40.9M pairs) | 1 h 08 min (4,107 / 4,134 s) | **4 min 58 s (306 / 290 s)** | **13.8×** |

- **Identical output.** nomics: 41/41 artifacts at 30× (703M mapped records, 4,773,797 VCF
  records) and 18/18 on the exome.
- **Which sarek.** The goldens are sarek on arm64 and, for the 30× genome, sarek on x86 with
  HaplotypeCaller in Java mode, which matches arm64 record for record. sarek's default x86
  configuration (Intel GKL) differs in about 2% of 30× records
  ([docs/REFERENCE-VARIABILITY.md](docs/REFERENCE-VARIABILITY.md)).
- **Same box, one job at a time,** page caches evicted between runs. The sarek times and the
  exome nomics times are from the 2026-10-17/18 and 2026-10-23 bake-off; the 30× nomics time
  is the 2026-10-27 run, the first after a reboot (three more runs that session, on
  fragmented memory, took 3,930–4,024 s).
- **Peak memory** at 30×: nomics 72.5 GB in the 2026-10-23 bake-off (77 GB in the
  2026-10-27 build).

### Where the 30× hour goes

| Stage | nomics | runs on |
|---|---:|---|
| Prep (fastp filtering and chunking) | 323 s | CPU, overlapped with alignment |
| Alignment (BWA-MEM + sort) | 2,172 s | GPU seeding, extension and mate rescue; CPU for the rest |
| MarkDuplicates | 286 s | CPU |
| BQSR (BaseRecalibrator + ApplyBQSR) | 201 s | CPU |
| HaplotypeCaller | 938 s | CPU, with a share of PairHMM on the GPU |
| **Total** | **3,619 s** | |

Stage times overlap slightly, so they do not sum exactly to the total. Alignment is now 60%
of the run. Its main remaining limit is the FM-index lookups BWA's algorithm requires, which
nomics must reproduce one for one (see [docs/LEARNINGS.md](docs/LEARNINGS.md)).

## How it got here

The 30× genome, FASTQ → VCF, as the work progressed. Every row is 41/41 identical.

| Date | Wall time | What changed |
|---|---:|---|
| 2026-10-11 | 9,994 s (2 h 47 min) | First whole-genome run, every stage streamed through disk |
| 2026-10-13/14 | 7,135 s | GPU mate rescue; prep overlapped with alignment |
| 2026-10-16 | 6,467 s | Parallel GVCF, NEON Smith-Waterman, mimalloc |
| 2026-10-20 | 5,240 s | Fused GPU seed extension (kernel v3, priority stream), 2.6× faster seed kernel |
| 2026-10-22 | 4,572 s | Alignment pipeline restructured: persistent engines, asynchronous GPU work, faster region merge, compact seeds |
| 2026-10-23 | 4,003 s | Three batches in flight, faster BAM compression; first run after a reboot |
| 2026-10-25 | 3,818 s | BQSR tables built during duplicate marking; fused extension the default |
| 2026-10-27 | **3,619 s** | Mate-rescue kernel 2.2× faster, seed-kernel trims, two-pass seeding |

The exome went from 785 s (2026-10-09, separate commands per stage) to 704 s as one process
the same day, then 459/462 s (10-17/18 bake-off), 379 s (10-19) and 298 s (10-23, mean of
two runs).

## Other platforms

| Platform | 30× FASTQ → VCF | Exome |
|---|---:|---:|
| DGX Spark, CUDA | 1 h 00 min | 4 min 58 s |
| Apple M4 Max (16 cores, 40-core GPU), Metal | 1 h 17 min (4,592 s) | 5 min 20 s (320 s) |
| Apple M4 Max, CPU only | 2 h 49 min (10,161 s) | 11 min 33 s (693 s) |

Every run on every platform is identical to the same sarek goldens.

Metal is experimental: seed search, mate rescue and seed extension run on the GPU;
PairHMM stays on the CPU, since Metal has no double precision
([bench/metal.md](bench/metal.md)).

## Caveats

- **Memory state moves times by 5–10%.** On a machine that has been up for a while, the
  index and GPU buffers may not get 2 MiB pages; the same 30× build then takes 3,930–4,024 s,
  still with identical output. The best numbers above are from runs right after a boot.
  [REPRODUCE.md](REPRODUCE.md) explains how runs flag this.
- **Single runs at 30×.** The 30× genome takes an hour per run, so most 30× rows are one run;
  chunk-level and exome measurements are repeated. Differences under ~30 s are not claimed.
- **One sample, one reference.** All gates are on HG002 with GRCh38. Other samples go
  through the same code, but have not been gated against their own sarek goldens.
- **What is not compared:** fastp's QC reports and `samtools merge`'s side output
  (`mapped/*.sorted.bam`), whose read-group suffixes are random in the reference itself
  ([docs/GATES.md](docs/GATES.md)).

Per-stage detail: [bench/](bench/). Design notes behind each step of the timeline:
[docs/strategy/](docs/strategy/).
