# Results

Every number on this page comes from a run whose outputs were gated identical to sarek
3.10.0 at every stage boundary (41 artifacts for a 30× genome, 18 for an exome) unless it
says otherwise. Hardware: one NVIDIA DGX Spark (GB10 GPU, 10 Cortex-X925 + 10 Cortex-A725
cores, 128 GB unified memory). Timings are on HG002, GRCh38 (iGenomes GATK bundle);
[nine more datasets](#exactness-across-samples) pass the same gates.

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

## Exactness across samples

The timings above are HG002. The same gates were then run on nine more inputs
(2026-11-05 to 11-09), each against its own pair of sarek 3.10.0 runs, r1 and r2, made on the
same Spark. Every artifact matched on every sample.

| | Sample | Data | Gates | VCF records |
|---|---|---|---|---:|
| E1 | HG001 | exome, Agilent v7, NovaSeq 6000 | 18/18 | 42,189 |
| E2 | HG003 | exome, IDT xGen, NovaSeq 6000 | 18/18 | 22,436 |
| E3 | HG004 | exome, IDT xGen, HiSeq 4000 (2×150) | 18/18 | 22,151 |
| E4 | HG005 | exome, Agilent v7, HiSeq 4000 (2×150) | 18/18 | 39,285 |
| E5 | HCC1395BL | exome, Agilent v6 + UTR, HiSeq 1500 (2×126), 88.5M pairs | 18/18 | 137,147 |
| E6 | HG006 | exome, IDT xGen, NovaSeq 6000 | 18/18 | 22,388 |
| G1 | HG002 | genome, 10×, NovaSeq X | 41/41 against r1, 41/41 against r2 | 4,283,337 |
| G2 | HG001 | genome, 20×, PCR-free, NovaSeq 6000 | 41/41 against r1, 41/41 against r2 | 4,680,788 |
| G3 | HG003 | genome, 20×, PCR library, NovaSeq 6000 | 41/41 against r1, 41/41 against r2 | 4,635,128 |

- **What the rows add to HG002:** six more individuals, two more exome kits, three more
  instruments, two more read lengths (2×150 and 2×126), a PCR genome library, and 10× and 20×
  depths.
- **Instrument read names.** E1–E4, E6, G2 and G3 keep the sequencer's read names, so
  MarkDuplicates' optical-duplicate detection is on, and on the NovaSeq 6000 samples so is
  fastp's poly-G trimming. Reads downloaded from SRA are renamed and switch both off; HG002
  had only exercised them on a small two-lane fixture.
- **sarek against itself.** For the exomes, r1 and r2 agree at every gated boundary. For the
  genomes they differ in one artifact, the recalibrated BAM, whose read-group IDs carry
  `samtools merge`'s clock-seeded suffixes. nomics matched each run from that run's own
  recorded scatter order ([docs/REFERENCE-VARIABILITY.md](docs/REFERENCE-VARIABILITY.md)).
- **One bug found.** HG006 was the first input with optical-duplicate coordinates above
  32,767. Picard keeps them as 16-bit values, so they wrap and change which duplicate is kept;
  nomics had not reproduced that. It does now ([docs/LEARNINGS.md](docs/LEARNINGS.md)). The
  HG002 results were unaffected: its SRA read names keep optical detection off.

The data plan, sources and md5s are in
[docs/strategy/sample-matrix.md](docs/strategy/sample-matrix.md).

**Synthetic edge cases.** Fourteen generated read sets aim at inputs real samples rarely
contain: homopolymer indels, long indels, fragments shorter than the reads, N runs,
duplicates around Picard's pixel distance, chimeras, unmapped mates, contig edges, quality
extremes, read lengths from 1 to 151, pileups of 3,000–8,000×, clustered variants, fully
overlapping mates, and coordinates past 16 bits. Each went through sarek twice and passes
18/18, with 435 to 690 VCF records apiece (2026-11-11). The reads are a function of a table
row, so any failure is a ready reproducer
([docs/strategy/fuzz-goldens.md](docs/strategy/fuzz-goldens.md)).

## Repeatability

Twenty full HG002 30× runs were queued back to back on one build (2026-11-09/10), each gated
41/41 by digest, with internal region checks on (`NOMICS_DIAG=regions`). The first stopped
when the disk filled. The other nineteen were clean: every artifact identical, no GPU
fallback.

| Runs | Configuration | Pipeline wall time |
|---:|---|---:|
| 13 | default | 3,869–3,948 s (64.5–65.8 min) |
| 1 | default, right after 739 GB of old runs were moved to network storage | 5,639 s |
| 5 | seed extension on the CPU (`NOMICS_GPU_EXTEND=0`) | 4,909–4,959 s |

The Spark had been up for eleven days when the runs began. So about 65 minutes is what a
machine in ordinary use gives; the 3,619 s headline is the best run, the first after a reboot.

Peak disk at 30× is 288 GB in the run directory, with the gate's intermediates kept
(2026-11-11). Each chunk's uncompressed spill is deleted once that chunk has aligned; before
that change the same run would have peaked near 500 GB.

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
  index and GPU buffers may not get 2 MiB pages; a 30× run then takes 3,869–4,024 s, still
  with identical output ([Repeatability](#repeatability)). The best numbers above are from
  runs right after a boot. [REPRODUCE.md](REPRODUCE.md) explains how runs flag this.
- **Single runs in the timeline.** Each row of "How it got here" is one 30× run, an hour
  each; chunk-level and exome measurements are repeated. Differences under ~30 s are not
  claimed.
- **Seven individuals, one reference, one pipeline version.** The gates cover HG002 (exome
  and 30× genome) and the nine datasets above, all on GRCh38 against sarek 3.10.0. All are
  Illumina short reads, and the only genome above 20× is HG002.
- **What is not compared:** `samtools merge`'s side output (`mapped/*.sorted.bam`), whose
  read-group suffixes are random in the reference itself, and FastQC's rendered pages and
  images ([docs/GATES.md](docs/GATES.md)). The timings above are of runs without the QC
  reports; `nomics run --reports` writes them and is not yet tuned (an exome takes 373 s
  with them against 278 s without; at 30× the reports made at the end take 370 s).

Per-stage detail: [bench/](bench/). Design notes behind each step of the timeline:
[docs/strategy/](docs/strategy/).
