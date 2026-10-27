# BQSR stage (BaseRecalibrator + ApplyBQSR)

Input: the golden `md.bam`, the GATK GRCh38 reference, three known-sites VCFs and the Agilent v7 BED (one interval group, as sarek ran it). Gates: `nomics dev gate-report` on the recal table (bytes) and `nomics dev gate-bam` on `recal.bam` (records, header excluding `@PG`).

| Run | Table | recal.bam | nomics wall | GATK wall (trace) |
|---|---|---|---|---|
| tier0, no cache (2026-10-08, run 01) | **IDENTICAL (bytes)** | **IDENTICAL** (873,049) | 38.6 s (36.8 s parsing dbSNP) | 15.6 s + 9.5 s |
| tier0, known-sites cache (run 03) | **IDENTICAL (bytes)** | **IDENTICAL** | **2.38 s** | 25.1 s |
| tier1, no cache (run 01) | **IDENTICAL (bytes)** | **IDENTICAL** (40,074,075) | 74.1 s | 518 s + 355 s = 873 s |
| tier1, known-sites cache (run 02) | **IDENTICAL (bytes)** | **IDENTICAL** | **37.8 s** (336 s user, 36.8 GB RSS) | 873 s (**23× faster**) |

tier1 stages with the cache: load+decode BAM 16.4 s, known sites 0.6 s, reference 1.0 s, recalibrator pass 4.2 s (parallel, per-thread tables merged exactly), apply pass 14.3 s (QUAL patched in place, `ParallelBgzfWriter`).

The known-sites cache (`--known-sites-cache`, 1.18 GB for these three VCFs) holds deduplicated per-contig site extents. It is keyed by source path, size and mtime and by the BAM contig names, and is rebuilt automatically when any of them changes.

Spec: `docs/spec/bqsr.md`. Before implementing, I checked its load-bearing claims on the goldens: the `EstimatedQReported` recurrence; empirical quality for all RecalTable0/1/2 rows at both tiers; the table layout (100-byte Arguments lines, padded headers, `\n\n` ending); and that `recal.bam` is an ordered subset of `md.bam` with only QUAL changed and Q<6 untouched. Evidence: `~/omics/runs/nomics-bqsr-tier{0,1}-20261008-0*/`.

## Chained: `nomics preprocess` (markdup → BQSR in one process)

MarkDuplicates hands its output records to BQSR in memory (`RecordSet`). md.bam is still written, because it is a published output, but on a thread that overlaps BQSR. Nothing is read back.

| Run | md.bam | metrics | recal table | recal.bam | wall | peak RSS |
|---|---|---|---|---|---|---|
| tier0 (2026-10-08) | IDENTICAL (977,738) | identical | IDENTICAL (bytes) | IDENTICAL (873,049) | 3.17 s | 8.0 GB |
| tier1 (2026-10-08) | IDENTICAL (81,578,107) | identical | IDENTICAL (bytes) | IDENTICAL (40,074,075) | **77.6 s** | 56.6 GB |

tier1 comparison: GATK MarkDuplicates + BaseRecalibrator + ApplyBQSR = 358 + 518 + 355 = 1,231 s (15.9× faster). Separate nomics commands = 56.7 + 37.8 = 94.5 s.
Costs: (1) peak RSS rises because the merged input arena and the output stream coexist during pass 2 (fixable by streaming); (2) the overlapped md.bam write (26.9 s) competes for cores with the apply pass (20.4 s vs 14.3 s alone).
