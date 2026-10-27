# MarkDuplicates stage

Inputs: the reference's 12 per-chunk BAMs (`preprocessing/mapped/*.NNNN.bam`). Gate: `nomics dev gate-bam` against golden `md.bam`, plus a metrics diff that skips lines 2 and 4 (command line, start time).

| Run | Records | Gate | Metrics | nomics wall | GATK wall (trace) |
|---|---|---|---|---|---|
| tier0 (2026-10-08) | 977,738 | **IDENTICAL** | identical | 2.5 s, 0.6 GB RSS | 8.8 s, 10.8 GB |
| tier1 (2026-10-08) | 81,578,107 | **IDENTICAL** | identical | 232 s, 301 s user, 32.6 GB RSS | 358 s, 1.5 cores, 21.9 GB |

Evidence: `~/omics/runs/nomics-markdup-tier{0,1}-20261008-01/` (`run.log`, outputs).

## Parallelization (2026-10-08)

Every row was gated: all records IDENTICAL to the golden `md.bam`, metrics identical.

| tier1 version | decode+merge | pass 1 | select | pass 2 | wall | user |
|---|---|---|---|---|---|---|
| 01 sequential | 19.6 s | 36.5 s | 7.4 s | 174.5 s (encode 40.0, BGZF 132.2) | 232 s | 301 s |
| 03 + parallel encode and `ParallelBgzfWriter` | 19.3 s | 36.4 s | 7.4 s | 21.2 s | 85.5 s | 378 s |
| 04 + parallel pass 1 via `RecordView` | 19.1 s | 7.9 s | 7.3 s | 21.2 s | **56.7 s** | 357 s |

That is 6.3× faster than GATK's 358 s. tier0: 0.68 s, against GATK's 8.8 s.

`ParallelBgzfWriter` compresses full blocks with fresh deflater state and the last block with a primed deflater (see the argument in `nomics-io/src/bgzf.rs`). Rewriting the tier1 `md.bam` and `recal.bam` through it is byte-identical. Remaining serial work: the k-way merge, the pair map, sorting in select, and the in-order write loop. Peak RSS is 32.6 GB, because every record is resident.

Original note: the first implementation is effectively single-threaded: only input decoding runs per-file in parallel; merge, re-encode and BGZF compression run on one thread. Next: per-stage timing, then parallel re-encode and compression (fresh-state deflate per block, measured safe on 668k golden blocks).
