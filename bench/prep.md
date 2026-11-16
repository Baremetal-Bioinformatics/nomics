# Read prep stage (fastp 1.1.0 equivalent)

Gate: `nomics dev gate-chunks` compares the 24 chunk FASTQs (`NNNN.HG002_WES-L1_R{1,2}.fastp.fastq.gz`) from the reference run's fastp work dir by decompressed content. The filter counts are also compared with fastp's JSON `filtering_result`.

| Run | Chunks | Filter counts vs fastp JSON | nomics wall | fastp wall (trace) |
|---|---|---|---|---|
| tier0 (2026-10-08) | **24/24 IDENTICAL** | identical (975,636 passed, 15,172 low quality, 310 N) | 0.78 s | — |
| tier1 (2026-10-08) | **24/24 IDENTICAL** (~3.37M pairs each) | identical (80,987,784 / 750,384 / 16,388) | 51.2 s (filter 29.1 s + chunk write 22 s) | 44.9 s (12 threads) |

Spec, read from fastp v1.1.0 source (`peprocessor.cpp`, `filter.cpp`, `evaluator.cpp`, `common.h`):
- Packs of 256 pairs go to worker `pack % 12`, and each worker writes its own split files.
- A read fails on quality first (> 40% of bases below Q15), then N bases (> 5), then length (< 15). A pair needs both reads to pass.
- Failed pairs are counted as 2 reads under the worse fastp code.
- Poly-G auto-trimming is triggered by two-colour read-name prefixes and is refused rather than implemented.

Not faster yet: the filter pass is bound by single-threaded inflate of ENA's single-member gzip, two streams in parallel. Writing chunk files disappears once prep feeds BWA in memory. Faster inflate is a later item.
