# Does GATK 4.0.12.0 reproduce itself? (2026-10-17/18)

*The cloud and kit scripts this document names (`scripts/cloud/`, `scripts/kit/`) drive one owner's AWS account and are not part of the public repository; the measurements and code changes they produced are.*

Question from the 4.0.12.0 port scope (`docs/strategy/port-4.0.12-scope.md`): GATK
4.0.12.0's assembly graph iterates HashSets keyed by `System.identityHashCode` (fixed in
4.1.4.0, #6195). If that made HaplotypeCaller non-deterministic, no port could be exact.

## Setup

- `scripts/cloud/gatk4012_determinism.sh` on a Verda CPU.32V.128G (AMD EPYC 7C13, AVX2,
  no AVX-512).
- Image: `broadinstitute/gatk:4.0.12.0`. FASTA, `.fai` and `.dict`: Broad hg38 bundle.
- Input: the wgs30s slice's recal.bam (HG002, 1/12 of the 30× reads, whole-genome).
- `HaplotypeCaller -ERC GVCF` with defaults, scattered by chromosome as deCODE did for UK
  Biobank, over the whole-genome calling regions.
- Each chromosome is one JVM, with up to 28 at once in Java mode and 8 in GKL mode.
- **Java mode:** `--pair-hmm-implementation LOGLESS_CACHING --smith-waterman JAVA` plus the JDK
  deflater and inflater. The log reports "the non-hardware-accelerated Java LOGLESS_CACHING
  PairHMM implementation".
- **GKL mode:** GATK's defaults. The log reports "the OpenMP multi-threaded AVX-accelerated
  native PairHMM implementation", FTZ enabled, 4 threads, IntelDeflater.

Each mode ran twice (r1, r2). Records were compared with headers excluded, per chromosome
(`compare.tsv`; run logs in `logs.tgz`, excerpt in `run-summary.txt`).

## Result

| Comparison | gVCF records per run | Differing lines |
|---|---:|---:|
| Java r1 vs Java r2 | 140,277,074 | **0** |
| GKL r1 vs GKL r2 | 140,277,074 | **0** |
| Java r1 vs GKL r1 | 140,277,074 | **1,786** |

Differing lines count both sides, so about 900 records (~0.0006%). They are spread over every
chromosome (4–362 per chromosome), which looks like float-vs-double PairHMM rounding rather
than a single localized desync.

## What it means

- **GATK 4.0.12.0 reproduced itself exactly here**, in both modes. On this workload and JVM,
  the identity-hash iteration did not change any output. That removes the scope's first
  blocker: an exact 4.0.12.0 port is a coherent target. Caveat: two runs per mode on one
  machine type and JVM. Identity hashes can differ with thread scheduling, so this is
  evidence, not proof. Every future 4.0.12.0 golden gets an r1/r2 pair anyway.
- **Java mode and GKL mode disagree on ~900 records.** A lab running GATK 4.0.12.0 on x86 with
  defaults holds GKL-mode output. An exact port of what labs validated therefore needs GKL's
  PairHMM (float first with double fallback, FTZ, M and X columns summed separately), as the
  scope assumed. This box is AVX2; the AVX-512 path may round differently again (unverified).
- PairHMM compute time on chr1: GKL 20 s vs Java 361 s.

Cost: $7.51. The run itself was ~$2; the rest was idle time after it finished, while the
teardown waited on a `--yes` confirmation.
