# Scope: an exact BWA 0.7.17 + GATK 4.0.12.0 port (UK Biobank deCODE profile)

Sources: `port-gatk-4.0.12-diff.md` (GATK, Picard, htsjdk and GKL internals) and
`port-ukb-pipeline-and-goldens.md` (pipeline steps, BWA, goldens), both 2026-10-17, plus
an internal survey of which versions are validated in production (why this version).

## What would be ported

The UK Biobank 500k WGS deCODE profile (247,522 samples). The single-sample gVCF is the
GATK boundary; cohort calling is GraphTyper, not GATK.

1. `bwa mem` 0.7.17, `-K 100000000 -Y`, alt-aware against GRCh38 + alts + decoys + HLA
   (no bwa-postalt).
2. samblaster 0.1.24 `-a --addMateTags`, which adds mate tags only (no duplicate marking).
3. samtools 1.9 sort.
4. Picard MarkDuplicates 2.20.3 in queryname mode, per read group.
5. GATK 4.0.12.0 BaseRecalibrator and ApplyBQSR, with FE 4-bin quality binning
   (SQQ 10/20/30) on the BAM fed to HaplotypeCaller.
6. GATK 4.0.12.0 `HaplotypeCaller -ERC GVCF`, scattered by chromosome.

The Sanger profile is blocked: its scatter interval list is not public, and output depends on
the scatter boundaries.

## Two findings that change the decision

1. **4.0.12.0's HaplotypeCaller is not deterministic in places.** `SeqVertex` hashes by
   `System.identityHashCode`, and graph code iterates HashSets of vertices; GATK fixed this in
   4.1.4.0 (#6195). Where that order matters, no port can be identical by construction. The
   reference cannot even reproduce itself. Before any port work, run 4.0.12.0 twice on the same
   input and measure how often it disagrees with itself. If the residual is not tiny, "exact
   4.0.12.0" is not a promise anyone can keep, and the target should move to the first
   deterministic version labs use (4.1.4.0 or later).
2. **Production UKB gVCFs were almost certainly computed in GKL mode on x86** (AVX PairHMM in
   float with a double fallback, FTZ, M and X columns summed separately). A Java-mode port
   matches a run nobody shipped. Matching what labs actually hold means porting GKL 0.8.6's
   PairHMM too. That looks feasible by construction, but it is a second numeric model to
   maintain, and GKL 0.8.6 and 0.8.11 already differ from each other (a hard-coded LN10).

Most of the HaplotypeCaller we ported was rewritten after 4.0.12.0, so the delta is a large
partial re-port rather than a patch:
- haplotype path finding;
- likelihood capping;
- which reads feed PL, AD and annotations;
- window padding;
- left-align and realign;
- overlapping-pair handling;
- the old exact AF calculator, which also drives isActive (about 1,500 lines deleted by 4.6.2.0);
- ExcessHet's old mid-p value;
- rank sums over deletions;
- NON_REF GT and PL rules;
- `qualToErrorProbLog10` arithmetic.

These carry over:
- the gVCF annotation set, GQ bands, PCR indel model, downsampling and its RNG seed;
- kmer sizes, read filters, Smith-Waterman parameters and the Java PairHMM;
- the BQSR covariates, quantizer and report format;
- the SQQ binning in ApplyBQSR;
- MarkDuplicates' duplicate-flagging logic.

## Effort

Engineer-days are for a person working without our agent workflow; see the calibration note
below.

| Component | Low | High | Notes |
|---|---:|---:|---|
| HaplotypeCaller deltas (Java mode), incl. the old exact AF calculator | 45 | 95 | Most of the cost; some spots may be unmatchable (finding 1) |
| BQSR deltas: empirical-quality posterior, known-site masking | 4 | 8 | |
| GenotypeGVCFs / GVCF writer deltas | 4 | 8 | Only if per-sample GenotypeGVCFs is in scope |
| htsjdk 2.18 output formatting (VCF encoder, headers) | 3 | 8 | |
| Picard MarkDuplicates, queryname input (2.18/2.20 semantics) | 6 | 17 | Our port takes coordinate-sorted input only |
| BWA 0.7.17 mode (MQ tag, `@HD`, `@PG`) | 1 | 1 | Engine unchanged |
| Alt-aware alignment (`.alt` loading, `AH:*` headers, gating) | 4 | 10 | Logic mostly present in nomics-align, disabled |
| samblaster mate tags (possible duplicate MC tags), samtools 1.9 sort and merge order, raw FASTQ entry, CRAM gated on records | 8 | 20 | |
| **Java-mode subtotal** | **75** | **167** | |
| GKL 0.8.6 PairHMM (float/double fallback, FTZ, column order) | 16 | 33 | Needed to match x86 production output |
| Goldens: x86 drivers, r1/r2 pairs on wgs30s then 30×, a GKL re-run, arm64 cross-check | 5.5 | 12 | Plus about $65–100 of cloud compute |
| **Total (Java mode + GKL + goldens)** | **~97** | **~212** | Roughly 5–10 engineer-months |

Calibration: the existing exact sarek germline port (fastp through HaplotypeCaller, GVCF and
joint calling) took about nine calendar days with this repo's agent workflow (worker pool,
oracle dumps, gates). These tables are written in conventional engineer-days, so treat them
as an upper bound on calendar time with the same workflow. The relative sizes are the
reliable part: HaplotypeCaller dominates, and GKL is a meaningful add-on.

## Recommendation

Do not start the port yet. In order:

1. **Measure 4.0.12.0's own non-determinism first.** About 2 days and about $30: r1 vs r2 on
   wgs30s in Java mode and in GKL mode. If 4.0.12.0 disagrees with itself in more than a
   handful of records, retarget to 4.1.4.0+ before spending anything else.
2. **Ask labs which version they validated** (the cheapest way to raise
   confidence). A named buyer version beats UK Biobank as a proxy, because UKB is a research
   pin.
3. **If the answer is still 4.0.x or 4.1.x, port HaplotypeCaller's deltas first,** gated on
   the gVCF from a golden recal.bam. HC is the long pole and the riskiest piece. BQSR,
   MarkDuplicates and BWA are comparatively cheap and can follow.
