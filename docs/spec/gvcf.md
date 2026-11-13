# HaplotypeCaller GVCF (GATK 4.6.2.0)

Scope: the single-sample diploid sarek 3.10.0 `--joint_germline` HC invocation,
`-ERC GVCF --pcr-indel-model CONSERVATIVE`, default annotations and GQ bands.
This emits the input to joint genotyping (`docs/spec/joint.md`).

## Source contract

Pinned source: [GATK 4.6.2.0](https://github.com/broadinstitute/gatk/tree/4.6.2.0).
Paths below are under `src/main/java/org/broadinstitute/hellbender/`.

- `tools/walkers/haplotypecaller/HaplotypeCallerEngine.java`,
  `validateAndInitializeArgs`, `callRegion`, `containsCalls`,
  `referenceModelForNoVariation`: confidence becomes **negative zero**, PLs are
  forced on, physical phasing stays enabled. Inactive/no-variation regions use
  finalized filtered reads. With calls, the left/right flanks use those reads
  clipped to each flank's padded span; the variant span uses realigned,
  likelihood-filtered evidence. If no genotype is called, use the original
  region's reference model. Reference-confidence records are generated inside
  deletions, then suppressed by the block writer.
- `tools/walkers/haplotypecaller/ReferenceConfidenceModel.java`: SNP likelihoods
  use Q > 6 (deletions use Q30), cap alternate GLs at hom-ref, compare unrounded
  hom-ref GQ with the indel model (Q45, at most 40 informative reads), then
  integerize the chosen PL vector. Indel informativeness tests offsets 1–10,
  caches from the first eligible pileup, and preserves the source's original
  read-length versus CIGAR-expanded-length boundary conditions. Each reference
  model invocation has a fresh cache.
- `utils/genotyper/AlleleLikelihoods.java`, `updateNonRefAlleleLikelihoods`:
  `<NON_REF>` is the median of concrete likelihoods **strictly below** the best;
  all ties use the best (NaN with at most one concrete allele). Append it before
  genotype likelihoods, AF calculation, allele subsetting and annotation.
- `tools/walkers/genotyper/GenotypingEngine.java`: `forceKeepAllele` in the HC
  subclass retains alleles in GVCF, but the monomorphic emit threshold still
  applies. `utils/variant/AlleleSubsettingUtils.java` keeps NON_REF **in addition
  to** six concrete alternates. htsjdk 4.2.0 `SimpleAllele.length()` returns zero
  for symbolic alleles, including in AF priors and upstream-deletion bookkeeping.
- `utils/variant/writers/{GVCFBlockCombiner,HomRefBlock,GVCFBlock}.java`:
  boundaries 1 through 60, 70, 80, 90, 99, 100; membership uses GQ capped at 99;
  componentwise minimum PLs, GQ from that vector, minimum DP, and rounded median
  DP. Blocks break at gaps, contigs, variants and GQ bands. DP is stored as an
  exact histogram so long uncovered intervals do not require a per-base buffer.
  Each scatter group has a fresh combiner; gather never recombines blocks.
- `tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java`, `phaseCalls`:
  use the haplotypes mapped to emitted alleles; connect equal/subset-on-all or
  disjoint complementary haplotype sets in call order, and abandon inconsistent
  phase groups. PGT describes alternate haplotype phase even for hom-var calls.
  PID is `POS_REF_firstALT`; PS is the first call position. GT is phased and
  reversed where required. Set iteration affects only membership/cardinality,
  not pair order or phase-group assignment.

## Annotations and random draw order

`HaplotypeCallerEngine.filterReferenceConfidenceAnnotations` (lines 264–278)
**removes QualByDepth, ChromosomeCounts, FisherStrand and StrandOddsRatio**, and
adds StrandBiasBySample. Thus the default GVCF command emits no QD, AC/AF/AN,
FS or SOR. RMSMappingQuality remains raw (`RAW_MQandDP`); the normal-mode tail
calls `finalizeRawMQ`, the GVCF tail does not. MLEAC/MLEAF include NON_REF, with
AD, DP, GQ, PL and SB in FORMAT. Default annotations do **not** include the
AS_ group; adding AS_ fields would fail this fixture's header and records.

This establishes the shared RNG order without inferring it from output:

1. Each HC invocation starts `Utils`' shared `Random(47382911)`.
2. The positional reservoir downsampler consumes `nextInt(t)` for each read
   after the first 50 at an alignment start, in input/traversal order. The
   existing `reads::downsample` and `activity::traverse` reproduce its reservoir
   slots and traversal lookahead (see `haplotypecaller.md` §1.5). The pinned
   `utils/downsampling/ReservoirDownsampler.java:102–130` increments
   `totalReadsSeen` before drawing; `PositionalDownsampler.finalizeReservoir`
   resets the reservoir at each coordinate without reseeding the generator.
3. Normal-mode QualByDepth can consume `nextGaussian` when QD >= 35.
   In GVCF this annotation is removed **before** calls are annotated, so there
   are **zero Gaussian draws**, regardless of the larger set of genotyped events.
   NON_REF likelihoods, AF calculation, reference confidence, physical phasing
   and block compression do not draw from the shared RNG.
4. Consequently GVCF's RNG is exactly the uninterrupted sequence of reservoir
   draws, reset only at scatter-group boundaries. The GVCF annotation test
   checks that a call which consumes QD jitter in normal mode leaves the GVCF
   RNG state unchanged (subsequent integer draws compared against a fresh RNG).

The tier0 fixture has **zero downsampled alignment starts**. It gates the
changed event set and annotation output, but does not empirically exercise
reservoir/QD interleaving. The RNG claim above is source-derived; the existing
Java RNG oracle tests still apply. No normal-mode random state is reused.

## Reproduce the gate

The read-only fixture is `$HOME/omics-fixtures/gvcf-tier0`. `REF-READY` must
exist before end-to-end execution. Its recorded HC command has one interval
file, `chr1_12146-12310.bed`, spanning all Agilent v7 targets; passing the full
`Agilent_v7.bed` as one group reproduces all output loci and the single RNG reset.

```sh
scripts/gvcf_tier0_gate.sh "$HOME/omics-fixtures/gvcf-tier0" /tmp/tier0.g.vcf.gz
```

Or use `nomics hc --emit-gvcf --input X.recal.bam --intervals group.bed
--reference ref.fasta --dbsnp dbsnp.vcf.gz --output X.haplotypecaller.g.vcf.gz`.
Repeat `--intervals` for already-scattered groups. `nomics run --emit-gvcf`
uses the pipeline's existing scatter/gather and chooses the `.g.vcf.gz` suffix.
`nomics dev gate-gvcf` is an alias for the bytewise VCF gate, not a relaxed gate.
No new normalization is introduced. Output is BGZF through nomics-io.

## Tier0 measured result (2026-10-15)

- Normalized header identical; **251,885/251,885 records identical in order**.
- 250,733 hom-ref blocks; 1,152 variant records, including 105 phased records.
- No residual record differences. Initial multiallelic residuals at
  chr20:46128304 and chr20:63279199 were resolved by preserving NON_REF outside
  the six-concrete-alt limit; a regression test covers this rule.
- First exact example run: 10.2 s calling on this Mac with 8 CPU threads,
  excluding input loading and output compression. No performance claim against
  GATK or the saturated normal-mode CPU path is made.
- Tier1 **GVCF** is separately gated below. WGS GVCF is not gated by these
  fixtures; existing normal-mode WES/WGS gates are distinct. GVCF output
  indexing remains outside this change.

## Validation and remaining coverage

- `scripts/gvcf_tier0_gate.sh`: public CLI gate passes with zero mismatches.
- Normal mode on the same BAM/BED: 1,047/1,047 records and the entire compressed
  VCF are byte-identical to a clean build of pre-change HEAD
  `8dc0be903eae9aacd2c21f8f758c04f3f455fd2e`.
- `cargo test -p nomics-hc -p nomics-gate -p nomics`: the existing
  `tier0_gatk_default_annotations` test fails because its external
  fixture (`$NOMICS_FIXTURES/hc_finalize/t0_calls.tsv.gz`) is absent.
  Re-running with `-- --skip tier0_gatk_default_annotations` passes: 40 tests,
  nine already-ignored external/heavy oracle tests, one explicitly filtered
  missing-fixture test. No test was changed to hide this absence.
- `cargo clippy -p nomics-hc -p nomics-gate -p nomics --all-targets` passes with
  existing workspace warnings. The pre-existing non-CUDA empty GPU read-length
  range is now selected with `cfg`, avoiding the new toolchain's denied lint;
  CUDA eligibility behavior is unchanged.
- Formatting checks pass for touched Rust files; `git diff --check` passes.

## Tier1 follow-up: source diagnosis

The full exome reproduction initially had 1,713,849 records against
1,713,854 golden records. Aligning unchanged records isolated 112 differing
nomics lines and 117 golden lines: 91 versus 96 reference blocks, 19 no-call
records on each side, and two called variants on each side. The positional gate correctly reported
many more mismatches after the first block-count shift; no gate normalization
was changed to accommodate it. The deletion at chr1:84574341 is present and
identical in both supplied files; the neighboring reference blocks differed.

Three omissions in the GVCF paths account for these differences:

1. **Adaptor filtering after realignment.**
   [AssemblyBasedCallerUtils.getPileupsOverReference, lines 579–600][tier1-pileups]
   constructs a new `LocusIteratorByState`, which applies
   [ReadUtils.isBaseInsideAdaptor, lines 1066–1072][tier1-adaptor] at each locus.
   Finalization's earlier hard clipping does not replace this filter: local
   realignment changes read coordinates and can put bases beyond the fragment
   boundary. The cutoff is deliberately signed (`fragmentLength > 100`), so
   negative TLENs still undergo filtering. The reference-confidence pileup now
   excludes these elements from both SNP likelihoods/depth and cached indel
   informativeness. This explains the reference-block DP/PL and boundary
   differences, including chr1:84574325, chr7:29566635, chr9:109119660,
   chr13:45596538 and chr16:4785779. The indel-size limit, call selection and
   block-combining rules needed no changes.
2. **No-call and empty-evidence annotation guards.**
   [DepthPerAlleleBySample][tier1-ad], [DepthPerSampleHC][tier1-dp] and
   [StrandBiasBySample][tier1-sb] return without annotating an uncalled genotype;
   [ExcessHet][tier1-eh] emits nothing when there are no called diploid samples.
   [Coverage][tier1-coverage] and [RMSMappingQuality][tier1-mq] also emit nothing
   when there is no evidence. For example chr2:144388936 needs `GT:GQ:PL`,
   without AD, sample DP, SB or ExcessHet; chr3:23921573 additionally lacks
   INFO DP and RAW_MQandDP. Missing sample DP uses htsjdk's -1 sentinel, and
   the writer omits absent fields rather than writing zero measurements.
3. **Annotation-time NON_REF likelihood refresh.**
   [HaplotypeCallerGenotypingEngine.prepareReadAlleleLikelihoodsForAnnotation,
   lines 584–623][tier1-refresh] calls `updateNonRefAlleleLikelihoods` when the
   emitted allele count differs from the likelihood matrix. This happens
   after genotyping and before zero-likelihood filtered reads are restored.
   [AlleleLikelihoods, lines 977–1007][tier1-nonref] restricts the median to
   retained alleles while still finding the best likelihood across all
   concrete alleles. Omitting the refresh changed informative-read assignment
   and thus AD/sample DP at chr16:22534296 and chr16:22534345, despite identical
   genotype PLs and QUAL. The refresh is now performed at the same boundary.

The changes to shared files are confined to the GVCF branch of `call_region`
and `annotate_mode`, plus presence-aware AD/DP serialization in `vcf::encode`.
Normal calls still have nonempty AD and nonnegative DP, and take the same
serialization path. `refconfidence.rs` is GVCF-only. Assembly, PairHMM,
genotyping, activity, random draws, phasing and the block combiner are unchanged.
Regression tests cover adaptor boundaries and signed TLEN, no-call annotation
and FORMAT omission, and the retained-allele median with a global best allele.

### Tier1 follow-up validation (2026-10-15)

- Full Agilent v7 BED in **one HC task**, exact tier1 recal BAM:
  `nomics dev gate-gvcf` reports **IDENTICAL**, normalized header identical,
  **1,713,854/1,713,854 records, zero mismatches**. Output:
  `/tmp/nomics-gvcf-tier1-fixed.g.vcf.gz`; gate log:
  `/tmp/nomics-gvcf-tier1-fixed-gate.log`.
- Tier0 rerun: **IDENTICAL**, normalized header identical,
  **251,885/251,885 records, zero mismatches**. Output:
  `/tmp/nomics-gvcf-tier0-fixed.g.vcf.gz`.
- Normal-mode tier0 against the clean pre-GVCF baseline: **IDENTICAL**,
  normalized header and **1,047/1,047 records**. The user's full normal-mode
  tier1/WGS gates were already 18/18 each before this follow-up; those full
  normal-mode runs were not repeated here.
- `cargo test -p nomics-hc -p nomics-gate -p nomics -- --skip tier0_gatk_default_annotations`:
  **43 passed**, nine existing ignored tests, one explicitly skipped external
  fixture test (same absent external fixture described above).
- `cargo clippy -p nomics-hc -p nomics-gate -p nomics --all-targets`: passes
  with existing warnings. Touched Rust files pass rustfmt; `git diff --check`
  passes. No new dependencies or unsafe code.
- Follow-up files: `crates/nomics-hc/src/{refconfidence,call,annotate,vcf}.rs`,
  `crates/nomics-hc/tests/gvcf.rs`, `docs/spec/gvcf.md`, `docs/GATES.md`.
  The stack map was verified byte-identical to its state at follow-up start;
  no publication was attempted.

The small diagnostic BED reproduced the reported loci and exercised each fix,
but also changed downsampling in high-depth regions by omitting earlier RNG
draws. Only the full-BED run above determines the verdict. Its exact match
covers the uninterrupted GVCF reservoir stream; no changes to `jrandom.rs` or
the source-derived draw-order proof were necessary.

[tier1-pileups]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java#L579-L600
[tier1-adaptor]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/utils/read/ReadUtils.java#L1066-L1072
[tier1-ad]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/DepthPerAlleleBySample.java#L39-L57
[tier1-dp]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/DepthPerSampleHC.java#L50-L91
[tier1-sb]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandBiasBySample.java#L62-L86
[tier1-eh]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/ExcessHet.java#L77-L104
[tier1-coverage]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/Coverage.java#L40-L52
[tier1-mq]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/RMSMappingQuality.java#L106-L119
[tier1-refresh]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java#L584-L623
[tier1-nonref]: https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/utils/genotyper/AlleleLikelihoods.java#L977-L1007
