# HaplotypeCaller 4.0.12.0 vs 4.6.2.0: walker, engine and assembly (detail)

Companion to `port-gatk-4.0.12-diff.md` and `port-4.0.12-scope.md`. Scope: `-ERC GVCF` with
default arguments, on the Java code paths. "A" is the 4.0.12.0 tree, "B" is 4.6.2.0; paths are
under `hellbender/`, and PR numbers are from `git log 4.0.12.0..4.6.2.0` (none of them is in
4.0.12.0). Written 2026-10-17 from a source read; items marked [inferred] or [unverified] were
not traced to the end.

## Argument defaults

The same in both tags:
- kmer-size 10,25; min-pruning 2; num-pruning-samples 1; min-dangling-branch-length 4;
  max-num-haplotypes-in-population 128;
- assembly regions 50/300 with padding 100; max-reads-per-alignment-start 50;
- active-probability-threshold 0.002; max-prob-propagation-distance 50; padding-around-snps 20;
- min-base-quality-score 10; pcr-indel-model CONSERVATIVE; gcp 10; mismapping rate 45;
  base-quality-score-threshold 18;
- native-pair-hmm-threads 4; pair-hmm-implementation FASTEST_AVAILABLE;
- GQ bands 1–60, 70, 80, 90, 99; ref-model indel size 10 and deletion quality 30;
  max-mnp-distance 0; heterozygosity priors;
- the read filter list.

Different, or present in only one tag:

| Argument | 4.0.12.0 | 4.6.2.0 | Affects the GVCF default run? |
|---|---|---|---|
| `--smith-waterman` | JAVA | FASTEST_AVAILABLE (#8485) | Yes: 4.0.12 always uses Java SW, even with AVX |
| `--padding-around-indels` | 150 | 75 | Yes |
| `--max-disc-ar-extension` | 25 | renamed, legacy trimmer only | Yes in 4.0.12 |
| `--padding-around-strs` | absent | 75 | Yes (new trimmer) |
| `--standard-min-confidence-threshold-for-calling` | 10.0 | 30.0 | GVCF forces -0.0 in both, but it feeds activity in A |
| `--use-new-qual-calculator` | false | true (#5484) | Yes (genotyping) |
| HMM normalisation | asymmetric | symmetric (#6196) | Yes |
| Pre-filter read transformer | identity | `IUPACReadTransformer` (#6868) | Only for non-ACGTN bases |

Both tags force, in GVCF mode: calling confidence -0.0, PLs at every site, physical phasing,
the GVCF writer with those bands. QD/FS/SOR/ChromosomeCounts are dropped and
StrandBiasBySample is added, so the shared `Random(47382911)` is drawn only by the downsampler.

## Output-affecting changes

Port size: S = small, M = medium, L = large.

| # | Change | A vs B | Affects output | Port |
|---|---|---|---|---|
| 1 | Activity probability (isActive), #6343 | A: heterozygosity priors with Jacobian log-sum and a `min(4, conf)` cutoff (0.0 in GVCF). B: Dirichlet AF calculator. | Yes, region boundaries everywhere | S |
| 2 | finalizeRegion order | A hard-clips low-quality ends first, evaluates fragment size on the trimmed read, keeps reads overlapping the region after clipping. B tests overlap before region clipping. Feeds the reference model. | Conditional (low-quality tails in soft clips, edge reads) | S |
| 3 | Overlapping-mate quality capping (#6886, #5794) | Different coordinate mapping when a mate starts in a deletion; A keeps the given order on equal soft starts. | Conditional | S |
| 4 | Trimming and genotyping window (#6358, #6507) | A: one padding per region (150 if any non-SNP event, else 20), symmetric extension in GVCF capped at 100. B: per-event padding (20/75/75+STR), asymmetric. B's `trimLegacy` is not A. | Yes, reads, haplotype span and GVCF block boundaries | M |
| 5 | Haplotype.trim (#6518, #6544) | A drops haplotypes whose trimmed CIGAR starts or ends with I/D; B strips leading/trailing I and keeps D. | Conditional | S |
| 6 | Graph acceptance order (#5562) | A checks cycles and low complexity before pruning; B prunes first. | Yes, which kmer sizes yield graphs and when kmer escalation triggers | S |
| 7 | Dangling-end recovery (#7086) | B adds a CIGAR offset when choosing the ref node, aborts on loops, and takes the heaviest outgoing edge. | Conditional | S |
| 8 | K-best haplotype finder (#5462, #5554, #5786, #5952, #6104) | A: lazy recursive per-vertex finder, edge score from `max(0.5, m)`, ranking in a `java.util.PriorityQueue` on score only, so ties follow heap order and edge insertion order. B breaks ties by bases. | Yes, with more than 128 paths or tied scores | L |
| 9 | findBestPaths / calculateCigar / left-align (#6885, #6427, #6510, #6403) | A has an equal-length ≤2-mismatch all-M shortcut and per-indel `leftAlignSingleIndel` with the old trim/consolidate helpers. | Yes, multi-indel and adjacent I/D haplotypes | M |
| 10 | EventMap (#8567, #6661, #6696) | B trims every event to minimal representation, including compound events. | Conditional | S |
| 12 | Likelihood capping (#6196) | A caps alleles at best non-ref − 4.5; B at overall best − 4.5. PairHMM kernels are byte-identical. | Yes, PL and GQ | S |
| 13 | Realignment to best haplotype (#6461, #6498) | A: single-indel left-align, old CIGAR helpers, leading-deletion shift. | Yes, annotations and ref-confidence pileups | M |
| 14 | Clipping engine (#6280, #6403) | A shifts starts wrongly for indels inside the clipped part (the #6280 bug). | Conditional | M |
| 15 | SW CIGAR consolidation | A merges adjacent equal ops; equivalent downstream [inferred]. | No [inferred] | S |
| 16 | Graph set iteration (#6195) | A returns `HashSet`s and `SeqVertex.hashCode` is `System.identityHashCode`; they feed the diamond/tail/suffix mergers and k-best edge order. | Conditional and possibly **non-deterministic** in A | Policy decision |

No output change from downsampling: same seed and draws, and end-of-input finalisation
(#5594) does not change which reads are emitted or their order. The rest is refactoring:
- region iteration, LIBS and span renames;
- `CigarUtils.isGood`;
- soft-clip reverts;
- identical threading, pruners and SeqGraph transformers;
- `ReadLikelihoods` renamed to `AlleleLikelihoods`;
- identical MathUtils tables.

## Unknowns

- **Determinism of 4.0.12.0** (item 16 plus #6104/#6195). Byte-exact output may be impossible on
  regions where identity-hash order decides graph simplification or k-best ties. Measure r1 vs
  r2 before committing scope.
- **jgrapht 0.9.1** (A) vs 1.1.0 (B) edge-set ordering, which decides k-best tie insertion order
  [unverified].
- **Not read here:** the GVCF emission path. That covers `ReferenceConfidenceModel.calculateRefConfidence`
  (a 529-line diff), indel-informative read counting (#5469, #5607), `getIndelPLs` and the
  GVCF writer. Genotyping and annotations are covered in `port-gatk-4.0.12-diff.md`.
- `PriorityQueue` tie behaviour in A's aggregated finder was not traced line by line
  [unverified].

The raw `@Argument` dumps and commit list were kept outside the repository (`hcdiff/`), which is
not durable.
