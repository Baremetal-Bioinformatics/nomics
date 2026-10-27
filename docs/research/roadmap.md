# After the 30× WGS gate: speed and scope (research pass, 2026-10-11)

Research only; nothing here is implemented. Sources: our own measurements (`docs/research/alignment.md`, `bench/`), GATK 4.6.2.0 and bwa 0.7.19 source, bwa-mem2 / bwa-mem3 code and docs, nf-core/sarek 3.10.0 config. Claims are tagged **[measured]** (ours), **[source]** (read in the reference code), **[claimed]** (paper or vendor, unverified), **[estimate]**.

## Where a 30× genome spends its time

Projection from the 10M-pair genome-wide rehearsal and the 30× chr20 slice (2026-10-11), Spark, whole-genome mode:

| Stage | Projected | Notes |
|---|---:|---|
| prep (spill) | ~6 min | single-threaded inflate + filter |
| alignment | ~1.9 h | per WGS read, single thread: rescue 91 µs, extend 81, chain 22, chain filter 7, seed 60 (on the GPU) **[measured]** |
| markdup, BQSR (streamed) | ~25 min | |
| HaplotypeCaller | ~1.6 h | PairHMM is ~80 % of thread time; centromeric pileups dominate **[measured]** |
| total | ~4 h | |

## Speed: exact levers, ranked

Every lever must leave output bit-identical; each is gated 18/18 on tier1 and wgs20, plus the 30× digests.

### 1. PairHMM: skip shared-prefix columns (HaplotypeCaller, est. 1.6–2× on HC)

- **[source]** GATK 4.6.2.0 has prefix-reuse machinery (`hapStartIndex`, `findFirstPositionWhereHaplotypesDiffer`) but it is dead code: `computeLog10Likelihoods` passes `recacheReadValues = true` for every haplotype, so GATK recomputes every cell. Our port does the same.
- Cell (i, j) depends only on columns ≤ j, the haplotype bases before j, and per-read constants; row 0 of D is `2^1020 / hapLen`. So for **equal-length** haplotypes sharing a prefix of length P, columns 1..P and the running final sum up to P are bit-identical. Skipping them is exact, not a reordering.
- **[measured]** on the PairHMM sample dump: **59 % of cells** sit in shared prefixes of equal-length haplotypes. (Exact duplicate reads: only 1.5 %; not worth it.)
- Shape: lanes become 8 *reads* (each with its own row parameters) instead of 8 haplotypes; walk a per-length haplotype trie column by column, saving the column state (all rows) at branch points. Each cell is still computed with Java's operand order.
- Risk: medium (new traversal); the bitwise PairHMM bench against the current lanes is a ready oracle.

### 2. Mate rescue: exact pruning and banding (alignment, est. 20–30 % of WGS alignment)

- **[source/claimed]** bwa-mem3 (Fulcrum) computes fewer cells per rescue call without skipping calls: exact 5-mer hits between mate and window bound the best possible score, only rows that can reach `min_seed_len·a` are computed, the first pass runs in diagonal bands around the hits, and the start-recovery pass runs in a band from the first pass's score and end. Stated preconditions cover bwa's default scoring. Reported 10–22 % wall on Graviton4 NEON for the scoring extensions in PR #541 (no default-scoring number found); gains are largest on narrow SIMD, which is our case.
- Our windows are ~1,339 rows for a 151 bp mate: **[estimate]** 3–6× fewer cells, against 91 µs/read of rescue.
- The proof must also reproduce `score2` (bwa's per-row `b[]` plateau merging in `ksw_u8`, which bwa-mem2's batched kernels got wrong: XS ±1–4, MAPQ ±2–18) and the qb/tb start pass.
- Risk: medium-high. Gate every call against the full-window kernel on the 62k exome and 70k WGS dumps before any pipeline gate.

### 3. HaplotypeCaller: run interval groups concurrently (est. ~10 % on HC)

- Each of sarek's 21 groups starts with fresh RNG state, so groups are independent **[source, ours]**. Today they run one after another; within a group the pool is ~91 % busy on 20 threads (chr20 centromere group: 1,834 user-s over 101 s) **[measured]**, so the gain is the tails, not throughput.
- In whole-genome mode this needs two groups' records in memory (~14 GB each at 30×).
- Risk: low.

### 4. Small certified alignment shortcuts (est. 5–10 %)

From bwa-mem3, each exact by construction or with a replay fallback: repeat-job reuse in rescue (byte-identical repeats of the thread's last few jobs), skipping the post-rescue dedup when provably a no-op, a memo for `cal_max_gap`, an index over kept seeds for the containment sweep, a certified band for the final CIGAR alignment, lockstep SMEM batching for the CPU-seeded overflow.

### 5. Chaining: dense SA interval 1 (chain 22 → ~10 µs/read)

- **[measured]** interval 1 halves chaining (15.5 → 6.9 µs/read on the exome), exact, but GPU seeding failed with CUDA error 400 while the 30 GB SA was resident (`alignment.md` idea 7).
- **Update 2026-10-13 [measured]:** interval 2 is now the run scripts' default. On the 30× chunk, chaining went 42 → 27.5 µs/record and the chunk 430 → 414.5 s. Interval 1 ran (no CUDA error) at 18.7 µs chaining but 475.6 s overall and a 77.5 GB peak. That run had only 79% huge-page coverage, so call it slower under those conditions rather than closed; a 38 GB index is impractical anyway.
- Memory: fine at whole-genome peaks (~65 GB projected) only if the CUDA 400 cause is found.

### 6. GPU: keep it to seeding

- GB10's FP64 rate is unpublished; if consumer-class (1:64), its ceiling is about our 20-core NEON PairHMM rate, while sharing power and LPDDR bandwidth with the CPU **[estimate]**. Measure a pure FP64 kernel once; drop the GPU PairHMM share unless it delivers >25 % of the CPU rate.
- **[claimed]** GPU SW/extension speedups are 2–2.5× even on an A100 against 64 cores. Our first GPU rescue measured 3e9 cells/s against 2.5e10 on the CPU. **Update 2026-10-13 [measured]:** a warp-wavefront `ksw_align2` with shared-memory staging runs at 2.4–2.7e11 cells/s, exact, and is now the default (30× align −17%). A row-per-warp `ksw_extend2` is exact but neutral end to end (opt-in; `alignment.md`). Seeding (random access, latency-bound) stays the GPU's job; batched SA lookups are the only other candidate.

### Ruled out or to avoid

| Idea | Why |
|---|---|
| NEON extension within a row | exact, but no faster than scalar on these cores **[measured]** |
| Interleaving two rescue alignments | exact, 1.45× on Apple M but 1.08× on the Spark **[measured]** |
| SIMD extension across reads; k-mer seed table | measured earlier (`alignment.md`); GPU mate rescue was here and is now kept (2026-10-13) |
| Float-first PairHMM with double fallback (GKL) | changes bits **[source]** |
| Seed reordering, tie-break changes, bwa-mem2's dense row-max `score2`, heuristic rescue skips, a different `-K` | change output |

## Scope: beyond the 18 gates

Our BAMs and VCF are identical to sarek's, so any downstream tool run unmodified on them (VEP, snpEff, MultiQC, bcftools) gives sarek's result; porting such a tool buys speed only. Every item below needs a new sarek r1/r2 golden pair first.

| # | Item | Effort | Why |
|---|---|---|---|
| 1 | fastp poly-G trimming for two-colour instruments | S–M | sarek never disables it and fastp enables it from NovaSeq/NextSeq read names; `nomics-prep` refuses such input today, so most users' own FASTQs fail on day one. Golden: HG002 re-headered to instrument names. |
| 2 | Multiple lanes / read groups, Picard optical duplicates | M | `rg_id` is hard-coded `{prefix}.L1` (`pipeline.rs:417`); optical duplicates are never counted (real read names switch them on in Picard). |
| 3 | `.bai` (then `.crai`) | S–M | already a stated, unproduced gate in `GATES.md`; every viewer needs it. |
| 4 | CRAM output (sarek's default), record-identical | M (L byte-exact) | what users archive; byte-exact needs htslib's codec choices. |
| 5 | QC files: samtools stats, mosdepth, bcftools stats, vcftools, FastQC data | M total | on by default in sarek, deterministic, cheap to gate; emit MultiQC-ingestible files, not MultiQC's HTML. |
| 6 | GVCF emit, then single-sample GenotypeGVCFs | L + M | cohorts and reanalysis. New subsystems: hom-ref blocks with GQ bands, `<NON_REF>`, raw annotation forms, and physical phasing (forced on in GVCF mode). VQSR (seeded Gaussian mixture, `Random(47382911)`): deterministic, so a bit-exact port is feasible. It is not deprecated in GATK; Illumina dropped it from DRAGEN. |
| 7 | Reference flexibility: no-BED interval generation, BQSR-skip path, ALT contigs | M | GRCh37, T2T (no known sites), GRCh38 with `.alt` (refused today, `index.rs`). |
| 8 | Mutect2 tumor–normal | L+ | separate track; reuses assembly, PairHMM and read prep; needs tumor–normal truth data. |

Non-goals: DeepVariant, the CNN filter (`haplotypecaller_filter`, TensorFlow), VEP/snpEff internals. GATK hard filters could be offered, labelled as not sarek.

## Suggested order after the 30× gate

1. **PairHMM shared-prefix traversal**: largest exact speed lever left, with a ready bitwise oracle.
2. **Poly-G and multi-lane/read-group support**: what stands between real users' FASTQs and the pipeline.
3. **Exact rescue pruning**: largest alignment lever, higher proof burden.
4. **`.bai`, CRAM, QC files**: completes sarek's default outputs.
5. **GVCF**: opens cohort use; the biggest new subsystem.

## End of the 2026-10-12/14 iteration window

30× genome: 2 h 47 min → 1 h 59 min (final run 7,135 s; earlier near-final builds 7,186 and 7,281 s), 41/41 identical. What was kept, what was rejected and the measurements behind each are in `bench/` and `docs/strategy/`.

Landed:
- levers 1, 3 and 5;
- GPU mate rescue;
- pipelined BGZF reader and writer;
- markdup merge/sort/classify parallelism;
- prep overlapped with alignment (inputs memory-mapped);
- multi-lane input, poly-G trimming, optical duplicates.

Next session starts from the velocity cards (V-02..V-09): digest caching, a 30× one-chunk slice tier, batched gates, more dump oracles, an x86 second lane, parallel candidates, timing-noise fixes, and surviving a full page cache. Open performance items:
- GPU extension kernel 3×;
- HC overlap with three groups at lower memory (−7% more HC);
- markdup pass 1's serial pairing sweep (now 133 s);
- the index-paging root cause.
