# Porting the UK Biobank 500k WGS pipeline (BWA 0.7.17 + GATK 4.0.12.0): scope and goldens

Scoping note, 2026-10-17. Analysis only. It covers the pipeline steps, BWA 0.7.17 vs 0.7.19, and how to
make goldens. The GATK 4.0.12.0 vs 4.6.2.0 diff (BQSR, ApplyBQSR, HaplotypeCaller, GenotypeGVCFs
internals) is being done by another worker and is **excluded** here. "Unverified" means we could not
retrieve a primary text, or the claim is our inference.

Sources:
- **[SI]**: UKB WGS Nature 2025 supplement ([docx](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09272-9/MediaObjects/41586_2025_9272_MOESM1_ESM.docx)), Supplementary Note 3.
- **[SB]**: Seven Bridges tables for the Sanger workflows, which the SI links ([UKBseq500k-methods/sb-ukb-tables](https://github.com/UKBseq500k-methods/sb-ukb-tables), PDF, Tables S1, S4, S5).
- **[FE]**: CCDG functional-equivalence standard ([PipelineStandard.md](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md), last commit 2017-03-25), and [Regier 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6168605).
- **[bwa]**: `git diff v0.7.17 v0.7.19` of github.com/lh3/bwa.

## 1. What the UKB WGS pipeline runs

There are three pipelines: deCODE, Sanger Vanguard and Sanger Main Phase. 247,522 samples were sequenced
by deCODE, 193,093 by Sanger, and 49,934 in the Vanguard phase [SI §"Sequence data processing"]. All
three claim FE compliance [SI]. Every centre's per-sample product is:
- a CRAM "containing original instrument basecall quality values", with its CRAI;
- a GATK gVCF;
- the BQSR table [SI].

Joint calling across the cohort is **not GATK**. It is GraphTyper v2.7.5 `genotype`, run on all CRAMs
against `GRCh38_full_analysis_set_plus_decoy_hla` [SI, "SNP and indel calling with GraphTyper"]. A
separate DRAGEN 3.7.8 release exists, which we do not port. GATK GenotypeGVCFs runs **per sample only**:
- deCODE uses it for QC and sex checks [SI].
- Sanger writes a single-sample VCF with it [SB S1/S5].

### 1a. Step by step

| Step | FE standard [FE] | deCODE [SI "deCODE pipeline"] | Sanger Main Phase [SI, SB S4/S5] (Vanguard [SB S1] where it differs) |
|---|---|---|---|
| Input | — | bcl2fastq; one FASTQ pair per (lane, index) read group | bambi CRAM → biobambam2 `bamtofastq` 2.0.144, per read group |
| Reference | GRCh38DH (1000G): primary + alts + decoys + HLA | "GRCh38 with alt contigs plus decoy contigs and HLA"; GraphTyper uses `GRCh38_full_analysis_set_plus_decoy_hla` | `GRCh38_primary_assembly_plus_ebv_alt_decoy_hla.fasta`, ".alt file supplied alongside it" [SB] |
| Align | BWA-MEM 0.7.15, `-K 100000000 -Y`, no `-M`; ".alt file … do not perform post-processing of alternate alignments" | `bwa mem` 0.7.17 `-K 100000000 -Y -t 24` | `bwa mem` 0.7.17 `-K 100000000 -Y -R … -t $NT` piped to `sambamba view -f bam` (0.6.7) |
| Mate tags | Keep RG, MQ, MC, SA; suggests `samblaster -a --addMateTags` | samblaster 0.1.24 `-a --addMateTags`, "to add MC and MQ tags". The step outline also lists "Fix mate pair information" between align and mark duplicates; whether a separate fixmate tool ran is **unverified** | **Omitted** "in agreement with WSI" [SB S4]. Vanguard S1 doesn't list samblaster either |
| Duplicates | Match Picard; Picard ≥2.4.1 on queryname input | Picard MarkDuplicates 2.20.3, `ASSUME_SORT_ORDER=queryname READ_NAME_REGEX=<custom>`, **per read group** | Picard MarkDuplicates 2.18.26, all RG BAMs as `INPUT`s, `OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname VALIDATION_STRINGENCY=SILENT` |
| Sort | — | samtools 1.9 sort | sambamba 0.5.9 sort |
| Merge, CRAM | Lossless CRAM; samtools ≥1.3.1; no OQ | QC per RG, then `samtools merge` 1.9 → CRAM | `samtools view --output-fmt CRAM` 1.9; `samtools index -c` |
| BQSR table | dbSNP138, Mills, known_indels; `-L chr1..chr22` optional | BaseRecalibrator 4.0.12, same three sites, `--preserve-qscores-less-than 6 -L chr1 … -L chr22` (reconstructed from garbled text) | BaseRecalibrator 4.0.12.0, same sites, `--intervals interval.bed --use-original-qualities`, then GatherBQSRReports |
| ApplyBQSR | `-SQQ 10 -SQQ 20 -SQQ 30 --disable_indel_quals` (4-bin: 2–6, 10, 20, 30) | per chromosome chr1–22, X, Y. Flags as extracted: "--preserve-qscores-less-than --static-quantized-quals 1… 2… 3… --create-output-bam-ind…", **garbled** (see below) | per interval: `--static-quantized-quals 10/20/30 --add-output-sam-program-record --use-original-qualities` |
| HC | (calling was centralised for FE) | HaplotypeCaller 4.0.12 `-ERC GVCF` per chromosome (24), then Picard MergeVcfs 2.20.3 | HaplotypeCaller 4.0.12.0 `--emit-ref-confidence GVCF --intervals interval.bed`; tabix 0.2.6 |
| Per-sample genotyping | — | GenotypeGVCFs 4.0.12 (QC) | GenotypeGVCFs `--dbsnp dbsnp138 --only-output-calls-starting-in-intervals -G StandardAnnotation --use-new-qual-calculator`, then MergeVcfs |
| Joint calling | — | GraphTyper 2.7.5 (cohort) | same cohort set |

Answers to the specific questions:

- **`-K 100000000 -Y`?** Yes, at both centres. No `-M` [SI, SB].
- **Alt-aware?** Yes. Sanger says outright that the `.alt` file was supplied with the index [SB S1/S4]. deCODE says only "with alt contigs", and FE requires the `.alt`. deCODE having it is **unverified but likely**. **No `bwa-postalt`** anywhere: FE forbids it, and neither centre lists it.
- **Duplicates:** Picard MarkDuplicates in queryname mode at both centres (2.20.3 at deCODE, 2.18.26 at Sanger). samblaster only adds mate tags (`-a` means `--acceptDupMarks`: it marks nothing; samblaster.cpp v.0.1.24 L1274). No bamsormadup.
- **BQSR:** GATK 4.0.12.0, not 3.x. Known sites are the FE trio from the GATK hg38 bundle.
- **Binning:** the 4-bin SQQ 10/20/30 scheme. In the SI's account, the binned recal BAM is an **intermediate that feeds HC**, and the delivered CRAM keeps the instrument qualities [SI]. **Unverified inconsistency:** Vanguard Table S1 orders ApplyBQSR → GatherBamFiles → `samtools view` CRAM, which would make the Vanguard CRAM recalibrated and binned. This changes what a port must emit for Vanguard samples.
- **Garbled extraction (unverified):** the SI's `READ_NAME_REGEX` and the deCODE ApplyBQSR flags are mangled in the docx; the same text came from a direct XML extract and from Firecrawl. From the Sanger table and FE, the ApplyBQSR line is presumably `--preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --create-output-bam-index`. The regex decides which reads count as optical duplicates (which feeds the duplicate metrics), so the exact string must come from deCODE.
- **Where deCODE and Sanger differ:**
  - samblaster (deCODE only);
  - Picard version (2.20.3 vs 2.18.26);
  - MarkDuplicates per lane (deCODE) vs per sample (Sanger);
  - optical pixel distance (default 100 vs 2500);
  - sorter (samtools 1.9 vs sambamba 0.5.9);
  - scatter (24 chromosomes vs an SB interval file);
  - `--use-original-qualities` (Sanger only). It reads the OQ tag, but Sanger realigns from `bamtofastq` output, which has no OQ, so the flag should be a no-op (inferred) and there is nothing to port;
  - reference FASTA label.

### 1b. What a nomics port must reproduce, boundary by boundary

| Boundary | deCODE | Sanger | Have it? |
|---|---|---|---|
| Raw `bwa mem` SAM per RG FASTQ (no fastp; `-K` batches run over the whole RG file) | yes | yes | Port is 0.7.19, needs a 0.7.17 mode (§2). `nomics run`/`prep` assume fastp chunks (`crates/nomics/src/main.rs` `Prep`, `Run`), so we need a raw-FASTQ align entry |
| Alt-aware alignment | yes | yes | **No.** `index.rs` L521 refuses `<prefix>.alt` |
| samblaster mate tags | yes | no | **No** |
| Picard MarkDuplicates, queryname mode | 2.20.3, per RG | 2.18.26, per sample | **No.** `nomics-markdup` is Picard 3.4.0, coordinate path only (`docs/spec/markdup.md` L553) |
| Coordinate sort tie order | samtools 1.9 | sambamba 0.5.9 | **No.** We have samtools 1.22.1 sort. Tie order feeds HC downsampling |
| `samtools merge` 1.9 of RG BAMs | yes | no | Partly (1.24 behaviour is documented in REFERENCE-VARIABILITY) |
| CRAM 3.0 output | yes | yes | **No CRAM code anywhere** in `crates/` |
| BQSR / ApplyBQSR (SQQ, preserve-qscores, OQ) / HC GVCF / GenotypeGVCFs | yes | yes | Other worker (GATK 4.0.12.0) |
| Scatter: 24 chromosomes | yes | — | Easy; note there is no chrM, alt or decoy in HC |
| Scatter: SB `interval.bed` from `sbg_prepare_intervals.py --mode 4` on `Homo_sapiens_primary_assembly_38_81_intervals.bed` | — | yes | **Blocked.** The file and the script are not public (searched; not found) |
| GraphTyper 2.7.5 joint genotyping | cohort | cohort | **Not ported**; recommend out of scope |

Tools we don't have, and our recommendations:
- **samblaster:** port it (small).
- **Picard 2.x queryname MarkDuplicates:** port it (large).
- **sambamba sort:** port only if we take on the Sanger profile.
- **CRAM:** gate on records via `samtools view` as a test oracle. GATES.md allows that category. We would not write a CRAM reader or writer yet.
- **GraphTyper:** do not port.
- **bwa-postalt:** not needed.

**Recommendation: port the deCODE profile first.** It is fully specified apart from the regex and has the
largest sample share. The Sanger profile is blocked on its interval list, and HC and BQSR depend on
scatter boundaries.

## 2. BWA 0.7.17 vs 0.7.19

Tags `v0.7.17` (`0.7.17-r1188`) and `v0.7.19` (`0.7.19-r1273`) exist. Between them are 83 commits and 32
files (+679/−88). Most of that is licence headers, CI, typos, `bwa aln` and `xa2multi.pl`. The changes
that can alter `bwa mem -K 100000000 -Y` output:

| Change | Commit | Effect on our flags |
|---|---|---|
| `MQ:i:<mate mapq>` written whenever a mate exists, after `MC:Z`, before `AS:i` (`bwamem.c` `mem_aln2sam`) | d8dd308 (2021-07-21) | **Changes every paired record.** 0.7.17 writes no MQ. Our port writes it (`sam.rs` L457-458), so 0.7.17 mode drops it |
| `@HD\tVN:1.5\tSO:unsorted\tGO:query` printed when `-H` has none (`bwa.c` `bwa_print_sam_hdr`) | 6b18630, a2cf60d, c232375 | Header only. 0.7.17 prints no `@HD`. Any later sort writes its own |
| `-u` (XB tag) and `-z` (XA drop ratio) options | e050ec9, 50122be | None. Both are off by default, and `XA_drop_ratio = 0.80` was already in 0.7.17 (`bwamem.c` L71) |
| `ksw.c` NEON / scalar SSE2 emulation | 165e524, ac612b6, ab01ab4 | x86 output unchanged. **0.7.17 has no arm64 path** |
| `CFLAGS -O2` → `-O3` | Makefile | No FP reassociation without `-ffast-math`. x86 baseline has no FMA. On arm64, contraction is the existing 58-site story (REFERENCE-VARIABILITY) |
| `kt_for` int → long; `smem_aux_t*` malloc size fix | 49e854e, c45d772 | None (memory and ABI only) |
| `bns_restore`: abort on `.alt` names ≥1023 chars | 20d0a13 | None for valid `.alt` files |
| `bwa-postalt.js` skips empty lines | 27dd1da | Not run |
| `bwaseqio.c` strlen → `name.l` | 973ac5f | `bwa aln` only |
| `@PG VN:` string | release commits | Header only |

`bwamem_pair.c` and `bntseq.c` changes, apart from the `.alt` check, are licence headers. There are no
changes to seeding, chaining, extension, pairing, MAPQ, XA/SA or `-K` batching. **So a 0.7.17 mode
is: no MQ, no `@HD`, a different `@PG` version.**

bwa 0.7.17 already writes `MC:Z` (it is a context line in the diff); only `MQ:i` is new. deCODE's samblaster then appends `MC:Z` and `MQ:i`
unconditionally, with no check for existing tags (samblaster.cpp v.0.1.24 L903-928, `addTag` L340).
deCODE records may therefore carry **two MC tags**. How Picard 2.20.3 and htsjdk treat that is
**unverified**, and is the first check to run on the deCODE golden.

**Alt-aware mode inside bwa.** With `<prefix>.alt` present, `bns_restore` sets `anns[i].is_alt`. The
0.7.17 source then changes these behaviours (`git grep is_alt v0.7.17`; about 25 sites):
- chain overlap filtering (`bwamem.c` L351);
- the `mem_ars_hash`/`hash2` sort keys (L397/400);
- `mem_mark_primary_se`. Primaries and `sub`, and with them **MAPQ**, are computed among non-ALT hits only, and ALT hits get `secondary = INT_MAX` (L523-557);
- `alt_sc`, giving a `pa:f` tag (L932);
- XA: up to 200 hits when any hit is ALT (`bwamem_extra.c` L110-118);
- `AH:*` on ALT `@SQ` lines (`bwa.c` L394);
- the soft-clip and supplementary-MAPQ exceptions (L818, L875, L1029);
- ALT handling in pairing (`bwamem_pair.c` L326-377).

Without `.alt` all of these paths are dead, which is the case for our current reference
(`docs/spec/bwa.md` §0.2). **nomics already has most of them in code:**
- `is_alt`, `n_pri`, `alt_sc`, `has_alt`/`MAX_XA_HITS_ALT` and `pa:f` are in `chain.rs`, `extend.rs`, `pair.rs` and `sam.rs`;
- but `Bns::load` hard-codes `is_alt: false` (`index.rs` L387) and the loader refuses `.alt` files (L521);
- there is no `AH:*` output.

The remaining work is: parse `.alt` the way `bns_restore` does, emit `AH:*`, and gate the newly live
paths against a 0.7.17 oracle on reads that hit alternate contigs (MHC/HLA, chr6/chr19 alts).

**Reference.** Our iGenomes `Homo_sapiens_assembly38.fasta.64` index has the same 3366 contigs
(`bwa.md` L672) but no `.alt`. The GATK hg38 bundle ships `Homo_sapiens_assembly38.fasta.64.alt`
(**unverified** that it matches deCODE's `.alt`). The FASTA names differ:
- deCODE: `GRCh38_full_analysis_set_plus_decoy_hla`, the 1000G GRCh38DH;
- Sanger: `GRCh38_primary_assembly_plus_ebv_alt_decoy_hla`.

Contig order and per-contig MD5 must be checked against each, because `.ann` order sets `rid`
tie-breaks and `@SQ` order. Whether the 1000G analysis set and Broad's assembly38 are identical
(chrY PAR masking) is **unverified**.

**Effort, given the exact 0.7.19 port:**
- 0.7.17 mode: about 1 day, including building a 0.7.17 oracle with `reference/oracle/instrument_bwa.py` and re-running the tier0 gates.
- Alt mode: 4–10 days, mostly oracle stage-dump gating on alt-heavy reads.

For scale: the 0.7.19 port went from scaffold (2026-10-08) to gated within days (git log of
`crates/nomics-align`).

## 3. Goldens

### 3a. Binaries and containers

| Tool | x86-64 | arm64 |
|---|---|---|
| bwa 0.7.17 | `quay.io/biocontainers/bwa:0.7.17--he4a0461_11` (also `hed695b0_7`). Pin the build string, because compilers differ by build | **No bioconda linux-aarch64 build** (0.7.18 is the first). Build from source with sse2neon or a backport of 0.7.18's `ksw.c` shim. That makes it "0.7.17-modified", which must be shown record-identical to x86 on the slice |
| samblaster 0.1.24 | biocontainers `0.1.24--hc9558a2_3` | No aarch64 build; build from source (plain C++) |
| samtools 1.9 | biocontainers (linux-64 only) | Build from source |
| sambamba 0.5.9 / 0.6.7 | biocontainers `0.5.9--1` | None (D toolchain); x86 only |
| Picard 2.20.3 / 2.18.26 | bioconda noarch jar | Same jar |
| GATK 4.0.12.0 | `broadinstitute/gatk:4.0.12.0` (amd64 only, 1.42 GB, OpenJDK 8 base); bioconda `gatk4-4.0.12.0-0` (noarch) | Bioconda noarch jar on a JDK 8 aarch64 (e.g. Temurin 8) |

Verified from `build.gradle` at tag 4.0.12.0: `requiredJavaVersion = "8"`, GKL 0.8.6, bundled Picard
2.18.16, htsjdk 2.18.1. The UKB centres ran standalone Picard 2.18.26/2.20.3, not the bundled version.

**Does GATK 4.0.12.0 run on arm64?** The Java code is portable, but GKL 0.8.6 ships x86-64 native
libraries only. The current arm64 setup already depends on GATK falling back to Java code when those
libraries are missing (REFERENCE-SPEC, DATA-MANIFEST), and that fallback is in place in the 4.0 line.
Whether it is clean in 4.0.12.0 is **unverified**: in particular, whether `--smith-waterman` exists in
4.0.12.0 and which deflater is used. Test it on the Spark with `--pair-hmm-implementation
LOGLESS_CACHING` and check the log.

**JDK 8 risk.** The repo's x86-Java = arm64-Java equivalence relies on `-XX:-UseLibmIntrinsic` under
JDK 17 (REFERENCE-VARIABILITY, `reference/cloud-hcjava.config`). JDK 8 x86 HotSpot has different `Math`
intrinsics, and the flag may not exist. Check with `java -XX:+PrintFlagsFinal -version | grep -i libm`
inside the 4.0.12.0 image. Until then the "portable Java-mode golden" is **unverified for JDK 8**.

### 3b. Sample

Use HG002 30× (SRR14724532; the `wgs30` goldens) and its 1/12 slice (`wgs30s`, `samplesheet.cloud-wgs30s.csv`).
- **Read groups.** UKB aligns per (lane, index) read group with no fastp. Each lane becomes one raw FASTQ pair, and with `-K 100000000` the batch boundaries depend on that file. To exercise deCODE's per-RG MarkDuplicates followed by `samtools merge`, split by lane if the names allow it; otherwise use one RG.
- **Read names.** SRA names carry no `tile:x:y`, so optical-duplicate detection never fires (`docs/spec/markdup.md` L226-228). Neither deCODE's regex nor Sanger's pixel distance of 2500 would be tested. To cover that path, add a small FASTQ with original Illumina names. Candidates are GIAB HG002 NovaSeq runs; availability is **unverified**.

### 3c. Plan

1. **Write `reference/ukb/run_ukb_decode.sh`**, a plain shell driver (not sarek), using pinned image digests and the commands in §1a. It should:
   - `tee` every intermediate that production discards: raw bwa SAM; post-samblaster SAM; MarkDuplicates BAM and metrics; sorted BAM; merged pre-CRAM BAM; CRAM; per-chromosome recal BAMs; BQSR table; 24 gVCFs and the merged gVCF; per-sample VCF;
   - record argv, `@PG` lines, image digests, start seconds (samtools merge 1.9 seeds its ID suffixes from the clock) and `/usr/bin/time -v`.
   - Also write `reference/ukb/build_ref.sh` to fetch `Homo_sapiens_assembly38.fasta.64.alt` and record contig MD5s, and extend `compare_runs.sh` with CRAM-by-records (samtools as an oracle).
2. **Slice r1 and r2** on the x86 box, HC in Java mode (the same flags as `cloud-hcjava.config`, adapted for JDK 8). Measure r1-vs-r2 at each boundary first, then write the gates (GATES.md method).
3. **GKL-mode HC-only rerun** on the same recal BAMs, with the default PairHMM (AVX). This records what labs actually ran, as a concordance number. It is not a gate, because GKL output is not portable (REFERENCE-VARIABILITY, 30× native vs Java). Optionally repeat on an AVX2-only instance type to see whether the GKL output depends on the CPU.
4. **Full 30× r1 and r2**, staggered on one box. Start r2 when r1 reaches HC: per-chromosome HC makes chr1 a long single-threaded tail, which leaves idle cores for r2's alignment.
5. **arm64 cross-check on the Spark** (under `flock`): the bioconda GATK jar on JDK 8, and the from-source 0.7.17-modified bwa. Run the slice plus the two heaviest chromosomes, compared against x86 Java mode.
6. **Sanger profile** (`run_ukb_sanger.sh`, with sambamba) only after the interval list is in hand.

**Compute and cost.** The brief's "32-vCPU box, ~4.8 h, ~$5" figure is the **nomics** CPU-only 30×
run: 17,187 s, of which 13,587 s was alignment (`docs/research/x86-lane-r1/stages.tsv`). It is not a
sarek run. No sarek 30× x86 wall time is recorded in the repo docs. That implies about $1/h for the box.

Projections (**unverified**):
- bwa 0.7.17: 6–8 h on 31 threads. This is an extrapolation. nomics CPU alignment of the same 700M reads on the same box took 13,587 s. `bench/align.md` puts nomics at ~13k reads/s/core and bwa at ~7.4k (tier0), so bwa should take about 1.76× as long, ≈ 6.6 h.
- Picard queryname MarkDuplicates: 2–4 h.
- Sort and merge: about 1 h.
- BQSR and ApplyBQSR: about 2 h.
- Per-chromosome Java-mode HC: dominated by chr1 at about 8% of the genome, single-threaded, plausibly 8–16 h.

| Item | Wall | Cost at ~$1/h (32 vCPU) |
|---|---|---|
| Slice r1+r2 | 3–5 h | $3–5 |
| Full 30× r1+r2 (staggered on one box; one run ≈ 20–30 h) | 30–45 h | $30–45 |
| GKL HC rerun (full) | 3–6 h | $3–6 |
| arm64 cross-check | Spark time | $0 |
| Sanger profile, later (same shape) | 30–45 h | $30–45 |
| **Total** | | **≈ $65–100**, plus storage (~0.5 TB of intermediates per 30× run, unverified) |

Instance prices are not re-checked.

## Risks

1. **There is no UKB output to compare against.** HG002 is not a UKB sample, so the goldens reproduce the *published description* of the pipeline, never a UKB file. Every unverified item below compounds this.
2. The garbled SI flags (deCODE `READ_NAME_REGEX`, ApplyBQSR) need confirmation from deCODE.
3. The Sanger interval list is unpublished, which blocks exactness for the Sanger profile.
4. The reference identity is uncertain: deCODE vs Sanger FASTAs vs the GATK bundle `.alt` file (contig order, MD5s, chrY PAR).
5. Duplicate MC tags from bwa 0.7.17 + samblaster, and how htsjdk handles them. An unlisted fixmate step at deCODE ("Fix mate pair information") would rewrite the mate fields and change this answer.
6. The SI's "original qualities" CRAM statement conflicts with Vanguard Table S1.
7. JDK 8 math intrinsics may break x86 = arm64 equivalence in Java mode.
8. Sorter tie order (samtools 1.9, sambamba 0.5.9) changes HC downsampling input.
9. bwa 0.7.17 on arm64 is a modified build.
10. SRA read names leave the optical-duplicate path untested.
11. The per-chromosome HC tail makes golden runs long. Spot or preemptible interruptions would cost reruns.

## Effort

Engineer-days. These exclude the GATK internals, which the other worker covers.

| Item | Low | High |
|---|---:|---:|
| bwa 0.7.17 mode (no MQ, no `@HD`, `@PG`) + 0.7.17 oracle + tier0 regate | 0.5 | 1.5 |
| Raw per-RG FASTQ align entry (no fastp chunking) | 1 | 2 |
| Alt-aware mode (`.alt` parse, `AH:*`, live paths, oracle gating on alt-heavy reads) | 4 | 10 |
| samblaster 0.1.24 `-a --addMateTags` | 1 | 3 |
| Picard MarkDuplicates 2.20.3 / 2.18.26, queryname mode, regex/optical, per-RG vs per-sample | 8 | 20 |
| samtools sort 1.9 tie order (deCODE) | 1 | 3 |
| samtools merge 1.9 (header injection, ID suffixes) | 1 | 3 |
| sambamba 0.5.9 sort (Sanger only) | 2 | 5 |
| CRAM: gate via samtools oracle (no codec) | 0.5 | 1 |
| CRAM 3.0 writer for drop-in output (deferred) | 15 | 40 |
| Scatter: per-chromosome (deCODE) | 0.5 | 1 |
| GraphTyper 2.7.5 | out of scope | |
| **Pipeline side, deCODE profile** (no CRAM writer) | **18** | **45** |
| BQSR/ApplyBQSR/HC/GenotypeGVCFs 4.0.12.0 deltas | excluded (other worker) | |

Golden generation, separately:

| Item | Low | High |
|---|---:|---:|
| Drivers, reference/`.alt` fetch, containers, oracle builds (x86 and arm64 bwa/samtools/samblaster) | 2 | 4 |
| Slice r1/r2 + variability write-up + gates in GATES.md | 1 | 2 |
| Full 30× r1/r2 + GKL rerun + arm64 cross-check | 1 | 3 (mostly waiting) |
| Sanger profile goldens (after the interval list) | 1.5 | 3 |
| **Total** | **5.5** | **12** |
