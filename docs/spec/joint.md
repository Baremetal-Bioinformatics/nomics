# Joint genotyping (GATK 4.6.2.0 GenotypeGVCFs)

Scope: sarek 3.10.0 `--joint_germline` on **one** diploid sample, i.e. the chain
`GATK4_GENOMICSDBIMPORT` (`--merge-input-intervals`, one BED group) →
`GATK4_GENOTYPEGVCFS` (`--dbsnp`, no `--intervals`, default
`StandardAnnotation` group, VQSR skipped with `haplotypecaller_filter`) →
`BCFTOOLS_SORT` (bcftools/htslib 1.23.1, `--output-type z`) →
`MERGE_GENOTYPEGVCFS` (GATK `MergeVcfs`) → `joint_germline.vcf.gz`.
`nomics joint-genotype --gvcf <g.vcf.gz> --dbsnp <dbsnp.vcf.gz> -o joint_germline.vcf.gz`
reproduces that file from the HaplotypeCaller gVCF (`crates/nomics-joint`).

## Source contract

Pinned source: [GATK 4.6.2.0](https://github.com/broadinstitute/gatk/tree/4.6.2.0),
paths under `src/main/java/org/broadinstitute/hellbender/`;
[htslib 1.23.1](https://github.com/samtools/htslib/tree/1.23.1);
sarek `modules/nf-core/gatk4/genotypegvcfs`, `conf/modules/joint_germline.config`.

- **Traversal** (`tools/walkers/GenotypeGVCFs.java`, `engine/VariantLocusWalker.java`):
  without `--include-non-variant-sites` or force-output intervals the walker
  runs by variant, one record per `apply`. Reference blocks merge to a
  reference-only allele list and are never variants, so only gVCF records with
  an alt besides `<NON_REF>` matter.
- **GenomicsDB order.** The driving iterator walks the workspace arrays
  (`<contig>$<start>$<end>`, one per contig under `--merge-input-intervals`) in
  sorted name order, so contigs are visited `chr1, chr10 … chr19, chr2, chr20 …`.
  `bcftools sort` restores dictionary order in the output, but the visit order
  fixes which QD-jitter draw each record gets (below). Evidence: under
  dictionary order chr1 (603 draws) matched and the first chr2 record matched
  draw 3,627 of `Random(47382911)`. That number is exactly the chr1 + chr10–19
  draw count, and the whole file matches in sorted array order.
- **GenomicsDB drops the gVCF ID column.** The output ID comes from the dbSNP
  overlap alone. Carrying the gVCF IDs through made 1,228 more records
  differ, each by its ID column only.
- **Merge** (`tools/walkers/ReferenceConfidenceVariantContextMerger.java`, `merge`):
  target alleles are the record's alleles without `<NON_REF>` (single sample:
  no remapping, no spanning events). PLs are re-indexed through
  `newToOldGenotypeMap` with no renormalization, and AD is subset by index.
  INFO values are parsed by header type (`Integer.parseInt` /
  `Double.parseDouble`). Flags (`DB`) and multi-valued keys are dropped.
  `RAW_MQandDP` is a reducible raw key (single sample: unchanged). `DP` becomes
  the record's INFO DP as a string. `AC, AF, AN, MLEAC, MLEAF, END` are removed.
- **Regenotype** (`GenotypeGVCFsEngine.regenotypeVC`, `GenotypingEngine.calculateGenotypes`):
  only records with DP > 0. Same AF calculator as HaplotypeCaller in VCF mode
  (`nomics_hc::genotype::GenotypingEngine::new`: standard confidence 30,
  `EMIT_VARIANTS_ONLY`, upstream-deletion bookkeeping for `*`). `PREFER_PLS`
  equals `USE_PLS_TO_ASSIGN` when PLs exist. The phasing flag of GT is kept.
  The engine returning nothing covers `isProperlyPolymorphic` and
  spanning-deletion-only sites.
- **Finalize** (`RMSMappingQuality.finalizeRawData`): `MQ = %.2f` of
  `sqrt(sumSq / n)` from `RAW_MQandDP`. The raw key is not kept.
  `reverseTrimAlleles` then removes shared trailing bases (`*` and symbolic
  alleles do not take part, and each allele keeps one base).
- **Annotate** (`VariantAnnotatorEngine.annotateContext`, no read likelihoods):
  - `ChromosomeCounts`: AC/AF/AN over called alleles.
  - `ExcessHet`: `"0.0000"`. With one diploid sample the exact test's p-value is 1.
  - `FisherStrand`/`StrandOddsRatio` from the genotype's SB via
    `getTableFromSamples`, with minimum counts 2 and 0 (`nomics_hc::annotate::{fisher_strand, strand_odds_ratio}`).
  - `QualByDepth`: `getDepth` over het/hom-var genotypes, with the AD sum and the
    AD-restricted depth when alt reads > 1. QD ≥ 35 is replaced by
    `30 + 3·nextGaussian()` from the tool-wide `Random(47382911)`.
  - Rank sums and Coverage need reads, so the merged values stay.
    `InbreedingCoeff` needs ≥ 10 samples.
- **dbSNP** (`utils/variant/VariantOverlapAnnotator`): `getRsID` uses dbSNP
  records starting at POS whose split, trimmed events match the final alleles.
  IDs are joined with `;`, and `DB` is set when any match
  (`nomics_hc::dbsnp::DbsnpSlice::id_for`).
- **Genotype cleanup** (`cleanupGenotypeAnnotations`, `keepSB = false`): SB is
  dropped, MIN_DP moves to DP, and PGT becomes `1|1` for hom-var.
- **Text** (`htsfmt.rs`): GATK formats values (htsjdk `formatVCFDouble` for
  Double attributes, QUAL via `formatQualValue`). `bcftools sort` then parses
  each Float INFO field with `hts_str2dbl`, stores it as a 32-bit float and
  prints it with `kputd` (6 significant digits, trailing zeros culled; for
  example `3.20` becomes `3.2` and `0.500` becomes `0.5`). QUAL goes through
  `atof` → float → `kputd`. `MergeVcfs` keeps INFO strings verbatim but
  re-prints QUAL with `formatQualValue` (`12345.7` becomes `12345.70`).
  INFO keys are sorted. FORMAT is GT followed by the sorted keys.
- **Header**: the gVCF header minus `GVCFBlock*`, HaplotypeCaller's
  `GATKCommandLine` and `source`, plus GenotypeGVCFs' lines (AC, AF, AN, DB,
  DP, FS, MQ, QD, SOR, RGQ, LowQual), bcftools' `FILTER=PASS`, the
  GenomicsDBImport command line and `source=GenomicsDBImport`, in htsjdk
  `VCFHeader` order.

## Verified on

Tier1 HG002 WES, x86 golden pair (`reference/run_tier1_joint_cloud_pair.sh`,
r1 == r2): **42,362 / 42,362 records identical**, and the header is identical
apart from `##GATKCommandLine`. The result is the same from the golden gVCF and
from nomics' own HaplotypeCaller gVCF (`hcgate2-20261015`), so FASTQ → joint
VCF is exact. The arm64 gVCF golden equals the x86 gVCF golden record for
record. Exercised in that data: 1,869 records with dropped alleles, 1,043
reverse-trimmed, 146 with `*` alleles, 2,036 with multiple dbSNP IDs, 6,241
QD jitter draws, and 3,535 gVCF variant records correctly not emitted.

**Not exercised** (ported from source, unverified by any golden):
no-call output genotypes, `excludeFromAnnotations` (GQ0 hom-ref with no DP,
which is not implemented), records without SB, more than six alternate
alleles, multiple samples, `--include-non-variant-sites`. A gVCF variant record
without PL is an error.

## Not produced

`joint_germline.vcf.gz.tbi`. The gate compares records and header, not BGZF
bytes. MergeVcfs' compression is not reproduced.
