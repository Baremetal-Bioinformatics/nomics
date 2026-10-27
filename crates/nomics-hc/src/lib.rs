//! Germline small-variant calling, equivalent to GATK 4.6.2.0
//! `HaplotypeCaller` (single sample, one run over the whole exome BED) on
//! aarch64, where PairHMM and Smith-Waterman run GATK's Java implementations.
//! Spec: `docs/spec/haplotypecaller.md`.
//!
//! Gate: `variant_calling/haplotypecaller/*.haplotypecaller.vcf.gz`, records
//! byte-identical, header identical except `##GATKCommandLine`.

pub mod activity;
pub mod annotate;
pub mod assembly;
pub mod call;
pub mod caller;
pub mod dbsnp;
pub mod finalize;
pub mod genotype;
pub mod gvcf;
pub mod javafmt;
pub mod jmath;
pub mod jrandom;
pub mod likelihood;
pub mod pairhmm;
pub mod pairhmm_trie;
pub mod readprep;
pub mod reads;
pub mod realign;
pub mod refconfidence;
pub mod tables;
pub mod trim;
pub mod vcf;
