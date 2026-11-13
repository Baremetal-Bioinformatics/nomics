//! The Mutect2 tumour-normal chain of sarek 3.10.0 (GATK 4.6.2.0), one module per tool.
//! Specs: docs/spec/mutect2.md (the caller) and docs/spec/mutect2-filtering.md (the
//! orientation model, pileup summaries, contamination and FilterMutectCalls).

pub mod contamination;
pub mod get_pileup_summaries;
mod math;
pub mod pileup_summary;
mod segmentation;

mod java_hashmap;
pub mod orientation;

pub mod filtering;
pub mod merge;

pub mod mutect2_regions;

pub mod mutect2_assembly;

pub mod mutect2_likelihoods;

pub mod mutect2_calls;

pub mod mutect2_f1r2;
