//! Paired-end short-read alignment, equivalent to `bwa mem -K 100000000 -Y`
//! (BWA 0.7.19-r1273).
//!
//! The FM-index lives in a `nomics_mem::HostBuf` so CPU and GPU share it in
//! place. Determinism contract: insert-size statistics are estimated per
//! 100 Mbp input batch (`-K`), and batches restart at each fastp chunk.
//!
//! Gate: `preprocessing/mapped/*.sorted.bam` records.

pub mod chain;
pub mod dense_sa;
pub mod diag;
pub mod extend;
pub mod fused;
pub mod glibc_math;
pub mod index;
pub mod kbtree;
pub mod kmer_table;
pub mod ksort;
pub mod ksw;
pub mod ksw_lanes;
pub mod libm;
pub mod mem;
pub mod pair;
pub mod sam;
pub mod samsort;
pub mod seed;
pub mod slice;
