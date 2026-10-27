//! Bioinformatics file formats, written from the specs (SAMv1, VCFv4.2, BGZF,
//! CRAMv3).
//!
//! Byte-identical BAM/VCF output is a gate, so the encoders must reproduce the
//! reference writers exactly, including BGZF block boundaries, deflate level and
//! implementation, header text and tag order. BAMs are gated on records and the
//! uncompressed stream, so by default they keep htsjdk's block boundaries but
//! use a faster deflater (`bgzf::BamCompression`); htsjdk's exact bytes stay
//! available with `NOMICS_BAM_COMPRESSION=htsjdk`. Normalizations allowed by
//! `docs/GATES.md` (e.g. `@PG`, `@SQ UR:` paths) are applied by `nomics-gate`,
//! never by the writers.

pub mod bam;
pub mod bamstream;
pub mod bgzf;
pub mod cram;
pub mod crc32;
pub mod deflate;
pub mod fasta;
pub mod gzip;
pub mod htsjdk_index;
pub mod inflate;
pub mod spill;
pub mod tabix;
pub mod vcf;
