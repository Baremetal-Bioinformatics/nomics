# Third-party notices

nomics links no bioinformatics or compression library and contains no code copied from the
tools below. It does reproduce their behavior, read from their source, and its tests use a
few of their files as oracles. This file lists both.

## Tools whose behavior nomics reproduces

| Tool | Version | License | nomics crate |
|---|---|---|---|
| [BWA](https://github.com/lh3/bwa) | 0.7.19-r1273 | GPL-3.0 | `nomics-align` |
| [samtools / htslib](https://github.com/samtools) | 1.22.1 (sort), 1.24 (merge), 1.21 (CRAM writing) | MIT/Expat | `nomics-align`, `nomics-io` |
| [fastp](https://github.com/OpenGene/fastp) | 1.1.0 | MIT | `nomics-prep` |
| [GATK](https://github.com/broadinstitute/gatk) | 4.6.2.0 | BSD-3-Clause | `nomics-markdup`, `nomics-bqsr`, `nomics-hc`, `nomics-joint` |
| [Picard](https://github.com/broadinstitute/picard) | 3.4.0 (as bundled in GATK) | MIT | `nomics-markdup` |
| [htsjdk](https://github.com/samtools/htsjdk) | 4.2.0 | MIT | `nomics-io` |
| [GenomicsDB](https://github.com/GenomicsDB/GenomicsDB) | as bundled in GATK 4.6.2.0 | MIT | `nomics-joint` |
| [bcftools](https://github.com/samtools/bcftools) | as run by sarek 3.10.0 | MIT/Expat | `nomics-joint` |
| [zlib](https://zlib.net) | deflate output reproduced exactly | zlib | `nomics-io` |
| [Salmon](https://github.com/COMBINE-lab/salmon) (experimental) | 2.7.0 | GPL-3.0 | `nomics-quant` |
| [Subread featureCounts](https://subread.sourceforge.net) (experimental) | 2.1.1 | GPL-3.0 | `nomics-count` |
| [STAR / STARsolo](https://github.com/alexdobin/STAR) (experimental) | 2.7.10b | MIT | `nomics-solo` |
| [nf-core/sarek](https://github.com/nf-core/sarek), [nf-core/rnaseq](https://github.com/nf-core/rnaseq) | 3.10.0, 3.27.0 | MIT | pipeline structure, `scripts/nfcore/` |

nomics-align follows BWA's algorithms as written, so it is distributed under BWA's license,
GPL-3.0; the whole project is GPL-3.0-only.

## Files copied into this repository

- `crates/nomics-align/tests/c/`: unmodified files from bwa 0.7.19. `ksort.h`, `ksw.c` and
  `kbtree.h` carry their own license headers (klib, MIT; `kbtree.h` BSD-2-Clause). `ksw.h`
  (klib), `malloc_wrap.h`, `neon_sse.h` and `scalar_sse.h` have no header and are covered by
  bwa's own license, GPL-3.0. Test oracles only; never linked into `nomics`.
- `crates/nomics-count/src/biotypes_header.txt`: unmodified from nf-core/rnaseq 3.27.0
  (MIT); see `crates/nomics-count/THIRD_PARTY.md`.
- `reference/oracle/hc/*.java`: our own test programs that call GATK 4.6.2.0 classes to dump
  internal state. They are run against a GATK jar you supply and are not part of `nomics`.

## Data

Test fixtures under `crates/*/tests/fixtures/` are outputs we generated from the tools above
(for example GATK's quality tables, Java formatting and random-number streams, glibc
`log`/`exp` samples) or from public HG002 data (GIAB; SRA SRR14724482 and SRR14724532).
No reference-pipeline output is redistributed beyond small derived fixtures and digest
manifests (`reference/*.digest.tsv`).

## Rust dependencies

General-purpose crates from crates.io, under their own licenses: clap, rayon, libc, serde,
serde_json, statrs and mimalloc (MIT or MIT/Apache-2.0). `cargo tree` lists the full set.
