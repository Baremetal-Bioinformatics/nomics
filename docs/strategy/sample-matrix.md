# Sample matrix: exactness across samples

Which samples the exactness gates cover beyond HG002, where the reads come from and what each
adds. Every row was run through sarek 3.10.0 twice on the DGX Spark (r1 and r2) and through
nomics, gated at every stage boundary; the results are in
[RESULTS.md](../../RESULTS.md#exactness-across-samples). The drivers are
`scripts/goldens/exome_golden.sh` and `scripts/goldens/genome_golden.sh`, reading
`reference/matrix/*.tsv`.

**Sources.** ENA served the Spark at about 0.15 MB/s (`scripts/goldens/probe_sources.sh`,
2026-11-06; Google's bucket gave 40 MB/s, S3 20–30 MB/s), so E1–E4 come from Google's copies of
the same Baid et al. runs (same individual, kit and instrument) and E5 from nf-core's S3 mirror.
Google's copies keep the instrument read names, so poly-G trimming is exercised on E1 and E2 as
well as E6, and optical-duplicate detection on every Google exome. HG002 and E5 keep SRA-style
names, so that path stays covered.

**Genome depths.** A whole-genome sarek run needs about 5x its FASTQ size in scratch even with
results published as hard links, and the Spark had 266 GB free, so the genomes are 10x and 20x
copies from Google's bucket rather than the 34x and 60x rows first planned. G2 is HG001 PCR-free
(a female genome, native names: poly-G and optical duplicates on). G3 is HG003 with a PCR
library. Each genome reduces each sarek run to a digest manifest as soon as it finishes, and
nomics is gated against r1 and against r2.

## The matrix

Six exomes and three genomes. Every accession, size and md5 below was read from the source
listing on 2026-10-30 (ENA filereport API, GCS JSON API, S3 listings); sizes are bytes of
R1 / R2 `.fastq.gz`.

| # | Sample | Run / file | Source | Instrument, reads | Kit or prep | Pairs | Size (R1 / R2) | md5 (R1 / R2) | New axes |
|---|---|---|---|---|---|---:|---|---|---|
| E1 | HG001 | `HG001.novaseq.wes_agilent.50x.R{1,2}.fastq.gz` | `gs://brain-genomics-public/research/sequencing/fastq/novaseq/wes_agilent/50x/` | NovaSeq 6000, 2×151, **native read names** | Agilent SureSelect v7 | ~20 M (est.) | 2,873,556,366 / 2,898,394,171 | 1ef717eebc6bbed2b2c160596cd7e060 / 0a3ba208464eaca4fffd54f5ac1146fc | individual, poly-G and optical dups on |
| E2 | HG003 | `HG003.novaseq.wes_idt.100x.R{1,2}.fastq.gz` | `gs://brain-genomics-public/research/sequencing/fastq/novaseq/wes_idt/100x/` | NovaSeq 6000, 2×151, **native read names** | IDT xGen | ~20 M (est.) | 2,610,112,136 / 2,681,020,186 | 08b57492faae928f2660a66eac72c0ea / 1419932c6371af5f6a1086e885346d0e | individual, kit, poly-G and optical dups on |
| E3 | HG004 | `HG004.hiseq4000.wes_idt.50x.R{1,2}.fastq.gz` | `gs://brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_idt/50x/` | HiSeq 4000, 2×150, **native read names** | IDT xGen | ~15 M (est.) | 2,013,085,701 / 2,577,847,567 | 9e55947f665c6718ed4faa279116a8a2 / 3ea5fd972ba8b0fb4a86e03c767419d9 | individual, instrument, read length, optical dups on |
| E4 | HG005 | `HG005.hiseq4000.wes_agilent.50x.R{1,2}.fastq.gz` | `gs://brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/` | HiSeq 4000, 2×150, **native read names** | Agilent SureSelect v7 | ~25 M (est.) | 3,377,752,690 / 4,017,940,282 | 25ee18a0a61eeab421be2968fc3a7963 / 077e73edd8d309a8f83462d00e979eb0 | individual, optical dups on |
| E5 | HCC1395BL | SRR7890919 | `s3://ngi-igenomes/test-data/sarek/` (nf-core's mirror of ENA SRP162370; same sizes and md5s) | HiSeq 1500, 2×126 | Agilent SureSelect XT v6 + UTR | 88,498,907 | 6,630,757,329 / 6,788,404,175 | 092cbd4c3a126509a3a6d06f14dca748 / 9a5a26a33a62cbcfbaa170150c5417b5 | instrument, read length, kit, depth (2× the HG002 exome) |
| E6 | HG006 | `HG006.novaseq.wes_idt.100x.R{1,2}.fastq.gz` | `s3://genomics-benchmark-datasets/google-brain/fastq/novaseq/wes_idt/100x/` (copy of `gs://brain-genomics-public/research/sequencing/fastq/`) | NovaSeq 6000, 2×151, **native read names** | IDT xGen | ~20 M (est.) | 3,099,073,237 / 3,282,613,139 | 642cb502836be6e726b0916c77bd23dc / 100645485ae4a4de63b1c142bc275a33 | individual, **fastp poly-G and optical dups on** |
| G1 | HG002 | `HG002.novaseqX.10x.R{1,2}.fastq.gz` (subset of SRR37356338) | `gs://brain-genomics-public/research/sequencing/fastq/novaseqx/` | NovaSeq X 25B, 2×150 | WGS | ~100 M (est.) | 7,624,745,528 / 7,565,366,293 | c1369eea32bceb4876dcd307e262e4fb / 44c72edb94351370db83ae7071cef850 | instrument, low depth (10×) |
| G2 | HG001 | `HG001.novaseq.pcr-free.20x.R{1,2}.fastq.gz` (Baid et al., PRJNA734598) | `gs://brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/20x/` | NovaSeq 6000, 2×151, **native read names** | WGS, PCR-free | ~230 M (est.) | 17,420,820,941 / 18,067,877,906 | aea0647f25a2aa86f5593aa303cbe728 / 32384c23dab61a296f6c1b86b4d43f00 | female genome, **poly-G and optical dups on in a genome** |
| G3 | HG003 | `HG003.novaseq.pcr-plus.20x.R{1,2}.fastq.gz` (Baid et al., PRJNA734598) | `gs://brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_plus/20x/` | NovaSeq 6000, 2×151, native read names | WGS, **PCR** prep | ~230 M (est.) | 17,562,264,097 / 18,171,557,418 | 4ffc529a85ce200ea44556bfe87e541a / 06a8a800119b48707d679d59f0379436 | individual, duplicate-heavy PCR library |

Spare, if a row fails for a reason outside nomics: HG007 Agilent v7 exome with native names
(`google-brain/fastq/novaseq/wes_agilent/100x/HG007…`, md5 5a8d5436dde7e7dfd7485afcfdb03233 /
2fca705605776367f661a9e26eb1f6c3). G2 fallback: NovaSeq X runs SRR37356339 + SRR37356338 as two
lanes of one sample (~85×, NIST, public domain), which loses the poly-G axis.

### Axis coverage

| Axis | Today | After |
|---|---|---|
| Individual | HG002 | + HG001, HG003, HG004, HG005, HG006, HCC1395BL; a female genome (G2) |
| Exome kit | Agilent v7 | + IDT xGen (E2, E3, E6), Agilent v6 + UTR (E5) |
| Instrument | NovaSeq 6000 | + NovaSeq X (G1), HiSeq 4000 (E3, E4), HiSeq 1500 (E5) |
| Read length | 2×151 | + 2×150 (E3, E4, G1), 2×126 (E5) |
| Genome depth | 30× | + 10× (G1), 20× (G2, G3) |
| Library | PCR-free genome | + PCR genome (G3) |
| fastp poly-G | never on for real input | on for E1, E2, E6 (NovaSeq), G2 and G3 |
| Optical duplicates | never on for real input | on for E1–E4, E6, G2 and G3 |

## Two things the brief assumed that the data says otherwise

- **NovaSeq X from SRA does not switch poly-G on.** fastp enables poly-G trimming from the first
  R1 name's prefix (`isTwoColorSystem`; our port in `crates/nomics-prep/src/lib.rs`). Every
  SRA/ENA FASTQ, including the NovaSeq X runs, is renamed to `@SRR…`, so the prefix never
  matches. Our HG002 inputs are SRA-named too (`@SRR14724482.1 1/1`). The same names also
  keep MarkDuplicates' optical-duplicate detection off, since it needs tile and x/y fields
  (`crates/nomics-markdup/src/optical.rs`). Real instrument names survive only in Google's
  copies of the Baid et al. data (`@A00744:…`), so E6 and G2 turn on both poly-G trimming and
  optical-duplicate detection on real reads.
- **The SEQC2 normal is not 2×126 alone.** It is also the heaviest exome here (88.5 M pairs,
  2.2× SRR14724482) and a lymphoblastoid line from a cancer study. It is a fine germline sample;
  it is not a GIAB reference sample, so it adds no hap.py record.

## Target BEDs

| Kit | BED | Source |
|---|---|---|
| Agilent v7 | `Agilent_v7.bed` (already on the Spark) | `s3://ngi-igenomes/test-data/sarek/` |
| IDT xGen | `idt_capture.grch38.bed`, 4,721,168 bytes, md5 fd44b2ee8dca744dd734a42045da3baf, chr-prefixed | `gs://brain-genomics-public/research/sequencing/grch38/bed/` |
| Agilent v6 + UTR | `S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed`, 64,015,319 bytes | `s3://ngi-igenomes/test-data/sarek/` (sarek's own `test_full`) |

Genomes use the bundle's WGS calling regions, as today.

## Terms

| Source | Terms |
|---|---|
| PRJNA734598 (E1–E4, G3) and the Google copies (E6, G1, G2) | Baid et al. 2020 data, "released under a CC-0 license" (bucket README; AWS Open Data registry entry lists CC0 1.0). The `novaseqx/` folder postdates that README and is not in the AWS copy, so its CC0 status is inferred, not stated. G1 is a subset of NIST's SRR37356338 (PRJNA1427896), which is US-government data. |
| SRP162370 (E5) | Public SRA data from the SEQC2 consortium (Fang et al. 2021, Nat Biotech); no access controls. Redistribution of results is normal practice; no explicit licence found. |
| BEDs | Agilent BEDs as redistributed by nf-core; IDT BED as redistributed by Google. Neither carries a stated licence. We publish digests and counts, not the BEDs. |

We publish digests, counts and verdicts, never read data, so every row is publishable on these
terms.
