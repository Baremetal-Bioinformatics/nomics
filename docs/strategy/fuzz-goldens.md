# Fuzz goldens: exactness on the inputs real samples rarely hit

The sample-matrix goldens show that nomics matches sarek 3.10.0 byte for byte on real exomes and
genomes. Real samples, though, exercise the common paths: most reads are clean, most variants are
simple SNPs, and coordinates stay in range. The bugs that break exactness live at the edges, as the
Picard short-wrap bug found on HG006 showed. The fuzz goldens aim synthetic reads at those edges and
gate nomics against sarek on each one, with the same checks as a real exome.

## How a fuzz golden runs

Each row of `reference/matrix/fuzz.tsv` names a scenario, a seed and a number of background pairs.
`scripts/goldens/exome_golden.sh` treats an `F*` id like an exome with three differences:

1. Its reads come from `scripts/fuzz/make_fuzz.py SCENARIO SEED PAIRS FASTA TARGETS OUTDIR`, run
   over the Spark's GATK GRCh38 FASTA and the Agilent v7 targets. The generator is deterministic:
   the same row gives the same bytes on any machine (seeded `random.Random`, gzip with mtime 0), so
   a golden can be rebuilt from its row alone. The md5s of the reads and the BED go in the artifact.
2. Its targets are the windows the reads came from (`targets.bed`, sorted in FASTA order and
   merged), so sarek's interval scatter and HaplotypeCaller cover every edge case.
3. It needs about 40 GB free instead of 220.

Everything else is the exome flow: sarek r1 and r2, r1 compared with r2, then a full nomics run
gated at every stage boundary (fastp chunks, mapped chunk BAMs, MarkDuplicates BAM and metrics,
recal table, recal BAM, VCF) against r1. A golden passes only when r1 equals r2 and every gate
passes. The background pairs (ordinary 2x151 reads at about 30x over target windows) keep all of
fastp's chunks non-empty and give BQSR a normal-looking table to learn from.

**Every window carries planted variants (2026-11-10).** The first fuzz goldens drew their
background and most scenario reads straight from the reference, so in nine of the fourteen
(F03-F10, F14) sarek's VCF had no records and the VCF gate compared nothing; F13 had one record.
The other stage gates were real. The generator now makes the sample diploid everywhere: each
background window gets one to three het or hom variants (SNVs, and deletions and insertions of
1-8 bases), and each scenario's own windows get SNVs, so edge-case reads sit next to calls.
Scenarios that put one read in a window (duplicates, chimera, unmapped, coord_wrap) add ordinary
30x coverage to some of those windows. A golden whose sarek VCF has fewer than 100 records
(`NOMICS_MIN_VCF_RECORDS`) now fails.

**Result, 2026-11-11.** All fourteen pass on the Spark: sarek r1 equals r2 at every judged
boundary, and nomics matches r1 at 18 of 18 gates.

| id | scenario | VCF records | | id | scenario | VCF records |
|---|---|---:|---|---|---|---:|
| F01 | homopolymer | 497 | | F08 | contig_edges | 546 |
| F02 | long_indel | 524 | | F09 | quality | 515 |
| F03 | short_insert | 525 | | F10 | read_length | 561 |
| F04 | n_bases | 507 | | F11 | pileup | 435 |
| F05 | duplicates | 606 | | F12 | complex_var | 690 |
| F06 | chimera | 474 | | F13 | overlap | 588 |
| F07 | unmapped | 685 | | F14 | coord_wrap | 631 |

F14 needed a second attempt: its first sarek run was stopped in a QC task when Seqera Wave
refused a container request (HTTP 429, a pull rate limit) after 27 sarek runs in under three
hours. The goldens still resolve containers through Wave on every task; `nomics verify` pins
image digests instead.

## The scenarios

| id  | scenario       | what it aims at |
|-----|----------------|-----------------|
| F01 | homopolymer    | indel errors and true indels in homopolymer runs (BWA gap placement, HC left-alignment) |
| F02 | long_indel     | het and hom deletions of 15-120 bp, insertions of 10-60 bp (soft clips, HC assembly) |
| F03 | short_insert   | fragments shorter than the reads, dovetails, inserts under fastp's 15 bp minimum |
| F04 | n_bases        | N runs at read ends, inside reads, at fastp's 5-N limit and over it, all-N reads |
| F05 | duplicates     | PCR and optical duplicates around Picard's 100-pixel distance, ties on base quality |
| F06 | chimera        | mates on different contigs, split reads (SA tags), inversions, 10 kb and 1 Mb inserts |
| F07 | unmapped       | random and low-complexity mates, both mates unmapped, adapter dimers |
| F08 | contig_edges   | chrM across its origin, contig starts and ends, N-gap edges, alt, decoy and EBV contigs |
| F09 | quality        | Q0-Q41 extremes, fastp's 40% low-quality limit, poly-G tails of 5-100 bp |
| F10 | read_length    | lengths 1-151, R1 and R2 of different lengths, around fastp's 15 bp minimum |
| F11 | pileup         | three loci at 3,000-8,000x with het SNPs and indels (HC downsampling, BQSR counts) |
| F12 | complex_var    | SNP clusters, MNPs, two alt alleles at one site, indels next to SNPs, STRs |
| F13 | overlap        | fragment equal to read length, start-tie pairs, mates that fully overlap |
| F14 | coord_wrap     | tile, x and y above 32767 and 65535, which Picard keeps as Java shorts |

Most scenarios use NovaSeq read names (`A00744:...`), which turn on fastp's poly-G trimming as on
the real NovaSeq samples. `unmapped` and `read_length` use HiSeq names, so the other fastp path is
covered too.

## What a pass means, and what it doesn't

A pass says nomics gives sarek's bytes on that scenario's reads at every stage. It does not say the
calls are right: synthetic reads have no truth set worth scoring, and accuracy is sarek's question,
not nomics'. A scenario that fails is a bug report with a reproducer: the row, the generator and the
gate log are all it takes to rebuild the case.

New scenarios go in as new functions in `make_fuzz.py` and new rows in `fuzz.tsv`. Changing an
existing scenario changes its reads, so the script notices that the saved sarek runs were made from
other inputs and runs sarek again.
