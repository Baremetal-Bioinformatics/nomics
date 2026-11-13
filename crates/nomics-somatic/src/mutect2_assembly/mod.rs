//! Mutect2 local assembly, default sarek tumour-normal mode. Experimental until
//! the proposed stage-two gates pass; see the oracle verification record.
pub mod pruner;
use crate::mutect2_regions::{self, Inputs, Read, Stats};
use nomics_hc::{
    activity::Region,
    assembly::{self, AsmRead},
    finalize::finalize_region_with_overlap,
    trim::{self, TrimHap, TrimResult},
};
use nomics_io::{bam::CIGAR_OPS, bam_query::BamQuery, fasta::read_fai};
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufWriter, Read as _, Seek, SeekFrom, Write},
    path::Path,
};

#[derive(Debug, Default)]
pub struct AssemblyStats {
    pub assembled: usize,
    pub non_reference: usize,
    pub before: usize,
    pub after: usize,
    pub artificial: usize,
}
fn cigar(c: &[u32]) -> String {
    c.iter()
        .map(|&c| format!("{}{}", c >> 4, CIGAR_OPS[(c & 15) as usize] as char))
        .collect()
}
fn span(contig: &str, start: i64, end: i64) -> String {
    format!("{contig}:{start}-{end}")
}
fn dump(w: &mut impl Write, marker: &str, haps: &[assembly::Haplotype]) -> io::Result<()> {
    writeln!(w, "{marker}S {}", haps.len())?;
    for h in haps {
        writeln!(
            w,
            "{marker} {} {} {} {} {} {} {}",
            h.is_ref,
            h.alignment_start_hap_wrt_ref,
            cigar(&h.cigar),
            span(&h.contig, h.start, h.end),
            crate::math::javafmt::double_to_string(h.score),
            h.kmer_size,
            String::from_utf8_lossy(&h.bases)
        )?;
    }
    Ok(())
}
fn remove_duplicates(reads: &[&Read]) -> Vec<usize> {
    let mut groups: HashMap<(&str, i64), Vec<usize>> = HashMap::new();
    for (i, r) in reads.iter().enumerate() {
        let a = &r.alignment;
        if a.flag & 1 != 0
            && a.flag & 8 == 0
            && (a.mate_contig != a.contig || i64::from(a.tlen).abs() > 1_000_000)
        {
            let unclipped = a.start
                - a.cigar
                    .iter()
                    .take_while(|c| matches!(**c & 15, 4 | 5))
                    .map(|c| i64::from(c >> 4))
                    .sum::<i64>();
            groups
                .entry((
                    &r.sample,
                    if a.flag & 64 != 0 {
                        unclipped
                    } else {
                        -unclipped
                    },
                ))
                .or_default()
                .push(i);
        }
    }
    let mut keep = vec![true; reads.len()];
    for group in groups.values() {
        let mut mates: HashMap<i32, Vec<usize>> = HashMap::new();
        for &i in group {
            mates
                .entry(reads[i].alignment.mate_contig)
                .or_default()
                .push(i);
        }
        for m in mates.values() {
            for &i in m.iter().skip(usize::from(m.len() > group.len() / 2)) {
                keep[i] = false;
            }
        }
    }
    (0..reads.len()).filter(|&i| keep[i]).collect()
}
/// Run stages one and two. `output` is the ordered internal dump; `artificial`
/// is the stage-determined projection of artificial BAM records in file order.
pub fn run(
    inputs: &Inputs<'_>,
    output: &Path,
    artificial: &Path,
) -> io::Result<(Stats, AssemblyStats)> {
    let mut out = BufWriter::new(File::create(output)?);
    let (traversal, stats, records) = run_with_sink(inputs, &mut out, true, &mut |_| Ok(()))?;
    let mut out = BufWriter::new(File::create(artificial)?);
    for (_, _, _, record) in records {
        out.write_all(record.as_bytes())?;
    }
    out.flush()?;
    Ok((traversal, stats))
}

/// Borrowed ordered stage-two boundary; original reads carry overlap-edited BI/BD.
pub struct AssembledRegion<'a> {
    pub ordinal: usize,
    pub reads: &'a [&'a Read],
    pub finalized: &'a [nomics_hc::finalize::HcRead],
    pub haplotypes: &'a [assembly::Haplotype],
    pub padded_span: (i64, i64),
    pub callable_span: (i64, i64),
    pub padded_reference: &'a [u8],
    pub reference_start: i64,
}
type ArtificialRecord = (usize, i64, String, String);
/// Consume each varying trimmed region immediately, preserving shared-read timing.
pub fn run_with_sink(
    inputs: &Inputs<'_>,
    out: &mut impl Write,
    artificial: bool,
    sink: &mut impl FnMut(AssembledRegion<'_>) -> io::Result<()>,
) -> io::Result<(Stats, AssemblyStats, Vec<ArtificialRecord>)> {
    run_with_sinks(inputs, out, artificial, sink, &mut mutect2_regions::NoLoci)
}
/// `run_with_sink`, with the traversal's per-locus pileups also delivered to `loci`
/// (the f1r2 collector), in the same single pass over the reads.
pub fn run_with_sinks(
    inputs: &Inputs<'_>,
    out: &mut impl Write,
    artificial: bool,
    sink: &mut impl FnMut(AssembledRegion<'_>) -> io::Result<()>,
    loci: &mut dyn mutect2_regions::LocusSink,
) -> io::Result<(Stats, AssemblyStats, Vec<ArtificialRecord>)> {
    if inputs.input.len() != 2 {
        return Err(io::Error::other(
            "only two BAM tumour-normal mode is supported",
        ));
    }
    let header = BamQuery::open(&inputs.input[0], &[])?.header;
    let fai = read_fai(Path::new(&format!("{}.fai", inputs.reference.display())))?;
    let mut reference = File::open(inputs.reference)?;
    let mut stats = AssemblyStats::default();
    let mut ordinal = 0;
    let mut records = Vec::new();
    let mut regions = |r: &Region, reads: &[&Read]| {
        ordinal += 1;
        if !r.active || reads.is_empty() {
            return Ok(());
        }
        stats.assembled += 1;
        let contig = String::from_utf8_lossy(&header.refs[r.contig].name);
        let length = i64::from(header.refs[r.contig].len);
        let (ps, pe) = ((r.start - 100).max(1), (r.end + 100).min(length));
        let (rs, re) = ((ps - 500).max(1), (pe + 500).min(length));
        let e = fai
            .iter()
            .find(|e| e.name == contig)
            .ok_or_else(|| io::Error::other("missing FASTA contig"))?;
        let mut bases = Vec::with_capacity((re - rs + 1) as usize);
        let mut pos = (rs - 1) as u64;
        while pos < re as u64 {
            let n = (e.line_bases - pos % e.line_bases).min(re as u64 - pos) as usize;
            reference.seek(SeekFrom::Start(
                e.offset + pos / e.line_bases * e.line_width + pos % e.line_bases,
            ))?;
            let at = bases.len();
            bases.resize(at + n, 0);
            reference.read_exact(&mut bases[at..])?;
            pos += n as u64;
        }
        bases.make_ascii_uppercase();
        let raw: Vec<_> = remove_duplicates(reads)
            .into_iter()
            .map(|i| {
                let mut a = reads[i].alignment.clone();
                a.record = i as u32;
                a
            })
            .collect();
        let finalized = finalize_region_with_overlap(length, r.start, r.end, &raw, false);
        writeln!(
            out,
            "R {ordinal} {} {} {}",
            span(&contig, r.start, r.end),
            span(&contig, ps, pe),
            finalized.len()
        )?;
        writeln!(
            out,
            "REF {} {}",
            span(&contig, rs, re),
            String::from_utf8_lossy(&bases)
        )?;
        for a in &finalized {
            let original = reads[a.record as usize];
            let left_hard = |c: &[u32]| {
                c.iter()
                    .take_while(|c| **c & 15 == 5)
                    .map(|c| (c >> 4) as usize)
                    .sum::<usize>()
            };
            let offset = left_hard(&a.cigar) - left_hard(&original.alignment.cigar);
            let hex = |bytes: &[u8]| bytes.iter().map(|q| format!("{q:02x}")).collect::<String>();
            writeln!(
                out,
                "D {} {} {} {} {} {} {} {} {}",
                original.sample,
                String::from_utf8_lossy(&a.name),
                a.flag,
                a.start,
                cigar(&a.cigar),
                String::from_utf8_lossy(&a.bases),
                hex(&a.quals),
                hex(&original.insertion_quals[offset..offset + a.bases.len()]),
                hex(&original.deletion_quals[offset..offset + a.bases.len()])
            )?;
        }
        let assembly_reads: Vec<_> = finalized
            .iter()
            .map(|a| AsmRead {
                name: a.name.clone(),
                flag: a.flag,
                start: a.start,
                cigar: a.cigar.clone(),
                bases: a.bases.clone(),
                quals: a.quals.clone(),
            })
            .collect();
        let samples: Vec<_> = finalized
            .iter()
            .map(|a| usize::from(reads[a.record as usize].is_normal))
            .collect();
        let haps = assembly::assemble_with_pruner(
            &contig,
            r.start,
            r.end,
            rs,
            &bases,
            &assembly_reads,
            Some(&pruner::chains_to_remove),
            Some(&samples),
        );
        stats.before += haps.len();
        stats.non_reference += usize::from(haps.iter().any(|h| !h.is_ref));
        dump(out, "H", &haps)?;
        let trim_haps: Vec<_> = haps
            .iter()
            .map(|h| TrimHap {
                bases: h.bases.clone(),
                is_ref: h.is_ref,
                alignment_start_hap_wrt_ref: i64::from(h.alignment_start_hap_wrt_ref),
                cigar: h.cigar.clone(),
                contig: h.contig.clone(),
                start: h.start,
                end: h.end,
            })
            .collect();
        match trim::trim_region_with_mnp(
            &contig,
            r.start,
            r.end,
            length,
            rs,
            &bases,
            &trim_haps,
            &[],
            1,
        ) {
            TrimResult::NoVariation => writeln!(out, "T none")?,
            TrimResult::Trimmed {
                span: s,
                padded_span: p,
                variation_present,
                haplotypes,
                ..
            } => {
                writeln!(
                    out,
                    "T {} {} {variation_present}",
                    span(&contig, s.0, s.1),
                    span(&contig, p.0, p.1)
                )?;
                let trimmed: Vec<_> = haplotypes
                    .into_iter()
                    .map(|h| {
                        let original = trim_haps
                            .iter()
                            .enumerate()
                            .find(|(_, o)| {
                                o.is_ref == h.is_ref
                                    && trim::trim_haplotype(o, p.0, p.1)
                                        .is_some_and(|t| t.bases == h.bases)
                            })
                            .map(|(i, _)| i)
                            .expect("trimmed original");
                        assembly::Haplotype {
                            bases: h.bases,
                            is_ref: h.is_ref,
                            alignment_start_hap_wrt_ref: h.alignment_start_hap_wrt_ref as i32,
                            cigar: h.cigar,
                            contig: h.contig,
                            start: h.start,
                            end: h.end,
                            score: haps[original].score,
                            kmer_size: if h.is_ref {
                                0
                            } else {
                                haps[original].kmer_size
                            },
                        }
                    })
                    .collect();
                stats.after += trimmed.len();
                dump(out, "h", &trimmed)?;
                if variation_present {
                    sink(AssembledRegion {
                        ordinal,
                        reads,
                        finalized: &finalized,
                        haplotypes: &trimmed,
                        padded_span: p,
                        callable_span: s,
                        padded_reference: &bases,
                        reference_start: rs,
                    })?;
                }
                if variation_present && artificial {
                    for (i, h) in trimmed.iter().enumerate() {
                        let counter = records.len() + 1;
                        let name = format!("HC_{counter}_{i}");
                        let hash = h.bases.iter().fold(1i32, |v, &b| {
                            v.wrapping_mul(31).wrapping_add(i32::from(b as i8))
                        });
                        records.push((r.contig,h.start,name,format!("{}\t{}\t{}\t{}\tHC:i:{hash}\tRG:Z:ArtificialHaplotypeRG\tCR:Z:{}\n",contig,h.start,cigar(&h.cigar),String::from_utf8_lossy(&h.bases),span(&contig,s.0,s.1))));
                    }
                }
            }
        }
        Ok(())
    };
    let traversal = mutect2_regions::run_with_sinks(inputs, &mut regions, loci)?;
    out.flush()?;
    records.sort_by(|a, b| (&a.0, &a.1, &a.2).cmp(&(&b.0, &b.1, &b.2)));
    stats.artificial = records.len();
    Ok((traversal, stats, records))
}
