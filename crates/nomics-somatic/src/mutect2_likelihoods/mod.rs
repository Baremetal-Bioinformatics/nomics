//! Mutect2's Java PairHMM and read realignment boundary, before callMutations.
//! The standalone projection invokes somatic genotyping to include XA tags.
//! Likelihood matrices retain the pre-genotyping boundary.
mod header;
pub mod records;
use crate::{
    mutect2_assembly::{self, AssembledRegion},
    mutect2_regions::{Inputs, Stats},
};
use nomics_hc::{
    finalize::{HcRead, coordinate_cmp, hard_clip_to_region},
    likelihood,
    pairhmm_trie::{self, HmmRead},
    readprep::hmm_qualities_with_indels,
    trim::TrimHap,
};
use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

#[derive(Default, Debug)]
pub struct LikelihoodStats {
    pub regions: usize,
    pub stubs: usize,
    pub evaluated: usize,
    pub mapq_capped_reads: usize,
    pub base_quality_floored_reads: usize,
    pub pcr_capped_reads: usize,
    pub indel_quality_floored_reads: usize,
    pub capped_reads: usize,
    pub capped_values: usize,
    pub poorly_modelled: usize,
    pub realigned: usize,
}
pub(crate) struct Aligned {
    pub source: (usize, u64),
    pub read: HcRead,
    pub ins: Vec<u8>,
    pub del: Vec<u8>,
    pub hc: Option<i32>,
    pub xa: Option<String>,
    pub indel_tags: bool,
    pub alignment_region: String,
    pub callable_region: String,
}
pub(crate) fn cigar(c: &[u32]) -> String {
    c.iter()
        .map(|c| {
            format!(
                "{}{}",
                c >> 4,
                nomics_io::bam::CIGAR_OPS[(c & 15) as usize] as char
            )
        })
        .collect()
}
fn hex(q: &[u8]) -> String {
    q.iter().map(|q| format!("{q:02x}")).collect()
}
fn left_hard(r: &HcRead) -> usize {
    r.cigar
        .iter()
        .take_while(|c| **c & 15 == 5)
        .map(|c| (c >> 4) as usize)
        .sum()
}
fn matrix(
    w: &mut impl Write,
    mark: &str,
    sample: &str,
    reads: &[HcRead],
    v: &[Vec<f64>],
) -> io::Result<()> {
    writeln!(w, "{mark} {sample} {}", reads.len())?;
    for r in reads {
        writeln!(w, "E {}", String::from_utf8_lossy(&r.name))?;
    }
    for (h, v) in v.iter().enumerate() {
        write!(w, "L {h}")?;
        for x in v {
            write!(w, " {:016x}", x.to_bits())?;
        }
        writeln!(w)?;
    }
    Ok(())
}
/// Execute stages one through three; outputs are diagnostic boundary projections,
/// not BAM/VCF. Existing stage-one and stage-two commands retain their paths.
pub fn run(
    inputs: &Inputs<'_>,
    output: &Path,
    realigned_records: &Path,
) -> io::Result<(Stats, LikelihoodStats)> {
    run_with_sink(inputs, output, realigned_records, &mut |region| {
        crate::mutect2_calls::supported_alleles(region, inputs.normal_sample)
    })
}

/// Read-haplotype natural-log likelihoods after realignment, before genotyping.
pub struct SampleLikelihoods {
    pub sample: String,
    pub reads: Vec<HcRead>,
    pub values: Vec<Vec<f64>>,
}
pub struct LikelihoodRegion<'a> {
    pub(crate) xa: &'a std::cell::RefCell<std::collections::BTreeMap<u32, String>>,
    pub haplotypes: &'a [TrimHap],
    pub samples: &'a [SampleLikelihoods],
    pub padded_span: (i64, i64),
    pub callable_span: (i64, i64),
    pub padded_reference: &'a [u8],
    pub reference_start: i64,
}
/// A synchronous consumer preserves traversal's shared-read mutation timing.
pub fn run_with_sink(
    inputs: &Inputs<'_>,
    output: &Path,
    realigned_records: &Path,
    sink: &mut impl FnMut(LikelihoodRegion<'_>) -> io::Result<()>,
) -> io::Result<(Stats, LikelihoodStats)> {
    run_with_sinks(
        inputs,
        output,
        realigned_records,
        sink,
        &mut crate::mutect2_regions::NoLoci,
    )
}
/// `run_with_sink`, with the traversal's per-locus pileups also delivered to `loci`.
pub fn run_with_sinks(
    inputs: &Inputs<'_>,
    output: &Path,
    realigned_records: &Path,
    sink: &mut impl FnMut(LikelihoodRegion<'_>) -> io::Result<()>,
    loci: &mut dyn crate::mutect2_regions::LocusSink,
) -> io::Result<(Stats, LikelihoodStats)> {
    let mut samples = Vec::new();
    for path in inputs.input {
        let header = nomics_io::bam_query::BamQuery::open(path, &[])?.header;
        for line in header
            .text
            .split(|&b| b == b'\n')
            .filter(|l| l.starts_with(b"@RG\t"))
        {
            if line.split(|&b| b == b'\t').any(|f| f.starts_with(b"FO:")) {
                return Err(io::Error::other("flow-order read groups are not supported"));
            }
            if let Some(sm) = line
                .split(|&b| b == b'\t')
                .find_map(|f| f.strip_prefix(b"SM:"))
            {
                samples.push(String::from_utf8_lossy(sm).into_owned());
            }
        }
    }
    samples.sort();
    samples.dedup();
    let mut w = BufWriter::new(File::create(output)?);
    let mut stats = LikelihoodStats::default();
    let mut aligned = Vec::new();
    let (traversal, _, _) = mutect2_assembly::run_with_sinks(
        inputs,
        &mut io::sink(),
        false,
        &mut |region| process(region, &samples, &mut w, &mut stats, &mut aligned, sink),
        loci,
    )?;
    w.flush()?;
    records::write(inputs, aligned, realigned_records)?;
    Ok((traversal, stats))
}
fn process(
    region: AssembledRegion<'_>,
    samples: &[String],
    w: &mut impl Write,
    stats: &mut LikelihoodStats,
    aligned: &mut Vec<Aligned>,
    sink: &mut impl FnMut(LikelihoodRegion<'_>) -> io::Result<()>,
) -> io::Result<()> {
    let (start, end) = region.padded_span;
    let mut reads: Vec<_> = region
        .finalized
        .iter()
        .cloned()
        .map(|r| hard_clip_to_region(r, start, end))
        .filter(|r| !r.bases.is_empty() && r.start <= end && r.end() >= start)
        .collect();
    reads.sort_by(coordinate_cmp);
    let stubs = reads.iter().filter(|r| r.bases.len() < 10).count();
    reads.retain(|r| r.bases.len() >= 10);
    stats.regions += 1;
    stats.stubs += stubs;
    stats.evaluated += reads.len();
    let haps: Vec<_> = region
        .haplotypes
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
    let contig = &haps[0].contig;
    writeln!(w, "R {} {contig}:{start}-{end} {stubs}", region.ordinal)?;
    writeln!(w, "HS {}", haps.len())?;
    for h in region.haplotypes {
        writeln!(
            w,
            "H {} {} {} {}:{}-{} {} {} {}",
            h.is_ref,
            h.alignment_start_hap_wrt_ref,
            cigar(&h.cigar),
            h.contig,
            h.start,
            h.end,
            crate::math::javafmt::double_to_string(h.score),
            h.kmer_size,
            String::from_utf8_lossy(&h.bases)
        )?;
    }
    let mut batches = Vec::new();
    for sample in samples {
        let rs: Vec<_> = reads
            .iter()
            .filter(|r| region.reads[r.record as usize].sample == *sample)
            .cloned()
            .collect();
        let indels: Vec<_> = rs
            .iter()
            .map(|r| {
                let original = region.reads[r.record as usize];
                let offset = left_hard(r) - left_hard(&original.alignment);
                (
                    original.insertion_quals[offset..offset + r.bases.len()].to_vec(),
                    original.deletion_quals[offset..offset + r.bases.len()].to_vec(),
                )
            })
            .collect();
        for (r, (ins, del)) in rs.iter().zip(&indels) {
            stats.mapq_capped_reads += usize::from(r.quals.iter().any(|&q| q > r.mapq));
            stats.base_quality_floored_reads +=
                usize::from(r.quals.iter().any(|&q| (q.min(r.mapq) as i8) < 18));
            stats.indel_quality_floored_reads +=
                usize::from(ins.iter().chain(del).any(|&q| (q as i8) < 6));
            let cache = nomics_hc::readprep::pcr_cache();
            stats.pcr_capped_reads += usize::from((0..r.bases.len().saturating_sub(1)).any(|i| {
                let cap = cache[nomics_hc::readprep::tandem_repeat_length(&r.bases, i)] as i8;
                ins[i] as i8 > cap || del[i] as i8 > cap
            }));
        }
        let prepared: Vec<_> = rs
            .iter()
            .zip(&indels)
            .map(|(r, (i, d))| hmm_qualities_with_indels(&r.bases, &r.quals, r.mapq, i, d))
            .collect();
        writeln!(w, "S {sample} {}", rs.len())?;
        for (r, (q, i, d)) in rs.iter().zip(&prepared) {
            writeln!(
                w,
                "D {} {} {} {} {} {} {} {} {}",
                String::from_utf8_lossy(&r.name),
                r.flag,
                r.start,
                cigar(&r.cigar),
                String::from_utf8_lossy(&r.bases),
                hex(q),
                hex(i),
                hex(d),
                "0a".repeat(r.bases.len())
            )?;
        }
        let gcps: Vec<_> = rs.iter().map(|r| vec![10; r.bases.len()]).collect();
        let hmm_reads: Vec<_> = rs
            .iter()
            .zip(&prepared)
            .zip(&gcps)
            .map(|((r, (q, i, d)), g)| HmmRead {
                bases: &r.bases,
                quals: q,
                ins: i,
                del: d,
                gcp: g,
            })
            .collect();
        let hb: Vec<_> = haps.iter().map(|h| h.bases.as_slice()).collect();
        // This shared prefix kernel is portable Rust (no nomics-simd calls).
        let raw = pairhmm_trie::likelihoods(&hb, &hmm_reads);
        batches.push((sample.to_owned(), rs, indels, raw));
    }
    for (s, rs, _, raw) in &batches {
        matrix(w, "RAW", s, rs, raw)?;
    }
    // Keep intermediate capping before HC's shared finish performs the same cap
    // and removes poorly modelled reads, so every adjustment has its own gate.
    for (s, rs, _, raw) in &mut batches {
        for r in 0..rs.len() {
            let best = raw.iter().map(|v| v[r]).fold(f64::NEG_INFINITY, f64::max);
            let cap = best - 4.5;
            let mut capped = false;
            for v in raw.iter_mut() {
                if v[r] < cap {
                    v[r] = cap;
                    stats.capped_values += 1;
                    capped = true;
                }
            }
            stats.capped_reads += usize::from(capped);
        }
        matrix(w, "CAP", s, rs, raw)?;
    }
    let mut kept = Vec::new();
    for (s, rs, indels, raw) in batches {
        let lens: Vec<_> = rs.iter().map(|r| r.bases.len()).collect();
        let result = likelihood::finish(raw, &lens);
        stats.poorly_modelled += result.disqualified.len();
        let reads: Vec<_> = result.kept.iter().map(|&i| rs[i].clone()).collect();
        let indels: Vec<_> = result.kept.iter().map(|&i| indels[i].clone()).collect();
        matrix(w, "KEEP", &s, &reads, &result.values)?;
        kept.push((s, reads, indels, result.values));
    }
    for (s, rs, _, v) in &mut kept {
        for row in v.iter_mut() {
            for x in row {
                *x *= nomics_hc::jmath::log(10.0);
            }
        }
        matrix(w, "LN", s, rs, v)?;
    }
    let aligned_start = aligned.len();
    let mut sample_likelihoods = Vec::new();
    for (s, rs, indels, v) in kept {
        let realigned = nomics_hc::realign::realign_reads_natural_log(&rs, &haps, start, &v);
        let mut reads = Vec::new();
        for ((r, hc), (ins, del)) in realigned.into_iter().zip(indels) {
            writeln!(
                w,
                "A {s} {} {} {} {} {}",
                String::from_utf8_lossy(&r.name),
                r.flag,
                r.start,
                cigar(&r.cigar),
                hc.map_or_else(|| "null".into(), |v| v.to_string())
            )?;
            stats.realigned += 1;
            reads.push(r.clone());
            let indel_tags = region.reads[r.record as usize].overlap_indels;
            aligned.push(Aligned {
                source: region.reads[r.record as usize].source,
                read: r,
                ins,
                del,
                hc,
                xa: None,
                indel_tags,
                alignment_region: format!("{contig}:{start}-{end}"),
                callable_region: format!(
                    "{contig}:{}-{}",
                    region.callable_span.0, region.callable_span.1
                ),
            });
        }
        sample_likelihoods.push(SampleLikelihoods {
            sample: s,
            reads,
            values: v,
        });
    }
    let xa = std::cell::RefCell::new(std::collections::BTreeMap::new());
    sink(LikelihoodRegion {
        xa: &xa,
        haplotypes: &haps,
        samples: &sample_likelihoods,
        padded_span: region.padded_span,
        callable_span: region.callable_span,
        padded_reference: region.padded_reference,
        reference_start: region.reference_start,
    })?;
    let xa = xa.into_inner();
    for a in &mut aligned[aligned_start..] {
        a.xa = xa.get(&a.read.record).cloned();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn bytes(s: &str) -> Vec<u8> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
            .collect()
    }
    #[test]
    fn independently_generated_gatk_quality_and_likelihood_bits() {
        for line in include_str!("numerics.tsv").lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            let b = f[0].as_bytes();
            let (q, i, d) = hmm_qualities_with_indels(
                b,
                &bytes(f[2]),
                f[1].parse().unwrap(),
                &bytes(f[3]),
                &bytes(f[4]),
            );
            assert_eq!(q, bytes(f[5]), "BQ {}", f[0]);
            assert_eq!(i, bytes(f[6]), "BI {}", f[0]);
            assert_eq!(d, bytes(f[7]), "BD {}", f[0]);
            let hmm = HmmRead {
                bases: b,
                quals: &q,
                ins: &i,
                del: &d,
                gcp: &vec![10; b.len()],
            };
            let v = pairhmm_trie::likelihoods(&[f[8].as_bytes()], &[hmm])[0][0];
            assert_eq!(v.to_bits(), u64::from_str_radix(f[9], 16).unwrap());
            assert_eq!(
                (v * nomics_hc::jmath::log(10.0)).to_bits(),
                u64::from_str_radix(f[10], 16).unwrap()
            );
        }
    }
}
