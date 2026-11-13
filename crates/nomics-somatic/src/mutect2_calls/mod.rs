//! Stage-four somatic evidence, annotations, phasing, and raw VCF output.
//! The diagnostic retains its original evidence and annotation boundaries.
//! F1R2 archive collection is not wired here.
mod annotation_dump;
pub mod annotations;
pub mod evidence;
pub mod fragments;
pub mod phasing;
pub mod vcf;

use crate::{
    mutect2_likelihoods::{self, LikelihoodRegion},
    mutect2_regions::{Inputs, Stats},
};
use std::{
    collections::BTreeSet,
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

#[derive(Default, Debug)]
pub struct GenotypingStats {
    pub regions: usize,
    pub events: usize,
    /// Events passing default tumour and normal thresholds, before annotation.
    pub passing_events: usize,
}

fn bits(w: &mut impl Write, name: &str, values: &[f64]) -> io::Result<()> {
    write!(w, "{name}")?;
    for v in values {
        write!(w, " {:016x}", v.to_bits())?;
    }
    writeln!(w)
}

fn matrix(w: &mut impl Write, name: &str, m: &[Vec<f64>]) -> io::Result<()> {
    writeln!(w, "{name} {} {}", m.len(), m[0].len())?;
    for row in m {
        bits(w, "L", row)?;
    }
    Ok(())
}

struct AnnotationWriters {
    raw: Vec<u8>,
    formatted: Vec<u8>,
    phasing: Vec<u8>,
}
fn dump(
    region: LikelihoodRegion<'_>,
    normal: &str,
    w: &mut impl Write,
    stats: &mut GenotypingStats,
    annotations: &mut AnnotationWriters,
    mut vcf: Option<&mut vcf::Writer>,
) -> io::Result<()> {
    stats.regions += 1;
    let haps = region.haplotypes;
    writeln!(
        w,
        "R {}:{}-{}",
        haps[0].contig, region.callable_span.0, region.callable_span.1
    )?;
    let mut samples = Vec::new();
    for sample in region.samples {
        let mut values = sample.values.clone();
        fragments::cap(&mut values);
        writeln!(w, "S {}", sample.sample)?;
        matrix(w, "C", &values)?;
        let (fragments, values) = fragments::group(&sample.reads, &values);
        for f in &fragments {
            writeln!(
                w,
                "F {} {} {} {}",
                f.start,
                f.end,
                f.evidence
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                f.reads
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
        }
        matrix(w, "G", &values)?;
        samples.push((sample.sample.as_str(), fragments, values));
    }
    let events: Vec<_> = haps
        .iter()
        .map(|h| {
            nomics_hc::trim::events_for_haplotype_with_mnp(
                h,
                region.padded_reference,
                region.reference_start,
                1,
            )
        })
        .collect();
    let starts: BTreeSet<_> = events
        .iter()
        .flatten()
        .map(|e| e.start)
        .filter(|&s| region.callable_span.0 <= s && s <= region.callable_span.1)
        .collect();
    let mut calls = Vec::new();
    let mut extra = Vec::new();
    let mut somatic = Vec::new();
    for start in starts {
        let (alleles, mapper) = nomics_hc::call::event_allele_mapping(
            &events,
            start,
            region.padded_reference[(start - region.reference_start) as usize],
            false,
        );
        if alleles.is_empty() {
            continue;
        }
        stats.events += 1;
        writeln!(
            w,
            "V {start} {}",
            alleles
                .iter()
                .map(|a| String::from_utf8_lossy(a))
                .collect::<Vec<_>>()
                .join(",")
        )?;
        let mut tumor = Vec::new();
        let mut normal_matrix = Vec::new();
        for (sample, fragments, values) in &samples {
            let end = start + alleles[0].len() as i64 - 1;
            let keep: Vec<_> = fragments
                .iter()
                .enumerate()
                .filter(|(_, f)| f.start <= end + 2 && f.end >= start - 2)
                .map(|(i, _)| i)
                .collect();
            let m: Vec<Vec<_>> = mapper
                .iter()
                .map(|hs| {
                    keep.iter()
                        .map(|&f| {
                            hs.iter()
                                .map(|&h| values[h][f])
                                .fold(f64::NEG_INFINITY, f64::max)
                        })
                        .collect()
                })
                .collect();
            writeln!(
                w,
                "M {sample} {}",
                keep.iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
            matrix(w, "A", &m)?;
            if *sample == normal {
                normal_matrix = m;
            } else {
                tumor = m;
            }
        }
        let tlod = evidence::somatic_log_odds(&tumor);
        let nlod = evidence::diploid_alt_log_odds(&normal_matrix);
        let nalod = evidence::somatic_log_odds(&normal_matrix);
        bits(w, "TLOD", &tlod)?;
        bits(w, "NLOD", &nlod)?;
        bits(w, "NALOD", &nalod)?;
        for a in 1..alleles.len() {
            if alleles[a] != b"*"
                && tlod[a - 1] > 3.0 * nomics_hc::jmath::log(10.0)
                && nlod[a - 1] > 2.2 * nomics_hc::jmath::log(10.0)
            {
                let e = vcf::event(start, &alleles[0], &alleles[a]);
                if !somatic.contains(&e) {
                    somatic.push(e);
                }
            }
        }
        let passing = (1..alleles.len()).any(|a| {
            alleles[a] != b"*"
                && tlod[a - 1] > 3.0 * nomics_hc::jmath::log(10.0)
                && nlod[a - 1] > 2.2 * nomics_hc::jmath::log(10.0)
        });
        writeln!(w, "PASS {passing}")?;
        stats.passing_events += usize::from(passing);
        if passing {
            calls.push(annotation_dump::write(
                &region,
                annotation_dump::EventInput {
                    start,
                    alleles: &alleles,
                    mapper: &mapper,
                    tlod: &tlod,
                    samples: &samples,
                },
                &mut annotations.raw,
                &mut annotations.formatted,
            )?);
            if let Some(writer) = vcf.as_deref_mut() {
                extra.push(writer.extra(&region, start, &alleles, &tlod, &nlod, &nalod)?);
            }
        }
    }
    phasing::write(
        &calls,
        &events,
        &haps[0].contig,
        &region
            .samples
            .iter()
            .map(|s| s.sample.clone())
            .collect::<Vec<_>>(),
        normal,
        &mut annotations.phasing,
    )?;
    if let Some(writer) = vcf {
        vcf::event_counts(&region, &events, &calls, &somatic, &mut extra);
        writer.region(&annotations.formatted, &annotations.phasing, extra)?;
        annotations.raw.clear();
        annotations.formatted.clear();
        annotations.phasing.clear();
    }
    Ok(())
}

/// Ordered diagnostic only: all fragment matrices, per-event marginal matrices
/// and unrounded tumour/normal LODs. Appended `.annotations.tsv` and
/// `.annotation-fields.tsv` sidecars contain an explicit subset of annotations
/// as raw attributes and VCF field text respectively; neither is a VCF.
/// The other files retain stage-three meaning.
pub fn run(
    inputs: &Inputs<'_>,
    output: &Path,
    likelihoods: &Path,
    realigned_records: &Path,
) -> io::Result<(Stats, GenotypingStats)> {
    let mut w = BufWriter::new(File::create(output)?);
    let mut stats = GenotypingStats::default();
    let sidecar = |suffix: &str| {
        let mut path = output.as_os_str().to_os_string();
        path.push(suffix);
        File::create(path)
    };
    let mut annotations = AnnotationWriters {
        raw: Vec::new(),
        formatted: Vec::new(),
        phasing: Vec::new(),
    };
    let (traversal, _) = mutect2_likelihoods::run_with_sink(
        inputs,
        likelihoods,
        realigned_records,
        &mut |region| {
            dump(
                region,
                inputs.normal_sample,
                &mut w,
                &mut stats,
                &mut annotations,
                None,
            )
        },
    )?;
    w.flush()?;
    sidecar(".annotations.tsv")?.write_all(&annotations.raw)?;
    sidecar(".annotation-fields.tsv")?.write_all(&annotations.formatted)?;
    sidecar(".phasing.tsv")?.write_all(&annotations.phasing)?;
    Ok((traversal, stats))
}

pub(crate) fn supported_alleles(region: LikelihoodRegion<'_>, normal: &str) -> io::Result<()> {
    dump(
        region,
        normal,
        &mut io::sink(),
        &mut GenotypingStats::default(),
        &mut AnnotationWriters {
            raw: Vec::new(),
            formatted: Vec::new(),
            phasing: Vec::new(),
        },
        None,
    )
}
