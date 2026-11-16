//! Explicit *subset* of the post-annotation boundary. This is not a VCF and
//! does not stand in for the full internal gate. Every supported field is also
//! written in VCF form to make seven-contig annotation checks Java-independent.
use super::{
    annotations::{self, SampleInput},
    fragments::Fragment,
};
use crate::mutect2_likelihoods::LikelihoodRegion;
use std::io::{self, Write};

type GroupedSample<'a> = (&'a str, Vec<Fragment>, Vec<Vec<f64>>);
pub(super) struct EventInput<'a> {
    pub start: i64,
    pub alleles: &'a [Vec<u8>],
    pub mapper: &'a [Vec<usize>],
    pub tlod: &'a [f64],
    pub samples: &'a [GroupedSample<'a>],
}

/// trimAlleles(true,true): reverse first, then forward; restore one base
/// when normalization would empty an allele. No symbolic alleles are emitted
/// by this default tumour-normal diagnostic.
pub(super) fn trim(mut start: i64, mut alleles: Vec<Vec<u8>>) -> (i64, Vec<Vec<u8>>) {
    if alleles.iter().any(|a| a.len() == 1) {
        return (start, alleles);
    }
    let mut end = 0;
    while alleles
        .iter()
        .all(|a| a.len() > end && a[a.len() - 1 - end] == alleles[0][alleles[0].len() - 1 - end])
    {
        end += 1;
    }
    let mut front = 0;
    while alleles
        .iter()
        .all(|a| a.len() > end + front && a[front] == alleles[0][front])
    {
        front += 1;
    }
    if alleles.iter().any(|a| a.len() == end + front) {
        if front == 0 {
            end -= 1;
        } else {
            front -= 1;
        }
    }
    for a in &mut alleles {
        *a = a[front..a.len() - end].to_vec();
    }
    start += front as i64;
    (start, alleles)
}
fn ints(xs: &[i32]) -> String {
    format!(
        "[{}]",
        xs.iter()
            .map(|x| format!("Integer:{x}"))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn doubles(xs: &[f64]) -> String {
    format!(
        "[{}]",
        xs.iter()
            .map(|x| format!("d:{:016x}", x.to_bits()))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn csv(xs: &[i32]) -> String {
    xs.iter().map(i32::to_string).collect::<Vec<_>>().join(",")
}

pub(super) fn write(
    region: &LikelihoodRegion<'_>,
    event: EventInput<'_>,
    raw: &mut impl Write,
    formatted: &mut impl Write,
) -> io::Result<super::phasing::Call> {
    let EventInput {
        start,
        alleles,
        mapper,
        tlod,
        samples,
    } = event;
    let emit: Vec<_> = std::iter::once(0)
        .chain(
            (1..alleles.len())
                .filter(|&a| alleles[a] != b"*" && tlod[a - 1] > 3.0 * nomics_hc::jmath::log(10.0)),
        )
        .collect();
    let emitted: Vec<_> = emit.iter().map(|&a| alleles[a].clone()).collect();
    let (pos, emitted) = trim(start, emitted);
    let order = annotations::allele_map_order(&emitted)?;
    let end = start + alleles[0].len() as i64 - 1;
    let mut rs = Vec::new();
    let mut rms = Vec::new();
    let mut fs = Vec::new();
    let mut fms = Vec::new();
    for (sn, (_, fragments, values)) in samples.iter().enumerate() {
        let sample = &region.samples[sn];
        let mut capped = sample.values.clone();
        super::fragments::cap(&mut capped);
        let read_keep: Vec<_> = sample
            .reads
            .iter()
            .enumerate()
            .filter(|(_, r)| r.start <= end + 2 && r.end() >= start - 2)
            .map(|(i, _)| i)
            .collect();
        let fragment_keep: Vec<_> = fragments
            .iter()
            .enumerate()
            .filter(|(_, f)| f.start <= end + 2 && f.end >= start - 2)
            .map(|(i, _)| i)
            .collect();
        let marginal = |matrix: &[Vec<f64>], keep: &[usize]| {
            emit.iter()
                .map(|&a| {
                    keep.iter()
                        .map(|&r| {
                            mapper[a]
                                .iter()
                                .map(|&h| matrix[h][r])
                                .fold(f64::NEG_INFINITY, f64::max)
                        })
                        .collect()
                })
                .collect::<Vec<Vec<f64>>>()
        };
        rs.push(
            read_keep
                .iter()
                .map(|&i| sample.reads[i].clone())
                .collect::<Vec<_>>(),
        );
        rms.push(marginal(&capped, &read_keep));
        fs.push(
            fragment_keep
                .iter()
                .map(|&i| &fragments[i])
                .collect::<Vec<_>>(),
        );
        fms.push(marginal(values, &fragment_keep));
    }
    let inputs: Vec<_> = (0..samples.len())
        .map(|sn| SampleInput {
            reads: &rs[sn],
            fragment_reads: &region.samples[sn].reads,
            read_matrix: &rms[sn],
            fragments: &fs[sn],
            fragment_matrix: &fms[sn],
        })
        .collect();
    for sn in 0..samples.len() {
        for (column, fragment) in fs[sn].iter().enumerate() {
            let (allele, confidence) = annotations::best(&fms[sn], column, &order);
            if confidence > 0.2 {
                for &read in &fragment.reads {
                    let record = region.samples[sn].reads[read].record;
                    let mut tags = region.xa.borrow_mut();
                    let value = tags.entry(record).or_default();
                    if !value.is_empty() {
                        value.push_str(", ");
                    }
                    value.push_str(&format!("{}:{pos}={allele}", region.haplotypes[0].contig));
                }
            }
        }
    }
    let (info, genotypes) = annotations::annotate(pos, &order, &inputs);
    let bases = emitted
        .iter()
        .map(|a| String::from_utf8_lossy(a))
        .collect::<Vec<_>>()
        .join(",");
    let chrom = &region.haplotypes[0].contig;
    writeln!(raw, "V {chrom}:{pos} {bases}")?;
    let strand = info
        .strands
        .iter()
        .map(|s| csv(s))
        .collect::<Vec<_>>()
        .join("|");
    writeln!(raw, "I AS_SB_TABLE String:{strand}")?;
    writeln!(raw, "I DP String:{}", info.dp)?;
    writeln!(raw, "I MBQ {}", ints(&info.mbq))?;
    writeln!(raw, "I MFRL {}", ints(&info.mfrl))?;
    writeln!(raw, "I MMQ {}", ints(&info.mmq))?;
    writeln!(raw, "I MPOS {}", ints(&info.mpos))?;
    writeln!(formatted, "V {chrom}:{pos} {bases}")?;
    writeln!(formatted, "I AS_SB_TABLE {strand}")?;
    writeln!(formatted, "I DP {}", info.dp)?;
    for (k, v) in [
        ("MBQ", &info.mbq),
        ("MFRL", &info.mfrl),
        ("MMQ", &info.mmq),
        ("MPOS", &info.mpos),
    ] {
        writeln!(formatted, "I {k} {}", csv(v))?;
    }
    for ((name, _, _), g) in samples.iter().zip(genotypes) {
        writeln!(raw, "S {name}")?;
        writeln!(raw, "AD {}", ints(&g.ad))?;
        writeln!(raw, "DP Integer:{}", g.dp)?;
        writeln!(raw, "G AF {}", doubles(&g.af))?;
        for (key, values) in [
            ("F1R2", g.f1r2.as_slice()),
            ("F2R1", g.f2r1.as_slice()),
            ("FAD", g.fad.as_slice()),
            ("SB", g.sb.as_slice()),
        ] {
            writeln!(raw, "G {key} {}", ints(values))?;
        }
        writeln!(formatted, "S {name}")?;
        writeln!(formatted, "AD {}", csv(&g.ad))?;
        writeln!(formatted, "DP {}", g.dp)?;
        writeln!(
            formatted,
            "G AF {}",
            g.af.iter()
                .map(|&x| nomics_hc::vcf::format_vcf_double(x))
                .collect::<Vec<_>>()
                .join(",")
        )?;
        for (key, values) in [
            ("F1R2", g.f1r2.as_slice()),
            ("F2R1", g.f2r1.as_slice()),
            ("FAD", g.fad.as_slice()),
            ("SB", g.sb.as_slice()),
        ] {
            writeln!(formatted, "G {key} {}", csv(values))?;
        }
    }
    Ok(super::phasing::Call {
        start: pos,
        alleles: emitted,
        haplotypes: emit
            .iter()
            .flat_map(|&a| mapper[a].iter().copied())
            .collect(),
    })
}
