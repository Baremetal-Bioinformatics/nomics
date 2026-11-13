//! Raw Mutect2 VCF serialization. Schema is GATK 4.6.2.0's default header.
use super::{AnnotationWriters, GenotypingStats, dump};
use crate::{
    mutect2_likelihoods::{self, LikelihoodRegion},
    mutect2_regions::Inputs,
};
use nomics_hc::vcf::format_vcf_double;
use nomics_io::{bgzf::ParallelBgzfWriter, tabix::Tabix, vcf::parse_line};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, Write},
    path::Path,
};
type Fields = BTreeMap<String, String>;
pub(super) struct Writer {
    output: ParallelBgzfWriter<File>,
    pon: Tabix,
    population: Tabix,
    meta: Vec<String>,
    contigs: Vec<(String, u64)>,
    normal: String,
}
#[derive(Default)]
struct Record {
    locus: String,
    alleles: String,
    info: Fields,
    samples: BTreeMap<String, Fields>,
}
fn parse(data: &[u8]) -> Vec<Record> {
    let mut records: Vec<Record> = Vec::new();
    let mut sample = String::new();
    for line in std::str::from_utf8(data).unwrap().lines() {
        let f: Vec<_> = line.split(' ').collect();
        match f[0] {
            "V" => records.push(Record {
                locus: f[1].into(),
                alleles: f[2].into(),
                ..Record::default()
            }),
            "S" => {
                sample = f[1].into();
            }
            "I" => {
                records
                    .last_mut()
                    .unwrap()
                    .info
                    .insert(f[1].into(), f[2].into());
            }
            _ => {
                let (k, v) = if f[0] == "G" {
                    (f[1], f[2])
                } else {
                    (f[0], f[1])
                };
                records
                    .last_mut()
                    .unwrap()
                    .samples
                    .entry(sample.clone())
                    .or_default()
                    .insert(k.into(), v.into());
            }
        }
    }
    records
}
// Mutect2Engine.getAttributeAsDoubleList: only the missing-value marker
// becomes the default; malformed numeric text is an error, even if unmatched.
fn resource_af_values(values: &[u8]) -> io::Result<Vec<f64>> {
    values
        .split(|&b| b == b',')
        .map(|v| {
            if v == b"." {
                return Ok(1e-6);
            }
            std::str::from_utf8(v)
                .ok()
                .and_then(|s| s.trim().parse().ok())
                .ok_or_else(|| io::Error::other("invalid germline resource AF value"))
        })
        .collect()
}
impl Writer {
    pub(super) fn extra(
        &self,
        region: &LikelihoodRegion<'_>,
        start: i64,
        alleles: &[Vec<u8>],
        tlod: &[f64],
        nlod: &[f64],
        nalod: &[f64],
    ) -> io::Result<Fields> {
        let emit: Vec<_> = (1..alleles.len())
            .filter(|&a| alleles[a] != b"*" && tlod[a - 1] > 3.0 * nomics_hc::jmath::log(10.0))
            .collect();
        let chrom = region.haplotypes[0].contig.as_bytes();
        let mut fields = Fields::new();
        for (key, values) in [("TLOD", tlod), ("NLOD", nlod), ("NALOD", nalod)] {
            fields.insert(
                key.into(),
                emit.iter()
                    .map(|&a| format_vcf_double(values[a - 1] * std::f64::consts::LOG10_E))
                    .collect::<Vec<_>>()
                    .join(","),
            );
        }
        if !self.pon.lines_starting_at(chrom, start as u64)?.is_empty() {
            fields.insert("PON".into(), String::new());
        }
        let mut af = vec![1e-6; emit.len()];
        for line in self.population.lines_starting_at(chrom, start as u64)? {
            let rec = parse_line(&line)?;
            let Some(values) = rec
                .info
                .split(|&b| b == b';')
                .find_map(|s| s.strip_prefix(b"AF="))
            else {
                continue;
            };
            // determineReferenceAllele checks equal-length disagreement only;
            // createAlleleMapping extends shorter alleles using the longer REF.
            if alleles[0].len() == rec.ref_allele.len() && alleles[0] != rec.ref_allele {
                return Err(io::Error::other(
                    "inconsistent equal-length germline resource REF",
                ));
            }
            let alts: Vec<_> = rec.alt.split(|&b| b == b',').collect();
            let values = resource_af_values(values)?;
            if alts.len() != values.len() {
                continue;
            }
            for (i, &a) in emit.iter().enumerate() {
                for (j, alt) in alts.iter().enumerate() {
                    let matches = if alleles[0].len() >= rec.ref_allele.len() {
                        alleles[a] == [alt, &alleles[0][rec.ref_allele.len()..]].concat()
                    } else {
                        [alleles[a].as_slice(), &rec.ref_allele[alleles[0].len()..]].concat()
                            == *alt
                    };
                    if matches {
                        af[i] = values[j];
                    }
                }
            }
        }
        fields.insert(
            "POPAF".into(),
            af.iter()
                .map(|&x| format_vcf_double(-nomics_hc::jmath::log10(x)))
                .collect::<Vec<_>>()
                .join(","),
        );
        let emitted = std::iter::once(alleles[0].clone())
            .chain(emit.iter().map(|&a| alleles[a].clone()))
            .collect();
        let (pos, emitted) = super::annotation_dump::trim(start, emitted);
        tandem(
            &emitted,
            &region.padded_reference[(pos + 1 - region.reference_start) as usize..],
            &mut fields,
        );
        Ok(fields)
    }
    pub(super) fn region(
        &mut self,
        formatted: &[u8],
        phase: &[u8],
        extras: Vec<Fields>,
    ) -> io::Result<()> {
        let mut records = parse(formatted);
        for ((r, p), extra) in records.iter_mut().zip(parse(phase)).zip(extras) {
            r.info.extend(extra);
            for (s, fields) in p.samples {
                r.samples.get_mut(&s).unwrap().extend(fields);
            }
        }
        for r in records {
            let (chrom, pos) = r.locus.rsplit_once(':').unwrap();
            let (reference, alts) = r.alleles.split_once(',').unwrap();
            let info = r
                .info
                .iter()
                .map(|(k, v)| {
                    if v.is_empty() {
                        k.clone()
                    } else {
                        format!("{k}={v}")
                    }
                })
                .collect::<Vec<_>>()
                .join(";");
            let mut keys: Vec<_> = r.samples.values().flat_map(|s| s.keys().cloned()).collect();
            keys.sort();
            keys.dedup();
            keys.retain(|k| k != "GT");
            keys.insert(0, "GT".into());
            write!(
                self.output,
                "{chrom}\t{pos}\t.\t{reference}\t{alts}\t.\t.\t{info}\t{}",
                keys.join(":")
            )?;
            for s in r.samples.values() {
                write!(
                    self.output,
                    "\t{}",
                    keys.iter()
                        .map(|k| s.get(k).map_or(".", String::as_str))
                        .collect::<Vec<_>>()
                        .join(":")
                )?;
            }
            writeln!(self.output)?;
        }
        Ok(())
    }
}
pub fn run(inputs: &Inputs<'_>, output: &Path, f1r2: Option<&Path>) -> io::Result<()> {
    let mut writer = Writer {
        output: ParallelBgzfWriter::new(File::create(output)?, 1),
        pon: Tabix::open(inputs.panel_of_normals)?,
        population: Tabix::open(inputs.germline_resource)?,
        normal: inputs.normal_sample.into(),
        meta: include_str!("header.txt")
            .lines()
            .map(|l| l[2..].into())
            .collect(),
        contigs: nomics_io::fasta::read_fai(Path::new(&format!(
            "{}.fai",
            inputs.reference.display()
        )))?
        .into_iter()
        .map(|e| (e.name, e.len))
        .collect(),
    };
    let mut samples = std::collections::BTreeSet::new();
    for path in inputs.input {
        let h = nomics_io::bam_query::BamQuery::open(path, &[])?.header;
        for line in h
            .text
            .split(|&b| b == b'\n')
            .filter(|l| l.starts_with(b"@RG\t"))
        {
            if let Some(sm) = line
                .split(|&b| b == b'\t')
                .find_map(|f| f.strip_prefix(b"SM:"))
            {
                samples.insert(String::from_utf8_lossy(sm).into_owned());
            }
        }
    }
    if samples.len() != 2 || !samples.contains(inputs.normal_sample) {
        return Err(io::Error::other(
            "only one named normal and one tumour sample are supported",
        ));
    }
    writer.meta.push(format!("normal_sample={}", writer.normal));
    writer.meta.push(format!(
        "tumor_sample={}",
        samples.iter().find(|s| **s != writer.normal).unwrap()
    ));
    let header = nomics_hc::vcf::Header {
        meta: writer.meta.clone(),
        contigs: writer.contigs.clone(),
        sample: samples.into_iter().collect::<Vec<_>>().join("\t"),
    };
    writer.output.write_all(header.text().as_bytes())?;
    let mut annotations = AnnotationWriters {
        raw: Vec::new(),
        formatted: Vec::new(),
        phasing: Vec::new(),
    };
    let mut stats = GenotypingStats::default();
    // The f1r2 counts come from the same pass over the reads as the calls, as in GATK.
    let mut collector = f1r2.map(|_| crate::mutect2_f1r2::Collector::default());
    let mut no_loci = crate::mutect2_regions::NoLoci;
    let loci: &mut dyn crate::mutect2_regions::LocusSink = match collector.as_mut() {
        Some(c) => c,
        None => &mut no_loci,
    };
    let (traversal, _) = mutect2_likelihoods::run_with_sinks(
        inputs,
        Path::new("/dev/null"),
        Path::new("/dev/null"),
        &mut |region| {
            dump(
                region,
                inputs.normal_sample,
                &mut io::sink(),
                &mut stats,
                &mut annotations,
                Some(&mut writer),
            )
        },
        loci,
    )?;
    writer.output.finish()?;
    if let (Some(collector), Some(path)) = (&collector, f1r2) {
        collector.write_archive(path)?;
    }
    eprint!("{}", traversal.filter_summary());
    Ok(())
}

pub(super) fn event(start: i64, reference: &[u8], alternate: &[u8]) -> nomics_hc::trim::Event {
    let mut r = reference.to_vec();
    let mut a = alternate.to_vec();
    if r.len() != 1 && a.len() != 1 {
        let n = r
            .iter()
            .rev()
            .zip(a.iter().rev())
            .take_while(|(x, y)| x == y)
            .count();
        r.truncate(r.len() - n);
        a.truncate(a.len() - n);
    }
    nomics_hc::trim::Event {
        start,
        reference: r,
        alternate: a,
    }
}
pub(super) fn event_counts(
    region: &LikelihoodRegion<'_>,
    events: &[Vec<nomics_hc::trim::Event>],
    calls: &[super::phasing::Call],
    somatic: &[nomics_hc::trim::Event],
    extras: &mut [Fields],
) {
    let mut support = vec![0; events.len()];
    let order: Vec<_> = (0..events.len()).collect();
    for sample in region.samples {
        let mut values = sample.values.clone();
        super::fragments::cap(&mut values);
        for r in 0..sample.reads.len() {
            support[super::annotations::best_with_reference(
                &values,
                r,
                &order,
                region.haplotypes.iter().position(|h| h.is_ref).unwrap(),
            )
            .0] += 1;
        }
    }
    for (call, info) in calls.iter().zip(extras) {
        info.insert("ECNT".into(), somatic.len().to_string());
        let counts: Vec<_> = call.alleles[1..]
            .iter()
            .filter_map(|a| {
                let e = event(call.start, &call.alleles[0], a);
                let mut best = None;
                for (h, es) in events.iter().enumerate() {
                    if es.contains(&e) && best.is_none_or(|b| support[h] > support[b]) {
                        best = Some(h);
                    }
                }
                best.map(|h| {
                    events[h]
                        .iter()
                        .filter(|e| somatic.contains(e))
                        .count()
                        .to_string()
                })
            })
            .collect();
        if !counts.is_empty() {
            info.insert("ECNTH".into(), counts.join(","));
        }
    }
}
fn tandem(alleles: &[Vec<u8>], context: &[u8], fields: &mut Fields) {
    // VariantContext.isIndel excludes MIXED indel/SNP or indel/MNP sites.
    if alleles[1..].iter().any(|a| a.len() == alleles[0].len()) {
        return;
    }
    let reference = &alleles[0][1..];
    let mut counts = Vec::new();
    let mut unit = Vec::new();
    for a in &alleles[1..] {
        let alt = &a[1..];
        let long = if alt.len() > reference.len() {
            alt
        } else {
            reference
        };
        let n = (1..=long.len())
            .find(|&n| long.len() % n == 0 && long.chunks(n).all(|c| c == &long[..n]))
            .unwrap_or(1);
        unit = long.get(..n).unwrap_or(&[0]).to_vec();
        let reps = |s: &[u8]| s.chunks(unit.len()).take_while(|c| *c == unit).count();
        let base = reps(reference) as i32;
        let nr = reps(&[reference, context].concat()) as i32 - base;
        let na = reps(&[alt, context].concat()) as i32 - base;
        if nr == 0 || na == 0 {
            return;
        }
        if counts.is_empty() {
            counts.push(nr.to_string());
        }
        counts.push(na.to_string());
    }
    fields.insert("RPA".into(), counts.join(","));
    fields.insert("RU".into(), String::from_utf8_lossy(&unit).into());
    fields.insert("STR".into(), String::new());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resource_missing_values_do_not_mask_malformed_numbers() {
        assert_eq!(resource_af_values(b".,0.25").unwrap(), vec![1e-6, 0.25]);
        assert!(resource_af_values(b"0.25,invalid").is_err());
    }

    #[test]
    fn mixed_mnp_indel_is_not_annotated_as_a_tandem_repeat() {
        // htsjdk VariantContext.determinePolymorphicType: INDEL + MNP = MIXED.
        // TandemRepeat.annotate rejects this before computing repeat units.
        let mut fields = Fields::new();
        tandem(
            &[b"AT".to_vec(), b"A".to_vec(), b"GT".to_vec()],
            b"TTTTG",
            &mut fields,
        );
        assert!(fields.is_empty());
    }
}
