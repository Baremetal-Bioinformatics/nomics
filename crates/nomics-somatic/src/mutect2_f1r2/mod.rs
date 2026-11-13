//! GATK 4.6.2.0 F1R2CountsCollector, default single-tumour sarek mode.
//!
//! Final-caller integration: construct `Collector::default()`, pass it together
//! with the caller's RegionSink to `mutect2_regions::run_with_sinks`, then call
//! `Collector::write_archive`. No second BAM traversal is needed. The traversal
//! initializes the sample from BAM headers, including when no reads survive.
//! Histograms use depth 1..200; member order is ref_histogram, alt_histogram,
//! alt_table. Singleton columns have deterministic context/base/orientation order
//! (the reference's JVM enum identity-hash order is explicitly normalized).
use crate::{
    mutect2_regions::{self, LocusSink, Read, activity},
    orientation::{hash, kmer, tar, url_encode},
};
use nomics_hc::activity::Region;
use std::{
    io::{self, Write},
    path::{Path, PathBuf},
};

const DEPTH: usize = 200;

/// The lossless AltSiteRecord boundary. The learner's internal `Record` retains
/// only n/k/x/alt and cannot represent ref_f1r2 or the original reference context.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AltSiteRecord {
    pub context: usize,
    pub ref_count: u32,
    pub alt_count: u32,
    pub ref_f1r2: u32,
    pub alt_f1r2: u32,
    pub alt: usize,
}

/// Counts from the activity pileup, not the later assembly/annotation pileup.
pub struct Collector {
    sample: Option<String>,
    reference: [Box<[u64; DEPTH]>; 64],
    singleton: [Box<[u64; DEPTH]>; 64 * 8],
    records: Vec<Vec<AltSiteRecord>>,
    mapping_qualities: Vec<u8>,
}
impl Default for Collector {
    fn default() -> Self {
        Self {
            sample: None,
            reference: std::array::from_fn(|_| Box::new([0; DEPTH])),
            singleton: std::array::from_fn(|_| Box::new([0; DEPTH])),
            records: Vec::new(),
            mapping_qualities: Vec::new(),
        }
    }
}

fn base(b: u8) -> Option<usize> {
    // BaseUtils.simpleBaseToBaseIndex includes lowercase and '*' as A.
    match b {
        b'A' | b'a' | b'*' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}
fn context_order() -> Vec<usize> {
    // HashMap<>(64) grows to 128 buckets at insertion 49. Stable collision
    // chains preserve ALL_KMERS' lexicographic insertion order across resize.
    let mut contexts: Vec<_> = (0..64).collect();
    contexts.sort_by_key(|&c| hash(&kmer(c)) & 127);
    contexts
}
impl LocusSink for Collector {
    fn on_samples(&mut self, samples: &[String]) -> io::Result<()> {
        if samples.len() != 1 || self.sample.is_some() {
            return Err(io::Error::other(
                "f1r2 requires one tumour and one traversal per collector",
            ));
        }
        self.sample = Some(samples[0].clone());
        Ok(())
    }
    fn on_locus(
        &mut self,
        position: i64,
        reference: &[u8],
        reads: &[Read],
        active: &[usize],
    ) -> io::Result<()> {
        let Some(sample) = &self.sample else {
            return Err(io::Error::other(
                "f1r2 collector has not received sample headers",
            ));
        };
        // ReferenceContext.getKmerAround requires the entire 3-base window.
        if position < 2 || position as usize >= reference.len() {
            return Ok(());
        }
        let context = &reference[position as usize - 2..=position as usize];
        let Some(c) = base(context[0])
            .zip(base(context[1]))
            .zip(base(context[2]))
            .map(|((a, b), c)| a * 16 + b * 4 + c)
        else {
            return Ok(());
        };
        let ref_base = c / 4 % 4;
        let (mut counts, mut f1r2) = ([0u32; 4], [0u32; 4]);
        let mut indels = 0u32;
        self.mapping_qualities.clear();
        for &id in active {
            let r = &reads[id];
            if r.is_normal || r.sample != *sample || !activity::covered(r, position) {
                continue;
            }
            let (i, off, within) = activity::element(r, position).unwrap();
            let a = &r.alignment;
            let cigar = a.cigar[i];
            // PileupElement.getQual returns 16 for deletions, which fail >20.
            if cigar & 15 == 2 || (a.quals[off] as i8) <= 20 {
                continue;
            }
            self.mapping_qualities.push(a.mapq);
            let at_end = within == i64::from(cigar >> 4) - 1;
            let before_deletion = at_end
                && a.cigar[i + 1..]
                    .iter()
                    .find(|&&c| matches!(c & 15, 0 | 2 | 7 | 8))
                    .is_some_and(|&c| c & 15 == 2);
            let after_insertion = within == 0 && i > 0 && a.cigar[i - 1] & 15 == 1;
            if before_deletion || after_insertion {
                indels += 1;
            }
            if let Some(b) = base(a.bases[off]) {
                counts[b] += 1;
                // GATKRead.isFirstOfPair is false for unpaired reads.
                let first = a.flag & 1 != 0 && a.flag & 64 != 0;
                if (a.flag & 16 != 0) != first {
                    // Nucleotide.decode('*') is not A: wildcard participates
                    // in base counts but never in the orientation count.
                    if a.bases[off] != b'*' {
                        f1r2[b] += 1;
                    }
                }
            }
        }
        let depth: u32 = counts.iter().sum();
        if depth == 0 || indels > depth / 100 {
            return Ok(());
        }
        self.mapping_qualities.sort_unstable();
        let n = self.mapping_qualities.len();
        // MathUtils.median uses the average of the middle two for even size.
        // Comparison to integer 50 is exactly equivalent without floating math.
        if u16::from(self.mapping_qualities[(n - 1) / 2]) + u16::from(self.mapping_qualities[n / 2])
            < 100
        {
            return Ok(());
        }
        let mut alt = if ref_base == 0 { 1 } else { 0 };
        for b in 0..4 {
            if b != ref_base && counts[b] > counts[alt] {
                alt = b;
            }
        }
        let bin = depth.min(DEPTH as u32) as usize - 1;
        match counts[alt] {
            0 => self.reference[c][bin] += 1,
            1 => self.singleton[c * 8 + alt * 2 + usize::from(f1r2[alt] != 1)][bin] += 1,
            alt_count => {
                if self.records.last().is_none_or(|chunk| chunk.len() == 256) {
                    self.records.push(Vec::with_capacity(256));
                }
                self.records.last_mut().unwrap().push(AltSiteRecord {
                    context: c,
                    ref_count: counts[ref_base],
                    alt_count,
                    ref_f1r2: f1r2[ref_base],
                    alt_f1r2: f1r2[alt],
                    alt,
                });
            }
        }
        Ok(())
    }
}
impl Collector {
    fn histogram(&self, alt: bool) -> io::Result<Vec<u8>> {
        let sample = self
            .sample
            .as_ref()
            .ok_or_else(|| io::Error::other("missing f1r2 sample"))?;
        let mut out = Vec::new();
        write!(
            out,
            "## htsjdk.samtools.metrics.StringHeader\n# {sample}\n\n\n## HISTOGRAM\tjava.lang.Integer\ndepth"
        )?;
        let contexts = context_order();
        let mut columns = Vec::new();
        for &c in &contexts {
            if alt {
                for b in 0..4 {
                    if b == c / 4 % 4 {
                        continue;
                    }
                    for (o, name) in ["F1R2", "F2R1"].iter().enumerate() {
                        write!(out, "\t{}_{}_{name}", kmer(c), b"ACGT"[b] as char)?;
                        columns.push(&self.singleton[c * 8 + b * 2 + o]);
                    }
                }
            } else {
                write!(out, "\t{}", kmer(c))?;
                columns.push(&self.reference[c]);
            }
        }
        writeln!(out)?;
        for d in 0..DEPTH {
            write!(out, "{}", d + 1)?;
            for col in &columns {
                write!(out, "\t{}", col[d])?;
            }
            writeln!(out)?;
        }
        writeln!(out)?;
        Ok(out)
    }
    /// Serialize all three payloads, using the orientation module's checked tar
    /// writer and native deflater. The wrapper metadata is deterministic.
    pub fn write_archive(&self, output: &Path) -> io::Result<()> {
        let sample = self
            .sample
            .as_ref()
            .ok_or_else(|| io::Error::other("missing f1r2 sample"))?;
        let stem = url_encode(sample);
        let mut table = Vec::new();
        writeln!(
            table,
            "#<METADATA>SAMPLE={sample}\ncontext\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt"
        )?;
        for r in self.records.iter().flatten() {
            writeln!(
                table,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                kmer(r.context),
                r.ref_count,
                r.alt_count,
                r.ref_f1r2,
                r.alt_f1r2,
                r.ref_count + r.alt_count,
                b"ACGT"[r.alt] as char
            )?;
        }
        tar::write(
            output,
            &[
                (format!("./{stem}.ref_histogram"), self.histogram(false)?),
                (format!("./{stem}.alt_histogram"), self.histogram(true)?),
                (format!("./{stem}.alt_table"), table),
            ],
        )
    }
}

pub struct Inputs<'a> {
    pub input: &'a [PathBuf],
    pub normal_sample: &'a str,
    pub reference: &'a Path,
    pub intervals: &'a Path,
    pub panel_of_normals: &'a Path,
    pub germline_resource: &'a Path,
    pub output: &'a Path,
}
/// Standalone collector. The final caller should use `Collector` directly with
/// `run_with_sinks` instead, sharing its own traversal and region calls.
pub fn run(inputs: &Inputs<'_>) -> io::Result<mutect2_regions::Stats> {
    let temporary = std::env::temp_dir().join(format!(
        "nomics-f1r2-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(io::Error::other)?
            .as_nanos()
    ));
    std::fs::create_dir(&temporary)?;
    let result = (|| {
        let mut collector = Collector::default();
        let stats = mutect2_regions::run_with_sinks(
            &mutect2_regions::Inputs {
                input: inputs.input,
                normal_sample: inputs.normal_sample,
                reference: inputs.reference,
                intervals: inputs.intervals,
                panel_of_normals: inputs.panel_of_normals,
                germline_resource: inputs.germline_resource,
                regions: &temporary.join("regions"),
                stats: &temporary.join("stats"),
                dump_prefix: None,
            },
            &mut |_: &Region, _: &[&Read]| Ok(()),
            &mut collector,
        )?;
        collector.write_archive(inputs.output)?;
        Ok(stats)
    })();
    std::fs::remove_dir_all(temporary)?;
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use nomics_hc::finalize::HcRead;

    fn read(base: u8, qual: u8, mapq: u8, cigar: Vec<u32>) -> Read {
        Read {
            alignment: HcRead {
                name: b"test".to_vec(),
                flag: 65,
                contig: 0,
                start: 2,
                cigar,
                bases: vec![base; 3],
                quals: vec![qual; 3],
                mapq,
                mate_contig: -1,
                mate_start: 0,
                tlen: 0,
                record: 0,
            },
            source: (0, 0),
            sample: "tumour".into(),
            is_normal: false,
            overlap_indels: false,
            insertion_quals: vec![45; 3],
            deletion_quals: vec![45; 3],
        }
    }
    fn collect(reads: &[Read], position: i64) -> Collector {
        let mut collector = Collector::default();
        collector.on_samples(&["tumour".into()]).unwrap();
        collector
            .on_locus(
                position,
                b"AAAAAA",
                reads,
                &(0..reads.len()).collect::<Vec<_>>(),
            )
            .unwrap();
        collector
    }
    #[test]
    fn quality_filtered_median_signed_quality_and_sample_selection() {
        let mut reads = vec![
            read(b'A', 21, 49, vec![48]),
            read(b'A', 21, 51, vec![48]),
            read(b'C', 20, 60, vec![48]),
            read(b'C', 255, 60, vec![48]),
            read(b'N', 21, 0, vec![48]),
        ];
        let mut normal = read(b'C', 40, 60, vec![48]);
        normal.is_normal = true;
        reads.push(normal);
        // N contributes to median MQ but not depth; lower-quality bases and
        // normal reads contribute to neither. Median [0,49,51] is 49: reject.
        assert_eq!(collect(&reads, 3).reference[0][1], 0);
        reads.remove(4);
        // Even-sized median is exactly 50, which passes; Q20/255 are excluded.
        assert_eq!(collect(&reads, 3).reference[0][1], 1);
    }
    #[test]
    fn dominant_alt_ties_table_depth_and_orientation() {
        let reads: Vec<_> = b"ACCTT"
            .iter()
            .copied()
            .map(|b| read(b, 30, 60, vec![48]))
            .collect();
        let collector = collect(&reads, 3);
        assert_eq!(
            collector.records[0],
            vec![AltSiteRecord {
                context: 0,
                ref_count: 1,
                alt_count: 2,
                ref_f1r2: 1,
                alt_f1r2: 2,
                alt: 1,
            }]
        ); // C beats T; serialized depth is 3 although pileup depth is 5.
        let mut unpaired = read(b'C', 30, 60, vec![48]);
        unpaired.alignment.flag = 64; // stale first-of-pair flag is ignored.
        assert_eq!(collect(&[unpaired], 3).singleton[3][0], 1); // F2R1
    }
    #[test]
    fn insertion_side_and_integer_indel_threshold() {
        let inserted = read(b'A', 30, 60, vec![16, 17, 16]); // 1M1I1M, insertion AFTER pos2.
        assert_eq!(
            collect(std::slice::from_ref(&inserted), 2).reference[0][0],
            1
        );
        assert_eq!(
            collect(std::slice::from_ref(&inserted), 3).reference[0][0],
            0
        );
        let mut reads = vec![read(b'A', 30, 60, vec![48]); 99];
        reads[0] = inserted;
        assert_eq!(collect(&reads, 3).reference[0][98], 0); // 1 > 99/100
        reads.push(read(b'A', 30, 60, vec![48]));
        assert_eq!(collect(&reads, 3).reference[0][99], 1); // 1 == 100/100
        let deleted = read(b'A', 30, 60, vec![16, 18, 16]);
        assert_eq!(
            collect(std::slice::from_ref(&deleted), 2).reference[0][0],
            0
        ); // before deletion
        assert_eq!(
            collect(std::slice::from_ref(&deleted), 3).reference[0][0],
            0
        ); // deletion Q16
        assert_eq!(collect(&[deleted], 4).reference[0][0], 1); // after deletion is allowed
    }
}
