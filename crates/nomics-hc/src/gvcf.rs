//! GATK 4.6.2.0 GVCFBlockCombiner / HomRefBlock and physical phasing.
use crate::trim::{TrimHap, events_for_haplotype};
use crate::vcf::{InfoValue, Record};
use std::collections::{BTreeMap, BTreeSet};

pub const NON_REF: &str = "<NON_REF>";

pub fn bands() -> Vec<i32> {
    (1..=60).chain([70, 80, 90, 99, 100]).collect()
}

fn band(gq: i32) -> usize {
    match gq {
        0..=59 => gq as usize,
        60..=69 => 60,
        70..=79 => 61,
        80..=89 => 62,
        90..=98 => 63,
        _ => 64,
    }
}

struct Block {
    record: Record,
    end: i64,
    depths: BTreeMap<i32, u64>,
    band: usize,
}

/// One instance per scattered HC invocation. Never recombine gathered groups.
#[derive(Default)]
pub struct BlockCombiner {
    block: Option<Block>,
    next_available: Option<(usize, i64)>,
    records: Vec<Record>,
}

impl BlockCombiner {
    /// `submit(site.to_record(contig))` without building the record unless it
    /// starts a block: a reference site is hom-ref `<NON_REF>` with 3 PLs.
    pub fn submit_ref(&mut self, contig: usize, site: &crate::refconfidence::RefSite) {
        if self
            .block
            .as_ref()
            .is_some_and(|b| b.record.contig != contig || site.pos > b.end + 1)
        {
            self.flush();
        }
        if let Some((c, end)) = self.next_available {
            if contig == c && site.pos <= end {
                return;
            }
            self.next_available = None;
        }
        let b = band(site.gq);
        if let Some(block) = self.block.as_mut()
            && block.band == b
            && block.record.genotype.pl.len() == 3
        {
            assert_eq!(site.pos, block.end + 1, "overlapping reference loci");
            for (old, new) in block.record.genotype.pl.iter_mut().zip(&site.pl) {
                *old = (*old).min(*new);
            }
            block.end = site.pos;
            *block.depths.entry(site.dp.max(0)).or_default() += 1;
        } else {
            self.flush();
            self.block = Some(Block {
                end: site.pos,
                depths: BTreeMap::from([(site.dp.max(0), 1)]),
                band: b,
                record: site.to_record(contig),
            });
        }
    }

    pub fn submit(&mut self, record: Record) {
        if self
            .block
            .as_ref()
            .is_some_and(|b| b.record.contig != record.contig || record.pos > b.end + 1)
        {
            self.flush();
        }
        let g = &record.genotype;
        let hom_ref = g.alleles.iter().all(|&a| a == 0)
            || (g
                .alleles
                .iter()
                .all(|&a| a == crate::genotype::NO_CALL_ALLELE)
                && g.pl.first() == Some(&0));
        if hom_ref && record.alts == [NON_REF] {
            if let Some((c, end)) = self.next_available {
                if record.contig == c && record.pos <= end {
                    return;
                }
                self.next_available = None;
            }
            let b = band(g.gq);
            if let Some(block) = self.block.as_mut()
                && block.band == b
                && block.record.genotype.pl.len() == g.pl.len()
            {
                assert_eq!(record.pos, block.end + 1, "overlapping reference loci");
                for (old, new) in block.record.genotype.pl.iter_mut().zip(&g.pl) {
                    *old = (*old).min(*new);
                }
                block.end = record.pos;
                *block.depths.entry(g.dp.max(0)).or_default() += 1;
            } else {
                self.flush();
                self.block = Some(Block {
                    end: record.pos,
                    depths: BTreeMap::from([(g.dp.max(0), 1)]),
                    band: b,
                    record,
                });
            }
        } else {
            self.flush();
            self.next_available = Some((
                record.contig,
                record.pos + record.ref_allele.len() as i64 - 1,
            ));
            self.records.push(record);
        }
    }

    fn flush(&mut self) {
        if let Some(mut b) = self.block.take() {
            let n: u64 = b.depths.values().sum();
            let at = |index| {
                let mut seen = 0;
                for (&depth, &count) in &b.depths {
                    seen += count;
                    if seen > index {
                        return i64::from(depth);
                    }
                }
                unreachable!("nonempty block depth histogram")
            };
            let dp = (at((n - 1) / 2) + at(n / 2) + 1) / 2;
            let g = &mut b.record.genotype;
            g.alleles = vec![0, 0];
            g.ad.clear();
            g.dp = dp as i32;
            let mut pl = g.pl.clone();
            pl.sort_unstable();
            g.gq = pl[1] - pl[0];
            g.gvcf = Some(crate::vcf::GvcfGenotype {
                min_dp: b.depths.first_key_value().map(|(&depth, _)| depth),
                ..Default::default()
            });
            b.record.info = vec![("END".into(), InfoValue::Ints(vec![b.end]))];
            self.records.push(b.record);
        }
    }

    pub fn finish(mut self) -> Vec<Record> {
        self.flush();
        self.records
    }
}

/// AssemblyBasedCallerUtils.phaseCalls. Sets are used only for membership,
/// cardinality and intersection; Java HashSet iteration cannot affect the result.
pub(crate) fn phase_calls(
    calls: &mut [Record],
    haplotypes: &[TrimHap],
    reference: &[u8],
    reference_start: i64,
) {
    let events: Vec<_> = haplotypes
        .iter()
        .map(|h| events_for_haplotype(h, reference, reference_start))
        .collect();
    let map: Vec<BTreeSet<usize>> = calls
        .iter()
        .map(|c| {
            let alts: Vec<_> = c
                .alts
                .iter()
                .filter(|a| a.as_str() != NON_REF && a.as_str() != "*")
                .collect();
            if alts.len() != 1 {
                return BTreeSet::new();
            }
            events
                .iter()
                .enumerate()
                .filter(|(_, ev)| {
                    ev.iter()
                        .any(|e| e.start == c.pos && e.alternate == alts[0].as_bytes())
                })
                .map(|(i, _)| i)
                .collect()
        })
        .collect();
    let total = map.iter().flatten().copied().collect::<BTreeSet<_>>().len();
    let mut phase: Vec<Option<(usize, bool)>> = vec![None; calls.len()];
    let mut groups = 0;
    for i in 0..calls.len().saturating_sub(1) {
        if map[i].is_empty() {
            continue;
        }
        let all = map[i].len() == total;
        let mut available = map[i].clone();
        for j in i + 1..calls.len() {
            if map[j].is_empty() {
                continue;
            }
            let together =
                map[i] == map[j] || (all && map[j].is_subset(&available)) || map[j].len() == total;
            let apart =
                !together && map[i].len() + map[j].len() == total && map[i].is_disjoint(&map[j]);
            if !together && !apart {
                continue;
            }
            if let Some((group, reverse)) = phase[i] {
                if phase[j].is_none() {
                    phase[j] = Some((group, reverse ^ apart));
                }
            } else {
                if phase[j].is_some() {
                    return;
                }
                phase[i] = Some((groups, false));
                phase[j] = Some((groups, apart));
                groups += 1;
                if together {
                    available.retain(|h| map[j].contains(h));
                }
            }
        }
    }
    for group in 0..groups {
        let first = phase
            .iter()
            .position(|p| p.is_some_and(|p| p.0 == group))
            .unwrap();
        let c = &calls[first];
        let ps = c.pos;
        let pid = format!("{}_{}_{}", ps, c.ref_allele, c.alts[0]);
        for (c, p) in calls.iter_mut().zip(&phase) {
            let Some((id, reverse)) = p else {
                continue;
            };
            if *id != group {
                continue;
            }
            let index = usize::from(!reverse);
            let g = &mut c.genotype;
            let a = g.alleles[index];
            let site_specific = a != 0
                && a != crate::genotype::NO_CALL_ALLELE
                && c.alts[a - 1] != NON_REF
                && c.alts[a - 1] != "*";
            if g.alleles[0] != g.alleles[1] && !site_specific {
                g.alleles.reverse();
            }
            g.gvcf.as_mut().unwrap().phase =
                Some((if *reverse { "1|0" } else { "0|1" }.into(), pid.clone(), ps));
        }
    }
}

pub fn header_meta(command: &str) -> Vec<String> {
    let mut meta = crate::caller::header_meta(command);
    meta.retain(|m| {
        !["AC", "AF", "AN", "FS", "SOR", "QD", "MQ"]
            .iter()
            .any(|k| m.starts_with(&format!("INFO=<ID={k},")))
    });
    meta.extend([
        r#"ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">"#,
        r#"FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">"#,
        r#"FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">"#,
        r#"FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">"#,
        r#"FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">"#,
        r#"FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">"#,
        r#"INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">"#,
        r#"INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">"#,
    ].into_iter().map(str::to_owned));
    let mut lower = 0;
    for upper in bands() {
        meta.push(format!(
            "GVCFBlock{lower}-{upper}=minGQ={lower}(inclusive),maxGQ={upper}(exclusive)"
        ));
        lower = upper;
    }
    meta
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{Genotype, GvcfGenotype, encode};

    fn hom_ref(pos: i64, dp: i32, pl: [i32; 3]) -> Record {
        Record {
            contig: 0,
            pos,
            id: None,
            ref_allele: "A".into(),
            alts: vec![NON_REF.into()],
            qual: f64::NAN,
            info: vec![],
            genotype: Genotype {
                alleles: vec![0, 0],
                ad: vec![dp, 0],
                dp,
                gq: pl[1].min(pl[2]),
                pl: pl.to_vec(),
                gvcf: Some(GvcfGenotype::default()),
            },
        }
    }

    fn text(record: &Record) -> String {
        let mut bytes = Vec::new();
        encode(record, &[("chr1".into(), 1000)], &mut bytes);
        String::from_utf8(bytes).unwrap()
    }

    #[test]
    fn blocks_take_component_minima_and_round_median_depth_up() {
        let mut writer = BlockCombiner::default();
        writer.submit(hom_ref(10, 10, [0, 61, 100]));
        writer.submit(hom_ref(11, 3, [0, 68, 95]));
        writer.submit(hom_ref(12, 7, [0, 70, 100])); // new GQ band
        writer.submit(hom_ref(14, 0, [0, 0, 0])); // gap
        let records = writer.finish();
        assert_eq!(records.len(), 3);
        assert_eq!(
            text(&records[0]),
            "chr1\t10\t.\tA\t<NON_REF>\t.\t.\tEND=11\tGT:DP:GQ:MIN_DP:PL\t0/0:7:61:3:0,61,95\n"
        );
        assert_eq!(
            text(&records[2]),
            "chr1\t14\t.\tA\t<NON_REF>\t.\t.\tEND=14\tGT:DP:GQ:MIN_DP:PL\t0/0:0:0:0:0,0,0\n"
        );
        assert_eq!(band(99), band(1000));
        assert_ne!(band(98), band(99));
    }

    #[test]
    fn deletion_suppresses_covered_reference_but_not_overlapping_variants() {
        let mut writer = BlockCombiner::default();
        let mut deletion = hom_ref(10, 10, [0, 30, 90]);
        deletion.ref_allele = "AAA".into();
        deletion.alts.insert(0, "A".into());
        deletion.genotype.alleles = vec![0, 1];
        writer.submit(deletion);
        writer.submit(hom_ref(11, 10, [0, 30, 90]));
        writer.submit(hom_ref(12, 10, [0, 30, 90]));
        writer.submit(hom_ref(13, 10, [0, 30, 90]));
        let records = writer.finish();
        assert_eq!(records.iter().map(|r| r.pos).collect::<Vec<_>>(), [10, 13]);
    }

    #[test]
    fn physical_phase_cis_and_trans() {
        let reference = b"AACCGG";
        let hap = |bases: &[u8], is_ref| TrimHap {
            bases: bases.to_vec(),
            is_ref,
            alignment_start_hap_wrt_ref: 0,
            cigar: vec![6 << 4],
            contig: "chr1".into(),
            start: 1,
            end: 6,
        };
        let calls = || {
            [2, 4]
                .into_iter()
                .map(|pos| {
                    let mut c = hom_ref(pos, 10, [0, 30, 90]);
                    c.ref_allele = (reference[(pos - 1) as usize] as char).to_string();
                    c.alts.insert(0, "T".into());
                    c.genotype.alleles = vec![0, 1];
                    c
                })
                .collect::<Vec<_>>()
        };
        let mut cis = calls();
        phase_calls(
            &mut cis,
            &[hap(reference, true), hap(b"ATCTGG", false)],
            reference,
            1,
        );
        for c in &cis {
            assert_eq!(
                c.genotype.gvcf.as_ref().unwrap().phase,
                Some(("0|1".into(), "2_A_T".into(), 2))
            );
        }
        let mut trans = calls();
        phase_calls(
            &mut trans,
            &[
                hap(reference, true),
                hap(b"ATCCGG", false),
                hap(b"AACTGG", false),
            ],
            reference,
            1,
        );
        assert_eq!(trans[0].genotype.alleles, [0, 1]);
        assert_eq!(trans[1].genotype.alleles, [1, 0]);
        assert_eq!(
            trans[1].genotype.gvcf.as_ref().unwrap().phase,
            Some(("1|0".into(), "2_A_T".into(), 2))
        );
    }
}
