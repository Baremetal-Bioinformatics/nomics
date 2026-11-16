use super::{cluster::Datum, invalid};
use std::{collections::BTreeMap, io};
pub(super) type Attributes = BTreeMap<String, String>;
#[derive(Clone)]
pub(super) struct Record {
    pub fields: Vec<String>,
    pub info: Attributes,
    pub gt: Vec<Attributes>,
    pub pos: u64,
    pub alts: Vec<String>,
    pub ad: Vec<u32>,
    pub normal_ad: Vec<u32>,
    pub lod: Vec<f64>,
    pub orientation: f64,
    pub maf: f64,
    pub phase: Option<String>,
}
pub(super) fn floats(s: &str) -> Vec<f64> {
    s.split(',')
        .map(|v| v.parse().expect("validated numeric VCF attribute"))
        .collect()
}
pub(super) fn ints(s: &str) -> Vec<u32> {
    s.split(',')
        .map(|v| v.parse().expect("validated integer VCF attribute"))
        .collect()
}
impl Record {
    pub fn parse(line: &str, tumor: usize, normal: usize) -> io::Result<Self> {
        let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if fields.len() != 11 {
            return Err(invalid("only two-sample tumour-normal VCF supported"));
        }
        let info: Attributes = fields[7]
            .split(';')
            .map(|s| s.split_once('=').unwrap_or((s, "")))
            .map(|(k, v)| (k.to_owned(), v.to_owned()))
            .collect();
        let keys: Vec<_> = fields[8].split(':').collect();
        let gt: Vec<Attributes> = fields[9..]
            .iter()
            .map(|g| {
                keys.iter()
                    .copied()
                    .zip(g.split(':'))
                    .map(|(k, v)| (k.to_owned(), v.to_owned()))
                    .collect()
            })
            .collect();
        let alts: Vec<_> = fields[4].split(',').map(str::to_owned).collect();
        if alts
            .iter()
            .any(|s| s.starts_with('<') || s == "*" || s.contains(['[', ']']))
        {
            return Err(invalid("symbolic alleles are not supported"));
        }
        let n = alts.len();
        for (key, offset) in [
            ("TLOD", 0),
            ("POPAF", 0),
            ("NLOD", 0),
            ("NALOD", 0),
            ("MPOS", 0),
            ("AS_UNIQ_ALT_READ_COUNT", 0),
            ("MBQ", 1),
            ("MMQ", 1),
            ("MFRL", 1),
            ("ECNT", usize::MAX),
            ("ECNTH", usize::MAX),
            ("NCount", usize::MAX),
        ] {
            if let Some(value) = info.get(key) {
                let parts: Vec<_> = value.split(',').collect();
                if parts
                    .iter()
                    .any(|v| !v.parse::<f64>().is_ok_and(f64::is_finite))
                    || (offset != usize::MAX && parts.len() != n + offset)
                {
                    return Err(invalid(format!(
                        "unsupported malformed INFO/{key} at {}:{}",
                        fields[0], fields[1]
                    )));
                }
            }
        }
        let valid_ints = |value: &str| {
            value
                .split(',')
                .all(|v| v.parse::<u32>().is_ok_and(|n| n <= i32::MAX as u32))
        };
        for g in &gt {
            for key in ["AD", "F1R2", "F2R1"] {
                if let Some(value) = g.get(key)
                    && (!valid_ints(value)
                        || value.split(',').count() != n + 1
                        || ints(value).iter().map(|&v| v as u64).sum::<u64>() > i32::MAX as u64)
                {
                    return Err(invalid(format!("unsupported malformed FORMAT/{key}")));
                }
            }
            if let Some(value) = g.get("AF")
                && (value.split(',').count() != n
                    || value
                        .split(',')
                        .any(|v| !v.parse::<f64>().is_ok_and(f64::is_finite)))
            {
                return Err(invalid("unsupported malformed FORMAT/AF"));
            }
        }
        if let Some(value) = info.get("RPA")
            && !valid_ints(value)
        {
            return Err(invalid("unsupported malformed RPA"));
        }
        if let Some(value) = info.get("AS_SB_TABLE")
            && (value.split('|').count() != n + 1
                || value
                    .split('|')
                    .any(|v| v.split(',').count() != 2 || !valid_ints(v)))
        {
            return Err(invalid("unsupported malformed AS_SB_TABLE"));
        }
        let ad = ints(
            gt[tumor]
                .get("AD")
                .ok_or_else(|| invalid("missing tumour AD"))?,
        );
        let normal_ad = ints(
            gt[normal]
                .get("AD")
                .ok_or_else(|| invalid("missing normal AD"))?,
        );
        if ad.len() != alts.len() + 1 || normal_ad.len() != ad.len() {
            return Err(invalid("AD cardinality mismatch"));
        }
        let pos = fields[1].parse().map_err(|_| invalid("bad VCF position"))?;
        let phase = gt[tumor]
            .get("PGT")
            .zip(gt[tumor].get("PID"))
            .map(|(a, b)| format!("{a}{b}"));
        let mut r = Self {
            fields,
            info,
            gt,
            pos,
            alts,
            ad,
            normal_ad,
            lod: Vec::new(),
            orientation: 0.0,
            maf: 0.5,
            phase,
        };
        r.lod = r
            .values("TLOD")
            .iter()
            .map(|v| v * std::f64::consts::LN_10)
            .collect();
        if r.lod.len() != r.alts.len() {
            return Err(invalid("missing TLOD or cardinality mismatch"));
        }
        Ok(r)
    }
    pub fn values(&self, k: &str) -> Vec<f64> {
        self.info.get(k).map_or_else(Vec::new, |s| floats(s))
    }
    pub fn depth(&self) -> u32 {
        self.ad.iter().sum()
    }
    pub fn len(&self, i: usize) -> i32 {
        self.alts[i].len() as i32 - self.fields[3].len() as i32
    }
    pub fn max_lod(&self) -> usize {
        let mut i = 0;
        for n in 1..self.lod.len() {
            if self.lod[n] > self.lod[i] {
                i = n;
            }
        }
        i
    }
    pub fn datum(&self, i: usize, artifact: f64, non_somatic: f64) -> Datum {
        Datum {
            lod: self.lod[i],
            artifact,
            non_seq: 1.0 - (1.0 - artifact) * (1.0 - non_somatic),
            k: self.ad[i + 1],
            n: self.depth(),
            len: self.len(i),
        }
    }
}
