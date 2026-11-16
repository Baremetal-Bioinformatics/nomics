//! Salmon 2.7.0 library formats and deterministic alignment-mode detection.
use std::io;
use std::str::FromStr;

use crate::fragment::{Hit, Status};

/// Canonical Salmon formats, in the stable histogram order. `Auto` has no
/// observed-format histogram slot; it resolves from unique proper pairs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum LibraryType {
    Iu,
    Isf,
    Isr,
    Ou,
    Osf,
    Osr,
    Mu,
    Msf,
    Msr,
    U,
    Sf,
    Sr,
    Auto,
}

impl FromStr for LibraryType {
    type Err = io::Error;
    fn from_str(text: &str) -> io::Result<Self> {
        Self::ALL.into_iter().chain([Self::Auto]).find(|f| f.as_str().eq_ignore_ascii_case(text))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput,
                format!("invalid library type {text:?}; expected A, IU, ISF, ISR, OU, OSF, OSR, MU, MSF, MSR, U, SF, or SR")))
    }
}

impl LibraryType {
    pub const ALL: [Self; 12] = [
        Self::Iu,
        Self::Isf,
        Self::Isr,
        Self::Ou,
        Self::Osf,
        Self::Osr,
        Self::Mu,
        Self::Msf,
        Self::Msr,
        Self::U,
        Self::Sf,
        Self::Sr,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Iu => "IU",
            Self::Isf => "ISF",
            Self::Isr => "ISR",
            Self::Ou => "OU",
            Self::Osf => "OSF",
            Self::Osr => "OSR",
            Self::Mu => "MU",
            Self::Msf => "MSF",
            Self::Msr => "MSR",
            Self::U => "U",
            Self::Sf => "SF",
            Self::Sr => "SR",
            Self::Auto => "A",
        }
    }

    pub(crate) fn compatible(self, hit: &Hit) -> bool {
        debug_assert_ne!(self, Self::Auto);
        let group = self as usize / 3;
        let strand = self as usize % 3;
        if hit.status == Status::Pair {
            return group == hit.format() / 3 && (strand == 0 || self as usize == hit.format());
        }
        if strand == 0 {
            return true;
        }
        // Single-end and matching use S/A; opposed pairs use SA/AS. They
        // are distinct encodings even where the same F/R suffix is printed.
        match hit.status {
            Status::Single => (group >= 2) && (hit.forward == (strand == 1)),
            Status::Left => (group != 3) && (hit.forward == (strand == 1)),
            Status::Right => {
                (group != 3) && (hit.forward == (strand == if group == 2 { 1 } else { 2 }))
            }
            Status::Pair => unreachable!(),
        }
    }

    /// Concordant histogram count, inconsistent/orphan count, and sense share.
    pub(crate) fn summary(self, counts: &[u64; 12]) -> (u64, u64, f64) {
        if self == Self::Auto {
            return (0, 0, 0.0);
        }
        let id = self as usize;
        let base = id / 3 * 3;
        let (sense, antisense) = (counts[base + 1], counts[base + 2]);
        let agree = if id.is_multiple_of(3) {
            sense + antisense
        } else {
            counts[id]
        };
        let disagree = counts.iter().sum::<u64>() - agree;
        let bias = if agree > 0 && sense + antisense > 0 {
            sense as f64 / (sense + antisense) as f64
        } else {
            0.0
        };
        (agree, disagree, bias)
    }
}

/// The BAM writer bakes only the unique proper-pair tally. No 50k prefix
/// sampling, no single-end fallback: with no such evidence A stays unresolved.
pub(crate) fn detect(counts: &[u64; 12]) -> Option<LibraryType> {
    if counts.iter().sum::<u64>() == 0 {
        return None;
    }
    let inward: u64 = counts[..3].iter().sum();
    let outward: u64 = counts[3..6].iter().sum();
    let matching: u64 = counts[6..9].iter().sum();
    let total = (inward + outward + matching) as f64;
    let forward = counts[1] + counts[4] + counts[7] + counts[10];
    let reverse = counts[2] + counts[5] + counts[8] + counts[11];
    if total == 0.0 || forward + reverse == 0 {
        return Some(LibraryType::Iu);
    }
    let (ri, ro, rm) = (
        inward as f64 / total,
        outward as f64 / total,
        matching as f64 / total,
    );
    let base = if ri >= ro && ri >= rm {
        0
    } else if ro >= ri && ro >= rm {
        3
    } else {
        6
    };
    let ratio = forward as f64 / (forward + reverse) as f64;
    let strand = if ratio < 0.3 {
        2
    } else if ratio < 0.7 {
        0
    } else {
        1
    };
    Some(LibraryType::ALL[base + strand])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(status: Status, forward: bool, mate_forward: bool) -> Hit {
        Hit {
            tid: 0,
            pos: 0,
            forward,
            mate_forward,
            length: 100,
            score: 0,
            status,
        }
    }

    #[test]
    fn complete_compatibility_matrix_including_orphan_encodings() {
        // Columns: inward F/R, matching F/R, single F/R, left F/R, right F/R.
        let observations = [
            hit(Status::Pair, true, false),
            hit(Status::Pair, false, true),
            hit(Status::Pair, true, true),
            hit(Status::Pair, false, false),
            hit(Status::Single, true, false),
            hit(Status::Single, false, false),
            hit(Status::Left, true, false),
            hit(Status::Left, false, false),
            hit(Status::Right, true, false),
            hit(Status::Right, false, false),
        ];
        let expected = [
            "1100111111",
            "1000001001",
            "0100000110", // I
            "0000111111",
            "0000001001",
            "0000000110", // O: RAD can't observe outward pairs
            "0011111111",
            "0010101010",
            "0001010101", // M
            "0000111111",
            "0000100000",
            "0000010000", // single
        ];
        for (format, expected) in LibraryType::ALL.into_iter().zip(expected) {
            let actual: String = observations
                .iter()
                .map(|h| if format.compatible(h) { '1' } else { '0' })
                .collect();
            assert_eq!(actual, expected, "{}", format.as_str());
            assert_eq!(
                format
                    .as_str()
                    .to_ascii_lowercase()
                    .parse::<LibraryType>()
                    .unwrap(),
                format
            );
        }
        assert!("S".parse::<LibraryType>().is_err());
        assert!("XYZ".parse::<LibraryType>().is_err());
    }

    #[test]
    fn detection_thresholds_ties_and_no_evidence() {
        assert_eq!(detect(&[0; 12]), None);
        for (fwd, rev, result) in [
            (29, 71, LibraryType::Isr),
            (30, 70, LibraryType::Iu),
            (69, 31, LibraryType::Iu),
            (70, 30, LibraryType::Isf),
        ] {
            let mut counts = [0; 12];
            counts[1] = fwd;
            counts[2] = rev;
            assert_eq!(detect(&counts), Some(result));
        }
        let mut counts = [0; 12];
        counts[1] = 10;
        counts[7] = 10;
        assert_eq!(detect(&counts), Some(LibraryType::Isf));
        counts[7] = 11;
        assert_eq!(detect(&counts), Some(LibraryType::Msf));
    }

    #[test]
    fn summaries_count_each_observed_format_and_zero_bias_without_agreement() {
        let mut counts = [0; 12];
        counts[1] = 30;
        counts[2] = 20;
        counts[10] = 5;
        assert_eq!(LibraryType::Iu.summary(&counts), (50, 5, 0.6));
        assert_eq!(LibraryType::Isr.summary(&counts), (20, 35, 0.6));
        assert_eq!(LibraryType::Ou.summary(&counts), (0, 55, 0.0));
        counts[2] = 0;
        assert_eq!(LibraryType::Isr.summary(&counts), (0, 35, 0.0));
        assert_eq!(LibraryType::Auto.summary(&counts), (0, 0, 0.0));
    }
}
