//! Supported BWA-MEM scoring. `-B` alone does not rescale other options.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Options {
    mismatch_penalty: i32,
    pub(crate) mat: [i8; 25],
}
impl Default for Options {
    fn default() -> Self {
        Self::new(4).unwrap()
    }
}
impl Options {
    pub fn new(mismatch_penalty: i32) -> Result<Self, String> {
        if !matches!(mismatch_penalty, 3 | 4) {
            return Err(
                "BWA mismatch penalty must be 3 or 4; other values are not verified exact".into(),
            );
        }
        Ok(Self {
            mismatch_penalty,
            mat: matrix(mismatch_penalty),
        })
    }
    pub fn mismatch_penalty(self) -> i32 {
        self.mismatch_penalty
    }
    /// bwamem.c:522 and bwamem_pair.c:255, including gap-cost maxima.
    pub(crate) fn suboptimal_window(self) -> i32 {
        use crate::extend::{A, E_DEL, E_INS, O_DEL, O_INS};
        (A + self.mismatch_penalty)
            .max(O_DEL + E_DEL)
            .max(O_INS + E_INS)
    }
    /// One line in the run log when the penalty is not BWA's default. The GPU
    /// scoring kernels take the matrix as a parameter, so routing does not
    /// depend on the penalty.
    pub fn log_routing(self) {
        if self.mismatch_penalty != crate::extend::B {
            eprintln!(
                "BWA -B {}: scoring matrix selected for CPU and GPU kernels",
                self.mismatch_penalty
            );
        }
    }
}
impl std::str::FromStr for Options {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::new(
            s.parse()
                .map_err(|_| "BWA mismatch penalty must be 3 or 4".to_owned())?,
        )
    }
}
/// bwa.c:136-145: N always scores -1, including N/N.
pub(crate) const fn matrix(b: i32) -> [i8; 25] {
    let mut mat = [-1; 25];
    let mut i = 0;
    while i < 4 {
        let mut j = 0;
        while j < 4 {
            mat[i * 5 + j] = if i == j { 1 } else { -b as i8 };
            j += 1;
        }
        i += 1;
    }
    mat
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bwa_scoring_and_secondary_window() {
        for b in [3, 4] {
            let opt = Options::new(b).unwrap();
            for i in 0..5 {
                for j in 0..5 {
                    assert_eq!(
                        opt.mat[i * 5 + j],
                        if i == 4 || j == 4 {
                            -1
                        } else if i == j {
                            1
                        } else {
                            -b as i8
                        }
                    );
                }
            }
            assert_eq!(opt.suboptimal_window(), 7);
        }
        for b in [-1, 0, 1, 2, 5, 127, 128] {
            assert!(Options::new(b).is_err());
        }
    }
}
