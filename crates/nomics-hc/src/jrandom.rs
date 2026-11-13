//! `java.util.Random` exactly as OpenJDK 17 implements it.
//!
//! GATK's `Utils.getRandomGenerator()` is one `Random(47382911)` shared by
//! the positional downsampler (`nextInt`) and `QualByDepth`'s QD jitter
//! (`nextGaussian`), so both consumers must draw from one instance in
//! traversal order (spec §1.5, §8.4).

/// The seed of GATK's shared generator (`Utils.java`).
pub const GATK_SEED: i64 = 47382911;

const MULTIPLIER: i64 = 0x5DEECE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;
const DOUBLE_UNIT: f64 = 1.0 / (1u64 << 53) as f64; // 0x1.0p-53

#[derive(Debug, Clone)]
pub struct JavaRandom {
    seed: i64,
    next_next_gaussian: f64,
    have_next_next_gaussian: bool,
}

impl JavaRandom {
    /// `new Random(seed)`.
    pub fn new(seed: i64) -> Self {
        Self {
            seed: (seed ^ MULTIPLIER) & MASK,
            next_next_gaussian: 0.0,
            have_next_next_gaussian: false,
        }
    }

    /// GATK's shared generator, freshly created.
    pub fn gatk() -> Self {
        Self::new(GATK_SEED)
    }

    /// `next(bits)`.
    #[inline]
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND)) & MASK;
        ((self.seed as u64) >> (48 - bits)) as i32
    }

    /// `nextInt(bound)` (bound > 0).
    pub fn next_int(&mut self, bound: i32) -> i32 {
        assert!(bound > 0, "bound must be positive");
        let mut r = self.next(31);
        let m = bound - 1;
        if bound & m == 0 {
            ((bound as i64 * r as i64) >> 31) as i32
        } else {
            let mut u = r;
            loop {
                r = u % bound;
                // `u - r + m < 0` in Java int arithmetic (overflow intended).
                if u.wrapping_sub(r).wrapping_add(m) >= 0 {
                    break;
                }
                u = self.next(31);
            }
            r
        }
    }

    /// `nextDouble()`.
    pub fn next_double(&mut self) -> f64 {
        let hi = (self.next(26) as i64) << 27;
        let lo = self.next(27) as i64;
        (hi + lo) as f64 * DOUBLE_UNIT
    }

    /// `nextGaussian()`: Marsaglia polar method with a cached second value,
    /// `StrictMath.log` and `StrictMath.sqrt`.
    pub fn next_gaussian(&mut self) -> f64 {
        if self.have_next_next_gaussian {
            self.have_next_next_gaussian = false;
            return self.next_next_gaussian;
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            if s < 1.0 && s != 0.0 {
                let multiplier = (-2.0 * crate::jmath::log(s) / s).sqrt();
                self.next_next_gaussian = v2 * multiplier;
                self.have_next_next_gaussian = true;
                return v1 * multiplier;
            }
        }
    }
}
