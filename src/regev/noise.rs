//! Rounded Gaussian error sampler.
//!
//! Samples are drawn from a centered continuous Gaussian with the given
//! standard deviation, rounded to the nearest integer, and reduced mod `q`.
//! This approximates the true discrete Gaussian.

use crate::SecureRng;
use crate::regev::arith::{Modulus, ZqElement};
use crate::regev::params::LweParams;
use rand_distr::{Distribution, Normal};

pub struct RoundedGaussian {
    normal: Normal<f64>,
    modulus: Modulus,
}

impl RoundedGaussian {
    pub fn new(sigma: f64, modulus: Modulus) -> Self {
        let normal = Normal::new(0.0, sigma).expect("creating normal distribution failed");
        Self { normal, modulus }
    }

    pub fn for_lwe_params(params: &LweParams) -> Self {
        Self::new(params.noise_sigma(), params.modulus())
    }

    /// Gaussian sample, rounded to i64. NaN/inf are caught by `abs < q`.
    pub(crate) fn sample_i64<R: SecureRng>(&self, rng: &mut R) -> i64 {
        let rounded = self.normal.sample(rng).round();
        let q = self.modulus.raw() as f64;
        assert!(
            rounded.abs() < q,
            "Gaussian sample {rounded} exceeds modulus {q}; sigma is too large for this q"
        );
        rounded as i64
    }

    pub fn sample<R: SecureRng>(&self, rng: &mut R) -> ZqElement {
        ZqElement::from_i64(self.sample_i64(rng), self.modulus)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn samples_are_in_range() {
        let q = Modulus::new(101);
        let s = RoundedGaussian::new(2.5, q);
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..1000 {
            assert!(s.sample(&mut rng).raw() < q.raw());
        }
    }

    #[test]
    fn sample_mean_near_zero() {
        let s = RoundedGaussian::new(50.0, Modulus::new(10_007));
        let mut rng = StdRng::seed_from_u64(7);
        let n = 10_000;
        let sum: i64 = (0..n).map(|_| s.sample_i64(&mut rng)).sum();
        let mean = sum as f64 / n as f64;
        assert!(mean.abs() < 3.0, "sample mean {mean} not near 0");
    }

    #[test]
    fn sample_std_dev_near_target() {
        let sigma = 100.0;
        let s = RoundedGaussian::new(sigma, Modulus::new(1_000_003));
        let mut rng = StdRng::seed_from_u64(99);
        let n = 10_000;
        let sum_sq: f64 = (0..n)
            .map(|_| {
                let e = s.sample_i64(&mut rng) as f64;
                e * e
            })
            .sum();
        let observed = (sum_sq / n as f64).sqrt();
        assert!(
            (observed - sigma).abs() / sigma < 0.05,
            "observed sigma {observed} far from target {sigma}"
        );
    }
}
