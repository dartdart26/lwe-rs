//! Parameter set for the Regev '05 LWE cryptosystem.
//!
//! See Regev's LWE survey, §4. A parameter choice that guarantees both
//! security and correctness is:
//!
//! ```text
//! modulus is a prime between dimension^2 and 2 * dimension^2
//! sample_count = 1.1 * dimension * log modulus
//! alpha = 1 / (sqrt(dimension) * (log dimension)^2)
//! ```
//!
//! The paper writes `sample_count` as a real-valued formula; we take the
//! ceiling in code.

use std::fmt;

pub use crate::regev::arith::Modulus;

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "LweParamsRaw")]
pub struct LweParams {
    dimension: Dimension,
    sample_count: SampleCount,
    modulus: Modulus,
    alpha: Alpha,
}

#[derive(serde::Deserialize)]
struct LweParamsRaw {
    dimension: Dimension,
    sample_count: SampleCount,
    modulus: Modulus,
    alpha: Alpha,
}

impl TryFrom<LweParamsRaw> for LweParams {
    type Error = ParamError;
    fn try_from(raw: LweParamsRaw) -> Result<Self, ParamError> {
        Self::new(raw.dimension, raw.sample_count, raw.modulus, raw.alpha)
    }
}

#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct Dimension(pub usize);

#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct SampleCount(pub usize);

/// Gaussian noise parameter. The error has standard deviation `alpha * modulus`.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct Alpha(pub f64);

#[derive(Debug)]
pub struct ParamError(&'static str);

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid LWE parameters: {}", self.0)
    }
}

impl std::error::Error for ParamError {}

impl From<&'static str> for ParamError {
    fn from(s: &'static str) -> Self {
        Self(s)
    }
}

impl LweParams {
    /// This function doesn't check for security, just for correct construction.
    /// It's up to callers to ensure input parameters are secure.
    pub fn new(
        dimension: Dimension,
        sample_count: SampleCount,
        modulus: Modulus,
        alpha: Alpha,
    ) -> Result<Self, ParamError> {
        if dimension.0 < 2 {
            return Err(ParamError("dimension must be at least 2"));
        }
        if sample_count.0 < dimension.0 {
            return Err(ParamError("sample_count must be at least dimension"));
        }
        let dim_sq = u32::try_from(dimension.0)
            .ok()
            .and_then(|d| d.checked_mul(d))
            .ok_or(ParamError("dimension^2 does not fit in u32"))?;
        let two_dim_sq = dim_sq
            .checked_mul(2)
            .ok_or(ParamError("2 * dimension^2 does not fit in u32"))?;
        let q = modulus.raw();
        if q <= dim_sq || q >= two_dim_sq {
            return Err(ParamError(
                "modulus must lie in (dimension^2, 2 * dimension^2)",
            ));
        }
        if !is_prime(q) {
            return Err(ParamError("modulus must be prime"));
        }
        if !alpha.0.is_finite() || alpha.0 <= 0.0 {
            return Err(ParamError("alpha must be finite and > 0"));
        }
        Ok(Self {
            dimension,
            sample_count,
            modulus,
            alpha,
        })
    }

    /// Parameter set from section 4 of the survey for a given dimension.
    pub fn regev_default(dimension: usize) -> Result<Self, ParamError> {
        let dim_u32 =
            u32::try_from(dimension).map_err(|_| ParamError("dimension does not fit in u32"))?;
        let dim_sq = dim_u32
            .checked_mul(dim_u32)
            .ok_or(ParamError("dimension^2 does not fit in u32"))?;
        let modulus = next_prime_above(dim_sq);
        let log_modulus = (modulus as f64).log2();
        let sample_count = (1.1 * dimension as f64 * log_modulus).ceil() as usize;
        let log_dimension = (dimension as f64).log2();
        let alpha = 1.0 / ((dimension as f64).sqrt() * log_dimension * log_dimension);
        Self::new(
            Dimension(dimension),
            SampleCount(sample_count),
            Modulus::try_from(modulus)?,
            Alpha(alpha),
        )
    }

    pub fn dimension(&self) -> Dimension {
        self.dimension
    }

    pub fn sample_count(&self) -> SampleCount {
        self.sample_count
    }

    pub fn modulus(&self) -> Modulus {
        self.modulus
    }

    pub fn alpha(&self) -> Alpha {
        self.alpha
    }

    /// Standard deviation of the per-sample error, `alpha * modulus`.
    pub fn noise_sigma(&self) -> f64 {
        self.alpha.0 * self.modulus.raw() as f64
    }
}

fn next_prime_above(lo: u32) -> u32 {
    let mut k = lo + 1;
    while !is_prime(k) {
        k += 1;
    }
    k
}

fn is_prime(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n.is_multiple_of(2) {
        return n == 2;
    }
    let mut d = 3u32;
    while d.saturating_mul(d) <= n {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 2;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regev_default_512_matches_paper() {
        let p = LweParams::regev_default(512).unwrap();
        assert_eq!(p.dimension().0, 512);
        assert_eq!(p.modulus().raw(), 262147);
        let dim_sq = (p.dimension().0 as u32).pow(2);
        assert!(p.modulus().raw() > dim_sq && p.modulus().raw() < 2 * dim_sq);
        assert_eq!(p.sample_count().0, 10138);
    }

    #[test]
    fn new_rejects_invalid() {
        // modulus composite (non-prime)
        assert!(
            LweParams::new(
                Dimension(16),
                SampleCount(100),
                Modulus::new(258),
                Alpha(0.01)
            )
            .is_err()
        );
        // modulus < dimension^2
        assert!(
            LweParams::new(
                Dimension(16),
                SampleCount(100),
                Modulus::new(251),
                Alpha(0.01)
            )
            .is_err()
        );
        // modulus >= 2 * dimension^2
        assert!(
            LweParams::new(
                Dimension(16),
                SampleCount(100),
                Modulus::new(521),
                Alpha(0.01)
            )
            .is_err()
        );
        // sample_count < dimension
        assert!(
            LweParams::new(
                Dimension(16),
                SampleCount(8),
                Modulus::new(257),
                Alpha(0.01)
            )
            .is_err()
        );
        // alpha = 0
        assert!(
            LweParams::new(
                Dimension(16),
                SampleCount(100),
                Modulus::new(257),
                Alpha(0.0)
            )
            .is_err()
        );
        // dimension too small
        assert!(
            LweParams::new(
                Dimension(1),
                SampleCount(100),
                Modulus::new(257),
                Alpha(0.01)
            )
            .is_err()
        );
    }

    #[test]
    fn primes_are_prime() {
        for &p in &[2u32, 3, 5, 7, 11, 65537, 262147] {
            assert!(is_prime(p), "{p} should be prime");
        }
        for &c in &[0u32, 1, 4, 9, 65536, 262144] {
            assert!(!is_prime(c), "{c} should be composite");
        }
    }
}
