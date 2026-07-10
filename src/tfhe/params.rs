//! GLWE parameter sets.

use crate::tfhe::modq::Modulus;

/// Parameters for GLWE over `Z_q[X]/(X^N + 1)`.
///
/// Plaintexts live in `Z_p` with `p = 2^plaintext_bits`, encoded into the top
/// bits of `Z_q` by scaling with `Δ = q/p`.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(from = "GlweParamsRaw", into = "GlweParamsRaw")]
pub struct GlweParams {
    k: usize,
    n: usize,
    q: Modulus,
    plaintext_bits: u32,
    sigma: f64,
}

impl GlweParams {
    pub fn new(k: usize, n: usize, q: Modulus, plaintext_bits: u32, sigma: f64) -> Self {
        assert!(k >= 1);
        assert!(n.is_power_of_two());
        assert!(plaintext_bits >= 1 && plaintext_bits < q.bits());
        assert!(sigma.is_finite() && sigma >= 0.0);
        Self {
            k,
            n,
            q,
            plaintext_bits,
            sigma,
        }
    }

    /// Number of secret polynomials.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Ring degree `N`.
    pub fn n(&self) -> usize {
        self.n
    }

    pub fn q(&self) -> Modulus {
        self.q
    }

    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Plaintext modulus `p = 2^plaintext_bits`.
    pub fn p(&self) -> u64 {
        1 << self.plaintext_bits
    }

    /// `Δ = q/p`, the scaling factor between plaintext and ciphertext space.
    pub fn delta(&self) -> u64 {
        1 << (self.q.bits() - self.plaintext_bits)
    }
}

/// Toy size: `q = 2^6`, `p = 2^2`, `N = 4`, `k = 2`, no noise, so results can be traced by hand.
pub fn toy() -> GlweParams {
    GlweParams::new(2, 4, Modulus::new(6), 2, 0.0)
}

/// Realistically shaped and, most probably, insecure parameters.
pub fn test_default() -> GlweParams {
    GlweParams::new(2, 256, Modulus::new(64), 8, (1u64 << 36) as f64)
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
struct GlweParamsRaw {
    k: usize,
    n: usize,
    q: Modulus,
    plaintext_bits: u32,
    sigma: f64,
}

impl From<GlweParamsRaw> for GlweParams {
    fn from(r: GlweParamsRaw) -> Self {
        Self::new(r.k, r.n, r.q, r.plaintext_bits, r.sigma)
    }
}

impl From<GlweParams> for GlweParamsRaw {
    fn from(p: GlweParams) -> Self {
        Self {
            k: p.k,
            n: p.n,
            q: p.q,
            plaintext_bits: p.plaintext_bits,
            sigma: p.sigma,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn n_must_be_power_of_two() {
        GlweParams::new(2, 3, Modulus::new(6), 2, 0.0);
    }

    #[test]
    #[should_panic]
    fn plaintext_must_be_smaller_than_q() {
        GlweParams::new(2, 4, Modulus::new(6), 6, 0.0);
    }

    #[test]
    fn serde_round_trip() {
        let params = test_default();
        let bytes = bincode::serialize(&params).unwrap();
        assert_eq!(bincode::deserialize::<GlweParams>(&bytes).unwrap(), params);
    }

    #[test]
    #[should_panic]
    fn deserializing_non_power_of_two_n_panics() {
        let raw = (2usize, 3usize, 6u32, 2u32, 0.0f64);
        let bytes = bincode::serialize(&raw).unwrap();
        let _ = bincode::deserialize::<GlweParams>(&bytes);
    }
}
