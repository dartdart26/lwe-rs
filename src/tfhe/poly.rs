//! Polynomials in `R_q = Z_q[X]/(X^N + 1)`.

use crate::SecureRng;
use crate::tfhe::modq::{Modulus, ZqElement, add_mod, mul_mod, neg_mod, sub_mod};

/// A polynomial with `N >= 1` coefficients in `Z_q`; index `i` holds the
/// coefficient of `X^i`.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(from = "Vec<ZqElement>", into = "Vec<ZqElement>")]
pub struct Poly {
    coeffs: Vec<ZqElement>,
}

impl Poly {
    pub fn new(coeffs: Vec<ZqElement>) -> Self {
        assert!(!coeffs.is_empty());
        Self { coeffs }
    }

    pub fn zero(n: usize) -> Self {
        Self::new(vec![ZqElement::ZERO; n])
    }

    pub fn from_signed(coeffs: &[i64], q: Modulus) -> Self {
        Self::new(coeffs.iter().map(|&x| ZqElement::from_i64(x, q)).collect())
    }

    pub fn random<R: SecureRng>(rng: &mut R, n: usize, q: Modulus) -> Self {
        Self::new((0..n).map(|_| ZqElement::random(rng, q)).collect())
    }

    /// The number of coefficients.
    pub fn n(&self) -> usize {
        self.coeffs.len()
    }

    pub fn coeff(&self, i: usize) -> ZqElement {
        self.coeffs[i]
    }

    pub fn coeffs(&self) -> &[ZqElement] {
        &self.coeffs
    }

    pub fn to_centered(&self, q: Modulus) -> Vec<i64> {
        self.coeffs.iter().map(|c| c.to_centered(q)).collect()
    }

    pub fn add(&self, other: &Self, q: Modulus) -> Self {
        assert_eq!(self.n(), other.n());
        Self::new(
            self.coeffs
                .iter()
                .zip(&other.coeffs)
                .map(|(&a, &b)| add_mod(a, b, q))
                .collect(),
        )
    }

    pub fn add_assign(&mut self, other: &Self, q: Modulus) {
        assert_eq!(self.n(), other.n());
        for (a, &b) in self.coeffs.iter_mut().zip(&other.coeffs) {
            *a = add_mod(*a, b, q);
        }
    }

    pub fn sub(&self, other: &Self, q: Modulus) -> Self {
        assert_eq!(self.n(), other.n());
        Self::new(
            self.coeffs
                .iter()
                .zip(&other.coeffs)
                .map(|(&a, &b)| sub_mod(a, b, q))
                .collect(),
        )
    }

    pub fn neg(&self, q: Modulus) -> Self {
        Self::new(self.coeffs.iter().map(|&a| neg_mod(a, q)).collect())
    }

    pub fn scalar_mul(&self, c: ZqElement, q: Modulus) -> Self {
        Self::new(self.coeffs.iter().map(|&a| mul_mod(a, c, q)).collect())
    }

    /// Schoolbook product in `Z_q[X]/(X^N + 1)`: terms that overflow degree
    /// `N - 1` wrap around negated, since `X^N = -1`.
    pub fn mul(&self, other: &Self, q: Modulus) -> Self {
        assert_eq!(self.n(), other.n());
        let n = self.n();
        let mut out = Self::zero(n);
        for i in 0..n {
            for j in 0..n {
                let p = mul_mod(self.coeffs[i], other.coeffs[j], q);
                if i + j < n {
                    out.coeffs[i + j] = add_mod(out.coeffs[i + j], p, q);
                } else {
                    out.coeffs[i + j - n] = sub_mod(out.coeffs[i + j - n], p, q);
                }
            }
        }
        out
    }
}

impl From<Vec<ZqElement>> for Poly {
    fn from(coeffs: Vec<ZqElement>) -> Self {
        Self::new(coeffs)
    }
}

impl From<Poly> for Vec<ZqElement> {
    fn from(p: Poly) -> Self {
        p.coeffs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    const Q6: Modulus = Modulus::new(6);

    fn p(coeffs: &[i64]) -> Poly {
        Poly::from_signed(coeffs, Q6)
    }

    #[test]
    #[should_panic]
    fn empty_poly_panics() {
        Poly::new(vec![]);
    }

    #[test]
    #[should_panic]
    fn mismatched_lengths_panic() {
        p(&[1, 2]).add(&p(&[1, 2, 3]), Q6);
    }

    #[test]
    fn signed_round_trip() {
        let a = p(&[28, -5, -30, 17]);
        assert_eq!(a.to_centered(Q6), vec![28, -5, -30, 17]);
    }

    #[test]
    fn add_sub_neg_match_naive() {
        let a = p(&[1, 2, 3, 4]);
        let b = p(&[5, -6, 7, -8]);
        assert_eq!(a.add(&b, Q6), p(&[6, -4, 10, -4]));
        assert_eq!(a.sub(&b, Q6), p(&[-4, 8, -4, 12]));
        assert_eq!(b.neg(Q6), p(&[-5, 6, -7, 8]));
    }

    #[test]
    fn add_assign_matches_add() {
        let a = p(&[1, 2, 3, 4]);
        let b = p(&[5, -6, 7, -8]);
        let mut c = a.clone();
        c.add_assign(&b, Q6);
        assert_eq!(c, a.add(&b, Q6));
    }

    #[test]
    fn scalar_mul_matches_naive() {
        let a = p(&[1, -2, 0, 31]);
        let c = ZqElement::from_i64(3, Q6);
        assert_eq!(a.scalar_mul(c, Q6), p(&[3, -6, 0, 93]));
    }

    #[test]
    fn x_times_x_cubed_wraps_to_minus_one() {
        let x = p(&[0, 1, 0, 0]);
        let x3 = p(&[0, 0, 0, 1]);
        assert_eq!(x.mul(&x3, Q6), p(&[-1, 0, 0, 0]));
    }

    #[test]
    fn mul_matches_naive() {
        let mut rng = StdRng::seed_from_u64(1);
        let n = 8;
        for _ in 0..20 {
            let a = Poly::random(&mut rng, n, Q6);
            let b = Poly::random(&mut rng, n, Q6);

            let mut full = vec![0i64; 2 * n - 1];
            for i in 0..n {
                for j in 0..n {
                    full[i + j] += (a.coeff(i).raw() * b.coeff(j).raw()) as i64;
                }
            }
            let mut folded = vec![0i64; n];
            for (k, &c) in full.iter().enumerate() {
                if k < n {
                    folded[k] += c;
                } else {
                    folded[k - n] -= c;
                }
            }

            assert_eq!(a.mul(&b, Q6), Poly::from_signed(&folded, Q6));
            assert_eq!(a.mul(&b, Q6), b.mul(&a, Q6));
        }
    }

    #[test]
    fn serde_round_trip() {
        let a = p(&[28, -5, -30, 17]);
        let bytes = bincode::serialize(&a).unwrap();
        assert_eq!(bincode::deserialize::<Poly>(&bytes).unwrap(), a);
    }

    #[test]
    #[should_panic]
    fn deserializing_empty_panics() {
        let bytes = bincode::serialize(&Vec::<ZqElement>::new()).unwrap();
        let _ = bincode::deserialize::<Poly>(&bytes);
    }
}