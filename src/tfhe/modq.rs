//! Arithmetic in `Z_q` for a power-of-two modulus `q = 2^bits`.

use crate::SecureRng;
use rand::RngExt;
use std::fmt;

/// `q = 2^bits` with `bits` in `[2, 64]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "u32", into = "u32")]
pub struct Modulus {
    bits: u32,
    mask: u64,
}

#[derive(Debug)]
pub struct ModulusError(&'static str);

impl fmt::Display for ModulusError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid modulus: {}", self.0)
    }
}

impl std::error::Error for ModulusError {}

impl Modulus {
    pub const fn new(bits: u32) -> Result<Self, ModulusError> {
        if bits < 2 || bits > 64 {
            return Err(ModulusError("bits must be in [2, 64]"));
        }
        Ok(Self {
            bits,
            mask: u64::MAX >> (64 - bits),
        })
    }

    pub const fn bits(self) -> u32 {
        self.bits
    }

    pub const fn q(self) -> u128 {
        1u128 << self.bits
    }
}

pub trait Modulo {
    /// `self mod q`.
    fn modulo(self, q: Modulus) -> u64;
}

impl Modulo for u64 {
    fn modulo(self, q: Modulus) -> u64 {
        self & q.mask
    }
}

impl TryFrom<u32> for Modulus {
    type Error = ModulusError;
    fn try_from(bits: u32) -> Result<Self, Self::Error> {
        Self::new(bits)
    }
}

impl From<Modulus> for u32 {
    fn from(m: Modulus) -> Self {
        m.bits
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ZqElement(u64);

impl ZqElement {
    pub const ZERO: Self = Self(0);

    pub fn from_u64(x: u64, q: Modulus) -> Self {
        Self(x.modulo(q))
    }

    pub fn from_i64(x: i64, q: Modulus) -> Self {
        Self((x as u64).modulo(q))
    }

    pub fn raw(self) -> u64 {
        self.0
    }

    /// Centered lift: the representative of this element in `[-q/2, q/2)`.
    pub fn to_centered(self, q: Modulus) -> i64 {
        let x = self.0 as i128;
        let q = q.q() as i128;
        if x < q / 2 { x as i64 } else { (x - q) as i64 }
    }

    pub fn random<R: SecureRng>(rng: &mut R, q: Modulus) -> Self {
        Self(rng.random::<u64>().modulo(q))
    }
}

pub fn add_mod(a: ZqElement, b: ZqElement, q: Modulus) -> ZqElement {
    ZqElement(a.0.wrapping_add(b.0).modulo(q))
}

pub fn sub_mod(a: ZqElement, b: ZqElement, q: Modulus) -> ZqElement {
    ZqElement(a.0.wrapping_sub(b.0).modulo(q))
}

pub fn neg_mod(a: ZqElement, q: Modulus) -> ZqElement {
    ZqElement(a.0.wrapping_neg().modulo(q))
}

pub fn mul_mod(a: ZqElement, b: ZqElement, q: Modulus) -> ZqElement {
    ZqElement(a.0.wrapping_mul(b.0).modulo(q))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    const Q6: Modulus = match Modulus::new(6) {
        Ok(q) => q,
        Err(_) => panic!("invalid modulus"),
    };
    const Q64: Modulus = match Modulus::new(64) {
        Ok(q) => q,
        Err(_) => panic!("invalid modulus"),
    };

    fn z6(x: u64) -> ZqElement {
        ZqElement::from_u64(x, Q6)
    }

    #[test]
    fn modulus_bounds() {
        assert!(Modulus::new(1).is_err());
        assert!(Modulus::new(65).is_err());
        assert!(Modulus::new(64).is_ok());
        assert_eq!(Q6.q(), 64);
        assert_eq!(Q64.q(), 1u128 << 64);
    }

    #[test]
    fn modulo_matches_rem() {
        for x in [0u64, 1, 63, 64, 65, 200, u64::MAX] {
            assert_eq!(x.modulo(Q6), x % 64);
        }
        assert_eq!(u64::MAX.modulo(Q64), u64::MAX);
    }

    #[test]
    fn ops_match_naive_small_q() {
        for a in 0..64 {
            for b in 0..64 {
                assert_eq!(add_mod(z6(a), z6(b), Q6).raw(), (a + b) % 64);
                assert_eq!(sub_mod(z6(a), z6(b), Q6).raw(), (64 + a - b) % 64);
                assert_eq!(mul_mod(z6(a), z6(b), Q6).raw(), (a * b) % 64);
            }
        }
    }

    #[test]
    fn ops_wrap_at_native_width() {
        let a = ZqElement::from_u64(u64::MAX, Q64);
        let b = ZqElement::from_u64(1, Q64);
        assert_eq!(add_mod(a, b, Q64), ZqElement::ZERO);
        assert_eq!(sub_mod(ZqElement::ZERO, b, Q64).raw(), u64::MAX);
        let c = ZqElement::from_u64(1u64 << 63, Q64);
        assert_eq!(mul_mod(c, ZqElement::from_u64(2, Q64), Q64), ZqElement::ZERO);
    }

    #[test]
    fn from_i64_handles_negatives() {
        assert_eq!(ZqElement::from_i64(-1, Q6).raw(), 63);
        assert_eq!(ZqElement::from_i64(-64, Q6).raw(), 0);
        assert_eq!(ZqElement::from_i64(5, Q6).raw(), 5);
        assert_eq!(ZqElement::from_i64(-1, Q64).raw(), u64::MAX);
    }

    #[test]
    fn to_centered_basics() {
        assert_eq!(z6(0).to_centered(Q6), 0);
        assert_eq!(z6(31).to_centered(Q6), 31);
        assert_eq!(z6(32).to_centered(Q6), -32);
        assert_eq!(z6(63).to_centered(Q6), -1);
        assert_eq!(ZqElement::from_u64(u64::MAX, Q64).to_centered(Q64), -1);
        assert_eq!(ZqElement::from_u64(1u64 << 63, Q64).to_centered(Q64), i64::MIN);
    }

    #[test]
    fn centered_round_trip() {
        for x in 0..64 {
            let e = z6(x);
            assert_eq!(ZqElement::from_i64(e.to_centered(Q6), Q6), e);
        }
    }

    #[test]
    fn random_is_in_range() {
        let mut rng = StdRng::seed_from_u64(0);
        for _ in 0..1000 {
            assert!(ZqElement::random(&mut rng, Q6).raw() < 64);
        }
    }
}
