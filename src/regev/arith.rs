//! Modular arithmetic primitives for `Z_q`.
use subtle::{Choice, ConditionallySelectable, ConstantTimeLess};

/// `Modulus` wraps a `u32` in `[2, 2^32)` together with a precomputed
/// Barrett constant `mu = floor(2^64 / q)`. The 32-bit cap keeps `a*b` in
/// `u64`, which lets `mul_mod` reduce via Barrett — integer multiplications
/// and a conditional subtract, no `%`, no data-dependent branch in source.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "u32", into = "u32")]
pub struct Modulus {
    q: u32,
    mu: u64,
}

impl Modulus {
    pub const fn new(q: u32) -> Self {
        assert!(q >= 2, "modulus must be at least 2");
        let mu = ((1u128 << 64) / q as u128) as u64;
        Self { q, mu }
    }

    pub const fn raw(self) -> u32 {
        self.q
    }
}

impl TryFrom<u32> for Modulus {
    type Error = &'static str;
    fn try_from(q: u32) -> Result<Self, Self::Error> {
        if q < 2 {
            return Err("modulus must be at least 2");
        }
        Ok(Self::new(q))
    }
}

impl From<Modulus> for u32 {
    fn from(m: Modulus) -> Self {
        m.q
    }
}

#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    zeroize::Zeroize,
    serde::Serialize,
    serde::Deserialize,
)]
pub struct ZqElement(u32);

impl ZqElement {
    pub const ZERO: Self = Self(0);

    /// Wraps a `u32` already in `[0, q)`.
    pub fn new(x: u32, q: Modulus) -> Self {
        assert!(x < q.raw());
        Self(x)
    }

    /// Reduces a signed integer into `[0, q)`.
    /// Not constant-time.
    pub fn from_i64(x: i64, q: Modulus) -> Self {
        Self(x.rem_euclid(q.raw() as i64) as u32)
    }

    pub fn raw(self) -> u32 {
        self.0
    }
}

impl ConditionallySelectable for ZqElement {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self(u32::conditional_select(&a.0, &b.0, choice))
    }
}

pub fn add_mod(a: ZqElement, b: ZqElement, q: Modulus) -> ZqElement {
    assert!(a.0 < q.raw() && b.0 < q.raw());
    let qq = q.raw() as u64;
    let sum = a.0 as u64 + b.0 as u64; // <= 2(q-1), fits in u64
    let reduced = sum.wrapping_sub(qq);
    let needs_reduce: Choice = !sum.ct_lt(&qq);
    ZqElement(u64::conditional_select(&sum, &reduced, needs_reduce) as u32)
}

pub fn sub_mod(a: ZqElement, b: ZqElement, q: Modulus) -> ZqElement {
    assert!(a.0 < q.raw() && b.0 < q.raw());
    let diff = a.0.wrapping_sub(b.0);
    let corrected = diff.wrapping_add(q.raw());
    let underflow: Choice = a.0.ct_lt(&b.0);
    ZqElement(u32::conditional_select(&diff, &corrected, underflow))
}

/// `(a * b) mod q` via Barrett reduction.
///
/// `quotient_est = ⌊product · mu / 2^64⌋` approximates `⌊product / q⌋`, where
/// `mu = ⌊2^64 / q⌋` was precomputed in `Modulus::new`. The estimate is always
/// within 1 of the true quotient, so `product - quotient_est * q` lands in
/// `[0, 2q)` — one conditional subtract finishes the reduction.
///
/// Avoids hardware division, which is variable-time on every modern CPU.
pub fn mul_mod(a: ZqElement, b: ZqElement, q: Modulus) -> ZqElement {
    assert!(a.0 < q.raw() && b.0 < q.raw());
    let qq = q.raw() as u64;

    let product = a.0 as u64 * b.0 as u64;
    let quotient_est = ((product as u128 * q.mu as u128) >> 64) as u64;
    let pre_reduced = product.wrapping_sub(quotient_est.wrapping_mul(qq)); // in [0, 2q)

    let reduced = pre_reduced.wrapping_sub(qq);
    let needs_reduce: Choice = !pre_reduced.ct_lt(&qq);
    ZqElement(u64::conditional_select(&pre_reduced, &reduced, needs_reduce) as u32)
}

pub fn inner_product_mod(a: &[ZqElement], b: &[ZqElement], q: Modulus) -> ZqElement {
    assert_eq!(a.len(), b.len(), "inner product length mismatch");
    let mut acc = ZqElement::ZERO;
    for (&x, &y) in a.iter().zip(b.iter()) {
        acc = add_mod(acc, mul_mod(x, y, q), q);
    }
    acc
}

pub fn add_assign_vec_mod(out: &mut [ZqElement], rhs: &[ZqElement], q: Modulus) {
    assert_eq!(out.len(), rhs.len(), "vector length mismatch");
    for (o, &r) in out.iter_mut().zip(rhs.iter()) {
        *o = add_mod(*o, r, q);
    }
}

/// Distance from `x` to `0` in `Z_q`, i.e. `min(x, q - x)`. Returns a value in `[0, ⌊q/2⌋]`.
pub fn distance_to_zero(x: ZqElement, q: Modulus) -> u32 {
    assert!(x.0 < q.raw());
    let other_way = q.raw() - x.0;
    let take_other_way: Choice = other_way.ct_lt(&x.0);
    u32::conditional_select(&x.0, &other_way, take_other_way)
}

#[cfg(test)]
mod tests {
    use super::*;

    const Q: Modulus = Modulus::new(17);

    fn z(x: u32) -> ZqElement {
        ZqElement::new(x, Q)
    }

    #[test]
    fn from_i64_handles_negatives() {
        assert_eq!(ZqElement::from_i64(0, Q).raw(), 0);
        assert_eq!(ZqElement::from_i64(5, Q).raw(), 5);
        assert_eq!(ZqElement::from_i64(17, Q).raw(), 0);
        assert_eq!(ZqElement::from_i64(18, Q).raw(), 1);
        assert_eq!(ZqElement::from_i64(-1, Q).raw(), 16);
        assert_eq!(ZqElement::from_i64(-17, Q).raw(), 0);
        assert_eq!(ZqElement::from_i64(-18, Q).raw(), 16);
    }

    #[test]
    fn add_sub_round_trip() {
        for a in 0..Q.raw() {
            for b in 0..Q.raw() {
                let s = add_mod(z(a), z(b), Q);
                assert!(s.raw() < Q.raw());
                assert_eq!(sub_mod(s, z(b), Q), z(a));
            }
        }
    }

    #[test]
    fn mul_matches_naive() {
        for a in 0..Q.raw() {
            for b in 0..Q.raw() {
                assert_eq!(mul_mod(z(a), z(b), Q).raw(), (a * b) % Q.raw());
            }
        }
    }

    #[test]
    fn mul_near_invariant_boundary() {
        // q just under 2^32. Exercises the Barrett path at the worst case where
        // (q-1)^2 is just under 2^64.
        let q = Modulus::new(u32::MAX - 4);
        let a = q.raw() - 1;
        let b = q.raw() - 1;
        let expected = ((a as u64 * b as u64) % q.raw() as u64) as u32;
        let result = mul_mod(ZqElement::new(a, q), ZqElement::new(b, q), q);
        assert_eq!(result.raw(), expected);
    }

    #[test]
    fn mul_mod_matches_naive_random() {
        use rand::rngs::StdRng;
        use rand::{RngExt, SeedableRng};
        const TEST_PRIMES: [u32; 8] = [
            7,
            17,
            257,
            1_009,
            65_537,
            1_000_003,
            2_147_483_647,
            4_294_967_291,
        ];
        let mut rng = StdRng::seed_from_u64(42);
        for &q_raw in &TEST_PRIMES {
            let q = Modulus::new(q_raw);
            for _ in 0..5000 {
                let a = rng.random_range(0..q_raw);
                let b = rng.random_range(0..q_raw);
                let got = mul_mod(ZqElement::new(a, q), ZqElement::new(b, q), q).raw();
                let want = ((a as u64 * b as u64) % q_raw as u64) as u32;
                assert_eq!(got, want, "q={q_raw} a={a} b={b}");
            }
        }
    }

    #[test]
    fn inner_product_simple() {
        let a = [z(1), z(2), z(3)];
        let b = [z(4), z(5), z(6)];
        // 1*4 + 2*5 + 3*6 = 32 = 17 + 15
        assert_eq!(inner_product_mod(&a, &b, Q).raw(), 15);
    }

    #[test]
    fn inner_product_empty_is_zero() {
        let a: [ZqElement; 0] = [];
        let b: [ZqElement; 0] = [];
        assert_eq!(inner_product_mod(&a, &b, Q), ZqElement::ZERO);
    }

    #[test]
    fn add_assign_vec_works() {
        let mut a = [z(1), z(2), z(3)];
        let b = [z(10), z(15), z(16)];
        add_assign_vec_mod(&mut a, &b, Q);
        assert_eq!(a, [z(11), z(0), z(2)]);
    }

    #[test]
    fn distance_to_zero_basics() {
        assert_eq!(distance_to_zero(z(0), Q), 0);
        assert_eq!(distance_to_zero(z(1), Q), 1);
        assert_eq!(distance_to_zero(z(8), Q), 8);
        assert_eq!(distance_to_zero(z(9), Q), 8);
        assert_eq!(distance_to_zero(z(16), Q), 1);
    }

    #[test]
    fn conditional_select_works() {
        let a = z(3);
        let b = z(7);
        assert_eq!(ZqElement::conditional_select(&a, &b, Choice::from(0)), a);
        assert_eq!(ZqElement::conditional_select(&a, &b, Choice::from(1)), b);
    }
}
