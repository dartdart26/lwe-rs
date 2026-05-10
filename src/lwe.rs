use crate::SecureRng;
use crate::arith::{Modulus, ZqElement, add_mod, distance_to_zero, inner_product_mod, sub_mod};
use crate::noise::RoundedGaussian;
use crate::params::LweParams;
use rand::RngExt;
use subtle::{Choice, ConditionallySelectable};

#[derive(
    Clone, Debug, serde::Serialize, serde::Deserialize, zeroize::Zeroize, zeroize::ZeroizeOnDrop,
)]
pub struct SecretKey {
    #[zeroize(skip)]
    params: LweParams,
    s: Vec<ZqElement>,
}

#[derive(Clone, Debug, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct LweSample {
    a: Vec<ZqElement>,
    b: ZqElement,
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct PublicKey {
    params: LweParams,
    samples: Vec<LweSample>,
}

pub type Ciphertext = LweSample;

impl SecretKey {
    pub fn params(&self) -> &LweParams {
        &self.params
    }

    /// Bytes in the secret data, not counting Vec/struct overhead or params.
    /// Use a serde format (e.g. bincode) for the actual on-the-wire size.
    pub fn data_bytes(&self) -> usize {
        self.s.len() * std::mem::size_of::<ZqElement>()
    }
}

impl PublicKey {
    pub fn params(&self) -> &LweParams {
        &self.params
    }

    /// Bytes in the LWE samples, not counting Vec/struct overhead or params.
    pub fn data_bytes(&self) -> usize {
        self.samples
            .iter()
            .map(|s| (s.a.len() + 1) * std::mem::size_of::<ZqElement>())
            .sum()
    }
}

impl Ciphertext {
    /// Bytes in the ciphertext data, not counting Vec/struct overhead.
    pub fn data_bytes(&self) -> usize {
        (self.a.len() + 1) * std::mem::size_of::<ZqElement>()
    }
}

fn uniform_zq<R: SecureRng>(rng: &mut R, q: Modulus) -> ZqElement {
    ZqElement::new(rng.random_range(0..q.raw()), q)
}

fn uniform_zq_vec<R: SecureRng>(rng: &mut R, n: usize, q: Modulus) -> Vec<ZqElement> {
    (0..n).map(|_| uniform_zq(rng, q)).collect()
}

pub fn keygen<R: SecureRng>(params: &LweParams, rng: &mut R) -> (SecretKey, PublicKey) {
    let n = params.dimension().0;
    let m = params.sample_count().0;
    let q = params.modulus();
    let noise = RoundedGaussian::for_lwe_params(params);

    let s = uniform_zq_vec(rng, n, q);

    let samples: Vec<LweSample> = (0..m)
        .map(|_| {
            let a = uniform_zq_vec(rng, n, q);
            let b = add_mod(inner_product_mod(&a, &s, q), noise.sample(rng), q);
            LweSample { a, b }
        })
        .collect();

    (
        SecretKey {
            params: params.clone(),
            s,
        },
        PublicKey {
            params: params.clone(),
            samples,
        },
    )
}

pub fn encrypt<R: SecureRng>(pk: &PublicKey, bit: bool, rng: &mut R) -> Ciphertext {
    let n = pk.params.dimension().0;
    let q = pk.params.modulus();

    let mut a = vec![ZqElement::ZERO; n];
    let mut b = ZqElement::ZERO;

    for sample in &pk.samples {
        let take: Choice = (rng.random::<bool>() as u8).into();
        for (a_i, &sample_a_i) in a.iter_mut().zip(sample.a.iter()) {
            let masked_a_i = ZqElement::conditional_select(&ZqElement::ZERO, &sample_a_i, take);
            *a_i = add_mod(*a_i, masked_a_i, q);
        }
        let masked_b = ZqElement::conditional_select(&ZqElement::ZERO, &sample.b, take);
        b = add_mod(b, masked_b, q);
    }

    // Encode the plaintext bit into `b`: add ⌊q/2⌋ when bit is 1, else 0.
    // `q.raw() / 2` is exactly ⌊q/2⌋.
    let half = ZqElement::new(q.raw() / 2, q);
    let bit_choice: Choice = (bit as u8).into();
    let to_add = ZqElement::conditional_select(&ZqElement::ZERO, &half, bit_choice);
    b = add_mod(b, to_add, q);

    Ciphertext { a, b }
}

pub fn decrypt(sk: &SecretKey, ct: &Ciphertext) -> bool {
    let q = sk.params.modulus();
    let inner = inner_product_mod(&ct.a, &sk.s, q);
    let diff = sub_mod(ct.b, inner, q);
    // `diff` lands near `0` (bit 0) or near `⌊q/2⌋` (bit 1). `⌊q/4⌋` is the
    // midpoint between those clusters: distance > ⌊q/4⌋ ⇒ closer to ⌊q/2⌋ ⇒ bit 1.
    distance_to_zero(diff, q) > q.raw() / 4
}

/// "Encrypts" a known plaintext bit as `(0, bit * ⌊q/2⌋).
///
/// Not secret as anyone can decrypt without the secret key.
pub fn trivial_encrypt(params: &LweParams, bit: bool) -> Ciphertext {
    let n = params.dimension().0;
    let q = params.modulus();
    let b = if bit {
        ZqElement::new(q.raw() / 2, q)
    } else {
        ZqElement::ZERO
    };
    Ciphertext {
        a: vec![ZqElement::ZERO; n],
        b,
    }
}

pub fn add(ct1: &Ciphertext, ct2: &Ciphertext, q: Modulus) -> Ciphertext {
    assert_eq!(ct1.a.len(), ct2.a.len(), "ciphertext dimension mismatch");
    let a = ct1
        .a
        .iter()
        .zip(ct2.a.iter())
        .map(|(&x, &y)| add_mod(x, y, q))
        .collect();
    let b = add_mod(ct1.b, ct2.b, q);
    Ciphertext { a, b }
}

pub fn sub(ct1: &Ciphertext, ct2: &Ciphertext, q: Modulus) -> Ciphertext {
    assert_eq!(ct1.a.len(), ct2.a.len(), "ciphertext dimension mismatch");
    let a = ct1
        .a
        .iter()
        .zip(ct2.a.iter())
        .map(|(&x, &y)| sub_mod(x, y, q))
        .collect();
    let b = sub_mod(ct1.b, ct2.b, q);
    Ciphertext { a, b }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn encrypt_decrypt_roundtrip() {
        let params = LweParams::regev_default(256).unwrap();
        let mut rng = StdRng::seed_from_u64(0);
        let (sk, pk) = keygen(&params, &mut rng);

        for bit in [true, false] {
            let ct = encrypt(&pk, bit, &mut rng);
            assert_eq!(decrypt(&sk, &ct), bit);
        }
    }

    #[test]
    fn encryptions_of_same_bit_differ() {
        let params = LweParams::regev_default(256).unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        let (_sk, pk) = keygen(&params, &mut rng);

        let ct1 = encrypt(&pk, true, &mut rng);
        let ct2 = encrypt(&pk, true, &mut rng);

        assert_ne!(ct1, ct2);
    }

    #[test]
    fn homomorphic_add_is_xor() {
        let params = LweParams::regev_default(256).unwrap();
        let q = params.modulus();
        let mut rng = StdRng::seed_from_u64(3);
        let (sk, pk) = keygen(&params, &mut rng);

        for b1 in [false, true] {
            for b2 in [false, true] {
                let ct1 = encrypt(&pk, b1, &mut rng);
                let ct2 = encrypt(&pk, b2, &mut rng);
                let ct_sum = add(&ct1, &ct2, q);
                assert_eq!(decrypt(&sk, &ct_sum), b1 ^ b2, "{b1} XOR {b2}");
            }
        }
    }

    #[test]
    fn homomorphic_sub_is_xor() {
        let params = LweParams::regev_default(256).unwrap();
        let q = params.modulus();
        let mut rng = StdRng::seed_from_u64(6);
        let (sk, pk) = keygen(&params, &mut rng);

        for b1 in [false, true] {
            for b2 in [false, true] {
                let ct1 = encrypt(&pk, b1, &mut rng);
                let ct2 = encrypt(&pk, b2, &mut rng);
                let ct_diff = sub(&ct1, &ct2, q);
                assert_eq!(decrypt(&sk, &ct_diff), b1 ^ b2, "{b1} - {b2}");
            }
        }
    }

    #[test]
    fn trivial_encrypt_decrypts_to_input() {
        let params = LweParams::regev_default(256).unwrap();
        let mut rng = StdRng::seed_from_u64(5);
        let (sk, _pk) = keygen(&params, &mut rng);

        for bit in [false, true] {
            let ct = trivial_encrypt(&params, bit);
            assert_eq!(decrypt(&sk, &ct), bit);
        }
    }

    #[test]
    fn add_with_trivial_is_plaintext_xor() {
        let params = LweParams::regev_default(256).unwrap();
        let q = params.modulus();
        let mut rng = StdRng::seed_from_u64(4);
        let (sk, pk) = keygen(&params, &mut rng);

        for ct_bit in [false, true] {
            for pt_bit in [false, true] {
                let ct = encrypt(&pk, ct_bit, &mut rng);
                let trivial = trivial_encrypt(&params, pt_bit);
                let combined = add(&ct, &trivial, q);
                assert_eq!(
                    decrypt(&sk, &combined),
                    ct_bit ^ pt_bit,
                    "{ct_bit} XOR {pt_bit}"
                );
            }
        }
    }

    #[test]
    fn serde_roundtrip() {
        let params = LweParams::regev_default(256).unwrap();
        let mut rng = StdRng::seed_from_u64(2);
        let (sk, pk) = keygen(&params, &mut rng);
        let ct = encrypt(&pk, true, &mut rng);

        let pk_bytes = bincode::serialize(&pk).unwrap();
        let pk2: PublicKey = bincode::deserialize(&pk_bytes).unwrap();

        let ct_bytes = bincode::serialize(&ct).unwrap();
        let ct2: Ciphertext = bincode::deserialize(&ct_bytes).unwrap();

        assert_eq!(ct, ct2);
        assert!(decrypt(&sk, &ct2));

        let ct3 = encrypt(&pk2, false, &mut rng);
        assert!(!decrypt(&sk, &ct3));
    }
}
