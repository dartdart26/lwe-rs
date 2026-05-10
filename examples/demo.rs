//! End-to-end LWE example: keygen, encrypt a bit, decrypt, print sizes.

use humansize::{BINARY, format_size};
use lwe::LweParams;
use lwe::lwe::{add, decrypt, encrypt, keygen, trivial_encrypt};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn main() {
    let params = LweParams::regev_default(2048).expect("regev_default(2048)");
    let mut rng = StdRng::seed_from_u64(0);

    let (sk, pk) = keygen(&params, &mut rng);

    println!(
        "n={}, q={}, m={}",
        params.dimension().0,
        params.modulus().raw(),
        params.sample_count().0,
    );
    println!(
        "sk={}, pk={}",
        format_size(sk.data_bytes(), BINARY),
        format_size(pk.data_bytes(), BINARY),
    );

    for bit in [true, false] {
        let ct = encrypt(&pk, bit, &mut rng);
        let recovered = decrypt(&sk, &ct);
        println!(
            "encrypted {bit} -> ct={} -> decrypted {recovered}",
            format_size(ct.data_bytes(), BINARY),
        );
        assert_eq!(bit, recovered);
    }

    println!();
    println!("homomorphic add (decrypts to XOR):");
    let q = params.modulus();
    for b1 in [false, true] {
        for b2 in [false, true] {
            let ct1 = encrypt(&pk, b1, &mut rng);
            let ct2 = encrypt(&pk, b2, &mut rng);
            let ct_sum = add(&ct1, &ct2, q);
            let recovered = decrypt(&sk, &ct_sum);
            let expected = b1 ^ b2;
            let mark = if recovered == expected { "✓" } else { "✗" };
            println!("  {mark} {b1} XOR {b2} = {recovered} (expected {expected})");
            assert_eq!(recovered, expected);
        }
    }

    println!();
    println!("trivial_encrypt decrypts to input:");
    for bit in [false, true] {
        let ct = trivial_encrypt(&params, bit);
        let recovered = decrypt(&sk, &ct);
        let mark = if recovered == bit { "✓" } else { "✗" };
        println!("  {mark} trivial({bit}) -> decrypted {recovered}");
        assert_eq!(recovered, bit);
    }

    println!();
    println!("XOR with plaintext bit (via trivial_encrypt + add):");
    for ct_bit in [false, true] {
        for pt_bit in [false, true] {
            let ct = encrypt(&pk, ct_bit, &mut rng);
            let trivial = trivial_encrypt(&params, pt_bit);
            let combined = add(&ct, &trivial, q);
            let recovered = decrypt(&sk, &combined);
            let expected = ct_bit ^ pt_bit;
            let mark = if recovered == expected { "✓" } else { "✗" };
            println!("  {mark} enc({ct_bit}) XOR {pt_bit} = {recovered} (expected {expected})");
            assert_eq!(recovered, expected);
        }
    }
}
