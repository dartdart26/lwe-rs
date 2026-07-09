# lwe-rs

Rust implementations of two LWE-based cryptosystems, written to learn the schemes by building them. The two modules are independent of each other.

This is a reference / educational implementation, not a production crypto library. Parameter sets are not security-evaluated, and the noise code could be insecure and needs a more careful look.

## Regev (`src/regev/`)

Regev's LWE public key cryptosystem (Section 4 of [Regev's LWE survey](https://cims.nyu.edu/~regev/papers/lwesurvey.pdf)) — the barebones `Z_q` variant as close to "production" as rational for a toy project. Prime modulus, uniform-`Z_q` secrets.

- Keygen, encrypt, decrypt for single-bit plaintexts
- Homomorphic `add` and `sub`
- `trivial_encrypt` for mixing known plaintexts into homomorphic computations
- Constant-time on the secret key paths (unverified)

```bash
cargo run --example regev_demo --release
```

TODO:

- **Noise tracking on `Ciphertext`**
  - Carry an estimation so callers can see how close they are to decryption failure
- **Vetted parameter sets**
- **Verify constant-time**
- **Plaintext modulus `t > 2`**
  - Encrypt elements of `Z_t` instead of just bits
- **Property tests and benchmarks**

## TFHE (`src/tfhe/`)

The constructions from Zama's [TFHE Deep Dive](https://www.zama.org/post/tfhe-deep-dive-part-1) blog series. Power-of-two modulus, binary secrets.

- GLWE encrypt/decrypt over `Z_q[X]/(X^N + 1)` with `Δ = q/p` encoding
- Homomorphic `add`, `sub`, and small-constant multiplication
- GLev ciphertexts and big-constant multiplication via decomposition

```bash
cargo run --example tfhe_demo --release
```

TODO:

- **Rest of deep dive part 3**: key switching, GGSW, external product, CMux
- **Programmable bootstrapping** (part 4)
- **NTT/FFT polynomial multiplication** (currently schoolbook `O(N²)`)
