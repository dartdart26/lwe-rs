# lwe-rs

Rust implementation of Regev's LWE public key cryptosystem (Section 4 of [Regev's LWE survey](https://cims.nyu.edu/~regev/papers/lwesurvey.pdf)).

Rationale for writing it is to learn more about LWE and implement the barebones Z_q variant as close to "production" as rational for a toy project.

## What's here

- Keygen, encrypt, decrypt for single-bit plaintexts
- Homomorphic `add` and `sub`
- `trivial_encrypt` for mixing known plaintexts into homomorphic computations
- Constant-time on the secret key paths

## Try it

```bash
cargo run --example demo --release
cargo test
```

## Caveats

This is a reference / educational implementation, not a production crypto library.

- **Noise code could be insecure**
  - Needs a more careful look
- **No bootstrapping**
  - Noise grows with each homomorphic op; after enough ops, decryption fails
- **Constant-time is unverified**

## TODO

- **Noise tracking on `Ciphertext`** 
  - Carry an estimation so callers can see how close they are to decryption failure
- **Vetted parameter sets**
- **Verify constant-time**
- **Plaintext modulus**
  - `t > 2`** — encrypt elements of `Z_t` instead of just bits
- **Ring-LWE**
- **Property tests and benchmarks**
