//! TFHE-style constructions following Zama's "TFHE Deep Dive" blog series
//! at https://www.zama.org/post/tfhe-deep-dive-part-1 .
//! GLWE and GLev ciphertexts over a power-of-two modulus, signed
//! decomposition, and multiplication by big constants.

pub mod modq;
pub mod poly;
