//! TFHE-style constructions following Zama's "TFHE Deep Dive" blog series
//! (parts 1-3): GLWE and GLev ciphertexts over a power-of-two modulus, signed
//! decomposition, and multiplication by big constants.
//!
//! Independent of the [`crate::regev`] module: power-of-two modulus (not
//! prime), binary secret keys (not uniform), and no constant-time discipline —
//! this module optimizes for readability against the blog posts.

pub mod modq;