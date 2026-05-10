pub mod arith;
pub mod lwe;
pub mod noise;
pub mod params;

pub use params::{Alpha, Dimension, LweParams, Modulus, ParamError, SampleCount};

pub trait SecureRng: rand::Rng + rand::CryptoRng {}
impl<T: rand::Rng + rand::CryptoRng> SecureRng for T {}
