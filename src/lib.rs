pub mod regev;
pub mod tfhe;

pub trait SecureRng: rand::Rng + rand::CryptoRng {}
impl<T: rand::Rng + rand::CryptoRng> SecureRng for T {}