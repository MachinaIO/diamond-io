use std::fmt::Debug;
pub trait PolyParams: Clone + Debug + PartialEq + Eq + Send + Sync {
    type Modulus: Debug + Clone;
    /// Returns the modulus value `q` used for polynomial coefficients in the ring `Z_q[x]/(x^n -
    /// 1)`.
    fn modulus(&self) -> Self::Modulus;
    /// A size of the base value used for a gadget vector and decomposition, i.e., `base =
    /// 2^base_bits`.
    fn base_bits(&self) -> u32;
    /// Fewest bits necessary to represent the modulus value `q`.
    fn modulus_bits(&self) -> usize;
    /// Fewest digits necessary to represent the modulus value `q` in the given base.
    fn modulus_digits(&self) -> usize;
    /// Returns the integer `n` that specifies the size of the polynomial ring used in this
    /// polynomial. Specifically, this is the degree parameter for the ring `Z_q[x]/(x^n - 1)`.
    fn ring_dimension(&self) -> u32;
}
