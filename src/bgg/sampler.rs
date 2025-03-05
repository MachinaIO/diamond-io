use super::{BggEncoding, BggPublicKey};
use crate::poly::{matrix::*, sampler::*, *};
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::{marker::PhantomData, sync::Arc};

/// A sampler of a public key A in the BGG+ RLWE encoding scheme
#[derive(Clone)]
pub struct BGGPublicKeySampler<K: AsRef<[u8]>, S: PolyHashSampler<K>> {
    pub sampler: Arc<S>,
    _k: PhantomData<K>,
}

impl<K: AsRef<[u8]>, S> BGGPublicKeySampler<K, S>
where
    S: PolyHashSampler<K>,
    <S as PolySampler>::M: Send + Sync,
{
    /// Create a new public key sampler
    /// # Arguments
    /// * `sampler`: The sampler to generate the public key matrix
    /// # Returns
    /// A new public key sampler
    pub fn new(sampler: Arc<S>) -> Self {
        Self { sampler, _k: PhantomData }
    }

    /// Sample a public key matrix
    /// # Arguments
    /// * `tag`: The tag to sample the public key matrix
    /// * `packed_input_size`: The packed input size, i.e., the number of necessary polynomials when
    ///   n bits are encoded into a single polynomial
    /// # Returns
    /// A vector of public key matrices
    pub fn sample(
        &self,
        tag: &[u8],
        packed_input_size: usize,
    ) -> Vec<BggPublicKey<<S as PolySampler>::M>> {
        let log_q = self.sampler.get_params().modulus_bits();
        let columns = 2 * log_q;
        let all_matrix =
            self.sampler.sample_hash(tag, 2, columns * packed_input_size, DistType::FinRingDist);
        (0..packed_input_size)
            .into_par_iter()
            .map(|idx| {
                BggPublicKey::new(all_matrix.slice_columns(columns * idx, columns * (idx + 1)))
            })
            .collect()
    }
}

/// A sampler of an encoding in the BGG+ RLWE encoding scheme
///
/// # Fields
/// * `secret`: The secret vector.
/// * `error_sampler`: The sampler to generate LWE errors.
#[derive(Clone)]
pub struct BGGEncodingSampler<S: PolyUniformSampler> {
    pub(crate) secret_vec: S::M,
    pub error_sampler: Arc<S>,
    pub gauss_sigma: f64,
}

impl<S> BGGEncodingSampler<S>
where
    S: PolyUniformSampler,
    <S as PolySampler>::M: Send + Sync,
    <<S as PolySampler>::M as PolyMatrix>::P: Send + Sync,
{
    /// Create a new encoding sampler
    /// # Arguments
    /// * `secret`: The secret polynomial
    /// * `error_sampler`: The sampler to generate LWE errors
    /// * `gauss_sigma`: The standard deviation of the Gaussian distribution
    /// # Returns
    /// A new encoding sampler
    pub fn new(secret: &<S::M as PolyMatrix>::P, error_sampler: Arc<S>, gauss_sigma: f64) -> Self {
        let params = error_sampler.get_params();
        let minus_one_poly = <S::M as PolyMatrix>::P::const_minus_one(&params);
        // 1*2 row vector
        let secret_vec = S::M::from_poly_vec_row(&params, vec![secret.clone(), minus_one_poly]);
        Self { secret_vec, error_sampler, gauss_sigma }
    }

    pub fn sample(
        &self,
        public_keys: &[BggPublicKey<S::M>],
        plaintexts: &[<S::M as PolyMatrix>::P],
        reveal_plaintexts: bool,
    ) -> Vec<BggEncoding<S::M>> {
        let params = self.error_sampler.get_params();
        let log_q = params.modulus_bits();
        let packed_input_size = plaintexts.len();
        let columns = 2 * log_q * packed_input_size;
        let error = self.error_sampler.sample_uniform(
            1,
            columns,
            DistType::GaussDist { sigma: self.gauss_sigma },
        );
        // first term sA
        // [TODO] Avoid memory cloning here.
        let all_public_key_matrix: S::M = public_keys[0]
            .matrix
            .concat_columns(&public_keys[1..].iter().map(|pk| pk.matrix.clone()).collect_vec());
        let first_term = self.secret_vec.clone() * all_public_key_matrix;
        // second term x \tensor sG
        let gadget = S::M::gadget_matrix(&params, 2);
        let sg = self.secret_vec.clone() * gadget;
        let encoded_polys_vec = S::M::from_poly_vec_row(&params, plaintexts.to_vec());
        let second_term = encoded_polys_vec.tensor(&sg);

        // all_vector = sA + x \tensor sG + e
        let all_vector = first_term + second_term + error;

        let encoding: Vec<BggEncoding<S::M>> = plaintexts
            .to_vec() // Convert &[T] to Vec<T> if needed
            .into_par_iter()
            .enumerate()
            .map(|(idx, plaintext)| {
                let vector = all_vector.slice_columns(2 * log_q * idx, 2 * log_q * (idx + 1));
                BggEncoding {
                    vector,
                    pubkey: public_keys[idx].clone(),
                    plaintext: if reveal_plaintexts { Some(plaintext.clone()) } else { None },
                }
            })
            .collect();
        encoding
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt::{
        DCRTPoly, DCRTPolyHashSampler, DCRTPolyParams, DCRTPolyUniformSampler,
    };
    use keccak_asm::Keccak256;

    #[test]
    fn test_bgg_pub_key_sampler() {
        let input_size = 10_usize;
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
        let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key, params);
        let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into());
        let sampled_pub_keys = bgg_sampler.sample(&tag_bytes, packed_input_size);
        let log_q = bgg_sampler.sampler.get_params().modulus_bits();
        let columns = 2 * log_q;
        assert_eq!(sampled_pub_keys.len(), packed_input_size);
        for m in sampled_pub_keys {
            assert_eq!(m.matrix.row_size(), 2);
            assert_eq!(m.matrix.col_size(), columns);
        }
    }

    #[test]
    fn test_bgg_encoding_sampler() {
        let input_size = 10_usize;
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
        let bgg_sampler = BGGPublicKeySampler::new(
            DCRTPolyHashSampler::<Keccak256>::new(key, params.clone()).into(),
        );
        let sampled_pub_keys = bgg_sampler.sample(&tag_bytes, packed_input_size);
        let uniform_sampler = DCRTPolyUniformSampler::new(params.clone());
        let secret = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let plaintexts = vec![DCRTPoly::const_one(&params); packed_input_size];
        let bgg_sampler = BGGEncodingSampler::new(&secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&sampled_pub_keys, &plaintexts, false);
        assert_eq!(bgg_encodings.len(), packed_input_size);
    }
}
