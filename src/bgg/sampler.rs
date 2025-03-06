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
    <S as PolyHashSampler<K>>::M: Send + Sync,
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
        params: &<<<S as PolyHashSampler<K>>::M as PolyMatrix>::P as Poly>::Params,
        tag: &[u8],
        packed_input_size: usize,
    ) -> Vec<BggPublicKey<<S as PolyHashSampler<K>>::M>> {
        let log_q = params.modulus_bits();
        let columns = 2 * log_q;
        let all_matrix = self.sampler.sample_hash(
            params,
            tag,
            2,
            columns * packed_input_size,
            DistType::FinRingDist,
        );
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
    <S as PolyUniformSampler>::M: Send + Sync,
    <<S as PolyUniformSampler>::M as PolyMatrix>::P: Send + Sync,
{
    /// Create a new encoding sampler
    /// # Arguments
    /// * `secret`: The secret polynomial
    /// * `error_sampler`: The sampler to generate LWE errors
    /// * `gauss_sigma`: The standard deviation of the Gaussian distribution
    /// # Returns
    /// A new encoding sampler
    pub fn new(
        params: &<<<S as PolyUniformSampler>::M as PolyMatrix>::P as Poly>::Params,
        secret: &<S::M as PolyMatrix>::P,
        error_sampler: Arc<S>,
        gauss_sigma: f64,
    ) -> Self {
        let minus_one_poly = <S::M as PolyMatrix>::P::const_minus_one(params);
        // 1*2 row vector
        let secret_vec = S::M::from_poly_vec_row(params, vec![secret.clone(), minus_one_poly]);
        Self { secret_vec, error_sampler, gauss_sigma }
    }

    pub fn sample(
        &self,
        params: &<<<S as PolyUniformSampler>::M as PolyMatrix>::P as Poly>::Params,
        public_keys: &[BggPublicKey<S::M>],
        plaintexts: &[<S::M as PolyMatrix>::P],
        reveal_plaintexts: bool,
    ) -> Vec<BggEncoding<S::M>> {
<<<<<<< HEAD
=======
        let secret_vec = &self.secret_vec;
>>>>>>> main
        let log_q = params.modulus_bits();
        let packed_input_size = plaintexts.len();
        let columns = 2 * log_q * packed_input_size;
        let error: S::M = self.error_sampler.sample_uniform(
            params,
            1,
            columns,
            DistType::GaussDist { sigma: self.gauss_sigma },
        );
        println!("err {:?}, ", error);
        // first term sA
        // [TODO] Avoid memory cloning here.
        let all_public_key_matrix: S::M = public_keys[0]
            .matrix
            .concat_columns(&public_keys[1..].iter().map(|pk| pk.matrix.clone()).collect_vec());
        let first_term = secret_vec.clone() * all_public_key_matrix;
        // second term x \tensor sG
        let gadget = S::M::gadget_matrix(params, 2);
        let encoded_polys_vec = S::M::from_poly_vec_row(params, plaintexts.to_vec());
        let second_term = encoded_polys_vec.tensor(&(secret_vec.clone() * gadget));
        // all_vector = sA + x \tensor sG + e
        let all_vector = first_term - second_term + error;
        let encoding: Vec<BggEncoding<S::M>> = plaintexts
            .to_vec()
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
        DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyUniformSampler,
    };
    use keccak_asm::Keccak256;

    #[test]
    fn test_bgg_pub_key_sampling() {
        let input_size = 10_usize;
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
        let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
        let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into());
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, packed_input_size);
        assert_eq!(sampled_pub_keys.len(), packed_input_size);
    }

    #[test]
    fn test_bgg_pub_key_addition() {
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = 2;
        let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
        let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into());
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, packed_input_size);
        let log_q = params.modulus_bits();
        let columns = 2 * log_q;

        for pair in sampled_pub_keys.chunks(2) {
            if let [a, b] = pair {
                let addition = a.clone() + b.clone();
                assert_eq!(addition.matrix.row_size(), 2);
                assert_eq!(addition.matrix.col_size(), columns);
                assert_eq!(addition.matrix, a.matrix.clone() + b.matrix.clone());
            }
        }
    }

    #[test]
    fn test_bgg_pub_key_multiplication() {
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = 2;
        let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
        let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into());
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, packed_input_size);
        let log_q = params.modulus_bits();
        let columns = 2 * log_q;

        for pair in sampled_pub_keys.chunks(2) {
            if let [a, b] = pair {
                let multiplication = a.clone() * b.clone();
                assert_eq!(multiplication.matrix.row_size(), 2);
                assert_eq!(multiplication.matrix.col_size(), columns);
                assert_eq!(multiplication.matrix, (a.matrix.clone() * b.matrix.decompose().clone()))
            }
        }
    }

    #[test]
    fn test_bgg_encoding_sampling() {
        let input_size = 10_usize;
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
        let bgg_sampler =
            BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into());
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, packed_input_size);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let plaintexts = vec![DCRTPoly::const_one(&params); packed_input_size];
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts, true);
        let g = DCRTPolyMatrix::gadget_matrix(&params, 2);
        assert_eq!(bgg_encodings.len(), packed_input_size);
        assert_eq!(
            bgg_encodings[0].vector,
            bgg_sampler.secret_vec.clone() * bgg_encodings[0].pubkey.matrix.clone()
                - bgg_sampler.secret_vec.clone()
                    * (g * bgg_encodings[0].plaintext.clone().unwrap())
        )
    }

    #[test]
    fn test_bgg_encoding_addition() {
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = 2;
        let bgg_sampler =
            BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into());
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, packed_input_size);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let plaintexts = vec![DCRTPoly::const_one(&params); packed_input_size];
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts, true);

        for pair in bgg_encodings.chunks(2) {
            if let [a, b] = pair {
                let addition = a.clone() + b.clone();
                assert_eq!(addition.pubkey, a.pubkey.clone() + b.pubkey.clone());
                assert_eq!(
                    addition.clone().plaintext.unwrap(),
                    a.plaintext.clone().unwrap() + b.plaintext.clone().unwrap()
                );
                let g = DCRTPolyMatrix::gadget_matrix(&params, 2);
                assert_eq!(addition.vector, a.clone().vector + b.clone().vector);
                assert_eq!(
                    addition.vector,
                    bgg_sampler.secret_vec.clone()
                        * (addition.pubkey.matrix - (g * addition.plaintext.unwrap()))
                )
            }
        }
    }

    #[test]
    fn test_bgg_encoding_multiplication() {
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = 2;
        let bgg_sampler =
            BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into());
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, packed_input_size);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let plaintexts = vec![DCRTPoly::const_one(&params); packed_input_size];
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts, true);

        for pair in bgg_encodings.chunks(2) {
            if let [a, b] = pair {
                let multiplication = a.clone() * b.clone();
                assert_eq!(multiplication.pubkey, (a.clone().pubkey * b.clone().pubkey));
                assert_eq!(
                    multiplication.clone().plaintext.unwrap(),
                    a.clone().plaintext.unwrap() * b.clone().plaintext.unwrap()
                );
                let g = DCRTPolyMatrix::gadget_matrix(&params, 2);
                assert_eq!(
                    multiplication.vector,
                    (bgg_sampler.secret_vec.clone()
                        * (multiplication.pubkey.matrix - (g * multiplication.plaintext.unwrap())))
                )
            }
        }
    }
}
