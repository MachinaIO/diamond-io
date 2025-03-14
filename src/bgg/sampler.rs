use super::{BggEncoding, BggPublicKey};
use crate::{
    parallel_iter,
    poly::{
        params::PolyParams,
        polynomial::Poly,
        sampler::{DistType, PolyHashSampler, PolyUniformSampler},
        PolyMatrix,
    },
};
use itertools::Itertools;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
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
    /// * `reveal_plaintexts`: A vector of booleans indicating whether the plaintexts associated to
    ///   the public keys should be revealed
    /// # Returns
    /// A vector of public key matrices
    pub fn sample(
        &self,
        params: &<<<S as PolyHashSampler<K>>::M as PolyMatrix>::P as Poly>::Params,
        tag: &[u8],
        reveal_plaintexts: &[bool],
    ) -> Vec<BggPublicKey<<S as PolyHashSampler<K>>::M>> {
        let log_q = params.modulus_bits();
        let columns = 2 * log_q;
        let packed_input_size = 1 + reveal_plaintexts.len(); // first slot is allocated to the constant 1 polynomial plaintext
        let all_matrix = self.sampler.sample_hash(
            params,
            tag,
            2,
            columns * packed_input_size,
            DistType::FinRingDist,
        );

        parallel_iter!(0..packed_input_size)
            .map(|idx| {
                let reveal_plaintext = if idx == 0 { true } else { reveal_plaintexts[idx - 1] };
                BggPublicKey::new(
                    all_matrix.slice_columns(columns * idx, columns * (idx + 1)),
                    reveal_plaintext,
                )
            })
            .collect()
    }
}

/// A sampler of an encoding in the BGG+ RLWE encoding scheme
///
/// # Fields
/// * `secret`: The secret vector.
/// * `error_sampler`: The sampler to generate RLWE errors.
/// * `gauss_sigma`: The standard deviation of the Gaussian distribution.
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
    /// * `error_sampler`: The sampler to generate RLWE errors
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
    ) -> Vec<BggEncoding<S::M>> {
        let secret_vec = &self.secret_vec;
        let log_q = params.modulus_bits();
        let packed_input_size = 1 + plaintexts.len(); // first slot is allocated to the constant 1 polynomial plaintext
        let plaintexts: Vec<<S::M as PolyMatrix>::P> =
            [&[<<S as PolyUniformSampler>::M as PolyMatrix>::P::const_one(params)], plaintexts]
                .concat();
        let columns = 2 * log_q * packed_input_size;
        let error: S::M = self.error_sampler.sample_uniform(
            params,
            1,
            columns,
            DistType::GaussDist { sigma: self.gauss_sigma },
        );

        let all_public_key_matrix: S::M = public_keys[0]
            .matrix
            .concat_columns(&public_keys[1..].iter().map(|pk| pk.matrix.clone()).collect_vec());
        let first_term = secret_vec.clone() * all_public_key_matrix;

        let gadget = S::M::gadget_matrix(params, 2);
        let encoded_polys_vec = S::M::from_poly_vec_row(params, plaintexts.to_vec());
        let second_term = encoded_polys_vec.tensor(&(secret_vec.clone() * gadget));

        let all_vector = first_term - second_term + error;

        parallel_iter!(plaintexts)
            .enumerate()
            .map(|(idx, plaintext)| {
                let vector = all_vector.slice_columns(2 * log_q * idx, 2 * log_q * (idx + 1));
                BggEncoding {
                    vector,
                    pubkey: public_keys[idx].clone(),
                    plaintext: if public_keys[idx].reveal_plaintext {
                        Some(plaintext.clone())
                    } else {
                        None
                    },
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bgg::eval::Evaluable,
        poly::dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyUniformSampler,
        },
        utils::{create_bit_random_poly, create_random_poly},
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
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        assert_eq!(sampled_pub_keys.len(), packed_input_size + 1);
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
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let log_q = params.modulus_bits();
        let columns = 2 * log_q;

        for pair in sampled_pub_keys[1..].chunks(2) {
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
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let log_q = params.modulus_bits();
        let columns = 2 * log_q;

        for pair in sampled_pub_keys[1..].chunks(2) {
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
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![DCRTPoly::const_one(&params); packed_input_size];
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);
        let g = DCRTPolyMatrix::gadget_matrix(&params, 2);
        assert_eq!(bgg_encodings.len(), packed_input_size + 1);
        assert_eq!(
            bgg_encodings[0].vector,
            bgg_sampler.secret_vec.clone() * bgg_encodings[0].pubkey.matrix.clone() -
                bgg_sampler.secret_vec.clone() *
                    (g.clone() * bgg_encodings[0].plaintext.clone().unwrap())
        );
        assert_eq!(
            bgg_encodings[1].vector,
            bgg_sampler.secret_vec.clone() * bgg_encodings[1].pubkey.matrix.clone() -
                bgg_sampler.secret_vec.clone() *
                    (g * bgg_encodings[1].plaintext.clone().unwrap())
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
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params); packed_input_size];
        // TODO: set the standard deviation to a non-zero value
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);

        for pair in bgg_encodings[1..].chunks(2) {
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
                    bgg_sampler.secret_vec.clone() *
                        (addition.pubkey.matrix - (g * addition.plaintext.unwrap()))
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
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params); packed_input_size];
        // TODO: set the standard deviation to a non-zero value
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);

        for pair in bgg_encodings[1..].chunks(2) {
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
                    (bgg_sampler.secret_vec.clone() *
                        (multiplication.pubkey.matrix - (g * multiplication.plaintext.unwrap())))
                )
            }
        }
    }

    #[test]
    fn test_bgg_encoding_scalar_multiplication() {
        let key: [u8; 32] = rand::random();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let params = DCRTPolyParams::default();
        let packed_input_size = 1;
        let bgg_sampler =
            BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into());
        let reveal_plaintexts = vec![true; packed_input_size];
        let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params); packed_input_size];
        // TODO: set the standard deviation to a non-zero value
        let bgg_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler.into(), 0.0);
        let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);

        // Create a scalar (polynomial) for scalar multiplication
        let scalar = create_random_poly(&params);

        // Test scalar multiplication for each encoding
        for encoding in &bgg_encodings[1..] {
            // Perform scalar multiplication
            let scalar_mul = encoding.scalar_mul(&params, &scalar);

            // Verify the pubkey is correctly transformed
            assert_eq!(scalar_mul.pubkey, encoding.pubkey.scalar_mul(&params, &scalar));

            // Verify the plaintext is correctly multiplied by the scalar
            assert_eq!(
                scalar_mul.plaintext.as_ref().unwrap(),
                &(encoding.plaintext.as_ref().unwrap().clone() * scalar.clone())
            );

            // Verify the vector is correctly transformed
            let g = DCRTPolyMatrix::gadget_matrix(&params, 2);
            let gadget_scalar = g.clone() * &scalar;
            let decomposed_gadget_scalar = gadget_scalar.decompose();

            // The vector should be the original vector multiplied by the decomposed gadget scalar
            assert_eq!(scalar_mul.vector, encoding.vector.clone() * decomposed_gadget_scalar);

            // Alternative verification: check that the vector satisfies the BGG encoding relation
            assert_eq!(
                scalar_mul.vector,
                bgg_sampler.secret_vec.clone() *
                    (scalar_mul.pubkey.matrix.clone() -
                        (g * scalar_mul.plaintext.as_ref().unwrap().clone()))
            );
        }
    }
}
