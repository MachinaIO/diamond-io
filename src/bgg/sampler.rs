// use super::{BggEncoding, BggPublicKey};
// use crate::{
//     parallel_iter,
//     poly::{
//         sampler::{DistType, PolyHashSampler, PolyUniformSampler},
//         Poly, PolyMatrix, PolyParams,
//     },
//     utils::debug_mem,
// };
// use rayon::prelude::*;
// use std::{marker::PhantomData, sync::Arc};

// /// A sampler of a public key A in the BGG+ RLWE encoding scheme
// #[derive(Clone)]
// pub struct BGGPublicKeySampler<K: AsRef<[u8]>, S: PolyHashSampler<K>> {
//     pub sampler: Arc<S>,
//     pub d: usize,
//     _k: PhantomData<K>,
// }

// impl<K: AsRef<[u8]>, S> BGGPublicKeySampler<K, S>
// where
//     S: PolyHashSampler<K>,
// {
//     /// Create a new public key sampler
//     /// # Arguments
//     /// * `sampler`: The sampler to generate the public key matrix
//     /// * `d`: The number of secret polynomials used with the sampled public key matrices.
//     /// # Returns
//     /// A new public key sampler
//     pub fn new(sampler: Arc<S>, d: usize) -> Self {
//         Self { sampler, d, _k: PhantomData }
//     }

//     /// Sample a public key matrix
//     /// # Arguments
//     /// * `tag`: The tag to sample the public key matrix
//     /// * `reveal_plaintexts`: A vector of booleans indicating whether the plaintexts associated
// to     ///   the public keys should be revealed
//     /// # Returns
//     /// A vector of public key matrices
//     pub fn sample(
//         &self,
//         params: &<<<S as PolyHashSampler<K>>::M as PolyMatrix>::P as Poly>::Params,
//         tag: &[u8],
//         reveal_plaintexts: &[bool],
//     ) -> Vec<BggPublicKey<<S as PolyHashSampler<K>>::M>> {
//         let log_base_q = params.modulus_digits();
//         let secret_vec_size = self.d + 1;
//         let columns = secret_vec_size * log_base_q;
//         let packed_input_size = 1 + reveal_plaintexts.len(); // first slot is allocated to the
// constant 1 polynomial plaintext         let all_matrix = self.sampler.sample_hash(
//             params,
//             tag,
//             secret_vec_size,
//             columns * packed_input_size,
//             DistType::FinRingDist,
//         );
//         parallel_iter!(0..packed_input_size)
//             .map(|idx| {
//                 let reveal_plaintext = if idx == 0 { true } else { reveal_plaintexts[idx - 1] };
//                 BggPublicKey::new(
//                     all_matrix.slice_columns(columns * idx, columns * (idx + 1)),
//                     reveal_plaintext,
//                 )
//             })
//             .collect()
//     }
// }

// /// A sampler of an encoding in the BGG+ RLWE encoding scheme
// ///
// /// # Fields
// /// * `secret`: The secret vector.
// /// * `error_sampler`: The sampler to generate RLWE errors.
// /// * `gauss_sigma`: The standard deviation of the Gaussian distribution.
// #[derive(Clone)]
// pub struct BGGEncodingSampler<S: PolyUniformSampler> {
//     pub(crate) secret_vec: S::M,
//     pub error_sampler: Arc<S>,
//     pub gauss_sigma: f64,
// }

// impl<S> BGGEncodingSampler<S>
// where
//     S: PolyUniformSampler,
// {
//     /// Create a new encoding sampler
//     /// # Arguments
//     /// * `secrets`: The secret polynomials
//     /// * `error_sampler`: The sampler to generate RLWE errors
//     /// * `gauss_sigma`: The standard deviation of the Gaussian distribution
//     /// # Returns
//     /// A new encoding sampler
//     pub fn new(
//         params: &<<<S as PolyUniformSampler>::M as PolyMatrix>::P as Poly>::Params,
//         secrets: &[<S::M as PolyMatrix>::P],
//         error_sampler: Arc<S>,
//         gauss_sigma: f64,
//     ) -> Self {
//         let minus_one_poly = <S::M as PolyMatrix>::P::const_minus_one(params);
//         // 1*(d+1) row vector
//         let mut secrets = secrets.to_vec();
//         secrets.push(minus_one_poly);
//         let secret_vec = S::M::from_poly_vec_row(params, secrets);
//         Self { secret_vec, error_sampler, gauss_sigma }
//     }

//     pub fn sample(
//         &self,
//         params: &<<<S as PolyUniformSampler>::M as PolyMatrix>::P as Poly>::Params,
//         public_keys: &[BggPublicKey<S::M>],
//         plaintexts: &[<S::M as PolyMatrix>::P],
//     ) -> Vec<BggEncoding<S::M>> {
//         let secret_vec = &self.secret_vec;
//         let log_base_q = params.modulus_digits();
//         let packed_input_size = 1 + plaintexts.len(); // first slot is allocated to the constant
// 1 polynomial plaintext         let plaintexts: Vec<<S::M as PolyMatrix>::P> =
//             [&[<<S as PolyUniformSampler>::M as PolyMatrix>::P::const_one(params)], plaintexts]
//                 .concat();
//         let secret_vec_size = self.secret_vec.col_size();
//         let columns = secret_vec_size * log_base_q * packed_input_size;
//         let error: S::M = self.error_sampler.sample_uniform(
//             params,
//             1,
//             columns,
//             DistType::GaussDist { sigma: self.gauss_sigma },
//         );
//         let all_public_key_matrix: S::M = public_keys[0]
//             .matrix
//             .concat_columns(&public_keys[1..].par_iter().map(|pk|
// &pk.matrix).collect::<Vec<_>>());         let first_term = secret_vec.clone() *
// all_public_key_matrix;

//         let gadget = S::M::gadget_matrix(params, secret_vec_size);
//         let encoded_polys_vec = S::M::from_poly_vec_row(params, plaintexts.to_vec());
//         let second_term = encoded_polys_vec.tensor(&(secret_vec.clone() * gadget));
//         let all_vector = first_term - second_term + error;

//         let m = secret_vec_size * log_base_q;
//         parallel_iter!(plaintexts)
//             .enumerate()
//             .map(|(idx, plaintext)| {
//                 let vector = all_vector.slice_columns(m * idx, m * (idx + 1));
//                 debug_mem("before constructing BggEncoding");
//                 BggEncoding {
//                     vector,
//                     pubkey: public_keys[idx].clone(),
//                     plaintext: if public_keys[idx].reveal_plaintext {
//                         Some(plaintext.clone())
//                     } else {
//                         None
//                     },
//                 }
//             })
//             .collect()
//     }
// }

// #[cfg(test)]
// #[cfg(feature = "test")]
// mod tests {
//     use super::*;
//     use crate::{
//         poly::dcrt::{
//             DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams,
// DCRTPolyUniformSampler,         },
//         utils::{create_bit_random_poly, create_random_poly},
//     };
//     use keccak_asm::Keccak256;

//     #[test]
//     fn test_bgg_pub_key_sampling() {
//         let input_size = 10_usize;
//         let key: [u8; 32] = rand::random();
//         let tag: u64 = rand::random();
//         let tag_bytes = tag.to_le_bytes();
//         let params = DCRTPolyParams::default();
//         let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
//         let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
//         let d = 3;
//         let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into(), d);
//         let reveal_plaintexts = vec![true; packed_input_size];
//         let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
//         assert_eq!(sampled_pub_keys.len(), packed_input_size + 1);
//     }

//     #[test]
//     fn test_bgg_pub_key_addition() {
//         let key: [u8; 32] = rand::random();
//         let tag: u64 = rand::random();
//         let tag_bytes = tag.to_le_bytes();
//         let params = DCRTPolyParams::default();
//         let packed_input_size = 2;
//         let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
//         let d = 3;
//         let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into(), d);
//         let reveal_plaintexts = vec![true; packed_input_size];
//         let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
//         let log_base_q = params.modulus_digits();
//         let columns = (d + 1) * log_base_q;

//         for pair in sampled_pub_keys[1..].chunks(2) {
//             if let [a, b] = pair {
//                 let addition = a.clone() + b.clone();
//                 assert_eq!(addition.matrix.row_size(), d + 1);
//                 assert_eq!(addition.matrix.col_size(), columns);
//                 assert_eq!(addition.matrix, a.matrix.clone() + b.matrix.clone());
//             }
//         }
//     }

//     #[test]
//     fn test_bgg_pub_key_multiplication() {
//         let key: [u8; 32] = rand::random();
//         let tag: u64 = rand::random();
//         let tag_bytes = tag.to_le_bytes();
//         let params = DCRTPolyParams::default();
//         let packed_input_size = 2;
//         let poly_hash_sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
//         let d = 3;
//         let bgg_sampler = BGGPublicKeySampler::new(poly_hash_sampler.into(), d);
//         let reveal_plaintexts = vec![true; packed_input_size];
//         let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
//         let log_base_q = params.modulus_digits();
//         let columns = (d + 1) * log_base_q;

//         for pair in sampled_pub_keys[1..].chunks(2) {
//             if let [a, b] = pair {
//                 let multiplication = a.clone() * b.clone();
//                 assert_eq!(multiplication.matrix.row_size(), d + 1);
//                 assert_eq!(multiplication.matrix.col_size(), columns);
//                 assert_eq!(multiplication.matrix, (a.matrix.clone() *
// b.matrix.decompose().clone()))             }
//         }
//     }

//     #[test]
//     fn test_bgg_encoding_sampling() {
//         let input_size = 10_usize;
//         let key: [u8; 32] = rand::random();
//         let tag: u64 = rand::random();
//         let tag_bytes = tag.to_le_bytes();
//         let params = DCRTPolyParams::default();
//         let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
//         let d = 3;
//         let bgg_sampler =
//             BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into(), d);
//         let reveal_plaintexts = vec![true; packed_input_size];
//         let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
//         let uniform_sampler = DCRTPolyUniformSampler::new();
//         let secrets = vec![create_bit_random_poly(&params); d];
//         let plaintexts = vec![DCRTPoly::const_one(&params); packed_input_size];
//         let bgg_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler.into(),
// 0.0);         let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);
//         let g = DCRTPolyMatrix::gadget_matrix(&params, d + 1);
//         assert_eq!(bgg_encodings.len(), packed_input_size + 1);
//         assert_eq!(
//             bgg_encodings[0].vector,
//             bgg_sampler.secret_vec.clone() * bgg_encodings[0].pubkey.matrix.clone() -
//                 bgg_sampler.secret_vec.clone() *
//                     (g.clone() * bgg_encodings[0].plaintext.clone().unwrap())
//         );
//         assert_eq!(
//             bgg_encodings[1].vector,
//             bgg_sampler.secret_vec.clone() * bgg_encodings[1].pubkey.matrix.clone() -
//                 bgg_sampler.secret_vec.clone() *
//                     (g * bgg_encodings[1].plaintext.clone().unwrap())
//         )
//     }

//     #[test]
//     fn test_bgg_encoding_addition() {
//         let key: [u8; 32] = rand::random();
//         let tag: u64 = rand::random();
//         let tag_bytes = tag.to_le_bytes();
//         let params = DCRTPolyParams::default();
//         let packed_input_size = 2;
//         let d = 3;
//         let bgg_sampler =
//             BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into(), d);
//         let reveal_plaintexts = vec![true; packed_input_size];
//         let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
//         let uniform_sampler = DCRTPolyUniformSampler::new();
//         let secrets = vec![create_bit_random_poly(&params); d];
//         let plaintexts = vec![create_random_poly(&params); packed_input_size];
//         // TODO: set the standard deviation to a non-zero value
//         let bgg_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler.into(),
// 0.0);         let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);

//         for pair in bgg_encodings[1..].chunks(2) {
//             if let [a, b] = pair {
//                 let addition = a.clone() + b.clone();
//                 assert_eq!(addition.pubkey, a.pubkey.clone() + b.pubkey.clone());
//                 assert_eq!(
//                     addition.clone().plaintext.unwrap(),
//                     a.plaintext.clone().unwrap() + b.plaintext.clone().unwrap()
//                 );
//                 let g = DCRTPolyMatrix::gadget_matrix(&params, d + 1);
//                 assert_eq!(addition.vector, a.clone().vector + b.clone().vector);
//                 assert_eq!(
//                     addition.vector,
//                     bgg_sampler.secret_vec.clone() *
//                         (addition.pubkey.matrix - (g * addition.plaintext.unwrap()))
//                 )
//             }
//         }
//     }

//     #[test]
//     fn test_bgg_encoding_multiplication() {
//         let key: [u8; 32] = rand::random();
//         let tag: u64 = rand::random();
//         let tag_bytes = tag.to_le_bytes();
//         let params = DCRTPolyParams::default();
//         let packed_input_size = 2;
//         let d = 3;
//         let bgg_sampler =
//             BGGPublicKeySampler::new(DCRTPolyHashSampler::<Keccak256>::new(key).into(), d);
//         let reveal_plaintexts = vec![true; packed_input_size];
//         let sampled_pub_keys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
//         let uniform_sampler = DCRTPolyUniformSampler::new();
//         let secrets = vec![create_bit_random_poly(&params); d];
//         let plaintexts = vec![create_random_poly(&params); packed_input_size];
//         let bgg_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler.into(),
// 0.0);         let bgg_encodings = bgg_sampler.sample(&params, &sampled_pub_keys, &plaintexts);

//         for pair in bgg_encodings[1..].chunks(2) {
//             if let [a, b] = pair {
//                 let multiplication = a.clone() * b.clone();
//                 assert_eq!(multiplication.pubkey, (a.clone().pubkey * b.clone().pubkey));
//                 assert_eq!(
//                     multiplication.clone().plaintext.unwrap(),
//                     a.clone().plaintext.unwrap() * b.clone().plaintext.unwrap()
//                 );
//                 let g = DCRTPolyMatrix::gadget_matrix(&params, d + 1);
//                 assert_eq!(
//                     multiplication.vector,
//                     (bgg_sampler.secret_vec.clone() *
//                         (multiplication.pubkey.matrix - (g * multiplication.plaintext.unwrap())))
//                 )
//             }
//         }
//     }
// }
