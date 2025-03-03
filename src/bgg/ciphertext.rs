// use crate::bgg::{
//     operations::{poly_add, poly_mul, poly_sub},
//     pub_key::PublicKey,
// };
// use crate::utils::empty_matrix_ring;
// use phantom_zone_crypto::util::distribution::NoiseDistribution;
// use phantom_zone_math::{
//     prelude::{Gaussian, ModulusOps, Sampler},
//     ring::RingOps,
// };
// use rand::{thread_rng, Rng};

// /// A ciphertext in the BGG+ RLWE encoding scheme
// ///
// /// # Fields
// ///
// /// * `inner`: The inner component matrix of size `(ell + 1) × m` where the rows are equal to the vectors `G+b_0, x_1G+b_1, ..., x_ellG+b_ell`
// /// * `secret`: The secret vector which is a polynomial in the ring
// /// * `error`: The error matrix of size `(ell + 1) × m` where each row is equal to `e_aS_0, e_aS_1, ..., e_aS_ell`
// #[derive(Clone, Debug)]
// pub struct Ciphertext {
//     inner: Vec<Vec<Vec<u64>>>,
//     secret: Vec<u64>,
//     error: Vec<Vec<Vec<u64>>>,
// }

// impl Ciphertext {
//     /// Create a new ciphertext by encoding an attribute vector
//     ///
//     /// # Arguments
//     /// * `public_key`: The public key matrix
//     /// * `x`: Attribute vector to encode
//     pub fn new(public_key: &PublicKey, x: &Vec<u64>) -> Self {
//         let params = public_key.params();
//         assert!(x.len() == params.ell() + 1);

//         let mut rng = thread_rng();
//         let ring = params.ring();
//         let m = *params.m();
//         let ell = *params.ell();
//         let s = ring.sample_uniform_vec(ring.ring_size(), &mut rng);
//         let mut ct_inner = public_key.b().clone();

//         let mut err_a = vec![vec![ring.zero(); ring.ring_size()]; m];
//         let gaussian: NoiseDistribution = Gaussian(3.19).into();

//         for err in err_a.iter_mut() {
//             ring.sample_into::<i64>(err, gaussian, rng.clone());
//         }

//         let mut error = empty_matrix_ring(ring, ell + 1, m);

//         for i in 0..ell + 1 {
//             for si in 0..m {
//                 for err_aj in err_a.iter() {
//                     let random_bit = if rng.gen_bool(0.5) { 1 } else { -1 };
//                     if random_bit == 1 {
//                         error[i][si] = poly_add(ring, &error[i][si], err_aj);
//                     }
//                     if random_bit == -1 {
//                         error[i][si] = poly_sub(ring, &error[i][si], err_aj);
//                     }
//                 }
//             }

//             let g = params.g();

//             // add gadget vector to ct_inner if x[i] = 1
//             if x[i] == 1 {
//                 for (ct, g_j) in ct_inner[i].iter_mut().zip(g.iter()) {
//                     *ct = poly_add(ring, ct, g_j);
//                 }
//             }
//         }

//         Self {
//             inner: ct_inner,
//             secret: s,
//             error,
//         }
//     }

//     /// Compute the ct_full as `inner * secret + error`
//     pub fn compute_ct_full(&self, pub_key: &PublicKey) -> Vec<Vec<Vec<u64>>> {
//         let params = pub_key.params();
//         let ring = params.ring();
//         let mut ct_full = self.inner().clone();

//         for i in 0..self.inner().len() {
//             for j in 0..self.inner()[0].len() {
//                 let c = poly_mul(ring, &self.inner()[i][j], self.secret());
//                 ct_full[i][j] = poly_add(ring, &c.to_vec(), &self.error()[i][j]);
//             }
//         }

//         ct_full
//     }

//     pub fn inner(&self) -> &Vec<Vec<Vec<u64>>> {
//         &self.inner
//     }

//     pub fn inner_concat(&self) -> Vec<Vec<u64>> {
//         let ct_inner = self.inner();
//         let mut concat_vec = Vec::new();

//         for ct_inner_i in ct_inner.iter() {
//             concat_vec.extend(ct_inner_i.clone());
//         }

//         concat_vec
//     }

//     pub fn secret(&self) -> &Vec<u64> {
//         &self.secret
//     }

//     pub fn error(&self) -> &Vec<Vec<Vec<u64>>> {
//         &self.error
//     }
// }
