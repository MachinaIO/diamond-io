// #[derive(Debug, Clone)]
// pub struct BGGEvaluator<
//     T: PolyElemOps,
//     P: PolyOps<T>,
//     M: PolyMatrixOps<T, P>,
//     G: PolyGadgetOps<T, P, M>,
// > {
//     pub poly_op: Arc<P>,
//     pub matrix_op: Arc<M>,
//     pub gadget_op: Arc<G>,
//     _t: PhantomData<T>,
// }

// impl<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, G: PolyGadgetOps<T, P, M>>
//     BGGEvaluator<T, P, M, G>
// {
//     pub fn new(poly_op: Arc<P>, matrix_op: Arc<M>, gadget_op: Arc<G>) -> Self {
//         Self { poly_op, matrix_op, gadget_op, _t: PhantomData }
//     }

//     pub fn add_public_keys(
//         &self,
//         a: &BggPublicKey<T, P, M>,
//         b: &BggPublicKey<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggPublicKey<T, P, M>, BggError> {
//         let new_matrix = self
//             .matrix_op
//             .add(&a.matrix, &b.matrix)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         Ok(BggPublicKey { matrix: new_matrix, index: new_index })
//     }

//     pub fn neg_public_keys(
//         &self,
//         a: &BggPublicKey<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggPublicKey<T, P, M>, BggError> {
//         let new_matrix =
//             self.matrix_op.neg(&a.matrix).map_err(|e| BggError::MatrixError(e.to_string()))?;
//         Ok(BggPublicKey { matrix: new_matrix, index: new_index })
//     }

//     pub fn sub_public_keys(
//         &self,
//         a: &BggPublicKey<T, P, M>,
//         b: &BggPublicKey<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggPublicKey<T, P, M>, BggError> {
//         let new_matrix = self
//             .matrix_op
//             .sub(&a.matrix, &b.matrix)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         Ok(BggPublicKey { matrix: new_matrix, index: new_index })
//     }

//     pub fn mul_public_keys(
//         &self,
//         a: &BggPublicKey<T, P, M>,
//         b: &BggPublicKey<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggPublicKey<T, P, M>, BggError> {
//         let decomposed_b = self
//             .gadget_op
//             .decompose(&b.matrix)
//             .map_err(|e| BggError::GadgetError(e.to_string()))?;
//         let new_matrix = self
//             .matrix_op
//             .mul(&a.matrix, &decomposed_b)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         Ok(BggPublicKey { matrix: new_matrix, index: new_index })
//     }

//     pub fn add_encodings(
//         &self,
//         a: &BggEncoding<T, P, M>,
//         b: &BggEncoding<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggEncoding<T, P, M>, BggError> {
//         let new_vector = self
//             .matrix_op
//             .add(&a.vector, &b.vector)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         // let new_plaintext = if let Some()
//         let new_plaintext = match (a.plaintext.as_ref(), b.plaintext.as_ref()) {
//             (Some(a_plain), Some(b_plain)) => Some(
//                 self.poly_op
//                     .add(a_plain, b_plain)
//                     .map_err(|e| BggError::PolyError(e.to_string()))?,
//             ),
//             _ => None,
//         };
//         Ok(BggEncoding { vector: new_vector, plaintext: new_plaintext, index: new_index })
//     }

//     pub fn neg_encoding(
//         &self,
//         a: &BggEncoding<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggEncoding<T, P, M>, BggError> {
//         let new_vector =
//             self.matrix_op.neg(&a.vector).map_err(|e| BggError::MatrixError(e.to_string()))?;
//         let new_plaintext = match a.plaintext.as_ref() {
//             Some(a_plain) => {
//                 Some(self.poly_op.neg(a_plain).map_err(|e| BggError::PolyError(e.to_string()))?)
//             }
//             _ => None,
//         };
//         Ok(BggEncoding { vector: new_vector, plaintext: new_plaintext, index: new_index })
//     }

//     pub fn sub_encodings(
//         &self,
//         a: &BggEncoding<T, P, M>,
//         b: &BggEncoding<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggEncoding<T, P, M>, BggError> {
//         let new_vector = self
//             .matrix_op
//             .sub(&a.vector, &b.vector)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         let new_plaintext = match (a.plaintext.as_ref(), b.plaintext.as_ref()) {
//             (Some(a_plain), Some(b_plain)) => Some(
//                 self.poly_op
//                     .sub(a_plain, b_plain)
//                     .map_err(|e| BggError::PolyError(e.to_string()))?,
//             ),
//             _ => None,
//         };
//         Ok(BggEncoding { vector: new_vector, plaintext: new_plaintext, index: new_index })
//     }

//     pub fn mul_encodings(
//         &self,
//         a: &BggEncoding<T, P, M>,
//         b: &BggEncoding<T, P, M>,
//         new_index: usize,
//     ) -> Result<BggEncoding<T, P, M>, BggError> {
//         if a.plaintext.is_none() {
//             return Err(BggError::UnknownPlaintextForMul(a.index, new_index));
//         }
//         let decomposed_b = self
//             .gadget_op
//             .decompose(&b.vector)
//             .map_err(|e| BggError::GadgetError(e.to_string()))?;
//         let first_term = self
//             .matrix_op
//             .mul(&a.vector, &decomposed_b)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         let second_term = self
//             .matrix_op
//             .scalar_mul(&b.vector, a.plaintext.as_ref().unwrap())
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         let new_vector = self
//             .matrix_op
//             .add(&first_term, &second_term)
//             .map_err(|e| BggError::MatrixError(e.to_string()))?;
//         let new_plaintext = match (a.plaintext.as_ref(), b.plaintext.as_ref()) {
//             (Some(a_plain), Some(b_plain)) => Some(
//                 self.poly_op
//                     .mul(a_plain, b_plain)
//                     .map_err(|e| BggError::PolyError(e.to_string()))?,
//             ),
//             _ => None,
//         };
//         Ok(BggEncoding { vector: new_vector, plaintext: new_plaintext, index: new_index })
//     }
// }

// use phantom_zone_math::{prelude::ModulusOps, ring::RingOps};

// use crate::{
//     operations::{
//         bit_decompose, gen_identity_matrix_to_scalar, poly_add, poly_mul, vec_mat_mul,
// vec_vec_add,     },
//     parameters::Parameters,
//     utils::empty_vector_ring,
// };

// /// `out = b_left + b_right`
// pub fn m_eval_add(
//     params: &Parameters,
//     b_left: &Vec<Vec<u64>>,
//     b_right: &Vec<Vec<u64>>,
// ) -> Vec<Vec<u64>> {
//     let ring = params.ring();
//     let m = *params.m();
//     let mut out = empty_vector_ring(ring, m);
//     for i in 0..m {
//         out[i] = poly_add(ring, &b_left[i], &b_right[i]);
//     }

//     out
// }

// /// `out = b_right * tau(-b_left)`
// pub fn m_eval_mul(
//     params: &Parameters,
//     b_left: &Vec<Vec<u64>>,
//     b_right: &Vec<Vec<u64>>,
// ) -> Vec<Vec<u64>> {
//     let ring = params.ring();
//     let m = *params.m();
//     let mut out = empty_vector_ring(ring, m);

//     let mut minus_b_left = vec![vec![ring.zero(); ring.ring_size()]; m];
//     for i in 0..m {
//         for j in 0..ring.ring_size() {
//             // To get -1 * coefficient in the ring, we subtract the coefficient from 0
//             minus_b_left[i][j] = ring.sub(&ring.zero(), &b_left[i][j]);
//         }
//     }

//     let tau = bit_decompose(params, &minus_b_left);

//     for i in 0..m {
//         for h in 0..m {
//             let product = poly_mul(ring, &b_right[h], &tau[h][i]);
//             out[i] = poly_add(ring, &out[i], &product.to_vec());
//         }
//     }

//     out
// }

// /// `out = ct_left + ct_right`
// pub fn m_eval_add_x(
//     params: &Parameters,
//     ct_left: &Vec<Vec<u64>>,
//     ct_right: &Vec<Vec<u64>>,
// ) -> Vec<Vec<u64>> {
//     let ring = params.ring();
//     vec_vec_add(ring, ct_left, ct_right)
// }

// /// `out = ct_left * identity_matrix_scaled_by_x_right + ct_right * tau(-b_left)`
// pub fn m_eval_mul_x(
//     params: &Parameters,
//     ct_left: &Vec<Vec<u64>>,
//     b_left: &Vec<Vec<u64>>,
//     ct_right: &Vec<Vec<u64>>,
//     x_right: &u64,
// ) -> Vec<Vec<u64>> {
//     let ring = params.ring();
//     let m = *params.m();

//     let part_left = gen_identity_matrix_to_scalar(ring, m, *x_right);

//     let mut minus_b_left = vec![vec![ring.zero(); ring.ring_size()]; m];
//     for i in 0..m {
//         for j in 0..ring.ring_size() {
//             minus_b_left[i][j] = ring.sub(&ring.zero(), &b_left[i][j]);
//         }
//     }

//     let part_right = bit_decompose(params, &minus_b_left);

//     let out_left = vec_mat_mul(ring, ct_left, &part_left); // TODO: make it more efficient

//     let out_right = vec_mat_mul(ring, ct_right, &part_right);
//     let out = vec_vec_add(ring, &out_left, &out_right);

//     out
// }

// #[cfg(test)]
// mod tests {
//     use crate::{
//         ciphertext::Ciphertext,
//         eval::{m_eval_add, m_eval_add_x, m_eval_mul, m_eval_mul_x},
//         operations::{poly_add, poly_mul},
//         parameters::Parameters,
//         pub_key::PublicKey,
//     };
//     use phantom_zone_math::{
//         prelude::{ElemFrom, ModulusOps},
//         ring::RingOps,
//     };
//     use rand::{thread_rng, Rng};

//     #[test]
//     fn test_matrix_encoding_homomorphism_add_gate() {
//         let params = Parameters::new(12, 51, 7);
//         let pub_key = PublicKey::new(params);
//         let mut rng = thread_rng();
//         let ring = pub_key.params().ring();
//         let m = *pub_key.params().m();
//         let ell = *pub_key.params().ell();
//         let g = pub_key.params().g();
//         let mut x = (0..ell + 1)
//             .map(|_| rng.gen_range(0..2))
//             .collect::<Vec<_>>();
//         x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just
// for easier arithmetic during encoding

//         let ciphertext = Ciphertext::new(&pub_key, &x);

//         // perform add gate of b[1] and b[2]
//         let b_out = m_eval_add(pub_key.params(), &pub_key.b()[1], &pub_key.b()[2]);

//         // perform add gate of ct_inner[1] and ct_inner[2]
//         let ct_out = m_eval_add_x(
//             pub_key.params(),
//             &ciphertext.inner()[1],
//             &ciphertext.inner()[2],
//         );

//         let lhs = ct_out.clone();

//         // rhs = b_out + g * (x1 + x2)
//         let mut rhs = b_out.clone();
//         let mut fx = vec![ring.zero(); ring.ring_size()];
//         fx[0] = ring.elem_from(x[1] + x[2]);

//         for i in 0..m {
//             let gi_times_fx = poly_mul(ring, &g[i], &fx);
//             rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx);
//         }
//         for i in 0..m {
//             assert_eq!(lhs[i], rhs[i]);
//         }
//     }

//     #[test]
//     fn test_matrix_encoding_homomorphism_mul_gate() {
//         let params = Parameters::new(12, 51, 7);
//         let pub_key = PublicKey::new(params);
//         let mut rng = thread_rng();
//         let ring = pub_key.params().ring();
//         let m = *pub_key.params().m();
//         let ell = *pub_key.params().ell();
//         let g = pub_key.params().g();
//         let mut x = (0..ell + 1)
//             .map(|_| rng.gen_range(0..2))
//             .collect::<Vec<_>>();
//         x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just
// for easier arithmetic during encoding

//         let ciphertext = Ciphertext::new(&pub_key, &x);

//         // perform mul gate of b[1] and b[2]
//         let b_out = m_eval_mul(pub_key.params(), &pub_key.b()[1], &pub_key.b()[2]);

//         // perform mul gate of ct_inner[1] and ct_inner[2]
//         let ct_out = m_eval_mul_x(
//             pub_key.params(),
//             &ciphertext.inner()[1],
//             &pub_key.b()[1],
//             &ciphertext.inner()[2],
//             &x[2],
//         );

//         let lhs = ct_out.clone();

//         // rhs = b_out + g * (x1 * x2)
//         let mut rhs = b_out.clone();
//         let mut fx = vec![ring.zero(); ring.ring_size()];
//         fx[0] = ring.elem_from(x[1] * x[2]);

//         for i in 0..m {
//             let gi_times_fx = poly_mul(ring, &g[i], &fx);
//             rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx);
//         }
//         for i in 0..m {
//             assert_eq!(lhs[i], rhs[i]);
//         }
//     }
// }
