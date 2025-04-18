#[cfg(feature = "bgm")]
use super::bgm::Player;

use super::{params::ObfuscationParams, utils::*, Obfuscation};
use crate::{
    bgg::{sampler::BGGPublicKeySampler, BggEncoding, DigitsToInt},
    parallel_iter,
    poly::{sampler::*, Poly, PolyElem, PolyMatrix, PolyParams},
    utils::log_mem,
};
use itertools::Itertools;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use std::{path::Path, sync::Arc};

impl<M> Obfuscation<M>
where
    M: PolyMatrix,
{
    pub fn read_dir<P: AsRef<Path> + Send + Sync>(
        obf_params: &ObfuscationParams<M>,
        dir_path: P,
    ) -> Self {
        let dir_path = dir_path.as_ref().to_path_buf();
        let b = M::read_from_files(&obf_params.params, 1, 1, &dir_path, "b");

        let dim = obf_params.params.ring_dimension() as usize;
        let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
        let d = obf_params.d;
        let d1 = d + 1;
        let log_base_q = obf_params.params.modulus_digits();
        #[cfg(feature = "test")]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
        #[cfg(not(feature = "test"))]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();
        let params = Arc::new(obf_params.params.clone());
        let encodings_init = parallel_iter!(0..packed_input_size)
            .map(|idx| {
                BggEncoding::read_from_files(
                    params.as_ref(),
                    d1,
                    log_base_q,
                    &dir_path,
                    &format!("encoding_init_{idx}"),
                    reveal_plaintexts[idx],
                )
            })
            .collect::<Vec<_>>();

        let m_b = (2 * d1) * (2 + log_base_q);
        let p_init = M::read_from_files(&obf_params.params, 1, m_b, &dir_path, "p_init");

        let level_size = (1u64 << obf_params.level_width) as usize;
        let depth = obf_params.input_size / obf_params.level_width;
        let m_preimages = parallel_iter!(0..depth)
            .map(|level| {
                parallel_iter!(0..level_size)
                    .map(|num| {
                        M::read_from_files(
                            params.as_ref(),
                            m_b,
                            m_b,
                            &dir_path,
                            &format!("m_preimage_{level}_{num}"),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let n_preimages = parallel_iter!(0..depth)
            .map(|level| {
                parallel_iter!(0..level_size)
                    .map(|num| {
                        M::read_from_files(
                            params.as_ref(),
                            m_b,
                            m_b,
                            &dir_path,
                            &format!("n_preimage_{level}_{num}"),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let k_columns = (1 + packed_input_size) * d1 * log_base_q;
        let k_preimages = parallel_iter!(0..depth)
            .map(|level| {
                parallel_iter!(0..level_size)
                    .map(|num| {
                        M::read_from_files(
                            params.as_ref(),
                            m_b,
                            k_columns,
                            &dir_path,
                            &format!("k_preimage_{level}_{num}"),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let packed_output_size = obf_params.public_circuit.num_output() / (2 * log_base_q);
        let final_preimage = M::read_from_files(
            &obf_params.params,
            m_b,
            packed_output_size,
            &dir_path,
            "final_preimage",
        );
        let hash_key = {
            let mut path = dir_path.clone();
            path.push("hash_key");
            let bytes = std::fs::read(&path).expect("Failed to read hash key file");
            let mut hash_key = [0u8; 32];
            hash_key.copy_from_slice(&bytes);
            hash_key
        };
        #[cfg(feature = "test")]
        let s_init = M::read_from_files(&obf_params.params, 1, d1, &dir_path, "s_init");

        #[cfg(feature = "test")]
        let minus_t_bar = <<M as PolyMatrix>::P as Poly>::read_from_file(
            &obf_params.params,
            &dir_path,
            "minus_t_bar",
        );

        #[cfg(feature = "test")]
        let hardcoded_key = <<M as PolyMatrix>::P as Poly>::read_from_file(
            &obf_params.params,
            &dir_path,
            "hardcoded_key",
        );

        #[cfg(feature = "test")]
        let bs = parallel_iter!(0..depth + 1)
            .map(|level| {
                let mut b_nums = parallel_iter!(0..level_size)
                    .map(|num| {
                        M::read_from_files(
                            params.as_ref(),
                            2 * d1,
                            m_b,
                            &dir_path,
                            &format!("b_{level}_{num}"),
                        )
                    })
                    .collect::<Vec<_>>();
                let b_star = M::read_from_files(
                    params.as_ref(),
                    2 * d1,
                    m_b,
                    &dir_path,
                    &format!("b_star_{level}"),
                );
                b_nums.push(b_star);
                b_nums
            })
            .collect::<Vec<_>>();

        Self {
            b,
            encodings_init,
            p_init,
            m_preimages,
            n_preimages,
            k_preimages,
            hash_key,
            final_preimage,
            #[cfg(feature = "test")]
            s_init,
            #[cfg(feature = "test")]
            minus_t_bar,
            #[cfg(feature = "test")]
            bs,
            #[cfg(feature = "test")]
            hardcoded_key,
        }
    }

    pub fn eval<SH, ST>(
        &self,
        obf_params: &ObfuscationParams<M>,
        mut sampler_hash: SH,
        inputs: &[bool],
    ) -> Vec<bool>
    where
        SH: PolyHashSampler<[u8; 32], M = M>,
        ST: PolyTrapdoorSampler<M = M>,
    {
        #[cfg(feature = "bgm")]
        let player = Player::new();

        #[cfg(feature = "bgm")]
        {
            player.play_music("bgm/eval_bgm1.mp3");
        }

        sampler_hash.set_key(self.hash_key);
        let params = Arc::new(obf_params.params.clone());
        let d = obf_params.d;
        let d1 = d + 1;
        let sampler = Arc::new(sampler_hash);
        debug_assert_eq!(inputs.len(), obf_params.input_size);
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone(), d);
        let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
        log_mem("Sampled public data");
        let packed_input_size = public_data.packed_input_size;
        let packed_output_size = public_data.packed_output_size;
        let (mut ps, mut encodings) = (vec![], vec![]);
        ps.push(self.p_init.clone());
        encodings.push(self.encodings_init.clone());

        let level_width = obf_params.level_width;
        #[cfg(feature = "test")]
        let level_size = (1u64 << obf_params.level_width) as usize;
        assert!(inputs.len() % level_width == 0);
        let depth = obf_params.input_size / level_width;

        #[cfg(feature = "test")]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
        #[cfg(not(feature = "test"))]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();
        let pub_key_init =
            sample_public_key_by_id(&bgg_pubkey_sampler, &obf_params.params, 0, &reveal_plaintexts);
        log_mem("Sampled pub_key_init");

        let mut pub_key_cur = pub_key_init;

        #[cfg(feature = "test")]
        if obf_params.encoding_sigma == 0.0 &&
            obf_params.hardcoded_key_sigma == 0.0 &&
            obf_params.p_sigma == 0.0
        {
            let expected_p_init = {
                let s_connect = self.s_init.concat_columns(&[&self.s_init]);
                s_connect * &self.bs[0][level_size]
            };
            assert_eq!(self.p_init, expected_p_init);

            let zero = <M::P as Poly>::const_zero(&params);
            let one = <M::P as Poly>::const_one(&params);
            let inserted_poly_gadget = {
                let mut polys = vec![];
                polys.push(one.clone());
                for _ in 0..(packed_input_size - 1) {
                    polys.push(zero.clone());
                }
                polys.push(self.minus_t_bar.clone());
                let gadget_d1 = M::gadget_matrix(&params, d1);
                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_d1)
            };
            let expected_encoding_init = self.s_init.clone() *
                &(pub_key_cur[0].concat_matrix(&pub_key_cur[1..]) - inserted_poly_gadget);
            assert_eq!(encodings[0][0].concat_vector(&encodings[0][1..]), expected_encoding_init);
        }
        let log_base_q = params.as_ref().modulus_digits();
        let dim = params.as_ref().ring_dimension() as usize;
        let nums: Vec<u64> = inputs
            .chunks(level_width)
            .map(|chunk| {
                chunk.iter().enumerate().fold(0u64, |acc, (i, &bit)| acc + ((bit as u64) << i))
            })
            .collect();
        debug_assert_eq!(nums.len(), depth);
        for (level, num) in nums.iter().enumerate() {
            let m = &self.m_preimages[level][*num as usize];
            let q = ps[level].clone() * m;
            log_mem(format!("q at {} computed", level));
            let n = &self.n_preimages[level][*num as usize];
            let p = q.clone() * n;
            log_mem(format!("p at {} computed", level));
            let k = &self.k_preimages[level][*num as usize];
            let v = q.clone() * k;
            log_mem(format!("v at {} computed", level));
            let new_encode_vec = {
                let rg = &public_data.rgs[*num as usize];
                let encode_vec = encodings[level][0].concat_vector(&encodings[level][1..]);
                let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
                encode_vec.mul_tensor_identity_decompose(rg, packed_input_size + 1) + v
            };
            log_mem(format!("new_encode_vec at {} computed", level));
            let mut new_encodings = vec![];
            let inserted_poly_index = 1 + (level * level_width) / dim;
            let pub_key_level = sample_public_key_by_id(
                &bgg_pubkey_sampler,
                &obf_params.params,
                level + 1,
                &reveal_plaintexts,
            );
            log_mem(format!("pub_key_level at {} computed", level));
            for (j, encode) in encodings[level].iter().enumerate() {
                let m = d1 * log_base_q;
                let new_vec = new_encode_vec.slice_columns(j * m, (j + 1) * m);
                log_mem(format!("new_vec at {}, {} computed", level, j));
                let plaintext = if j == inserted_poly_index {
                    let inserted_coeff_indices =
                        (0..level_width).map(|i| (i + (level * level_width)) % dim).collect_vec();
                    let mut coeffs = encode.plaintext.as_ref().unwrap().coeffs().clone();
                    let num_bits: Vec<bool> =
                        (0..level_width).map(|i| (num >> i) & 1 == 1).collect();
                    debug_assert_eq!(num_bits.len(), level_width);
                    for (i, coeff_idx) in inserted_coeff_indices.iter().enumerate() {
                        let bit = num_bits[i];
                        if bit {
                            coeffs[*coeff_idx] = <M::P as Poly>::Elem::one(&params.modulus());
                        }
                    }
                    Some(M::P::from_coeffs(params.as_ref(), &coeffs))
                } else {
                    encode.plaintext.clone()
                };
                log_mem(format!("plaintext at {}, {} computed", level, j));
                let new_encode: BggEncoding<M> =
                    BggEncoding::new(new_vec, pub_key_level[j].clone(), plaintext);
                log_mem(format!("new_encode at {}, {} computed", level, j));
                new_encodings.push(new_encode);
            }
            ps.push(p.clone());
            pub_key_cur = pub_key_level;
            encodings.push(new_encodings);
            #[cfg(feature = "test")]
            if obf_params.encoding_sigma == 0.0 &&
                obf_params.hardcoded_key_sigma == 0.0 &&
                obf_params.p_sigma == 0.0
            {
                let mut cur_s = self.s_init.clone();
                for prev_num in nums[0..level].iter() {
                    let r = public_data.rs[*prev_num as usize].clone();
                    cur_s = cur_s * r;
                }
                let new_s = cur_s.clone() * &public_data.rs[*num as usize];
                let b_next_bit = self.bs[level + 1][*num as usize].clone();
                let expected_q = cur_s.concat_columns(&[&new_s]) * &b_next_bit;
                assert_eq!(q, expected_q);
                let expected_p = new_s.concat_columns(&[&new_s]) * &self.bs[level + 1][level_size];
                assert_eq!(p, expected_p);
                let expcted_new_encode = {
                    let dim = params.ring_dimension() as usize;
                    let one = <M::P as Poly>::const_one(&params);
                    let gadget_d1 = M::gadget_matrix(&params, d1);
                    let inserted_poly_gadget = {
                        let mut polys = vec![];
                        polys.push(one.clone());
                        let mut coeffs = vec![];
                        for bit in inputs[0..(level_width * (level + 1))].iter() {
                            if *bit {
                                coeffs.push(<M::P as Poly>::Elem::one(&params.modulus()));
                            } else {
                                coeffs.push(<M::P as Poly>::Elem::zero(&params.modulus()));
                            }
                        }
                        for _ in 0..(obf_params.input_size - level_width * (level + 1)) {
                            coeffs.push(<M::P as Poly>::Elem::zero(&params.modulus()));
                        }
                        let input_polys = coeffs
                            .chunks(dim)
                            .map(|coeffs| M::P::from_coeffs(&params, coeffs))
                            .collect_vec();
                        polys.extend(input_polys);
                        polys.push(self.minus_t_bar.clone());
                        M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_d1)
                    };
                    let pubkey = pub_key_cur[0].concat_matrix(&pub_key_cur[1..]);
                    new_s * (pubkey - inserted_poly_gadget)
                };
                assert_eq!(new_encode_vec, expcted_new_encode);
            }
        }

        #[cfg(feature = "bgm")]
        {
            player.play_music("bgm/eval_bgm2.mp3");
        }

        let a_decomposed = public_data.a_rlwe_bar.entry(0, 0).decompose_base(params.as_ref());
        let b_decomposed = &self.b.entry(0, 0).decompose_base(params.as_ref());
        log_mem("a,b decomposed");
        let final_circuit = build_final_digits_circuit::<M::P, BggEncoding<M>>(
            &a_decomposed,
            b_decomposed,
            obf_params.public_circuit.clone(),
        );
        log_mem("final_circuit built");
        let last_input_encodings = encodings.last().unwrap();
        let output_encodings = final_circuit.eval::<BggEncoding<M>>(
            params.as_ref(),
            &last_input_encodings[0],
            &last_input_encodings[1..],
        );
        log_mem("final_circuit evaluated");
        let output_encoding_ints = output_encodings
            .par_chunks(log_base_q)
            .map(|digits| BggEncoding::digits_to_int(digits, &params))
            .collect::<Vec<_>>();
        let output_encodings_vec =
            output_encoding_ints[0].concat_vector(&output_encoding_ints[1..]);
        log_mem("final_circuit evaluated and recomposed");
        let final_preimage = &self.final_preimage;
        let final_v = ps.last().unwrap().clone() * final_preimage;
        log_mem("final_v computed");
        let z = output_encodings_vec.clone() - final_v.clone();
        log_mem("z computed");
        debug_assert_eq!(z.size(), (1, packed_output_size));
        #[cfg(feature = "test")]
        if obf_params.encoding_sigma == 0.0 &&
            obf_params.hardcoded_key_sigma == 0.0 &&
            obf_params.p_sigma == 0.0
        {
            let mut last_s = self.s_init.clone();
            for num in nums.iter() {
                let r = public_data.rs[*num as usize].clone();
                last_s = last_s * r;
            }
            {
                let expected = last_s *
                    (output_encoding_ints[0].pubkey.matrix.clone() -
                        M::unit_column_vector(params.as_ref(), d1, d1 - 1) *
                            output_encoding_ints[0].plaintext.clone().unwrap());
                assert_eq!(output_encoding_ints[0].vector, expected);
            }
            assert_eq!(z.size(), (1, packed_output_size));
            if inputs[0] {
                assert_eq!(
                    output_encoding_ints[0]
                        .plaintext
                        .clone()
                        .unwrap()
                        .extract_bits_with_threshold(&params),
                    self.hardcoded_key.to_bool_vec()
                );
            }
        }
        z.get_row(0).into_iter().flat_map(|p| p.extract_bits_with_threshold(&params)).collect_vec()
    }
}
