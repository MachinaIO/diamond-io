use super::utils::*;
use super::{Obfuscation, ObfuscationParams};
use crate::bgg::sampler::{BGGEncodingSampler, BGGPublicKeySampler};
use crate::bgg::BggPublicKey;
use crate::poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams};
use rand::{Rng, RngCore};
use std::sync::Arc;

pub fn obfuscate<M, S, R>(
    obf_params: ObfuscationParams<M>,
    mut sampler: S,
    rng: &mut R,
) -> Obfuscation<M>
where
    M: PolyMatrix,
    S: PolyUniformSampler<M = M> + PolyHashSampler<[u8; 32], M = M> + PolyTrapdoorSampler<M = M>,
    R: RngCore,
{
    let public_circuit = &obf_params.public_circuit;
    debug_assert_eq!(public_circuit.num_input(), obf_params.input_size);
    let hash_key = rng.random::<[u8; 32]>();
    sampler.set_key(hash_key);
    let sampler = Arc::new(sampler);
    let params = Arc::new(obf_params.params.clone());
    let log_q = params.as_ref().modulus_bits();
    let dim = params.as_ref().ring_dimension() as usize;
    // let packed_input_size = obf_params.input_size.div_ceil(dim);
    let packed_output_size = public_circuit.num_output() / 2 / log_q;
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone());
    let public_data =
        PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler, packed_output_size);
    let packed_input_size = public_data.packed_input_size;

    let s_bar = sampler.sample_uniform(&params, 1, 1, DistType::BitDist).entry(1, 1).clone();
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bar,
        sampler.clone(),
        obf_params.error_gauss_sigma,
    );
    let s_init = &bgg_encode_sampler.secret_vec;
    let t_bar = sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    let minus_one_poly =
        M::from_poly_vec_row(params.as_ref(), vec![M::P::const_minus_one(params.as_ref())]);
    let t = t_bar.concat_columns(&[minus_one_poly]);

    let hardcoded_key = sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    let enc_hardcoded_key = {
        let e = sampler.sample_uniform(
            &params,
            1,
            1,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        t_bar.clone() * &public_data.a_rlwe_bar.clone() + &e + &hardcoded_key
    };
    let enc_hardcoded_key_polys = enc_hardcoded_key.decompose().get_column(0);

    let mut plaintexts = vec![];
    plaintexts.extend(enc_hardcoded_key_polys);
    let zero_plaintexts: Vec<M::P> =
        (0..packed_input_size).map(|_| M::P::const_zero(params.as_ref())).collect();
    plaintexts.extend(zero_plaintexts);
    plaintexts.push(t.entry(0, 0).clone());
    // let mut input_encoded_polys: Vec<<M as PolyMatrix>::P> =
    //     vec![<M::P as Poly>::const_one(params.as_ref())];
    // input_encoded_polys.extend(enc_hardcoded_key_polys);
    // input_encoded_polys.extend(zero_plaintexts);
    let encodings_init = bgg_encode_sampler.sample(&params, &public_data.pubkeys[0], &plaintexts);
    // let encode_fhe_key =
    //     bgg_encode_sampler.sample(&params, &public_data.pubkeys_fhe_key[0], &t.get_row(0), false);

    let mut bs = vec![];
    let mut b_trapdoors = vec![];
    for _ in 0..=obf_params.input_size {
        let (b_0_trapdoor, b_0) = sampler.trapdoor(&params);
        let (b_1_trapdoor, b_1) = sampler.trapdoor(&params);
        let (b_star_trapdoor, b_star) = sampler.trapdoor(&params);
        bs.push((b_0, b_1, b_star));
        b_trapdoors.push((b_0_trapdoor, b_1_trapdoor, b_star_trapdoor));
    }
    let m_b = 2 + log_q;
    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init.clone()]);
        let s_b = s_connect * &bs[0].2;
        let error = sampler.sample_uniform(
            &params,
            1,
            m_b,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        s_b + error
    };
    let identity_2 = M::identity(params.as_ref(), 2, None);
    let u_0 = identity_2.concat_diag(&[public_data.r_0.clone()]);
    let u_1 = identity_2.concat_diag(&[public_data.r_1.clone()]);
    let u_star = {
        let zeros = M::zero(params.as_ref(), 2, 4);
        let identities = identity_2.concat_columns(&[identity_2.clone()]);
        zeros.concat_rows(&[identities])
    };
    let gadget_2 = M::gadget_matrix(params.as_ref(), 2);

    let (mut m_preimages, mut n_preimages, mut k_preimages) = (vec![], vec![], vec![]);
    for idx in 0..obf_params.input_size {
        let (_, _, b_cur_star) = &bs[idx];
        let (b_next_0, b_next_1, b_next_star) = &bs[idx + 1];
        let (_, _, b_cur_star_trapdoor) = &b_trapdoors[idx];
        let (b_next_0_trapdoor, b_next_1_trapdoor, _) = &b_trapdoors[idx + 1];
        let m_0 = {
            let ub = u_0.clone() * b_next_0;
            sampler.preimage(params.as_ref(), b_cur_star_trapdoor, b_cur_star, &ub)
        };
        let m_1 = {
            let ub = u_1.clone() * b_next_1;
            sampler.preimage(params.as_ref(), b_cur_star_trapdoor, b_cur_star, &ub)
        };
        m_preimages.push((m_0, m_1));

        let ub_star = u_star.clone() * b_next_star;
        let n_0 = sampler.preimage(&params, b_next_0_trapdoor, b_next_0, &ub_star);
        let n_1 = sampler.preimage(&params, b_next_1_trapdoor, b_next_1, &ub_star);
        n_preimages.push((n_0, n_1));

        let mut ks = vec![];
        for bit in 0..1 {
            let t = &public_data.ts[bit];
            let top = public_data.pubkeys[idx][0].concat_matrix(&public_data.pubkeys[idx][1..]) * t;
            let inserted_poly_index = 1 + log_q + idx / dim;
            let inserted_coeff_index = idx % dim;
            let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            let mut coeffs = vec![zero_coeff; dim];
            coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus());
            let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
            let inserted_poly_gadget = {
                let zero = <M::P as Poly>::const_zero(params.as_ref());
                let mut polys = vec![];
                for _ in 0..(inserted_poly_index) {
                    polys.push(zero.clone());
                }
                polys.push(inserted_poly);
                for _ in (inserted_poly_index + 1)..(packed_input_size + 1) {
                    polys.push(zero.clone());
                }
                M::from_poly_vec_row(params.as_ref(), polys) * &gadget_2
            };
            let bottom = public_data.pubkeys[idx + 1][0]
                .concat_matrix(&public_data.pubkeys[idx + 1][1..])
                - &inserted_poly_gadget;
            let k_target = top.concat_rows(&[bottom]);
            let b_matrix = if bit == 0 { b_next_0 } else { b_next_1 };
            let trapdoor = if bit == 0 { b_next_0_trapdoor } else { b_next_1_trapdoor };
            let k = sampler.preimage(&params, trapdoor, b_matrix, &k_target);
            ks.push(k);

            // let (t_input, t_fhe_key) = if bit == 0 { &public_data.t_0 } else { &public_data.t_1 };
            // let at_input = public_data.pubkeys_input[idx][0]
            //     .concat_matrix(&public_data.pubkeys_input[idx][1..])
            //     * t_input;
            // let at_fhe_key = public_data.pubkeys_fhe_key[idx][0]
            //     .concat_matrix(&public_data.pubkeys_fhe_key[idx][1..])
            //     * t_fhe_key;
            // let former = at_input.concat_columns(&[at_fhe_key]);
            // let inserted_poly_index = 1 + log_q + idx / dim;
            // let inserted_coeff_index = idx % dim;
            // let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            // let mut coeffs = vec![zero_coeff; dim];
            // coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus());
            // let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
            // let inserted_poly_gadget = {
            //     let zero = <M::P as Poly>::const_zero(params.as_ref());
            //     let mut polys = vec![];
            //     for _ in 0..(inserted_poly_index) {
            //         polys.push(zero.clone());
            //     }
            //     polys.push(inserted_poly);
            //     for _ in inserted_poly_index + 1..packed_input_size {
            //         polys.push(zero.clone());
            //     }
            //     M::from_poly_vec_row(params.as_ref(), polys) * &gadget_2
            // };
            // let a_input_next = public_data.pubkeys_input[idx + 1][0]
            //     .concat_matrix(&public_data.pubkeys_input[idx + 1][1..])
            //     - &inserted_poly_gadget;
            // let latter = a_input_next.concat_columns(&[public_data.pubkeys_fhe_key[idx + 1][0]
            //     .concat_matrix(&public_data.pubkeys_fhe_key[idx + 1][1..])]);
            // let k_target = former.concat_rows(&[latter]);
            // let b_matrix = if bit == 0 { b_next_0 } else { b_next_1 };
            // let trapdoor = if bit == 0 { b_next_0_trapdoor } else { b_next_1_trapdoor };
            // let k = sampler.preimage(&params, trapdoor, b_matrix, &k_target);
            // ks.push(k);
        }
        k_preimages.push((ks[0].clone(), ks[1].clone()));
    }

    let final_circuit =
        build_final_step_circuit::<_, BggPublicKey<M>>(&params, public_circuit.clone());
    let final_preimage_target = {
        let one = public_data.pubkeys[obf_params.input_size][0].clone();
        let input = &public_data.pubkeys[obf_params.input_size][1..];
        let eval_outputs = final_circuit.eval_poly_circuit(params.as_ref(), one, input);
        let eval_outputs_matrix = eval_outputs[0].concat_matrix(&eval_outputs[1..]);
        debug_assert_eq!(eval_outputs_matrix.col_size(), packed_output_size);
        (eval_outputs_matrix + &public_data.a_prf).concat_rows(&[M::zero(
            params.as_ref(),
            2,
            packed_output_size,
        )])
    };

    // // here we support only inner product between the fhe secret key t and the input x
    // let mut ip_pubkey = None;
    // for idx in 0..packed_input_size {
    //     let muled = public_data.pubkeys_input[obf_params.input_size][idx].clone()
    //         * &public_data.pubkeys_fhe_key[obf_params.input_size][0];
    //     match ip_pubkey {
    //         None => {
    //             ip_pubkey = Some(muled);
    //         }
    //         Some(ip) => {
    //             ip_pubkey = Some(ip + muled);
    //         }
    //     }
    // }
    // let ip_pubkey = ip_pubkey.unwrap().matrix + &public_data.a_prf;
    // let final_preimage_target =
    //     ip_pubkey.concat_rows(&[M::zero(params.as_ref(), 2, ip_pubkey.col_size())]);
    let (_, _, b_final) = &bs[obf_params.input_size];
    let (_, _, b_final_trapdoor) = &b_trapdoors[obf_params.input_size];
    let final_preimage =
        sampler.preimage(&params, b_final_trapdoor, b_final, &final_preimage_target);

    Obfuscation {
        hash_key,
        encodings_init,
        p_init,
        m_preimages,
        n_preimages,
        k_preimages,
        final_preimage,
    }
}
