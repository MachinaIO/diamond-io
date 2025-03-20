use super::{utils::*, Obfuscation, ObfuscationParams};
use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggPublicKey, BitToInt,
    },
    poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams},
};
use itertools::Itertools;
use rand::{Rng, RngCore};
use std::{ops::Mul, sync::Arc};

const TAG_BGG_PUBKEY_INPUT_PREFIX: &[u8] = b"BGG_PUBKEY_INPUT:";

pub fn obfuscate<M, SU, SH, ST, R>(
    obf_params: ObfuscationParams<M>,
    sampler_uniform: SU,
    mut sampler_hash: SH,
    sampler_trapdoor: ST,
    rng: &mut R,
) -> Obfuscation<M>
where
    M: PolyMatrix,
    SU: PolyUniformSampler<M = M>,
    SH: PolyHashSampler<[u8; 32], M = M>,
    ST: PolyTrapdoorSampler<M = M>,
    R: RngCore,
    for<'a> &'a M: Mul<&'a <M as PolyMatrix>::P, Output = M>,
    for<'a> &'a M: Mul<&'a M, Output = M>,
{
    let public_circuit = &obf_params.public_circuit;
    let dim = obf_params.params.ring_dimension() as usize;
    let log_q = obf_params.params.modulus_bits();
    debug_assert_eq!(public_circuit.num_input(), log_q + obf_params.input_size);
    let d = obf_params.d;
    let hash_key = rng.random::<[u8; 32]>();
    sampler_hash.set_key(hash_key);
    let sampler_uniform = Arc::new(sampler_uniform);
    let sampler_trapdoor = Arc::new(sampler_trapdoor);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(Arc::new(sampler_hash), d);
    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    let params = Arc::new(obf_params.params);
    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let s_bars =
        sampler_uniform.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0).clone();
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bars,
        sampler_uniform.clone(),
        obf_params.encoding_sigma,
    );
    let s_init = &bgg_encode_sampler.secret_vec;
    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist); // TODO: add sample for poly
    let hardcoded_key_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
    let enc_hardcoded_key = {
        let e = sampler_uniform.sample_uniform(
            &params,
            1,
            1,
            DistType::GaussDist { sigma: obf_params.hardcoded_key_sigma },
        );
        let scale = M::P::from_const(&params, &<M::P as Poly>::Elem::half_q(&params.modulus()));
        &t_bar_matrix * &public_data.a_rlwe_bar + &e - &(&hardcoded_key_matrix * &scale)
    };

    let enc_hardcoded_key_polys = enc_hardcoded_key.decompose().get_column(0); // TODO: first get column then decompose
    let t_bar = t_bar_matrix.entry(0, 0).clone();
    #[cfg(test)]
    let hardcoded_key = hardcoded_key_matrix.entry(0, 0).clone();

    let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
    #[cfg(not(test))]
    let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();

    let pub_keys_init = bgg_pubkey_sampler.sample(
        &params,
        &[TAG_BGG_PUBKEY_INPUT_PREFIX, &(0 as u64).to_le_bytes()].concat(),
        &reveal_plaintexts,
    );

    let mut plaintexts = (0..obf_params.input_size.div_ceil(dim))
        .map(|_| M::P::const_zero(params.as_ref()))
        .collect_vec();
    plaintexts.push(t_bar.clone());
    let encodings_init = bgg_encode_sampler.sample(&params, &pub_keys_init, &plaintexts);

    let (b_star_trapdoor_init, b_star_init) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));

    let m_b = (2 * (d + 1)) * (2 + log_q);
    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init]);
        let s_b = s_connect * &b_star_init;
        let error = sampler_uniform.sample_uniform(
            &params,
            1,
            m_b,
            DistType::GaussDist { sigma: obf_params.p_sigma },
        );
        s_b + error
    };
    let identity_d_plus_1 = M::identity(params.as_ref(), d + 1, None);
    let u_0 = identity_d_plus_1.concat_diag(&[&public_data.r_0]);
    let u_1 = identity_d_plus_1.concat_diag(&[&public_data.r_1]);
    let u_bits = [u_0, u_1];
    let u_star = {
        let zeros = M::zero(params.as_ref(), d + 1, 2 * (d + 1));
        let identities = identity_d_plus_1.concat_columns(&[&identity_d_plus_1]);
        zeros.concat_rows(&[&identities])
    };
    let gadget_d_plus_1 = M::gadget_matrix(params.as_ref(), d + 1); // TODO: fetch it from public data

    let mut b_star_trapdoor_cur = b_star_trapdoor_init;
    let mut b_star_cur = b_star_init;
    let mut pub_keys_cur = pub_keys_init;

    let mut m_preimages = Vec::<Vec<M>>::with_capacity(obf_params.input_size);
    let mut n_preimages = Vec::<Vec<M>>::with_capacity(obf_params.input_size);
    let mut k_preimages = Vec::<Vec<M>>::with_capacity(obf_params.input_size);

    for idx in 0..obf_params.input_size {
        // Sample pub keys
        let pub_keys_idx = bgg_pubkey_sampler.sample(
            &params,
            &[TAG_BGG_PUBKEY_INPUT_PREFIX, &(idx as u64).to_le_bytes()].concat(),
            &reveal_plaintexts,
        );

        // Sample B and B trapdoor
        let (b_star_trapdoor_idx, b_star_idx) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));

        // Precomputation for k_preimage that are not bit dependent
        let lhs = -pub_keys_cur[0].concat_matrix(&pub_keys_cur[1..]);
        let inserted_poly_index = 1 + idx / dim;
        let inserted_coeff_index = idx % dim;

        let mut mp = Vec::with_capacity(2);
        let mut np = Vec::with_capacity(2);
        let mut kp = Vec::with_capacity(2);

        for bit in 0..=1 {
            let (b_bit_trapdoor_idx, b_bit_idx) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
            let m_preimage_bit = sampler_trapdoor.preimage(
                &params,
                &b_star_trapdoor_cur,
                &b_star_cur,
                &(&u_bits[bit] * &b_bit_idx),
            );

            mp.push(m_preimage_bit);

            let n_preimage_bit = sampler_trapdoor.preimage(
                &params,
                &b_bit_trapdoor_idx,
                &b_bit_idx,
                &(&u_star * &b_star_idx.clone()),
            );

            np.push(n_preimage_bit);

            // compute k_preimage
            let rg = &public_data.rgs[bit];
            let top = lhs.mul_tensor_identity_decompose(rg, 1 + packed_input_size);

            // TO DO: how is the following part bit dependent?
            let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            let mut coeffs = vec![zero_coeff; dim];
            if bit != 0 {
                coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus())
            };
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
                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_d_plus_1)
            };
            let bottom = pub_keys_idx[0].concat_matrix(&pub_keys_idx[1..]) - &inserted_poly_gadget;
            let k_target = top.concat_rows(&[&bottom]);
            let k_preimage_bit =
                sampler_trapdoor.preimage(&params, &b_bit_trapdoor_idx, &b_bit_idx, &k_target);

            kp.push(k_preimage_bit);
        }

        m_preimages.push(mp);
        n_preimages.push(np);
        k_preimages.push(kp);

        b_star_trapdoor_cur = b_star_trapdoor_idx;
        b_star_cur = b_star_idx;
        pub_keys_cur = pub_keys_idx;
    }

    let a_decomposed_polys = public_data.a_rlwe_bar.decompose().get_column(0); // TODO: first get column then decompose
    let final_circuit = build_final_bits_circuit::<M::P, BggPublicKey<M>>(
        &a_decomposed_polys,
        &enc_hardcoded_key_polys,
        public_circuit.clone(),
    );
    let final_preimage_target = {
        let one = pub_keys_cur[0].clone();
        let input = &pub_keys_cur[1..];
        let eval_outputs = final_circuit.eval(params.as_ref(), one, input);
        assert_eq!(eval_outputs.len(), log_q * packed_output_size);
        let output_ints = eval_outputs
            .chunks(log_q)
            .map(|bits| BggPublicKey::bits_to_int(bits, &params))
            .collect_vec();
        let eval_outputs_matrix = output_ints[0].concat_matrix(&output_ints[1..]);
        debug_assert_eq!(eval_outputs_matrix.col_size(), packed_output_size);
        (eval_outputs_matrix + public_data.a_prf).concat_rows(&[&M::zero(
            params.as_ref(),
            d + 1,
            packed_output_size,
        )])
    };
    let final_preimage = sampler_trapdoor.preimage(
        &params,
        &b_star_trapdoor_cur,
        &b_star_cur,
        &final_preimage_target,
    );
    Obfuscation {
        hash_key,
        enc_hardcoded_key,
        encodings_init,
        p_init,
        m_preimages,
        n_preimages,
        k_preimages,
        final_preimage,
        #[cfg(test)]
        s_init: s_init.clone(),
        #[cfg(test)]
        t_bar: t_bar.clone(),
        // #[cfg(test)]
        // bs,
        #[cfg(test)]
        hardcoded_key: hardcoded_key.clone(),
        #[cfg(test)]
        final_preimage_target,
    }
}
