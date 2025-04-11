use super::{params::ObfuscationParams, utils::sample_public_key_by_idx, Obfuscation};
use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggPublicKey, BitToInt,
    },
    io::utils::{build_final_bits_circuit, PublicSampledData},
    poly::{
        enc::rlwe_encrypt,
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyElem, PolyMatrix, PolyParams,
    },
    utils::log_mem,
};
use itertools::Itertools;
use rand::{Rng, RngCore};
use std::sync::Arc;

pub fn obfuscate<M, SU, SH, ST, R>(
    obf_params: ObfuscationParams<M>,
    sampler_uniform: SU,
    mut sampler_hash: SH,
    sampler_trapdoor: ST,
    hardcoded_key: M::P,
    rng: &mut R,
) -> Obfuscation<M>
where
    M: PolyMatrix,
    SU: PolyUniformSampler<M = M>,
    SH: PolyHashSampler<[u8; 32], M = M>,
    ST: PolyTrapdoorSampler<M = M>,
    R: RngCore,
{
    let public_circuit = &obf_params.public_circuit;
    let dim = obf_params.params.ring_dimension() as usize;
    let log_q = obf_params.params.modulus_bits();
    let log_base_q = obf_params.params.modulus_digits();
    debug_assert_eq!(public_circuit.num_input(), (2 * log_q) + obf_params.input_size);
    let d = obf_params.d;
    let hash_key = rng.random::<[u8; 32]>();
    sampler_hash.set_key(hash_key);
    let sampler_uniform = Arc::new(sampler_uniform);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(Arc::new(sampler_hash), d);
    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    log_mem("Sampled public data");

    let packed_input_size = public_data.packed_input_size;
    #[cfg(feature = "test")]
    let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
    #[cfg(not(feature = "test"))]
    let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();

    let pub_key_init =
        sample_public_key_by_idx(&bgg_pubkey_sampler, &obf_params.params, 0, &reveal_plaintexts);
    log_mem("Sampled pub key init");

    let params = Arc::new(obf_params.params);
    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let s_bars = sampler_uniform.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
    log_mem("Sampled s_bars");
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bars,
        sampler_uniform.clone(),
        obf_params.encoding_sigma,
    );

    let s_init = &bgg_encode_sampler.secret_vec;
    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    log_mem("Sampled t_bar_matrix");

    let hardcoded_key_matrix = M::from_poly_vec_row(&params, vec![hardcoded_key.clone()]);
    log_mem("Sampled hardcoded_key_matrix");

    let a = public_data.a_rlwe_bar;

    let b = rlwe_encrypt(
        params.as_ref(),
        sampler_uniform.as_ref(),
        &t_bar_matrix,
        &a,
        &hardcoded_key_matrix,
        obf_params.hardcoded_key_sigma,
    );

    log_mem("Generated RLWE ciphertext {a, b}");

    let a_decomposed = a.entry(0, 0).decompose_bits(params.as_ref());
    let b_decomposed = b.entry(0, 0).decompose_bits(params.as_ref());

    log_mem("Decomposed RLWE ciphertext into {bits(a), bits(b)}");

    let t_bar = t_bar_matrix.entry(0, 0);

    let mut plaintexts = (0..obf_params.input_size.div_ceil(dim))
        .map(|_| M::P::const_zero(params.as_ref()))
        .collect_vec();
    plaintexts.push(t_bar.clone());

    let encodings_init = bgg_encode_sampler.sample(&params, &pub_key_init, &plaintexts);
    log_mem("Sampled initial encodings");

    let (mut b_star_trapdoor_cur, mut b_star_cur) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
    log_mem("b star trapdoor init sampled");

    let p_init = {
        let m_b = (2 * (d + 1)) * (2 + log_base_q);
        let s_connect = s_init.concat_columns(&[s_init]);
        let s_b = s_connect * &b_star_cur;
        let error = sampler_uniform.sample_uniform(
            &params,
            1,
            m_b,
            DistType::GaussDist { sigma: obf_params.p_sigma },
        );
        s_b + error
    };
    log_mem("Computed p_init");

    let identity_d_plus_1 = M::identity(params.as_ref(), d + 1, None);
    let u_0 = identity_d_plus_1.concat_diag(&[&public_data.r_0]);
    let u_1 = identity_d_plus_1.concat_diag(&[&public_data.r_1]);
    let u_bits = [u_0, u_1];
    let u_star = {
        let zeros = M::zero(params.as_ref(), d + 1, 2 * (d + 1));
        let identities = identity_d_plus_1.concat_columns(&[&identity_d_plus_1]);
        zeros.concat_rows(&[&identities])
    };
    log_mem("Computed u_0, u_1, u_star");

    let (mut m_preimages, mut n_preimages, mut k_preimages) = (
        vec![Vec::with_capacity(2); obf_params.input_size],
        vec![Vec::with_capacity(2); obf_params.input_size],
        vec![Vec::with_capacity(2); obf_params.input_size],
    );

    #[cfg(feature = "test")]
    let mut bs: Vec<Vec<M>> =
        vec![vec![M::zero(params.as_ref(), 0, 0); 3]; obf_params.input_size + 1];

    #[cfg(feature = "test")]
    {
        bs[0][2] = b_star_cur.clone();
    }

    let mut pub_key_cur = pub_key_init;

    for idx in 0..obf_params.input_size {
        let (b_star_trapdoor_idx, b_star_idx) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
        log_mem("Sampled b_star trapdoor for idx");

        let pub_key_idx =
            sample_public_key_by_idx(&bgg_pubkey_sampler, &params, idx + 1, &reveal_plaintexts);
        log_mem("Sampled pub key idx");

        #[cfg(feature = "test")]
        {
            bs[idx + 1][2] = b_star_idx.clone();
        }

        // Precomputation for k_preimage that are not bit dependent
        let lhs = -pub_key_cur[0].concat_matrix(&pub_key_cur[1..]);
        let inserted_poly_index = 1 + idx / dim;
        let inserted_coeff_index = idx % dim;
        let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
        let mut coeffs = vec![zero_coeff; dim];

        for bit in 0..=1 {
            let (b_bit_trapdoor_idx, b_bit_idx) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
            log_mem("Sampled b trapdoor for idx and bit");

            #[cfg(feature = "test")]
            {
                bs[idx + 1][bit] = b_bit_idx.clone();
            }

            let m_preimage_bit = sampler_trapdoor.preimage(
                &params,
                &b_star_trapdoor_cur,
                &b_star_cur,
                &(u_bits[bit].clone() * &b_bit_idx),
            );

            log_mem("Computed m_preimage_bit");

            m_preimages[idx].push(m_preimage_bit);

            let n_preimage_bit = sampler_trapdoor.preimage(
                &params,
                &b_bit_trapdoor_idx,
                &b_bit_idx,
                &(u_star.clone() * &b_star_idx.clone()),
            );
            log_mem("Computed n_preimage_bit");

            n_preimages[idx].push(n_preimage_bit);

            let rg = &public_data.rgs[bit];
            let top = lhs.mul_tensor_identity_decompose(rg, 1 + packed_input_size);
            if bit != 0 {
                coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus())
            };
            let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
            let inserted_poly_gadget = {
                let gadget_d_plus_1 = M::gadget_matrix(&params, d + 1);
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
            let bottom = pub_key_idx[0].concat_matrix(&pub_key_idx[1..]) - &inserted_poly_gadget;
            let k_target = top.concat_rows(&[&bottom]);
            let k_preimage_bit =
                sampler_trapdoor.preimage(&params, &b_bit_trapdoor_idx, &b_bit_idx, &k_target);
            log_mem("Computed k_preimage_bit");

            k_preimages[idx].push(k_preimage_bit);
        }

        b_star_trapdoor_cur = b_star_trapdoor_idx;
        b_star_cur = b_star_idx;
        pub_key_cur = pub_key_idx;
    }

    let final_preimage_target = {
        let final_circuit = build_final_bits_circuit::<M::P, BggPublicKey<M>>(
            &a_decomposed,
            &b_decomposed,
            public_circuit.clone(),
        );
        log_mem("Computed final_circuit");
        let eval_outputs = final_circuit.eval(params.as_ref(), &pub_key_cur[0], &pub_key_cur[1..]);
        log_mem("Evaluated outputs");
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
    log_mem("Computed final_preimage_target");

    let final_preimage = sampler_trapdoor.preimage(
        &params,
        &b_star_trapdoor_cur,
        &b_star_cur,
        &final_preimage_target,
    );
    log_mem("Sampled final_preimage");

    Obfuscation {
        hash_key,
        ct_b: b,
        encodings_init,
        p_init,
        m_preimages,
        n_preimages,
        k_preimages,
        final_preimage,
        #[cfg(feature = "test")]
        s_init: s_init.clone(),
        #[cfg(feature = "test")]
        t_bar,
        #[cfg(feature = "test")]
        bs,
        #[cfg(feature = "test")]
        hardcoded_key,
        #[cfg(feature = "test")]
        final_preimage_target,
    }
}
