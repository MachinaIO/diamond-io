#[cfg(feature = "bgm")]
use super::bgm::Player;

use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggPublicKey, DigitsToInt,
    },
    io::{
        params::ObfuscationParams,
        utils::{build_final_digits_circuit, sample_public_key_by_id, PublicSampledData},
        Obfuscation,
    },
    poly::{
        enc::rlwe_encrypt,
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyElem, PolyMatrix, PolyParams,
    },
    utils::log_mem,
};
use itertools::Itertools;
use rand::{Rng, RngCore};
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use std::sync::Arc;

pub fn obfuscate<M, SU, SH, ST, R>(
    obf_params: ObfuscationParams<M>,
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
    #[cfg(feature = "bgm")]
    let player = Player::new();
    #[cfg(feature = "bgm")]
    player.play_music("bgm/obf_bgm1.mp3");

    let public_circuit = &obf_params.public_circuit;
    let dim = obf_params.params.ring_dimension() as usize;
    let log_base_q = obf_params.params.modulus_digits();
    let d = obf_params.d;
    let hash_key = rng.random::<[u8; 32]>();
    let sampler_uniform = SU::new();
    let sampler_trapdoor = ST::new(&obf_params.params, obf_params.trapdoor_sigma);
    let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(hash_key, d);
    let public_data = PublicSampledData::<SH>::sample(&obf_params, hash_key);
    log_mem("Sampled public data");
    let packed_input_size = public_data.packed_input_size;
    assert_eq!(public_circuit.num_input(), (2 * log_base_q) + (packed_input_size - 1));
    #[cfg(feature = "test")]
    let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
    #[cfg(not(feature = "test"))]
    let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();

    let pub_key_init =
        sample_public_key_by_id(&bgg_pubkey_sampler, &obf_params.params, 0, &reveal_plaintexts);
    log_mem("Sampled pub key init");

    let params = Arc::new(obf_params.params);
    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let s_bars = sampler_uniform.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
    log_mem("Sampled s_bars");
    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    log_mem("Sampled t_bar_matrix");
    let hardcoded_key_matrix = M::from_poly_vec_row(&params, vec![hardcoded_key.clone()]);
    log_mem("Sampled hardcoded_key_matrix");

    let a = public_data.a_rlwe_bar;

    let b = rlwe_encrypt(
        params.as_ref(),
        &sampler_uniform,
        &t_bar_matrix,
        &a,
        &hardcoded_key_matrix,
        obf_params.hardcoded_key_sigma,
    );

    log_mem("Generated RLWE ciphertext {a, b}");
    let m_b = (2 * (d + 1)) * (2 + log_base_q);
    let error = sampler_uniform.sample_uniform(
        &params,
        1,
        m_b,
        DistType::GaussDist { sigma: obf_params.p_sigma },
    );

    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bars,
        sampler_uniform,
        obf_params.encoding_sigma,
    );

    let s_init = &bgg_encode_sampler.secret_vec;

    let a_decomposed = a.entry(0, 0).decompose_base(params.as_ref());
    let b_decomposed = b.entry(0, 0).decompose_base(params.as_ref());

    log_mem("Decomposed RLWE ciphertext into {BaseDecompose(a), BaseDecompose(b)}");

    let minus_t_bar = -t_bar_matrix.entry(0, 0);

    let mut plaintexts = (0..obf_params.input_size.div_ceil(dim))
        .map(|_| M::P::const_zero(params.as_ref()))
        .collect_vec();
    plaintexts.push(minus_t_bar.clone());

    let encodings_init = bgg_encode_sampler.sample(&params, &pub_key_init, &plaintexts);
    log_mem("Sampled initial encodings");

    let (mut b_star_trapdoor_cur, mut b_star_cur) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
    log_mem("b star trapdoor init sampled");

    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init]);
        let s_b = s_connect * &b_star_cur;

        s_b + error
    };
    log_mem("Computed p_init");

    let identity_d_plus_1 = M::identity(params.as_ref(), d + 1, None);
    let level_width = obf_params.level_width; // number of bits to be inserted at each level
    assert_eq!(obf_params.input_size % level_width, 0);
    assert_eq!(dim % level_width, 0);
    let level_size = (1u64 << obf_params.level_width) as usize;
    let depth = obf_params.input_size / level_width; // number of levels necessary to encode the input
    let mut u_nums = Vec::with_capacity(level_size);
    for i in 0..level_size {
        let u_i = identity_d_plus_1.concat_diag(&[&public_data.rs[i]]);
        u_nums.push(u_i);
    }
    let u_star = {
        let zeros = M::zero(params.as_ref(), d + 1, 2 * (d + 1));
        let identities = identity_d_plus_1.concat_columns(&[&identity_d_plus_1]);
        zeros.concat_rows(&[&identities])
    };
    log_mem("Computed u_0, u_1, u_star");

    let (mut m_preimages, mut n_preimages, mut k_preimages) = (
        vec![Vec::with_capacity(level_size); depth],
        vec![Vec::with_capacity(level_size); depth],
        vec![Vec::with_capacity(level_size); depth],
    );

    #[cfg(feature = "test")]
    let mut bs: Vec<Vec<M>> = vec![vec![M::zero(params.as_ref(), 0, 0); level_size + 1]; depth + 1];

    #[cfg(feature = "test")]
    {
        bs[0][level_size] = b_star_cur.clone();
    }

    let mut pub_key_cur = pub_key_init;

    for level in 0..depth {
        let (b_star_trapdoor_level, b_star_level) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
        log_mem("Sampled b_star trapdoor for level");

        let pub_key_level =
            sample_public_key_by_id(&bgg_pubkey_sampler, &params, level + 1, &reveal_plaintexts);
        log_mem("Sampled pub key level");

        #[cfg(feature = "test")]
        {
            bs[level + 1][level_size] = b_star_level.clone();
        }

        // Precomputation for k_preimage that are not num dependent
        let lhs = -pub_key_cur[0].concat_matrix(&pub_key_cur[1..]);
        let inserted_poly_index = 1 + (level * level_width) / dim;
        let inserted_coeff_indices =
            (0..level_width).map(|i| (i + (level * level_width)) % dim).collect_vec();
        debug_assert_eq!(inserted_coeff_indices.len(), level_width);
        let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
        let mut coeffs = vec![zero_coeff; dim];

        for num in 0..level_size {
            #[cfg(feature = "bgm")]
            {
                player.play_music(format!("bgm/obf_bgm{}.mp3", (2 * level + num) % 3 + 2));
            }

            let (b_num_trapdoor_level, b_num_level) =
                sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
            log_mem("Sampled b trapdoor for level and num");

            #[cfg(feature = "test")]
            {
                bs[level + 1][num] = b_num_level.clone();
            }

            let m_preimage_num = sampler_trapdoor.preimage(
                &params,
                &b_star_trapdoor_cur,
                &b_star_cur,
                &(u_nums[num].clone() * &b_num_level),
            );

            log_mem("Computed m_preimage_num");

            m_preimages[level].push(m_preimage_num);

            let n_preimage_num = sampler_trapdoor.preimage(
                &params,
                &b_num_trapdoor_level,
                &b_num_level,
                &(u_star.clone() * &b_star_level.clone()),
            );
            log_mem("Computed n_preimage_num");

            n_preimages[level].push(n_preimage_num);

            let rg = &public_data.rgs[num];
            let top = lhs.mul_tensor_identity_decompose(rg, 1 + packed_input_size);
            log_mem("Computed top");
            // bit decompose num over level_width bits
            let num_bits: Vec<bool> = (0..level_width).map(|i| (num >> i) & 1 == 1).collect();
            debug_assert_eq!(num_bits.len(), level_width);
            for (i, coeff_idx) in inserted_coeff_indices.iter().enumerate() {
                let bit = num_bits[i];
                if bit {
                    coeffs[*coeff_idx] = <M::P as Poly>::Elem::one(&params.modulus());
                }
            }
            let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
            log_mem("Computed inserted_poly");
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
            log_mem("Computed inserted_poly_gadget");
            let bottom =
                pub_key_level[0].concat_matrix(&pub_key_level[1..]) - &inserted_poly_gadget;
            log_mem("Computed bottom");
            let k_target = top.concat_rows(&[&bottom]);
            log_mem("Computed k_target");
            let k_preimage_num =
                sampler_trapdoor.preimage(&params, &b_num_trapdoor_level, &b_num_level, &k_target);
            log_mem("Computed k_preimage_num");

            k_preimages[level].push(k_preimage_num);
        }

        b_star_trapdoor_cur = b_star_trapdoor_level;
        b_star_cur = b_star_level;
        pub_key_cur = pub_key_level;
    }
    #[cfg(feature = "bgm")]
    {
        player.play_music("bgm/obf_bgm5.mp3");
    }

    let final_preimage_target = {
        let final_circuit = build_final_digits_circuit::<M::P, BggPublicKey<M>>(
            &a_decomposed,
            &b_decomposed,
            public_circuit.clone(),
        );
        log_mem("Computed final_circuit");
        let eval_outputs = final_circuit.eval(params.as_ref(), &pub_key_cur[0], &pub_key_cur[1..]);
        log_mem("Evaluated outputs");
        debug_assert_eq!(eval_outputs.len(), log_base_q * packed_output_size);
        let output_ints = eval_outputs
            .par_chunks(log_base_q)
            .map(|digits| BggPublicKey::digits_to_int(digits, &params))
            .collect::<Vec<_>>();
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
        minus_t_bar,
        #[cfg(feature = "test")]
        bs,
        #[cfg(feature = "test")]
        hardcoded_key,
        #[cfg(feature = "test")]
        final_preimage_target,
    }
}
