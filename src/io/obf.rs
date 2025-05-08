#[cfg(feature = "bgm")]
use super::bgm::Player;
use crate::{
    bgg::{sampler::BGGPublicKeySampler, BggPublicKey, DigitsToInt},
    io::{
        params::ObfuscationParams,
        utils::{
            build_final_digits_circuit, build_u_mask_multi, sample_public_key_by_id,
            PublicSampledData,
        },
    },
    poly::{
        enc::rlwe_encrypt,
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
    utils::log_mem,
};
use futures::future::join_all;
use itertools::Itertools;
use rand::{Rng, RngCore};
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use std::{path::Path, sync::Arc};
use tokio::runtime::Handle;

pub async fn obfuscate<M, SU, SH, ST, R, P>(
    obf_params: ObfuscationParams<M>,
    hardcoded_key: M::P,
    rng: &mut R,
    dir_path: P,
) where
    M: PolyMatrix + 'static,
    SU: PolyUniformSampler<M = M>,
    SH: PolyHashSampler<[u8; 32], M = M>,
    ST: PolyTrapdoorSampler<M = M>,
    R: RngCore,
    P: AsRef<Path>,
{
    #[cfg(feature = "bgm")]
    let player = Player::new();
    #[cfg(feature = "bgm")]
    player.play_music("bgm/obf_bgm1.mp3");

    let mut handles: Vec<tokio::task::JoinHandle<()>> = Vec::new();
    let dir_path = dir_path.as_ref().to_path_buf();
    if !dir_path.exists() {
        std::fs::create_dir_all(&dir_path).expect("Failed to create directory");
    }

    let hash_key = rng.random::<[u8; 32]>();
    let public_data = PublicSampledData::<SH>::sample(&obf_params, hash_key);
    log_mem("Sampled public data");
    let params = Arc::new(obf_params.params);
    let public_circuit = obf_params.public_circuit;
    let packed_input_size = public_data.packed_input_size;
    let log_base_q = params.modulus_digits();
    assert_eq!(public_circuit.num_input(), (2 * log_base_q) + (packed_input_size - 1));
    let n = params.ring_dimension() as usize;
    let d = obf_params.d;
    let sampler_uniform = SU::new();
    let sampler_trapdoor = ST::new(&params, obf_params.trapdoor_sigma);
    let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(hash_key, d);
    let m_b = (1 + packed_input_size) * (d + 1) * (2 + log_base_q);
    let packed_output_size = public_data.packed_output_size;

    /*
    =============================================================================
    We reveal all input slots so the evaluator can multiply encodings,
    but we keep the secret-key slot (last slot) hidden unless we're in debug mode.
    In the paper, the last slot corresponds to the FHE secret-key vector t.

    Sample the initial public key (level 0) with our reveal flags.

    The length of encodings_init is (1 + packed_input_size):
       - 1 the encoding of constant 1
       - packed_input_size for the packed evaluator inputs and for FHE secret key t
    =============================================================================
    */

    // Sample BGG+ encoding secret key
    let s_bars = sampler_uniform.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
    log_mem("Sampled s_bars");
    // Sample FHE secret key t
    let t_bar = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
    log_mem("Sampled t_bar");
    let minus_t_bar = -t_bar.entry(0, 0);
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_poly(minus_t_bar.clone(), &dir_path, "minus_t_bar"));
    // This plaintexts is (1, 0_L, t), total length is L + 2, packed_input_size - 1 is L
    let one = M::P::const_one(&params);
    let mut plaintexts = vec![one];
    plaintexts.extend((0..packed_input_size - 1).map(|_| M::P::const_zero(&params)).collect_vec());
    plaintexts.push(minus_t_bar);

    let mut reveal_plaintexts = vec![true; plaintexts.len()];
    // We reveal FHE secret key, only in the debug mode
    reveal_plaintexts[packed_input_size - 1] = cfg!(feature = "debug");
    let s_init = {
        let minus_one_poly = <SU::M as PolyMatrix>::P::const_minus_one(&params);
        let mut secrets = s_bars.to_vec();
        secrets.push(minus_one_poly);
        SU::M::from_poly_vec_row(&params, secrets)
    };
    log_mem(format!("s_init ({},{}) 1x(d+1)", s_init.row_size(), s_init.col_size()));

    /*
    =============================================================================
    Pre‐loop initialization:

    1) Sample input‐dependent random basis B_* (trapdoor & public matrix).
    2) Compute initial secret key p_init = ((1, 0_L, t) ⊗ s_init)·B_* + error.
    where 0_L is dummy mask for input space and t is fhe secret key,
    3) Create I_{d+1}, derive level_width, level_size, and depth.
    4) For each i in 0..level_size,
        For each j in 0..depth,
            U_{j,0} = I_{1 + packed_input_size}
            U_(j,i) = Matrix to update input bit to 1

    These values (B_*, p_init, I_{d+1}, U_{j, i}) are all set up
    before entering the main preimage generation loop.
    =============================================================================
    */

    // Sample input dependent random matrix B_*
    let (mut b_star_trapdoor_cur, mut b_star_cur) =
        sampler_trapdoor.trapdoor(&params, (1 + packed_input_size) * (d + 1));
    log_mem(format!(
        "b star epsilon ({},{}) (d+1)xm_B and trapdoor epsilon sampled",
        b_star_cur.row_size(),
        b_star_cur.col_size()
    ));

    // Compute p_init
    let encoded_bits = M::from_poly_vec_row(&params, plaintexts);
    log_mem(format!("encoded_bits ({},{})", encoded_bits.row_size(), encoded_bits.col_size()));
    let s_connect = encoded_bits.tensor(&s_init);
    log_mem(format!("s_connect ({},{})", s_connect.row_size(), s_connect.col_size()));
    let s_b = s_connect * &b_star_cur;
    let p_init_error = sampler_uniform.sample_uniform(
        &params,
        1,
        m_b,
        DistType::GaussDist { sigma: obf_params.p_sigma },
    );
    let p_init = s_b + p_init_error;
    log_mem(format!("Computed p_epsilon ({},{})", p_init.row_size(), p_init.col_size()));
    handles.push(store_and_drop_matrix(p_init, &dir_path, "p_init"));
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_matrix(s_init, &dir_path, "s_init"));
    let identity_1_plus_packed_input_size =
        M::identity(params.as_ref(), 1 + packed_input_size, None);
    log_mem(format!(
        "Computed identity_1_plus_packed_input_size {} {}",
        identity_1_plus_packed_input_size.row_size(),
        identity_1_plus_packed_input_size.col_size()
    ));
    // number of bits to be inserted at each level
    let level_width = obf_params.level_width;
    assert_eq!(obf_params.input_size % level_width, 0);
    // otherwise we need >1 polynomial to insert the bits for each level
    assert!(level_width <= n);
    // otherwise we get to a point in which the inserted bits have to be split between two
    // polynomials
    if obf_params.input_size > n {
        assert_eq!(n % level_width, 0);
    }
    let level_size = (1u64 << obf_params.level_width) as usize;
    // number of levels necessary to encode the input
    let depth = obf_params.input_size / level_width;
    let mut u_nums = Vec::with_capacity(depth);

    for j in 0..depth {
        let mut masks_at_level = Vec::with_capacity(level_size);
        masks_at_level.push(identity_1_plus_packed_input_size.clone());
        if level_size > 0 {
            for i in 1..level_size {
                masks_at_level.push(build_u_mask_multi::<M>(
                    &params,
                    packed_input_size,
                    level_width,
                    j,
                    i,
                ));
            }
        }

        u_nums.push(masks_at_level);
    }

    log_mem(format!("Computed u_0, .. u_{}", depth));
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_matrix(b_star_cur.clone(), &dir_path, "b_star_0"));

    /*
    Trapdoor preimage generation for the input insertion step.
    For each depth, sample K preimage at the corresponding level size.
    Level starts from 1 because we already have p_init at level 0.
    */
    let one_identity = M::identity(params.as_ref(), 1, None);
    let hash_sampler = SH::new();
    for level in 1..(depth + 1) {
        let (b_star_trapdoor_level, b_star_level) =
            sampler_trapdoor.trapdoor(&params, (1 + packed_input_size) * (d + 1));
        log_mem(format!(
            "Sampled b_star ({},{}) and trapdoor for level {}",
            b_star_level.row_size(),
            b_star_level.col_size(),
            level
        ));

        #[cfg(feature = "debug")]
        handles.push(store_and_drop_matrix(
            b_star_level.clone(),
            &dir_path,
            &format!("b_star_{}", level),
        ));

        for num in 0..level_size {
            let mut handles_per_level: Vec<tokio::task::JoinHandle<()>> = Vec::new();
            #[cfg(feature = "bgm")]
            {
                player.play_music(format!("bgm/obf_bgm{}.mp3", (2 * level + num) % 3 + 2));
            }

            let tag = format!("S_{}_{}", level, num).into_bytes();
            let s_i_bar =
                hash_sampler.sample_hash(&params, hash_key, &tag, d, d, DistType::BitDist);
            let s_i_num = s_i_bar.concat_diag(&[&one_identity]);

            #[cfg(feature = "debug")]
            handles.push(store_and_drop_matrix(
                s_i_num.clone(),
                &dir_path,
                &format!("s_{}_{}", level, num),
            ));

            log_mem(format!(
                "Computed S ({},{}) (d+1)x(d+1)",
                s_i_num.row_size(),
                s_i_num.col_size()
            ));
            let u = &u_nums[level - 1][num];
            log_mem(format!("Get U ({},{}) L'xL'", u.row_size(), u.col_size()));
            let u_tensor_s = u.tensor(&s_i_num);
            let k_target_error = sampler_uniform.sample_uniform(
                &params,
                (1 + packed_input_size) * (d + 1),
                m_b,
                DistType::GaussDist { sigma: obf_params.p_sigma },
            );
            log_mem(format!(
                "Computed U ⊗ S ({},{})",
                u_tensor_s.row_size(),
                u_tensor_s.col_size()
            ));
            let k_target = (u_tensor_s * &b_star_level) + k_target_error;
            let k_preimage_num =
                sampler_trapdoor.preimage(&params, &b_star_trapdoor_cur, &b_star_cur, &k_target);
            log_mem(format!(
                "Computed k_preimage_num ({},{})",
                k_preimage_num.row_size(),
                k_preimage_num.col_size()
            ));
            handles_per_level.push(store_and_drop_matrix(
                k_preimage_num,
                &dir_path,
                &format!("k_preimage_{level}_{num}"),
            ));
            join_all(handles_per_level).await;
        }

        b_star_trapdoor_cur = b_star_trapdoor_level;
        b_star_cur = b_star_level;
    }

    /*
    =============================================================================
     Preimage‐generation for final BGG+ encoding evaluation step (after all bit insertion)

     1) Build the "final digits" circuit f[x_L] from a_decomposed, b_decomposed and the public circuit.
     2) Evaluate f[x_L] on (C = {a_decomposed, b_decomposed}, public key bits…).
     3) Form the target matrix Y = (eval_outputs_matrix + a_prf) ∥ 0.
     4) Sample a trapdoor preimage of Y under B*_star_basis → final_preimage.
    =============================================================================
    */

    #[cfg(feature = "bgm")]
    {
        player.play_music("bgm/obf_bgm5.mp3");
    }

    // Sample A_att
    let pub_key_att = sample_public_key_by_id(&bgg_pubkey_sampler, &params, 0, &reveal_plaintexts);
    let (first_key, other_keys) =
        pub_key_att.split_first().expect("pub_key_att must contain at least one key");
    let pub_key_att_matrix: M = first_key.concat_matrix(other_keys);
    log_mem(format!(
        "Sampled pub_key_att_matrix ({},{})",
        pub_key_att_matrix.row_size(),
        pub_key_att_matrix.col_size()
    ));
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_poly(hardcoded_key.clone(), &dir_path, "hardcoded_key"));
    let hardcoded_key_matrix = M::from_poly_vec_row(&params, vec![hardcoded_key]);

    // Generate RLWE ciphertext for the hardcoded key
    let sampler_uniform = SU::new();
    let a = public_data.a_rlwe_bar;
    let b = rlwe_encrypt(
        params.as_ref(),
        &sampler_uniform,
        &t_bar,
        &a,
        &hardcoded_key_matrix,
        obf_params.hardcoded_key_sigma,
    );
    log_mem("Generated RLWE ciphertext {a, b}");

    // Decompose RLWE ciphertext. This is required to lower down error growth for BGG+ encoding
    // evaluation.
    let a_decomposed = a.entry(0, 0).decompose_base(params.as_ref());
    let b_decomposed = b.entry(0, 0).decompose_base(params.as_ref());
    log_mem("Decomposed RLWE ciphertext into {BaseDecompose(a), BaseDecompose(b)}");
    handles.push(store_and_drop_matrix(b, &dir_path, "b"));

    // P_att := u_1_L' ⊗ A_att - I_L' ⊗ G_d+1
    // computing u_1_L' ⊗ A_att
    let params = &params;
    let zeros = M::zero(
        params.as_ref(),
        packed_input_size * pub_key_att_matrix.row_size(),
        pub_key_att_matrix.col_size(),
    );
    let k_lhs = pub_key_att_matrix.concat_rows(&[&zeros]);

    // computing I_L' ⊗ G_d+1
    let gadget = M::gadget_matrix(params, d + 1);
    let other_blocks: Vec<&M> = std::iter::repeat_n(&gadget, packed_input_size).collect();
    let k_rhs = gadget.concat_diag(&other_blocks);
    let final_preimage_target_att = k_lhs - k_rhs;
    let final_preimage_att = sampler_trapdoor.preimage(
        params,
        &b_star_trapdoor_cur,
        &b_star_cur,
        &final_preimage_target_att,
    );
    log_mem("Sampled final_preimage_att");
    handles.push(store_and_drop_matrix(final_preimage_att, &dir_path, "final_preimage_att"));
    handles.push(store_and_drop_matrix(pub_key_att_matrix, &dir_path, "pub_key_att"));

    // P_F := u_1_L' ⊗ (A_F + A_p)
    let final_preimage_target_f = {
        let final_circuit = build_final_digits_circuit::<M::P, BggPublicKey<M>>(
            &a_decomposed,
            &b_decomposed,
            public_circuit,
        );
        log_mem("Computed final_circuit");
        let eval_outputs = final_circuit.eval(params.as_ref(), &pub_key_att[0], &pub_key_att[1..]);
        log_mem("Evaluated outputs");
        debug_assert_eq!(eval_outputs.len(), log_base_q * packed_output_size);
        let output_ints = eval_outputs
            .par_chunks(log_base_q)
            .map(|digits| BggPublicKey::digits_to_int(digits, params))
            .collect::<Vec<_>>();
        let eval_outputs_matrix = output_ints[0].concat_matrix(&output_ints[1..]);
        #[cfg(feature = "debug")]
        assert_eq!(eval_outputs_matrix.col_size(), packed_output_size);
        #[cfg(feature = "debug")]
        handles.push(store_and_drop_matrix(
            eval_outputs_matrix.clone() + public_data.a_prf.clone(),
            &dir_path,
            "eval_outputs_matrix_plus_a_prf",
        ));
        let summat = eval_outputs_matrix + public_data.a_prf;
        let extra_blocks = packed_input_size;
        let zero_rows = extra_blocks * summat.row_size();
        let zero_cols = summat.col_size();
        let zeros = M::zero(params.as_ref(), zero_rows, zero_cols);
        summat.concat_rows(&[&zeros])
    };
    log_mem("Computed final_preimage_target_f");

    let final_preimage_f = sampler_trapdoor.preimage(
        params,
        &b_star_trapdoor_cur,
        &b_star_cur,
        &final_preimage_target_f,
    );
    log_mem("Sampled final_preimage_f");
    handles.push(store_and_drop_matrix(final_preimage_f, &dir_path, "final_preimage_f"));

    let store_hash_key = tokio::task::spawn_blocking(move || {
        let path = dir_path.join("hash_key");
        std::fs::write(&path, hash_key).expect("Failed to write hash_key file");
        log_mem("Stored hash_key");
    });
    handles.push(store_hash_key);

    join_all(handles).await;
}

fn store_and_drop_matrix<M: PolyMatrix + 'static>(
    matrix: M,
    dir_path: &Path,
    id: &str,
) -> tokio::task::JoinHandle<()> {
    let dir_path = dir_path.to_path_buf();
    let id_str = id.to_string();

    tokio::task::spawn_blocking(move || {
        log_mem(format!("Storing {id_str}"));
        Handle::current().block_on(async {
            matrix.write_to_files(&dir_path, &id_str).await;
        });
        drop(matrix);
        log_mem(format!("Stored {id_str}"));
    })
}

#[cfg(feature = "debug")]
fn store_and_drop_poly<P: Poly + 'static>(
    poly: P,
    dir_path: &Path,
    id: &str,
) -> tokio::task::JoinHandle<()> {
    let dir_path = dir_path.to_path_buf();
    let id_str = id.to_string();
    tokio::task::spawn_blocking(move || {
        log_mem(format!("Storing {id_str}"));
        Handle::current().block_on(async {
            poly.write_to_file(&dir_path, &id_str).await;
        });
        drop(poly);
        log_mem(format!("Stored {id_str}"));
    })
}
