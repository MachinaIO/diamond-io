#[cfg(feature = "bgm")]
use super::bgm::Player;

use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggEncoding, BggPublicKey, DigitsToInt,
    },
    io::{
        params::ObfuscationParams,
        utils::{build_final_digits_circuit, sample_public_key_by_id, PublicSampledData},
    },
    poly::{
        enc::rlwe_encrypt,
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyElem, PolyMatrix, PolyParams,
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
    let dim = params.ring_dimension() as usize;
    let d = obf_params.d;
    let sampler_uniform = SU::new();
    let sampler_trapdoor = ST::new(&params, obf_params.trapdoor_sigma);
    let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(hash_key, d);
    let packed_output_size = public_data.packed_output_size;

    /*
    =============================================================================
    We reveal all input slots so the evaluator can multiply encodings,
    but we keep the secret-key slot (last slot) hidden unless we're in debug mode.
    In the paper, the last slot corresponds to the FHE secret-key vector t.

    Sample the initial public key (level 0) with our reveal flags.

    Sample the initial BGG+ encodings.
    **NOTE:** the paper treats
        - c_att = encoding of (1, bits(X)), where X represents the evaluator's inputs.
        - c_t = encoding of the FHE secret-key
    as two distinct vectors.

    The length of encodings_init is (1 + packed_input_size + 1):
       - 1 for the encoding of 1
       - packed_input_size for the packed evaluator inputs,
       - 1 for the encoding of t (secret key).

    encodings_init is bundled into a single Vec<bool> of length L,
    where:
       - indices 0..L-2 → c_att
       - index   L-1   → c_t
    =============================================================================
    */

    // Sample FHE secret key t
    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    log_mem("Sampled t_bar_matrix");
    let minus_t_bar = -t_bar_matrix.entry(0, 0);
    let mut plaintexts =
        (0..(packed_input_size - 1)).map(|_| M::P::const_zero(&params)).collect_vec();
    plaintexts.push(minus_t_bar.clone());
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_poly(minus_t_bar, &dir_path, "minus_t_bar"));
    let mut reveal_plaintexts = vec![true; packed_input_size];
    reveal_plaintexts[packed_input_size - 1] = cfg!(feature = "debug");

    // Sample public key and secret key
    let pub_key_init = sample_public_key_by_id(&bgg_pubkey_sampler, &params, 0, &reveal_plaintexts);
    log_mem("Sampled pub key init");
    let s_bars = sampler_uniform.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
    log_mem("Sampled s_bars");
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bars,
        sampler_uniform,
        obf_params.encoding_sigma,
    );
    let s_init = bgg_encode_sampler.secret_vec;

    /*
    =============================================================================
    Pre‐loop initialization:

    1) Sample input‐dependent random basis B_* (trapdoor & public matrix).
    2) Compute initial secret key p_init = ((x_init , 1_L) ⊗ s_init)·B_* + error.
    where x_init (1, bits(FHE encryption), t), t is fhe secret key, 1_L is dummy mask for input space
    3) Create I_{d+1}, derive level_width, level_size, and depth.
    4) For each i in 0..level_size, build
            U_i = [ [ I_{d+1}, 0 ],
                    [ 0,       R_i ] ]
    5) Build wildcard matrix
            U_* = [ [ 0,       0     ],
                    [ I_{d+1}, I_{d+1} ] ]

    These values (B_*, p_init, I_{d+1}, U_i, U_*) are all set up
    before entering the main preimage generation loop.
    =============================================================================
    */

    // Sample input dependent random matrix B_*
    let (mut b_star_trapdoor_cur, mut b_star_cur) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
    log_mem("b star trapdoor init sampled");

    // Compute Initial secret key p_init := (((x_init , 1_L) ⊗ s_init)·B_*
    let hardcoded_key_matrix = M::from_poly_vec_row(&params, vec![hardcoded_key.clone()]);
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_poly(hardcoded_key, &dir_path, "hardcoded_key"));
    let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
    let mut coeffs = vec![zero_coeff; dim];
    let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
    let inserted_poly_matrix = M::from_poly_vec_row(&params, vec![inserted_poly]);
    let encoded_bits = hardcoded_key_matrix.concat_columns(&[&inserted_poly_matrix]);
    let s_connect = encoded_bits.tensor(&s_init);
    let s_b = s_connect * &b_star_cur;
    let p_init_error = bgg_encode_sampler.error_sampler.sample_uniform(
        &params,
        1,
        (2 * (d + 1)) * (2 + log_base_q),
        DistType::GaussDist { sigma: obf_params.p_sigma },
    );
    let p_init = s_b + p_init_error;
    log_mem("Computed p_init");
    handles.push(store_and_drop_matrix(p_init, &dir_path, "p_init"));
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_matrix(s_init, &dir_path, "s_init"));

    let identity_d_plus_1 = M::identity(params.as_ref(), d + 1, None);
    // number of bits to be inserted at each level
    let level_width = obf_params.level_width;
    assert_eq!(obf_params.input_size % level_width, 0);
    // otherwise we need >1 polynomial to insert the bits for each level
    assert!(level_width <= dim);
    // otherwise we get to a point in which the inserted bits have to be split between two
    // polynomials
    if obf_params.input_size > dim {
        assert_eq!(dim % level_width, 0);
    }
    let level_size = (1u64 << obf_params.level_width) as usize;
    // number of levels necessary to encode the input
    let depth = obf_params.input_size / level_width;
    let mut u_nums = Vec::with_capacity(level_size);
    for i in 0..level_size {
        let u_i = identity_d_plus_1.concat_diag(&[&public_data.rs[i]]);
        u_nums.push(u_i);
    }
    log_mem("Computed u_0, u_1");
    #[cfg(feature = "debug")]
    handles.push(store_and_drop_matrix(b_star_cur.clone(), &dir_path, "b_star_0"));

    /*
    Trapdoor preimage generation for the input insertion step.
    For each depth, sample K preimage at the corresponding level size.
    */

    for level in 0..depth {
        let (b_star_trapdoor_level, b_star_level) = sampler_trapdoor.trapdoor(&params, 2 * (d + 1));
        log_mem("Sampled b_star trapdoor for level");

        #[cfg(feature = "debug")]
        handles.push(store_and_drop_matrix(
            b_star_level.clone(),
            &dir_path,
            &format!("b_star_{}", level + 1),
        ));

        // Precomputation for k_preimage that are not num dependent
        let lhs = -pub_key_init[0].concat_matrix(&pub_key_init[1..]);
        let inserted_poly_index = 1 + (level * level_width) / dim;
        let inserted_coeff_indices =
            (0..level_width).map(|i| (i + (level * level_width)) % dim).collect_vec();
        debug_assert_eq!(inserted_coeff_indices.len(), level_width);

        for num in 0..level_size {
            let mut handles_per_level: Vec<tokio::task::JoinHandle<()>> = Vec::new();
            #[cfg(feature = "bgm")]
            {
                player.play_music(format!("bgm/obf_bgm{}.mp3", (2 * level + num) % 3 + 2));
            }
            let rg = &public_data.rgs[num];
            let top = lhs.mul_tensor_identity_decompose(rg, 1 + packed_input_size);
            log_mem("Computed top");
            // bit decompose num over level_width bits
            let num_bits: Vec<bool> = (0..level_width).map(|i| (num >> i) & 1 == 1).collect();
            debug_assert_eq!(num_bits.len(), level_width);
            // update input by updating coefficient
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
            let bottom = pub_key_init[0].concat_matrix(&pub_key_init[1..]) - &inserted_poly_gadget;
            log_mem("Computed bottom");
            let k_target = top.concat_rows(&[&bottom]);
            log_mem("Computed k_target");
            let k_preimage_num =
                sampler_trapdoor.preimage(&params, &b_star_trapdoor_cur, &b_star_level, &k_target);
            log_mem("Computed k_preimage_num");
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

     1) Build the “final digits” circuit f[x_L] from a_decomposed, b_decomposed and the public circuit.
     2) Evaluate f[x_L] on (C = {a_decomposed, b_decomposed}, public key bits…).
     3) Form the target matrix Y = (eval_outputs_matrix + a_prf) ∥ 0.
     4) Sample a trapdoor preimage of Y under B*_star_basis → final_preimage.
    =============================================================================
    */

    #[cfg(feature = "bgm")]
    {
        player.play_music("bgm/obf_bgm5.mp3");
    }

    // Generate RLWE ciphertext for the hardcoded key
    let sampler_uniform = SU::new();
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

    // Decompose RLWE ciphertext. This is required to lower down error growth for BGG+ encoding
    // evaluation.
    let a_decomposed = a.entry(0, 0).decompose_base(params.as_ref());
    let b_decomposed = b.entry(0, 0).decompose_base(params.as_ref());
    log_mem("Decomposed RLWE ciphertext into {BaseDecompose(a), BaseDecompose(b)}");
    handles.push(store_and_drop_matrix(b, &dir_path, "b"));

    let final_preimage_target = {
        let final_circuit = build_final_digits_circuit::<M::P, BggPublicKey<M>>(
            &a_decomposed,
            &b_decomposed,
            public_circuit,
        );
        log_mem("Computed final_circuit");
        let eval_outputs =
            final_circuit.eval(params.as_ref(), &pub_key_init[0], &pub_key_init[1..]);
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
    handles.push(store_and_drop_matrix(final_preimage, &dir_path, "final_preimage"));

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

fn store_and_drop_bgg_encoding<M: PolyMatrix + 'static>(
    encoding: BggEncoding<M>,
    dir_path: &Path,
    id: &str,
) -> tokio::task::JoinHandle<()> {
    let dir_path = dir_path.to_path_buf();
    let id_str = id.to_string();
    tokio::task::spawn_blocking(move || {
        log_mem(format!("Storing {id_str}"));
        Handle::current().block_on(async {
            encoding.write_to_files(&dir_path, &id_str).await;
        });
        drop(encoding);
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
