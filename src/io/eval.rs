#[cfg(feature = "bgm")]
use super::bgm::Player;
use super::{params::ObfuscationParams, utils::build_poly_vec};
#[cfg(feature = "debug")]
use crate::parallel_iter;
use crate::{
    bgg::{
        encoding::BggEncodingPltEvaluator, sampler::BGGPublicKeySampler, BggEncoding, DigitsToInt,
    },
    io::utils::{build_final_digits_circuit, sample_public_key_by_id, PublicSampledData},
    poly::{
        sampler::{PolyHashSampler, PolyTrapdoorSampler},
        Poly, PolyMatrix, PolyParams,
    },
    utils::{log_mem, timed_read},
};
use itertools::Itertools;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use std::{path::Path, sync::Arc, time::Duration};

pub fn evaluate<M, SH, ST, P>(
    obf_params: ObfuscationParams<M>,
    inputs: &[bool],
    dir_path: P,
) -> Vec<bool>
where
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M>,
    P: AsRef<Path>,
{
    #[cfg(feature = "bgm")]
    let player = Player::new();

    #[cfg(feature = "bgm")]
    {
        player.play_music("bgm/eval_bgm1.mp3");
    }
    let d = obf_params.d;
    let params = Arc::new(obf_params.params.clone());
    let log_base_q = params.modulus_digits();
    let dir_path = dir_path.as_ref().to_path_buf();
    assert_eq!(inputs.len(), obf_params.input_size);

    let mut total_load = Duration::ZERO;

    let hash_key = {
        let mut path = dir_path.clone();
        path.push("hash_key");
        let bytes = timed_read(
            "hash_key",
            || std::fs::read(&path).expect("Failed to read hash key file"),
            &mut total_load,
        );
        let mut hash_key = [0u8; 32];
        hash_key.copy_from_slice(&bytes);
        hash_key
    };
    log_mem("hash_key loaded");

    let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(hash_key, d);
    let public_data = PublicSampledData::<SH>::sample(&obf_params, hash_key);
    log_mem("Sampled public data");
    let packed_input_size = public_data.packed_input_size;
    let m_b = (1 + packed_input_size) * (d + 1) * (2 + log_base_q);
    let packed_output_size = public_data.packed_output_size;
    let mut p_cur = timed_read(
        "p_cur",
        || M::read_from_files(&obf_params.params, 1, m_b, &dir_path, "p_init"),
        &mut total_load,
    );
    log_mem(format!("p_init ({},{}) loaded", p_cur.row_size(), p_cur.col_size()));

    let reveal_plaintexts = {
        let mut reveals = vec![true; packed_input_size + 1];
        // Only reveal FHE secret key in debug mode.
        reveals[packed_input_size] = cfg!(feature = "debug");
        reveals
    };
    let params = Arc::new(obf_params.params.clone());
    let level_width = obf_params.level_width;
    assert!(inputs.len() % level_width == 0);
    let depth = obf_params.input_size / level_width;
    #[cfg(feature = "debug")]
    let s_init = timed_read(
        "s_init",
        || M::read_from_files(&obf_params.params, 1, d + 1, &dir_path, "s_init"),
        &mut total_load,
    );
    #[cfg(feature = "debug")]
    let minus_t_bar = timed_read(
        "minus_t_bar",
        || {
            <<M as PolyMatrix>::P as Poly>::read_from_file(
                &obf_params.params,
                &dir_path,
                "minus_t_bar",
            )
        },
        &mut total_load,
    );

    #[cfg(feature = "debug")]
    let (b_stars, durations) = parallel_iter!(0..(depth + 1))
        .map(|level| {
            let mut local_duration = Duration::ZERO;
            let b_star = timed_read(
                "b_star",
                || {
                    M::read_from_files(
                        params.as_ref(),
                        (1 + packed_input_size) * (d + 1),
                        m_b,
                        &dir_path,
                        &format!("b_star_{level}"),
                    )
                },
                &mut local_duration,
            );
            (b_star, local_duration)
        })
        .unzip::<_, _, Vec<_>, Vec<_>>();

    #[cfg(feature = "debug")]
    {
        for duration in durations {
            total_load += duration;
        }
    }

    #[cfg(feature = "debug")]
    if obf_params.hardcoded_key_sigma == 0.0 && obf_params.p_sigma == 0.0 {
        let one = M::P::const_one(&params);
        let mut plaintexts = vec![one];
        plaintexts
            .extend((0..packed_input_size - 1).map(|_| M::P::const_zero(&params)).collect_vec());
        plaintexts.push(minus_t_bar.clone());
        let encoded_bits = M::from_poly_vec_row(&params, plaintexts);
        let s_connect = encoded_bits.tensor(&s_init);
        let expected_p_init = s_connect * &b_stars[0];
        assert_eq!(p_cur, expected_p_init);
    }
    // Pack the input bits into little-endian limbs:
    // for every `level_width`-sized chunk, interpret the Booleans as a binary
    // number (bit 0 = LSB) and collect the resulting integers.
    let nums: Vec<u64> = inputs
        .chunks(level_width)
        .map(|chunk| {
            chunk.iter().enumerate().fold(0u64, |acc, (i, &bit)| acc + ((bit as u64) << i))
        })
        .collect();
    debug_assert_eq!(nums.len(), depth);
    #[cfg(feature = "debug")]
    let mut s_cur = s_init;

    for (level, num) in nums.iter().enumerate() {
        let level = level + 1;
        let k = timed_read(
            "k",
            || {
                M::read_from_files(
                    params.as_ref(),
                    m_b,
                    m_b,
                    &dir_path,
                    &format!("k_preimage_{level}_{num}"),
                )
            },
            &mut total_load,
        );
        log_mem(format!("k_{}_{} loaded ({},{})", level, num, k.row_size(), k.col_size()));
        let p = p_cur * k;
        log_mem(format!("p at {} computed ({},{})", level, p.row_size(), p.col_size()));

        #[cfg(feature = "debug")]
        if obf_params.hardcoded_key_sigma == 0.0 && obf_params.p_sigma == 0.0 {
            let s_matrix = timed_read(
                "s",
                || {
                    M::read_from_files(
                        params.as_ref(),
                        d + 1,
                        d + 1,
                        &dir_path,
                        &format!("s_{level}_{num}"),
                    )
                },
                &mut total_load,
            );
            s_cur = s_cur * s_matrix;
            let plaintexts = build_poly_vec::<M>(
                &params,
                inputs,
                level_width,
                level,
                obf_params.input_size,
                Some(minus_t_bar.clone()),
            );
            let encoded_bits = M::from_poly_vec_row(&params, plaintexts);
            let s_connect = encoded_bits.tensor(&s_cur);
            let expected_p = s_connect * &b_stars[level];
            assert_eq!(p, expected_p, "debug check failed at level {level}");
        }

        p_cur = p;
    }
    let m_b_small = (d + 1) * (2 + log_base_q);
    let k_plus_one = timed_read(
        "k_l_plus_one",
        || M::read_from_files(params.as_ref(), m_b, m_b_small, &dir_path, &format!("k_l_plus_one")),
        &mut total_load,
    );
    let p_l_plus_one = p_cur.clone() * k_plus_one;
    log_mem(format!(
        "p_l_plus_one computed ({},{})",
        p_l_plus_one.row_size(),
        p_l_plus_one.col_size()
    ));

    #[cfg(feature = "bgm")]
    {
        player.play_music("bgm/eval_bgm2.mp3");
    }

    let b = timed_read(
        "b",
        || M::read_from_files(&obf_params.params, 1, 1, &dir_path, "b"),
        &mut total_load,
    );
    log_mem("b loaded");

    let a_decomposed = public_data.a_rlwe_bar.entry(0, 0).decompose_base(&params);
    let b_decomposed = &b.entry(0, 0).decompose_base(&params);
    log_mem("a,b decomposed");
    let final_circuit = build_final_digits_circuit::<M::P, BggEncoding<M>>(
        &a_decomposed,
        b_decomposed,
        obf_params.public_circuit,
    );
    log_mem("final_circuit built");
    let final_preimage_f = timed_read(
        "final_preimage_f",
        || {
            M::read_from_files(
                &obf_params.params,
                m_b_small,
                packed_output_size,
                &dir_path,
                "final_preimage_f",
            )
        },
        &mut total_load,
    );
    log_mem("final_preimage loaded");
    // v := p_l_plus_one * K_F
    let final_v = p_l_plus_one.clone() * final_preimage_f;
    log_mem("final_v computed");

    #[cfg(feature = "debug")]
    if obf_params.hardcoded_key_sigma == 0.0 && obf_params.p_sigma == 0.0 {
        let eval_outputs_matrix_plus_a_prf = timed_read(
            "eval_outputs_matrix_plus_a_prf",
            || {
                M::read_from_files(
                    &obf_params.params,
                    d + 1,
                    packed_output_size,
                    &dir_path,
                    "eval_outputs_matrix_plus_a_prf",
                )
            },
            &mut total_load,
        );
        let expected_final_v = s_cur.clone() * eval_outputs_matrix_plus_a_prf;
        assert_eq!(final_v, expected_final_v);
        log_mem("final_v debug check passed");
    }
    let final_preimage_att = timed_read(
        "final_preimage_att",
        || {
            M::read_from_files(
                &obf_params.params,
                m_b,
                (1 + packed_input_size) * (d + 1) * log_base_q,
                &dir_path,
                "final_preimage_att",
            )
        },
        &mut total_load,
    );
    // c_att := p * K_att
    let c_att = p_cur.clone() * final_preimage_att;
    log_mem(format!("Computed c_att ({}, {})", c_att.row_size(), c_att.col_size()));
    let pub_key_att = sample_public_key_by_id(&bgg_pubkey_sampler, &params, 0, &reveal_plaintexts);
    log_mem(format!("Sampled pub_key_att {} ", pub_key_att.len()));

    let m = (d + 1) * log_base_q;
    #[cfg(not(feature = "debug"))]
    let polys =
        build_poly_vec::<M>(&params, inputs, level_width, nums.len(), obf_params.input_size, None);
    #[cfg(feature = "debug")]
    let polys = build_poly_vec::<M>(
        &params,
        inputs,
        level_width,
        nums.len(),
        obf_params.input_size,
        Some(minus_t_bar.clone()),
    );

    #[cfg(feature = "debug")]
    if obf_params.hardcoded_key_sigma == 0.0 && obf_params.p_sigma == 0.0 {
        let gadget = M::gadget_matrix(&params, d + 1);
        let plaintexts = M::from_poly_vec_row(&params, polys.clone());
        let pubkey = pub_key_att[0].concat_matrix(&pub_key_att[1..]);
        let inner = pubkey - plaintexts.tensor(&gadget);
        let expected_c_att = s_cur.clone() * inner;
        assert_eq!(c_att, expected_c_att);
        log_mem("c_att debug check passed");
    }

    let mut new_encodings = vec![];
    #[cfg(not(feature = "debug"))]
    let plaintexts_len = pub_key_att.len();
    for (j, pub_key) in pub_key_att.into_iter().enumerate() {
        let new_vec = c_att.slice_columns(j * m, (j + 1) * m);
        #[cfg(feature = "debug")]
        {
            let new_encode: BggEncoding<M> =
                BggEncoding::new(new_vec, pub_key, Some(polys[j].clone()));
            new_encodings.push(new_encode);
        }

        #[cfg(not(feature = "debug"))]
        {
            let new_encode: BggEncoding<M> = if j == plaintexts_len - 1 {
                BggEncoding::new(new_vec, pub_key, None)
            } else {
                BggEncoding::new(new_vec, pub_key, Some(polys[j].clone()))
            };
            new_encodings.push(new_encode);
        }
    }
    let bgg_encoding_plt_evaluator =
        BggEncodingPltEvaluator::<M, SH>::new(hash_key, dir_path.clone(), p_l_plus_one);
    let output_encodings = final_circuit.eval(
        params.as_ref(),
        &new_encodings[0],
        &new_encodings[1..],
        Some(bgg_encoding_plt_evaluator),
    );
    log_mem("final_circuit evaluated");
    let output_encoding_ints = output_encodings
        .par_chunks(log_base_q)
        .map(|digits| BggEncoding::digits_to_int(digits, &params))
        .collect::<Vec<_>>();
    let output_encodings_vec = output_encoding_ints[0].concat_vector(&output_encoding_ints[1..]);
    log_mem("final_circuit evaluated and recomposed");
    let z = output_encodings_vec - final_v;
    log_mem("z computed");
    debug_assert_eq!(z.size(), (1, packed_output_size));
    #[cfg(feature = "debug")]
    if obf_params.hardcoded_key_sigma == 0.0 && obf_params.p_sigma == 0.0 {
        let hardcoded_key = timed_read(
            "hardcoded_key",
            || {
                <<M as PolyMatrix>::P as Poly>::read_from_file(
                    &obf_params.params,
                    &dir_path,
                    "hardcoded_key",
                )
            },
            &mut total_load,
        );
        {
            let expected = s_cur *
                (output_encoding_ints[0].pubkey.matrix.clone() -
                    M::unit_column_vector(&params, d + 1, d) *
                        output_encoding_ints[0].plaintext.clone().unwrap());
            assert_eq!(output_encoding_ints[0].vector, expected);
        }
        if inputs[0] {
            assert_eq!(
                output_encoding_ints[0]
                    .plaintext
                    .clone()
                    .unwrap()
                    .extract_bits_with_threshold(&params),
                hardcoded_key.to_bool_vec()
            );
        }
        assert_eq!(z.size(), (1, packed_output_size));
    }
    log_mem(format!("total loading time {total_load:?}"));
    z.get_row(0).into_iter().flat_map(|p| p.extract_bits_with_threshold(&params)).collect_vec()
}
