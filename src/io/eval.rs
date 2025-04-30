#[cfg(feature = "bgm")]
use super::bgm::Player;
use super::params::ObfuscationParams;
use crate::{
    bgg::{sampler::BGGPublicKeySampler, BggEncoding, DigitsToInt},
    io::utils::{build_final_digits_circuit, sample_public_key_by_id, PublicSampledData},
    parallel_iter,
    poly::{
        element::PolyElem,
        sampler::{PolyHashSampler, PolyTrapdoorSampler},
        Poly, PolyMatrix, PolyParams,
    },
    utils::log_mem,
};
use itertools::Itertools;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use std::{path::Path, sync::Arc};

pub fn evaluate<M, SH, ST, P>(
    obf_params: ObfuscationParams<M>,
    inputs: &[bool],
    dir_path: P,
) -> Vec<bool>
where
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M>,
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
    let d_plus_1 = d + 1;
    let params = Arc::new(obf_params.params.clone());
    let log_base_q = params.modulus_digits();
    let dim = params.ring_dimension() as usize;
    let m_b = (2 * d_plus_1) * (2 + log_base_q);
    let dir_path = dir_path.as_ref().to_path_buf();
    assert_eq!(inputs.len(), obf_params.input_size);

    let hash_key = {
        let mut path = dir_path.clone();
        path.push("hash_key");
        let bytes = std::fs::read(&path).expect("Failed to read hash key file");
        let mut hash_key = [0u8; 32];
        hash_key.copy_from_slice(&bytes);
        hash_key
    };
    log_mem("hash_key loaded");

    let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(hash_key, d);
    let public_data = PublicSampledData::<SH>::sample(&obf_params, hash_key);
    log_mem("Sampled public data");

    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let mut p_cur = M::read_from_files(&obf_params.params, 1, m_b, &dir_path, "p_init");
    log_mem("p_init loaded");

    #[cfg(feature = "debug")]
    let reveal_plaintexts = [vec![true; packed_input_size], vec![true; 1]].concat();
    #[cfg(not(feature = "debug"))]
    let reveal_plaintexts = [vec![true; packed_input_size], vec![false; 1]].concat();
    let params = Arc::new(obf_params.params.clone());
    let level_width = obf_params.level_width;
    #[cfg(feature = "debug")]
    let level_size = (1u64 << obf_params.level_width) as usize;
    assert!(inputs.len() % level_width == 0);
    let depth = obf_params.input_size / level_width;

    // drop the first element from `reveal_plaintexts`` (this is set to true for `one` encofing
    // within the sample logic)
    let reveal_plaintexts = reveal_plaintexts[1..].to_vec();
    let pub_key_init = sample_public_key_by_id(&bgg_pubkey_sampler, &params, 0, &reveal_plaintexts);
    log_mem("Sampled pub_key_init");
    #[cfg(feature = "debug")]
    let s_init = M::read_from_files(&obf_params.params, 1, d_plus_1, &dir_path, "s_init");
    #[cfg(feature = "debug")]
    let minus_t_bar = <<M as PolyMatrix>::P as Poly>::read_from_file(
        &obf_params.params,
        &dir_path,
        "minus_t_bar",
    );

    #[cfg(feature = "debug")]
    let b_stars = parallel_iter!(0..depth + 1)
        .map(|level| {
            let b_star = M::read_from_files(
                params.as_ref(),
                2 * d_plus_1,
                m_b,
                &dir_path,
                &format!("b_star_{level}"),
            );
            b_star
        })
        .collect::<Vec<_>>();

    #[cfg(feature = "debug")]
    if obf_params.encoding_sigma == 0.0 &&
        obf_params.hardcoded_key_sigma == 0.0 &&
        obf_params.p_sigma == 0.0
    {
        let hardcoded_key = <<M as PolyMatrix>::P as Poly>::read_from_file(
            &obf_params.params,
            &dir_path,
            "hardcoded_key",
        );
        let hardcoded_key_matrix = M::from_poly_vec_row(&params, vec![hardcoded_key]);
        let zero = <M::P as Poly>::const_zero(&params);
        let one = <M::P as Poly>::const_one(&params);
        let mut polys = vec![];
        polys.push(one);
        for _ in 0..(packed_input_size - 1) {
            polys.push(zero.clone());
        }
        polys.push(minus_t_bar.clone());
        let inserted_poly_matrix = M::from_poly_vec_row(&params, polys);
        let encoded_bits = hardcoded_key_matrix.concat_columns(&[&inserted_poly_matrix]);
        let s_connect = encoded_bits.tensor(&s_init);
        let expected_p_init = s_connect * &b_stars[0];
        assert_eq!(p_cur, expected_p_init);
    }
    let nums: Vec<u64> = inputs
        .chunks(level_width)
        .map(|chunk| {
            chunk.iter().enumerate().fold(0u64, |acc, (i, &bit)| acc + ((bit as u64) << i))
        })
        .collect();
    debug_assert_eq!(nums.len(), depth);

    for (level, num) in nums.iter().enumerate() {
        let k_columns = (1 + packed_input_size) * d_plus_1 * log_base_q;
        let k = M::read_from_files(
            params.as_ref(),
            m_b,
            k_columns,
            &dir_path,
            &format!("k_preimage_{level}_{num}"),
        );
        log_mem(format!("k at {} loaded", level));
        let p = p_cur * k;
        log_mem(format!("p at {} computed", level));
        let inserted_poly_index = 1 + (level * level_width) / dim;
        p_cur = p.clone();
        #[cfg(feature = "debug")]
        if obf_params.encoding_sigma == 0.0 &&
            obf_params.hardcoded_key_sigma == 0.0 &&
            obf_params.p_sigma == 0.0
        {
            let dim = params.ring_dimension() as usize;
            let one = <M::P as Poly>::const_one(&params);
            let mut polys = vec![];
            polys.push(one);
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
            let input_polys =
                coeffs.chunks(dim).map(|coeffs| M::P::from_coeffs(&params, coeffs)).collect_vec();
            polys.extend(input_polys);
            polys.push(minus_t_bar.clone());
            let inserted_poly_matrix = M::from_poly_vec_row(&params, polys);
            let hardcoded_key = <<M as PolyMatrix>::P as Poly>::read_from_file(
                &obf_params.params,
                &dir_path,
                "hardcoded_key",
            );
            let hardcoded_key_matrix = M::from_poly_vec_row(&params, vec![hardcoded_key]);
            let encoded_bits = hardcoded_key_matrix.concat_columns(&[&inserted_poly_matrix]);
            let s_connect = encoded_bits.tensor(&s_init);
            let expected_p = s_connect * &b_stars[level + 1];
            assert_eq!(p, expected_p);
        }
    }

    #[cfg(feature = "bgm")]
    {
        player.play_music("bgm/eval_bgm2.mp3");
    }

    let b = M::read_from_files(&obf_params.params, 1, 1, &dir_path, "b");
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
    // todo: build last_input_encodings
    let last_input_encodings = encodings_cur;
    let output_encodings = final_circuit.eval::<BggEncoding<M>>(
        &params,
        &last_input_encodings[0],
        &last_input_encodings[1..],
    );
    log_mem("final_circuit evaluated");
    let output_encoding_ints = output_encodings
        .par_chunks(log_base_q)
        .map(|digits| BggEncoding::digits_to_int(digits, &params))
        .collect::<Vec<_>>();
    let output_encodings_vec = output_encoding_ints[0].concat_vector(&output_encoding_ints[1..]);
    log_mem("final_circuit evaluated and recomposed");
    let final_preimage = M::read_from_files(
        &obf_params.params,
        m_b,
        packed_output_size,
        &dir_path,
        "final_preimage",
    );
    log_mem("final_preimage loaded");
    let final_v = p_cur * final_preimage;
    log_mem("final_v computed");
    let z = output_encodings_vec - final_v;
    log_mem("z computed");
    debug_assert_eq!(z.size(), (1, packed_output_size));
    #[cfg(feature = "debug")]
    if obf_params.encoding_sigma == 0.0 &&
        obf_params.hardcoded_key_sigma == 0.0 &&
        obf_params.p_sigma == 0.0
    {
        let hardcoded_key = <<M as PolyMatrix>::P as Poly>::read_from_file(
            &obf_params.params,
            &dir_path,
            "hardcoded_key",
        );

        let mut last_s = s_init.clone();
        for num in nums.iter() {
            let r = public_data.rs[*num as usize].clone();
            last_s = last_s * r;
        }
        {
            let expected = last_s *
                (output_encoding_ints[0].pubkey.matrix.clone() -
                    M::unit_column_vector(&params, d_plus_1, d_plus_1 - 1) *
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
                hardcoded_key.to_bool_vec()
            );
        }
    }
    z.get_row(0).into_iter().flat_map(|p| p.extract_bits_with_threshold(&params)).collect_vec()
}
