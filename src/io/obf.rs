use super::{utils::*, Obfuscation, ObfuscationParams};
use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggPublicKey,
    },
    join,
    poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams},
    utils::log_mem,
};
use rand::{Rng, RngCore};
use std::{ops::Mul, path::Path, sync::Arc};
use tracing::info;

pub fn obfuscate<M, SU, SH, ST, R>(
    obf_params: ObfuscationParams<M>,
    sampler_uniform: SU,
    mut sampler_hash: SH,
    sampler_trapdoor: ST,
    rng: &mut R,
    fs_dir_path: &Path,
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
    info!("Obfuscating circuit");
    log_mem();
    let public_circuit = &obf_params.public_circuit;
    let dim = obf_params.params.ring_dimension() as usize;
    let log_q = obf_params.params.modulus_bits();
    debug_assert_eq!(public_circuit.num_input(), log_q + obf_params.input_size);
    let hash_key = rng.random::<[u8; 32]>();
    info!("Gen hash key");
    log_mem();
    sampler_hash.set_key(hash_key);
    info!("Set hash key");
    log_mem();
    let sampler_uniform = Arc::new(sampler_uniform);
    info!("Created sampler_uniform");
    log_mem();
    let sampler_trapdoor = Arc::new(sampler_trapdoor);
    info!("Created sampler_trapdoor");
    log_mem();
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(Arc::new(sampler_hash));
    info!("Created bgg_pubkey_sampler");
    log_mem();
    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    info!("Sampled public data");
    log_mem();
    let params = Arc::new(obf_params.params);
    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let s_bar =
        sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist).entry(0, 0).clone();
    info!("Sampled s_bar");
    log_mem();
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bar,
        sampler_uniform.clone(),
        obf_params.error_gauss_sigma,
    );
    info!("Created bgg_encode_sampler");
    log_mem();
    let s_init = &bgg_encode_sampler.secret_vec;
    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    info!("Sampled t_bar_matrix");
    log_mem();
    let hardcoded_key_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
    info!("Sampled hardcoded_key_matrix");
    log_mem();
    let enc_hardcoded_key = {
        let e = sampler_uniform.sample_uniform(
            &params,
            1,
            1,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        info!("Sampled enc_hardcoded_key");
        log_mem();
        let scale = M::P::from_const(&params, &<M::P as Poly>::Elem::half_q(&params.modulus()));
        &t_bar_matrix * &public_data.a_rlwe_bar + &e - &(&hardcoded_key_matrix * &scale)
    };
    info!("Computed enc_hardcoded_key");
    log_mem();

    let enc_hardcoded_key_polys = enc_hardcoded_key.get_column_matrix(0).decompose().get_column(0);
    info!("Decomposed enc_hardcoded_key");
    log_mem();
    let t_bar = t_bar_matrix.entry(0, 0).clone();
    #[cfg(test)]
    let hardcoded_key = hardcoded_key_matrix.entry(0, 0).clone();

    let mut plaintexts: Vec<M::P> = (0..obf_params.input_size.div_ceil(dim))
        .map(|_| M::P::const_zero(params.as_ref()))
        .collect();
    plaintexts.push(t_bar.clone());

    info!("Generated plaintexts");
    log_mem();

    // let mut input_encoded_polys: Vec<<M as PolyMatrix>::P> =
    //     vec![<M::P as Poly>::const_one(params.as_ref())];
    // input_encoded_polys.extend(enc_hardcoded_key_polys);
    // input_encoded_polys.extend(zero_plaintexts);
    let encodings_init = bgg_encode_sampler.sample(&params, &public_data.pubkeys[0], &plaintexts);
    // let encode_fhe_key =
    //     bgg_encode_sampler.sample(&params, &public_data.pubkeys_fhe_key[0], &t.get_row(0),
    // false);

    info!("Sampled encodings_init");
    log_mem();

    let mut bs_path = Vec::with_capacity(obf_params.input_size);
    let mut b_trapdoors = Vec::with_capacity(obf_params.input_size);
    for _i in 0..=obf_params.input_size {
        let (b_0_trapdoor, b_0) = sampler_trapdoor.trapdoor(&params, 4);
        info!("Sampled b_0 trapdoor for input size {}", _i);
        log_mem();
        let b_0_matrix_path = fs_dir_path.join(format!("b_0_{}", _i));
        b_0.store(&b_0_matrix_path);
        drop(b_0);
        info!("Drop b_0 for input size {}", _i);
        log_mem();
        let (b_1_trapdoor, b_1) = sampler_trapdoor.trapdoor(&params, 4);
        info!("Sampled b_1 trapdoor for input size {}", _i);
        log_mem();
        let b_1_matrix_path = fs_dir_path.join(format!("b_1_{}", _i));
        b_1.store(&b_1_matrix_path);
        drop(b_1);
        info!("Dropped b_1 for input size {}", _i);
        log_mem();
        let (b_star_trapdoor, b_star) = sampler_trapdoor.trapdoor(&params, 4);
        info!("Sampled b_star trapdoor for input size {}", _i);
        log_mem();
        let b_star_matrix_path = fs_dir_path.join(format!("b_star_{}", _i));
        b_star.store(&b_star_matrix_path);
        drop(b_star);
        info!("Dropped b_star for input size {}", _i);
        log_mem();
        bs_path.push((b_0_matrix_path, b_1_matrix_path, b_star_matrix_path));
        b_trapdoors.push((b_0_trapdoor, b_1_trapdoor, b_star_trapdoor));
    }
    let m_b = 4 * (2 + log_q);

    let identity_2 = M::identity(params.as_ref(), 2, None);
    info!("Generate identity_2");
    log_mem();
    let u_0 = identity_2.concat_diag(&[&public_data.r_0]);
    let u_1 = identity_2.concat_diag(&[&public_data.r_1]);
    let u_star = {
        let zeros = M::zero(params.as_ref(), 2, 4);
        let identities = identity_2.concat_columns(&[&identity_2]);
        zeros.concat_rows(&[&identities])
    };
    info!("Generated u_0, u_1, u_star by concatenations");
    log_mem();
    let gadget_2 = M::gadget_matrix(params.as_ref(), 2);
    let (mut m_preimages_paths, mut n_preimages_paths, mut k_preimages_paths) = (
        Vec::with_capacity(obf_params.input_size),
        Vec::with_capacity(obf_params.input_size),
        Vec::with_capacity(obf_params.input_size),
    );
    info!("Starting preimage calculations");
    log_mem();
    for idx in 0..obf_params.input_size {
        info!("Processing preimage for input {}/{}", idx + 1, obf_params.input_size);
        log_mem();
        let (_, _, b_cur_star_path) = &bs_path[idx];
        let (b_next_0_path, b_next_1_path, b_next_star_path) = &bs_path[idx + 1];
        let (_, _, b_cur_star_trapdoor) = &b_trapdoors[idx];
        let (b_next_0_trapdoor, b_next_1_trapdoor, _) = &b_trapdoors[idx + 1];
        let m_preimage = |a, m_i| {
            info!("Computed m_preimage for input {} bit {}", idx + 1, m_i);
            let b_cur_star = M::load(b_cur_star_path);
            let m: M =
                sampler_trapdoor.preimage(params.as_ref(), b_cur_star_trapdoor, &b_cur_star, &a);
            drop(b_cur_star);
            drop(a);
            log_mem();
            let m_path = fs_dir_path.join(format!("m_preimage_{}_{}", m_i, idx));
            m.store(&m_path);
            drop(m);
            info!("Dropped m_preimage for input {} bit {}", idx + 1, m_i);
            log_mem();
            m_path
        };
        info!("Computed ub_star for input {}", idx + 1);
        log_mem();
        let b_next_0 = M::load(b_next_0_path);
        let b_next_1 = M::load(b_next_1_path);
        let mp = || join!(|| m_preimage(&u_0 * &b_next_0, 0), || m_preimage(&u_1 * &b_next_1, 1));
        let b_next_star = M::load(b_next_star_path);
        let ub_star = &u_star * &b_next_star;
        info!("Computed ub_star for input {}", idx + 1);
        log_mem();

        let n_preimage = |t, n, n_idx| {
            info!("Before computing n_preimage for input {} bit {}", idx + 1, n_idx);
            log_mem();
            let matrix_n = sampler_trapdoor.preimage(&params, t, n, &ub_star);
            info!("Computed n_preimage for input {} bit {}", idx + 1, n_idx);
            log_mem();
            let n_path = fs_dir_path.join(format!("n_preimage_{}_{}", n_idx, idx));
            matrix_n.store(&n_path);
            drop(matrix_n);
            info!("Dropped n_preimage for input {} bit {}", idx + 1, n_idx);
            log_mem();
            n_path
        };

        let np = || {
            join!(|| n_preimage(b_next_0_trapdoor, &b_next_0, 0), || n_preimage(
                b_next_1_trapdoor,
                &b_next_1,
                1
            ))
        };
        let k_preimage = |bit: usize| {
            info!("Before computing k_preimage for input {} bit {}", idx + 1, bit);
            log_mem();

            let rg = &public_data.rgs[bit];
            let lhs = -public_data.pubkeys[idx][0].concat_matrix(&public_data.pubkeys[idx][1..]);

            info!("Performed lhs computation for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();

            let top = lhs.mul_tensor_identity_decompose(rg, 1 + packed_input_size);

            info!("Computed top matrix for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();

            let inserted_poly_index = 1 + idx / dim;
            let inserted_coeff_index = idx % dim;
            let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            let mut coeffs = vec![zero_coeff; dim];
            if bit != 0 {
                coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus())
            };
            let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);

            info!("Created inserted_poly for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();

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
                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_2)
            };

            info!("Created inserted_poly_gadget for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();

            let bottom = public_data.pubkeys[idx + 1][0]
                .concat_matrix(&public_data.pubkeys[idx + 1][1..]) -
                &inserted_poly_gadget;
            let k_target = top.concat_rows(&[&bottom]);

            info!("Computed k_target for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();

            let b_matrix = if bit == 0 { b_next_0.clone() } else { b_next_1.clone() };
            let trapdoor = if bit == 0 { b_next_0_trapdoor } else { b_next_1_trapdoor };
            let k_matrix = sampler_trapdoor.preimage(&params, trapdoor, &b_matrix, &k_target);
            info!("Computed k_matrix for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();
            let k_path = fs_dir_path.join(format!("k_preimage_{}_{}", bit, idx));
            k_matrix.store(&k_path);
            drop(k_matrix);
            info!("Dropped k_matrix for k_preimage input {} bit {}", idx + 1, bit);
            log_mem();
            k_path
        };

        let kp = || join!(|| k_preimage(0), || k_preimage(1));

        let (mp, (np, kp)) = join!(mp, || join!(np, kp));

        n_preimages_paths.push(np);
        m_preimages_paths.push(mp);
        k_preimages_paths.push(kp);
    }

    info!("Finished all preimage calculations");
    log_mem();

    let a_decomposed_polys = public_data.a_rlwe_bar.get_column_matrix(0).decompose().get_column(0);

    info!("Decomposed a_rlwe_bar");
    log_mem();

    let final_circuit = build_final_step_circuit::<_, BggPublicKey<M>>(
        &params,
        &a_decomposed_polys,
        &enc_hardcoded_key_polys,
        public_circuit.clone(),
    );

    info!("Built final step circuit");
    log_mem();

    let final_preimage_target = {
        let one = public_data.pubkeys[obf_params.input_size][0].clone();
        let input = &public_data.pubkeys[obf_params.input_size][1..];
        let eval_outputs = final_circuit.eval(params.as_ref(), one, input);
        info!("Evaluated final circuit");
        log_mem();
        let mut eval_outputs_matrix = eval_outputs[0].concat_matrix(&eval_outputs[1..]);
        let unit_vector = identity_2.slice_columns(1, 2);
        eval_outputs_matrix = eval_outputs_matrix * unit_vector.decompose();
        info!("Computed eval_outputs_matrix");
        log_mem();
        debug_assert_eq!(eval_outputs_matrix.col_size(), packed_output_size);
        (eval_outputs_matrix + public_data.a_prf).concat_rows(&[&M::zero(
            params.as_ref(),
            2,
            packed_output_size,
        )])
    };
    let (_, _, b_final_path) = &bs_path[obf_params.input_size];
    let (_, _, b_final_trapdoor) = &b_trapdoors[obf_params.input_size];
    let b_final = M::load(b_final_path);
    let final_preimage =
        sampler_trapdoor.preimage(&params, b_final_trapdoor, &b_final, &final_preimage_target);
    let final_preimage_path = fs_dir_path.join("final_preimage");
    final_preimage.store(&final_preimage_path);
    drop(final_preimage);
    info!("Computed final preimage");
    log_mem();

    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init]);
        let bs_0_2 = M::load(&bs_path[0].2);
        let s_b = s_connect * &bs_0_2;
        let error = sampler_uniform.sample_uniform(
            &params,
            1,
            m_b,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        s_b + error
    };
    let p_init_path = fs_dir_path.join("p_init");
    p_init.store(&p_init_path);
    drop(p_init);
    info!("Computed p_init");
    log_mem();

    info!("OBFUSCATION COMPLETED");
    Obfuscation {
        hash_key,
        enc_hardcoded_key,
        encodings_init,
        p_init_path,
        m_preimages_paths,
        n_preimages_paths,
        k_preimages_paths,
        final_preimage_path,
        #[cfg(test)]
        s_init: s_init.clone(),
        #[cfg(test)]
        t_bar: t_bar.clone(),
        #[cfg(test)]
        bs_path,
        #[cfg(test)]
        hardcoded_key: hardcoded_key.clone(),
        #[cfg(test)]
        final_preimage_target,
    }
}
