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
    sampler_hash.set_key(hash_key);
    log_mem();

    let sampler_uniform = Arc::new(sampler_uniform);
    let sampler_trapdoor = Arc::new(sampler_trapdoor);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(Arc::new(sampler_hash));

    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    info!("Sampled public data");
    log_mem();

    let params = Arc::new(obf_params.params);
    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;

    let s_bar =
        sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist).entry(0, 0).clone();
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bar,
        sampler_uniform.clone(),
        obf_params.error_gauss_sigma,
    );
    let s_init = &bgg_encode_sampler.secret_vec;

    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    let hardcoded_key_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

    let t_bar = t_bar_matrix.entry(0, 0).clone();

    #[cfg(test)]
    let hardcoded_key = hardcoded_key_matrix.entry(0, 0).clone();

    let mut plaintexts: Vec<M::P> = Vec::with_capacity(obf_params.input_size.div_ceil(dim) + 1);
    for _ in 0..obf_params.input_size.div_ceil(dim) {
        plaintexts.push(M::P::const_zero(params.as_ref()));
    }
    plaintexts.push(t_bar.clone());

    let mut bs_path = Vec::with_capacity(obf_params.input_size + 1);
    let mut b_trapdoors = Vec::with_capacity(obf_params.input_size + 1);

    // Generate trapdoors with immediate file storage to reduce memory usage
    for i in 0..=obf_params.input_size {
        info!("Generating trapdoors for input {}/{}", i, obf_params.input_size);

        // Generate b_0 and its trapdoor
        let (b_0, b_0_trapdoor_first, b_0_trapdoor_second) = sampler_trapdoor.trapdoor(&params, 4);
        let b_0_trapdoor_first_path = fs_dir_path.join(format!("b_0_trapdoor_first_{}", i));
        b_0_trapdoor_first.store(&b_0_trapdoor_first_path);
        drop(b_0_trapdoor_first);
        let b_0_trapdoor_second_path = fs_dir_path.join(format!("b_0_trapdoor_second_{}", i));
        b_0_trapdoor_second.store(&b_0_trapdoor_second_path);
        drop(b_0_trapdoor_second);
        let b_0_matrix_path = fs_dir_path.join(format!("b_0_{}", i));
        b_0.store(&b_0_matrix_path);
        drop(b_0);

        // Generate b_1 and its trapdoor
        let (b_1, b_1_trapdoor_first, b_1_trapdoor_second) = sampler_trapdoor.trapdoor(&params, 4);
        let b_1_trapdoor_first_path = fs_dir_path.join(format!("b_1_trapdoor_first_{}", i));
        b_1_trapdoor_first.store(&b_1_trapdoor_first_path);
        drop(b_1_trapdoor_first);
        let b_1_trapdoor_second_path = fs_dir_path.join(format!("b_1_trapdoor_second_{}", i));
        b_1_trapdoor_second.store(&b_1_trapdoor_second_path);
        drop(b_1_trapdoor_second);
        let b_1_matrix_path = fs_dir_path.join(format!("b_1_{}", i));
        b_1.store(&b_1_matrix_path);
        drop(b_1);

        // Generate b_star and its trapdoor
        let (b_star, b_star_trapdoor_first, b_star_trapdoor_second) =
            sampler_trapdoor.trapdoor(&params, 4);
        let b_star_trapdoor_first_path = fs_dir_path.join(format!("b_star_trapdoor_first_{}", i));
        b_star_trapdoor_first.store(&b_star_trapdoor_first_path);
        drop(b_star_trapdoor_first);
        let b_star_trapdoor_second_path = fs_dir_path.join(format!("b_star_trapdoor_second_{}", i));
        b_star_trapdoor_second.store(&b_star_trapdoor_second_path);
        drop(b_star_trapdoor_second);
        let b_star_matrix_path = fs_dir_path.join(format!("b_star_{}", i));
        b_star.store(&b_star_matrix_path);
        drop(b_star);

        // Store paths
        bs_path.push((b_0_matrix_path, b_1_matrix_path, b_star_matrix_path));
        b_trapdoors.push((
            b_0_trapdoor_first_path,
            b_0_trapdoor_second_path,
            b_1_trapdoor_first_path,
            b_1_trapdoor_second_path,
            b_star_trapdoor_first_path,
            b_star_trapdoor_second_path,
        ));

        log_mem();
    }

    let m_b = 4 * (2 + log_q);

    let identity_2 = M::identity(params.as_ref(), 2, None);
    let u_0 = identity_2.concat_diag(&[&public_data.r_0]);
    let u_1 = identity_2.concat_diag(&[&public_data.r_1]);
    let u_star = {
        let zeros = M::zero(params.as_ref(), 2, 4);
        let identities = identity_2.concat_columns(&[&identity_2]);
        zeros.concat_rows(&[&identities])
    };

    let gadget_2 = M::gadget_matrix(params.as_ref(), 2);

    // Prepare storage for preimage paths
    let mut m_preimages_paths = Vec::with_capacity(obf_params.input_size);
    let mut n_preimages_paths = Vec::with_capacity(obf_params.input_size);
    let mut k_preimages_paths = Vec::with_capacity(obf_params.input_size);

    // Process preimages with careful memory management
    for idx in 0..obf_params.input_size {
        info!("Processing preimage for input {}/{}", idx + 1, obf_params.input_size);
        log_mem();

        // Get paths for current and next b matrices
        let (_, _, b_cur_star_path) = &bs_path[idx];
        let (b_next_0_path, b_next_1_path, b_next_star_path) = &bs_path[idx + 1];
        let (_, _, _, _, b_cur_star_trapdoor_first_path, b_cur_star_trapdoor_second_path) =
            &b_trapdoors[idx];
        let (
            b_next_0_trapdoor_first,
            b_next_0_trapdoor_second,
            b_next_1_trapdoor_first,
            b_next_1_trapdoor_second,
            _,
            _,
        ) = &b_trapdoors[idx + 1];

        // Create m_preimage function that handles memory cleanup
        let m_preimage = |a, m_i| {
            info!("Computing m_preimage for input {} bit {}", idx + 1, m_i);
            log_mem();

            let b_cur_star = M::load(b_cur_star_path);
            let b_cur_star_trapdoor_first = M::load(b_cur_star_trapdoor_first_path);
            let b_cur_star_trapdoor_second = M::load(b_cur_star_trapdoor_second_path);

            info!("Loaded m_preimage for input {} bit {}", idx + 1, m_i);
            log_mem();

            let m: M = sampler_trapdoor.preimage(
                params.as_ref(),
                &b_cur_star,
                &b_cur_star_trapdoor_first,
                &b_cur_star_trapdoor_second,
                &a,
            );

            info!("Sampled m_preimage for input {} bit {}", idx + 1, m_i);
            log_mem();

            drop(b_cur_star);
            drop(b_cur_star_trapdoor_first);
            drop(b_cur_star_trapdoor_second);
            drop(a);

            // Store result
            let m_path = fs_dir_path.join(format!("m_preimage_{}_{}", m_i, idx));
            m.store(&m_path);
            drop(m);

            m_path
        };

        // Load matrices one at a time to minimize memory usage
        let b_next_0 = M::load(b_next_0_path);
        let b_next_1 = M::load(b_next_1_path);

        // Compute m preimages in parallel
        let mp = || join!(|| m_preimage(&u_0 * &b_next_0, 0), || m_preimage(&u_1 * &b_next_1, 1));

        let b_next_star = M::load(b_next_star_path);
        let ub_star = &u_star * &b_next_star;

        // Define n_preimage function
        let n_preimage = |n, t_first_path, t_second_path, n_idx| {
            info!("Computing n_preimage for input {} bit {}", idx + 1, n_idx);
            log_mem();

            let t_first = M::load(t_first_path);
            let t_second = M::load(t_second_path);
            info!("Loaded n_preimage for input {} bit {}", idx + 1, n_idx);
            log_mem();

            let matrix_n = sampler_trapdoor.preimage(&params, n, &t_first, &t_second, &ub_star);
            info!("Sampled n_preimage for input {} bit {}", idx + 1, n_idx);
            log_mem();

            drop(t_first);
            drop(t_second);

            let n_path = fs_dir_path.join(format!("n_preimage_{}_{}", n_idx, idx));
            matrix_n.store(&n_path);
            drop(matrix_n);

            n_path
        };

        let np = || {
            join!(
                || n_preimage(&b_next_0, b_next_0_trapdoor_first, b_next_0_trapdoor_second, 0),
                || n_preimage(&b_next_1, b_next_1_trapdoor_first, b_next_1_trapdoor_second, 1)
            )
        };

        let k_preimage = |bit: usize| {
            info!("Computing k_preimage for input {} bit {}", idx + 1, bit);
            log_mem();

            let rg = &public_data.rgs[bit];
            let lhs = -public_data.pubkeys[idx][0].concat_matrix(&public_data.pubkeys[idx][1..]);

            let top = lhs.mul_tensor_identity_decompose(rg, 1 + packed_input_size);
            drop(lhs);

            let inserted_poly_index = 1 + idx / dim;
            let inserted_coeff_index = idx % dim;
            let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            let mut coeffs = vec![zero_coeff; dim];
            if bit != 0 {
                coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus())
            };
            let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);

            let inserted_poly_gadget = {
                let zero = <M::P as Poly>::const_zero(params.as_ref());
                let mut polys = Vec::with_capacity(packed_input_size + 1);

                for _ in 0..(inserted_poly_index) {
                    polys.push(zero.clone());
                }
                polys.push(inserted_poly);
                for _ in (inserted_poly_index + 1)..(packed_input_size + 1) {
                    polys.push(zero.clone());
                }

                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_2)
            };

            let bottom = public_data.pubkeys[idx + 1][0]
                .concat_matrix(&public_data.pubkeys[idx + 1][1..]) -
                &inserted_poly_gadget;

            drop(inserted_poly_gadget);

            let k_target = top.concat_rows(&[&bottom]);

            drop(top);
            drop(bottom);

            let (b_matrix, trapdoor_first_path, trapdoor_second_path) = if bit == 0 {
                (&b_next_0, b_next_0_trapdoor_first, b_next_0_trapdoor_second)
            } else {
                (&b_next_1, b_next_1_trapdoor_first, b_next_1_trapdoor_second)
            };

            let trapdoor_first = M::load(trapdoor_first_path);
            let trapdoor_second = M::load(trapdoor_second_path);

            info!("Loaded k_preimage for input {} bit {}", idx + 1, bit);
            log_mem();

            let k_matrix = sampler_trapdoor.preimage(
                &params,
                b_matrix,
                &trapdoor_first,
                &trapdoor_second,
                &k_target,
            );

            info!("Sampled k_preimage for input {} bit {}", idx + 1, bit);
            log_mem();

            drop(trapdoor_first);
            drop(trapdoor_second);
            drop(k_target);

            let k_path = fs_dir_path.join(format!("k_preimage_{}_{}", bit, idx));
            k_matrix.store(&k_path);
            drop(k_matrix);

            k_path
        };

        let kp = || join!(|| k_preimage(0), || k_preimage(1));

        let (mp, (np, kp)) = join!(mp, || join!(np, kp));

        n_preimages_paths.push(np);
        m_preimages_paths.push(mp);
        k_preimages_paths.push(kp);

        drop(b_next_0);
        drop(b_next_1);
        drop(b_next_star);
        drop(ub_star);

        log_mem();
    }

    info!("Finished all preimage calculations");
    log_mem();

    let a_decomposed_polys = public_data.a_rlwe_bar.get_column_matrix(0).decompose().get_column(0);
    let enc_hardcoded_key = {
        let e = sampler_uniform.sample_uniform(
            &params,
            1,
            1,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        let scale = M::P::from_const(&params, &<M::P as Poly>::Elem::half_q(&params.modulus()));
        let result =
            &t_bar_matrix * &public_data.a_rlwe_bar + &e - &(&hardcoded_key_matrix * &scale);
        drop(e);
        result
    };
    let enc_hardcoded_key_polys = enc_hardcoded_key.get_column_matrix(0).decompose().get_column(0);
    let final_circuit = build_final_step_circuit::<_, BggPublicKey<M>>(
        &params,
        &a_decomposed_polys,
        &enc_hardcoded_key_polys,
        public_circuit.clone(),
    );

    let final_preimage_target = {
        let one = public_data.pubkeys[obf_params.input_size][0].clone();
        let input = &public_data.pubkeys[obf_params.input_size][1..];

        let eval_outputs = final_circuit.eval(params.as_ref(), one, input);

        let mut eval_outputs_matrix = eval_outputs[0].concat_matrix(&eval_outputs[1..]);
        let unit_vector = identity_2.slice_columns(1, 2);
        eval_outputs_matrix = eval_outputs_matrix * unit_vector.decompose();

        drop(unit_vector);

        (eval_outputs_matrix + public_data.a_prf).concat_rows(&[&M::zero(
            params.as_ref(),
            2,
            packed_output_size,
        )])
    };

    let (_, _, b_final_path) = &bs_path[obf_params.input_size];
    let (_, _, _, _, b_final_trapdoor_first, b_final_trapdoor_second) =
        &b_trapdoors[obf_params.input_size];
    let b_final = M::load(b_final_path);
    let b_final_trapdoor_first = M::load(b_final_trapdoor_first);
    let b_final_trapdoor_second = M::load(b_final_trapdoor_second);

    let final_preimage = sampler_trapdoor.preimage(
        &params,
        &b_final,
        &b_final_trapdoor_first,
        &b_final_trapdoor_second,
        &final_preimage_target,
    );

    drop(b_final);
    drop(b_final_trapdoor_first);
    drop(b_final_trapdoor_second);

    let final_preimage_path = fs_dir_path.join("final_preimage");
    final_preimage.store(&final_preimage_path);
    drop(final_preimage);

    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init]);
        let bs_0_2 = M::load(&bs_path[0].2);
        let s_b = &s_connect * &bs_0_2;

        drop(s_connect);
        drop(bs_0_2);

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
    let encodings_init = bgg_encode_sampler.sample(&params, &public_data.pubkeys[0], &plaintexts);

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
        hardcoded_key,
        #[cfg(test)]
        final_preimage_target,
    }
}
