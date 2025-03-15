use super::{utils::*, Obfuscation, ObfuscationParams};
use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggPublicKey,
    },
    join,
    poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams},
};
use rand::{Rng, RngCore};
use std::{ops::Mul, sync::Arc};
use tracing::info;

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
    info!("gm");
    let public_circuit = &obf_params.public_circuit;
    let dim = obf_params.params.ring_dimension() as usize;
    let log_q = obf_params.params.modulus_bits();
    debug_assert_eq!(public_circuit.num_input(), log_q + obf_params.input_size);
    let hash_key = rng.random::<[u8; 32]>();
    sampler_hash.set_key(hash_key);
    let sampler_uniform = Arc::new(sampler_uniform);
    let sampler_trapdoor = Arc::new(sampler_trapdoor);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(Arc::new(sampler_hash));
    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    let params = Arc::new(obf_params.params);
    let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let s_bar =
        sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist).entry(0, 0).clone();
    info!("s_bar computed");
    let bgg_encode_sampler = BGGEncodingSampler::new(
        params.as_ref(),
        &s_bar,
        sampler_uniform.clone(),
        obf_params.error_gauss_sigma,
    );
    let s_init = &bgg_encode_sampler.secret_vec;
    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    info!("t_bar computed");
    let hardcoded_key_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
    let enc_hardcoded_key = {
        let e = sampler_uniform.sample_uniform(
            &params,
            1,
            1,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        let scale = M::P::from_const(&params, &<M::P as Poly>::Elem::half_q(&params.modulus()));
        &t_bar_matrix * &public_data.a_rlwe_bar + &e - &(&hardcoded_key_matrix * &scale)
    };

    let enc_hardcoded_key_polys = enc_hardcoded_key.decompose().get_column(0);
    info!("enc_hardcoded_key computed");
    let t_bar = t_bar_matrix.entry(0, 0).clone();
    #[cfg(test)]
    let hardcoded_key = hardcoded_key_matrix.entry(0, 0).clone();

    let mut plaintexts = vec![];
    let zero_plaintexts: Vec<M::P> = (0..obf_params.input_size.div_ceil(dim))
        .map(|_| M::P::const_zero(params.as_ref()))
        .collect();
    plaintexts.extend(zero_plaintexts);
    plaintexts.push(t_bar.clone());
    // let mut input_encoded_polys: Vec<<M as PolyMatrix>::P> =
    //     vec![<M::P as Poly>::const_one(params.as_ref())];
    // input_encoded_polys.extend(enc_hardcoded_key_polys);
    // input_encoded_polys.extend(zero_plaintexts);
    let encodings_init = bgg_encode_sampler.sample(&params, &public_data.pubkeys[0], &plaintexts);
    info!("encodings_init computed");
    // let encode_fhe_key =
    //     bgg_encode_sampler.sample(&params, &public_data.pubkeys_fhe_key[0], &t.get_row(0),
    // false);

    let mut bs = vec![];
    let mut b_trapdoors = vec![];
    for _ in 0..=obf_params.input_size {
        let (b_0_trapdoor, b_0) = sampler_trapdoor.trapdoor(&params, 4);
        let (b_1_trapdoor, b_1) = sampler_trapdoor.trapdoor(&params, 4);
        let (b_star_trapdoor, b_star) = sampler_trapdoor.trapdoor(&params, 4);
        bs.push((b_0, b_1, b_star));
        b_trapdoors.push((b_0_trapdoor, b_1_trapdoor, b_star_trapdoor));
        info!("bs computed");
    }
    let m_b = 4 * (2 + log_q);
    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init]);
        let s_b = s_connect * &bs[0].2;
        let error = sampler_uniform.sample_uniform(
            &params,
            1,
            m_b,
            DistType::GaussDist { sigma: obf_params.error_gauss_sigma },
        );
        s_b + error
    };
    info!("p_init computed");
    let identity_2 = M::identity(params.as_ref(), 2, None);
    let u_0 = identity_2.concat_diag(&[&public_data.r_0]);
    let u_1 = identity_2.concat_diag(&[&public_data.r_1]);
    let u_star = {
        let zeros = M::zero(params.as_ref(), 2, 4);
        let identities = identity_2.concat_columns(&[&identity_2]);
        zeros.concat_rows(&[&identities])
    };
    let gadget_2 = M::gadget_matrix(params.as_ref(), 2);

    let (mut m_preimages, mut n_preimages, mut k_preimages) = (vec![], vec![], vec![]);
    for idx in 0..obf_params.input_size {
        let (_, _, b_cur_star) = &bs[idx];
        let (b_next_0, b_next_1, b_next_star) = &bs[idx + 1];
        let (_, _, b_cur_star_trapdoor) = &b_trapdoors[idx];
        let (b_next_0_trapdoor, b_next_1_trapdoor, _) = &b_trapdoors[idx + 1];

        let m_preimage =
            |a| sampler_trapdoor.preimage(params.as_ref(), b_cur_star_trapdoor, b_cur_star, &a);
        let mp =
            || join!(|| m_preimage(u_0.clone() * b_next_0), || m_preimage(u_1.clone() * b_next_1));

        let ub_star = u_star.clone() * b_next_star;
        let n_preimage = |t, n| sampler_trapdoor.preimage(&params, t, n, &ub_star);
        let np = || {
            join!(|| n_preimage(b_next_0_trapdoor, b_next_0), || n_preimage(
                b_next_1_trapdoor,
                b_next_1
            ))
        };

        let k_preimage = |bit: usize| {
            let rg_decomposed = &public_data.rgs_decomposed[bit];
            let lhs = -public_data.pubkeys[idx][0].concat_matrix(&public_data.pubkeys[idx][1..]);
            let top = lhs.mul_tensor_identity(rg_decomposed, 1 + packed_input_size);
            let inserted_poly_index = 1 + idx / dim;
            let inserted_coeff_index = idx % dim;
            let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            let mut coeffs = vec![zero_coeff; dim];
            coeffs[inserted_coeff_index] = if bit == 0 {
                <M::P as Poly>::Elem::zero(&params.modulus())
            } else {
                <M::P as Poly>::Elem::one(&params.modulus())
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
                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_2)
            };
            let bottom = public_data.pubkeys[idx + 1][0]
                .concat_matrix(&public_data.pubkeys[idx + 1][1..]) -
                &inserted_poly_gadget;
            let k_target = top.concat_rows(&[&bottom]);
            let b_matrix = if bit == 0 { b_next_0 } else { b_next_1 };
            let trapdoor = if bit == 0 { b_next_0_trapdoor } else { b_next_1_trapdoor };
            sampler_trapdoor.preimage(&params, trapdoor, b_matrix, &k_target)
        };
        let kp = || join!(|| k_preimage(0), || k_preimage(1));

        let (mp, (np, kp)) = join!(mp, || join!(np, kp));

        m_preimages.push(mp);
        n_preimages.push(np);
        k_preimages.push(kp);
        info!("m_preimages, n_preimages, k_preimages computed");
    }

    let a_decomposed_polys = public_data.a_rlwe_bar.decompose().get_column(0);
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
        debug_assert_eq!(eval_outputs_matrix.col_size(), packed_output_size);
        (eval_outputs_matrix + public_data.a_prf).concat_rows(&[&M::zero(
            params.as_ref(),
            2,
            packed_output_size,
        )])
    };
    let (_, _, b_final) = &bs[obf_params.input_size];
    let (_, _, b_final_trapdoor) = &b_trapdoors[obf_params.input_size];
    let final_preimage =
        sampler_trapdoor.preimage(&params, b_final_trapdoor, b_final, &final_preimage_target);
    info!("final_preimage computed");
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
        #[cfg(test)]
        bs,
        #[cfg(test)]
        hardcoded_key: hardcoded_key.clone(),
        #[cfg(test)]
        final_preimage_target,
    }
}
