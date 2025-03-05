use super::utils::*;
use super::Obfuscation;
use crate::bgg::sampler::{BGGEncodingSampler, BGGPublicKeySampler};
use crate::poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams};
use rand::{Rng, RngCore};
use std::sync::Arc;

pub fn obfuscate<M, S, R>(
    params: <M::P as Poly>::Params,
    mut sampler: S,
    input_size: usize,
    error_gauss_sigma: f64,
    rng: &mut R,
) -> Obfuscation<M>
where
    M: PolyMatrix,
    S: PolyUniformSampler<M = M> + PolyHashSampler<[u8; 32], M = M> + PolyTrapdoorSampler<M = M>,
    R: RngCore,
{
    let hash_key = rng.random::<[u8; 32]>();
    sampler.set_key(hash_key);
    let sampler = Arc::new(sampler);
    let params = Arc::new(params);
    let dim = params.as_ref().ring_dimension() as usize;
    let packed_input_size = input_size.div_ceil(dim);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone());
    let public_data = PublicSampledData::sample(
        params.as_ref(),
        sampler.clone(),
        &bgg_pubkey_sampler,
        input_size,
        packed_input_size,
    );
    let s_bar = sampler.sample_uniform(&params, 1, 1, DistType::BitDist).entry(1, 1).clone();
    let bgg_encode_sampler =
        BGGEncodingSampler::new(params.as_ref(), &s_bar, sampler.clone(), error_gauss_sigma);
    let s_init = &bgg_encode_sampler.secret_vec;
    let t_bar = sampler.sample_uniform(&params, 1, 1, DistType::BitDist);
    let minus_one_poly =
        M::from_poly_vec_row(params.as_ref(), vec![M::P::const_minus_one(params.as_ref())]);
    let t = t_bar.concat_columns(&[minus_one_poly]);
    let log_q = params.as_ref().modulus_bits();
    let a_fhe_bar = &public_data.a_fhe_bar;
    let b_fhe = {
        let muled = t_bar * a_fhe_bar;
        let error = sampler.sample_uniform(
            &params,
            1,
            2 * log_q,
            DistType::GaussDist { sigma: error_gauss_sigma },
        );
        muled + error
    };
    // let a_fhe = a_fhe_bar.concat_rows(&[b_fhe]);

    // // let bgg_pubkeys_input =
    // //     bgg_pubkey_sampler.sample(TAG_BGG_PUBKEY_INPUT, packed_input_size + 1)?;
    // // let bgg_pubkeys_fhe_key = bgg_pubkey_sampler.sample(TAG_BGG_PUBKEY_FHEKEY, 2)?;
    let zero_plaintexts: Vec<M::P> =
        (0..packed_input_size).map(|_| M::P::const_zero(params.as_ref())).collect();
    let mut one_and_zeros = vec![<M::P as Poly>::const_one(params.as_ref())];
    one_and_zeros.extend(zero_plaintexts);
    let encode_input =
        bgg_encode_sampler.sample(&params, &public_data.pubkeys_input[0], &one_and_zeros, true);
    let encode_fhe_key =
        bgg_encode_sampler.sample(&params, &public_data.pubkeys_fhe_key[0], &t.get_row(1), false);
    let mut bs = vec![];
    let mut b_trapdoors = vec![];
    for _ in 0..=input_size {
        let (b_0_trapdoor, b_0) = sampler.trapdoor(&params);
        let (b_1_trapdoor, b_1) = sampler.trapdoor(&params);
        let (b_star_trapdoor, b_star) = sampler.trapdoor(&params);
        bs.push((b_0, b_1, b_star));
        b_trapdoors.push((b_0_trapdoor, b_1_trapdoor, b_star_trapdoor));
    }
    let m_b = 2 + log_q;
    let p_init = {
        let s_connect = s_init.concat_columns(&[s_init.clone()]);
        let s_b = s_connect * &bs[0].2;
        let error = sampler.sample_uniform(
            &params,
            1,
            m_b,
            DistType::GaussDist { sigma: error_gauss_sigma },
        );
        s_b + error
    };
    let identity_2 = M::identity(params.as_ref(), 2, None);
    let u_0 = identity_2.concat_diag(&[public_data.r_0.clone()]);
    let u_1 = identity_2.concat_diag(&[public_data.r_1.clone()]);
    let u_star = {
        let zeros = M::zero(params.as_ref(), 2, 4);
        let identities = identity_2.concat_columns(&[identity_2.clone()]);
        zeros.concat_rows(&[identities])
    };
    let gadget_2 = M::gadget_matrix(params.as_ref(), 2);
    let (mut m_preimages, mut n_preimages, mut k_preimages) = (vec![], vec![], vec![]);
    for idx in 0..input_size {
        let (_, _, b_cur_star) = &bs[idx];
        let (b_next_0, b_next_1, b_next_star) = &bs[idx + 1];
        let (_, _, b_cur_star_trapdoor) = &b_trapdoors[idx];
        let (b_next_0_trapdoor, b_next_1_trapdoor, _) = &b_trapdoors[idx + 1];
        let m_0 = {
            let ub = u_0.clone() * b_next_0;
            sampler.preimage(params.as_ref(), b_cur_star_trapdoor, b_cur_star, &ub)
        };
        let m_1 = {
            let ub = u_1.clone() * b_next_1;
            sampler.preimage(params.as_ref(), b_cur_star_trapdoor, b_cur_star, &ub)
        };
        m_preimages.push((m_0, m_1));

        let ub_star = u_star.clone() * b_next_star;
        let n_0 = sampler.preimage(&params, b_next_0_trapdoor, b_next_0, &ub_star);
        let n_1 = sampler.preimage(&params, b_next_1_trapdoor, b_next_1, &ub_star);
        n_preimages.push((n_0, n_1));

        let mut ks = vec![];
        for bit in 0..1 {
            let (t_input, t_fhe_key) = if bit == 0 { &public_data.t_0 } else { &public_data.t_1 };
            let at_input = public_data.pubkeys_input[idx][0]
                .concat_matrix(&public_data.pubkeys_input[idx][1..])
                * t_input;
            let at_fhe_key = public_data.pubkeys_fhe_key[idx][0]
                .concat_matrix(&public_data.pubkeys_fhe_key[idx][1..])
                * t_fhe_key;
            let former = at_input.concat_columns(&[at_fhe_key]);
            let inserted_poly_index = idx / dim;
            let inserted_coeff_index = idx % dim;
            let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
            let mut coeffs = vec![zero_coeff; dim];
            coeffs[inserted_coeff_index] = <M::P as Poly>::Elem::one(&params.modulus());
            let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
            let inserted_poly_gadget = {
                let zero = <M::P as Poly>::const_zero(params.as_ref());
                let mut polys = vec![];
                for _ in 0..inserted_poly_index {
                    polys.push(zero.clone());
                }
                polys.push(inserted_poly);
                for _ in inserted_poly_index + 1..packed_input_size {
                    polys.push(zero.clone());
                }
                M::from_poly_vec_row(params.as_ref(), polys) * &gadget_2
            };
            let a_input_next = public_data.pubkeys_input[idx + 1][0]
                .concat_matrix(&public_data.pubkeys_input[idx + 1][1..])
                - &inserted_poly_gadget;
            let latter = a_input_next.concat_columns(&[public_data.pubkeys_fhe_key[idx + 1][0]
                .concat_matrix(&public_data.pubkeys_fhe_key[idx + 1][1..])]);
            let k_target = former.concat_rows(&[latter]);
            let b_matrix = if bit == 0 { b_next_0 } else { b_next_1 };
            let trapdoor = if bit == 0 { b_next_0_trapdoor } else { b_next_1_trapdoor };
            let k = sampler.preimage(&params, trapdoor, b_matrix, &k_target);
            ks.push(k);
        }
        k_preimages.push((ks[0].clone(), ks[1].clone()));
    }

    // here we support only inner product between the fhe secret key t and the input x
    let mut ip_pubkey = None;
    for idx in 0..packed_input_size {
        let muled = public_data.pubkeys_input[input_size][idx].clone()
            * &public_data.pubkeys_fhe_key[input_size][0];
        match ip_pubkey {
            None => {
                ip_pubkey = Some(muled);
            }
            Some(ip) => {
                ip_pubkey = Some(ip + muled);
            }
        }
    }
    let ip_pubkey = ip_pubkey.unwrap();
    let final_preimage_target =
        ip_pubkey.matrix.concat_rows(&[M::zero(params.as_ref(), 2, ip_pubkey.matrix.col_size())]);
    let (_, _, b_final) = &bs[input_size];
    let (_, _, b_final_trapdoor) = &b_trapdoors[input_size];
    let final_preimage =
        sampler.preimage(&params, b_final_trapdoor, b_final, &final_preimage_target);

    Obfuscation {
        hash_key,
        b_fhe,
        encode_input,
        encode_fhe_key,
        p_init,
        m_preimages,
        n_preimages,
        k_preimages,
        final_preimage,
    }
}
