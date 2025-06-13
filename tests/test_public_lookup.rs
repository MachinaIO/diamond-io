use diamond_io::poly::{
    dcrt::{
        DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler,
    },
    sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
    Poly, PolyMatrix, PolyParams,
};
use itertools::Itertools;

#[test]
fn test_public_lookup() {
    /*
       Get p vector
    */
    let packed_input_size = 3;
    let d = 2;
    let params = DCRTPolyParams::default();
    let log_base_q = params.modulus_digits();
    let sampler_uniform = DCRTPolyUniformSampler::new();
    let trapdoor_sigma = 4.578;
    let p_sigma = 1.0;
    let one = DCRTPoly::const_one(&params);
    let t_bar = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
    let minus_t_bar = -t_bar.entry(0, 0);
    let s_bars = sampler_uniform.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
    let mut plaintexts = vec![one];
    plaintexts
        .extend((0..packed_input_size - 1).map(|_| DCRTPoly::const_zero(&params)).collect_vec());
    plaintexts.push(minus_t_bar);
    let encoded_bits = DCRTPolyMatrix::from_poly_vec_row(&params, plaintexts);
    let s_init = {
        let minus_one_poly = DCRTPoly::const_minus_one(&params);
        let mut secrets = s_bars.to_vec();
        secrets.push(minus_one_poly);
        DCRTPolyMatrix::from_poly_vec_row(&params, secrets)
    };
    let s_connect = encoded_bits.tensor(&s_init);
    let sampler_trapdoor = DCRTPolyTrapdoorSampler::new(&params, trapdoor_sigma);
    let (mut b_star_trapdoor_cur, mut b_star_cur) =
        sampler_trapdoor.trapdoor(&params, (1 + packed_input_size) * (d + 1));
    let m_b = (1 + packed_input_size) * (d + 1) * (2 + log_base_q);
    let s_b = s_connect * &b_star_cur;
    let p_init_error =
        sampler_uniform.sample_uniform(&params, 1, m_b, DistType::GaussDist { sigma: p_sigma });
    let p_init = s_b + p_init_error;

    /*
       Public lookup table preimage
    */
}
