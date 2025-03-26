#[cfg(test)]
mod test {
    use diamond_io::{
        poly::{
            dcrt::{DCRTPolyParams, DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler},
            sampler::{PolyTrapdoorSampler, PolyUniformSampler},
        },
        utils::{init_tracing, log_mem},
    };

    #[test]
    fn test_preimage_sampling_real_params() {
        init_tracing();
        let params = DCRTPolyParams::new(8192, 7, 51);
        let d = 4;

        // gen trapdoor
        let sampler_trapdoor = DCRTPolyTrapdoorSampler::new();
        log_mem("Before trapdoor gen");
        let (trapdoor, public_matrix) = sampler_trapdoor.trapdoor(&params, d);
        log_mem("After trapdoor gen");

        // gen a dxd target matrix from uniform distribution
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(
            &params,
            d,
            d,
            diamond_io::poly::sampler::DistType::BitDist,
        );

        // sample preimage
        log_mem("Before preimage gen");
        let preimage = sampler_trapdoor.preimage(&params, &trapdoor, &public_matrix, &target);
        log_mem("After preimage gen");
    }
}
