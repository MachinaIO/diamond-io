use super::{circuit::Evaluable, BggEncoding, BggPublicKey};
use crate::poly::{Poly, PolyMatrix, PolyParams};
use rayon::prelude::*;
pub trait DigitsToInt<P: Poly>: Evaluable {
    fn power_of_base(&self, params: &P::Params, k: usize) -> Self;
    fn digits_to_int(digits: &[Self], params: &P::Params) -> Self {
        let log_base_q = params.modulus_digits();
        debug_assert_eq!(digits.len(), log_base_q);

        digits
            .par_iter()
            .enumerate()
            .map(|(i, digit)| digit.power_of_base(params, i))
            .reduce_with(|a, b| a + b)
            .unwrap()
    }
}

impl<P: Poly> DigitsToInt<P> for P {
    fn power_of_base(&self, params: &P::Params, k: usize) -> Self {
        let power_of_base = P::const_power_of_base(params, k);
        self.clone() * power_of_base
    }
}

impl<M: PolyMatrix> DigitsToInt<M::P> for BggPublicKey<M> {
    fn power_of_base(&self, params: &<M::P as crate::poly::Poly>::Params, k: usize) -> Self {
        let scalar = M::P::const_power_of_base(params, k);
        // d+1
        let d1 = self.matrix.row_size();
        let unit_vector = M::unit_column_vector(params, d1, d1 - 1);
        let scalared = unit_vector * scalar;
        let decomposed = scalared.decompose();
        let matrix = self.matrix.clone() * decomposed;
        Self { matrix, reveal_plaintext: self.reveal_plaintext }
    }
}

impl<M: PolyMatrix> DigitsToInt<M::P> for BggEncoding<M> {
    fn power_of_base(&self, params: &<M::P as crate::poly::Poly>::Params, k: usize) -> Self {
        let scalar = M::P::const_power_of_base(params, k);
        // d+1
        let d1 = self.pubkey.matrix.row_size();
        let unit_vector = M::unit_column_vector(params, d1, d1 - 1);
        let scalared = unit_vector * &scalar;
        let decomposed = scalared.decompose();
        let vector = self.vector.clone() * decomposed;
        let pubkey = self.pubkey.power_of_base(params, k);
        let plaintext = self.plaintext.clone().map(|plaintext| plaintext * scalar);
        Self { vector, pubkey, plaintext }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bgg::{
            sampler::{BGGEncodingSampler, BGGPublicKeySampler},
            BggEncoding, BggPublicKey,
        },
        poly::{
            dcrt::{
                element::FinRingElem,
                matrix::DCRTPolyMatrix,
                params::DCRTPolyParams,
                poly::DCRTPoly,
                sampler::{hash::DCRTPolyHashSampler, uniform::DCRTPolyUniformSampler},
            },
            sampler::PolyUniformSampler,
            PolyElem,
        },
        utils::{create_bit_random_poly, create_random_poly},
    };
    use keccak_asm::Keccak256;

    #[test]
    fn test_dcrtpoly_digits_to_int_random() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Generate digit polynomials
        let random_poly = create_random_poly(&params);
        let digit_polys = random_poly.decompose_base(&params);

        // Compute the integer representation
        let result = DCRTPoly::digits_to_int(&digit_polys, &params);

        assert_eq!(result, random_poly, "digits_to_int result does not match expected value");
    }

    #[test]
    fn test_dcrtpoly_digits_to_int_static() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create digit polynomials with known values
        let digit_polys =
            DCRTPoly::from_const(&params, &FinRingElem::constant(&params.modulus(), 13))
                .decompose_base(&params);

        // Compute the integer representation
        let result = DCRTPoly::digits_to_int(&digit_polys, &params);

        // Expected result: 1 + 2 + 0 + 8 = 11
        // In polynomial form, this is a constant polynomial with value 11
        let expected = DCRTPoly::from_const(&params, &FinRingElem::new(13u32, params.modulus()));
        assert_eq!(result, expected, "digits_to_int result does not match expected value 13");
    }

    #[test]
    fn test_bggpublickey_digits_to_int() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys with reveal_plaintext set to true
        let reveal_plaintexts = vec![true; params.modulus_digits() + 1];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Extract digit public keys
        let digit_pubkeys: Vec<BggPublicKey<DCRTPolyMatrix>> =
            pubkeys.into_iter().skip(1).collect();

        // Compute the integer representation
        let result = BggPublicKey::<DCRTPolyMatrix>::digits_to_int(&digit_pubkeys, &params);

        // Verify the result by manually computing the expected value
        let mut expected = digit_pubkeys[0].power_of_base(&params, 0);
        for (i, digit_pubkey) in digit_pubkeys.iter().enumerate().skip(1) {
            expected = expected + digit_pubkey.power_of_base(&params, i);
        }

        assert_eq!(
            result.matrix, expected.matrix,
            "bits_to_int matrix does not match expected value"
        );
        assert_eq!(
            result.reveal_plaintext, expected.reveal_plaintext,
            "bits_to_int reveal_plaintext does not match expected value"
        );
    }

    #[test]
    fn test_bggencoding_digits_to_int_static() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_pubkey_sampler =
            BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        let uniform_sampler = DCRTPolyUniformSampler::new();

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = vec![true; params.modulus_digits() + 1]; // +1 for the one encoding
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts (digit polynomials)
        let secrets = vec![create_bit_random_poly(&params); d];
        let int_poly = DCRTPoly::from_const(&params, &FinRingElem::constant(&params.modulus(), 13));
        let plaintexts = int_poly.decompose_base(&params);

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);

        // Extract digit encodings (skip the first one which is the one encoding)
        let digit_encodings: Vec<BggEncoding<DCRTPolyMatrix>> =
            encodings.into_iter().skip(1).collect();

        // Compute the integer representation
        let result = BggEncoding::<DCRTPolyMatrix>::digits_to_int(&digit_encodings, &params);
        assert_eq!(result.plaintext.unwrap(), int_poly);
    }

    #[test]
    fn test_bggencoding_digits_to_int_random() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_pubkey_sampler =
            BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        let uniform_sampler = DCRTPolyUniformSampler::new();

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = vec![true; params.modulus_digits() + 1]; // +1 for the one encoding
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts (digit polynomials)
        let secrets = vec![create_bit_random_poly(&params); d];
        let int_poly = create_bit_random_poly(&params);
        let plaintexts = int_poly.decompose_base(&params);

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);

        // Extract digit encodings (skip the first one which is the one encoding)
        let digit_encodings: Vec<BggEncoding<DCRTPolyMatrix>> =
            encodings.into_iter().skip(1).collect();

        // Compute the integer representation
        let result = BggEncoding::<DCRTPolyMatrix>::digits_to_int(&digit_encodings, &params);

        // Verify the result by manually computing the expected value
        let mut expected = digit_encodings[0].power_of_base(&params, 0);
        for (i, digit_encoding) in digit_encodings.iter().enumerate().skip(1) {
            expected = expected + digit_encoding.power_of_base(&params, i);
        }

        let expected_pubkey = BggPublicKey::digits_to_int(&pubkeys[1..], &params);

        assert_eq!(
            result.vector, expected.vector,
            "bits_to_int vector does not match expected value"
        );
        assert_eq!(
            result.pubkey, expected_pubkey,
            "bits_to_int pubkey.matrix does not match expected value"
        );

        assert_eq!(result.plaintext.unwrap(), int_poly);
    }
}
