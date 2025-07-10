use super::{circuit::Evaluable, BggPublicKey};
use crate::{
    bgg::lut::public_lut::PublicLut,
    poly::{Poly, PolyMatrix},
    utils::timed_read,
};
use rayon::prelude::*;
use std::{
    ops::{Add, Mul, Sub},
    path::PathBuf,
    time::Duration,
};
use tracing::{info, warn};

#[derive(Debug, Clone)]
pub struct BggEncoding<M: PolyMatrix> {
    pub vector: M,
    pub pubkey: BggPublicKey<M>,
    pub plaintext: Option<<M as PolyMatrix>::P>,
}

impl<M: PolyMatrix> BggEncoding<M> {
    pub fn new(
        vector: M,
        pubkey: BggPublicKey<M>,
        plaintext: Option<<M as PolyMatrix>::P>,
    ) -> Self {
        Self { vector, pubkey, plaintext }
    }

    pub fn concat_vector(&self, others: &[Self]) -> M {
        self.vector.concat_columns(&others.par_iter().map(|x| &x.vector).collect::<Vec<_>>()[..])
    }

    /// Writes the encoding with id to files under the given directory.
    pub async fn write_to_files<P: AsRef<std::path::Path> + Send + Sync>(
        &self,
        dir_path: P,
        id: &str,
    ) {
        // Write the vector
        self.vector.write_to_files(&dir_path, &format!("{id}_vector")).await;

        // Write the pubkey
        self.pubkey.write_to_files(&dir_path, &format!("{id}_pubkey")).await;

        // Write the plaintext component if it exists
        if let Some(plaintext) = &self.plaintext {
            plaintext.write_to_file(&dir_path, &format!("{id}_plaintext")).await;
        }
    }

    /// Reads an encoding with id from files under the given directory.
    pub fn read_from_files<P: AsRef<std::path::Path> + Send + Sync>(
        params: &<M::P as Poly>::Params,
        d1: usize,
        log_base_q: usize,
        dir_path: P,
        id: &str,
        reveal_plaintext: bool,
    ) -> Self {
        let ncol = d1 * log_base_q;

        // Read the vector
        let vector = M::read_from_files(params, 1, ncol, &dir_path, &format!("{id}_vector"));

        // Read the pubkey
        let pubkey = BggPublicKey::read_from_files(
            params,
            d1,
            ncol,
            &dir_path,
            &format!("{id}_pubkey"),
            reveal_plaintext,
        );

        // If reveal_plaintext is true, read the plaintext
        let plaintext = if reveal_plaintext {
            Some(M::P::read_from_file(params, &dir_path, &format!("{id}_plaintext")))
        } else {
            None
        };

        Self { vector, pubkey, plaintext }
    }
}

impl<M: PolyMatrix> Add for BggEncoding<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggEncoding<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        let vector = self.vector + &other.vector;
        let pubkey = self.pubkey + &other.pubkey;
        let plaintext = match (self.plaintext, other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a + b),
            _ => None,
        };
        Self { vector, pubkey, plaintext }
    }
}

impl<M: PolyMatrix> Sub for BggEncoding<M> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self - &other
    }
}

impl<M: PolyMatrix> Sub<&Self> for BggEncoding<M> {
    type Output = Self;
    fn sub(self, other: &Self) -> Self {
        let vector = self.vector - &other.vector;
        let pubkey = self.pubkey - &other.pubkey;
        let plaintext = match (self.plaintext, other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a - b),
            _ => None,
        };
        Self { vector, pubkey, plaintext }
    }
}

impl<M: PolyMatrix> Mul for BggEncoding<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggEncoding<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        if self.plaintext.is_none() {
            panic!("Unknown plaintext for the left-hand input of multiplication");
        }
        let decomposed_b = other.pubkey.matrix.decompose();
        let first_term = self.vector * decomposed_b.clone();
        let second_term = other.vector.clone() * self.plaintext.as_ref().unwrap();
        let new_vector = first_term + second_term;
        let new_plaintext = match (self.plaintext, other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a * b),
            _ => None,
        };

        let new_pubkey = BggPublicKey {
            matrix: self.pubkey.matrix * decomposed_b,
            reveal_plaintext: self.pubkey.reveal_plaintext & other.pubkey.reveal_plaintext,
        };
        Self { vector: new_vector, pubkey: new_pubkey, plaintext: new_plaintext }
    }
}

impl<M: PolyMatrix> Evaluable for BggEncoding<M> {
    type Params = <M::P as Poly>::Params;
    type Matrix = M;

    fn rotate(self, params: &Self::Params, shift: usize) -> Self {
        let rotate_poly = <M::P>::const_rotate_poly(params, shift);
        let vector = self.vector.clone() * &rotate_poly;
        let pubkey = self.pubkey.rotate(params, shift);
        let plaintext = self.plaintext.clone().map(|plaintext| plaintext * rotate_poly);
        Self { vector, pubkey, plaintext }
    }

    fn from_digits(params: &Self::Params, one: &Self, digits: &[u32]) -> Self {
        let const_poly =
            <M::P as Evaluable>::from_digits(params, &<M::P>::const_one(params), digits);
        let vector = one.vector.clone() * &const_poly;
        let pubkey = BggPublicKey::from_digits(params, &one.pubkey, digits);
        let plaintext = one.plaintext.clone().map(|plaintext| plaintext * const_poly);
        Self { vector, pubkey, plaintext }
    }

    fn public_lookup(
        self,
        params: &Self::Params,
        plt: &mut PublicLut<M>,
        helper_lookup: Option<(M, PathBuf, usize, usize)>,
    ) -> Self {
        let c_z = &self.plaintext.clone().expect("the BGG encoding should revealed plaintext");
        info!("c_z {:?}", c_z.coeffs());
        let k = c_z.to_const_int();
        info!("k is {}", k);
        let (p_x_l, dir_path, m, m_b) = helper_lookup.expect("BGG encoding's helper needed");
        if let Some((x_k, y_k)) = plt.f.get(&k) {
            let r_k = plt.r_k_s.slice_columns(k * m, (k + 1) * m);
            let l_k = timed_read(
                &format!("L_{k}"),
                || M::read_from_files(params, m_b, m, &dir_path, &format!("L_{k}")),
                &mut Duration::default(),
            );
            let c_lt_k = p_x_l * l_k;
            let pubkey = BggPublicKey::new(plt.a_lt.clone(), self.pubkey.reveal_plaintext);
            info!("x_k {:?}", x_k.coeffs());
            let vector = self.vector * &r_k.decompose() + c_lt_k;
            Self { vector, pubkey, plaintext: Some(y_k.clone()) }
        } else {
            warn!("Public lookup: no row for index k = {}", k);
            self
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bgg::{
            circuit::PolyCircuit,
            lut::{public_lut::PublicLut, utils::p_vector_for_inputs},
            sampler::{BGGEncodingSampler, BGGPublicKeySampler},
            BggEncoding, BggPublicKey,
        },
        poly::{
            dcrt::{
                matrix::base::BaseMatrix,
                params::DCRTPolyParams,
                sampler::{hash::DCRTPolyHashSampler, uniform::DCRTPolyUniformSampler},
                DCRTPoly, DCRTPolyMatrix, DCRTPolyTrapdoorSampler,
            },
            sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
            Poly, PolyMatrix, PolyParams,
        },
        utils::{create_bit_random_poly, create_random_poly, init_tracing},
    };
    use futures::future::join_all;
    use keccak_asm::Keccak256;
    use rand::Rng;
    use serial_test::serial;
    use std::{collections::HashMap, fs, path::Path};
    use tempfile::tempdir;
    use tokio;
    use tracing::info;

    const SIGMA: f64 = 4.578;

    #[tokio::test]
    #[ignore = "file cannot be read"]
    async fn test_encoding_plt_for_dio() {
        init_tracing();

        let tmp_dir = tempdir().unwrap().into_path();

        /* Setup */
        let params = DCRTPolyParams::default();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);

        /* constant Lookup mapping k => (x_k, y_k) */
        let mut f = HashMap::new();
        f.insert(0, (DCRTPoly::const_int(&params, 0), DCRTPoly::const_int(&params, 7)));
        f.insert(1, (DCRTPoly::const_int(&params, 1), DCRTPoly::const_int(&params, 5)));
        f.insert(2, (DCRTPoly::const_int(&params, 2), DCRTPoly::const_int(&params, 6)));
        f.insert(3, (DCRTPoly::const_int(&params, 3), DCRTPoly::const_int(&params, 1)));
        f.insert(4, (DCRTPoly::const_int(&params, 4), DCRTPoly::const_int(&params, 0)));
        f.insert(5, (DCRTPoly::const_int(&params, 5), DCRTPoly::const_int(&params, 3)));
        f.insert(6, (DCRTPoly::const_int(&params, 6), DCRTPoly::const_int(&params, 4)));
        f.insert(7, (DCRTPoly::const_int(&params, 7), DCRTPoly::const_int(&params, 2)));

        /* Obfuscation Step */
        let key: [u8; 32] = rand::random();
        let d = 1;
        let mut handles: Vec<tokio::task::JoinHandle<()>> = Vec::new();
        let input_size = 1;
        let uni = DCRTPolyUniformSampler::new();
        let bgg_pubkey_sampler =
            BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let (b_l_trapdoor, b_l) = trapdoor_sampler.trapdoor(&params, (d + 1) * (1 + input_size));
        info!("b_l ({},{})", b_l.row_size(), b_l.col_size());
        let m = (1 + d) * params.modulus_digits();
        let m_b = (1 + input_size) * (d + 1) * (2 + params.modulus_digits());
        /* BGG+ encoding setup */
        let secrets = uni.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
        // in reality there should be input insertion step that updates the secret s_init to
        let s_x_l = {
            let minus_one_poly = DCRTPoly::const_minus_one(&params);
            let mut secrets = secrets.to_vec();
            secrets.push(minus_one_poly);
            DCRTPolyMatrix::from_poly_vec_row(&params, secrets)
        };
        info!("s_x_L ({},{})", s_x_l.row_size(), s_x_l.col_size());
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let lut = PublicLut::<DCRTPolyMatrix>::new::<
            DCRTPolyUniformSampler,
            DCRTPolyHashSampler<Keccak256>,
        >(&params, d, f, rand::random());

        let a_lt = lut.clone().a_lt;

        // Create a simple circuit with an plt operation
        let mut circuit = PolyCircuit::new();
        {
            let inputs = circuit.input(1);
            let plt_id = circuit.register_public_lookup(lut);
            let plt_gate = circuit.public_lookup_gate(inputs[0], plt_id);
            circuit.output(vec![plt_gate]);
        }

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        assert_eq!(pubkeys.len(), 2);
        let _ = circuit.eval(&params, &pubkeys[0], &[pubkeys[1].clone()], None);
        // after update the target matrix above while evaluation, now can sample preimage L_k
        circuit.preimage_sample_all_lookups(
            &params,
            &b_l,
            &trapdoor_sampler,
            &b_l_trapdoor,
            input_size + 1,
            &tmp_dir,
            &mut handles,
        );
        join_all(handles).await;

        // Create secret and plaintexts
        let plaintexts = vec![DCRTPoly::const_int(&params, 2)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();

        assert_eq!(*enc1.plaintext.as_ref().unwrap(), DCRTPoly::const_int(&params, 2));

        // Evaluate the circuit
        let p_x_l = p_vector_for_inputs(&b_l, plaintexts, &params, 0.0, &s_x_l);
        let result = circuit.eval(&params, &enc_one, &[enc1], Some((p_x_l, tmp_dir, m, m_b)));

        let expected_encodings = bgg_encoding_sampler.sample(
            &params,
            &[BggPublicKey::new(a_lt.clone(), true), BggPublicKey::new(a_lt.clone(), true)],
            &[DCRTPoly::const_int(&params, 6)],
        );
        let expected_enc1 = expected_encodings[1].clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected_enc1.vector);
        assert_eq!(result[0].pubkey.matrix, a_lt.clone());
        assert_eq!(*result[0].plaintext.as_ref().unwrap(), DCRTPoly::const_int(&params, 6));
    }

    #[test]
    fn test_encoding_add() {
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
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &enc_one.clone(), &[enc1.clone(), enc2.clone()], None);

        // Expected result
        let expected = enc1.clone() + enc2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_sub() {
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
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with a Sub operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &enc_one, &[enc1.clone(), enc2.clone()], None);

        // Expected result
        let expected = enc1.clone() - enc2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_mul() {
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
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with a Mul operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &enc_one, &[enc1.clone(), enc2.clone()], None);

        // Expected result
        let expected = enc1.clone() * enc2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_circuit_operations() {
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
        let reveal_plaintexts = [true; 4];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();
        let enc3 = encodings[3].clone();

        // Create a circuit: ((enc1 + enc2)^2) - enc3
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(3);

        // enc1 + enc2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (enc1 + enc2)^2
        let square_gate = circuit.mul_gate(add_gate, add_gate);

        // ((enc1 + enc2)^2) - enc3
        let sub_gate = circuit.sub_gate(square_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval(&params, &enc_one, &[enc1.clone(), enc2.clone(), enc3.clone()], None);

        // Expected result: ((enc1 + enc2)^2) - enc3
        let expected =
            ((enc1.clone() + enc2.clone()) * (enc1.clone() + enc2.clone())) - enc3.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_complex_circuit() {
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
        let reveal_plaintexts = [true; 5];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();
        let enc3 = encodings[3].clone();
        let enc4 = encodings[4].clone();

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = enc1 + enc2, b = enc3 * enc4
        // Level 2: c = a * b, d = enc1 - enc3
        // Level 3: e = c + d
        // Level 4: f = e * e
        // Output: f
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(4);

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // enc1 + enc2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // enc3 * enc4

        // Level 2
        let c = circuit.mul_gate(a, b); // (enc1 + enc2) * (enc3 * enc4)
        let d = circuit.sub_gate(inputs[0], inputs[2]); // enc1 - enc3

        // Level 3
        let e = circuit.add_gate(c, d); // ((enc1 + enc2) * (enc3 * enc4)) + (enc1 - enc3)

        // Level 4
        let f = circuit.mul_gate(e, e); // (((enc1 + enc2) * (enc3 * enc4)) + (enc1 - enc3))^2

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval(
            &params,
            &enc_one,
            &[enc1.clone(), enc2.clone(), enc3.clone(), enc4.clone()],
            None,
        );

        // Expected result: (((enc1 + enc2) * (enc3 * enc4)) + (enc1 - enc3)) * scalar
        let sum1 = enc1.clone() + enc2.clone();
        let prod1 = enc3.clone() * enc4.clone();
        let prod2 = sum1.clone() * prod1;
        let diff = enc1.clone() - enc3.clone();
        let sum2 = prod2 + diff;
        let expected = sum2.clone() * sum2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_register_and_call_sub_circuit() {
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
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a sub-circuit that performs addition and multiplication
        let mut sub_circuit = PolyCircuit::new();
        let sub_inputs = sub_circuit.input(2);

        // Add operation: enc1 + enc2
        let add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);

        // Mul operation: enc1 * enc2
        let mul_gate = sub_circuit.mul_gate(sub_inputs[0], sub_inputs[1]);

        // Set the outputs of the sub-circuit
        sub_circuit.output(vec![add_gate, mul_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::new();
        let main_inputs = main_circuit.input(2);

        // Register the sub-circuit and get its ID
        let sub_circuit_id = main_circuit.register_sub_circuit(sub_circuit);

        // Call the sub-circuit with the main circuit's inputs
        let sub_outputs =
            main_circuit.call_sub_circuit(sub_circuit_id, &[main_inputs[0], main_inputs[1]]);

        // Verify we got two outputs from the sub-circuit
        assert_eq!(sub_outputs.len(), 2);

        // Use the sub-circuit outputs for further computation
        // For example, subtract the multiplication result from the addition result
        let final_gate = main_circuit.sub_gate(sub_outputs[0], sub_outputs[1]);

        // Set the output of the main circuit
        main_circuit.output(vec![final_gate]);

        // Evaluate the main circuit
        let result = main_circuit.eval(&params, &enc_one, &[enc1.clone(), enc2.clone()], None);

        // Expected result: (enc1 + enc2) - (enc1 * enc2)
        let expected = (enc1.clone() + enc2.clone()) - (enc1.clone() * enc2.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_nested_sub_circuits() {
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
        let reveal_plaintexts = [true; 4];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();
        let enc3 = encodings[3].clone();

        // Create the innermost sub-circuit that performs multiplication
        let mut inner_circuit = PolyCircuit::new();
        let inner_inputs = inner_circuit.input(2);
        let mul_gate = inner_circuit.mul_gate(inner_inputs[0], inner_inputs[1]);
        inner_circuit.output(vec![mul_gate]);

        // Create a middle sub-circuit that uses the inner sub-circuit
        let mut middle_circuit = PolyCircuit::new();
        let middle_inputs = middle_circuit.input(3);

        // Register the inner circuit
        let inner_circuit_id = middle_circuit.register_sub_circuit(inner_circuit);

        // Call the inner circuit with the first two inputs
        let inner_outputs = middle_circuit
            .call_sub_circuit(inner_circuit_id, &[middle_inputs[0], middle_inputs[1]]);

        // Add the result of the inner circuit with the third input
        let add_gate = middle_circuit.add_gate(inner_outputs[0], middle_inputs[2]);
        middle_circuit.output(vec![add_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::new();
        let main_inputs = main_circuit.input(3);

        // Register the middle circuit
        let middle_circuit_id = main_circuit.register_sub_circuit(middle_circuit);

        // Call the middle circuit with all inputs
        let middle_outputs = main_circuit
            .call_sub_circuit(middle_circuit_id, &[main_inputs[0], main_inputs[1], main_inputs[2]]);

        // Use the output for square
        let scalar_mul_gate = main_circuit.mul_gate(middle_outputs[0], middle_outputs[0]);

        // Set the output of the main circuit
        main_circuit.output(vec![scalar_mul_gate]);

        // Evaluate the main circuit
        let result =
            main_circuit.eval(&params, &enc_one, &[enc1.clone(), enc2.clone(), enc3.clone()], None);

        // Expected result: ((enc1 * enc2) + enc3)^2
        let expected = ((enc1.clone() * enc2.clone()) + enc3.clone()) *
            ((enc1.clone() * enc2.clone()) + enc3.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[tokio::test]
    #[serial]
    async fn test_encoding_write_read() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        // let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let d1 = d + 1;
        let bgg_pubkey_sampler =
            BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let log_base_q = params.modulus_digits();

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys with different reveal_plaintext flags
        let mut rng = rand::rng();
        let reveal_plaintexts = [
            rng.random::<bool>(),
            rng.random::<bool>(),
            rng.random::<bool>(),
            rng.random::<bool>(),
            rng.random::<bool>(),
        ];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secrets = vec![create_bit_random_poly(&params); d];
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secrets, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);

        // Create a temporary directory for testing
        let test_dir = Path::new("test_encoding_write_read");
        if !test_dir.exists() {
            fs::create_dir(test_dir).unwrap();
        } else {
            // Clean it first to ensure no old files interfere
            fs::remove_dir_all(test_dir).unwrap();
            fs::create_dir(test_dir).unwrap();
        }

        for (idx, encoding) in encodings.iter().enumerate() {
            // Write the encoding to files
            let id = format!("test_encoding_{}", idx);
            encoding.write_to_files(test_dir, &id).await;

            // Read the encoding from files and verify it matches the original
            if idx == 0 || (idx > 0 && reveal_plaintexts[idx - 1]) {
                // first encoding (for one) is always true
                let read_enc: BggEncoding<BaseMatrix<DCRTPoly>> =
                    BggEncoding::read_from_files(&params, d1, log_base_q, test_dir, &id, true);
                assert_eq!(read_enc.vector, encoding.vector);
                assert_eq!(read_enc.pubkey.matrix, encoding.pubkey.matrix);
                assert_eq!(
                    read_enc.plaintext.as_ref().unwrap(),
                    encoding.plaintext.as_ref().unwrap()
                );
            } else {
                let read_enc: BggEncoding<BaseMatrix<DCRTPoly>> =
                    BggEncoding::read_from_files(&params, d1, log_base_q, test_dir, &id, false);
                assert_eq!(read_enc.vector, encoding.vector);
                assert_eq!(read_enc.pubkey.matrix, encoding.pubkey.matrix);
                assert!(read_enc.plaintext.is_none());
            }
        }
        // Clean up the test directory
        std::fs::remove_dir_all(test_dir).unwrap();
    }
}
