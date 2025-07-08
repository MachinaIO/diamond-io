use super::circuit::Evaluable;
use crate::{
    bgg::lut::public_lut::PublicLut,
    poly::{Poly, PolyMatrix},
    utils::debug_mem,
};
use rayon::prelude::*;
use std::{
    ops::{Add, Mul, Sub},
    path::PathBuf,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BggPublicKey<M: PolyMatrix> {
    pub matrix: M,
    pub reveal_plaintext: bool,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn new(matrix: M, reveal_plaintext: bool) -> Self {
        Self { matrix, reveal_plaintext }
    }

    pub fn concat_matrix(&self, others: &[Self]) -> M {
        self.matrix.concat_columns(&others.par_iter().map(|x| &x.matrix).collect::<Vec<_>>()[..])
    }

    /// Writes the public key with id to files under the given directory.
    pub async fn write_to_files<P: AsRef<std::path::Path> + Send + Sync>(
        &self,
        dir_path: P,
        id: &str,
    ) {
        self.matrix.write_to_files(dir_path, id).await;
    }

    /// Reads a public of given rows and cols with id from files under the given directory.
    pub fn read_from_files<P: AsRef<std::path::Path> + Send + Sync>(
        params: &<M::P as Poly>::Params,
        nrow: usize,
        ncol: usize,
        dir_path: P,
        id: &str,
        reveal_plaintext: bool,
    ) -> Self {
        let matrix = M::read_from_files(params, nrow, ncol, dir_path, id);
        Self { matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Add for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix: self.matrix + &other.matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Sub for BggPublicKey<M> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self - &other
    }
}

impl<M: PolyMatrix> Sub<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn sub(self, other: &Self) -> Self {
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix: self.matrix - &other.matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Mul for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        debug_mem(format!("BGGPublicKey::mul {:?}, {:?}", self.matrix.size(), other.matrix.size()));
        let decomposed = other.matrix.decompose();
        debug_mem("BGGPublicKey::mul decomposed");
        let matrix = self.matrix * decomposed;
        debug_mem("BGGPublicKey::mul matrix multiplied");
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        debug_mem("BGGPublicKey::mul reveal_plaintext");
        Self { matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Evaluable for BggPublicKey<M> {
    type Params = <M::P as Poly>::Params;
    type Matrix = M;

    fn rotate(self, params: &Self::Params, shift: usize) -> Self {
        debug_mem(format!("BGGPublicKey::rotate {:?}, {:?}", self.matrix.size(), shift));
        let rotate_poly = <M::P>::const_rotate_poly(params, shift);
        debug_mem("BGGPublicKey::rotate rotate_poly");
        let matrix = self.matrix.clone() * rotate_poly;
        debug_mem("BGGPublicKey::rotate matrix multiplied");
        Self { matrix, reveal_plaintext: self.reveal_plaintext }
    }

    fn from_digits(params: &Self::Params, one: &Self, digits: &[u32]) -> Self {
        debug_mem(format!("BGGPublicKey::from_digits {:?}, {:?}", one.matrix.size(), digits.len()));
        let const_poly =
            <M::P as Evaluable>::from_digits(params, &<M::P>::const_one(params), digits);
        debug_mem("BGGPublicKey::from_digits const_poly");
        let matrix = one.matrix.clone() * const_poly;
        debug_mem("BGGPublicKey::from_digits matrix multiplied");
        Self { matrix, reveal_plaintext: one.reveal_plaintext }
    }

    fn public_lookup(
        self,
        _: &Self::Params,
        plt: &PublicLut<Self::Matrix>,
        _: Option<Self::Matrix>,
        _: usize,
        _: PathBuf,
    ) -> Self {
        Self { matrix: plt.a_lt.clone(), reveal_plaintext: self.reveal_plaintext }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bgg::{
            circuit::PolyCircuit, lut::public_lut::PublicLut, sampler::BGGPublicKeySampler,
            BggPublicKey,
        },
        poly::{
            dcrt::{
                params::DCRTPolyParams, DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix,
                DCRTPolyUniformSampler,
            },
            Poly,
        },
    };
    use keccak_asm::Keccak256;
    use rand::Rng;
    use serial_test::serial;
    use std::{collections::HashMap, fs, path::Path};
    use tokio;

    #[test]
    fn test_pubkey_plt() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();
        /* Lookup mapping k => (x_k, y_k) */
        let mut f = HashMap::new();
        f.insert(0, (DCRTPoly::const_int(&params, 0), DCRTPoly::const_int(&params, 7)));
        f.insert(1, (DCRTPoly::const_int(&params, 1), DCRTPoly::const_int(&params, 5)));
        f.insert(2, (DCRTPoly::const_int(&params, 2), DCRTPoly::const_int(&params, 6)));
        f.insert(3, (DCRTPoly::const_int(&params, 3), DCRTPoly::const_int(&params, 1)));
        f.insert(4, (DCRTPoly::const_int(&params, 4), DCRTPoly::const_int(&params, 0)));
        f.insert(5, (DCRTPoly::const_int(&params, 5), DCRTPoly::const_int(&params, 3)));
        f.insert(6, (DCRTPoly::const_int(&params, 6), DCRTPoly::const_int(&params, 4)));
        f.insert(7, (DCRTPoly::const_int(&params, 7), DCRTPoly::const_int(&params, 2)));

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(1);
        let lut = PublicLut::<DCRTPolyMatrix>::new::<
            DCRTPolyUniformSampler,
            DCRTPolyHashSampler<Keccak256>,
        >(&params, d, f);
        let a_lt = lut.a_lt.clone();
        let plt_id = circuit.register_public_lookup(lut);
        let plt_gate = circuit.public_lookup_gate(inputs[0], plt_id);
        circuit.output(vec![plt_gate]);

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1], None, Some("io_plt".into()));

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, a_lt);
    }

    #[test]
    fn test_pubkey_add() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], None, None);

        // Expected result
        let expected = pk1.clone() + pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_sub() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with a Sub operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], None, None);

        // Expected result
        let expected = pk1.clone() - pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with a Mul operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], None, None);

        // Expected result
        let expected = pk1.clone() * pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_circuit_operations() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 4];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();

        // Create a circuit: ((pk1 + pk2)^2) - pk3
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(3);

        // pk1 + pk2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (pk1 + pk2)^2
        let square_gate = circuit.mul_gate(add_gate, add_gate);

        // ((pk1 + pk2)^2) - pk3
        let sub_gate = circuit.sub_gate(square_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone(), pk3.clone()], None, None);

        // Expected result: ((pk1 + pk2)-2) - pk3
        let expected = ((pk1.clone() + pk2.clone()) * (pk1.clone() + pk2.clone())) - pk3.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_complex_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 5];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();
        let pk4 = pubkeys[4].clone();

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = pk1 + pk2, b = pk3 * pk4
        // Level 2: c = a * b, d = pk1 - pk3
        // Level 3: e = c + d
        // Level 4: f = e * e
        // Output: f
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(4);

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // pk1 + pk2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // pk3 * pk4

        // Level 2
        let c = circuit.mul_gate(a, b); // (pk1 + pk2) * (pk3 * pk4)
        let d = circuit.sub_gate(inputs[0], inputs[2]); // pk1 - pk3

        // Level 3
        let e = circuit.add_gate(c, d); // ((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3)

        // Level 4
        let f = circuit.mul_gate(e, e); // (((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3))^2

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval(
            &params,
            &pk_one,
            &[pk1.clone(), pk2.clone(), pk3.clone(), pk4.clone()],
            None,
            None,
        );

        // Expected result: (((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3))^2
        let sum1 = pk1.clone() + pk2.clone();
        let prod1 = pk3.clone() * pk4.clone();
        let prod2 = sum1.clone() * prod1;
        let diff = pk1.clone() - pk3.clone();
        let sum2 = prod2 + diff;
        let expected = sum2.clone() * sum2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_register_and_call_sub_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a sub-circuit that performs addition and multiplication
        let mut sub_circuit = PolyCircuit::new();
        let sub_inputs = sub_circuit.input(2);

        // Add operation: pk1 + pk2
        let add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);

        // Mul operation: pk1 * pk2
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
        let result = main_circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], None, None);

        // Expected result: (pk1 + pk2) - (pk1 * pk2)
        let expected = (pk1.clone() + pk2.clone()) - (pk1.clone() * pk2.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_nested_sub_circuits() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 4];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();

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
        let square_gate = main_circuit.mul_gate(middle_outputs[0], middle_outputs[0]);

        // Set the output of the main circuit
        main_circuit.output(vec![square_gate]);

        // Evaluate the main circuit
        let result = main_circuit.eval(
            &params,
            &pk_one,
            &[pk1.clone(), pk2.clone(), pk3.clone()],
            None,
            None,
        );

        // Expected result: ((pk1 * pk2) + pk3)^2
        let expected = ((pk1.clone() * pk2.clone()) + pk3.clone()) *
            ((pk1.clone() * pk2.clone()) + pk3.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[tokio::test]
    #[serial]
    async fn test_pubkey_write_read() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let mut rng = rand::rng();
        let reveal_plaintexts = [
            rng.random::<bool>(),
            rng.random::<bool>(),
            rng.random::<bool>(),
            rng.random::<bool>(),
        ];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create a temporary directory for testing
        let test_dir = Path::new("test_pubkey_write_read");
        if !test_dir.exists() {
            fs::create_dir(test_dir).unwrap();
        } else {
            // Clean it first to ensure no old files interfere
            fs::remove_dir_all(test_dir).unwrap();
            fs::create_dir(test_dir).unwrap();
        }

        for (idx, pubkey) in pubkeys.iter().enumerate() {
            // Write the public key to files
            let id = format!("test_pubkey_{}", idx);
            pubkey.write_to_files(test_dir, &id).await;

            // Get the size of the original matrix
            let (nrow, ncol) = pubkey.matrix.size();

            let read_pk;

            // Read the public key from files and verify it matches the original
            if idx == 0 {
                read_pk = BggPublicKey::read_from_files(&params, nrow, ncol, test_dir, &id, true);
                assert_eq!(pubkey, &read_pk);
            } else {
                read_pk = BggPublicKey::read_from_files(
                    &params,
                    nrow,
                    ncol,
                    test_dir,
                    &id,
                    reveal_plaintexts[idx - 1],
                );
                assert_eq!(pubkey.matrix, read_pk.matrix);
            }
        }

        std::fs::remove_dir_all(test_dir).unwrap();
    }
}
