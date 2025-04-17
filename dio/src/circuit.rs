use diamond_io::bgg::circuit::PolyCircuit;

pub enum BenchCircuit {
    AddMul(PolyCircuit),
    AddMulVerify(PolyCircuit),
}

impl BenchCircuit {
    pub fn new_add_mul(add_n: usize, mul_n: usize) -> Self {
        let mut public_circuit = PolyCircuit::new();
        let n = add_n + mul_n + 1;

        {
            let inputs = public_circuit.input(n);
            let mut outputs = vec![];
            let eval_input = inputs[n - 1];
            for ct_input in inputs[0..add_n].iter() {
                let added = public_circuit.add_gate(*ct_input, eval_input);
                outputs.push(added);
            }
            for ct_input in inputs[add_n..add_n + mul_n].iter() {
                let muled = public_circuit.mul_gate(*ct_input, eval_input);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        Self::AddMul(public_circuit)
    }

    pub fn new_add_mul_verify(add_n: usize, mul_n: usize) -> Self {
        let mut public_circuit = PolyCircuit::new();
        {
            let inputs = public_circuit.input(2);
            let mut outputs = vec![];
            let input_poly = inputs[1];
            let hardcoded_key = inputs[0];
            // ? if add_n is 0 all gates are mul this works but how can i perform verification if
            // add_n and mul_n is mixed?
            let muled = public_circuit.mul_gate(hardcoded_key, input_poly);
            outputs.push(muled);
            public_circuit.output(outputs);
        }

        Self::AddMulVerify(public_circuit)
    }

    pub fn as_poly_circuit(self) -> PolyCircuit {
        match self {
            BenchCircuit::AddMul(poly_circuit) => poly_circuit,
            BenchCircuit::AddMulVerify(poly_circuit) => poly_circuit,
        }
    }
}
