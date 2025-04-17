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
            let eval_input = inputs[1];
            let hardcoded_key = inputs[0];
            let num = add_n.checked_div(mul_n).unwrap_or(1);
            for _ in 0..num {
                let added = public_circuit.add_gate(hardcoded_key, eval_input);
                outputs.push(added);
            }
            let num = mul_n.checked_div(add_n).unwrap_or(1);
            for _ in 0..num {
                let muled = public_circuit.mul_gate(hardcoded_key, eval_input);
                outputs.push(muled);
            }
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
