use mxx::{circuit::PolyCircuit, lookup::PublicLut, poly::dcrt::poly::DCRTPoly};

pub enum BenchCircuit {
    AddMul(PolyCircuit<DCRTPoly>),
    AddMulVerify(PolyCircuit<DCRTPoly>),
    Plt(PolyCircuit<DCRTPoly>),
    PltVerify(PolyCircuit<DCRTPoly>),
}

impl BenchCircuit {
    pub fn new_add_mul(add_n: usize, mul_n: usize, log_base_q: usize) -> Self {
        let mut public_circuit = PolyCircuit::new();

        // inputs: BaseDecompose(ct), eval_input
        // outputs: BaseDecompose(ct) * acc
        {
            let inputs = public_circuit.input((2 * log_base_q) + 1);
            let mut outputs = vec![];
            let eval_input = inputs[2 * log_base_q];

            // compute acc according to add_n and mul_n logic
            let mut acc = if add_n == 0 {
                public_circuit.const_one_gate()
            } else {
                public_circuit.const_zero_gate()
            };
            for _ in 0..add_n {
                acc = public_circuit.add_gate(acc, eval_input);
            }
            for _ in 0..mul_n {
                acc = public_circuit.mul_gate(acc, eval_input);
            }

            // compute the output
            for ct_input in inputs[0..2 * log_base_q].iter() {
                let muled = public_circuit.mul_gate(*ct_input, acc);
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

            // compute acc according to add_n and mul_n logic
            let mut acc = if add_n == 0 {
                public_circuit.const_one_gate()
            } else {
                public_circuit.const_zero_gate()
            };
            for _ in 0..add_n {
                acc = public_circuit.add_gate(acc, eval_input);
            }
            for _ in 0..mul_n {
                acc = public_circuit.mul_gate(acc, eval_input);
            }

            // compute the output
            for ct_input in inputs[0..1].iter() {
                let muled = public_circuit.mul_gate(*ct_input, acc);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        Self::AddMulVerify(public_circuit)
    }

    pub fn new_plt(log_base_q: usize, lut: PublicLut<DCRTPoly>) -> Self {
        let mut public_circuit = PolyCircuit::new();

        // inputs: BaseDecompose(ct), eval_input
        // outputs: BaseDecompose(ct) * acc
        {
            let inputs = public_circuit.input((2 * log_base_q) + 1);
            let mut outputs = vec![];
            let eval_input = inputs[2 * log_base_q];

            // register and call lookup
            let plt_id = public_circuit.register_public_lookup(lut);
            let plt_out = public_circuit.public_lookup_gate(eval_input, plt_id);

            // compute the output
            for ct_input in inputs[0..2 * log_base_q].iter() {
                let muled = public_circuit.mul_gate(*ct_input, plt_out);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        Self::Plt(public_circuit)
    }

    pub fn new_plt_verify() -> Self {
        let mut public_circuit = PolyCircuit::new();

        // inputs: BaseDecompose(ct), eval_input
        // outputs: BaseDecompose(ct) * acc
        {
            let inputs = public_circuit.input(2);
            let mut outputs = vec![];
            let eval_input = inputs[1];

            // compute the output
            for ct_input in inputs[0..1].iter() {
                let muled = public_circuit.mul_gate(*ct_input, eval_input);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        Self::PltVerify(public_circuit)
    }

    pub fn as_poly_circuit(self) -> PolyCircuit<DCRTPoly> {
        match self {
            BenchCircuit::AddMul(poly_circuit) => poly_circuit,
            BenchCircuit::AddMulVerify(poly_circuit) => poly_circuit,
            BenchCircuit::Plt(poly_circuit) => poly_circuit,
            BenchCircuit::PltVerify(poly_circuit) => poly_circuit,
        }
    }
}
