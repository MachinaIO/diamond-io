use clap::ValueEnum;
use diamond_io::bgg::circuit::PolyCircuit;

#[derive(Debug, Clone, ValueEnum)]
pub enum BenchCircuitType {
    And,
    AndVerify,
}

pub enum BenchCircuit {
    AND(PolyCircuit),
    AndVerify(PolyCircuit),
}

impl BenchCircuit {
    pub fn new(circuit_type: BenchCircuitType, n: usize) -> Self {
        let mut public_circuit = PolyCircuit::new();
        match circuit_type {
            BenchCircuitType::And => {
                // inputs: binary polynomial, hardcoded_key, eval_input
                // outputs: binary AND eval_input | hardcoded_key AND eval_input
                {
                    let inputs = public_circuit.input(n);
                    let mut outputs = vec![];
                    let eval_input = inputs[n - 1];
                    for ct_input in inputs[0..(n - 1)].iter() {
                        let muled = public_circuit.and_gate(*ct_input, eval_input);
                        outputs.push(muled);
                    }
                    for ct_input in inputs[0..(n - 1)].iter() {
                        let muled = public_circuit.and_gate(*ct_input, eval_input);
                        outputs.push(muled);
                    }
                    public_circuit.output(outputs);
                }

                Self::AND(public_circuit)
            }
            BenchCircuitType::AndVerify => {
                let inputs = public_circuit.input(n);
                let mut outputs = vec![];
                let eval_input = inputs[n - 1];
                for ct_input in inputs[0..(n - 1)].iter() {
                    let muled = public_circuit.and_gate(*ct_input, eval_input);
                    outputs.push(muled);
                }

                public_circuit.output(outputs);
                Self::AndVerify(public_circuit)
            }
        }
    }

    pub fn as_poly_circuit(&self) -> &PolyCircuit {
        match self {
            BenchCircuit::AND(poly_circuit) => poly_circuit,
            BenchCircuit::AndVerify(poly_circuit) => poly_circuit,
        }
    }
}
