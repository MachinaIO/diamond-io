use crate::poly::Poly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyGate<P: Poly> {
    pub gate_id: usize,
    pub gate_type: PolyGateType<P>,
    pub input_gates: Vec<usize>,
}

impl<P: Poly> PolyGate<P> {
    pub fn new(gate_id: usize, gate_type: PolyGateType<P>, input_gates: Vec<usize>) -> Self {
        Self { gate_id, gate_type, input_gates }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolyGateType<P: Poly> {
    Input,
    Const { bits: Vec<bool> },
    Add,
    Sub,
    Mul,
    ScalarMul { scalar: P },
    Rotate { shift: usize },
    Call { circuit_id: usize, num_input: usize, output_id: usize },
}

impl<P: Poly> PolyGateType<P> {
    pub fn num_input(&self) -> usize {
        match self {
            PolyGateType::Input | PolyGateType::Const { .. } => 0,
            PolyGateType::Rotate { .. } | PolyGateType::ScalarMul { .. } => 1,
            PolyGateType::Add | PolyGateType::Sub | PolyGateType::Mul => 2,
            PolyGateType::Call { num_input, .. } => *num_input,
        }
    }
}
