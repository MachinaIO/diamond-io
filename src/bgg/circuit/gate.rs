use crate::poly::Poly;

#[derive(Debug, Clone)]
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
    Add,
    Sub,
    ScalarMul(P),
    Mul,
    Call { circuit_id: usize, num_input: usize, output_id: usize },
}

impl<P: Poly> PolyGateType<P> {
    pub fn num_input(&self) -> usize {
        match self {
            PolyGateType::Input => 0,
            PolyGateType::ScalarMul(_) => 1,
            PolyGateType::Add | PolyGateType::Sub | PolyGateType::Mul => 2,
            PolyGateType::Call { num_input, .. } => *num_input,
        }
    }
}
