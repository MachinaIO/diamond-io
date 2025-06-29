#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyGate {
    pub gate_id: usize,
    pub gate_type: PolyGateType,
    pub input_gates: Vec<usize>,
}

impl PolyGate {
    pub fn new(gate_id: usize, gate_type: PolyGateType, input_gates: Vec<usize>) -> Self {
        Self { gate_id, gate_type, input_gates }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum PolyGateType {
    Input,
    Const { digits: Vec<u32> },
    Add,
    Sub,
    Mul,
    Rotate { shift: usize },
    Call { circuit_id: usize, num_input: usize, output_id: usize },
    PubLut { lookup_id: usize },
}

impl PolyGateType {
    pub fn num_input(&self) -> usize {
        match self {
            PolyGateType::Input | PolyGateType::Const { .. } => 0,
            PolyGateType::Rotate { .. } | PolyGateType::PubLut { .. } => 1,
            PolyGateType::Add | PolyGateType::Sub | PolyGateType::Mul => 2,
            PolyGateType::Call { num_input, .. } => *num_input,
        }
    }
}
