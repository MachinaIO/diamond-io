use serde::{Deserialize, Serialize};
use std::fmt::Display;
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
)]
pub struct GateId(pub(crate) usize);

impl Display for GateId {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyGate {
    pub gate_id: GateId,
    pub gate_type: PolyGateType,
    pub input_gates: Vec<GateId>,
}

impl PolyGate {
    pub fn new(gate_id: GateId, gate_type: PolyGateType, input_gates: Vec<GateId>) -> Self {
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
