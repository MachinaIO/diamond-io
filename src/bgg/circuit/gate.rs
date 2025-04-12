use crate::poly::Poly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyGate<P: Poly> {
    pub gate_id: usize,
    pub gate_type: PolyGateType<P>,
    pub input_gates: Vec<usize>,
    pub input_gates_scalar: Vec<usize>,
}

impl<P: Poly> PolyGate<P> {
    pub fn new(
        gate_id: usize,
        gate_type: PolyGateType<P>,
        input_gates: Vec<usize>,
        input_gates_scalar: Vec<usize>,
    ) -> Self {
        Self { gate_id, gate_type, input_gates, input_gates_scalar }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolyGateType<P: Poly> {
    Input,
    InputScalar,
    Const { bits: Vec<bool> },
    Add,
    Sub,
    Mul,
    ScalarMul,
    Rotate { shift: usize },
    Call { circuit_id: usize, num_input: usize, num_scalar_input: usize, output_id: usize },
    _Phantom(std::marker::PhantomData<P>),
}

impl<P: Poly> PolyGateType<P> {
    pub fn num_input(&self) -> usize {
        match self {
            PolyGateType::Input | PolyGateType::InputScalar | PolyGateType::Const { .. } => 0,
            PolyGateType::Rotate { .. } | PolyGateType::ScalarMul => 1,
            PolyGateType::Add | PolyGateType::Sub | PolyGateType::Mul => 2,
            PolyGateType::Call { num_input, .. } => *num_input,
            PolyGateType::_Phantom(_) => 0,
        }
    }

    pub fn num_scalar_input(&self) -> usize {
        match self {
            PolyGateType::Input | PolyGateType::InputScalar | PolyGateType::Const { .. } => 0,
            PolyGateType::Rotate { .. } => 0,
            PolyGateType::Add | PolyGateType::Sub | PolyGateType::Mul => 0,
            PolyGateType::ScalarMul => 1,
            PolyGateType::Call { num_scalar_input, .. } => *num_scalar_input,
            PolyGateType::_Phantom(_) => 0,
        }
    }
}
