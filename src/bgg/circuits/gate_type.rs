use crate::poly::{Poly, PolyMatrix};

#[derive(Debug, Clone)]
pub enum PolyGateType<M: PolyMatrix> {
    Input,
    Output,
    Add,
    Sub,
    ScalarMul(M::P),
    MatrixMul(M),
    Mul,
}

impl<M: PolyMatrix> PolyGateType<M> {
    pub fn num_input(&self) -> usize {
        match self {
            PolyGateType::Input => 1,
            PolyGateType::Output => 0,
            PolyGateType::Add => 2,
            PolyGateType::Sub => 2,
            PolyGateType::ScalarMul(_) => 1,
            PolyGateType::MatrixMul(m) => m.col_size(),
            PolyGateType::Mul => 2,
        }
    }

    pub fn num_output(&self) -> usize {
        match self {
            PolyGateType::Input => 0,
            PolyGateType::Output => 1,
            PolyGateType::Add => 1,
            PolyGateType::Sub => 1,
            PolyGateType::ScalarMul(_) => 1,
            PolyGateType::Mul => 1,
        }
    }
}
