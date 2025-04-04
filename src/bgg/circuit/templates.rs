use super::PolyCircuit;
/// Returns a circuit that computes multiple inner products between
/// an input vector (`lhs`) and segments of a vector (`rhs`).
/// The `rhs` vector must be a multiple of the `lhs` vector in length.
/// Each segment of `rhs` is paired with `lhs` to compute an inner product,
/// and the circuit outputs all resulting inner products.
pub fn ip_circuit(lhs_len: usize, rhs_len: usize) -> PolyCircuit {
    assert!(rhs_len % lhs_len == 0, "rhs length must be a multiple of lhs length");
    let mut ip_circuit = PolyCircuit::new();
    let ip_inputs = ip_circuit.input(lhs_len + rhs_len);
    let lhs_inputs = &ip_inputs[0..lhs_len];
    let rhs_inputs = &ip_inputs[lhs_len..];
    let num_ip_outputs = rhs_len / lhs_len;

    let mut ip_outputs = Vec::with_capacity(num_ip_outputs);
    for out_idx in 0..num_ip_outputs {
        let mut ip_output = 0;
        for (priv_input, pub_output) in
            lhs_inputs.iter().zip(rhs_inputs[(lhs_len * out_idx)..(lhs_len * (out_idx + 1))].iter())
        {
            let mul = ip_circuit.mul_gate(*pub_output, *priv_input);
            if ip_output == 0 {
                ip_output = mul;
            } else {
                ip_output = ip_circuit.add_gate(ip_output, mul);
            }
        }
        ip_outputs.push(ip_output);
    }

    ip_circuit.output(ip_outputs);
    ip_circuit
}
