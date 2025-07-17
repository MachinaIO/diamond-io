use crate::{bgg::BggEncoding, poly::PolyMatrix};

pub struct CRTBggEncoding<M: PolyMatrix> {
    pub inner: Vec<BggEncoding<M>>,
}

impl<M: PolyMatrix> CRTBggEncoding<M> {
    // pub fn crt_decompose(
    //     params: <M::P as Poly>::Params,
    //     public_key: BggPublicKey<M>,
    //     encoding_q: BggEncoding<M>,
    // ) -> Self {
    //     let crt_params = params.to_crt();
    //     let mut inner = Vec::new();
    //     for modulus in crt_params {
    //         // q_i is param.modulus
    //         let new_vector = encoding_q.vector.modulus_switch(&modulus);
    //         let new_pubkey = encoding_q.pubkey.modulus_switch(&modulus);
    //         let new_plaintext =
    //             encoding_q.plaintext.clone().map(|poly| poly.modulus_switch(&params, modulus));
    //         let encoding = BggEncoding::new(new_vector, public_key.clone(), new_plaintext);
    //         inner.push(encoding);
    //     }
    //     Self { inner }
    // }

    // pub fn to_encoding(self) -> BggEncoding<M> {
    //     for encoding in self.inner {}
    // }
}
