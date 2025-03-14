use bytes::Bytes;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    bgg::serde::SerializableBggEncoding,
    poly::{Poly, PolyMatrix},
};

use super::Obfuscation;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SerializableObfuscation {
    pub hash_key: [u8; 32],
    pub encodings_init: Vec<SerializableBggEncoding>,
    pub p_init: Vec<Bytes>,
    pub m_preimages: Vec<(Vec<Bytes>, Vec<Bytes>)>,
    pub n_preimages: Vec<(Vec<Bytes>, Vec<Bytes>)>,
    pub k_preimages: Vec<(Vec<Bytes>, Vec<Bytes>)>,
    pub final_preimage: Vec<Bytes>,
    #[cfg(test)]
    pub s_init: Bytes,
    #[cfg(test)]
    pub t_bar: Bytes,
    #[cfg(test)]
    pub bs: Vec<(Bytes, Bytes, Bytes)>,
    #[cfg(test)]
    pub hardcoded_key: Bytes,
    #[cfg(test)]
    pub enc_hardcoded_key: Bytes,
    #[cfg(test)]
    pub final_preimage_target: Bytes,
}

impl<M: PolyMatrix> Obfuscation<M> {
    pub fn to_compact_bytes(&self, byte_size: usize) -> SerializableObfuscation {
        SerializableObfuscation {
            hash_key: self.hash_key,
            encodings_init: self
                .encodings_init
                .iter()
                .map(|e| e.to_compact_bytes(byte_size))
                .collect_vec(),
            p_init: self.p_init.to_compact_bytes(byte_size),
            m_preimages: self
                .m_preimages
                .iter()
                .map(|(m_0, m_1)| {
                    (m_0.to_compact_bytes(byte_size), m_1.to_compact_bytes(byte_size))
                })
                .collect_vec(),
            n_preimages: self
                .n_preimages
                .iter()
                .map(|(n_0, n_1)| {
                    (n_0.to_compact_bytes(byte_size), n_1.to_compact_bytes(byte_size))
                })
                .collect_vec(),
            k_preimages: self
                .k_preimages
                .iter()
                .map(|(k_0, k_1)| {
                    (k_0.to_compact_bytes(byte_size), k_1.to_compact_bytes(byte_size))
                })
                .collect_vec(),
            final_preimage: self.final_preimage.to_compact_bytes(byte_size),
            #[cfg(test)]
            s_init: todo!(),
            #[cfg(test)]
            t_bar: todo!(),
            #[cfg(test)]
            bs: todo!(),
            #[cfg(test)]
            hardcoded_key: todo!(),
            #[cfg(test)]
            enc_hardcoded_key: todo!(),
            #[cfg(test)]
            final_preimage_target: todo!(),
        }
    }
}

impl SerializableObfuscation {
    pub fn from_compact_bytes<M: PolyMatrix>(
        self,
        params: &<M::P as Poly>::Params,
        byte_size: usize,
    ) -> Obfuscation<M> {
        Obfuscation {
            hash_key: self.hash_key,
            encodings_init: self
                .encodings_init
                .into_iter()
                .map(|e| e.from_compact_bytes(params, byte_size))
                .collect(),
            p_init: M::from_compact_bytes(params, self.p_init),
            m_preimages: self
                .m_preimages
                .into_iter()
                .map(|(m_0, m_1)| {
                    (M::from_compact_bytes(params, m_0), M::from_compact_bytes(params, m_1))
                })
                .collect(),
            n_preimages: self
                .n_preimages
                .into_iter()
                .map(|(n_0, n_1)| {
                    (M::from_compact_bytes(params, n_0), M::from_compact_bytes(params, n_1))
                })
                .collect(),
            k_preimages: self
                .k_preimages
                .into_iter()
                .map(|(k_0, k_1)| {
                    (M::from_compact_bytes(params, k_0), M::from_compact_bytes(params, k_1))
                })
                .collect(),
            final_preimage: M::from_compact_bytes(params, self.final_preimage),
            #[cfg(test)]
            s_init: todo!(),
            #[cfg(test)]
            t_bar: todo!(),
            #[cfg(test)]
            bs: todo!(),
            #[cfg(test)]
            hardcoded_key: todo!(),
            #[cfg(test)]
            enc_hardcoded_key: todo!(),
            #[cfg(test)]
            final_preimage_target: todo!(),
        }
    }
}
