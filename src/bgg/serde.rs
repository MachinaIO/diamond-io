use bytes::Bytes;
use serde::{Deserialize, Serialize};

use crate::poly::{Poly, PolyMatrix};

use super::{BggEncoding, BggPublicKey};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SerializableBggPublicKey {
    pub matrix: Vec<Bytes>,
    pub reveal_plaintext: bool,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn to_compact_bytes(&self, byte_size: usize, offset: usize) -> SerializableBggPublicKey {
        SerializableBggPublicKey {
            matrix: self.matrix.to_compact_bytes(byte_size, offset),
            reveal_plaintext: self.reveal_plaintext,
        }
    }
}

impl SerializableBggPublicKey {
    pub fn from_compact_bytes<M: PolyMatrix>(
        self,
        params: &<M::P as Poly>::Params,
        byte_size: usize,
    ) -> BggPublicKey<M> {
        BggPublicKey {
            matrix: M::from_compact_bytes(params, self.matrix),
            reveal_plaintext: self.reveal_plaintext,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SerializableBggEncoding {
    pub vector: Vec<Bytes>,
    pub pubkey: SerializableBggPublicKey,
    // todo: we don't know yet how to turn poly into bytes
    pub plaintext: Option<Bytes>,
}

impl<M: PolyMatrix> BggEncoding<M> {
    pub fn to_compact_bytes(&self, byte_size: usize) -> SerializableBggEncoding {
        SerializableBggEncoding {
            vector: self.vector.to_compact_bytes(byte_size, 0),
            pubkey: self.pubkey.to_compact_bytes(byte_size, 0),
            // todo: we don't know yet how to turn poly into bytes
            plaintext: None,
        }
    }
}

impl SerializableBggEncoding {
    pub fn from_compact_bytes<M: PolyMatrix>(
        self,
        params: &<M::P as Poly>::Params,
        byte_size: usize,
    ) -> BggEncoding<M> {
        BggEncoding {
            vector: M::from_compact_bytes(params, self.vector),
            pubkey: self.pubkey.from_compact_bytes(params, byte_size),
            // todo: we don't know yet how to turn poly into bytes
            plaintext: None,
        }
    }
}
