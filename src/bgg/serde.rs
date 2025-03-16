use super::{BggEncoding, BggPublicKey};
use crate::poly::{Poly, PolyMatrix};
use bytes::Bytes;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SerializableBggPublicKey {
    pub matrix: Vec<Bytes>,
    pub reveal_plaintext: bool,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn to_compact_bytes(&self) -> SerializableBggPublicKey {
        SerializableBggPublicKey {
            matrix: self.matrix.to_compact_bytes(),
            reveal_plaintext: self.reveal_plaintext,
        }
    }
}

impl SerializableBggPublicKey {
    pub fn from_compact_bytes<M: PolyMatrix>(
        self,
        params: &<M::P as Poly>::Params,
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
    pub fn to_compact_bytes(&self) -> SerializableBggEncoding {
        SerializableBggEncoding {
            vector: self.vector.to_compact_bytes(),
            pubkey: self.pubkey.to_compact_bytes(),
            // todo: we don't know yet how to turn poly into bytes
            plaintext: None,
        }
    }
}

impl SerializableBggEncoding {
    pub fn from_compact_bytes<M: PolyMatrix>(
        self,
        params: &<M::P as Poly>::Params,
    ) -> BggEncoding<M> {
        BggEncoding {
            vector: M::from_compact_bytes(params, self.vector),
            pubkey: self.pubkey.from_compact_bytes(params),
            // todo: we don't know yet how to turn poly into bytes
            plaintext: None,
        }
    }
}
