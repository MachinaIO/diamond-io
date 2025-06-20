//! BGG+ partitioning

use crate::{bgg::BggEncoding, poly::PolyMatrix};
use rayon::prelude::*;

/// Immutable view returned by [`iter_pairs`].
#[derive(Clone, Copy)]
pub struct Pair<'a, M>
where
    M: PolyMatrix,
{
    pub c_one: &'a BggEncoding<M>,
    pub c_attr: &'a BggEncoding<M>,
}

impl<'a, M> Pair<'a, M>
where
    M: PolyMatrix + Clone,
{
    /// Return A_i: (d+1) x 2m
    pub fn a_matrix(&self) -> M {
        let a_one = &self.c_one.pubkey.matrix;
        let a_attr = &self.c_attr.pubkey.matrix;
        a_one.concat_columns(&[&a_attr])
    }

    /// Return c_i = c_one | c_attr : (1 x 2m)
    pub fn combined_matrix(&self) -> M {
        let v_one = &self.c_one.vector;
        let v_attr = &self.c_attr.vector;
        v_one.concat_columns(&[v_attr])
    }
}

pub fn iter_pairs<M>(encs: &[BggEncoding<M>]) -> impl IndexedParallelIterator<Item = Pair<'_, M>>
where
    M: PolyMatrix,
{
    assert!(encs.len() >= 2, "need at least one attribute");
    let (c_one, attrs) = encs.split_first().unwrap();
    attrs.par_iter().map(move |attr| Pair { c_one, c_attr: attr })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{bgg::lut::utils::random_bgg_encodings_for_bits, poly::dcrt::DCRTPolyParams};

    #[test]
    fn pair_count_matches_attrs() {
        let params = DCRTPolyParams::default();
        let encs = random_bgg_encodings_for_bits(8, 3, &params);
        let l = encs.len() - 1;

        let cnt = iter_pairs(&encs).count();
        assert_eq!(cnt, l);
    }
}
