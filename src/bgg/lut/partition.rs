//! split a (L+1)-wide BGG⁺ row into L attribute pairs.

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

pub fn iter_pairs<M>(
    encs: &[BggEncoding<M>],
) -> impl IndexedParallelIterator<Item = (&BggEncoding<M>, &BggEncoding<M>)>
where
    M: PolyMatrix,
{
    assert!(encs.len() >= 2, "need at least one attribute");
    let (c_one, attrs) = encs.split_first().unwrap();
    attrs.par_iter().map(move |attr| (c_one, attr))
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
