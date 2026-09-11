#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

use crate::maybe_parallelize;
use crate::scalar::morton::hierarchy as scalar;

/// Vectorized [`scalar::depth`].
pub fn depth(ipix: &[u64], nthreads: usize) -> Vec<u8> {
    let mut result = Vec::<u8>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, scalar::depth);

    result
}

/// Vectorized [`scalar::parent`]; panics on a word that does not decode.
pub fn parent(ipix: &[u64], depth: &u8, nthreads: usize) -> Vec<u64> {
    let mut result = Vec::<u64>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, |hash| scalar::parent(hash, depth)
        .unwrap_or_else(|| panic!("{} is not a valid morton cell id", hash)));

    result
}

/// Whether `ancestor` contains each of `ipix` (containment-is-truncation).
pub fn contains(ancestor: &u64, ipix: &[u64], nthreads: usize) -> Vec<bool> {
    let mut result = Vec::<bool>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, |hash| scalar::contains(
        ancestor, hash
    ));

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::morton::conversion::from_nested;

    #[test]
    fn test_parent_and_contains() {
        let words: Vec<u64> = vec![from_nested(&164, &3), from_nested(&165, &3)];

        assert_eq!(depth(&words, 1), vec![3, 3]);

        let parents = parent(&words, &2, 1);
        assert_eq!(parents, vec![from_nested(&41, &2), from_nested(&41, &2)]);

        let ancestor = from_nested(&41, &2);
        assert_eq!(contains(&ancestor, &words, 1), vec![true, true]);

        let elsewhere = from_nested(&40, &2);
        assert_eq!(contains(&elsewhere, &words, 1), vec![false, false]);
    }
}
