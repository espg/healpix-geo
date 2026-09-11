#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

use crate::maybe_parallelize;
use crate::scalar::morton::conversion as scalar;
use crate::vectorized::depth::Depth;

/// Vectorized [`scalar::to_nested`]; panics on a word that does not decode.
pub fn to_nested(ipix: &[u64], nthreads: usize) -> (Vec<u64>, Vec<u8>) {
    let mut result = Vec::<(u64, u8)>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, |hash| scalar::to_nested(hash)
        .unwrap_or_else(|| panic!("{} is not a valid morton cell id", hash)));

    result.into_iter().unzip()
}

pub fn from_nested(ipix: &[u64], depth: Depth, nthreads: usize) -> Vec<u64> {
    let mut result = Vec::<u64>::with_capacity(ipix.len());

    match depth {
        Depth::Scalar(d) => {
            maybe_parallelize!(nthreads, ipix, result, |hash| scalar::from_nested(hash, d));
        }
        Depth::Array(d) => {
            let zipped: Vec<_> = ipix.iter().zip(d.iter()).collect();
            maybe_parallelize!(nthreads, zipped, result, |(hash, depth)| {
                scalar::from_nested(hash, depth)
            });
        }
    };

    result
}

/// Vectorized [`scalar::to_ring`]; panics on a word that does not decode.
pub fn to_ring(ipix: &[u64], nthreads: usize) -> (Vec<u64>, Vec<u8>) {
    let mut result = Vec::<(u64, u8)>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, |hash| scalar::to_ring(hash)
        .unwrap_or_else(|| panic!("{} is not a valid morton cell id", hash)));

    result.into_iter().unzip()
}

pub fn from_ring(ipix: &[u64], depth: Depth, nthreads: usize) -> Vec<u64> {
    let mut result = Vec::<u64>::with_capacity(ipix.len());

    match depth {
        Depth::Scalar(d) => {
            maybe_parallelize!(nthreads, ipix, result, |hash| scalar::from_ring(hash, d));
        }
        Depth::Array(d) => {
            let zipped: Vec<_> = ipix.iter().zip(d.iter()).collect();
            maybe_parallelize!(nthreads, zipped, result, |(hash, depth)| {
                scalar::from_ring(hash, depth)
            });
        }
    };

    result
}

/// Vectorized [`scalar::to_zuniq`]; panics on a word that does not decode.
pub fn to_zuniq(ipix: &[u64], nthreads: usize) -> Vec<u64> {
    let mut result = Vec::<u64>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, |hash| scalar::to_zuniq(hash)
        .unwrap_or_else(|| panic!("{} is not a valid morton cell id", hash)));

    result
}

pub fn from_zuniq(ipix: &[u64], nthreads: usize) -> Vec<u64> {
    let mut result = Vec::<u64>::with_capacity(ipix.len());
    maybe_parallelize!(nthreads, ipix, result, scalar::from_zuniq);

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round_trip_matches_scalar() {
        let nested: Vec<u64> = vec![0, 164, 3071];
        let depths: Vec<u8> = vec![0, 3, 4];

        let words = from_nested(&nested, Depth::Array(&depths), 1);
        assert_eq!(
            words,
            vec![
                1152921504606846976,
                4107282860161892355,
                14983475960261640196
            ]
        );

        let (back, back_depths) = to_nested(&words, 1);
        assert_eq!(back, nested);
        assert_eq!(back_depths, depths);

        let zuniq = to_zuniq(&words, 1);
        assert_eq!(from_zuniq(&zuniq, 1), words);
    }

    #[test]
    fn test_scalar_depth() {
        let nested: Vec<u64> = vec![0, 164];
        let words = from_nested(&nested, Depth::Scalar(&3), 1);
        let (ring, depths) = to_ring(&words, 1);

        assert_eq!(from_ring(&ring, Depth::Array(&depths), 1), words);
    }
}
