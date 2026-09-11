use mortie_core::decimal_morton as mortie;

/// The refinement level encoded in a Morton word's suffix.
///
/// Total over every bit pattern (the suffix numbering has no invalid values);
/// combine with [`super::conversion::is_canonical`] when the word is
/// untrusted.
pub fn depth(hash: &u64) -> u8 {
    mortie::order_of(*hash)
}

/// The ancestor of a Morton word at `depth`, discarding finer detail.
///
/// `depth` is the **absolute** depth to land at, not a number of levels to
/// drop — the opposite reading of the `delta_depth` taken by
/// [`crate::vectorized::nested::hierarchy::parents`], which is why this is
/// not called `parent`. Returns the word unchanged when `depth >=` its own
/// depth, and `None` if the word does not decode.
pub fn ancestor(hash: &u64, depth: &u8) -> Option<u64> {
    mortie::coarsen(*hash, *depth)
}

/// Containment-is-truncation: whether `ancestor` contains `descendant`.
///
/// Decided on the decoded cell paths, so that it answers the question a tiled
/// store asks — does this word fall inside that tile — rather than a question
/// about bit patterns:
///
/// - `ancestor` must be a canonical **area** word: a point word claims no
///   area, so it contains nothing, not even itself, and a non-canonical word
///   is not an id at all.
/// - `descendant` may be any canonical word, area or point.
/// - true iff `ancestor` is no deeper than `descendant` and truncating
///   `descendant`'s nested path to `ancestor`'s depth yields `ancestor`'s
///   cell — no coordinate math involved.
///
/// So it is reflexive exactly on canonical area words, and a depth-29 area
/// cell contains the point words falling inside it — which comparing the
/// truncated *words* instead of the cells they name would miss.
pub fn contains(ancestor: &u64, descendant: &u64) -> bool {
    use super::conversion::{is_canonical, is_point, to_nested};

    if is_point(ancestor) || !is_canonical(ancestor) || !is_canonical(descendant) {
        return false;
    }

    // canonicality implies both words decode
    let (nested_a, depth_a) = to_nested(ancestor).unwrap();
    let (nested_b, depth_b) = to_nested(descendant).unwrap();

    depth_a <= depth_b && nested_b >> (2 * u32::from(depth_b - depth_a)) == nested_a
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::morton::conversion::from_nested;

    #[test]
    fn test_depth() {
        assert_eq!(depth(&from_nested(&0, &0)), 0);
        assert_eq!(depth(&from_nested(&164, &3)), 3);
        assert_eq!(depth(&from_nested(&0, &29)), 29);
    }

    #[test]
    fn test_ancestor_golden() {
        // generated with the published mortie-core 0.1.0
        let word = from_nested(&164, &3);
        assert_eq!(ancestor(&word, &1), Some(4035225266123964417));
        assert_eq!(ancestor(&word, &2), Some(4107282860161892354));
        assert_eq!(ancestor(&word, &3), Some(word));
        assert_eq!(ancestor(&word, &29), Some(word));
        assert_eq!(ancestor(&0, &1), None);
    }

    #[test]
    fn test_ancestor_matches_nested_truncation() {
        // nested 164 at depth 3 sits under nested 164 >> 2 at depth 2
        let word = from_nested(&164, &3);
        assert_eq!(ancestor(&word, &2), Some(from_nested(&(164 >> 2), &2)));
    }

    #[test]
    fn test_contains() {
        let child = from_nested(&164, &3);
        let parent_ = from_nested(&(164 >> 2), &2);
        let other = from_nested(&((164 >> 2) + 1), &2);

        assert!(contains(&parent_, &child));
        assert!(!contains(&other, &child));
        assert!(!contains(&child, &parent_)); // not symmetric
        assert!(contains(&child, &child)); // reflexive

        // invalid words never contain or get contained
        assert!(!contains(&0, &child));
        assert!(!contains(&child, &0));
    }

    #[test]
    fn test_contains_requires_a_canonical_ancestor() {
        let word = from_nested(&164, &3);
        let junk = word | (1 << 30); // decodes to the same cell, not canonical

        assert!(!contains(&junk, &junk));
        assert!(!contains(&junk, &word));
        assert!(!contains(&word, &junk));
        assert!(!contains(&0, &0));
        assert!(!contains(&u64::MAX, &u64::MAX));
    }

    #[test]
    fn test_contains_covers_point_words() {
        use crate::scalar::morton::conversion::from_nested_point;

        let nested = 164u64 << (2 * 26);
        let point = from_nested_point(&nested);
        let area = from_nested(&nested, &29);

        // the depth-29 area cell covering a point contains it, as does every
        // proper ancestor
        assert!(contains(&area, &point));
        assert!(contains(&from_nested(&164, &3), &point));
        assert!(contains(&from_nested(&(164 >> 6), &0), &point));
        assert!(!contains(&from_nested(&165, &3), &point));

        // a point claims no area, so it contains nothing — not even itself
        assert!(!contains(&point, &point));
        assert!(!contains(&point, &area));
    }

    #[test]
    fn test_base_cell_contains_everything_under_it() {
        let base = from_nested(&2, &0);
        for depth_ in [1u8, 5, 29] {
            let first = from_nested(&(2u64 << (2 * depth_ as u32)), &depth_);
            assert!(contains(&base, &first));
        }
        let elsewhere = from_nested(&(3u64 << 4), &2);
        assert!(!contains(&base, &elsewhere));
    }
}
