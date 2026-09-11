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
/// Returns the word unchanged when `depth >=` its own depth, and `None` if
/// the word does not decode.
pub fn parent(hash: &u64, depth: &u8) -> Option<u64> {
    mortie::coarsen(*hash, *depth)
}

/// Containment-is-truncation: whether `ancestor` contains `descendant`.
///
/// True iff truncating `descendant` to `ancestor`'s depth yields exactly
/// `ancestor` — no coordinate math involved. Reflexive for canonical area
/// words (`contains(w, w)` is true); false whenever either word does not
/// decode or `ancestor` is not canonical.
pub fn contains(ancestor: &u64, descendant: &u64) -> bool {
    depth(ancestor) <= depth(descendant)
        && mortie::coarsen(*descendant, depth(ancestor)) == Some(*ancestor)
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
    fn test_parent_golden() {
        // generated with the published mortie-core 0.1.0
        let word = from_nested(&164, &3);
        assert_eq!(parent(&word, &1), Some(4035225266123964417));
        assert_eq!(parent(&word, &2), Some(4107282860161892354));
        assert_eq!(parent(&word, &3), Some(word));
        assert_eq!(parent(&word, &29), Some(word));
        assert_eq!(parent(&0, &1), None);
    }

    #[test]
    fn test_parent_matches_nested_truncation() {
        // nested 164 at depth 3 sits under nested 164 >> 2 at depth 2
        let word = from_nested(&164, &3);
        assert_eq!(parent(&word, &2), Some(from_nested(&(164 >> 2), &2)));
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
