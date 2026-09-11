use healpix_geo::scalar::morton::hierarchy;
use wasm_bindgen::prelude::*;

use crate::grid::{from_morton_checked, to_level};

/// Codec-only operations on `morton` cell ids.
///
/// A `morton` id is a packed-`u64` Morton word (the `mortie-core` codec): it
/// encodes its own refinement level, a raw unsigned sort is a Z-order
/// traversal with parents immediately before their children, and containment
/// testing is plain prefix truncation. None of these operations need a
/// reference body, so they live here rather than on `Grid`; coordinate math
/// on morton ids goes through `new Grid({ scheme: "morton", ... })`.
///
/// Every method validates its ids the way `Grid` does, so all of them take
/// canonical area words only: max-encoded point words (a coordinate cast to
/// level 29 with no area claim) are rejected here too and stay a codec-level
/// concept of the core crate.
#[wasm_bindgen(js_name = morton)]
pub struct Morton;

impl Morton {
    pub(crate) fn level_impl(cell: u64) -> Result<u8, String> {
        let (level, _) = from_morton_checked(cell)?;

        Ok(level)
    }

    pub(crate) fn parent_impl(cell: u64, level: f64) -> Result<u64, String> {
        from_morton_checked(cell)?;
        let level = to_level(level)?;

        // the id decodes (checked above), so `parent` cannot return `None`
        Ok(hierarchy::parent(&cell, &level).unwrap())
    }

    pub(crate) fn contains_impl(ancestor: u64, descendant: u64) -> Result<bool, String> {
        from_morton_checked(ancestor)?;
        from_morton_checked(descendant)?;

        Ok(hierarchy::contains(&ancestor, &descendant))
    }
}

#[wasm_bindgen(js_class = morton)]
impl Morton {
    /// The refinement level encoded in the given cell id
    #[wasm_bindgen(js_name = level)]
    pub fn level(cell: u64) -> Result<u8, JsValue> {
        Morton::level_impl(cell).map_err(|message| JsError::new(&message).into())
    }

    /// The ancestor of the given cell id at `level`
    ///
    /// Truncates the id to `level`, discarding finer detail; a `level` at or
    /// above the id's own level returns the id unchanged.
    #[wasm_bindgen(js_name = parent)]
    pub fn parent(cell: u64, level: f64) -> Result<u64, JsValue> {
        Morton::parent_impl(cell, level).map_err(|message| JsError::new(&message).into())
    }

    /// Whether `ancestor` contains `descendant` (containment-is-truncation)
    ///
    /// True iff truncating `descendant` to `ancestor`'s level yields exactly
    /// `ancestor`; a cell contains itself.
    #[wasm_bindgen(js_name = contains)]
    pub fn contains(ancestor: u64, descendant: u64) -> Result<bool, JsValue> {
        Morton::contains_impl(ancestor, descendant).map_err(|message| JsError::new(&message).into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use healpix_geo::scalar::morton::conversion::from_nested;

    #[test]
    fn test_level() {
        assert_eq!(Morton::level_impl(from_nested(&164, &3)).unwrap(), 3);
        assert_eq!(Morton::level_impl(from_nested(&2, &0)).unwrap(), 0);
        assert!(Morton::level_impl(0).is_err());
        assert!(Morton::level_impl(u64::MAX).is_err());
    }

    #[test]
    fn test_parent() {
        let cell = from_nested(&164, &3);

        assert_eq!(
            Morton::parent_impl(cell, 2.0).unwrap(),
            from_nested(&41, &2)
        );
        assert_eq!(Morton::parent_impl(cell, 3.0).unwrap(), cell);
        assert_eq!(Morton::parent_impl(cell, 29.0).unwrap(), cell);
        assert!(Morton::parent_impl(0, 2.0).is_err());
        assert!(Morton::parent_impl(cell, 2.5).is_err());
        assert!(Morton::parent_impl(cell, 30.0).is_err());
        assert!(Morton::parent_impl(cell, -1.0).is_err());
    }

    #[test]
    fn test_contains() {
        let child = from_nested(&164, &3);
        let parent = from_nested(&41, &2);

        assert!(Morton::contains_impl(parent, child).unwrap());
        assert!(!Morton::contains_impl(child, parent).unwrap());
        assert!(Morton::contains_impl(child, child).unwrap());
        assert!(Morton::contains_impl(0, child).is_err());
        assert!(Morton::contains_impl(child, u64::MAX).is_err());
    }

    #[test]
    fn test_rejects_point_words() {
        use healpix_geo::scalar::morton::conversion::from_nested_point;

        let nested = 164u64 << (2 * 26);
        let point = from_nested_point(&nested);
        let area = from_nested(&nested, &29);

        assert!(Morton::level_impl(point).is_err());
        assert!(Morton::parent_impl(point, 3.0).is_err());
        assert!(Morton::contains_impl(point, point).is_err());
        assert!(Morton::contains_impl(area, point).is_err());

        // the area cell of the same body is unaffected
        assert_eq!(Morton::level_impl(area).unwrap(), 29);
    }
}
