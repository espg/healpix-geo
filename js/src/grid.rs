use cdshealpix as healpix;
use geodesy::prelude::EllipsoidBase;
use healpix_geo_core::ellipsoid::{Ellipsoid as RustEllipsoid, ReferenceBody};
use healpix_geo_core::scalar::nested::coordinates as scalar;
use serde::Deserialize;
use serde_wasm_bindgen::from_value;
use wasm_bindgen::prelude::*;

use crate::coordinates::Coordinate;
use crate::ellipsoid::EllipsoidLike;
use crate::geometry::spherical_vertex;

const MAX_DEPTH: u8 = 29;

/// Number of cells of a full HEALPix grid at `depth` (12 · 4^depth).
fn n_cells(depth: u8) -> u64 {
    12u64 << (2 * depth)
}

/// Decode a `zuniq` cell id.
///
/// `cdshealpix`'s `from_zuniq` is infallible in the type system but panics —
/// an unrecoverable wasm trap — for values that are not valid zuniq
/// encodings: `0` underflows the depth computation, and ids whose sentinel bit
/// sits at an odd position or whose hash is out of range at the encoded depth
/// (`u64::MAX`, for instance) blow up further down. Validate first.
fn from_zuniq_checked(cell: u64) -> Result<(u8, u64), String> {
    // a well-formed zuniq id has its sentinel bit at 2·(29 − depth); note
    // that `0` has 64 trailing zeros and is rejected by the upper bound
    let trailing = cell.trailing_zeros();
    if trailing % 2 != 0 || trailing > 2 * u32::from(MAX_DEPTH) {
        return Err(format!("{} is not a valid zuniq cell id", cell));
    }

    let (depth, hash) = healpix::nested::from_zuniq(cell);
    if hash >= n_cells(depth) {
        return Err(format!(
            "{} is not a valid zuniq cell id: hash {} is out of range at depth {}",
            cell, hash, depth
        ));
    }

    Ok((depth, hash))
}

#[derive(Deserialize, Debug, Clone, Copy, PartialEq)]
#[serde(rename_all = "lowercase")]
pub(crate) enum Scheme {
    Nested,
    Ring,
    Zuniq,
}

impl Scheme {
    fn name(&self) -> &'static str {
        match self {
            Self::Nested => "nested",
            Self::Ring => "ring",
            Self::Zuniq => "zuniq",
        }
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct GridOptions {
    scheme: Scheme,
    #[serde(default)]
    depth: Option<u8>,
    #[serde(default)]
    nside: Option<u64>,
    #[serde(default)]
    ellipsoid: Option<EllipsoidLike>,
}

/// The options object accepted by `createGrid`.
#[wasm_bindgen(typescript_custom_section)]
const GRID_OPTIONS: &'static str = r#"
export type GridOptions = {
    scheme: "nested" | "ring" | "zuniq";
    depth?: number;
    nside?: number;
    ellipsoid?: EllipsoidInput | null;
};
"#;

/// A HEALPix grid at a fixed depth on a fixed reference body.
///
/// Constructed through [`create_grid`] (`createGrid` in JS). The scheme
/// dispatch, the depth, and the parsed ellipsoid state (including the
/// authalic-latitude coefficients) are stored once, so per-call overhead is
/// limited to the actual coordinate math.
///
/// Every method validates its inputs and reports misuse as a catchable JS
/// `Error` rather than trapping the wasm instance.
#[wasm_bindgen]
pub struct Grid {
    scheme: Scheme,
    depth: u8,
    pub(crate) ellipsoid: RustEllipsoid,
}

impl Grid {
    pub(crate) fn from_options(options: GridOptions) -> Result<Grid, String> {
        let depth = match (options.depth, options.nside) {
            (Some(_), Some(_)) => {
                return Err("pass either `depth` or `nside`, not both".to_string());
            }
            (Some(depth), None) => depth,
            (None, Some(nside)) => {
                if !nside.is_power_of_two() {
                    return Err(format!("`nside` must be a power of two, got {}", nside));
                }
                nside.trailing_zeros() as u8
            }
            (None, None) => {
                return Err("one of `depth` or `nside` is required".to_string());
            }
        };
        if depth > MAX_DEPTH {
            return Err(format!(
                "`depth` must be at most {}, got {}",
                MAX_DEPTH, depth
            ));
        }

        let ellipsoid = options
            .ellipsoid
            .map(|e| e.into_ellipsoid())
            .unwrap_or_default();

        Ok(Grid {
            scheme: options.scheme,
            depth,
            ellipsoid,
        })
    }

    /// Number of cells at the grid's depth.
    fn n_cells(&self) -> u64 {
        n_cells(self.depth)
    }

    fn check_cell(&self, cell: u64) -> Result<(), String> {
        if cell >= self.n_cells() {
            return Err(format!(
                "cell id {} is out of range at depth {} ({} cells)",
                cell,
                self.depth,
                self.n_cells()
            ));
        }

        Ok(())
    }

    /// The cell in the nested scheme plus the depth it lives at.
    ///
    /// The `ring` branch converts to nested and then uses the nested
    /// coordinate implementation (the same route `js/src/ring.rs` takes),
    /// rather than `healpix_geo_core::scalar::ring::coordinates`; the two
    /// agree everywhere tested, but only one of them is exercised.
    pub(crate) fn to_nested(&self, cell: u64) -> Result<(u8, u64), String> {
        match self.scheme {
            Scheme::Nested => {
                self.check_cell(cell)?;
                Ok((self.depth, cell))
            }
            Scheme::Ring => {
                self.check_cell(cell)?;
                Ok((self.depth, healpix::nested::get(self.depth).from_ring(cell)))
            }
            Scheme::Zuniq => from_zuniq_checked(cell),
        }
    }

    fn from_nested(&self, hash: u64) -> u64 {
        match self.scheme {
            Scheme::Nested => hash,
            Scheme::Ring => healpix::nested::get(self.depth).to_ring(hash),
            Scheme::Zuniq => healpix::nested::to_zuniq(self.depth, hash),
        }
    }

    pub(crate) fn center_impl(&self, cell: u64) -> Result<Coordinate, String> {
        let (depth, hash) = self.to_nested(cell)?;
        let layer = healpix::nested::get(depth);

        let (lon, lat) = scalar::healpix_to_lonlat(&hash, layer, &self.ellipsoid);

        Ok(Coordinate { lon, lat })
    }

    pub(crate) fn vertex_impl(&self, cell: u64, u: f64, v: f64) -> Result<Coordinate, String> {
        let (depth, hash) = self.to_nested(cell)?;
        let layer = healpix::nested::get(depth);

        let center = layer.center_of_projected_cell(hash);
        let (lon, lat) = spherical_vertex(center, depth, (u, v));

        Ok(Coordinate {
            lon: lon.to_degrees().rem_euclid(360.0),
            lat: self
                .ellipsoid
                .latitude_authalic_to_geographic(lat)
                .to_degrees(),
        })
    }

    pub(crate) fn bit_combine_impl(&self, i: u32, j: u32) -> Result<u64, String> {
        let nside = self.nside();
        if i >= nside || j >= nside {
            return Err(format!(
                "z-order coordinates ({}, {}) are out of range for nside {}",
                i, j, nside
            ));
        }

        let zoc = healpix::nested::zordercurve::get_zoc(self.depth);

        Ok(self.from_nested(zoc.ij2h(i, j)))
    }
}

#[wasm_bindgen]
impl Grid {
    /// The indexing scheme: "nested", "ring" or "zuniq"
    #[wasm_bindgen(getter)]
    pub fn scheme(&self) -> String {
        self.scheme.name().to_string()
    }

    /// The depth (refinement level) of the grid
    #[wasm_bindgen(getter)]
    pub fn depth(&self) -> u8 {
        self.depth
    }

    /// The nside (2^depth) of the grid
    #[wasm_bindgen(getter)]
    pub fn nside(&self) -> u32 {
        1u32 << self.depth
    }

    /// Semi-major axis of the reference body (the radius, for a sphere), in
    /// meters
    #[wasm_bindgen(getter, js_name = semiMajorAxis)]
    pub fn semi_major_axis(&self) -> f64 {
        self.ellipsoid.ellipsoid().semimajor_axis()
    }

    /// Flattening of the reference body (0 for a sphere)
    #[wasm_bindgen(getter)]
    pub fn flattening(&self) -> f64 {
        self.ellipsoid.ellipsoid().flattening()
    }

    /// Whether the reference body is a sphere
    #[wasm_bindgen(getter, js_name = isSphere)]
    pub fn is_sphere(&self) -> bool {
        matches!(self.ellipsoid, RustEllipsoid::Sphere(_))
    }

    /// The reference body of the grid, as a **new** handle
    ///
    /// Every read clones the parsed ellipsoid into a freshly allocated
    /// `Ellipsoid` (wasm allocation + JS wrapper + finalizer registration), so
    /// this is not a free property access: hoist it out of hot paths and
    /// `free()` it when done, or use the `semiMajorAxis` / `flattening` /
    /// `isSphere` getters above, which return plain numbers.
    #[wasm_bindgen(getter)]
    pub fn ellipsoid(&self) -> crate::ellipsoid::Ellipsoid {
        crate::ellipsoid::Ellipsoid {
            inner: self.ellipsoid.clone(),
        }
    }

    /// Center coordinates of the given cell
    pub fn center(&self, cell: u64) -> Result<Coordinate, JsValue> {
        self.center_impl(cell)
            .map_err(|message| JsError::new(&message).into())
    }

    /// The cell containing the given coordinate
    #[wasm_bindgen(js_name = cellAt)]
    pub fn cell_at(&self, lon: f64, lat: f64) -> u64 {
        let layer = healpix::nested::get(self.depth);
        let hash = scalar::lonlat_to_healpix(&lon, &lat, layer, &self.ellipsoid);

        self.from_nested(hash)
    }

    /// Single vertex of the given cell
    ///
    /// `u` and `v` are offsets from the southern vertex of the cell. For the
    /// `zuniq` scheme, the depth encoded in the cell id is used.
    pub fn vertex(&self, cell: u64, u: f64, v: f64) -> Result<Coordinate, JsValue> {
        self.vertex_impl(cell, u, v)
            .map_err(|message| JsError::new(&message).into())
    }

    /// Cell index at the given z-order coordinates
    ///
    /// Interleaves the bits of `i` and `j` (the two axes of the nested
    /// z-order numbering within a base-resolution pixel) into the cell index
    /// at the grid's depth, expressed in the grid's scheme. Both coordinates
    /// have to be smaller than `nside`.
    #[wasm_bindgen(js_name = bitCombine)]
    pub fn bit_combine(&self, i: u32, j: u32) -> Result<u64, JsValue> {
        self.bit_combine_impl(i, j)
            .map_err(|message| JsError::new(&message).into())
    }
}

/// Create a [`Grid`] handle from plain options.
///
/// Options:
/// - `scheme`: `"nested"`, `"ring"` or `"zuniq"` (required)
/// - `depth` or `nside`: the refinement level (exactly one required)
/// - `ellipsoid`: a plain object as accepted by `Ellipsoid.from`, or
///   `null`/absent for the default sphere
#[wasm_bindgen(js_name = createGrid)]
pub fn create_grid(
    #[wasm_bindgen(unchecked_param_type = "GridOptions")] options: JsValue,
) -> Result<Grid, JsValue> {
    let options: GridOptions = from_value(options)?;

    Grid::from_options(options).map_err(|message| JsError::new(&message).into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options(scheme: Scheme, depth: Option<u8>, nside: Option<u64>) -> GridOptions {
        GridOptions {
            scheme,
            depth,
            nside,
            ellipsoid: None,
        }
    }

    fn grid(scheme: Scheme, depth: u8) -> Grid {
        Grid::from_options(options(scheme, Some(depth), None)).unwrap()
    }

    #[test]
    fn test_depth_nside_validation() {
        assert!(Grid::from_options(options(Scheme::Nested, Some(2), Some(4))).is_err());
        assert!(Grid::from_options(options(Scheme::Nested, None, None)).is_err());
        assert!(Grid::from_options(options(Scheme::Nested, Some(30), None)).is_err());
        assert!(Grid::from_options(options(Scheme::Nested, None, Some(3))).is_err());

        let from_nside = Grid::from_options(options(Scheme::Nested, None, Some(8))).unwrap();
        assert_eq!(from_nside.depth(), 3);
        assert_eq!(from_nside.nside(), 8);
    }

    #[test]
    fn test_vertex_matches_scheme_statics() {
        let grid = grid(Scheme::Nested, 0);

        let expected: Vec<(f64, f64)> = vec![
            (45.0, 0.0),
            (67.5, 19.47122063),
            (90.0, 41.8103149),
            (45.0, 90.0),
        ];
        let uv: Vec<(f64, f64)> = vec![(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.0, 1.0)];

        for ((u, v), (lon, lat)) in uv.into_iter().zip(expected) {
            let actual = grid.vertex_impl(0, u, v).unwrap();
            assert!((actual.lon - lon).abs() < 1e-4);
            assert!((actual.lat - lat).abs() < 1e-4);
        }
    }

    #[test]
    fn test_zuniq_uses_encoded_depth() {
        // depth 0, nested cell 0 encoded as zuniq
        let cell = healpix::nested::to_zuniq(0, 0);
        let grid = grid(Scheme::Zuniq, 4);

        let vertex = grid.vertex_impl(cell, 0.0, 0.0).unwrap();
        assert!((vertex.lon - 45.0).abs() < 1e-4);
        assert!(vertex.lat.abs() < 1e-4);
    }

    #[test]
    fn test_ring_roundtrip() {
        let grid = grid(Scheme::Ring, 4);

        let center = grid.center_impl(164).unwrap();
        let cell = grid.cell_at(center.lon, center.lat);
        assert_eq!(cell, 164);
    }

    #[test]
    fn test_flattened_ellipsoid_getters_match_the_handle() {
        let grid = grid(Scheme::Nested, 4);

        assert_eq!(grid.semi_major_axis(), grid.ellipsoid().semi_major_axis());
        assert_eq!(grid.flattening(), grid.ellipsoid().flattening());
        assert_eq!(grid.is_sphere(), grid.ellipsoid().is_sphere());
        assert!(grid.is_sphere());
    }

    #[test]
    fn test_rejects_out_of_range_cells() {
        // 12 · 4^4 == 3072 is the first invalid id at depth 4; cdshealpix
        // panics on it ("Wrong hash value: too large")
        for scheme in [Scheme::Nested, Scheme::Ring] {
            let grid = grid(scheme, 4);
            assert_eq!(grid.n_cells(), 3072);

            assert!(grid.center_impl(3072).is_err());
            assert!(grid.vertex_impl(3072, 0.5, 0.5).is_err());
            assert!(grid.center_impl(u64::MAX).is_err());
            assert!(grid.center_impl(3071).is_ok());
        }
    }

    #[test]
    fn test_rejects_invalid_zuniq_cells() {
        let grid = grid(Scheme::Zuniq, 4);

        // `0` used to underflow the depth computation
        assert!(grid.center_impl(0).is_err());
        // sentinel bit at an odd position
        assert!(grid.center_impl(2).is_err());
        // hash out of range at the encoded depth
        assert!(grid.center_impl(u64::MAX).is_err());

        let valid = healpix::nested::to_zuniq(4, 164);
        assert!(grid.center_impl(valid).is_ok());
    }

    #[test]
    fn test_bit_combine_rejects_out_of_range_coordinates() {
        // ring is the dangerous one: an out-of-range nested hash used to
        // panic inside `to_ring` ("assertion failed: j_d0h <= 2")
        for scheme in [Scheme::Nested, Scheme::Ring, Scheme::Zuniq] {
            let grid = grid(scheme, 1);
            assert_eq!(grid.nside(), 2);

            assert!(grid.bit_combine_impl(0, 2).is_err());
            assert!(grid.bit_combine_impl(256, 0).is_err());
            assert!(grid.bit_combine_impl(1, 1).is_ok());
        }
    }

    #[test]
    fn test_bit_combine_matches_statics() {
        assert_eq!(
            grid(Scheme::Zuniq, 1).bit_combine_impl(0, 1).unwrap(),
            360287970189639680
        );
        assert_eq!(grid(Scheme::Ring, 1).bit_combine_impl(0, 1).unwrap(), 4);
        assert_eq!(grid(Scheme::Nested, 1).bit_combine_impl(0, 1).unwrap(), 2);
    }
}
