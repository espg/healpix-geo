use cdshealpix as healpix;
use healpix_geo_core::ellipsoid::{Ellipsoid as RustEllipsoid, ReferenceBody};
use healpix_geo_core::scalar::nested::coordinates as scalar;
use serde::Deserialize;
use serde_wasm_bindgen::from_value;
use wasm_bindgen::prelude::*;

use crate::coordinates::Coordinate;
use crate::ellipsoid::EllipsoidLike;
use crate::geometry::spherical_vertex;

const MAX_DEPTH: u8 = 29;

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

/// A HEALPix grid at a fixed depth on a fixed reference body.
///
/// Constructed through [`create_grid`] (`createGrid` in JS). The scheme
/// dispatch, the depth, and the parsed ellipsoid state (including the
/// authalic-latitude coefficients) are stored once, so per-call overhead is
/// limited to the actual coordinate math.
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
            return Err(format!("`depth` must be at most {}, got {}", MAX_DEPTH, depth));
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

    /// The cell in the nested scheme plus the depth it lives at.
    fn to_nested(&self, cell: u64) -> (u8, u64) {
        match self.scheme {
            Scheme::Nested => (self.depth, cell),
            Scheme::Ring => (self.depth, healpix::nested::get(self.depth).from_ring(cell)),
            Scheme::Zuniq => healpix::nested::from_zuniq(cell),
        }
    }

    fn from_nested(&self, hash: u64) -> u64 {
        match self.scheme {
            Scheme::Nested => hash,
            Scheme::Ring => healpix::nested::get(self.depth).to_ring(hash),
            Scheme::Zuniq => healpix::nested::to_zuniq(self.depth, hash),
        }
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

    /// The reference body of the grid, as a new handle
    #[wasm_bindgen(getter)]
    pub fn ellipsoid(&self) -> crate::ellipsoid::Ellipsoid {
        crate::ellipsoid::Ellipsoid {
            inner: self.ellipsoid.clone(),
        }
    }

    /// Center coordinates of the given cell
    pub fn center(&self, cell: u64) -> Coordinate {
        let (depth, hash) = self.to_nested(cell);
        let layer = healpix::nested::get(depth);

        let (lon, lat) = scalar::healpix_to_lonlat(&hash, layer, &self.ellipsoid);

        Coordinate { lon, lat }
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
    pub fn vertex(&self, cell: u64, u: f64, v: f64) -> Coordinate {
        let (depth, hash) = self.to_nested(cell);
        let layer = healpix::nested::get(depth);

        let center = layer.center_of_projected_cell(hash);
        let (lon, lat) = spherical_vertex(center, depth, (u, v));

        Coordinate {
            lon: lon.to_degrees().rem_euclid(360.0),
            lat: self
                .ellipsoid
                .latitude_authalic_to_geographic(lat)
                .to_degrees(),
        }
    }

    /// Cell index at the given z-order coordinates
    ///
    /// Interleaves the bits of `i` and `j` (the two axes of the nested
    /// z-order numbering within a base-resolution pixel) into the cell index
    /// at the grid's depth, expressed in the grid's scheme.
    #[wasm_bindgen(js_name = bitCombine)]
    pub fn bit_combine(&self, i: u32, j: u32) -> u64 {
        let zoc = healpix::nested::zordercurve::get_zoc(self.depth);

        self.from_nested(zoc.ij2h(i, j))
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
pub fn create_grid(options: JsValue) -> Result<Grid, JsValue> {
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
        let grid = Grid::from_options(options(Scheme::Nested, Some(0), None)).unwrap();

        let expected: Vec<(f64, f64)> = vec![
            (45.0, 0.0),
            (67.5, 19.47122063),
            (90.0, 41.8103149),
            (45.0, 90.0),
        ];
        let uv: Vec<(f64, f64)> = vec![(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.0, 1.0)];

        for ((u, v), (lon, lat)) in uv.into_iter().zip(expected) {
            let actual = grid.vertex(0, u, v);
            assert!((actual.lon - lon).abs() < 1e-4);
            assert!((actual.lat - lat).abs() < 1e-4);
        }
    }

    #[test]
    fn test_zuniq_uses_encoded_depth() {
        // depth 0, nested cell 0 encoded as zuniq
        let cell = healpix::nested::to_zuniq(0, 0);
        let grid = Grid::from_options(options(Scheme::Zuniq, Some(4), None)).unwrap();

        let vertex = grid.vertex(cell, 0.0, 0.0);
        assert!((vertex.lon - 45.0).abs() < 1e-4);
        assert!(vertex.lat.abs() < 1e-4);
    }

    #[test]
    fn test_ring_roundtrip() {
        let grid = Grid::from_options(options(Scheme::Ring, Some(4), None)).unwrap();

        let center = grid.center(164);
        let cell = grid.cell_at(center.lon, center.lat);
        assert_eq!(cell, 164);
    }

    #[test]
    fn test_bit_combine_matches_statics() {
        let grid = Grid::from_options(options(Scheme::Zuniq, Some(1), None)).unwrap();
        assert_eq!(grid.bit_combine(0, 1), 360287970189639680);

        let grid = Grid::from_options(options(Scheme::Ring, Some(1), None)).unwrap();
        assert_eq!(grid.bit_combine(0, 1), 4);

        let grid = Grid::from_options(options(Scheme::Nested, Some(1), None)).unwrap();
        assert_eq!(grid.bit_combine(0, 1), 2);
    }
}
