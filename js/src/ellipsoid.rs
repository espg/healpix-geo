use geodesy::ellps::Ellipsoid as GeodesyEllipsoid;
use geodesy::prelude::EllipsoidBase;
use healpix_geo_core::ellipsoid::{
    Ellipsoid as RustEllipsoid, ReferenceBody, ReferenceEllipsoid, ReferenceSphere,
};
use serde::Deserialize;
use serde_wasm_bindgen::from_value;
use wasm_bindgen::prelude::*;

/// Plain-object shapes accepted as ellipsoid definitions.
///
/// This is an input (deserialization) type only; the parsed, reusable state
/// lives in [`Ellipsoid`].
#[derive(Deserialize, Debug, PartialEq)]
pub enum EllipsoidLike {
    #[serde(untagged)]
    EllipsoidInverseFlattening(EllipsoidInverseFlattening),
    #[serde(untagged)]
    EllipsoidSemiMinorAxis(EllipsoidSemiMinorAxis),
    #[serde(untagged)]
    Sphere(Sphere),
}

#[derive(Deserialize, Debug, PartialEq)]
pub struct EllipsoidInverseFlattening {
    pub semi_major_axis: f64,
    pub inverse_flattening: f64,
}

#[derive(Deserialize, Debug, PartialEq)]
pub struct EllipsoidSemiMinorAxis {
    pub semi_major_axis: f64,
    pub semi_minor_axis: f64,
}

impl From<EllipsoidSemiMinorAxis> for EllipsoidInverseFlattening {
    fn from(val: EllipsoidSemiMinorAxis) -> EllipsoidInverseFlattening {
        let a = val.semi_major_axis;
        let b = val.semi_minor_axis;

        EllipsoidInverseFlattening {
            semi_major_axis: a,
            inverse_flattening: a / (a - b),
        }
    }
}

#[derive(Deserialize, Debug, PartialEq)]
pub struct Sphere {
    pub radius: f64,
}

impl Default for Sphere {
    fn default() -> Sphere {
        Sphere { radius: 6370997.0 }
    }
}

impl EllipsoidLike {
    pub fn into_ellipsoid(self) -> RustEllipsoid {
        match self {
            Self::EllipsoidInverseFlattening(ell) => {
                let ellipsoid =
                    GeodesyEllipsoid::new(ell.semi_major_axis, 1.0f64 / ell.inverse_flattening);

                RustEllipsoid::Ellipsoid(ReferenceEllipsoid::new(ellipsoid))
            }
            Self::EllipsoidSemiMinorAxis(ell) => {
                let a = ell.semi_major_axis;
                let b = ell.semi_minor_axis;
                let flattening = (a - b) / a;
                let ellipsoid = GeodesyEllipsoid::new(a, flattening);

                RustEllipsoid::Ellipsoid(ReferenceEllipsoid::new(ellipsoid))
            }
            Self::Sphere(sphere) => {
                let ellipsoid = GeodesyEllipsoid::new(sphere.radius, 0.0f64);

                RustEllipsoid::Sphere(ReferenceSphere::new(ellipsoid))
            }
        }
    }
}

/// A parsed, reusable reference body.
///
/// Unlike the previous `parseEllipsoid` result (a dynamic union that was
/// *consumed* — moved into WASM — on first use), this handle is passed to
/// functions by reference and can be reused for any number of calls. Parsing
/// and the authalic-latitude Fourier coefficient computation happen once, at
/// construction.
#[wasm_bindgen]
pub struct Ellipsoid {
    pub(crate) inner: RustEllipsoid,
}

#[wasm_bindgen]
impl Ellipsoid {
    /// Parse a plain object (or `null` for the default sphere) into a
    /// reusable ellipsoid handle.
    ///
    /// Accepted shapes:
    /// - `null` — the default sphere (radius 6370997 m)
    /// - `{ radius }` — a sphere
    /// - `{ semi_major_axis, inverse_flattening }`
    /// - `{ semi_major_axis, semi_minor_axis }`
    #[wasm_bindgen(js_name = from)]
    pub fn from_value(obj: JsValue) -> Result<Ellipsoid, JsValue> {
        let inner = if obj.is_null() || obj.is_undefined() {
            RustEllipsoid::default()
        } else {
            let parsed: EllipsoidLike = from_value(obj)?;
            parsed.into_ellipsoid()
        };

        Ok(Ellipsoid { inner })
    }

    /// Semi-major axis (the radius, for a sphere) in meters
    #[wasm_bindgen(getter, js_name = semiMajorAxis)]
    pub fn semi_major_axis(&self) -> f64 {
        self.inner.ellipsoid().semimajor_axis()
    }

    /// Flattening (0 for a sphere)
    #[wasm_bindgen(getter)]
    pub fn flattening(&self) -> f64 {
        self.inner.ellipsoid().flattening()
    }

    /// Whether this reference body is a sphere
    #[wasm_bindgen(getter, js_name = isSphere)]
    pub fn is_sphere(&self) -> bool {
        matches!(self.inner, RustEllipsoid::Sphere(_))
    }
}

/// Parse a plain object into a reusable [`Ellipsoid`] handle.
///
/// Alias of `Ellipsoid.from`, kept for compatibility with the 0.2.x API.
/// Note: the returned value is now a reusable handle instead of a one-shot
/// union value.
#[wasm_bindgen(js_name = parseEllipsoid)]
pub fn parse_ellipsoid(obj: JsValue) -> Result<Ellipsoid, JsValue> {
    Ellipsoid::from_value(obj)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ellipsoidlike_to_ellipsoid() {
        let a: f64 = 6378137.0;
        let if_: f64 = 298.257223563;
        let f: f64 = 1.0 / if_;

        let obj = EllipsoidLike::EllipsoidInverseFlattening(EllipsoidInverseFlattening {
            semi_major_axis: a,
            inverse_flattening: if_,
        });

        let actual = obj.into_ellipsoid();
        match actual {
            RustEllipsoid::Ellipsoid(ell) => {
                let unpacked = ell.ellipsoid();

                assert_eq!(unpacked.semimajor_axis(), a);
                assert_eq!(unpacked.flattening(), f);
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn test_deserialize_ellipsoid() {
        let a = 6378137.0;
        let if_ = 298.257223563;

        let json = r#"{"name": "WGS84", "semi_major_axis": 6378137.0, "inverse_flattening": 298.257223563}"#;

        let actual: EllipsoidLike = serde_json::from_str(json).unwrap();
        let expected = EllipsoidLike::EllipsoidInverseFlattening(EllipsoidInverseFlattening {
            semi_major_axis: a,
            inverse_flattening: if_,
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_deserialize_ellipsoid_semiminor() {
        let a = 6378137.0;
        let b = 6356752.314245179;

        let json = r#"{"name": "WGS84", "semi_major_axis": 6378137.0, "semi_minor_axis": 6356752.314245179}"#;

        let actual: EllipsoidLike = serde_json::from_str(json).unwrap();
        let expected = EllipsoidLike::EllipsoidSemiMinorAxis(EllipsoidSemiMinorAxis {
            semi_major_axis: a,
            semi_minor_axis: b,
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_deserialize_sphere() {
        let json = r#"{"name": "sphere", "radius": 6370997.0}"#;

        let actual: EllipsoidLike = serde_json::from_str(json).unwrap();
        let expected = EllipsoidLike::Sphere(Sphere { radius: 6370997.0 });

        assert_eq!(actual, expected);
    }
}

#[cfg(all(test, target_arch = "wasm32"))]
pub mod tests_wasm32 {
    use super::*;
    use serde_json::json;
    use serde_wasm_bindgen::to_value;
    use wasm_bindgen_test::*;

    wasm_bindgen_test_configure!(run_in_browser);

    #[wasm_bindgen_test]
    fn test_parse_ellipsoid_ellipsoid() {
        let a: f64 = 6378137.0;
        let if_: f64 = 298.257223563;

        let data = json!({"name": "WGS84", "semi_major_axis": a, "inverse_flattening": if_});
        let obj = to_value(&data).map_err(JsValue::from).unwrap();

        let actual = Ellipsoid::from_value(obj).map_err(JsValue::from).unwrap();
        assert_eq!(actual.semi_major_axis(), a);
        assert!(!actual.is_sphere());
    }
}
