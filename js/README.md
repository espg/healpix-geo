# javascript and typescript bindings for `healpix-geo`

This module provides javascript / typescript bindings to the `healpix_geo_core::scalar` crate.

Usage:

```typescript
import init, * as healpixGeo from "healpix-geo";

await init(); // or use a bundler

// parse once, reuse for any number of calls
const ellipsoid = healpixGeo.Ellipsoid.from({
  semi_major_axis: 6378137.0,
  inverse_flattening: 298.257223563,
}); // or Ellipsoid.from(null) for the default sphere

const cellId: bigint = 10;
const level: number = 2;
const { lon, lat } = healpixGeo.nested.healpixToLonLat(
  cellId,
  level,
  ellipsoid,
);
```

### Migrating from 0.2.x

- The ellipsoid argument of the `nested` / `ring` / `zuniq` functions is now
  required and takes an `Ellipsoid` handle: pass `Ellipsoid.from(null)` where
  `null` used to be passed.
- `parseEllipsoid` keeps its name and input shapes, but returns the handle
  instead of a `Sphere` / `EllipsoidInverseFlattening` /
  `EllipsoidSemiMinorAxis` value — a breaking return-type change. Code that
  only parses and passes the result is unaffected (and gains reuse); code that
  read `radius` / `semi_major_axis` / `inverse_flattening` /
  `semi_minor_axis` off the result reads `semiMajorAxis`, `flattening` and
  `isSphere` instead. Those three classes are no longer exported.
