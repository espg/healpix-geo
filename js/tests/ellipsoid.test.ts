import init, * as healpixGeo from "../pkg/index.js";
import { Ellipsoid } from "../pkg/index.js";
import { describe, expect, test } from "vitest";

describe("Ellipsoid.from", () => {
  test("null gives the default sphere", () => {
    const result = Ellipsoid.from(null);
    expect(result.isSphere).to.equal(true);
    expect(result.semiMajorAxis).to.equal(6370997.0);
    expect(result.flattening).to.equal(0);
  });

  test("undefined gives the default sphere", () => {
    // upstream's parseEllipsoid rejected `undefined`; documented in the
    // doc comment of `Ellipsoid.from`
    const explicit = Ellipsoid.from(undefined);
    expect(explicit.isSphere).to.equal(true);
    expect(explicit.semiMajorAxis).to.equal(6370997.0);

    // a missing argument behaves the same way
    const missing = (Ellipsoid.from as () => Ellipsoid)();
    expect(missing.isSphere).to.equal(true);
  });

  test("sphere", () => {
    const result = Ellipsoid.from({ name: "sphere", radius: 6371000.0 });
    expect(result.isSphere).to.equal(true);
    expect(result.semiMajorAxis).to.equal(6371000);
  });

  test("ellipsoid from semi-major axis and inverse flattening", () => {
    const result = Ellipsoid.from({
      name: "WGS84",
      semi_major_axis: 6378137.0,
      inverse_flattening: 298.257223563,
    });
    expect(result.isSphere).to.equal(false);
    expect(result.semiMajorAxis).to.equal(6378137.0);
    expect(result.flattening).to.be.closeTo(1 / 298.257223563, 1e-15);
  });

  test("ellipsoid from semi-major axis and semi-minor axis", () => {
    const result = Ellipsoid.from({
      name: "WGS84",
      semi_major_axis: 6378137.0,
      semi_minor_axis: 6356752.314245179,
    });
    expect(result.isSphere).to.equal(false);
    expect(result.semiMajorAxis).to.equal(6378137.0);
    expect(result.flattening).to.be.closeTo(1 / 298.257223563, 1e-9);
  });

  test("invalid input throws", () => {
    expect(() => Ellipsoid.from({ bogus: 1 })).to.throw();
  });
});

describe("parseEllipsoid (deprecated alias, breaking return type)", () => {
  test("returns a reusable handle instead of the input-shaped union", () => {
    const result = healpixGeo.parseEllipsoid({ radius: 6371000.0 });
    expect(result).to.be.instanceOf(Ellipsoid);
    expect(result.semiMajorAxis).to.equal(6371000);

    // the 0.2.x property surface is gone — this is the breaking part of the
    // change, not a compatibility shim
    expect(result).to.not.have.property("radius");
    expect(result).to.not.have.property("semi_major_axis");
    expect(result).to.not.have.property("inverse_flattening");
  });
});

describe("handle reuse", () => {
  test("the same handle survives repeated calls", () => {
    const ellipsoid = Ellipsoid.from({
      semi_major_axis: 6378137.0,
      inverse_flattening: 298.257223563,
    });

    // before the handle existed, the second call would throw
    // "invalid dynamic union value" because the union argument was
    // consumed (moved into WASM) by the first call
    const first = healpixGeo.nested.healpixToLonLat(164n, 4, ellipsoid);
    const second = healpixGeo.nested.healpixToLonLat(164n, 4, ellipsoid);

    expect(first.lon).to.equal(second.lon);
    expect(first.lat).to.equal(second.lat);
    expect(ellipsoid.semiMajorAxis).to.equal(6378137.0);
  });
});
