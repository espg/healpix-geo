import init, * as healpixGeo from "../pkg/index.js";
import { Coordinate, Ellipsoid, Grid } from "../pkg/index.js";
import { describe, expect, test } from "vitest";

const wgs84 = {
  semi_major_axis: 6378137.0,
  inverse_flattening: 298.257223563,
};

describe("createGrid", () => {
  test("basic construction", () => {
    const grid = healpixGeo.createGrid({ scheme: "nested", depth: 4 });
    expect(grid.scheme).to.equal("nested");
    expect(grid.depth).to.equal(4);
    expect(grid.nside).to.equal(16);
    expect(grid.ellipsoid.isSphere).to.equal(true);
  });

  test("accepts nside instead of depth", () => {
    const grid = healpixGeo.createGrid({ scheme: "ring", nside: 8 });
    expect(grid.depth).to.equal(3);
    expect(grid.nside).to.equal(8);
  });

  test("accepts an ellipsoid definition", () => {
    const grid = healpixGeo.createGrid({
      scheme: "nested",
      depth: 4,
      ellipsoid: wgs84,
    });
    expect(grid.ellipsoid.isSphere).to.equal(false);
    expect(grid.ellipsoid.semiMajorAxis).to.equal(wgs84.semi_major_axis);
  });

  test("flattened ellipsoid getters avoid the handle allocation", () => {
    // `grid.ellipsoid` mints a new handle (allocation + finalizer) on every
    // read; these return plain numbers
    const grid = healpixGeo.createGrid({
      scheme: "nested",
      depth: 4,
      ellipsoid: wgs84,
    });

    const handle = grid.ellipsoid;
    expect(grid.semiMajorAxis).to.equal(handle.semiMajorAxis);
    expect(grid.flattening).to.equal(handle.flattening);
    expect(grid.isSphere).to.equal(handle.isSphere);
    expect(grid.isSphere).to.equal(false);
    handle.free();
  });

  test("rejects invalid options", () => {
    expect(() => healpixGeo.createGrid({ scheme: "nested" })).to.throw(
      /depth.*nside/,
    );
    expect(() =>
      healpixGeo.createGrid({ scheme: "nested", depth: 4, nside: 16 }),
    ).to.throw(/not both/);
    expect(() =>
      healpixGeo.createGrid({ scheme: "nested", nside: 3 }),
    ).to.throw(/power of two/);
    expect(() =>
      healpixGeo.createGrid({ scheme: "nested", depth: 30 }),
    ).to.throw(/at most 29/);
    expect(() =>
      healpixGeo.createGrid({ scheme: "bogus", depth: 4 }),
    ).to.throw();
  });
});

describe("Grid input validation", () => {
  // every one of these used to be an unrecoverable wasm trap, which poisons
  // the whole module instance rather than throwing
  test("rejects cell ids out of range at the grid's depth", () => {
    for (const scheme of ["nested", "ring"] as const) {
      const grid = healpixGeo.createGrid({ scheme, depth: 4 });
      // 12 * 4^4 == 3072
      expect(() => grid.center(3072n)).to.throw(/out of range/);
      expect(() => grid.vertex(3072n, 0.5, 0.5)).to.throw(/out of range/);
      expect(() => grid.center(2n ** 64n - 1n)).to.throw(/out of range/);
      expect(grid.center(3071n).lon).to.be.a("number");
    }
  });

  test("rejects malformed zuniq cell ids", () => {
    const grid = healpixGeo.createGrid({ scheme: "zuniq", depth: 4 });

    expect(() => grid.center(0n)).to.throw(/zuniq/);
    expect(() => grid.center(2n)).to.throw(/zuniq/);
    expect(() => grid.center(2n ** 64n - 1n)).to.throw(/zuniq/);
  });

  test("rejects z-order coordinates outside the base cell", () => {
    for (const scheme of ["nested", "ring", "zuniq"] as const) {
      const grid = healpixGeo.createGrid({ scheme, depth: 1 });
      expect(() => grid.bitCombine(0, 2)).to.throw(/out of range/);
      expect(() => grid.bitCombine(256, 0)).to.throw(/out of range/);
      expect(grid.bitCombine(1, 1)).to.be.a("bigint");
    }
  });

  test("the module stays usable after a rejected call", () => {
    const grid = healpixGeo.createGrid({ scheme: "ring", depth: 1 });
    expect(() => grid.bitCombine(256, 0)).to.throw();
    expect(grid.center(0n).lon).to.be.a("number");
  });
});

describe("Grid parity with the scheme statics", () => {
  test("nested", () => {
    const grid = healpixGeo.createGrid({
      scheme: "nested",
      depth: 4,
      ellipsoid: wgs84,
    });
    const ellipsoid = Ellipsoid.from(wgs84);

    const center: Coordinate = grid.center(164n);
    const staticCenter = healpixGeo.nested.healpixToLonLat(164n, 4, ellipsoid);
    expect(center.lon).to.equal(staticCenter.lon);
    expect(center.lat).to.equal(staticCenter.lat);

    const vertex = grid.vertex(164n, 0.25, 0.75);
    const staticVertex = healpixGeo.nested.vertex(
      164n,
      4,
      0.25,
      0.75,
      ellipsoid,
    );
    expect(vertex.lon).to.equal(staticVertex.lon);
    expect(vertex.lat).to.equal(staticVertex.lat);

    const cell = grid.cellAt(center.lon, center.lat);
    expect(cell).to.equal(164n);

    expect(grid.bitCombine(3, 5)).to.equal(
      healpixGeo.nested.bitCombine(4, 3, 5),
    );
  });

  test("ring", () => {
    const grid = healpixGeo.createGrid({ scheme: "ring", depth: 4 });
    const sphere = Ellipsoid.from(null);

    const center = grid.center(164n);
    const staticCenter = healpixGeo.ring.healpixToLonLat(164n, 4, sphere);
    expect(center.lon).to.equal(staticCenter.lon);
    expect(center.lat).to.equal(staticCenter.lat);

    const vertex = grid.vertex(164n, 0.5, 0.5);
    const staticVertex = healpixGeo.ring.vertex(164n, 4, 0.5, 0.5, sphere);
    expect(vertex.lon).to.equal(staticVertex.lon);
    expect(vertex.lat).to.equal(staticVertex.lat);
  });

  test("zuniq hides the signature difference", () => {
    // the zuniq statics take no depth for vertex/center; the grid object
    // exposes the exact same call shape for all three schemes
    const grid = healpixGeo.createGrid({ scheme: "zuniq", depth: 0 });
    const sphere = Ellipsoid.from(null);

    const cell = healpixGeo.zuniq.bitCombine(0, 0, 0);

    const vertex = grid.vertex(cell, 1.0, 1.0);
    const staticVertex = healpixGeo.zuniq.vertex(cell, 1.0, 1.0, sphere);
    expect(vertex.lon).to.equal(staticVertex.lon);
    expect(vertex.lat).to.equal(staticVertex.lat);

    const center = grid.center(cell);
    const staticCenter = healpixGeo.zuniq.healpixToLonLat(cell, sphere);
    expect(center.lon).to.equal(staticCenter.lon);
    expect(center.lat).to.equal(staticCenter.lat);
  });

  test("uniform loop over all schemes", () => {
    // the motivating gridlook use case: one code path for every scheme
    for (const scheme of ["nested", "ring", "zuniq"] as const) {
      const grid = healpixGeo.createGrid({ scheme, depth: 2 });
      const cell = grid.bitCombine(1, 2);
      const vertex = grid.vertex(cell, 0.5, 0.5);
      const center = grid.center(cell);
      expect(vertex.lon).to.equal(center.lon);
      expect(vertex.lat).to.equal(center.lat);
    }
  });
});

describe("Grid ellipsoid state reuse", () => {
  test("repeated calls do not need re-parsing", () => {
    const grid = healpixGeo.createGrid({
      scheme: "nested",
      depth: 4,
      ellipsoid: wgs84,
    });

    const first = grid.vertex(164n, 0.5, 0.5);
    const second = grid.vertex(164n, 0.5, 0.5);
    expect(first.lon).to.equal(second.lon);
    expect(first.lat).to.equal(second.lat);
  });

  test("cells above 2^53 stay exact", () => {
    // depth 29 nested ids exceed Number.MAX_SAFE_INTEGER
    const grid = healpixGeo.createGrid({ scheme: "nested", depth: 29 });
    const cell = 11n * 4n ** 29n - 1n; // last cell at depth 29
    expect(cell > BigInt(Number.MAX_SAFE_INTEGER)).to.equal(true);

    const center = grid.center(cell);
    const roundtrip = grid.cellAt(center.lon, center.lat);
    expect(roundtrip).to.equal(cell);
  });
});
