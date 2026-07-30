import init, * as healpixGeo from "../pkg/index.js";
import { describe, expect, test } from "vitest";

const wgs84 = {
  semi_major_axis: 6378137.0,
  inverse_flattening: 298.257223563,
};

describe("Grid.vertices", () => {
  test("matches the scalar vertex loop", () => {
    const steps = 5;
    for (const ellipsoid of [null, wgs84]) {
      const grid = healpixGeo.createGrid({
        scheme: "nested",
        depth: 4,
        ellipsoid,
      });

      const bulk = grid.vertices(164n, steps);
      expect(bulk).to.be.instanceOf(Float64Array);
      expect(bulk.length).to.equal(steps * steps * 2);

      let k = 0;
      for (let i = 0; i < steps; ++i) {
        const u = i / (steps - 1);
        for (let j = 0; j < steps; ++j) {
          const v = j / (steps - 1);
          const scalar = grid.vertex(164n, u, v);
          expect(bulk[k++]).to.equal(scalar.lon);
          expect(bulk[k++]).to.equal(scalar.lat);
        }
      }
    }
  });

  test("works for zuniq cells (depth from the id)", () => {
    const grid = healpixGeo.createGrid({ scheme: "zuniq", depth: 4 });
    const cell = healpixGeo.zuniq.bitCombine(4, 3, 5);

    const bulk = grid.vertices(cell, 3);
    const corner = grid.vertex(cell, 0, 0);
    expect(bulk[0]).to.equal(corner.lon);
    expect(bulk[1]).to.equal(corner.lat);
  });

  test("rejects steps < 2", () => {
    const grid = healpixGeo.createGrid({ scheme: "nested", depth: 4 });
    expect(() => grid.vertices(164n, 1)).to.throw(/at least 2/);
  });

  test("cells above 2^53 stay exact", () => {
    const grid = healpixGeo.createGrid({ scheme: "nested", depth: 29 });
    const cell = 11n * 4n ** 29n - 1n;
    expect(cell > BigInt(Number.MAX_SAFE_INTEGER)).to.equal(true);

    const bulk = grid.vertices(cell, 2);
    const corner = grid.vertex(cell, 0, 0);
    expect(bulk[0]).to.equal(corner.lon);
    expect(bulk[1]).to.equal(corner.lat);
  });
});

describe("Grid.cellsAt and Grid.centers", () => {
  test("roundtrip through typed arrays", () => {
    const grid = healpixGeo.createGrid({
      scheme: "ring",
      depth: 4,
      ellipsoid: wgs84,
    });
    const cells = new BigUint64Array([0n, 164n, 700n]);

    const centers = grid.centers(cells);
    expect(centers).to.be.instanceOf(Float64Array);
    expect(centers.length).to.equal(cells.length * 2);

    const roundtrip = grid.cellsAt(centers);
    expect(roundtrip).to.be.instanceOf(BigUint64Array);
    expect([...roundtrip]).to.deep.equal([...cells]);
  });

  test("cellsAt matches the scalar static", () => {
    const grid = healpixGeo.createGrid({ scheme: "nested", depth: 4 });
    const sphere = healpixGeo.Ellipsoid.from(null);

    const cells = grid.cellsAt(new Float64Array([16.875, 38.68, 45.0, 0.0]));
    expect(cells[0]).to.equal(
      healpixGeo.nested.lonLatToHealpix(16.875, 38.68, 4, sphere),
    );
    expect(cells[1]).to.equal(
      healpixGeo.nested.lonLatToHealpix(45.0, 0.0, 4, sphere),
    );
  });

  test("cellsAt rejects odd-length input", () => {
    const grid = healpixGeo.createGrid({ scheme: "nested", depth: 4 });
    expect(() => grid.cellsAt(new Float64Array([1, 2, 3]))).to.throw(
      /even length/,
    );
  });
});

describe("Grid.bitCombineTable", () => {
  test("matches the scalar bitCombine (gridlook unshuffle layout)", () => {
    const size = 8;
    const grid = healpixGeo.createGrid({ scheme: "nested", nside: size });

    const table = grid.bitCombineTable(size);
    expect(table).to.be.instanceOf(BigUint64Array);
    expect(table.length).to.equal(size * size);

    // gridlook's loop: for i { for j { temp[idx++] = bitCombine(j, i) } }
    let idx = 0;
    for (let i = 0; i < size; ++i) {
      for (let j = 0; j < size; ++j) {
        expect(table[idx++]).to.equal(healpixGeo.nested.bitCombine(3, j, i));
      }
    }
  });
});
