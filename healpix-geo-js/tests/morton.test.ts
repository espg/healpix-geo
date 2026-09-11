import * as healpixGeo from "../pkg/healpix_geo.js";
import { Grid } from "../pkg/healpix_geo.js";
import { describe, expect, test } from "vitest";

// golden words generated with the published mortie-core 0.1.0:
// from_nested(nested, level)
const BASE0_L0 = 1152921504606846976n; // from_nested(0, 0)
const CELL164_L3 = 4107282860161892355n; // from_nested(164, 3)
const CELL41_L2 = 4107282860161892354n; // from_nested(41, 2) == parent of 164@3
const CELL41_L1 = 4035225266123964417n; // from_nested(10, 1)

describe("morton statics", () => {
  test("level reads the embedded refinement level", () => {
    expect(healpixGeo.morton.level(BASE0_L0)).to.equal(0);
    expect(healpixGeo.morton.level(CELL164_L3)).to.equal(3);
  });

  test("level rejects invalid words", () => {
    expect(() => healpixGeo.morton.level(0n)).to.throw(
      "not a valid morton cell id",
    );
    expect(() => healpixGeo.morton.level(2n ** 64n - 1n)).to.throw();
  });

  test("parent truncates to the requested level", () => {
    expect(healpixGeo.morton.parent(CELL164_L3, 2)).to.equal(CELL41_L2);
    expect(healpixGeo.morton.parent(CELL164_L3, 1)).to.equal(CELL41_L1);
    // at or above the embedded level, the id is unchanged
    expect(healpixGeo.morton.parent(CELL164_L3, 3)).to.equal(CELL164_L3);
    expect(healpixGeo.morton.parent(CELL164_L3, 29)).to.equal(CELL164_L3);
  });

  test("parent validates its inputs", () => {
    expect(() => healpixGeo.morton.parent(0n, 2)).to.throw();
    expect(() => healpixGeo.morton.parent(CELL164_L3, 2.5)).to.throw();
    expect(() => healpixGeo.morton.parent(CELL164_L3, 30)).to.throw();
    expect(() => healpixGeo.morton.parent(CELL164_L3, -1)).to.throw();
  });

  test("contains is containment-is-truncation", () => {
    expect(healpixGeo.morton.contains(CELL41_L2, CELL164_L3)).to.equal(true);
    expect(healpixGeo.morton.contains(CELL41_L1, CELL164_L3)).to.equal(true);
    // not symmetric; reflexive
    expect(healpixGeo.morton.contains(CELL164_L3, CELL41_L2)).to.equal(false);
    expect(healpixGeo.morton.contains(CELL164_L3, CELL164_L3)).to.equal(true);
    // a different level-2 cell does not contain it
    expect(healpixGeo.morton.contains(BASE0_L0, CELL164_L3)).to.equal(false);
  });

  test("contains validates its inputs", () => {
    expect(() => healpixGeo.morton.contains(0n, CELL164_L3)).to.throw();
    expect(() =>
      healpixGeo.morton.contains(CELL164_L3, 2n ** 64n - 1n),
    ).to.throw();
  });
});

describe("morton grid", () => {
  const grid = new Grid({ scheme: "morton", level: 3 });
  const nested = new Grid({ scheme: "nested", level: 3 });

  test("coordinates match the nested scheme", () => {
    const expected = nested.healpixToLonLat(new BigUint64Array([164n]));
    const actual = grid.healpixToLonLat(new BigUint64Array([CELL164_L3]));
    expect(Array.from(actual)).to.deep.equal(Array.from(expected));
  });

  test("lonLatToHealpix encodes at the grid's level", () => {
    const centers = grid.healpixToLonLat(new BigUint64Array([CELL164_L3]));
    const roundtrip = grid.lonLatToHealpix(centers);
    expect(Array.from(roundtrip)).to.deep.equal([CELL164_L3]);
  });

  test("ids are read at their embedded level", () => {
    // a level-0 id in a level-3 grid: the base-cell vertex, as in the zuniq
    // scheme
    const vertex = grid.vertex(BASE0_L0, 0.0, 0.0);
    expect(vertex.lon).to.be.closeTo(45.0, 1e-4);
    expect(vertex.lat).to.be.closeTo(0.0, 1e-4);
  });

  test("invalid and non-canonical ids throw instead of trapping", () => {
    expect(() => grid.vertex(0n, 0.5, 0.5)).to.throw(
      "not a valid morton cell id",
    );
    expect(() => grid.vertex(2n ** 64n - 1n, 0.5, 0.5)).to.throw();
    // junk below the encoded level makes the word non-canonical
    expect(() => grid.vertex(CELL164_L3 | (1n << 30n), 0.5, 0.5)).to.throw();
  });

  test("toScheme round-trips through nested and zuniq", () => {
    expect(grid.toScheme(CELL164_L3, "nested")).to.equal(164n);
    expect(nested.toScheme(164n, "morton")).to.equal(CELL164_L3);

    const zuniq = grid.toScheme(CELL164_L3, "zuniq");
    const zuniqGrid = new Grid({ scheme: "zuniq", level: 3 });
    expect(zuniqGrid.toScheme(zuniq, "morton")).to.equal(CELL164_L3);

    // identity, as for zuniq
    expect(grid.toScheme(CELL164_L3, "morton")).to.equal(CELL164_L3);
  });

  test("toScheme accepts a level override only from level-free schemes", () => {
    expect(nested.toScheme(164n, "morton", 3)).to.equal(CELL164_L3);
    expect(() => grid.toScheme(CELL164_L3, "morton", 3)).to.throw();
    expect(() => nested.toScheme(164n, "ring", 3)).to.throw();
  });

  test("raw unsigned order sorts parents immediately before children", () => {
    // the property gridlook's hive addressing relies on
    expect(CELL41_L1 < CELL41_L2).to.equal(true);
    expect(CELL41_L2 < CELL164_L3).to.equal(true);
  });
});
