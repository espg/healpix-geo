import { Grid } from "../pkg/healpix_geo.js";
import { describe, expect, test } from "vitest";

const nested = (level: number) => new Grid({ scheme: "nested", level });

describe("Grid.neighbours", () => {
  test("returns eight positional ids per cell by default", () => {
    // level 2, cell 42: the edge neighbours are SW 111, NW 21, NE 43, SE 40
    const all = nested(2).neighbours(new BigUint64Array([42n]));
    expect(all).to.be.instanceOf(BigInt64Array);
    expect(all.length).to.equal(8);
    // S, SW, W, NW, N, NE, E, SE
    expect([all[1], all[3], all[5], all[7]]).to.deep.equal([
      111n,
      21n,
      43n,
      40n,
    ]);
  });

  test("edge and vertex connectivity split the eight", () => {
    const grid = nested(2);
    const cells = new BigUint64Array([42n, 7n]);
    const all = grid.neighbours(cells, "all");
    const edge = grid.neighbours(cells, "edge");
    const vertex = grid.neighbours(cells, "vertex");

    expect(edge.length).to.equal(8);
    expect(vertex.length).to.equal(8);
    for (let cell = 0; cell < cells.length; ++cell) {
      const row = Array.from(all.subarray(cell * 8, cell * 8 + 8));
      expect(Array.from(edge.subarray(cell * 4, cell * 4 + 4))).to.deep.equal([
        row[1],
        row[3],
        row[5],
        row[7],
      ]);
      expect(Array.from(vertex.subarray(cell * 4, cell * 4 + 4))).to.deep.equal(
        [row[0], row[2], row[4], row[6]],
      );
    }
  });

  test("crosses base cell boundaries and is symmetric", () => {
    const level = 3;
    const grid = nested(level);
    const perBase = 4n ** BigInt(level);
    const cells = BigUint64Array.from(
      { length: Number(12n * perBase) },
      (_, index) => BigInt(index),
    );
    const table = grid.neighbours(cells);

    let missing = 0;
    let crossings = 0;
    for (let cell = 0; cell < cells.length; ++cell) {
      for (let k = 0; k < 8; ++k) {
        const other = table[cell * 8 + k];
        if (other === -1n) {
          missing += 1;
          continue;
        }
        if (other / perBase !== BigInt(cell) / perBase) {
          crossings += 1;
        }
        // every neighbour lists this cell back
        const back = table.subarray(Number(other) * 8, Number(other) * 8 + 8);
        expect(Array.from(back)).to.include(BigInt(cell));
      }
    }

    // the 24 cells at the corners where only three base cells meet
    expect(missing).to.equal(24);
    expect(crossings).to.be.greaterThan(0);
  });

  test("answers in the grid's scheme", () => {
    const level = 3;
    const reference = nested(level);
    const cell = 300n;
    const expected = reference.neighbours(new BigUint64Array([cell]));

    for (const scheme of ["ring", "zuniq"] as const) {
      const grid = new Grid({ scheme, level });
      const actual = grid.neighbours(
        new BigUint64Array([reference.toScheme(cell, scheme)]),
      );
      expect(Array.from(actual)).to.deep.equal(
        Array.from(expected, (id) =>
          id === -1n ? -1n : reference.toScheme(id, scheme),
        ),
      );
    }
  });

  test("a zuniq cell's neighbours live at its own level", () => {
    const fine = nested(5);
    const cell = 1234n;
    const zuniq = new Grid({ scheme: "zuniq", level: 1 });
    const actual = zuniq.neighbours(
      new BigUint64Array([fine.toScheme(cell, "zuniq")]),
    );
    const expected = fine.neighbours(new BigUint64Array([cell]));
    expect(Array.from(actual)).to.deep.equal(
      Array.from(expected, (id) => fine.toScheme(id, "zuniq")),
    );
  });

  test("rejects bad input as a catchable error", () => {
    const grid = nested(1);
    expect(() => grid.neighbours(new BigUint64Array([48n]))).to.throw(
      /cells\[0\].*out of range/,
    );
    expect(() =>
      // @ts-expect-error not a Connectivity
      grid.neighbours(new BigUint64Array([0n]), "diagonal"),
    ).to.throw(/unknown connectivity/);
  });
});

describe("Grid.kthNeighbourhood and Grid.kthNeighbours", () => {
  test("ring 1 is the cell plus its neighbours", () => {
    const grid = nested(4);
    const cell = 777n;
    const hood = grid.kthNeighbourhood(new BigUint64Array([cell]), 1);
    const ring = grid.kthNeighbours(new BigUint64Array([cell]), 1);
    const direct = grid.neighbours(new BigUint64Array([cell]));

    expect(hood.length).to.equal(9);
    expect(ring.length).to.equal(8);
    const sorted = (ids: Iterable<bigint>) =>
      Array.from(ids)
        .filter((id) => id !== -1n)
        .sort((a, b) => (a < b ? -1 : 1));
    expect(sorted(ring)).to.deep.equal(sorted(direct));
    expect(sorted(hood)).to.deep.equal(sorted([cell, ...direct]));
  });

  test("the neighbourhood is the union of the rings", () => {
    const grid = nested(4);
    const cells = new BigUint64Array([0n, 777n, 3071n]);
    const ring = 3;
    const side = 2 * ring + 1;
    const hood = grid.kthNeighbourhood(cells, ring);
    expect(hood.length).to.equal(cells.length * side * side);

    for (let index = 0; index < cells.length; ++index) {
      const one = new BigUint64Array([cells[index]]);
      const union = new Set<bigint>([cells[index]]);
      for (let k = 1; k <= ring; ++k) {
        const ids = grid.kthNeighbours(one, k);
        expect(ids.length).to.equal(8 * k);
        ids.forEach((id) => id !== -1n && union.add(id));
      }
      const row = hood.subarray(index * side * side, (index + 1) * side * side);
      expect(new Set(Array.from(row).filter((id) => id !== -1n))).to.deep.equal(
        union,
      );
    }
  });

  test("ring 0 is the cell itself", () => {
    const hood = nested(3).kthNeighbourhood(new BigUint64Array([5n, 9n]), 0);
    expect(Array.from(hood)).to.deep.equal([5n, 9n]);
  });

  test("rejects a ring the grid cannot hold", () => {
    const grid = nested(2);
    const cells = new BigUint64Array([0n]);
    expect(() => grid.kthNeighbours(cells, 0)).to.throw(/at least 1/);
    expect(() => grid.kthNeighbours(cells, 1.5)).to.throw(/integer/);
    expect(() => grid.kthNeighbourhood(cells, 5)).to.throw(/at most `nside`/);
    // the widest ring the grid accepts does not abort the module
    expect(grid.kthNeighbourhood(cells, 4).length).to.equal(81);
  });
});
