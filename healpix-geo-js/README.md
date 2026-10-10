# javascript and typescript bindings for `healpix-geo`

This module provides javascript / typescript bindings to the `healpix_geo_core::scalar` crate.

## `Grid`

`Grid` is the recommended entry point: it holds the scheme, the refinement
level and the parsed ellipsoid state (including the authalic-latitude
coefficients), so the per-call overhead is the coordinate math alone, and
every scheme gets the same call shape.

```typescript
import init, { Grid } from "healpix-geo";

await init(); // or use a bundler

const grid = new Grid({
  scheme: "nested", // or "ring" or "zuniq"
  level: 4, // 0-indexed; level 0 is the 12 base cells
  ellipsoid: { semi_major_axis: 6378137.0, inverse_flattening: 298.257223563 },
}); // `ellipsoid` is optional; omit it for the default sphere

const corner = grid.vertex(164n, 0, 0); // u, v are cell-local offsets in [0, 1]
const cell = grid.bitCombine(3, 5);
const ring = grid.toScheme(cell, "ring");
const coarse = grid.toScheme(5n, "zuniq", 0); // read (and encode) 5n at level 0
```

Every method reports misuse — a cell id out of range at the grid's level, a
vertex offset outside `[0, 1]`, a non-integer z-order coordinate — as a
catchable JS `Error` rather than trapping the wasm instance.

Three things worth knowing:

- For the `zuniq` scheme, methods that read a cell id use the level **embedded
  in the id**, not the grid's level; methods that produce one encode the
  grid's level. So `zuniq` → `zuniq` is the identity, and a coarse id survives
  a read unchanged.
- `toScheme` takes an optional `level` that overrides the level the cell is
  read and encoded at. It is only valid when converting to a scheme that
  encodes the level in its cell ids (`"zuniq"` today) and from one that does
  not; anything else throws.
- `grid.ellipsoid` allocates a **new** handle on every read (wasm allocation +
  JS wrapper + finalizer). Hoist it out of hot paths and `free()` it, or use
  the `semiMajorAxis` / `flattening` / `isSphere` getters, which return plain
  numbers.

### bulk methods

Four methods run the whole loop inside WASM and hand back a single typed
array, which is what makes a per-vertex mesh build practical:

```typescript
const grid = new Grid({ scheme: "nested", level: 4 });

// steps x steps vertices of one cell, [lon0, lat0, lon1, lat1, ...]
// (i outer, j inner, matching vertex(cell, i / (steps - 1), j / (steps - 1)))
const mesh: Float64Array = grid.vertices(164n, 65);

// centers of many cells, same interleaving
const centers: Float64Array = grid.healpixToLonLat(
  new BigUint64Array([0n, 164n]),
);

// and back: one cell id per [lon, lat] pair
const cells: BigUint64Array = grid.lonLatToHealpix(centers);

// entry `row * size + col` is `bitCombine(col, row)` — the layout for
// unshuffling a z-order-flattened size x size chunk into row-major order
const table: BigUint64Array = grid.bitCombineTable(8);
```

A batch is rejected as a whole, naming the offending index
(`cells[2]: ...`, `lonlats[2]: ...`), rather than trapping on one bad element.

### neighbours

Three more bulk methods look up neighbouring cells, across base cell
boundaries like anywhere else. Each returns a flat `BigInt64Array` with a
fixed number of ids per input cell, in the grid's scheme, with `-1n` where
there is no cell:

```typescript
const cells = new BigUint64Array([42n, 164n]);

// 8 per cell, positional: S, SW, W, NW, N, NE, E, SE
const around: BigInt64Array = grid.neighbours(cells);
// 4 per cell: "edge" is SW, NW, NE, SE; "vertex" is S, W, N, E
const sharingAnEdge: BigInt64Array = grid.neighbours(cells, "edge");

// (2 * ring + 1) ** 2 per cell: everything within `ring` cells, itself included
const window: BigInt64Array = grid.kthNeighbourhood(cells, 2);
// 8 * ring per cell: only the cells at distance `ring`
const border: BigInt64Array = grid.kthNeighbours(cells, 2);
```

24 cells of every grid, at the corners where only three base cells meet, have
seven neighbours. `neighbours` keeps its directional positions and marks the
missing one; the k-ring methods are not positional and pad at the end. `ring`
is limited to `nside`.

Note that `healpixToLonLat` and `lonLatToHealpix` are only exact inverses
within one level: on a `zuniq` grid `healpixToLonLat` reads each id at its
embedded level while `lonLatToHealpix` encodes at the grid's level, so a
coarse id in a mixed-level batch comes back refined.

## Low-level scheme functions

The `nested` / `ring` / `zuniq` namespaces expose the per-scheme functions
`Grid` is built on (their `depth` parameter is the same 0-indexed quantity as
the grid's `level`). Reach for them when a single call is all you need;
otherwise prefer `Grid`, which parses the ellipsoid once instead of per call.

```typescript
import init, * as healpixGeo from "healpix-geo";

await init(); // or use a bundler

// parse once, reuse for any number of calls
const ellipsoid = healpixGeo.Ellipsoid.from({
  semi_major_axis: 6378137.0,
  inverse_flattening: 298.257223563,
}); // or Ellipsoid.from(null) for the default sphere

const cellId: bigint = 10n;
const level: number = 2;
const { lon, lat } = healpixGeo.nested.healpixToLonLat(
  cellId,
  level,
  ellipsoid,
);
```
