import { readFileSync } from "node:fs";
import { describe, expect, test } from "vitest";

// The generated `.d.ts` is the API surface TypeScript consumers actually see.
// `JsValue` parameters default to `any` there, which silently drops the input
// contract, so the hand-written `typescript_custom_section` types are pinned
// here.
const declarations = readFileSync(
  new URL("../pkg/index.d.ts", import.meta.url),
  "utf8",
);

describe("generated index.d.ts", () => {
  test("declares the accepted ellipsoid input shapes", () => {
    expect(declarations).to.include("export type EllipsoidInput");
    expect(declarations).to.include("{ radius: number }");
    expect(declarations).to.include(
      "{ semi_major_axis: number; inverse_flattening: number }",
    );
    expect(declarations).to.include(
      "{ semi_major_axis: number; semi_minor_axis: number }",
    );
  });

  test("ellipsoid entry points are typed, not `any`", () => {
    expect(declarations).to.include(
      "parseEllipsoid(obj: EllipsoidInput | null | undefined)",
    );
    expect(declarations).to.include(
      "from(obj: EllipsoidInput | null | undefined)",
    );
    expect(declarations).to.not.match(/parseEllipsoid\(obj: any\)/);
  });

  test("declares the createGrid options shape", () => {
    expect(declarations).to.include("export type GridOptions");
    expect(declarations).to.include('scheme: "nested" | "ring" | "zuniq"');
    expect(declarations).to.include("createGrid(options: GridOptions)");
    expect(declarations).to.not.match(/createGrid\(options: any\)/);
  });
});
