import assert from "node:assert/strict";
import test from "node:test";
import { sourceScfUsesPrimitiveCell } from "../../src/lib/kPathTransforms";

test("detects primitive SCF cell representations directly", () => {
  assert.equal(
    sourceScfUsesPrimitiveCell({ cell_representation: "primitive_spglib" }, 4),
    true,
  );
  assert.equal(
    sourceScfUsesPrimitiveCell({ cell_representation: "primitive_input" }, 4),
    true,
  );
});

test("treats optimized-source SCF as primitive when atom count shrinks vs conventional", () => {
  const scfParams = {
    cell_representation: "optimized_source",
    source_structure: {
      atoms: [{ symbol: "Cu", position: [0, 0, 0] }],
    },
  };
  assert.equal(sourceScfUsesPrimitiveCell(scfParams, 4), true);
});

test("keeps optimized-source SCF non-primitive when atom count matches conventional", () => {
  const scfParams = {
    cell_representation: "optimized_source",
    source_structure: {
      atoms: [
        { symbol: "Cu", position: [0, 0, 0] },
        { symbol: "Cu", position: [0, 0.5, 0.5] },
        { symbol: "Cu", position: [0.5, 0, 0.5] },
        { symbol: "Cu", position: [0.5, 0.5, 0] },
      ],
    },
  };
  assert.equal(sourceScfUsesPrimitiveCell(scfParams, 4), false);
});

test("does not infer primitive optimized-source mode without conventional atom-count context", () => {
  const scfParams = {
    cell_representation: "optimized_source",
    source_structure: {
      atoms: [{ symbol: "Cu", position: [0, 0, 0] }],
    },
  };
  assert.equal(sourceScfUsesPrimitiveCell(scfParams), false);
});

