import assert from "node:assert/strict";
import test from "node:test";
import { balancedPoolCounts, bandPathPointCount, savedPointStructure } from "../../src/lib/qePointEstimate";

test("balanced pools divide both ranks per image and distinct points", () => {
  assert.deepEqual(balancedPoolCounts(64, 8, 12), [1, 2, 4]);
  assert.deepEqual(balancedPoolCounts(32, 1, 10), [1, 2]);
  assert.deepEqual(balancedPoolCounts(16, 1, 1), [1]);
  assert.deepEqual(balancedPoolCounts(63, 8, 12), []);
  assert.deepEqual(balancedPoolCounts(64, 0, 12), []);
});

test("saved Cartesian positions use the saved cell, including Bohr units", () => {
  const structure = savedPointStructure({
    cell_parameters: [[4, 0, 0], [0, 6, 0], [0, 0, 8]], cell_units: "bohr",
    position_units: "bohr", atoms: [{ symbol: "Si", position: [1, 3, 6] }],
  });
  assert.ok(structure);
  for (const [i, expected] of [0.25, 0.5, 0.75].entries()) {
    assert.ok(Math.abs(structure.atoms[0].position[i] - expected) < 1e-12);
  }
  assert.ok(Math.abs(structure.lattice[0][0] - 4 * 0.529177210903) < 1e-12);
});

test("an unknown alat scale is not silently estimated", () => {
  assert.equal(savedPointStructure({ cell_parameters: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    cell_units: "alat", position_units: "crystal", atoms: [{ symbol: "Si", position: [0, 0, 0] }] }), null);
});

test("explicit band-path counts include the final endpoint and disconnected samples", () => {
  assert.equal(bandPathPointCount([{ npoints: 20 }, { npoints: 0 }, { npoints: 10 }, { npoints: 0 }]), 32);
  assert.equal(bandPathPointCount([]), 0);
});
