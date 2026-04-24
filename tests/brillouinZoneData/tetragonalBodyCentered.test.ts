import assert from "node:assert/strict";
import test from "node:test";
import {
  findHighSymmetryPoint,
  getBrillouinZoneData,
  getHighSymmetryPointId,
} from "../../src/lib/brillouinZoneData";

const tetragonalAngles = {
  alpha: 90,
  beta: 90,
  gamma: 90,
};

function assertVecClose(actual: number[], expected: number[], tolerance = 1e-12) {
  assert.equal(actual.length, expected.length);
  for (let i = 0; i < actual.length; i += 1) {
    assert.ok(
      Math.abs(actual[i] - expected[i]) <= tolerance,
      `component ${i}: expected ${expected[i]}, received ${actual[i]}`,
    );
  }
}

test("body-centered tetragonal tI2 uses SeeK-path labels and coordinates", () => {
  const a = 3.87;
  const c = 12.74;
  const data = getBrillouinZoneData("tI", {
    a,
    b: a,
    c,
    ...tetragonalAngles,
  });

  const eta = (1 + (a * a) / (c * c)) / 4;
  const zeta = (a * a) / (2 * c * c);

  assert.equal(data.latticeType, "tI2");
  assert.deepEqual(data.points.map((point) => point.label), [
    "Γ",
    "M",
    "X",
    "P",
    "N",
    "S₀",
    "S",
    "R",
    "G",
  ]);
  assertVecClose(findHighSymmetryPoint(data, "M")?.coords ?? [], [0.5, 0.5, -0.5]);
  assertVecClose(findHighSymmetryPoint(data, "S_0")?.coords ?? [], [-eta, eta, eta]);
  assertVecClose(findHighSymmetryPoint(data, "S")?.coords ?? [], [eta, 1 - eta, -eta]);
  assertVecClose(findHighSymmetryPoint(data, "R")?.coords ?? [], [-zeta, zeta, 0.5]);
  assertVecClose(findHighSymmetryPoint(data, "G")?.coords ?? [], [0.5, 0.5, -zeta]);

  assert.equal(findHighSymmetryPoint(data, "G")?.label, "G");
  assert.equal(findHighSymmetryPoint(data, "gG")?.label, "Γ");
  assert.equal(getHighSymmetryPointId(findHighSymmetryPoint(data, "G")!), "G_BCT2");
  assert.deepEqual(data.recommendedPath, [
    ["Γ", "X"],
    ["X", "P"],
    ["P", "N"],
    ["N", "Γ"],
    ["Γ", "M"],
    ["M", "S"],
    ["S₀", "Γ"],
    ["X", "R"],
    ["G", "M"],
  ]);
});

test("body-centered tetragonal tI1 uses SeeK-path eta, Z0 label, and path", () => {
  const a = 4;
  const c = 3;
  const data = getBrillouinZoneData("tI", {
    a,
    b: a,
    c,
    ...tetragonalAngles,
  });

  const eta = (1 + (c * c) / (a * a)) / 4;

  assert.equal(data.latticeType, "tI1");
  assertVecClose(findHighSymmetryPoint(data, "Z")?.coords ?? [], [eta, eta, -eta]);
  assertVecClose(findHighSymmetryPoint(data, "Z_0")?.coords ?? [], [-eta, 1 - eta, eta]);
  assert.deepEqual(data.recommendedPath, [
    ["Γ", "X"],
    ["X", "M"],
    ["M", "Γ"],
    ["Γ", "Z"],
    ["Z₀", "M"],
    ["X", "P"],
    ["P", "N"],
    ["N", "Γ"],
  ]);
});
