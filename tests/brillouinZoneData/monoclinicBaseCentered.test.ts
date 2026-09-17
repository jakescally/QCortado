import assert from "node:assert/strict";
import test from "node:test";
import { getBrillouinZoneData } from "../../src/lib/brillouinZoneData";
import { getBrillouinZoneLattices } from "../../src/lib/brillouinZoneLattice";
import {
  calculateBrillouinZone, conventionalToPrimitive, dot, fractionalToCartesian, realSpaceLatticeVectors,
  reciprocalLatticeVectors, type Matrix3x3, type Vec3,
} from "../../src/lib/reciprocalLattice";
import { createPathCoordinateConverters, resolvePathTransformContext } from "../../src/lib/kPathTransforms";
import {
  transformWien2kKPathForAcceptedStructure,
  transformWien2kKPathForKlistBand,
} from "../../src/lib/wien2kBandsWizard";
import type { CrystalData } from "../../src/lib/types";
import type { SymmetryTransformResult } from "../../src/lib/symmetryTransform";

// Actual cell from Fall 2026 / XPn2 / TaP2, CIF ICSD 648187 (SG 5).
const tap2 = { a: 8.87, b: 3.267, c: 7.497, alpha: 90, beta: 119.4, gamma: 90 };
type Params = typeof tap2;
function crystal(params: Params = tap2): CrystalData {
  return {
    cell_length_a: { value: params.a }, cell_length_b: { value: params.b }, cell_length_c: { value: params.c },
    cell_angle_alpha: { value: params.alpha }, cell_angle_beta: { value: params.beta }, cell_angle_gamma: { value: params.gamma },
    space_group_IT_number: 5, space_group_HM: "C 1 2 1",
    atom_sites: [], symmetry_operations: [], anisotropic_params: [],
  };
}
function lattice(p: Params): Matrix3x3 {
  return realSpaceLatticeVectors(p.a, p.b, p.c, p.alpha, p.beta, p.gamma);
}
function closeVector(actual: Vec3, expected: Vec3) {
  actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-9, `${actual} != ${expected}`));
}

// Independent Wigner-Seitz definition, not the wireframe's hull algorithm:
// k is in the first BZ iff k.G <= |G|²/2 for every reciprocal lattice vector G.
function assertInside(points: { label: string; coords: Vec3 }[], basis: Matrix3x3) {
  const planes: Vec3[] = [];
  for (let h = -4; h <= 4; h++) for (let k = -4; k <= 4; k++) for (let l = -4; l <= 4; l++) {
    if (h || k || l) planes.push(fractionalToCartesian([h, k, l], basis));
  }
  for (const point of points) {
    const cart = fractionalToCartesian(point.coords, basis);
    for (const g of planes) {
      assert.ok(2 * dot(cart, g) <= dot(g, g) + 1e-9, `${point.label} outside BZ: ${point.coords}`);
    }
  }
}

test("TaP2 selects mC1 and all table points lie in the actual input-cell Brillouin zone", () => {
  const data = getBrillouinZoneData("mC", tap2);
  const { primitiveLattice } = getBrillouinZoneLattices(crystal(), "C", true, null);
  assert.equal(data.latticeType, "mC1");
  assert.equal(data.points.length, 16);
  assertInside(data.points, reciprocalLatticeVectors(primitiveLattice));
  // Independently evaluated SC Table 17, then transformed to CIF axes.
  closeVector(data.points.find(p => p.label === "F")!.coords,
    [0.5325537153939529, -0.5325537153939529, -0.1120980598552377]);
  closeVector(data.points.find(p => p.label === "Y")!.coords, [0.5, -0.5, 0]);
  assert.ok(data.recommendedPath.some(([a, b]) => a === "F" && b === "L"));
});

// Reduced conventional cells: c*cos(alpha) < b/2. Non-reduced test
// cells can describe a different BZ variation after lattice reduction.
const branchAlpha = 75;
const branchC = 1.2;
const gammaBoundary = Math.sin(branchAlpha * Math.PI / 180);
const criterionBoundary = gammaBoundary / Math.sqrt(1 - Math.cos(branchAlpha * Math.PI / 180) / branchC);
const branches = [
  { type: "mC1", a: 0.7 },
  { type: "mC2", a: gammaBoundary },
  { type: "mC3", a: 1.13 },
  { type: "mC4", a: criterionBoundary },
  { type: "mC5", a: 1 },
];
for (const branch of branches) {
  for (const uniqueB of [false, true]) for (const obtuse of [false, true]) {
    test(`${branch.type} has valid points and paths, unique-${uniqueB ? "b" : "a"}, ${obtuse ? "obtuse" : "acute"}`, () => {
      const angle = obtuse ? 180 - branchAlpha : branchAlpha;
      const p = uniqueB
        ? { a: 1, b: branch.a, c: branchC, alpha: 90, beta: angle, gamma: 90 }
        : { a: branch.a, b: 1, c: branchC, alpha: angle, beta: 90, gamma: 90 };
      const data = getBrillouinZoneData("mC", p);
      assert.equal(data.latticeType, branch.type);
      assertInside(data.points, reciprocalLatticeVectors(conventionalToPrimitive(lattice(p), "C")));
      for (const [from, to] of data.recommendedPath) {
        assert.ok(data.points.some(p => p.label === from));
        assert.ok(data.points.some(p => p.label === to));
      }
      if (branch.type === "mC3" || branch.type === "mC4") {
        assert.ok(!data.points.some(p => p.label === "L")); // Not in SC Table 18.
      }
      if (branch.type === "mC5") assert.ok(data.recommendedPath.some(([a, b]) => a === "F" && b === "L"));
    });
  }
}

test("mC branch tolerances keep both exact boundaries and their neighboring branches finite", () => {
  for (const [a, expected] of [
    [gammaBoundary - 1e-5, "mC1"], [gammaBoundary + 1e-5, "mC5"],
    [criterionBoundary - 1e-5, "mC5"], [criterionBoundary + 1e-5, "mC3"],
  ] as const) {
    const p = { a, b: 1, c: branchC, alpha: branchAlpha, beta: 90, gamma: 90 };
    const data = getBrillouinZoneData("mC", p);
    assert.equal(data.latticeType, expected);
    assertInside(data.points, reciprocalLatticeVectors(conventionalToPrimitive(lattice(p), "C")));
  }
});

test("TaP2 displayed points and QE/WIEN2k exports have identical Cartesian k-vectors", () => {
  const input = lattice(tap2);
  // Saved TaP2 SCF uses spglib's (a-b)/2, (a+b)/2, c primitive cell,
  // distinct from the viewer's (a+b)/2, (-a+b)/2, c cell.
  const spglibPrimitive: Matrix3x3 = [
    [4.435, -1.6335, 0], [4.435, 1.6335, 0], input[2],
  ];
  const identity: Matrix3x3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const symmetry: SymmetryTransformResult = {
    spacegroupNumber: 12, hallNumber: 63, internationalSymbol: "C2/m", choice: "b1",
    inputLattice: input, standardizedConventionalLattice: input, standardizedPrimitiveLattice: spglibPrimitive,
    standardizedConventionalAtoms: [], standardizedPrimitiveAtoms: [],
    inputToPrimitiveReciprocal: [[0.5, -0.5, 0], [0.5, 0.5, 0], [0, 0, 1]],
    primitiveToInputReciprocal: [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
    primitiveToStandardizedConventionalReciprocal: [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
    standardizedConventionalToPrimitiveReciprocal: [[0.5, -0.5, 0], [0.5, 0.5, 0], [0, 0, 1]],
    transformationMatrix: identity, originShift: [0, 0, 0],
  };
  const c = crystal();
  const { primitiveLattice } = getBrillouinZoneLattices(c, "C", true, symmetry);
  const basis = reciprocalLatticeVectors(primitiveLattice);
  const wireframe = calculateBrillouinZone(basis);
  const oldWireframe = calculateBrillouinZone(reciprocalLatticeVectors(spglibPrimitive));
  assert.equal(wireframe.vertices.length, oldWireframe.vertices.length);
  assert.equal(wireframe.edges.length, oldWireframe.edges.length);
  for (const vertex of wireframe.vertices) {
    assert.ok(oldWireframe.vertices.some(v => v.every((x, i) => Math.abs(x - vertex[i]) < 1e-9)));
  }
  const data = getBrillouinZoneData("mC", tap2);
  assertInside(data.points, basis);
  const path = data.points.map(p => ({ ...p, npoints: 20 }));
  const wienPath = transformWien2kKPathForKlistBand(path, c);
  const acceptedStructurePath = transformWien2kKPathForAcceptedStructure(path, c, symmetry);
  assert.deepEqual(acceptedStructurePath, wienPath);
  for (const backend of [null, symmetry]) {
    const converters = createPathCoordinateConverters(resolvePathTransformContext(c, backend), backend);
    for (const [i, point] of data.points.entries()) {
      const displayed = fractionalToCartesian(point.coords, basis);
      closeVector(fractionalToCartesian(converters.toInputConventionalCoords(point.coords), reciprocalLatticeVectors(input)), displayed);
      closeVector(fractionalToCartesian(wienPath[i].coords, reciprocalLatticeVectors(input)), displayed);
      if (backend) closeVector(fractionalToCartesian(converters.toSymmetryPrimitiveCoords(point.coords), reciprocalLatticeVectors(spglibPrimitive)), displayed);
    }
  }
  // A different backend conventional frame must not change mC's input-based contract.
  const permuted = { ...symmetry, standardizedConventionalLattice: [input[1], input[0], input[2]] as Matrix3x3 };
  assert.deepEqual(getBrillouinZoneLattices(c, "C", true, permuted).primitiveLattice, primitiveLattice);
});
