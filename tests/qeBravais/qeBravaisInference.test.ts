import assert from "node:assert/strict";
import test from "node:test";
import { inferQeBravaisCellFromCif, transformKPointToQeBasis } from "../../src/lib/engines/qe/bravaisInference";
import { Matrix3x3, Vec3, reciprocalLatticeVectors, fractionalToCartesian } from "../../src/lib/reciprocalLattice";
import { SymmetryTransformResult } from "../../src/lib/symmetryTransform";
import { CrystalData } from "../../src/lib/types";

const BOHR_TO_ANGSTROM = 0.529177210903;

function identityMatrix(): Matrix3x3 {
  return [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
  ];
}

function createCrystalData(overrides: Partial<CrystalData> = {}): CrystalData {
  return {
    cell_length_a: { value: 5 },
    cell_length_b: { value: 5 },
    cell_length_c: { value: 5 },
    cell_angle_alpha: { value: 90 },
    cell_angle_beta: { value: 90 },
    cell_angle_gamma: { value: 90 },
    atom_sites: [{
      label: "X1",
      type_symbol: "X",
      fract_x: 0,
      fract_y: 0,
      fract_z: 0,
      occupancy: 1,
    }],
    symmetry_operations: [],
    anisotropic_params: [],
    ...overrides,
  };
}

function createSymmetryResult(
  spacegroupNumber: number,
  conventional: Matrix3x3,
  primitive: Matrix3x3,
  primitiveAtoms: { symbol: string; position: Vec3 }[],
): SymmetryTransformResult {
  const atoms = primitiveAtoms.map((atom, index) => ({
    symbol: atom.symbol,
    position: atom.position,
    typeIndex: index,
  }));
  return {
    spacegroupNumber,
    hallNumber: 0,
    internationalSymbol: "",
    choice: "",
    inputLattice: conventional,
    standardizedConventionalLattice: conventional,
    standardizedPrimitiveLattice: primitive,
    standardizedConventionalAtoms: atoms,
    standardizedPrimitiveAtoms: atoms,
    primitiveToInputReciprocal: identityMatrix(),
    inputToPrimitiveReciprocal: identityMatrix(),
    primitiveToStandardizedConventionalReciprocal: identityMatrix(),
    standardizedConventionalToPrimitiveReciprocal: identityMatrix(),
    transformationMatrix: identityMatrix(),
    originShift: [0, 0, 0],
  };
}

function almostEqual(actual: number, expected: number, tolerance = 1e-9): void {
  assert.ok(
    Math.abs(actual - expected) <= tolerance,
    `Expected ${actual} to be within ${tolerance} of ${expected}`,
  );
}

test("infers cubic primitive ibrav from CIF space group", () => {
  const crystalData = createCrystalData({
    space_group_IT_number: 221,
    space_group_HM: "Pm-3m",
  });
  const lattice: Matrix3x3 = [
    [5, 0, 0],
    [0, 5, 0],
    [0, 0, 5],
  ];
  const symmetry = createSymmetryResult(221, lattice, lattice, [
    { symbol: "X", position: [0.125, 0.25, 0.5] },
  ]);

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "cubic_p");
  almostEqual(inferred.celldm[0], 5 / BOHR_TO_ANGSTROM, 1e-8);
  assert.deepEqual(inferred.atoms[0].position, [0.125, 0.25, 0.5]);
});

test("infers fcc ibrav with primitive basis from space group 225", () => {
  const a = 6;
  const conventional: Matrix3x3 = [
    [a, 0, 0],
    [0, a, 0],
    [0, 0, a],
  ];
  const primitive: Matrix3x3 = [
    [-a / 2, 0, a / 2],
    [0, a / 2, a / 2],
    [-a / 2, a / 2, 0],
  ];
  const symmetry = createSymmetryResult(225, conventional, primitive, [
    { symbol: "X", position: [0, 0, 0] },
  ]);
  const crystalData = createCrystalData({
    cell_length_a: { value: a },
    cell_length_b: { value: a },
    cell_length_c: { value: a },
    space_group_IT_number: 225,
    space_group_HM: "Fm-3m",
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "cubic_f");
  almostEqual(inferred.celldm[0], a / BOHR_TO_ANGSTROM, 1e-8);
  assert.deepEqual(inferred.atoms[0].position, [0, 0, 0]);
});

test("returns null when primitive basis is incompatible with inferred ibrav", () => {
  const crystalData = createCrystalData({
    cell_length_a: { value: 4 },
    cell_length_b: { value: 5 },
    cell_length_c: { value: 6 },
    space_group_IT_number: 221,
    space_group_HM: "Pm-3m",
  });
  const lattice: Matrix3x3 = [
    [4, 0, 0],
    [0, 5, 0],
    [0, 0, 6],
  ];
  const symmetry = createSymmetryResult(221, lattice, lattice, [
    { symbol: "X", position: [0, 0, 0] },
  ]);

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.equal(inferred, null);
});

test("infers triclinic ibrav and celldm values from low-symmetry CIF", () => {
  const conventional: Matrix3x3 = [
    [4, 0, 0],
    [1.0, 5, 0],
    [0.8, 0.7, 6],
  ];
  const symmetry = createSymmetryResult(2, conventional, conventional, [
    { symbol: "X", position: [0.1, 0.2, 0.3] },
  ]);
  const crystalData = createCrystalData({
    cell_length_a: { value: 4 },
    cell_length_b: { value: Math.sqrt(26) },
    cell_length_c: { value: Math.sqrt(37.13) },
    cell_angle_alpha: { value: 81 },
    cell_angle_beta: { value: 82 },
    cell_angle_gamma: { value: 79 },
    space_group_IT_number: 2,
    space_group_HM: "P-1",
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "triclinic");
  assert.ok(inferred.celldm[0] > 0);
  assert.ok(inferred.celldm[1] > 0);
  assert.ok(inferred.celldm[2] > 0);
  assert.ok(Math.abs(inferred.celldm[3]) < 1);
  assert.ok(Math.abs(inferred.celldm[4]) < 1);
  assert.ok(Math.abs(inferred.celldm[5]) < 1);
});

test("falls back to symmetry space group when CIF metadata is missing", () => {
  const lattice: Matrix3x3 = [
    [5, 0, 0],
    [0, 5, 0],
    [0, 0, 5],
  ];
  const symmetry = createSymmetryResult(221, lattice, lattice, [
    { symbol: "X", position: [0, 0, 0] },
  ]);
  const crystalData = createCrystalData({
    space_group_IT_number: undefined,
    space_group_HM: undefined,
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "cubic_p");
});

test("infers tetragonal body-centered ibrav from space group 139", () => {
  const a = 4.5;
  const c = 7.2;
  const conventional: Matrix3x3 = [
    [a, 0, 0],
    [0, a, 0],
    [0, 0, c],
  ];
  const primitive: Matrix3x3 = [
    [a / 2, -a / 2, c / 2],
    [a / 2, a / 2, c / 2],
    [-a / 2, -a / 2, c / 2],
  ];
  const symmetry = createSymmetryResult(139, conventional, primitive, [
    { symbol: "X", position: [0.2, 0.3, 0.4] },
  ]);
  const crystalData = createCrystalData({
    cell_length_a: { value: a },
    cell_length_b: { value: a },
    cell_length_c: { value: c },
    space_group_IT_number: 139,
    space_group_HM: "I4/mmm",
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "tetragonal_i");
  almostEqual(inferred.celldm[0], a / BOHR_TO_ANGSTROM, 1e-8);
  almostEqual(inferred.celldm[2], c / a, 1e-12);
});

test("infers orthorhombic face-centered ibrav from space group 69", () => {
  const a = 4;
  const b = 5;
  const c = 6;
  const conventional: Matrix3x3 = [
    [a, 0, 0],
    [0, b, 0],
    [0, 0, c],
  ];
  const primitive: Matrix3x3 = [
    [a / 2, 0, c / 2],
    [a / 2, b / 2, 0],
    [0, b / 2, c / 2],
  ];
  const symmetry = createSymmetryResult(69, conventional, primitive, [
    { symbol: "X", position: [0, 0, 0] },
  ]);
  const crystalData = createCrystalData({
    cell_length_a: { value: a },
    cell_length_b: { value: b },
    cell_length_c: { value: c },
    space_group_IT_number: 69,
    space_group_HM: "Fmmm",
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "orthorhombic_fc");
  almostEqual(inferred.celldm[1], b / a, 1e-12);
  almostEqual(inferred.celldm[2], c / a, 1e-12);
});

test("infers trigonal-r ibrav from rhombohedral space group", () => {
  const a = 5.1;
  const cosAlpha = 0.2;
  const term1 = Math.sqrt(1 + 2 * cosAlpha);
  const term2 = Math.sqrt(1 - cosAlpha);
  const primitive: Matrix3x3 = [
    [(a * term2) / Math.sqrt(2), -(a * term2) / Math.sqrt(6), (a * term1) / Math.sqrt(3)],
    [0, (Math.sqrt(2) * a * term2) / Math.sqrt(3), (a * term1) / Math.sqrt(3)],
    [-(a * term2) / Math.sqrt(2), -(a * term2) / Math.sqrt(6), (a * term1) / Math.sqrt(3)],
  ];
  const symmetry = createSymmetryResult(166, primitive, primitive, [
    { symbol: "X", position: [0.11, 0.22, 0.33] },
  ]);
  const crystalData = createCrystalData({
    space_group_IT_number: 166,
    space_group_HM: "R-3m",
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "trigonal_r");
  almostEqual(inferred.celldm[3], cosAlpha, 1e-10);
});

test("handles common monoclinic unique-b conventional setting via axis remapping", () => {
  const a = 5;
  const b = 6;
  const c = 7;
  const betaDeg = 110;
  const betaRad = (betaDeg * Math.PI) / 180;
  const conventional: Matrix3x3 = [
    [a, 0, 0],
    [0, b, 0],
    [c * Math.cos(betaRad), 0, c * Math.sin(betaRad)],
  ];
  // QE positive monoclinic (ibrav=12) expects unique axis c:
  // use remapped parameters a'=a, b'=c, c'=b, gamma'=beta.
  const primitive: Matrix3x3 = [
    [a, 0, 0],
    [c * Math.cos(betaRad), c * Math.sin(betaRad), 0],
    [0, 0, b],
  ];
  const symmetry = createSymmetryResult(14, conventional, primitive, [
    { symbol: "X", position: [0.2, 0.4, 0.6] },
  ]);
  const crystalData = createCrystalData({
    cell_length_a: { value: a },
    cell_length_b: { value: b },
    cell_length_c: { value: c },
    cell_angle_alpha: { value: 90 },
    cell_angle_beta: { value: betaDeg },
    cell_angle_gamma: { value: 90 },
    space_group_IT_number: 14,
    space_group_HM: "P21/c",
  });

  const inferred = inferQeBravaisCellFromCif(crystalData, symmetry);

  assert.ok(inferred);
  assert.equal(inferred.ibrav, "monoclinic_p");
  almostEqual(inferred.celldm[1], c / a, 1e-10);
  almostEqual(inferred.celldm[2], b / a, 1e-10);
});

// Independent direct-vector definitions from QE's latgen (positive ibrav),
// compared against spglib's primitive cells in a fixed Cartesian frame.
const centeredQeCases = [
  { name: "fcc", sg: 225, a: 4, b: 4, c: 4, source: [[0, 1, 1], [1, 0, 1], [1, 1, 0]], target: [[-1, 0, 1], [0, 1, 1], [-1, 1, 0]] },
  { name: "bcc", sg: 229, a: 4, b: 4, c: 4, source: [[-1, 1, 1], [1, -1, 1], [1, 1, -1]], target: [[1, 1, 1], [-1, 1, 1], [-1, -1, 1]] },
  { name: "tI1", sg: 139, a: 4, b: 4, c: 2.5, source: [[-1, 1, 1], [1, -1, 1], [1, 1, -1]], target: [[1, -1, 1], [1, 1, 1], [-1, -1, 1]] },
  { name: "Sr2RuO4 tI2", sg: 139, a: 3.871, b: 3.871, c: 12.702, source: [[-1, 1, 1], [1, -1, 1], [1, 1, -1]], target: [[1, -1, 1], [1, 1, 1], [-1, -1, 1]] },
  { name: "oC", sg: 65, a: 3, b: 4, c: 5, source: [[1, -1, 0], [1, 1, 0], [0, 0, 2]], target: [[1, 1, 0], [-1, 1, 0], [0, 0, 2]] },
  { name: "oI", sg: 71, a: 3, b: 4, c: 5, source: [[-1, 1, 1], [1, -1, 1], [1, 1, -1]], target: [[1, 1, 1], [-1, 1, 1], [-1, -1, 1]] },
  { name: "oF", sg: 69, a: 3, b: 4, c: 5, source: [[0, 1, 1], [1, 0, 1], [1, 1, 0]], target: [[1, 0, 1], [1, 1, 0], [0, 1, 1]] },
];
for (const fixture of centeredQeCases) {
  test(`${fixture.name}: exported k-points preserve Cartesian positions throughout segments`, () => {
    const lengths = [fixture.a, fixture.b, fixture.c];
    const basis = (rows: number[][]) => rows.map(row => row.map((v, i) => v * lengths[i] / 2)) as Matrix3x3;
    const conventional: Matrix3x3 = [[fixture.a, 0, 0], [0, fixture.b, 0], [0, 0, fixture.c]];
    const source = basis(fixture.source);
    const target = basis(fixture.target);
    const crystal = createCrystalData({ space_group_IT_number: fixture.sg });
    const symmetry = createSymmetryResult(fixture.sg, conventional, source, [{ symbol: "Si", position: [0, 0, 0] }]);
    const cell = inferQeBravaisCellFromCif(crystal, symmetry);
    assert.ok(cell);
    const sourceReciprocal = reciprocalLatticeVectors(source);
    const targetReciprocal = reciprocalLatticeVectors(target);
    const endpoints: Vec3[] = [[0, 0, 0], [0.5, 0.5, -0.5], [0, 0, 0.5], [0.25, 0.25, 0.25], [-0.23, 0.42, 0.71]];
    for (const endpoint of endpoints) {
      // Include interior points: equivalent endpoints alone cannot certify a path.
      for (const t of [0, 0.17, 0.5, 0.83, 1]) {
        const q = endpoint.map(v => v * t) as Vec3;
        const expected = fractionalToCartesian(q, sourceReciprocal);
        const actual = fractionalToCartesian(transformKPointToQeBasis(q, cell), targetReciprocal);
        actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-10, `${fixture.name}: ${q}, ${actual} != ${expected}`));
      }
    }
    if (fixture.name === "Sr2RuO4 tI2") {
      assert.deepEqual(transformKPointToQeBasis([0.5, 0.5, -0.5], cell), [0.5, 0.5, 0.5]);
    }
  });
}

test("explicit spglib cell needs no additional QE basis conversion", () => {
  const q: Vec3 = [-0.4, 0.2, 0.3];
  assert.deepEqual(transformKPointToQeBasis(q, null), q);
});
