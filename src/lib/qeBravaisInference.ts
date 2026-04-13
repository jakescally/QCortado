import { BravaisLattice, detectBravaisLattice } from "./brillouinZone";
import { Matrix3x3, Vec3 } from "./reciprocalLattice";
import { SymmetryTransformResult } from "./symmetryTransform";
import { CrystalData } from "./types";

const BOHR_TO_ANGSTROM = 0.529177210903;
const ANGSTROM_TO_BOHR = 1 / BOHR_TO_ANGSTROM;
const MONOCLINIC_ANGLE_THRESHOLD_DEG = 1e-2;
const UNIMODULAR_ENTRY_TOL = 2e-2;

type Matrix = [[number, number, number], [number, number, number], [number, number, number]];

export type QeBravaisIbrav =
  | "cubic_p"
  | "cubic_f"
  | "cubic_i"
  | "hexagonal"
  | "trigonal_r"
  | "tetragonal_p"
  | "tetragonal_i"
  | "orthorhombic_p"
  | "orthorhombic_bc"
  | "orthorhombic_fc"
  | "orthorhombic_i"
  | "monoclinic_p"
  | "monoclinic_bc"
  | "triclinic";

interface QeIbravDefinition {
  ibrav: QeBravaisIbrav;
  numeric: number;
}

interface LatticeParameters {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

interface PrimitiveAtom {
  symbol: string;
  position: Vec3;
}

export interface InferredQeBravaisCell {
  ibrav: QeBravaisIbrav;
  celldm: [number, number, number, number, number, number];
  atoms: PrimitiveAtom[];
}

function coerceSpaceGroupNumber(value: unknown): number | null {
  if (value == null) return null;
  const parsed = typeof value === "number"
    ? value
    : Number.parseInt(String(value).trim(), 10);
  if (!Number.isInteger(parsed) || parsed < 1 || parsed > 230) {
    return null;
  }
  return parsed;
}

function extractSpaceGroupNumberFromHM(hm: string | undefined): number | null {
  if (!hm) return null;
  const match = hm.match(/#\s*(\d{1,3})/);
  if (!match) return null;
  return coerceSpaceGroupNumber(match[1]);
}

function bravaisToQeIbrav(bravais: BravaisLattice): QeIbravDefinition {
  switch (bravais) {
    case "cubic-P":
      return { ibrav: "cubic_p", numeric: 1 };
    case "cubic-F":
      return { ibrav: "cubic_f", numeric: 2 };
    case "cubic-I":
      return { ibrav: "cubic_i", numeric: 3 };
    case "hexagonal":
      return { ibrav: "hexagonal", numeric: 4 };
    case "trigonal-R":
      return { ibrav: "trigonal_r", numeric: 5 };
    case "tetragonal-P":
      return { ibrav: "tetragonal_p", numeric: 6 };
    case "tetragonal-I":
      return { ibrav: "tetragonal_i", numeric: 7 };
    case "orthorhombic-P":
      return { ibrav: "orthorhombic_p", numeric: 8 };
    case "orthorhombic-C":
      return { ibrav: "orthorhombic_bc", numeric: 9 };
    case "orthorhombic-F":
      return { ibrav: "orthorhombic_fc", numeric: 10 };
    case "orthorhombic-I":
      return { ibrav: "orthorhombic_i", numeric: 11 };
    case "monoclinic-P":
      return { ibrav: "monoclinic_p", numeric: 12 };
    case "monoclinic-C":
      return { ibrav: "monoclinic_bc", numeric: 13 };
    case "triclinic":
      return { ibrav: "triclinic", numeric: 14 };
    default: {
      const _exhaustive: never = bravais;
      throw new Error(`Unsupported Bravais lattice: ${_exhaustive}`);
    }
  }
}

function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function magnitude(v: Vec3): number {
  return Math.sqrt(dot(v, v));
}

function angleDegrees(a: Vec3, b: Vec3): number {
  const denom = magnitude(a) * magnitude(b);
  if (denom <= 0) return 0;
  const cosine = Math.max(-1, Math.min(1, dot(a, b) / denom));
  return (Math.acos(cosine) * 180) / Math.PI;
}

function latticeParameters(lattice: Matrix3x3): LatticeParameters {
  const [aVec, bVec, cVec] = lattice;
  return {
    a: magnitude(aVec),
    b: magnitude(bVec),
    c: magnitude(cVec),
    alpha: angleDegrees(bVec, cVec),
    beta: angleDegrees(aVec, cVec),
    gamma: angleDegrees(aVec, bVec),
  };
}

function monoclinicToUniqueCAxes(params: LatticeParameters): {
  a: number;
  b: number;
  c: number;
  gamma: number;
} | null {
  const deviations = [
    Math.abs(params.alpha - 90),
    Math.abs(params.beta - 90),
    Math.abs(params.gamma - 90),
  ];
  const maxDeviation = Math.max(...deviations);
  if (maxDeviation < MONOCLINIC_ANGLE_THRESHOLD_DEG) return null;

  // QE positive ibrav=12/13 expects unique axis c (gamma is the non-90 angle).
  if (deviations[2] === maxDeviation) {
    return { a: params.a, b: params.b, c: params.c, gamma: params.gamma };
  }
  if (deviations[1] === maxDeviation) {
    // Common CIF setting (unique axis b) -> swap b and c.
    return { a: params.a, b: params.c, c: params.b, gamma: params.beta };
  }
  // Unique axis a -> rotate axes so gamma carries the non-90 angle.
  return { a: params.b, b: params.c, c: params.a, gamma: params.alpha };
}

function buildCelldm(
  qeIbrav: QeIbravDefinition,
  conventional: Matrix3x3,
  primitive: Matrix3x3,
): [number, number, number, number, number, number] | null {
  const conv = latticeParameters(conventional);
  const prim = latticeParameters(primitive);
  const celldm: [number, number, number, number, number, number] = [0, 0, 0, 0, 0, 0];

  switch (qeIbrav.numeric) {
    case 1:
    case 2:
    case 3:
      celldm[0] = conv.a * ANGSTROM_TO_BOHR;
      return celldm;
    case 4:
    case 6:
    case 7:
      celldm[0] = conv.a * ANGSTROM_TO_BOHR;
      celldm[2] = conv.c / conv.a;
      return celldm;
    case 5: {
      const lengths = [prim.a, prim.b, prim.c];
      const minLen = Math.min(...lengths);
      const maxLen = Math.max(...lengths);
      if (minLen <= 0 || (maxLen - minLen) / minLen > 0.05) return null;
      const cosines = [
        Math.cos((prim.alpha * Math.PI) / 180),
        Math.cos((prim.beta * Math.PI) / 180),
        Math.cos((prim.gamma * Math.PI) / 180),
      ];
      const cosMin = Math.min(...cosines);
      const cosMax = Math.max(...cosines);
      if (cosMax - cosMin > 0.05) return null;
      celldm[0] = ((prim.a + prim.b + prim.c) / 3) * ANGSTROM_TO_BOHR;
      celldm[3] = (cosines[0] + cosines[1] + cosines[2]) / 3;
      return celldm;
    }
    case 8:
    case 9:
    case 10:
    case 11:
      celldm[0] = conv.a * ANGSTROM_TO_BOHR;
      celldm[1] = conv.b / conv.a;
      celldm[2] = conv.c / conv.a;
      return celldm;
    case 12:
    case 13: {
      const mono = monoclinicToUniqueCAxes(conv);
      if (!mono) return null;
      celldm[0] = mono.a * ANGSTROM_TO_BOHR;
      celldm[1] = mono.b / mono.a;
      celldm[2] = mono.c / mono.a;
      celldm[3] = Math.cos((mono.gamma * Math.PI) / 180);
      return celldm;
    }
    case 14:
      celldm[0] = prim.a * ANGSTROM_TO_BOHR;
      celldm[1] = prim.b / prim.a;
      celldm[2] = prim.c / prim.a;
      celldm[3] = Math.cos((prim.alpha * Math.PI) / 180); // cos(b,c)
      celldm[4] = Math.cos((prim.beta * Math.PI) / 180); // cos(a,c)
      celldm[5] = Math.cos((prim.gamma * Math.PI) / 180); // cos(a,b)
      return celldm;
    default:
      return null;
  }
}

function qeBasisFromCelldm(
  qeIbrav: QeIbravDefinition,
  celldm: [number, number, number, number, number, number],
): Matrix3x3 | null {
  const a = celldm[0];
  const bOverA = celldm[1];
  const cOverA = celldm[2];
  const c4 = celldm[3];
  const c5 = celldm[4];
  const c6 = celldm[5];
  const sr3 = Math.sqrt(3);
  const basisBohr: Matrix3x3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

  switch (qeIbrav.numeric) {
    case 1:
      basisBohr[0] = [a, 0, 0];
      basisBohr[1] = [0, a, 0];
      basisBohr[2] = [0, 0, a];
      break;
    case 2: {
      const term = a / 2;
      basisBohr[0] = [-term, 0, term];
      basisBohr[1] = [0, term, term];
      basisBohr[2] = [-term, term, 0];
      break;
    }
    case 3: {
      const term = a / 2;
      basisBohr[0] = [term, term, term];
      basisBohr[1] = [-term, term, term];
      basisBohr[2] = [-term, -term, term];
      break;
    }
    case 4:
      basisBohr[0] = [a, 0, 0];
      basisBohr[1] = [-a / 2, (a * sr3) / 2, 0];
      basisBohr[2] = [0, 0, a * cOverA];
      break;
    case 5: {
      const term1 = Math.sqrt(1 + 2 * c4);
      const term2 = Math.sqrt(1 - c4);
      basisBohr[1] = [0, (Math.sqrt(2) * a * term2) / sr3, (a * term1) / sr3];
      basisBohr[0] = [(a * term2) / Math.sqrt(2), -(a * term2) / Math.sqrt(6), basisBohr[1][2]];
      basisBohr[2] = [-basisBohr[0][0], basisBohr[0][1], basisBohr[1][2]];
      break;
    }
    case 6:
      basisBohr[0] = [a, 0, 0];
      basisBohr[1] = [0, a, 0];
      basisBohr[2] = [0, 0, a * cOverA];
      break;
    case 7: {
      const term = a / 2;
      const z = (a * cOverA) / 2;
      basisBohr[0] = [term, -term, z];
      basisBohr[1] = [term, term, z];
      basisBohr[2] = [-term, -term, z];
      break;
    }
    case 8:
      basisBohr[0] = [a, 0, 0];
      basisBohr[1] = [0, a * bOverA, 0];
      basisBohr[2] = [0, 0, a * cOverA];
      break;
    case 9: {
      const halfA = 0.5 * a;
      basisBohr[0] = [halfA, halfA * bOverA, 0];
      basisBohr[1] = [-halfA, halfA * bOverA, 0];
      basisBohr[2] = [0, 0, a * cOverA];
      break;
    }
    case 10: {
      const halfA = 0.5 * a;
      basisBohr[1] = [halfA, halfA * bOverA, 0];
      basisBohr[0] = [halfA, 0, halfA * cOverA];
      basisBohr[2] = [0, halfA * bOverA, halfA * cOverA];
      break;
    }
    case 11: {
      const halfA = 0.5 * a;
      basisBohr[0] = [halfA, halfA * bOverA, halfA * cOverA];
      basisBohr[1] = [-halfA, halfA * bOverA, halfA * cOverA];
      basisBohr[2] = [-halfA, -halfA * bOverA, halfA * cOverA];
      break;
    }
    case 12: {
      const sinGamma = Math.sqrt(1 - c4 * c4);
      basisBohr[0] = [a, 0, 0];
      basisBohr[1] = [a * bOverA * c4, a * bOverA * sinGamma, 0];
      basisBohr[2] = [0, 0, a * cOverA];
      break;
    }
    case 13: {
      const sinGamma = Math.sqrt(1 - c4 * c4);
      basisBohr[0] = [0.5 * a, 0, -0.5 * a * cOverA];
      basisBohr[1] = [a * bOverA * c4, a * bOverA * sinGamma, 0];
      basisBohr[2] = [0.5 * a, 0, 0.5 * a * cOverA];
      break;
    }
    case 14: {
      const sinGamma = Math.sqrt(1 - c6 * c6);
      const termRaw = 1 + 2 * c4 * c5 * c6 - c4 * c4 - c5 * c5 - c6 * c6;
      if (sinGamma <= 0 || termRaw <= 0) return null;
      const term = Math.sqrt(termRaw / (1 - c6 * c6));
      basisBohr[0] = [a, 0, 0];
      basisBohr[1] = [a * bOverA * c6, a * bOverA * sinGamma, 0];
      basisBohr[2] = [
        a * cOverA * c5,
        a * cOverA * (c4 - c5 * c6) / sinGamma,
        a * cOverA * term,
      ];
      break;
    }
    default:
      return null;
  }

  return basisBohr.map((row) => row.map((value) => value * BOHR_TO_ANGSTROM) as Vec3) as Matrix3x3;
}

function determinant3(matrix: Matrix): number {
  const [r0, r1, r2] = matrix;
  return (
    r0[0] * (r1[1] * r2[2] - r1[2] * r2[1]) -
    r0[1] * (r1[0] * r2[2] - r1[2] * r2[0]) +
    r0[2] * (r1[0] * r2[1] - r1[1] * r2[0])
  );
}

function inverse3(matrix: Matrix): Matrix | null {
  const det = determinant3(matrix);
  if (!Number.isFinite(det) || Math.abs(det) < 1e-10) return null;
  const invDet = 1 / det;
  const [r0, r1, r2] = matrix;
  return [
    [
      (r1[1] * r2[2] - r1[2] * r2[1]) * invDet,
      (r0[2] * r2[1] - r0[1] * r2[2]) * invDet,
      (r0[1] * r1[2] - r0[2] * r1[1]) * invDet,
    ],
    [
      (r1[2] * r2[0] - r1[0] * r2[2]) * invDet,
      (r0[0] * r2[2] - r0[2] * r2[0]) * invDet,
      (r0[2] * r1[0] - r0[0] * r1[2]) * invDet,
    ],
    [
      (r1[0] * r2[1] - r1[1] * r2[0]) * invDet,
      (r0[1] * r2[0] - r0[0] * r2[1]) * invDet,
      (r0[0] * r1[1] - r0[1] * r1[0]) * invDet,
    ],
  ];
}

function basisColumnMatrix(basis: Matrix3x3): Matrix {
  return [
    [basis[0][0], basis[1][0], basis[2][0]],
    [basis[0][1], basis[1][1], basis[2][1]],
    [basis[0][2], basis[1][2], basis[2][2]],
  ];
}

function expressVectorInBasis(vector: Vec3, basis: Matrix3x3): Vec3 | null {
  const inv = inverse3(basisColumnMatrix(basis));
  if (!inv) return null;
  return [
    inv[0][0] * vector[0] + inv[0][1] * vector[1] + inv[0][2] * vector[2],
    inv[1][0] * vector[0] + inv[1][1] * vector[1] + inv[1][2] * vector[2],
    inv[2][0] * vector[0] + inv[2][1] * vector[1] + inv[2][2] * vector[2],
  ];
}

function deriveUnimodularRowTransform(sourceBasis: Matrix3x3, targetBasis: Matrix3x3): Matrix | null {
  const rows: Vec3[] = [];
  for (const targetVector of targetBasis) {
    const coeffs = expressVectorInBasis(targetVector, sourceBasis);
    if (!coeffs) return null;
    rows.push(coeffs);
  }

  const rounded: Matrix = [
    [Math.round(rows[0][0]), Math.round(rows[0][1]), Math.round(rows[0][2])],
    [Math.round(rows[1][0]), Math.round(rows[1][1]), Math.round(rows[1][2])],
    [Math.round(rows[2][0]), Math.round(rows[2][1]), Math.round(rows[2][2])],
  ];

  for (let i = 0; i < 3; i += 1) {
    for (let j = 0; j < 3; j += 1) {
      if (Math.abs(rows[i][j] - rounded[i][j]) > UNIMODULAR_ENTRY_TOL) {
        return null;
      }
    }
  }

  const det = determinant3(rounded);
  if (Math.abs(Math.abs(det) - 1) > 1e-6) {
    return null;
  }

  return rounded;
}

function rowVectorTimesMatrix(vector: Vec3, matrix: Matrix): Vec3 {
  return [
    vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
    vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
    vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
  ];
}

function wrapFractionalCoordinate(value: number): number {
  let wrapped = value - Math.floor(value);
  if (Math.abs(wrapped) < 1e-10 || Math.abs(wrapped - 1) < 1e-10) {
    wrapped = 0;
  }
  return Number(wrapped.toFixed(12));
}

function wrapFractionalVector(vector: Vec3): Vec3 {
  return [
    wrapFractionalCoordinate(vector[0]),
    wrapFractionalCoordinate(vector[1]),
    wrapFractionalCoordinate(vector[2]),
  ];
}

export function inferQeBravaisCellFromCif(
  crystalData: CrystalData,
  symmetryTransform: SymmetryTransformResult | null,
): InferredQeBravaisCell | null {
  if (!symmetryTransform) return null;
  if (symmetryTransform.standardizedPrimitiveAtoms.length === 0) return null;

  const spaceGroup = coerceSpaceGroupNumber(crystalData.space_group_IT_number)
    ?? extractSpaceGroupNumberFromHM(crystalData.space_group_HM)
    ?? coerceSpaceGroupNumber(symmetryTransform.spacegroupNumber);
  if (spaceGroup == null) return null;

  let bravais: BravaisLattice;
  try {
    bravais = detectBravaisLattice(spaceGroup);
  } catch {
    return null;
  }

  const qeIbrav = bravaisToQeIbrav(bravais);
  const sourcePrimitive = symmetryTransform.standardizedPrimitiveLattice;
  const sourceConventional = symmetryTransform.standardizedConventionalLattice;
  const celldm = buildCelldm(qeIbrav, sourceConventional, sourcePrimitive);
  if (!celldm) return null;

  const qePrimitive = qeBasisFromCelldm(qeIbrav, celldm);
  if (!qePrimitive) return null;

  const transform = deriveUnimodularRowTransform(sourcePrimitive, qePrimitive);
  if (!transform) return null;
  const transformInverse = inverse3(transform);
  if (!transformInverse) return null;

  const atoms = symmetryTransform.standardizedPrimitiveAtoms.map((atom) => ({
    symbol: atom.symbol,
    position: wrapFractionalVector(rowVectorTimesMatrix(atom.position, transformInverse)),
  }));

  return {
    ibrav: qeIbrav.ibrav,
    celldm,
    atoms,
  };
}
