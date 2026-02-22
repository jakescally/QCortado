/**
 * Brillouin zone geometry and high-symmetry points for all Bravais lattices.
 *
 * Coordinates follow Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010).
 * All k-point coordinates are in fractional (crystal) coordinates of the
 * reciprocal lattice vectors, i.e., k = k1*b1 + k2*b2 + k3*b3.
 *
 * This is the convention used by Quantum ESPRESSO's crystal_b k-point format.
 */

import {
  Vec3,
  realSpaceLatticeVectors,
  reciprocalLatticeVectors,
  conventionalToPrimitive,
  magnitude,
  dot,
} from "./reciprocalLattice";

export type RhombohedralConvention = "sc_primitive" | "bilbao_hex";

export interface BrillouinZoneDataOptions {
  rhombohedralConvention?: RhombohedralConvention;
}

export interface HighSymmetryPoint {
  /**
   * Convention-stable point identifier used to remap paths when labels change
   * across naming schemes (e.g., rhombohedral primitive vs HEX labels).
   */
  id?: string;
  label: string;
  /** Optional alias labels used in other conventions (Bilbao/papers/QE snippets). */
  aliases?: string[];
  /** Fractional coordinates in reciprocal lattice basis [k1, k2, k3] */
  coords: Vec3;
  /** Description of the point location */
  description: string;
}

export interface BrillouinZoneData {
  /** Bravais lattice type */
  latticeType: string;
  /** Human-readable name */
  name: string;
  /** High-symmetry points */
  points: HighSymmetryPoint[];
  /** Recommended path segments (pairs of point labels) */
  recommendedPath: [string, string][];
  /**
   * BZ vertices in fractional reciprocal coordinates.
   * For visualization purposes.
   */
  vertices: Vec3[];
  /**
   * BZ edges as pairs of vertex indices.
   * For visualization purposes.
   */
  edges: [number, number][];
}

const NUMERIC_TOLERANCE = 1e-9;
const SUBSCRIPT_TO_ASCII: Record<string, string> = {
  "₀": "0",
  "₁": "1",
  "₂": "2",
  "₃": "3",
  "₄": "4",
  "₅": "5",
  "₆": "6",
  "₇": "7",
  "₈": "8",
  "₉": "9",
};

function toAsciiLabel(label: string): string {
  return label
    .trim()
    .replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (ch) => SUBSCRIPT_TO_ASCII[ch] ?? ch)
    .replace(/Γ/gi, "G")
    .replace(/[^A-Za-z0-9]/g, "")
    .toUpperCase();
}

function autoAliasesForLabel(label: string): string[] {
  const aliases = new Set<string>();
  if (label === "Γ") {
    aliases.add("G");
    aliases.add("Gamma");
    aliases.add("gG");
  }
  const ascii = label.replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (ch) => SUBSCRIPT_TO_ASCII[ch] ?? ch);
  if (ascii !== label) {
    aliases.add(ascii);
  }
  return [...aliases];
}

export function getHighSymmetryPointId(point: Pick<HighSymmetryPoint, "id" | "label">): string {
  const explicitId = point.id?.trim();
  if (explicitId && explicitId.length > 0) {
    return explicitId;
  }
  return normalizeHighSymmetryLabel(point.label);
}

export function defaultRhombohedralConventionForSetting(setting: "hexagonal" | "rhombohedral" | null | undefined): RhombohedralConvention {
  return setting === "hexagonal" ? "bilbao_hex" : "sc_primitive";
}

function angleDegrees(u: Vec3, v: Vec3): number {
  const denom = magnitude(u) * magnitude(v);
  if (denom <= 0) return 0;
  const cosine = Math.max(-1, Math.min(1, dot(u, v) / denom));
  return (Math.acos(cosine) * 180) / Math.PI;
}

function nearlyEqual(a: number, b: number, tolerance = NUMERIC_TOLERANCE): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tolerance * scale;
}

/**
 * Normalize high-symmetry labels across common variants:
 * - Γ / G / gamma / gG
 * - subscript digits (X₁ <-> X1)
 */
export function normalizeHighSymmetryLabel(label: string): string {
  const ascii = toAsciiLabel(label);
  if (ascii === "GG" || ascii === "GAMMA" || ascii === "G") {
    return "Γ";
  }
  return ascii;
}

/**
 * Resolve a point by label using both canonical labels and aliases.
 */
export function findHighSymmetryPoint(
  data: BrillouinZoneData,
  label: string,
): HighSymmetryPoint | null {
  const needle = normalizeHighSymmetryLabel(label);

  // Pass 1: canonical labels and explicit ids.
  for (const point of data.points) {
    if (normalizeHighSymmetryLabel(point.label) === needle) {
      return point;
    }
    if (normalizeHighSymmetryLabel(getHighSymmetryPointId(point)) === needle) {
      return point;
    }
  }

  // Pass 2: generated aliases.
  for (const point of data.points) {
    for (const alias of autoAliasesForLabel(point.label)) {
      if (normalizeHighSymmetryLabel(alias) === needle) {
        return point;
      }
    }
  }

  // Pass 3: convention-specific aliases.
  for (const point of data.points) {
    for (const alias of point.aliases ?? []) {
      if (normalizeHighSymmetryLabel(alias) === needle) {
        return point;
      }
    }
  }
  return null;
}

export function findHighSymmetryPointById(
  data: BrillouinZoneData,
  id: string,
): HighSymmetryPoint | null {
  const needle = normalizeHighSymmetryLabel(id);
  for (const point of data.points) {
    if (normalizeHighSymmetryLabel(getHighSymmetryPointId(point)) === needle) {
      return point;
    }
  }
  return null;
}

// ============================================================================
// Cubic Lattices
// ============================================================================

/**
 * Simple Cubic (cP) - Primitive cubic
 * BZ is a cube
 */
export function getCubicPrimitiveBZ(): BrillouinZoneData {
  return {
    latticeType: "cP",
    name: "Cubic (Primitive)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "X", coords: [0, 0.5, 0], description: "Face center" },
      { label: "M", coords: [0.5, 0.5, 0], description: "Edge center" },
      { label: "R", coords: [0.5, 0.5, 0.5], description: "Corner" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "M"], ["M", "Γ"], ["Γ", "R"], ["R", "X"], ["M", "R"]
    ],
    vertices: [
      [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
      [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5],
    ],
    edges: [
      [0, 1], [1, 2], [2, 3], [3, 0], // bottom
      [4, 5], [5, 6], [6, 7], [7, 4], // top
      [0, 4], [1, 5], [2, 6], [3, 7], // sides
    ],
  };
}

/**
 * Face-Centered Cubic (cF) - FCC
 * BZ is a truncated octahedron
 */
export function getCubicFaceCenteredBZ(): BrillouinZoneData {
  return {
    latticeType: "cF",
    name: "Cubic (Face-Centered)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "X", coords: [0.5, 0, 0.5], description: "Square face center" },
      { label: "L", coords: [0.5, 0.5, 0.5], description: "Hexagonal face center" },
      { label: "W", coords: [0.5, 0.25, 0.75], description: "Corner (square-hex edge)" },
      { label: "U", coords: [0.625, 0.25, 0.625], description: "Edge center (square-hex)" },
      { label: "K", coords: [0.375, 0.375, 0.75], description: "Edge center (hex-hex)" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "W"], ["W", "K"], ["K", "Γ"], ["Γ", "L"], ["L", "U"], ["U", "W"], ["W", "L"], ["L", "K"], ["U", "X"]
    ],
    // Truncated octahedron vertices (24 vertices)
    vertices: [
      // Square faces perpendicular to x
      [0.5, 0.25, 0.25], [0.5, -0.25, 0.25], [0.5, -0.25, -0.25], [0.5, 0.25, -0.25],
      [-0.5, 0.25, 0.25], [-0.5, -0.25, 0.25], [-0.5, -0.25, -0.25], [-0.5, 0.25, -0.25],
      // Square faces perpendicular to y
      [0.25, 0.5, 0.25], [-0.25, 0.5, 0.25], [-0.25, 0.5, -0.25], [0.25, 0.5, -0.25],
      [0.25, -0.5, 0.25], [-0.25, -0.5, 0.25], [-0.25, -0.5, -0.25], [0.25, -0.5, -0.25],
      // Square faces perpendicular to z
      [0.25, 0.25, 0.5], [-0.25, 0.25, 0.5], [-0.25, -0.25, 0.5], [0.25, -0.25, 0.5],
      [0.25, 0.25, -0.5], [-0.25, 0.25, -0.5], [-0.25, -0.25, -0.5], [0.25, -0.25, -0.5],
    ],
    edges: [
      // This is approximate - full edge connectivity for truncated octahedron
      // For visualization, we'll generate this properly
      [0, 1], [1, 2], [2, 3], [3, 0],
      [4, 5], [5, 6], [6, 7], [7, 4],
      [8, 9], [9, 10], [10, 11], [11, 8],
      [12, 13], [13, 14], [14, 15], [15, 12],
      [16, 17], [17, 18], [18, 19], [19, 16],
      [20, 21], [21, 22], [22, 23], [23, 20],
    ],
  };
}

/**
 * Body-Centered Cubic (cI) - BCC
 * BZ is a rhombic dodecahedron
 */
export function getCubicBodyCenteredBZ(): BrillouinZoneData {
  return {
    latticeType: "cI",
    name: "Cubic (Body-Centered)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "H", coords: [0.5, -0.5, 0.5], description: "Corner vertex" },
      { label: "P", coords: [0.25, 0.25, 0.25], description: "Edge center" },
      { label: "N", coords: [0, 0, 0.5], description: "Face center" },
    ],
    recommendedPath: [
      ["Γ", "H"], ["H", "N"], ["N", "Γ"], ["Γ", "P"], ["P", "H"], ["P", "N"]
    ],
    // Rhombic dodecahedron vertices (14 vertices)
    vertices: [
      [0.5, 0, 0], [-0.5, 0, 0],
      [0, 0.5, 0], [0, -0.5, 0],
      [0, 0, 0.5], [0, 0, -0.5],
      [0.25, 0.25, 0.25], [0.25, 0.25, -0.25], [0.25, -0.25, 0.25], [0.25, -0.25, -0.25],
      [-0.25, 0.25, 0.25], [-0.25, 0.25, -0.25], [-0.25, -0.25, 0.25], [-0.25, -0.25, -0.25],
    ],
    edges: [
      [0, 6], [0, 7], [0, 8], [0, 9],
      [1, 10], [1, 11], [1, 12], [1, 13],
      [2, 6], [2, 7], [2, 10], [2, 11],
      [3, 8], [3, 9], [3, 12], [3, 13],
      [4, 6], [4, 8], [4, 10], [4, 12],
      [5, 7], [5, 9], [5, 11], [5, 13],
    ],
  };
}

// ============================================================================
// Tetragonal Lattices
// ============================================================================

/**
 * Simple Tetragonal (tP)
 * BZ is a rectangular prism
 */
export function getTetragonalPrimitiveBZ(_c_over_a: number): BrillouinZoneData {
  // Note: eta = (1 + c/a^2) / (4*(c/a)^2) is used for body-centered tetragonal

  return {
    latticeType: "tP",
    name: "Tetragonal (Primitive)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "X", coords: [0, 0.5, 0], description: "Face center (ab plane)" },
      { label: "M", coords: [0.5, 0.5, 0], description: "Edge center (ab plane)" },
      { label: "Z", coords: [0, 0, 0.5], description: "Face center (c axis)" },
      { label: "R", coords: [0, 0.5, 0.5], description: "Edge center" },
      { label: "A", coords: [0.5, 0.5, 0.5], description: "Corner" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "M"], ["M", "Γ"], ["Γ", "Z"], ["Z", "R"], ["R", "A"], ["A", "Z"], ["X", "R"], ["M", "A"]
    ],
    vertices: [
      [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
      [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5],
    ],
    edges: [
      [0, 1], [1, 2], [2, 3], [3, 0],
      [4, 5], [5, 6], [6, 7], [7, 4],
      [0, 4], [1, 5], [2, 6], [3, 7],
    ],
  };
}

/**
 * Body-Centered Tetragonal (tI)
 * Two cases depending on c/a ratio
 */
export function getTetragonalBodyCenteredBZ(c_over_a: number): BrillouinZoneData {
  if (c_over_a < 1) {
    // tI1 (c < a)
    const eta = (1 + c_over_a * c_over_a) / (4 * c_over_a * c_over_a);
    return {
      latticeType: "tI1",
      name: "Tetragonal I (c < a)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "M", coords: [-0.5, 0.5, 0.5], description: "Corner" },
        { label: "N", coords: [0, 0.5, 0], description: "Face center" },
        { label: "P", coords: [0.25, 0.25, 0.25], description: "Edge midpoint" },
        { label: "X", coords: [0, 0, 0.5], description: "Face center" },
        { label: "Z", coords: [eta, eta, -eta], description: "On Σ line" },
        { label: "Z₁", coords: [-eta, 1 - eta, eta], description: "On Σ₁ line" },
      ],
      recommendedPath: [
        ["Γ", "X"], ["X", "M"], ["M", "Γ"], ["Γ", "Z"], ["Z", "P"], ["P", "N"], ["N", "Z₁"], ["Z₁", "M"], ["X", "P"]
      ],
      vertices: [], // Complex shape
      edges: [],
    };
  } else {
    // tI2 (c > a)
    const eta = (1 + c_over_a * c_over_a) / (4 * c_over_a * c_over_a);
    const zeta = c_over_a * c_over_a / (2 + 2 * c_over_a * c_over_a);
    return {
      latticeType: "tI2",
      name: "Tetragonal I (c > a)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "N", coords: [0, 0.5, 0], description: "Face center" },
        { label: "P", coords: [0.25, 0.25, 0.25], description: "Vertex" },
        { label: "S", coords: [-eta, eta, eta], description: "On Σ line" },
        { label: "S₁", coords: [eta, 1 - eta, -eta], description: "On Σ₁ line" },
        { label: "X", coords: [0, 0, 0.5], description: "Face center" },
        { label: "Y", coords: [-zeta, zeta, 0.5], description: "Edge point" },
        { label: "Y₁", coords: [0.5, 0.5, -zeta], description: "Edge point" },
        { label: "Z", coords: [0.5, 0.5, -0.5], description: "Corner" },
      ],
      recommendedPath: [
        ["Γ", "X"], ["X", "Y"], ["Y", "S"], ["S", "Γ"], ["Γ", "Z"], ["Z", "S₁"], ["S₁", "N"], ["N", "P"], ["P", "Y₁"], ["Y₁", "Z"], ["X", "P"]
      ],
      vertices: [],
      edges: [],
    };
  }
}

// ============================================================================
// Orthorhombic Lattices
// ============================================================================

/**
 * Simple Orthorhombic (oP)
 */
export function getOrthorhombicPrimitiveBZ(): BrillouinZoneData {
  return {
    latticeType: "oP",
    name: "Orthorhombic (Primitive)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "X", coords: [0.5, 0, 0], description: "Face center" },
      { label: "Y", coords: [0, 0.5, 0], description: "Face center" },
      { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
      { label: "S", coords: [0.5, 0.5, 0], description: "Edge center" },
      { label: "T", coords: [0, 0.5, 0.5], description: "Edge center" },
      { label: "U", coords: [0.5, 0, 0.5], description: "Edge center" },
      { label: "R", coords: [0.5, 0.5, 0.5], description: "Corner" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "S"], ["S", "Y"], ["Y", "Γ"], ["Γ", "Z"], ["Z", "U"], ["U", "R"], ["R", "T"], ["T", "Z"], ["Y", "T"], ["U", "X"], ["S", "R"]
    ],
    vertices: [
      [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
      [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5],
    ],
    edges: [
      [0, 1], [1, 2], [2, 3], [3, 0],
      [4, 5], [5, 6], [6, 7], [7, 4],
      [0, 4], [1, 5], [2, 6], [3, 7],
    ],
  };
}

/**
 * Face-centered orthorhombic (oF)
 * Three branches depending on reciprocal metric.
 */
export function getOrthorhombicFaceCenteredBZ(a: number, b: number, c: number): BrillouinZoneData {
  const invA2 = 1 / (a * a);
  const invB2 = 1 / (b * b);
  const invC2 = 1 / (c * c);
  const metricSplit = invB2 + invC2;

  if (invA2 > metricSplit && !nearlyEqual(invA2, metricSplit)) {
    // oF1
    const eta = (1 + (a * a) / (b * b) + (a * a) / (c * c)) / 4;
    const zeta = (1 + (a * a) / (b * b) - (a * a) / (c * c)) / 4;
    return {
      latticeType: "oF1",
      name: "Orthorhombic (Face-Centered, oF1)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "A", coords: [0.5, 0.5 + zeta, zeta], description: "Face point" },
        { label: "A₁", coords: [0.5, 0.5 - zeta, 1 - zeta], description: "Face point" },
        { label: "L", coords: [0.5, 0.5, 0.5], description: "Zone boundary point" },
        { label: "T", coords: [1, 0.5, 0.5], description: "Edge point" },
        { label: "X", coords: [0, eta, eta], description: "Edge point" },
        { label: "X₁", coords: [1, 1 - eta, 1 - eta], description: "Edge point" },
        { label: "Y", coords: [0.5, 0, 0.5], description: "Face center" },
        { label: "Z", coords: [0.5, 0.5, 0], description: "Face center" },
      ],
      recommendedPath: [
        ["Γ", "Y"], ["Y", "T"], ["T", "Z"], ["Z", "Γ"], ["Γ", "X"], ["X", "A₁"], ["A₁", "Y"],
        ["T", "X₁"],
        ["X", "A"], ["A", "Z"],
        ["L", "Γ"],
      ],
      vertices: [],
      edges: [],
    };
  }

  if (invA2 < metricSplit && !nearlyEqual(invA2, metricSplit)) {
    // oF2
    const eta = (1 + (a * a) / (b * b) - (a * a) / (c * c)) / 4;
    const phi = (1 + (c * c) / (b * b) - (c * c) / (a * a)) / 4;
    const delta = (1 + (b * b) / (a * a) - (b * b) / (c * c)) / 4;
    return {
      latticeType: "oF2",
      name: "Orthorhombic (Face-Centered, oF2)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "C", coords: [0.5, 0.5 - eta, 1 - eta], description: "Face point" },
        { label: "C₁", coords: [0.5, 0.5 + eta, eta], description: "Face point" },
        { label: "D", coords: [0.5 - delta, 0.5, 1 - delta], description: "Edge point" },
        { label: "D₁", coords: [0.5 + delta, 0.5, delta], description: "Edge point" },
        { label: "L", coords: [0.5, 0.5, 0.5], description: "Zone boundary point" },
        { label: "H", coords: [1 - phi, 0.5 - phi, 0.5], description: "Edge point" },
        { label: "H₁", coords: [phi, 0.5 + phi, 0.5], description: "Edge point" },
        { label: "X", coords: [0, 0.5, 0.5], description: "Face center" },
        { label: "Y", coords: [0.5, 0, 0.5], description: "Face center" },
        { label: "Z", coords: [0.5, 0.5, 0], description: "Face center" },
      ],
      recommendedPath: [
        ["Γ", "Y"], ["Y", "C"], ["C", "D"], ["D", "X"], ["X", "Γ"], ["Γ", "Z"], ["Z", "D₁"], ["D₁", "H"], ["H", "C"],
        ["C₁", "Z"],
        ["X", "H₁"],
        ["H", "Y"],
        ["L", "Γ"],
      ],
      vertices: [],
      edges: [],
    };
  }

  // oF3: boundary condition invA2 == invB2 + invC2
  const zeta = (1 + (a * a) / (b * b) - (a * a) / (c * c)) / 4;
  const eta = (1 + (a * a) / (b * b) + (a * a) / (c * c)) / 4;
  return {
    latticeType: "oF3",
    name: "Orthorhombic (Face-Centered, oF3)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "A", coords: [0.5, 0.5 + zeta, zeta], description: "Face point" },
      { label: "A₁", coords: [0.5, 0.5 - zeta, 1 - zeta], description: "Face point" },
      { label: "L", coords: [0.5, 0.5, 0.5], description: "Zone boundary point" },
      { label: "T", coords: [1, 0.5, 0.5], description: "Edge point" },
      { label: "X", coords: [0, eta, eta], description: "Edge point" },
      { label: "Y", coords: [0.5, 0, 0.5], description: "Face center" },
      { label: "Z", coords: [0.5, 0.5, 0], description: "Face center" },
    ],
    recommendedPath: [
      ["Γ", "Y"], ["Y", "T"], ["T", "Z"], ["Z", "Γ"], ["Γ", "X"], ["X", "A₁"], ["A₁", "Y"],
      ["X", "A"], ["A", "Z"],
      ["L", "Γ"],
    ],
    vertices: [],
    edges: [],
  };
}

/**
 * Body-centered orthorhombic (oI).
 */
export function getOrthorhombicBodyCenteredBZ(a: number, b: number, c: number): BrillouinZoneData {
  const zeta = (1 + (a * a) / (c * c)) / 4;
  const eta = (1 + (b * b) / (c * c)) / 4;
  const delta = ((b * b) - (a * a)) / (4 * c * c);
  const mu = ((a * a) + (b * b)) / (4 * c * c);

  return {
    latticeType: "oI",
    name: "Orthorhombic (Body-Centered)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "L", coords: [-mu, mu, 0.5 - delta], description: "Edge point" },
      { label: "L₁", coords: [mu, -mu, 0.5 + delta], description: "Edge point" },
      { label: "L₂", coords: [0.5 - delta, 0.5 + delta, -mu], description: "Edge point" },
      { label: "R", coords: [0, 0.5, 0], description: "Face center" },
      { label: "S", coords: [0.5, 0, 0], description: "Face center" },
      { label: "T", coords: [0, 0, 0.5], description: "Face center" },
      { label: "W", coords: [0.25, 0.25, 0.25], description: "Body point" },
      { label: "X", coords: [-zeta, zeta, zeta], description: "Edge point" },
      { label: "X₁", coords: [zeta, 1 - zeta, -zeta], description: "Edge point" },
      { label: "Y", coords: [eta, -eta, eta], description: "Edge point" },
      { label: "Y₁", coords: [1 - eta, eta, -eta], description: "Edge point" },
      { label: "Z", coords: [0.5, 0.5, -0.5], description: "Corner" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "L"], ["L", "T"], ["T", "W"], ["W", "R"], ["R", "X₁"], ["X₁", "Z"], ["Z", "Γ"], ["Γ", "Y"], ["Y", "S"], ["S", "W"],
      ["L₁", "Y"],
      ["Y₁", "Z"],
    ],
    vertices: [],
    edges: [],
  };
}

/**
 * Base-centered orthorhombic (oC).
 */
export function getOrthorhombicBaseCenteredBZ(a: number, b: number): BrillouinZoneData {
  const zeta = (1 + (a * a) / (b * b)) / 4;
  return {
    latticeType: "oC",
    name: "Orthorhombic (Base-Centered)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "A", coords: [zeta, zeta, 0.5], description: "Face point" },
      { label: "A₁", coords: [-zeta, 1 - zeta, 0.5], description: "Face point" },
      { label: "R", coords: [0, 0.5, 0.5], description: "Edge center" },
      { label: "S", coords: [0, 0.5, 0], description: "Face center" },
      { label: "T", coords: [-0.5, 0.5, 0.5], description: "Corner" },
      { label: "X", coords: [zeta, zeta, 0], description: "Face point" },
      { label: "X₁", coords: [-zeta, 1 - zeta, 0], description: "Face point" },
      { label: "Y", coords: [-0.5, 0.5, 0], description: "Edge center" },
      { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "S"], ["S", "R"], ["R", "A"], ["A", "Z"], ["Z", "Γ"], ["Γ", "Y"], ["Y", "X₁"], ["X₁", "A₁"], ["A₁", "T"], ["T", "Y"],
      ["Z", "T"],
    ],
    vertices: [],
    edges: [],
  };
}

// ============================================================================
// Hexagonal Lattice
// ============================================================================

/**
 * Hexagonal (hP)
 * BZ is a hexagonal prism
 */
export function getHexagonalBZ(): BrillouinZoneData {
  return {
    latticeType: "hP",
    name: "Hexagonal",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "A", coords: [0, 0, 0.5], description: "Top/bottom face center" },
      { label: "H", coords: [1/3, 1/3, 0.5], description: "Top corner" },
      { label: "K", coords: [1/3, 1/3, 0], description: "Hexagonal face corner" },
      { label: "L", coords: [0.5, 0, 0.5], description: "Rectangular face center" },
      { label: "M", coords: [0.5, 0, 0], description: "Edge center" },
    ],
    recommendedPath: [
      ["Γ", "M"], ["M", "K"], ["K", "Γ"], ["Γ", "A"], ["A", "L"], ["L", "H"], ["H", "A"], ["L", "M"], ["K", "H"]
    ],
    // Hexagonal prism vertices (12 vertices)
    vertices: [
      // Bottom hexagon
      [0.5, 0, -0.5], [0.25, 0.433, -0.5], [-0.25, 0.433, -0.5],
      [-0.5, 0, -0.5], [-0.25, -0.433, -0.5], [0.25, -0.433, -0.5],
      // Top hexagon
      [0.5, 0, 0.5], [0.25, 0.433, 0.5], [-0.25, 0.433, 0.5],
      [-0.5, 0, 0.5], [-0.25, -0.433, 0.5], [0.25, -0.433, 0.5],
    ],
    edges: [
      // Bottom hexagon
      [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
      // Top hexagon
      [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 6],
      // Vertical edges
      [0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11],
    ],
  };
}

// ============================================================================
// Rhombohedral (Trigonal) Lattice
// ============================================================================

/**
 * Rhombohedral (hR)
 * Two cases depending on rhombohedral angle
 */
export function getRhombohedralBZ(
  alpha: number,
  convention: RhombohedralConvention = "sc_primitive",
): BrillouinZoneData {
  const alphaRad = (alpha * Math.PI) / 180;

  interface RhombohedralPointTemplate {
    id: string;
    coords: Vec3;
    description: string;
  }

  const hR1PrimitiveLabels: Record<string, string> = {
    gamma: "Γ",
    b: "B",
    b1: "B₁",
    f_face: "F",
    l_edge: "L",
    l1_edge: "L₁",
    p: "P",
    p1: "P₁",
    p2: "P₂",
    q: "Q",
    x: "X",
    z_corner: "Z",
  };
  const hR2PrimitiveLabels: Record<string, string> = {
    gamma: "Γ",
    f_face: "F",
    l_edge: "L",
    p: "P",
    p1: "P₁",
    q: "Q",
    q1: "Q₁",
    z_corner: "Z",
  };
  const hR1BilbaoLabels: Record<string, string> = {
    gamma: "Γ",
    b: "B",
    b1: "B₁",
    f_face: "L",
    l_edge: "F",
    l1_edge: "F₁",
    p: "P",
    p1: "P₁",
    p2: "P₂",
    q: "Q",
    x: "X",
    z_corner: "T",
  };
  const hR2BilbaoLabels: Record<string, string> = {
    gamma: "Γ",
    f_face: "L",
    l_edge: "F",
    p: "P",
    p1: "P₁",
    q: "Q",
    q1: "Q₁",
    z_corner: "T",
  };

  const pathFromIds = (
    ids: [string, string][],
    labels: Record<string, string>,
  ): [string, string][] => ids.map(([fromId, toId]) => [labels[fromId], labels[toId]]);

  const makePoints = (
    templates: RhombohedralPointTemplate[],
    activeLabels: Record<string, string>,
    alternateLabels: Record<string, string>,
  ): HighSymmetryPoint[] => {
    return templates.map((point) => {
      const label = activeLabels[point.id];
      const alternate = alternateLabels[point.id];
      const aliases = alternate && alternate !== label ? [alternate] : undefined;
      return {
        id: point.id,
        label,
        aliases,
        coords: point.coords,
        description: point.description,
      };
    });
  };

  if (alpha < 90) {
    // hR1 (alpha < 90°)
    const eta = (1 + 4 * Math.cos(alphaRad)) / (2 + 4 * Math.cos(alphaRad));
    const nu = 0.75 - eta / 2;
    const hR1Points: RhombohedralPointTemplate[] = [
      { id: "gamma", coords: [0, 0, 0], description: "Zone center" },
      { id: "b", coords: [eta, 0.5, 1 - eta], description: "Edge point" },
      { id: "b1", coords: [0.5, 1 - eta, eta - 1], description: "Edge point" },
      { id: "f_face", coords: [0.5, 0.5, 0], description: "Face center" },
      { id: "l_edge", coords: [0.5, 0, 0], description: "Edge center" },
      { id: "l1_edge", coords: [0, 0, -0.5], description: "Edge center" },
      { id: "p", coords: [eta, nu, nu], description: "On Λ line" },
      { id: "p1", coords: [1 - nu, 1 - nu, 1 - eta], description: "On P line" },
      { id: "p2", coords: [nu, nu, eta - 1], description: "On P line" },
      { id: "q", coords: [1 - nu, nu, 0], description: "On Q line" },
      { id: "x", coords: [nu, 0, -nu], description: "On X line" },
      { id: "z_corner", coords: [0.5, 0.5, 0.5], description: "Corner" },
    ];
    const labels = convention === "bilbao_hex" ? hR1BilbaoLabels : hR1PrimitiveLabels;
    const alternateLabels = convention === "bilbao_hex" ? hR1PrimitiveLabels : hR1BilbaoLabels;
    const pathIds: [string, string][] = [
      ["gamma", "l_edge"],
      ["l_edge", "b1"],
      ["b1", "b"],
      ["b", "z_corner"],
      ["z_corner", "gamma"],
      ["gamma", "x"],
      ["x", "q"],
      ["q", "f_face"],
      ["f_face", "p1"],
      ["p1", "z_corner"],
      ["l_edge", "p"],
    ];

    return {
      latticeType: "hR1",
      name: convention === "bilbao_hex"
        ? "Rhombohedral (α < 90°, Bilbao/CDML HEX labels)"
        : "Rhombohedral (α < 90°)",
      points: makePoints(hR1Points, labels, alternateLabels),
      recommendedPath: pathFromIds(pathIds, labels),
      vertices: [],
      edges: [],
    };
  } else {
    // hR2 (alpha > 90°)
    const eta = 1 / (2 * Math.tan(alphaRad / 2) * Math.tan(alphaRad / 2));
    const nu = 0.75 - eta / 2;
    const hR2Points: RhombohedralPointTemplate[] = [
      { id: "gamma", coords: [0, 0, 0], description: "Zone center" },
      { id: "f_face", coords: [0.5, -0.5, 0], description: "Face center" },
      { id: "l_edge", coords: [0.5, 0, 0], description: "Edge center" },
      { id: "p", coords: [1 - nu, -nu, 1 - nu], description: "Vertex" },
      { id: "p1", coords: [nu, nu - 1, nu - 1], description: "Vertex" },
      { id: "q", coords: [eta, eta, eta], description: "On Λ line" },
      { id: "q1", coords: [1 - eta, -eta, -eta], description: "On Q line" },
      { id: "z_corner", coords: [0.5, -0.5, 0.5], description: "Corner" },
    ];
    const labels = convention === "bilbao_hex" ? hR2BilbaoLabels : hR2PrimitiveLabels;
    const alternateLabels = convention === "bilbao_hex" ? hR2PrimitiveLabels : hR2BilbaoLabels;
    const pathIds: [string, string][] = [
      ["gamma", "p"],
      ["p", "z_corner"],
      ["z_corner", "q"],
      ["q", "gamma"],
      ["gamma", "f_face"],
      ["f_face", "p1"],
      ["p1", "q1"],
      ["q1", "l_edge"],
      ["l_edge", "z_corner"],
    ];

    return {
      latticeType: "hR2",
      name: convention === "bilbao_hex"
        ? "Rhombohedral (α > 90°, Bilbao/CDML HEX labels)"
        : "Rhombohedral (α > 90°)",
      points: makePoints(hR2Points, labels, alternateLabels),
      recommendedPath: pathFromIds(pathIds, labels),
      vertices: [],
      edges: [],
    };
  }
}

// ============================================================================
// Monoclinic Lattices
// ============================================================================

/**
 * Simple Monoclinic (mP)
 * Unique axis b (standard setting)
 */
export function getMonoclinicPrimitiveBZ(b_over_a: number, c_over_a: number, beta: number): BrillouinZoneData {
  const betaRad = (beta * Math.PI) / 180;
  const eta = (1 - b_over_a * Math.cos(betaRad)) / (2 * Math.sin(betaRad) * Math.sin(betaRad));
  const nu = 0.5 - eta * c_over_a * Math.cos(betaRad) / b_over_a;

  return {
    latticeType: "mP",
    name: "Monoclinic (Primitive)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "A", coords: [0.5, 0.5, 0], description: "Edge point" },
      { label: "C", coords: [0, 0.5, 0.5], description: "Edge point" },
      { label: "D", coords: [0.5, 0, 0.5], description: "Edge point" },
      { label: "D₁", coords: [0.5, 0, -0.5], description: "Edge point" },
      { label: "E", coords: [0.5, 0.5, 0.5], description: "Corner" },
      { label: "H", coords: [0, eta, 1 - nu], description: "Face point" },
      { label: "H₁", coords: [0, 1 - eta, nu], description: "Face point" },
      { label: "H₂", coords: [0, eta, -nu], description: "Face point" },
      { label: "M", coords: [0.5, eta, 1 - nu], description: "Edge point" },
      { label: "M₁", coords: [0.5, 1 - eta, nu], description: "Edge point" },
      { label: "M₂", coords: [0.5, eta, -nu], description: "Edge point" },
      { label: "X", coords: [0, 0.5, 0], description: "Face center" },
      { label: "Y", coords: [0, 0, 0.5], description: "Face center" },
      { label: "Y₁", coords: [0, 0, -0.5], description: "Face center" },
      { label: "Z", coords: [0.5, 0, 0], description: "Face center" },
    ],
    recommendedPath: [
      ["Γ", "Z"], ["Z", "D"], ["D", "Y"], ["Y", "Γ"], ["Γ", "A"], ["A", "E"], ["E", "C"], ["C", "Γ"]
    ],
    vertices: [],
    edges: [],
  };
}

/**
 * Base-centered monoclinic (mC).
 * Branch selection and formulas follow Setyawan-Curtarolo conventions.
 */
export function getMonoclinicBaseCenteredBZ(
  a: number,
  b: number,
  c: number,
  alpha: number,
  beta: number,
  gamma: number,
): BrillouinZoneData {
  // Convert conventional mC vectors to primitive vectors so branch logic
  // and point formulas use the correct primitive metric.
  const conventional = realSpaceLatticeVectors(a, b, c, alpha, beta, gamma);
  const primitive = conventionalToPrimitive(conventional, "C");
  const aP = magnitude(primitive[0]);
  const bP = magnitude(primitive[1]);
  const cP = magnitude(primitive[2]);
  const alphaP = angleDegrees(primitive[1], primitive[2]);
  const alphaRad = (alphaP * Math.PI) / 180;
  const sinAlpha = Math.sin(alphaRad);
  const cosAlpha = Math.cos(alphaRad);

  const primitiveReciprocal = reciprocalLatticeVectors(primitive);
  const kGamma = angleDegrees(primitiveReciprocal[0], primitiveReciprocal[1]);
  const criterion =
    (bP * cosAlpha) / cP + ((bP * bP) * (sinAlpha * sinAlpha)) / (aP * aP);

  const isKGamma90 = nearlyEqual(kGamma, 90, 1e-7);
  const isCriterion1 = nearlyEqual(criterion, 1, 1e-7);

  if (kGamma > 90 && !isKGamma90) {
    if (criterion < 1 && !isCriterion1) {
      // mC1
      const zeta = (2 - (bP * cosAlpha) / cP) / (4 * sinAlpha * sinAlpha);
      const eta = 0.5 + (2 * zeta * cP * cosAlpha) / bP;
      const psi = 0.75 - (aP * aP) / (4 * bP * bP * sinAlpha * sinAlpha);
      const phi = psi + (0.75 - psi) * ((bP * cosAlpha) / cP);
      return {
        latticeType: "mC1",
        name: "Monoclinic (Base-Centered, mC1)",
        points: [
          { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
          { label: "N", coords: [0.5, 0, 0], description: "Face center" },
          { label: "N₁", coords: [0, -0.5, 0], description: "Face center" },
          { label: "F", coords: [1 - zeta, 1 - zeta, 1 - eta], description: "Face point" },
          { label: "F₁", coords: [zeta, zeta, eta], description: "Face point" },
          { label: "F₂", coords: [zeta, zeta - 1, eta], description: "Face point" },
          { label: "I", coords: [phi, 1 - phi, 0.5], description: "Edge point" },
          { label: "I₁", coords: [1 - phi, phi - 1, 0.5], description: "Edge point" },
          { label: "L", coords: [0.5, 0.5, 0.5], description: "Corner" },
          { label: "M", coords: [0.5, 0, 0.5], description: "Edge center" },
          { label: "X", coords: [1 - psi, psi - 1, 0], description: "Edge point" },
          { label: "X₁", coords: [psi, 1 - psi, 0], description: "Edge point" },
          { label: "X₂", coords: [psi - 1, -psi, 0], description: "Edge point" },
          { label: "Y", coords: [0.5, 0.5, 0], description: "Face center" },
          { label: "Y₁", coords: [-0.5, -0.5, 0], description: "Face center" },
          { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
        ],
        recommendedPath: [
          ["Γ", "Y"], ["Y", "F"], ["F", "L"], ["L", "I"],
          ["I₁", "Z"], ["Z", "F₁"],
          ["Y", "X₁"],
          ["X", "Γ"], ["Γ", "N"],
          ["M", "Γ"],
        ],
        vertices: [],
        edges: [],
      };
    }

    if (isCriterion1) {
      // mC2
      const zeta = (2 - (bP * cosAlpha) / cP) / (4 * sinAlpha * sinAlpha);
      const eta = 0.5 + (2 * zeta * cP * cosAlpha) / bP;
      const psi = 0.75 - (aP * aP) / (4 * bP * bP * sinAlpha * sinAlpha);
      const phi = psi + (0.75 - psi) * ((bP * cosAlpha) / cP);
      return {
        latticeType: "mC2",
        name: "Monoclinic (Base-Centered, mC2)",
        points: [
          { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
          { label: "N", coords: [0.5, 0, 0], description: "Face center" },
          { label: "N₁", coords: [0, -0.5, 0], description: "Face center" },
          { label: "F", coords: [1 - zeta, 1 - zeta, 1 - eta], description: "Face point" },
          { label: "F₁", coords: [zeta, zeta, eta], description: "Face point" },
          { label: "F₂", coords: [zeta, zeta - 1, eta], description: "Face point" },
          { label: "F₃", coords: [1 - zeta, -zeta, 1 - eta], description: "Face point" },
          { label: "I", coords: [phi, 1 - phi, 0.5], description: "Edge point" },
          { label: "I₁", coords: [1 - phi, phi - 1, 0.5], description: "Edge point" },
          { label: "L", coords: [0.5, 0.5, 0.5], description: "Corner" },
          { label: "M", coords: [0.5, 0, 0.5], description: "Edge center" },
          { label: "X", coords: [1 - psi, psi - 1, 0], description: "Edge point" },
          { label: "Y", coords: [0.5, 0.5, 0], description: "Face center" },
          { label: "Y₁", coords: [-0.5, -0.5, 0], description: "Face center" },
          { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
        ],
        recommendedPath: [
          ["Γ", "Y"], ["Y", "F"], ["F", "L"], ["L", "I"],
          ["I₁", "Z"], ["Z", "F₁"],
          ["N", "Γ"], ["Γ", "M"],
        ],
        vertices: [],
        edges: [],
      };
    }

    // mC3
    const mu = (1 + (bP * bP) / (aP * aP)) / 4;
    const delta = (bP * cP * cosAlpha) / (2 * aP * aP);
    const zeta =
      mu -
      0.25 +
      (1 - (bP * cosAlpha) / cP) / (4 * sinAlpha * sinAlpha);
    const eta = 0.5 + (2 * zeta * cP * cosAlpha) / bP;
    const phi = 1 + zeta - 2 * mu;
    const psi = eta - 2 * delta;
    return {
      latticeType: "mC3",
      name: "Monoclinic (Base-Centered, mC3)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "F", coords: [1 - phi, 1 - phi, 1 - psi], description: "Face point" },
        { label: "F₁", coords: [phi, phi - 1, psi], description: "Face point" },
        { label: "F₂", coords: [1 - phi, -phi, 1 - psi], description: "Face point" },
        { label: "H", coords: [zeta, zeta, eta], description: "Edge point" },
        { label: "H₁", coords: [1 - zeta, -zeta, 1 - eta], description: "Edge point" },
        { label: "H₂", coords: [-zeta, -zeta, 1 - eta], description: "Edge point" },
        { label: "I", coords: [0.5, -0.5, 0.5], description: "Edge center" },
        { label: "I₁", coords: [0.5, 0.5, -0.5], description: "Edge center" },
        { label: "L", coords: [0.5, 0.5, 0.5], description: "Corner" },
        { label: "M", coords: [0.5, 0, 0.5], description: "Edge center" },
        { label: "N", coords: [0.5, 0, 0], description: "Face center" },
        { label: "N₁", coords: [0, -0.5, 0], description: "Face center" },
        { label: "X", coords: [0.5, -0.5, 0], description: "Edge center" },
        { label: "Y", coords: [mu, mu, delta], description: "Edge point" },
        { label: "Y₁", coords: [1 - mu, -mu, -delta], description: "Edge point" },
        { label: "Y₂", coords: [-mu, -mu, -delta], description: "Edge point" },
        { label: "Y₃", coords: [mu, mu - 1, delta], description: "Edge point" },
        { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
      ],
      recommendedPath: [
        ["Γ", "Y"], ["Y", "F"], ["F", "H"], ["H", "Z"], ["Z", "I"], ["I", "F₁"],
        ["H₁", "Y₁"], ["Y₁", "X"], ["X", "Γ"], ["Γ", "N"],
        ["M", "Γ"],
      ],
      vertices: [],
      edges: [],
    };
  }

  if (isKGamma90) {
    // mC4
    const mu = (1 + (bP * bP) / (aP * aP)) / 4;
    const delta = (bP * cP * cosAlpha) / (2 * aP * aP);
    const zeta =
      mu -
      0.25 +
      (1 - (bP * cosAlpha) / cP) / (4 * sinAlpha * sinAlpha);
    const eta = 0.5 + (2 * zeta * cP * cosAlpha) / bP;
    const phi = 1 + zeta - 2 * mu;
    const psi = eta - 2 * delta;
    return {
      latticeType: "mC4",
      name: "Monoclinic (Base-Centered, mC4)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "F", coords: [1 - phi, 1 - phi, 1 - psi], description: "Face point" },
        { label: "F₁", coords: [phi, phi - 1, psi], description: "Face point" },
        { label: "F₂", coords: [1 - phi, -phi, 1 - psi], description: "Face point" },
        { label: "H", coords: [zeta, zeta, eta], description: "Edge point" },
        { label: "H₁", coords: [1 - zeta, -zeta, 1 - eta], description: "Edge point" },
        { label: "H₂", coords: [-zeta, -zeta, 1 - eta], description: "Edge point" },
        { label: "I", coords: [phi, 1 - phi, 0.5], description: "Edge point" },
        { label: "I₁", coords: [1 - phi, phi - 1, 0.5], description: "Edge point" },
        { label: "L", coords: [0.5, 0.5, 0.5], description: "Corner" },
        { label: "M", coords: [0.5, 0, 0.5], description: "Edge center" },
        { label: "N", coords: [0.5, 0, 0], description: "Face center" },
        { label: "N₁", coords: [0, -0.5, 0], description: "Face center" },
        { label: "X", coords: [0.5, -0.5, 0], description: "Edge center" },
        { label: "Y", coords: [mu, mu, delta], description: "Edge point" },
        { label: "Y₁", coords: [1 - mu, -mu, -delta], description: "Edge point" },
        { label: "Y₂", coords: [-mu, -mu, -delta], description: "Edge point" },
        { label: "Y₃", coords: [mu, mu - 1, delta], description: "Edge point" },
        { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
      ],
      recommendedPath: [
        ["Γ", "Y"], ["Y", "F"], ["F", "H"], ["H", "Z"], ["Z", "I"],
        ["F₁", "H₁"], ["H₁", "Y₁"], ["Y₁", "X"], ["X", "Γ"], ["Γ", "N"],
        ["M", "Γ"],
      ],
      vertices: [],
      edges: [],
    };
  }

  // mC5: kGamma < 90
  const zeta =
    ((bP * bP) / (aP * aP) + (1 - (bP * cosAlpha) / cP) / (sinAlpha * sinAlpha)) / 4;
  const eta = 0.5 + (2 * zeta * cP * cosAlpha) / bP;
  const mu = eta / 2 + (bP * bP) / (4 * aP * aP) - (bP * cP * cosAlpha) / (2 * aP * aP);
  const nu = 2 * mu - zeta;
  const rho = 1 - (zeta * aP * aP) / (bP * bP);
  const omega =
    ((4 * nu - 1 - ((bP * bP) * (sinAlpha * sinAlpha)) / (aP * aP)) * cP) /
    (2 * bP * cosAlpha);
  const delta = (zeta * cP * cosAlpha) / bP + omega / 2 - 0.25;

  return {
    latticeType: "mC5",
    name: "Monoclinic (Base-Centered, mC5)",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "F", coords: [nu, nu, omega], description: "Face point" },
      { label: "F₁", coords: [1 - nu, 1 - nu, 1 - omega], description: "Face point" },
      { label: "F₂", coords: [nu, nu - 1, omega], description: "Face point" },
      { label: "H", coords: [zeta, zeta, eta], description: "Edge point" },
      { label: "H₁", coords: [1 - zeta, -zeta, 1 - eta], description: "Edge point" },
      { label: "H₂", coords: [-zeta, -zeta, 1 - eta], description: "Edge point" },
      { label: "I", coords: [rho, 1 - rho, 0.5], description: "Edge point" },
      { label: "I₁", coords: [1 - rho, rho - 1, 0.5], description: "Edge point" },
      { label: "L", coords: [0.5, 0.5, 0.5], description: "Corner" },
      { label: "M", coords: [0.5, 0, 0.5], description: "Edge center" },
      { label: "N", coords: [0.5, 0, 0], description: "Face center" },
      { label: "N₁", coords: [0, -0.5, 0], description: "Face center" },
      { label: "X", coords: [0.5, -0.5, 0], description: "Edge center" },
      { label: "Y", coords: [mu, mu, delta], description: "Edge point" },
      { label: "Y₁", coords: [1 - mu, -mu, -delta], description: "Edge point" },
      { label: "Y₂", coords: [-mu, -mu, -delta], description: "Edge point" },
      { label: "Y₃", coords: [mu, mu - 1, delta], description: "Edge point" },
      { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
    ],
    recommendedPath: [
      ["Γ", "Y"], ["Y", "F"], ["F", "H"], ["H", "Z"], ["Z", "I"], ["I", "F₁"],
      ["H₁", "Y₁"], ["Y₁", "X"], ["X", "Γ"], ["Γ", "N"],
      ["M", "Γ"],
    ],
    vertices: [],
    edges: [],
  };
}

// ============================================================================
// Triclinic Lattice
// ============================================================================

/**
 * Triclinic (aP)
 * Most general case - BZ shape depends on all 6 lattice parameters
 */
export function getTriclinicBZ(): BrillouinZoneData {
  return {
    latticeType: "aP",
    name: "Triclinic",
    points: [
      { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
      { label: "L", coords: [0.5, 0.5, 0], description: "Face center" },
      { label: "M", coords: [0, 0.5, 0.5], description: "Face center" },
      { label: "N", coords: [0.5, 0, 0.5], description: "Face center" },
      { label: "R", coords: [0.5, 0.5, 0.5], description: "Corner" },
      { label: "X", coords: [0.5, 0, 0], description: "Face center" },
      { label: "Y", coords: [0, 0.5, 0], description: "Face center" },
      { label: "Z", coords: [0, 0, 0.5], description: "Face center" },
    ],
    recommendedPath: [
      ["Γ", "X"], ["X", "L"], ["L", "Y"], ["Y", "Γ"], ["Γ", "Z"], ["Z", "M"], ["M", "N"], ["N", "R"], ["R", "Z"]
    ],
    vertices: [
      [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
      [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5],
    ],
    edges: [
      [0, 1], [1, 2], [2, 3], [3, 0],
      [4, 5], [5, 6], [6, 7], [7, 4],
      [0, 4], [1, 5], [2, 6], [3, 7],
    ],
  };
}

// ============================================================================
// Helper to get BZ data for a Bravais lattice type
// ============================================================================

export type BravaisLatticeType =
  | "cP" | "cF" | "cI"  // Cubic
  | "tP" | "tI"         // Tetragonal
  | "oP" | "oI" | "oF" | "oC"  // Orthorhombic
  | "hP" | "hR"         // Hexagonal, Rhombohedral
  | "mP" | "mC"         // Monoclinic
  | "aP";               // Triclinic

/**
 * Get Brillouin zone data for a given Bravais lattice type and parameters.
 */
export function getBrillouinZoneData(
  latticeType: BravaisLatticeType,
  params: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  },
  options?: BrillouinZoneDataOptions,
): BrillouinZoneData {
  const { a, b, c, alpha, beta, gamma } = params;
  const c_over_a = c / a;
  const b_over_a = b / a;

  switch (latticeType) {
    case "cP":
      return getCubicPrimitiveBZ();
    case "cF":
      return getCubicFaceCenteredBZ();
    case "cI":
      return getCubicBodyCenteredBZ();
    case "tP":
      return getTetragonalPrimitiveBZ(c_over_a);
    case "tI":
      return getTetragonalBodyCenteredBZ(c_over_a);
    case "oP":
      return getOrthorhombicPrimitiveBZ();
    case "oC":
      return getOrthorhombicBaseCenteredBZ(a, b);
    case "oI":
      return getOrthorhombicBodyCenteredBZ(a, b, c);
    case "oF":
      return getOrthorhombicFaceCenteredBZ(a, b, c);
    case "hP":
      return getHexagonalBZ();
    case "hR":
      return getRhombohedralBZ(alpha, options?.rhombohedralConvention);
    case "mP":
      return getMonoclinicPrimitiveBZ(b_over_a, c_over_a, beta);
    case "mC":
      return getMonoclinicBaseCenteredBZ(a, b, c, alpha, beta, gamma);
    case "aP":
    default:
      return getTriclinicBZ();
  }
}
