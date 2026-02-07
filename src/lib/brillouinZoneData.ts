/**
 * Brillouin zone geometry and high-symmetry points for all Bravais lattices.
 *
 * Coordinates follow Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010).
 * All k-point coordinates are in fractional (crystal) coordinates of the
 * reciprocal lattice vectors, i.e., k = k1*b1 + k2*b2 + k3*b3.
 *
 * This is the convention used by Quantum ESPRESSO's crystal_b k-point format.
 */

import { Vec3 } from "./reciprocalLattice";

export interface HighSymmetryPoint {
  label: string;
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
export function getRhombohedralBZ(alpha: number): BrillouinZoneData {
  const alphaRad = (alpha * Math.PI) / 180;

  if (alpha < 90) {
    // hR1 (alpha < 90°)
    const eta = (1 + 4 * Math.cos(alphaRad)) / (2 + 4 * Math.cos(alphaRad));
    const nu = 0.75 - eta / 2;

    return {
      latticeType: "hR1",
      name: "Rhombohedral (α < 90°)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "B", coords: [eta, 0.5, 1 - eta], description: "Edge point" },
        { label: "B₁", coords: [0.5, 1 - eta, eta - 1], description: "Edge point" },
        { label: "F", coords: [0.5, 0.5, 0], description: "Face center" },
        { label: "L", coords: [0.5, 0, 0], description: "Edge center" },
        { label: "L₁", coords: [0, 0, -0.5], description: "Edge center" },
        { label: "P", coords: [eta, nu, nu], description: "On Λ line" },
        { label: "P₁", coords: [1 - nu, 1 - nu, 1 - eta], description: "On P line" },
        { label: "P₂", coords: [nu, nu, eta - 1], description: "On P line" },
        { label: "Q", coords: [1 - nu, nu, 0], description: "On Q line" },
        { label: "X", coords: [nu, 0, -nu], description: "On X line" },
        { label: "Z", coords: [0.5, 0.5, 0.5], description: "Corner" },
      ],
      recommendedPath: [
        ["Γ", "L"], ["L", "B₁"], ["B₁", "B"], ["B", "Z"], ["Z", "Γ"], ["Γ", "X"], ["X", "Q"], ["Q", "F"], ["F", "P₁"], ["P₁", "Z"], ["L", "P"]
      ],
      vertices: [],
      edges: [],
    };
  } else {
    // hR2 (alpha > 90°)
    const eta = 1 / (2 * Math.tan(alphaRad / 2) * Math.tan(alphaRad / 2));
    const nu = 0.75 - eta / 2;

    return {
      latticeType: "hR2",
      name: "Rhombohedral (α > 90°)",
      points: [
        { label: "Γ", coords: [0, 0, 0], description: "Zone center" },
        { label: "F", coords: [0.5, -0.5, 0], description: "Face center" },
        { label: "L", coords: [0.5, 0, 0], description: "Edge center" },
        { label: "P", coords: [1 - nu, -nu, 1 - nu], description: "Vertex" },
        { label: "P₁", coords: [nu, nu - 1, nu - 1], description: "Vertex" },
        { label: "Q", coords: [eta, eta, eta], description: "On Λ line" },
        { label: "Q₁", coords: [1 - eta, -eta, -eta], description: "On Q line" },
        { label: "Z", coords: [0.5, -0.5, 0.5], description: "Corner" },
      ],
      recommendedPath: [
        ["Γ", "P"], ["P", "Z"], ["Z", "Q"], ["Q", "Γ"], ["Γ", "F"], ["F", "P₁"], ["P₁", "Q₁"], ["Q₁", "L"], ["L", "Z"]
      ],
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
  }
): BrillouinZoneData {
  const { a, b, c, alpha, beta } = params;
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
    case "oC":
    case "oI":
    case "oF":
      // All orthorhombic types use similar BZ with different points
      // For now, use primitive orthorhombic as base
      return getOrthorhombicPrimitiveBZ();
    case "hP":
      return getHexagonalBZ();
    case "hR":
      return getRhombohedralBZ(alpha);
    case "mP":
    case "mC":
      return getMonoclinicPrimitiveBZ(b_over_a, c_over_a, beta);
    case "aP":
    default:
      return getTriclinicBZ();
  }
}
