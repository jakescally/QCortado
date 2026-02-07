/**
 * Reciprocal lattice calculations for Brillouin zone visualization.
 *
 * Conventions follow Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010).
 * All reciprocal space coordinates are in units of 2π/a (crystallographic convention).
 */

export type Vec3 = [number, number, number];
export type Matrix3x3 = [Vec3, Vec3, Vec3];

/**
 * Compute the cross product of two 3D vectors.
 */
export function cross(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

/**
 * Compute the dot product of two 3D vectors.
 */
export function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Scale a vector by a scalar.
 */
export function scale(v: Vec3, s: number): Vec3 {
  return [v[0] * s, v[1] * s, v[2] * s];
}

/**
 * Add two vectors.
 */
export function add(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

/**
 * Compute the magnitude of a vector.
 */
export function magnitude(v: Vec3): number {
  return Math.sqrt(dot(v, v));
}

/**
 * Normalize a vector.
 */
export function normalize(v: Vec3): Vec3 {
  const mag = magnitude(v);
  if (mag === 0) return [0, 0, 0];
  return scale(v, 1 / mag);
}

/**
 * Compute real-space lattice vectors from crystallographic parameters.
 * Uses the standard crystallographic convention where:
 * - a is along x
 * - b is in the xy plane
 * - c is determined by angles
 *
 * @param a - lattice parameter a (Angstrom)
 * @param b - lattice parameter b (Angstrom)
 * @param c - lattice parameter c (Angstrom)
 * @param alpha - angle between b and c (degrees)
 * @param beta - angle between a and c (degrees)
 * @param gamma - angle between a and b (degrees)
 * @returns [a_vec, b_vec, c_vec] - real space lattice vectors
 */
export function realSpaceLatticeVectors(
  a: number,
  b: number,
  c: number,
  alpha: number,
  beta: number,
  gamma: number
): Matrix3x3 {
  // Convert angles to radians
  const alphaRad = (alpha * Math.PI) / 180;
  const betaRad = (beta * Math.PI) / 180;
  const gammaRad = (gamma * Math.PI) / 180;

  // a vector along x
  const a_vec: Vec3 = [a, 0, 0];

  // b vector in xy plane
  const b_vec: Vec3 = [
    b * Math.cos(gammaRad),
    b * Math.sin(gammaRad),
    0,
  ];

  // c vector - general direction
  const cx = c * Math.cos(betaRad);
  const cy = c * (Math.cos(alphaRad) - Math.cos(betaRad) * Math.cos(gammaRad)) / Math.sin(gammaRad);
  const cz = Math.sqrt(c * c - cx * cx - cy * cy);
  const c_vec: Vec3 = [cx, cy, cz];

  return [a_vec, b_vec, c_vec];
}

/**
 * Compute the volume of the unit cell (scalar triple product).
 */
export function cellVolume(lattice: Matrix3x3): number {
  const [a, b, c] = lattice;
  return Math.abs(dot(a, cross(b, c)));
}

/**
 * Compute reciprocal lattice vectors from real-space lattice vectors.
 *
 * The reciprocal lattice vectors are defined as:
 * b1 = 2π (a2 × a3) / V
 * b2 = 2π (a3 × a1) / V
 * b3 = 2π (a1 × a2) / V
 *
 * where V = a1 · (a2 × a3) is the cell volume.
 *
 * Note: We include the 2π factor (physicist convention).
 *
 * @param realLattice - [a1, a2, a3] real space lattice vectors
 * @returns [b1, b2, b3] reciprocal lattice vectors
 */
export function reciprocalLatticeVectors(realLattice: Matrix3x3): Matrix3x3 {
  const [a1, a2, a3] = realLattice;
  const V = dot(a1, cross(a2, a3));
  const twoPiOverV = (2 * Math.PI) / V;

  const b1 = scale(cross(a2, a3), twoPiOverV);
  const b2 = scale(cross(a3, a1), twoPiOverV);
  const b3 = scale(cross(a1, a2), twoPiOverV);

  return [b1, b2, b3];
}

/**
 * Transform a point from fractional (crystal) coordinates to Cartesian coordinates.
 *
 * @param fractional - point in fractional coordinates [k1, k2, k3]
 * @param basis - basis vectors (either real or reciprocal lattice)
 * @returns Cartesian coordinates
 */
export function fractionalToCartesian(fractional: Vec3, basis: Matrix3x3): Vec3 {
  const [b1, b2, b3] = basis;
  return [
    fractional[0] * b1[0] + fractional[1] * b2[0] + fractional[2] * b3[0],
    fractional[0] * b1[1] + fractional[1] * b2[1] + fractional[2] * b3[1],
    fractional[0] * b1[2] + fractional[1] * b2[2] + fractional[2] * b3[2],
  ];
}

/**
 * Compute the metric tensor for a lattice.
 * G_ij = a_i · a_j
 */
export function metricTensor(lattice: Matrix3x3): Matrix3x3 {
  const [a, b, c] = lattice;
  return [
    [dot(a, a), dot(a, b), dot(a, c)],
    [dot(b, a), dot(b, b), dot(b, c)],
    [dot(c, a), dot(c, b), dot(c, c)],
  ];
}

/**
 * Centering types for Bravais lattices.
 * P = primitive, F = face-centered, I = body-centered,
 * C = base-centered (C-face), A = base-centered (A-face), B = base-centered (B-face)
 * R = rhombohedral
 */
export type CenteringType = "P" | "F" | "I" | "C" | "A" | "B" | "R";

/**
 * Convert conventional cell vectors to primitive cell vectors.
 *
 * The Brillouin zone must be computed from the primitive cell, not the conventional cell.
 * For centered lattices (F, I, C, etc.), this transformation is essential.
 *
 * @param conventional - Conventional cell vectors [a, b, c]
 * @param centering - Centering type
 * @returns Primitive cell vectors
 */
export function conventionalToPrimitive(
  conventional: Matrix3x3,
  centering: CenteringType
): Matrix3x3 {
  const [a, b, c] = conventional;

  switch (centering) {
    case "P":
      // Primitive - no transformation needed
      return conventional;

    case "F":
      // Face-centered: primitive vectors are half face diagonals
      // a' = (b + c) / 2
      // b' = (a + c) / 2
      // c' = (a + b) / 2
      return [
        scale(add(b, c), 0.5),
        scale(add(a, c), 0.5),
        scale(add(a, b), 0.5),
      ];

    case "I":
      // Body-centered: primitive vectors involve the body diagonal
      // a' = (-a + b + c) / 2
      // b' = (a - b + c) / 2
      // c' = (a + b - c) / 2
      return [
        scale(add(add(scale(a, -1), b), c), 0.5),
        scale(add(add(a, scale(b, -1)), c), 0.5),
        scale(add(add(a, b), scale(c, -1)), 0.5),
      ];

    case "C":
      // C-centered (ab-face centered)
      // a' = (a + b) / 2
      // b' = (-a + b) / 2
      // c' = c
      return [
        scale(add(a, b), 0.5),
        scale(add(scale(a, -1), b), 0.5),
        c,
      ];

    case "A":
      // A-centered (bc-face centered)
      // a' = a
      // b' = (b + c) / 2
      // c' = (-b + c) / 2
      return [
        a,
        scale(add(b, c), 0.5),
        scale(add(scale(b, -1), c), 0.5),
      ];

    case "B":
      // B-centered (ac-face centered)
      // a' = (a + c) / 2
      // b' = b
      // c' = (-a + c) / 2
      return [
        scale(add(a, c), 0.5),
        b,
        scale(add(scale(a, -1), c), 0.5),
      ];

    case "R":
      // Rhombohedral - already primitive in rhombohedral setting
      // If given in hexagonal setting, would need different transformation
      return conventional;

    default:
      return conventional;
  }
}

// ============================================================================
// Brillouin Zone Calculation (Wigner-Seitz cell of reciprocal lattice)
// ============================================================================

/**
 * Represents a plane in 3D space: normal · r = d
 */
interface Plane {
  normal: Vec3;  // The G vector (not normalized)
  d: number;     // |G|²/2
  hkl: [number, number, number];  // Miller indices for identification
}

/**
 * A corner of the Brillouin zone
 */
interface BZCorner {
  position: Vec3;
  planes: Plane[];  // The 3 planes that intersect to form this corner
}

/**
 * Result of BZ calculation
 */
export interface BrillouinZoneGeometry {
  vertices: Vec3[];
  edges: [number, number][];
}

/**
 * Subtract two vectors.
 */
export function subtract(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

/**
 * Solve a 3x3 linear system Ax = b using Cramer's rule.
 * Returns null if the system is singular (determinant ≈ 0).
 */
function solve3x3(A: Matrix3x3, b: Vec3): Vec3 | null {
  const [r0, r1, r2] = A;

  // Compute determinant of A
  const detA =
    r0[0] * (r1[1] * r2[2] - r1[2] * r2[1]) -
    r0[1] * (r1[0] * r2[2] - r1[2] * r2[0]) +
    r0[2] * (r1[0] * r2[1] - r1[1] * r2[0]);

  if (Math.abs(detA) < 1e-10) {
    return null; // Singular matrix
  }

  // Cramer's rule
  const detX =
    b[0] * (r1[1] * r2[2] - r1[2] * r2[1]) -
    r0[1] * (b[1] * r2[2] - r1[2] * b[2]) +
    r0[2] * (b[1] * r2[1] - r1[1] * b[2]);

  const detY =
    r0[0] * (b[1] * r2[2] - r1[2] * b[2]) -
    b[0] * (r1[0] * r2[2] - r1[2] * r2[0]) +
    r0[2] * (r1[0] * b[2] - b[1] * r2[0]);

  const detZ =
    r0[0] * (r1[1] * b[2] - b[1] * r2[1]) -
    r0[1] * (r1[0] * b[2] - b[1] * r2[0]) +
    b[0] * (r1[0] * r2[1] - r1[1] * r2[0]);

  return [detX / detA, detY / detA, detZ / detA];
}

/**
 * Check if two planes are the same (same G vector direction).
 */
function planesEqual(p1: Plane, p2: Plane): boolean {
  return p1.hkl[0] === p2.hkl[0] &&
         p1.hkl[1] === p2.hkl[1] &&
         p1.hkl[2] === p2.hkl[2];
}

/**
 * Count how many planes two corners share.
 */
function countSharedPlanes(c1: BZCorner, c2: BZCorner): number {
  let count = 0;
  for (const p1 of c1.planes) {
    for (const p2 of c2.planes) {
      if (planesEqual(p1, p2)) {
        count++;
        break;
      }
    }
  }
  return count;
}

/**
 * Check if two Vec3 positions are approximately equal.
 */
function positionsEqual(a: Vec3, b: Vec3, tol: number = 1e-6): boolean {
  return Math.abs(a[0] - b[0]) < tol &&
         Math.abs(a[1] - b[1]) < tol &&
         Math.abs(a[2] - b[2]) < tol;
}

/**
 * Calculate the first Brillouin zone geometry from reciprocal lattice vectors.
 *
 * Algorithm (Wigner-Seitz construction):
 * 1. Generate reciprocal lattice vectors G_hkl = h*b1 + k*b2 + l*b3
 * 2. Each G defines a plane perpendicular to G, passing through G/2
 *    Plane equation: G · k = |G|²/2
 * 3. Find corners by solving intersections of 3 planes
 * 4. Keep only corners closer to Γ than to any other lattice point
 * 5. Connect corners that share exactly 2 planes to form edges
 *
 * Reference: https://lampz.tugraz.at/~hadley/ss1/bzones/
 */
export function calculateBrillouinZone(reciprocalBasis: Matrix3x3): BrillouinZoneGeometry {
  const [b1, b2, b3] = reciprocalBasis;

  // Step 1: Generate all G vectors for h,k,l ∈ {-1, 0, 1}, excluding (0,0,0)
  // We also include h,k,l ∈ {-2, -1, 0, 1, 2} to handle edge cases
  const gVectors: Plane[] = [];

  for (let h = -2; h <= 2; h++) {
    for (let k = -2; k <= 2; k++) {
      for (let l = -2; l <= 2; l++) {
        if (h === 0 && k === 0 && l === 0) continue;

        const G: Vec3 = [
          h * b1[0] + k * b2[0] + l * b3[0],
          h * b1[1] + k * b2[1] + l * b3[1],
          h * b1[2] + k * b2[2] + l * b3[2],
        ];

        const G2 = dot(G, G);

        gVectors.push({
          normal: G,
          d: G2 / 2,
          hkl: [h, k, l],
        });
      }
    }
  }

  // Sort by |G|² to process shorter G vectors first (more likely to be boundaries)
  gVectors.sort((a, b) => dot(a.normal, a.normal) - dot(b.normal, b.normal));

  // Step 2: Find candidate boundary planes
  // A plane is a boundary if G/2 is closer to Γ than to any other lattice point
  // This is equivalent to checking that |G| is minimal in its direction
  const boundaryPlanes: Plane[] = [];

  for (const plane of gVectors) {
    const G = plane.normal;

    // Check if this G/2 point is closer to origin than to any other G
    let isBoundary = true;
    const halfG: Vec3 = [G[0] / 2, G[1] / 2, G[2] / 2];

    for (const other of gVectors) {
      if (other === plane) continue;

      const otherG = other.normal;
      const diff = subtract(halfG, otherG);
      const distToOther = dot(diff, diff);
      const distToOrigin = dot(halfG, halfG);

      if (distToOther < distToOrigin - 1e-10) {
        isBoundary = false;
        break;
      }
    }

    if (isBoundary) {
      boundaryPlanes.push(plane);
    }
  }

  // Step 3: Find corners by solving intersections of 3 boundary planes
  const corners: BZCorner[] = [];

  for (let i = 0; i < boundaryPlanes.length; i++) {
    for (let j = i + 1; j < boundaryPlanes.length; j++) {
      for (let k = j + 1; k < boundaryPlanes.length; k++) {
        const p1 = boundaryPlanes[i];
        const p2 = boundaryPlanes[j];
        const p3 = boundaryPlanes[k];

        // Solve the system:
        // G1 · r = |G1|²/2
        // G2 · r = |G2|²/2
        // G3 · r = |G3|²/2
        const A: Matrix3x3 = [p1.normal, p2.normal, p3.normal];
        const b: Vec3 = [p1.d, p2.d, p3.d];

        const solution = solve3x3(A, b);
        if (!solution) continue; // Planes don't intersect at a point

        // Check if this corner is inside the BZ (closer to Γ than to any G)
        let isInside = true;
        for (const plane of gVectors) {
          const G = plane.normal;
          const diff = subtract(solution, G);
          const distToG = dot(diff, diff);
          const distToOrigin = dot(solution, solution);

          if (distToG < distToOrigin - 1e-10) {
            isInside = false;
            break;
          }
        }

        if (!isInside) continue;

        // Check if this corner already exists
        let isDuplicate = false;
        for (const existing of corners) {
          if (positionsEqual(existing.position, solution)) {
            isDuplicate = true;
            break;
          }
        }

        if (!isDuplicate) {
          corners.push({
            position: solution,
            planes: [p1, p2, p3],
          });
        }
      }
    }
  }

  // Step 4: Build edges by connecting corners that share exactly 2 planes
  const edges: [number, number][] = [];

  for (let i = 0; i < corners.length; i++) {
    for (let j = i + 1; j < corners.length; j++) {
      const sharedPlanes = countSharedPlanes(corners[i], corners[j]);
      if (sharedPlanes === 2) {
        edges.push([i, j]);
      }
    }
  }

  // Extract vertex positions
  const vertices: Vec3[] = corners.map(c => c.position);

  return { vertices, edges };
}
