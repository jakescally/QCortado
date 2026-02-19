// K-path compatibility adapter.
//
// Canonical source of high-symmetry points/paths now lives in
// `brillouinZoneData.ts`, which is also used by the BZ viewer.

import { BravaisLattice } from "./brillouinZone";
import { BravaisLatticeType, getBrillouinZoneData } from "./brillouinZoneData";

export interface HighSymmetryPoint {
  label: string; // Display label (e.g., "Γ", "X", "M")
  qeLabel: string; // Label for QE text output (e.g., "gG", "X", "M")
  coords: [number, number, number]; // Crystal coordinates
  description: string; // Human-readable description
}

export interface KPathSegment {
  from: string;
  to: string;
}

export interface StandardKPath {
  lattice: BravaisLattice;
  displayName: string;
  points: HighSymmetryPoint[];
  defaultPath: KPathSegment[];
  description: string;
}

export interface KPathLatticeParams {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

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

const BRAVAIS_TO_BZ: Record<BravaisLattice, BravaisLatticeType> = {
  "cubic-P": "cP",
  "cubic-F": "cF",
  "cubic-I": "cI",
  "tetragonal-P": "tP",
  "tetragonal-I": "tI",
  "orthorhombic-P": "oP",
  "orthorhombic-C": "oC",
  "orthorhombic-I": "oI",
  "orthorhombic-F": "oF",
  "hexagonal": "hP",
  "trigonal-R": "hR",
  "monoclinic-P": "mP",
  "monoclinic-C": "mC",
  "triclinic": "aP",
};

const BRAVAIS_DISPLAY_NAME: Record<BravaisLattice, string> = {
  "cubic-P": "Simple Cubic",
  "cubic-F": "Face-Centered Cubic",
  "cubic-I": "Body-Centered Cubic",
  "tetragonal-P": "Simple Tetragonal",
  "tetragonal-I": "Body-Centered Tetragonal",
  "orthorhombic-P": "Simple Orthorhombic",
  "orthorhombic-C": "Base-Centered Orthorhombic",
  "orthorhombic-I": "Body-Centered Orthorhombic",
  "orthorhombic-F": "Face-Centered Orthorhombic",
  "hexagonal": "Hexagonal",
  "trigonal-R": "Rhombohedral",
  "monoclinic-P": "Simple Monoclinic",
  "monoclinic-C": "Base-Centered Monoclinic",
  "triclinic": "Triclinic",
};

const DEFAULT_PARAMS: Record<BravaisLattice, KPathLatticeParams> = {
  "cubic-P": { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 },
  "cubic-F": { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 },
  "cubic-I": { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 },
  "tetragonal-P": { a: 1, b: 1, c: 1.3, alpha: 90, beta: 90, gamma: 90 },
  "tetragonal-I": { a: 1, b: 1, c: 1.3, alpha: 90, beta: 90, gamma: 90 },
  "orthorhombic-P": { a: 1, b: 1.15, c: 1.4, alpha: 90, beta: 90, gamma: 90 },
  "orthorhombic-C": { a: 1, b: 1.15, c: 1.4, alpha: 90, beta: 90, gamma: 90 },
  "orthorhombic-I": { a: 1, b: 1.15, c: 1.4, alpha: 90, beta: 90, gamma: 90 },
  "orthorhombic-F": { a: 1, b: 1.15, c: 1.4, alpha: 90, beta: 90, gamma: 90 },
  "hexagonal": { a: 1, b: 1, c: 1.6, alpha: 90, beta: 90, gamma: 120 },
  // hR defaults to alpha < 90 so the richer hR1 point set is exposed by default.
  "trigonal-R": { a: 1, b: 1, c: 1, alpha: 78, beta: 78, gamma: 78 },
  "monoclinic-P": { a: 1, b: 1.1, c: 1.3, alpha: 90, beta: 105, gamma: 90 },
  "monoclinic-C": { a: 1, b: 1.1, c: 1.3, alpha: 90, beta: 105, gamma: 90 },
  "triclinic": { a: 1, b: 1.07, c: 1.21, alpha: 82, beta: 96, gamma: 108 },
};

function toQeLabel(label: string): string {
  if (label === "Γ") {
    return "gG";
  }
  return label
    .replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (ch) => SUBSCRIPT_TO_ASCII[ch] ?? ch)
    .replace(/[^A-Za-z0-9]/g, "");
}

function resolveParams(lattice: BravaisLattice, overrides?: Partial<KPathLatticeParams>): KPathLatticeParams {
  return {
    ...DEFAULT_PARAMS[lattice],
    ...overrides,
  };
}

// Helper to create a path string for display
export function pathToString(segments: KPathSegment[]): string {
  if (segments.length === 0) return "";

  let result = segments[0].from;
  let lastTo = segments[0].from;

  for (const segment of segments) {
    if (segment.from === lastTo) {
      result += ` → ${segment.to}`;
    } else {
      result += ` | ${segment.from} → ${segment.to}`;
    }
    lastTo = segment.to;
  }

  return result;
}

// Helper to get points in order from segments
export function getOrderedPoints(
  segments: KPathSegment[],
  allPoints: HighSymmetryPoint[]
): { point: HighSymmetryPoint; npoints: number }[] {
  if (segments.length === 0) return [];

  const pointMap = new Map(allPoints.map((p) => [p.label, p]));
  const result: { point: HighSymmetryPoint; npoints: number }[] = [];

  let lastTo = "";

  for (let i = 0; i < segments.length; i++) {
    const segment = segments[i];
    const isNewPath = segment.from !== lastTo;

    // Add "from" point if it's a new path or first segment
    if (isNewPath || i === 0) {
      const fromPoint = pointMap.get(segment.from);
      if (fromPoint) {
        result.push({ point: fromPoint, npoints: 20 });
      }
    }

    // Add "to" point
    const toPoint = pointMap.get(segment.to);
    if (toPoint) {
      const isLastInPath = i === segments.length - 1 || segments[i + 1].from !== segment.to;
      result.push({ point: toPoint, npoints: isLastInPath ? 0 : 20 });
    }

    lastTo = segment.to;
  }

  return result;
}

function buildStandardKPath(lattice: BravaisLattice, params: KPathLatticeParams): StandardKPath {
  const bzType = BRAVAIS_TO_BZ[lattice];
  const bzData = getBrillouinZoneData(bzType, params);
  return {
    lattice,
    displayName: BRAVAIS_DISPLAY_NAME[lattice],
    description: `Canonical path for ${BRAVAIS_DISPLAY_NAME[lattice]} (${bzData.latticeType}).`,
    points: bzData.points.map((p) => ({
      label: p.label,
      qeLabel: toQeLabel(p.label),
      coords: [p.coords[0], p.coords[1], p.coords[2]],
      description: p.description,
    })),
    defaultPath: bzData.recommendedPath.map(([from, to]) => ({ from, to })),
  };
}

/**
 * Compatibility export for previous API shape.
 *
 * This material-agnostic table uses representative lattice parameters per
 * Bravais family. Use `getStandardKPath(lattice, params)` when branch-sensitive
 * families (hR, oF, mC, tI) require exact metric-dependent paths.
 */
export const KPATH_DATABASE: Record<BravaisLattice, StandardKPath> = {
  "cubic-P": buildStandardKPath("cubic-P", DEFAULT_PARAMS["cubic-P"]),
  "cubic-F": buildStandardKPath("cubic-F", DEFAULT_PARAMS["cubic-F"]),
  "cubic-I": buildStandardKPath("cubic-I", DEFAULT_PARAMS["cubic-I"]),
  "tetragonal-P": buildStandardKPath("tetragonal-P", DEFAULT_PARAMS["tetragonal-P"]),
  "tetragonal-I": buildStandardKPath("tetragonal-I", DEFAULT_PARAMS["tetragonal-I"]),
  "orthorhombic-P": buildStandardKPath("orthorhombic-P", DEFAULT_PARAMS["orthorhombic-P"]),
  "orthorhombic-C": buildStandardKPath("orthorhombic-C", DEFAULT_PARAMS["orthorhombic-C"]),
  "orthorhombic-I": buildStandardKPath("orthorhombic-I", DEFAULT_PARAMS["orthorhombic-I"]),
  "orthorhombic-F": buildStandardKPath("orthorhombic-F", DEFAULT_PARAMS["orthorhombic-F"]),
  "hexagonal": buildStandardKPath("hexagonal", DEFAULT_PARAMS["hexagonal"]),
  "trigonal-R": buildStandardKPath("trigonal-R", DEFAULT_PARAMS["trigonal-R"]),
  "monoclinic-P": buildStandardKPath("monoclinic-P", DEFAULT_PARAMS["monoclinic-P"]),
  "monoclinic-C": buildStandardKPath("monoclinic-C", DEFAULT_PARAMS["monoclinic-C"]),
  triclinic: buildStandardKPath("triclinic", DEFAULT_PARAMS["triclinic"]),
};

/**
 * Get the standard k-path for a Bravais lattice.
 *
 * Optional lattice parameters let callers resolve branch-dependent conventions
 * (e.g., hR1/hR2, oF1/oF2/oF3, mC1-5).
 */
export function getStandardKPath(
  lattice: BravaisLattice,
  params?: Partial<KPathLatticeParams>
): StandardKPath {
  return buildStandardKPath(lattice, resolveParams(lattice, params));
}

/**
 * Format k-path for QE bands input.
 */
export function formatKPathForQE(
  points: { point: HighSymmetryPoint; npoints: number }[]
): string {
  const lines: string[] = [];
  lines.push("K_POINTS {crystal_b}");
  lines.push(`${points.length}`);

  for (const { point, npoints } of points) {
    const [x, y, z] = point.coords;
    lines.push(`${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)} ${npoints}  ${point.qeLabel}`);
  }

  return lines.join("\n");
}
