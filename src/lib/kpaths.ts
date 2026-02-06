// K-path database for band structure calculations
// Based on: Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010)

import { BravaisLattice } from './brillouinZone';

export interface HighSymmetryPoint {
  label: string;           // Display label (e.g., "Γ", "X", "M")
  qeLabel: string;         // Label for QE input (e.g., "gG", "X", "M")
  coords: [number, number, number];  // Crystal coordinates
  description: string;     // Human-readable description
}

export interface KPathSegment {
  from: string;  // Label of starting point
  to: string;    // Label of ending point
}

export interface StandardKPath {
  lattice: BravaisLattice;
  displayName: string;
  points: HighSymmetryPoint[];
  defaultPath: KPathSegment[];  // Default path as segments
  description: string;
}

// Helper to create a path string for display
export function pathToString(segments: KPathSegment[]): string {
  if (segments.length === 0) return '';

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

  const pointMap = new Map(allPoints.map(p => [p.label, p]));
  const result: { point: HighSymmetryPoint; npoints: number }[] = [];

  let lastTo = '';

  for (let i = 0; i < segments.length; i++) {
    const segment = segments[i];
    const isNewPath = segment.from !== lastTo;

    // Add "from" point if it's a new path or first segment
    if (isNewPath || i === 0) {
      const fromPoint = pointMap.get(segment.from);
      if (fromPoint) {
        result.push({ point: fromPoint, npoints: 20 }); // default 20 points per segment
      }
    }

    // Add "to" point
    const toPoint = pointMap.get(segment.to);
    if (toPoint) {
      // Last point in a disconnected segment or final point gets 0
      const isLastInPath = i === segments.length - 1 || segments[i + 1].from !== segment.to;
      result.push({ point: toPoint, npoints: isLastInPath ? 0 : 20 });
    }

    lastTo = segment.to;
  }

  return result;
}

// ============================================================================
// K-PATH DATABASE
// ============================================================================

const KPATH_CUBIC_P: StandardKPath = {
  lattice: 'cubic-P',
  displayName: 'Simple Cubic',
  description: 'Standard path for primitive cubic lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'X', qeLabel: 'X', coords: [0, 0.5, 0], description: 'Face center' },
    { label: 'M', qeLabel: 'M', coords: [0.5, 0.5, 0], description: 'Edge center' },
    { label: 'R', qeLabel: 'R', coords: [0.5, 0.5, 0.5], description: 'Corner' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'X' },
    { from: 'X', to: 'M' },
    { from: 'M', to: 'Γ' },
    { from: 'Γ', to: 'R' },
    { from: 'R', to: 'X' },
    { from: 'M', to: 'R' },
  ],
};

const KPATH_CUBIC_F: StandardKPath = {
  lattice: 'cubic-F',
  displayName: 'Face-Centered Cubic',
  description: 'Standard path for FCC lattices (e.g., Cu, Al, Si)',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'X', qeLabel: 'X', coords: [0.5, 0, 0.5], description: 'Face center' },
    { label: 'W', qeLabel: 'W', coords: [0.5, 0.25, 0.75], description: 'Corner of square face' },
    { label: 'K', qeLabel: 'K', coords: [0.375, 0.375, 0.75], description: 'Middle of edge' },
    { label: 'L', qeLabel: 'L', coords: [0.5, 0.5, 0.5], description: 'Center of hexagonal face' },
    { label: 'U', qeLabel: 'U', coords: [0.625, 0.25, 0.625], description: 'On square face edge' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'X' },
    { from: 'X', to: 'W' },
    { from: 'W', to: 'K' },
    { from: 'K', to: 'Γ' },
    { from: 'Γ', to: 'L' },
    { from: 'L', to: 'U' },
    { from: 'U', to: 'W' },
    { from: 'L', to: 'K' },
    { from: 'U', to: 'X' },
  ],
};

const KPATH_CUBIC_I: StandardKPath = {
  lattice: 'cubic-I',
  displayName: 'Body-Centered Cubic',
  description: 'Standard path for BCC lattices (e.g., Fe, W, Na)',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'H', qeLabel: 'H', coords: [0.5, -0.5, 0.5], description: 'Corner' },
    { label: 'N', qeLabel: 'N', coords: [0, 0, 0.5], description: 'Face center' },
    { label: 'P', qeLabel: 'P', coords: [0.25, 0.25, 0.25], description: 'Edge midpoint' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'H' },
    { from: 'H', to: 'N' },
    { from: 'N', to: 'Γ' },
    { from: 'Γ', to: 'P' },
    { from: 'P', to: 'H' },
    { from: 'P', to: 'N' },
  ],
};

const KPATH_HEXAGONAL: StandardKPath = {
  lattice: 'hexagonal',
  displayName: 'Hexagonal',
  description: 'Standard path for hexagonal lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'M', qeLabel: 'M', coords: [0.5, 0, 0], description: 'Edge midpoint' },
    { label: 'K', qeLabel: 'K', coords: [1/3, 1/3, 0], description: 'Corner' },
    { label: 'A', qeLabel: 'A', coords: [0, 0, 0.5], description: 'Top center' },
    { label: 'L', qeLabel: 'L', coords: [0.5, 0, 0.5], description: 'Top edge midpoint' },
    { label: 'H', qeLabel: 'H', coords: [1/3, 1/3, 0.5], description: 'Top corner' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'M' },
    { from: 'M', to: 'K' },
    { from: 'K', to: 'Γ' },
    { from: 'Γ', to: 'A' },
    { from: 'A', to: 'L' },
    { from: 'L', to: 'H' },
    { from: 'H', to: 'A' },
    { from: 'L', to: 'M' },
    { from: 'K', to: 'H' },
  ],
};

const KPATH_TETRAGONAL_P: StandardKPath = {
  lattice: 'tetragonal-P',
  displayName: 'Simple Tetragonal',
  description: 'Standard path for primitive tetragonal lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'X', qeLabel: 'X', coords: [0, 0.5, 0], description: 'Face center' },
    { label: 'M', qeLabel: 'M', coords: [0.5, 0.5, 0], description: 'Edge center' },
    { label: 'Z', qeLabel: 'Z', coords: [0, 0, 0.5], description: 'Top center' },
    { label: 'R', qeLabel: 'R', coords: [0, 0.5, 0.5], description: 'Top face center' },
    { label: 'A', qeLabel: 'A', coords: [0.5, 0.5, 0.5], description: 'Top corner' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'X' },
    { from: 'X', to: 'M' },
    { from: 'M', to: 'Γ' },
    { from: 'Γ', to: 'Z' },
    { from: 'Z', to: 'R' },
    { from: 'R', to: 'A' },
    { from: 'A', to: 'Z' },
    { from: 'X', to: 'R' },
    { from: 'M', to: 'A' },
  ],
};

const KPATH_TETRAGONAL_I: StandardKPath = {
  lattice: 'tetragonal-I',
  displayName: 'Body-Centered Tetragonal',
  description: 'Standard path for body-centered tetragonal lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'M', qeLabel: 'M', coords: [-0.5, 0.5, 0.5], description: 'Edge center' },
    { label: 'X', qeLabel: 'X', coords: [0, 0, 0.5], description: 'Face center' },
    { label: 'P', qeLabel: 'P', coords: [0.25, 0.25, 0.25], description: 'Body point' },
    { label: 'N', qeLabel: 'N', coords: [0, 0.5, 0], description: 'Face point' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'X' },
    { from: 'X', to: 'M' },
    { from: 'M', to: 'Γ' },
    { from: 'Γ', to: 'N' },
    { from: 'N', to: 'P' },
    { from: 'P', to: 'M' },
  ],
};

const KPATH_ORTHORHOMBIC_P: StandardKPath = {
  lattice: 'orthorhombic-P',
  displayName: 'Simple Orthorhombic',
  description: 'Standard path for primitive orthorhombic lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'X', qeLabel: 'X', coords: [0.5, 0, 0], description: 'Face center (a)' },
    { label: 'Y', qeLabel: 'Y', coords: [0, 0.5, 0], description: 'Face center (b)' },
    { label: 'Z', qeLabel: 'Z', coords: [0, 0, 0.5], description: 'Face center (c)' },
    { label: 'S', qeLabel: 'S', coords: [0.5, 0.5, 0], description: 'Edge (ab)' },
    { label: 'T', qeLabel: 'T', coords: [0, 0.5, 0.5], description: 'Edge (bc)' },
    { label: 'U', qeLabel: 'U', coords: [0.5, 0, 0.5], description: 'Edge (ac)' },
    { label: 'R', qeLabel: 'R', coords: [0.5, 0.5, 0.5], description: 'Corner' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'X' },
    { from: 'X', to: 'S' },
    { from: 'S', to: 'Y' },
    { from: 'Y', to: 'Γ' },
    { from: 'Γ', to: 'Z' },
    { from: 'Z', to: 'U' },
    { from: 'U', to: 'R' },
    { from: 'R', to: 'T' },
    { from: 'T', to: 'Z' },
  ],
};

const KPATH_TRIGONAL_R: StandardKPath = {
  lattice: 'trigonal-R',
  displayName: 'Rhombohedral',
  description: 'Standard path for rhombohedral lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'L', qeLabel: 'L', coords: [0.5, 0, 0], description: 'Face center' },
    { label: 'F', qeLabel: 'F', coords: [0.5, 0.5, 0], description: 'Edge point' },
    { label: 'Z', qeLabel: 'Z', coords: [0.5, 0.5, 0.5], description: 'Body point' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'L' },
    { from: 'L', to: 'F' },
    { from: 'F', to: 'Γ' },
    { from: 'Γ', to: 'Z' },
  ],
};

const KPATH_MONOCLINIC_P: StandardKPath = {
  lattice: 'monoclinic-P',
  displayName: 'Simple Monoclinic',
  description: 'Standard path for primitive monoclinic lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'Y', qeLabel: 'Y', coords: [0.5, 0, 0], description: 'Face center' },
    { label: 'A', qeLabel: 'A', coords: [0, 0.5, 0], description: 'Face center' },
    { label: 'B', qeLabel: 'B', coords: [0, 0, 0.5], description: 'Face center' },
    { label: 'C', qeLabel: 'C', coords: [0.5, 0.5, 0], description: 'Edge point' },
    { label: 'D', qeLabel: 'D', coords: [0.5, 0, 0.5], description: 'Edge point' },
    { label: 'E', qeLabel: 'E', coords: [0, 0.5, 0.5], description: 'Edge point' },
    { label: 'Z', qeLabel: 'Z', coords: [0.5, 0.5, 0.5], description: 'Corner' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'Y' },
    { from: 'Y', to: 'C' },
    { from: 'C', to: 'A' },
    { from: 'A', to: 'Γ' },
    { from: 'Γ', to: 'B' },
    { from: 'B', to: 'D' },
    { from: 'D', to: 'Z' },
    { from: 'Z', to: 'E' },
    { from: 'E', to: 'B' },
  ],
};

const KPATH_TRICLINIC: StandardKPath = {
  lattice: 'triclinic',
  displayName: 'Triclinic',
  description: 'Simple path for triclinic lattices',
  points: [
    { label: 'Γ', qeLabel: 'gG', coords: [0, 0, 0], description: 'Zone center' },
    { label: 'X', qeLabel: 'X', coords: [0.5, 0, 0], description: 'a* direction' },
    { label: 'Y', qeLabel: 'Y', coords: [0, 0.5, 0], description: 'b* direction' },
    { label: 'Z', qeLabel: 'Z', coords: [0, 0, 0.5], description: 'c* direction' },
  ],
  defaultPath: [
    { from: 'Γ', to: 'X' },
    { from: 'X', to: 'Γ' },
    { from: 'Γ', to: 'Y' },
    { from: 'Y', to: 'Γ' },
    { from: 'Γ', to: 'Z' },
  ],
};

// Fallback paths for lattice types not fully specified above
const KPATH_ORTHORHOMBIC_C: StandardKPath = {
  ...KPATH_ORTHORHOMBIC_P,
  lattice: 'orthorhombic-C',
  displayName: 'Base-Centered Orthorhombic',
  description: 'Standard path for C-centered orthorhombic lattices',
};

const KPATH_ORTHORHOMBIC_I: StandardKPath = {
  ...KPATH_ORTHORHOMBIC_P,
  lattice: 'orthorhombic-I',
  displayName: 'Body-Centered Orthorhombic',
  description: 'Standard path for I-centered orthorhombic lattices',
};

const KPATH_ORTHORHOMBIC_F: StandardKPath = {
  ...KPATH_ORTHORHOMBIC_P,
  lattice: 'orthorhombic-F',
  displayName: 'Face-Centered Orthorhombic',
  description: 'Standard path for F-centered orthorhombic lattices',
};

const KPATH_MONOCLINIC_C: StandardKPath = {
  ...KPATH_MONOCLINIC_P,
  lattice: 'monoclinic-C',
  displayName: 'Base-Centered Monoclinic',
  description: 'Standard path for C-centered monoclinic lattices',
};

// ============================================================================
// EXPORTS
// ============================================================================

export const KPATH_DATABASE: Record<BravaisLattice, StandardKPath> = {
  'cubic-P': KPATH_CUBIC_P,
  'cubic-F': KPATH_CUBIC_F,
  'cubic-I': KPATH_CUBIC_I,
  'tetragonal-P': KPATH_TETRAGONAL_P,
  'tetragonal-I': KPATH_TETRAGONAL_I,
  'orthorhombic-P': KPATH_ORTHORHOMBIC_P,
  'orthorhombic-C': KPATH_ORTHORHOMBIC_C,
  'orthorhombic-I': KPATH_ORTHORHOMBIC_I,
  'orthorhombic-F': KPATH_ORTHORHOMBIC_F,
  'hexagonal': KPATH_HEXAGONAL,
  'trigonal-R': KPATH_TRIGONAL_R,
  'monoclinic-P': KPATH_MONOCLINIC_P,
  'monoclinic-C': KPATH_MONOCLINIC_C,
  'triclinic': KPATH_TRICLINIC,
};

/**
 * Get the standard k-path for a given Bravais lattice
 */
export function getStandardKPath(lattice: BravaisLattice): StandardKPath {
  return KPATH_DATABASE[lattice];
}

/**
 * Format k-path for QE bands input
 */
export function formatKPathForQE(
  points: { point: HighSymmetryPoint; npoints: number }[]
): string {
  const lines: string[] = [];
  lines.push(`K_POINTS {crystal_b}`);
  lines.push(`${points.length}`);

  for (const { point, npoints } of points) {
    const [x, y, z] = point.coords;
    lines.push(`${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)} ${npoints}  ${point.qeLabel}`);
  }

  return lines.join('\n');
}
