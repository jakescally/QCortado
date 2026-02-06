// Crystal system and Bravais lattice detection for k-path selection

export type BravaisLattice =
  | 'cubic-P'        // Primitive cubic (simple cubic)
  | 'cubic-F'        // Face-centered cubic (fcc)
  | 'cubic-I'        // Body-centered cubic (bcc)
  | 'tetragonal-P'   // Primitive tetragonal
  | 'tetragonal-I'   // Body-centered tetragonal
  | 'orthorhombic-P' // Primitive orthorhombic
  | 'orthorhombic-C' // Base-centered orthorhombic
  | 'orthorhombic-I' // Body-centered orthorhombic
  | 'orthorhombic-F' // Face-centered orthorhombic
  | 'hexagonal'      // Hexagonal
  | 'trigonal-R'     // Rhombohedral (trigonal)
  | 'monoclinic-P'   // Primitive monoclinic
  | 'monoclinic-C'   // Base-centered monoclinic
  | 'triclinic';     // Triclinic

export type CrystalSystem =
  | 'triclinic'
  | 'monoclinic'
  | 'orthorhombic'
  | 'tetragonal'
  | 'trigonal'
  | 'hexagonal'
  | 'cubic';

/**
 * Determines the crystal system from space group number
 */
export function getCrystalSystem(spaceGroupNumber: number): CrystalSystem {
  if (spaceGroupNumber >= 1 && spaceGroupNumber <= 2) {
    return 'triclinic';
  } else if (spaceGroupNumber >= 3 && spaceGroupNumber <= 15) {
    return 'monoclinic';
  } else if (spaceGroupNumber >= 16 && spaceGroupNumber <= 74) {
    return 'orthorhombic';
  } else if (spaceGroupNumber >= 75 && spaceGroupNumber <= 142) {
    return 'tetragonal';
  } else if (spaceGroupNumber >= 143 && spaceGroupNumber <= 167) {
    return 'trigonal';
  } else if (spaceGroupNumber >= 168 && spaceGroupNumber <= 194) {
    return 'hexagonal';
  } else if (spaceGroupNumber >= 195 && spaceGroupNumber <= 230) {
    return 'cubic';
  }
  throw new Error(`Invalid space group number: ${spaceGroupNumber}`);
}

// Face-centered cubic space groups
const CUBIC_F_GROUPS = new Set([
  196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228
]);

// Body-centered cubic space groups
const CUBIC_I_GROUPS = new Set([
  197, 199, 204, 206, 211, 214, 217, 220, 229, 230
]);

// Body-centered tetragonal space groups
const TETRAGONAL_I_GROUPS = new Set([
  79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110,
  119, 120, 121, 122, 139, 140, 141, 142
]);

// C-centered orthorhombic (also includes A-centered, which are equivalent)
const ORTHORHOMBIC_C_GROUPS = new Set([
  20, 21, 35, 36, 37, 38, 39, 40, 41, 63, 64, 65, 66, 67, 68
]);

// Face-centered orthorhombic
const ORTHORHOMBIC_F_GROUPS = new Set([
  22, 42, 43, 69, 70
]);

// Body-centered orthorhombic
const ORTHORHOMBIC_I_GROUPS = new Set([
  23, 24, 44, 45, 46, 71, 72, 73, 74
]);

// C-centered monoclinic
const MONOCLINIC_C_GROUPS = new Set([
  5, 8, 9, 12, 15
]);

// Rhombohedral (R) trigonal space groups
const TRIGONAL_R_GROUPS = new Set([
  146, 148, 155, 160, 161, 166, 167
]);

/**
 * Determines the Bravais lattice type from space group number.
 * Uses the International Tables for Crystallography conventions.
 */
export function detectBravaisLattice(spaceGroupNumber: number): BravaisLattice {
  const crystalSystem = getCrystalSystem(spaceGroupNumber);

  switch (crystalSystem) {
    case 'triclinic':
      return 'triclinic';

    case 'monoclinic':
      if (MONOCLINIC_C_GROUPS.has(spaceGroupNumber)) {
        return 'monoclinic-C';
      }
      return 'monoclinic-P';

    case 'orthorhombic':
      if (ORTHORHOMBIC_F_GROUPS.has(spaceGroupNumber)) {
        return 'orthorhombic-F';
      }
      if (ORTHORHOMBIC_I_GROUPS.has(spaceGroupNumber)) {
        return 'orthorhombic-I';
      }
      if (ORTHORHOMBIC_C_GROUPS.has(spaceGroupNumber)) {
        return 'orthorhombic-C';
      }
      return 'orthorhombic-P';

    case 'tetragonal':
      if (TETRAGONAL_I_GROUPS.has(spaceGroupNumber)) {
        return 'tetragonal-I';
      }
      return 'tetragonal-P';

    case 'trigonal':
      if (TRIGONAL_R_GROUPS.has(spaceGroupNumber)) {
        return 'trigonal-R';
      }
      // Non-rhombohedral trigonal uses hexagonal lattice
      return 'hexagonal';

    case 'hexagonal':
      return 'hexagonal';

    case 'cubic':
      if (CUBIC_F_GROUPS.has(spaceGroupNumber)) {
        return 'cubic-F';
      }
      if (CUBIC_I_GROUPS.has(spaceGroupNumber)) {
        return 'cubic-I';
      }
      return 'cubic-P';

    default:
      return 'triclinic';
  }
}

/**
 * Returns a human-readable name for the Bravais lattice
 */
export function getBravaisLatticeName(lattice: BravaisLattice): string {
  const names: Record<BravaisLattice, string> = {
    'cubic-P': 'Simple Cubic (cP)',
    'cubic-F': 'Face-Centered Cubic (cF)',
    'cubic-I': 'Body-Centered Cubic (cI)',
    'tetragonal-P': 'Simple Tetragonal (tP)',
    'tetragonal-I': 'Body-Centered Tetragonal (tI)',
    'orthorhombic-P': 'Simple Orthorhombic (oP)',
    'orthorhombic-C': 'Base-Centered Orthorhombic (oC)',
    'orthorhombic-I': 'Body-Centered Orthorhombic (oI)',
    'orthorhombic-F': 'Face-Centered Orthorhombic (oF)',
    'hexagonal': 'Hexagonal (hP)',
    'trigonal-R': 'Rhombohedral (hR)',
    'monoclinic-P': 'Simple Monoclinic (mP)',
    'monoclinic-C': 'Base-Centered Monoclinic (mC)',
    'triclinic': 'Triclinic (aP)',
  };
  return names[lattice];
}

/**
 * Returns the crystal system name
 */
export function getCrystalSystemName(system: CrystalSystem): string {
  const names: Record<CrystalSystem, string> = {
    triclinic: 'Triclinic',
    monoclinic: 'Monoclinic',
    orthorhombic: 'Orthorhombic',
    tetragonal: 'Tetragonal',
    trigonal: 'Trigonal',
    hexagonal: 'Hexagonal',
    cubic: 'Cubic',
  };
  return names[system];
}
