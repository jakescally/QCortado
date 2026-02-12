/**
 * Primitive cell conversion for band structure calculations.
 *
 * For accurate band structure calculations that match literature,
 * we need to use the primitive cell with the appropriate ibrav value,
 * not the conventional cell with ibrav=0.
 *
 * This module detects lattice types and converts to primitive cells.
 */

import { CrystalData } from "./types";

/**
 * Strip oxidation state and other suffixes from element symbol.
 * e.g., "Si0+" -> "Si", "Fe3+" -> "Fe", "O2-" -> "O"
 */
function getBaseElement(symbol: string): string {
  return symbol.replace(/[\d+-]+$/, "");
}

/**
 * Normalize Hermann-Mauguin space group strings so common formatting variants
 * (spaces, underscores, unicode minus) can be matched reliably.
 */
function normalizeSpaceGroupHM(value: string): string {
  return value
    .toLowerCase()
    .replace(/[−–—]/g, "-")
    .replace(/[\s_]/g, "");
}

/**
 * Check if a structure is FCC diamond (like silicon, germanium, diamond).
 * Space group 227 (Fd-3m) is the diamond structure.
 */
export function isFCCDiamond(crystalData: CrystalData): boolean {
  const sgHM = crystalData.space_group_HM || "";
  const sgNum = crystalData.space_group_IT_number;
  const normalizedHM = normalizeSpaceGroupHM(sgHM);

  // Check by space group number first (most reliable)
  // 227: Fd-3m (diamond), 225: Fm-3m (FCC metals), 216: F-43m (zincblende)
  if (sgNum === 227 || sgNum === 225 || sgNum === 216) {
    return true;
  }

  // Check common HM notations for FCC/diamond/zincblende families.
  if (normalizedHM === "fd-3m" || normalizedHM === "fd3m") {
    return true;
  }
  if (normalizedHM === "fm-3m" || normalizedHM === "fm3m") {
    return true;
  }
  if (normalizedHM === "f-43m" || normalizedHM === "f43m") {
    return true;
  }

  // Also check if the string contains the number
  if (normalizedHM.includes("227") || normalizedHM.includes("225") || normalizedHM.includes("216")) {
    return true;
  }

  return false;
}

/**
 * Check if a structure is BCC (body-centered cubic).
 * Space group 229 (Im-3m) is common for BCC metals.
 */
export function isBCC(crystalData: CrystalData): boolean {
  const sgHM = crystalData.space_group_HM || "";
  const sgNum = crystalData.space_group_IT_number;
  const normalizedHM = normalizeSpaceGroupHM(sgHM);

  // Check by space group number first
  if (sgNum === 229) {
    return true;
  }

  if (normalizedHM === "im-3m" || normalizedHM === "im3m") {
    return true;
  }

  return false;
}

/**
 * Primitive cell data for QE calculations.
 */
export interface PrimitiveCell {
  /** ibrav enum value for QE (matches Rust BravaisLattice enum) */
  ibrav: "cubic_f" | "cubic_i" | "cubic_p";
  /** ibrav numeric value for logging (2 = FCC, 3 = BCC, etc.) */
  ibravNumeric: number;
  /** Lattice parameter in Bohr */
  celldm1: number;
  /** Atoms in primitive cell (crystal coordinates for the primitive cell) */
  atoms: { symbol: string; position: [number, number, number] }[];
  /** Number of atoms */
  nat: number;
}

/**
 * Convert conventional FCC diamond cell to primitive cell.
 *
 * For diamond structure (Si, Ge, C):
 * - Conventional cell: 8 atoms, cubic with lattice parameter a
 * - Primitive cell: 2 atoms at (0,0,0) and (1/4, 1/4, 1/4)
 *
 * The primitive FCC lattice vectors in terms of conventional a are:
 *   a1 = a/2 * (0, 1, 1)
 *   a2 = a/2 * (1, 0, 1)
 *   a3 = a/2 * (1, 1, 0)
 *
 * QE uses ibrav=2 with celldm(1) = a (in Bohr) for this.
 */
export function convertToFCCPrimitive(crystalData: CrystalData): PrimitiveCell | null {
  // Get the lattice parameter (a = b = c for cubic)
  const a = crystalData.cell_length_a.value;

  // Verify it's cubic
  const alpha = crystalData.cell_angle_alpha.value;
  const beta = crystalData.cell_angle_beta.value;
  const gamma = crystalData.cell_angle_gamma.value;

  if (Math.abs(alpha - 90) > 0.1 || Math.abs(beta - 90) > 0.1 || Math.abs(gamma - 90) > 0.1) {
    console.warn("Not a cubic cell, cannot convert to FCC primitive");
    return null;
  }

  // Convert a from Angstrom to Bohr (1 Bohr = 0.529177 Å)
  const BOHR_TO_ANGSTROM = 0.529177;
  const celldm1 = a / BOHR_TO_ANGSTROM;

  // Get unique element symbols (strip oxidation states)
  const elements = new Set(crystalData.atom_sites.map(site => getBaseElement(site.type_symbol)));

  // For diamond structure with one element type
  if (elements.size === 1) {
    const symbol = getBaseElement(crystalData.atom_sites[0].type_symbol);

    return {
      ibrav: "cubic_f",
      ibravNumeric: 2,
      celldm1,
      atoms: [
        { symbol, position: [0, 0, 0] },
        { symbol, position: [0.25, 0.25, 0.25] },
      ],
      nat: 2,
    };
  }

  // For zincblende structure (GaAs, ZnS, etc.) with two element types
  if (elements.size === 2) {
    const symbols = Array.from(elements);

    // Determine which element is at which site based on the CIF data
    // In zincblende, one element at (0,0,0) and one at (1/4,1/4,1/4)
    let symbol1 = symbols[0];
    let symbol2 = symbols[1];

    // Try to find which is which from the atom sites
    for (const site of crystalData.atom_sites) {
      const x = site.fract_x;
      const y = site.fract_y;
      const z = site.fract_z;

      // Check if this atom is at origin (0,0,0) or equivalent
      if (Math.abs(x) < 0.01 && Math.abs(y) < 0.01 && Math.abs(z) < 0.01) {
        symbol1 = getBaseElement(site.type_symbol);
      }
      // Check if at (0.25, 0.25, 0.25) or equivalent
      if (Math.abs(x - 0.25) < 0.01 && Math.abs(y - 0.25) < 0.01 && Math.abs(z - 0.25) < 0.01) {
        symbol2 = getBaseElement(site.type_symbol);
      }
    }

    return {
      ibrav: "cubic_f",
      ibravNumeric: 2,
      celldm1,
      atoms: [
        { symbol: symbol1, position: [0, 0, 0] },
        { symbol: symbol2, position: [0.25, 0.25, 0.25] },
      ],
      nat: 2,
    };
  }

  // More than 2 element types - not a simple diamond/zincblende structure
  console.warn("Complex structure with >2 elements, cannot auto-convert to primitive");
  return null;
}

/**
 * Convert conventional BCC cell to primitive cell.
 *
 * For BCC structure:
 * - Conventional cell: 2 atoms (corners + body center)
 * - Primitive cell: 1 atom
 *
 * QE uses ibrav=3 with celldm(1) = a (in Bohr) for BCC.
 */
export function convertToBCCPrimitive(crystalData: CrystalData): PrimitiveCell | null {
  const a = crystalData.cell_length_a.value;

  // Verify it's cubic
  const alpha = crystalData.cell_angle_alpha.value;
  const beta = crystalData.cell_angle_beta.value;
  const gamma = crystalData.cell_angle_gamma.value;

  if (Math.abs(alpha - 90) > 0.1 || Math.abs(beta - 90) > 0.1 || Math.abs(gamma - 90) > 0.1) {
    console.warn("Not a cubic cell, cannot convert to BCC primitive");
    return null;
  }

  // Convert a from Angstrom to Bohr
  const BOHR_TO_ANGSTROM = 0.529177;
  const celldm1 = a / BOHR_TO_ANGSTROM;

  // Get the element (strip oxidation state)
  const symbol = getBaseElement(crystalData.atom_sites[0].type_symbol);

  return {
    ibrav: "cubic_i",
    ibravNumeric: 3,
    celldm1,
    atoms: [
      { symbol, position: [0, 0, 0] },
    ],
    nat: 1,
  };
}

/**
 * Get the primitive cell for a crystal structure if it's a known lattice type.
 * Returns null if the structure should use the conventional cell (ibrav=0).
 */
export function getPrimitiveCell(crystalData: CrystalData): PrimitiveCell | null {
  if (isFCCDiamond(crystalData)) {
    return convertToFCCPrimitive(crystalData);
  }

  if (isBCC(crystalData)) {
    return convertToBCCPrimitive(crystalData);
  }

  return null;
}

/**
 * Standard k-point paths for different lattice types.
 * These are in crystal (reciprocal lattice) coordinates for the primitive cell.
 *
 * References:
 * - Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010)
 * - QE examples (example01 for silicon)
 */
export const STANDARD_KPATHS = {
  /** FCC (ibrav=2): L-Γ-X-U|K-Γ */
  FCC: [
    { label: "L", coords: [0.5, 0.5, 0.5] as [number, number, number] },
    { label: "Γ", coords: [0.0, 0.0, 0.0] as [number, number, number] },
    { label: "X", coords: [0.0, 0.5, 0.5] as [number, number, number] },
    { label: "U", coords: [0.25, 0.625, 0.625] as [number, number, number] },
    { label: "Γ", coords: [0.0, 0.0, 0.0] as [number, number, number] },
  ],

  /** Alternative FCC path: Γ-X-W-K-Γ-L-U-W-L-K */
  FCC_FULL: [
    { label: "Γ", coords: [0.0, 0.0, 0.0] as [number, number, number] },
    { label: "X", coords: [0.0, 0.5, 0.5] as [number, number, number] },
    { label: "W", coords: [0.25, 0.5, 0.75] as [number, number, number] },
    { label: "K", coords: [0.375, 0.375, 0.75] as [number, number, number] },
    { label: "Γ", coords: [0.0, 0.0, 0.0] as [number, number, number] },
    { label: "L", coords: [0.5, 0.5, 0.5] as [number, number, number] },
  ],

  /** BCC (ibrav=3): Γ-H-N-Γ-P-H */
  BCC: [
    { label: "Γ", coords: [0.0, 0.0, 0.0] as [number, number, number] },
    { label: "H", coords: [0.5, -0.5, 0.5] as [number, number, number] },
    { label: "N", coords: [0.0, 0.0, 0.5] as [number, number, number] },
    { label: "Γ", coords: [0.0, 0.0, 0.0] as [number, number, number] },
    { label: "P", coords: [0.25, 0.25, 0.25] as [number, number, number] },
  ],
};
