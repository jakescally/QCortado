import { invoke } from "@tauri-apps/api/core";
import { CrystalData } from "./types";
import { Matrix3x3, Vec3 } from "./reciprocalLattice";

export interface SymmetryAtomInput {
  symbol: string;
  position: Vec3;
}

export interface SymmetryAtomOutput {
  symbol: string;
  position: Vec3;
  typeIndex: number;
}

export interface SymmetryTransformResult {
  spacegroupNumber: number;
  hallNumber: number;
  internationalSymbol: string;
  choice: string;
  inputLattice: Matrix3x3;
  standardizedConventionalLattice: Matrix3x3;
  standardizedPrimitiveLattice: Matrix3x3;
  standardizedConventionalAtoms: SymmetryAtomOutput[];
  standardizedPrimitiveAtoms: SymmetryAtomOutput[];
  primitiveToInputReciprocal: Matrix3x3;
  inputToPrimitiveReciprocal: Matrix3x3;
  primitiveToStandardizedConventionalReciprocal: Matrix3x3;
  standardizedConventionalToPrimitiveReciprocal: Matrix3x3;
  transformationMatrix: Matrix3x3;
  originShift: Vec3;
}

export interface SymmetryAnalyzeRequest {
  lattice: Matrix3x3;
  atoms: SymmetryAtomInput[];
  symprec?: number;
  angleTolerance?: number;
}

const RELAXED_SYMPREC = 1e-3;

function getBaseElement(symbol: string): string {
  return symbol.replace(/[\d+-]+$/, "");
}

function coerceSpaceGroupNumber(value: unknown): number | null {
  if (value == null) return null;
  const parsed = typeof value === "number" ? value : Number.parseInt(String(value).trim(), 10);
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

function crystalSystemRank(spaceGroupNumber: number): number {
  if (spaceGroupNumber >= 1 && spaceGroupNumber <= 2) return 0; // triclinic
  if (spaceGroupNumber >= 3 && spaceGroupNumber <= 15) return 1; // monoclinic
  if (spaceGroupNumber >= 16 && spaceGroupNumber <= 74) return 2; // orthorhombic
  if (spaceGroupNumber >= 75 && spaceGroupNumber <= 142) return 3; // tetragonal
  if (spaceGroupNumber >= 143 && spaceGroupNumber <= 167) return 4; // trigonal
  if (spaceGroupNumber >= 168 && spaceGroupNumber <= 194) return 5; // hexagonal
  if (spaceGroupNumber >= 195 && spaceGroupNumber <= 230) return 6; // cubic
  return -1;
}

function inferSpaceGroupHint(crystalData: CrystalData): number | null {
  return (
    coerceSpaceGroupNumber(crystalData.space_group_IT_number) ??
    extractSpaceGroupNumberFromHM(crystalData.space_group_HM)
  );
}

function mismatchScore(spaceGroup: number, hintedSpaceGroup: number): number {
  const rankDelta = Math.abs(
    crystalSystemRank(spaceGroup) - crystalSystemRank(hintedSpaceGroup),
  );
  const numberDelta = Math.abs(spaceGroup - hintedSpaceGroup);
  return rankDelta * 1000 + numberDelta;
}

function shouldRetryWithRelaxedSymprec(
  strictResult: SymmetryTransformResult,
  hintedSpaceGroup: number | null,
  symprec: number,
): boolean {
  if (hintedSpaceGroup == null) return false;
  if (symprec >= RELAXED_SYMPREC) return false;
  if (strictResult.spacegroupNumber === hintedSpaceGroup) return false;

  const strictRank = crystalSystemRank(strictResult.spacegroupNumber);
  const hintRank = crystalSystemRank(hintedSpaceGroup);

  // Retry when strict analysis collapses to much lower symmetry than CIF metadata.
  return strictRank + 1 < hintRank || strictResult.spacegroupNumber <= 15;
}

function shouldPreferRelaxedResult(
  strictResult: SymmetryTransformResult,
  relaxedResult: SymmetryTransformResult,
  hintedSpaceGroup: number | null,
): boolean {
  if (hintedSpaceGroup != null) {
    const strictScore = mismatchScore(strictResult.spacegroupNumber, hintedSpaceGroup);
    const relaxedScore = mismatchScore(relaxedResult.spacegroupNumber, hintedSpaceGroup);
    return relaxedScore < strictScore;
  }
  return crystalSystemRank(relaxedResult.spacegroupNumber) > crystalSystemRank(strictResult.spacegroupNumber);
}

export function buildConventionalLatticeFromCrystalData(crystalData: CrystalData): Matrix3x3 {
  const a = crystalData.cell_length_a.value;
  const b = crystalData.cell_length_b.value;
  const c = crystalData.cell_length_c.value;
  const alpha = (crystalData.cell_angle_alpha.value * Math.PI) / 180;
  const beta = (crystalData.cell_angle_beta.value * Math.PI) / 180;
  const gamma = (crystalData.cell_angle_gamma.value * Math.PI) / 180;

  const aVec: Vec3 = [a, 0, 0];
  const bVec: Vec3 = [b * Math.cos(gamma), b * Math.sin(gamma), 0];
  const cx = c * Math.cos(beta);
  const cy = (c * (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma))) / Math.sin(gamma);
  const cz = Math.sqrt(Math.max(0, c * c - cx * cx - cy * cy));
  const cVec: Vec3 = [cx, cy, cz];

  return [aVec, bVec, cVec];
}

export function buildSymmetryRequestFromCrystalData(
  crystalData: CrystalData,
  symprec = 1e-5,
  angleTolerance = -1,
): SymmetryAnalyzeRequest {
  return {
    lattice: buildConventionalLatticeFromCrystalData(crystalData),
    atoms: crystalData.atom_sites.map((site) => ({
      symbol: getBaseElement(site.type_symbol),
      position: [site.fract_x, site.fract_y, site.fract_z] as Vec3,
    })),
    symprec,
    angleTolerance,
  };
}

const symmetryCache = new Map<string, Promise<SymmetryTransformResult>>();

export async function analyzeCrystalSymmetry(
  crystalData: CrystalData,
  symprec = 1e-5,
  angleTolerance = -1,
): Promise<SymmetryTransformResult> {
  const request = buildSymmetryRequestFromCrystalData(crystalData, symprec, angleTolerance);
  const cacheKey = JSON.stringify(request);
  const cached = symmetryCache.get(cacheKey);
  if (cached) return cached;

  const pending = (async () => {
    const strictResult = await invoke<SymmetryTransformResult>("analyze_structure_symmetry", {
      input: request,
    });

    const hintedSpaceGroup = inferSpaceGroupHint(crystalData);
    if (!shouldRetryWithRelaxedSymprec(strictResult, hintedSpaceGroup, symprec)) {
      return strictResult;
    }

    const relaxedSymprec = Math.max(symprec, RELAXED_SYMPREC);
    const relaxedRequest = buildSymmetryRequestFromCrystalData(
      crystalData,
      relaxedSymprec,
      angleTolerance,
    );

    try {
      const relaxedResult = await invoke<SymmetryTransformResult>("analyze_structure_symmetry", {
        input: relaxedRequest,
      });
      if (shouldPreferRelaxedResult(strictResult, relaxedResult, hintedSpaceGroup)) {
        console.info(
          `Symmetry analysis switched to relaxed symprec=${relaxedSymprec} (strict sg=${strictResult.spacegroupNumber}, relaxed sg=${relaxedResult.spacegroupNumber}, hinted sg=${hintedSpaceGroup}).`,
        );
        return relaxedResult;
      }
    } catch (error) {
      console.warn(
        `Relaxed symmetry analysis failed at symprec=${relaxedSymprec}; keeping strict result.`,
        error,
      );
    }

    return strictResult;
  })().catch((error) => {
    symmetryCache.delete(cacheKey);
    throw error;
  });
  symmetryCache.set(cacheKey, pending);
  return pending;
}

export function multiplyMatrixVector(matrix: Matrix3x3, vector: Vec3): Vec3 {
  return [
    matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
    matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
    matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
  ];
}
