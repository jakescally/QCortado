import { CrystalData } from "./types";
import { detectBravaisLattice, BravaisLattice } from "./brillouinZone";
import {
  CenteringType,
  RhombohedralSetting,
  Vec3,
  detectRhombohedralSettingFromLattice,
  kPointPrimitiveToConventional,
} from "./reciprocalLattice";
import {
  SymmetryTransformResult,
  buildConventionalLatticeFromCrystalData,
  multiplyMatrixVector,
} from "./symmetryTransform";

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

function normalizeSpaceGroupHM(value: string | undefined): string {
  if (!value) return "";
  return value
    .toLowerCase()
    .replace(/[−–—]/g, "-")
    .replace(/[\s_:'"]/g, "");
}

function detectBravaisFromHMSymbol(value: string | undefined): BravaisLattice | null {
  const normalized = normalizeSpaceGroupHM(value);
  if (!normalized) return null;
  if (normalized.startsWith("r")) {
    return "trigonal-R";
  }
  if (
    normalized.startsWith("p6") ||
    normalized.startsWith("p-6") ||
    normalized.startsWith("p63")
  ) {
    return "hexagonal";
  }
  return null;
}

function crystalSystemRankFromSpaceGroup(spaceGroupNumber: number): number {
  if (spaceGroupNumber >= 1 && spaceGroupNumber <= 2) return 0; // triclinic
  if (spaceGroupNumber >= 3 && spaceGroupNumber <= 15) return 1; // monoclinic
  if (spaceGroupNumber >= 16 && spaceGroupNumber <= 74) return 2; // orthorhombic
  if (spaceGroupNumber >= 75 && spaceGroupNumber <= 142) return 3; // tetragonal
  if (spaceGroupNumber >= 143 && spaceGroupNumber <= 167) return 4; // trigonal
  if (spaceGroupNumber >= 168 && spaceGroupNumber <= 194) return 5; // hexagonal
  if (spaceGroupNumber >= 195 && spaceGroupNumber <= 230) return 6; // cubic
  return -1;
}

export interface PathTransformContext {
  centering: CenteringType;
  rhombohedralSetting?: RhombohedralSetting;
}

export function resolvePathTransformContext(
  crystalData: CrystalData,
  symmetryTransform: SymmetryTransformResult | null,
): PathTransformContext {
  const cifSpaceGroup =
    coerceSpaceGroupNumber(crystalData.space_group_IT_number) ??
    extractSpaceGroupNumberFromHM(crystalData.space_group_HM);
  const symmetrySpaceGroup = coerceSpaceGroupNumber(symmetryTransform?.spacegroupNumber);

  const preferCifSpaceGroup =
    cifSpaceGroup != null &&
    symmetrySpaceGroup != null &&
    crystalSystemRankFromSpaceGroup(symmetrySpaceGroup) <
      crystalSystemRankFromSpaceGroup(cifSpaceGroup);

  const spaceGroup = preferCifSpaceGroup
    ? cifSpaceGroup
    : (symmetrySpaceGroup ?? cifSpaceGroup);

  const hmFallback = preferCifSpaceGroup
    ? (detectBravaisFromHMSymbol(crystalData.space_group_HM) ??
      detectBravaisFromHMSymbol(symmetryTransform?.internationalSymbol))
    : (detectBravaisFromHMSymbol(symmetryTransform?.internationalSymbol) ??
      detectBravaisFromHMSymbol(crystalData.space_group_HM));

  let bravaisType: BravaisLattice = hmFallback ?? "triclinic";
  if (spaceGroup != null) {
    try {
      bravaisType = detectBravaisLattice(spaceGroup);
    } catch {
      if (hmFallback) {
        bravaisType = hmFallback;
      }
    }
  }

  const centeringMap: Record<BravaisLattice, CenteringType> = {
    "cubic-P": "P",
    "cubic-F": "F",
    "cubic-I": "I",
    "tetragonal-P": "P",
    "tetragonal-I": "I",
    "orthorhombic-P": "P",
    "orthorhombic-C": "C",
    "orthorhombic-I": "I",
    "orthorhombic-F": "F",
    "hexagonal": "P",
    "trigonal-R": "R",
    "monoclinic-P": "P",
    "monoclinic-C": "C",
    "triclinic": "P",
  };

  const centering = centeringMap[bravaisType] ?? "P";
  const rhombohedralSetting = centering === "R"
    ? detectRhombohedralSettingFromLattice(buildConventionalLatticeFromCrystalData(crystalData))
    : undefined;

  return {
    centering,
    rhombohedralSetting,
  };
}

function roundCoord(value: number): number {
  const rounded = Math.abs(value) < 1e-12 ? 0 : Number(value.toFixed(12));
  return Math.abs(rounded) < 1e-12 ? 0 : rounded;
}

export function roundVec3(coords: Vec3): Vec3 {
  return [
    roundCoord(coords[0]),
    roundCoord(coords[1]),
    roundCoord(coords[2]),
  ];
}

export interface PathCoordinateConverters {
  toInputConventionalCoords: (coords: Vec3) => Vec3;
  toSymmetryPrimitiveCoords: (coords: Vec3) => Vec3;
}

export function createPathCoordinateConverters(
  context: PathTransformContext,
  symmetryTransform: SymmetryTransformResult | null,
): PathCoordinateConverters {
  const toInputConventionalCoords = (coords: Vec3): Vec3 => {
    return roundVec3(
      kPointPrimitiveToConventional(
        coords,
        context.centering,
        context.centering === "R"
          ? { rhombohedralSetting: context.rhombohedralSetting }
          : undefined,
      ),
    );
  };

  const toSymmetryPrimitiveCoords = (coords: Vec3): Vec3 => {
    const inputCoords = toInputConventionalCoords(coords);
    if (!symmetryTransform) {
      return inputCoords;
    }
    return roundVec3(
      multiplyMatrixVector(
        symmetryTransform.inputToPrimitiveReciprocal,
        inputCoords,
      ),
    );
  };

  return {
    toInputConventionalCoords,
    toSymmetryPrimitiveCoords,
  };
}

export function sourceScfUsesPrimitiveCell(scfParameters: unknown): boolean {
  const params = (scfParameters ?? {}) as Record<string, unknown>;
  const sourceCellRepresentation = String(params.cell_representation || "").toLowerCase();
  return sourceCellRepresentation.startsWith("primitive");
}

type Vec3WithCoords = { coords: Vec3 };

export function mapPathCoordinates<T extends Vec3WithCoords>(
  path: T[],
  transformCoords: (coords: Vec3) => Vec3,
): T[] {
  return path.map((point) => ({
    ...point,
    coords: transformCoords(point.coords),
  }));
}
