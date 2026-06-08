import type { CrystalData } from "./types";
import type { KPathPoint } from "../components/BrillouinZoneViewer";
import { getLeadingElementSymbol } from "./elements";
import { resolvePathTransformContext, roundVec3 } from "./kPathTransforms";
import { kPointPrimitiveToConventional } from "./reciprocalLattice";

export type Wien2kBandProjectionKind = "all" | "atom" | "orbital";

export interface Wien2kBandProjectionOption {
  value: string;
  label: string;
  kind: Wien2kBandProjectionKind;
  atomIndex: number;
  characterAtom: number;
  characterL: number;
}

export interface Wien2kBandProjectionSite {
  symbol?: string | null;
  positions?: unknown[] | null;
}

const WIEN2K_ORBITAL_LABELS: Array<{ l: number; label: string }> = [
  { l: 0, label: "s" },
  { l: 1, label: "p" },
  { l: 2, label: "d" },
  { l: 3, label: "f" },
];

function clampInt(value: number, min: number, max: number): number {
  if (!Number.isFinite(value)) return min;
  return Math.min(max, Math.max(min, Math.round(value)));
}

function getConnectedSegmentIndices(path: KPathPoint[]): number[] {
  const indices: number[] = [];
  for (let i = 0; i < path.length - 1; i += 1) {
    if (path[i].npoints > 0) {
      indices.push(i);
    }
  }
  return indices;
}

export function applyWien2kBandsTotalKPoints(path: KPathPoint[], totalKPoints: number): KPathPoint[] {
  const connectedSegmentIndices = getConnectedSegmentIndices(path);
  if (path.length === 0 || connectedSegmentIndices.length === 0) {
    return path.map((point) => ({ ...point, npoints: 0 }));
  }

  const safeTotal = clampInt(totalKPoints, connectedSegmentIndices.length, 5000);
  const remainingAfterBaseline = safeTotal - connectedSegmentIndices.length;
  const lengths = connectedSegmentIndices.map((segmentIndex) => {
    const from = path[segmentIndex];
    const to = path[segmentIndex + 1];
    const dx = to.coords[0] - from.coords[0];
    const dy = to.coords[1] - from.coords[1];
    const dz = to.coords[2] - from.coords[2];
    const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
    return Number.isFinite(distance) && distance > 1e-9 ? distance : 1e-9;
  });
  const totalLength = lengths.reduce((sum, len) => sum + len, 0);
  const rawExtras = lengths.map((length) => (
    totalLength > 0
      ? (length / totalLength) * remainingAfterBaseline
      : remainingAfterBaseline / connectedSegmentIndices.length
  ));
  const extraPoints = rawExtras.map((value) => Math.floor(value));
  let leftovers = remainingAfterBaseline - extraPoints.reduce((sum, value) => sum + value, 0);
  const order = rawExtras
    .map((value, idx) => ({ idx, frac: value - Math.floor(value), len: lengths[idx] }))
    .sort((a, b) => {
      if (b.frac !== a.frac) return b.frac - a.frac;
      if (b.len !== a.len) return b.len - a.len;
      return a.idx - b.idx;
    });

  for (let i = 0; i < order.length && leftovers > 0; i += 1) {
    extraPoints[order[i].idx] += 1;
    leftovers -= 1;
  }

  const segmentPoints = new Map<number, number>();
  for (let i = 0; i < connectedSegmentIndices.length; i += 1) {
    segmentPoints.set(connectedSegmentIndices[i], 1 + extraPoints[i]);
  }

  return path.map((point, index) => {
    if (index >= path.length - 1) {
      return { ...point, npoints: 0 };
    }
    if (!segmentPoints.has(index)) {
      return { ...point, npoints: 0 };
    }
    return { ...point, npoints: segmentPoints.get(index)! };
  });
}

export function transformWien2kKPathForKlistBand(
  path: KPathPoint[],
  crystalData: CrystalData,
): KPathPoint[] {
  const context = resolvePathTransformContext(crystalData, null);

  return path.map((point) => ({
    ...point,
    // WIEN2k klist_band uses conventional reciprocal components for centered
    // non-R lattices, while R band paths are supplied in primitive reciprocal coordinates.
    coords: context.centering === "R"
      ? roundVec3(point.coords)
      : roundVec3(kPointPrimitiveToConventional(point.coords, context.centering)),
  }));
}

function normalizeElementLabel(value: string): string {
  return value.trim().replace(/\s+/g, " ");
}

function projectionElementLabel(site: CrystalData["atom_sites"][number], fallback: string): string {
  return (
    getLeadingElementSymbol(site.type_symbol)
    ?? getLeadingElementSymbol(site.label)
    ?? normalizeElementLabel(site.type_symbol || site.label || fallback)
  );
}

export function getWien2kBandProjectionOptions(crystalData: CrystalData | null): Wien2kBandProjectionOption[] {
  const options: Wien2kBandProjectionOption[] = [
    {
      value: "all",
      label: "All",
      kind: "all",
      atomIndex: 0,
      characterAtom: 0,
      characterL: 0,
    },
  ];

  const elementGroups: Array<{ label: string }> = [];
  const groupIndices = new Map<string, number>();
  for (const [index, site] of (crystalData?.atom_sites ?? []).entries()) {
    const label = projectionElementLabel(site, `Atom ${index + 1}`);
    const key = label.toLowerCase();
    const existingIndex = groupIndices.get(key);
    if (existingIndex != null) {
      continue;
    }
    groupIndices.set(key, elementGroups.length);
    elementGroups.push({ label });
  }

  for (const [index, group] of elementGroups.entries()) {
    const atomIndex = index + 1;
    const atomLabel = group.label;
    const atomValue = `atom:${atomIndex}`;
    options.push({
      value: atomValue,
      label: atomLabel,
      kind: "atom",
      atomIndex,
      characterAtom: atomIndex,
      characterL: 0,
    });
    for (const orbital of WIEN2K_ORBITAL_LABELS) {
      options.push({
        value: `${atomValue}:orbital:${orbital.l}`,
        label: `${atomLabel} ${orbital.label}`,
        kind: "orbital",
        atomIndex,
        characterAtom: atomIndex,
        characterL: orbital.l,
      });
    }
  }

  return options;
}

export function getWien2kBandProjectionOptionsFromSites(
  sites: Wien2kBandProjectionSite[] | null | undefined,
  fallbackCrystalData: CrystalData | null = null,
): Wien2kBandProjectionOption[] {
  if (!Array.isArray(sites) || sites.length === 0) {
    return getWien2kBandProjectionOptions(fallbackCrystalData);
  }

  const options = getWien2kBandProjectionOptions(null);
  for (const [index, site] of sites.entries()) {
    const atomIndex = index + 1;
    const atomLabel = (
      getLeadingElementSymbol(site.symbol ?? "")
      ?? normalizeElementLabel(site.symbol || `Atom ${atomIndex}`)
    );
    const atomValue = `jatom:${atomIndex}`;
    options.push({
      value: atomValue,
      label: atomLabel,
      kind: "atom",
      atomIndex,
      characterAtom: atomIndex,
      characterL: 0,
    });
    for (const orbital of WIEN2K_ORBITAL_LABELS) {
      options.push({
        value: `${atomValue}:orbital:${orbital.l}`,
        label: `${atomLabel} ${orbital.label}`,
        kind: "orbital",
        atomIndex,
        characterAtom: atomIndex,
        characterL: orbital.l,
      });
    }
  }

  return options;
}

export function resolveWien2kBandProjectionOption(
  options: Wien2kBandProjectionOption[],
  selection: string,
): Wien2kBandProjectionOption {
  return options.find((option) => option.value === selection) ?? options[0] ?? {
    value: "all",
    label: "All",
    kind: "all",
    atomIndex: 0,
    characterAtom: 0,
    characterL: 0,
  };
}
