import { useEffect, useMemo, useState } from "react";
import type { QeSourceLayout } from "./qeAutoParallelization";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData, SavedStructureData } from "./types";
import { dot, Matrix3x3, reciprocalLatticeVectors, Vec3 } from "./reciprocalLattice";
import { analyzeCrystalSymmetry, buildSymmetryRequestFromCrystalData, SymmetryAnalyzeRequest, SymmetryTransformResult } from "./symmetryTransform";
import { inferQeBravaisCellFromCif } from "./engines/qe/bravaisInference";

export interface PointEstimate { fullPoints: number; irreduciblePoints: number; }
export interface PointEstimateState { result: PointEstimate | null; loading: boolean; error: string | null; conservative?: boolean; explicit?: boolean; }
export interface QePointEstimates { sourceLayout?: QeSourceLayout; k: PointEstimateState; q?: PointEstimateState; kLabel?: string; }

export function savedPointStructure(saved: SavedStructureData | null): SymmetryAnalyzeRequest | null {
  if (!saved?.cell_parameters || !saved.atoms.length) return null;
  if (saved.cell_units && !["angstrom", "bohr"].includes(saved.cell_units)) return null;
  if (!["crystal", "angstrom", "bohr"].includes(saved.position_units)) return null;
  const factor = saved.cell_units === "bohr" ? 0.529177210903 : 1;
  const lattice = saved.cell_parameters.map((row) => row.map((v) => v * factor)) as Matrix3x3;
  const reciprocal = reciprocalLatticeVectors(lattice);
  return { lattice, atoms: saved.atoms.map((atom) => {
    const position = saved.position_units === "crystal" ? atom.position : reciprocal.map((row) =>
      dot(row, atom.position) * (saved.position_units === "bohr" ? 0.529177210903 : 1) / (2 * Math.PI)) as Vec3;
    return { symbol: atom.symbol, position };
  }) };
}

/** Use the same direct basis as the automatic QE grid, including inferred ibrav. */
export function crystalPointStructure(crystal: CrystalData, symmetry: SymmetryTransformResult | null, primitive: boolean): SymmetryAnalyzeRequest | null {
  if (!primitive) return buildSymmetryRequestFromCrystalData(crystal);
  if (!symmetry) return null;
  const inferred = inferQeBravaisCellFromCif(crystal, symmetry);
  const primitiveLattice = symmetry.standardizedPrimitiveLattice;
  const lattice = inferred ? inferred.primitiveToQeReciprocal.map((row) =>
    [0, 1, 2].map((j) => row.reduce((sum, value, i) => sum + value * primitiveLattice[i][j], 0))) as Matrix3x3 : primitiveLattice;
  return { lattice, atoms: (inferred?.atoms ?? symmetry.standardizedPrimitiveAtoms).map(({ symbol, position }) => ({ symbol, position })) };
}

export function usePointStructure(crystal: CrystalData | null, symmetry: SymmetryTransformResult | null, primitive: boolean, saved: SavedStructureData | null = null) {
  const [resolved, setResolved] = useState<{ key: string; symmetry: SymmetryTransformResult } | null>(null);
  const key = JSON.stringify(crystal);
  useEffect(() => {
    if (!crystal || !primitive || symmetry) return;
    let cancelled = false;
    void analyzeCrystalSymmetry(crystal).then((value) => {
      if (!cancelled) setResolved({ key, symmetry: value });
    }).catch(() => { /* Leave the count unavailable if the matching cell cannot be resolved. */ });
    return () => { cancelled = true; };
  }, [key, primitive, symmetry]);
  return useMemo(() => {
    if (!primitive && saved) return savedPointStructure(saved);
    if (!crystal) return null;
    return crystalPointStructure(crystal, symmetry ?? (resolved?.key === key ? resolved.symmetry : null), primitive);
  }, [key, symmetry, primitive, JSON.stringify(saved), resolved]);
}

export function usePointMeshEstimate(structure: SymmetryAnalyzeRequest | null, mesh: number[] | null,
  shift: number[] = [0, 0, 0], conservative = false): PointEstimateState {
  const request = structure && mesh ? { structure, mesh, shift, timeReversal: !conservative, useSymmetry: !conservative } : null;
  const key = JSON.stringify(request);
  const [state, setState] = useState<{ key: string; value: PointEstimateState } | null>(null);
  useEffect(() => {
    if (!request) return;
    let cancelled = false;
    const timer = window.setTimeout(() => {
      void invoke<PointEstimate>("estimate_point_mesh", { input: request }).then((result) => {
        if (!cancelled) setState({ key, value: { result, loading: false, error: null, conservative } });
      }).catch((error) => {
        if (!cancelled) setState({ key, value: { result: null, loading: false, error: String(error), conservative } });
      });
    }, 250);
    return () => { cancelled = true; window.clearTimeout(timer); };
  }, [key]);
  return state?.key === key ? state.value : { result: null, loading: !!request, error: request ? null : "Matching calculation cell is unavailable.", conservative };
}

export function balancedPoolCounts(ranks: number, images: number, points: number): number[] {
  if (!Number.isSafeInteger(ranks) || !Number.isSafeInteger(images) || ranks < 1 || images < 1 || ranks % images) return [];
  const perImage = ranks / images;
  const counts = new Set<number>();
  for (let divisor = 1; divisor * divisor <= perImage; divisor++) {
    if (perImage % divisor) continue;
    for (const candidate of [divisor, perImage / divisor]) {
      if (candidate <= points && points % candidate === 0) counts.add(candidate);
    }
  }
  return [...counts].sort((a, b) => a - b);
}

/** QE crystal_b includes a final endpoint and one sample at each path break. */
export function bandPathPointCount(path: { npoints: number }[]): number {
  return path.length ? path.slice(0, -1).reduce((sum, point) => sum + Math.max(1, point.npoints), 1) : 0;
}
