// Project Dashboard - Main view for working with a project's structures and calculations

import { useState, useEffect, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { parseCIF } from "../lib/cifParser";
import { CrystalData, SCFPreset, OptimizedStructureOption, SavedCellSummary } from "../lib/types";
import { getPrimitiveCell } from "../lib/primitiveCell";
import { getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { isPhononReadyScf } from "../lib/phononReady";
import { sourceScfUsesPrimitiveCell } from "../lib/kPathTransforms";
import { extractOptimizedStructure, isSavedStructureData, summarizeCell } from "../lib/optimizedStructure";
import { downloadHpcCalculationArtifacts } from "../lib/hpcConfig";
import { detectBravaisLattice } from "../lib/brillouinZone";
import { detectRhombohedralSettingFromLattice } from "../lib/reciprocalLattice";
import type { BravaisLattice } from "../lib/brillouinZone";
import type { CenteringType, RhombohedralSetting } from "../lib/reciprocalLattice";
import { buildConventionalLatticeFromCrystalData } from "../lib/symmetryTransform";
import { formatWannierConvergenceFlag, getWannierIssueCounts, getWannierQualityIssues } from "../lib/wannierQuality";
import { EditProjectDialog } from "./EditProjectDialog";
import type { TransportResult } from "../lib/transport";

interface QEResult {
  converged: boolean;
  total_energy: number | null;
  fermi_energy: number | null;
  n_scf_steps: number | null;
  wall_time_seconds: number | null;
  raw_output: string;
  band_data?: any;  // Band structure data for bands calculations
  dos_data?: any;  // Electronic DOS data for DOS calculations
  phonon_data?: any;  // Phonon data (DOS + dispersion) for phonon calculations
  wannier_data?: any;  // Wannier90 data payload for Wannier calculations
  transport_data?: TransportResult;
}

export interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: QEResult | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
  storage_bytes?: number | null;
}

export interface WannierBandOverlayOption {
  id: string;
  label: string;
  data: any;
  fermiEnergy: number | null;
}

interface CifVariant {
  id: string;
  filename: string;
  formula: string;
  added_at: string;
  calculations: CalculationRun[];
}

interface Project {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  cif_variants: CifVariant[];
  last_opened_cif_id: string | null;
}

interface ProjectDashboardProps {
  projectId: string;
  readOnly?: boolean;
  refreshToken?: number;
  onBack: () => void;
  onDeleted: () => void;
  onRunSCF: (
    cifId: string,
    crystalData: CrystalData,
    cifContent: string,
    filename: string,
    preset?: SCFPreset,
    presetLock?: boolean,
    optimizedStructures?: OptimizedStructureOption[],
  ) => void;
  onRunBands: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewBands: (bandData: any, scfFermiEnergy: number | null) => void;
  onRunDos: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewDos: (dosData: any, scfFermiEnergy: number | null) => void;
  onRunWannier: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewWannier: (
    wannierData: any,
    scfFermiEnergy: number | null,
    overlayOptions?: WannierBandOverlayOption[],
  ) => void;
  onRunTransport: (cifId: string, crystalData: CrystalData, wannierCalculations: CalculationRun[]) => void;
  onViewTransport: (transportData: TransportResult) => void;
  onRunFermiSurface: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onRunPhonons: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewPhonons: (phononData: any, viewMode: "bands" | "dos") => void;
}

type CalcTagType = "info" | "feature" | "special" | "geometry";
type CellViewMode = "conventional" | "primitive";
type CalculationSortMode = "recent" | "best";
type CalculationCategory = "scf" | "bands" | "dos" | "wannier" | "transport" | "fermi_surface" | "phonon" | "optimization";
type CalculationRuntimeKind = "wall" | "cpu";

interface CalculationRuntimeDisplay {
  kind: CalculationRuntimeKind;
  seconds: number;
}

type HpcDownloadPhase = "connecting" | "collecting" | "downloading" | "saved" | "error";
type LudwigExportMode = "strict_2d" | "quasi_2d_fixed_slice";
type LudwigChemicalPotentialMode = "saved" | "override";

interface HpcDownloadProgress {
  phase: HpcDownloadPhase;
  processedFiles: number;
  totalFiles: number;
  downloadedBytes: number;
  totalBytes: number;
  skippedFiles: number;
}

interface LudwigExportResult {
  bundle_path: string;
  band_count: number;
  chemical_potential_ev: number;
  grid_shape: [number, number];
}

interface DisplayCellMetrics {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

type CellMatrix = [[number, number, number], [number, number, number], [number, number, number]];

const DELETE_CONFIRM_TEXT = "DELETE";
const SOC_PRIORITY_BOOST = 250;
const PINNED_TAG = "pinned";
const RECIPROCAL_AXIS_OPTIONS = [
  { value: 0, label: "a*" },
  { value: 1, label: "b*" },
  { value: 2, label: "c*" },
] as const;

const DISPLAY_CENTERING_BY_BRAVAIS: Record<BravaisLattice, CenteringType> = {
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
  triclinic: "P",
};

interface OptimizationDisplayCellContext {
  centering: CenteringType;
  rhombohedralSetting?: RhombohedralSetting;
}

function findSliceAxis(primaryAxis: number, secondaryAxis: number): number | null {
  const selected = new Set([primaryAxis, secondaryAxis]);
  const remaining = RECIPROCAL_AXIS_OPTIONS.find((option) => !selected.has(option.value));
  return remaining ? remaining.value : null;
}

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

function normalizeSpaceGroupHM(value: string | undefined): string {
  if (!value) return "";
  return value
    .toLowerCase()
    .replace(/[−–—]/g, "-")
    .replace(/[\s_:'"]/g, "");
}

function resolveOptimizationDisplayCellContext(
  calc: CalculationRun,
  crystalData: CrystalData | null,
): OptimizationDisplayCellContext | null {
  const storedSpaceGroup = coerceSpaceGroupNumber(calc.parameters?.symmetry_spacegroup);
  const cifSpaceGroup = coerceSpaceGroupNumber(crystalData?.space_group_IT_number);
  let centering: CenteringType | null = null;

  const spaceGroup = storedSpaceGroup ?? cifSpaceGroup;
  if (spaceGroup != null) {
    centering = DISPLAY_CENTERING_BY_BRAVAIS[detectBravaisLattice(spaceGroup)] ?? "P";
  } else {
    const normalizedHM = normalizeSpaceGroupHM(crystalData?.space_group_HM);
    if (normalizedHM.startsWith("r")) {
      centering = "R";
    }
  }

  if (!centering) return null;

  return {
    centering,
    rhombohedralSetting:
      centering === "R" && crystalData
        ? detectRhombohedralSettingFromLattice(buildConventionalLatticeFromCrystalData(crystalData))
        : undefined,
  };
}

function isHpcCalculation(calc: CalculationRun): boolean {
  const params = calc.parameters || {};
  const backend = String(params.execution_backend || "").trim().toLowerCase();
  if (backend === "hpc") {
    return true;
  }
  if (params.remote_job_id || params.remote_workdir || params.remote_project_path) {
    return true;
  }
  const rawOutput = typeof calc.result?.raw_output === "string" ? calc.result.raw_output : "";
  return rawOutput.includes("HPC_STAGE|") || rawOutput.includes("HPC_CMD|");
}

function getCalcTagClass(tag: { label: string; type: string }): string {
  const normalizedLabel = tag.label.trim().toUpperCase();
  const isHpcTag = normalizedLabel === "HPC";
  const isDownloadedTag = normalizedLabel === "DOWNLOADED";
  return `calc-tag calc-tag-${tag.type}${isHpcTag ? " calc-tag-hpc" : ""}${isDownloadedTag ? " calc-tag-downloaded" : ""}`;
}

function asNonNegativeInteger(value: unknown): number | null {
  const numeric = Number(value);
  if (!Number.isFinite(numeric) || numeric < 0) return null;
  return Math.floor(numeric);
}

function getLocalStorageBytes(calc: CalculationRun): number | null {
  const fromStorageField = asNonNegativeInteger(calc.storage_bytes);
  if (fromStorageField != null) return fromStorageField;
  return asNonNegativeInteger(calc.parameters?.local_storage_bytes);
}

function getRemoteStorageBytes(calc: CalculationRun): number | null {
  const params = calc.parameters || {};
  const explicitRemote = asNonNegativeInteger(params.remote_storage_bytes);
  if (explicitRemote != null) return explicitRemote;
  const fallbackRemote = asNonNegativeInteger(params.remote_artifact_bytes ?? params.remote_total_bytes);
  if (fallbackRemote != null) return fallbackRemote;
  const fermiPayloadBytes = asNonNegativeInteger(params.total_frmsf_bytes ?? params.total_bxsf_bytes);
  if (fermiPayloadBytes != null) return fermiPayloadBytes;
  return null;
}

function isHpcArtifactsDownloaded(calc: CalculationRun): boolean {
  const params = calc.parameters || {};
  if (params.artifacts_downloaded_full === true) return true;
  const syncMode = String(params.artifact_sync_mode || "").trim().toLowerCase();
  return syncMode === "full";
}

function isWannierReadyScf(calc: CalculationRun): boolean {
  if (calc.calc_type !== "scf" || !calc.result?.converged) {
    return false;
  }
  const params = calc.parameters || {};
  if (!sourceScfUsesPrimitiveCell(params)) {
    return false;
  }
  if (Number(params.nspin) !== 1) {
    return false;
  }
  if (Boolean(params.noncolin) || Boolean(params.lspinorb)) {
    return false;
  }
  if (Boolean(params.lda_plus_u)) {
    return false;
  }
  const vdw = String(params.vdw_corr || "").trim().toLowerCase();
  return !vdw || vdw === "none";
}

// Helper to generate calculation feature tags from parameters
function getCalculationTags(calc: CalculationRun): { label: string; type: CalcTagType }[] {
  const tags: { label: string; type: CalcTagType }[] = [];
  const params = calc.parameters || {};
  const phononReady = calc.calc_type === "scf" && isPhononReadyScf(params, calc.tags);
  const wannierReady = isWannierReadyScf(calc);
  let hasGeometryTag = false;
  const pushTag = (label: string, type: CalcTagType) => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  // Special tags from stored tags array (phonon-ready, structure-optimized).
  if (calc.tags) {
    for (const tag of calc.tags) {
      if (tag === "phonon-ready") {
        if (phononReady) {
          pushTag("Phonon-Ready", "special");
        }
      } else if (tag === "structure-optimized") {
        pushTag("Optimized", "special");
      } else if (tag === "geometry") {
        pushTag("Geometry", "geometry");
        hasGeometryTag = true;
      }
    }
  }

  if (!hasGeometryTag && params.structure_source?.type === "optimization") {
    pushTag("Geometry", "geometry");
  }

  if (phononReady) {
    pushTag("Phonon-Ready", "special");
  }
  if (wannierReady) {
    pushTag("Wannier-Ready", "special");
  }

  // K-points grid
  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }

  // Convergence threshold
  if (params.conv_thr) {
    const thr = params.conv_thr;
    // Format as scientific notation if small
    const label = thr < 0.001 ? thr.toExponential(0) : thr.toString();
    pushTag(label, "info");
  }

  // Feature tags
  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }

  if (params.nspin === 4) {
    pushTag("Non-collinear", "feature");
  } else if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }

  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }

  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
    if (isHpcArtifactsDownloaded(calc)) {
      pushTag("Downloaded", "feature");
    }
  }

  return tags;
}

function isOptimizationCalculation(calc: CalculationRun): boolean {
  return calc.calc_type === "optimization" || calc.calc_type === "relax" || calc.calc_type === "vcrelax";
}

function getOptimizationMode(calc: CalculationRun): "relax" | "vcrelax" {
  const mode = calc.parameters?.optimization_mode;
  if (mode === "relax" || calc.calc_type === "relax") return "relax";
  if (mode === "vcrelax" || calc.calc_type === "vcrelax") return "vcrelax";
  return "vcrelax";
}

function getOptimizationTags(calc: CalculationRun): { label: string; type: CalcTagType }[] {
  const tags: { label: string; type: CalcTagType }[] = [];
  const params = calc.parameters || {};
  tags.push({ label: "Geometry", type: "geometry" });
  tags.push({ label: getOptimizationMode(calc) === "vcrelax" ? "VC-Relax" : "Relax", type: "info" });

  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    tags.push({ label: `${k1}×${k2}×${k3}`, type: "info" });
  }

  const forceConv = formatThreshold(params.forc_conv_thr);
  if (forceConv) {
    tags.push({ label: `F ${forceConv}`, type: "info" });
  }

  const energyConv = formatThreshold(params.etot_conv_thr);
  if (energyConv) {
    tags.push({ label: `E ${energyConv}`, type: "info" });
  }

  if (isHpcCalculation(calc)) {
    tags.push({ label: "HPC", type: "feature" });
  }

  return tags;
}

// Helper to generate bands-specific tags
function getBandsTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  // K-points info
  if (params.total_k_points) {
    pushTag(`${params.total_k_points} k-pts`, "info");
  }

  // Inherited feature tags from SCF
  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }
  if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }
  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  // Orbital projection/fat-band marker
  const hasProjectionTag = Array.isArray(calc.tags)
    && calc.tags.some((tag) => {
      const normalized = String(tag).trim().toLowerCase();
      return normalized === "proj" || normalized === "orb";
    });
  if (hasProjectionTag || params.fat_bands_requested) {
    // Backward compatibility: legacy entries may still use the older `orb` marker.
    pushTag("Proj", "feature");
  }

  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
  }

  return tags;
}

// Helper to generate electronic DOS-specific tags
function getDosTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.dos_k_grid) {
    const [k1, k2, k3] = params.dos_k_grid;
    pushTag(`${k1}×${k2}×${k3} K`, "info");
  }

  if (params.n_points) {
    pushTag(`${params.n_points} pts`, "info");
  }

  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }
  if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }
  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
  }

  return tags;
}

function getWannierTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};
  const issueCounts = getWannierIssueCounts(
    calc.result?.wannier_data ?? null,
    calc.result?.raw_output ?? null,
    calc.result?.fermi_energy ?? null,
  );
  const pushTag = (label: string, type: "info" | "feature") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.k_grid) {
    const [k1, k2, k3] = params.k_grid;
    pushTag(`${k1}×${k2}×${k3} K`, "info");
  }
  if (params.n_wann) {
    pushTag(`${params.n_wann} WF`, "info");
  }
  if (params.n_bands) {
    pushTag(`${params.n_bands} bands`, "info");
  }
  if (Array.isArray(params.exclude_bands) && params.exclude_bands.length > 0) {
    pushTag("Disentangled", "feature");
  }
  if (params.total_spread != null && Number.isFinite(Number(params.total_spread))) {
    pushTag(`Ω ${Number(params.total_spread).toFixed(3)}`, "info");
  }
  if (issueCounts.errors > 0) {
    pushTag("Needs Review", "feature");
  } else if (issueCounts.warnings > 0) {
    pushTag("Warning", "feature");
  }
  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
  }

  return tags;
}

function getTransportTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.boltz_kmesh) {
    const [k1, k2, k3] = params.boltz_kmesh;
    pushTag(`${k1}×${k2}×${k3} K`, "info");
  }
  if (params.mu_points) {
    pushTag(`${params.mu_points} μ`, "info");
  }
  if (params.temperature_points) {
    pushTag(`${params.temperature_points} T`, "info");
  }
  if (Number.isFinite(Number(params.relaxation_time_fs))) {
    pushTag(`τ ${Number(params.relaxation_time_fs).toFixed(1)} fs`, "info");
  }
  if (params.engine) {
    pushTag(String(params.engine), "feature");
  }
  if (params.is_2d) {
    pushTag("2D", "feature");
  }
  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
  }

  return tags;
}

function getFermiSurfaceTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.fermi_k_grid) {
    const [k1, k2, k3] = params.fermi_k_grid;
    pushTag(`${k1}×${k2}×${k3} K`, "info");
  }

  if (params.n_frmsf_files) {
    pushTag(`${params.n_frmsf_files} FRMSF`, "info");
  } else if (params.n_bxsf_files) {
    pushTag(`${params.n_bxsf_files} BXSF`, "info");
  }

  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }
  if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }
  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
  }

  return tags;
}

// Helper to generate phonon-specific tags
function getPhononTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};

  if (params.q_grid) {
    const [q1, q2, q3] = params.q_grid;
    tags.push({ label: `${q1}×${q2}×${q3} Q`, type: "info" });
  }

  if (params.n_modes) {
    tags.push({ label: `${params.n_modes} modes`, type: "info" });
  }

  if (params.calculate_dos) {
    tags.push({ label: "DOS", type: "feature" });
  }

  if (params.calculate_dispersion) {
    tags.push({ label: "Dispersion", type: "feature" });
  }

  if (isHpcCalculation(calc)) {
    tags.push({ label: "HPC", type: "feature" });
  }

  return tags;
}

function getRecencyTimestamp(calc: CalculationRun): number {
  const completed = calc.completed_at ? Date.parse(calc.completed_at) : Number.NaN;
  if (Number.isFinite(completed)) return completed;
  const started = Date.parse(calc.started_at);
  if (Number.isFinite(started)) return started;
  return 0;
}

function getMeshProduct(mesh: unknown): number {
  if (!Array.isArray(mesh) || mesh.length !== 3) return 0;
  const values = mesh.map((entry) => Number(entry));
  if (!values.every((value) => Number.isFinite(value) && value > 0)) return 0;
  return values[0] * values[1] * values[2];
}

function getThresholdTightness(value: unknown, maxScore = 20): number {
  const numeric = Number(value);
  if (!Number.isFinite(numeric) || numeric <= 0) return 0;
  return Math.max(0, Math.min(maxScore, -Math.log10(numeric)));
}

function formatThreshold(value: unknown): string | null {
  const numeric = Number(value);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  return numeric < 0.001 ? numeric.toExponential(1) : numeric.toString();
}

function parsePositiveInt(value: string | undefined): number {
  const parsed = Number(value ?? "");
  if (!Number.isFinite(parsed) || parsed < 0) return 0;
  return Math.floor(parsed);
}

function defaultHpcDownloadProgress(phase: HpcDownloadPhase): HpcDownloadProgress {
  return {
    phase,
    processedFiles: 0,
    totalFiles: 0,
    downloadedBytes: 0,
    totalBytes: 0,
    skippedFiles: 0,
  };
}

function getHpcDownloadPercent(progress: HpcDownloadProgress): number | null {
  if (progress.totalBytes > 0) {
    return Math.min(
      100,
      Math.max(0, Math.round((progress.downloadedBytes / progress.totalBytes) * 100)),
    );
  }
  if (progress.totalFiles > 0) {
    return Math.min(
      100,
      Math.max(0, Math.round((progress.processedFiles / progress.totalFiles) * 100)),
    );
  }
  return null;
}

function parseQeTimeToken(token: string): number | null {
  const value = token.trim();
  if (!value) return null;

  let totalSeconds = 0;
  let matched = false;

  const hourMatch = value.match(/([\d.]+)h/);
  if (hourMatch) {
    const hours = Number(hourMatch[1]);
    if (!Number.isFinite(hours)) return null;
    totalSeconds += hours * 3600;
    matched = true;
  }

  const minuteMatch = value.match(/([\d.]+)m/);
  if (minuteMatch) {
    const minutes = Number(minuteMatch[1]);
    if (!Number.isFinite(minutes)) return null;
    totalSeconds += minutes * 60;
    matched = true;
  }

  const secondMatch = value.match(/([\d.]+)s/);
  if (secondMatch) {
    const seconds = Number(secondMatch[1]);
    if (!Number.isFinite(seconds)) return null;
    totalSeconds += seconds;
    matched = true;
  }

  if (matched) {
    return totalSeconds;
  }

  const fallbackSeconds = Number(value);
  return Number.isFinite(fallbackSeconds) ? fallbackSeconds : null;
}

function parseRuntimeFromRawOutput(rawOutput: string): CalculationRuntimeDisplay | null {
  let totalCpuSeconds = 0;
  let totalWallSeconds = 0;
  let foundTimingLine = false;

  const timingPattern = /^\s*[A-Za-z0-9_.+-]+\s*:\s*([0-9.hms]+)\s+CPU(?:\s+([0-9.hms]+)\s+WALL)?/gim;
  for (const match of rawOutput.matchAll(timingPattern)) {
    foundTimingLine = true;
    const cpuSeconds = parseQeTimeToken(match[1]);
    const wallSeconds = match[2] ? parseQeTimeToken(match[2]) : null;

    if (typeof cpuSeconds === "number" && cpuSeconds > 0) {
      totalCpuSeconds += cpuSeconds;
    }
    if (typeof wallSeconds === "number" && wallSeconds > 0) {
      totalWallSeconds += wallSeconds;
    }
  }

  if (!foundTimingLine) return null;
  if (totalWallSeconds > 0) return { kind: "wall", seconds: totalWallSeconds };
  if (totalCpuSeconds > 0) return { kind: "cpu", seconds: totalCpuSeconds };
  return null;
}

function getCalculationBestScore(calc: CalculationRun, category: CalculationCategory): number {
  const params = calc.parameters || {};
  const convergedBonus = calc.result?.converged ? 100 : 0;
  const socBonus = params.lspinorb ? SOC_PRIORITY_BOOST : 0;

  if (category === "scf") {
    const kScore = Math.log2(Math.max(1, getMeshProduct(params.kgrid)));
    const convScore = getThresholdTightness(params.conv_thr);
    const ecutScore = Math.log2(Math.max(1, Number(params.ecutwfc) || 1));
    return convergedBonus + (4 * convScore) + (3 * kScore) + socBonus + ecutScore;
  }

  if (category === "phonon") {
    const qScore = Math.log2(Math.max(1, getMeshProduct(params.q_grid)));
    const tr2Score = getThresholdTightness(params.tr2_ph);
    return convergedBonus + (6 * qScore) + tr2Score;
  }

  if (category === "optimization") {
    const kScore = Math.log2(Math.max(1, getMeshProduct(params.kgrid)));
    const convScore = getThresholdTightness(params.conv_thr);
    const forceScore = getThresholdTightness(params.forc_conv_thr);
    const energyScore = getThresholdTightness(params.etot_conv_thr);
    return convergedBonus + (2 * kScore) + (2 * convScore) + (4 * forceScore) + (3 * energyScore) + socBonus;
  }

  if (category === "fermi_surface") {
    const kScore = Math.log2(Math.max(1, getMeshProduct(params.fermi_k_grid)));
    const bandScore = Math.log2(Math.max(1, Number(params.fermi_n_bands_requested) || 1));
    return convergedBonus + (4 * kScore) + bandScore + socBonus;
  }

  if (category === "dos") {
    const kScore = Math.log2(Math.max(1, getMeshProduct(params.dos_k_grid)));
    const resolutionScore = getThresholdTightness(params.dos_delta_e, 12);
    return convergedBonus + (4 * kScore) + (2 * resolutionScore) + socBonus;
  }

  if (category === "wannier") {
    const meshScore = Math.log2(Math.max(1, getMeshProduct(params.k_grid)));
    const wannScore = Math.log2(Math.max(1, Number(params.n_wann) || 1));
    const bandScore = Math.log2(Math.max(1, Number(params.n_bands) || 1));
    const interpolationScore = Math.log2(Math.max(1, Number(params.total_k_points) || 1));
    return convergedBonus + (4 * meshScore) + wannScore + bandScore + interpolationScore + socBonus;
  }

  if (category === "transport") {
    const meshScore = Math.log2(Math.max(1, getMeshProduct(params.boltz_kmesh)));
    const muScore = Math.log2(Math.max(1, Number(params.mu_points) || 1));
    const tempScore = Math.log2(Math.max(1, Number(params.temperature_points) || 1));
    return convergedBonus + (4 * meshScore) + (2 * muScore) + (2 * tempScore);
  }

  // Bands: prioritize denser path sampling, then inherited SCF settings when present.
  const pathScore = Math.log2(Math.max(1, Number(params.total_k_points) || 0));
  const bandCountScore = Math.log2(Math.max(1, Number(params.n_bands) || 0));
  return convergedBonus + (3 * pathScore) + bandCountScore + socBonus;
}

function sortCalculations(
  calculations: CalculationRun[],
  sortMode: CalculationSortMode,
  category: CalculationCategory,
  pinnedIds: Set<string>,
): CalculationRun[] {
  const sorted = [...calculations];
  sorted.sort((a, b) => {
    const aPinned = pinnedIds.has(a.id);
    const bPinned = pinnedIds.has(b.id);
    if (aPinned !== bPinned) {
      return aPinned ? -1 : 1;
    }

    if (sortMode === "best") {
      const diff = getCalculationBestScore(b, category) - getCalculationBestScore(a, category);
      if (Math.abs(diff) > 1e-9) return diff;
    }
    return getRecencyTimestamp(b) - getRecencyTimestamp(a);
  });
  return sorted;
}

function calculateMetricsFromVectors(
  v1: [number, number, number],
  v2: [number, number, number],
  v3: [number, number, number],
): DisplayCellMetrics {
  const norm = (v: [number, number, number]) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  const dot = (u: [number, number, number], v: [number, number, number]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const angle = (u: [number, number, number], v: [number, number, number]) => {
    const denom = norm(u) * norm(v);
    if (denom === 0) return 0;
    const cosine = dot(u, v) / denom;
    const clamped = Math.max(-1, Math.min(1, cosine));
    return (Math.acos(clamped) * 180) / Math.PI;
  };

  return {
    a: norm(v1),
    b: norm(v2),
    c: norm(v3),
    alpha: angle(v2, v3),
    beta: angle(v1, v3),
    gamma: angle(v1, v2),
  };
}

function calculateVolumeFromMetrics(metrics: DisplayCellMetrics): number {
  const alpha = (metrics.alpha * Math.PI) / 180;
  const beta = (metrics.beta * Math.PI) / 180;
  const gamma = (metrics.gamma * Math.PI) / 180;
  const cosAlpha = Math.cos(alpha);
  const cosBeta = Math.cos(beta);
  const cosGamma = Math.cos(gamma);
  const factor = 1 + (2 * cosAlpha * cosBeta * cosGamma)
    - (cosAlpha * cosAlpha)
    - (cosBeta * cosBeta)
    - (cosGamma * cosGamma);
  const safeFactor = Math.max(0, factor);
  return metrics.a * metrics.b * metrics.c * Math.sqrt(safeFactor);
}

function calculateVolumeFromVectors(
  v1: [number, number, number],
  v2: [number, number, number],
  v3: [number, number, number],
): number {
  const cross: [number, number, number] = [
    v2[1] * v3[2] - v2[2] * v3[1],
    v2[2] * v3[0] - v2[0] * v3[2],
    v2[0] * v3[1] - v2[1] * v3[0],
  ];
  return Math.abs(v1[0] * cross[0] + v1[1] * cross[1] + v1[2] * cross[2]);
}

function relativeApproximatelyEqual(a: number, b: number, tolerance = 0.35): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tolerance * scale;
}

function getCenteringMultiplicity(centering: CenteringType): number {
  switch (centering) {
    case "F":
      return 4;
    case "I":
    case "A":
    case "B":
    case "C":
      return 2;
    case "R":
      return 3;
    case "P":
    default:
      return 1;
  }
}

function isLikelyPrimitiveOptimizationCell(
  calc: CalculationRun,
  optimizedCell: CellMatrix,
  displayContext: OptimizationDisplayCellContext | null,
  crystalData: CrystalData | null,
): boolean {
  if (!displayContext || displayContext.centering === "P") {
    return false;
  }

  const cellRepresentation = String(calc.parameters?.cell_representation || "").trim().toLowerCase();
  if (cellRepresentation.startsWith("primitive")) {
    return true;
  }

  if (displayContext.centering === "R" && displayContext.rhombohedralSetting === "hexagonal") {
    if (detectRhombohedralSettingFromLattice(optimizedCell) === "rhombohedral") {
      return true;
    }
  }

  if (!crystalData) {
    return false;
  }

  const conventionalCell = buildConventionalLatticeFromCrystalData(crystalData);
  const conventionalVolume = calculateVolumeFromVectors(
    conventionalCell[0],
    conventionalCell[1],
    conventionalCell[2],
  );
  const optimizedVolume = calculateVolumeFromVectors(
    optimizedCell[0],
    optimizedCell[1],
    optimizedCell[2],
  );
  if (conventionalVolume <= 0 || optimizedVolume <= 0) {
    return false;
  }

  return relativeApproximatelyEqual(
    optimizedVolume * getCenteringMultiplicity(displayContext.centering),
    conventionalVolume,
  );
}

function convertPrimitiveToConventionalCell(
  primitiveCell: CellMatrix,
  centering: CenteringType,
  rhombohedralSetting?: RhombohedralSetting,
): CellMatrix | null {
  const [p1, p2, p3] = primitiveCell;

  switch (centering) {
    case "P":
      return primitiveCell;
    case "F":
      return [
        [-p1[0] + p2[0] + p3[0], -p1[1] + p2[1] + p3[1], -p1[2] + p2[2] + p3[2]],
        [p1[0] - p2[0] + p3[0], p1[1] - p2[1] + p3[1], p1[2] - p2[2] + p3[2]],
        [p1[0] + p2[0] - p3[0], p1[1] + p2[1] - p3[1], p1[2] + p2[2] - p3[2]],
      ];
    case "I":
      return [
        [p2[0] + p3[0], p2[1] + p3[1], p2[2] + p3[2]],
        [p1[0] + p3[0], p1[1] + p3[1], p1[2] + p3[2]],
        [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]],
      ];
    case "C":
      return [
        [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]],
        [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]],
        p3,
      ];
    case "A":
      return [
        p1,
        [p2[0] - p3[0], p2[1] - p3[1], p2[2] - p3[2]],
        [p2[0] + p3[0], p2[1] + p3[1], p2[2] + p3[2]],
      ];
    case "B":
      return [
        [p1[0] - p3[0], p1[1] - p3[1], p1[2] - p3[2]],
        p2,
        [p1[0] + p3[0], p1[1] + p3[1], p1[2] + p3[2]],
      ];
    case "R":
      if (rhombohedralSetting !== "hexagonal") {
        return null;
      }
      return [
        [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]],
        [p2[0] - p3[0], p2[1] - p3[1], p2[2] - p3[2]],
        [p1[0] + p2[0] + p3[0], p1[1] + p2[1] + p3[1], p1[2] + p2[2] + p3[2]],
      ];
    default:
      return null;
  }
}

function asCellMatrix(value: unknown): CellMatrix | null {
  if (!Array.isArray(value) || value.length !== 3) return null;
  const rows: [number, number, number][] = [];
  for (const row of value) {
    if (!Array.isArray(row) || row.length !== 3) return null;
    const x = Number(row[0]);
    const y = Number(row[1]);
    const z = Number(row[2]);
    if (![x, y, z].every((entry) => Number.isFinite(entry))) return null;
    rows.push([x, y, z]);
  }
  return [rows[0], rows[1], rows[2]];
}

export function ProjectDashboard({
  projectId,
  readOnly = false,
  refreshToken,
  onBack,
  onDeleted,
  onRunSCF,
  onRunBands,
  onViewBands,
  onRunDos,
  onViewDos,
  onRunWannier,
  onViewWannier,
  onRunTransport,
  onViewTransport,
  onRunFermiSurface,
  onRunPhonons,
  onViewPhonons,
}: ProjectDashboardProps) {
  const [cellViewMode, setCellViewMode] = useState<CellViewMode>("conventional");
  const [project, setProject] = useState<Project | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // Current CIF selection
  const [selectedCifId, setSelectedCifId] = useState<string | null>(null);
  const [crystalData, setCrystalData] = useState<CrystalData | null>(null);
  const [cifContent, setCifContent] = useState<string>("");

  // Delete project confirmation dialog state
  const [showDeleteDialog, setShowDeleteDialog] = useState(false);
  const [showEditProjectDialog, setShowEditProjectDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeleting, setIsDeleting] = useState(false);

  // Delete calculation confirmation dialog state
  const [showDeleteCalcDialog, setShowDeleteCalcDialog] = useState(false);
  const [calcToDelete, setCalcToDelete] = useState<{ calcId: string; calcType: string } | null>(null);
  const [isDeletingCalc, setIsDeletingCalc] = useState(false);
  const [showLudwigExportDialog, setShowLudwigExportDialog] = useState(false);
  const [calcToExportLudwig, setCalcToExportLudwig] = useState<CalculationRun | null>(null);
  const [isExportingLudwig, setIsExportingLudwig] = useState(false);
  const [ludwigExportMode, setLudwigExportMode] = useState<LudwigExportMode>("strict_2d");
  const [ludwigPrimaryAxis, setLudwigPrimaryAxis] = useState<number>(0);
  const [ludwigSecondaryAxis, setLudwigSecondaryAxis] = useState<number>(1);
  const [ludwigSliceCoordinateInput, setLudwigSliceCoordinateInput] = useState("0.0");
  const [ludwigNkxInput, setLudwigNkxInput] = useState("161");
  const [ludwigNkyInput, setLudwigNkyInput] = useState("161");
  const [ludwigChemicalPotentialMode, setLudwigChemicalPotentialMode] =
    useState<LudwigChemicalPotentialMode>("saved");
  const [ludwigChemicalPotentialInput, setLudwigChemicalPotentialInput] = useState("");

  // Import state
  const [isImporting, setIsImporting] = useState(false);
  const [infoMessage, setInfoMessage] = useState<string | null>(null);
  const [isRefreshingProject, setIsRefreshingProject] = useState(false);
  const [launchingFermiCalcId, setLaunchingFermiCalcId] = useState<string | null>(null);
  const [downloadingCalcId, setDownloadingCalcId] = useState<string | null>(null);
  const [downloadProgressByCalcId, setDownloadProgressByCalcId] = useState<Record<string, HpcDownloadProgress>>({});
  const [calculationDetailsById, setCalculationDetailsById] = useState<Record<string, CalculationRun>>({});

  // Expanded calculation
  const [expandedCalc, setExpandedCalc] = useState<string | null>(null);
  const [calculationSortMode, setCalculationSortMode] = useState<CalculationSortMode>(() => getStoredSortMode());

  useEffect(() => {
    loadProject();
  }, [projectId, refreshToken]);

  useEffect(() => {
    setCalculationDetailsById({});
  }, [projectId]);

  function formatBytes(bytes: number): string {
    if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
    const units = ["B", "KB", "MB", "GB", "TB"];
    let value = bytes;
    let unitIdx = 0;
    while (value >= 1024 && unitIdx < units.length - 1) {
      value /= 1024;
      unitIdx += 1;
    }
    const precision = value >= 10 || unitIdx === 0 ? 0 : 1;
    return `${value.toFixed(precision)} ${units[unitIdx]}`;
  }

  function renderStorageDetailItems(calc: CalculationRun) {
    const localBytes = getLocalStorageBytes(calc);
    const remoteBytes = getRemoteStorageBytes(calc);
    if (localBytes == null && remoteBytes == null) return null;

    return (
      <>
        {localBytes != null && (
          <div className="detail-item">
            <label>Local Size</label>
            <span>{formatBytes(localBytes)}</span>
          </div>
        )}
        {remoteBytes != null && (
          <div className="detail-item">
            <label>Remote Size</label>
            <span>{formatBytes(remoteBytes)}</span>
          </div>
        )}
      </>
    );
  }

  async function loadProject(options: { showLoading?: boolean; refreshSelectedCif?: boolean } = {}): Promise<boolean> {
    const { showLoading = true, refreshSelectedCif = true } = options;
    if (showLoading) {
      setIsLoading(true);
    }
    setError(null);
    try {
      const proj = await invoke<Project>("get_project", { projectId });
      setProject(proj);
      setCalculationDetailsById({});

      // Determine which CIF to show
      if (proj.cif_variants.length > 0) {
        const selectedStillExists =
          selectedCifId !== null && proj.cif_variants.some((variant) => variant.id === selectedCifId);
        const cifToOpen =
          selectedStillExists
            ? selectedCifId
            : proj.last_opened_cif_id &&
                proj.cif_variants.some((variant) => variant.id === proj.last_opened_cif_id)
              ? proj.last_opened_cif_id
              : proj.cif_variants[0].id;

        if (!cifToOpen) return true;

        if (refreshSelectedCif || cifToOpen !== selectedCifId || !crystalData) {
          await selectCif(cifToOpen);
        }
      } else {
        setSelectedCifId(null);
        setCrystalData(null);
        setCifContent("");
        setExpandedCalc(null);
      }
      return true;
    } catch (e) {
      console.error("Failed to load project:", e);
      setError(String(e));
      return false;
    } finally {
      if (showLoading) {
        setIsLoading(false);
      }
    }
  }

  async function handleRefreshProject() {
    if (isRefreshingProject) return;
    setIsRefreshingProject(true);
    setInfoMessage(null);
    const refreshed = await loadProject({ showLoading: false, refreshSelectedCif: false });
    if (refreshed) {
      setInfoMessage("Project refreshed.");
    }
    setIsRefreshingProject(false);
  }

  function getCalculationRecord(calc: CalculationRun): CalculationRun {
    return calculationDetailsById[calc.id] ?? calc;
  }

  async function ensureCalculationDetails(calc: CalculationRun): Promise<CalculationRun> {
    const cached = calculationDetailsById[calc.id];
    if (cached) {
      return cached;
    }

    const detail = await invoke<CalculationRun>("get_project_calculation", {
      projectId,
      calcId: calc.id,
    });
    setCalculationDetailsById((prev) => (
      prev[calc.id]
        ? prev
        : {
          ...prev,
          [calc.id]: detail,
        }
    ));
    return detail;
  }

  function updateCalculationRecords(
    calcId: string,
    updater: (calc: CalculationRun) => CalculationRun,
  ) {
    setProject((prev) => {
      if (!prev) return prev;
      return {
        ...prev,
        cif_variants: prev.cif_variants.map((variant) => ({
          ...variant,
          calculations: variant.calculations.map((calc) => (
            calc.id === calcId ? updater(calc) : calc
          )),
        })),
      };
    });
    setCalculationDetailsById((prev) => {
      const existing = prev[calcId];
      if (!existing) return prev;
      return {
        ...prev,
        [calcId]: updater(existing),
      };
    });
  }

  async function selectCif(cifId: string) {
    setSelectedCifId(cifId);
    setCellViewMode("conventional");
    try {
      // Load crystal data
      const data = await invoke<CrystalData>("get_cif_crystal_data", {
        projectId,
        cifId,
      });
      setCrystalData(data);

      // Load CIF content
      const content = await invoke<string>("get_cif_content", {
        projectId,
        cifId,
      });
      setCifContent(content);

      // Update last opened
      await invoke("set_last_opened_cif", { projectId, cifId });
    } catch (e) {
      console.error("Failed to load CIF data:", e);
      setError(`Failed to load structure: ${e}`);
    }
  }

  async function handleImportCIF() {
    if (readOnly) return;
    try {
      const selected = await open({
        multiple: false,
        filters: [{ name: "CIF Files", extensions: ["cif"] }],
        title: "Select CIF File",
      });

      if (selected && typeof selected === "string") {
        setIsImporting(true);
        setError(null);

        const content = await readTextFile(selected);
        const parsedData = parseCIF(content);
        const filename = selected.split("/").pop() || "structure.cif";
        const formula = parsedData.formula_sum || parsedData.formula_structural || "Unknown";

        const newVariant = await invoke<CifVariant>("add_cif_to_project", {
          projectId,
          cifData: {
            filename,
            formula,
            content,
            crystal_data: parsedData,
          },
        });

        // Reload project and select the new CIF
        await loadProject();
        await selectCif(newVariant.id);
      }
    } catch (e) {
      console.error("Failed to import CIF:", e);
      setError(`Failed to import CIF: ${e}`);
    } finally {
      setIsImporting(false);
    }
  }

  function openEditProjectDialog() {
    if (readOnly) return;
    setShowEditProjectDialog(true);
  }

  function handleProjectMetadataSaved(updatedProject: {
    id: string;
    name: string;
    description: string | null;
  }) {
    setProject((currentProject) => {
      if (!currentProject || currentProject.id !== updatedProject.id) {
        return currentProject;
      }
      return {
        ...currentProject,
        name: updatedProject.name,
        description: updatedProject.description,
      };
    });
    setInfoMessage(`Updated "${updatedProject.name}".`);
    setError(null);
  }

  async function handleConfirmDelete() {
    if (readOnly) return;
    if (deleteConfirmText !== DELETE_CONFIRM_TEXT) return;

    setIsDeleting(true);
    try {
      await invoke("delete_project", { projectId });
      onDeleted();
    } catch (e) {
      console.error("Failed to delete project:", e);
      setError(String(e));
      setIsDeleting(false);
      setShowDeleteDialog(false);
    }
  }

  function openDeleteCalcDialog(calcId: string, calcType: string) {
    setCalcToDelete({ calcId, calcType });
    setShowDeleteCalcDialog(true);
  }

  async function handleConfirmDeleteCalc() {
    if (readOnly) return;
    if (!calcToDelete || !selectedCifId) return;

    setIsDeletingCalc(true);
    const scrollY = window.scrollY;
    try {
      await invoke("delete_calculation", {
        projectId,
        cifId: selectedCifId,
        calcId: calcToDelete.calcId,
      });
      // Refresh the data in-place so the viewport position is preserved.
      await loadProject({ showLoading: false, refreshSelectedCif: false });
      if (expandedCalc === calcToDelete.calcId) {
        setExpandedCalc(null);
      }
      setShowDeleteCalcDialog(false);
      setCalcToDelete(null);
      requestAnimationFrame(() => {
        window.scrollTo({ top: scrollY, left: 0, behavior: "auto" });
      });
    } catch (e) {
      console.error("Failed to delete calculation:", e);
      setError(String(e));
    } finally {
      setIsDeletingCalc(false);
    }
  }

  function closeLudwigExportDialog() {
    if (isExportingLudwig) return;
    setShowLudwigExportDialog(false);
    setCalcToExportLudwig(null);
  }

  function openLudwigExportDialog(calc: CalculationRun) {
    if (readOnly) return;
    const savedChemicalPotential = calc.result?.fermi_energy;
    setCalcToExportLudwig(calc);
    setLudwigExportMode("strict_2d");
    setLudwigPrimaryAxis(0);
    setLudwigSecondaryAxis(1);
    setLudwigSliceCoordinateInput("0.0");
    setLudwigNkxInput("161");
    setLudwigNkyInput("161");
    if (savedChemicalPotential != null && Number.isFinite(savedChemicalPotential)) {
      setLudwigChemicalPotentialMode("saved");
      setLudwigChemicalPotentialInput(savedChemicalPotential.toFixed(6));
    } else {
      setLudwigChemicalPotentialMode("override");
      setLudwigChemicalPotentialInput("0.000000");
    }
    setError(null);
    setInfoMessage(null);
    setShowLudwigExportDialog(true);
  }

  async function handleExportLudwigBundle() {
    if (readOnly || !calcToExportLudwig) return;

    const sliceAxis = findSliceAxis(ludwigPrimaryAxis, ludwigSecondaryAxis);
    if (sliceAxis == null) {
      setError("Choose two distinct in-plane reciprocal axes.");
      return;
    }

    const nkx = Number.parseInt(ludwigNkxInput, 10);
    const nky = Number.parseInt(ludwigNkyInput, 10);
    if (!Number.isFinite(nkx) || nkx < 2 || !Number.isFinite(nky) || nky < 2) {
      setError("Ludwig export requires nkx and nky to both be at least 2.");
      return;
    }

    const sliceCoordinate =
      ludwigExportMode === "strict_2d" ? 0.0 : Number.parseFloat(ludwigSliceCoordinateInput);
    if (!Number.isFinite(sliceCoordinate)) {
      setError("Enter a valid fractional slice coordinate.");
      return;
    }

    let chemicalPotentialEv: number | null = null;
    if (ludwigChemicalPotentialMode === "override") {
      chemicalPotentialEv = Number.parseFloat(ludwigChemicalPotentialInput);
      if (!Number.isFinite(chemicalPotentialEv)) {
        setError("Enter a valid chemical potential in eV.");
        return;
      }
    } else if (!Number.isFinite(calcToExportLudwig.result?.fermi_energy ?? NaN)) {
      setError("This Wannier calculation has no saved Fermi energy. Switch to a manual chemical potential override.");
      return;
    }

    const selectedDestination = await open({
      title: "Select Ludwig Export Folder",
      directory: true,
      multiple: false,
    });
    if (!selectedDestination || Array.isArray(selectedDestination)) {
      return;
    }

    setIsExportingLudwig(true);
    setError(null);
    setInfoMessage(null);
    try {
      const result = await invoke<LudwigExportResult>("export_wannier_for_ludwig", {
        config: {
          projectId,
          calculationId: calcToExportLudwig.id,
          destinationRoot: selectedDestination,
          mode: ludwigExportMode,
          inPlaneAxes: [ludwigPrimaryAxis, ludwigSecondaryAxis],
          sliceAxis,
          sliceCoordinate,
          nkx,
          nky,
          chemicalPotentialEv,
        },
      });
      setShowLudwigExportDialog(false);
      setCalcToExportLudwig(null);
      setInfoMessage(
        `Exported Ludwig bundle (${result.band_count} bands, ${result.grid_shape[0]}×${result.grid_shape[1]} grid) to ${result.bundle_path}.`,
      );
    } catch (e) {
      console.error("Failed to export Ludwig bundle:", e);
      setError(`Failed to export Ludwig bundle: ${e}`);
    } finally {
      setIsExportingLudwig(false);
    }
  }

  async function getOptimizedStructureOptions(calculations: CalculationRun[]): Promise<OptimizedStructureOption[]> {
    const optimizationCalcs = calculations.filter((calc) => isOptimizationCalculation(calc));
    const resolvedCalcs = await Promise.all(optimizationCalcs.map(async (calc) => {
      if (
        isSavedStructureData(calc.parameters?.optimized_structure)
        || isSavedStructureData(calc.parameters?.source_structure)
        || (calc.result?.raw_output && calc.result.raw_output.trim().length > 0)
      ) {
        return calc;
      }
      try {
        return await ensureCalculationDetails(calc);
      } catch (e) {
        console.warn("Failed to load optimization detail:", e);
        return calc;
      }
    }));

    return resolvedCalcs
      .map((calc) => {
        const storedOptimized = isSavedStructureData(calc.parameters?.optimized_structure)
          ? calc.parameters.optimized_structure
          : null;
        const storedSource = isSavedStructureData(calc.parameters?.source_structure)
          ? calc.parameters.source_structure
          : null;
        const parsedStructure = extractOptimizedStructure(
          calc.result?.raw_output || "",
          storedOptimized || storedSource,
        );
        const structure = parsedStructure || storedOptimized;
        if (!structure || !Array.isArray(structure.atoms) || structure.atoms.length === 0) {
          return null;
        }
        const derivedSummary = summarizeCell(structure);

        const mode = getOptimizationMode(calc);
        const completedAt = calc.completed_at;
        const dateLabel = completedAt
          ? new Date(completedAt).toLocaleDateString(undefined, { month: "short", day: "numeric", year: "numeric" })
          : "In progress";

        return {
          calcId: calc.id,
          label: `${mode === "vcrelax" ? "VC-Relax" : "Relax"} (${dateLabel})`,
          mode,
          completedAt,
          structure,
          cellSummary: derivedSummary ?? calc.parameters?.optimized_cell_summary ?? null,
        } as OptimizedStructureOption;
      })
      .filter((opt): opt is OptimizedStructureOption => opt !== null);
  }

  async function handleRunSCF() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    const optimizedStructures = await getOptimizedStructureOptions(variant.calculations);
    onRunSCF(
      selectedCifId,
      crystalData,
      cifContent,
      variant.filename,
      undefined,
      undefined,
      optimizedStructures,
    );
  }

  function handleRunOptimization() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunSCF(selectedCifId, crystalData, cifContent, variant.filename, "relax", true);
  }

  function handleRunBands() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    // Pass all calculations for this CIF - the wizard will filter for SCF
    onRunBands(selectedCifId, crystalData, variant.calculations);
  }

  function handleRunDos() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunDos(selectedCifId, crystalData, variant.calculations);
  }

  function handleRunWannier() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunWannier(selectedCifId, crystalData, variant.calculations);
  }

  function handleRunTransport() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunTransport(selectedCifId, crystalData, variant.calculations);
  }

  function handleRunFermiSurface() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunFermiSurface(selectedCifId, crystalData, variant.calculations);
  }

  async function handleViewBands(calc: CalculationRun) {
    try {
      const detail = await ensureCalculationDetails(calc);
      const bandData = detail.result?.band_data ?? calc.result?.band_data ?? null;
      if (!bandData) {
        setError("Saved band data is unavailable for this calculation.");
        return;
      }
      onViewBands(bandData, detail.result?.fermi_energy ?? calc.result?.fermi_energy ?? null);
    } catch (e) {
      console.error("Failed to load band data:", e);
      setError(`Failed to load band data: ${e}`);
    }
  }

  async function handleViewDos(calc: CalculationRun) {
    try {
      const detail = await ensureCalculationDetails(calc);
      const dosData = detail.result?.dos_data ?? calc.result?.dos_data ?? null;
      if (!dosData) {
        setError("Saved DOS data is unavailable for this calculation.");
        return;
      }
      onViewDos(dosData, detail.result?.fermi_energy ?? calc.result?.fermi_energy ?? null);
    } catch (e) {
      console.error("Failed to load DOS data:", e);
      setError(`Failed to load DOS data: ${e}`);
    }
  }

  async function loadWannierBandOverlayOptions(
    calc: CalculationRun,
    detail: CalculationRun,
  ): Promise<WannierBandOverlayOption[]> {
    const variant = getSelectedVariant();
    if (!variant) {
      return [];
    }

    const sourceScfId = String(detail.parameters?.source_scf_id ?? calc.parameters?.source_scf_id ?? "").trim();
    const kPath = normalizeSavedKPath(detail.parameters?.k_path ?? calc.parameters?.k_path);
    if (!sourceScfId || !kPath) {
      return [];
    }

    const matchingBandRuns = variant.calculations.filter((candidate) => {
      if (candidate.calc_type !== "bands") {
        return false;
      }
      const candidateSourceScfId = String(candidate.parameters?.source_scf_id ?? "").trim();
      const candidateKPath = normalizeSavedKPath(candidate.parameters?.k_path);
      return candidateSourceScfId === sourceScfId && candidateKPath === kPath;
    });

    const settled = await Promise.allSettled(
      matchingBandRuns.map(async (candidate) => {
        const candidateDetail = await ensureCalculationDetails(candidate);
        const bandData = candidateDetail.result?.band_data ?? candidate.result?.band_data ?? null;
        if (!bandData) {
          return null;
        }
        const startedAt = candidateDetail.started_at ?? candidate.started_at;
        return {
          id: candidate.id,
          label: `Bands · ${new Date(startedAt).toLocaleString()}`,
          data: bandData,
          fermiEnergy: candidateDetail.result?.fermi_energy ?? candidate.result?.fermi_energy ?? null,
        } satisfies WannierBandOverlayOption;
      }),
    );

    return settled
      .flatMap((entry) => (entry.status === "fulfilled" && entry.value ? [entry.value] : []))
      .sort((a, b) => {
        const left = matchingBandRuns.find((calcRun) => calcRun.id === a.id)?.started_at ?? "";
        const right = matchingBandRuns.find((calcRun) => calcRun.id === b.id)?.started_at ?? "";
        return right.localeCompare(left);
      });
  }

  async function handleViewWannier(calc: CalculationRun) {
    try {
      const detail = await ensureCalculationDetails(calc);
      const wannierData = detail.result?.wannier_data ?? calc.result?.wannier_data ?? null;
      if (!wannierData) {
        setError("Saved Wannier data is unavailable for this calculation.");
        return;
      }
      const overlayOptions = await loadWannierBandOverlayOptions(calc, detail);
      onViewWannier(
        wannierData,
        detail.result?.fermi_energy ?? calc.result?.fermi_energy ?? null,
        overlayOptions,
      );
    } catch (e) {
      console.error("Failed to load Wannier data:", e);
      setError(`Failed to load Wannier data: ${e}`);
    }
  }

  async function handleViewTransport(calc: CalculationRun) {
    try {
      const detail = await ensureCalculationDetails(calc);
      const transportData = detail.result?.transport_data ?? calc.result?.transport_data ?? null;
      if (!transportData) {
        setError("Saved transport data is unavailable for this calculation.");
        return;
      }
      onViewTransport(transportData);
    } catch (e) {
      console.error("Failed to load transport data:", e);
      setError(`Failed to load transport data: ${e}`);
    }
  }

  async function handleViewFermiSurface(calc: CalculationRun, surfaceFile: string | null) {
    if (readOnly) return;
    setLaunchingFermiCalcId(calc.id);
    setError(null);
    try {
      await invoke("launch_fermi_surface_viewer", {
        projectId,
        calculationId: calc.id,
        surfaceFile,
      });
      setInfoMessage(
        `Opened ${surfaceFile ?? "primary Fermi-surface file"} in FermiSurfer.`,
      );
    } catch (e) {
      console.error("Failed to launch FermiSurfer:", e);
      setError(`Failed to launch FermiSurfer: ${e}`);
    } finally {
      setLaunchingFermiCalcId((current) => (current === calc.id ? null : current));
    }
  }

  async function handleDownloadHpcArtifacts(calc: CalculationRun) {
    if (!isHpcCalculation(calc)) return;
    setDownloadingCalcId(calc.id);
    setError(null);
    setInfoMessage(null);
    setDownloadProgressByCalcId((prev) => ({
      ...prev,
      [calc.id]: defaultHpcDownloadProgress("connecting"),
    }));
    let unlisten: UnlistenFn | null = null;
    let downloadCompleted = false;
    try {
      unlisten = await listen<string>(`task-output:${calc.id}`, (event) => {
        const line = String(event.payload || "");
        if (line.startsWith("HPC_STAGE|Connecting|")) {
          setDownloadProgressByCalcId((prev) => ({
            ...prev,
            [calc.id]: {
              ...(prev[calc.id] ?? defaultHpcDownloadProgress("connecting")),
              phase: "connecting",
            },
          }));
          return;
        }
        if (line.startsWith("HPC_STAGE|Collecting|")) {
          setDownloadProgressByCalcId((prev) => ({
            ...prev,
            [calc.id]: {
              ...(prev[calc.id] ?? defaultHpcDownloadProgress("collecting")),
              phase: "collecting",
            },
          }));
          return;
        }
        if (line.startsWith("HPC_TRANSFER|start|")) {
          const parts = line.split("|");
          setDownloadProgressByCalcId((prev) => ({
            ...prev,
            [calc.id]: {
              ...(prev[calc.id] ?? defaultHpcDownloadProgress("downloading")),
              phase: "downloading",
              totalFiles: parsePositiveInt(parts[3]),
              totalBytes: parsePositiveInt(parts[4]),
            },
          }));
          return;
        }
        if (line.startsWith("HPC_TRANSFER|progress|")) {
          const parts = line.split("|");
          setDownloadProgressByCalcId((prev) => ({
            ...prev,
            [calc.id]: {
              ...(prev[calc.id] ?? defaultHpcDownloadProgress("downloading")),
              phase: "downloading",
              processedFiles: parsePositiveInt(parts[3]),
              totalFiles: parsePositiveInt(parts[4]),
              downloadedBytes: parsePositiveInt(parts[5]),
              totalBytes: parsePositiveInt(parts[6]),
              skippedFiles: parsePositiveInt(parts[7]),
            },
          }));
          return;
        }
        if (line.startsWith("HPC_TRANSFER|done|")) {
          const parts = line.split("|");
          const downloadedFiles = parsePositiveInt(parts[3]);
          const totalFiles = parsePositiveInt(parts[4]);
          const skippedFiles = parsePositiveInt(parts[7]);
          setDownloadProgressByCalcId((prev) => ({
            ...prev,
            [calc.id]: {
              ...(prev[calc.id] ?? defaultHpcDownloadProgress("downloading")),
              phase: "downloading",
              processedFiles:
                downloadedFiles + skippedFiles > 0
                  ? downloadedFiles + skippedFiles
                  : totalFiles,
              totalFiles,
              downloadedBytes: parsePositiveInt(parts[5]),
              totalBytes: parsePositiveInt(parts[6]),
              skippedFiles,
            },
          }));
          return;
        }
        if (line.startsWith("HPC_STAGE|Saved|")) {
          setDownloadProgressByCalcId((prev) => ({
            ...prev,
            [calc.id]: {
              ...(prev[calc.id] ?? defaultHpcDownloadProgress("saved")),
              phase: "saved",
            },
          }));
        }
      });

      const report = await downloadHpcCalculationArtifacts(projectId, calc.id, null, true);
      const remoteStorageBytes = report.downloaded_bytes + report.skipped_bytes;
      setDownloadProgressByCalcId((prev) => {
        const current = prev[calc.id] ?? defaultHpcDownloadProgress("saved");
        return {
          ...prev,
          [calc.id]: {
            ...current,
            phase: "saved",
            processedFiles:
              current.totalFiles > 0 ? current.totalFiles : report.downloaded_files + report.skipped_files,
            totalFiles:
              current.totalFiles > 0 ? current.totalFiles : report.downloaded_files + report.skipped_files,
            downloadedBytes: report.downloaded_bytes,
            totalBytes:
              current.totalBytes > 0 ? current.totalBytes : report.downloaded_bytes + report.skipped_bytes,
            skippedFiles: report.skipped_files,
          },
        };
      });
      const skippedSuffix = report.skipped_files > 0 ? ` (${report.skipped_files} skipped)` : "";
      setInfoMessage(
        `Downloaded ${report.downloaded_files} files (${formatBytes(report.downloaded_bytes)})${skippedSuffix}.`,
      );
      updateCalculationRecords(calc.id, (entry) => {
        const nextParams = {
          ...(entry.parameters || {}),
          artifacts_downloaded_full: true,
          artifact_sync_mode: "full",
          artifact_synced_at: new Date().toISOString(),
          remote_storage_bytes: remoteStorageBytes,
        };
        return {
          ...entry,
          parameters: nextParams,
        };
      });
      downloadCompleted = true;
      await loadProject({ showLoading: false, refreshSelectedCif: false });
    } catch (e) {
      console.error("Failed to download HPC artifacts:", e);
      setDownloadProgressByCalcId((prev) => ({
        ...prev,
        [calc.id]: {
          ...(prev[calc.id] ?? defaultHpcDownloadProgress("error")),
          phase: "error",
        },
      }));
      setError(`Failed to download HPC artifacts: ${e}`);
    } finally {
      if (unlisten) {
        unlisten();
      }
      if (downloadCompleted) {
        setDownloadProgressByCalcId((prev) => {
          if (!(calc.id in prev)) return prev;
          const next = { ...prev };
          delete next[calc.id];
          return next;
        });
      }
      setDownloadingCalcId((current) => (current === calc.id ? null : current));
    }
  }

  function renderHpcDownloadProgress(calc: CalculationRun) {
    const progress = downloadProgressByCalcId[calc.id];
    if (!progress) return null;

    const percent = getHpcDownloadPercent(progress);
    const statusLabel =
      progress.phase === "connecting"
        ? "Connecting..."
        : progress.phase === "collecting"
          ? "Collecting..."
          : progress.phase === "downloading"
            ? "Downloading..."
            : progress.phase === "saved"
              ? "Saved"
              : "Failed";

    return (
      <div className="calc-download-progress">
        <div className="calc-download-progress-row">
          <span>{statusLabel}</span>
          <span>{percent == null ? "..." : `${percent}%`}</span>
        </div>
        <div className="calc-download-progress-bar">
          <div
            className={`calc-download-progress-fill ${percent == null ? "indeterminate" : ""} ${progress.phase === "error" ? "error" : ""}`}
            style={percent == null ? undefined : { width: `${percent}%` }}
          />
        </div>
        <div className="calc-download-progress-meta">
          <span>{`${progress.processedFiles}/${progress.totalFiles || "?"} files`}</span>
          <span>{`${formatBytes(progress.downloadedBytes)} / ${formatBytes(progress.totalBytes)}`}</span>
        </div>
      </div>
    );
  }

  function renderHpcDownloadButton(calc: CalculationRun) {
    if (readOnly) return null;
    if (!isHpcCalculation(calc)) return null;
    const isDownloading = downloadingCalcId === calc.id;
    const isDownloaded = isHpcArtifactsDownloaded(calc) && !isDownloading;
    return (
      <button
        className={`download-calc-btn${isDownloaded ? " downloaded" : ""}`}
        onClick={(e) => {
          e.stopPropagation();
          void handleDownloadHpcArtifacts(calc);
        }}
        disabled={isDownloading || isDownloaded}
      >
        {isDownloading ? "Downloading..." : isDownloaded ? "Downloaded" : "Download Full"}
      </button>
    );
  }

  function handleRunPhonons() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunPhonons(selectedCifId, crystalData, variant.calculations);
  }

  async function handleViewPhonon(
    calc: CalculationRun,
    viewMode: "bands" | "dos",
    basePhononData: { dos_data: any; dispersion_data: any },
  ) {
    let phononData = basePhononData;
    try {
      const detail = await ensureCalculationDetails(calc);
      const detailResult = detail.result as any;
      const detailPhononData = detailResult?.phonon_data
        ?? ((detailResult?.dos_data != null || detailResult?.dispersion_data != null)
          ? {
            dos_data: detailResult?.dos_data ?? null,
            dispersion_data: detailResult?.dispersion_data ?? null,
          }
          : null);
      if (detailPhononData) {
        phononData = {
          dos_data: detailPhononData.dos_data ?? phononData.dos_data ?? null,
          dispersion_data: detailPhononData.dispersion_data ?? phononData.dispersion_data ?? null,
        };
      }
    } catch (e) {
      console.warn("Failed to load phonon detail:", e);
    }
    const needsDispersion = viewMode === "bands" && phononData.dispersion_data == null;
    const needsDos = viewMode === "dos" && phononData.dos_data == null;

    if (needsDispersion || needsDos) {
      try {
        const recovered = await invoke<{ dos_data: any | null; dispersion_data: any | null }>(
          "get_saved_phonon_data",
          {
            projectId,
            calcId: calc.id,
          },
        );
        phononData = {
          dos_data: phononData.dos_data ?? recovered?.dos_data ?? null,
          dispersion_data: phononData.dispersion_data ?? recovered?.dispersion_data ?? null,
        };
      } catch (e) {
        console.warn("Failed to recover saved phonon data:", e);
      }
    }

    onViewPhonons(phononData, viewMode);
  }

  function hasConvergedSCF(): boolean {
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return false;
    return variant.calculations.some(c => c.calc_type === "scf" && c.result?.converged);
  }

  function hasWannierReadyScf(): boolean {
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return false;
    return variant.calculations.some((calc) => isWannierReadyScf(calc));
  }

  function hasSavedWannierCalculation(): boolean {
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return false;
    return variant.calculations.some(
      (calc) => calc.calc_type === "wannier" && (Boolean(calc.completed_at) || Boolean(calc.parameters?.seedname) || Boolean(calc.result?.wannier_data)),
    );
  }

  function getCalculationRuntime(calc: CalculationRun): CalculationRuntimeDisplay | null {
    const wallSeconds = calc.result?.wall_time_seconds;
    if (typeof wallSeconds === "number" && Number.isFinite(wallSeconds) && wallSeconds > 0) {
      return { kind: "wall", seconds: wallSeconds };
    }

    const rawOutput = calc.result?.raw_output;
    if (typeof rawOutput === "string" && rawOutput.trim().length > 0) {
      const parsedRuntime = parseRuntimeFromRawOutput(rawOutput);
      if (parsedRuntime) return parsedRuntime;
    }

    if (!calc.completed_at) return null;
    const startMs = Date.parse(calc.started_at);
    const endMs = Date.parse(calc.completed_at);
    if (Number.isFinite(startMs) && Number.isFinite(endMs) && endMs > startMs) {
      return { kind: "wall", seconds: (endMs - startMs) / 1000 };
    }

    return null;
  }

  function formatRuntimeDuration(seconds: number): string {
    const clamped = Math.max(0, seconds);
    if (clamped < 60) {
      return `${clamped.toFixed(1)} s`;
    }

    const rounded = Math.round(clamped);
    const hours = Math.floor(rounded / 3600);
    const minutes = Math.floor((rounded % 3600) / 60);
    const secs = rounded % 60;

    if (hours > 0) {
      return `${hours}h ${String(minutes).padStart(2, "0")}m ${String(secs).padStart(2, "0")}s`;
    }
    return `${minutes}m ${String(secs).padStart(2, "0")}s`;
  }

  function formatDate(isoString: string): string {
    try {
      const date = new Date(isoString);
      return date.toLocaleDateString(undefined, {
        year: "numeric",
        month: "short",
        day: "numeric",
        hour: "2-digit",
        minute: "2-digit",
      });
    } catch {
      return isoString;
    }
  }

function formatEnergy(energy: number): string {
  return `${energy.toFixed(6)} Ry`;
}

function normalizeSavedKPath(value: unknown): string {
  return String(value || "")
    .replace(/\s*→\s*/g, "→")
    .replace(/\s+/g, " ")
    .trim();
}

  function getSelectedVariant(): CifVariant | undefined {
    return project?.cif_variants.find(v => v.id === selectedCifId);
  }

  function handleCalculationSortModeChange(mode: CalculationSortMode) {
    setCalculationSortMode(mode);
    setStoredSortMode(mode);
  }

  async function togglePinnedCalculation(calcId: string, isPinned: boolean) {
    if (readOnly) return;
    if (!selectedCifId) return;
    const cifId = selectedCifId;
    const shouldBePinned = !isPinned;

    try {
      await invoke("set_calculation_tag", {
        projectId,
        cifId,
        calcId,
        tag: PINNED_TAG,
        enabled: shouldBePinned,
      });

      // Keep local UI state in sync without reloading the whole dashboard.
      updateCalculationRecords(calcId, (calc) => {
        const tags = Array.isArray(calc.tags) ? calc.tags.filter((tag) => tag !== PINNED_TAG) : [];
        if (shouldBePinned) {
          tags.push(PINNED_TAG);
        }
        return {
          ...calc,
          tags,
        };
      });
    } catch (e) {
      console.error("Failed to update pin status:", e);
      setError(`Failed to update pin status: ${e}`);
    }
  }

  const selectedVariant = getSelectedVariant();
  const pinnedCalcIds = useMemo<Set<string>>(() => {
    if (!selectedVariant) return new Set<string>();
    return new Set(
      selectedVariant.calculations
        .filter((calc) => Array.isArray(calc.tags) && calc.tags.includes(PINNED_TAG))
        .map((calc) => calc.id),
    );
  }, [selectedVariant]);
  const scfCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "scf") || [],
      calculationSortMode,
      "scf",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const bandCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "bands") || [],
      calculationSortMode,
      "bands",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const dosCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "dos") || [],
      calculationSortMode,
      "dos",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const wannierCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "wannier") || [],
      calculationSortMode,
      "wannier",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const transportCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "transport") || [],
      calculationSortMode,
      "transport",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const fermiSurfaceCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "fermi_surface") || [],
      calculationSortMode,
      "fermi_surface",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const phononCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => calc.calc_type === "phonon") || [],
      calculationSortMode,
      "phonon",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const optimizationCalculations = useMemo<CalculationRun[]>(
    () => sortCalculations(
      selectedVariant?.calculations.filter((calc) => isOptimizationCalculation(calc)) || [],
      calculationSortMode,
      "optimization",
      pinnedCalcIds,
    ),
    [selectedVariant, calculationSortMode, pinnedCalcIds],
  );
  const primitiveCell = useMemo(() => {
    if (!crystalData) return null;
    return getPrimitiveCell(crystalData);
  }, [crystalData]);
  const conventionalCellMetrics = useMemo<DisplayCellMetrics | null>(() => {
    if (!crystalData) return null;
    return {
      a: crystalData.cell_length_a.value,
      b: crystalData.cell_length_b.value,
      c: crystalData.cell_length_c.value,
      alpha: crystalData.cell_angle_alpha.value,
      beta: crystalData.cell_angle_beta.value,
      gamma: crystalData.cell_angle_gamma.value,
    };
  }, [crystalData]);

  const primitiveCellMetrics = useMemo<DisplayCellMetrics | null>(() => {
    if (!primitiveCell) return null;

    const BOHR_TO_ANGSTROM = 0.529177;
    const a = primitiveCell.celldm1 * BOHR_TO_ANGSTROM;

    if (primitiveCell.ibrav === "cubic_f") {
      const v1: [number, number, number] = [0, a / 2, a / 2];
      const v2: [number, number, number] = [a / 2, 0, a / 2];
      const v3: [number, number, number] = [a / 2, a / 2, 0];
      return calculateMetricsFromVectors(v1, v2, v3);
    }

    if (primitiveCell.ibrav === "cubic_i") {
      const v1: [number, number, number] = [a / 2, a / 2, -a / 2];
      const v2: [number, number, number] = [-a / 2, a / 2, a / 2];
      const v3: [number, number, number] = [a / 2, -a / 2, a / 2];
      return calculateMetricsFromVectors(v1, v2, v3);
    }

    return {
      a,
      b: a,
      c: a,
      alpha: 90,
      beta: 90,
      gamma: 90,
    };
  }, [primitiveCell]);

  const hasPrimitiveDisplay = primitiveCellMetrics !== null;
  const displayedCellMetrics = cellViewMode === "primitive" && primitiveCellMetrics
    ? primitiveCellMetrics
    : conventionalCellMetrics;
  const displayedCellVolume = useMemo<number | null>(() => {
    if (!displayedCellMetrics) {
      return crystalData?.cell_volume ?? null;
    }
    return calculateVolumeFromMetrics(displayedCellMetrics);
  }, [displayedCellMetrics, crystalData]);

  useEffect(() => {
    if (!hasPrimitiveDisplay && cellViewMode === "primitive") {
      setCellViewMode("conventional");
    }
  }, [hasPrimitiveDisplay, cellViewMode]);

  if (isLoading) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Loading...</h2>
        </div>
      </div>
    );
  }

  if (error && !project) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Error</h2>
        </div>
        {infoMessage && <div className="info-banner">{infoMessage}</div>}
        <div className="error-banner">{error}</div>
      </div>
    );
  }

  if (!project) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Project not found</h2>
        </div>
      </div>
    );
  }

  // Empty state - no CIF files yet
  if (project.cif_variants.length === 0) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <div className="dashboard-title">
            <div className="dashboard-title-row">
              <h2>{project.name}</h2>
              {!readOnly && (
                <button
                  className="project-title-edit-btn"
                  type="button"
                  onClick={openEditProjectDialog}
                  title="Edit project"
                  aria-label="Edit project"
                >
                  <svg viewBox="0 0 24 24" fill="none" aria-hidden="true">
                    <path
                      d="M4 20h4l10-10a2.12 2.12 0 0 0-3-3L5 17v3z"
                      stroke="currentColor"
                      strokeWidth="2"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                    <path
                      d="M13.5 6.5l4 4"
                      stroke="currentColor"
                      strokeWidth="2"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                  </svg>
                </button>
              )}
            </div>
            {project.description && (
              <p className="dashboard-description">{project.description}</p>
            )}
          </div>
          <div className="dashboard-header-actions">
            <button
              className="dashboard-refresh-btn"
              onClick={() => void handleRefreshProject()}
              disabled={isRefreshingProject}
              title="Reload project data"
            >
              {isRefreshingProject ? "Refreshing..." : "Refresh"}
            </button>
          </div>
        </div>

        {infoMessage && <div className="info-banner">{infoMessage}</div>}
        {error && <div className="error-banner">{error}</div>}

        <div className="dashboard-content">
          <div className="empty-state">
            <div className="empty-icon">📄</div>
            <h3>No Structure File</h3>
            {readOnly ? (
              <p>No structures are available in this synced project snapshot.</p>
            ) : (
              <>
                <p>Import a CIF file to get started with your calculations</p>
                <button
                  className="add-structure-btn primary"
                  onClick={handleImportCIF}
                  disabled={isImporting}
                >
                  {isImporting ? "Importing..." : "Import CIF File"}
                </button>
              </>
            )}
          </div>
        </div>

        {/* Delete Dialog */}
        {!readOnly && (
          <>
            <EditProjectDialog
              isOpen={showEditProjectDialog}
              projectId={project.id}
              initialName={project.name}
              initialDescription={project.description}
              onClose={() => setShowEditProjectDialog(false)}
              onSaved={handleProjectMetadataSaved}
            />
            {showDeleteDialog && renderDeleteDialog()}
          </>
        )}
      </div>
    );
  }

  return (
    <div className="dashboard-container">
      <div className="dashboard-header">
        <button className="back-btn" onClick={onBack}>
          ← Back
        </button>
        <div className="dashboard-title">
          <div className="dashboard-title-row">
            <h2>{project.name}</h2>
            {!readOnly && (
              <button
                className="project-title-edit-btn"
                type="button"
                onClick={openEditProjectDialog}
                title="Edit project"
                aria-label="Edit project"
              >
                <svg viewBox="0 0 24 24" fill="none" aria-hidden="true">
                  <path
                    d="M4 20h4l10-10a2.12 2.12 0 0 0-3-3L5 17v3z"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                  <path
                    d="M13.5 6.5l4 4"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </svg>
              </button>
            )}
          </div>
          {project.description && (
            <p className="dashboard-description">{project.description}</p>
          )}
        </div>
        <div className="dashboard-header-actions">
          <button
            className="dashboard-refresh-btn"
            onClick={() => void handleRefreshProject()}
            disabled={isRefreshingProject}
            title="Reload project data"
          >
            {isRefreshingProject ? "Refreshing..." : "Refresh"}
          </button>
          <div className="structure-selector">
            <label className="structure-selector-label">Structure</label>
            <select
              value={selectedCifId || ""}
              onChange={(e) => selectCif(e.target.value)}
            >
              {project.cif_variants.map((variant) => (
                <option key={variant.id} value={variant.id}>
                  {variant.formula} ({variant.filename})
                </option>
              ))}
            </select>
            {!readOnly && (
              <button
                className="add-structure-inline-btn"
                onClick={handleImportCIF}
                disabled={isImporting}
                title="Add new structure"
              >
                +
              </button>
            )}
          </div>
        </div>
      </div>

      {infoMessage && <div className="info-banner">{infoMessage}</div>}
      {error && <div className="error-banner">{error}</div>}

      <div className="dashboard-content">
        {/* Structure Info Hero */}
        {crystalData && (
          <section className="structure-hero">
            <div className="hero-formula">
              {crystalData.formula_sum || crystalData.formula_structural || "Unknown"}
            </div>
            {hasPrimitiveDisplay && (
              <div className="hero-cell-toggle">
                <button
                  type="button"
                  className={`hero-cell-toggle-btn ${cellViewMode === "conventional" ? "active" : ""}`}
                  onClick={() => setCellViewMode("conventional")}
                  title="Show conventional CIF lattice parameters"
                >
                  Conventional
                </button>
                <button
                  type="button"
                  className={`hero-cell-toggle-btn ${cellViewMode === "primitive" ? "active" : ""}`}
                  onClick={() => setCellViewMode("primitive")}
                  title="Show primitive lattice parameters used by QE (when detected)"
                >
                  Primitive
                </button>
              </div>
            )}
            <div className="hero-details">
              <div className="hero-detail-item">
                <label>Space Group</label>
                <span>
                  {crystalData.space_group_HM || "N/A"}
                  {crystalData.space_group_IT_number && ` (#${crystalData.space_group_IT_number})`}
                </span>
              </div>
              <div className="hero-detail-item">
                <label>Lattice Parameters</label>
                <span>
                  a = {(displayedCellMetrics?.a ?? crystalData.cell_length_a.value).toFixed(4)} A,{" "}
                  b = {(displayedCellMetrics?.b ?? crystalData.cell_length_b.value).toFixed(4)} A,{" "}
                  c = {(displayedCellMetrics?.c ?? crystalData.cell_length_c.value).toFixed(4)} A
                </span>
              </div>
              <div className="hero-detail-item">
                <label>Angles</label>
                <span>
                  alpha = {(displayedCellMetrics?.alpha ?? crystalData.cell_angle_alpha.value).toFixed(2)} deg,{" "}
                  beta = {(displayedCellMetrics?.beta ?? crystalData.cell_angle_beta.value).toFixed(2)} deg,{" "}
                  gamma = {(displayedCellMetrics?.gamma ?? crystalData.cell_angle_gamma.value).toFixed(2)} deg
                </span>
              </div>
              {hasPrimitiveDisplay && (
                <div className="hero-detail-item">
                  <label>Cell View</label>
                  <span>{cellViewMode === "primitive" ? "Primitive (QE)" : "Conventional (CIF)"}</span>
                </div>
              )}
              {displayedCellVolume !== null && (
                <div className="hero-detail-item">
                  <label>Volume</label>
                  <span>{displayedCellVolume.toFixed(2)} A^3</span>
                </div>
              )}
              <div className="hero-detail-item">
                <label>Atoms</label>
                <span>{crystalData.atom_sites.length} sites</span>
              </div>
            </div>
          </section>
        )}

        {/* Calculation Actions */}
        <section className="actions-section">
          <div className="actions-section-header">
            <h3>Calculations</h3>
            <div className="history-sort-control">
              <label htmlFor="dashboard-sort-mode">Sort Entries</label>
              <select
                id="dashboard-sort-mode"
                value={calculationSortMode}
                onChange={(e) => handleCalculationSortModeChange(e.target.value as CalculationSortMode)}
              >
                <option value="recent">Most Recent</option>
                <option value="best">Best</option>
              </select>
            </div>
          </div>
          {!readOnly && (
            <div className="calc-action-grid">
              <button className="calc-action-btn" onClick={handleRunSCF}>
                <span className="calc-action-icon">SCF</span>
                <span className="calc-action-label">Self-Consistent Field</span>
                <span className="calc-action-hint">Ground state energy</span>
              </button>
              <button
                className="calc-action-btn"
                onClick={handleRunBands}
                disabled={!hasConvergedSCF()}
              >
                <span className="calc-action-icon">Band</span>
                <span className="calc-action-label">Band Structure</span>
                <span className="calc-action-hint">
                  {hasConvergedSCF() ? "Electronic bands" : "Requires SCF"}
                </span>
              </button>
              <button
                className="calc-action-btn"
                onClick={handleRunDos}
                disabled={!hasConvergedSCF()}
              >
                <span className="calc-action-icon">DOS</span>
                <span className="calc-action-label">Electronic DOS</span>
                <span className="calc-action-hint">
                  {hasConvergedSCF() ? "Total density of states" : "Requires SCF"}
                </span>
              </button>
              <button
                className="calc-action-btn"
                onClick={handleRunWannier}
                disabled={!hasWannierReadyScf()}
              >
                <span className="calc-action-icon">W90</span>
                <span className="calc-action-label">Wannier90</span>
                <span className="calc-action-hint">
                  {hasWannierReadyScf() ? "MLWFs + interpolated bands" : "Requires primitive scalar SCF"}
                </span>
              </button>
              <button
                className="calc-action-btn"
                onClick={handleRunTransport}
                disabled={!hasSavedWannierCalculation()}
              >
                <span className="calc-action-icon">BW</span>
                <span className="calc-action-label">BoltzWann Transport</span>
                <span className="calc-action-hint">
                  {hasSavedWannierCalculation() ? "Transport via postw90.x" : "Requires saved Wannier"}
                </span>
              </button>
              <button
                className="calc-action-btn"
                onClick={handleRunFermiSurface}
                disabled={!hasConvergedSCF()}
              >
                <span className="calc-action-icon">FS</span>
                <span className="calc-action-label">Fermi Surface</span>
                <span className="calc-action-hint">
                  {hasConvergedSCF() ? "Generate FRMSF" : "Requires SCF"}
                </span>
              </button>
              <button
                className="calc-action-btn"
                onClick={handleRunPhonons}
                disabled={!hasConvergedSCF()}
              >
                <span className="calc-action-icon">Ph</span>
                <span className="calc-action-label">Phonons</span>
                <span className="calc-action-hint">
                  {hasConvergedSCF() ? "DOS & Dispersion" : "Requires SCF"}
                </span>
              </button>
              <button className="calc-action-btn" onClick={handleRunOptimization}>
                <span className="calc-action-icon">Opt</span>
                <span className="calc-action-label">Geometry Optimization</span>
                <span className="calc-action-hint">VC-Relax preset</span>
              </button>
            </div>
          )}
        </section>

        {/* SCF Calculations */}
        {scfCalculations.length > 0 && (
          <section className="history-section">
            <h3>SCFs</h3>
            <div className="calculations-list">
              {scfCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                return (
                  <div key={calc.id} className="calculation-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">SCF</span>
                        {calc.result && (
                          <span
                            className={`calc-status ${
                              calc.result.converged ? "converged" : "failed"
                            }`}
                          >
                            {calc.result.converged ? "Converged" : "Not converged"}
                          </span>
                        )}
                        {calc.result?.total_energy && (
                          <span className="calc-energy">
                            E = {formatEnergy(calc.result.total_energy)}
                          </span>
                        )}
                        <div className="calc-tags">
                          {getCalculationTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && calc.result && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          {calc.result.total_energy && (
                            <div className="detail-item">
                              <label>Total Energy</label>
                              <span>{formatEnergy(calc.result.total_energy)}</span>
                            </div>
                          )}
                          {calc.result.fermi_energy && (
                            <div className="detail-item">
                              <label>Fermi Energy</label>
                              <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                            </div>
                          )}
                          {calc.result.n_scf_steps && (
                            <div className="detail-item">
                              <label>SCF Steps</label>
                              <span>{calc.result.n_scf_steps}</span>
                            </div>
                          )}
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        <div className="detail-item parameters">
                          <label>Parameters</label>
                          <pre>{JSON.stringify(calc.parameters, null, 2)}</pre>
                        </div>
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete Calculation
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {/* Band Structure Calculations */}
        {bandCalculations.length > 0 && (
          <section className="history-section bands-section">
            <h3>Bands</h3>
            <div className="calculations-list">
              {bandCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                return (
                  <div key={calc.id} className="calculation-item bands-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">BANDS</span>
                        {calc.parameters?.k_path && (
                          <span className="calc-kpath">{calc.parameters.k_path}</span>
                        )}
                        <div className="calc-tags">
                          {getBandsTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>K-Path</label>
                            <span>{calc.parameters?.k_path || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Total K-Points</label>
                            <span>{calc.parameters?.total_k_points || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Number of Bands</label>
                            <span>{calc.parameters?.n_bands || "N/A"}</span>
                          </div>
                          {calc.result?.fermi_energy && (
                            <div className="detail-item">
                              <label>Fermi Energy</label>
                              <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                            </div>
                          )}
                          <div className="detail-item">
                            <label>Source SCF</label>
                            <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                          </div>
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {calc.result && (
                            <button
                              className="view-bands-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                void handleViewBands(calc);
                              }}
                            >
                              View Bands
                            </button>
                          )}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {/* Electronic DOS Calculations */}
        {dosCalculations.length > 0 && (
          <section className="history-section dos-section">
            <h3>Electronic DOS</h3>
            <div className="calculations-list">
              {dosCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                const energyMin = Number(calc.parameters?.dos_emin);
                const energyMax = Number(calc.parameters?.dos_emax);
                return (
                  <div key={calc.id} className="calculation-item dos-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">DOS</span>
                        <div className="calc-tags">
                          {getDosTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>DOS K-Grid</label>
                            <span>
                              {calc.parameters?.dos_k_grid
                                ? `${calc.parameters.dos_k_grid[0]}×${calc.parameters.dos_k_grid[1]}×${calc.parameters.dos_k_grid[2]}`
                                : "N/A"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>DeltaE</label>
                            <span>{calc.parameters?.dos_delta_e ?? "N/A"} eV</span>
                          </div>
                          <div className="detail-item">
                            <label>Energy Window</label>
                            <span>
                              {Number.isFinite(energyMin) && Number.isFinite(energyMax)
                                ? `${energyMin.toFixed(2)} to ${energyMax.toFixed(2)} eV`
                                : "Auto"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>DOS Points</label>
                            <span>{calc.parameters?.n_points || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Source SCF</label>
                            <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                          </div>
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {calc.result?.fermi_energy != null && (
                            <div className="detail-item">
                              <label>Fermi Energy</label>
                              <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {calc.result && (
                            <button
                              className="view-dos-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                void handleViewDos(calc);
                              }}
                            >
                              View DOS
                            </button>
                          )}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {wannierCalculations.length > 0 && (
          <section className="history-section wannier-section">
            <h3>Wannier90</h3>
            <div className="calculations-list">
              {wannierCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                const wannierData = calcData.result?.wannier_data ?? null;
                const bandData = wannierData?.band_data ?? calcData.result?.band_data ?? null;
                const wannierIssues = getWannierQualityIssues(
                  wannierData,
                  calcData.result?.raw_output ?? calc.result?.raw_output ?? null,
                  calcData.result?.fermi_energy ?? calc.result?.fermi_energy ?? null,
                );
                const totalSpread = Number(wannierData?.total_spread ?? calc.parameters?.total_spread);
                return (
                  <div key={calc.id} className="calculation-item bands-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">WANNIER</span>
                        {calc.parameters?.k_path && (
                          <span className="calc-kpath">{calc.parameters.k_path}</span>
                        )}
                        <div className="calc-tags">
                          {getWannierTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        {wannierIssues.length > 0 && (
                          <div className="warning-banner">
                            {wannierIssues.map((issue) => issue.message).join(" ")}
                          </div>
                        )}
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>K-Mesh</label>
                            <span>
                              {calc.parameters?.k_grid
                                ? `${calc.parameters.k_grid[0]}×${calc.parameters.k_grid[1]}×${calc.parameters.k_grid[2]}`
                                : "N/A"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>num_wann</label>
                            <span>{wannierData?.num_wann ?? calc.parameters?.n_wann ?? "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>num_bands</label>
                            <span>{wannierData?.num_bands ?? calc.parameters?.n_bands ?? "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Interpolated K-Points</label>
                            <span>{bandData?.n_kpoints ?? calc.parameters?.total_k_points ?? "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Total Spread</label>
                            <span>{Number.isFinite(totalSpread) ? `${totalSpread.toFixed(6)} A^2` : "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Source SCF</label>
                            <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Minimization</label>
                            <span>{formatWannierConvergenceFlag(wannierData?.convergence?.minimization_converged)}</span>
                          </div>
                          <div className="detail-item">
                            <label>Disentanglement</label>
                            <span>{formatWannierConvergenceFlag(wannierData?.convergence?.disentanglement_converged)}</span>
                          </div>
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {calc.result?.fermi_energy != null && (
                            <div className="detail-item">
                              <label>Fermi Energy</label>
                              <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        {Array.isArray(calc.parameters?.projection_summary) && calc.parameters.projection_summary.length > 0 && (
                          <div className="detail-item parameters">
                            <label>Projections</label>
                            <pre>{calc.parameters.projection_summary.join("\n")}</pre>
                          </div>
                        )}
                        {(wannierData?.convergence?.failure_reasons?.length || wannierData?.convergence?.warnings?.length) ? (
                          <div className="detail-item parameters">
                            <label>Quality Checks</label>
                            <pre>
                              {[
                                ...(wannierData?.convergence?.failure_reasons || []).map((entry: string) => `Error: ${entry}`),
                                ...(wannierData?.convergence?.warnings || []).map((entry: string) => `Warning: ${entry}`),
                              ].join("\n")}
                            </pre>
                          </div>
                        ) : null}
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {calc.result && (
                            <button
                              className="view-bands-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                void handleViewWannier(calc);
                              }}
                            >
                              View Wannier
                            </button>
                          )}
                          {!readOnly && (
                            <button
                              className="download-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openLudwigExportDialog(calc);
                              }}
                              disabled={isExportingLudwig}
                            >
                              {isExportingLudwig && calcToExportLudwig?.id === calc.id
                                ? "Exporting..."
                                : "Export to Ludwig"}
                            </button>
                          )}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {transportCalculations.length > 0 && (
          <section className="history-section dos-section">
            <h3>BoltzWann Transport</h3>
            <div className="calculations-list">
              {transportCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                const transportData = calcData.result?.transport_data ?? null;
                return (
                  <div key={calc.id} className="calculation-item dos-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">BOLTZWANN</span>
                        <div className="calc-tags">
                          {getTransportTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>Source Wannier</label>
                            <span>{calc.parameters?.source_wannier_calc_id?.slice(0, 8) || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Boltz k-mesh</label>
                            <span>
                              {calc.parameters?.boltz_kmesh
                                ? `${calc.parameters.boltz_kmesh[0]}×${calc.parameters.boltz_kmesh[1]}×${calc.parameters.boltz_kmesh[2]}`
                                : "N/A"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>μ points</label>
                            <span>{calc.parameters?.mu_points ?? transportData?.mu_values_ev?.length ?? "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>T points</label>
                            <span>{calc.parameters?.temperature_points ?? transportData?.temperature_values_k?.length ?? "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>τ</label>
                            <span>
                              {Number.isFinite(Number(calc.parameters?.relaxation_time_fs))
                                ? `${Number(calc.parameters.relaxation_time_fs).toFixed(2)} fs`
                                : "N/A"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>Engine</label>
                            <span>{calc.parameters?.engine || transportData?.engine || "N/A"}</span>
                          </div>
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {calc.result?.fermi_energy != null && (
                            <div className="detail-item">
                              <label>Reference Fermi Energy</label>
                              <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {(transportData || calc.result) && (
                            <button
                              className="view-dos-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                void handleViewTransport(calc);
                              }}
                            >
                              View Transport
                            </button>
                          )}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {/* Fermi Surface Calculations */}
        {!readOnly && fermiSurfaceCalculations.length > 0 && (
          <section className="history-section fermi-surface-section">
            <h3>Fermi Surface</h3>
            <div className="calculations-list">
              {fermiSurfaceCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                const frmsfFiles = Array.isArray(calc.parameters?.frmsf_files)
                  ? calc.parameters.frmsf_files
                  : [];
                const legacyBxsfFiles = Array.isArray(calc.parameters?.bxsf_files)
                  ? calc.parameters.bxsf_files
                  : [];
                const surfaceFiles = frmsfFiles.length > 0 ? frmsfFiles : legacyBxsfFiles;
                const primaryFile = calc.parameters?.primary_frmsf_file
                  ?? calc.parameters?.primary_bxsf_file
                  ?? (surfaceFiles.length > 0 ? surfaceFiles[0] : null);
                const totalBytes = Number(
                  calc.parameters?.total_frmsf_bytes
                  ?? calc.parameters?.total_bxsf_bytes,
                );
                const fileKind = (frmsfFiles.length > 0 || calc.parameters?.n_frmsf_files) ? "FRMSF" : "BXSF";
                return (
                  <div key={calc.id} className="calculation-item fermi-surface-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">FERMI</span>
                        <div className="calc-tags">
                          {getFermiSurfaceTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>Dense K-Grid</label>
                            <span>
                              {calc.parameters?.fermi_k_grid
                                ? `${calc.parameters.fermi_k_grid[0]}×${calc.parameters.fermi_k_grid[1]}×${calc.parameters.fermi_k_grid[2]}`
                                : "N/A"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>File Format</label>
                            <span>{fileKind}</span>
                          </div>
                          <div className="detail-item">
                            <label>FermiSurfer Files</label>
                            <span>{calc.parameters?.n_frmsf_files ?? calc.parameters?.n_bxsf_files ?? surfaceFiles.length ?? "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Primary File</label>
                            <span>{primaryFile || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Total File Size</label>
                            <span>{Number.isFinite(totalBytes) ? formatBytes(totalBytes) : "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Source SCF</label>
                            <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                          </div>
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {calc.result?.fermi_energy != null && (
                            <div className="detail-item">
                              <label>Fermi Energy</label>
                              <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        {surfaceFiles.length > 0 && (
                          <div className="detail-item parameters">
                            <label>{fileKind} Files</label>
                            <pre>{JSON.stringify(surfaceFiles, null, 2)}</pre>
                          </div>
                        )}
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          <button
                            className="view-fermi-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              void handleViewFermiSurface(calc, primaryFile);
                            }}
                            disabled={launchingFermiCalcId === calc.id || (!primaryFile && surfaceFiles.length === 0)}
                          >
                            {launchingFermiCalcId === calc.id ? "Launching..." : "View Fermi Surface"}
                          </button>
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {/* Phonon Calculations */}
        {phononCalculations.length > 0 && (
          <section className="history-section phonon-section">
            <h3>Phonons</h3>
            <div className="calculations-list">
              {phononCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
                const calcData = getCalculationRecord(calc);
                const runtime = getCalculationRuntime(calcData);
                const resultData = calcData.result as any;
                const phononData = resultData?.phonon_data
                  ?? ((resultData?.dos_data != null || resultData?.dispersion_data != null)
                    ? {
                      dos_data: resultData?.dos_data ?? null,
                      dispersion_data: resultData?.dispersion_data ?? null,
                    }
                    : null);
                const hasDispersion = phononData?.dispersion_data != null || calc.parameters?.calculate_dispersion === true;
                const hasDos = phononData?.dos_data != null || calc.parameters?.calculate_dos === true;
                const fallbackPhononData = phononData ?? {
                  dos_data: null,
                  dispersion_data: null,
                };
                return (
                  <div key={calc.id} className="calculation-item phonon-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
                    >
                      <div className="calculation-info">
                        <span className="calc-type">PHONON</span>
                        {calc.parameters?.q_path && (
                          <span className="calc-kpath">{calc.parameters.q_path}</span>
                        )}
                        <div className="calc-tags">
                          {getPhononTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>Q-Grid</label>
                            <span>
                              {calc.parameters?.q_grid
                                ? `${calc.parameters.q_grid[0]}×${calc.parameters.q_grid[1]}×${calc.parameters.q_grid[2]}`
                                : "N/A"}
                            </span>
                          </div>
                          <div className="detail-item">
                            <label>Modes</label>
                            <span>{calc.parameters?.n_modes || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Q-Points</label>
                            <span>{calc.parameters?.n_qpoints || "N/A"}</span>
                          </div>
                          <div className="detail-item">
                            <label>DOS</label>
                            <span>{calc.parameters?.calculate_dos ? "Yes" : "No"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Dispersion</label>
                            <span>{calc.parameters?.calculate_dispersion ? "Yes" : "No"}</span>
                          </div>
                          <div className="detail-item">
                            <label>Source SCF</label>
                            <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                          </div>
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {hasDispersion && (
                            <button
                              className="view-phonon-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                void handleViewPhonon(calc, "bands", fallbackPhononData);
                              }}
                            >
                              View Bands
                            </button>
                          )}
                          {hasDos && (
                            <button
                              className="view-phonon-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                void handleViewPhonon(calc, "dos", fallbackPhononData);
                              }}
                            >
                              View DOS
                            </button>
                          )}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}

        {/* Geometry Optimization Calculations */}
        {optimizationCalculations.length > 0 && (
          <section className="history-section">
            <h3>Optimizations</h3>
            <div className="calculations-list">
              {optimizationCalculations.map((calc) => {
                const calcData = getCalculationRecord(calc);
                const storedOptimizedStructure = isSavedStructureData(calc.parameters?.optimized_structure)
                  ? calc.parameters.optimized_structure
                  : null;
                const storedSourceStructure = isSavedStructureData(calc.parameters?.source_structure)
                  ? calc.parameters.source_structure
                  : null;
                const parsedOptimizedStructure = extractOptimizedStructure(
                  calcData.result?.raw_output || "",
                  storedOptimizedStructure || storedSourceStructure,
                );
                const optimizedStructure = parsedOptimizedStructure || storedOptimizedStructure;
                const summary = (optimizedStructure ? summarizeCell(optimizedStructure) : null)
                  || (calc.parameters?.optimized_cell_summary as SavedCellSummary | null)
                  || null;
                const optimizedCell = asCellMatrix(optimizedStructure?.cell_parameters);
                const units = summary?.units || optimizedStructure?.cell_units || "angstrom";
                const displayCellContext = resolveOptimizationDisplayCellContext(calc, crystalData);
                const likelyPrimitiveCell = optimizedCell
                  ? isLikelyPrimitiveOptimizationCell(calc, optimizedCell, displayCellContext, crystalData)
                  : false;
                let displaySummary: SavedCellSummary | null = summary ? { ...summary } : null;
                let cellBasisLabel = likelyPrimitiveCell
                  ? (displayCellContext?.centering === "R" ? "Rhombohedral primitive (QE output)" : "Primitive (QE output)")
                  : "Stored (QE output)";

                if (!displaySummary && optimizedCell) {
                  const metrics = calculateMetricsFromVectors(optimizedCell[0], optimizedCell[1], optimizedCell[2]);
                  displaySummary = {
                    a: metrics.a,
                    b: metrics.b,
                    c: metrics.c,
                    alpha: metrics.alpha,
                    beta: metrics.beta,
                    gamma: metrics.gamma,
                    volume: calculateVolumeFromVectors(optimizedCell[0], optimizedCell[1], optimizedCell[2]),
                    units,
                  };
                }

                if (
                  cellViewMode === "conventional"
                  && likelyPrimitiveCell
                  && optimizedCell
                  && displayCellContext
                ) {
                  const convertedCell = convertPrimitiveToConventionalCell(
                    optimizedCell,
                    displayCellContext.centering,
                    displayCellContext.rhombohedralSetting,
                  );
                  if (convertedCell) {
                    const metrics = calculateMetricsFromVectors(convertedCell[0], convertedCell[1], convertedCell[2]);
                    displaySummary = {
                      a: metrics.a,
                      b: metrics.b,
                      c: metrics.c,
                      alpha: metrics.alpha,
                      beta: metrics.beta,
                      gamma: metrics.gamma,
                      volume: calculateVolumeFromVectors(convertedCell[0], convertedCell[1], convertedCell[2]),
                      units,
                    };
                    cellBasisLabel = "Conventional (derived)";
                  }
                }

                const mode = getOptimizationMode(calc);
                const modeLabel = mode === "vcrelax" ? "VC-Relax" : "Relax";
                const sourceLabel =
                  calc.parameters?.structure_source?.type === "optimization"
                    ? `Saved optimization ${String(calc.parameters?.structure_source?.calc_id || "").slice(0, 8)}`
                    : "From CIF";
                const electronConvLabel = formatThreshold(calc.parameters?.conv_thr);
                const forceConvLabel = formatThreshold(calc.parameters?.forc_conv_thr);
                const energyConvLabel = formatThreshold(calc.parameters?.etot_conv_thr);
                const pressValue = Number(calc.parameters?.press);
                const isPinned = pinnedCalcIds.has(calc.id);
                const runtime = getCalculationRuntime(calcData);

                return (
                  <div key={calc.id} className="calculation-item">
                    <div
                      className="calculation-header"
                      onClick={() => {
                        if (expandedCalc !== calc.id) {
                          void ensureCalculationDetails(calc).catch((e) => {
                            console.warn("Failed to load optimization detail:", e);
                          });
                        }
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id);
                      }}
                    >
                      <div className="calculation-info">
                        <span className="calc-type">OPT</span>
                        {calc.result && (
                          <span
                            className={`calc-status ${
                              calc.result.converged ? "converged" : "failed"
                            }`}
                          >
                            {calc.result.converged ? "Converged" : "Not converged"}
                          </span>
                        )}
                        {calc.result?.total_energy && (
                          <span className="calc-energy">
                            E = {formatEnergy(calc.result.total_energy)}
                          </span>
                        )}
                        <div className="calc-tags">
                          {getOptimizationTags(calc).map((tag, i) => (
                            <span key={i} className={getCalcTagClass(tag)}>
                              {tag.label}
                            </span>
                          ))}
                        </div>
                      </div>
                      <div className="calculation-meta">
                        <button
                          type="button"
                          className={`pin-calc-btn ${isPinned ? "pinned" : ""}`}
                          onClick={(e) => {
                            e.stopPropagation();
                            void togglePinnedCalculation(calc.id, isPinned);
                          }}
                          disabled={readOnly}
                          title={isPinned ? "Unpin calculation" : "Pin calculation"}
                          aria-label={isPinned ? "Unpin calculation" : "Pin calculation"}
                        >
                          <svg viewBox="0 0 24 24" aria-hidden="true">
                            <path d="M12 2.5L14.9 8.38L21.4 9.33L16.7 13.91L17.81 20.38L12 17.33L6.19 20.38L7.3 13.91L2.6 9.33L9.1 8.38L12 2.5Z" />
                          </svg>
                        </button>
                        <span className="calc-date">
                          {calc.completed_at
                            ? formatDate(calc.completed_at)
                            : "In progress..."}
                        </span>
                        {runtime && (
                          <span className="calc-runtime">
                            {formatRuntimeDuration(runtime.seconds)}
                          </span>
                        )}
                        {calc.storage_bytes != null && (
                          <span className="calc-size">{formatBytes(calc.storage_bytes)}</span>
                        )}
                        <span className="expand-icon">
                          {expandedCalc === calc.id ? "▼" : "▶"}
                        </span>
                      </div>
                    </div>

                    {expandedCalc === calc.id && (
                      <div className="calculation-details">
                        <div className="details-grid">
                          <div className="detail-item">
                            <label>Mode</label>
                            <span>{modeLabel}</span>
                          </div>
                          <div className="detail-item">
                            <label>Source Structure</label>
                            <span>{sourceLabel}</span>
                          </div>
                          <div className="detail-item">
                            <label>Atoms</label>
                            <span>{optimizedStructure?.atoms?.length || "N/A"}</span>
                          </div>
                          {electronConvLabel && (
                            <div className="detail-item">
                              <label>Electron Conv.</label>
                              <span>{electronConvLabel}</span>
                            </div>
                          )}
                          {forceConvLabel && (
                            <div className="detail-item">
                              <label>Force Conv.</label>
                              <span>{forceConvLabel}</span>
                            </div>
                          )}
                          {energyConvLabel && (
                            <div className="detail-item">
                              <label>Energy Conv.</label>
                              <span>{energyConvLabel}</span>
                            </div>
                          )}
                          {mode === "vcrelax" && Number.isFinite(pressValue) && (
                            <div className="detail-item">
                              <label>Target Pressure</label>
                              <span>{pressValue.toFixed(2)} kbar</span>
                            </div>
                          )}
                          {displaySummary && (
                            <div className="detail-item">
                              <label>Cell Basis</label>
                              <span>{cellBasisLabel}</span>
                            </div>
                          )}
                          {displaySummary && (
                            <div className="detail-item">
                              <label>Lattice ({displaySummary.units})</label>
                              <span>{`${displaySummary.a.toFixed(4)} / ${displaySummary.b.toFixed(4)} / ${displaySummary.c.toFixed(4)}`}</span>
                            </div>
                          )}
                          {displaySummary && (
                            <div className="detail-item">
                              <label>Angles (deg)</label>
                              <span>{`${displaySummary.alpha.toFixed(2)} / ${displaySummary.beta.toFixed(2)} / ${displaySummary.gamma.toFixed(2)}`}</span>
                            </div>
                          )}
                          {displaySummary && (
                            <div className="detail-item">
                              <label>Volume</label>
                              <span>{`${displaySummary.volume.toFixed(4)} ${displaySummary.units === "angstrom" ? "A^3" : `${displaySummary.units}^3`}`}</span>
                            </div>
                          )}
                          {runtime && (
                            <div className="detail-item">
                              <label>Time</label>
                              <span>{formatRuntimeDuration(runtime.seconds)}</span>
                            </div>
                          )}
                          {renderStorageDetailItems(calc)}
                        </div>
                        <div className="calc-actions">
                          {renderHpcDownloadProgress(calc)}
                          {renderHpcDownloadButton(calc)}
                          {!readOnly && (
                            <button
                              className="delete-calc-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                openDeleteCalcDialog(calc.id, calc.calc_type);
                              }}
                            >
                              Delete
                            </button>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </section>
        )}
      </div>

      {!readOnly && (
        <>
          <EditProjectDialog
            isOpen={showEditProjectDialog}
            projectId={project.id}
            initialName={project.name}
            initialDescription={project.description}
            onClose={() => setShowEditProjectDialog(false)}
            onSaved={handleProjectMetadataSaved}
          />

          {showLudwigExportDialog && calcToExportLudwig && (
            <div className="dialog-overlay" onClick={closeLudwigExportDialog}>
              <div className="dialog-content" onClick={(e) => e.stopPropagation()}>
                <div className="dialog-header">
                  <h2>Export to Ludwig</h2>
                  <button
                    className="dialog-close"
                    onClick={closeLudwigExportDialog}
                    disabled={isExportingLudwig}
                  >
                    &times;
                  </button>
                </div>
                <div className="dialog-body">
                  <div className="save-form">
                    <div className="form-group">
                      <label>Source Calculation</label>
                      <input
                        type="text"
                        value={`Wannier ${calcToExportLudwig.id.slice(0, 8)}`}
                        disabled
                      />
                    </div>
                    <div className="form-group">
                      <label>Export Mode</label>
                      <select
                        value={ludwigExportMode}
                        onChange={(e) => {
                          const nextMode = e.target.value as LudwigExportMode;
                          setLudwigExportMode(nextMode);
                          if (nextMode === "strict_2d") {
                            setLudwigSliceCoordinateInput("0.0");
                          }
                        }}
                        disabled={isExportingLudwig}
                      >
                        <option value="strict_2d">Strict 2D</option>
                        <option value="quasi_2d_fixed_slice">Quasi-2D Fixed Slice</option>
                      </select>
                    </div>
                    <div className="form-group">
                      <label>Primary In-Plane Axis</label>
                      <select
                        value={String(ludwigPrimaryAxis)}
                        onChange={(e) => setLudwigPrimaryAxis(Number.parseInt(e.target.value, 10))}
                        disabled={isExportingLudwig}
                      >
                        {RECIPROCAL_AXIS_OPTIONS.map((option) => (
                          <option key={option.value} value={option.value}>
                            {option.label}
                          </option>
                        ))}
                      </select>
                    </div>
                    <div className="form-group">
                      <label>Secondary In-Plane Axis</label>
                      <select
                        value={String(ludwigSecondaryAxis)}
                        onChange={(e) => setLudwigSecondaryAxis(Number.parseInt(e.target.value, 10))}
                        disabled={isExportingLudwig}
                      >
                        {RECIPROCAL_AXIS_OPTIONS.map((option) => (
                          <option key={option.value} value={option.value}>
                            {option.label}
                          </option>
                        ))}
                      </select>
                    </div>
                    <div className="form-group">
                      <label>Slice Axis</label>
                      <input
                        type="text"
                        value={
                          RECIPROCAL_AXIS_OPTIONS.find(
                            (option) => option.value === findSliceAxis(ludwigPrimaryAxis, ludwigSecondaryAxis),
                          )?.label ?? "Choose two distinct in-plane axes"
                        }
                        disabled
                      />
                    </div>
                    <div className="form-group">
                      <label>Slice Coordinate (fractional)</label>
                      <input
                        type="number"
                        step="0.01"
                        value={ludwigSliceCoordinateInput}
                        onChange={(e) => setLudwigSliceCoordinateInput(e.target.value)}
                        disabled={isExportingLudwig || ludwigExportMode === "strict_2d"}
                      />
                    </div>
                    <div className="form-group">
                      <label>nkx</label>
                      <input
                        type="number"
                        min={2}
                        step={1}
                        value={ludwigNkxInput}
                        onChange={(e) => setLudwigNkxInput(e.target.value)}
                        disabled={isExportingLudwig}
                      />
                    </div>
                    <div className="form-group">
                      <label>nky</label>
                      <input
                        type="number"
                        min={2}
                        step={1}
                        value={ludwigNkyInput}
                        onChange={(e) => setLudwigNkyInput(e.target.value)}
                        disabled={isExportingLudwig}
                      />
                    </div>
                    <div className="form-group">
                      <label>Chemical Potential</label>
                      <select
                        value={ludwigChemicalPotentialMode}
                        onChange={(e) => setLudwigChemicalPotentialMode(e.target.value as LudwigChemicalPotentialMode)}
                        disabled={isExportingLudwig}
                      >
                        <option value="saved">Use Saved Fermi Energy</option>
                        <option value="override">Override Manually</option>
                      </select>
                    </div>
                    <div className="form-group">
                      <label>Chemical Potential Override (eV)</label>
                      <input
                        type="number"
                        step="0.001"
                        value={ludwigChemicalPotentialInput}
                        onChange={(e) => setLudwigChemicalPotentialInput(e.target.value)}
                        disabled={isExportingLudwig || ludwigChemicalPotentialMode !== "override"}
                      />
                    </div>
                    <p className="project-archive-hint">
                      QCortado will export a Ludwig-ready band bundle from the saved Wannier artifacts, then Ludwig can build
                      its own mesh and collision operator externally.
                    </p>
                    <p className="warning-text">
                      Ludwig remains a 2D code. Quasi-2D export here means a fixed reciprocal-space slice, not a full
                      kz-integrated transport calculation.
                    </p>
                  </div>
                </div>
                <div className="dialog-footer">
                  <button
                    className="dialog-btn cancel"
                    onClick={closeLudwigExportDialog}
                    disabled={isExportingLudwig}
                  >
                    Cancel
                  </button>
                  <button
                    className="dialog-btn save"
                    onClick={() => {
                      void handleExportLudwigBundle();
                    }}
                    disabled={isExportingLudwig}
                  >
                    {isExportingLudwig ? "Exporting..." : "Export Bundle"}
                  </button>
                </div>
              </div>
            </div>
          )}

          {/* Delete Project Dialog */}
          {showDeleteDialog && renderDeleteDialog()}

          {/* Delete Calculation Dialog */}
          {showDeleteCalcDialog && (
            <div className="dialog-overlay" onClick={() => !isDeletingCalc && setShowDeleteCalcDialog(false)}>
              <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
                <div className="dialog-header">
                  <h2>Delete Calculation</h2>
                  <button
                    className="dialog-close"
                    onClick={() => setShowDeleteCalcDialog(false)}
                    disabled={isDeletingCalc}
                  >
                    &times;
                  </button>
                </div>

                <div className="dialog-body">
                  <p className="exit-warning">
                    Are you sure you want to delete this {calcToDelete?.calcType.toUpperCase()} calculation?
                  </p>
                  <p className="exit-hint">
                    This will permanently remove the calculation results and all associated input/output files.
                  </p>
                </div>

                <div className="dialog-footer">
                  <button
                    className="dialog-btn cancel"
                    onClick={() => setShowDeleteCalcDialog(false)}
                    disabled={isDeletingCalc}
                  >
                    Cancel
                  </button>
                  <button
                    className="dialog-btn delete width-lock"
                    onClick={handleConfirmDeleteCalc}
                    disabled={isDeletingCalc}
                  >
                    {isDeletingCalc ? "Deleting..." : "Delete Calculation"}
                  </button>
                </div>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );

  function renderDeleteDialog() {
    return (
      <div className="dialog-overlay" onClick={() => !isDeleting && setShowDeleteDialog(false)}>
        <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
          <div className="dialog-header">
            <h2>Delete Project</h2>
            <button
              className="dialog-close"
              onClick={() => setShowDeleteDialog(false)}
              disabled={isDeleting}
            >
              &times;
            </button>
          </div>

          <div className="dialog-body">
            <div className="delete-warning">
              <p>
                You are about to permanently delete <strong>{project?.name}</strong> and all of its data:
              </p>
              <ul>
                <li>{project?.cif_variants.length} structure{project?.cif_variants.length !== 1 ? "s" : ""}</li>
                <li>
                  {project?.cif_variants.reduce((sum, v) => sum + v.calculations.length, 0)} calculation
                  {project?.cif_variants.reduce((sum, v) => sum + v.calculations.length, 0) !== 1 ? "s" : ""}
                </li>
                <li>All input/output files</li>
              </ul>
              <p className="delete-warning-emphasis">
                This action cannot be undone.
              </p>
            </div>

            <div className="form-group">
              <label>
                Type <code>{DELETE_CONFIRM_TEXT}</code> to confirm:
              </label>
              <input
                type="text"
                value={deleteConfirmText}
                onChange={(e) => setDeleteConfirmText(e.target.value)}
                placeholder={DELETE_CONFIRM_TEXT}
                disabled={isDeleting}
                autoFocus
              />
            </div>
          </div>

          <div className="dialog-footer">
            <button
              className="dialog-btn cancel"
              onClick={() => setShowDeleteDialog(false)}
              disabled={isDeleting}
            >
              Cancel
            </button>
            <button
              className="dialog-btn delete width-lock"
              onClick={handleConfirmDelete}
              disabled={deleteConfirmText !== DELETE_CONFIRM_TEXT || isDeleting}
            >
              {isDeleting ? "Deleting..." : "Delete Project"}
            </button>
          </div>
        </div>
      </div>
    );
  }
}
