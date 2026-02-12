// Project Dashboard - Main view for working with a project's structures and calculations

import { useState, useEffect, useRef, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { parseCIF } from "../lib/cifParser";
import { CrystalData, SCFPreset, OptimizedStructureOption, SavedCellSummary, SavedStructureData } from "../lib/types";
import { getPrimitiveCell } from "../lib/primitiveCell";
import { getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";

interface QEResult {
  converged: boolean;
  total_energy: number | null;
  fermi_energy: number | null;
  n_scf_steps: number | null;
  wall_time_seconds: number | null;
  raw_output: string;
  band_data?: any;  // Band structure data for bands calculations
  phonon_data?: any;  // Phonon data (DOS + dispersion) for phonon calculations
}

export interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: QEResult | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
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
  onRunPhonons: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewPhonons: (phononData: any, viewMode: "bands" | "dos") => void;
}

type CalcTagType = "info" | "feature" | "special" | "geometry";
type CellViewMode = "conventional" | "primitive";
type CalculationSortMode = "recent" | "best";
type CalculationCategory = "scf" | "bands" | "phonon" | "optimization";

interface DisplayCellMetrics {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

type CellMatrix = [[number, number, number], [number, number, number], [number, number, number]];

const CONFIRM_TEXT = "delete my project for good";
const SOC_PRIORITY_BOOST = 250;
const PINNED_TAG = "pinned";

// Helper to generate calculation feature tags from parameters
function getCalculationTags(calc: CalculationRun): { label: string; type: CalcTagType }[] {
  const tags: { label: string; type: CalcTagType }[] = [];
  const params = calc.parameters || {};
  let hasGeometryTag = false;

  // Special tags from stored tags array (phonon-ready, structure-optimized)
  if (calc.tags) {
    for (const tag of calc.tags) {
      if (tag === "phonon-ready") {
        tags.push({ label: "Phonon-Ready", type: "special" });
      } else if (tag === "structure-optimized") {
        tags.push({ label: "Optimized", type: "special" });
      } else if (tag === "geometry") {
        tags.push({ label: "Geometry", type: "geometry" });
        hasGeometryTag = true;
      }
    }
  }

  if (!hasGeometryTag && params.structure_source?.type === "optimization") {
    tags.push({ label: "Geometry", type: "geometry" });
  }

  // K-points grid
  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    tags.push({ label: `${k1}×${k2}×${k3}`, type: "info" });
  }

  // Convergence threshold
  if (params.conv_thr) {
    const thr = params.conv_thr;
    // Format as scientific notation if small
    const label = thr < 0.001 ? thr.toExponential(0) : thr.toString();
    tags.push({ label: label, type: "info" });
  }

  // Feature tags
  if (params.lspinorb) {
    tags.push({ label: "SOC", type: "feature" });
  }

  if (params.nspin === 4) {
    tags.push({ label: "Non-collinear", type: "feature" });
  } else if (params.nspin === 2) {
    tags.push({ label: "Magnetic", type: "feature" });
  }

  if (params.lda_plus_u) {
    tags.push({ label: "DFT+U", type: "feature" });
  }

  if (params.vdw_corr && params.vdw_corr !== "none") {
    tags.push({ label: "vdW", type: "feature" });
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

  return tags;
}

// Helper to generate bands-specific tags
function getBandsTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};

  // K-points info
  if (params.total_k_points) {
    tags.push({ label: `${params.total_k_points} k-pts`, type: "info" });
  }

  // Inherited feature tags from SCF
  if (params.lspinorb) {
    tags.push({ label: "SOC", type: "feature" });
  }
  if (params.nspin === 2) {
    tags.push({ label: "Magnetic", type: "feature" });
  }
  if (params.lda_plus_u) {
    tags.push({ label: "DFT+U", type: "feature" });
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    tags.push({ label: "vdW", type: "feature" });
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

function convertPrimitiveToConventionalCell(
  primitiveCell: CellMatrix,
  primitiveLatticeKind: "cubic_f" | "cubic_i" | "cubic_p" | null,
): CellMatrix | null {
  if (!primitiveLatticeKind) return null;

  const [p1, p2, p3] = primitiveCell;

  if (primitiveLatticeKind === "cubic_f") {
    // Conventional vectors from FCC primitive basis:
    // a = -p1 + p2 + p3, b = p1 - p2 + p3, c = p1 + p2 - p3
    return [
      [-p1[0] + p2[0] + p3[0], -p1[1] + p2[1] + p3[1], -p1[2] + p2[2] + p3[2]],
      [p1[0] - p2[0] + p3[0], p1[1] - p2[1] + p3[1], p1[2] - p2[2] + p3[2]],
      [p1[0] + p2[0] - p3[0], p1[1] + p2[1] - p3[1], p1[2] + p2[2] - p3[2]],
    ];
  }

  if (primitiveLatticeKind === "cubic_i") {
    // Conventional vectors from BCC primitive basis:
    // a = p1 + p3, b = p1 + p2, c = p2 + p3
    return [
      [p1[0] + p3[0], p1[1] + p3[1], p1[2] + p3[2]],
      [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]],
      [p2[0] + p3[0], p2[1] + p3[1], p2[2] + p3[2]],
    ];
  }

  if (primitiveLatticeKind === "cubic_p") {
    return primitiveCell;
  }

  return null;
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
  onBack,
  onDeleted,
  onRunSCF,
  onRunBands,
  onViewBands,
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

  // Settings menu state
  const [showSettingsMenu, setShowSettingsMenu] = useState(false);
  const settingsRef = useRef<HTMLDivElement>(null);

  // Delete project confirmation dialog state
  const [showDeleteDialog, setShowDeleteDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeleting, setIsDeleting] = useState(false);

  // Delete calculation confirmation dialog state
  const [showDeleteCalcDialog, setShowDeleteCalcDialog] = useState(false);
  const [calcToDelete, setCalcToDelete] = useState<{ calcId: string; calcType: string } | null>(null);
  const [isDeletingCalc, setIsDeletingCalc] = useState(false);

  // Import state
  const [isImporting, setIsImporting] = useState(false);
  const [isRecoveringPhonon, setIsRecoveringPhonon] = useState(false);
  const [infoMessage, setInfoMessage] = useState<string | null>(null);

  // Expanded calculation
  const [expandedCalc, setExpandedCalc] = useState<string | null>(null);
  const [calculationSortMode, setCalculationSortMode] = useState<CalculationSortMode>(() => getStoredSortMode());

  useEffect(() => {
    loadProject();
  }, [projectId]);

  // Close settings menu when clicking outside
  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (settingsRef.current && !settingsRef.current.contains(event.target as Node)) {
        setShowSettingsMenu(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  async function loadProject() {
    setIsLoading(true);
    setError(null);
    try {
      const proj = await invoke<Project>("get_project", { projectId });
      setProject(proj);

      // Determine which CIF to show
      if (proj.cif_variants.length > 0) {
        const cifToOpen = proj.last_opened_cif_id &&
          proj.cif_variants.some(v => v.id === proj.last_opened_cif_id)
          ? proj.last_opened_cif_id
          : proj.cif_variants[0].id;

        await selectCif(cifToOpen);
      }
    } catch (e) {
      console.error("Failed to load project:", e);
      setError(String(e));
    } finally {
      setIsLoading(false);
    }
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

  function openDeleteDialog() {
    setShowSettingsMenu(false);
    setDeleteConfirmText("");
    setShowDeleteDialog(true);
  }

  async function handleConfirmDelete() {
    if (deleteConfirmText !== CONFIRM_TEXT) return;

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
    if (!calcToDelete || !selectedCifId) return;

    setIsDeletingCalc(true);
    try {
      await invoke("delete_calculation", {
        projectId,
        cifId: selectedCifId,
        calcId: calcToDelete.calcId,
      });
      // Reload project to reflect changes
      await loadProject();
      setShowDeleteCalcDialog(false);
      setCalcToDelete(null);
    } catch (e) {
      console.error("Failed to delete calculation:", e);
      setError(String(e));
    } finally {
      setIsDeletingCalc(false);
    }
  }

  function getOptimizedStructureOptions(calculations: CalculationRun[]): OptimizedStructureOption[] {
    return calculations
      .filter((calc) => isOptimizationCalculation(calc))
      .map((calc) => {
        const structure = calc.parameters?.optimized_structure;
        if (!structure || !Array.isArray(structure.atoms) || structure.atoms.length === 0) {
          return null;
        }

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
          cellSummary: calc.parameters?.optimized_cell_summary ?? null,
        } as OptimizedStructureOption;
      })
      .filter((opt): opt is OptimizedStructureOption => opt !== null);
  }

  function handleRunSCF() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunSCF(
      selectedCifId,
      crystalData,
      cifContent,
      variant.filename,
      undefined,
      undefined,
      getOptimizedStructureOptions(variant.calculations),
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

  function handleRunPhonons() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunPhonons(selectedCifId, crystalData, variant.calculations);
  }

  async function handleRecoverPhonon() {
    if (!selectedCifId) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;

    const defaultTmpDir = "/tmp/qcortado_phonon";

    const fallbackScf = variant.calculations
      .filter((calc) => calc.calc_type === "scf" && calc.result?.converged)
      .sort((a, b) => {
        const aTime = a.completed_at ?? a.started_at;
        const bTime = b.completed_at ?? b.started_at;
        return bTime.localeCompare(aTime);
      })[0];

    setIsRecoveringPhonon(true);
    setError(null);
    setInfoMessage(null);
    try {
      // First try the default pipeline scratch path automatically.
      try {
        await invoke("recover_phonon_calculation", {
          projectId,
          cifId: selectedCifId,
          workingDir: defaultTmpDir,
          sourceScfId: fallbackScf?.id ?? null,
        });
        await loadProject();
        setInfoMessage(`Recovered phonon calculation from ${defaultTmpDir}.`);
        return;
      } catch {
        // If default location is unavailable, let the user pick a folder.
      }

      const selected = await open({
        multiple: false,
        directory: true,
        defaultPath: defaultTmpDir,
        title: "Select phonon scratch directory",
      });

      if (!selected || Array.isArray(selected)) {
        setInfoMessage("Phonon recovery canceled.");
        return;
      }

      await invoke("recover_phonon_calculation", {
        projectId,
        cifId: selectedCifId,
        workingDir: selected,
        sourceScfId: fallbackScf?.id ?? null,
      });
      await loadProject();
      setInfoMessage(`Recovered phonon calculation from ${selected}.`);
    } catch (e) {
      console.error("Failed to recover phonon calculation:", e);
      setError(`Phonon recovery failed: ${e}`);
    } finally {
      setIsRecoveringPhonon(false);
    }
  }

  async function handleViewPhonon(
    calc: CalculationRun,
    viewMode: "bands" | "dos",
    basePhononData: { dos_data: any; dispersion_data: any },
  ) {
    let phononData = basePhononData;
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

  function getSelectedVariant(): CifVariant | undefined {
    return project?.cif_variants.find(v => v.id === selectedCifId);
  }

  function handleCalculationSortModeChange(mode: CalculationSortMode) {
    setCalculationSortMode(mode);
    setStoredSortMode(mode);
  }

  async function togglePinnedCalculation(calcId: string, isPinned: boolean) {
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
      setProject((prev) => {
        if (!prev) return prev;
        return {
          ...prev,
          cif_variants: prev.cif_variants.map((variant) => {
            if (variant.id !== cifId) return variant;
            return {
              ...variant,
              calculations: variant.calculations.map((calc) => {
                if (calc.id !== calcId) return calc;
                const tags = Array.isArray(calc.tags) ? calc.tags.filter((tag) => tag !== PINNED_TAG) : [];
                if (shouldBePinned) {
                  tags.push(PINNED_TAG);
                }
                return {
                  ...calc,
                  tags,
                };
              }),
            };
          }),
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
  const primitiveLatticeKind = primitiveCell?.ibrav ?? null;
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
            <h2>{project.name}</h2>
            {project.description && (
              <p className="dashboard-description">{project.description}</p>
            )}
          </div>
        </div>

        {infoMessage && <div className="info-banner">{infoMessage}</div>}
        {error && <div className="error-banner">{error}</div>}

        <div className="dashboard-content">
          <div className="empty-state">
            <div className="empty-icon">📄</div>
            <h3>No Structure File</h3>
            <p>Import a CIF file to get started with your calculations</p>
            <button
              className="add-structure-btn primary"
              onClick={handleImportCIF}
              disabled={isImporting}
            >
              {isImporting ? "Importing..." : "Import CIF File"}
            </button>
          </div>
        </div>

        {/* Floating Settings Button */}
        <div className="floating-settings" ref={settingsRef}>
          <button
            className="floating-settings-btn"
            onClick={() => setShowSettingsMenu(!showSettingsMenu)}
            title="Project settings"
          >
            <svg width="24" height="24" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.061 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.533 1.533 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.533 1.533 0 01.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd" />
            </svg>
          </button>

          {showSettingsMenu && (
            <div className="floating-settings-menu">
              <button className="settings-menu-item danger" onClick={openDeleteDialog}>
                Delete Project
              </button>
            </div>
          )}
        </div>

        {/* Delete Dialog */}
        {showDeleteDialog && renderDeleteDialog()}
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
          <h2>{project.name}</h2>
        </div>
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
          <button
            className="add-structure-inline-btn"
            onClick={handleImportCIF}
            disabled={isImporting}
            title="Add new structure"
          >
            +
          </button>
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
            <button
              className="calc-action-btn"
              onClick={handleRecoverPhonon}
              disabled={isRecoveringPhonon}
            >
              <span className="calc-action-icon">Rec</span>
              <span className="calc-action-label">Recover Phonon</span>
              <span className="calc-action-hint">
                {isRecoveringPhonon ? "Recovering..." : "Import completed tmp run"}
              </span>
            </button>
          </div>
        </section>

        {/* SCF Calculations */}
        {scfCalculations.length > 0 && (
          <section className="history-section">
            <h3>SCFs</h3>
            <div className="calculations-list">
              {scfCalculations.map((calc) => {
                const isPinned = pinnedCalcIds.has(calc.id);
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
                            <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
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
                          {calc.result.wall_time_seconds && (
                            <div className="detail-item">
                              <label>Wall Time</label>
                              <span>{calc.result.wall_time_seconds.toFixed(1)} s</span>
                            </div>
                          )}
                        </div>
                        <div className="detail-item parameters">
                          <label>Parameters</label>
                          <pre>{JSON.stringify(calc.parameters, null, 2)}</pre>
                        </div>
                        <div className="calc-actions">
                          <button
                            className="delete-calc-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              openDeleteCalcDialog(calc.id, calc.calc_type);
                            }}
                          >
                            Delete Calculation
                          </button>
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
                            <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
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
                        </div>
                        <div className="calc-actions">
                          {calc.result?.band_data && (
                            <button
                              className="view-bands-btn"
                              onClick={(e) => {
                                e.stopPropagation();
                                onViewBands(calc.result?.band_data, calc.result?.fermi_energy ?? null);
                              }}
                            >
                              View Bands
                            </button>
                          )}
                          <button
                            className="delete-calc-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              openDeleteCalcDialog(calc.id, calc.calc_type);
                            }}
                          >
                            Delete
                          </button>
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
                const resultData = calc.result as any;
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
                            <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
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
                        </div>
                        <div className="calc-actions">
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
                          <button
                            className="delete-calc-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              openDeleteCalcDialog(calc.id, calc.calc_type);
                            }}
                          >
                            Delete
                          </button>
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
                const summary = calc.parameters?.optimized_cell_summary as SavedCellSummary | undefined;
                const optimizedStructure = calc.parameters?.optimized_structure as SavedStructureData | undefined;
                const optimizedCell = asCellMatrix(optimizedStructure?.cell_parameters);
                const units = summary?.units || optimizedStructure?.cell_units || "angstrom";
                const supportsConventionalTransform =
                  primitiveLatticeKind === "cubic_f" ||
                  primitiveLatticeKind === "cubic_i" ||
                  primitiveLatticeKind === "cubic_p";
                let displaySummary: SavedCellSummary | null = summary
                  ? { ...summary }
                  : null;
                let cellBasisLabel = supportsConventionalTransform
                  ? "Primitive (QE)"
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

                if (cellViewMode === "conventional" && supportsConventionalTransform && optimizedCell) {
                  const convertedCell = convertPrimitiveToConventionalCell(optimizedCell, primitiveLatticeKind);
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
                  } else {
                    cellBasisLabel = "Primitive (QE)";
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

                return (
                  <div key={calc.id} className="calculation-item">
                    <div
                      className="calculation-header"
                      onClick={() =>
                        setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                      }
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
                            <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
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
                            <span>{calc.parameters?.optimized_structure?.atoms?.length || "N/A"}</span>
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
                          {calc.result?.wall_time_seconds && (
                            <div className="detail-item">
                              <label>Wall Time</label>
                              <span>{calc.result.wall_time_seconds.toFixed(1)} s</span>
                            </div>
                          )}
                        </div>
                        <div className="calc-actions">
                          <button
                            className="delete-calc-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              openDeleteCalcDialog(calc.id, calc.calc_type);
                            }}
                          >
                            Delete
                          </button>
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

      {/* Floating Settings Button */}
      <div className="floating-settings" ref={settingsRef}>
        <button
          className="floating-settings-btn"
          onClick={() => setShowSettingsMenu(!showSettingsMenu)}
          title="Project settings"
        >
          <svg width="24" height="24" viewBox="0 0 20 20" fill="currentColor">
            <path fillRule="evenodd" d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.061 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.533 1.533 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.533 1.533 0 01.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd" />
          </svg>
        </button>

        {showSettingsMenu && (
          <div className="floating-settings-menu">
            <button className="settings-menu-item danger" onClick={openDeleteDialog}>
              Delete Project
            </button>
          </div>
        )}
      </div>

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
                className="dialog-btn delete"
                onClick={handleConfirmDeleteCalc}
                disabled={isDeletingCalc}
              >
                {isDeletingCalc ? "Deleting..." : "Delete Calculation"}
              </button>
            </div>
          </div>
        </div>
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
                Type <code>{CONFIRM_TEXT}</code> to confirm:
              </label>
              <input
                type="text"
                value={deleteConfirmText}
                onChange={(e) => setDeleteConfirmText(e.target.value)}
                placeholder={CONFIRM_TEXT}
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
              className="dialog-btn delete"
              onClick={handleConfirmDelete}
              disabled={deleteConfirmText !== CONFIRM_TEXT || isDeleting}
            >
              {isDeleting ? "Deleting..." : "Delete Project"}
            </button>
          </div>
        </div>
      </div>
    );
  }
}
