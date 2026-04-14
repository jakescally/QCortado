// Band Structure Wizard - Calculate and visualize electronic band structures

import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { save } from "@tauri-apps/plugin-dialog";
import { writeTextFile } from "@tauri-apps/plugin-fs";
import {
  CrystalData,
  ELEMENT_MASSES,
  ExecutionMode,
  HpcProfile,
  SlurmResourceRequest,
} from "../lib/types";
import { BandData } from "./BandPlot";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";
import { Vec3 } from "../lib/reciprocalLattice";
import {
  analyzeCrystalSymmetry,
  buildConventionalLatticeFromCrystalData,
  SymmetryTransformResult,
} from "../lib/symmetryTransform";
import { inferQeBravaisCellFromCif } from "../lib/qeBravaisInference";
import {
  createPathCoordinateConverters,
  mapPathCoordinates,
  roundVec3,
  resolvePathTransformContext,
  sourceScfUsesPrimitiveCell,
} from "../lib/kPathTransforms";
import { resolveSavedScfStructure } from "../lib/optimizedStructure";
import {
  RhombohedralConvention,
  defaultRhombohedralConventionForSetting,
} from "../lib/brillouinZoneData";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { isPhononReadyScf } from "../lib/phononReady";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import {
  buildExecutionTarget,
  buildHpcQeInputCommandLine,
  downloadHpcCalculationArtifacts,
  defaultResourcesForProfile,
  listRemotePseudopotentials,
  resolveProfileRemoteQeBinDir,
  sampleHpcUtilization,
  saveExecutionMode,
} from "../lib/hpcConfig";
import { HpcRunSettings } from "./HpcRunSettings";

interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    total_energy: number | null;
    fermi_energy: number | null;
    raw_output?: string | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
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

function hasFullScfBundle(calc: CalculationRun, downloadedIds: Set<string>): boolean {
  if (downloadedIds.has(calc.id)) return true;
  const params = calc.parameters || {};
  if (params.artifacts_downloaded_full === true) return true;
  const mode = String(params.artifact_sync_mode || "").trim().toLowerCase();
  return mode === "full";
}

function getScfProfileId(calc: CalculationRun): string | null {
  const value = calc.parameters?.hpc_profile_id;
  if (typeof value !== "string") return null;
  const trimmed = value.trim();
  return trimmed.length > 0 ? trimmed : null;
}

// Helper to generate calculation feature tags from parameters
function getCalculationTags(
  calc: CalculationRun,
  downloadedIds?: Set<string>,
): { label: string; type: "info" | "feature" | "special" | "geometry" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" | "geometry" }[] = [];
  const params = calc.parameters || {};
  const phononReady = calc.calc_type === "scf" && isPhononReadyScf(params, calc.tags);
  const pushTag = (label: string, type: "info" | "feature" | "special" | "geometry") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

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
      }
    }
  }

  if (params.structure_source?.type === "optimization") {
    pushTag("Geometry", "geometry");
  }

  // K-points grid
  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }

  // Convergence threshold
  if (params.conv_thr) {
    const thr = params.conv_thr;
    const label = thr < 0.001 ? thr.toExponential(0) : thr.toString();
    pushTag(label, "info");
  }

  if (phononReady) {
    pushTag("Phonon-Ready", "special");
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
    if (hasFullScfBundle(calc, downloadedIds ?? new Set<string>())) {
      pushTag("Downloaded", "feature");
    }
  }

  return tags;
}

function Tooltip({ text }: { text: string }) {
  return (
    <span className="tooltip-container">
      <span className="tooltip-icon">?</span>
      <span className="tooltip-text">{text}</span>
    </span>
  );
}

const BANDS_WIZARD_SETTINGS_STORAGE_KEY = "qcortado-bands-wizard-settings-v1";

interface StoredBandWizardSettings {
  nbnd: number | "auto";
  nscfConvThrInput: string;
  nscfMixingBetaInput: string;
  nscfOccupations: "fixed" | "smearing" | "from_input" | "tetrahedra";
  nscfSmearing: "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac";
  nscfDegaussInput: string;
  nscfVerbosity: "low" | "high" | "debug";
  bandsFilbandInput: string;
  bandsLsym: boolean;
  bandsNoOverlap: boolean;
  enableProjections: boolean;
  projectionLsym: boolean;
  projectionDiagBasis: boolean;
  projectionPawproj: boolean;
  projectionFilprojInput: string;
  autoSaveLogEnabled: boolean;
  autoSaveLogPath: string;
}

function readStoredBandWizardSettings(): Partial<StoredBandWizardSettings> | null {
  if (typeof window === "undefined") return null;
  try {
    const raw = window.localStorage.getItem(BANDS_WIZARD_SETTINGS_STORAGE_KEY);
    if (!raw) return null;
    const parsed = JSON.parse(raw);
    if (!parsed || typeof parsed !== "object") return null;
    return parsed as Partial<StoredBandWizardSettings>;
  } catch {
    return null;
  }
}

function writeStoredBandWizardSettings(settings: StoredBandWizardSettings): void {
  if (typeof window === "undefined") return;
  try {
    window.localStorage.setItem(BANDS_WIZARD_SETTINGS_STORAGE_KEY, JSON.stringify(settings));
  } catch {
    // Ignore persistence failures and keep in-memory behavior.
  }
}

function parseOptionalNumber(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (trimmed.length === 0) return null;
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed)) {
    throw new Error(`Invalid ${label}.`);
  }
  return parsed;
}

function parseOptionalPositiveNumber(input: string, label: string): number | null {
  const parsed = parseOptionalNumber(input, label);
  if (parsed == null) return null;
  if (parsed <= 0) {
    throw new Error(`${label} must be positive.`);
  }
  return parsed;
}

function parseOptionalPositiveInt(input: string, label: string): number | null {
  const parsed = parseOptionalNumber(input, label);
  if (parsed == null) return null;
  if (!Number.isInteger(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive integer.`);
  }
  return parsed;
}

function sanitizeOutputFilename(raw: string, fallback: string): string {
  const trimmed = raw.trim();
  if (trimmed.length === 0) return fallback;
  const sanitized = trimmed
    .replace(/[^a-zA-Z0-9._-]/g, "_")
    .replace(/^_+|_+$/g, "");
  return sanitized.length > 0 ? sanitized : fallback;
}

function normalizeOccupations(raw: unknown): "fixed" | "smearing" | "from_input" | "tetrahedra" {
  const lowered = String(raw || "smearing").toLowerCase();
  if (lowered === "fixed") return "fixed";
  if (lowered === "from_input") return "from_input";
  if (lowered === "tetrahedra") return "tetrahedra";
  return "smearing";
}

function normalizeSmearing(
  raw: unknown,
  fallback: "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" = "marzari-vanderbilt",
): "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" {
  const lowered = String(raw || fallback).toLowerCase();
  if (lowered === "methfessel-paxton") return "methfessel-paxton";
  if (lowered === "marzari-vanderbilt") return "marzari-vanderbilt";
  if (lowered === "fermi-dirac") return "fermi-dirac";
  if (lowered === "gaussian") return "gaussian";
  return fallback;
}

const MAX_VIEWER_POINTS_PER_SEGMENT = 400;
const MAX_TOTAL_K_POINTS = 5000;

function clampInt(value: number, min: number, max: number): number {
  if (!Number.isFinite(value)) return min;
  return Math.min(max, Math.max(min, Math.round(value)));
}

function getConnectedSegmentIndices(path: KPathPoint[]): number[] {
  const indices: number[] = [];
  for (let i = 0; i < path.length - 1; i++) {
    if (path[i].npoints > 0) {
      indices.push(i);
    }
  }
  return indices;
}

function applyTotalKPoints(path: KPathPoint[], totalKPoints: number): KPathPoint[] {
  const connectedSegmentIndices = getConnectedSegmentIndices(path);
  if (path.length === 0 || connectedSegmentIndices.length === 0) {
    return path.map((point) => ({
      ...point,
      npoints: 0,
    }));
  }

  const safeTotal = clampInt(
    totalKPoints,
    connectedSegmentIndices.length,
    MAX_TOTAL_K_POINTS,
  );
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
  const rawExtras = lengths.map((length) =>
    totalLength > 0
      ? (length / totalLength) * remainingAfterBaseline
      : remainingAfterBaseline / connectedSegmentIndices.length,
  );
  const extraPoints = rawExtras.map((value) => Math.floor(value));
  const assignedExtra = extraPoints.reduce((sum, value) => sum + value, 0);
  let leftovers = remainingAfterBaseline - assignedExtra;
  const order = rawExtras
    .map((value, idx) => ({
      idx,
      frac: value - Math.floor(value),
      len: lengths[idx],
    }))
    .sort((a, b) => {
      if (b.frac !== a.frac) return b.frac - a.frac;
      if (b.len !== a.len) return b.len - a.len;
      return a.idx - b.idx;
    });

  for (let i = 0; i < order.length && leftovers > 0; i++) {
    extraPoints[order[i].idx] += 1;
    leftovers -= 1;
  }

  const segmentPoints = new Map<number, number>();
  for (let i = 0; i < connectedSegmentIndices.length; i++) {
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

interface BandStructureWizardProps {
  onBack: () => void;
  onExecutionModeChange?: (mode: ExecutionMode) => Promise<void> | void;
  onViewBands: (bandData: BandData, scfFermiEnergy: number | null) => void;
  qePath: string;
  defaultSmearing?: "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac";
  executionMode?: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "kpath" | "parameters" | "run" | "results";
const BANDS_WORK_DIR = "/tmp/qcortado_bands";

interface BandTaskPlan {
  taskLabel: string;
  taskParams: Record<string, any>;
  saveParameters: Record<string, any>;
  saveTags: string[];
}

export function BandStructureWizard({
  onBack,
  onExecutionModeChange,
  onViewBands,
  qePath,
  defaultSmearing = "marzari-vanderbilt",
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId: _cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: BandStructureWizardProps) {
  const storedBandWizardSettingsRef = useRef<Partial<StoredBandWizardSettings> | null>(readStoredBandWizardSettings());
  const shouldPreserveStoredConfigRef = useRef(Boolean(storedBandWizardSettingsRef.current));
  const storedBandWizardSettings = storedBandWizardSettingsRef.current;
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  const hasBlockingExternalTask = !isHpcMode && hasExternalRunningTask;
  const resolvedDefaultSmearing = normalizeSmearing(defaultSmearing);
  // Wizard state
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);

  // Step 1: Source SCF
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [downloadedDependencyScfIds, setDownloadedDependencyScfIds] = useState<Set<string>>(new Set());
  const [isResolvingDependency, setIsResolvingDependency] = useState(false);
  const [dependencyStatus, setDependencyStatus] = useState<string | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());

  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);

  // Step 2: K-Path (using BrillouinZoneViewer)
  const [kPath, setKPath] = useState<KPathPoint[]>([]);
  const [kPathRhombohedralConvention, setKPathRhombohedralConvention] = useState<RhombohedralConvention | undefined>(undefined);
  const [totalKPointsTarget, setTotalKPointsTarget] = useState(120);
  const [totalKPointsInput, setTotalKPointsInput] = useState("120");

  // Step 3: Parameters
  const [nbnd, setNbnd] = useState<number | "auto">(() => {
    const storedNbnd = storedBandWizardSettings?.nbnd;
    if (storedNbnd === "auto") return "auto";
    if (Number.isInteger(storedNbnd) && Number(storedNbnd) > 0) {
      return Number(storedNbnd);
    }
    return "auto";
  });
  const [nscfConvThrInput, setNscfConvThrInput] = useState(
    () => storedBandWizardSettings?.nscfConvThrInput ?? "1e-8",
  );
  const [nscfMixingBetaInput, setNscfMixingBetaInput] = useState(
    () => storedBandWizardSettings?.nscfMixingBetaInput ?? "0.7",
  );
  const [nscfOccupations, setNscfOccupations] = useState<"fixed" | "smearing" | "from_input" | "tetrahedra">(
    () => normalizeOccupations(storedBandWizardSettings?.nscfOccupations),
  );
  const [nscfSmearing, setNscfSmearing] = useState<"gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac">(
    () => normalizeSmearing(storedBandWizardSettings?.nscfSmearing, resolvedDefaultSmearing),
  );
  const [nscfDegaussInput, setNscfDegaussInput] = useState(
    () => storedBandWizardSettings?.nscfDegaussInput ?? "0.02",
  );
  const [nscfVerbosity, setNscfVerbosity] = useState<"low" | "high" | "debug">(() => {
    const storedVerbosity = storedBandWizardSettings?.nscfVerbosity;
    if (storedVerbosity === "low" || storedVerbosity === "debug") {
      return storedVerbosity;
    }
    return "high";
  });
  const [bandsFilbandInput, setBandsFilbandInput] = useState(
    () => storedBandWizardSettings?.bandsFilbandInput ?? "bands.dat",
  );
  const [bandsLsym, setBandsLsym] = useState(
    () => storedBandWizardSettings?.bandsLsym ?? true,
  );
  const [bandsNoOverlap, setBandsNoOverlap] = useState(
    () => storedBandWizardSettings?.bandsNoOverlap ?? true,
  );
  const [enableProjections, setEnableProjections] = useState(
    () => storedBandWizardSettings?.enableProjections ?? true,
  );
  const [projectionLsym, setProjectionLsym] = useState(
    () => storedBandWizardSettings?.projectionLsym ?? false,
  );
  const [projectionDiagBasis, setProjectionDiagBasis] = useState(
    () => storedBandWizardSettings?.projectionDiagBasis ?? false,
  );
  const [projectionPawproj, setProjectionPawproj] = useState(
    () => storedBandWizardSettings?.projectionPawproj ?? false,
  );
  const [projectionFilprojInput, setProjectionFilprojInput] = useState(
    () => storedBandWizardSettings?.projectionFilprojInput ?? "bands.projwfc.dat",
  );
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({
    core: true,
    nscf: true,
    post: true,
    projections: true,
    mpi: false,
  });

  // Step 4: Running
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const [outputLineCount, setOutputLineCount] = useState(0);
  const [autoSaveLogEnabled, setAutoSaveLogEnabled] = useState(
    () => storedBandWizardSettings?.autoSaveLogEnabled ?? false,
  );
  const [autoSaveLogPath, setAutoSaveLogPath] = useState(
    () => storedBandWizardSettings?.autoSaveLogPath ?? "",
  );
  const [isSavingLog, setIsSavingLog] = useState(false);
  const [logSaveStatus, setLogSaveStatus] = useState<string | null>(null);
  const visibleOutputLineCount = useMemo(() => countVisibleOutputLines(output), [output]);
  useViewportScrollLock(step === "run");
  useEffect(() => {
    if (visibleOutputLineCount > outputLineCount) {
      setOutputLineCount(visibleOutputLineCount);
    }
  }, [outputLineCount, visibleOutputLineCount]);

  const toggleSection = (section: string) => {
    setExpandedSections((prev) => ({ ...prev, [section]: !prev[section] }));
  };
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Band structure",
  });

  // Step 5: Results
  const [bandData, setBandData] = useState<BandData | null>(null);
  // Store SCF Fermi energy separately to ensure it persists
  const [scfFermiEnergy, setScfFermiEnergy] = useState<number | null>(null);
  // Show calculation output in results
  const [showOutput, setShowOutput] = useState(true);
  // Track if calculation was saved
  const [isSaved, setIsSaved] = useState(false);
  // Track calculation timing
  const [calcStartTime, setCalcStartTime] = useState<string>("");

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(
    defaultResourcesForProfile(activeHpcProfile),
  );
  const [hpcTelemetryOutput, setHpcTelemetryOutput] = useState<string>(
    "Waiting for remote job allocation...",
  );
  const [hpcTelemetrySource, setHpcTelemetrySource] = useState<string>("pending");
  const [hpcTelemetryError, setHpcTelemetryError] = useState<string | null>(null);
  const [hpcTelemetryUpdatedAt, setHpcTelemetryUpdatedAt] = useState<string | null>(null);
  const [hpcTelemetryLoading, setHpcTelemetryLoading] = useState(false);

  // Pseudopotentials (pseudopotentials list used internally for auto-selection)
  const [, setPseudopotentials] = useState<string[]>([]);
  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});
  const [symmetryTransform, setSymmetryTransform] = useState<SymmetryTransformResult | null>(null);
  const [symmetryError, setSymmetryError] = useState<string | null>(null);

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [isHpcMode, activeHpcProfile?.id, activeHpcProfile?.resource_mode]);

  useEffect(() => {
    writeStoredBandWizardSettings({
      nbnd,
      nscfConvThrInput,
      nscfMixingBetaInput,
      nscfOccupations,
      nscfSmearing,
      nscfDegaussInput,
      nscfVerbosity,
      bandsFilbandInput,
      bandsLsym,
      bandsNoOverlap,
      enableProjections,
      projectionLsym,
      projectionDiagBasis,
      projectionPawproj,
      projectionFilprojInput,
      autoSaveLogEnabled,
      autoSaveLogPath,
    });
  }, [
    nbnd,
    nscfConvThrInput,
    nscfMixingBetaInput,
    nscfOccupations,
    nscfSmearing,
    nscfDegaussInput,
    nscfVerbosity,
    bandsFilbandInput,
    bandsLsym,
    bandsNoOverlap,
    enableProjections,
    projectionLsym,
    projectionDiagBasis,
    projectionPawproj,
    projectionFilprojInput,
    autoSaveLogEnabled,
    autoSaveLogPath,
  ]);

  // Helper functions
  function getBaseElement(symbol: string): string {
    return symbol.replace(/[\d+-]+$/, "");
  }

  function pseudoMatchesElement(filename: string, element: string): boolean {
    const lowerFile = filename.toLowerCase();
    const lowerEl = getBaseElement(element).toLowerCase();
    return (
      lowerFile.startsWith(lowerEl + ".") ||
      lowerFile.startsWith(lowerEl + "_") ||
      lowerFile.startsWith(lowerEl + "-")
    );
  }

  const applyScfDefaults = useCallback((scf: CalculationRun) => {
    const params = scf.parameters || {};
    const convThr = Number(params.conv_thr);
    if (Number.isFinite(convThr) && convThr > 0) {
      setNscfConvThrInput(String(convThr));
    }

    const mixingBeta = Number(params.mixing_beta);
    if (Number.isFinite(mixingBeta) && mixingBeta > 0) {
      setNscfMixingBetaInput(String(mixingBeta));
    }

    setNscfOccupations(normalizeOccupations(params.occupations));
    setNscfSmearing(normalizeSmearing(params.smearing, resolvedDefaultSmearing));

    const degauss = Number(params.degauss);
    if (Number.isFinite(degauss) && degauss > 0) {
      setNscfDegaussInput(String(degauss));
    } else {
      setNscfDegaussInput("");
    }

    const verbosityRaw = String(params.verbosity || "high").toLowerCase();
    if (verbosityRaw === "low" || verbosityRaw === "debug") {
      setNscfVerbosity(verbosityRaw);
    } else {
      setNscfVerbosity("high");
    }

    const sourceNbnd = Number(params.nbnd);
    if (Number.isInteger(sourceNbnd) && sourceNbnd > 0) {
      setNbnd(sourceNbnd);
    } else {
      setNbnd("auto");
    }
  }, [resolvedDefaultSmearing]);

  const selectSourceScf = useCallback((scf: CalculationRun) => {
    setSelectedScf(scf);
    setDependencyStatus(null);
    setScfFermiEnergy(scf.result?.fermi_energy ?? null);
    if (!shouldPreserveStoredConfigRef.current) {
      applyScfDefaults(scf);
    }
  }, [applyScfDefaults]);

  const selectedScfDependencyBlocked = useMemo(() => {
    if (isHpcMode || !selectedScf) return false;
    if (!isHpcCalculation(selectedScf)) return false;
    return !hasFullScfBundle(selectedScf, downloadedDependencyScfIds);
  }, [isHpcMode, selectedScf, downloadedDependencyScfIds]);

  async function handleSwitchToHpcMode() {
    setIsResolvingDependency(true);
    setDependencyStatus("Switching execution mode to HPC...");
    setError(null);
    try {
      if (onExecutionModeChange) {
        await onExecutionModeChange("hpc");
      } else {
        await saveExecutionMode("hpc");
      }
      setDependencyStatus("Execution mode switched to HPC. Review HPC settings and submit remotely.");
    } catch (e) {
      setError(`Failed to switch execution mode: ${e}`);
      setDependencyStatus(null);
    } finally {
      setIsResolvingDependency(false);
    }
  }

  async function handleDownloadDependencyBundle() {
    if (!selectedScf) return;
    setIsResolvingDependency(true);
    setDependencyStatus("Downloading full SCF bundle from remote...");
    setError(null);
    try {
      const report = await downloadHpcCalculationArtifacts(
        projectId,
        selectedScf.id,
        getScfProfileId(selectedScf),
        true,
      );
      setDownloadedDependencyScfIds((prev) => {
        const next = new Set(prev);
        next.add(selectedScf.id);
        return next;
      });
      setSelectedScf((prev) => {
        if (!prev || prev.id !== selectedScf.id) return prev;
        return {
          ...prev,
          parameters: {
            ...(prev.parameters || {}),
            artifacts_downloaded_full: true,
            artifact_sync_mode: "full",
            remote_storage_bytes: report.downloaded_bytes + report.skipped_bytes,
          },
        };
      });
      setDependencyStatus("Full SCF bundle downloaded. Local run is now available.");
    } catch (e) {
      setError(`Failed to download full SCF bundle: ${e}`);
      setDependencyStatus(null);
    } finally {
      setIsResolvingDependency(false);
    }
  }

  // Load MPI info and pseudopotentials
  useEffect(() => {
    async function init() {
      try {
        // Check MPI
        const count = await invoke<number>("get_cpu_count");
        const safeCount = Math.max(1, Math.floor(count));
        setCpuCount(safeCount);
        const defaults = await loadGlobalMpiDefaults(safeCount);

        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);
        setMpiEnabled(available ? defaults.enabled : false);
        setMpiProcs(defaults.nprocs);

        // Load pseudopotentials
        const pseudoDir = isHpcMode
          ? (activeHpcProfile?.remote_pseudo_dir || "")
          : qePath.replace(/\/bin\/?$/, "/pseudo");
        const pseudos = isHpcMode
          ? await listRemotePseudopotentials(pseudoDir, activeHpcProfile?.id ?? null)
          : await invoke<string[]>("list_pseudopotentials", { pseudoDir });
        setPseudopotentials(pseudos);

        // Auto-select pseudopotentials for elements in the structure
        const elements = [...new Set(crystalData.atom_sites.map(s => getBaseElement(s.type_symbol)))];
        const selected: Record<string, string> = {};
        for (const el of elements) {
          const match = pseudos.find(p => pseudoMatchesElement(p, el));
          if (match) {
            selected[el] = match;
          }
        }
        setSelectedPseudos(selected);
      } catch (e) {
        console.error("Failed to initialize:", e);
      }
    }
    init();
  }, [qePath, crystalData, isHpcMode, activeHpcProfile?.id, activeHpcProfile?.remote_pseudo_dir]);

  useEffect(() => {
    let cancelled = false;
    setSymmetryError(null);
    analyzeCrystalSymmetry(crystalData)
      .then((result) => {
        if (cancelled) return;
        setSymmetryTransform(result);
      })
      .catch((err) => {
        if (cancelled) return;
        console.error("Failed to analyze symmetry with spglib:", err);
        setSymmetryTransform(null);
        setSymmetryError(String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [crystalData]);

  useEffect(() => {
    setKPathRhombohedralConvention(undefined);
  }, [crystalData]);

  useEffect(() => {
    if (!isHpcMode || step !== "run") {
      return;
    }

    const taskIsRunning = isRunning || activeTask?.status === "running";
    if (!taskIsRunning) {
      return;
    }

    const profileId = activeHpcProfile?.id ?? null;
    if (!profileId) {
      setHpcTelemetryError("No active HPC profile selected.");
      setHpcTelemetrySource("unavailable");
      return;
    }

    if (activeTask?.hpc?.backend && activeTask.hpc.backend !== "hpc") {
      return;
    }

    const remoteJobId = activeTask?.hpc?.remote_job_id?.trim() || "";
    const remoteNode = activeTask?.hpc?.remote_node?.trim() || "";
    if (!remoteJobId) {
      setHpcTelemetryLoading(false);
      setHpcTelemetryError(null);
      setHpcTelemetrySource("pending");
      setHpcTelemetryOutput("Waiting for remote job allocation...");
      return;
    }

    let cancelled = false;
    let timeoutId: number | null = null;

    const pollTelemetry = async () => {
      if (cancelled) return;
      setHpcTelemetryLoading(true);
      try {
        const sample = await sampleHpcUtilization(profileId, remoteJobId, remoteNode || null);
        if (cancelled) return;
        setHpcTelemetryOutput(sample.output || "No telemetry output received from remote host.");
        setHpcTelemetrySource(sample.source || "unknown");
        setHpcTelemetryUpdatedAt(sample.captured_at || new Date().toISOString());
        setHpcTelemetryError(null);
      } catch (e) {
        if (cancelled) return;
        setHpcTelemetrySource("error");
        setHpcTelemetryError(String(e));
      } finally {
        if (cancelled) return;
        setHpcTelemetryLoading(false);
        timeoutId = window.setTimeout(() => {
          void pollTelemetry();
        }, 5000);
      }
    };

    void pollTelemetry();

    return () => {
      cancelled = true;
      if (timeoutId !== null) {
        window.clearTimeout(timeoutId);
      }
    };
  }, [
    isHpcMode,
    step,
    isRunning,
    activeTask?.status,
    activeTask?.hpc?.backend,
    activeTask?.hpc?.remote_job_id,
    activeTask?.hpc?.remote_node,
    activeHpcProfile?.id,
  ]);

  // Reconnect to a running/completed background task
  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) {
      taskContext.reconnectToTask(activeTaskId);
      return;
    }

    setIsRunning(task.status === "running");
    if (task.outputLineCount > 0 || task.status !== "running") {
      setOutput(task.outputText);
      setOutputLineCount(task.outputLineCount);
    }
    setProgress(task.progress);
    setCalcStartTime(task.startedAt);

    if (task.status === "completed" && task.result) {
      setBandData(task.result as any);
      setStep("results");
    } else if (task.status === "failed" || task.status === "cancelled") {
      setError(task.error || "Task failed");
    } else {
      setStep("run");
    }
  }, [activeTaskId, taskContext.getTask(activeTaskId ?? "")?.status]);

  // Sync output/progress from active task in real-time
  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) return;

    if (task.outputLineCount > 0) {
      setOutput(task.outputText);
      setOutputLineCount(task.outputLineCount);
    }
    setProgress(task.progress);
    setIsRunning(task.status === "running");
  }, [
    activeTaskId,
    taskContext.getTask(activeTaskId ?? "")?.output.length,
    taskContext.getTask(activeTaskId ?? "")?.status,
  ]);

  // Handle k-path changes from the BZ viewer
  const handleKPathChange = useCallback((newPath: KPathPoint[]) => {
    setKPath(
      applyTotalKPoints(newPath, totalKPointsTarget),
    );
  }, [totalKPointsTarget]);

  useEffect(() => {
    setKPath((prevPath) =>
      applyTotalKPoints(prevPath, totalKPointsTarget),
    );
  }, [totalKPointsTarget]);

  const kPathSegmentCount = useMemo(
    () => getConnectedSegmentIndices(kPath).length,
    [kPath],
  );
  const totalKPoints = useMemo(
    () => kPath.reduce((sum, point) => sum + point.npoints, 0),
    [kPath],
  );
  const minimumTotalKPoints = Math.max(1, kPathSegmentCount);
  const viewerPointsPerSegment = clampInt(
    totalKPointsTarget / minimumTotalKPoints,
    1,
    MAX_VIEWER_POINTS_PER_SEGMENT,
  );
  const hpcCommandLines = useMemo(() => {
    const qeBinDir = resolveProfileRemoteQeBinDir(activeHpcProfile, hpcResources.resource_type);
    const lines = [
      "cd \"$SLURM_SUBMIT_DIR\"",
      `QE_BIN="${qeBinDir}"`,
      buildHpcQeInputCommandLine(activeHpcProfile, "pw.x", "bands.in", "bands.out"),
      buildHpcQeInputCommandLine(activeHpcProfile, "bands.x", "bands_pp.in", "bands_pp.out"),
    ];
    if (enableProjections) {
      lines.push(buildHpcQeInputCommandLine(activeHpcProfile, "projwfc.x", "projwfc.in", "projwfc.out"));
    }
    return lines;
  }, [activeHpcProfile, enableProjections, hpcResources.resource_type]);

  useEffect(() => {
    setTotalKPointsInput(String(totalKPointsTarget));
  }, [totalKPointsTarget]);

  const commitTotalKPointsInput = useCallback(() => {
    const parsed = Number.parseInt(totalKPointsInput.trim(), 10);
    const fallback = Number.isFinite(totalKPointsTarget)
      ? totalKPointsTarget
      : Math.max(120, minimumTotalKPoints);
    const committed = Number.isFinite(parsed)
      ? clampInt(parsed, minimumTotalKPoints, MAX_TOTAL_K_POINTS)
      : clampInt(fallback, minimumTotalKPoints, MAX_TOTAL_K_POINTS);
    setTotalKPointsTarget(committed);
    setTotalKPointsInput(String(committed));
  }, [minimumTotalKPoints, totalKPointsInput, totalKPointsTarget]);

  const buildBandTaskPlan = async (): Promise<BandTaskPlan> => {
    if (!selectedScf?.result) {
      throw new Error("No source SCF calculation selected");
    }

    if (kPath.length < 2) {
      throw new Error("Please select at least 2 points for the k-path");
    }

    const parsedConvThr = parseOptionalPositiveNumber(nscfConvThrInput, "NSCF convergence threshold");
    if (parsedConvThr == null) {
      throw new Error("Please enter a valid positive NSCF convergence threshold.");
    }

    const parsedMixingBeta = parseOptionalPositiveNumber(nscfMixingBetaInput, "NSCF mixing beta");
    if (parsedMixingBeta == null) {
      throw new Error("Please enter a valid positive NSCF mixing beta.");
    }
    if (parsedMixingBeta > 1.0) {
      throw new Error("NSCF mixing beta should typically be in the range (0, 1].");
    }

    const parsedDegauss = nscfOccupations === "smearing"
      ? parseOptionalPositiveNumber(nscfDegaussInput, "NSCF degauss")
      : null;
    if (nscfOccupations === "smearing" && parsedDegauss == null) {
      throw new Error("Please provide a positive degauss value when using smearing occupations.");
    }

    const manualNbnd = nbnd === "auto"
      ? null
      : parseOptionalPositiveInt(String(nbnd), "number of bands");

    const bandsFilband = sanitizeOutputFilename(bandsFilbandInput, "bands.dat");
    const projectionFilproj = sanitizeOutputFilename(projectionFilprojInput, "bands.projwfc.dat");

    // Get the SCF parameters for inherited defaults and tags.
    const scfParams = selectedScf.parameters || {};
    const savedPseudoMap = (scfParams.selected_pseudos && typeof scfParams.selected_pseudos === "object")
      ? scfParams.selected_pseudos as Record<string, string>
      : {};

    // Get unique elements from crystal data
    const elements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];

    // Resolve pseudopotentials, preferring the source SCF mapping.
    const resolvedPseudos: Record<string, string> = {};
    for (const el of elements) {
      const savedPseudo = savedPseudoMap[el];
      const detectedPseudo = selectedPseudos[el];
      const resolvedPseudo = (typeof savedPseudo === "string" && savedPseudo.length > 0)
        ? savedPseudo
        : detectedPseudo;
      if (!resolvedPseudo) {
        throw new Error(`No pseudopotential selected for element ${el}`);
      }
      resolvedPseudos[el] = resolvedPseudo;
    }

    // Build the full calculation config from crystalData
    // Use the same prefix as the source SCF so we can read its .save directory
    // SCFWizard uses "qcortado_scf" as the prefix
    const scfPrefix = scfParams.prefix || "qcortado_scf";
    const pseudoDir = isHpcMode
      ? (activeHpcProfile?.remote_pseudo_dir || "")
      : qePath.replace(/\/bin\/?$/, "/pseudo");
    if (!pseudoDir.trim()) {
      throw new Error(
        isHpcMode
          ? "Remote pseudopotential directory is not configured for the active HPC profile."
          : "Local pseudopotential directory is not configured.",
      );
    }

    const ecutwfcValue = Number(scfParams.ecutwfc);
    const ecutwfc = Number.isFinite(ecutwfcValue) && ecutwfcValue > 0 ? ecutwfcValue : 40;
    const ecutrhoValue = Number(scfParams.ecutrho);
    const ecutrho = Number.isFinite(ecutrhoValue) && ecutrhoValue > 0
      ? ecutrhoValue
      : ecutwfc * 8;
    const sourceNspin = Number(scfParams.nspin);
    const nspin = Number.isFinite(sourceNspin) && sourceNspin > 0 ? sourceNspin : 1;
    const lspinorb = Boolean(scfParams.lspinorb);
    const noncolin = nspin === 4 || Boolean(scfParams.noncolin) || lspinorb;

    let resolvedSymmetry = symmetryTransform;
    let resolvedSymmetryError = symmetryError;
    if (!resolvedSymmetry) {
      try {
        resolvedSymmetry = await analyzeCrystalSymmetry(crystalData);
        setSymmetryTransform(resolvedSymmetry);
        setSymmetryError(null);
        resolvedSymmetryError = null;
      } catch (err) {
        resolvedSymmetryError = String(err);
        setSymmetryError(resolvedSymmetryError);
      }
    }

    const structureForNscf = resolveSavedScfStructure(scfParams);
    const isOptimizedSourceScf =
      String(scfParams.cell_representation || "").toLowerCase() === "optimized_source";
    if (isOptimizedSourceScf && !structureForNscf) {
      throw new Error(
        "Selected SCF was run from an optimized structure, but its saved structure metadata is missing. Re-run the SCF from the optimized source and try again.",
      );
    }
    const sourceUsesPrimitive = sourceScfUsesPrimitiveCell(
      scfParams,
      crystalData.atom_sites.length,
    );
    const canUseSymmetryPrimitive =
      sourceUsesPrimitive &&
      resolvedSymmetry !== null &&
      resolvedSymmetry.standardizedPrimitiveAtoms.length > 0;
    if (sourceUsesPrimitive && !canUseSymmetryPrimitive) {
      throw new Error(
        "Selected SCF was run in a primitive cell, but symmetry conversion data is unavailable. Re-run the SCF or refresh the structure metadata.",
      );
    }

    const conventionalLattice = buildConventionalLatticeFromCrystalData(crystalData);

    const species = elements.map((el) => ({
      symbol: el,
      mass: ELEMENT_MASSES[el] || 1.0,
      pseudopotential: resolvedPseudos[el],
    }));

    const context = resolvePathTransformContext(crystalData, resolvedSymmetry);
    const effectiveRhombohedralConvention = kPathRhombohedralConvention ??
      defaultRhombohedralConventionForSetting(context.rhombohedralSetting ?? null);
    const converters = createPathCoordinateConverters(context, resolvedSymmetry);

    let baseCalculation;
    let transformedKPath: KPathPoint[];

    if (isOptimizedSourceScf && structureForNscf) {
      baseCalculation = {
        calculation: "scf",
        prefix: scfPrefix,
        outdir: "./tmp",
        pseudo_dir: pseudoDir,
        system: {
          ibrav: "free",
          celldm: null,
          cell_parameters: structureForNscf.cell_parameters,
          cell_units: structureForNscf.cell_units || "angstrom",
          species,
          atoms: structureForNscf.atoms.map((atom) => ({
            symbol: atom.symbol,
            position: roundVec3(atom.position),
            if_pos: [true, true, true],
          })),
          position_units: "crystal",
          ecutwfc,
          ecutrho,
          nspin,
          noncolin,
          lspinorb,
          occupations: nscfOccupations,
          smearing: nscfSmearing,
          degauss: parsedDegauss,
        },
        kpoints: { type: "gamma" },
        conv_thr: parsedConvThr,
        mixing_beta: parsedMixingBeta,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: nscfVerbosity,
      };
      transformedKPath = mapPathCoordinates(
        kPath,
        canUseSymmetryPrimitive
          ? converters.toSymmetryPrimitiveCoords
          : converters.toInputConventionalCoords,
      ).map((point) => ({
        label: point.label,
        coords: point.coords as Vec3,
        npoints: point.npoints,
      }));
    } else if (canUseSymmetryPrimitive && resolvedSymmetry) {
      const inferredBravais = inferQeBravaisCellFromCif(crystalData, resolvedSymmetry);
      baseCalculation = {
        calculation: "scf",
        prefix: scfPrefix,
        outdir: "./tmp",
        pseudo_dir: pseudoDir,
        system: {
          ibrav: inferredBravais?.ibrav ?? "free",
          celldm: inferredBravais?.celldm ?? null,
          cell_parameters: inferredBravais ? null : resolvedSymmetry.standardizedPrimitiveLattice,
          cell_units: inferredBravais ? null : "angstrom",
          species,
          atoms: (inferredBravais?.atoms ?? resolvedSymmetry.standardizedPrimitiveAtoms).map((atom) => ({
            symbol: atom.symbol,
            position: roundVec3(atom.position),
            if_pos: [true, true, true],
          })),
          position_units: "crystal",
          ecutwfc,
          ecutrho,
          nspin,
          noncolin,
          lspinorb,
          occupations: nscfOccupations,
          smearing: nscfSmearing,
          degauss: parsedDegauss,
        },
        kpoints: { type: "gamma" },
        conv_thr: parsedConvThr,
        mixing_beta: parsedMixingBeta,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: nscfVerbosity,
      };
      transformedKPath = mapPathCoordinates(kPath, converters.toSymmetryPrimitiveCoords).map((point) => ({
        label: point.label,
        coords: point.coords as Vec3,
        npoints: point.npoints,
      }));
    } else {
      baseCalculation = {
        calculation: "scf",
        prefix: scfPrefix,
        outdir: "./tmp",
        pseudo_dir: pseudoDir,
        system: {
          ibrav: "free",
          celldm: null,
          cell_parameters: conventionalLattice,
          cell_units: "angstrom",
          species,
          atoms: crystalData.atom_sites.map((site) => ({
            symbol: getBaseElement(site.type_symbol),
            position: roundVec3([site.fract_x, site.fract_y, site.fract_z]),
            if_pos: [true, true, true],
          })),
          position_units: "crystal",
          ecutwfc,
          ecutrho,
          nspin,
          noncolin,
          lspinorb,
          occupations: nscfOccupations,
          smearing: nscfSmearing,
          degauss: parsedDegauss,
        },
        kpoints: { type: "gamma" },
        conv_thr: parsedConvThr,
        mixing_beta: parsedMixingBeta,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: nscfVerbosity,
      };

      transformedKPath = mapPathCoordinates(kPath, converters.toInputConventionalCoords).map((point) => ({
        label: point.label,
        coords: point.coords as Vec3,
        npoints: point.npoints,
      }));
    }

    const bandCellRepresentation = isOptimizedSourceScf && structureForNscf
      ? "optimized_source"
      : canUseSymmetryPrimitive
        ? "primitive_spglib"
        : "conventional_input";

    const taskLabel = `Bands - ${crystalData?.formula_sum || ""}`;
    const taskParams = {
      config: {
        base_calculation: baseCalculation,
        k_path: transformedKPath,
        nbnd: manualNbnd,
        project_id: projectId,
        scf_calc_id: selectedScf.id,
        bands_x: {
          filband: bandsFilband,
          lsym: bandsLsym,
          no_overlap: bandsNoOverlap,
        },
        projections: {
          enabled: enableProjections,
          filproj: projectionFilproj,
          lsym: projectionLsym,
          diag_basis: projectionDiagBasis,
          pawproj: projectionPawproj,
        },
      },
      workingDir: BANDS_WORK_DIR,
      mpiConfig: !isHpcMode && mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      executionTarget: buildExecutionTarget(
        executionMode,
        activeHpcProfile?.id ?? null,
        isHpcMode ? hpcResources : null,
        false,
      ),
    };

    const pathString = kPath.map((point) => point.label).join(" → ");
    const saveParameters = {
      source_scf_id: selectedScf.id,
      k_path: pathString,
      k_path_sampling_mode: "total",
      total_k_points_target: totalKPoints,
      total_k_points: null,
      n_bands: manualNbnd,
      n_bands_requested: manualNbnd,
      nscf_conv_thr: parsedConvThr,
      nscf_mixing_beta: parsedMixingBeta,
      nscf_occupations: nscfOccupations,
      nscf_smearing: nscfSmearing,
      nscf_degauss: parsedDegauss,
      nscf_verbosity: nscfVerbosity,
      bands_x_filband: bandsFilband,
      bands_x_lsym: bandsLsym,
      bands_x_no_overlap: bandsNoOverlap,
      // Inherit SCF parameters for tags
      ecutwfc: scfParams.ecutwfc,
      nspin: scfParams.nspin,
      lspinorb: scfParams.lspinorb,
      lda_plus_u: scfParams.lda_plus_u,
      vdw_corr: scfParams.vdw_corr,
      fat_bands_requested: enableProjections,
      projection_filproj: projectionFilproj,
      projection_lsym: projectionLsym,
      projection_diag_basis: projectionDiagBasis,
      projection_pawproj: projectionPawproj,
      log_auto_save_enabled: autoSaveLogEnabled,
      log_auto_save_path: autoSaveLogEnabled ? autoSaveLogPath.trim() || null : null,
      scf_fermi_energy: selectedScf.result?.fermi_energy ?? scfFermiEnergy,
      cell_representation: bandCellRepresentation,
      k_path_convention: context.centering === "R" ? effectiveRhombohedralConvention : null,
      k_path_rhombohedral_setting: context.centering === "R" ? (context.rhombohedralSetting ?? null) : null,
      symmetry_spacegroup: resolvedSymmetry?.spacegroupNumber ?? null,
      symmetry_hall_number: resolvedSymmetry?.hallNumber ?? null,
      primitive_to_input_reciprocal: resolvedSymmetry?.primitiveToInputReciprocal ?? null,
      symmetry_error: resolvedSymmetryError,
    };

    return {
      taskLabel,
      taskParams,
      saveParameters,
      saveTags: enableProjections ? ["Proj"] : [],
    };
  };

  // Run the calculation
  const runCalculation = async () => {
    if (selectedScfDependencyBlocked) {
      setError("Selected SCF was computed remotely and needs a full local bundle for local execution.");
      return;
    }
    if (autoSaveLogEnabled && autoSaveLogPath.trim().length === 0) {
      setError("Set a local log file path or disable auto-save log before running.");
      return;
    }
    if (hasBlockingExternalTask) {
      setError("Another local task is currently running. Queue this task or wait for completion.");
      return;
    }
    setIsRunning(true);
    setOutput("");
    setOutputLineCount(0);
    setError(null);
    setBandData(null);
    setIsSaved(false);
    setLogSaveStatus(null);
    setProgress(defaultProgressState("Band structure"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");
    if (isHpcMode) {
      setHpcTelemetryOutput("Waiting for remote job allocation...");
      setHpcTelemetrySource("pending");
      setHpcTelemetryError(null);
      setHpcTelemetryUpdatedAt(null);
      setHpcTelemetryLoading(false);
    }

    try {
      const plan = await buildBandTaskPlan();

      const taskId = await taskContext.startTask("bands", plan.taskParams, plan.taskLabel);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Calculation failed");
      }

      const result = finalTask.result as BandData;
      const outputContent = finalTask.output.join("\n");
      const endTime = new Date().toISOString();
      await persistLogToConfiguredPath(outputContent);
      const hpcSaveParams = (isHpcMode || finalTask.hpc.backend === "hpc")
        ? {
          execution_backend: "hpc",
          hpc_profile_id: activeHpcProfile?.id ?? null,
          remote_job_id: finalTask.hpc.remote_job_id ?? null,
          scheduler_state: finalTask.hpc.scheduler_state ?? null,
          remote_node: finalTask.hpc.remote_node ?? null,
          remote_workdir: finalTask.hpc.remote_workdir ?? null,
          remote_project_path: finalTask.hpc.remote_project_path ?? null,
          remote_storage_bytes: finalTask.hpc.remote_storage_bytes ?? null,
        }
        : {};
      setBandData(result);
      setStep("results");
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      // Auto-save the band calculation to the project
      try {
        await invoke("save_calculation", {
          projectId,
          cifId: _cifId,
          calcData: {
            calc_type: "bands",
            parameters: {
              ...plan.saveParameters,
              total_k_points: result.n_kpoints,
              n_bands: result.n_bands,
              ...hpcSaveParams,
            },
            result: {
              converged: true,
              total_energy: null,
              fermi_energy: selectedScf?.result?.fermi_energy ?? scfFermiEnergy,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: outputContent,
              // Store band data for later viewing
              band_data: result,
            },
            started_at: startTime,
            completed_at: endTime,
            input_content: "", // TODO: store bands input
            output_content: outputContent,
            tags: plan.saveTags,
          },
          workingDir: finalTask.hpc.local_sync_dir ?? BANDS_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save band calculation:", saveError);
        setError(`Failed to auto-save band calculation: ${saveError}`);
      }
    } catch (e) {
      const errorText = String(e);
      setError(errorText);
      const fallbackLog = output.trim().length > 0
        ? `${output}\nError: ${errorText}\n`
        : `Error: ${errorText}\n`;
      await persistLogToConfiguredPath(fallbackLog);
      setOutput((prev) => prev + `\nError: ${errorText}\n`);
      setProgress((prev) => ({
        ...prev,
        status: "error",
        percent: null,
        phase: "Error",
      }));
    } finally {
      setIsRunning(false);
    }
  };

  const queueCalculation = async () => {
    if (isHpcMode) {
      setError("Queueing is unavailable in HPC mode. Submit directly to Andromeda.");
      return;
    }
    if (autoSaveLogEnabled) {
      setError("Auto-save log currently applies to direct runs only. Disable it before queueing.");
      return;
    }
    if (selectedScfDependencyBlocked) {
      setError("Selected SCF was computed remotely and needs a full local bundle for local execution.");
      return;
    }
    try {
      const plan = await buildBandTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "bands",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId: _cifId,
          workingDir: BANDS_WORK_DIR,
          calcType: "bands",
          parameters: plan.saveParameters,
          tags: plan.saveTags,
          inputContent: "",
        },
      );
    } catch (e) {
      setError(String(e));
    }
  };

  const handleChooseAutoLogPath = useCallback(async () => {
    setLogSaveStatus(null);
    try {
      const timestamp = new Date().toISOString().replace(/:/g, "-").replace(/\..+$/, "");
      const defaultPath = autoSaveLogPath.trim().length > 0
        ? autoSaveLogPath.trim()
        : `bands-log-${isHpcMode ? "hpc" : "local"}-${timestamp}.log`;
      const destinationPath = await save({
        title: "Choose Band Structure Log Path",
        defaultPath,
        filters: [{ name: "Log file", extensions: ["log", "txt"] }],
      });
      if (!destinationPath || Array.isArray(destinationPath)) {
        return;
      }
      setAutoSaveLogPath(destinationPath);
      setAutoSaveLogEnabled(true);
    } catch (e) {
      console.error("Failed to choose band-structure log path:", e);
      setLogSaveStatus(`Failed to choose log path: ${e}`);
    }
  }, [autoSaveLogPath, isHpcMode]);

  const persistLogToConfiguredPath = useCallback(async (logText: string) => {
    if (!autoSaveLogEnabled) return;
    const destinationPath = autoSaveLogPath.trim();
    if (!destinationPath) {
      setLogSaveStatus("Auto-save log is enabled, but no log file path is set.");
      return;
    }

    setIsSavingLog(true);
    setLogSaveStatus(null);
    try {
      const content = logText.endsWith("\n")
        ? logText
        : `${logText}\n`;
      await writeTextFile(destinationPath, content);
      setLogSaveStatus(`Run log saved to ${destinationPath}.`);
    } catch (e) {
      console.error("Failed to save configured band-structure log:", e);
      const errorText = String(e);
      const permissionHint = /scope|denied|forbidden|not allowed|permission/i.test(errorText)
        ? " Re-select the path with 'Choose...' to refresh file access permissions."
        : "";
      const statusText = `Failed to save log: ${errorText}${permissionHint}`.trim();
      setLogSaveStatus(statusText);
      setError(`Run completed, but log save failed. ${errorText}${permissionHint}`.trim());
    } finally {
      setIsSavingLog(false);
    }
  }, [autoSaveLogEnabled, autoSaveLogPath]);

  // Render current step
  const renderStep = () => {
    switch (step) {
      case "source":
        return renderSourceStep();
      case "kpath":
        return renderKPathStep();
      case "parameters":
        return renderParametersStep();
      case "run":
        return renderRunStep();
      case "results":
        return renderResultsStep();
    }
  };

  // Step 1: Select source SCF
  const renderSourceStep = () => {
    const validScfs = scfCalculations.filter(c => c.calc_type === "scf" && c.result?.converged);
    const sortedScfs = sortScfByMode(validScfs, scfSortMode);

    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Band structure calculations require a completed SCF calculation.
            Please run an SCF calculation first, then return here.
          </p>
          <button className="primary-button" onClick={onBack}>
            Back to Dashboard
          </button>
        </div>
      );
    }

    return (
      <div className="wizard-step source-step">
        <div className="source-step-header">
          <h3>Select Source SCF Calculation</h3>
          <div className="source-sort-control">
            <label htmlFor="bands-scf-sort">Sort SCFs</label>
            <select
              id="bands-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose the SCF calculation to use as the starting point for band structure.
        </p>

        <div className="scf-list">
          {sortedScfs.map((scf) => (
            <div
              key={scf.id}
              className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""}`}
              onClick={() => selectSourceScf(scf)}
            >
              <div className="scf-option-header">
                <input
                  type="radio"
                  checked={selectedScf?.id === scf.id}
                  onChange={() => selectSourceScf(scf)}
                />
                <span className="scf-date">
                  {new Date(scf.started_at).toLocaleDateString()}
                </span>
              </div>
              <div className="scf-details">
                <span>E = {scf.result?.total_energy?.toFixed(6)} Ry</span>
                {scf.result?.fermi_energy && (
                  <span>E_F = {scf.result.fermi_energy.toFixed(3)} eV</span>
                )}
              </div>
              <div className="calc-tags">
                {getCalculationTags(scf, downloadedDependencyScfIds).map((tag, i) => (
                  <span
                    key={i}
                    className={`calc-tag calc-tag-${tag.type}${tag.label.trim().toUpperCase() === "HPC" ? " calc-tag-hpc" : ""}${tag.label.trim().toUpperCase() === "DOWNLOADED" ? " calc-tag-downloaded" : ""}`}
                  >
                    {tag.label}
                  </span>
                ))}
              </div>
            </div>
          ))}
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Cancel
          </button>
          <button
            className="primary-button"
            disabled={!selectedScf}
            onClick={() => setStep("kpath")}
          >
            Next: K-Path
          </button>
        </div>
      </div>
    );
  };

  // Step 2: K-Path selection with 3D Brillouin Zone viewer
  const renderKPathStep = () => {
    return (
      <div className="wizard-step kpath-step">
        <h3>Select K-Path</h3>

        <div className="crystal-info">
          {crystalData.space_group_HM && (
            <div className="info-row">
              <span className="label">Space Group:</span>
              <span className="value">
                {crystalData.space_group_HM}
                {crystalData.space_group_IT_number && ` (#${crystalData.space_group_IT_number})`}
              </span>
            </div>
          )}
        </div>

        <BrillouinZoneViewer
          crystalData={crystalData}
          onPathChange={handleKPathChange}
          initialPath={kPath}
          pointsPerSegment={viewerPointsPerSegment}
          symmetryTransform={symmetryTransform}
          rhombohedralConvention={kPathRhombohedralConvention}
          onRhombohedralConventionChange={setKPathRhombohedralConvention}
        />

        <div className="kpath-sampling-panel">
          <div className="kpath-sampling-header">
            <div>
              <h4>K-Path Sampling</h4>
              <p>K-points are distributed along the full path by segment length.</p>
            </div>
            <span className="kpath-sampling-summary">
              {totalKPoints} total k-points
            </span>
          </div>

          <label className="kpath-sampling-input">
            <span>Total k-points</span>
            <input
              type="number"
              min={minimumTotalKPoints}
              max={MAX_TOTAL_K_POINTS}
              value={totalKPointsInput}
              onChange={(e) => {
                setTotalKPointsInput(e.target.value);
              }}
              onBlur={commitTotalKPointsInput}
              onKeyDown={(e) => {
                if (e.key === "Enter") {
                  e.preventDefault();
                  commitTotalKPointsInput();
                }
              }}
            />
          </label>

          <p className="kpath-sampling-note">
            {kPathSegmentCount > 0
              ? `Evenly distributed by segment length across ${kPathSegmentCount} path segment${kPathSegmentCount === 1 ? "" : "s"}.`
              : "Add at least 2 points to distribute k-points along the path."}
          </p>
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button
            className="primary-button"
            disabled={kPath.length < 2}
            onClick={() => setStep("parameters")}
          >
            Next: Parameters
          </button>
        </div>
      </div>
    );
  };

  // Step 3: Parameters
  const renderParametersStep = () => {
    // Calculate total k-points from the path
    // Each point's npoints indicates k-points in the segment TO the next point
    // Last point has npoints=0, so sum gives total k-points
    const totalKPoints = kPath.reduce((sum, p) => sum + p.npoints, 0);
    // Format path for display
    const pathString = kPath.map((p) => p.label).join(" → ");
    const safeParsePositive = (value: string): number | null => {
      const trimmed = value.trim();
      if (!trimmed) return null;
      const parsed = Number(trimmed);
      if (!Number.isFinite(parsed) || parsed <= 0) return null;
      return parsed;
    };
    const parsedConvThr = safeParsePositive(nscfConvThrInput);
    const parsedMixingBeta = safeParsePositive(nscfMixingBetaInput);
    const parsedDegauss = safeParsePositive(nscfDegaussInput);
    const degaussRequired = nscfOccupations === "smearing";
    const isNbndValid = nbnd === "auto" || (Number.isInteger(nbnd) && nbnd > 0);
    const isConvThrValid = parsedConvThr !== null;
    const isMixingBetaValid = parsedMixingBeta !== null && parsedMixingBeta <= 1.0;
    const isDegaussValid = !degaussRequired || parsedDegauss !== null;
    const hasValidLogPath = !autoSaveLogEnabled || autoSaveLogPath.trim().length > 0;
    const bandsFilband = sanitizeOutputFilename(bandsFilbandInput, "bands.dat");
    const projectionFilproj = sanitizeOutputFilename(projectionFilprojInput, "bands.projwfc.dat");
    const canRun = isNbndValid && isConvThrValid && isMixingBetaValid && isDegaussValid && hasValidLogPath;

    return (
      <div className="wizard-step parameters-step">
        <h3>Band Structure Parameters</h3>
        <p className="step-description">
          Configure NSCF electronic settings and post-processing options for `bands.x` and optional `projwfc.x`.
        </p>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("core")} className="section-header">
            <span className={`collapse-icon ${expandedSections.core ? "expanded" : ""}`}>▶</span>
            Core Band Sampling
            <Tooltip text="Primary controls for the NSCF band path run, including explicit band count selection." />
          </h4>
          {expandedSections.core && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    Number of bands
                    <Tooltip text="QE variable: `nbnd` in `pw.x` NSCF run. Use manual mode to force more conduction bands for crossings/high-energy features." />
                  </label>
                  <div className="nbnd-input">
                    <select
                      value={nbnd === "auto" ? "auto" : "manual"}
                      onChange={(e) => setNbnd(e.target.value === "auto" ? "auto" : 20)}
                    >
                      <option value="auto">Auto (QE default)</option>
                      <option value="manual">Manual</option>
                    </select>
                    {nbnd !== "auto" && (
                      <input
                        type="number"
                        min={1}
                        value={nbnd}
                        onChange={(e) => setNbnd(Math.max(1, parseInt(e.target.value, 10) || 1))}
                      />
                    )}
                  </div>
                  {!isNbndValid && (
                    <span className="param-hint input-error">Use a positive integer for manual `nbnd`.</span>
                  )}
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("nscf")} className="section-header">
            <span className={`collapse-icon ${expandedSections.nscf ? "expanded" : ""}`}>▶</span>
            NSCF Electronic Controls
            <Tooltip text="Advanced `pw.x` NSCF settings: convergence (`conv_thr`), mixing (`mixing_beta`), occupations/smearing, and verbosity." />
          </h4>
          {expandedSections.nscf && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    NSCF convergence threshold
                    <Tooltip text="QE variable: `conv_thr` in `&ELECTRONS`. Tighter thresholds reduce numerical noise in near-degenerate regions at higher cost." />
                  </label>
                  <input
                    type="text"
                    value={nscfConvThrInput}
                    onChange={(e) => setNscfConvThrInput(e.target.value)}
                    placeholder="1e-8"
                    spellCheck={false}
                  />
                  {!isConvThrValid && (
                    <span className="param-hint input-error">Use a positive number (e.g. `1e-8`).</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Mixing beta
                    <Tooltip text="QE variable: `mixing_beta` in `&ELECTRONS`. Typical stable range is 0.2-0.8; values above 1 are generally unstable." />
                  </label>
                  <input
                    type="text"
                    value={nscfMixingBetaInput}
                    onChange={(e) => setNscfMixingBetaInput(e.target.value)}
                    placeholder="0.7"
                    spellCheck={false}
                  />
                  {!isMixingBetaValid && (
                    <span className="param-hint input-error">Use a positive number in the range (0, 1].</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Occupations
                    <Tooltip text="QE variable: `occupations` in `&SYSTEM`. Use smearing for metals and challenging Fermi-level crossings." />
                  </label>
                  <select
                    value={nscfOccupations}
                    onChange={(e) => setNscfOccupations(e.target.value as "fixed" | "smearing" | "from_input" | "tetrahedra")}
                  >
                    <option value="smearing">smearing</option>
                    <option value="fixed">fixed</option>
                    <option value="tetrahedra">tetrahedra</option>
                    <option value="from_input">from_input</option>
                  </select>
                </div>

                <div className="phonon-field">
                  <label>
                    Smearing type
                    <Tooltip text="QE variable: `smearing`. Only used when occupations = smearing." />
                  </label>
                  <select
                    value={nscfSmearing}
                    onChange={(e) => setNscfSmearing(e.target.value as "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac")}
                    disabled={nscfOccupations !== "smearing"}
                  >
                    <option value="gaussian">gaussian</option>
                    <option value="methfessel-paxton">methfessel-paxton</option>
                    <option value="marzari-vanderbilt">marzari-vanderbilt</option>
                    <option value="fermi-dirac">fermi-dirac</option>
                  </select>
                </div>

                <div className="phonon-field">
                  <label>
                    Degauss (Ry)
                    <Tooltip text="QE variable: `degauss`. Required when using smearing occupations; controls broadening width in Ry." />
                  </label>
                  <input
                    type="text"
                    value={nscfDegaussInput}
                    onChange={(e) => setNscfDegaussInput(e.target.value)}
                    placeholder={nscfOccupations === "smearing" ? "0.02" : "unused"}
                    spellCheck={false}
                    disabled={nscfOccupations !== "smearing"}
                  />
                  {!isDegaussValid && (
                    <span className="param-hint input-error">Provide a positive degauss when occupations = smearing.</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Output detail level
                    <Tooltip text="QE variable: `verbosity` in `&CONTROL`. Higher verbosity helps debugging but increases output volume." />
                  </label>
                  <select
                    value={nscfVerbosity}
                    onChange={(e) => setNscfVerbosity(e.target.value as "low" | "high" | "debug")}
                  >
                    <option value="low">low</option>
                    <option value="high">high</option>
                    <option value="debug">debug</option>
                  </select>
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("post")} className="section-header">
            <span className={`collapse-icon ${expandedSections.post ? "expanded" : ""}`}>▶</span>
            bands.x Post-Processing
            <Tooltip text="Configure `bands.x` output file and ordering controls (`filband`, `lsym`, `no_overlap`)." />
          </h4>
          {expandedSections.post && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    bands.x output file
                    <span className="band-control-tech-name">filband</span>
                    <Tooltip text="QE variable: `filband` in `&BANDS`. The parsed file is `<filband>.gnu`." />
                  </label>
                  <input
                    type="text"
                    value={bandsFilbandInput}
                    onChange={(e) => setBandsFilbandInput(e.target.value)}
                    placeholder="bands.dat"
                  />
                  <span className="param-hint">Resolved name: `{bandsFilband}`</span>
                </div>
              </div>
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={bandsLsym}
                  onChange={(e) => setBandsLsym(e.target.checked)}
                />
                <span>
                  Enable symmetry handling in bands.x
                  <span className="band-control-tech-name">lsym</span>
                  <Tooltip text="QE variable: `lsym` in `&BANDS`. Reorders/classifies bands using symmetry information; this can change band index assignment across crossings." />
                </span>
              </label>
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={bandsNoOverlap}
                  onChange={(e) => setBandsNoOverlap(e.target.checked)}
                  disabled={bandsLsym}
                />
                <span>
                  Keep default energy-order indexing
                  <span className="band-control-tech-name">no_overlap</span>
                  <Tooltip text="QE variable: `no_overlap` in `&BANDS`. Only used when `lsym = .false.`. Keep enabled for default ordering at each k-point; disable to order by maximal overlap with neighboring k-points." />
                </span>
              </label>
              {bandsLsym && (
                <span className="param-hint">`no_overlap` is ignored while `lsym` is enabled.</span>
              )}
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("projections")} className="section-header">
            <span className={`collapse-icon ${expandedSections.projections ? "expanded" : ""}`}>▶</span>
            Projection Controls (projwfc.x)
            <Tooltip text="Fat-band projection settings and projwfc.x output controls." />
          </h4>
          {expandedSections.projections && (
            <div className="option-params">
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={enableProjections}
                  onChange={(e) => setEnableProjections(e.target.checked)}
                />
                <span>
                  Enable orbital projections (fat bands)
                  <Tooltip text="Runs `projwfc.x` after bands.x and attaches parsed projection groups to the result." />
                </span>
              </label>

              {enableProjections && (
                <div className="option-params">
                  <div className="phonon-grid">
                    <div className="phonon-field">
                      <label>
                        projwfc output file
                        <span className="band-control-tech-name">filproj</span>
                        <Tooltip text="QE variable: `filproj` in `&PROJWFC`. Parsed from this file when available." />
                      </label>
                      <input
                        type="text"
                        value={projectionFilprojInput}
                        onChange={(e) => setProjectionFilprojInput(e.target.value)}
                        placeholder="bands.projwfc.dat"
                      />
                      <span className="param-hint">Resolved name: `{projectionFilproj}`</span>
                    </div>
                  </div>

                  <label className="option-checkbox">
                    <input
                      type="checkbox"
                      checked={projectionLsym}
                      onChange={(e) => setProjectionLsym(e.target.checked)}
                    />
                    <span>
                      Symmetrize projection weights
                      <span className="band-control-tech-name">lsym</span>
                      <Tooltip text="projwfc.x `lsym`: symmetry-averages projection weights." />
                    </span>
                  </label>

                  <label className="option-checkbox">
                    <input
                      type="checkbox"
                      checked={projectionDiagBasis}
                      onChange={(e) => setProjectionDiagBasis(e.target.checked)}
                    />
                    <span>
                      Local crystal-field basis
                      <span className="band-control-tech-name">diag_basis</span>
                      <Tooltip text="projwfc.x `diag_basis`: rotates local orbital basis into a diagonal crystal-field frame." />
                    </span>
                  </label>

                  <label className="option-checkbox">
                    <input
                      type="checkbox"
                      checked={projectionPawproj}
                      onChange={(e) => setProjectionPawproj(e.target.checked)}
                    />
                    <span>
                      PAW projection correction
                      <span className="band-control-tech-name">pawproj</span>
                      <Tooltip text="projwfc.x `pawproj`: PAW-specific projection treatment." />
                    </span>
                  </label>
                </div>
              )}
            </div>
          )}
        </div>

        <div className="option-section config-section">
          <h4>Run Logging</h4>
          <div className="option-params">
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={autoSaveLogEnabled}
                onChange={(e) => setAutoSaveLogEnabled(e.target.checked)}
              />
              <span>
                Auto-save calculation log to local file
                <Tooltip text="Writes the full run log to the configured local file path after completion (success or failure)." />
              </span>
            </label>

            {autoSaveLogEnabled && (
              <div className="log-path-controls">
                <label>
                  Log file path
                </label>
                <div className="log-path-row">
                  <input
                    className="config-path-input"
                    type="text"
                    value={autoSaveLogPath}
                    onChange={(e) => setAutoSaveLogPath(e.target.value)}
                    placeholder="/tmp/bands-run.log"
                    spellCheck={false}
                  />
                  <button
                    className="secondary-button"
                    onClick={() => void handleChooseAutoLogPath()}
                    disabled={isSavingLog}
                  >
                    Choose...
                  </button>
                </div>
                <span className="param-hint">
                  The file is overwritten each run.
                </span>
                {autoSaveLogPath.trim().length === 0 && (
                  <span className="param-hint input-error">
                    Set a local path before running.
                  </span>
                )}
              </div>
            )}
            {logSaveStatus && <span className="param-hint">{logSaveStatus}</span>}
          </div>
        </div>

        {isHpcMode ? (
          <HpcRunSettings
            profileId={activeHpcProfile?.id ?? null}
            profileName={activeHpcProfile?.name ?? "Andromeda"}
            taskKind="bands"
            commandLines={hpcCommandLines}
            resources={hpcResources}
            resourceMode={activeHpcProfile?.resource_mode ?? "both"}
            defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
            defaultGpuResources={activeHpcProfile?.default_gpu_resources ?? null}
            onResourcesChange={setHpcResources}
            disabled={isRunning}
          />
        ) : (
          <div className="option-section mpi-section config-section collapsible">
            <h4 onClick={() => toggleSection("mpi")} className="section-header">
              <span className={`collapse-icon ${expandedSections.mpi ? "expanded" : ""}`}>▶</span>
              Parallelization
              <Tooltip text="Process-level MPI controls for the NSCF stage." />
            </h4>
            {expandedSections.mpi && (
              <div className="option-params">
                {mpiAvailable ? (
                  <div className="mpi-toggle">
                    <label>
                      <input
                        type="checkbox"
                        checked={mpiEnabled}
                        onChange={(e) => setMpiEnabled(e.target.checked)}
                      />
                      Enable MPI ({cpuCount} cores available)
                    </label>
                    {mpiEnabled && (
                      <div className="mpi-procs">
                        <label>
                          Number of processes:
                          <input
                            type="number"
                            min={1}
                            max={cpuCount}
                            value={mpiProcs}
                            onChange={(e) => setMpiProcs(Math.max(1, parseInt(e.target.value, 10) || 1))}
                          />
                        </label>
                      </div>
                    )}
                  </div>
                ) : (
                  <p className="mpi-unavailable">
                    MPI not available. Running in serial mode.
                  </p>
                )}
              </div>
            )}
          </div>
        )}

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>K-path:</span>
            <span>{pathString || "Not selected"}</span>
          </div>
          <div className="summary-row">
            <span>Total k-points:</span>
            <span>{totalKPoints}</span>
          </div>
          <div className="summary-row">
            <span>Requested `nbnd`:</span>
            <span>{nbnd === "auto" ? "Auto" : nbnd}</span>
          </div>
          <div className="summary-row">
            <span>bands.x output:</span>
            <span>{bandsFilband}.gnu</span>
          </div>
          <div className="summary-row">
            <span>Fat-band projections:</span>
            <span>{enableProjections ? `Enabled (${projectionFilproj})` : "Disabled"}</span>
          </div>
          <div className="summary-row">
            <span>Log auto-save:</span>
            <span>
              {autoSaveLogEnabled
                ? (autoSaveLogPath.trim().length > 0 ? autoSaveLogPath : "Enabled (path not set)")
                : "Disabled"}
            </span>
          </div>
        </div>

        {selectedScfDependencyBlocked && (
          <div className="warning-banner">
            <strong>Remote SCF dependency detected.</strong>{" "}
            This source SCF was run on HPC, and local bands requires a full local SCF bundle.
            <div className="run-error-actions">
              <button
                className="secondary-button"
                onClick={() => void handleSwitchToHpcMode()}
                disabled={isResolvingDependency}
              >
                Switch to HPC Mode
              </button>
              <button
                className="primary-button"
                onClick={() => void handleDownloadDependencyBundle()}
                disabled={isResolvingDependency}
              >
                {isResolvingDependency ? "Downloading..." : "Download Full Bundle"}
              </button>
            </div>
            {dependencyStatus && <div className="param-hint">{dependencyStatus}</div>}
          </div>
        )}

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("kpath")}>
            Back
          </button>
          {!isHpcMode && (
            <button
              className="secondary-button"
              onClick={queueCalculation}
              disabled={!canRun || selectedScfDependencyBlocked || isResolvingDependency}
            >
              Queue Task
            </button>
          )}
          <button
            className="primary-button"
            onClick={runCalculation}
            disabled={!canRun || hasBlockingExternalTask || selectedScfDependencyBlocked || isResolvingDependency}
          >
            {isHpcMode ? "Submit Bands to Andromeda" : "Run Calculation"}
          </button>
        </div>
      </div>
    );
  };

  // Step 4: Running
  const renderRunStep = () => {
    return (
      <div className="wizard-step run-step run-step-focused">
        <div className="run-step-headline">
          <h3>{isRunning ? "Running Band Structure Calculation" : "Band Structure Output"}</h3>
          <span className={`run-step-status-pill ${isRunning ? "running" : error ? "error" : "idle"}`}>
            {isRunning ? "Live output" : error ? "Run failed" : "Output"}
          </span>
        </div>

        <div className="run-status-rail">
          <ProgressBar
            status={progress.status}
            percent={progress.percent}
            phase={progress.phase}
            detail={progress.detail}
            compact
          />
          <div className="run-status-meta">
            <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />
          </div>
        </div>

        {error && (
          <>
            <div className="run-inline-error">{error}</div>
            <div className="run-error-actions">
              <button className="secondary-button" onClick={() => setStep("parameters")}>
                Back to Parameters
              </button>
            </div>
          </>
        )}

        <div className={`run-layout ${isHpcMode ? "run-layout-hpc-telemetry" : "run-layout-single"}`}>
          <LiveOutputPanel
            title={isRunning ? "Running..." : "Output"}
            output={output}
            placeholder="Starting calculation..."
            totalLineCount={outputLineCount}
            visibleLineCount={visibleOutputLineCount}
          />

          {isHpcMode && (
            <div className="telemetry-panel">
              <div className="telemetry-header">
                <h3>Remote Utilization</h3>
                <span className="telemetry-meta">
                  {hpcTelemetryLoading
                    ? "Refreshing..."
                    : hpcTelemetryUpdatedAt
                    ? `Updated ${new Date(hpcTelemetryUpdatedAt).toLocaleTimeString()}`
                    : "Waiting for first sample..."}
                </span>
              </div>
              <div className="telemetry-meta-row">
                <span>Job: {activeTask?.hpc?.remote_job_id || "pending allocation"}</span>
                <span>Node: {activeTask?.hpc?.remote_node || "pending"}</span>
                <span>Source: {hpcTelemetrySource}</span>
              </div>
              <pre className="telemetry-output">{hpcTelemetryOutput}</pre>
              {hpcTelemetryError && <p className="telemetry-error">{hpcTelemetryError}</p>}
            </div>
          )}
        </div>
      </div>
    );
  };

  // Step 5: Results
  const renderResultsStep = () => {
    if (!bandData) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("parameters")}>
            Try Again
          </button>
        </div>
      );
    }

    return (
      <div className="wizard-step results-step">
        <h3>Band Structure Results</h3>
        <p className="step-description">
          Calculation complete. Use the main Bands viewer for full plotting controls.
        </p>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Bands:</span>
              <span className="value">{bandData.n_bands}</span>
            </div>
            <div className="summary-item">
              <span className="label">K-points:</span>
              <span className="value">{bandData.n_kpoints}</span>
            </div>
            <div className="summary-item">
              <span className="label">Fermi Energy:</span>
              <span className="value">{bandData.fermi_energy.toFixed(3)} eV</span>
            </div>
            {/* Band gap info - disabled for now
            {bandData.band_gap ? (
              <>
                <div className="summary-item">
                  <span className="label">Band Gap:</span>
                  <span className="value gap-value">{bandData.band_gap.value.toFixed(3)} eV</span>
                </div>
                <div className="summary-item">
                  <span className="label">Gap Type:</span>
                  <span className="value">{bandData.band_gap.is_direct ? "Direct" : "Indirect"}</span>
                </div>
              </>
            ) : (
              <div className="summary-item">
                <span className="label">Character:</span>
                <span className="value metal-value">Metallic</span>
              </div>
            )}
            */}
          </div>
        </div>

        {/* Collapsible calculation output */}
        <div className="output-section">
          <button
            className="output-toggle"
            onClick={() => setShowOutput(!showOutput)}
          >
            {showOutput ? "▼" : "▶"} Calculation Output
          </button>
          {showOutput && (
            <pre className="calculation-output results-output">
              {output || "No output available"}
            </pre>
          )}
        </div>

        {isSaved && (
          <div className="save-status save-status-inline">
            <span className="saved">Saved to project</span>
          </div>
        )}

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Back to Dashboard
          </button>
          <button
            className="primary-button"
            onClick={() => onViewBands(bandData, scfFermiEnergy)}
          >
            View Bands
          </button>
          <button className="primary-button" onClick={() => setStep("parameters")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className={`band-structure-wizard wizard-step-${step}`}>
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Band Structure Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "kpath" ? "active" : ["parameters", "run", "results"].includes(step) ? "completed" : ""}>
            2. K-Path
          </span>
          <span className={step === "parameters" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            3. Parameters
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "completed" : ""}>
            4. Run
          </span>
          <span className={step === "results" ? "active" : ""}>
            5. Results
          </span>
        </div>
      </div>

      <div className="wizard-content">
        {renderStep()}
      </div>
    </div>
  );
}
