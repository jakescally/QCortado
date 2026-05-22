import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import {
  CrystalData,
  ELEMENT_MASSES,
  ExecutionMode,
  HpcProfile,
  SlurmResourceRequest,
} from "../lib/types";
import {
  BandData,
  BandPlot,
  BandPlotWindowOverlay,
  BandPlotWindowOverlayUpdate,
  buildBandPlotProjectionOptions,
  getDefaultBandPlotEnergyRange,
} from "./BandPlot";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";
import {
  analyzeCrystalSymmetry,
  buildConventionalLatticeFromCrystalData,
  SymmetryTransformResult,
} from "../lib/symmetryTransform";
import { inferQeBravaisCellFromCif } from "../lib/engines/qe/bravaisInference";
import {
  createPathCoordinateConverters,
  mapPathCoordinates,
  resolvePathTransformContext,
  roundVec3,
  sourceScfUsesPrimitiveCell,
} from "../lib/kPathTransforms";
import {
  RhombohedralConvention,
  defaultRhombohedralConventionForSetting,
} from "../lib/brillouinZoneData";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/engines/qe/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { InfoTooltip } from "./InfoTooltip";
import { useTaskContext } from "../lib/TaskContext";
import { defaultProgressState, ProgressState } from "../lib/engines/qe/progress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import {
  buildExecutionTarget,
  defaultResourcesForProfile,
} from "../lib/hpcConfig";
import {
  buildHpcQeInputCommandLine,
  resolveProfileRemotePseudoDir,
} from "../lib/engines/qe/hpc";
import { validateHpcTasksWithinBandCount } from "../lib/hpcBandLimits";
import { HpcRunSettings } from "./HpcRunSettings";
import {
  formatWannierConvergenceFlag,
  getWannierQualityIssues,
  WannierQualityIssue,
} from "../lib/engines/qe/wannierQuality";
import {
  getNeutralElectronConfiguration,
  getOutermostOccupiedOrbital,
} from "../lib/electronConfigurations";
import { resolveSavedScfStructure } from "../lib/optimizedStructure";
import { getMagneticSpeciesFields } from "../lib/magnetism";
import { formatCalculationSourceLabel, getCalculationName } from "../lib/calculationNames";
import { readProjectWizardSettings, writeProjectWizardSettings } from "../lib/projectWizardSettings";
import type { EngineId } from "../lib/engines/types";

interface CalculationRun {
  id: string;
  engine_id?: EngineId | null;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    fermi_energy: number | null;
    band_data?: BandData | null;
    raw_output?: string | null;
    wannier_data?: WannierResult | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

interface WannierBandOverlayOption {
  id: string;
  label: string;
  data: BandData;
  fermiEnergy: number | null;
}

interface WannierSpread {
  index: number;
  centre: [number, number, number];
  spread: number;
}

interface WannierCentreRecord {
  label: string;
  position: [number, number, number];
}

interface WannierConvergenceData {
  converged: boolean;
  iterations?: number | null;
  minimization_converged?: boolean | null;
  disentanglement_converged?: boolean | null;
  max_iterations_reached?: boolean | null;
  omega_i?: number | null;
  omega_d?: number | null;
  omega_od?: number | null;
  omega_total?: number | null;
  warnings?: string[];
  failure_reasons?: string[];
}

interface WannierFermiAlignment {
  source_brackets_fermi?: boolean | null;
  wannier_brackets_fermi: boolean;
  source_min_distance_ev?: number | null;
  wannier_min_distance_ev: number;
  source_energy_range_ev?: [number, number] | null;
  wannier_energy_range_ev?: [number, number] | null;
}

interface WannierArtifact {
  file_name: string;
  size_bytes: number;
}

interface WannierResult {
  seedname: string;
  num_wann: number;
  num_bands: number;
  k_grid: [number, number, number];
  band_data: BandData;
  total_spread?: number | null;
  spreads: WannierSpread[];
  centres: WannierCentreRecord[];
  convergence: WannierConvergenceData;
  fermi_alignment?: WannierFermiAlignment | null;
  quality_issues?: WannierQualityIssue[];
  artifact_manifest: WannierArtifact[];
}

type ProjectionTargetType = "element" | "site";
type ProjectionOrbital = "s" | "p" | "d" | "f" | "sp" | "sp2" | "sp3" | "sp3d" | "sp3d2";
type WizardStep = "source" | "mesh" | "projections" | "run" | "results";
type NumBandsMode = "auto" | "manual";
type ReferenceEnergyMode = "relative" | "absolute";
const WANNIER_WIZARD_SETTINGS_ID = "wannier";

interface ProjectionDraft {
  id: string;
  targetType: ProjectionTargetType;
  symbol: string;
  orbital: ProjectionOrbital;
  siteIndex: number | null;
}

interface StoredWannierWizardSettings {
  kGridInput: [string, string, string];
  totalKPointsTarget: number;
  totalKPointsInput: string;
  numBandsInput: string;
  numBandsMode: NumBandsMode;
  seednameInput: string;
  projectionDrafts: ProjectionDraft[];
  numWannOverrideInput: string;
  excludeBandsInput: string;
  disWinMinInput: string;
  disWinMaxInput: string;
  disFrozMinInput: string;
  disFrozMaxInput: string;
  referenceEnergyMode: ReferenceEnergyMode;
  referenceShowBandGapOverlay: boolean;
  referenceProjectionSelection: string;
}

interface WannierWizardProps {
  onBack: () => void;
  onExecutionModeChange?: (mode: ExecutionMode) => Promise<void> | void;
  onViewWannier: (
    result: WannierResult,
    scfFermiEnergy: number | null,
    overlayOptions?: WannierBandOverlayOption[],
  ) => void;
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

const WANNIER_WORK_DIR = "/tmp/qcortado_wannier";
const ORBITAL_OPTIONS: ProjectionOrbital[] = ["s", "p", "d", "f", "sp", "sp2", "sp3", "sp3d", "sp3d2"];

function getBaseElement(symbol: string): string {
  return symbol.replace(/[\d+-]+$/, "");
}

function normalizeOccupations(raw: unknown): "fixed" | "smearing" | "from_input" | "tetrahedra" | "tetrahedra_lin" | "tetrahedra_opt" {
  const lowered = String(raw || "fixed").toLowerCase();
  if (lowered === "smearing") return "smearing";
  if (lowered === "from_input") return "from_input";
  if (lowered === "tetrahedra") return "tetrahedra";
  if (lowered === "tetrahedra_lin") return "tetrahedra_lin";
  if (lowered === "tetrahedra_opt") return "tetrahedra_opt";
  return "fixed";
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

function orbitalTemplateCount(orbital: ProjectionOrbital): number {
  switch (orbital) {
    case "s":
      return 1;
    case "p":
      return 3;
    case "d":
      return 5;
    case "f":
      return 7;
    case "sp":
      return 4;
    case "sp2":
      return 3;
    case "sp3":
      return 4;
    case "sp3d":
      return 5;
    case "sp3d2":
      return 6;
    default:
      return 1;
  }
}

function renderElectronConfiguration(config: string) {
  const tokens = config.split(/\s+/).filter(Boolean);
  return (
    <>
      {tokens.map((token, index) => {
        const subshellMatch = token.match(/^(\d)([spdf])(\d+)$/);
        const spacer = index < tokens.length - 1 ? " " : "";
        if (!subshellMatch) {
          return <span key={`${token}-${index}`}>{token}{spacer}</span>;
        }
        return (
          <span key={`${token}-${index}`}>
            {subshellMatch[1]}
            {subshellMatch[2]}
            <sup>{subshellMatch[3]}</sup>
            {spacer}
          </span>
        );
      })}
    </>
  );
}

function getNeutralElectronConfigSummary(symbol: string) {
  return getNeutralElectronConfiguration(symbol);
}

function getDefaultProjectionOrbital(symbol: string): ProjectionOrbital {
  return getOutermostOccupiedOrbital(symbol) ?? "s";
}

function countMatchingElementSites(symbol: string, atomSymbols: string[]): number {
  const trimmed = getBaseElement(symbol).trim().toLowerCase();
  if (!trimmed) return 0;
  return atomSymbols.filter((entry) => getBaseElement(entry).trim().toLowerCase() === trimmed).length;
}

function projectionContributionCount(
  draft: ProjectionDraft,
  atomSymbols: string[],
): number {
  const orbitalCount = orbitalTemplateCount(draft.orbital);
  if (draft.targetType === "site") {
    return orbitalCount;
  }
  return orbitalCount * countMatchingElementSites(draft.symbol, atomSymbols);
}

function parseElectronCount(rawOutput: string | null | undefined): number | null {
  if (!rawOutput) return null;
  const match = rawOutput.match(/number of electrons\s*=\s*([0-9.+-Ee]+)/i);
  if (!match) return null;
  const value = Number(match[1]);
  return Number.isFinite(value) && value > 0 ? value : null;
}

function parseOptionalNumber(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (!trimmed) return null;
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed)) {
    throw new Error(`Invalid ${label}.`);
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

function sanitizeSeedname(raw: string): string {
  const trimmed = raw.trim();
  if (!trimmed) return "qcortado_wannier";
  const sanitized = trimmed
    .replace(/[^a-zA-Z0-9._-]/g, "_")
    .replace(/^_+|_+$/g, "");
  return sanitized || "qcortado_wannier";
}

function parseExcludeBandsInput(raw: string): number[] {
  const values: number[] = [];
  for (const chunk of raw.split(",")) {
    const token = chunk.trim();
    if (!token) continue;
    const range = token.split("-").map((part) => Number.parseInt(part.trim(), 10));
    if (range.length === 2 && Number.isFinite(range[0]) && Number.isFinite(range[1])) {
      const start = Math.min(range[0], range[1]);
      const end = Math.max(range[0], range[1]);
      for (let value = start; value <= end; value += 1) {
        if (value > 0) values.push(value);
      }
      continue;
    }
    const single = Number.parseInt(token, 10);
    if (Number.isFinite(single) && single > 0) {
      values.push(single);
    }
  }
  return [...new Set(values)].sort((a, b) => a - b);
}

function parseFiniteInput(input: string): number | null {
  const trimmed = input.trim();
  if (!trimmed) return null;
  const parsed = Number.parseFloat(trimmed);
  return Number.isFinite(parsed) ? parsed : null;
}

function formatEnergyInputValue(value: number): string {
  if (!Number.isFinite(value)) return "";
  return Number.parseFloat(value.toFixed(6)).toString();
}

function toAscendingPair(first: number, second: number): [number, number] {
  return first <= second ? [first, second] : [second, first];
}

function clampNumber(value: number, min: number, max: number): number {
  return Math.min(max, Math.max(min, value));
}

const MAX_VIEWER_POINTS_PER_SEGMENT = 400;
const MAX_TOTAL_K_POINTS = 5000;

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

function applyTotalKPoints(path: KPathPoint[], totalKPoints: number): KPathPoint[] {
  const connectedSegmentIndices = getConnectedSegmentIndices(path);
  if (path.length === 0 || connectedSegmentIndices.length === 0) {
    return path.map((point) => ({
      ...point,
      npoints: 0,
    }));
  }

  const safeTotal = clampInt(totalKPoints, connectedSegmentIndices.length, MAX_TOTAL_K_POINTS);
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

function getWannierSourceIssues(calc: CalculationRun): string[] {
  const params = calc.parameters || {};
  const issues: string[] = [];
  const nspin = Number(params.nspin ?? 1);
  if (nspin !== 1) {
    issues.push("Scalar nspin=1 required");
  }
  if (Boolean(params.noncolin) || Boolean(params.lspinorb)) {
    issues.push("No noncollinear / SOC");
  }
  if (Boolean(params.lda_plus_u)) {
    issues.push("DFT+U not supported in v1");
  }
  const vdw = String(params.vdw_corr || "").trim().toLowerCase();
  if (vdw && vdw !== "none") {
    issues.push("vdW source not certified in v1");
  }
  if (normalizeOccupations(params.occupations) === "from_input") {
    issues.push("occupations='from_input' not supported");
  }
  return issues;
}

function getCalculationTags(calc: CalculationRun): { label: string; type: "info" | "feature" | "special" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" }[] = [];
  const params = calc.parameters || {};
  const issues = getWannierSourceIssues(calc);
  const pushTag = (label: string, type: "info" | "feature" | "special") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }
  if (params.conv_thr) {
    const thr = params.conv_thr;
    pushTag(thr < 0.001 ? thr.toExponential(0) : thr.toString(), "info");
  }
  if (issues.length === 0) {
    pushTag("Wannier-Ready", "special");
  }
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
  if (String(params.execution_backend || "").trim().toLowerCase() === "hpc") {
    pushTag("HPC", "feature");
  }

  return tags;
}

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes < 1024) return `${bytes || 0} B`;
  const kb = bytes / 1024;
  if (kb < 1024) return `${kb.toFixed(1)} KB`;
  const mb = kb / 1024;
  if (mb < 1024) return `${mb.toFixed(1)} MB`;
  return `${(mb / 1024).toFixed(1)} GB`;
}

function makeProjectionId(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") {
    return crypto.randomUUID();
  }
  return `proj_${Date.now()}_${Math.floor(Math.random() * 1_000_000)}`;
}

function formatReferenceBandsLabel(calc: CalculationRun): string {
  const dateLabel = new Date(calc.completed_at ?? calc.started_at).toLocaleDateString();
  const sourceScfId = String(calc.parameters?.source_scf_id || "").trim();
  const sourceLabel = sourceScfId ? sourceScfId.slice(0, 8) : "unknown";
  const hasProjections = Boolean(calc.result?.band_data?.projections);
  return `${dateLabel} · ${calc.id.slice(0, 8)} · source ${sourceLabel} · ${hasProjections ? "fat bands" : "line bands"}`;
}

function normalizeSavedKPath(value: unknown): string {
  return String(value || "")
    .replace(/\s*→\s*/g, "→")
    .replace(/\s+/g, " ")
    .trim();
}

export function WannierWizard({
  onBack,
  onViewWannier,
  qePath,
  defaultSmearing = "marzari-vanderbilt",
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: WannierWizardProps) {
  const resolvedDefaultSmearing = normalizeSmearing(defaultSmearing);
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const storedWizardSettings = useMemo(
    () => readProjectWizardSettings<StoredWannierWizardSettings>(projectId, WANNIER_WIZARD_SETTINGS_ID),
    [projectId],
  );
  const shouldPreserveStoredSettingsRef = useRef(Boolean(storedWizardSettings));
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  const hasBlockingExternalTask = !isHpcMode && hasExternalRunningTask;
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());
  const [symmetryTransform, setSymmetryTransform] = useState<SymmetryTransformResult | null>(null);
  const [symmetryError, setSymmetryError] = useState<string | null>(null);
  const [kGridInput, setKGridInput] = useState<[string, string, string]>(() => storedWizardSettings?.kGridInput ?? ["4", "4", "4"]);
  const [totalKPointsTarget, setTotalKPointsTarget] = useState(() => storedWizardSettings?.totalKPointsTarget ?? 120);
  const [totalKPointsInput, setTotalKPointsInput] = useState(() => storedWizardSettings?.totalKPointsInput ?? "120");
  const [numBandsInput, setNumBandsInput] = useState(() => storedWizardSettings?.numBandsInput ?? "");
  const [numBandsMode, setNumBandsMode] = useState<NumBandsMode>(() => storedWizardSettings?.numBandsMode ?? "auto");
  const [seednameInput, setSeednameInput] = useState(() => storedWizardSettings?.seednameInput ?? "qcortado_wannier");
  const [kPath, setKPath] = useState<KPathPoint[]>([]);
  const [kPathRhombohedralConvention, setKPathRhombohedralConvention] = useState<RhombohedralConvention | undefined>(undefined);
  const [projectionDrafts, setProjectionDrafts] = useState<ProjectionDraft[]>(() => storedWizardSettings?.projectionDrafts ?? []);
  const [numWannOverrideInput, setNumWannOverrideInput] = useState(() => storedWizardSettings?.numWannOverrideInput ?? "");
  const [excludeBandsInput, setExcludeBandsInput] = useState(() => storedWizardSettings?.excludeBandsInput ?? "");
  const [disWinMinInput, setDisWinMinInput] = useState(() => storedWizardSettings?.disWinMinInput ?? "");
  const [disWinMaxInput, setDisWinMaxInput] = useState(() => storedWizardSettings?.disWinMaxInput ?? "");
  const [disFrozMinInput, setDisFrozMinInput] = useState(() => storedWizardSettings?.disFrozMinInput ?? "");
  const [disFrozMaxInput, setDisFrozMaxInput] = useState(() => storedWizardSettings?.disFrozMaxInput ?? "");
  const [referenceBandsCalcId, setReferenceBandsCalcId] = useState("");
  const [referenceEnergyMode, setReferenceEnergyMode] = useState<ReferenceEnergyMode>(() => storedWizardSettings?.referenceEnergyMode ?? "absolute");
  const [referenceShowBandGapOverlay, setReferenceShowBandGapOverlay] = useState(() => storedWizardSettings?.referenceShowBandGapOverlay ?? true);
  const [referenceProjectionSelection, setReferenceProjectionSelection] = useState(() => storedWizardSettings?.referenceProjectionSelection ?? "none");
  const [detailedBandsById, setDetailedBandsById] = useState<Record<string, CalculationRun>>({});
  const [referenceBandsLoading, setReferenceBandsLoading] = useState(false);
  const [progress, setProgress] = useState<ProgressState>(defaultProgressState("Wannier90"));
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const [outputLineCount, setOutputLineCount] = useState(0);
  const [result, setResult] = useState<WannierResult | null>(null);
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [isSaved, setIsSaved] = useState(false);
  const [showOutput, setShowOutput] = useState(true);
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(
    defaultResourcesForProfile(activeHpcProfile),
  );
  const visibleOutputLineCount = useMemo(() => countVisibleOutputLines(output), [output]);
  useViewportScrollLock(step === "run");
  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);
  const sourcePseudoMap = useMemo<Record<string, string>>(() => {
    const params = selectedScf?.parameters || {};
    return (params.selected_pseudos && typeof params.selected_pseudos === "object")
      ? params.selected_pseudos as Record<string, string>
      : {};
  }, [selectedScf]);
  const uniqueCrystalElements = useMemo(
    () => [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))],
    [crystalData],
  );
  const conventionalAtomCount = crystalData.atom_sites.length;
  const structureElectronConfigs = useMemo(
    () => uniqueCrystalElements.map((element) => ({
      element,
      config: getNeutralElectronConfigSummary(element),
    })),
    [uniqueCrystalElements],
  );

  function getAllowedOrbitalOptionsForDraft(draft: ProjectionDraft): ProjectionOrbital[] {
    void draft;
    return ORBITAL_OPTIONS;
  }

  function getProjectionElement(draft: ProjectionDraft): string {
    return draft.targetType === "site"
      ? getBaseElement(crystalData.atom_sites[draft.siteIndex ?? 0]?.type_symbol || draft.symbol)
      : getBaseElement(draft.symbol);
  }

  async function prepareSelectedScfSource(scf: CalculationRun): Promise<void> {
    const params = scf.parameters || {};
    if (sourceScfUsesPrimitiveCell(params, conventionalAtomCount)) {
      await ensureSymmetryTransform();
    }
  }

  useEffect(() => {
    if (visibleOutputLineCount > outputLineCount) {
      setOutputLineCount(visibleOutputLineCount);
    }
  }, [outputLineCount, visibleOutputLineCount]);

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [isHpcMode, activeHpcProfile?.id, activeHpcProfile?.resource_mode]);

  useEffect(() => {
    writeProjectWizardSettings(projectId, WANNIER_WIZARD_SETTINGS_ID, {
      kGridInput,
      totalKPointsTarget,
      totalKPointsInput,
      numBandsInput,
      numBandsMode,
      seednameInput,
      projectionDrafts,
      numWannOverrideInput,
      excludeBandsInput,
      disWinMinInput,
      disWinMaxInput,
      disFrozMinInput,
      disFrozMaxInput,
      referenceEnergyMode,
      referenceShowBandGapOverlay,
      referenceProjectionSelection,
    });
  }, [
    projectId,
    kGridInput,
    totalKPointsTarget,
    totalKPointsInput,
    numBandsInput,
    numBandsMode,
    seednameInput,
    projectionDrafts,
    numWannOverrideInput,
    excludeBandsInput,
    disWinMinInput,
    disWinMaxInput,
    disFrozMinInput,
    disFrozMaxInput,
    referenceEnergyMode,
    referenceShowBandGapOverlay,
    referenceProjectionSelection,
  ]);

  useEffect(() => {
    async function initLocalExecutionDefaults() {
      if (isHpcMode) return;
      try {
        const count = await invoke<number>("get_cpu_count");
        const safeCount = Math.max(1, Math.floor(count));
        setCpuCount(safeCount);
        const defaults = await loadGlobalMpiDefaults(safeCount);
        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);
        setMpiEnabled(available ? defaults.enabled : false);
        setMpiProcs(defaults.nprocs);
      } catch {
        setCpuCount(1);
        setMpiAvailable(false);
        setMpiEnabled(false);
        setMpiProcs(1);
      }
    }

    void initLocalExecutionDefaults();
  }, [isHpcMode]);

  useEffect(() => {
    const task = activeTask;
    if (!task) return;
    if (task.outputLineCount > 0 || task.status !== "running") {
      setOutput(task.outputText);
      setOutputLineCount(task.outputLineCount);
    }
    setProgress(task.progress);
    if (task.status === "completed" && task.result) {
      const nextResult = task.result as WannierResult;
      setResult(nextResult);
      setIsRunning(false);
      setStep("results");
      setProgress({
        status: "complete",
        percent: 100,
        phase: "Complete",
      });
    } else if (task.status === "failed") {
      setIsRunning(false);
      setError(task.error || "Wannier calculation failed.");
    } else if (task.status === "cancelled") {
      setIsRunning(false);
      setError("Wannier calculation cancelled.");
    }
  }, [activeTask]);

  useEffect(() => {
    if (!selectedScf) return;
    const params = selectedScf.parameters || {};
    if (!shouldPreserveStoredSettingsRef.current) {
      const storedGrid = Array.isArray(params.kgrid) && params.kgrid.length === 3
        ? params.kgrid.map((value: unknown) => String(value || "4")) as [string, string, string]
        : ["4", "4", "4"];
      setKGridInput([storedGrid[0], storedGrid[1], storedGrid[2]]);
    }

    if (projectionDrafts.length === 0) {
      const uniqueElements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
      if (uniqueElements.length === 1) {
        const defaultOrbital = getDefaultProjectionOrbital(uniqueElements[0]);
        setProjectionDrafts([
          {
            id: makeProjectionId(),
            targetType: "element",
            symbol: uniqueElements[0],
            orbital: defaultOrbital,
            siteIndex: null,
          },
        ]);
      }
    }
  }, [selectedScf, crystalData, projectionDrafts.length]);

  useEffect(() => {
    if (shouldPreserveStoredSettingsRef.current) return;
    setNumBandsMode("auto");
    setNumBandsInput("");
  }, [selectedScf?.id]);

  const wannierCellAtomSymbols = useMemo(() => {
    const params = selectedScf?.parameters || {};
    if (
      sourceScfUsesPrimitiveCell(params, conventionalAtomCount)
      && symmetryTransform
      && symmetryTransform.standardizedPrimitiveAtoms.length > 0
    ) {
      return symmetryTransform.standardizedPrimitiveAtoms.map((atom) => getBaseElement(atom.symbol));
    }
    return crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol));
  }, [conventionalAtomCount, crystalData, selectedScf, symmetryTransform]);

  const validScfs = useMemo(
    () => scfCalculations.filter((calc) => calc.calc_type === "scf" && calc.result?.converged),
    [scfCalculations],
  );
  const readyScfs = useMemo(
    () => validScfs.filter((calc) => getWannierSourceIssues(calc).length === 0),
    [validScfs],
  );
  const sortedReadyScfs = useMemo(
    () => sortScfByMode(readyScfs, scfSortMode),
    [readyScfs, scfSortMode],
  );

  const selectedIssues = useMemo(
    () => (selectedScf ? getWannierSourceIssues(selectedScf) : []),
    [selectedScf],
  );
  const bandCalculationSummaries = useMemo(
    () => scfCalculations.filter((calc) => calc.calc_type === "bands"),
    [scfCalculations],
  );
  const hydratedBandCalculations = useMemo(
    () => bandCalculationSummaries.map((calc) => detailedBandsById[calc.id] ?? calc),
    [bandCalculationSummaries, detailedBandsById],
  );
  const referenceBandsCalculations = useMemo(
    () => hydratedBandCalculations
      .filter((calc) => Boolean(calc.result?.band_data))
      .slice()
      .sort((left, right) => {
        const leftTime = new Date(left.completed_at ?? left.started_at).getTime();
        const rightTime = new Date(right.completed_at ?? right.started_at).getTime();
        return rightTime - leftTime;
      }),
    [hydratedBandCalculations],
  );
  const selectedReferenceBandsCalculation = useMemo(
    () => referenceBandsCalculations.find((calc) => calc.id === referenceBandsCalcId)
      ?? referenceBandsCalculations[0]
      ?? null,
    [referenceBandsCalcId, referenceBandsCalculations],
  );
  const selectedReferenceBandData = selectedReferenceBandsCalculation?.result?.band_data ?? null;
  const selectedReferenceFermiEnergy = useMemo(() => {
    const sourceScfEf = Number(selectedScf?.result?.fermi_energy);
    if (Number.isFinite(sourceScfEf)) {
      return sourceScfEf;
    }
    const sourceBandsEf = Number(selectedReferenceBandsCalculation?.result?.fermi_energy);
    if (Number.isFinite(sourceBandsEf)) {
      return sourceBandsEf;
    }
    return 0;
  }, [selectedReferenceBandsCalculation?.result?.fermi_energy, selectedScf?.result?.fermi_energy]);
  const referencePlotFermi = referenceEnergyMode === "relative" ? selectedReferenceFermiEnergy : 0;
  const referenceProjectionOptions = useMemo(
    () => selectedReferenceBandData
      ? buildBandPlotProjectionOptions(selectedReferenceBandData)
      : [{ value: "none", label: "none" }],
    [selectedReferenceBandData],
  );

  useEffect(() => {
    if (!selectedScf) return;
    if (readyScfs.some((calc) => calc.id === selectedScf.id)) {
      return;
    }
    setSelectedScf(null);
  }, [readyScfs, selectedScf]);

  useEffect(() => {
    if (step !== "projections") return;
    const missingBandDetails = bandCalculationSummaries.filter((calc) => {
      if (calc.result?.band_data) return false;
      if (detailedBandsById[calc.id]?.result?.band_data) return false;
      return true;
    });
    if (missingBandDetails.length === 0) {
      setReferenceBandsLoading(false);
      return;
    }

    let cancelled = false;
    setReferenceBandsLoading(true);
    void (async () => {
      const details = await Promise.all(
        missingBandDetails.map(async (calc) => {
          try {
            return await invoke<CalculationRun>("get_project_calculation", {
              projectId,
              calcId: calc.id,
            });
          } catch (err) {
            console.warn(`Failed to hydrate bands calculation ${calc.id}:`, err);
            return null;
          }
        }),
      );
      if (cancelled) return;
      setDetailedBandsById((prev) => {
        const next = { ...prev };
        let changed = false;
        for (const detail of details) {
          if (!detail || !detail.result?.band_data) continue;
          if (!next[detail.id]?.result?.band_data) {
            next[detail.id] = detail;
            changed = true;
          }
        }
        return changed ? next : prev;
      });
      setReferenceBandsLoading(false);
    })();

    return () => {
      cancelled = true;
    };
  }, [bandCalculationSummaries, detailedBandsById, projectId, step]);

  useEffect(() => {
    if (!selectedReferenceBandsCalculation) {
      if (referenceBandsCalcId) {
        setReferenceBandsCalcId("");
      }
      return;
    }
    if (selectedReferenceBandsCalculation.id !== referenceBandsCalcId) {
      setReferenceBandsCalcId(selectedReferenceBandsCalculation.id);
    }
  }, [referenceBandsCalcId, selectedReferenceBandsCalculation]);

  useEffect(() => {
    if (!referenceProjectionOptions.some((option) => option.value === referenceProjectionSelection)) {
      setReferenceProjectionSelection("none");
    }
  }, [referenceProjectionOptions, referenceProjectionSelection]);

  const referenceDefaultDisplayRange = useMemo<[number, number] | null>(() => {
    if (!selectedReferenceBandData) {
      return null;
    }
    return getDefaultBandPlotEnergyRange(
      selectedReferenceBandData,
      referencePlotFermi,
      "scf",
    );
  }, [referencePlotFermi, selectedReferenceBandData]);

  const displayEnergyToAbsolute = useCallback((value: number) => {
    if (!Number.isFinite(value)) {
      return value;
    }
    return referenceEnergyMode === "absolute"
      ? value
      : value + selectedReferenceFermiEnergy;
  }, [referenceEnergyMode, selectedReferenceFermiEnergy]);

  const absoluteEnergyToDisplay = useCallback((value: number) => {
    if (!Number.isFinite(value)) {
      return value;
    }
    return referenceEnergyMode === "absolute"
      ? value
      : value - selectedReferenceFermiEnergy;
  }, [referenceEnergyMode, selectedReferenceFermiEnergy]);

  const defaultDisWindowAbsolute = useMemo<[number, number]>(() => {
    if (!referenceDefaultDisplayRange) {
      return [selectedReferenceFermiEnergy - 4, selectedReferenceFermiEnergy + 4];
    }
    const converted = toAscendingPair(
      displayEnergyToAbsolute(referenceDefaultDisplayRange[0]),
      displayEnergyToAbsolute(referenceDefaultDisplayRange[1]),
    );
    if (converted[1] - converted[0] < 1e-4) {
      return [converted[0] - 1, converted[1] + 1];
    }
    return converted;
  }, [displayEnergyToAbsolute, referenceDefaultDisplayRange, selectedReferenceFermiEnergy]);

  const defaultFrozenWindowAbsolute = useMemo<[number, number]>(() => {
    const [disMin, disMax] = defaultDisWindowAbsolute;
    const disSpan = Math.max(disMax - disMin, 1e-4);
    const halfSpan = Math.min(1, disSpan * 0.2);
    let min = selectedReferenceFermiEnergy - halfSpan;
    let max = selectedReferenceFermiEnergy + halfSpan;
    min = clampNumber(min, disMin, disMax);
    max = clampNumber(max, disMin, disMax);
    if (max <= min) {
      min = disMin + disSpan * 0.25;
      max = disMax - disSpan * 0.25;
      if (max <= min) {
        min = disMin;
        max = disMax;
      }
    }
    return [min, max];
  }, [defaultDisWindowAbsolute, selectedReferenceFermiEnergy]);

  const clampFrozenInsideDisWindow = useCallback(
    (frozenCandidate: [number, number], disWindow: [number, number]): [number, number] => {
      const [disMin, disMax] = toAscendingPair(disWindow[0], disWindow[1]);
      const minSpan = Math.max(1e-4, (disMax - disMin) * 1e-6);
      let [frozenMin, frozenMax] = toAscendingPair(frozenCandidate[0], frozenCandidate[1]);
      frozenMin = clampNumber(frozenMin, disMin, disMax);
      frozenMax = clampNumber(frozenMax, disMin, disMax);

      if (frozenMax - frozenMin < minSpan) {
        if (disMax - disMin <= minSpan) {
          return [disMin, disMax];
        }
        const center = clampNumber((frozenMin + frozenMax) / 2, disMin, disMax);
        frozenMin = clampNumber(center - minSpan / 2, disMin, disMax);
        frozenMax = clampNumber(center + minSpan / 2, disMin, disMax);
        if (frozenMax - frozenMin < minSpan) {
          frozenMin = disMin;
          frozenMax = disMax;
        }
      }

      return [frozenMin, frozenMax];
    },
    [],
  );

  const clampDisWindowAroundFrozen = useCallback(
    (disCandidate: [number, number], frozenWindow: [number, number]): [number, number] => {
      const minSpan = 1e-4;
      let [disMin, disMax] = toAscendingPair(disCandidate[0], disCandidate[1]);
      if (disMax - disMin < minSpan) {
        disMax = disMin + minSpan;
      }
      const [frozenMin, frozenMax] = toAscendingPair(frozenWindow[0], frozenWindow[1]);
      if (disMin > frozenMin) {
        disMin = frozenMin;
      }
      if (disMax < frozenMax) {
        disMax = frozenMax;
      }
      if (disMax - disMin < minSpan) {
        disMax = disMin + minSpan;
      }
      return [disMin, disMax];
    },
    [],
  );

  const disWinMinParsed = parseFiniteInput(disWinMinInput);
  const disWinMaxParsed = parseFiniteInput(disWinMaxInput);
  const disFrozMinParsed = parseFiniteInput(disFrozMinInput);
  const disFrozMaxParsed = parseFiniteInput(disFrozMaxInput);

  const effectiveFrozenWindowAbsoluteBase = useMemo<[number, number]>(() => {
    return toAscendingPair(
      disFrozMinParsed ?? defaultFrozenWindowAbsolute[0],
      disFrozMaxParsed ?? defaultFrozenWindowAbsolute[1],
    );
  }, [defaultFrozenWindowAbsolute, disFrozMaxParsed, disFrozMinParsed]);

  const effectiveDisWindowAbsolute = useMemo<[number, number]>(() => {
    const candidate: [number, number] = [
      disWinMinParsed ?? defaultDisWindowAbsolute[0],
      disWinMaxParsed ?? defaultDisWindowAbsolute[1],
    ];
    return clampDisWindowAroundFrozen(candidate, effectiveFrozenWindowAbsoluteBase);
  }, [
    clampDisWindowAroundFrozen,
    defaultDisWindowAbsolute,
    disWinMaxParsed,
    disWinMinParsed,
    effectiveFrozenWindowAbsoluteBase,
  ]);

  const effectiveFrozenWindowAbsolute = useMemo<[number, number]>(() => {
    return clampFrozenInsideDisWindow(
      effectiveFrozenWindowAbsoluteBase,
      effectiveDisWindowAbsolute,
    );
  }, [
    clampFrozenInsideDisWindow,
    effectiveDisWindowAbsolute,
    effectiveFrozenWindowAbsoluteBase,
  ]);

  const setDisentanglementInputsFromAbsolute = useCallback(
    (disWindow: [number, number], frozenWindow: [number, number]) => {
      setDisWinMinInput(formatEnergyInputValue(disWindow[0]));
      setDisWinMaxInput(formatEnergyInputValue(disWindow[1]));
      setDisFrozMinInput(formatEnergyInputValue(frozenWindow[0]));
      setDisFrozMaxInput(formatEnergyInputValue(frozenWindow[1]));
    },
    [],
  );

  useEffect(() => {
    if (!selectedReferenceBandData) {
      setReferenceProjectionSelection("none");
      return;
    }
    setReferenceProjectionSelection("none");
  }, [selectedReferenceBandData, selectedReferenceBandsCalculation?.id]);

  useEffect(() => {
    if (!selectedReferenceBandData) {
      return;
    }
    const hasAnyDisWindowInput = [
      disWinMinInput,
      disWinMaxInput,
      disFrozMinInput,
      disFrozMaxInput,
    ].some((value) => value.trim().length > 0);
    if (hasAnyDisWindowInput) {
      return;
    }
    setDisentanglementInputsFromAbsolute(defaultDisWindowAbsolute, defaultFrozenWindowAbsolute);
  }, [
    defaultDisWindowAbsolute,
    defaultFrozenWindowAbsolute,
    disFrozMaxInput,
    disFrozMinInput,
    disWinMaxInput,
    disWinMinInput,
    selectedReferenceBandData,
    selectedReferenceBandsCalculation?.id,
    setDisentanglementInputsFromAbsolute,
  ]);

  const commitDisentanglementInputRanges = useCallback(() => {
    const disCandidate: [number, number] = [
      parseFiniteInput(disWinMinInput) ?? effectiveDisWindowAbsolute[0],
      parseFiniteInput(disWinMaxInput) ?? effectiveDisWindowAbsolute[1],
    ];
    const frozenCandidate: [number, number] = [
      parseFiniteInput(disFrozMinInput) ?? effectiveFrozenWindowAbsolute[0],
      parseFiniteInput(disFrozMaxInput) ?? effectiveFrozenWindowAbsolute[1],
    ];
    const nextDisWindow = clampDisWindowAroundFrozen(disCandidate, frozenCandidate);
    const nextFrozenWindow = clampFrozenInsideDisWindow(frozenCandidate, nextDisWindow);
    setDisentanglementInputsFromAbsolute(nextDisWindow, nextFrozenWindow);
  }, [
    clampDisWindowAroundFrozen,
    clampFrozenInsideDisWindow,
    disFrozMaxInput,
    disFrozMinInput,
    disWinMaxInput,
    disWinMinInput,
    effectiveDisWindowAbsolute,
    effectiveFrozenWindowAbsolute,
    setDisentanglementInputsFromAbsolute,
  ]);

  const handleDisentanglementInputKeyDown = useCallback((event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key !== "Enter") return;
    event.preventDefault();
    commitDisentanglementInputRanges();
  }, [commitDisentanglementInputRanges]);

  const referenceWindowOverlays = useMemo<BandPlotWindowOverlay[]>(() => {
    if (!selectedReferenceBandData) {
      return [];
    }
    return [
      {
        id: "disentanglement",
        min: absoluteEnergyToDisplay(effectiveDisWindowAbsolute[0]),
        max: absoluteEnergyToDisplay(effectiveDisWindowAbsolute[1]),
        color: "#0ea5e9",
        side: "left",
        label: "Disentanglement",
      },
      {
        id: "frozen",
        min: absoluteEnergyToDisplay(effectiveFrozenWindowAbsolute[0]),
        max: absoluteEnergyToDisplay(effectiveFrozenWindowAbsolute[1]),
        color: "#f59e0b",
        side: "right",
        label: "Frozen",
      },
    ];
  }, [
    absoluteEnergyToDisplay,
    effectiveDisWindowAbsolute,
    effectiveFrozenWindowAbsolute,
    selectedReferenceBandData,
  ]);

  const handleReferenceWindowOverlayChange = useCallback((update: BandPlotWindowOverlayUpdate) => {
    const candidateAbsoluteWindow = toAscendingPair(
      displayEnergyToAbsolute(update.min),
      displayEnergyToAbsolute(update.max),
    );
    if (update.id === "disentanglement") {
      const nextDisWindow = clampDisWindowAroundFrozen(
        candidateAbsoluteWindow,
        effectiveFrozenWindowAbsolute,
      );
      setDisentanglementInputsFromAbsolute(nextDisWindow, effectiveFrozenWindowAbsolute);
      return;
    }
    if (update.id === "frozen") {
      const nextFrozenWindow = clampFrozenInsideDisWindow(
        candidateAbsoluteWindow,
        effectiveDisWindowAbsolute,
      );
      setDisentanglementInputsFromAbsolute(effectiveDisWindowAbsolute, nextFrozenWindow);
    }
  }, [
    clampDisWindowAroundFrozen,
    clampFrozenInsideDisWindow,
    displayEnergyToAbsolute,
    effectiveDisWindowAbsolute,
    effectiveFrozenWindowAbsolute,
    setDisentanglementInputsFromAbsolute,
  ]);

  const buildViewerOverlayOptions = useCallback(async (): Promise<WannierBandOverlayOption[]> => {
    const sourceScfId = String(selectedScf?.id ?? "").trim();
    const selectedKPath = normalizeSavedKPath(kPath.map((point) => point.label).join(" → "));
    if (!sourceScfId || !selectedKPath) {
      return [];
    }

    const matchingBandRuns = bandCalculationSummaries.filter((candidate) => {
      if (candidate.calc_type !== "bands") {
        return false;
      }
      const candidateSourceScfId = String(candidate.parameters?.source_scf_id ?? "").trim();
      const candidateKPath = normalizeSavedKPath(candidate.parameters?.k_path);
      return candidateSourceScfId === sourceScfId && candidateKPath === selectedKPath;
    });
    if (matchingBandRuns.length === 0) {
      return [];
    }

    const hydratedBandDetails: CalculationRun[] = [];
    const settled = await Promise.allSettled(
      matchingBandRuns.map(async (candidate) => {
        let detail = detailedBandsById[candidate.id] ?? candidate;
        if (!detail.result?.band_data) {
          try {
            detail = await invoke<CalculationRun>("get_project_calculation", {
              projectId,
              calcId: candidate.id,
            });
          } catch {
            return null;
          }
        }
        if (detail.result?.band_data) {
          hydratedBandDetails.push(detail);
        }
        const bandData = detail.result?.band_data ?? candidate.result?.band_data ?? null;
        if (!bandData) {
          return null;
        }
        const startedAt = detail.started_at ?? candidate.started_at;
        return {
          id: candidate.id,
          label: `Bands · ${new Date(startedAt).toLocaleString()}`,
          data: bandData,
          fermiEnergy: detail.result?.fermi_energy ?? candidate.result?.fermi_energy ?? null,
          startedAt,
        };
      }),
    );

    if (hydratedBandDetails.length > 0) {
      setDetailedBandsById((prev) => {
        const next = { ...prev };
        let changed = false;
        for (const detail of hydratedBandDetails) {
          if (!next[detail.id]?.result?.band_data && detail.result?.band_data) {
            next[detail.id] = detail;
            changed = true;
          }
        }
        return changed ? next : prev;
      });
    }

    return settled
      .flatMap((entry) => (entry.status === "fulfilled" && entry.value ? [entry.value] : []))
      .sort((left, right) => right.startedAt.localeCompare(left.startedAt))
      .map(({ startedAt: _startedAt, ...overlay }) => overlay);
  }, [bandCalculationSummaries, detailedBandsById, kPath, projectId, selectedScf?.id]);

  const derivedNumWann = useMemo(
    () => projectionDrafts.reduce(
      (sum, item) => sum + projectionContributionCount(item, wannierCellAtomSymbols),
      0,
    ),
    [projectionDrafts, wannierCellAtomSymbols],
  );

  const resolvedNumWann = useMemo(() => {
    const override = Number.parseInt(numWannOverrideInput.trim(), 10);
    if (Number.isInteger(override) && override > 0) {
      return override;
    }
    return derivedNumWann;
  }, [derivedNumWann, numWannOverrideInput]);

  const projectionSummary = useMemo(
    () => projectionDrafts.map((item) => {
      if (item.targetType === "site" && item.siteIndex != null) {
        const site = crystalData.atom_sites[item.siteIndex];
        return `${site?.label || `site ${item.siteIndex + 1}`}: ${item.orbital}`;
      }
      const siteCount = countMatchingElementSites(item.symbol, wannierCellAtomSymbols);
      const totalCount = projectionContributionCount(item, wannierCellAtomSymbols);
      return `${item.symbol}: ${item.orbital} (${siteCount} site${siteCount === 1 ? "" : "s"}, ${totalCount} functions)`;
    }),
    [projectionDrafts, crystalData, wannierCellAtomSymbols],
  );
  const sourceElectronCount = useMemo(
    () => parseElectronCount(selectedScf?.result?.raw_output),
    [selectedScf],
  );
  const sourceNbnd = useMemo(() => {
    const value = Number(selectedScf?.parameters?.nbnd);
    return Number.isInteger(value) && value > 0 ? value : null;
  }, [selectedScf]);
  const minOccupiedBands = useMemo(
    () => sourceElectronCount != null ? Math.ceil(sourceElectronCount / 2) : null,
    [sourceElectronCount],
  );
  const minimumNumBands = useMemo(
    () => Math.max(resolvedNumWann, minOccupiedBands ?? 0, 1),
    [minOccupiedBands, resolvedNumWann],
  );
  const recommendedNumBands = useMemo(() => {
    if (!selectedScf) return null;
    return Math.max(minimumNumBands, sourceNbnd ?? 0);
  }, [minimumNumBands, selectedScf, sourceNbnd]);
  const effectiveNumBandsInput = useMemo(
    () => (numBandsMode === "auto" ? (recommendedNumBands != null ? String(recommendedNumBands) : "") : numBandsInput),
    [numBandsInput, numBandsMode, recommendedNumBands],
  );
  const manualNumBandsForHpc = useMemo(() => {
    if (numBandsMode !== "manual") {
      return null;
    }
    const parsed = Number(effectiveNumBandsInput);
    return Number.isInteger(parsed) && parsed > 0 ? parsed : null;
  }, [effectiveNumBandsInput, numBandsMode]);
  const numBandsHelpText = useMemo(() => {
    if (!selectedScf || recommendedNumBands == null) {
      return "Select a source SCF to derive num_bands.";
    }
    if (sourceNbnd != null && sourceNbnd > minimumNumBands) {
      return `Auto uses ${recommendedNumBands} bands because the source SCF already set nbnd=${sourceNbnd}, which is larger than the minimum requirement ${minimumNumBands}.`;
    }
    if (minOccupiedBands != null && sourceElectronCount != null) {
      return `Auto uses ${recommendedNumBands} bands = max(${minOccupiedBands} occupied bands from ${sourceElectronCount} electrons, num_wann=${resolvedNumWann}). Increase it only if you want a larger disentanglement space.`;
    }
    return `Auto uses ${recommendedNumBands} bands = max(num_wann=${resolvedNumWann}${sourceNbnd != null ? `, source nbnd=${sourceNbnd}` : ""}). Increase it only if you want a larger disentanglement space.`;
  }, [
    minOccupiedBands,
    minimumNumBands,
    recommendedNumBands,
    resolvedNumWann,
    selectedScf,
    sourceElectronCount,
    sourceNbnd,
  ]);
  const handleKPathChange = useCallback((newPath: KPathPoint[]) => {
    setKPath(applyTotalKPoints(newPath, totalKPointsTarget));
  }, [totalKPointsTarget]);
  useEffect(() => {
    setKPath((prevPath) => applyTotalKPoints(prevPath, totalKPointsTarget));
  }, [totalKPointsTarget]);
  const kPathSegmentCount = useMemo(
    () => getConnectedSegmentIndices(kPath).length,
    [kPath],
  );
  const totalInterpolatedPoints = useMemo(
    () => kPath.reduce((sum, point) => sum + (Number(point.npoints) || 0), 0),
    [kPath],
  );
  const minimumTotalKPoints = Math.max(1, kPathSegmentCount);
  const viewerPointsPerSegment = clampInt(
    totalKPointsTarget / minimumTotalKPoints,
    1,
    MAX_VIEWER_POINTS_PER_SEGMENT,
  );
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
  const parameterValidationError = useMemo(() => {
    if (!selectedScf) {
      return "Select a source SCF first.";
    }
    if (selectedIssues.length > 0) {
      return selectedIssues.join(". ");
    }
    if (projectionDrafts.length === 0) {
      return "Add at least one Wannier projection.";
    }
    if (resolvedNumWann !== derivedNumWann) {
      return "num_wann override must match the rendered projection count.";
    }
    if (kPath.length < 2) {
      return "Select an interpolation k-path.";
    }

    try {
      const kx = parseOptionalPositiveInt(kGridInput[0], "k-grid x");
      const ky = parseOptionalPositiveInt(kGridInput[1], "k-grid y");
      const kz = parseOptionalPositiveInt(kGridInput[2], "k-grid z");
      const totalKPoints = parseOptionalPositiveInt(totalKPointsInput, "total interpolation k-points");
      const numBands = parseOptionalPositiveInt(effectiveNumBandsInput, "num_bands");
      if (kx == null || ky == null || kz == null || totalKPoints == null || numBands == null) {
        return "Use positive integers for the mesh, total interpolation k-points, and num_bands.";
      }
      if (totalKPoints < minimumTotalKPoints) {
        return `Total interpolation k-points must be at least ${minimumTotalKPoints} for the selected path.`;
      }
      if (minOccupiedBands != null && numBands < minOccupiedBands) {
        return `num_bands must be at least ${minOccupiedBands} for this scalar source (${sourceElectronCount} electrons => ${minOccupiedBands} occupied bands).`;
      }
      if (numBands < resolvedNumWann) {
        return "num_bands must be greater than or equal to num_wann.";
      }
      if (isHpcMode && numBandsMode === "manual") {
        const taskLimitError = validateHpcTasksWithinBandCount(hpcResources, numBands, "num_bands");
        if (taskLimitError) {
          return taskLimitError;
        }
      }

      if (numBands > resolvedNumWann) {
        const disWinMin = parseOptionalNumber(disWinMinInput, "dis_win_min");
        const disWinMax = parseOptionalNumber(disWinMaxInput, "dis_win_max");
        const disFrozMin = parseOptionalNumber(disFrozMinInput, "dis_froz_min");
        const disFrozMax = parseOptionalNumber(disFrozMaxInput, "dis_froz_max");
        if (disWinMin != null && disWinMax != null && disWinMin > disWinMax) {
          return "dis_win_min must be less than or equal to dis_win_max.";
        }
        if (disFrozMin != null && disFrozMax != null && disFrozMin > disFrozMax) {
          return "dis_froz_min must be less than or equal to dis_froz_max.";
        }
        if (disWinMin != null && disFrozMin != null && disFrozMin < disWinMin) {
          return "Frozen window must lie inside the disentanglement window.";
        }
        if (disWinMax != null && disFrozMax != null && disFrozMax > disWinMax) {
          return "Frozen window must lie inside the disentanglement window.";
        }
      }
    } catch (err) {
      return err instanceof Error ? err.message : String(err);
    }

    return null;
  }, [
    derivedNumWann,
    disFrozMaxInput,
    disFrozMinInput,
    disWinMaxInput,
    disWinMinInput,
    kGridInput,
    kPath,
    effectiveNumBandsInput,
    hpcResources,
    isHpcMode,
    numBandsMode,
    projectionDrafts.length,
    resolvedNumWann,
    selectedIssues,
    selectedScf,
    minimumTotalKPoints,
    minOccupiedBands,
    sourceElectronCount,
    totalKPointsInput,
  ]);

  const hpcCommandLines = useMemo(() => {
    const remoteWannier = (activeHpcProfile?.remote_wannier90_path || "wannier90.x").trim() || "wannier90.x";
    const seedname = sanitizeSeedname(seednameInput);
    return [
      `"${remoteWannier}" -pp ${seedname} > wannier90_pre.out 2>&1`,
      buildHpcQeInputCommandLine(activeHpcProfile, "pw.x", "nscf.in", "nscf.out", undefined, hpcResources.resource_type),
      buildHpcQeInputCommandLine(activeHpcProfile, "pw2wannier90.x", "pw2wan.in", "pw2wan.out", undefined, hpcResources.resource_type),
      `"${remoteWannier}" ${seedname} > wannier90.out 2>&1`,
    ];
  }, [activeHpcProfile, hpcResources.resource_type, seednameInput]);

  const resultQualityIssues = useMemo(
    () => getWannierQualityIssues(result, output, selectedScf?.result?.fermi_energy ?? null),
    [output, result, selectedScf],
  );

  const addProjection = useCallback(() => {
    const firstElement = getBaseElement(crystalData.atom_sites[0]?.type_symbol || "X");
    const defaultOrbital = getDefaultProjectionOrbital(firstElement);
    setProjectionDrafts((prev) => [
      ...prev,
      {
        id: makeProjectionId(),
        targetType: "element",
        symbol: firstElement,
        orbital: defaultOrbital,
        siteIndex: null,
      },
    ]);
  }, [crystalData]);

  const updateProjection = useCallback((id: string, patch: Partial<ProjectionDraft>) => {
    setProjectionDrafts((prev) => prev.map((item) => (item.id === id ? { ...item, ...patch } : item)));
  }, []);

  const removeProjection = useCallback((id: string) => {
    setProjectionDrafts((prev) => prev.filter((item) => item.id !== id));
  }, []);

  async function ensureSymmetryTransform(): Promise<SymmetryTransformResult | null> {
    if (symmetryTransform) {
      return symmetryTransform;
    }
    try {
      const transform = await analyzeCrystalSymmetry(crystalData);
      setSymmetryTransform(transform);
      setSymmetryError(null);
      return transform;
    } catch (err) {
      const message = String(err);
      setSymmetryError(message);
      return null;
    }
  }

  async function buildTaskPlan() {
    if (!selectedScf) {
      throw new Error("Select a source SCF first.");
    }
    if (selectedIssues.length > 0) {
      throw new Error(selectedIssues.join(". "));
    }
    if (projectionDrafts.length === 0) {
      throw new Error("Add at least one Wannier projection.");
    }
    if (resolvedNumWann !== derivedNumWann) {
      throw new Error("num_wann override must match the rendered projection count.");
    }
    if (kPath.length < 2) {
      throw new Error("Select an interpolation k-path.");
    }

    const kGrid = [
      parseOptionalPositiveInt(kGridInput[0], "k-grid x") ?? 0,
      parseOptionalPositiveInt(kGridInput[1], "k-grid y") ?? 0,
      parseOptionalPositiveInt(kGridInput[2], "k-grid z") ?? 0,
    ] as [number, number, number];
    const totalKPoints = parseOptionalPositiveInt(totalKPointsInput, "total interpolation k-points") ?? 0;
    const numBands = parseOptionalPositiveInt(effectiveNumBandsInput, "num_bands") ?? 0;
    const seedname = sanitizeSeedname(seednameInput);
    if (totalKPoints < minimumTotalKPoints) {
      throw new Error(`Total interpolation k-points must be at least ${minimumTotalKPoints} for the selected path.`);
    }
    if (minOccupiedBands != null && numBands < minOccupiedBands) {
      throw new Error(
        `num_bands must be at least ${minOccupiedBands} for this scalar source (${sourceElectronCount} electrons => ${minOccupiedBands} occupied bands).`,
      );
    }
    if (numBands < resolvedNumWann) {
      throw new Error("num_bands must be greater than or equal to num_wann.");
    }
    if (isHpcMode && numBandsMode === "manual") {
      const taskLimitError = validateHpcTasksWithinBandCount(hpcResources, numBands, "num_bands");
      if (taskLimitError) {
        throw new Error(taskLimitError);
      }
    }

    const params = selectedScf.parameters || {};
    const ecutwfc = Number(params.ecutwfc);
    const ecutrho = Number(params.ecutrho);
    const nspin = Number(params.nspin) || 1;
    const lspinorb = Boolean(params.lspinorb);
    const noncolin = nspin === 4 || Boolean(params.noncolin) || lspinorb;
    const occupations = normalizeOccupations(params.occupations);
    const smearing = normalizeSmearing(params.smearing, resolvedDefaultSmearing);
    const degauss = Number(params.degauss);
    const convThr = Number(params.conv_thr);
    const mixingBeta = Number(params.mixing_beta);
    const baseElements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
    const species = baseElements.map((element) => ({
      symbol: element,
      mass: ELEMENT_MASSES[element] || 1,
      pseudopotential: sourcePseudoMap[element],
      ...getMagneticSpeciesFields(params, element),
    }));
    if (species.some((entry) => !entry.pseudopotential)) {
      throw new Error("Source SCF pseudopotential metadata is incomplete. Re-run the SCF before starting Wannier.");
    }

    const resolvedSymmetry = await ensureSymmetryTransform();
    const structureForNscf = resolveSavedScfStructure(params);
    const isOptimizedSourceScf =
      String(params.cell_representation || "").toLowerCase() === "optimized_source";
    if (isOptimizedSourceScf && !structureForNscf) {
      throw new Error(
        "Selected SCF was run from an optimized structure, but its saved structure metadata is missing. Re-run the SCF from the optimized source and try again.",
      );
    }
    const sourceUsesPrimitive = sourceScfUsesPrimitiveCell(params, conventionalAtomCount);
    const canUseSymmetryPrimitive =
      sourceUsesPrimitive &&
      resolvedSymmetry !== null &&
      resolvedSymmetry.standardizedPrimitiveAtoms.length > 0;
    if (sourceUsesPrimitive && !canUseSymmetryPrimitive) {
      throw new Error(
        "Selected SCF was run in a primitive cell, but symmetry conversion data is unavailable.",
      );
    }

    const conventionalLattice = buildConventionalLatticeFromCrystalData(crystalData);
    const context = resolvePathTransformContext(crystalData, resolvedSymmetry);
    const effectiveRhombohedralConvention = kPathRhombohedralConvention
      ?? defaultRhombohedralConventionForSetting(context.rhombohedralSetting ?? null);
    const converters = createPathCoordinateConverters(context, resolvedSymmetry);
    const inferredBravais =
      canUseSymmetryPrimitive && resolvedSymmetry
        ? inferQeBravaisCellFromCif(crystalData, resolvedSymmetry)
        : null;

    let baseCalculation;
    let transformedPath;
    const pseudoDir = isHpcMode
      ? resolveProfileRemotePseudoDir(activeHpcProfile, hpcResources.resource_type)
      : qePath.replace(/\/bin\/?$/, "/pseudo");

    if (isOptimizedSourceScf && structureForNscf) {
      baseCalculation = {
        calculation: "scf",
        prefix: params.prefix || "qcortado_scf",
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
          position_units: structureForNscf.position_units || "crystal",
          ecutwfc,
          ecutrho,
          nbnd: null,
          nspin,
          noncolin,
          lspinorb,
          occupations,
          smearing,
          degauss: Number.isFinite(degauss) && degauss > 0 ? degauss : null,
          nosym: false,
          noinv: false,
        },
        kpoints: { type: "gamma" },
        conv_thr: Number.isFinite(convThr) && convThr > 0 ? convThr : 1e-8,
        mixing_beta: Number.isFinite(mixingBeta) && mixingBeta > 0 ? mixingBeta : 0.7,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: "high",
      };
      transformedPath = mapPathCoordinates(
        kPath,
        canUseSymmetryPrimitive && resolvedSymmetry
          ? converters.toSymmetryPrimitiveCoords
          : converters.toInputConventionalCoords,
      ).map((point) => ({
        label: point.label,
        coords: point.coords,
        npoints: point.npoints,
      }));
    } else if (canUseSymmetryPrimitive && resolvedSymmetry) {
      baseCalculation = {
        calculation: "scf",
        prefix: params.prefix || "qcortado_scf",
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
          nbnd: null,
          nspin,
          noncolin,
          lspinorb,
          occupations,
          smearing,
          degauss: Number.isFinite(degauss) && degauss > 0 ? degauss : null,
          nosym: false,
          noinv: false,
        },
        kpoints: { type: "gamma" },
        conv_thr: Number.isFinite(convThr) && convThr > 0 ? convThr : 1e-8,
        mixing_beta: Number.isFinite(mixingBeta) && mixingBeta > 0 ? mixingBeta : 0.7,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: "high",
      };
      transformedPath = mapPathCoordinates(
        kPath,
        converters.toSymmetryPrimitiveCoords,
      ).map((point) => ({
        label: point.label,
        coords: point.coords,
        npoints: point.npoints,
      }));
    } else {
      baseCalculation = {
        calculation: "scf",
        prefix: params.prefix || "qcortado_scf",
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
          nbnd: null,
          nspin,
          noncolin,
          lspinorb,
          occupations,
          smearing,
          degauss: Number.isFinite(degauss) && degauss > 0 ? degauss : null,
          nosym: false,
          noinv: false,
        },
        kpoints: { type: "gamma" },
        conv_thr: Number.isFinite(convThr) && convThr > 0 ? convThr : 1e-8,
        mixing_beta: Number.isFinite(mixingBeta) && mixingBeta > 0 ? mixingBeta : 0.7,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: "high",
      };
      transformedPath = mapPathCoordinates(
        kPath,
        converters.toInputConventionalCoords,
      ).map((point) => ({
        label: point.label,
        coords: point.coords,
        npoints: point.npoints,
      }));
    }

    const projectionSpecs = projectionDrafts.map((item) => {
      if (item.targetType === "site" && item.siteIndex != null) {
        const site = crystalData.atom_sites[item.siteIndex];
        return {
          target_type: "site",
          symbol: getBaseElement(site.type_symbol),
          orbital: item.orbital,
          site_index: item.siteIndex + 1,
          fractional_position: [site.fract_x, site.fract_y, site.fract_z] as [number, number, number],
        };
      }
      return {
        target_type: "element",
        symbol: item.symbol,
        orbital: item.orbital,
      };
    });

    const disentanglement =
      numBands > resolvedNumWann
        ? {
          exclude_bands: parseExcludeBandsInput(excludeBandsInput),
          dis_win_min: parseOptionalNumber(disWinMinInput, "dis_win_min"),
          dis_win_max: parseOptionalNumber(disWinMaxInput, "dis_win_max"),
          dis_froz_min: parseOptionalNumber(disFrozMinInput, "dis_froz_min"),
          dis_froz_max: parseOptionalNumber(disFrozMaxInput, "dis_froz_max"),
        }
        : null;
    if (disentanglement) {
      if (
        disentanglement.dis_win_min != null &&
        disentanglement.dis_win_max != null &&
        disentanglement.dis_win_min > disentanglement.dis_win_max
      ) {
        throw new Error("dis_win_min must be less than or equal to dis_win_max.");
      }
      if (
        disentanglement.dis_froz_min != null &&
        disentanglement.dis_froz_max != null &&
        disentanglement.dis_froz_min > disentanglement.dis_froz_max
      ) {
        throw new Error("dis_froz_min must be less than or equal to dis_froz_max.");
      }
      if (
        disentanglement.dis_win_min != null &&
        disentanglement.dis_froz_min != null &&
        disentanglement.dis_froz_min < disentanglement.dis_win_min
      ) {
        throw new Error("Frozen window must lie inside the disentanglement window.");
      }
      if (
        disentanglement.dis_win_max != null &&
        disentanglement.dis_froz_max != null &&
        disentanglement.dis_froz_max > disentanglement.dis_win_max
      ) {
        throw new Error("Frozen window must lie inside the disentanglement window.");
      }
    }

    return {
      taskLabel: `Wannier - ${crystalData.formula_sum || crystalData.formula_structural || "structure"}`,
      taskParams: {
        config: {
          base_calculation: baseCalculation,
          k_grid: kGrid,
          num_wann: resolvedNumWann,
          num_bands: numBands,
          seedname,
          projections: projectionSpecs,
          band_path: {
            k_path: transformedPath,
            bands_num_points: 0,
            total_k_points_target: totalKPoints,
          },
          disentanglement,
          project_id: projectId,
          scf_calc_id: selectedScf.id,
          source_metadata: {
            cell_representation: params.cell_representation || null,
            nspin,
            noncolin,
            lspinorb,
            lda_plus_u: Boolean(params.lda_plus_u),
            vdw_corr: params.vdw_corr || null,
            electron_count: sourceElectronCount,
          },
        },
        workingDir: WANNIER_WORK_DIR,
        mpiConfig: !isHpcMode ? { enabled: mpiEnabled, nprocs: mpiProcs } : null,
        executionTarget: buildExecutionTarget(
          executionMode,
          activeHpcProfile?.id ?? null,
          isHpcMode ? hpcResources : null,
          false,
        ),
      },
      saveParameters: {
        source_scf_id: selectedScf.id,
        seedname,
        k_grid: kGrid,
        num_wann: resolvedNumWann,
        num_bands: numBands,
        projection_summary: projectionSummary,
        exclude_bands: disentanglement?.exclude_bands ?? [],
        dis_win_min: disentanglement?.dis_win_min ?? null,
        dis_win_max: disentanglement?.dis_win_max ?? null,
        dis_froz_min: disentanglement?.dis_froz_min ?? null,
        dis_froz_max: disentanglement?.dis_froz_max ?? null,
        k_path: kPath.map((point) => point.label).join(" → "),
        total_k_points_target: totalKPoints,
        scf_fermi_energy: selectedScf.result?.fermi_energy ?? null,
        cell_representation: isOptimizedSourceScf
          ? "optimized_source"
          : (canUseSymmetryPrimitive ? "primitive_spglib" : "conventional_input"),
        ecutwfc: params.ecutwfc,
        nspin: params.nspin,
        noncolin: params.noncolin,
        lspinorb: params.lspinorb,
        starting_magnetization: params.starting_magnetization ?? null,
        starting_magnetization_theta: params.starting_magnetization_theta ?? params.starting_magnetization_angle1 ?? params.theta ?? params.angle1 ?? null,
        starting_magnetization_phi: params.starting_magnetization_phi ?? params.starting_magnetization_angle2 ?? params.phi ?? params.angle2 ?? null,
        lda_plus_u: params.lda_plus_u,
        vdw_corr: params.vdw_corr,
        symmetry_spacegroup: resolvedSymmetry?.spacegroupNumber ?? null,
        symmetry_hall_number: resolvedSymmetry?.hallNumber ?? null,
        k_path_convention: context.centering === "R" ? effectiveRhombohedralConvention : null,
      },
    };
  }

  async function runCalculation() {
    if (hasBlockingExternalTask) {
      setError("Another local task is running. Wait for it to finish or switch to HPC mode.");
      return;
    }

    setError(null);
    setIsRunning(true);
    setResult(null);
    setIsSaved(false);
    setShowOutput(true);
    setOutput("");
    setOutputLineCount(0);
    setProgress(defaultProgressState("Wannier90"));
    setStep("run");
    const startedAt = new Date().toISOString();
    setCalcStartTime(startedAt);

    try {
      const plan = await buildTaskPlan();
      const hpcSaveSpec = isHpcMode
        ? {
          projectId,
          cifId,
          workingDir: WANNIER_WORK_DIR,
          calcType: "wannier" as const,
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        }
        : null;
      const taskId = await taskContext.startTask("wannier", plan.taskParams, plan.taskLabel, hpcSaveSpec);
      setActiveTaskId(taskId);
      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask.error || "Wannier calculation failed.");
      }
      const nextResult = finalTask.result as WannierResult;
      const outputText = finalTask.output.join("\n");
      const visibleTask = taskContext.getTask(taskId);
      setResult(nextResult);
      setOutput(visibleTask?.outputText ?? outputText);
      setOutputLineCount(visibleTask?.outputLineCount ?? countVisibleOutputLines(outputText));
      setIsRunning(false);
      setStep("results");
      const hpcSaveParams = (isHpcMode || finalTask.hpc.backend === "hpc")
        ? {
          execution_backend: "hpc",
          hpc_profile_id: activeHpcProfile?.id ?? null,
          hpc_resource_type: finalTask.hpc.hpc_resource_type ?? hpcResources.resource_type,
          remote_job_id: finalTask.hpc.remote_job_id ?? null,
          scheduler_state: finalTask.hpc.scheduler_state ?? null,
          remote_node: finalTask.hpc.remote_node ?? null,
          remote_workdir: finalTask.hpc.remote_workdir ?? null,
          remote_project_path: finalTask.hpc.remote_project_path ?? null,
          remote_storage_bytes: finalTask.hpc.remote_storage_bytes ?? null,
        }
        : {};
      await invoke("save_calculation", {
        projectId,
        cifId,
        calcData: {
          calc_type: "wannier",
          parameters: {
            ...plan.saveParameters,
            total_spread: nextResult.total_spread ?? null,
            n_wann: nextResult.num_wann,
            n_bands: nextResult.num_bands,
            total_k_points: nextResult.band_data.n_kpoints,
            ...hpcSaveParams,
          },
          result: {
            converged: nextResult.convergence?.converged ?? true,
            total_energy: null,
            fermi_energy: selectedScf?.result?.fermi_energy ?? null,
            n_scf_steps: null,
            wall_time_seconds: null,
            raw_output: outputText,
            band_data: nextResult.band_data,
            wannier_data: nextResult,
          },
          started_at: startedAt,
          completed_at: new Date().toISOString(),
          input_content: "",
          output_content: outputText,
          tags: [],
        },
        workingDir: finalTask.hpc.local_sync_dir ?? WANNIER_WORK_DIR,
      });
      setIsSaved(true);
    } catch (err) {
      setError(String(err));
      setIsRunning(false);
    }
  }

  const renderSourceStep = () => {
    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Wannier90 requires a converged scalar SCF calculation.
            Please run an SCF calculation first.
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
            <label htmlFor="wannier-scf-sort">Sort SCFs</label>
            <select
              id="wannier-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose the converged scalar SCF that will supply the restart density and source metadata for Wannierization.
        </p>

        {readyScfs.length === 0 && (
          <div className="warning-banner">
            No Wannier-ready SCFs found. Wannier v1 requires a scalar source with no SOC, noncollinearity, DFT+U, or certified vdW corrections.
          </div>
        )}

        <div className="scf-list">
          {sortedReadyScfs.map((scf) => {
            const scfName = getCalculationName(scf);
            return (
              <div
                key={scf.id}
                className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""}`}
                onClick={() => {
                  setSelectedScf(scf);
                  setError(null);
                }}
              >
                <div className="scf-option-header">
                  <input
                    type="radio"
                    checked={selectedScf?.id === scf.id}
                    onChange={() => {
                      setSelectedScf(scf);
                      setError(null);
                    }}
                  />
                  {scfName && (
                    <span className="scf-name" title={formatCalculationSourceLabel(scf)}>
                      {scfName}
                    </span>
                  )}
                  <span className="scf-date">
                    {new Date(scf.started_at).toLocaleDateString()}
                  </span>
                </div>
                <div className="scf-details">
                  {typeof scf.result?.fermi_energy === "number" && (
                    <span>EF = {scf.result.fermi_energy.toFixed(3)} eV</span>
                  )}
                </div>
                <div className="calc-tags">
                  {getCalculationTags(scf).map((tag, i) => (
                    <span
                      key={`${tag.label}-${i}`}
                      className={`calc-tag calc-tag-${tag.type}${tag.label.trim().toUpperCase() === "HPC" ? " calc-tag-hpc" : ""}`}
                    >
                      {tag.label}
                    </span>
                  ))}
                </div>
              </div>
            );
          })}
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Cancel
          </button>
          <button
            className="primary-button"
            disabled={!selectedScf}
            onClick={() => {
              if (!selectedScf) return;
              void prepareSelectedScfSource(selectedScf)
                .then(() => setStep("mesh"))
                .catch((err) => setError(String(err)));
            }}
          >
            Next: K-Mesh
          </button>
        </div>
      </div>
    );
  };

  const renderMeshStep = () => {
    return (
      <div className="wizard-step kpath-step">
        <h3>K-Mesh and Interpolation Path</h3>
        <p className="step-description">
          Define the explicit zero-shift Monkhorst-Pack mesh used for the Wannier NSCF and the interpolation path passed to Wannier90.
        </p>

        <div className="param-section">
          <h4>NSCF + Interpolation Controls</h4>
          <div className="param-grid">
            <div className="param-row">
              <label>Uniform NSCF mesh</label>
              <div className="kgrid-inputs">
                {kGridInput.map((value, index) => (
                  <input
                    key={index}
                    type="number"
                    min={1}
                    value={value}
                    onChange={(event) => {
                      const next = [...kGridInput] as [string, string, string];
                      next[index] = event.target.value;
                      setKGridInput(next);
                    }}
                  />
                ))}
              </div>
            </div>
            <div className="param-row">
              <label>
                Total interpolated k-points
                <InfoTooltip text="QCortado target for the full interpolation path. QCortado converts this to Wannier90 input `bands_num_points` internally; because Wannier90 scales each segment by relative length, the final plotted total may differ slightly from the requested target." />
              </label>
              <input
                type="number"
                min={minimumTotalKPoints}
                max={MAX_TOTAL_K_POINTS}
                value={totalKPointsInput}
                onChange={(event) => setTotalKPointsInput(event.target.value)}
                onBlur={commitTotalKPointsInput}
                onKeyDown={(event) => {
                  if (event.key === "Enter") {
                    event.preventDefault();
                    commitTotalKPointsInput();
                  }
                }}
              />
            </div>
            <div className="param-row">
              <label>
                Calculation file prefix
                <InfoTooltip text="Wannier90 input: `seedname`. QCortado uses this prefix when writing the `.win`, `.wout`, `_band.dat`, and related Wannier files." />
              </label>
              <input
                type="text"
                value={seednameInput}
                onChange={(event) => setSeednameInput(event.target.value)}
                spellCheck={false}
              />
            </div>
          </div>
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

        {symmetryError && <div className="warning-banner">{symmetryError}</div>}
        {selectedIssues.length > 0 && (
          <div className="warning-banner">{selectedIssues.join(". ")}</div>
        )}

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>Source SCF:</span>
            <span>{formatCalculationSourceLabel(selectedScf)}</span>
          </div>
          <div className="summary-row">
            <span>NSCF mesh:</span>
            <span>{kGridInput.join("×")}</span>
          </div>
          <div className="summary-row">
            <span>Total k-point target:</span>
            <span>{totalKPointsTarget}</span>
          </div>
          <div className="summary-row">
            <span>Path nodes:</span>
            <span>{kPath.map((point) => point.label).join(" → ") || "Not selected"}</span>
          </div>
          <div className="summary-row">
            <span>Distributed path points:</span>
            <span>{totalInterpolatedPoints}</span>
          </div>
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button
            className="primary-button"
            disabled={!selectedScf}
            onClick={() => setStep("projections")}
          >
            Next: Projections
          </button>
        </div>
      </div>
    );
  };

  const renderProjectionStep = () => {
    return (
      <div className="wizard-step parameters-step parameters-step-fullscreen">
        <div className="wannier-projection-workbench">
          <div className="wannier-projection-controls-pane">
            <h3>Projections and Disentanglement</h3>
            <p className="step-description">
              Define the initial orbital projections and tune disentanglement/frozen windows against the saved bands preview.
            </p>
            <div className="info-banner">
              Element projections apply to every matching atom in the Wannier cell. Use site-targeted projections if you want only one atom.
            </div>
            <details className="wannier-guidance">
              <summary>How to choose transport-ready projections</summary>
              <div className="guidance-body">
                <p>Start from the orbitals whose source bands actually cross or sit close to E_F. If the interpolated Wannier manifold misses E_F, BoltzWann transport will collapse toward zero.</p>
                <p>For simple isolated manifolds, keep `num_bands = num_wann`. For entangled metallic manifolds, increase `num_bands`, then use `dis_win_*` to capture the full target space and `dis_froz_*` to pin the bands you must reproduce near E_F.</p>
                <p>QCortado shows each element's neutral-atom electron configuration as a reference only. It does not restrict the orbital templates you can choose.</p>
              </div>
            </details>
            {structureElectronConfigs.length > 0 && (
              <div className="info-banner">
                Neutral atom references:&nbsp;
                {structureElectronConfigs.map(({ element, config }, index) => (
                  <span key={element}>
                    <strong>{element}</strong>: {config ? renderElectronConfiguration(config.compact) : "unavailable"}
                    {index < structureElectronConfigs.length - 1 ? "; " : ""}
                  </span>
                ))}
              </div>
            )}

            <div className="param-section">
              <div className="source-step-header">
                <h4>
                  Initial Projections
                  <InfoTooltip text="Choose orbitals that span the bands you want the Wannier interpolation to reproduce. For transport, prioritize orbitals contributing near the source SCF Fermi level." />
                </h4>
                <button className="secondary-button" type="button" onClick={addProjection}>
                  Add Projection
                </button>
              </div>

              {projectionDrafts.length === 0 && (
                <div className="warning-banner">
                  Add at least one element- or site-targeted orbital template.
                </div>
              )}

              {projectionDrafts.map((item) => {
                const allowedOrbitalOptions = getAllowedOrbitalOptionsForDraft(item);
                const projectionElement = getProjectionElement(item);
                const electronConfig = getNeutralElectronConfigSummary(projectionElement);
                return (
                  <div key={item.id} className="param-grid" style={{ marginBottom: "1rem" }}>
                    <div className="param-row">
                      <label>
                        Target
                        <InfoTooltip text="Element targets apply the chosen template to every matching atom. Site targets apply it to one specific atomic site." />
                      </label>
                      <select
                        value={item.targetType}
                        onChange={(event) => updateProjection(item.id, {
                          targetType: event.target.value === "site" ? "site" : "element",
                          siteIndex: event.target.value === "site" ? 0 : null,
                        })}
                      >
                        <option value="element">Element</option>
                        <option value="site">Site</option>
                      </select>
                    </div>
                    {item.targetType === "element" ? (
                      <div className="param-row">
                        <label>Element</label>
                        <select
                          value={item.symbol}
                          onChange={(event) => updateProjection(item.id, { symbol: event.target.value })}
                        >
                          {uniqueCrystalElements.map((element) => (
                            <option key={element} value={element}>{element}</option>
                          ))}
                        </select>
                      </div>
                    ) : (
                      <div className="param-row">
                        <label>Site</label>
                        <select
                          value={item.siteIndex ?? 0}
                          onChange={(event) => updateProjection(item.id, { siteIndex: Number.parseInt(event.target.value, 10) })}
                        >
                          {crystalData.atom_sites.map((site, index) => (
                            <option key={`${site.label}-${index}`} value={index}>
                              {site.label} ({getBaseElement(site.type_symbol)}) [{site.fract_x.toFixed(3)}, {site.fract_y.toFixed(3)}, {site.fract_z.toFixed(3)}]
                            </option>
                          ))}
                        </select>
                      </div>
                    )}
                    <div className="param-row">
                      <label>
                        Orbital template
                        <InfoTooltip text="QCortado leaves all Wannier projection templates available. Use the neutral-atom reference below as a reminder of the element's baseline shell filling, not as a hard rule for transport-ready choices." />
                      </label>
                      <select
                        value={item.orbital}
                        onChange={(event) => updateProjection(item.id, { orbital: event.target.value as ProjectionOrbital })}
                      >
                        {allowedOrbitalOptions.map((orbital) => (
                          <option key={orbital} value={orbital}>{orbital}</option>
                        ))}
                      </select>
                      <span className="param-hint">
                        Neutral atom reference for {projectionElement}:{" "}
                        {electronConfig ? renderElectronConfiguration(electronConfig.compact) : "unavailable"}
                      </span>
                    </div>
                    <div className="param-row">
                      <label>Contribution</label>
                      <span>{projectionContributionCount(item, wannierCellAtomSymbols)} Wannier functions</span>
                    </div>
                    <div className="param-row">
                      <label>&nbsp;</label>
                      <button className="delete-calc-btn" type="button" onClick={() => removeProjection(item.id)}>
                        Remove
                      </button>
                    </div>
                  </div>
                );
              })}
            </div>

            <div className="param-section">
              <h4>Wannier Subspace</h4>
              <div className="param-grid">
                <div className="param-row">
                  <label>
                    Total Wannier functions
                    <InfoTooltip text="Wannier90 input: `num_wann`. QCortado derives this from the selected projections. Override it only if you are deliberately matching that same total by hand." />
                  </label>
                  <input
                    type="number"
                    min={1}
                    value={numWannOverrideInput}
                    placeholder={`${derivedNumWann}`}
                    onChange={(event) => setNumWannOverrideInput(event.target.value)}
                  />
                </div>
                <div className="param-row">
                  <label>
                    Bands included in Wannierization
                    <InfoTooltip text="Wannier90 input: `num_bands`. Increase this above `num_wann` when the target bands are entangled. This enlarges the space that Wannier90 may disentangle before constructing the final subspace." />
                  </label>
                  <div style={{ display: "flex", alignItems: "center", gap: "0.75rem", flexWrap: "wrap" }}>
                    <input
                      type="number"
                      min={1}
                      value={effectiveNumBandsInput}
                      onChange={(event) => {
                        setNumBandsMode("manual");
                        setNumBandsInput(event.target.value);
                      }}
                    />
                    <label className="toggle-label" style={{ minWidth: "auto", fontSize: "0.85rem" }}>
                      <input
                        type="checkbox"
                        checked={numBandsMode === "auto"}
                        onChange={(event) => {
                          if (event.target.checked) {
                            setNumBandsMode("auto");
                            return;
                          }
                          setNumBandsInput(effectiveNumBandsInput);
                          setNumBandsMode("manual");
                        }}
                      />
                      Auto
                    </label>
                  </div>
                </div>
                <div className="param-row" style={{ alignItems: "flex-start" }}>
                  <label>Auto-fill logic</label>
                  <span className="help-text" style={{ margin: 0, maxWidth: "720px" }}>
                    {numBandsHelpText}
                  </span>
                </div>
              </div>
            </div>

            <div className="param-section">
              <h4>
                Disentanglement Controls
                <InfoTooltip text="Wannier90 inputs: `exclude_bands`, `dis_win_min`, `dis_win_max`, `dis_froz_min`, and `dis_froz_max`. Use the outer window to include all bands that may mix into the target space, and the frozen window to enforce an exact match for the bands you care about most, usually near E_F." />
              </h4>
              <div className="param-grid">
                <div className="param-row">
                  <label>
                    Excluded source bands
                    <InfoTooltip text="Wannier90 input: `exclude_bands`. Exclude low-lying semicore or otherwise irrelevant bands from the disentanglement pool. Use QE/Wannier band plots to identify them first." />
                  </label>
                  <input
                    type="text"
                    value={excludeBandsInput}
                    placeholder="e.g. 1-4, 12"
                    onChange={(event) => setExcludeBandsInput(event.target.value)}
                  />
                </div>
                <div className="param-row">
                  <label>
                    Outer window minimum (eV)
                    <InfoTooltip text="Wannier90 input: `dis_win_min`. Lower edge of the outer disentanglement window in absolute eV." />
                  </label>
                  <input
                    value={disWinMinInput}
                    onChange={(event) => setDisWinMinInput(event.target.value)}
                    onBlur={commitDisentanglementInputRanges}
                    onKeyDown={handleDisentanglementInputKeyDown}
                  />
                </div>
                <div className="param-row">
                  <label>
                    Outer window maximum (eV)
                    <InfoTooltip text="Wannier90 input: `dis_win_max`. Upper edge of the outer disentanglement window in absolute eV." />
                  </label>
                  <input
                    value={disWinMaxInput}
                    onChange={(event) => setDisWinMaxInput(event.target.value)}
                    onBlur={commitDisentanglementInputRanges}
                    onKeyDown={handleDisentanglementInputKeyDown}
                  />
                </div>
                <div className="param-row">
                  <label>
                    Frozen window minimum (eV)
                    <InfoTooltip text="Wannier90 input: `dis_froz_min`. Lower edge of the frozen window in absolute eV." />
                  </label>
                  <input
                    value={disFrozMinInput}
                    onChange={(event) => setDisFrozMinInput(event.target.value)}
                    onBlur={commitDisentanglementInputRanges}
                    onKeyDown={handleDisentanglementInputKeyDown}
                  />
                </div>
                <div className="param-row">
                  <label>
                    Frozen window maximum (eV)
                    <InfoTooltip text="Wannier90 input: `dis_froz_max`. Upper edge of the frozen window in absolute eV." />
                  </label>
                  <input
                    value={disFrozMaxInput}
                    onChange={(event) => setDisFrozMaxInput(event.target.value)}
                    onBlur={commitDisentanglementInputRanges}
                    onKeyDown={handleDisentanglementInputKeyDown}
                  />
                </div>
              </div>
              <p className="help-text" style={{ marginBottom: 0 }}>
                Drag the cyan/orange windows in the reference plot or use left/right side sliders; numeric fields stay synchronized in absolute eV.
              </p>
            </div>

            <div className="calculation-summary">
              <h4>Summary</h4>
              <div className="summary-row">
                <span>Source SCF:</span>
                <span>{formatCalculationSourceLabel(selectedScf)}</span>
              </div>
              <div className="summary-row">
                <span>NSCF mesh:</span>
                <span>{kGridInput.join("×")}</span>
              </div>
              <div className="summary-row">
                <span>num_wann:</span>
                <span>{resolvedNumWann}</span>
              </div>
              <div className="summary-row">
                <span>num_bands:</span>
                <span>{effectiveNumBandsInput || "N/A"}{numBandsMode === "auto" ? " (auto)" : ""}</span>
              </div>
              <div className="summary-row">
                <span>Projections:</span>
                <span>{projectionSummary.length > 0 ? projectionSummary.join("; ") : "None"}</span>
              </div>
            </div>

            {isHpcMode ? (
              <HpcRunSettings
                profileId={activeHpcProfile?.id ?? null}
                profileName={activeHpcProfile?.name ?? "Andromeda"}
                taskKind="wannier"
                commandLines={hpcCommandLines}
                resources={hpcResources}
                onResourcesChange={setHpcResources}
                resourceMode={activeHpcProfile?.resource_mode ?? "both"}
                defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
                defaultGpuResources={activeHpcProfile?.default_gpu_resources ?? null}
                maxTasks={manualNumBandsForHpc != null ? { value: manualNumBandsForHpc, reason: "manual num_bands is active" } : null}
                disabled={isRunning}
              />
            ) : (
              <div className="mpi-section">
                <h4>Parallelization</h4>
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
                            onChange={(e) => setMpiProcs(Math.max(1, Number.parseInt(e.target.value, 10) || 1))}
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

            {error && <div className="error-message">{error}</div>}
            {parameterValidationError && <div className="error-message">{parameterValidationError}</div>}

            <div className="step-actions step-actions-sticky">
              <button className="secondary-button" onClick={() => setStep("mesh")}>
                Back
              </button>
              <button
                className="primary-button"
                onClick={() => void runCalculation()}
                disabled={Boolean(parameterValidationError) || hasBlockingExternalTask}
              >
                {isHpcMode ? "Submit Wannier to Andromeda" : "Run Calculation"}
              </button>
            </div>
          </div>

          <div className="wannier-projection-plot-pane">
            <div className="source-step-header">
              <h4>
                Reference Bands Plotter
                <InfoTooltip text="Load a saved bands calculation and interactively drag disentanglement/frozen windows." />
              </h4>
            </div>

            {referenceBandsCalculations.length === 0 ? (
              <div className="info-banner">
                {referenceBandsLoading
                  ? "Loading saved bands calculations..."
                  : "No saved bands calculations with plot data are available yet. Run a bands calculation first to use this plotter."}
              </div>
            ) : (
              <>
                {referenceBandsLoading && (
                  <div className="info-banner">
                    Loading additional saved bands calculations...
                  </div>
                )}
                <div className="param-grid">
                  <div className="param-row">
                    <label>
                      Bands calculation
                      <InfoTooltip text="Choose which saved bands result to preview for window tuning." />
                    </label>
                    <select
                      value={selectedReferenceBandsCalculation?.id ?? ""}
                      onChange={(event) => setReferenceBandsCalcId(event.target.value)}
                    >
                      {referenceBandsCalculations.map((calc) => (
                        <option key={calc.id} value={calc.id}>
                          {formatReferenceBandsLabel(calc)}
                        </option>
                      ))}
                    </select>
                  </div>
                  <div className="param-row">
                    <label>
                      Energy mode
                      <InfoTooltip text="Toggle between E − E_F and absolute eV display. `dis_*` values remain absolute eV." />
                    </label>
                    <select
                      value={referenceEnergyMode}
                      onChange={(event) => setReferenceEnergyMode(event.target.value as ReferenceEnergyMode)}
                    >
                      <option value="relative">E − E_F (eV)</option>
                      <option value="absolute">Absolute energy (eV)</option>
                    </select>
                  </div>
                  <div className="param-row">
                    <label>
                      Fat-band projection
                      <InfoTooltip text="Select the saved projection group to render as fat bands. Disabled if the selected run has no projection data." />
                    </label>
                    <select
                      value={referenceProjectionSelection}
                      onChange={(event) => setReferenceProjectionSelection(event.target.value)}
                      disabled={referenceProjectionOptions.length <= 1}
                    >
                      {referenceProjectionOptions.map((option) => (
                        <option key={option.value} value={option.value}>
                          {option.label}
                        </option>
                      ))}
                    </select>
                  </div>
                  <div className="param-row">
                    <label>
                      Band-gap overlay
                      <InfoTooltip text="Show or hide the detected band-gap shading and labels on the plot." />
                    </label>
                    <label className="toggle-label" style={{ minWidth: "auto", fontSize: "0.85rem" }}>
                      <input
                        type="checkbox"
                        checked={referenceShowBandGapOverlay}
                        onChange={(event) => setReferenceShowBandGapOverlay(event.target.checked)}
                      />
                      Show
                    </label>
                  </div>
                </div>
                {selectedReferenceBandData && (
                  <div className="wannier-reference-bands-plot wannier-reference-bands-plot-full">
                    <BandPlot
                      data={selectedReferenceBandData}
                      scfFermiEnergy={referencePlotFermi}
                      projectionSelection={referenceProjectionSelection}
                      viewerType="electronic"
                      yAxisLabel={referenceEnergyMode === "absolute" ? "Energy (eV)" : "E − E_F (eV)"}
                      valueLabel={referenceEnergyMode === "absolute" ? "Energy" : "E − E_F"}
                      showFermiLevel={referenceEnergyMode !== "absolute"}
                      showSidebar={false}
                      enableWheelRangeControl
                      enableHoverScrollLock
                      yClampRange={null}
                      height={640}
                      windowOverlays={referenceWindowOverlays}
                      onWindowOverlayChange={handleReferenceWindowOverlayChange}
                      windowOverlayHint="Drag window edges or side sliders for fine absolute-eV tuning."
                      showBandGapOverlayOverride={referenceShowBandGapOverlay}
                      scrollHint="Scroll: zoom Y | Shift+Scroll: pan energy"
                    />
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      </div>
    );
  };

  const renderRunStep = () => {
    return (
      <div className="wizard-step run-step run-step-focused">
        <div className="run-step-headline">
          <h3>{isRunning ? "Running Wannier90 Workflow" : "Wannier90 Output"}</h3>
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
              <button className="secondary-button" onClick={() => setStep("projections")}>
                Back to Parameters
              </button>
            </div>
          </>
        )}

        <div className="run-layout run-layout-single">
          <LiveOutputPanel
            title={isRunning ? "Running..." : "Output"}
            output={output}
            placeholder="Starting calculation..."
            totalLineCount={outputLineCount}
            visibleLineCount={visibleOutputLineCount}
          />
        </div>
      </div>
    );
  };

  const renderResultsStep = () => {
    if (!result) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("projections")}>
            Back to Parameters
          </button>
        </div>
      );
    }

    return (
      <div className="wizard-step results-step">
        <h3>Wannier90 Results</h3>
        <p className="step-description">
          Calculation complete. Open the main Wannier viewer for interactive band plotting.
        </p>

        {resultQualityIssues.length > 0 && (
          <div className="warning-banner">
            {resultQualityIssues.map((issue) => issue.message).join(" ")}
          </div>
        )}

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">num_wann:</span>
              <span className="value">{result.num_wann}</span>
            </div>
            <div className="summary-item">
              <span className="label">num_bands:</span>
              <span className="value">{result.num_bands}</span>
            </div>
            <div className="summary-item">
              <span className="label">K-points:</span>
              <span className="value">{result.band_data.n_kpoints}</span>
            </div>
            <div className="summary-item">
              <span className="label">Total Spread:</span>
              <span className="value">
                {result.total_spread != null ? `${result.total_spread.toFixed(6)} A^2` : "N/A"}
              </span>
            </div>
            <div className="summary-item">
              <span className="label">Converged:</span>
              <span className="value">{result.convergence?.converged ? "Yes" : "No"}</span>
            </div>
            <div className="summary-item">
              <span className="label">Iterations:</span>
              <span className="value">{result.convergence?.iterations ?? "N/A"}</span>
            </div>
            <div className="summary-item">
              <span className="label">Minimization:</span>
              <span className="value">{formatWannierConvergenceFlag(result.convergence?.minimization_converged)}</span>
            </div>
            <div className="summary-item">
              <span className="label">Disentanglement:</span>
              <span className="value">{formatWannierConvergenceFlag(result.convergence?.disentanglement_converged)}</span>
            </div>
          </div>
        </div>

        {(result.convergence?.failure_reasons?.length || result.convergence?.warnings?.length) ? (
          <div className="detail-item parameters">
            <label>Quality Checks</label>
            <pre>
              {[
                ...(result.convergence?.failure_reasons || []).map((entry) => `Error: ${entry}`),
                ...(result.convergence?.warnings || []).map((entry) => `Warning: ${entry}`),
              ].join("\n")}
            </pre>
          </div>
        ) : null}

        {result.spreads.length > 0 && (
          <div className="detail-item parameters">
            <label>Spreads</label>
            <pre>
              {result.spreads.map((entry) =>
                `WF ${entry.index}: center=(${entry.centre.map((value) => value.toFixed(5)).join(", ")}) spread=${entry.spread.toFixed(6)}`
              ).join("\n")}
            </pre>
          </div>
        )}

        {result.artifact_manifest.length > 0 && (
          <div className="detail-item parameters">
            <label>Artifacts</label>
            <pre>
              {result.artifact_manifest.map((entry) => `${entry.file_name} (${formatBytes(entry.size_bytes)})`).join("\n")}
            </pre>
          </div>
        )}

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
            onClick={() => {
              void (async () => {
                const overlays = await buildViewerOverlayOptions();
                onViewWannier(result, selectedScf?.result?.fermi_energy ?? null, overlays);
              })();
            }}
          >
            View Wannier
          </button>
          <button className="primary-button" onClick={() => setStep("projections")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className={`wannier-wizard wizard-step-${step}`}>
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Wannier90 Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "mesh" ? "active" : ["projections", "run", "results"].includes(step) ? "completed" : ""}>
            2. K-Mesh
          </span>
          <span className={step === "projections" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
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
        {step === "source" && renderSourceStep()}
        {step === "mesh" && renderMeshStep()}
        {step === "projections" && renderProjectionStep()}
        {step === "run" && renderRunStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
