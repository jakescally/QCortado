import { useCallback, useEffect, useMemo, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import {
  CrystalData,
  ELEMENT_MASSES,
  ExecutionMode,
  HpcProfile,
  SlurmResourceRequest,
} from "../lib/types";
import { BandData } from "./BandPlot";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";
import {
  analyzeCrystalSymmetry,
  buildConventionalLatticeFromCrystalData,
  SymmetryTransformResult,
} from "../lib/symmetryTransform";
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
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { useTaskContext } from "../lib/TaskContext";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import {
  buildExecutionTarget,
  buildHpcLauncherCommand,
  buildHpcQeInputCommandLine,
  defaultResourcesForProfile,
} from "../lib/hpcConfig";
import { HpcRunSettings } from "./HpcRunSettings";

interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    fermi_energy: number | null;
    raw_output?: string | null;
    wannier_data?: WannierResult | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
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
  omega_i?: number | null;
  omega_d?: number | null;
  omega_od?: number | null;
  omega_total?: number | null;
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
  artifact_manifest: WannierArtifact[];
}

type ProjectionTargetType = "element" | "site";
type ProjectionOrbital = "s" | "p" | "d" | "f" | "sp" | "sp2" | "sp3" | "sp3d" | "sp3d2";
type WizardStep = "source" | "mesh" | "projections" | "run" | "results";
type NumBandsMode = "auto" | "manual";

interface ProjectionDraft {
  id: string;
  targetType: ProjectionTargetType;
  symbol: string;
  orbital: ProjectionOrbital;
  siteIndex: number | null;
}

interface WannierWizardProps {
  onBack: () => void;
  onExecutionModeChange?: (mode: ExecutionMode) => Promise<void> | void;
  onViewWannier: (result: WannierResult, scfFermiEnergy: number | null) => void;
  qePath: string;
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

function normalizeSmearing(raw: unknown): "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" {
  const lowered = String(raw || "gaussian").toLowerCase();
  if (lowered === "methfessel-paxton") return "methfessel-paxton";
  if (lowered === "marzari-vanderbilt") return "marzari-vanderbilt";
  if (lowered === "fermi-dirac") return "fermi-dirac";
  return "gaussian";
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

function getWannierSourceIssues(calc: CalculationRun): string[] {
  const params = calc.parameters || {};
  const issues: string[] = [];
  if (!sourceScfUsesPrimitiveCell(params)) {
    issues.push("Primitive-cell SCF required");
  }
  if (Number(params.nspin) !== 1) {
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

export function WannierWizard({
  onBack,
  onViewWannier,
  qePath,
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: WannierWizardProps) {
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
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
  const [kGridInput, setKGridInput] = useState<[string, string, string]>(["4", "4", "4"]);
  const [bandsNumPointsInput, setBandsNumPointsInput] = useState("100");
  const [numBandsInput, setNumBandsInput] = useState("");
  const [numBandsMode, setNumBandsMode] = useState<NumBandsMode>("auto");
  const [seednameInput, setSeednameInput] = useState("qcortado_wannier");
  const [kPath, setKPath] = useState<KPathPoint[]>([]);
  const [kPathRhombohedralConvention, setKPathRhombohedralConvention] = useState<RhombohedralConvention | undefined>(undefined);
  const [projectionDrafts, setProjectionDrafts] = useState<ProjectionDraft[]>([]);
  const [numWannOverrideInput, setNumWannOverrideInput] = useState("");
  const [excludeBandsInput, setExcludeBandsInput] = useState("");
  const [disWinMinInput, setDisWinMinInput] = useState("");
  const [disWinMaxInput, setDisWinMaxInput] = useState("");
  const [disFrozMinInput, setDisFrozMinInput] = useState("");
  const [disFrozMaxInput, setDisFrozMaxInput] = useState("");
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
    const storedGrid = Array.isArray(params.kgrid) && params.kgrid.length === 3
      ? params.kgrid.map((value: unknown) => String(value || "4")) as [string, string, string]
      : ["4", "4", "4"];
    setKGridInput([storedGrid[0], storedGrid[1], storedGrid[2]]);

    if (projectionDrafts.length === 0) {
      const uniqueElements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
      if (uniqueElements.length === 1) {
        setProjectionDrafts([
          {
            id: makeProjectionId(),
            targetType: "element",
            symbol: uniqueElements[0],
            orbital: "sp3",
            siteIndex: null,
          },
        ]);
      }
    }
  }, [selectedScf, crystalData, projectionDrafts.length]);

  useEffect(() => {
    setNumBandsMode("auto");
    setNumBandsInput("");
  }, [selectedScf?.id]);

  useEffect(() => {
    if (!selectedScf) return;
    const params = selectedScf.parameters || {};
    if (!sourceScfUsesPrimitiveCell(params)) return;
    if (symmetryTransform || symmetryError) return;
    void ensureSymmetryTransform();
  }, [selectedScf, symmetryError, symmetryTransform]);

  const wannierCellAtomSymbols = useMemo(() => {
    const params = selectedScf?.parameters || {};
    if (
      sourceScfUsesPrimitiveCell(params)
      && symmetryTransform
      && symmetryTransform.standardizedPrimitiveAtoms.length > 0
    ) {
      return symmetryTransform.standardizedPrimitiveAtoms.map((atom) => getBaseElement(atom.symbol));
    }
    return crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol));
  }, [crystalData, selectedScf, symmetryTransform]);

  const validScfs = useMemo(
    () => scfCalculations.filter((calc) => calc.calc_type === "scf" && calc.result?.converged),
    [scfCalculations],
  );
  const sortedScfs = useMemo(
    () => sortScfByMode(validScfs, scfSortMode),
    [validScfs, scfSortMode],
  );
  const readyScfs = useMemo(
    () => validScfs.filter((calc) => getWannierSourceIssues(calc).length === 0),
    [validScfs],
  );

  const selectedIssues = useMemo(
    () => (selectedScf ? getWannierSourceIssues(selectedScf) : []),
    [selectedScf],
  );

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
  const totalInterpolatedPoints = useMemo(
    () => kPath.reduce((sum, point) => sum + (Number(point.npoints) || 0), 0),
    [kPath],
  );
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
      const bandsNumPoints = parseOptionalPositiveInt(bandsNumPointsInput, "bands_num_points");
      const numBands = parseOptionalPositiveInt(effectiveNumBandsInput, "num_bands");
      if (kx == null || ky == null || kz == null || bandsNumPoints == null || numBands == null) {
        return "Use positive integers for the mesh, bands_num_points, and num_bands.";
      }
      if (minOccupiedBands != null && numBands < minOccupiedBands) {
        return `num_bands must be at least ${minOccupiedBands} for this scalar source (${sourceElectronCount} electrons => ${minOccupiedBands} occupied bands).`;
      }
      if (numBands < resolvedNumWann) {
        return "num_bands must be greater than or equal to num_wann.";
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
    bandsNumPointsInput,
    derivedNumWann,
    disFrozMaxInput,
    disFrozMinInput,
    disWinMaxInput,
    disWinMinInput,
    kGridInput,
    kPath,
    effectiveNumBandsInput,
    projectionDrafts.length,
    resolvedNumWann,
    selectedIssues,
    selectedScf,
    minOccupiedBands,
    sourceElectronCount,
  ]);

  const hpcCommandLines = useMemo(() => {
    const launcher = buildHpcLauncherCommand(activeHpcProfile);
    const remoteWannier = (activeHpcProfile?.remote_wannier90_path || "wannier90.x").trim() || "wannier90.x";
    const seedname = sanitizeSeedname(seednameInput);
    return [
      `${launcher} "${remoteWannier}" -pp ${seedname} > wannier90_pre.out 2>&1`,
      buildHpcQeInputCommandLine(activeHpcProfile, "pw.x", "nscf.in", "nscf.out"),
      buildHpcQeInputCommandLine(activeHpcProfile, "pw2wannier90.x", "pw2wan.in", "pw2wan.out"),
      `${launcher} "${remoteWannier}" ${seedname} > wannier90.out 2>&1`,
    ];
  }, [activeHpcProfile, seednameInput]);

  const sourcePseudoMap = useMemo<Record<string, string>>(() => {
    const params = selectedScf?.parameters || {};
    return (params.selected_pseudos && typeof params.selected_pseudos === "object")
      ? params.selected_pseudos as Record<string, string>
      : {};
  }, [selectedScf]);

  const addProjection = useCallback(() => {
    const firstElement = getBaseElement(crystalData.atom_sites[0]?.type_symbol || "X");
    setProjectionDrafts((prev) => [
      ...prev,
      {
        id: makeProjectionId(),
        targetType: "element",
        symbol: firstElement,
        orbital: "s",
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
    const bandsNumPoints = parseOptionalPositiveInt(bandsNumPointsInput, "bands_num_points") ?? 0;
    const numBands = parseOptionalPositiveInt(effectiveNumBandsInput, "num_bands") ?? 0;
    const seedname = sanitizeSeedname(seednameInput);
    if (minOccupiedBands != null && numBands < minOccupiedBands) {
      throw new Error(
        `num_bands must be at least ${minOccupiedBands} for this scalar source (${sourceElectronCount} electrons => ${minOccupiedBands} occupied bands).`,
      );
    }
    if (numBands < resolvedNumWann) {
      throw new Error("num_bands must be greater than or equal to num_wann.");
    }

    const params = selectedScf.parameters || {};
    const ecutwfc = Number(params.ecutwfc);
    const ecutrho = Number(params.ecutrho);
    const nspin = Number(params.nspin) || 1;
    const noncolin = Boolean(params.noncolin);
    const lspinorb = Boolean(params.lspinorb);
    const occupations = normalizeOccupations(params.occupations);
    const smearing = normalizeSmearing(params.smearing);
    const degauss = Number(params.degauss);
    const convThr = Number(params.conv_thr);
    const mixingBeta = Number(params.mixing_beta);
    const baseElements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
    const species = baseElements.map((element) => ({
      symbol: element,
      mass: ELEMENT_MASSES[element] || 1,
      pseudopotential: sourcePseudoMap[element],
    }));
    if (species.some((entry) => !entry.pseudopotential)) {
      throw new Error("Source SCF pseudopotential metadata is incomplete. Re-run the SCF before starting Wannier.");
    }

    const resolvedSymmetry = await ensureSymmetryTransform();
    const sourceUsesPrimitive = sourceScfUsesPrimitiveCell(params);
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

    const baseCalculation = canUseSymmetryPrimitive && resolvedSymmetry
      ? {
        calculation: "scf",
        prefix: params.prefix || "qcortado_scf",
        outdir: "./tmp",
        pseudo_dir: isHpcMode ? (activeHpcProfile?.remote_pseudo_dir || "") : qePath.replace(/\/bin\/?$/, "/pseudo"),
        system: {
          ibrav: "free",
          celldm: null,
          cell_parameters: resolvedSymmetry.standardizedPrimitiveLattice,
          cell_units: "angstrom",
          species,
          atoms: resolvedSymmetry.standardizedPrimitiveAtoms.map((atom) => ({
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
      }
      : {
        calculation: "scf",
        prefix: params.prefix || "qcortado_scf",
        outdir: "./tmp",
        pseudo_dir: isHpcMode ? (activeHpcProfile?.remote_pseudo_dir || "") : qePath.replace(/\/bin\/?$/, "/pseudo"),
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

    const transformedPath = mapPathCoordinates(
      kPath,
      canUseSymmetryPrimitive && resolvedSymmetry
        ? converters.toSymmetryPrimitiveCoords
        : converters.toInputConventionalCoords,
    ).map((point) => ({
      label: point.label,
      coords: point.coords,
      npoints: point.npoints,
    }));

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
            bands_num_points: bandsNumPoints,
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
        bands_num_points: bandsNumPoints,
        scf_fermi_energy: selectedScf.result?.fermi_energy ?? null,
        cell_representation: canUseSymmetryPrimitive ? "primitive_spglib" : "conventional_input",
        ecutwfc: params.ecutwfc,
        nspin: params.nspin,
        lspinorb: params.lspinorb,
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
      const taskId = await taskContext.startTask("wannier", plan.taskParams, plan.taskLabel);
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
          Choose the converged scalar SCF that will supply the restart density and primitive-cell metadata for Wannierization.
        </p>

        {readyScfs.length === 0 && (
          <div className="warning-banner">
            No Wannier-ready SCFs found. Wannier v1 requires a primitive-cell scalar source with no SOC, noncollinearity, DFT+U, or certified vdW corrections.
          </div>
        )}

        <div className="scf-list">
          {sortedScfs.map((scf) => {
            const issues = getWannierSourceIssues(scf);
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
                  {issues.map((issue) => (
                    <span key={issue} className="calc-tag calc-tag-feature">
                      {issue}
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
            onClick={() => setStep("mesh")}
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
              <label>bands_num_points</label>
              <input
                type="number"
                min={1}
                value={bandsNumPointsInput}
                onChange={(event) => setBandsNumPointsInput(event.target.value)}
              />
            </div>
            <div className="param-row">
              <label>seedname</label>
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
          onPathChange={setKPath}
          initialPath={kPath}
          pointsPerSegment={Math.max(1, Number.parseInt(bandsNumPointsInput, 10) || 100)}
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
            <span>{selectedScf?.id.slice(0, 8) || "N/A"}</span>
          </div>
          <div className="summary-row">
            <span>NSCF mesh:</span>
            <span>{kGridInput.join("×")}</span>
          </div>
          <div className="summary-row">
            <span>bands_num_points:</span>
            <span>{bandsNumPointsInput || "100"}</span>
          </div>
          <div className="summary-row">
            <span>Path nodes:</span>
            <span>{kPath.map((point) => point.label).join(" → ") || "Not selected"}</span>
          </div>
          <div className="summary-row">
            <span>Interpolated k-points:</span>
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
      <div className="wizard-step parameters-step">
        <h3>Projections and Disentanglement</h3>
        <p className="step-description">
          Define the initial orbital projections and any disentanglement windows used when `num_bands` exceeds `num_wann`.
        </p>
        <div className="info-banner">
          Element projections apply to every matching atom in the Wannier cell. Use site-targeted projections if you want only one atom.
        </div>

        <div className="param-section">
          <div className="source-step-header">
            <h4>Initial Projections</h4>
            <button className="secondary-button" type="button" onClick={addProjection}>
              Add Projection
            </button>
          </div>

          {projectionDrafts.length === 0 && (
            <div className="warning-banner">
              Add at least one element- or site-targeted orbital template.
            </div>
          )}

          {projectionDrafts.map((item) => (
            <div key={item.id} className="param-grid" style={{ marginBottom: "1rem" }}>
              <div className="param-row">
                <label>Target</label>
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
                    {[...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))].map((element) => (
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
                <label>Orbital template</label>
                <select
                  value={item.orbital}
                  onChange={(event) => updateProjection(item.id, { orbital: event.target.value as ProjectionOrbital })}
                >
                  {ORBITAL_OPTIONS.map((orbital) => (
                    <option key={orbital} value={orbital}>{orbital}</option>
                  ))}
                </select>
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
          ))}
        </div>

        <div className="param-section">
          <h4>Wannier Subspace</h4>
          <div className="param-grid">
            <div className="param-row">
              <label>num_wann override</label>
              <input
                type="number"
                min={1}
                value={numWannOverrideInput}
                placeholder={`${derivedNumWann}`}
                onChange={(event) => setNumWannOverrideInput(event.target.value)}
              />
            </div>
            <div className="param-row">
              <label>num_bands</label>
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
          <h4>Disentanglement Controls</h4>
          <div className="param-grid">
            <div className="param-row">
              <label>exclude_bands</label>
              <input
                type="text"
                value={excludeBandsInput}
                placeholder="e.g. 1-4, 12"
                onChange={(event) => setExcludeBandsInput(event.target.value)}
              />
            </div>
            <div className="param-row">
              <label>dis_win_min (eV)</label>
              <input value={disWinMinInput} onChange={(event) => setDisWinMinInput(event.target.value)} />
            </div>
            <div className="param-row">
              <label>dis_win_max (eV)</label>
              <input value={disWinMaxInput} onChange={(event) => setDisWinMaxInput(event.target.value)} />
            </div>
            <div className="param-row">
              <label>dis_froz_min (eV)</label>
              <input value={disFrozMinInput} onChange={(event) => setDisFrozMinInput(event.target.value)} />
            </div>
            <div className="param-row">
              <label>dis_froz_max (eV)</label>
              <input value={disFrozMaxInput} onChange={(event) => setDisFrozMaxInput(event.target.value)} />
            </div>
          </div>
        </div>

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>Source SCF:</span>
            <span>{selectedScf?.id.slice(0, 8) || "N/A"}</span>
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

        <div className="step-actions">
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

        <div className="run-layout">
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
          </div>
        </div>

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
            onClick={() => onViewWannier(result, selectedScf?.result?.fermi_energy ?? null)}
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
