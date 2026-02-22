// Phonon Calculation Wizard - Calculate phonon DOS and dispersion

import { useState, useEffect, useRef, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData } from "../lib/types";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { EstimatedRemainingTime } from "./EstimatedRemainingTime";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { isPhononReadyScf } from "../lib/phononReady";
import { analyzeCrystalSymmetry, SymmetryTransformResult } from "../lib/symmetryTransform";
import {
  createPathCoordinateConverters,
  mapPathCoordinates,
  resolvePathTransformContext,
  sourceScfUsesPrimitiveCell,
} from "../lib/kPathTransforms";
import {
  RhombohedralConvention,
  defaultRhombohedralConventionForSetting,
} from "../lib/brillouinZoneData";

function Tooltip({ text }: { text: string }) {
  return (
    <span className="tooltip-container">
      <span className="tooltip-icon">?</span>
      <span className="tooltip-text">{text}</span>
    </span>
  );
}

function parseOptionalPositiveNumber(value: string): number | null {
  const trimmed = value.trim();
  if (!trimmed) return null;
  const parsed = Number(trimmed);
  return Number.isFinite(parsed) && parsed > 0 ? parsed : null;
}

function parseOptionalPositiveInt(value: string): number | null {
  const trimmed = value.trim();
  if (!trimmed) return null;
  const parsed = Number(trimmed);
  return Number.isInteger(parsed) && parsed > 0 ? parsed : null;
}

function isOptionalPositiveNumberValid(value: string): boolean {
  return value.trim() === "" || parseOptionalPositiveNumber(value) !== null;
}

function isOptionalPositiveIntValid(value: string): boolean {
  return value.trim() === "" || parseOptionalPositiveInt(value) !== null;
}

type QPathSamplingMode = "segment" | "total";

const MIN_POINTS_PER_SEGMENT = 5;
const MAX_POINTS_PER_SEGMENT = 400;
const MAX_TOTAL_Q_POINTS = 5000;

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

function applyPointsPerSegment(path: KPathPoint[], pointsPerSegment: number): KPathPoint[] {
  const connectedSegmentIndices = new Set(getConnectedSegmentIndices(path));
  if (path.length === 0 || connectedSegmentIndices.size === 0) {
    return path.map((point) => ({
      ...point,
      npoints: 0,
    }));
  }

  const safePointsPerSegment = clampInt(
    pointsPerSegment,
    MIN_POINTS_PER_SEGMENT,
    MAX_POINTS_PER_SEGMENT,
  );
  return path.map((point, index) => {
    if (index >= path.length - 1) {
      return { ...point, npoints: 0 };
    }
    if (!connectedSegmentIndices.has(index)) {
      return { ...point, npoints: 0 };
    }
    return { ...point, npoints: safePointsPerSegment };
  });
}

function applyTotalQPoints(path: KPathPoint[], totalQPoints: number): KPathPoint[] {
  const connectedSegmentIndices = getConnectedSegmentIndices(path);
  if (path.length === 0 || connectedSegmentIndices.length === 0) {
    return path.map((point) => ({
      ...point,
      npoints: 0,
    }));
  }

  const safeTotal = clampInt(
    totalQPoints,
    connectedSegmentIndices.length,
    MAX_TOTAL_Q_POINTS,
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

function normalizeQPathSampling(
  path: KPathPoint[],
  samplingMode: QPathSamplingMode,
  pointsPerSegment: number,
  totalQPoints: number,
): KPathPoint[] {
  if (samplingMode === "total") {
    return applyTotalQPoints(path, totalQPoints);
  }
  return applyPointsPerSegment(path, pointsPerSegment);
}

interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    total_energy: number | null;
    fermi_energy: number | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

// Helper to generate calculation feature tags from parameters
function getCalculationTags(calc: CalculationRun): { label: string; type: "info" | "feature" | "special" | "geometry" }[] {
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

  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }

  if (params.conv_thr) {
    const thr = params.conv_thr;
    const label = thr < 0.001 ? thr.toExponential(0) : thr.toString();
    pushTag(label, "info");
  }

  // Mark phonon-ready calculations.
  if (phononReady) {
    pushTag("Phonon-Ready", "special");
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

  return tags;
}

interface PhononDOS {
  frequencies: number[];
  dos: number[];
  omega_max: number;
  omega_min: number;
}

interface PhononHighSymmetryMarker {
  q_distance: number;
  label: string;
}

interface PhononDispersion {
  q_points: number[];
  frequencies: number[][];
  high_symmetry_points: PhononHighSymmetryMarker[];
  n_modes: number;
  n_qpoints: number;
  frequency_range: [number, number];
}

interface PhononResult {
  converged: boolean;
  n_qpoints: number;
  n_modes: number;
  dos_data: PhononDOS | null;
  dispersion_data: PhononDispersion | null;
  raw_output: string;
}

interface PhononViewerData {
  dos_data: PhononDOS | null;
  dispersion_data: PhononDispersion | null;
}

interface PhononWizardProps {
  onBack: () => void;
  onViewPhonons: (phononData: PhononViewerData, viewMode: "bands" | "dos") => void;
  qePath: string;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "qgrid" | "options" | "run" | "results";
const PHONON_WORK_DIR = "/tmp/qcortado_phonon";

interface PhononTaskPlan {
  taskLabel: string;
  taskParams: Record<string, any>;
  saveParameters: Record<string, any>;
  progressMeta: {
    hasDos: boolean;
    hasDispersion: boolean;
  };
}

export function PhononWizard({
  onBack,
  onViewPhonons,
  qePath: _qePath,
  projectId,
  cifId: _cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: PhononWizardProps) {
  const taskContext = useTaskContext();
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  // Wizard state
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);

  // Step 1: Source SCF
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());
  const [symmetryTransform, setSymmetryTransform] = useState<SymmetryTransformResult | null>(null);
  const [symmetryError, setSymmetryError] = useState<string | null>(null);

  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);

  // Step 2: Q-Grid
  const [qGrid, setQGrid] = useState<[number, number, number]>([4, 4, 4]);
  const [tr2Ph, setTr2Ph] = useState<number>(1e-12);
  const [tr2PhInput, setTr2PhInput] = useState<string>("1e-12");
  const [asr, setAsr] = useState<"none" | "simple" | "crystal" | "one-dim" | "zero-dim">("crystal");
  const [recover, setRecover] = useState(false);
  const [trans, setTrans] = useState(true);
  const [epsil, setEpsil] = useState(false);
  const [fpol, setFpol] = useState(false);
  const [lraman, setLraman] = useState(false);
  const [verbosity, setVerbosity] = useState<"default" | "high" | "debug">("default");
  const [nmixPhInput, setNmixPhInput] = useState<string>("");
  const [niterPhInput, setNiterPhInput] = useState<string>("");
  const [alphaMixInput, setAlphaMixInput] = useState<string>("");
  const [startQInput, setStartQInput] = useState<string>("");
  const [lastQInput, setLastQInput] = useState<string>("");
  const [startIrrInput, setStartIrrInput] = useState<string>("");
  const [lastIrrInput, setLastIrrInput] = useState<string>("");
  const [electronPhononMode, setElectronPhononMode] = useState<"none" | "interpolated" | "lambda_tetra" | "gamma_tetra">("none");
  const [elPhSigmaInput, setElPhSigmaInput] = useState<string>("0.02");
  const [elPhNsigmaInput, setElPhNsigmaInput] = useState<string>("10");
  const [fildvscfEnabled, setFildvscfEnabled] = useState(false);
  const [fildvscf, setFildvscf] = useState("dvscf");
  const [epwPreparationEnabled, setEpwPreparationEnabled] = useState(false);
  const [preserveFullArtifacts, setPreserveFullArtifacts] = useState(true);

  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({
    solver: false,
    longRange: false,
    electronPhonon: false,
    epw: true,
    execution: false,
    dos: true,
    dispersion: true,
    mpi: false,
  });

  const toggleSection = (section: string) => {
    setExpandedSections((prev) => ({ ...prev, [section]: !prev[section] }));
  };

  // Step 3: Options
  const [calculateDos, setCalculateDos] = useState(true);
  const [dosGrid, setDosGrid] = useState<[number, number, number]>([20, 20, 20]);
  const [dosDeltaEInput, setDosDeltaEInput] = useState<string>("1.0");
  const [calculateDispersion, setCalculateDispersion] = useState(true);
  const [qPath, setQPath] = useState<KPathPoint[]>([]);
  const [qPathRhombohedralConvention, setQPathRhombohedralConvention] = useState<RhombohedralConvention | undefined>(undefined);
  const [qPathSamplingMode, setQPathSamplingMode] = useState<QPathSamplingMode>("segment");
  const [pointsPerSegment, setPointsPerSegment] = useState(20);
  const [totalQPointsTarget, setTotalQPointsTarget] = useState(120);
  const [totalQPointsInput, setTotalQPointsInput] = useState("120");

  // Step 4: Running
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);

  const handleOutputScroll = () => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followOutputRef.current = distanceToBottom <= 16;
  };
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Phonon calculation",
  });

  // Step 5: Results
  const [phononResult, setPhononResult] = useState<PhononResult | null>(null);
  const [showOutput, setShowOutput] = useState(false);
  const [isSaved, setIsSaved] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string>("");

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);

  // Load MPI info
  useEffect(() => {
    async function init() {
      try {
        const count = await invoke<number>("get_cpu_count");
        const safeCount = Math.max(1, Math.floor(count));
        setCpuCount(safeCount);
        const defaults = await loadGlobalMpiDefaults(safeCount);

        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);
        setMpiEnabled(available ? defaults.enabled : false);
        setMpiProcs(defaults.nprocs);
      } catch (e) {
        console.error("Failed to initialize:", e);
      }
    }
    init();
  }, []);

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
    setQPathRhombohedralConvention(undefined);
  }, [crystalData]);

  useEffect(() => {
    if (!epwPreparationEnabled) return;
    if (!fildvscfEnabled) setFildvscfEnabled(true);
    if (electronPhononMode === "none") {
      setElectronPhononMode("interpolated");
    }
  }, [epwPreparationEnabled, fildvscfEnabled, electronPhononMode]);

  // Auto-scroll output only if user is at the bottom
  useEffect(() => {
    const el = outputRef.current;
    if (!el || !followOutputRef.current) return;
    el.scrollTop = el.scrollHeight;
  }, [output]);

  // Reconnect to a running/completed background task
  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) {
      taskContext.reconnectToTask(activeTaskId);
      return;
    }

    setIsRunning(task.status === "running");
    setOutput(task.output.join("\n") + "\n");
    setProgress(task.progress);
    setCalcStartTime(task.startedAt);

    if (task.status === "completed" && task.result) {
      setPhononResult(task.result as any);
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
    if (!task || task.status !== "running") return;

    setOutput(task.output.join("\n") + "\n");
    setProgress(task.progress);
  }, [activeTaskId, taskContext.getTask(activeTaskId ?? "")?.output.length]);

  // Handle Q-path changes from the BZ viewer
  const handleQPathChange = useCallback((newPath: KPathPoint[]) => {
    setQPath(
      normalizeQPathSampling(newPath, qPathSamplingMode, pointsPerSegment, totalQPointsTarget),
    );
  }, [qPathSamplingMode, pointsPerSegment, totalQPointsTarget]);

  useEffect(() => {
    setQPath((prevPath) =>
      normalizeQPathSampling(prevPath, qPathSamplingMode, pointsPerSegment, totalQPointsTarget),
    );
  }, [qPathSamplingMode, pointsPerSegment, totalQPointsTarget]);

  const qPathSegmentCount = getConnectedSegmentIndices(qPath).length;
  const totalDispersionQPoints = qPath.reduce((sum, point) => sum + point.npoints, 0);
  const minimumTotalQPoints = Math.max(1, qPathSegmentCount);
  const viewerPointsPerSegment = qPathSamplingMode === "total"
    ? clampInt(
      totalQPointsTarget / minimumTotalQPoints,
      1,
      MAX_POINTS_PER_SEGMENT,
    )
    : pointsPerSegment;

  useEffect(() => {
    setTotalQPointsInput(String(totalQPointsTarget));
  }, [totalQPointsTarget]);

  const commitTotalQPointsInput = useCallback(() => {
    const parsed = Number.parseInt(totalQPointsInput.trim(), 10);
    const fallback = Number.isFinite(totalQPointsTarget)
      ? totalQPointsTarget
      : Math.max(120, minimumTotalQPoints);
    const committed = Number.isFinite(parsed)
      ? clampInt(parsed, minimumTotalQPoints, MAX_TOTAL_Q_POINTS)
      : clampInt(fallback, minimumTotalQPoints, MAX_TOTAL_Q_POINTS);
    setTotalQPointsTarget(committed);
    setTotalQPointsInput(String(committed));
  }, [minimumTotalQPoints, totalQPointsInput, totalQPointsTarget]);

  async function buildPhononTaskPlan(): Promise<PhononTaskPlan> {
    if (!selectedScf?.result) {
      throw new Error("No source SCF calculation selected");
    }

    const scfParams = selectedScf.parameters || {};
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

    const shouldCalculateDispersion = calculateDispersion && qPath.length >= 2;
    if (calculateDispersion && !shouldCalculateDispersion) {
      throw new Error("Select at least 2 Q-path points for dispersion, or disable dispersion.");
    }

    const sourceUsesPrimitive = sourceScfUsesPrimitiveCell(scfParams);
    const canUseSymmetryPrimitive =
      sourceUsesPrimitive &&
      resolvedSymmetry !== null &&
      resolvedSymmetry.standardizedPrimitiveAtoms.length > 0;
    if (shouldCalculateDispersion && sourceUsesPrimitive && !canUseSymmetryPrimitive) {
      throw new Error(
        "Selected SCF was run in a primitive cell, but symmetry conversion data is unavailable. Re-run the SCF or refresh the structure metadata.",
      );
    }

    const context = resolvePathTransformContext(crystalData, resolvedSymmetry);
    const effectiveRhombohedralConvention = qPathRhombohedralConvention ??
      defaultRhombohedralConventionForSetting(context.rhombohedralSetting ?? null);
    const converters = createPathCoordinateConverters(context, resolvedSymmetry);

    const parsedDosDeltaE = parseOptionalPositiveNumber(dosDeltaEInput);
    if (calculateDos && parsedDosDeltaE === null) {
      throw new Error("Invalid DOS frequency step. Use a positive number (e.g. 1.0).");
    }

    const nmixPh = parseOptionalPositiveInt(nmixPhInput);
    if (nmixPhInput.trim() && nmixPh === null) {
      throw new Error("Response mixing history must be a positive integer (`nmix_ph`).");
    }
    const niterPh = parseOptionalPositiveInt(niterPhInput);
    if (niterPhInput.trim() && niterPh === null) {
      throw new Error("Max DFPT iterations must be a positive integer (`niter_ph`).");
    }
    const alphaMix = parseOptionalPositiveNumber(alphaMixInput);
    if (alphaMixInput.trim() && alphaMix === null) {
      throw new Error("Initial response mixing factor must be a positive number (`alpha_mix(1)`).");
    }

    const startQ = parseOptionalPositiveInt(startQInput);
    if (startQInput.trim() && startQ === null) {
      throw new Error("Q-range start index must be a positive integer (`start_q`).");
    }
    const lastQ = parseOptionalPositiveInt(lastQInput);
    if (lastQInput.trim() && lastQ === null) {
      throw new Error("Q-range end index must be a positive integer (`last_q`).");
    }
    if (startQ !== null && lastQ !== null && startQ > lastQ) {
      throw new Error("Q-range start index (`start_q`) cannot be greater than end index (`last_q`).");
    }

    const startIrr = parseOptionalPositiveInt(startIrrInput);
    if (startIrrInput.trim() && startIrr === null) {
      throw new Error("Irreducible-representation start index must be a positive integer (`start_irr`).");
    }
    const lastIrr = parseOptionalPositiveInt(lastIrrInput);
    if (lastIrrInput.trim() && lastIrr === null) {
      throw new Error("Irreducible-representation end index must be a positive integer (`last_irr`).");
    }
    if (startIrr !== null && lastIrr !== null && startIrr > lastIrr) {
      throw new Error("Irreducible-representation start index (`start_irr`) cannot be greater than end index (`last_irr`).");
    }

    const electronPhonon = electronPhononMode === "none" ? null : electronPhononMode;
    const elPhSigma = electronPhonon ? parseOptionalPositiveNumber(elPhSigmaInput) : null;
    const elPhNsigma = electronPhonon ? parseOptionalPositiveInt(elPhNsigmaInput) : null;
    if (electronPhonon && elPhSigma === null) {
      throw new Error("E-ph smearing width must be a positive number when electron-phonon mode is enabled (`el_ph_sigma`).");
    }
    if (electronPhonon && elPhNsigma === null) {
      throw new Error("Number of smearing samples must be a positive integer when electron-phonon mode is enabled (`el_ph_nsigma`).");
    }

    const fildvscfPath = fildvscfEnabled ? (fildvscf.trim() || "dvscf") : null;
    const verbosityValue = verbosity === "default" ? null : verbosity;
    const shouldPreserveFullArtifacts = epwPreparationEnabled && preserveFullArtifacts;
    const transformedQPath = shouldCalculateDispersion
      ? mapPathCoordinates(
        qPath,
        canUseSymmetryPrimitive
          ? converters.toSymmetryPrimitiveCoords
          : converters.toInputConventionalCoords,
      ).map((point) => ({
        label: point.label,
        coords: point.coords,
        npoints: point.npoints,
      }))
      : null;
    const phononCellRepresentation = canUseSymmetryPrimitive
      ? "primitive_spglib"
      : "conventional_input";
    const phononConfig = {
      phonon: {
        prefix: "qcortado_scf",
        outdir: "./tmp",
        fildyn: "dynmat",
        nq: qGrid,
        tr2_ph: tr2Ph,
        ldisp: true,
        recover,
        trans,
        epsil,
        fpol,
        lraman,
        fildvscf: fildvscfPath,
        electron_phonon: electronPhonon,
        el_ph_sigma: elPhSigma,
        el_ph_nsigma: elPhNsigma,
        nmix_ph: nmixPh,
        niter_ph: niterPh,
        alpha_mix: alphaMix,
        start_q: startQ,
        last_q: lastQ,
        start_irr: startIrr,
        last_irr: lastIrr,
        verbosity: verbosityValue,
        asr,
      },
      calculate_dos: calculateDos,
      dos_grid: calculateDos ? dosGrid : null,
      dos_delta_e: calculateDos ? parsedDosDeltaE : null,
      calculate_dispersion: shouldCalculateDispersion,
      q_path: transformedQPath,
      points_per_segment: qPathSamplingMode === "segment" ? pointsPerSegment : viewerPointsPerSegment,
      project_id: projectId,
      scf_calc_id: selectedScf.id,
    };

    const taskLabel = `Phonon - ${crystalData?.formula_sum || ""}`;
    const taskParams = {
      config: phononConfig,
      workingDir: PHONON_WORK_DIR,
      mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
    };
    const pathString = shouldCalculateDispersion ? qPath.map((point) => point.label).join(" -> ") : "";
    const saveParameters = {
      source_scf_id: selectedScf.id,
      q_grid: qGrid,
      tr2_ph: tr2Ph,
      asr,
      recover,
      trans,
      epsil,
      fpol,
      lraman,
      nmix_ph: nmixPh,
      niter_ph: niterPh,
      alpha_mix: alphaMix,
      start_q: startQ,
      last_q: lastQ,
      start_irr: startIrr,
      last_irr: lastIrr,
      verbosity: verbosityValue,
      electron_phonon: electronPhonon,
      el_ph_sigma: elPhSigma,
      el_ph_nsigma: elPhNsigma,
      fildvscf: fildvscfPath,
      calculate_dos: calculateDos,
      dos_grid: calculateDos ? dosGrid : null,
      dos_delta_e: calculateDos ? parsedDosDeltaE : null,
      calculate_dispersion: shouldCalculateDispersion,
      q_path: pathString,
      q_path_sampling_mode: qPathSamplingMode,
      points_per_segment: qPathSamplingMode === "segment" ? pointsPerSegment : null,
      total_q_points_target: qPathSamplingMode === "total" ? totalDispersionQPoints : null,
      q_path_convention: context.centering === "R" ? effectiveRhombohedralConvention : null,
      q_path_rhombohedral_setting: context.centering === "R" ? (context.rhombohedralSetting ?? null) : null,
      cell_representation: phononCellRepresentation,
      symmetry_spacegroup: resolvedSymmetry?.spacegroupNumber ?? null,
      symmetry_hall_number: resolvedSymmetry?.hallNumber ?? null,
      primitive_to_input_reciprocal: resolvedSymmetry?.primitiveToInputReciprocal ?? null,
      symmetry_error: resolvedSymmetryError,
      epw_preparation: {
        enabled: epwPreparationEnabled,
        preserve_full_artifacts: shouldPreserveFullArtifacts,
        fildvscf: fildvscfPath,
        electron_phonon: electronPhonon,
      },
      n_qpoints: null,
      n_modes: null,
    };

    return {
      taskLabel,
      taskParams,
      saveParameters,
      progressMeta: {
        hasDos: calculateDos,
        hasDispersion: shouldCalculateDispersion,
      },
    };
  }

  // Run the phonon calculation
  const runCalculation = async () => {
    if (hasExternalRunningTask) {
      setError("Another task is currently running. Queue this task or wait for completion.");
      return;
    }
    let plan: PhononTaskPlan;
    try {
      plan = await buildPhononTaskPlan();
    } catch (e) {
      setError(String(e));
      return;
    }

    setIsRunning(true);
    followOutputRef.current = true;
    setOutput("");
    setError(null);
    setPhononResult(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Phonon calculation", {
      phonon: {
        hasDos: plan.progressMeta.hasDos,
        hasDispersion: plan.progressMeta.hasDispersion,
      },
    }));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

    try {
      const taskId = await taskContext.startTask("phonon", plan.taskParams, plan.taskLabel);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Calculation failed");
      }

      const result = finalTask.result as PhononResult;
      const outputContent = finalTask.output.join("\n");
      const endTime = new Date().toISOString();
      setPhononResult(result);
      setStep("results");
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      // Auto-save the phonon calculation to the project
      try {
        await invoke("save_calculation", {
          projectId,
          cifId: _cifId,
          calcData: {
            calc_type: "phonon",
            parameters: {
              ...plan.saveParameters,
              n_qpoints: result.n_qpoints,
              n_modes: result.n_modes,
            },
            result: {
              converged: result.converged,
              total_energy: null,
              fermi_energy: null,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: result.raw_output,
              phonon_data: {
                dos_data: result.dos_data,
                dispersion_data: result.dispersion_data,
              },
            },
            started_at: startTime,
            completed_at: endTime,
            input_content: "",
            output_content: outputContent,
          },
          workingDir: PHONON_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save phonon calculation:", saveError);
        setError(`Failed to auto-save phonon calculation: ${saveError}`);
      }
    } catch (e) {
      setError(String(e));
      setOutput((prev) => prev + `\nError: ${e}\n`);
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
    try {
      const plan = await buildPhononTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "phonon",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId: _cifId,
          workingDir: PHONON_WORK_DIR,
          calcType: "phonon",
          parameters: plan.saveParameters,
          inputContent: "",
        },
      );
    } catch (e) {
      setError(String(e));
    }
  };

  // Check if SCF is phonon-ready (optimized structure + conv_thr <= 1e-12).
  const isPhononReady = (calc: CalculationRun): boolean => {
    return isPhononReadyScf(calc.parameters, calc.tags);
  };

  // Render current step
  const renderStep = () => {
    switch (step) {
      case "source":
        return renderSourceStep();
      case "qgrid":
        return renderQGridStep();
      case "options":
        return renderOptionsStep();
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
    const phononReadyScfs = validScfs.filter(isPhononReady);

    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Phonon calculations require a converged SCF calculation.
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
            <label htmlFor="phonon-scf-sort">Sort SCFs</label>
            <select
              id="phonon-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose an SCF calculation to use as the starting point for phonons.
          For best results, use a phonon-ready SCF from an optimized structure
          with conv_thr &lt;= 1e-12.
        </p>

        {phononReadyScfs.length === 0 && (
          <div className="warning-banner">
            No phonon-ready SCF calculations found. Run an SCF on an optimized
            structure with the "Phonon-Ready" preset (conv_thr = 1e-12).
          </div>
        )}

        <div className="scf-list">
          {sortedScfs.map((scf) => {
            const ready = isPhononReady(scf);
            return (
              <div
                key={scf.id}
                className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""} ${ready ? "phonon-ready" : ""}`}
                onClick={() => setSelectedScf(scf)}
              >
                <div className="scf-option-header">
                  <input
                    type="radio"
                    checked={selectedScf?.id === scf.id}
                    onChange={() => setSelectedScf(scf)}
                  />
                  <span className="scf-date">
                    {new Date(scf.started_at).toLocaleDateString()}
                  </span>
                  {ready && <span className="phonon-ready-badge">Phonon-ready</span>}
                </div>
                <div className="scf-details">
                  <span>E = {scf.result?.total_energy?.toFixed(6)} Ry</span>
                  {scf.result?.fermi_energy && (
                    <span>E_F = {scf.result.fermi_energy.toFixed(3)} eV</span>
                  )}
                </div>
                <div className="calc-tags">
                  {getCalculationTags(scf).map((tag, i) => (
                    <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
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
            onClick={() => setStep("qgrid")}
          >
            Next: ph.x Setup
          </button>
        </div>
      </div>
    );
  };

  // Step 2: ph.x + EPW configuration
  const renderQGridStep = () => {
    const totalQPoints = qGrid[0] * qGrid[1] * qGrid[2];
    const parsedTr2Ph = Number(tr2PhInput);
    const isTr2PhValid = Number.isFinite(parsedTr2Ph) && parsedTr2Ph > 0;
    const isNmixValid = isOptionalPositiveIntValid(nmixPhInput);
    const isNiterValid = isOptionalPositiveIntValid(niterPhInput);
    const isAlphaMixValid = isOptionalPositiveNumberValid(alphaMixInput);
    const isStartQValid = isOptionalPositiveIntValid(startQInput);
    const isLastQValid = isOptionalPositiveIntValid(lastQInput);
    const isStartIrrValid = isOptionalPositiveIntValid(startIrrInput);
    const isLastIrrValid = isOptionalPositiveIntValid(lastIrrInput);
    const parsedStartQ = parseOptionalPositiveInt(startQInput);
    const parsedLastQ = parseOptionalPositiveInt(lastQInput);
    const parsedStartIrr = parseOptionalPositiveInt(startIrrInput);
    const parsedLastIrr = parseOptionalPositiveInt(lastIrrInput);
    const hasQRangeError =
      parsedStartQ !== null && parsedLastQ !== null && parsedStartQ > parsedLastQ;
    const hasIrrRangeError =
      parsedStartIrr !== null && parsedLastIrr !== null && parsedStartIrr > parsedLastIrr;
    const electronPhEnabled = electronPhononMode !== "none";
    const isElPhSigmaValid =
      !electronPhEnabled || parseOptionalPositiveNumber(elPhSigmaInput) !== null;
    const isElPhNsigmaValid =
      !electronPhEnabled || parseOptionalPositiveInt(elPhNsigmaInput) !== null;
    const isFildvscfValid = !fildvscfEnabled || fildvscf.trim().length > 0;
    const canContinue =
      isTr2PhValid &&
      isNmixValid &&
      isNiterValid &&
      isAlphaMixValid &&
      isStartQValid &&
      isLastQValid &&
      isStartIrrValid &&
      isLastIrrValid &&
      !hasQRangeError &&
      !hasIrrRangeError &&
      isElPhSigmaValid &&
      isElPhNsigmaValid &&
      isFildvscfValid;

    return (
      <div className="wizard-step qgrid-step">
        <h3>ph.x + EPW Setup</h3>
        <p className="step-description">
          Configure `ph.x` input options and optional EPW-ready outputs. Denser grids and tighter thresholds improve accuracy but increase compute time.
        </p>

        <div className="option-section">
          <h4>Core ph.x Controls</h4>
          <div className="phonon-grid">
            <div className="phonon-field">
              <label>
                Q-point grid
                <Tooltip text="Phonon dynamical matrices are computed on this coarse reciprocal-space mesh (`nq1`, `nq2`, `nq3`). Larger meshes improve interpolation and thermodynamics, but cost scales quickly." />
              </label>
              <div className="qgrid-inputs">
                <input
                  type="number"
                  min={1}
                  max={24}
                  value={qGrid[0]}
                  onChange={(e) => setQGrid([parseInt(e.target.value) || 1, qGrid[1], qGrid[2]])}
                />
                <span>x</span>
                <input
                  type="number"
                  min={1}
                  max={24}
                  value={qGrid[1]}
                  onChange={(e) => setQGrid([qGrid[0], parseInt(e.target.value) || 1, qGrid[2]])}
                />
                <span>x</span>
                <input
                  type="number"
                  min={1}
                  max={24}
                  value={qGrid[2]}
                  onChange={(e) => setQGrid([qGrid[0], qGrid[1], parseInt(e.target.value) || 1])}
                />
              </div>
              <span className="param-hint">{totalQPoints} total q-points</span>
            </div>

            <div className="phonon-field">
              <label>
                DFPT convergence threshold
                <Tooltip text="QE variable: `tr2_ph`. Linear-response convergence threshold. Tighter values reduce force-constant noise but increase SCF iterations per perturbation." />
              </label>
              <input
                type="text"
                value={tr2PhInput}
                onChange={(e) => {
                  const value = e.target.value.trim();
                  setTr2PhInput(value);
                  const parsed = Number(value);
                  if (Number.isFinite(parsed) && parsed > 0) {
                    setTr2Ph(parsed);
                  }
                }}
                placeholder="e.g. 1e-12"
                spellCheck={false}
              />
              {!isTr2PhValid && (
                <span className="param-hint input-error">Enter a positive number (e.g. 1e-12).</span>
              )}
            </div>

            <div className="phonon-field">
              <label>
                Acoustic sum rule
                <Tooltip text="QE variable: `asr` (used by `q2r.x`/`matdyn.x`). Enforces translational invariance of acoustic modes. `crystal` is typically best for periodic bulk systems." />
              </label>
              <select
                value={asr}
                onChange={(e) => setAsr(e.target.value as "none" | "simple" | "crystal" | "one-dim" | "zero-dim")}
              >
                <option value="crystal">crystal</option>
                <option value="simple">simple</option>
                <option value="none">none</option>
                <option value="one-dim">one-dim</option>
                <option value="zero-dim">zero-dim</option>
              </select>
            </div>
          </div>
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("solver")} className="section-header">
            <span className={`collapse-icon ${expandedSections.solver ? "expanded" : ""}`}>▶</span>
            Convergence & Solver Settings
            <Tooltip text="Advanced `ph.x` solver controls (`nmix_ph`, `niter_ph`, `alpha_mix(1)`, `verbosity`)." />
          </h4>
          {expandedSections.solver && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    Response mixing history
                    <Tooltip text="QE variable: `nmix_ph`. Number of mixed iterations retained in linear-response SCF. Larger values can stabilize difficult cases but cost more memory/time." />
                  </label>
                  <input
                    type="text"
                    value={nmixPhInput}
                    onChange={(e) => setNmixPhInput(e.target.value)}
                    placeholder="default"
                    spellCheck={false}
                  />
                  {!isNmixValid && <span className="param-hint input-error">Use a positive integer.</span>}
                </div>

                <div className="phonon-field">
                  <label>
                    Max DFPT iterations
                    <Tooltip text="QE variable: `niter_ph`. Maximum DFPT self-consistency iterations per perturbation. Raise for hard-to-converge metallic or magnetic systems." />
                  </label>
                  <input
                    type="text"
                    value={niterPhInput}
                    onChange={(e) => setNiterPhInput(e.target.value)}
                    placeholder="default"
                    spellCheck={false}
                  />
                  {!isNiterValid && <span className="param-hint input-error">Use a positive integer.</span>}
                </div>

                <div className="phonon-field">
                  <label>
                    Initial response mixing factor
                    <Tooltip text="QE variable: `alpha_mix(1)`. Lower values can help convergence in difficult cases, but usually require more iterations." />
                  </label>
                  <input
                    type="text"
                    value={alphaMixInput}
                    onChange={(e) => setAlphaMixInput(e.target.value)}
                    placeholder="default"
                    spellCheck={false}
                  />
                  {!isAlphaMixValid && <span className="param-hint input-error">Use a positive number.</span>}
                </div>

                <div className="phonon-field">
                  <label>
                    Output detail level
                    <Tooltip text="QE variable: `verbosity`. Higher verbosity helps debugging but increases output size." />
                  </label>
                  <select
                    value={verbosity}
                    onChange={(e) => setVerbosity(e.target.value as "default" | "high" | "debug")}
                  >
                    <option value="default">default</option>
                    <option value="high">high</option>
                    <option value="debug">debug</option>
                  </select>
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("longRange")} className="section-header">
            <span className={`collapse-icon ${expandedSections.longRange ? "expanded" : ""}`}>▶</span>
            Long-Range & Spectroscopy
            <Tooltip text="Optional long-range electrostatics and spectroscopy outputs (`epsil`, `fpol`, `lraman`)." />
          </h4>
          {expandedSections.longRange && (
            <div className="option-params">
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={epsil}
                  onChange={(e) => setEpsil(e.target.checked)}
                />
                <span>
                  Compute dielectric tensor and Born charges
                  <Tooltip text="QE variable: `epsil`. Computes dielectric tensor and Born effective charges for non-metallic systems; required for LO-TO splitting analysis at Γ." />
                </span>
              </label>
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={fpol}
                  onChange={(e) => setFpol(e.target.checked)}
                />
                <span>
                  Compute polarizability derivatives
                  <Tooltip text="QE variable: `fpol`. Enables additional response quantities used by optical analyses; adds DFPT workload." />
                </span>
              </label>
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={lraman}
                  onChange={(e) => setLraman(e.target.checked)}
                />
                <span>
                  Compute Raman tensors
                  <Tooltip text="QE variable: `lraman`. Computes Raman tensors from DFPT response. Useful for spectroscopy, but significantly increases runtime." />
                </span>
              </label>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("electronPhonon")} className="section-header">
            <span className={`collapse-icon ${expandedSections.electronPhonon ? "expanded" : ""}`}>▶</span>
            Electron-Phonon & Partial Run
            <Tooltip text="Electron-phonon mode and optional partial-range execution (`electron_phonon`, `el_ph_sigma`, `el_ph_nsigma`, `start_q`, `last_q`, `start_irr`, `last_irr`)." />
          </h4>
          {expandedSections.electronPhonon && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    Electron-phonon mode
                    <Tooltip text="QE variable: `electron_phonon`. Keep `none` for pure vibrational runs; enable for coupling workflows." />
                  </label>
                  <select
                    value={electronPhononMode}
                    onChange={(e) => setElectronPhononMode(e.target.value as "none" | "interpolated" | "lambda_tetra" | "gamma_tetra")}
                  >
                    <option value="none">none</option>
                    <option value="interpolated">interpolated</option>
                    <option value="lambda_tetra">lambda_tetra</option>
                    <option value="gamma_tetra">gamma_tetra</option>
                  </select>
                </div>

                <div className="phonon-field">
                  <label>
                    E-ph smearing width
                    <Tooltip text="QE variable: `el_ph_sigma`. Smearing width used in electron-phonon integrations. Larger broadening smooths noisy integrals but can blur fine features." />
                  </label>
                  <input
                    type="text"
                    value={elPhSigmaInput}
                    onChange={(e) => setElPhSigmaInput(e.target.value)}
                    placeholder="e.g. 0.02"
                    spellCheck={false}
                  />
                  {!isElPhSigmaValid && (
                    <span className="param-hint input-error">Use a positive number when electron-phonon mode is enabled.</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Number of smearing samples
                    <Tooltip text="QE variable: `el_ph_nsigma`. Number of sigma values sampled in electron-phonon integrations. More values improve convergence checks but increase cost." />
                  </label>
                  <input
                    type="text"
                    value={elPhNsigmaInput}
                    onChange={(e) => setElPhNsigmaInput(e.target.value)}
                    placeholder="e.g. 10"
                    spellCheck={false}
                  />
                  {!isElPhNsigmaValid && (
                    <span className="param-hint input-error">Use a positive integer when electron-phonon mode is enabled.</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Q-point index range
                    <Tooltip text="QE variables: `start_q`, `last_q`. Run only a subset of irreducible q-points (1-based indices). Useful for restart/debug or distributed job splitting." />
                  </label>
                  <div className="dual-input">
                    <input
                      type="text"
                      value={startQInput}
                      onChange={(e) => setStartQInput(e.target.value)}
                      placeholder="start_q"
                      spellCheck={false}
                    />
                    <input
                      type="text"
                      value={lastQInput}
                      onChange={(e) => setLastQInput(e.target.value)}
                      placeholder="last_q"
                      spellCheck={false}
                    />
                  </div>
                  {(!isStartQValid || !isLastQValid) && (
                    <span className="param-hint input-error">Use positive integers.</span>
                  )}
                  {hasQRangeError && (
                    <span className="param-hint input-error">Start index must be less than or equal to end index.</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Irreducible-representation range
                    <Tooltip text="QE variables: `start_irr`, `last_irr`. Restrict irreducible representation range for each q-point. Primarily for advanced restart and debugging workflows." />
                  </label>
                  <div className="dual-input">
                    <input
                      type="text"
                      value={startIrrInput}
                      onChange={(e) => setStartIrrInput(e.target.value)}
                      placeholder="start_irr"
                      spellCheck={false}
                    />
                    <input
                      type="text"
                      value={lastIrrInput}
                      onChange={(e) => setLastIrrInput(e.target.value)}
                      placeholder="last_irr"
                      spellCheck={false}
                    />
                  </div>
                  {(!isStartIrrValid || !isLastIrrValid) && (
                    <span className="param-hint input-error">Use positive integers.</span>
                  )}
                  {hasIrrRangeError && (
                    <span className="param-hint input-error">Start index must be less than or equal to end index.</span>
                  )}
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="option-section epw-prep-section config-section collapsible">
          <h4 onClick={() => toggleSection("epw")} className="section-header">
            <span className={`collapse-icon ${expandedSections.epw ? "expanded" : ""}`}>▶</span>
            EPW Preparation
            <Tooltip text="Prepare additional artifacts for later EPW workflows (`fildvscf`, artifact retention)." />
          </h4>
          {expandedSections.epw && (
            <div className="option-params">
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={epwPreparationEnabled}
                  onChange={(e) => setEpwPreparationEnabled(e.target.checked)}
                />
                <span>
                  Enable EPW-ready artifact preparation
                  <Tooltip text="Turns on options commonly needed before EPW (dvscf outputs and preservation controls). This does not run `epw.x`; it prepares prerequisite phonon data." />
                </span>
              </label>
              <div className="option-params">
                <label className="option-checkbox">
                  <input
                    type="checkbox"
                    checked={fildvscfEnabled}
                    onChange={(e) => setFildvscfEnabled(e.target.checked)}
                  />
                  <span>
                    Write perturbation potentials
                    <Tooltip text="QE variable: `fildvscf`. Stores dvscf files used by electron-phonon workflows. Increases on-disk data substantially." />
                  </span>
                </label>
                {fildvscfEnabled && (
                  <div className="phonon-field inline-field">
                    <label>
                      Perturbation potential filename root
                      <Tooltip text="QE variable: `fildvscf`. Base filename used for perturbation-potential files." />
                    </label>
                    <input
                      type="text"
                      value={fildvscf}
                      onChange={(e) => setFildvscf(e.target.value)}
                      placeholder="dvscf"
                      spellCheck={false}
                    />
                    {!isFildvscfValid && (
                      <span className="param-hint input-error">Filename root cannot be empty when enabled.</span>
                    )}
                  </div>
                )}
                <label className="option-checkbox">
                  <input
                    type="checkbox"
                    checked={preserveFullArtifacts}
                    onChange={(e) => setPreserveFullArtifacts(e.target.checked)}
                  />
                  <span>
                    Preserve full phonon scratch artifacts on save
                    <Tooltip text="Keeps complete phonon temporary data in the project instead of compact copy. Needed for maximum EPW flexibility, but can consume many GB." />
                  </span>
                </label>
                {epwPreparationEnabled && !preserveFullArtifacts && (
                  <div className="epw-warning">
                    EPW prep is enabled, but compact-save mode may omit intermediate files needed by some EPW post-processing paths.
                  </div>
                )}
              </div>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("execution")} className="section-header">
            <span className={`collapse-icon ${expandedSections.execution ? "expanded" : ""}`}>▶</span>
            Execution Behavior
            <Tooltip text="Runtime behavior flags (`recover`, `trans`)." />
          </h4>
          {expandedSections.execution && (
            <div className="option-params">
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={recover}
                  onChange={(e) => setRecover(e.target.checked)}
                />
                <span>
                  Resume from checkpoints
                  <Tooltip text="QE variable: `recover`. Attempts to resume from existing phonon checkpoints in the working directory. Helpful after interrupted long runs." />
                </span>
              </label>
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={trans}
                  onChange={(e) => setTrans(e.target.checked)}
                />
                <span>
                  Use Cartesian-transform output
                  <Tooltip text="QE variable: `trans`. Controls representation transform in `ph.x` output. Leave enabled unless you need specific low-level output conventions." />
                </span>
              </label>
            </div>
          )}
        </div>

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button
            className="primary-button"
            disabled={!canContinue}
            onClick={() => {
              setTr2Ph(parsedTr2Ph);
              setStep("options");
            }}
          >
            Next: Post-Processing
          </button>
        </div>
      </div>
    );
  };

  // Step 3: post-processing and run
  const renderOptionsStep = () => {
    const parsedDosDeltaE = parseOptionalPositiveNumber(dosDeltaEInput);
    const isDosDeltaEValid = !calculateDos || parsedDosDeltaE !== null;
    const dispersionReady = !calculateDispersion || qPath.length >= 2;
    const canRun = (calculateDos || calculateDispersion) && isDosDeltaEValid && dispersionReady;

    return (
      <div className="wizard-step options-step">
        <h3>Post-Processing + Run</h3>
        <p className="step-description">
          Configure `matdyn.x` post-processing (DOS and/or dispersion) and parallel execution.
        </p>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("dos")} className="section-header">
            <span className={`collapse-icon ${expandedSections.dos ? "expanded" : ""}`}>▶</span>
            Phonon Density of States (matdyn.x)
            <Tooltip text="Configure DOS post-processing with `matdyn.x` (`dos=.true.`, `nk1/nk2/nk3`, `deltaE`)." />
          </h4>
          {expandedSections.dos && (
            <div className="option-params">
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={calculateDos}
                  onChange={(e) => setCalculateDos(e.target.checked)}
                />
                <span>
                  Compute phonon density of states
                  <Tooltip text="Computes vibrational density of states from force constants. Useful for thermodynamics and spectral comparisons." />
                </span>
              </label>

              {calculateDos && (
                <div className="option-params">
                  <div className="phonon-grid">
                    <div className="phonon-field">
                      <label>
                        DOS sampling grid
                        <Tooltip text="`matdyn.x` integration grid (`nk1`, `nk2`, `nk3`). Dense q sampling smooths DOS and improves integrated thermodynamic quantities, at additional cost." />
                      </label>
                      <div className="qgrid-inputs">
                        <input
                          type="number"
                          min={5}
                          max={60}
                          value={dosGrid[0]}
                          onChange={(e) => setDosGrid([parseInt(e.target.value) || 10, dosGrid[1], dosGrid[2]])}
                        />
                        <span>x</span>
                        <input
                          type="number"
                          min={5}
                          max={60}
                          value={dosGrid[1]}
                          onChange={(e) => setDosGrid([dosGrid[0], parseInt(e.target.value) || 10, dosGrid[2]])}
                        />
                        <span>x</span>
                        <input
                          type="number"
                          min={5}
                          max={60}
                          value={dosGrid[2]}
                          onChange={(e) => setDosGrid([dosGrid[0], dosGrid[1], parseInt(e.target.value) || 10])}
                        />
                      </div>
                    </div>
                    <div className="phonon-field">
                      <label>
                        DOS frequency bin width (cm^-1)
                        <Tooltip text="`matdyn.x` variable: `deltaE`. Smaller values give finer spectra but larger files and potentially noisier curves." />
                      </label>
                      <input
                        type="text"
                        value={dosDeltaEInput}
                        onChange={(e) => setDosDeltaEInput(e.target.value)}
                        placeholder="1.0"
                        spellCheck={false}
                      />
                      {!isDosDeltaEValid && (
                        <span className="param-hint input-error">Use a positive number.</span>
                      )}
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("dispersion")} className="section-header">
            <span className={`collapse-icon ${expandedSections.dispersion ? "expanded" : ""}`}>▶</span>
            Phonon Dispersion (matdyn.x)
            <Tooltip text="Configure band-structure-style phonon interpolation along a high-symmetry path." />
          </h4>
          {expandedSections.dispersion && (
            <div className="option-params">
              <label className="option-checkbox">
                <input
                  type="checkbox"
                  checked={calculateDispersion}
                  onChange={(e) => setCalculateDispersion(e.target.checked)}
                />
                <span>
                  Compute phonon dispersion
                  <Tooltip text="Computes phonon branches along a high-symmetry path. Essential for diagnosing dynamical stability and mode character." />
                </span>
              </label>

              {calculateDispersion && (
                <div className="option-params">
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
                    onPathChange={handleQPathChange}
                    initialPath={qPath}
                    pointsPerSegment={viewerPointsPerSegment}
                    symmetryTransform={symmetryTransform}
                    rhombohedralConvention={qPathRhombohedralConvention}
                    onRhombohedralConventionChange={setQPathRhombohedralConvention}
                  />

                  <div className="kpath-sampling-panel">
                    <div className="kpath-sampling-header">
                      <div>
                        <h4>Q-Path Sampling</h4>
                        <p>Set interpolation density along the selected high-symmetry path.</p>
                      </div>
                      <span className="kpath-sampling-summary">
                        {totalDispersionQPoints} total q-points
                      </span>
                    </div>

                    <div className="phonon-unit-toggle kpath-sampling-toggle" role="group" aria-label="Q-path sampling mode">
                      <button
                        type="button"
                        className={`phonon-unit-btn ${qPathSamplingMode === "segment" ? "active" : ""}`}
                        onClick={() => setQPathSamplingMode("segment")}
                      >
                        Points / segment
                      </button>
                      <button
                        type="button"
                        className={`phonon-unit-btn ${qPathSamplingMode === "total" ? "active" : ""}`}
                        onClick={() => {
                          if (totalDispersionQPoints > 0) {
                            setTotalQPointsTarget(totalDispersionQPoints);
                            setTotalQPointsInput(String(totalDispersionQPoints));
                          }
                          setQPathSamplingMode("total");
                        }}
                      >
                        Total q-points
                      </button>
                    </div>

                    {qPathSamplingMode === "segment" ? (
                      <label className="kpath-sampling-input">
                        <span>
                          Points per path segment
                          <Tooltip text="Number of interpolated points between adjacent high-symmetry nodes in `matdyn.x` path generation. Higher values produce smoother curves with longer runs." />
                        </span>
                        <input
                          type="number"
                          min={MIN_POINTS_PER_SEGMENT}
                          max={MAX_POINTS_PER_SEGMENT}
                          value={pointsPerSegment}
                          onChange={(e) => {
                            const parsed = Number.parseInt(e.target.value, 10);
                            setPointsPerSegment(
                              Number.isFinite(parsed)
                                ? clampInt(parsed, MIN_POINTS_PER_SEGMENT, MAX_POINTS_PER_SEGMENT)
                                : 20,
                            );
                          }}
                        />
                      </label>
                    ) : (
                      <label className="kpath-sampling-input">
                        <span>Total q-points</span>
                        <input
                          type="number"
                          min={minimumTotalQPoints}
                          max={MAX_TOTAL_Q_POINTS}
                          value={totalQPointsInput}
                          onChange={(e) => {
                            setTotalQPointsInput(e.target.value);
                          }}
                          onBlur={commitTotalQPointsInput}
                          onKeyDown={(e) => {
                            if (e.key === "Enter") {
                              e.preventDefault();
                              commitTotalQPointsInput();
                            }
                          }}
                        />
                      </label>
                    )}

                    {qPathSamplingMode === "total" && (
                      <p className="kpath-sampling-note">
                        {qPathSegmentCount > 0
                          ? `Evenly distributed by segment length across ${qPathSegmentCount} path segment${qPathSegmentCount === 1 ? "" : "s"}.`
                          : "Add at least 2 points to distribute q-points along the path."}
                      </p>
                    )}

                    {!dispersionReady && (
                      <span className="param-hint input-error">Select at least 2 Q-path points.</span>
                    )}
                  </div>
                </div>
              )}
            </div>
          )}
        </div>

        <div className="option-section mpi-section config-section collapsible">
          <h4 onClick={() => toggleSection("mpi")} className="section-header">
            <span className={`collapse-icon ${expandedSections.mpi ? "expanded" : ""}`}>▶</span>
            Parallelization
            <Tooltip text="Controls process-level parallel execution. More MPI processes often reduce wall time, but scaling depends on system size and workload balance." />
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
                          onChange={(e) => setMpiProcs(parseInt(e.target.value) || 1)}
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

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("qgrid")}>
            Back
          </button>
          <button
            className="secondary-button"
            disabled={!canRun}
            onClick={queueCalculation}
          >
            Queue Task
          </button>
          <button
            className="primary-button"
            disabled={!canRun || hasExternalRunningTask}
            onClick={runCalculation}
          >
            Run Calculation
          </button>
        </div>
      </div>
    );
  };

  // Step 4: Running
  const renderRunStep = () => {
    const phononMeta = progress.meta?.phonon;
    const completedQPoints = phononMeta?.completedQPoints
      ? new Set(phononMeta.completedQPoints).size
      : 0;
    const totalQPoints = phononMeta?.totalQPoints;

    return (
      <div className="wizard-step run-step">
        <h3>{isRunning ? "Running Phonon Calculation" : "Phonon Output"}</h3>
        <p className="step-description">
          This may take a while depending on system size and q-grid density.
        </p>

        <ProgressBar
          status={progress.status}
          percent={progress.percent}
          phase={progress.phase}
          detail={progress.detail}
        />
        <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />
        <EstimatedRemainingTime
          startedAt={calcStartTime}
          isRunning={isRunning}
          completedUnits={completedQPoints}
          totalUnits={totalQPoints}
          label="ETA (q-points)"
        />

        <pre ref={outputRef} className="calculation-output" onScroll={handleOutputScroll}>
          {output || "Starting calculation..."}
        </pre>

        {error && (
          <div className="error-actions">
            <button className="secondary-button" onClick={() => setStep("options")}>
              Back to Post-Processing
            </button>
          </div>
        )}
      </div>
    );
  };

  // Step 5: Results
  const renderResultsStep = () => {
    if (!phononResult) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("options")}>
            Back to Post-Processing
          </button>
        </div>
      );
    }

    const hasDos = phononResult.dos_data !== null;
    const hasDispersion = phononResult.dispersion_data !== null;

    return (
      <div className="wizard-step results-step">
        <h3>Phonon Calculation Results</h3>
        <p className="step-description">
          Calculation complete. Use the main Phonon viewer for full plotting controls.
        </p>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Q-points:</span>
              <span className="value">{phononResult.n_qpoints}</span>
            </div>
            <div className="summary-item">
              <span className="label">Modes:</span>
              <span className="value">{phononResult.n_modes}</span>
            </div>
            {hasDispersion && phononResult.dispersion_data && (
              <div className="summary-item">
                <span className="label">Frequency Range:</span>
                <span className="value">
                  {phononResult.dispersion_data.frequency_range[0].toFixed(1)} to{" "}
                  {phononResult.dispersion_data.frequency_range[1].toFixed(1)} cm^-1
                </span>
              </div>
            )}
          </div>
        </div>

        {/* Collapsible calculation output */}
        <div className="output-section">
          <button
            className="output-toggle"
            onClick={() => setShowOutput(!showOutput)}
          >
            {showOutput ? "v" : ">"} Calculation Output
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
          {hasDispersion && (
            <button
              className="primary-button"
              onClick={() =>
                onViewPhonons(
                  {
                    dos_data: phononResult.dos_data,
                    dispersion_data: phononResult.dispersion_data,
                  },
                  "bands",
                )
              }
            >
              View Phonon Bands
            </button>
          )}
          {hasDos && (
            <button
              className="primary-button"
              onClick={() =>
                onViewPhonons(
                  {
                    dos_data: phononResult.dos_data,
                    dispersion_data: phononResult.dispersion_data,
                  },
                  "dos",
                )
              }
            >
              View Phonon DOS
            </button>
          )}
          <button className="primary-button" onClick={() => setStep("options")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className="phonon-wizard">
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Phonon Calculation Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "qgrid" ? "active" : ["options", "run", "results"].includes(step) ? "completed" : ""}>
            2. ph.x + EPW
          </span>
          <span className={step === "options" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            3. Post-Processing
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
