// Band Structure Wizard - Calculate and visualize electronic band structures

import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData, ELEMENT_MASSES } from "../lib/types";
import { BandData } from "./BandPlot";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";
import { Vec3 } from "../lib/reciprocalLattice";
import {
  analyzeCrystalSymmetry,
  buildConventionalLatticeFromCrystalData,
  SymmetryTransformResult,
} from "../lib/symmetryTransform";
import {
  createPathCoordinateConverters,
  mapPathCoordinates,
  roundVec3,
  resolvePathTransformContext,
  sourceScfUsesPrimitiveCell,
} from "../lib/kPathTransforms";
import {
  RhombohedralConvention,
  defaultRhombohedralConventionForSetting,
} from "../lib/brillouinZoneData";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { isPhononReadyScf } from "../lib/phononReady";

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

function normalizeSmearing(raw: unknown): "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" {
  const lowered = String(raw || "gaussian").toLowerCase();
  if (lowered === "methfessel-paxton") return "methfessel-paxton";
  if (lowered === "marzari-vanderbilt") return "marzari-vanderbilt";
  if (lowered === "fermi-dirac") return "fermi-dirac";
  return "gaussian";
}

type KPathSamplingMode = "segment" | "total";

const MIN_POINTS_PER_SEGMENT = 5;
const MAX_POINTS_PER_SEGMENT = 400;
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

function normalizeKPathSampling(
  path: KPathPoint[],
  samplingMode: KPathSamplingMode,
  pointsPerSegment: number,
  totalKPoints: number,
): KPathPoint[] {
  if (samplingMode === "total") {
    return applyTotalKPoints(path, totalKPoints);
  }
  return applyPointsPerSegment(path, pointsPerSegment);
}

interface BandStructureWizardProps {
  onBack: () => void;
  onViewBands: (bandData: BandData, scfFermiEnergy: number | null) => void;
  qePath: string;
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
  onViewBands,
  qePath,
  projectId,
  cifId: _cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: BandStructureWizardProps) {
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

  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);

  // Step 2: K-Path (using BrillouinZoneViewer)
  const [kPath, setKPath] = useState<KPathPoint[]>([]);
  const [kPathRhombohedralConvention, setKPathRhombohedralConvention] = useState<RhombohedralConvention | undefined>(undefined);
  const [kPathSamplingMode, setKPathSamplingMode] = useState<KPathSamplingMode>("segment");
  const [pointsPerSegment, setPointsPerSegment] = useState(20);
  const [totalKPointsTarget, setTotalKPointsTarget] = useState(120);
  const [totalKPointsInput, setTotalKPointsInput] = useState("120");

  // Step 3: Parameters
  const [nbnd, setNbnd] = useState<number | "auto">("auto");
  const [nscfConvThrInput, setNscfConvThrInput] = useState("1e-8");
  const [nscfMixingBetaInput, setNscfMixingBetaInput] = useState("0.7");
  const [nscfOccupations, setNscfOccupations] = useState<"fixed" | "smearing" | "from_input" | "tetrahedra">("smearing");
  const [nscfSmearing, setNscfSmearing] = useState<"gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac">("gaussian");
  const [nscfDegaussInput, setNscfDegaussInput] = useState("0.02");
  const [nscfVerbosity, setNscfVerbosity] = useState<"low" | "high" | "debug">("high");
  const [bandsFilbandInput, setBandsFilbandInput] = useState("bands.dat");
  const [bandsLsym, setBandsLsym] = useState(true);
  const [bandsNoOverlap, setBandsNoOverlap] = useState(true);
  const [enableProjections, setEnableProjections] = useState(true);
  const [projectionLsym, setProjectionLsym] = useState(false);
  const [projectionDiagBasis, setProjectionDiagBasis] = useState(false);
  const [projectionPawproj, setProjectionPawproj] = useState(false);
  const [projectionFilprojInput, setProjectionFilprojInput] = useState("bands.projwfc.dat");
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
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);

  const handleOutputScroll = () => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followOutputRef.current = distanceToBottom <= 16;
  };
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
  const [showOutput, setShowOutput] = useState(false);
  // Track if calculation was saved
  const [isSaved, setIsSaved] = useState(false);
  // Track calculation timing
  const [calcStartTime, setCalcStartTime] = useState<string>("");

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);

  // Pseudopotentials (pseudopotentials list used internally for auto-selection)
  const [, setPseudopotentials] = useState<string[]>([]);
  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});
  const [symmetryTransform, setSymmetryTransform] = useState<SymmetryTransformResult | null>(null);
  const [symmetryError, setSymmetryError] = useState<string | null>(null);

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
    setNscfSmearing(normalizeSmearing(params.smearing));

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
  }, []);

  const selectSourceScf = useCallback((scf: CalculationRun) => {
    setSelectedScf(scf);
    setScfFermiEnergy(scf.result?.fermi_energy ?? null);
    applyScfDefaults(scf);
  }, [applyScfDefaults]);

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
        const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
        const pseudos = await invoke<string[]>("list_pseudopotentials", { pseudoDir });
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
  }, [qePath, crystalData]);

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
    if (!task || task.status !== "running") return;

    setOutput(task.output.join("\n") + "\n");
    setProgress(task.progress);
  }, [activeTaskId, taskContext.getTask(activeTaskId ?? "")?.output.length]);

  // Handle k-path changes from the BZ viewer
  const handleKPathChange = useCallback((newPath: KPathPoint[]) => {
    setKPath(
      normalizeKPathSampling(newPath, kPathSamplingMode, pointsPerSegment, totalKPointsTarget),
    );
  }, [kPathSamplingMode, pointsPerSegment, totalKPointsTarget]);

  useEffect(() => {
    setKPath((prevPath) =>
      normalizeKPathSampling(prevPath, kPathSamplingMode, pointsPerSegment, totalKPointsTarget),
    );
  }, [kPathSamplingMode, pointsPerSegment, totalKPointsTarget]);

  const kPathSegmentCount = useMemo(
    () => getConnectedSegmentIndices(kPath).length,
    [kPath],
  );
  const totalKPoints = useMemo(
    () => kPath.reduce((sum, point) => sum + point.npoints, 0),
    [kPath],
  );
  const minimumTotalKPoints = Math.max(1, kPathSegmentCount);
  const viewerPointsPerSegment = kPathSamplingMode === "total"
    ? clampInt(
      totalKPointsTarget / minimumTotalKPoints,
      1,
      MAX_POINTS_PER_SEGMENT,
    )
    : pointsPerSegment;

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
    const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");

    const ecutwfcValue = Number(scfParams.ecutwfc);
    const ecutwfc = Number.isFinite(ecutwfcValue) && ecutwfcValue > 0 ? ecutwfcValue : 40;
    const ecutrhoValue = Number(scfParams.ecutrho);
    const ecutrho = Number.isFinite(ecutrhoValue) && ecutrhoValue > 0
      ? ecutrhoValue
      : ecutwfc * 8;
    const sourceNspin = Number(scfParams.nspin);
    const nspin = Number.isFinite(sourceNspin) && sourceNspin > 0 ? sourceNspin : 1;

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

    const sourceUsesPrimitive = sourceScfUsesPrimitiveCell(scfParams);
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

    if (canUseSymmetryPrimitive && resolvedSymmetry) {
      baseCalculation = {
        calculation: "scf",
        prefix: scfPrefix,
        outdir: "./tmp",
        pseudo_dir: pseudoDir,
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
          nspin,
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

    const bandCellRepresentation = canUseSymmetryPrimitive
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
      mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
    };

    const pathString = kPath.map((point) => point.label).join(" → ");
    const saveParameters = {
      source_scf_id: selectedScf.id,
      k_path: pathString,
      k_path_sampling_mode: kPathSamplingMode,
      points_per_segment: kPathSamplingMode === "segment" ? pointsPerSegment : null,
      total_k_points_target: kPathSamplingMode === "total" ? totalKPoints : null,
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
    if (hasExternalRunningTask) {
      setError("Another task is currently running. Queue this task or wait for completion.");
      return;
    }
    setIsRunning(true);
    followOutputRef.current = true;
    setOutput("");
    setError(null);
    setBandData(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Band structure"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

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
          workingDir: BANDS_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save band calculation:", saveError);
        setError(`Failed to auto-save band calculation: ${saveError}`);
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
                {getCalculationTags(scf).map((tag, i) => (
                  <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
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
              <p>Choose how k-points are distributed along the selected path.</p>
            </div>
            <span className="kpath-sampling-summary">
              {totalKPoints} total k-points
            </span>
          </div>

          <div className="phonon-unit-toggle kpath-sampling-toggle" role="group" aria-label="K-path sampling mode">
            <button
              type="button"
              className={`phonon-unit-btn ${kPathSamplingMode === "segment" ? "active" : ""}`}
              onClick={() => setKPathSamplingMode("segment")}
            >
              Points / segment
            </button>
            <button
              type="button"
              className={`phonon-unit-btn ${kPathSamplingMode === "total" ? "active" : ""}`}
              onClick={() => {
                if (totalKPoints > 0) {
                  setTotalKPointsTarget(totalKPoints);
                  setTotalKPointsInput(String(totalKPoints));
                }
                setKPathSamplingMode("total");
              }}
            >
              Total k-points
            </button>
          </div>

          {kPathSamplingMode === "segment" ? (
            <label className="kpath-sampling-input">
              <span>Points per segment</span>
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
          )}

          {kPathSamplingMode === "total" && (
            <p className="kpath-sampling-note">
              {kPathSegmentCount > 0
                ? `Evenly distributed by segment length across ${kPathSegmentCount} path segment${kPathSegmentCount === 1 ? "" : "s"}.`
                : "Add at least 2 points to distribute k-points along the path."}
            </p>
          )}
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
    const bandsFilband = sanitizeOutputFilename(bandsFilbandInput, "bands.dat");
    const projectionFilproj = sanitizeOutputFilename(projectionFilprojInput, "bands.projwfc.dat");
    const canRun = isNbndValid && isConvThrValid && isMixingBetaValid && isDegaussValid;

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
        </div>

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("kpath")}>
            Back
          </button>
          <button className="secondary-button" onClick={queueCalculation} disabled={!canRun}>
            Queue Task
          </button>
          <button className="primary-button" onClick={runCalculation} disabled={!canRun || hasExternalRunningTask}>
            Run Calculation
          </button>
        </div>
      </div>
    );
  };

  // Step 4: Running
  const renderRunStep = () => {
    return (
      <div className="wizard-step run-step">
        <h3>{isRunning ? "Running Band Structure Calculation" : "Band Structure Output"}</h3>

        <ProgressBar
          status={progress.status}
          percent={progress.percent}
          phase={progress.phase}
          detail={progress.detail}
        />
        <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />

        <pre ref={outputRef} className="calculation-output" onScroll={handleOutputScroll}>
          {output || "Starting calculation..."}
        </pre>

        {error && (
          <div className="error-actions">
            <button className="secondary-button" onClick={() => setStep("parameters")}>
              Back to Parameters
            </button>
          </div>
        )}
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
    <div className="band-structure-wizard">
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
