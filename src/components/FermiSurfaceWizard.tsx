import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import {
  CrystalData,
  ELEMENT_MASSES,
  ExecutionMode,
  HpcProfile,
  SlurmResourceRequest,
} from "../lib/types";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { getPrimitiveCell, PrimitiveCell } from "../lib/primitiveCell";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import {
  buildExecutionTarget,
  buildHpcQeInputCommandLine,
  downloadHpcCalculationArtifacts,
  defaultResourcesForProfile,
  listRemotePseudopotentials,
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

interface FermiSurfaceWizardProps {
  onBack: () => void;
  onExecutionModeChange?: (mode: ExecutionMode) => Promise<void> | void;
  qePath: string;
  executionMode?: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "parameters" | "run" | "results";
const FERMI_WORK_DIR = "/tmp/qcortado_fermi_surface";

interface FrmsfFileData {
  file_name: string;
  size_bytes: number;
}

interface FermiSurfaceData {
  k_grid: [number, number, number];
  fermi_energy: number | null;
  primary_file: string;
  frmsf_files: FrmsfFileData[];
}

interface FermiSurfaceTaskPlan {
  taskLabel: string;
  taskParams: Record<string, any>;
  saveParameters: Record<string, any>;
}

function Tooltip({ text }: { text: string }) {
  return (
    <span className="tooltip-container">
      <span className="tooltip-icon">?</span>
      <span className="tooltip-text">{text}</span>
    </span>
  );
}

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

function getCalculationTags(
  calc: CalculationRun,
  downloadedIds?: Set<string>,
): { label: string; type: "info" | "feature" | "special" | "geometry" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" | "geometry" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature" | "special" | "geometry") => {
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

function getBaseElement(symbol: string): string {
  return symbol.replace(/[\d+-]+$/, "");
}

function pseudoMatchesElement(filename: string, element: string): boolean {
  const lowerFile = filename.toLowerCase();
  const lowerEl = getBaseElement(element).toLowerCase();
  return (
    lowerFile.startsWith(`${lowerEl}.`)
    || lowerFile.startsWith(`${lowerEl}_`)
    || lowerFile.startsWith(`${lowerEl}-`)
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

export function FermiSurfaceWizard({
  onBack,
  onExecutionModeChange,
  qePath,
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: FermiSurfaceWizardProps) {
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );

  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);

  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [downloadedDependencyScfIds, setDownloadedDependencyScfIds] = useState<Set<string>>(new Set());
  const [isResolvingDependency, setIsResolvingDependency] = useState(false);
  const [dependencyStatus, setDependencyStatus] = useState<string | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());
  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);

  const [fermiKGrid, setFermiKGrid] = useState<[number, number, number]>([24, 24, 24]);
  const [fermiKOffset, setFermiKOffset] = useState<[0 | 1, 0 | 1, 0 | 1]>([0, 0, 0]);
  const [fermiNbnd, setFermiNbnd] = useState<number | "auto">("auto");
  const [nscfConvThrInput, setNscfConvThrInput] = useState("1e-8");
  const [nscfMixingBetaInput, setNscfMixingBetaInput] = useState("0.7");
  const [nscfOccupations, setNscfOccupations] = useState<"fixed" | "smearing" | "from_input" | "tetrahedra">("smearing");
  const [nscfSmearing, setNscfSmearing] = useState<"gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac">("gaussian");
  const [nscfDegaussInput, setNscfDegaussInput] = useState("0.02");
  const [nscfVerbosity, setNscfVerbosity] = useState<"low" | "high" | "debug">("high");
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({
    mesh: true,
    nscf: true,
    mpi: false,
  });

  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Fermi Surface",
  });
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [fermiData, setFermiData] = useState<FermiSurfaceData | null>(null);
  const [scfFermiEnergy, setScfFermiEnergy] = useState<number | null>(null);
  const [showOutput, setShowOutput] = useState(true);
  const [isSaved, setIsSaved] = useState(false);

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

  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});
  const hpcCommandLines = useMemo(
    () => [
      "cd \"$SLURM_SUBMIT_DIR\"",
      "QE_BIN=\"$HOME/qe/bin\"",
      buildHpcQeInputCommandLine(
        activeHpcProfile,
        "fermi_velocity.x",
        "fermi_velocity.in",
        "fermi_velocity.out",
        "-npool 1",
      ),
    ],
    [activeHpcProfile],
  );

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [isHpcMode, activeHpcProfile?.id, activeHpcProfile?.resource_mode]);

  const handleOutputScroll = () => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followOutputRef.current = distanceToBottom <= 16;
  };

  const toggleSection = (section: string) => {
    setExpandedSections((prev) => ({ ...prev, [section]: !prev[section] }));
  };

  const applyScfDefaults = useCallback((scf: CalculationRun) => {
    const params = scf.parameters || {};

    if (Array.isArray(params.kgrid) && params.kgrid.length === 3) {
      const parsedGrid = params.kgrid.map((value: unknown) => Number(value));
      if (parsedGrid.every((value: number) => Number.isInteger(value) && value > 0)) {
        setFermiKGrid([parsedGrid[0], parsedGrid[1], parsedGrid[2]]);
      }
    }

    if (Array.isArray(params.kgrid_offset) && params.kgrid_offset.length === 3) {
      const parsedOffset = params.kgrid_offset.map((value: unknown) => Number(value));
      if (parsedOffset.every((value: number) => value === 0 || value === 1)) {
        setFermiKOffset([parsedOffset[0] as 0 | 1, parsedOffset[1] as 0 | 1, parsedOffset[2] as 0 | 1]);
      }
    }

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

    const sourceNbnd = Number(params.nbnd);
    if (Number.isInteger(sourceNbnd) && sourceNbnd > 0) {
      setFermiNbnd(sourceNbnd);
    } else {
      setFermiNbnd("auto");
    }

    const verbosityRaw = String(params.verbosity || "high").toLowerCase();
    if (verbosityRaw === "low" || verbosityRaw === "debug") {
      setNscfVerbosity(verbosityRaw);
    } else {
      setNscfVerbosity("high");
    }
  }, []);

  const selectSourceScf = useCallback((scf: CalculationRun) => {
    setSelectedScf(scf);
    setDependencyStatus(null);
    setScfFermiEnergy(scf.result?.fermi_energy ?? null);
    applyScfDefaults(scf);
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

        const pseudoDir = isHpcMode
          ? (activeHpcProfile?.remote_pseudo_dir || "")
          : qePath.replace(/\/bin\/?$/, "/pseudo");
        const pseudos = isHpcMode
          ? await listRemotePseudopotentials(pseudoDir, activeHpcProfile?.id ?? null)
          : await invoke<string[]>("list_pseudopotentials", { pseudoDir });
        const elements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
        const selected: Record<string, string> = {};
        for (const el of elements) {
          const match = pseudos.find((pseudo) => pseudoMatchesElement(pseudo, el));
          if (match) selected[el] = match;
        }
        setSelectedPseudos(selected);
      } catch (e) {
        console.error("Failed to initialize Fermi-surface wizard:", e);
      }
    }

    void init();
  }, [crystalData, qePath, isHpcMode, activeHpcProfile?.id, activeHpcProfile?.remote_pseudo_dir]);

  useEffect(() => {
    const el = outputRef.current;
    if (!el || !followOutputRef.current) return;
    el.scrollTop = el.scrollHeight;
  }, [output]);

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

  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) {
      void taskContext.reconnectToTask(activeTaskId);
      return;
    }

    setIsRunning(task.status === "running");
    setOutput(task.output.join("\n") + (task.output.length > 0 ? "\n" : ""));
    setProgress(task.progress);
    setCalcStartTime(task.startedAt);

    if (task.status === "completed" && task.result) {
      setFermiData(task.result as FermiSurfaceData);
      setStep("results");
    } else if (task.status === "failed" || task.status === "cancelled") {
      setError(task.error || "Task failed");
    } else {
      setStep("run");
    }
  }, [activeTaskId, taskContext, taskContext.getTask(activeTaskId ?? "")?.status]);

  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) return;

    setOutput(task.output.join("\n") + (task.output.length > 0 ? "\n" : ""));
    setProgress(task.progress);
    setIsRunning(task.status === "running");
  }, [
    activeTaskId,
    taskContext,
    taskContext.getTask(activeTaskId ?? "")?.output.length,
    taskContext.getTask(activeTaskId ?? "")?.status,
  ]);

  function buildTaskPlan(): FermiSurfaceTaskPlan {
    if (!selectedScf) {
      throw new Error("No source SCF calculation selected");
    }

    const parsedConvThr = parseOptionalPositiveNumber(nscfConvThrInput, "convergence threshold");
    if (parsedConvThr == null) {
      throw new Error("Please provide a valid positive convergence threshold.");
    }

    const parsedMixingBeta = parseOptionalPositiveNumber(nscfMixingBetaInput, "mixing beta");
    if (parsedMixingBeta == null || parsedMixingBeta > 1.0) {
      throw new Error("Mixing beta must be in the range (0, 1].");
    }

    const parsedDegauss = nscfOccupations === "smearing"
      ? parseOptionalPositiveNumber(nscfDegaussInput, "degauss")
      : null;
    if (nscfOccupations === "smearing" && parsedDegauss == null) {
      throw new Error("Please provide a positive degauss value when using smearing occupations.");
    }

    const parsedNbnd = fermiNbnd === "auto"
      ? null
      : parseOptionalPositiveInt(String(fermiNbnd), "number of bands");

    const scfParams = selectedScf.parameters || {};
    const savedPseudoMap = (scfParams.selected_pseudos && typeof scfParams.selected_pseudos === "object")
      ? scfParams.selected_pseudos as Record<string, string>
      : {};
    const elements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
    const resolvedPseudos: Record<string, string> = {};
    for (const element of elements) {
      const savedPseudo = savedPseudoMap[element];
      const detectedPseudo = selectedPseudos[element];
      const resolvedPseudo = (typeof savedPseudo === "string" && savedPseudo.length > 0)
        ? savedPseudo
        : detectedPseudo;
      if (!resolvedPseudo) {
        throw new Error(`No pseudopotential selected for element ${element}`);
      }
      resolvedPseudos[element] = resolvedPseudo;
    }

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
    const prefix = scfParams.prefix || "qcortado_scf";

    const ecutwfcValue = Number(scfParams.ecutwfc);
    const ecutwfc = Number.isFinite(ecutwfcValue) && ecutwfcValue > 0 ? ecutwfcValue : 40;
    const ecutrhoValue = Number(scfParams.ecutrho);
    const ecutrho = Number.isFinite(ecutrhoValue) && ecutrhoValue > 0
      ? ecutrhoValue
      : ecutwfc * 8;

    const primitiveCell: PrimitiveCell | null = getPrimitiveCell(crystalData);
    const commonSystemFields = {
      species: elements.map((element) => ({
        symbol: element,
        mass: ELEMENT_MASSES[element] || 1.0,
        pseudopotential: resolvedPseudos[element],
      })),
      position_units: "crystal",
      ecutwfc,
      ecutrho,
      nspin: Number(scfParams.nspin) || 1,
      occupations: nscfOccupations,
      smearing: nscfSmearing,
      degauss: parsedDegauss,
    };

    const systemConfig = primitiveCell
      ? {
        ...commonSystemFields,
        ibrav: primitiveCell.ibrav,
        celldm: [primitiveCell.celldm1, 0, 0, 0, 0, 0],
        cell_parameters: null,
        cell_units: null,
        atoms: primitiveCell.atoms.map((atom) => ({
          symbol: atom.symbol,
          position: atom.position,
          if_pos: [true, true, true] as [boolean, boolean, boolean],
        })),
      }
      : (() => {
        const gammaRadians = (crystalData.cell_angle_gamma.value * Math.PI) / 180;
        const betaRadians = (crystalData.cell_angle_beta.value * Math.PI) / 180;
        const alphaRadians = (crystalData.cell_angle_alpha.value * Math.PI) / 180;
        const c = crystalData.cell_length_c.value;
        const cVector: [number, number, number] = [
          c * Math.cos(betaRadians),
          (c * (Math.cos(alphaRadians) - Math.cos(betaRadians) * Math.cos(gammaRadians))) / Math.sin(gammaRadians),
          Math.sqrt(
            Math.max(
              0,
              c * c
                - (c * Math.cos(betaRadians)) ** 2
                - (
                  (c * (Math.cos(alphaRadians) - Math.cos(betaRadians) * Math.cos(gammaRadians)))
                  / Math.sin(gammaRadians)
                ) ** 2,
            ),
          ),
        ];

        return {
          ...commonSystemFields,
          ibrav: "free" as const,
          celldm: null,
          cell_parameters: [
            [crystalData.cell_length_a.value, 0, 0],
            [
              crystalData.cell_length_b.value * Math.cos(gammaRadians),
              crystalData.cell_length_b.value * Math.sin(gammaRadians),
              0,
            ],
            cVector,
          ] as [[number, number, number], [number, number, number], [number, number, number]],
          cell_units: "angstrom" as const,
          atoms: crystalData.atom_sites.map((site) => ({
            symbol: getBaseElement(site.type_symbol),
            position: [site.fract_x, site.fract_y, site.fract_z] as [number, number, number],
            if_pos: [true, true, true] as [boolean, boolean, boolean],
          })),
        };
      })();

    const baseCalculation = {
      calculation: "scf",
      prefix,
      outdir: "./tmp",
      pseudo_dir: pseudoDir,
      system: systemConfig,
      kpoints: { type: "gamma" },
      conv_thr: parsedConvThr,
      mixing_beta: parsedMixingBeta,
      tprnfor: false,
      tstress: false,
      forc_conv_thr: null,
      etot_conv_thr: null,
      verbosity: nscfVerbosity,
    };

    const taskLabel = `Fermi Surface - ${crystalData?.formula_sum || ""}`;
    const taskParams = {
      config: {
        base_calculation: baseCalculation,
        k_grid: fermiKGrid,
        k_offset: fermiKOffset,
        nbnd: parsedNbnd,
        project_id: projectId,
        scf_calc_id: selectedScf.id,
      },
      workingDir: FERMI_WORK_DIR,
      mpiConfig: !isHpcMode && mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      executionTarget: buildExecutionTarget(
        executionMode,
        activeHpcProfile?.id ?? null,
        isHpcMode ? hpcResources : null,
        false,
      ),
    };

    const saveParameters = {
      source_scf_id: selectedScf.id,
      fermi_k_grid: fermiKGrid,
      fermi_k_offset: fermiKOffset,
      fermi_n_bands_requested: parsedNbnd,
      nscf_conv_thr: parsedConvThr,
      nscf_mixing_beta: parsedMixingBeta,
      nscf_occupations: nscfOccupations,
      nscf_smearing: nscfSmearing,
      nscf_degauss: parsedDegauss,
      nscf_verbosity: nscfVerbosity,
      fermi_surface_tool: "fermi_velocity.x",
      n_frmsf_files: null,
      frmsf_files: [],
      primary_frmsf_file: null,
      total_frmsf_bytes: null,
      nspin: scfParams.nspin,
      lspinorb: scfParams.lspinorb,
      lda_plus_u: scfParams.lda_plus_u,
      vdw_corr: scfParams.vdw_corr,
      scf_fermi_energy: selectedScf.result?.fermi_energy ?? scfFermiEnergy,
    };

    return {
      taskLabel,
      taskParams,
      saveParameters,
    };
  }

  async function runCalculation() {
    if (selectedScfDependencyBlocked) {
      setError("Selected SCF was computed remotely and needs a full local bundle for local execution.");
      return;
    }
    if (hasExternalRunningTask) {
      setError("Another task is currently running. Queue this task or wait for completion.");
      return;
    }

    setIsRunning(true);
    followOutputRef.current = true;
    setOutput("");
    setError(null);
    setFermiData(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Fermi Surface"));
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
      const plan = buildTaskPlan();
      const taskId = await taskContext.startTask("fermi_surface", plan.taskParams, plan.taskLabel);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Calculation failed");
      }

      const result = finalTask.result as FermiSurfaceData;
      const outputContent = finalTask.output.join("\n");
      const endTime = new Date().toISOString();
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
      setFermiData(result);
      setStep("results");
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      try {
        const frmsfFiles = result.frmsf_files.map((file) => file.file_name);
        const totalFrmsfBytes = result.frmsf_files.reduce((sum, file) => sum + (Number(file.size_bytes) || 0), 0);

        await invoke("save_calculation", {
          projectId,
          cifId,
          calcData: {
            calc_type: "fermi_surface",
            parameters: {
              ...plan.saveParameters,
              n_frmsf_files: result.frmsf_files.length,
              frmsf_files: frmsfFiles,
              primary_frmsf_file: result.primary_file,
              total_frmsf_bytes: totalFrmsfBytes,
              ...hpcSaveParams,
            },
            result: {
              converged: true,
              total_energy: null,
              fermi_energy: result.fermi_energy ?? selectedScf?.result?.fermi_energy ?? scfFermiEnergy,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: outputContent,
            },
            started_at: startTime,
            completed_at: endTime,
            input_content: "",
            output_content: outputContent,
            tags: [],
          },
          workingDir: FERMI_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save Fermi-surface calculation:", saveError);
        setError(`Failed to auto-save Fermi-surface calculation: ${saveError}`);
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
  }

  function queueCalculation() {
    if (selectedScfDependencyBlocked) {
      setError("Selected SCF was computed remotely and needs a full local bundle for local execution.");
      return;
    }
    try {
      const plan = buildTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "fermi_surface",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId,
          workingDir: FERMI_WORK_DIR,
          calcType: "fermi_surface",
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        },
      );
    } catch (e) {
      setError(String(e));
    }
  }

  const renderSourceStep = () => {
    const validScfs = scfCalculations.filter((calc) => calc.calc_type === "scf" && calc.result?.converged);
    const sortedScfs = sortScfByMode(validScfs, scfSortMode);

    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Fermi-surface generation requires a completed SCF calculation.
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
            <label htmlFor="fermi-scf-sort">Sort SCFs</label>
            <select
              id="fermi-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose the SCF calculation used to run QE `fermi_velocity.x` (FermiSurfer tutorial workflow).
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
                {scf.result?.fermi_energy != null && (
                  <span>EF = {scf.result.fermi_energy.toFixed(3)} eV</span>
                )}
              </div>
              <div className="calc-tags">
                {getCalculationTags(scf, downloadedDependencyScfIds).map((tag, i) => (
                  <span
                    key={`${tag.label}-${i}`}
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
            onClick={() => setStep("parameters")}
          >
            Next: Parameters
          </button>
        </div>
      </div>
    );
  };

  const renderParametersStep = () => {
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
    const isNbndValid = fermiNbnd === "auto" || (Number.isInteger(fermiNbnd) && fermiNbnd > 0);
    const isConvThrValid = parsedConvThr !== null;
    const isMixingBetaValid = parsedMixingBeta !== null && parsedMixingBeta <= 1.0;
    const isDegaussValid = !degaussRequired || parsedDegauss !== null;
    const canRun = isNbndValid && isConvThrValid && isMixingBetaValid && isDegaussValid;

    return (
      <div className="wizard-step parameters-step">
        <h3>Fermi-Surface Parameters</h3>
        <p className="step-description">
          Configure QE `fermi_velocity.x` input controls (Section 7 tutorial workflow).
        </p>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("mesh")} className="section-header">
            <span className={`collapse-icon ${expandedSections.mesh ? "expanded" : ""}`}>▶</span>
            K-Grid & Band Count
            <Tooltip text="Controls the `K_POINTS automatic` mesh and optional `nbnd` override used to run `fermi_velocity.x`." />
          </h4>
          {expandedSections.mesh && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    K-point grid
                    <Tooltip text="Finer meshes improve Fermi-surface detail, but increase runtime and memory usage." />
                  </label>
                  <div className="qgrid-inputs">
                    <input
                      type="number"
                      min={1}
                      value={fermiKGrid[0]}
                      onChange={(e) => setFermiKGrid([Math.max(1, parseInt(e.target.value, 10) || 1), fermiKGrid[1], fermiKGrid[2]])}
                    />
                    <span>×</span>
                    <input
                      type="number"
                      min={1}
                      value={fermiKGrid[1]}
                      onChange={(e) => setFermiKGrid([fermiKGrid[0], Math.max(1, parseInt(e.target.value, 10) || 1), fermiKGrid[2]])}
                    />
                    <span>×</span>
                    <input
                      type="number"
                      min={1}
                      value={fermiKGrid[2]}
                      onChange={(e) => setFermiKGrid([fermiKGrid[0], fermiKGrid[1], Math.max(1, parseInt(e.target.value, 10) || 1)])}
                    />
                  </div>
                </div>

                <div className="phonon-field">
                  <label>
                    K-point offset (0/1)
                    <Tooltip text="Monkhorst-Pack offset per reciprocal direction. Use 1 to shift half a grid step in that direction." />
                  </label>
                  <div className="qgrid-inputs">
                    <select
                      value={fermiKOffset[0]}
                      onChange={(e) => setFermiKOffset([Number(e.target.value) as 0 | 1, fermiKOffset[1], fermiKOffset[2]])}
                    >
                      <option value={0}>0</option>
                      <option value={1}>1</option>
                    </select>
                    <span>×</span>
                    <select
                      value={fermiKOffset[1]}
                      onChange={(e) => setFermiKOffset([fermiKOffset[0], Number(e.target.value) as 0 | 1, fermiKOffset[2]])}
                    >
                      <option value={0}>0</option>
                      <option value={1}>1</option>
                    </select>
                    <span>×</span>
                    <select
                      value={fermiKOffset[2]}
                      onChange={(e) => setFermiKOffset([fermiKOffset[0], fermiKOffset[1], Number(e.target.value) as 0 | 1])}
                    >
                      <option value={0}>0</option>
                      <option value={1}>1</option>
                    </select>
                  </div>
                </div>

                <div className="phonon-field">
                  <label>
                    Number of bands
                    <Tooltip text="QE variable: `nbnd` in the pw-style input used by `fermi_velocity.x`." />
                  </label>
                  <div className="nbnd-input">
                    <select
                      value={fermiNbnd === "auto" ? "auto" : "manual"}
                      onChange={(e) => setFermiNbnd(e.target.value === "auto" ? "auto" : 20)}
                    >
                      <option value="auto">Auto (QE default)</option>
                      <option value="manual">Manual</option>
                    </select>
                    {fermiNbnd !== "auto" && (
                      <input
                        type="number"
                        min={1}
                        value={fermiNbnd}
                        onChange={(e) => setFermiNbnd(Math.max(1, parseInt(e.target.value, 10) || 1))}
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
            Electronic Controls
            <Tooltip text="Advanced controls for the generated pw-style input: `conv_thr`, `mixing_beta`, occupations/smearing/degauss, and verbosity." />
          </h4>
          {expandedSections.nscf && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    Convergence threshold
                    <Tooltip text="QE variable: `conv_thr` in `&ELECTRONS`." />
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
                    <Tooltip text="QE variable: `mixing_beta`. Stable values are typically in (0, 1]." />
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
                    <Tooltip text="QE variable: `occupations`. Metals typically require smearing." />
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
                    <Tooltip text="QE variable: `smearing`. Used only when occupations = smearing." />
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
                    <Tooltip text="QE variable: `degauss`. Required with smearing occupations." />
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
                    <Tooltip text="QE variable: `verbosity` in `&CONTROL`." />
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

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>Source SCF:</span>
            <span>{selectedScf?.id.slice(0, 8) || "N/A"}</span>
          </div>
          <div className="summary-row">
            <span>K-grid:</span>
            <span>{fermiKGrid[0]}×{fermiKGrid[1]}×{fermiKGrid[2]}</span>
          </div>
          <div className="summary-row">
            <span>K-offset:</span>
            <span>{fermiKOffset[0]} {fermiKOffset[1]} {fermiKOffset[2]}</span>
          </div>
          <div className="summary-row">
            <span>Requested `nbnd`:</span>
            <span>{fermiNbnd === "auto" ? "Auto" : fermiNbnd}</span>
          </div>
          <div className="summary-row">
            <span>Expected output:</span>
            <span>`vfermi.frmsf` (or spin-resolved `vfermi1/2.frmsf`)</span>
          </div>
        </div>

        {selectedScfDependencyBlocked && (
          <div className="warning-banner">
            <strong>Remote SCF dependency detected.</strong>{" "}
            This source SCF was run on HPC, and local Fermi-surface requires a full local SCF bundle.
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

        {isHpcMode ? (
          <HpcRunSettings
            profileId={activeHpcProfile?.id ?? null}
            profileName={activeHpcProfile?.name ?? "Andromeda"}
            taskKind="fermi_surface"
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
              <Tooltip text="Process-level MPI controls for `fermi_velocity.x` (`npool` is fixed to 1 per tutorial guidance)." />
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

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button
            className="secondary-button"
            onClick={() => void queueCalculation()}
            disabled={!canRun || selectedScfDependencyBlocked || isResolvingDependency}
          >
            Queue Task
          </button>
          <button
            className="primary-button"
            onClick={() => void runCalculation()}
            disabled={!canRun || hasExternalRunningTask || selectedScfDependencyBlocked || isResolvingDependency}
          >
            {isHpcMode ? "Submit Fermi Surface to Andromeda" : "Run Calculation"}
          </button>
        </div>
      </div>
    );
  };

  const renderRunStep = () => (
    <div className="wizard-step run-step run-step-focused">
      <div className="run-step-headline">
        <h3>{isRunning ? "Running Fermi-Surface Generation" : "Fermi-Surface Output"}</h3>
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

      <div className={`run-layout ${isHpcMode ? "run-layout-hpc-telemetry" : ""}`}>
        <div className="output-panel">
          <h3>{isRunning ? "Running..." : "Output"}</h3>
          <pre ref={outputRef} className="output-text" onScroll={handleOutputScroll}>
            {output || "Starting calculation..."}
          </pre>
        </div>

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

  const renderResultsStep = () => {
    if (!fermiData) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("parameters")}>
            Try Again
          </button>
        </div>
      );
    }

    const totalBytes = fermiData.frmsf_files.reduce((sum, file) => sum + (Number(file.size_bytes) || 0), 0);

    return (
      <div className="wizard-step results-step">
        <h3>Fermi Surface Results</h3>
        <p className="step-description">
          FRMSF generation complete.
        </p>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Dense k-grid:</span>
              <span className="value">{fermiData.k_grid[0]}×{fermiData.k_grid[1]}×{fermiData.k_grid[2]}</span>
            </div>
            <div className="summary-item">
              <span className="label">FRMSF Files:</span>
              <span className="value">{fermiData.frmsf_files.length}</span>
            </div>
            <div className="summary-item">
              <span className="label">Primary File:</span>
              <span className="value">{fermiData.primary_file}</span>
            </div>
            <div className="summary-item">
              <span className="label">Total Size:</span>
              <span className="value">{formatBytes(totalBytes)}</span>
            </div>
            <div className="summary-item">
              <span className="label">Fermi Energy:</span>
              <span className="value">
                {fermiData.fermi_energy != null ? `${fermiData.fermi_energy.toFixed(3)} eV` : "Unavailable"}
              </span>
            </div>
          </div>
        </div>

        <div className="detail-item parameters">
          <label>Generated FRMSF Metadata</label>
          <pre>{JSON.stringify(fermiData.frmsf_files, null, 2)}</pre>
        </div>

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
          <button className="primary-button" onClick={() => setStep("parameters")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className={`electronic-dos-wizard fermi-surface-wizard wizard-step-${step}`}>
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Fermi Surface Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "parameters" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            2. Parameters
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "completed" : ""}>
            3. Run
          </span>
          <span className={step === "results" ? "active" : ""}>
            4. Results
          </span>
        </div>
      </div>

      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "parameters" && renderParametersStep()}
        {step === "run" && renderRunStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
