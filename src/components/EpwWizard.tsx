import { useEffect, useMemo, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { AppHeaderPortal } from "./AppHeaderPortal";
import {
  CrystalData,
  ExecutionMode,
  HpcProfile,
  SlurmResourceRequest,
} from "../lib/types";
import type { EpwViewerData } from "./EpwViewer";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { InfoTooltip } from "./InfoTooltip";
import { defaultProgressState, ProgressState } from "../lib/engines/qe/progress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import {
  buildExecutionTarget,
  buildHpcLauncherCommand,
  defaultResourcesForProfile,
} from "../lib/hpcConfig";
import {
  buildHpcQeRuntimeSetupLines,
  qeEngineUsesModuleMode,
} from "../lib/engines/qe/hpc";
import { HpcRunSettings } from "./HpcRunSettings";
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
    raw_output?: string | null;
    epw_data?: EpwCalculationResult | null;
    wannier_data?: any;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

interface EpwArtifactManifestEntry {
  source_calc_id: string;
  source_calc_type: string;
  rel_path: string;
  size_bytes: number;
}

interface EpwPrerequisiteValidation {
  ok: boolean;
  errors: string[];
  warnings: string[];
  remediation_hints: string[];
  manifests: EpwArtifactManifestEntry[];
}

interface EpwSourceRef {
  calc_id: string;
  calc_type: string;
}

interface EpwGoalSelection {
  coupling: boolean;
  phonon_linewidth_a2f: boolean;
  electron_self_energy: boolean;
  transport_mobility: boolean;
  superconductivity: boolean;
}

interface EpwCalculationResult extends EpwViewerData {
  schema_version: number;
  sources: {
    phonon: EpwSourceRef;
    wannier?: EpwSourceRef | null;
    scf?: EpwSourceRef | null;
    manifests: EpwArtifactManifestEntry[];
  };
  input: {
    prefix: string;
    outdir: string;
    dvscf_dir: string;
    wannier_dir: string;
    k_mesh: [number, number, number];
    q_mesh: [number, number, number];
    coarse_k_mesh?: [number, number, number] | null;
    fine_k_mesh?: [number, number, number] | null;
    coarse_q_mesh?: [number, number, number] | null;
    fine_q_mesh?: [number, number, number] | null;
    epbwrite: boolean;
    epbread: boolean;
    epwwrite: boolean;
    epwread: boolean;
    wannierize: boolean;
    fsthick_ev?: number | null;
    degaussw_ev?: number | null;
    nbndsub?: number | null;
  };
  goals?: EpwGoalSelection;
  runtime: {
    pools?: number | null;
    max_seconds?: number | null;
    artifact_sync_mode?: string | null;
  };
  artifacts: Array<{ file_name: string; size_bytes: number }>;
  result_summary: {
    completed: boolean;
    elapsed_seconds?: number | null;
    core_metrics: Record<string, number>;
    generated_outputs: Array<{ file_name: string; size_bytes: number }>;
    unknown_metrics: Record<string, unknown>;
    parse_partial: boolean;
    notes: string[];
    warnings?: string[];
    parse_coverage?: string[];
  };
  errors: Array<{ code: string; message: string; hint?: string | null }>;
}

interface EpwWizardProps {
  onBack: () => void;
  onViewEPW: (epwData: EpwCalculationResult, rawOutput?: string | null) => void;
  qePath: string;
  executionMode?: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  calculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "controls" | "run" | "results";

const EPW_WORK_DIR = "/tmp/qcortado_epw";
const EPW_WIZARD_SETTINGS_ID = "epw";

interface StoredEpwWizardSettings {
  prefix: string;
  outdir: string;
  dvscfDir: string;
  wannierDir: string;
  kMesh: [number, number, number];
  qMesh: [number, number, number];
  fineQMesh: [number, number, number];
  goals: EpwGoalSelection;
  epbwrite: boolean;
  epbread: boolean;
  epwwrite: boolean;
  epwread: boolean;
  wannierize: boolean;
  fsthickInput: string;
  degausswInput: string;
  nbndsubInput: string;
  mode: string;
  runtimePoolsInput: string;
  runtimeMaxSecondsInput: string;
  artifactSyncMode: "minimal" | "epw-ready" | "full";
  advancedOverridesText: string;
}

function parseOptionalPositiveNumber(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (!trimmed) return null;
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive number.`);
  }
  return parsed;
}

function parseOptionalPositiveInt(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (!trimmed) return null;
  const parsed = Number.parseInt(trimmed, 10);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive integer.`);
  }
  return parsed;
}

function parseOverrides(text: string): Record<string, string> {
  const overrides: Record<string, string> = {};
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line || line.startsWith("#") || line.startsWith("!")) {
      continue;
    }
    const [rawKey, ...rest] = line.split("=");
    const key = rawKey?.trim().toLowerCase();
    if (!key) {
      continue;
    }
    const value = rest.join("=").trim();
    overrides[key] = value;
  }
  return overrides;
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

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value.toFixed(value >= 10 || unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
}

function prettyOverridesPreview(overrides: Record<string, string>): string {
  const entries = Object.entries(overrides).sort(([left], [right]) => left.localeCompare(right));
  if (entries.length === 0) {
    return "No advanced overrides.";
  }
  return entries.map(([key, value]) => `${key} = ${value}`).join("\n");
}

function extractQMeshCompatibilityMessage(validation: EpwPrerequisiteValidation | null): { text: string; isError: boolean } | null {
  if (!validation) return null;
  const qMeshRegex = /\bq-?mesh\b/i;
  const qGridRegex = /\bq-?grid\b/i;

  const errorMessage = validation.errors.find((entry) => qMeshRegex.test(entry) && qGridRegex.test(entry));
  if (errorMessage) {
    return { text: errorMessage, isError: true };
  }

  const warningMessage = validation.warnings.find((entry) => qMeshRegex.test(entry) && qGridRegex.test(entry));
  if (warningMessage) {
    return { text: warningMessage, isError: false };
  }

  return null;
}

function extractKMeshCompatibilityMessage(validation: EpwPrerequisiteValidation | null): { text: string; isError: boolean } | null {
  if (!validation) return null;
  const kMeshRegex = /\bk-?mesh\b/i;
  const kGridRegex = /\bk-?grid\b/i;

  const errorMessage = validation.errors.find((entry) => kMeshRegex.test(entry) && kGridRegex.test(entry));
  if (errorMessage) {
    return { text: errorMessage, isError: true };
  }

  const warningMessage = validation.warnings.find((entry) => kMeshRegex.test(entry) && kGridRegex.test(entry));
  if (warningMessage) {
    return { text: warningMessage, isError: false };
  }

  return null;
}

function extractSavedPhononQGrid(calc: CalculationRun | null): [number, number, number] | null {
  const rawGrid = calc?.parameters?.q_grid;
  if (!Array.isArray(rawGrid) || rawGrid.length !== 3) {
    return null;
  }
  const parsed = rawGrid.map((entry) => Number.parseInt(String(entry), 10));
  if (parsed.some((value) => !Number.isFinite(value) || value <= 0)) {
    return null;
  }
  return [parsed[0], parsed[1], parsed[2]];
}

function extractSavedWannierKGrid(calc: CalculationRun | null): [number, number, number] | null {
  const rawGrid = calc?.parameters?.k_grid;
  if (!Array.isArray(rawGrid) || rawGrid.length !== 3) {
    return null;
  }
  const parsed = rawGrid.map((entry) => Number.parseInt(String(entry), 10));
  if (parsed.some((value) => !Number.isFinite(value) || value <= 0)) {
    return null;
  }
  return [parsed[0], parsed[1], parsed[2]];
}

export function EpwWizard({
  onBack,
  onViewEPW,
  qePath: _qePath,
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId,
  crystalData,
  calculations,
  reconnectTaskId,
}: EpwWizardProps) {
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const storedWizardSettings = useMemo(
    () => readProjectWizardSettings<StoredEpwWizardSettings>(projectId, EPW_WIZARD_SETTINGS_ID),
    [projectId],
  );
  const shouldPreserveStoredSettingsRef = useRef(Boolean(storedWizardSettings));
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  const hasBlockingExternalTask = !isHpcMode && hasExternalRunningTask;

  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);

  const [selectedPhononId, setSelectedPhononId] = useState<string>("");
  const [selectedWannierId, setSelectedWannierId] = useState<string>("");

  const [prefix, setPrefix] = useState(() => storedWizardSettings?.prefix ?? "qcortado_scf");
  const [outdir, setOutdir] = useState(() => storedWizardSettings?.outdir ?? "./tmp");
  const [dvscfDir, setDvscfDir] = useState(() => storedWizardSettings?.dvscfDir ?? "./save");
  const [wannierDir, setWannierDir] = useState(() => storedWizardSettings?.wannierDir ?? "./wannier");
  const [kMesh, setKMesh] = useState<[number, number, number]>(() => storedWizardSettings?.kMesh ?? [24, 24, 24]);
  const [qMesh, setQMesh] = useState<[number, number, number]>(() => storedWizardSettings?.qMesh ?? [6, 6, 6]);
  const [fineQMesh, setFineQMesh] = useState<[number, number, number]>(() => storedWizardSettings?.fineQMesh ?? [24, 24, 24]);
  const [goals, setGoals] = useState<EpwGoalSelection>(() => storedWizardSettings?.goals ?? {
    coupling: true,
    phonon_linewidth_a2f: true,
    electron_self_energy: false,
    transport_mobility: false,
    superconductivity: false,
  });
  const [epbwrite, setEpbwrite] = useState(() => storedWizardSettings?.epbwrite ?? true);
  const [epbread, setEpbread] = useState(() => storedWizardSettings?.epbread ?? false);
  const [epwwrite, setEpwwrite] = useState(() => storedWizardSettings?.epwwrite ?? true);
  const [epwread, setEpwread] = useState(() => storedWizardSettings?.epwread ?? false);
  const [wannierize, setWannierize] = useState(() => storedWizardSettings?.wannierize ?? true);
  const [fsthickInput, setFsthickInput] = useState(() => storedWizardSettings?.fsthickInput ?? "0.4");
  const [degausswInput, setDegausswInput] = useState(() => storedWizardSettings?.degausswInput ?? "0.02");
  const [nbndsubInput, setNbndsubInput] = useState(() => storedWizardSettings?.nbndsubInput ?? "");

  const [mode, setMode] = useState(() => storedWizardSettings?.mode ?? "standard");
  const [runtimePoolsInput, setRuntimePoolsInput] = useState(() => storedWizardSettings?.runtimePoolsInput ?? "");
  const [runtimeMaxSecondsInput, setRuntimeMaxSecondsInput] = useState(() => storedWizardSettings?.runtimeMaxSecondsInput ?? "");
  const [artifactSyncMode, setArtifactSyncMode] = useState<"minimal" | "epw-ready" | "full">(() => storedWizardSettings?.artifactSyncMode ?? "epw-ready");

  const [advancedOverridesText, setAdvancedOverridesText] = useState(() => storedWizardSettings?.advancedOverridesText ?? "");

  const [prereqValidation, setPrereqValidation] = useState<EpwPrerequisiteValidation | null>(null);
  const [isValidatingPrerequisites, setIsValidatingPrerequisites] = useState(false);

  const [isRunning, setIsRunning] = useState(false);
  const [progress, setProgress] = useState<ProgressState>(defaultProgressState("EPW"));
  const [output, setOutput] = useState("");
  const [outputLineCount, setOutputLineCount] = useState(0);
  const visibleOutputLineCount = useMemo(() => countVisibleOutputLines(output), [output]);
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [epwResult, setEpwResult] = useState<EpwCalculationResult | null>(null);
  const [isSaved, setIsSaved] = useState(false);

  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(
    defaultResourcesForProfile(activeHpcProfile),
  );

  useViewportScrollLock(step === "run");

  const phononCalculations = useMemo(
    () => calculations
      .filter((calc) => calc.calc_type === "phonon")
      .sort((left, right) => (right.completed_at || right.started_at).localeCompare(left.completed_at || left.started_at)),
    [calculations],
  );

  const wannierCalculations = useMemo(
    () => calculations
      .filter((calc) => calc.calc_type === "wannier")
      .sort((left, right) => (right.completed_at || right.started_at).localeCompare(left.completed_at || left.started_at)),
    [calculations],
  );

  const selectedPhononCalculation = useMemo(
    () => phononCalculations.find((calc) => calc.id === selectedPhononId) ?? null,
    [phononCalculations, selectedPhononId],
  );

  const selectedWannierCalculation = useMemo(
    () => wannierCalculations.find((calc) => calc.id === selectedWannierId) ?? null,
    [wannierCalculations, selectedWannierId],
  );

  const hpcCommandLines = useMemo(
    () => {
      const remoteEpw = (activeHpcProfile?.remote_epw_path || "").trim();
      const epwExecutable = qeEngineUsesModuleMode(activeHpcProfile)
        ? "epw.x"
        : remoteEpw.length > 0 ? `"${remoteEpw}"` : "\"$QE_BIN/epw.x\"";
      return [
        "cd \"$SLURM_SUBMIT_DIR\"",
        ...buildHpcQeRuntimeSetupLines(activeHpcProfile, hpcResources.resource_type),
        `${buildHpcLauncherCommand(activeHpcProfile, hpcResources.resource_type)} ${epwExecutable} -in epw.in > epw.out 2> epw.err`,
      ];
    },
    [activeHpcProfile, hpcResources.resource_type],
  );

  const overrideMap = useMemo(
    () => parseOverrides(advancedOverridesText),
    [advancedOverridesText],
  );
  const qMeshCompatibilityMessage = useMemo(
    () => extractQMeshCompatibilityMessage(prereqValidation),
    [prereqValidation],
  );
  const kMeshCompatibilityMessage = useMemo(
    () => extractKMeshCompatibilityMessage(prereqValidation),
    [prereqValidation],
  );

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.resource_mode, isHpcMode]);

  useEffect(() => {
    writeProjectWizardSettings(projectId, EPW_WIZARD_SETTINGS_ID, {
      prefix,
      outdir,
      dvscfDir,
      wannierDir,
      kMesh,
      qMesh,
      fineQMesh,
      goals,
      epbwrite,
      epbread,
      epwwrite,
      epwread,
      wannierize,
      fsthickInput,
      degausswInput,
      nbndsubInput,
      mode,
      runtimePoolsInput,
      runtimeMaxSecondsInput,
      artifactSyncMode,
      advancedOverridesText,
    });
  }, [
    projectId,
    prefix,
    outdir,
    dvscfDir,
    wannierDir,
    kMesh,
    qMesh,
    fineQMesh,
    goals,
    epbwrite,
    epbread,
    epwwrite,
    epwread,
    wannierize,
    fsthickInput,
    degausswInput,
    nbndsubInput,
    mode,
    runtimePoolsInput,
    runtimeMaxSecondsInput,
    artifactSyncMode,
    advancedOverridesText,
  ]);

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
      } catch (initError) {
        console.error("Failed to initialize EPW wizard:", initError);
      }
    }

    void init();
  }, []);

  useEffect(() => {
    if (step !== "controls") return;
    if (!selectedPhononId) {
      setPrereqValidation(null);
      return;
    }
    void validatePrerequisites(true);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [step, selectedPhononId, selectedWannierId, qMesh[0], qMesh[1], qMesh[2], kMesh[0], kMesh[1], kMesh[2], fineQMesh[0], fineQMesh[1], fineQMesh[2]]);

  useEffect(() => {
    setFineQMesh(kMesh);
  }, [kMesh[0], kMesh[1], kMesh[2]]);

  useEffect(() => {
    if (!activeTask) return;

    setProgress(activeTask.progress);
    const outputText = activeTask.output.join("\n");
    setOutput(outputText);
    setOutputLineCount(activeTask.outputLineCount || countVisibleOutputLines(outputText));
    if (activeTask.startedAt) {
      setCalcStartTime(activeTask.startedAt);
    }

    if (activeTask.status === "running") {
      setIsRunning(true);
      setStep("run");
      return;
    }

    setIsRunning(false);
    if (activeTask.status === "completed" && activeTask.result) {
      setEpwResult(activeTask.result as EpwCalculationResult);
      setStep("results");
      setProgress({ status: "complete", percent: 100, phase: "Complete" });
    } else if (activeTask.status === "failed" || activeTask.status === "cancelled") {
      setError(activeTask.error || "EPW run failed");
      setStep("run");
      setProgress((prev) => ({
        ...prev,
        status: "error",
        phase: activeTask.status === "cancelled" ? "Cancelled" : "Failed",
      }));
    }
  }, [activeTask]);

  function handleSelectPhononSource(calc: CalculationRun) {
    setSelectedPhononId(calc.id);
    const savedQGrid = extractSavedPhononQGrid(calc);
    if (savedQGrid && !shouldPreserveStoredSettingsRef.current) {
      setQMesh(savedQGrid);
    }
    setError(null);
  }

  function handleSelectWannierSource(calc: CalculationRun | null) {
    if (!calc) {
      setSelectedWannierId("");
      setError(null);
      return;
    }
    setSelectedWannierId(calc.id);
    const savedKGrid = extractSavedWannierKGrid(calc);
    if (savedKGrid && !shouldPreserveStoredSettingsRef.current) {
      setKMesh(savedKGrid);
    }
    setError(null);
  }

  useEffect(() => {
    if (visibleOutputLineCount > outputLineCount) {
      setOutputLineCount(visibleOutputLineCount);
    }
  }, [outputLineCount, visibleOutputLineCount]);

  function updateMesh(
    kind: "k" | "q" | "fine-q",
    index: 0 | 1 | 2,
    value: string,
  ) {
    const parsed = Number.parseInt(value, 10);
    const safe = Number.isFinite(parsed) && parsed > 0 ? parsed : 1;
    if (kind === "k") {
      const next: [number, number, number] = [...kMesh] as [number, number, number];
      next[index] = safe;
      setKMesh(next);
    } else if (kind === "q") {
      const next: [number, number, number] = [...qMesh] as [number, number, number];
      next[index] = safe;
      setQMesh(next);
    } else {
      const next: [number, number, number] = [...fineQMesh] as [number, number, number];
      next[index] = safe;
      setFineQMesh(next);
    }
  }

  function updateGoal(key: keyof EpwGoalSelection, checked: boolean) {
    setGoals((prev) => ({ ...prev, [key]: checked }));
  }

  function buildConfig() {
    if (!selectedPhononCalculation) {
      throw new Error("Select a phonon source calculation.");
    }
    const coarseQGrid = extractSavedPhononQGrid(selectedPhononCalculation);
    if (!coarseQGrid) {
      throw new Error("Selected phonon source does not include a saved q-grid. Re-run phonons before EPW.");
    }
    const coarseKGrid = extractSavedWannierKGrid(selectedWannierCalculation) ?? kMesh;

    const fsthick = parseOptionalPositiveNumber(fsthickInput, "fsthick");
    const degaussw = parseOptionalPositiveNumber(degausswInput, "degaussw");
    const nbndsub = parseOptionalPositiveInt(nbndsubInput, "nbndsub");
    const pools = parseOptionalPositiveInt(runtimePoolsInput, "MPI pools");
    const maxSeconds = parseOptionalPositiveInt(runtimeMaxSecondsInput, "Max seconds");

    return {
      project_id: projectId,
      source_phonon_calc_id: selectedPhononCalculation.id,
      source_wannier_calc_id: selectedWannierCalculation?.id || null,
      source_scf_calc_id: null,
      mode,
      goals,
      input: {
        prefix: prefix.trim() || "qcortado_scf",
        outdir: outdir.trim() || "./tmp",
        dvscf_dir: dvscfDir.trim() || "./save",
        wannier_dir: wannierDir.trim() || "./wannier",
        k_mesh: coarseKGrid,
        q_mesh: coarseQGrid,
        coarse_k_mesh: coarseKGrid,
        fine_k_mesh: kMesh,
        coarse_q_mesh: coarseQGrid,
        fine_q_mesh: fineQMesh,
        epbwrite,
        epbread,
        epwwrite,
        epwread,
        wannierize,
        fsthick_ev: fsthick,
        degaussw_ev: degaussw,
        nbndsub,
      },
      runtime: {
        pools,
        max_seconds: maxSeconds,
        artifact_sync_mode: artifactSyncMode,
      },
      advanced_overrides: overrideMap,
    };
  }

  function buildTaskPlan() {
    const config = buildConfig();
    const taskLabel = `EPW - ${crystalData?.formula_sum || ""}`;
    const taskParams = {
      config,
      workingDir: EPW_WORK_DIR,
      mpiConfig: !isHpcMode && mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      executionTarget: buildExecutionTarget(
        executionMode,
        activeHpcProfile?.id ?? null,
        isHpcMode ? hpcResources : null,
        false,
      ),
    };

    const saveParameters = {
      mode,
      source_phonon_calc_id: config.source_phonon_calc_id,
      source_wannier_calc_id: config.source_wannier_calc_id,
      epw_goals: goals,
      prefix: config.input.prefix,
      outdir: config.input.outdir,
      dvscf_dir: config.input.dvscf_dir,
      wannier_dir: config.input.wannier_dir,
      k_mesh: config.input.fine_k_mesh,
      q_mesh: config.input.coarse_q_mesh,
      coarse_k_grid: config.input.coarse_k_mesh,
      fine_k_grid: config.input.fine_k_mesh,
      coarse_q_grid: config.input.coarse_q_mesh,
      fine_q_grid: config.input.fine_q_mesh,
      epbwrite: config.input.epbwrite,
      epbread: config.input.epbread,
      epwwrite: config.input.epwwrite,
      epwread: config.input.epwread,
      wannierize: config.input.wannierize,
      fsthick_ev: config.input.fsthick_ev,
      degaussw_ev: config.input.degaussw_ev,
      nbndsub: config.input.nbndsub,
      runtime_pools: config.runtime.pools,
      runtime_max_seconds: config.runtime.max_seconds,
      artifact_sync_mode: config.runtime.artifact_sync_mode,
      advanced_override_count: Object.keys(config.advanced_overrides || {}).length,
      parse_partial: null,
      artifact_count: null,
    };

    return {
      config,
      taskLabel,
      taskParams,
      saveParameters,
    };
  }

  async function validatePrerequisites(silent = false) {
    if (!selectedPhononId) {
      if (!silent) setError("Select a phonon source first.");
      setPrereqValidation(null);
      return;
    }

    setIsValidatingPrerequisites(true);
    if (!silent) setError(null);
    try {
      const config = buildConfig();
      const validation = await invoke<EpwPrerequisiteValidation>("validate_epw_prerequisites", {
        config,
      });
      setPrereqValidation(validation);
    } catch (validationError) {
      if (!silent) {
        setError(String(validationError));
      } else {
        console.warn("EPW prerequisite validation refresh failed:", validationError);
      }
      setPrereqValidation(null);
    } finally {
      setIsValidatingPrerequisites(false);
    }
  }

  async function runCalculation() {
    if (hasBlockingExternalTask) {
      setError("Another local task is currently running. Queue this task or wait for completion.");
      return;
    }

    setError(null);
    setIsRunning(true);
    setOutput("");
    setOutputLineCount(0);
    setEpwResult(null);
    setIsSaved(false);
    setProgress(defaultProgressState("EPW"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

    try {
      const plan = buildTaskPlan();

      const validation = await invoke<EpwPrerequisiteValidation>("validate_epw_prerequisites", {
        config: plan.config,
      });
      setPrereqValidation(validation);
      if (!validation.ok) {
        throw new Error("EPW prerequisite validation failed. Resolve the reported issues before running.");
      }

      const hpcSaveSpec = isHpcMode
        ? {
          projectId,
          cifId,
          workingDir: EPW_WORK_DIR,
          calcType: "epw" as const,
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        }
        : null;
      const taskId = await taskContext.startTask("epw", plan.taskParams, plan.taskLabel, hpcSaveSpec);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "EPW calculation failed");
      }

      const result = finalTask.result as EpwCalculationResult;
      const outputContent = finalTask.output.join("\n");
      setOutput(outputContent);
      setOutputLineCount(countVisibleOutputLines(outputContent));
      const completedAt = new Date().toISOString();
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

      setEpwResult(result);
      setProgress({ status: "complete", percent: 100, phase: "Complete" });
      setStep("results");

      try {
        await invoke("save_calculation", {
          projectId,
          cifId,
          calcData: {
            calc_type: "epw",
            parameters: {
              ...plan.saveParameters,
              parse_partial: result?.result_summary?.parse_partial ?? null,
              artifact_count: Array.isArray(result?.artifacts) ? result.artifacts.length : null,
              ...hpcSaveParams,
            },
            result: {
              converged: Boolean(result?.result_summary?.completed),
              total_energy: null,
              fermi_energy: null,
              n_scf_steps: null,
              wall_time_seconds: result?.result_summary?.elapsed_seconds ?? null,
              raw_output: outputContent,
              epw_data: result,
            },
            started_at: startTime,
            completed_at: completedAt,
            input_content: "",
            output_content: outputContent,
            tags: [],
          },
          workingDir: finalTask.hpc.local_sync_dir ?? EPW_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save EPW calculation:", saveError);
        setError(`Failed to auto-save EPW calculation: ${saveError}`);
      }
    } catch (runError) {
      setError(String(runError));
      setProgress({
        status: "error",
        percent: null,
        phase: "Error",
      });
    } finally {
      setIsRunning(false);
    }
  }

  function queueCalculation() {
    if (isHpcMode) {
      setError("Queueing is unavailable in HPC mode. Submit directly to HPC.");
      return;
    }

    try {
      const plan = buildTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "epw",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId,
          workingDir: EPW_WORK_DIR,
          calcType: "epw",
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        },
      );
      setStep("run");
    } catch (queueError) {
      setError(String(queueError));
    }
  }

  function renderSourceStep() {
    if (phononCalculations.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No Phonon Calculations Available</h3>
          <p className="warning-text">
            EPW requires at least one completed phonon calculation with retained artifacts.
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
          <h3>Select EPW Sources</h3>
        </div>
        <p className="step-description">
          Pick a saved phonon calculation and, optionally, a compatible Wannier calculation.
        </p>

        <div className="option-section">
          <h4>
            Phonon Source
            <InfoTooltip text="Choose the completed phonon run that provides dyn*/dvscf prerequisites for EPW staging." />
          </h4>
          <div className="scf-list">
            {phononCalculations.map((calc) => {
              const isSelected = selectedPhononId === calc.id;
              const calcName = getCalculationName(calc);
              return (
                <div
                  key={calc.id}
                  className={`scf-option ${isSelected ? "selected" : ""}`}
                  onClick={() => handleSelectPhononSource(calc)}
                >
                  <div className="scf-option-header">
                    <input
                      type="radio"
                      checked={isSelected}
                      onChange={() => handleSelectPhononSource(calc)}
                    />
                    {calcName && (
                      <span className="scf-name" title={formatCalculationSourceLabel(calc, "Phonon")}>
                        {calcName}
                      </span>
                    )}
                    <span className="scf-date">
                      {new Date(calc.started_at).toLocaleDateString()}
                    </span>
                  </div>
                  <div className="scf-details">
                    <span>ID = {calc.id.slice(0, 8)}</span>
                    <span>{calc.completed_at ? "Completed" : "In progress"}</span>
                  </div>
                  <div className="calc-tags">
                    {isHpcCalculation(calc) && <span className="calc-tag calc-tag-feature calc-tag-hpc">HPC</span>}
                    {calc.result?.converged && <span className="calc-tag calc-tag-special">Converged</span>}
                  </div>
                </div>
              );
            })}
          </div>
        </div>

        <div className="option-section">
          <h4>
            Wannier Source (Optional)
            <InfoTooltip text="Optional reusable Wannier outputs. Leave as none when running EPW without an upstream Wannier stage." />
          </h4>
          <div className="scf-list">
            <div
              className={`scf-option ${selectedWannierId === "" ? "selected" : ""}`}
              onClick={() => handleSelectWannierSource(null)}
            >
              <div className="scf-option-header">
                <input
                  type="radio"
                  checked={selectedWannierId === ""}
                  onChange={() => handleSelectWannierSource(null)}
                />
                <span className="scf-date">Optional</span>
              </div>
              <div className="scf-details">
                <span>No Wannier source</span>
              </div>
            </div>

            {wannierCalculations.map((calc) => {
              const isSelected = selectedWannierId === calc.id;
              const calcName = getCalculationName(calc);
              return (
                <div
                  key={calc.id}
                  className={`scf-option ${isSelected ? "selected" : ""}`}
                  onClick={() => handleSelectWannierSource(calc)}
                >
                  <div className="scf-option-header">
                    <input
                      type="radio"
                      checked={isSelected}
                      onChange={() => handleSelectWannierSource(calc)}
                    />
                    {calcName && (
                      <span className="scf-name" title={formatCalculationSourceLabel(calc, "Wannier")}>
                        {calcName}
                      </span>
                    )}
                    <span className="scf-date">
                      {new Date(calc.started_at).toLocaleDateString()}
                    </span>
                  </div>
                  <div className="scf-details">
                    <span>ID = {calc.id.slice(0, 8)}</span>
                    <span>{calc.completed_at ? "Completed" : "In progress"}</span>
                  </div>
                  <div className="calc-tags">
                    {isHpcCalculation(calc) && <span className="calc-tag calc-tag-feature calc-tag-hpc">HPC</span>}
                    {calc.result?.converged && <span className="calc-tag calc-tag-special">Converged</span>}
                  </div>
                </div>
              );
            })}
          </div>
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>Cancel</button>
          <button
            className="primary-button"
            disabled={!selectedPhononId}
            onClick={() => setStep("controls")}
          >
            Next: Controls
          </button>
        </div>
      </div>
    );
  }

  function renderControlsStep() {
    const selectedQGrid = selectedPhononCalculation ? extractSavedPhononQGrid(selectedPhononCalculation) : null;
    const selectedCoarseKGrid = extractSavedWannierKGrid(selectedWannierCalculation);

    return (
      <div className="wizard-step">
        <h3>EPW Controls</h3>
        <p className="step-description">
          Choose the EPW result families and runtime scale. Source-compatible settings are inferred from the selected calculations.
        </p>

        <div className="option-section epw-goals-section">
          <h4>
            Results
            <InfoTooltip text="QCortado maps these goals to the EPW flags needed for the selected calculation." />
          </h4>
          <div className="epw-goals-grid">
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={goals.coupling}
                onChange={(event) => updateGoal("coupling", event.target.checked)}
              />
              <span>
                Electron-phonon coupling
                <InfoTooltip text="Enables EPW electron-phonon coupling setup and reusable matrix outputs for downstream analyses." />
              </span>
            </label>
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={goals.phonon_linewidth_a2f}
                onChange={(event) => updateGoal("phonon_linewidth_a2f", event.target.checked)}
              />
              <span>
                Phonon linewidths / alpha2F
                <InfoTooltip text="Requests phonon self-energies and Eliashberg spectral function output (`phonselfen` and `a2f`)." />
              </span>
            </label>
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={goals.electron_self_energy}
                onChange={(event) => updateGoal("electron_self_energy", event.target.checked)}
              />
              <span>
                Electron self-energy / spectra
                <InfoTooltip text="Requests electron self-energy and spectral-function style outputs. This can add runtime and output volume." />
              </span>
            </label>
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={goals.transport_mobility}
                onChange={(event) => updateGoal("transport_mobility", event.target.checked)}
              />
              <span>
                Transport mobility
                <InfoTooltip text="Enables EPW scattering and mobility flags for carrier transport analysis. Carrier/scissor details remain in advanced overrides." />
              </span>
            </label>
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={goals.superconductivity}
                onChange={(event) => updateGoal("superconductivity", event.target.checked)}
              />
              <span>
                Superconductivity
                <InfoTooltip text="Enables isotropic Eliashberg superconductivity outputs. Temperature, Coulomb parameter, and convergence details can be overridden in advanced settings." />
              </span>
            </label>
          </div>
        </div>

        <div className="phonon-grid">
          <div className="phonon-field">
            <label>
              Fine K Mesh
              <InfoTooltip text="Fine electron interpolation grid (`nkf1 nkf2 nkf3`). The coarse NSCF/Wannier grid is inferred from the selected source." />
            </label>
            <div className="qgrid-inputs">
              <input type="number" min={1} value={kMesh[0]} onChange={(event) => updateMesh("k", 0, event.target.value)} />
              <span>x</span>
              <input type="number" min={1} value={kMesh[1]} onChange={(event) => updateMesh("k", 1, event.target.value)} />
              <span>x</span>
              <input type="number" min={1} value={kMesh[2]} onChange={(event) => updateMesh("k", 2, event.target.value)} />
            </div>
            {!isValidatingPrerequisites && kMeshCompatibilityMessage && (
              <span className={`param-hint ${kMeshCompatibilityMessage.isError ? "input-error" : ""}`}>
                {kMeshCompatibilityMessage.text}
              </span>
            )}
          </div>

          <div className="phonon-field">
            <label>
              Phonon Q Mesh
              <InfoTooltip text="Coarse phonon grid (`nq1 nq2 nq3`) read from the selected phonon run." />
            </label>
            <div className="qgrid-inputs">
              <input type="number" value={selectedQGrid?.[0] ?? qMesh[0]} readOnly aria-readonly="true" />
              <span>x</span>
              <input type="number" value={selectedQGrid?.[1] ?? qMesh[1]} readOnly aria-readonly="true" />
              <span>x</span>
              <input type="number" value={selectedQGrid?.[2] ?? qMesh[2]} readOnly aria-readonly="true" />
            </div>
            {!selectedQGrid && (
              <span className="param-hint input-error">
                The selected phonon run does not expose q-grid metadata.
              </span>
            )}
            {isValidatingPrerequisites && (
              <span className="param-hint">Checking EPW source compatibility...</span>
            )}
            {!isValidatingPrerequisites && qMeshCompatibilityMessage && (
              <span className={`param-hint ${qMeshCompatibilityMessage.isError ? "input-error" : ""}`}>
                {qMeshCompatibilityMessage.text}
              </span>
            )}
          </div>
        </div>

        {!isHpcMode && (
          <div className="option-section epw-controls-section">
            <h4>Local MPI</h4>
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={mpiEnabled}
                onChange={(event) => setMpiEnabled(event.target.checked)}
                disabled={!mpiAvailable}
              />
              <span>
                Enable MPI
                <InfoTooltip text="Launch `epw.x` through MPI locally. Process count is controlled below." />
              </span>
            </label>
            {!mpiAvailable && (
              <p className="param-hint">
                MPI launcher not detected. Install/configure `mpirun` to enable local MPI execution.
              </p>
            )}
            {mpiEnabled && (
              <div className="phonon-field inline-field">
                <label>
                  MPI Processes
                  <InfoTooltip text="Local MPI rank count. Higher values can reduce runtime but increase resource contention." />
                </label>
                <input
                  type="number"
                  min={1}
                  max={cpuCount}
                  value={mpiProcs}
                  onChange={(event) => {
                    const parsed = Number.parseInt(event.target.value, 10);
                    if (!Number.isFinite(parsed)) return;
                    setMpiProcs(Math.max(1, Math.min(cpuCount, parsed)));
                  }}
                />
              </div>
            )}
          </div>
        )}

        <details className="option-section epw-advanced-details">
          <summary>
            <span className="epw-summary-label">
              <span className="epw-summary-chevron" aria-hidden="true">›</span>
              Advanced
            </span>
            <span onClick={(event) => event.stopPropagation()}>
              <InfoTooltip text="Overrides for EPW experts. Defaults are inferred from selected sources and goals." />
            </span>
          </summary>
          <div className="epw-controls-section">
            <div className="phonon-grid">
              <div className="phonon-field">
                <label>
                  Run Mode
                  <InfoTooltip text="Preset EPW mode label retained for compatibility with saved runs." />
                </label>
                <select value={mode} onChange={(event) => setMode(event.target.value)}>
                  <option value="standard">standard</option>
                  <option value="epb-build">epb-build</option>
                  <option value="epw-read">epw-read</option>
                </select>
              </div>
              <div className="phonon-field">
                <label>prefix</label>
                <input value={prefix} onChange={(event) => setPrefix(event.target.value)} />
              </div>
              <div className="phonon-field">
                <label>outdir</label>
                <input value={outdir} onChange={(event) => setOutdir(event.target.value)} />
              </div>
              <div className="phonon-field">
                <label>dvscf_dir</label>
                <input value={dvscfDir} onChange={(event) => setDvscfDir(event.target.value)} />
              </div>
              <div className="phonon-field">
                <label>wannier_dir</label>
                <input value={wannierDir} onChange={(event) => setWannierDir(event.target.value)} />
              </div>
              <div className="phonon-field">
                <label>
                  Coarse K Mesh
                  <InfoTooltip text="Inferred from the selected Wannier source when present; EPW maps this to `nk1 nk2 nk3`." />
                </label>
                <input
                  value={selectedCoarseKGrid ? `${selectedCoarseKGrid[0]} x ${selectedCoarseKGrid[1]} x ${selectedCoarseKGrid[2]}` : "inferred from source or fine K mesh"}
                  readOnly
                  aria-readonly="true"
                />
              </div>
              <div className="phonon-field">
                <label>
                  Fine Q Mesh
                  <InfoTooltip text="Fine phonon interpolation grid (`nqf1 nqf2 nqf3`). Defaults to the fine K mesh." />
                </label>
                <div className="qgrid-inputs">
                  <input type="number" min={1} value={fineQMesh[0]} onChange={(event) => updateMesh("fine-q", 0, event.target.value)} />
                  <span>x</span>
                  <input type="number" min={1} value={fineQMesh[1]} onChange={(event) => updateMesh("fine-q", 1, event.target.value)} />
                  <span>x</span>
                  <input type="number" min={1} value={fineQMesh[2]} onChange={(event) => updateMesh("fine-q", 2, event.target.value)} />
                </div>
              </div>
              <div className="phonon-field">
                <label>fsthick (eV)</label>
                <input value={fsthickInput} onChange={(event) => setFsthickInput(event.target.value)} placeholder="0.4" />
              </div>
              <div className="phonon-field">
                <label>degaussw (eV)</label>
                <input value={degausswInput} onChange={(event) => setDegausswInput(event.target.value)} placeholder="0.02" />
              </div>
              <div className="phonon-field">
                <label>nbndsub</label>
                <input value={nbndsubInput} onChange={(event) => setNbndsubInput(event.target.value)} placeholder="optional" />
              </div>
            </div>

            <div className="epw-raw-flag-grid">
              <label className="option-checkbox">
                <input type="checkbox" checked={epbwrite} onChange={(event) => setEpbwrite(event.target.checked)} />
                <span>epbwrite</span>
              </label>
              <label className="option-checkbox">
                <input type="checkbox" checked={epbread} onChange={(event) => setEpbread(event.target.checked)} />
                <span>epbread</span>
              </label>
              <label className="option-checkbox">
                <input type="checkbox" checked={epwwrite} onChange={(event) => setEpwwrite(event.target.checked)} />
                <span>epwwrite</span>
              </label>
              <label className="option-checkbox">
                <input type="checkbox" checked={epwread} onChange={(event) => setEpwread(event.target.checked)} />
                <span>epwread</span>
              </label>
              <label className="option-checkbox">
                <input type="checkbox" checked={wannierize} onChange={(event) => setWannierize(event.target.checked)} />
                <span>wannierize</span>
              </label>
            </div>

            {isHpcMode && (
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>Artifact Sync Mode</label>
                  <select
                    value={artifactSyncMode}
                    onChange={(event) => setArtifactSyncMode(event.target.value as "minimal" | "epw-ready" | "full")}
                  >
                    <option value="minimal">minimal</option>
                    <option value="epw-ready">epw-ready</option>
                    <option value="full">full</option>
                  </select>
                </div>
                <div className="phonon-field">
                  <label>MPI Pools</label>
                  <input
                    value={runtimePoolsInput}
                    onChange={(event) => setRuntimePoolsInput(event.target.value)}
                    placeholder="optional"
                  />
                </div>
                <div className="phonon-field">
                  <label>Max Seconds</label>
                  <input
                    value={runtimeMaxSecondsInput}
                    onChange={(event) => setRuntimeMaxSecondsInput(event.target.value)}
                    placeholder="optional"
                  />
                </div>
              </div>
            )}

            <h4>
              Keyword Overrides
              <InfoTooltip text="Free-form `key=value` EPW keywords merged last; use for unsupported or new upstream options." />
            </h4>
            <p className="param-hint">
              One keyword per line. Blank lines and lines starting with `#` or `!` are ignored.
            </p>
            <textarea
              className="epw-overrides-textarea"
              value={advancedOverridesText}
              onChange={(event) => setAdvancedOverridesText(event.target.value)}
              rows={8}
              placeholder={"# one per line\n\nfsthick = 0.35\nnbndsub = 24"}
            />
            <div className="epw-overrides-meta">
              <span>Parsed overrides: {Object.keys(overrideMap).length}</span>
            </div>
            <pre className="epw-overrides-preview">{prettyOverridesPreview(overrideMap)}</pre>
          </div>
        </details>

        {isHpcMode && (
          <HpcRunSettings
            profileId={activeHpcProfile?.id ?? null}
            profileName={activeHpcProfile?.name ?? "Andromeda"}
            taskKind="epw"
            commandLines={hpcCommandLines}
            resources={hpcResources}
            resourceMode={activeHpcProfile?.resource_mode ?? null}
            defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
            defaultGpuResources={activeHpcProfile?.default_gpu_resources ?? null}
            onResourcesChange={setHpcResources}
          />
        )}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>Back</button>
          {!isHpcMode && (
            <button className="secondary-button" onClick={queueCalculation}>
              Queue
            </button>
          )}
          <button
            className="primary-button"
            onClick={() => void runCalculation()}
            disabled={!selectedPhononId || !selectedQGrid || hasBlockingExternalTask}
          >
            {isHpcMode ? "Submit EPW to Andromeda" : "Run EPW"}
          </button>
        </div>
      </div>
    );
  }

  function renderRunStep() {
    return (
      <div className="wizard-step run-step run-step-focused">
        <div className="run-step-headline">
          <h3>{isRunning ? "Running EPW Workflow" : "EPW Output"}</h3>
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
              <button className="secondary-button" onClick={() => setStep("controls")}>
                Back to Controls
              </button>
            </div>
          </>
        )}

        <div className="run-layout run-layout-single">
          <LiveOutputPanel
            title={isRunning ? "Running..." : "Output"}
            output={output}
            placeholder="Starting EPW calculation..."
            totalLineCount={outputLineCount}
            visibleLineCount={visibleOutputLineCount}
          />
        </div>
      </div>
    );
  }

  function renderResultsStep() {
    if (!epwResult) {
      return (
        <div className="wizard-step">
          <h3>No EPW Results</h3>
          <p className="warning-text">No EPW result payload is available.</p>
          <div className="step-actions">
            <button className="secondary-button" onClick={() => setStep("run")}>Back to Run</button>
          </div>
        </div>
      );
    }

    const artifacts = epwResult.artifacts || [];
    const metrics = epwResult.result_summary?.core_metrics || {};
    const parseLimited = Boolean(epwResult.result_summary?.parse_partial);
    const parseNotes = [
      ...(epwResult.result_summary?.notes || []),
      ...(epwResult.result_summary?.warnings || []),
    ];

    return (
      <div className="wizard-step results-step">
        <h3>EPW Results</h3>
        <p className="step-description">
          {epwResult.result_summary?.completed ? "EPW run completed." : "EPW run did not complete."}
          {parseLimited ? " Parsed results are limited." : " Parsed results are complete."}
        </p>

        <div className="details-grid">
          <div className="detail-item">
            <label>Completed</label>
            <span>{epwResult.result_summary?.completed ? "Yes" : "No"}</span>
          </div>
          <div className="detail-item">
            <label>Parsed Results</label>
            <span>{parseLimited ? "Limited" : "Complete"}</span>
          </div>
          <div className="detail-item">
            <label>Elapsed (s)</label>
            <span>{epwResult.result_summary?.elapsed_seconds ?? "N/A"}</span>
          </div>
          <div className="detail-item">
            <label>Artifacts</label>
            <span>{artifacts.length}</span>
          </div>
        </div>

        {Object.keys(metrics).length > 0 && (
          <div className="option-section">
            <h4>Core Metrics</h4>
            <div className="details-grid">
              {Object.entries(metrics).map(([key, value]) => (
                <div key={key} className="detail-item">
                  <label>{key}</label>
                  <span>{value}</span>
                </div>
              ))}
            </div>
          </div>
        )}

        {epwResult.errors?.length > 0 && (
          <div className="option-section">
            <h4>Run Notes</h4>
            <ul className="warning-list">
              {epwResult.errors.map((entry, index) => (
                <li key={`${entry.code}-${index}`}>
                  [{entry.code}] {entry.message} {entry.hint ? `Hint: ${entry.hint}` : ""}
                </li>
              ))}
            </ul>
          </div>
        )}

        {parseLimited && parseNotes.length > 0 && (
          <div className="option-section">
            <h4>Parse Notes</h4>
            <ul className="warning-list">
              {parseNotes.slice(0, 8).map((entry, index) => (
                <li key={`${entry}-${index}`}>{entry}</li>
              ))}
            </ul>
          </div>
        )}

        <div className="option-section">
          <h4>Generated Outputs</h4>
          <div className="calculations-list">
            {artifacts.map((artifact) => (
              <div key={artifact.file_name} className="calculation-item">
                <div className="calculation-header">
                  <div className="calculation-info">
                    <span className="calc-type">{artifact.file_name}</span>
                  </div>
                  <div className="calculation-meta">
                    <span className="calc-size">{formatBytes(Number(artifact.size_bytes) || 0)}</span>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("run")}>
            Back to Run
          </button>
          <button className="secondary-button" onClick={() => onViewEPW(epwResult, output || null)}>
            Open Viewer
          </button>
          <button className="primary-button" onClick={onBack}>
            Back to Dashboard
          </button>
        </div>

        {!isSaved && (
          <p className="warning-text">Result is not saved yet. Re-run or save may have failed.</p>
        )}
      </div>
    );
  }

  return (
    <div className={`epw-wizard wizard-step-${step}`}>
      <AppHeaderPortal className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>EPW Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "controls" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            2. Controls
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "completed" : ""}>
            3. Run
          </span>
          <span className={step === "results" ? "active" : ""}>
            4. Results
          </span>
        </div>
      </AppHeaderPortal>

      {error && <div className="error-banner">{error}</div>}

      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "controls" && renderControlsStep()}
        {step === "run" && renderRunStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
