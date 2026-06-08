import { useEffect, useMemo, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { AppHeaderPortal } from "./AppHeaderPortal";
import {
  CrystalData,
  ExecutionMode,
  HpcProfile,
  SlurmResourceRequest,
} from "../lib/types";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { defaultProgressState, ProgressState } from "../lib/engines/qe/progress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import {
  buildExecutionTarget,
  buildHpcLauncherCommand,
  defaultResourcesForProfile,
  saveExecutionMode,
} from "../lib/hpcConfig";
import {
  buildHpcQeRuntimeSetupLines,
  resolveProfileRemoteQeAuxiliaryExecutable,
} from "../lib/engines/qe/hpc";
import { HpcRunSettings } from "./HpcRunSettings";
import { TransportPlot } from "./TransportPlot";
import { getWannierIssueCounts, getWannierQualityIssues } from "../lib/engines/qe/wannierQuality";
import { formatCalculationSourceLabel, getCalculationName } from "../lib/calculationNames";
import { readProjectWizardSettings, writeProjectWizardSettings } from "../lib/projectWizardSettings";
import type { TransportResult } from "../lib/transport";
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
    wannier_data?: any;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

interface TransportWizardProps {
  onBack: () => void;
  onExecutionModeChange?: (mode: ExecutionMode) => Promise<void> | void;
  qePath: string;
  executionMode?: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  wannierCalculations: CalculationRun[];
  reconnectTaskId?: string;
  onViewTransport: (transportData: TransportResult) => void;
}

type WizardStep = "source" | "parameters" | "run" | "results";
type TransportSourceSortMode = "recent" | "best";

const TRANSPORT_WORK_DIR = "/tmp/qcortado_transport";
const TRANSPORT_WIZARD_SETTINGS_ID = "transport";

interface StoredTransportWizardSettings {
  boltzKMesh: [number, number, number];
  muOffsetMinInput: string;
  muOffsetMaxInput: string;
  muOffsetStepInput: string;
  tempMinInput: string;
  tempMaxInput: string;
  tempStepInput: string;
  relaxationTimeInput: string;
  is2d: boolean;
  boltz2dDir: string;
}

function formatNumber(value: number | null | undefined, digits = 4): string {
  if (!Number.isFinite(Number(value))) {
    return "N/A";
  }
  return Number(value).toFixed(digits);
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

function parsePositiveNumber(value: string, label: string): number {
  const parsed = Number(value.trim());
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive number.`);
  }
  return parsed;
}

function parseInteger(value: string, label: string): number {
  const parsed = Number.parseInt(value.trim(), 10);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive integer.`);
  }
  return parsed;
}

function deriveRemotePostw90Path(profile: HpcProfile | null | undefined): string {
  const remoteWannier90 = (profile?.remote_wannier90_path || "").trim();
  if (remoteWannier90.includes("/") || remoteWannier90.startsWith("~")) {
    const segments = remoteWannier90.split("/");
    segments[segments.length - 1] = "postw90.x";
    return segments.join("/");
  }
  return "postw90.x";
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

function sourceSeedname(calc: CalculationRun): string {
  return String(
    calc.parameters?.seedname
      || calc.result?.wannier_data?.seedname
      || "qcortado_wannier",
  );
}

function sourceWannierSpread(calc: CalculationRun): number | null {
  const spread = Number(calc.parameters?.total_spread ?? calc.result?.wannier_data?.total_spread);
  return Number.isFinite(spread) ? spread : null;
}

function sourceWannierNumWann(calc: CalculationRun): number | null {
  const value = Number(calc.parameters?.n_wann ?? calc.result?.wannier_data?.num_wann);
  return Number.isFinite(value) && value > 0 ? value : null;
}

function sourceWannierNumBands(calc: CalculationRun): number | null {
  const value = Number(calc.parameters?.n_bands ?? calc.result?.wannier_data?.num_bands);
  return Number.isFinite(value) && value > 0 ? value : null;
}

function sourceWannierConverged(calc: CalculationRun): boolean {
  const issueCounts = getWannierIssueCounts(
    calc.result?.wannier_data ?? null,
    calc.result?.raw_output ?? null,
    calc.result?.fermi_energy ?? null,
  );
  if (issueCounts.errors > 0) {
    return false;
  }
  return Boolean(calc.result?.converged || calc.result?.wannier_data?.convergence?.converged);
}

function defaultBoltzKMesh(calc: CalculationRun | null | undefined): [number, number, number] {
  const parameterGrid = Array.isArray(calc?.parameters?.k_grid)
    ? calc?.parameters?.k_grid
    : Array.isArray(calc?.result?.wannier_data?.k_grid)
      ? calc?.result?.wannier_data?.k_grid
      : [];
  const values = parameterGrid.map((value: unknown) => Number(value));
  if (values.length === 3 && values.every((value: number) => Number.isInteger(value) && value > 0)) {
    return [values[0], values[1], values[2]];
  }
  return [24, 24, 24];
}

function sourceWannierSortTimestamp(calc: CalculationRun): number {
  const iso = calc.completed_at || calc.started_at;
  const timestamp = new Date(iso).getTime();
  return Number.isFinite(timestamp) ? timestamp : 0;
}

function sortWannierCalculations(calculations: CalculationRun[], mode: TransportSourceSortMode): CalculationRun[] {
  return [...calculations].sort((left, right) => {
    if (mode === "best") {
      const leftConverged = sourceWannierConverged(left) ? 1 : 0;
      const rightConverged = sourceWannierConverged(right) ? 1 : 0;
      if (leftConverged !== rightConverged) {
        return rightConverged - leftConverged;
      }

      const leftSpread = sourceWannierSpread(left);
      const rightSpread = sourceWannierSpread(right);
      if (leftSpread != null && rightSpread != null && leftSpread !== rightSpread) {
        return leftSpread - rightSpread;
      }
      if (leftSpread != null && rightSpread == null) return -1;
      if (leftSpread == null && rightSpread != null) return 1;
    }

    return sourceWannierSortTimestamp(right) - sourceWannierSortTimestamp(left);
  });
}

function getWannierSourceTags(calc: CalculationRun): Array<{ label: string; type: "info" | "feature" | "special" }> {
  const tags: Array<{ label: string; type: "info" | "feature" | "special" }> = [];
  const issueCounts = getWannierIssueCounts(
    calc.result?.wannier_data ?? null,
    calc.result?.raw_output ?? null,
    calc.result?.fermi_energy ?? null,
  );
  const pushTag = (label: string, type: "info" | "feature" | "special") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  const spread = sourceWannierSpread(calc);
  if (spread != null) {
    pushTag(`Ω ${spread.toFixed(3)}`, "info");
  }
  const numWann = sourceWannierNumWann(calc);
  if (numWann != null) {
    pushTag(`${numWann} WF`, "info");
  }
  if (sourceWannierConverged(calc)) {
    pushTag("Converged", "special");
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

interface TransportTaskPlan {
  taskLabel: string;
  taskParams: Record<string, any>;
  saveParameters: Record<string, any>;
}

export function TransportWizard({
  onBack,
  onExecutionModeChange,
  qePath,
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId,
  crystalData,
  wannierCalculations,
  reconnectTaskId,
  onViewTransport,
}: TransportWizardProps) {
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const storedWizardSettings = useMemo(
    () => readProjectWizardSettings<StoredTransportWizardSettings>(projectId, TRANSPORT_WIZARD_SETTINGS_ID),
    [projectId],
  );
  const shouldPreserveStoredSettingsRef = useRef(Boolean(storedWizardSettings));
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  const hasBlockingExternalTask = !isHpcMode && hasExternalRunningTask;

  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);
  const [selectedWannierCalculation, setSelectedWannierCalculation] = useState<CalculationRun | null>(null);
  const [sourceSortMode, setSourceSortMode] = useState<TransportSourceSortMode>("recent");
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const [outputLineCount, setOutputLineCount] = useState(0);
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Transport",
  });
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [transportData, setTransportData] = useState<TransportResult | null>(null);
  const [isSaved, setIsSaved] = useState(false);

  const [boltzKMesh, setBoltzKMesh] = useState<[number, number, number]>(() => storedWizardSettings?.boltzKMesh ?? [24, 24, 24]);
  const [muOffsetMinInput, setMuOffsetMinInput] = useState(() => storedWizardSettings?.muOffsetMinInput ?? "-0.5");
  const [muOffsetMaxInput, setMuOffsetMaxInput] = useState(() => storedWizardSettings?.muOffsetMaxInput ?? "0.5");
  const [muOffsetStepInput, setMuOffsetStepInput] = useState(() => storedWizardSettings?.muOffsetStepInput ?? "0.05");
  const [tempMinInput, setTempMinInput] = useState(() => storedWizardSettings?.tempMinInput ?? "100");
  const [tempMaxInput, setTempMaxInput] = useState(() => storedWizardSettings?.tempMaxInput ?? "800");
  const [tempStepInput, setTempStepInput] = useState(() => storedWizardSettings?.tempStepInput ?? "50");
  const [relaxationTimeInput, setRelaxationTimeInput] = useState(() => storedWizardSettings?.relaxationTimeInput ?? "10");
  const [is2d, setIs2d] = useState(() => storedWizardSettings?.is2d ?? false);
  const [boltz2dDir, setBoltz2dDir] = useState(() => storedWizardSettings?.boltz2dDir ?? "z");

  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(
    defaultResourcesForProfile(activeHpcProfile),
  );

  const visibleOutputLineCount = useMemo(() => countVisibleOutputLines(output), [output]);
  useViewportScrollLock(step === "run");

  const sourceFermiEnergy = selectedWannierCalculation?.result?.fermi_energy ?? null;
  const seedname = selectedWannierCalculation ? sourceSeedname(selectedWannierCalculation) : "qcortado_wannier";
  const sourceSpread = selectedWannierCalculation ? sourceWannierSpread(selectedWannierCalculation) : null;
  const sourceNumWann = selectedWannierCalculation ? sourceWannierNumWann(selectedWannierCalculation) : null;
  const sourceNumBands = selectedWannierCalculation ? sourceWannierNumBands(selectedWannierCalculation) : null;
  const selectedWannierIssues = useMemo(
    () => selectedWannierCalculation
      ? getWannierQualityIssues(
        selectedWannierCalculation.result?.wannier_data ?? null,
        selectedWannierCalculation.result?.raw_output ?? null,
        selectedWannierCalculation.result?.fermi_energy ?? null,
      )
      : [],
    [selectedWannierCalculation],
  );
  const availableWannierCalculations = useMemo(
    () => wannierCalculations.filter((calc) => calc.calc_type === "wannier"),
    [wannierCalculations],
  );
  const sortedWannierCalculations = useMemo(
    () => sortWannierCalculations(availableWannierCalculations, sourceSortMode),
    [availableWannierCalculations, sourceSortMode],
  );

  const hpcCommandLines = useMemo(
    () => {
      const remotePostw90 = resolveProfileRemoteQeAuxiliaryExecutable(
        activeHpcProfile,
        deriveRemotePostw90Path(activeHpcProfile),
        "postw90.x",
      );
      return [
        "cd \"$SLURM_SUBMIT_DIR\"",
        ...buildHpcQeRuntimeSetupLines(activeHpcProfile, hpcResources.resource_type),
        `${buildHpcLauncherCommand(activeHpcProfile, hpcResources.resource_type)} ${remotePostw90} ${seedname} > ${seedname}.wpout 2> ${seedname}.werr`,
      ];
    },
    [activeHpcProfile, hpcResources.resource_type, seedname],
  );

  function selectSourceWannier(calc: CalculationRun) {
    setSelectedWannierCalculation(calc);
    if (!shouldPreserveStoredSettingsRef.current) {
      setBoltzKMesh(defaultBoltzKMesh(calc));
    }
    setError(null);
  }

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.resource_mode, isHpcMode]);

  useEffect(() => {
    writeProjectWizardSettings(projectId, TRANSPORT_WIZARD_SETTINGS_ID, {
      boltzKMesh,
      muOffsetMinInput,
      muOffsetMaxInput,
      muOffsetStepInput,
      tempMinInput,
      tempMaxInput,
      tempStepInput,
      relaxationTimeInput,
      is2d,
      boltz2dDir,
    });
  }, [
    projectId,
    boltzKMesh,
    muOffsetMinInput,
    muOffsetMaxInput,
    muOffsetStepInput,
    tempMinInput,
    tempMaxInput,
    tempStepInput,
    relaxationTimeInput,
    is2d,
    boltz2dDir,
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
        console.error("Failed to initialize transport wizard:", initError);
      }
    }

    void init();
  }, [qePath]);

  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) {
      void taskContext.reconnectToTask(activeTaskId);
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
      setTransportData(task.result as TransportResult);
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
    if (task.outputLineCount > 0) {
      setOutput(task.outputText);
      setOutputLineCount(task.outputLineCount);
    }
    setProgress(task.progress);
    setIsRunning(task.status === "running");
  }, [
    activeTaskId,
    taskContext,
    taskContext.getTask(activeTaskId ?? "")?.outputLineCount,
    taskContext.getTask(activeTaskId ?? "")?.status,
  ]);

  async function handleSwitchToHpcMode() {
    setError(null);
    try {
      if (onExecutionModeChange) {
        await onExecutionModeChange("hpc");
      } else {
        await saveExecutionMode("hpc");
      }
    } catch (switchError) {
      setError(`Failed to switch execution mode: ${switchError}`);
    }
  }

  async function buildTaskPlan(): Promise<TransportTaskPlan> {
    if (!selectedWannierCalculation?.id) {
      throw new Error("No saved Wannier calculation selected.");
    }
    if (!Number.isFinite(Number(sourceFermiEnergy))) {
      throw new Error("The selected Wannier calculation has no saved Fermi energy.");
    }

    const muOffsetMin = Number(muOffsetMinInput.trim());
    const muOffsetMax = Number(muOffsetMaxInput.trim());
    const muOffsetStep = parsePositiveNumber(muOffsetStepInput, "Chemical-potential step");
    const tempMin = parsePositiveNumber(tempMinInput, "Minimum temperature");
    const tempMax = parsePositiveNumber(tempMaxInput, "Maximum temperature");
    const tempStep = parsePositiveNumber(tempStepInput, "Temperature step");
    const relaxationTimeFs = parsePositiveNumber(relaxationTimeInput, "Relaxation time");

    if (!Number.isFinite(muOffsetMin) || !Number.isFinite(muOffsetMax)) {
      throw new Error("Chemical-potential limits must be finite numbers.");
    }
    if (muOffsetMax < muOffsetMin) {
      throw new Error("Maximum Δμ must be greater than or equal to minimum Δμ.");
    }
    if (tempMax < tempMin) {
      throw new Error("Maximum temperature must be greater than or equal to minimum temperature.");
    }

    const taskLabel = `BoltzWann Transport - ${crystalData?.formula_sum || seedname}`;
    const taskParams = {
      config: {
        project_id: projectId,
        source_wannier_calc_id: selectedWannierCalculation.id,
        boltz_kmesh: boltzKMesh,
        mu_offset_min: muOffsetMin,
        mu_offset_max: muOffsetMax,
        mu_offset_step: muOffsetStep,
        temp_min: tempMin,
        temp_max: tempMax,
        temp_step: tempStep,
        relaxation_time_fs: relaxationTimeFs,
        is_2d: is2d,
        boltz_2d_dir: is2d ? boltz2dDir : null,
      },
      workingDir: TRANSPORT_WORK_DIR,
      mpiConfig: !isHpcMode && mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      executionTarget: buildExecutionTarget(
        executionMode,
        activeHpcProfile?.id ?? null,
        isHpcMode ? hpcResources : null,
        false,
      ),
    };

    return {
      taskLabel,
      taskParams,
      saveParameters: {
        engine: "boltzwann",
        seedname,
        source_wannier_calc_id: selectedWannierCalculation.id,
        boltz_kmesh: boltzKMesh,
        mu_offset_min: muOffsetMin,
        mu_offset_max: muOffsetMax,
        mu_offset_step: muOffsetStep,
        temp_min: tempMin,
        temp_max: tempMax,
        temp_step: tempStep,
        relaxation_time_fs: relaxationTimeFs,
        reference_fermi_energy_ev: sourceFermiEnergy,
        is_2d: is2d,
        boltz_2d_dir: is2d ? boltz2dDir : null,
        mu_points: null,
        temperature_points: null,
      },
    };
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
    setTransportData(null);
    setIsSaved(false);
    setProgress(defaultProgressState("BoltzWann Transport"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

    try {
      const plan = await buildTaskPlan();
      const hpcSaveSpec = isHpcMode
        ? {
          projectId,
          cifId,
          workingDir: TRANSPORT_WORK_DIR,
          calcType: "transport" as const,
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        }
        : null;
      const taskId = await taskContext.startTask("transport", plan.taskParams, plan.taskLabel, hpcSaveSpec);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Transport calculation failed");
      }

      const result = finalTask.result as TransportResult;
      const outputContent = finalTask.output.join("\n");
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

      setTransportData(result);
      setStep("results");
      setProgress({
        status: "complete",
        percent: 100,
        phase: "Complete",
      });

      try {
        await invoke("save_calculation", {
          projectId,
          cifId,
          calcData: {
            calc_type: "transport",
            parameters: {
              ...plan.saveParameters,
              mu_points: Array.isArray(result?.mu_values_ev) ? result.mu_values_ev.length : null,
              temperature_points: Array.isArray(result?.temperature_values_k) ? result.temperature_values_k.length : null,
              ...hpcSaveParams,
            },
            result: {
              converged: true,
              total_energy: null,
              fermi_energy: result?.reference_fermi_energy_ev ?? selectedWannierCalculation?.result?.fermi_energy ?? sourceFermiEnergy,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: outputContent,
              transport_data: result,
            },
            started_at: startTime,
            completed_at: completedAt,
            input_content: "",
            output_content: outputContent,
            tags: [],
          },
          workingDir: finalTask.hpc.local_sync_dir ?? TRANSPORT_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save transport calculation:", saveError);
        setError(`Failed to auto-save transport calculation: ${saveError}`);
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

  async function queueCalculation() {
    if (isHpcMode) {
      setError("Queueing is unavailable in HPC mode. Submit directly to Andromeda.");
      return;
    }

    try {
      const plan = await buildTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "transport",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId,
          workingDir: TRANSPORT_WORK_DIR,
          calcType: "transport",
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        },
      );
    } catch (queueError) {
      setError(String(queueError));
    }
  }

  const renderSourceStep = () => {
    if (availableWannierCalculations.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No Wannier Calculations Available</h3>
          <p className="warning-text">
            BoltzWann transport requires a saved Wannier calculation. Run and save Wannier90 first, then return here.
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
          <h3>Select Source Wannier Calculation</h3>
          <div className="source-sort-control">
            <label htmlFor="transport-wannier-sort">Sort Wannier</label>
            <select
              id="transport-wannier-sort"
              value={sourceSortMode}
              onChange={(event) => setSourceSortMode(event.target.value as TransportSourceSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose the saved Wannier calculation that BoltzWann should reuse as its transport source.
        </p>

        <div className="scf-list">
          {sortedWannierCalculations.map((calc) => {
            const spread = sourceWannierSpread(calc);
            const fermiEnergy = calc.result?.fermi_energy;
            const numWann = sourceWannierNumWann(calc);
            const numBands = sourceWannierNumBands(calc);
            const calcName = getCalculationName(calc);
            const sourceIssues = getWannierQualityIssues(
              calc.result?.wannier_data ?? null,
              calc.result?.raw_output ?? null,
              calc.result?.fermi_energy ?? null,
            );

            return (
              <div
                key={calc.id}
                className={`scf-option ${selectedWannierCalculation?.id === calc.id ? "selected" : ""}`}
                onClick={() => selectSourceWannier(calc)}
              >
                <div className="scf-option-header">
                  <input
                    type="radio"
                    checked={selectedWannierCalculation?.id === calc.id}
                    onChange={() => selectSourceWannier(calc)}
                  />
                  {calcName && (
                    <span className="scf-name" title={formatCalculationSourceLabel(calc, "Wannier")}>
                      {calcName}
                    </span>
                  )}
                  <span className="scf-date">
                    {new Date(calc.completed_at || calc.started_at).toLocaleDateString()}
                  </span>
                </div>
                <div className="scf-details">
                  <span>seed = {sourceSeedname(calc)}</span>
                  {numWann != null && (
                    <span>num_wann / num_bands = {numWann} / {numBands ?? "N/A"}</span>
                  )}
                  {fermiEnergy != null && (
                    <span>EF = {fermiEnergy.toFixed(3)} eV</span>
                  )}
                  {spread != null && (
                    <span>Ω = {spread.toFixed(6)} A^2</span>
                  )}
                </div>
                <div className="calc-tags">
                  {getWannierSourceTags(calc).map((tag, index) => (
                    <span
                      key={`${tag.label}-${index}`}
                      className={`calc-tag calc-tag-${tag.type}${tag.label.trim().toUpperCase() === "HPC" ? " calc-tag-hpc" : ""}`}
                    >
                      {tag.label}
                    </span>
                  ))}
                </div>
                {sourceIssues.length > 0 && (
                  <div className="param-hint" style={{ marginTop: "0.5rem" }}>
                    {sourceIssues[0].message}
                  </div>
                )}
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
            disabled={!selectedWannierCalculation}
            onClick={() => setStep("parameters")}
          >
            Next: Parameters
          </button>
        </div>
      </div>
    );
  };

  const renderParametersStep = () => (
    <div className="wizard-step parameters-step">
      <h3>BoltzWann Transport Parameters</h3>
      <p className="step-description">
        QCortado will copy the saved Wannier artifacts, write a transport-specific `.win`, and run
        `postw90.x` with `boltzwann = true`.
      </p>
      <details className="wannier-guidance">
        <summary>How to choose a BoltzWann-ready Wannier source</summary>
        <div className="guidance-body">
          <p>Prefer Wannier runs whose interpolated bands reach the source SCF Fermi level. If the Wannier manifold misses E_F, BoltzWann transport will usually come back nearly zero.</p>
          <p>Choose projections that span the orbitals crossing E_F, then use a larger `num_bands` plus `dis_win_*` and `dis_froz_*` when the transport manifold is entangled or metallic.</p>
          <p>QCortado flags saved Wannier runs that failed disentanglement or appear to miss E_F. Treat those warnings as blockers for transport-quality results.</p>
        </div>
      </details>

      {!selectedWannierCalculation ? (
        <div className="warning-banner">
          Select a source Wannier calculation before configuring BoltzWann transport.
        </div>
      ) : (
      <div className="calculation-summary">
        <h4>Source Wannier Run</h4>
        <div className="summary-row">
          <span>Calculation ID</span>
          <span>{selectedWannierCalculation.id.slice(0, 8)}</span>
        </div>
        <div className="summary-row">
          <span>seedname</span>
          <span>{seedname}</span>
        </div>
        <div className="summary-row">
          <span>Reference Fermi Energy</span>
          <span>{formatNumber(sourceFermiEnergy, 4)} eV</span>
        </div>
        <div className="summary-row">
          <span>num_wann / num_bands</span>
          <span>{sourceNumWann ?? "N/A"} / {sourceNumBands ?? "N/A"}</span>
        </div>
        <div className="summary-row">
          <span>Total Spread</span>
          <span>{sourceSpread != null ? `${sourceSpread.toFixed(6)} A^2` : "N/A"}</span>
        </div>
      </div>
      )}

      {selectedWannierIssues.length > 0 && (
        <div className="warning-banner">
          {selectedWannierIssues.map((issue) => issue.message).join(" ")}
        </div>
      )}

      <div className="option-section config-section">
        <h4>BoltzWann Scan</h4>
        <div className="phonon-grid">
          <div className="phonon-field">
            <label>Boltz k-mesh</label>
            <div className="qgrid-inputs">
              {[0, 1, 2].map((index) => (
                <input
                  key={index}
                  type="number"
                  min={1}
                  value={boltzKMesh[index]}
                  onChange={(event) => {
                    const next = [...boltzKMesh] as [number, number, number];
                    next[index] = Math.max(1, parseInteger(event.target.value || "1", "Boltz k-mesh"));
                    setBoltzKMesh(next);
                  }}
                />
              ))}
            </div>
          </div>

          <div className="phonon-field">
            <label>Δμ min / max / step (eV)</label>
            <div className="qgrid-inputs">
              <input type="text" value={muOffsetMinInput} onChange={(event) => setMuOffsetMinInput(event.target.value)} spellCheck={false} />
              <input type="text" value={muOffsetMaxInput} onChange={(event) => setMuOffsetMaxInput(event.target.value)} spellCheck={false} />
              <input type="text" value={muOffsetStepInput} onChange={(event) => setMuOffsetStepInput(event.target.value)} spellCheck={false} />
            </div>
          </div>

          <div className="phonon-field">
            <label>T min / max / step (K)</label>
            <div className="qgrid-inputs">
              <input type="text" value={tempMinInput} onChange={(event) => setTempMinInput(event.target.value)} spellCheck={false} />
              <input type="text" value={tempMaxInput} onChange={(event) => setTempMaxInput(event.target.value)} spellCheck={false} />
              <input type="text" value={tempStepInput} onChange={(event) => setTempStepInput(event.target.value)} spellCheck={false} />
            </div>
          </div>

          <div className="phonon-field">
            <label>Relaxation time τ (fs)</label>
            <input
              type="text"
              value={relaxationTimeInput}
              onChange={(event) => setRelaxationTimeInput(event.target.value)}
              spellCheck={false}
            />
            <span className="param-hint">This is a user assumption, not a quantity inferred from Wannier data.</span>
          </div>

          <div className="phonon-field">
            <label className="option-checkbox">
              <input
                type="checkbox"
                checked={is2d}
                onChange={(event) => setIs2d(event.target.checked)}
              />
              <span>Treat as 2D transport</span>
            </label>
          </div>

          {is2d && (
            <div className="phonon-field">
              <label>Boltz 2D Direction</label>
              <select value={boltz2dDir} onChange={(event) => setBoltz2dDir(event.target.value)}>
                <option value="x">x</option>
                <option value="y">y</option>
                <option value="z">z</option>
              </select>
            </div>
          )}
        </div>
      </div>

      {isHpcMode ? (
        <HpcRunSettings
          profileId={activeHpcProfile?.id ?? null}
          profileName={activeHpcProfile?.name ?? "Andromeda"}
          taskKind="transport"
          commandLines={hpcCommandLines}
          resources={hpcResources}
          resourceMode={activeHpcProfile?.resource_mode ?? "both"}
          defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
          defaultGpuResources={activeHpcProfile?.default_gpu_resources ?? null}
          onResourcesChange={setHpcResources}
          disabled={isRunning}
        />
      ) : (
        <div className="option-section mpi-section config-section">
          <h4>Parallelization</h4>
          {mpiAvailable ? (
            <div className="mpi-toggle">
              <label>
                <input
                  type="checkbox"
                  checked={mpiEnabled}
                  onChange={(event) => setMpiEnabled(event.target.checked)}
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
                      onChange={(event) => setMpiProcs(Math.max(1, parseInt(event.target.value, 10) || 1))}
                    />
                  </label>
                </div>
              )}
            </div>
          ) : (
            <p className="mpi-unavailable">MPI not available. Running in serial mode.</p>
          )}
        </div>
      )}

      <div className="warning-banner">
        EPW transport is intentionally out of scope here. This workflow uses saved Wannier artifacts and upstream
        BoltzWann semantics, including BoltzWann `K`.
      </div>

      {error && <div className="error-message">{error}</div>}

      <div className="step-actions">
        <button className="secondary-button" onClick={() => setStep("source")}>
          Back to Source
        </button>
        {!isHpcMode && (
          <button className="secondary-button" onClick={() => void queueCalculation()}>
            Queue Task
          </button>
        )}
        {!isHpcMode && (
          <button className="secondary-button" onClick={() => void handleSwitchToHpcMode()}>
            Switch to HPC
          </button>
        )}
        <button
          className="primary-button"
          onClick={() => void runCalculation()}
          disabled={!selectedWannierCalculation || hasBlockingExternalTask}
        >
          {isHpcMode ? "Submit BoltzWann Transport" : "Run BoltzWann Transport"}
        </button>
      </div>
    </div>
  );

  const renderRunStep = () => (
    <div className="wizard-step run-step run-step-focused">
      <div className="run-step-headline">
        <h3>{isRunning ? "Running BoltzWann Transport" : "BoltzWann Transport Output"}</h3>
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
            <button className="secondary-button" onClick={() => setStep(selectedWannierCalculation ? "parameters" : "source")}>
              {selectedWannierCalculation ? "Back to Parameters" : "Back to Source"}
            </button>
          </div>
        </>
      )}

      <div className="run-layout run-layout-single">
        <LiveOutputPanel
          title={isRunning ? "Running..." : "Output"}
          output={output}
          placeholder="Starting transport calculation..."
          totalLineCount={outputLineCount}
          visibleLineCount={visibleOutputLineCount}
        />
      </div>
    </div>
  );

  const renderResultsStep = () => {
    if (!transportData) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep(selectedWannierCalculation ? "parameters" : "source")}>
            Try Again
          </button>
        </div>
      );
    }

    const artifactBytes = Array.isArray(transportData?.artifact_manifest)
      ? transportData.artifact_manifest.reduce((sum: number, item: any) => sum + (Number(item?.size_bytes) || 0), 0)
      : 0;

    return (
      <div className="wizard-step results-step">
        <h3>BoltzWann Transport Results</h3>
        <p className="step-description">
          BoltzWann finished successfully. QCortado preserved the upstream transport tables and provenance.
        </p>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Engine:</span>
              <span className="value">{transportData.engine || "boltzwann"}</span>
            </div>
            <div className="summary-item">
              <span className="label">μ points:</span>
              <span className="value">{transportData.mu_values_ev?.length ?? "N/A"}</span>
            </div>
            <div className="summary-item">
              <span className="label">T points:</span>
              <span className="value">{transportData.temperature_values_k?.length ?? "N/A"}</span>
            </div>
            <div className="summary-item">
              <span className="label">Reference EF:</span>
              <span className="value">{formatNumber(transportData.reference_fermi_energy_ev, 4)} eV</span>
            </div>
            <div className="summary-item">
              <span className="label">τ:</span>
              <span className="value">{formatNumber(transportData.relaxation_time_fs, 2)} fs</span>
            </div>
            <div className="summary-item">
              <span className="label">Artifacts:</span>
              <span className="value">{formatBytes(artifactBytes)}</span>
            </div>
          </div>
        </div>

        <TransportPlot data={transportData} />

        {isSaved && (
          <div className="save-status save-status-inline">
            <span className="saved">Saved to project</span>
          </div>
        )}

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Back to Dashboard
          </button>
          <button className="secondary-button" onClick={() => onViewTransport(transportData)}>
            Open Viewer
          </button>
          <button
            className="primary-button"
            onClick={() => {
              setActiveTaskId(null);
              setTransportData(null);
              setError(null);
              setIsSaved(false);
              setStep("source");
            }}
          >
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className={`transport-wizard wizard-step-${step}`}>
      <AppHeaderPortal className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>BoltzWann Transport Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : ["parameters", "run", "results"].includes(step) ? "completed" : ""}>
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
      </AppHeaderPortal>

      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "parameters" && renderParametersStep()}
        {step === "run" && renderRunStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
