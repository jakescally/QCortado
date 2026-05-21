import { createContext, useContext, useEffect, useMemo, useState, useCallback, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { ProgressState, progressReducer, defaultProgressState } from "./qeProgress";
import { extractOptimizedStructure, isSavedStructureData, summarizeCell } from "./optimizedStructure";
import { buildVisibleOutputWindow } from "./liveOutput";
import { HpcTaskMeta } from "./types";

export type TaskStatus = "running" | "completed" | "failed" | "cancelled";
export type TaskType = "scf" | "bands" | "dos" | "fermi_surface" | "hubbard_lrt" | "phonon" | "epw" | "wannier" | "transport";
export type QueueItemStatus = "queued" | "running" | "saving" | "completed" | "failed" | "cancelled";

export interface TaskState {
  taskId: string;
  taskType: TaskType;
  label: string;
  startedAt: string;
  status: TaskStatus;
  progress: ProgressState;
  output: string[];
  outputText: string;
  outputLineCount: number;
  result: any | null;
  error: string | null;
  hpc: HpcTaskMeta;
}

interface TaskSummary {
  task_id: string;
  task_type: string;
  label: string;
  started_at: string;
  status: TaskStatus;
  backend?: string | null;
  hpc_resource_type?: "cpu" | "gpu" | null;
  remote_job_id?: string | null;
  scheduler_state?: string | null;
  remote_node?: string | null;
  remote_workdir?: string | null;
  remote_project_path?: string | null;
  remote_storage_bytes?: number | null;
  hpc_profile_id?: string | null;
  local_sync_dir?: string | null;
  recovery_save?: HpcRecoverySaveSpec | null;
  headless_attached?: boolean;
}

interface TaskInfo {
  task_id: string;
  task_type: string;
  label: string;
  started_at: string;
  status: TaskStatus;
  result: any | null;
  error: string | null;
  backend?: string | null;
  hpc_resource_type?: "cpu" | "gpu" | null;
  remote_job_id?: string | null;
  scheduler_state?: string | null;
  remote_node?: string | null;
  remote_workdir?: string | null;
  remote_project_path?: string | null;
  remote_storage_bytes?: number | null;
  hpc_profile_id?: string | null;
  local_sync_dir?: string | null;
  recovery_save?: HpcRecoverySaveSpec | null;
  headless_attached?: boolean;
}

interface QueueSaveSpec {
  projectId: string;
  cifId: string;
  workingDir?: string | null;
  calcType: "scf" | "bands" | "dos" | "fermi_surface" | "hubbard_lrt" | "phonon" | "epw" | "wannier" | "transport" | "optimization";
  parameters: Record<string, any>;
  tags?: string[];
  inputContent?: string;
}

interface HpcRecoverySaveSpec {
  project_id: string;
  cif_id: string;
  calc_type: string;
  parameters: Record<string, any>;
  tags?: string[];
  input_content?: string;
}

export interface QueueItem {
  id: string;
  taskType: TaskType;
  label: string;
  params: Record<string, any>;
  status: QueueItemStatus;
  createdAt: string;
  startedAt: string | null;
  finishedAt: string | null;
  taskId: string | null;
  error: string | null;
  saveSpec: QueueSaveSpec | null;
}

export interface QueueSummary {
  total: number;
  activeIndex: number | null;
}

interface TaskContextValue {
  tasks: Map<string, TaskState>;
  activeTasks: TaskState[];
  startTask: (
    type: TaskType,
    params: Record<string, any>,
    label: string,
    saveSpec?: QueueSaveSpec | null,
  ) => Promise<string>;
  enqueueTask: (
    type: TaskType,
    params: Record<string, any>,
    label: string,
    saveSpec?: QueueSaveSpec | null,
  ) => string;
  queueItems: QueueItem[];
  queueSummary: QueueSummary;
  cancelQueueItem: (queueItemId: string) => Promise<void>;
  removeQueueItem: (queueItemId: string) => void;
  moveQueueItem: (queueItemId: string, direction: "up" | "down") => void;
  clearFinishedQueueItems: () => void;
  cancelTask: (taskId: string) => Promise<void>;
  dismissTask: (taskId: string) => Promise<void>;
  getTask: (taskId: string) => TaskState | undefined;
  waitForTaskCompletion: (taskId: string) => Promise<TaskState>;
  waitForQueueItemCompletion: (taskId: string) => Promise<QueueItem | null>;
  reconnectToTask: (taskId: string) => Promise<void>;
}

const TaskContext = createContext<TaskContextValue | null>(null);

const COMMAND_MAP: Record<TaskType, string> = {
  scf: "start_scf_calculation",
  bands: "start_bands_calculation",
  dos: "start_dos_calculation",
  fermi_surface: "start_fermi_surface_calculation",
  hubbard_lrt: "start_hubbard_lrt_calculation",
  phonon: "start_phonon_calculation",
  epw: "start_epw_calculation",
  wannier: "start_wannier_calculation",
  transport: "start_transport_calculation",
};

function sleep(ms: number): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function generateQueueItemId(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") {
    return crypto.randomUUID();
  }
  return `queue_${Date.now()}_${Math.floor(Math.random() * 1_000_000)}`;
}

function isBusyError(error: unknown): boolean {
  const text = String(error).toLowerCase();
  return text.includes("already running");
}

function isHpcStartParams(params: Record<string, any>): boolean {
  const mode = String(params?.executionTarget?.mode || "").toLowerCase();
  return mode === "hpc";
}

function normalizeTaskType(taskType: string): TaskType {
  if (taskType === "scf" || taskType === "bands" || taskType === "dos" || taskType === "fermi_surface" || taskType === "hubbard_lrt" || taskType === "phonon" || taskType === "epw" || taskType === "wannier" || taskType === "transport") {
    return taskType;
  }
  return "scf";
}

function taskInfoToHpcMeta(info: Partial<TaskInfo> | Partial<TaskSummary>): HpcTaskMeta {
  return {
    backend: info.backend ?? null,
    hpc_resource_type: info.hpc_resource_type === "gpu" ? "gpu" : info.hpc_resource_type === "cpu" ? "cpu" : null,
    remote_job_id: info.remote_job_id ?? null,
    scheduler_state: info.scheduler_state ?? null,
    remote_node: info.remote_node ?? null,
    remote_workdir: info.remote_workdir ?? null,
    remote_project_path: info.remote_project_path ?? null,
    remote_storage_bytes: info.remote_storage_bytes ?? null,
    hpc_profile_id: info.hpc_profile_id ?? null,
    local_sync_dir: info.local_sync_dir ?? null,
    recovery_save: info.recovery_save ?? null,
    headless_attached: Boolean(info.headless_attached),
  };
}

function withHpcRecoverySave(
  params: Record<string, any>,
  saveSpec?: QueueSaveSpec | null,
): Record<string, any> {
  if (!saveSpec || !isHpcStartParams(params)) {
    return params;
  }
  return {
    ...params,
    executionTarget: {
      ...(params.executionTarget || {}),
      hpc: {
        ...(params.executionTarget?.hpc || {}),
        recovery_save: {
          project_id: saveSpec.projectId,
          cif_id: saveSpec.cifId,
          calc_type: saveSpec.calcType,
          parameters: saveSpec.parameters || {},
          tags: saveSpec.tags || [],
          input_content: saveSpec.inputContent || "",
        },
      },
    },
  };
}

function progressFromOutput(taskType: TaskType, output: string[], status: TaskStatus): ProgressState {
  let progress = defaultProgressState("Starting...");
  for (const line of output) {
    progress = progressReducer(taskType, line, progress);
  }

  if (status === "completed") {
    return {
      status: "complete",
      percent: 100,
      phase: "Complete",
    };
  }

  if (status === "cancelled") {
    return {
      status: "error",
      percent: progress.percent,
      phase: "Cancelled",
    };
  }

  if (status === "failed") {
    return {
      status: "error",
      percent: progress.percent,
      phase: "Failed",
    };
  }

  return progress;
}

function isPslibraryAtomicOrbitalLabelError(text: string): boolean {
  const normalized = text.toLowerCase();
  return normalized.includes("offset_atom_wfc")
    && (
      normalized.includes("does not contain labels for atomic orbitals")
      || normalized.includes("please add them by hand in the pseudo")
    );
}

function taskUsesRemotePseudos(hpc: HpcTaskMeta): boolean {
  return hpc.backend === "hpc"
    || Boolean(hpc.hpc_profile_id)
    || Boolean(hpc.remote_job_id)
    || Boolean(hpc.remote_workdir)
    || Boolean(hpc.remote_project_path);
}

function addPslibraryPseudoRepairHint(
  error: string | null,
  output: string[],
  hpc: HpcTaskMeta,
): string | null {
  const baseError = error ?? "Calculation failed.";
  const combinedText = `${baseError}\n${output.join("\n")}`;
  if (!isPslibraryAtomicOrbitalLabelError(combinedText)) {
    return error;
  }
  if (
    baseError.includes("Repair Local PSLibrary Pseudos")
    || baseError.includes("Repair Remote PSLibrary Pseudos")
  ) {
    return baseError;
  }

  const remote = taskUsesRemotePseudos(hpc);
  const settingsPath = remote ? "Settings > HPC" : "Settings > General";
  const buttonLabel = remote ? "Repair Remote PSLibrary Pseudos" : "Repair Local PSLibrary Pseudos";
  return `${baseError}\n\nThis matches a known PSLibrary UPF defect where the tenth atomic wavefunction tag is written as PP_CHI.1 instead of PP_CHI.10. Open ${settingsPath} and run ${buttonLabel}, then rerun the calculation.`;
}

function buildTaskState(
  taskId: string,
  taskType: TaskType,
  label: string,
  startedAt: string,
  status: TaskStatus,
  output: string[],
  result: any | null,
  error: string | null,
  hpc: HpcTaskMeta,
  visibleOnly = true,
): TaskState {
  const visibleOutput = visibleOnly
    ? buildVisibleOutputWindow(output)
    : {
      lines: output,
      text: output.join("\n"),
      totalLineCount: output.length,
    };

  return {
    taskId,
    taskType,
    label,
    startedAt,
    status,
    progress: progressFromOutput(taskType, output, status),
    output: visibleOutput.lines,
    outputText: visibleOutput.text,
    outputLineCount: visibleOutput.totalLineCount,
    result,
    error: addPslibraryPseudoRepairHint(error, output, hpc),
    hpc,
  };
}

function appendVisibleTaskOutput(task: TaskState, taskType: TaskType, line: string): TaskState {
  const visibleOutput = buildVisibleOutputWindow([...task.output, line]);
  return {
    ...task,
    progress: progressReducer(taskType, line, task.progress),
    output: visibleOutput.lines,
    outputText: visibleOutput.text,
    outputLineCount: task.outputLineCount + 1,
  };
}

function buildQueuedResult(
  taskType: TaskType,
  taskResult: any,
  outputText: string,
  parameters: Record<string, any>,
  taskStatus: TaskStatus = "completed",
  taskError: string | null = null,
) {
  const failed = taskStatus === "failed";
  const failureOutput = failed && taskError && !outputText.includes(taskError)
    ? `${outputText}${outputText.endsWith("\n") || outputText.length === 0 ? "" : "\n"}Error: ${taskError}`
    : outputText;

  if (taskType === "scf") {
    if (failed) {
      return {
        converged: false,
        total_energy: null,
        fermi_energy: null,
        n_scf_steps: null,
        wall_time_seconds: null,
        raw_output: failureOutput,
      };
    }
    if (taskResult && typeof taskResult === "object") {
      return taskResult;
    }
    return {
      converged: false,
      raw_output: outputText,
    };
  }

  if (taskType === "bands") {
    const fallbackEf = Number(parameters?.scf_fermi_energy);
    const fermiEnergy = Number.isFinite(fallbackEf) ? fallbackEf : null;
    return {
      converged: true,
      total_energy: null,
      fermi_energy: fermiEnergy,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failureOutput,
      band_data: failed ? null : taskResult,
    };
  }

  if (taskType === "dos") {
    const resultEf = Number(taskResult?.fermi_energy);
    const fallbackEf = Number(parameters?.scf_fermi_energy);
    const fermiEnergy = Number.isFinite(resultEf)
      ? resultEf
      : Number.isFinite(fallbackEf)
        ? fallbackEf
        : null;
    return {
      converged: true,
      total_energy: null,
      fermi_energy: fermiEnergy,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failureOutput,
      dos_data: failed ? null : taskResult,
    };
  }

  if (taskType === "fermi_surface") {
    const resultEf = Number(taskResult?.fermi_energy);
    const fallbackEf = Number(parameters?.scf_fermi_energy);
    const fermiEnergy = Number.isFinite(resultEf)
      ? resultEf
      : Number.isFinite(fallbackEf)
        ? fallbackEf
        : null;
    return {
      converged: true,
      total_energy: null,
      fermi_energy: fermiEnergy,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failureOutput,
    };
  }

  if (taskType === "hubbard_lrt") {
    return {
      converged: failed ? false : Boolean(taskResult?.converged ?? true),
      total_energy: null,
      fermi_energy: null,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failed ? failureOutput : taskResult?.raw_output ?? outputText,
      hubbard_lrt_data: failed ? null : taskResult,
    };
  }

  if (taskType === "wannier") {
    const resultEf = Number(taskResult?.band_data?.fermi_energy);
    const fallbackEf = Number(parameters?.scf_fermi_energy);
    const fermiEnergy = Number.isFinite(resultEf)
      ? resultEf
      : Number.isFinite(fallbackEf)
        ? fallbackEf
        : null;
    return {
      converged: taskResult?.convergence?.converged ?? true,
      total_energy: null,
      fermi_energy: fermiEnergy,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failureOutput,
      band_data: taskResult?.band_data ?? null,
      wannier_data: failed ? null : taskResult,
    };
  }

  if (taskType === "transport") {
    const resultEf = Number(taskResult?.reference_fermi_energy_ev);
    const fallbackEf = Number(parameters?.reference_fermi_energy_ev);
    const fermiEnergy = Number.isFinite(resultEf)
      ? resultEf
      : Number.isFinite(fallbackEf)
        ? fallbackEf
        : null;
    return {
      converged: true,
      total_energy: null,
      fermi_energy: fermiEnergy,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failureOutput,
      transport_data: failed ? null : taskResult,
    };
  }

  if (taskType === "epw") {
    return {
      converged: failed ? false : Boolean(taskResult?.result_summary?.completed),
      total_energy: null,
      fermi_energy: null,
      n_scf_steps: null,
      wall_time_seconds: null,
      raw_output: failureOutput,
      epw_data: failed ? null : taskResult,
    };
  }

  const converged = failed ? false : taskResult?.converged ?? true;
  return {
    converged,
    total_energy: null,
    fermi_energy: null,
    n_scf_steps: null,
    wall_time_seconds: null,
    raw_output: failed ? failureOutput : taskResult?.raw_output ?? outputText,
    phonon_data: {
      dos_data: failed ? null : taskResult?.dos_data ?? null,
      dispersion_data: failed ? null : taskResult?.dispersion_data ?? null,
    },
  };
}

function augmentQueuedParameters(taskType: TaskType, baseParameters: Record<string, any>, taskResult: any): Record<string, any> {
  const next = { ...baseParameters };

  if (taskType === "scf") {
    const mode = String(next.optimization_mode ?? next.calculation_mode ?? "").toLowerCase();
    const isOptimization = mode === "relax" || mode === "vcrelax";
    if (isOptimization) {
      const rawOutput = typeof taskResult?.raw_output === "string" ? taskResult.raw_output : "";
      const fallbackSource = isSavedStructureData(next.source_structure) ? next.source_structure : null;
      const fallbackOptimized = isSavedStructureData(next.optimized_structure) ? next.optimized_structure : null;
      const optimizedStructure = extractOptimizedStructure(rawOutput, fallbackOptimized || fallbackSource);
      if (optimizedStructure) {
        next.optimized_structure = optimizedStructure;
        const optimizedSummary = summarizeCell(optimizedStructure);
        if (optimizedSummary) {
          next.optimized_cell_summary = optimizedSummary;
        }
      }
    }
  } else if (taskType === "bands") {
    if (next.total_k_points == null && Number.isFinite(Number(taskResult?.n_kpoints))) {
      next.total_k_points = Number(taskResult.n_kpoints);
    }
    if (next.n_bands == null && Number.isFinite(Number(taskResult?.n_bands))) {
      next.n_bands = Number(taskResult.n_bands);
    }
  } else if (taskType === "dos") {
    if (next.n_points == null && Number.isFinite(Number(taskResult?.points))) {
      next.n_points = Number(taskResult.points);
    }
  } else if (taskType === "fermi_surface") {
    if (next.n_frmsf_files == null && Array.isArray(taskResult?.frmsf_files)) {
      next.n_frmsf_files = taskResult.frmsf_files.length;
    }
    if ((!Array.isArray(next.frmsf_files) || next.frmsf_files.length === 0) && Array.isArray(taskResult?.frmsf_files)) {
      next.frmsf_files = taskResult.frmsf_files.map((entry: any) => entry?.file_name).filter(Boolean);
    }
    if (next.primary_frmsf_file == null && typeof taskResult?.primary_file === "string") {
      next.primary_frmsf_file = taskResult.primary_file;
    }
    if (next.total_frmsf_bytes == null && Array.isArray(taskResult?.frmsf_files)) {
      next.total_frmsf_bytes = taskResult.frmsf_files.reduce(
        (sum: number, entry: any) => sum + (Number(entry?.size_bytes) || 0),
        0,
      );
    }
    // Backward compatibility for older fs.x / BXSF task payloads.
    if (next.n_bxsf_files == null && Array.isArray(taskResult?.bxsf_files)) {
      next.n_bxsf_files = taskResult.bxsf_files.length;
    }
    if ((!Array.isArray(next.bxsf_files) || next.bxsf_files.length === 0) && Array.isArray(taskResult?.bxsf_files)) {
      next.bxsf_files = taskResult.bxsf_files.map((entry: any) => entry?.file_name).filter(Boolean);
    }
    if (next.primary_bxsf_file == null && typeof taskResult?.primary_file === "string") {
      next.primary_bxsf_file = taskResult.primary_file;
    }
    if (next.total_bxsf_bytes == null && Array.isArray(taskResult?.bxsf_files)) {
      next.total_bxsf_bytes = taskResult.bxsf_files.reduce(
        (sum: number, entry: any) => sum + (Number(entry?.size_bytes) || 0),
        0,
      );
    }
  } else if (taskType === "phonon") {
    if (next.n_qpoints == null && Number.isFinite(Number(taskResult?.n_qpoints))) {
      next.n_qpoints = Number(taskResult.n_qpoints);
    }
    if (next.n_modes == null && Number.isFinite(Number(taskResult?.n_modes))) {
      next.n_modes = Number(taskResult.n_modes);
    }
  } else if (taskType === "hubbard_lrt") {
    if (next.u_value_count == null && Array.isArray(taskResult?.u_values)) {
      next.u_value_count = taskResult.u_values.length;
    }
    if (next.q_mesh == null && Array.isArray(taskResult?.q_mesh)) {
      next.q_mesh = taskResult.q_mesh;
    }
  } else if (taskType === "wannier") {
    if (next.total_spread == null && Number.isFinite(Number(taskResult?.total_spread))) {
      next.total_spread = Number(taskResult.total_spread);
    }
    if (next.n_wann == null && Number.isFinite(Number(taskResult?.num_wann))) {
      next.n_wann = Number(taskResult.num_wann);
    }
    if (next.n_bands == null && Number.isFinite(Number(taskResult?.num_bands))) {
      next.n_bands = Number(taskResult.num_bands);
    }
    if (next.total_k_points == null && Number.isFinite(Number(taskResult?.band_data?.n_kpoints))) {
      next.total_k_points = Number(taskResult.band_data.n_kpoints);
    }
  } else if (taskType === "transport") {
    if (next.seedname == null && typeof taskResult?.seedname === "string") {
      next.seedname = taskResult.seedname;
    }
    if (next.reference_fermi_energy_ev == null && Number.isFinite(Number(taskResult?.reference_fermi_energy_ev))) {
      next.reference_fermi_energy_ev = Number(taskResult.reference_fermi_energy_ev);
    }
    if (next.mu_points == null && Array.isArray(taskResult?.mu_values_ev)) {
      next.mu_points = taskResult.mu_values_ev.length;
    }
    if (next.temperature_points == null && Array.isArray(taskResult?.temperature_values_k)) {
      next.temperature_points = taskResult.temperature_values_k.length;
    }
    if (next.engine == null && typeof taskResult?.engine === "string") {
      next.engine = taskResult.engine;
    }
  } else if (taskType === "epw") {
    if (next.source_phonon_calc_id == null && typeof taskResult?.sources?.phonon?.calc_id === "string") {
      next.source_phonon_calc_id = taskResult.sources.phonon.calc_id;
    }
    if (next.source_wannier_calc_id == null && typeof taskResult?.sources?.wannier?.calc_id === "string") {
      next.source_wannier_calc_id = taskResult.sources.wannier.calc_id;
    }
    if (next.artifact_count == null && Array.isArray(taskResult?.artifacts)) {
      next.artifact_count = taskResult.artifacts.length;
    }
    if (next.parse_partial == null && typeof taskResult?.result_summary?.parse_partial === "boolean") {
      next.parse_partial = taskResult.result_summary.parse_partial;
    }
    if (next.epw_goals == null && taskResult?.goals && typeof taskResult.goals === "object") {
      next.epw_goals = taskResult.goals;
    }
    if (next.coarse_k_grid == null && Array.isArray(taskResult?.input?.coarse_k_mesh)) {
      next.coarse_k_grid = taskResult.input.coarse_k_mesh;
    }
    if (next.fine_k_grid == null && Array.isArray(taskResult?.input?.fine_k_mesh)) {
      next.fine_k_grid = taskResult.input.fine_k_mesh;
    }
    if (next.coarse_q_grid == null && Array.isArray(taskResult?.input?.coarse_q_mesh)) {
      next.coarse_q_grid = taskResult.input.coarse_q_mesh;
    }
    if (next.fine_q_grid == null && Array.isArray(taskResult?.input?.fine_q_mesh)) {
      next.fine_q_grid = taskResult.input.fine_q_mesh;
    }
  }

  return next;
}

async function saveRecoveredTaskResult(task: TaskState, recoverySave: HpcRecoverySaveSpec): Promise<void> {
  const outputLines = await invoke<string[]>("get_task_output", { taskId: task.taskId })
    .catch(() => task.output);
  const outputText = outputLines.join("\n");
  const failed = task.status === "failed";
  const parameters = augmentQueuedParameters(task.taskType, recoverySave.parameters || {}, task.result);
  if (failed) {
    parameters.run_status = "failed";
    parameters.failure_reason = task.error ?? "Task failed";
    parameters.failed_at = new Date().toISOString();
  }
  if (task.hpc.backend === "hpc") {
    parameters.execution_backend = "hpc";
    parameters.hpc_resource_type = task.hpc.hpc_resource_type ?? null;
    parameters.remote_job_id = task.hpc.remote_job_id ?? null;
    parameters.scheduler_state = task.hpc.scheduler_state ?? null;
    parameters.remote_node = task.hpc.remote_node ?? null;
    parameters.remote_workdir = task.hpc.remote_workdir ?? null;
    parameters.remote_project_path = task.hpc.remote_project_path ?? null;
    parameters.remote_storage_bytes = task.hpc.remote_storage_bytes ?? null;
    parameters.hpc_profile_id = task.hpc.hpc_profile_id ?? null;
    parameters.local_sync_dir = task.hpc.local_sync_dir ?? null;
  }
  const resultPayload = buildQueuedResult(task.taskType, task.result, outputText, parameters, task.status, task.error);
  const tags = [...(recoverySave.tags || [])];
  if (failed && !tags.includes("failed")) {
    tags.push("failed");
  }
  await invoke("save_calculation", {
    projectId: recoverySave.project_id,
    cifId: recoverySave.cif_id,
    calcData: {
      calc_type: recoverySave.calc_type,
      parameters,
      result: resultPayload,
      started_at: task.startedAt || new Date().toISOString(),
      completed_at: new Date().toISOString(),
      input_content: recoverySave.input_content || "",
      output_content: outputText,
      tags,
    },
    workingDir: task.hpc.local_sync_dir ?? null,
  });
}

export function useTaskContext(): TaskContextValue {
  const ctx = useContext(TaskContext);
  if (!ctx) {
    throw new Error("useTaskContext must be used within a TaskProvider");
  }
  return ctx;
}

export function TaskProvider({ children }: { children: React.ReactNode }) {
  const [tasks, setTasks] = useState<Map<string, TaskState>>(new Map());
  const [queueItems, setQueueItems] = useState<QueueItem[]>([]);
  const unlistenRefs = useRef<Map<string, UnlistenFn[]>>(new Map());
  const queueProcessingRef = useRef(false);
  const tasksRef = useRef(tasks);
  const queueRef = useRef(queueItems);
  const autoSavedRecoveryTasksRef = useRef<Set<string>>(new Set());

  useEffect(() => {
    tasksRef.current = tasks;
  }, [tasks]);

  useEffect(() => {
    queueRef.current = queueItems;
  }, [queueItems]);

  const maybeAutoSaveRecoveredTask = useCallback(async (task: TaskState) => {
    const recoverySave = task.hpc.recovery_save;
    const shouldSaveCompletedRecovery = task.hpc.headless_attached && task.status === "completed";
    const shouldSaveFailedArtifacts =
      task.status === "failed"
      && task.hpc.backend === "hpc"
      && Boolean(task.hpc.local_sync_dir);
    if (
      !recoverySave
      || (!shouldSaveCompletedRecovery && !shouldSaveFailedArtifacts)
      || autoSavedRecoveryTasksRef.current.has(task.taskId)
    ) {
      return;
    }
    autoSavedRecoveryTasksRef.current.add(task.taskId);
    try {
      await saveRecoveredTaskResult(task, recoverySave);
    } catch (e) {
      console.error("Failed to auto-save recovered HPC task:", e);
    }
  }, []);

  // Sync with backend on mount (handles app reload)
  useEffect(() => {
    void syncWithBackend();
    const intervalId = window.setInterval(() => {
      void syncWithBackend();
    }, 2000);
    return () => {
      window.clearInterval(intervalId);
      for (const fns of unlistenRefs.current.values()) {
        for (const fn of fns) fn();
      }
    };
  }, []);

  const subscribeToTask = useCallback((taskId: string, taskType: TaskType) => {
    const fns: UnlistenFn[] = [];
    const refreshTaskSnapshot = async (): Promise<TaskState | null> => {
      try {
        const [info, output] = await Promise.all([
          invoke<TaskInfo>("get_task_info", { taskId }).catch(() => null),
          invoke<string[]>("get_task_output", { taskId }).catch(() => null),
        ]);
        if (!output) return null;

        let refreshedTask: TaskState | null = null;
        setTasks((prev) => {
          const task = prev.get(taskId);
          if (!task) return prev;
          const next = new Map(prev);
          refreshedTask = buildTaskState(
            taskId,
            task.taskType,
            task.label,
            task.startedAt,
            info?.status ?? task.status,
            output,
            info?.result ?? task.result,
            info?.error ?? task.error,
            info ? taskInfoToHpcMeta(info) : task.hpc,
          );
          next.set(taskId, refreshedTask);
          return next;
        });
        return refreshedTask;
      } catch (e) {
        console.error("Failed to refresh task snapshot:", e);
        return null;
      }
    };

    listen<string>(`task-output:${taskId}`, (event) => {
      setTasks((prev) => {
        const task = prev.get(taskId);
        if (!task) return prev;
        const next = new Map(prev);
        next.set(taskId, appendVisibleTaskOutput(task, taskType, event.payload));
        return next;
      });
    }).then((fn) => fns.push(fn));

    listen<string>(`task-complete:${taskId}`, async () => {
      try {
        const info = await invoke<TaskInfo>("get_task_info", { taskId });
        let completedTask: TaskState | null = null;
        setTasks((prev) => {
          const task = prev.get(taskId);
          if (!task) return prev;
          const next = new Map(prev);
          completedTask = {
            ...task,
            status: "completed",
            result: info.result,
            hpc: taskInfoToHpcMeta(info),
            progress: {
              status: "complete",
              percent: 100,
              phase: "Complete",
            },
          };
          next.set(taskId, completedTask);
          return next;
        });
        if (completedTask) {
          void maybeAutoSaveRecoveredTask(completedTask);
        }
        void refreshTaskSnapshot();
      } catch (e) {
        console.error("Failed to get task info on completion:", e);
      }
    }).then((fn) => fns.push(fn));

    listen<string>(`task-status:${taskId}`, async (event) => {
      const payload = event.payload;
      if (payload.startsWith("failed:")) {
        const errorMsg = payload.slice(7);
        setTasks((prev) => {
          const task = prev.get(taskId);
          if (!task) return prev;
          const next = new Map(prev);
          next.set(taskId, {
            ...task,
            status: "failed",
            error: addPslibraryPseudoRepairHint(errorMsg, task.output, task.hpc),
            progress: {
              status: "error",
              percent: task.progress.percent,
              phase: "Failed",
            },
          });
          return next;
        });
        const refreshedTask = await refreshTaskSnapshot();
        if (refreshedTask) {
          void maybeAutoSaveRecoveredTask(refreshedTask);
        }
      }
    }).then((fn) => fns.push(fn));

    unlistenRefs.current.set(taskId, fns);
  }, [maybeAutoSaveRecoveredTask]);

  const startTaskInternal = useCallback(
    async (
      type: TaskType,
      params: Record<string, any>,
      label: string,
      saveSpec?: QueueSaveSpec | null,
    ): Promise<string> => {
      const effectiveParams = withHpcRecoverySave(params, saveSpec);
      const taskId = await invoke<string>(COMMAND_MAP[type], {
        ...effectiveParams,
        label,
      });

      const initialState: TaskState = {
        taskId,
        taskType: type,
        label,
        startedAt: new Date().toISOString(),
        status: "running",
        progress: defaultProgressState("Starting..."),
        output: [],
        outputText: "",
        outputLineCount: 0,
        result: null,
        error: null,
        hpc: taskInfoToHpcMeta({ recovery_save: effectiveParams.executionTarget?.hpc?.recovery_save ?? null }),
      };

      setTasks((prev) => {
        const next = new Map(prev);
        next.set(taskId, initialState);
        return next;
      });

      subscribeToTask(taskId, type);
      return taskId;
    },
    [subscribeToTask],
  );

  const enqueueTaskInternal = useCallback(
    (
      type: TaskType,
      params: Record<string, any>,
      label: string,
      saveSpec?: QueueSaveSpec | null,
    ): string => {
      const queueItemId = generateQueueItemId();
      const nextItem: QueueItem = {
        id: queueItemId,
        taskType: type,
        label,
        params,
        status: "queued",
        createdAt: new Date().toISOString(),
        startedAt: null,
        finishedAt: null,
        taskId: null,
        error: null,
        saveSpec: saveSpec ?? null,
      };
      setQueueItems((prev) => [...prev, nextItem]);
      return queueItemId;
    },
    [],
  );

  const enqueueTask = useCallback(
    (
      type: TaskType,
      params: Record<string, any>,
      label: string,
      saveSpec?: QueueSaveSpec | null,
    ): string => enqueueTaskInternal(type, params, label, saveSpec),
    [enqueueTaskInternal],
  );

  const cancelTask = useCallback(async (taskId: string) => {
    await invoke("cancel_task", { taskId });
    setTasks((prev) => {
      const task = prev.get(taskId);
      if (!task) return prev;
      const next = new Map(prev);
      next.set(taskId, {
        ...task,
        status: "cancelled",
        error: "Cancelled by user",
        progress: {
          status: "error",
          percent: task.progress.percent,
          phase: "Cancelled",
        },
      });
      return next;
    });

    const fns = unlistenRefs.current.get(taskId);
    if (fns) {
      for (const fn of fns) fn();
      unlistenRefs.current.delete(taskId);
    }
  }, []);

  const cancelQueueItem = useCallback(async (queueItemId: string) => {
    const item = queueRef.current.find((entry) => entry.id === queueItemId);
    if (!item) return;

    if ((item.status === "running" || item.status === "saving") && item.taskId) {
      try {
        await cancelTask(item.taskId);
      } catch (e) {
        console.error("Failed to cancel running queued task:", e);
      }
    }

    setQueueItems((prev) => prev.map((entry) => (
      entry.id === queueItemId
        ? {
          ...entry,
          status: "cancelled",
          error: entry.error ?? "Cancelled by user",
          finishedAt: new Date().toISOString(),
        }
        : entry
    )));
  }, [cancelTask]);

  const dismissTask = useCallback(async (taskId: string) => {
    await invoke("dismiss_task", { taskId });
    setTasks((prev) => {
      const next = new Map(prev);
      next.delete(taskId);
      return next;
    });
    const fns = unlistenRefs.current.get(taskId);
    if (fns) {
      for (const fn of fns) fn();
      unlistenRefs.current.delete(taskId);
    }
  }, []);

  const getTask = useCallback(
    (taskId: string) => tasks.get(taskId),
    [tasks],
  );

  const reconnectToTask = useCallback(async (taskId: string) => {
    try {
      const info = await invoke<TaskInfo>("get_task_info", { taskId });
      const taskType = normalizeTaskType(info.task_type);
      const hasSubscription = unlistenRefs.current.has(taskId);
      const localTask = tasksRef.current.get(taskId);
      const shouldLoadFullOutput = info.status !== "running" || !localTask || localTask.status !== "running" || !hasSubscription;

      if (shouldLoadFullOutput) {
        const output = await invoke<string[]>("get_task_output", { taskId });
        const loadedTask = buildTaskState(
          taskId,
          taskType,
          info.label,
          info.started_at,
          info.status,
          output,
          info.result,
          info.error,
          taskInfoToHpcMeta(info),
        );
        setTasks((prev) => {
          const next = new Map(prev);
          next.set(taskId, loadedTask);
          return next;
        });
        void maybeAutoSaveRecoveredTask(loadedTask);
      } else {
        setTasks((prev) => {
          const current = prev.get(taskId);
          if (!current) return prev;
          const next = new Map(prev);
          next.set(taskId, {
            ...current,
            label: info.label,
            startedAt: info.started_at,
            status: info.status,
            result: info.result ?? current.result,
            error: info.error ?? current.error,
            hpc: taskInfoToHpcMeta(info),
          });
          return next;
        });
      }

      if (info.status === "running" && !hasSubscription) {
        subscribeToTask(taskId, taskType);
      }
    } catch (e) {
      console.error("Failed to reconnect to task:", e);
    }
  }, [maybeAutoSaveRecoveredTask, subscribeToTask]);

  async function syncWithBackend() {
    try {
      const summaries = await invoke<TaskSummary[]>("list_running_tasks");
      for (const summary of summaries) {
        if (summary.status === "running") {
          await reconnectToTask(summary.task_id);
        } else {
          const info = await invoke<TaskInfo>("get_task_info", { taskId: summary.task_id });
          const output = await invoke<string[]>("get_task_output", { taskId: summary.task_id })
            .catch(() => []);
          const taskType = normalizeTaskType(summary.task_type);

          setTasks((prev) => {
            const next = new Map(prev);
            next.set(summary.task_id, buildTaskState(
              summary.task_id,
              taskType,
              summary.label,
              summary.started_at,
              info.status,
              output,
              info.result,
              info.error,
              taskInfoToHpcMeta(info),
            ));
            return next;
          });
        }
      }
    } catch (e) {
      console.error("Failed to sync tasks with backend:", e);
    }
  }

  const removeQueueItem = useCallback((queueItemId: string) => {
    setQueueItems((prev) => prev.filter((item) => item.id !== queueItemId));
  }, []);

  const moveQueueItem = useCallback((queueItemId: string, direction: "up" | "down") => {
    setQueueItems((prev) => {
      const index = prev.findIndex((item) => item.id === queueItemId);
      if (index < 0) return prev;
      const item = prev[index];
      if (item.status !== "queued") return prev;

      const queuedIndices = prev
        .map((entry, idx) => ({ entry, idx }))
        .filter(({ entry }) => entry.status === "queued")
        .map(({ idx }) => idx);
      const queuedPosition = queuedIndices.indexOf(index);
      if (queuedPosition < 0) return prev;

      const neighborPosition = direction === "up" ? queuedPosition - 1 : queuedPosition + 1;
      if (neighborPosition < 0 || neighborPosition >= queuedIndices.length) return prev;

      const neighborIndex = queuedIndices[neighborPosition];
      const next = [...prev];
      [next[index], next[neighborIndex]] = [next[neighborIndex], next[index]];
      return next;
    });
  }, []);

  const clearFinishedQueueItems = useCallback(() => {
    setQueueItems((prev) => prev.filter((item) => item.status === "queued" || item.status === "running" || item.status === "saving"));
  }, []);

  const waitForTaskCompletion = useCallback(async (taskId: string): Promise<TaskState> => {
    while (true) {
      const local = tasksRef.current.get(taskId);
      if (local && local.status !== "running") {
        try {
          const [info, output] = await Promise.all([
            invoke<TaskInfo>("get_task_info", { taskId }),
            invoke<string[]>("get_task_output", { taskId }),
          ]);
          const taskType = normalizeTaskType(info.task_type);
          const fullTask = buildTaskState(
            taskId,
            taskType,
            info.label,
            info.started_at,
            info.status,
            output,
            info.result,
            info.error,
            taskInfoToHpcMeta(info),
            false,
          );
          setTasks((prev) => {
            const next = new Map(prev);
            next.set(taskId, buildTaskState(
              taskId,
              taskType,
              info.label,
              info.started_at,
              info.status,
              output,
              info.result,
              info.error,
              taskInfoToHpcMeta(info),
            ));
            return next;
          });
          void maybeAutoSaveRecoveredTask(fullTask);
          return fullTask;
        } catch (e) {
          console.error("Failed to hydrate completed task output:", e);
          void maybeAutoSaveRecoveredTask(local);
          return local;
        }
      }

      try {
        const info = await invoke<TaskInfo>("get_task_info", { taskId });
        if (info.status !== "running") {
          const output = await invoke<string[]>("get_task_output", { taskId });
          const taskType = normalizeTaskType(info.task_type);
          const reconstructed = buildTaskState(
            taskId,
            taskType,
            info.label,
            info.started_at,
            info.status,
            output,
            info.result,
            info.error,
            taskInfoToHpcMeta(info),
            false,
          );

          setTasks((prev) => {
            const next = new Map(prev);
            next.set(taskId, buildTaskState(
              taskId,
              taskType,
              info.label,
              info.started_at,
              info.status,
              output,
              info.result,
              info.error,
              taskInfoToHpcMeta(info),
            ));
            return next;
          });

          void maybeAutoSaveRecoveredTask(reconstructed);
          return reconstructed;
        }
      } catch (e) {
        console.error("Failed to poll task status:", e);
      }

      await sleep(500);
    }
  }, [maybeAutoSaveRecoveredTask]);

  const saveQueuedTaskResult = useCallback(async (item: QueueItem, task: TaskState) => {
    const spec = item.saveSpec;
    if (!spec) return;

    const outputLines = await invoke<string[]>("get_task_output", { taskId: task.taskId })
      .catch(() => task.output);
    const outputText = outputLines.join("\n");
    const nowIso = new Date().toISOString();
    const taskResult = task.result;
    const parameters = augmentQueuedParameters(item.taskType, spec.parameters || {}, taskResult);
    if (task.hpc.backend === "hpc") {
      parameters.execution_backend = "hpc";
      parameters.hpc_resource_type = task.hpc.hpc_resource_type ?? null;
      parameters.remote_job_id = task.hpc.remote_job_id ?? null;
      parameters.scheduler_state = task.hpc.scheduler_state ?? null;
      parameters.remote_node = task.hpc.remote_node ?? null;
      parameters.remote_workdir = task.hpc.remote_workdir ?? null;
      parameters.remote_project_path = task.hpc.remote_project_path ?? null;
      parameters.remote_storage_bytes = task.hpc.remote_storage_bytes ?? null;
      parameters.hpc_profile_id = task.hpc.hpc_profile_id ?? null;
    }
    const resultPayload = buildQueuedResult(item.taskType, taskResult, outputText, parameters);

    await invoke("save_calculation", {
      projectId: spec.projectId,
      cifId: spec.cifId,
      calcData: {
        calc_type: spec.calcType,
        parameters,
        result: resultPayload,
        started_at: item.startedAt || task.startedAt || nowIso,
        completed_at: nowIso,
        input_content: spec.inputContent ?? "",
        output_content: outputText,
        tags: spec.tags ?? [],
      },
      workingDir: task.hpc.local_sync_dir ?? spec.workingDir ?? item.params.workingDir ?? null,
    });
  }, []);

  const processQueue = useCallback(async () => {
    if (queueProcessingRef.current) {
      return;
    }

    const nextQueued = queueRef.current.find((item) => item.status === "queued");
    if (!nextQueued) {
      return;
    }

    const hasRunningTask = Array.from(tasksRef.current.values()).some((task) => task.status === "running");
    if (hasRunningTask) {
      return;
    }

    queueProcessingRef.current = true;

    try {
      setQueueItems((prev) => prev.map((item) => (
        item.id === nextQueued.id
          ? {
            ...item,
            status: "running",
            startedAt: new Date().toISOString(),
            error: null,
          }
          : item
      )));

      let taskId: string;
      try {
        taskId = await startTaskInternal(nextQueued.taskType, nextQueued.params, nextQueued.label);
      } catch (e) {
        if (isBusyError(e)) {
          setQueueItems((prev) => prev.map((item) => (
            item.id === nextQueued.id
              ? { ...item, status: "queued", error: null }
              : item
          )));
          return;
        }

        setQueueItems((prev) => prev.map((item) => (
          item.id === nextQueued.id
            ? {
              ...item,
              status: "failed",
              finishedAt: new Date().toISOString(),
              error: String(e),
            }
            : item
        )));
        return;
      }

      setQueueItems((prev) => prev.map((item) => (
        item.id === nextQueued.id
          ? { ...item, taskId }
          : item
      )));

      const finalTask = await waitForTaskCompletion(taskId);
      if (finalTask.status === "completed") {
        const queueItem = queueRef.current.find((item) => item.id === nextQueued.id);
        if (queueItem?.saveSpec) {
          setQueueItems((prev) => prev.map((item) => (
            item.id === nextQueued.id
              ? { ...item, status: "saving" }
              : item
          )));

          try {
            await saveQueuedTaskResult(queueItem, finalTask);
          } catch (saveError) {
            setQueueItems((prev) => prev.map((item) => (
              item.id === nextQueued.id
                ? {
                  ...item,
                  status: "failed",
                  finishedAt: new Date().toISOString(),
                  error: `Failed to save result: ${saveError}`,
                }
                : item
            )));
            return;
          }
        }

        setQueueItems((prev) => prev.map((item) => (
          item.id === nextQueued.id
            ? {
              ...item,
              status: "completed",
              finishedAt: new Date().toISOString(),
              error: null,
            }
            : item
        )));
      } else if (finalTask.status === "cancelled") {
        setQueueItems((prev) => prev.map((item) => (
          item.id === nextQueued.id
            ? {
              ...item,
              status: "cancelled",
              finishedAt: new Date().toISOString(),
              error: finalTask.error ?? "Cancelled",
            }
            : item
        )));
      } else {
        setQueueItems((prev) => prev.map((item) => (
          item.id === nextQueued.id
            ? {
              ...item,
              status: "failed",
              finishedAt: new Date().toISOString(),
              error: finalTask.error ?? "Task failed",
            }
            : item
        )));
      }
    } finally {
      queueProcessingRef.current = false;
      if (queueRef.current.some((item) => item.status === "queued")) {
        window.setTimeout(() => {
          void processQueue();
        }, 0);
      }
    }
  }, [saveQueuedTaskResult, startTaskInternal, waitForTaskCompletion]);

  useEffect(() => {
    let unlisten: UnlistenFn | null = null;
    void listen("app-resumed", () => {
      void syncWithBackend();
      window.setTimeout(() => {
        void processQueue();
      }, 0);
    }).then((fn) => {
      unlisten = fn;
    });

    return () => {
      if (unlisten) {
        unlisten();
      }
    };
  }, [processQueue]);

  const waitForQueueTaskStart = useCallback(async (queueItemId: string): Promise<string> => {
    const startMs = Date.now();
    const maxWaitMs = 30_000;

    while (true) {
      const queueItem = queueRef.current.find((item) => item.id === queueItemId);
      if (!queueItem) {
        if (Date.now() - startMs > maxWaitMs) {
          throw new Error("Queued task was removed before it could start.");
        }
        await sleep(50);
        continue;
      }

      if (queueItem.taskId) {
        return queueItem.taskId;
      }

      if (queueItem.status === "failed" || queueItem.status === "cancelled") {
        throw new Error(queueItem.error ?? "Queued task failed before starting.");
      }

      await sleep(100);
    }
  }, []);

  const waitForQueueItemCompletion = useCallback(async (taskId: string): Promise<QueueItem | null> => {
    const startMs = Date.now();
    const maxWaitMs = 120_000;

    while (true) {
      const queueItem = queueRef.current.find((item) => item.taskId === taskId);
      if (!queueItem) {
        if (Date.now() - startMs > maxWaitMs) {
          return null;
        }
        await sleep(100);
        continue;
      }

      if (queueItem.status === "completed" || queueItem.status === "failed" || queueItem.status === "cancelled") {
        return queueItem;
      }

      await sleep(100);
    }
  }, []);

  const startTask = useCallback(
    async (
      type: TaskType,
      params: Record<string, any>,
      label: string,
      saveSpec?: QueueSaveSpec | null,
    ): Promise<string> => {
      if (isHpcStartParams(params)) {
        return startTaskInternal(type, params, label, saveSpec ?? null);
      }
      const queueItemId = enqueueTaskInternal(type, params, label, saveSpec ?? null);
      window.setTimeout(() => {
        void processQueue();
      }, 0);
      return waitForQueueTaskStart(queueItemId);
    },
    [enqueueTaskInternal, processQueue, startTaskInternal, waitForQueueTaskStart],
  );

  const hasQueuedItems = queueItems.some((item) => item.status === "queued");
  const hasRunningTask = useMemo(
    () => Array.from(tasks.values()).some((task) => task.status === "running"),
    [tasks],
  );

  useEffect(() => {
    if (hasQueuedItems && !hasRunningTask) {
      void processQueue();
    }
  }, [hasQueuedItems, hasRunningTask, processQueue]);

  const activeTasks = Array.from(tasks.values()).filter(
    (task) => task.status === "running" || task.status === "completed" || task.status === "failed" || task.status === "cancelled",
  );

  const queueSummary = useMemo<QueueSummary>(() => {
    const pending = queueItems.filter((item) => item.status === "queued" || item.status === "running" || item.status === "saving");
    if (pending.length === 0) {
      return { total: 0, activeIndex: null };
    }
    const runningIndex = pending.findIndex((item) => item.status === "running" || item.status === "saving");
    return {
      total: pending.length,
      activeIndex: runningIndex >= 0 ? runningIndex + 1 : 1,
    };
  }, [queueItems]);

  const value: TaskContextValue = {
    tasks,
    activeTasks,
    startTask,
    enqueueTask,
    queueItems,
    queueSummary,
    cancelQueueItem,
    removeQueueItem,
    moveQueueItem,
    clearFinishedQueueItems,
    cancelTask,
    dismissTask,
    getTask,
    waitForTaskCompletion,
    waitForQueueItemCompletion,
    reconnectToTask,
  };

  return (
    <TaskContext.Provider value={value}>
      {children}
    </TaskContext.Provider>
  );
}
