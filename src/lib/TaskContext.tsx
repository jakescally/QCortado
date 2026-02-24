import { createContext, useContext, useEffect, useMemo, useState, useCallback, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { ProgressState, progressReducer, defaultProgressState } from "./qeProgress";
import { extractOptimizedStructure, isSavedStructureData, summarizeCell } from "./optimizedStructure";
import { HpcTaskMeta } from "./types";

export type TaskStatus = "running" | "completed" | "failed" | "cancelled";
export type TaskType = "scf" | "bands" | "dos" | "fermi_surface" | "phonon";
export type QueueItemStatus = "queued" | "running" | "saving" | "completed" | "failed" | "cancelled";

export interface TaskState {
  taskId: string;
  taskType: TaskType;
  label: string;
  startedAt: string;
  status: TaskStatus;
  progress: ProgressState;
  output: string[];
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
  remote_job_id?: string | null;
  scheduler_state?: string | null;
  remote_node?: string | null;
  remote_workdir?: string | null;
  remote_project_path?: string | null;
  remote_storage_bytes?: number | null;
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
  remote_job_id?: string | null;
  scheduler_state?: string | null;
  remote_node?: string | null;
  remote_workdir?: string | null;
  remote_project_path?: string | null;
  remote_storage_bytes?: number | null;
}

interface QueueSaveSpec {
  projectId: string;
  cifId: string;
  workingDir?: string | null;
  calcType: "scf" | "bands" | "dos" | "fermi_surface" | "phonon" | "optimization";
  parameters: Record<string, any>;
  tags?: string[];
  inputContent?: string;
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
  phonon: "start_phonon_calculation",
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

function normalizeTaskType(taskType: string): TaskType {
  if (taskType === "scf" || taskType === "bands" || taskType === "dos" || taskType === "fermi_surface" || taskType === "phonon") {
    return taskType;
  }
  return "scf";
}

function taskInfoToHpcMeta(info: Partial<TaskInfo> | Partial<TaskSummary>): HpcTaskMeta {
  return {
    backend: info.backend ?? null,
    remote_job_id: info.remote_job_id ?? null,
    scheduler_state: info.scheduler_state ?? null,
    remote_node: info.remote_node ?? null,
    remote_workdir: info.remote_workdir ?? null,
    remote_project_path: info.remote_project_path ?? null,
    remote_storage_bytes: info.remote_storage_bytes ?? null,
  };
}

function buildQueuedResult(taskType: TaskType, taskResult: any, outputText: string, parameters: Record<string, any>) {
  if (taskType === "scf") {
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
      raw_output: outputText,
      band_data: taskResult,
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
      raw_output: outputText,
      dos_data: taskResult,
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
      raw_output: outputText,
    };
  }

  const converged = taskResult?.converged ?? true;
  return {
    converged,
    total_energy: null,
    fermi_energy: null,
    n_scf_steps: null,
    wall_time_seconds: null,
    raw_output: taskResult?.raw_output ?? outputText,
    phonon_data: {
      dos_data: taskResult?.dos_data ?? null,
      dispersion_data: taskResult?.dispersion_data ?? null,
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
  }

  return next;
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

  useEffect(() => {
    tasksRef.current = tasks;
  }, [tasks]);

  useEffect(() => {
    queueRef.current = queueItems;
  }, [queueItems]);

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
    const refreshTaskSnapshot = async () => {
      try {
        const [info, output] = await Promise.all([
          invoke<TaskInfo>("get_task_info", { taskId }).catch(() => null),
          invoke<string[]>("get_task_output", { taskId }).catch(() => null),
        ]);
        if (!output) return;

        setTasks((prev) => {
          const task = prev.get(taskId);
          if (!task) return prev;

          let nextProgress = defaultProgressState("Starting...");
          for (const line of output) {
            nextProgress = progressReducer(taskType, line, nextProgress);
          }

          const nextStatus = info?.status ?? task.status;
          if (nextStatus !== "running") {
            nextProgress = {
              status: nextStatus === "completed" ? "complete" : "error",
              percent: nextStatus === "completed" ? 100 : nextProgress.percent,
              phase: nextStatus === "completed"
                ? "Complete"
                : nextStatus === "cancelled"
                  ? "Cancelled"
                  : "Failed",
            };
          }

          const next = new Map(prev);
          next.set(taskId, {
            ...task,
            status: nextStatus,
            output,
            progress: nextProgress,
            result: info?.result ?? task.result,
            error: info?.error ?? task.error,
            hpc: info ? taskInfoToHpcMeta(info) : task.hpc,
          });
          return next;
        });
      } catch (e) {
        console.error("Failed to refresh task snapshot:", e);
      }
    };

    listen<string>(`task-output:${taskId}`, (event) => {
      setTasks((prev) => {
        const task = prev.get(taskId);
        if (!task) return prev;
        const next = new Map(prev);
        const newProgress = progressReducer(taskType, event.payload, task.progress);
        next.set(taskId, {
          ...task,
          output: [...task.output, event.payload],
          progress: newProgress,
        });
        return next;
      });
    }).then((fn) => fns.push(fn));

    listen<string>(`task-complete:${taskId}`, async () => {
      try {
        const info = await invoke<TaskInfo>("get_task_info", { taskId });
        setTasks((prev) => {
          const task = prev.get(taskId);
          if (!task) return prev;
          const next = new Map(prev);
          next.set(taskId, {
            ...task,
            status: "completed",
            result: info.result,
            hpc: taskInfoToHpcMeta(info),
            progress: {
              status: "complete",
              percent: 100,
              phase: "Complete",
            },
          });
          return next;
        });
        void refreshTaskSnapshot();
      } catch (e) {
        console.error("Failed to get task info on completion:", e);
      }
    }).then((fn) => fns.push(fn));

    listen<string>(`task-status:${taskId}`, (event) => {
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
            error: errorMsg,
            progress: {
              status: "error",
              percent: task.progress.percent,
              phase: "Failed",
            },
          });
          return next;
        });
        void refreshTaskSnapshot();
      }
    }).then((fn) => fns.push(fn));

    unlistenRefs.current.set(taskId, fns);
  }, []);

  const startTaskInternal = useCallback(
    async (
      type: TaskType,
      params: Record<string, any>,
      label: string,
    ): Promise<string> => {
      const taskId = await invoke<string>(COMMAND_MAP[type], {
        ...params,
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
        result: null,
        error: null,
        hpc: {},
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
      const output = await invoke<string[]>("get_task_output", { taskId });

      const taskType = normalizeTaskType(info.task_type);
      let progress = defaultProgressState("Starting...");
      for (const line of output) {
        progress = progressReducer(taskType, line, progress);
      }

      if (info.status !== "running") {
        progress = {
          status: info.status === "completed" ? "complete" : "error",
          percent: info.status === "completed" ? 100 : progress.percent,
          phase: info.status === "completed" ? "Complete" : "Failed",
        };
      }

      setTasks((prev) => {
        const next = new Map(prev);
        next.set(taskId, {
          taskId,
          taskType,
          label: info.label,
          startedAt: info.started_at,
          status: info.status,
          progress,
          output,
          result: info.result,
          error: info.error,
          hpc: taskInfoToHpcMeta(info),
        });
        return next;
      });

      if (info.status === "running" && !unlistenRefs.current.has(taskId)) {
        subscribeToTask(taskId, taskType);
      }
    } catch (e) {
      console.error("Failed to reconnect to task:", e);
    }
  }, [subscribeToTask]);

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
          let progress = defaultProgressState("Starting...");
          for (const line of output) {
            progress = progressReducer(taskType, line, progress);
          }
          if (info.status === "completed") {
            progress = {
              status: "complete",
              percent: 100,
              phase: "Complete",
            };
          } else if (info.status === "cancelled") {
            progress = {
              status: "error",
              percent: progress.percent,
              phase: "Cancelled",
            };
          } else {
            progress = {
              status: "error",
              percent: progress.percent,
              phase: "Failed",
            };
          }

          setTasks((prev) => {
            const next = new Map(prev);
            next.set(summary.task_id, {
              taskId: summary.task_id,
              taskType,
              label: summary.label,
              startedAt: summary.started_at,
              status: info.status,
              progress,
              output,
              result: info.result,
              error: info.error,
              hpc: taskInfoToHpcMeta(info),
            });
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
        return local;
      }

      try {
        const info = await invoke<TaskInfo>("get_task_info", { taskId });
        if (info.status !== "running") {
          const output = await invoke<string[]>("get_task_output", { taskId });
          const taskType = normalizeTaskType(info.task_type);
          let progress = defaultProgressState("Starting...");
          for (const line of output) {
            progress = progressReducer(taskType, line, progress);
          }
          if (info.status === "completed") {
            progress = {
              status: "complete",
              percent: 100,
              phase: "Complete",
            };
          } else {
            progress = {
              status: "error",
              percent: progress.percent,
              phase: info.status === "cancelled" ? "Cancelled" : "Failed",
            };
          }

          const reconstructed: TaskState = {
            taskId,
            taskType,
            label: info.label,
            startedAt: info.started_at,
            status: info.status,
            progress,
            output,
            result: info.result,
            error: info.error,
            hpc: taskInfoToHpcMeta(info),
          };

          setTasks((prev) => {
            const next = new Map(prev);
            next.set(taskId, reconstructed);
            return next;
          });

          return reconstructed;
        }
      } catch (e) {
        console.error("Failed to poll task status:", e);
      }

      await sleep(500);
    }
  }, []);

  const saveQueuedTaskResult = useCallback(async (item: QueueItem, task: TaskState) => {
    const spec = item.saveSpec;
    if (!spec) return;

    const outputLines = task.output.length > 0
      ? task.output
      : await invoke<string[]>("get_task_output", { taskId: task.taskId });
    const outputText = outputLines.join("\n");
    const nowIso = new Date().toISOString();
    const taskResult = task.result;
    const parameters = augmentQueuedParameters(item.taskType, spec.parameters || {}, taskResult);
    if (task.hpc.backend === "hpc") {
      parameters.execution_backend = "hpc";
      parameters.remote_job_id = task.hpc.remote_job_id ?? null;
      parameters.scheduler_state = task.hpc.scheduler_state ?? null;
      parameters.remote_node = task.hpc.remote_node ?? null;
      parameters.remote_workdir = task.hpc.remote_workdir ?? null;
      parameters.remote_project_path = task.hpc.remote_project_path ?? null;
      parameters.remote_storage_bytes = task.hpc.remote_storage_bytes ?? null;
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
      workingDir: spec.workingDir ?? item.params.workingDir ?? null,
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
      const queueItemId = enqueueTaskInternal(type, params, label, saveSpec ?? null);
      window.setTimeout(() => {
        void processQueue();
      }, 0);
      return waitForQueueTaskStart(queueItemId);
    },
    [enqueueTaskInternal, processQueue, waitForQueueTaskStart],
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
