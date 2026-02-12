import { createContext, useContext, useEffect, useState, useCallback, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { ProgressState, progressReducer, defaultProgressState } from "./qeProgress";

export type TaskStatus = "running" | "completed" | "failed" | "cancelled";
export type TaskType = "scf" | "bands" | "dos" | "phonon";

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
}

interface TaskSummary {
  task_id: string;
  task_type: string;
  label: string;
  started_at: string;
  status: TaskStatus;
}

interface TaskInfo {
  task_id: string;
  task_type: string;
  label: string;
  started_at: string;
  status: TaskStatus;
  result: any | null;
  error: string | null;
}

interface TaskContextValue {
  tasks: Map<string, TaskState>;
  activeTasks: TaskState[];
  startTask: (
    type: TaskType,
    params: Record<string, any>,
    label: string,
  ) => Promise<string>;
  cancelTask: (taskId: string) => Promise<void>;
  dismissTask: (taskId: string) => Promise<void>;
  getTask: (taskId: string) => TaskState | undefined;
  reconnectToTask: (taskId: string) => Promise<void>;
}

const TaskContext = createContext<TaskContextValue | null>(null);

export function useTaskContext(): TaskContextValue {
  const ctx = useContext(TaskContext);
  if (!ctx) {
    throw new Error("useTaskContext must be used within a TaskProvider");
  }
  return ctx;
}

export function TaskProvider({ children }: { children: React.ReactNode }) {
  const [tasks, setTasks] = useState<Map<string, TaskState>>(new Map());
  const unlistenRefs = useRef<Map<string, UnlistenFn[]>>(new Map());

  // Sync with backend on mount (handles app reload)
  useEffect(() => {
    syncWithBackend();
    return () => {
      // Cleanup all listeners
      for (const fns of unlistenRefs.current.values()) {
        for (const fn of fns) fn();
      }
    };
  }, []);

  async function syncWithBackend() {
    try {
      const summaries = await invoke<TaskSummary[]>("list_running_tasks");
      for (const summary of summaries) {
        if (summary.status === "running") {
          await reconnectToTask(summary.task_id);
        } else {
          // Load completed/failed task info
          const info = await invoke<TaskInfo>("get_task_info", { taskId: summary.task_id });
          setTasks((prev) => {
            const next = new Map(prev);
            next.set(summary.task_id, {
              taskId: summary.task_id,
              taskType: summary.task_type as TaskType,
              label: summary.label,
              startedAt: summary.started_at,
              status: info.status,
              progress: {
                status: info.status === "completed" ? "complete" : "error",
                percent: info.status === "completed" ? 100 : null,
                phase: info.status === "completed" ? "Complete" : "Failed",
              },
              output: [],
              result: info.result,
              error: info.error,
            });
            return next;
          });
        }
      }
    } catch (e) {
      console.error("Failed to sync tasks with backend:", e);
    }
  }

  function subscribeToTask(taskId: string, taskType: TaskType) {
    const fns: UnlistenFn[] = [];

    // Listen for output lines
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

    // Listen for completion
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
            progress: {
              status: "complete",
              percent: 100,
              phase: "Complete",
            },
          });
          return next;
        });
      } catch (e) {
        console.error("Failed to get task info on completion:", e);
      }
    }).then((fn) => fns.push(fn));

    // Listen for failure/status changes
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
      }
    }).then((fn) => fns.push(fn));

    unlistenRefs.current.set(taskId, fns);
  }

  const startTask = useCallback(
    async (
      type: TaskType,
      params: Record<string, any>,
      label: string,
    ): Promise<string> => {
      const commandMap: Record<TaskType, string> = {
        scf: "start_scf_calculation",
        bands: "start_bands_calculation",
        dos: "start_dos_calculation",
        phonon: "start_phonon_calculation",
      };

      const taskId = await invoke<string>(commandMap[type], {
        ...params,
        label,
      });

      // Create initial task state
      const initialState: TaskState = {
        taskId,
        taskType: type,
        label,
        startedAt: new Date().toISOString(),
        status: "running",
        progress: defaultProgressState(
          type === "scf"
            ? "Starting..."
            : type === "bands"
              ? "Starting..."
              : type === "dos"
                ? "Starting..."
              : "Starting...",
        ),
        output: [],
        result: null,
        error: null,
      };

      setTasks((prev) => {
        const next = new Map(prev);
        next.set(taskId, initialState);
        return next;
      });

      subscribeToTask(taskId, type);
      return taskId;
    },
    [],
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
    // Clean up listeners
    const fns = unlistenRefs.current.get(taskId);
    if (fns) {
      for (const fn of fns) fn();
      unlistenRefs.current.delete(taskId);
    }
  }, []);

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

      // Replay progress from buffered output
      const taskType = info.task_type as TaskType;
      let progress = defaultProgressState("Starting...");
      for (const line of output) {
        progress = progressReducer(taskType, line, progress);
      }

      // If task already completed/failed, set final state
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
        });
        return next;
      });

      // Only subscribe to events if the task is still running
      if (info.status === "running") {
        subscribeToTask(taskId, taskType);
      }
    } catch (e) {
      console.error("Failed to reconnect to task:", e);
    }
  }, []);

  const activeTasks = Array.from(tasks.values()).filter(
    (t) => t.status === "running" || t.status === "completed" || t.status === "failed" || t.status === "cancelled",
  );

  const value: TaskContextValue = {
    tasks,
    activeTasks,
    startTask,
    cancelTask,
    dismissTask,
    getTask,
    reconnectToTask,
  };

  return (
    <TaskContext.Provider value={value}>
      {children}
    </TaskContext.Provider>
  );
}
