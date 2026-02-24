import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { useTaskContext } from "../lib/TaskContext";
import { openHpcActivityWindow } from "../lib/hpcConfig";

interface HpcActivityPanelProps {
  standalone?: boolean;
}

interface TimelineState {
  [key: string]: boolean;
}

interface TaskSummarySnapshot {
  task_id: string;
  status: "running" | "completed" | "failed" | "cancelled";
  backend?: string | null;
}

const STAGE_ORDER = [
  "Connecting",
  "Uploading",
  "Submitting",
  "Submitted",
  "Pending",
  "Running",
  "Completing",
  "Collecting",
  "Saved",
];

function deriveTimeline(output: string[]): TimelineState {
  const timeline: TimelineState = {};
  for (const stage of STAGE_ORDER) {
    timeline[stage] = false;
  }
  for (const line of output) {
    if (line.startsWith("HPC_STAGE|Connecting")) timeline.Connecting = true;
    if (line.startsWith("HPC_STAGE|Uploading")) timeline.Uploading = true;
    if (line.startsWith("HPC_STAGE|Submitting")) timeline.Submitting = true;
    if (line.startsWith("HPC_STAGE|Submitted")) timeline.Submitted = true;
    if (line.startsWith("HPC_STAGE|Collecting")) timeline.Collecting = true;
    if (line.startsWith("HPC_STAGE|Saved")) timeline.Saved = true;
    if (line.startsWith("HPC_SCHED|PENDING")) timeline.Pending = true;
    if (line.startsWith("HPC_SCHED|RUNNING")) timeline.Running = true;
    if (line.startsWith("HPC_SCHED|COMPLETING")) timeline.Completing = true;
  }
  return timeline;
}

function filterReadableOutput(lines: string[]): string[] {
  const rendered: string[] = [];
  for (const line of lines) {
    if (line.startsWith("HPC_CMD|") || line.startsWith("HPC_SCRIPT|") || line.startsWith("HPC_SCRIPT_BEGIN") || line.startsWith("HPC_SCRIPT_END")) {
      continue;
    }
    if (line.startsWith("HPC_STAGE|")) {
      rendered.push(`[stage] ${line.slice("HPC_STAGE|".length)}`);
      continue;
    }
    if (line.startsWith("HPC_SCHED|")) {
      rendered.push(`[scheduler] ${line.slice("HPC_SCHED|".length)}`);
      continue;
    }
    if (line.startsWith("HPC_WARNING|")) {
      rendered.push(`[warning] ${line.slice("HPC_WARNING|".length)}`);
      continue;
    }
    if (line.startsWith("HPC_TRANSFER|")) {
      continue;
    }
    rendered.push(line);
  }
  return rendered;
}

function extractScriptData(lines: string[]): { submitCommand: string | null; script: string | null } {
  let submitCommand: string | null = null;
  const scriptLines: string[] = [];

  for (const line of lines) {
    if (line.startsWith("HPC_CMD|")) {
      submitCommand = line.slice("HPC_CMD|".length).trim() || null;
    } else if (line.startsWith("HPC_SCRIPT|")) {
      scriptLines.push(line.slice("HPC_SCRIPT|".length));
    }
  }

  return {
    submitCommand,
    script: scriptLines.length > 0 ? scriptLines.join("\n") : null,
  };
}

export function HpcActivityPanel({ standalone = false }: HpcActivityPanelProps) {
  const { tasks, reconnectToTask } = useTaskContext();
  const tasksRef = useRef(tasks);
  const [isSyncing, setIsSyncing] = useState(false);
  const [lastSyncAt, setLastSyncAt] = useState<number | null>(null);
  const logRef = useRef<HTMLPreElement>(null);
  const followLogRef = useRef(true);
  const syncInFlightRef = useRef(false);

  useEffect(() => {
    tasksRef.current = tasks;
  }, [tasks]);

  const syncFromBackend = useCallback(async () => {
    if (syncInFlightRef.current) return;
    syncInFlightRef.current = true;
    setIsSyncing(true);
    try {
      const snapshots = await invoke<TaskSummarySnapshot[]>("list_running_tasks");
      const hpcSnapshots = snapshots.filter((item) => item.backend === "hpc");
      for (const snapshot of hpcSnapshots) {
        const local = tasksRef.current.get(snapshot.task_id);
        if (snapshot.status === "running" || !local || local.status !== snapshot.status) {
          await reconnectToTask(snapshot.task_id);
        }
      }
      setLastSyncAt(Date.now());
    } catch (e) {
      console.error("Failed to sync HPC activity panel:", e);
    } finally {
      setIsSyncing(false);
      syncInFlightRef.current = false;
    }
  }, [reconnectToTask]);

  useEffect(() => {
    void syncFromBackend();
    const interval = window.setInterval(() => {
      void syncFromBackend();
    }, 2000);
    return () => {
      window.clearInterval(interval);
    };
  }, [syncFromBackend]);

  const handleLogScroll = () => {
    const el = logRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followLogRef.current = distanceToBottom <= 16;
  };

  const hpcTasks = useMemo(
    () =>
      Array.from(tasks.values())
        .filter((task) => task.hpc.backend === "hpc")
        .sort((a, b) => b.startedAt.localeCompare(a.startedAt)),
    [tasks],
  );
  const [selectedTaskId, setSelectedTaskId] = useState<string | null>(null);

  const selectedTask = useMemo(() => {
    if (hpcTasks.length === 0) return null;
    if (selectedTaskId) {
      const found = hpcTasks.find((task) => task.taskId === selectedTaskId);
      if (found) return found;
    }
    return hpcTasks[0];
  }, [hpcTasks, selectedTaskId]);

  const timeline = useMemo(
    () => (selectedTask ? deriveTimeline(selectedTask.output) : deriveTimeline([])),
    [selectedTask],
  );

  const readableOutput = useMemo(
    () => (selectedTask ? filterReadableOutput(selectedTask.output).slice(-1600) : []),
    [selectedTask],
  );
  const scriptData = useMemo(
    () => (selectedTask ? extractScriptData(selectedTask.output) : { submitCommand: null, script: null }),
    [selectedTask],
  );
  const selectedStatus = selectedTask?.hpc.scheduler_state || selectedTask?.status || "idle";
  const selectedIsRunning = selectedTask?.status === "running";
  const syncLabel = isSyncing
    ? "Syncing..."
    : lastSyncAt
      ? `Updated ${new Date(lastSyncAt).toLocaleTimeString()}`
      : "Waiting for first sync";

  useEffect(() => {
    const el = logRef.current;
    if (!el || !followLogRef.current) return;
    el.scrollTop = el.scrollHeight;
  }, [selectedTask?.taskId, readableOutput.length]);

  if (hpcTasks.length === 0) {
    return (
      <section className={`hpc-activity-panel ${standalone ? "standalone" : "docked"}`}>
        <div className="hpc-activity-header">
          <div className="hpc-activity-header-main">
            <h3>Cluster Activity</h3>
            <div className="hpc-activity-header-actions">
              <button type="button" className="hpc-activity-refresh" onClick={() => void syncFromBackend()}>
                Refresh
              </button>
            </div>
          </div>
          <div className="hpc-activity-title-row">
            <span className="hpc-activity-sync">{syncLabel}</span>
          </div>
        </div>
        <p className="hpc-activity-empty">No HPC tasks are active yet.</p>
      </section>
    );
  }

  return (
    <section className={`hpc-activity-panel ${standalone ? "standalone" : "docked"}`}>
      <div className="hpc-activity-header">
        <div className="hpc-activity-header-main">
          <h3>Cluster Activity</h3>
          <div className="hpc-activity-header-actions">
            <button type="button" className="hpc-activity-refresh" onClick={() => void syncFromBackend()}>
              Refresh
            </button>
            {!standalone && (
              <button
                type="button"
                className="hpc-activity-popout"
                onClick={() => {
                  void openHpcActivityWindow();
                }}
              >
                Popout
              </button>
            )}
          </div>
        </div>
        <div className="hpc-activity-title-row">
          <span className={`hpc-activity-live ${selectedIsRunning ? "running" : "idle"}`}>
            {selectedIsRunning ? `Live · ${selectedStatus}` : selectedStatus}
          </span>
          <span className="hpc-activity-sync">{syncLabel}</span>
        </div>
      </div>

      <div className="hpc-activity-task-list">
        {hpcTasks.map((task) => (
          <button
            key={task.taskId}
            type="button"
            className={`hpc-activity-task-item ${selectedTask?.taskId === task.taskId ? "active" : ""}`}
            onClick={() => setSelectedTaskId(task.taskId)}
          >
            <span>{task.label}</span>
            <span>{task.hpc.remote_job_id ? `Job ${task.hpc.remote_job_id}` : task.taskType.toUpperCase()}</span>
          </button>
        ))}
      </div>

      {selectedTask && (
        <div className="hpc-activity-detail">
          <div className="hpc-activity-meta">
            <span>Job: {selectedTask.hpc.remote_job_id || "—"}</span>
            <span>State: {selectedTask.hpc.scheduler_state || selectedTask.status}</span>
            <span>Node: {selectedTask.hpc.remote_node || "—"}</span>
            <span>Task: {selectedTask.taskType.toUpperCase()}</span>
            <span>Started: {new Date(selectedTask.startedAt).toLocaleTimeString()}</span>
          </div>

          <pre ref={logRef} className="hpc-activity-log" onScroll={handleLogScroll}>
            {readableOutput.join("\n") || "Waiting for output..."}
          </pre>

          <div className="hpc-activity-foot">
            <div className="hpc-activity-timeline">
              {STAGE_ORDER.map((stage) => (
                <div key={stage} className={`hpc-activity-stage ${timeline[stage] ? "done" : ""}`}>
                  <span>{timeline[stage] ? "●" : "○"}</span>
                  <span>{stage}</span>
                </div>
              ))}
            </div>

            {(scriptData.submitCommand || scriptData.script) && (
              <div className="hpc-activity-command">
                {scriptData.submitCommand && (
                  <p>
                    <strong>Submit:</strong> <code>{scriptData.submitCommand}</code>
                  </p>
                )}
                {scriptData.script && (
                  <details>
                    <summary>Generated Slurm script</summary>
                    <pre>{scriptData.script}</pre>
                  </details>
                )}
              </div>
            )}
          </div>

        </div>
      )}
    </section>
  );
}
