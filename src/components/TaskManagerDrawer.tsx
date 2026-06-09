import { useEffect, useMemo, useState } from "react";
import {
  ArrowDown,
  ArrowUp,
  ExternalLink,
  ListRestart,
  RefreshCw,
  Square,
  Trash2,
  X,
} from "lucide-react";
import { openHpcActivityWindow } from "../lib/hpcConfig";
import {
  buildTaskManagerEntries,
  filterTaskManagerEntries,
  summarizeTaskManagerEntries,
} from "../lib/taskManager";
import type { TaskManagerEntry, TaskManagerFilter } from "../lib/taskManager";
import { useTaskContext } from "../lib/TaskContext";
import { useOverlayDrawer } from "../lib/useOverlayDrawer";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";

interface TaskManagerDrawerProps {
  isOpen: boolean;
  requestedFilter: TaskManagerFilter;
  focusedTaskId: string | null;
  onClose: () => void;
  onNavigateToTask?: (taskId: string, taskType: string) => boolean;
}

const FILTERS: Array<{ id: TaskManagerFilter; label: string }> = [
  { id: "active", label: "Active" },
  { id: "all", label: "All" },
  { id: "finished", label: "Finished" },
  { id: "hpc", label: "HPC" },
];

function formatDate(value: string | null): string {
  if (!value) return "Not yet";
  const parsed = new Date(value);
  return Number.isNaN(parsed.getTime()) ? value : parsed.toLocaleString();
}

function statusLabel(entry: TaskManagerEntry): string {
  if (entry.status === "saving") return "Saving";
  if (entry.status === "completed") return "Complete";
  return entry.status.charAt(0).toUpperCase() + entry.status.slice(1);
}

function extractScriptData(output: string[]): { submitCommand: string | null; script: string | null } {
  let submitCommand: string | null = null;
  const script: string[] = [];
  for (const line of output) {
    if (line.startsWith("HPC_CMD|")) submitCommand = line.slice("HPC_CMD|".length).trim() || null;
    if (line.startsWith("HPC_SCRIPT|")) script.push(line.slice("HPC_SCRIPT|".length));
  }
  return { submitCommand, script: script.length > 0 ? script.join("\n") : null };
}

export function TaskManagerDrawer({
  isOpen,
  requestedFilter,
  focusedTaskId,
  onClose,
  onNavigateToTask,
}: TaskManagerDrawerProps) {
  const {
    tasks,
    queueItems,
    moveQueueItem,
    removeQueueItem,
    cancelQueueItem,
    cancelTask,
    dismissTask,
    clearFinishedTasks,
    reconnectToTask,
  } = useTaskContext();
  const drawerRef = useOverlayDrawer<HTMLElement>(isOpen, onClose);
  const [filter, setFilter] = useState<TaskManagerFilter>(requestedFilter);
  const [selectedEntryId, setSelectedEntryId] = useState<string | null>(null);
  const [confirmingCancelId, setConfirmingCancelId] = useState<string | null>(null);
  const [navigationError, setNavigationError] = useState<string | null>(null);

  const entries = useMemo(
    () => buildTaskManagerEntries(tasks.values(), queueItems),
    [queueItems, tasks],
  );
  const summary = useMemo(() => summarizeTaskManagerEntries(entries), [entries]);
  const filteredEntries = useMemo(() => filterTaskManagerEntries(entries, filter), [entries, filter]);
  const selectedEntry = useMemo(
    () => filteredEntries.find((entry) => entry.id === selectedEntryId)
      ?? filteredEntries.find((entry) => entry.taskId === focusedTaskId)
      ?? filteredEntries[0]
      ?? null,
    [filteredEntries, focusedTaskId, selectedEntryId],
  );

  useEffect(() => {
    if (isOpen) setFilter(requestedFilter);
  }, [isOpen, requestedFilter]);

  useEffect(() => {
    if (!isOpen || !focusedTaskId) return;
    const focused = entries.find((entry) => entry.taskId === focusedTaskId);
    if (focused) setSelectedEntryId(focused.id);
  }, [entries, focusedTaskId, isOpen]);

  useEffect(() => {
    setConfirmingCancelId(null);
    setNavigationError(null);
  }, [selectedEntry?.id]);

  const queuedEntries = entries.filter((entry) => entry.group === "queued");
  const selectedQueuedIndex = selectedEntry
    ? queuedEntries.findIndex((entry) => entry.id === selectedEntry.id)
    : -1;
  const selectedScript = selectedEntry?.task
    ? extractScriptData(selectedEntry.task.output)
    : { submitCommand: null, script: null };

  async function handleCancel(entry: TaskManagerEntry) {
    if (confirmingCancelId !== entry.id) {
      setConfirmingCancelId(entry.id);
      window.setTimeout(() => setConfirmingCancelId((current) => current === entry.id ? null : current), 3000);
      return;
    }
    setConfirmingCancelId(null);
    if (entry.queueItemId) {
      await cancelQueueItem(entry.queueItemId);
    } else if (entry.taskId) {
      await cancelTask(entry.taskId);
    }
  }

  async function handleRemove(entry: TaskManagerEntry) {
    if (entry.queueItemId) removeQueueItem(entry.queueItemId);
    if (entry.taskId && entry.task?.status !== "running") await dismissTask(entry.taskId);
  }

  function handleOpenRun(entry: TaskManagerEntry) {
    if (!entry.taskId || !onNavigateToTask) return;
    const opened = onNavigateToTask(entry.taskId, entry.taskType);
    if (!opened) {
      setNavigationError("This run no longer has enough project context to reopen its workflow. Its status, log, and artifacts remain available here.");
    }
  }

  return (
    <div className={`app-drawer-layer right ${isOpen ? "open" : ""}`} aria-hidden={!isOpen} inert={!isOpen}>
      <aside ref={drawerRef} className="app-drawer task-manager-drawer" role="dialog" aria-modal="false" aria-label="Task manager">
        <div className="app-drawer-header task-manager-header">
          <div>
            <span className="app-drawer-eyebrow">Session activity</span>
            <h2>Task Manager</h2>
          </div>
          <div className="task-manager-header-actions">
            <button
              type="button"
              className="chrome-icon-btn"
              onClick={() => void clearFinishedTasks()}
              disabled={summary.finished === 0}
              aria-label="Clear finished tasks"
              title="Clear finished tasks"
            >
              <ListRestart size={18} />
            </button>
            <button type="button" className="chrome-icon-btn" onClick={onClose} aria-label="Close task manager">
              <X size={19} />
            </button>
          </div>
        </div>

        <div className="task-manager-filters" role="tablist" aria-label="Task filters">
          {FILTERS.map((option) => (
            <button
              key={option.id}
              type="button"
              className={filter === option.id ? "active" : ""}
              onClick={() => setFilter(option.id)}
              role="tab"
              aria-selected={filter === option.id}
            >
              {option.label}
              <span>{
                option.id === "active" ? summary.running + summary.queued
                  : option.id === "finished" ? summary.finished
                    : option.id === "hpc" ? summary.hpc
                      : summary.total
              }</span>
            </button>
          ))}
        </div>

        <div className="task-manager-layout">
          <div className="task-manager-list" aria-label="Tasks">
            {filteredEntries.length === 0 ? (
              <div className="task-manager-empty">
                <strong>No tasks in this view</strong>
                <span>New calculations and queued work will appear here.</span>
              </div>
            ) : filteredEntries.map((entry) => (
              <button
                key={entry.id}
                type="button"
                className={`task-manager-list-item ${selectedEntry?.id === entry.id ? "active" : ""}`}
                onClick={() => setSelectedEntryId(entry.id)}
              >
                <span className={`task-state-dot ${entry.status}`} aria-hidden="true" />
                <span className="task-manager-list-copy">
                  <strong>{entry.label}</strong>
                  <span>{entry.taskType.toUpperCase()} · {statusLabel(entry)}{entry.isHpc ? " · HPC" : ""}</span>
                </span>
                {entry.task?.progress.percent != null && entry.group === "running" && (
                  <span className="task-manager-list-percent">{Math.round(entry.task.progress.percent)}%</span>
                )}
              </button>
            ))}
          </div>

          <div className="task-manager-detail">
            {!selectedEntry ? (
              <div className="task-manager-empty">
                <strong>No task selected</strong>
              </div>
            ) : (
              <>
                <div className="task-detail-heading">
                  <div>
                    <span className={`task-detail-status ${selectedEntry.status}`}>{statusLabel(selectedEntry)}</span>
                    <h3>{selectedEntry.label}</h3>
                    <p>{selectedEntry.taskType.toUpperCase()} · {selectedEntry.isHpc ? "HPC" : "Local"}</p>
                  </div>
                  <div className="task-detail-actions">
                    {selectedEntry.taskId && (
                      <button type="button" onClick={() => handleOpenRun(selectedEntry)}>
                        <ExternalLink size={15} /> Open run
                      </button>
                    )}
                    {selectedEntry.taskId && (
                      <button type="button" onClick={() => void reconnectToTask(selectedEntry.taskId!)}>
                        <RefreshCw size={15} /> Refresh
                      </button>
                    )}
                    {(selectedEntry.group === "running" || selectedEntry.group === "queued") && (
                      <button type="button" className="danger" onClick={() => void handleCancel(selectedEntry)}>
                        <Square size={14} />
                        {confirmingCancelId === selectedEntry.id ? "Confirm cancel" : "Cancel"}
                      </button>
                    )}
                    {selectedEntry.group === "finished" && (
                      <button type="button" onClick={() => void handleRemove(selectedEntry)}>
                        <Trash2 size={15} /> Dismiss
                      </button>
                    )}
                  </div>
                </div>

                <div className="task-detail-summary">
                  <span>Created<strong>{formatDate(selectedEntry.createdAt)}</strong></span>
                  <span>Started<strong>{formatDate(selectedEntry.startedAt)}</strong></span>
                  {selectedEntry.group === "running" && selectedEntry.startedAt && (
                    <span>Elapsed<ElapsedTimer startedAt={selectedEntry.startedAt} isRunning label="" /></span>
                  )}
                  {selectedEntry.group === "queued" && selectedQueuedIndex >= 0 && (
                    <span>Queue position<strong>{selectedQueuedIndex + 1} of {queuedEntries.length}</strong></span>
                  )}
                  {selectedEntry.isHpc && selectedEntry.task?.hpc.remote_job_id && (
                    <span>Job<strong>{selectedEntry.task.hpc.remote_job_id}</strong></span>
                  )}
                </div>

                {selectedEntry.task?.progress.phase && selectedEntry.group === "running" && (
                  <div className="task-detail-progress">
                    <div className={`task-detail-progress-track ${selectedEntry.task.progress.percent == null ? "indeterminate" : ""}`}>
                      <span style={selectedEntry.task.progress.percent == null ? undefined : { width: `${selectedEntry.task.progress.percent}%` }} />
                    </div>
                    <div><strong>{selectedEntry.task.progress.phase}</strong><span>{selectedEntry.task.progress.detail}</span></div>
                  </div>
                )}

                {selectedEntry.error && <div className="task-detail-error">{selectedEntry.error}</div>}
                {navigationError && <div className="task-detail-error">{navigationError}</div>}

                {selectedEntry.queueItemId && selectedEntry.group === "queued" && (
                  <div className="task-queue-controls">
                    <button
                      type="button"
                      disabled={selectedQueuedIndex <= 0}
                      onClick={() => moveQueueItem(selectedEntry.queueItemId!, "up")}
                    >
                      <ArrowUp size={15} /> Move up
                    </button>
                    <button
                      type="button"
                      disabled={selectedQueuedIndex < 0 || selectedQueuedIndex >= queuedEntries.length - 1}
                      onClick={() => moveQueueItem(selectedEntry.queueItemId!, "down")}
                    >
                      <ArrowDown size={15} /> Move down
                    </button>
                    <button type="button" onClick={() => void handleRemove(selectedEntry)}>
                      <Trash2 size={15} /> Remove
                    </button>
                  </div>
                )}

                {selectedEntry.isHpc && selectedEntry.task && (
                  <div className="task-hpc-detail">
                    <div className="task-hpc-toolbar">
                      <span>{selectedEntry.task.hpc.scheduler_state ?? selectedEntry.status}</span>
                      <button type="button" onClick={() => void openHpcActivityWindow()}>
                        <ExternalLink size={15} /> Popout
                      </button>
                    </div>
                    <div className="task-hpc-meta">
                      <span>Node<strong>{selectedEntry.task.hpc.remote_node ?? "Not assigned"}</strong></span>
                      <span>Remote directory<strong>{selectedEntry.task.hpc.remote_workdir ?? "Not available"}</strong></span>
                      <span>Local sync<strong>{selectedEntry.task.hpc.local_sync_dir ?? "Not available"}</strong></span>
                    </div>
                    <LiveOutputPanel
                      title="Live Log"
                      output={selectedEntry.task.outputText}
                      totalLineCount={selectedEntry.task.outputLineCount}
                      visibleLineCount={selectedEntry.task.output.length}
                      panelClassName="task-live-output"
                      outputClassName="task-live-output-text"
                      showTools={false}
                    />
                    {(selectedScript.submitCommand || selectedScript.script) && (
                      <details className="task-hpc-script">
                        <summary>Generated Slurm script</summary>
                        {selectedScript.submitCommand && <code>{selectedScript.submitCommand}</code>}
                        {selectedScript.script && <pre>{selectedScript.script}</pre>}
                      </details>
                    )}
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      </aside>
    </div>
  );
}
