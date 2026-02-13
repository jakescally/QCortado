import { useMemo } from "react";
import { QueueItem, useTaskContext } from "../lib/TaskContext";

interface TaskQueuePageProps {
  onBack: () => void;
}

function formatDate(iso: string | null): string {
  if (!iso) return "—";
  try {
    return new Date(iso).toLocaleString();
  } catch {
    return iso;
  }
}

function statusLabel(item: QueueItem): string {
  switch (item.status) {
    case "queued":
      return "Queued";
    case "running":
      return "Running";
    case "saving":
      return "Saving";
    case "completed":
      return "Complete";
    case "failed":
      return "Failed";
    case "cancelled":
      return "Cancelled";
    default:
      return item.status;
  }
}

export function TaskQueuePage({ onBack }: TaskQueuePageProps) {
  const {
    queueItems,
    moveQueueItem,
    removeQueueItem,
    cancelQueueItem,
    clearFinishedQueueItems,
    getTask,
  } = useTaskContext();

  const pendingCount = useMemo(
    () => queueItems.filter((item) => item.status === "queued" || item.status === "running" || item.status === "saving").length,
    [queueItems],
  );
  const completedCount = useMemo(
    () => queueItems.filter((item) => item.status === "completed").length,
    [queueItems],
  );
  const failedCount = useMemo(
    () => queueItems.filter((item) => item.status === "failed" || item.status === "cancelled").length,
    [queueItems],
  );
  const queuedOnly = useMemo(
    () => queueItems.filter((item) => item.status === "queued").map((item) => item.id),
    [queueItems],
  );

  return (
    <div className="queue-page-container">
      <div className="queue-page-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Task Queue</h2>
      </div>

      <div className="queue-page-summary">
        <div className="queue-summary-card">
          <span className="queue-summary-label">Pending</span>
          <span className="queue-summary-value">{pendingCount}</span>
        </div>
        <div className="queue-summary-card">
          <span className="queue-summary-label">Completed</span>
          <span className="queue-summary-value">{completedCount}</span>
        </div>
        <div className="queue-summary-card">
          <span className="queue-summary-label">Failed/Cancelled</span>
          <span className="queue-summary-value">{failedCount}</span>
        </div>
        <button
          className="queue-clear-finished-btn"
          onClick={clearFinishedQueueItems}
          disabled={completedCount + failedCount === 0}
        >
          Clear Finished
        </button>
      </div>

      {queueItems.length === 0 ? (
        <div className="queue-empty-state">
          <h3>No queued tasks</h3>
          <p>Queued calculations will appear here.</p>
        </div>
      ) : (
        <div className="queue-list">
          {queueItems.map((item) => {
            const task = item.taskId ? getTask(item.taskId) : undefined;
            const queuedIndex = queuedOnly.indexOf(item.id);
            const canMoveUp = item.status === "queued" && queuedIndex > 0;
            const canMoveDown = item.status === "queued" && queuedIndex >= 0 && queuedIndex < queuedOnly.length - 1;
            const canCancel = item.status === "queued" || item.status === "running" || item.status === "saving";
            const canRemove = item.status !== "running" && item.status !== "saving";

            return (
              <div key={item.id} className="queue-item">
                <div className="queue-item-main">
                  <div className="queue-item-title-row">
                    <span className="queue-item-type">{item.taskType.toUpperCase()}</span>
                    <span className={`queue-item-status status-${item.status}`}>
                      {statusLabel(item)}
                    </span>
                  </div>
                  <div className="queue-item-label">{item.label}</div>
                  <div className="queue-item-meta">
                    <span>Queued: {formatDate(item.createdAt)}</span>
                    <span>Started: {formatDate(item.startedAt)}</span>
                    <span>Finished: {formatDate(item.finishedAt)}</span>
                  </div>
                  {task?.progress?.phase && (item.status === "running" || item.status === "saving") && (
                    <div className="queue-item-progress">
                      <span>{task.progress.phase}</span>
                      {task.progress.percent != null && <span>{Math.round(task.progress.percent)}%</span>}
                    </div>
                  )}
                  {item.error && (
                    <div className="queue-item-error">{item.error}</div>
                  )}
                </div>

                <div className="queue-item-actions">
                  <button onClick={() => moveQueueItem(item.id, "up")} disabled={!canMoveUp}>
                    ↑
                  </button>
                  <button onClick={() => moveQueueItem(item.id, "down")} disabled={!canMoveDown}>
                    ↓
                  </button>
                  <button onClick={() => void cancelQueueItem(item.id)} disabled={!canCancel}>
                    Cancel
                  </button>
                  <button onClick={() => removeQueueItem(item.id)} disabled={!canRemove}>
                    Remove
                  </button>
                </div>
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}

