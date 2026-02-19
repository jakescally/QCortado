import { useEffect, useState } from "react";
import { useTaskContext } from "../lib/TaskContext";
import { ElapsedTimer } from "./ElapsedTimer";

interface ProcessIndicatorProps {
  onNavigateToTask?: (taskId: string, taskType: string) => void;
}

export function ProcessIndicator({ onNavigateToTask }: ProcessIndicatorProps) {
  const { activeTasks, cancelTask, dismissTask, queueSummary } = useTaskContext();
  const [confirmingCancel, setConfirmingCancel] = useState<string | null>(null);
  const [isHidden, setIsHidden] = useState(false);

  // Show the most relevant task: prefer running, then most recent non-running
  const task = activeTasks.find((t) => t.status === "running") ??
    activeTasks[activeTasks.length - 1] ?? null;

  useEffect(() => {
    if (activeTasks.length === 0) {
      setIsHidden(false);
      setConfirmingCancel(null);
    }
  }, [activeTasks.length]);

  if (!task) return null;

  const isRunning = task.status === "running";
  const isComplete = task.status === "completed";
  const isFailed = task.status === "failed" || task.status === "cancelled";

  const statusClass = isRunning ? "running" : isComplete ? "completed" : "failed";
  const percent = task.progress.percent;
  const clamped = percent === null ? 0 : Math.max(0, Math.min(100, percent));

  function handleClick() {
    if (onNavigateToTask) {
      onNavigateToTask(task.taskId, task.taskType);
    }
  }

  async function handleCancel() {
    if (confirmingCancel === task.taskId) {
      await cancelTask(task.taskId);
      setConfirmingCancel(null);
    } else {
      setConfirmingCancel(task.taskId);
      // Auto-dismiss confirmation after 3 seconds
      setTimeout(() => setConfirmingCancel(null), 3000);
    }
  }

  async function handleDismiss() {
    await dismissTask(task.taskId);
  }

  function handleHide() {
    setIsHidden(true);
    setConfirmingCancel(null);
  }

  const typeLabels: Record<string, string> = {
    scf: "SCF",
    bands: "Band Structure",
    dos: "Electronic DOS",
    fermi_surface: "Fermi Surface",
    phonon: "Phonon",
  };
  const typeLabel = typeLabels[task.taskType] || task.taskType.toUpperCase();

  if (isHidden) {
    const summary = isRunning
      ? `${typeLabel} in progress`
      : isComplete
        ? `${typeLabel} complete`
        : `${typeLabel} ${task.status === "cancelled" ? "cancelled" : "failed"}`;
    return (
      <button
        type="button"
        className={`process-indicator-collapsed ${statusClass}`}
        onClick={() => setIsHidden(false)}
        title={`Show status: ${summary}`}
        aria-label={`Show status box: ${summary}`}
      >
        <svg
          className="process-indicator-collapsed-ring"
          viewBox="0 0 86 34"
          preserveAspectRatio="none"
          aria-hidden="true"
        >
          <rect
            className="process-indicator-collapsed-ring-head"
            x="1"
            y="1"
            width="84"
            height="32"
            rx="16"
            ry="16"
            pathLength="100"
          />
        </svg>
        <span className="process-indicator-collapsed-dot" aria-hidden="true" />
        <span className="process-indicator-collapsed-text">Status</span>
      </button>
    );
  }

  return (
    <div className={`process-indicator ${statusClass}`}>
      <div className="process-indicator-content" onClick={handleClick}>
        <div className="process-indicator-header">
          <span className="process-indicator-type">{typeLabel}</span>
          <span className="process-indicator-label">{task.label}</span>
          {queueSummary.total > 0 && queueSummary.activeIndex !== null && (
            <span className="process-indicator-queue">
              {queueSummary.activeIndex} of {queueSummary.total} in queue
            </span>
          )}
          {isRunning && (
            <ElapsedTimer
              startedAt={task.startedAt}
              isRunning={true}
              label=""
            />
          )}
        </div>

        {isRunning && (
          <div className="process-indicator-progress">
            <div className={`process-indicator-track ${percent === null ? "indeterminate" : ""}`}>
              <div
                className="process-indicator-fill"
                style={percent === null ? undefined : { width: `${clamped}%` }}
              />
            </div>
            <span className="process-indicator-phase">{task.progress.phase}</span>
            {task.progress.detail && (
              <span className="process-indicator-detail">{task.progress.detail}</span>
            )}
          </div>
        )}

        {isComplete && (
          <div className="process-indicator-status">Complete</div>
        )}

        {isFailed && (
          <div className="process-indicator-status process-indicator-error">
            {task.status === "cancelled" ? "Cancelled" : "Failed"}
          </div>
        )}
      </div>

      <div className="process-indicator-actions">
        <button
          type="button"
          className="process-indicator-btn process-indicator-hide"
          onClick={(e) => { e.stopPropagation(); handleHide(); }}
          title="Hide status box"
        >
          Hide
        </button>
        {isRunning && (
          <button
            type="button"
            className="process-indicator-btn process-indicator-cancel"
            onClick={(e) => { e.stopPropagation(); handleCancel(); }}
            title={confirmingCancel === task.taskId ? "Click again to confirm" : "Cancel calculation"}
          >
            {confirmingCancel === task.taskId ? "Confirm cancel" : "Cancel"}
          </button>
        )}
        {!isRunning && (
          <button
            type="button"
            className="process-indicator-btn process-indicator-dismiss"
            onClick={(e) => { e.stopPropagation(); handleDismiss(); }}
            title="Dismiss"
          >
            Dismiss
          </button>
        )}
      </div>
    </div>
  );
}
