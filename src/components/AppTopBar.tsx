import { Activity, Menu } from "lucide-react";
import type { TaskManagerEntry, TaskManagerSummary } from "../lib/taskManager";

interface AppTopBarProps {
  contextLabel: string;
  task: TaskManagerEntry | null;
  summary: TaskManagerSummary;
  onOpenMenu: () => void;
  onOpenTasks: (taskId?: string | null) => void;
  embeddedHeader?: boolean;
}

function statusLabel(task: TaskManagerEntry): string {
  if (task.status === "running" && task.task?.progress.phase) return task.task.progress.phase;
  if (task.status === "saving") return "Saving";
  if (task.status === "queued") return "Queued";
  if (task.status === "completed") return "Complete";
  if (task.status === "cancelled") return "Cancelled";
  return "Failed";
}

export function AppTopBar({
  contextLabel,
  task,
  summary,
  onOpenMenu,
  onOpenTasks,
  embeddedHeader = false,
}: AppTopBarProps) {
  return (
    <header className={`app-top-bar ${embeddedHeader ? "app-top-bar-shell" : ""}`}>
      <button type="button" className="chrome-icon-btn" onClick={onOpenMenu} aria-label="Open application menu">
        <Menu size={19} />
      </button>
      {embeddedHeader ? (
        <div
          id="app-dynamic-header-slot"
          className="app-dynamic-header-slot"
          aria-label={`${contextLabel} controls`}
        />
      ) : (
        <div className="app-top-bar-context">{contextLabel}</div>
      )}
      <button
        type="button"
        className={`task-status-pill ${embeddedHeader ? "task-status-pill-shell" : ""} ${task?.group ?? "idle"} ${task?.status ?? ""}`.trim()}
        onClick={() => onOpenTasks(task?.taskId)}
        aria-label="Open task manager"
      >
        <Activity size={16} />
        <span className="task-status-pill-copy">
          <strong>{task ? task.label : "No active tasks"}</strong>
          <span>{task ? statusLabel(task) : "Task manager"}</span>
        </span>
        {summary.running + summary.queued > 0 && (
          <span className="task-status-count">{summary.running + summary.queued}</span>
        )}
      </button>
    </header>
  );
}
