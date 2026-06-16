import type { QueueItem, QueueItemStatus, TaskState, TaskStatus } from "./TaskContext";

export type TaskManagerFilter = "active" | "all" | "finished" | "hpc";
export type TaskManagerGroup = "running" | "queued" | "finished";

export interface TaskManagerEntry {
  id: string;
  taskId: string | null;
  queueItemId: string | null;
  taskType: string;
  label: string;
  status: QueueItemStatus | TaskStatus;
  group: TaskManagerGroup;
  createdAt: string;
  startedAt: string | null;
  finishedAt: string | null;
  error: string | null;
  isHpc: boolean;
  metadata: Record<string, any> | null;
  task: TaskState | null;
  queueItem: QueueItem | null;
}

export interface TaskManagerSummary {
  running: number;
  queued: number;
  finished: number;
  failed: number;
  hpc: number;
  hpcRunning: number;
  hpcQueued: number;
  total: number;
}

function groupForStatus(status: QueueItemStatus | TaskStatus): TaskManagerGroup {
  if (status === "running" || status === "saving") return "running";
  if (status === "queued") return "queued";
  return "finished";
}

function entryTimestamp(entry: TaskManagerEntry): number {
  const timestamp = entry.finishedAt ?? entry.startedAt ?? entry.createdAt;
  const parsed = Date.parse(timestamp);
  return Number.isFinite(parsed) ? parsed : 0;
}

function groupRank(group: TaskManagerGroup): number {
  if (group === "running") return 0;
  if (group === "queued") return 1;
  return 2;
}

export function buildTaskManagerEntries(
  tasks: Iterable<TaskState>,
  queueItems: QueueItem[],
): TaskManagerEntry[] {
  const tasksById = new Map(Array.from(tasks, (task) => [task.taskId, task]));
  const linkedTaskIds = new Set<string>();
  const entries: TaskManagerEntry[] = queueItems.map((queueItem) => {
    const task = queueItem.taskId ? tasksById.get(queueItem.taskId) ?? null : null;
    if (task) linkedTaskIds.add(task.taskId);
    const queueIsActive = queueItem.status === "queued"
      || queueItem.status === "running"
      || queueItem.status === "saving";
    const status = queueIsActive ? queueItem.status : task?.status ?? queueItem.status;

    return {
      id: `queue:${queueItem.id}`,
      taskId: task?.taskId ?? queueItem.taskId,
      queueItemId: queueItem.id,
      taskType: task?.taskType ?? queueItem.taskType,
      label: task?.label ?? queueItem.label,
      status,
      group: groupForStatus(status),
      createdAt: queueItem.createdAt,
      startedAt: task?.startedAt ?? queueItem.startedAt,
      finishedAt: queueItem.finishedAt,
      error: task?.error ?? queueItem.error,
      isHpc: task?.hpc.backend === "hpc"
        || String(queueItem.params?.executionTarget?.mode ?? "").toLowerCase() === "hpc",
      metadata: task?.metadata ?? null,
      task,
      queueItem,
    };
  });

  for (const task of tasksById.values()) {
    if (linkedTaskIds.has(task.taskId)) continue;
    entries.push({
      id: `task:${task.taskId}`,
      taskId: task.taskId,
      queueItemId: null,
      taskType: task.taskType,
      label: task.label,
      status: task.status,
      group: groupForStatus(task.status),
      createdAt: task.startedAt,
      startedAt: task.startedAt,
      finishedAt: null,
      error: task.error,
      isHpc: task.hpc.backend === "hpc",
      metadata: task.metadata,
      task,
      queueItem: null,
    });
  }

  return entries.sort((left, right) => {
    const groupDelta = groupRank(left.group) - groupRank(right.group);
    if (groupDelta !== 0) return groupDelta;
    if (left.group === "queued" && right.group === "queued") {
      return queueItems.findIndex((item) => item.id === left.queueItemId)
        - queueItems.findIndex((item) => item.id === right.queueItemId);
    }
    return entryTimestamp(right) - entryTimestamp(left);
  });
}

export function filterTaskManagerEntries(
  entries: TaskManagerEntry[],
  filter: TaskManagerFilter,
): TaskManagerEntry[] {
  if (filter === "active") return entries.filter((entry) => entry.group !== "finished");
  if (filter === "finished") return entries.filter((entry) => entry.group === "finished");
  if (filter === "hpc") return entries.filter((entry) => entry.isHpc);
  return entries;
}

export function summarizeTaskManagerEntries(entries: TaskManagerEntry[]): TaskManagerSummary {
  return entries.reduce<TaskManagerSummary>((summary, entry) => {
    summary.total += 1;
    summary[entry.group] += 1;
    if (entry.status === "failed" || entry.status === "cancelled") summary.failed += 1;
    if (entry.isHpc) {
      summary.hpc += 1;
      if (entry.group === "running") summary.hpcRunning += 1;
      if (entry.group === "queued") summary.hpcQueued += 1;
    }
    return summary;
  }, {
    running: 0,
    queued: 0,
    finished: 0,
    failed: 0,
    hpc: 0,
    hpcRunning: 0,
    hpcQueued: 0,
    total: 0,
  });
}

export function findRelevantTaskManagerEntry(entries: TaskManagerEntry[]): TaskManagerEntry | null {
  return entries.find((entry) => entry.group === "running")
    ?? entries.find((entry) => entry.group === "queued")
    ?? entries.find((entry) => entry.group === "finished")
    ?? null;
}
