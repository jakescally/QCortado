import assert from "node:assert/strict";
import test from "node:test";
import {
  buildTaskManagerEntries,
  filterTaskManagerEntries,
  findRelevantTaskManagerEntry,
  summarizeTaskManagerEntries,
} from "../../src/lib/taskManager";
import type { QueueItem, TaskState } from "../../src/lib/TaskContext";

function task(overrides: Partial<TaskState> & Pick<TaskState, "taskId" | "status">): TaskState {
  return {
    taskId: overrides.taskId,
    taskType: "scf",
    label: overrides.taskId,
    startedAt: "2026-06-08T12:00:00.000Z",
    status: overrides.status,
    progress: { status: "running", percent: 25, phase: "Running" },
    output: [],
    outputText: "",
    outputLineCount: 0,
    result: null,
    error: null,
    hpc: {
      backend: null,
      hpc_resource_type: null,
      remote_job_id: null,
      scheduler_state: null,
      remote_node: null,
      remote_workdir: null,
      remote_project_path: null,
      remote_storage_bytes: null,
      hpc_profile_id: null,
      local_sync_dir: null,
      recovery_save: null,
      headless_attached: false,
    },
    metadata: null,
    ...overrides,
  };
}

function queue(overrides: Partial<QueueItem> & Pick<QueueItem, "id" | "status">): QueueItem {
  return {
    id: overrides.id,
    taskType: "scf",
    label: overrides.id,
    params: {},
    status: overrides.status,
    createdAt: "2026-06-08T11:00:00.000Z",
    startedAt: null,
    finishedAt: null,
    taskId: null,
    error: null,
    saveSpec: null,
    ...overrides,
  };
}

test("merges a started queue item with its task and keeps direct tasks", () => {
  const linked = task({ taskId: "linked", status: "running", label: "Linked task" });
  const direct = task({
    taskId: "direct",
    status: "running",
    label: "Direct HPC task",
    hpc: { ...linked.hpc, backend: "hpc" },
  });
  const entries = buildTaskManagerEntries(
    [linked, direct],
    [
      queue({ id: "queued", status: "queued" }),
      queue({ id: "linked-queue", status: "running", taskId: "linked" }),
    ],
  );

  assert.equal(entries.length, 3);
  assert.equal(entries.filter((entry) => entry.taskId === "linked").length, 1);
  assert.equal(entries.find((entry) => entry.taskId === "linked")?.queueItemId, "linked-queue");
  assert.equal(entries.find((entry) => entry.taskId === "direct")?.isHpc, true);
});

test("orders running before queued before finished and reports summaries", () => {
  const entries = buildTaskManagerEntries(
    [
      task({ taskId: "finished", status: "completed" }),
      task({ taskId: "running", status: "running" }),
      task({ taskId: "failed", status: "failed" }),
    ],
    [queue({ id: "queued", status: "queued" })],
  );

  assert.deepEqual(entries.map((entry) => entry.group), ["running", "queued", "finished", "finished"]);
  assert.equal(findRelevantTaskManagerEntry(entries)?.taskId, "running");
  assert.deepEqual(summarizeTaskManagerEntries(entries), {
    running: 1,
    queued: 1,
    finished: 2,
    failed: 1,
    hpc: 0,
    hpcRunning: 0,
    hpcQueued: 0,
    total: 4,
  });
});

test("filters active, finished, and HPC entries", () => {
  const local = task({ taskId: "local", status: "completed" });
  const hpc = task({
    taskId: "hpc",
    status: "running",
    hpc: { ...local.hpc, backend: "hpc" },
  });
  const entries = buildTaskManagerEntries([local, hpc], [queue({ id: "queued", status: "queued" })]);

  assert.deepEqual(filterTaskManagerEntries(entries, "active").map((entry) => entry.id), ["task:hpc", "queue:queued"]);
  assert.deepEqual(filterTaskManagerEntries(entries, "finished").map((entry) => entry.id), ["task:local"]);
  assert.deepEqual(filterTaskManagerEntries(entries, "hpc").map((entry) => entry.id), ["task:hpc"]);
});

test("recognizes queued HPC work before a backend task exists", () => {
  const entries = buildTaskManagerEntries([], [
    queue({
      id: "queued-hpc",
      status: "queued",
      params: { executionTarget: { mode: "hpc" } },
    }),
  ]);

  assert.equal(entries[0]?.isHpc, true);
  assert.deepEqual(summarizeTaskManagerEntries(entries), {
    running: 0,
    queued: 1,
    finished: 0,
    failed: 0,
    hpc: 1,
    hpcRunning: 0,
    hpcQueued: 1,
    total: 1,
  });
});

test("preserves resumable utility metadata on task manager entries", () => {
  const utility = task({
    taskId: "init-task",
    taskType: "wien2k_scf_initialize",
    status: "completed",
    metadata: {
      wizardResume: {
        view: "wien2k-scf-wizard",
        utility: "initialization",
        context: { projectId: "project-1", cifId: "cif-1" },
      },
    },
  });
  const entries = buildTaskManagerEntries([utility], []);

  assert.equal(entries[0]?.taskType, "wien2k_scf_initialize");
  assert.equal(entries[0]?.metadata?.wizardResume?.view, "wien2k-scf-wizard");
});

test("preserves WIEN2k structure stage resume metadata", () => {
  const utility = task({
    taskId: "structure-task",
    taskType: "wien2k_structure_stage",
    status: "running",
    metadata: {
      wizardResume: {
        view: "wien2k-structure-wizard",
        utility: "structure_stage",
        stage: "rmt",
        context: { projectId: "project-1", cifId: "cif-1" },
      },
    },
  });
  const entries = buildTaskManagerEntries([utility], []);

  assert.equal(entries[0]?.taskType, "wien2k_structure_stage");
  assert.equal(entries[0]?.metadata?.wizardResume?.stage, "rmt");
});
