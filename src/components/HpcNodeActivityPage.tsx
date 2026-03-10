import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  getHpcClusterSnapshot,
  HpcClusterSnapshot,
  HpcNodeSnapshot,
  HpcQueueScope,
  HpcQueueSnapshot,
} from "../lib/hpcConfig";
import { ExecutionMode } from "../lib/types";

interface HpcNodeActivityPageProps {
  onBack: () => void;
  executionMode: ExecutionMode;
  activeProfileId: string | null;
  activeProfileName?: string | null;
}

type ActivityTab = "nodes" | "queue";

const AUTO_REFRESH_MS = 15_000;
const NODE_ACTIVITY_FILTERS_STORAGE_KEY = "qcortado.hpc.node_activity.filters.v1";

interface PersistedNodeActivityFilters {
  tab: ActivityTab;
  queueScope: HpcQueueScope;
  nodeSearch: string;
  nodeTypeFilter: "all" | "cpu" | "gpu";
  nodePartitionFilter: string;
  nodeStateFilter: string;
  nodeMinCpu: string;
  nodeMinGpu: string;
  nodeMinMemoryGb: string;
  queueSearch: string;
  queueStateFilter: string;
}

function loadPersistedNodeActivityFilters(): Partial<PersistedNodeActivityFilters> {
  try {
    const raw = window.localStorage.getItem(NODE_ACTIVITY_FILTERS_STORAGE_KEY);
    if (!raw) return {};
    const parsed = JSON.parse(raw) as Partial<PersistedNodeActivityFilters>;
    return parsed && typeof parsed === "object" ? parsed : {};
  } catch {
    return {};
  }
}

function persistNodeActivityFilters(value: PersistedNodeActivityFilters): void {
  try {
    window.localStorage.setItem(NODE_ACTIVITY_FILTERS_STORAGE_KEY, JSON.stringify(value));
  } catch {
    // Ignore storage errors (private mode, quota, etc).
  }
}

function parseNonNegativeInt(value: string): number {
  const parsed = Number.parseInt(value, 10);
  if (!Number.isFinite(parsed) || parsed < 0) return 0;
  return parsed;
}

function formatGbFromMb(valueMb?: number | null): string {
  if (!Number.isFinite(valueMb) || valueMb === null || valueMb === undefined || valueMb <= 0) {
    return "0";
  }
  const gb = valueMb / 1024;
  return gb >= 100 ? gb.toFixed(0) : gb.toFixed(1);
}

function formatTimestamp(timestamp?: string | null): string {
  if (!timestamp) return "Not yet queried";
  const date = new Date(timestamp);
  if (Number.isNaN(date.getTime())) return timestamp;
  return date.toLocaleTimeString();
}

function nodeStateTone(state: string): "good" | "busy" | "warn" | "bad" | "neutral" {
  const normalized = state.trim().toUpperCase();
  if (normalized === "IDLE") return "good";
  if (normalized === "ALLOC" || normalized === "MIXED" || normalized === "RUNNING") return "busy";
  if (normalized.includes("DRAIN") || normalized.includes("MAINT")) return "warn";
  if (normalized.includes("DOWN") || normalized.includes("FAIL") || normalized.includes("UNKNOWN")) return "bad";
  return "neutral";
}

function queueStateTone(state: string): "good" | "busy" | "warn" | "bad" | "neutral" {
  const normalized = state.trim().toUpperCase();
  if (normalized === "RUNNING" || normalized === "COMPLETING") return "busy";
  if (normalized === "PENDING" || normalized === "CONFIGURING") return "warn";
  if (normalized === "COMPLETED") return "good";
  if (normalized.includes("FAIL") || normalized.includes("CANCEL") || normalized === "TIMEOUT") return "bad";
  return "neutral";
}

function normalizeText(value: string): string {
  return value.trim().toLocaleLowerCase();
}

function matchesNodeSearch(node: HpcNodeSnapshot, query: string): boolean {
  if (!query) return true;
  const haystack = [
    node.node_name,
    node.partitions.join(","),
    node.state,
    node.raw_state,
    node.cpu_name || "",
    node.gpu_name || "",
    node.features.join(","),
  ]
    .join(" ")
    .toLocaleLowerCase();
  return haystack.includes(query);
}

function matchesQueueSearch(job: HpcQueueSnapshot, query: string): boolean {
  if (!query) return true;
  const haystack = [
    job.job_id,
    job.user,
    job.state,
    job.raw_state,
    job.partition,
    job.nodelist,
    job.reason,
    job.name,
  ]
    .join(" ")
    .toLocaleLowerCase();
  return haystack.includes(query);
}

export function HpcNodeActivityPage({
  onBack,
  executionMode,
  activeProfileId,
  activeProfileName,
}: HpcNodeActivityPageProps) {
  const persistedFiltersRef = useRef<Partial<PersistedNodeActivityFilters>>(
    loadPersistedNodeActivityFilters(),
  );

  const [tab, setTab] = useState<ActivityTab>(
    persistedFiltersRef.current.tab === "queue" ? "queue" : "nodes",
  );
  const [snapshot, setSnapshot] = useState<HpcClusterSnapshot | null>(null);
  const [isInitialLoading, setIsInitialLoading] = useState(true);
  const [isRefreshing, setIsRefreshing] = useState(false);
  const [refreshError, setRefreshError] = useState<string | null>(null);
  const [queueScope, setQueueScope] = useState<HpcQueueScope>(
    persistedFiltersRef.current.queueScope === "mine" ? "mine" : "all",
  );
  const [refreshCycle, setRefreshCycle] = useState(0);

  const [nodeSearch, setNodeSearch] = useState(persistedFiltersRef.current.nodeSearch || "");
  const [nodeTypeFilter, setNodeTypeFilter] = useState<"all" | "cpu" | "gpu">(
    persistedFiltersRef.current.nodeTypeFilter === "cpu"
      ? "cpu"
      : persistedFiltersRef.current.nodeTypeFilter === "gpu"
        ? "gpu"
        : "all",
  );
  const [nodePartitionFilter, setNodePartitionFilter] = useState<string>(
    persistedFiltersRef.current.nodePartitionFilter || "all",
  );
  const [nodeStateFilter, setNodeStateFilter] = useState<string>(
    persistedFiltersRef.current.nodeStateFilter || "all",
  );
  const [nodeMinCpu, setNodeMinCpu] = useState(persistedFiltersRef.current.nodeMinCpu || "0");
  const [nodeMinGpu, setNodeMinGpu] = useState(persistedFiltersRef.current.nodeMinGpu || "0");
  const [nodeMinMemoryGb, setNodeMinMemoryGb] = useState(
    persistedFiltersRef.current.nodeMinMemoryGb || "0",
  );

  const [queueSearch, setQueueSearch] = useState(persistedFiltersRef.current.queueSearch || "");
  const [queueStateFilter, setQueueStateFilter] = useState<string>(
    persistedFiltersRef.current.queueStateFilter || "all",
  );

  const inFlightRef = useRef(false);
  const snapshotRef = useRef<HpcClusterSnapshot | null>(null);
  const timerRef = useRef<number | null>(null);

  useEffect(() => {
    snapshotRef.current = snapshot;
  }, [snapshot]);

  const clearRefreshTimer = useCallback(() => {
    if (timerRef.current !== null) {
      window.clearTimeout(timerRef.current);
      timerRef.current = null;
    }
  }, []);

  const resetAutoRefresh = useCallback(() => {
    setRefreshCycle((prev) => prev + 1);
  }, []);

  const fetchSnapshot = useCallback(async () => {
    if (executionMode !== "hpc") return;
    if (!activeProfileId) {
      setIsInitialLoading(false);
      setRefreshError("No active HPC profile selected. Open Settings > HPC and choose a profile.");
      return;
    }
    if (inFlightRef.current) return;

    inFlightRef.current = true;
    if (!snapshotRef.current) {
      setIsInitialLoading(true);
    }
    setIsRefreshing(true);
    try {
      const next = await getHpcClusterSnapshot(activeProfileId, queueScope, true, 1500);
      setSnapshot(next);
      setRefreshError(null);
    } catch (error) {
      setRefreshError(String(error));
    } finally {
      setIsInitialLoading(false);
      setIsRefreshing(false);
      inFlightRef.current = false;
    }
  }, [activeProfileId, executionMode, queueScope]);

  useEffect(() => {
    if (executionMode !== "hpc") return;
    void fetchSnapshot().finally(() => {
      resetAutoRefresh();
    });
  }, [executionMode, fetchSnapshot, resetAutoRefresh]);

  useEffect(() => {
    clearRefreshTimer();
    if (executionMode !== "hpc") return;
    timerRef.current = window.setTimeout(() => {
      void fetchSnapshot().finally(() => {
        setRefreshCycle((prev) => prev + 1);
      });
    }, AUTO_REFRESH_MS);
    return () => {
      clearRefreshTimer();
    };
  }, [clearRefreshTimer, executionMode, fetchSnapshot, refreshCycle]);

  useEffect(() => {
    persistNodeActivityFilters({
      tab,
      queueScope,
      nodeSearch,
      nodeTypeFilter,
      nodePartitionFilter,
      nodeStateFilter,
      nodeMinCpu,
      nodeMinGpu,
      nodeMinMemoryGb,
      queueSearch,
      queueStateFilter,
    });
  }, [
    nodeMinCpu,
    nodeMinGpu,
    nodeMinMemoryGb,
    nodePartitionFilter,
    nodeSearch,
    nodeStateFilter,
    nodeTypeFilter,
    queueScope,
    queueSearch,
    queueStateFilter,
    tab,
  ]);

  const partitionOptions = useMemo(() => {
    if (!snapshot) return [];
    const set = new Set<string>();
    for (const node of snapshot.nodes) {
      for (const partition of node.partitions) {
        if (partition.trim()) set.add(partition.trim());
      }
    }
    return Array.from(set).sort((a, b) => a.localeCompare(b));
  }, [snapshot]);

  const nodeStateOptions = useMemo(() => {
    if (!snapshot) return [];
    const set = new Set<string>();
    for (const node of snapshot.nodes) {
      if (node.state.trim()) set.add(node.state.trim());
    }
    return Array.from(set).sort((a, b) => a.localeCompare(b));
  }, [snapshot]);

  const queueStateOptions = useMemo(() => {
    if (!snapshot) return [];
    const set = new Set<string>();
    for (const job of snapshot.queue) {
      if (job.state.trim()) set.add(job.state.trim());
    }
    return Array.from(set).sort((a, b) => a.localeCompare(b));
  }, [snapshot]);

  const filteredNodes = useMemo(() => {
    if (!snapshot) return [];

    const query = normalizeText(nodeSearch);
    const minCpu = parseNonNegativeInt(nodeMinCpu);
    const minGpu = parseNonNegativeInt(nodeMinGpu);
    const minMemoryGb = parseNonNegativeInt(nodeMinMemoryGb);

    return snapshot.nodes.filter((node) => {
      if (nodeTypeFilter !== "all" && node.node_type.toLowerCase() !== nodeTypeFilter) {
        return false;
      }
      if (nodePartitionFilter !== "all" && !node.partitions.includes(nodePartitionFilter)) {
        return false;
      }
      if (nodeStateFilter !== "all" && node.state !== nodeStateFilter) {
        return false;
      }
      if (node.cpus < minCpu) {
        return false;
      }
      if (node.gpus < minGpu) {
        return false;
      }
      const memoryGb = node.memory_mb / 1024;
      if (memoryGb < minMemoryGb) {
        return false;
      }
      if (!matchesNodeSearch(node, query)) {
        return false;
      }
      return true;
    });
  }, [
    nodeMinCpu,
    nodeMinGpu,
    nodeMinMemoryGb,
    nodePartitionFilter,
    nodeSearch,
    nodeStateFilter,
    nodeTypeFilter,
    snapshot,
  ]);

  const filteredQueue = useMemo(() => {
    if (!snapshot) return [];
    const query = normalizeText(queueSearch);
    return snapshot.queue.filter((job) => {
      if (queueStateFilter !== "all" && job.state !== queueStateFilter) {
        return false;
      }
      if (!matchesQueueSearch(job, query)) {
        return false;
      }
      return true;
    });
  }, [queueSearch, queueStateFilter, snapshot]);

  const nodeSummary = useMemo(() => {
    const nodes = snapshot?.nodes || [];
    let idle = 0;
    let busy = 0;
    let unavailable = 0;
    for (const node of nodes) {
      const tone = nodeStateTone(node.state);
      if (tone === "good") idle += 1;
      else if (tone === "busy") busy += 1;
      else if (tone === "warn" || tone === "bad") unavailable += 1;
    }
    return {
      total: nodes.length,
      idle,
      busy,
      unavailable,
    };
  }, [snapshot?.nodes]);

  const queueSummary = useMemo(() => {
    const jobs = snapshot?.queue || [];
    let running = 0;
    let pending = 0;
    for (const job of jobs) {
      const state = job.state.toUpperCase();
      if (state === "RUNNING" || state === "COMPLETING") running += 1;
      if (state === "PENDING") pending += 1;
    }
    return {
      total: jobs.length,
      running,
      pending,
    };
  }, [snapshot?.queue]);

  const statusMessage = useMemo(() => {
    if (!snapshot) return "No snapshot loaded";
    const stamp = formatTimestamp(snapshot.captured_at);
    return `Updated ${stamp}`;
  }, [snapshot]);

  const showLoadingBlock = isInitialLoading && !snapshot;
  const showErrorBanner = refreshError !== null;

  return (
    <div className="queue-page-container node-activity-page">
      <div className="queue-page-header node-activity-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <div className="node-activity-heading">
          <h2>Node Activity</h2>
          <p>
            Cluster: {activeProfileName || "Andromeda"}
          </p>
        </div>
        <div className="node-activity-header-actions">
          <span className="node-activity-status">{statusMessage}</span>
          <button
            className="queue-clear-finished-btn"
            type="button"
            onClick={() => {
              void fetchSnapshot().finally(() => {
                resetAutoRefresh();
              });
            }}
            disabled={isRefreshing || !activeProfileId}
          >
            {isRefreshing ? "Refreshing..." : "Refresh"}
          </button>
        </div>
      </div>

      {!activeProfileId && (
        <div className="queue-empty-state">
          <h3>No active HPC profile</h3>
          <p>Select an Andromeda profile in Settings to query node activity.</p>
        </div>
      )}

      {showLoadingBlock && (
        <div className="node-activity-loading-state">
          <p className="pseudo-querying">
            <span className="pseudo-querying-spinner" />
            Querying the database...
          </p>
        </div>
      )}

      {showErrorBanner && (
        <div className="node-activity-error-banner" role="alert">
          <strong>Refresh failed.</strong> {refreshError}
          {snapshot && " Showing last successful snapshot."}
        </div>
      )}

      {snapshot && snapshot.warnings.length > 0 && (
        <div className="node-activity-warning-banner">
          {snapshot.warnings.map((warning) => (
            <p key={warning}>{warning}</p>
          ))}
        </div>
      )}

      {snapshot && (
        <>
          <div className="node-activity-tab-row">
            <button
              type="button"
              className={`node-activity-tab ${tab === "nodes" ? "active" : ""}`}
              onClick={() => setTab("nodes")}
            >
              Nodes
            </button>
            <button
              type="button"
              className={`node-activity-tab ${tab === "queue" ? "active" : ""}`}
              onClick={() => setTab("queue")}
            >
              Queue
            </button>
          </div>

          {tab === "nodes" && (
            <div className="node-activity-panel">
              <div className="queue-page-summary node-activity-summary">
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Total Nodes</span>
                  <span className="queue-summary-value">{nodeSummary.total}</span>
                </div>
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Idle</span>
                  <span className="queue-summary-value">{nodeSummary.idle}</span>
                </div>
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Allocated/Mixed</span>
                  <span className="queue-summary-value">{nodeSummary.busy}</span>
                </div>
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Unavailable</span>
                  <span className="queue-summary-value">{nodeSummary.unavailable}</span>
                </div>
              </div>

              <div className="node-activity-filter-grid">
                <label>
                  Search
                  <input
                    value={nodeSearch}
                    onChange={(event) => setNodeSearch(event.target.value)}
                    placeholder="node, partition, feature, reason"
                  />
                </label>
                <label>
                  Node Type
                  <select
                    value={nodeTypeFilter}
                    onChange={(event) =>
                      setNodeTypeFilter(
                        event.target.value === "cpu" || event.target.value === "gpu"
                          ? event.target.value
                          : "all",
                      )
                    }
                  >
                    <option value="all">All</option>
                    <option value="cpu">CPU</option>
                    <option value="gpu">GPU</option>
                  </select>
                </label>
                <label>
                  Partition
                  <select
                    value={nodePartitionFilter}
                    onChange={(event) => setNodePartitionFilter(event.target.value)}
                  >
                    <option value="all">All</option>
                    {partitionOptions.map((partition) => (
                      <option key={partition} value={partition}>
                        {partition}
                      </option>
                    ))}
                  </select>
                </label>
                <label>
                  State
                  <select
                    value={nodeStateFilter}
                    onChange={(event) => setNodeStateFilter(event.target.value)}
                  >
                    <option value="all">All</option>
                    {nodeStateOptions.map((state) => (
                      <option key={state} value={state}>
                        {state}
                      </option>
                    ))}
                  </select>
                </label>
                <label>
                  Min CPU
                  <input
                    type="number"
                    min={0}
                    value={nodeMinCpu}
                    onChange={(event) => setNodeMinCpu(event.target.value)}
                  />
                </label>
                <label>
                  Min GPU
                  <input
                    type="number"
                    min={0}
                    value={nodeMinGpu}
                    onChange={(event) => setNodeMinGpu(event.target.value)}
                  />
                </label>
                <label>
                  Min Memory (GB)
                  <input
                    type="number"
                    min={0}
                    value={nodeMinMemoryGb}
                    onChange={(event) => setNodeMinMemoryGb(event.target.value)}
                  />
                </label>
              </div>

              <div className="node-activity-table-shell">
                <table className="node-activity-table">
                  <thead>
                    <tr>
                      <th>Node</th>
                      <th>Type</th>
                      <th>State</th>
                      <th>Partition(s)</th>
                      <th>CPU Cores</th>
                      <th>CPU Name</th>
                      <th>GPU Count</th>
                      <th>GPU Name</th>
                      <th>Mem (GB)</th>
                      <th>Free (GB)</th>
                      <th>Features</th>
                    </tr>
                  </thead>
                  <tbody>
                    {filteredNodes.length === 0 ? (
                      <tr>
                        <td colSpan={11} className="node-activity-empty-cell">
                          No nodes matched current filters.
                        </td>
                      </tr>
                    ) : (
                      filteredNodes.map((node) => (
                        <tr key={node.node_name}>
                          <td>{node.node_name}</td>
                          <td>{node.node_type.toUpperCase()}</td>
                          <td>
                            <span className={`node-activity-pill tone-${nodeStateTone(node.state)}`}>
                              {node.state}
                            </span>
                          </td>
                          <td>{node.partitions.join(", ") || "—"}</td>
                          <td>{node.cpus}</td>
                          <td>{node.cpu_name || "—"}</td>
                          <td>{node.gpus}</td>
                          <td>{node.gpu_name || "—"}</td>
                          <td>{formatGbFromMb(node.memory_mb)}</td>
                          <td>{formatGbFromMb(node.free_memory_mb)}</td>
                          <td>{node.features.join(", ") || "—"}</td>
                        </tr>
                      ))
                    )}
                  </tbody>
                </table>
              </div>
            </div>
          )}

          {tab === "queue" && (
            <div className="node-activity-panel">
              <div className="queue-page-summary node-activity-summary">
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Total Jobs</span>
                  <span className="queue-summary-value">{queueSummary.total}</span>
                </div>
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Running</span>
                  <span className="queue-summary-value">{queueSummary.running}</span>
                </div>
                <div className="queue-summary-card">
                  <span className="queue-summary-label">Pending</span>
                  <span className="queue-summary-value">{queueSummary.pending}</span>
                </div>
              </div>

              <div className="node-activity-scope-toggle" role="group" aria-label="Queue scope">
                <button
                  type="button"
                  className={queueScope === "all" ? "active" : ""}
                  onClick={() => setQueueScope("all")}
                >
                  All Jobs
                </button>
                <button
                  type="button"
                  className={queueScope === "mine" ? "active" : ""}
                  onClick={() => setQueueScope("mine")}
                >
                  My Jobs
                </button>
              </div>

              <div className="node-activity-filter-grid">
                <label>
                  Search
                  <input
                    value={queueSearch}
                    onChange={(event) => setQueueSearch(event.target.value)}
                    placeholder="job, user, name, node, reason"
                  />
                </label>
                <label>
                  State
                  <select
                    value={queueStateFilter}
                    onChange={(event) => setQueueStateFilter(event.target.value)}
                  >
                    <option value="all">All</option>
                    {queueStateOptions.map((state) => (
                      <option key={state} value={state}>
                        {state}
                      </option>
                    ))}
                  </select>
                </label>
              </div>

              <div className="node-activity-table-shell">
                <table className="node-activity-table">
                  <thead>
                    <tr>
                      <th>Job</th>
                      <th>User</th>
                      <th>State</th>
                      <th>Partition</th>
                      <th>Nodes</th>
                      <th>Elapsed</th>
                      <th>Limit</th>
                      <th>Node / Reason</th>
                      <th>Name</th>
                    </tr>
                  </thead>
                  <tbody>
                    {filteredQueue.length === 0 ? (
                      <tr>
                        <td colSpan={9} className="node-activity-empty-cell">
                          No queue entries matched current filters.
                        </td>
                      </tr>
                    ) : (
                      filteredQueue.map((job) => (
                        <tr key={`${job.job_id}-${job.name}-${job.partition}`}>
                          <td>{job.job_id}</td>
                          <td>{job.user}</td>
                          <td>
                            <span className={`node-activity-pill tone-${queueStateTone(job.state)}`}>
                              {job.state}
                            </span>
                          </td>
                          <td>{job.partition || "—"}</td>
                          <td>{job.nodes}</td>
                          <td>{job.elapsed || "—"}</td>
                          <td>{job.time_limit || "—"}</td>
                          <td>{job.nodelist && job.nodelist !== "unknown" ? job.nodelist : job.reason || "—"}</td>
                          <td>{job.name || "—"}</td>
                        </tr>
                      ))
                    )}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
}
