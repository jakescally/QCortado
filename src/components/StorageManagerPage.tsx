import { useEffect, useMemo, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { ProgressBar } from "./ProgressBar";
import type {
  ExecutionMode,
  HpcProfile,
  StorageCalculationEntry,
  StorageCalculationTarget,
  StorageDeleteEntryTarget,
  StorageEntryKind,
  StorageInventoryMode,
  StorageInventoryProgress,
  StorageInventoryResult,
  StorageMutationResult,
  StoragePathEntry,
} from "../lib/types";

interface StorageManagerPageProps {
  onBack: () => void;
  executionMode: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
}

type SortMode = "size_desc" | "size_asc" | "label";
type ConfirmAction = "delete_entries" | "delete_selected" | "lighten_calculations";

interface GroupConfig {
  key: "projects" | "calculations" | "temp" | "orphans" | "other";
  title: string;
  description: string;
}

const GROUPS: GroupConfig[] = [
  { key: "calculations", title: "Calculations", description: "Saved calculation artifacts and outputs." },
  { key: "projects", title: "Project Data", description: "Project metadata, structures, and non-calculation files." },
  { key: "temp", title: "Temporary Storage", description: "QCortado temp folders under system temp roots." },
  { key: "orphans", title: "Orphans", description: "Remote QCortado directories not linked to saved calculations." },
  { key: "other", title: "Other", description: "Additional QCortado-owned metadata or support directories." },
];

function createProgressEventId(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") {
    return crypto.randomUUID();
  }
  return `storage-${Date.now()}-${Math.random().toString(36).slice(2)}`;
}

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIdx = 0;
  while (value >= 1024 && unitIdx < units.length - 1) {
    value /= 1024;
    unitIdx += 1;
  }
  const precision = value >= 10 || unitIdx === 0 ? 0 : 1;
  return `${value.toFixed(precision)} ${units[unitIdx]}`;
}

function formatDateTime(value: string | null | undefined): string {
  if (!value) return "Never";
  try {
    return new Date(value).toLocaleString();
  } catch {
    return value;
  }
}

function sortByMode<T extends { bytes: number; label: string }>(items: T[], sortMode: SortMode): T[] {
  const next = [...items];
  if (sortMode === "size_asc") {
    next.sort((a, b) => a.bytes - b.bytes || a.label.localeCompare(b.label));
    return next;
  }
  if (sortMode === "label") {
    next.sort((a, b) => a.label.localeCompare(b.label) || b.bytes - a.bytes);
    return next;
  }
  next.sort((a, b) => b.bytes - a.bytes || a.label.localeCompare(b.label));
  return next;
}

function progressPercent(progress: StorageInventoryProgress | null): number | null {
  if (!progress) return null;
  if (progress.total_bytes_estimate > 0) {
    return Math.min(100, Math.max(0, Math.round((progress.bytes_seen / progress.total_bytes_estimate) * 100)));
  }
  if (progress.total_items > 0) {
    return Math.min(100, Math.max(0, Math.round((progress.scanned_items / progress.total_items) * 100)));
  }
  return null;
}

function progressPhaseLabel(progress: StorageInventoryProgress | null): string {
  if (!progress) return "Preparing scan";
  switch (progress.phase) {
    case "preparing":
      return "Preparing local scan";
    case "scanning-calculations":
      return "Scanning calculations";
    case "scanning-projects":
      return "Scanning project data";
    case "scanning-temp":
      return "Scanning temp storage";
    case "scanning-other":
      return "Scanning other QCortado data";
    case "discovering-remote":
      return "Discovering remote paths";
    case "scanning-remote":
      return "Scanning remote storage";
    case "deleting-storage":
      return "Deleting storage";
    case "deleting-calculations":
      return "Deleting calculations";
    case "deleting-selected":
      return "Deleting selected items";
    case "lightening-calculations":
      return "Lightening calculations";
    case "done":
      return "Operation complete";
    default:
      return progress.phase;
  }
}

function progressDetail(progress: StorageInventoryProgress | null): string | undefined {
  if (!progress) return undefined;
  if (progress.total_items > 0) {
    return `${progress.scanned_items}/${progress.total_items} items • ${formatBytes(progress.bytes_seen)} scanned`;
  }
  if (progress.bytes_seen > 0) {
    return `${formatBytes(progress.bytes_seen)} scanned`;
  }
  return undefined;
}

function entryKindLabel(kind: StorageEntryKind): string {
  switch (kind) {
    case "project":
      return "Project";
    case "temp":
      return "Temp";
    case "orphan":
      return "Orphan";
    case "other":
      return "Other";
    default:
      return kind;
  }
}

export function StorageManagerPage({
  onBack,
  executionMode,
  activeHpcProfile = null,
}: StorageManagerPageProps) {
  const [mode, setMode] = useState<StorageInventoryMode>(executionMode === "hpc" ? "hpc" : "local");
  const [inventory, setInventory] = useState<StorageInventoryResult | null>(null);
  const [progress, setProgress] = useState<StorageInventoryProgress | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [isMutating, setIsMutating] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [statusMessage, setStatusMessage] = useState<string | null>(null);
  const [sortMode, setSortMode] = useState<SortMode>("size_desc");
  const [selectedIds, setSelectedIds] = useState<string[]>([]);
  const [confirmAction, setConfirmAction] = useState<ConfirmAction | null>(null);
  const [scanNonce, setScanNonce] = useState(0);
  const activeProgressIdRef = useRef<string | null>(null);

  useEffect(() => {
    setMode(executionMode === "hpc" ? "hpc" : "local");
  }, [executionMode]);

  useEffect(() => {
    const unlistenPromise = listen<StorageInventoryProgress>("storage-manager-progress", (event) => {
      const payload = event.payload;
      if (!activeProgressIdRef.current || payload.progress_event_id !== activeProgressIdRef.current) {
        return;
      }
      setProgress(payload);
    });

    return () => {
      activeProgressIdRef.current = null;
      void unlistenPromise.then((unlisten) => unlisten());
    };
  }, []);

  useEffect(() => {
    let cancelled = false;
    const progressEventId = createProgressEventId();
    activeProgressIdRef.current = progressEventId;

    async function loadInventory() {
      if (mode === "hpc" && !activeHpcProfile) {
        setInventory(null);
        setProgress(null);
        setIsLoading(false);
        setError("No active HPC profile is configured.");
        return;
      }

      setIsLoading(true);
      setError(null);
      setStatusMessage(null);
      setProgress({
        progress_event_id: progressEventId,
        phase: "preparing",
        scanned_items: 0,
        total_items: 0,
        bytes_seen: 0,
        total_bytes_estimate: 0,
      });

      try {
        const result = await invoke<StorageInventoryResult>("scan_storage_inventory", {
          mode,
          profileId: mode === "hpc" ? activeHpcProfile?.id ?? null : null,
          progressEventId,
        });
        if (cancelled || activeProgressIdRef.current !== progressEventId) return;
        setInventory(result);
        setSelectedIds([]);
      } catch (scanError) {
        if (cancelled || activeProgressIdRef.current !== progressEventId) return;
        console.error("Failed to scan storage inventory:", scanError);
        setInventory(null);
        setError(String(scanError));
      } finally {
        if (!cancelled && activeProgressIdRef.current === progressEventId) {
          setIsLoading(false);
        }
      }
    }

    void loadInventory();
    return () => {
      cancelled = true;
    };
  }, [mode, activeHpcProfile?.id, scanNonce]);

  const calculationEntries = useMemo(
    () => sortByMode(inventory?.calculations ?? [], sortMode),
    [inventory?.calculations, sortMode],
  );
  const projectEntries = useMemo(
    () => sortByMode(inventory?.projects ?? [], sortMode),
    [inventory?.projects, sortMode],
  );
  const tempEntries = useMemo(
    () => sortByMode(inventory?.temp ?? [], sortMode),
    [inventory?.temp, sortMode],
  );
  const orphanEntries = useMemo(
    () => sortByMode(inventory?.orphans ?? [], sortMode),
    [inventory?.orphans, sortMode],
  );
  const otherEntries = useMemo(
    () => sortByMode(inventory?.other ?? [], sortMode),
    [inventory?.other, sortMode],
  );

  const pathEntriesById = useMemo(() => {
    const map = new Map<string, StoragePathEntry>();
    for (const entry of [...projectEntries, ...tempEntries, ...orphanEntries, ...otherEntries]) {
      map.set(entry.id, entry);
    }
    return map;
  }, [projectEntries, tempEntries, orphanEntries, otherEntries]);

  const calculationsById = useMemo(() => {
    const map = new Map<string, StorageCalculationEntry>();
    for (const calc of calculationEntries) {
      map.set(calc.id, calc);
    }
    return map;
  }, [calculationEntries]);

  const selectedCalculations = useMemo(
    () => selectedIds
      .map((id) => calculationsById.get(id))
      .filter((entry): entry is StorageCalculationEntry => Boolean(entry)),
    [calculationsById, selectedIds],
  );
  const selectedPathEntries = useMemo(
    () => selectedIds
      .map((id) => pathEntriesById.get(id))
      .filter((entry): entry is StoragePathEntry => Boolean(entry)),
    [pathEntriesById, selectedIds],
  );

  const selectedTotalBytes = useMemo(
    () => selectedCalculations.reduce((sum, entry) => sum + entry.bytes, 0)
      + selectedPathEntries.reduce((sum, entry) => sum + entry.bytes, 0),
    [selectedCalculations, selectedPathEntries],
  );
  const selectedContainsMixedTypes = selectedCalculations.length > 0 && selectedPathEntries.length > 0;
  const canDeleteSelected = selectedIds.length > 0
    && selectedCalculations.every((entry) => entry.can_delete)
    && selectedPathEntries.every((entry) => entry.delete_supported);
  const canLightenCalculations = selectedCalculations.length > 0
    && selectedPathEntries.length === 0
    && selectedCalculations.every((entry) => entry.can_lighten)
    && mode === "local";
  const canDeleteEntries = selectedPathEntries.length > 0
    && selectedCalculations.length === 0
    && selectedPathEntries.every((entry) => entry.delete_supported);

  function queueRescan(message?: string) {
    if (message) {
      setStatusMessage(message);
    }
    setScanNonce((value) => value + 1);
  }

  function toggleSelection(id: string) {
    setSelectedIds((current) => (
      current.includes(id) ? current.filter((entryId) => entryId !== id) : [...current, id]
    ));
  }

  function toggleGroupSelection(ids: string[]) {
    if (ids.length === 0) return;
    setSelectedIds((current) => {
      const allSelected = ids.every((id) => current.includes(id));
      if (allSelected) {
        return current.filter((id) => !ids.includes(id));
      }
      const merged = new Set(current);
      for (const id of ids) merged.add(id);
      return [...merged];
    });
  }

  function areAllSelected(ids: string[]): boolean {
    return ids.length > 0 && ids.every((id) => selectedIds.includes(id));
  }

  async function runDeleteSelected() {
    const calculations: StorageCalculationTarget[] = selectedCalculations.map((entry) => ({
      project_id: entry.project_id,
      cif_id: entry.cif_id,
      calc_id: entry.calc_id,
    }));
    const entries: StorageDeleteEntryTarget[] = selectedPathEntries.map((entry) => ({
      id: entry.id,
      kind: entry.kind,
      path: entry.path,
    }));
    const progressEventId = createProgressEventId();
    activeProgressIdRef.current = progressEventId;
    setProgress({
      progress_event_id: progressEventId,
      phase: "deleting-selected",
      scanned_items: 0,
      total_items: calculations.length + entries.length,
      bytes_seen: 0,
      total_bytes_estimate: 0,
    });
    setIsMutating(true);
    setError(null);
    try {
      const result = await invoke<StorageMutationResult>("delete_storage_selection", {
        request: {
          mode,
          profileId: mode === "hpc" ? activeHpcProfile?.id ?? null : null,
          progressEventId,
          calculations,
          entries,
        },
      });
      if (result.failures.length > 0) {
        setError(result.failures.map((failure) => failure.message).join("\n"));
      }
      queueRescan(
        `Deleted ${result.succeeded.length} item${result.succeeded.length === 1 ? "" : "s"} and freed ${formatBytes(result.bytes_freed)}.`,
      );
    } catch (mutationError) {
      console.error("Failed to delete selected storage:", mutationError);
      setError(String(mutationError));
    } finally {
      setIsMutating(false);
      setConfirmAction(null);
    }
  }

  async function runDeleteEntries() {
    const entries: StorageDeleteEntryTarget[] = selectedPathEntries.map((entry) => ({
      id: entry.id,
      kind: entry.kind,
      path: entry.path,
    }));
    const progressEventId = createProgressEventId();
    activeProgressIdRef.current = progressEventId;
    setProgress({
      progress_event_id: progressEventId,
      phase: "deleting-storage",
      scanned_items: 0,
      total_items: entries.length,
      bytes_seen: 0,
      total_bytes_estimate: 0,
    });
    setIsMutating(true);
    setError(null);
    try {
      const result = await invoke<StorageMutationResult>("delete_storage_entries", {
        request: {
          mode,
          profileId: mode === "hpc" ? activeHpcProfile?.id ?? null : null,
          progressEventId,
          entries,
        },
      });
      if (result.failures.length > 0) {
        setError(result.failures.map((failure) => failure.message).join("\n"));
      }
      queueRescan(`Deleted ${result.succeeded.length} storage entr${result.succeeded.length === 1 ? "y" : "ies"} and freed ${formatBytes(result.bytes_freed)}.`);
    } catch (mutationError) {
      console.error("Failed to delete storage entries:", mutationError);
      setError(String(mutationError));
    } finally {
      setIsMutating(false);
      setConfirmAction(null);
    }
  }

  async function runLightenCalculations() {
    const calculations: StorageCalculationTarget[] = selectedCalculations.map((entry) => ({
      project_id: entry.project_id,
      cif_id: entry.cif_id,
      calc_id: entry.calc_id,
    }));
    const progressEventId = createProgressEventId();
    activeProgressIdRef.current = progressEventId;
    setProgress({
      progress_event_id: progressEventId,
      phase: "lightening-calculations",
      scanned_items: 0,
      total_items: calculations.length,
      bytes_seen: 0,
      total_bytes_estimate: 0,
    });
    setIsMutating(true);
    setError(null);
    try {
      const result = await invoke<StorageMutationResult>("lighten_storage_calculations", {
        request: { calculations, progressEventId },
      });
      if (result.failures.length > 0) {
        setError(result.failures.map((failure) => failure.message).join("\n"));
      }
      queueRescan(`Lightened ${result.succeeded.length} calculation${result.succeeded.length === 1 ? "" : "s"} and freed ${formatBytes(result.bytes_freed)}.`);
    } catch (mutationError) {
      console.error("Failed to lighten calculations:", mutationError);
      setError(String(mutationError));
    } finally {
      setIsMutating(false);
      setConfirmAction(null);
    }
  }

  function renderPathGroup(group: GroupConfig, entries: StoragePathEntry[]) {
    if (entries.length === 0) return null;
    const selectableEntries = entries.filter((entry) => entry.delete_supported);
    const allSelected = areAllSelected(selectableEntries.map((entry) => entry.id));
    return (
      <section key={group.key} className="storage-section-card">
        <div className="storage-group-header">
          <div>
            <h3>{group.title}</h3>
            <p className="settings-menu-hint">{group.description}</p>
          </div>
          <div className="storage-group-actions">
            <span className="storage-group-total">{formatBytes(entries.reduce((sum, entry) => sum + entry.bytes, 0))}</span>
            {selectableEntries.length > 0 && (
              <button className="settings-page-tab" type="button" onClick={() => toggleGroupSelection(selectableEntries.map((entry) => entry.id))}>
                {allSelected ? "Deselect All" : "Select All"}
              </button>
            )}
          </div>
        </div>
        <div className="storage-entry-list">
          {entries.map((entry) => (
            <label
              key={entry.id}
              className={`queue-item storage-entry-row ${selectedIds.includes(entry.id) ? "selected" : ""} ${!entry.delete_supported ? "readonly" : ""}`}
            >
              <input
                type="checkbox"
                checked={selectedIds.includes(entry.id)}
                disabled={!entry.delete_supported}
                onChange={() => toggleSelection(entry.id)}
              />
              <div className="storage-entry-main">
                <div className="storage-entry-title-row">
                  <span className="storage-entry-title">{entry.label}</span>
                  <div className="calc-tags">
                    <span className="calc-tag calc-tag-info">{entryKindLabel(entry.kind)}</span>
                  </div>
                </div>
                {entry.description && <div className="queue-item-meta storage-entry-description"><span>{entry.description}</span></div>}
              </div>
              <div className="storage-entry-size">{formatBytes(entry.bytes)}</div>
            </label>
          ))}
        </div>
      </section>
    );
  }

  function renderCalculationGroup() {
    if (calculationEntries.length === 0) return null;
    const selectableCalculationIds = calculationEntries
      .filter((entry) => entry.can_delete || entry.can_lighten)
      .map((entry) => entry.id);
    const allSelected = areAllSelected(selectableCalculationIds);
    return (
      <section className="storage-section-card">
        <div className="storage-group-header">
          <div>
            <h3>Calculations</h3>
            <p className="settings-menu-hint">Saved calculation directories across all projects.</p>
          </div>
          <div className="storage-group-actions">
            <span className="storage-group-total">{formatBytes(calculationEntries.reduce((sum, entry) => sum + entry.bytes, 0))}</span>
            <button className="settings-page-tab" type="button" onClick={() => toggleGroupSelection(selectableCalculationIds)}>
              {allSelected ? "Deselect All" : "Select All"}
            </button>
          </div>
        </div>
        <div className="storage-entry-list">
          {calculationEntries.map((entry) => (
            <label
              key={entry.id}
              className={`queue-item storage-entry-row ${selectedIds.includes(entry.id) ? "selected" : ""}`}
            >
              <input
                type="checkbox"
                checked={selectedIds.includes(entry.id)}
                onChange={() => toggleSelection(entry.id)}
              />
              <div className="storage-entry-main">
                <div className="storage-entry-title-row">
                  <span className="storage-entry-title">{entry.label}</span>
                  <div className="calc-tags storage-entry-badges">
                    <span className="calc-tag calc-tag-info">{entry.calc_type.toUpperCase()}</span>
                    {entry.has_remote_artifacts && <span className="calc-tag calc-tag-hpc">Remote</span>}
                    {entry.can_lighten && <span className="calc-tag calc-tag-geometry">Lightenable</span>}
                  </div>
                </div>
                <div className="queue-item-label storage-entry-subtitle">{entry.project_name} • {entry.cif_label}</div>
                {(entry.local_path || entry.remote_paths.length > 0) && (
                  <div className="queue-item-meta storage-entry-meta">
                    {entry.local_path && <span>Local: {entry.local_path}</span>}
                    {entry.remote_paths.length > 0 && <span>Remote: {entry.remote_paths.length} linked path{entry.remote_paths.length === 1 ? "" : "s"}</span>}
                  </div>
                )}
              </div>
              <div className="storage-entry-size">{formatBytes(entry.bytes)}</div>
            </label>
          ))}
        </div>
      </section>
    );
  }

  const confirmModal = confirmAction ? (
    <div className="storage-confirm-overlay" onClick={() => !isMutating && setConfirmAction(null)}>
      <div className="storage-confirm-dialog" onClick={(event) => event.stopPropagation()}>
        <h3>
          {confirmAction === "delete_entries" && "Delete selected storage?"}
          {confirmAction === "delete_selected" && "Delete selected items?"}
          {confirmAction === "lighten_calculations" && "Lighten selected calculations?"}
        </h3>
        {confirmAction === "delete_entries" && (
          <p>
            This will permanently remove {selectedPathEntries.length} storage entr{selectedPathEntries.length === 1 ? "y" : "ies"}.
          </p>
        )}
        {confirmAction === "delete_selected" && (
          <>
            <p>
              This will permanently delete {[
                selectedCalculations.length > 0 ? `${selectedCalculations.length} calculation${selectedCalculations.length === 1 ? "" : "s"}` : null,
                selectedPathEntries.length > 0 ? `${selectedPathEntries.length} storage entr${selectedPathEntries.length === 1 ? "y" : "ies"}` : null,
              ].filter(Boolean).join(" and ")}.
            </p>
            {selectedCalculations.some((entry) => entry.has_remote_artifacts) && (
              <p className="storage-confirm-warning">
                One or more selected calculations also has linked remote HPC artifacts. Those remote files will be deleted too.
              </p>
            )}
          </>
        )}
        {confirmAction === "lighten_calculations" && (
          <p>
            This removes heavy QE wavefunction archives from {selectedCalculations.length} saved calculation{selectedCalculations.length === 1 ? "" : "s"} while keeping the rest of the saved results.
          </p>
        )}
        <div className="storage-confirm-actions">
          <button type="button" className="dialog-btn cancel" onClick={() => setConfirmAction(null)} disabled={isMutating}>
            Cancel
          </button>
          {confirmAction === "delete_entries" && (
            <button type="button" className="dialog-btn delete" onClick={() => void runDeleteEntries()} disabled={isMutating}>
              {isMutating ? "Deleting..." : "Delete Storage"}
            </button>
          )}
          {confirmAction === "delete_selected" && (
            <button type="button" className="dialog-btn delete" onClick={() => void runDeleteSelected()} disabled={isMutating}>
              {isMutating ? "Deleting..." : "Delete Selected"}
            </button>
          )}
          {confirmAction === "lighten_calculations" && (
            <button type="button" className="dialog-btn save" onClick={() => void runLightenCalculations()} disabled={isMutating}>
              {isMutating ? "Lightening..." : "Lighten Calculations"}
            </button>
          )}
        </div>
      </div>
    </div>
  ) : null;

  return (
    <div className="queue-page-container storage-manager-page">
      <div className="queue-page-header storage-manager-header">
        <div className="storage-manager-header-main">
          <button className="back-button" onClick={onBack}>
            ← Back
          </button>
          <h2>Storage Manager</h2>
        </div>
        <div className="storage-manager-toolbar">
          <div className="theme-toggle-group" role="group" aria-label="Storage mode">
            <button
              type="button"
              className={`theme-toggle-btn ${mode === "local" ? "active" : ""}`}
              onClick={() => setMode("local")}
            >
              Local
            </button>
            <button
              type="button"
              className={`theme-toggle-btn ${mode === "hpc" ? "active" : ""}`}
              onClick={() => setMode("hpc")}
            >
              HPC
            </button>
          </div>
          <select className="settings-menu-input storage-sort-select" value={sortMode} onChange={(event) => setSortMode(event.target.value as SortMode)}>
            <option value="size_desc">Biggest First</option>
            <option value="size_asc">Smallest First</option>
            <option value="label">Label</option>
          </select>
          <button className="dialog-btn cancel" type="button" onClick={() => queueRescan()} disabled={isLoading || isMutating}>
            {isLoading ? "Scanning..." : "Refresh"}
          </button>
        </div>
      </div>

      <div className="storage-manager-subhead">
        <span>Showing {mode === "hpc" ? "remote" : "local"} QCortado storage.</span>
        {mode === "hpc" && (
          <span>Profile: {activeHpcProfile ? `${activeHpcProfile.name} (${activeHpcProfile.username}@${activeHpcProfile.host})` : "None"}</span>
        )}
        <span>Last scan: {formatDateTime(inventory?.scanned_at)}</span>
      </div>

      {(isLoading || isMutating || progress) && (
        <div className="storage-progress-card">
          <ProgressBar
            status={error ? "error" : progress?.phase === "done" ? "complete" : "running"}
            percent={progressPercent(progress)}
            phase={progressPhaseLabel(progress)}
            detail={progressDetail(progress)}
          />
        </div>
      )}

      {error && <div className="storage-status error">{error}</div>}
      {statusMessage && <div className="storage-status success">{statusMessage}</div>}

      <div className="queue-page-summary storage-summary-grid">
        <div className="queue-summary-card">
          <span className="queue-summary-label">Total</span>
          <span className="queue-summary-value">{formatBytes(inventory?.total_bytes ?? 0)}</span>
        </div>
        <div className="queue-summary-card">
          <span className="queue-summary-label">Calculations</span>
          <span className="queue-summary-value">{formatBytes(inventory?.calculation_bytes ?? 0)}</span>
        </div>
        <div className="queue-summary-card">
          <span className="queue-summary-label">{mode === "local" ? "Temporary" : "Orphans"}</span>
          <span className="queue-summary-value">{formatBytes(mode === "local" ? (inventory?.temp_bytes ?? 0) : (inventory?.orphan_bytes ?? 0))}</span>
        </div>
        <div className="queue-summary-card">
          <span className="queue-summary-label">Other</span>
          <span className="queue-summary-value">{formatBytes((inventory?.project_bytes ?? 0) + (inventory?.other_bytes ?? 0))}</span>
        </div>
      </div>

      <div className="storage-manager-layout">
        <div className="storage-manager-main">
          {renderCalculationGroup()}
          {renderPathGroup(GROUPS[1], projectEntries)}
          {renderPathGroup(GROUPS[2], tempEntries)}
          {renderPathGroup(GROUPS[3], orphanEntries)}
          {renderPathGroup(GROUPS[4], otherEntries)}
          {!isLoading && inventory && calculationEntries.length === 0 && projectEntries.length === 0 && tempEntries.length === 0 && orphanEntries.length === 0 && otherEntries.length === 0 && (
            <div className="queue-empty-state">
              <h3>No QCortado storage found</h3>
              <p>Open a project or run a calculation to populate this view.</p>
            </div>
          )}
        </div>

        <aside className="storage-manager-sidebar">
          <div className="storage-section-card storage-sidebar-card">
            <h3>Selection</h3>
            {selectedIds.length === 0 ? (
              <p className="settings-menu-hint">Select calculations or deletable storage entries to manage them.</p>
            ) : (
              <>
                <div className="storage-selection-metrics">
                  <span>{selectedIds.length} selected</span>
                  <strong>{formatBytes(selectedTotalBytes)}</strong>
                </div>
                {selectedContainsMixedTypes && (
                  <p className="storage-selection-warning">
                    Mixed selections can still be deleted together. Lightening only applies to saved calculations.
                  </p>
                )}
                <div className="storage-selection-list">
                  {selectedCalculations.map((entry) => (
                    <div key={entry.id} className="storage-selection-item">
                      <span>{entry.label}</span>
                      <small>{entry.project_name}</small>
                    </div>
                  ))}
                  {selectedPathEntries.map((entry) => (
                    <div key={entry.id} className="storage-selection-item">
                      <span>{entry.label}</span>
                      <small>{entry.description ?? entry.path}</small>
                    </div>
                  ))}
                </div>
                <div className="storage-sidebar-actions">
                  <button className="settings-menu-item danger" type="button" onClick={() => setConfirmAction("delete_selected")} disabled={!canDeleteSelected || isMutating}>
                    Delete Selected
                  </button>
                  <button className="settings-menu-item" type="button" onClick={() => setConfirmAction("lighten_calculations")} disabled={!canLightenCalculations || isMutating}>
                    Lighten Calculations
                  </button>
                  <button className="settings-menu-item warning" type="button" onClick={() => setConfirmAction("delete_entries")} disabled={!canDeleteEntries || isMutating}>
                    Delete Temp/Orphans
                  </button>
                  <button className="settings-menu-item" type="button" onClick={() => setSelectedIds([])} disabled={isMutating}>
                    Clear Selection
                  </button>
                </div>
              </>
            )}
          </div>
        </aside>
      </div>

      {confirmModal}
    </div>
  );
}
