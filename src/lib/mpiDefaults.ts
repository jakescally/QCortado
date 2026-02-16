import { invoke } from "@tauri-apps/api/core";

export interface MpiDefaults {
  enabled: boolean;
  nprocs: number;
}

const MPI_DEFAULTS_STORAGE_KEY = "qcortado_global_mpi_defaults";

export function getDefaultMpiProcs(cpuCount: number): number {
  return Math.max(1, Math.floor(Math.max(1, cpuCount) * 0.75));
}

export function clampMpiProcs(rawValue: number, cpuCount: number): number {
  if (!Number.isFinite(rawValue)) {
    return getDefaultMpiProcs(cpuCount);
  }
  const safeCpuCount = Math.max(1, Math.floor(cpuCount));
  const rounded = Math.floor(rawValue);
  return Math.max(1, Math.min(safeCpuCount, rounded));
}

function normalizeMpiDefaults(raw: Partial<MpiDefaults> | null | undefined, cpuCount: number): MpiDefaults {
  return {
    enabled: Boolean(raw?.enabled),
    nprocs: clampMpiProcs(Number(raw?.nprocs ?? getDefaultMpiProcs(cpuCount)), cpuCount),
  };
}

function loadMpiDefaultsFromLocalStorage(cpuCount: number): MpiDefaults | null {
  try {
    const raw = localStorage.getItem(MPI_DEFAULTS_STORAGE_KEY);
    if (!raw) return null;
    const parsed = JSON.parse(raw) as Partial<MpiDefaults>;
    return normalizeMpiDefaults(parsed, cpuCount);
  } catch {
    return null;
  }
}

function saveMpiDefaultsToLocalStorage(defaults: MpiDefaults): void {
  try {
    localStorage.setItem(MPI_DEFAULTS_STORAGE_KEY, JSON.stringify(defaults));
  } catch {
    // Ignore local storage errors; backend persistence can still succeed.
  }
}

export async function loadGlobalMpiDefaults(cpuCount: number): Promise<MpiDefaults> {
  const localFallback = loadMpiDefaultsFromLocalStorage(cpuCount);
  try {
    const defaults = await invoke<MpiDefaults | null>("get_mpi_defaults");
    if (defaults) {
      const normalized = normalizeMpiDefaults(defaults, cpuCount);
      saveMpiDefaultsToLocalStorage(normalized);
      return normalized;
    }
    if (localFallback) {
      // Backend has no value yet; repopulate from local app-wide cache.
      try {
        await invoke("set_mpi_defaults", { defaults: localFallback });
      } catch {
        // Best-effort sync only.
      }
      return localFallback;
    }
  } catch {
    if (localFallback) {
      return localFallback;
    }
  }
  const normalizedDefault = normalizeMpiDefaults(null, cpuCount);
  saveMpiDefaultsToLocalStorage(normalizedDefault);
  return normalizedDefault;
}

export async function saveGlobalMpiDefaults(defaults: MpiDefaults, cpuCount: number): Promise<MpiDefaults> {
  const normalized = normalizeMpiDefaults(defaults, cpuCount);
  saveMpiDefaultsToLocalStorage(normalized);
  try {
    await invoke("set_mpi_defaults", { defaults: normalized });
  } catch {
    // Keep local persistence even if backend command/config is unavailable.
  }
  return normalized;
}
