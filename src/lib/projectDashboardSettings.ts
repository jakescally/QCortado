import type { EngineId } from "./engines/types";

const ENGINE_FILTER_STORAGE_KEY = "qcortado.project-dashboard.engine-filter.v1";
const ACTIVE_ENGINE_STORAGE_KEY = "qcortado.project-dashboard.active-engine.v1";

export function getStoredProjectDashboardEngineFilter(): string {
  if (typeof window === "undefined") return "all";

  try {
    const stored = window.localStorage.getItem(ENGINE_FILTER_STORAGE_KEY);
    if (typeof stored === "string" && stored.trim().length > 0) {
      return stored;
    }
  } catch {
    // Ignore persistence failures.
  }
  return "all";
}

export function setStoredProjectDashboardEngineFilter(
  filter: string,
): void {
  if (typeof window === "undefined") return;

  try {
    window.localStorage.setItem(ENGINE_FILTER_STORAGE_KEY, filter);
  } catch {
    // Ignore persistence failures.
  }
}

export function getStoredProjectDashboardActiveEngineId(): EngineId | null {
  if (typeof window === "undefined") return null;

  try {
    const stored = window.localStorage.getItem(ACTIVE_ENGINE_STORAGE_KEY);
    if (stored === "qe" || stored === "wien2k") {
      return stored;
    }
  } catch {
    // Ignore persistence failures.
  }
  return null;
}

export function setStoredProjectDashboardActiveEngineId(engineId: EngineId): void {
  if (typeof window === "undefined") return;

  try {
    window.localStorage.setItem(ACTIVE_ENGINE_STORAGE_KEY, engineId);
  } catch {
    // Ignore persistence failures.
  }
}
