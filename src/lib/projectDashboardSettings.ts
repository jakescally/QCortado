const STORAGE_KEY = "qcortado.project-dashboard.engine-filter.v1";

export function getStoredProjectDashboardEngineFilter(): string {
  if (typeof window === "undefined") return "all";

  try {
    const stored = window.localStorage.getItem(STORAGE_KEY);
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
    window.localStorage.setItem(STORAGE_KEY, filter);
  } catch {
    // Ignore persistence failures.
  }
}
