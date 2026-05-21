const PROJECT_WIZARD_SETTINGS_PREFIX = "qcortado.project_wizard_settings.v1";

function safeProjectId(projectId: string | null | undefined): string | null {
  const trimmed = String(projectId || "").trim();
  return trimmed.length > 0 ? trimmed : null;
}

function storageKey(projectId: string, wizardId: string): string {
  return `${PROJECT_WIZARD_SETTINGS_PREFIX}.${encodeURIComponent(projectId)}.${wizardId}`;
}

export function readProjectWizardSettings<T extends object>(
  projectId: string | null | undefined,
  wizardId: string,
): Partial<T> | null {
  if (typeof window === "undefined") return null;
  const resolvedProjectId = safeProjectId(projectId);
  if (!resolvedProjectId) return null;

  try {
    const raw = window.localStorage.getItem(storageKey(resolvedProjectId, wizardId));
    if (!raw) return null;
    const parsed = JSON.parse(raw);
    if (!parsed || typeof parsed !== "object" || Array.isArray(parsed)) return null;
    return parsed as Partial<T>;
  } catch {
    return null;
  }
}

export function writeProjectWizardSettings<T extends object>(
  projectId: string | null | undefined,
  wizardId: string,
  settings: T,
): void {
  if (typeof window === "undefined") return;
  const resolvedProjectId = safeProjectId(projectId);
  if (!resolvedProjectId) return;

  try {
    window.localStorage.setItem(storageKey(resolvedProjectId, wizardId), JSON.stringify(settings));
  } catch {
    // Ignore persistence failures; wizard defaults still work.
  }
}
