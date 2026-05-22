const SCF_RUN_SETTINGS_STORAGE_KEY = "qcortado.scf_run_settings_clipboard_text";
export const SCF_RUN_SETTINGS_UPDATED_EVENT = "qcortado:scf-run-settings-updated";
const SCF_RUN_SETTINGS_FORMAT = "qcortado-scf-settings-v1";

export interface ScfRunSettingsSource {
  calc_id: string;
  calc_type: string;
  calc_name: string | null;
  element_symbols: string[];
}

export interface ScfRunSettingsClipboardPayload {
  format: string;
  exported_at: string;
  source: ScfRunSettingsSource;
  settings: Record<string, unknown>;
}

interface ScfRunLike {
  id: string;
  calc_type: string;
  name?: string | null;
  parameters?: unknown;
}

function asRecord(value: unknown): Record<string, unknown> | null {
  if (!value || typeof value !== "object" || Array.isArray(value)) {
    return null;
  }
  return value as Record<string, unknown>;
}

function collectElementSymbolsFromSettings(settings: Record<string, unknown>): string[] {
  const sourceKeys: string[] = [];
  const selectedPseudos = asRecord(settings.selected_pseudos);
  if (selectedPseudos) {
    sourceKeys.push(...Object.keys(selectedPseudos));
  }
  const startingMagnetization = asRecord(settings.starting_magnetization);
  if (startingMagnetization) {
    sourceKeys.push(...Object.keys(startingMagnetization));
  }
  const startingMagnetizationTheta = asRecord(
    settings.starting_magnetization_theta ?? settings.starting_magnetization_angle1 ?? settings.theta ?? settings.angle1,
  );
  if (startingMagnetizationTheta) {
    sourceKeys.push(...Object.keys(startingMagnetizationTheta));
  }
  const startingMagnetizationPhi = asRecord(
    settings.starting_magnetization_phi ?? settings.starting_magnetization_angle2 ?? settings.phi ?? settings.angle2,
  );
  if (startingMagnetizationPhi) {
    sourceKeys.push(...Object.keys(startingMagnetizationPhi));
  }
  const hubbardU = asRecord(settings.hubbard_u);
  if (hubbardU) {
    sourceKeys.push(...Object.keys(hubbardU));
  }
  const hubbardManifold = asRecord(settings.hubbard_manifold);
  if (hubbardManifold) {
    sourceKeys.push(...Object.keys(hubbardManifold));
  }
  return [...new Set(sourceKeys.map((key) => key.trim()).filter(Boolean))].sort();
}

export function serializeScfRunSettings(run: ScfRunLike): string {
  const settings = asRecord(run.parameters) ?? {};
  const payload: ScfRunSettingsClipboardPayload = {
    format: SCF_RUN_SETTINGS_FORMAT,
    exported_at: new Date().toISOString(),
    source: {
      calc_id: run.id,
      calc_type: run.calc_type,
      calc_name: run.name?.trim() || null,
      element_symbols: collectElementSymbolsFromSettings(settings),
    },
    settings,
  };

  return [
    "QCortado SCF Run Settings",
    `Format: ${SCF_RUN_SETTINGS_FORMAT}`,
    "",
    JSON.stringify(payload, null, 2),
  ].join("\n");
}

function normalizeClipboardPayload(raw: unknown): ScfRunSettingsClipboardPayload | null {
  const parsed = asRecord(raw);
  if (!parsed) return null;
  if (parsed.format !== SCF_RUN_SETTINGS_FORMAT) return null;

  const sourceRecord = asRecord(parsed.source);
  const settingsRecord = asRecord(parsed.settings);
  if (!sourceRecord || !settingsRecord) return null;

  const calcId = String(sourceRecord.calc_id ?? "").trim();
  if (!calcId) return null;

  const elementSymbolsRaw = Array.isArray(sourceRecord.element_symbols)
    ? sourceRecord.element_symbols
    : [];
  const elementSymbols = [...new Set(
    elementSymbolsRaw
      .map((entry) => String(entry).trim())
      .filter(Boolean),
  )].sort();

  return {
    format: SCF_RUN_SETTINGS_FORMAT,
    exported_at: String(parsed.exported_at ?? ""),
    source: {
      calc_id: calcId,
      calc_type: String(sourceRecord.calc_type ?? ""),
      calc_name: sourceRecord.calc_name == null ? null : String(sourceRecord.calc_name),
      element_symbols: elementSymbols,
    },
    settings: settingsRecord,
  };
}

function parseJsonCandidate(candidate: string): ScfRunSettingsClipboardPayload | null {
  const trimmed = candidate.trim();
  if (!trimmed) return null;
  try {
    const parsed = JSON.parse(trimmed);
    return normalizeClipboardPayload(parsed);
  } catch {
    return null;
  }
}

export function parseScfRunSettingsClipboardText(text: string | null | undefined): ScfRunSettingsClipboardPayload | null {
  if (!text) return null;

  const candidates = new Set<string>();
  candidates.add(text);

  const jsonStart = text.indexOf("{");
  const jsonEnd = text.lastIndexOf("}");
  if (jsonStart >= 0 && jsonEnd > jsonStart) {
    candidates.add(text.slice(jsonStart, jsonEnd + 1));
  }

  for (const candidate of candidates) {
    const parsed = parseJsonCandidate(candidate);
    if (parsed) {
      return parsed;
    }
  }
  return null;
}

export function getStoredScfRunSettingsClipboardText(): string | null {
  try {
    return window.localStorage.getItem(SCF_RUN_SETTINGS_STORAGE_KEY);
  } catch {
    return null;
  }
}

export function rememberScfRunSettingsClipboardText(text: string): void {
  try {
    window.localStorage.setItem(SCF_RUN_SETTINGS_STORAGE_KEY, text);
  } catch {
    // Ignore localStorage failures; clipboard still carries the payload.
  }
  window.dispatchEvent(new Event(SCF_RUN_SETTINGS_UPDATED_EVENT));
}

export function hasStoredScfRunSettingsClipboardText(): boolean {
  return parseScfRunSettingsClipboardText(getStoredScfRunSettingsClipboardText()) !== null;
}
