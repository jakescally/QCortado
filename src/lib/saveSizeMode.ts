import { invoke } from "@tauri-apps/api/core";

export type SaveSizeMode = "large" | "small";

function normalizeSaveSizeMode(raw: unknown): SaveSizeMode {
  return raw === "small" ? "small" : "large";
}

export async function loadGlobalSaveSizeMode(): Promise<SaveSizeMode> {
  try {
    const mode = await invoke<SaveSizeMode>("get_save_size_mode");
    return normalizeSaveSizeMode(mode);
  } catch {
    return "large";
  }
}

export async function saveGlobalSaveSizeMode(mode: SaveSizeMode): Promise<SaveSizeMode> {
  const normalized = normalizeSaveSizeMode(mode);
  try {
    await invoke("set_save_size_mode", { mode: normalized });
  } catch {
    // Keep UI state even if backend persistence is unavailable.
  }
  return normalized;
}
