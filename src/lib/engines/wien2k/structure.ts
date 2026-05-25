import { invoke } from "@tauri-apps/api/core";

export type Wien2kStructureStage = "rmt" | "sgroup" | "symmetry";
export type Wien2kStructureSessionPhase = "staged" | "rmt_ready" | "sgroup_ready" | "symmetry_ready";

export interface Wien2kStructureSiteOverride {
  siteIndex: number;
  npt?: number | null;
  r0?: number | null;
  rmt?: number | null;
}

export interface Wien2kStructureControls {
  rmtReductionPercent: number;
  sgroupTolerance: number;
  siteOverrides: Wien2kStructureSiteOverride[];
}

export interface Wien2kStructureSite {
  siteIndex: number;
  symbol: string;
  atomicNumber: number;
  positions: Array<[number, number, number]>;
  npt: number;
  r0: number;
  rmt: number;
}

export interface Wien2kStructureDraft {
  projectId: string;
  cifId: string;
  caseName: string;
  latticeType: string;
  spacegroupNumber: number;
  internationalSymbol: string;
  cellParameters: [number, number, number, number, number, number];
  standardizedLattice: number[][];
  sites: Wien2kStructureSite[];
  controls: Wien2kStructureControls;
  structContent: string;
}

export interface Wien2kStructureSession {
  sessionId: string;
  draft: Wien2kStructureDraft;
  remoteCaseDir: string;
  remoteInstallRoot: string;
  hpcProfileId: string;
  phase: Wien2kStructureSessionPhase;
  currentStruct: string;
  artifacts: Record<string, string>;
  transcript: string;
  startedAt: string;
}

export interface Wien2kStructureStageResult {
  sessionId: string;
  stage: Wien2kStructureStage;
  phase: Wien2kStructureSessionPhase;
  candidateStruct: string;
  sites: Wien2kStructureSite[];
  nativeOutput: string;
  saveAllowed: boolean;
  diagnostics: string[];
}

export interface Wien2kStructureSourceRecord {
  id: string;
  engine_id?: string | null;
  calc_type: string;
  parameters?: Record<string, unknown> | null;
  completed_at?: string | null;
}

export const DEFAULT_WIEN2K_STRUCTURE_CONTROLS: Wien2kStructureControls = {
  rmtReductionPercent: 0,
  sgroupTolerance: 1e-5,
  siteOverrides: [],
};

export function validateWien2kStructureControls(controls: Wien2kStructureControls): string | null {
  if (!Number.isFinite(controls.rmtReductionPercent)
    || controls.rmtReductionPercent < 0
    || controls.rmtReductionPercent >= 100) {
    return "RMT reduction must be at least 0 and less than 100 percent.";
  }
  if (!Number.isFinite(controls.sgroupTolerance)
    || controls.sgroupTolerance < 1e-7
    || controls.sgroupTolerance > 1e-3) {
    return "SGROUP tolerance must be between 1e-7 and 1e-3.";
  }
  for (const override of controls.siteOverrides) {
    if (override.npt != null && (!Number.isInteger(override.npt) || override.npt <= 0 || override.npt % 2 !== 1)) {
      return `NPT for site ${override.siteIndex + 1} must be a positive odd integer.`;
    }
    if (override.r0 != null && (!Number.isFinite(override.r0) || override.r0 <= 0)) {
      return `R0 for site ${override.siteIndex + 1} must be positive.`;
    }
    if (override.rmt != null && (!Number.isFinite(override.rmt) || override.rmt <= 0)) {
      return `RMT for site ${override.siteIndex + 1} must be positive.`;
    }
    if (override.r0 != null && override.rmt != null && override.r0 >= override.rmt) {
      return `R0 for site ${override.siteIndex + 1} must be smaller than RMT.`;
    }
  }
  return null;
}

export function listWien2kStructureSources<T extends Wien2kStructureSourceRecord>(calculations: T[]): T[] {
  return calculations.filter((calculation) =>
    calculation.engine_id === "wien2k"
      && calculation.calc_type === "engine_setup"
      && calculation.parameters?.setup_kind === "structure",
  );
}

export function prepareWien2kStructureDraft(
  projectId: string,
  cifId: string,
  caseName: string,
  controls: Wien2kStructureControls,
): Promise<Wien2kStructureDraft> {
  return invoke("wien2k_prepare_structure_draft", { projectId, cifId, caseName, controls });
}

export function startWien2kStructureSession(draft: Wien2kStructureDraft): Promise<Wien2kStructureSession> {
  return invoke("wien2k_start_structure_session", { draft });
}

export function runWien2kStructureStage(
  sessionId: string,
  stage: Wien2kStructureStage,
  controls: Wien2kStructureControls,
): Promise<Wien2kStructureStageResult> {
  return invoke("wien2k_run_structure_stage", { sessionId, stage, controls });
}

export function saveWien2kStructureSource(sessionId: string): Promise<Wien2kStructureSourceRecord> {
  return invoke("wien2k_save_structure_source", { sessionId });
}

export function discardWien2kStructureSession(sessionId: string): Promise<void> {
  return invoke("wien2k_discard_structure_session", { sessionId });
}
