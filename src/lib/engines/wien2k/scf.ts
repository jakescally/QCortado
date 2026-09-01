import { invoke } from "@tauri-apps/api/core";
import type { SlurmResourceRequest } from "../../types";
import type { Wien2kInitializationSettings, Wien2kScfRunSettings } from "./types";
import type { Wien2kStructureSite } from "./structure";

export type Wien2kScfSessionPhase = "staged" | "initialized" | "scf_complete" | "failed";
export type ScfConvergenceState = "converged" | "not_converged" | "failed" | "cancelled" | "unknown";

export interface NormalizedScfSummary {
  schema: "cortado.scf_summary.v1";
  convergence: ScfConvergenceState;
  totalEnergy?: { value: number; unit: "Ry" | "eV" | "Ha" } | null;
  fermiEnergyEv?: number | null;
  scfSteps?: number | null;
  wallTimeSeconds?: number | null;
  totalMagnetization?: number | null;
  metadata?: Record<string, unknown> | null;
}

export interface Wien2kScfSession {
  sessionId: string;
  projectId: string;
  cifId: string;
  sourceStructureCalculationId: string;
  sourceStructureSites: Wien2kStructureSite[];
  caseName: string;
  remoteCaseDir: string;
  remoteInstallRoot: string;
  hpcProfileId: string;
  phase: Wien2kScfSessionPhase;
  initialization?: Wien2kInitializationSettings | null;
  latestRun?: Wien2kScfRunSettings | null;
  latestCalculationId?: string | null;
  artifacts: Record<string, string>;
  transcript: string[];
  startedAt: string;
}

export interface Wien2kLstartCoreLeakSuggestion {
  suggestedCutoffRy: number;
  referenceEnergyRy: number;
  bufferRy: number;
  atom?: string | null;
  orbital?: string | null;
  leakCharge?: number | null;
}

export interface Wien2kInitializationResult {
  sessionId: string;
  phase: Wien2kScfSessionPhase;
  nativeOutput: string;
  diagnostics: string[];
  lstartCoreLeakSuggestion?: Wien2kLstartCoreLeakSuggestion | null;
  artifacts: Record<string, string>;
}

export interface Wien2kScfExecutionResult {
  sessionId: string;
  phase: Wien2kScfSessionPhase;
  nativeOutput: string;
  diagnostics: string[];
  summary: NormalizedScfSummary;
  calculationId: string;
  continuation: boolean;
}

export function startWien2kScfSession(
  projectId: string,
  cifId: string,
  structureCalculationId: string,
): Promise<Wien2kScfSession> {
  return invoke("wien2k_start_scf_session", { projectId, cifId, structureCalculationId });
}

export function startWien2kScfContinuationSession(
  projectId: string,
  cifId: string,
  calculationId: string,
): Promise<Wien2kScfSession> {
  return invoke("wien2k_start_scf_continuation_session", { projectId, cifId, calculationId });
}

export function initializeWien2kScfSession(
  sessionId: string,
  settings: Wien2kInitializationSettings,
): Promise<Wien2kInitializationResult> {
  return invoke("wien2k_initialize_scf_session", { sessionId, settings });
}

export function runWien2kScfSession(
  sessionId: string,
  settings: Wien2kScfRunSettings,
  continuation: boolean,
  parentCalculationId?: string | null,
  resources?: SlurmResourceRequest | null,
): Promise<Wien2kScfExecutionResult> {
  return invoke("wien2k_run_scf_session", {
    sessionId,
    settings,
    continuation,
    parentCalculationId: parentCalculationId ?? null,
    resources: resources ?? null,
  });
}

export function discardWien2kScfSession(sessionId: string): Promise<void> {
  return invoke("wien2k_discard_scf_session", { sessionId });
}
