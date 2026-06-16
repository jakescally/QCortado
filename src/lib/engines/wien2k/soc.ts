import { invoke } from "@tauri-apps/api/core";
import type { SlurmResourceRequest } from "../../types";
import type { NormalizedScfSummary } from "./scf";
import type { Wien2kScfRunSettings, Wien2kSpinMode } from "./types";

export type Wien2kSocSessionPhase =
  | "staged"
  | "prepared"
  | "symmetry_ready"
  | "soc_complete"
  | "failed";

export interface Wien2kSocRlo {
  atomIndex: number;
  energyRy: number;
  de: number;
  switch: "STOP" | "CONT";
}

export interface Wien2kSocRloCandidate extends Wien2kSocRlo {
  sourceFile: string;
}

export interface Wien2kSocPrepareSettings {
  magnetizationDirection: [number, number, number];
  disabledAtomIndices: number[];
  lapw1EmaxRy: number;
  outputEnergyMinRy: number;
  outputEnergyMaxRy: number;
  rloAtoms: Wien2kSocRlo[];
}

export interface Wien2kSocRunSettings {
  scf: Wien2kScfRunSettings;
  diagnosticLog: boolean;
}

export interface Wien2kSocSession {
  sessionId: string;
  projectId: string;
  cifId: string;
  sourceScfCalculationId: string;
  sourceStructureCalculationId: string;
  caseName: string;
  remoteCaseDir: string;
  sourceRemoteCaseDir: string;
  remoteInstallRoot: string;
  hpcProfileId: string;
  spinMode: Wien2kSpinMode;
  sourceKMesh?: [number, number, number] | null;
  sourceDftU: boolean;
  phase: Wien2kSocSessionPhase;
  latestPrepare?: Wien2kSocPrepareSettings | null;
  latestRun?: Wien2kSocRunSettings | null;
  latestCalculationId?: string | null;
  artifacts: Record<string, string>;
  transcript: string[];
  startedAt: string;
}

export interface Wien2kSocPrepareResult {
  sessionId: string;
  phase: Wien2kSocSessionPhase;
  symmetryReviewRequired: boolean;
  nativeOutput: string;
  diagnostics: string[];
  artifacts: Record<string, string>;
}

export interface Wien2kSocExecutionResult {
  sessionId: string;
  phase: Wien2kSocSessionPhase;
  nativeOutput: string;
  diagnostics: string[];
  summary: NormalizedScfSummary;
  calculationId: string;
  continuation: boolean;
}

export const DEFAULT_WIEN2K_SOC_PREPARE_SETTINGS: Wien2kSocPrepareSettings = {
  magnetizationDirection: [0, 0, 1],
  disabledAtomIndices: [],
  lapw1EmaxRy: 5,
  outputEnergyMinRy: -10,
  outputEnergyMaxRy: 1.9,
  rloAtoms: [],
};

export function getWien2kSocRloCandidates(
  projectId: string,
  sourceScfCalculationId: string,
): Promise<Wien2kSocRloCandidate[]> {
  return invoke("wien2k_get_soc_rlo_candidates", { projectId, sourceScfCalculationId });
}

export function startWien2kSocSession(
  projectId: string,
  cifId: string,
  sourceScfCalculationId: string,
): Promise<Wien2kSocSession> {
  return invoke("wien2k_start_soc_session", { projectId, cifId, sourceScfCalculationId });
}

export function prepareWien2kSocSession(
  sessionId: string,
  settings: Wien2kSocPrepareSettings,
): Promise<Wien2kSocPrepareResult> {
  return invoke("wien2k_prepare_soc_session", { sessionId, settings });
}

export function acceptWien2kSocSymmetry(sessionId: string): Promise<Wien2kSocPrepareResult> {
  return invoke("wien2k_accept_soc_symmetry", { sessionId });
}

export function runWien2kSocSession(
  sessionId: string,
  settings: Wien2kSocRunSettings,
  continuation: boolean,
  parentCalculationId?: string | null,
  resources?: SlurmResourceRequest | null,
): Promise<Wien2kSocExecutionResult> {
  return invoke("wien2k_run_soc_session", {
    sessionId,
    settings,
    continuation,
    parentCalculationId: parentCalculationId ?? null,
    resources: resources ?? null,
  });
}

export function discardWien2kSocSession(sessionId: string): Promise<void> {
  return invoke("wien2k_discard_soc_session", { sessionId });
}
