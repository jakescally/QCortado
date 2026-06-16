import { invoke } from "@tauri-apps/api/core";
import type { BandData } from "../../../components/BandPlot";
import type { SlurmResourceRequest } from "../../types";
import type { Wien2kSpinMode } from "./types";

export type Wien2kBandsSessionPhase = "staged" | "prepared" | "bands_complete" | "failed";
export type Wien2kBandsSpinChannel = "none" | "up" | "down";

export interface Wien2kBandsKPathPoint {
  label: string;
  coords: [number, number, number];
  npoints: number;
}

export interface Wien2kBandsPrepareSettings {
  kPath: Wien2kBandsKPathPoint[];
  energyMinEv: number;
  energyMaxEv: number;
  characterAtom: number;
  characterL: number;
  characterScale: number;
  runLapw2Qtl: boolean;
  runIrrep: boolean;
  spinChannel?: Wien2kBandsSpinChannel | null;
}

export interface Wien2kBandsRunSettings {
  spinChannel: Wien2kBandsSpinChannel;
  runLapw2Qtl: boolean;
  runIrrep: boolean;
  spinOrbit: boolean;
  diagnosticLog: boolean;
}

export interface Wien2kBandsSession {
  sessionId: string;
  projectId: string;
  cifId: string;
  sourceScfCalculationId: string;
  caseName: string;
  remoteCaseDir: string;
  sourceRemoteCaseDir: string;
  remoteInstallRoot: string;
  hpcProfileId: string;
  spinMode: Wien2kSpinMode;
  sourceSpinOrbit: boolean;
  fermiEnergyEv?: number | null;
  phase: Wien2kBandsSessionPhase;
  latestPrepare?: Wien2kBandsPrepareSettings | null;
  artifacts: Record<string, string>;
  transcript: string[];
  startedAt: string;
}

export interface Wien2kBandsPrepareResult {
  sessionId: string;
  phase: Wien2kBandsSessionPhase;
  nativeOutput: string;
  artifacts: Record<string, string>;
}

export interface Wien2kBandsExecutionResult {
  sessionId: string;
  phase: Wien2kBandsSessionPhase;
  nativeOutput: string;
  diagnostics: string[];
  bandData: BandData;
  bandDataset: unknown;
  calculationId: string;
}

export function startWien2kBandsSession(
  projectId: string,
  cifId: string,
  sourceScfCalculationId: string,
): Promise<Wien2kBandsSession> {
  return invoke("wien2k_start_bands_session", { projectId, cifId, sourceScfCalculationId });
}

export function prepareWien2kBandsSession(
  sessionId: string,
  settings: Wien2kBandsPrepareSettings,
): Promise<Wien2kBandsPrepareResult> {
  return invoke("wien2k_prepare_bands_session", { sessionId, settings });
}

export function runWien2kBandsSession(
  sessionId: string,
  settings: Wien2kBandsRunSettings,
  resources?: SlurmResourceRequest | null,
): Promise<Wien2kBandsExecutionResult> {
  return invoke("wien2k_run_bands_session", {
    sessionId,
    settings,
    resources: resources ?? null,
  });
}

export function discardWien2kBandsSession(sessionId: string): Promise<void> {
  return invoke("wien2k_discard_bands_session", { sessionId });
}
