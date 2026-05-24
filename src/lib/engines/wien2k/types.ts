/**
 * Hidden WIEN2k skeleton types.
 *
 * These model WIEN2k's remote case-directory workflow. They are not wired to
 * UI routes or Tauri commands, and they deliberately do not reuse QE input
 * concepts such as pseudopotentials or ecutwfc/ecutrho.
 */

export type Wien2kSpinMode = "non_spin_polarized" | "spin_polarized";

export type Wien2kCasePhase =
  | "unstaged"
  | "struct_staged"
  | "initialization_running"
  | "initialized"
  | "scf_running"
  | "scf_complete"
  | "bands_running"
  | "bands_complete"
  | "failed";

export type Wien2kCaseArtifactRole =
  | "struct"
  | "initialization_input"
  | "initialization_output"
  | "scf_output"
  | "dayfile"
  | "density"
  | "vector"
  | "bands_output"
  | "scratch";

export type Wien2kCommandProgram = "init_lapw" | "run_lapw" | "runsp_lapw" | "x";

export interface Wien2kCaseReference {
  caseName: string;
  remoteCaseDir: string;
  remoteScratchDir?: string | null;
  remoteArchiveDir?: string | null;
  localShadowDir?: string | null;
  projectId?: string | null;
  cifId?: string | null;
}

export interface Wien2kCaseArtifact {
  role: Wien2kCaseArtifactRole;
  basename: string;
  requiredForResume: boolean;
}

export interface Wien2kCaseState {
  case: Wien2kCaseReference;
  phase: Wien2kCasePhase;
  artifacts: Wien2kCaseArtifact[];
  remoteJobId?: string | null;
}

export interface Wien2kInitializationSettings {
  rmtReductionPercent?: number | null;
  rkmax: number;
  gmax: number;
  lmax: number;
  kMesh: [number, number, number];
  exchangeCorrelation: number;
  lstartEnergyCutoffRy: number;
  spinMode: Wien2kSpinMode;
}

export interface Wien2kScfRunSettings {
  spinMode: Wien2kSpinMode;
  maxIterations: number;
  energyConvergenceRy: number;
  chargeConvergence: number;
  forceConvergenceMryBohr?: number | null;
  parallel: boolean;
}

export interface Wien2kRemoteRuntimeProfile {
  profileId: string;
  wienroot: string;
  remoteWorkspaceRoot: string;
  remoteProjectRoot: string;
  remoteScratchRoot?: string | null;
}

export interface Wien2kCommandPlan {
  program: Wien2kCommandProgram;
  argv: string[];
  workingDirectory: string;
  environment: Array<[string, string]>;
  phase: Wien2kCasePhase;
  expectedArtifacts: Wien2kCaseArtifact[];
}

export const DEFAULT_WIEN2K_INITIALIZATION_SETTINGS: Wien2kInitializationSettings = {
  rmtReductionPercent: null,
  rkmax: 7,
  gmax: 12,
  lmax: 10,
  kMesh: [6, 6, 6],
  exchangeCorrelation: 13,
  lstartEnergyCutoffRy: -6,
  spinMode: "non_spin_polarized",
};

export const DEFAULT_WIEN2K_SCF_RUN_SETTINGS: Wien2kScfRunSettings = {
  spinMode: "non_spin_polarized",
  maxIterations: 40,
  energyConvergenceRy: 0.0001,
  chargeConvergence: 0.0001,
  forceConvergenceMryBohr: null,
  parallel: false,
};
