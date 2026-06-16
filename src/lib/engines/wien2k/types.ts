/** WIEN2k-native remote case-directory workflow types. */

export type Wien2kSpinMode = "non_spin_polarized" | "spin_polarized";
export type Wien2kParallelMode = "openmp" | "kpoint";

export type Wien2kInitialSpinConfiguration = "up" | "down" | "non_magnetic";
export type Wien2kDftUDoubleCounting = "amf" | "sic" | "hmf";
export type Wien2kFermiMethod = "tetra" | "temp" | "temps";
export type Wien2kMixerMode = "MSR1" | "MSEC3" | "MSEC4" | "MSR2" | "PRATT" | "PRAT0";
export type Wien2kMixerTrust = "default" | "STIFF" | "STIFFER" | "FAST";
export type Wien2kDispersionCorrection = "none" | "dftd3" | "dftd4";

export interface Wien2kStartingMagnetization {
  siteIndex: number;
  element: string;
  configuration: Wien2kInitialSpinConfiguration;
  momentBohrMagneton: number;
}

export interface Wien2kHubbardTarget {
  siteIndex: number;
  element: string;
  manifold: string;
  orbitalL: number;
  uEv: number;
  jEv: number;
  recommended?: boolean;
  reason?: string | null;
}

export interface Wien2kDftUSettings {
  enabled: boolean;
  doubleCounting: Wien2kDftUDoubleCounting;
  targets: Wien2kHubbardTarget[];
}

export interface Wien2kMixerSettings {
  mode: Wien2kMixerMode;
  greed: number;
  history: number;
  trust: Wien2kMixerTrust;
}

export interface Wien2kExchangeCorrelationOption {
  value: number;
  label: string;
}

export const WIEN2K_EXCHANGE_CORRELATION_OPTIONS: readonly Wien2kExchangeCorrelationOption[] = [
  { value: 13, label: "PBE (vxc=13)" },
  { value: 5, label: "LDA (vxc=5)" },
  { value: 11, label: "WC (vxc=11)" },
  { value: 19, label: "PBEsol (vxc=19)" },
];

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
  rkmax: number;
  gmax: number;
  lmax: number;
  kMesh: [number, number, number];
  exchangeCorrelation: number;
  lstartEnergyCutoffRy: number;
  spinMode: Wien2kSpinMode;
  fermiMethod: Wien2kFermiMethod;
  fermiSmearingRy?: number | null;
  startingMagnetization: Wien2kStartingMagnetization[];
}

export interface Wien2kScfRunSettings {
  spinMode: Wien2kSpinMode;
  parallelMode: Wien2kParallelMode;
  maxIterations: number;
  energyConvergenceRy: number;
  chargeConvergence: number;
  forceConvergenceMryBohr?: number | null;
  dftU: Wien2kDftUSettings;
  dispersionCorrection: Wien2kDispersionCorrection;
  iterativeDiagonalization: boolean;
  forceMinimization: boolean;
  mixer: Wien2kMixerSettings;
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
  rkmax: 7,
  gmax: 12,
  lmax: 10,
  kMesh: [6, 6, 6],
  exchangeCorrelation: 13,
  lstartEnergyCutoffRy: -6,
  spinMode: "non_spin_polarized",
  fermiMethod: "tetra",
  fermiSmearingRy: null,
  startingMagnetization: [],
};

export const DEFAULT_WIEN2K_SCF_RUN_SETTINGS: Wien2kScfRunSettings = {
  spinMode: "non_spin_polarized",
  parallelMode: "openmp",
  maxIterations: 40,
  energyConvergenceRy: 0.0001,
  chargeConvergence: 0.0001,
  forceConvergenceMryBohr: null,
  dftU: {
    enabled: false,
    doubleCounting: "sic",
    targets: [],
  },
  dispersionCorrection: "none",
  iterativeDiagonalization: false,
  forceMinimization: false,
  mixer: {
    mode: "MSR1",
    greed: 0.2,
    history: 8,
    trust: "default",
  },
};
