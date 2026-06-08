import type {
  Wien2kCaseArtifact,
  Wien2kCaseArtifactRole,
  Wien2kCasePhase,
  Wien2kCaseReference,
  Wien2kCommandPlan,
  Wien2kInitializationSettings,
  Wien2kScfRunSettings,
} from "./types";
import { WIEN2K_EXCHANGE_CORRELATION_OPTIONS } from "./types";

export function normalizeWien2kCaseName(raw: string): string | null {
  const trimmed = raw.trim();
  if (trimmed.length === 0) return null;
  return /^[A-Za-z0-9_.-]+$/.test(trimmed) ? trimmed : null;
}

export function wien2kCaseFile(caseName: string, suffix: string): string {
  return `${caseName}.${suffix}`;
}

export function isTerminalWien2kCasePhase(phase: Wien2kCasePhase): boolean {
  return phase === "scf_complete" || phase === "bands_complete" || phase === "failed";
}

export function coreWien2kCaseArtifacts(caseName: string): Wien2kCaseArtifact[] {
  return [
    artifact("struct", caseName, "struct", true),
    artifact("initialization_output", caseName, "outputnn", false),
    artifact("initialization_output", caseName, "outputst", false),
    artifact("initialization_output", caseName, "outputs", false),
    artifact("density", caseName, "clmsum", true),
    artifact("scf_output", caseName, "scf", true),
    artifact("dayfile", caseName, "dayfile", false),
  ];
}

export function initializedWien2kCaseArtifacts(caseName: string): Wien2kCaseArtifact[] {
  return [
    artifact("struct", caseName, "struct", true),
    artifact("initialization_output", caseName, "outputnn", false),
    artifact("initialization_output", caseName, "outputs", false),
    artifact("initialization_output", caseName, "outputst", false),
    artifact("density", caseName, "clmsum", true),
    artifact("dayfile", caseName, "dayfile", false),
    artifact("initialization_input", caseName, "in0", true),
    artifact("initialization_input", caseName, "in1", true),
    artifact("initialization_input", caseName, "in2", true),
    artifact("initialization_input", caseName, "inc", true),
    artifact("initialization_input", caseName, "inm", true),
    artifact("initialization_input", caseName, "klist", true),
  ];
}

export function buildWien2kInitCommandPlan(
  caseRef: Wien2kCaseReference,
  settings: Wien2kInitializationSettings,
): Wien2kCommandPlan {
  const argv = [
    "-b",
    "-rkmax",
    formatNumber(settings.rkmax),
    "-gmax",
    formatNumber(settings.gmax),
    "-lmax",
    String(settings.lmax),
    "-numk",
    "0",
    String(settings.kMesh[0]),
    String(settings.kMesh[1]),
    String(settings.kMesh[2]),
    "-vxc",
    String(settings.exchangeCorrelation),
    "-ecut",
    formatNumber(settings.lstartEnergyCutoffRy),
  ];

  if (settings.spinMode === "spin_polarized") {
    argv.push("-sp");
  }

  return {
    program: "init_lapw",
    argv,
    workingDirectory: caseRef.remoteCaseDir,
    environment: caseRef.remoteScratchDir ? [["SCRATCH", caseRef.remoteScratchDir]] : [],
    phase: "initialization_running",
    expectedArtifacts: initializedWien2kCaseArtifacts(caseRef.caseName),
  };
}

export function buildWien2kScfCommandPlan(
  caseRef: Wien2kCaseReference,
  settings: Wien2kScfRunSettings,
  continuation = false,
): Wien2kCommandPlan {
  const argv = [
    "-ec",
    formatNumber(settings.energyConvergenceRy),
    "-cc",
    formatNumber(settings.chargeConvergence),
    "-i",
    String(settings.maxIterations),
  ];

  if (settings.forceConvergenceMryBohr != null) {
    argv.push("-fc", formatNumber(settings.forceConvergenceMryBohr));
  }
  if (settings.iterativeDiagonalization) {
    argv.push("-it");
  }
  if (settings.forceMinimization) {
    argv.push("-min");
  }
  if (settings.dispersionCorrection !== "none") {
    argv.push(`-${settings.dispersionCorrection}`);
  }
  if (settings.dftU.enabled) {
    argv.push("-orb");
  }
  if (continuation) {
    argv.push("-NI");
  }

  return {
    program: settings.spinMode === "spin_polarized" ? "runsp_lapw" : "run_lapw",
    argv,
    workingDirectory: caseRef.remoteCaseDir,
    environment: caseRef.remoteScratchDir ? [["SCRATCH", caseRef.remoteScratchDir]] : [],
    phase: "scf_running",
    expectedArtifacts: coreWien2kCaseArtifacts(caseRef.caseName),
  };
}

export function validateWien2kInitializationSettings(settings: Wien2kInitializationSettings): string | null {
  if (!Number.isFinite(settings.rkmax) || settings.rkmax <= 0) return "RKMAX must be positive.";
  if (!Number.isFinite(settings.gmax) || settings.gmax <= 0) return "GMAX must be positive.";
  if (!Number.isInteger(settings.lmax) || settings.lmax <= 0) return "LMAX must be a positive integer.";
  if (settings.kMesh.some((value) => !Number.isInteger(value) || value <= 0)) {
    return "All k-mesh dimensions must be positive integers.";
  }
  if (!WIEN2K_EXCHANGE_CORRELATION_OPTIONS.some((option) => option.value === settings.exchangeCorrelation)) {
    return "XC functional must be one of WIEN2k's native initialization options: LDA, PBE, WC, or PBEsol.";
  }
  if (!Number.isFinite(settings.lstartEnergyCutoffRy)) return "The LSTART energy cutoff must be finite.";
  if (!["tetra", "temp", "temps"].includes(settings.fermiMethod)) return "Fermi integration method is not supported.";
  if (settings.fermiMethod !== "tetra"
    && (settings.fermiSmearingRy == null || !Number.isFinite(settings.fermiSmearingRy) || settings.fermiSmearingRy <= 0)) {
    return "Fermi smearing must be positive when TEMP or TEMPS is selected.";
  }
  for (const entry of settings.startingMagnetization ?? []) {
    if (!Number.isInteger(entry.siteIndex) || entry.siteIndex <= 0) return "Starting magnetization site indices must be positive integers.";
    if (!["up", "down", "non_magnetic"].includes(entry.configuration)) return "Starting spin configuration is not supported.";
    if (!Number.isFinite(entry.momentBohrMagneton) || entry.momentBohrMagneton < 0) return "Starting magnetization moments must be non-negative.";
  }
  return null;
}

export function validateWien2kScfRunSettings(settings: Wien2kScfRunSettings): string | null {
  if (!Number.isInteger(settings.maxIterations) || settings.maxIterations <= 0) {
    return "Maximum SCF iterations must be a positive integer.";
  }
  if (!Number.isFinite(settings.energyConvergenceRy) || settings.energyConvergenceRy <= 0) {
    return "Energy convergence must be positive.";
  }
  if (!Number.isFinite(settings.chargeConvergence) || settings.chargeConvergence <= 0) {
    return "Charge convergence must be positive.";
  }
  if (settings.forceConvergenceMryBohr != null
    && (!Number.isFinite(settings.forceConvergenceMryBohr) || settings.forceConvergenceMryBohr <= 0)) {
    return "Force convergence must be positive when set.";
  }
  if (settings.dftU.enabled) {
    if (settings.spinMode !== "spin_polarized") return "WIEN2k DFT+U requires spin-polarized SCF.";
    if (!["amf", "sic", "hmf"].includes(settings.dftU.doubleCounting)) return "DFT+U double-counting mode is not supported.";
    if (settings.dftU.targets.length === 0) return "Enable at least one DFT+U target.";
    for (const target of settings.dftU.targets) {
      if (!Number.isInteger(target.siteIndex) || target.siteIndex <= 0) return "DFT+U site indices must be positive integers.";
      if (!/^\d+[spdf]$/.test(target.manifold)) return "DFT+U manifolds must look like 3d, 4f, etc.";
      if (!Number.isInteger(target.orbitalL) || target.orbitalL < 0 || target.orbitalL > 3) return "DFT+U orbital l must be 0, 1, 2, or 3.";
      if (!Number.isFinite(target.uEv) || target.uEv <= 0) return "DFT+U U values must be positive eV values.";
      if (!Number.isFinite(target.jEv) || target.jEv < 0) return "DFT+U J values must be non-negative eV values.";
    }
  }
  if (!["none", "dftd3", "dftd4"].includes(settings.dispersionCorrection)) {
    return "Dispersion correction is not supported.";
  }
  if (!["MSR1", "MSEC3", "MSEC4", "MSR2", "PRATT", "PRAT0"].includes(settings.mixer.mode)) {
    return "Mixer mode is not supported.";
  }
  if (!Number.isFinite(settings.mixer.greed) || settings.mixer.greed <= 0 || settings.mixer.greed > 1) {
    return "Mixer greed must be in the range (0, 1].";
  }
  if (!Number.isInteger(settings.mixer.history) || settings.mixer.history <= 0) {
    return "Mixer history must be a positive integer.";
  }
  if (!["default", "STIFF", "STIFFER", "FAST"].includes(settings.mixer.trust)) {
    return "Mixer trust setting is not supported.";
  }
  return null;
}

function artifact(
  role: Wien2kCaseArtifactRole,
  caseName: string,
  suffix: string,
  requiredForResume: boolean,
): Wien2kCaseArtifact {
  return {
    role,
    basename: wien2kCaseFile(caseName, suffix),
    requiredForResume,
  };
}

function formatNumber(value: number): string {
  if (!Number.isFinite(value)) return String(value);
  return String(Number(value.toFixed(8)));
}
