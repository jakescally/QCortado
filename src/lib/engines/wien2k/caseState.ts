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
