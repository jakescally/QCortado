import type {
  Wien2kCaseArtifact,
  Wien2kCaseArtifactRole,
  Wien2kCasePhase,
  Wien2kCaseReference,
  Wien2kCommandPlan,
  Wien2kInitializationSettings,
  Wien2kScfRunSettings,
} from "./types";

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
    artifact("initialization_output", caseName, "outputsgroup", false),
    artifact("initialization_output", caseName, "outputs", false),
    artifact("density", caseName, "clmsum", true),
    artifact("scf_output", caseName, "scf", true),
    artifact("dayfile", caseName, "dayfile", false),
  ];
}

export function initializedWien2kCaseArtifacts(caseName: string): Wien2kCaseArtifact[] {
  return [
    ...coreWien2kCaseArtifacts(caseName),
    artifact("initialization_input", caseName, "in0", true),
    artifact("initialization_input", caseName, "in1", true),
    artifact("initialization_input", caseName, "in2", true),
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

  if (settings.rmtReductionPercent != null) {
    argv.push("-red", formatNumber(settings.rmtReductionPercent));
  }
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
  if (settings.parallel) {
    argv.push("-p");
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
