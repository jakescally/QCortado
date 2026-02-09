export type ProgressStatus = "idle" | "running" | "error" | "complete";

export interface ProgressState {
  status: ProgressStatus;
  percent: number | null;
  phase: string;
  detail?: string;
}

const SCF_ITER_RE = /iteration\s+#\s*(\d+)/i;
const SCF_ACCURACY_RE = /estimated\s+scf\s+accuracy\s*[=<]\s*([-\d.Ee+]+)\s*(\S+)?/i;

const SCF_DONE_MARKERS = [
  "convergence has been achieved",
  "End of self-consistent calculation",
  "End of BFGS Geometry Optimization",
];

function updateWithDetail(state: ProgressState, detail?: string): ProgressState {
  if (!detail) return state;
  return { ...state, detail };
}

export function updateScfProgress(line: string, state: ProgressState): ProgressState {
  let next: ProgressState = { ...state, status: "running", phase: "SCF iterations", percent: null };

  let iteration: string | null = null;
  const iterMatch = line.match(SCF_ITER_RE);
  if (iterMatch) {
    iteration = iterMatch[1];
  }

  let accuracy: string | null = null;
  const accMatch = line.match(SCF_ACCURACY_RE);
  if (accMatch) {
    const value = accMatch[1];
    const unit = accMatch[2] ? ` ${accMatch[2]}` : "";
    accuracy = `${value}${unit}`;
  }

  if (iteration || accuracy) {
    const pieces: string[] = [];
    if (iteration) pieces.push(`Iteration ${iteration}`);
    if (accuracy) pieces.push(`accuracy = ${accuracy}`);
    next = updateWithDetail(next, pieces.join(" • "));
  }

  for (const marker of SCF_DONE_MARKERS) {
    if (line.includes(marker)) {
      return {
        ...next,
        status: "complete",
        percent: 100,
        phase: "Complete",
      };
    }
  }

  return next;
}

export function updateBandsProgress(line: string, state: ProgressState): ProgressState {
  let next = { ...state, status: "running" };

  if (line.includes("Step 1/2: Running NSCF")) {
    return { ...next, percent: 10, phase: "NSCF along k-path" };
  }
  if (line.includes("Step 2/2: Running bands.x")) {
    return { ...next, percent: 80, phase: "bands.x post-processing" };
  }
  if (line.includes("Parsing band structure data")) {
    return { ...next, percent: 95, phase: "Parsing bands" };
  }
  if (line.includes("=== Band Structure Complete ===")) {
    return { ...next, percent: 100, status: "complete", phase: "Complete" };
  }

  return next;
}

export function updatePhononProgress(line: string, state: ProgressState): ProgressState {
  let next = { ...state, status: "running" };

  if (line.includes("=== Step 1/4")) {
    return { ...next, percent: 5, phase: "ph.x" };
  }
  if (line.includes("=== Step 2/4")) {
    return { ...next, percent: 35, phase: "q2r.x" };
  }
  if (line.includes("=== Step 3/4")) {
    return { ...next, percent: 65, phase: "matdyn.x (DOS)" };
  }
  if (line.includes("=== Step 4/4")) {
    return { ...next, percent: 85, phase: "matdyn.x (dispersion)" };
  }
  if (line.includes("=== Phonon Calculation Complete ===")) {
    return { ...next, percent: 100, status: "complete", phase: "Complete" };
  }

  return next;
}

type ProgressKind = "scf" | "bands" | "phonon";

export function progressReducer(
  kind: ProgressKind,
  line: string,
  state: ProgressState,
): ProgressState {
  switch (kind) {
    case "scf":
      return updateScfProgress(line, state);
    case "bands":
      return updateBandsProgress(line, state);
    case "phonon":
      return updatePhononProgress(line, state);
    default:
      return state;
  }
}

export function defaultProgressState(phase: string): ProgressState {
  return {
    status: "running",
    percent: null,
    phase,
  };
}
