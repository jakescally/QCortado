export type ProgressStatus = "idle" | "running" | "error" | "complete";

export interface PhononProgressMeta {
  totalQPoints?: number;
  completedQPoints?: string[];
  inProgressQ?: string;
  hasDos?: boolean;
  hasDispersion?: boolean;
}

export interface ProgressMeta {
  phonon?: PhononProgressMeta;
}

export interface ProgressState {
  status: ProgressStatus;
  percent: number | null;
  phase: string;
  detail?: string;
  meta?: ProgressMeta;
}

const SCF_ITER_RE = /iteration\s+#\s*(\d+)/i;
const SCF_ACCURACY_RE = /estimated\s+scf\s+accuracy\s*[=<]\s*([-\d.Ee+]+)\s*(\S+)?/i;

const QPOINTS_TOTAL_RE = /\(\s*(\d+)\s*q-points?\s*\)/i;
const CALC_Q_RE = /Calculation of q\s*=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)/i;

const SCF_DONE_MARKERS = [
  "convergence has been achieved",
  "End of self-consistent calculation",
  "End of BFGS Geometry Optimization",
];

function updateWithDetail(state: ProgressState, detail?: string): ProgressState {
  if (!detail) return state;
  return { ...state, detail };
}

function normalizeQPoint(x: string, y: string, z: string): string {
  const nums = [x, y, z].map((value) => Number.parseFloat(value));
  if (nums.some((value) => Number.isNaN(value))) {
    return `${x.trim()} ${y.trim()} ${z.trim()}`;
  }
  return nums.map((value) => value.toFixed(6)).join(" ");
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
  let next: ProgressState = { ...state, status: "running" };

  const bandsStepMatch = line.match(/Step\s+(\d+)\/(\d+):\s*(.+)/i);
  if (bandsStepMatch) {
    const step = Number.parseInt(bandsStepMatch[1], 10);
    const total = Number.parseInt(bandsStepMatch[2], 10);
    if (step === 1) {
      return { ...next, percent: 10, phase: "NSCF along k-path" };
    }
    if (step === 2) {
      return {
        ...next,
        percent: total >= 3 ? 60 : 80,
        phase: "bands.x post-processing",
      };
    }
    if (step === 3) {
      return { ...next, percent: 90, phase: "projwfc.x projections" };
    }
  }
  if (line.includes("Parsing band structure data")) {
    return { ...next, percent: 95, phase: "Parsing bands" };
  }
  if (line.includes("=== Band Structure Complete ===")) {
    return { ...next, percent: 100, status: "complete", phase: "Complete" };
  }

  return next;
}

export function updateDosProgress(line: string, state: ProgressState): ProgressState {
  const next: ProgressState = { ...state, status: "running" };

  const dosStepMatch = line.match(/Step\s+(\d+)\/(\d+):\s*(.+)/i);
  if (dosStepMatch) {
    const step = Number.parseInt(dosStepMatch[1], 10);
    if (step === 1) {
      return { ...next, percent: 15, phase: "NSCF on dense k-grid" };
    }
    if (step === 2) {
      return { ...next, percent: 70, phase: "dos.x post-processing" };
    }
  }

  if (line.includes("Parsing DOS data")) {
    return { ...next, percent: 92, phase: "Parsing DOS" };
  }
  if (line.includes("=== Electronic DOS Complete ===")) {
    return { ...next, percent: 100, status: "complete", phase: "Complete" };
  }

  return next;
}

export function updateFermiSurfaceProgress(line: string, state: ProgressState): ProgressState {
  const next: ProgressState = { ...state, status: "running" };

  const stepMatch = line.match(/Step\s+(\d+)\/(\d+):\s*(.+)/i);
  if (stepMatch) {
    const step = Number.parseInt(stepMatch[1], 10);
    if (step === 1) {
      return { ...next, percent: 15, phase: "NSCF on dense k-grid" };
    }
    if (step === 2) {
      return { ...next, percent: 72, phase: "fs.x post-processing" };
    }
  }

  if (line.includes("Collecting BXSF artifacts")) {
    return { ...next, percent: 92, phase: "Collecting BXSF files" };
  }
  if (line.includes("=== Fermi Surface Generation Complete ===")) {
    return { ...next, percent: 100, status: "complete", phase: "Complete" };
  }

  return next;
}

export function updatePhononProgress(line: string, state: ProgressState): ProgressState {
  let next: ProgressState = { ...state, status: "running" };
  const meta = { ...(state.meta?.phonon ?? {}) };
  let metaUpdated = false;

  const totalMatch = line.match(QPOINTS_TOTAL_RE);
  if (totalMatch) {
    const total = Number.parseInt(totalMatch[1], 10);
    if (!Number.isNaN(total) && total > 0) {
      meta.totalQPoints = total;
      metaUpdated = true;
    }
  }

  const calcMatch = line.match(CALC_Q_RE);
  if (calcMatch) {
    const key = normalizeQPoint(calcMatch[1], calcMatch[2], calcMatch[3]);
    const completed = meta.completedQPoints ? [...meta.completedQPoints] : [];
    if (meta.inProgressQ && meta.inProgressQ !== key && !completed.includes(meta.inProgressQ)) {
      completed.push(meta.inProgressQ);
    }
    meta.inProgressQ = key;
    meta.completedQPoints = completed;
    metaUpdated = true;

    const total = meta.totalQPoints ?? 0;
    const completedCount = meta.completedQPoints?.length ?? 0;
    if (total > 0) {
      const phMax = meta.hasDos ? 95 : 99;
      const ratio = Math.min(1, completedCount / total);
      const percent = Math.max(next.percent ?? 0, ratio * phMax);
      const runningIndex = Math.min(completedCount + 1, total);
      next = updateWithDetail(
        { ...next, percent, phase: "ph.x" },
        `Q-point ${runningIndex}/${total}`
      );
    } else {
      next = { ...next, phase: "ph.x" };
    }
  }

  const attachMeta = (stateToReturn: ProgressState): ProgressState => {
    if (metaUpdated || state.meta?.phonon) {
      return { ...stateToReturn, meta: { ...(state.meta ?? {}), phonon: meta } };
    }
    return stateToReturn;
  };

  if (line.includes("=== Step 1/4")) {
    return attachMeta({ ...next, percent: next.percent ?? 0, phase: "ph.x" });
  }
  if (line.includes("=== Step 2/4")) {
    const phMax = meta.hasDos ? 95 : 99;
    const completed = meta.completedQPoints ? [...meta.completedQPoints] : [];
    if (meta.inProgressQ && !completed.includes(meta.inProgressQ)) {
      completed.push(meta.inProgressQ);
      meta.completedQPoints = completed;
      meta.inProgressQ = undefined;
      metaUpdated = true;
    }
    return attachMeta({ ...next, percent: Math.max(next.percent ?? 0, phMax), phase: "q2r.x" });
  }
  if (line.includes("=== Step 3/4")) {
    const target = meta.hasDos ? 95 : 99;
    return attachMeta({ ...next, percent: Math.max(next.percent ?? 0, target), phase: "matdyn.x (DOS)" });
  }
  if (line.includes("=== Step 4/4")) {
    const target = 100;
    return attachMeta({ ...next, percent: Math.max(next.percent ?? 0, target), phase: "matdyn.x (dispersion)" });
  }
  if (line.includes("=== Phonon Calculation Complete ===")) {
    return attachMeta({ ...next, percent: 100, status: "complete", phase: "Complete" });
  }

  return attachMeta(next);
}

type ProgressKind = "scf" | "bands" | "dos" | "fermi_surface" | "phonon";

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
    case "dos":
      return updateDosProgress(line, state);
    case "fermi_surface":
      return updateFermiSurfaceProgress(line, state);
    case "phonon":
      return updatePhononProgress(line, state);
    default:
      return state;
  }
}

export function defaultProgressState(phase: string, meta?: ProgressMeta): ProgressState {
  return {
    status: "running",
    percent: null,
    phase,
    meta,
  };
}
