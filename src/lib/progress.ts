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

export function defaultProgressState(phase: string, meta?: ProgressMeta): ProgressState {
  return {
    status: "running",
    percent: null,
    phase,
    meta,
  };
}
