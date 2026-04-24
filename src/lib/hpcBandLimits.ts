import { SlurmResourceRequest } from "./types";

export function requestedHpcTaskCount(resources: SlurmResourceRequest | null | undefined): number {
  const parsed = Number(resources?.ntasks);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    return 1;
  }
  return Math.max(1, Math.floor(parsed));
}

export function validateHpcTasksWithinBandCount(
  resources: SlurmResourceRequest | null | undefined,
  bandCount: number | null | undefined,
  bandLabel = "bands",
): string | null {
  if (bandCount == null || !Number.isFinite(bandCount) || bandCount <= 0) {
    return null;
  }

  const resolvedBandCount = Math.floor(bandCount);
  const taskCount = requestedHpcTaskCount(resources);
  if (taskCount <= resolvedBandCount) {
    return null;
  }

  return `HPC tasks (${taskCount}) cannot exceed manually requested ${bandLabel} (${resolvedBandCount}). Lower Tasks or increase ${bandLabel}.`;
}
