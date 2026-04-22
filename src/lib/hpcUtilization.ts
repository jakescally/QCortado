import { HpcUtilizationSample } from "./hpcConfig";
import { HpcResourceType } from "./types";

export interface CpuUtilizationPoint {
  jobId?: string | null;
  capturedAt?: string | null;
  totalCpuSeconds?: number | null;
  allocatedCpus?: number | null;
}

export function formatTelemetryBytes(value?: number | null): string {
  if (!Number.isFinite(value) || value == null || value <= 0) {
    return "Unavailable";
  }
  const units = ["B", "KB", "MB", "GB", "TB"];
  let scaled = value;
  let unitIndex = 0;
  while (scaled >= 1024 && unitIndex < units.length - 1) {
    scaled /= 1024;
    unitIndex += 1;
  }
  const precision = scaled >= 100 || unitIndex === 0 ? 0 : 1;
  return `${scaled.toFixed(precision)} ${units[unitIndex]}`;
}

export function formatTelemetryPercent(value?: number | null): string {
  if (!Number.isFinite(value) || value == null) {
    return "Unavailable";
  }
  return `${Math.max(0, value).toFixed(0)}%`;
}

export function clampTelemetryPercent(value?: number | null): number {
  if (!Number.isFinite(value) || value == null) {
    return 0;
  }
  return Math.max(0, Math.min(100, value));
}

export function formatTelemetryTime(timestamp?: string | null): string {
  if (!timestamp) return "Waiting for first sample";
  const date = new Date(timestamp);
  if (Number.isNaN(date.getTime())) return timestamp;
  return `Updated ${date.toLocaleTimeString()}`;
}

export function getCpuUtilizationPoint(sample: HpcUtilizationSample | null): CpuUtilizationPoint | null {
  if (!sample?.cpu?.total_cpu_seconds || !sample.captured_at) {
    return null;
  }
  return {
    jobId: sample.job_id,
    capturedAt: sample.captured_at,
    totalCpuSeconds: sample.cpu.total_cpu_seconds,
    allocatedCpus: sample.cpu.allocated_cpus ?? sample.scheduler?.allocated_cpus ?? null,
  };
}

export function calculateCpuUtilizationPercent(
  previous: CpuUtilizationPoint | null,
  current: CpuUtilizationPoint | null,
): number | null {
  if (!previous || !current || previous.jobId !== current.jobId) {
    return null;
  }
  const previousTime = previous.capturedAt ? new Date(previous.capturedAt).getTime() : NaN;
  const currentTime = current.capturedAt ? new Date(current.capturedAt).getTime() : NaN;
  const previousCpu = previous.totalCpuSeconds;
  const currentCpu = current.totalCpuSeconds;
  const cpus = current.allocatedCpus;
  if (
    !Number.isFinite(previousTime)
    || !Number.isFinite(currentTime)
    || !Number.isFinite(previousCpu)
    || !Number.isFinite(currentCpu)
    || !Number.isFinite(cpus)
    || previousCpu == null
    || currentCpu == null
    || cpus == null
    || cpus <= 0
  ) {
    return null;
  }
  const wallSeconds = (currentTime - previousTime) / 1000;
  const cpuDelta = currentCpu - previousCpu;
  if (wallSeconds <= 0 || cpuDelta < 0) {
    return null;
  }
  return (cpuDelta / (wallSeconds * cpus)) * 100;
}

export function resolveTelemetryResourceType(
  taskResourceType?: HpcResourceType | null,
  fallbackResourceType?: HpcResourceType | null,
): HpcResourceType {
  if (taskResourceType === "gpu" || taskResourceType === "cpu") {
    return taskResourceType;
  }
  return fallbackResourceType === "gpu" ? "gpu" : "cpu";
}
