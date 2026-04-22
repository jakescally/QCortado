import { useEffect, useRef, useState } from "react";
import {
  HpcCpuStepTelemetry,
  HpcUtilizationSample,
  sampleHpcUtilization,
} from "../lib/hpcConfig";
import { HpcResourceType } from "../lib/types";
import {
  calculateCpuUtilizationPercent,
  clampTelemetryPercent,
  formatTelemetryBytes,
  formatTelemetryPercent,
  formatTelemetryTime,
  getCpuUtilizationPoint,
} from "../lib/hpcUtilization";

interface UseHpcUtilizationArgs {
  enabled: boolean;
  profileId?: string | null;
  remoteJobId?: string | null;
  remoteNode?: string | null;
  resourceType: HpcResourceType;
}

interface HpcUtilizationState {
  sample: HpcUtilizationSample | null;
  loading: boolean;
  error: string | null;
  cpuUtilizationPercent: number | null;
}

interface RemoteUtilizationPanelProps {
  profileId?: string | null;
  remoteJobId?: string | null;
  remoteNode?: string | null;
  resourceType: HpcResourceType;
  enabled: boolean;
}

const POLL_INTERVAL_MS = 5000;

export function useHpcUtilization({
  enabled,
  profileId,
  remoteJobId,
  remoteNode,
  resourceType,
}: UseHpcUtilizationArgs): HpcUtilizationState {
  const [sample, setSample] = useState<HpcUtilizationSample | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [cpuUtilizationPercent, setCpuUtilizationPercent] = useState<number | null>(null);
  const previousCpuPointRef = useRef<ReturnType<typeof getCpuUtilizationPoint>>(null);

  useEffect(() => {
    previousCpuPointRef.current = null;
    setCpuUtilizationPercent(null);
    setSample(null);
    setError(null);
  }, [remoteJobId, resourceType]);

  useEffect(() => {
    if (!enabled) {
      setLoading(false);
      return;
    }
    if (!profileId) {
      setError("No active HPC profile selected.");
      setLoading(false);
      return;
    }
    if (!remoteJobId) {
      setError(null);
      setLoading(false);
      return;
    }

    let cancelled = false;
    let timeoutId: number | null = null;

    const poll = async () => {
      if (cancelled) return;
      setLoading(true);
      try {
        const next = await sampleHpcUtilization(profileId, remoteJobId, remoteNode ?? null, resourceType);
        if (cancelled) return;
        const currentPoint = getCpuUtilizationPoint(next);
        const nextCpuPercent = calculateCpuUtilizationPercent(previousCpuPointRef.current, currentPoint);
        previousCpuPointRef.current = currentPoint;
        setSample(next);
        setCpuUtilizationPercent(nextCpuPercent);
        setError(null);
      } catch (pollError) {
        if (cancelled) return;
        setError(String(pollError));
      } finally {
        if (cancelled) return;
        setLoading(false);
        timeoutId = window.setTimeout(() => {
          void poll();
        }, POLL_INTERVAL_MS);
      }
    };

    void poll();

    return () => {
      cancelled = true;
      if (timeoutId !== null) {
        window.clearTimeout(timeoutId);
      }
    };
  }, [enabled, profileId, remoteJobId, remoteNode, resourceType]);

  return { sample, loading, error, cpuUtilizationPercent };
}

function MetricCard({
  label,
  value,
  detail,
}: {
  label: string;
  value: string;
  detail?: string | null;
}) {
  return (
    <div className="telemetry-metric-card">
      <span>{label}</span>
      <strong>{value}</strong>
      {detail && <small>{detail}</small>}
    </div>
  );
}

function UtilizationBar({ value, label }: { value?: number | null; label?: string }) {
  const percent = clampTelemetryPercent(value);
  return (
    <div className="telemetry-bar" aria-label={label}>
      <div className="telemetry-bar-fill" style={{ width: `${percent}%` }} />
    </div>
  );
}

function formatRatio(used?: number | null, total?: number | null): string {
  if (!used || !total) {
    return "Unavailable";
  }
  return `${formatTelemetryBytes(used)} / ${formatTelemetryBytes(total)}`;
}

function stepMemory(step: HpcCpuStepTelemetry): string {
  return formatTelemetryBytes(step.peak_rss_bytes ?? step.average_rss_bytes ?? null);
}

export function RemoteUtilizationPanel({
  profileId,
  remoteJobId,
  remoteNode,
  resourceType,
  enabled,
}: RemoteUtilizationPanelProps) {
  const { sample, loading, error, cpuUtilizationPercent } = useHpcUtilization({
    enabled,
    profileId,
    remoteJobId,
    remoteNode,
    resourceType,
  });
  const warnings = sample?.warnings ?? [];
  const sources = sample?.sources?.length ? sample.sources.join(", ") : "pending";
  const schedulerNode = sample?.node || remoteNode || "pending";
  const currentRam = sample?.memory?.current_rss_bytes ?? null;
  const peakRam = sample?.memory?.peak_rss_bytes ?? null;
  const requestedRam = sample?.memory?.requested_bytes ?? sample?.scheduler?.requested_memory_bytes ?? null;
  const ramPercent = currentRam && requestedRam ? (currentRam / requestedRam) * 100 : null;

  return (
    <div className="telemetry-panel">
      <div className="telemetry-header">
        <h3>Remote Utilization</h3>
        <span className="telemetry-meta">
          {loading ? "Refreshing..." : formatTelemetryTime(sample?.captured_at)}
        </span>
      </div>

      <div className="telemetry-meta-row">
        <span>Job: {remoteJobId || "pending allocation"}</span>
        <span>Node: {schedulerNode}</span>
        <span>Type: {resourceType.toUpperCase()}</span>
        <span>Source: {sources}</span>
      </div>

      {!remoteJobId ? (
        <div className="telemetry-empty-state">Waiting for remote job allocation.</div>
      ) : error ? (
        <div className="telemetry-empty-state telemetry-error-state">Telemetry probe failed: {error}</div>
      ) : sample ? (
        <>
          {resourceType === "gpu" ? (
            <>
              <div className="telemetry-metric-grid">
                <MetricCard
                  label="GPU"
                  value={formatTelemetryPercent(sample.gpu?.average_utilization_percent)}
                  detail={`${sample.gpu?.devices?.length ?? 0} device${(sample.gpu?.devices?.length ?? 0) === 1 ? "" : "s"}`}
                />
                <MetricCard
                  label="VRAM"
                  value={formatTelemetryPercent(sample.gpu?.memory_used_percent)}
                  detail={formatRatio(sample.gpu?.memory_used_bytes, sample.gpu?.memory_total_bytes)}
                />
                <MetricCard
                  label="RAM"
                  value={formatTelemetryBytes(currentRam)}
                  detail={requestedRam ? `of ${formatTelemetryBytes(requestedRam)}` : sample.memory?.source ?? null}
                />
                <MetricCard label="Peak RAM" value={formatTelemetryBytes(peakRam)} detail={sample.memory?.source ?? null} />
              </div>
              <div className="telemetry-device-list">
                {(sample.gpu?.devices ?? []).map((device) => (
                  <div className="telemetry-device-row" key={device.index}>
                    <div className="telemetry-device-heading">
                      <strong>GPU {device.index}</strong>
                      <span>{device.name}</span>
                    </div>
                    <div className="telemetry-device-bars">
                      <span>Compute {formatTelemetryPercent(device.utilization_percent)}</span>
                      <UtilizationBar value={device.utilization_percent} label={`GPU ${device.index} compute`} />
                      <span>VRAM {formatRatio(device.memory_used_bytes, device.memory_total_bytes)}</span>
                      <UtilizationBar value={device.memory_utilization_percent} label={`GPU ${device.index} memory`} />
                    </div>
                    <span className="telemetry-device-temp">
                      {device.temperature_c != null ? `${device.temperature_c.toFixed(0)} C` : "Temp unavailable"}
                    </span>
                  </div>
                ))}
                {(sample.gpu?.devices?.length ?? 0) === 0 && (
                  <div className="telemetry-empty-state">GPU metrics unavailable inside the allocation.</div>
                )}
              </div>
            </>
          ) : (
            <>
              <div className="telemetry-metric-grid">
                <MetricCard
                  label="CPU"
                  value={formatTelemetryPercent(cpuUtilizationPercent)}
                  detail={cpuUtilizationPercent == null ? "Needs two accounting samples" : "Live job usage"}
                />
                <MetricCard
                  label="Allocated CPUs"
                  value={String(sample.cpu?.allocated_cpus ?? sample.scheduler?.allocated_cpus ?? "Unknown")}
                  detail={sample.scheduler?.state ?? null}
                />
                <MetricCard
                  label="RAM"
                  value={formatTelemetryBytes(currentRam)}
                  detail={requestedRam ? `of ${formatTelemetryBytes(requestedRam)}` : sample.memory?.source ?? null}
                />
                <MetricCard label="Peak RAM" value={formatTelemetryBytes(peakRam)} detail={sample.memory?.source ?? null} />
              </div>
              <div className="telemetry-bar-group">
                <div>
                  <span>CPU utilization</span>
                  <span>{formatTelemetryPercent(cpuUtilizationPercent)}</span>
                </div>
                <UtilizationBar value={cpuUtilizationPercent} label="CPU utilization" />
              </div>
              <div className="telemetry-bar-group">
                <div>
                  <span>RAM utilization</span>
                  <span>{requestedRam ? formatTelemetryPercent(ramPercent) : "Request unknown"}</span>
                </div>
                <UtilizationBar value={ramPercent} label="RAM utilization" />
              </div>
              <div className="telemetry-step-table-wrap">
                <table className="telemetry-step-table">
                  <thead>
                    <tr>
                      <th>Step</th>
                      <th>Tasks</th>
                      <th>Ave CPU</th>
                      <th>Peak RSS</th>
                    </tr>
                  </thead>
                  <tbody>
                    {(sample.cpu?.steps ?? []).slice(0, 6).map((step) => (
                      <tr key={step.job_step}>
                        <td>{step.job_step}</td>
                        <td>{step.n_tasks ?? "Unknown"}</td>
                        <td>{step.average_cpu_display ?? "Unavailable"}</td>
                        <td>{stepMemory(step)}</td>
                      </tr>
                    ))}
                    {(sample.cpu?.steps?.length ?? 0) === 0 && (
                      <tr>
                        <td colSpan={4}>No Slurm accounting rows are available yet.</td>
                      </tr>
                    )}
                  </tbody>
                </table>
              </div>
            </>
          )}

          {warnings.length > 0 && (
            <div className="telemetry-warning-list">
              {warnings.map((warning) => (
                <span key={warning}>{warning}</span>
              ))}
            </div>
          )}

          {sample.raw && (
            <details className="telemetry-raw-details">
              <summary>Raw telemetry</summary>
              <pre>{sample.raw}</pre>
            </details>
          )}
        </>
      ) : (
        <div className="telemetry-empty-state">Waiting for first telemetry sample.</div>
      )}
    </div>
  );
}
