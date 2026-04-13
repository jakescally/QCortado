import { useMemo } from "react";

export interface EpwSourceRef {
  calc_id?: string;
  calc_type?: string;
}

export interface EpwArtifactManifestEntry {
  source_calc_id?: string;
  source_calc_type?: string;
  rel_path?: string;
  size_bytes?: number;
}

export interface EpwNamedArtifact {
  file_name?: string;
  size_bytes?: number;
}

export interface EpwResultSummary {
  completed?: boolean;
  elapsed_seconds?: number | null;
  core_metrics?: Record<string, number>;
  generated_outputs?: EpwNamedArtifact[];
  unknown_metrics?: Record<string, unknown>;
  parse_partial?: boolean;
  notes?: string[];
}

export interface EpwErrorRecord {
  code?: string;
  message?: string;
  hint?: string | null;
}

export interface EpwViewerData {
  schema_version?: number;
  sources?: {
    phonon?: EpwSourceRef | null;
    wannier?: EpwSourceRef | null;
    scf?: EpwSourceRef | null;
    manifests?: EpwArtifactManifestEntry[];
  };
  extensions?: {
    superconductivity?: unknown;
  };
  runtime?: {
    pools?: number | null;
    max_seconds?: number | null;
    artifact_sync_mode?: string | null;
  };
  artifacts?: EpwNamedArtifact[];
  result_summary?: EpwResultSummary;
  errors?: EpwErrorRecord[];
}

export interface EpwViewerPayload {
  data: EpwViewerData;
  rawOutput?: string | null;
}

interface EpwViewerProps {
  payload: EpwViewerPayload;
  onBack: () => void;
}

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value.toFixed(value >= 10 || unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
}

function formatMetricValue(value: unknown): string {
  if (value === null || value === undefined) return "null";
  if (typeof value === "number") {
    if (!Number.isFinite(value)) return String(value);
    return Number.isInteger(value) ? String(value) : value.toFixed(6).replace(/0+$/, "").replace(/\.$/, "");
  }
  if (typeof value === "string") return value;
  if (typeof value === "boolean") return value ? "true" : "false";
  try {
    return JSON.stringify(value);
  } catch {
    return String(value);
  }
}

function sourceIdLabel(source: EpwSourceRef | null | undefined): string {
  const raw = String(source?.calc_id || "").trim();
  return raw ? raw.slice(0, 8) : "N/A";
}

export function EpwViewer({ payload, onBack }: EpwViewerProps) {
  const data = payload.data || {};
  const summary = data.result_summary || {};
  const runtime = data.runtime || {};
  const rawOutput = typeof payload.rawOutput === "string" ? payload.rawOutput : "";
  const coreMetricsEntries = Object.entries(summary.core_metrics || {});
  const unknownMetricEntries = Object.entries(summary.unknown_metrics || {});
  const notes = Array.isArray(summary.notes) ? summary.notes : [];
  const errors = Array.isArray(data.errors) ? data.errors : [];
  const manifests = Array.isArray(data.sources?.manifests) ? data.sources!.manifests! : [];

  const outputFiles = useMemo(() => {
    const merged = new Map<string, EpwNamedArtifact>();
    const push = (entry: EpwNamedArtifact | null | undefined) => {
      const name = String(entry?.file_name || "").trim();
      if (!name) return;
      merged.set(name, {
        file_name: name,
        size_bytes: Number(entry?.size_bytes) || 0,
      });
    };

    for (const entry of summary.generated_outputs || []) push(entry);
    for (const entry of data.artifacts || []) push(entry);

    return Array.from(merged.values()).sort((left, right) => String(left.file_name).localeCompare(String(right.file_name)));
  }, [data.artifacts, summary.generated_outputs]);

  return (
    <div className="bands-viewer-container">
      <div className="bands-viewer-header">
        <button className="back-button" onClick={onBack}>
          ← Back to Dashboard
        </button>
        <h2>EPW Results</h2>
      </div>

      <div className="bands-viewer-content bands-viewer-content-stacked">
        <div className="bands-viewer-details-region">
          <div className="details-grid">
            <div className="detail-item">
              <label>Completed</label>
              <span>{summary.completed ? "Yes" : "No"}</span>
            </div>
            <div className="detail-item">
              <label>Parse Partial</label>
              <span>{summary.parse_partial ? "Yes" : "No"}</span>
            </div>
            <div className="detail-item">
              <label>Elapsed (s)</label>
              <span>{summary.elapsed_seconds ?? "N/A"}</span>
            </div>
            <div className="detail-item">
              <label>Schema</label>
              <span>{data.schema_version ?? "N/A"}</span>
            </div>
            <div className="detail-item">
              <label>MPI Pools</label>
              <span>{runtime.pools ?? "N/A"}</span>
            </div>
            <div className="detail-item">
              <label>Sync Mode</label>
              <span>{runtime.artifact_sync_mode || "N/A"}</span>
            </div>
          </div>

          <div className="option-section">
            <h4>Prerequisite Provenance</h4>
            <div className="details-grid">
              <div className="detail-item">
                <label>Source Phonon</label>
                <span>{sourceIdLabel(data.sources?.phonon)}</span>
              </div>
              <div className="detail-item">
                <label>Source Wannier</label>
                <span>{sourceIdLabel(data.sources?.wannier)}</span>
              </div>
              <div className="detail-item">
                <label>Source SCF</label>
                <span>{sourceIdLabel(data.sources?.scf)}</span>
              </div>
              <div className="detail-item">
                <label>Manifest Entries</label>
                <span>{manifests.length}</span>
              </div>
            </div>
          </div>

          {coreMetricsEntries.length > 0 && (
            <div className="option-section">
              <h4>Core Metrics</h4>
              <div className="details-grid">
                {coreMetricsEntries.map(([key, value]) => (
                  <div key={key} className="detail-item">
                    <label>{key}</label>
                    <span>{formatMetricValue(value)}</span>
                  </div>
                ))}
              </div>
            </div>
          )}

          {unknownMetricEntries.length > 0 && (
            <div className="option-section">
              <h4>Unknown Metrics</h4>
              <div className="details-grid">
                {unknownMetricEntries.map(([key, value]) => (
                  <div key={key} className="detail-item">
                    <label>{key}</label>
                    <span>{formatMetricValue(value)}</span>
                  </div>
                ))}
              </div>
            </div>
          )}

          {notes.length > 0 && (
            <div className="option-section">
              <h4>Parser Notes</h4>
              <ul className="warning-list">
                {notes.map((entry, index) => <li key={`${entry}-${index}`}>{entry}</li>)}
              </ul>
            </div>
          )}

          {errors.length > 0 && (
            <div className="option-section">
              <h4>Errors</h4>
              <ul className="warning-list">
                {errors.map((entry, index) => (
                  <li key={`${entry.code || "error"}-${index}`}>
                    [{entry.code || "error"}] {entry.message || "No message"}
                    {entry.hint ? ` Hint: ${entry.hint}` : ""}
                  </li>
                ))}
              </ul>
            </div>
          )}

          <div className="option-section">
            <h4>Output Files</h4>
            {outputFiles.length > 0 ? (
              <div className="calculations-list">
                {outputFiles.map((artifact) => (
                  <div key={artifact.file_name} className="calculation-item">
                    <div className="calculation-header">
                      <div className="calculation-info">
                        <span className="calc-type">{artifact.file_name}</span>
                      </div>
                      <div className="calculation-meta">
                        <span className="calc-size">{formatBytes(Number(artifact.size_bytes) || 0)}</span>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            ) : (
              <p className="param-hint">No output artifacts recorded.</p>
            )}
          </div>

          <div className="option-section">
            <h4>Run Log</h4>
            <pre className="output-text">{rawOutput || "No saved log output available for this calculation."}</pre>
          </div>
        </div>
      </div>
    </div>
  );
}
