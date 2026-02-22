import { ProgressState } from "../lib/qeProgress";

interface ProgressBarProps extends ProgressState {
  compact?: boolean;
}

export function ProgressBar({ status, percent, phase, detail, compact = false }: ProgressBarProps) {
  const isIndeterminate = percent === null;
  const clamped = percent === null ? 0 : Math.max(0, Math.min(100, percent));

  return (
    <div className={`progress-bar ${status} ${isIndeterminate ? "indeterminate" : ""} ${compact ? "compact" : ""}`}>
      <div className="progress-bar-header">
        <span className="progress-bar-phase">{phase}</span>
        {compact && detail && (
          <span className="progress-bar-detail-inline" title={detail}>{detail}</span>
        )}
        {percent !== null && (
          <span className="progress-bar-percent">{Math.round(clamped)}%</span>
        )}
      </div>
      <div className="progress-bar-track">
        <div
          className="progress-bar-fill"
          style={isIndeterminate ? undefined : { width: `${clamped}%` }}
        ></div>
      </div>
      {!compact && detail && <div className="progress-bar-detail">{detail}</div>}
    </div>
  );
}
