import { ProgressState } from "../lib/qeProgress";

interface ProgressBarProps extends ProgressState {}

export function ProgressBar({ status, percent, phase, detail }: ProgressBarProps) {
  const isIndeterminate = percent === null;
  const clamped = percent === null ? 0 : Math.max(0, Math.min(100, percent));

  return (
    <div className={`progress-bar ${status} ${isIndeterminate ? "indeterminate" : ""}`}>
      <div className="progress-bar-header">
        <span className="progress-bar-phase">{phase}</span>
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
      {detail && <div className="progress-bar-detail">{detail}</div>}
    </div>
  );
}
