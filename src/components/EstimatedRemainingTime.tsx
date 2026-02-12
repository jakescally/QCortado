import { useEffect, useState } from "react";

interface EstimatedRemainingTimeProps {
  startedAt?: string | null;
  isRunning: boolean;
  completedUnits: number;
  totalUnits?: number;
  label?: string;
}

function getElapsedSeconds(startedAt: string): number {
  const startMs = Date.parse(startedAt);
  if (Number.isNaN(startMs)) return 0;
  return Math.max(0, Math.floor((Date.now() - startMs) / 1000));
}

function formatDuration(seconds: number): string {
  const total = Math.max(0, Math.floor(seconds));
  const hours = Math.floor(total / 3600);
  const minutes = Math.floor((total % 3600) / 60);
  const secs = total % 60;
  return `${String(hours).padStart(2, "0")}:${String(minutes).padStart(2, "0")}:${String(secs).padStart(2, "0")}`;
}

function computeRemainingSeconds(
  startedAt: string,
  completedUnits: number,
  totalUnits: number,
): number | null {
  if (completedUnits < 1 || totalUnits <= completedUnits) return null;
  const elapsed = getElapsedSeconds(startedAt);
  if (elapsed <= 0) return null;
  const avgPerUnit = elapsed / completedUnits;
  const remainingUnits = totalUnits - completedUnits;
  return Math.max(0, Math.round(avgPerUnit * remainingUnits));
}

export function EstimatedRemainingTime({
  startedAt,
  isRunning,
  completedUnits,
  totalUnits,
  label = "ETA",
}: EstimatedRemainingTimeProps) {
  const [remainingSeconds, setRemainingSeconds] = useState<number | null>(null);
  const hasStart = typeof startedAt === "string" && startedAt.trim().length > 0;
  const hasTotals = typeof totalUnits === "number" && totalUnits > 0;

  useEffect(() => {
    if (!hasStart || !startedAt || !hasTotals || totalUnits === undefined) {
      setRemainingSeconds(null);
      return;
    }
    setRemainingSeconds(computeRemainingSeconds(startedAt, completedUnits, totalUnits));
  }, [hasStart, startedAt, hasTotals, completedUnits, totalUnits]);

  useEffect(() => {
    if (!hasStart || !startedAt || !hasTotals || totalUnits === undefined || !isRunning) {
      return;
    }

    const intervalId = window.setInterval(() => {
      setRemainingSeconds(computeRemainingSeconds(startedAt, completedUnits, totalUnits));
    }, 1000);

    return () => {
      window.clearInterval(intervalId);
    };
  }, [hasStart, startedAt, hasTotals, completedUnits, totalUnits, isRunning]);

  if (!isRunning || !hasStart || !hasTotals || totalUnits === undefined || remainingSeconds === null) {
    return null;
  }

  return (
    <div className="eta-time" aria-live="polite">
      <span className="eta-time-label">{label}:</span>
      <span className="eta-time-value">{formatDuration(remainingSeconds)}</span>
    </div>
  );
}
