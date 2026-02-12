import { useEffect, useState } from "react";

interface ElapsedTimerProps {
  startedAt?: string | null;
  isRunning: boolean;
  label?: string;
}

function getElapsedSeconds(startedAt: string): number {
  const startMs = Date.parse(startedAt);
  if (Number.isNaN(startMs)) return 0;
  return Math.max(0, Math.floor((Date.now() - startMs) / 1000));
}

function formatElapsed(seconds: number): string {
  const hours = Math.floor(seconds / 3600);
  const minutes = Math.floor((seconds % 3600) / 60);
  const secs = seconds % 60;
  return `${String(hours).padStart(2, "0")}:${String(minutes).padStart(2, "0")}:${String(secs).padStart(2, "0")}`;
}

export function ElapsedTimer({ startedAt, isRunning, label = "Elapsed" }: ElapsedTimerProps) {
  const [elapsedSeconds, setElapsedSeconds] = useState(0);
  const hasStart = typeof startedAt === "string" && startedAt.trim().length > 0;

  useEffect(() => {
    if (!hasStart || !startedAt) {
      setElapsedSeconds(0);
      return;
    }
    setElapsedSeconds(getElapsedSeconds(startedAt));
  }, [hasStart, startedAt]);

  useEffect(() => {
    if (!hasStart || !startedAt) return;

    if (!isRunning) {
      setElapsedSeconds(getElapsedSeconds(startedAt));
      return;
    }

    const intervalId = window.setInterval(() => {
      setElapsedSeconds(getElapsedSeconds(startedAt));
    }, 1000);

    return () => {
      window.clearInterval(intervalId);
    };
  }, [hasStart, isRunning, startedAt]);

  if (!hasStart) return null;

  return (
    <div className="elapsed-time" aria-live="polite">
      <span className="elapsed-time-label">{label}:</span>
      <span className="elapsed-time-value">{formatElapsed(elapsedSeconds)}</span>
    </div>
  );
}
