import { useCallback, useEffect, useMemo, useState } from "react";
import { ArchiveRestore, RefreshCw, ShieldCheck } from "lucide-react";
import { AppHeaderPortal } from "./AppHeaderPortal";
import type { CalculationRun } from "./ProjectDashboard";
import type { HpcProfile } from "../lib/types";
import {
  listHpcCalculationBackups,
  restoreHpcCalculationBackup,
  type HpcCalculationBackupSummary,
} from "../lib/hpcConfig";

interface BackupManagerPageProps {
  onBack: () => void;
  onRestored: (projectId: string) => void;
  activeHpcProfile?: HpcProfile | null;
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

function formatDateTime(value: string): string {
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? value : date.toLocaleString();
}

function derivedBackupRoot(profile: HpcProfile | null | undefined): string {
  const configured = profile?.remote_backup_root?.trim();
  if (configured) return configured;
  const projectRoot = profile?.remote_project_root?.replace(/\/+$/, "");
  if (!projectRoot) return "Sibling backups directory";
  const slash = projectRoot.lastIndexOf("/");
  return slash >= 0 ? `${projectRoot.slice(0, slash)}/backups` : "backups";
}

export function BackupManagerPage({
  onBack,
  onRestored,
  activeHpcProfile = null,
}: BackupManagerPageProps) {
  const [backups, setBackups] = useState<HpcCalculationBackupSummary[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [restoringId, setRestoringId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [statusMessage, setStatusMessage] = useState<string | null>(null);

  const loadBackups = useCallback(async () => {
    if (!activeHpcProfile) {
      setBackups([]);
      setError("No active HPC profile is configured.");
      return;
    }
    setIsLoading(true);
    setError(null);
    try {
      setBackups(await listHpcCalculationBackups(activeHpcProfile.id));
    } catch (loadError) {
      console.error("Failed to list calculation backups:", loadError);
      setBackups([]);
      setError(String(loadError));
    } finally {
      setIsLoading(false);
    }
  }, [activeHpcProfile]);

  useEffect(() => {
    void loadBackups();
  }, [loadBackups]);

  const totalBytes = useMemo(
    () => backups.reduce((sum, backup) => sum + backup.bytes, 0),
    [backups],
  );

  async function restoreBackup(backup: HpcCalculationBackupSummary) {
    setRestoringId(backup.backup_id);
    setError(null);
    setStatusMessage(null);
    try {
      const restored = await restoreHpcCalculationBackup<CalculationRun>(
        backup.backup_id,
        activeHpcProfile?.id,
      );
      setStatusMessage(
        `Restored ${backup.calc_type.toUpperCase()} as a new calculation (${restored.id.slice(0, 8)}).`,
      );
      onRestored(backup.project_id);
    } catch (restoreError) {
      console.error("Failed to restore calculation backup:", restoreError);
      setError(String(restoreError));
    } finally {
      setRestoringId(null);
    }
  }

  return (
    <div className="queue-page-container storage-manager-page backup-manager-page">
      <AppHeaderPortal className="queue-page-header storage-manager-header">
        <div className="storage-manager-header-main">
          <button className="back-button" onClick={onBack}>← Back</button>
          <h2>Calculation Backups</h2>
        </div>
        <button
          className="dialog-btn cancel"
          type="button"
          onClick={() => void loadBackups()}
          disabled={isLoading || restoringId !== null}
        >
          <RefreshCw size={15} />
          {isLoading ? "Checking..." : "Refresh"}
        </button>
      </AppHeaderPortal>

      <div className="storage-manager-subhead">
        <span>Profile: {activeHpcProfile?.name ?? "None"}</span>
        <span>Protected remote root: {derivedBackupRoot(activeHpcProfile)}</span>
      </div>

      {error && <div className="storage-status error">{error}</div>}
      {statusMessage && <div className="storage-status success">{statusMessage}</div>}

      <div className="queue-page-summary storage-summary-grid backup-summary-grid">
        <div className="queue-summary-card">
          <span className="queue-summary-label">Backups</span>
          <span className="queue-summary-value">{backups.length}</span>
        </div>
        <div className="queue-summary-card">
          <span className="queue-summary-label">Protected Storage</span>
          <span className="queue-summary-value">{formatBytes(totalBytes)}</span>
        </div>
      </div>

      <section className="storage-section-card">
        <div className="storage-group-header">
          <div>
            <h3>Recoverable calculations</h3>
            <p className="settings-menu-hint">
              Restoring creates a new calculation entry and a new normal project copy; the backup remains unchanged.
            </p>
          </div>
          <ShieldCheck size={22} aria-label="Write protected" />
        </div>

        <div className="storage-entry-list">
          {backups.map((backup) => (
            <div key={backup.backup_id} className="queue-item storage-entry-row backup-entry-row">
              <ArchiveRestore size={19} />
              <div className="storage-entry-main">
                <div className="storage-entry-title-row">
                  <span className="storage-entry-title">{backup.project_name} • {backup.cif_label}</span>
                  <div className="calc-tags storage-entry-badges">
                    <span className="calc-tag calc-tag-info">{backup.calc_type.toUpperCase()}</span>
                    <span className="calc-tag calc-tag-hpc">Protected</span>
                  </div>
                </div>
                <div className="queue-item-meta storage-entry-meta">
                  <span>Created {formatDateTime(backup.created_at)}</span>
                  <span>{formatBytes(backup.bytes)}</span>
                  <span>{backup.remote_backup_path}</span>
                </div>
              </div>
              <button
                type="button"
                className="dialog-btn save"
                onClick={() => void restoreBackup(backup)}
                disabled={restoringId !== null}
              >
                {restoringId === backup.backup_id ? "Restoring..." : "Restore"}
              </button>
            </div>
          ))}
        </div>

        {!isLoading && backups.length === 0 && !error && (
          <div className="queue-empty-state">
            <h3>No backups yet</h3>
            <p>Open a project and use Create Backup on a saved HPC calculation.</p>
          </div>
        )}
      </section>
    </div>
  );
}
