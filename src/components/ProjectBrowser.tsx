// Project Browser - Grid view of all projects

import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { open, save } from "@tauri-apps/plugin-dialog";
import { CreateProjectDialog } from "./CreateProjectDialog";

interface ProjectSummary {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  formula: string | null;
  calculation_count: number;
  last_activity: string;
}

interface ProjectBrowserProps {
  onBack: () => void;
  onSelectProject?: (projectId: string) => void;
  onProjectsChanged?: () => void;
}

interface ProjectArchiveExportResult {
  archive_path: string;
  project_id: string;
  project_name: string;
  archive_size_bytes: number;
}

interface ProjectArchiveImportResult {
  project_id: string;
  project_name: string;
  imported_with_new_id: boolean;
}

interface ProjectArchiveExportProgress {
  export_id: string;
  project_id: string;
  phase: "compressing" | "done" | string;
  processed_bytes: number;
  total_bytes: number;
  progress_percent: number;
}

interface ProjectArchiveImportProgress {
  import_id: string;
  phase: "extracting" | "installing" | "done" | string;
  processed_bytes: number;
  total_bytes: number;
  progress_percent: number;
}

export function ProjectBrowser({
  onBack,
  onSelectProject,
  onProjectsChanged,
}: ProjectBrowserProps) {
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [statusMessage, setStatusMessage] = useState<string | null>(null);
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [showExportDialog, setShowExportDialog] = useState(false);
  const [selectedExportProjectId, setSelectedExportProjectId] = useState<string>("");
  const [isExporting, setIsExporting] = useState(false);
  const [isCancelingExport, setIsCancelingExport] = useState(false);
  const [isImporting, setIsImporting] = useState(false);
  const [exportProgress, setExportProgress] = useState<ProjectArchiveExportProgress | null>(null);
  const [importProgress, setImportProgress] = useState<ProjectArchiveImportProgress | null>(null);
  const activeExportIdRef = useRef<string | null>(null);
  const activeImportIdRef = useRef<string | null>(null);

  useEffect(() => {
    loadProjects();
  }, []);

  useEffect(() => {
    const unlistenPromise = listen<ProjectArchiveExportProgress>("project-export-progress", (event) => {
      const payload = event.payload;
      if (!activeExportIdRef.current || payload.export_id !== activeExportIdRef.current) {
        return;
      }
      setExportProgress(payload);
    });

    return () => {
      void unlistenPromise.then((unlisten) => unlisten());
    };
  }, []);

  useEffect(() => {
    const unlistenPromise = listen<ProjectArchiveImportProgress>("project-import-progress", (event) => {
      const payload = event.payload;
      if (!activeImportIdRef.current || payload.import_id !== activeImportIdRef.current) {
        return;
      }
      setImportProgress(payload);
    });

    return () => {
      void unlistenPromise.then((unlisten) => unlisten());
    };
  }, []);

  async function loadProjects() {
    setIsLoading(true);
    setError(null);
    try {
      const projectList = await invoke<ProjectSummary[]>("list_projects");
      setProjects(projectList);
      setSelectedExportProjectId((currentSelection) => {
        if (projectList.length === 0) return "";
        if (currentSelection && projectList.some((project) => project.id === currentSelection)) {
          return currentSelection;
        }
        return projectList[0].id;
      });
    } catch (e) {
      console.error("Failed to load projects:", e);
      setError(String(e));
    } finally {
      setIsLoading(false);
    }
  }

  function formatDate(isoString: string): string {
    try {
      const date = new Date(isoString);
      return date.toLocaleDateString(undefined, {
        year: "numeric",
        month: "short",
        day: "numeric",
      });
    } catch {
      return isoString;
    }
  }

  function formatRelativeTime(isoString: string): string {
    try {
      const date = new Date(isoString);
      const now = new Date();
      const diffMs = now.getTime() - date.getTime();
      const diffDays = Math.floor(diffMs / (1000 * 60 * 60 * 24));

      if (diffDays === 0) return "Today";
      if (diffDays === 1) return "Yesterday";
      if (diffDays < 7) return `${diffDays} days ago`;
      if (diffDays < 30) return `${Math.floor(diffDays / 7)} weeks ago`;
      if (diffDays < 365) return `${Math.floor(diffDays / 30)} months ago`;
      return `${Math.floor(diffDays / 365)} years ago`;
    } catch {
      return isoString;
    }
  }

  function formatBytes(bytes: number): string {
    if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
    const units = ["B", "KB", "MB", "GB", "TB"];
    let size = bytes;
    let unitIndex = 0;
    while (size >= 1024 && unitIndex < units.length - 1) {
      size /= 1024;
      unitIndex += 1;
    }
    const decimals = size >= 10 || unitIndex === 0 ? 0 : 1;
    return `${size.toFixed(decimals)} ${units[unitIndex]}`;
  }

  function getProgressText(progress: ProjectArchiveExportProgress | null): string {
    if (!progress) return "Preparing export...";
    if (progress.phase === "done") return "Export complete";
    if (progress.phase === "scanning") return "Scanning project files...";
    if (progress.phase === "finalizing") return "Finalizing archive...";
    if (progress.total_bytes <= 0) return "Compressing project files...";
    return `Compressing ${Math.round(progress.progress_percent)}%`;
  }

  function getImportProgressText(progress: ProjectArchiveImportProgress | null): string {
    if (!progress) return "Preparing import...";
    if (progress.phase === "done") return "Import complete";
    if (progress.phase === "extracting") {
      if (progress.total_bytes <= 0) return "Extracting archive...";
      return `Extracting archive ${Math.round(progress.progress_percent)}%`;
    }
    if (progress.phase === "installing") {
      if (progress.total_bytes <= 0) return "Installing project files...";
      return `Installing files ${Math.round(progress.progress_percent)}%`;
    }
    return "Importing project...";
  }

  function sanitizeFilenamePart(value: string): string {
    const cleaned = value
      .trim()
      .replace(/[^a-zA-Z0-9._-]+/g, "-")
      .replace(/-+/g, "-")
      .replace(/^-|-$/g, "");
    return cleaned || "project";
  }

  function openExportDialog() {
    if (projects.length === 0) {
      setError("No projects are available to export.");
      return;
    }

    setError(null);
    setStatusMessage(null);
    setSelectedExportProjectId((currentSelection) => {
      if (currentSelection && projects.some((project) => project.id === currentSelection)) {
        return currentSelection;
      }
      return projects[0].id;
    });
    setExportProgress(null);
    setShowExportDialog(true);
  }

  function closeExportDialog() {
    if (isExporting) return;
    setExportProgress(null);
    setShowExportDialog(false);
  }

  function isExportCanceledError(error: unknown): boolean {
    return String(error).toLowerCase().includes("export canceled by user");
  }

  async function handleCancelExport() {
    if (!isExporting) {
      closeExportDialog();
      return;
    }

    const activeExportId = activeExportIdRef.current;
    if (!activeExportId) {
      closeExportDialog();
      return;
    }

    setIsCancelingExport(true);
    try {
      await invoke<boolean>("cancel_project_export", {
        exportId: activeExportId,
      });
    } catch (e) {
      console.error("Failed to cancel export:", e);
      setError(`Failed to cancel export: ${e}`);
      setIsCancelingExport(false);
    }
  }

  async function handleExportProject() {
    if (!selectedExportProjectId) {
      setError("Select a project to export.");
      return;
    }

    const selectedProject = projects.find((project) => project.id === selectedExportProjectId);
    if (!selectedProject) {
      setError("Selected project could not be found.");
      return;
    }

    setIsExporting(true);
    setIsCancelingExport(false);
    setError(null);
    setStatusMessage(null);
    setExportProgress(null);

    try {
      const defaultArchiveName = `${sanitizeFilenamePart(selectedProject.name)}-${selectedProject.id}.qcproj`;
      const destinationPath = await save({
        title: "Export Project Archive",
        defaultPath: defaultArchiveName,
        filters: [{ name: "QCortado Project Archive", extensions: ["qcproj"] }],
      });

      if (!destinationPath) {
        return;
      }

      const exportId = `export_${Date.now()}_${Math.random().toString(16).slice(2)}`;
      activeExportIdRef.current = exportId;
      const result = await invoke<ProjectArchiveExportResult>("export_project_archive", {
        projectId: selectedExportProjectId,
        destinationPath,
        exportId,
      });
      setShowExportDialog(false);
      setExportProgress(null);
      setStatusMessage(
        `Exported "${result.project_name}" (${formatBytes(result.archive_size_bytes)}) to ${result.archive_path}`,
      );
    } catch (e) {
      if (isExportCanceledError(e)) {
        setShowExportDialog(false);
        setExportProgress(null);
        setStatusMessage("Export canceled.");
      } else {
        console.error("Failed to export project:", e);
        setError(String(e));
      }
    } finally {
      activeExportIdRef.current = null;
      setIsExporting(false);
      setIsCancelingExport(false);
    }
  }

  async function handleImportProject() {
    if (isImporting || isExporting) {
      return;
    }

    setError(null);
    setStatusMessage(null);

    try {
      const selectedArchivePath = await open({
        title: "Import Project Archive",
        directory: false,
        multiple: false,
        filters: [{ name: "QCortado Project Archive", extensions: ["qcproj"] }],
      });

      if (!selectedArchivePath || Array.isArray(selectedArchivePath)) {
        return;
      }

      const importId = `import_${Date.now()}_${Math.random().toString(16).slice(2)}`;
      activeImportIdRef.current = importId;
      setImportProgress(null);
      setIsImporting(true);

      const result = await invoke<ProjectArchiveImportResult>("import_project_archive", {
        archivePath: selectedArchivePath,
        importId,
      });
      await loadProjects();
      onProjectsChanged?.();
      setStatusMessage(
        result.imported_with_new_id
          ? `Imported "${result.project_name}" with a new project ID to avoid a conflict.`
          : `Imported "${result.project_name}".`,
      );
    } catch (e) {
      console.error("Failed to import project archive:", e);
      setError(String(e));
    } finally {
      activeImportIdRef.current = null;
      setImportProgress(null);
      setIsImporting(false);
    }
  }

  return (
    <div className="browser-container">
      <div className="browser-header">
        <button className="back-btn" onClick={onBack}>
          ← Back
        </button>
        <h2>Projects</h2>
        <div className="browser-actions">
          <button
            className="secondary-project-btn"
            onClick={handleImportProject}
            disabled={isImporting || isExporting}
          >
            {isImporting ? "Importing..." : "Import"}
          </button>
          <button
            className="secondary-project-btn"
            onClick={openExportDialog}
            disabled={projects.length === 0 || isImporting || isExporting}
          >
            Export
          </button>
          <button className="new-project-btn" onClick={() => setShowCreateDialog(true)}>
            + New Project
          </button>
        </div>
      </div>

      {error && <div className="error-banner">{error}</div>}
      {statusMessage && <div className="success-banner">{statusMessage}</div>}

      <div className="browser-content">
        {isLoading ? (
          <div className="loading-state">Loading projects...</div>
        ) : projects.length === 0 ? (
          <div className="empty-state">
            <div className="empty-icon">📁</div>
            <h3>No Projects Yet</h3>
            <p>Create a project to organize your calculations</p>
            <button className="new-project-btn" onClick={() => setShowCreateDialog(true)}>
              Create Your First Project
            </button>
          </div>
        ) : (
          <div className="project-grid">
            {projects.map((project) => (
              <div
                key={project.id}
                className="project-card"
                onClick={() => onSelectProject?.(project.id)}
              >
                <div className="project-card-header">
                  <h3 className="project-name">{project.name}</h3>
                </div>

                {project.description && (
                  <p className="project-description">{project.description}</p>
                )}

                <div className="project-meta">
                  {project.formula && (
                    <span className="project-formula">{project.formula}</span>
                  )}
                  <span className="project-stats">
                    {project.calculation_count} calculation
                    {project.calculation_count !== 1 ? "s" : ""}
                  </span>
                </div>

                <div className="project-footer">
                  <span className="project-date" title={formatDate(project.created_at)}>
                    Created {formatDate(project.created_at)}
                  </span>
                  <span className="project-activity" title={formatDate(project.last_activity)}>
                    Active {formatRelativeTime(project.last_activity)}
                  </span>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>

      <CreateProjectDialog
        isOpen={showCreateDialog}
        onClose={() => setShowCreateDialog(false)}
        onCreated={() => {
          setShowCreateDialog(false);
          void loadProjects();
          onProjectsChanged?.();
        }}
      />

      {showExportDialog && (
        <div className="dialog-overlay" onClick={closeExportDialog}>
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>Export Project</h2>
              <button
                className="dialog-close"
                onClick={() => {
                  void handleCancelExport();
                }}
                disabled={isCancelingExport}
              >
                &times;
              </button>
            </div>
            <div className="dialog-body">
              <div className="save-form">
                <div className="form-group">
                  <label>Project</label>
                  <select
                    value={selectedExportProjectId}
                    onChange={(e) => setSelectedExportProjectId(e.target.value)}
                    disabled={isExporting}
                  >
                    {projects.map((project) => (
                      <option key={project.id} value={project.id}>
                        {project.name}
                      </option>
                    ))}
                  </select>
                </div>
                <p className="project-archive-hint">
                  The full project directory will be compressed into a `.qcproj` archive.
                </p>
                {isExporting && (
                  <div className="export-progress">
                    <div className="export-progress-labels">
                      <span>{getProgressText(exportProgress)}</span>
                      {exportProgress && exportProgress.total_bytes > 0 && (
                        <span>
                          {formatBytes(exportProgress.processed_bytes)} / {formatBytes(exportProgress.total_bytes)}
                        </span>
                      )}
                    </div>
                    <div
                      className={`export-progress-track ${
                        exportProgress && exportProgress.total_bytes > 0 ? "" : "indeterminate"
                      }`}
                    >
                      <div
                        className="export-progress-fill"
                        style={{
                          width:
                            exportProgress && exportProgress.total_bytes > 0
                              ? `${Math.min(100, Math.max(0, exportProgress.progress_percent))}%`
                              : "35%",
                        }}
                      />
                    </div>
                  </div>
                )}
              </div>
            </div>
            <div className="dialog-footer">
              <button
                className="dialog-btn cancel"
                onClick={() => {
                  void handleCancelExport();
                }}
                disabled={isCancelingExport}
              >
                {isExporting ? (isCancelingExport ? "Canceling..." : "Cancel Export") : "Cancel"}
              </button>
              <button
                className="dialog-btn save"
                onClick={() => {
                  void handleExportProject();
                }}
                disabled={isExporting || isCancelingExport || !selectedExportProjectId}
              >
                {isExporting ? "Exporting..." : "Export Project"}
              </button>
            </div>
          </div>
        </div>
      )}

      {isImporting && (
        <div className="dialog-overlay">
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>Import Project</h2>
            </div>
            <div className="dialog-body">
              <div className="save-form">
                <p className="project-archive-hint">
                  The archive is being extracted and copied into your projects directory.
                </p>
                <div className="export-progress">
                  <div className="export-progress-labels">
                    <span>{getImportProgressText(importProgress)}</span>
                    {importProgress && importProgress.total_bytes > 0 && (
                      <span>
                        {formatBytes(importProgress.processed_bytes)} / {formatBytes(importProgress.total_bytes)}
                      </span>
                    )}
                  </div>
                  <div
                    className={`export-progress-track ${
                      importProgress && importProgress.total_bytes > 0 ? "" : "indeterminate"
                    }`}
                  >
                    <div
                      className="export-progress-fill"
                      style={{
                        width:
                          importProgress && importProgress.total_bytes > 0
                            ? `${Math.min(100, Math.max(0, importProgress.progress_percent))}%`
                            : "35%",
                      }}
                    />
                  </div>
                </div>
              </div>
            </div>
            <div className="dialog-footer">
              <button className="dialog-btn save" disabled>
                Importing...
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
