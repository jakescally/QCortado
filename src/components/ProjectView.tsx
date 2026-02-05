// Project View - Display project details, CIF variants, and calculations

import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";

interface QEResult {
  converged: boolean;
  total_energy: number | null;
  fermi_energy: number | null;
  n_scf_steps: number | null;
  wall_time_seconds: number | null;
  raw_output: string;
}

interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: unknown;
  result: QEResult | null;
  started_at: string;
  completed_at: string | null;
}

interface CifVariant {
  id: string;
  filename: string;
  formula: string;
  added_at: string;
  calculations: CalculationRun[];
}

interface Project {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  cif_variants: CifVariant[];
}

interface ProjectViewProps {
  projectId: string;
  onBack: () => void;
  onDeleted: () => void;
}

const CONFIRM_TEXT = "delete my project for good";

export function ProjectView({ projectId, onBack, onDeleted }: ProjectViewProps) {
  const [project, setProject] = useState<Project | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [expandedCalc, setExpandedCalc] = useState<string | null>(null);

  // Settings menu state
  const [showSettingsMenu, setShowSettingsMenu] = useState(false);
  const settingsRef = useRef<HTMLDivElement>(null);

  // Delete confirmation dialog state
  const [showDeleteDialog, setShowDeleteDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeleting, setIsDeleting] = useState(false);

  useEffect(() => {
    loadProject();
  }, [projectId]);

  // Close settings menu when clicking outside
  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (settingsRef.current && !settingsRef.current.contains(event.target as Node)) {
        setShowSettingsMenu(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  async function loadProject() {
    setIsLoading(true);
    setError(null);
    try {
      const proj = await invoke<Project>("get_project", { projectId });
      setProject(proj);
    } catch (e) {
      console.error("Failed to load project:", e);
      setError(String(e));
    } finally {
      setIsLoading(false);
    }
  }

  function openDeleteDialog() {
    setShowSettingsMenu(false);
    setDeleteConfirmText("");
    setShowDeleteDialog(true);
  }

  async function handleConfirmDelete() {
    if (deleteConfirmText !== CONFIRM_TEXT) return;

    setIsDeleting(true);
    try {
      await invoke("delete_project", { projectId });
      onDeleted();
    } catch (e) {
      console.error("Failed to delete project:", e);
      setError(String(e));
      setIsDeleting(false);
      setShowDeleteDialog(false);
    }
  }

  function formatDate(isoString: string): string {
    try {
      const date = new Date(isoString);
      return date.toLocaleDateString(undefined, {
        year: "numeric",
        month: "short",
        day: "numeric",
        hour: "2-digit",
        minute: "2-digit",
      });
    } catch {
      return isoString;
    }
  }

  function formatEnergy(energy: number): string {
    return `${energy.toFixed(6)} Ry`;
  }

  function getTotalCalculations(): number {
    if (!project) return 0;
    return project.cif_variants.reduce((sum, v) => sum + v.calculations.length, 0);
  }

  if (isLoading) {
    return (
      <div className="project-view-container">
        <div className="project-view-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Loading...</h2>
        </div>
      </div>
    );
  }

  if (error || !project) {
    return (
      <div className="project-view-container">
        <div className="project-view-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Error</h2>
        </div>
        <div className="error-banner">{error || "Project not found"}</div>
      </div>
    );
  }

  return (
    <div className="project-view-container">
      <div className="project-view-header">
        <button className="back-btn" onClick={onBack}>
          ← Back
        </button>
        <div className="project-view-title">
          <h2>{project.name}</h2>
          {project.description && (
            <p className="project-view-description">{project.description}</p>
          )}
        </div>
        <div className="project-view-meta">
          <span className="meta-item">
            Created {formatDate(project.created_at)}
          </span>
          <span className="meta-item">
            {project.cif_variants.length} structure{project.cif_variants.length !== 1 ? "s" : ""}
          </span>
          <span className="meta-item">
            {getTotalCalculations()} calculation{getTotalCalculations() !== 1 ? "s" : ""}
          </span>
        </div>
      </div>

      <div className="project-view-content">
        {project.cif_variants.length === 0 ? (
          <div className="empty-state">
            <div className="empty-icon">📄</div>
            <h3>No Structures Yet</h3>
            <p>Run a calculation and save it to this project to add structures</p>
          </div>
        ) : (
          <div className="cif-variants-list">
            {project.cif_variants.map((variant) => (
              <div key={variant.id} className="cif-variant-card">
                <div className="cif-variant-header">
                  <div className="cif-variant-info">
                    <span className="cif-formula">{variant.formula}</span>
                    <span className="cif-filename">{variant.filename}</span>
                  </div>
                  <span className="cif-added">
                    Added {formatDate(variant.added_at)}
                  </span>
                </div>

                {variant.calculations.length === 0 ? (
                  <div className="no-calculations">
                    No calculations yet
                  </div>
                ) : (
                  <div className="calculations-list">
                    {variant.calculations.map((calc) => (
                      <div key={calc.id} className="calculation-item">
                        <div
                          className="calculation-header"
                          onClick={() =>
                            setExpandedCalc(
                              expandedCalc === calc.id ? null : calc.id
                            )
                          }
                        >
                          <div className="calculation-info">
                            <span className="calc-type">
                              {calc.calc_type.toUpperCase()}
                            </span>
                            {calc.result && (
                              <span
                                className={`calc-status ${
                                  calc.result.converged ? "converged" : "failed"
                                }`}
                              >
                                {calc.result.converged ? "Converged" : "Not converged"}
                              </span>
                            )}
                            {calc.result?.total_energy && (
                              <span className="calc-energy">
                                E = {formatEnergy(calc.result.total_energy)}
                              </span>
                            )}
                          </div>
                          <div className="calculation-meta">
                            <span className="calc-date">
                              {calc.completed_at
                                ? formatDate(calc.completed_at)
                                : "In progress..."}
                            </span>
                            <span className="expand-icon">
                              {expandedCalc === calc.id ? "▼" : "▶"}
                            </span>
                          </div>
                        </div>

                        {expandedCalc === calc.id && calc.result && (
                          <div className="calculation-details">
                            <div className="details-grid">
                              {calc.result.total_energy && (
                                <div className="detail-item">
                                  <label>Total Energy</label>
                                  <span>{formatEnergy(calc.result.total_energy)}</span>
                                </div>
                              )}
                              {calc.result.fermi_energy && (
                                <div className="detail-item">
                                  <label>Fermi Energy</label>
                                  <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                                </div>
                              )}
                              {calc.result.n_scf_steps && (
                                <div className="detail-item">
                                  <label>SCF Steps</label>
                                  <span>{calc.result.n_scf_steps}</span>
                                </div>
                              )}
                              {calc.result.wall_time_seconds && (
                                <div className="detail-item">
                                  <label>Wall Time</label>
                                  <span>{calc.result.wall_time_seconds.toFixed(1)} s</span>
                                </div>
                              )}
                            </div>
                            <div className="detail-item parameters">
                              <label>Parameters</label>
                              <pre>{JSON.stringify(calc.parameters, null, 2)}</pre>
                            </div>
                          </div>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Floating Settings Button */}
      <div className="floating-settings" ref={settingsRef}>
        <button
          className="floating-settings-btn"
          onClick={() => setShowSettingsMenu(!showSettingsMenu)}
          title="Project settings"
        >
          <svg width="24" height="24" viewBox="0 0 20 20" fill="currentColor">
            <path fillRule="evenodd" d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.061 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.533 1.533 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.533 1.533 0 01.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd" />
          </svg>
        </button>

        {showSettingsMenu && (
          <div className="floating-settings-menu">
            <button className="settings-menu-item danger" onClick={openDeleteDialog}>
              Delete Project
            </button>
          </div>
        )}
      </div>

      {/* Delete Confirmation Dialog */}
      {showDeleteDialog && (
        <div className="dialog-overlay" onClick={() => !isDeleting && setShowDeleteDialog(false)}>
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>Delete Project</h2>
              <button
                className="dialog-close"
                onClick={() => setShowDeleteDialog(false)}
                disabled={isDeleting}
              >
                &times;
              </button>
            </div>

            <div className="dialog-body">
              <div className="delete-warning">
                <p>
                  You are about to permanently delete <strong>{project.name}</strong> and all of its data:
                </p>
                <ul>
                  <li>{project.cif_variants.length} structure{project.cif_variants.length !== 1 ? "s" : ""}</li>
                  <li>{getTotalCalculations()} calculation{getTotalCalculations() !== 1 ? "s" : ""}</li>
                  <li>All input/output files</li>
                </ul>
                <p className="delete-warning-emphasis">
                  This action cannot be undone.
                </p>
              </div>

              <div className="form-group">
                <label>
                  Type <code>{CONFIRM_TEXT}</code> to confirm:
                </label>
                <input
                  type="text"
                  value={deleteConfirmText}
                  onChange={(e) => setDeleteConfirmText(e.target.value)}
                  placeholder={CONFIRM_TEXT}
                  disabled={isDeleting}
                  autoFocus
                />
              </div>
            </div>

            <div className="dialog-footer">
              <button
                className="dialog-btn cancel"
                onClick={() => setShowDeleteDialog(false)}
                disabled={isDeleting}
              >
                Cancel
              </button>
              <button
                className="dialog-btn delete"
                onClick={handleConfirmDelete}
                disabled={deleteConfirmText !== CONFIRM_TEXT || isDeleting}
              >
                {isDeleting ? "Deleting..." : "Delete Project"}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
