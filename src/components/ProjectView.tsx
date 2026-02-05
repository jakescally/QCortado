// Project View - Display project details, CIF variants, and calculations

import { useState, useEffect } from "react";
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
}

export function ProjectView({ projectId, onBack }: ProjectViewProps) {
  const [project, setProject] = useState<Project | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [expandedCalc, setExpandedCalc] = useState<string | null>(null);

  useEffect(() => {
    loadProject();
  }, [projectId]);

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
    </div>
  );
}
