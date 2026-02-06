// Project Dashboard - Main view for working with a project's structures and calculations

import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { parseCIF } from "../lib/cifParser";
import { CrystalData } from "../lib/types";

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
  last_opened_cif_id: string | null;
}

interface ProjectDashboardProps {
  projectId: string;
  onBack: () => void;
  onDeleted: () => void;
  onRunSCF: (cifId: string, crystalData: CrystalData, cifContent: string, filename: string) => void;
}

const CONFIRM_TEXT = "delete my project for good";

export function ProjectDashboard({
  projectId,
  onBack,
  onDeleted,
  onRunSCF,
}: ProjectDashboardProps) {
  const [project, setProject] = useState<Project | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // Current CIF selection
  const [selectedCifId, setSelectedCifId] = useState<string | null>(null);
  const [crystalData, setCrystalData] = useState<CrystalData | null>(null);
  const [cifContent, setCifContent] = useState<string>("");

  // Settings menu state
  const [showSettingsMenu, setShowSettingsMenu] = useState(false);
  const settingsRef = useRef<HTMLDivElement>(null);

  // Delete confirmation dialog state
  const [showDeleteDialog, setShowDeleteDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeleting, setIsDeleting] = useState(false);

  // Import state
  const [isImporting, setIsImporting] = useState(false);

  // Expanded calculation
  const [expandedCalc, setExpandedCalc] = useState<string | null>(null);

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

      // Determine which CIF to show
      if (proj.cif_variants.length > 0) {
        const cifToOpen = proj.last_opened_cif_id &&
          proj.cif_variants.some(v => v.id === proj.last_opened_cif_id)
          ? proj.last_opened_cif_id
          : proj.cif_variants[0].id;

        await selectCif(cifToOpen);
      }
    } catch (e) {
      console.error("Failed to load project:", e);
      setError(String(e));
    } finally {
      setIsLoading(false);
    }
  }

  async function selectCif(cifId: string) {
    setSelectedCifId(cifId);
    try {
      // Load crystal data
      const data = await invoke<CrystalData>("get_cif_crystal_data", {
        projectId,
        cifId,
      });
      setCrystalData(data);

      // Load CIF content
      const content = await invoke<string>("get_cif_content", {
        projectId,
        cifId,
      });
      setCifContent(content);

      // Update last opened
      await invoke("set_last_opened_cif", { projectId, cifId });
    } catch (e) {
      console.error("Failed to load CIF data:", e);
      setError(`Failed to load structure: ${e}`);
    }
  }

  async function handleImportCIF() {
    try {
      const selected = await open({
        multiple: false,
        filters: [{ name: "CIF Files", extensions: ["cif"] }],
        title: "Select CIF File",
      });

      if (selected && typeof selected === "string") {
        setIsImporting(true);
        setError(null);

        const content = await readTextFile(selected);
        const parsedData = parseCIF(content);
        const filename = selected.split("/").pop() || "structure.cif";
        const formula = parsedData.formula_sum || parsedData.formula_structural || "Unknown";

        const newVariant = await invoke<CifVariant>("add_cif_to_project", {
          projectId,
          cifData: {
            filename,
            formula,
            content,
            crystal_data: parsedData,
          },
        });

        // Reload project and select the new CIF
        await loadProject();
        await selectCif(newVariant.id);
      }
    } catch (e) {
      console.error("Failed to import CIF:", e);
      setError(`Failed to import CIF: ${e}`);
    } finally {
      setIsImporting(false);
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

  function handleRunSCF() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunSCF(selectedCifId, crystalData, cifContent, variant.filename);
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

  function getSelectedVariant(): CifVariant | undefined {
    return project?.cif_variants.find(v => v.id === selectedCifId);
  }

  if (isLoading) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Loading...</h2>
        </div>
      </div>
    );
  }

  if (error && !project) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Error</h2>
        </div>
        <div className="error-banner">{error}</div>
      </div>
    );
  }

  if (!project) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <h2>Project not found</h2>
        </div>
      </div>
    );
  }

  // Empty state - no CIF files yet
  if (project.cif_variants.length === 0) {
    return (
      <div className="dashboard-container">
        <div className="dashboard-header">
          <button className="back-btn" onClick={onBack}>
            ← Back
          </button>
          <div className="dashboard-title">
            <h2>{project.name}</h2>
            {project.description && (
              <p className="dashboard-description">{project.description}</p>
            )}
          </div>
        </div>

        {error && <div className="error-banner">{error}</div>}

        <div className="dashboard-content">
          <div className="empty-state">
            <div className="empty-icon">📄</div>
            <h3>No Structure File</h3>
            <p>Import a CIF file to get started with your calculations</p>
            <button
              className="add-structure-btn primary"
              onClick={handleImportCIF}
              disabled={isImporting}
            >
              {isImporting ? "Importing..." : "Import CIF File"}
            </button>
          </div>
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

        {/* Delete Dialog */}
        {showDeleteDialog && renderDeleteDialog()}
      </div>
    );
  }

  const selectedVariant = getSelectedVariant();

  return (
    <div className="dashboard-container">
      <div className="dashboard-header">
        <button className="back-btn" onClick={onBack}>
          ← Back
        </button>
        <div className="dashboard-title">
          <h2>{project.name}</h2>
        </div>
        <div className="structure-selector">
          <select
            value={selectedCifId || ""}
            onChange={(e) => selectCif(e.target.value)}
          >
            {project.cif_variants.map((variant) => (
              <option key={variant.id} value={variant.id}>
                {variant.formula} ({variant.filename})
              </option>
            ))}
          </select>
          <button
            className="add-structure-inline-btn"
            onClick={handleImportCIF}
            disabled={isImporting}
            title="Add another structure"
          >
            +
          </button>
        </div>
      </div>

      {error && <div className="error-banner">{error}</div>}

      <div className="dashboard-content">
        {/* Structure Info Hero */}
        {crystalData && (
          <section className="structure-hero">
            <div className="hero-formula">
              {crystalData.formula_sum || crystalData.formula_structural || "Unknown"}
            </div>
            <div className="hero-details">
              <div className="hero-detail-item">
                <label>Space Group</label>
                <span>
                  {crystalData.space_group_HM || "N/A"}
                  {crystalData.space_group_IT_number && ` (#${crystalData.space_group_IT_number})`}
                </span>
              </div>
              <div className="hero-detail-item">
                <label>Lattice Parameters</label>
                <span>
                  a = {crystalData.cell_length_a.value.toFixed(4)} A,{" "}
                  b = {crystalData.cell_length_b.value.toFixed(4)} A,{" "}
                  c = {crystalData.cell_length_c.value.toFixed(4)} A
                </span>
              </div>
              <div className="hero-detail-item">
                <label>Angles</label>
                <span>
                  alpha = {crystalData.cell_angle_alpha.value.toFixed(2)} deg,{" "}
                  beta = {crystalData.cell_angle_beta.value.toFixed(2)} deg,{" "}
                  gamma = {crystalData.cell_angle_gamma.value.toFixed(2)} deg
                </span>
              </div>
              {crystalData.cell_volume && (
                <div className="hero-detail-item">
                  <label>Volume</label>
                  <span>{crystalData.cell_volume.toFixed(2)} A^3</span>
                </div>
              )}
              <div className="hero-detail-item">
                <label>Atoms</label>
                <span>{crystalData.atom_sites.length} sites</span>
              </div>
            </div>
          </section>
        )}

        {/* Calculation Actions */}
        <section className="actions-section">
          <h3>Calculations</h3>
          <div className="calc-action-grid">
            <button className="calc-action-btn" onClick={handleRunSCF}>
              <span className="calc-action-icon">SCF</span>
              <span className="calc-action-label">Self-Consistent Field</span>
              <span className="calc-action-hint">Ground state energy</span>
            </button>
            <button className="calc-action-btn" disabled>
              <span className="calc-action-icon">Band</span>
              <span className="calc-action-label">Band Structure</span>
              <span className="calc-action-hint">Coming soon</span>
            </button>
            <button className="calc-action-btn" disabled>
              <span className="calc-action-icon">DOS</span>
              <span className="calc-action-label">Density of States</span>
              <span className="calc-action-hint">Coming soon</span>
            </button>
            <button className="calc-action-btn" disabled>
              <span className="calc-action-icon">Opt</span>
              <span className="calc-action-label">Geometry Optimization</span>
              <span className="calc-action-hint">Coming soon</span>
            </button>
          </div>
        </section>

        {/* Calculation History */}
        {selectedVariant && selectedVariant.calculations.length > 0 && (
          <section className="history-section">
            <h3>Calculation History</h3>
            <div className="calculations-list">
              {selectedVariant.calculations.map((calc) => (
                <div key={calc.id} className="calculation-item">
                  <div
                    className="calculation-header"
                    onClick={() =>
                      setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                    }
                  >
                    <div className="calculation-info">
                      <span className="calc-type">{calc.calc_type.toUpperCase()}</span>
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
                        {expandedCalc === calc.id ? "v" : ">"}
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
          </section>
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

      {/* Delete Dialog */}
      {showDeleteDialog && renderDeleteDialog()}
    </div>
  );

  function renderDeleteDialog() {
    return (
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
                You are about to permanently delete <strong>{project?.name}</strong> and all of its data:
              </p>
              <ul>
                <li>{project?.cif_variants.length} structure{project?.cif_variants.length !== 1 ? "s" : ""}</li>
                <li>
                  {project?.cif_variants.reduce((sum, v) => sum + v.calculations.length, 0)} calculation
                  {project?.cif_variants.reduce((sum, v) => sum + v.calculations.length, 0) !== 1 ? "s" : ""}
                </li>
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
    );
  }
}
