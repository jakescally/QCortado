// Project Dashboard - Main view for working with a project's structures and calculations

import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { parseCIF } from "../lib/cifParser";
import { CrystalData, SCFPreset } from "../lib/types";

interface QEResult {
  converged: boolean;
  total_energy: number | null;
  fermi_energy: number | null;
  n_scf_steps: number | null;
  wall_time_seconds: number | null;
  raw_output: string;
  band_data?: any;  // Band structure data for bands calculations
}

export interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: QEResult | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
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
  onRunSCF: (
    cifId: string,
    crystalData: CrystalData,
    cifContent: string,
    filename: string,
    preset?: SCFPreset,
    presetLock?: boolean,
  ) => void;
  onRunBands: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewBands: (bandData: any, scfFermiEnergy: number | null) => void;
  onRunPhonons: (cifId: string, crystalData: CrystalData, scfCalculations: CalculationRun[]) => void;
  onViewPhonons: (phononData: any) => void;
}

const CONFIRM_TEXT = "delete my project for good";

// Helper to generate calculation feature tags from parameters
function getCalculationTags(calc: CalculationRun): { label: string; type: "info" | "feature" | "special" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" }[] = [];
  const params = calc.parameters || {};

  // Special tags from stored tags array (phonon-ready, structure-optimized)
  if (calc.tags) {
    for (const tag of calc.tags) {
      if (tag === "phonon-ready") {
        tags.push({ label: "Phonon-Ready", type: "special" });
      } else if (tag === "structure-optimized") {
        tags.push({ label: "Optimized", type: "special" });
      }
    }
  }

  // K-points grid
  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    tags.push({ label: `${k1}×${k2}×${k3}`, type: "info" });
  }

  // Convergence threshold
  if (params.conv_thr) {
    const thr = params.conv_thr;
    // Format as scientific notation if small
    const label = thr < 0.001 ? thr.toExponential(0) : thr.toString();
    tags.push({ label: label, type: "info" });
  }

  // Feature tags
  if (params.lspinorb) {
    tags.push({ label: "SOC", type: "feature" });
  }

  if (params.nspin === 4) {
    tags.push({ label: "Non-collinear", type: "feature" });
  } else if (params.nspin === 2) {
    tags.push({ label: "Magnetic", type: "feature" });
  }

  if (params.lda_plus_u) {
    tags.push({ label: "DFT+U", type: "feature" });
  }

  if (params.vdw_corr && params.vdw_corr !== "none") {
    tags.push({ label: "vdW", type: "feature" });
  }

  return tags;
}

// Helper to generate bands-specific tags
function getBandsTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};

  // K-points info
  if (params.total_k_points) {
    tags.push({ label: `${params.total_k_points} k-pts`, type: "info" });
  }

  // Inherited feature tags from SCF
  if (params.lspinorb) {
    tags.push({ label: "SOC", type: "feature" });
  }
  if (params.nspin === 2) {
    tags.push({ label: "Magnetic", type: "feature" });
  }
  if (params.lda_plus_u) {
    tags.push({ label: "DFT+U", type: "feature" });
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    tags.push({ label: "vdW", type: "feature" });
  }

  return tags;
}

// Helper to generate phonon-specific tags
function getPhononTags(calc: CalculationRun): { label: string; type: "info" | "feature" }[] {
  const tags: { label: string; type: "info" | "feature" }[] = [];
  const params = calc.parameters || {};

  if (params.q_grid) {
    const [q1, q2, q3] = params.q_grid;
    tags.push({ label: `${q1}×${q2}×${q3} Q`, type: "info" });
  }

  if (params.n_modes) {
    tags.push({ label: `${params.n_modes} modes`, type: "info" });
  }

  if (params.calculate_dos) {
    tags.push({ label: "DOS", type: "feature" });
  }

  if (params.calculate_dispersion) {
    tags.push({ label: "Dispersion", type: "feature" });
  }

  return tags;
}

export function ProjectDashboard({
  projectId,
  onBack,
  onDeleted,
  onRunSCF,
  onRunBands,
  onViewBands,
  onRunPhonons,
  onViewPhonons,
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

  // Delete project confirmation dialog state
  const [showDeleteDialog, setShowDeleteDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeleting, setIsDeleting] = useState(false);

  // Delete calculation confirmation dialog state
  const [showDeleteCalcDialog, setShowDeleteCalcDialog] = useState(false);
  const [calcToDelete, setCalcToDelete] = useState<{ calcId: string; calcType: string } | null>(null);
  const [isDeletingCalc, setIsDeletingCalc] = useState(false);

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

  function openDeleteCalcDialog(calcId: string, calcType: string) {
    setCalcToDelete({ calcId, calcType });
    setShowDeleteCalcDialog(true);
  }

  async function handleConfirmDeleteCalc() {
    if (!calcToDelete || !selectedCifId) return;

    setIsDeletingCalc(true);
    try {
      await invoke("delete_calculation", {
        projectId,
        cifId: selectedCifId,
        calcId: calcToDelete.calcId,
      });
      // Reload project to reflect changes
      await loadProject();
      setShowDeleteCalcDialog(false);
      setCalcToDelete(null);
    } catch (e) {
      console.error("Failed to delete calculation:", e);
      setError(String(e));
    } finally {
      setIsDeletingCalc(false);
    }
  }

  function handleRunSCF() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunSCF(selectedCifId, crystalData, cifContent, variant.filename);
  }

  function handleRunOptimization() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunSCF(selectedCifId, crystalData, cifContent, variant.filename, "relax", true);
  }

  function handleRunBands() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    // Pass all calculations for this CIF - the wizard will filter for SCF
    onRunBands(selectedCifId, crystalData, variant.calculations);
  }

  function handleRunPhonons() {
    if (!selectedCifId || !crystalData) return;
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return;
    onRunPhonons(selectedCifId, crystalData, variant.calculations);
  }

  function hasConvergedSCF(): boolean {
    const variant = project?.cif_variants.find(v => v.id === selectedCifId);
    if (!variant) return false;
    return variant.calculations.some(c => c.calc_type === "scf" && c.result?.converged);
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
          <label className="structure-selector-label">Structure</label>
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
            title="Add new structure"
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
            <button
              className="calc-action-btn"
              onClick={handleRunBands}
              disabled={!hasConvergedSCF()}
            >
              <span className="calc-action-icon">Band</span>
              <span className="calc-action-label">Band Structure</span>
              <span className="calc-action-hint">
                {hasConvergedSCF() ? "Electronic bands" : "Requires SCF"}
              </span>
            </button>
            <button
              className="calc-action-btn"
              onClick={handleRunPhonons}
              disabled={!hasConvergedSCF()}
            >
              <span className="calc-action-icon">Ph</span>
              <span className="calc-action-label">Phonons</span>
              <span className="calc-action-hint">
                {hasConvergedSCF() ? "DOS & Dispersion" : "Requires SCF"}
              </span>
            </button>
            <button className="calc-action-btn" onClick={handleRunOptimization}>
              <span className="calc-action-icon">Opt</span>
              <span className="calc-action-label">Geometry Optimization</span>
              <span className="calc-action-hint">VC-Relax preset</span>
            </button>
          </div>
        </section>

        {/* SCF Calculations */}
        {selectedVariant && selectedVariant.calculations.filter(c => c.calc_type === "scf").length > 0 && (
          <section className="history-section">
            <h3>SCFs</h3>
            <div className="calculations-list">
              {selectedVariant.calculations.filter(c => c.calc_type === "scf").map((calc) => (
                <div key={calc.id} className="calculation-item">
                  <div
                    className="calculation-header"
                    onClick={() =>
                      setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                    }
                  >
                    <div className="calculation-info">
                      <span className="calc-type">SCF</span>
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
                      <div className="calc-tags">
                        {getCalculationTags(calc).map((tag, i) => (
                          <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
                            {tag.label}
                          </span>
                        ))}
                      </div>
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
                      <div className="calc-actions">
                        <button
                          className="delete-calc-btn"
                          onClick={(e) => {
                            e.stopPropagation();
                            openDeleteCalcDialog(calc.id, calc.calc_type);
                          }}
                        >
                          Delete Calculation
                        </button>
                      </div>
                    </div>
                  )}
                </div>
              ))}
            </div>
          </section>
        )}

        {/* Band Structure Calculations */}
        {selectedVariant && selectedVariant.calculations.filter(c => c.calc_type === "bands").length > 0 && (
          <section className="history-section bands-section">
            <h3>Bands</h3>
            <div className="calculations-list">
              {selectedVariant.calculations.filter(c => c.calc_type === "bands").map((calc) => (
                <div key={calc.id} className="calculation-item bands-item">
                  <div
                    className="calculation-header"
                    onClick={() =>
                      setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                    }
                  >
                    <div className="calculation-info">
                      <span className="calc-type">BANDS</span>
                      {calc.parameters?.k_path && (
                        <span className="calc-kpath">{calc.parameters.k_path}</span>
                      )}
                      <div className="calc-tags">
                        {getBandsTags(calc).map((tag, i) => (
                          <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
                            {tag.label}
                          </span>
                        ))}
                      </div>
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

                  {expandedCalc === calc.id && (
                    <div className="calculation-details">
                      <div className="details-grid">
                        <div className="detail-item">
                          <label>K-Path</label>
                          <span>{calc.parameters?.k_path || "N/A"}</span>
                        </div>
                        <div className="detail-item">
                          <label>Total K-Points</label>
                          <span>{calc.parameters?.total_k_points || "N/A"}</span>
                        </div>
                        <div className="detail-item">
                          <label>Number of Bands</label>
                          <span>{calc.parameters?.n_bands || "N/A"}</span>
                        </div>
                        {calc.result?.fermi_energy && (
                          <div className="detail-item">
                            <label>Fermi Energy</label>
                            <span>{calc.result.fermi_energy.toFixed(4)} eV</span>
                          </div>
                        )}
                        <div className="detail-item">
                          <label>Source SCF</label>
                          <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                        </div>
                      </div>
                      <div className="calc-actions">
                        {calc.result?.band_data && (
                          <button
                            className="view-bands-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              onViewBands(calc.result?.band_data, calc.result?.fermi_energy ?? null);
                            }}
                          >
                            View Bands
                          </button>
                        )}
                        <button
                          className="delete-calc-btn"
                          onClick={(e) => {
                            e.stopPropagation();
                            openDeleteCalcDialog(calc.id, calc.calc_type);
                          }}
                        >
                          Delete
                        </button>
                      </div>
                    </div>
                  )}
                </div>
              ))}
            </div>
          </section>
        )}

        {/* Phonon Calculations */}
        {selectedVariant && selectedVariant.calculations.filter(c => c.calc_type === "phonon").length > 0 && (
          <section className="history-section phonon-section">
            <h3>Phonons</h3>
            <div className="calculations-list">
              {selectedVariant.calculations.filter(c => c.calc_type === "phonon").map((calc) => (
                <div key={calc.id} className="calculation-item phonon-item">
                  <div
                    className="calculation-header"
                    onClick={() =>
                      setExpandedCalc(expandedCalc === calc.id ? null : calc.id)
                    }
                  >
                    <div className="calculation-info">
                      <span className="calc-type">PHONON</span>
                      {calc.parameters?.q_path && (
                        <span className="calc-kpath">{calc.parameters.q_path}</span>
                      )}
                      <div className="calc-tags">
                        {getPhononTags(calc).map((tag, i) => (
                          <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
                            {tag.label}
                          </span>
                        ))}
                      </div>
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

                  {expandedCalc === calc.id && (
                    <div className="calculation-details">
                      <div className="details-grid">
                        <div className="detail-item">
                          <label>Q-Grid</label>
                          <span>
                            {calc.parameters?.q_grid
                              ? `${calc.parameters.q_grid[0]}×${calc.parameters.q_grid[1]}×${calc.parameters.q_grid[2]}`
                              : "N/A"}
                          </span>
                        </div>
                        <div className="detail-item">
                          <label>Modes</label>
                          <span>{calc.parameters?.n_modes || "N/A"}</span>
                        </div>
                        <div className="detail-item">
                          <label>Q-Points</label>
                          <span>{calc.parameters?.n_qpoints || "N/A"}</span>
                        </div>
                        <div className="detail-item">
                          <label>DOS</label>
                          <span>{calc.parameters?.calculate_dos ? "Yes" : "No"}</span>
                        </div>
                        <div className="detail-item">
                          <label>Dispersion</label>
                          <span>{calc.parameters?.calculate_dispersion ? "Yes" : "No"}</span>
                        </div>
                        <div className="detail-item">
                          <label>Source SCF</label>
                          <span>{calc.parameters?.source_scf_id?.slice(0, 8) || "N/A"}</span>
                        </div>
                      </div>
                      <div className="calc-actions">
                        {(calc.result as any)?.phonon_data && (
                          <button
                            className="view-phonon-btn"
                            onClick={(e) => {
                              e.stopPropagation();
                              onViewPhonons((calc.result as any)?.phonon_data);
                            }}
                          >
                            View Phonons
                          </button>
                        )}
                        <button
                          className="delete-calc-btn"
                          onClick={(e) => {
                            e.stopPropagation();
                            openDeleteCalcDialog(calc.id, calc.calc_type);
                          }}
                        >
                          Delete
                        </button>
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

      {/* Delete Project Dialog */}
      {showDeleteDialog && renderDeleteDialog()}

      {/* Delete Calculation Dialog */}
      {showDeleteCalcDialog && (
        <div className="dialog-overlay" onClick={() => !isDeletingCalc && setShowDeleteCalcDialog(false)}>
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>Delete Calculation</h2>
              <button
                className="dialog-close"
                onClick={() => setShowDeleteCalcDialog(false)}
                disabled={isDeletingCalc}
              >
                &times;
              </button>
            </div>

            <div className="dialog-body">
              <p className="exit-warning">
                Are you sure you want to delete this {calcToDelete?.calcType.toUpperCase()} calculation?
              </p>
              <p className="exit-hint">
                This will permanently remove the calculation results and all associated input/output files.
              </p>
            </div>

            <div className="dialog-footer">
              <button
                className="dialog-btn cancel"
                onClick={() => setShowDeleteCalcDialog(false)}
                disabled={isDeletingCalc}
              >
                Cancel
              </button>
              <button
                className="dialog-btn delete"
                onClick={handleConfirmDeleteCalc}
                disabled={isDeletingCalc}
              >
                {isDeletingCalc ? "Deleting..." : "Delete Calculation"}
              </button>
            </div>
          </div>
        </div>
      )}
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
