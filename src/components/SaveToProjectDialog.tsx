// Save to Project Dialog - Modal for saving calculation results to projects

import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData } from "../lib/types";

interface ProjectSummary {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  formula: string | null;
  calculation_count: number;
  last_activity: string;
}

interface CifVariant {
  id: string;
  filename: string;
  formula: string;
  added_at: string;
  calculations: unknown[];
}

interface Project {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  cif_variants: CifVariant[];
}

interface CalculationData {
  calc_type: string;
  parameters: unknown;
  result: unknown;
  started_at: string;
  completed_at: string;
  input_content: string;
  output_content: string;
  tags?: string[];
}

interface CifData {
  filename: string;
  formula: string;
  content: string;
  crystal_data: CrystalData;
}

interface SaveToProjectDialogProps {
  isOpen: boolean;
  onClose: () => void;
  onSaved: () => void;
  calculationData: CalculationData;
  cifData: CifData;
  workingDir?: string;
  /** Pre-selected project/CIF context (from dashboard) */
  projectContext?: { projectId: string; cifId: string };
}

type SaveMode = "new" | "existing";

export function SaveToProjectDialog({
  isOpen,
  onClose,
  onSaved,
  calculationData,
  cifData,
  workingDir,
  projectContext,
}: SaveToProjectDialogProps) {
  // If we have a project context, default to existing mode
  const [mode, setMode] = useState<SaveMode>(projectContext ? "existing" : "new");
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [selectedProjectId, setSelectedProjectId] = useState<string>(projectContext?.projectId || "");
  const [selectedProject, setSelectedProject] = useState<Project | null>(null);
  const [selectedCifId, setSelectedCifId] = useState<string>(projectContext?.cifId || "");
  const [addNewCif, setAddNewCif] = useState(!projectContext);

  const [newProjectName, setNewProjectName] = useState("");
  const [newProjectDescription, setNewProjectDescription] = useState("");

  const [isSaving, setIsSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Load projects on open
  useEffect(() => {
    if (isOpen) {
      loadProjects();
      // Reset state based on context
      if (projectContext) {
        setMode("existing");
        setSelectedProjectId(projectContext.projectId);
        setSelectedCifId(projectContext.cifId);
        setAddNewCif(false);
      }
    }
  }, [isOpen, projectContext]);

  // Load full project when selection changes
  useEffect(() => {
    if (selectedProjectId) {
      loadProject(selectedProjectId);
    } else {
      setSelectedProject(null);
      setSelectedCifId("");
    }
  }, [selectedProjectId]);

  async function loadProjects() {
    try {
      const projectList = await invoke<ProjectSummary[]>("list_projects");
      setProjects(projectList);
      if (projectList.length === 0) {
        setMode("new");
      }
    } catch (e) {
      console.error("Failed to load projects:", e);
      setError(String(e));
    }
  }

  async function loadProject(projectId: string) {
    try {
      const project = await invoke<Project>("get_project", { projectId });
      setSelectedProject(project);
      // Check if this CIF already exists in the project
      const existingCif = project.cif_variants.find(
        (v) => v.formula === cifData.formula
      );
      if (existingCif) {
        setSelectedCifId(existingCif.id);
        setAddNewCif(false);
      } else {
        setSelectedCifId("");
        setAddNewCif(true);
      }
    } catch (e) {
      console.error("Failed to load project:", e);
      setError(String(e));
    }
  }

  async function handleSave() {
    setIsSaving(true);
    setError(null);

    try {
      let projectId: string;
      let cifId: string;

      if (mode === "new") {
        // Create new project
        if (!newProjectName.trim()) {
          throw new Error("Project name is required");
        }

        const project = await invoke<Project>("create_project", {
          name: newProjectName.trim(),
          description: newProjectDescription.trim() || null,
        });
        projectId = project.id;

        // Add CIF to new project
        const variant = await invoke<CifVariant>("add_cif_to_project", {
          projectId,
          cifData: {
            filename: cifData.filename,
            formula: cifData.formula,
            content: cifData.content,
            crystal_data: cifData.crystal_data,
          },
        });
        cifId = variant.id;
      } else {
        // Add to existing project
        if (!selectedProjectId) {
          throw new Error("Please select a project");
        }
        projectId = selectedProjectId;

        if (addNewCif || !selectedCifId) {
          // Add new CIF to existing project
          const variant = await invoke<CifVariant>("add_cif_to_project", {
            projectId,
            cifData: {
              filename: cifData.filename,
              formula: cifData.formula,
              content: cifData.content,
              crystal_data: cifData.crystal_data,
            },
          });
          cifId = variant.id;
        } else {
          cifId = selectedCifId;
        }
      }

      // Save the calculation
      await invoke("save_calculation", {
        projectId,
        cifId,
        calcData: calculationData,
        workingDir: workingDir || null,
      });

      onSaved();
      onClose();
    } catch (e) {
      console.error("Failed to save:", e);
      setError(String(e));
    } finally {
      setIsSaving(false);
    }
  }

  if (!isOpen) return null;

  return (
    <div className="dialog-overlay" onClick={onClose}>
      <div className="dialog-content" onClick={(e) => e.stopPropagation()}>
        <div className="dialog-header">
          <h2>Save to Project</h2>
          <button className="dialog-close" onClick={onClose}>
            &times;
          </button>
        </div>

        <div className="dialog-body">
          {error && <div className="dialog-error">{error}</div>}

          <div className="save-mode-tabs">
            <button
              className={`mode-tab ${mode === "new" ? "active" : ""}`}
              onClick={() => setMode("new")}
            >
              Create New Project
            </button>
            <button
              className={`mode-tab ${mode === "existing" ? "active" : ""}`}
              onClick={() => setMode("existing")}
              disabled={projects.length === 0}
            >
              Add to Existing
            </button>
          </div>

          {mode === "new" && (
            <div className="save-form">
              <div className="form-group">
                <label>Project Name *</label>
                <input
                  type="text"
                  value={newProjectName}
                  onChange={(e) => setNewProjectName(e.target.value)}
                  placeholder="e.g., Silicon Band Structure Study"
                  autoFocus
                />
              </div>
              <div className="form-group">
                <label>Description (optional)</label>
                <textarea
                  value={newProjectDescription}
                  onChange={(e) => setNewProjectDescription(e.target.value)}
                  placeholder="Brief description of the project..."
                  rows={3}
                />
              </div>
            </div>
          )}

          {mode === "existing" && (
            <div className="save-form">
              <div className="form-group">
                <label>Select Project</label>
                <select
                  value={selectedProjectId}
                  onChange={(e) => setSelectedProjectId(e.target.value)}
                >
                  <option value="">Choose a project...</option>
                  {projects.map((p) => (
                    <option key={p.id} value={p.id}>
                      {p.name} {p.formula && `(${p.formula})`}
                    </option>
                  ))}
                </select>
              </div>

              {selectedProject && selectedProject.cif_variants.length > 0 && (
                <div className="form-group">
                  <label>CIF Structure</label>
                  <div className="cif-options">
                    <label className="radio-option">
                      <input
                        type="radio"
                        checked={addNewCif}
                        onChange={() => {
                          setAddNewCif(true);
                          setSelectedCifId("");
                        }}
                      />
                      <span>Add as new structure ({cifData.formula})</span>
                    </label>
                    {selectedProject.cif_variants.map((v) => (
                      <label key={v.id} className="radio-option">
                        <input
                          type="radio"
                          checked={!addNewCif && selectedCifId === v.id}
                          onChange={() => {
                            setAddNewCif(false);
                            setSelectedCifId(v.id);
                          }}
                        />
                        <span>
                          {v.formula} ({v.calculations.length} calculations)
                        </span>
                      </label>
                    ))}
                  </div>
                </div>
              )}
            </div>
          )}

          <div className="calculation-summary">
            <h4>Calculation to Save</h4>
            <div className="summary-grid">
              <span className="summary-label">Type:</span>
              <span className="summary-value">{calculationData.calc_type.toUpperCase()}</span>
              <span className="summary-label">Structure:</span>
              <span className="summary-value">{cifData.formula}</span>
              <span className="summary-label">Completed:</span>
              <span className="summary-value">
                {new Date(calculationData.completed_at).toLocaleString()}
              </span>
            </div>
          </div>
        </div>

        <div className="dialog-footer">
          <button className="dialog-btn cancel" onClick={onClose}>
            Cancel
          </button>
          <button
            className="dialog-btn save"
            onClick={handleSave}
            disabled={isSaving || (mode === "new" && !newProjectName.trim())}
          >
            {isSaving ? "Saving..." : "Save to Project"}
          </button>
        </div>
      </div>
    </div>
  );
}
