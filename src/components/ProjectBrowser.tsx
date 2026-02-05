// Project Browser - Grid view of all projects

import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
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
}

export function ProjectBrowser({
  onBack,
  onSelectProject,
}: ProjectBrowserProps) {
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [showCreateDialog, setShowCreateDialog] = useState(false);

  useEffect(() => {
    loadProjects();
  }, []);

  async function loadProjects() {
    setIsLoading(true);
    setError(null);
    try {
      const projectList = await invoke<ProjectSummary[]>("list_projects");
      setProjects(projectList);
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

  return (
    <div className="browser-container">
      <div className="browser-header">
        <button className="back-btn" onClick={onBack}>
          ← Back
        </button>
        <h2>Projects</h2>
        <button className="new-project-btn" onClick={() => setShowCreateDialog(true)}>
          + New Project
        </button>
      </div>

      {error && <div className="error-banner">{error}</div>}

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
          loadProjects();
        }}
      />
    </div>
  );
}
