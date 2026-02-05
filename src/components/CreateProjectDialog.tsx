// Create Project Dialog - Modal for creating a new project

import { useState } from "react";
import { invoke } from "@tauri-apps/api/core";

interface Project {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  cif_variants: unknown[];
}

interface CreateProjectDialogProps {
  isOpen: boolean;
  onClose: () => void;
  onCreated: (project: Project) => void;
}

export function CreateProjectDialog({
  isOpen,
  onClose,
  onCreated,
}: CreateProjectDialogProps) {
  const [name, setName] = useState("");
  const [description, setDescription] = useState("");
  const [isCreating, setIsCreating] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function handleCreate() {
    if (!name.trim()) {
      setError("Project name is required");
      return;
    }

    setIsCreating(true);
    setError(null);

    try {
      const project = await invoke<Project>("create_project", {
        name: name.trim(),
        description: description.trim() || null,
      });

      onCreated(project);
      handleClose();
    } catch (e) {
      console.error("Failed to create project:", e);
      setError(String(e));
    } finally {
      setIsCreating(false);
    }
  }

  function handleClose() {
    setName("");
    setDescription("");
    setError(null);
    onClose();
  }

  function handleKeyDown(e: React.KeyboardEvent) {
    if (e.key === "Enter" && !e.shiftKey && name.trim()) {
      e.preventDefault();
      handleCreate();
    }
    if (e.key === "Escape") {
      handleClose();
    }
  }

  if (!isOpen) return null;

  return (
    <div className="dialog-overlay" onClick={handleClose}>
      <div
        className="dialog-content dialog-small"
        onClick={(e) => e.stopPropagation()}
        onKeyDown={handleKeyDown}
      >
        <div className="dialog-header">
          <h2>New Project</h2>
          <button className="dialog-close" onClick={handleClose}>
            &times;
          </button>
        </div>

        <div className="dialog-body">
          {error && <div className="dialog-error">{error}</div>}

          <div className="save-form">
            <div className="form-group">
              <label>Project Name *</label>
              <input
                type="text"
                value={name}
                onChange={(e) => setName(e.target.value)}
                placeholder="e.g., Silicon Band Structure Study"
                autoFocus
              />
            </div>
            <div className="form-group">
              <label>Description (optional)</label>
              <textarea
                value={description}
                onChange={(e) => setDescription(e.target.value)}
                placeholder="Brief description of the project..."
                rows={3}
              />
            </div>
          </div>
        </div>

        <div className="dialog-footer">
          <button className="dialog-btn cancel" onClick={handleClose}>
            Cancel
          </button>
          <button
            className="dialog-btn save"
            onClick={handleCreate}
            disabled={isCreating || !name.trim()}
          >
            {isCreating ? "Creating..." : "Create Project"}
          </button>
        </div>
      </div>
    </div>
  );
}
