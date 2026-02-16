import { useEffect, useState } from "react";
import type { KeyboardEvent } from "react";
import { invoke } from "@tauri-apps/api/core";

interface UpdatedProjectMetadata {
  id: string;
  name: string;
  description: string | null;
}

interface EditProjectDialogProps {
  isOpen: boolean;
  projectId: string | null;
  initialName: string;
  initialDescription: string | null;
  onClose: () => void;
  onSaved: (updatedProject: UpdatedProjectMetadata) => void;
}

export function EditProjectDialog({
  isOpen,
  projectId,
  initialName,
  initialDescription,
  onClose,
  onSaved,
}: EditProjectDialogProps) {
  const [name, setName] = useState(initialName);
  const [description, setDescription] = useState(initialDescription ?? "");
  const [isSaving, setIsSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!isOpen) return;
    setName(initialName);
    setDescription(initialDescription ?? "");
    setError(null);
  }, [isOpen, initialName, initialDescription]);

  if (!isOpen || !projectId) return null;

  async function handleSave() {
    if (!name.trim()) {
      setError("Project name is required");
      return;
    }

    setIsSaving(true);
    setError(null);

    try {
      const updatedProject = await invoke<UpdatedProjectMetadata>("update_project_metadata", {
        projectId,
        name: name.trim(),
        description: description.trim() || null,
      });
      onSaved(updatedProject);
      onClose();
    } catch (e) {
      console.error("Failed to update project metadata:", e);
      setError(String(e));
    } finally {
      setIsSaving(false);
    }
  }

  function handleKeyDown(e: KeyboardEvent<HTMLDivElement>) {
    if (e.key === "Enter" && !e.shiftKey && name.trim()) {
      e.preventDefault();
      void handleSave();
    }
    if (e.key === "Escape" && !isSaving) {
      onClose();
    }
  }

  return (
    <div className="dialog-overlay" onClick={() => !isSaving && onClose()}>
      <div
        className="dialog-content dialog-small"
        onClick={(e) => e.stopPropagation()}
        onKeyDown={handleKeyDown}
      >
        <div className="dialog-header">
          <h2>Edit Project</h2>
          <button className="dialog-close" onClick={onClose} disabled={isSaving}>
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
          <button className="dialog-btn cancel" onClick={onClose} disabled={isSaving}>
            Cancel
          </button>
          <button
            className="dialog-btn save"
            onClick={() => {
              void handleSave();
            }}
            disabled={isSaving || !name.trim()}
          >
            {isSaving ? "Saving..." : "Save Changes"}
          </button>
        </div>
      </div>
    </div>
  );
}
