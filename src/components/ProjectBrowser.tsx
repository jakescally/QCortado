// Project Browser - Grid view of all projects

import { useState, useEffect, useRef, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { open, save } from "@tauri-apps/plugin-dialog";
import { CreateProjectDialog } from "./CreateProjectDialog";
import { EditProjectDialog } from "./EditProjectDialog";

interface ProjectSummary {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  folder_id: string | null;
  formula: string | null;
  calculation_count: number;
  calculation_types?: string[];
  last_activity: string;
}

interface ProjectFolder {
  id: string;
  name: string;
  created_at: string;
}

type ProjectCalculationType =
  | "scf"
  | "bands"
  | "dos"
  | "phonon"
  | "optimization"
  | "fermi_surface";

const PROJECT_CALCULATION_TYPE_ORDER: ProjectCalculationType[] = [
  "scf",
  "bands",
  "dos",
  "phonon",
  "optimization",
  "fermi_surface",
];

const PROJECT_CALCULATION_TYPE_LABELS: Record<ProjectCalculationType, string> = {
  scf: "SCF",
  bands: "Bands",
  dos: "DOS",
  phonon: "Phonon",
  optimization: "Geometry Optimization",
  fermi_surface: "Fermi Surface",
};

type ProjectFilter = "all" | ProjectCalculationType;

interface ProjectWithCalculationTypes {
  project: ProjectSummary;
  calculationTypes: ProjectCalculationType[];
}

interface FolderWithProjects {
  folder: ProjectFolder;
  projects: ProjectWithCalculationTypes[];
  calculationTypeCounts: Record<ProjectCalculationType, number>;
}

const PROJECT_FILTER_ORDER: ProjectFilter[] = ["all", ...PROJECT_CALCULATION_TYPE_ORDER];

function createEmptyCalculationTypeCounts(): Record<ProjectCalculationType, number> {
  return {
    scf: 0,
    bands: 0,
    dos: 0,
    phonon: 0,
    optimization: 0,
    fermi_surface: 0,
  };
}

function normalizeProjectCalculationType(calcType: string): ProjectCalculationType | null {
  const normalized = calcType.trim().toLowerCase();
  if (normalized === "scf") return "scf";
  if (normalized === "band" || normalized === "bands") return "bands";
  if (normalized === "dos") return "dos";
  if (normalized === "phonon") return "phonon";
  if (
    normalized === "optimization"
    || normalized === "geometry_optimization"
    || normalized === "geometry optimization"
    || normalized === "relax"
    || normalized === "vcrelax"
  ) {
    return "optimization";
  }
  if (
    normalized === "fermi_surface"
    || normalized === "fermisurface"
    || normalized === "fermi-surface"
    || normalized === "fermi surface"
  ) {
    return "fermi_surface";
  }
  return null;
}

function getProjectCalculationTypes(project: ProjectSummary): ProjectCalculationType[] {
  const presentTypes = new Set<ProjectCalculationType>();
  for (const rawType of project.calculation_types ?? []) {
    const normalizedType = normalizeProjectCalculationType(rawType);
    if (normalizedType) {
      presentTypes.add(normalizedType);
    }
  }
  return PROJECT_CALCULATION_TYPE_ORDER.filter((calcType) => presentTypes.has(calcType));
}

interface ProjectBrowserProps {
  onBack: () => void;
  onSelectProject?: (projectId: string, folderId: string | null) => void;
  onProjectsChanged?: () => void;
  initialActiveFolderId?: string | null;
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

interface CreatedProject {
  id: string;
  name: string;
}

export function ProjectBrowser({
  onBack,
  onSelectProject,
  onProjectsChanged,
  initialActiveFolderId = null,
}: ProjectBrowserProps) {
  const [projects, setProjects] = useState<ProjectSummary[]>([]);
  const [folders, setFolders] = useState<ProjectFolder[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [statusMessage, setStatusMessage] = useState<string | null>(null);
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [editingProjectId, setEditingProjectId] = useState<string | null>(null);
  const [showCreateFolderDialog, setShowCreateFolderDialog] = useState(false);
  const [newFolderName, setNewFolderName] = useState("");
  const [renamingFolder, setRenamingFolder] = useState<ProjectFolder | null>(null);
  const [renameFolderName, setRenameFolderName] = useState("");
  const [isSavingFolder, setIsSavingFolder] = useState(false);
  const [activeFolderId, setActiveFolderId] = useState<string | null>(initialActiveFolderId);
  const [openProjectMenuId, setOpenProjectMenuId] = useState<string | null>(null);
  const [movingProjectId, setMovingProjectId] = useState<string | null>(null);
  const [showExportDialog, setShowExportDialog] = useState(false);
  const [selectedExportProjectId, setSelectedExportProjectId] = useState<string>("");
  const [isExporting, setIsExporting] = useState(false);
  const [isCancelingExport, setIsCancelingExport] = useState(false);
  const [isImporting, setIsImporting] = useState(false);
  const [activeProjectFilters, setActiveProjectFilters] = useState<ProjectCalculationType[]>([]);
  const [exportProgress, setExportProgress] = useState<ProjectArchiveExportProgress | null>(null);
  const [importProgress, setImportProgress] = useState<ProjectArchiveImportProgress | null>(null);
  const activeExportIdRef = useRef<string | null>(null);
  const activeImportIdRef = useRef<string | null>(null);

  const folderById = useMemo<Map<string, ProjectFolder>>(
    () => new Map(folders.map((folder) => [folder.id, folder])),
    [folders],
  );

  const activeFolder = useMemo<ProjectFolder | null>(() => {
    if (!activeFolderId) return null;
    return folderById.get(activeFolderId) ?? null;
  }, [activeFolderId, folderById]);

  const projectsWithCalculationTypes = useMemo<ProjectWithCalculationTypes[]>(
    () => projects.map((project) => ({
      project,
      calculationTypes: getProjectCalculationTypes(project),
    })),
    [projects],
  );

  const projectsForActiveFolder = useMemo<ProjectWithCalculationTypes[]>(() => {
    if (!activeFolderId) {
      return projectsWithCalculationTypes.filter(({ project }) => !project.folder_id);
    }
    return projectsWithCalculationTypes.filter(
      ({ project }) => project.folder_id === activeFolderId,
    );
  }, [projectsWithCalculationTypes, activeFolderId]);

  const hasActiveFilters = activeProjectFilters.length > 0;

  const calculationTypeProjectCounts = useMemo<Record<ProjectCalculationType, number>>(() => {
    const counts = createEmptyCalculationTypeCounts();
    for (const { calculationTypes } of projectsForActiveFolder) {
      for (const calcType of calculationTypes) {
        counts[calcType] += 1;
      }
    }
    return counts;
  }, [projectsForActiveFolder]);

  const foldersWithProjects = useMemo<FolderWithProjects[]>(() => {
    const projectsByFolderId = new Map<string, ProjectWithCalculationTypes[]>();
    for (const projectWithTypes of projectsWithCalculationTypes) {
      const folderId = projectWithTypes.project.folder_id;
      if (!folderId) continue;
      const existing = projectsByFolderId.get(folderId);
      if (existing) {
        existing.push(projectWithTypes);
      } else {
        projectsByFolderId.set(folderId, [projectWithTypes]);
      }
    }

    return folders.map((folder) => {
      const folderProjects = [...(projectsByFolderId.get(folder.id) ?? [])].sort((a, b) => (
        a.project.name.localeCompare(b.project.name, undefined, { sensitivity: "base" })
      ));
      const calculationTypeCounts = createEmptyCalculationTypeCounts();
      for (const entry of folderProjects) {
        for (const calcType of entry.calculationTypes) {
          calculationTypeCounts[calcType] += 1;
        }
      }
      return {
        folder,
        projects: folderProjects,
        calculationTypeCounts,
      };
    });
  }, [projectsWithCalculationTypes, folders]);

  const filteredProjects = useMemo<ProjectWithCalculationTypes[]>(() => {
    if (activeProjectFilters.length === 0) {
      return projectsForActiveFolder;
    }
    return projectsForActiveFolder.filter(({ calculationTypes }) => (
      activeProjectFilters.some((filterType) => calculationTypes.includes(filterType))
    ));
  }, [projectsForActiveFolder, activeProjectFilters]);

  function toggleProjectFilter(filterType: ProjectCalculationType) {
    setActiveProjectFilters((currentFilters) => {
      if (currentFilters.includes(filterType)) {
        return currentFilters.filter((current) => current !== filterType);
      }
      return [...currentFilters, filterType];
    });
  }

  useEffect(() => {
    void loadProjectsAndFolders();
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

  useEffect(() => {
    if (!activeFolderId) return;
    if (!folderById.has(activeFolderId)) {
      setActiveFolderId(null);
    }
  }, [activeFolderId, folderById]);

  useEffect(() => {
    setActiveFolderId(initialActiveFolderId ?? null);
  }, [initialActiveFolderId]);

  useEffect(() => {
    function handleDocumentClick(event: MouseEvent) {
      const target = event.target as HTMLElement | null;
      if (!target?.closest(".project-card-menu-wrapper")) {
        setOpenProjectMenuId(null);
      }
    }

    document.addEventListener("mousedown", handleDocumentClick);
    return () => {
      document.removeEventListener("mousedown", handleDocumentClick);
    };
  }, []);

  async function loadProjectsAndFolders(showLoadingState = true) {
    setError(null);
    if (showLoadingState) {
      setIsLoading(true);
    }
    try {
      const [projectList, folderList] = await Promise.all([
        invoke<ProjectSummary[]>("list_projects"),
        invoke<ProjectFolder[]>("list_project_folders"),
      ]);
      setProjects(projectList);
      setFolders(folderList);
      setSelectedExportProjectId((currentSelection) => {
        if (projectList.length === 0) return "";
        if (currentSelection && projectList.some((project) => project.id === currentSelection)) {
          return currentSelection;
        }
        return projectList[0].id;
      });
    } catch (e) {
      console.error("Failed to load projects/folders:", e);
      setError(String(e));
    } finally {
      if (showLoadingState) {
        setIsLoading(false);
      }
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
      await loadProjectsAndFolders(false);
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

  function openEditProjectDialog(projectId: string) {
    setOpenProjectMenuId(null);
    setError(null);
    setStatusMessage(null);
    setEditingProjectId(projectId);
  }

  function closeEditProjectDialog() {
    setEditingProjectId(null);
  }

  function handleProjectMetadataSaved(updatedProject: {
    id: string;
    name: string;
    description: string | null;
  }) {
    setProjects((currentProjects) => currentProjects.map((project) => (
      project.id === updatedProject.id
        ? { ...project, name: updatedProject.name, description: updatedProject.description }
        : project
    )));
    setStatusMessage(`Updated "${updatedProject.name}".`);
  }

  function openCreateFolderModal() {
    setError(null);
    setStatusMessage(null);
    setNewFolderName("");
    setShowCreateFolderDialog(true);
  }

  function openRenameFolderModal(folder: ProjectFolder) {
    setError(null);
    setStatusMessage(null);
    setRenamingFolder(folder);
    setRenameFolderName(folder.name);
  }

  function closeCreateFolderModal() {
    if (isSavingFolder) return;
    setShowCreateFolderDialog(false);
  }

  function closeRenameFolderModal() {
    if (isSavingFolder) return;
    setRenamingFolder(null);
  }

  async function handleCreateFolder() {
    if (!newFolderName.trim()) {
      setError("Folder name is required.");
      return;
    }

    setIsSavingFolder(true);
    setError(null);
    setStatusMessage(null);
    try {
      const folder = await invoke<ProjectFolder>("create_project_folder", {
        name: newFolderName.trim(),
      });
      await loadProjectsAndFolders(false);
      setShowCreateFolderDialog(false);
      setStatusMessage(`Created folder "${folder.name}".`);
    } catch (e) {
      console.error("Failed to create folder:", e);
      setError(String(e));
    } finally {
      setIsSavingFolder(false);
    }
  }

  async function handleRenameFolder() {
    if (!renamingFolder) return;
    if (!renameFolderName.trim()) {
      setError("Folder name is required.");
      return;
    }

    setIsSavingFolder(true);
    setError(null);
    setStatusMessage(null);
    try {
      const updatedFolder = await invoke<ProjectFolder>("rename_project_folder", {
        folderId: renamingFolder.id,
        name: renameFolderName.trim(),
      });
      await loadProjectsAndFolders(false);
      setRenamingFolder(null);
      setStatusMessage(`Renamed folder to "${updatedFolder.name}".`);
    } catch (e) {
      console.error("Failed to rename folder:", e);
      setError(String(e));
    } finally {
      setIsSavingFolder(false);
    }
  }

  async function handleMoveProject(projectId: string, targetFolderId: string | null) {
    const destinationFolderId = targetFolderId?.trim() || null;
    const project = projects.find((item) => item.id === projectId);
    if (!project) return;
    if ((project.folder_id ?? null) === destinationFolderId) return;

    setOpenProjectMenuId(null);
    setMovingProjectId(projectId);
    setError(null);
    setStatusMessage(null);
    try {
      await invoke("move_project_to_folder", {
        projectId,
        folderId: destinationFolderId,
      });
      await loadProjectsAndFolders(false);
      onProjectsChanged?.();
      if (!destinationFolderId) {
        setStatusMessage(`Moved "${project.name}" to root.`);
      } else {
        const destinationFolderName = folderById.get(destinationFolderId)?.name ?? "folder";
        setStatusMessage(`Moved "${project.name}" to "${destinationFolderName}".`);
      }
    } catch (e) {
      console.error("Failed to move project:", e);
      setError(String(e));
    } finally {
      setMovingProjectId(null);
    }
  }

  async function handleProjectCreated(project: CreatedProject) {
    setShowCreateDialog(false);
    setError(null);
    setStatusMessage(null);

    const targetFolderId = activeFolderId;
    try {
      if (targetFolderId) {
        await invoke("move_project_to_folder", {
          projectId: project.id,
          folderId: targetFolderId,
        });
        const targetFolderName = folderById.get(targetFolderId)?.name ?? "folder";
        setStatusMessage(`Created "${project.name}" in "${targetFolderName}".`);
      } else {
        setStatusMessage(`Created "${project.name}".`);
      }
    } catch (e) {
      console.error("Project created but move failed:", e);
      setError(`Project was created, but moving it failed: ${e}`);
    } finally {
      await loadProjectsAndFolders(false);
      onProjectsChanged?.();
    }
  }

  const editingProject = editingProjectId
    ? projects.find((project) => project.id === editingProjectId) ?? null
    : null;

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
          <button
            className="secondary-project-btn"
            onClick={openCreateFolderModal}
            disabled={isImporting || isExporting}
          >
            + New Folder
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
        ) : projects.length === 0 && folders.length === 0 ? (
          <div className="empty-state">
            <div className="empty-icon">📁</div>
            <h3>No Projects Yet</h3>
            <p>Create a project to organize your calculations</p>
            <button className="new-project-btn" onClick={() => setShowCreateDialog(true)}>
              Create Your First Project
            </button>
          </div>
        ) : (
          <>
            {!activeFolder && foldersWithProjects.length > 0 && (
              <div className="folder-section">
                <div className="folder-section-header">
                  <h3>Folders</h3>
                  <span>
                    {foldersWithProjects.length} folder{foldersWithProjects.length !== 1 ? "s" : ""}
                  </span>
                </div>
                <div className="folder-grid">
                  {foldersWithProjects.map(({ folder, projects: folderProjects, calculationTypeCounts }) => (
                    <div
                      key={folder.id}
                      className={`folder-card ${activeFolderId === folder.id ? "active" : ""}`}
                      onClick={() => {
                        setOpenProjectMenuId(null);
                        setActiveFolderId(folder.id);
                      }}
                      onKeyDown={(e) => {
                        if (e.key === "Enter" || e.key === " ") {
                          e.preventDefault();
                          setOpenProjectMenuId(null);
                          setActiveFolderId(folder.id);
                        }
                      }}
                      role="button"
                      tabIndex={0}
                    >
                      <div className="folder-card-header">
                        <h4>{folder.name}</h4>
                        <button
                          className="folder-rename-btn"
                          type="button"
                          onClick={(e) => {
                            e.stopPropagation();
                            openRenameFolderModal(folder);
                          }}
                        >
                          Rename
                        </button>
                      </div>

                      {folderProjects.length === 0 ? (
                        <p className="folder-empty-text">No projects in this folder yet.</p>
                      ) : (
                        <div className="folder-project-list">
                          {folderProjects.map(({ project }) => (
                            <span key={project.id} className="folder-project-title">
                              {project.name}
                            </span>
                          ))}
                        </div>
                      )}

                      <div className="folder-tags">
                        {PROJECT_CALCULATION_TYPE_ORDER.filter(
                          (calcType) => calculationTypeCounts[calcType] > 0,
                        ).map((calcType) => (
                          <span
                            key={calcType}
                            className={`project-calc-tag project-calc-tag-${calcType.replace(/_/g, "-")}`}
                          >
                            {PROJECT_CALCULATION_TYPE_LABELS[calcType]} {calculationTypeCounts[calcType]}
                          </span>
                        ))}
                        {PROJECT_CALCULATION_TYPE_ORDER.every(
                          (calcType) => calculationTypeCounts[calcType] === 0,
                        ) && <span className="folder-tag-empty">No tags yet</span>}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {activeFolder && (
              <div className="folder-browse-bar">
                <button
                  className="secondary-project-btn folder-back-btn"
                  type="button"
                  onClick={() => setActiveFolderId(null)}
                >
                  ← All Projects
                </button>
                <div className="folder-browse-meta">
                  <h3>{activeFolder.name}</h3>
                  <span>
                    {projectsForActiveFolder.length} project
                    {projectsForActiveFolder.length !== 1 ? "s" : ""}
                  </span>
                </div>
              </div>
            )}

            <div className="project-filter-bar">
              {PROJECT_FILTER_ORDER.map((filterType) => {
                const filterLabel = filterType === "all"
                  ? "All"
                  : PROJECT_CALCULATION_TYPE_LABELS[filterType];
                const filterCount = filterType === "all"
                  ? projectsForActiveFolder.length
                  : calculationTypeProjectCounts[filterType];
                const variantClass = filterType === "all"
                  ? "project-filter-tab-all"
                  : `project-filter-tab-${filterType.replace(/_/g, "-")}`;
                const isActive = filterType === "all"
                  ? activeProjectFilters.length === 0
                  : activeProjectFilters.includes(filterType);
                return (
                  <button
                    key={filterType}
                    className={`project-filter-tab ${variantClass} ${isActive ? "active" : ""}`}
                    onClick={() => {
                      if (filterType === "all") {
                        setActiveProjectFilters([]);
                        return;
                      }
                      toggleProjectFilter(filterType);
                    }}
                    type="button"
                  >
                    <span>{filterLabel}</span>
                    <span className="project-filter-count">{filterCount}</span>
                  </button>
                );
              })}
            </div>

            {filteredProjects.length === 0 ? (
              <div className="empty-state project-filter-empty-state">
                {projectsForActiveFolder.length === 0 && activeFolder ? (
                  <>
                    <h3>No Projects In This Folder</h3>
                    <p>Move a project into {activeFolder.name} using a project menu.</p>
                    <button
                      className="secondary-project-btn"
                      type="button"
                      onClick={() => setActiveFolderId(null)}
                    >
                      View All Projects
                    </button>
                  </>
                ) : projectsForActiveFolder.length === 0 && !activeFolder && !hasActiveFilters ? (
                  <>
                    <h3>No Root Projects</h3>
                    <p>Projects moved into folders are hidden from root.</p>
                  </>
                ) : projects.length === 0 ? (
                  <>
                    <h3>No Projects Yet</h3>
                    <p>Create a project and place it in root or inside a folder.</p>
                    <button className="new-project-btn" onClick={() => setShowCreateDialog(true)}>
                      Create Project
                    </button>
                  </>
                ) : (
                  <>
                    <h3>No Matching Projects</h3>
                    <p>No projects include the selected tags.</p>
                    <button
                      className="secondary-project-btn"
                      type="button"
                      onClick={() => setActiveProjectFilters([])}
                    >
                      Clear Filters
                    </button>
                  </>
                )}
              </div>
            ) : (
              <div className="project-grid">
                {filteredProjects.map(({ project, calculationTypes }) => (
                  <div
                    key={project.id}
                    className="project-card"
                    onClick={() => onSelectProject?.(project.id, project.folder_id ?? null)}
                  >
                    <div className="project-card-header">
                      <h3 className="project-name">{project.name}</h3>
                      <div className="project-card-menu-wrapper" onClick={(e) => e.stopPropagation()}>
                        <button
                          className="project-card-menu-btn"
                          type="button"
                          onClick={() => {
                            setOpenProjectMenuId((current) => (
                              current === project.id ? null : project.id
                            ));
                          }}
                          title="Project options"
                          aria-label={`Open options for ${project.name}`}
                        >
                          ⋮
                        </button>
                        {openProjectMenuId === project.id && (
                          <div className="project-card-menu">
                            <button
                              type="button"
                              className="project-card-menu-item"
                              onClick={() => {
                                openEditProjectDialog(project.id);
                              }}
                            >
                              Edit Project
                            </button>
                            <div className="project-card-menu-divider" />
                            <button
                              type="button"
                              className="project-card-menu-item"
                              onClick={() => {
                                void handleMoveProject(project.id, null);
                              }}
                              disabled={movingProjectId === project.id || !project.folder_id}
                            >
                              Move to Root
                            </button>
                            {folders.map((folder) => (
                              <button
                                key={`${project.id}-${folder.id}`}
                                type="button"
                                className="project-card-menu-item"
                                onClick={() => {
                                  void handleMoveProject(project.id, folder.id);
                                }}
                                disabled={movingProjectId === project.id || project.folder_id === folder.id}
                              >
                                Move to {folder.name}
                              </button>
                            ))}
                          </div>
                        )}
                      </div>
                    </div>

                    {project.description && (
                      <p className="project-description">{project.description}</p>
                    )}

                    <div className="project-meta">
                      {project.formula && (
                        <span className="project-formula">{project.formula}</span>
                      )}
                      <span className="project-folder-pill">
                        {project.folder_id
                          ? (folderById.get(project.folder_id)?.name ?? "Unknown folder")
                          : "Root"}
                      </span>
                      <span className="project-stats">
                        {project.calculation_count} calculation
                        {project.calculation_count !== 1 ? "s" : ""}
                      </span>
                    </div>

                    {calculationTypes.length > 0 && (
                      <div className="project-calculation-tags">
                        {calculationTypes.map((calcType) => (
                          <span
                            key={calcType}
                            className={`project-calc-tag project-calc-tag-${calcType.replace(/_/g, "-")}`}
                          >
                            {PROJECT_CALCULATION_TYPE_LABELS[calcType]}
                          </span>
                        ))}
                      </div>
                    )}

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
          </>
        )}
      </div>

      <CreateProjectDialog
        isOpen={showCreateDialog}
        onClose={() => setShowCreateDialog(false)}
        onCreated={(project) => {
          void handleProjectCreated(project);
        }}
      />

      <EditProjectDialog
        isOpen={Boolean(editingProject)}
        projectId={editingProject?.id ?? null}
        initialName={editingProject?.name ?? ""}
        initialDescription={editingProject?.description ?? null}
        onClose={closeEditProjectDialog}
        onSaved={handleProjectMetadataSaved}
      />

      {showCreateFolderDialog && (
        <div className="dialog-overlay" onClick={closeCreateFolderModal}>
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>New Folder</h2>
              <button className="dialog-close" onClick={closeCreateFolderModal} disabled={isSavingFolder}>
                &times;
              </button>
            </div>
            <div className="dialog-body">
              <div className="save-form">
                <div className="form-group">
                  <label>Folder Name *</label>
                  <input
                    type="text"
                    value={newFolderName}
                    onChange={(e) => setNewFolderName(e.target.value)}
                    placeholder="e.g., Catalyst studies"
                    autoFocus
                  />
                </div>
              </div>
            </div>
            <div className="dialog-footer">
              <button className="dialog-btn cancel" onClick={closeCreateFolderModal} disabled={isSavingFolder}>
                Cancel
              </button>
              <button
                className="dialog-btn save"
                onClick={() => {
                  void handleCreateFolder();
                }}
                disabled={isSavingFolder || !newFolderName.trim()}
              >
                {isSavingFolder ? "Creating..." : "Create Folder"}
              </button>
            </div>
          </div>
        </div>
      )}

      {renamingFolder && (
        <div className="dialog-overlay" onClick={closeRenameFolderModal}>
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>Rename Folder</h2>
              <button className="dialog-close" onClick={closeRenameFolderModal} disabled={isSavingFolder}>
                &times;
              </button>
            </div>
            <div className="dialog-body">
              <div className="save-form">
                <div className="form-group">
                  <label>Folder Name *</label>
                  <input
                    type="text"
                    value={renameFolderName}
                    onChange={(e) => setRenameFolderName(e.target.value)}
                    autoFocus
                  />
                </div>
              </div>
            </div>
            <div className="dialog-footer">
              <button className="dialog-btn cancel" onClick={closeRenameFolderModal} disabled={isSavingFolder}>
                Cancel
              </button>
              <button
                className="dialog-btn save"
                onClick={() => {
                  void handleRenameFolder();
                }}
                disabled={isSavingFolder || !renameFolderName.trim()}
              >
                {isSavingFolder ? "Saving..." : "Save Name"}
              </button>
            </div>
          </div>
        </div>
      )}

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
