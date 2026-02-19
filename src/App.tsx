import { useEffect, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { open } from "@tauri-apps/plugin-dialog";
import "./App.css";
import { SCFWizard } from "./components/SCFWizard";
import { BandStructureWizard } from "./components/BandStructureWizard";
import { BandData, BandPlot } from "./components/BandPlot";
import { ElectronicDOSWizard } from "./components/ElectronicDOSWizard";
import { ElectronicDOSData, ElectronicDOSPlot } from "./components/ElectronicDOSPlot";
import { FermiSurfaceWizard } from "./components/FermiSurfaceWizard";
import { PhononWizard } from "./components/PhononWizard";
import { PhononDOSPlot } from "./components/PhononPlot";
import { ProjectBrowser } from "./components/ProjectBrowser";
import { ProjectDashboard, CalculationRun } from "./components/ProjectDashboard";
import { CreateProjectDialog } from "./components/CreateProjectDialog";
import { ProcessIndicator } from "./components/ProcessIndicator";
import { TaskQueuePage } from "./components/TaskQueuePage";
import { TaskProvider } from "./lib/TaskContext";
import { ThemeProvider, useTheme } from "./lib/ThemeContext";
import { useWindowSize } from "./lib/useWindowSize";
import { clampMpiProcs, loadGlobalMpiDefaults, saveGlobalMpiDefaults } from "./lib/mpiDefaults";
import { SaveSizeMode, loadGlobalSaveSizeMode, saveGlobalSaveSizeMode } from "./lib/saveSizeMode";
import { CrystalData, SCFPreset, OptimizedStructureOption } from "./lib/types";

interface ProjectSummary {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  formula: string | null;
  calculation_count: number;
  last_activity: string;
}

interface TempCleanupResult {
  removed_paths: string[];
  failed_paths: string[];
  bytes_freed: number;
}

interface RecoveryCalculation {
  id: string;
  calc_type: string;
  started_at: string;
  completed_at: string | null;
  result: {
    converged: boolean;
  } | null;
}

interface RecoveryCifVariant {
  id: string;
  calculations: RecoveryCalculation[];
}

interface SettingsProjectSnapshot {
  id: string;
  name: string;
  last_opened_cif_id: string | null;
  cif_variants: RecoveryCifVariant[];
}

const DELETE_CONFIRM_TEXT = "delete my project for good";
const DEFAULT_FERMI_SURFER_PATH = "/usr/local/bin/fermisurfer";

type AppView = "home" | "scf-wizard" | "bands-wizard" | "bands-viewer" | "dos-wizard" | "dos-viewer" | "fermi-surface-wizard" | "phonon-wizard" | "phonon-viewer" | "project-browser" | "project-dashboard" | "task-queue";

interface SCFContext {
  cifId: string;
  crystalData: CrystalData;
  cifContent: string;
  filename: string;
  projectId: string;
  initialPreset?: SCFPreset;
  presetLock?: boolean;
  optimizedStructures?: OptimizedStructureOption[];
}

interface BandsContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface DosContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface FermiSurfaceContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface PhononsContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface PhononData {
  dos_data: any | null;
  dispersion_data: any | null;
}

type PhononViewMode = "bands" | "dos";
type PhononFrequencyUnit = "cm-1" | "thz";
type PhononBandFocus = "full" | "acoustic" | "optical";
const CM1_TO_THZ = 0.0299792458;

function toBandDataFromPhononDispersion(phononDispersion: any): BandData {
  return {
    k_points: phononDispersion.q_points || [],
    energies: phononDispersion.frequencies || [],
    fermi_energy: 0,
    high_symmetry_points: (phononDispersion.high_symmetry_points || []).map((point: any) => ({
      k_distance: point.q_distance,
      label: point.label,
    })),
    n_bands: phononDispersion.n_modes || 0,
    n_kpoints: phononDispersion.n_qpoints || 0,
    band_gap: null,
    energy_range: phononDispersion.frequency_range || [0, 0],
  };
}

function convertPhononBandDataUnit(data: BandData, unit: PhononFrequencyUnit): BandData {
  if (unit === "cm-1") {
    return data;
  }

  const scale = CM1_TO_THZ;
  return {
    ...data,
    energies: data.energies.map((band) => band.map((freq) => freq * scale)),
    energy_range: [
      data.energy_range[0] * scale,
      data.energy_range[1] * scale,
    ],
  };
}

function getPhononViewerRange(phononBandData: BandData): [number, number] {
  const rawMin = Number(phononBandData.energy_range?.[0]);
  const rawMax = Number(phononBandData.energy_range?.[1]);
  const safeMax = Number.isFinite(rawMax) ? rawMax : 0;
  const upperBase = Math.max(safeMax, 1);
  const padding = Math.max(0.2, upperBase * 0.08);
  const lower = 0;
  const upper = upperBase + padding;
  if (Number.isFinite(rawMin) && rawMin > upper) {
    return [0, rawMin + padding];
  }
  return [lower, upper];
}

function getFrequencyBounds(modes: number[][]): [number, number] | null {
  let min = Infinity;
  let max = -Infinity;

  for (const mode of modes) {
    for (const value of mode) {
      if (!Number.isFinite(value)) continue;
      if (value < min) min = value;
      if (value > max) max = value;
    }
  }

  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return null;
  }

  return [min, max];
}

function getPhononFocusRanges(phononBandData: BandData): {
  full: [number, number];
  acoustic: [number, number] | null;
  optical: [number, number] | null;
} {
  const full = getPhononViewerRange(phononBandData);
  const acousticCount = Math.min(3, phononBandData.energies.length);

  let acoustic: [number, number] | null = null;
  if (acousticCount > 0) {
    const acousticBounds = getFrequencyBounds(phononBandData.energies.slice(0, acousticCount));
    if (acousticBounds) {
      const upperBase = Math.max(acousticBounds[1], 1);
      const span = Math.max(upperBase - Math.max(acousticBounds[0], 0), 1);
      const padding = Math.max(0.2, span * 0.08);
      acoustic = [0, upperBase + padding];
    }
  }

  let optical: [number, number] | null = null;
  if (phononBandData.energies.length > acousticCount) {
    const opticalBounds = getFrequencyBounds(phononBandData.energies.slice(acousticCount));
    if (opticalBounds) {
      const span = Math.max(opticalBounds[1] - opticalBounds[0], 1);
      const padding = Math.max(0.2, span * 0.08);
      optical = [Math.max(0, opticalBounds[0] - padding * 0.6), opticalBounds[1] + padding];
    }
  }

  return { full, acoustic, optical };
}

function AppInner() {
  const { theme, setTheme } = useTheme();
  const windowSize = useWindowSize();
  const plotHeight = Math.max(400, windowSize.height - 160);

  const [qePath, setQePath] = useState<string | null>(null);
  const [qePathInput, setQePathInput] = useState("");
  const [fermiSurferPathInput, setFermiSurferPathInput] = useState(DEFAULT_FERMI_SURFER_PATH);
  const [isSavingQePath, setIsSavingQePath] = useState(false);
  const [isSavingFermiSurferPath, setIsSavingFermiSurferPath] = useState(false);
  const [availableExecutables, setAvailableExecutables] = useState<string[]>([]);
  const [qeStatus, setQeStatus] = useState<"Found" | "Not configured" | "Not found">("Not configured");
  const [fermiSurferStatus, setFermiSurferStatus] = useState<"Found" | "Not configured" | "Not found">("Not configured");
  const [error, setError] = useState<string | null>(null);
  const [currentView, setCurrentView] = useState<AppView>("home");
  const [projectCount, setProjectCount] = useState<number>(0);
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [selectedProjectId, setSelectedProjectId] = useState<string | null>(null);
  const [projectBrowserFolderId, setProjectBrowserFolderId] = useState<string | null>(null);
  const [showCloseConfirm, setShowCloseConfirm] = useState(false);
  const [showQueueMenu, setShowQueueMenu] = useState(false);
  const [lastNonQueueView, setLastNonQueueView] = useState<AppView>("home");
  const queueMenuRef = useRef<HTMLDivElement | null>(null);
  const [showSettingsMenu, setShowSettingsMenu] = useState(false);
  const settingsMenuRef = useRef<HTMLDivElement | null>(null);
  const [executionPrefixInput, setExecutionPrefixInput] = useState("");
  const [isSavingExecutionPrefix, setIsSavingExecutionPrefix] = useState(false);
  const [prefixStatus, setPrefixStatus] = useState<string | null>(null);
  const [globalMpiEnabled, setGlobalMpiEnabled] = useState(false);
  const [globalMpiProcs, setGlobalMpiProcs] = useState(1);
  const [globalMpiCpuCount, setGlobalMpiCpuCount] = useState(1);
  const [isSavingGlobalMpi, setIsSavingGlobalMpi] = useState(false);
  const [globalMpiStatus, setGlobalMpiStatus] = useState<string | null>(null);
  const [saveSizeMode, setSaveSizeMode] = useState<SaveSizeMode>("large");
  const [isSavingSaveSizeMode, setIsSavingSaveSizeMode] = useState(false);
  const [saveSizeStatus, setSaveSizeStatus] = useState<string | null>(null);
  const [isClearingTempStorage, setIsClearingTempStorage] = useState(false);
  const [tempStorageStatus, setTempStorageStatus] = useState<string | null>(null);
  const [isRecoveringPhonon, setIsRecoveringPhonon] = useState(false);
  const [recoveryStatus, setRecoveryStatus] = useState<string | null>(null);
  const [showDeleteProjectDialog, setShowDeleteProjectDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeletingProject, setIsDeletingProject] = useState(false);
  const [deleteProjectSnapshot, setDeleteProjectSnapshot] = useState<SettingsProjectSnapshot | null>(null);
  const [projectDashboardRefreshToken, setProjectDashboardRefreshToken] = useState(0);

  // Active task ID for reconnection when navigating to wizard from indicator
  const [reconnectTaskId, setReconnectTaskId] = useState<string | null>(null);

  // Context for running SCF from a project
  const [scfContext, setScfContext] = useState<SCFContext | null>(null);

  // Context for running Bands from a project
  const [bandsContext, setBandsContext] = useState<BandsContext | null>(null);

  // Context for viewing saved band data
  const [viewBandsData, setViewBandsData] = useState<{ bandData: any; fermiEnergy: number | null } | null>(null);

  // Context for running Electronic DOS from a project
  const [dosContext, setDosContext] = useState<DosContext | null>(null);

  // Context for viewing saved DOS data
  const [viewDosData, setViewDosData] = useState<{ dosData: ElectronicDOSData; fermiEnergy: number | null } | null>(null);

  // Context for running Fermi-surface generation from a project
  const [fermiSurfaceContext, setFermiSurfaceContext] = useState<FermiSurfaceContext | null>(null);

  // Context for running Phonons from a project
  const [phononsContext, setPhononsContext] = useState<PhononsContext | null>(null);

  // Context for viewing saved phonon data
  const [viewPhononData, setViewPhononData] = useState<{ data: PhononData; mode: PhononViewMode } | null>(null);
  const [phononBandsUnit, setPhononBandsUnit] = useState<PhononFrequencyUnit>("cm-1");
  const [phononBandFocus, setPhononBandFocus] = useState<PhononBandFocus>("full");

  // Check for existing QE configuration on startup
  useEffect(() => {
    checkQEPath();
    loadProjectCount();
    void loadFermiSurferPath();
    void loadExecutionPrefix();
    void loadGlobalMpiSettings();
    void loadGlobalSaveSizeSetting();
  }, []);

  // Listen for close confirmation events from backend
  useEffect(() => {
    const unlisten = listen("confirm-close", () => {
      setShowCloseConfirm(true);
    });
    return () => { unlisten.then((fn) => fn()); };
  }, []);

  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (queueMenuRef.current && !queueMenuRef.current.contains(event.target as Node)) {
        setShowQueueMenu(false);
      }
      if (settingsMenuRef.current && !settingsMenuRef.current.contains(event.target as Node)) {
        setShowSettingsMenu(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  async function loadProjectCount() {
    try {
      const projects = await invoke<ProjectSummary[]>("list_projects");
      setProjectCount(projects.length);
    } catch (e) {
      console.log("Failed to load project count:", e);
    }
  }

  async function checkQEPath() {
    try {
      const path = await invoke<string | null>("get_qe_path");
      if (path) {
        setQePath(path);
        setQePathInput(path);
        await loadExecutables();
      } else {
        setQePath(null);
        setQePathInput("");
        setAvailableExecutables([]);
        setQeStatus("Not configured");
      }
    } catch (e) {
      console.log("No QE path configured yet");
      setQePath(null);
      setQePathInput("");
      setAvailableExecutables([]);
      setQeStatus("Not configured");
    }
  }

  async function loadExecutables() {
    try {
      const exes = await invoke<string[]>("check_qe_executables");
      setAvailableExecutables(exes);
      setQeStatus(exes.includes("pw.x") ? "Found" : "Not found");
      setError(null);
    } catch (e) {
      setAvailableExecutables([]);
      setQeStatus("Not found");
      setError(String(e));
    }
  }

  async function selectQEPath() {
    try {
      const selected = await open({
        directory: true,
        multiple: false,
        title: "Select Quantum ESPRESSO bin directory",
      });

      if (selected && typeof selected === "string") {
        setQePathInput(selected);
        setError(null);
      }
    } catch (e) {
      setError(String(e));
    }
  }

  async function saveQEPath() {
    const normalized = qePathInput.trim();
    if (!normalized) {
      setError("QE path cannot be empty.");
      return;
    }

    setIsSavingQePath(true);
    try {
      await invoke("set_qe_path", { path: normalized });
      setQePath(normalized);
      setQePathInput(normalized);
      await loadExecutables();
      setError(null);
    } catch (e) {
      setError(String(e));
    } finally {
      setIsSavingQePath(false);
    }
  }

  async function loadFermiSurferPath() {
    try {
      const path = await invoke<string | null>("get_fermi_surfer_path");
      if (path) {
        setFermiSurferPathInput(path);
        setFermiSurferStatus("Found");
      } else {
        setFermiSurferPathInput(DEFAULT_FERMI_SURFER_PATH);
        setFermiSurferStatus("Not configured");
      }
    } catch (e) {
      console.error("Failed to load FermiSurfer path:", e);
      setFermiSurferPathInput(DEFAULT_FERMI_SURFER_PATH);
      setFermiSurferStatus("Not found");
    }
  }

  async function selectFermiSurferPath() {
    try {
      const selected = await open({
        directory: false,
        multiple: false,
        defaultPath: fermiSurferPathInput || qePath || "/usr/local/bin",
        title: "Select FermiSurfer executable",
      });

      if (selected && typeof selected === "string") {
        setFermiSurferPathInput(selected);
        setError(null);
      }
    } catch (e) {
      setError(String(e));
    }
  }

  async function saveFermiSurferPath() {
    const normalized = fermiSurferPathInput.trim();
    setIsSavingFermiSurferPath(true);
    try {
      await invoke("set_fermi_surfer_path", {
        path: normalized.length > 0 ? normalized : null,
      });
      if (normalized.length > 0) {
        setFermiSurferPathInput(normalized);
        setFermiSurferStatus("Found");
      } else {
        setFermiSurferPathInput(DEFAULT_FERMI_SURFER_PATH);
        setFermiSurferStatus("Not configured");
      }
      setError(null);
    } catch (e) {
      setFermiSurferStatus("Not found");
      setError(String(e));
    } finally {
      setIsSavingFermiSurferPath(false);
    }
  }

  async function loadExecutionPrefix() {
    try {
      const prefix = await invoke<string | null>("get_execution_prefix");
      setExecutionPrefixInput(prefix ?? "");
    } catch (e) {
      console.error("Failed to load execution prefix:", e);
    }
  }

  async function saveExecutionPrefix() {
    const normalized = executionPrefixInput.trim();
    setIsSavingExecutionPrefix(true);
    setPrefixStatus(null);
    try {
      await invoke("set_execution_prefix", {
        prefix: normalized.length > 0 ? normalized : null,
      });
      setExecutionPrefixInput(normalized);
      setPrefixStatus("Saved");
    } catch (e) {
      console.error("Failed to save execution prefix:", e);
      setPrefixStatus("Failed to save");
    } finally {
      setIsSavingExecutionPrefix(false);
    }
  }

  async function loadGlobalMpiSettings() {
    try {
      const cores = await invoke<number>("get_cpu_count");
      const safeCores = Math.max(1, Math.floor(cores));
      setGlobalMpiCpuCount(safeCores);
      const defaults = await loadGlobalMpiDefaults(safeCores);
      setGlobalMpiEnabled(defaults.enabled);
      setGlobalMpiProcs(defaults.nprocs);
    } catch (e) {
      console.error("Failed to load global MPI defaults:", e);
    }
  }

  async function saveGlobalMpiSettings() {
    setIsSavingGlobalMpi(true);
    setGlobalMpiStatus(null);
    try {
      const saved = await saveGlobalMpiDefaults(
        { enabled: globalMpiEnabled, nprocs: globalMpiProcs },
        globalMpiCpuCount,
      );
      setGlobalMpiEnabled(saved.enabled);
      setGlobalMpiProcs(saved.nprocs);
      setGlobalMpiStatus("Saved");
    } catch (e) {
      console.error("Failed to save global MPI defaults:", e);
      setGlobalMpiStatus("Failed to save");
    } finally {
      setIsSavingGlobalMpi(false);
    }
  }

  async function loadGlobalSaveSizeSetting() {
    const mode = await loadGlobalSaveSizeMode();
    setSaveSizeMode(mode);
  }

  async function saveGlobalSaveSizeSetting() {
    setIsSavingSaveSizeMode(true);
    setSaveSizeStatus(null);
    try {
      const saved = await saveGlobalSaveSizeMode(saveSizeMode);
      setSaveSizeMode(saved);
      setSaveSizeStatus("Saved");
    } catch (e) {
      console.error("Failed to save global save-size mode:", e);
      setSaveSizeStatus("Failed to save");
    } finally {
      setIsSavingSaveSizeMode(false);
    }
  }

  function formatBytes(bytes: number): string {
    if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
    const units = ["B", "KB", "MB", "GB", "TB"];
    let value = bytes;
    let unitIdx = 0;
    while (value >= 1024 && unitIdx < units.length - 1) {
      value /= 1024;
      unitIdx += 1;
    }
    const precision = value >= 10 || unitIdx === 0 ? 0 : 1;
    return `${value.toFixed(precision)} ${units[unitIdx]}`;
  }

  async function clearTempStorage() {
    setIsClearingTempStorage(true);
    setTempStorageStatus(null);
    try {
      const result = await invoke<TempCleanupResult>("clear_temp_storage");
      if (result.failed_paths.length > 0) {
        setTempStorageStatus(
          `Removed ${result.removed_paths.length} item(s), but ${result.failed_paths.length} item(s) could not be removed.`,
        );
      } else if (result.removed_paths.length > 0) {
        setTempStorageStatus(
          `Cleared ${formatBytes(result.bytes_freed)} from ${result.removed_paths.length} item(s).`,
        );
      } else {
        setTempStorageStatus("No QCortado temporary data found.");
      }
    } catch (e) {
      console.error("Failed to clear temp storage:", e);
      setTempStorageStatus(`Failed to clear temporary storage: ${e}`);
    } finally {
      setIsClearingTempStorage(false);
    }
  }

  async function recoverPhononFromSettings() {
    if (!selectedProjectId) {
      setRecoveryStatus("Open a project to recover phonon data.");
      return;
    }

    setIsRecoveringPhonon(true);
    setRecoveryStatus(null);
    try {
      const project = await invoke<SettingsProjectSnapshot>("get_project", { projectId: selectedProjectId });
      const selectedVariant = (project.last_opened_cif_id && project.cif_variants.some((variant) => variant.id === project.last_opened_cif_id))
        ? project.cif_variants.find((variant) => variant.id === project.last_opened_cif_id) ?? null
        : project.cif_variants[0] ?? null;
      if (!selectedVariant) {
        setRecoveryStatus("No structure found in this project.");
        return;
      }

      const fallbackScf = selectedVariant.calculations
        .filter((calc) => calc.calc_type === "scf" && calc.result?.converged)
        .sort((a, b) => {
          const aTime = a.completed_at ?? a.started_at;
          const bTime = b.completed_at ?? b.started_at;
          return bTime.localeCompare(aTime);
        })[0];

      const defaultTmpDir = "/tmp/qcortado_phonon";
      try {
        await invoke("recover_phonon_calculation", {
          projectId: selectedProjectId,
          cifId: selectedVariant.id,
          workingDir: defaultTmpDir,
          sourceScfId: fallbackScf?.id ?? null,
        });
        setRecoveryStatus(`Recovered phonon calculation from ${defaultTmpDir}.`);
        setProjectDashboardRefreshToken((prev) => prev + 1);
        return;
      } catch {
        // Fall back to directory picker below.
      }

      const selected = await open({
        multiple: false,
        directory: true,
        defaultPath: defaultTmpDir,
        title: "Select phonon scratch directory",
      });

      if (!selected || Array.isArray(selected)) {
        setRecoveryStatus("Phonon recovery canceled.");
        return;
      }

      await invoke("recover_phonon_calculation", {
        projectId: selectedProjectId,
        cifId: selectedVariant.id,
        workingDir: selected,
        sourceScfId: fallbackScf?.id ?? null,
      });
      setRecoveryStatus(`Recovered phonon calculation from ${selected}.`);
      setProjectDashboardRefreshToken((prev) => prev + 1);
    } catch (e) {
      console.error("Failed to recover phonon calculation:", e);
      setRecoveryStatus(`Phonon recovery failed: ${e}`);
    } finally {
      setIsRecoveringPhonon(false);
    }
  }

  async function openDeleteProjectDialog() {
    if (!selectedProjectId) return;
    try {
      const project = await invoke<SettingsProjectSnapshot>("get_project", { projectId: selectedProjectId });
      setDeleteProjectSnapshot(project);
      setDeleteConfirmText("");
      setShowDeleteProjectDialog(true);
      setShowSettingsMenu(false);
    } catch (e) {
      console.error("Failed to load project for deletion:", e);
      setPrefixStatus("Failed to open delete dialog");
    }
  }

  async function handleConfirmDeleteProject() {
    if (!selectedProjectId || deleteConfirmText !== DELETE_CONFIRM_TEXT) return;
    setIsDeletingProject(true);
    try {
      await invoke("delete_project", { projectId: selectedProjectId });
      setShowDeleteProjectDialog(false);
      setDeleteProjectSnapshot(null);
      setDeleteConfirmText("");
      setScfContext(null);
      setBandsContext(null);
      setDosContext(null);
      setFermiSurfaceContext(null);
      setPhononsContext(null);
      setViewBandsData(null);
      setViewDosData(null);
      setViewPhononData(null);
      setReconnectTaskId(null);
      setSelectedProjectId(null);
      setCurrentView("project-browser");
      await loadProjectCount();
    } catch (e) {
      console.error("Failed to delete project:", e);
      setPrefixStatus("Failed to delete project");
    } finally {
      setIsDeletingProject(false);
    }
  }

  function handleNavigateToTask(taskId: string, taskType: string) {
    setReconnectTaskId(taskId);
    const viewMap: Record<string, AppView> = {
      scf: "scf-wizard",
      bands: "bands-wizard",
      dos: "dos-wizard",
      fermi_surface: "fermi-surface-wizard",
      phonon: "phonon-wizard",
    };
    const view = viewMap[taskType];
    if (view) {
      setCurrentView(view);
    }
  }

  function navigateToQueue() {
    setShowQueueMenu(false);
    if (currentView !== "task-queue") {
      setLastNonQueueView(currentView);
      setCurrentView("task-queue");
    }
  }

  function returnFromQueue() {
    const fallback: AppView = selectedProjectId ? "project-dashboard" : "home";
    const destination = lastNonQueueView === "task-queue" ? fallback : lastNonQueueView;
    setCurrentView(destination);
  }

  async function handleConfirmClose() {
    setShowCloseConfirm(false);
    await invoke("shutdown_and_close");
  }

  // Render the close confirmation modal overlay
  const closeConfirmModal = showCloseConfirm ? (
    <div className="close-confirm-overlay" onClick={() => setShowCloseConfirm(false)}>
      <div className="close-confirm-dialog" onClick={(e) => e.stopPropagation()}>
        <h3>Calculations Still Running</h3>
        <p>
          One or more calculations are still running. Closing the app will terminate them. Are you sure?
        </p>
        <div className="close-confirm-actions">
          <button onClick={() => setShowCloseConfirm(false)}>Cancel</button>
          <button className="close-confirm-danger" onClick={handleConfirmClose}>
            Close Anyway
          </button>
        </div>
      </div>
    </div>
  ) : null;

  // The process indicator is always rendered
  const processIndicator = <ProcessIndicator onNavigateToTask={handleNavigateToTask} />;
  const queueLauncher = (
    <div className="floating-queue" ref={queueMenuRef}>
      <button
        className="floating-queue-btn"
        onClick={() => {
          setShowSettingsMenu(false);
          setShowQueueMenu((prev) => !prev);
        }}
        title="Task queue menu"
      >
        ☰
      </button>
      {showQueueMenu && (
        <div className="floating-queue-menu">
          <button onClick={navigateToQueue}>View Queue</button>
        </div>
      )}
    </div>
  );

  const settingsLauncher = (
    <div className="floating-settings" ref={settingsMenuRef}>
      <button
        className="floating-settings-btn"
        onClick={() => {
          setShowQueueMenu(false);
          setShowSettingsMenu((prev) => !prev);
        }}
        title="Settings"
      >
        <svg width="24" height="24" viewBox="0 0 20 20" fill="currentColor">
          <path fillRule="evenodd" d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.061 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.533 1.533 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.533 1.533 0 01.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd" />
        </svg>
      </button>

      {showSettingsMenu && (
        <div className="settings-window-overlay" onClick={() => setShowSettingsMenu(false)}>
          <div className="floating-settings-menu" onClick={(event) => event.stopPropagation()} role="dialog" aria-modal="true" aria-label="Settings">
            <div className="settings-window-header">
              <h3>Settings</h3>
              <button
                className="settings-window-close"
                onClick={() => setShowSettingsMenu(false)}
                aria-label="Close settings"
              >
                &times;
              </button>
            </div>
            <div className="settings-window-content">
              <div className="settings-menu-section">
                <label className="settings-menu-label" htmlFor="global-execution-prefix-input">
                  MPI Command Path
                </label>
                <input
                  id="global-execution-prefix-input"
                  className="settings-menu-input"
                  value={executionPrefixInput}
                  onChange={(e) => {
                    setExecutionPrefixInput(e.target.value);
                    setPrefixStatus(null);
                  }}
                  placeholder="e.g. /opt/homebrew/bin/mpirun"
                />
                <p className="settings-menu-hint">
                  Path to MPI launcher command (recommended: absolute path). Prepended before every QE executable launch.
                </p>
                <button
                  className="settings-menu-item"
                  onClick={() => void saveExecutionPrefix()}
                  disabled={isSavingExecutionPrefix}
                >
                  {isSavingExecutionPrefix ? "Saving..." : "Save MPI Command"}
                </button>
                {prefixStatus && <div className="settings-menu-status">{prefixStatus}</div>}
              </div>

              <div className="settings-menu-divider" />
              <div className="settings-menu-section">
                <label className="settings-menu-label">MPI Defaults</label>
                <label className="toggle-label">
                  <input
                    type="checkbox"
                    checked={globalMpiEnabled}
                    onChange={(e) => {
                      setGlobalMpiEnabled(e.target.checked);
                      setGlobalMpiStatus(null);
                    }}
                  />
                  <span>Enable MPI by default</span>
                </label>
                <label className="settings-menu-label" htmlFor="global-mpi-procs-input">
                  Default MPI Processes
                </label>
                <input
                  id="global-mpi-procs-input"
                  type="number"
                  min={1}
                  max={globalMpiCpuCount}
                  className="settings-menu-input"
                  value={globalMpiProcs}
                  onChange={(e) => {
                    setGlobalMpiProcs(clampMpiProcs(Number.parseInt(e.target.value, 10), globalMpiCpuCount));
                    setGlobalMpiStatus(null);
                  }}
                />
                <p className="settings-menu-hint">
                  Used as the initial MPI option in all calculation wizards ({globalMpiCpuCount} cores available).
                </p>
                <button
                  className="settings-menu-item"
                  onClick={() => void saveGlobalMpiSettings()}
                  disabled={isSavingGlobalMpi}
                >
                  {isSavingGlobalMpi ? "Saving..." : "Save MPI Defaults"}
                </button>
                {globalMpiStatus && <div className="settings-menu-status">{globalMpiStatus}</div>}
              </div>

              <div className="settings-menu-divider" />
              <div className="settings-menu-section">
                <label className="settings-menu-label" htmlFor="global-save-size-mode">
                  Calculation Save Size
                </label>
                <select
                  id="global-save-size-mode"
                  className="settings-menu-input"
                  value={saveSizeMode}
                  onChange={(event) => {
                    const value = event.target.value === "small" ? "small" : "large";
                    setSaveSizeMode(value);
                    setSaveSizeStatus(null);
                  }}
                >
                  <option value="large">Large (full restart data)</option>
                  <option value="small">Small (strip wavefunction archives)</option>
                </select>
                <p className="settings-menu-hint">
                  `Small` keeps useful outputs while removing large `wfc*` scratch files from saved calculation folders.
                </p>
                <button
                  className="settings-menu-item"
                  onClick={() => void saveGlobalSaveSizeSetting()}
                  disabled={isSavingSaveSizeMode}
                >
                  {isSavingSaveSizeMode ? "Saving..." : "Save Size Mode"}
                </button>
                {saveSizeStatus && <div className="settings-menu-status">{saveSizeStatus}</div>}
              </div>

              <div className="settings-menu-divider" />
              <div className="settings-menu-section">
                <label className="settings-menu-label">Temporary Storage</label>
                <p className="settings-menu-hint">
                  Remove `/tmp` and system temp QCortado working folders.
                </p>
                <button
                  className="settings-menu-item warning"
                  onClick={() => void clearTempStorage()}
                  disabled={isClearingTempStorage}
                >
                  {isClearingTempStorage ? "Clearing..." : "Clear Temp Storage"}
                </button>
                {tempStorageStatus && <div className="settings-menu-status">{tempStorageStatus}</div>}
              </div>

              <div className="settings-menu-divider" />
              <div className="settings-menu-section">
                <label className="settings-menu-label">Recovery</label>
                <p className="settings-menu-hint">
                  Import a completed phonon scratch run into the active project history.
                </p>
                <button
                  className="settings-menu-item"
                  onClick={() => void recoverPhononFromSettings()}
                  disabled={isRecoveringPhonon || !selectedProjectId}
                >
                  {isRecoveringPhonon ? "Recovering..." : "Recover Phonon"}
                </button>
                {recoveryStatus && <div className="settings-menu-status">{recoveryStatus}</div>}
              </div>

              <div className="settings-menu-divider" />
              <div className="settings-menu-section">
                <label className="settings-menu-label">Theme</label>
                <div className="theme-toggle-group" role="group" aria-label="Theme">
                  <button
                    type="button"
                    className={`theme-toggle-btn ${theme === "system" ? "active" : ""}`}
                    onClick={() => setTheme("system")}
                  >
                    System
                  </button>
                  <button
                    type="button"
                    className={`theme-toggle-btn ${theme === "light" ? "active" : ""}`}
                    onClick={() => setTheme("light")}
                  >
                    Light
                  </button>
                  <button
                    type="button"
                    className={`theme-toggle-btn ${theme === "dark" ? "active" : ""}`}
                    onClick={() => setTheme("dark")}
                  >
                    Dark
                  </button>
                </div>
              </div>

              {selectedProjectId && (
                <>
                  <div className="settings-menu-divider" />
                  <button className="settings-menu-item danger" onClick={() => void openDeleteProjectDialog()}>
                    Delete Project
                  </button>
                </>
              )}
            </div>
          </div>
        </div>
      )}
    </div>
  );

  const deleteProjectModal = showDeleteProjectDialog && deleteProjectSnapshot ? (
    <div className="dialog-overlay" onClick={() => !isDeletingProject && setShowDeleteProjectDialog(false)}>
      <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
        <div className="dialog-header">
          <h2>Delete Project</h2>
          <button
            className="dialog-close"
            onClick={() => setShowDeleteProjectDialog(false)}
            disabled={isDeletingProject}
          >
            &times;
          </button>
        </div>

        <div className="dialog-body">
          <div className="delete-warning">
            <p>
              You are about to permanently delete <strong>{deleteProjectSnapshot.name}</strong> and all of its data:
            </p>
            <ul>
              <li>{deleteProjectSnapshot.cif_variants.length} structure{deleteProjectSnapshot.cif_variants.length !== 1 ? "s" : ""}</li>
              <li>
                {deleteProjectSnapshot.cif_variants.reduce((sum, variant) => sum + variant.calculations.length, 0)} calculation
                {deleteProjectSnapshot.cif_variants.reduce((sum, variant) => sum + variant.calculations.length, 0) !== 1 ? "s" : ""}
              </li>
              <li>All input/output files</li>
            </ul>
            <p className="delete-warning-emphasis">
              This action cannot be undone.
            </p>
          </div>

          <div className="form-group">
            <label>
              Type <code>{DELETE_CONFIRM_TEXT}</code> to confirm:
            </label>
            <input
              type="text"
              value={deleteConfirmText}
              onChange={(e) => setDeleteConfirmText(e.target.value)}
              placeholder={DELETE_CONFIRM_TEXT}
              disabled={isDeletingProject}
              autoFocus
            />
          </div>
        </div>

        <div className="dialog-footer">
          <button
            className="dialog-btn cancel"
            onClick={() => setShowDeleteProjectDialog(false)}
            disabled={isDeletingProject}
          >
            Cancel
          </button>
          <button
            className="dialog-btn delete"
            onClick={() => void handleConfirmDeleteProject()}
            disabled={deleteConfirmText !== DELETE_CONFIRM_TEXT || isDeletingProject}
          >
            {isDeletingProject ? "Deleting..." : "Delete Project"}
          </button>
        </div>
      </div>
    </div>
  ) : null;

  const appChrome = (
    <>
      {queueLauncher}
      {settingsLauncher}
      {processIndicator}
      {closeConfirmModal}
      {deleteProjectModal}
    </>
  );

  if (currentView === "task-queue") {
    return (
      <>
        <TaskQueuePage onBack={returnFromQueue} />
        {appChrome}
      </>
    );
  }

  if (currentView === "bands-wizard" && qePath && (bandsContext || reconnectTaskId)) {
    return (
      <>
        <BandStructureWizard
          qePath={qePath}
          onViewBands={(bandData, fermiEnergy) => {
            setViewBandsData({ bandData, fermiEnergy });
            setCurrentView("bands-viewer");
            setReconnectTaskId(null);
          }}
          onBack={() => {
            setCurrentView("project-dashboard");
            setBandsContext(null);
            setReconnectTaskId(null);
          }}
          projectId={bandsContext?.projectId ?? ""}
          cifId={bandsContext?.cifId ?? ""}
          crystalData={bandsContext?.crystalData ?? { a: 0, b: 0, c: 0, alpha: 0, beta: 0, gamma: 0, spaceGroup: "", formula: "", atoms: [], species: [] } as any}
          scfCalculations={bandsContext?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskId ?? undefined}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "bands-viewer" && viewBandsData) {
    return (
      <>
        <div className="bands-viewer-container">
          <div className="bands-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewBandsData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>Band Structure</h2>
          </div>
          <div className="bands-viewer-content">
            <BandPlot
              data={viewBandsData.bandData}
              scfFermiEnergy={viewBandsData.fermiEnergy ?? undefined}
              viewerType="electronic"
            />
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "dos-wizard" && qePath && (dosContext || reconnectTaskId)) {
    return (
      <>
        <ElectronicDOSWizard
          qePath={qePath}
          onViewDos={(dosData, fermiEnergy) => {
            setViewDosData({ dosData, fermiEnergy });
            setCurrentView("dos-viewer");
            setReconnectTaskId(null);
          }}
          onBack={() => {
            setCurrentView("project-dashboard");
            setDosContext(null);
            setReconnectTaskId(null);
          }}
          projectId={dosContext?.projectId ?? ""}
          cifId={dosContext?.cifId ?? ""}
          crystalData={dosContext?.crystalData ?? { a: 0, b: 0, c: 0, alpha: 0, beta: 0, gamma: 0, spaceGroup: "", formula: "", atoms: [], species: [] } as any}
          scfCalculations={dosContext?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskId ?? undefined}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "dos-viewer" && viewDosData) {
    return (
      <>
        <div className="bands-viewer-container">
          <div className="bands-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewDosData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>Electronic DOS</h2>
          </div>
          <div className="bands-viewer-content">
            <ElectronicDOSPlot
              data={{
                ...viewDosData.dosData,
                fermi_energy: viewDosData.dosData.fermi_energy ?? viewDosData.fermiEnergy,
              }}
            />
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "fermi-surface-wizard" && qePath && (fermiSurfaceContext || reconnectTaskId)) {
    return (
      <>
        <FermiSurfaceWizard
          qePath={qePath}
          onBack={() => {
            setCurrentView("project-dashboard");
            setFermiSurfaceContext(null);
            setReconnectTaskId(null);
          }}
          projectId={fermiSurfaceContext?.projectId ?? ""}
          cifId={fermiSurfaceContext?.cifId ?? ""}
          crystalData={fermiSurfaceContext?.crystalData ?? { a: 0, b: 0, c: 0, alpha: 0, beta: 0, gamma: 0, spaceGroup: "", formula: "", atoms: [], species: [] } as any}
          scfCalculations={fermiSurfaceContext?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskId ?? undefined}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "phonon-wizard" && qePath && (phononsContext || reconnectTaskId)) {
    return (
      <>
        <PhononWizard
          qePath={qePath}
          onViewPhonons={(phononData, viewMode) => {
            setViewPhononData({
              data: phononData,
              mode: viewMode,
            });
            setCurrentView("phonon-viewer");
            setReconnectTaskId(null);
          }}
          onBack={() => {
            setCurrentView("project-dashboard");
            setPhononsContext(null);
            setReconnectTaskId(null);
          }}
          projectId={phononsContext?.projectId ?? ""}
          cifId={phononsContext?.cifId ?? ""}
          crystalData={phononsContext?.crystalData ?? { a: 0, b: 0, c: 0, alpha: 0, beta: 0, gamma: 0, spaceGroup: "", formula: "", atoms: [], species: [] } as any}
          scfCalculations={phononsContext?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskId ?? undefined}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "phonon-viewer" && viewPhononData) {
    const phononData = viewPhononData.data;
    const showingBands = viewPhononData.mode === "bands";
    const showingDos = viewPhononData.mode === "dos";
    const hasDos = phononData.dos_data !== null;
    const hasDispersion = phononData.dispersion_data !== null;
    const phononBandData = hasDispersion
      ? toBandDataFromPhononDispersion(phononData.dispersion_data)
      : null;
    const displayPhononBandData = phononBandData
      ? convertPhononBandDataUnit(phononBandData, phononBandsUnit)
      : null;
    const phononUnitLabel = phononBandsUnit === "thz" ? "THz" : "cm^-1";
    const focusRanges = displayPhononBandData ? getPhononFocusRanges(displayPhononBandData) : null;
    const canShowAcoustic = focusRanges?.acoustic !== null;
    const canShowOptical = focusRanges?.optical !== null;
    const resolvedFocus: PhononBandFocus =
      phononBandFocus === "acoustic" && !canShowAcoustic
        ? "full"
        : phononBandFocus === "optical" && !canShowOptical
          ? "full"
          : phononBandFocus;
    const activePhononRange =
      !focusRanges
        ? null
        : resolvedFocus === "acoustic"
          ? focusRanges.acoustic
          : resolvedFocus === "optical"
            ? focusRanges.optical
            : focusRanges.full;
    const plotKey = activePhononRange
      ? `${phononBandsUnit}-${resolvedFocus}-${activePhononRange[0].toFixed(6)}-${activePhononRange[1].toFixed(6)}`
      : `${phononBandsUnit}-${resolvedFocus}`;

    return (
      <>
        <div className="phonon-viewer-container">
          <div className="phonon-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewPhononData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>{showingBands ? "Phonon Bands" : "Phonon DOS"}</h2>
            {showingBands && displayPhononBandData && (
              <div className="phonon-view-controls">
                <div className="phonon-unit-toggle" role="group" aria-label="Phonon frequency units">
                  <button
                    type="button"
                    className={`phonon-unit-btn ${phononBandsUnit === "cm-1" ? "active" : ""}`}
                    onClick={() => setPhononBandsUnit("cm-1")}
                  >
                    cm^-1
                  </button>
                  <button
                    type="button"
                    className={`phonon-unit-btn ${phononBandsUnit === "thz" ? "active" : ""}`}
                    onClick={() => setPhononBandsUnit("thz")}
                  >
                    THz
                  </button>
                </div>
                <div className="phonon-focus-toggle" role="group" aria-label="Phonon frequency focus">
                  <button
                    type="button"
                    className={`phonon-focus-btn ${resolvedFocus === "full" ? "active" : ""}`}
                    onClick={() => setPhononBandFocus("full")}
                  >
                    Full
                  </button>
                  <button
                    type="button"
                    className={`phonon-focus-btn ${resolvedFocus === "acoustic" ? "active" : ""}`}
                    onClick={() => setPhononBandFocus("acoustic")}
                    disabled={!canShowAcoustic}
                  >
                    Acoustic
                  </button>
                  <button
                    type="button"
                    className={`phonon-focus-btn ${resolvedFocus === "optical" ? "active" : ""}`}
                    onClick={() => setPhononBandFocus("optical")}
                    disabled={!canShowOptical}
                  >
                    Optical
                  </button>
                </div>
              </div>
            )}
          </div>
          <div className="phonon-viewer-content">
            {showingBands && displayPhononBandData && activePhononRange ? (
              <BandPlot
                key={plotKey}
                data={displayPhononBandData}
                energyRange={activePhononRange}
                showFermiLevel={false}
                yAxisLabel={`Frequency (${phononUnitLabel})`}
                pointLabel="Mode"
                valueLabel="Frequency"
                valueUnit={phononUnitLabel}
                valueDecimals={phononBandsUnit === "thz" ? 3 : 1}
                primaryCountLabel="modes"
                secondaryCountLabel="q-points"
                scrollHint="Scroll: zoom frequency | Shift+Scroll: pan"
                yClampRange={null}
                viewerType="phonon"
              />
            ) : showingDos && hasDos ? (
              <PhononDOSPlot
                data={phononData.dos_data}
                width={Math.min(500, windowSize.width - 80)}
                height={plotHeight}
              />
            ) : (
              <p>{showingBands ? "No phonon dispersion data available" : "No phonon DOS data available"}</p>
            )}
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "scf-wizard" && qePath) {
    return (
      <>
        <SCFWizard
          qePath={qePath}
          onBack={() => {
            if (scfContext) {
              setCurrentView("project-dashboard");
              setScfContext(null);
            } else {
              setCurrentView("home");
            }
            setReconnectTaskId(null);
          }}
          initialCif={scfContext || undefined}
          initialPreset={scfContext?.initialPreset}
          presetLock={scfContext?.presetLock}
          optimizedStructures={scfContext?.optimizedStructures}
          reconnectTaskId={reconnectTaskId ?? undefined}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "project-browser") {
    return (
      <>
        <ProjectBrowser
          initialActiveFolderId={projectBrowserFolderId}
          onBack={() => {
            setCurrentView("home");
            loadProjectCount();
          }}
          onProjectsChanged={() => {
            void loadProjectCount();
          }}
          onSelectProject={(projectId, folderId) => {
            setSelectedProjectId(projectId);
            setProjectBrowserFolderId(folderId);
            setCurrentView("project-dashboard");
          }}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "project-dashboard" && selectedProjectId) {
    return (
      <>
        <ProjectDashboard
          projectId={selectedProjectId}
          showFloatingSettings={false}
          refreshToken={projectDashboardRefreshToken}
          onBack={() => {
            setCurrentView("project-browser");
            setSelectedProjectId(null);
          }}
          onDeleted={() => {
            setCurrentView("project-browser");
            setSelectedProjectId(null);
            loadProjectCount();
          }}
          onRunSCF={(cifId, crystalData, cifContent, filename, preset, presetLock, optimizedStructures) => {
            setScfContext({
              cifId,
              crystalData,
              cifContent,
              filename,
              projectId: selectedProjectId,
              initialPreset: preset,
              presetLock,
              optimizedStructures,
            });
            setCurrentView("scf-wizard");
          }}
          onRunBands={(cifId, crystalData, scfCalculations) => {
            setBandsContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            setCurrentView("bands-wizard");
          }}
          onViewBands={(bandData, fermiEnergy) => {
            setViewBandsData({ bandData, fermiEnergy });
            setCurrentView("bands-viewer");
          }}
          onRunDos={(cifId, crystalData, scfCalculations) => {
            setDosContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            setCurrentView("dos-wizard");
          }}
          onViewDos={(dosData, fermiEnergy) => {
            setViewDosData({ dosData, fermiEnergy });
            setCurrentView("dos-viewer");
          }}
          onRunFermiSurface={(cifId, crystalData, scfCalculations) => {
            setFermiSurfaceContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            setCurrentView("fermi-surface-wizard");
          }}
          onRunPhonons={(cifId, crystalData, scfCalculations) => {
            setPhononsContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            setCurrentView("phonon-wizard");
          }}
          onViewPhonons={(phononData, viewMode) => {
            setViewPhononData({
              data: phononData,
              mode: viewMode,
            });
            setCurrentView("phonon-viewer");
          }}
        />
        {appChrome}
      </>
    );
  }

  const availablePrograms: Array<{ name: string; type: "qe" | "fermisurfer" }> = [
    ...availableExecutables.map((name) => ({ name, type: "qe" as const })),
    ...(fermiSurferStatus === "Found"
      ? [{ name: "fermisurfer", type: "fermisurfer" as const }]
      : []),
  ];

  const qeStatusClass = qeStatus === "Found" ? "ready" : qeStatus === "Not found" ? "error" : "pending";
  const fermiStatusClass =
    fermiSurferStatus === "Found" ? "ready" : fermiSurferStatus === "Not found" ? "error" : "pending";

  return (
    <>
      <main className="container">
        <header className="header">
          <h1>QCortado</h1>
          <p className="subtitle">A Modern Interface for Quantum ESPRESSO</p>
        </header>

        <section className="config-section">
          <h2>Configuration</h2>
          <div className="config-row">
            <label>QE Installation:</label>
            <input
              type="text"
              className="config-path-input"
              value={qePathInput}
              onChange={(e) => {
                setQePathInput(e.target.value);
                setError(null);
              }}
              placeholder="/path/to/qe/bin"
              spellCheck={false}
            />
            <div className="config-row-actions">
              <button onClick={selectQEPath}>Browse</button>
              <button
                onClick={() => void saveQEPath()}
                disabled={isSavingQePath || qePathInput.trim().length === 0}
              >
                {isSavingQePath ? "Saving..." : "Save"}
              </button>
            </div>
          </div>

          <div className="config-row">
            <label>FermiSurfer:</label>
            <input
              type="text"
              className="config-path-input"
              value={fermiSurferPathInput}
              onChange={(e) => {
                setFermiSurferPathInput(e.target.value);
                setError(null);
              }}
              placeholder={DEFAULT_FERMI_SURFER_PATH}
              spellCheck={false}
            />
            <div className="config-row-actions">
              <button onClick={selectFermiSurferPath}>Browse</button>
              <button
                onClick={() => void saveFermiSurferPath()}
                disabled={isSavingFermiSurferPath}
              >
                {isSavingFermiSurferPath ? "Saving..." : "Save"}
              </button>
            </div>
          </div>

          <div className="status-row">
            <label>QE Status:</label>
            <span className={`status ${qeStatusClass}`}>
              {qeStatus}
            </span>
          </div>

          <div className="status-row">
            <label>FermiSurfer Status:</label>
            <span className={`status ${fermiStatusClass}`}>
              {fermiSurferStatus}
            </span>
          </div>

          {error && <div className="error">{error}</div>}

          {availablePrograms.length > 0 && (
            <div className="executables">
              <label>Available programs:</label>
              <div className="exe-list">
                {availablePrograms.map((program) => (
                  <span
                    key={program.name}
                    className={`exe-tag ${program.type === "fermisurfer" ? "exe-tag-fermisurfer" : ""}`}
                  >
                    {program.name}
                  </span>
                ))}
              </div>
            </div>
          )}
        </section>

        {qePath && (
          <>
            <section className="actions-section">
              <h2>Projects</h2>
              <div className="action-grid">
                <button className="action-btn" onClick={() => setShowCreateDialog(true)}>
                  <span className="action-icon">+</span>
                  <span className="action-label">New Project</span>
                  <span className="action-hint">Create a project</span>
                </button>
                <button
                  className="action-btn"
                  onClick={() => {
                    setProjectBrowserFolderId(null);
                    setCurrentView("project-browser");
                  }}
                  disabled={projectCount === 0}
                >
                  <span className="action-icon">{projectCount}</span>
                  <span className="action-label">
                    {projectCount === 1 ? "View Project" : "View Projects"}
                  </span>
                  <span className="action-hint">
                    {projectCount === 0 ? "No projects yet" : "Browse & manage"}
                  </span>
                </button>
              </div>
            </section>

          </>
        )}

        <footer className="footer">
          <p>QCortado v0.1.0 - Quantum ESPRESSO 7.5 Interface</p>
          <div className="theme-toggle-group" role="group" aria-label="Theme">
            <button
              type="button"
              className={`theme-toggle-btn ${theme === "system" ? "active" : ""}`}
              onClick={() => setTheme("system")}
            >
              System
            </button>
            <button
              type="button"
              className={`theme-toggle-btn ${theme === "light" ? "active" : ""}`}
              onClick={() => setTheme("light")}
            >
              Light
            </button>
            <button
              type="button"
              className={`theme-toggle-btn ${theme === "dark" ? "active" : ""}`}
              onClick={() => setTheme("dark")}
            >
              Dark
            </button>
          </div>
        </footer>

        <CreateProjectDialog
          isOpen={showCreateDialog}
          onClose={() => setShowCreateDialog(false)}
          onCreated={() => {
            setShowCreateDialog(false);
            loadProjectCount();
          }}
        />
      </main>
        {appChrome}
    </>
  );
}

function App() {
  return (
    <ThemeProvider>
      <TaskProvider>
        <AppInner />
      </TaskProvider>
    </ThemeProvider>
  );
}

export default App;
