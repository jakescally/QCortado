import { useEffect, useMemo, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { getCurrentWindow } from "@tauri-apps/api/window";
import { open, save } from "@tauri-apps/plugin-dialog";
import "./App.css";
import { EpwViewer } from "./components/qe";
import type { EpwViewerPayload } from "./components/qe";
import {
  canRenderEngineWorkflowHost,
  EngineWorkflowHost,
} from "./components/engine";
import type { EngineWorkflowBackDestination } from "./components/engine";
import { BandPlot } from "./components/BandPlot";
import type { BandData, BandPlotData } from "./components/BandPlot";
import { ElectronicDOSData, ElectronicDOSPlot } from "./components/ElectronicDOSPlot";
import { PhononDOSPlot } from "./components/PhononPlot";
import { TransportPlot } from "./components/TransportPlot";
import { ProjectBrowser } from "./components/ProjectBrowser";
import {
  ProjectDashboard,
  CalculationRun,
  SavedBandsCalculationContext,
  WannierBandOverlayOption,
} from "./components/ProjectDashboard";
import { HpcActivityPanel } from "./components/HpcActivityPanel";
import { AppNavigationDrawer } from "./components/AppNavigationDrawer";
import { AppHeaderPortal } from "./components/AppHeaderPortal";
import { AppTopBar } from "./components/AppTopBar";
import { TaskManagerDrawer } from "./components/TaskManagerDrawer";
import { HpcSetupWizard } from "./components/HpcSetupWizard";
import { HpcProfileEditor } from "./components/HpcProfileEditor";
import { HpcNodeActivityPage } from "./components/HpcNodeActivityPage";
import { StorageManagerPage } from "./components/StorageManagerPage";
import { BandsMultiview } from "./components/BandsMultiview";
import { InfoTooltip } from "./components/InfoTooltip";
import type { BandsMultiviewCalculation } from "./components/BandsMultiview";
import { TaskProvider, useTaskContext } from "./lib/TaskContext";
import { ThemeProvider, useTheme } from "./lib/ThemeContext";
import { useWindowSize } from "./lib/useWindowSize";
import {
  buildTaskManagerEntries,
  findRelevantTaskManagerEntry,
  summarizeTaskManagerEntries,
} from "./lib/taskManager";
import type { TaskManagerFilter } from "./lib/taskManager";
import { clampMpiProcs, loadGlobalMpiDefaults, saveGlobalMpiDefaults } from "./lib/mpiDefaults";
import { SaveSizeMode, loadGlobalSaveSizeMode, saveGlobalSaveSizeMode } from "./lib/saveSizeMode";
import {
  listEngineInstallations,
  resolveEngineWorkflowHostRoute,
} from "./lib/engines";
import type { EngineInstallation, EngineWorkflowHostRoute } from "./lib/engines";
import type { EngineWorkflowView } from "./lib/engines/plugin";
import { formatWannierConvergenceFlag, getWannierQualityIssues } from "./lib/engines/qe/wannierQuality";
import {
  CrystalData,
  SCFPreset,
  OptimizedStructureOption,
  ExecutionMode,
  HpcLauncher,
  HpcProfile,
  HpcResourceMode,
  QeDefaults,
  QeSmearingType,
  SlurmResourceRequest,
} from "./lib/types";
import {
  cleanHpcRemoteOrphans,
  defaultCpuResources,
  defaultGpuResources,
  deleteHpcProfile,
  exportHpcPresetBundle,
  getActiveHpcProfileId,
  importHpcPresetBundle,
  attachHeadlessHpcJob,
  listHeadlessHpcJobs,
  listHpcProfiles,
  loadExecutionMode,
  migrateHpcRemoteRoots,
  normalizeCliDashText,
  saveExecutionMode,
  setActiveHpcProfile,
  updateHpcProfileDefaults,
} from "./lib/hpcConfig";
import { resolveProfileRemotePseudoDir } from "./lib/engines/qe/hpc";
import type { HpcHeadlessJobCandidate } from "./lib/hpcConfig";
import type { TransportResult } from "./lib/transport";
import type { CalculationKind, EngineId } from "./lib/engines/types";

interface TempCleanupResult {
  removed_paths: string[];
  failed_paths: string[];
  bytes_freed: number;
}

interface PslibraryPseudoRepairResult {
  pseudo_dir: string;
  scanned: number;
  candidates: number;
  patched: number;
  already_clean: number;
  patched_files: string[];
  clean_files: string[];
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

interface HpcRecoverableRemotePhononRun {
  profile_id: string;
  remote_workdir: string;
  location: string;
  modified_at_epoch: number;
}

interface HpcRemotePhononRecoveryDebugReport {
  profile_id: string;
  workspace_root: string;
  project_phonon_root: string;
  workspace_probe_output: string;
  project_probe_output: string;
  recoverable_runs: HpcRecoverableRemotePhononRun[];
}

interface RemotePhononRecoveryContext {
  projectId: string;
  cifId: string;
  sourceScfId: string | null;
}

const DELETE_CONFIRM_TEXT = "DELETE";
const DEFAULT_FERMI_SURFER_PATH = "/usr/local/bin/fermisurfer";
const DEFAULT_WANNIER90_PATH = "qe-7.5/external/wannier90/wannier90.x";
const DEFAULT_POSTW90_PATH = "qe-7.5/external/wannier90/postw90.x";
const DEFAULT_QE_DEFAULTS: QeDefaults = {
  smearing: "marzari-vanderbilt",
};

type AppView = "scf-wizard" | "bands-wizard" | "bands-viewer" | "bands-multiview" | "dos-wizard" | "dos-viewer" | "wannier-wizard" | "wannier-viewer" | "transport-wizard" | "transport-viewer" | "fermi-surface-wizard" | "hubbard-lrt-wizard" | "phonon-wizard" | "phonon-viewer" | "epw-wizard" | "epw-viewer" | "wien2k-structure-wizard" | "wien2k-scf-wizard" | "project-browser" | "project-dashboard" | "settings" | "node-activity" | "storage-manager";

interface OpenTaskViewRequest {
  taskId: string;
  taskType: string;
}

interface SCFContext {
  cifId: string;
  crystalData: CrystalData;
  cifContent: string;
  filename: string;
  projectId: string;
  initialPreset?: SCFPreset;
  presetLock?: boolean;
  optimizedStructures?: OptimizedStructureOption[];
  calculations?: CalculationRun[];
  continuationCalculationId?: string | null;
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

interface WannierContext {
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

interface TransportContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  wannierCalculations: CalculationRun[];
}

interface PhononsContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface EpwContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  calculations: CalculationRun[];
}

interface HubbardLrtContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface Wien2kStructureContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
}

interface PhononData {
  dos_data: any | null;
  dispersion_data: any | null;
}

type PhononViewMode = "bands" | "dos";
type PhononFrequencyUnit = "cm-1" | "thz";
type PhononBandFocus = "full" | "acoustic" | "optical";
const CM1_TO_THZ = 0.0299792458;

function formatEpochSeconds(epochSeconds: number): string {
  if (!Number.isFinite(epochSeconds) || epochSeconds <= 0) {
    return "unknown time";
  }
  return new Date(epochSeconds * 1000).toLocaleString();
}

function formatRemotePhononCandidateLocation(location: string): string {
  if (location === "project_archive") {
    return "archive";
  }
  return "workspace";
}

function normalizePositiveIntInput(input: string, fallback: number): number {
  const parsed = Number.parseInt(input, 10);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    return fallback;
  }
  return parsed;
}

function derivePostw90PathFromWannier90Path(path: string | null | undefined): string {
  const trimmed = String(path || "").trim();
  if (!trimmed) return DEFAULT_POSTW90_PATH;
  if (!trimmed.includes("/") && !trimmed.startsWith("~")) {
    return "postw90.x";
  }
  const segments = trimmed.split("/");
  segments[segments.length - 1] = "postw90.x";
  const derived = segments.join("/");
  return derived.trim().length > 0 ? derived : DEFAULT_POSTW90_PATH;
}

function deriveBundledWannier90PathFromQeBinDir(path: string | null | undefined): string | null {
  const trimmed = String(path || "").trim();
  if (!trimmed) return null;

  const normalized = trimmed.replace(/\\/g, "/").replace(/\/+$/, "");
  if (!normalized) return null;

  const segments = normalized.split("/");
  if (segments.length === 0) return null;
  if (segments[segments.length - 1] === "bin") {
    segments.pop();
  }
  const qeRoot = segments.join("/").trim();
  if (!qeRoot) return null;
  return `${qeRoot}/external/wannier90/wannier90.x`;
}

function resolveDefaultWannier90Path(qeBinDir: string | null | undefined): string {
  return deriveBundledWannier90PathFromQeBinDir(qeBinDir) || DEFAULT_WANNIER90_PATH;
}

function normalizeQeSmearing(raw: unknown): QeSmearingType {
  const lowered = String(raw || "").toLowerCase();
  if (lowered === "gaussian") return "gaussian";
  if (lowered === "methfessel-paxton") return "methfessel-paxton";
  if (lowered === "fermi-dirac") return "fermi-dirac";
  return "marzari-vanderbilt";
}

function cloneResourceDefaults(
  resourceType: "cpu" | "gpu",
  resources?: SlurmResourceRequest | null,
): SlurmResourceRequest {
  const fallback = resourceType === "gpu" ? defaultGpuResources() : defaultCpuResources();
  const merged = resources ? { ...fallback, ...resources } : fallback;
  return {
    ...merged,
    resource_type: resourceType,
    additional_sbatch: [...(merged.additional_sbatch || [])],
  };
}

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
  const taskContext = useTaskContext();
  const windowSize = useWindowSize();
  const plotHeight = Math.max(400, windowSize.height - 160);

  const [qePath, setQePath] = useState<string | null>(null);
  const [qePathInput, setQePathInput] = useState("");
  const [fermiSurferPathInput, setFermiSurferPathInput] = useState(DEFAULT_FERMI_SURFER_PATH);
  const [wannier90PathInput, setWannier90PathInput] = useState(() => resolveDefaultWannier90Path(null));
  const [isSavingQePath, setIsSavingQePath] = useState(false);
  const [isSavingFermiSurferPath, setIsSavingFermiSurferPath] = useState(false);
  const [isSavingWannier90Path, setIsSavingWannier90Path] = useState(false);
  const [isSavingQeDefaults, setIsSavingQeDefaults] = useState(false);
  const [availableExecutables, setAvailableExecutables] = useState<string[]>([]);
  const [qeStatus, setQeStatus] = useState<"Found" | "Not configured" | "Not found">("Not configured");
  const [fermiSurferStatus, setFermiSurferStatus] = useState<"Found" | "Not configured" | "Not found">("Not configured");
  const [wannier90Status, setWannier90Status] = useState<"Found" | "Not configured" | "Not found">("Not configured");
  const [error, setError] = useState<string | null>(null);
  const [currentView, setCurrentView] = useState<AppView>("project-browser");
  const [activeWorkflowRoute, setActiveWorkflowRoute] = useState<EngineWorkflowHostRoute | null>(null);
  const [selectedProjectId, setSelectedProjectId] = useState<string | null>(null);
  const [projectBrowserFolderId, setProjectBrowserFolderId] = useState<string | null>(null);
  const [bandsMultiviewInitialCalculations, setBandsMultiviewInitialCalculations] =
    useState<BandsMultiviewCalculation[] | undefined>(undefined);
  const [showCloseConfirm, setShowCloseConfirm] = useState(false);
  const [showNavigationDrawer, setShowNavigationDrawer] = useState(false);
  const [showTaskDrawer, setShowTaskDrawer] = useState(false);
  const [taskDrawerFilter, setTaskDrawerFilter] = useState<TaskManagerFilter>("active");
  const [taskDrawerFocusId, setTaskDrawerFocusId] = useState<string | null>(null);
  const [lastNonUtilityView, setLastNonUtilityView] = useState<AppView>("project-browser");
  const [showSettingsMenu, setShowSettingsMenu] = useState(false);
  const [settingsPage, setSettingsPage] = useState<"general" | "hpc" | "hpc-profile">("general");
  const [hpcProfileEditorDirty, setHpcProfileEditorDirty] = useState(false);
  const settingsMenuRef = useRef<HTMLDivElement | null>(null);
  const [executionMode, setExecutionMode] = useState<ExecutionMode>("local");
  const [hpcProfiles, setHpcProfiles] = useState<HpcProfile[]>([]);
  const [activeHpcProfileId, setActiveHpcProfileId] = useState<string | null>(null);
  const [engineInstallations, setEngineInstallations] = useState<EngineInstallation[]>([]);
  const [showHpcSetupWizard, setShowHpcSetupWizard] = useState(false);
  const [editingHpcProfileId, setEditingHpcProfileId] = useState<string | null>(null);
  const [hpcStatus, setHpcStatus] = useState<string | null>(null);
  const [qeDefaults, setQeDefaults] = useState<QeDefaults>(DEFAULT_QE_DEFAULTS);
  const [qeDefaultsStatus, setQeDefaultsStatus] = useState<string | null>(null);
  const [isExportingHpcPresetBundle, setIsExportingHpcPresetBundle] = useState(false);
  const [isImportingHpcPresetBundle, setIsImportingHpcPresetBundle] = useState(false);
  const [isCleaningHpcRemote, setIsCleaningHpcRemote] = useState(false);
  const [isMigratingHpcRoots, setIsMigratingHpcRoots] = useState(false);
  const [isRepairingLocalPseudos, setIsRepairingLocalPseudos] = useState(false);
  const [isRepairingRemotePseudos, setIsRepairingRemotePseudos] = useState(false);
  const [pseudoRepairStatus, setPseudoRepairStatus] = useState<string | null>(null);
  const [showHeadlessRecoveryDialog, setShowHeadlessRecoveryDialog] = useState(false);
  const [headlessCandidates, setHeadlessCandidates] = useState<HpcHeadlessJobCandidate[]>([]);
  const [isLoadingHeadlessCandidates, setIsLoadingHeadlessCandidates] = useState(false);
  const [isAttachingHeadlessJob, setIsAttachingHeadlessJob] = useState(false);
  const [selectedHeadlessJobId, setSelectedHeadlessJobId] = useState<string | null>(null);
  const [headlessRecoveryStatus, setHeadlessRecoveryStatus] = useState<string | null>(null);
  const [hpcDefaultCpuDraft, setHpcDefaultCpuDraft] = useState<SlurmResourceRequest>(
    defaultCpuResources(),
  );
  const [hpcDefaultGpuDraft, setHpcDefaultGpuDraft] = useState<SlurmResourceRequest>(
    defaultGpuResources(),
  );
  const [hpcResourceModeDraft, setHpcResourceModeDraft] = useState<HpcResourceMode>("both");
  const [hpcLauncherDraft, setHpcLauncherDraft] = useState<HpcLauncher>("srun");
  const [hpcLauncherCpuExtraArgsDraft, setHpcLauncherCpuExtraArgsDraft] = useState("");
  const [hpcLauncherGpuExtraArgsDraft, setHpcLauncherGpuExtraArgsDraft] = useState("");
  const [isSavingHpcDefaults, setIsSavingHpcDefaults] = useState(false);
  const [hpcDefaultsSaved, setHpcDefaultsSaved] = useState(false);
  const [hpcDefaultsStatus, setHpcDefaultsStatus] = useState<string | null>(null);
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
  const [isLoadingRemotePhononCandidates, setIsLoadingRemotePhononCandidates] = useState(false);
  const [isRecoveringRemotePhonon, setIsRecoveringRemotePhonon] = useState(false);
  const [showRemotePhononSelectionDialog, setShowRemotePhononSelectionDialog] = useState(false);
  const [remotePhononRecoveryContext, setRemotePhononRecoveryContext] = useState<RemotePhononRecoveryContext | null>(null);
  const [remotePhononCandidates, setRemotePhononCandidates] = useState<HpcRecoverableRemotePhononRun[]>([]);
  const [selectedRemotePhononWorkdir, setSelectedRemotePhononWorkdir] = useState("");
  const [recoveryStatus, setRecoveryStatus] = useState<string | null>(null);
  const [showDeleteProjectDialog, setShowDeleteProjectDialog] = useState(false);
  const [deleteConfirmText, setDeleteConfirmText] = useState("");
  const [isDeletingProject, setIsDeletingProject] = useState(false);
  const [deleteProjectSnapshot, setDeleteProjectSnapshot] = useState<SettingsProjectSnapshot | null>(null);
  const [showCleanRemoteConfirmDialog, setShowCleanRemoteConfirmDialog] = useState(false);
  const [showMigrateHpcRootsDialog, setShowMigrateHpcRootsDialog] = useState(false);
  const [migrateWorkspaceRootDraft, setMigrateWorkspaceRootDraft] = useState("");
  const [migrateProjectRootDraft, setMigrateProjectRootDraft] = useState("");
  const [projectDashboardRefreshToken, setProjectDashboardRefreshToken] = useState(0);

  function openEngineWorkflow(engineId: EngineId, kind: CalculationKind, fallbackView: EngineWorkflowView) {
    const route = resolveEngineWorkflowHostRoute(engineId, kind);
    if (!route) {
      console.warn(`No frontend workflow route registered for ${engineId}:${kind}`);
      setActiveWorkflowRoute(null);
      setCurrentView(fallbackView);
      return;
    }
    setActiveWorkflowRoute(route);
    setCurrentView(route.view);
  }

  // Active task ID for reconnection when navigating to wizard from indicator
  const [reconnectTaskId, setReconnectTaskId] = useState<string | null>(null);

  // Context for running SCF from a project
  const [scfContext, setScfContext] = useState<SCFContext | null>(null);

  // Context for running Bands from a project
  const [bandsContext, setBandsContext] = useState<BandsContext | null>(null);

  // Context for viewing saved band data
  const [viewBandsData, setViewBandsData] = useState<{
    bandData: BandPlotData;
    fermiEnergy: number | null;
    calculationParameters?: Record<string, unknown> | null;
    calculationContext?: SavedBandsCalculationContext | null;
  } | null>(null);

  // Context for running Electronic DOS from a project
  const [dosContext, setDosContext] = useState<DosContext | null>(null);

  // Context for viewing saved DOS data
  const [viewDosData, setViewDosData] = useState<{ dosData: ElectronicDOSData; fermiEnergy: number | null } | null>(null);

  // Context for running Wannier90 from a project
  const [wannierContext, setWannierContext] = useState<WannierContext | null>(null);

  // Context for viewing saved Wannier90 data
  const [viewWannierData, setViewWannierData] = useState<{
    result: any;
    fermiEnergy: number | null;
    overlayOptions: WannierBandOverlayOption[];
  } | null>(null);

  // Context for running Fermi-surface generation from a project
  const [fermiSurfaceContext, setFermiSurfaceContext] = useState<FermiSurfaceContext | null>(null);

  // Context for running transport from a saved Wannier calculation
  const [transportContext, setTransportContext] = useState<TransportContext | null>(null);

  // Context for viewing saved transport data
  const [viewTransportData, setViewTransportData] = useState<{ data: TransportResult } | null>(null);

  // Context for running Phonons from a project
  const [phononsContext, setPhononsContext] = useState<PhononsContext | null>(null);

  // Context for running Hubbard linear response from a project
  const [hubbardLrtContext, setHubbardLrtContext] = useState<HubbardLrtContext | null>(null);

  // Context for running EPW from a project
  const [epwContext, setEpwContext] = useState<EpwContext | null>(null);

  // Context for preparing a WIEN2k structure source from a project CIF
  const [wien2kStructureContext, setWien2kStructureContext] = useState<Wien2kStructureContext | null>(null);

  // Context for viewing saved phonon data
  const [viewPhononData, setViewPhononData] = useState<{ data: PhononData; mode: PhononViewMode } | null>(null);
  const [viewEpwData, setViewEpwData] = useState<EpwViewerPayload | null>(null);
  const [phononBandsUnit, setPhononBandsUnit] = useState<PhononFrequencyUnit>("cm-1");
  const [phononBandFocus, setPhononBandFocus] = useState<PhononBandFocus>("full");
  const isHpcActivityPopout = new URLSearchParams(window.location.search).get("hpc_activity") === "1";
  const activeHpcProfile = useMemo(
    () => hpcProfiles.find((profile) => profile.id === activeHpcProfileId) ?? null,
    [hpcProfiles, activeHpcProfileId],
  );
  const taskManagerEntries = useMemo(
    () => buildTaskManagerEntries(taskContext.tasks.values(), taskContext.queueItems),
    [taskContext.queueItems, taskContext.tasks],
  );
  const taskManagerSummary = useMemo(
    () => summarizeTaskManagerEntries(taskManagerEntries),
    [taskManagerEntries],
  );
  const relevantTaskManagerEntry = useMemo(
    () => findRelevantTaskManagerEntry(taskManagerEntries),
    [taskManagerEntries],
  );

  useEffect(() => {
    const hasTopBar = !isHpcActivityPopout;
    document.documentElement.setAttribute("data-app-chrome", hasTopBar ? "topbar" : "none");
    return () => {
      document.documentElement.removeAttribute("data-app-chrome");
    };
  }, [currentView, isHpcActivityPopout]);
  const editingHpcProfile = useMemo(
    () => hpcProfiles.find((profile) => profile.id === editingHpcProfileId) ?? null,
    [hpcProfiles, editingHpcProfileId],
  );

  function confirmDiscardHpcProfileChanges(): boolean {
    return settingsPage !== "hpc-profile"
      || !hpcProfileEditorDirty
      || window.confirm("Discard unsaved changes to this HPC profile?");
  }

  function closeSettingsMenu() {
    if (!confirmDiscardHpcProfileChanges()) return false;
    setHpcProfileEditorDirty(false);
    setShowSettingsMenu(false);
    if (currentView === "settings") {
      returnFromUtilityView();
    }
    return true;
  }

  function leaveHpcProfileEditor() {
    if (!confirmDiscardHpcProfileChanges()) return;
    setHpcProfileEditorDirty(false);
    setSettingsPage("hpc");
  }

  function openNavigationDrawer() {
    if (windowSize.width <= 900) setShowTaskDrawer(false);
    setShowNavigationDrawer(true);
  }

  function openTaskDrawer(filter: TaskManagerFilter = "active", taskId: string | null = null) {
    if (windowSize.width <= 900) setShowNavigationDrawer(false);
    const focusedEntry = taskId
      ? taskManagerEntries.find((entry) => entry.taskId === taskId)
      : null;
    setTaskDrawerFilter(filter === "active" && focusedEntry?.group === "finished" ? "finished" : filter);
    setTaskDrawerFocusId(taskId);
    setShowTaskDrawer(true);
  }

  function openSettingsWorkspace() {
    setShowNavigationDrawer(false);
    if (currentView !== "settings") setLastNonUtilityView(currentView);
    setShowSettingsMenu(true);
    setCurrentView("settings");
  }

  useEffect(() => {
    if (!activeHpcProfile) {
      setHpcDefaultCpuDraft(defaultCpuResources());
      setHpcDefaultGpuDraft(defaultGpuResources());
      setHpcResourceModeDraft("both");
      setHpcLauncherDraft("srun");
      setHpcLauncherCpuExtraArgsDraft("");
      setHpcLauncherGpuExtraArgsDraft("");
      clearHpcDefaultsSaveFeedback();
      return;
    }
    setHpcDefaultCpuDraft(cloneResourceDefaults("cpu", activeHpcProfile.default_cpu_resources));
    setHpcDefaultGpuDraft(cloneResourceDefaults("gpu", activeHpcProfile.default_gpu_resources));
    setHpcResourceModeDraft(activeHpcProfile.resource_mode ?? "both");
    setHpcLauncherDraft(activeHpcProfile.launcher ?? "srun");
    setHpcLauncherCpuExtraArgsDraft(
      activeHpcProfile.launcher_cpu_extra_args
        ?? activeHpcProfile.launcher_extra_args
        ?? "",
    );
    setHpcLauncherGpuExtraArgsDraft(
      activeHpcProfile.launcher_gpu_extra_args
        ?? activeHpcProfile.launcher_extra_args
        ?? "",
    );
    clearHpcDefaultsSaveFeedback();
  }, [activeHpcProfile?.id, activeHpcProfile?.updated_at]);

  function clearHpcDefaultsSaveFeedback() {
    setHpcDefaultsSaved(false);
    setHpcDefaultsStatus(null);
  }

  // Check for existing QE configuration on startup
  useEffect(() => {
    checkQEPath();
    void loadQeDefaults();
    void loadFermiSurferPath();
    void loadWannier90Path();
    void loadExecutionPrefix();
    void loadHpcExecutionSettings();
    void loadGlobalMpiSettings();
    void loadGlobalSaveSizeSetting();
  }, []);

  // Listen for explicit quit confirmation events from backend
  useEffect(() => {
    const unlisten = listen("confirm-quit", () => {
      setShowCloseConfirm(true);
    });
    return () => { unlisten.then((fn) => fn()); };
  }, []);

  useEffect(() => {
    if (executionMode !== "hpc" && (settingsPage === "hpc" || settingsPage === "hpc-profile")) {
      setSettingsPage("general");
    }
  }, [executionMode, settingsPage]);

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
      if (wannier90Status === "Not configured") {
        setWannier90PathInput(resolveDefaultWannier90Path(normalized));
      }
      await loadExecutables();
      setError(null);
    } catch (e) {
      setError(String(e));
    } finally {
      setIsSavingQePath(false);
    }
  }

  function getLocalPseudoDir() {
    const basePath = (qePath || qePathInput).trim();
    if (!basePath) return "";
    return basePath.replace(/\/bin\/?$/, "/pseudo");
  }

  function formatPseudoRepairResult(scope: "Local" | "Remote", result: PslibraryPseudoRepairResult) {
    const changed = result.patched_files.slice(0, 4).join(", ");
    const suffix = result.patched_files.length > 4
      ? `, +${result.patched_files.length - 4} more`
      : "";
    const patchedDetail = result.patched > 0
      ? ` Patched: ${changed}${suffix}.`
      : " No changes needed.";
    return `${scope} PSLibrary repair scanned ${result.scanned} UPFs, found ${result.candidates} affected candidates, patched ${result.patched}, already clean ${result.already_clean}.${patchedDetail}`;
  }

  async function repairLocalPslibraryPseudos() {
    const pseudoDir = getLocalPseudoDir();
    if (!pseudoDir) {
      setPseudoRepairStatus("Set the local QE bin path first.");
      return;
    }

    setIsRepairingLocalPseudos(true);
    setPseudoRepairStatus(null);
    try {
      const result = await invoke<PslibraryPseudoRepairResult>("repair_local_pslibrary_pseudopotentials", {
        pseudoDir,
      });
      setPseudoRepairStatus(formatPseudoRepairResult("Local", result));
      setError(null);
    } catch (e) {
      const message = String(e);
      setPseudoRepairStatus(message);
      setError(message);
    } finally {
      setIsRepairingLocalPseudos(false);
    }
  }

  async function repairRemotePslibraryPseudos() {
    if (!activeHpcProfile) {
      setPseudoRepairStatus("Select an active HPC profile first.");
      return;
    }

    setIsRepairingRemotePseudos(true);
    setPseudoRepairStatus(null);
    try {
      const result = await invoke<PslibraryPseudoRepairResult>("hpc_repair_remote_pslibrary_pseudopotentials", {
        profileId: activeHpcProfile.id,
        pseudoDir: resolveProfileRemotePseudoDir(activeHpcProfile),
      });
      setPseudoRepairStatus(formatPseudoRepairResult("Remote", result));
      setError(null);
    } catch (e) {
      const message = String(e);
      setPseudoRepairStatus(message);
      setError(message);
    } finally {
      setIsRepairingRemotePseudos(false);
    }
  }

  async function loadQeDefaults() {
    try {
      const defaults = await invoke<Partial<QeDefaults> | null>("get_qe_defaults");
      setQeDefaults({
        smearing: normalizeQeSmearing(defaults?.smearing),
      });
      setQeDefaultsStatus(null);
    } catch (e) {
      console.error("Failed to load QE defaults:", e);
      setQeDefaults(DEFAULT_QE_DEFAULTS);
      setQeDefaultsStatus("Failed to load QE defaults");
    }
  }

  async function saveQeDefaults() {
    setIsSavingQeDefaults(true);
    setQeDefaultsStatus(null);
    try {
      const normalized: QeDefaults = {
        smearing: normalizeQeSmearing(qeDefaults.smearing),
      };
      await invoke("set_qe_defaults", { defaults: normalized });
      setQeDefaults(normalized);
      setQeDefaultsStatus("Saved");
      setError(null);
    } catch (e) {
      setQeDefaultsStatus("Failed to save QE defaults");
      setError(String(e));
    } finally {
      setIsSavingQeDefaults(false);
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

  async function loadWannier90Path() {
    try {
      const [path, configuredQePath] = await Promise.all([
        invoke<string | null>("get_wannier90_path"),
        invoke<string | null>("get_qe_path"),
      ]);
      const defaultWannier90Path = resolveDefaultWannier90Path(configuredQePath);
      if (path) {
        setWannier90PathInput(path);
        setWannier90Status("Found");
      } else {
        setWannier90PathInput(defaultWannier90Path);
        setWannier90Status("Not configured");
      }
    } catch (e) {
      console.error("Failed to load Wannier90 path:", e);
      setWannier90PathInput(resolveDefaultWannier90Path(qePath));
      setWannier90Status("Not found");
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

  async function selectWannier90Path() {
    try {
      const selected = await open({
        directory: false,
        multiple: false,
        defaultPath: wannier90PathInput || qePath || "/usr/local/bin",
        title: "Select Wannier90 executable",
      });

      if (selected && typeof selected === "string") {
        setWannier90PathInput(selected);
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

  async function saveWannier90Path() {
    const normalized = wannier90PathInput.trim();
    const defaultWannier90Path = resolveDefaultWannier90Path(qePath);
    setIsSavingWannier90Path(true);
    try {
      await invoke("set_wannier90_path", {
        path: normalized.length > 0 ? normalized : null,
      });
      if (normalized.length > 0) {
        setWannier90PathInput(normalized);
        setWannier90Status("Found");
      } else {
        setWannier90PathInput(defaultWannier90Path);
        setWannier90Status("Not configured");
      }
      setError(null);
    } catch (e) {
      setWannier90Status("Not found");
      setError(String(e));
    } finally {
      setIsSavingWannier90Path(false);
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

  async function loadHpcExecutionSettings() {
    try {
      const [mode, profiles, activeProfileId, installations] = await Promise.all([
        loadExecutionMode(),
        listHpcProfiles(),
        getActiveHpcProfileId(),
        listEngineInstallations(),
      ]);
      setExecutionMode(mode);
      setHpcProfiles(profiles);
      setActiveHpcProfileId(activeProfileId);
      setEngineInstallations(installations);
    } catch (e) {
      console.error("Failed to load HPC settings:", e);
      setExecutionMode("local");
      setHpcProfiles([]);
      setActiveHpcProfileId(null);
      setEngineInstallations([]);
    }
  }

  async function handleExecutionModeChange(mode: ExecutionMode) {
    setHpcStatus(null);
    try {
      await saveExecutionMode(mode);
      setExecutionMode(mode);
      if (mode === "hpc" && hpcProfiles.length === 0) {
        setShowHpcSetupWizard(true);
      }
      setHpcStatus("Execution mode saved.");
    } catch (e) {
      console.error("Failed to save execution mode:", e);
      setHpcStatus(`Failed to save execution mode: ${e}`);
    }
  }

  async function handleSelectHpcProfile(profileId: string) {
    setHpcStatus(null);
    try {
      await setActiveHpcProfile(profileId);
      setActiveHpcProfileId(profileId);
      setHpcStatus("Active profile updated.");
    } catch (e) {
      console.error("Failed to set active HPC profile:", e);
      setHpcStatus(`Failed to set active profile: ${e}`);
    }
  }

  async function handleDeleteHpcProfile(profileId: string) {
    setHpcStatus(null);
    try {
      await deleteHpcProfile(profileId);
      await loadHpcExecutionSettings();
      setHpcStatus("Profile deleted.");
    } catch (e) {
      console.error("Failed to delete HPC profile:", e);
      setHpcStatus(`Failed to delete profile: ${e}`);
    }
  }

  async function handleExportHpcPresetBundle() {
    if (hpcProfiles.length === 0) {
      setHpcStatus("No HPC profiles to export.");
      return;
    }

    setIsExportingHpcPresetBundle(true);
    setHpcStatus(null);
    try {
      const destinationPath = await save({
        title: "Export HPC Presets + Defaults",
        defaultPath: `qcortado-hpc-presets-${new Date().toISOString().slice(0, 10)}.qchpc`,
        filters: [{ name: "QCortado HPC Preset Bundle", extensions: ["qchpc", "json"] }],
      });
      if (!destinationPath) {
        return;
      }

      const result = await exportHpcPresetBundle(destinationPath);
      setHpcStatus(
        `Exported ${result.profile_count} profile preset(s) to ${result.bundle_path}. Usernames and credentials were excluded.`,
      );
    } catch (e) {
      console.error("Failed to export HPC preset bundle:", e);
      setHpcStatus(`Failed to export presets: ${e}`);
    } finally {
      setIsExportingHpcPresetBundle(false);
    }
  }

  async function handleImportHpcPresetBundle() {
    setIsImportingHpcPresetBundle(true);
    setHpcStatus(null);
    try {
      const selectedPath = await open({
        title: "Import HPC Presets + Defaults",
        directory: false,
        multiple: false,
        filters: [{ name: "QCortado HPC Preset Bundle", extensions: ["qchpc", "json"] }],
      });
      if (!selectedPath || Array.isArray(selectedPath)) {
        return;
      }

      const result = await importHpcPresetBundle(selectedPath);
      await loadHpcExecutionSettings();

      const summary = `Imported ${result.imported_profile_count} profile preset(s): ${result.updated_profile_count} updated, ${result.created_profile_count} created.`;
      if (result.profiles_requiring_username.length > 0) {
        setHpcStatus(
          `${summary} Set usernames before connecting for: ${result.profiles_requiring_username.join(", ")}.`,
        );
      } else {
        setHpcStatus(summary);
      }
    } catch (e) {
      console.error("Failed to import HPC preset bundle:", e);
      setHpcStatus(`Failed to import presets: ${e}`);
    } finally {
      setIsImportingHpcPresetBundle(false);
    }
  }

  function openCleanRemoteConfirmDialog() {
    if (!activeHpcProfile) {
      setHpcStatus("Select an active HPC profile first.");
      return;
    }
    setShowCleanRemoteConfirmDialog(true);
  }

  async function handleConfirmCleanRemoteOrphans() {
    if (!activeHpcProfile) {
      setShowCleanRemoteConfirmDialog(false);
      setHpcStatus("Select an active HPC profile first.");
      return;
    }

    setIsCleaningHpcRemote(true);
    setHpcStatus(null);
    try {
      const result = await cleanHpcRemoteOrphans(activeHpcProfile.id);
      if (result.failed_paths.length > 0) {
        setHpcStatus(
          `Removed ${result.removed_paths.length} orphan path(s), but ${result.failed_paths.length} could not be removed.`,
        );
      } else if (result.removed_paths.length > 0) {
        setHpcStatus(
          `Removed ${result.removed_paths.length} orphan remote path(s) after scanning ${result.scanned_paths}.`,
        );
      } else {
        setHpcStatus(`No orphaned QCortado remote paths found (${result.scanned_paths} scanned).`);
      }
      setShowCleanRemoteConfirmDialog(false);
    } catch (e) {
      console.error("Failed to clean remote HPC artifacts:", e);
      setHpcStatus(`Failed to clean remote artifacts: ${e}`);
    } finally {
      setIsCleaningHpcRemote(false);
    }
  }

  function openMigrateHpcRootsDialog() {
    if (!activeHpcProfile) {
      setHpcStatus("Select an active HPC profile first.");
      return;
    }
    setMigrateWorkspaceRootDraft(activeHpcProfile.remote_workspace_root || "");
    setMigrateProjectRootDraft(activeHpcProfile.remote_project_root || "");
    setShowMigrateHpcRootsDialog(true);
  }

  async function handleConfirmMigrateHpcRoots() {
    if (!activeHpcProfile) {
      setHpcStatus("Select an active HPC profile first.");
      setShowMigrateHpcRootsDialog(false);
      return;
    }

    const newWorkspaceRoot = migrateWorkspaceRootDraft.trim();
    const newProjectRoot = migrateProjectRootDraft.trim();
    if (!newWorkspaceRoot || !newProjectRoot) {
      setHpcStatus("Both new remote roots are required.");
      return;
    }

    if (
      newWorkspaceRoot === activeHpcProfile.remote_workspace_root
      && newProjectRoot === activeHpcProfile.remote_project_root
    ) {
      setHpcStatus("Remote roots are unchanged.");
      return;
    }

    setIsMigratingHpcRoots(true);
    setHpcStatus(null);
    try {
      const updated = await migrateHpcRemoteRoots(
        activeHpcProfile.id,
        newWorkspaceRoot,
        newProjectRoot,
      );
      setHpcProfiles((prev) => prev.map((profile) => (profile.id === updated.id ? updated : profile)));
      setHpcStatus("Remote roots migrated and profile updated.");
      setShowMigrateHpcRootsDialog(false);
    } catch (e) {
      console.error("Failed to migrate HPC remote roots:", e);
      setHpcStatus(`Failed to migrate remote roots: ${e}`);
    } finally {
      setIsMigratingHpcRoots(false);
    }
  }

  async function openHeadlessRecoveryDialog() {
    if (!activeHpcProfile) {
      setHpcStatus("Select an active HPC profile first.");
      return;
    }
    setShowNavigationDrawer(false);
    setShowHeadlessRecoveryDialog(true);
    setHeadlessRecoveryStatus(null);
    setIsLoadingHeadlessCandidates(true);
    try {
      const candidates = await listHeadlessHpcJobs(activeHpcProfile.id, 50);
      setHeadlessCandidates(candidates);
      setSelectedHeadlessJobId(candidates[0]?.remote_job_id ?? null);
      setHeadlessRecoveryStatus(
        candidates.length > 0
          ? `Found ${candidates.length} recoverable QCortado job(s).`
          : "No recoverable QCortado jobs were found for the active profile.",
      );
    } catch (e) {
      console.error("Failed to list headless HPC jobs:", e);
      setHeadlessCandidates([]);
      setSelectedHeadlessJobId(null);
      setHeadlessRecoveryStatus(`Failed to query headless jobs: ${e}`);
    } finally {
      setIsLoadingHeadlessCandidates(false);
    }
  }

  async function handleAttachHeadlessJob() {
    const candidate = headlessCandidates.find((item) => item.remote_job_id === selectedHeadlessJobId);
    if (!candidate) {
      setHeadlessRecoveryStatus("Select a job to recover.");
      return;
    }
    setIsAttachingHeadlessJob(true);
    setHeadlessRecoveryStatus(null);
    try {
      const result = await attachHeadlessHpcJob(
        candidate.remote_job_id,
        candidate.remote_workdir ?? null,
        candidate.profile_id,
      );
      await taskContext.reconnectToTask(result.task_id);
      setShowHeadlessRecoveryDialog(false);
      openTaskDrawer("hpc", result.task_id);
    } catch (e) {
      console.error("Failed to attach headless HPC job:", e);
      setHeadlessRecoveryStatus(`Failed to attach job: ${e}`);
    } finally {
      setIsAttachingHeadlessJob(false);
    }
  }

  async function saveHpcDefaultRunSettings() {
    if (!activeHpcProfile) {
      setHpcDefaultsStatus("Select an active profile first.");
      return;
    }

    setIsSavingHpcDefaults(true);
    clearHpcDefaultsSaveFeedback();
    try {
      const saved = await updateHpcProfileDefaults(
        activeHpcProfile.id,
        cloneResourceDefaults("cpu", hpcDefaultCpuDraft),
        cloneResourceDefaults("gpu", hpcDefaultGpuDraft),
        hpcResourceModeDraft,
        hpcLauncherDraft,
        hpcLauncherCpuExtraArgsDraft.trim() || null,
        hpcLauncherGpuExtraArgsDraft.trim() || null,
      );
      setHpcProfiles((prev) => prev.map((profile) => (profile.id === saved.id ? saved : profile)));
      setHpcDefaultCpuDraft(cloneResourceDefaults("cpu", saved.default_cpu_resources));
      setHpcDefaultGpuDraft(cloneResourceDefaults("gpu", saved.default_gpu_resources));
      setHpcResourceModeDraft(saved.resource_mode ?? "both");
      setHpcLauncherDraft(saved.launcher ?? "srun");
      setHpcLauncherCpuExtraArgsDraft(saved.launcher_cpu_extra_args || "");
      setHpcLauncherGpuExtraArgsDraft(saved.launcher_gpu_extra_args || "");
      setHpcDefaultsSaved(true);
      setHpcDefaultsStatus(null);
    } catch (e) {
      console.error("Failed to save HPC defaults:", e);
      setHpcDefaultsSaved(false);
      setHpcDefaultsStatus(`Failed to save defaults: ${e}`);
    } finally {
      setIsSavingHpcDefaults(false);
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

  async function recoverRemotePhononFromSettings() {
    if (!selectedProjectId) {
      setRecoveryStatus("Open a project to recover remote phonon data.");
      return;
    }

    setRecoveryStatus(null);
    setShowSettingsMenu(false);
    returnFromUtilityView();
    closeRemotePhononSelectionDialog(true);
    setShowRemotePhononSelectionDialog(true);
    setIsLoadingRemotePhononCandidates(true);
    try {
      const project = await invoke<SettingsProjectSnapshot>("get_project", { projectId: selectedProjectId });
      const selectedVariant = (project.last_opened_cif_id && project.cif_variants.some((variant) => variant.id === project.last_opened_cif_id))
        ? project.cif_variants.find((variant) => variant.id === project.last_opened_cif_id) ?? null
        : project.cif_variants[0] ?? null;
      if (!selectedVariant) {
        setRecoveryStatus("No structure found in this project.");
        closeRemotePhononSelectionDialog(true);
        return;
      }

      const fallbackScf = selectedVariant.calculations
        .filter((calc) => calc.calc_type === "scf" && calc.result?.converged)
        .sort((a, b) => {
          const aTime = a.completed_at ?? a.started_at;
          const bTime = b.completed_at ?? b.started_at;
          return bTime.localeCompare(aTime);
        })[0];

      const candidates = await invoke<HpcRecoverableRemotePhononRun[]>(
        "hpc_list_recoverable_remote_phonon_runs",
        { profileId: null, limit: 20 },
      );
      console.info("[remote-phonon-recovery] recoverable candidates", candidates);

      if (candidates.length === 0) {
        try {
          const debugReport = await invoke<HpcRemotePhononRecoveryDebugReport>(
            "hpc_debug_remote_phonon_recovery",
            { profileId: null },
          );
          console.info("[remote-phonon-recovery] debug probe", debugReport);
          const workspaceScanned = debugReport.workspace_probe_output
            .split(/\r?\n/)
            .filter((line) => line.trim().length > 0).length;
          const projectScanned = debugReport.project_probe_output
            .split(/\r?\n/)
            .filter((line) => line.trim().length > 0).length;
          setRecoveryStatus(
            `No recoverable remote phonon runs found. Debug scanned ${workspaceScanned} workspace and ${projectScanned} archive directories; see DevTools console for details.`,
          );
        } catch (debugError) {
          console.error("[remote-phonon-recovery] debug probe failed", debugError);
          setRecoveryStatus(
            `No recoverable remote phonon runs found. Debug probe failed: ${String(debugError)}`,
          );
        }
        setRemotePhononCandidates([]);
        setSelectedRemotePhononWorkdir("");
        setRemotePhononRecoveryContext(null);
        return;
      }

      const sortedCandidates = [...candidates].sort((a, b) => {
        if (b.modified_at_epoch !== a.modified_at_epoch) {
          return b.modified_at_epoch - a.modified_at_epoch;
        }
        return a.remote_workdir.localeCompare(b.remote_workdir);
      });
      const defaultCandidate = sortedCandidates[0];
      setRemotePhononCandidates(sortedCandidates);
      setSelectedRemotePhononWorkdir(defaultCandidate.remote_workdir);
      setRemotePhononRecoveryContext({
        projectId: selectedProjectId,
        cifId: selectedVariant.id,
        sourceScfId: fallbackScf?.id ?? null,
      });
      setRecoveryStatus(
        `Found ${sortedCandidates.length} recoverable remote phonon runs. Select a run by date/time.`,
      );
    } catch (e) {
      console.error("Failed to recover remote phonon calculation:", e);
      const message = e instanceof Error ? e.message : String(e);
      setRecoveryStatus(`Remote phonon recovery failed: ${message}`);
    } finally {
      setIsLoadingRemotePhononCandidates(false);
    }
  }

  function closeRemotePhononSelectionDialog(force = false) {
    if (isRecoveringRemotePhonon && !force) {
      return;
    }
    setIsLoadingRemotePhononCandidates(false);
    setShowRemotePhononSelectionDialog(false);
    setRemotePhononRecoveryContext(null);
    setRemotePhononCandidates([]);
    setSelectedRemotePhononWorkdir("");
  }

  async function confirmRemotePhononRecoverySelection() {
    if (isLoadingRemotePhononCandidates) {
      return;
    }
    if (!remotePhononRecoveryContext) {
      setRecoveryStatus("Remote recovery context is missing. Re-open recovery and select a run.");
      return;
    }
    if (!selectedRemotePhononWorkdir.trim()) {
      setRecoveryStatus("Select a remote phonon run to recover.");
      return;
    }

    setIsRecoveringRemotePhonon(true);
    try {
      await invoke("hpc_recover_remote_phonon_calculation", {
        projectId: remotePhononRecoveryContext.projectId,
        cifId: remotePhononRecoveryContext.cifId,
        remoteWorkdir: selectedRemotePhononWorkdir,
        profileId: null,
        sourceScfId: remotePhononRecoveryContext.sourceScfId,
      });
      setRecoveryStatus(`Recovered remote phonon calculation from ${selectedRemotePhononWorkdir}.`);
      setProjectDashboardRefreshToken((prev) => prev + 1);
      closeRemotePhononSelectionDialog(true);
    } catch (e) {
      console.error("Failed to recover selected remote phonon calculation:", e);
      const message = e instanceof Error ? e.message : String(e);
      setRecoveryStatus(`Remote phonon recovery failed: ${message}`);
    } finally {
      setIsRecoveringRemotePhonon(false);
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
      setHubbardLrtContext(null);
      setPhononsContext(null);
      setEpwContext(null);
      setViewBandsData(null);
      setViewDosData(null);
      setViewPhononData(null);
      setViewEpwData(null);
      setReconnectTaskId(null);
      setActiveWorkflowRoute(null);
      setSelectedProjectId(null);
      setCurrentView("project-browser");
    } catch (e) {
      console.error("Failed to delete project:", e);
      setPrefixStatus("Failed to delete project");
    } finally {
      setIsDeletingProject(false);
    }
  }

  function handleNavigateToTask(taskId: string, taskType: string) {
    setReconnectTaskId(taskId);
    const kindMap: Record<string, CalculationKind> = {
      scf: "scf",
      bands: "bands",
      wien2k_scf: "scf",
      wien2k_bands: "bands",
      wien2k_fermi_surface: "fermi_surface",
      dos: "dos",
      wannier: "wannier",
      transport: "transport",
      fermi_surface: "fermi_surface",
      hubbard_lrt: "hubbard_lrt",
      phonon: "phonon",
      epw: "epw",
    };
    const fallbackViewMap: Record<string, EngineWorkflowView> = {
      scf: "scf-wizard",
      bands: "bands-wizard",
      wien2k_scf: "wien2k-scf-wizard",
      wien2k_bands: "bands-wizard",
      wien2k_fermi_surface: "fermi-surface-wizard",
      dos: "dos-wizard",
      wannier: "wannier-wizard",
      transport: "transport-wizard",
      fermi_surface: "fermi-surface-wizard",
      hubbard_lrt: "hubbard-lrt-wizard",
      phonon: "phonon-wizard",
      epw: "epw-wizard",
    };
    const kind = kindMap[taskType];
    const fallbackView = fallbackViewMap[taskType];
    if (kind && fallbackView) {
      openEngineWorkflow(taskType.startsWith("wien2k_") ? "wien2k" : "qe", kind, fallbackView);
    }
  }

  function clearWorkflowContext(view: EngineWorkflowView) {
    switch (view) {
      case "scf-wizard":
        setScfContext(null);
        break;
      case "bands-wizard":
        setBandsContext(null);
        break;
      case "dos-wizard":
        setDosContext(null);
        break;
      case "wannier-wizard":
        setWannierContext(null);
        break;
      case "transport-wizard":
        setTransportContext(null);
        break;
      case "fermi-surface-wizard":
        setFermiSurfaceContext(null);
        break;
      case "hubbard-lrt-wizard":
        setHubbardLrtContext(null);
        break;
      case "phonon-wizard":
        setPhononsContext(null);
        break;
      case "epw-wizard":
        setEpwContext(null);
        break;
      case "wien2k-structure-wizard":
        setWien2kStructureContext(null);
        break;
      case "wien2k-scf-wizard":
        setScfContext(null);
        break;
    }
  }

  function handleWorkflowBack(view: EngineWorkflowView, destination: EngineWorkflowBackDestination) {
    if (activeWorkflowRoute?.engineId !== "wien2k") {
      clearWorkflowContext(view);
    }
    setReconnectTaskId(null);
    setActiveWorkflowRoute(null);
    setCurrentView(destination);
  }

  useEffect(() => {
    if (isHpcActivityPopout) return;
    let cancelled = false;
    let cleanup: (() => void) | null = null;

    listen<OpenTaskViewRequest>("hpc-open-task-view", async (event) => {
      const { taskId, taskType } = event.payload;
      if (!taskId || !taskType) return;
      handleNavigateToTask(taskId, taskType);

      try {
        const currentWindow = getCurrentWindow();
        await currentWindow.show();
        await currentWindow.unminimize();
        await currentWindow.setFocus();
      } catch (e) {
        console.error("Failed to focus main window for task view:", e);
      }
    }).then((unlisten) => {
      if (cancelled) {
        unlisten();
      } else {
        cleanup = unlisten;
      }
    });

    return () => {
      cancelled = true;
      if (cleanup) cleanup();
    };
  }, [isHpcActivityPopout]);

  function navigateToProjects() {
    setShowNavigationDrawer(false);
    setCurrentView("project-browser");
  }

  function navigateToNodeActivity() {
    setShowNavigationDrawer(false);
    if (currentView !== "node-activity") {
      setLastNonUtilityView(currentView);
      setCurrentView("node-activity");
    }
  }

  function navigateToStorageManager() {
    setShowNavigationDrawer(false);
    if (currentView !== "storage-manager") {
      setLastNonUtilityView(currentView);
      setCurrentView("storage-manager");
    }
  }

  function returnFromUtilityView() {
    const fallback: AppView = selectedProjectId ? "project-dashboard" : "project-browser";
    const destination = (lastNonUtilityView === "settings" || lastNonUtilityView === "node-activity" || lastNonUtilityView === "storage-manager")
      ? fallback
      : lastNonUtilityView;
    setCurrentView(destination);
  }

  const derivedPostw90Path = derivePostw90PathFromWannier90Path(wannier90PathInput);
  const availablePrograms: Array<{ name: string; type: "qe" | "fermisurfer" | "wannier90" | "postw90" }> = [
    ...availableExecutables.map((name) => ({ name, type: "qe" as const })),
    ...(fermiSurferStatus === "Found"
      ? [{ name: "fermisurfer", type: "fermisurfer" as const }]
      : []),
    ...(wannier90Status === "Found"
      ? [{ name: "wannier90.x", type: "wannier90" as const }]
      : []),
    ...(wannier90Status === "Found"
      ? [{ name: "postw90.x", type: "postw90" as const }]
      : []),
  ];
  const qeStatusClass = qeStatus === "Found" ? "ready" : qeStatus === "Not found" ? "error" : "pending";
  const fermiStatusClass =
    fermiSurferStatus === "Found" ? "ready" : fermiSurferStatus === "Not found" ? "error" : "pending";
  const wannierStatusClass =
    wannier90Status === "Found" ? "ready" : wannier90Status === "Not found" ? "error" : "pending";

  const hpcCpuAdditionalSbatchText = (hpcDefaultCpuDraft.additional_sbatch || []).join("\n");
  const hpcGpuAdditionalSbatchText = (hpcDefaultGpuDraft.additional_sbatch || []).join("\n");

  async function handleConfirmQuit() {
    try {
      setShowCloseConfirm(false);
      await invoke("shutdown_and_close");
    } catch (e) {
      setShowCloseConfirm(true);
      setError(`Quit failed: ${e}`);
    }
  }

  // Render the quit confirmation modal overlay
  const closeConfirmModal = showCloseConfirm ? (
    <div className="close-confirm-overlay" onClick={() => setShowCloseConfirm(false)}>
      <div className="close-confirm-dialog" onClick={(e) => e.stopPropagation()}>
        <h3>Quit QCortado?</h3>
        <p>
          Quit will cancel running jobs and clear pending queue.
        </p>
        <div className="close-confirm-actions">
          <button onClick={() => setShowCloseConfirm(false)}>Cancel</button>
          <button className="close-confirm-danger" onClick={handleConfirmQuit}>
            Quit Anyway
          </button>
        </div>
      </div>
    </div>
  ) : null;

  const settingsWorkspace = showSettingsMenu ? (
    <div className="settings-workspace-shell" ref={settingsMenuRef}>
      <div className="settings-workspace" role="region" aria-label="Settings">
            <AppHeaderPortal className="settings-window-header">
              {settingsPage === "hpc-profile" && (
                <button className="settings-header-back" onClick={leaveHpcProfileEditor} aria-label="Back to HPC settings">
                  Back
                </button>
              )}
              <h3>Settings</h3>
              <div className="settings-header-actions">
                {settingsPage === "hpc-profile" && (
                  <button type="submit" form="hpc-profile-editor-form" className="settings-header-save">
                    Save
                  </button>
                )}
                <button
                  className="settings-window-close"
                  onClick={() => closeSettingsMenu()}
                  aria-label="Close settings"
                >
                  &times;
                </button>
              </div>
            </AppHeaderPortal>
            {settingsPage !== "hpc-profile" && <div className="settings-page-nav">
              <button
                className={`settings-page-tab ${settingsPage === "general" ? "active" : ""}`}
                onClick={() => setSettingsPage("general")}
              >
                General
              </button>
              <button
                className={`settings-page-tab ${settingsPage === "hpc" ? "active" : ""}`}
                onClick={() => setSettingsPage("hpc")}
                disabled={executionMode !== "hpc"}
              >
                HPC
              </button>
            </div>}
            <div className="settings-window-content">
              {settingsPage === "hpc-profile" && editingHpcProfile ? (
                <HpcProfileEditor
                  profile={editingHpcProfile}
                  installations={engineInstallations}
                  onDirtyChange={setHpcProfileEditorDirty}
                  onSaved={(_profile, message) => {
                    setHpcStatus(message);
                    void loadHpcExecutionSettings();
                  }}
                />
              ) : (
              <div className="settings-pane">
              <div className="settings-menu-section">
                <label className="settings-menu-label" htmlFor="execution-mode-select">
                  Execution Mode
                </label>
                <select
                  id="execution-mode-select"
                  className="settings-menu-input"
                  value={executionMode}
                  onChange={(event) => {
                    const nextMode = event.target.value === "hpc" ? "hpc" : "local";
                    void handleExecutionModeChange(nextMode);
                  }}
                >
                  <option value="local">Local</option>
                  <option value="hpc">HPC (Andromeda)</option>
                </select>
                <p className="settings-menu-hint">
                  HPC target host: <code>andromeda.bc.edu</code>. Access usually requires BC network or VPN.
                </p>
                {executionMode === "hpc" && settingsPage === "hpc" && (
                  <div className="settings-hpc-profile-shell">
                    <label className="settings-menu-label" htmlFor="hpc-profile-select">
                      Active HPC Profile
                    </label>
                    <select
                      id="hpc-profile-select"
                      className="settings-menu-input"
                      value={activeHpcProfileId ?? ""}
                      onChange={(event) => {
                        const value = event.target.value;
                        if (value) {
                          void handleSelectHpcProfile(value);
                        }
                      }}
                    >
                      {hpcProfiles.length === 0 ? (
                        <option value="">No HPC profile configured</option>
                      ) : (
                        hpcProfiles.map((profile) => (
                          <option key={profile.id} value={profile.id}>
                            {profile.name} ({profile.username}@{profile.host})
                          </option>
                        ))
                      )}
                    </select>
                    <div className="settings-hpc-actions">
                      <button
                        className="settings-menu-item"
                        onClick={() => {
                          setEditingHpcProfileId(null);
                          setShowHpcSetupWizard(true);
                        }}
                      >
                        New HPC Profile
                      </button>
                      <button
                        className="settings-menu-item"
                        onClick={() => {
                          setEditingHpcProfileId(activeHpcProfileId);
                          setHpcProfileEditorDirty(false);
                          setSettingsPage("hpc-profile");
                        }}
                        disabled={!activeHpcProfile}
                      >
                        Edit Active Profile
                      </button>
                      {activeHpcProfileId && (
                        <button
                          className="settings-menu-item warning"
                          onClick={() => void handleDeleteHpcProfile(activeHpcProfileId)}
                        >
                          Delete Active Profile
                        </button>
                      )}
                      <button
                        className="settings-menu-item"
                        onClick={openMigrateHpcRootsDialog}
                        disabled={!activeHpcProfile || isMigratingHpcRoots}
                      >
                        {isMigratingHpcRoots ? "Migrating Remote..." : "Migrate Remote Roots"}
                      </button>
                      <button
                        className="settings-menu-item warning"
                        onClick={openCleanRemoteConfirmDialog}
                        disabled={!activeHpcProfile || isCleaningHpcRemote}
                      >
                        {isCleaningHpcRemote ? "Cleaning Remote..." : "Clean Remote Orphans"}
                      </button>
                      <div className="settings-repair-button-wrap">
                        <InfoTooltip text="Scans the active HPC profile's remote pseudo directory and fixes only PSLibrary UPFs where the tenth atomic wavefunction tag was written as PP_CHI.1 instead of PP_CHI.10.">
                          <button
                            className="settings-menu-item"
                            onClick={() => void repairRemotePslibraryPseudos()}
                            disabled={!activeHpcProfile || isRepairingRemotePseudos}
                          >
                            {isRepairingRemotePseudos ? "Repairing Remote Pseudos..." : "Repair Remote PSLibrary Pseudos"}
                          </button>
                        </InfoTooltip>
                      </div>
                      <button
                        className="settings-menu-item"
                        onClick={() => void handleExportHpcPresetBundle()}
                        disabled={isExportingHpcPresetBundle || isImportingHpcPresetBundle || hpcProfiles.length === 0}
                      >
                        {isExportingHpcPresetBundle ? "Exporting..." : "Export Presets + Defaults"}
                      </button>
                      <button
                        className="settings-menu-item"
                        onClick={() => void handleImportHpcPresetBundle()}
                        disabled={isImportingHpcPresetBundle || isExportingHpcPresetBundle}
                      >
                        {isImportingHpcPresetBundle ? "Importing..." : "Import Presets + Defaults"}
                      </button>
                    </div>
                  </div>
                )}
                {executionMode === "hpc" && settingsPage === "general" && (
                  <button
                    className="settings-menu-item"
                    onClick={() => setSettingsPage("hpc")}
                  >
                    Open HPC Settings
                  </button>
                )}
                {settingsPage === "hpc" && hpcStatus && <div className="settings-menu-status">{hpcStatus}</div>}
                {settingsPage === "hpc" && pseudoRepairStatus && <div className="settings-menu-status">{pseudoRepairStatus}</div>}
              </div>

              {executionMode === "hpc" && settingsPage === "hpc" && activeHpcProfile && (
                <>
                  <div className="settings-menu-divider" />
                  <div className="settings-menu-section settings-hpc-defaults-section">
                    <div className="settings-menu-label settings-field-label">
                      HPC Default Run Settings
                      <InfoTooltip text="In HPC mode, these become SLURM #SBATCH requests first. After allocation, QCortado launches QE inside that job with the selected launcher." />
                    </div>
                    <p className="settings-menu-hint">
                      These values prefill the HPC run block in all calculation wizards. They remain editable per run.
                    </p>
                    <div className="settings-hpc-launcher-grid">
                      <label>
                        <span className="settings-field-label">
                          Supported Resource Types
                          <InfoTooltip text="Choose whether this profile can submit CPU jobs, GPU jobs, or both. This also controls which QE bin path (CPU or GPU) is used." />
                        </span>
                        <select
                          value={hpcResourceModeDraft}
                          onChange={(event) => {
                            const value = event.target.value;
                            const mode = value === "cpu_only" || value === "gpu_only" ? value : "both";
                            setHpcResourceModeDraft(mode);
                            clearHpcDefaultsSaveFeedback();
                          }}
                        >
                          <option value="both">CPU + GPU</option>
                          <option value="cpu_only">CPU only</option>
                          <option value="gpu_only">GPU only</option>
                        </select>
                      </label>
                      <label>
                        <span className="settings-field-label">
                          MPI Launcher
                          <InfoTooltip text="How QE is started inside the SLURM allocation. With mpirun, QCortado uses the SLURM task count as MPI ranks." />
                        </span>
                        <select
                          value={hpcLauncherDraft}
                          onChange={(event) => {
                            const launcher = event.target.value === "mpirun" ? "mpirun" : "srun";
                            setHpcLauncherDraft(launcher);
                            clearHpcDefaultsSaveFeedback();
                          }}
                        >
                          <option value="srun">srun</option>
                          <option value="mpirun">mpirun</option>
                        </select>
                      </label>
                      <label>
                        <span className="settings-field-label">
                          CPU Launcher Extra Args
                          <InfoTooltip text="Optional CPU-run flags appended to srun/mpirun before the QE executable command." />
                        </span>
                        <input
                          value={hpcLauncherCpuExtraArgsDraft}
                          onChange={(event) => {
                            setHpcLauncherCpuExtraArgsDraft(normalizeCliDashText(event.target.value));
                            clearHpcDefaultsSaveFeedback();
                          }}
                          placeholder="e.g. --bind-to none"
                          autoCorrect="off"
                          autoCapitalize="off"
                          spellCheck={false}
                        />
                      </label>
                      <label>
                        <span className="settings-field-label">
                          GPU Launcher Extra Args
                          <InfoTooltip text="Optional GPU-run flags appended to srun/mpirun before the QE executable command." />
                        </span>
                        <input
                          value={hpcLauncherGpuExtraArgsDraft}
                          onChange={(event) => {
                            setHpcLauncherGpuExtraArgsDraft(normalizeCliDashText(event.target.value));
                            clearHpcDefaultsSaveFeedback();
                          }}
                          placeholder="e.g. --gpu-bind=closest"
                          autoCorrect="off"
                          autoCapitalize="off"
                          spellCheck={false}
                        />
                      </label>
                    </div>
                    <p className="settings-menu-hint">
                      Launcher controls execution inside the allocated Slurm job; it does not change resource requests.
                    </p>

                    <div className="settings-hpc-defaults-grid">
                      <div className="settings-hpc-default-card">
                        <h4>CPU Defaults</h4>
                        <div className="hpc-grid">
                          <label>
                            Partition
                            <input
                              value={hpcDefaultCpuDraft.partition || ""}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  partition: event.target.value,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder="short"
                            />
                          </label>
                          <label>
                            Walltime (HH:MM:SS)
                            <input
                              value={hpcDefaultCpuDraft.walltime || ""}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  walltime: event.target.value,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder="02:00:00"
                            />
                          </label>
                          <label>
                            Nodes
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultCpuDraft.nodes ?? 1}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  nodes: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            <span className="settings-field-label">
                              Tasks
                              <InfoTooltip text="SLURM --ntasks. This is the MPI rank count for the job." />
                            </span>
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultCpuDraft.ntasks ?? 1}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  ntasks: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            <span className="settings-field-label">
                              CPUs / Task
                              <InfoTooltip text="SLURM --cpus-per-task. CPUs per MPI rank. MPI-only runs often use 1; hybrid MPI+OpenMP runs use higher values." />
                            </span>
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultCpuDraft.cpus_per_task ?? 1}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  cpus_per_task: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            Memory (GB)
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultCpuDraft.memory_gb ?? 16}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  memory_gb: normalizePositiveIntInput(event.target.value, 16),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                        </div>
                        <div className="hpc-advanced-grid">
                          <label>
                            QoS
                            <input
                              value={hpcDefaultCpuDraft.qos || ""}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  qos: event.target.value || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            Account
                            <input
                              value={hpcDefaultCpuDraft.account || ""}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  account: event.target.value || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            Constraint
                            <input
                              value={hpcDefaultCpuDraft.constraint || ""}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  constraint: event.target.value || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label className="hpc-advanced-wide">
                            Module / Preamble
                            <textarea
                              rows={3}
                              value={hpcDefaultCpuDraft.module_preamble || ""}
                              onChange={(event) => {
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  module_preamble: normalizeCliDashText(event.target.value) || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder={"module purge\nmodule load qe"}
                              autoCorrect="off"
                              autoCapitalize="off"
                              spellCheck={false}
                            />
                          </label>
                          <label className="hpc-advanced-wide">
                            Extra SBATCH Lines
                            <textarea
                              rows={3}
                              value={hpcCpuAdditionalSbatchText}
                              onChange={(event) => {
                                const parsed = normalizeCliDashText(event.target.value)
                                  .split(/\r?\n/)
                                  .map((line) => line.trim())
                                  .filter((line) => line.length > 0);
                                setHpcDefaultCpuDraft((prev) => ({
                                  ...prev,
                                  additional_sbatch: parsed,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder={"--mail-type=END\n--mail-user=you@example.edu"}
                              autoCorrect="off"
                              autoCapitalize="off"
                              spellCheck={false}
                            />
                          </label>
                        </div>
                      </div>

                      <div className="settings-hpc-default-card">
                        <h4>GPU Defaults</h4>
                        <div className="hpc-grid">
                          <label>
                            Partition
                            <input
                              value={hpcDefaultGpuDraft.partition || ""}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  partition: event.target.value,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder="short"
                            />
                          </label>
                          <label>
                            Walltime (HH:MM:SS)
                            <input
                              value={hpcDefaultGpuDraft.walltime || ""}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  walltime: event.target.value,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder="02:00:00"
                            />
                          </label>
                          <label>
                            Nodes
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultGpuDraft.nodes ?? 1}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  nodes: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            <span className="settings-field-label">
                              Tasks
                              <InfoTooltip text="SLURM --ntasks. This is the MPI rank count for the GPU job." />
                            </span>
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultGpuDraft.ntasks ?? 1}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  ntasks: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            <span className="settings-field-label">
                              CPUs / Task
                              <InfoTooltip text="SLURM --cpus-per-task for each MPI rank in the GPU job. Increase if host-side work is heavy." />
                            </span>
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultGpuDraft.cpus_per_task ?? 1}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  cpus_per_task: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            Memory (GB)
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultGpuDraft.memory_gb ?? 32}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  memory_gb: normalizePositiveIntInput(event.target.value, 32),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            <span className="settings-field-label">
                              GPUs
                              <InfoTooltip text="SLURM --gres=gpu:N. Number of GPUs requested. A common starting point is ntasks approximately equal to gpus." />
                            </span>
                            <input
                              type="number"
                              min={1}
                              value={hpcDefaultGpuDraft.gpus ?? 1}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  gpus: normalizePositiveIntInput(event.target.value, 1),
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                        </div>
                        <div className="hpc-advanced-grid">
                          <label>
                            QoS
                            <input
                              value={hpcDefaultGpuDraft.qos || ""}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  qos: event.target.value || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            Account
                            <input
                              value={hpcDefaultGpuDraft.account || ""}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  account: event.target.value || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label>
                            Constraint
                            <input
                              value={hpcDefaultGpuDraft.constraint || ""}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  constraint: event.target.value || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                            />
                          </label>
                          <label className="hpc-advanced-wide">
                            Module / Preamble
                            <textarea
                              rows={3}
                              value={hpcDefaultGpuDraft.module_preamble || ""}
                              onChange={(event) => {
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  module_preamble: normalizeCliDashText(event.target.value) || null,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder={"module purge\nmodule load qe"}
                              autoCorrect="off"
                              autoCapitalize="off"
                              spellCheck={false}
                            />
                          </label>
                          <label className="hpc-advanced-wide">
                            Extra SBATCH Lines
                            <textarea
                              rows={3}
                              value={hpcGpuAdditionalSbatchText}
                              onChange={(event) => {
                                const parsed = normalizeCliDashText(event.target.value)
                                  .split(/\r?\n/)
                                  .map((line) => line.trim())
                                  .filter((line) => line.length > 0);
                                setHpcDefaultGpuDraft((prev) => ({
                                  ...prev,
                                  additional_sbatch: parsed,
                                }));
                                clearHpcDefaultsSaveFeedback();
                              }}
                              placeholder={"--mail-type=END\n--mail-user=you@example.edu"}
                              autoCorrect="off"
                              autoCapitalize="off"
                              spellCheck={false}
                            />
                          </label>
                        </div>
                      </div>
                    </div>

                    <button
                      className={`settings-menu-item ${hpcDefaultsSaved ? "saved" : ""}`.trim()}
                      onClick={() => void saveHpcDefaultRunSettings()}
                      disabled={isSavingHpcDefaults}
                    >
                      {isSavingHpcDefaults
                        ? "Saving..."
                        : hpcDefaultsSaved
                          ? "✓ Saved"
                          : "Save HPC Default Run Settings"}
                    </button>
                    {hpcDefaultsStatus && <div className="settings-menu-status">{hpcDefaultsStatus}</div>}
                  </div>
                </>
              )}

              {settingsPage === "general" && (
                <>
              {executionMode === "local" && (
                <>
                  <div className="settings-menu-divider" />
                  <div className="settings-menu-section settings-local-tools">
                    <label className="settings-menu-label">Local Executable Paths</label>
                    <p className="settings-menu-hint">
                      These are only used when execution mode is set to Local.
                    </p>

                    <div className="settings-local-config-grid">
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

          <div className="config-row">
            <label>Wannier90:</label>
            <input
              type="text"
                          className="config-path-input"
                          value={wannier90PathInput}
                          onChange={(e) => {
                            setWannier90PathInput(e.target.value);
                            setError(null);
                          }}
                          placeholder={resolveDefaultWannier90Path(qePath)}
                          spellCheck={false}
                        />
                        <div className="config-row-actions">
                          <button onClick={selectWannier90Path}>Browse</button>
                          <button
                            onClick={() => void saveWannier90Path()}
                            disabled={isSavingWannier90Path}
                          >
                            {isSavingWannier90Path ? "Saving..." : "Save"}
              </button>
            </div>
          </div>

          <div className="config-row">
            <label>Default Smearing:</label>
            <select
              className="config-path-input"
              value={qeDefaults.smearing}
              onChange={(e) => {
                setQeDefaults((prev) => ({ ...prev, smearing: normalizeQeSmearing(e.target.value) }));
                setQeDefaultsStatus(null);
              }}
            >
              <option value="marzari-vanderbilt">marzari-vanderbilt</option>
              <option value="gaussian">gaussian</option>
              <option value="methfessel-paxton">methfessel-paxton</option>
              <option value="fermi-dirac">fermi-dirac</option>
            </select>
            <div className="config-row-actions">
              <button
                onClick={() => void saveQeDefaults()}
                disabled={isSavingQeDefaults}
              >
                {isSavingQeDefaults ? "Saving..." : "Save"}
              </button>
            </div>
          </div>
          {qeDefaultsStatus && <p className="settings-menu-hint">{qeDefaultsStatus}</p>}

          <div className="config-row">
            <label>postw90 (auto):</label>
                        <input
                          type="text"
                          className="config-path-input"
                          value={derivedPostw90Path}
                          readOnly
                          spellCheck={false}
                        />
                        <div className="config-row-actions">
                          <InfoTooltip text="Derived from the configured Wannier90 path.">
                            <button disabled aria-label="Derived from the configured Wannier90 path.">Auto</button>
                          </InfoTooltip>
                        </div>
                      </div>
                    </div>

                    <p className="settings-menu-hint">
                      QCortado derives <code>postw90.x</code> from the Wannier90 executable path by using the same directory.
                    </p>

                    <div className="settings-repair-button-wrap">
                      <InfoTooltip text="Scans the pseudo directory next to the configured local QE bin path and fixes only PSLibrary UPFs where the tenth atomic wavefunction tag was written as PP_CHI.1 instead of PP_CHI.10.">
                        <button
                          className="settings-menu-item"
                          onClick={() => void repairLocalPslibraryPseudos()}
                          disabled={isRepairingLocalPseudos || getLocalPseudoDir().length === 0}
                        >
                          {isRepairingLocalPseudos ? "Repairing Local Pseudos..." : "Repair Local PSLibrary Pseudos"}
                        </button>
                      </InfoTooltip>
                    </div>
                    {pseudoRepairStatus && <div className="settings-menu-status">{pseudoRepairStatus}</div>}

                    <div className="settings-local-status-grid">
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

                      <div className="status-row">
                        <label>Wannier90 Status:</label>
                        <span className={`status ${wannierStatusClass}`}>
                          {wannier90Status}
                        </span>
                      </div>
                    </div>

                    {error && <div className="error">{error}</div>}

                    {availablePrograms.length > 0 && (
                      <div className="executables">
                        <label>Available programs:</label>
                        <div className="exe-list">
                          {availablePrograms.map((program) => (
                            <span
                              key={program.name}
                              className={`exe-tag ${program.type === "fermisurfer" || program.type === "wannier90" || program.type === "postw90" ? "exe-tag-fermisurfer" : ""}`}
                            >
                              {program.name}
                            </span>
                          ))}
                        </div>
                      </div>
                    )}
                  </div>
                </>
              )}

              <div className="settings-menu-divider" />
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
                  Non-SCF Save Size
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
                  <option value="large">Large (keep full non-SCF restart files)</option>
                  <option value="small">Small (compact non-SCF saves)</option>
                </select>
                <p className="settings-menu-hint">
                  This only affects non-SCF saves. SCF always keeps restart files required for phonon workflows.
                </p>
                <button
                  className="settings-menu-item"
                  onClick={() => void saveGlobalSaveSizeSetting()}
                  disabled={isSavingSaveSizeMode}
                >
                  {isSavingSaveSizeMode ? "Saving..." : "Save Non-SCF Size Mode"}
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
                  disabled={isRecoveringPhonon || isRecoveringRemotePhonon || !selectedProjectId}
                >
                  {isRecoveringPhonon ? "Recovering..." : "Recover Phonon"}
                </button>
                <button
                  className="settings-menu-item"
                  onClick={() => void recoverRemotePhononFromSettings()}
                  disabled={isRecoveringPhonon || isRecoveringRemotePhonon || showRemotePhononSelectionDialog || !selectedProjectId}
                >
                  {isRecoveringRemotePhonon ? "Recovering..." : "Recover Remote Phonon"}
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
                </>
              )}
              </div>
              )}
            </div>
      </div>
    </div>
  ) : null;

  const remotePhononSelectionModal = showRemotePhononSelectionDialog ? (
    <div
      className="remote-recovery-overlay"
      onPointerDown={(e) => e.stopPropagation()}
      onMouseDown={(e) => e.stopPropagation()}
      onClick={() => closeRemotePhononSelectionDialog()}
    >
      <div
        className="remote-recovery-dialog"
        onPointerDown={(e) => e.stopPropagation()}
        onMouseDown={(e) => e.stopPropagation()}
        onClick={(e) => e.stopPropagation()}
      >
        <div className="dialog-header">
          <h2>Recover Remote Phonon</h2>
          <button
            className="dialog-close"
            onClick={() => closeRemotePhononSelectionDialog()}
            disabled={isRecoveringRemotePhonon}
          >
            &times;
          </button>
        </div>

        <div className="remote-recovery-body">
          <p className="exit-hint">
            Select the remote phonon run to import. Entries are sorted by newest timestamp first.
          </p>
          {isLoadingRemotePhononCandidates ? (
            <div className="remote-recovery-state">Loading recoverable remote phonon runs...</div>
          ) : remotePhononCandidates.length === 0 ? (
            <div className="remote-recovery-state">
              No recoverable remote phonon runs found. Check the recovery status message for debug details.
            </div>
          ) : (
            <div className="remote-recovery-list" role="radiogroup" aria-label="Recoverable remote phonon runs">
              {remotePhononCandidates.map((candidate) => {
                const isSelected = selectedRemotePhononWorkdir === candidate.remote_workdir;
                return (
                  <label
                    key={`${candidate.remote_workdir}-${candidate.modified_at_epoch}`}
                    className={`radio-option remote-recovery-option ${isSelected ? "selected" : ""}`.trim()}
                  >
                    <input
                      type="radio"
                      name="remote-phonon-recovery-choice"
                      checked={isSelected}
                      onChange={() => setSelectedRemotePhononWorkdir(candidate.remote_workdir)}
                      disabled={isRecoveringRemotePhonon}
                    />
                    <div className="remote-recovery-option-details">
                      <strong>{formatEpochSeconds(candidate.modified_at_epoch)}</strong>
                      <div className="remote-recovery-option-meta">
                        <span>Source: {formatRemotePhononCandidateLocation(candidate.location)}</span>
                      </div>
                      <code className="remote-recovery-option-path">{candidate.remote_workdir}</code>
                    </div>
                  </label>
                );
              })}
            </div>
          )}
        </div>

        <div className="remote-recovery-footer">
          <button
            className="remote-recovery-btn remote-recovery-btn-cancel"
            onClick={() => closeRemotePhononSelectionDialog()}
            disabled={isRecoveringRemotePhonon}
          >
            Cancel
          </button>
          <button
            className="remote-recovery-btn remote-recovery-btn-primary"
            onClick={() => void confirmRemotePhononRecoverySelection()}
            disabled={
              isLoadingRemotePhononCandidates
              || isRecoveringRemotePhonon
              || !selectedRemotePhononWorkdir
            }
          >
            {isRecoveringRemotePhonon ? "Recovering..." : "Recover Selected Run"}
          </button>
        </div>
      </div>
    </div>
  ) : null;

  const headlessRecoveryModal = showHeadlessRecoveryDialog ? (
    <div
      className="remote-recovery-overlay"
      onClick={() => !isAttachingHeadlessJob && setShowHeadlessRecoveryDialog(false)}
    >
      <div className="remote-recovery-dialog" onClick={(e) => e.stopPropagation()}>
        <div className="dialog-header">
          <h2>Recover Headless Job</h2>
          <button
            className="dialog-close"
            onClick={() => setShowHeadlessRecoveryDialog(false)}
            disabled={isAttachingHeadlessJob}
          >
            &times;
          </button>
        </div>

        <div className="remote-recovery-body">
          {isLoadingHeadlessCandidates ? (
            <div className="remote-recovery-state">Querying cluster jobs...</div>
          ) : headlessCandidates.length === 0 ? (
            <div className="remote-recovery-state">
              {headlessRecoveryStatus || "No recoverable QCortado jobs found."}
            </div>
          ) : (
            <div className="remote-recovery-list" role="radiogroup" aria-label="Recoverable headless HPC jobs">
              {headlessCandidates.map((candidate) => {
                const isSelected = selectedHeadlessJobId === candidate.remote_job_id;
                return (
                  <label
                    key={`${candidate.remote_job_id}-${candidate.remote_workdir || ""}`}
                    className={`radio-option remote-recovery-option ${isSelected ? "selected" : ""}`.trim()}
                  >
                    <input
                      type="radio"
                      name="headless-hpc-recovery-choice"
                      checked={isSelected}
                      onChange={() => setSelectedHeadlessJobId(candidate.remote_job_id)}
                      disabled={isAttachingHeadlessJob}
                    />
                    <div className="remote-recovery-option-details">
                      <strong>
                        Job {candidate.remote_job_id} · {candidate.task_kind.toUpperCase()}
                      </strong>
                      <div className="remote-recovery-option-meta">
                        <span>{candidate.scheduler_state}</span>
                        <span>{candidate.remote_node || "No node yet"}</span>
                        <span>{candidate.auto_save_available ? "Auto-save ready" : "Limited legacy recovery"}</span>
                      </div>
                      <code className="remote-recovery-option-path">{candidate.remote_workdir || "Remote workdir unavailable"}</code>
                    </div>
                  </label>
                );
              })}
            </div>
          )}
          {headlessRecoveryStatus && headlessCandidates.length > 0 && (
            <div className="remote-recovery-state">{headlessRecoveryStatus}</div>
          )}
        </div>

        <div className="remote-recovery-footer">
          <button
            className="remote-recovery-btn remote-recovery-btn-cancel"
            onClick={() => setShowHeadlessRecoveryDialog(false)}
            disabled={isAttachingHeadlessJob}
          >
            Cancel
          </button>
          <button
            className="remote-recovery-btn remote-recovery-btn-primary"
            onClick={() => void handleAttachHeadlessJob()}
            disabled={isLoadingHeadlessCandidates || isAttachingHeadlessJob || !selectedHeadlessJobId}
          >
            {isAttachingHeadlessJob ? "Attaching..." : "Attach Job"}
          </button>
        </div>
      </div>
    </div>
  ) : null;

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

  const cleanRemoteConfirmModal = showCleanRemoteConfirmDialog ? (
    <div className="dialog-overlay" onClick={() => !isCleaningHpcRemote && setShowCleanRemoteConfirmDialog(false)}>
      <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
        <div className="dialog-header">
          <h2>Clean Remote Orphans</h2>
          <button
            className="dialog-close"
            onClick={() => setShowCleanRemoteConfirmDialog(false)}
            disabled={isCleaningHpcRemote}
          >
            &times;
          </button>
        </div>

        <div className="dialog-body">
          <p className="exit-warning">
            Delete remote QCortado bundle directories not referenced by local projects?
          </p>
          <p className="exit-hint">
            Referenced calculations are preserved. Only orphaned remote directories are removed.
          </p>
        </div>

        <div className="dialog-footer">
          <button
            className="dialog-btn cancel"
            onClick={() => setShowCleanRemoteConfirmDialog(false)}
            disabled={isCleaningHpcRemote}
          >
            Cancel
          </button>
          <button
            className="dialog-btn delete"
            onClick={() => void handleConfirmCleanRemoteOrphans()}
            disabled={isCleaningHpcRemote}
          >
            {isCleaningHpcRemote ? "Cleaning..." : "Clean Orphans"}
          </button>
        </div>
      </div>
    </div>
  ) : null;

  const migrateHpcRootsModal = showMigrateHpcRootsDialog ? (
    <div className="dialog-overlay" onClick={() => !isMigratingHpcRoots && setShowMigrateHpcRootsDialog(false)}>
      <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
        <div className="dialog-header">
          <h2>Migrate Remote Roots</h2>
          <button
            className="dialog-close"
            onClick={() => setShowMigrateHpcRootsDialog(false)}
            disabled={isMigratingHpcRoots}
          >
            &times;
          </button>
        </div>

        <div className="dialog-body">
          <p className="exit-warning">
            This will copy all data to new roots and then remove old roots.
          </p>
          <p className="exit-hint">
            QCortado will copy all contents from the current remote workspace/project roots to the new roots, then remove the old roots and update this profile.
          </p>
          <div className="form-group">
            <label>
              New Remote Workspace Root
            </label>
            <input
              type="text"
              value={migrateWorkspaceRootDraft}
              onChange={(e) => setMigrateWorkspaceRootDraft(e.target.value)}
              placeholder="e.g. ~/scratch/qcortado"
              disabled={isMigratingHpcRoots}
              autoFocus
              autoCorrect="off"
              autoCapitalize="off"
              spellCheck={false}
            />
          </div>
          <div className="form-group">
            <label>
              New Remote Project Root
            </label>
            <input
              type="text"
              value={migrateProjectRootDraft}
              onChange={(e) => setMigrateProjectRootDraft(e.target.value)}
              placeholder="e.g. ~/projects/qcortado"
              disabled={isMigratingHpcRoots}
              autoCorrect="off"
              autoCapitalize="off"
              spellCheck={false}
            />
          </div>
        </div>

        <div className="dialog-footer">
          <button
            className="dialog-btn cancel"
            onClick={() => setShowMigrateHpcRootsDialog(false)}
            disabled={isMigratingHpcRoots}
          >
            Cancel
          </button>
          <button
            className="dialog-btn save"
            onClick={() => void handleConfirmMigrateHpcRoots()}
            disabled={isMigratingHpcRoots || migrateWorkspaceRootDraft.trim().length === 0 || migrateProjectRootDraft.trim().length === 0}
          >
            {isMigratingHpcRoots ? "Migrating..." : "Migrate Roots"}
          </button>
        </div>
      </div>
    </div>
  ) : null;

  const activeDialogModal =
    closeConfirmModal
    || remotePhononSelectionModal
    || headlessRecoveryModal
    || deleteProjectModal
    || cleanRemoteConfirmModal
    || migrateHpcRootsModal;

  const appContextLabel: Record<AppView, string> = {
    "project-browser": "Projects",
    "project-dashboard": "Project",
    "bands-multiview": "Band Comparison",
    "settings": "Settings",
    "node-activity": "Node Activity",
    "storage-manager": "Storage Manager",
    "scf-wizard": "SCF Workflow",
    "bands-wizard": "Band Structure Workflow",
    "bands-viewer": "Band Structure",
    "dos-wizard": "Electronic DOS Workflow",
    "dos-viewer": "Electronic DOS",
    "wannier-wizard": "Wannier90 Workflow",
    "wannier-viewer": "Wannier90",
    "transport-wizard": "Transport Workflow",
    "transport-viewer": "Transport",
    "fermi-surface-wizard": "Fermi Surface Workflow",
    "hubbard-lrt-wizard": "Hubbard LRT Workflow",
    "phonon-wizard": "Phonon Workflow",
    "phonon-viewer": "Phonons",
    "epw-wizard": "EPW Workflow",
    "epw-viewer": "EPW",
    "wien2k-structure-wizard": "WIEN2k Structure Workflow",
    "wien2k-scf-wizard": "WIEN2k SCF Workflow",
  };

  const appChrome = (
    <>
      {!isHpcActivityPopout && (
        <AppTopBar
          contextLabel={appContextLabel[currentView]}
          task={relevantTaskManagerEntry}
          summary={taskManagerSummary}
          onOpenMenu={openNavigationDrawer}
          onOpenTasks={(taskId) => openTaskDrawer("active", taskId ?? null)}
          embeddedHeader
        />
      )}
      <AppNavigationDrawer
        isOpen={showNavigationDrawer}
        executionMode={executionMode}
        activeHpcProfile={activeHpcProfile}
        taskSummary={taskManagerSummary}
        onClose={() => setShowNavigationDrawer(false)}
        onProjects={navigateToProjects}
        onStorage={navigateToStorageManager}
        onNodeActivity={navigateToNodeActivity}
        onRecoverHeadless={() => {
          setShowNavigationDrawer(false);
          void openHeadlessRecoveryDialog();
        }}
        onSettings={openSettingsWorkspace}
        onHpcTasks={() => {
          setShowNavigationDrawer(false);
          openTaskDrawer("hpc");
        }}
      />
      <TaskManagerDrawer
        isOpen={showTaskDrawer}
        requestedFilter={taskDrawerFilter}
        focusedTaskId={taskDrawerFocusId}
        onClose={() => setShowTaskDrawer(false)}
        onNavigateToTask={(taskId, taskType) => {
          setShowTaskDrawer(false);
          handleNavigateToTask(taskId, taskType);
        }}
      />
      {activeDialogModal}
      <HpcSetupWizard
        isOpen={showHpcSetupWizard}
        initialProfile={null}
        onClose={() => {
          setShowHpcSetupWizard(false);
          setEditingHpcProfileId(null);
        }}
        onSaved={() => {
          setEditingHpcProfileId(null);
          setHpcStatus("HPC profile saved.");
          void loadHpcExecutionSettings();
        }}
      />
    </>
  );

  if (isHpcActivityPopout) {
    return (
      <div className="hpc-activity-standalone-shell">
        <HpcActivityPanel standalone />
      </div>
    );
  }

  if (currentView === "settings") {
    return (
      <>
        {settingsWorkspace}
        {appChrome}
      </>
    );
  }

  if (currentView === "node-activity") {
    return (
      <>
        <HpcNodeActivityPage
          onBack={returnFromUtilityView}
          executionMode={executionMode}
          activeProfileId={activeHpcProfileId}
          activeProfileName={activeHpcProfile?.name ?? "Andromeda"}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "storage-manager") {
    return (
      <>
        <StorageManagerPage
          onBack={returnFromUtilityView}
          executionMode={executionMode}
          activeHpcProfile={activeHpcProfile}
        />
        {appChrome}
      </>
    );
  }

  const workflowHostRuntime = {
    qePath,
    defaultSmearing: qeDefaults.smearing,
    executionMode,
    onExecutionModeChange: handleExecutionModeChange,
    activeHpcProfile,
  };
  const workflowHostContexts = {
    scf: scfContext,
    bands: bandsContext,
    dos: dosContext,
    wannier: wannierContext,
    fermiSurface: fermiSurfaceContext,
    hubbardLrt: hubbardLrtContext,
    phonons: phononsContext,
    transport: transportContext,
    epw: epwContext,
    structureSetup: wien2kStructureContext,
    wien2kScf: scfContext,
  };

  if (
    activeWorkflowRoute &&
    currentView === activeWorkflowRoute.view &&
    canRenderEngineWorkflowHost({
      route: activeWorkflowRoute,
      runtime: workflowHostRuntime,
      contexts: workflowHostContexts,
      reconnectTaskId,
    })
  ) {
    return (
      <>
        <EngineWorkflowHost
          route={activeWorkflowRoute}
          runtime={workflowHostRuntime}
          contexts={workflowHostContexts}
          reconnectTaskId={reconnectTaskId}
          onBack={handleWorkflowBack}
          onViewBands={(bandData, fermiEnergy, calculationParameters, calculationContext) => {
            setViewBandsData({
              bandData,
              fermiEnergy,
              calculationParameters,
              calculationContext,
            });
            setActiveWorkflowRoute(null);
            setCurrentView("bands-viewer");
            setReconnectTaskId(null);
          }}
          onViewDos={(dosData, fermiEnergy) => {
            setViewDosData({ dosData, fermiEnergy });
            setActiveWorkflowRoute(null);
            setCurrentView("dos-viewer");
            setReconnectTaskId(null);
          }}
          onViewWannier={(result, fermiEnergy, overlayOptions = []) => {
            setViewWannierData({ result, fermiEnergy, overlayOptions });
            setActiveWorkflowRoute(null);
            setCurrentView("wannier-viewer");
            setReconnectTaskId(null);
          }}
          onViewTransport={(transportData) => {
            setViewTransportData({ data: transportData });
            setActiveWorkflowRoute(null);
            setCurrentView("transport-viewer");
            setReconnectTaskId(null);
          }}
          onViewPhonons={(phononData, viewMode) => {
            setViewPhononData({
              data: phononData,
              mode: viewMode,
            });
            setActiveWorkflowRoute(null);
            setCurrentView("phonon-viewer");
            setReconnectTaskId(null);
          }}
          onViewEpw={(epwData, rawOutput) => {
            setViewEpwData({
              data: epwData,
              rawOutput: rawOutput ?? null,
            });
            setActiveWorkflowRoute(null);
            setCurrentView("epw-viewer");
            setReconnectTaskId(null);
          }}
          onStructureSourceSaved={() => {
            setProjectDashboardRefreshToken((previous) => previous + 1);
            handleWorkflowBack("wien2k-structure-wizard", "project-dashboard");
          }}
          onWien2kScfSaved={() => {
            setProjectDashboardRefreshToken((previous) => previous + 1);
            handleWorkflowBack("wien2k-scf-wizard", "project-dashboard");
          }}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "bands-viewer" && viewBandsData) {
    return (
      <>
        <div className="bands-viewer-container">
          <AppHeaderPortal className="bands-viewer-header">
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
          </AppHeaderPortal>
          <div className="bands-viewer-content">
            <BandPlot
              data={viewBandsData.bandData}
              scfFermiEnergy={viewBandsData.fermiEnergy ?? undefined}
              calculationParameters={viewBandsData.calculationParameters ?? null}
              onPersistSelectedValenceBandIndex={
                viewBandsData.calculationContext
                  ? async (bandIndex) => {
                    const updatedCalculation = await invoke<CalculationRun>(
                      "update_calculation_band_viewer_metadata",
                      {
                        projectId: viewBandsData.calculationContext?.projectId,
                        cifId: viewBandsData.calculationContext?.cifId,
                        calcId: viewBandsData.calculationContext?.calcId,
                        selectedValenceBandIndex: bandIndex,
                      },
                    );
                    setViewBandsData((prev) => prev
                      ? {
                        ...prev,
                        calculationParameters: (updatedCalculation.parameters ?? null) as Record<string, unknown> | null,
                      }
                      : prev);
                  }
                  : undefined
              }
              viewerType="electronic"
            />
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "dos-viewer" && viewDosData) {
    return (
      <>
        <div className="bands-viewer-container">
          <AppHeaderPortal className="bands-viewer-header">
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
          </AppHeaderPortal>
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

  if (currentView === "wannier-viewer" && viewWannierData) {
    const result = viewWannierData.result;
    const wannierIssues = getWannierQualityIssues(result, null, viewWannierData.fermiEnergy ?? null);
    return (
      <>
        <div className="bands-viewer-container">
          <AppHeaderPortal className="bands-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewWannierData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>Wannier90</h2>
          </AppHeaderPortal>
          <div className="bands-viewer-content bands-viewer-content-stacked">
            <div className="bands-viewer-plot-region">
              <BandPlot
                data={result.band_data}
                scfFermiEnergy={viewWannierData.fermiEnergy ?? undefined}
                viewerType="electronic"
                comparisonOptions={viewWannierData.overlayOptions}
                comparisonTitle="Saved Band Overlay"
                comparisonNoneLabel="No overlay"
              />
            </div>
            <div className="bands-viewer-details-region">
              {wannierIssues.length > 0 && (
                <div className="warning-banner">
                  {wannierIssues.map((issue) => issue.message).join(" ")}
                </div>
              )}
              <div className="details-grid">
              <div className="detail-item">
                <label>seedname</label>
                <span>{result.seedname || "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>num_wann</label>
                <span>{result.num_wann ?? "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>num_bands</label>
                <span>{result.num_bands ?? "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>Total Spread</label>
                <span>{result.total_spread != null ? `${Number(result.total_spread).toFixed(6)} A^2` : "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>Converged</label>
                <span>{result.convergence?.converged ? "Yes" : "No"}</span>
              </div>
              <div className="detail-item">
                <label>Iterations</label>
                <span>{result.convergence?.iterations ?? "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>Minimization</label>
                <span>{formatWannierConvergenceFlag(result.convergence?.minimization_converged)}</span>
              </div>
              <div className="detail-item">
                <label>Disentanglement</label>
                <span>{formatWannierConvergenceFlag(result.convergence?.disentanglement_converged)}</span>
              </div>
            </div>
            </div>
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "transport-viewer" && viewTransportData) {
    return (
      <>
        <div className="bands-viewer-container transport-viewer-container">
          <AppHeaderPortal className="bands-viewer-header transport-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewTransportData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>BoltzWann Transport</h2>
          </AppHeaderPortal>
          <div className="bands-viewer-content transport-viewer-content">
            <TransportPlot data={viewTransportData.data} />
          </div>
        </div>
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
          <AppHeaderPortal className="phonon-viewer-header">
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
          </AppHeaderPortal>
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

  if (currentView === "epw-viewer" && viewEpwData) {
    return (
      <>
        <EpwViewer
          payload={viewEpwData}
          onBack={() => {
            setCurrentView("project-dashboard");
            setViewEpwData(null);
          }}
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
          onSelectProject={(projectId, folderId) => {
            setSelectedProjectId(projectId);
            setProjectBrowserFolderId(folderId);
            setCurrentView("project-dashboard");
          }}
          onOpenBandsMultiview={(initialCalculations) => {
            setBandsMultiviewInitialCalculations(initialCalculations);
            setCurrentView("bands-multiview");
          }}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "bands-multiview") {
    return (
      <>
        <BandsMultiview
          onBack={() => {
            setBandsMultiviewInitialCalculations(undefined);
            setCurrentView("project-browser");
          }}
          initialCalculations={bandsMultiviewInitialCalculations}
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
          refreshToken={projectDashboardRefreshToken}
          onBack={() => {
            setCurrentView("project-browser");
            setSelectedProjectId(null);
          }}
          onDeleted={() => {
            setCurrentView("project-browser");
            setSelectedProjectId(null);
          }}
          onRunSCF={(engineId, cifId, crystalData, cifContent, filename, preset, presetLock, optimizedStructures, calculations) => {
            setScfContext({
              cifId,
              crystalData,
              cifContent,
              filename,
              projectId: selectedProjectId,
              initialPreset: preset,
              presetLock,
              optimizedStructures,
              calculations,
            });
            openEngineWorkflow(engineId, preset === "relax" ? "structure_optimization" : "scf", "scf-wizard");
          }}
          onContinueWien2kScf={(cifId, crystalData, cifContent, filename, calculations, calculationId) => {
            setScfContext({
              cifId,
              crystalData,
              cifContent,
              filename,
              projectId: selectedProjectId,
              calculations,
              continuationCalculationId: calculationId,
            });
            openEngineWorkflow("wien2k", "scf", "wien2k-scf-wizard");
          }}
          onRunEngineSetup={(engineId, cifId, crystalData) => {
            setWien2kStructureContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
            });
            openEngineWorkflow(engineId, "engine_setup", "wien2k-structure-wizard");
          }}
          onRunBands={(engineId, cifId, crystalData, scfCalculations) => {
            setBandsContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            openEngineWorkflow(engineId, "bands", "bands-wizard");
          }}
          onViewBands={(bandData, fermiEnergy, calculationParameters, calculationContext) => {
            setViewBandsData({
              bandData,
              fermiEnergy,
              calculationParameters,
              calculationContext,
            });
            setCurrentView("bands-viewer");
          }}
          onRunDos={(engineId, cifId, crystalData, scfCalculations) => {
            setDosContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            openEngineWorkflow(engineId, "dos", "dos-wizard");
          }}
          onViewDos={(dosData, fermiEnergy) => {
            setViewDosData({ dosData, fermiEnergy });
            setCurrentView("dos-viewer");
          }}
          onRunWannier={(engineId, cifId, crystalData, scfCalculations) => {
            setWannierContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            openEngineWorkflow(engineId, "wannier", "wannier-wizard");
          }}
          onViewWannier={(wannierData, fermiEnergy, overlayOptions = []) => {
            setViewWannierData({ result: wannierData, fermiEnergy, overlayOptions });
            setCurrentView("wannier-viewer");
          }}
          onRunTransport={(engineId, cifId, crystalData, wannierCalculations) => {
            setTransportContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              wannierCalculations,
            });
            openEngineWorkflow(engineId, "transport", "transport-wizard");
          }}
          onViewTransport={(transportData) => {
            setViewTransportData({ data: transportData });
            setCurrentView("transport-viewer");
          }}
          onRunFermiSurface={(engineId, cifId, crystalData, scfCalculations) => {
            setFermiSurfaceContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            openEngineWorkflow(engineId, "fermi_surface", "fermi-surface-wizard");
          }}
          onRunPhonons={(engineId, cifId, crystalData, scfCalculations) => {
            setPhononsContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            openEngineWorkflow(engineId, "phonon", "phonon-wizard");
          }}
          onRunHubbardLrt={(engineId, cifId, crystalData, scfCalculations) => {
            setHubbardLrtContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              scfCalculations,
            });
            openEngineWorkflow(engineId, "hubbard_lrt", "hubbard-lrt-wizard");
          }}
          onRunEPW={(engineId, cifId, crystalData, calculations) => {
            setEpwContext({
              cifId,
              crystalData,
              projectId: selectedProjectId,
              calculations,
            });
            openEngineWorkflow(engineId, "epw", "epw-wizard");
          }}
          onViewPhonons={(phononData, viewMode) => {
            setViewPhononData({
              data: phononData,
              mode: viewMode,
            });
            setCurrentView("phonon-viewer");
          }}
          onViewEPW={(epwData, rawOutput) => {
            setViewEpwData({
              data: epwData,
              rawOutput: rawOutput ?? null,
            });
            setCurrentView("epw-viewer");
          }}
        />
        {appChrome}
      </>
    );
  }

  return (
    <>
      <ProjectBrowser
        initialActiveFolderId={projectBrowserFolderId}
        onSelectProject={(projectId, folderId) => {
          setSelectedProjectId(projectId);
          setProjectBrowserFolderId(folderId);
          setCurrentView("project-dashboard");
        }}
        onOpenBandsMultiview={(initialCalculations) => {
          setBandsMultiviewInitialCalculations(initialCalculations);
          setCurrentView("bands-multiview");
        }}
      />
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
