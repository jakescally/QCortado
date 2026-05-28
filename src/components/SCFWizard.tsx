// SCF Calculation Wizard - Import CIF, configure, and run SCF calculations

import { useState, useEffect, useCallback, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import {
  CrystalData,
  ELEMENT_MASSES,
  ExecutionMode,
  HpcProfile,
  OptimizedStructureOption,
  PseudopotentialMetadata,
  QePositionUnit,
  SavedCellSummary,
  SavedStructureData,
  SCFPreset,
  SlurmResourceRequest,
} from "../lib/types";
import { parseCIF } from "../lib/cifParser";
import { UnitCellViewer } from "./UnitCellViewer";
import { SaveToProjectDialog } from "./SaveToProjectDialog";
import {
  analyzeCrystalSymmetry,
  buildConventionalLatticeFromCrystalData,
  SymmetryTransformResult,
} from "../lib/symmetryTransform";
import { inferQeBravaisCellFromCif } from "../lib/qeBravaisInference";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { InfoTooltip } from "./InfoTooltip";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { isPhononReadyScf } from "../lib/phononReady";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import {
  buildExecutionTarget,
  buildHpcQeInputCommandLine,
  defaultResourcesForProfile,
  getRemotePseudopotentialMetadata,
  loadRemoteSsspData,
  listRemotePseudopotentialInventory,
  listRemotePseudopotentials,
  listRemotePseudopotentialMetadata,
  resolveProfileRemoteQeBinDir,
  resolveProfileRemotePseudoDir,
} from "../lib/hpcConfig";
import { HpcRunSettings } from "./HpcRunSettings";
import { RemoteUtilizationPanel } from "./RemoteUtilizationPanel";
import {
  getDefaultHubbardManifold,
  getOutermostOccupiedOrbitalManifold,
} from "../lib/electronConfigurations";
import {
  getHubbardRecommendations,
  getLatestHubbardLrtValue,
  getHundJDefaultEv,
  GENERAL_HUBBARD_U_GUESS_EV,
  HUBBARD_J_SOURCE,
} from "../lib/hubbard";
import {
  CutoffDerivation,
  CutoffProvenance,
  CutoffStatus,
  SSSPElementData,
  summarizeSelectedPseudoCutoffs,
} from "../lib/pseudopotentialCutoffs";
import {
  getStoredScfRunSettingsClipboardText,
  hasStoredScfRunSettingsClipboardText,
  parseScfRunSettingsClipboardText,
  rememberScfRunSettingsClipboardText,
  SCF_RUN_SETTINGS_UPDATED_EVENT,
  ScfRunSettingsClipboardPayload,
} from "../lib/scfRunSettingsClipboard";
import { readProjectWizardSettings, writeProjectWizardSettings } from "../lib/projectWizardSettings";
import {
  buildPseudopotentialFilenameSignature,
  buildPseudopotentialInventorySignature,
  isPseudopotentialMetadataCacheFresh,
  readPseudopotentialMetadataCache,
  updateCachedPseudopotentialMetadata,
  writePseudopotentialMetadataCache,
} from "../lib/pseudopotentialMetadataCache";
import type { PseudopotentialInventoryEntry } from "../lib/pseudopotentialMetadataCache";
interface InitialCifData {
  cifId: string;
  crystalData: CrystalData;
  cifContent: string;
  filename: string;
  projectId: string;
  calculations?: any[];
}

interface SCFWizardProps {
  onBack: () => void;
  qePath: string;
  defaultSmearing?: "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac";
  executionMode?: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
  /** Pre-loaded CIF data from project dashboard */
  initialCif?: InitialCifData;
  initialPreset?: SCFPreset;
  presetLock?: boolean;
  optimizedStructures?: OptimizedStructureOption[];
  /** If provided, reconnect to a running/completed background task */
  reconnectTaskId?: string;
}

interface SCFConfig {
  // Calculation type
  calculation: "scf" | "relax" | "vcrelax";

  // Basic parameters
  ecutwfc: number;
  ecutrho: number;
  kgrid: [number, number, number];
  kgrid_offset: [number, number, number];

  // Relaxation parameters (for relax/vcrelax)
  forc_conv_thr: number;
  etot_conv_thr: number;
  press: number;  // Target pressure for vcrelax (kbar)

  // Electronic structure
  occupations: "smearing" | "tetrahedra" | "tetrahedra_lin" | "tetrahedra_opt" | "fixed" | "from_input";
  smearing: "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" | "cold";
  degauss: number;
  nbnd: number | null;  // null = auto
  tot_charge: number;

  // Magnetism & Spin
  nspin: 1 | 2 | 4;
  noncolin: boolean;
  lspinorb: boolean;
  starting_magnetization: Record<string, number>;  // per element
  starting_magnetization_theta: Record<string, number>;  // polar angle from z, degrees
  starting_magnetization_phi: Record<string, number>;  // azimuthal angle in xy, degrees
  tot_magnetization: number | null;  // null = auto, only for nspin=2
  constrained_magnetization: "none" | "total" | "atomic" | "total direction" | "atomic direction";

  // SCF Convergence
  conv_thr: number;
  electron_maxstep: number;
  mixing_mode: "plain" | "TF" | "local-TF";
  mixing_beta: number;
  mixing_ndim: number;
  diagonalization: "david" | "cg" | "ppcg" | "paro" | "rmm-davidson";
  startingpot: "atomic" | "file";
  startingwfc: "atomic" | "atomic+random" | "random" | "file";

  // DFT+U
  lda_plus_u: boolean;
  lda_plus_u_kind: 0 | 1 | 2;
  hubbard_projector: HubbardProjector;
  hubbard_manifold: Record<string, string>;  // per element, e.g. 3d
  hubbard_u: Record<string, number>;  // per element
  hubbard_j: Record<string, number>;  // per element (for kind=1,2)

  // Van der Waals
  vdw_corr: "none" | "grimme-d2" | "grimme-d3" | "ts-vdw" | "xdm" | "dft-d";

  // Isolated systems
  assume_isolated: "none" | "makov-payne" | "martyna-tuckerman" | "esm" | "2D";

  // XC functional override
  input_dft: string;  // empty = use pseudopotential default

  // Output control
  verbosity: "low" | "high";
  tprnfor: boolean;
  tstress: boolean;
  disk_io: "low" | "medium" | "high" | "nowf";
}

type WizardStep = "import" | "configure" | "run" | "results";
type CellViewMode = "conventional" | "primitive";
type HubbardProjector = "atomic" | "ortho-atomic" | "norm-atomic" | "wf" | "pseudo";
type HubbardParameter = {
  parameter: "U" | "J" | "J0";
  manifold: string;
  value: number;
};

interface DisplayCellMetrics {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

interface ScfTaskPlan {
  calculation: any;
  inputText: string;
  taskLabel: string;
  taskParams: Record<string, any>;
  sourceStructure: SavedStructureData;
  magnetismViewerStructure: SavedStructureData;
  sourceDescriptor: { type: "cif" | "optimization"; calc_id?: string };
}

function normalizeDefaultSmearing(
  raw: unknown,
): "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" {
  const lowered = String(raw || "").toLowerCase();
  if (lowered === "gaussian") return "gaussian";
  if (lowered === "methfessel-paxton") return "methfessel-paxton";
  if (lowered === "fermi-dirac") return "fermi-dirac";
  return "marzari-vanderbilt";
}

function normalizeHubbardManifold(raw: string): string {
  return raw.trim().replace(/\s+/g, "");
}

function asRecord(value: unknown): Record<string, unknown> | null {
  if (!value || typeof value !== "object" || Array.isArray(value)) {
    return null;
  }
  return value as Record<string, unknown>;
}

function asFiniteNumber(value: unknown): number | null {
  const parsed = typeof value === "number" ? value : Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function asBoolean(value: unknown): boolean | null {
  if (typeof value === "boolean") return value;
  if (value === "true") return true;
  if (value === "false") return false;
  return null;
}

function asNumberTriplet(value: unknown): [number, number, number] | null {
  if (!Array.isArray(value) || value.length !== 3) return null;
  const parsed = value.map((entry) => asFiniteNumber(entry));
  if (parsed.some((entry) => entry === null)) return null;
  return [parsed[0]!, parsed[1]!, parsed[2]!];
}

function toIntegerTriplet(value: [number, number, number]): [number, number, number] {
  return value.map((entry) => Math.max(1, Math.round(entry))) as [number, number, number];
}

function toOffsetTriplet(value: [number, number, number]): [number, number, number] {
  return value.map((entry) => (entry >= 0.5 ? 1 : 0)) as [number, number, number];
}

type PseudopotentialPreset = "sssp" | "paw" | "uspp" | "ncpp";
const SCF_WIZARD_SETTINGS_ID = "scf";
const OPTIMIZATION_WIZARD_SETTINGS_ID = "optimization";
const PSEUDO_CACHE_REMOTE_CHECK_TTL_MS = 5 * 60 * 1000;

interface StoredScfWizardSettings {
  config: SCFConfig;
  selectedPseudos: Record<string, string>;
  selectedPseudoPreset: PseudopotentialPreset;
  structureSource: string;
}

export function SCFWizard({
  onBack,
  qePath,
  defaultSmearing = "marzari-vanderbilt",
  executionMode = "local",
  activeHpcProfile = null,
  initialCif,
  initialPreset,
  presetLock,
  optimizedStructures = [],
  reconnectTaskId,
}: SCFWizardProps) {
  const resolvedDefaultSmearing = normalizeDefaultSmearing(defaultSmearing);
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const lockedPreset = presetLock && initialPreset ? initialPreset : null;
  const isOptimizationWizard = lockedPreset === "relax";
  const wizardSettingsId = isOptimizationWizard ? OPTIMIZATION_WIZARD_SETTINGS_ID : SCF_WIZARD_SETTINGS_ID;
  const storedWizardSettings = useMemo(
    () => readProjectWizardSettings<StoredScfWizardSettings>(initialCif?.projectId, wizardSettingsId),
    [initialCif?.projectId, wizardSettingsId],
  );
  // Track background task for this wizard
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  const hasBlockingExternalTask = !isHpcMode && hasExternalRunningTask;
  // If we have initial CIF data, skip to configure step
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : (initialCif ? "configure" : "import"));
  const [crystalData, setCrystalData] = useState<CrystalData | null>(initialCif?.crystalData || null);
  const [cifFilename, setCifFilename] = useState<string>(initialCif?.filename || "");
  const [cifContent, setCifContent] = useState<string>(initialCif?.cifContent || "");
  const [error, setError] = useState<string | null>(null);
  const [copiedRunSettingsAvailable, setCopiedRunSettingsAvailable] = useState<boolean>(
    hasStoredScfRunSettingsClipboardText(),
  );
  const [copiedRunSettingsMessage, setCopiedRunSettingsMessage] = useState<string | null>(null);

  // Track project context for saving
  const [projectContext, setProjectContext] = useState<{ projectId: string; cifId: string } | null>(
    initialCif ? { projectId: initialCif.projectId, cifId: initialCif.cifId } : null
  );
  const [pseudopotentials, setPseudopotentials] = useState<PseudopotentialMetadata[]>([]);
  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>(
    () => storedWizardSettings?.selectedPseudos ?? {},
  );
  const selectedPseudoMetadata = useMemo(() => {
    const next: Record<string, PseudopotentialMetadata> = {};
    for (const [element, filename] of Object.entries(selectedPseudos)) {
      const metadata = pseudopotentials.find((pseudo) => pseudo.filename === filename);
      if (metadata) {
        next[element] = metadata;
      }
    }
    return next;
  }, [pseudopotentials, selectedPseudos]);

  const WORK_DIR = "/tmp/qcortado_work";

  // Save dialog state
  const [showSaveDialog, setShowSaveDialog] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [calcEndTime, setCalcEndTime] = useState<string>("");
  const [generatedInput, setGeneratedInput] = useState<string>("");
  const [resultSaved, setResultSaved] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [runSourceStructure, setRunSourceStructure] = useState<SavedStructureData | null>(null);
  const [runMagnetismViewerStructure, setRunMagnetismViewerStructure] = useState<SavedStructureData | null>(null);
  const [runSourceDescriptor, setRunSourceDescriptor] = useState<{ type: "cif" | "optimization"; calc_id?: string } | null>(null);
  const [cellViewMode, setCellViewMode] = useState<CellViewMode>("conventional");
  const [symmetryTransform, setSymmetryTransform] = useState<SymmetryTransformResult | null>(null);
  const [symmetryError, setSymmetryError] = useState<string | null>(null);
  const activeQueueItem = useMemo(
    () => (activeTaskId ? taskContext.queueItems.find((item) => item.taskId === activeTaskId) ?? null : null),
    [activeTaskId, taskContext.queueItems],
  );
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const hasTaskLinkedAutosave = Boolean(activeQueueItem?.saveSpec);
  const autoSaveExpected = Boolean(projectContext || hasTaskLinkedAutosave);
  const initialConfig = useMemo<SCFConfig>(() => {
    const storedConfig = (storedWizardSettings?.config ?? {}) as Partial<SCFConfig>;
    const storedCalculation = storedConfig.calculation;
    const normalizedCalculation = isOptimizationWizard
      ? (storedCalculation === "scf" || storedCalculation === "relax" || storedCalculation === "vcrelax"
        ? storedCalculation
        : "vcrelax")
      : "scf";
    return {
      // Basic parameters
      ecutwfc: 40,
      ecutrho: 320,
      kgrid: [4, 4, 4],
      kgrid_offset: [0, 0, 0],

      // Relaxation parameters
      forc_conv_thr: 1e-4,
      etot_conv_thr: 1e-5,
      press: 0.0,

      // Electronic structure
      occupations: "smearing",
      smearing: resolvedDefaultSmearing,
      degauss: 0.01,
      nbnd: null,
      tot_charge: 0,

      // Magnetism & Spin
      nspin: 1,
      noncolin: false,
      lspinorb: false,
      starting_magnetization: {},
      starting_magnetization_theta: {},
      starting_magnetization_phi: {},
      tot_magnetization: null,
      constrained_magnetization: "none",

      // SCF Convergence
      conv_thr: 1e-12,
      electron_maxstep: 1000,
      mixing_mode: "plain",
      mixing_beta: 0.7,
      mixing_ndim: 8,
      diagonalization: "david",
      startingpot: "atomic",
      startingwfc: "atomic",

      // DFT+U
      lda_plus_u: false,
      lda_plus_u_kind: 0,
      hubbard_projector: "ortho-atomic",
      hubbard_manifold: {},
      hubbard_u: {},
      hubbard_j: {},

      // Van der Waals
      vdw_corr: "none",

      // Isolated systems
      assume_isolated: "none",

      // XC functional override
      input_dft: "",

      // Output control
      verbosity: "high",
      tprnfor: true,
      tstress: true,
      disk_io: "low",
      ...storedConfig,
      calculation: normalizedCalculation,
    };
  }, [isOptimizationWizard, resolvedDefaultSmearing, storedWizardSettings?.config]);

  const [config, setConfig] = useState<SCFConfig>(() => ({
    ...initialConfig,
  }));
  const [manuallyEditedHubbardU, setManuallyEditedHubbardU] = useState<Record<string, boolean>>({});
  const [manuallyEditedHubbardJ, setManuallyEditedHubbardJ] = useState<Record<string, boolean>>({});
  const [hubbardJDefaultLabels, setHubbardJDefaultLabels] = useState<Record<string, string>>({});
  const [selectedPreset, setSelectedPreset] = useState<SCFPreset | null>(null);
  const [structureSource, setStructureSource] = useState<string>(() => storedWizardSettings?.structureSource ?? "cif");
  const selectedOptimizedStructure =
    structureSource === "cif"
      ? null
      : optimizedStructures.find((option) => option.calcId === structureSource) || null;
  const phononPresetDisabled = structureSource === "cif";

  const applyPreset = useCallback((preset: SCFPreset) => {
    if (preset === "phonon" && phononPresetDisabled) {
      return;
    }
    setSelectedPreset(preset);
    setConfig((prev) => {
      switch (preset) {
        case "standard":
          return {
            ...prev,
            calculation: "scf",
            conv_thr: 1e-12,
            nspin: prev.nspin === 4 ? 1 : prev.nspin,
            noncolin: prev.nspin === 4 ? false : prev.noncolin,
            lspinorb: false,
            disk_io: "low",
          };
        case "phonon":
          return {
            ...prev,
            calculation: "scf",
            conv_thr: 1e-12,
            nspin: prev.nspin === 4 ? 1 : prev.nspin,
            noncolin: prev.nspin === 4 ? false : prev.noncolin,
            lspinorb: false,
            disk_io: "medium",
          };
        case "relax":
          return {
            ...prev,
            calculation: "vcrelax",
            conv_thr: 1e-12,
            nspin: prev.nspin === 4 ? 1 : prev.nspin,
            noncolin: prev.nspin === 4 ? false : prev.noncolin,
            lspinorb: false,
            forc_conv_thr: 1e-4,
            etot_conv_thr: 1e-5,
            press: 0,
            tprnfor: true,
            tstress: true,
          };
        case "soc":
          return {
            ...prev,
            calculation: "scf",
            conv_thr: 1e-12,
            noncolin: true,
            nspin: 4,
            lspinorb: true,
          };
        default:
          return prev;
      }
    });
    if (preset === "phonon" || preset === "relax" || preset === "standard" || preset === "soc") {
      setConvThrInput("1e-12");
    }
  }, [phononPresetDisabled]);

  useEffect(() => {
    if (initialPreset) {
      applyPreset(initialPreset);
    }
  }, [initialPreset, applyPreset]);

  useEffect(() => {
    if (projectContext || !activeQueueItem?.saveSpec) return;
    setProjectContext({
      projectId: activeQueueItem.saveSpec.projectId,
      cifId: activeQueueItem.saveSpec.cifId,
    });
  }, [activeQueueItem, projectContext]);

  useEffect(() => {
    const refresh = () => {
      void refreshCopiedRunSettingsAvailability();
    };
    refresh();
    window.addEventListener(SCF_RUN_SETTINGS_UPDATED_EVENT, refresh);
    window.addEventListener("focus", refresh);
    return () => {
      window.removeEventListener(SCF_RUN_SETTINGS_UPDATED_EVENT, refresh);
      window.removeEventListener("focus", refresh);
    };
  }, []);

  useEffect(() => {
    if (step !== "configure") return;
    void refreshCopiedRunSettingsAvailability();
  }, [step]);

  // Reconnect to a running/completed background task
  useEffect(() => {
    if (!activeTaskId) return;
    if (!activeTask) {
      // Task may not be loaded yet, try reconnecting
      taskContext.reconnectToTask(activeTaskId);
      return;
    }

    // Sync wizard state from task
    setIsRunning(activeTask.status === "running");
    if (activeTask.outputLineCount > 0 || activeTask.status !== "running") {
      setOutput(activeTask.outputText);
      setOutputLineCount(activeTask.outputLineCount);
    }
    if (activeTask.status === "completed" && activeQueueItem?.saveSpec) {
      if (activeQueueItem.status === "completed") {
        setProgress({
          ...activeTask.progress,
          status: "complete",
          percent: 100,
          phase: "Complete",
        });
      } else if (activeQueueItem.status === "failed" || activeQueueItem.status === "cancelled") {
        setProgress({
          ...activeTask.progress,
          status: "error",
          percent: null,
          phase: "Save Error",
        });
      } else {
        setProgress({
          ...activeTask.progress,
          status: "running",
          percent: null,
          phase: "Saving to project",
        });
      }
    } else {
      setProgress(activeTask.progress);
    }
    setCalcStartTime(activeTask.startedAt);

    if (activeTask.status === "completed" && activeTask.result) {
      setResult(activeTask.result);
      setCalcEndTime(new Date().toISOString());
      setStep("results");
    } else if (activeTask.status === "failed" || activeTask.status === "cancelled") {
      setError(activeTask.error || "Task failed");
    } else {
      setStep("run");
    }
  }, [activeQueueItem, activeTask, activeTaskId, taskContext.reconnectToTask]);

  // When active task updates, sync output/progress
  useEffect(() => {
    if (!activeTask || activeTask.status !== "running") return;
    if (activeTask.outputLineCount > 0) {
      setOutput(activeTask.outputText);
      setOutputLineCount(activeTask.outputLineCount);
    }
    setProgress(activeTask.progress);
  }, [activeTask]);

  useEffect(() => {
    if (!autoSaveExpected || !activeTask?.result) return;
    if (!activeQueueItem?.saveSpec || !activeTaskId) return;

    if (activeQueueItem.status === "completed") {
      setIsSaving(false);
      setResultSaved(true);
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));
      return;
    }

    if (activeQueueItem.status === "failed" || activeQueueItem.status === "cancelled") {
      setIsSaving(false);
      setProgress((prev) => ({
        ...prev,
        status: "error",
        percent: null,
        phase: "Save Error",
      }));
      setError(activeQueueItem.error || "Failed to auto-save calculation.");
      return;
    }

    if (activeTask?.status === "completed") {
      setIsSaving(true);
      setProgress((prev) => ({
        ...prev,
        status: "running",
        percent: null,
        phase: "Saving to project",
      }));
    }
  }, [activeQueueItem, activeTask, activeTaskId, autoSaveExpected]);

  const wizardTitle = isOptimizationWizard ? "Structure Optimization Wizard" : "SCF Calculation Wizard";
  const structureSourceTooltip = isOptimizationWizard
    ? "For complicated runs, you may want to run a relaxation with a coarse k-grid, then use that result as a structure source and run a finer relaxation."
    : "Choose the geometry used for this SCF. 'From CIF' uses the original imported structure; optimized structures come from saved geometry optimization runs.";
  const showPresetRow = !lockedPreset || lockedPreset === "relax" || lockedPreset === "soc";
  const showRelaxPreset = lockedPreset === "relax";
  const showSocPreset = !lockedPreset || lockedPreset === "soc";
  const socPresetSelected = selectedPreset === "soc";
  const conventionalCellMetrics = useMemo<DisplayCellMetrics | null>(() => {
    if (!crystalData) return null;
    return {
      a: crystalData.cell_length_a.value,
      b: crystalData.cell_length_b.value,
      c: crystalData.cell_length_c.value,
      alpha: crystalData.cell_angle_alpha.value,
      beta: crystalData.cell_angle_beta.value,
      gamma: crystalData.cell_angle_gamma.value,
    };
  }, [crystalData]);

  const primitiveCellMetrics = useMemo<DisplayCellMetrics | null>(() => {
    if (!symmetryTransform) return null;
    const [v1, v2, v3] = symmetryTransform.standardizedPrimitiveLattice;
    return calculateMetricsFromVectors(v1, v2, v3);
  }, [symmetryTransform]);

  const hasPrimitiveDisplay = primitiveCellMetrics !== null;
  const displayedCellMetrics = cellViewMode === "primitive" && primitiveCellMetrics
    ? primitiveCellMetrics
    : conventionalCellMetrics;

  useEffect(() => {
    if (structureSource === "cif") return;
    if (!selectedOptimizedStructure) {
      setStructureSource("cif");
    }
  }, [structureSource, selectedOptimizedStructure]);

  useEffect(() => {
    if (phononPresetDisabled && selectedPreset === "phonon") {
      setSelectedPreset(null);
    }
  }, [phononPresetDisabled, selectedPreset]);

  useEffect(() => {
    if (!hasPrimitiveDisplay && cellViewMode === "primitive") {
      setCellViewMode("conventional");
    }
  }, [hasPrimitiveDisplay, cellViewMode]);

  // Collapsed section states
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({
    basic: true,
    electronic: false,
    magnetism: false,
    convergence: false,
    dftu: false,
    vdw: false,
    advanced: false,
  });

  const toggleSection = (section: string) => {
    setExpandedSections(prev => ({ ...prev, [section]: !prev[section] }));
  };

  // Local string state for conv_thr input (to allow typing scientific notation)
  const [convThrInput, setConvThrInput] = useState("1e-12");

  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState<string>("");
  const [outputLineCount, setOutputLineCount] = useState(0);
  const [result, setResult] = useState<any>(null);
  const [pseudoDir, setPseudoDir] = useState<string>("");
  const [pseudoFilenames, setPseudoFilenames] = useState<string[]>([]);
  const [pseudoError, setPseudoError] = useState<string | null>(null);
  const [pseudoPresetWarnings, setPseudoPresetWarnings] = useState<Record<string, string>>({});
  const [selectedPseudoPreset, setSelectedPseudoPreset] = useState<PseudopotentialPreset>(
    () => storedWizardSettings?.selectedPseudoPreset ?? "sssp",
  );
  const [isPseudoListLoading, setIsPseudoListLoading] = useState(false);
  const [isPseudoMetadataLoading, setIsPseudoMetadataLoading] = useState(false);
  const [isSsspLoading, setIsSsspLoading] = useState(false);
  const [pseudoListLoaded, setPseudoListLoaded] = useState(false);
  const [pseudoMetadataLoaded, setPseudoMetadataLoaded] = useState(false);
  const [autoCutoffFromPseudosEnabled, setAutoCutoffFromPseudosEnabled] = useState(
    () => !storedWizardSettings?.config,
  );
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "SCF iterations",
  });

  const [ssspData, setSsspData] = useState<Record<string, SSSPElementData> | null>(null);
  const [ssspMissing, setSsspMissing] = useState(false);
  const isPseudoLoading = isPseudoListLoading || isPseudoMetadataLoading || isSsspLoading;
  const pseudoMetadataByFilename = useMemo(
    () => Object.fromEntries(pseudopotentials.map((pseudo) => [pseudo.filename, pseudo])),
    [pseudopotentials],
  );
  const hasPseudoCutoffMetadata = useMemo(
    () => pseudopotentials.some((pseudo) => pseudo.cutoff_wfc != null || pseudo.cutoff_rho != null)
      || Boolean(ssspData && Object.keys(ssspData).length > 0),
    [pseudopotentials, ssspData],
  );

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(
    defaultResourcesForProfile(activeHpcProfile),
  );
  const resolvedPseudoDir = useMemo(
    () => (
      isHpcMode
        ? resolveProfileRemotePseudoDir(activeHpcProfile, hpcResources.resource_type)
        : (qePath ? qePath.replace(/\/bin\/?$/, "/pseudo") : "")
    ),
    [
      activeHpcProfile?.remote_cpu_pseudo_dir,
      activeHpcProfile?.remote_gpu_pseudo_dir,
      activeHpcProfile?.remote_pseudo_dir,
      hpcResources.resource_type,
      isHpcMode,
      qePath,
    ],
  );
  const pseudoCacheScope = useMemo(
    () => ({
      kind: isHpcMode ? "remote" as const : "local" as const,
      pseudoDir: resolvedPseudoDir,
      profileId: isHpcMode ? activeHpcProfile?.id ?? null : null,
    }),
    [activeHpcProfile?.id, isHpcMode, resolvedPseudoDir],
  );
  const visibleOutputLineCount = useMemo(() => countVisibleOutputLines(output), [output]);
  useViewportScrollLock(step === "run");

  useEffect(() => {
    const persistedConfig = isOptimizationWizard
      ? config
      : { ...config, calculation: "scf" as const };
    writeProjectWizardSettings(projectContext?.projectId, wizardSettingsId, {
      config: persistedConfig,
      selectedPseudos,
      selectedPseudoPreset,
      structureSource,
    });
  }, [config, isOptimizationWizard, projectContext?.projectId, selectedPseudos, selectedPseudoPreset, structureSource, wizardSettingsId]);

  useEffect(() => {
    if (isOptimizationWizard || config.calculation === "scf") return;
    setConfig((prev) => (prev.calculation === "scf" ? prev : { ...prev, calculation: "scf" }));
  }, [config.calculation, isOptimizationWizard]);

  useEffect(() => {
    if (visibleOutputLineCount > outputLineCount) {
      setOutputLineCount(visibleOutputLineCount);
    }
  }, [outputLineCount, visibleOutputLineCount]);

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [isHpcMode, activeHpcProfile?.id, activeHpcProfile?.resource_mode]);

  const hpcCommandLines = useMemo(
    () => [
      "cd \"$SLURM_SUBMIT_DIR\"",
      `QE_BIN="${resolveProfileRemoteQeBinDir(activeHpcProfile, hpcResources.resource_type)}"`,
      buildHpcQeInputCommandLine(activeHpcProfile, "pw.x", "pw.in", "pw.out", undefined, hpcResources.resource_type),
    ],
    [activeHpcProfile, hpcResources.resource_type],
  );

  useEffect(() => {
    setPseudoDir(resolvedPseudoDir);
    const cached = readPseudopotentialMetadataCache(pseudoCacheScope);
    setPseudoFilenames(cached?.filenames ?? []);
    setPseudopotentials(cached?.metadata ?? []);
    setSsspData(cached?.ssspData ?? null);
    setPseudoListLoaded(Boolean(cached));
    setPseudoMetadataLoaded(Boolean(cached?.metadata.length));
    setSsspMissing(cached?.ssspMissing ?? false);
    setPseudoPresetWarnings({});
    setIsPseudoListLoading(false);
    setIsPseudoMetadataLoading(false);
    setIsSsspLoading(false);

    if (!resolvedPseudoDir.trim()) {
      setPseudoError(
        isHpcMode
          ? "Remote pseudopotential directory is not configured in the active HPC profile."
          : "Local pseudopotential directory is not configured.",
      );
      return;
    }

    setPseudoError(null);
  }, [isHpcMode, pseudoCacheScope, resolvedPseudoDir]);

  const ensurePseudoFilenamesLoaded = useCallback(async (): Promise<string[]> => {
    if (!resolvedPseudoDir.trim()) {
      const message = isHpcMode
        ? "Remote pseudopotential directory is not configured in the active HPC profile."
        : "Local pseudopotential directory is not configured.";
      setPseudoError(message);
      throw new Error(message);
    }
    if (pseudoListLoaded) {
      return pseudoFilenames;
    }

    setIsPseudoListLoading(true);
    try {
      const filenames = isHpcMode
        ? await listRemotePseudopotentials(resolvedPseudoDir, activeHpcProfile?.id ?? null)
        : await invoke<string[]>("list_pseudopotentials", { pseudoDir: resolvedPseudoDir });
      setPseudoFilenames(filenames);
      setPseudoListLoaded(true);
      setPseudoError(null);
      return filenames;
    } catch (e) {
      console.error("Failed to list pseudopotentials:", e);
      const message = String(e);
      setPseudoError(message);
      setPseudoFilenames([]);
      throw new Error(message);
    } finally {
      setIsPseudoListLoading(false);
    }
  }, [activeHpcProfile?.id, isHpcMode, pseudoFilenames, pseudoListLoaded, resolvedPseudoDir]);

  const ensurePseudoMetadataLoaded = useCallback(async (): Promise<PseudopotentialMetadata[]> => {
    if (!resolvedPseudoDir.trim()) {
      const message = isHpcMode
        ? "Remote pseudopotential directory is not configured in the active HPC profile."
        : "Local pseudopotential directory is not configured.";
      setPseudoError(message);
      throw new Error(message);
    }
    if (pseudoMetadataLoaded) {
      return pseudopotentials;
    }

    setIsPseudoMetadataLoading(true);
    try {
      const pseudos = isHpcMode
        ? await listRemotePseudopotentialMetadata(resolvedPseudoDir, activeHpcProfile?.id ?? null)
        : await invoke<PseudopotentialMetadata[]>("list_pseudopotential_metadata", { pseudoDir: resolvedPseudoDir });
      setPseudopotentials(pseudos);
      setPseudoMetadataLoaded(true);
      setPseudoListLoaded(true);
      setPseudoFilenames(pseudos.map((pseudo) => pseudo.filename));
      writePseudopotentialMetadataCache(pseudoCacheScope, {
        filenames: pseudos.map((pseudo) => pseudo.filename),
        signature: buildPseudopotentialFilenameSignature(pseudos.map((pseudo) => pseudo.filename)),
        metadata: pseudos,
        ssspData,
        ssspMissing,
        updatedAt: new Date().toISOString(),
        checkedAt: new Date().toISOString(),
      });
      setPseudoError(null);
      return pseudos;
    } catch (e) {
      console.error("Failed to load pseudopotential metadata:", e);
      const message = String(e);
      setPseudoError(message);
      setPseudopotentials([]);
      throw new Error(message);
    } finally {
      setIsPseudoMetadataLoading(false);
    }
  }, [activeHpcProfile?.id, isHpcMode, pseudoCacheScope, pseudoMetadataLoaded, pseudopotentials, resolvedPseudoDir, ssspData, ssspMissing]);

  const ensureSsspDataLoaded = useCallback(async (): Promise<Record<string, SSSPElementData> | null> => {
    if (!resolvedPseudoDir.trim()) {
      const message = isHpcMode
        ? "Remote pseudopotential directory is not configured in the active HPC profile."
        : "Local pseudopotential directory is not configured.";
      setPseudoError(message);
      throw new Error(message);
    }
    if (ssspData) {
      return ssspData;
    }
    if (ssspMissing) {
      return null;
    }

    setIsSsspLoading(true);
    try {
      const sssp = isHpcMode
        ? await loadRemoteSsspData<Record<string, SSSPElementData>>(
          resolvedPseudoDir,
          activeHpcProfile?.id ?? null,
        )
        : await invoke<Record<string, SSSPElementData>>("load_sssp_data", {
          pseudoDir: resolvedPseudoDir,
        });
      setSsspData(sssp);
      setSsspMissing(false);
      const cached = readPseudopotentialMetadataCache(pseudoCacheScope);
      if (cached) {
        writePseudopotentialMetadataCache(pseudoCacheScope, {
          ...cached,
          ssspData: sssp,
          ssspMissing: false,
          updatedAt: new Date().toISOString(),
          checkedAt: new Date().toISOString(),
        });
      }
      return sssp;
    } catch (e) {
      console.error("Failed to load SSSP data:", e);
      setSsspData(null);
      setSsspMissing(true);
      const cached = readPseudopotentialMetadataCache(pseudoCacheScope);
      if (cached) {
        writePseudopotentialMetadataCache(pseudoCacheScope, {
          ...cached,
          ssspData: null,
          ssspMissing: true,
          updatedAt: new Date().toISOString(),
          checkedAt: new Date().toISOString(),
        });
      }
      return null;
    } finally {
      setIsSsspLoading(false);
    }
  }, [activeHpcProfile?.id, isHpcMode, pseudoCacheScope, resolvedPseudoDir, ssspData, ssspMissing]);

  useEffect(() => {
    if (!resolvedPseudoDir.trim()) return;
    let cancelled = false;

    async function refreshPseudoRecordIfChanged() {
      const cachedBeforeRefresh = readPseudopotentialMetadataCache(pseudoCacheScope);
      if (isPseudopotentialMetadataCacheFresh(cachedBeforeRefresh, PSEUDO_CACHE_REMOTE_CHECK_TTL_MS)) {
        return;
      }

      const hasUsableCachedRecord = Boolean(cachedBeforeRefresh?.filenames.length);
      if (!hasUsableCachedRecord) {
        setIsPseudoListLoading(true);
      }
      try {
        const inventory = isHpcMode
          ? await listRemotePseudopotentialInventory(resolvedPseudoDir, activeHpcProfile?.id ?? null)
          : await invoke<PseudopotentialInventoryEntry[]>("list_pseudopotential_inventory", {
            pseudoDir: resolvedPseudoDir,
          });
        if (cancelled) return;

        const filenames = inventory.map((entry) => entry.filename);
        const signature = buildPseudopotentialInventorySignature(inventory);
        const cached = readPseudopotentialMetadataCache(pseudoCacheScope);
        setPseudoFilenames(filenames);
        setPseudoListLoaded(true);
        setPseudoError(null);

        if (cached && cached.signature === signature) {
          const checkedAt = new Date().toISOString();
          setPseudopotentials(cached.metadata);
          setPseudoMetadataLoaded(cached.metadata.length > 0);
          setSsspData(cached.ssspData);
          setSsspMissing(cached.ssspMissing);
          writePseudopotentialMetadataCache(pseudoCacheScope, {
            ...cached,
            checkedAt,
          });
          return;
        }

        const metadata = isHpcMode
          ? await listRemotePseudopotentialMetadata(resolvedPseudoDir, activeHpcProfile?.id ?? null)
          : await invoke<PseudopotentialMetadata[]>("list_pseudopotential_metadata", { pseudoDir: resolvedPseudoDir });
        if (cancelled) return;

        let nextSsspData: Record<string, SSSPElementData> | null = null;
        let nextSsspMissing = false;
        try {
          nextSsspData = isHpcMode
            ? await loadRemoteSsspData<Record<string, SSSPElementData>>(resolvedPseudoDir, activeHpcProfile?.id ?? null)
            : await invoke<Record<string, SSSPElementData>>("load_sssp_data", { pseudoDir: resolvedPseudoDir });
        } catch {
          nextSsspMissing = true;
        }
        if (cancelled) return;

        setPseudopotentials(metadata);
        setPseudoMetadataLoaded(true);
        setSsspData(nextSsspData);
        setSsspMissing(nextSsspMissing);
        writePseudopotentialMetadataCache(pseudoCacheScope, {
          filenames,
          signature,
          metadata,
          ssspData: nextSsspData,
          ssspMissing: nextSsspMissing,
          updatedAt: new Date().toISOString(),
          checkedAt: new Date().toISOString(),
        });
      } catch (e) {
        if (!cancelled) {
          console.error("Failed to refresh pseudopotential metadata cache:", e);
          setPseudoError(String(e));
        }
      } finally {
        if (!cancelled) {
          setIsPseudoListLoading(false);
          setIsPseudoMetadataLoading(false);
        }
      }
    }

    void refreshPseudoRecordIfChanged();
    return () => {
      cancelled = true;
    };
  }, [activeHpcProfile?.id, isHpcMode, pseudoCacheScope, resolvedPseudoDir]);

  const ensurePseudoMetadataForFilename = useCallback(async (filename: string): Promise<void> => {
    if (!filename || pseudoMetadataByFilename[filename]) return;
    if (!resolvedPseudoDir.trim()) return;

    setIsPseudoMetadataLoading(true);
    try {
      const metadata = isHpcMode
        ? await getRemotePseudopotentialMetadata(filename, resolvedPseudoDir, activeHpcProfile?.id ?? null)
        : await invoke<PseudopotentialMetadata>("get_pseudopotential_metadata", {
          pseudoDir: resolvedPseudoDir,
          filename,
        });
      setPseudopotentials((prev) => {
        const next = prev.filter((pseudo) => pseudo.filename !== filename);
        next.push(metadata);
        next.sort((left, right) => left.filename.localeCompare(right.filename));
        return next;
      });
      setPseudoMetadataLoaded(true);
      updateCachedPseudopotentialMetadata(pseudoCacheScope, filename, metadata);
      setPseudoError(null);
    } catch (e) {
      console.error("Failed to parse selected pseudopotential metadata:", e);
      setPseudoError(String(e));
    } finally {
      setIsPseudoMetadataLoading(false);
    }
  }, [activeHpcProfile?.id, isHpcMode, pseudoCacheScope, pseudoMetadataByFilename, resolvedPseudoDir]);

  // Load CPU count and check MPI availability
  useEffect(() => {
    async function loadMpiInfo() {
      try {
        const cores = await invoke<number>("get_cpu_count");
        const safeCores = Math.max(1, Math.floor(cores));
        setCpuCount(safeCores);
        const defaults = await loadGlobalMpiDefaults(safeCores);
        const mpiOk = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(mpiOk);
        setMpiEnabled(mpiOk ? defaults.enabled : false);
        setMpiProcs(defaults.nprocs);
      } catch (e) {
        console.error("Failed to load MPI info:", e);
      }
    }
    loadMpiInfo();
  }, []);

  // Strip oxidation state from element symbol (e.g., "Ni0+" -> "Ni", "Fe3+" -> "Fe")
  function getBaseElement(symbol: string): string {
    return symbol.replace(/[\d+-]+$/, "");
  }

  // Helper to check if a pseudopotential file matches an element
  // Matches patterns like: Si.pbe-..., Si_ONCV_..., Si-pbe-..., etc.
  function pseudoMatchesElement(filename: string, element: string): boolean {
    const lowerFile = filename.toLowerCase();
    const lowerEl = getBaseElement(element).toLowerCase();
    // Match element followed by dot, underscore, or hyphen
    return (
      lowerFile.startsWith(lowerEl + ".") ||
      lowerFile.startsWith(lowerEl + "_") ||
      lowerFile.startsWith(lowerEl + "-")
    );
  }

  function isPslibraryPseudoFilename(filename: string): boolean {
    const lowerFile = filename.toLowerCase();
    return lowerFile.includes("pslibrary") || /(^|[._-])psl([._-]|$)/.test(lowerFile);
  }

  function matchesPseudoPresetType(
    pseudo: PseudopotentialMetadata | undefined,
    preset: Exclude<PseudopotentialPreset, "sssp">,
  ): boolean {
    const pseudoType = (pseudo?.pseudo_type || "").toUpperCase();
    switch (preset) {
      case "paw":
        return pseudoType.includes("PAW");
      case "uspp":
        return pseudoType.includes("US");
      case "ncpp":
        return pseudoType.includes("NC");
      default:
        return false;
    }
  }

  function getPseudoPresetDisplayName(preset: PseudopotentialPreset): string {
    switch (preset) {
      case "paw":
        return "PAW";
      case "uspp":
        return "USPP";
      case "ncpp":
        return "NCPP";
      default:
        return "SSSP";
    }
  }

  function choosePreferredPseudoFilename(
    element: string,
    candidates: string[],
    availableSsspData?: Record<string, SSSPElementData> | null,
  ): string | null {
    if (candidates.length === 0) {
      return null;
    }
    const ssspFilename = availableSsspData?.[element]?.filename?.toLowerCase();
    if (ssspFilename) {
      const ssspMatch = candidates.find((filename) => filename.toLowerCase() === ssspFilename);
      if (ssspMatch) {
        return ssspMatch;
      }
    }
    const pbeMatch = candidates.find((filename) => filename.toLowerCase().includes("pbe"));
    return pbeMatch || candidates[0];
  }

  function findFilenameCaseInsensitive(candidates: string[], target: string | null | undefined): string | null {
    if (!target) {
      return null;
    }
    return candidates.find((candidate) => candidate.toLowerCase() === target.toLowerCase()) || null;
  }

  function clearPseudoPresetWarning(element: string) {
    setPseudoPresetWarnings((prev) => {
      if (!prev[element]) {
        return prev;
      }
      const next = { ...prev };
      delete next[element];
      return next;
    });
  }

  function cutoffProvenanceLabel(provenance: CutoffProvenance): string {
    switch (provenance) {
      case "sssp":
        return "SSSP";
      case "djrepo":
        return "PseudoDojo";
      case "upf_info":
        return "UPF info";
      case "upf":
        return "UPF";
      case "upf_fallback":
        return "UPF fallback";
      case "unknown":
        return "metadata";
      default:
        return "";
    }
  }

  function cutoffStatusLabel(
    status: CutoffStatus,
    provenance: CutoffProvenance,
    derivation: CutoffDerivation,
    ratio: number | null,
  ): string {
    switch (status) {
      case "parsed":
        return provenance === "sssp"
          ? "SSSP"
          : `Parsed • ${cutoffProvenanceLabel(provenance)}`;
      case "inferred":
        if (derivation === "from_wfc" && ratio != null) {
          return provenance === "missing"
            ? `Auto x${ratio}`
            : `Auto x${ratio} • ${cutoffProvenanceLabel(provenance)}`;
        }
        if (derivation === "from_rho" && ratio != null) {
          return provenance === "missing"
            ? `Auto /${ratio}`
            : `Auto /${ratio} • ${cutoffProvenanceLabel(provenance)}`;
        }
        return provenance === "missing"
          ? "Autocalculated"
          : `Auto • ${cutoffProvenanceLabel(provenance)}`;
      case "missing":
        return "Missing";
      default:
        return "";
    }
  }

  const selectedPseudoCutoffSummary = useMemo(() => {
    if (pseudopotentials.length === 0 && !ssspData) {
      return {
        maxWfc: 0,
        maxRho: 0,
        wfcStatus: "idle" as const,
        rhoStatus: "idle" as const,
        wfcProvenance: "missing" as const,
        rhoProvenance: "missing" as const,
        wfcDerivation: "missing" as const,
        rhoDerivation: "missing" as const,
        wfcRatio: null,
        rhoRatio: null,
        hasInferredCutoff: false,
        hasMissingCutoff: false,
      };
    }

    const elements = getUniqueElements();
    if (elements.length === 0) {
      return {
        maxWfc: 0,
        maxRho: 0,
        wfcStatus: "idle" as const,
        rhoStatus: "idle" as const,
        wfcProvenance: "missing" as const,
        rhoProvenance: "missing" as const,
        wfcDerivation: "missing" as const,
        rhoDerivation: "missing" as const,
        wfcRatio: null,
        rhoRatio: null,
        hasInferredCutoff: false,
        hasMissingCutoff: false,
      };
    }

    return summarizeSelectedPseudoCutoffs(
      elements,
      selectedPseudos,
      pseudoMetadataByFilename,
      ssspData,
    );
  }, [crystalData, optimizedStructures, pseudoMetadataByFilename, selectedPseudos, ssspData, structureSource]);

  const applyPseudopotentialPresetSelection = useCallback(async (preset: PseudopotentialPreset) => {
    const elements = getUniqueElements();
    setPseudoPresetWarnings({});
    if (elements.length === 0) {
      setSelectedPseudos({});
      return;
    }

    const filenames = await ensurePseudoFilenamesLoaded();
    const nextSelected: Record<string, string> = {};

    if (preset === "sssp") {
      const nextSsspData = await ensureSsspDataLoaded();
      const metadata = socPresetSelected ? await ensurePseudoMetadataLoaded() : null;
      const metadataByFilename = metadata
        ? Object.fromEntries(metadata.map((pseudo) => [pseudo.filename, pseudo]))
        : null;
      let hasCompleteSsspCoverage = Boolean(nextSsspData);
      if (nextSsspData) {
        for (const element of elements) {
          const ssspFilename = nextSsspData[element]?.filename;
          const matchedSsspFilename = findFilenameCaseInsensitive(filenames, ssspFilename);
          const ssspExists = Boolean(matchedSsspFilename);
          const ssspSupportsSoc = !socPresetSelected || (
            Boolean(matchedSsspFilename)
            && Boolean(metadataByFilename?.[matchedSsspFilename!]?.supports_soc)
          );
          if (!ssspExists || !ssspSupportsSoc) {
            hasCompleteSsspCoverage = false;
            break;
          }
        }
      }

      if (!hasCompleteSsspCoverage) {
        setPseudoError("SSSP pseudos not found");
        for (const element of elements) {
          const candidates = filenames.filter((filename) => pseudoMatchesElement(filename, element));
          const compatibleCandidates = socPresetSelected
            ? candidates.filter((filename) => metadataByFilename?.[filename]?.supports_soc)
            : candidates;
          const fallback = choosePreferredPseudoFilename(element, compatibleCandidates, nextSsspData);
          if (fallback) {
            nextSelected[element] = fallback;
          }
        }
        setSelectedPseudos(nextSelected);
        return;
      }

      for (const element of elements) {
        const ssspFilename = nextSsspData?.[element]?.filename;
        const matchedSsspFilename = findFilenameCaseInsensitive(filenames, ssspFilename);
        if (matchedSsspFilename) {
          nextSelected[element] = matchedSsspFilename;
        }
      }
      setPseudoError(null);
      setSelectedPseudos(nextSelected);
      return;
    }

    const metadata = await ensurePseudoMetadataLoaded();
    const metadataByFilename = Object.fromEntries(metadata.map((pseudo) => [pseudo.filename, pseudo]));
    const nextWarnings: Record<string, string> = {};

    for (const element of elements) {
      const matchingFilenames = filenames.filter((filename) => pseudoMatchesElement(filename, element));
      const socCompatibleFilenames = socPresetSelected
        ? matchingFilenames.filter((filename) => metadataByFilename[filename]?.supports_soc)
        : matchingFilenames;

      const typedMatches = socCompatibleFilenames.filter((filename) =>
        matchesPseudoPresetType(metadataByFilename[filename], preset),
      );
      const chosenTypedMatch = choosePreferredPseudoFilename(element, typedMatches, ssspData);
      if (chosenTypedMatch) {
        nextSelected[element] = chosenTypedMatch;
        continue;
      }

      const fallback = choosePreferredPseudoFilename(element, socCompatibleFilenames, ssspData);
      if (fallback) {
        nextSelected[element] = fallback;
        nextWarnings[element] = `No ${getPseudoPresetDisplayName(preset)} pseudo found. The calculation may still work.`;
      }
    }

    setPseudoError(null);
    setPseudoPresetWarnings(nextWarnings);
    setSelectedPseudos(nextSelected);
  }, [
    ensurePseudoFilenamesLoaded,
    ensurePseudoMetadataLoaded,
    ensureSsspDataLoaded,
    socPresetSelected,
    ssspData,
  ]);

  const handlePseudopotentialPresetClick = useCallback((preset: PseudopotentialPreset) => {
    setAutoCutoffFromPseudosEnabled(true);
    setSelectedPseudoPreset(preset);
    void applyPseudopotentialPresetSelection(preset).catch((e) => {
      console.error("Failed to apply pseudopotential preset:", e);
    });
  }, [applyPseudopotentialPresetSelection]);

  useEffect(() => {
    if (storedWizardSettings?.selectedPseudos && Object.keys(storedWizardSettings.selectedPseudos).length > 0) {
      return;
    }
    if (!resolvedPseudoDir.trim()) {
      return;
    }
    void applyPseudopotentialPresetSelection(selectedPseudoPreset).catch((e) => {
      console.error("Failed to apply pseudopotential preset:", e);
    });
  }, [applyPseudopotentialPresetSelection, resolvedPseudoDir, structureSource, crystalData, selectedOptimizedStructure?.calcId, socPresetSelected]);

  useEffect(() => {
    if (!resolvedPseudoDir.trim() || pseudoMetadataLoaded || isPseudoMetadataLoading) {
      return;
    }

    const selectedPslibraryPseudo = Object.values(selectedPseudos).some((filename) =>
      isPslibraryPseudoFilename(filename),
    );
    if (!selectedPslibraryPseudo) {
      return;
    }

    void Promise.all(
      Object.values(selectedPseudos)
        .filter((filename) => isPslibraryPseudoFilename(filename))
        .map((filename) => ensurePseudoMetadataForFilename(filename)),
    ).catch((e) => {
      console.error("Failed to load PSLibrary pseudopotential metadata:", e);
    });
  }, [
    ensurePseudoMetadataForFilename,
    isPseudoMetadataLoading,
    pseudoMetadataLoaded,
    resolvedPseudoDir,
    selectedPseudos,
  ]);

  // Keep cutoffs aligned with the currently selected pseudopotentials, including manual changes.
  useEffect(() => {
    if (!autoCutoffFromPseudosEnabled) return;
    if (pseudopotentials.length === 0 && !ssspData) return;

    setConfig((prev) => {
      let nextEcutwfc = prev.ecutwfc;
      let nextEcutrho = prev.ecutrho;
      let changed = false;

      if (selectedPseudoCutoffSummary.maxWfc > 0 && prev.ecutwfc !== selectedPseudoCutoffSummary.maxWfc) {
        nextEcutwfc = selectedPseudoCutoffSummary.maxWfc;
        changed = true;
      }

      if (selectedPseudoCutoffSummary.maxRho > 0 && prev.ecutrho !== selectedPseudoCutoffSummary.maxRho) {
        nextEcutrho = selectedPseudoCutoffSummary.maxRho;
        changed = true;
      }

      if (changed) {
        return {
          ...prev,
          ecutwfc: nextEcutwfc,
          ecutrho: nextEcutrho,
        };
      }
      return prev;
    });
  }, [autoCutoffFromPseudosEnabled, pseudopotentials, selectedPseudoCutoffSummary]);

  useEffect(() => {
    if (!crystalData || structureSource !== "cif") {
      setSymmetryTransform(null);
      setSymmetryError(null);
      return;
    }

    let cancelled = false;
    setSymmetryError(null);
    analyzeCrystalSymmetry(crystalData)
      .then((result) => {
        if (cancelled) return;
        setSymmetryTransform(result);
      })
      .catch((err) => {
        if (cancelled) return;
        console.error("Failed to analyze symmetry with spglib:", err);
        setSymmetryTransform(null);
        setSymmetryError(String(err));
      });

    return () => {
      cancelled = true;
    };
  }, [crystalData, structureSource]);

  async function handleImportCIF() {
    try {
      const selected = await open({
        multiple: false,
        filters: [{ name: "CIF Files", extensions: ["cif"] }],
        title: "Select CIF File",
      });

      if (selected && typeof selected === "string") {
        const content = await readTextFile(selected);
        const parsed = parseCIF(content);
        setCrystalData(parsed);
        setCifFilename(selected.split("/").pop() || "structure.cif");
        setCifContent(content);
        setError(null);
        setStep("configure");
      }
    } catch (e) {
      setError(`Failed to import CIF: ${e}`);
    }
  }

  function getUniqueElements(): string[] {
    if (selectedOptimizedStructure) {
      return [...new Set(selectedOptimizedStructure.structure.atoms.map((a) => getBaseElement(a.symbol)))];
    }
    if (!crystalData) return [];
    return [...new Set(crystalData.atom_sites.map((a) => getBaseElement(a.type_symbol)))];
  }

  function isSameMaterialAsClipboard(payload: ScfRunSettingsClipboardPayload): boolean {
    const currentElements = new Set(getUniqueElements().map((symbol) => getBaseElement(symbol)));
    if (currentElements.size === 0) return false;
    const copiedElements = new Set(payload.source.element_symbols.map((symbol) => getBaseElement(symbol)));
    if (copiedElements.size !== currentElements.size) return false;
    for (const symbol of copiedElements) {
      if (!currentElements.has(symbol)) {
        return false;
      }
    }
    return true;
  }

  async function refreshCopiedRunSettingsAvailability() {
    let available = hasStoredScfRunSettingsClipboardText();
    if (available) {
      setCopiedRunSettingsAvailable(true);
      return;
    }

    try {
      if (navigator.clipboard?.readText) {
        const clipboardText = await navigator.clipboard.readText();
        const parsed = parseScfRunSettingsClipboardText(clipboardText);
        if (parsed) {
          rememberScfRunSettingsClipboardText(clipboardText);
          available = true;
        }
      }
    } catch {
      // Ignore clipboard-read failures; availability falls back to in-app copied state.
    }
    setCopiedRunSettingsAvailable(available);
  }

  function applyCopiedScfSettings(payload: ScfRunSettingsClipboardPayload): {
    applied: boolean;
    skippedElementSpecific: boolean;
  } {
    const settings = asRecord(payload.settings);
    if (!settings) {
      return { applied: false, skippedElementSpecific: false };
    }

    const currentElements = new Set(getUniqueElements().map((symbol) => getBaseElement(symbol)));
    const sameMaterial = isSameMaterialAsClipboard(payload);
    const selectedPseudosForCurrentMaterial: Record<string, string> | null = sameMaterial
      ? (() => {
        const selectedPseudos = asRecord(settings.selected_pseudos);
        if (!selectedPseudos) return null;
        const nextPseudos: Record<string, string> = {};
        for (const [element, filename] of Object.entries(selectedPseudos)) {
          const normalized = getBaseElement(element);
          if (!currentElements.has(normalized)) continue;
          const text = String(filename ?? "").trim();
          if (!text) continue;
          nextPseudos[normalized] = text;
        }
        return nextPseudos;
      })()
      : null;
    const filterByCurrentElements = (value: unknown): Record<string, number> => {
      const raw = asRecord(value);
      if (!raw) return {};
      const next: Record<string, number> = {};
      for (const [element, amount] of Object.entries(raw)) {
        const normalized = getBaseElement(element);
        if (!currentElements.has(normalized)) continue;
        const numeric = asFiniteNumber(amount);
        if (numeric == null) continue;
        next[normalized] = numeric;
      }
      return next;
    };
    const filterManifoldsByCurrentElements = (value: unknown): Record<string, string> => {
      const raw = asRecord(value);
      if (!raw) return {};
      const next: Record<string, string> = {};
      for (const [element, manifold] of Object.entries(raw)) {
        const normalized = getBaseElement(element);
        if (!currentElements.has(normalized)) continue;
        const text = String(manifold ?? "").trim();
        if (!text) continue;
        next[normalized] = text;
      }
      return next;
    };

    setConfig((prev) => {
      const next: SCFConfig = { ...prev };

      const calculationMode = settings.calculation_mode;
      if (
        isOptimizationWizard
        && (calculationMode === "scf" || calculationMode === "relax" || calculationMode === "vcrelax")
      ) {
        next.calculation = calculationMode;
      } else {
        next.calculation = "scf";
      }

      const occupations = settings.occupations;
      if (
        occupations === "smearing" || occupations === "tetrahedra" || occupations === "tetrahedra_lin"
        || occupations === "tetrahedra_opt" || occupations === "fixed" || occupations === "from_input"
      ) {
        next.occupations = occupations;
      }

      const smearing = settings.smearing;
      if (
        smearing === "gaussian" || smearing === "methfessel-paxton"
        || smearing === "marzari-vanderbilt" || smearing === "fermi-dirac" || smearing === "cold"
      ) {
        next.smearing = smearing;
      }

      const simpleNumberFields: Array<[keyof SCFConfig, unknown]> = [
        ["ecutwfc", settings.ecutwfc],
        ["ecutrho", settings.ecutrho],
        ["degauss", settings.degauss],
        ["conv_thr", settings.conv_thr],
        ["mixing_beta", settings.mixing_beta],
        ["forc_conv_thr", settings.forc_conv_thr],
        ["etot_conv_thr", settings.etot_conv_thr],
        ["press", settings.press],
        ["electron_maxstep", settings.electron_maxstep],
        ["mixing_ndim", settings.mixing_ndim],
        ["tot_charge", settings.tot_charge],
      ];
      for (const [field, value] of simpleNumberFields) {
        const parsed = asFiniteNumber(value);
        if (parsed != null) {
          (next[field] as number) = parsed;
        }
      }

      const nbndValue = settings.nbnd;
      if (nbndValue == null || nbndValue === "") {
        next.nbnd = null;
      } else {
        const parsedNbnd = asFiniteNumber(nbndValue);
        if (parsedNbnd != null) {
          next.nbnd = Math.max(1, Math.round(parsedNbnd));
        }
      }

      const kgrid = asNumberTriplet(settings.kgrid);
      if (kgrid) {
        next.kgrid = toIntegerTriplet(kgrid);
      }
      const kgridOffset = asNumberTriplet(settings.kgrid_offset);
      if (kgridOffset) {
        next.kgrid_offset = toOffsetTriplet(kgridOffset);
      }

      const nspin = asFiniteNumber(settings.nspin);
      if (nspin === 1 || nspin === 2 || nspin === 4) {
        next.nspin = nspin;
      }

      const boolFields: Array<[keyof SCFConfig, unknown]> = [
        ["noncolin", settings.noncolin],
        ["lspinorb", settings.lspinorb],
        ["tprnfor", settings.tprnfor],
        ["tstress", settings.tstress],
      ];
      for (const [field, value] of boolFields) {
        const parsed = asBoolean(value);
        if (parsed != null) {
          (next[field] as boolean) = parsed;
        }
      }

      const totMag = settings.tot_magnetization;
      if (totMag == null || totMag === "") {
        next.tot_magnetization = null;
      } else {
        const parsedTotMag = asFiniteNumber(totMag);
        if (parsedTotMag != null) {
          next.tot_magnetization = parsedTotMag;
        }
      }

      const constrainedMag = settings.constrained_magnetization;
      if (
        constrainedMag === "none" || constrainedMag === "total"
        || constrainedMag === "atomic" || constrainedMag === "total direction"
        || constrainedMag === "atomic direction"
      ) {
        next.constrained_magnetization = constrainedMag;
      }

      const vdw = settings.vdw_corr;
      if (
        vdw === "none" || vdw === "grimme-d2" || vdw === "grimme-d3"
        || vdw === "ts-vdw" || vdw === "xdm" || vdw === "dft-d"
      ) {
        next.vdw_corr = vdw;
      }

      const assumeIsolated = settings.assume_isolated;
      if (
        assumeIsolated === "none" || assumeIsolated === "makov-payne"
        || assumeIsolated === "martyna-tuckerman" || assumeIsolated === "esm"
        || assumeIsolated === "2D"
      ) {
        next.assume_isolated = assumeIsolated;
      }

      const verbosity = settings.verbosity;
      if (verbosity === "low" || verbosity === "high") {
        next.verbosity = verbosity;
      }

      const diskIo = settings.disk_io;
      if (diskIo === "low" || diskIo === "medium" || diskIo === "high" || diskIo === "nowf") {
        next.disk_io = diskIo;
      }

      const mixingMode = settings.mixing_mode;
      if (mixingMode === "plain" || mixingMode === "TF" || mixingMode === "local-TF") {
        next.mixing_mode = mixingMode;
      }

      const diagonalization = settings.diagonalization;
      if (
        diagonalization === "david" || diagonalization === "cg" || diagonalization === "ppcg"
        || diagonalization === "paro" || diagonalization === "rmm-davidson"
      ) {
        next.diagonalization = diagonalization;
      }

      const startingpot = settings.startingpot;
      if (startingpot === "atomic" || startingpot === "file") {
        next.startingpot = startingpot;
      }

      const startingwfc = settings.startingwfc;
      if (
        startingwfc === "atomic" || startingwfc === "atomic+random"
        || startingwfc === "random" || startingwfc === "file"
      ) {
        next.startingwfc = startingwfc;
      }

      const inputDft = settings.input_dft;
      if (typeof inputDft === "string") {
        next.input_dft = inputDft;
      }

      if (sameMaterial) {
        const startingMag = filterByCurrentElements(settings.starting_magnetization);
        if (Object.keys(startingMag).length > 0) {
          next.starting_magnetization = startingMag;
        }
        const startingMagTheta = filterByCurrentElements(
          settings.starting_magnetization_theta ?? settings.starting_magnetization_angle1 ?? settings.theta ?? settings.angle1,
        );
        if (Object.keys(startingMagTheta).length > 0) {
          next.starting_magnetization_theta = startingMagTheta;
        }
        const startingMagPhi = filterByCurrentElements(
          settings.starting_magnetization_phi ?? settings.starting_magnetization_angle2 ?? settings.phi ?? settings.angle2,
        );
        if (Object.keys(startingMagPhi).length > 0) {
          next.starting_magnetization_phi = startingMagPhi;
        }

        const ldaPlusU = asBoolean(settings.lda_plus_u);
        if (ldaPlusU != null) {
          next.lda_plus_u = ldaPlusU;
        }

        const hubbardProjector = settings.hubbard_projector;
        if (
          hubbardProjector === "atomic" || hubbardProjector === "ortho-atomic"
          || hubbardProjector === "norm-atomic" || hubbardProjector === "wf"
          || hubbardProjector === "pseudo"
        ) {
          next.hubbard_projector = hubbardProjector;
        }

        const hubbardKind = asFiniteNumber(settings.hubbard_formulation);
        if (hubbardKind === 0 || hubbardKind === 1 || hubbardKind === 2) {
          next.lda_plus_u_kind = hubbardKind;
        }

        const hubbardManifold = filterManifoldsByCurrentElements(settings.hubbard_manifold);
        if (Object.keys(hubbardManifold).length > 0) {
          next.hubbard_manifold = hubbardManifold;
        }

        const hubbardU = filterByCurrentElements(settings.hubbard_u);
        if (Object.keys(hubbardU).length > 0) {
          next.hubbard_u = hubbardU;
        }

        const hubbardJ = filterByCurrentElements(settings.hubbard_j);
        if (Object.keys(hubbardJ).length > 0) {
          next.hubbard_j = hubbardJ;
        }
      }

      return next;
    });

    const copiedConvThr = asFiniteNumber(settings.conv_thr);
    setConvThrInput(copiedConvThr != null ? String(copiedConvThr) : String(config.conv_thr));
    if (selectedPseudosForCurrentMaterial) {
      setSelectedPseudos(selectedPseudosForCurrentMaterial);
      setPseudoPresetWarnings({});
    }
    setSelectedPreset(null);

    return {
      applied: true,
      skippedElementSpecific: !sameMaterial,
    };
  }

  async function handlePasteFromCopiedRun() {
    setError(null);
    setCopiedRunSettingsMessage(null);
    let payload: ScfRunSettingsClipboardPayload | null = null;
    let clipboardText: string | null = null;

    try {
      if (navigator.clipboard?.readText) {
        clipboardText = await navigator.clipboard.readText();
        payload = parseScfRunSettingsClipboardText(clipboardText);
      }
    } catch {
      // Ignore clipboard-read failures and fall back to in-app stored copy.
    }

    if (!payload) {
      const storedText = getStoredScfRunSettingsClipboardText();
      payload = parseScfRunSettingsClipboardText(storedText);
      clipboardText = storedText;
    }

    if (!payload) {
      setCopiedRunSettingsAvailable(false);
      setCopiedRunSettingsMessage("No copied run settings were found. Copy an existing run first.");
      return;
    }

    if (clipboardText) {
      rememberScfRunSettingsClipboardText(clipboardText);
    }

    const applyOutcome = applyCopiedScfSettings(payload);
    if (!applyOutcome.applied) {
      setCopiedRunSettingsMessage("Copied settings could not be parsed.");
      return;
    }
    setCopiedRunSettingsAvailable(true);
    if (applyOutcome.skippedElementSpecific) {
      setCopiedRunSettingsMessage("Pasted settings. Element-specific options were skipped because the material differs.");
    } else {
      setCopiedRunSettingsMessage("Pasted settings from copied run.");
    }
  }

  const hubbardDefaultElements = getUniqueElements();
  const hubbardDefaultElementKey = hubbardDefaultElements.join("|");
  const hubbardRecommendations = useMemo(
    () => getHubbardRecommendations(hubbardDefaultElements),
    [hubbardDefaultElementKey],
  );
  const hubbardRecommendationText = hubbardRecommendations
    .map((entry) => `${entry.element}-${entry.manifold}`)
    .join(", ");
  const isHubbardRecommended = hubbardRecommendations.length > 0;

  useEffect(() => {
    if (hubbardDefaultElements.length === 0) return;

    setConfig((prev) => {
      let changed = false;
      const nextManifolds = { ...prev.hubbard_manifold };

      for (const element of hubbardDefaultElements) {
        const current = normalizeHubbardManifold(nextManifolds[element] || "");
        const defaultManifold = getDefaultHubbardManifold(element);
        const previousOuterShellDefault = getOutermostOccupiedOrbitalManifold(element);

        if (!defaultManifold) {
          if (current && current === previousOuterShellDefault) {
            delete nextManifolds[element];
            changed = true;
          }
          continue;
        }

        if (!current || current === previousOuterShellDefault) {
          nextManifolds[element] = defaultManifold;
          changed = true;
        }
      }

      return changed ? { ...prev, hubbard_manifold: nextManifolds } : prev;
    });
  }, [hubbardDefaultElementKey]);

  useEffect(() => {
    if (!isHubbardRecommended) return;
    setExpandedSections((prev) => (prev.dftu ? prev : { ...prev, dftu: true }));
  }, [isHubbardRecommended, hubbardDefaultElementKey]);

  useEffect(() => {
    if (hubbardRecommendations.length === 0) return;
    setConfig((prev) => {
      let changed = false;
      const nextU = { ...prev.hubbard_u };

      for (const recommendation of hubbardRecommendations) {
        if (manuallyEditedHubbardU[recommendation.element]) continue;
        const current = nextU[recommendation.element];
        if (current != null && current !== 0) continue;
        nextU[recommendation.element] = GENERAL_HUBBARD_U_GUESS_EV;
        changed = true;
      }

      if (!changed) return prev;
      return { ...prev, hubbard_u: nextU };
    });
  }, [hubbardDefaultElementKey, manuallyEditedHubbardU]);

  useEffect(() => {
    if (config.lda_plus_u_kind !== 1 || hubbardDefaultElements.length === 0) return;
    const nextJLabels: Record<string, string> = {};

    setConfig((prev) => {
      let changed = false;
      const nextJ = { ...prev.hubbard_j };

      for (const element of hubbardDefaultElements) {
        if (manuallyEditedHubbardJ[element]) continue;
        const defaultJ = getHundJDefaultEv(element);
        if (defaultJ == null) continue;
        const current = nextJ[element];
        if (current != null && current !== 0) continue;
        nextJ[element] = defaultJ;
        nextJLabels[element] = `${HUBBARD_J_SOURCE}: ${defaultJ.toFixed(3)} eV`;
        changed = true;
      }

      if (!changed) return prev;
      return { ...prev, hubbard_j: nextJ };
    });

    if (Object.keys(nextJLabels).length > 0) {
      setHubbardJDefaultLabels((prev) => ({ ...prev, ...nextJLabels }));
    }
  }, [config.lda_plus_u_kind, hubbardDefaultElementKey, manuallyEditedHubbardJ]);

  function canRun(): boolean {
    if (!crystalData) return false;
    const elements = getUniqueElements();
    return elements.every((el) => selectedPseudos[el]);
  }

  function getHubbardManifoldForElement(element: string): string | null {
    const explicit = normalizeHubbardManifold(config.hubbard_manifold[element] || "");
    return explicit || getDefaultHubbardManifold(element);
  }

  function applyCalculatedHubbardU(element: string) {
    const manifold = getHubbardManifoldForElement(element);
    if (!manifold) return;

    const lrtValue = getLatestHubbardLrtValue(element, manifold, (initialCif?.calculations ?? []) as any[]);
    if (!lrtValue) return;

    setConfig((prev) => ({
      ...prev,
      hubbard_u: { ...prev.hubbard_u, [element]: lrtValue.value },
    }));
    setManuallyEditedHubbardU((prev) => ({ ...prev, [element]: true }));
  }

  function buildHubbardSettings(elements: string[]) {
    if (!config.lda_plus_u) return undefined;

    const parameters: HubbardParameter[] = [];
    for (const element of elements) {
      const uValue = config.hubbard_u[element] ?? 0;
      const jValue = config.hubbard_j[element] ?? 0;
      const hasJValue = config.lda_plus_u_kind > 0 && jValue !== 0;
      if (uValue === 0 && !hasJValue) continue;

      const manifold = getHubbardManifoldForElement(element);
      if (!manifold) {
        throw new Error(`Missing Hubbard manifold for ${element}. Enter a value such as 3d or 4f.`);
      }

      const target = `${element}-${manifold}`;
      parameters.push({ parameter: "U", manifold: target, value: uValue });
      if (hasJValue) {
        parameters.push({
          parameter: config.lda_plus_u_kind === 2 ? "J0" : "J",
          manifold: target,
          value: jValue,
        });
      }
    }

    if (parameters.length === 0) return undefined;
    return {
      projector: config.hubbard_projector,
      parameters,
    };
  }

  async function buildScfTaskPlan(): Promise<ScfTaskPlan> {
    if (!crystalData || !canRun()) {
      throw new Error("Missing structure or pseudopotential selections.");
    }

    const elements = getUniqueElements();
    const cifStructure = buildCifConventionalStructure(crystalData);
    const sourceStructure: SavedStructureData = selectedOptimizedStructure
      ? {
          ...selectedOptimizedStructure.structure,
          cell_parameters: selectedOptimizedStructure.structure.cell_parameters || cifStructure.cell_parameters,
          cell_units: selectedOptimizedStructure.structure.cell_units || cifStructure.cell_units,
        }
      : cifStructure;
    const sourceDescriptor = selectedOptimizedStructure
      ? { type: "optimization" as const, calc_id: selectedOptimizedStructure.calcId }
      : { type: "cif" as const };
    let magnetismViewerStructure: SavedStructureData = sourceStructure;

    let resolvedSymmetry = symmetryTransform;
    if (!selectedOptimizedStructure && !resolvedSymmetry) {
      try {
        resolvedSymmetry = await analyzeCrystalSymmetry(crystalData);
        setSymmetryTransform(resolvedSymmetry);
        setSymmetryError(null);
      } catch (err) {
        setSymmetryError(String(err));
      }
    }

    const canUseSymmetryPrimitive =
      !selectedOptimizedStructure &&
      resolvedSymmetry !== null &&
      resolvedSymmetry.standardizedPrimitiveAtoms.length > 0;

    if (canUseSymmetryPrimitive && resolvedSymmetry) {
      console.log(
        `Using spglib primitive cell: nat=${resolvedSymmetry.standardizedPrimitiveAtoms.length}, spacegroup=${resolvedSymmetry.spacegroupNumber}`,
      );
    }

    // Build species list (same for both primitive and conventional)
    const speciesList = elements.map((el) => ({
      symbol: el,
      mass: ELEMENT_MASSES[el] || 1.0,
      pseudopotential: selectedPseudos[el],
      starting_magnetization: config.starting_magnetization[el] || 0,
      ...(config.nspin === 4
        ? {
          theta: config.starting_magnetization_theta[el] ?? 0,
          phi: config.starting_magnetization_phi[el] ?? 0,
        }
        : {}),
    }));

    // Build system configuration based on cell type
    const normalizedInputDft = config.input_dft.trim();
    const hubbardSettings = buildHubbardSettings(elements);

    const systemConfig: any = {
      // Common properties
      species: speciesList,
      position_units: sourceStructure.position_units || "crystal",
      ecutwfc: config.ecutwfc,
      ecutrho: config.ecutrho,
      // Electronic structure
      occupations: config.occupations,
      smearing: config.occupations === "smearing" ? config.smearing : undefined,
      degauss: config.occupations === "smearing" ? config.degauss : undefined,
      nbnd: config.nbnd,
      tot_charge: config.tot_charge !== 0 ? config.tot_charge : undefined,
      // Magnetism
      nspin: config.nspin,
      noncolin: config.noncolin || undefined,
      lspinorb: config.lspinorb || undefined,
      tot_magnetization: config.nspin === 2 && config.tot_magnetization !== null ? config.tot_magnetization : undefined,
      constrained_magnetization: config.constrained_magnetization !== "none" ? config.constrained_magnetization : undefined,
      // DFT+U
      hubbard: hubbardSettings,
      // Van der Waals
      vdw_corr: config.vdw_corr !== "none" ? config.vdw_corr : undefined,
      // Isolated systems
      assume_isolated: config.assume_isolated !== "none" ? config.assume_isolated : undefined,
      // XC functional override
      input_dft: normalizedInputDft || undefined,
    };

    // Add cell-specific properties
    if (selectedOptimizedStructure) {
      systemConfig.ibrav = "free";
      systemConfig.celldm = null;
      systemConfig.cell_parameters = sourceStructure.cell_parameters;
      systemConfig.cell_units = sourceStructure.cell_units || "angstrom";
      systemConfig.atoms = sourceStructure.atoms.map((atom) => ({
        symbol: getBaseElement(atom.symbol),
        position: atom.position,
        if_pos: [true, true, true],
      }));
      magnetismViewerStructure = {
        ...sourceStructure,
        atoms: systemConfig.atoms.map((atom: { symbol: string; position: [number, number, number] }) => ({
          symbol: atom.symbol,
          position: atom.position,
        })),
      };
    } else if (canUseSymmetryPrimitive && resolvedSymmetry) {
      const inferredBravais = inferQeBravaisCellFromCif(crystalData, resolvedSymmetry);
      if (inferredBravais) {
        systemConfig.ibrav = inferredBravais.ibrav;
        systemConfig.celldm = inferredBravais.celldm;
        systemConfig.cell_parameters = null;
        systemConfig.cell_units = null;
        systemConfig.atoms = inferredBravais.atoms.map((atom) => ({
          symbol: atom.symbol,
          position: atom.position,
          if_pos: [true, true, true],
        }));
        magnetismViewerStructure = {
          position_units: "crystal",
          atoms: inferredBravais.atoms.map((atom) => ({
            symbol: atom.symbol,
            position: atom.position,
          })),
          cell_parameters: resolvedSymmetry.standardizedPrimitiveLattice,
          cell_units: "angstrom",
        };
      } else {
        systemConfig.ibrav = "free";
        systemConfig.celldm = null;
        systemConfig.cell_parameters = resolvedSymmetry.standardizedPrimitiveLattice;
        systemConfig.cell_units = "angstrom";
        systemConfig.atoms = resolvedSymmetry.standardizedPrimitiveAtoms.map((atom) => ({
          symbol: atom.symbol,
          position: atom.position,
          if_pos: [true, true, true],
        }));
        magnetismViewerStructure = {
          position_units: "crystal",
          atoms: resolvedSymmetry.standardizedPrimitiveAtoms.map((atom) => ({
            symbol: atom.symbol,
            position: atom.position,
          })),
          cell_parameters: resolvedSymmetry.standardizedPrimitiveLattice,
          cell_units: "angstrom",
        };
      }
    } else {
      // Conventional cell with ibrav=0 fallback when primitive extraction is unavailable.
      systemConfig.ibrav = "free";
      systemConfig.celldm = null;
      systemConfig.cell_parameters = buildConventionalLatticeFromCrystalData(crystalData);
      systemConfig.cell_units = "angstrom";
      systemConfig.atoms = crystalData.atom_sites.map((site) => ({
        symbol: getBaseElement(site.type_symbol),
        position: [site.fract_x, site.fract_y, site.fract_z],
        if_pos: [true, true, true],
      }));
      magnetismViewerStructure = {
        position_units: "crystal",
        atoms: systemConfig.atoms.map((atom: { symbol: string; position: [number, number, number] }) => ({
          symbol: atom.symbol,
          position: atom.position,
        })),
        cell_parameters: systemConfig.cell_parameters,
        cell_units: "angstrom",
      };
    }

    // Build the calculation configuration with all options
    const pseudoDir = resolvedPseudoDir;
    if (!pseudoDir.trim()) {
      throw new Error(
        isHpcMode
          ? "Remote pseudopotential directory is not configured for the active HPC profile."
          : "Local pseudopotential directory is not configured.",
      );
    }

    const calculation: any = {
      calculation: config.calculation,
      prefix: "qcortado_scf",
      outdir: "./tmp",
      pseudo_dir: pseudoDir,
      verbosity: config.verbosity,
      tprnfor: config.tprnfor,
      tstress: config.tstress,
      disk_io: config.disk_io,
      system: systemConfig,
      kpoints: {
        type: "automatic",
        grid: config.kgrid,
        offset: config.kgrid_offset,
      },
      // SCF convergence settings
      conv_thr: config.conv_thr,
      electron_maxstep: config.electron_maxstep,
      mixing_mode: config.mixing_mode,
      mixing_beta: config.mixing_beta,
      mixing_ndim: config.mixing_ndim,
      diagonalization: config.diagonalization,
      startingpot: config.startingpot,
      startingwfc: config.startingwfc,
      // Relaxation settings (only used for relax/vcrelax)
      forc_conv_thr: (config.calculation === "relax" || config.calculation === "vcrelax") ? config.forc_conv_thr : null,
      etot_conv_thr: (config.calculation === "relax" || config.calculation === "vcrelax") ? config.etot_conv_thr : null,
      // vcrelax specific
      press: config.calculation === "vcrelax" ? config.press : null,
    };

    // Generate and display the input file
    const inputText = await invoke<string>("generate_input", { calculation });
    const taskLabel = crystalData?.formula_sum || cifFilename || "SCF";

    return {
      calculation,
      inputText,
      taskLabel,
      taskParams: {
        calculation,
        workingDir: WORK_DIR,
        mpiConfig: !isHpcMode && mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
        executionTarget: buildExecutionTarget(
          executionMode,
          activeHpcProfile?.id ?? null,
          isHpcMode ? hpcResources : null,
          false,
        ),
      },
      sourceStructure,
      magnetismViewerStructure,
      sourceDescriptor,
    };
  }

  async function handleRun() {
    if (!crystalData || !canRun()) return;
    if (hasBlockingExternalTask) {
      setError("Another local task is currently running. Queue this task or wait for completion.");
      return;
    }

    setIsRunning(true);
    setOutput("");
    setOutputLineCount(0);
    setResult(null);
    setResultSaved(false);
    setProgress(defaultProgressState("SCF iterations"));
    setStep("run");

    // Track calculation start time
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);

    // Run the calculation as a background task

    try {
      const plan = await buildScfTaskPlan();
      const { inputText, taskLabel, sourceStructure, magnetismViewerStructure, sourceDescriptor } = plan;
      const draftCalcData = buildCalculationData(
        {
          converged: false,
          total_energy: null,
          fermi_energy: null,
          n_scf_steps: null,
          wall_time_seconds: null,
          raw_output: "",
        },
        startTime,
        startTime,
        inputText,
        sourceStructure,
        sourceDescriptor,
        magnetismViewerStructure,
      );
      const queueCalcType: "scf" | "optimization" = draftCalcData.calc_type === "optimization" ? "optimization" : "scf";
      const runSaveSpec = projectContext
        ? {
            projectId: projectContext.projectId,
            cifId: projectContext.cifId,
            workingDir: WORK_DIR,
            calcType: queueCalcType,
            parameters: draftCalcData.parameters,
            tags: draftCalcData.tags,
            inputContent: inputText,
          }
        : null;

      setRunSourceStructure(sourceStructure);
      setRunMagnetismViewerStructure(magnetismViewerStructure);
      setRunSourceDescriptor(sourceDescriptor);
      setGeneratedInput(inputText);
      setOutput(`=== Generated Input ===\n${inputText}\n\n=== Running pw.x ===\n`);

      // Start as a background task
      const taskId = await taskContext.startTask("scf", plan.taskParams, `SCF - ${taskLabel}`, runSaveSpec);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Calculation failed");
      }

      const calcResult = finalTask.result;

      // Track calculation end time
      const endTime = new Date().toISOString();
      setCalcEndTime(endTime);

      setResult(calcResult);
      setOutput((prev) => prev + "\n=== Calculation Complete ===\n");
      setStep("results");
      if (runSaveSpec && !isHpcMode && finalTask.hpc.backend !== "hpc") {
        setIsSaving(true);
        setProgress((prev) => ({
          ...prev,
          status: "running",
          percent: null,
          phase: "Saving to project",
        }));
        const saveQueueItem = await taskContext.waitForQueueItemCompletion(taskId);
        setIsSaving(false);
        if (saveQueueItem?.status === "completed") {
          setResultSaved(true);
          setProgress((prev) => ({
            ...prev,
            status: "complete",
            percent: 100,
            phase: "Complete",
          }));
        } else {
          const queueError = saveQueueItem?.error || "Unknown save failure";
          setError(`Failed to auto-save calculation: ${queueError}`);
          setProgress((prev) => ({
            ...prev,
            status: "error",
            percent: null,
            phase: "Save Error",
          }));
        }
      } else if (projectContext && (isHpcMode || finalTask.hpc.backend === "hpc")) {
        setIsSaving(true);
        setProgress((prev) => ({
          ...prev,
          status: "running",
          percent: null,
          phase: "Saving to project",
        }));
        try {
          const hpcCalcData = buildCalculationData(
            calcResult,
            startTime,
            endTime,
            inputText,
            sourceStructure,
            sourceDescriptor,
            magnetismViewerStructure,
            finalTask.hpc,
          );
          await invoke("save_calculation", {
            projectId: projectContext.projectId,
            cifId: projectContext.cifId,
            calcData: hpcCalcData,
            workingDir: finalTask.hpc.local_sync_dir ?? WORK_DIR,
          });
          setResultSaved(true);
          setProgress((prev) => ({
            ...prev,
            status: "complete",
            percent: 100,
            phase: "Complete",
          }));
        } catch (saveError) {
          setError(`Failed to auto-save calculation: ${saveError}`);
          setProgress((prev) => ({
            ...prev,
            status: "error",
            percent: null,
            phase: "Save Error",
          }));
        } finally {
          setIsSaving(false);
        }
      } else {
        setProgress((prev) => ({
          ...prev,
          status: "complete",
          percent: 100,
          phase: "Complete",
        }));
      }
    } catch (e) {
      setError(`Calculation failed: ${e}`);
      setOutput((prev) => prev + `\n\nERROR: ${e}`);
      setProgress((prev) => ({
        ...prev,
        status: "error",
        percent: null,
        phase: "Error",
      }));
    } finally {
      setIsRunning(false);
    }
  }

  async function handleQueue() {
    if (!crystalData || !canRun()) return;
    if (isHpcMode) {
      setError("Queueing is unavailable in HPC mode. Submit directly to Andromeda.");
      return;
    }

    if (!projectContext) {
      setError("Queueing is available when running from a project.");
      return;
    }

    try {
      const plan = await buildScfTaskPlan();
      const nowIso = new Date().toISOString();
      const draftCalcData = buildCalculationData(
        {
          converged: false,
          total_energy: null,
          fermi_energy: null,
          n_scf_steps: null,
          wall_time_seconds: null,
          raw_output: "",
        },
        nowIso,
        nowIso,
        plan.inputText,
        plan.sourceStructure,
        plan.sourceDescriptor,
        plan.magnetismViewerStructure,
      );
      const queueCalcType = draftCalcData.calc_type === "optimization" ? "optimization" : "scf";
      const queueLabel = queueCalcType === "optimization"
        ? `${config.calculation === "vcrelax" ? "VC-Relax" : "Relax"} - ${plan.taskLabel}`
        : `SCF - ${plan.taskLabel}`;

      taskContext.enqueueTask(
        "scf",
        plan.taskParams,
        queueLabel,
        {
          projectId: projectContext.projectId,
          cifId: projectContext.cifId,
          workingDir: WORK_DIR,
          calcType: queueCalcType,
          parameters: draftCalcData.parameters,
          tags: draftCalcData.tags,
          inputContent: plan.inputText,
        },
      );
      setError(null);
    } catch (e) {
      setError(`Failed to queue calculation: ${e}`);
    }
  }

  function calculateMetricsFromVectors(
    v1: [number, number, number],
    v2: [number, number, number],
    v3: [number, number, number],
  ): DisplayCellMetrics {
    const norm = (v: [number, number, number]) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    const dot = (u: [number, number, number], v: [number, number, number]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const angle = (u: [number, number, number], v: [number, number, number]) => {
      const denom = norm(u) * norm(v);
      if (denom === 0) return 0;
      const cosine = dot(u, v) / denom;
      const clamped = Math.max(-1, Math.min(1, cosine));
      return (Math.acos(clamped) * 180) / Math.PI;
    };

    return {
      a: norm(v1),
      b: norm(v2),
      c: norm(v3),
      alpha: angle(v2, v3),
      beta: angle(v1, v3),
      gamma: angle(v1, v2),
    };
  }

  function isQeUnit(value: string): value is QePositionUnit {
    return value === "alat" || value === "bohr" || value === "angstrom" || value === "crystal";
  }

  function parseUnitFromHeader(line: string): QePositionUnit | null {
    const match = line.match(/[\(\{]\s*([a-zA-Z]+)/);
    if (!match) return null;
    const unit = match[1].toLowerCase();
    return isQeUnit(unit) ? unit : null;
  }

  function parseAlatFromHeader(line: string): number | null {
    const match = line.match(/[\(\{]\s*alat\s*=\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)/i);
    if (!match) return null;
    const alatValue = Number.parseFloat(match[1]);
    if (!Number.isFinite(alatValue) || alatValue <= 0) return null;
    return alatValue;
  }

  function parseTriplet(line: string): [number, number, number] | null {
    const numbers = line.match(/[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?/g);
    if (!numbers || numbers.length < 3) return null;
    const first = Number.parseFloat(numbers[0]);
    const second = Number.parseFloat(numbers[1]);
    const third = Number.parseFloat(numbers[2]);
    if ([first, second, third].some((value) => Number.isNaN(value))) return null;
    return [first, second, third];
  }

  function parseCellBlock(lines: string[], index: number): {
    cellParameters: [[number, number, number], [number, number, number], [number, number, number]] | null;
    cellUnits: QePositionUnit | null;
    cellAlat: number | null;
  } {
    const vectors: [number, number, number][] = [];
    for (let i = index + 1; i < Math.min(lines.length, index + 8); i += 1) {
      const vector = parseTriplet(lines[i]);
      if (!vector) break;
      vectors.push(vector);
      if (vectors.length === 3) break;
    }
    return {
      cellParameters: vectors.length === 3
        ? [vectors[0], vectors[1], vectors[2]]
        : null,
      cellUnits: parseUnitFromHeader(lines[index]),
      cellAlat: parseAlatFromHeader(lines[index]),
    };
  }

  function parseAtomsBlock(lines: string[], index: number): {
    atoms: { symbol: string; position: [number, number, number] }[];
    positionUnits: QePositionUnit;
  } {
    const atoms: { symbol: string; position: [number, number, number] }[] = [];
    const parsedUnit = parseUnitFromHeader(lines[index]);
    const positionUnits: QePositionUnit = parsedUnit || "crystal";

    for (let i = index + 1; i < lines.length; i += 1) {
      const line = lines[i].trim();
      if (!line) break;
      if (
        line.startsWith("CELL_PARAMETERS") ||
        line.startsWith("ATOMIC_POSITIONS") ||
        line.startsWith("End final coordinates")
      ) {
        break;
      }
      const match = line.match(/^([A-Za-z]{1,3}[A-Za-z0-9]*)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)/);
      if (!match) continue;
      const x = Number.parseFloat(match[2]);
      const y = Number.parseFloat(match[3]);
      const z = Number.parseFloat(match[4]);
      if ([x, y, z].some((value) => Number.isNaN(value))) continue;
      atoms.push({
        symbol: getBaseElement(match[1]),
        position: [x, y, z],
      });
    }

    return { atoms, positionUnits };
  }

  function buildCifConventionalStructure(data: CrystalData): SavedStructureData {
    return {
      position_units: "crystal",
      cell_units: "angstrom",
      cell_parameters: buildConventionalLatticeFromCrystalData(data),
      atoms: data.atom_sites.map((site) => ({
        symbol: getBaseElement(site.type_symbol),
        position: [site.fract_x, site.fract_y, site.fract_z] as [number, number, number],
      })),
    };
  }

  function convertCellToAngstrom(
    cellParameters: [[number, number, number], [number, number, number], [number, number, number]] | null,
    cellUnits: QePositionUnit | null,
    alatInBohr: number | null,
  ): {
    cellParameters: [[number, number, number], [number, number, number], [number, number, number]] | null;
    cellUnits: QePositionUnit | null;
  } {
    if (!cellParameters || !cellUnits) {
      return { cellParameters, cellUnits };
    }

    const BOHR_TO_ANGSTROM = 0.529177210903;

    if (cellUnits === "angstrom") {
      return { cellParameters, cellUnits };
    }

    let scale: number | null = null;
    if (cellUnits === "bohr") {
      scale = BOHR_TO_ANGSTROM;
    } else if (cellUnits === "alat" && alatInBohr && Number.isFinite(alatInBohr) && alatInBohr > 0) {
      // QE reports "alat" in Bohr, so convert vector components to Angstrom.
      scale = alatInBohr * BOHR_TO_ANGSTROM;
    }

    if (scale === null) {
      return { cellParameters, cellUnits };
    }

    const scaleValue = scale;
    const scaled = cellParameters.map((vector) => [
      vector[0] * scaleValue,
      vector[1] * scaleValue,
      vector[2] * scaleValue,
    ]) as [[number, number, number], [number, number, number], [number, number, number]];

    return {
      cellParameters: scaled,
      cellUnits: "angstrom",
    };
  }

  function extractOptimizedStructure(rawOutput: string, fallback: SavedStructureData): SavedStructureData {
    const lines = rawOutput.split(/\r?\n/);
    const begin = lines.findIndex((line) => line.includes("Begin final coordinates"));
    const end = lines.findIndex((line) => line.includes("End final coordinates"));
    const scope = begin >= 0 && end > begin ? lines.slice(begin, end + 1) : lines;

    let lastCell = fallback.cell_parameters;
    let lastCellUnits = fallback.cell_units;
    let lastCellAlat: number | null = null;
    let lastAtoms = fallback.atoms;
    let lastPositionUnits = fallback.position_units;

    for (let i = 0; i < scope.length; i += 1) {
      const line = scope[i];
      if (line.includes("CELL_PARAMETERS")) {
        const parsed = parseCellBlock(scope, i);
        if (parsed.cellParameters) lastCell = parsed.cellParameters;
        if (parsed.cellUnits) {
          lastCellUnits = parsed.cellUnits;
          lastCellAlat = parsed.cellUnits === "alat" ? parsed.cellAlat : null;
        }
        if (parsed.cellAlat !== null) {
          lastCellAlat = parsed.cellAlat;
        }
      }
      if (line.includes("ATOMIC_POSITIONS")) {
        const parsed = parseAtomsBlock(scope, i);
        if (parsed.atoms.length > 0) {
          lastAtoms = parsed.atoms;
          lastPositionUnits = parsed.positionUnits;
        }
      }
    }

    const normalizedCell = convertCellToAngstrom(lastCell, lastCellUnits, lastCellAlat);

    return {
      position_units: lastPositionUnits,
      cell_units: normalizedCell.cellUnits,
      cell_parameters: normalizedCell.cellParameters,
      atoms: lastAtoms,
    };
  }

  function summarizeCell(structure: SavedStructureData): SavedCellSummary | null {
    if (!structure.cell_parameters) return null;
    const [aVec, bVec, cVec] = structure.cell_parameters;
    const norm = (v: [number, number, number]) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    const dot = (u: [number, number, number], v: [number, number, number]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const cross = (u: [number, number, number], v: [number, number, number]): [number, number, number] => [
      u[1] * v[2] - u[2] * v[1],
      u[2] * v[0] - u[0] * v[2],
      u[0] * v[1] - u[1] * v[0],
    ];
    const angle = (u: [number, number, number], v: [number, number, number]) => {
      const cosine = dot(u, v) / (norm(u) * norm(v));
      const clamped = Math.max(-1, Math.min(1, cosine));
      return (Math.acos(clamped) * 180) / Math.PI;
    };

    const a = norm(aVec);
    const b = norm(bVec);
    const c = norm(cVec);
    const alpha = angle(bVec, cVec);
    const beta = angle(aVec, cVec);
    const gamma = angle(aVec, bVec);
    const volume = Math.abs(dot(aVec, cross(bVec, cVec)));

    return {
      a,
      b,
      c,
      alpha,
      beta,
      gamma,
      volume,
      units: structure.cell_units || "angstrom",
    };
  }

  function buildCalculationData(
    calcResult: any,
    startedAt: string,
    completedAt: string,
    inputContent: string,
    sourceStructure: SavedStructureData,
    sourceDescriptor: { type: "cif" | "optimization"; calc_id?: string },
    magnetismViewerStructure?: SavedStructureData | null,
    hpcMeta?: {
      backend?: string | null;
      hpc_resource_type?: "cpu" | "gpu" | null;
      remote_job_id?: string | null;
      scheduler_state?: string | null;
      remote_node?: string | null;
      remote_workdir?: string | null;
      remote_project_path?: string | null;
      remote_storage_bytes?: number | null;
    },
  ) {
    const isOptimization = config.calculation === "relax" || config.calculation === "vcrelax";
    const usesOptimizedSource = sourceDescriptor.type === "optimization";
    const isPhononReady = !isOptimization && isPhononReadyScf({
      conv_thr: config.conv_thr,
      structure_source: sourceDescriptor,
    });
    const optimizedStructure = isOptimization
      ? extractOptimizedStructure(calcResult?.raw_output || "", sourceStructure)
      : null;
    const optimizedCellSummary = optimizedStructure ? summarizeCell(optimizedStructure) : null;
    const isRemoteRun = isHpcMode || hpcMeta?.backend === "hpc";
    const hubbardManifolds = config.lda_plus_u
      ? Object.fromEntries(getUniqueElements().map((element) => [element, getHubbardManifoldForElement(element)]))
      : null;
    const hubbardSettings = config.lda_plus_u ? buildHubbardSettings(getUniqueElements()) : null;

    return {
      calc_type: isOptimization ? "optimization" : "scf",
      parameters: {
        prefix: "qcortado_scf",
        calculation_mode: config.calculation,
        occupations: config.occupations,
        smearing: config.occupations === "smearing" ? config.smearing : null,
        degauss: config.occupations === "smearing" ? config.degauss : null,
        ecutwfc: config.ecutwfc,
        ecutrho: config.ecutrho,
        kgrid: config.kgrid,
        kgrid_offset: config.kgrid_offset,
        conv_thr: config.conv_thr,
        mixing_beta: config.mixing_beta,
        selected_pseudos: selectedPseudos,
        selected_pseudo_metadata: selectedPseudoMetadata,
        structure_source: sourceDescriptor,
        source_structure: sourceStructure,
        magnetism_viewer_structure: magnetismViewerStructure ?? sourceStructure,
        cell_representation:
          sourceDescriptor.type === "optimization"
            ? "optimized_source"
            : (symmetryTransform ? "primitive_spglib" : "conventional_input"),
        symmetry_spacegroup: symmetryTransform?.spacegroupNumber ?? null,
        symmetry_hall_number: symmetryTransform?.hallNumber ?? null,
        symmetry_error: symmetryError,
        // Feature flags for tags
        nspin: config.nspin,
        noncolin: config.noncolin,
        lspinorb: config.lspinorb,
        starting_magnetization: config.nspin === 1
          ? {}
          : Object.fromEntries(
            getUniqueElements().map((element) => [
              element,
              config.starting_magnetization[element] ?? 0,
            ]),
          ),
        starting_magnetization_theta: config.nspin === 4
          ? Object.fromEntries(
            getUniqueElements().map((element) => [
              element,
              config.starting_magnetization_theta[element] ?? 0,
            ]),
          )
          : {},
        starting_magnetization_phi: config.nspin === 4
          ? Object.fromEntries(
            getUniqueElements().map((element) => [
              element,
              config.starting_magnetization_phi[element] ?? 0,
            ]),
          )
          : {},
        tot_magnetization: config.nspin === 2 ? config.tot_magnetization : null,
        constrained_magnetization: config.constrained_magnetization,
        lda_plus_u: config.lda_plus_u,
        hubbard_projector: config.lda_plus_u ? config.hubbard_projector : null,
        hubbard_formulation: config.lda_plus_u ? config.lda_plus_u_kind : null,
        hubbard_manifold: hubbardManifolds,
        hubbard_parameters: hubbardSettings?.parameters ?? null,
        hubbard_u: config.lda_plus_u ? config.hubbard_u : null,
        hubbard_j: config.lda_plus_u && config.lda_plus_u_kind > 0 ? config.hubbard_j : null,
        vdw_corr: config.vdw_corr,
        // Relaxation parameters
        forc_conv_thr: config.forc_conv_thr,
        etot_conv_thr: config.etot_conv_thr,
        press: config.press,
        optimization_mode: isOptimization ? config.calculation : null,
        optimized_structure: optimizedStructure,
        optimized_cell_summary: optimizedCellSummary,
        ...(isRemoteRun
          ? {
            execution_backend: "hpc",
            hpc_profile_id: activeHpcProfile?.id ?? null,
            hpc_resource_type: hpcMeta?.hpc_resource_type ?? hpcResources.resource_type,
            remote_job_id: hpcMeta?.remote_job_id ?? null,
            scheduler_state: hpcMeta?.scheduler_state ?? null,
            remote_node: hpcMeta?.remote_node ?? null,
            remote_workdir: hpcMeta?.remote_workdir ?? null,
            remote_project_path: hpcMeta?.remote_project_path ?? null,
            remote_storage_bytes: hpcMeta?.remote_storage_bytes ?? null,
          }
          : {}),
      },
      result: calcResult,
      started_at: startedAt,
      completed_at: completedAt,
      input_content: inputContent,
      output_content: calcResult?.raw_output || "",
      tags: [
        ...(isOptimization || usesOptimizedSource ? ["geometry"] : []),
        // Tag as structure-optimized for vcrelax or relax
        ...(config.calculation === "vcrelax" || config.calculation === "relax" ? ["structure-optimized"] : []),
        // Tag as phonon-ready only for optimized-source SCFs with tight convergence.
        ...(isPhononReady ? ["phonon-ready"] : []),
      ],
    };
  }

  const currentSourceStructure: SavedStructureData | null = crystalData
    ? (
        selectedOptimizedStructure
          ? {
              ...selectedOptimizedStructure.structure,
              cell_parameters: selectedOptimizedStructure.structure.cell_parameters || buildCifConventionalStructure(crystalData).cell_parameters,
              cell_units: selectedOptimizedStructure.structure.cell_units || buildCifConventionalStructure(crystalData).cell_units,
            }
          : buildCifConventionalStructure(crystalData)
      )
    : null;
  const currentSourceDescriptor = runSourceDescriptor || (
    selectedOptimizedStructure
      ? { type: "optimization" as const, calc_id: selectedOptimizedStructure.calcId }
      : { type: "cif" as const }
  );

  const calculationData = result && (runSourceStructure || currentSourceStructure)
    ? buildCalculationData(
        result,
        calcStartTime,
        calcEndTime,
        generatedInput,
        runSourceStructure || currentSourceStructure!,
        currentSourceDescriptor,
        runMagnetismViewerStructure || runSourceStructure || currentSourceStructure,
        activeTask?.hpc,
      )
    : null;

  return (
    <div className={`wizard-container wizard-step-${step}`}>
      <div className="wizard-header">
        <button className="back-btn" onClick={onBack}>
          ← Back
        </button>
        <h2>{wizardTitle}</h2>
        <div className="step-indicator">
          <span className={step === "import" ? "active" : "done"}>
            1. Import
          </span>
          <span className={step === "configure" ? "active" : step === "run" || step === "results" ? "done" : ""}>
            2. Configure
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "done" : ""}>
            3. Run
          </span>
          <span className={step === "results" ? "active" : ""}>4. Results</span>
        </div>
      </div>

      {error && <div className="error-banner">{error}</div>}
      {copiedRunSettingsMessage && <div className="info-banner">{copiedRunSettingsMessage}</div>}

      <div className="wizard-content">
        {step === "import" && (
          <div className="import-step">
            <div className="import-zone" onClick={handleImportCIF}>
              <div className="import-icon">📁</div>
              <h3>Import CIF File</h3>
              <p>Click to select a Crystallographic Information File (.cif)</p>
            </div>
          </div>
        )}

        {step === "configure" && crystalData && (
          <div className="configure-step">
            <div className="config-layout">
              <div className="config-left">
                {/* Crystal Info */}
                <section className="config-section">
                  <h3>Crystal Structure</h3>
                  {hasPrimitiveDisplay && (
                    <div className="preset-buttons">
                      <InfoTooltip text="Show conventional CIF lattice parameters">
                        <button
                          type="button"
                          className={`preset-btn ${cellViewMode === "conventional" ? "active" : ""}`}
                          onClick={() => setCellViewMode("conventional")}
                          aria-label="Show conventional CIF lattice parameters"
                        >
                          Conventional
                        </button>
                      </InfoTooltip>
                      <InfoTooltip text="Show primitive lattice parameters used by QE (when detected)">
                        <button
                          type="button"
                          className={`preset-btn ${cellViewMode === "primitive" ? "active" : ""}`}
                          onClick={() => setCellViewMode("primitive")}
                          aria-label="Show primitive lattice parameters used by QE (when detected)"
                        >
                          Primitive
                        </button>
                      </InfoTooltip>
                    </div>
                  )}
                  <div className="info-grid">
                    <div className="info-item">
                      <label>Formula</label>
                      <span>{crystalData.formula_sum || crystalData.formula_structural || "N/A"}</span>
                    </div>
                    <div className="info-item">
                      <label>Space Group</label>
                      <span>{crystalData.space_group_HM || "N/A"}</span>
                    </div>
                    <div className="info-item">
                      <label>Lattice (Å)</label>
                      <span>
                        a={(displayedCellMetrics?.a ?? crystalData.cell_length_a.value).toFixed(3)}, b=
                        {(displayedCellMetrics?.b ?? crystalData.cell_length_b.value).toFixed(3)}, c=
                        {(displayedCellMetrics?.c ?? crystalData.cell_length_c.value).toFixed(3)}
                      </span>
                    </div>
                    <div className="info-item">
                      <label>Angles (°)</label>
                      <span>
                        α={(displayedCellMetrics?.alpha ?? crystalData.cell_angle_alpha.value).toFixed(1)}, β=
                        {(displayedCellMetrics?.beta ?? crystalData.cell_angle_beta.value).toFixed(1)}, γ=
                        {(displayedCellMetrics?.gamma ?? crystalData.cell_angle_gamma.value).toFixed(1)}
                      </span>
                    </div>
                    {hasPrimitiveDisplay && (
                      <div className="info-item">
                        <label>Cell View</label>
                        <span>{cellViewMode === "primitive" ? "Primitive (QE)" : "Conventional (CIF)"}</span>
                      </div>
                    )}
                    <div className="info-item">
                      <label>Atoms</label>
                      <span>{crystalData.atom_sites.length} sites</span>
                    </div>
                  </div>
                </section>

                {/* Pseudopotentials */}
                <section className="config-section">
                  <h3>
                    Pseudopotentials
                    <InfoTooltip text="Pseudopotentials approximate the effect of core electrons (those tightly bound to the nucleus) so the calculation only needs to explicitly handle valence electrons. This dramatically reduces computation cost while maintaining accuracy for chemical bonding and material properties. Each element needs its own pseudopotential file, typically generated with a specific exchange-correlation functional (like PBE)." />
                  </h3>
                  {socPresetSelected && (
                    <p className="pseudo-hint">
                      SOC-compatible mode filters this list to fully relativistic pseudopotentials only. Files without
                      `has_so = true` or `relativistic = "full"` in the UPF header are hidden.
                    </p>
                  )}
                  <p className="pseudo-dir-info">
                    Directory: <code>{pseudoDir}</code>
                    {pseudoFilenames.length > 0 && (
                      <span className="pseudo-count"> ({pseudoFilenames.length} files)</span>
                    )}
                  </p>
                  {isPseudoLoading && (
                    <p className="pseudo-querying">
                      <span className="pseudo-querying-spinner" />
                      Querying pseudopotentials{isHpcMode ? " on cluster" : ""}...
                    </p>
                  )}
                  {pseudoError && (
                    <div className="pseudo-error">{pseudoError}</div>
                  )}
                  <div className="pseudo-list">
                    {getUniqueElements().map((el) => {
                      const matchingPseudoFiles = pseudoFilenames.filter((filename) =>
                        pseudoMatchesElement(filename, el)
                      );
                      const availablePseudoFiles = socPresetSelected
                        ? matchingPseudoFiles.filter((filename) => pseudoMetadataByFilename[filename]?.supports_soc)
                        : matchingPseudoFiles;
                      const selectedValue = selectedPseudos[el] || "";
                      const dropdownOptions = !socPresetSelected
                        && selectedValue
                        && !availablePseudoFiles.includes(selectedValue)
                        ? [selectedValue, ...availablePseudoFiles]
                        : availablePseudoFiles;
                      const warningText = pseudoPresetWarnings[el];
                      return (
                        <div key={el} className="pseudo-row">
                          <label className="pseudo-row-label">
                            <span>{el}</span>
                            {warningText && (
                              <InfoTooltip text={warningText}>
                                <span className="pseudo-warning-icon" aria-label={warningText}>!</span>
                              </InfoTooltip>
                            )}
                          </label>
                          <select
                            value={selectedValue}
                            disabled={isPseudoLoading || (socPresetSelected && !pseudoMetadataLoaded)}
                            onChange={(e) =>
                              {
                                const nextValue = e.target.value;
                                clearPseudoPresetWarning(el);
                                setAutoCutoffFromPseudosEnabled(true);
                                setSelectedPseudos((prev) => ({
                                  ...prev,
                                  [el]: nextValue,
                                }));
                                if (nextValue) {
                                  void ensurePseudoMetadataForFilename(nextValue);
                                }
                              }
                            }
                          >
                            <option value="">
                              {dropdownOptions.length === 0
                                ? socPresetSelected
                                  ? `No SOC-compatible ${el} pseudopotentials found`
                                  : `No ${el} pseudopotentials found`
                                : socPresetSelected && !pseudoMetadataLoaded
                                  ? "Loading SOC-compatible pseudos..."
                                  : "Select..."}
                            </option>
                            {dropdownOptions.map((filename) => (
                              <option key={filename} value={filename}>
                                {filename}
                              </option>
                            ))}
                          </select>
                        </div>
                      );
                    })}
                  </div>
                  {(!socPresetSelected || pseudoMetadataLoaded) && getUniqueElements().some((el) => {
                    const matches = pseudoFilenames.filter((filename) => pseudoMatchesElement(filename, el));
                    return socPresetSelected
                      ? !matches.some((filename) => pseudoMetadataByFilename[filename]?.supports_soc)
                      : matches.length === 0;
                  }) && (
                    <p className="pseudo-hint">
                      Missing pseudopotentials? Download them from the{" "}
                      <a
                        href="https://www.materialscloud.org/discover/sssp/table/precision"
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                        SSSP Precision Library
                      </a>{" "}
                      and place them in your QE pseudo directory. SOC-compatible runs require fully relativistic files.
                    </p>
                  )}
                </section>

                {/* Basic Parameters */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("basic")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.basic ? "expanded" : ""}`}>▶</span>
                    Basic Parameters
                    <InfoTooltip text="Essential parameters for any SCF calculation: energy cutoffs and k-point sampling." />
                  </h3>
                  {expandedSections.basic && (
                    <div className="param-grid">
                      {showPresetRow && (
                        <div className="param-row full-width">
                          <label>Calculation Presets</label>
                          <div className="preset-buttons">
                            {showRelaxPreset && (
                              <InfoTooltip text="Variable-cell relaxation for structure optimization">
                                <button
                                  type="button"
                                  className={`preset-btn preset-relax ${selectedPreset === "relax" ? "active" : ""}`}
                                  onClick={() => applyPreset("relax")}
                                  aria-label="Variable-cell relaxation for structure optimization"
                                >
                                  Optimize Structure
                                </button>
                              </InfoTooltip>
                            )}
                            {showSocPreset && (
                              <InfoTooltip text="Enable the current SOC-compatible setup: noncollinear spin, `nspin = 4`, and spin-orbit coupling.">
                                <button
                                  type="button"
                                  className={`preset-btn ${selectedPreset === "soc" ? "active" : ""}`}
                                  onClick={() => applyPreset("soc")}
                                  aria-label="SOC-compatible SCF calculation"
                                >
                                  SOC-Compatible
                                </button>
                              </InfoTooltip>
                            )}
                          </div>
                        </div>
                      )}

                      <div className="param-row full-width">
                        <label>
                          Pseudopotential Presets
                          <InfoTooltip text="Auto-select pseudopotentials by library or pseudopotential type. SSSP is the default. PAW, USPP, and NCPP try to match that type for each element and fall back to another pseudo with a warning when needed." />
                        </label>
                        <div className="preset-buttons">
                          <InfoTooltip text="Prefer the SSSP-recommended pseudopotentials for each element.">
                            <button
                              type="button"
                              className={`preset-btn ${selectedPseudoPreset === "sssp" ? "active" : ""}`}
                              onClick={() => handlePseudopotentialPresetClick("sssp")}
                              aria-label="SSSP pseudopotential preset"
                            >
                              SSSP
                            </button>
                          </InfoTooltip>
                          <InfoTooltip text="Prefer PAW pseudopotentials when available.">
                            <button
                              type="button"
                              className={`preset-btn ${selectedPseudoPreset === "paw" ? "active" : ""}`}
                              onClick={() => handlePseudopotentialPresetClick("paw")}
                              aria-label="PAW pseudopotential preset"
                            >
                              PAW
                            </button>
                          </InfoTooltip>
                          <InfoTooltip text="Prefer ultrasoft pseudopotentials when available.">
                            <button
                              type="button"
                              className={`preset-btn ${selectedPseudoPreset === "uspp" ? "active" : ""}`}
                              onClick={() => handlePseudopotentialPresetClick("uspp")}
                              aria-label="USPP pseudopotential preset"
                            >
                              USPP
                            </button>
                          </InfoTooltip>
                          <InfoTooltip text="Prefer norm-conserving pseudopotentials when available.">
                            <button
                              type="button"
                              className={`preset-btn ${selectedPseudoPreset === "ncpp" ? "active" : ""}`}
                              onClick={() => handlePseudopotentialPresetClick("ncpp")}
                              aria-label="NCPP pseudopotential preset"
                            >
                              NCPP
                            </button>
                          </InfoTooltip>
                        </div>
                      </div>

                      <div className="param-row">
                        <label>
                          Structure Source
                          <InfoTooltip text={structureSourceTooltip} />
                        </label>
                        <select
                          value={structureSource}
                          onChange={(e) => setStructureSource(e.target.value)}
                        >
                          <option value="cif">From CIF (original structure)</option>
                          {optimizedStructures.map((option) => (
                            <option key={option.calcId} value={option.calcId}>
                              {option.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      {isOptimizationWizard && (
                        <div className="param-row">
                          <label>
                            Calculation Type
                            <InfoTooltip text="scf: ground state energy only. relax: optimize atomic positions. vcrelax: optimize both positions and cell shape/size." />
                          </label>
                          <select
                            value={config.calculation}
                            onChange={(e) => setConfig((prev) => ({ ...prev, calculation: e.target.value as "scf" | "relax" | "vcrelax" }))}
                            disabled={Boolean(lockedPreset)}
                          >
                            <option value="scf">SCF (ground state)</option>
                            <option value="relax">Relax (fixed cell)</option>
                            <option value="vcrelax">VC-Relax (variable cell)</option>
                          </select>
                        </div>
                      )}

                      {/* Relaxation parameters - shown only for relax/vcrelax */}
                      {(config.calculation === "relax" || config.calculation === "vcrelax") && (
                        <>
                          <div className="param-row">
                            <label>
                              Force Convergence
                              <InfoTooltip text="Convergence threshold for forces (Ry/Bohr). Relaxation stops when all forces are below this value. Use 1e-4 for rough optimization, 1e-5 for phonon calculations." />
                            </label>
                            <div className="param-input-group">
                              <input
                                type="text"
                                value={config.forc_conv_thr}
                                onChange={(e) => {
                                  const val = parseFloat(e.target.value);
                                  if (!isNaN(val)) {
                                    setConfig((prev) => ({ ...prev, forc_conv_thr: val }));
                                  }
                                }}
                              />
                              <span className="param-unit">Ry/Bohr</span>
                            </div>
                          </div>
                          <div className="param-row">
                            <label>
                              Energy Convergence
                              <InfoTooltip text="Convergence threshold for total energy change between ionic steps (Ry)." />
                            </label>
                            <div className="param-input-group">
                              <input
                                type="text"
                                value={config.etot_conv_thr}
                                onChange={(e) => {
                                  const val = parseFloat(e.target.value);
                                  if (!isNaN(val)) {
                                    setConfig((prev) => ({ ...prev, etot_conv_thr: val }));
                                  }
                                }}
                              />
                              <span className="param-unit">Ry</span>
                            </div>
                          </div>
                        </>
                      )}

                      {/* VC-Relax specific: target pressure */}
                      {config.calculation === "vcrelax" && (
                        <div className="param-row">
                          <label>
                            Target Pressure
                            <InfoTooltip text="Target external pressure in kbar. Use 0 for ambient pressure, positive for compression." />
                          </label>
                          <div className="param-input-group">
                            <input
                              type="number"
                              value={config.press}
                              onChange={(e) => setConfig((prev) => ({ ...prev, press: parseFloat(e.target.value) || 0 }))}
                            />
                            <span className="param-unit">kbar</span>
                          </div>
                        </div>
                      )}

                      <div className="param-row">
                        <label>
                          Wavefunction Cutoff
                          <InfoTooltip text="Energy cutoff for plane-wave expansion of wavefunctions (in Rydberg). Higher = more accurate but slower. Typical: 30-80 Ry." />
                        </label>
                        <div className="param-input-group">
                          <input
                            type="number"
                            className={selectedPseudoCutoffSummary.wfcStatus === "idle" ? "" : `cutoff-input cutoff-input-${selectedPseudoCutoffSummary.wfcStatus}`}
                            value={config.ecutwfc}
                            onChange={(e) => setConfig((prev) => ({ ...prev, ecutwfc: parseFloat(e.target.value) }))} />
                          <span className="param-unit">Ry</span>
                          {selectedPseudoCutoffSummary.wfcStatus !== "idle" && (
                            <span className={`cutoff-status cutoff-status-${selectedPseudoCutoffSummary.wfcStatus}`}>
                              {cutoffStatusLabel(
                                selectedPseudoCutoffSummary.wfcStatus,
                                selectedPseudoCutoffSummary.wfcProvenance,
                                selectedPseudoCutoffSummary.wfcDerivation,
                                selectedPseudoCutoffSummary.wfcRatio,
                              )}
                            </span>
                          )}
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          Charge Density Cutoff
                          <InfoTooltip text="Energy cutoff for charge density (in Rydberg). Typical fallback is 4x ecutwfc for NC/PAW and 8x for US pseudopotentials." />
                        </label>
                        <div className="param-input-group">
                          <input
                            type="number"
                            className={selectedPseudoCutoffSummary.rhoStatus === "idle" ? "" : `cutoff-input cutoff-input-${selectedPseudoCutoffSummary.rhoStatus}`}
                            value={config.ecutrho}
                            onChange={(e) => setConfig((prev) => ({ ...prev, ecutrho: parseFloat(e.target.value) }))} />
                          <span className="param-unit">Ry</span>
                          {selectedPseudoCutoffSummary.rhoStatus !== "idle" && (
                            <span className={`cutoff-status cutoff-status-${selectedPseudoCutoffSummary.rhoStatus}`}>
                              {cutoffStatusLabel(
                                selectedPseudoCutoffSummary.rhoStatus,
                                selectedPseudoCutoffSummary.rhoProvenance,
                                selectedPseudoCutoffSummary.rhoDerivation,
                                selectedPseudoCutoffSummary.rhoRatio,
                              )}
                            </span>
                          )}
                        </div>
                      </div>
                      {selectedPseudoCutoffSummary.hasInferredCutoff && (
                        <div className="cutoff-metadata-warning">
                          QCortado autocalculated missing cutoff values from the companion cutoff when needed using 4x for NC/PAW and 8x for US pseudopotentials. Verify these values with a convergence test.
                        </div>
                      )}
                      {selectedPseudoCutoffSummary.hasMissingCutoff && (
                        <div className="cutoff-metadata-warning cutoff-metadata-warning-missing">
                          Some selected pseudopotentials do not embed enough cutoff metadata for QCortado to determine every recommended value.
                        </div>
                      )}
                      <div className="param-row">
                        <label>
                          K-point Grid
                          <InfoTooltip text="Monkhorst-Pack grid for Brillouin zone sampling. Denser = more accurate. Metals need denser grids." />
                        </label>
                        <div className="kgrid-inputs">
                          <input type="number" value={config.kgrid[0]}
                            onChange={(e) => setConfig((prev) => ({ ...prev, kgrid: [parseInt(e.target.value), prev.kgrid[1], prev.kgrid[2]] }))} />
                          <span>×</span>
                          <input type="number" value={config.kgrid[1]}
                            onChange={(e) => setConfig((prev) => ({ ...prev, kgrid: [prev.kgrid[0], parseInt(e.target.value), prev.kgrid[2]] }))} />
                          <span>×</span>
                          <input type="number" value={config.kgrid[2]}
                            onChange={(e) => setConfig((prev) => ({ ...prev, kgrid: [prev.kgrid[0], prev.kgrid[1], parseInt(e.target.value)] }))} />
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          K-point Offset
                          <InfoTooltip text="Shift the k-point grid. Use (1,1,1) for metals to avoid high-symmetry points, (0,0,0) for insulators." />
                        </label>
                        <div className="kgrid-inputs">
                          <input type="number" value={config.kgrid_offset[0]} min={0} max={1}
                            onChange={(e) => setConfig((prev) => ({ ...prev, kgrid_offset: [parseInt(e.target.value), prev.kgrid_offset[1], prev.kgrid_offset[2]] as [number, number, number] }))} />
                          <span>,</span>
                          <input type="number" value={config.kgrid_offset[1]} min={0} max={1}
                            onChange={(e) => setConfig((prev) => ({ ...prev, kgrid_offset: [prev.kgrid_offset[0], parseInt(e.target.value), prev.kgrid_offset[2]] as [number, number, number] }))} />
                          <span>,</span>
                          <input type="number" value={config.kgrid_offset[2]} min={0} max={1}
                            onChange={(e) => setConfig((prev) => ({ ...prev, kgrid_offset: [prev.kgrid_offset[0], prev.kgrid_offset[1], parseInt(e.target.value)] as [number, number, number] }))} />
                        </div>
                      </div>
                    </div>
                  )}
                  {ssspMissing && expandedSections.basic && (
                    <p className="sssp-hint">
                      Cutoff values are defaults unless the pseudopotential header supplies recommendations. For optimized values, download the{" "}
                      <a href="https://www.materialscloud.org/discover/sssp/table/precision" target="_blank" rel="noopener noreferrer">SSSP library</a>.
                    </p>
                  )}
                  {hasPseudoCutoffMetadata && expandedSections.basic && (
                    <p className="sssp-success">Cutoffs auto-configured from pseudopotential metadata when available, with SSSP used as a fallback.</p>
                  )}
                </section>

                {/* Electronic Structure */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("electronic")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.electronic ? "expanded" : ""}`}>▶</span>
                    Electronic Structure
                    <InfoTooltip text="Control how electronic states are occupied and how the Fermi level is determined." />
                  </h3>
                  {expandedSections.electronic && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Occupations
                          <InfoTooltip text="How to fill electronic states. 'smearing' for metals, 'tetrahedra' for accurate DOS, 'fixed' for insulators with gap." />
                        </label>
                        <select value={config.occupations}
                          onChange={(e) => setConfig((prev) => ({ ...prev, occupations: e.target.value as any }))}>
                          <option value="smearing">Smearing (metals/default)</option>
                          <option value="tetrahedra">Tetrahedra</option>
                          <option value="tetrahedra_lin">Tetrahedra (linear)</option>
                          <option value="tetrahedra_opt">Tetrahedra (optimized)</option>
                          <option value="fixed">Fixed</option>
                        </select>
                      </div>
                      {config.occupations === "smearing" && (
                        <>
                          <div className="param-row">
                            <label>
                              Smearing Type
                              <InfoTooltip text="gaussian: simple, good for most cases. methfessel-paxton: better for metals. marzari-vanderbilt: 'cold smearing', good for forces. fermi-dirac: physical but needs higher temp." />
                            </label>
                            <select value={config.smearing}
                              onChange={(e) => setConfig((prev) => ({ ...prev, smearing: e.target.value as any }))}>
                              <option value="gaussian">Gaussian</option>
                              <option value="methfessel-paxton">Methfessel-Paxton</option>
                              <option value="marzari-vanderbilt">Marzari-Vanderbilt (cold)</option>
                              <option value="fermi-dirac">Fermi-Dirac</option>
                            </select>
                          </div>
                          <div className="param-row">
                            <label>
                              Smearing Width (degauss)
                              <InfoTooltip text="Smearing width in Ry. Smaller = more accurate but harder to converge. Typical: 0.005-0.02 Ry for metals." />
                            </label>
                            <div className="param-input-group">
                              <input type="number" step="0.001" value={config.degauss}
                                onChange={(e) => setConfig((prev) => ({ ...prev, degauss: parseFloat(e.target.value) }))} />
                              <span className="param-unit">Ry</span>
                            </div>
                          </div>
                        </>
                      )}
                      <div className="param-row">
                        <label>
                          Number of Bands
                          <InfoTooltip text="Total number of Kohn-Sham states. Leave empty for automatic (occupied + some empty). Increase for accurate DOS or optical properties." />
                        </label>
                        <div className="param-input-group">
                          <input type="number" value={config.nbnd ?? ""} placeholder="auto"
                            onChange={(e) => setConfig((prev) => ({ ...prev, nbnd: e.target.value ? parseInt(e.target.value) : null }))} />
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          Total Charge
                          <InfoTooltip text="Total charge of the system in units of e. Positive = remove electrons, negative = add electrons. 0 = neutral." />
                        </label>
                        <div className="param-input-group">
                          <input type="number" step="0.1" value={config.tot_charge}
                            onChange={(e) => setConfig((prev) => ({ ...prev, tot_charge: parseFloat(e.target.value) || 0 }))} />
                          <span className="param-unit">e</span>
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          XC Functional Override
                          <InfoTooltip text="Override the exchange-correlation functional from pseudopotentials. Leave empty to use pseudopotential default. Examples: PBE, PBEsol, SCAN, HSE, B3LYP." />
                        </label>
                        <input type="text" value={config.input_dft} placeholder="(use pseudopotential default)"
                          onChange={(e) => setConfig((prev) => ({ ...prev, input_dft: e.target.value }))} />
                      </div>
                    </div>
                  )}
                </section>

                {/* Magnetism & Spin */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("magnetism")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.magnetism ? "expanded" : ""}`}>▶</span>
                    Magnetism & Spin
                    <InfoTooltip text="Enable spin-polarized calculations for magnetic materials and spin-orbit coupling for heavy elements." />
                  </h3>
                  {expandedSections.magnetism && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Spin Polarization
                          <InfoTooltip text="nspin=1: non-magnetic. nspin=2: collinear magnetic (spin up/down). nspin=4: non-collinear (required for spin-orbit)." />
                        </label>
                        <select value={config.nspin}
                          onChange={(e) => {
                            const nspin = parseInt(e.target.value) as 1 | 2 | 4;
                            setConfig((prev) => ({
                              ...prev,
                              nspin,
                              noncolin: nspin === 4,
                              lspinorb: nspin === 4 ? prev.lspinorb : false,
                            }));
                          }}>
                          <option value={1}>Non-polarized (nspin=1)</option>
                          <option value={2}>Collinear (nspin=2)</option>
                          <option value={4}>Non-collinear (nspin=4)</option>
                        </select>
                      </div>
                      {config.nspin === 4 && (
                        <div className="param-row">
                          <label>
                            Spin-Orbit Coupling
                            <InfoTooltip text="Include spin-orbit interaction. Required for topological properties and heavy elements. Needs fully-relativistic pseudopotentials." />
                          </label>
                          <label className="toggle-label">
                            <input type="checkbox" checked={config.lspinorb}
                              onChange={(e) => setConfig((prev) => ({ ...prev, lspinorb: e.target.checked }))} />
                            <span>Enable spin-orbit coupling</span>
                          </label>
                        </div>
                      )}
                      {config.nspin === 2 && (
                        <div className="param-row">
                          <label>
                            Total Magnetization
                            <InfoTooltip text="Fix total magnetization (N_up - N_down). Leave empty for unconstrained." />
                          </label>
                          <div className="param-input-group">
                            <input type="number" step="0.1" value={config.tot_magnetization ?? ""} placeholder="auto"
                              onChange={(e) => setConfig((prev) => ({ ...prev, tot_magnetization: e.target.value ? parseFloat(e.target.value) : null }))} />
                            <span className="param-unit">Bohr mag</span>
                          </div>
                        </div>
                      )}
                      {(config.nspin === 2 || config.nspin === 4) && (
                        <>
                          <div className="param-row full-width">
                            <label>
                              Starting Magnetization (per element)
                              <InfoTooltip text="Initial spin polarization for each element (-1 to +1). Helps SCF converge to correct magnetic state." />
                            </label>
                          </div>
                          <div className="magnetization-grid">
                            {getUniqueElements().map((el) => (
                              <div key={el} className="mag-row">
                                <label>{el}</label>
                                <input type="number" step="0.1" min={-1} max={1}
                                  value={config.starting_magnetization[el] ?? 0}
                                  onChange={(e) => setConfig((prev) => ({
                                    ...prev,
                                    starting_magnetization: { ...prev.starting_magnetization, [el]: parseFloat(e.target.value) || 0 }
                                  }))} />
                              </div>
                            ))}
                          </div>
                          {config.nspin === 4 && (
                            <>
                              <div className="param-row full-width">
                                <label>
                                  Starting Magnetization Direction (per element)
                                  <InfoTooltip text="Noncollinear magnetization direction in degrees: theta is the polar angle from z, phi is the azimuthal angle in the xy plane." />
                                </label>
                              </div>
                              <div className="magnetization-grid magnetization-angle-grid">
                                {getUniqueElements().map((el) => (
                                  <div key={el} className="mag-row mag-angle-row">
                                    <label>{el}</label>
                                    <span>theta</span>
                                    <input type="number" step="1" min={0} max={180}
                                      value={config.starting_magnetization_theta[el] ?? 0}
                                      onChange={(e) => setConfig((prev) => ({
                                        ...prev,
                                        starting_magnetization_theta: {
                                          ...prev.starting_magnetization_theta,
                                          [el]: parseFloat(e.target.value) || 0,
                                        },
                                      }))} />
                                    <span>phi</span>
                                    <input type="number" step="1" min={0} max={360}
                                      value={config.starting_magnetization_phi[el] ?? 0}
                                      onChange={(e) => setConfig((prev) => ({
                                        ...prev,
                                        starting_magnetization_phi: {
                                          ...prev.starting_magnetization_phi,
                                          [el]: parseFloat(e.target.value) || 0,
                                        },
                                      }))} />
                                  </div>
                                ))}
                              </div>
                            </>
                          )}
                        </>
                      )}
                    </div>
                  )}
                </section>

                {/* SCF Convergence */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("convergence")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.convergence ? "expanded" : ""}`}>▶</span>
                    SCF Convergence
                    <InfoTooltip text="Parameters controlling how the self-consistent field iteration converges." />
                  </h3>
                  {expandedSections.convergence && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Convergence Threshold
                          <InfoTooltip text="SCF stops when energy change falls below this value (Ry). 1e-12 is the default for the built-in presets, and 1e-12+ remains appropriate for phonons." />
                        </label>
                        <div className="param-input-group">
                          <input type="text" value={convThrInput}
                            onChange={(e) => setConvThrInput(e.target.value)}
                            onBlur={() => {
                              const parsed = parseFloat(convThrInput);
                              if (!isNaN(parsed) && parsed > 0) {
                                setConfig((prev) => ({ ...prev, conv_thr: parsed }));
                              } else {
                                setConvThrInput(config.conv_thr.toString());
                              }
                            }}
                            onKeyDown={(e) => { if (e.key === "Enter") e.currentTarget.blur(); }} />
                          <span className="param-unit">Ry</span>
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          Max SCF Iterations
                          <InfoTooltip text="Maximum number of SCF iterations before giving up." />
                        </label>
                        <input type="number" value={config.electron_maxstep}
                          onChange={(e) => setConfig((prev) => ({ ...prev, electron_maxstep: parseInt(e.target.value) }))} />
                      </div>
                      <div className="param-row">
                        <label>
                          Mixing Mode
                          <InfoTooltip text="How to mix old and new charge density. 'plain' is default, 'TF' or 'local-TF' can help metals converge." />
                        </label>
                        <select value={config.mixing_mode}
                          onChange={(e) => setConfig((prev) => ({ ...prev, mixing_mode: e.target.value as any }))}>
                          <option value="plain">Plain</option>
                          <option value="TF">Thomas-Fermi</option>
                          <option value="local-TF">Local Thomas-Fermi</option>
                        </select>
                      </div>
                      <div className="param-row">
                        <label>
                          Mixing Beta
                          <InfoTooltip text="Fraction of new density to mix in (0-1). Lower = more stable but slower. 0.7 default, use 0.1-0.3 for difficult cases." />
                        </label>
                        <input type="number" step="0.05" min={0} max={1} value={config.mixing_beta}
                          onChange={(e) => setConfig((prev) => ({ ...prev, mixing_beta: parseFloat(e.target.value) }))} />
                      </div>
                      <div className="param-row">
                        <label>
                          Mixing Dimensions
                          <InfoTooltip text="Number of previous iterations used in Broyden mixing. Higher can help difficult cases converge." />
                        </label>
                        <input type="number" min={2} max={20} value={config.mixing_ndim}
                          onChange={(e) => setConfig((prev) => ({ ...prev, mixing_ndim: parseInt(e.target.value) }))} />
                      </div>
                      <div className="param-row">
                        <label>
                          Diagonalization
                          <InfoTooltip text="Algorithm for solving eigenvalue problem. 'david' is fast for most cases, 'cg' for difficult cases, 'ppcg'/'paro' for large parallel runs." />
                        </label>
                        <select value={config.diagonalization}
                          onChange={(e) => setConfig((prev) => ({ ...prev, diagonalization: e.target.value as any }))}>
                          <option value="david">Davidson</option>
                          <option value="cg">Conjugate Gradient</option>
                          <option value="ppcg">PPCG</option>
                          <option value="paro">ParO</option>
                          <option value="rmm-davidson">RMM-Davidson</option>
                        </select>
                      </div>
                      <div className="param-row">
                        <label>
                          Starting Potential
                          <InfoTooltip text="'atomic': superposition of atomic potentials (default). 'file': read from previous calculation." />
                        </label>
                        <select value={config.startingpot}
                          onChange={(e) => setConfig((prev) => ({ ...prev, startingpot: e.target.value as any }))}>
                          <option value="atomic">Atomic</option>
                          <option value="file">From file</option>
                        </select>
                      </div>
                      <div className="param-row">
                        <label>
                          Starting Wavefunctions
                          <InfoTooltip text="Initial guess for wavefunctions. 'atomic': atomic orbitals. 'atomic+random': with randomization. 'random': fully random." />
                        </label>
                        <select value={config.startingwfc}
                          onChange={(e) => setConfig((prev) => ({ ...prev, startingwfc: e.target.value as any }))}>
                          <option value="atomic">Atomic</option>
                          <option value="atomic+random">Atomic + Random</option>
                          <option value="random">Random</option>
                          <option value="file">From file</option>
                        </select>
                      </div>
                    </div>
                  )}
                </section>

                {/* DFT+U */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("dftu")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.dftu ? "expanded" : ""}`}>▶</span>
                    DFT+U (Hubbard Correction)
                    <InfoTooltip text="Add Hubbard U correction to improve description of localized d and f electrons in transition metals and lanthanides." />
                    {isHubbardRecommended && (
                      <span className="section-recommendation">
                        Recommended for this system
                        <InfoTooltip text={`QCortado detected localized Hubbard manifolds (${hubbardRecommendationText}). ${hubbardRecommendations.map((entry) => entry.reason).join(" ")}`} />
                      </span>
                    )}
                  </h3>
                  {expandedSections.dftu && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Enable DFT+U
                          <InfoTooltip text="Apply Hubbard U correction. Essential for correlated materials like transition metal oxides." />
                        </label>
                        <label className="toggle-label">
                          <input type="checkbox" checked={config.lda_plus_u}
                            onChange={(e) => setConfig((prev) => ({ ...prev, lda_plus_u: e.target.checked }))} />
                          <span>Enable Hubbard correction</span>
                        </label>
                        {isHubbardRecommended && (
                          <p className="field-hint">
                            Recommended for {hubbardRecommendationText}. U keeps any saved run value, otherwise uses a general 6.0 eV guess, and the button lets you pull in a saved LRT value on demand.
                          </p>
                        )}
                      </div>
                      {config.lda_plus_u && (
                        <>
                          <div className="param-row">
                            <label>
                              Hubbard Projectors
                              <InfoTooltip text="QE 7.3.1+ requires the projector type on the HUBBARD card. Ortho-atomic is recommended when supported by the pseudopotentials." />
                            </label>
                            <select value={config.hubbard_projector}
                              onChange={(e) => setConfig((prev) => ({ ...prev, hubbard_projector: e.target.value as HubbardProjector }))}>
                              <option value="ortho-atomic">ortho-atomic</option>
                              <option value="atomic">atomic</option>
                              <option value="norm-atomic">norm-atomic</option>
                              <option value="wf">wf</option>
                              <option value="pseudo">pseudo</option>
                            </select>
                          </div>
                          <div className="param-row">
                            <label>
                              Hubbard Formulation
                              <InfoTooltip text="Dudarev DFT+U writes U lines. Liechtenstein DFT+U+J writes U and J lines. Dudarev DFT+U+J0 writes U and J0 lines." />
                            </label>
                            <select value={config.lda_plus_u_kind}
                              onChange={(e) => setConfig((prev) => ({ ...prev, lda_plus_u_kind: parseInt(e.target.value) as 0 | 1 | 2 }))}>
                              <option value={0}>DFT+U (Dudarev)</option>
                              <option value={1}>DFT+U+J (Liechtenstein)</option>
                              <option value={2}>DFT+U+J0 (Dudarev)</option>
                            </select>
                          </div>
                          <div className="param-row full-width">
                            <label>Hubbard manifolds and values (per element, in eV)</label>
                          </div>
                          <div className="hubbard-grid">
                            {getUniqueElements().map((el) => {
                              const manifold = getHubbardManifoldForElement(el);
                              const lrtValue = manifold
                                ? getLatestHubbardLrtValue(el, manifold, (initialCif?.calculations ?? []) as any[])
                                : null;
                              return (
                              <div key={el} className="hubbard-row">
                                <label className="hubbard-element">{el}</label>
                                <label className="hubbard-field">
                                  <span>Manifold</span>
                                  <input type="text"
                                    value={config.hubbard_manifold[el] ?? ""}
                                    placeholder={getDefaultHubbardManifold(el) || "manual"}
                                    onChange={(e) => setConfig((prev) => ({
                                      ...prev,
                                      hubbard_manifold: { ...prev.hubbard_manifold, [el]: e.target.value }
                                    }))} />
                                </label>
                                <label className="hubbard-field hubbard-value-field">
                                  <span>U</span>
                                  <span className="hubbard-value-inline">
                                    <input type="number" step="0.1" min={0}
                                      value={config.hubbard_u[el] ?? 0}
                                      onChange={(e) => setConfig((prev) => ({
                                        ...prev,
                                        hubbard_u: { ...prev.hubbard_u, [el]: parseFloat(e.target.value) || 0 }
                                      }))}
                                      onInput={() => setManuallyEditedHubbardU((prev) => ({ ...prev, [el]: true }))} />
                                    {lrtValue && (
                                      <button
                                        type="button"
                                        className="secondary-button hubbard-lrt-use-btn"
                                        onClick={() => applyCalculatedHubbardU(el)}
                                        title="Populate this element's U from the latest completed Hubbard LRT calculation saved in the project."
                                      >
                                        Auto
                                      </button>
                                    )}
                                  </span>
                                </label>
                                {config.lda_plus_u_kind > 0 && (
                                  <label className="hubbard-field hubbard-value-field">
                                    <span>{config.lda_plus_u_kind === 2 ? "J0" : "J"}</span>
                                    <input type="number" step="0.1" min={0}
                                      value={config.hubbard_j[el] ?? 0}
                                      onChange={(e) => setConfig((prev) => ({
                                        ...prev,
                                        hubbard_j: { ...prev.hubbard_j, [el]: parseFloat(e.target.value) || 0 }
                                      }))}
                                      onInput={() => setManuallyEditedHubbardJ((prev) => ({ ...prev, [el]: true }))} />
                                    {config.lda_plus_u_kind === 1 && hubbardJDefaultLabels[el] && (
                                      <small className="hubbard-default-label">{hubbardJDefaultLabels[el]}</small>
                                    )}
                                  </label>
                                )}
                              </div>
                              );
                            })}
                          </div>
                        </>
                      )}
                    </div>
                  )}
                </section>

                {/* Van der Waals */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("vdw")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.vdw ? "expanded" : ""}`}>▶</span>
                    Van der Waals Corrections
                    <InfoTooltip text="Add dispersion corrections for systems where van der Waals interactions are important (layered materials, molecules, etc.)." />
                  </h3>
                  {expandedSections.vdw && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          vdW Correction Method
                          <InfoTooltip text="None: standard DFT. DFT-D2/D3: Grimme's empirical corrections. TS-vdW: Tkatchenko-Scheffler. XDM: exchange-hole dipole moment." />
                        </label>
                        <select value={config.vdw_corr}
                          onChange={(e) => setConfig((prev) => ({ ...prev, vdw_corr: e.target.value as any }))}>
                          <option value="none">None</option>
                          <option value="grimme-d2">Grimme DFT-D2</option>
                          <option value="grimme-d3">Grimme DFT-D3</option>
                          <option value="ts-vdw">Tkatchenko-Scheffler</option>
                          <option value="xdm">XDM</option>
                        </select>
                      </div>
                    </div>
                  )}
                </section>

                {/* Advanced Options */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("advanced")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.advanced ? "expanded" : ""}`}>▶</span>
                    Advanced Options
                    <InfoTooltip text="Additional settings for special cases: isolated systems, output control, etc." />
                  </h3>
                  {expandedSections.advanced && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Isolated System
                          <InfoTooltip text="For molecules or slabs, remove spurious interactions with periodic images. 'martyna-tuckerman' for molecules, 'esm' or '2D' for slabs." />
                        </label>
                        <select value={config.assume_isolated}
                          onChange={(e) => setConfig((prev) => ({ ...prev, assume_isolated: e.target.value as any }))}>
                          <option value="none">None (3D periodic)</option>
                          <option value="makov-payne">Makov-Payne (charged)</option>
                          <option value="martyna-tuckerman">Martyna-Tuckerman (molecules)</option>
                          <option value="esm">ESM (slabs)</option>
                          <option value="2D">2D (slabs)</option>
                        </select>
                      </div>
                      <div className="param-row">
                        <label>
                          Verbosity
                          <InfoTooltip text="Amount of output. 'high' gives more detail for debugging." />
                        </label>
                        <select value={config.verbosity}
                          onChange={(e) => setConfig((prev) => ({ ...prev, verbosity: e.target.value as any }))}>
                          <option value="low">Low</option>
                          <option value="high">High</option>
                        </select>
                      </div>
                      <div className="param-row">
                        <label>
                          Disk I/O
                          <InfoTooltip text="How much to write to disk. 'low' saves space, 'nowf' doesn't save wavefunctions (can't restart)." />
                        </label>
                        <select value={config.disk_io}
                          onChange={(e) => setConfig((prev) => ({ ...prev, disk_io: e.target.value as any }))}>
                          <option value="low">Low</option>
                          <option value="medium">Medium</option>
                          <option value="high">High</option>
                          <option value="nowf">No wavefunctions</option>
                        </select>
                      </div>
                      <div className="param-row">
                        <label>Calculate Forces</label>
                        <label className="toggle-label">
                          <input type="checkbox" checked={config.tprnfor}
                            onChange={(e) => setConfig((prev) => ({ ...prev, tprnfor: e.target.checked }))} />
                          <span>Print forces on atoms</span>
                        </label>
                      </div>
                      <div className="param-row">
                        <label>Calculate Stress</label>
                        <label className="toggle-label">
                          <input type="checkbox" checked={config.tstress}
                            onChange={(e) => setConfig((prev) => ({ ...prev, tstress: e.target.checked }))} />
                          <span>Print stress tensor</span>
                        </label>
                      </div>
                    </div>
                  )}
                </section>

                {isHpcMode ? (
                  <HpcRunSettings
                    profileId={activeHpcProfile?.id ?? null}
                    profileName={activeHpcProfile?.name ?? "Andromeda"}
                    taskKind="scf"
                    commandLines={hpcCommandLines}
                    resources={hpcResources}
                    resourceMode={activeHpcProfile?.resource_mode ?? "both"}
                    defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
                    defaultGpuResources={activeHpcProfile?.default_gpu_resources ?? null}
                    onResourcesChange={setHpcResources}
                    disabled={isRunning}
                  />
                ) : (
                  <section className="config-section">
                    <h3>
                      Parallelization
                      <InfoTooltip text="MPI (Message Passing Interface) allows running calculations in parallel across multiple CPU cores, significantly speeding up large calculations. Your Quantum ESPRESSO installation must be compiled with MPI support for this to work." />
                    </h3>

                    <div className="mpi-toggle-row">
                      <label className="toggle-label">
                        <input
                          type="checkbox"
                          checked={mpiEnabled}
                          onChange={(e) => setMpiEnabled(e.target.checked)}
                          disabled={!mpiAvailable}
                        />
                        <span>Enable MPI parallel execution</span>
                      </label>
                      {!mpiAvailable && (
                        <span className="mpi-unavailable">
                          (mpirun not found on system)
                        </span>
                      )}
                    </div>

                    {mpiEnabled && (
                      <div className="mpi-settings">
                        <div className="param-row">
                          <label>
                            Number of processes
                            <InfoTooltip text="Number of parallel MPI processes to use. More processes can speed up large calculations but have overhead for small systems. Using too many processes for a small system may actually slow things down." />
                          </label>
                          <div className="mpi-input-group">
                            <input
                              type="number"
                              min={1}
                              max={cpuCount}
                              value={mpiProcs}
                              onChange={(e) => setMpiProcs(Math.max(1, Math.min(cpuCount, parseInt(e.target.value) || 1)))}
                            />
                            <span className="mpi-core-info">
                              of {cpuCount} available cores
                            </span>
                          </div>
                        </div>

                        <div className="mpi-warning">
                          <strong>Note:</strong> MPI only works if Quantum ESPRESSO was compiled with MPI support.
                          If the calculation fails, try disabling MPI or check your QE installation.
                        </div>
                      </div>
                    )}
                  </section>
                )}

                <div className="run-btn-group">
                  <button
                    className={`secondary-button paste-copied-run-btn ${copiedRunSettingsAvailable ? "" : "is-disabled"}`}
                    onClick={() => void handlePasteFromCopiedRun()}
                    title={copiedRunSettingsAvailable
                      ? "Paste settings from copied run"
                      : "Copy an Existing Run's Settings First"}
                  >
                    Paste from copied run
                  </button>
                  {projectContext && !isHpcMode && (
                    <button
                      className="secondary-button"
                      onClick={() => void handleQueue()}
                      disabled={!canRun()}
                    >
                      Queue Task
                    </button>
                  )}
                  <button
                    className="run-btn"
                    onClick={handleRun}
                    disabled={!canRun() || hasBlockingExternalTask}
                  >
                    {isHpcMode
                      ? "Submit SCF to Andromeda"
                      : mpiEnabled && mpiProcs > 1
                      ? `Run SCF Calculation (${mpiProcs} cores)`
                      : "Run SCF Calculation"}
                  </button>
                </div>
              </div>

              <div className="config-right">
                <div className="viewer-container">
                  <UnitCellViewer crystalData={crystalData} />
                </div>
              </div>
            </div>
          </div>
        )}

        {(step === "run" || step === "results") && (
          <div className="wizard-step run-step run-step-focused scf-run-step">
            <div className="run-step-headline">
              <h3>
                {isOptimizationWizard
                  ? (isRunning ? "Running Structure Optimization" : "Structure Optimization Output")
                  : (isRunning ? "Running SCF Calculation" : "SCF Output")}
              </h3>
              <span className={`run-step-status-pill ${isRunning ? "running" : error ? "error" : "idle"}`}>
                {isRunning ? "Live output" : error ? "Run failed" : "Output"}
              </span>
            </div>

            <div className="run-status-rail scf-run-status">
              <ProgressBar
                status={progress.status}
                percent={progress.percent}
                phase={progress.phase}
                detail={progress.detail}
                compact
              />
              <div className="run-status-meta">
                <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />
              </div>
            </div>
            <div className={`run-layout ${isHpcMode && !result ? "run-layout-hpc-telemetry" : result ? "" : "run-layout-single"}`}>
              <LiveOutputPanel
                title={isRunning ? "Running..." : "Output"}
                output={output}
                placeholder="Starting calculation..."
                totalLineCount={outputLineCount}
                visibleLineCount={visibleOutputLineCount}
              />

              {isHpcMode && !result && (
                <RemoteUtilizationPanel
                  enabled={isRunning || activeTask?.status === "running"}
                  profileId={activeHpcProfile?.id ?? null}
                  remoteJobId={activeTask?.hpc?.remote_job_id ?? null}
                  remoteNode={activeTask?.hpc?.remote_node ?? null}
                  resourceType={activeTask?.hpc?.hpc_resource_type ?? hpcResources.resource_type}
                />
              )}

              {result && (
                <div className="results-panel">
                  <h3>Results</h3>
                  <div className="result-grid">
                    <div className="result-item">
                      <label>Status</label>
                      <span className={result.converged ? "success" : "error"}>
                        {result.converged ? "Converged" : "Not Converged"}
                      </span>
                    </div>
                    {result.total_energy && (
                      <div className="result-item">
                        <label>Total Energy</label>
                        <span>{result.total_energy.toFixed(8)} Ry</span>
                      </div>
                    )}
                    {result.n_scf_steps && (
                      <div className="result-item">
                        <label>SCF Iterations</label>
                        <span>{result.n_scf_steps}</span>
                      </div>
                    )}
                    {result.wall_time_seconds && (
                      <div className="result-item">
                        <label>Wall Time</label>
                        <span>{result.wall_time_seconds.toFixed(1)} s</span>
                      </div>
                    )}
                  </div>
                  {resultSaved && (
                    <div className="save-status save-status-inline">
                      <span className="saved">Saved to project</span>
                    </div>
                  )}
                  {autoSaveExpected && isSaving && (
                    <div className="save-status save-status-inline">
                      <span>Saving to project...</span>
                    </div>
                  )}
                </div>
              )}
            </div>

            <div className="run-actions">
              <button onClick={() => setStep("configure")}>
                ← Back to Configure
              </button>
              {result && (
                <>
                  {!autoSaveExpected && (
                    <button
                      onClick={() => setShowSaveDialog(true)}
                      className="save-project-btn"
                      disabled={isSaving}
                    >
                      {resultSaved ? "Saved" : "Save to Project"}
                    </button>
                  )}
                  <button onClick={() => setStep("import")} className="new-calc-btn">
                    New Calculation
                  </button>
                </>
              )}
            </div>
          </div>
        )}
      </div>

      {/* Save to Project Dialog */}
      {!autoSaveExpected && showSaveDialog && crystalData && calculationData && (
        <SaveToProjectDialog
          isOpen={showSaveDialog}
          onClose={() => setShowSaveDialog(false)}
          onSaved={() => {
            setShowSaveDialog(false);
            setResultSaved(true);
          }}
          calculationData={calculationData}
          cifData={{
            filename: cifFilename,
            formula: crystalData.formula_sum || crystalData.formula_structural || "Unknown",
            content: cifContent,
            crystal_data: crystalData,
          }}
          workingDir={WORK_DIR}
          projectContext={projectContext || undefined}
        />
      )}
    </div>
  );
}
