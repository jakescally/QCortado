// SCF Calculation Wizard - Import CIF, configure, and run SCF calculations

import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import {
  CrystalData,
  ELEMENT_MASSES,
  OptimizedStructureOption,
  QePositionUnit,
  SavedCellSummary,
  SavedStructureData,
  SCFPreset,
} from "../lib/types";
import { parseCIF } from "../lib/cifParser";
import { UnitCellViewer } from "./UnitCellViewer";
import { SaveToProjectDialog } from "./SaveToProjectDialog";
import { getPrimitiveCell, PrimitiveCell } from "../lib/primitiveCell";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { isPhononReadyScf } from "../lib/phononReady";

// Tooltip component for help icons
function Tooltip({ text }: { text: string }) {
  return (
    <span className="tooltip-container">
      <span className="tooltip-icon">?</span>
      <span className="tooltip-text">{text}</span>
    </span>
  );
}

interface InitialCifData {
  cifId: string;
  crystalData: CrystalData;
  cifContent: string;
  filename: string;
  projectId: string;
}

interface SCFWizardProps {
  onBack: () => void;
  qePath: string;
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
  sourceDescriptor: { type: "cif" | "optimization"; calc_id?: string };
}

export function SCFWizard({
  onBack,
  qePath,
  initialCif,
  initialPreset,
  presetLock,
  optimizedStructures = [],
  reconnectTaskId,
}: SCFWizardProps) {
  const taskContext = useTaskContext();
  // Track background task for this wizard
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  // If we have initial CIF data, skip to configure step
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : (initialCif ? "configure" : "import"));
  const [crystalData, setCrystalData] = useState<CrystalData | null>(initialCif?.crystalData || null);
  const [cifFilename, setCifFilename] = useState<string>(initialCif?.filename || "");
  const [cifContent, setCifContent] = useState<string>(initialCif?.cifContent || "");
  const [error, setError] = useState<string | null>(null);

  // Track project context for saving
  const [projectContext] = useState<{ projectId: string; cifId: string } | null>(
    initialCif ? { projectId: initialCif.projectId, cifId: initialCif.cifId } : null
  );
  const [pseudopotentials, setPseudopotentials] = useState<string[]>([]);
  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});

  const WORK_DIR = "/tmp/qcortado_work";

  // Save dialog state
  const [showSaveDialog, setShowSaveDialog] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [calcEndTime, setCalcEndTime] = useState<string>("");
  const [generatedInput, setGeneratedInput] = useState<string>("");
  const [resultSaved, setResultSaved] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [runSourceStructure, setRunSourceStructure] = useState<SavedStructureData | null>(null);
  const [runSourceDescriptor, setRunSourceDescriptor] = useState<{ type: "cif" | "optimization"; calc_id?: string } | null>(null);
  const [cellViewMode, setCellViewMode] = useState<CellViewMode>("conventional");

  // Exit confirmation dialog
  const [showExitConfirm, setShowExitConfirm] = useState(false);
  const [pendingExitAction, setPendingExitAction] = useState<(() => void) | null>(null);

  // Ref for output auto-follow (only when user is at the bottom)
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);

  const handleOutputScroll = () => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followOutputRef.current = distanceToBottom <= 16;
  };

  const [config, setConfig] = useState<SCFConfig>({
    // Calculation type
    calculation: "scf",

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
    smearing: "gaussian",
    degauss: 0.01,
    nbnd: null,
    tot_charge: 0,

    // Magnetism & Spin
    nspin: 1,
    noncolin: false,
    lspinorb: false,
    starting_magnetization: {},
    tot_magnetization: null,
    constrained_magnetization: "none",

    // SCF Convergence
    conv_thr: 1e-6,
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
  });
  const [selectedPreset, setSelectedPreset] = useState<SCFPreset | null>(null);
  const [structureSource, setStructureSource] = useState<string>("cif");
  const selectedOptimizedStructure =
    structureSource === "cif"
      ? null
      : optimizedStructures.find((option) => option.calcId === structureSource) || null;
  const phononPresetDisabled = structureSource === "cif";
  const phononPresetDisabledMessage = "Select an optimized structure (or run a structure optimization) to start a phonon-ready calculation.";

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
            conv_thr: 1e-6,
            disk_io: "low",
          };
        case "phonon":
          return {
            ...prev,
            calculation: "scf",
            conv_thr: 1e-12,
            disk_io: "medium",
          };
        case "relax":
          return {
            ...prev,
            calculation: "vcrelax",
            forc_conv_thr: 1e-4,
            etot_conv_thr: 1e-5,
            press: 0,
            tprnfor: true,
            tstress: true,
          };
        default:
          return prev;
      }
    });
    if (preset === "phonon") {
      setConvThrInput("1e-12");
    } else if (preset === "standard") {
      setConvThrInput("1e-6");
    }
  }, [phononPresetDisabled]);

  useEffect(() => {
    if (initialPreset) {
      applyPreset(initialPreset);
    }
  }, [initialPreset, applyPreset]);

  // Reconnect to a running/completed background task
  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) {
      // Task may not be loaded yet, try reconnecting
      taskContext.reconnectToTask(activeTaskId);
      return;
    }

    // Sync wizard state from task
    setIsRunning(task.status === "running");
    setOutput(task.output.join("\n") + "\n");
    setProgress(task.progress);
    setCalcStartTime(task.startedAt);

    if (task.status === "completed" && task.result) {
      setResult(task.result);
      setCalcEndTime(new Date().toISOString());
      setStep("results");
    } else if (task.status === "failed" || task.status === "cancelled") {
      setError(task.error || "Task failed");
    } else {
      setStep("run");
    }
  }, [activeTaskId, taskContext.getTask(activeTaskId ?? "")?.status]);

  // When active task updates, sync output/progress
  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task || task.status !== "running") return;

    setOutput(task.output.join("\n") + "\n");
    setProgress(task.progress);

    if (followOutputRef.current && outputRef.current) {
      outputRef.current.scrollTop = outputRef.current.scrollHeight;
    }
  }, [activeTaskId, taskContext.getTask(activeTaskId ?? "")?.output.length]);

  const lockedPreset = presetLock && initialPreset ? initialPreset : null;
  const isOptimizationWizard = lockedPreset === "relax";
  const wizardTitle = isOptimizationWizard ? "Structure Optimization Wizard" : "SCF Calculation Wizard";
  const showPresetRow = !lockedPreset || lockedPreset === "relax";
  const showStandardPreset = !lockedPreset || lockedPreset === "standard";
  const showPhononPreset = !lockedPreset || lockedPreset === "phonon";
  const showRelaxPreset = lockedPreset === "relax";
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
    if (!crystalData) return null;
    const primitive = getPrimitiveCell(crystalData);
    if (!primitive) return null;
    const BOHR_TO_ANGSTROM = 0.529177;
    const a = primitive.celldm1 * BOHR_TO_ANGSTROM;

    if (primitive.ibrav === "cubic_f") {
      const v1: [number, number, number] = [0, a / 2, a / 2];
      const v2: [number, number, number] = [a / 2, 0, a / 2];
      const v3: [number, number, number] = [a / 2, a / 2, 0];
      return calculateMetricsFromVectors(v1, v2, v3);
    }

    if (primitive.ibrav === "cubic_i") {
      const v1: [number, number, number] = [a / 2, a / 2, -a / 2];
      const v2: [number, number, number] = [-a / 2, a / 2, a / 2];
      const v3: [number, number, number] = [a / 2, -a / 2, a / 2];
      return calculateMetricsFromVectors(v1, v2, v3);
    }

    return {
      a,
      b: a,
      c: a,
      alpha: 90,
      beta: 90,
      gamma: 90,
    };
  }, [crystalData]);

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
  const [convThrInput, setConvThrInput] = useState("1e-6");

  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState<string>("");
  const [result, setResult] = useState<any>(null);
  const [pseudoDir, setPseudoDir] = useState<string>("");
  const [pseudoError, setPseudoError] = useState<string | null>(null);
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "SCF iterations",
  });

  // SSSP data for recommended cutoffs and pseudopotentials
  interface SSSPElementData {
    filename: string;
    md5?: string;
    pseudopotential?: string;
    cutoff_wfc: number;
    cutoff_rho: number;
  }
  const [ssspData, setSsspData] = useState<Record<string, SSSPElementData> | null>(null);
  const [ssspMissing, setSsspMissing] = useState(false);

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);

  // Load available pseudopotentials and SSSP data
  useEffect(() => {
    async function loadPseudos() {
      try {
        // Compute pseudo directory from QE path
        const computedPseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
        setPseudoDir(computedPseudoDir);

        const pseudos = await invoke<string[]>("list_pseudopotentials", {
          pseudoDir: computedPseudoDir,
        });
        setPseudopotentials(pseudos);
        setPseudoError(null);

        // Try to load SSSP data
        try {
          const sssp = await invoke<Record<string, SSSPElementData>>("load_sssp_data", {
            pseudoDir: computedPseudoDir,
          });
          setSsspData(sssp);
          setSsspMissing(false);
        } catch {
          setSsspData(null);
          setSsspMissing(true);
        }
      } catch (e) {
        console.error("Failed to load pseudopotentials:", e);
        setPseudoError(String(e));
        setPseudopotentials([]);
      }
    }
    loadPseudos();
  }, [qePath]);

  // Auto-scroll output only if user is at the bottom
  useEffect(() => {
    const el = outputRef.current;
    if (!el || !followOutputRef.current) return;
    el.scrollTop = el.scrollHeight;
  }, [output]);

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

  // Auto-select pseudopotentials and set cutoffs when structure source or SSSP data changes
  useEffect(() => {
    if (pseudopotentials.length > 0) {
      const selected: Record<string, string> = {};
      const elements = getUniqueElements();
      if (elements.length === 0) return;

      let maxWfc = 0;
      let maxRho = 0;

      for (const element of elements) {
        // First, try to use SSSP recommended pseudopotential
        if (ssspData && ssspData[element]) {
          const ssspEntry = ssspData[element];
          // Check if the SSSP recommended file exists in our pseudopotentials
          if (pseudopotentials.some((p) => p.toLowerCase() === ssspEntry.filename.toLowerCase())) {
            selected[element] = ssspEntry.filename;
            maxWfc = Math.max(maxWfc, ssspEntry.cutoff_wfc);
            maxRho = Math.max(maxRho, ssspEntry.cutoff_rho);
            continue;
          }
        }

        // Fallback: find a matching pseudopotential by element name
        const matches = pseudopotentials.filter((p) =>
          pseudoMatchesElement(p, element)
        );
        if (matches.length > 0) {
          // Prefer PBE pseudopotentials
          const pbe = matches.find((m) => m.toLowerCase().includes("pbe"));
          selected[element] = pbe || matches[0];
        }
      }

      setSelectedPseudos(selected);

      // Update cutoffs if we got SSSP values
      if (maxWfc > 0 && maxRho > 0) {
        setConfig((prev) => ({
          ...prev,
          ecutwfc: maxWfc,
          ecutrho: maxRho,
        }));
      }
    }
  }, [crystalData, pseudopotentials, ssspData, structureSource, optimizedStructures]);

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

  function canRun(): boolean {
    if (!crystalData) return false;
    const elements = getUniqueElements();
    return elements.every((el) => selectedPseudos[el]);
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

    // Check if this is an FCC structure that should use primitive cell
    // Using primitive cell with ibrav gives correct symmetry and faster calculations
    const primitiveCell: PrimitiveCell | null = selectedOptimizedStructure ? null : getPrimitiveCell(crystalData);
    const usePrimitiveCell = !selectedOptimizedStructure && primitiveCell !== null;

    if (usePrimitiveCell && primitiveCell) {
      console.log(`Using primitive cell: ibrav=${primitiveCell.ibravNumeric} (${primitiveCell.ibrav}), celldm(1)=${primitiveCell.celldm1.toFixed(4)} Bohr, ${primitiveCell.nat} atoms`);
    }

    // Build species list (same for both primitive and conventional)
    const speciesList = elements.map((el) => ({
      symbol: el,
      mass: ELEMENT_MASSES[el] || 1.0,
      pseudopotential: selectedPseudos[el],
      starting_magnetization: config.starting_magnetization[el] || 0,
      hubbard_u: config.lda_plus_u ? (config.hubbard_u[el] || 0) : undefined,
      hubbard_j: config.lda_plus_u && config.lda_plus_u_kind > 0 ? (config.hubbard_j[el] || 0) : undefined,
    }));

    // Build system configuration based on cell type
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
      lda_plus_u: config.lda_plus_u || undefined,
      lda_plus_u_kind: config.lda_plus_u ? config.lda_plus_u_kind : undefined,
      // Van der Waals
      vdw_corr: config.vdw_corr !== "none" ? config.vdw_corr : undefined,
      // Isolated systems
      assume_isolated: config.assume_isolated !== "none" ? config.assume_isolated : undefined,
      // XC functional override
      input_dft: config.input_dft || undefined,
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
    } else if (usePrimitiveCell && primitiveCell) {
      // Primitive cell with ibrav (e.g., ibrav=2 for FCC)
      systemConfig.ibrav = primitiveCell.ibrav;
      // celldm array: [a, b/a, c/a, cos(alpha), cos(beta), cos(gamma)]
      // For cubic, only celldm(1) = a is needed
      systemConfig.celldm = [primitiveCell.celldm1, 0, 0, 0, 0, 0];
      systemConfig.cell_parameters = null;
      systemConfig.cell_units = null;
      systemConfig.atoms = primitiveCell.atoms.map((atom) => ({
        symbol: atom.symbol,
        position: atom.position,
        if_pos: [true, true, true],
      }));
    } else {
      // Conventional cell with ibrav=0 (fallback for non-FCC structures)
      systemConfig.ibrav = "free";
      systemConfig.celldm = null;
      systemConfig.cell_parameters = [
        [crystalData.cell_length_a.value, 0, 0],
        [
          crystalData.cell_length_b.value *
            Math.cos((crystalData.cell_angle_gamma.value * Math.PI) / 180),
          crystalData.cell_length_b.value *
            Math.sin((crystalData.cell_angle_gamma.value * Math.PI) / 180),
          0,
        ],
        calculateCVector(crystalData),
      ];
      systemConfig.cell_units = "angstrom";
      systemConfig.atoms = crystalData.atom_sites.map((site) => ({
        symbol: getBaseElement(site.type_symbol),
        position: [site.fract_x, site.fract_y, site.fract_z],
        if_pos: [true, true, true],
      }));
    }

    // Build the calculation configuration with all options
    const calculation: any = {
      calculation: config.calculation,
      prefix: "qcortado_scf",
      outdir: "./tmp",
      pseudo_dir: qePath.replace("/bin", "/pseudo"),
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
        mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      },
      sourceStructure,
      sourceDescriptor,
    };
  }

  async function handleRun() {
    if (!crystalData || !canRun()) return;
    if (hasExternalRunningTask) {
      setError("Another task is currently running. Queue this task or wait for completion.");
      return;
    }

    setIsRunning(true);
    followOutputRef.current = true;
    setOutput("");
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
      const { inputText, taskLabel, sourceStructure, sourceDescriptor } = plan;
      setRunSourceStructure(sourceStructure);
      setRunSourceDescriptor(sourceDescriptor);
      setGeneratedInput(inputText);
      setOutput(`=== Generated Input ===\n${inputText}\n\n=== Running pw.x ===\n`);

      // Start as a background task
      const taskId = await taskContext.startTask("scf", plan.taskParams, `SCF - ${taskLabel}`);
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
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      // Auto-save to project if we have project context
      if (projectContext) {
        try {
          setIsSaving(true);
          const calcData = buildCalculationData(
            calcResult,
            startTime,
            endTime,
            inputText,
            sourceStructure,
            sourceDescriptor,
          );
          await invoke("save_calculation", {
            projectId: projectContext.projectId,
            cifId: projectContext.cifId,
            calcData,
            workingDir: WORK_DIR,
          });
          setResultSaved(true);
        } catch (saveError) {
          console.error("Failed to auto-save calculation:", saveError);
          setError(`Failed to auto-save calculation: ${saveError}`);
        } finally {
          setIsSaving(false);
        }
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

  // Handle exit attempts - show confirmation if results not saved
  function handleExitAttempt(action: () => void) {
    if (result && !resultSaved) {
      setPendingExitAction(() => action);
      setShowExitConfirm(true);
    } else {
      action();
    }
  }

  function confirmExit() {
    if (pendingExitAction) {
      pendingExitAction();
    }
    setShowExitConfirm(false);
    setPendingExitAction(null);
  }

  function cancelExit() {
    setShowExitConfirm(false);
    setPendingExitAction(null);
  }

  function calculateCVector(data: CrystalData): [number, number, number] {
    const c = data.cell_length_c.value;
    const alpha = (data.cell_angle_alpha.value * Math.PI) / 180;
    const beta = (data.cell_angle_beta.value * Math.PI) / 180;
    const gamma = (data.cell_angle_gamma.value * Math.PI) / 180;

    const cx = c * Math.cos(beta);
    const cy = (c * (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma))) / Math.sin(gamma);
    const cz = Math.sqrt(c * c - cx * cx - cy * cy);

    return [cx, cy, cz];
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
      cell_parameters: [
        [data.cell_length_a.value, 0, 0],
        [
          data.cell_length_b.value * Math.cos((data.cell_angle_gamma.value * Math.PI) / 180),
          data.cell_length_b.value * Math.sin((data.cell_angle_gamma.value * Math.PI) / 180),
          0,
        ],
        calculateCVector(data),
      ],
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
        conv_thr: config.conv_thr,
        mixing_beta: config.mixing_beta,
        selected_pseudos: selectedPseudos,
        structure_source: sourceDescriptor,
        source_structure: sourceStructure,
        // Feature flags for tags
        nspin: config.nspin,
        lspinorb: config.lspinorb,
        lda_plus_u: config.lda_plus_u,
        vdw_corr: config.vdw_corr,
        // Relaxation parameters
        forc_conv_thr: config.forc_conv_thr,
        etot_conv_thr: config.etot_conv_thr,
        press: config.press,
        optimization_mode: isOptimization ? config.calculation : null,
        optimized_structure: optimizedStructure,
        optimized_cell_summary: optimizedCellSummary,
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
      )
    : null;

  return (
    <div className={`wizard-container wizard-step-${step}`}>
      <div className="wizard-header">
        <button className="back-btn" onClick={() => handleExitAttempt(onBack)}>
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
                      <button
                        type="button"
                        className={`preset-btn ${cellViewMode === "conventional" ? "active" : ""}`}
                        onClick={() => setCellViewMode("conventional")}
                        title="Show conventional CIF lattice parameters"
                      >
                        Conventional
                      </button>
                      <button
                        type="button"
                        className={`preset-btn ${cellViewMode === "primitive" ? "active" : ""}`}
                        onClick={() => setCellViewMode("primitive")}
                        title="Show primitive lattice parameters used by QE (when detected)"
                      >
                        Primitive
                      </button>
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
                    <Tooltip text="Pseudopotentials approximate the effect of core electrons (those tightly bound to the nucleus) so the calculation only needs to explicitly handle valence electrons. This dramatically reduces computation cost while maintaining accuracy for chemical bonding and material properties. Each element needs its own pseudopotential file, typically generated with a specific exchange-correlation functional (like PBE)." />
                  </h3>
                  <p className="pseudo-dir-info">
                    Directory: <code>{pseudoDir}</code>
                    {pseudopotentials.length > 0 && (
                      <span className="pseudo-count"> ({pseudopotentials.length} files)</span>
                    )}
                  </p>
                  {pseudoError && (
                    <div className="pseudo-error">{pseudoError}</div>
                  )}
                  <div className="pseudo-list">
                    {getUniqueElements().map((el) => {
                      const matchingPseudos = pseudopotentials.filter((p) =>
                        pseudoMatchesElement(p, el)
                      );
                      return (
                        <div key={el} className="pseudo-row">
                          <label>{el}</label>
                          <select
                            value={selectedPseudos[el] || ""}
                            onChange={(e) =>
                              setSelectedPseudos((prev) => ({
                                ...prev,
                                [el]: e.target.value,
                              }))
                            }
                          >
                            <option value="">
                              {matchingPseudos.length === 0
                                ? `No ${el} pseudopotentials found`
                                : "Select..."}
                            </option>
                            {matchingPseudos.map((p) => (
                              <option key={p} value={p}>
                                {p}
                              </option>
                            ))}
                          </select>
                        </div>
                      );
                    })}
                  </div>
                  {getUniqueElements().some(
                    (el) =>
                      !pseudopotentials.some((p) => pseudoMatchesElement(p, el))
                  ) && (
                    <p className="pseudo-hint">
                      Missing pseudopotentials? Download them from the{" "}
                      <a
                        href="https://www.materialscloud.org/discover/sssp/table/precision"
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                        SSSP Precision Library
                      </a>{" "}
                      and place them in your QE pseudo directory.
                    </p>
                  )}
                </section>

                {/* Basic Parameters */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("basic")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.basic ? "expanded" : ""}`}>▶</span>
                    Basic Parameters
                    <Tooltip text="Essential parameters for any SCF calculation: energy cutoffs and k-point sampling." />
                  </h3>
                  {expandedSections.basic && (
                    <div className="param-grid">
                      {showPresetRow && (
                        <div className="param-row full-width">
                          <label>Presets</label>
                          <div className="preset-buttons">
                            {showStandardPreset && (
                              <button
                                type="button"
                                className={`preset-btn ${selectedPreset === "standard" ? "active" : ""}`}
                                onClick={() => applyPreset("standard")}
                                title="Standard SCF calculation"
                              >
                                Standard
                              </button>
                            )}
                            {showPhononPreset && (
                              <span className="preset-tooltip-container">
                                <button
                                  type="button"
                                  className={`preset-btn preset-phonon ${selectedPreset === "phonon" && !phononPresetDisabled ? "active" : ""} ${phononPresetDisabled ? "preset-disabled" : ""}`}
                                  onClick={() => applyPreset("phonon")}
                                  aria-disabled={phononPresetDisabled}
                                >
                                  Phonon-Ready
                                </button>
                                {phononPresetDisabled && (
                                  <span className="tooltip-text preset-disabled-tooltip">
                                    {phononPresetDisabledMessage}
                                  </span>
                                )}
                              </span>
                            )}
                            {showRelaxPreset && (
                              <button
                                type="button"
                                className={`preset-btn preset-relax ${selectedPreset === "relax" ? "active" : ""}`}
                                onClick={() => applyPreset("relax")}
                                title="Variable-cell relaxation for structure optimization"
                              >
                                Optimize Structure
                              </button>
                            )}
                          </div>
                        </div>
                      )}

                      {!isOptimizationWizard && (
                        <div className="param-row">
                          <label>
                            Structure Source
                            <Tooltip text="Choose the geometry used for this SCF. 'From CIF' uses the original imported structure; optimized structures come from saved geometry optimization runs." />
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
                      )}

                      <div className="param-row">
                        <label>
                          Calculation Type
                          <Tooltip text="scf: ground state energy only. relax: optimize atomic positions. vcrelax: optimize both positions and cell shape/size." />
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

                      {/* Relaxation parameters - shown only for relax/vcrelax */}
                      {(config.calculation === "relax" || config.calculation === "vcrelax") && (
                        <>
                          <div className="param-row">
                            <label>
                              Force Convergence
                              <Tooltip text="Convergence threshold for forces (Ry/Bohr). Relaxation stops when all forces are below this value. Use 1e-4 for rough optimization, 1e-5 for phonon calculations." />
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
                              <Tooltip text="Convergence threshold for total energy change between ionic steps (Ry)." />
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
                            <Tooltip text="Target external pressure in kbar. Use 0 for ambient pressure, positive for compression." />
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
                          <Tooltip text="Energy cutoff for plane-wave expansion of wavefunctions (in Rydberg). Higher = more accurate but slower. Typical: 30-80 Ry." />
                        </label>
                        <div className="param-input-group">
                          <input type="number" value={config.ecutwfc}
                            onChange={(e) => setConfig((prev) => ({ ...prev, ecutwfc: parseFloat(e.target.value) }))} />
                          <span className="param-unit">Ry</span>
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          Charge Density Cutoff
                          <Tooltip text="Energy cutoff for charge density (in Rydberg). Use 4x ecutwfc for NC, 8-12x for US/PAW pseudopotentials." />
                        </label>
                        <div className="param-input-group">
                          <input type="number" value={config.ecutrho}
                            onChange={(e) => setConfig((prev) => ({ ...prev, ecutrho: parseFloat(e.target.value) }))} />
                          <span className="param-unit">Ry</span>
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          K-point Grid
                          <Tooltip text="Monkhorst-Pack grid for Brillouin zone sampling. Denser = more accurate. Metals need denser grids." />
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
                          <Tooltip text="Shift the k-point grid. Use (1,1,1) for metals to avoid high-symmetry points, (0,0,0) for insulators." />
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
                      Cutoff values are defaults. For optimized values, download the{" "}
                      <a href="https://www.materialscloud.org/discover/sssp/table/precision" target="_blank" rel="noopener noreferrer">SSSP library</a>.
                    </p>
                  )}
                  {!ssspMissing && ssspData && expandedSections.basic && (
                    <p className="sssp-success">Cutoffs auto-configured from SSSP library.</p>
                  )}
                </section>

                {/* Electronic Structure */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("electronic")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.electronic ? "expanded" : ""}`}>▶</span>
                    Electronic Structure
                    <Tooltip text="Control how electronic states are occupied and how the Fermi level is determined." />
                  </h3>
                  {expandedSections.electronic && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Occupations
                          <Tooltip text="How to fill electronic states. 'smearing' for metals, 'tetrahedra' for accurate DOS, 'fixed' for insulators with gap." />
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
                              <Tooltip text="gaussian: simple, good for most cases. methfessel-paxton: better for metals. marzari-vanderbilt: 'cold smearing', good for forces. fermi-dirac: physical but needs higher temp." />
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
                              <Tooltip text="Smearing width in Ry. Smaller = more accurate but harder to converge. Typical: 0.005-0.02 Ry for metals." />
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
                          <Tooltip text="Total number of Kohn-Sham states. Leave empty for automatic (occupied + some empty). Increase for accurate DOS or optical properties." />
                        </label>
                        <div className="param-input-group">
                          <input type="number" value={config.nbnd ?? ""} placeholder="auto"
                            onChange={(e) => setConfig((prev) => ({ ...prev, nbnd: e.target.value ? parseInt(e.target.value) : null }))} />
                        </div>
                      </div>
                      <div className="param-row">
                        <label>
                          Total Charge
                          <Tooltip text="Total charge of the system in units of e. Positive = remove electrons, negative = add electrons. 0 = neutral." />
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
                          <Tooltip text="Override the exchange-correlation functional from pseudopotentials. Leave empty to use pseudopotential default. Examples: PBE, PBEsol, SCAN, HSE, B3LYP." />
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
                    <Tooltip text="Enable spin-polarized calculations for magnetic materials and spin-orbit coupling for heavy elements." />
                  </h3>
                  {expandedSections.magnetism && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Spin Polarization
                          <Tooltip text="nspin=1: non-magnetic. nspin=2: collinear magnetic (spin up/down). nspin=4: non-collinear (required for spin-orbit)." />
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
                            <Tooltip text="Include spin-orbit interaction. Required for topological properties and heavy elements. Needs fully-relativistic pseudopotentials." />
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
                            <Tooltip text="Fix total magnetization (N_up - N_down). Leave empty for unconstrained." />
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
                              <Tooltip text="Initial spin polarization for each element (-1 to +1). Helps SCF converge to correct magnetic state." />
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
                    <Tooltip text="Parameters controlling how the self-consistent field iteration converges." />
                  </h3>
                  {expandedSections.convergence && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Convergence Threshold
                          <Tooltip text="SCF stops when energy change falls below this value (Ry). 1e-6 for most cases, 1e-8+ for phonons." />
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
                          <Tooltip text="Maximum number of SCF iterations before giving up." />
                        </label>
                        <input type="number" value={config.electron_maxstep}
                          onChange={(e) => setConfig((prev) => ({ ...prev, electron_maxstep: parseInt(e.target.value) }))} />
                      </div>
                      <div className="param-row">
                        <label>
                          Mixing Mode
                          <Tooltip text="How to mix old and new charge density. 'plain' is default, 'TF' or 'local-TF' can help metals converge." />
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
                          <Tooltip text="Fraction of new density to mix in (0-1). Lower = more stable but slower. 0.7 default, use 0.1-0.3 for difficult cases." />
                        </label>
                        <input type="number" step="0.05" min={0} max={1} value={config.mixing_beta}
                          onChange={(e) => setConfig((prev) => ({ ...prev, mixing_beta: parseFloat(e.target.value) }))} />
                      </div>
                      <div className="param-row">
                        <label>
                          Mixing Dimensions
                          <Tooltip text="Number of previous iterations used in Broyden mixing. Higher can help difficult cases converge." />
                        </label>
                        <input type="number" min={2} max={20} value={config.mixing_ndim}
                          onChange={(e) => setConfig((prev) => ({ ...prev, mixing_ndim: parseInt(e.target.value) }))} />
                      </div>
                      <div className="param-row">
                        <label>
                          Diagonalization
                          <Tooltip text="Algorithm for solving eigenvalue problem. 'david' is fast for most cases, 'cg' for difficult cases, 'ppcg'/'paro' for large parallel runs." />
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
                          <Tooltip text="'atomic': superposition of atomic potentials (default). 'file': read from previous calculation." />
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
                          <Tooltip text="Initial guess for wavefunctions. 'atomic': atomic orbitals. 'atomic+random': with randomization. 'random': fully random." />
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
                    <Tooltip text="Add Hubbard U correction to improve description of localized d and f electrons in transition metals and lanthanides." />
                  </h3>
                  {expandedSections.dftu && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Enable DFT+U
                          <Tooltip text="Apply Hubbard U correction. Essential for correlated materials like transition metal oxides." />
                        </label>
                        <label className="toggle-label">
                          <input type="checkbox" checked={config.lda_plus_u}
                            onChange={(e) => setConfig((prev) => ({ ...prev, lda_plus_u: e.target.checked }))} />
                          <span>Enable Hubbard correction</span>
                        </label>
                      </div>
                      {config.lda_plus_u && (
                        <>
                          <div className="param-row">
                            <label>
                              DFT+U Type
                              <Tooltip text="0: simplified rotationally invariant. 1: rotationally invariant with J. 2: DFT+U+J (full)." />
                            </label>
                            <select value={config.lda_plus_u_kind}
                              onChange={(e) => setConfig((prev) => ({ ...prev, lda_plus_u_kind: parseInt(e.target.value) as 0 | 1 | 2 }))}>
                              <option value={0}>Simplified (kind=0)</option>
                              <option value={1}>Rotationally invariant (kind=1)</option>
                              <option value={2}>Full DFT+U+J (kind=2)</option>
                            </select>
                          </div>
                          <div className="param-row full-width">
                            <label>Hubbard U values (per element, in eV)</label>
                          </div>
                          <div className="hubbard-grid">
                            {getUniqueElements().map((el) => (
                              <div key={el} className="hubbard-row">
                                <label>{el}</label>
                                <span>U:</span>
                                <input type="number" step="0.1" min={0}
                                  value={config.hubbard_u[el] ?? 0}
                                  onChange={(e) => setConfig((prev) => ({
                                    ...prev,
                                    hubbard_u: { ...prev.hubbard_u, [el]: parseFloat(e.target.value) || 0 }
                                  }))} />
                                {config.lda_plus_u_kind > 0 && (
                                  <>
                                    <span>J:</span>
                                    <input type="number" step="0.1" min={0}
                                      value={config.hubbard_j[el] ?? 0}
                                      onChange={(e) => setConfig((prev) => ({
                                        ...prev,
                                        hubbard_j: { ...prev.hubbard_j, [el]: parseFloat(e.target.value) || 0 }
                                      }))} />
                                  </>
                                )}
                              </div>
                            ))}
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
                    <Tooltip text="Add dispersion corrections for systems where van der Waals interactions are important (layered materials, molecules, etc.)." />
                  </h3>
                  {expandedSections.vdw && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          vdW Correction Method
                          <Tooltip text="None: standard DFT. DFT-D2/D3: Grimme's empirical corrections. TS-vdW: Tkatchenko-Scheffler. XDM: exchange-hole dipole moment." />
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
                    <Tooltip text="Additional settings for special cases: isolated systems, output control, etc." />
                  </h3>
                  {expandedSections.advanced && (
                    <div className="param-grid">
                      <div className="param-row">
                        <label>
                          Isolated System
                          <Tooltip text="For molecules or slabs, remove spurious interactions with periodic images. 'martyna-tuckerman' for molecules, 'esm' or '2D' for slabs." />
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
                          <Tooltip text="Amount of output. 'high' gives more detail for debugging." />
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
                          <Tooltip text="How much to write to disk. 'low' saves space, 'nowf' doesn't save wavefunctions (can't restart)." />
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

                {/* MPI Parallelization */}
                <section className="config-section">
                  <h3>
                    Parallelization
                    <Tooltip text="MPI (Message Passing Interface) allows running calculations in parallel across multiple CPU cores, significantly speeding up large calculations. Your Quantum ESPRESSO installation must be compiled with MPI support for this to work." />
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
                          <Tooltip text="Number of parallel MPI processes to use. More processes can speed up large calculations but have overhead for small systems. Using too many processes for a small system may actually slow things down." />
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

                <div className="run-btn-group">
                  {projectContext && (
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
                    disabled={!canRun() || hasExternalRunningTask}
                  >
                    {mpiEnabled && mpiProcs > 1
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
          <div className="run-step">
            <ProgressBar
              status={progress.status}
              percent={progress.percent}
              phase={progress.phase}
              detail={progress.detail}
            />
            <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />
            <div className="run-layout">
              <div className="output-panel">
                <h3>{isRunning ? "Running..." : "Output"}</h3>
                <pre className="output-text" ref={outputRef} onScroll={handleOutputScroll}>{output}</pre>
              </div>

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
                </div>
              )}
            </div>

            <div className="run-actions">
              <button onClick={() => handleExitAttempt(() => setStep("configure"))}>
                ← Back to Configure
              </button>
              {result && (
                <>
                  {!projectContext && (
                    <button
                      onClick={() => setShowSaveDialog(true)}
                      className="save-project-btn"
                      disabled={isSaving}
                    >
                      {resultSaved ? "Saved" : "Save to Project"}
                    </button>
                  )}
                  <button onClick={() => handleExitAttempt(() => setStep("import"))} className="new-calc-btn">
                    New Calculation
                  </button>
                </>
              )}
            </div>
          </div>
        )}
      </div>

      {/* Save to Project Dialog */}
      {!projectContext && showSaveDialog && crystalData && calculationData && (
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

      {/* Exit Confirmation Dialog */}
      {showExitConfirm && (
        <div className="dialog-overlay" onClick={cancelExit}>
          <div className="dialog-content dialog-small" onClick={(e) => e.stopPropagation()}>
            <div className="dialog-header">
              <h2>Unsaved Results</h2>
              <button className="dialog-close" onClick={cancelExit}>
                &times;
              </button>
            </div>

            <div className="dialog-body">
              <p className="exit-warning">
                {projectContext
                  ? "Your calculation results have not finished auto-saving to this project yet."
                  : "Your calculation results have not been saved to a project."}
                {" "}
                If you leave now, the results may be lost.
              </p>
              <p className="exit-hint">
                {projectContext
                  ? "Wait for auto-save to complete, or choose \"Leave Anyway\" to discard this run."
                  : "Click \"Save to Project\" to keep your results, or \"Leave Anyway\" to discard them."}
              </p>
            </div>

            <div className="dialog-footer">
              <button className="dialog-btn cancel" onClick={cancelExit}>
                Stay Here
              </button>
              <button className="dialog-btn delete" onClick={confirmExit}>
                Leave Anyway
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
