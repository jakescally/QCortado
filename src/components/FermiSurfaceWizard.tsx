import { useState, useEffect, useRef, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData, ELEMENT_MASSES } from "../lib/types";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { getPrimitiveCell, PrimitiveCell } from "../lib/primitiveCell";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { useTaskContext } from "../lib/TaskContext";

interface CalculationRun {
  id: string;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    total_energy: number | null;
    fermi_energy: number | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

interface FermiSurfaceWizardProps {
  onBack: () => void;
  qePath: string;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "parameters" | "run" | "results";
const FERMI_WORK_DIR = "/tmp/qcortado_fermi_surface";

interface BxsfFileData {
  file_name: string;
  size_bytes: number;
}

interface FermiSurfaceData {
  k_grid: [number, number, number];
  fermi_energy: number | null;
  delta_e: number | null;
  fil_fermi: string;
  primary_file: string;
  bxsf_files: BxsfFileData[];
}

interface FermiSurfaceTaskPlan {
  taskLabel: string;
  taskParams: Record<string, any>;
  saveParameters: Record<string, any>;
}

function Tooltip({ text }: { text: string }) {
  return (
    <span className="tooltip-container">
      <span className="tooltip-icon">?</span>
      <span className="tooltip-text">{text}</span>
    </span>
  );
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

function getCalculationTags(calc: CalculationRun): { label: string; type: "info" | "feature" | "special" | "geometry" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" | "geometry" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature" | "special" | "geometry") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }
  if (params.conv_thr) {
    const thr = params.conv_thr;
    pushTag(thr < 0.001 ? thr.toExponential(0) : thr.toString(), "info");
  }
  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }
  if (params.nspin === 4) {
    pushTag("Non-collinear", "feature");
  } else if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }
  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  return tags;
}

function getBaseElement(symbol: string): string {
  return symbol.replace(/[\d+-]+$/, "");
}

function pseudoMatchesElement(filename: string, element: string): boolean {
  const lowerFile = filename.toLowerCase();
  const lowerEl = getBaseElement(element).toLowerCase();
  return (
    lowerFile.startsWith(`${lowerEl}.`)
    || lowerFile.startsWith(`${lowerEl}_`)
    || lowerFile.startsWith(`${lowerEl}-`)
  );
}

function parseOptionalNumber(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (trimmed.length === 0) return null;
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed)) {
    throw new Error(`Invalid ${label}.`);
  }
  return parsed;
}

function parseOptionalPositiveNumber(input: string, label: string): number | null {
  const parsed = parseOptionalNumber(input, label);
  if (parsed == null) return null;
  if (parsed <= 0) {
    throw new Error(`${label} must be positive.`);
  }
  return parsed;
}

function parseOptionalPositiveInt(input: string, label: string): number | null {
  const parsed = parseOptionalNumber(input, label);
  if (parsed == null) return null;
  if (!Number.isInteger(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive integer.`);
  }
  return parsed;
}

function sanitizeOutputStem(raw: string): string {
  const trimmed = raw.trim();
  if (trimmed.length === 0) return "fermi_surface";
  const sanitized = trimmed
    .replace(/[^a-zA-Z0-9_-]/g, "_")
    .replace(/^_+|_+$/g, "");
  return sanitized.length > 0 ? sanitized : "fermi_surface";
}

function normalizeOccupations(raw: unknown): "fixed" | "smearing" | "from_input" | "tetrahedra" {
  const lowered = String(raw || "smearing").toLowerCase();
  if (lowered === "fixed") return "fixed";
  if (lowered === "from_input") return "from_input";
  if (lowered === "tetrahedra") return "tetrahedra";
  return "smearing";
}

function normalizeSmearing(raw: unknown): "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac" {
  const lowered = String(raw || "gaussian").toLowerCase();
  if (lowered === "methfessel-paxton") return "methfessel-paxton";
  if (lowered === "marzari-vanderbilt") return "marzari-vanderbilt";
  if (lowered === "fermi-dirac") return "fermi-dirac";
  return "gaussian";
}

export function FermiSurfaceWizard({
  onBack,
  qePath,
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: FermiSurfaceWizardProps) {
  const taskContext = useTaskContext();
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );

  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);

  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());
  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);

  const [fermiKGrid, setFermiKGrid] = useState<[number, number, number]>([24, 24, 24]);
  const [fermiKOffset, setFermiKOffset] = useState<[0 | 1, 0 | 1, 0 | 1]>([0, 0, 0]);
  const [fermiNbnd, setFermiNbnd] = useState<number | "auto">("auto");
  const [nscfConvThrInput, setNscfConvThrInput] = useState("1e-8");
  const [nscfMixingBetaInput, setNscfMixingBetaInput] = useState("0.7");
  const [nscfOccupations, setNscfOccupations] = useState<"fixed" | "smearing" | "from_input" | "tetrahedra">("smearing");
  const [nscfSmearing, setNscfSmearing] = useState<"gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac">("gaussian");
  const [nscfDegaussInput, setNscfDegaussInput] = useState("0.02");
  const [nscfVerbosity, setNscfVerbosity] = useState<"low" | "high" | "debug">("high");
  const [deltaEInput, setDeltaEInput] = useState("0.02");
  const [outputStemInput, setOutputStemInput] = useState("fermi_surface");
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({
    mesh: true,
    nscf: true,
    fsx: true,
    mpi: false,
  });

  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Fermi Surface",
  });
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [fermiData, setFermiData] = useState<FermiSurfaceData | null>(null);
  const [scfFermiEnergy, setScfFermiEnergy] = useState<number | null>(null);
  const [showOutput, setShowOutput] = useState(false);
  const [isSaved, setIsSaved] = useState(false);

  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);

  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});

  const handleOutputScroll = () => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followOutputRef.current = distanceToBottom <= 16;
  };

  const toggleSection = (section: string) => {
    setExpandedSections((prev) => ({ ...prev, [section]: !prev[section] }));
  };

  const applyScfDefaults = useCallback((scf: CalculationRun) => {
    const params = scf.parameters || {};

    const convThr = Number(params.conv_thr);
    if (Number.isFinite(convThr) && convThr > 0) {
      setNscfConvThrInput(String(convThr));
    }

    const mixingBeta = Number(params.mixing_beta);
    if (Number.isFinite(mixingBeta) && mixingBeta > 0) {
      setNscfMixingBetaInput(String(mixingBeta));
    }

    setNscfOccupations(normalizeOccupations(params.occupations));
    setNscfSmearing(normalizeSmearing(params.smearing));

    const degauss = Number(params.degauss);
    if (Number.isFinite(degauss) && degauss > 0) {
      setNscfDegaussInput(String(degauss));
    } else {
      setNscfDegaussInput("");
    }

    const sourceNbnd = Number(params.nbnd);
    if (Number.isInteger(sourceNbnd) && sourceNbnd > 0) {
      setFermiNbnd(sourceNbnd);
    } else {
      setFermiNbnd("auto");
    }

    const verbosityRaw = String(params.verbosity || "high").toLowerCase();
    if (verbosityRaw === "low" || verbosityRaw === "debug") {
      setNscfVerbosity(verbosityRaw);
    } else {
      setNscfVerbosity("high");
    }
  }, []);

  const selectSourceScf = useCallback((scf: CalculationRun) => {
    setSelectedScf(scf);
    setScfFermiEnergy(scf.result?.fermi_energy ?? null);
    applyScfDefaults(scf);
  }, [applyScfDefaults]);

  useEffect(() => {
    async function init() {
      try {
        const count = await invoke<number>("get_cpu_count");
        setCpuCount(count);
        setMpiProcs(Math.max(1, Math.floor(count * 0.75)));

        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);

        const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
        const pseudos = await invoke<string[]>("list_pseudopotentials", { pseudoDir });
        const elements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
        const selected: Record<string, string> = {};
        for (const el of elements) {
          const match = pseudos.find((pseudo) => pseudoMatchesElement(pseudo, el));
          if (match) selected[el] = match;
        }
        setSelectedPseudos(selected);
      } catch (e) {
        console.error("Failed to initialize Fermi-surface wizard:", e);
      }
    }

    void init();
  }, [crystalData, qePath]);

  useEffect(() => {
    const el = outputRef.current;
    if (!el || !followOutputRef.current) return;
    el.scrollTop = el.scrollHeight;
  }, [output]);

  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task) {
      void taskContext.reconnectToTask(activeTaskId);
      return;
    }

    setIsRunning(task.status === "running");
    setOutput(task.output.join("\n") + (task.output.length > 0 ? "\n" : ""));
    setProgress(task.progress);
    setCalcStartTime(task.startedAt);

    if (task.status === "completed" && task.result) {
      setFermiData(task.result as FermiSurfaceData);
      setStep("results");
    } else if (task.status === "failed" || task.status === "cancelled") {
      setError(task.error || "Task failed");
    } else {
      setStep("run");
    }
  }, [activeTaskId, taskContext, taskContext.getTask(activeTaskId ?? "")?.status]);

  useEffect(() => {
    if (!activeTaskId) return;
    const task = taskContext.getTask(activeTaskId);
    if (!task || task.status !== "running") return;

    setOutput(task.output.join("\n") + (task.output.length > 0 ? "\n" : ""));
    setProgress(task.progress);
  }, [activeTaskId, taskContext, taskContext.getTask(activeTaskId ?? "")?.output.length]);

  function buildTaskPlan(): FermiSurfaceTaskPlan {
    if (!selectedScf) {
      throw new Error("No source SCF calculation selected");
    }

    const parsedDeltaE = parseOptionalPositiveNumber(deltaEInput, "DeltaE");
    const parsedConvThr = parseOptionalPositiveNumber(nscfConvThrInput, "NSCF convergence threshold");
    if (parsedConvThr == null) {
      throw new Error("Please provide a valid positive NSCF convergence threshold.");
    }

    const parsedMixingBeta = parseOptionalPositiveNumber(nscfMixingBetaInput, "NSCF mixing beta");
    if (parsedMixingBeta == null || parsedMixingBeta > 1.0) {
      throw new Error("NSCF mixing beta must be in the range (0, 1].");
    }

    const parsedDegauss = nscfOccupations === "smearing"
      ? parseOptionalPositiveNumber(nscfDegaussInput, "NSCF degauss")
      : null;
    if (nscfOccupations === "smearing" && parsedDegauss == null) {
      throw new Error("Please provide a positive degauss value when using smearing occupations.");
    }

    const parsedNbnd = fermiNbnd === "auto"
      ? null
      : parseOptionalPositiveInt(String(fermiNbnd), "number of bands");

    const outputStem = sanitizeOutputStem(outputStemInput);

    const scfParams = selectedScf.parameters || {};
    const savedPseudoMap = (scfParams.selected_pseudos && typeof scfParams.selected_pseudos === "object")
      ? scfParams.selected_pseudos as Record<string, string>
      : {};
    const elements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
    const resolvedPseudos: Record<string, string> = {};
    for (const element of elements) {
      const savedPseudo = savedPseudoMap[element];
      const detectedPseudo = selectedPseudos[element];
      const resolvedPseudo = (typeof savedPseudo === "string" && savedPseudo.length > 0)
        ? savedPseudo
        : detectedPseudo;
      if (!resolvedPseudo) {
        throw new Error(`No pseudopotential selected for element ${element}`);
      }
      resolvedPseudos[element] = resolvedPseudo;
    }

    const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
    const prefix = scfParams.prefix || "qcortado_scf";

    const ecutwfcValue = Number(scfParams.ecutwfc);
    const ecutwfc = Number.isFinite(ecutwfcValue) && ecutwfcValue > 0 ? ecutwfcValue : 40;
    const ecutrhoValue = Number(scfParams.ecutrho);
    const ecutrho = Number.isFinite(ecutrhoValue) && ecutrhoValue > 0
      ? ecutrhoValue
      : ecutwfc * 8;

    const primitiveCell: PrimitiveCell | null = getPrimitiveCell(crystalData);
    const commonSystemFields = {
      species: elements.map((element) => ({
        symbol: element,
        mass: ELEMENT_MASSES[element] || 1.0,
        pseudopotential: resolvedPseudos[element],
      })),
      position_units: "crystal",
      ecutwfc,
      ecutrho,
      nspin: Number(scfParams.nspin) || 1,
      occupations: nscfOccupations,
      smearing: nscfSmearing,
      degauss: parsedDegauss,
    };

    const systemConfig = primitiveCell
      ? {
        ...commonSystemFields,
        ibrav: primitiveCell.ibrav,
        celldm: [primitiveCell.celldm1, 0, 0, 0, 0, 0],
        cell_parameters: null,
        cell_units: null,
        atoms: primitiveCell.atoms.map((atom) => ({
          symbol: atom.symbol,
          position: atom.position,
          if_pos: [true, true, true] as [boolean, boolean, boolean],
        })),
      }
      : (() => {
        const gammaRadians = (crystalData.cell_angle_gamma.value * Math.PI) / 180;
        const betaRadians = (crystalData.cell_angle_beta.value * Math.PI) / 180;
        const alphaRadians = (crystalData.cell_angle_alpha.value * Math.PI) / 180;
        const c = crystalData.cell_length_c.value;
        const cVector: [number, number, number] = [
          c * Math.cos(betaRadians),
          (c * (Math.cos(alphaRadians) - Math.cos(betaRadians) * Math.cos(gammaRadians))) / Math.sin(gammaRadians),
          Math.sqrt(
            Math.max(
              0,
              c * c
                - (c * Math.cos(betaRadians)) ** 2
                - (
                  (c * (Math.cos(alphaRadians) - Math.cos(betaRadians) * Math.cos(gammaRadians)))
                  / Math.sin(gammaRadians)
                ) ** 2,
            ),
          ),
        ];

        return {
          ...commonSystemFields,
          ibrav: "free" as const,
          celldm: null,
          cell_parameters: [
            [crystalData.cell_length_a.value, 0, 0],
            [
              crystalData.cell_length_b.value * Math.cos(gammaRadians),
              crystalData.cell_length_b.value * Math.sin(gammaRadians),
              0,
            ],
            cVector,
          ] as [[number, number, number], [number, number, number], [number, number, number]],
          cell_units: "angstrom" as const,
          atoms: crystalData.atom_sites.map((site) => ({
            symbol: getBaseElement(site.type_symbol),
            position: [site.fract_x, site.fract_y, site.fract_z] as [number, number, number],
            if_pos: [true, true, true] as [boolean, boolean, boolean],
          })),
        };
      })();

    const baseCalculation = {
      calculation: "scf",
      prefix,
      outdir: "./tmp",
      pseudo_dir: pseudoDir,
      system: systemConfig,
      kpoints: { type: "gamma" },
      conv_thr: parsedConvThr,
      mixing_beta: parsedMixingBeta,
      tprnfor: false,
      tstress: false,
      forc_conv_thr: null,
      etot_conv_thr: null,
      verbosity: nscfVerbosity,
    };

    const taskLabel = `Fermi Surface - ${crystalData?.formula_sum || ""}`;
    const taskParams = {
      config: {
        base_calculation: baseCalculation,
        k_grid: fermiKGrid,
        k_offset: fermiKOffset,
        nbnd: parsedNbnd,
        delta_e: parsedDeltaE,
        fil_fermi: outputStem,
        project_id: projectId,
        scf_calc_id: selectedScf.id,
      },
      workingDir: FERMI_WORK_DIR,
      mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
    };

    const saveParameters = {
      source_scf_id: selectedScf.id,
      fermi_k_grid: fermiKGrid,
      fermi_k_offset: fermiKOffset,
      fermi_n_bands_requested: parsedNbnd,
      fermi_delta_e: parsedDeltaE,
      fermi_output_stem: outputStem,
      nscf_conv_thr: parsedConvThr,
      nscf_mixing_beta: parsedMixingBeta,
      nscf_occupations: nscfOccupations,
      nscf_smearing: nscfSmearing,
      nscf_degauss: parsedDegauss,
      nscf_verbosity: nscfVerbosity,
      n_bxsf_files: null,
      bxsf_files: [],
      primary_bxsf_file: null,
      total_bxsf_bytes: null,
      nspin: scfParams.nspin,
      lspinorb: scfParams.lspinorb,
      lda_plus_u: scfParams.lda_plus_u,
      vdw_corr: scfParams.vdw_corr,
      scf_fermi_energy: selectedScf.result?.fermi_energy ?? scfFermiEnergy,
    };

    return {
      taskLabel,
      taskParams,
      saveParameters,
    };
  }

  async function runCalculation() {
    if (hasExternalRunningTask) {
      setError("Another task is currently running. Queue this task or wait for completion.");
      return;
    }

    setIsRunning(true);
    followOutputRef.current = true;
    setOutput("");
    setError(null);
    setFermiData(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Fermi Surface"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

    try {
      const plan = buildTaskPlan();
      const taskId = await taskContext.startTask("fermi_surface", plan.taskParams, plan.taskLabel);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Calculation failed");
      }

      const result = finalTask.result as FermiSurfaceData;
      const outputContent = finalTask.output.join("\n");
      const endTime = new Date().toISOString();
      setFermiData(result);
      setStep("results");
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      try {
        const bxsfFiles = result.bxsf_files.map((file) => file.file_name);
        const totalBxsfBytes = result.bxsf_files.reduce((sum, file) => sum + (Number(file.size_bytes) || 0), 0);

        await invoke("save_calculation", {
          projectId,
          cifId,
          calcData: {
            calc_type: "fermi_surface",
            parameters: {
              ...plan.saveParameters,
              fermi_output_stem: result.fil_fermi,
              n_bxsf_files: result.bxsf_files.length,
              bxsf_files: bxsfFiles,
              primary_bxsf_file: result.primary_file,
              total_bxsf_bytes: totalBxsfBytes,
            },
            result: {
              converged: true,
              total_energy: null,
              fermi_energy: result.fermi_energy ?? selectedScf?.result?.fermi_energy ?? scfFermiEnergy,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: outputContent,
            },
            started_at: startTime,
            completed_at: endTime,
            input_content: "",
            output_content: outputContent,
            tags: [],
          },
          workingDir: FERMI_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save Fermi-surface calculation:", saveError);
        setError(`Failed to auto-save Fermi-surface calculation: ${saveError}`);
      }
    } catch (e) {
      setError(String(e));
      setOutput((prev) => prev + `\nError: ${e}\n`);
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

  function queueCalculation() {
    try {
      const plan = buildTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "fermi_surface",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId,
          workingDir: FERMI_WORK_DIR,
          calcType: "fermi_surface",
          parameters: plan.saveParameters,
          tags: [],
          inputContent: "",
        },
      );
    } catch (e) {
      setError(String(e));
    }
  }

  const renderSourceStep = () => {
    const validScfs = scfCalculations.filter((calc) => calc.calc_type === "scf" && calc.result?.converged);
    const sortedScfs = sortScfByMode(validScfs, scfSortMode);

    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Fermi-surface generation requires a completed SCF calculation.
            Please run an SCF calculation first, then return here.
          </p>
          <button className="primary-button" onClick={onBack}>
            Back to Dashboard
          </button>
        </div>
      );
    }

    return (
      <div className="wizard-step source-step">
        <div className="source-step-header">
          <h3>Select Source SCF Calculation</h3>
          <div className="source-sort-control">
            <label htmlFor="fermi-scf-sort">Sort SCFs</label>
            <select
              id="fermi-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose the SCF calculation used to build a dense NSCF mesh for `fs.x`.
        </p>

        <div className="scf-list">
          {sortedScfs.map((scf) => (
            <div
              key={scf.id}
              className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""}`}
              onClick={() => selectSourceScf(scf)}
            >
              <div className="scf-option-header">
                <input
                  type="radio"
                  checked={selectedScf?.id === scf.id}
                  onChange={() => selectSourceScf(scf)}
                />
                <span className="scf-date">
                  {new Date(scf.started_at).toLocaleDateString()}
                </span>
              </div>
              <div className="scf-details">
                <span>E = {scf.result?.total_energy?.toFixed(6)} Ry</span>
                {scf.result?.fermi_energy != null && (
                  <span>EF = {scf.result.fermi_energy.toFixed(3)} eV</span>
                )}
              </div>
              <div className="calc-tags">
                {getCalculationTags(scf).map((tag, i) => (
                  <span key={`${tag.label}-${i}`} className={`calc-tag calc-tag-${tag.type}`}>
                    {tag.label}
                  </span>
                ))}
              </div>
            </div>
          ))}
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Cancel
          </button>
          <button
            className="primary-button"
            disabled={!selectedScf}
            onClick={() => setStep("parameters")}
          >
            Next: Parameters
          </button>
        </div>
      </div>
    );
  };

  const renderParametersStep = () => {
    const outputStem = sanitizeOutputStem(outputStemInput);
    const safeParsePositive = (value: string): number | null => {
      const trimmed = value.trim();
      if (!trimmed) return null;
      const parsed = Number(trimmed);
      if (!Number.isFinite(parsed) || parsed <= 0) return null;
      return parsed;
    };
    const parsedConvThr = safeParsePositive(nscfConvThrInput);
    const parsedMixingBeta = safeParsePositive(nscfMixingBetaInput);
    const parsedDegauss = safeParsePositive(nscfDegaussInput);
    const parsedDeltaE = safeParsePositive(deltaEInput);
    const degaussRequired = nscfOccupations === "smearing";
    const isNbndValid = fermiNbnd === "auto" || (Number.isInteger(fermiNbnd) && fermiNbnd > 0);
    const isConvThrValid = parsedConvThr !== null;
    const isMixingBetaValid = parsedMixingBeta !== null && parsedMixingBeta <= 1.0;
    const isDegaussValid = !degaussRequired || parsedDegauss !== null;
    const isDeltaEValid = deltaEInput.trim().length === 0 || parsedDeltaE !== null;
    const canRun = isNbndValid && isConvThrValid && isMixingBetaValid && isDegaussValid && isDeltaEValid;

    return (
      <div className="wizard-step parameters-step">
        <h3>Fermi-Surface Parameters</h3>
        <p className="step-description">
          Configure dense NSCF controls and `fs.x` options for Fermi-surface BXSF generation.
        </p>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("mesh")} className="section-header">
            <span className={`collapse-icon ${expandedSections.mesh ? "expanded" : ""}`}>▶</span>
            NSCF Mesh & Band Count
            <Tooltip text="Controls dense NSCF sampling mesh (`K_POINTS automatic`) and optional manual `nbnd` override." />
          </h4>
          {expandedSections.mesh && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    Dense k-point grid
                    <Tooltip text="Finer NSCF meshes improve Fermi-surface topology fidelity but increase runtime and memory usage significantly." />
                  </label>
                  <div className="qgrid-inputs">
                    <input
                      type="number"
                      min={1}
                      value={fermiKGrid[0]}
                      onChange={(e) => setFermiKGrid([Math.max(1, parseInt(e.target.value, 10) || 1), fermiKGrid[1], fermiKGrid[2]])}
                    />
                    <span>×</span>
                    <input
                      type="number"
                      min={1}
                      value={fermiKGrid[1]}
                      onChange={(e) => setFermiKGrid([fermiKGrid[0], Math.max(1, parseInt(e.target.value, 10) || 1), fermiKGrid[2]])}
                    />
                    <span>×</span>
                    <input
                      type="number"
                      min={1}
                      value={fermiKGrid[2]}
                      onChange={(e) => setFermiKGrid([fermiKGrid[0], fermiKGrid[1], Math.max(1, parseInt(e.target.value, 10) || 1)])}
                    />
                  </div>
                </div>

                <div className="phonon-field">
                  <label>
                    K-point offset (0/1)
                    <Tooltip text="Monkhorst-Pack offset per reciprocal direction. Use 1 to shift half a grid step in that direction." />
                  </label>
                  <div className="qgrid-inputs">
                    <select
                      value={fermiKOffset[0]}
                      onChange={(e) => setFermiKOffset([Number(e.target.value) as 0 | 1, fermiKOffset[1], fermiKOffset[2]])}
                    >
                      <option value={0}>0</option>
                      <option value={1}>1</option>
                    </select>
                    <span>×</span>
                    <select
                      value={fermiKOffset[1]}
                      onChange={(e) => setFermiKOffset([fermiKOffset[0], Number(e.target.value) as 0 | 1, fermiKOffset[2]])}
                    >
                      <option value={0}>0</option>
                      <option value={1}>1</option>
                    </select>
                    <span>×</span>
                    <select
                      value={fermiKOffset[2]}
                      onChange={(e) => setFermiKOffset([fermiKOffset[0], fermiKOffset[1], Number(e.target.value) as 0 | 1])}
                    >
                      <option value={0}>0</option>
                      <option value={1}>1</option>
                    </select>
                  </div>
                </div>

                <div className="phonon-field">
                  <label>
                    Number of bands
                    <Tooltip text="QE variable: `nbnd` in NSCF `pw.x`. Increase to include enough unoccupied states near Fermi crossings." />
                  </label>
                  <div className="nbnd-input">
                    <select
                      value={fermiNbnd === "auto" ? "auto" : "manual"}
                      onChange={(e) => setFermiNbnd(e.target.value === "auto" ? "auto" : 20)}
                    >
                      <option value="auto">Auto (QE default)</option>
                      <option value="manual">Manual</option>
                    </select>
                    {fermiNbnd !== "auto" && (
                      <input
                        type="number"
                        min={1}
                        value={fermiNbnd}
                        onChange={(e) => setFermiNbnd(Math.max(1, parseInt(e.target.value, 10) || 1))}
                      />
                    )}
                  </div>
                  {!isNbndValid && (
                    <span className="param-hint input-error">Use a positive integer for manual `nbnd`.</span>
                  )}
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("nscf")} className="section-header">
            <span className={`collapse-icon ${expandedSections.nscf ? "expanded" : ""}`}>▶</span>
            NSCF Electronic Controls
            <Tooltip text="Advanced NSCF settings: `conv_thr`, `mixing_beta`, occupations/smearing/degauss, and output verbosity." />
          </h4>
          {expandedSections.nscf && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    NSCF convergence threshold
                    <Tooltip text="QE variable: `conv_thr` in `&ELECTRONS`." />
                  </label>
                  <input
                    type="text"
                    value={nscfConvThrInput}
                    onChange={(e) => setNscfConvThrInput(e.target.value)}
                    placeholder="1e-8"
                    spellCheck={false}
                  />
                  {!isConvThrValid && (
                    <span className="param-hint input-error">Use a positive number (e.g. `1e-8`).</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Mixing beta
                    <Tooltip text="QE variable: `mixing_beta`. Stable values are typically in (0, 1]." />
                  </label>
                  <input
                    type="text"
                    value={nscfMixingBetaInput}
                    onChange={(e) => setNscfMixingBetaInput(e.target.value)}
                    placeholder="0.7"
                    spellCheck={false}
                  />
                  {!isMixingBetaValid && (
                    <span className="param-hint input-error">Use a positive number in the range (0, 1].</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Occupations
                    <Tooltip text="QE variable: `occupations` for NSCF. Metals typically require smearing." />
                  </label>
                  <select
                    value={nscfOccupations}
                    onChange={(e) => setNscfOccupations(e.target.value as "fixed" | "smearing" | "from_input" | "tetrahedra")}
                  >
                    <option value="smearing">smearing</option>
                    <option value="fixed">fixed</option>
                    <option value="tetrahedra">tetrahedra</option>
                    <option value="from_input">from_input</option>
                  </select>
                </div>

                <div className="phonon-field">
                  <label>
                    Smearing type
                    <Tooltip text="QE variable: `smearing`. Used only when occupations = smearing." />
                  </label>
                  <select
                    value={nscfSmearing}
                    onChange={(e) => setNscfSmearing(e.target.value as "gaussian" | "methfessel-paxton" | "marzari-vanderbilt" | "fermi-dirac")}
                    disabled={nscfOccupations !== "smearing"}
                  >
                    <option value="gaussian">gaussian</option>
                    <option value="methfessel-paxton">methfessel-paxton</option>
                    <option value="marzari-vanderbilt">marzari-vanderbilt</option>
                    <option value="fermi-dirac">fermi-dirac</option>
                  </select>
                </div>

                <div className="phonon-field">
                  <label>
                    Degauss (Ry)
                    <Tooltip text="QE variable: `degauss`. Required with smearing occupations." />
                  </label>
                  <input
                    type="text"
                    value={nscfDegaussInput}
                    onChange={(e) => setNscfDegaussInput(e.target.value)}
                    placeholder={nscfOccupations === "smearing" ? "0.02" : "unused"}
                    spellCheck={false}
                    disabled={nscfOccupations !== "smearing"}
                  />
                  {!isDegaussValid && (
                    <span className="param-hint input-error">Provide a positive degauss when occupations = smearing.</span>
                  )}
                </div>

                <div className="phonon-field">
                  <label>
                    Output detail level
                    <Tooltip text="QE variable: `verbosity` in `&CONTROL`." />
                  </label>
                  <select
                    value={nscfVerbosity}
                    onChange={(e) => setNscfVerbosity(e.target.value as "low" | "high" | "debug")}
                  >
                    <option value="low">low</option>
                    <option value="high">high</option>
                    <option value="debug">debug</option>
                  </select>
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="option-section config-section collapsible">
          <h4 onClick={() => toggleSection("fsx")} className="section-header">
            <span className={`collapse-icon ${expandedSections.fsx ? "expanded" : ""}`}>▶</span>
            fs.x Output Controls
            <Tooltip text="Configure `fs.x` inputs (`filfermi`, `DeltaE`) used to generate BXSF output." />
          </h4>
          {expandedSections.fsx && (
            <div className="option-params">
              <div className="phonon-grid">
                <div className="phonon-field">
                  <label>
                    Output stem
                    <Tooltip text="QE variable: `filfermi` in `&fermi`. BXSF output is typically `<stem>_fs.bxsf`." />
                  </label>
                  <input
                    type="text"
                    value={outputStemInput}
                    onChange={(e) => setOutputStemInput(e.target.value)}
                    placeholder="fermi_surface"
                  />
                  <span className="param-hint">Resolved stem: `{outputStem}`</span>
                </div>

                <div className="phonon-field">
                  <label>
                    DeltaE (eV, optional)
                    <Tooltip text="QE variable: `DeltaE` in `fs.x`. Leave blank to use fs.x defaults." />
                  </label>
                  <input
                    type="text"
                    value={deltaEInput}
                    onChange={(e) => setDeltaEInput(e.target.value)}
                    placeholder="0.02"
                    spellCheck={false}
                  />
                  {!isDeltaEValid && (
                    <span className="param-hint input-error">Use a positive number, or leave blank for default behavior.</span>
                  )}
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>Source SCF:</span>
            <span>{selectedScf?.id.slice(0, 8) || "N/A"}</span>
          </div>
          <div className="summary-row">
            <span>Dense k-grid:</span>
            <span>{fermiKGrid[0]}×{fermiKGrid[1]}×{fermiKGrid[2]}</span>
          </div>
          <div className="summary-row">
            <span>K-offset:</span>
            <span>{fermiKOffset[0]} {fermiKOffset[1]} {fermiKOffset[2]}</span>
          </div>
          <div className="summary-row">
            <span>Requested `nbnd`:</span>
            <span>{fermiNbnd === "auto" ? "Auto" : fermiNbnd}</span>
          </div>
          <div className="summary-row">
            <span>Output stem:</span>
            <span>{outputStem}</span>
          </div>
          <div className="summary-row">
            <span>Expected BXSF:</span>
            <span>{outputStem}_fs.bxsf</span>
          </div>
          <div className="summary-row">
            <span>DeltaE:</span>
            <span>{deltaEInput.trim().length > 0 ? `${deltaEInput.trim()} eV` : "Default"}</span>
          </div>
        </div>

        <div className="option-section mpi-section config-section collapsible">
          <h4 onClick={() => toggleSection("mpi")} className="section-header">
            <span className={`collapse-icon ${expandedSections.mpi ? "expanded" : ""}`}>▶</span>
            Parallelization
            <Tooltip text="Process-level MPI controls for the NSCF stage." />
          </h4>
          {expandedSections.mpi && (
            <div className="option-params">
              {mpiAvailable ? (
                <div className="mpi-toggle">
                  <label>
                    <input
                      type="checkbox"
                      checked={mpiEnabled}
                      onChange={(e) => setMpiEnabled(e.target.checked)}
                    />
                    Enable MPI for NSCF ({cpuCount} cores available)
                  </label>
                  {mpiEnabled && (
                    <div className="mpi-procs">
                      <label>
                        Number of processes:
                        <input
                          type="number"
                          min={1}
                          max={cpuCount}
                          value={mpiProcs}
                          onChange={(e) => setMpiProcs(Math.max(1, parseInt(e.target.value, 10) || 1))}
                        />
                      </label>
                    </div>
                  )}
                </div>
              ) : (
                <p className="mpi-unavailable">
                  MPI not available. Running in serial mode.
                </p>
              )}
            </div>
          )}
        </div>

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button className="secondary-button" onClick={() => void queueCalculation()} disabled={!canRun}>
            Queue Task
          </button>
          <button className="primary-button" onClick={() => void runCalculation()} disabled={!canRun || hasExternalRunningTask}>
            Run Calculation
          </button>
        </div>
      </div>
    );
  };

  const renderRunStep = () => (
    <div className="wizard-step run-step">
      <h3>{isRunning ? "Running Fermi-Surface Generation" : "Fermi-Surface Output"}</h3>

      <ProgressBar
        status={progress.status}
        percent={progress.percent}
        phase={progress.phase}
        detail={progress.detail}
      />
      <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />

      <pre ref={outputRef} className="calculation-output" onScroll={handleOutputScroll}>
        {output || "Starting calculation..."}
      </pre>

      {error && (
        <div className="error-actions">
          <button className="secondary-button" onClick={() => setStep("parameters")}>
            Back to Parameters
          </button>
        </div>
      )}
    </div>
  );

  const renderResultsStep = () => {
    if (!fermiData) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("parameters")}>
            Try Again
          </button>
        </div>
      );
    }

    const totalBytes = fermiData.bxsf_files.reduce((sum, file) => sum + (Number(file.size_bytes) || 0), 0);

    return (
      <div className="wizard-step results-step">
        <h3>Fermi Surface Results</h3>
        <p className="step-description">
          BXSF generation complete. Viewer launch wiring is the next step.
        </p>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Dense k-grid:</span>
              <span className="value">{fermiData.k_grid[0]}×{fermiData.k_grid[1]}×{fermiData.k_grid[2]}</span>
            </div>
            <div className="summary-item">
              <span className="label">BXSF Files:</span>
              <span className="value">{fermiData.bxsf_files.length}</span>
            </div>
            <div className="summary-item">
              <span className="label">Primary File:</span>
              <span className="value">{fermiData.primary_file}</span>
            </div>
            <div className="summary-item">
              <span className="label">Total Size:</span>
              <span className="value">{formatBytes(totalBytes)}</span>
            </div>
            <div className="summary-item">
              <span className="label">Fermi Energy:</span>
              <span className="value">
                {fermiData.fermi_energy != null ? `${fermiData.fermi_energy.toFixed(3)} eV` : "Unavailable"}
              </span>
            </div>
          </div>
        </div>

        <div className="detail-item parameters">
          <label>Generated BXSF Metadata</label>
          <pre>{JSON.stringify(fermiData.bxsf_files, null, 2)}</pre>
        </div>

        <div className="output-section">
          <button
            className="output-toggle"
            onClick={() => setShowOutput(!showOutput)}
          >
            {showOutput ? "▼" : "▶"} Calculation Output
          </button>
          {showOutput && (
            <pre className="calculation-output results-output">
              {output || "No output available"}
            </pre>
          )}
        </div>

        {isSaved && (
          <div className="save-status save-status-inline">
            <span className="saved">Saved to project</span>
          </div>
        )}

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Back to Dashboard
          </button>
          <button className="primary-button" onClick={() => setStep("parameters")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className="electronic-dos-wizard fermi-surface-wizard">
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Fermi Surface Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "parameters" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            2. Parameters
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "completed" : ""}>
            3. Run
          </span>
          <span className={step === "results" ? "active" : ""}>
            4. Results
          </span>
        </div>
      </div>

      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "parameters" && renderParametersStep()}
        {step === "run" && renderRunStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
