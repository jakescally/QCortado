import { useState, useEffect, useRef, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData, ELEMENT_MASSES } from "../lib/types";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { defaultProgressState, ProgressState } from "../lib/qeProgress";
import { useTaskContext } from "../lib/TaskContext";
import { ElectronicDOSData } from "./ElectronicDOSPlot";

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

interface ElectronicDOSWizardProps {
  onBack: () => void;
  onViewDos: (dosData: ElectronicDOSData, scfFermiEnergy: number | null) => void;
  qePath: string;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "parameters" | "run" | "results";
const DOS_WORK_DIR = "/tmp/qcortado_dos";

interface DosTaskPlan {
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

export function ElectronicDOSWizard({
  onBack,
  onViewDos,
  qePath,
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: ElectronicDOSWizardProps) {
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

  const [dosKGrid, setDosKGrid] = useState<[number, number, number]>([16, 16, 16]);
  const [deltaEInput, setDeltaEInput] = useState("0.02");
  const [degaussInput, setDegaussInput] = useState("0.02");
  const [eminInput, setEminInput] = useState("");
  const [emaxInput, setEmaxInput] = useState("");

  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Electronic DOS",
  });
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [dosData, setDosData] = useState<ElectronicDOSData | null>(null);
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
        console.error("Failed to initialize DOS wizard:", e);
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
      setDosData(task.result as ElectronicDOSData);
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

  function buildDosTaskPlan(): DosTaskPlan {
    if (!selectedScf) {
      throw new Error("No source SCF calculation selected");
    }

    const parsedDeltaE = parseOptionalPositiveNumber(deltaEInput, "DOS energy step");
    const parsedDegauss = parseOptionalPositiveNumber(degaussInput, "DOS degauss");
    const parsedEmin = parseOptionalNumber(eminInput, "minimum energy");
    const parsedEmax = parseOptionalNumber(emaxInput, "maximum energy");

    if (parsedEmin != null && parsedEmax != null && parsedEmin >= parsedEmax) {
      throw new Error("Energy window is invalid. Emin must be smaller than Emax.");
    }

    const scfParams = selectedScf.parameters || {};
    const elements = [...new Set(crystalData.atom_sites.map((site) => getBaseElement(site.type_symbol)))];
    for (const element of elements) {
      if (!selectedPseudos[element]) {
        throw new Error(`No pseudopotential selected for element ${element}`);
      }
    }

    const occupationRaw = String(scfParams.occupations || "smearing").toLowerCase();
    const occupations = occupationRaw === "fixed"
      ? "fixed"
      : occupationRaw === "tetrahedra"
        ? "tetrahedra"
        : occupationRaw === "from_input"
          ? "from_input"
          : "smearing";
    const smearingRaw = String(scfParams.smearing || "gaussian").toLowerCase();
    const smearing = smearingRaw === "methfessel-paxton"
      || smearingRaw === "marzari-vanderbilt"
      || smearingRaw === "fermi-dirac"
      ? smearingRaw
      : "gaussian";
    const inheritedDegaussValue = Number(scfParams.degauss);
    const inheritedDegauss = Number.isFinite(inheritedDegaussValue) ? inheritedDegaussValue : null;
    const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
    const prefix = scfParams.prefix || "qcortado_scf";

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

    const baseCalculation = {
      calculation: "scf",
      prefix,
      outdir: "./tmp",
      pseudo_dir: pseudoDir,
      system: {
        ibrav: "free",
        celldm: null,
        cell_parameters: [
          [crystalData.cell_length_a.value, 0, 0],
          [
            crystalData.cell_length_b.value * Math.cos(gammaRadians),
            crystalData.cell_length_b.value * Math.sin(gammaRadians),
            0,
          ],
          cVector,
        ],
        cell_units: "angstrom",
        species: elements.map((element) => ({
          symbol: element,
          mass: ELEMENT_MASSES[element] || 1.0,
          pseudopotential: selectedPseudos[element],
        })),
        atoms: crystalData.atom_sites.map((site) => ({
          symbol: getBaseElement(site.type_symbol),
          position: [site.fract_x, site.fract_y, site.fract_z],
          if_pos: [true, true, true],
        })),
        position_units: "crystal",
        ecutwfc: Number(scfParams.ecutwfc) || 40,
        ecutrho: Number(scfParams.ecutrho) || 320,
        nspin: Number(scfParams.nspin) || 1,
        occupations,
        smearing,
        degauss: parsedDegauss ?? inheritedDegauss ?? 0.02,
      },
      kpoints: { type: "gamma" },
      conv_thr: Number(scfParams.conv_thr) || 1e-8,
      mixing_beta: Number(scfParams.mixing_beta) || 0.7,
      tprnfor: false,
      tstress: false,
      forc_conv_thr: null,
      etot_conv_thr: null,
      verbosity: "high",
    };

    const taskLabel = `DOS - ${crystalData?.formula_sum || ""}`;
    const taskParams = {
      config: {
        base_calculation: baseCalculation,
        k_grid: dosKGrid,
        degauss: parsedDegauss,
        emin: parsedEmin,
        emax: parsedEmax,
        delta_e: parsedDeltaE,
        project_id: projectId,
        scf_calc_id: selectedScf.id,
      },
      workingDir: DOS_WORK_DIR,
      mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
    };

    const saveParameters = {
      source_scf_id: selectedScf.id,
      dos_k_grid: dosKGrid,
      dos_delta_e: parsedDeltaE,
      dos_degauss: parsedDegauss,
      dos_emin: parsedEmin,
      dos_emax: parsedEmax,
      n_points: null,
      // Inherit SCF feature parameters for dashboard tags
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
    setDosData(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Electronic DOS"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

    try {
      const plan = buildDosTaskPlan();
      const taskId = await taskContext.startTask("dos", plan.taskParams, plan.taskLabel);
      setActiveTaskId(taskId);

      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask?.error || "Calculation failed");
      }

      const result = finalTask.result as ElectronicDOSData;
      const outputContent = finalTask.output.join("\n");
      const endTime = new Date().toISOString();
      setDosData(result);
      setStep("results");
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      try {
        await invoke("save_calculation", {
          projectId,
          cifId,
          calcData: {
            calc_type: "dos",
            parameters: {
              ...plan.saveParameters,
              n_points: result.points,
            },
            result: {
              converged: true,
              total_energy: null,
              fermi_energy: result.fermi_energy ?? selectedScf?.result?.fermi_energy ?? scfFermiEnergy,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: outputContent,
              dos_data: result,
            },
            started_at: startTime,
            completed_at: endTime,
            input_content: "",
            output_content: outputContent,
            tags: [],
          },
          workingDir: DOS_WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save DOS calculation:", saveError);
        setError(`Failed to auto-save DOS calculation: ${saveError}`);
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
      const plan = buildDosTaskPlan();
      setError(null);
      taskContext.enqueueTask(
        "dos",
        plan.taskParams,
        plan.taskLabel,
        {
          projectId,
          cifId,
          workingDir: DOS_WORK_DIR,
          calcType: "dos",
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
            Electronic DOS calculations require a completed SCF calculation.
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
            <label htmlFor="dos-scf-sort">Sort SCFs</label>
            <select
              id="dos-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose the SCF calculation used to build the dense NSCF sampling for DOS.
        </p>

        <div className="scf-list">
          {sortedScfs.map((scf) => (
            <div
              key={scf.id}
              className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""}`}
              onClick={() => {
                setSelectedScf(scf);
                setScfFermiEnergy(scf.result?.fermi_energy ?? null);
              }}
            >
              <div className="scf-option-header">
                <input
                  type="radio"
                  checked={selectedScf?.id === scf.id}
                  onChange={() => {
                    setSelectedScf(scf);
                    setScfFermiEnergy(scf.result?.fermi_energy ?? null);
                  }}
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
    return (
      <div className="wizard-step parameters-step">
        <h3>DOS Parameters</h3>
        <p className="step-description">
          Configure a dense NSCF k-grid and `dos.x` energy sampling.
        </p>

        <div className="param-section">
          <h4>
            NSCF Sampling Grid
            <Tooltip text="Dense k-point sampling controls DOS quality. Denser grids reduce numerical noise and improve peak positions, but increase NSCF runtime and memory use." />
          </h4>
          <div className="qgrid-inputs">
            <input
              type="number"
              min={1}
              value={dosKGrid[0]}
              onChange={(e) => setDosKGrid([Math.max(1, parseInt(e.target.value, 10) || 1), dosKGrid[1], dosKGrid[2]])}
            />
            <span>×</span>
            <input
              type="number"
              min={1}
              value={dosKGrid[1]}
              onChange={(e) => setDosKGrid([dosKGrid[0], Math.max(1, parseInt(e.target.value, 10) || 1), dosKGrid[2]])}
            />
            <span>×</span>
            <input
              type="number"
              min={1}
              value={dosKGrid[2]}
              onChange={(e) => setDosKGrid([dosKGrid[0], dosKGrid[1], Math.max(1, parseInt(e.target.value, 10) || 1)])}
            />
          </div>
        </div>

        <div className="param-section">
          <h4>DOS Energy Window</h4>
          <div className="param-grid">
            <div className="param-row">
              <label>
                DeltaE (eV)
                <Tooltip text="Energy bin width used by dos.x. Smaller DeltaE increases resolution of fine features, but also increases point count and can expose sampling noise." />
              </label>
              <input
                type="text"
                value={deltaEInput}
                onChange={(e) => setDeltaEInput(e.target.value)}
                placeholder="0.02"
              />
            </div>
            <div className="param-row">
              <label>
                Degauss (Ry)
                <Tooltip text="Broadening width for DOS integration. Larger values smooth the DOS; smaller values preserve sharper features but require denser k-point sampling." />
              </label>
              <input
                type="text"
                value={degaussInput}
                onChange={(e) => setDegaussInput(e.target.value)}
                placeholder="0.02"
              />
            </div>
            <div className="param-row">
              <label>
                Emin (eV, optional)
                <Tooltip text="Lower energy bound for DOS output. Set this to focus on valence-region features and reduce unnecessary output outside the range of interest." />
              </label>
              <input
                type="text"
                value={eminInput}
                onChange={(e) => setEminInput(e.target.value)}
                placeholder="Auto"
              />
            </div>
            <div className="param-row">
              <label>
                Emax (eV, optional)
                <Tooltip text="Upper energy bound for DOS output. Set this to focus on the relevant conduction-energy region instead of writing the full default range." />
              </label>
              <input
                type="text"
                value={emaxInput}
                onChange={(e) => setEmaxInput(e.target.value)}
                placeholder="Auto"
              />
            </div>
          </div>
        </div>

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>Source SCF:</span>
            <span>{selectedScf?.id.slice(0, 8) || "N/A"}</span>
          </div>
          <div className="summary-row">
            <span>DOS k-grid:</span>
            <span>{dosKGrid[0]}×{dosKGrid[1]}×{dosKGrid[2]}</span>
          </div>
          <div className="summary-row">
            <span>DeltaE:</span>
            <span>{deltaEInput || "Auto"} eV</span>
          </div>
          <div className="summary-row">
            <span>Energy window:</span>
            <span>{eminInput || "Auto"} to {emaxInput || "Auto"} eV</span>
          </div>
        </div>

        <div className="mpi-section">
          <h4>Parallelization</h4>
          {mpiAvailable ? (
            <div className="mpi-toggle">
              <label>
                <input
                  type="checkbox"
                  checked={mpiEnabled}
                  onChange={(e) => setMpiEnabled(e.target.checked)}
                />
                Enable MPI ({cpuCount} cores available)
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

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button className="secondary-button" onClick={() => void queueCalculation()}>
            Queue Task
          </button>
          <button className="primary-button" onClick={() => void runCalculation()} disabled={hasExternalRunningTask}>
            Run Calculation
          </button>
        </div>
      </div>
    );
  };

  const renderRunStep = () => (
    <div className="wizard-step run-step">
      <h3>{isRunning ? "Running Electronic DOS Calculation" : "Electronic DOS Output"}</h3>

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
    if (!dosData) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("parameters")}>
            Try Again
          </button>
        </div>
      );
    }

    return (
      <div className="wizard-step results-step">
        <h3>Electronic DOS Results</h3>
        <p className="step-description">
          DOS calculation complete. Open the DOS viewer for interactive plotting.
        </p>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Points:</span>
              <span className="value">{dosData.points}</span>
            </div>
            <div className="summary-item">
              <span className="label">Energy Min:</span>
              <span className="value">{dosData.energy_range[0].toFixed(3)} eV</span>
            </div>
            <div className="summary-item">
              <span className="label">Energy Max:</span>
              <span className="value">{dosData.energy_range[1].toFixed(3)} eV</span>
            </div>
            <div className="summary-item">
              <span className="label">Fermi Energy:</span>
              <span className="value">
                {dosData.fermi_energy != null ? `${dosData.fermi_energy.toFixed(3)} eV` : "Unavailable"}
              </span>
            </div>
          </div>
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
          <button
            className="primary-button"
            onClick={() => onViewDos(dosData, scfFermiEnergy)}
          >
            View DOS
          </button>
          <button className="primary-button" onClick={() => setStep("parameters")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className="electronic-dos-wizard">
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Electronic DOS Wizard</h2>
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
