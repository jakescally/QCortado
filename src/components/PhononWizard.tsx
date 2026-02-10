// Phonon Calculation Wizard - Calculate phonon DOS and dispersion

import { useState, useEffect, useRef, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { CrystalData } from "../lib/types";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";
import { PhononPlot, PhononDOSPlot, PhononDispersionPlot } from "./PhononPlot";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/scfSorting";
import { ProgressBar } from "./ProgressBar";
import { defaultProgressState, progressReducer, ProgressState } from "../lib/qeProgress";

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

// Helper to generate calculation feature tags from parameters
function getCalculationTags(calc: CalculationRun): { label: string; type: "info" | "feature" | "special" | "geometry" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" | "geometry" }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: "info" | "feature" | "special" | "geometry") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (calc.tags) {
    for (const tag of calc.tags) {
      if (tag === "phonon-ready") {
        pushTag("Phonon-Ready", "special");
      } else if (tag === "structure-optimized") {
        pushTag("Optimized", "special");
      } else if (tag === "geometry") {
        pushTag("Geometry", "geometry");
      }
    }
  }

  if (params.structure_source?.type === "optimization") {
    pushTag("Geometry", "geometry");
  }

  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }

  if (params.conv_thr) {
    const thr = params.conv_thr;
    const label = thr < 0.001 ? thr.toExponential(0) : thr.toString();
    pushTag(label, "info");
  }

  // Mark phonon-ready calculations
  if (params.conv_thr && params.conv_thr <= 1e-10) {
    pushTag("Phonon-Ready", "special");
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

interface PhononDOS {
  frequencies: number[];
  dos: number[];
  omega_max: number;
  omega_min: number;
}

interface PhononHighSymmetryMarker {
  q_distance: number;
  label: string;
}

interface PhononDispersion {
  q_points: number[];
  frequencies: number[][];
  high_symmetry_points: PhononHighSymmetryMarker[];
  n_modes: number;
  n_qpoints: number;
  frequency_range: [number, number];
}

interface PhononResult {
  converged: boolean;
  n_qpoints: number;
  n_modes: number;
  dos_data: PhononDOS | null;
  dispersion_data: PhononDispersion | null;
  raw_output: string;
}

interface PhononWizardProps {
  onBack: () => void;
  qePath: string;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
}

type WizardStep = "source" | "qgrid" | "options" | "run" | "results";

export function PhononWizard({
  onBack,
  qePath: _qePath,
  projectId,
  cifId: _cifId,
  crystalData,
  scfCalculations,
}: PhononWizardProps) {
  // Wizard state
  const [step, setStep] = useState<WizardStep>("source");
  const [error, setError] = useState<string | null>(null);

  // Step 1: Source SCF
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());

  const handleScfSortModeChange = useCallback((mode: ScfSortMode) => {
    setScfSortMode(mode);
    setStoredSortMode(mode);
  }, []);

  // Step 2: Q-Grid
  const [qGrid, setQGrid] = useState<[number, number, number]>([4, 4, 4]);
  const [tr2Ph, setTr2Ph] = useState<number>(1e-12);

  // Step 3: Options
  const [calculateDos, setCalculateDos] = useState(true);
  const [dosGrid, setDosGrid] = useState<[number, number, number]>([20, 20, 20]);
  const [calculateDispersion, setCalculateDispersion] = useState(true);
  const [qPath, setQPath] = useState<KPathPoint[]>([]);
  const [pointsPerSegment, setPointsPerSegment] = useState(20);

  // Step 4: Running
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const outputRef = useRef<HTMLPreElement>(null);
  const followOutputRef = useRef(true);

  const handleOutputScroll = () => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    followOutputRef.current = distanceToBottom <= 16;
  };
  const [progress, setProgress] = useState<ProgressState>({
    status: "idle",
    percent: null,
    phase: "Phonon calculation",
  });

  // Step 5: Results
  const [phononResult, setPhononResult] = useState<PhononResult | null>(null);
  const [showOutput, setShowOutput] = useState(false);
  const [isSaved, setIsSaved] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string>("");

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);

  // Load MPI info
  useEffect(() => {
    async function init() {
      try {
        const count = await invoke<number>("get_cpu_count");
        setCpuCount(count);
        setMpiProcs(Math.max(1, Math.floor(count * 0.75)));

        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);
      } catch (e) {
        console.error("Failed to initialize:", e);
      }
    }
    init();
  }, []);

  // Auto-scroll output only if user is at the bottom
  useEffect(() => {
    const el = outputRef.current;
    if (!el || !followOutputRef.current) return;
    el.scrollTop = el.scrollHeight;
  }, [output]);

  // Handle Q-path changes from the BZ viewer
  const handleQPathChange = useCallback((newPath: KPathPoint[]) => {
    setQPath(newPath);
  }, []);

  // Run the phonon calculation
  const runCalculation = async () => {
    if (!selectedScf?.result) {
      setError("No source SCF calculation selected");
      return;
    }

    const shouldCalculateDispersion = calculateDispersion && qPath.length >= 2;
    if (calculateDispersion && !shouldCalculateDispersion) {
      setError("Select at least 2 Q-path points for dispersion, or disable dispersion.");
      return;
    }

    setIsRunning(true);
    followOutputRef.current = true;
    setOutput("");
    setError(null);
    setPhononResult(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Phonon calculation", {
      phonon: {
        hasDos: calculateDos,
        hasDispersion: shouldCalculateDispersion,
      },
    }));
    setCalcStartTime(new Date().toISOString());
    setStep("run");

    let unlisten: UnlistenFn | null = null;

    try {
      // Listen for output events
      unlisten = await listen<string>("qe-output-line", (event) => {
        setOutput((prev) => prev + event.payload + "\n");
        setProgress((prev) => progressReducer("phonon", event.payload, prev));
      });

      // Build the phonon configuration
      const phononConfig = {
        phonon: {
          prefix: "qcortado_scf",
          outdir: "./tmp",
          fildyn: "dynmat",
          nq: qGrid,
          tr2_ph: tr2Ph,
          ldisp: true,
          recover: false,
          asr: "crystal",
        },
        calculate_dos: calculateDos,
        dos_grid: calculateDos ? dosGrid : null,
        calculate_dispersion: shouldCalculateDispersion,
        q_path: shouldCalculateDispersion ? qPath.map(p => ({
          label: p.label,
          coords: p.coords,
          npoints: p.npoints,
        })) : null,
        points_per_segment: pointsPerSegment,
        project_id: projectId,
        scf_calc_id: selectedScf.id,
      };

      // Call the backend
      const result = await invoke<PhononResult>("run_phonon_calculation", {
        config: phononConfig,
        workingDir: "/tmp/qcortado_phonon",
        mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      });

      const endTime = new Date().toISOString();
      setPhononResult(result);
      setStep("results");
      setProgress((prev) => ({
        ...prev,
        status: "complete",
        percent: 100,
        phase: "Complete",
      }));

      // Auto-save the phonon calculation to the project
      try {
        setIsSaving(true);
        const pathString = shouldCalculateDispersion ? qPath.map((p) => p.label).join(" -> ") : "";

        await invoke("save_calculation", {
          projectId,
          cifId: _cifId,
          calcData: {
            calc_type: "phonon",
            parameters: {
              source_scf_id: selectedScf.id,
              q_grid: qGrid,
              tr2_ph: tr2Ph,
              calculate_dos: calculateDos,
              dos_grid: calculateDos ? dosGrid : null,
              calculate_dispersion: shouldCalculateDispersion,
              q_path: pathString,
              points_per_segment: pointsPerSegment,
              n_qpoints: result.n_qpoints,
              n_modes: result.n_modes,
            },
            result: {
              converged: result.converged,
              total_energy: null,
              fermi_energy: null,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: result.raw_output,
              phonon_data: {
                dos_data: result.dos_data,
                dispersion_data: result.dispersion_data,
              },
            },
            started_at: calcStartTime,
            completed_at: endTime,
            input_content: "",
            output_content: output,
          },
          workingDir: "/tmp/qcortado_phonon",
        });
        setIsSaved(true);
      } catch (saveError) {
        console.error("Failed to save phonon calculation:", saveError);
      } finally {
        setIsSaving(false);
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
      if (unlisten) unlisten();
      setIsRunning(false);
    }
  };

  // Check if SCF is phonon-ready (conv_thr <= 1e-10)
  const isPhononReady = (calc: CalculationRun): boolean => {
    const convThr = calc.parameters?.conv_thr;
    return convThr && convThr <= 1e-10;
  };

  // Render current step
  const renderStep = () => {
    switch (step) {
      case "source":
        return renderSourceStep();
      case "qgrid":
        return renderQGridStep();
      case "options":
        return renderOptionsStep();
      case "run":
        return renderRunStep();
      case "results":
        return renderResultsStep();
    }
  };

  // Step 1: Select source SCF
  const renderSourceStep = () => {
    const validScfs = scfCalculations.filter(c => c.calc_type === "scf" && c.result?.converged);
    const sortedScfs = sortScfByMode(validScfs, scfSortMode);
    const phononReadyScfs = validScfs.filter(isPhononReady);

    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Phonon calculations require a converged SCF calculation.
            Please run an SCF calculation first.
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
            <label htmlFor="phonon-scf-sort">Sort SCFs</label>
            <select
              id="phonon-scf-sort"
              value={scfSortMode}
              onChange={(e) => handleScfSortModeChange(e.target.value as ScfSortMode)}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose an SCF calculation to use as the starting point for phonons.
          For best results, use a phonon-ready SCF (conv_thr &lt;= 1e-10).
        </p>

        {phononReadyScfs.length === 0 && (
          <div className="warning-banner">
            No phonon-ready SCF calculations found. Consider running an SCF with
            the "Phonon-Ready" preset (conv_thr = 1e-12) for accurate phonons.
          </div>
        )}

        <div className="scf-list">
          {sortedScfs.map((scf) => {
            const ready = isPhononReady(scf);
            return (
              <div
                key={scf.id}
                className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""} ${ready ? "phonon-ready" : ""}`}
                onClick={() => setSelectedScf(scf)}
              >
                <div className="scf-option-header">
                  <input
                    type="radio"
                    checked={selectedScf?.id === scf.id}
                    onChange={() => setSelectedScf(scf)}
                  />
                  <span className="scf-date">
                    {new Date(scf.started_at).toLocaleDateString()}
                  </span>
                  {ready && <span className="phonon-ready-badge">Phonon-ready</span>}
                </div>
                <div className="scf-details">
                  <span>E = {scf.result?.total_energy?.toFixed(6)} Ry</span>
                  {scf.result?.fermi_energy && (
                    <span>E_F = {scf.result.fermi_energy.toFixed(3)} eV</span>
                  )}
                </div>
                <div className="calc-tags">
                  {getCalculationTags(scf).map((tag, i) => (
                    <span key={i} className={`calc-tag calc-tag-${tag.type}`}>
                      {tag.label}
                    </span>
                  ))}
                </div>
              </div>
            );
          })}
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Cancel
          </button>
          <button
            className="primary-button"
            disabled={!selectedScf}
            onClick={() => setStep("qgrid")}
          >
            Next: Q-Grid
          </button>
        </div>
      </div>
    );
  };

  // Step 2: Q-Grid configuration
  const renderQGridStep = () => {
    const totalQPoints = qGrid[0] * qGrid[1] * qGrid[2];

    return (
      <div className="wizard-step qgrid-step">
        <h3>Q-Point Grid Configuration</h3>
        <p className="step-description">
          Configure the q-point grid for phonon calculations.
          More q-points = more accurate but slower.
        </p>

        <div className="param-section">
          <label>Q-point grid:</label>
          <div className="qgrid-inputs">
            <input
              type="number"
              min={1}
              max={20}
              value={qGrid[0]}
              onChange={(e) => setQGrid([parseInt(e.target.value) || 1, qGrid[1], qGrid[2]])}
            />
            <span>x</span>
            <input
              type="number"
              min={1}
              max={20}
              value={qGrid[1]}
              onChange={(e) => setQGrid([qGrid[0], parseInt(e.target.value) || 1, qGrid[2]])}
            />
            <span>x</span>
            <input
              type="number"
              min={1}
              max={20}
              value={qGrid[2]}
              onChange={(e) => setQGrid([qGrid[0], qGrid[1], parseInt(e.target.value) || 1])}
            />
          </div>
          <span className="param-hint">{totalQPoints} total q-points</span>
        </div>

        <div className="param-section">
          <label>
            Convergence threshold (tr2_ph):
            <select
              value={tr2Ph}
              onChange={(e) => setTr2Ph(parseFloat(e.target.value))}
            >
              <option value={1e-12}>1e-12 (standard)</option>
              <option value={1e-14}>1e-14 (high precision)</option>
              <option value={1e-16}>1e-16 (very high precision)</option>
            </select>
          </label>
        </div>

        {error && <div className="error-message">{error}</div>}

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button className="primary-button" onClick={() => setStep("options")}>
            Next: Options
          </button>
        </div>
      </div>
    );
  };

  // Step 3: Options (DOS + Dispersion)
  const renderOptionsStep = () => {
    return (
      <div className="wizard-step options-step">
        <h3>Calculation Options</h3>

        {/* DOS Options */}
        <div className="option-section">
          <label className="option-checkbox">
            <input
              type="checkbox"
              checked={calculateDos}
              onChange={(e) => setCalculateDos(e.target.checked)}
            />
            <span>Calculate Phonon DOS</span>
          </label>

          {calculateDos && (
            <div className="option-params">
              <label>DOS sampling grid:</label>
              <div className="qgrid-inputs">
                <input
                  type="number"
                  min={5}
                  max={50}
                  value={dosGrid[0]}
                  onChange={(e) => setDosGrid([parseInt(e.target.value) || 10, dosGrid[1], dosGrid[2]])}
                />
                <span>x</span>
                <input
                  type="number"
                  min={5}
                  max={50}
                  value={dosGrid[1]}
                  onChange={(e) => setDosGrid([dosGrid[0], parseInt(e.target.value) || 10, dosGrid[2]])}
                />
                <span>x</span>
                <input
                  type="number"
                  min={5}
                  max={50}
                  value={dosGrid[2]}
                  onChange={(e) => setDosGrid([dosGrid[0], dosGrid[1], parseInt(e.target.value) || 10])}
                />
              </div>
            </div>
          )}
        </div>

        {/* Dispersion Options */}
        <div className="option-section">
          <label className="option-checkbox">
            <input
              type="checkbox"
              checked={calculateDispersion}
              onChange={(e) => setCalculateDispersion(e.target.checked)}
            />
            <span>Calculate Phonon Dispersion</span>
          </label>

          {calculateDispersion && (
            <div className="option-params">
              <div className="crystal-info">
                {crystalData.space_group_HM && (
                  <div className="info-row">
                    <span className="label">Space Group:</span>
                    <span className="value">
                      {crystalData.space_group_HM}
                      {crystalData.space_group_IT_number && ` (#${crystalData.space_group_IT_number})`}
                    </span>
                  </div>
                )}
              </div>

              <BrillouinZoneViewer
                crystalData={crystalData}
                onPathChange={handleQPathChange}
                initialPath={qPath}
                pointsPerSegment={pointsPerSegment}
              />

              <div className="points-per-segment">
                <label>
                  Points per segment:
                  <input
                    type="number"
                    min={5}
                    max={100}
                    value={pointsPerSegment}
                    onChange={(e) => setPointsPerSegment(parseInt(e.target.value) || 20)}
                  />
                </label>
              </div>
            </div>
          )}
        </div>

        {/* MPI Settings */}
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
                      onChange={(e) => setMpiProcs(parseInt(e.target.value) || 1)}
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
          <button className="secondary-button" onClick={() => setStep("qgrid")}>
            Back
          </button>
          <button
            className="primary-button"
            disabled={!calculateDos && !calculateDispersion}
            onClick={runCalculation}
          >
            Run Calculation
          </button>
        </div>
      </div>
    );
  };

  // Step 4: Running
  const renderRunStep = () => {
    return (
      <div className="wizard-step run-step">
        <h3>{isRunning ? "Running Phonon Calculation" : "Phonon Output"}</h3>
        <p className="step-description">
          This may take a while depending on system size and q-grid density.
        </p>

        <ProgressBar
          status={progress.status}
          percent={progress.percent}
          phase={progress.phase}
          detail={progress.detail}
        />

        <pre ref={outputRef} className="calculation-output" onScroll={handleOutputScroll}>
          {output || "Starting calculation..."}
        </pre>

        {error && (
          <div className="error-actions">
            <button className="secondary-button" onClick={() => setStep("options")}>
              Back to Options
            </button>
          </div>
        )}
      </div>
    );
  };

  // Step 5: Results
  const renderResultsStep = () => {
    if (!phononResult) {
      return (
        <div className="wizard-step results-step">
          <h3>No Results</h3>
          <button className="secondary-button" onClick={() => setStep("options")}>
            Try Again
          </button>
        </div>
      );
    }

    const hasDos = phononResult.dos_data !== null;
    const hasDispersion = phononResult.dispersion_data !== null;

    return (
      <div className="wizard-step results-step">
        <h3>Phonon Calculation Results</h3>

        {/* Combined plot if both DOS and dispersion are available */}
        {hasDos && hasDispersion && phononResult.dos_data && phononResult.dispersion_data && (
          <div className="phonon-plot-wrapper combined">
            <PhononPlot
              dos={phononResult.dos_data}
              dispersion={phononResult.dispersion_data}
              width={900}
              height={500}
            />
          </div>
        )}

        {/* DOS only */}
        {hasDos && !hasDispersion && phononResult.dos_data && (
          <div className="phonon-plot-wrapper dos-only">
            <PhononDOSPlot
              data={phononResult.dos_data}
              width={400}
              height={500}
            />
          </div>
        )}

        {/* Dispersion only */}
        {hasDispersion && !hasDos && phononResult.dispersion_data && (
          <div className="phonon-plot-wrapper dispersion-only">
            <PhononDispersionPlot
              data={phononResult.dispersion_data}
              width={700}
              height={500}
            />
          </div>
        )}

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Q-points:</span>
              <span className="value">{phononResult.n_qpoints}</span>
            </div>
            <div className="summary-item">
              <span className="label">Modes:</span>
              <span className="value">{phononResult.n_modes}</span>
            </div>
            {hasDispersion && phononResult.dispersion_data && (
              <div className="summary-item">
                <span className="label">Frequency Range:</span>
                <span className="value">
                  {phononResult.dispersion_data.frequency_range[0].toFixed(1)} to{" "}
                  {phononResult.dispersion_data.frequency_range[1].toFixed(1)} cm^-1
                </span>
              </div>
            )}
          </div>
        </div>

        {/* Collapsible calculation output */}
        <div className="output-section">
          <button
            className="output-toggle"
            onClick={() => setShowOutput(!showOutput)}
          >
            {showOutput ? "v" : ">"} Calculation Output
          </button>
          {showOutput && (
            <pre className="calculation-output results-output">
              {output || "No output available"}
            </pre>
          )}
        </div>

        {/* Save status */}
        <div className="save-status">
          {isSaving && <span className="saving">Saving to project...</span>}
          {isSaved && <span className="saved">Saved to project</span>}
        </div>

        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>
            Back to Dashboard
          </button>
          <button className="primary-button" onClick={() => setStep("options")}>
            Run Another
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className="phonon-wizard">
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Phonon Calculation Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "qgrid" ? "active" : ["options", "run", "results"].includes(step) ? "completed" : ""}>
            2. Q-Grid
          </span>
          <span className={step === "options" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            3. Options
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "completed" : ""}>
            4. Run
          </span>
          <span className={step === "results" ? "active" : ""}>
            5. Results
          </span>
        </div>
      </div>

      <div className="wizard-content">
        {renderStep()}
      </div>
    </div>
  );
}
