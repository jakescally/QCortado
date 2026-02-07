// Band Structure Wizard - Calculate and visualize electronic band structures

import { useState, useEffect, useRef, useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { CrystalData, ELEMENT_MASSES } from "../lib/types";
import { BandPlot, BandData } from "./BandPlot";
import { BrillouinZoneViewer, KPathPoint } from "./BrillouinZoneViewer";

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
}

interface BandStructureWizardProps {
  onBack: () => void;
  qePath: string;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
}

type WizardStep = "source" | "kpath" | "parameters" | "run" | "results";

export function BandStructureWizard({
  onBack,
  qePath,
  projectId,
  cifId: _cifId,
  crystalData,
  scfCalculations,
}: BandStructureWizardProps) {
  // Wizard state
  const [step, setStep] = useState<WizardStep>("source");
  const [error, setError] = useState<string | null>(null);

  // Step 1: Source SCF
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);

  // Step 2: K-Path (using BrillouinZoneViewer)
  const [kPath, setKPath] = useState<KPathPoint[]>([]);
  const [pointsPerSegment, setPointsPerSegment] = useState(20);

  // Step 3: Parameters
  const [nbnd, setNbnd] = useState<number | "auto">("auto");

  // Step 4: Running
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const outputRef = useRef<HTMLPreElement>(null);

  // Step 5: Results
  const [bandData, setBandData] = useState<BandData | null>(null);
  // Store SCF Fermi energy separately to ensure it persists
  const [scfFermiEnergy, setScfFermiEnergy] = useState<number | null>(null);

  // MPI settings
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [cpuCount, setCpuCount] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);

  // Pseudopotentials (pseudopotentials list used internally for auto-selection)
  const [, setPseudopotentials] = useState<string[]>([]);
  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});

  // Helper functions
  function getBaseElement(symbol: string): string {
    return symbol.replace(/[\d+-]+$/, "");
  }

  function pseudoMatchesElement(filename: string, element: string): boolean {
    const lowerFile = filename.toLowerCase();
    const lowerEl = getBaseElement(element).toLowerCase();
    return (
      lowerFile.startsWith(lowerEl + ".") ||
      lowerFile.startsWith(lowerEl + "_") ||
      lowerFile.startsWith(lowerEl + "-")
    );
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

  // Load MPI info and pseudopotentials
  useEffect(() => {
    async function init() {
      try {
        // Check MPI
        const count = await invoke<number>("get_cpu_count");
        setCpuCount(count);
        setMpiProcs(Math.max(1, Math.floor(count * 0.75)));

        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);

        // Load pseudopotentials
        const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
        const pseudos = await invoke<string[]>("list_pseudopotentials", { pseudoDir });
        setPseudopotentials(pseudos);

        // Auto-select pseudopotentials for elements in the structure
        const elements = [...new Set(crystalData.atom_sites.map(s => getBaseElement(s.type_symbol)))];
        const selected: Record<string, string> = {};
        for (const el of elements) {
          const match = pseudos.find(p => pseudoMatchesElement(p, el));
          if (match) {
            selected[el] = match;
          }
        }
        setSelectedPseudos(selected);
      } catch (e) {
        console.error("Failed to initialize:", e);
      }
    }
    init();
  }, [qePath, crystalData]);

  // Auto-scroll output
  useEffect(() => {
    if (outputRef.current) {
      outputRef.current.scrollTop = outputRef.current.scrollHeight;
    }
  }, [output]);

  // Handle k-path changes from the BZ viewer
  const handleKPathChange = useCallback((newPath: KPathPoint[]) => {
    setKPath(newPath);
  }, []);

  // Run the calculation
  const runCalculation = async () => {
    if (!selectedScf?.result) {
      setError("No source SCF calculation selected");
      return;
    }

    setIsRunning(true);
    setOutput("");
    setError(null);
    setBandData(null);
    setStep("run");

    let unlisten: UnlistenFn | null = null;

    try {
      // Listen for output events
      unlisten = await listen<string>("qe-output-line", (event) => {
        setOutput((prev) => prev + event.payload + "\n");
      });

      // Validate k-path
      if (kPath.length < 2) {
        throw new Error("Please select at least 2 points for the k-path");
      }

      // Get the SCF parameters for cutoffs, etc.
      const scfParams = selectedScf.parameters || {};

      // Get unique elements from crystal data
      const elements = [...new Set(crystalData.atom_sites.map(s => getBaseElement(s.type_symbol)))];

      // Check we have pseudopotentials for all elements
      for (const el of elements) {
        if (!selectedPseudos[el]) {
          throw new Error(`No pseudopotential selected for element ${el}`);
        }
      }

      // Build the full calculation config from crystalData
      // Use the same prefix as the source SCF so we can read its .save directory
      // SCFWizard uses "qcortado_scf" as the prefix
      const scfPrefix = scfParams.prefix || "qcortado_scf";
      const pseudoDir = qePath.replace(/\/bin\/?$/, "/pseudo");
      const baseCalculation = {
        calculation: "scf",
        prefix: scfPrefix,
        outdir: "./tmp",
        pseudo_dir: pseudoDir,
        system: {
          ibrav: "free",
          celldm: null,
          cell_parameters: [
            [crystalData.cell_length_a.value, 0, 0],
            [
              crystalData.cell_length_b.value *
                Math.cos((crystalData.cell_angle_gamma.value * Math.PI) / 180),
              crystalData.cell_length_b.value *
                Math.sin((crystalData.cell_angle_gamma.value * Math.PI) / 180),
              0,
            ],
            calculateCVector(crystalData),
          ],
          cell_units: "angstrom",
          species: elements.map((el) => ({
            symbol: el,
            mass: ELEMENT_MASSES[el] || 1.0,
            pseudopotential: selectedPseudos[el],
          })),
          atoms: crystalData.atom_sites.map((site) => ({
            symbol: getBaseElement(site.type_symbol),
            position: [site.fract_x, site.fract_y, site.fract_z],
            if_pos: [true, true, true],
          })),
          position_units: "crystal",
          ecutwfc: scfParams.ecutwfc || 40,
          ecutrho: scfParams.ecutrho || 320,
          nspin: 1,
          occupations: "smearing",
          smearing: "gaussian",
          degauss: 0.01,
        },
        kpoints: { type: "gamma" }, // Will be replaced by backend
        conv_thr: scfParams.conv_thr || 1e-6,
        mixing_beta: scfParams.mixing_beta || 0.7,
        tprnfor: false,
        tstress: false,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: "high",
      };

      // Call the backend
      const result = await invoke<BandData>("run_bands_calculation", {
        config: {
          base_calculation: baseCalculation,
          k_path: kPath,
          nbnd: nbnd === "auto" ? null : nbnd,
          project_id: projectId,
          scf_calc_id: selectedScf.id,
        },
        workingDir: "/tmp/qcortado_bands",
        mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      });

      setBandData(result);
      setStep("results");
    } catch (e) {
      setError(String(e));
      setOutput((prev) => prev + `\nError: ${e}\n`);
    } finally {
      if (unlisten) unlisten();
      setIsRunning(false);
    }
  };

  // Render current step
  const renderStep = () => {
    switch (step) {
      case "source":
        return renderSourceStep();
      case "kpath":
        return renderKPathStep();
      case "parameters":
        return renderParametersStep();
      case "run":
        return renderRunStep();
      case "results":
        return renderResultsStep();
    }
  };

  // Step 1: Select source SCF
  const renderSourceStep = () => {
    const validScfs = scfCalculations.filter(c => c.calc_type === "scf" && c.result?.converged);

    if (validScfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">
            Band structure calculations require a completed SCF calculation.
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
        <h3>Select Source SCF Calculation</h3>
        <p className="step-description">
          Choose the SCF calculation to use as the starting point for band structure.
        </p>

        <div className="scf-list">
          {validScfs.map((scf) => (
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
                {scf.result?.fermi_energy && (
                  <span>E_F = {scf.result.fermi_energy.toFixed(3)} eV</span>
                )}
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
            onClick={() => setStep("kpath")}
          >
            Next: K-Path
          </button>
        </div>
      </div>
    );
  };

  // Step 2: K-Path selection with 3D Brillouin Zone viewer
  const renderKPathStep = () => {
    return (
      <div className="wizard-step kpath-step">
        <h3>Select K-Path</h3>

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
          onPathChange={handleKPathChange}
          initialPath={kPath}
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

        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>
            Back
          </button>
          <button
            className="primary-button"
            disabled={kPath.length < 2}
            onClick={() => setStep("parameters")}
          >
            Next: Parameters
          </button>
        </div>
      </div>
    );
  };

  // Step 3: Parameters
  const renderParametersStep = () => {
    // Calculate total k-points from the path
    const totalKPoints = kPath.reduce((sum, p) => sum + p.npoints, 0) + kPath.filter(p => p.npoints === 0).length;
    // Format path for display
    const pathString = kPath.map((p) => p.label).join(" → ");

    return (
      <div className="wizard-step parameters-step">
        <h3>Calculation Parameters</h3>

        <div className="param-section">
          <label>
            Number of bands:
            <div className="nbnd-input">
              <select
                value={nbnd === "auto" ? "auto" : "manual"}
                onChange={(e) => setNbnd(e.target.value === "auto" ? "auto" : 20)}
              >
                <option value="auto">Auto (recommended)</option>
                <option value="manual">Manual</option>
              </select>
              {nbnd !== "auto" && (
                <input
                  type="number"
                  min={1}
                  value={nbnd}
                  onChange={(e) => setNbnd(parseInt(e.target.value) || 20)}
                />
              )}
            </div>
          </label>
        </div>

        <div className="calculation-summary">
          <h4>Summary</h4>
          <div className="summary-row">
            <span>K-path:</span>
            <span>{pathString || "Not selected"}</span>
          </div>
          <div className="summary-row">
            <span>Total k-points:</span>
            <span>{totalKPoints}</span>
          </div>
          <div className="summary-row">
            <span>Source SCF energy:</span>
            <span>{selectedScf?.result?.total_energy?.toFixed(6)} Ry</span>
          </div>
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
          <button className="secondary-button" onClick={() => setStep("kpath")}>
            Back
          </button>
          <button className="primary-button" onClick={runCalculation}>
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
        <h3>Running Band Structure Calculation</h3>

        <div className="progress-indicator">
          {isRunning ? (
            <div className="spinner"></div>
          ) : error ? (
            <span className="error-icon">Error</span>
          ) : (
            <span className="success-icon">Complete</span>
          )}
        </div>

        <pre ref={outputRef} className="calculation-output">
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
  };

  // Step 5: Results
  const renderResultsStep = () => {
    if (!bandData) {
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
        <h3>Band Structure Results</h3>

        <div className="band-plot-wrapper">
          <BandPlot
            data={bandData}
            width={800}
            height={500}
            scfFermiEnergy={scfFermiEnergy ?? undefined}
          />
        </div>

        <div className="results-summary">
          <div className="summary-grid">
            <div className="summary-item">
              <span className="label">Bands:</span>
              <span className="value">{bandData.n_bands}</span>
            </div>
            <div className="summary-item">
              <span className="label">K-points:</span>
              <span className="value">{bandData.n_kpoints}</span>
            </div>
            <div className="summary-item">
              <span className="label">Fermi Energy:</span>
              <span className="value">{bandData.fermi_energy.toFixed(3)} eV</span>
            </div>
            {/* Band gap info - disabled for now
            {bandData.band_gap ? (
              <>
                <div className="summary-item">
                  <span className="label">Band Gap:</span>
                  <span className="value gap-value">{bandData.band_gap.value.toFixed(3)} eV</span>
                </div>
                <div className="summary-item">
                  <span className="label">Gap Type:</span>
                  <span className="value">{bandData.band_gap.is_direct ? "Direct" : "Indirect"}</span>
                </div>
              </>
            ) : (
              <div className="summary-item">
                <span className="label">Character:</span>
                <span className="value metal-value">Metallic</span>
              </div>
            )}
            */}
          </div>
        </div>

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
    <div className="band-structure-wizard">
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Band Structure Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "kpath" ? "active" : ["parameters", "run", "results"].includes(step) ? "completed" : ""}>
            2. K-Path
          </span>
          <span className={step === "parameters" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            3. Parameters
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
