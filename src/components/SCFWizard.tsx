// SCF Calculation Wizard - Import CIF, configure, and run SCF calculations

import { useState, useEffect, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { CrystalData, ELEMENT_MASSES } from "../lib/types";
import { parseCIF } from "../lib/cifParser";
import { UnitCellViewer } from "./UnitCellViewer";
import { SaveToProjectDialog } from "./SaveToProjectDialog";

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
}

interface SCFConfig {
  ecutwfc: number;
  ecutrho: number;
  kgrid: [number, number, number];
  conv_thr: number;
  mixing_beta: number;
}

type WizardStep = "import" | "configure" | "run" | "results";

export function SCFWizard({ onBack, qePath, initialCif }: SCFWizardProps) {
  // If we have initial CIF data, skip to configure step
  const [step, setStep] = useState<WizardStep>(initialCif ? "configure" : "import");
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

  // Save dialog state
  const [showSaveDialog, setShowSaveDialog] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string>("");
  const [calcEndTime, setCalcEndTime] = useState<string>("");
  const [generatedInput, setGeneratedInput] = useState<string>("");
  const [resultSaved, setResultSaved] = useState(false);

  // Exit confirmation dialog
  const [showExitConfirm, setShowExitConfirm] = useState(false);
  const [pendingExitAction, setPendingExitAction] = useState<(() => void) | null>(null);

  // Ref for auto-scrolling output
  const outputRef = useRef<HTMLPreElement>(null);

  const [config, setConfig] = useState<SCFConfig>({
    ecutwfc: 40,
    ecutrho: 320,
    kgrid: [4, 4, 4],
    conv_thr: 1e-6,
    mixing_beta: 0.7,
  });

  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState<string>("");
  const [result, setResult] = useState<any>(null);
  const [pseudoDir, setPseudoDir] = useState<string>("");
  const [pseudoError, setPseudoError] = useState<string | null>(null);

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

  // Auto-scroll output when it changes
  useEffect(() => {
    if (outputRef.current) {
      outputRef.current.scrollTop = outputRef.current.scrollHeight;
    }
  }, [output]);

  // Load CPU count and check MPI availability
  useEffect(() => {
    async function loadMpiInfo() {
      try {
        const cores = await invoke<number>("get_cpu_count");
        setCpuCount(cores);
        // Default to ~75% of cores, minimum 1
        setMpiProcs(Math.max(1, Math.floor(cores * 0.75)));

        const mpiOk = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(mpiOk);
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

  // Auto-select pseudopotentials and set cutoffs when crystal data or SSSP data changes
  useEffect(() => {
    if (crystalData && pseudopotentials.length > 0) {
      const selected: Record<string, string> = {};
      const elements = [...new Set(crystalData.atom_sites.map((a) => getBaseElement(a.type_symbol)))];

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
  }, [crystalData, pseudopotentials, ssspData]);

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
    if (!crystalData) return [];
    // Strip oxidation states and get unique base elements
    return [...new Set(crystalData.atom_sites.map((a) => getBaseElement(a.type_symbol)))];
  }

  function canRun(): boolean {
    if (!crystalData) return false;
    const elements = getUniqueElements();
    return elements.every((el) => selectedPseudos[el]);
  }

  async function handleRun() {
    if (!crystalData || !canRun()) return;

    setIsRunning(true);
    setOutput("");
    setResult(null);
    setStep("run");

    // Track calculation start time
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);

    // Set up event listener for streaming output
    let unlisten: UnlistenFn | null = null;

    try {
      const elements = getUniqueElements();

      // Build the calculation configuration
      const calculation = {
        calculation: "scf",
        prefix: "qcortado_scf",
        outdir: "./tmp",
        pseudo_dir: qePath.replace("/bin", "/pseudo"),
        system: {
          ibrav: "free", // Free lattice - use CELL_PARAMETERS
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
          ecutwfc: config.ecutwfc,
          ecutrho: config.ecutrho,
          nspin: 1,
          occupations: "smearing",
          smearing: "gaussian",
          degauss: 0.01,
        },
        kpoints: {
          type: "automatic",
          grid: config.kgrid,
          offset: [0, 0, 0],
        },
        conv_thr: config.conv_thr,
        mixing_beta: config.mixing_beta,
        tprnfor: true,
        tstress: true,
        forc_conv_thr: null,
        etot_conv_thr: null,
        verbosity: "high",
      };

      // Generate and display the input file
      const inputText = await invoke<string>("generate_input", { calculation });
      setGeneratedInput(inputText);
      setOutput(`=== Generated Input ===\n${inputText}\n\n=== Running pw.x ===\n`);

      // Listen for streaming output
      unlisten = await listen<string>("qe-output-line", (event) => {
        setOutput((prev) => prev + event.payload + "\n");
      });

      // Run the calculation with streaming
      const calcResult = await invoke<any>("run_calculation_streaming", {
        calculation,
        workingDir: "/tmp/qcortado_work",
        mpiConfig: mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
      });

      // Track calculation end time
      const endTime = new Date().toISOString();
      setCalcEndTime(endTime);

      setResult(calcResult);
      setOutput((prev) => prev + "\n=== Calculation Complete ===\n");
      setStep("results");
    } catch (e) {
      setError(`Calculation failed: ${e}`);
      setOutput((prev) => prev + `\n\nERROR: ${e}`);
    } finally {
      // Clean up event listener
      if (unlisten) {
        unlisten();
      }
      setIsRunning(false);
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

  return (
    <div className="wizard-container">
      <div className="wizard-header">
        <button className="back-btn" onClick={() => handleExitAttempt(onBack)}>
          ← Back
        </button>
        <h2>SCF Calculation Wizard</h2>
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
                        a={crystalData.cell_length_a.value.toFixed(3)}, b=
                        {crystalData.cell_length_b.value.toFixed(3)}, c=
                        {crystalData.cell_length_c.value.toFixed(3)}
                      </span>
                    </div>
                    <div className="info-item">
                      <label>Angles (°)</label>
                      <span>
                        α={crystalData.cell_angle_alpha.value.toFixed(1)}, β=
                        {crystalData.cell_angle_beta.value.toFixed(1)}, γ=
                        {crystalData.cell_angle_gamma.value.toFixed(1)}
                      </span>
                    </div>
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

                {/* SCF Parameters */}
                <section className="config-section">
                  <h3>
                    SCF Parameters
                    <Tooltip text="Self-Consistent Field (SCF) calculations iteratively solve the Kohn-Sham equations until the electron density converges. These parameters control the accuracy and computational cost of the calculation." />
                  </h3>
                  <div className="param-grid">
                    <div className="param-row">
                      <label>
                        Wavefunction Cutoff
                        <Tooltip text="Energy cutoff for plane-wave expansion of wavefunctions (in Rydberg). Higher values = more accurate but slower. Typical range: 30-80 Ry. Start with your pseudopotential's recommended value. Doubling this roughly 8x the computation time." />
                      </label>
                      <div className="param-input-group">
                        <input
                          type="number"
                          value={config.ecutwfc}
                          onChange={(e) =>
                            setConfig((prev) => ({
                              ...prev,
                              ecutwfc: parseFloat(e.target.value),
                            }))
                          }
                        />
                        <span className="param-unit">Ry</span>
                      </div>
                    </div>
                    <div className="param-row">
                      <label>
                        Charge Density Cutoff
                        <Tooltip text="Energy cutoff for charge density and potential (in Rydberg). Should be 4x wavefunction cutoff for norm-conserving pseudopotentials, or 8-12x for ultrasoft/PAW. Higher values improve accuracy of forces and stress." />
                      </label>
                      <div className="param-input-group">
                        <input
                          type="number"
                          value={config.ecutrho}
                          onChange={(e) =>
                            setConfig((prev) => ({
                              ...prev,
                              ecutrho: parseFloat(e.target.value),
                            }))
                          }
                        />
                        <span className="param-unit">Ry</span>
                      </div>
                    </div>
                    <div className="param-row">
                      <label>
                        K-point Grid
                        <Tooltip text="Sampling grid for the Brillouin zone (reciprocal space). Denser grids = more accurate but slower. Metals need denser grids than insulators. For cubic cells, use equal values. Doubling each dimension = 8x more k-points = ~8x slower." />
                      </label>
                      <div className="kgrid-inputs">
                        <input
                          type="number"
                          value={config.kgrid[0]}
                          onChange={(e) =>
                            setConfig((prev) => ({
                              ...prev,
                              kgrid: [parseInt(e.target.value), prev.kgrid[1], prev.kgrid[2]],
                            }))
                          }
                        />
                        <span>×</span>
                        <input
                          type="number"
                          value={config.kgrid[1]}
                          onChange={(e) =>
                            setConfig((prev) => ({
                              ...prev,
                              kgrid: [prev.kgrid[0], parseInt(e.target.value), prev.kgrid[2]],
                            }))
                          }
                        />
                        <span>×</span>
                        <input
                          type="number"
                          value={config.kgrid[2]}
                          onChange={(e) =>
                            setConfig((prev) => ({
                              ...prev,
                              kgrid: [prev.kgrid[0], prev.kgrid[1], parseInt(e.target.value)],
                            }))
                          }
                        />
                      </div>
                    </div>
                    <div className="param-row">
                      <label>
                        Convergence Threshold
                        <Tooltip text="SCF iteration stops when energy change falls below this value (in Rydberg). Smaller = more accurate but more iterations. 1e-6 is good for most calculations. Use 1e-8 or smaller for phonons or precise forces." />
                      </label>
                      <div className="param-input-group">
                        <input
                          type="text"
                          value={config.conv_thr}
                          onChange={(e) =>
                            setConfig((prev) => ({
                              ...prev,
                              conv_thr: parseFloat(e.target.value),
                            }))
                          }
                        />
                        <span className="param-unit">Ry</span>
                      </div>
                    </div>
                  </div>
                  {ssspMissing && (
                    <p className="sssp-hint">
                      Cutoff values are defaults. For optimized values, download the{" "}
                      <a
                        href="https://www.materialscloud.org/discover/sssp/table/precision"
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                        SSSP library
                      </a>{" "}
                      and place the JSON file in your pseudo directory.
                    </p>
                  )}
                  {!ssspMissing && ssspData && (
                    <p className="sssp-success">
                      Cutoffs auto-configured from SSSP library (max values for all elements).
                    </p>
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

                <button
                  className="run-btn"
                  onClick={handleRun}
                  disabled={!canRun()}
                >
                  {mpiEnabled && mpiProcs > 1
                    ? `Run SCF Calculation (${mpiProcs} cores)`
                    : "Run SCF Calculation"}
                </button>
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
            <div className="run-layout">
              <div className="output-panel">
                <h3>{isRunning ? "Running..." : "Output"}</h3>
                <pre className="output-text" ref={outputRef}>{output}</pre>
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
                </div>
              )}
            </div>

            <div className="run-actions">
              <button onClick={() => handleExitAttempt(() => setStep("configure"))} disabled={isRunning}>
                ← Back to Configure
              </button>
              {result && (
                <>
                  <button
                    onClick={() => setShowSaveDialog(true)}
                    className="save-project-btn"
                  >
                    {resultSaved ? "Saved" : "Save to Project"}
                  </button>
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
      {showSaveDialog && crystalData && result && (
        <SaveToProjectDialog
          isOpen={showSaveDialog}
          onClose={() => setShowSaveDialog(false)}
          onSaved={() => {
            setShowSaveDialog(false);
            setResultSaved(true);
          }}
          calculationData={{
            calc_type: "scf",
            parameters: {
              ecutwfc: config.ecutwfc,
              ecutrho: config.ecutrho,
              kgrid: config.kgrid,
              conv_thr: config.conv_thr,
              mixing_beta: config.mixing_beta,
            },
            result: result,
            started_at: calcStartTime,
            completed_at: calcEndTime,
            input_content: generatedInput,
            output_content: result.raw_output || "",
          }}
          cifData={{
            filename: cifFilename,
            formula: crystalData.formula_sum || crystalData.formula_structural || "Unknown",
            content: cifContent,
            crystal_data: crystalData,
          }}
          workingDir="/tmp/qcortado_work"
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
                Your calculation results have not been saved to a project.
                If you leave now, the results will be lost.
              </p>
              <p className="exit-hint">
                Click "Save to Project" to keep your results, or "Leave Anyway" to discard them.
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
