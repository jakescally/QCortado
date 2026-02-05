// SCF Calculation Wizard - Import CIF, configure, and run SCF calculations

import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { CrystalData, ELEMENT_MASSES } from "../lib/types";
import { parseCIF } from "../lib/cifParser";
import { UnitCellViewer } from "./UnitCellViewer";

interface SCFWizardProps {
  onBack: () => void;
  qePath: string;
}

interface SCFConfig {
  ecutwfc: number;
  ecutrho: number;
  kgrid: [number, number, number];
  conv_thr: number;
  mixing_beta: number;
}

type WizardStep = "import" | "configure" | "run" | "results";

export function SCFWizard({ onBack, qePath }: SCFWizardProps) {
  const [step, setStep] = useState<WizardStep>("import");
  const [crystalData, setCrystalData] = useState<CrystalData | null>(null);
  const [_cifFilename, setCifFilename] = useState<string>("");
  const [error, setError] = useState<string | null>(null);
  const [pseudopotentials, setPseudopotentials] = useState<string[]>([]);
  const [selectedPseudos, setSelectedPseudos] = useState<Record<string, string>>({});

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

  // Load available pseudopotentials
  useEffect(() => {
    async function loadPseudos() {
      try {
        const pseudoDir = qePath.replace("/bin", "/pseudo");
        const pseudos = await invoke<string[]>("list_pseudopotentials", {
          pseudoDir,
        });
        setPseudopotentials(pseudos);
      } catch (e) {
        console.error("Failed to load pseudopotentials:", e);
      }
    }
    loadPseudos();
  }, [qePath]);

  // Auto-select pseudopotentials when crystal data changes
  useEffect(() => {
    if (crystalData && pseudopotentials.length > 0) {
      const selected: Record<string, string> = {};
      const elements = [...new Set(crystalData.atom_sites.map((a) => a.type_symbol))];

      for (const element of elements) {
        // Try to find a matching pseudopotential
        const matches = pseudopotentials.filter((p) =>
          p.toLowerCase().startsWith(element.toLowerCase() + ".")
        );
        if (matches.length > 0) {
          // Prefer PBE pseudopotentials
          const pbe = matches.find((m) => m.toLowerCase().includes("pbe"));
          selected[element] = pbe || matches[0];
        }
      }
      setSelectedPseudos(selected);
    }
  }, [crystalData, pseudopotentials]);

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
        setError(null);
        setStep("configure");
      }
    } catch (e) {
      setError(`Failed to import CIF: ${e}`);
    }
  }

  function getUniqueElements(): string[] {
    if (!crystalData) return [];
    return [...new Set(crystalData.atom_sites.map((a) => a.type_symbol))];
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

    try {
      const elements = getUniqueElements();

      // Build the calculation configuration
      const calculation = {
        calculation: "scf",
        prefix: "qcortado_scf",
        outdir: "./tmp",
        pseudo_dir: qePath.replace("/bin", "/pseudo"),
        system: {
          ibrav: 0, // Free lattice - use CELL_PARAMETERS
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
            symbol: site.type_symbol,
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
      setOutput(`=== Generated Input ===\n${inputText}\n\n=== Running... ===\n`);

      // Run the calculation
      const calcResult = await invoke<any>("run_calculation", {
        calculation,
        workingDir: "/tmp/qcortado_work",
      });

      setResult(calcResult);
      setOutput((prev) => prev + "\n=== Calculation Complete ===\n" + calcResult.raw_output);
      setStep("results");
    } catch (e) {
      setError(`Calculation failed: ${e}`);
      setOutput((prev) => prev + `\n\nERROR: ${e}`);
    } finally {
      setIsRunning(false);
    }
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
        <button className="back-btn" onClick={onBack}>
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
                  <h3>Pseudopotentials</h3>
                  <div className="pseudo-list">
                    {getUniqueElements().map((el) => (
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
                          <option value="">Select...</option>
                          {pseudopotentials
                            .filter((p) =>
                              p.toLowerCase().startsWith(el.toLowerCase() + ".")
                            )
                            .map((p) => (
                              <option key={p} value={p}>
                                {p}
                              </option>
                            ))}
                        </select>
                      </div>
                    ))}
                  </div>
                </section>

                {/* SCF Parameters */}
                <section className="config-section">
                  <h3>SCF Parameters</h3>
                  <div className="param-grid">
                    <div className="param-row">
                      <label>ecutwfc (Ry)</label>
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
                    </div>
                    <div className="param-row">
                      <label>ecutrho (Ry)</label>
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
                    </div>
                    <div className="param-row">
                      <label>K-grid</label>
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
                      <label>conv_thr</label>
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
                    </div>
                  </div>
                </section>

                <button
                  className="run-btn"
                  onClick={handleRun}
                  disabled={!canRun()}
                >
                  Run SCF Calculation
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
                <pre className="output-text">{output}</pre>
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
              <button onClick={() => setStep("configure")} disabled={isRunning}>
                ← Back to Configure
              </button>
              {result && (
                <button onClick={() => setStep("import")} className="new-calc-btn">
                  New Calculation
                </button>
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
