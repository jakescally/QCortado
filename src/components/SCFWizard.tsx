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
  // Basic parameters
  ecutwfc: number;
  ecutrho: number;
  kgrid: [number, number, number];
  kgrid_offset: [number, number, number];

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
    // Basic parameters
    ecutwfc: 40,
    ecutrho: 320,
    kgrid: [4, 4, 4],
    kgrid_offset: [0, 0, 0],

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

      // Build the calculation configuration with all options
      const calculation: any = {
        calculation: "scf",
        prefix: "qcortado_scf",
        outdir: "./tmp",
        pseudo_dir: qePath.replace("/bin", "/pseudo"),
        verbosity: config.verbosity,
        tprnfor: config.tprnfor,
        tstress: config.tstress,
        disk_io: config.disk_io,
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
            starting_magnetization: config.starting_magnetization[el] || 0,
            hubbard_u: config.lda_plus_u ? (config.hubbard_u[el] || 0) : undefined,
            hubbard_j: config.lda_plus_u && config.lda_plus_u_kind > 0 ? (config.hubbard_j[el] || 0) : undefined,
          })),
          atoms: crystalData.atom_sites.map((site) => ({
            symbol: getBaseElement(site.type_symbol),
            position: [site.fract_x, site.fract_y, site.fract_z],
            if_pos: [true, true, true],
          })),
          position_units: "crystal",
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
        },
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
        forc_conv_thr: null,
        etot_conv_thr: null,
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

                {/* Basic Parameters */}
                <section className="config-section collapsible">
                  <h3 onClick={() => toggleSection("basic")} className="section-header">
                    <span className={`collapse-icon ${expandedSections.basic ? "expanded" : ""}`}>▶</span>
                    Basic Parameters
                    <Tooltip text="Essential parameters for any SCF calculation: energy cutoffs and k-point sampling." />
                  </h3>
                  {expandedSections.basic && (
                    <div className="param-grid">
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
              prefix: "qcortado_scf",
              ecutwfc: config.ecutwfc,
              ecutrho: config.ecutrho,
              kgrid: config.kgrid,
              conv_thr: config.conv_thr,
              mixing_beta: config.mixing_beta,
              // Feature flags for tags
              nspin: config.nspin,
              lspinorb: config.lspinorb,
              lda_plus_u: config.lda_plus_u,
              vdw_corr: config.vdw_corr,
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
