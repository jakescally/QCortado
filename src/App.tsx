import { useEffect, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import "./App.css";
import { SCFWizard } from "./components/SCFWizard";
import { BandStructureWizard } from "./components/BandStructureWizard";
import { BandData, BandPlot } from "./components/BandPlot";
import { PhononWizard } from "./components/PhononWizard";
import { PhononDOSPlot } from "./components/PhononPlot";
import { ProjectBrowser } from "./components/ProjectBrowser";
import { ProjectDashboard, CalculationRun } from "./components/ProjectDashboard";
import { CreateProjectDialog } from "./components/CreateProjectDialog";
import { CrystalData, SCFPreset, OptimizedStructureOption } from "./lib/types";

interface ProjectSummary {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  formula: string | null;
  calculation_count: number;
  last_activity: string;
}

type AppView = "home" | "scf-wizard" | "bands-wizard" | "bands-viewer" | "phonon-wizard" | "phonon-viewer" | "project-browser" | "project-dashboard";

interface SCFContext {
  cifId: string;
  crystalData: CrystalData;
  cifContent: string;
  filename: string;
  projectId: string;
  initialPreset?: SCFPreset;
  presetLock?: boolean;
  optimizedStructures?: OptimizedStructureOption[];
}

interface BandsContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface PhononsContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

interface PhononData {
  dos_data: any | null;
  dispersion_data: any | null;
}

type PhononViewMode = "bands" | "dos";
type PhononFrequencyUnit = "cm-1" | "thz";
type PhononBandFocus = "full" | "acoustic" | "optical";
const CM1_TO_THZ = 0.0299792458;

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

function App() {
  const [qePath, setQePath] = useState<string | null>(null);
  const [availableExecutables, setAvailableExecutables] = useState<string[]>([]);
  const [status, setStatus] = useState<string>("Not configured");
  const [error, setError] = useState<string | null>(null);
  const [currentView, setCurrentView] = useState<AppView>("home");
  const [projectCount, setProjectCount] = useState<number>(0);
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [selectedProjectId, setSelectedProjectId] = useState<string | null>(null);

  // Context for running SCF from a project
  const [scfContext, setScfContext] = useState<SCFContext | null>(null);

  // Context for running Bands from a project
  const [bandsContext, setBandsContext] = useState<BandsContext | null>(null);

  // Context for viewing saved band data
  const [viewBandsData, setViewBandsData] = useState<{ bandData: any; fermiEnergy: number | null } | null>(null);

  // Context for running Phonons from a project
  const [phononsContext, setPhononsContext] = useState<PhononsContext | null>(null);

  // Context for viewing saved phonon data
  const [viewPhononData, setViewPhononData] = useState<{ data: PhononData; mode: PhononViewMode } | null>(null);
  const [phononBandsUnit, setPhononBandsUnit] = useState<PhononFrequencyUnit>("cm-1");
  const [phononBandFocus, setPhononBandFocus] = useState<PhononBandFocus>("full");

  // Check for existing QE configuration on startup
  useEffect(() => {
    checkQEPath();
    loadProjectCount();
  }, []);

  async function loadProjectCount() {
    try {
      const projects = await invoke<ProjectSummary[]>("list_projects");
      setProjectCount(projects.length);
    } catch (e) {
      console.log("Failed to load project count:", e);
    }
  }

  async function checkQEPath() {
    try {
      const path = await invoke<string | null>("get_qe_path");
      if (path) {
        setQePath(path);
        await loadExecutables();
        setStatus("Ready");
      }
    } catch (e) {
      console.log("No QE path configured yet");
    }
  }

  async function loadExecutables() {
    try {
      const exes = await invoke<string[]>("check_qe_executables");
      setAvailableExecutables(exes);
      setError(null);
    } catch (e) {
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
        await invoke("set_qe_path", { path: selected });
        setQePath(selected);
        await loadExecutables();
        setStatus("Ready");
        setError(null);
      }
    } catch (e) {
      setError(String(e));
    }
  }

  if (currentView === "bands-wizard" && qePath && bandsContext) {
    return (
      <BandStructureWizard
        qePath={qePath}
        onBack={() => {
          setCurrentView("project-dashboard");
          setBandsContext(null);
        }}
        projectId={bandsContext.projectId}
        cifId={bandsContext.cifId}
        crystalData={bandsContext.crystalData}
        scfCalculations={bandsContext.scfCalculations}
      />
    );
  }

  if (currentView === "bands-viewer" && viewBandsData) {
    return (
      <div className="bands-viewer-container">
        <div className="bands-viewer-header">
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
        </div>
        <div className="bands-viewer-content">
          <BandPlot
            data={viewBandsData.bandData}
            width={900}
            height={600}
            scfFermiEnergy={viewBandsData.fermiEnergy ?? undefined}
          />
        </div>
      </div>
    );
  }

  if (currentView === "phonon-wizard" && qePath && phononsContext) {
    return (
      <PhononWizard
        qePath={qePath}
        onBack={() => {
          setCurrentView("project-dashboard");
          setPhononsContext(null);
        }}
        projectId={phononsContext.projectId}
        cifId={phononsContext.cifId}
        crystalData={phononsContext.crystalData}
        scfCalculations={phononsContext.scfCalculations}
      />
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
      <div className="phonon-viewer-container">
        <div className="phonon-viewer-header">
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
        </div>
        <div className="phonon-viewer-content">
          {showingBands && displayPhononBandData && activePhononRange ? (
            <BandPlot
              key={plotKey}
              data={displayPhononBandData}
              width={900}
              height={600}
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
            />
          ) : showingDos && hasDos ? (
            <PhononDOSPlot
              data={phononData.dos_data}
              width={400}
              height={600}
            />
          ) : (
            <p>{showingBands ? "No phonon dispersion data available" : "No phonon DOS data available"}</p>
          )}
        </div>
      </div>
    );
  }

  if (currentView === "scf-wizard" && qePath) {
    return (
      <SCFWizard
        qePath={qePath}
        onBack={() => {
          // If we came from a project dashboard, go back there
          if (scfContext) {
            setCurrentView("project-dashboard");
            setScfContext(null);
          } else {
            setCurrentView("home");
          }
        }}
        initialCif={scfContext || undefined}
        initialPreset={scfContext?.initialPreset}
        presetLock={scfContext?.presetLock}
        optimizedStructures={scfContext?.optimizedStructures}
      />
    );
  }

  if (currentView === "project-browser") {
    return (
      <ProjectBrowser
        onBack={() => {
          setCurrentView("home");
          loadProjectCount(); // Refresh count in case projects were added/deleted
        }}
        onSelectProject={(projectId) => {
          setSelectedProjectId(projectId);
          setCurrentView("project-dashboard");
        }}
      />
    );
  }

  if (currentView === "project-dashboard" && selectedProjectId) {
    return (
      <ProjectDashboard
        projectId={selectedProjectId}
        onBack={() => {
          setCurrentView("project-browser");
          setSelectedProjectId(null);
        }}
        onDeleted={() => {
          setCurrentView("project-browser");
          setSelectedProjectId(null);
          loadProjectCount();
        }}
        onRunSCF={(cifId, crystalData, cifContent, filename, preset, presetLock, optimizedStructures) => {
          setScfContext({
            cifId,
            crystalData,
            cifContent,
            filename,
            projectId: selectedProjectId,
            initialPreset: preset,
            presetLock,
            optimizedStructures,
          });
          setCurrentView("scf-wizard");
        }}
        onRunBands={(cifId, crystalData, scfCalculations) => {
          setBandsContext({
            cifId,
            crystalData,
            projectId: selectedProjectId,
            scfCalculations,
          });
          setCurrentView("bands-wizard");
        }}
        onViewBands={(bandData, fermiEnergy) => {
          setViewBandsData({ bandData, fermiEnergy });
          setCurrentView("bands-viewer");
        }}
        onRunPhonons={(cifId, crystalData, scfCalculations) => {
          setPhononsContext({
            cifId,
            crystalData,
            projectId: selectedProjectId,
            scfCalculations,
          });
          setCurrentView("phonon-wizard");
        }}
        onViewPhonons={(phononData, viewMode) => {
          setViewPhononData({
            data: phononData,
            mode: viewMode,
          });
          setCurrentView("phonon-viewer");
        }}
      />
    );
  }

  return (
    <main className="container">
      <header className="header">
        <h1>QCortado</h1>
        <p className="subtitle">A Modern Interface for Quantum ESPRESSO</p>
      </header>

      <section className="config-section">
        <h2>Configuration</h2>
        <div className="config-row">
          <label>QE Installation:</label>
          {qePath ? (
            <span className="path">{qePath}</span>
          ) : (
            <span className="not-set">Not configured</span>
          )}
          <button onClick={selectQEPath}>
            {qePath ? "Change" : "Configure"}
          </button>
        </div>

        <div className="status-row">
          <label>Status:</label>
          <span className={`status ${status === "Ready" ? "ready" : "pending"}`}>
            {status}
          </span>
        </div>

        {error && <div className="error">{error}</div>}

        {availableExecutables.length > 0 && (
          <div className="executables">
            <label>Available programs:</label>
            <div className="exe-list">
              {availableExecutables.map((exe) => (
                <span key={exe} className="exe-tag">{exe}</span>
              ))}
            </div>
          </div>
        )}
      </section>

      {qePath && (
        <>
          <section className="actions-section">
            <h2>Projects</h2>
            <div className="action-grid">
              <button className="action-btn" onClick={() => setShowCreateDialog(true)}>
                <span className="action-icon">+</span>
                <span className="action-label">New Project</span>
                <span className="action-hint">Create a project</span>
              </button>
              <button
                className="action-btn"
                onClick={() => setCurrentView("project-browser")}
                disabled={projectCount === 0}
              >
                <span className="action-icon">{projectCount}</span>
                <span className="action-label">
                  {projectCount === 1 ? "View Project" : "View Projects"}
                </span>
                <span className="action-hint">
                  {projectCount === 0 ? "No projects yet" : "Browse & manage"}
                </span>
              </button>
            </div>
          </section>

          <section className="actions-section">
            <h2>Quick Calculations</h2>
            <div className="action-grid">
              <button className="action-btn" onClick={() => setCurrentView("scf-wizard")}>
                <span className="action-icon">SCF</span>
                <span className="action-label">SCF Calculation</span>
                <span className="action-hint">Import CIF & run</span>
              </button>
              <button className="action-btn" disabled>
                <span className="action-icon">Band</span>
                <span className="action-label">Band Structure</span>
                <span className="action-hint">Coming soon</span>
              </button>
              <button className="action-btn" disabled>
                <span className="action-icon">DOS</span>
                <span className="action-label">Density of States</span>
                <span className="action-hint">Coming soon</span>
              </button>
            </div>
          </section>
        </>
      )}

      <footer className="footer">
        <p>QCortado v0.1.0 - Quantum ESPRESSO 7.5 Interface</p>
      </footer>

      <CreateProjectDialog
        isOpen={showCreateDialog}
        onClose={() => setShowCreateDialog(false)}
        onCreated={() => {
          setShowCreateDialog(false);
          loadProjectCount();
        }}
      />
    </main>
  );
}

export default App;
