import { useEffect, useMemo, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import "./App.css";
import { BandPlot } from "./components/BandPlot";
import type { BandData, BandPlotData } from "./components/BandPlot";
import { ElectronicDOSData, ElectronicDOSPlot } from "./components/ElectronicDOSPlot";
import { EpwViewer } from "./components/qe";
import type { EpwViewerPayload } from "./components/qe";
import { HpcSetupWizard } from "./components/HpcSetupWizard";
import { PhononDOSPlot } from "./components/PhononPlot";
import { ProjectBrowser } from "./components/ProjectBrowser";
import {
  CalculationRun,
  ProjectDashboard,
  SavedBandsCalculationContext,
  WannierBandOverlayOption,
} from "./components/ProjectDashboard";
import { TransportPlot } from "./components/TransportPlot";
import { TaskProvider } from "./lib/TaskContext";
import { ThemeProvider } from "./lib/ThemeContext";
import { getActiveHpcProfileId, listHpcProfiles } from "./lib/hpcConfig";
import { HpcProfile } from "./lib/types";
import { useWindowSize } from "./lib/useWindowSize";
import { BandsMultiview } from "./components/BandsMultiview";
import type { BandsMultiviewCalculation } from "./components/BandsMultiview";
import { formatWannierConvergenceFlag, getWannierQualityIssues } from "./lib/engines/qe/wannierQuality";
import type { TransportResult } from "./lib/transport";

interface ViewerSyncStatus {
  last_synced_at?: string | null;
  last_error?: string | null;
  local_project_count: number;
}

interface ViewerSyncResult {
  synced_at: string;
  downloaded_projects: number;
  removed_projects: number;
  skipped_projects: number;
  total_projects: number;
}

type ViewerView = "home" | "project-browser" | "project-dashboard" | "bands-viewer" | "bands-multiview" | "dos-viewer" | "wannier-viewer" | "transport-viewer" | "phonon-viewer" | "epw-viewer";
type PhononViewMode = "bands" | "dos";

interface PhononData {
  dos_data: any | null;
  dispersion_data: any | null;
}

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

function formatSyncTime(isoString: string | null | undefined): string {
  if (!isoString) return "Never";
  const parsed = new Date(isoString);
  if (Number.isNaN(parsed.getTime())) return isoString;
  return parsed.toLocaleString();
}

function ViewerAppInner() {
  const windowSize = useWindowSize();
  const plotHeight = Math.max(400, windowSize.height - 180);

  const [currentView, setCurrentView] = useState<ViewerView>("home");
  const [selectedProjectId, setSelectedProjectId] = useState<string | null>(null);
  const [projectBrowserFolderId, setProjectBrowserFolderId] = useState<string | null>(null);
  const [bandsMultiviewInitialCalculations, setBandsMultiviewInitialCalculations] =
    useState<BandsMultiviewCalculation[] | undefined>(undefined);
  const [bandsMultiviewReturnView, setBandsMultiviewReturnView] =
    useState<"project-browser" | "project-dashboard">("project-browser");
  const [openBandsMultiviewInitialCalculationsInPreview, setOpenBandsMultiviewInitialCalculationsInPreview] =
    useState(false);
  const [showHpcSetupWizard, setShowHpcSetupWizard] = useState(false);
  const [editingHpcProfileId, setEditingHpcProfileId] = useState<string | null>(null);

  const [hpcProfiles, setHpcProfiles] = useState<HpcProfile[]>([]);
  const [activeHpcProfileId, setActiveHpcProfileId] = useState<string | null>(null);

  const [syncStatus, setSyncStatus] = useState<ViewerSyncStatus>({
    last_synced_at: null,
    last_error: null,
    local_project_count: 0,
  });
  const [isSyncing, setIsSyncing] = useState(false);
  const [syncMessage, setSyncMessage] = useState<string | null>(null);

  const [viewBandsData, setViewBandsData] = useState<{
    bandData: BandPlotData;
    fermiEnergy: number | null;
    calculationParameters?: Record<string, unknown> | null;
    calculationContext?: SavedBandsCalculationContext | null;
  } | null>(null);
  const [viewDosData, setViewDosData] = useState<{ dosData: ElectronicDOSData; fermiEnergy: number | null } | null>(null);
  const [viewWannierData, setViewWannierData] = useState<{
    result: any;
    fermiEnergy: number | null;
    overlayOptions: WannierBandOverlayOption[];
  } | null>(null);
  const [viewTransportData, setViewTransportData] = useState<{ data: TransportResult } | null>(null);
  const [viewPhononData, setViewPhononData] = useState<{ data: PhononData; mode: PhononViewMode } | null>(null);
  const [viewEpwData, setViewEpwData] = useState<EpwViewerPayload | null>(null);

  const activeHpcProfile = useMemo(
    () => hpcProfiles.find((profile) => profile.id === activeHpcProfileId) ?? null,
    [hpcProfiles, activeHpcProfileId],
  );
  const editingHpcProfile = useMemo(
    () => hpcProfiles.find((profile) => profile.id === editingHpcProfileId) ?? null,
    [hpcProfiles, editingHpcProfileId],
  );

  useEffect(() => {
    void initializeViewer();
  }, []);

  useEffect(() => {
    document.documentElement.setAttribute("data-current-view", currentView);
    return () => {
      document.documentElement.removeAttribute("data-current-view");
    };
  }, [currentView]);

  useEffect(() => {
    if (!activeHpcProfileId) return;
    void runSync(false);
    const handle = window.setInterval(() => {
      void runSync(false);
    }, 120_000);
    return () => window.clearInterval(handle);
  }, [activeHpcProfileId]);

  async function initializeViewer() {
    await refreshProfiles();
    await refreshSyncStatus();
  }

  async function refreshProfiles() {
    try {
      const [profiles, activeId] = await Promise.all([
        listHpcProfiles(),
        getActiveHpcProfileId(),
      ]);
      setHpcProfiles(profiles);
      setActiveHpcProfileId(activeId);
      return { profiles, activeId };
    } catch (e) {
      setSyncMessage(`Failed to load HPC profiles: ${e}`);
      setHpcProfiles([]);
      setActiveHpcProfileId(null);
      return { profiles: [], activeId: null };
    }
  }

  async function refreshSyncStatus() {
    try {
      const status = await invoke<ViewerSyncStatus>("viewer_get_sync_status");
      setSyncStatus(status);
    } catch (e) {
      setSyncStatus((current) => ({
        ...current,
        last_error: String(e),
      }));
    }
  }

  async function runSync(manual: boolean) {
    if (isSyncing || !activeHpcProfileId) return;
    setIsSyncing(true);
    if (manual) {
      setSyncMessage(null);
    }
    try {
      const result = await invoke<ViewerSyncResult>("viewer_sync_remote_library", {
        profileId: activeHpcProfileId,
      });
      if (manual) {
        setSyncMessage(
          `Synced ${result.total_projects} projects (${result.downloaded_projects} updated, ${result.removed_projects} removed).`,
        );
      }
    } catch (e) {
      setSyncMessage(`Sync failed: ${e}`);
    } finally {
      setIsSyncing(false);
      await refreshSyncStatus();
    }
  }

  const viewOnlyNoopScf = (
    _engineId: any,
    _cifId: string,
    _crystalData: any,
    _cifContent: string,
    _filename: string,
    _preset?: any,
    _presetLock?: boolean,
    _optimizedStructures?: any[],
  ) => undefined;
  const viewOnlyNoopWien2kContinue = (
    _cifId: string,
    _crystalData: any,
    _cifContent: string,
    _filename: string,
    _calculations: CalculationRun[],
    _calculationId: string,
  ) => undefined;
  const viewOnlyNoopEngineSetup = (_engineId: any, _cifId: string, _crystalData: any) => undefined;
  const viewOnlyNoopCalc = (_engineId: any, _cifId: string, _crystalData: any, _scfCalculations: CalculationRun[]) => undefined;
  const viewOnlyNoopTransport = (_engineId: any, _cifId: string, _crystalData: any, _wannierCalculations: CalculationRun[]) => undefined;

  const appChrome = (
    <HpcSetupWizard
      isOpen={showHpcSetupWizard}
      initialProfile={editingHpcProfile}
      onClose={() => {
        setShowHpcSetupWizard(false);
        setEditingHpcProfileId(null);
      }}
      onSaved={(profile) => {
        setEditingHpcProfileId(profile.id);
        setShowHpcSetupWizard(false);
        void refreshProfiles();
        void runSync(true);
      }}
    />
  );

  if (currentView === "project-browser") {
    return (
      <>
        <ProjectBrowser
          readOnly
          initialActiveFolderId={projectBrowserFolderId}
          onBack={() => setCurrentView("home")}
          onSelectProject={(projectId, folderId) => {
            setSelectedProjectId(projectId);
            setProjectBrowserFolderId(folderId);
            setCurrentView("project-dashboard");
          }}
          onOpenBandsMultiview={(initialCalculations) => {
            setBandsMultiviewInitialCalculations(initialCalculations);
            setBandsMultiviewReturnView("project-browser");
            setOpenBandsMultiviewInitialCalculationsInPreview(false);
            setCurrentView("bands-multiview");
          }}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "bands-multiview") {
    return (
      <>
        <BandsMultiview
          backLabel={bandsMultiviewReturnView === "project-dashboard" ? "Back to Project" : "Back to Projects"}
          onBack={() => {
            setBandsMultiviewInitialCalculations(undefined);
            setOpenBandsMultiviewInitialCalculationsInPreview(false);
            setCurrentView(bandsMultiviewReturnView);
          }}
          initialCalculations={bandsMultiviewInitialCalculations}
          openInitialCalculationsInPreview={openBandsMultiviewInitialCalculationsInPreview}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "project-dashboard" && selectedProjectId) {
    return (
      <>
        <ProjectDashboard
          projectId={selectedProjectId}
          readOnly
          onBack={() => {
            setCurrentView("project-browser");
            setSelectedProjectId(null);
          }}
          onDeleted={() => {
            setCurrentView("project-browser");
            setSelectedProjectId(null);
          }}
          onRunSCF={viewOnlyNoopScf}
          onContinueWien2kScf={viewOnlyNoopWien2kContinue}
          onRunEngineSetup={viewOnlyNoopEngineSetup}
          onRunBands={viewOnlyNoopCalc}
          onRunSoc={(_cifId, _crystalData, _scfCalculations) => undefined}
          onViewBands={(bandData, fermiEnergy, calculationParameters, calculationContext) => {
            setViewBandsData({
              bandData,
              fermiEnergy,
              calculationParameters,
              calculationContext,
            });
            setCurrentView("bands-viewer");
          }}
          onRunDos={viewOnlyNoopCalc}
          onViewDos={(dosData, fermiEnergy) => {
            setViewDosData({ dosData, fermiEnergy });
            setCurrentView("dos-viewer");
          }}
          onRunWannier={viewOnlyNoopCalc}
          onViewWannier={(wannierData, fermiEnergy, overlayOptions = []) => {
            setViewWannierData({ result: wannierData, fermiEnergy, overlayOptions });
            setCurrentView("wannier-viewer");
          }}
          onRunTransport={viewOnlyNoopTransport}
          onViewTransport={(transportData) => {
            setViewTransportData({ data: transportData });
            setCurrentView("transport-viewer");
          }}
          onRunFermiSurface={viewOnlyNoopCalc}
          onRunHubbardLrt={viewOnlyNoopCalc}
          onRunPhonons={viewOnlyNoopCalc}
          onRunEPW={viewOnlyNoopCalc}
          onViewPhonons={(phononData, viewMode) => {
            setViewPhononData({ data: phononData, mode: viewMode });
            setCurrentView("phonon-viewer");
          }}
          onViewEPW={(epwData, rawOutput) => {
            setViewEpwData({
              data: epwData,
              rawOutput: rawOutput ?? null,
            });
            setCurrentView("epw-viewer");
          }}
          onOpenBandsMultiview={(initialCalculations) => {
            setBandsMultiviewInitialCalculations(initialCalculations);
            setBandsMultiviewReturnView("project-dashboard");
            setOpenBandsMultiviewInitialCalculationsInPreview(true);
            setCurrentView("bands-multiview");
          }}
        />
        {appChrome}
      </>
    );
  }

  if (currentView === "bands-viewer" && viewBandsData) {
    return (
      <>
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
              scfFermiEnergy={viewBandsData.fermiEnergy ?? undefined}
              calculationParameters={viewBandsData.calculationParameters ?? null}
              viewerType="electronic"
            />
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "dos-viewer" && viewDosData) {
    return (
      <>
        <div className="bands-viewer-container">
          <div className="bands-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewDosData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>Electronic DOS</h2>
          </div>
          <div className="bands-viewer-content">
            <ElectronicDOSPlot
              data={{
                ...viewDosData.dosData,
                fermi_energy: viewDosData.dosData.fermi_energy ?? viewDosData.fermiEnergy,
              }}
            />
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "wannier-viewer" && viewWannierData) {
    const result = viewWannierData.result;
    const wannierIssues = getWannierQualityIssues(result, null, viewWannierData.fermiEnergy ?? null);
    return (
      <>
        <div className="bands-viewer-container">
          <div className="bands-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewWannierData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>Wannier90</h2>
          </div>
          <div className="bands-viewer-content bands-viewer-content-stacked">
            <div className="bands-viewer-plot-region">
              <BandPlot
                data={result.band_data}
                scfFermiEnergy={viewWannierData.fermiEnergy ?? undefined}
                viewerType="electronic"
                comparisonOptions={viewWannierData.overlayOptions}
                comparisonTitle="Saved Band Overlay"
                comparisonNoneLabel="No overlay"
              />
            </div>
            <div className="bands-viewer-details-region">
              {wannierIssues.length > 0 && (
                <div className="warning-banner">
                  {wannierIssues.map((issue) => issue.message).join(" ")}
                </div>
              )}
              <div className="details-grid">
              <div className="detail-item">
                <label>seedname</label>
                <span>{result.seedname || "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>num_wann</label>
                <span>{result.num_wann ?? "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>num_bands</label>
                <span>{result.num_bands ?? "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>Total Spread</label>
                <span>{result.total_spread != null ? `${Number(result.total_spread).toFixed(6)} A^2` : "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>Converged</label>
                <span>{result.convergence?.converged ? "Yes" : "No"}</span>
              </div>
              <div className="detail-item">
                <label>Iterations</label>
                <span>{result.convergence?.iterations ?? "N/A"}</span>
              </div>
              <div className="detail-item">
                <label>Minimization</label>
                <span>{formatWannierConvergenceFlag(result.convergence?.minimization_converged)}</span>
              </div>
              <div className="detail-item">
                <label>Disentanglement</label>
                <span>{formatWannierConvergenceFlag(result.convergence?.disentanglement_converged)}</span>
              </div>
            </div>
            </div>
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "transport-viewer" && viewTransportData) {
    return (
      <>
        <div className="bands-viewer-container transport-viewer-container">
          <div className="bands-viewer-header transport-viewer-header">
            <button
              className="back-button"
              onClick={() => {
                setCurrentView("project-dashboard");
                setViewTransportData(null);
              }}
            >
              ← Back to Dashboard
            </button>
            <h2>BoltzWann Transport</h2>
          </div>
          <div className="bands-viewer-content transport-viewer-content">
            <TransportPlot data={viewTransportData.data} />
          </div>
        </div>
        {appChrome}
      </>
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

    return (
      <>
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
          </div>
          <div className="phonon-viewer-content">
            {showingBands && phononBandData ? (
              <BandPlot
                data={phononBandData}
                showFermiLevel={false}
                yAxisLabel="Frequency (cm^-1)"
                pointLabel="Mode"
                valueLabel="Frequency"
                valueUnit="cm^-1"
                valueDecimals={1}
                primaryCountLabel="modes"
                secondaryCountLabel="q-points"
                scrollHint="Scroll: zoom frequency | Shift+Scroll: pan"
                yClampRange={null}
                viewerType="phonon"
              />
            ) : showingDos && hasDos ? (
              <PhononDOSPlot
                data={phononData.dos_data}
                width={Math.min(500, windowSize.width - 80)}
                height={plotHeight}
              />
            ) : (
              <p>{showingBands ? "No phonon dispersion data available" : "No phonon DOS data available"}</p>
            )}
          </div>
        </div>
        {appChrome}
      </>
    );
  }

  if (currentView === "epw-viewer" && viewEpwData) {
    return (
      <>
        <EpwViewer
          payload={viewEpwData}
          onBack={() => {
            setCurrentView("project-dashboard");
            setViewEpwData(null);
          }}
        />
        {appChrome}
      </>
    );
  }

  return (
    <>
      <main className="container">
        <header className="header">
          <h1>QCortado Viewer</h1>
          <p className="subtitle">Read-only access to synced HPC project results</p>
        </header>

        <section className="config-section">
          <h2>Viewer Sync</h2>
          <div className="status-row">
            <label>Execution Mode:</label>
            <span className="status ready">HPC</span>
          </div>
          <div className="status-row">
            <label>Active Profile:</label>
            <span className={`status ${activeHpcProfile ? "ready" : "pending"}`}>
              {activeHpcProfile ? `${activeHpcProfile.name} (${activeHpcProfile.username}@${activeHpcProfile.host})` : "Not configured"}
            </span>
          </div>
          <div className="status-row">
            <label>Last Sync:</label>
            <span className={`status ${syncStatus.last_error ? "error" : syncStatus.last_synced_at ? "ready" : "pending"}`}>
              {formatSyncTime(syncStatus.last_synced_at)}
            </span>
          </div>
          <div className="status-row">
            <label>Local Projects:</label>
            <span className="status">{syncStatus.local_project_count}</span>
          </div>

          {syncStatus.last_error && <div className="error">{syncStatus.last_error}</div>}
          {syncMessage && <div className="success-banner">{syncMessage}</div>}

          <div className="action-grid">
            <button
              className="action-btn"
              onClick={() => {
                setEditingHpcProfileId(activeHpcProfileId);
                setShowHpcSetupWizard(true);
              }}
            >
              <span className="action-icon">⚙</span>
              <span className="action-label">{activeHpcProfile ? "Edit HPC Profile" : "Configure HPC Profile"}</span>
              <span className="action-hint">Andromeda access and shared paths</span>
            </button>

            <button
              className="action-btn"
              onClick={() => {
                void runSync(true);
              }}
              disabled={!activeHpcProfileId || isSyncing}
            >
              <span className="action-icon">{isSyncing ? "…" : "↻"}</span>
              <span className="action-label">{isSyncing ? "Syncing..." : "Sync Now"}</span>
              <span className="action-hint">Pull latest published snapshot</span>
            </button>

            <button
              className="action-btn"
              onClick={() => setCurrentView("project-browser")}
              disabled={syncStatus.local_project_count === 0}
            >
              <span className="action-icon">{syncStatus.local_project_count}</span>
              <span className="action-label">Open Projects</span>
              <span className="action-hint">
                {syncStatus.local_project_count > 0 ? "Browse results" : "Sync to load projects"}
              </span>
            </button>
          </div>
        </section>
      </main>
      {appChrome}
    </>
  );
}

function ViewerApp() {
  return (
    <ThemeProvider>
      <TaskProvider>
        <ViewerAppInner />
      </TaskProvider>
    </ThemeProvider>
  );
}

export default ViewerApp;
