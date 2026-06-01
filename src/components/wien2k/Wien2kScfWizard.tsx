import { useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import {
  DEFAULT_WIEN2K_INITIALIZATION_SETTINGS,
  DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
  WIEN2K_EXCHANGE_CORRELATION_OPTIONS,
  discardWien2kScfSession,
  initializeWien2kScfSession,
  listWien2kStructureSources,
  runWien2kScfSession,
  startWien2kScfSession,
  validateWien2kInitializationSettings,
  validateWien2kScfRunSettings,
} from "../../lib/engines/wien2k";
import type {
  NormalizedScfSummary,
  Wien2kInitializationSettings,
  Wien2kScfExecutionResult,
  Wien2kScfRunSettings,
  Wien2kScfSession,
  Wien2kStructureSourceRecord,
} from "../../lib/engines/wien2k";
import { defaultCpuResources } from "../../lib/hpcConfig";
import type { HpcProfile, SlurmResourceRequest } from "../../lib/types";
import { useViewportScrollLock } from "../../lib/useViewportScrollLock";
import { ElapsedTimer } from "../ElapsedTimer";
import { HpcRunSettings } from "../HpcRunSettings";
import { LiveOutputPanel } from "../LiveOutputPanel";
import { ProgressBar } from "../ProgressBar";
import { RemoteUtilizationPanel } from "../RemoteUtilizationPanel";
import { Wien2kFieldLabel } from "./Wien2kFieldLabel";

interface Wien2kScfWizardProps {
  projectId: string;
  cifId: string;
  calculations: Wien2kStructureSourceRecord[];
  activeHpcProfile?: HpcProfile | null;
  onBack: () => void;
  onSaved: () => void;
}

type ScfWizardStep = "source" | "initialize" | "scf" | "results";
const SCF_STEPS: Array<{ id: ScfWizardStep; label: string }> = [
  { id: "source", label: "Source" },
  { id: "initialize", label: "Initialize" },
  { id: "scf", label: "SCF" },
  { id: "results", label: "Results" },
];

function convergenceLabel(summary: NormalizedScfSummary): string {
  switch (summary.convergence) {
    case "converged": return "Converged";
    case "not_converged": return "Not converged";
    case "failed": return "Failed";
    case "cancelled": return "Cancelled";
    default: return "Unknown";
  }
}

function numberField(value: string, fallback: number): number {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : fallback;
}

type SectionKey = "source" | "radii" | "initialization" | "scf" | "hpc";

function cloneWien2kCpuResources(profile: HpcProfile | null | undefined): SlurmResourceRequest {
  const source = profile?.default_cpu_resources ?? defaultCpuResources();
  const openMpThreads = Math.max(1, source.ntasks ?? 1) * Math.max(1, source.cpus_per_task ?? 1);
  return {
    ...source,
    resource_type: "cpu",
    nodes: 1,
    ntasks: 1,
    cpus_per_task: openMpThreads,
    gpus: 0,
    additional_sbatch: [...(source.additional_sbatch ?? [])],
  };
}

function shellQuote(value: string): string {
  return `'${value.replace(/'/g, `'\"'\"'`)}'`;
}

function moduleSetupLines(profile: HpcProfile | null | undefined): string[] {
  if (profile?.wien2k_path_mode !== "module") {
    return profile?.remote_wien2k_install_root
      ? [`export WIENROOT=${shellQuote(profile.remote_wien2k_install_root)}`, `export PATH="$WIENROOT:$PATH"`]
      : [];
  }
  const lines: string[] = [];
  const moduleUse = (profile.wien2k_module_use ?? "").trim();
  const moduleLoad = (profile.wien2k_module_load ?? "").trim();
  if (moduleUse) lines.push(`module use ${shellQuote(moduleUse)}`);
  if (moduleLoad) lines.push(`module load ${shellQuote(moduleLoad)}`);
  return lines;
}

function buildWien2kMachinesPreview(resources: SlurmResourceRequest): string {
  const cpusPerTask = Math.max(1, resources.cpus_per_task ?? 1);
  return `rm -f .machines .processes && export OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK:-${cpusPerTask}}"`;
}

function buildRunPreviewCommand(profile: HpcProfile | null | undefined, settings: Wien2kScfRunSettings): string {
  const program = settings.spinMode === "spin_polarized" ? "runsp_lapw" : "run_lapw";
  const executable = profile?.wien2k_path_mode === "module" ? program : `"$WIENROOT/${program}"`;
  const args = [
    "-ec", String(settings.energyConvergenceRy),
    "-cc", String(settings.chargeConvergence),
    "-i", String(settings.maxIterations),
    ...(settings.forceConvergenceMryBohr != null ? ["-fc", String(settings.forceConvergenceMryBohr)] : []),
  ];
  return `${executable} ${args.map(shellQuote).join(" ")}`;
}

function parseRemoteJobId(line: string): string | null {
  const match = line.match(/Scheduled (?:utility )?job submitted:\s*(\S+)/);
  return match?.[1] ?? null;
}

function parseRemoteNode(line: string): string | null {
  const match = line.match(/^\[QCortado\] Host:\s*(.+)$/);
  return match?.[1]?.trim() ?? null;
}

export function Wien2kScfWizard({
  projectId,
  cifId,
  calculations,
  activeHpcProfile = null,
  onBack,
  onSaved,
}: Wien2kScfWizardProps) {
  const sources = useMemo(() => listWien2kStructureSources(calculations), [calculations]);
  const [sourceId, setSourceId] = useState(() => sources[0]?.id ?? "");
  const [initialization, setInitialization] = useState<Wien2kInitializationSettings>(
    DEFAULT_WIEN2K_INITIALIZATION_SETTINGS,
  );
  const [runSettings, setRunSettings] = useState<Wien2kScfRunSettings>(
    DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
  );
  const [session, setSession] = useState<Wien2kScfSession | null>(null);
  const [result, setResult] = useState<Wien2kScfExecutionResult | null>(null);
  const [outputLines, setOutputLines] = useState<string[]>([]);
  const [isInitializing, setIsInitializing] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [expandedSections, setExpandedSections] = useState<Record<SectionKey, boolean>>({
    source: true,
    radii: true,
    initialization: true,
    scf: false,
    hpc: false,
  });
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(() => cloneWien2kCpuResources(activeHpcProfile));
  const [scfRunStarted, setScfRunStarted] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string | null>(null);
  const [remoteJobId, setRemoteJobId] = useState<string | null>(null);
  const [remoteNode, setRemoteNode] = useState<string | null>(null);
  const outputUnlistenRef = useRef<UnlistenFn | null>(null);
  const selectedSource = sources.find((source) => source.id === sourceId) ?? sources[0] ?? null;
  const initError = validateWien2kInitializationSettings(initialization);
  const effectiveRunSettings = { ...runSettings, spinMode: initialization.spinMode };
  const runError = validateWien2kScfRunSettings(effectiveRunSettings);
  const initialized = session?.phase === "initialized" || session?.phase === "scf_complete" || session?.phase === "failed";
  const currentStep: ScfWizardStep = !session
    ? "source"
    : !initialized
      ? "initialize"
      : result
        ? "results"
        : "scf";
  const currentStepIndex = SCF_STEPS.findIndex((step) => step.id === currentStep);
  const sourceParameters = selectedSource?.parameters ?? {};
  const caseName = typeof sourceParameters.case_name === "string" ? sourceParameters.case_name : "case";
  const finalSummary = sourceParameters.final_structure_summary as {
    spacegroupSymbol?: string | null;
    spacegroupNumber?: number | null;
  } | undefined;
  const sites = Array.isArray(sourceParameters.sites)
    ? sourceParameters.sites as Array<{ symbol?: string; rmt?: number; positions?: unknown[] }>
    : [];
  const hpcCommandLines = useMemo(
    () => [
      "cd \"$SLURM_SUBMIT_DIR\"",
      ...moduleSetupLines(activeHpcProfile),
      buildWien2kMachinesPreview(hpcResources),
      buildRunPreviewCommand(activeHpcProfile, effectiveRunSettings),
    ],
    [activeHpcProfile, effectiveRunSettings, hpcResources],
  );

  useViewportScrollLock(scfRunStarted);

  useEffect(() => () => outputUnlistenRef.current?.(), []);

  useEffect(() => {
    setHpcResources(cloneWien2kCpuResources(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.updated_at]);

  function toggleSection(section: SectionKey) {
    if ((section === "scf" || section === "hpc") && !initialized) return;
    setExpandedSections((current) => ({ ...current, [section]: !current[section] }));
  }

  async function attachOutputListener(sessionId: string) {
    outputUnlistenRef.current?.();
    outputUnlistenRef.current = await listen<string>(`wien2k-scf-output:${sessionId}`, (event) => {
      const jobId = parseRemoteJobId(event.payload);
      if (jobId) setRemoteJobId(jobId);
      const node = parseRemoteNode(event.payload);
      if (node) setRemoteNode(node);
      setOutputLines((current) => [...current, event.payload].slice(-1500));
    });
  }

  async function beginInitialization() {
    if (!selectedSource || initError) return;
    setIsInitializing(true);
    setError(null);
    setResult(null);
    try {
      let activeSession = session;
      if (!activeSession) {
        activeSession = await startWien2kScfSession(projectId, cifId, selectedSource.id);
        setSession(activeSession);
        setOutputLines([`Staged accepted ${activeSession.caseName}.struct in ${activeSession.remoteCaseDir}.`]);
        await attachOutputListener(activeSession.sessionId);
      }
      const initializedResult = await initializeWien2kScfSession(activeSession.sessionId, initialization);
      setSession((current) => current ? {
        ...current,
        phase: initializedResult.phase,
        initialization,
        artifacts: { ...current.artifacts, ...initializedResult.artifacts },
      } : current);
      setOutputLines((current) => [
        ...current,
        `[native initialization validated: ${Object.keys(initializedResult.artifacts).length} text artifacts retrieved]`,
      ].slice(-1500));
      setExpandedSections({
        source: false,
        radii: false,
        initialization: false,
        scf: true,
        hpc: true,
      });
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsInitializing(false);
    }
  }

  async function submitScf(continuation: boolean) {
    if (!session || runError) return;
    setIsRunning(true);
    setError(null);
    setScfRunStarted(true);
    setCalcStartTime(new Date().toISOString());
    setRemoteJobId(null);
    setRemoteNode(null);
    try {
      const next = await runWien2kScfSession(
        session.sessionId,
        effectiveRunSettings,
        continuation,
        continuation ? result?.calculationId ?? null : null,
        hpcResources,
      );
      setResult(next);
      setSession((current) => current ? {
        ...current,
        phase: next.phase,
        latestRun: effectiveRunSettings,
        latestCalculationId: next.calculationId,
      } : current);
      setOutputLines((current) => [
        ...current,
        `[saved calculation ${next.calculationId}: ${convergenceLabel(next.summary)}]`,
      ].slice(-1500));
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsRunning(false);
    }
  }

  async function resetInitialization() {
    if (!session) return;
    setError(null);
    try {
      await discardWien2kScfSession(session.sessionId);
      outputUnlistenRef.current?.();
      outputUnlistenRef.current = null;
      setSession(null);
      setResult(null);
      setOutputLines([]);
      setScfRunStarted(false);
      setCalcStartTime(null);
      setRemoteJobId(null);
      setRemoteNode(null);
      setExpandedSections({
        source: true,
        radii: true,
        initialization: true,
        scf: false,
        hpc: false,
      });
    } catch (reason) {
      setError(String(reason));
    }
  }

  async function leaveWizard(destination: "back" | "saved") {
    setError(null);
    try {
      if (session) {
        await discardWien2kScfSession(session.sessionId);
      }
      outputUnlistenRef.current?.();
      outputUnlistenRef.current = null;
      if (destination === "saved") onSaved();
      else onBack();
    } catch (reason) {
      setError(String(reason));
    }
  }

  function renderSection(
    key: SectionKey,
    title: string,
    body: ReactNode,
    options: { locked?: boolean; status?: string } = {},
  ) {
    const expanded = expandedSections[key] && !options.locked;
    return (
      <section className={`config-section collapsible wien2k-scf-section${options.locked ? " is-locked" : ""}`}>
        <h4
          className="section-header"
          onClick={() => toggleSection(key)}
          aria-expanded={expanded}
        >
          <span className={`collapse-icon ${expanded ? "expanded" : ""}`}>▶</span>
          {title}
          {options.status && <span className="wien2k-section-status">{options.status}</span>}
        </h4>
        {expanded && body}
      </section>
    );
  }

  function renderResultSummary() {
    if (!result) return null;
    return (
      <div className="wien2k-summary wien2k-scf-results">
        <h3>SCF Result: {convergenceLabel(result.summary)}</h3>
        <p>
          {result.summary.totalEnergy ? `E = ${result.summary.totalEnergy.value.toFixed(8)} ${result.summary.totalEnergy.unit}; ` : ""}
          {result.summary.scfSteps != null ? `${result.summary.scfSteps} iterations` : "Iterations unavailable"}
        </p>
        {result.summary.fermiEnergyEv != null && <p>Fermi energy: {result.summary.fermiEnergyEv.toFixed(6)} eV</p>}
        {result.summary.totalMagnetization != null && <p>Total magnetization: {result.summary.totalMagnetization.toFixed(6)}</p>}
        {result.diagnostics.map((diagnostic) => <p key={diagnostic} className="wien2k-validation">{diagnostic}</p>)}
      </div>
    );
  }

  function renderRunStep() {
    const progressStatus: "idle" | "running" | "error" | "complete" =
      error ? "error" : result ? "complete" : isRunning ? "running" : "idle";
    return (
      <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard wizard-step-run">
        <div className="wizard-header">
          <button className="back-btn" type="button" disabled={isRunning} onClick={() => void leaveWizard("back")}>
            Back
          </button>
          <h2>WIEN2k SCF</h2>
          <div className="step-indicator">
            {SCF_STEPS.map((step, index) => (
              <span
                key={step.id}
                className={index === currentStepIndex ? "active" : index < currentStepIndex ? "completed" : ""}
              >
                {index + 1}. {step.label}
              </span>
            ))}
          </div>
        </div>
        {error && <div className="error-banner">{error}</div>}
        <div className="wizard-content">
          <div className="wizard-step run-step run-step-focused scf-run-step">
            <div className="run-step-headline">
              <h3>{isRunning ? "Running WIEN2k SCF" : "WIEN2k SCF Output"}</h3>
              <span className={`run-step-status-pill ${isRunning ? "running" : error ? "error" : "idle"}`}>
                {isRunning ? "Live output" : error ? "Run failed" : "Output"}
              </span>
            </div>

            <div className="run-status-rail scf-run-status">
              <ProgressBar
                status={progressStatus}
                percent={null}
                phase={isRunning ? "SCF cycle" : result ? convergenceLabel(result.summary) : "SCF output"}
                detail={remoteJobId ? `Slurm job ${remoteJobId}` : "Waiting for remote allocation"}
                compact
              />
              <div className="run-status-meta">
                <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />
              </div>
            </div>

            <div className="run-layout run-layout-hpc-telemetry">
              <LiveOutputPanel
                title={isRunning ? "Running..." : "Output"}
                output={outputLines.join("\n")}
                placeholder="Starting WIEN2k SCF..."
                totalLineCount={outputLines.length}
                visibleLineCount={outputLines.length}
              />
              <RemoteUtilizationPanel
                enabled={isRunning}
                profileId={activeHpcProfile?.id ?? null}
                remoteJobId={remoteJobId}
                remoteNode={remoteNode}
                resourceType="cpu"
              />
            </div>
            {renderResultSummary()}
            <div className="run-actions">
              <button type="button" disabled={isRunning} onClick={() => setScfRunStarted(false)}>
                Back to SCF Setup
              </button>
              {result?.summary.convergence === "not_converged" && (
                <button type="button" className="primary-button" disabled={Boolean(runError) || isRunning} onClick={() => void submitScf(true)}>
                  {isRunning ? "Continuing..." : "Continue SCF"}
                </button>
              )}
              {result && (
                <button type="button" className="secondary-button" disabled={isRunning} onClick={() => void leaveWizard("saved")}>
                  Return to Project
                </button>
              )}
            </div>
          </div>
        </div>
      </div>
    );
  }

  if (scfRunStarted) {
    return renderRunStep();
  }

  return (
    <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
      <div className="wizard-header">
        <button className="back-btn" type="button" disabled={isInitializing || isRunning} onClick={() => void leaveWizard("back")}>
          Back
        </button>
        <h2>WIEN2k SCF</h2>
        <div className="step-indicator">
          {SCF_STEPS.map((step, index) => (
            <span
              key={step.id}
              className={index === currentStepIndex ? "active" : index < currentStepIndex ? "completed" : ""}
            >
              {index + 1}. {step.label}
            </span>
          ))}
        </div>
      </div>
      {error && <div className="error-banner">{error}</div>}
      <div className="wien2k-structure-content">
        <section className="wien2k-structure-controls">
          {renderSection("source", "Accepted Structure Source", (
            <div className="wien2k-summary">
              <label>
                <Wien2kFieldLabel tooltip="Selects the accepted WIEN2k structure, symmetry, and muffin-tin radii used as the SCF basis. Recommended choice: the latest reviewed structure for the intended material; change geometry or RMT in the Structure workflow.">
                  Saved case.struct
                </Wien2kFieldLabel>
                <select value={sourceId} disabled={Boolean(session)} onChange={(event) => setSourceId(event.target.value)}>
                  {sources.map((source) => (
                    <option key={source.id} value={source.id}>
                      {String(source.parameters?.case_name ?? "case")}.struct
                      {source.completed_at ? ` - ${new Date(source.completed_at).toLocaleString()}` : ""}
                    </option>
                  ))}
                </select>
              </label>
              {selectedSource ? (
                <p>
                  {caseName}.struct;{" "}
                  {finalSummary?.spacegroupSymbol
                    ? `${finalSummary.spacegroupSymbol} (#${finalSummary.spacegroupNumber ?? "?"})`
                    : "native symmetry accepted"}.
                  RMT values are fixed by this structure source.
                </p>
              ) : (
                <p className="wien2k-validation">No accepted WIEN2k Structure source is available.</p>
              )}
            </div>
          ), { status: session ? "Locked" : undefined })}

          {renderSection("radii", "Accepted Muffin-Tin Radii", sites.length > 0 ? (
            <table className="wien2k-site-table">
              <thead><tr><th>Site</th><th>Multiplicity</th><th>RMT</th></tr></thead>
              <tbody>
                {sites.map((site, index) => (
                  <tr key={`${site.symbol ?? "site"}-${index}`}>
                    <td>{site.symbol ?? `Site ${index + 1}`}</td>
                    <td>{Array.isArray(site.positions) ? site.positions.length : "-"}</td>
                    <td>{typeof site.rmt === "number" ? site.rmt.toFixed(4) : "-"}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          ) : (
            <p className="wien2k-validation">No accepted muffin-tin radii were saved with this structure source.</p>
          ), { status: initialized ? "Locked" : undefined })}

          {renderSection("initialization", "Initialization", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="RMT(min) times KMAX controls the LAPW plane-wave basis size; larger values improve basis accuracy while increasing memory and runtime. Typical starting range: 6 to 9, followed by convergence testing.">
                    RKMAX
                  </Wien2kFieldLabel>
                  <input type="number" min="0.1" step="0.1" disabled={initialized} value={initialization.rkmax} onChange={(event) => setInitialization((current) => ({ ...current, rkmax: numberField(event.target.value, current.rkmax) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Sets the reciprocal-space cutoff for representing density and potential components; raising it resolves sharper density features at higher cost. Typical starting range: 10 to 16.">
                    GMAX
                  </Wien2kFieldLabel>
                  <input type="number" min="0.1" step="0.1" disabled={initialized} value={initialization.gmax} onChange={(event) => setInitialization((current) => ({ ...current, gmax: numberField(event.target.value, current.gmax) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Maximum angular momentum in the muffin-tin expansion; higher values capture more aspherical character but add work. Typical starting range: 8 to 12.">
                    LMAX
                  </Wien2kFieldLabel>
                  <input type="number" min="1" step="1" disabled={initialized} value={initialization.lmax} onChange={(event) => setInitialization((current) => ({ ...current, lmax: numberField(event.target.value, current.lmax) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Chooses the exchange-correlation approximation used to initialize and run the density. WIEN2k initialization provides LDA, PBE, WC, and PBEsol; PBE is a common general-purpose default, and comparisons should keep XC fixed.">
                    XC functional
                  </Wien2kFieldLabel>
                  <select disabled={initialized} value={initialization.exchangeCorrelation} onChange={(event) => setInitialization((current) => ({ ...current, exchangeCorrelation: Number(event.target.value) }))}>
                    {WIEN2K_EXCHANGE_CORRELATION_OPTIONS.map((option) => (
                      <option key={option.value} value={option.value}>{option.label}</option>
                    ))}
                  </select>
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Separates core and valence states during LSTART, affecting which electrons participate in the variational SCF basis. Typical starting range: -8 to -4 Ry; check for core leakage warnings.">
                    LSTART cutoff (Ry)
                  </Wien2kFieldLabel>
                  <input type="number" step="0.1" disabled={initialized} value={initialization.lstartEnergyCutoffRy} onChange={(event) => setInitialization((current) => ({ ...current, lstartEnergyCutoffRy: numberField(event.target.value, current.lstartEnergyCutoffRy) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Enables collinear spin-resolved densities for magnetic calculations, adding separate spin channels and computational work. Recommended choice: non-spin-polarized unless magnetism is physically expected.">
                    Spin mode
                  </Wien2kFieldLabel>
                  <select disabled={initialized} value={initialization.spinMode} onChange={(event) => setInitialization((current) => ({ ...current, spinMode: event.target.value as Wien2kInitializationSettings["spinMode"] }))}>
                    <option value="non_spin_polarized">Non-spin-polarized</option>
                    <option value="spin_polarized">Spin-polarized</option>
                  </select>
                </label>
              </div>
              <div className="wien2k-control-grid wien2k-kmesh-grid">
                {initialization.kMesh.map((value, index) => (
                  <label key={index}>
                    <Wien2kFieldLabel tooltip={`Number of reciprocal-space sampling divisions along lattice direction ${index + 1}; denser meshes improve Brillouin-zone integration while multiplying runtime. Typical starting range: 4 to 12 per direction, scaled for cell shape.`}>
                      k{index + 1}
                    </Wien2kFieldLabel>
                    <input type="number" min="1" step="1" disabled={initialized} value={value} onChange={(event) => setInitialization((current) => {
                      const next = [...current.kMesh] as [number, number, number];
                      next[index] = numberField(event.target.value, value);
                      return { ...current, kMesh: next };
                    })} />
                  </label>
                ))}
              </div>
              {initError && <p className="wien2k-validation">{initError}</p>}
            </>
          ), { status: initialized ? "Locked" : undefined })}

          {renderSection("scf", "SCF Cycle", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="Stops SCF when successive total-energy changes are sufficiently small; tighter tolerances improve reproducibility but may require more iterations. Typical range: 1e-5 to 1e-4 Ry.">
                    Energy convergence (Ry)
                  </Wien2kFieldLabel>
                  <input type="number" min="0" step="0.00001" value={runSettings.energyConvergenceRy} onChange={(event) => setRunSettings((current) => ({ ...current, energyConvergenceRy: numberField(event.target.value, current.energyConvergenceRy) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Stops SCF when the charge-density residual is small, directly measuring density self-consistency. Typical range: 1e-4 to 1e-3 electrons.">
                    Charge convergence
                  </Wien2kFieldLabel>
                  <input type="number" min="0" step="0.00001" value={runSettings.chargeConvergence} onChange={(event) => setRunSettings((current) => ({ ...current, chargeConvergence: numberField(event.target.value, current.chargeConvergence) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Caps SCF update cycles before the retained case is reported as unconverged and may be continued. Typical range: 30 to 100 iterations.">
                    Maximum iterations
                  </Wien2kFieldLabel>
                  <input type="number" min="1" step="1" value={runSettings.maxIterations} onChange={(event) => setRunSettings((current) => ({ ...current, maxIterations: numberField(event.target.value, current.maxIterations) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Optionally requires Hellmann-Feynman forces to fall below a threshold in addition to electronic convergence. Typical range when used: 1 to 5 mRy/Bohr.">
                    Force convergence (mRy/Bohr)
                  </Wien2kFieldLabel>
                  <input type="number" min="0" step="0.1" placeholder="Optional" value={runSettings.forceConvergenceMryBohr ?? ""} onChange={(event) => setRunSettings((current) => ({ ...current, forceConvergenceMryBohr: event.target.value.trim() ? numberField(event.target.value, 0) : null }))} />
                </label>
              </div>
              {runError && <p className="wien2k-validation">{runError}</p>}
            </>
          ), { locked: !initialized, status: initialized ? undefined : "Locked - run initialization first" })}

          {renderSection("hpc", "HPC Run Settings", (
              <HpcRunSettings
                profileId={activeHpcProfile?.id ?? null}
                profileName={activeHpcProfile?.name ?? "Andromeda"}
                taskKind="wien2k-scf"
                commandLines={hpcCommandLines}
                resources={hpcResources}
                resourceMode="cpu_only"
                defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
                defaultGpuResources={null}
                resourceModeMessage="WIEN2k SCF currently uses one CPU Slurm task with OpenMP threads."
                lockResourceTypeWhenSingleMode={false}
                maxTasks={{
                  value: 1,
                  reason: "WIEN2k k-point .machines mode requires passwordless compute-node SSH and is not enabled by this wizard yet",
                }}
                onResourcesChange={(next) => setHpcResources({
                  ...next,
                  resource_type: "cpu",
                  nodes: 1,
                  ntasks: 1,
                  gpus: 0,
                })}
                disabled={isRunning}
              />
          ), { locked: !initialized, status: initialized ? undefined : "Locked - run initialization first" })}
        </section>
        <section className="wien2k-output-column">
          <LiveOutputPanel
            title="WIEN2k initialization and SCF"
            output={outputLines.join("\n")}
            totalLineCount={outputLines.length}
            visibleLineCount={outputLines.length}
            panelClassName="output-panel wien2k-inline-output"
            outputClassName="output-text wien2k-inline-output-text"
          />
        </section>
      </div>
      <div className="step-actions">
        {session && !result && initialized && (
          <button className="secondary-button" type="button" disabled={isRunning} onClick={() => void resetInitialization()}>
            Restart Initialization
          </button>
        )}
        {!initialized && (
          <button type="button" className="primary-button" disabled={!selectedSource || Boolean(initError) || isInitializing} onClick={() => void beginInitialization()}>
            {isInitializing ? "Initializing..." : "Run Initialization"}
          </button>
        )}
        {initialized && !result && (
          <button type="button" className="primary-button" disabled={Boolean(runError) || isRunning} onClick={() => void submitScf(false)}>
            {isRunning ? "Running SCF..." : "Run SCF"}
          </button>
        )}
        {result?.summary.convergence === "not_converged" && (
          <button type="button" className="primary-button" disabled={Boolean(runError) || isRunning} onClick={() => void submitScf(true)}>
            {isRunning ? "Continuing..." : "Continue SCF"}
          </button>
        )}
        {result && (
          <button type="button" className="secondary-button" disabled={isRunning} onClick={() => void leaveWizard("saved")}>
            Return to Project
          </button>
        )}
      </div>
    </div>
  );
}
