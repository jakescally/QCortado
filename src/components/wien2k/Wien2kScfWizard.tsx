import { useEffect, useMemo, useRef, useState } from "react";
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
import { LiveOutputPanel } from "../LiveOutputPanel";
import { Wien2kFieldLabel } from "./Wien2kFieldLabel";

interface Wien2kScfWizardProps {
  projectId: string;
  cifId: string;
  calculations: Wien2kStructureSourceRecord[];
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

export function Wien2kScfWizard({
  projectId,
  cifId,
  calculations,
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

  useEffect(() => () => outputUnlistenRef.current?.(), []);

  async function attachOutputListener(sessionId: string) {
    outputUnlistenRef.current?.();
    outputUnlistenRef.current = await listen<string>(`wien2k-scf-output:${sessionId}`, (event) => {
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
    try {
      const next = await runWien2kScfSession(
        session.sessionId,
        effectiveRunSettings,
        continuation,
        continuation ? result?.calculationId ?? null : null,
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
          <div className="wien2k-summary">
            <h3>Accepted Structure Source</h3>
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

          {sites.length > 0 && (
            <>
              <h3 className="wien2k-site-heading">Accepted muffin-tin radii</h3>
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
            </>
          )}

          <h3 className="wien2k-site-heading">Initialization</h3>
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

          {initialized && (
            <>
              <h3 className="wien2k-site-heading">SCF Cycle</h3>
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
          )}

          {result && (
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
          )}
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
