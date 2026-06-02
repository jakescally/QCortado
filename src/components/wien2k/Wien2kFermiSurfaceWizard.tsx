import { useEffect, useMemo, useState, type ReactNode } from "react";
import { invoke } from "@tauri-apps/api/core";
import { ElapsedTimer } from "../ElapsedTimer";
import { HpcRunSettings } from "../HpcRunSettings";
import { InfoTooltip } from "../InfoTooltip";
import { LiveOutputPanel } from "../LiveOutputPanel";
import { ProgressBar } from "../ProgressBar";
import { RemoteUtilizationPanel } from "../RemoteUtilizationPanel";
import { Wien2kFieldLabel } from "./Wien2kFieldLabel";
import { defaultCpuResources } from "../../lib/hpcConfig";
import { formatCalculationSourceLabel, getCalculationName } from "../../lib/calculationNames";
import { useTaskContext } from "../../lib/TaskContext";
import type { CrystalData, HpcProfile, SlurmResourceRequest } from "../../lib/types";
import type {
  Wien2kBandsSpinChannel,
  Wien2kFermiSurfaceResult,
  Wien2kFermiSurfaceSettings,
} from "../../lib/engines/wien2k";
import { useViewportScrollLock } from "../../lib/useViewportScrollLock";

interface CalculationRun {
  id: string;
  name?: string | null;
  engine_id?: string | null;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    total_energy: number | null;
    fermi_energy: number | null;
    raw_output?: string | null;
  } | null;
  scf_summary?: {
    convergence: "converged" | "not_converged" | "failed" | "cancelled" | "unknown";
    totalEnergy?: { value: number; unit: string } | null;
    fermiEnergyEv?: number | null;
    scfSteps?: number | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

interface Wien2kFermiSurfaceWizardProps {
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  activeHpcProfile?: HpcProfile | null;
  reconnectTaskId?: string;
  onBack: () => void;
}

type WizardStep = "source" | "parameters" | "run" | "results";
type SectionKey = "source" | "mesh" | "xcrysden" | "logging" | "hpc";

const STEPS: Array<{ id: WizardStep; label: string }> = [
  { id: "source", label: "Source" },
  { id: "parameters", label: "Parameters" },
  { id: "run", label: "Run" },
  { id: "results", label: "Results" },
];

function cloneResources(profile: HpcProfile | null | undefined): SlurmResourceRequest {
  const source = profile?.default_cpu_resources ?? defaultCpuResources();
  return {
    ...source,
    resource_type: "cpu",
    gpus: 0,
    additional_sbatch: [...(source.additional_sbatch ?? [])],
  };
}

function numberField(value: string, fallback: number): number {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : fallback;
}

function clampKMesh(value: number): number {
  if (!Number.isFinite(value)) return 1;
  return Math.min(500, Math.max(1, Math.round(value)));
}

function isConvergedWien2kScf(calc: CalculationRun): boolean {
  return calc.engine_id === "wien2k"
    && calc.calc_type === "scf"
    && calc.scf_summary?.convergence === "converged"
    && typeof calc.parameters?.remote_case_dir === "string";
}

function sourceFermi(calc: CalculationRun | null): number | null {
  if (!calc) return null;
  const summaryValue = Number(calc.scf_summary?.fermiEnergyEv);
  if (Number.isFinite(summaryValue)) return summaryValue;
  const resultValue = Number(calc.result?.fermi_energy);
  return Number.isFinite(resultValue) ? resultValue : null;
}

function sourceSpinMode(calc: CalculationRun | null): "non_spin_polarized" | "spin_polarized" {
  const value =
    calc?.parameters?.run?.spinMode
    ?? calc?.parameters?.run?.spin_mode
    ?? calc?.parameters?.initialization?.spinMode
    ?? calc?.parameters?.initialization?.spin_mode;
  return value === "spin_polarized" ? "spin_polarized" : "non_spin_polarized";
}

function parseRemoteJobId(line: string): string | null {
  const match = line.match(/Scheduled (?:utility )?job submitted:\s*(\S+)/);
  return match?.[1] ?? null;
}

function parseRemoteNode(line: string): string | null {
  const match = line.match(/^\[QCortado\]\s+Host:\s*(.+)$/);
  return match?.[1]?.trim() ?? null;
}

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value.toFixed(value >= 10 || unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
}

export function Wien2kFermiSurfaceWizard({
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  activeHpcProfile = null,
  reconnectTaskId,
  onBack,
}: Wien2kFermiSurfaceWizardProps) {
  const taskContext = useTaskContext();
  const sources = useMemo(() => scfCalculations.filter(isConvergedWien2kScf), [scfCalculations]);
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(sources[0] ?? null);
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [kMesh, setKMesh] = useState<[number, number, number]>(() => {
    const sourceMesh = sources[0]?.parameters?.initialization?.kMesh ?? sources[0]?.parameters?.initialization?.k_mesh;
    return Array.isArray(sourceMesh) && sourceMesh.length === 3
      ? [clampKMesh(Number(sourceMesh[0]) * 2), clampKMesh(Number(sourceMesh[1]) * 2), clampKMesh(Number(sourceMesh[2]) * 2)]
      : [24, 24, 24];
  });
  const [spinChannel, setSpinChannel] = useState<Wien2kBandsSpinChannel>("none");
  const [spinOrbit, setSpinOrbit] = useState(false);
  const [diagnosticLog, setDiagnosticLog] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(() => cloneResources(activeHpcProfile));
  const [expandedSections, setExpandedSections] = useState<Record<SectionKey, boolean>>({
    source: true,
    mesh: true,
    xcrysden: true,
    logging: true,
    hpc: true,
  });
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const [result, setResult] = useState<Wien2kFermiSurfaceResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [calcStartTime, setCalcStartTime] = useState<string | null>(null);
  const [remoteJobId, setRemoteJobId] = useState<string | null>(null);
  const [remoteNode, setRemoteNode] = useState<string | null>(null);
  const [isLaunching, setIsLaunching] = useState(false);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const taskResult = activeTask?.result as Wien2kFermiSurfaceResult | null | undefined;
  const displayedResult = result ?? taskResult ?? null;
  const runIsActive = activeTask?.status === "running";
  const runHasFailed = activeTask?.status === "failed" || activeTask?.status === "cancelled" || Boolean(error);
  const taskOutputText = activeTask ? activeTask.outputText : "";
  const taskOutputLines = activeTask ? activeTask.output : [];
  const taskOutputLineCount = activeTask ? activeTask.outputLineCount : 0;
  const currentStepIndex = STEPS.findIndex((entry) => entry.id === step);
  const selectedSpinMode = sourceSpinMode(selectedScf);
  const fermiEnergy = sourceFermi(selectedScf);
  const totalBytes = displayedResult?.bxsfFiles.reduce((sum, file) => sum + (Number(file.sizeBytes) || 0), 0) ?? 0;
  useViewportScrollLock(step === "run" && runIsActive);

  useEffect(() => {
    setHpcResources(cloneResources(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.updated_at]);

  useEffect(() => {
    if (selectedSpinMode === "spin_polarized" && spinChannel === "none") {
      setSpinChannel("up");
    }
    if (selectedSpinMode === "non_spin_polarized" && spinChannel !== "none") {
      setSpinChannel("none");
    }
  }, [selectedSpinMode, spinChannel]);

  useEffect(() => {
    if (!reconnectTaskId) return;
    setActiveTaskId(reconnectTaskId);
    setStep("run");
    void taskContext.reconnectToTask(reconnectTaskId);
  }, [reconnectTaskId, taskContext.reconnectToTask]);

  useEffect(() => {
    if (!activeTaskId || activeTask) return;
    void taskContext.reconnectToTask(activeTaskId);
  }, [activeTaskId, activeTask, taskContext.reconnectToTask]);

  useEffect(() => {
    if (!activeTask) return;
    for (const line of activeTask.output) {
      const jobId = parseRemoteJobId(line);
      if (jobId) setRemoteJobId(jobId);
      const node = parseRemoteNode(line);
      if (node) setRemoteNode(node);
    }
    if (activeTask.status === "completed" && taskResult) {
      setResult(taskResult);
      setStep("results");
    } else if (activeTask.status === "failed" || activeTask.status === "cancelled") {
      setError(activeTask.error ?? "WIEN2k Fermi-surface calculation failed.");
    }
  }, [activeTask, taskResult]);

  function toggleSection(section: SectionKey) {
    setExpandedSections((current) => ({ ...current, [section]: !current[section] }));
  }

  function commandLines(): string[] {
    const spinArg = spinChannel === "none" ? "" : ` -${spinChannel === "up" ? "up" : "dn"}`;
    const spinOrbitArg = spinOrbit ? " -so" : "";
    return [
      "cd \"$SLURM_SUBMIT_DIR\"",
      `# write case.klist for ${kMesh[0]} ${kMesh[1]} ${kMesh[2]}`,
      `x lapw1${spinArg}${spinOrbitArg}`,
      `x lapw2 -fermi${spinArg}${spinOrbitArg}`,
      "xcrysden --wien_fermisurface .",
      "find . -maxdepth 3 -iname '*.bxsf'",
    ];
  }

  function buildSettings(): Wien2kFermiSurfaceSettings {
    if (selectedSpinMode === "spin_polarized" && spinChannel === "none") {
      throw new Error("Spin-polarized WIEN2k Fermi surfaces require an up or down spin channel.");
    }
    return {
      kMesh,
      spinChannel,
      spinOrbit,
      diagnosticLog,
    };
  }

  async function runFermiSurface() {
    if (!selectedScf) return;
    setError(null);
    setResult(null);
    setCalcStartTime(new Date().toISOString());
    setRemoteJobId(null);
    setRemoteNode(null);
    setStep("run");
    try {
      const taskId = await taskContext.startTask(
        "wien2k_fermi_surface",
        {
          projectId,
          cifId,
          sourceScfCalculationId: selectedScf.id,
          settings: buildSettings(),
          resources: hpcResources,
        },
        `WIEN2k Fermi Surface - ${String(selectedScf.parameters?.case_name ?? crystalData.formula_sum ?? "case")}`,
      );
      setActiveTaskId(taskId);
    } catch (reason) {
      setError(String(reason));
    }
  }

  async function launchFermiSurfer() {
    if (!displayedResult) return;
    setIsLaunching(true);
    setError(null);
    try {
      await invoke("launch_fermi_surface_viewer", {
        projectId,
        calculationId: displayedResult.calculationId,
        surfaceFile: displayedResult.primaryFile,
      });
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsLaunching(false);
    }
  }

  function renderHeader() {
    return (
      <div className="wizard-header">
        <button className="back-btn" type="button" disabled={runIsActive} onClick={onBack}>
          ← Exit
        </button>
        <h2>WIEN2k Fermi Surface</h2>
        <div className="step-indicator">
          {STEPS.map((entry, index) => (
            <span key={entry.id} className={index === currentStepIndex ? "active" : index < currentStepIndex ? "completed" : ""}>
              {index + 1}. {entry.label}
            </span>
          ))}
        </div>
      </div>
    );
  }

  function renderSection(key: SectionKey, title: string, body: ReactNode, status?: string) {
    const expanded = expandedSections[key];
    return (
      <section className="config-section collapsible wien2k-scf-section">
        <h4 className="section-header" onClick={() => toggleSection(key)} aria-expanded={expanded}>
          <span className={`collapse-icon ${expanded ? "expanded" : ""}`}>▶</span>
          {title}
          {status && <span className="wien2k-section-status">{status}</span>}
        </h4>
        {expanded && body}
      </section>
    );
  }

  function renderSourceStep() {
    if (sources.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No WIEN2k SCF Calculations Available</h3>
          <p className="warning-text">Fermi-surface calculations require a converged WIEN2k SCF with a retained native case directory.</p>
          <button className="primary-button" type="button" onClick={onBack}>Back to Dashboard</button>
        </div>
      );
    }
    return (
      <div className="wizard-step source-step">
        <div className="source-step-header">
          <h3>Select WIEN2k Source SCF</h3>
        </div>
        <p className="step-description">Choose the converged native case to clone for the dense Fermi-surface mesh.</p>
        <div className="scf-list">
          {sources.map((scf) => {
            const name = getCalculationName(scf);
            const summary = scf.scf_summary;
            return (
              <div
                key={scf.id}
                className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""}`}
                onClick={() => setSelectedScf(scf)}
              >
                <div className="scf-option-header">
                  <input type="radio" checked={selectedScf?.id === scf.id} onChange={() => setSelectedScf(scf)} />
                  {name && <span className="scf-name" title={formatCalculationSourceLabel(scf)}>{name}</span>}
                  <span className="scf-date">{new Date(scf.started_at).toLocaleDateString()}</span>
                </div>
                <div className="scf-details">
                  <span>{String(scf.parameters?.case_name ?? "case")}</span>
                  {summary?.totalEnergy && <span>E = {summary.totalEnergy.value.toFixed(6)} {summary.totalEnergy.unit}</span>}
                  {summary?.fermiEnergyEv != null && <span>E_F = {summary.fermiEnergyEv.toFixed(3)} eV</span>}
                </div>
                <div className="calc-tags">
                  <span className="calc-tag calc-tag-feature calc-tag-hpc">HPC</span>
                  <span className="calc-tag calc-tag-special">WIEN2k</span>
                </div>
              </div>
            );
          })}
        </div>
        <div className="step-actions">
          <button className="secondary-button" type="button" onClick={onBack}>Cancel</button>
          <button className="primary-button" type="button" disabled={!selectedScf} onClick={() => setStep("parameters")}>
            Next: Parameters
          </button>
        </div>
      </div>
    );
  }

  function renderParametersStep() {
    const canRun = Boolean(selectedScf) && (selectedSpinMode !== "spin_polarized" || spinChannel !== "none");
    return (
      <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
        {renderHeader()}
        {error && <div className="error-banner">{error}</div>}
        <div className="wien2k-structure-content">
          <section className="wien2k-structure-controls">
            {renderSection("source", "Source Case", (
              <div className="wien2k-summary">
                <p>{String(selectedScf?.parameters?.case_name ?? "case")} from SCF {selectedScf?.id.slice(0, 8)}.</p>
                {fermiEnergy != null && <p>Fermi energy: {fermiEnergy.toFixed(6)} eV</p>}
                <p>Spin mode: {selectedSpinMode === "spin_polarized" ? "spin-polarized" : "non-spin-polarized"}</p>
              </div>
            ))}
            {renderSection("mesh", "Fermi Surface Mesh", (
              <>
                <div className="wien2k-control-grid wien2k-kmesh-grid">
                  {(["k1", "k2", "k3"] as const).map((label, index) => (
                    <label key={label}>
                      <Wien2kFieldLabel tooltip="Regular unshifted k mesh written to case.klist for the WIEN2k lapw1/lapw2 Fermi-surface pass. XCrySDen's WIEN2k Fermi-surface flow does not support shifted meshes.">
                        {label.toUpperCase()}
                      </Wien2kFieldLabel>
                      <input
                        type="number"
                        min="1"
                        max="500"
                        step="1"
                        value={kMesh[index]}
                        onChange={(event) => {
                          const next = [...kMesh] as [number, number, number];
                          next[index] = clampKMesh(numberField(event.target.value, next[index]));
                          setKMesh(next);
                        }}
                      />
                    </label>
                  ))}
                  <label>
                    <Wien2kFieldLabel tooltip="Spin channel passed to WIEN2k lapw1/lapw2 for spin-polarized source cases.">
                      Spin channel
                    </Wien2kFieldLabel>
                    <select value={spinChannel} onChange={(event) => setSpinChannel(event.target.value as Wien2kBandsSpinChannel)}>
                      <option value="none" disabled={selectedSpinMode === "spin_polarized"}>None</option>
                      <option value="up">Up</option>
                      <option value="down">Down</option>
                    </select>
                  </label>
                </div>
                {selectedSpinMode === "spin_polarized" && spinChannel === "none" && (
                  <p className="wien2k-validation">Spin-polarized WIEN2k sources must run an up or down channel.</p>
                )}
              </>
            ), `${kMesh[0]}×${kMesh[1]}×${kMesh[2]}`)}
            {renderSection("xcrysden", "XCrySDen Conversion", (
              <>
                <p className="warning-text">
                  Fermi surface calculation from WIEN2k SCFs relies on xcrysden. <InfoTooltip text="XCrySDen must be in the remote path." />
                </p>
                <label className="option-checkbox">
                  <input type="checkbox" checked={spinOrbit} onChange={(event) => setSpinOrbit(event.target.checked)} />
                  <span>
                    Enable spin-orbit coupling
                    <Wien2kFieldLabel tooltip="Passes -so to the WIEN2k Fermi-surface lapw1/lapw2 steps. Use only when the source case contains the required spin-orbit files."> </Wien2kFieldLabel>
                  </span>
                </label>
              </>
            ), "BXSF")}
            {renderSection("logging", "Run Logging", (
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="Normal streams Slurm stdout and stderr. Diagnostic tails also prints WIEN2k dayfile, lapw outputs, XCrySDen logs, and error-file tails after native commands.">
                    Verbosity
                  </Wien2kFieldLabel>
                  <select value={diagnosticLog ? "diagnostic" : "normal"} onChange={(event) => setDiagnosticLog(event.target.value === "diagnostic")}>
                    <option value="normal">Normal</option>
                    <option value="diagnostic">Diagnostic tails</option>
                  </select>
                </label>
              </div>
            ), diagnosticLog ? "Diagnostic" : "Normal")}
            {renderSection("hpc", "HPC Run Settings", (
              <HpcRunSettings
                profileId={activeHpcProfile?.id ?? null}
                profileName={activeHpcProfile?.name ?? "Andromeda"}
                taskKind="fermi-surface"
                commandLines={commandLines()}
                resources={hpcResources}
                resourceMode="cpu_only"
                defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
                defaultGpuResources={null}
                onResourcesChange={setHpcResources}
                disabled={runIsActive}
              />
            ), "Ready")}
            <div className="run-actions">
              <button type="button" disabled={runIsActive} onClick={() => setStep("source")}>Back</button>
              <button type="button" className="primary-button" disabled={!canRun || runIsActive} onClick={() => void runFermiSurface()}>
                Run WIEN2k Fermi Surface
              </button>
            </div>
          </section>
          <section className="wien2k-structure-preview">
            <LiveOutputPanel
              title="Run preview"
              output={commandLines().join("\n")}
              placeholder="Run commands will appear here."
              totalLineCount={commandLines().length}
              visibleLineCount={commandLines().length}
            />
          </section>
        </div>
      </div>
    );
  }

  function renderRunStep() {
    const progressStatus: "idle" | "running" | "error" | "complete" =
      runHasFailed ? "error" : displayedResult ? "complete" : runIsActive ? "running" : "idle";
    return (
      <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard wizard-step-run">
        {renderHeader()}
        {error && <div className="error-banner">{error}</div>}
        <div className="wizard-content">
          <div className="wizard-step run-step run-step-focused scf-run-step">
            <div className="run-step-headline">
              <h3>{runIsActive ? "Running WIEN2k Fermi Surface" : "WIEN2k Fermi Surface Output"}</h3>
              <span className={`run-step-status-pill ${runIsActive ? "running" : runHasFailed ? "error" : "idle"}`}>
                {runIsActive ? "Live output" : runHasFailed ? "Run failed" : "Ready"}
              </span>
            </div>
            <div className="run-status-rail scf-run-status">
              <ProgressBar
                status={progressStatus}
                percent={null}
                phase={runIsActive ? "lapw1/lapw2/xcrysden" : displayedResult ? "Complete" : "Ready"}
                detail={remoteJobId ? `Slurm job ${remoteJobId}` : "Waiting for remote allocation"}
                compact
              />
              <div className="run-status-meta">
                <ElapsedTimer startedAt={activeTask?.startedAt ?? calcStartTime} isRunning={runIsActive} />
              </div>
            </div>
            <div className="run-layout run-layout-hpc-telemetry">
              <LiveOutputPanel
                title={runIsActive ? "Running..." : "Output"}
                output={taskOutputText}
                placeholder="Start the WIEN2k Fermi-surface run to see output."
                totalLineCount={taskOutputLineCount}
                visibleLineCount={taskOutputLines.length}
              />
              <RemoteUtilizationPanel
                enabled={runIsActive}
                profileId={activeHpcProfile?.id ?? null}
                remoteJobId={activeTask?.hpc.remote_job_id ?? remoteJobId}
                remoteNode={activeTask?.hpc.remote_node ?? remoteNode}
                resourceType="cpu"
              />
            </div>
            {displayedResult && (
              <div className="wien2k-summary wien2k-scf-results">
                <h3>Fermi Surface Result: Complete</h3>
                <p>Generated {displayedResult.bxsfFiles.length} BXSF file{displayedResult.bxsfFiles.length === 1 ? "" : "s"}.</p>
                <p>Primary file: {displayedResult.primaryFile}</p>
                {displayedResult.fermiEnergy != null && <p>Fermi energy: {displayedResult.fermiEnergy.toFixed(6)} eV</p>}
                {displayedResult.diagnostics.map((diagnostic) => (
                  <p key={diagnostic} className="wien2k-validation">{diagnostic}</p>
                ))}
              </div>
            )}
            <div className="run-actions">
              <button type="button" disabled={runIsActive} onClick={() => setStep("parameters")}>Back</button>
              {!displayedResult && !runIsActive && (
                <button type="button" className="primary-button" disabled={!selectedScf} onClick={() => void runFermiSurface()}>
                  Run Fermi Surface
                </button>
              )}
              {displayedResult && (
                <>
                  <button type="button" className="secondary-button" disabled={runIsActive} onClick={onBack}>
                    Return to Project
                  </button>
                  <button type="button" className="primary-button" disabled={isLaunching} onClick={() => void launchFermiSurfer()}>
                    {isLaunching ? "Launching..." : "View in FermiSurfer"}
                  </button>
                </>
              )}
            </div>
          </div>
        </div>
      </div>
    );
  }

  function renderResultsStep() {
    return (
      <div className="wizard-step results-step">
        <h3>WIEN2k Fermi Surface Complete</h3>
        {displayedResult ? (
          <div className="results-summary">
            <p>Generated {displayedResult.bxsfFiles.length} BXSF file{displayedResult.bxsfFiles.length === 1 ? "" : "s"} ({formatBytes(totalBytes)}).</p>
            <p>Primary file: {displayedResult.primaryFile}</p>
          </div>
        ) : (
          <p>No result is available.</p>
        )}
        <LiveOutputPanel
          title="Output"
          output={taskOutputText}
          placeholder="No output captured."
          totalLineCount={taskOutputLineCount}
          visibleLineCount={taskOutputLines.length}
        />
        <div className="step-actions">
          <button className="secondary-button" type="button" onClick={onBack}>Return to Project</button>
          <button className="primary-button" type="button" disabled={!displayedResult || isLaunching} onClick={() => void launchFermiSurfer()}>
            {isLaunching ? "Launching..." : "View in FermiSurfer"}
          </button>
        </div>
      </div>
    );
  }

  if (step === "parameters") return renderParametersStep();
  if (step === "run") return renderRunStep();

  return (
    <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
      {renderHeader()}
      {error && <div className="error-banner">{error}</div>}
      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
