import { useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import { AppHeaderPortal } from "../AppHeaderPortal";
import { ElapsedTimer } from "../ElapsedTimer";
import { HpcRunSettings } from "../HpcRunSettings";
import { LiveOutputPanel } from "../LiveOutputPanel";
import { ProgressBar } from "../ProgressBar";
import { RemoteUtilizationPanel } from "../RemoteUtilizationPanel";
import { Wien2kFieldLabel } from "./Wien2kFieldLabel";
import { defaultCpuResources } from "../../lib/hpcConfig";
import { formatCalculationSourceLabel, getCalculationName } from "../../lib/calculationNames";
import { getCalculationTagBadges, getCalcTagClass } from "../../lib/calculationTags";
import { useTaskContext } from "../../lib/TaskContext";
import type { CalculationRun } from "../ProjectDashboard";
import type { HpcProfile, SlurmResourceRequest } from "../../lib/types";
import {
  DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
  DEFAULT_WIEN2K_SOC_PREPARE_SETTINGS,
  acceptWien2kSocSymmetry,
  discardWien2kSocSession,
  getWien2kSocRloCandidates,
  prepareWien2kSocSession,
  startWien2kSocSession,
  validateWien2kScfRunSettings,
} from "../../lib/engines/wien2k";
import type {
  Wien2kScfRunSettings,
  Wien2kSocExecutionResult,
  Wien2kSocPrepareSettings,
  Wien2kSocRlo,
  Wien2kSocRloCandidate,
  Wien2kSocSession,
  Wien2kSpinMode,
} from "../../lib/engines/wien2k";
import { useViewportScrollLock } from "../../lib/useViewportScrollLock";

interface Wien2kSocWizardProps {
  projectId: string;
  cifId: string;
  scfCalculations: CalculationRun[];
  activeHpcProfile?: HpcProfile | null;
  reconnectTaskId?: string;
  onBack: () => void;
  onSaved: () => void;
}

type WizardStep = "source" | "prepare" | "symmetry" | "run" | "results";
type SectionKey = "source" | "initialization" | "rlo" | "symmetry" | "cycle" | "hpc";

const STEPS: Array<{ id: WizardStep; label: string }> = [
  { id: "source", label: "Source" },
  { id: "prepare", label: "Initialize SOC" },
  { id: "symmetry", label: "Review" },
  { id: "run", label: "SOC SCF" },
  { id: "results", label: "Results" },
];

function cloneResources(profile: HpcProfile | null | undefined): SlurmResourceRequest {
  const source = profile?.default_cpu_resources ?? defaultCpuResources();
  const threads = Math.max(1, source.ntasks ?? 1) * Math.max(1, source.cpus_per_task ?? 1);
  return {
    ...source,
    resource_type: "cpu",
    nodes: 1,
    ntasks: 1,
    cpus_per_task: threads,
    gpus: 0,
    additional_sbatch: [...(source.additional_sbatch ?? [])],
  };
}

function isSocCalculation(calc: CalculationRun): boolean {
  return Boolean(
    calc.parameters?.spin_orbit
    || calc.parameters?.spinOrbit
    || calc.parameters?.soc
    || calc.tags?.includes("soc"),
  );
}

function isEligibleSource(calc: CalculationRun): boolean {
  return calc.engine_id === "wien2k"
    && calc.calc_type === "scf"
    && calc.scf_summary?.convergence === "converged"
    && typeof calc.parameters?.remote_case_dir === "string"
    && !isSocCalculation(calc);
}

function sourceSpinMode(calc: CalculationRun | null): Wien2kSpinMode {
  const value =
    calc?.parameters?.run?.spinMode
    ?? calc?.parameters?.run?.spin_mode
    ?? calc?.parameters?.initialization?.spinMode
    ?? calc?.parameters?.initialization?.spin_mode;
  return value === "spin_polarized" ? "spin_polarized" : "non_spin_polarized";
}

function inheritedRunSettings(calc: CalculationRun | null): Wien2kScfRunSettings {
  const source = calc?.parameters?.run ?? {};
  const spinMode = sourceSpinMode(calc);
  return {
    ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
    ...source,
    spinMode,
    parallelMode: "openmp",
    forceConvergenceMryBohr: null,
    forceMinimization: false,
    dftU: {
      ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS.dftU,
      ...(source.dftU ?? source.dft_u ?? {}),
    },
    mixer: {
      ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS.mixer,
      ...(source.mixer ?? {}),
    },
  };
}

function numberField(value: string, fallback: number): number {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : fallback;
}

function rloCandidateKey(candidate: Wien2kSocRloCandidate): string {
  return [
    candidate.atomIndex,
    candidate.energyRy,
    candidate.de,
    candidate.switch,
    candidate.sourceFile,
  ].join("|");
}

function matchingRloCandidate(
  rlo: Wien2kSocRlo,
  candidates: Wien2kSocRloCandidate[],
): Wien2kSocRloCandidate | undefined {
  return candidates.find((candidate) => (
    candidate.atomIndex === rlo.atomIndex
    && candidate.energyRy === rlo.energyRy
    && candidate.de === rlo.de
    && candidate.switch === rlo.switch
  ));
}

function parseAtomIndices(value: string): number[] {
  return Array.from(new Set(
    value
      .split(/[\s,]+/)
      .map((entry) => Number.parseInt(entry, 10))
      .filter((entry) => Number.isInteger(entry) && entry > 0),
  )).sort((left, right) => left - right);
}

function rloError(settings: Wien2kSocPrepareSettings): string | null {
  const disabled = new Set(settings.disabledAtomIndices);
  const rloAtoms = new Set<number>();
  for (const rlo of settings.rloAtoms) {
    if (!Number.isInteger(rlo.atomIndex) || rlo.atomIndex <= 0) {
      return "Each RLO must target a positive WIEN2k inequivalent atom index.";
    }
    if (rloAtoms.has(rlo.atomIndex)) {
      return `Atom ${rlo.atomIndex} has more than one RLO entry.`;
    }
    if (disabled.has(rlo.atomIndex)) {
      return `Atom ${rlo.atomIndex} cannot have an RLO while SOC is disabled for that atom.`;
    }
    if (!Number.isFinite(rlo.energyRy) || !Number.isFinite(rlo.de) || rlo.de < 0) {
      return "RLO energy parameters must be finite and DE must be non-negative.";
    }
    rloAtoms.add(rlo.atomIndex);
  }
  return null;
}

function prepareError(settings: Wien2kSocPrepareSettings): string | null {
  if (settings.magnetizationDirection.every((entry) => entry === 0)) {
    return "Magnetization direction cannot be the zero vector.";
  }
  if (!Number.isFinite(settings.lapw1EmaxRy) || settings.lapw1EmaxRy <= 0) {
    return "LAPW1 EMAX must be positive.";
  }
  if (settings.outputEnergyMinRy >= settings.outputEnergyMaxRy) {
    return "The case.inso energy window must be ordered.";
  }
  return rloError(settings);
}

function convergenceLabel(result: Wien2kSocExecutionResult | null): string {
  if (!result) return "SOC SCF output";
  return result.summary.convergence === "converged" ? "Converged" : result.summary.convergence.replace("_", " ");
}

export function Wien2kSocWizard({
  projectId,
  cifId,
  scfCalculations,
  activeHpcProfile = null,
  reconnectTaskId,
  onBack,
  onSaved,
}: Wien2kSocWizardProps) {
  const taskContext = useTaskContext();
  const sources = useMemo(() => scfCalculations.filter(isEligibleSource), [scfCalculations]);
  const [selectedId, setSelectedId] = useState(() => sources[0]?.id ?? "");
  const selectedSource = sources.find((calc) => calc.id === selectedId) ?? sources[0] ?? null;
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [session, setSession] = useState<Wien2kSocSession | null>(null);
  const [prepareSettings, setPrepareSettings] = useState<Wien2kSocPrepareSettings>({
    ...DEFAULT_WIEN2K_SOC_PREPARE_SETTINGS,
    magnetizationDirection: [...DEFAULT_WIEN2K_SOC_PREPARE_SETTINGS.magnetizationDirection],
    disabledAtomIndices: [],
    rloAtoms: [],
  });
  const [disabledAtomsText, setDisabledAtomsText] = useState("");
  const [rloCandidates, setRloCandidates] = useState<Wien2kSocRloCandidate[]>([]);
  const [rloCandidateStatus, setRloCandidateStatus] = useState<"idle" | "loading" | "loaded" | "failed">("idle");
  const [runSettings, setRunSettings] = useState<Wien2kScfRunSettings>(() => inheritedRunSettings(sources[0] ?? null));
  const [diagnosticLog, setDiagnosticLog] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(() => cloneResources(activeHpcProfile));
  const [expandedSections, setExpandedSections] = useState<Record<SectionKey, boolean>>({
    source: true,
    initialization: true,
    rlo: false,
    symmetry: false,
    cycle: false,
    hpc: false,
  });
  const [outputLines, setOutputLines] = useState<string[]>([]);
  const [result, setResult] = useState<Wien2kSocExecutionResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [isPreparing, setIsPreparing] = useState(false);
  const [isAccepting, setIsAccepting] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const [calcStartTime, setCalcStartTime] = useState<string | null>(null);
  const outputUnlistenRef = useRef<UnlistenFn | null>(null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const taskResult = activeTask?.result as Wien2kSocExecutionResult | null | undefined;
  const displayedResult = result ?? taskResult ?? null;
  const runIsActive = activeTask?.status === "running" || isRunning;
  const runHasFailed = activeTask?.status === "failed" || activeTask?.status === "cancelled" || Boolean(error);
  const selectedSpinMode = sourceSpinMode(selectedSource);
  const effectiveRunSettings = useMemo(() => ({
    ...runSettings,
    spinMode: selectedSpinMode,
    forceConvergenceMryBohr: null,
    forceMinimization: false,
  }), [runSettings, selectedSpinMode]);
  const runError = validateWien2kScfRunSettings(effectiveRunSettings);
  const currentStepIndex = STEPS.findIndex((entry) => entry.id === step);
  const taskOutput = activeTask?.outputText ?? outputLines.join("\n");
  const taskOutputLines = activeTask?.output ?? outputLines;
  const taskOutputLineCount = activeTask?.outputLineCount ?? outputLines.length;
  const remoteJobId = activeTask?.hpc.remote_job_id ?? null;
  const remoteNode = activeTask?.hpc.remote_node ?? null;
  const symmetryOutput = session
    ? Object.entries(session.artifacts ?? {}).find(([name]) => name.endsWith(".outsymso"))?.[1] ?? ""
    : "";
  const preparationLocked = Boolean(session) || isPreparing;
  const symmetryReviewRequired = session?.phase === "symmetry_ready";
  const prepared = session?.phase === "prepared" || session?.phase === "soc_complete";
  const sourceSites = useMemo(() => {
    const raw = selectedSource?.parameters?.source_structure_sites ?? selectedSource?.parameters?.sites;
    if (!Array.isArray(raw)) return [];
    return raw
      .map((site, index) => ({
        siteIndex: Number(site?.siteIndex ?? site?.site_index ?? index + 1),
        symbol: String(site?.symbol ?? `Site ${index + 1}`),
      }))
      .filter((site) => Number.isInteger(site.siteIndex) && site.siteIndex > 0);
  }, [selectedSource]);
  useViewportScrollLock(step === "run" && runIsActive);

  useEffect(() => () => outputUnlistenRef.current?.(), []);

  useEffect(() => {
    setHpcResources(cloneResources(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.updated_at]);

  useEffect(() => {
    setRunSettings(inheritedRunSettings(selectedSource));
  }, [selectedSource?.id]);

  useEffect(() => {
    let cancelled = false;
    setRloCandidates([]);
    if (!selectedSource) {
      setRloCandidateStatus("idle");
      return () => {
        cancelled = true;
      };
    }
    setRloCandidateStatus("loading");
    void getWien2kSocRloCandidates(projectId, selectedSource.id)
      .then((candidates) => {
        if (cancelled) return;
        setRloCandidates(candidates);
        setRloCandidateStatus("loaded");
      })
      .catch(() => {
        if (cancelled) return;
        setRloCandidateStatus("failed");
      });
    return () => {
      cancelled = true;
    };
  }, [projectId, selectedSource?.id]);

  useEffect(() => {
    if (step === "symmetry" && symmetryReviewRequired) {
      setExpandedSections((current) => ({
        ...current,
        source: false,
        initialization: false,
        symmetry: true,
        cycle: false,
        hpc: false,
      }));
    }
  }, [step, symmetryReviewRequired]);

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
    if (activeTask.status === "failed" || activeTask.status === "cancelled") {
      setError(activeTask.error ?? "WIEN2k SOC calculation failed.");
      setIsRunning(false);
      return;
    }
    if (activeTask.status === "completed" && taskResult) {
      setResult(taskResult);
      setIsRunning(false);
      setStep("results");
      setSession((current) => current ? {
        ...current,
        phase: taskResult.phase,
        latestCalculationId: taskResult.calculationId,
      } : current);
    }
  }, [activeTask, taskResult]);

  async function attachOutputListener(sessionId: string) {
    outputUnlistenRef.current?.();
    outputUnlistenRef.current = await listen<string>(`wien2k-soc-output:${sessionId}`, (event) => {
      setOutputLines((current) => [...current, event.payload].slice(-1500));
    });
  }

  async function prepareSoc() {
    if (!selectedSource || prepareError(prepareSettings)) return;
    setIsPreparing(true);
    setError(null);
    try {
      let activeSession = session;
      if (!activeSession) {
        const staged = await startWien2kSocSession(projectId, cifId, selectedSource.id);
        activeSession = { ...staged, artifacts: staged.artifacts ?? {} };
        setSession(activeSession);
        setOutputLines([`Copied converged ${activeSession.caseName} case for SOC preparation.`]);
        await attachOutputListener(activeSession.sessionId);
      }
      const prepared = await prepareWien2kSocSession(activeSession.sessionId, prepareSettings);
      setSession((current) => current ? {
        ...current,
        phase: prepared.phase,
        latestPrepare: prepareSettings,
        artifacts: { ...(current.artifacts ?? {}), ...(prepared.artifacts ?? {}) },
      } : current);
      setOutputLines((current) => [...current, prepared.nativeOutput].slice(-1500));
      if (prepared.symmetryReviewRequired) {
        setStep("symmetry");
      } else {
        setExpandedSections((current) => ({
          ...current,
          source: false,
          initialization: false,
          rlo: false,
          symmetry: false,
          cycle: true,
          hpc: true,
        }));
        setStep("symmetry");
      }
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsPreparing(false);
    }
  }

  async function acceptSymmetry() {
    if (!session) return;
    setIsAccepting(true);
    setError(null);
    try {
      const accepted = await acceptWien2kSocSymmetry(session.sessionId);
      setSession((current) => current ? {
        ...current,
        phase: accepted.phase,
        artifacts: { ...(current.artifacts ?? {}), ...(accepted.artifacts ?? {}) },
      } : current);
      setOutputLines((current) => [...current, accepted.nativeOutput].slice(-1500));
      setExpandedSections((current) => ({
        ...current,
        symmetry: false,
        cycle: true,
        hpc: true,
      }));
      setStep("symmetry");
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsAccepting(false);
    }
  }

  function addRlo() {
    const used = new Set(prepareSettings.rloAtoms.map((entry) => entry.atomIndex));
    const candidate = rloCandidates.find((entry) => !used.has(entry.atomIndex))
      ?? null;
    const atomIndex = candidate?.atomIndex
      ?? sourceSites.find((site) => !used.has(site.siteIndex))?.siteIndex
      ?? Array.from({ length: used.size + 1 }, (_, index) => index + 1).find((index) => !used.has(index))
      ?? used.size + 1;
    setPrepareSettings((current) => ({
      ...current,
      rloAtoms: [
        ...current.rloAtoms,
        {
          atomIndex,
          energyRy: candidate?.energyRy ?? -3,
          de: candidate?.de ?? 0.005,
          switch: candidate?.switch ?? "STOP",
        },
      ],
    }));
  }

  function selectRloSite(index: number, atomIndex: number) {
    const candidate = rloCandidates.find((entry) => entry.atomIndex === atomIndex);
    updateRlo(index, candidate
      ? {
          atomIndex,
          energyRy: candidate.energyRy,
          de: candidate.de,
          switch: candidate.switch,
        }
      : { atomIndex });
  }

  function selectRloCandidate(index: number, key: string) {
    const candidate = rloCandidates.find((entry) => rloCandidateKey(entry) === key);
    if (!candidate) return;
    updateRlo(index, {
      atomIndex: candidate.atomIndex,
      energyRy: candidate.energyRy,
      de: candidate.de,
      switch: candidate.switch,
    });
  }

  function updateRlo(index: number, changes: Partial<Wien2kSocRlo>) {
    setPrepareSettings((current) => ({
      ...current,
      rloAtoms: current.rloAtoms.map((entry, entryIndex) => (
        entryIndex === index ? { ...entry, ...changes } : entry
      )),
    }));
  }

  function removeRlo(index: number) {
    setPrepareSettings((current) => ({
      ...current,
      rloAtoms: current.rloAtoms.filter((_, entryIndex) => entryIndex !== index),
    }));
  }

  async function submitSoc(continuation: boolean) {
    if (!session || runError) return;
    setError(null);
    setIsRunning(true);
    setCalcStartTime(new Date().toISOString());
    setStep("run");
    try {
      const taskId = await taskContext.startTask(
        "wien2k_soc",
        {
          sessionId: session.sessionId,
          settings: { scf: effectiveRunSettings, diagnosticLog },
          continuation,
          parentCalculationId: continuation ? displayedResult?.calculationId ?? session.latestCalculationId ?? null : selectedSource?.id ?? null,
          resources: hpcResources,
        },
        `WIEN2k SOC - ${session.caseName}`,
      );
      setActiveTaskId(taskId);
    } catch (reason) {
      setError(String(reason));
      setIsRunning(false);
    }
  }

  async function leaveWizard(saved: boolean) {
    setError(null);
    try {
      if (session && !runIsActive) {
        await discardWien2kSocSession(session.sessionId);
      }
      outputUnlistenRef.current?.();
      outputUnlistenRef.current = null;
      if (saved) onSaved();
      else onBack();
    } catch (reason) {
      setError(String(reason));
    }
  }

  const commandLines = [
    "cd \"$SLURM_SUBMIT_DIR\"",
    `${selectedSpinMode === "spin_polarized" ? "runsp_lapw" : "run_lapw"} -so -ec ${effectiveRunSettings.energyConvergenceRy} -cc ${effectiveRunSettings.chargeConvergence} -i ${effectiveRunSettings.maxIterations}${effectiveRunSettings.dftU.enabled ? " -orb" : ""}`,
  ];

  function toggleSection(key: SectionKey) {
    setExpandedSections((current) => ({ ...current, [key]: !current[key] }));
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
          onClick={() => {
            if (!options.locked) toggleSection(key);
          }}
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

  function renderHeader() {
    return (
      <AppHeaderPortal className="wizard-header">
        <button className="back-btn" type="button" disabled={isPreparing || isAccepting} onClick={() => void leaveWizard(false)}>← Exit</button>
        <h2>WIEN2k SOC</h2>
        <div className="step-indicator">
          {STEPS.map((entry, index) => <span key={entry.id} className={index === currentStepIndex ? "active" : index < currentStepIndex ? "completed" : ""}>{index + 1}. {entry.label}</span>)}
        </div>
      </AppHeaderPortal>
    );
  }

  function renderSourceStep() {
    if (sources.length === 0) {
      return (
        <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
          {renderHeader()}
          {error && <div className="error-banner">{error}</div>}
          <div className="wizard-content">
            <div className="wizard-step source-step">
              <h3>No WIEN2k SCF Calculations Available</h3>
              <p className="warning-text">SOC requires a converged scalar-relativistic WIEN2k SCF with a retained native case directory.</p>
              <button className="primary-button" type="button" onClick={onBack}>Back to Dashboard</button>
            </div>
          </div>
        </div>
      );
    }
    return (
      <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
        {renderHeader()}
        {error && <div className="error-banner">{error}</div>}
        <div className="wizard-content">
          <div className="wizard-step source-step">
            <div className="source-step-header">
              <h3>Select WIEN2k Source SCF</h3>
            </div>
            <p className="step-description">Choose the converged scalar-relativistic native case that will be cloned and prepared for SOC.</p>
            <div className="scf-list">
              {sources.map((scf) => {
                const name = getCalculationName(scf);
                const summary = scf.scf_summary;
                return (
                  <div
                    key={scf.id}
                    className={`scf-option ${selectedSource?.id === scf.id ? "selected" : ""}`}
                    onClick={() => setSelectedId(scf.id)}
                  >
                    <div className="scf-option-header">
                      <input type="radio" checked={selectedSource?.id === scf.id} onChange={() => setSelectedId(scf.id)} />
                      {name && <span className="scf-name" title={formatCalculationSourceLabel(scf)}>{name}</span>}
                      <span className="scf-date">{new Date(scf.started_at).toLocaleDateString()}</span>
                    </div>
                    <div className="scf-details">
                      <span>{String(scf.parameters?.case_name ?? "case")}</span>
                      {summary?.totalEnergy && <span>E = {summary.totalEnergy.value.toFixed(6)} {summary.totalEnergy.unit}</span>}
                      {summary?.fermiEnergyEv != null && <span>E_F = {summary.fermiEnergyEv.toFixed(3)} eV</span>}
                    </div>
                    <div className="calc-tags">
                      {getCalculationTagBadges(
                        { ...scf, tags: [...(scf.tags ?? []), "scalar-relativistic"] },
                        { legacyFallback: true, calcId: scf.id },
                      ).map((tag, i) => (
                        <span key={`${tag.label}-${i}`} className={getCalcTagClass(tag)}>
                          {tag.label}
                        </span>
                      ))}
                    </div>
                  </div>
                );
              })}
            </div>
            <div className="step-actions">
              <button className="secondary-button" type="button" onClick={onBack}>Cancel</button>
              <button className="primary-button" type="button" disabled={!selectedSource} onClick={() => setStep("prepare")}>Next: Configure SOC</button>
            </div>
          </div>
        </div>
      </div>
    );
  }

  if (step === "source") {
    return renderSourceStep();
  }

  if (step === "run" || step === "results") {
    const progressStatus: "idle" | "running" | "error" | "complete" =
      runHasFailed ? "error" : displayedResult ? "complete" : runIsActive ? "running" : "idle";
    return (
      <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard wizard-step-run">
        {renderHeader()}
        {error && <div className="error-banner">{error}</div>}
        <div className="wizard-content">
          <div className="wizard-step run-step run-step-focused scf-run-step">
            <div className="run-step-headline">
              <h3>{runIsActive ? "Running self-consistent SOC" : "WIEN2k SOC Output"}</h3>
              <span className={`run-step-status-pill ${runIsActive ? "running" : runHasFailed ? "error" : "idle"}`}>
                {runIsActive ? "Live output" : runHasFailed ? "Run failed" : convergenceLabel(displayedResult)}
              </span>
            </div>
            <div className="run-status-rail scf-run-status">
              <ProgressBar status={progressStatus} percent={null} phase={runIsActive ? "SOC SCF cycle" : convergenceLabel(displayedResult)} detail={remoteJobId ? `Slurm job ${remoteJobId}` : "WIEN2k run with -so"} compact />
              <div className="run-status-meta"><ElapsedTimer startedAt={activeTask?.startedAt ?? calcStartTime} isRunning={runIsActive} /></div>
            </div>
            <div className="run-layout run-layout-hpc-telemetry">
              <LiveOutputPanel title={runIsActive ? "Running..." : "Output"} output={taskOutput} placeholder="Starting WIEN2k SOC..." totalLineCount={taskOutputLineCount} visibleLineCount={taskOutputLines.length} />
              <RemoteUtilizationPanel enabled={runIsActive} profileId={activeTask?.hpc.hpc_profile_id ?? activeHpcProfile?.id ?? null} remoteJobId={remoteJobId} remoteNode={remoteNode} resourceType={activeTask?.hpc.hpc_resource_type ?? hpcResources.resource_type} />
            </div>
            {displayedResult && (
              <div className="wien2k-summary wien2k-scf-results">
                <h3>SOC Result: {convergenceLabel(displayedResult)}</h3>
                {displayedResult.summary.totalEnergy && <p>Total energy: {displayedResult.summary.totalEnergy.value.toFixed(8)} {displayedResult.summary.totalEnergy.unit}</p>}
                {displayedResult.summary.fermiEnergyEv != null && <p>Fermi energy: {displayedResult.summary.fermiEnergyEv.toFixed(6)} eV</p>}
                {displayedResult.diagnostics.map((diagnostic) => <p key={diagnostic} className="wien2k-validation">{diagnostic}</p>)}
              </div>
            )}
            <div className="run-actions">
              {!runIsActive && !displayedResult && <button type="button" onClick={() => setStep("symmetry")}>Back</button>}
              {displayedResult?.summary.convergence === "not_converged" && <button className="primary-button" type="button" disabled={Boolean(runError) || runIsActive} onClick={() => void submitSoc(true)}>Continue SOC SCF</button>}
              {displayedResult && <button className="secondary-button" type="button" disabled={runIsActive} onClick={() => void leaveWizard(true)}>Return to Project</button>}
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
      {renderHeader()}
      {error && <div className="error-banner">{error}</div>}
      <div className="wien2k-structure-content">
        <section className="wien2k-structure-controls">
          {renderSection("source", "Source Case", (
            <div className="wien2k-summary">
              <p>{String(selectedSource?.parameters?.case_name ?? "case")} from SCF {selectedSource?.id.slice(0, 8)}.</p>
              {selectedSource?.scf_summary?.totalEnergy && <p>Total energy: {selectedSource.scf_summary.totalEnergy.value.toFixed(8)} {selectedSource.scf_summary.totalEnergy.unit}</p>}
              {selectedSource?.scf_summary?.fermiEnergyEv != null && <p>Fermi energy: {selectedSource.scf_summary.fermiEnergyEv.toFixed(6)} eV</p>}
              <p>Spin mode: {selectedSpinMode.replace(/_/g, " ")}. Existing DFT+U settings are inherited.</p>
            </div>
          ), { status: session || isPreparing ? "Locked" : undefined })}

          {renderSection("initialization", "SOC Initialization", (
            <>
              <div className="wien2k-control-grid">
                {[0, 1, 2].map((index) => (
                  <label key={index}>
                    <Wien2kFieldLabel tooltip="Integer h, k, l component of the SOC magnetization direction written to case.inso.">Direction {["h", "k", "l"][index]}</Wien2kFieldLabel>
                    <input type="number" step="1" disabled={preparationLocked} value={prepareSettings.magnetizationDirection[index]} onChange={(event) => setPrepareSettings((current) => {
                      const direction: [number, number, number] = [...current.magnetizationDirection];
                      direction[index] = numberField(event.target.value, direction[index]);
                      return { ...current, magnetizationDirection: direction };
                    })} />
                  </label>
                ))}
                <label>
                  <Wien2kFieldLabel tooltip="Raises the LAPW1 upper energy limit before lapwso. WIEN2k commonly uses 5 Ry.">LAPW1 EMAX (Ry)</Wien2kFieldLabel>
                  <input type="number" min="0" step="0.5" disabled={preparationLocked} value={prepareSettings.lapw1EmaxRy} onChange={(event) => setPrepareSettings((current) => ({ ...current, lapw1EmaxRy: numberField(event.target.value, current.lapw1EmaxRy) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Optional positive atom indices that should not receive spin-orbit coupling.">Atoms without SOC</Wien2kFieldLabel>
                  <input type="text" disabled={preparationLocked} placeholder="e.g. 2, 4" value={disabledAtomsText} onChange={(event) => {
                    setDisabledAtomsText(event.target.value);
                    setPrepareSettings((current) => ({ ...current, disabledAtomIndices: parseAtomIndices(event.target.value) }));
                  }} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Lower energy bound written to case.inso.">Output Emin (Ry)</Wien2kFieldLabel>
                  <input type="number" step="0.1" disabled={preparationLocked} value={prepareSettings.outputEnergyMinRy} onChange={(event) => setPrepareSettings((current) => ({ ...current, outputEnergyMinRy: numberField(event.target.value, current.outputEnergyMinRy) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Upper energy bound written to case.inso.">Output Emax (Ry)</Wien2kFieldLabel>
                  <input type="number" step="0.1" disabled={preparationLocked} value={prepareSettings.outputEnergyMaxRy} onChange={(event) => setPrepareSettings((current) => ({ ...current, outputEnergyMaxRy: numberField(event.target.value, current.outputEnergyMaxRy) }))} />
                </label>
              </div>
              {prepareError(prepareSettings) && <p className="wien2k-validation">{prepareError(prepareSettings)}</p>}
            </>
          ), { status: preparationLocked ? "Locked" : "Ready" })}

          {renderSection("rlo", "Relativistic Local Orbitals (RLO)", (
            <>
              <p className="wien2k-hint">
                Improve the SOC basis for heavy atoms with p semicore states.{" "}
                <Wien2kFieldLabel tooltip="RLOs add relativistic p local orbitals to reduce linearization error for heavy-element semicore states. QCortado lists every l=1 entry from the selected source case.in1c or case.in1; choose the appropriate p semicore entry.">Details</Wien2kFieldLabel>
              </p>
              {rloCandidateStatus === "loading" && <p className="wien2k-hint">Loading l=1 candidates from the source case...</p>}
              {rloCandidateStatus === "loaded" && rloCandidates.length === 0 && <p className="wien2k-hint">No l=1 candidates were found in the source case. RLO values can still be entered manually.</p>}
              {rloCandidateStatus === "failed" && <p className="wien2k-validation">Could not read l=1 candidates from the source case. RLO values can still be entered manually.</p>}
              {prepareSettings.rloAtoms.length > 0 ? (
                <table className="wien2k-site-table wien2k-rlo-table">
                  <thead>
                    <tr>
                      <th>WIEN2k site</th>
                      <th><Wien2kFieldLabel tooltip="Every l=1 entry found for this site in the selected source case.in1c or case.in1. Selecting one fills all RLO values together.">Source p candidate</Wien2kFieldLabel></th>
                      <th><Wien2kFieldLabel tooltip="The l=1 linearization energy in Ry used to construct this relativistic p local orbital. It should match the selected p energy parameter in the site's case.in1 or case.in1c block.">p energy (Ry)</Wien2kFieldLabel></th>
                      <th><Wien2kFieldLabel tooltip="Energy derivative step used when constructing the local-orbital radial basis. Use the DE value associated with the selected p energy parameter in case.in1 or case.in1c.">DE</Wien2kFieldLabel></th>
                      <th><Wien2kFieldLabel tooltip="Controls radial-function continuation for this energy parameter. STOP terminates after the local-orbital solution; CONT continues processing subsequent energy parameters. Match the switch from case.in1 or case.in1c.">Switch</Wien2kFieldLabel></th>
                      <th />
                    </tr>
                  </thead>
                  <tbody>
                    {prepareSettings.rloAtoms.map((rlo, index) => (
                      <tr key={`${index}-${rlo.atomIndex}`}>
                        <td>
                          {sourceSites.length > 0 ? (
                            <select disabled={preparationLocked} value={rlo.atomIndex} onChange={(event) => selectRloSite(index, Number(event.target.value))}>
                              {!sourceSites.some((site) => site.siteIndex === rlo.atomIndex) && <option value={rlo.atomIndex}>Site {rlo.atomIndex}</option>}
                              {sourceSites.map((site) => <option key={site.siteIndex} value={site.siteIndex}>{site.siteIndex}. {site.symbol}</option>)}
                            </select>
                          ) : (
                            <input type="number" min="1" step="1" disabled={preparationLocked} value={rlo.atomIndex} onChange={(event) => selectRloSite(index, Math.max(1, Math.round(numberField(event.target.value, rlo.atomIndex))))} />
                          )}
                        </td>
                        <td>
                          <select
                            disabled={preparationLocked}
                            value={matchingRloCandidate(rlo, rloCandidates) ? rloCandidateKey(matchingRloCandidate(rlo, rloCandidates)!) : ""}
                            onChange={(event) => selectRloCandidate(index, event.target.value)}
                          >
                            <option value="">Manual values</option>
                            {rloCandidates
                              .filter((candidate) => candidate.atomIndex === rlo.atomIndex)
                              .map((candidate) => (
                                <option key={rloCandidateKey(candidate)} value={rloCandidateKey(candidate)}>
                                  E={candidate.energyRy}, DE={candidate.de}, {candidate.switch} ({candidate.sourceFile})
                                </option>
                              ))}
                          </select>
                        </td>
                        <td><input type="number" step="0.01" disabled={preparationLocked} value={rlo.energyRy} onChange={(event) => updateRlo(index, { energyRy: numberField(event.target.value, rlo.energyRy) })} /></td>
                        <td><input type="number" min="0" step="0.001" disabled={preparationLocked} value={rlo.de} onChange={(event) => updateRlo(index, { de: numberField(event.target.value, rlo.de) })} /></td>
                        <td>
                          <select disabled={preparationLocked} value={rlo.switch} onChange={(event) => updateRlo(index, { switch: event.target.value as Wien2kSocRlo["switch"] })}>
                            <option value="STOP">STOP</option>
                            <option value="CONT">CONT</option>
                          </select>
                        </td>
                        <td><button className="secondary-button" type="button" disabled={preparationLocked} onClick={() => removeRlo(index)}>Remove</button></td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              ) : (
                <p className="wien2k-hint">No RLOs configured.</p>
              )}
              <button className="secondary-button" type="button" disabled={preparationLocked || (sourceSites.length > 0 && prepareSettings.rloAtoms.length >= sourceSites.length)} onClick={addRlo}>
                Add RLO
              </button>
              {rloError(prepareSettings) && <p className="wien2k-validation">{rloError(prepareSettings)}</p>}
            </>
          ), {
            status: preparationLocked
              ? "Locked"
              : prepareSettings.rloAtoms.length > 0
                ? `${prepareSettings.rloAtoms.length} configured`
                : "None",
          })}

          {renderSection("symmetry", "Review symmetso Proposal", (
            <>
              <p className="wien2k-hint">Spin-polarized SOC may lower symmetry and split equivalent atoms. Accept the generated `_so` files only after reviewing this output.</p>
              <pre className="output-text wien2k-inline-output-text">{symmetryOutput || "No case.outsymso text was retrieved."}</pre>
              <button className="primary-button" type="button" disabled={isAccepting} onClick={() => void acceptSymmetry()}>{isAccepting ? "Accepting..." : "Accept Symmetry Proposal"}</button>
            </>
          ), {
            locked: !symmetryReviewRequired,
            status: symmetryReviewRequired ? "Review required" : prepared ? "Accepted / not required" : "Locked - run preparation first",
          })}

          {renderSection("cycle", "Self-consistent SOC Cycle", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="Maximum self-consistent SOC cycles.">Maximum iterations</Wien2kFieldLabel>
                  <input type="number" min="1" step="1" value={runSettings.maxIterations} onChange={(event) => setRunSettings((current) => ({ ...current, maxIterations: numberField(event.target.value, current.maxIterations) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Total-energy convergence threshold.">Energy convergence (Ry)</Wien2kFieldLabel>
                  <input type="number" min="0" step="0.00001" value={runSettings.energyConvergenceRy} onChange={(event) => setRunSettings((current) => ({ ...current, energyConvergenceRy: numberField(event.target.value, current.energyConvergenceRy) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Charge-density convergence threshold.">Charge convergence</Wien2kFieldLabel>
                  <input type="number" min="0" step="0.00001" value={runSettings.chargeConvergence} onChange={(event) => setRunSettings((current) => ({ ...current, chargeConvergence: numberField(event.target.value, current.chargeConvergence) }))} />
                </label>
                <label className="checkbox-label">
                  <input type="checkbox" checked={diagnosticLog} onChange={(event) => setDiagnosticLog(event.target.checked)} />
                  <span>Diagnostic command logging</span>
                </label>
              </div>
              <p className="wien2k-hint">Force convergence and minimization are disabled because WIEN2k forces are incorrect for atoms treated with SOC.</p>
              {runError && <p className="wien2k-validation">{runError}</p>}
            </>
          ), { locked: !prepared, status: prepared ? "Ready" : "Locked - prepare SOC first" })}

          {renderSection("hpc", "HPC Run Settings", (
            <HpcRunSettings profileId={activeHpcProfile?.id ?? null} profileName={activeHpcProfile?.name ?? "Andromeda"} taskKind="wien2k-soc" commandLines={commandLines} resources={hpcResources} resourceMode="cpu_only" defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null} defaultGpuResources={null} resourceModeMessage="WIEN2k SOC uses one CPU Slurm task with OpenMP threads." maxTasks={{ value: 1, reason: "Distributed WIEN2k .machines mode is not enabled by this wizard." }} onResourcesChange={(next) => setHpcResources({ ...next, resource_type: "cpu", nodes: 1, ntasks: 1, gpus: 0 })} disabled={isPreparing || isAccepting} />
          ), { locked: !prepared, status: prepared ? "Ready" : "Locked - prepare SOC first" })}
        </section>
        <section className="wien2k-output-column">
          <LiveOutputPanel title="WIEN2k SOC preparation" output={outputLines.join("\n")} placeholder="SOC preparation output will appear here." totalLineCount={outputLines.length} visibleLineCount={outputLines.length} panelClassName="output-panel wien2k-inline-output" outputClassName="output-text wien2k-inline-output-text" />
        </section>
      </div>
      <div className="step-actions">
        {!session && <button className="secondary-button" type="button" disabled={isPreparing} onClick={() => setStep("source")}>Back</button>}
        {(!session || session.phase === "staged") && <button className="primary-button" type="button" disabled={!selectedSource || Boolean(prepareError(prepareSettings)) || isPreparing} onClick={() => void prepareSoc()}>{isPreparing ? "Preparing..." : "Prepare SOC"}</button>}
        {prepared && <button className="primary-button" type="button" disabled={Boolean(runError) || isRunning} onClick={() => void submitSoc(false)}>Run self-consistent SOC</button>}
      </div>
    </div>
  );
}
