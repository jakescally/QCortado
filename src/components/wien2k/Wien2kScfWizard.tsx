import { useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import { AppHeaderPortal } from "../AppHeaderPortal";
import {
  DEFAULT_WIEN2K_INITIALIZATION_SETTINGS,
  DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
  WIEN2K_EXCHANGE_CORRELATION_OPTIONS,
  discardWien2kScfSession,
  listWien2kStructureSources,
  startWien2kScfContinuationSession,
  startWien2kScfSession,
  validateWien2kInitializationSettings,
  validateWien2kScfRunSettings,
} from "../../lib/engines/wien2k";
import type {
  NormalizedScfSummary,
  Wien2kHubbardTarget,
  Wien2kInitializationResult,
  Wien2kInitializationSettings,
  Wien2kScfExecutionResult,
  Wien2kScfRunSettings,
  Wien2kScfSession,
  Wien2kStructureSourceRecord,
} from "../../lib/engines/wien2k";
import { defaultCpuResources } from "../../lib/hpcConfig";
import type { HpcProfile, SlurmResourceRequest } from "../../lib/types";
import { getHubbardRecommendations, resolveHubbardUDefault } from "../../lib/engines/qe/hubbard";
import { useTaskContext } from "../../lib/TaskContext";
import { useViewportScrollLock } from "../../lib/useViewportScrollLock";
import { ElapsedTimer } from "../ElapsedTimer";
import { HpcRunSettings } from "../HpcRunSettings";
import { InfoTooltip } from "../InfoTooltip";
import { LiveOutputPanel } from "../LiveOutputPanel";
import { ProgressBar } from "../ProgressBar";
import { RemoteUtilizationPanel } from "../RemoteUtilizationPanel";
import { Wien2kFieldLabel } from "./Wien2kFieldLabel";

interface Wien2kScfWizardProps {
  projectId: string;
  cifId: string;
  calculations: Wien2kStructureSourceRecord[];
  activeHpcProfile?: HpcProfile | null;
  reconnectTaskId?: string;
  continuationCalculationId?: string | null;
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

type SectionKey = "source" | "radii" | "initialization" | "magnetism" | "dftu" | "scf" | "advanced" | "hpc";

const ORBITAL_L_BY_MANIFOLD: Record<string, number> = { s: 0, p: 1, d: 2, f: 3 };

function normalizeElementSymbol(symbol: string): string {
  const trimmed = String(symbol || "").trim().replace(/[\d+-]+$/, "");
  if (!trimmed) return "";
  return trimmed.charAt(0).toUpperCase() + trimmed.slice(1).toLowerCase();
}

function orbitalLFromManifold(manifold: string): number {
  return ORBITAL_L_BY_MANIFOLD[manifold.trim().slice(-1)] ?? 2;
}

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

function cloneWien2kKpointResources(profile: HpcProfile | null | undefined): SlurmResourceRequest {
  const source = profile?.default_cpu_resources ?? defaultCpuResources();
  const nodes = Math.max(1, source.nodes ?? 1);
  return {
    ...source,
    resource_type: "cpu",
    nodes,
    ntasks: Math.max(nodes, source.ntasks ?? 1),
    cpus_per_task: Math.max(1, source.cpus_per_task ?? 1),
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

function buildWien2kParallelPreview(
  resources: SlurmResourceRequest,
  parallelMode: Wien2kScfRunSettings["parallelMode"],
): string {
  const cpusPerTask = Math.max(1, resources.cpus_per_task ?? 1);
  if (parallelMode === "kpoint") {
    return `rm -f .machines .processes && export OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK:-${cpusPerTask}}" && printf 'granularity:1\\n' > .machines && srun --nodes="\${SLURM_JOB_NUM_NODES:-1}" --ntasks="\${SLURM_NTASKS:-1}" hostname | awk 'NF { print "1:" $1 }' >> .machines && printf 'extrafine:1\\n' >> .machines`;
  }
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
    ...(settings.iterativeDiagonalization ? ["-it"] : []),
    ...(settings.forceMinimization ? ["-min"] : []),
    ...(settings.dispersionCorrection !== "none" ? [`-${settings.dispersionCorrection}`] : []),
    ...(settings.dftU.enabled ? ["-orb"] : []),
    ...(settings.parallelMode === "kpoint" ? ["-p"] : []),
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
  reconnectTaskId,
  continuationCalculationId = null,
  onBack,
  onSaved,
}: Wien2kScfWizardProps) {
  const taskContext = useTaskContext();
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
  const [isPreparingContinuation, setIsPreparingContinuation] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [expandedSections, setExpandedSections] = useState<Record<SectionKey, boolean>>({
    source: true,
    radii: true,
    initialization: true,
    magnetism: false,
    dftu: false,
    scf: false,
    advanced: false,
    hpc: false,
  });
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(() => cloneWien2kCpuResources(activeHpcProfile));
  const [scfRunStarted, setScfRunStarted] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState<string | null>(null);
  const [remoteJobId, setRemoteJobId] = useState<string | null>(null);
  const [remoteNode, setRemoteNode] = useState<string | null>(null);
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const [activeInitializationTaskId, setActiveInitializationTaskId] = useState<string | null>(null);
  const [continuationParentCalculationId, setContinuationParentCalculationId] = useState<string | null>(
    continuationCalculationId,
  );
  const outputUnlistenRef = useRef<UnlistenFn | null>(null);
  const reconnectTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const activeTask = reconnectTask?.taskType === "wien2k_scf" ? reconnectTask : undefined;
  const activeInitializationTask = activeInitializationTaskId
    ? taskContext.getTask(activeInitializationTaskId)
    : reconnectTask?.taskType === "wien2k_scf_initialize"
      ? reconnectTask
      : undefined;
  const selectedSource = sources.find((source) => source.id === sourceId) ?? sources[0] ?? null;
  const initError = validateWien2kInitializationSettings(initialization);
  const effectiveRunSettings = useMemo(
    () => ({ ...runSettings, spinMode: initialization.spinMode }),
    [runSettings, initialization.spinMode],
  );
  const runError = validateWien2kScfRunSettings(effectiveRunSettings);
  const initializationRunning = activeInitializationTask?.status === "running" || isInitializing;
  const initialized = session?.phase === "initialized" || session?.phase === "scf_complete" || session?.phase === "failed";
  const initializationLocked = initialized || initializationRunning;
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
  const hubbardRecommendations = useMemo(() => {
    const recommendations = getHubbardRecommendations(sites.map((site) => site.symbol ?? ""));
    return new Map(recommendations.map((entry) => [entry.element, entry]));
  }, [sites]);
  const recommendedHubbardTargets = useMemo<Wien2kHubbardTarget[]>(() => {
    const targets: Wien2kHubbardTarget[] = [];
    sites.forEach((site, index) => {
      const element = normalizeElementSymbol(site.symbol ?? "");
      const recommendation = hubbardRecommendations.get(element);
      if (!element || !recommendation) return;
      const defaultU = resolveHubbardUDefault(element, recommendation.manifold, []);
      targets.push({
        siteIndex: index + 1,
        element,
        manifold: recommendation.manifold,
        orbitalL: orbitalLFromManifold(recommendation.manifold),
        uEv: defaultU.value,
        jEv: 0,
        recommended: true,
        reason: recommendation.reason,
      });
    });
    return targets;
  }, [sites, hubbardRecommendations]);
  const hpcCommandLines = useMemo(
    () => [
      "cd \"$SLURM_SUBMIT_DIR\"",
      ...moduleSetupLines(activeHpcProfile),
      buildWien2kParallelPreview(hpcResources, effectiveRunSettings.parallelMode),
      buildRunPreviewCommand(activeHpcProfile, effectiveRunSettings),
    ],
    [activeHpcProfile, effectiveRunSettings, hpcResources],
  );

  const taskResult = activeTask?.result as Wien2kScfExecutionResult | null | undefined;
  const displayedResult = result ?? taskResult ?? null;
  const taskOutputLines = activeTask ? activeTask.output : outputLines;
  const taskOutputText = activeTask ? activeTask.outputText : outputLines.join("\n");
  const taskOutputLineCount = activeTask ? activeTask.outputLineCount : outputLines.length;
  const initializationOutputLines = activeInitializationTask ? activeInitializationTask.output : outputLines;
  const initializationOutputText = activeInitializationTask ? activeInitializationTask.outputText : outputLines.join("\n");
  const initializationOutputLineCount = activeInitializationTask ? activeInitializationTask.outputLineCount : outputLines.length;
  const runIsActive = activeTask?.status === "running" || isRunning;
  const runHasFailed = activeTask?.status === "failed" || activeTask?.status === "cancelled" || Boolean(error);
  const runRemoteJobId = activeTask?.hpc.remote_job_id ?? remoteJobId;
  const runRemoteNode = activeTask?.hpc.remote_node ?? remoteNode;
  const isContinuationSession = Boolean(continuationParentCalculationId);

  useViewportScrollLock(scfRunStarted && runIsActive);

  useEffect(() => () => outputUnlistenRef.current?.(), []);

  useEffect(() => {
    if (!continuationCalculationId) return;
    let cancelled = false;
    setIsPreparingContinuation(true);
    setError(null);
    setResult(null);
    setContinuationParentCalculationId(continuationCalculationId);
    void startWien2kScfContinuationSession(projectId, cifId, continuationCalculationId)
      .then(async (continuationSession) => {
        if (cancelled) return;
        setSession(continuationSession);
        setSourceId(continuationSession.sourceStructureCalculationId);
        if (continuationSession.initialization) {
          setInitialization(continuationSession.initialization);
        }
        if (continuationSession.latestRun) {
          setRunSettings(continuationSession.latestRun);
        }
        setOutputLines([
          `Reopened retained ${continuationSession.caseName} case for continuation.`,
          `Remote case: ${continuationSession.remoteCaseDir}`,
        ]);
        setExpandedSections({
          source: false,
          radii: false,
          initialization: false,
          magnetism: false,
          dftu: Boolean(continuationSession.latestRun?.dftU.enabled),
          scf: true,
          advanced: false,
          hpc: true,
        });
        await attachOutputListener(continuationSession.sessionId);
      })
      .catch((reason) => {
        if (!cancelled) setError(String(reason));
      })
      .finally(() => {
        if (!cancelled) setIsPreparingContinuation(false);
      });
    return () => {
      cancelled = true;
    };
  }, [continuationCalculationId, projectId, cifId]);

  useEffect(() => {
    if (!reconnectTaskId) return;
    setActiveTaskId(reconnectTaskId);
    void taskContext.reconnectToTask(reconnectTaskId);
  }, [reconnectTaskId, taskContext.reconnectToTask]);

  useEffect(() => {
    if (!activeTaskId || reconnectTask) return;
    void taskContext.reconnectToTask(activeTaskId);
  }, [activeTaskId, reconnectTask, taskContext.reconnectToTask]);

  useEffect(() => {
    if (!reconnectTask) return;
    if (reconnectTask.taskType === "wien2k_scf_initialize") {
      setActiveInitializationTaskId(reconnectTask.taskId);
      setScfRunStarted(false);
    } else if (reconnectTask.taskType === "wien2k_scf") {
      setScfRunStarted(true);
    }
  }, [reconnectTask]);

  useEffect(() => {
    if (!activeTask) return;
    if (activeTask.status === "failed" || activeTask.status === "cancelled") {
      setError(activeTask.error ?? (activeTask.status === "cancelled" ? "Cancelled by user" : "WIEN2k SCF failed."));
      setIsRunning(false);
      return;
    }
    if (activeTask.status !== "completed" || !taskResult) return;
    setResult(taskResult);
    setSession((current) => current ? {
      ...current,
      phase: taskResult.phase,
      latestRun: effectiveRunSettings,
      latestCalculationId: taskResult.calculationId,
    } : current);
    setContinuationParentCalculationId(taskResult.calculationId);
    setIsRunning(false);
  }, [activeTask, taskResult, effectiveRunSettings]);

  useEffect(() => {
    if (!activeInitializationTask) return;
    const resume = activeInitializationTask.metadata?.wizardResume;
    if (resume?.sourceId) setSourceId(String(resume.sourceId));
    if (resume?.initialization) setInitialization(resume.initialization as Wien2kInitializationSettings);
    if (resume?.runSettings) setRunSettings(resume.runSettings as Wien2kScfRunSettings);
    if (resume?.session) setSession(resume.session as Wien2kScfSession);
    if (activeInitializationTask.status === "failed" || activeInitializationTask.status === "cancelled") {
      setError(activeInitializationTask.error ?? "WIEN2k initialization failed.");
      setIsInitializing(false);
      return;
    }
    if (activeInitializationTask.status !== "completed") return;
    const initResult = activeInitializationTask.result as (Wien2kInitializationResult & { session?: Wien2kScfSession }) | null;
    if (!initResult) return;
    setSession(initResult.session ?? ((current) => current ? {
      ...current,
      phase: initResult.phase,
      initialization,
      artifacts: { ...current.artifacts, ...initResult.artifacts },
    } : current));
    setOutputLines(activeInitializationTask.output);
    setExpandedSections({
      source: false,
      radii: false,
      initialization: false,
      magnetism: false,
      dftu: runSettings.dftU.enabled,
      scf: true,
      advanced: false,
      hpc: true,
    });
    setIsInitializing(false);
  }, [activeInitializationTask, initialization, runSettings]);

  useEffect(() => {
    setHpcResources(runSettings.parallelMode === "kpoint"
      ? cloneWien2kKpointResources(activeHpcProfile)
      : cloneWien2kCpuResources(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.updated_at, runSettings.parallelMode]);

  useEffect(() => {
    if (sites.length === 0) return;
    setInitialization((current) => {
      if (current.startingMagnetization.length > 0) return current;
      return {
        ...current,
        startingMagnetization: sites.map((site, index) => ({
          siteIndex: index + 1,
          element: normalizeElementSymbol(site.symbol ?? `Site ${index + 1}`),
          configuration: "up",
          momentBohrMagneton: 1,
        })),
      };
    });
  }, [sites]);

  function toggleSection(section: SectionKey) {
    if ((section === "scf" || section === "advanced" || section === "hpc") && !initialized) return;
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
      const taskId = await taskContext.startTask(
        "wien2k_scf_initialize",
        {
          sessionId: activeSession.sessionId,
          settings: initialization,
          workflowMetadata: {
            wizardResume: {
              view: "wien2k-scf-wizard",
              utility: "initialization",
              context: { projectId, cifId, calculations },
              sourceId: selectedSource.id,
              initialization,
              runSettings,
              session: activeSession,
            },
          },
        },
        `WIEN2k Initialization - ${activeSession.caseName}`,
      );
      setActiveInitializationTaskId(taskId);
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
      const taskId = await taskContext.startTask(
        "wien2k_scf",
        {
          sessionId: session.sessionId,
          settings: effectiveRunSettings,
          continuation,
          parentCalculationId: continuation
            ? continuationParentCalculationId ?? result?.calculationId ?? session.latestCalculationId ?? null
            : null,
          resources: hpcResources,
        },
        `WIEN2k SCF - ${caseName}`,
      );
      setActiveTaskId(taskId);
    } catch (reason) {
      setError(String(reason));
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
      setActiveTaskId(null);
      setActiveInitializationTaskId(null);
      setContinuationParentCalculationId(null);
      setExpandedSections({
        source: true,
        radii: true,
        initialization: true,
        magnetism: false,
        dftu: false,
        scf: false,
        advanced: false,
        hpc: false,
      });
    } catch (reason) {
      setError(String(reason));
    }
  }

  async function leaveWizard(destination: "back" | "saved") {
    setError(null);
    try {
      if (session && !runIsActive && !initializationRunning) {
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

  function setSpinMode(spinMode: Wien2kInitializationSettings["spinMode"]) {
    setInitialization((current) => ({ ...current, spinMode }));
    setRunSettings((current) => ({
      ...current,
      spinMode,
      dftU: spinMode === "spin_polarized" ? current.dftU : { ...current.dftU, enabled: false },
    }));
  }

  function setParallelMode(parallelMode: Wien2kScfRunSettings["parallelMode"]) {
    setRunSettings((current) => ({ ...current, parallelMode }));
    setHpcResources((current) => {
      if (parallelMode === "kpoint") {
        return cloneWien2kKpointResources(activeHpcProfile);
      }
      const openMpThreads = Math.max(1, current.ntasks ?? 1) * Math.max(1, current.cpus_per_task ?? 1);
      return {
        ...current,
        nodes: 1,
        ntasks: 1,
        cpus_per_task: openMpThreads,
      };
    });
  }

  function enableRecommendedDftU(enabled: boolean) {
    setRunSettings((current) => ({
      ...current,
      spinMode: enabled ? "spin_polarized" : current.spinMode,
      dftU: {
        ...current.dftU,
        enabled,
        targets: enabled && current.dftU.targets.length === 0 ? recommendedHubbardTargets : current.dftU.targets,
      },
    }));
    if (enabled) {
      setInitialization((current) => ({ ...current, spinMode: "spin_polarized" }));
      setExpandedSections((current) => ({ ...current, magnetism: true, dftu: true }));
    }
  }

  function toggleHubbardTarget(target: Wien2kScfRunSettings["dftU"]["targets"][number], enabled: boolean) {
    setRunSettings((current) => {
      const targets = enabled
        ? [...current.dftU.targets.filter((entry) => entry.siteIndex !== target.siteIndex || entry.manifold !== target.manifold), target]
        : current.dftU.targets.filter((entry) => entry.siteIndex !== target.siteIndex || entry.manifold !== target.manifold);
      return { ...current, dftU: { ...current.dftU, targets } };
    });
  }

  function updateHubbardTarget(
    siteIndex: number,
    changes: Partial<Wien2kScfRunSettings["dftU"]["targets"][number]>,
  ) {
    setRunSettings((current) => ({
      ...current,
      dftU: {
        ...current.dftU,
        targets: current.dftU.targets.map((target) => (
          target.siteIndex === siteIndex ? { ...target, ...changes } : target
        )),
      },
    }));
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
    if (!displayedResult) return null;
    return (
      <div className="wien2k-summary wien2k-scf-results">
        <h3>SCF Result: {convergenceLabel(displayedResult.summary)}</h3>
        <p>
          {displayedResult.summary.totalEnergy ? `E = ${displayedResult.summary.totalEnergy.value.toFixed(8)} ${displayedResult.summary.totalEnergy.unit}; ` : ""}
          {displayedResult.summary.scfSteps != null ? `${displayedResult.summary.scfSteps} iterations` : "Iterations unavailable"}
        </p>
        {displayedResult.summary.fermiEnergyEv != null && <p>Fermi energy: {displayedResult.summary.fermiEnergyEv.toFixed(6)} eV</p>}
        {displayedResult.summary.totalMagnetization != null && <p>Total magnetization: {displayedResult.summary.totalMagnetization.toFixed(6)}</p>}
        {displayedResult.diagnostics.map((diagnostic) => <p key={diagnostic} className="wien2k-validation">{diagnostic}</p>)}
      </div>
    );
  }

  function renderRunStep() {
    const progressStatus: "idle" | "running" | "error" | "complete" =
      runHasFailed ? "error" : displayedResult ? "complete" : runIsActive ? "running" : "idle";
    return (
      <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard wizard-step-run">
        <AppHeaderPortal className="wizard-header">
          <button className="back-btn" type="button" onClick={() => void leaveWizard("back")}>
            ← Exit
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
        </AppHeaderPortal>
        {error && <div className="error-banner">{error}</div>}
        {isPreparingContinuation && (
          <div className="info-banner">Loading saved WIEN2k continuation state...</div>
        )}
        <div className="wizard-content">
          <div className="wizard-step run-step run-step-focused scf-run-step">
            <div className="run-step-headline">
              <h3>{runIsActive ? "Running WIEN2k SCF" : "WIEN2k SCF Output"}</h3>
              <span className={`run-step-status-pill ${runIsActive ? "running" : runHasFailed ? "error" : "idle"}`}>
                {runIsActive ? "Live output" : runHasFailed ? "Run failed" : "Output"}
              </span>
            </div>

            <div className="run-status-rail scf-run-status">
              <ProgressBar
                status={progressStatus}
                percent={null}
                phase={runIsActive ? "SCF cycle" : displayedResult ? convergenceLabel(displayedResult.summary) : "SCF output"}
                detail={runRemoteJobId ? `Slurm job ${runRemoteJobId}` : "Waiting for remote allocation"}
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
                placeholder="Starting WIEN2k SCF..."
                totalLineCount={taskOutputLineCount}
                visibleLineCount={taskOutputLines.length}
              />
              <RemoteUtilizationPanel
                enabled={runIsActive}
                profileId={activeTask?.hpc.hpc_profile_id ?? activeHpcProfile?.id ?? null}
                remoteJobId={runRemoteJobId}
                remoteNode={runRemoteNode}
                resourceType={activeTask?.hpc.hpc_resource_type ?? hpcResources.resource_type}
              />
            </div>
            {renderResultSummary()}
            <div className="run-actions">
              <button type="button" disabled={runIsActive} onClick={() => setScfRunStarted(false)}>
                Back
              </button>
              {displayedResult?.summary.convergence === "not_converged" && (
                <button type="button" className="primary-button" disabled={Boolean(runError) || runIsActive} onClick={() => void submitScf(true)}>
                  {runIsActive ? "Continuing..." : "Continue SCF"}
                </button>
              )}
              {displayedResult && (
                <button type="button" className="secondary-button" disabled={runIsActive} onClick={() => void leaveWizard("saved")}>
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
      <AppHeaderPortal className="wizard-header">
        <button className="back-btn" type="button" disabled={isRunning} onClick={() => void leaveWizard("back")}>
          ← Exit
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
      </AppHeaderPortal>
      {error && <div className="error-banner">{error}</div>}
      {isPreparingContinuation && (
        <div className="info-banner">Loading saved WIEN2k continuation state...</div>
      )}
      <div className="wien2k-structure-content">
        <section className="wien2k-structure-controls">
          {renderSection("source", "Accepted Structure Source", (
            <div className="wien2k-summary">
              <label>
                <Wien2kFieldLabel tooltip="Selects the accepted WIEN2k structure, symmetry, and muffin-tin radii used as the SCF basis. Recommended choice: the latest reviewed structure for the intended material; change geometry or RMT in the Structure workflow.">
                  Saved case.struct
                </Wien2kFieldLabel>
                <select value={sourceId} disabled={Boolean(session) || isInitializing} onChange={(event) => setSourceId(event.target.value)}>
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
          ), { status: session || isInitializing ? "Locked" : undefined })}

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
          ), { status: initializationLocked ? "Locked" : undefined })}

          {renderSection("initialization", "Initialization", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="RMT(min) times KMAX controls the LAPW plane-wave basis size; larger values improve basis accuracy while increasing memory and runtime. Typical starting range: 6 to 9, followed by convergence testing.">
                    RKMAX
                  </Wien2kFieldLabel>
                  <input type="number" min="0.1" step="0.1" disabled={initializationLocked} value={initialization.rkmax} onChange={(event) => setInitialization((current) => ({ ...current, rkmax: numberField(event.target.value, current.rkmax) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Sets the reciprocal-space cutoff for representing density and potential components; raising it resolves sharper density features at higher cost. Typical starting range: 10 to 16.">
                    GMAX
                  </Wien2kFieldLabel>
                  <input type="number" min="0.1" step="0.1" disabled={initializationLocked} value={initialization.gmax} onChange={(event) => setInitialization((current) => ({ ...current, gmax: numberField(event.target.value, current.gmax) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Maximum angular momentum in the muffin-tin expansion; higher values capture more aspherical character but add work. Typical starting range: 8 to 12.">
                    LMAX
                  </Wien2kFieldLabel>
                  <input type="number" min="1" step="1" disabled={initializationLocked} value={initialization.lmax} onChange={(event) => setInitialization((current) => ({ ...current, lmax: numberField(event.target.value, current.lmax) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Chooses the exchange-correlation approximation used to initialize and run the density. WIEN2k initialization provides LDA, PBE, WC, and PBEsol; PBE is a common general-purpose default, and comparisons should keep XC fixed.">
                    XC functional
                  </Wien2kFieldLabel>
                  <select disabled={initializationLocked} value={initialization.exchangeCorrelation} onChange={(event) => setInitialization((current) => ({ ...current, exchangeCorrelation: Number(event.target.value) }))}>
                    {WIEN2K_EXCHANGE_CORRELATION_OPTIONS.map((option) => (
                      <option key={option.value} value={option.value}>{option.label}</option>
                    ))}
                  </select>
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Separates core and valence states during LSTART, affecting which electrons participate in the variational SCF basis. Typical starting range: -8 to -4 Ry; check for core leakage warnings.">
                    LSTART cutoff (Ry)
                  </Wien2kFieldLabel>
                  <input type="number" step="0.1" disabled={initializationLocked} value={initialization.lstartEnergyCutoffRy} onChange={(event) => setInitialization((current) => ({ ...current, lstartEnergyCutoffRy: numberField(event.target.value, current.lstartEnergyCutoffRy) }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Fermi integration method written to case.in2 after initialization. TETRA is the default for well-sampled solids; TEMP/TEMPS are useful for metallic or hard-to-converge occupations.">
                    Fermi method
                  </Wien2kFieldLabel>
                  <select disabled={initializationLocked} value={initialization.fermiMethod} onChange={(event) => setInitialization((current) => ({ ...current, fermiMethod: event.target.value as Wien2kInitializationSettings["fermiMethod"] }))}>
                    <option value="tetra">TETRA</option>
                    <option value="temp">TEMP</option>
                    <option value="temps">TEMPS</option>
                  </select>
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Smearing value in Ry used with TEMP or TEMPS Fermi integration. Ignored for TETRA.">
                    Smearing (Ry)
                  </Wien2kFieldLabel>
                  <input type="number" min="0" step="0.0001" disabled={initializationLocked || initialization.fermiMethod === "tetra"} placeholder="TETRA" value={initialization.fermiSmearingRy ?? ""} onChange={(event) => setInitialization((current) => ({ ...current, fermiSmearingRy: event.target.value.trim() ? numberField(event.target.value, 0.002) : null }))} />
                </label>
              </div>
              <div className="wien2k-control-grid wien2k-kmesh-grid">
                {initialization.kMesh.map((value, index) => (
                  <label key={index}>
                    <Wien2kFieldLabel tooltip={`Number of reciprocal-space sampling divisions along lattice direction ${index + 1}; denser meshes improve Brillouin-zone integration while multiplying runtime. Typical starting range: 4 to 12 per direction, scaled for cell shape.`}>
                      k{index + 1}
                    </Wien2kFieldLabel>
                    <input type="number" min="1" step="1" disabled={initializationLocked} value={value} onChange={(event) => setInitialization((current) => {
                      const next = [...current.kMesh] as [number, number, number];
                      next[index] = numberField(event.target.value, value);
                      return { ...current, kMesh: next };
                    })} />
                  </label>
                ))}
              </div>
              {initError && <p className="wien2k-validation">{initError}</p>}
            </>
          ), { status: initializationLocked ? "Locked" : undefined })}

          {renderSection("magnetism", "Magnetism and Spin", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="Selects nonmagnetic run_lapw or collinear spin-polarized initialization plus runsp_lapw. DFT+U requires the spin-polarized path in WIEN2k.">
                    Spin mode
                  </Wien2kFieldLabel>
                  <select disabled={initializationLocked || runSettings.dftU.enabled} value={initialization.spinMode} onChange={(event) => setSpinMode(event.target.value as Wien2kInitializationSettings["spinMode"])}>
                    <option value="non_spin_polarized">Non-spin-polarized</option>
                    <option value="spin_polarized">Spin-polarized</option>
                  </select>
                </label>
              </div>
              {initialization.spinMode === "spin_polarized" && (
                <table className="wien2k-site-table">
                  <thead><tr><th>Site</th><th>Start</th><th>Moment (uB)</th></tr></thead>
                  <tbody>
                    {initialization.startingMagnetization.map((entry, index) => (
                      <tr key={`${entry.siteIndex}-${entry.element}`}>
                        <td>{entry.siteIndex}. {entry.element || `Site ${entry.siteIndex}`}</td>
                        <td>
                          <select disabled={initializationLocked} value={entry.configuration} onChange={(event) => setInitialization((current) => {
                            const next = [...current.startingMagnetization];
                            next[index] = { ...entry, configuration: event.target.value as typeof entry.configuration };
                            return { ...current, startingMagnetization: next };
                          })}>
                            <option value="up">Up</option>
                            <option value="down">Down</option>
                            <option value="non_magnetic">Non-magnetic</option>
                          </select>
                        </td>
                        <td>
                          <input type="number" min="0" step="0.1" disabled={initializationLocked} value={entry.momentBohrMagneton} onChange={(event) => setInitialization((current) => {
                            const next = [...current.startingMagnetization];
                            next[index] = { ...entry, momentBohrMagneton: numberField(event.target.value, entry.momentBohrMagneton) };
                            return { ...current, startingMagnetization: next };
                          })} />
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              )}
              {initialization.spinMode === "spin_polarized" && (
                <p className="wien2k-validation">
                  Spin-polarized initialization is applied natively. Site-level AFM/down/non-magnetic occupation editing is saved with the run settings, but automatic case.inst rewriting is not enabled in this workflow yet.
                </p>
              )}
            </>
          ), { status: initializationLocked ? "Locked" : undefined })}

          {renderSection("dftu", "DFT+U Hubbard Corrections", (
            <>
              <div className="wien2k-dftu-header-grid">
                <label className="wien2k-inline-toggle">
                  <input type="checkbox" disabled={initializationLocked} checked={runSettings.dftU.enabled} onChange={(event) => enableRecommendedDftU(event.target.checked)} />
                  <span className="wien2k-field-label-row">
                    Enable DFT+U
                    <InfoTooltip text="WIEN2k applies DFT+U through ORB. QCortado uses the supported runsp_lapw -orb path, so enabling this switches the case to spin-polarized initialization." />
                  </span>
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="WIEN2k ORB double-counting mode. SIC is the common fully localized choice; AMF may be more appropriate for itinerant metallic states.">
                    Double counting
                  </Wien2kFieldLabel>
                  <select disabled={initializationLocked || !runSettings.dftU.enabled} value={runSettings.dftU.doubleCounting} onChange={(event) => setRunSettings((current) => ({ ...current, dftU: { ...current.dftU, doubleCounting: event.target.value as typeof current.dftU.doubleCounting } }))}>
                    <option value="sic">SIC</option>
                    <option value="amf">AMF</option>
                    <option value="hmf">HMF</option>
                  </select>
                </label>
              </div>
              {runSettings.dftU.enabled && (
                <table className="wien2k-site-table">
                  <thead>
                    <tr>
                      <th>Use</th>
                      <th>Site</th>
                      <th>Manifold</th>
                      <th>U (eV) <InfoTooltip text="On-site Hubbard U for this WIEN2k orbital channel." /></th>
                      <th>J (eV) <InfoTooltip text="Hund exchange J; keep 0 for a U-only correction." /></th>
                    </tr>
                  </thead>
                  <tbody>
                    {recommendedHubbardTargets.map((target) => {
                      const active = runSettings.dftU.targets.find((entry) => entry.siteIndex === target.siteIndex && entry.manifold === target.manifold);
                      return (
                        <tr key={`${target.siteIndex}-${target.manifold}`}>
                          <td><input type="checkbox" disabled={initializationLocked} checked={Boolean(active)} onChange={(event) => toggleHubbardTarget(target, event.target.checked)} /></td>
                          <td title={target.reason ?? undefined}>{target.siteIndex}. {target.element}</td>
                          <td>
                            <input value={active?.manifold ?? target.manifold} disabled={initializationLocked || !active} onChange={(event) => updateHubbardTarget(target.siteIndex, {
                              manifold: event.target.value,
                              orbitalL: orbitalLFromManifold(event.target.value),
                            })} />
                          </td>
                          <td><input type="number" min="0" step="0.1" disabled={initializationLocked || !active} value={active?.uEv ?? target.uEv} onChange={(event) => updateHubbardTarget(target.siteIndex, { uEv: numberField(event.target.value, target.uEv) })} /></td>
                          <td><input type="number" min="0" step="0.1" disabled={initializationLocked || !active} value={active?.jEv ?? target.jEv} onChange={(event) => updateHubbardTarget(target.siteIndex, { jEv: numberField(event.target.value, target.jEv) })} /></td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              )}
              {runSettings.dftU.enabled && recommendedHubbardTargets.length === 0 && (
                <p className="wien2k-validation">No transition-metal d or lanthanide/actinide f manifolds were detected for this structure.</p>
              )}
              {runError && runSettings.dftU.enabled && <p className="wien2k-validation">{runError}</p>}
            </>
          ), { status: initializationLocked ? "Locked" : undefined })}

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

          {renderSection("advanced", "Advanced SCF Options", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="Selects a WIEN2k dispersion correction switch for run_lapw/runsp_lapw. DFT-D3/D4 require the corresponding WIEN2k-side executable to be installed.">
                    Dispersion
                  </Wien2kFieldLabel>
                  <select value={runSettings.dispersionCorrection} onChange={(event) => setRunSettings((current) => ({ ...current, dispersionCorrection: event.target.value as typeof current.dispersionCorrection }))}>
                    <option value="none">None</option>
                    <option value="dftd3">DFT-D3</option>
                    <option value="dftd4">DFT-D4</option>
                  </select>
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="MIXER mode written to case.inm. MSR1 is the robust default; PRATT is mainly useful for difficult starts with small greed.">
                    Mixer
                  </Wien2kFieldLabel>
                  <select value={runSettings.mixer.mode} onChange={(event) => setRunSettings((current) => ({ ...current, mixer: { ...current.mixer, mode: event.target.value as typeof current.mixer.mode } }))}>
                    <option value="MSR1">MSR1</option>
                    <option value="MSEC3">MSEC3</option>
                    <option value="MSEC4">MSEC4</option>
                    <option value="MSR2">MSR2</option>
                    <option value="PRATT">PRATT</option>
                    <option value="PRAT0">PRAT0</option>
                  </select>
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Mixing greed Q. Keep 0.2 for most MSR1/MSEC runs; lower values may help unstable magnetic or correlated starts.">
                    Mixer greed
                  </Wien2kFieldLabel>
                  <input type="number" min="0.01" max="1" step="0.01" value={runSettings.mixer.greed} onChange={(event) => setRunSettings((current) => ({ ...current, mixer: { ...current.mixer, greed: numberField(event.target.value, current.mixer.greed) } }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Number of previous iterations used by multisecant mixing. 6 to 10 is typical; larger cells sometimes benefit from more history.">
                    Mixer history
                  </Wien2kFieldLabel>
                  <input type="number" min="1" step="1" value={runSettings.mixer.history} onChange={(event) => setRunSettings((current) => ({ ...current, mixer: { ...current.mixer, history: numberField(event.target.value, current.mixer.history) } }))} />
                </label>
                <label>
                  <Wien2kFieldLabel tooltip="Optional MIXER trust hint. STIFF/STIFFER can help difficult oscillating cases; FAST can accelerate easy cases.">
                    Trust
                  </Wien2kFieldLabel>
                  <select value={runSettings.mixer.trust} onChange={(event) => setRunSettings((current) => ({ ...current, mixer: { ...current.mixer, trust: event.target.value as typeof current.mixer.trust } }))}>
                    <option value="default">Default</option>
                    <option value="STIFF">STIFF</option>
                    <option value="STIFFER">STIFFER</option>
                    <option value="FAST">FAST</option>
                  </select>
                </label>
              </div>
              <div className="wien2k-control-grid">
                <label className="checkbox-label">
                  <input type="checkbox" checked={runSettings.iterativeDiagonalization} onChange={(event) => setRunSettings((current) => ({ ...current, iterativeDiagonalization: event.target.checked }))} />
                  <span>Iterative diagonalization</span>
                </label>
                <label className="checkbox-label">
                  <input type="checkbox" checked={runSettings.forceMinimization} onChange={(event) => setRunSettings((current) => ({ ...current, forceMinimization: event.target.checked }))} />
                  <span>MSR1a force minimization</span>
                </label>
              </div>
            </>
          ), { locked: !initialized, status: initialized ? undefined : "Locked - run initialization first" })}

          {renderSection("hpc", "HPC Run Settings", (
            <>
              <div className="wien2k-control-grid">
                <label>
                  <span className="wien2k-field-label-row">
                    Parallelization mode
                    <InfoTooltip text="OpenMP runs one WIEN2k process on one node and uses CPUs / Task as threads. K-point mode runs WIEN2k with -p, creates .machines from the Slurm allocation, and distributes independent k points across Tasks and Nodes; CPUs / Task controls OpenMP threads within each worker. Speedup is limited by the number of irreducible k points." />
                  </span>
                  <select
                    value={runSettings.parallelMode}
                    onChange={(event) => setParallelMode(event.target.value === "kpoint" ? "kpoint" : "openmp")}
                    disabled={isRunning}
                  >
                    <option value="openmp">OpenMP (single node)</option>
                    <option value="kpoint">K-point parallel (-p)</option>
                  </select>
                </label>
              </div>
              <HpcRunSettings
                profileId={activeHpcProfile?.id ?? null}
                profileName={activeHpcProfile?.name ?? "Andromeda"}
                taskKind="wien2k-scf"
                commandLines={hpcCommandLines}
                resources={hpcResources}
                resourceMode="cpu_only"
                defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
                defaultGpuResources={null}
                resourceModeMessage={runSettings.parallelMode === "kpoint"
                  ? "K-point mode distributes WIEN2k workers across Slurm tasks and nodes using .machines."
                  : "OpenMP mode uses one CPU Slurm task on one node; CPUs / Task sets the thread count."}
                lockResourceTypeWhenSingleMode={false}
                maxTasks={runSettings.parallelMode === "openmp" ? {
                  value: 1,
                  reason: "OpenMP mode uses one Slurm task",
                } : null}
                onResourcesChange={(next) => {
                  const nodes = runSettings.parallelMode === "openmp" ? 1 : Math.max(1, next.nodes ?? 1);
                  setHpcResources({
                    ...next,
                    resource_type: "cpu",
                    nodes,
                    ntasks: runSettings.parallelMode === "openmp" ? 1 : Math.max(nodes, next.ntasks ?? 1),
                    gpus: 0,
                  });
                }}
                disabled={isRunning}
                nodesDisabled={runSettings.parallelMode === "openmp"}
                tasksDisabled={runSettings.parallelMode === "openmp"}
              />
            </>
          ), { locked: !initialized, status: initialized ? undefined : "Locked - run initialization first" })}
        </section>
        <section className="wien2k-output-column">
          <LiveOutputPanel
            title="WIEN2k initialization and SCF"
            output={initializationOutputText}
            totalLineCount={initializationOutputLineCount}
            visibleLineCount={initializationOutputLines.length}
            panelClassName="output-panel wien2k-inline-output"
            outputClassName="output-text wien2k-inline-output-text"
          />
        </section>
      </div>
      <div className="step-actions">
        {session && !result && initialized && !isContinuationSession && (
          <button className="secondary-button" type="button" disabled={isRunning || initializationRunning} onClick={() => void resetInitialization()}>
            Restart Initialization
          </button>
        )}
        {!initialized && (
          <button type="button" className="primary-button" disabled={!selectedSource || Boolean(initError) || initializationRunning} onClick={() => void beginInitialization()}>
            {initializationRunning ? "Initializing..." : "Run Initialization"}
          </button>
        )}
        {initialized && !result && (
          <button
            type="button"
            className="primary-button"
            disabled={Boolean(runError) || isRunning || isPreparingContinuation}
            onClick={() => void submitScf(isContinuationSession)}
          >
            {isRunning ? (isContinuationSession ? "Continuing..." : "Running SCF...") : isContinuationSession ? "Continue SCF" : "Run SCF"}
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
