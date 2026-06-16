import { useCallback, useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import { AppHeaderPortal } from "../AppHeaderPortal";
import type { BandData } from "../BandPlot";
import { BrillouinZoneViewer, type KPathPoint } from "../BrillouinZoneViewer";
import { ElapsedTimer } from "../ElapsedTimer";
import { HpcRunSettings } from "../HpcRunSettings";
import { LiveOutputPanel } from "../LiveOutputPanel";
import { ProgressBar } from "../ProgressBar";
import { RemoteUtilizationPanel } from "../RemoteUtilizationPanel";
import { Wien2kFieldLabel } from "./Wien2kFieldLabel";
import {
  discardWien2kBandsSession,
  prepareWien2kBandsSession,
  startWien2kBandsSession,
  type Wien2kBandsExecutionResult,
  type Wien2kBandsRunSettings,
  type Wien2kBandsSession,
  type Wien2kBandsSpinChannel,
} from "../../lib/engines/wien2k";
import {
  applyWien2kBandsTotalKPoints,
  getWien2kBandProjectionOptions,
  getWien2kBandProjectionOptionsFromSites,
  transformWien2kKPathForKlistBand,
  type Wien2kBandProjectionOption,
  type Wien2kBandProjectionSite,
} from "../../lib/wien2kBandsWizard";
import { defaultCpuResources } from "../../lib/hpcConfig";
import type { CrystalData, HpcProfile, SlurmResourceRequest } from "../../lib/types";
import { useViewportScrollLock } from "../../lib/useViewportScrollLock";
import { formatCalculationSourceLabel, getCalculationName } from "../../lib/calculationNames";
import { getCalculationTagBadges, getCalcTagClass } from "../../lib/calculationTags";
import { useTaskContext } from "../../lib/TaskContext";

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

interface Wien2kBandsWizardProps {
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  activeHpcProfile?: HpcProfile | null;
  reconnectTaskId?: string;
  onBack: () => void;
  onViewBands: (
    bandData: BandData,
    scfFermiEnergy: number | null,
    calculationParameters?: Record<string, unknown> | null,
    calculationContext?: { projectId: string; cifId: string; calcId: string } | null,
  ) => void;
}

type BandsWizardStep = "source" | "kpath" | "prepare" | "run" | "results";
type SectionKey = "source" | "files" | "projections" | "spaghetti" | "logging" | "hpc";

const STEPS: Array<{ id: BandsWizardStep; label: string }> = [
  { id: "source", label: "Source" },
  { id: "kpath", label: "K-path" },
  { id: "prepare", label: "Prepare" },
  { id: "run", label: "Run" },
  { id: "results", label: "Results" },
];
const MAX_VIEWER_POINTS_PER_SEGMENT = 400;
const MAX_TOTAL_K_POINTS = 5000;

function cloneWien2kBandsResources(profile: HpcProfile | null | undefined): SlurmResourceRequest {
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

function buildWien2kBandsOpenMpPreview(resources: SlurmResourceRequest): string {
  const openMpThreads = Math.max(1, resources.ntasks ?? 1) * Math.max(1, resources.cpus_per_task ?? 1);
  return `rm -f .machines .processes && export OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK:-${openMpThreads}}"`;
}

function parseRemoteJobId(line: string): string | null {
  const match = line.match(/Scheduled (?:utility )?job submitted:\s*(\S+)/);
  return match?.[1] ?? null;
}

function parseRemoteNode(line: string): string | null {
  const match = line.match(/^\[QCortado\] Host:\s*(.+)$/);
  return match?.[1]?.trim() ?? null;
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

function sourceHasSpinOrbit(calc: CalculationRun | null): boolean {
  return Boolean(
    calc?.parameters?.spin_orbit
    || calc?.parameters?.spinOrbit
    || calc?.parameters?.soc
    || calc?.tags?.includes("soc"),
  );
}

function clampInt(value: number, min: number, max: number): number {
  if (!Number.isFinite(value)) return min;
  return Math.min(max, Math.max(min, Math.round(value)));
}

function sourceStructureSites(calc: CalculationRun | null): Wien2kBandProjectionSite[] {
  const sourceSites = calc?.parameters?.source_structure_sites;
  if (Array.isArray(sourceSites)) {
    return sourceSites as Wien2kBandProjectionSite[];
  }
  const legacySites = calc?.parameters?.sites;
  return Array.isArray(legacySites) ? legacySites as Wien2kBandProjectionSite[] : [];
}

export function Wien2kBandsWizard({
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  activeHpcProfile = null,
  reconnectTaskId,
  onBack,
  onViewBands,
}: Wien2kBandsWizardProps) {
  const taskContext = useTaskContext();
  const sources = useMemo(
    () => scfCalculations.filter(isConvergedWien2kScf),
    [scfCalculations],
  );
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(sources[0] ?? null);
  const [step, setStep] = useState<BandsWizardStep>("source");
  const [kPath, setKPath] = useState<KPathPoint[]>([]);
  const [session, setSession] = useState<Wien2kBandsSession | null>(null);
  const [result, setResult] = useState<Wien2kBandsExecutionResult | null>(null);
  const [outputLines, setOutputLines] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [isPreparing, setIsPreparing] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [totalKPointsTarget, setTotalKPointsTarget] = useState(120);
  const [totalKPointsInput, setTotalKPointsInput] = useState("120");
  const [energyMinEv, setEnergyMinEv] = useState(-8);
  const [energyMaxEv, setEnergyMaxEv] = useState(6);
  const [characterAtom, setCharacterAtom] = useState(0);
  const [characterL, setCharacterL] = useState(0);
  const [characterScale, setCharacterScale] = useState(0.2);
  const [runLapw2Qtl, setRunLapw2Qtl] = useState(false);
  const [runIrrep, setRunIrrep] = useState(false);
  const [spinOrbit, setSpinOrbit] = useState(false);
  const [diagnosticLog, setDiagnosticLog] = useState(false);
  const [spinChannel, setSpinChannel] = useState<Wien2kBandsSpinChannel>("none");
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(() => cloneWien2kBandsResources(activeHpcProfile));
  const [expandedSections, setExpandedSections] = useState<Record<SectionKey, boolean>>({
    source: true,
    files: true,
    projections: true,
    spaghetti: true,
    logging: true,
    hpc: false,
  });
  const [calcStartTime, setCalcStartTime] = useState<string | null>(null);
  const [remoteJobId, setRemoteJobId] = useState<string | null>(null);
  const [remoteNode, setRemoteNode] = useState<string | null>(null);
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const outputUnlistenRef = useRef<UnlistenFn | null>(null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const taskResult = activeTask?.result as Wien2kBandsExecutionResult | null | undefined;
  const displayedResult = result ?? taskResult ?? null;
  const taskOutputLines = activeTask ? activeTask.output : outputLines;
  const taskOutputText = activeTask ? activeTask.outputText : outputLines.join("\n");
  const taskOutputLineCount = activeTask ? activeTask.outputLineCount : outputLines.length;
  const runIsActive = activeTask?.status === "running" || isRunning;
  const runHasFailed = activeTask?.status === "failed" || activeTask?.status === "cancelled" || Boolean(error);
  const runRemoteJobId = activeTask?.hpc.remote_job_id ?? remoteJobId;
  const runRemoteNode = activeTask?.hpc.remote_node ?? remoteNode;
  const fermiEnergy = sourceFermi(selectedScf);
  const selectedSpinMode = session?.spinMode ?? sourceSpinMode(selectedScf);
  const selectedSourceHasSpinOrbit = sourceHasSpinOrbit(selectedScf);
  const currentStepIndex = STEPS.findIndex((entry) => entry.id === step);
  const runStepActive = step === "run" && runIsActive;
  const projectionOptions = useMemo(
    () => {
      const sites = sourceStructureSites(selectedScf);
      return sites.length > 0
        ? getWien2kBandProjectionOptionsFromSites(sites, crystalData)
        : getWien2kBandProjectionOptions(crystalData);
    },
    [crystalData, selectedScf],
  );
  const projectionGroups = useMemo(() => {
    const groups = new Map<number, { atomIndex: number; label: string; orbitals: Wien2kBandProjectionOption[] }>();
    for (const option of projectionOptions) {
      if (option.kind === "atom") {
        groups.set(option.atomIndex, {
          atomIndex: option.atomIndex,
          label: option.label,
          orbitals: [],
        });
      }
    }
    for (const option of projectionOptions) {
      if (option.kind !== "orbital") continue;
      const group = groups.get(option.atomIndex) ?? {
        atomIndex: option.atomIndex,
        label: `Atom ${option.atomIndex}`,
        orbitals: [],
      };
      group.orbitals.push(option);
      groups.set(option.atomIndex, group);
    }
    return Array.from(groups.values());
  }, [projectionOptions]);
  const projectionSummary = characterAtom > 0
    ? `atom ${characterAtom} / L=${characterL}`
    : "plain bands";
  const kPathSegmentCount = useMemo(
    () => kPath.filter((point, index) => index < kPath.length - 1 && point.npoints > 0).length,
    [kPath],
  );
  const totalKPoints = useMemo(
    () => kPath.reduce((sum, point) => sum + point.npoints, 0) + (kPath.length > 0 ? 1 : 0),
    [kPath],
  );
  const minimumTotalKPoints = Math.max(1, kPathSegmentCount);
  const viewerPointsPerSegment = clampInt(
    totalKPointsTarget / minimumTotalKPoints,
    1,
    MAX_VIEWER_POINTS_PER_SEGMENT,
  );
  useViewportScrollLock(runStepActive);

  useEffect(() => () => outputUnlistenRef.current?.(), []);

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
      setError(activeTask.error ?? (activeTask.status === "cancelled" ? "Cancelled by user" : "WIEN2k bands failed."));
      setIsRunning(false);
      return;
    }
    if (activeTask.status !== "completed" || !taskResult) return;
    setResult(taskResult);
    setSession((current) => current ? { ...current, phase: taskResult.phase } : current);
    setIsRunning(false);
  }, [activeTask, taskResult]);

  useEffect(() => {
    setHpcResources(cloneWien2kBandsResources(activeHpcProfile));
  }, [activeHpcProfile?.id, activeHpcProfile?.updated_at]);

  useEffect(() => {
    if (selectedSpinMode === "spin_polarized" && spinChannel === "none") {
      setSpinChannel("up");
    }
  }, [selectedSpinMode, spinChannel]);

  useEffect(() => {
    if (!selectedSourceHasSpinOrbit) return;
    setSpinOrbit(true);
    setSpinChannel(selectedSpinMode === "spin_polarized" ? "up" : "none");
  }, [selectedSourceHasSpinOrbit, selectedSpinMode]);

  useEffect(() => {
    if (characterAtom > 0 && !runLapw2Qtl) {
      setRunLapw2Qtl(true);
    }
  }, [characterAtom, runLapw2Qtl]);

  useEffect(() => {
    setKPath((prevPath) => applyWien2kBandsTotalKPoints(prevPath, totalKPointsTarget));
  }, [totalKPointsTarget]);

  useEffect(() => {
    setTotalKPointsInput(String(totalKPointsTarget));
  }, [totalKPointsTarget]);

  const handleKPathChange = useCallback((newPath: KPathPoint[]) => {
    setKPath(applyWien2kBandsTotalKPoints(newPath, totalKPointsTarget));
  }, [totalKPointsTarget]);

  const commitTotalKPointsInput = useCallback(() => {
    const parsed = Number.parseInt(totalKPointsInput.trim(), 10);
    const fallback = Number.isFinite(totalKPointsTarget) ? totalKPointsTarget : 120;
    const committed = Number.isFinite(parsed)
      ? clampInt(parsed, minimumTotalKPoints, MAX_TOTAL_K_POINTS)
      : clampInt(fallback, minimumTotalKPoints, MAX_TOTAL_K_POINTS);
    setTotalKPointsTarget(committed);
    setTotalKPointsInput(String(committed));
  }, [minimumTotalKPoints, totalKPointsInput, totalKPointsTarget]);

  async function attachOutputListener(sessionId: string) {
    outputUnlistenRef.current?.();
    outputUnlistenRef.current = await listen<string>(`wien2k-bands-output:${sessionId}`, (event) => {
      const jobId = parseRemoteJobId(event.payload);
      if (jobId) setRemoteJobId(jobId);
      const node = parseRemoteNode(event.payload);
      if (node) setRemoteNode(node);
      setOutputLines((current) => [...current, event.payload].slice(-1500));
    });
  }

  function toggleSection(section: SectionKey) {
    setExpandedSections((current) => ({ ...current, [section]: !current[section] }));
  }

  function selectProjection(option: Wien2kBandProjectionOption) {
    setCharacterAtom(option.characterAtom);
    setCharacterL(option.characterL);
    if (option.characterAtom > 0) setRunLapw2Qtl(true);
  }

  function bandsCommandLines(): string[] {
    const spinArg = spinChannel === "none" ? "" : ` -${spinChannel === "up" ? "up" : "dn"}`;
    const spinOrbitArg = spinOrbit ? " -so" : "";
    return [
      "cd \"$SLURM_SUBMIT_DIR\"",
      ...moduleSetupLines(activeHpcProfile),
      buildWien2kBandsOpenMpPreview(hpcResources),
      "rm -f lapw1.error lapw2.error irrep.error spaghetti.error *.error",
      ...(diagnosticLog ? ["# QCortado diagnostic tails enabled"] : []),
      `x lapw1 -band${spinArg}`,
      ...(spinOrbit ? [`x lapwso${spinChannel === "none" ? "" : " -up"}`] : []),
      ...(runLapw2Qtl ? [`x lapw2 -qtl -band${spinArg}${spinOrbitArg}`] : []),
      ...(runIrrep ? [`x irrep -band${spinArg}`] : []),
      `x spaghetti${spinArg}${spinOrbitArg}`,
    ];
  }

  async function ensureSession(): Promise<Wien2kBandsSession> {
    if (session) return session;
    if (!selectedScf) throw new Error("Select a WIEN2k source SCF first.");
    const next = await startWien2kBandsSession(projectId, cifId, selectedScf.id);
    setSession(next);
    setSpinChannel(next.spinMode === "spin_polarized" ? "up" : "none");
    setSpinOrbit(next.sourceSpinOrbit);
    setOutputLines([`Staged ${next.caseName} from source SCF ${selectedScf.id.slice(0, 8)}.`]);
    await attachOutputListener(next.sessionId);
    return next;
  }

  async function prepareBands() {
    if (kPath.length < 2) return;
    setIsPreparing(true);
    setError(null);
    try {
      const activeSession = await ensureSession();
      const klistPath = transformWien2kKPathForKlistBand(kPath, crystalData);
      const prepared = await prepareWien2kBandsSession(activeSession.sessionId, {
        kPath: klistPath.map((point) => ({
          label: point.label,
          coords: point.coords,
          npoints: point.npoints,
        })),
        energyMinEv,
        energyMaxEv,
        characterAtom,
        characterL,
        characterScale,
        runLapw2Qtl,
        runIrrep,
        spinChannel,
      });
      setSession((current) => current ? { ...current, phase: prepared.phase } : current);
      setOutputLines((current) => [
        ...current,
        `[prepared ${Object.keys(prepared.artifacts).join(", ")}]`,
      ].slice(-1500));
      setExpandedSections((current) => ({ ...current, hpc: true }));
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsPreparing(false);
    }
  }

  async function runBands() {
    if (!session) return;
    setIsRunning(true);
    setError(null);
    setCalcStartTime(new Date().toISOString());
    setRemoteJobId(null);
    setRemoteNode(null);
    try {
      const runSettings: Wien2kBandsRunSettings = {
        spinChannel,
        runLapw2Qtl,
        runIrrep,
        spinOrbit,
        diagnosticLog,
      };
      const taskId = await taskContext.startTask(
        "wien2k_bands",
        {
          sessionId: session.sessionId,
          settings: runSettings,
          resources: hpcResources,
        },
        `WIEN2k Bands - ${session.caseName}`,
      );
      setActiveTaskId(taskId);
    } catch (reason) {
      setError(String(reason));
      setIsRunning(false);
    }
  }

  function startPreparedRun() {
    if (!session || runIsActive) return;
    setStep("run");
    void runBands();
  }

  async function leaveWizard(destination: "back" | "view" = "back") {
    setError(null);
    try {
      if (session && !displayedResult && !runIsActive) {
        await discardWien2kBandsSession(session.sessionId);
      }
      outputUnlistenRef.current?.();
      outputUnlistenRef.current = null;
      if (destination === "view" && displayedResult) {
        onViewBands(
          displayedResult.bandData,
          fermiEnergy,
          { engine_id: "wien2k", source_scf_id: selectedScf?.id ?? null },
          { projectId, cifId, calcId: displayedResult.calculationId },
        );
      } else {
        onBack();
      }
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
        <button className="back-btn" type="button" disabled={isPreparing} onClick={() => void leaveWizard("back")}>
          ← Exit
        </button>
        <h2>WIEN2k Bands</h2>
        <div className="step-indicator">
          {STEPS.map((entry, index) => (
            <span key={entry.id} className={index === currentStepIndex ? "active" : index < currentStepIndex ? "completed" : ""}>
              {index + 1}. {entry.label}
            </span>
          ))}
        </div>
      </AppHeaderPortal>
    );
  }

  function renderSourceStep() {
    if (sources.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No WIEN2k SCF Calculations Available</h3>
          <p className="warning-text">Band calculations require a converged WIEN2k SCF with a retained native case directory.</p>
          <button className="primary-button" type="button" onClick={onBack}>Back to Dashboard</button>
        </div>
      );
    }
    return (
      <div className="wizard-step source-step">
        <div className="source-step-header">
          <h3>Select WIEN2k Source SCF</h3>
        </div>
        <p className="step-description">Choose the converged native case to clone for the band path.</p>
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
                  {getCalculationTagBadges(scf as any, { legacyFallback: true, calcId: scf.id }).map((tag, i) => (
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
          <button className="primary-button" type="button" disabled={!selectedScf} onClick={() => setStep("kpath")}>
            Next: K-path
          </button>
        </div>
      </div>
    );
  }

  function renderKPathStep() {
    return (
      <div className="wizard-step kpath-step">
        <h3>Select K-path</h3>
        <div className="crystal-info">
          {crystalData.space_group_HM && (
            <div className="info-row">
              <span className="label">Space Group:</span>
              <span className="value">{crystalData.space_group_HM}{crystalData.space_group_IT_number && ` (#${crystalData.space_group_IT_number})`}</span>
            </div>
          )}
        </div>
        <BrillouinZoneViewer
          crystalData={crystalData}
          onPathChange={handleKPathChange}
          initialPath={kPath}
          pointsPerSegment={viewerPointsPerSegment}
        />
        <div className="kpath-sampling-panel">
          <div className="kpath-sampling-header">
            <div>
              <h4>K-Path Sampling</h4>
              <p>K-points are distributed along the full path by segment length.</p>
            </div>
            <span className="kpath-sampling-summary">
              {totalKPoints} total k-points
            </span>
          </div>

          <label className="kpath-sampling-input">
            <span>Total k-points</span>
            <input
              type="number"
              min={minimumTotalKPoints}
              max={MAX_TOTAL_K_POINTS}
              value={totalKPointsInput}
              onChange={(event) => {
                setTotalKPointsInput(event.target.value);
              }}
              onBlur={commitTotalKPointsInput}
              onKeyDown={(event) => {
                if (event.key === "Enter") {
                  event.preventDefault();
                  commitTotalKPointsInput();
                }
              }}
            />
          </label>

          <p className="kpath-sampling-note">
            {kPathSegmentCount > 0
              ? `Evenly distributed by segment length across ${kPathSegmentCount} path segment${kPathSegmentCount === 1 ? "" : "s"}.`
              : "Add at least 2 points to distribute k-points along the path."}
          </p>
        </div>
        <div className="step-actions">
          <button className="secondary-button" type="button" onClick={() => setStep("source")}>Back</button>
          <button className="primary-button" type="button" disabled={kPath.length < 2} onClick={() => setStep("prepare")}>
            Next: Prepare
          </button>
        </div>
      </div>
    );
  }

  function renderPrepareStep() {
    const pathString = kPath.map((point) => point.label).join(" → ");
    const totalKPoints = kPath.reduce((sum, point) => sum + point.npoints, 0) + (kPath.length > 0 ? 1 : 0);
    const prepared = session?.phase === "prepared" || session?.phase === "bands_complete";
    const canPrepare = kPath.length >= 2 && energyMinEv < energyMaxEv && characterScale >= 0 && !isPreparing;
    const canRun = prepared && Boolean(session) && !runIsActive && (selectedSpinMode !== "spin_polarized" || spinChannel !== "none");
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
                <p>K-path: {pathString}; {totalKPoints} points written to `case.klist_band`.</p>
                <p>Projection: {projectionSummary}.</p>
              </div>
            ))}
            {renderSection("files", "Band Preparation Files", (
              <>
                <div className="wien2k-control-grid">
                  <label>
                    <Wien2kFieldLabel tooltip="Lower plotted-energy bound in eV for WIEN2k `case.insp`. This sets the spaghetti/viewer range relative to the Fermi level; it does not change how many eigenvalues lapw1 calculates.">
                      Energy min (eV)
                    </Wien2kFieldLabel>
                    <input type="number" step="0.5" value={energyMinEv} onChange={(event) => setEnergyMinEv(numberField(event.target.value, energyMinEv))} />
                  </label>
                  <label>
                    <Wien2kFieldLabel tooltip="Upper plotted-energy bound in eV for WIEN2k `case.insp`. This controls the plotted window relative to the Fermi level, not the raw eigenvalue count.">
                      Energy max (eV)
                    </Wien2kFieldLabel>
                    <input type="number" step="0.5" value={energyMaxEv} onChange={(event) => setEnergyMaxEv(numberField(event.target.value, energyMaxEv))} />
                  </label>
                  <label>
                    <Wien2kFieldLabel tooltip="Spin channel passed to `x lapw1`, `x lapw2`, `x irrep`, and `x spaghetti`. Spin-polarized SCFs usually need separate up/down band runs.">
                      Spin channel
                    </Wien2kFieldLabel>
                    <select value={spinChannel} onChange={(event) => setSpinChannel(event.target.value as Wien2kBandsSpinChannel)}>
                      <option value="none" disabled={selectedSpinMode === "spin_polarized"}>None</option>
                      <option value="up">Up</option>
                      <option value="down">Down</option>
                    </select>
                  </label>
                </div>
                {energyMinEv >= energyMaxEv && <p className="wien2k-validation">Energy min must be below energy max.</p>}
                {selectedSpinMode === "spin_polarized" && spinChannel === "none" && (
                  <p className="wien2k-validation">Spin-polarized WIEN2k sources must run an up or down band channel.</p>
                )}
              </>
            ), { status: prepared ? "Prepared" : undefined })}
            {renderSection("projections", "Projections", (
              <div className="wien2k-projection-dropdown">
                <div className="wien2k-projection-exact-grid">
                  <label>
                    <Wien2kFieldLabel tooltip="Atom index used by spaghetti character plotting. Keep 0 for plain line bands. Set with L below when preparing character/fat-band output.">
                      Character atom
                    </Wien2kFieldLabel>
                    <input
                      type="number"
                      min="0"
                      step="1"
                      value={characterAtom}
                      onChange={(event) => setCharacterAtom(Math.max(0, Math.round(numberField(event.target.value, characterAtom))))}
                    />
                  </label>
                  <label>
                    <Wien2kFieldLabel tooltip="Angular momentum channel for spaghetti character plotting: 0=s, 1=p, 2=d, 3=f. Only used when Character atom is nonzero.">
                      Character L
                    </Wien2kFieldLabel>
                    <input
                      type="number"
                      min="0"
                      max="3"
                      step="1"
                      value={characterL}
                      onChange={(event) => setCharacterL(Math.min(3, Math.max(0, Math.round(numberField(event.target.value, characterL)))))}
                    />
                  </label>
                  <label>
                    <Wien2kFieldLabel tooltip="Symbol-size scale for spaghetti character plotting. Set 0 for plain lines.">
                      Character scale
                    </Wien2kFieldLabel>
                    <input
                      type="number"
                      min="0"
                      step="0.05"
                      value={characterScale}
                      onChange={(event) => setCharacterScale(Math.max(0, numberField(event.target.value, characterScale)))}
                    />
                  </label>
                </div>
                <div className="wien2k-projection-choice-row">
                  <label className="wien2k-projection-checkbox">
                    <input
                      type="checkbox"
                      checked={characterAtom === 0}
                      onChange={() => {
                        setCharacterAtom(0);
                        setCharacterL(0);
                      }}
                    />
                    <span>Plain bands</span>
                  </label>
                  <span className="wien2k-projection-summary">{projectionSummary}</span>
                </div>
                <div className="wien2k-projection-groups">
                  {projectionGroups.map((group) => (
                    <div className="wien2k-projection-group" key={group.atomIndex}>
                      <h5>{group.atomIndex}. {group.label}</h5>
                      <div className="wien2k-projection-checkbox-grid">
                        {group.orbitals.map((option) => {
                          const orbitalLabel = option.label.startsWith(group.label)
                            ? option.label.slice(group.label.length).trim()
                            : option.label;
                          const checked = characterAtom === option.characterAtom && characterL === option.characterL;
                          return (
                            <label className="wien2k-projection-checkbox" key={option.value}>
                              <input
                                type="checkbox"
                                checked={checked}
                                onChange={(event) => {
                                  if (event.target.checked) {
                                    selectProjection(option);
                                  } else if (checked) {
                                    setCharacterAtom(0);
                                    setCharacterL(0);
                                  }
                                }}
                              />
                              <span>{orbitalLabel || option.label}</span>
                            </label>
                          );
                        })}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            ), { status: projectionSummary })}
            {renderSection("spaghetti", "Optional Setup Scripts", (
              <>
                <label className="option-checkbox">
                  <input
                    type="checkbox"
                    checked={runLapw2Qtl}
                    onChange={(event) => {
                      setRunLapw2Qtl(event.target.checked);
                      if (!event.target.checked) {
                        setCharacterAtom(0);
                        setCharacterL(0);
                      }
                    }}
                  />
                  <span>
                    Generate character data
                    <Wien2kFieldLabel tooltip="Runs `x lapw2 -qtl -band` to generate QTL/character data for spaghetti. Enable this before using character/fat-band style output."> </Wien2kFieldLabel>
                  </span>
                </label>
                <label className="option-checkbox">
                  <input type="checkbox" checked={runIrrep} onChange={(event) => setRunIrrep(event.target.checked)} />
                  <span>
                    Compute symmetry labels
                    <Wien2kFieldLabel tooltip="Runs `x irrep -band` to compute irreducible-representation data along the band path when WIEN2k can classify states by symmetry."> </Wien2kFieldLabel>
                  </span>
                </label>
                <label className="option-checkbox">
                  <input type="checkbox" checked={spinOrbit} disabled={selectedSourceHasSpinOrbit} onChange={(event) => setSpinOrbit(event.target.checked)} />
                  <span>
                    Enable spin-orbit coupling
                    <Wien2kFieldLabel tooltip="Passes `-so` to the band commands. Use only when the source case was prepared for spin-orbit coupling and contains the required SO files."> </Wien2kFieldLabel>
                  </span>
                </label>
                {selectedSourceHasSpinOrbit && <p className="wien2k-hint">SOC is required because the selected source is a converged SOC calculation.</p>}
              </>
            ), { status: prepared ? "Prepared" : undefined })}
            {renderSection("logging", "Run Logging", (
              <div className="wien2k-control-grid">
                <label>
                  <Wien2kFieldLabel tooltip="Normal streams Slurm stdout and stderr. Diagnostic tails also prints WIEN2k dayfile, output, QTL, irrep, spaghetti, and error-file tails after each native command.">
                    Verbosity
                  </Wien2kFieldLabel>
                  <select value={diagnosticLog ? "diagnostic" : "normal"} onChange={(event) => setDiagnosticLog(event.target.value === "diagnostic")}>
                    <option value="normal">Normal</option>
                    <option value="diagnostic">Diagnostic tails</option>
                  </select>
                </label>
              </div>
            ), { status: diagnosticLog ? "Diagnostic" : "Normal" })}
            {renderSection("hpc", "HPC Run Settings", (
              <HpcRunSettings
                profileId={activeHpcProfile?.id ?? null}
                profileName={activeHpcProfile?.name ?? "Andromeda"}
                taskKind="bands"
                commandLines={bandsCommandLines()}
                resources={hpcResources}
                resourceMode="cpu_only"
                defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
                defaultGpuResources={null}
                onResourcesChange={setHpcResources}
                disabled={!prepared || isPreparing || runIsActive}
              />
            ), {
              locked: !prepared,
              status: prepared ? "Ready" : "Locked - prepare files first",
            })}
            <div className="run-actions">
              <button type="button" disabled={isPreparing} onClick={() => setStep("kpath")}>Back</button>
              <button
                type="button"
                className="primary-button"
                disabled={prepared ? !canRun : !canPrepare}
                onClick={() => {
                  if (prepared) startPreparedRun();
                  else void prepareBands();
                }}
              >
                {prepared ? "Run WIEN2k Bands" : isPreparing ? "Preparing..." : "Prepare WIEN2k Files"}
              </button>
            </div>
          </section>
          <section className="wien2k-structure-preview">
            <LiveOutputPanel
              title={isPreparing ? "Preparing..." : "Preparation log"}
              output={outputLines.join("\n")}
              placeholder="Preparation output will appear here."
              totalLineCount={outputLines.length}
              visibleLineCount={outputLines.length}
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
              <h3>{runIsActive ? "Running WIEN2k Bands" : "WIEN2k Bands Output"}</h3>
              <span className={`run-step-status-pill ${runIsActive ? "running" : runHasFailed ? "error" : "idle"}`}>
                {runIsActive ? "Live output" : runHasFailed ? "Run failed" : "Ready"}
              </span>
            </div>
            <div className="run-status-rail scf-run-status">
              <ProgressBar
                status={progressStatus}
                percent={null}
                phase={runIsActive ? "lapw1/spaghetti" : displayedResult ? "Complete" : "Ready"}
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
                placeholder="Start the WIEN2k bands run to see output."
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
            {displayedResult && (
              <div className="wien2k-summary wien2k-scf-results">
                <h3>Band Result: Complete</h3>
                <p>Parsed {displayedResult.bandData.n_bands} bands across {displayedResult.bandData.n_kpoints} k-points.</p>
                {displayedResult.diagnostics.map((diagnostic) => (
                  <p key={diagnostic} className="wien2k-validation">{diagnostic}</p>
                ))}
              </div>
            )}
            <div className="run-actions">
              <button type="button" disabled={runIsActive} onClick={() => setStep("prepare")}>Back</button>
              {!displayedResult && !runIsActive && (
                <button type="button" className="primary-button" disabled={!session} onClick={() => void runBands()}>
                  Run Band Calculation
                </button>
              )}
              {displayedResult && (
                <>
                  <button type="button" className="secondary-button" disabled={runIsActive} onClick={() => void leaveWizard("back")}>
                    Return to Project
                  </button>
                  <button type="button" className="primary-button" disabled={runIsActive} onClick={() => void leaveWizard("view")}>
                    View Bands
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
        <h3>WIEN2k Bands Complete</h3>
        {displayedResult ? (
          <div className="results-summary">
            <p>Parsed {displayedResult.bandData.n_bands} bands across {displayedResult.bandData.n_kpoints} k-points.</p>
            {displayedResult.diagnostics.map((diagnostic) => <p key={diagnostic} className="wien2k-validation">{diagnostic}</p>)}
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
          <button className="secondary-button" type="button" onClick={() => void leaveWizard("back")}>Return to Project</button>
          <button className="primary-button" type="button" disabled={!displayedResult} onClick={() => void leaveWizard("view")}>
            View Bands
          </button>
        </div>
      </div>
    );
  }

  if (step === "prepare") return renderPrepareStep();
  if (step === "run") return renderRunStep();

  return (
    <div className="wizard-container wien2k-structure-wizard wien2k-scf-wizard">
      {renderHeader()}
      {error && <div className="error-banner">{error}</div>}
      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "kpath" && renderKPathStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
