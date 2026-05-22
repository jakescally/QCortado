import { useEffect, useMemo, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { CrystalData, ExecutionMode, HpcProfile, SlurmResourceRequest } from "../lib/types";
import { sortScfByMode, ScfSortMode, getStoredSortMode, setStoredSortMode } from "../lib/engines/qe/scfSorting";
import { defaultProgressState, ProgressState } from "../lib/engines/qe/progress";
import { countVisibleOutputLines } from "../lib/liveOutput";
import { useTaskContext } from "../lib/TaskContext";
import { loadGlobalMpiDefaults } from "../lib/mpiDefaults";
import { useViewportScrollLock } from "../lib/useViewportScrollLock";
import { getHubbardEligibilityReason, getScfHubbardUDisplayValues, isDudarevDftUScf, normalizeHubbardLrtUValues } from "../lib/engines/qe/hubbard";
import { isPhononReadyScf } from "../lib/engines/qe/phononReady";
import { formatCalculationSourceLabel, getCalculationName } from "../lib/calculationNames";
import { readProjectWizardSettings, writeProjectWizardSettings } from "../lib/projectWizardSettings";
import {
  buildExecutionTarget,
  defaultResourcesForProfile,
  downloadHpcCalculationArtifacts,
  saveExecutionMode,
} from "../lib/hpcConfig";
import {
  buildHpcQeInputCommandLine,
  resolveProfileRemoteQeBinDir,
} from "../lib/engines/qe/hpc";
import { ProgressBar } from "./ProgressBar";
import { ElapsedTimer } from "./ElapsedTimer";
import { LiveOutputPanel } from "./LiveOutputPanel";
import { HpcRunSettings } from "./HpcRunSettings";
import { RemoteUtilizationPanel } from "./RemoteUtilizationPanel";
import { InfoTooltip } from "./InfoTooltip";
import type { EngineId } from "../lib/engines/types";

interface CalculationRun {
  id: string;
  engine_id?: EngineId | null;
  calc_type: string;
  parameters: any;
  result: {
    converged: boolean;
    total_energy: number | null;
    fermi_energy: number | null;
    raw_output?: string | null;
  } | null;
  started_at: string;
  completed_at: string | null;
  tags?: string[];
}

interface HubbardLrtWizardProps {
  onBack: () => void;
  onExecutionModeChange?: (mode: ExecutionMode) => Promise<void> | void;
  qePath: string;
  executionMode?: ExecutionMode;
  activeHpcProfile?: HpcProfile | null;
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  scfCalculations: CalculationRun[];
  reconnectTaskId?: string;
}

type WizardStep = "source" | "parameters" | "run" | "results";

interface StoredSettings {
  qMesh: [number, number, number];
  findAtpert: number;
  convThrChi: string;
  niterMax: number;
  iverbosity: number;
  noMetq0: boolean;
  skipEquivalenceQ: boolean;
  doccThr: string;
  ethrNscf: string;
  threshInit: string;
  alphaMix: string;
  nmix: number;
  maxSeconds: string;
}

const HUBBARD_LRT_WIZARD_SETTINGS_ID = "hubbard-lrt";
const WORK_DIR = "/tmp/qcortado_hubbard_lrt";
const FIND_ATPERT_TOOLTIP = [
  "Atomic perturbation selector. QE variable: `find_atpert`.",
  "1: distinguish inequivalent Hubbard atoms by unperturbed occupations.",
  "2: group by Hubbard atomic type only; fastest, but can miss symmetry or occupation inequivalence.",
  "3: distinguish inequivalent Hubbard atoms by symmetry.",
  "4: perturb all Hubbard atoms; most expensive.",
].join(" ");
const IVERBOSITY_TOOLTIP = [
  "Controls hp.x output detail. QE variable: `iverbosity`.",
  "1 is low detail: minimal output.",
  "2 is medium detail: includes symmetry matrices, response matrices, inverse matrices, and full U matrix.",
  "Higher QE levels exist, but the wizard keeps the UI to the commonly useful range.",
].join(" ");

function qeVariableTooltip(variable: string, description: string): string {
  return `${description} QE variable: \`${variable}\`.`;
}

function isHpcCalculation(calc: CalculationRun): boolean {
  const params = calc.parameters || {};
  const backend = String(params.execution_backend || "").trim().toLowerCase();
  if (backend === "hpc") return true;
  if (params.remote_job_id || params.remote_workdir || params.remote_project_path) return true;
  const rawOutput = typeof calc.result?.raw_output === "string" ? calc.result.raw_output : "";
  return rawOutput.includes("HPC_STAGE|") || rawOutput.includes("HPC_CMD|");
}

function hasFullScfBundle(calc: CalculationRun, downloadedIds: Set<string>): boolean {
  if (downloadedIds.has(calc.id)) return true;
  const params = calc.parameters || {};
  if (params.artifacts_downloaded_full === true) return true;
  return String(params.artifact_sync_mode || "").trim().toLowerCase() === "full";
}

function getScfProfileId(calc: CalculationRun): string | null {
  const value = calc.parameters?.hpc_profile_id;
  return typeof value === "string" && value.trim() ? value.trim() : null;
}

function formatHubbardUDisplay(calc: CalculationRun): string | null {
  const values = getScfHubbardUDisplayValues(calc);
  if (values.length === 0) return null;
  return values
    .map((entry) => `${entry.target} = ${entry.value_ev.toFixed(3)} eV`)
    .join(", ");
}

function getCalculationTags(
  calc: CalculationRun,
  downloadedIds?: Set<string>,
): { label: string; type: "info" | "feature" | "special" | "geometry" }[] {
  const tags: { label: string; type: "info" | "feature" | "special" | "geometry" }[] = [];
  const params = calc.parameters || {};
  const phononReady = calc.calc_type === "scf" && isPhononReadyScf(params, calc.tags);
  const pushTag = (label: string, type: "info" | "feature" | "special" | "geometry") => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (calc.tags) {
    for (const tag of calc.tags) {
      if (tag === "phonon-ready") {
        if (phononReady) pushTag("Phonon-Ready", "special");
      } else if (tag === "structure-optimized") {
        pushTag("Optimized", "special");
      } else if (tag === "geometry") {
        pushTag("Geometry", "geometry");
      }
    }
  }

  if (params.structure_source?.type === "optimization") {
    pushTag("Geometry", "geometry");
  }

  if (params.kgrid) {
    const [k1, k2, k3] = params.kgrid;
    pushTag(`${k1}×${k2}×${k3}`, "info");
  }

  if (params.conv_thr) {
    const threshold = Number(params.conv_thr);
    const label = Number.isFinite(threshold) && threshold < 0.001
      ? threshold.toExponential(0)
      : String(params.conv_thr);
    pushTag(label, "info");
  }

  if (phononReady) {
    pushTag("Phonon-Ready", "special");
  }

  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }

  if (params.nspin === 4) {
    pushTag("Non-collinear", "feature");
  } else if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }

  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }

  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
    if (hasFullScfBundle(calc, downloadedIds ?? new Set<string>())) {
      pushTag("Downloaded", "feature");
    }
  }

  return tags;
}

function parsePositiveNumber(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (!trimmed) return null;
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive number.`);
  }
  return parsed;
}

function parseOptionalNumber(input: string, label: string): number | null {
  const trimmed = input.trim();
  if (!trimmed) return null;
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed)) {
    throw new Error(`${label} must be numeric.`);
  }
  return parsed;
}

export function HubbardLrtWizard({
  onBack,
  onExecutionModeChange,
  executionMode = "local",
  activeHpcProfile = null,
  projectId,
  cifId,
  crystalData,
  scfCalculations,
  reconnectTaskId,
}: HubbardLrtWizardProps) {
  const stored = useMemo(
    () => readProjectWizardSettings<StoredSettings>(projectId, HUBBARD_LRT_WIZARD_SETTINGS_ID),
    [projectId],
  );
  const taskContext = useTaskContext();
  const isHpcMode = executionMode === "hpc";
  const [activeTaskId, setActiveTaskId] = useState<string | null>(reconnectTaskId ?? null);
  const activeTask = activeTaskId ? taskContext.getTask(activeTaskId) : undefined;
  const [step, setStep] = useState<WizardStep>(reconnectTaskId ? "run" : "source");
  const [error, setError] = useState<string | null>(null);
  const [selectedScf, setSelectedScf] = useState<CalculationRun | null>(null);
  const [downloadedDependencyScfIds, setDownloadedDependencyScfIds] = useState<Set<string>>(new Set());
  const [isResolvingDependency, setIsResolvingDependency] = useState(false);
  const [dependencyStatus, setDependencyStatus] = useState<string | null>(null);
  const [scfSortMode, setScfSortMode] = useState<ScfSortMode>(() => getStoredSortMode());
  const [qMesh, setQMesh] = useState<[number, number, number]>(() => stored?.qMesh ?? [2, 2, 2]);
  const [findAtpert, setFindAtpert] = useState(() => stored?.findAtpert ?? 1);
  const [convThrChi, setConvThrChi] = useState(() => stored?.convThrChi ?? "1e-5");
  const [niterMax, setNiterMax] = useState(() => stored?.niterMax ?? 100);
  const [iverbosity, setIverbosity] = useState(() => stored?.iverbosity ?? 1);
  const [noMetq0, setNoMetq0] = useState(() => stored?.noMetq0 ?? false);
  const [skipEquivalenceQ, setSkipEquivalenceQ] = useState(() => stored?.skipEquivalenceQ ?? false);
  const [doccThr, setDoccThr] = useState(() => stored?.doccThr ?? "1e-5");
  const [ethrNscf, setEthrNscf] = useState(() => stored?.ethrNscf ?? "1e-11");
  const [threshInit, setThreshInit] = useState(() => stored?.threshInit ?? "1e-14");
  const [alphaMix, setAlphaMix] = useState(() => stored?.alphaMix ?? "0.3");
  const [nmix, setNmix] = useState(() => stored?.nmix ?? 4);
  const [maxSeconds, setMaxSeconds] = useState(() => stored?.maxSeconds ?? "");
  const [expandedAdvanced, setExpandedAdvanced] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState("");
  const [outputLineCount, setOutputLineCount] = useState(0);
  const [progress, setProgress] = useState<ProgressState>(defaultProgressState("Hubbard LRT"));
  const [result, setResult] = useState<any | null>(null);
  const [isSaved, setIsSaved] = useState(false);
  const [calcStartTime, setCalcStartTime] = useState("");
  const [mpiEnabled, setMpiEnabled] = useState(false);
  const [mpiProcs, setMpiProcs] = useState(1);
  const [mpiAvailable, setMpiAvailable] = useState(false);
  const [hpcResources, setHpcResources] = useState<SlurmResourceRequest>(
    defaultResourcesForProfile(activeHpcProfile),
  );
  const visibleOutputLineCount = useMemo(() => countVisibleOutputLines(output), [output]);
  useViewportScrollLock(step === "run");

  const hasExternalRunningTask = taskContext.activeTasks.some(
    (task) => task.status === "running" && task.taskId !== activeTaskId,
  );
  const hasBlockingExternalTask = !isHpcMode && hasExternalRunningTask;

  useEffect(() => {
    writeProjectWizardSettings(projectId, HUBBARD_LRT_WIZARD_SETTINGS_ID, {
      qMesh,
      findAtpert,
      convThrChi,
      niterMax,
      iverbosity,
      noMetq0,
      skipEquivalenceQ,
      doccThr,
      ethrNscf,
      threshInit,
      alphaMix,
      nmix,
      maxSeconds,
    });
  }, [projectId, qMesh, findAtpert, convThrChi, niterMax, iverbosity, noMetq0, skipEquivalenceQ, doccThr, ethrNscf, threshInit, alphaMix, nmix, maxSeconds]);

  useEffect(() => {
    if (visibleOutputLineCount > outputLineCount) {
      setOutputLineCount(visibleOutputLineCount);
    }
  }, [outputLineCount, visibleOutputLineCount]);

  useEffect(() => {
    async function initMpi() {
      try {
        const count = await invoke<number>("get_cpu_count");
        const defaults = await loadGlobalMpiDefaults(Math.max(1, Math.floor(count)));
        const available = await invoke<boolean>("check_mpi_available");
        setMpiAvailable(available);
        setMpiEnabled(available ? defaults.enabled : false);
        setMpiProcs(defaults.nprocs);
      } catch {
        setMpiAvailable(false);
      }
    }
    void initMpi();
  }, []);

  useEffect(() => {
    if (!isHpcMode) return;
    setHpcResources(defaultResourcesForProfile(activeHpcProfile));
  }, [isHpcMode, activeHpcProfile?.id, activeHpcProfile?.resource_mode]);

  useEffect(() => {
    if (!activeTask) return;
    setIsRunning(activeTask.status === "running");
    setOutput(activeTask.outputText);
    setOutputLineCount(activeTask.outputLineCount);
    setProgress(activeTask.progress);
    setCalcStartTime(activeTask.startedAt);
    if (activeTask.status === "completed" && activeTask.result) {
      setResult(activeTask.result);
      setStep("results");
    } else if (activeTask.status === "failed" || activeTask.status === "cancelled") {
      setError(activeTask.error || "Task failed");
    } else {
      setStep("run");
    }
  }, [activeTask]);

  const hpcCommandLines = useMemo(
    () => [
      "cd \"$SLURM_SUBMIT_DIR\"",
      `QE_BIN="${resolveProfileRemoteQeBinDir(activeHpcProfile, hpcResources.resource_type)}"`,
      buildHpcQeInputCommandLine(activeHpcProfile, "hp.x", "hp.in", "hp.out", undefined, hpcResources.resource_type),
    ],
    [activeHpcProfile, hpcResources.resource_type],
  );

  const selectedScfDependencyBlocked = useMemo(() => {
    if (isHpcMode || !selectedScf) return false;
    if (!isHpcCalculation(selectedScf)) return false;
    return !hasFullScfBundle(selectedScf, downloadedDependencyScfIds);
  }, [isHpcMode, selectedScf, downloadedDependencyScfIds]);

  async function handleSwitchToHpcMode() {
    setIsResolvingDependency(true);
    setDependencyStatus("Switching execution mode to HPC...");
    setError(null);
    try {
      if (onExecutionModeChange) {
        await onExecutionModeChange("hpc");
      } else {
        await saveExecutionMode("hpc");
      }
      setDependencyStatus("Execution mode switched to HPC.");
    } catch (e) {
      setError(`Failed to switch execution mode: ${e}`);
      setDependencyStatus(null);
    } finally {
      setIsResolvingDependency(false);
    }
  }

  async function handleDownloadDependencyBundle() {
    if (!selectedScf) return;
    setIsResolvingDependency(true);
    setDependencyStatus("Downloading full SCF bundle from remote...");
    setError(null);
    try {
      const report = await downloadHpcCalculationArtifacts(
        projectId,
        selectedScf.id,
        getScfProfileId(selectedScf),
        true,
      );
      setDownloadedDependencyScfIds((prev) => new Set(prev).add(selectedScf.id));
      setSelectedScf((prev) => prev && prev.id === selectedScf.id
        ? {
          ...prev,
          parameters: {
            ...(prev.parameters || {}),
            artifacts_downloaded_full: true,
            artifact_sync_mode: "full",
            remote_storage_bytes: report.downloaded_bytes + report.skipped_bytes,
          },
        }
        : prev);
      setDependencyStatus("Full SCF bundle downloaded. Local run is now available.");
    } catch (e) {
      setError(`Failed to download full SCF bundle: ${e}`);
      setDependencyStatus(null);
    } finally {
      setIsResolvingDependency(false);
    }
  }

  function buildConfig() {
    if (!selectedScf) throw new Error("Select a Dudarev DFT+U source SCF.");
    const hp = {
      prefix: selectedScf.parameters?.prefix || "qcortado_scf",
      outdir: "./tmp",
      nq: qMesh,
      find_atpert: findAtpert,
      conv_thr_chi: parsePositiveNumber(convThrChi, "conv_thr_chi"),
      niter_max: niterMax,
      iverbosity,
      no_metq0: noMetq0,
      skip_equivalence_q: skipEquivalenceQ,
      docc_thr: parsePositiveNumber(doccThr, "docc_thr"),
      ethr_nscf: parsePositiveNumber(ethrNscf, "ethr_nscf"),
      thresh_init: parsePositiveNumber(threshInit, "thresh_init"),
      alpha_mix: parsePositiveNumber(alphaMix, "alpha_mix(1)"),
      nmix,
      max_seconds: parseOptionalNumber(maxSeconds, "max_seconds"),
    };
    return {
      hp,
      project_id: projectId,
      scf_calc_id: selectedScf.id,
    };
  }

  function buildSaveParameters(config: any) {
    const scfParams = selectedScf?.parameters || {};
    return {
      source_scf_id: selectedScf?.id,
      prefix: config.hp.prefix,
      q_mesh: qMesh,
      find_atpert: findAtpert,
      conv_thr_chi: config.hp.conv_thr_chi,
      niter_max: niterMax,
      iverbosity,
      no_metq0: noMetq0,
      skip_equivalence_q: skipEquivalenceQ,
      docc_thr: config.hp.docc_thr,
      ethr_nscf: config.hp.ethr_nscf,
      thresh_init: config.hp.thresh_init,
      alpha_mix: config.hp.alpha_mix,
      nmix,
      max_seconds: config.hp.max_seconds,
      source_hubbard_parameters: scfParams.hubbard_parameters ?? null,
      hubbard_projector: scfParams.hubbard_projector ?? null,
      hubbard_formulation: scfParams.hubbard_formulation ?? null,
    };
  }

  async function runCalculation() {
    if (selectedScfDependencyBlocked) {
      setError("Selected SCF was computed remotely and needs a full local bundle for local execution.");
      return;
    }
    if (hasBlockingExternalTask) {
      setError("Another local task is currently running. Queue this task or wait for completion.");
      return;
    }
    setIsRunning(true);
    setOutput("");
    setOutputLineCount(0);
    setError(null);
    setResult(null);
    setIsSaved(false);
    setProgress(defaultProgressState("Hubbard LRT"));
    const startTime = new Date().toISOString();
    setCalcStartTime(startTime);
    setStep("run");

    try {
      const config = buildConfig();
      const saveParameters = buildSaveParameters(config);
      const taskParams = {
        config,
        workingDir: WORK_DIR,
        mpiConfig: !isHpcMode && mpiEnabled ? { enabled: true, nprocs: mpiProcs } : null,
        executionTarget: buildExecutionTarget(
          executionMode,
          activeHpcProfile?.id ?? null,
          isHpcMode ? hpcResources : null,
          false,
        ),
      };
      const hpcSaveSpec = isHpcMode
        ? {
          projectId,
          cifId,
          workingDir: WORK_DIR,
          calcType: "hubbard_lrt" as const,
          parameters: saveParameters,
          tags: ["Hubbard-LRT"],
          inputContent: "",
        }
        : null;
      const label = `Hubbard LRT - ${crystalData?.formula_sum || ""}`;
      const taskId = await taskContext.startTask("hubbard_lrt", taskParams, label, hpcSaveSpec);
      setActiveTaskId(taskId);
      const finalTask = await taskContext.waitForTaskCompletion(taskId);
      if (finalTask.status !== "completed" || !finalTask.result) {
        throw new Error(finalTask.error || "Hubbard LRT failed");
      }
      const taskResult = finalTask.result;
      const outputContent = finalTask.output.join("\n");
      const endTime = new Date().toISOString();
      const hpcSaveParams = (isHpcMode || finalTask.hpc.backend === "hpc")
        ? {
          execution_backend: "hpc",
          hpc_profile_id: finalTask.hpc.hpc_profile_id ?? activeHpcProfile?.id ?? null,
          hpc_resource_type: finalTask.hpc.hpc_resource_type ?? hpcResources.resource_type,
          remote_job_id: finalTask.hpc.remote_job_id ?? null,
          scheduler_state: finalTask.hpc.scheduler_state ?? null,
          remote_node: finalTask.hpc.remote_node ?? null,
          remote_workdir: finalTask.hpc.remote_workdir ?? null,
          remote_project_path: finalTask.hpc.remote_project_path ?? null,
          remote_storage_bytes: finalTask.hpc.remote_storage_bytes ?? null,
        }
        : {};
      setResult(taskResult);
      setStep("results");
      setProgress({ ...finalTask.progress, status: "complete", percent: 100, phase: "Complete" });

      try {
        await invoke("save_calculation", {
          projectId,
          cifId,
          calcData: {
            calc_type: "hubbard_lrt",
            parameters: {
              ...saveParameters,
              u_value_count: Array.isArray(taskResult.u_values) ? taskResult.u_values.length : 0,
              ...hpcSaveParams,
            },
            result: {
              converged: Boolean(taskResult.converged),
              total_energy: null,
              fermi_energy: null,
              n_scf_steps: null,
              wall_time_seconds: null,
              raw_output: taskResult.raw_output || outputContent,
              hubbard_lrt_data: taskResult,
            },
            started_at: startTime,
            completed_at: endTime,
            input_content: "",
            output_content: outputContent,
            tags: ["Hubbard-LRT"],
          },
          workingDir: finalTask.hpc.local_sync_dir ?? WORK_DIR,
        });
        setIsSaved(true);
      } catch (saveError) {
        setError(`Failed to auto-save Hubbard LRT calculation: ${saveError}`);
      }
    } catch (e) {
      setError(String(e));
      setProgress((prev) => ({ ...prev, status: "error", phase: "Error" }));
    } finally {
      setIsRunning(false);
    }
  }

  function updateQMesh(index: number, value: number) {
    setQMesh((prev) => {
      const next = [...prev] as [number, number, number];
      next[index] = Math.max(1, Math.round(value || 1));
      return next;
    });
  }

  function renderSourceStep() {
    const scfs = sortScfByMode(
      scfCalculations.filter((calc) => calc.calc_type === "scf"),
      scfSortMode,
    );
    const eligibleCount = scfs.filter(isDudarevDftUScf).length;
    if (scfs.length === 0) {
      return (
        <div className="wizard-step source-step">
          <h3>No SCF Calculations Available</h3>
          <p className="warning-text">Hubbard LRT requires a completed Dudarev DFT+U SCF calculation.</p>
          <button className="primary-button" onClick={onBack}>Back to Dashboard</button>
        </div>
      );
    }

    return (
      <div className="wizard-step source-step">
        <div className="source-step-header">
          <h3>Select Source SCF Calculation</h3>
          <div className="source-sort-control">
            <label htmlFor="hubbard-lrt-scf-sort">Sort SCFs</label>
            <select
              id="hubbard-lrt-scf-sort"
              value={scfSortMode}
              onChange={(event) => {
                const mode = event.target.value as ScfSortMode;
                setScfSortMode(mode);
                setStoredSortMode(mode);
              }}
            >
              <option value="recent">Most Recent</option>
              <option value="best">Best</option>
            </select>
          </div>
        </div>
        <p className="step-description">
          Choose a completed Dudarev DFT+U SCF. DFT+U+J and DFT+U+J0 are intentionally disabled for hp.x.
        </p>
        {eligibleCount === 0 && (
          <p className="warning-text">No eligible Dudarev DFT+U SCF is available yet.</p>
        )}
        <div className="scf-list">
          {scfs.map((scf) => {
            const reason = getHubbardEligibilityReason(scf);
            const disabled = Boolean(reason);
            const hubbardUDisplay = formatHubbardUDisplay(scf);
            const scfName = getCalculationName(scf);
            return (
              <div
                key={scf.id}
                className={`scf-option ${selectedScf?.id === scf.id ? "selected" : ""} ${disabled ? "disabled" : ""}`}
                onClick={() => {
                  if (!disabled) {
                    setSelectedScf(scf);
                    setDependencyStatus(null);
                  }
                }}
              >
                <div className="scf-option-header">
                  <input
                    type="radio"
                    checked={selectedScf?.id === scf.id}
                    disabled={disabled}
                    onChange={() => {
                      if (!disabled) setSelectedScf(scf);
                    }}
                  />
                  {scfName && (
                    <span className="scf-name" title={formatCalculationSourceLabel(scf)}>
                      {scfName}
                    </span>
                  )}
                  <span className="scf-date">{new Date(scf.started_at).toLocaleDateString()}</span>
                </div>
                <div className="scf-details">
                  <span>E = {scf.result?.total_energy?.toFixed(6) ?? "N/A"} Ry</span>
                  {scf.result?.fermi_energy != null && (
                    <span>E_F = {scf.result.fermi_energy.toFixed(3)} eV</span>
                  )}
                  {hubbardUDisplay && (
                    <span>Hubbard U: {hubbardUDisplay}</span>
                  )}
                  {scf.parameters?.hubbard_formulation != null && (
                    <span>Formulation = {Number(scf.parameters.hubbard_formulation) === 0 ? "DFT+U" : Number(scf.parameters.hubbard_formulation) === 1 ? "DFT+U+J" : "DFT+U+J0"}</span>
                  )}
                </div>
                <div className="calc-tags">
                  {getCalculationTags(scf, downloadedDependencyScfIds).map((tag, i) => (
                    <span
                      key={i}
                      className={`calc-tag calc-tag-${tag.type}${tag.label.trim().toUpperCase() === "HPC" ? " calc-tag-hpc" : ""}${tag.label.trim().toUpperCase() === "DOWNLOADED" ? " calc-tag-downloaded" : ""}`}
                    >
                      {tag.label}
                    </span>
                  ))}
                </div>
                {reason && <p className="param-hint">{reason}</p>}
              </div>
            );
          })}
        </div>
        <div className="step-actions">
          <button className="secondary-button" onClick={onBack}>Cancel</button>
          <button className="primary-button" disabled={!selectedScf} onClick={() => setStep("parameters")}>
            Next: Parameters
          </button>
        </div>
      </div>
    );
  }

  function renderParametersStep() {
    return (
      <div className="wizard-step parameters-step">
        <h3>Configure hp.x</h3>
        <p className="step-description">Set the q mesh and convergence controls for QE linear-response Hubbard U.</p>
        <div className="param-grid">
          <div className="param-row">
            <label>
              Q Mesh <InfoTooltip text={qeVariableTooltip("nq1, nq2, nq3", "Uniform q mesh for the linear-response perturbations.")} />
            </label>
            <div className="kgrid-inputs">
              <input type="number" min={1} value={qMesh[0]} onChange={(e) => updateQMesh(0, Number(e.target.value))} />
              <span>×</span>
              <input type="number" min={1} value={qMesh[1]} onChange={(e) => updateQMesh(1, Number(e.target.value))} />
              <span>×</span>
              <input type="number" min={1} value={qMesh[2]} onChange={(e) => updateQMesh(2, Number(e.target.value))} />
            </div>
          </div>
          <div className="param-row">
            <label>
              Perturbation Selection <InfoTooltip text={FIND_ATPERT_TOOLTIP} />
            </label>
            <input
              type="number"
              min={1}
              max={4}
              value={findAtpert}
              onChange={(e) => setFindAtpert(Math.min(4, Math.max(1, Math.round(Number(e.target.value) || 1))))}
            />
          </div>
          <div className="param-row">
            <label>
              Response Convergence <InfoTooltip text={qeVariableTooltip("conv_thr_chi", "Convergence threshold for the self-consistent response function.")} />
            </label>
            <input type="text" value={convThrChi} onChange={(e) => setConvThrChi(e.target.value)} />
          </div>
          <div className="param-row">
            <label>
              Max Iterations <InfoTooltip text={qeVariableTooltip("niter_max", "Maximum number of response self-consistency iterations.")} />
            </label>
            <input type="number" min={1} value={niterMax} onChange={(e) => setNiterMax(Math.max(1, Math.round(Number(e.target.value) || 1)))} />
          </div>
          <div className="param-row">
            <label>
              Output Detail <InfoTooltip text={IVERBOSITY_TOOLTIP} />
            </label>
            <select value={iverbosity} onChange={(e) => setIverbosity(Number(e.target.value))}>
              <option value={1}>Low - minimal output</option>
              <option value={2}>Medium - response matrices and full U matrix</option>
            </select>
          </div>
          <div className="param-row">
            <label className="toggle-label">
              <input type="checkbox" checked={noMetq0} onChange={(e) => setNoMetq0(e.target.checked)} />
              <span>
                Skip Metallic q=0 Term <InfoTooltip text={qeVariableTooltip("no_metq0", "Controls whether hp.x skips the special metallic q=0 contribution.")} />
              </span>
            </label>
          </div>
        </div>
        <section className="option-section config-section collapsible">
          <h4 className="section-header" onClick={() => setExpandedAdvanced((prev) => !prev)}>
            <span className={`collapse-icon ${expandedAdvanced ? "expanded" : ""}`}>▶</span>
            Advanced
          </h4>
          {expandedAdvanced && (
            <div className="option-params param-grid">
              <div className="param-row">
                <label className="toggle-label">
                  <input type="checkbox" checked={skipEquivalenceQ} onChange={(e) => setSkipEquivalenceQ(e.target.checked)} />
                  <span>
                    Skip Q-Point Equivalence <InfoTooltip text={qeVariableTooltip("skip_equivalence_q", "Disables symmetry equivalence reduction among q points.")} />
                  </span>
                </label>
              </div>
              <div className="param-row">
                <label>
                  Occupation Change Threshold <InfoTooltip text={qeVariableTooltip("docc_thr", "Threshold for changes in Hubbard occupations during response iterations.")} />
                </label>
                <input type="text" value={doccThr} onChange={(e) => setDoccThr(e.target.value)} />
              </div>
              <div className="param-row">
                <label>
                  NSCF Energy Threshold <InfoTooltip text={qeVariableTooltip("ethr_nscf", "Electronic convergence threshold for the non-self-consistent response calculations.")} />
                </label>
                <input type="text" value={ethrNscf} onChange={(e) => setEthrNscf(e.target.value)} />
              </div>
              <div className="param-row">
                <label>
                  Initial Response Threshold <InfoTooltip text={qeVariableTooltip("thresh_init", "Initial threshold used by hp.x before response convergence begins.")} />
                </label>
                <input type="text" value={threshInit} onChange={(e) => setThreshInit(e.target.value)} />
              </div>
              <div className="param-row">
                <label>
                  Mixing Factor <InfoTooltip text={qeVariableTooltip("alpha_mix(1)", "First response mixing factor.")} />
                </label>
                <input type="text" value={alphaMix} onChange={(e) => setAlphaMix(e.target.value)} />
              </div>
              <div className="param-row">
                <label>
                  Mixing History <InfoTooltip text={qeVariableTooltip("nmix", "Number of previous response iterations kept for mixing.")} />
                </label>
                <input type="number" min={1} value={nmix} onChange={(e) => setNmix(Math.max(1, Math.round(Number(e.target.value) || 1)))} />
              </div>
              <div className="param-row">
                <label>
                  Time Limit <InfoTooltip text={qeVariableTooltip("max_seconds", "Optional wall-time limit passed to hp.x.")} />
                </label>
                <input type="text" placeholder="none" value={maxSeconds} onChange={(e) => setMaxSeconds(e.target.value)} />
              </div>
            </div>
          )}
        </section>
        {selectedScfDependencyBlocked && (
          <div className="dependency-warning">
            <strong>Remote SCF dependency detected.</strong> This source SCF needs a full local bundle for local Hubbard LRT.
            <div className="dependency-actions">
              <button className="secondary-button" onClick={() => void handleSwitchToHpcMode()} disabled={isResolvingDependency}>
                Switch to HPC
              </button>
              <button className="primary-button" onClick={() => void handleDownloadDependencyBundle()} disabled={isResolvingDependency}>
                {isResolvingDependency ? "Downloading..." : "Download Full Bundle"}
              </button>
            </div>
            {dependencyStatus && <div className="param-hint">{dependencyStatus}</div>}
          </div>
        )}
        {isHpcMode ? (
          <HpcRunSettings
            profileId={activeHpcProfile?.id ?? null}
            profileName={activeHpcProfile?.name ?? "Andromeda"}
            taskKind="hubbard_lrt"
            commandLines={hpcCommandLines}
            resources={hpcResources}
            resourceMode={activeHpcProfile?.resource_mode ?? "both"}
            defaultCpuResources={activeHpcProfile?.default_cpu_resources ?? null}
            defaultGpuResources={activeHpcProfile?.default_gpu_resources ?? null}
            onResourcesChange={setHpcResources}
            disabled={isRunning}
          />
        ) : (
          <div className="mpi-section">
            <h4>Parallelization</h4>
            <label className="toggle-label">
              <input type="checkbox" checked={mpiEnabled} disabled={!mpiAvailable} onChange={(e) => setMpiEnabled(e.target.checked)} />
              <span>{mpiAvailable ? "Use MPI" : "MPI unavailable"}</span>
            </label>
            {mpiEnabled && (
              <input type="number" min={1} value={mpiProcs} onChange={(e) => setMpiProcs(Math.max(1, Math.round(Number(e.target.value) || 1)))} />
            )}
          </div>
        )}
        <div className="step-actions">
          <button className="secondary-button" onClick={() => setStep("source")}>Back</button>
          <button className="primary-button" disabled={!selectedScf || selectedScfDependencyBlocked || hasBlockingExternalTask} onClick={() => void runCalculation()}>
            Run Hubbard LRT
          </button>
        </div>
      </div>
    );
  }

  function renderRunStep() {
    return (
      <div className="wizard-step run-step run-step-focused">
        <div className="run-step-headline">
          <h3>{isRunning ? "Running Hubbard LRT" : "Hubbard LRT Output"}</h3>
          <span className={`run-step-status-pill ${isRunning ? "running" : error ? "error" : "idle"}`}>
            {isRunning ? "Live output" : error ? "Run failed" : "Output"}
          </span>
        </div>

        <div className="run-status-rail">
          <ProgressBar
            status={progress.status}
            percent={progress.percent}
            phase={progress.phase}
            detail={progress.detail}
            compact
          />
          <div className="run-status-meta">
            {calcStartTime && <ElapsedTimer startedAt={calcStartTime} isRunning={isRunning} />}
          </div>
        </div>

        {error && (
          <>
            <div className="run-inline-error">{error}</div>
            <div className="run-error-actions">
              <button className="secondary-button" onClick={() => setStep("parameters")}>
                Back to Parameters
              </button>
            </div>
          </>
        )}

        <div className={`run-layout ${isHpcMode ? "run-layout-hpc-telemetry" : "run-layout-single"}`}>
          <LiveOutputPanel
            title={isRunning ? "Running..." : "hp.x Output"}
            output={output}
            placeholder="Starting calculation..."
            totalLineCount={outputLineCount}
            visibleLineCount={visibleOutputLineCount}
          />

          {isHpcMode && (
            <RemoteUtilizationPanel
              enabled={isRunning || activeTask?.status === "running"}
              profileId={activeHpcProfile?.id ?? null}
              remoteJobId={activeTask?.hpc?.remote_job_id ?? null}
              remoteNode={activeTask?.hpc?.remote_node ?? null}
              resourceType={activeTask?.hpc?.hpc_resource_type ?? hpcResources.resource_type}
            />
          )}
        </div>
      </div>
    );
  }

  function renderResultsStep() {
    const uValues = normalizeHubbardLrtUValues(result?.u_values);
    const artifactCount = Array.isArray(result?.artifacts) ? result.artifacts.length : 0;
    return (
      <div className="wizard-step results-step hubbard-lrt-results-step">
        <div className="results-hero">
          <div>
            <h3>Hubbard LRT Complete</h3>
            <p className="step-description">
              {isSaved ? "The calculation has been saved to this project." : "The calculation completed and is ready to save."}
            </p>
          </div>
          <span className={`calc-tag ${result?.converged ? "calc-tag-special" : "calc-tag-muted"}`}>
            {result?.converged ? "Converged" : "Completed"}
          </span>
        </div>
        {error && <div className="error-message">{error}</div>}

        <section className="lrt-results-card lrt-u-card">
          <div className="lrt-results-card-header">
            <h4>Calculated Hubbard U</h4>
            <span>{uValues.length} value{uValues.length === 1 ? "" : "s"}</span>
          </div>
          {uValues.length ? (
            <div className="lrt-u-grid">
              {uValues.map((entry) => (
                <div key={entry.target} className="lrt-u-value">
                  <span className="lrt-u-target">{entry.target}</span>
                  <strong>{entry.value_ev.toFixed(3)} eV</strong>
                </div>
              ))}
            </div>
          ) : (
            <p className="param-hint">No Hubbard U values were parsed from the hp.x artifacts.</p>
          )}
        </section>

        <section className="lrt-results-card">
          <div className="lrt-results-card-header">
            <h4>Run Summary</h4>
          </div>
          <div className="details-grid lrt-summary-grid">
            <div className="detail-item">
              <label>Q Mesh</label>
              <span>{qMesh.join(" × ")}</span>
            </div>
            <div className="detail-item">
              <label>Source SCF</label>
              <span>{formatCalculationSourceLabel(selectedScf)}</span>
            </div>
            <div className="detail-item">
              <label>Artifacts</label>
              <span>{artifactCount || "N/A"}</span>
            </div>
          </div>
        </section>
        <div className="step-actions">
          <button className="primary-button" onClick={onBack}>Back to Dashboard</button>
        </div>
      </div>
    );
  }

  return (
    <div className={`hubbard-lrt-wizard wizard-step-${step}`}>
      <div className="wizard-header">
        <button className="back-button" onClick={onBack}>
          ← Back
        </button>
        <h2>Hubbard LRT Wizard</h2>
        <div className="step-indicator">
          <span className={step === "source" ? "active" : "completed"}>
            1. Source
          </span>
          <span className={step === "parameters" ? "active" : ["run", "results"].includes(step) ? "completed" : ""}>
            2. Parameters
          </span>
          <span className={step === "run" ? "active" : step === "results" ? "completed" : ""}>
            3. Run
          </span>
          <span className={step === "results" ? "active" : ""}>
            4. Results
          </span>
        </div>
      </div>

      <div className="wizard-content">
        {step === "source" && renderSourceStep()}
        {step === "parameters" && renderParametersStep()}
        {step === "run" && renderRunStep()}
        {step === "results" && renderResultsStep()}
      </div>
    </div>
  );
}
