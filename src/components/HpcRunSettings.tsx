import { useEffect, useMemo, useRef, useState } from "react";
import { applyQeAutoPlan, planQeAutoAllocation } from "../lib/qeAutoParallelization";
import { normalizeCliDashText, previewSlurmScript } from "../lib/hpcConfig";
import { HpcResourceMode, SlurmResourceRequest } from "../lib/types";
import { InfoTooltip } from "./InfoTooltip";
import type { QePointEstimates } from "../lib/qePointEstimate";
import { QePointParallelization } from "./QePointParallelization";

interface HpcRunSettingsProps {
  profileId: string | null;
  profileName?: string | null;
  taskKind: string;
  pointEstimates?: QePointEstimates;
  commandLines: string[];
  resources: SlurmResourceRequest;
  resourceMode?: HpcResourceMode | null;
  defaultCpuResources?: SlurmResourceRequest | null;
  defaultGpuResources?: SlurmResourceRequest | null;
  resourceModeMessage?: string | null;
  lockResourceTypeWhenSingleMode?: boolean;
  maxTasks?: {
    value: number;
    reason: string;
  } | null;
  onResourcesChange: (next: SlurmResourceRequest) => void;
  disabled?: boolean;
  nodesDisabled?: boolean;
  tasksDisabled?: boolean;
}

interface ScriptPreviewState {
  script: string;
  sbatch: string;
  warnings: string[];
  errors: string[];
}

/** Keep intermediate keystrokes separate from the validated resource request. */
function HpcIntegerInput({ value, onCommit, min = 1, max = 4294967295, optional = false,
  placeholder, disabled = false }: {
  value: number | null; onCommit: (value: number | null) => void; min?: number; max?: number;
  optional?: boolean; placeholder?: string; disabled?: boolean;
}) {
  const [draft, setDraft] = useState(value == null ? "" : String(value));
  const editing = useRef(false);
  useEffect(() => {
    if (!editing.current) setDraft(value == null ? "" : String(value));
  }, [value]);
  const commit = () => {
    const parsed = Number(draft);
    const next = optional && draft.trim() === "" ? null
      : Math.min(max, Math.max(min, draft.trim() !== "" && Number.isFinite(parsed)
        ? Math.trunc(parsed) : value ?? min));
    setDraft(next == null ? "" : String(next));
    onCommit(next);
  };
  return <input type="number" min={min} max={max} step={1} required={!optional}
    value={draft} placeholder={placeholder} disabled={disabled}
    onFocus={() => { editing.current = true; }}
    onChange={(event) => setDraft(event.target.value)}
    onBlur={() => { editing.current = false; commit(); }}
    onKeyDown={(event) => { if (event.key === "Enter") { event.preventDefault(); commit(); } }} />;
}

function cloneResourceRequest(source: SlurmResourceRequest, resourceType: "cpu" | "gpu"): SlurmResourceRequest {
  return {
    ...source,
    resource_type: resourceType,
    additional_sbatch: [...(source.additional_sbatch || [])],
  };
}

export function HpcRunSettings({
  profileId,
  profileName,
  taskKind,
  pointEstimates,
  commandLines,
  resources,
  resourceMode = "both",
  defaultCpuResources = null,
  defaultGpuResources = null,
  resourceModeMessage = null,
  lockResourceTypeWhenSingleMode = true,
  maxTasks = null,
  onResourcesChange,
  disabled = false,
  nodesDisabled = false,
  tasksDisabled = false,
}: HpcRunSettingsProps) {
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [previewState, setPreviewState] = useState<ScriptPreviewState>({
    script: "",
    sbatch: "",
    warnings: [],
    errors: [],
  });
  const [previewError, setPreviewError] = useState<string | null>(null);
  const [loadingPreview, setLoadingPreview] = useState(false);
  const allowCpu = resourceMode !== "gpu_only";
  const allowGpu = resourceMode !== "cpu_only";
  const resourceModeLabel =
    resourceMode === "cpu_only"
      ? "CPU-only profile"
      : resourceMode === "gpu_only"
        ? "GPU-only profile"
        : null;
  const gpuCount = Math.max(1, resources.gpus || 1);
  const taskCount = Math.max(1, resources.ntasks || 1);
  const maxTaskCount = maxTasks && Number.isFinite(maxTasks.value) && maxTasks.value > 0
    ? Math.max(1, Math.floor(maxTasks.value))
    : null;
  const gpuOversubscriptionWarning = resources.resource_type === "gpu" && taskCount > gpuCount
    ? `Tasks (${taskCount}) exceed GPUs (${gpuCount}). This can oversubscribe GPUs and hurt QE performance.`
    : null;
  const resourceTypeSelectDisabled = disabled
    || (!allowCpu && !allowGpu)
    || (lockResourceTypeWhenSingleMode && (!allowCpu || !allowGpu));

  const commandLinesKey = useMemo(() => commandLines.join("\n"), [commandLines]);
  const pointParallelizationAvailable = ["scf", "bands", "dos", "phonon"].includes(taskKind) && /(?:pw|ph)\.x/.test(commandLinesKey);
  const autoAvailable = pointParallelizationAvailable && resources.resource_type === "cpu" && !nodesDisabled && !tasksDisabled;
  const defaultParallelMode = (resources.qe_kpoint_pools ?? 1) > 1 || (resources.qe_qpoint_images ?? 1) > 1 || (resources.cpus_per_task ?? 1) > 1 ? "manual" : "auto";
  const autoActive = autoAvailable && (resources.qe_parallel_mode ?? defaultParallelMode) === "auto";
  const totalCpus = resources.qe_total_cpus ?? taskCount * Math.max(1, resources.cpus_per_task ?? 1);
  const autoState = useMemo(() => {
    if (!autoActive) return { plan: null, error: null };
    try {
      return { plan: planQeAutoAllocation({
        nodes: resources.nodes,
        totalCpus,
        phonon: taskKind === "phonon",
        kPoints: pointEstimates?.k.conservative ? null : pointEstimates?.k.result?.irreduciblePoints,
        qPoints: pointEstimates?.q?.conservative ? null : pointEstimates?.qPointsToCompute ?? pointEstimates?.q?.result?.irreduciblePoints,
        maxRanks: maxTaskCount,
        sourceLayout: pointEstimates?.sourceLayout,
        irreps: resources.qe_phonon_irreps,
        fftPlanes: pointEstimates?.fftPlanes,
        bands: pointEstimates?.bands,
        pencilDecomposition: resources.pencil_decomposition === true,
      }), error: null };
    } catch (error) {
      return { plan: null, error: String(error) };
    }
  }, [autoActive, resources.nodes, totalCpus, taskKind, maxTaskCount,
    pointEstimates?.k.conservative, pointEstimates?.k.result?.irreduciblePoints,
    pointEstimates?.q?.conservative, pointEstimates?.q?.result?.irreduciblePoints,
    pointEstimates?.sourceLayout?.ranks, pointEstimates?.sourceLayout?.pools,
    resources.qe_phonon_irreps, pointEstimates?.fftPlanes, pointEstimates?.bands, pointEstimates?.qPointsToCompute, resources.pencil_decomposition]);

  useEffect(() => {
    if (disabled || !autoState.plan) return;
    const next = applyQeAutoPlan(resources, autoState.plan);
    if (next.ntasks !== resources.ntasks || next.cpus_per_task !== resources.cpus_per_task
      || next.ntasks_per_node !== resources.ntasks_per_node || next.qe_kpoint_pools !== resources.qe_kpoint_pools
      || next.qe_qpoint_images !== resources.qe_qpoint_images || next.qe_parallel_mode !== resources.qe_parallel_mode
      || next.qe_total_cpus !== resources.qe_total_cpus || next.pencil_decomposition !== resources.pencil_decomposition) {
      onResourcesChange(next);
    }
  }, [autoState.plan, disabled, resources, onResourcesChange]);

  const pencilDecompositionAvailable = /(^|\s)-pd\s+\.true\.(?=\s|$)/im.test(commandLinesKey);
  const previewCommandLines = useMemo(
    () => resources.pencil_decomposition !== true
      ? commandLines.map((line) => line.replace(/\s+-pd\s+\.true\./gi, ""))
      : commandLines,
    [commandLines, commandLinesKey, resources.pencil_decomposition, autoActive],
  );

  useEffect(() => {
    if (allowCpu && allowGpu) {
      return;
    }
    if (!allowCpu && resources.resource_type !== "gpu") {
      const next: SlurmResourceRequest = defaultGpuResources
        ? cloneResourceRequest(defaultGpuResources, "gpu")
        : { ...resources, resource_type: "gpu", gpus: Math.max(1, resources.gpus || 1) };
      onResourcesChange(next);
      return;
    }
    if (!allowGpu && resources.resource_type !== "cpu") {
      const next: SlurmResourceRequest = defaultCpuResources
        ? cloneResourceRequest(defaultCpuResources, "cpu")
        : { ...resources, resource_type: "cpu", gpus: 0 };
      onResourcesChange(next);
    }
  }, [
    allowCpu,
    allowGpu,
    resources,
    defaultCpuResources,
    defaultGpuResources,
    onResourcesChange,
  ]);

  useEffect(() => {
    if (maxTaskCount == null || taskCount <= maxTaskCount) {
      return;
    }
    onResourcesChange({
      ...resources,
      ntasks: maxTaskCount,
    });
  }, [maxTaskCount, onResourcesChange, resources, taskCount]);

  useEffect(() => {
    let cancelled = false;

    async function loadPreview() {
      if (!profileId) {
        setPreviewState({ script: "", sbatch: "", warnings: [], errors: [] });
        setPreviewError("Select an active HPC profile in Settings before submitting.");
        return;
      }

      setLoadingPreview(true);
      setPreviewError(null);
      try {
        const previewResources = autoState.plan ? applyQeAutoPlan(resources, autoState.plan) : resources;
        const result = await previewSlurmScript(taskKind, previewCommandLines, previewResources, profileId);
        if (cancelled) return;
        setPreviewState({
          script: result.script,
          sbatch: result.sbatch_preview,
          warnings: result.validation.warnings,
          errors: result.validation.errors,
        });
      } catch (error) {
        if (cancelled) return;
        setPreviewError(String(error));
        setPreviewState((prev) => ({
          ...prev,
          script: "",
          sbatch: "",
        }));
      } finally {
        if (!cancelled) {
          setLoadingPreview(false);
        }
      }
    }

    void loadPreview();
    return () => {
      cancelled = true;
    };
  }, [profileId, taskKind, resources, previewCommandLines, commandLinesKey, autoState.plan]);

  const update = (patch: Partial<SlurmResourceRequest>) => {
    onResourcesChange({
      ...resources,
      ...patch,
    });
  };
  const additionalSbatchText = (resources.additional_sbatch || []).join("\n");

  return (
    <section className="hpc-run-settings">
      <div className="hpc-run-settings-header">
        <h4>HPC Run Settings</h4>
        <span className="hpc-run-settings-cluster">
          Target: {profileName || "Andromeda"}
        </span>
      </div>

      <p className="hpc-run-settings-hint">
        Andromeda limits: interactive (12h), short (12h), medium (48h), long (5d). GPU jobs require
        <code> --gres=gpu:N </code>.
      </p>
      {resourceModeLabel && (
        <p className="hpc-run-settings-hint">
          {resourceModeMessage || `${resourceModeLabel}. Change this in HPC profile settings if needed.`}
        </p>
      )}
      {gpuOversubscriptionWarning && (
        <div className="hpc-script-warnings">
          <p>{gpuOversubscriptionWarning}</p>
        </div>
      )}
      {maxTaskCount != null && (
        <p className="hpc-run-settings-hint">
          Tasks are capped at {maxTaskCount} because {maxTasks?.reason}.
        </p>
      )}

      <div className="hpc-grid">
        <label>
          Resource Type
          <select
            value={resources.resource_type}
            onChange={(event) => {
              const nextType = event.target.value === "gpu" ? "gpu" : "cpu";
              const template =
                nextType === "gpu"
                  ? (defaultGpuResources ? cloneResourceRequest(defaultGpuResources, "gpu") : null)
                  : (defaultCpuResources ? cloneResourceRequest(defaultCpuResources, "cpu") : null);
              const next: SlurmResourceRequest = template || {
                ...resources,
                resource_type: nextType,
              };
              if (nextType === "gpu") {
                next.gpus = Math.max(1, next.gpus || 1);
                next.nodes = Math.max(1, next.nodes || 1);
              } else {
                next.gpus = 0;
              }
              next.qe_parallel_mode = resources.qe_parallel_mode;
              next.qe_total_cpus = (next.ntasks ?? 1) * (next.cpus_per_task ?? 1);
              next.qe_kpoint_pools = resources.qe_kpoint_pools;
              next.qe_qpoint_images = resources.qe_qpoint_images;
              next.qe_phonon_irreps = resources.qe_phonon_irreps;
              next.pencil_decomposition = resources.pencil_decomposition ?? false;
              next.ntasks_per_node = resources.ntasks_per_node ?? null;
              onResourcesChange(next);
            }}
            disabled={resourceTypeSelectDisabled}
          >
            <option value="cpu" disabled={!allowCpu}>CPU</option>
            <option value="gpu" disabled={!allowGpu}>GPU</option>
          </select>
        </label>

        <label>
          Partition
          <input
            value={resources.partition || ""}
            onChange={(event) => update({ partition: event.target.value })}
            placeholder="short"
            disabled={disabled}
          />
        </label>

        <label>
          Walltime (HH:MM:SS)
          <input
            value={resources.walltime || ""}
            onChange={(event) => update({ walltime: event.target.value })}
            placeholder="02:00:00"
            disabled={disabled}
          />
        </label>

        {autoAvailable && <label>
          <span className="hpc-checkbox-text">Parallelization <InfoTooltip text="Auto compares layouts using one, two and four OpenMP threads per MPI rank, favoring groups within a node and balanced point work. Saved SCF FFT dimensions guide group sizes when available. Known irreps can enable more phonon images than q-points. Source restart compatibility takes priority. This benchmark-informed heuristic does not measure performance or predict memory; use Manual for cluster-specific tuning." /></span>
          <select value={autoActive ? "auto" : "manual"} disabled={disabled}
            onChange={(event) => update({ qe_parallel_mode: event.target.value === "auto" ? "auto" : "manual", qe_total_cpus: taskCount * Math.max(1, resources.cpus_per_task ?? 1) })}>
            <option value="auto">Auto</option>
            <option value="manual">Manual</option>
          </select>
        </label>}

        <label>
          <span className="hpc-checkbox-text">Nodes <InfoTooltip text="Optional Slurm node count. Leave blank to let Slurm choose enough nodes for the requested CPUs. Auto cannot optimize for node boundaries without an explicit count." /></span>
          <HpcIntegerInput value={resources.nodes ?? null} optional placeholder="Auto"
            max={autoActive ? (maxTaskCount ?? 4294967295) : 4294967295}
            onCommit={(nodes) => update({ nodes })} disabled={disabled || nodesDisabled} />
        </label>

        {autoActive ? <label>
          <span className="hpc-checkbox-text">Total allocated CPUs <InfoTooltip text="Total CPU cores across all requested nodes, not CPUs per node. Auto balances MPI ranks with OpenMP threads based on allocation size and available point parallelism. Slurm receives the displayed ranks × CPUs per rank, preserving your full CPU allocation." /></span>
          <HpcIntegerInput value={totalCpus} min={resources.nodes ?? 1}
            disabled={disabled} onCommit={(value) => update({ qe_total_cpus: value! })} />
        </label> : <>
        <label>
          Tasks
          <HpcIntegerInput value={resources.ntasks ?? 1} max={maxTaskCount ?? 4294967295}
            onCommit={(value) => update({ ntasks: value! })} disabled={disabled || tasksDisabled} />
        </label>

        <label>
          CPUs / Task
          <HpcIntegerInput value={resources.cpus_per_task ?? 1}
            onCommit={(value) => update({ cpus_per_task: value! })} disabled={disabled} />
        </label>
        </>}

        <label>
          Memory per node (GiB)
          <HpcIntegerInput value={resources.memory_gb ?? 16}
            onCommit={(value) => update({ memory_gb: value! })} disabled={disabled} />
        </label>

        {resources.resource_type === "gpu" && (
          <label>
            GPUs
            <HpcIntegerInput value={resources.gpus ?? 1}
            onCommit={(value) => update({ gpus: value! })} disabled={disabled} />
          </label>
        )}
      </div>

      {autoState.error && <p className="hpc-script-error">{autoState.error}</p>}
      {pointParallelizationAvailable && (
        <QePointParallelization ranks={autoState.plan?.ranks ?? taskCount} pools={autoState.plan?.pools ?? resources.qe_kpoint_pools ?? 1}
          estimates={pointEstimates} images={autoState.plan?.images ?? resources.qe_qpoint_images ?? 1}
          autoPlan={autoActive ? autoState.plan ?? undefined : undefined}
          phonon={taskKind === "phonon"} disabled={disabled || autoActive}
          onChange={(pools, images) => update({ qe_parallel_mode: "manual", qe_kpoint_pools: pools, qe_qpoint_images: images })} />
      )}

      <button
        type="button"
        className="hpc-advanced-toggle"
        onClick={() => setShowAdvanced((prev) => !prev)}
        aria-expanded={showAdvanced}
        aria-controls="hpc-advanced-slurm-fields"
      >
        {showAdvanced
          ? "Hide Advanced Slurm Fields"
          : "Show Advanced Slurm Fields"}
      </button>

      {showAdvanced && (
        <div className="hpc-advanced-grid" id="hpc-advanced-slurm-fields">
          {taskKind === "phonon" && autoAvailable && <label>
            <span className="hpc-checkbox-text">Known phonon irreps <InfoTooltip text="Optional total number of irreps to compute across the selected q-points, from QE initialization for this exact task. Do not enter the number of displacement modes or the irrep count at just one q-point. Auto may use more images than q-points, but never more than this count. Leave blank when unknown; update it when changing the q-grid or representation range." /></span>
            <HpcIntegerInput value={resources.qe_phonon_irreps ?? null} optional placeholder="Optional"
              disabled={disabled} onCommit={(value) => update({ qe_phonon_irreps: value })} />
          </label>}
          <label>
            <span className="hpc-checkbox-text">
              Tasks / Node
              <InfoTooltip text="Slurm `--ntasks-per-node`. Limits how many MPI tasks are placed on each node; leave blank to let Slurm distribute the total task count." />
            </span>
            <HpcIntegerInput value={resources.ntasks_per_node ?? null} optional placeholder="Auto"
              disabled={disabled || autoActive} onCommit={(value) => update({ ntasks_per_node: value })} />
          </label>
          <label className="hpc-checkbox-label">
            <input
              type="checkbox"
              checked={resources.contiguous ?? false}
              onChange={(event) => update({ contiguous: event.target.checked })}
              disabled={disabled}
            />
            <span className="hpc-checkbox-text">
              Contiguous nodes
              <InfoTooltip text="Adds Slurm `--contiguous` to request a contiguous set of compute nodes." />
            </span>
          </label>
          {pencilDecompositionAvailable && (
            <label className="hpc-checkbox-label">
              <input
                type="checkbox"
                checked={resources.pencil_decomposition ?? false}
                onChange={(event) => update({ pencil_decomposition: event.target.checked })}
                disabled={disabled}
              />
              <span className="hpc-checkbox-text">
                Pencil decomposition
                <InfoTooltip text="Adds Quantum ESPRESSO's `-pd .true.` option to supported MPI commands. Off by default for reliability across multiple nodes. You can enable it in Auto or Manual when your QE build and cluster MPI setup support it reliably." />
              </span>
            </label>
          )}
          <label>
            QoS
            <input
              value={resources.qos || ""}
              onChange={(event) => update({ qos: event.target.value || null })}
              disabled={disabled}
            />
          </label>
          <label>
            Account
            <input
              value={resources.account || ""}
              onChange={(event) => update({ account: event.target.value || null })}
              disabled={disabled}
            />
          </label>
          <label>
            Constraint
            <input
              value={resources.constraint || ""}
              onChange={(event) => update({ constraint: event.target.value || null })}
              disabled={disabled}
            />
          </label>

          <label className="hpc-advanced-wide">
            Module / Preamble
            <textarea
              rows={3}
              value={resources.module_preamble || ""}
              onChange={(event) =>
                update({ module_preamble: normalizeCliDashText(event.target.value) || null })
              }
              placeholder={"module purge\nmodule load qe"}
              disabled={disabled}
              autoCorrect="off"
              autoCapitalize="off"
              spellCheck={false}
            />
          </label>

          <label className="hpc-advanced-wide">
            Extra SBATCH Lines
            <textarea
              rows={3}
              value={additionalSbatchText}
              onChange={(event) => {
                const parsed = normalizeCliDashText(event.target.value)
                  .split(/\r?\n/)
                  .map((line) => line.trim())
                  .filter((line) => line.length > 0);
                update({ additional_sbatch: parsed });
              }}
              placeholder={"--mail-type=END\n--mail-user=you@example.edu"}
              disabled={disabled}
              autoCorrect="off"
              autoCapitalize="off"
              spellCheck={false}
            />
          </label>
        </div>
      )}

      <div className="hpc-script-preview">
        <h5>Generated Slurm Script (read-only)</h5>
        {loadingPreview && <p className="hpc-script-loading">Refreshing preview...</p>}
        {previewError && <p className="hpc-script-error">{previewError}</p>}
        {!previewError && previewState.sbatch && (
          <p className="hpc-script-cmd">
            <strong>Submit command:</strong> <code>{previewState.sbatch}</code>
          </p>
        )}
        {previewState.errors.length > 0 && (
          <div className="hpc-script-errors">
            {previewState.errors.map((msg) => (
              <p key={msg}>{msg}</p>
            ))}
          </div>
        )}
        {previewState.warnings.length > 0 && (
          <div className="hpc-script-warnings">
            {previewState.warnings.map((msg) => (
              <p key={msg}>{msg}</p>
            ))}
          </div>
        )}
        <pre className="hpc-script-text">{previewState.script || "Script preview unavailable."}</pre>
      </div>
    </section>
  );
}
