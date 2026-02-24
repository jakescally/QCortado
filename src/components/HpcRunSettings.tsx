import { useEffect, useMemo, useState } from "react";
import { normalizeCliDashText, previewSlurmScript } from "../lib/hpcConfig";
import { HpcResourceMode, SlurmResourceRequest } from "../lib/types";

interface HpcRunSettingsProps {
  profileId: string | null;
  profileName?: string | null;
  taskKind: string;
  commandLines: string[];
  resources: SlurmResourceRequest;
  resourceMode?: HpcResourceMode | null;
  defaultCpuResources?: SlurmResourceRequest | null;
  defaultGpuResources?: SlurmResourceRequest | null;
  onResourcesChange: (next: SlurmResourceRequest) => void;
  disabled?: boolean;
}

interface ScriptPreviewState {
  script: string;
  sbatch: string;
  warnings: string[];
  errors: string[];
}

function normalizePositiveInt(input: string, fallback: number): number {
  const parsed = Number.parseInt(input, 10);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    return fallback;
  }
  return parsed;
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
  commandLines,
  resources,
  resourceMode = "both",
  defaultCpuResources = null,
  defaultGpuResources = null,
  onResourcesChange,
  disabled = false,
}: HpcRunSettingsProps) {
  const [showAdvanced, setShowAdvanced] = useState(true);
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
  const gpuOversubscriptionWarning = resources.resource_type === "gpu" && taskCount > gpuCount
    ? `Tasks (${taskCount}) exceed GPUs (${gpuCount}). This can oversubscribe GPUs and hurt QE performance.`
    : null;

  const commandLinesKey = useMemo(() => commandLines.join("\n"), [commandLines]);

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
        const result = await previewSlurmScript(taskKind, commandLines, resources, profileId);
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
  }, [profileId, taskKind, resources, commandLines, commandLinesKey]);

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
        Andromeda limits: interactive (12h), short (24h), medium (72h), long (7d). GPU jobs require
        <code> --gres=gpu:N </code>.
      </p>
      {resourceModeLabel && (
        <p className="hpc-run-settings-hint">
          {resourceModeLabel}. Change this in HPC profile settings if needed.
        </p>
      )}
      {gpuOversubscriptionWarning && (
        <div className="hpc-script-warnings">
          <p>{gpuOversubscriptionWarning}</p>
        </div>
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
              onResourcesChange(next);
            }}
            disabled={disabled || !allowCpu || !allowGpu}
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

        <label>
          Nodes
          <input
            type="number"
            min={1}
            value={resources.nodes ?? 1}
            onChange={(event) => update({ nodes: normalizePositiveInt(event.target.value, 1) })}
            disabled={disabled}
          />
        </label>

        <label>
          Tasks
          <input
            type="number"
            min={1}
            value={resources.ntasks ?? 1}
            onChange={(event) => update({ ntasks: normalizePositiveInt(event.target.value, 1) })}
            disabled={disabled}
          />
        </label>

        <label>
          CPUs / Task
          <input
            type="number"
            min={1}
            value={resources.cpus_per_task ?? 1}
            onChange={(event) =>
              update({ cpus_per_task: normalizePositiveInt(event.target.value, 1) })
            }
            disabled={disabled}
          />
        </label>

        <label>
          Memory (GB)
          <input
            type="number"
            min={1}
            value={resources.memory_gb ?? 16}
            onChange={(event) => update({ memory_gb: normalizePositiveInt(event.target.value, 16) })}
            disabled={disabled}
          />
        </label>

        {resources.resource_type === "gpu" && (
          <label>
            GPUs
            <input
              type="number"
              min={1}
              value={resources.gpus ?? 1}
              onChange={(event) => update({ gpus: normalizePositiveInt(event.target.value, 1) })}
              disabled={disabled}
            />
          </label>
        )}
      </div>

      <button
        type="button"
        className="hpc-advanced-toggle"
        onClick={() => setShowAdvanced((prev) => !prev)}
      >
        {showAdvanced
          ? "Hide Advanced Slurm Fields"
          : "Show Advanced Slurm Fields (qos/account/constraint/module preamble)"}
      </button>

      {showAdvanced && (
        <div className="hpc-advanced-grid">
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
