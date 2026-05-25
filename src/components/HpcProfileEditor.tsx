import { useEffect, useMemo, useState, type FormEvent } from "react";
import { addEngineInstallation, type EngineInstallation } from "../lib/engines";
import { normalizeCliDashText, saveHpcProfile } from "../lib/hpcConfig";
import type { EngineId } from "../lib/engines/types";
import type { HpcProfile, SlurmResourceRequest } from "../lib/types";

interface HpcProfileEditorProps {
  profile: HpcProfile;
  installations: EngineInstallation[];
  onSaved: (profile: HpcProfile, message: string) => void;
}

function normalizedOptional(value: string | null | undefined): string | null {
  const trimmed = String(value ?? "").trim();
  return trimmed.length > 0 ? trimmed : null;
}

function wien2kRootForProfile(profile: HpcProfile, installations: EngineInstallation[]): string | null {
  return normalizedOptional(profile.remote_wien2k_install_root)
    ?? installations.find(
      (installation) => installation.engineId === "wien2k" && installation.hpcProfileId === profile.id,
    )?.remoteInstallRoot
    ?? null;
}

function withInstallationFallback(profile: HpcProfile, installations: EngineInstallation[]): HpcProfile {
  return {
    ...profile,
    remote_wien2k_install_root: wien2kRootForProfile(profile, installations),
  };
}

function positiveInt(raw: string, fallback: number): number {
  const parsed = Number.parseInt(raw, 10);
  return Number.isFinite(parsed) && parsed >= 1 ? parsed : fallback;
}

function optionalText(raw: string): string | null {
  return raw.trim().length > 0 ? raw : null;
}

function deriveRemotePostw90Path(remoteWannier90Path: string | null | undefined): string {
  const remoteWannier90 = String(remoteWannier90Path ?? "").trim();
  if (remoteWannier90.includes("/") || remoteWannier90.startsWith("~")) {
    const segments = remoteWannier90.split("/");
    segments[segments.length - 1] = "postw90.x";
    return segments.join("/").trim() || "postw90.x";
  }
  return "postw90.x";
}

function resourceText(resources: SlurmResourceRequest, key: "module_preamble"): string {
  return resources[key] ?? "";
}

export function HpcProfileEditor({ profile, installations, onSaved }: HpcProfileEditorProps) {
  const [draft, setDraft] = useState<HpcProfile>(() => withInstallationFallback(profile, installations));
  const [credential, setCredential] = useState("");
  const [persistCredential, setPersistCredential] = useState(profile.credential_persisted);
  const [selectedEngine, setSelectedEngine] = useState<EngineId>("qe");
  const [isSaving, setIsSaving] = useState(false);
  const [status, setStatus] = useState<string | null>(null);
  const savedWien2kRoot = useMemo(
    () => wien2kRootForProfile(profile, installations),
    [installations, profile],
  );
  const hasWien2k = installations.some(
    (installation) => installation.engineId === "wien2k" && installation.hpcProfileId === profile.id,
  );

  useEffect(() => {
    setDraft(withInstallationFallback(profile, installations));
    setCredential("");
    setPersistCredential(profile.credential_persisted);
  }, [installations, profile]);

  function updateResource(
    resourceType: "cpu" | "gpu",
    patch: Partial<SlurmResourceRequest>,
  ) {
    setDraft((current) => ({
      ...current,
      [resourceType === "cpu" ? "default_cpu_resources" : "default_gpu_resources"]: {
        ...(resourceType === "cpu" ? current.default_cpu_resources : current.default_gpu_resources),
        ...patch,
      },
    }));
  }

  async function handleSave(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (isSaving) return;

    const desiredWien2kRoot = normalizedOptional(draft.remote_wien2k_install_root);
    if (hasWien2k && !desiredWien2kRoot) {
      setStatus("A configured WIEN2k engine must retain a remote WIENROOT.");
      return;
    }

    setIsSaving(true);
    setStatus(null);
    let savedBase: HpcProfile | null = null;
    try {
      // Keep the last verified WIEN2k root until a changed root has been
      // successfully checked on the remote host.
      savedBase = await saveHpcProfile(
        { ...draft, remote_wien2k_install_root: savedWien2kRoot },
        credential.trim().length > 0 ? credential : null,
        persistCredential,
      );

      if (hasWien2k && desiredWien2kRoot !== savedWien2kRoot) {
        setStatus("Verifying updated WIEN2k installation...");
        await addEngineInstallation({
          engineId: "wien2k",
          hpcProfileId: savedBase.id,
          remoteInstallRoot: desiredWien2kRoot!,
        });
        const verified = { ...savedBase, remote_wien2k_install_root: desiredWien2kRoot };
        setDraft(verified);
        onSaved(verified, "Profile saved and WIEN2k installation verified.");
        setStatus("Saved. Updated WIEN2k installation verified.");
        return;
      }

      setDraft(savedBase);
      onSaved(savedBase, "HPC profile saved.");
      setStatus("Saved.");
    } catch (reason) {
      if (savedBase) {
        onSaved(savedBase, "Profile settings saved, but engine verification failed.");
      }
      setStatus(`Failed to save profile: ${reason}`);
    } finally {
      setIsSaving(false);
    }
  }

  function renderResourceCard(resourceType: "cpu" | "gpu") {
    const resources = resourceType === "cpu" ? draft.default_cpu_resources : draft.default_gpu_resources;
    const label = resourceType.toUpperCase();
    return (
      <div className="settings-hpc-default-card">
        <h4>{label} Defaults</h4>
        <div className="hpc-grid">
          <label>
            Partition
            <input value={resources.partition ?? ""} onChange={(event) => updateResource(resourceType, { partition: optionalText(event.target.value) })} />
          </label>
          <label>
            Walltime
            <input value={resources.walltime ?? ""} onChange={(event) => updateResource(resourceType, { walltime: optionalText(event.target.value) })} placeholder="02:00:00" />
          </label>
          <label>
            Nodes
            <input type="number" min={1} value={resources.nodes ?? 1} onChange={(event) => updateResource(resourceType, { nodes: positiveInt(event.target.value, 1) })} />
          </label>
          <label>
            Tasks
            <input type="number" min={1} value={resources.ntasks ?? 1} onChange={(event) => updateResource(resourceType, { ntasks: positiveInt(event.target.value, 1) })} />
          </label>
          <label>
            CPUs / Task
            <input type="number" min={1} value={resources.cpus_per_task ?? 1} onChange={(event) => updateResource(resourceType, { cpus_per_task: positiveInt(event.target.value, 1) })} />
          </label>
          <label>
            Memory (GB)
            <input type="number" min={1} value={resources.memory_gb ?? 16} onChange={(event) => updateResource(resourceType, { memory_gb: positiveInt(event.target.value, 16) })} />
          </label>
          {resourceType === "gpu" && (
            <label>
              GPUs
              <input type="number" min={1} value={resources.gpus ?? 1} onChange={(event) => updateResource(resourceType, { gpus: positiveInt(event.target.value, 1) })} />
            </label>
          )}
        </div>
        <div className="hpc-advanced-grid">
          <label>
            QoS
            <input value={resources.qos ?? ""} onChange={(event) => updateResource(resourceType, { qos: optionalText(event.target.value) })} />
          </label>
          <label>
            Account
            <input value={resources.account ?? ""} onChange={(event) => updateResource(resourceType, { account: optionalText(event.target.value) })} />
          </label>
          <label>
            Constraint
            <input value={resources.constraint ?? ""} onChange={(event) => updateResource(resourceType, { constraint: optionalText(event.target.value) })} />
          </label>
          <label className="hpc-advanced-wide">
            Module / Preamble
            <textarea
              rows={3}
              value={resourceText(resources, "module_preamble")}
              onChange={(event) => updateResource(resourceType, { module_preamble: optionalText(normalizeCliDashText(event.target.value)) })}
            />
          </label>
          <label className="hpc-advanced-wide">
            Additional SBATCH Directives
            <textarea
              rows={3}
              value={(resources.additional_sbatch ?? []).join("\n")}
              onChange={(event) => updateResource(resourceType, {
                additional_sbatch: normalizeCliDashText(event.target.value)
                  .split("\n")
                  .map((line) => line.trim())
                  .filter(Boolean),
              })}
            />
          </label>
        </div>
      </div>
    );
  }

  return (
    <form id="hpc-profile-editor-form" className="settings-pane hpc-profile-editor" onSubmit={(event) => void handleSave(event)}>
      <div className="settings-menu-section">
        <label className="settings-menu-label">Connection</label>
        <div className="hpc-profile-editor-grid">
          <label>
            Profile Name
            <input className="settings-menu-input" value={draft.name} onChange={(event) => setDraft((current) => ({ ...current, name: event.target.value }))} />
          </label>
          <label>
            Cluster
            <input className="settings-menu-input" value={draft.cluster} onChange={(event) => setDraft((current) => ({ ...current, cluster: event.target.value }))} />
          </label>
          <label>
            SSH Host
            <input className="settings-menu-input" value={draft.host} onChange={(event) => setDraft((current) => ({ ...current, host: event.target.value }))} />
          </label>
          <label>
            SSH Port
            <input className="settings-menu-input" type="number" min={1} max={65535} value={draft.port} onChange={(event) => setDraft((current) => ({ ...current, port: positiveInt(event.target.value, 22) }))} />
          </label>
          <label>
            Username
            <input className="settings-menu-input" value={draft.username} onChange={(event) => setDraft((current) => ({ ...current, username: event.target.value }))} />
          </label>
          <label>
            Authentication Method
            <select className="settings-menu-input" value={draft.auth_method} onChange={(event) => setDraft((current) => ({ ...current, auth_method: event.target.value === "password" ? "password" : "ssh_key" }))}>
              <option value="ssh_key">SSH Key</option>
              <option value="password">Password</option>
            </select>
          </label>
          {draft.auth_method === "ssh_key" && (
            <label>
              SSH Key Path
              <input className="settings-menu-input" value={draft.ssh_key_path ?? ""} onChange={(event) => setDraft((current) => ({ ...current, ssh_key_path: optionalText(event.target.value) }))} />
            </label>
          )}
          <label>
            {draft.auth_method === "password" ? "Password" : "Key Passphrase (optional)"}
            <input className="settings-menu-input" type="password" value={credential} onChange={(event) => setCredential(event.target.value)} placeholder="Unchanged if left blank" />
          </label>
        </div>
        <label className="toggle-label">
          <input type="checkbox" checked={persistCredential} onChange={(event) => setPersistCredential(event.target.checked)} />
          <span>Store credential in encrypted keychain</span>
        </label>
      </div>

      <div className="settings-menu-section">
        <label className="settings-menu-label">Remote Project Storage</label>
        <p className="settings-menu-hint">
          These roots are shared by all engines in this profile. Engine runs stage beneath the project and engine namespace.
        </p>
        <div className="hpc-profile-editor-grid">
          <label>
            Remote Workspace Root
            <input className="settings-menu-input" value={draft.remote_workspace_root} onChange={(event) => setDraft((current) => ({ ...current, remote_workspace_root: event.target.value }))} />
          </label>
          <label>
            Remote Project Root
            <input className="settings-menu-input" value={draft.remote_project_root} onChange={(event) => setDraft((current) => ({ ...current, remote_project_root: event.target.value }))} />
          </label>
        </div>
      </div>

      <div className="settings-menu-section settings-profile-engines">
        <div className="settings-engine-header">
          <label className="settings-menu-label">Engines</label>
          <div className="engine-switcher" role="group" aria-label="Profile computation engine settings">
            <button type="button" className={selectedEngine === "qe" ? "active" : ""} aria-pressed={selectedEngine === "qe"} onClick={() => setSelectedEngine("qe")}>QE</button>
            {hasWien2k && (
              <button type="button" className={selectedEngine === "wien2k" ? "active" : ""} aria-pressed={selectedEngine === "wien2k"} onClick={() => setSelectedEngine("wien2k")}>W2K</button>
            )}
          </div>
        </div>

        {selectedEngine === "qe" && (
          <>
            <p className="settings-menu-hint">Quantum ESPRESSO is the default engine. Pseudopotentials apply only to this engine.</p>
            <div className="hpc-profile-editor-grid">
              <label>
                CPU Bin Directory
                <input className="settings-menu-input" value={draft.remote_qe_cpu_bin_dir ?? draft.remote_qe_bin_dir} onChange={(event) => setDraft((current) => ({ ...current, remote_qe_bin_dir: event.target.value, remote_qe_cpu_bin_dir: optionalText(event.target.value) }))} />
              </label>
              <label>
                GPU Bin Directory
                <input className="settings-menu-input" value={draft.remote_qe_gpu_bin_dir ?? ""} onChange={(event) => setDraft((current) => ({ ...current, remote_qe_gpu_bin_dir: optionalText(event.target.value) }))} />
              </label>
              <label>
                CPU Pseudopotential Directory
                <input className="settings-menu-input" value={draft.remote_cpu_pseudo_dir ?? draft.remote_pseudo_dir} onChange={(event) => setDraft((current) => ({ ...current, remote_pseudo_dir: event.target.value, remote_cpu_pseudo_dir: optionalText(event.target.value) }))} />
              </label>
              <label>
                GPU Pseudopotential Directory
                <input className="settings-menu-input" value={draft.remote_gpu_pseudo_dir ?? ""} onChange={(event) => setDraft((current) => ({ ...current, remote_gpu_pseudo_dir: optionalText(event.target.value) }))} />
              </label>
              <label>
                EPW Executable Override
                <input className="settings-menu-input" value={draft.remote_epw_path ?? ""} onChange={(event) => setDraft((current) => ({ ...current, remote_epw_path: optionalText(event.target.value) }))} />
              </label>
              <label>
                Wannier90 Executable
                <input className="settings-menu-input" value={draft.remote_wannier90_path ?? ""} onChange={(event) => setDraft((current) => ({ ...current, remote_wannier90_path: optionalText(event.target.value) }))} />
              </label>
              <label>
                Postw90 Executable (derived)
                <input className="settings-menu-input" value={deriveRemotePostw90Path(draft.remote_wannier90_path)} readOnly />
              </label>
            </div>
          </>
        )}
        {selectedEngine === "wien2k" && (
          <>
            <p className="settings-menu-hint">
              WIEN2k uses a remote WIENROOT installation and does not use Quantum ESPRESSO pseudopotential paths. Changing this directory is verified when saved.
            </p>
            <label>
              Remote WIENROOT
              <input
                className="settings-menu-input"
                value={draft.remote_wien2k_install_root ?? ""}
                onChange={(event) => setDraft((current) => ({ ...current, remote_wien2k_install_root: optionalText(event.target.value) }))}
                placeholder="/opt/WIEN2k"
              />
            </label>
          </>
        )}
      </div>

      <div className="settings-menu-section settings-hpc-defaults-section">
        <label className="settings-menu-label">Scheduler Defaults</label>
        <div className="settings-hpc-launcher-grid">
          <label>
            Supported Resource Types
            <select value={draft.resource_mode} onChange={(event) => setDraft((current) => ({ ...current, resource_mode: event.target.value === "cpu_only" || event.target.value === "gpu_only" ? event.target.value : "both" }))}>
              <option value="both">CPU + GPU</option>
              <option value="cpu_only">CPU only</option>
              <option value="gpu_only">GPU only</option>
            </select>
          </label>
          <label>
            MPI Launcher
            <select value={draft.launcher} onChange={(event) => setDraft((current) => ({ ...current, launcher: event.target.value === "mpirun" ? "mpirun" : "srun" }))}>
              <option value="srun">srun</option>
              <option value="mpirun">mpirun</option>
            </select>
          </label>
          <label>
            CPU Launcher Extra Args
            <input value={draft.launcher_cpu_extra_args ?? draft.launcher_extra_args ?? ""} onChange={(event) => setDraft((current) => ({ ...current, launcher_cpu_extra_args: optionalText(normalizeCliDashText(event.target.value)), launcher_extra_args: null }))} />
          </label>
          <label>
            GPU Launcher Extra Args
            <input value={draft.launcher_gpu_extra_args ?? draft.launcher_extra_args ?? ""} onChange={(event) => setDraft((current) => ({ ...current, launcher_gpu_extra_args: optionalText(normalizeCliDashText(event.target.value)), launcher_extra_args: null }))} />
          </label>
        </div>
        <div className="settings-hpc-defaults-grid">
          {renderResourceCard("cpu")}
          {renderResourceCard("gpu")}
        </div>
      </div>
      {status && <div className="settings-menu-status hpc-profile-editor-status">{status}</div>}
      {isSaving && <p className="settings-menu-hint">Saving profile...</p>}
    </form>
  );
}
