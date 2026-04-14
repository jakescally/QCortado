import { useEffect, useMemo, useRef, useState } from "react";
import { open } from "@tauri-apps/plugin-dialog";
import {
  defaultCpuResources,
  defaultGpuResources,
  getActiveHpcProfileId,
  importHpcPresetBundle,
  listHpcProfiles,
  normalizeCliDashText,
  saveHpcProfile,
  setActiveHpcProfile,
  testHpcConnection,
  validateHpcEnvironment,
} from "../lib/hpcConfig";
import { HpcProfile } from "../lib/types";

interface HpcSetupWizardProps {
  isOpen: boolean;
  initialProfile?: HpcProfile | null;
  onClose: () => void;
  onSaved: (profile: HpcProfile) => void;
}

const STEP_TITLES = ["Connection", "Authentication", "Paths", "Scheduler", "Validation"];

function deriveRemotePostw90Path(remoteWannier90Path: string | null | undefined): string {
  const remoteWannier90 = String(remoteWannier90Path || "").trim();
  if (remoteWannier90.includes("/") || remoteWannier90.startsWith("~")) {
    const segments = remoteWannier90.split("/");
    segments[segments.length - 1] = "postw90.x";
    const derived = segments.join("/");
    return derived.trim().length > 0 ? derived : "postw90.x";
  }
  return "postw90.x";
}

function makeDefaultProfile(): HpcProfile {
  return {
    id: "",
    name: "BC Andromeda",
    cluster: "andromeda",
    host: "andromeda.bc.edu",
    port: 22,
    username: "",
    auth_method: "ssh_key",
    ssh_key_path: "~/.ssh/id_rsa",
    remote_qe_bin_dir: "~/qe/bin",
    remote_qe_cpu_bin_dir: "~/qe/bin",
    remote_qe_gpu_bin_dir: "~/qe-gpu/bin",
    remote_wannier90_path: "wannier90.x",
    remote_pseudo_dir: "~/qe/pseudo",
    remote_workspace_root: "~/qcortado/work",
    remote_project_root: "~/qcortado/projects",
    resource_mode: "both",
    launcher: "srun",
    launcher_extra_args: null,
    default_cpu_resources: defaultCpuResources(),
    default_gpu_resources: defaultGpuResources(),
    credential_persisted: false,
    warning_acknowledged: false,
    created_at: "",
    updated_at: "",
  };
}

function normalizeQeProfilePaths(profile: HpcProfile): HpcProfile {
  const legacy = (profile.remote_qe_bin_dir || "").trim();
  const fallback = legacy.length > 0 ? legacy : "~/qe/bin";
  const cpu = (profile.remote_qe_cpu_bin_dir || "").trim();
  const gpu = (profile.remote_qe_gpu_bin_dir || "").trim();
  const primary = cpu.length > 0
    ? cpu
    : gpu.length > 0
      ? gpu
      : fallback;
  const normalizedCpu = cpu.length > 0 ? cpu : primary;
  const normalizedGpu = gpu.length > 0 ? gpu : primary;

  return {
    ...profile,
    remote_qe_bin_dir: primary,
    remote_qe_cpu_bin_dir: normalizedCpu,
    remote_qe_gpu_bin_dir: normalizedGpu,
  };
}

export function HpcSetupWizard({
  isOpen,
  initialProfile = null,
  onClose,
  onSaved,
}: HpcSetupWizardProps) {
  const wasOpenRef = useRef(false);
  const [step, setStep] = useState(0);
  const [profile, setProfile] = useState<HpcProfile>(normalizeQeProfilePaths(initialProfile ?? makeDefaultProfile()));
  const [credential, setCredential] = useState("");
  const [persistCredential, setPersistCredential] = useState(initialProfile?.credential_persisted ?? false);
  const [isSaving, setIsSaving] = useState(false);
  const [isImportingPresetBundle, setIsImportingPresetBundle] = useState(false);
  const [importMessage, setImportMessage] = useState<string | null>(null);
  const [validationMessage, setValidationMessage] = useState<string | null>(null);
  const [validationDetails, setValidationDetails] = useState<string[]>([]);
  const isBusy = isSaving || isImportingPresetBundle;

  useEffect(() => {
    const isNewlyOpened = isOpen && !wasOpenRef.current;
    if (isNewlyOpened) {
      setStep(0);
      setProfile(normalizeQeProfilePaths(initialProfile ?? makeDefaultProfile()));
      setCredential("");
      setPersistCredential(initialProfile?.credential_persisted ?? false);
      setImportMessage(null);
      setValidationMessage(null);
      setValidationDetails([]);
    }
    wasOpenRef.current = isOpen;
  }, [initialProfile, isOpen]);

  const canContinue = useMemo(() => {
    if (step === 0) {
      return profile.host.trim().length > 0 && profile.username.trim().length > 0 && profile.port > 0;
    }
    if (step === 1) {
      if (profile.auth_method === "ssh_key") {
        return (profile.ssh_key_path || "").trim().length > 0;
      }
      return credential.trim().length > 0 || profile.id.trim().length > 0;
    }
    if (step === 2) {
      const legacyPath = (profile.remote_qe_bin_dir || "").trim();
      const cpuPath = (profile.remote_qe_cpu_bin_dir || "").trim();
      const gpuPath = (profile.remote_qe_gpu_bin_dir || "").trim();
      const resolvedCpuPath = cpuPath.length > 0 ? cpuPath : legacyPath;
      const resolvedGpuPath = gpuPath.length > 0 ? gpuPath : legacyPath;
      const hasCpuPath = profile.resource_mode !== "gpu_only" ? resolvedCpuPath.length > 0 : true;
      const hasGpuPath = profile.resource_mode !== "cpu_only" ? resolvedGpuPath.length > 0 : true;
      return (
        hasCpuPath
        && hasGpuPath
        && profile.remote_pseudo_dir.trim().length > 0
        && profile.remote_workspace_root.trim().length > 0
        && profile.remote_project_root.trim().length > 0
      );
    }
    return true;
  }, [credential, profile, step]);

  async function handlePersistProfile(runValidation: boolean) {
    setIsSaving(true);
    setValidationMessage(null);
    try {
      const normalizedProfile = normalizeQeProfilePaths(profile);
      const saved = await saveHpcProfile(
        normalizedProfile,
        credential.trim().length > 0 ? credential : null,
        persistCredential,
      );
      await setActiveHpcProfile(saved.id);
      setProfile(saved);
      onSaved(saved);

      if (runValidation) {
        const connection = await testHpcConnection(saved.id);
        const validation = await validateHpcEnvironment(saved.id);
        const detail: string[] = [];
        detail.push(connection.success ? "SSH connection successful." : `SSH connection failed: ${connection.message}`);
        if (!validation.sbatch_available) detail.push("`sbatch` is not available.");
        if (!validation.squeue_available) detail.push("`squeue` is not available.");
        if (!validation.sacct_available) detail.push("`sacct` is not available.");
        if (!validation.qe_pw_available) detail.push("`pw.x` is not executable at configured QE path.");
        if (!validation.qe_epw_available) detail.push("`epw.x` is not executable at configured QE path.");
        if (!validation.qe_pw2wannier_available) detail.push("`pw2wannier90.x` is not executable at configured QE path.");
        if (!validation.wannier90_available) detail.push("`wannier90.x` is not executable at configured path or on `PATH`.");
        if (!validation.postw90_available) detail.push("Derived `postw90.x` is not executable from the configured `wannier90.x` location.");
        if (!validation.workspace_writable) detail.push("Remote workspace is not writable.");
        detail.push(...validation.messages);

        if (
          connection.success
          && validation.sbatch_available
          && validation.squeue_available
          && validation.sacct_available
          && validation.qe_pw_available
          && validation.qe_epw_available
          && validation.qe_pw2wannier_available
          && validation.wannier90_available
          && validation.postw90_available
          && validation.workspace_writable
        ) {
          setValidationMessage("Validation succeeded.");
        } else {
          setValidationMessage("Validation completed with issues.");
        }
        setValidationDetails(detail);
      }
    } catch (e) {
      setValidationMessage(`Failed to save profile: ${e}`);
    } finally {
      setIsSaving(false);
    }
  }

  async function handleImportExistingProfile() {
    if (isBusy) {
      return;
    }

    setIsImportingPresetBundle(true);
    setImportMessage(null);
    setValidationMessage(null);
    setValidationDetails([]);
    try {
      const selectedPath = await open({
        title: "Import HPC Presets + Defaults",
        directory: false,
        multiple: false,
        filters: [{ name: "QCortado HPC Preset Bundle", extensions: ["qchpc", "json"] }],
      });
      if (!selectedPath || Array.isArray(selectedPath)) {
        return;
      }

      const result = await importHpcPresetBundle(selectedPath);
      const [profiles, activeProfileId] = await Promise.all([
        listHpcProfiles(),
        getActiveHpcProfileId(),
      ]);
      const importedProfile = profiles.find((candidate) => candidate.id === activeProfileId) ?? profiles[0] ?? null;

      const summaryBase = `Imported ${result.imported_profile_count} profile preset(s): ${result.updated_profile_count} updated, ${result.created_profile_count} created.`;
      const summary = result.profiles_requiring_username.length > 0
        ? `${summaryBase} Set usernames before connecting for: ${result.profiles_requiring_username.join(", ")}.`
        : summaryBase;
      setImportMessage(summary);

      if (!importedProfile) {
        return;
      }

      setStep(0);
      setProfile(normalizeQeProfilePaths(importedProfile));
      setCredential("");
      setPersistCredential(importedProfile.credential_persisted);
      onSaved(importedProfile);
    } catch (e) {
      setImportMessage(`Failed to import presets: ${e}`);
    } finally {
      setIsImportingPresetBundle(false);
    }
  }

  if (!isOpen) return null;

  return (
    <div className="settings-window-overlay" onClick={() => !isBusy && onClose()}>
      <div
        className="floating-settings-menu hpc-setup-wizard"
        onClick={(event) => event.stopPropagation()}
        role="dialog"
        aria-modal="true"
        aria-label="HPC setup wizard"
      >
        <div className="settings-window-header">
          <h3>HPC Setup ({STEP_TITLES[step]})</h3>
          <button className="settings-window-close" onClick={onClose} disabled={isBusy} aria-label="Close setup">
            &times;
          </button>
        </div>

        <div className="hpc-step-indicator">
          {STEP_TITLES.map((title, index) => (
            <span key={title} className={index === step ? "active" : index < step ? "done" : ""}>
              {index + 1}. {title}
            </span>
          ))}
        </div>

        <div className="settings-window-content">
          {step === 0 && (
            <div className="settings-menu-section">
              <button
                className="settings-menu-item"
                onClick={() => void handleImportExistingProfile()}
                disabled={isBusy}
              >
                {isImportingPresetBundle ? "Importing..." : "Import Existing Profile"}
              </button>
              <p className="settings-menu-hint">
                Load a saved HPC preset bundle (<code>.qchpc</code> or <code>.json</code>) instead of entering all fields manually.
              </p>
              {importMessage && <div className="settings-menu-status">{importMessage}</div>}
              <label className="settings-menu-label">Profile Name</label>
              <input
                className="settings-menu-input"
                value={profile.name}
                onChange={(event) => setProfile((prev) => ({ ...prev, name: event.target.value }))}
              />
              <label className="settings-menu-label">SSH Host</label>
              <input
                className="settings-menu-input"
                value={profile.host}
                onChange={(event) => setProfile((prev) => ({ ...prev, host: event.target.value }))}
              />
              <p className="settings-menu-hint">
                Andromeda host: <code>andromeda.bc.edu</code> (requires BC network or VPN).
              </p>
              <label className="settings-menu-label">SSH Port</label>
              <input
                className="settings-menu-input"
                type="number"
                min={1}
                max={65535}
                value={profile.port}
                onChange={(event) => {
                  const value = Number.parseInt(event.target.value, 10);
                  setProfile((prev) => ({ ...prev, port: Number.isFinite(value) ? value : 22 }));
                }}
              />
              <label className="settings-menu-label">Username</label>
              <input
                className="settings-menu-input"
                value={profile.username}
                onChange={(event) => setProfile((prev) => ({ ...prev, username: event.target.value }))}
              />
            </div>
          )}

          {step === 1 && (
            <div className="settings-menu-section">
              <label className="settings-menu-label">Authentication Method</label>
              <select
                className="settings-menu-input"
                value={profile.auth_method}
                onChange={(event) => {
                  const method = event.target.value === "password" ? "password" : "ssh_key";
                  setProfile((prev) => ({ ...prev, auth_method: method }));
                }}
              >
                <option value="ssh_key">SSH Key</option>
                <option value="password">Password</option>
              </select>
              {profile.auth_method === "ssh_key" && (
                <>
                  <label className="settings-menu-label">SSH Key Path</label>
                  <input
                    className="settings-menu-input"
                    value={profile.ssh_key_path || ""}
                    onChange={(event) => setProfile((prev) => ({ ...prev, ssh_key_path: event.target.value }))}
                  />
                </>
              )}
              <label className="settings-menu-label">
                {profile.auth_method === "password" ? "Password" : "Key Passphrase (optional)"}
              </label>
              <input
                className="settings-menu-input"
                type="password"
                value={credential}
                onChange={(event) => setCredential(event.target.value)}
                placeholder={profile.auth_method === "password" ? "Enter password" : "Leave empty if ssh-agent handles auth"}
              />
              <label className="toggle-label">
                <input
                  type="checkbox"
                  checked={persistCredential}
                  onChange={(event) => setPersistCredential(event.target.checked)}
                />
                <span>Store credential in encrypted keychain</span>
              </label>
            </div>
          )}

          {step === 2 && (
            <div className="settings-menu-section">
              <div className="hpc-path-grid">
                <label className="settings-menu-label">Remote QE CPU Bin Directory</label>
                <input
                  className="settings-menu-input"
                  value={profile.remote_qe_cpu_bin_dir || ""}
                  onChange={(event) =>
                    setProfile((prev) => ({ ...prev, remote_qe_cpu_bin_dir: event.target.value || null }))
                  }
                  placeholder="~/qe/bin"
                />
                <label className="settings-menu-label">Remote QE GPU Bin Directory</label>
                <input
                  className="settings-menu-input"
                  value={profile.remote_qe_gpu_bin_dir || ""}
                  onChange={(event) =>
                    setProfile((prev) => ({ ...prev, remote_qe_gpu_bin_dir: event.target.value || null }))
                  }
                  placeholder="~/qe-gpu/bin"
                />
              </div>
              <p className="settings-menu-hint">
                CPU runs use the CPU path and GPU runs use the GPU path. Legacy fallback path is preserved automatically.
              </p>
              <label className="settings-menu-label">Remote Wannier90 Executable</label>
              <input
                className="settings-menu-input"
                value={profile.remote_wannier90_path || ""}
                onChange={(event) => setProfile((prev) => ({ ...prev, remote_wannier90_path: event.target.value || null }))}
                placeholder="wannier90.x or /path/to/wannier90.x"
              />
              <p className="settings-menu-hint">
                Leave as <code>wannier90.x</code> to resolve from the remote <code>PATH</code>, or set an absolute path.
              </p>
              <label className="settings-menu-label">Remote postw90 Executable (auto-derived)</label>
              <input
                className="settings-menu-input"
                value={deriveRemotePostw90Path(profile.remote_wannier90_path)}
                readOnly
              />
              <p className="settings-menu-hint">
                QCortado derives this automatically by replacing the executable name in the Wannier90 path with <code>postw90.x</code>.
              </p>
              <label className="settings-menu-label">Remote Pseudopotential Directory</label>
              <input
                className="settings-menu-input"
                value={profile.remote_pseudo_dir}
                onChange={(event) => setProfile((prev) => ({ ...prev, remote_pseudo_dir: event.target.value }))}
              />
              <label className="settings-menu-label">Remote Workspace Root</label>
              <input
                className="settings-menu-input"
                value={profile.remote_workspace_root}
                onChange={(event) => setProfile((prev) => ({ ...prev, remote_workspace_root: event.target.value }))}
              />
              <label className="settings-menu-label">Remote Project Root</label>
              <input
                className="settings-menu-input"
                value={profile.remote_project_root}
                onChange={(event) => setProfile((prev) => ({ ...prev, remote_project_root: event.target.value }))}
              />
            </div>
          )}

          {step === 3 && (
            <div className="settings-menu-section">
              <label className="settings-menu-label">Supported Resource Types</label>
              <select
                className="settings-menu-input"
                value={profile.resource_mode}
                onChange={(event) => {
                  const value = event.target.value;
                  const resourceMode = value === "cpu_only" || value === "gpu_only" ? value : "both";
                  setProfile((prev) => ({ ...prev, resource_mode: resourceMode }));
                }}
              >
                <option value="both">CPU + GPU</option>
                <option value="cpu_only">CPU only</option>
                <option value="gpu_only">GPU only</option>
              </select>
              <p className="settings-menu-hint">
                Controls which resource types appear in HPC run settings.
              </p>
              <label className="settings-menu-label">MPI Launcher</label>
              <select
                className="settings-menu-input"
                value={profile.launcher}
                onChange={(event) => {
                  const value = event.target.value === "mpirun" ? "mpirun" : "srun";
                  setProfile((prev) => ({ ...prev, launcher: value }));
                }}
              >
                <option value="srun">srun</option>
                <option value="mpirun">mpirun</option>
              </select>
              <label className="settings-menu-label">Launcher Extra Args</label>
              <input
                className="settings-menu-input"
                value={profile.launcher_extra_args || ""}
                onChange={(event) =>
                  setProfile((prev) => ({
                    ...prev,
                    launcher_extra_args: normalizeCliDashText(event.target.value) || null,
                  }))
                }
                placeholder="e.g. --bind-to none"
                autoCorrect="off"
                autoCapitalize="off"
                spellCheck={false}
              />
              <p className="settings-menu-hint">
                Appended to launcher command before QE executable.
              </p>
              <label className="settings-menu-label">Default CPU Partition</label>
              <input
                className="settings-menu-input"
                value={profile.default_cpu_resources.partition || "short"}
                onChange={(event) =>
                  setProfile((prev) => ({
                    ...prev,
                    default_cpu_resources: { ...prev.default_cpu_resources, partition: event.target.value },
                  }))
                }
              />
              <label className="settings-menu-label">Default CPU Walltime (HH:MM:SS)</label>
              <input
                className="settings-menu-input"
                value={profile.default_cpu_resources.walltime || "02:00:00"}
                onChange={(event) =>
                  setProfile((prev) => ({
                    ...prev,
                    default_cpu_resources: { ...prev.default_cpu_resources, walltime: event.target.value },
                  }))
                }
              />
              <label className="settings-menu-label">Default GPU Walltime (HH:MM:SS)</label>
              <input
                className="settings-menu-input"
                value={profile.default_gpu_resources.walltime || "02:00:00"}
                onChange={(event) =>
                  setProfile((prev) => ({
                    ...prev,
                    default_gpu_resources: { ...prev.default_gpu_resources, walltime: event.target.value },
                  }))
                }
              />
              <label className="settings-menu-label">Default GPU Count</label>
              <input
                className="settings-menu-input"
                type="number"
                min={1}
                value={profile.default_gpu_resources.gpus || 1}
                onChange={(event) =>
                  setProfile((prev) => ({
                    ...prev,
                    default_gpu_resources: {
                      ...prev.default_gpu_resources,
                      gpus: Number.parseInt(event.target.value, 10) || 1,
                    },
                  }))
                }
              />
            </div>
          )}

          {step === 4 && (
            <div className="settings-menu-section">
              <p className="settings-menu-hint">
                Runs checks for SSH login, Slurm commands (`sbatch`/`squeue`/`sacct`), `pw.x`, `epw.x`, `pw2wannier90.x`, `wannier90.x`, and remote workspace write permissions.
              </p>
              <button
                className="settings-menu-item"
                onClick={() => void handlePersistProfile(true)}
                disabled={isBusy}
              >
                {isSaving ? "Running checks..." : "Save Profile + Run Checks"}
              </button>
              {validationMessage && <div className="settings-menu-status">{validationMessage}</div>}
              {validationDetails.length > 0 && (
                <ul className="hpc-validation-list">
                  {validationDetails.map((detail) => (
                    <li key={detail}>{detail}</li>
                  ))}
                </ul>
              )}
            </div>
          )}
        </div>

        <div className="dialog-footer">
          <button
            className="dialog-btn cancel"
            onClick={() => setStep((prev) => Math.max(0, prev - 1))}
            disabled={step === 0 || isBusy}
          >
            Back
          </button>
          {step < STEP_TITLES.length - 1 ? (
            <button
              className="dialog-btn primary"
              onClick={() => setStep((prev) => Math.min(STEP_TITLES.length - 1, prev + 1))}
              disabled={!canContinue || isBusy}
            >
              Next
            </button>
          ) : (
            <button
              className="dialog-btn primary"
              onClick={async () => {
                if (isBusy) return;
                await handlePersistProfile(false);
                onClose();
              }}
              disabled={isBusy}
            >
              Done
            </button>
          )}
        </div>
      </div>
    </div>
  );
}
