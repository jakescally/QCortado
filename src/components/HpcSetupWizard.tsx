import { useEffect, useMemo, useRef, useState } from "react";
import {
  defaultCpuResources,
  defaultGpuResources,
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

export function HpcSetupWizard({
  isOpen,
  initialProfile = null,
  onClose,
  onSaved,
}: HpcSetupWizardProps) {
  const wasOpenRef = useRef(false);
  const [step, setStep] = useState(0);
  const [profile, setProfile] = useState<HpcProfile>(initialProfile ?? makeDefaultProfile());
  const [credential, setCredential] = useState("");
  const [persistCredential, setPersistCredential] = useState(initialProfile?.credential_persisted ?? false);
  const [isSaving, setIsSaving] = useState(false);
  const [validationMessage, setValidationMessage] = useState<string | null>(null);
  const [validationDetails, setValidationDetails] = useState<string[]>([]);

  useEffect(() => {
    const isNewlyOpened = isOpen && !wasOpenRef.current;
    if (isNewlyOpened) {
      setStep(0);
      setProfile(initialProfile ?? makeDefaultProfile());
      setCredential("");
      setPersistCredential(initialProfile?.credential_persisted ?? false);
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
      return (
        profile.remote_qe_bin_dir.trim().length > 0
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
      const saved = await saveHpcProfile(
        profile,
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
        if (!validation.workspace_writable) detail.push("Remote workspace is not writable.");
        detail.push(...validation.messages);

        if (connection.success && validation.sbatch_available && validation.squeue_available && validation.sacct_available && validation.qe_pw_available && validation.workspace_writable) {
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

  if (!isOpen) return null;

  return (
    <div className="settings-window-overlay" onClick={() => !isSaving && onClose()}>
      <div
        className="floating-settings-menu hpc-setup-wizard"
        onClick={(event) => event.stopPropagation()}
        role="dialog"
        aria-modal="true"
        aria-label="HPC setup wizard"
      >
        <div className="settings-window-header">
          <h3>HPC Setup ({STEP_TITLES[step]})</h3>
          <button className="settings-window-close" onClick={onClose} disabled={isSaving} aria-label="Close setup">
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
              <label className="settings-menu-label">Remote QE Bin Directory</label>
              <input
                className="settings-menu-input"
                value={profile.remote_qe_bin_dir}
                onChange={(event) => setProfile((prev) => ({ ...prev, remote_qe_bin_dir: event.target.value }))}
              />
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
                Runs checks for SSH login, Slurm commands (`sbatch`/`squeue`/`sacct`), `pw.x`, and remote workspace write permissions.
              </p>
              <button
                className="settings-menu-item"
                onClick={() => void handlePersistProfile(true)}
                disabled={isSaving}
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
            disabled={step === 0 || isSaving}
          >
            Back
          </button>
          {step < STEP_TITLES.length - 1 ? (
            <button
              className="dialog-btn primary"
              onClick={() => setStep((prev) => Math.min(STEP_TITLES.length - 1, prev + 1))}
              disabled={!canContinue || isSaving}
            >
              Next
            </button>
          ) : (
            <button
              className="dialog-btn primary"
              onClick={async () => {
                if (isSaving) return;
                await handlePersistProfile(false);
                onClose();
              }}
              disabled={isSaving}
            >
              Done
            </button>
          )}
        </div>
      </div>
    </div>
  );
}
