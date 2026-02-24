import { invoke } from "@tauri-apps/api/core";
import {
  ExecutionMode,
  ExecutionTarget,
  HpcLauncher,
  HpcProfile,
  HpcResourceMode,
  SlurmResourceRequest,
} from "./types";

export interface HpcConnectionTestResult {
  success: boolean;
  message: string;
}

export interface HpcEnvironmentValidation {
  reachable: boolean;
  sbatch_available: boolean;
  squeue_available: boolean;
  sacct_available: boolean;
  qe_pw_available: boolean;
  workspace_writable: boolean;
  messages: string[];
}

export interface HpcScriptPreviewResult {
  script: string;
  sbatch_preview: string;
  validation: {
    errors: string[];
    warnings: string[];
  };
}

export interface HpcUtilizationSample {
  captured_at: string;
  source: string;
  node?: string | null;
  output: string;
}

export interface HpcArtifactSyncReport {
  mode: string;
  downloaded_files: number;
  downloaded_bytes: number;
  skipped_files: number;
  skipped_bytes: number;
}

export interface HpcPresetBundleExportResult {
  bundle_path: string;
  profile_count: number;
}

export interface HpcPresetBundleImportResult {
  imported_profile_count: number;
  updated_profile_count: number;
  created_profile_count: number;
  profiles_requiring_username: string[];
  active_profile_id?: string | null;
}

export function normalizeCliDashText(input: string): string {
  return input
    .replace(/\u2014/g, "--")
    .replace(/[\u2013\u2010\u2011\u2012\u2212\uFE63\uFF0D]/g, "-");
}

export function defaultCpuResources(): SlurmResourceRequest {
  return {
    resource_type: "cpu",
    partition: "short",
    walltime: "02:00:00",
    nodes: 1,
    ntasks: 4,
    cpus_per_task: 1,
    memory_gb: 16,
    gpus: 0,
    additional_sbatch: [],
  };
}

export function defaultGpuResources(): SlurmResourceRequest {
  return {
    resource_type: "gpu",
    partition: "short",
    walltime: "02:00:00",
    nodes: 1,
    ntasks: 1,
    cpus_per_task: 8,
    memory_gb: 32,
    gpus: 1,
    additional_sbatch: [],
  };
}

function cloneResourceTemplate(
  source: SlurmResourceRequest,
  resourceType: "cpu" | "gpu",
): SlurmResourceRequest {
  return {
    ...source,
    resource_type: resourceType,
    additional_sbatch: [...(source.additional_sbatch || [])],
  };
}

export function supportsProfileResourceType(
  profile: HpcProfile | null | undefined,
  resourceType: "cpu" | "gpu",
): boolean {
  if (!profile) {
    return true;
  }
  if (profile.resource_mode === "both") {
    return true;
  }
  return profile.resource_mode === "cpu_only" ? resourceType === "cpu" : resourceType === "gpu";
}

export function defaultResourcesForProfile(profile: HpcProfile | null | undefined): SlurmResourceRequest {
  if (!profile) {
    return defaultCpuResources();
  }
  if (profile.resource_mode === "gpu_only") {
    return cloneResourceTemplate(profile.default_gpu_resources, "gpu");
  }
  return cloneResourceTemplate(profile.default_cpu_resources, "cpu");
}

export function buildHpcLauncherCommand(profile: HpcProfile | null | undefined): string {
  const launcherBase = profile?.launcher === "mpirun"
    ? "mpirun -np \"${SLURM_NTASKS:-1}\""
    : "srun";
  const extraArgs = (profile?.launcher_extra_args || "").trim();
  if (!extraArgs) {
    return launcherBase;
  }
  return `${launcherBase} ${extraArgs}`;
}

export function buildHpcQeInputCommandLine(
  profile: HpcProfile | null | undefined,
  executable: string,
  inputFile: string,
  outputFile: string,
  extraArgs?: string,
): string {
  const launcher = buildHpcLauncherCommand(profile);
  const args = (extraArgs || "").trim();
  const argSegment = args.length > 0 ? ` ${args}` : "";
  return `${launcher} "$QE_BIN/${executable}"${argSegment} -in ${inputFile} > ${outputFile} 2>&1`;
}

export async function loadExecutionMode(): Promise<ExecutionMode> {
  try {
    return await invoke<ExecutionMode>("get_execution_mode");
  } catch {
    return "local";
  }
}

export async function saveExecutionMode(mode: ExecutionMode): Promise<void> {
  await invoke("set_execution_mode", { mode });
}

export async function listHpcProfiles(): Promise<HpcProfile[]> {
  return invoke<HpcProfile[]>("hpc_list_profiles");
}

export async function getActiveHpcProfileId(): Promise<string | null> {
  return invoke<string | null>("hpc_get_active_profile_id");
}

export async function setActiveHpcProfile(profileId: string): Promise<void> {
  await invoke("hpc_set_active_profile", { profileId });
}

export async function saveHpcProfile(
  profile: HpcProfile,
  credential: string | null,
  persistCredential: boolean,
): Promise<HpcProfile> {
  return invoke<HpcProfile>("hpc_save_profile", {
    profile,
    credential,
    persistCredential,
  });
}

export async function updateHpcProfileDefaults(
  profileId: string,
  defaultCpuResources: SlurmResourceRequest,
  defaultGpuResources: SlurmResourceRequest,
  resourceMode?: HpcResourceMode,
  launcher?: HpcLauncher,
  launcherExtraArgs?: string | null,
): Promise<HpcProfile> {
  return invoke<HpcProfile>("hpc_update_profile_defaults", {
    profileId,
    resourceMode: resourceMode ?? null,
    launcher: launcher ?? null,
    launcherExtraArgs: launcherExtraArgs ?? null,
    defaultCpuResources,
    defaultGpuResources,
  });
}

export async function deleteHpcProfile(profileId: string): Promise<void> {
  await invoke("hpc_delete_profile", { profileId });
}

export async function exportHpcPresetBundle(
  destinationPath: string,
): Promise<HpcPresetBundleExportResult> {
  return invoke<HpcPresetBundleExportResult>("hpc_export_preset_bundle", {
    destinationPath,
  });
}

export async function importHpcPresetBundle(bundlePath: string): Promise<HpcPresetBundleImportResult> {
  return invoke<HpcPresetBundleImportResult>("hpc_import_preset_bundle", {
    bundlePath,
  });
}

export async function testHpcConnection(profileId?: string | null): Promise<HpcConnectionTestResult> {
  return invoke<HpcConnectionTestResult>("hpc_test_connection", {
    profileId: profileId ?? null,
  });
}

export async function validateHpcEnvironment(profileId?: string | null): Promise<HpcEnvironmentValidation> {
  return invoke<HpcEnvironmentValidation>("hpc_validate_environment", {
    profileId: profileId ?? null,
  });
}

export async function listRemotePseudopotentials(
  pseudoDir?: string | null,
  profileId?: string | null,
): Promise<string[]> {
  return invoke<string[]>("hpc_list_remote_pseudopotentials", {
    profileId: profileId ?? null,
    pseudoDir: pseudoDir ?? null,
  });
}

export async function loadRemoteSsspData<T = Record<string, unknown>>(
  pseudoDir?: string | null,
  profileId?: string | null,
): Promise<T> {
  return invoke<T>("hpc_load_remote_sssp_data", {
    profileId: profileId ?? null,
    pseudoDir: pseudoDir ?? null,
  });
}

export async function previewSlurmScript(
  taskKind: string,
  commandLines: string[],
  resources: SlurmResourceRequest,
  profileId?: string | null,
): Promise<HpcScriptPreviewResult> {
  return invoke<HpcScriptPreviewResult>("hpc_preview_slurm_script", {
    profileId: profileId ?? null,
    taskKind,
    commandLines,
    resources,
  });
}

export async function openHpcActivityWindow(): Promise<void> {
  await invoke("hpc_open_activity_window");
}

export async function sampleHpcUtilization(
  profileId?: string | null,
  remoteJobId?: string | null,
  remoteNode?: string | null,
): Promise<HpcUtilizationSample> {
  return invoke<HpcUtilizationSample>("hpc_sample_utilization", {
    profileId: profileId ?? null,
    remoteJobId: remoteJobId ?? null,
    remoteNode: remoteNode ?? null,
  });
}

export async function downloadHpcTaskArtifacts(
  taskId: string,
  full = true,
): Promise<HpcArtifactSyncReport> {
  return invoke<HpcArtifactSyncReport>("hpc_download_task_artifacts", {
    taskId,
    full,
  });
}

export async function downloadHpcCalculationArtifacts(
  projectId: string,
  calcId: string,
  profileId?: string | null,
  full = true,
): Promise<HpcArtifactSyncReport> {
  return invoke<HpcArtifactSyncReport>("hpc_download_calculation_artifacts", {
    projectId,
    calcId,
    profileId: profileId ?? null,
    full,
  });
}

export function buildExecutionTarget(
  mode: ExecutionMode,
  profileId: string | null,
  resources: SlurmResourceRequest | null,
  interactiveDebug = false,
): ExecutionTarget {
  if (mode !== "hpc") {
    return { mode: "local" };
  }
  return {
    mode: "hpc",
    hpc: {
      profile_id: profileId,
      resources,
      interactive_debug: interactiveDebug,
    },
  };
}
