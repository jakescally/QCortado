import { invoke } from "@tauri-apps/api/core";
import {
  ExecutionMode,
  ExecutionTarget,
  HpcLauncher,
  HpcProfile,
  HpcResourceType,
  HpcResourceMode,
  PseudopotentialMetadata,
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
  qe_epw_available: boolean;
  qe_pw2wannier_available: boolean;
  wannier90_available: boolean;
  postw90_available: boolean;
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
  resource_type: HpcResourceType;
  job_id?: string | null;
  node?: string | null;
  sources?: string[];
  warnings?: string[];
  scheduler?: HpcSchedulerTelemetry | null;
  cpu?: HpcCpuTelemetry | null;
  memory?: HpcMemoryTelemetry | null;
  gpu?: HpcGpuTelemetry | null;
  raw?: string | null;
}

export interface HpcSchedulerTelemetry {
  state?: string | null;
  node?: string | null;
  allocated_cpus?: number | null;
  requested_memory_bytes?: number | null;
  nodes?: number | null;
  elapsed?: string | null;
  time_limit?: string | null;
  reason?: string | null;
}

export interface HpcCpuTelemetry {
  allocated_cpus?: number | null;
  total_cpu_seconds?: number | null;
  steps?: HpcCpuStepTelemetry[];
}

export interface HpcCpuStepTelemetry {
  job_step: string;
  n_tasks?: number | null;
  average_cpu_seconds?: number | null;
  average_cpu_display?: string | null;
  average_rss_bytes?: number | null;
  peak_rss_bytes?: number | null;
  average_vm_bytes?: number | null;
  peak_vm_bytes?: number | null;
}

export interface HpcMemoryTelemetry {
  source: string;
  current_rss_bytes?: number | null;
  peak_rss_bytes?: number | null;
  average_vm_bytes?: number | null;
  peak_vm_bytes?: number | null;
  requested_bytes?: number | null;
}

export interface HpcGpuTelemetry {
  source: string;
  devices?: HpcGpuDeviceTelemetry[];
  average_utilization_percent?: number | null;
  memory_used_bytes?: number | null;
  memory_total_bytes?: number | null;
  memory_used_percent?: number | null;
}

export interface HpcGpuDeviceTelemetry {
  index: number;
  name: string;
  utilization_percent?: number | null;
  memory_utilization_percent?: number | null;
  memory_used_bytes?: number | null;
  memory_total_bytes?: number | null;
  temperature_c?: number | null;
}

export type HpcQueueScope = "all" | "mine";

export interface HpcNodeSnapshot {
  node_name: string;
  partitions: string[];
  state: string;
  raw_state: string;
  cpus: number;
  cpu_name?: string | null;
  memory_mb: number;
  free_memory_mb?: number | null;
  gpus: number;
  gpu_name?: string | null;
  gres?: string | null;
  features: string[];
  reason?: string | null;
  node_type: string;
}

export interface HpcQueueSnapshot {
  job_id: string;
  user: string;
  state: string;
  raw_state: string;
  partition: string;
  nodes: number;
  elapsed: string;
  time_limit: string;
  reason: string;
  nodelist: string;
  name: string;
}

export interface HpcClusterSnapshot {
  captured_at: string;
  cluster: string;
  host: string;
  queue_scope: HpcQueueScope | string;
  queue_included: boolean;
  queue_limit: number;
  nodes: HpcNodeSnapshot[];
  queue: HpcQueueSnapshot[];
  warnings: string[];
}

export interface HpcArtifactSyncReport {
  mode: string;
  downloaded_files: number;
  downloaded_bytes: number;
  skipped_files: number;
  skipped_bytes: number;
}

export interface HpcRemoteOrphanCleanupResult {
  profile_id: string;
  scanned_paths: number;
  referenced_paths: number;
  orphan_paths: number;
  removed_paths: string[];
  failed_paths: string[];
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

function normalizeRemoteQeBinDir(value: string | null | undefined): string | null {
  if (typeof value !== "string") {
    return null;
  }
  const trimmed = value.trim();
  return trimmed.length > 0 ? trimmed : null;
}

export function resolveProfileRemoteQeBinDir(
  profile: HpcProfile | null | undefined,
  resourceType?: HpcResourceType | null,
): string {
  const fallback = normalizeRemoteQeBinDir(profile?.remote_qe_bin_dir) || "$HOME/qe/bin";
  if (!profile) {
    return fallback;
  }
  const requestedType = resourceType
    || (profile.resource_mode === "gpu_only" ? "gpu" : "cpu");
  if (requestedType === "gpu") {
    return normalizeRemoteQeBinDir(profile.remote_qe_gpu_bin_dir) || fallback;
  }
  return normalizeRemoteQeBinDir(profile.remote_qe_cpu_bin_dir) || fallback;
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

const QE_PENCIL_DECOMPOSITION_EXECUTABLES = new Set([
  "pw.x",
  "bands.x",
  "projwfc.x",
  "dos.x",
  "fermi_velocity.x",
  "ph.x",
  "q2r.x",
  "matdyn.x",
  "pw2wannier90.x",
]);

function commandArgsIncludePencilDecomposition(rawArgs: string): boolean {
  const normalized = rawArgs
    .replace(/\u2014/g, "--")
    .replace(/[\u2013\u2010\u2011\u2012\u2212\uFE63\uFF0D]/g, "-");
  return normalized
    .split(/\s+/)
    .some((token) => {
      const lower = token.toLowerCase();
      return lower === "-pd"
        || lower === "-use_pd"
        || lower === "-pencil_decomposition"
        || lower === "-use_pencil_decomposition"
        || lower.startsWith("-pd=")
        || lower.startsWith("-use_pd=")
        || lower.startsWith("-pencil_decomposition=")
        || lower.startsWith("-use_pencil_decomposition=");
    });
}

function qeExecutableUsesPencilDecomposition(executable: string): boolean {
  const executableName = executable.trim().split(/[\\/]/).pop()?.toLowerCase() || "";
  return QE_PENCIL_DECOMPOSITION_EXECUTABLES.has(executableName);
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
  const hasPencilDecomposition = commandArgsIncludePencilDecomposition(args);
  const pdSegment = qeExecutableUsesPencilDecomposition(executable) && !hasPencilDecomposition
    ? " -pd .true."
    : "";
  return `${launcher} "$QE_BIN/${executable}"${argSegment}${pdSegment} -in ${inputFile} > ${outputFile} 2>&1`;
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

export async function migrateHpcRemoteRoots(
  profileId: string,
  newWorkspaceRoot: string,
  newProjectRoot: string,
): Promise<HpcProfile> {
  return invoke<HpcProfile>("hpc_migrate_remote_roots", {
    profileId,
    newWorkspaceRoot,
    newProjectRoot,
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

export async function listRemotePseudopotentialMetadata(
  pseudoDir?: string | null,
  profileId?: string | null,
): Promise<PseudopotentialMetadata[]> {
  return invoke<PseudopotentialMetadata[]>("hpc_list_remote_pseudopotential_metadata", {
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
  resourceType?: HpcResourceType | null,
): Promise<HpcUtilizationSample> {
  return invoke<HpcUtilizationSample>("hpc_sample_utilization", {
    profileId: profileId ?? null,
    remoteJobId: remoteJobId ?? null,
    remoteNode: remoteNode ?? null,
    resourceType: resourceType ?? null,
  });
}

export async function getHpcClusterSnapshot(
  profileId?: string | null,
  queueScope: HpcQueueScope = "all",
  includeQueue = true,
  queueLimit = 1500,
): Promise<HpcClusterSnapshot> {
  return invoke<HpcClusterSnapshot>("hpc_get_cluster_snapshot", {
    profileId: profileId ?? null,
    queueScope,
    includeQueue,
    queueLimit,
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

export async function cleanHpcRemoteOrphans(
  profileId?: string | null,
): Promise<HpcRemoteOrphanCleanupResult> {
  return invoke<HpcRemoteOrphanCleanupResult>("hpc_clean_remote_orphans", {
    profileId: profileId ?? null,
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
