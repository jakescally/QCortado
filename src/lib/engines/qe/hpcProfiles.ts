import {
  buildHpcLauncherCommand,
  supportsProfileResourceType,
} from "../../hpcConfig";
import {
  resolveProfileRemotePseudoDir,
  resolveProfileRemoteQeBinDir,
} from "./hpc";
import type {
  HpcLauncher,
  HpcProfile,
  HpcResourceMode,
  HpcResourceType,
  SlurmResourceRequest,
} from "../../types";
import type { EngineId } from "../types";

export interface HpcRuntimeProfileOptions {
  resourceType?: HpcResourceType | null;
}

export interface EngineHpcRuntimeProfileBase<Engine extends EngineId> {
  engineId: Engine;
  profileId: string;
  profileName: string;
  cluster: string;
  host: string;
  port: number;
  username: string;
  resourceMode: HpcResourceMode;
  resourceType: HpcResourceType;
  launcher: HpcLauncher;
  launcherCommand: string;
  remoteWorkspaceRoot: string;
  remoteProjectRoot: string;
  defaultResources: SlurmResourceRequest;
}

export interface QeHpcRuntimePaths {
  remoteBinDir: string;
  remotePseudoDir: string;
  remoteEpwPath: string | null;
  remoteWannier90Path: string | null;
  remotePostw90Path: string | null;
  usesPseudopotentials: true;
}

export interface QeHpcRuntimeProfile extends EngineHpcRuntimeProfileBase<"qe"> {
  supported: true;
  supportsResourceType: boolean;
  paths: QeHpcRuntimePaths;
}

export interface UnsupportedEngineHpcRuntimeProfile {
  engineId: Exclude<EngineId, "qe">;
  supported: false;
  profileId: string;
  profileName: string;
  reason: string;
}

export type EngineHpcRuntimeProfile =
  | QeHpcRuntimeProfile
  | UnsupportedEngineHpcRuntimeProfile;

function normalizeOptionalText(value: string | null | undefined): string | null {
  if (typeof value !== "string") {
    return null;
  }
  const trimmed = value.trim();
  return trimmed.length > 0 ? trimmed : null;
}

function preferredResourceType(profile: HpcProfile): HpcResourceType {
  return profile.resource_mode === "gpu_only" ? "gpu" : "cpu";
}

function cloneResourceRequest(
  source: SlurmResourceRequest,
  resourceType: HpcResourceType,
): SlurmResourceRequest {
  return {
    ...source,
    resource_type: resourceType,
    additional_sbatch: [...(source.additional_sbatch || [])],
  };
}

function defaultResourcesForResourceType(
  profile: HpcProfile,
  resourceType: HpcResourceType,
): SlurmResourceRequest {
  return resourceType === "gpu"
    ? cloneResourceRequest(profile.default_gpu_resources, "gpu")
    : cloneResourceRequest(profile.default_cpu_resources, "cpu");
}

export function qeHpcProfileToRuntimeProfile(
  profile: HpcProfile,
  options: HpcRuntimeProfileOptions = {},
): QeHpcRuntimeProfile {
  const resourceType = options.resourceType ?? preferredResourceType(profile);

  return {
    engineId: "qe",
    supported: true,
    profileId: profile.id,
    profileName: profile.name,
    cluster: profile.cluster,
    host: profile.host,
    port: profile.port,
    username: profile.username,
    resourceMode: profile.resource_mode,
    resourceType,
    supportsResourceType: supportsProfileResourceType(profile, resourceType),
    launcher: profile.launcher,
    launcherCommand: buildHpcLauncherCommand(profile, resourceType),
    remoteWorkspaceRoot: profile.remote_workspace_root,
    remoteProjectRoot: profile.remote_project_root,
    defaultResources: defaultResourcesForResourceType(profile, resourceType),
    paths: {
      remoteBinDir: resolveProfileRemoteQeBinDir(profile, resourceType),
      remotePseudoDir: resolveProfileRemotePseudoDir(profile, resourceType),
      remoteEpwPath: normalizeOptionalText(profile.remote_epw_path),
      remoteWannier90Path: normalizeOptionalText(profile.remote_wannier90_path),
      remotePostw90Path: normalizeOptionalText(profile.remote_postw90_path),
      usesPseudopotentials: true,
    },
  };
}

export function hpcProfileToEngineRuntimeProfile(
  profile: HpcProfile,
  engineId: EngineId = "qe",
  options: HpcRuntimeProfileOptions = {},
): EngineHpcRuntimeProfile {
  if (engineId === "qe") {
    return qeHpcProfileToRuntimeProfile(profile, options);
  }

  return {
    engineId,
    supported: false,
    profileId: profile.id,
    profileName: profile.name,
    reason: "This engine does not have an HPC runtime profile adapter yet.",
  };
}
