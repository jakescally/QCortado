import { invoke } from "@tauri-apps/api/core";

import { QE_ENGINE_DESCRIPTOR } from "./registry";
import type { EngineDescriptor, EngineId, EngineImplementationStatus } from "./types";
import { WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST } from "./wien2k";
import type { HpcProfile } from "../types";

export interface EngineInstallation {
  engineId: EngineId;
  hpcProfileId: string;
  remoteInstallRoot: string;
  remoteWorkspaceRoot: string;
  remoteProjectRoot: string;
  verifiedExecutables: string[];
  versionHint?: string | null;
  verifiedAt: string;
}

export interface EngineInstallationRequest {
  engineId: EngineId;
  hpcProfileId: string;
  remoteInstallRoot: string;
}

export interface EngineInstallationVerification {
  success: boolean;
  message: string;
  checkedExecutables: string[];
  versionHint?: string | null;
}

export interface AddEngineInstallationResult {
  installation: EngineInstallation;
  verification: EngineInstallationVerification;
}

export interface EngineInstallFormDefaults {
  hpcProfileId: string;
  remoteInstallRoot: string;
}

export const INSTALLABLE_ENGINE_DESCRIPTORS: readonly EngineDescriptor[] = [
  QE_ENGINE_DESCRIPTOR,
  WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.descriptor,
];

export function isSelectableEngineStatus(status: EngineImplementationStatus): boolean {
  return status === "implemented" || status === "configured";
}

export function isEngineAlreadyAvailable(
  engineId: EngineId,
  descriptors: readonly EngineDescriptor[],
  installations: readonly EngineInstallation[] = [],
): boolean {
  return (
    descriptors.some(
      (descriptor) => descriptor.id === engineId && isSelectableEngineStatus(descriptor.status),
    ) || installations.some((installation) => installation.engineId === engineId)
  );
}

export function buildDefaultEngineInstallForm(
  profile: HpcProfile | null | undefined,
): EngineInstallFormDefaults {
  return {
    hpcProfileId: profile?.id ?? "",
    remoteInstallRoot: "",
  };
}

export async function listEngineInstallations(): Promise<EngineInstallation[]> {
  return invoke<EngineInstallation[]>("list_engine_installations");
}

export async function addEngineInstallation(
  request: EngineInstallationRequest,
): Promise<AddEngineInstallationResult> {
  return invoke<AddEngineInstallationResult>("add_engine_installation", { request });
}
