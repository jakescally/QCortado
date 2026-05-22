import type { CalculationKind, EngineDescriptor, EngineExecutionMode, EngineId } from "./types";

export type SharedWorkflowPanelKind =
  | "structure_source"
  | "structure_viewer"
  | "source_calculation_selector"
  | "k_path_viewer"
  | "hpc_run_settings"
  | "live_output"
  | "result_summary"
  | "viewer";

export type WorkflowInputRequirementKind =
  | "crystal_structure"
  | "cif_content"
  | "source_scf_calculation"
  | "source_bands_calculation"
  | "source_wannier_calculation"
  | "source_phonon_calculation"
  | "engine_runtime_profile";

export type EngineResultDatasetKind =
  | "scf_summary"
  | "band_dataset"
  | "dos_dataset"
  | "phonon_bands"
  | "phonon_dos"
  | "wannier_result"
  | "transport_dataset"
  | "epw_dataset"
  | "native_artifacts";

export interface SharedWorkflowPanelDescriptor {
  kind: SharedWorkflowPanelKind;
  label: string;
  order: number;
  required: boolean;
}

export interface EngineWorkflowPanelDescriptor {
  id: string;
  label: string;
  componentKey: string;
  order: number;
  required: boolean;
  description?: string | null;
}

export interface WorkflowInputRequirement {
  kind: WorkflowInputRequirementKind;
  label: string;
  required: boolean;
}

export interface EngineWorkflowDescriptor {
  kind: CalculationKind;
  label: string;
  routeKey: string;
  executionModes: EngineExecutionMode[];
  sharedPanels: SharedWorkflowPanelDescriptor[];
  enginePanels: EngineWorkflowPanelDescriptor[];
  inputRequirements: WorkflowInputRequirement[];
  produces: EngineResultDatasetKind[];
}

export interface EngineHpcInterfaceDescriptor {
  supportsLocal: boolean;
  supportsHpc: boolean;
  supportsRemoteOnly: boolean;
  runtimeProfileKey: string;
}

export interface EnginePluginManifest {
  descriptor: EngineDescriptor;
  hpc: EngineHpcInterfaceDescriptor;
  workflows: EngineWorkflowDescriptor[];
}

export type EngineWorkflowView =
  | "scf-wizard"
  | "bands-wizard"
  | "dos-wizard"
  | "fermi-surface-wizard"
  | "hubbard-lrt-wizard"
  | "phonon-wizard"
  | "epw-wizard"
  | "wannier-wizard"
  | "transport-wizard";

export interface FrontendEnginePlugin {
  id: EngineId;
  manifest: EnginePluginManifest;
  workflowViews: Partial<Record<CalculationKind, EngineWorkflowView>>;
  getWorkflowView: (kind: CalculationKind) => EngineWorkflowView | null;
}

