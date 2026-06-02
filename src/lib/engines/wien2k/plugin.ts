import type {
  EnginePluginManifest,
  EngineWorkflowPanelDescriptor,
  SharedWorkflowPanelDescriptor,
  SharedWorkflowPanelKind,
  WorkflowInputRequirement,
  WorkflowInputRequirementKind,
} from "../plugin";

function sharedPanel(
  kind: SharedWorkflowPanelKind,
  label: string,
  order: number,
  required = true,
): SharedWorkflowPanelDescriptor {
  return {
    kind,
    label,
    order,
    required,
  };
}

function enginePanel(
  id: string,
  label: string,
  componentKey: string,
  order: number,
  required: boolean,
  description: string,
): EngineWorkflowPanelDescriptor {
  return {
    id,
    label,
    componentKey,
    order,
    required,
    description,
  };
}

function requirement(
  kind: WorkflowInputRequirementKind,
  label: string,
  required = true,
): WorkflowInputRequirement {
  return {
    kind,
    label,
    required,
  };
}

const structurePanels: SharedWorkflowPanelDescriptor[] = [
  sharedPanel("structure_source", "Project structure", 10),
  sharedPanel("structure_viewer", "Structure viewer", 20),
  sharedPanel("hpc_run_settings", "Remote execution settings", 90),
  sharedPanel("live_output", "Live output", 100),
  sharedPanel("result_summary", "Results", 110),
];

const sourceScfPanels: SharedWorkflowPanelDescriptor[] = [
  sharedPanel("source_calculation_selector", "Source Wien2k SCF case", 10),
  sharedPanel("hpc_run_settings", "Remote execution settings", 90),
  sharedPanel("live_output", "Live output", 100),
  sharedPanel("result_summary", "Results", 110),
];

/** Complete WIEN2k manifest, including workflows reserved for later phases. */
export const WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST: EnginePluginManifest = {
  descriptor: {
    id: "wien2k",
    label: "WIEN2k",
    status: "reserved",
    executionModes: ["remote"],
    calculationKinds: ["engine_setup", "scf", "bands", "fermi_surface", "dos"],
  },
  hpc: {
    supportsLocal: false,
    supportsHpc: true,
    supportsRemoteOnly: true,
    runtimeProfileKey: "wien2k.remote.caseRuntimeProfile",
  },
  workflows: [
    {
      kind: "engine_setup",
      label: "Wien2k case setup",
      routeKey: "wien2k.caseSetup",
      executionModes: ["remote"],
      sharedPanels: structurePanels,
      enginePanels: [
        enginePanel("wien2k-case-directory", "Case directory", "wien2k.case.directory", 30, true, "Remote WIEN2k case prefix and case-directory staging controls."),
        enginePanel("wien2k-struct-generation", "case.struct", "wien2k.case.struct", 40, true, "CIF-to-case.struct generation and validation owned by the WIEN2k engine."),
      ],
      inputRequirements: [
        requirement("crystal_structure", "Crystal structure"),
        requirement("cif_content", "CIF content"),
        requirement("engine_runtime_profile", "Wien2k remote runtime profile"),
      ],
      produces: ["native_artifacts"],
    },
    {
      kind: "scf",
      label: "Wien2k SCF",
      routeKey: "wien2k.scf",
      executionModes: ["remote"],
      sharedPanels: structurePanels,
      enginePanels: [
        enginePanel("wien2k-case-state", "Case state", "wien2k.case.state", 30, true, "Remote case-state selection, staging, and resume policy for WIEN2k."),
        enginePanel("wien2k-scf-basis", "Basis and k mesh", "wien2k.scf.basis", 40, true, "WIEN2k RKmax, Gmax, lmax, k-mesh, and LSTART controls. Accepted RMT values remain owned by case.struct."),
        enginePanel("wien2k-scf-cycle", "SCF cycle", "wien2k.scf.cycle", 50, true, "Serial WIEN2k run_lapw/runsp_lapw convergence and basic spin-mode controls."),
      ],
      inputRequirements: [
        requirement("crystal_structure", "Crystal structure"),
        requirement("cif_content", "CIF content"),
        requirement("engine_runtime_profile", "Wien2k remote runtime profile"),
        requirement("engine_case_state", "Wien2k remote case state", false),
      ],
      produces: ["scf_summary", "native_artifacts"],
    },
    {
      kind: "bands",
      label: "Wien2k bands",
      routeKey: "wien2k.bands",
      executionModes: ["remote"],
      sharedPanels: [
        sourceScfPanels[0],
        sharedPanel("k_path_viewer", "Band path", 20),
        ...sourceScfPanels.slice(1),
        sharedPanel("viewer", "Band viewer", 120),
      ],
      enginePanels: [
        enginePanel("wien2k-bands-spaghetti", "lapw1/lapw2/spaghetti", "wien2k.bands.spaghetti", 30, true, "WIEN2k band command sequence, including case.klist_band and spaghetti output parsing."),
      ],
      inputRequirements: [
        requirement("source_scf_calculation", "Source Wien2k SCF calculation"),
        requirement("engine_case_state", "Converged Wien2k case state"),
        requirement("engine_runtime_profile", "Wien2k remote runtime profile"),
      ],
      produces: ["band_dataset", "native_artifacts"],
    },
    {
      kind: "fermi_surface",
      label: "Wien2k Fermi surface",
      routeKey: "wien2k.fermiSurface",
      executionModes: ["remote"],
      sharedPanels: [
        sourceScfPanels[0],
        sharedPanel("hpc_run_settings", "Remote execution settings", 90),
        sharedPanel("live_output", "Live output", 100),
        sharedPanel("viewer", "FermiSurfer", 120),
      ],
      enginePanels: [
        enginePanel("wien2k-fermi-xcrysden", "XCrySDen BXSF", "wien2k.fermi.xcrysden", 30, true, "WIEN2k dense k mesh, lapw1/lapw2 -fermi, and XCrySDen BXSF conversion for FermiSurfer."),
      ],
      inputRequirements: [
        requirement("source_scf_calculation", "Source Wien2k SCF calculation"),
        requirement("engine_case_state", "Converged Wien2k case state"),
        requirement("engine_runtime_profile", "Wien2k remote runtime profile"),
      ],
      produces: ["native_artifacts"],
    },
    {
      kind: "dos",
      label: "Wien2k DOS",
      routeKey: "wien2k.dos",
      executionModes: ["remote"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("wien2k-dos", "DOS sequence", "wien2k.dos.sequence", 30, true, "Future WIEN2k DOS command sequence and DOS file parser boundary."),
      ],
      inputRequirements: [
        requirement("source_scf_calculation", "Source Wien2k SCF calculation"),
        requirement("engine_case_state", "Converged Wien2k case state"),
      ],
      produces: ["dos_dataset", "native_artifacts"],
    },
  ],
};

/** WIEN2k workflows currently available after remote installation. */
export const WIEN2K_STRUCTURE_ENGINE_PLUGIN_MANIFEST: EnginePluginManifest = {
  ...WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST,
  descriptor: {
    ...WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.descriptor,
    status: "configured",
    calculationKinds: ["engine_setup", "scf", "bands", "fermi_surface"],
  },
  workflows: WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.workflows.filter(
    (workflow) => workflow.kind === "engine_setup" || workflow.kind === "scf" || workflow.kind === "bands" || workflow.kind === "fermi_surface",
  ),
};
