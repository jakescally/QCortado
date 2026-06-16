import type { CalculationKind, EngineDescriptor, EngineId, ImplementedEngineId } from "./types";
import type {
  EnginePluginManifest,
  EngineWorkflowPanelDescriptor,
  EngineWorkflowView,
  FrontendEnginePlugin,
  SharedWorkflowPanelDescriptor,
  SharedWorkflowPanelKind,
  WorkflowInputRequirement,
  WorkflowInputRequirementKind,
} from "./plugin";
import { WIEN2K_STRUCTURE_ENGINE_PLUGIN_MANIFEST } from "./wien2k";

export const DEFAULT_ENGINE_ID: ImplementedEngineId = "qe";

export const QE_ENGINE_DESCRIPTOR: EngineDescriptor = {
  id: "qe",
  label: "Quantum ESPRESSO",
  status: "implemented",
  executionModes: ["local", "hpc"],
  calculationKinds: [
    "scf",
    "structure_optimization",
    "bands",
    "dos",
    "fermi_surface",
    "phonon",
    "hubbard_lrt",
    "wannier",
    "transport",
    "epw",
  ],
};

export const FALLBACK_ENGINE_DESCRIPTORS: readonly EngineDescriptor[] = [
  QE_ENGINE_DESCRIPTOR,
];

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
  sharedPanel("hpc_run_settings", "Execution settings", 90),
  sharedPanel("live_output", "Live output", 100),
  sharedPanel("result_summary", "Results", 110),
];

const sourceScfPanels: SharedWorkflowPanelDescriptor[] = [
  sharedPanel("source_calculation_selector", "Source SCF calculation", 10),
  sharedPanel("hpc_run_settings", "Execution settings", 90),
  sharedPanel("live_output", "Live output", 100),
  sharedPanel("result_summary", "Results", 110),
];

export const QE_ENGINE_PLUGIN_MANIFEST: EnginePluginManifest = {
  descriptor: QE_ENGINE_DESCRIPTOR,
  hpc: {
    supportsLocal: true,
    supportsHpc: true,
    supportsRemoteOnly: false,
    runtimeProfileKey: "qe.hpc.runtimeProfile",
  },
  workflows: [
    {
      kind: "scf",
      label: "SCF",
      routeKey: "qe.scf",
      executionModes: ["local", "hpc"],
      sharedPanels: structurePanels,
      enginePanels: [
        enginePanel("qe-scf-calculation-mode", "Calculation mode", "qe.scf.calculationMode", 30, true, "QE pw.x calculation mode, including SCF and relax/vc-relax options."),
        enginePanel("qe-scf-pseudopotentials", "Pseudopotentials", "qe.scf.pseudopotentials", 40, true, "QE pseudopotential selection and SSSP/UPF-derived cutoff hints."),
        enginePanel("qe-scf-electronic-controls", "Electronic controls", "qe.scf.electronicControls", 50, true, "QE ecutwfc, ecutrho, occupations, smearing, k-points, and convergence controls."),
        enginePanel("qe-scf-hubbard-magnetism", "Hubbard and magnetism", "qe.scf.hubbardMagnetism", 60, false, "QE DFT+U/HUBBARD syntax and spin initialization controls."),
      ],
      inputRequirements: [
        requirement("crystal_structure", "Crystal structure"),
        requirement("cif_content", "CIF content"),
        requirement("engine_runtime_profile", "QE runtime profile"),
      ],
      produces: ["scf_summary", "native_artifacts"],
    },
    {
      kind: "structure_optimization",
      label: "Structure optimization",
      routeKey: "qe.scf",
      executionModes: ["local", "hpc"],
      sharedPanels: structurePanels,
      enginePanels: [
        enginePanel("qe-relax-mode", "Relaxation mode", "qe.scf.calculationMode", 30, true, "QE relax/vc-relax mode and optimization constraints."),
        enginePanel("qe-relax-pseudopotentials", "Pseudopotentials", "qe.scf.pseudopotentials", 40, true, "QE pseudopotential selection and cutoff hints for relaxation."),
        enginePanel("qe-relax-electronic-controls", "Electronic controls", "qe.scf.electronicControls", 50, true, "QE electronic and ionic convergence controls."),
      ],
      inputRequirements: [
        requirement("crystal_structure", "Crystal structure"),
        requirement("cif_content", "CIF content"),
        requirement("engine_runtime_profile", "QE runtime profile"),
      ],
      produces: ["scf_summary", "native_artifacts"],
    },
    {
      kind: "bands",
      label: "Bands",
      routeKey: "qe.bands",
      executionModes: ["local", "hpc"],
      sharedPanels: [
        sourceScfPanels[0],
        sharedPanel("k_path_viewer", "Band path", 20),
        ...sourceScfPanels.slice(1),
        sharedPanel("viewer", "Band viewer", 120),
      ],
      enginePanels: [
        enginePanel("qe-bands-nscf", "NSCF and bands.x", "qe.bands.nscfBands", 30, true, "QE pw.x NSCF and bands.x controls."),
        enginePanel("qe-bands-projections", "Projections", "qe.bands.projections", 40, false, "Optional QE projwfc.x fat-band projection controls."),
      ],
      inputRequirements: [
        requirement("source_scf_calculation", "Source SCF calculation"),
        requirement("engine_runtime_profile", "QE runtime profile"),
      ],
      produces: ["band_dataset", "native_artifacts"],
    },
    {
      kind: "dos",
      label: "DOS",
      routeKey: "qe.dos",
      executionModes: ["local", "hpc"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("qe-dos-projwfc", "DOS and projections", "qe.dos.projwfc", 30, true, "QE NSCF, dos.x, and projwfc.x DOS controls."),
      ],
      inputRequirements: [
        requirement("source_scf_calculation", "Source SCF calculation"),
        requirement("engine_runtime_profile", "QE runtime profile"),
      ],
      produces: ["dos_dataset"],
    },
    {
      kind: "fermi_surface",
      label: "Fermi surface",
      routeKey: "qe.fermiSurface",
      executionModes: ["local", "hpc"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("qe-fermi-surface", "Fermi surface", "qe.fermiSurface.fs", 30, true, "QE dense NSCF and Fermi-surface artifact generation controls."),
      ],
      inputRequirements: [requirement("source_scf_calculation", "Source SCF calculation")],
      produces: ["native_artifacts"],
    },
    {
      kind: "phonon",
      label: "Phonons",
      routeKey: "qe.phonon",
      executionModes: ["local", "hpc"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("qe-phonon-pipeline", "Phonon pipeline", "qe.phonon.pipeline", 30, true, "QE ph.x, q2r.x, and matdyn.x controls."),
      ],
      inputRequirements: [requirement("source_scf_calculation", "Source SCF calculation")],
      produces: ["phonon_bands", "phonon_dos", "native_artifacts"],
    },
    {
      kind: "hubbard_lrt",
      label: "Hubbard LRT",
      routeKey: "qe.hubbardLrt",
      executionModes: ["local", "hpc"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("qe-hubbard-lrt", "hp.x linear response", "qe.hubbardLrt.hp", 30, true, "QE hp.x q-mesh, perturbation, and Hubbard parameter controls."),
      ],
      inputRequirements: [requirement("source_scf_calculation", "DFT+U source SCF calculation")],
      produces: ["native_artifacts"],
    },
    {
      kind: "wannier",
      label: "Wannier",
      routeKey: "qe.wannier",
      executionModes: ["local", "hpc"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("qe-wannier-nscf", "NSCF and interfaces", "qe.wannier.nscfInterfaces", 30, true, "QE NSCF, pw2wannier90.x, and wannier90.x setup controls."),
        enginePanel("qe-wannier-projections", "Wannier projections", "qe.wannier.projections", 40, true, "QE/Wannier90 projection and disentanglement controls."),
      ],
      inputRequirements: [requirement("source_scf_calculation", "Source SCF calculation")],
      produces: ["wannier_result", "band_dataset", "native_artifacts"],
    },
    {
      kind: "transport",
      label: "Transport",
      routeKey: "qe.transport",
      executionModes: ["local", "hpc"],
      sharedPanels: [
        sharedPanel("source_calculation_selector", "Source Wannier calculation", 10),
        sharedPanel("hpc_run_settings", "Execution settings", 90),
        sharedPanel("live_output", "Live output", 100),
        sharedPanel("viewer", "Transport viewer", 120),
      ],
      enginePanels: [
        enginePanel("qe-transport-boltz", "BoltzWann transport", "qe.transport.boltzWann", 30, true, "QE/Wannier90 BoltzWann and postw90 transport controls."),
      ],
      inputRequirements: [requirement("source_wannier_calculation", "Source Wannier calculation")],
      produces: ["transport_dataset"],
    },
    {
      kind: "epw",
      label: "EPW",
      routeKey: "qe.epw",
      executionModes: ["local", "hpc"],
      sharedPanels: sourceScfPanels,
      enginePanels: [
        enginePanel("qe-epw", "EPW", "qe.epw.workflow", 30, true, "QE epw.x source selection, mesh, transport, and superconductivity controls."),
      ],
      inputRequirements: [
        requirement("source_scf_calculation", "Source SCF/NSCF calculation"),
        requirement("source_phonon_calculation", "Source phonon calculation", false),
      ],
      produces: ["epw_dataset", "native_artifacts"],
    },
  ],
};

export const FALLBACK_ENGINE_PLUGIN_MANIFESTS: readonly EnginePluginManifest[] = [
  QE_ENGINE_PLUGIN_MANIFEST,
];

export const QE_WORKFLOW_VIEWS: Partial<Record<CalculationKind, EngineWorkflowView>> = {
  scf: "scf-wizard",
  structure_optimization: "scf-wizard",
  bands: "bands-wizard",
  dos: "dos-wizard",
  fermi_surface: "fermi-surface-wizard",
  hubbard_lrt: "hubbard-lrt-wizard",
  phonon: "phonon-wizard",
  epw: "epw-wizard",
  wannier: "wannier-wizard",
  transport: "transport-wizard",
};

export const QE_FRONTEND_ENGINE_PLUGIN: FrontendEnginePlugin = {
  id: "qe",
  manifest: QE_ENGINE_PLUGIN_MANIFEST,
  workflowViews: QE_WORKFLOW_VIEWS,
  getWorkflowView(kind) {
    return QE_WORKFLOW_VIEWS[kind] ?? null;
  },
};

export const WIEN2K_WORKFLOW_VIEWS: Partial<Record<CalculationKind, EngineWorkflowView>> = {
  engine_setup: "wien2k-structure-wizard",
  scf: "wien2k-scf-wizard",
  soc: "wien2k-soc-wizard",
  bands: "bands-wizard",
  fermi_surface: "fermi-surface-wizard",
};

export const WIEN2K_FRONTEND_ENGINE_PLUGIN: FrontendEnginePlugin = {
  id: "wien2k",
  manifest: WIEN2K_STRUCTURE_ENGINE_PLUGIN_MANIFEST,
  workflowViews: WIEN2K_WORKFLOW_VIEWS,
  getWorkflowView(kind) {
    return WIEN2K_WORKFLOW_VIEWS[kind] ?? null;
  },
};

export function isImplementedEngineId(engineId: EngineId): engineId is ImplementedEngineId {
  return engineId === DEFAULT_ENGINE_ID;
}

export function getFrontendEnginePlugin(engineId: EngineId): FrontendEnginePlugin | null {
  if (engineId === "qe") {
    return QE_FRONTEND_ENGINE_PLUGIN;
  }
  if (engineId === "wien2k") {
    return WIEN2K_FRONTEND_ENGINE_PLUGIN;
  }
  return null;
}

export function getEngineWorkflowView(
  engineId: EngineId,
  kind: CalculationKind,
): EngineWorkflowView | null {
  return getFrontendEnginePlugin(engineId)?.getWorkflowView(kind) ?? null;
}

export function getEngineLabel(
  descriptors: readonly EngineDescriptor[],
  engineId: EngineId,
): string {
  return descriptors.find((descriptor) => descriptor.id === engineId)?.label ?? engineId;
}

export function getEngineShortLabel(engine: EngineDescriptor): string {
  if (engine.id === "qe") {
    return "QE";
  }
  if (engine.id === "wien2k") {
    return "WIEN2k";
  }
  return engine.label;
}
