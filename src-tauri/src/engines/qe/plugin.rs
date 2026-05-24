//! QE implementation of the engine plugin contract.

use crate::engines::common::qe_engine_descriptor;
use crate::engines::plugin::{
    EngineHpcInterfaceDescriptor, EnginePlugin, EnginePluginManifest, EngineResultDatasetKind,
    EngineWorkflowDescriptor, EngineWorkflowPanelDescriptor, SharedWorkflowPanelDescriptor,
    SharedWorkflowPanelKind, WorkflowInputRequirement, WorkflowInputRequirementKind,
};
use crate::engines::types::{CalculationKind, EngineDescriptor, EngineExecutionMode, EngineId};

pub struct QeEnginePlugin;

pub static QE_ENGINE_PLUGIN: QeEnginePlugin = QeEnginePlugin;

fn shared_panel(
    kind: SharedWorkflowPanelKind,
    label: &str,
    order: u16,
    required: bool,
) -> SharedWorkflowPanelDescriptor {
    SharedWorkflowPanelDescriptor {
        kind,
        label: label.to_string(),
        order,
        required,
    }
}

fn engine_panel(
    id: &str,
    label: &str,
    component_key: &str,
    order: u16,
    required: bool,
    description: &str,
) -> EngineWorkflowPanelDescriptor {
    EngineWorkflowPanelDescriptor {
        id: id.to_string(),
        label: label.to_string(),
        component_key: component_key.to_string(),
        order,
        required,
        description: Some(description.to_string()),
    }
}

fn requirement(
    kind: WorkflowInputRequirementKind,
    label: &str,
    required: bool,
) -> WorkflowInputRequirement {
    WorkflowInputRequirement {
        kind,
        label: label.to_string(),
        required,
    }
}

fn common_structure_panels() -> Vec<SharedWorkflowPanelDescriptor> {
    vec![
        shared_panel(
            SharedWorkflowPanelKind::StructureSource,
            "Project structure",
            10,
            true,
        ),
        shared_panel(
            SharedWorkflowPanelKind::StructureViewer,
            "Structure viewer",
            20,
            true,
        ),
        shared_panel(
            SharedWorkflowPanelKind::HpcRunSettings,
            "Execution settings",
            90,
            true,
        ),
        shared_panel(
            SharedWorkflowPanelKind::LiveOutput,
            "Live output",
            100,
            true,
        ),
        shared_panel(SharedWorkflowPanelKind::ResultSummary, "Results", 110, true),
    ]
}

fn source_scf_panels() -> Vec<SharedWorkflowPanelDescriptor> {
    vec![
        shared_panel(
            SharedWorkflowPanelKind::SourceCalculationSelector,
            "Source SCF calculation",
            10,
            true,
        ),
        shared_panel(
            SharedWorkflowPanelKind::HpcRunSettings,
            "Execution settings",
            90,
            true,
        ),
        shared_panel(
            SharedWorkflowPanelKind::LiveOutput,
            "Live output",
            100,
            true,
        ),
        shared_panel(SharedWorkflowPanelKind::ResultSummary, "Results", 110, true),
    ]
}

fn qe_workflows() -> Vec<EngineWorkflowDescriptor> {
    vec![
        EngineWorkflowDescriptor {
            kind: CalculationKind::Scf,
            label: "SCF".to_string(),
            route_key: "qe.scf".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: common_structure_panels(),
            engine_panels: vec![
                engine_panel(
                    "qe-scf-calculation-mode",
                    "Calculation mode",
                    "qe.scf.calculationMode",
                    30,
                    true,
                    "QE pw.x calculation mode, including SCF and relax/vc-relax options.",
                ),
                engine_panel(
                    "qe-scf-pseudopotentials",
                    "Pseudopotentials",
                    "qe.scf.pseudopotentials",
                    40,
                    true,
                    "QE pseudopotential selection and SSSP/UPF-derived cutoff hints.",
                ),
                engine_panel(
                    "qe-scf-electronic-controls",
                    "Electronic controls",
                    "qe.scf.electronicControls",
                    50,
                    true,
                    "QE ecutwfc, ecutrho, occupations, smearing, k-points, and convergence controls.",
                ),
                engine_panel(
                    "qe-scf-hubbard-magnetism",
                    "Hubbard and magnetism",
                    "qe.scf.hubbardMagnetism",
                    60,
                    false,
                    "QE DFT+U/HUBBARD syntax and spin initialization controls.",
                ),
            ],
            input_requirements: vec![
                requirement(WorkflowInputRequirementKind::CrystalStructure, "Crystal structure", true),
                requirement(WorkflowInputRequirementKind::CifContent, "CIF content", true),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "QE runtime profile",
                    true,
                ),
            ],
            produces: vec![
                EngineResultDatasetKind::ScfSummary,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::StructureOptimization,
            label: "Structure optimization".to_string(),
            route_key: "qe.scf".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: common_structure_panels(),
            engine_panels: vec![
                engine_panel(
                    "qe-relax-mode",
                    "Relaxation mode",
                    "qe.scf.calculationMode",
                    30,
                    true,
                    "QE relax/vc-relax mode and optimization constraints.",
                ),
                engine_panel(
                    "qe-relax-pseudopotentials",
                    "Pseudopotentials",
                    "qe.scf.pseudopotentials",
                    40,
                    true,
                    "QE pseudopotential selection and cutoff hints for relaxation.",
                ),
                engine_panel(
                    "qe-relax-electronic-controls",
                    "Electronic controls",
                    "qe.scf.electronicControls",
                    50,
                    true,
                    "QE electronic and ionic convergence controls.",
                ),
            ],
            input_requirements: vec![
                requirement(WorkflowInputRequirementKind::CrystalStructure, "Crystal structure", true),
                requirement(WorkflowInputRequirementKind::CifContent, "CIF content", true),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "QE runtime profile",
                    true,
                ),
            ],
            produces: vec![
                EngineResultDatasetKind::ScfSummary,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Bands,
            label: "Bands".to_string(),
            route_key: "qe.bands".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: {
                let mut panels = source_scf_panels();
                panels.insert(
                    1,
                    shared_panel(SharedWorkflowPanelKind::KPathViewer, "Band path", 20, true),
                );
                panels.push(shared_panel(SharedWorkflowPanelKind::Viewer, "Band viewer", 120, true));
                panels
            },
            engine_panels: vec![
                engine_panel(
                    "qe-bands-nscf",
                    "NSCF and bands.x",
                    "qe.bands.nscfBands",
                    30,
                    true,
                    "QE pw.x NSCF and bands.x controls.",
                ),
                engine_panel(
                    "qe-bands-projections",
                    "Projections",
                    "qe.bands.projections",
                    40,
                    false,
                    "Optional QE projwfc.x fat-band projection controls.",
                ),
            ],
            input_requirements: vec![
                requirement(
                    WorkflowInputRequirementKind::SourceScfCalculation,
                    "Source SCF calculation",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "QE runtime profile",
                    true,
                ),
            ],
            produces: vec![
                EngineResultDatasetKind::BandDataset,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Dos,
            label: "DOS".to_string(),
            route_key: "qe.dos".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: source_scf_panels(),
            engine_panels: vec![engine_panel(
                "qe-dos-projwfc",
                "DOS and projections",
                "qe.dos.projwfc",
                30,
                true,
                "QE NSCF, dos.x, and projwfc.x DOS controls.",
            )],
            input_requirements: vec![
                requirement(
                    WorkflowInputRequirementKind::SourceScfCalculation,
                    "Source SCF calculation",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "QE runtime profile",
                    true,
                ),
            ],
            produces: vec![EngineResultDatasetKind::DosDataset],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::FermiSurface,
            label: "Fermi surface".to_string(),
            route_key: "qe.fermiSurface".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: source_scf_panels(),
            engine_panels: vec![engine_panel(
                "qe-fermi-surface",
                "Fermi surface",
                "qe.fermiSurface.fs",
                30,
                true,
                "QE dense NSCF and Fermi-surface artifact generation controls.",
            )],
            input_requirements: vec![requirement(
                WorkflowInputRequirementKind::SourceScfCalculation,
                "Source SCF calculation",
                true,
            )],
            produces: vec![EngineResultDatasetKind::NativeArtifacts],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Phonon,
            label: "Phonons".to_string(),
            route_key: "qe.phonon".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: source_scf_panels(),
            engine_panels: vec![engine_panel(
                "qe-phonon-pipeline",
                "Phonon pipeline",
                "qe.phonon.pipeline",
                30,
                true,
                "QE ph.x, q2r.x, and matdyn.x controls.",
            )],
            input_requirements: vec![requirement(
                WorkflowInputRequirementKind::SourceScfCalculation,
                "Source SCF calculation",
                true,
            )],
            produces: vec![
                EngineResultDatasetKind::PhononBands,
                EngineResultDatasetKind::PhononDos,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::HubbardLrt,
            label: "Hubbard LRT".to_string(),
            route_key: "qe.hubbardLrt".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: source_scf_panels(),
            engine_panels: vec![engine_panel(
                "qe-hubbard-lrt",
                "hp.x linear response",
                "qe.hubbardLrt.hp",
                30,
                true,
                "QE hp.x q-mesh, perturbation, and Hubbard parameter controls.",
            )],
            input_requirements: vec![requirement(
                WorkflowInputRequirementKind::SourceScfCalculation,
                "DFT+U source SCF calculation",
                true,
            )],
            produces: vec![EngineResultDatasetKind::NativeArtifacts],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Wannier,
            label: "Wannier".to_string(),
            route_key: "qe.wannier".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: source_scf_panels(),
            engine_panels: vec![
                engine_panel(
                    "qe-wannier-nscf",
                    "NSCF and interfaces",
                    "qe.wannier.nscfInterfaces",
                    30,
                    true,
                    "QE NSCF, pw2wannier90.x, and wannier90.x setup controls.",
                ),
                engine_panel(
                    "qe-wannier-projections",
                    "Wannier projections",
                    "qe.wannier.projections",
                    40,
                    true,
                    "QE/Wannier90 projection and disentanglement controls.",
                ),
            ],
            input_requirements: vec![requirement(
                WorkflowInputRequirementKind::SourceScfCalculation,
                "Source SCF calculation",
                true,
            )],
            produces: vec![
                EngineResultDatasetKind::WannierResult,
                EngineResultDatasetKind::BandDataset,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Transport,
            label: "Transport".to_string(),
            route_key: "qe.transport".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: vec![
                shared_panel(
                    SharedWorkflowPanelKind::SourceCalculationSelector,
                    "Source Wannier calculation",
                    10,
                    true,
                ),
                shared_panel(
                    SharedWorkflowPanelKind::HpcRunSettings,
                    "Execution settings",
                    90,
                    true,
                ),
                shared_panel(SharedWorkflowPanelKind::LiveOutput, "Live output", 100, true),
                shared_panel(
                    SharedWorkflowPanelKind::Viewer,
                    "Transport viewer",
                    120,
                    true,
                ),
            ],
            engine_panels: vec![engine_panel(
                "qe-transport-boltz",
                "BoltzWann transport",
                "qe.transport.boltzWann",
                30,
                true,
                "QE/Wannier90 BoltzWann and postw90 transport controls.",
            )],
            input_requirements: vec![requirement(
                WorkflowInputRequirementKind::SourceWannierCalculation,
                "Source Wannier calculation",
                true,
            )],
            produces: vec![EngineResultDatasetKind::TransportDataset],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Epw,
            label: "EPW".to_string(),
            route_key: "qe.epw".to_string(),
            execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
            shared_panels: source_scf_panels(),
            engine_panels: vec![engine_panel(
                "qe-epw",
                "EPW",
                "qe.epw.workflow",
                30,
                true,
                "QE epw.x source selection, mesh, transport, and superconductivity controls.",
            )],
            input_requirements: vec![
                requirement(
                    WorkflowInputRequirementKind::SourceScfCalculation,
                    "Source SCF/NSCF calculation",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::SourcePhononCalculation,
                    "Source phonon calculation",
                    false,
                ),
            ],
            produces: vec![
                EngineResultDatasetKind::EpwDataset,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
    ]
}

impl EnginePlugin for QeEnginePlugin {
    fn id(&self) -> EngineId {
        EngineId::Qe
    }

    fn descriptor(&self) -> EngineDescriptor {
        qe_engine_descriptor()
    }

    fn manifest(&self) -> EnginePluginManifest {
        EnginePluginManifest {
            descriptor: self.descriptor(),
            hpc: EngineHpcInterfaceDescriptor {
                supports_local: true,
                supports_hpc: true,
                supports_remote_only: false,
                runtime_profile_key: "qe.hpc.runtimeProfile".to_string(),
            },
            workflows: qe_workflows(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qe_plugin_manifest_exposes_workflows_without_shared_input_schema() {
        let manifest = QE_ENGINE_PLUGIN.manifest();

        assert_eq!(manifest.descriptor.id, EngineId::Qe);
        assert!(manifest.hpc.supports_local);
        assert!(manifest.hpc.supports_hpc);
        assert!(!manifest.hpc.supports_remote_only);
        assert!(manifest
            .workflows
            .iter()
            .any(|workflow| workflow.kind == CalculationKind::Scf
                && workflow.route_key == "qe.scf"
                && workflow
                    .engine_panels
                    .iter()
                    .any(|panel| panel.component_key == "qe.scf.pseudopotentials")));
        assert!(manifest
            .workflows
            .iter()
            .any(|workflow| workflow.kind == CalculationKind::Bands
                && workflow
                    .produces
                    .contains(&EngineResultDatasetKind::BandDataset)));
    }
}
