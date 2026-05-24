//! Hidden WIEN2k plugin skeleton.
//!
//! This manifest is deliberately not registered in `implemented_engine_plugins`.
//! It captures the future platform-to-engine contract while keeping the UI and
//! command API unchanged.

use crate::engines::plugin::{
    EngineHpcInterfaceDescriptor, EnginePlugin, EnginePluginManifest, EngineResultDatasetKind,
    EngineWorkflowDescriptor, EngineWorkflowPanelDescriptor, SharedWorkflowPanelDescriptor,
    SharedWorkflowPanelKind, WorkflowInputRequirement, WorkflowInputRequirementKind,
};
use crate::engines::types::{
    CalculationKind, EngineDescriptor, EngineExecutionMode, EngineId, EngineImplementationStatus,
};

pub struct Wien2kReservedEnginePlugin;

pub static WIEN2K_RESERVED_ENGINE_PLUGIN: Wien2kReservedEnginePlugin = Wien2kReservedEnginePlugin;

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

fn structure_panels() -> Vec<SharedWorkflowPanelDescriptor> {
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
            "Remote execution settings",
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
            "Source Wien2k SCF case",
            10,
            true,
        ),
        shared_panel(
            SharedWorkflowPanelKind::HpcRunSettings,
            "Remote execution settings",
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

fn workflows() -> Vec<EngineWorkflowDescriptor> {
    vec![
        EngineWorkflowDescriptor {
            kind: CalculationKind::EngineSetup,
            label: "Wien2k case setup".to_string(),
            route_key: "wien2k.caseSetup".to_string(),
            execution_modes: vec![EngineExecutionMode::Remote],
            shared_panels: structure_panels(),
            engine_panels: vec![
                engine_panel(
                    "wien2k-case-directory",
                    "Case directory",
                    "wien2k.case.directory",
                    30,
                    true,
                    "Remote WIEN2k case prefix and case-directory staging controls.",
                ),
                engine_panel(
                    "wien2k-struct-generation",
                    "case.struct",
                    "wien2k.case.struct",
                    40,
                    true,
                    "CIF-to-case.struct generation and validation owned by the WIEN2k engine.",
                ),
                engine_panel(
                    "wien2k-init-lapw",
                    "init_lapw",
                    "wien2k.init.controls",
                    50,
                    true,
                    "WIEN2k initialization controls for RMT, RKmax, Gmax, lmax, lstart, kgen, and dstart.",
                ),
            ],
            input_requirements: vec![
                requirement(WorkflowInputRequirementKind::CrystalStructure, "Crystal structure", true),
                requirement(WorkflowInputRequirementKind::CifContent, "CIF content", true),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "Wien2k remote runtime profile",
                    true,
                ),
            ],
            produces: vec![EngineResultDatasetKind::NativeArtifacts],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Scf,
            label: "Wien2k SCF".to_string(),
            route_key: "wien2k.scf".to_string(),
            execution_modes: vec![EngineExecutionMode::Remote],
            shared_panels: structure_panels(),
            engine_panels: vec![
                engine_panel(
                    "wien2k-case-state",
                    "Case state",
                    "wien2k.case.state",
                    30,
                    true,
                    "Remote case-state selection, staging, and resume policy for WIEN2k.",
                ),
                engine_panel(
                    "wien2k-scf-basis",
                    "Basis and spheres",
                    "wien2k.scf.basis",
                    40,
                    true,
                    "WIEN2k RMT, RKmax, Gmax, and lmax controls. No pseudopotentials are used.",
                ),
                engine_panel(
                    "wien2k-scf-cycle",
                    "SCF cycle",
                    "wien2k.scf.cycle",
                    50,
                    true,
                    "WIEN2k run_lapw/runsp_lapw convergence, spin mode, and parallel flags.",
                ),
            ],
            input_requirements: vec![
                requirement(WorkflowInputRequirementKind::CrystalStructure, "Crystal structure", true),
                requirement(WorkflowInputRequirementKind::CifContent, "CIF content", true),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "Wien2k remote runtime profile",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::EngineCaseState,
                    "Wien2k remote case state",
                    false,
                ),
            ],
            produces: vec![
                EngineResultDatasetKind::ScfSummary,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
        EngineWorkflowDescriptor {
            kind: CalculationKind::Bands,
            label: "Wien2k bands".to_string(),
            route_key: "wien2k.bands".to_string(),
            execution_modes: vec![EngineExecutionMode::Remote],
            shared_panels: {
                let mut panels = source_scf_panels();
                panels.insert(
                    1,
                    shared_panel(SharedWorkflowPanelKind::KPathViewer, "Band path", 20, true),
                );
                panels.push(shared_panel(SharedWorkflowPanelKind::Viewer, "Band viewer", 120, true));
                panels
            },
            engine_panels: vec![engine_panel(
                "wien2k-bands-spaghetti",
                "lapw1/lapw2/spaghetti",
                "wien2k.bands.spaghetti",
                30,
                true,
                "WIEN2k band command sequence, including case.klist_band and spaghetti output parsing.",
            )],
            input_requirements: vec![
                requirement(
                    WorkflowInputRequirementKind::SourceScfCalculation,
                    "Source Wien2k SCF calculation",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::EngineCaseState,
                    "Converged Wien2k case state",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::EngineRuntimeProfile,
                    "Wien2k remote runtime profile",
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
            label: "Wien2k DOS".to_string(),
            route_key: "wien2k.dos".to_string(),
            execution_modes: vec![EngineExecutionMode::Remote],
            shared_panels: source_scf_panels(),
            engine_panels: vec![engine_panel(
                "wien2k-dos",
                "DOS sequence",
                "wien2k.dos.sequence",
                30,
                true,
                "Future WIEN2k DOS command sequence and DOS file parser boundary.",
            )],
            input_requirements: vec![
                requirement(
                    WorkflowInputRequirementKind::SourceScfCalculation,
                    "Source Wien2k SCF calculation",
                    true,
                ),
                requirement(
                    WorkflowInputRequirementKind::EngineCaseState,
                    "Converged Wien2k case state",
                    true,
                ),
            ],
            produces: vec![
                EngineResultDatasetKind::DosDataset,
                EngineResultDatasetKind::NativeArtifacts,
            ],
        },
    ]
}

impl EnginePlugin for Wien2kReservedEnginePlugin {
    fn id(&self) -> EngineId {
        EngineId::Wien2k
    }

    fn descriptor(&self) -> EngineDescriptor {
        EngineDescriptor {
            id: EngineId::Wien2k,
            label: "WIEN2k".to_string(),
            status: EngineImplementationStatus::Reserved,
            execution_modes: vec![EngineExecutionMode::Remote],
            calculation_kinds: vec![
                CalculationKind::EngineSetup,
                CalculationKind::Scf,
                CalculationKind::Bands,
                CalculationKind::Dos,
            ],
        }
    }

    fn manifest(&self) -> EnginePluginManifest {
        EnginePluginManifest {
            descriptor: self.descriptor(),
            hpc: EngineHpcInterfaceDescriptor {
                supports_local: false,
                supports_hpc: true,
                supports_remote_only: true,
                runtime_profile_key: "wien2k.remote.caseRuntimeProfile".to_string(),
            },
            workflows: workflows(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engines::plugin::{get_implemented_engine_plugin, EnginePlugin};

    #[test]
    fn reserved_wien2k_plugin_is_hidden_from_implemented_registry() {
        assert!(get_implemented_engine_plugin(EngineId::Wien2k).is_none());

        let manifest = WIEN2K_RESERVED_ENGINE_PLUGIN.manifest();
        assert_eq!(manifest.descriptor.id, EngineId::Wien2k);
        assert_eq!(
            manifest.descriptor.status,
            EngineImplementationStatus::Reserved
        );
        assert!(manifest.hpc.supports_remote_only);
        assert!(!manifest.hpc.supports_local);
    }

    #[test]
    fn reserved_manifest_models_case_state_without_pseudopotentials() {
        let manifest = WIEN2K_RESERVED_ENGINE_PLUGIN.manifest();
        let scf = manifest
            .workflows
            .iter()
            .find(|workflow| workflow.kind == CalculationKind::Scf)
            .expect("SCF workflow should exist in reserved skeleton");

        assert!(scf.input_requirements.iter().any(|requirement| {
            requirement.kind == WorkflowInputRequirementKind::EngineCaseState
        }));
        assert!(scf
            .engine_panels
            .iter()
            .all(|panel| !panel.component_key.contains("pseudo")));
        assert!(scf.produces.contains(&EngineResultDatasetKind::ScfSummary));
    }
}
