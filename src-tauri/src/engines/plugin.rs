//! Engine plugin contract for platform-to-engine interaction.
//!
//! This is the backend API boundary for computation engines. It describes
//! engine identity, supported workflows, and frontend panel ownership without
//! defining shared calculation input schemas.

use serde::{Deserialize, Serialize};

use super::types::{CalculationKind, EngineDescriptor, EngineExecutionMode, EngineId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SharedWorkflowPanelKind {
    StructureSource,
    StructureViewer,
    SourceCalculationSelector,
    KPathViewer,
    HpcRunSettings,
    LiveOutput,
    ResultSummary,
    Viewer,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WorkflowInputRequirementKind {
    CrystalStructure,
    CifContent,
    SourceScfCalculation,
    SourceBandsCalculation,
    SourceWannierCalculation,
    SourcePhononCalculation,
    EngineRuntimeProfile,
    /// Engine-owned persisted run state, such as a future Wien2k remote case
    /// directory. This is metadata about where native state lives, not a
    /// shared calculation input schema.
    EngineCaseState,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineResultDatasetKind {
    ScfSummary,
    BandDataset,
    DosDataset,
    PhononBands,
    PhononDos,
    WannierResult,
    TransportDataset,
    EpwDataset,
    NativeArtifacts,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SharedWorkflowPanelDescriptor {
    pub kind: SharedWorkflowPanelKind,
    pub label: String,
    pub order: u16,
    #[serde(default)]
    pub required: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineWorkflowPanelDescriptor {
    pub id: String,
    pub label: String,
    pub component_key: String,
    pub order: u16,
    #[serde(default)]
    pub required: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WorkflowInputRequirement {
    pub kind: WorkflowInputRequirementKind,
    pub label: String,
    #[serde(default)]
    pub required: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineWorkflowDescriptor {
    pub kind: CalculationKind,
    pub label: String,
    /// Stable frontend route key. The frontend plugin registry maps this to
    /// concrete wizard/viewer components for the selected engine.
    pub route_key: String,
    pub execution_modes: Vec<EngineExecutionMode>,
    pub shared_panels: Vec<SharedWorkflowPanelDescriptor>,
    pub engine_panels: Vec<EngineWorkflowPanelDescriptor>,
    pub input_requirements: Vec<WorkflowInputRequirement>,
    pub produces: Vec<EngineResultDatasetKind>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineHpcInterfaceDescriptor {
    pub supports_local: bool,
    pub supports_hpc: bool,
    pub supports_remote_only: bool,
    pub runtime_profile_key: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EnginePluginManifest {
    pub descriptor: EngineDescriptor,
    pub hpc: EngineHpcInterfaceDescriptor,
    pub workflows: Vec<EngineWorkflowDescriptor>,
}

pub trait EnginePlugin: Send + Sync {
    fn id(&self) -> EngineId;
    fn descriptor(&self) -> EngineDescriptor;
    fn manifest(&self) -> EnginePluginManifest;

    fn workflow(&self, kind: CalculationKind) -> Option<EngineWorkflowDescriptor> {
        self.manifest()
            .workflows
            .into_iter()
            .find(|workflow| workflow.kind == kind)
    }

    fn supports_workflow(&self, kind: CalculationKind) -> bool {
        self.workflow(kind).is_some()
    }
}

pub fn implemented_engine_plugins() -> Vec<&'static dyn EnginePlugin> {
    vec![&crate::engines::qe::plugin::QE_ENGINE_PLUGIN]
}

pub fn implemented_engine_manifests() -> Vec<EnginePluginManifest> {
    implemented_engine_plugins()
        .into_iter()
        .map(EnginePlugin::manifest)
        .collect()
}

pub fn get_implemented_engine_plugin(engine_id: EngineId) -> Option<&'static dyn EnginePlugin> {
    implemented_engine_plugins()
        .into_iter()
        .find(|plugin| plugin.id() == engine_id)
}

pub fn get_implemented_engine_manifest(engine_id: EngineId) -> Option<EnginePluginManifest> {
    get_implemented_engine_plugin(engine_id).map(EnginePlugin::manifest)
}
