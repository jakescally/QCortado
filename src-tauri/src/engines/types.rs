//! Engine-neutral metadata types.
//!
//! These types are deliberately narrow. They identify engines and broad result
//! categories for platform metadata, but they are not calculation input
//! schemas. QE-specific inputs remain in [`crate::qe`], and future Wien2k
//! inputs should model Wien2k workflows directly.

use serde::{Deserialize, Serialize};

pub type MetadataMap = serde_json::Map<String, serde_json::Value>;
pub type Vector3 = [f64; 3];
pub type Matrix3x3 = [Vector3; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineId {
    /// The existing, implemented Quantum ESPRESSO engine.
    Qe,
    /// Reserved for a future remote-only Wien2k backend.
    ///
    /// This variant is a placeholder identity only; it does not mean Wien2k is
    /// implemented.
    Wien2k,
}

impl EngineId {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Qe => "qe",
            Self::Wien2k => "wien2k",
        }
    }
}

impl Default for EngineId {
    fn default() -> Self {
        Self::Qe
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineImplementationStatus {
    Implemented,
    Reserved,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineExecutionMode {
    Local,
    Hpc,
    /// Reserved for future engines whose execution is managed outside the
    /// local desktop process.
    Remote,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CalculationKind {
    Scf,
    StructureOptimization,
    Bands,
    Dos,
    FermiSurface,
    Phonon,
    HubbardLrt,
    Wannier,
    Transport,
    Epw,
    EngineSetup,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineTaskKind {
    /// A pw.x or engine-native self-consistent-field calculation.
    Scf,
    /// A relax/vc-relax or engine-native structural optimization task.
    StructureOptimization,
    Bands,
    Dos,
    FermiSurface,
    Phonon,
    HubbardLrt,
    Wannier,
    Transport,
    Epw,
    /// Engine setup or initialization work, such as a future Wien2k case setup.
    EngineSetup,
    /// Artifact upload/download/synchronization work.
    ArtifactSync,
    /// Engine installation, profile, or remote environment inspection.
    EnvironmentProbe,
    /// Parser-to-normalized-result conversion work.
    ResultNormalization,
    Other,
}

impl EngineTaskKind {
    pub const fn calculation_kind(self) -> Option<CalculationKind> {
        match self {
            Self::Scf => Some(CalculationKind::Scf),
            Self::StructureOptimization => Some(CalculationKind::StructureOptimization),
            Self::Bands => Some(CalculationKind::Bands),
            Self::Dos => Some(CalculationKind::Dos),
            Self::FermiSurface => Some(CalculationKind::FermiSurface),
            Self::Phonon => Some(CalculationKind::Phonon),
            Self::HubbardLrt => Some(CalculationKind::HubbardLrt),
            Self::Wannier => Some(CalculationKind::Wannier),
            Self::Transport => Some(CalculationKind::Transport),
            Self::Epw => Some(CalculationKind::Epw),
            Self::EngineSetup => Some(CalculationKind::EngineSetup),
            Self::Other => Some(CalculationKind::Other),
            Self::ArtifactSync | Self::EnvironmentProbe | Self::ResultNormalization => None,
        }
    }
}

impl From<CalculationKind> for EngineTaskKind {
    fn from(kind: CalculationKind) -> Self {
        match kind {
            CalculationKind::Scf => Self::Scf,
            CalculationKind::StructureOptimization => Self::StructureOptimization,
            CalculationKind::Bands => Self::Bands,
            CalculationKind::Dos => Self::Dos,
            CalculationKind::FermiSurface => Self::FermiSurface,
            CalculationKind::Phonon => Self::Phonon,
            CalculationKind::HubbardLrt => Self::HubbardLrt,
            CalculationKind::Wannier => Self::Wannier,
            CalculationKind::Transport => Self::Transport,
            CalculationKind::Epw => Self::Epw,
            CalculationKind::EngineSetup => Self::EngineSetup,
            CalculationKind::Other => Self::Other,
        }
    }
}

/// Platform-facing metadata for an engine.
///
/// This is suitable for registry and project metadata work. It intentionally
/// does not contain engine-specific input defaults such as QE pseudopotentials
/// or future Wien2k case setup fields.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineDescriptor {
    pub id: EngineId,
    pub label: String,
    pub status: EngineImplementationStatus,
    pub execution_modes: Vec<EngineExecutionMode>,
    pub calculation_kinds: Vec<CalculationKind>,
}

/// Platform-facing metadata for an engine task.
///
/// This is not an input schema. Engine-native input payloads should stay in the
/// corresponding engine module (`crate::qe` today, future engine modules later).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineTaskMetadata {
    pub engine_id: EngineId,
    pub task_kind: EngineTaskKind,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub calculation_kind: Option<CalculationKind>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub task_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub calculation_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub project_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cif_id: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_task_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_calculation_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub created_at: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<MetadataMap>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArtifactLocationKind {
    LocalPath,
    RemotePath,
    ProjectRelativePath,
    ArchiveRelativePath,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArtifactRole {
    Input,
    Output,
    Log,
    Scratch,
    Viewer,
    Provenance,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ArtifactReference {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub role: ArtifactRole,
    pub location_kind: ArtifactLocationKind,
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub engine_id: Option<EngineId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size_bytes: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub media_type: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checksum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<MetadataMap>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiagnosticSeverity {
    Info,
    Warning,
    Error,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiagnosticSource {
    Platform,
    Engine,
    Parser,
    Viewer,
    Hpc,
    Storage,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DiagnosticHint {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub severity: DiagnosticSeverity,
    pub source: DiagnosticSource,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub code: Option<String>,
    pub title: String,
    pub message: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub suggested_action: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub related_artifact_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<MetadataMap>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NormalizedResultSchema {
    #[serde(rename = "cortado.scf_summary.v1")]
    ScfSummaryV1,
    #[serde(rename = "cortado.band_result_metadata.v1")]
    BandResultMetadataV1,
}

fn default_scf_summary_schema() -> NormalizedResultSchema {
    NormalizedResultSchema::ScfSummaryV1
}

fn default_band_result_metadata_schema() -> NormalizedResultSchema {
    NormalizedResultSchema::BandResultMetadataV1
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineResultProvenance {
    pub engine_id: EngineId,
    pub task_kind: EngineTaskKind,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub calculation_kind: Option<CalculationKind>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub calculation_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub project_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cif_id: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_calculation_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub generated_at: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EnergyUnit {
    #[serde(rename = "eV")]
    Ev,
    #[serde(rename = "Ry")]
    Ry,
    #[serde(rename = "Ha")]
    Ha,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ForceUnit {
    #[serde(rename = "Ry/Bohr")]
    RyPerBohr,
    #[serde(rename = "eV/angstrom")]
    EvPerAngstrom,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum StressUnit {
    #[serde(rename = "kbar")]
    Kbar,
    #[serde(rename = "GPa")]
    Gpa,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EnergyQuantity {
    pub value: f64,
    pub unit: EnergyUnit,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ForceSet {
    pub values: Vec<Vector3>,
    pub unit: ForceUnit,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct StressTensor {
    pub values: Matrix3x3,
    pub unit: StressUnit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ScfConvergenceState {
    Converged,
    NotConverged,
    Failed,
    Cancelled,
    Unknown,
}

/// Engine-neutral SCF result summary for shared project/result metadata.
///
/// This is an output contract only. It does not replace [`crate::qe::QEResult`]
/// in this phase and it should not be used as an engine-neutral input model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NormalizedScfSummary {
    #[serde(default = "default_scf_summary_schema")]
    pub schema: NormalizedResultSchema,
    pub provenance: EngineResultProvenance,
    pub convergence: ScfConvergenceState,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub total_energy: Option<EnergyQuantity>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fermi_energy_ev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scf_steps: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_seconds: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub total_magnetization: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub atomic_magnetic_moments: Option<Vec<Vector3>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub forces: Option<ForceSet>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stress: Option<StressTensor>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub diagnostics: Vec<DiagnosticHint>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub artifacts: Vec<ArtifactReference>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<MetadataMap>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BandPathMarker {
    pub x: f64,
    pub label: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NormalizedBandGap {
    pub value_ev: f64,
    pub is_direct: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vbm_x: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cbm_x: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vbm_energy_ev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cbm_energy_ev: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BandProjectionGroupKind {
    Atom,
    Orbital,
    ElementOrbital,
    Other,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BandProjectionSummary {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub total_group_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub atom_group_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orbital_group_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub element_orbital_group_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub group_kinds: Vec<BandProjectionGroupKind>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<MetadataMap>,
}

/// Engine-neutral band-result metadata.
///
/// This intentionally stores metadata and summary values, not full band-energy
/// arrays. Full engine-native payloads such as QE [`crate::qe::BandData`] stay
/// unchanged until a later adapter phase wires normalized viewer datasets into
/// project storage or commands.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NormalizedBandResultMetadata {
    #[serde(default = "default_band_result_metadata_schema")]
    pub schema: NormalizedResultSchema,
    pub provenance: EngineResultProvenance,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub n_bands: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub n_kpoints: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_energy_ev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub energy_range_ev: Option<[f64; 2]>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub path_markers: Vec<BandPathMarker>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub band_gap: Option<NormalizedBandGap>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub projections: Option<BandProjectionSummary>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub diagnostics: Vec<DiagnosticHint>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub artifacts: Vec<ArtifactReference>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<MetadataMap>,
}
