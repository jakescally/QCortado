//! Engine-neutral metadata types.
//!
//! These types are deliberately narrow. They identify engines and broad result
//! categories for platform metadata, but they are not calculation input
//! schemas. QE-specific inputs remain in [`crate::qe`], and future Wien2k
//! inputs should model Wien2k workflows directly.

use serde::{Deserialize, Serialize};

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
