//! Shared engine metadata helpers.
//!
//! These helpers describe platform-facing engine identity only. They do not
//! replace the current QE implementation and they do not introduce shared
//! calculation input models.

use super::types::{
    CalculationKind, EngineDescriptor, EngineExecutionMode, EngineId, EngineImplementationStatus,
};

/// Existing projects and saved calculations without explicit engine metadata
/// should be treated as Quantum ESPRESSO during the migration.
pub const LEGACY_PROJECT_ENGINE_ID: EngineId = EngineId::Qe;

/// Returns the descriptor for the currently implemented QE engine.
///
/// The descriptor is metadata for future engine selection boundaries. Runtime
/// behavior still flows through [`crate::qe`] and the existing Tauri commands.
pub fn qe_engine_descriptor() -> EngineDescriptor {
    EngineDescriptor {
        id: EngineId::Qe,
        label: "Quantum ESPRESSO".to_string(),
        status: EngineImplementationStatus::Implemented,
        execution_modes: vec![EngineExecutionMode::Local, EngineExecutionMode::Hpc],
        calculation_kinds: vec![
            CalculationKind::Scf,
            CalculationKind::StructureOptimization,
            CalculationKind::Bands,
            CalculationKind::Dos,
            CalculationKind::FermiSurface,
            CalculationKind::Phonon,
            CalculationKind::HubbardLrt,
            CalculationKind::Wannier,
            CalculationKind::Transport,
            CalculationKind::Epw,
        ],
    }
}

/// Engine descriptors that are safe to present as implemented today.
pub fn implemented_engine_descriptors() -> Vec<EngineDescriptor> {
    vec![qe_engine_descriptor()]
}
