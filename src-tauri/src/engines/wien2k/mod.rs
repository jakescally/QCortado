//! Hidden WIEN2k engine skeleton.
//!
//! WIEN2k is reserved but not implemented. The skeleton exists so future work
//! can build against a real engine boundary without changing QE behavior or
//! exposing WIEN2k in the UI. Unlike QE, WIEN2k owns a remote case directory
//! whose files and lifecycle must stay behind the engine adapter.

use super::types::EngineId;

pub mod case_state;
pub mod commands;
pub mod plugin;
pub mod types;

pub const ENGINE_ID: EngineId = EngineId::Wien2k;

pub use case_state::{
    core_case_artifacts, initialized_case_artifacts, normalize_case_name, wien2k_case_file,
    Wien2kCaseState,
};
pub use commands::{build_init_lapw_plan, build_scf_run_plan};
pub use plugin::{Wien2kReservedEnginePlugin, WIEN2K_RESERVED_ENGINE_PLUGIN};
pub use types::*;
