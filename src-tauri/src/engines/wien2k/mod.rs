//! WIEN2k engine boundary.
//!
//! Structure-source setup and SCF are exposed WIEN2k workflows. Unlike QE,
//! WIEN2k owns a remote case directory whose files and lifecycle stay behind
//! this engine adapter.

use super::types::EngineId;

pub mod bands;
pub mod case_state;
pub mod commands;
pub mod plugin;
pub mod scf;
pub mod structure;
pub mod types;

pub const ENGINE_ID: EngineId = EngineId::Wien2k;

pub use bands::*;
pub use case_state::{
    core_case_artifacts, initialized_case_artifacts, normalize_case_name, wien2k_case_file,
    Wien2kCaseState,
};
pub use commands::{build_init_lapw_plan, build_scf_run_plan};
pub use plugin::{Wien2kReservedEnginePlugin, WIEN2K_RESERVED_ENGINE_PLUGIN};
pub use scf::*;
pub use structure::*;
pub use types::*;
