//! Engine boundary scaffolding for the Cortado platform migration.
//!
//! This module is intentionally additive. The existing [`crate::qe`] module
//! remains the source of truth for all Quantum ESPRESSO behavior, command
//! handlers, parsers, runners, and input generation. New code can use this
//! namespace to describe engine identity and migration boundaries without
//! moving QE files or introducing shared input schemas.
//!
//! Do not add engine-agnostic calculation input structs here. Engine inputs
//! must stay native to each engine; shared platform code should normalize only
//! metadata and result/viewer datasets where that is useful.

pub mod common;
pub mod qe;
pub mod types;

pub use types::*;
