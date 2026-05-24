//! Engine boundary scaffolding for the Cortado platform migration.
//!
//! Quantum ESPRESSO lives under [`crate::engines::qe`] as the first explicit
//! engine implementation. The legacy [`crate::qe`] module is a compatibility
//! shim that re-exports this implementation.
//!
//! Do not add engine-agnostic calculation input structs here. Engine inputs
//! must stay native to each engine; shared platform code should normalize only
//! metadata and result/viewer datasets where that is useful.

pub mod common;
pub mod plugin;
pub mod qe;
pub mod types;
pub mod wien2k;

pub use plugin::*;
pub use types::*;
