//! Compatibility shim for the Quantum ESPRESSO engine.
//!
//! New backend code should import QE functionality through
//! [`crate::engines::qe`]. This module remains so existing saved paths,
//! internal tests, and any downstream imports continue to compile during the
//! engine migration.

pub use crate::engines::qe::*;
