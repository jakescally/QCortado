//! Quantum ESPRESSO integration module.
//!
//! This module provides:
//! - Type definitions for QE calculations (`types`)
//! - Input file generation (`input`)
//! - Output parsing (`output`)
//! - Process execution (`runner`)

pub mod input;
pub mod output;
pub mod runner;
pub mod types;

// Re-export commonly used items
pub use input::generate_pw_input;
pub use output::{parse_dos_file, parse_pw_output};
pub use runner::{QERunner, RunnerError};
pub use types::*;
