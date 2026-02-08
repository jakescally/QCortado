//! Quantum ESPRESSO integration module.
//!
//! This module provides:
//! - Type definitions for QE calculations (`types`)
//! - Input file generation (`input`)
//! - Output parsing (`output`)
//! - Process execution (`runner`)
//! - Band structure support (`bands`)
//! - Phonon calculation support (`phonon`)

pub mod bands;
pub mod input;
pub mod output;
pub mod phonon;
pub mod runner;
pub mod types;

// Re-export commonly used items
pub use bands::{BandData, BandsXConfig, KPathPoint, generate_bands_x_input, parse_bands_gnu, read_bands_gnu_file};
pub use input::{generate_pw_input, generate_ph_input, generate_q2r_input, generate_matdyn_dos_input, generate_matdyn_bands_input};
pub use output::{parse_dos_file, parse_pw_output};
pub use phonon::{parse_ph_output, parse_phonon_dos, parse_phonon_dispersion, add_phonon_symmetry_markers, read_phonon_dos_file, read_phonon_dispersion_file};
pub use runner::{QERunner, RunnerError};
pub use types::*;
