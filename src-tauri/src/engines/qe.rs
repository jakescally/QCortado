//! Quantum ESPRESSO engine facade.
//!
//! This module exists to give future engine-aware backend code a stable
//! namespace while the existing [`crate::qe`] module remains the implementation
//! source of truth. Re-exporting selected QE items here must not change command
//! behavior or parser behavior.

pub use crate::qe::{
    generate_bands_x_input, generate_dos_input, generate_hp_input, generate_matdyn_bands_input,
    generate_matdyn_dos_input, generate_ph_input, generate_projwfc_input, generate_pw_input,
    generate_q2r_input, parse_dos_file, parse_ph_output, parse_projwfc_projection_groups_aligned,
    parse_pw_output, read_bands_gnu_file, BandData, BandsXConfig, DosCalculation, HpCalculation,
    MatdynCalculation, PhononCalculation, ProjwfcConfig, Q2RCalculation, QECalculation, QEResult,
    QERunner,
};

use super::types::EngineId;

/// Engine identity for the existing QE implementation.
pub const ENGINE_ID: EngineId = EngineId::Qe;
