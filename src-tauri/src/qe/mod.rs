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
pub mod epw;
pub mod input;
pub mod output;
pub mod phonon;
pub mod runner;
pub mod transport;
pub mod types;
pub mod wannier;

// Re-export commonly used items
pub use bands::{
    generate_bands_x_input, generate_projwfc_input, parse_bands_gnu,
    parse_projwfc_projection_groups, parse_projwfc_projection_groups_aligned, read_bands_gnu_file,
    BandData, BandProjectionData, BandProjectionGroup, BandsXConfig, KPathPoint, ProjwfcConfig,
};
pub use epw::{
    build_epw_input, build_epw_input_preview, build_epw_keyword_map, collect_epw_artifacts,
    parse_epw_result_summary, render_epw_input, validate_epw_config, EpwArtifactManifestEntry,
    EpwCalculationConfig, EpwCalculationV1, EpwErrorRecord, EpwInputConfig, EpwInputPreviewResult,
    EpwExtensionsV1, EpwPrerequisiteValidation, EpwResultSummaryV1, EpwRuntimeConfig,
    EpwSourceRef, EpwSourcesV1,
    EPW_SCHEMA_VERSION,
};
pub use input::{
    generate_dos_input, generate_matdyn_bands_input, generate_matdyn_dos_input, generate_ph_input,
    generate_pw_input, generate_q2r_input,
};
pub use output::{parse_dos_file, parse_pw_output};
pub use phonon::{
    add_phonon_symmetry_markers, parse_ph_output, parse_phonon_dispersion, parse_phonon_dos,
    read_phonon_dispersion_file, read_phonon_dos_file,
};
pub use runner::{QERunner, RunnerError};
pub use transport::{
    build_transport_win, collect_transport_artifacts, parse_transport_result,
    validate_transport_config, TransportCalculationConfig, TransportDataset, TransportResult,
    TransportTdfData,
};
pub use types::*;
pub use wannier::{
    collect_wannier_artifacts, export_ludwig_bundle, generate_pw2wannier90_input,
    generate_uniform_mp_kpoints, generate_wannier90_win, parse_wannier_band_data,
    parse_wannier_centres_xyz, parse_wannier_hamiltonian, parse_wannier_wout,
    prepare_wannier_nscf_calculation, read_wannier_result, validate_wannier_config,
    LudwigExportConfig, LudwigExportMetadata, LudwigExportMode, LudwigExportResult,
    Pw2Wannier90Config, WannierArtifact, WannierBandPathConfig, WannierCalculationConfig,
    WannierCentreRecord, WannierConvergenceData, WannierDisentanglementConfig, WannierHamiltonian,
    WannierProjectionSpec, WannierProjectionTargetType, WannierResult, WannierSourceMetadata,
    WannierSpread,
};
