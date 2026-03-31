//! Wannier90 input generation and result parsing helpers.

use super::bands::{add_symmetry_markers, parse_bands_gnu, BandData, KPathPoint};
use super::types::{
    CalculationType, CellMatrix, KPoint, KPoints, Occupations, PositionUnits, QECalculation,
    SmearingType, StartingPotential,
};
use nalgebra::DMatrix;
use num_complex::Complex64;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Write;
use std::path::{Path, PathBuf};

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;
const TWO_PI: f64 = std::f64::consts::PI * 2.0;

fn default_seedname() -> String {
    "qcortado_wannier".to_string()
}

fn default_guiding_centres() -> bool {
    true
}

fn default_use_ws_distance() -> bool {
    true
}

fn default_write_hr() -> bool {
    true
}

fn default_write_xyz() -> bool {
    true
}

fn default_bands_plot() -> bool {
    true
}

fn default_bands_num_points() -> u32 {
    100
}

fn default_conv_window() -> u32 {
    5
}

fn default_conv_tol() -> f64 {
    1.0e-10
}

fn default_num_iter() -> u32 {
    100
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WannierProjectionTargetType {
    Element,
    Site,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierProjectionSpec {
    pub target_type: WannierProjectionTargetType,
    #[serde(default)]
    pub symbol: Option<String>,
    pub orbital: String,
    #[serde(default)]
    pub site_index: Option<usize>,
    #[serde(default)]
    pub fractional_position: Option<[f64; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WannierDisentanglementConfig {
    #[serde(default)]
    pub exclude_bands: Vec<u32>,
    #[serde(default)]
    pub dis_win_min: Option<f64>,
    #[serde(default)]
    pub dis_win_max: Option<f64>,
    #[serde(default)]
    pub dis_froz_min: Option<f64>,
    #[serde(default)]
    pub dis_froz_max: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierBandPathConfig {
    pub k_path: Vec<KPathPoint>,
    #[serde(default = "default_bands_num_points")]
    pub bands_num_points: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pw2Wannier90Config {
    #[serde(default = "default_true")]
    pub write_amn: bool,
    #[serde(default = "default_true")]
    pub write_mmn: bool,
    #[serde(default)]
    pub write_spn: bool,
    #[serde(default)]
    pub write_unk: bool,
    #[serde(default)]
    pub write_dmn: bool,
    #[serde(default = "default_spin_component")]
    pub spin_component: String,
    #[serde(default)]
    pub scdm_proj: bool,
    #[serde(default)]
    pub atom_proj: bool,
}

fn default_true() -> bool {
    true
}

fn default_spin_component() -> String {
    "none".to_string()
}

impl Default for Pw2Wannier90Config {
    fn default() -> Self {
        Self {
            write_amn: true,
            write_mmn: true,
            write_spn: false,
            write_unk: false,
            write_dmn: false,
            spin_component: default_spin_component(),
            scdm_proj: false,
            atom_proj: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WannierSourceMetadata {
    #[serde(default)]
    pub cell_representation: Option<String>,
    #[serde(default)]
    pub electron_count: Option<f64>,
    #[serde(default)]
    pub nspin: Option<u32>,
    #[serde(default)]
    pub noncolin: Option<bool>,
    #[serde(default)]
    pub lspinorb: Option<bool>,
    #[serde(default)]
    pub lda_plus_u: Option<bool>,
    #[serde(default)]
    pub vdw_corr: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierCalculationConfig {
    pub base_calculation: QECalculation,
    pub k_grid: [u32; 3],
    pub num_wann: u32,
    pub num_bands: u32,
    #[serde(default = "default_seedname")]
    pub seedname: String,
    pub projections: Vec<WannierProjectionSpec>,
    pub band_path: WannierBandPathConfig,
    #[serde(default)]
    pub disentanglement: Option<WannierDisentanglementConfig>,
    #[serde(default)]
    pub pw2wannier90: Option<Pw2Wannier90Config>,
    #[serde(default)]
    pub project_id: Option<String>,
    #[serde(default)]
    pub scf_calc_id: Option<String>,
    #[serde(default)]
    pub source_metadata: Option<WannierSourceMetadata>,
    #[serde(default = "default_guiding_centres")]
    pub guiding_centres: bool,
    #[serde(default = "default_use_ws_distance")]
    pub use_ws_distance: bool,
    #[serde(default = "default_write_hr")]
    pub write_hr: bool,
    #[serde(default = "default_write_xyz")]
    pub write_xyz: bool,
    #[serde(default = "default_bands_plot")]
    pub bands_plot: bool,
    #[serde(default = "default_conv_window")]
    pub conv_window: u32,
    #[serde(default = "default_conv_tol")]
    pub conv_tol: f64,
    #[serde(default = "default_num_iter")]
    pub num_iter: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierSpread {
    pub index: u32,
    pub centre: [f64; 3],
    pub spread: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WannierConvergenceData {
    #[serde(default)]
    pub converged: bool,
    #[serde(default)]
    pub iterations: Option<u32>,
    #[serde(default)]
    pub omega_i: Option<f64>,
    #[serde(default)]
    pub omega_d: Option<f64>,
    #[serde(default)]
    pub omega_od: Option<f64>,
    #[serde(default)]
    pub omega_total: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierArtifact {
    pub file_name: String,
    pub size_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierCentreRecord {
    pub label: String,
    pub position: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierResult {
    pub seedname: String,
    pub num_wann: u32,
    pub num_bands: u32,
    pub k_grid: [u32; 3],
    pub band_data: BandData,
    #[serde(default)]
    pub total_spread: Option<f64>,
    #[serde(default)]
    pub spreads: Vec<WannierSpread>,
    #[serde(default)]
    pub centres: Vec<WannierCentreRecord>,
    #[serde(default)]
    pub convergence: WannierConvergenceData,
    #[serde(default)]
    pub artifact_manifest: Vec<WannierArtifact>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierHrTerm {
    pub r: [i32; 3],
    pub m: usize,
    pub n: usize,
    pub degeneracy: u32,
    pub real: f64,
    pub imag: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WannierHamiltonian {
    pub seedname: String,
    pub num_wann: u32,
    pub lattice_vectors_angstrom: CellMatrix,
    pub hr_terms: Vec<WannierHrTerm>,
    #[serde(default)]
    pub use_ws_distance: bool,
    #[serde(default)]
    pub ws_translations: HashMap<String, Vec<[i32; 3]>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum LudwigExportMode {
    Strict2d,
    Quasi2dFixedSlice,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LudwigExportConfig {
    pub project_id: String,
    pub calculation_id: String,
    pub destination_root: String,
    pub mode: LudwigExportMode,
    pub in_plane_axes: [u8; 2],
    pub slice_axis: u8,
    pub slice_coordinate: f64,
    pub nkx: u32,
    pub nky: u32,
    #[serde(default)]
    pub chemical_potential_ev: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LudwigExportMetadata {
    pub project_id: String,
    pub calculation_id: String,
    pub seedname: String,
    pub mode: LudwigExportMode,
    pub in_plane_axes: [u8; 2],
    pub slice_axis: u8,
    pub slice_coordinate: f64,
    pub chemical_potential_ev: f64,
    pub grid_shape: [u32; 2],
    pub band_count: u32,
    pub lattice_vectors_angstrom: CellMatrix,
    pub in_plane_lattice_angstrom: [[f64; 2]; 2],
    pub fractional_domain: String,
    pub energy_reference: String,
    pub provenance_files: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LudwigExportResult {
    pub bundle_path: String,
    pub band_count: u32,
    pub chemical_potential_ev: f64,
    pub grid_shape: [u32; 2],
}

fn orbital_template_count(orbital: &str) -> Option<u32> {
    match orbital.trim().to_ascii_lowercase().as_str() {
        "s" => Some(1),
        "p" => Some(3),
        "d" => Some(5),
        "f" => Some(7),
        "sp" => Some(4),
        "sp2" => Some(3),
        "sp3" => Some(4),
        "sp3d" => Some(5),
        "sp3d2" => Some(6),
        _ => None,
    }
}

fn parse_f64_token(token: &str, context: &str) -> Result<f64, String> {
    token
        .parse::<f64>()
        .map_err(|_| format!("Failed to parse {} from '{}'.", context, token))
}

fn parse_i32_token(token: &str, context: &str) -> Result<i32, String> {
    token
        .parse::<i32>()
        .map_err(|_| format!("Failed to parse {} from '{}'.", context, token))
}

fn parse_usize_token(token: &str, context: &str) -> Result<usize, String> {
    token
        .parse::<usize>()
        .map_err(|_| format!("Failed to parse {} from '{}'.", context, token))
}

fn convert_cell_matrix_to_angstrom(system: &super::types::QESystem) -> Result<CellMatrix, String> {
    let cell = system
        .cell_parameters
        .ok_or_else(|| "Wannier export requires explicit cell_parameters.".to_string())?;
    let cell_units = system.cell_units.unwrap_or(PositionUnits::Angstrom);

    let scale = match cell_units {
        PositionUnits::Angstrom => 1.0,
        PositionUnits::Bohr => BOHR_TO_ANGSTROM,
        PositionUnits::Alat => {
            let celldm = system.celldm.ok_or_else(|| {
                "CELL_PARAMETERS in alat units require celldm(1) to convert to angstrom."
                    .to_string()
            })?;
            let alat_bohr = celldm[0];
            if alat_bohr <= 0.0 {
                return Err("Invalid celldm(1) for alat cell conversion.".to_string());
            }
            alat_bohr * BOHR_TO_ANGSTROM
        }
        PositionUnits::Crystal => {
            return Err(
                "CELL_PARAMETERS in crystal units cannot be exported to Wannier90 unit_cell_cart."
                    .to_string(),
            );
        }
    };

    let mut converted = cell;
    for row in &mut converted {
        for value in row {
            *value *= scale;
        }
    }
    Ok(converted)
}

fn wsvec_key(r: [i32; 3], m: usize, n: usize) -> String {
    format!(
        "{}:{}:{}:{}:{}",
        r[0],
        r[1],
        r[2],
        m.saturating_add(1),
        n.saturating_add(1)
    )
}

fn parse_win_bool_assignment(content: &str, key: &str) -> Option<bool> {
    let prefix = format!("{} =", key);
    content.lines().find_map(|line| {
        let trimmed = line.trim();
        if !trimmed
            .to_ascii_lowercase()
            .starts_with(&prefix.to_ascii_lowercase())
        {
            return None;
        }
        let (_, raw_value) = trimmed.split_once('=')?;
        let normalized = raw_value
            .trim()
            .trim_matches('.')
            .trim_matches(',')
            .to_ascii_lowercase();
        match normalized.as_str() {
            "true" | "t" => Some(true),
            "false" | "f" => Some(false),
            _ => None,
        }
    })
}

fn parse_win_unit_cell_angstrom(content: &str) -> Result<CellMatrix, String> {
    let mut lines = content.lines().enumerate();
    while let Some((_, line)) = lines.next() {
        if line.trim().eq_ignore_ascii_case("begin unit_cell_cart") {
            let (_, unit_line) = lines
                .next()
                .ok_or_else(|| "Missing unit line in begin unit_cell_cart block.".to_string())?;
            let unit = unit_line.trim().to_ascii_lowercase();
            let scale = match unit.as_str() {
                "ang" | "angstrom" => 1.0,
                "bohr" => BOHR_TO_ANGSTROM,
                "alat" => return Err(
                    "Wannier .win files written in alat units are not supported for Ludwig export."
                        .to_string(),
                ),
                _ => {
                    return Err(format!(
                        "Unsupported unit_cell_cart unit '{}' in Wannier .win file.",
                        unit_line.trim()
                    ))
                }
            };

            let mut cell = [[0.0; 3]; 3];
            for row in &mut cell {
                let (_, row_line) = lines.next().ok_or_else(|| {
                    "Incomplete unit_cell_cart block in Wannier .win file.".to_string()
                })?;
                let parts: Vec<&str> = row_line.split_whitespace().collect();
                if parts.len() < 3 {
                    return Err("Malformed unit_cell_cart row in Wannier .win file.".to_string());
                }
                row[0] = parse_f64_token(parts[0], "unit cell x")? * scale;
                row[1] = parse_f64_token(parts[1], "unit cell y")? * scale;
                row[2] = parse_f64_token(parts[2], "unit cell z")? * scale;
            }
            return Ok(cell);
        }
    }

    Err("Could not find begin unit_cell_cart block in Wannier .win file.".to_string())
}

fn embed_in_plane_lattice(
    lattice_vectors_angstrom: &CellMatrix,
    in_plane_axes: [u8; 2],
) -> Result<[[f64; 2]; 2], String> {
    let axis0 = usize::from(in_plane_axes[0]);
    let axis1 = usize::from(in_plane_axes[1]);
    if axis0 >= 3 || axis1 >= 3 || axis0 == axis1 {
        return Err("In-plane lattice axes must be two distinct values from 0, 1, 2.".to_string());
    }

    let a1 = lattice_vectors_angstrom[axis0];
    let a2 = lattice_vectors_angstrom[axis1];

    let norm_a1 = (a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2]).sqrt();
    if norm_a1 <= 1.0e-12 {
        return Err("Selected first in-plane lattice vector has zero length.".to_string());
    }
    let e1 = [a1[0] / norm_a1, a1[1] / norm_a1, a1[2] / norm_a1];
    let dot_a2_e1 = a2[0] * e1[0] + a2[1] * e1[1] + a2[2] * e1[2];
    let v2 = [
        a2[0] - dot_a2_e1 * e1[0],
        a2[1] - dot_a2_e1 * e1[1],
        a2[2] - dot_a2_e1 * e1[2],
    ];
    let norm_v2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
    if norm_v2 <= 1.0e-12 {
        return Err("Selected in-plane lattice vectors are collinear.".to_string());
    }

    Ok([[norm_a1, 0.0], [dot_a2_e1, norm_v2]])
}

pub fn validate_wannier_config(config: &WannierCalculationConfig) -> Result<(), String> {
    if config.k_grid.iter().any(|value| *value == 0) {
        return Err("Wannier k-grid entries must all be positive.".to_string());
    }
    if config.num_wann == 0 {
        return Err("num_wann must be positive.".to_string());
    }
    if config.num_bands < config.num_wann {
        return Err("num_bands must be greater than or equal to num_wann.".to_string());
    }
    if config.projections.is_empty() {
        return Err("At least one Wannier projection is required.".to_string());
    }
    if config.band_path.k_path.len() < 2 {
        return Err("Wannier interpolation path requires at least two k-points.".to_string());
    }
    if config.band_path.bands_num_points == 0 {
        return Err("bands_num_points must be positive.".to_string());
    }

    let projected_wann = config.projections.iter().try_fold(0u32, |sum, spec| {
        let count = orbital_template_count(&spec.orbital).ok_or_else(|| {
            format!("Unsupported Wannier projection orbital template: {}", spec.orbital)
        })?;
        match spec.target_type {
            WannierProjectionTargetType::Element => {
                let symbol = spec.symbol.as_deref().unwrap_or("").trim();
                if symbol.is_empty() {
                    return Err("Element-targeted Wannier projections require a symbol.".to_string());
                }
                let atom_matches = config
                    .base_calculation
                    .system
                    .atoms
                    .iter()
                    .filter(|atom| atom.symbol.eq_ignore_ascii_case(symbol))
                    .count() as u32;
                if atom_matches == 0 {
                    return Err(format!(
                        "Element-targeted Wannier projection '{}' does not match any atom in the Wannier cell.",
                        symbol
                    ));
                }
                return Ok(sum.saturating_add(count.saturating_mul(atom_matches)));
            }
            WannierProjectionTargetType::Site => {
                if spec.fractional_position.is_none() {
                    return Err(
                        "Site-targeted Wannier projections require fractional_position."
                            .to_string(),
                    );
                }
            }
        }
        Ok(sum.saturating_add(count))
    })?;
    if projected_wann != config.num_wann {
        return Err(format!(
            "Projection set expands to {} Wannier functions, but num_wann = {}.",
            projected_wann, config.num_wann
        ));
    }

    let source = config.source_metadata.as_ref();
    let source_nspin = source
        .and_then(|value| value.nspin)
        .unwrap_or(config.base_calculation.system.nspin);
    let source_noncolin = source
        .and_then(|value| value.noncolin)
        .unwrap_or(config.base_calculation.system.noncolin);
    let source_lspinorb = source
        .and_then(|value| value.lspinorb)
        .unwrap_or(config.base_calculation.system.lspinorb);
    if source_nspin != 1 || source_noncolin || source_lspinorb {
        return Err(
            "Scalar Wannier v1 requires a source SCF with nspin = 1, noncolin = false, and lspinorb = false."
                .to_string(),
        );
    }
    if let Some(electron_count) = source.and_then(|value| value.electron_count) {
        if electron_count > 0.0 {
            let min_bands = (electron_count / 2.0).ceil() as u32;
            if config.num_bands < min_bands {
                return Err(format!(
                    "num_bands = {} is too small for this scalar source: {} electrons require at least {} occupied bands before adding empty/disentanglement states.",
                    config.num_bands, electron_count, min_bands
                ));
            }
        }
    }
    if source
        .and_then(|value| value.cell_representation.as_ref())
        .map(|value| value.to_ascii_lowercase())
        .filter(|value| !value.starts_with("primitive"))
        .is_some()
    {
        return Err(
            "Scalar Wannier v1 requires a source SCF saved in a primitive-cell representation."
                .to_string(),
        );
    }
    if source.and_then(|value| value.lda_plus_u).unwrap_or(false) {
        return Err("Scalar Wannier v1 does not support DFT+U source calculations.".to_string());
    }
    if source
        .and_then(|value| value.vdw_corr.as_ref())
        .map(|value| {
            let trimmed = value.trim();
            !trimmed.is_empty() && !trimmed.eq_ignore_ascii_case("none")
        })
        .unwrap_or(false)
    {
        return Err(
            "Scalar Wannier v1 does not certify vdW-corrected source calculations.".to_string(),
        );
    }

    if let Some(disentanglement) = config.disentanglement.as_ref() {
        if config.num_bands == config.num_wann {
            if !disentanglement.exclude_bands.is_empty()
                || disentanglement.dis_win_min.is_some()
                || disentanglement.dis_win_max.is_some()
                || disentanglement.dis_froz_min.is_some()
                || disentanglement.dis_froz_max.is_some()
            {
                return Err(
                    "Disentanglement controls are only valid when num_bands > num_wann."
                        .to_string(),
                );
            }
        }

        if let (Some(win_min), Some(win_max)) =
            (disentanglement.dis_win_min, disentanglement.dis_win_max)
        {
            if win_min >= win_max {
                return Err("dis_win_min must be smaller than dis_win_max.".to_string());
            }
        }
        if let (Some(froz_min), Some(froz_max)) =
            (disentanglement.dis_froz_min, disentanglement.dis_froz_max)
        {
            if froz_min >= froz_max {
                return Err("dis_froz_min must be smaller than dis_froz_max.".to_string());
            }
        }
        if let (Some(win_min), Some(froz_min)) =
            (disentanglement.dis_win_min, disentanglement.dis_froz_min)
        {
            if froz_min < win_min {
                return Err(
                    "Frozen window must lie inside the outer disentanglement window.".to_string(),
                );
            }
        }
        if let (Some(win_max), Some(froz_max)) =
            (disentanglement.dis_win_max, disentanglement.dis_froz_max)
        {
            if froz_max > win_max {
                return Err(
                    "Frozen window must lie inside the outer disentanglement window.".to_string(),
                );
            }
        }
    }

    Ok(())
}

pub fn generate_uniform_mp_kpoints(grid: [u32; 3]) -> Vec<KPoint> {
    let total = (grid[0] as f64) * (grid[1] as f64) * (grid[2] as f64);
    let weight = if total > 0.0 { 1.0 / total } else { 1.0 };
    let mut points = Vec::with_capacity((grid[0] * grid[1] * grid[2]) as usize);
    for i in 0..grid[0] {
        for j in 0..grid[1] {
            for k in 0..grid[2] {
                points.push(KPoint {
                    k: [
                        (i as f64) / (grid[0] as f64),
                        (j as f64) / (grid[1] as f64),
                        (k as f64) / (grid[2] as f64),
                    ],
                    weight,
                });
            }
        }
    }
    points
}

fn normalize_wannier_nscf_occupations(
    nscf_calc: &mut QECalculation,
    config: &WannierCalculationConfig,
    notes: &mut Vec<String>,
) -> Result<(), String> {
    match nscf_calc.system.occupations {
        Occupations::FromInput => Err(
            "Wannier v1 does not support source calculations with occupations='from_input'."
                .to_string(),
        ),
        Occupations::Tetrahedra | Occupations::TetrahedraLin | Occupations::TetrahedraOpt => {
            if config.num_bands == config.num_wann {
                nscf_calc.system.occupations = Occupations::Fixed;
                nscf_calc.system.degauss = None;
                notes.push(
                    "Using fixed occupations for the isolated Wannier NSCF because QE tetrahedra occupations require an automatically generated k-point grid, while the Wannier interface needs an explicit full K_POINTS crystal list."
                        .to_string(),
                );
            } else {
                nscf_calc.system.occupations = Occupations::Smearing;
                nscf_calc.system.smearing = SmearingType::Gaussian;
                if nscf_calc
                    .system
                    .degauss
                    .map(|value| value <= 0.0)
                    .unwrap_or(true)
                {
                    nscf_calc.system.degauss = Some(0.02);
                }
                notes.push(
                    "Converted tetrahedra occupations to gaussian smearing for the Wannier NSCF because QE tetrahedra occupations require an automatically generated k-point grid, while the Wannier interface needs an explicit full K_POINTS crystal list."
                        .to_string(),
                );
            }
            Ok(())
        }
        Occupations::Smearing => {
            if config.num_bands == config.num_wann {
                nscf_calc.system.occupations = Occupations::Fixed;
                nscf_calc.system.degauss = None;
                notes.push(
                    "Using fixed occupations for the isolated Wannier NSCF because num_bands = num_wann, avoiding unnecessary Fermi-level bracketing on an isolated manifold."
                        .to_string(),
                );
            } else if nscf_calc
                .system
                .degauss
                .map(|value| value <= 0.0)
                .unwrap_or(true)
            {
                nscf_calc.system.degauss = Some(0.02);
                notes.push(
                    "Applied default degauss = 0.02 Ry for the Wannier NSCF smearing run."
                        .to_string(),
                );
            }
            Ok(())
        }
        Occupations::Fixed => {
            nscf_calc.system.degauss = None;
            Ok(())
        }
    }
}

pub fn prepare_wannier_nscf_calculation(
    config: &WannierCalculationConfig,
) -> Result<(QECalculation, Vec<String>), String> {
    validate_wannier_config(config)?;

    let mut nscf_calc = config.base_calculation.clone();
    let mut notes = Vec::new();

    nscf_calc.calculation = CalculationType::Nscf;
    nscf_calc.system.nbnd = Some(config.num_bands);
    nscf_calc.system.nosym = true;
    nscf_calc.system.noinv = true;
    nscf_calc.kpoints = KPoints::Crystal {
        points: generate_uniform_mp_kpoints(config.k_grid),
    };
    nscf_calc.startingpot = Some(StartingPotential::File);
    nscf_calc.diago_full_acc = true;
    nscf_calc.verbosity = Some("high".to_string());

    normalize_wannier_nscf_occupations(&mut nscf_calc, config, &mut notes)?;

    Ok((nscf_calc, notes))
}

fn render_bool(value: bool) -> &'static str {
    if value {
        "true"
    } else {
        "false"
    }
}

fn format_projection(spec: &WannierProjectionSpec) -> Result<String, String> {
    let orbital = spec.orbital.trim().to_ascii_lowercase();
    if orbital_template_count(&orbital).is_none() {
        return Err(format!(
            "Unsupported Wannier projection orbital template: {}",
            spec.orbital
        ));
    }

    match spec.target_type {
        WannierProjectionTargetType::Element => {
            let symbol = spec
                .symbol
                .as_deref()
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .ok_or_else(|| {
                    "Element-targeted Wannier projections require a symbol.".to_string()
                })?;
            Ok(format!("{}:{}", symbol, orbital))
        }
        WannierProjectionTargetType::Site => {
            let position = spec.fractional_position.ok_or_else(|| {
                "Site-targeted Wannier projections require fractional_position.".to_string()
            })?;
            Ok(format!(
                "f={:.10},{:.10},{:.10}:{}",
                position[0], position[1], position[2], orbital
            ))
        }
    }
}

fn format_exclude_bands(exclude_bands: &[u32]) -> Option<String> {
    if exclude_bands.is_empty() {
        return None;
    }
    let mut values: Vec<u32> = exclude_bands
        .iter()
        .copied()
        .filter(|value| *value > 0)
        .collect();
    values.sort_unstable();
    values.dedup();
    if values.is_empty() {
        return None;
    }

    let mut chunks: Vec<String> = Vec::new();
    let mut start = values[0];
    let mut prev = values[0];
    for value in values.into_iter().skip(1) {
        if value == prev + 1 {
            prev = value;
            continue;
        }
        if start == prev {
            chunks.push(start.to_string());
        } else {
            chunks.push(format!("{}-{}", start, prev));
        }
        start = value;
        prev = value;
    }
    if start == prev {
        chunks.push(start.to_string());
    } else {
        chunks.push(format!("{}-{}", start, prev));
    }
    Some(chunks.join(", "))
}

pub fn generate_wannier90_win(
    config: &WannierCalculationConfig,
    kpoints: &[KPoint],
) -> Result<String, String> {
    validate_wannier_config(config)?;
    let system = &config.base_calculation.system;
    let cell = convert_cell_matrix_to_angstrom(system)?;

    let mut output = String::new();
    writeln!(output, "num_wann = {}", config.num_wann).unwrap();
    writeln!(output, "num_bands = {}", config.num_bands).unwrap();
    writeln!(
        output,
        "mp_grid = {} {} {}",
        config.k_grid[0], config.k_grid[1], config.k_grid[2]
    )
    .unwrap();
    writeln!(
        output,
        "guiding_centres = {}",
        render_bool(config.guiding_centres)
    )
    .unwrap();
    writeln!(
        output,
        "use_ws_distance = {}",
        render_bool(config.use_ws_distance)
    )
    .unwrap();
    writeln!(output, "write_hr = {}", render_bool(config.write_hr)).unwrap();
    writeln!(output, "write_xyz = {}", render_bool(config.write_xyz)).unwrap();
    writeln!(output, "bands_plot = {}", render_bool(config.bands_plot)).unwrap();
    writeln!(output, "bands_plot_format = gnuplot").unwrap();
    writeln!(
        output,
        "bands_num_points = {}",
        config.band_path.bands_num_points
    )
    .unwrap();
    writeln!(output, "conv_window = {}", config.conv_window).unwrap();
    writeln!(output, "conv_tol = {:.1e}", config.conv_tol).unwrap();
    writeln!(output, "num_iter = {}", config.num_iter).unwrap();

    if config.num_bands > config.num_wann {
        if let Some(disentanglement) = config.disentanglement.as_ref() {
            if let Some(exclude_bands) = format_exclude_bands(&disentanglement.exclude_bands) {
                writeln!(output, "exclude_bands = {}", exclude_bands).unwrap();
            }
            if let Some(value) = disentanglement.dis_win_min {
                writeln!(output, "dis_win_min = {}", value).unwrap();
            }
            if let Some(value) = disentanglement.dis_win_max {
                writeln!(output, "dis_win_max = {}", value).unwrap();
            }
            if let Some(value) = disentanglement.dis_froz_min {
                writeln!(output, "dis_froz_min = {}", value).unwrap();
            }
            if let Some(value) = disentanglement.dis_froz_max {
                writeln!(output, "dis_froz_max = {}", value).unwrap();
            }
        }
    }

    writeln!(output).unwrap();
    writeln!(output, "begin unit_cell_cart").unwrap();
    writeln!(output, "angstrom").unwrap();
    for row in &cell {
        writeln!(output, "{:16.10} {:16.10} {:16.10}", row[0], row[1], row[2]).unwrap();
    }
    writeln!(output, "end unit_cell_cart").unwrap();
    writeln!(output).unwrap();

    writeln!(output, "begin atoms_frac").unwrap();
    for atom in &system.atoms {
        writeln!(
            output,
            "{} {:16.10} {:16.10} {:16.10}",
            atom.symbol, atom.position[0], atom.position[1], atom.position[2]
        )
        .unwrap();
    }
    writeln!(output, "end atoms_frac").unwrap();
    writeln!(output).unwrap();

    writeln!(output, "begin projections").unwrap();
    for projection in &config.projections {
        writeln!(output, "{}", format_projection(projection)?).unwrap();
    }
    writeln!(output, "end projections").unwrap();
    writeln!(output).unwrap();

    writeln!(output, "begin kpoints").unwrap();
    for point in kpoints {
        writeln!(
            output,
            "{:16.10} {:16.10} {:16.10}",
            point.k[0], point.k[1], point.k[2]
        )
        .unwrap();
    }
    writeln!(output, "end kpoints").unwrap();
    writeln!(output).unwrap();

    writeln!(output, "begin kpoint_path").unwrap();
    for window in config.band_path.k_path.windows(2) {
        let from = &window[0];
        let to = &window[1];
        writeln!(
            output,
            "{} {:16.10} {:16.10} {:16.10} {} {:16.10} {:16.10} {:16.10}",
            from.label,
            from.coords[0],
            from.coords[1],
            from.coords[2],
            to.label,
            to.coords[0],
            to.coords[1],
            to.coords[2]
        )
        .unwrap();
    }
    writeln!(output, "end kpoint_path").unwrap();

    Ok(output)
}

pub fn generate_pw2wannier90_input(
    config: &WannierCalculationConfig,
    pw2_config: &Pw2Wannier90Config,
) -> String {
    let mut output = String::new();
    writeln!(output, "&inputpp").unwrap();
    writeln!(output, "  prefix = '{}',", config.base_calculation.prefix).unwrap();
    writeln!(output, "  outdir = '{}',", config.base_calculation.outdir).unwrap();
    writeln!(output, "  seedname = '{}',", config.seedname).unwrap();
    writeln!(
        output,
        "  write_amn = .{}.,",
        if pw2_config.write_amn {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(
        output,
        "  write_mmn = .{}.,",
        if pw2_config.write_mmn {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(
        output,
        "  write_spn = .{}.,",
        if pw2_config.write_spn {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(
        output,
        "  write_unk = .{}.,",
        if pw2_config.write_unk {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(
        output,
        "  write_dmn = .{}.,",
        if pw2_config.write_dmn {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(
        output,
        "  spin_component = '{}',",
        pw2_config.spin_component
    )
    .unwrap();
    writeln!(
        output,
        "  scdm_proj = .{}.,",
        if pw2_config.scdm_proj {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(
        output,
        "  atom_proj = .{}.,",
        if pw2_config.atom_proj {
            "true"
        } else {
            "false"
        }
    )
    .unwrap();
    writeln!(output, "/").unwrap();
    output
}

pub fn parse_wannier_wout(
    content: &str,
) -> Result<(Vec<WannierSpread>, WannierConvergenceData), String> {
    let spread_re = Regex::new(
        r"WF centre and spread\s+(\d+)\s+\(\s*([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*\)\s*([-\d.Ee+]+)",
    )
    .map_err(|e| format!("Failed to compile Wannier spread regex: {}", e))?;
    let omega_re = Regex::new(r"Omega\s+(I|D|OD|Total)\s*=\s*([-\d.Ee+]+)")
        .map_err(|e| format!("Failed to compile Wannier omega regex: {}", e))?;
    let iter_re = Regex::new(r"^\s*(\d+)\s+[-\d.Ee+]+\s+[-\d.Ee+]+\s+[-\d.Ee+]+")
        .map_err(|e| format!("Failed to compile Wannier iteration regex: {}", e))?;

    let mut spreads = Vec::new();
    let mut convergence = WannierConvergenceData::default();
    for line in content.lines() {
        if let Some(caps) = spread_re.captures(line) {
            let index = caps
                .get(1)
                .and_then(|m| m.as_str().parse::<u32>().ok())
                .unwrap_or(0);
            let centre = [
                caps.get(2)
                    .and_then(|m| m.as_str().parse::<f64>().ok())
                    .unwrap_or(0.0),
                caps.get(3)
                    .and_then(|m| m.as_str().parse::<f64>().ok())
                    .unwrap_or(0.0),
                caps.get(4)
                    .and_then(|m| m.as_str().parse::<f64>().ok())
                    .unwrap_or(0.0),
            ];
            let spread = caps
                .get(5)
                .and_then(|m| m.as_str().parse::<f64>().ok())
                .unwrap_or(0.0);
            spreads.push(WannierSpread {
                index,
                centre,
                spread,
            });
            continue;
        }

        if let Some(caps) = omega_re.captures(line) {
            let label = caps.get(1).map(|m| m.as_str()).unwrap_or("");
            let value = caps.get(2).and_then(|m| m.as_str().parse::<f64>().ok());
            match label {
                "I" => convergence.omega_i = value,
                "D" => convergence.omega_d = value,
                "OD" => convergence.omega_od = value,
                "Total" => convergence.omega_total = value,
                _ => {}
            }
            continue;
        }

        if line.contains("<-- CONV") {
            if let Some(caps) = iter_re.captures(line) {
                convergence.iterations = caps.get(1).and_then(|m| m.as_str().parse::<u32>().ok());
            }
            convergence.converged = true;
        }
    }

    if convergence.omega_total.is_none() && !spreads.is_empty() {
        convergence.omega_total = Some(spreads.iter().map(|entry| entry.spread).sum());
    }

    if spreads.is_empty() {
        return Err("No Wannier centre/spread entries found in .wout.".to_string());
    }

    Ok((spreads, convergence))
}

pub fn parse_wannier_centres_xyz(content: &str) -> Result<Vec<WannierCentreRecord>, String> {
    let mut records = Vec::new();
    for (index, line) in content.lines().enumerate() {
        if index < 2 {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let x = parts[1]
            .parse::<f64>()
            .map_err(|_| "Failed to parse centre x coordinate".to_string())?;
        let y = parts[2]
            .parse::<f64>()
            .map_err(|_| "Failed to parse centre y coordinate".to_string())?;
        let z = parts[3]
            .parse::<f64>()
            .map_err(|_| "Failed to parse centre z coordinate".to_string())?;
        records.push(WannierCentreRecord {
            label: parts[0].to_string(),
            position: [x, y, z],
        });
    }
    Ok(records)
}

pub fn parse_wannier_band_data(
    content: &str,
    fermi_energy: f64,
    k_path: &[KPathPoint],
) -> Result<BandData, String> {
    match parse_bands_gnu(content, fermi_energy) {
        Ok(mut data) => {
            add_symmetry_markers(&mut data, k_path);
            Ok(data)
        }
        Err(_) => {
            let rows: Vec<Vec<f64>> = content
                .lines()
                .filter_map(|line| {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        return None;
                    }
                    let parsed: Option<Vec<f64>> = trimmed
                        .split_whitespace()
                        .map(|token| token.parse::<f64>().ok())
                        .collect();
                    parsed
                })
                .collect();
            if rows.is_empty() {
                return Err("Empty Wannier band file.".to_string());
            }
            let width = rows[0].len();
            if width < 2 || rows.iter().any(|row| row.len() != width) {
                return Err("Unrecognized Wannier band file format.".to_string());
            }

            let k_points: Vec<f64> = rows.iter().map(|row| row[0]).collect();
            let n_kpoints = k_points.len();
            let n_bands = width - 1;
            let mut energies = vec![vec![0.0; n_kpoints]; n_bands];
            for (k_index, row) in rows.iter().enumerate() {
                for band_index in 0..n_bands {
                    energies[band_index][k_index] = row[band_index + 1];
                }
            }

            let mut e_min = f64::INFINITY;
            let mut e_max = f64::NEG_INFINITY;
            for band in &energies {
                for energy in band {
                    if *energy < e_min {
                        e_min = *energy;
                    }
                    if *energy > e_max {
                        e_max = *energy;
                    }
                }
            }

            let mut data = BandData {
                k_points,
                energies,
                fermi_energy,
                high_symmetry_points: Vec::new(),
                n_bands,
                n_kpoints,
                band_gap: None,
                energy_range: [e_min, e_max],
                projections: None,
            };
            add_symmetry_markers(&mut data, k_path);
            Ok(data)
        }
    }
}

pub fn collect_wannier_artifacts(work_path: &Path, seedname: &str) -> Vec<WannierArtifact> {
    let file_names = [
        format!("{}.win", seedname),
        format!("{}.nnkp", seedname),
        format!("{}.amn", seedname),
        format!("{}.mmn", seedname),
        format!("{}.eig", seedname),
        format!("{}_wsvec.dat", seedname),
        format!("{}.wout", seedname),
        format!("{}.chk", seedname),
        format!("{}_hr.dat", seedname),
        format!("{}_centres.xyz", seedname),
        format!("{}_band.dat", seedname),
        format!("{}_band.kpt", seedname),
        "nscf.in".to_string(),
        "nscf.out".to_string(),
        "pw2wan.in".to_string(),
        "pw2wan.out".to_string(),
        "wannier90_pre.out".to_string(),
        "wannier90.out".to_string(),
        "run.sbatch".to_string(),
        "slurm.out".to_string(),
        "slurm.err".to_string(),
    ];
    let mut artifacts = Vec::new();
    for file_name in file_names {
        let path = work_path.join(&file_name);
        if let Ok(metadata) = std::fs::metadata(&path) {
            if metadata.is_file() {
                artifacts.push(WannierArtifact {
                    file_name,
                    size_bytes: metadata.len(),
                });
            }
        }
    }
    artifacts
}

pub fn read_wannier_result(
    work_path: &Path,
    config: &WannierCalculationConfig,
    fermi_energy: f64,
) -> Result<WannierResult, String> {
    let wout_path = work_path.join(format!("{}.wout", config.seedname));
    let band_path = work_path.join(format!("{}_band.dat", config.seedname));
    if !wout_path.exists() {
        return Err(format!(
            "Wannier output file not found: {}",
            wout_path.display()
        ));
    }
    if !band_path.exists() {
        return Err(format!(
            "Interpolated Wannier band file not found: {}",
            band_path.display()
        ));
    }

    let wout_content = std::fs::read_to_string(&wout_path)
        .map_err(|e| format!("Failed to read {}: {}", wout_path.display(), e))?;
    let band_content = std::fs::read_to_string(&band_path)
        .map_err(|e| format!("Failed to read {}: {}", band_path.display(), e))?;

    let (spreads, convergence) = parse_wannier_wout(&wout_content)?;
    let band_data = parse_wannier_band_data(&band_content, fermi_energy, &config.band_path.k_path)?;
    let centres_path = work_path.join(format!("{}_centres.xyz", config.seedname));
    let centres = if centres_path.exists() {
        std::fs::read_to_string(&centres_path)
            .ok()
            .and_then(|content| parse_wannier_centres_xyz(&content).ok())
            .unwrap_or_default()
    } else {
        spreads
            .iter()
            .map(|entry| WannierCentreRecord {
                label: format!("WF{}", entry.index),
                position: entry.centre,
            })
            .collect()
    };

    Ok(WannierResult {
        seedname: config.seedname.clone(),
        num_wann: config.num_wann,
        num_bands: config.num_bands,
        k_grid: config.k_grid,
        band_data,
        total_spread: convergence.omega_total,
        spreads,
        centres,
        convergence,
        artifact_manifest: collect_wannier_artifacts(work_path, &config.seedname),
    })
}

fn parse_wannier_hr(content: &str) -> Result<(u32, Vec<WannierHrTerm>), String> {
    let mut lines = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty());
    let _header = lines
        .next()
        .ok_or_else(|| "Wannier hr.dat file is empty.".to_string())?;
    let num_wann = lines
        .next()
        .ok_or_else(|| "Wannier hr.dat missing num_wann line.".to_string())
        .and_then(|line| {
            line.parse::<u32>()
                .map_err(|_| format!("Invalid num_wann entry in hr.dat: {}", line))
        })?;
    let nrpts = lines
        .next()
        .ok_or_else(|| "Wannier hr.dat missing nrpts line.".to_string())
        .and_then(|line| {
            line.parse::<usize>()
                .map_err(|_| format!("Invalid nrpts entry in hr.dat: {}", line))
        })?;

    let mut degeneracies = Vec::with_capacity(nrpts);
    while degeneracies.len() < nrpts {
        let line = lines
            .next()
            .ok_or_else(|| "Unexpected end of hr.dat while reading degeneracies.".to_string())?;
        for token in line.split_whitespace() {
            degeneracies.push(
                token
                    .parse::<u32>()
                    .map_err(|_| format!("Invalid hr.dat degeneracy entry '{}'.", token))?,
            );
            if degeneracies.len() == nrpts {
                break;
            }
        }
    }

    let expected_terms = nrpts
        .checked_mul(num_wann as usize)
        .and_then(|value| value.checked_mul(num_wann as usize))
        .ok_or_else(|| "Wannier hr.dat term count overflow.".to_string())?;

    let mut terms = Vec::with_capacity(expected_terms);
    for term_index in 0..expected_terms {
        let line = lines.next().ok_or_else(|| {
            "Unexpected end of hr.dat while reading Hamiltonian terms.".to_string()
        })?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 7 {
            return Err(format!("Malformed hr.dat term line: '{}'.", line));
        }
        let r = [
            parse_i32_token(parts[0], "R1")?,
            parse_i32_token(parts[1], "R2")?,
            parse_i32_token(parts[2], "R3")?,
        ];
        let m = parse_usize_token(parts[3], "m")?
            .checked_sub(1)
            .ok_or_else(|| "Wannier hr.dat orbital indices must start at 1.".to_string())?;
        let n = parse_usize_token(parts[4], "n")?
            .checked_sub(1)
            .ok_or_else(|| "Wannier hr.dat orbital indices must start at 1.".to_string())?;
        let real = parse_f64_token(parts[5], "real hopping")?;
        let imag = parse_f64_token(parts[6], "imag hopping")?;
        let degeneracy = degeneracies[term_index / ((num_wann as usize) * (num_wann as usize))];
        if degeneracy == 0 {
            return Err("Encountered zero degeneracy in hr.dat.".to_string());
        }
        terms.push(WannierHrTerm {
            r,
            m,
            n,
            degeneracy,
            real,
            imag,
        });
    }

    Ok((num_wann, terms))
}

fn parse_wannier_wsvec(content: &str) -> Result<HashMap<String, Vec<[i32; 3]>>, String> {
    let mut entries = HashMap::new();
    let mut lines = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#') && !line.starts_with('!'));

    while let Some(header_line) = lines.next() {
        let header_parts: Vec<&str> = header_line.split_whitespace().collect();
        if header_parts.len() < 5 {
            return Err(format!("Malformed wsvec header line: '{}'.", header_line));
        }
        let r = [
            parse_i32_token(header_parts[0], "wsvec R1")?,
            parse_i32_token(header_parts[1], "wsvec R2")?,
            parse_i32_token(header_parts[2], "wsvec R3")?,
        ];
        let m = parse_usize_token(header_parts[3], "wsvec m")?
            .checked_sub(1)
            .ok_or_else(|| "wsvec orbital indices must start at 1.".to_string())?;
        let n = parse_usize_token(header_parts[4], "wsvec n")?
            .checked_sub(1)
            .ok_or_else(|| "wsvec orbital indices must start at 1.".to_string())?;
        let count_line = lines.next().ok_or_else(|| {
            "Unexpected end of wsvec.dat while reading translation count.".to_string()
        })?;
        let count = count_line
            .parse::<usize>()
            .map_err(|_| format!("Invalid wsvec multiplicity '{}'.", count_line))?;
        let mut translations = Vec::with_capacity(count);
        for _ in 0..count {
            let vector_line = lines.next().ok_or_else(|| {
                "Unexpected end of wsvec.dat while reading translation vectors.".to_string()
            })?;
            let parts: Vec<&str> = vector_line.split_whitespace().collect();
            if parts.len() < 3 {
                return Err(format!(
                    "Malformed wsvec translation line: '{}'.",
                    vector_line
                ));
            }
            translations.push([
                parse_i32_token(parts[0], "wsvec T1")?,
                parse_i32_token(parts[1], "wsvec T2")?,
                parse_i32_token(parts[2], "wsvec T3")?,
            ]);
        }
        entries.insert(wsvec_key(r, m, n), translations);
    }

    Ok(entries)
}

pub fn parse_wannier_hamiltonian(
    seedname: &str,
    win_content: &str,
    hr_content: &str,
    wsvec_content: Option<&str>,
) -> Result<WannierHamiltonian, String> {
    let lattice_vectors_angstrom = parse_win_unit_cell_angstrom(win_content)?;
    let use_ws_distance = parse_win_bool_assignment(win_content, "use_ws_distance")
        .unwrap_or_else(|| wsvec_content.is_some());
    let (num_wann, hr_terms) = parse_wannier_hr(hr_content)?;
    let ws_translations = if use_ws_distance {
        let content = wsvec_content.ok_or_else(|| {
            "Wannier export requires seedname_wsvec.dat because use_ws_distance = true.".to_string()
        })?;
        parse_wannier_wsvec(content)?
    } else {
        HashMap::new()
    };

    Ok(WannierHamiltonian {
        seedname: seedname.to_string(),
        num_wann,
        lattice_vectors_angstrom,
        hr_terms,
        use_ws_distance,
        ws_translations,
    })
}

impl WannierHamiltonian {
    fn phase_factor(&self, term: &WannierHrTerm, k_fractional: [f64; 3]) -> Complex64 {
        let translations = if self.use_ws_distance {
            self.ws_translations
                .get(&wsvec_key(term.r, term.m, term.n))
                .map(|values| values.as_slice())
        } else {
            None
        };
        let translations = translations.unwrap_or(&[[0, 0, 0]]);
        let mut phase_sum = Complex64::new(0.0, 0.0);
        for translation in translations {
            let r_total = [
                term.r[0] + translation[0],
                term.r[1] + translation[1],
                term.r[2] + translation[2],
            ];
            let angle = TWO_PI
                * (k_fractional[0] * f64::from(r_total[0])
                    + k_fractional[1] * f64::from(r_total[1])
                    + k_fractional[2] * f64::from(r_total[2]));
            phase_sum += Complex64::from_polar(1.0, angle);
        }
        phase_sum / (translations.len() as f64)
    }

    fn hamiltonian_matrix(&self, k_fractional: [f64; 3]) -> DMatrix<Complex64> {
        let size = self.num_wann as usize;
        let mut matrix = DMatrix::<Complex64>::zeros(size, size);
        for term in &self.hr_terms {
            let hopping = Complex64::new(term.real, term.imag);
            let phase = self.phase_factor(term, k_fractional);
            matrix[(term.m, term.n)] += hopping * phase / f64::from(term.degeneracy);
        }
        for row in 0..size {
            for col in 0..row {
                let averaged = (matrix[(row, col)] + matrix[(col, row)].conj()) * 0.5;
                matrix[(row, col)] = averaged;
                matrix[(col, row)] = averaged.conj();
            }
        }
        matrix
    }

    pub fn eigenvalues_at(&self, k_fractional: [f64; 3]) -> Result<Vec<f64>, String> {
        let matrix = self.hamiltonian_matrix(k_fractional);
        let eigenvalues = matrix
            .eigenvalues()
            .ok_or_else(|| "Failed to diagonalize Wannier Hamiltonian.".to_string())?;
        let mut values: Vec<f64> = eigenvalues.iter().map(|value| value.re).collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        Ok(values)
    }
}

fn sanitize_fragment(raw: &str) -> String {
    let sanitized: String = raw
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                ch
            } else {
                '_'
            }
        })
        .collect();
    let trimmed = sanitized.trim_matches('_');
    if trimmed.is_empty() {
        "export".to_string()
    } else {
        trimmed.to_string()
    }
}

fn unique_export_dir(root: &Path, base_name: &str) -> Result<PathBuf, String> {
    std::fs::create_dir_all(root)
        .map_err(|e| format!("Failed to create export root {}: {}", root.display(), e))?;
    let sanitized = sanitize_fragment(base_name);
    for suffix in 0..1000 {
        let candidate = if suffix == 0 {
            root.join(&sanitized)
        } else {
            root.join(format!("{}_{}", sanitized, suffix))
        };
        if !candidate.exists() {
            std::fs::create_dir_all(&candidate).map_err(|e| {
                format!(
                    "Failed to create Ludwig export bundle directory {}: {}",
                    candidate.display(),
                    e
                )
            })?;
            return Ok(candidate);
        }
    }
    Err("Could not allocate a unique Ludwig export directory.".to_string())
}

pub fn export_ludwig_bundle(
    hamiltonian: &WannierHamiltonian,
    config: &LudwigExportConfig,
    chemical_potential_ev: f64,
    provenance_files: &[PathBuf],
) -> Result<LudwigExportResult, String> {
    if config.nkx < 2 || config.nky < 2 {
        return Err("Ludwig export requires nkx and nky to both be at least 2.".to_string());
    }
    let axis_a = usize::from(config.in_plane_axes[0]);
    let axis_b = usize::from(config.in_plane_axes[1]);
    let slice_axis = usize::from(config.slice_axis);
    if axis_a >= 3 || axis_b >= 3 || slice_axis >= 3 {
        return Err("Ludwig export axes must each be one of 0, 1, 2.".to_string());
    }
    if axis_a == axis_b || axis_a == slice_axis || axis_b == slice_axis {
        return Err(
            "Ludwig export requires a permutation of three distinct reciprocal axes.".to_string(),
        );
    }

    let slice_coordinate = config.slice_coordinate.rem_euclid(1.0);
    let in_plane_lattice_angstrom =
        embed_in_plane_lattice(&hamiltonian.lattice_vectors_angstrom, config.in_plane_axes)?;
    let destination_root = PathBuf::from(config.destination_root.trim());
    if destination_root.as_os_str().is_empty() {
        return Err("Ludwig export destination root is required.".to_string());
    }
    let bundle_dir = unique_export_dir(
        &destination_root,
        &format!("{}_ludwig_{}", hamiltonian.seedname, config.calculation_id),
    )?;

    let mut csv = String::new();
    csv.push_str("ix,iy,kx_frac,ky_frac");
    for band_index in 0..hamiltonian.num_wann {
        let _ = write!(csv, ",band_{}", band_index + 1);
    }
    csv.push('\n');

    for iy in 0..config.nky {
        let ky_frac = f64::from(iy) / f64::from(config.nky);
        for ix in 0..config.nkx {
            let kx_frac = f64::from(ix) / f64::from(config.nkx);
            let mut k_fractional = [0.0; 3];
            k_fractional[axis_a] = kx_frac;
            k_fractional[axis_b] = ky_frac;
            k_fractional[slice_axis] = slice_coordinate;
            let eigenvalues = hamiltonian.eigenvalues_at(k_fractional)?;
            let _ = write!(csv, "{},{},{:.10},{:.10}", ix, iy, kx_frac, ky_frac);
            for energy in eigenvalues {
                let _ = write!(csv, ",{:.12}", energy - chemical_potential_ev);
            }
            csv.push('\n');
        }
    }

    let mut copied_provenance = Vec::new();
    for source in provenance_files {
        if !source.exists() {
            continue;
        }
        let Some(file_name) = source.file_name().and_then(|value| value.to_str()) else {
            continue;
        };
        let destination = bundle_dir.join(file_name);
        std::fs::copy(source, &destination).map_err(|e| {
            format!(
                "Failed to copy Ludwig provenance file {} to {}: {}",
                source.display(),
                destination.display(),
                e
            )
        })?;
        copied_provenance.push(file_name.to_string());
    }

    let metadata = LudwigExportMetadata {
        project_id: config.project_id.clone(),
        calculation_id: config.calculation_id.clone(),
        seedname: hamiltonian.seedname.clone(),
        mode: config.mode.clone(),
        in_plane_axes: config.in_plane_axes,
        slice_axis: config.slice_axis,
        slice_coordinate,
        chemical_potential_ev,
        grid_shape: [config.nkx, config.nky],
        band_count: hamiltonian.num_wann,
        lattice_vectors_angstrom: hamiltonian.lattice_vectors_angstrom,
        in_plane_lattice_angstrom,
        fractional_domain: "0 <= kx_frac, ky_frac < 1 in the selected reciprocal-basis directions"
            .to_string(),
        energy_reference:
            "Band energies are shifted by the exported chemical potential so Ludwig sees EF = 0."
                .to_string(),
        provenance_files: copied_provenance,
    };

    let metadata_path = bundle_dir.join("metadata.json");
    let metadata_json = serde_json::to_string_pretty(&metadata)
        .map_err(|e| format!("Failed to serialize Ludwig export metadata: {}", e))?;
    std::fs::write(&metadata_path, metadata_json).map_err(|e| {
        format!(
            "Failed to write Ludwig metadata file {}: {}",
            metadata_path.display(),
            e
        )
    })?;

    let bands_path = bundle_dir.join("bands.csv");
    std::fs::write(&bands_path, csv).map_err(|e| {
        format!(
            "Failed to write Ludwig band grid {}: {}",
            bands_path.display(),
            e
        )
    })?;

    Ok(LudwigExportResult {
        bundle_path: bundle_dir.display().to_string(),
        band_count: hamiltonian.num_wann,
        chemical_potential_ev,
        grid_shape: [config.nkx, config.nky],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qe::types::{
        Atom, AtomicSpecies, BravaisLattice, CalculationType, KPoints, Occupations, PositionUnits,
        QECalculation, QESystem, SmearingType,
    };
    use std::time::{SystemTime, UNIX_EPOCH};

    fn sample_base_calculation() -> QECalculation {
        QECalculation {
            calculation: CalculationType::Scf,
            prefix: "qcortado_scf".to_string(),
            outdir: "./tmp".to_string(),
            pseudo_dir: "./pseudo".to_string(),
            disk_io: None,
            system: QESystem {
                ibrav: BravaisLattice::Free,
                celldm: None,
                cell_parameters: Some([
                    [0.0, 2.715, 2.715],
                    [2.715, 0.0, 2.715],
                    [2.715, 2.715, 0.0],
                ]),
                cell_units: Some(PositionUnits::Angstrom),
                species: vec![AtomicSpecies {
                    symbol: "Si".to_string(),
                    mass: 28.0855,
                    pseudopotential: "Si.upf".to_string(),
                }],
                atoms: vec![
                    Atom {
                        symbol: "Si".to_string(),
                        position: [0.0, 0.0, 0.0],
                        if_pos: [true, true, true],
                    },
                    Atom {
                        symbol: "Si".to_string(),
                        position: [0.25, 0.25, 0.25],
                        if_pos: [true, true, true],
                    },
                ],
                position_units: PositionUnits::Crystal,
                ecutwfc: 40.0,
                ecutrho: Some(320.0),
                nbnd: Some(8),
                tot_charge: None,
                input_dft: None,
                nspin: 1,
                noncolin: false,
                lspinorb: false,
                occupations: Occupations::Fixed,
                smearing: SmearingType::Gaussian,
                degauss: None,
                nosym: false,
                noinv: false,
            },
            kpoints: KPoints::Gamma,
            conv_thr: 1.0e-8,
            electron_maxstep: None,
            mixing_mode: None,
            mixing_beta: 0.7,
            mixing_ndim: None,
            diagonalization: None,
            startingpot: None,
            startingwfc: None,
            diago_full_acc: false,
            tprnfor: false,
            tstress: false,
            forc_conv_thr: None,
            etot_conv_thr: None,
            press: None,
            verbosity: Some("high".to_string()),
        }
    }

    fn sample_site_projection() -> WannierProjectionSpec {
        WannierProjectionSpec {
            target_type: WannierProjectionTargetType::Site,
            symbol: Some("Si".to_string()),
            orbital: "sp3".to_string(),
            site_index: Some(1),
            fractional_position: Some([0.0, 0.0, 0.0]),
        }
    }

    #[test]
    fn uniform_grid_is_zero_shifted() {
        let grid = generate_uniform_mp_kpoints([2, 2, 1]);
        assert_eq!(grid.len(), 4);
        assert_eq!(grid[0].k, [0.0, 0.0, 0.0]);
        assert_eq!(grid[1].k, [0.0, 0.5, 0.0]);
        assert_eq!(grid[2].k, [0.5, 0.0, 0.0]);
        assert!((grid[0].weight - 0.25).abs() < 1.0e-12);
    }

    #[test]
    fn win_writer_renders_required_blocks() {
        let config = WannierCalculationConfig {
            base_calculation: sample_base_calculation(),
            k_grid: [2, 2, 2],
            num_wann: 4,
            num_bands: 8,
            seedname: "si_mlwf".to_string(),
            projections: vec![sample_site_projection()],
            band_path: WannierBandPathConfig {
                k_path: vec![
                    KPathPoint {
                        label: "G".to_string(),
                        coords: [0.0, 0.0, 0.0],
                        npoints: 20,
                    },
                    KPathPoint {
                        label: "X".to_string(),
                        coords: [0.5, 0.0, 0.5],
                        npoints: 0,
                    },
                ],
                bands_num_points: 100,
            },
            disentanglement: None,
            pw2wannier90: None,
            project_id: None,
            scf_calc_id: None,
            source_metadata: Some(WannierSourceMetadata {
                cell_representation: Some("primitive_spglib".to_string()),
                electron_count: Some(8.0),
                nspin: Some(1),
                noncolin: Some(false),
                lspinorb: Some(false),
                lda_plus_u: Some(false),
                vdw_corr: Some("none".to_string()),
            }),
            guiding_centres: true,
            use_ws_distance: true,
            write_hr: true,
            write_xyz: true,
            bands_plot: true,
            conv_window: 5,
            conv_tol: 1.0e-10,
            num_iter: 100,
        };

        let win = generate_wannier90_win(&config, &generate_uniform_mp_kpoints([2, 2, 2])).unwrap();
        assert!(win.contains("num_wann = 4"));
        assert!(win.contains("mp_grid = 2 2 2"));
        assert!(win.contains("begin unit_cell_cart\nangstrom"));
        assert!(win.contains("begin kpoints"));
        assert!(win.contains("begin kpoint_path"));
        assert!(win.contains("f=0.0000000000,0.0000000000,0.0000000000:sp3"));
    }

    #[test]
    fn wout_parser_extracts_spreads() {
        let content = r#"
 Final State
  WF centre and spread    1  (   0.000000,   0.000000,   0.000000 )     1.23456789
  WF centre and spread    2  (   1.000000,   1.000000,   1.000000 )     2.00000000
  Omega I     =      1.0
  Omega D     =      0.5
  Omega OD    =      1.73456789
  Omega Total =      3.23456789
        12   0.0000  0.0000  0.0000   <-- CONV
"#;
        let (spreads, convergence) = parse_wannier_wout(content).unwrap();
        assert_eq!(spreads.len(), 2);
        assert_eq!(spreads[0].index, 1);
        assert!(convergence.converged);
        assert_eq!(convergence.iterations, Some(12));
        assert_eq!(convergence.omega_total, Some(3.23456789));
    }

    #[test]
    fn prepare_wannier_nscf_uses_fixed_for_isolated_smearing_source() {
        let mut base = sample_base_calculation();
        base.system.occupations = Occupations::Smearing;
        base.system.smearing = SmearingType::MarzariVanderbilt;
        base.system.degauss = Some(0.01);

        let config = WannierCalculationConfig {
            base_calculation: base,
            k_grid: [4, 4, 4],
            num_wann: 4,
            num_bands: 4,
            seedname: "silicon".to_string(),
            projections: vec![sample_site_projection()],
            band_path: WannierBandPathConfig {
                k_path: vec![
                    KPathPoint {
                        label: "G".to_string(),
                        coords: [0.0, 0.0, 0.0],
                        npoints: 20,
                    },
                    KPathPoint {
                        label: "X".to_string(),
                        coords: [0.5, 0.0, 0.5],
                        npoints: 0,
                    },
                ],
                bands_num_points: 100,
            },
            disentanglement: None,
            pw2wannier90: None,
            project_id: None,
            scf_calc_id: None,
            source_metadata: Some(WannierSourceMetadata {
                cell_representation: Some("primitive_spglib".to_string()),
                electron_count: Some(8.0),
                nspin: Some(1),
                noncolin: Some(false),
                lspinorb: Some(false),
                lda_plus_u: Some(false),
                vdw_corr: Some("none".to_string()),
            }),
            guiding_centres: true,
            use_ws_distance: true,
            write_hr: true,
            write_xyz: true,
            bands_plot: true,
            conv_window: 5,
            conv_tol: 1.0e-10,
            num_iter: 100,
        };

        let (nscf, notes) = prepare_wannier_nscf_calculation(&config).unwrap();
        assert_eq!(nscf.calculation, CalculationType::Nscf);
        assert_eq!(nscf.system.nbnd, Some(4));
        assert_eq!(nscf.system.occupations, Occupations::Fixed);
        assert_eq!(nscf.system.degauss, None);
        assert!(nscf.system.nosym);
        assert!(nscf.system.noinv);
        assert!(matches!(nscf.kpoints, KPoints::Crystal { .. }));
        assert!(notes.iter().any(|note| note.contains("fixed occupations")));
    }

    #[test]
    fn prepare_wannier_nscf_converts_tetrahedra_to_smearing_for_entangled_run() {
        let mut base = sample_base_calculation();
        base.system.occupations = Occupations::Tetrahedra;
        base.system.degauss = None;

        let config = WannierCalculationConfig {
            base_calculation: base,
            k_grid: [4, 4, 4],
            num_wann: 4,
            num_bands: 8,
            seedname: "silicon".to_string(),
            projections: vec![sample_site_projection()],
            band_path: WannierBandPathConfig {
                k_path: vec![
                    KPathPoint {
                        label: "G".to_string(),
                        coords: [0.0, 0.0, 0.0],
                        npoints: 20,
                    },
                    KPathPoint {
                        label: "X".to_string(),
                        coords: [0.5, 0.0, 0.5],
                        npoints: 0,
                    },
                ],
                bands_num_points: 100,
            },
            disentanglement: Some(WannierDisentanglementConfig::default()),
            pw2wannier90: None,
            project_id: None,
            scf_calc_id: None,
            source_metadata: Some(WannierSourceMetadata {
                cell_representation: Some("primitive_spglib".to_string()),
                electron_count: Some(8.0),
                nspin: Some(1),
                noncolin: Some(false),
                lspinorb: Some(false),
                lda_plus_u: Some(false),
                vdw_corr: Some("none".to_string()),
            }),
            guiding_centres: true,
            use_ws_distance: true,
            write_hr: true,
            write_xyz: true,
            bands_plot: true,
            conv_window: 5,
            conv_tol: 1.0e-10,
            num_iter: 100,
        };

        let (nscf, notes) = prepare_wannier_nscf_calculation(&config).unwrap();
        assert_eq!(nscf.system.occupations, Occupations::Smearing);
        assert_eq!(nscf.system.smearing, SmearingType::Gaussian);
        assert_eq!(nscf.system.degauss, Some(0.02));
        assert!(notes
            .iter()
            .any(|note| note.contains("tetrahedra occupations")));
    }

    #[test]
    fn prepare_wannier_nscf_rejects_from_input_occupations() {
        let mut base = sample_base_calculation();
        base.system.occupations = Occupations::FromInput;

        let config = WannierCalculationConfig {
            base_calculation: base,
            k_grid: [4, 4, 4],
            num_wann: 4,
            num_bands: 4,
            seedname: "silicon".to_string(),
            projections: vec![sample_site_projection()],
            band_path: WannierBandPathConfig {
                k_path: vec![
                    KPathPoint {
                        label: "G".to_string(),
                        coords: [0.0, 0.0, 0.0],
                        npoints: 20,
                    },
                    KPathPoint {
                        label: "X".to_string(),
                        coords: [0.5, 0.0, 0.5],
                        npoints: 0,
                    },
                ],
                bands_num_points: 100,
            },
            disentanglement: None,
            pw2wannier90: None,
            project_id: None,
            scf_calc_id: None,
            source_metadata: Some(WannierSourceMetadata {
                cell_representation: Some("primitive_spglib".to_string()),
                electron_count: Some(8.0),
                nspin: Some(1),
                noncolin: Some(false),
                lspinorb: Some(false),
                lda_plus_u: Some(false),
                vdw_corr: Some("none".to_string()),
            }),
            guiding_centres: true,
            use_ws_distance: true,
            write_hr: true,
            write_xyz: true,
            bands_plot: true,
            conv_window: 5,
            conv_tol: 1.0e-10,
            num_iter: 100,
        };

        let error = prepare_wannier_nscf_calculation(&config).unwrap_err();
        assert!(error.contains("occupations='from_input'"));
    }

    #[test]
    fn validate_wannier_config_counts_element_projections_over_all_matching_atoms() {
        let config = WannierCalculationConfig {
            base_calculation: sample_base_calculation(),
            k_grid: [2, 2, 2],
            num_wann: 4,
            num_bands: 8,
            seedname: "si_mlwf".to_string(),
            projections: vec![WannierProjectionSpec {
                target_type: WannierProjectionTargetType::Element,
                symbol: Some("Si".to_string()),
                orbital: "sp3".to_string(),
                site_index: None,
                fractional_position: None,
            }],
            band_path: WannierBandPathConfig {
                k_path: vec![
                    KPathPoint {
                        label: "G".to_string(),
                        coords: [0.0, 0.0, 0.0],
                        npoints: 20,
                    },
                    KPathPoint {
                        label: "X".to_string(),
                        coords: [0.5, 0.0, 0.5],
                        npoints: 0,
                    },
                ],
                bands_num_points: 100,
            },
            disentanglement: None,
            pw2wannier90: None,
            project_id: None,
            scf_calc_id: None,
            source_metadata: Some(WannierSourceMetadata {
                cell_representation: Some("primitive_spglib".to_string()),
                electron_count: Some(8.0),
                nspin: Some(1),
                noncolin: Some(false),
                lspinorb: Some(false),
                lda_plus_u: Some(false),
                vdw_corr: Some("none".to_string()),
            }),
            guiding_centres: true,
            use_ws_distance: true,
            write_hr: true,
            write_xyz: true,
            bands_plot: true,
            conv_window: 5,
            conv_tol: 1.0e-10,
            num_iter: 100,
        };

        let error = validate_wannier_config(&config).unwrap_err();
        assert!(error.contains("expands to 8 Wannier functions"));
    }

    #[test]
    fn validate_wannier_config_rejects_num_bands_below_occupied_manifold() {
        let config = WannierCalculationConfig {
            base_calculation: sample_base_calculation(),
            k_grid: [2, 2, 2],
            num_wann: 4,
            num_bands: 4,
            seedname: "si_mlwf".to_string(),
            projections: vec![sample_site_projection()],
            band_path: WannierBandPathConfig {
                k_path: vec![
                    KPathPoint {
                        label: "G".to_string(),
                        coords: [0.0, 0.0, 0.0],
                        npoints: 20,
                    },
                    KPathPoint {
                        label: "X".to_string(),
                        coords: [0.5, 0.0, 0.5],
                        npoints: 0,
                    },
                ],
                bands_num_points: 100,
            },
            disentanglement: None,
            pw2wannier90: None,
            project_id: None,
            scf_calc_id: None,
            source_metadata: Some(WannierSourceMetadata {
                cell_representation: Some("primitive_spglib".to_string()),
                electron_count: Some(10.0),
                nspin: Some(1),
                noncolin: Some(false),
                lspinorb: Some(false),
                lda_plus_u: Some(false),
                vdw_corr: Some("none".to_string()),
            }),
            guiding_centres: true,
            use_ws_distance: true,
            write_hr: true,
            write_xyz: true,
            bands_plot: true,
            conv_window: 5,
            conv_tol: 1.0e-10,
            num_iter: 100,
        };

        let error = validate_wannier_config(&config).unwrap_err();
        assert!(error.contains("too small for this scalar source"));
    }

    fn sample_win_content(use_ws_distance: bool) -> String {
        format!(
            r#"
num_wann = 1
num_bands = 1
use_ws_distance = {}
begin unit_cell_cart
angstrom
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 5.0
end unit_cell_cart
"#,
            if use_ws_distance { "true" } else { "false" }
        )
    }

    fn sample_hr_content() -> &'static str {
        r#"
written by qcortado test
1
5
1 1 1 1 1
0 0 0 1 1 0.0 0.0
1 0 0 1 1 -1.0 0.0
-1 0 0 1 1 -1.0 0.0
0 1 0 1 1 -1.0 0.0
0 -1 0 1 1 -1.0 0.0
"#
    }

    fn sample_wsvec_content() -> &'static str {
        r#"
0 0 0 1 1
1
0 0 0
1 0 0 1 1
1
0 0 0
-1 0 0 1 1
1
0 0 0
0 1 0 1 1
1
0 0 0
0 -1 0 1 1
1
0 0 0
"#
    }

    fn unique_temp_dir(label: &str) -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "qcortado_{}_{}_{}",
            label,
            std::process::id(),
            nanos
        ))
    }

    #[test]
    fn parse_wannier_hamiltonian_reconstructs_square_lattice_dispersion() {
        let hamiltonian = parse_wannier_hamiltonian(
            "toy",
            &sample_win_content(true),
            sample_hr_content(),
            Some(sample_wsvec_content()),
        )
        .unwrap();

        let gamma = hamiltonian.eigenvalues_at([0.0, 0.0, 0.0]).unwrap();
        let m = hamiltonian.eigenvalues_at([0.5, 0.5, 0.0]).unwrap();
        let x = hamiltonian.eigenvalues_at([0.5, 0.0, 0.0]).unwrap();

        assert_eq!(gamma.len(), 1);
        assert!((gamma[0] + 4.0).abs() < 1.0e-9);
        assert!((m[0] - 4.0).abs() < 1.0e-9);
        assert!(x[0].abs() < 1.0e-9);
    }

    #[test]
    fn collect_wannier_artifacts_includes_wsvec() {
        let temp_dir = unique_temp_dir("wannier_artifacts");
        std::fs::create_dir_all(&temp_dir).unwrap();
        let seedname = "toy";
        std::fs::write(
            temp_dir.join(format!("{}_wsvec.dat", seedname)),
            sample_wsvec_content(),
        )
        .unwrap();
        std::fs::write(
            temp_dir.join(format!("{}_hr.dat", seedname)),
            sample_hr_content(),
        )
        .unwrap();

        let artifacts = collect_wannier_artifacts(&temp_dir, seedname);
        assert!(artifacts
            .iter()
            .any(|entry| entry.file_name == "toy_wsvec.dat"));

        let _ = std::fs::remove_dir_all(temp_dir);
    }

    #[test]
    fn export_ludwig_bundle_writes_metadata_and_band_grid() {
        let hamiltonian = parse_wannier_hamiltonian(
            "toy",
            &sample_win_content(true),
            sample_hr_content(),
            Some(sample_wsvec_content()),
        )
        .unwrap();
        let root_dir = unique_temp_dir("ludwig_export_root");
        std::fs::create_dir_all(&root_dir).unwrap();

        let provenance_dir = unique_temp_dir("ludwig_export_provenance");
        std::fs::create_dir_all(&provenance_dir).unwrap();
        let win_path = provenance_dir.join("toy.win");
        let hr_path = provenance_dir.join("toy_hr.dat");
        let wsvec_path = provenance_dir.join("toy_wsvec.dat");
        let wout_path = provenance_dir.join("toy.wout");
        std::fs::write(&win_path, sample_win_content(true)).unwrap();
        std::fs::write(&hr_path, sample_hr_content()).unwrap();
        std::fs::write(&wsvec_path, sample_wsvec_content()).unwrap();
        std::fs::write(&wout_path, "dummy").unwrap();

        let config = LudwigExportConfig {
            project_id: "project".to_string(),
            calculation_id: "calc123".to_string(),
            destination_root: root_dir.display().to_string(),
            mode: LudwigExportMode::Quasi2dFixedSlice,
            in_plane_axes: [0, 1],
            slice_axis: 2,
            slice_coordinate: 0.25,
            nkx: 4,
            nky: 3,
            chemical_potential_ev: None,
        };
        let result = export_ludwig_bundle(
            &hamiltonian,
            &config,
            1.5,
            &[
                win_path.clone(),
                hr_path.clone(),
                wsvec_path.clone(),
                wout_path.clone(),
            ],
        )
        .unwrap();

        let metadata =
            std::fs::read_to_string(PathBuf::from(&result.bundle_path).join("metadata.json"))
                .unwrap();
        let bands =
            std::fs::read_to_string(PathBuf::from(&result.bundle_path).join("bands.csv")).unwrap();

        assert!(metadata.contains("\"slice_coordinate\": 0.25"));
        assert!(metadata.contains("\"band_count\": 1"));
        assert!(bands.lines().next().unwrap().contains("band_1"));
        assert!(bands.contains("-5.500000000000"));
        assert!(PathBuf::from(&result.bundle_path)
            .join("toy_wsvec.dat")
            .exists());

        let _ = std::fs::remove_dir_all(root_dir);
        let _ = std::fs::remove_dir_all(provenance_dir);
    }
}
