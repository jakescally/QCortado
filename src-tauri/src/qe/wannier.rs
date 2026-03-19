//! Wannier90 input generation and result parsing helpers.

use super::bands::{add_symmetry_markers, parse_bands_gnu, BandData, KPathPoint};
use super::types::{
    CalculationType, KPoint, KPoints, Occupations, PositionUnits, QECalculation, SmearingType,
    StartingPotential,
};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fmt::Write;
use std::path::Path;

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
    let source_nspin = source.and_then(|value| value.nspin).unwrap_or(config.base_calculation.system.nspin);
    let source_noncolin =
        source.and_then(|value| value.noncolin).unwrap_or(config.base_calculation.system.noncolin);
    let source_lspinorb =
        source.and_then(|value| value.lspinorb).unwrap_or(config.base_calculation.system.lspinorb);
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
            "Scalar Wannier v1 does not certify vdW-corrected source calculations."
                .to_string(),
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
                return Err("Frozen window must lie inside the outer disentanglement window.".to_string());
            }
        }
        if let (Some(win_max), Some(froz_max)) =
            (disentanglement.dis_win_max, disentanglement.dis_froz_max)
        {
            if froz_max > win_max {
                return Err("Frozen window must lie inside the outer disentanglement window.".to_string());
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
                if nscf_calc.system.degauss.map(|value| value <= 0.0).unwrap_or(true) {
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
            } else if nscf_calc.system.degauss.map(|value| value <= 0.0).unwrap_or(true) {
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
                .ok_or_else(|| "Element-targeted Wannier projections require a symbol.".to_string())?;
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
    let mut values: Vec<u32> = exclude_bands.iter().copied().filter(|value| *value > 0).collect();
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

pub fn generate_wannier90_win(config: &WannierCalculationConfig, kpoints: &[KPoint]) -> Result<String, String> {
    validate_wannier_config(config)?;
    let system = &config.base_calculation.system;
    let cell = system
        .cell_parameters
        .as_ref()
        .ok_or_else(|| "Wannier v1 requires explicit cell_parameters in the base calculation.".to_string())?;
    let cell_units = system
        .cell_units
        .as_ref()
        .unwrap_or(&PositionUnits::Angstrom)
        .as_str()
        .to_string();

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
    writeln!(output, "{}", cell_units).unwrap();
    for row in cell {
        writeln!(
            output,
            "{:16.10} {:16.10} {:16.10}",
            row[0], row[1], row[2]
        )
        .unwrap();
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
        if pw2_config.write_amn { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  write_mmn = .{}.,",
        if pw2_config.write_mmn { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  write_spn = .{}.,",
        if pw2_config.write_spn { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  write_unk = .{}.,",
        if pw2_config.write_unk { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  write_dmn = .{}.,",
        if pw2_config.write_dmn { "true" } else { "false" }
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
        if pw2_config.scdm_proj { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  atom_proj = .{}.,",
        if pw2_config.atom_proj { "true" } else { "false" }
    )
    .unwrap();
    writeln!(output, "/").unwrap();
    output
}

pub fn parse_wannier_wout(content: &str) -> Result<(Vec<WannierSpread>, WannierConvergenceData), String> {
    let spread_re = Regex::new(
        r"WF centre and spread\s+(\d+)\s+\(\s*([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*\)\s*([-\d.Ee+]+)",
    )
    .map_err(|e| format!("Failed to compile Wannier spread regex: {}", e))?;
    let omega_re = Regex::new(
        r"Omega\s+(I|D|OD|Total)\s*=\s*([-\d.Ee+]+)",
    )
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
                convergence.iterations = caps
                    .get(1)
                    .and_then(|m| m.as_str().parse::<u32>().ok());
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
        return Err(format!("Wannier output file not found: {}", wout_path.display()));
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qe::types::{
        Atom, AtomicSpecies, BravaisLattice, CalculationType, KPoints, Occupations,
        PositionUnits, QESystem, QECalculation, SmearingType,
    };

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
        assert!(notes.iter().any(|note| note.contains("tetrahedra occupations")));
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
}
