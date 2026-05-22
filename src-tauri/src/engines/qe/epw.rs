//! EPW (`epw.x`) helpers: config validation, input generation, and lightweight output parsing.

use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write;
use std::path::{Path, PathBuf};

use super::wannier::WannierArtifact;

pub const EPW_SCHEMA_VERSION: u32 = 2;

const MAX_PARSED_TABLE_BYTES: u64 = 25 * 1024 * 1024;
const MAX_PARSED_TABLE_ROWS: usize = 50_000;

fn default_prefix() -> String {
    "qcortado_scf".to_string()
}

fn default_outdir() -> String {
    "./tmp".to_string()
}

fn default_phonon_dir() -> String {
    "./save".to_string()
}

fn default_wannier_dir() -> String {
    "./wannier".to_string()
}

fn default_k_mesh() -> [u32; 3] {
    [24, 24, 24]
}

fn default_q_mesh() -> [u32; 3] {
    [6, 6, 6]
}

fn default_true() -> bool {
    true
}

fn default_wannierize() -> bool {
    true
}

fn default_epw_goal_coupling() -> bool {
    true
}

fn default_epw_goal_linewidth() -> bool {
    true
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwGoalSelection {
    #[serde(default = "default_epw_goal_coupling")]
    pub coupling: bool,
    #[serde(default = "default_epw_goal_linewidth")]
    pub phonon_linewidth_a2f: bool,
    #[serde(default)]
    pub electron_self_energy: bool,
    #[serde(default)]
    pub transport_mobility: bool,
    #[serde(default)]
    pub superconductivity: bool,
}

impl Default for EpwGoalSelection {
    fn default() -> Self {
        Self {
            coupling: true,
            phonon_linewidth_a2f: true,
            electron_self_energy: false,
            transport_mobility: false,
            superconductivity: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwInputConfig {
    #[serde(default = "default_prefix")]
    pub prefix: String,
    #[serde(default = "default_outdir")]
    pub outdir: String,
    #[serde(default = "default_phonon_dir")]
    pub dvscf_dir: String,
    #[serde(default = "default_wannier_dir")]
    pub wannier_dir: String,
    #[serde(default = "default_k_mesh")]
    pub k_mesh: [u32; 3],
    #[serde(default = "default_q_mesh")]
    pub q_mesh: [u32; 3],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coarse_k_mesh: Option<[u32; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fine_k_mesh: Option<[u32; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coarse_q_mesh: Option<[u32; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fine_q_mesh: Option<[u32; 3]>,
    #[serde(default = "default_true")]
    pub epbwrite: bool,
    #[serde(default)]
    pub epbread: bool,
    #[serde(default = "default_true")]
    pub epwwrite: bool,
    #[serde(default)]
    pub epwread: bool,
    #[serde(default = "default_wannierize")]
    pub wannierize: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fsthick_ev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degaussw_ev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nbndsub: Option<u32>,
}

impl Default for EpwInputConfig {
    fn default() -> Self {
        Self {
            prefix: default_prefix(),
            outdir: default_outdir(),
            dvscf_dir: default_phonon_dir(),
            wannier_dir: default_wannier_dir(),
            k_mesh: default_k_mesh(),
            q_mesh: default_q_mesh(),
            coarse_k_mesh: None,
            fine_k_mesh: None,
            coarse_q_mesh: None,
            fine_q_mesh: None,
            epbwrite: true,
            epbread: false,
            epwwrite: true,
            epwread: false,
            wannierize: default_wannierize(),
            fsthick_ev: Some(0.4),
            degaussw_ev: Some(0.02),
            nbndsub: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwRuntimeConfig {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pools: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_seconds: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub artifact_sync_mode: Option<String>,
}

/// Reserved extension payloads for future EPW feature families (e.g. superconductivity).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwExtensionsV1 {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub superconductivity: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwCalculationConfig {
    pub project_id: String,
    pub source_phonon_calc_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_wannier_calc_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_scf_calc_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mode: Option<String>,
    #[serde(default)]
    pub input: EpwInputConfig,
    #[serde(default)]
    pub runtime: EpwRuntimeConfig,
    #[serde(default)]
    pub extensions: EpwExtensionsV1,
    #[serde(default)]
    pub goals: EpwGoalSelection,
    #[serde(default)]
    pub advanced_overrides: HashMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwSourceRef {
    pub calc_id: String,
    pub calc_type: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwArtifactManifestEntry {
    pub source_calc_id: String,
    pub source_calc_type: String,
    pub rel_path: String,
    pub size_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwSourcesV1 {
    pub phonon: EpwSourceRef,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wannier: Option<EpwSourceRef>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scf: Option<EpwSourceRef>,
    #[serde(default)]
    pub manifests: Vec<EpwArtifactManifestEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwErrorRecord {
    pub code: String,
    pub message: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hint: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwResultSummaryV1 {
    #[serde(default)]
    pub completed: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elapsed_seconds: Option<f64>,
    #[serde(default)]
    pub core_metrics: HashMap<String, f64>,
    #[serde(default)]
    pub generated_outputs: Vec<WannierArtifact>,
    #[serde(default)]
    pub unknown_metrics: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub parse_partial: bool,
    #[serde(default)]
    pub notes: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub parse_coverage: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwParsedTable {
    pub file_name: String,
    pub family: String,
    pub title: String,
    #[serde(default)]
    pub column_labels: Vec<String>,
    #[serde(default)]
    pub rows: Vec<Vec<Option<f64>>>,
    #[serde(default)]
    pub skipped: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skip_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwMobilityDataset {
    pub carrier_type: String,
    pub method: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub iteration: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub converged: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_error: Option<f64>,
    #[serde(default)]
    pub component_labels: Vec<String>,
    #[serde(default)]
    pub temperature_values_k: Vec<f64>,
    #[serde(default)]
    pub fermi_values_ev: Vec<f64>,
    #[serde(default)]
    pub density_values_cm3: Vec<f64>,
    #[serde(default)]
    pub population_values: Vec<Option<f64>>,
    #[serde(default)]
    pub mobility_values: Vec<Vec<Option<f64>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwTransportData {
    #[serde(default)]
    pub mobility: Vec<EpwMobilityDataset>,
    #[serde(default)]
    pub scattering_file_notices: u32,
    #[serde(default)]
    pub notes: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwEliashbergIteration {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature_k: Option<f64>,
    pub iteration: u32,
    pub ethr: f64,
    pub znorm: f64,
    pub delta_mev: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwGapSummary {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub free_energy_mev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap_min_mev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap_max_mev: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwSuperconductivityData {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lambda: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lambda_tr: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub electron_phonon_coupling: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tc_mcmillan_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tc_allen_dynes_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tc_sisso_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub w_log_mev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bcs_gap_mev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub muc: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frequency_cutoff_ev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frequency_points: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eliashberg_converged: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eliashberg_nsiter: Option<u32>,
    #[serde(default)]
    pub temperatures_k: Vec<f64>,
    #[serde(default)]
    pub eliashberg_iterations: Vec<EpwEliashbergIteration>,
    #[serde(default)]
    pub gap_summaries: Vec<EpwGapSummary>,
    #[serde(default)]
    pub spectral_tables: Vec<EpwParsedTable>,
    #[serde(default)]
    pub notes: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwSelfEnergyMode {
    pub mode_label: String,
    pub lambda: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lambda_tr: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gamma_mev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gamma_tr_mev: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub omega_mev: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwSpectralData {
    #[serde(default)]
    pub self_energy_modes: Vec<EpwSelfEnergyMode>,
    #[serde(default)]
    pub tables: Vec<EpwParsedTable>,
    #[serde(default)]
    pub notes: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct EpwParsedResultV2 {
    pub summary: EpwResultSummaryV1,
    pub transport: Option<EpwTransportData>,
    pub superconductivity: Option<EpwSuperconductivityData>,
    pub spectral: Option<EpwSpectralData>,
    pub parsed_tables: Vec<EpwParsedTable>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwCalculationV1 {
    pub schema_version: u32,
    pub sources: EpwSourcesV1,
    pub input: EpwInputConfig,
    #[serde(default)]
    pub goals: EpwGoalSelection,
    pub runtime: EpwRuntimeConfig,
    #[serde(default)]
    pub extensions: EpwExtensionsV1,
    #[serde(default)]
    pub artifacts: Vec<WannierArtifact>,
    pub result_summary: EpwResultSummaryV1,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transport: Option<EpwTransportData>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub superconductivity: Option<EpwSuperconductivityData>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spectral: Option<EpwSpectralData>,
    #[serde(default)]
    pub parsed_tables: Vec<EpwParsedTable>,
    #[serde(default)]
    pub errors: Vec<EpwErrorRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EpwPrerequisiteValidation {
    #[serde(default)]
    pub ok: bool,
    #[serde(default)]
    pub errors: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub remediation_hints: Vec<String>,
    #[serde(default)]
    pub manifests: Vec<EpwArtifactManifestEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwInputPreviewResult {
    pub schema_version: u32,
    pub input_text: String,
    pub merged_keywords: BTreeMap<String, String>,
}

fn quote_string(value: &str) -> String {
    let escaped = value.replace('\'', "''");
    format!("'{}'", escaped)
}

fn bool_keyword(value: bool) -> String {
    if value {
        ".true.".to_string()
    } else {
        ".false.".to_string()
    }
}

fn mesh_has_zero(mesh: &[u32; 3]) -> bool {
    mesh.iter().any(|entry| *entry == 0)
}

pub fn epw_coarse_k_mesh(input: &EpwInputConfig) -> [u32; 3] {
    input.coarse_k_mesh.unwrap_or(input.k_mesh)
}

pub fn epw_fine_k_mesh(input: &EpwInputConfig) -> [u32; 3] {
    input.fine_k_mesh.unwrap_or(input.k_mesh)
}

pub fn epw_coarse_q_mesh(input: &EpwInputConfig) -> [u32; 3] {
    input.coarse_q_mesh.unwrap_or(input.q_mesh)
}

pub fn epw_fine_q_mesh(input: &EpwInputConfig) -> [u32; 3] {
    input.fine_q_mesh.unwrap_or_else(|| epw_fine_k_mesh(input))
}

fn format_decimal(value: f64) -> String {
    let formatted = format!("{:.12}", value);
    formatted
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

pub fn validate_epw_config(config: &EpwCalculationConfig) -> Result<(), String> {
    if config.project_id.trim().is_empty() {
        return Err("EPW project_id is required.".to_string());
    }
    if config.source_phonon_calc_id.trim().is_empty() {
        return Err("EPW source_phonon_calc_id is required.".to_string());
    }
    if config.input.prefix.trim().is_empty() {
        return Err("EPW input prefix cannot be empty.".to_string());
    }
    if config.input.outdir.trim().is_empty() {
        return Err("EPW input outdir cannot be empty.".to_string());
    }
    if config.input.dvscf_dir.trim().is_empty() {
        return Err("EPW input dvscf_dir cannot be empty.".to_string());
    }
    if mesh_has_zero(&config.input.k_mesh) {
        return Err("EPW k_mesh values must all be positive.".to_string());
    }
    if mesh_has_zero(&config.input.q_mesh) {
        return Err("EPW q_mesh values must all be positive.".to_string());
    }
    for (label, mesh) in [
        ("coarse_k_mesh", config.input.coarse_k_mesh),
        ("fine_k_mesh", config.input.fine_k_mesh),
        ("coarse_q_mesh", config.input.coarse_q_mesh),
        ("fine_q_mesh", config.input.fine_q_mesh),
    ] {
        if let Some(mesh) = mesh {
            if mesh_has_zero(&mesh) {
                return Err(format!("EPW {} values must all be positive.", label));
            }
        }
    }
    if let Some(value) = config.input.fsthick_ev {
        if !value.is_finite() || value <= 0.0 {
            return Err("EPW fsthick_ev must be a positive finite number.".to_string());
        }
    }
    if let Some(value) = config.input.degaussw_ev {
        if !value.is_finite() || value <= 0.0 {
            return Err("EPW degaussw_ev must be a positive finite number.".to_string());
        }
    }
    if let Some(value) = config.input.nbndsub {
        if value == 0 {
            return Err("EPW nbndsub must be positive when set.".to_string());
        }
    }
    if let Some(value) = config.runtime.pools {
        if value == 0 {
            return Err("EPW runtime pools must be positive when set.".to_string());
        }
    }
    if let Some(value) = config.runtime.max_seconds {
        if value == 0 {
            return Err("EPW runtime max_seconds must be positive when set.".to_string());
        }
    }
    Ok(())
}

pub fn build_epw_keyword_map(
    config: &EpwCalculationConfig,
) -> Result<BTreeMap<String, String>, String> {
    validate_epw_config(config)?;

    let mut keywords: BTreeMap<String, String> = BTreeMap::new();
    keywords.insert(
        "prefix".to_string(),
        quote_string(config.input.prefix.trim()),
    );
    keywords.insert(
        "outdir".to_string(),
        quote_string(config.input.outdir.trim()),
    );
    keywords.insert(
        "dvscf_dir".to_string(),
        quote_string(config.input.dvscf_dir.trim()),
    );
    let coarse_k_mesh = epw_coarse_k_mesh(&config.input);
    let fine_k_mesh = epw_fine_k_mesh(&config.input);
    let coarse_q_mesh = epw_coarse_q_mesh(&config.input);
    let fine_q_mesh = epw_fine_q_mesh(&config.input);
    keywords.insert("nk1".to_string(), coarse_k_mesh[0].to_string());
    keywords.insert("nk2".to_string(), coarse_k_mesh[1].to_string());
    keywords.insert("nk3".to_string(), coarse_k_mesh[2].to_string());
    keywords.insert("nq1".to_string(), coarse_q_mesh[0].to_string());
    keywords.insert("nq2".to_string(), coarse_q_mesh[1].to_string());
    keywords.insert("nq3".to_string(), coarse_q_mesh[2].to_string());
    keywords.insert("nkf1".to_string(), fine_k_mesh[0].to_string());
    keywords.insert("nkf2".to_string(), fine_k_mesh[1].to_string());
    keywords.insert("nkf3".to_string(), fine_k_mesh[2].to_string());
    keywords.insert("nqf1".to_string(), fine_q_mesh[0].to_string());
    keywords.insert("nqf2".to_string(), fine_q_mesh[1].to_string());
    keywords.insert("nqf3".to_string(), fine_q_mesh[2].to_string());

    let goals = &config.goals;
    let any_epw_physics = goals.coupling
        || goals.phonon_linewidth_a2f
        || goals.electron_self_energy
        || goals.transport_mobility
        || goals.superconductivity;
    if any_epw_physics {
        keywords.insert("elph".to_string(), ".true.".to_string());
    }
    if goals.coupling || goals.superconductivity {
        keywords.insert("ep_coupling".to_string(), ".true.".to_string());
    }
    if goals.phonon_linewidth_a2f {
        keywords.insert("phonselfen".to_string(), ".true.".to_string());
        keywords.insert("a2f".to_string(), ".true.".to_string());
    }
    if goals.electron_self_energy {
        keywords.insert("elecselfen".to_string(), ".true.".to_string());
        keywords.insert("specfun_el".to_string(), ".true.".to_string());
    }
    if goals.transport_mobility {
        keywords.insert("scattering".to_string(), ".true.".to_string());
        keywords.insert("scattering_serta".to_string(), ".true.".to_string());
        keywords.insert("int_mob".to_string(), ".true.".to_string());
    }
    if goals.superconductivity {
        keywords.insert("ephwrite".to_string(), ".true.".to_string());
        keywords.insert("eliashberg".to_string(), ".true.".to_string());
        keywords.insert("liso".to_string(), ".true.".to_string());
        keywords.insert("limag".to_string(), ".true.".to_string());
    }

    keywords.insert("epbwrite".to_string(), bool_keyword(config.input.epbwrite));
    keywords.insert("epbread".to_string(), bool_keyword(config.input.epbread));
    keywords.insert("epwwrite".to_string(), bool_keyword(config.input.epwwrite));
    keywords.insert("epwread".to_string(), bool_keyword(config.input.epwread));
    keywords.insert(
        "wannierize".to_string(),
        bool_keyword(config.input.wannierize),
    );

    if let Some(value) = config.input.fsthick_ev {
        keywords.insert("fsthick".to_string(), format_decimal(value));
    }
    if let Some(value) = config.input.degaussw_ev {
        keywords.insert("degaussw".to_string(), format_decimal(value));
    }
    if let Some(value) = config.input.nbndsub {
        keywords.insert("nbndsub".to_string(), value.to_string());
    }
    if let Some(value) = config.runtime.pools {
        keywords.insert("npool".to_string(), value.to_string());
    }

    for (raw_key, raw_value) in &config.advanced_overrides {
        let key = raw_key.trim().to_ascii_lowercase();
        if key.is_empty() {
            continue;
        }
        let value = raw_value.trim();
        if value.is_empty() {
            keywords.remove(&key);
            continue;
        }
        keywords.insert(key, value.to_string());
    }

    Ok(keywords)
}

pub fn render_epw_input(keyword_map: &BTreeMap<String, String>) -> String {
    let mut output = String::new();
    output.push_str("&inputepw\n");
    for (key, value) in keyword_map {
        writeln!(output, "  {} = {}", key, value).unwrap();
    }
    output.push_str("/\n");
    output
}

pub fn build_epw_input(config: &EpwCalculationConfig) -> Result<String, String> {
    let keyword_map = build_epw_keyword_map(config)?;
    Ok(render_epw_input(&keyword_map))
}

pub fn build_epw_input_preview(
    config: &EpwCalculationConfig,
) -> Result<EpwInputPreviewResult, String> {
    let merged_keywords = build_epw_keyword_map(config)?;
    let input_text = render_epw_input(&merged_keywords);
    Ok(EpwInputPreviewResult {
        schema_version: EPW_SCHEMA_VERSION,
        input_text,
        merged_keywords,
    })
}

fn should_keep_epw_artifact(lower: &str) -> bool {
    if is_epw_scratch_artifact(lower) {
        return lower == "epw.in" || lower == "epw.out" || lower == "epw.err";
    }
    lower == "epw.in"
        || lower == "epw.out"
        || lower == "epw.err"
        || lower == "run.sbatch"
        || lower == "slurm.out"
        || lower == "slurm.err"
        || lower.starts_with("epw")
        || lower.ends_with(".epw")
        || lower.ends_with(".a2f")
        || lower.ends_with(".freq")
        || lower.ends_with(".frq")
        || lower.ends_with(".fmt")
        || lower.ends_with(".dat")
        || lower.ends_with(".dos")
        || lower.ends_with(".phdos")
        || lower.ends_with(".gnu")
        || lower.ends_with(".txt")
        || lower.contains("specfun")
        || lower.contains("selfen")
        || lower.contains("scat")
        || lower.contains("mob")
}

fn collect_epw_artifacts_recursive(
    root: &Path,
    current: &Path,
    artifacts: &mut Vec<WannierArtifact>,
) {
    let Ok(entries) = std::fs::read_dir(current) else {
        return;
    };

    for entry in entries.flatten() {
        let path = entry.path();
        let Ok(meta) = entry.metadata() else {
            continue;
        };
        if meta.is_dir() {
            collect_epw_artifacts_recursive(root, &path, artifacts);
            continue;
        }
        if !meta.is_file() {
            continue;
        }
        let rel_path = path
            .strip_prefix(root)
            .unwrap_or(path.as_path())
            .to_string_lossy()
            .replace('\\', "/");
        let lower = rel_path.to_ascii_lowercase();
        if should_keep_epw_artifact(&lower) {
            artifacts.push(WannierArtifact {
                file_name: rel_path,
                size_bytes: meta.len(),
            });
        }
    }
}

pub fn collect_epw_artifacts(work_path: &Path) -> Vec<WannierArtifact> {
    let mut artifacts: Vec<WannierArtifact> = Vec::new();
    collect_epw_artifacts_recursive(work_path, work_path, &mut artifacts);
    artifacts.sort_by(|a, b| a.file_name.cmp(&b.file_name));
    artifacts
}

fn is_epw_scratch_artifact(file_name: &str) -> bool {
    let lower = file_name.to_ascii_lowercase();
    lower.contains(".save/")
        || lower.starts_with("save/")
        || lower.starts_with("phonon/")
        || lower.starts_with("wannier/")
        || lower.contains("/_ph0/")
        || lower.contains("/wfc")
        || lower.ends_with("charge-density.dat")
        || lower.ends_with("data-file-schema.xml")
}

fn parse_fortran_float(raw: &str) -> Option<f64> {
    let cleaned = raw
        .trim()
        .trim_matches(',')
        .trim_matches(';')
        .replace(['D', 'd'], "E");
    cleaned
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
}

fn numeric_tokens(line: &str) -> Vec<f64> {
    line.split_whitespace()
        .filter_map(parse_fortran_float)
        .collect()
}

fn push_unique_number(values: &mut Vec<f64>, value: f64) {
    if values.iter().any(|existing| {
        let scale = existing.abs().max(value.abs()).max(1.0);
        (existing - value).abs() <= scale * 1.0e-9
    }) {
        return;
    }
    values.push(value);
}

fn push_unique_string(values: &mut Vec<String>, value: &str) {
    let trimmed = value.trim();
    if trimmed.is_empty() || values.iter().any(|entry| entry == trimmed) {
        return;
    }
    values.push(trimmed.to_string());
}

#[derive(Debug, Clone)]
struct ParsedMobilityRow {
    temperature_k: f64,
    fermi_ev: f64,
    density_cm3: f64,
    population: Option<f64>,
    values: HashMap<String, f64>,
}

fn detect_axis_component(line: &str, fallback: &str) -> String {
    let lower = line.to_ascii_lowercase();
    if lower.contains("x-axis") {
        "xx".to_string()
    } else if lower.contains("y-axis") {
        "yy".to_string()
    } else if lower.contains("z-axis") {
        "zz".to_string()
    } else if lower.contains("avg") {
        "avg".to_string()
    } else {
        fallback.to_string()
    }
}

fn finalize_mobility_dataset(
    carrier_type: &str,
    method: &str,
    iteration: Option<u32>,
    rows: Vec<ParsedMobilityRow>,
    converged: Option<bool>,
    max_error: Option<f64>,
) -> Option<EpwMobilityDataset> {
    if rows.is_empty() {
        return None;
    }

    let mut component_labels: Vec<String> = Vec::new();
    for preferred in ["xx", "yy", "zz", "avg", "xy", "xz", "yx", "yz", "zx", "zy"] {
        if rows.iter().any(|row| row.values.contains_key(preferred)) {
            component_labels.push(preferred.to_string());
        }
    }
    let mut extra_labels = rows
        .iter()
        .flat_map(|row| row.values.keys().cloned())
        .filter(|label| !component_labels.iter().any(|existing| existing == label))
        .collect::<Vec<_>>();
    extra_labels.sort();
    component_labels.extend(extra_labels);
    if component_labels.is_empty() {
        return None;
    }

    let temperature_values_k = rows.iter().map(|row| row.temperature_k).collect::<Vec<_>>();
    let fermi_values_ev = rows.iter().map(|row| row.fermi_ev).collect::<Vec<_>>();
    let density_values_cm3 = rows.iter().map(|row| row.density_cm3).collect::<Vec<_>>();
    let population_values = rows.iter().map(|row| row.population).collect::<Vec<_>>();
    let mobility_values = component_labels
        .iter()
        .map(|component| {
            rows.iter()
                .map(|row| row.values.get(component).copied())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    Some(EpwMobilityDataset {
        carrier_type: carrier_type.to_string(),
        method: method.to_string(),
        iteration,
        converged,
        max_error,
        component_labels,
        temperature_values_k,
        fermi_values_ev,
        density_values_cm3,
        population_values,
        mobility_values,
    })
}

fn parse_axis_mobility_block(
    lines: &[&str],
    start_index: usize,
    carrier_type: &str,
    method: &str,
    iteration: Option<u32>,
    converged: Option<bool>,
    max_error: Option<f64>,
) -> (Option<EpwMobilityDataset>, usize) {
    let mut rows: Vec<ParsedMobilityRow> = Vec::new();
    let mut index = start_index + 1;

    while index < lines.len() {
        let line = lines[index];
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.chars().all(|ch| ch == '=') {
            index += 1;
            continue;
        }
        let lower = trimmed.to_ascii_lowercase();
        if lower.contains("temp [k]")
            || lower.contains("total time")
            || lower.contains("writing scattering")
        {
            break;
        }
        let numbers = numeric_tokens(trimmed);
        if numbers.len() < 4 {
            if rows.is_empty() {
                index += 1;
                continue;
            }
            break;
        }

        let mut values = HashMap::new();
        values.insert(detect_axis_component(trimmed, "xx"), numbers[3]);
        index += 1;
        while index < lines.len() {
            let continuation = lines[index].trim();
            if continuation.is_empty() || continuation.chars().all(|ch| ch == '=') {
                index += 1;
                continue;
            }
            if continuation
                .split_whitespace()
                .next()
                .and_then(parse_fortran_float)
                .is_some()
                && numeric_tokens(continuation).len() >= 4
            {
                break;
            }
            let continuation_numbers = numeric_tokens(continuation);
            if continuation_numbers.is_empty() {
                break;
            }
            let component = detect_axis_component(continuation, "value");
            values.insert(component, continuation_numbers[0]);
            index += 1;
            if values.len() >= 4 {
                break;
            }
        }

        rows.push(ParsedMobilityRow {
            temperature_k: numbers[0],
            fermi_ev: numbers[1],
            density_cm3: numbers[2],
            population: None,
            values,
        });
    }

    (
        finalize_mobility_dataset(carrier_type, method, iteration, rows, converged, max_error),
        index,
    )
}

fn parse_tensor_mobility_block(
    lines: &[&str],
    start_index: usize,
    carrier_type: &str,
    method: &str,
    iteration: Option<u32>,
    converged: Option<bool>,
    max_error: Option<f64>,
) -> (Option<EpwMobilityDataset>, usize) {
    let mut rows: Vec<ParsedMobilityRow> = Vec::new();
    let mut index = start_index + 1;

    while index < lines.len() {
        let line = lines[index].trim();
        if line.is_empty() || line.chars().all(|ch| ch == '=') || line.starts_with("[K]") {
            index += 1;
            continue;
        }
        let lower = line.to_ascii_lowercase();
        if lower.contains("iteration number")
            || lower.contains("max error")
            || lower.contains("unfolding")
            || lower.contains("electron-phonon")
            || lower.contains("total program")
        {
            break;
        }
        let first_numbers = numeric_tokens(line);
        if first_numbers.len() < 7 {
            if rows.is_empty() {
                index += 1;
                continue;
            }
            break;
        }

        let mut tensor_rows = vec![[
            first_numbers[first_numbers.len() - 3],
            first_numbers[first_numbers.len() - 2],
            first_numbers[first_numbers.len() - 1],
        ]];
        let mut population = Some(first_numbers[3]);
        let temperature_k = first_numbers[0];
        let fermi_ev = first_numbers[1];
        let density_cm3 = first_numbers[2];
        index += 1;

        while tensor_rows.len() < 3 && index < lines.len() {
            let continuation = lines[index].trim();
            if continuation.is_empty() {
                index += 1;
                continue;
            }
            let continuation_numbers = numeric_tokens(continuation);
            if continuation_numbers.len() < 4 {
                break;
            }
            tensor_rows.push([
                continuation_numbers[continuation_numbers.len() - 3],
                continuation_numbers[continuation_numbers.len() - 2],
                continuation_numbers[continuation_numbers.len() - 1],
            ]);
            if population.is_none() {
                population = Some(continuation_numbers[0]);
            }
            index += 1;
        }

        if tensor_rows.len() == 3 {
            let mut values = HashMap::new();
            let labels = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"];
            for (label, value) in labels.iter().zip(tensor_rows.iter().flatten()) {
                values.insert((*label).to_string(), *value);
            }
            rows.push(ParsedMobilityRow {
                temperature_k,
                fermi_ev,
                density_cm3,
                population,
                values,
            });
        }
    }

    (
        finalize_mobility_dataset(carrier_type, method, iteration, rows, converged, max_error),
        index,
    )
}

fn parse_epw_transport_data(output: &str) -> Option<EpwTransportData> {
    let lines = output.lines().collect::<Vec<_>>();
    let mut mobility: Vec<EpwMobilityDataset> = Vec::new();
    let mut notes = Vec::new();
    let mut warnings = Vec::new();
    let mut scattering_file_notices = 0_u32;
    let mut method = "generic".to_string();
    let mut iteration: Option<u32> = None;
    let mut index = 0;
    let iteration_re = Regex::new(r"(?i)Iteration number:\s*(\d+)").expect("valid iteration regex");
    let max_error_re = Regex::new(r"([-\d.Ee+Dd+]+)\s+Max error").expect("valid max error regex");

    while index < lines.len() {
        let line = lines[index];
        let lower = line.to_ascii_lowercase();
        if lower.contains("bte in the self-energy relaxation time approximation") {
            method = "SERTA".to_string();
        } else if lower.contains("start solving iterative boltzmann transport equation") {
            method = "IBTE".to_string();
            iteration = None;
        } else if let Some(captures) = iteration_re.captures(line) {
            iteration = captures
                .get(1)
                .and_then(|entry| entry.as_str().parse::<u32>().ok());
        } else if lower.contains("writing scattering rate to file") {
            scattering_file_notices += 1;
        } else if lower.contains("mobility are sorted") {
            push_unique_string(&mut warnings, line.trim());
        } else if lower.contains("no intermediate mobility") {
            push_unique_string(&mut notes, line.trim());
        } else if let Some(captures) = max_error_re.captures(line) {
            if let Some(value) = captures
                .get(1)
                .and_then(|entry| parse_fortran_float(entry.as_str()))
            {
                if let Some(last) = mobility.last_mut() {
                    last.max_error = Some(value);
                }
            }
        }

        if lower.contains("temp")
            && lower.contains("fermi")
            && lower.contains("mobility")
            && (lower.contains("hole") || lower.contains("elec"))
        {
            let carrier_type = if lower.contains("hole") {
                "hole"
            } else {
                "electron"
            };
            let is_tensor = lower.contains("population") || lower.contains("drift");
            let (dataset, next_index) = if is_tensor {
                parse_tensor_mobility_block(
                    &lines,
                    index,
                    carrier_type,
                    &method,
                    iteration,
                    None,
                    None,
                )
            } else {
                parse_axis_mobility_block(
                    &lines,
                    index,
                    carrier_type,
                    &method,
                    iteration,
                    None,
                    None,
                )
            };
            if let Some(dataset) = dataset {
                mobility.push(dataset);
            }
            index = next_index;
            continue;
        }

        if lower.contains("the iteration reached the maximum but did not converge") {
            if let Some(last) = mobility.last_mut() {
                last.converged = Some(false);
            }
            push_unique_string(&mut warnings, line.trim());
        }

        index += 1;
    }

    if mobility.is_empty()
        && scattering_file_notices == 0
        && notes.is_empty()
        && warnings.is_empty()
    {
        return None;
    }

    Some(EpwTransportData {
        mobility,
        scattering_file_notices,
        notes,
        warnings,
    })
}

fn parse_first_capture(output: &str, pattern: &str) -> Option<f64> {
    let regex = Regex::new(pattern).ok()?;
    regex
        .captures_iter(output)
        .last()
        .and_then(|captures| captures.get(1))
        .and_then(|entry| parse_fortran_float(entry.as_str()))
}

fn parse_epw_superconductivity_data(
    output: &str,
    parsed_tables: &[EpwParsedTable],
) -> Option<EpwSuperconductivityData> {
    let mut data = EpwSuperconductivityData::default();
    data.electron_phonon_coupling = parse_first_capture(
        output,
        r"(?i)Electron-phonon coupling strength\s*=\s*([-\d.Ee+Dd+]+)",
    );
    data.lambda =
        parse_first_capture(output, r"(?i)\blambda\s*[:=]\s*([-\d.Ee+Dd+]+)").or_else(|| {
            parse_first_capture(output, r"(?i)lambda___\(\s*tot\s*\)\s*=\s*([-\d.Ee+Dd+]+)")
        });
    data.lambda_tr = parse_first_capture(output, r"(?i)\blambda_tr\s*[:=]\s*([-\d.Ee+Dd+]+)")
        .or_else(|| {
            parse_first_capture(output, r"(?i)lambda_tr\(\s*tot\s*\)\s*=\s*([-\d.Ee+Dd+]+)")
        });
    data.tc_mcmillan_k = parse_first_capture(
        output,
        r"(?i)Estimated Tc using McMillan expression\s*=\s*([-\d.Ee+Dd+]+)\s*K",
    );
    data.tc_allen_dynes_k = parse_first_capture(
        output,
        r"(?i)Estimated Tc using Allen-Dynes(?: modified McMillan expression)?\s*=\s*([-\d.Ee+Dd+]+)\s*K",
    );
    data.tc_sisso_k = parse_first_capture(
        output,
        r"(?i)Estimated Tc using SISSO machine learning model\s*=\s*([-\d.Ee+Dd+]+)\s*K",
    );
    data.w_log_mev =
        parse_first_capture(output, r"(?i)Estimated w_log\s*=\s*([-\d.Ee+Dd+]+)\s*meV");
    data.bcs_gap_mev = parse_first_capture(
        output,
        r"(?i)Estimated BCS superconducting gap using McMillan Tc\s*=\s*([-\d.Ee+Dd+]+)\s*meV",
    );
    data.muc = parse_first_capture(output, r"(?i)\bmuc\s*=\s*([-\d.Ee+Dd+]+)")
        .or_else(|| parse_first_capture(output, r"(?i)for muc\s*=\s*([-\d.Ee+Dd+]+)"));
    data.frequency_cutoff_ev =
        parse_first_capture(output, r"(?i)Cutoff frequency wscut\s*=\s*([-\d.Ee+Dd+]+)");
    data.frequency_points =
        Regex::new(r"(?i)Total number of frequency points nsiw\(\s*\d+\)\s*=\s*(\d+)")
            .ok()
            .and_then(|regex| regex.captures_iter(output).last())
            .and_then(|captures| captures.get(1))
            .and_then(|entry| entry.as_str().parse::<u32>().ok())
            .or_else(|| {
                Regex::new(r"(?i)Actual number of frequency points\s*\(\s*\d+\)\s*=\s*(\d+)")
                    .ok()
                    .and_then(|regex| regex.captures_iter(output).last())
                    .and_then(|captures| captures.get(1))
                    .and_then(|entry| entry.as_str().parse::<u32>().ok())
            });

    let temp_re =
        Regex::new(r"(?i)temp\(\s*\d+\)\s*=\s*([-\d.Ee+Dd+]+)\s*K").expect("valid EPW temp regex");
    let itemp_re = Regex::new(r"(?i)Temp\s*\(itemp\s*=\s*\d+\)\s*=\s*([-\d.Ee+Dd+]+)\s*K(?:\s+Free energy\s*=\s*([-\d.Ee+Dd+]+)\s*meV)?")
        .expect("valid EPW itemp regex");
    let iter_re =
        Regex::new(r"^\s*(\d+)\s+([-\d.Ee+Dd+]+)\s+([-\d.Ee+Dd+]+)\s+([-\d.Ee+Dd+]+)\s*$")
            .expect("valid Eliashberg iteration regex");
    let gap_re = Regex::new(r"(?i)Min\.\s*/\s*Max\.\s*values of superconducting gap\s*=\s*([-\d.Ee+Dd+]+)\s+([-\d.Ee+Dd+]+)\s*meV")
        .expect("valid gap regex");
    let convergence_re = Regex::new(r"(?i)Convergence was reached in nsiter\s*=\s*(\d+)")
        .expect("valid convergence regex");
    let mut current_temperature: Option<f64> = None;
    let mut active_gap_summary: Option<usize> = None;

    for line in output.lines() {
        if let Some(captures) = temp_re.captures(line) {
            current_temperature = captures
                .get(1)
                .and_then(|entry| parse_fortran_float(entry.as_str()));
            if let Some(temp) = current_temperature {
                push_unique_number(&mut data.temperatures_k, temp);
            }
        }
        if let Some(captures) = itemp_re.captures(line) {
            current_temperature = captures
                .get(1)
                .and_then(|entry| parse_fortran_float(entry.as_str()));
            if let Some(temp) = current_temperature {
                push_unique_number(&mut data.temperatures_k, temp);
            }
            data.gap_summaries.push(EpwGapSummary {
                temperature_k: current_temperature,
                free_energy_mev: captures
                    .get(2)
                    .and_then(|entry| parse_fortran_float(entry.as_str())),
                gap_min_mev: None,
                gap_max_mev: None,
            });
            active_gap_summary = Some(data.gap_summaries.len() - 1);
        }
        if let Some(captures) = iter_re.captures(line) {
            if let (Some(iteration), Some(ethr), Some(znorm), Some(delta_mev)) = (
                captures
                    .get(1)
                    .and_then(|entry| entry.as_str().parse::<u32>().ok()),
                captures
                    .get(2)
                    .and_then(|entry| parse_fortran_float(entry.as_str())),
                captures
                    .get(3)
                    .and_then(|entry| parse_fortran_float(entry.as_str())),
                captures
                    .get(4)
                    .and_then(|entry| parse_fortran_float(entry.as_str())),
            ) {
                data.eliashberg_iterations.push(EpwEliashbergIteration {
                    temperature_k: current_temperature,
                    iteration,
                    ethr,
                    znorm,
                    delta_mev,
                });
            }
        }
        if let Some(captures) = gap_re.captures(line) {
            let gap_min_mev = captures
                .get(1)
                .and_then(|entry| parse_fortran_float(entry.as_str()));
            let gap_max_mev = captures
                .get(2)
                .and_then(|entry| parse_fortran_float(entry.as_str()));
            if let Some(index) = active_gap_summary {
                if let Some(summary) = data.gap_summaries.get_mut(index) {
                    summary.gap_min_mev = gap_min_mev;
                    summary.gap_max_mev = gap_max_mev;
                }
            } else {
                data.gap_summaries.push(EpwGapSummary {
                    temperature_k: current_temperature,
                    free_energy_mev: None,
                    gap_min_mev,
                    gap_max_mev,
                });
            }
        }
        if let Some(captures) = convergence_re.captures(line) {
            data.eliashberg_converged = Some(true);
            data.eliashberg_nsiter = captures
                .get(1)
                .and_then(|entry| entry.as_str().parse::<u32>().ok());
        }
        let lower = line.to_ascii_lowercase();
        if lower.contains("eliashberg") || lower.contains("a2f file") {
            push_unique_string(&mut data.notes, line.trim());
        }
        if lower.contains("warning") && (lower.contains("tc") || lower.contains("eliashberg")) {
            push_unique_string(&mut data.warnings, line.trim());
        }
    }

    data.spectral_tables = parsed_tables
        .iter()
        .filter(|table| matches!(table.family.as_str(), "a2f" | "phdos" | "dos"))
        .cloned()
        .collect();

    let has_data = data.lambda.is_some()
        || data.lambda_tr.is_some()
        || data.electron_phonon_coupling.is_some()
        || data.tc_mcmillan_k.is_some()
        || data.tc_allen_dynes_k.is_some()
        || data.tc_sisso_k.is_some()
        || !data.eliashberg_iterations.is_empty()
        || !data.gap_summaries.is_empty()
        || !data.spectral_tables.is_empty();
    if has_data {
        Some(data)
    } else {
        None
    }
}

fn parse_epw_self_energy_modes(output: &str) -> Vec<EpwSelfEnergyMode> {
    let lambda_re = Regex::new(
        r"(?i)lambda___\(\s*([^)]+?)\s*\)\s*=\s*([-\d.Ee+Dd+]+).*?gamma___=\s*([-\d.Ee+Dd+]+)\s*meV\s+omega=\s*([-\d.Ee+Dd+]+)\s*meV",
    )
    .expect("valid lambda mode regex");
    let lambda_tr_re = Regex::new(
        r"(?i)lambda_tr\(\s*([^)]+?)\s*\)\s*=\s*([-\d.Ee+Dd+]+).*?gamma_tr=\s*([-\d.Ee+Dd+]+)\s*meV",
    )
    .expect("valid lambda_tr mode regex");
    let mut modes: Vec<EpwSelfEnergyMode> = Vec::new();

    for line in output.lines() {
        if let Some(captures) = lambda_re.captures(line) {
            let label = captures
                .get(1)
                .map(|entry| entry.as_str().trim().to_string())
                .unwrap_or_default();
            if let Some(lambda) = captures
                .get(2)
                .and_then(|entry| parse_fortran_float(entry.as_str()))
            {
                modes.push(EpwSelfEnergyMode {
                    mode_label: label,
                    lambda,
                    lambda_tr: None,
                    gamma_mev: captures
                        .get(3)
                        .and_then(|entry| parse_fortran_float(entry.as_str())),
                    gamma_tr_mev: None,
                    omega_mev: captures
                        .get(4)
                        .and_then(|entry| parse_fortran_float(entry.as_str())),
                });
            }
        } else if let Some(captures) = lambda_tr_re.captures(line) {
            let label = captures
                .get(1)
                .map(|entry| entry.as_str().trim().to_string())
                .unwrap_or_default();
            if let Some(mode) = modes.iter_mut().rev().find(|mode| mode.mode_label == label) {
                mode.lambda_tr = captures
                    .get(2)
                    .and_then(|entry| parse_fortran_float(entry.as_str()));
                mode.gamma_tr_mev = captures
                    .get(3)
                    .and_then(|entry| parse_fortran_float(entry.as_str()));
            }
        }
    }

    modes
}

fn infer_epw_table_family(file_name: &str) -> String {
    let lower = file_name.to_ascii_lowercase();
    if lower.ends_with(".a2f") || lower.contains("a2f") {
        "a2f".to_string()
    } else if lower.contains("phdos") || lower.ends_with(".phdos") {
        "phdos".to_string()
    } else if lower.contains("specfun") {
        "spectral_function".to_string()
    } else if lower.contains("selfen") {
        "self_energy".to_string()
    } else if lower.contains("scat") || lower.contains("scatter") {
        "scattering".to_string()
    } else if lower.contains("mob") {
        "mobility".to_string()
    } else if lower.ends_with(".dos") || lower.contains(".dos") {
        "dos".to_string()
    } else if lower.ends_with(".freq") || lower.ends_with(".frq") {
        "frequency".to_string()
    } else {
        "generic".to_string()
    }
}

fn is_chartable_epw_file(file_name: &str) -> bool {
    let lower = file_name.to_ascii_lowercase();
    if is_epw_scratch_artifact(&lower)
        || lower.ends_with(".xml")
        || lower.ends_with("wfc.dat")
        || lower.contains("/wfc")
        || lower.contains("charge-density")
    {
        return false;
    }
    lower.ends_with(".a2f")
        || lower.ends_with(".dos")
        || lower.ends_with(".phdos")
        || lower.ends_with(".freq")
        || lower.ends_with(".frq")
        || lower.ends_with(".fmt")
        || lower.ends_with(".gnu")
        || lower.ends_with(".txt")
        || lower.ends_with(".dat") && infer_epw_table_family(&lower) != "generic"
        || lower.contains("specfun")
        || lower.contains("selfen")
        || lower.contains("scat")
        || lower.contains("mob")
}

fn labels_from_header(header: Option<&str>, count: usize) -> Vec<String> {
    if let Some(raw) = header {
        let mut labels = raw
            .trim()
            .trim_start_matches('#')
            .trim_start_matches('!')
            .split_whitespace()
            .map(|entry| {
                entry.trim_matches(|ch: char| ch == '[' || ch == ']' || ch == ',' || ch == ':')
            })
            .filter(|entry| !entry.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        if labels.len() == count {
            return labels;
        }
        if labels.len() > count {
            labels.truncate(count);
            return labels;
        }
    }
    (0..count).map(|index| format!("c{}", index + 1)).collect()
}

fn skipped_table(file_name: &str, size_bytes: u64, reason: String) -> EpwParsedTable {
    EpwParsedTable {
        file_name: file_name.to_string(),
        family: infer_epw_table_family(file_name),
        title: format!("{} ({})", file_name, format!("{} bytes", size_bytes)),
        column_labels: Vec::new(),
        rows: Vec::new(),
        skipped: true,
        skip_reason: Some(reason),
    }
}

fn parse_numeric_table_file(
    work_path: &Path,
    artifact: &WannierArtifact,
) -> Option<EpwParsedTable> {
    if !is_chartable_epw_file(&artifact.file_name) {
        return None;
    }
    if artifact.size_bytes > MAX_PARSED_TABLE_BYTES {
        return Some(skipped_table(
            &artifact.file_name,
            artifact.size_bytes,
            format!(
                "File is larger than the {} MB parser preview limit.",
                MAX_PARSED_TABLE_BYTES / 1024 / 1024
            ),
        ));
    }
    let path = artifact
        .file_name
        .split('/')
        .fold(PathBuf::from(work_path), |acc, part| acc.join(part));
    let Ok(content) = std::fs::read_to_string(&path) else {
        return Some(skipped_table(
            &artifact.file_name,
            artifact.size_bytes,
            "File could not be read as UTF-8 text.".to_string(),
        ));
    };

    let mut header: Option<String> = None;
    let mut rows: Vec<Vec<Option<f64>>> = Vec::new();
    let mut column_count: Option<usize> = None;
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') || trimmed.starts_with('!') {
            header = Some(trimmed.to_string());
            continue;
        }
        let tokens = trimmed.split_whitespace().collect::<Vec<_>>();
        let parsed = tokens
            .iter()
            .map(|token| parse_fortran_float(token))
            .collect::<Vec<_>>();
        if parsed.iter().all(Option::is_some) && parsed.len() >= 2 {
            let count = parsed.len();
            if let Some(expected) = column_count {
                if expected != count {
                    continue;
                }
            } else {
                column_count = Some(count);
            }
            rows.push(parsed);
            if rows.len() > MAX_PARSED_TABLE_ROWS {
                return Some(skipped_table(
                    &artifact.file_name,
                    artifact.size_bytes,
                    format!("File has more than {} numeric rows.", MAX_PARSED_TABLE_ROWS),
                ));
            }
        }
    }

    let column_count = column_count?;
    if rows.is_empty() {
        return None;
    }
    Some(EpwParsedTable {
        file_name: artifact.file_name.clone(),
        family: infer_epw_table_family(&artifact.file_name),
        title: artifact.file_name.clone(),
        column_labels: labels_from_header(header.as_deref(), column_count),
        rows,
        skipped: false,
        skip_reason: None,
    })
}

fn parse_epw_numeric_tables(
    work_path: &Path,
    artifacts: &[WannierArtifact],
) -> Vec<EpwParsedTable> {
    artifacts
        .iter()
        .filter_map(|artifact| parse_numeric_table_file(work_path, artifact))
        .collect()
}

fn parse_epw_spectral_data(
    output: &str,
    parsed_tables: &[EpwParsedTable],
) -> Option<EpwSpectralData> {
    let self_energy_modes = parse_epw_self_energy_modes(output);
    let tables = parsed_tables
        .iter()
        .filter(|table| {
            matches!(
                table.family.as_str(),
                "spectral_function" | "self_energy" | "frequency" | "scattering"
            )
        })
        .cloned()
        .collect::<Vec<_>>();
    let mut notes = Vec::new();
    let mut warnings = Vec::new();
    for line in output.lines() {
        let lower = line.to_ascii_lowercase();
        if lower.contains("self-energy")
            || lower.contains("spectral function")
            || lower.contains("lambda___")
        {
            push_unique_string(&mut notes, line.trim());
        }
        if lower.contains("warning") && (lower.contains("self") || lower.contains("spectral")) {
            push_unique_string(&mut warnings, line.trim());
        }
    }
    if self_energy_modes.is_empty() && tables.is_empty() && notes.is_empty() {
        return None;
    }
    Some(EpwSpectralData {
        self_energy_modes,
        tables,
        notes,
        warnings,
    })
}

fn epw_output_reports_completion(output: &str) -> bool {
    if output.contains("JOB DONE") {
        return true;
    }
    let lower = output.to_ascii_lowercase();
    lower.contains("total program execution") && lower.contains("epw")
}

pub fn parse_epw_result_summary(
    output: &str,
    artifacts: Vec<WannierArtifact>,
    runner_completed: bool,
) -> EpwResultSummaryV1 {
    let mut summary = EpwResultSummaryV1 {
        completed: runner_completed || epw_output_reports_completion(output),
        elapsed_seconds: None,
        core_metrics: HashMap::new(),
        generated_outputs: artifacts,
        unknown_metrics: HashMap::new(),
        parse_partial: false,
        notes: Vec::new(),
        warnings: Vec::new(),
        parse_coverage: Vec::new(),
    };

    let elapsed_re =
        Regex::new(r"finished with exit=\d+\s+elapsed=(\d+)s").expect("valid elapsed regex");
    if let Some(captures) = elapsed_re.captures_iter(output).last() {
        if let Some(seconds) = captures
            .get(1)
            .and_then(|value| value.as_str().parse::<f64>().ok())
        {
            summary.elapsed_seconds = Some(seconds);
            summary
                .core_metrics
                .insert("elapsed_seconds".to_string(), seconds);
        }
    }

    let lambda_re =
        Regex::new(r"(?i)\blambda\s*[:=]\s*([-\d.Ee+Dd+]+)").expect("valid lambda regex");
    if let Some(captures) = lambda_re.captures_iter(output).last() {
        if let Some(value) = captures
            .get(1)
            .and_then(|entry| parse_fortran_float(entry.as_str()))
        {
            summary.core_metrics.insert("lambda".to_string(), value);
        }
    }
    let lambda_tot_re = Regex::new(r"(?i)lambda___\(\s*tot\s*\)\s*=\s*([-\d.Ee+Dd+]+)")
        .expect("valid lambda tot regex");
    if !summary.core_metrics.contains_key("lambda") {
        if let Some(captures) = lambda_tot_re.captures_iter(output).last() {
            if let Some(value) = captures
                .get(1)
                .and_then(|entry| parse_fortran_float(entry.as_str()))
            {
                summary.core_metrics.insert("lambda".to_string(), value);
            }
        }
    }
    let lambda_tr_re =
        Regex::new(r"(?i)\blambda_tr\s*[:=]\s*([-\d.Ee+Dd+]+)").expect("valid lambda_tr regex");
    if let Some(captures) = lambda_tr_re.captures_iter(output).last() {
        if let Some(value) = captures
            .get(1)
            .and_then(|entry| parse_fortran_float(entry.as_str()))
        {
            summary.core_metrics.insert("lambda_tr".to_string(), value);
        }
    }

    let tc_re = Regex::new(r"(?i)\btc\s*=\s*([-\d.Ee+Dd+]+)").expect("valid tc regex");
    if let Some(captures) = tc_re.captures_iter(output).last() {
        if let Some(value) = captures
            .get(1)
            .and_then(|entry| parse_fortran_float(entry.as_str()))
        {
            summary.core_metrics.insert("tc".to_string(), value);
        }
    }
    let tc_named_re = Regex::new(
        r"(?i)Estimated Tc using (?:Allen-Dynes|McMillan|SISSO)[^=]*=\s*([-\d.Ee+Dd+]+)\s*K",
    )
    .expect("valid named Tc regex");
    if !summary.core_metrics.contains_key("tc") {
        if let Some(captures) = tc_named_re.captures_iter(output).last() {
            if let Some(value) = captures
                .get(1)
                .and_then(|entry| parse_fortran_float(entry.as_str()))
            {
                summary.core_metrics.insert("tc".to_string(), value);
            }
        }
    }

    let warning_re = Regex::new(r"(?i)\b(warning|error|fatal)\b").expect("valid EPW warning regex");
    for line in output.lines() {
        if warning_re.is_match(line) {
            push_unique_string(&mut summary.warnings, line.trim());
        }
    }

    let temps_default_re =
        Regex::new(r"(?i)No temperature supplied\.\s*Setting temps\(:\)\s*to\s*([-\d.Ee+]+)\s*K")
            .expect("valid EPW default temperature regex");
    if let Some(captures) = temps_default_re.captures_iter(output).last() {
        if let Some(raw) = captures.get(1).map(|entry| entry.as_str()) {
            if let Ok(value) = raw.parse::<f64>() {
                summary.notes.push(format!(
                    "No EPW temperatures were supplied; EPW defaulted `temps(:)` to {} K.",
                    format_decimal(value)
                ));
                summary.unknown_metrics.insert(
                    "default_temperature_k".to_string(),
                    serde_json::Value::from(value),
                );
            } else {
                summary.notes.push(
                    "No EPW temperatures were supplied; EPW used an internal default for `temps(:)`."
                        .to_string(),
                );
            }
        }
    }

    if summary.completed && summary.core_metrics.is_empty() {
        summary.parse_partial = true;
        summary.notes.push(
            "Run completed but no known EPW scalar metrics were extracted from stdout.".to_string(),
        );
    }
    summary
}

pub fn parse_epw_result_v2(
    output: &str,
    work_path: &Path,
    artifacts: Vec<WannierArtifact>,
    runner_completed: bool,
) -> EpwParsedResultV2 {
    let parsed_tables = parse_epw_numeric_tables(work_path, &artifacts);
    let transport = parse_epw_transport_data(output);
    let superconductivity = parse_epw_superconductivity_data(output, &parsed_tables);
    let spectral = parse_epw_spectral_data(output, &parsed_tables);
    let mut summary = parse_epw_result_summary(output, artifacts, runner_completed);

    if let Some(transport) = &transport {
        if !transport.mobility.is_empty() {
            summary.core_metrics.insert(
                "epw_mobility_datasets".to_string(),
                transport.mobility.len() as f64,
            );
            summary.parse_coverage.push("transport".to_string());
        }
        if transport.scattering_file_notices > 0 {
            summary.core_metrics.insert(
                "scattering_file_notices".to_string(),
                transport.scattering_file_notices as f64,
            );
        }
        for warning in &transport.warnings {
            push_unique_string(&mut summary.warnings, warning);
        }
    }
    if let Some(superconductivity) = &superconductivity {
        if let Some(value) = superconductivity
            .lambda
            .or(superconductivity.electron_phonon_coupling)
        {
            summary.core_metrics.insert("lambda".to_string(), value);
        }
        if let Some(value) = superconductivity.lambda_tr {
            summary.core_metrics.insert("lambda_tr".to_string(), value);
        }
        if let Some(value) = superconductivity
            .tc_allen_dynes_k
            .or(superconductivity.tc_mcmillan_k)
            .or(superconductivity.tc_sisso_k)
        {
            summary.core_metrics.insert("tc".to_string(), value);
        }
        summary.parse_coverage.push("superconductivity".to_string());
        for warning in &superconductivity.warnings {
            push_unique_string(&mut summary.warnings, warning);
        }
    }
    if let Some(spectral) = &spectral {
        if !spectral.self_energy_modes.is_empty() {
            summary.core_metrics.insert(
                "self_energy_modes".to_string(),
                spectral.self_energy_modes.len() as f64,
            );
        }
        summary.parse_coverage.push("spectral".to_string());
        for warning in &spectral.warnings {
            push_unique_string(&mut summary.warnings, warning);
        }
    }
    let parsed_preview_tables = parsed_tables.iter().filter(|table| !table.skipped).count();
    if parsed_preview_tables > 0 {
        summary
            .core_metrics
            .insert("parsed_tables".to_string(), parsed_preview_tables as f64);
        summary.parse_coverage.push("numeric_tables".to_string());
    }
    for table in &parsed_tables {
        if table.skipped && !is_epw_scratch_artifact(&table.file_name) {
            summary.notes.push(format!(
                "Skipped preview for {}: {}",
                table.file_name,
                table.skip_reason.as_deref().unwrap_or("not previewed")
            ));
        }
    }

    if summary.completed
        && summary.parse_coverage.is_empty()
        && summary
            .core_metrics
            .keys()
            .all(|key| key == "elapsed_seconds")
    {
        summary.parse_partial = true;
        summary
            .notes
            .push("Run completed but no structured EPW result families were detected.".to_string());
    }

    EpwParsedResultV2 {
        summary,
        transport,
        superconductivity,
        spectral,
        parsed_tables,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_dir(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("qcortado_epw_{}_{}", name, uuid::Uuid::new_v4()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn parses_serta_mobility_blocks() {
        let dir = temp_dir("serta");
        let output = r#"
     ===================================================================
     Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]
     ===================================================================

      300.000      6.8411       0.100000E+14       0.148369E+03  x-axis
                                                   0.196063E+03  y-axis
                                                   0.207783E+03  z-axis
                                                   0.184071E+03     avg
      350.000      6.9357       0.100000E+14       0.107924E+03  x-axis
                                                   0.141331E+03  y-axis
                                                   0.149746E+03  z-axis
                                                   0.133001E+03     avg

     Note: Mobility are sorted by ascending values and might not correspond
                                              to the expected (x,y,z) axis.
     Writing scattering rate to file
     JOB DONE.
"#;
        let parsed = parse_epw_result_v2(output, &dir, Vec::new(), true);
        let transport = parsed.transport.expect("transport parsed");
        assert_eq!(transport.scattering_file_notices, 1);
        assert_eq!(transport.mobility.len(), 1);
        let dataset = &transport.mobility[0];
        assert_eq!(dataset.carrier_type, "hole");
        assert_eq!(dataset.method, "generic");
        assert_eq!(dataset.temperature_values_k, vec![300.0, 350.0]);
        assert_eq!(dataset.component_labels, vec!["xx", "yy", "zz", "avg"]);
        assert_eq!(dataset.mobility_values[0][0], Some(148.369));
        assert_eq!(dataset.mobility_values[3][1], Some(133.001));
        assert!(transport
            .warnings
            .iter()
            .any(|entry| entry.contains("Mobility are sorted")));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn parses_ibte_mobility_iterations() {
        let dir = temp_dir("ibte");
        let output = r#"
     BTE in the self-energy relaxation time approximation (SERTA)

       Temp     Fermi   Elec density  Population SR            Drift Elec mobility
        [K]      [eV]     [cm^-3]      [e per cell]                    [cm^2/Vs]

      400.000   7.0305   0.10000E+14  -0.35390E-11    0.290381E+02    0.661168E+01    0.661168E+01
                                       0.86631E-11    0.661168E+01    0.290382E+02    0.661169E+01
                                      -0.45474E-10    0.661168E+01    0.661169E+01    0.290381E+02
     Start solving iterative Boltzmann Transport Equation
     Iteration number:         2

       Temp     Fermi   Elec density  Population SR            Drift Elec mobility
        [K]      [eV]     [cm^-3]      [e per cell]                    [cm^2/Vs]

      400.000   7.0305   0.10000E+14  -0.60592E-11    0.252006E+02    0.520832E+01    0.520832E+01
                                       0.36870E-10    0.520832E+01    0.252006E+02    0.520833E+01
                                      -0.20633E-09    0.520832E+01    0.520833E+01    0.252005E+02

                                                      0.767385E+00    Max error
     The iteration reached the maximum but did not converge.
"#;
        let parsed = parse_epw_result_v2(output, &dir, Vec::new(), true);
        let datasets = parsed.transport.expect("transport parsed").mobility;
        assert_eq!(datasets.len(), 2);
        assert_eq!(datasets[0].method, "SERTA");
        assert_eq!(datasets[1].method, "IBTE");
        assert_eq!(datasets[1].iteration, Some(2));
        assert_eq!(datasets[1].converged, Some(false));
        assert_eq!(datasets[1].max_error, Some(0.767385));
        assert_eq!(
            datasets[1].component_labels,
            vec!["xx", "yy", "zz", "xy", "xz", "yx", "yz", "zx", "zy"]
        );
        assert_eq!(datasets[1].mobility_values[0][0], Some(25.2006));
        assert_eq!(datasets[1].mobility_values[2][0], Some(25.2005));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn parses_superconductivity_and_eliashberg_summary() {
        let dir = temp_dir("super");
        let output = r#"
     Electron-phonon coupling strength =    0.8714948
     Estimated Tc using McMillan expression =  26.4082 K for muc =   0.1600
     Estimated Tc using Allen-Dynes modified McMillan expression =  27.4569 K
     Estimated Tc using SISSO machine learning model =  31.5015 K
     Estimated w_log =  61.4704 meV
     Estimated BCS superconducting gap using McMillan Tc =   4.0052 meV
     temp(  1) =     15.00000 K
        iter      ethr          znormi      deltai [meV]
          1   2.720968E+00   1.843051E+00   4.851085E+00
          2   1.228282E-01   1.839523E+00   5.515340E+00
     Convergence was reached in nsiter =      2
     Temp (itemp =   1) =   15.000 K  Free energy =    -0.048960 meV
     Min. / Max. values of superconducting gap =     2.639256   14.001996 meV
     JOB DONE.
"#;
        let parsed = parse_epw_result_v2(output, &dir, Vec::new(), true);
        let supercon = parsed.superconductivity.expect("superconductivity parsed");
        assert_eq!(supercon.electron_phonon_coupling, Some(0.8714948));
        assert_eq!(supercon.tc_mcmillan_k, Some(26.4082));
        assert_eq!(supercon.tc_allen_dynes_k, Some(27.4569));
        assert_eq!(supercon.tc_sisso_k, Some(31.5015));
        assert_eq!(supercon.w_log_mev, Some(61.4704));
        assert_eq!(supercon.bcs_gap_mev, Some(4.0052));
        assert_eq!(supercon.muc, Some(0.16));
        assert_eq!(supercon.eliashberg_converged, Some(true));
        assert_eq!(supercon.eliashberg_nsiter, Some(2));
        assert_eq!(supercon.eliashberg_iterations.len(), 2);
        assert_eq!(supercon.gap_summaries[0].gap_max_mev, Some(14.001996));
        assert!(parsed
            .summary
            .parse_coverage
            .contains(&"superconductivity".to_string()));
        let json = serde_json::to_string(&EpwCalculationV1 {
            schema_version: EPW_SCHEMA_VERSION,
            sources: EpwSourcesV1::default(),
            input: EpwInputConfig::default(),
            goals: EpwGoalSelection::default(),
            runtime: EpwRuntimeConfig::default(),
            extensions: EpwExtensionsV1::default(),
            artifacts: Vec::new(),
            result_summary: parsed.summary,
            transport: parsed.transport,
            superconductivity: Some(supercon),
            spectral: parsed.spectral,
            parsed_tables: parsed.parsed_tables,
            errors: Vec::new(),
        })
        .unwrap();
        assert!(json.contains("\"schema_version\":2"));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn parses_numeric_artifacts_and_skips_large_tables() {
        let dir = temp_dir("tables");
        std::fs::write(
            dir.join("sample.a2f"),
            "# omega a2f lambda\n0.0 0.0 0.0\n1.0 2.5 3.5\nbad row\n",
        )
        .unwrap();
        let artifacts = vec![
            WannierArtifact {
                file_name: "sample.a2f".to_string(),
                size_bytes: std::fs::metadata(dir.join("sample.a2f")).unwrap().len(),
            },
            WannierArtifact {
                file_name: "huge.a2f".to_string(),
                size_bytes: MAX_PARSED_TABLE_BYTES + 1,
            },
        ];
        let parsed = parse_epw_result_v2("JOB DONE.", &dir, artifacts, true);
        assert_eq!(parsed.parsed_tables.len(), 2);
        let a2f = parsed
            .parsed_tables
            .iter()
            .find(|table| table.file_name == "sample.a2f")
            .unwrap();
        assert_eq!(a2f.family, "a2f");
        assert_eq!(a2f.column_labels, vec!["omega", "a2f", "lambda"]);
        assert_eq!(a2f.rows.len(), 2);
        let huge = parsed
            .parsed_tables
            .iter()
            .find(|table| table.file_name == "huge.a2f")
            .unwrap();
        assert!(huge.skipped);
        assert!(huge.skip_reason.as_deref().unwrap_or("").contains("larger"));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn treats_successful_epw_without_job_done_as_completed() {
        let dir = temp_dir("completed_no_job_done");
        let output = r#"
     Program EPW v.6.0 starts
     lambda = 0.42
     Total program execution
"#;
        let parsed = parse_epw_result_v2(output, &dir, Vec::new(), true);
        assert!(parsed.summary.completed);
        assert!(!parsed.summary.parse_partial);

        let failed = parse_epw_result_v2("lambda = 0.42", &dir, Vec::new(), false);
        assert!(!failed.summary.completed);
        assert!(!failed.summary.parse_partial);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn excludes_qe_scratch_files_from_parsed_tables() {
        let dir = temp_dir("scratch");
        std::fs::create_dir_all(dir.join("tmp/qcortado_scf.save")).unwrap();
        std::fs::create_dir_all(dir.join("save/qcortado_scf.phsave")).unwrap();
        std::fs::write(dir.join("tmp/qcortado_scf.save/wfc1.dat"), "1 2\n3 4\n").unwrap();
        std::fs::write(
            dir.join("tmp/qcortado_scf.save/charge-density.dat"),
            "1 2\n",
        )
        .unwrap();
        std::fs::write(dir.join("save/qcortado_scf.dvscf_q1"), "1 2\n").unwrap();
        std::fs::write(dir.join("linewidth.phself"), "1.0 2.0\n").unwrap();
        let artifacts = vec![
            WannierArtifact {
                file_name: "tmp/qcortado_scf.save/wfc1.dat".to_string(),
                size_bytes: 8,
            },
            WannierArtifact {
                file_name: "tmp/qcortado_scf.save/charge-density.dat".to_string(),
                size_bytes: 4,
            },
            WannierArtifact {
                file_name: "save/qcortado_scf.dvscf_q1".to_string(),
                size_bytes: 4,
            },
            WannierArtifact {
                file_name: "linewidth.phself".to_string(),
                size_bytes: 8,
            },
        ];
        let parsed = parse_epw_result_v2("JOB DONE.", &dir, artifacts, true);
        assert!(parsed.parsed_tables.iter().all(
            |table| !table.file_name.contains(".save") && !table.file_name.starts_with("save/")
        ));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn epw_goals_emit_expected_keywords_and_meshes() {
        let mut config = EpwCalculationConfig {
            project_id: "project".to_string(),
            source_phonon_calc_id: "phonon".to_string(),
            source_wannier_calc_id: None,
            source_scf_calc_id: None,
            mode: None,
            input: EpwInputConfig {
                k_mesh: [6, 6, 6],
                q_mesh: [3, 3, 3],
                fine_k_mesh: Some([24, 24, 24]),
                fine_q_mesh: Some([12, 12, 12]),
                ..EpwInputConfig::default()
            },
            runtime: EpwRuntimeConfig::default(),
            extensions: EpwExtensionsV1::default(),
            goals: EpwGoalSelection {
                electron_self_energy: true,
                transport_mobility: true,
                superconductivity: true,
                ..EpwGoalSelection::default()
            },
            advanced_overrides: HashMap::new(),
        };
        config.input.coarse_k_mesh = Some(config.input.k_mesh);
        config.input.coarse_q_mesh = Some(config.input.q_mesh);

        let keywords = build_epw_keyword_map(&config).unwrap();
        assert_eq!(keywords.get("nk1").map(String::as_str), Some("6"));
        assert_eq!(keywords.get("nkf1").map(String::as_str), Some("24"));
        assert_eq!(keywords.get("nq1").map(String::as_str), Some("3"));
        assert_eq!(keywords.get("nqf1").map(String::as_str), Some("12"));
        assert_eq!(keywords.get("elph").map(String::as_str), Some(".true."));
        assert_eq!(
            keywords.get("phonselfen").map(String::as_str),
            Some(".true.")
        );
        assert_eq!(keywords.get("a2f").map(String::as_str), Some(".true."));
        assert_eq!(
            keywords.get("elecselfen").map(String::as_str),
            Some(".true.")
        );
        assert_eq!(
            keywords.get("specfun_el").map(String::as_str),
            Some(".true.")
        );
        assert_eq!(
            keywords.get("scattering").map(String::as_str),
            Some(".true.")
        );
        assert_eq!(
            keywords.get("scattering_serta").map(String::as_str),
            Some(".true.")
        );
        assert_eq!(keywords.get("int_mob").map(String::as_str), Some(".true."));
        assert_eq!(keywords.get("ephwrite").map(String::as_str), Some(".true."));
        assert_eq!(
            keywords.get("eliashberg").map(String::as_str),
            Some(".true.")
        );
        assert_eq!(keywords.get("liso").map(String::as_str), Some(".true."));
        assert_eq!(keywords.get("limag").map(String::as_str), Some(".true."));
    }
}
