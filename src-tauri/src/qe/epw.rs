//! EPW (`epw.x`) helpers: config validation, input generation, and lightweight output parsing.

use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write;
use std::path::Path;

use super::wannier::WannierArtifact;

pub const EPW_SCHEMA_VERSION: u32 = 1;

fn default_prefix() -> String {
    "qcortado_scf".to_string()
}

fn default_outdir() -> String {
    "./tmp".to_string()
}

fn default_phonon_dir() -> String {
    "./phonon".to_string()
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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpwCalculationV1 {
    pub schema_version: u32,
    pub sources: EpwSourcesV1,
    pub input: EpwInputConfig,
    pub runtime: EpwRuntimeConfig,
    #[serde(default)]
    pub extensions: EpwExtensionsV1,
    #[serde(default)]
    pub artifacts: Vec<WannierArtifact>,
    pub result_summary: EpwResultSummaryV1,
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
    if config.input.k_mesh.iter().any(|entry| *entry == 0) {
        return Err("EPW k_mesh values must all be positive.".to_string());
    }
    if config.input.q_mesh.iter().any(|entry| *entry == 0) {
        return Err("EPW q_mesh values must all be positive.".to_string());
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
    keywords.insert("prefix".to_string(), quote_string(config.input.prefix.trim()));
    keywords.insert("outdir".to_string(), quote_string(config.input.outdir.trim()));
    keywords.insert(
        "dvscf_dir".to_string(),
        quote_string(config.input.dvscf_dir.trim()),
    );
    keywords.insert(
        "nk1".to_string(),
        config.input.k_mesh[0].to_string(),
    );
    keywords.insert(
        "nk2".to_string(),
        config.input.k_mesh[1].to_string(),
    );
    keywords.insert(
        "nk3".to_string(),
        config.input.k_mesh[2].to_string(),
    );
    keywords.insert(
        "nq1".to_string(),
        config.input.q_mesh[0].to_string(),
    );
    keywords.insert(
        "nq2".to_string(),
        config.input.q_mesh[1].to_string(),
    );
    keywords.insert(
        "nq3".to_string(),
        config.input.q_mesh[2].to_string(),
    );
    keywords.insert(
        "epbwrite".to_string(),
        if config.input.epbwrite {
            ".true.".to_string()
        } else {
            ".false.".to_string()
        },
    );
    keywords.insert(
        "epbread".to_string(),
        if config.input.epbread {
            ".true.".to_string()
        } else {
            ".false.".to_string()
        },
    );
    keywords.insert(
        "epwwrite".to_string(),
        if config.input.epwwrite {
            ".true.".to_string()
        } else {
            ".false.".to_string()
        },
    );
    keywords.insert(
        "epwread".to_string(),
        if config.input.epwread {
            ".true.".to_string()
        } else {
            ".false.".to_string()
        },
    );
    keywords.insert(
        "wannierize".to_string(),
        if config.input.wannierize {
            ".true.".to_string()
        } else {
            ".false.".to_string()
        },
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

pub fn build_epw_input_preview(config: &EpwCalculationConfig) -> Result<EpwInputPreviewResult, String> {
    let merged_keywords = build_epw_keyword_map(config)?;
    let input_text = render_epw_input(&merged_keywords);
    Ok(EpwInputPreviewResult {
        schema_version: EPW_SCHEMA_VERSION,
        input_text,
        merged_keywords,
    })
}

pub fn collect_epw_artifacts(work_path: &Path) -> Vec<WannierArtifact> {
    let mut artifacts: Vec<WannierArtifact> = Vec::new();
    let Ok(entries) = std::fs::read_dir(work_path) else {
        return artifacts;
    };

    for entry in entries.flatten() {
        let Ok(meta) = entry.metadata() else {
            continue;
        };
        if !meta.is_file() {
            continue;
        }
        let file_name = entry.file_name().to_string_lossy().to_string();
        let lower = file_name.to_ascii_lowercase();
        let keep = lower == "epw.in"
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
            || lower.ends_with(".xml")
            || lower.ends_with(".dat");
        if keep {
            artifacts.push(WannierArtifact {
                file_name,
                size_bytes: meta.len(),
            });
        }
    }

    artifacts.sort_by(|a, b| a.file_name.cmp(&b.file_name));
    artifacts
}

pub fn parse_epw_result_summary(
    output: &str,
    artifacts: Vec<WannierArtifact>,
) -> EpwResultSummaryV1 {
    let mut summary = EpwResultSummaryV1 {
        completed: output.contains("JOB DONE"),
        elapsed_seconds: None,
        core_metrics: HashMap::new(),
        generated_outputs: artifacts,
        unknown_metrics: HashMap::new(),
        parse_partial: false,
        notes: Vec::new(),
    };

    let elapsed_re =
        Regex::new(r"finished with exit=\d+\s+elapsed=(\d+)s").expect("valid elapsed regex");
    if let Some(captures) = elapsed_re.captures_iter(output).last() {
        if let Some(seconds) = captures.get(1).and_then(|value| value.as_str().parse::<f64>().ok()) {
            summary.elapsed_seconds = Some(seconds);
            summary
                .core_metrics
                .insert("elapsed_seconds".to_string(), seconds);
        }
    }

    let lambda_re = Regex::new(r"(?i)\blambda\s*=\s*([-\d.Ee+]+)").expect("valid lambda regex");
    if let Some(captures) = lambda_re.captures_iter(output).last() {
        if let Some(value) = captures.get(1).and_then(|entry| entry.as_str().parse::<f64>().ok()) {
            summary.core_metrics.insert("lambda".to_string(), value);
        }
    }

    let tc_re = Regex::new(r"(?i)\btc\s*=\s*([-\d.Ee+]+)").expect("valid tc regex");
    if let Some(captures) = tc_re.captures_iter(output).last() {
        if let Some(value) = captures.get(1).and_then(|entry| entry.as_str().parse::<f64>().ok()) {
            summary.core_metrics.insert("tc".to_string(), value);
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
                summary
                    .unknown_metrics
                    .insert("default_temperature_k".to_string(), serde_json::Value::from(value));
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
    if !summary.completed {
        summary.parse_partial = true;
        summary
            .notes
            .push("Run did not report JOB DONE; output metrics may be incomplete.".to_string());
    }

    summary
}
