//! Native WIEN2k SCF session and result helpers.

use std::collections::BTreeMap;

use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::engines::types::{
    CalculationKind, DiagnosticHint, DiagnosticSeverity, DiagnosticSource, EnergyQuantity,
    EnergyUnit, EngineId, EngineResultProvenance, EngineTaskKind, NormalizedResultSchema,
    NormalizedScfSummary, ScfConvergenceState,
};

use super::{Wien2kInitializationSettings, Wien2kScfRunSettings, Wien2kSpinMode};

const RY_TO_EV: f64 = 13.605_693_122_994;
const SUPPORTED_INITIALIZATION_VXC: [u16; 4] = [5, 11, 13, 19];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kScfSessionPhase {
    Staged,
    Initialized,
    ScfComplete,
    Failed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kScfSession {
    pub session_id: String,
    pub project_id: String,
    pub cif_id: String,
    pub source_structure_calculation_id: String,
    pub case_name: String,
    pub remote_case_dir: String,
    pub remote_install_root: String,
    pub hpc_profile_id: String,
    pub phase: Wien2kScfSessionPhase,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initialization: Option<Wien2kInitializationSettings>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latest_run: Option<Wien2kScfRunSettings>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latest_calculation_id: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub artifacts: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transcript: Vec<String>,
    pub started_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kInitializationResult {
    pub session_id: String,
    pub phase: Wien2kScfSessionPhase,
    pub native_output: String,
    pub diagnostics: Vec<String>,
    pub artifacts: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kScfExecutionResult {
    pub session_id: String,
    pub phase: Wien2kScfSessionPhase,
    pub native_output: String,
    pub diagnostics: Vec<String>,
    pub summary: NormalizedScfSummary,
    pub calculation_id: String,
    pub continuation: bool,
}

pub fn validate_initialization_settings(
    settings: &Wien2kInitializationSettings,
) -> Result<(), String> {
    if !settings.rkmax.is_finite() || settings.rkmax <= 0.0 {
        return Err("RKMAX must be positive.".to_string());
    }
    if !settings.gmax.is_finite() || settings.gmax <= 0.0 {
        return Err("GMAX must be positive.".to_string());
    }
    if settings.lmax == 0 {
        return Err("LMAX must be positive.".to_string());
    }
    if settings.k_mesh.iter().any(|point| *point == 0) {
        return Err("All k-mesh dimensions must be positive.".to_string());
    }
    if !SUPPORTED_INITIALIZATION_VXC.contains(&settings.exchange_correlation) {
        return Err(
            "XC functional must be one of WIEN2k's native initialization options: LDA, PBE, WC, or PBEsol."
                .to_string(),
        );
    }
    if !settings.lstart_energy_cutoff_ry.is_finite() {
        return Err("The LSTART energy cutoff must be finite.".to_string());
    }
    Ok(())
}

pub fn validate_run_settings(settings: &Wien2kScfRunSettings) -> Result<(), String> {
    if settings.max_iterations == 0 {
        return Err("Maximum SCF iterations must be positive.".to_string());
    }
    if !settings.energy_convergence_ry.is_finite() || settings.energy_convergence_ry <= 0.0 {
        return Err("Energy convergence must be positive.".to_string());
    }
    if !settings.charge_convergence.is_finite() || settings.charge_convergence <= 0.0 {
        return Err("Charge convergence must be positive.".to_string());
    }
    if settings
        .force_convergence_mry_bohr
        .is_some_and(|value| !value.is_finite() || value <= 0.0)
    {
        return Err("Force convergence must be positive when set.".to_string());
    }
    Ok(())
}

pub fn required_initialization_suffixes(spin_mode: Wien2kSpinMode) -> Vec<&'static str> {
    let mut required = vec!["struct", "in0", "inc", "inm", "klist", "clmsum"];
    if spin_mode == Wien2kSpinMode::SpinPolarized {
        required.extend(["clmup", "clmdn"]);
    }
    required
}

pub fn initialization_diagnostics(output: &str) -> Vec<String> {
    let upper = output.to_ascii_uppercase();
    let mut diagnostics = Vec::new();
    if initialization_has_error_marker(output) {
        diagnostics.push(
            "WIEN2k initialization reported an error; review the native output before retrying."
                .to_string(),
        );
    }
    if upper.contains("OVERLAP") {
        diagnostics.push("WIEN2k reported sphere overlap during initialization; return to the Structure workflow to change RMT values.".to_string());
    }
    diagnostics
}

fn initialization_has_error_marker(output: &str) -> bool {
    output.lines().any(|line| {
        let upper = line.trim().to_ascii_uppercase();
        if upper.is_empty() {
            return false;
        }
        if upper.contains("STOPPED") {
            return true;
        }
        upper
            .split(|character: char| !character.is_ascii_alphanumeric())
            .any(|token| token == "ERROR")
    })
}

pub fn parse_scf_summary(
    project_id: &str,
    cif_id: &str,
    source_structure_calculation_id: &str,
    settings: &Wien2kScfRunSettings,
    scf_output: &str,
    dayfile: &str,
    command_failed: bool,
) -> NormalizedScfSummary {
    let energies = tagged_values(scf_output, ":ENE");
    let total_energy = energies.last().copied();
    let charge_residual = tagged_values(scf_output, ":DIS").last().copied();
    let energy_delta = if energies.len() >= 2 {
        Some((energies[energies.len() - 1] - energies[energies.len() - 2]).abs())
    } else {
        None
    };
    let fermi_energy_ev = tagged_values(scf_output, ":FER")
        .last()
        .copied()
        .map(|value| value * RY_TO_EV);
    let total_magnetization = tagged_values(scf_output, ":MMTOT").last().copied();
    let combined_upper = format!("{}\n{}", scf_output, dayfile).to_ascii_uppercase();
    let native_converged = combined_upper.contains("CONVERGED")
        || combined_upper.contains("CONVERGENCE CRITERIA SATISFIED");
    let threshold_converged = charge_residual
        .is_some_and(|value| value <= settings.charge_convergence)
        && energy_delta.is_some_and(|value| value <= settings.energy_convergence_ry);
    let convergence = if command_failed || combined_upper.contains("ERROR") {
        ScfConvergenceState::Failed
    } else if native_converged || threshold_converged {
        ScfConvergenceState::Converged
    } else {
        ScfConvergenceState::NotConverged
    };

    let mut diagnostics = Vec::new();
    match convergence {
        ScfConvergenceState::Failed => diagnostics.push(diagnostic(
            DiagnosticSeverity::Error,
            "WIEN2k SCF failed.",
            "Review case.dayfile and case.scf for the native failure.",
        )),
        ScfConvergenceState::NotConverged => diagnostics.push(diagnostic(
            DiagnosticSeverity::Warning,
            "WIEN2k SCF did not converge.",
            "Continue the retained case with -NI or adjust the convergence settings.",
        )),
        _ => {}
    }

    let mut metadata = serde_json::Map::new();
    if let Some(value) = charge_residual {
        metadata.insert("chargeResidual".to_string(), serde_json::json!(value));
    }
    if let Some(value) = energy_delta {
        metadata.insert("energyDeltaRy".to_string(), serde_json::json!(value));
    }
    metadata.insert(
        "spinMode".to_string(),
        serde_json::to_value(settings.spin_mode).unwrap_or(serde_json::Value::Null),
    );

    NormalizedScfSummary {
        schema: NormalizedResultSchema::ScfSummaryV1,
        provenance: EngineResultProvenance {
            engine_id: EngineId::Wien2k,
            task_kind: EngineTaskKind::Scf,
            calculation_kind: Some(CalculationKind::Scf),
            calculation_id: None,
            project_id: Some(project_id.to_string()),
            cif_id: Some(cif_id.to_string()),
            source_calculation_ids: vec![source_structure_calculation_id.to_string()],
            generated_at: None,
        },
        convergence,
        total_energy: total_energy.map(|value| EnergyQuantity {
            value,
            unit: EnergyUnit::Ry,
        }),
        fermi_energy_ev,
        scf_steps: (!energies.is_empty()).then_some(energies.len() as u32),
        wall_time_seconds: None,
        total_magnetization,
        atomic_magnetic_moments: None,
        forces: None,
        stress: None,
        diagnostics,
        artifacts: Vec::new(),
        metadata: Some(metadata),
    }
}

fn tagged_values(content: &str, tag: &str) -> Vec<f64> {
    content
        .lines()
        .filter(|line| line.contains(tag))
        .filter_map(last_numeric_value)
        .collect()
}

fn last_numeric_value(line: &str) -> Option<f64> {
    let matcher = Regex::new(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?").ok()?;
    matcher.find_iter(line).last().and_then(|capture| {
        capture
            .as_str()
            .replace(['D', 'd'], "E")
            .parse::<f64>()
            .ok()
    })
}

fn diagnostic(severity: DiagnosticSeverity, title: &str, message: &str) -> DiagnosticHint {
    DiagnosticHint {
        id: None,
        severity,
        source: DiagnosticSource::Engine,
        code: None,
        title: title.to_string(),
        message: message.to_string(),
        suggested_action: None,
        related_artifact_ids: Vec::new(),
        metadata: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn initialization_accepts_native_lstart_xc_choices_only() {
        for value in SUPPORTED_INITIALIZATION_VXC {
            let settings = Wien2kInitializationSettings {
                exchange_correlation: value,
                ..Wien2kInitializationSettings::default()
            };
            assert!(validate_initialization_settings(&settings).is_ok());
        }

        let invalid = Wien2kInitializationSettings {
            exchange_correlation: 999,
            ..Wien2kInitializationSettings::default()
        };
        assert!(validate_initialization_settings(&invalid).is_err());
    }

    #[test]
    fn initialization_diagnostics_ignore_benign_error_plural_warnings() {
        let output = "\
 atom 1 has a large sphere and is a heavy element, consider setting HDLOs and/or larger LVNS\n\
 Atomic spheres .gt. 2.35 may lead to linearization errors\n\
 For more accuracy rerun with -hdlo switch\n\
  init_lapw finished ok\n";

        assert!(initialization_diagnostics(output).is_empty());
    }

    #[test]
    fn initialization_diagnostics_report_explicit_error_markers() {
        let output = "[QCortado] ERROR: init_lapw did not generate case.in0\n";

        assert!(initialization_diagnostics(output)
            .iter()
            .any(|diagnostic| diagnostic.contains("reported an error")));
    }

    #[test]
    fn initialization_diagnostics_report_stopped_runs() {
        let output = "LSTART STOPPED in case.outputst\n";

        assert!(initialization_diagnostics(output)
            .iter()
            .any(|diagnostic| diagnostic.contains("reported an error")));
    }

    #[test]
    fn parser_reports_converged_scalar_output() {
        let output = "\
:ENE  : TOTAL ENERGY IN Ry = -10.000500\n\
:DIS  : CHARGE DISTANCE = 0.001\n\
:ENE  : TOTAL ENERGY IN Ry = -10.000550\n\
:DIS  : CHARGE DISTANCE = 0.000050\n\
:FER  : F E R M I - ENERGY = 0.250000\n";
        let summary = parse_scf_summary(
            "project",
            "cif",
            "struct",
            &Wien2kScfRunSettings::default(),
            output,
            "CONVERGED",
            false,
        );

        assert_eq!(summary.convergence, ScfConvergenceState::Converged);
        assert_eq!(summary.scf_steps, Some(2));
        assert_eq!(summary.total_energy.expect("energy").value, -10.00055);
        assert!(summary.fermi_energy_ev.expect("fermi") > 3.4);
    }

    #[test]
    fn parser_preserves_nonconverged_spin_attempt_for_continuation() {
        let settings = Wien2kScfRunSettings {
            spin_mode: Wien2kSpinMode::SpinPolarized,
            ..Wien2kScfRunSettings::default()
        };
        let summary = parse_scf_summary(
            "project",
            "cif",
            "struct",
            &settings,
            ":ENE : TOTAL ENERGY IN Ry = -5.0\n:DIS : CHARGE DISTANCE = 0.4\n:MMTOT: 2.0\n",
            "",
            false,
        );

        assert_eq!(summary.convergence, ScfConvergenceState::NotConverged);
        assert_eq!(summary.total_magnetization, Some(2.0));
        assert_eq!(summary.diagnostics.len(), 1);
    }
}
