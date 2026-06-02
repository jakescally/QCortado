//! Native WIEN2k SCF session and result helpers.

use std::collections::BTreeMap;

use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::engines::types::{
    CalculationKind, DiagnosticHint, DiagnosticSeverity, DiagnosticSource, EnergyQuantity,
    EnergyUnit, EngineId, EngineResultProvenance, EngineTaskKind, NormalizedResultSchema,
    NormalizedScfSummary, ScfConvergenceState,
};

use super::structure::Wien2kStructureSite;
use super::{
    Wien2kDftUDoubleCounting, Wien2kFermiMethod, Wien2kInitializationSettings, Wien2kMixerTrust,
    Wien2kScfRunSettings, Wien2kSpinMode,
};

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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_structure_sites: Vec<Wien2kStructureSite>,
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
    if settings.fermi_method != Wien2kFermiMethod::Tetra
        && settings
            .fermi_smearing_ry
            .is_none_or(|value| !value.is_finite() || value <= 0.0)
    {
        return Err("Fermi smearing must be positive when TEMP or TEMPS is selected.".to_string());
    }
    for entry in &settings.starting_magnetization {
        if entry.site_index == 0 {
            return Err("Starting magnetization site indices must be positive.".to_string());
        }
        if !entry.moment_bohr_magneton.is_finite() || entry.moment_bohr_magneton < 0.0 {
            return Err("Starting magnetization moments must be non-negative.".to_string());
        }
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
    if settings.dft_u.enabled {
        if settings.spin_mode != Wien2kSpinMode::SpinPolarized {
            return Err("WIEN2k DFT+U requires spin-polarized SCF.".to_string());
        }
        if settings.dft_u.targets.is_empty() {
            return Err("Enable at least one DFT+U target.".to_string());
        }
        for target in &settings.dft_u.targets {
            if target.site_index == 0 {
                return Err("DFT+U site indices must be positive.".to_string());
            }
            if target.orbital_l > 3 {
                return Err("DFT+U orbital l must be 0, 1, 2, or 3.".to_string());
            }
            if !is_valid_manifold(&target.manifold) {
                return Err("DFT+U manifolds must look like 3d, 4f, etc.".to_string());
            }
            if !target.u_ev.is_finite() || target.u_ev <= 0.0 {
                return Err("DFT+U U values must be positive eV values.".to_string());
            }
            if !target.j_ev.is_finite() || target.j_ev < 0.0 {
                return Err("DFT+U J values must be non-negative eV values.".to_string());
            }
        }
    }
    if !settings.mixer.greed.is_finite()
        || settings.mixer.greed <= 0.0
        || settings.mixer.greed > 1.0
    {
        return Err("Mixer greed must be in the range (0, 1].".to_string());
    }
    if settings.mixer.history == 0 {
        return Err("Mixer history must be positive.".to_string());
    }
    Ok(())
}

pub fn build_dft_u_input_files(
    settings: &Wien2kScfRunSettings,
) -> Option<BTreeMap<String, String>> {
    if !settings.dft_u.enabled || settings.dft_u.targets.is_empty() {
        return None;
    }

    let mut files = BTreeMap::new();
    files.insert("inorb".to_string(), build_case_inorb(settings));
    files.insert("indm".to_string(), build_case_indm(settings));
    Some(files)
}

pub fn build_case_inm(settings: &Wien2kScfRunSettings) -> String {
    let trust = match settings.mixer.trust {
        Wien2kMixerTrust::Default => "#",
        Wien2kMixerTrust::STIFF => "STIFF",
        Wien2kMixerTrust::STIFFER => "STIFFER",
        Wien2kMixerTrust::FAST => "FAST",
    };
    format!(
        "{} 0.d0 YES\n{:.8}\n1.0 1.0\n999 {}\n{}\n",
        serde_plain_mixer_mode(settings),
        settings.mixer.greed,
        settings.mixer.history,
        trust,
    )
}

fn build_case_inorb(settings: &Wien2kScfRunSettings) -> String {
    let targets = &settings.dft_u.targets;
    let double_counting = match settings.dft_u.double_counting {
        Wien2kDftUDoubleCounting::Amf => 0,
        Wien2kDftUDoubleCounting::Sic => 1,
        Wien2kDftUDoubleCounting::Hmf => 2,
    };
    let mut lines = vec![format!("1 {} 0", targets.len()), "PRATT,1.0".to_string()];
    for target in targets {
        lines.push(format!("{} 1 {}", target.site_index, target.orbital_l));
    }
    lines.push(double_counting.to_string());
    for target in targets {
        lines.push(format!(
            "{:.8} {:.8}",
            ev_to_ry(target.u_ev),
            ev_to_ry(target.j_ev)
        ));
    }
    lines.push(String::new());
    lines.join("\n")
}

fn build_case_indm(settings: &Wien2kScfRunSettings) -> String {
    let targets = &settings.dft_u.targets;
    let mut lines = vec!["-9.0".to_string(), targets.len().to_string()];
    for target in targets {
        lines.push(format!("{} 1 {}", target.site_index, target.orbital_l));
    }
    lines.push("0 0".to_string());
    lines.push(String::new());
    lines.join("\n")
}

fn ev_to_ry(value: f64) -> f64 {
    value / RY_TO_EV
}

fn is_valid_manifold(value: &str) -> bool {
    let mut chars = value.chars().peekable();
    let mut digits = 0;
    while chars.peek().is_some_and(|ch| ch.is_ascii_digit()) {
        digits += 1;
        chars.next();
    }
    digits > 0 && matches!(chars.next(), Some('s' | 'p' | 'd' | 'f')) && chars.next().is_none()
}

fn serde_plain_mixer_mode(settings: &Wien2kScfRunSettings) -> &'static str {
    match settings.mixer.mode {
        super::types::Wien2kMixerMode::MSR1 => "MSR1",
        super::types::Wien2kMixerMode::MSEC3 => "MSEC3",
        super::types::Wien2kMixerMode::MSEC4 => "MSEC4",
        super::types::Wien2kMixerMode::MSR2 => "MSR2",
        super::types::Wien2kMixerMode::PRATT => "PRATT",
        super::types::Wien2kMixerMode::PRAT0 => "PRAT0",
    }
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
    fn dft_u_inputs_are_rendered_in_wien2k_native_units() {
        let settings = Wien2kScfRunSettings {
            spin_mode: Wien2kSpinMode::SpinPolarized,
            dft_u: crate::engines::wien2k::Wien2kDftUSettings {
                enabled: true,
                double_counting: crate::engines::wien2k::Wien2kDftUDoubleCounting::Sic,
                targets: vec![crate::engines::wien2k::Wien2kHubbardTarget {
                    site_index: 1,
                    element: "Ni".to_string(),
                    manifold: "3d".to_string(),
                    orbital_l: 2,
                    u_ev: 6.0,
                    j_ev: 0.0,
                    recommended: true,
                    reason: None,
                }],
            },
            ..Wien2kScfRunSettings::default()
        };

        let files = build_dft_u_input_files(&settings).expect("DFT+U files");
        assert!(validate_run_settings(&settings).is_ok());
        assert!(files["inorb"].contains("1 1 0"));
        assert!(files["inorb"].contains("1 1 2"));
        assert!(files["inorb"].contains("0.440991"));
        assert!(files["indm"].contains("1 1 2"));
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
