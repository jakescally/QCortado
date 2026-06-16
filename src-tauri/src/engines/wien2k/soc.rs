//! WIEN2k spin-orbit-coupled SCF workflow types and input helpers.

use std::collections::{BTreeMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::engines::types::NormalizedScfSummary;

use super::{scf::validate_run_settings, Wien2kScfRunSettings, Wien2kSpinMode};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kSocSessionPhase {
    Staged,
    Prepared,
    SymmetryReady,
    SocComplete,
    Failed,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocRlo {
    pub atom_index: u32,
    pub energy_ry: f64,
    pub de: f64,
    #[serde(default = "default_rlo_switch")]
    pub switch: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocRloCandidate {
    pub atom_index: u32,
    pub energy_ry: f64,
    pub de: f64,
    pub switch: String,
    pub source_file: String,
}

fn default_rlo_switch() -> String {
    "STOP".to_string()
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocPrepareSettings {
    #[serde(default = "default_direction")]
    pub magnetization_direction: [i32; 3],
    #[serde(default)]
    pub disabled_atom_indices: Vec<u32>,
    #[serde(default = "default_lapw1_emax")]
    pub lapw1_emax_ry: f64,
    #[serde(default = "default_output_emin")]
    pub output_energy_min_ry: f64,
    #[serde(default = "default_output_emax")]
    pub output_energy_max_ry: f64,
    #[serde(default)]
    pub rlo_atoms: Vec<Wien2kSocRlo>,
}

const fn default_direction() -> [i32; 3] {
    [0, 0, 1]
}

const fn default_lapw1_emax() -> f64 {
    5.0
}

const fn default_output_emin() -> f64 {
    -10.0
}

const fn default_output_emax() -> f64 {
    1.9
}

impl Default for Wien2kSocPrepareSettings {
    fn default() -> Self {
        Self {
            magnetization_direction: default_direction(),
            disabled_atom_indices: Vec::new(),
            lapw1_emax_ry: default_lapw1_emax(),
            output_energy_min_ry: default_output_emin(),
            output_energy_max_ry: default_output_emax(),
            rlo_atoms: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocRunSettings {
    pub scf: Wien2kScfRunSettings,
    #[serde(default)]
    pub diagnostic_log: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocSession {
    pub session_id: String,
    pub project_id: String,
    pub cif_id: String,
    pub source_scf_calculation_id: String,
    pub source_structure_calculation_id: String,
    pub case_name: String,
    pub remote_case_dir: String,
    pub source_remote_case_dir: String,
    pub remote_install_root: String,
    pub hpc_profile_id: String,
    pub spin_mode: Wien2kSpinMode,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_k_mesh: Option<[u32; 3]>,
    #[serde(default)]
    pub source_dft_u: bool,
    pub phase: Wien2kSocSessionPhase,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latest_prepare: Option<Wien2kSocPrepareSettings>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latest_run: Option<Wien2kSocRunSettings>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latest_calculation_id: Option<String>,
    #[serde(default)]
    pub artifacts: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transcript: Vec<String>,
    pub started_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocPrepareResult {
    pub session_id: String,
    pub phase: Wien2kSocSessionPhase,
    pub symmetry_review_required: bool,
    pub native_output: String,
    pub diagnostics: Vec<String>,
    pub artifacts: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSocExecutionResult {
    pub session_id: String,
    pub phase: Wien2kSocSessionPhase,
    pub native_output: String,
    pub diagnostics: Vec<String>,
    pub summary: NormalizedScfSummary,
    pub calculation_id: String,
    pub continuation: bool,
}

pub fn validate_soc_prepare_settings(settings: &Wien2kSocPrepareSettings) -> Result<(), String> {
    if settings.magnetization_direction == [0, 0, 0] {
        return Err("SOC magnetization direction cannot be the zero vector.".to_string());
    }
    if !settings.lapw1_emax_ry.is_finite() || settings.lapw1_emax_ry <= 0.0 {
        return Err("LAPW1 EMAX must be a positive finite value.".to_string());
    }
    if !settings.output_energy_min_ry.is_finite()
        || !settings.output_energy_max_ry.is_finite()
        || settings.output_energy_min_ry >= settings.output_energy_max_ry
    {
        return Err("The case.inso output-energy window must be finite and ordered.".to_string());
    }
    let mut disabled = HashSet::new();
    for atom in &settings.disabled_atom_indices {
        if *atom == 0 || !disabled.insert(*atom) {
            return Err("Atoms without SOC must be unique positive atom indices.".to_string());
        }
    }
    let mut rlo_atoms = HashSet::new();
    for rlo in &settings.rlo_atoms {
        if rlo.atom_index == 0 || !rlo_atoms.insert(rlo.atom_index) {
            return Err("RLO atoms must be unique positive atom indices.".to_string());
        }
        if disabled.contains(&rlo.atom_index) {
            return Err("An atom cannot both use an SOC RLO and have SOC disabled.".to_string());
        }
        if !rlo.energy_ry.is_finite() || !rlo.de.is_finite() || rlo.de < 0.0 {
            return Err(
                "SOC RLO energy and de values must be finite, with non-negative de.".to_string(),
            );
        }
        if !matches!(
            rlo.switch.trim().to_ascii_uppercase().as_str(),
            "STOP" | "CONT"
        ) {
            return Err("SOC RLO switch must be STOP or CONT.".to_string());
        }
    }
    Ok(())
}

pub fn validate_soc_run_settings(
    settings: &Wien2kSocRunSettings,
    spin_mode: Wien2kSpinMode,
) -> Result<(), String> {
    validate_run_settings(&settings.scf)?;
    if settings.scf.spin_mode != spin_mode {
        return Err("SOC run spin mode must match the selected source SCF.".to_string());
    }
    if settings.scf.force_convergence_mry_bohr.is_some() || settings.scf.force_minimization {
        return Err(
            "WIEN2k SOC workflow does not support force convergence or force minimization."
                .to_string(),
        );
    }
    Ok(())
}

pub fn build_case_inso(settings: &Wien2kSocPrepareSettings) -> String {
    let mut lines = vec![
        "WFFIL".to_string(),
        " 4 0 0                        llmax,ipr,kpot".to_string(),
        format!(
            " {:.8} {:.8}           Emin, Emax",
            settings.output_energy_min_ry, settings.output_energy_max_ry
        ),
        format!(
            "    {} {} {}                       h,k,l (direction of magnetization)",
            settings.magnetization_direction[0],
            settings.magnetization_direction[1],
            settings.magnetization_direction[2]
        ),
        format!(
            " {}                             number of atoms with RLO",
            settings.rlo_atoms.len()
        ),
    ];
    for rlo in &settings.rlo_atoms {
        lines.push(format!(
            " {} {:.8} {:.8} {}       atom-number, E-param for RLO",
            rlo.atom_index,
            rlo.energy_ry,
            rlo.de,
            rlo.switch.trim().to_ascii_uppercase()
        ));
    }
    let disabled = settings
        .disabled_atom_indices
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(" ");
    lines.push(format!(
        " {} {}     number of atoms without SO, atomnumbers",
        settings.disabled_atom_indices.len(),
        disabled
    ));
    lines.push(String::new());
    lines.join("\n")
}

fn parse_wien2k_float(value: &str) -> Option<f64> {
    value.replace(['D', 'd'], "E").parse().ok()
}

pub fn parse_soc_rlo_candidates(content: &str, source_file: &str) -> Vec<Wien2kSocRloCandidate> {
    let lines = content.lines().collect::<Vec<_>>();
    let mut candidates = Vec::new();
    let mut line_index = 2;
    let mut atom_index = 1_u32;

    while line_index < lines.len() {
        let line = lines[line_index];
        if line.to_ascii_uppercase().contains("K-VECTORS") {
            break;
        }
        let header = line.split_whitespace().collect::<Vec<_>>();
        let Some(entry_count) = header.get(1).and_then(|value| value.parse::<usize>().ok()) else {
            line_index += 1;
            continue;
        };
        if header
            .first()
            .and_then(|value| parse_wien2k_float(value))
            .is_none()
        {
            line_index += 1;
            continue;
        }

        for entry_line in lines.iter().skip(line_index + 1).take(entry_count) {
            let fields = entry_line.split_whitespace().collect::<Vec<_>>();
            if fields.first().and_then(|value| value.parse::<u32>().ok()) != Some(1) {
                continue;
            }
            let Some(energy_ry) = fields.get(1).and_then(|value| parse_wien2k_float(value)) else {
                continue;
            };
            let Some(de) = fields.get(2).and_then(|value| parse_wien2k_float(value)) else {
                continue;
            };
            let Some(switch) = fields.get(3).map(|value| value.to_ascii_uppercase()) else {
                continue;
            };
            if !matches!(switch.as_str(), "STOP" | "CONT") {
                continue;
            }
            candidates.push(Wien2kSocRloCandidate {
                atom_index,
                energy_ry,
                de,
                switch,
                source_file: source_file.to_string(),
            });
        }

        line_index += entry_count + 1;
        atom_index += 1;
    }

    candidates
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_wien2k_inso_with_direction_and_disabled_atoms() {
        let settings = Wien2kSocPrepareSettings {
            magnetization_direction: [1, 1, 0],
            disabled_atom_indices: vec![2, 4],
            rlo_atoms: vec![Wien2kSocRlo {
                atom_index: 1,
                energy_ry: -3.25,
                de: 0.005,
                switch: "CONT".to_string(),
            }],
            ..Wien2kSocPrepareSettings::default()
        };
        let inso = build_case_inso(&settings);
        assert!(inso.contains("1 1 0"));
        assert!(inso.contains("1 -3.25000000 0.00500000 CONT"));
        assert!(inso.contains("2 2 4"));
        assert!(inso.starts_with("WFFIL"));
    }

    #[test]
    fn rejects_zero_direction_and_overlapping_rlo() {
        let zero = Wien2kSocPrepareSettings {
            magnetization_direction: [0, 0, 0],
            ..Wien2kSocPrepareSettings::default()
        };
        assert!(validate_soc_prepare_settings(&zero).is_err());

        let overlap = Wien2kSocPrepareSettings {
            disabled_atom_indices: vec![1],
            rlo_atoms: vec![Wien2kSocRlo {
                atom_index: 1,
                energy_ry: -3.5,
                de: 0.005,
                switch: "STOP".to_string(),
            }],
            ..Wien2kSocPrepareSettings::default()
        };
        assert!(validate_soc_prepare_settings(&overlap).is_err());
    }

    #[test]
    fn parses_all_l1_rlo_candidates_by_site() {
        let content = "\
WFFIL EF=0.5
 8.50 10 4
 0.30 5 0
 0 -3.24 0.001 STOP 1
 1 -1.28 0.002 CONT 1
 1 0.30 0.000 CONT 1
 2 0.30 0.005 CONT 1
 3 0.30 0.005 CONT 1
 0.30 3 0
 0 -1.53 0.002 CONT 1
 0 0.30 0.000 CONT 1
 1 -2.5D+00 5.0d-03 STOP 1
K-VECTORS FROM UNIT:4 -9.0 4.5 28
";
        let candidates = parse_soc_rlo_candidates(content, "case.in1c");

        assert_eq!(candidates.len(), 3);
        assert_eq!(candidates[0].atom_index, 1);
        assert_eq!(candidates[0].energy_ry, -1.28);
        assert_eq!(candidates[1].energy_ry, 0.30);
        assert_eq!(candidates[2].atom_index, 2);
        assert_eq!(candidates[2].energy_ry, -2.5);
        assert_eq!(candidates[2].de, 0.005);
        assert_eq!(candidates[2].switch, "STOP");
        assert_eq!(candidates[2].source_file, "case.in1c");
    }

    #[test]
    fn staged_session_serializes_empty_artifacts_for_frontend_rendering() {
        let session = Wien2kSocSession {
            session_id: "session".to_string(),
            project_id: "project".to_string(),
            cif_id: "cif".to_string(),
            source_scf_calculation_id: "scf".to_string(),
            source_structure_calculation_id: "structure".to_string(),
            case_name: "case".to_string(),
            remote_case_dir: "/remote/case".to_string(),
            source_remote_case_dir: "/remote/source".to_string(),
            remote_install_root: "/remote/wien2k".to_string(),
            hpc_profile_id: "profile".to_string(),
            spin_mode: Wien2kSpinMode::NonSpinPolarized,
            source_k_mesh: None,
            source_dft_u: false,
            phase: Wien2kSocSessionPhase::Staged,
            latest_prepare: None,
            latest_run: None,
            latest_calculation_id: None,
            artifacts: BTreeMap::new(),
            transcript: Vec::new(),
            started_at: "2026-06-10T00:00:00Z".to_string(),
        };

        let serialized = serde_json::to_value(session).unwrap();
        assert_eq!(serialized["artifacts"], serde_json::json!({}));
    }
}
