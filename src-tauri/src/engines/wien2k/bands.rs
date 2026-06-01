//! Native WIEN2k band-structure session and parser helpers.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::engines::qe::bands::{BandData, BandGap, HighSymmetryMarker};
use crate::engines::types::EngineId;

use super::Wien2kSpinMode;

const EV_TO_RY: f64 = 1.0 / 13.605_693_122_994;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kBandsSessionPhase {
    Staged,
    Prepared,
    BandsComplete,
    Failed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kBandsSpinChannel {
    None,
    Up,
    Down,
}

impl Wien2kBandsSpinChannel {
    pub fn from_spin_mode(spin_mode: Wien2kSpinMode) -> Self {
        match spin_mode {
            Wien2kSpinMode::NonSpinPolarized => Self::None,
            Wien2kSpinMode::SpinPolarized => Self::Up,
        }
    }

    pub fn x_arg(self) -> Option<&'static str> {
        match self {
            Self::None => None,
            Self::Up => Some("-up"),
            Self::Down => Some("-dn"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kKPathPoint {
    pub label: String,
    pub coords: [f64; 3],
    pub npoints: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kBandsPrepareSettings {
    pub k_path: Vec<Wien2kKPathPoint>,
    pub energy_min_ev: f64,
    pub energy_max_ev: f64,
    #[serde(default)]
    pub character_atom: u16,
    #[serde(default)]
    pub character_l: u16,
    #[serde(default = "default_character_scale")]
    pub character_scale: f64,
    #[serde(default)]
    pub run_lapw2_qtl: bool,
    #[serde(default)]
    pub run_irrep: bool,
    #[serde(default)]
    pub spin_channel: Option<Wien2kBandsSpinChannel>,
}

fn default_character_scale() -> f64 {
    0.2
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kBandsRunSettings {
    pub spin_channel: Wien2kBandsSpinChannel,
    #[serde(default)]
    pub run_lapw2_qtl: bool,
    #[serde(default)]
    pub run_irrep: bool,
    #[serde(default)]
    pub spin_orbit: bool,
    #[serde(default)]
    pub diagnostic_log: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kBandsSession {
    pub session_id: String,
    pub project_id: String,
    pub cif_id: String,
    pub source_scf_calculation_id: String,
    pub case_name: String,
    pub remote_case_dir: String,
    pub source_remote_case_dir: String,
    pub remote_install_root: String,
    pub hpc_profile_id: String,
    pub spin_mode: Wien2kSpinMode,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fermi_energy_ev: Option<f64>,
    pub phase: Wien2kBandsSessionPhase,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latest_prepare: Option<Wien2kBandsPrepareSettings>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub artifacts: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transcript: Vec<String>,
    pub started_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kBandsPrepareResult {
    pub session_id: String,
    pub phase: Wien2kBandsSessionPhase,
    pub native_output: String,
    pub artifacts: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kBandsExecutionResult {
    pub session_id: String,
    pub phase: Wien2kBandsSessionPhase,
    pub native_output: String,
    pub diagnostics: Vec<String>,
    pub band_data: BandData,
    pub band_dataset: serde_json::Value,
    pub calculation_id: String,
}

pub fn validate_prepare_settings(settings: &Wien2kBandsPrepareSettings) -> Result<(), String> {
    if settings.k_path.len() < 2 {
        return Err("Select at least two k-path points.".to_string());
    }
    if settings.k_path.iter().all(|point| point.npoints == 0) {
        return Err("The k-path must contain at least one sampled segment.".to_string());
    }
    if !settings.energy_min_ev.is_finite()
        || !settings.energy_max_ev.is_finite()
        || settings.energy_min_ev >= settings.energy_max_ev
    {
        return Err("Energy window must be finite and ordered.".to_string());
    }
    if !settings.character_scale.is_finite() || settings.character_scale < 0.0 {
        return Err("Character scale must be non-negative.".to_string());
    }
    Ok(())
}

pub fn build_klist_band(_case_name: &str, path: &[Wien2kKPathPoint]) -> String {
    let expanded = expand_k_path(path);
    let denominator = 10_000_i32;
    let mut output = String::new();
    for point in expanded {
        let label = sanitize_k_label(&point.label);
        let kx = (point.coords[0] * denominator as f64).round() as i32;
        let ky = (point.coords[1] * denominator as f64).round() as i32;
        let kz = (point.coords[2] * denominator as f64).round() as i32;
        output.push_str(&format!(
            "{:<10}{:>5}{:>5}{:>5}{:>5}{:>5.1}\n",
            label, kx, ky, kz, denominator, 2.0_f64
        ));
    }
    output.push_str("END\n");
    output
}

pub fn build_insp(
    _case_name: &str,
    settings: &Wien2kBandsPrepareSettings,
    fermi_ev: Option<f64>,
) -> String {
    let fermi_ry = fermi_ev.unwrap_or(0.0) * EV_TO_RY;
    let line_switch = if settings.run_lapw2_qtl && settings.character_atom > 0 {
        0
    } else {
        1
    };
    let color_switch = if settings.run_irrep { 4 } else { 0 };
    let character_column = if settings.character_atom > 0 {
        settings.character_l.max(1)
    } else {
        1
    };
    format!(
        "### Figure configuration\n\
         5.0 3.0 # paper offset of plot\n\
         10.0 15.0 # xsize,ysize [cm]\n\
         1.0 4 # major ticks, minor ticks\n\
         1.0 1 # character height, font switch\n\
         1.1 {line_switch} {color_switch} # line width, line switch, color switch\n\
         ### Data configuration\n\
         {emin:.8} {emax:.8} 2 # energy range, energy switch (1:Ry, 2:eV)\n\
         1 {fermi:.10} # Fermi switch, Fermi-level (in Ry units)\n\
         1 999 # lower and upper band index for heavier plotting\n\
         {atom} {jcol} {scale:.6} # jatom, jtype, size of heavier plotting\n",
        emin = settings.energy_min_ev,
        emax = settings.energy_max_ev,
        fermi = fermi_ry,
        atom = settings.character_atom,
        jcol = character_column,
        scale = settings.character_scale,
    )
}

pub fn expand_k_path(path: &[Wien2kKPathPoint]) -> Vec<Wien2kKPathPoint> {
    let mut expanded = Vec::new();
    if path.is_empty() {
        return expanded;
    }

    for segment_index in 0..path.len().saturating_sub(1) {
        let from = &path[segment_index];
        let to = &path[segment_index + 1];
        let samples = from.npoints.max(1);
        for sample_index in 0..samples {
            let t = sample_index as f64 / samples as f64;
            let label = if sample_index == 0 {
                from.label.clone()
            } else {
                String::new()
            };
            expanded.push(Wien2kKPathPoint {
                label,
                coords: [
                    from.coords[0] + (to.coords[0] - from.coords[0]) * t,
                    from.coords[1] + (to.coords[1] - from.coords[1]) * t,
                    from.coords[2] + (to.coords[2] - from.coords[2]) * t,
                ],
                npoints: 0,
            });
        }
    }
    if let Some(last) = path.last() {
        expanded.push(Wien2kKPathPoint {
            label: last.label.clone(),
            coords: last.coords,
            npoints: 0,
        });
    }
    expanded
}

pub fn parse_spaghetti_xy(content: &str, fermi_energy_ev: f64) -> Result<BandData, String> {
    let mut bands: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut current: Vec<(f64, f64)> = Vec::new();
    for raw in content.lines() {
        let line = raw.trim();
        if line.is_empty() || line == "&" {
            if !current.is_empty() {
                bands.push(std::mem::take(&mut current));
            }
            continue;
        }
        if line.starts_with('#')
            || line.starts_with('@')
            || line.starts_with('*')
            || line.starts_with('!')
        {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 2 {
            continue;
        }
        let Ok(x) = fields[0].replace('D', "E").replace('d', "E").parse::<f64>() else {
            continue;
        };
        let Ok(y) = fields[1].replace('D', "E").replace('d', "E").parse::<f64>() else {
            continue;
        };
        current.push((x, y));
    }
    if !current.is_empty() {
        bands.push(current);
    }
    if bands.is_empty() {
        return Err("No x/y band traces were parsed from WIEN2k spaghetti output.".to_string());
    }
    let n_kpoints = bands[0].len();
    if n_kpoints == 0 {
        return Err("WIEN2k spaghetti output contained no k-points.".to_string());
    }
    for band in &bands {
        if band.len() != n_kpoints {
            return Err(
                "WIEN2k spaghetti output has bands with inconsistent k-point counts.".to_string(),
            );
        }
    }
    let k_points = bands[0].iter().map(|(x, _)| *x).collect::<Vec<_>>();
    let energies = bands
        .iter()
        .map(|band| band.iter().map(|(_, energy)| *energy).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut e_min = f64::INFINITY;
    let mut e_max = f64::NEG_INFINITY;
    for band in &energies {
        for energy in band {
            e_min = e_min.min(*energy);
            e_max = e_max.max(*energy);
        }
    }
    Ok(BandData {
        k_points,
        n_bands: energies.len(),
        n_kpoints,
        band_gap: calculate_band_gap(&energies, &bands[0], fermi_energy_ev),
        energy_range: [e_min, e_max],
        energies,
        fermi_energy: fermi_energy_ev,
        high_symmetry_points: Vec::new(),
        projections: None,
    })
}

pub fn add_symmetry_markers(data: &mut BandData, path: &[Wien2kKPathPoint]) {
    if path.is_empty() || data.k_points.is_empty() {
        return;
    }
    let mut markers = Vec::new();
    let mut k_index = 0_usize;
    for (index, point) in path.iter().enumerate() {
        let k_distance = data
            .k_points
            .get(k_index)
            .copied()
            .unwrap_or_else(|| data.k_points.last().copied().unwrap_or(0.0));
        markers.push(HighSymmetryMarker {
            k_distance,
            label: point.label.clone(),
        });
        if index < path.len() - 1 {
            k_index = k_index.saturating_add(point.npoints as usize);
        }
    }
    data.high_symmetry_points = markers;
}

pub fn band_dataset_json(
    band_data: &BandData,
    calculation_id: Option<&str>,
    project_id: &str,
    cif_id: &str,
    source_scf_id: &str,
    generated_at: &str,
) -> serde_json::Value {
    let series = band_data
        .energies
        .iter()
        .enumerate()
        .map(|(index, values)| {
            serde_json::json!({
                "id": format!("band-{}", index + 1),
                "label": format!("Band {}", index + 1),
                "values": values,
                "unit": "eV",
                "metadata": { "nativeBandIndex": index }
            })
        })
        .collect::<Vec<_>>();
    let markers = band_data
        .high_symmetry_points
        .iter()
        .map(|marker| serde_json::json!({ "x": marker.k_distance, "label": marker.label }))
        .collect::<Vec<_>>();
    serde_json::json!({
        "schema": "cortado.band_path.v1",
        "provenance": {
            "engineId": EngineId::Wien2k,
            "calculationKind": "bands",
            "calculationId": calculation_id,
            "projectId": project_id,
            "cifId": cif_id,
            "sourceCalculationIds": [source_scf_id],
            "generatedAt": generated_at
        },
        "quantity": "electronic_energy",
        "x": band_data.k_points,
        "series": series,
        "xUnit": "path_distance",
        "yUnit": "eV",
        "referenceEnergyEv": band_data.fermi_energy,
        "markers": markers,
        "bandGap": band_data.band_gap.as_ref().map(|gap| serde_json::json!({
            "valueEv": gap.value,
            "isDirect": gap.is_direct,
            "vbmX": gap.vbm_k,
            "cbmX": gap.cbm_k,
            "vbmEnergyEv": gap.vbm_energy,
            "cbmEnergyEv": gap.cbm_energy
        })),
        "projections": null,
        "metadata": {
            "nBands": band_data.n_bands,
            "nKpoints": band_data.n_kpoints,
            "energyRangeEv": band_data.energy_range,
            "sourceFormat": "wien2k-spaghetti"
        }
    })
}

fn sanitize_k_label(label: &str) -> String {
    let ascii = label
        .replace('Γ', "GAMMA")
        .replace('γ', "GAMMA")
        .replace('Δ', "DELTA")
        .replace('δ', "DELTA")
        .replace('Λ', "LAMBDA")
        .replace('λ', "LAMBDA")
        .replace('Σ', "SIGMA")
        .replace('σ', "SIGMA");
    ascii
        .chars()
        .filter(|character| character.is_ascii_alphanumeric())
        .take(10)
        .collect::<String>()
}

fn calculate_band_gap(
    energies: &[Vec<f64>],
    first_band_pairs: &[(f64, f64)],
    fermi_energy: f64,
) -> Option<BandGap> {
    let mut vbm_energy = f64::NEG_INFINITY;
    let mut cbm_energy = f64::INFINITY;
    let mut vbm_k = 0.0;
    let mut cbm_k = 0.0;
    for band in energies {
        for (k_index, energy) in band.iter().enumerate() {
            let k = first_band_pairs
                .get(k_index)
                .map(|(x, _)| *x)
                .unwrap_or(0.0);
            if *energy <= fermi_energy && *energy > vbm_energy {
                vbm_energy = *energy;
                vbm_k = k;
            }
            if *energy >= fermi_energy && *energy < cbm_energy {
                cbm_energy = *energy;
                cbm_k = k;
            }
        }
    }
    if !vbm_energy.is_finite() || !cbm_energy.is_finite() || cbm_energy < vbm_energy {
        return None;
    }
    Some(BandGap {
        value: cbm_energy - vbm_energy,
        is_direct: (vbm_k - cbm_k).abs() < 1e-6,
        vbm_k,
        cbm_k,
        vbm_energy,
        cbm_energy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn klist_band_expands_segments_and_writes_end() {
        let path = vec![
            Wien2kKPathPoint {
                label: "Γ".to_string(),
                coords: [0.0, 0.0, 0.0],
                npoints: 2,
            },
            Wien2kKPathPoint {
                label: "X".to_string(),
                coords: [0.5, 0.0, 0.0],
                npoints: 0,
            },
        ];
        let klist = build_klist_band("Si", &path);
        let lines = klist.lines().collect::<Vec<_>>();
        assert_eq!(lines[0], "GAMMA         0    0    010000  2.0");
        assert_eq!(lines[1], "           2500    0    010000  2.0");
        assert_eq!(lines[2], "X          5000    0    010000  2.0");
        assert!(klist.trim_end().ends_with("END"));
    }

    #[test]
    fn insp_uses_modern_spaghetti_template() {
        let settings = Wien2kBandsPrepareSettings {
            k_path: vec![
                Wien2kKPathPoint {
                    label: "G".to_string(),
                    coords: [0.0, 0.0, 0.0],
                    npoints: 2,
                },
                Wien2kKPathPoint {
                    label: "X".to_string(),
                    coords: [0.5, 0.0, 0.0],
                    npoints: 0,
                },
            ],
            energy_min_ev: -8.0,
            energy_max_ev: 6.0,
            character_atom: 0,
            character_l: 0,
            character_scale: 0.2,
            run_lapw2_qtl: false,
            run_irrep: false,
            spin_channel: None,
        };
        let insp = build_insp("Si", &settings, Some(13.605_693_122_994));

        assert!(insp.starts_with("### Figure configuration\n"));
        assert!(insp.contains("### Data configuration"));
        assert!(insp.contains("-8.00000000 6.00000000 2"));
        assert!(insp.contains("1 1.0000000000"));
        assert!(insp.contains("0 1 0.200000"));
    }

    #[test]
    fn spaghetti_xy_parser_accepts_agr_blocks() {
        let content = "\
@ title \"bands\"\n\
0.0 -1.0\n\
1.0 -0.5\n\
&\n\
0.0 0.2\n\
1.0 0.6\n";
        let parsed = parse_spaghetti_xy(content, 0.0).expect("parse bands");
        assert_eq!(parsed.n_bands, 2);
        assert_eq!(parsed.n_kpoints, 2);
        assert_eq!(parsed.energies[0][1], -0.5);
    }
}
