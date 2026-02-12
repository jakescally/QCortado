//! Band structure calculation support for Quantum ESPRESSO.
//!
//! This module provides:
//! - Input generation for bands.x post-processing
//! - Output parsing for bands.dat.gnu file
//! - Types for band structure data

use serde::{Deserialize, Serialize};
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write;
use std::path::Path;

/// A point in the k-path for band structure calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KPathPoint {
    pub label: String,
    pub coords: [f64; 3],
    pub npoints: u32, // 0 for last point in a segment
}

/// Band gap information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BandGap {
    /// Band gap value in eV
    pub value: f64,
    /// Whether the gap is direct
    pub is_direct: bool,
    /// K-point distance of valence band maximum
    pub vbm_k: f64,
    /// K-point distance of conduction band minimum
    pub cbm_k: f64,
    /// Energy of VBM in eV (relative to Fermi)
    pub vbm_energy: f64,
    /// Energy of CBM in eV (relative to Fermi)
    pub cbm_energy: f64,
}

/// High-symmetry point marker for the band plot
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HighSymmetryMarker {
    /// K-path distance
    pub k_distance: f64,
    /// Label (e.g., "Γ", "X", "M")
    pub label: String,
}

/// Complete band structure data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BandData {
    /// K-path distances for each k-point
    pub k_points: Vec<f64>,
    /// Energies for each band at each k-point: [band_index][k_index]
    pub energies: Vec<Vec<f64>>,
    /// Fermi energy in eV
    pub fermi_energy: f64,
    /// High-symmetry point markers
    pub high_symmetry_points: Vec<HighSymmetryMarker>,
    /// Number of bands
    pub n_bands: usize,
    /// Number of k-points
    pub n_kpoints: usize,
    /// Band gap information (if semiconductor/insulator)
    pub band_gap: Option<BandGap>,
    /// Energy range [min, max] in eV
    pub energy_range: [f64; 2],
    /// Optional orbital projection data for fat-band visualization
    #[serde(default)]
    pub projections: Option<BandProjectionData>,
}

/// Configuration for bands.x post-processing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BandsXConfig {
    pub prefix: String,
    pub outdir: String,
    pub filband: String,
    /// Whether to print symmetry labels
    pub lsym: bool,
}

/// projwfc.x configuration for orbital projection analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjwfcConfig {
    pub prefix: String,
    pub outdir: String,
    /// Output file containing k/band-resolved projection weights
    pub filproj: String,
    /// Symmetrize projections
    #[serde(default = "default_projwfc_lsym")]
    pub lsym: bool,
    /// Rotate wavefunctions to local crystal-field basis
    #[serde(default)]
    pub diag_basis: bool,
    /// Use PAW-specific projection treatment
    #[serde(default)]
    pub pawproj: bool,
}

fn default_projwfc_lsym() -> bool {
    false
}

impl Default for ProjwfcConfig {
    fn default() -> Self {
        Self {
            prefix: "qcortado_scf".to_string(),
            outdir: "./tmp".to_string(),
            filproj: "bands.projwfc.dat".to_string(),
            lsym: false,
            diag_basis: false,
            pawproj: false,
        }
    }
}

/// A projection channel for fat-band rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BandProjectionGroup {
    pub id: String,
    pub label: String,
    /// "atom" or "orbital"
    pub kind: String,
    /// Weight matrix [band][k-point]
    pub weights: Vec<Vec<f64>>,
}

/// Collection of projection groups parsed from projwfc.x outputs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BandProjectionData {
    /// Data source used to build these groups.
    pub source: String,
    pub atom_groups: Vec<BandProjectionGroup>,
    pub orbital_groups: Vec<BandProjectionGroup>,
}

impl Default for BandsXConfig {
    fn default() -> Self {
        Self {
            prefix: "qcortado".to_string(),
            outdir: "./tmp".to_string(),
            filband: "bands.dat".to_string(),
            lsym: true,
        }
    }
}

/// Generates input for bands.x post-processing
pub fn generate_bands_x_input(config: &BandsXConfig) -> String {
    let mut output = String::new();

    writeln!(output, "&BANDS").unwrap();
    writeln!(output, "  prefix = '{}',", config.prefix).unwrap();
    writeln!(output, "  outdir = '{}',", config.outdir).unwrap();
    writeln!(output, "  filband = '{}',", config.filband).unwrap();
    writeln!(
        output,
        "  lsym = .{}.,",
        if config.lsym { "true" } else { "false" }
    )
    .unwrap();
    writeln!(output, "/").unwrap();

    output
}

/// Generates input for projwfc.x projection analysis.
pub fn generate_projwfc_input(config: &ProjwfcConfig) -> String {
    let mut output = String::new();
    writeln!(output, "&PROJWFC").unwrap();
    writeln!(output, "  prefix = '{}',", config.prefix).unwrap();
    writeln!(output, "  outdir = '{}',", config.outdir).unwrap();
    writeln!(output, "  filproj = '{}',", config.filproj).unwrap();
    writeln!(
        output,
        "  lsym = .{}.,",
        if config.lsym { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  diag_basis = .{}.,",
        if config.diag_basis { "true" } else { "false" }
    )
    .unwrap();
    writeln!(
        output,
        "  pawproj = .{}.,",
        if config.pawproj { "true" } else { "false" }
    )
    .unwrap();
    writeln!(output, "/").unwrap();
    output
}

/// Parses the bands.dat.gnu file output from bands.x
///
/// The format is:
/// ```text
/// k_distance  energy  (band 1, k-point 1)
/// k_distance  energy  (band 1, k-point 2)
/// ...
/// k_distance  energy  (band 1, k-point N)
/// (blank line)
/// k_distance  energy  (band 2, k-point 1)
/// k_distance  energy  (band 2, k-point 2)
/// ...
/// ```
/// Each block separated by blank lines represents one band across all k-points.
pub fn parse_bands_gnu(content: &str, fermi_energy: f64) -> Result<BandData, String> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.is_empty() {
        return Err("Empty bands file".to_string());
    }

    // Parse data into bands, separated by blank lines
    let mut bands: Vec<Vec<(f64, f64)>> = Vec::new(); // Each band is a vec of (k, energy) pairs
    let mut current_band: Vec<(f64, f64)> = Vec::new();

    for line in lines.iter() {
        let line = line.trim();
        if line.is_empty() {
            // Blank line marks end of a band
            if !current_band.is_empty() {
                bands.push(current_band);
                current_band = Vec::new();
            }
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let k: f64 = parts[0].parse().map_err(|_| "Failed to parse k value")?;
            let e: f64 = parts[1]
                .parse()
                .map_err(|_| "Failed to parse energy value")?;
            current_band.push((k, e));
        }
    }

    // Don't forget the last band if file doesn't end with blank line
    if !current_band.is_empty() {
        bands.push(current_band);
    }

    if bands.is_empty() {
        return Err("No valid data in bands file".to_string());
    }

    // Extract k_values from the first band (all bands should have same k-points)
    let k_values: Vec<f64> = bands[0].iter().map(|(k, _)| *k).collect();
    let n_kpoints = k_values.len();
    let n_bands = bands.len();

    if n_kpoints == 0 {
        return Err("No k-points found".to_string());
    }

    // Convert to energies[band][kpoint] format
    let mut energies: Vec<Vec<f64>> = Vec::with_capacity(n_bands);
    for band in &bands {
        let band_energies: Vec<f64> = band.iter().map(|(_, e)| *e).collect();
        // Verify this band has the same number of k-points
        if band_energies.len() != n_kpoints {
            return Err(format!(
                "Band has {} points, expected {}",
                band_energies.len(),
                n_kpoints
            ));
        }
        energies.push(band_energies);
    }

    // Calculate energy range
    let mut e_min = f64::MAX;
    let mut e_max = f64::MIN;
    for band in &energies {
        for &e in band {
            if e < e_min {
                e_min = e;
            }
            if e > e_max {
                e_max = e;
            }
        }
    }

    // Calculate band gap (if present)
    let band_gap = calculate_band_gap(&energies, &k_values, fermi_energy);

    Ok(BandData {
        k_points: k_values,
        energies,
        fermi_energy,
        high_symmetry_points: Vec::new(), // Will be set separately
        n_bands,
        n_kpoints,
        band_gap,
        energy_range: [e_min, e_max],
        projections: None,
    })
}

#[derive(Debug, Clone)]
struct StateProjectionInfo {
    atom_index: usize,
    atom_symbol: String,
    orbital: String,
}

fn orbital_label_from_l(l: Option<usize>) -> String {
    match l {
        Some(0) => "s".to_string(),
        Some(1) => "p".to_string(),
        Some(2) => "d".to_string(),
        Some(3) => "f".to_string(),
        Some(4) => "g".to_string(),
        Some(value) => format!("l{}", value),
        None => "other".to_string(),
    }
}

fn commit_projection_buffer(
    current_k: Option<usize>,
    current_band: Option<usize>,
    buffer: &mut Vec<(usize, f64)>,
    state_to_groups: &HashMap<usize, (usize, usize)>,
    atom_groups: &mut [BandProjectionGroup],
    orbital_groups: &mut [BandProjectionGroup],
    n_bands: usize,
    n_kpoints: usize,
) -> bool {
    let (Some(k_idx), Some(band_idx)) = (current_k, current_band) else {
        buffer.clear();
        return false;
    };
    if band_idx >= n_bands || k_idx >= n_kpoints {
        buffer.clear();
        return false;
    }

    let mut committed = false;
    for (state_id, weight) in buffer.iter().copied() {
        if !weight.is_finite() || weight <= 0.0 {
            continue;
        }
        if let Some((atom_group_idx, orbital_group_idx)) = state_to_groups.get(&state_id).copied() {
            atom_groups[atom_group_idx].weights[band_idx][k_idx] += weight;
            orbital_groups[orbital_group_idx].weights[band_idx][k_idx] += weight;
            committed = true;
        }
    }
    buffer.clear();
    committed
}

/// Parse projwfc projection text (from `filproj` output) into atom and orbital groups.
///
/// The parser is intentionally tolerant to formatting differences across QE versions.
pub fn parse_projwfc_projection_groups(
    content: &str,
    n_bands: usize,
    n_kpoints: usize,
) -> Result<BandProjectionData, String> {
    if content.trim().is_empty() {
        return Err("Projection output is empty".to_string());
    }

    let state_re = Regex::new(
        r"state\s*#\s*(\d+)\s*:\s*atom\s*(\d+)\s*\(\s*([A-Za-z0-9_+\-]+)\s*\).*?(?:l\s*=\s*(\d+))?",
    )
    .map_err(|e| format!("Failed to build state regex: {}", e))?;
    let k_re = Regex::new(r"^\s*k\s*=\s*[-\d.Ee+]+\s+[-\d.Ee+]+\s+[-\d.Ee+]+")
        .map_err(|e| format!("Failed to build k-point regex: {}", e))?;
    let band_re = Regex::new(r"e\(\s*(\d+)\s*\)\s*=\s*[-\d.Ee+]+\s*eV")
        .map_err(|e| format!("Failed to build band regex: {}", e))?;
    let coeff_re = Regex::new(
        r"(?:\(\s*([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*\)|([-\d.Ee+]+))\s*\*\s*\[#\s*(\d+)\s*\]",
    )
    .map_err(|e| format!("Failed to build coefficient regex: {}", e))?;

    let mut states: HashMap<usize, StateProjectionInfo> = HashMap::new();
    for caps in state_re.captures_iter(content) {
        let Some(state_id) = caps.get(1).and_then(|m| m.as_str().parse::<usize>().ok()) else {
            continue;
        };
        let Some(atom_index) = caps.get(2).and_then(|m| m.as_str().parse::<usize>().ok()) else {
            continue;
        };
        let atom_symbol = caps
            .get(3)
            .map(|m| m.as_str().trim().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "X".to_string());
        let l_value = caps.get(4).and_then(|m| m.as_str().parse::<usize>().ok());
        states.insert(
            state_id,
            StateProjectionInfo {
                atom_index,
                atom_symbol,
                orbital: orbital_label_from_l(l_value),
            },
        );
    }

    if states.is_empty() {
        return Err("No projection state definitions found in projwfc output".to_string());
    }

    let mut atom_labels: BTreeMap<usize, String> = BTreeMap::new();
    let mut orbital_labels: BTreeSet<String> = BTreeSet::new();
    for info in states.values() {
        atom_labels
            .entry(info.atom_index)
            .or_insert_with(|| info.atom_symbol.clone());
        orbital_labels.insert(info.orbital.clone());
    }

    let mut atom_groups: Vec<BandProjectionGroup> = atom_labels
        .iter()
        .map(|(atom_idx, atom_symbol)| BandProjectionGroup {
            id: format!("atom-{}-{}", atom_idx, atom_symbol),
            label: format!("Atom {} ({})", atom_idx, atom_symbol),
            kind: "atom".to_string(),
            weights: vec![vec![0.0; n_kpoints]; n_bands],
        })
        .collect();
    let mut orbital_groups: Vec<BandProjectionGroup> = orbital_labels
        .iter()
        .map(|orbital| BandProjectionGroup {
            id: format!("orbital-{}", orbital),
            label: format!("{} orbitals", orbital.to_uppercase()),
            kind: "orbital".to_string(),
            weights: vec![vec![0.0; n_kpoints]; n_bands],
        })
        .collect();

    let atom_lookup: HashMap<usize, usize> = atom_labels
        .keys()
        .enumerate()
        .map(|(idx, atom_idx)| (*atom_idx, idx))
        .collect();
    let orbital_lookup: HashMap<String, usize> = orbital_labels
        .iter()
        .enumerate()
        .map(|(idx, orbital)| (orbital.clone(), idx))
        .collect();

    let mut state_to_groups: HashMap<usize, (usize, usize)> = HashMap::new();
    for (state_id, info) in &states {
        let Some(atom_group_idx) = atom_lookup.get(&info.atom_index).copied() else {
            continue;
        };
        let Some(orbital_group_idx) = orbital_lookup.get(&info.orbital).copied() else {
            continue;
        };
        state_to_groups.insert(*state_id, (atom_group_idx, orbital_group_idx));
    }

    let mut current_k: Option<usize> = None;
    let mut current_band: Option<usize> = None;
    let mut coefficient_buffer: Vec<(usize, f64)> = Vec::new();
    let mut committed_blocks = 0usize;

    for line in content.lines() {
        if k_re.is_match(line) {
            if commit_projection_buffer(
                current_k,
                current_band,
                &mut coefficient_buffer,
                &state_to_groups,
                &mut atom_groups,
                &mut orbital_groups,
                n_bands,
                n_kpoints,
            ) {
                committed_blocks += 1;
            }
            current_k = Some(current_k.map(|value| value + 1).unwrap_or(0));
            current_band = None;
            continue;
        }

        if let Some(caps) = band_re.captures(line) {
            if commit_projection_buffer(
                current_k,
                current_band,
                &mut coefficient_buffer,
                &state_to_groups,
                &mut atom_groups,
                &mut orbital_groups,
                n_bands,
                n_kpoints,
            ) {
                committed_blocks += 1;
            }
            current_band = caps
                .get(1)
                .and_then(|m| m.as_str().parse::<usize>().ok())
                .map(|value| value.saturating_sub(1));
            continue;
        }

        if current_k.is_none() || current_band.is_none() || !line.contains("[#") {
            continue;
        }

        for caps in coeff_re.captures_iter(line) {
            let Some(state_id) = caps.get(4).and_then(|m| m.as_str().parse::<usize>().ok()) else {
                continue;
            };
            let weight = if let (Some(re), Some(im)) = (
                caps.get(1).and_then(|m| m.as_str().parse::<f64>().ok()),
                caps.get(2).and_then(|m| m.as_str().parse::<f64>().ok()),
            ) {
                re * re + im * im
            } else {
                caps.get(3)
                    .and_then(|m| m.as_str().parse::<f64>().ok())
                    .map(|value| value.abs())
                    .unwrap_or(0.0)
            };
            if weight > 0.0 {
                coefficient_buffer.push((state_id, weight));
            }
        }
    }

    if commit_projection_buffer(
        current_k,
        current_band,
        &mut coefficient_buffer,
        &state_to_groups,
        &mut atom_groups,
        &mut orbital_groups,
        n_bands,
        n_kpoints,
    ) {
        committed_blocks += 1;
    }

    if committed_blocks == 0 {
        return Err("No k-point/band projection blocks were parsed from projwfc output".to_string());
    }

    Ok(BandProjectionData {
        source: "projwfc".to_string(),
        atom_groups,
        orbital_groups,
    })
}

/// Calculate band gap from band energies
fn calculate_band_gap(
    energies: &[Vec<f64>],
    k_points: &[f64],
    fermi_energy: f64,
) -> Option<BandGap> {
    // Find highest occupied and lowest unoccupied bands
    // Consider bands with max energy < fermi as occupied
    // Consider bands with min energy > fermi as unoccupied

    let mut vbm = f64::MIN;
    let mut vbm_k = 0.0;
    let mut cbm = f64::MAX;
    let mut cbm_k = 0.0;

    for band in energies.iter() {
        let band_max = band.iter().cloned().fold(f64::MIN, f64::max);
        let band_min = band.iter().cloned().fold(f64::MAX, f64::min);

        // Is this band occupied? (max energy below Fermi)
        if band_max <= fermi_energy + 0.01 {
            // This is an occupied band, check if it has the VBM
            for (k_idx, &e) in band.iter().enumerate() {
                if e > vbm {
                    vbm = e;
                    vbm_k = k_points[k_idx];
                }
            }
        }

        // Is this band unoccupied? (min energy above Fermi)
        if band_min >= fermi_energy - 0.01 && band_max > fermi_energy {
            // This is an unoccupied band, check if it has the CBM
            for (k_idx, &e) in band.iter().enumerate() {
                if e < cbm && e > fermi_energy - 0.01 {
                    cbm = e;
                    cbm_k = k_points[k_idx];
                }
            }
        }
    }

    // Only report band gap if both VBM and CBM were found
    if vbm > f64::MIN / 2.0 && cbm < f64::MAX / 2.0 && cbm > vbm {
        let gap = cbm - vbm;
        // Only report gap if it's meaningful (> 0.01 eV)
        if gap > 0.01 {
            let is_direct = (vbm_k - cbm_k).abs() < 0.01;
            return Some(BandGap {
                value: gap,
                is_direct,
                vbm_k,
                cbm_k,
                vbm_energy: vbm,
                cbm_energy: cbm,
            });
        }
    }

    None
}

/// Read bands.dat.gnu file from disk
pub fn read_bands_gnu_file(path: &Path, fermi_energy: f64) -> Result<BandData, String> {
    let content =
        std::fs::read_to_string(path).map_err(|e| format!("Failed to read bands file: {}", e))?;
    parse_bands_gnu(&content, fermi_energy)
}

/// Add high-symmetry point markers to band data
///
/// The k-path defines segments like: L(20) -> Γ(20) -> X(20) -> U(0)
/// where each npoints value indicates how many k-points in the segment TO that point.
/// High-symmetry points occur at cumulative indices: 0, 20, 40, 60.
/// We read the actual k-distance from the parsed band data at those indices.
pub fn add_symmetry_markers(data: &mut BandData, k_path: &[KPathPoint]) {
    if k_path.is_empty() || data.k_points.is_empty() {
        return;
    }

    let mut markers = Vec::new();
    let mut k_index: usize = 0;

    for (i, point) in k_path.iter().enumerate() {
        // Get the actual k-distance from the band data at this index
        let k_distance = if k_index < data.k_points.len() {
            data.k_points[k_index]
        } else {
            // Fallback to last k-point if we're past the end
            data.k_points.last().copied().unwrap_or(0.0)
        };

        markers.push(HighSymmetryMarker {
            k_distance,
            label: point.label.clone(),
        });

        // Move to the next high-symmetry point index
        // npoints indicates how many points in the segment AFTER this point
        if i < k_path.len() - 1 {
            k_index += point.npoints as usize;
        }
    }

    data.high_symmetry_points = markers;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_bands_x_input() {
        let config = BandsXConfig {
            prefix: "silicon".to_string(),
            outdir: "./tmp".to_string(),
            filband: "si_bands.dat".to_string(),
            lsym: true,
        };

        let input = generate_bands_x_input(&config);

        assert!(input.contains("prefix = 'silicon'"));
        assert!(input.contains("outdir = './tmp'"));
        assert!(input.contains("filband = 'si_bands.dat'"));
        assert!(input.contains("lsym = .true."));
    }

    #[test]
    fn test_parse_bands_gnu() {
        // Correct format: each band is a block of k-points, separated by blank lines
        let content = r#"
0.000000  -6.12
0.100000  -6.00
0.200000  -5.80

0.000000  -1.23
0.100000  -1.10
0.200000  -0.90

0.000000   4.56
0.100000   4.70
0.200000   4.85
"#;

        let result = parse_bands_gnu(content, 0.0);
        assert!(result.is_ok());

        let data = result.unwrap();
        assert_eq!(data.n_kpoints, 3);
        assert_eq!(data.n_bands, 3);
        assert_eq!(data.k_points.len(), 3);

        // Verify the energies are correctly organized
        assert!((data.energies[0][0] - (-6.12)).abs() < 0.01); // Band 0, k-point 0
        assert!((data.energies[1][0] - (-1.23)).abs() < 0.01); // Band 1, k-point 0
        assert!((data.energies[2][2] - 4.85).abs() < 0.01); // Band 2, k-point 2
    }

    #[test]
    fn test_generate_projwfc_input() {
        let input = generate_projwfc_input(&ProjwfcConfig {
            prefix: "silicon".to_string(),
            outdir: "./tmp".to_string(),
            filproj: "bands.projwfc.dat".to_string(),
            lsym: false,
            diag_basis: false,
            pawproj: false,
        });
        assert!(input.contains("&PROJWFC"));
        assert!(input.contains("prefix = 'silicon'"));
        assert!(input.contains("filproj = 'bands.projwfc.dat'"));
        assert!(input.contains("lsym = .false."));
    }

    #[test]
    fn test_parse_projwfc_projection_groups() {
        let content = r#"
state #   1: atom   1 (Si ) wfc  1 (l=0 m=1)
state #   2: atom   2 (B  ) wfc  2 (l=1 m=1)

 k = 0.0000 0.0000 0.0000
 e( 1 ) = -1.0000 eV
 psi = 0.8*[# 1] + 0.6*[# 2]
 e( 2 ) = 2.0000 eV
 psi = 1.0*[# 2]

 k = 0.5000 0.0000 0.0000
 e( 1 ) = -0.5000 eV
 psi = 1.0*[# 1]
 e( 2 ) = 1.5000 eV
 psi = (0.0,0.5)*[# 2]
"#;

        let parsed = parse_projwfc_projection_groups(content, 2, 2).unwrap();
        assert_eq!(parsed.atom_groups.len(), 2);
        assert_eq!(parsed.orbital_groups.len(), 2);

        // Atom 1 has weight on band 1 at both k-points.
        let atom1 = parsed
            .atom_groups
            .iter()
            .find(|group| group.id.contains("atom-1"))
            .unwrap();
        assert!(atom1.weights[0][0] > 0.0);
        assert!(atom1.weights[0][1] > 0.0);

        // Orbital p channel should have contributions on both bands.
        let p_orbital = parsed
            .orbital_groups
            .iter()
            .find(|group| group.id == "orbital-p")
            .unwrap();
        assert!(p_orbital.weights[0][0] > 0.0);
        assert!(p_orbital.weights[1][0] > 0.0);
    }
}
