//! Band structure calculation support for Quantum ESPRESSO.
//!
//! This module provides:
//! - Input generation for bands.x post-processing
//! - Output parsing for bands.dat.gnu file
//! - Types for band structure data

use serde::{Deserialize, Serialize};
use std::fmt::Write;
use std::path::Path;

/// A point in the k-path for band structure calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KPathPoint {
    pub label: String,
    pub coords: [f64; 3],
    pub npoints: u32,  // 0 for last point in a segment
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
    writeln!(output, "  lsym = .{}.,", if config.lsym { "true" } else { "false" }).unwrap();
    writeln!(output, "/").unwrap();

    output
}

/// Parses the bands.dat.gnu file output from bands.x
///
/// The format is:
/// ```text
/// k_distance  energy_band_1
/// k_distance  energy_band_2
/// ...
/// (blank line)
/// k_distance  energy_band_1  (next k-point)
/// ...
/// ```
pub fn parse_bands_gnu(content: &str, fermi_energy: f64) -> Result<BandData, String> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.is_empty() {
        return Err("Empty bands file".to_string());
    }

    // Parse all data points
    let mut k_values: Vec<f64> = Vec::new();
    let mut all_energies: Vec<(f64, f64)> = Vec::new(); // (k, energy) pairs

    for line in lines.iter() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let k: f64 = parts[0].parse().map_err(|_| "Failed to parse k value")?;
            let e: f64 = parts[1].parse().map_err(|_| "Failed to parse energy value")?;
            all_energies.push((k, e));

            // Collect unique k values in order
            if k_values.is_empty() || (k_values.last().unwrap() - k).abs() > 1e-10 {
                k_values.push(k);
            }
        }
    }

    if k_values.is_empty() || all_energies.is_empty() {
        return Err("No valid data in bands file".to_string());
    }

    let n_kpoints = k_values.len();
    let n_bands = all_energies.len() / n_kpoints;

    if n_bands == 0 {
        return Err("Could not determine number of bands".to_string());
    }

    // Organize into bands array: energies[band][kpoint]
    let mut energies: Vec<Vec<f64>> = vec![vec![0.0; n_kpoints]; n_bands];

    for (idx, (_, e)) in all_energies.iter().enumerate() {
        let band_idx = idx / n_kpoints;
        let k_idx = idx % n_kpoints;
        if band_idx < n_bands && k_idx < n_kpoints {
            energies[band_idx][k_idx] = *e;
        }
    }

    // Calculate energy range
    let mut e_min = f64::MAX;
    let mut e_max = f64::MIN;
    for band in &energies {
        for &e in band {
            if e < e_min { e_min = e; }
            if e > e_max { e_max = e; }
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
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read bands file: {}", e))?;
    parse_bands_gnu(&content, fermi_energy)
}

/// Add high-symmetry point markers to band data
pub fn add_symmetry_markers(
    data: &mut BandData,
    k_path: &[KPathPoint],
) {
    let mut markers = Vec::new();

    // Estimate total path length from data
    let total_k = data.k_points.last().copied().unwrap_or(0.0);

    for (i, point) in k_path.iter().enumerate() {
        if i == 0 || point.npoints == 0 || k_path.get(i.saturating_sub(1)).map_or(false, |p| p.npoints == 0) {
            // This is a high-symmetry point
            // Estimate its k-distance based on position in path
            let k_frac = if k_path.len() > 1 {
                i as f64 / (k_path.len() - 1) as f64
            } else {
                0.0
            };

            markers.push(HighSymmetryMarker {
                k_distance: k_frac * total_k,
                label: point.label.clone(),
            });
        }
    }

    // If we have markers, use the actual k_points data to place them more accurately
    if !markers.is_empty() && !data.k_points.is_empty() {
        // Simple heuristic: place first at k=0, last at k_max
        if let Some(first) = markers.first_mut() {
            first.k_distance = 0.0;
        }
        if let Some(last) = markers.last_mut() {
            last.k_distance = total_k;
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
        let content = r#"
0.000000  -6.12
0.000000  -1.23
0.000000   4.56
0.100000  -6.00
0.100000  -1.10
0.100000   4.70
0.200000  -5.80
0.200000  -0.90
0.200000   4.85
"#;

        let result = parse_bands_gnu(content, 0.0);
        assert!(result.is_ok());

        let data = result.unwrap();
        assert_eq!(data.n_kpoints, 3);
        assert_eq!(data.n_bands, 3);
        assert_eq!(data.k_points.len(), 3);
    }
}
