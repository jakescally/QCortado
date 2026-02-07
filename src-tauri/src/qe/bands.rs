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
            let e: f64 = parts[1].parse().map_err(|_| "Failed to parse energy value")?;
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
    if k_path.is_empty() {
        return;
    }

    let total_k = data.k_points.last().copied().unwrap_or(0.0);

    // Calculate total npoints to determine proportional positions
    // Each point's npoints indicates the segment AFTER it, so we sum all but use
    // cumulative sums to place markers
    let total_npoints: u32 = k_path.iter().map(|p| p.npoints).sum();

    if total_npoints == 0 {
        // Fallback: just place markers evenly
        let markers: Vec<HighSymmetryMarker> = k_path
            .iter()
            .enumerate()
            .map(|(i, point)| {
                let k_frac = if k_path.len() > 1 {
                    i as f64 / (k_path.len() - 1) as f64
                } else {
                    0.0
                };
                HighSymmetryMarker {
                    k_distance: k_frac * total_k,
                    label: point.label.clone(),
                }
            })
            .collect();
        data.high_symmetry_points = markers;
        return;
    }

    // Place markers based on cumulative npoints
    // First point is at k=0, subsequent points are placed proportionally
    let mut markers = Vec::new();
    let mut cumulative: u32 = 0;

    for (i, point) in k_path.iter().enumerate() {
        let k_distance = if i == 0 {
            0.0
        } else {
            (cumulative as f64 / total_npoints as f64) * total_k
        };

        markers.push(HighSymmetryMarker {
            k_distance,
            label: point.label.clone(),
        });

        cumulative += point.npoints;
    }

    // Ensure last marker is exactly at total_k
    if let Some(last) = markers.last_mut() {
        last.k_distance = total_k;
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
        assert!((data.energies[2][2] - 4.85).abs() < 0.01);    // Band 2, k-point 2
    }
}
