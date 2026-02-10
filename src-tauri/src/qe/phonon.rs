//! Phonon calculation support for Quantum ESPRESSO.
//!
//! This module provides:
//! - Output parsing for ph.x
//! - Phonon DOS parsing from matdyn.x
//! - Phonon dispersion parsing from matdyn.x
//! - High-symmetry point markers for dispersion plots

use super::types::{PhononDOS, PhononDispersion, PhononHighSymmetryMarker, QPathPoint};

/// Parse ph.x output to check for convergence
///
/// Returns (converged, n_qpoints_completed)
pub fn parse_ph_output(output: &str) -> (bool, u32) {
    let mut converged = false;
    let mut n_qpoints = 0;

    for line in output.lines() {
        // Check for successful completion
        if line.contains("JOB DONE") {
            converged = true;
        }

        // Count completed q-points
        // ph.x outputs "Calculation of q =" for each q-point
        if line.contains("Calculation of q =") {
            n_qpoints += 1;
        }

        // Also check for "Representation" completion
        if line.contains("Self-consistent Calculation") {
            // This indicates phonon SCF is running
        }
    }

    (converged, n_qpoints)
}

/// Parse the phonon DOS file from matdyn.x
///
/// The format is two columns: frequency (cm^-1) and DOS
pub fn parse_phonon_dos(content: &str) -> Result<PhononDOS, String> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.is_empty() {
        return Err("Empty DOS file".to_string());
    }

    let mut frequencies = Vec::new();
    let mut dos = Vec::new();

    for line in lines.iter() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let freq: f64 = parts[0].parse().map_err(|_| "Failed to parse frequency")?;
            let d: f64 = parts[1].parse().map_err(|_| "Failed to parse DOS value")?;
            frequencies.push(freq);
            dos.push(d);
        }
    }

    if frequencies.is_empty() {
        return Err("No valid data in DOS file".to_string());
    }

    let omega_min = frequencies.iter().cloned().fold(f64::MAX, f64::min);
    let omega_max = frequencies.iter().cloned().fold(f64::MIN, f64::max);

    Ok(PhononDOS {
        frequencies,
        dos,
        omega_max,
        omega_min,
    })
}

/// Parse the phonon dispersion file from matdyn.x
///
/// The format is similar to bands.dat.gnu:
/// - Each block (separated by blank lines) represents one phonon mode
/// - Within each block: q_distance  frequency
pub fn parse_phonon_dispersion(content: &str) -> Result<PhononDispersion, String> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.is_empty() {
        return Err("Empty dispersion file".to_string());
    }

    // Detect whether this file is in:
    // 1) mode-block format (blank-line separated, two columns: q freq), or
    // 2) q-row format (one row per q, columns: q f1 f2 ... fn)
    let mut detected_columns: Option<usize> = None;
    for line in &lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let numeric_values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|part| part.parse::<f64>().ok())
            .collect();
        if numeric_values.len() >= 2 {
            detected_columns = Some(numeric_values.len());
            break;
        }
    }
    let is_q_row_format = detected_columns.map(|cols| cols > 2).unwrap_or(false);

    if is_q_row_format {
        // Parse q-row format: q f1 f2 ... fn
        let mut q_points: Vec<f64> = Vec::new();
        let mut frequencies: Vec<Vec<f64>> = Vec::new();
        let mut mode_count: Option<usize> = None;

        for line in &lines {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let numeric_values: Vec<f64> = line
                .split_whitespace()
                .filter_map(|part| part.parse::<f64>().ok())
                .collect();

            if numeric_values.len() < 2 {
                continue;
            }

            let this_mode_count = numeric_values.len() - 1;
            if this_mode_count == 0 {
                continue;
            }

            if let Some(expected) = mode_count {
                if this_mode_count < expected {
                    return Err(format!(
                        "Q-row has {} modes, expected at least {}",
                        this_mode_count, expected
                    ));
                }
            } else {
                mode_count = Some(this_mode_count);
                frequencies = vec![Vec::new(); this_mode_count];
            }

            q_points.push(numeric_values[0]);
            let expected = mode_count.unwrap_or(0);
            for mode_idx in 0..expected {
                frequencies[mode_idx].push(numeric_values[mode_idx + 1]);
            }
        }

        if q_points.is_empty() || frequencies.is_empty() {
            return Err("No valid q-row dispersion data found".to_string());
        }

        let n_qpoints = q_points.len();
        let n_modes = frequencies.len();

        let mut f_min = f64::MAX;
        let mut f_max = f64::MIN;
        for mode in &frequencies {
            for &f in mode {
                if f < f_min {
                    f_min = f;
                }
                if f > f_max {
                    f_max = f;
                }
            }
        }

        return Ok(PhononDispersion {
            q_points,
            frequencies,
            high_symmetry_points: Vec::new(),
            n_modes,
            n_qpoints,
            frequency_range: [f_min, f_max],
        });
    }

    // Parse data into modes, separated by blank lines
    let mut modes: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut current_mode: Vec<(f64, f64)> = Vec::new();

    for line in lines.iter() {
        let line = line.trim();
        if line.is_empty() {
            if !current_mode.is_empty() {
                modes.push(current_mode);
                current_mode = Vec::new();
            }
            continue;
        }

        // Skip comment lines
        if line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let q: f64 = match parts[0].parse() {
                Ok(value) => value,
                Err(_) => continue,
            };
            let freq: f64 = match parts[1].parse() {
                Ok(value) => value,
                Err(_) => continue,
            };
            current_mode.push((q, freq));
        }
    }

    // Don't forget the last mode
    if !current_mode.is_empty() {
        modes.push(current_mode);
    }

    if modes.is_empty() {
        return Err("No valid data in dispersion file".to_string());
    }

    // Extract q-points from first mode (all modes should have same q-points)
    let q_points: Vec<f64> = modes[0].iter().map(|(q, _)| *q).collect();
    let n_qpoints = q_points.len();
    let n_modes = modes.len();

    if n_qpoints == 0 {
        return Err("No q-points found".to_string());
    }

    // Convert to frequencies[mode][qpoint] format
    let mut frequencies: Vec<Vec<f64>> = Vec::with_capacity(n_modes);
    for mode in &modes {
        let mode_freqs: Vec<f64> = mode.iter().map(|(_, f)| *f).collect();
        if mode_freqs.len() != n_qpoints {
            return Err(format!(
                "Mode has {} points, expected {}",
                mode_freqs.len(),
                n_qpoints
            ));
        }
        frequencies.push(mode_freqs);
    }

    // Calculate frequency range
    let mut f_min = f64::MAX;
    let mut f_max = f64::MIN;
    for mode in &frequencies {
        for &f in mode {
            if f < f_min {
                f_min = f;
            }
            if f > f_max {
                f_max = f;
            }
        }
    }

    Ok(PhononDispersion {
        q_points,
        frequencies,
        high_symmetry_points: Vec::new(),
        n_modes,
        n_qpoints,
        frequency_range: [f_min, f_max],
    })
}

/// Add high-symmetry point markers to phonon dispersion data
///
/// Uses the same logic as band structure markers:
/// - q_path defines segments like: Γ(20) -> X(20) -> L(0)
/// - High-symmetry points occur at cumulative indices
pub fn add_phonon_symmetry_markers(data: &mut PhononDispersion, q_path: &[QPathPoint]) {
    if q_path.is_empty() || data.q_points.is_empty() {
        return;
    }

    let mut markers = Vec::new();
    let mut q_index: usize = 0;

    for (i, point) in q_path.iter().enumerate() {
        // Get the actual q-distance from the dispersion data at this index
        let q_distance = if q_index < data.q_points.len() {
            data.q_points[q_index]
        } else {
            // Fallback to last q-point if we're past the end
            data.q_points.last().copied().unwrap_or(0.0)
        };

        markers.push(PhononHighSymmetryMarker {
            q_distance,
            label: point.label.clone(),
        });

        // Move to the next high-symmetry point index
        if i < q_path.len() - 1 {
            q_index += point.npoints as usize;
        }
    }

    data.high_symmetry_points = markers;
}

/// Read phonon DOS file from disk
pub fn read_phonon_dos_file(path: &std::path::Path) -> Result<PhononDOS, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read phonon DOS file: {}", e))?;
    parse_phonon_dos(&content)
}

/// Read phonon dispersion file from disk
pub fn read_phonon_dispersion_file(path: &std::path::Path) -> Result<PhononDispersion, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read phonon dispersion file: {}", e))?;
    parse_phonon_dispersion(&content)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ph_output_converged() {
        let output = r#"
     Program PHONON v.7.0 starts on 15Jan2024 at 12:00:00

     Calculation of q =    0.0000000   0.0000000   0.0000000
     ...
     Calculation of q =    0.5000000   0.0000000   0.0000000
     ...
     JOB DONE.
        "#;

        let (converged, n_qpoints) = parse_ph_output(output);
        assert!(converged);
        assert_eq!(n_qpoints, 2);
    }

    #[test]
    fn test_parse_phonon_dos() {
        let content = r#"
# Phonon DOS
0.0    0.0
50.0   0.5
100.0  1.2
150.0  0.8
200.0  0.3
        "#;

        let result = parse_phonon_dos(content);
        assert!(result.is_ok());

        let dos = result.unwrap();
        assert_eq!(dos.frequencies.len(), 5);
        assert_eq!(dos.dos.len(), 5);
        assert!((dos.omega_max - 200.0).abs() < 0.01);
        assert!((dos.omega_min - 0.0).abs() < 0.01);
    }

    #[test]
    fn test_parse_phonon_dispersion() {
        // Simulate matdyn output with 3 modes and 3 q-points
        let content = r#"
0.000000  0.00
0.100000  50.00
0.200000  100.00

0.000000  200.00
0.100000  220.00
0.200000  240.00

0.000000  400.00
0.100000  420.00
0.200000  450.00
        "#;

        let result = parse_phonon_dispersion(content);
        assert!(result.is_ok());

        let data = result.unwrap();
        assert_eq!(data.n_modes, 3);
        assert_eq!(data.n_qpoints, 3);
        assert_eq!(data.q_points.len(), 3);

        // Check first mode values
        assert!((data.frequencies[0][0] - 0.0).abs() < 0.01);
        assert!((data.frequencies[0][2] - 100.0).abs() < 0.01);

        // Check last mode values
        assert!((data.frequencies[2][0] - 400.0).abs() < 0.01);
        assert!((data.frequencies[2][2] - 450.0).abs() < 0.01);
    }

    #[test]
    fn test_parse_phonon_dispersion_q_row_format() {
        // Simulate matdyn output with q followed by all mode frequencies on one line
        let content = r#"
0.000000  0.00  200.00  400.00
0.100000 50.00  220.00  420.00
0.200000 100.00 240.00  450.00
        "#;

        let result = parse_phonon_dispersion(content);
        assert!(result.is_ok());

        let data = result.unwrap();
        assert_eq!(data.n_modes, 3);
        assert_eq!(data.n_qpoints, 3);
        assert_eq!(data.q_points.len(), 3);
        assert!((data.frequencies[0][2] - 100.0).abs() < 0.01);
        assert!((data.frequencies[2][1] - 420.0).abs() < 0.01);
    }

    #[test]
    fn test_add_phonon_symmetry_markers() {
        let mut data = PhononDispersion {
            q_points: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
            frequencies: vec![vec![0.0; 6], vec![100.0; 6]],
            high_symmetry_points: Vec::new(),
            n_modes: 2,
            n_qpoints: 6,
            frequency_range: [0.0, 100.0],
        };

        let q_path = vec![
            QPathPoint {
                label: "Γ".to_string(),
                coords: [0.0, 0.0, 0.0],
                npoints: 3,
            },
            QPathPoint {
                label: "X".to_string(),
                coords: [0.5, 0.0, 0.5],
                npoints: 3,
            },
            QPathPoint {
                label: "L".to_string(),
                coords: [0.5, 0.5, 0.5],
                npoints: 0,
            },
        ];

        add_phonon_symmetry_markers(&mut data, &q_path);

        assert_eq!(data.high_symmetry_points.len(), 3);
        assert_eq!(data.high_symmetry_points[0].label, "Γ");
        assert_eq!(data.high_symmetry_points[1].label, "X");
        assert_eq!(data.high_symmetry_points[2].label, "L");
    }
}
