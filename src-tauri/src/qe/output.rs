//! Output parsing for Quantum ESPRESSO.
//!
//! Parses pw.x stdout to extract energies, convergence, forces, etc.

use super::types::QEResult;
use regex::Regex;

/// Parses the output of a pw.x calculation.
pub fn parse_pw_output(output: &str) -> QEResult {
    let mut result = QEResult {
        raw_output: output.to_string(),
        ..Default::default()
    };

    // Check convergence
    result.converged = output.contains("convergence has been achieved")
        || output.contains("End of self-consistent calculation")
        || output.contains("End of BFGS Geometry Optimization");

    // Check for explicit non-convergence
    if output.contains("convergence NOT achieved") {
        result.converged = false;
    }

    // Parse total energy
    result.total_energy = parse_total_energy(output);

    // Parse Fermi energy
    result.fermi_energy = parse_fermi_energy(output);

    // Parse total magnetization
    result.total_magnetization = parse_magnetization(output);

    // Parse forces
    result.forces = parse_forces(output);

    // Parse stress tensor
    result.stress = parse_stress(output);

    // Parse number of SCF iterations
    result.n_scf_steps = parse_scf_iterations(output);

    // Parse wall time
    result.wall_time_seconds = parse_wall_time(output);

    result
}

/// Extracts the final total energy in Ry.
fn parse_total_energy(output: &str) -> Option<f64> {
    // Pattern: "!    total energy              =     -XXX.XXXXXXXX Ry"
    // The "!" prefix indicates the final (converged) value
    let re = Regex::new(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry").ok()?;

    // Find the last match (final energy)
    re.captures_iter(output)
        .last()
        .and_then(|cap| cap.get(1))
        .and_then(|m| m.as_str().parse::<f64>().ok())
}

/// Extracts the Fermi energy in eV.
fn parse_fermi_energy(output: &str) -> Option<f64> {
    // Pattern: "the Fermi energy is    XX.XXXX ev"
    let re = Regex::new(r"the Fermi energy is\s+([-\d.]+)\s+ev").ok()?;

    re.captures(output)
        .and_then(|cap| cap.get(1))
        .and_then(|m| m.as_str().parse::<f64>().ok())
}

/// Extracts total magnetization.
fn parse_magnetization(output: &str) -> Option<f64> {
    // Pattern: "total magnetization       =     X.XX Bohr mag/cell"
    let re = Regex::new(r"total magnetization\s+=\s+([-\d.]+)").ok()?;

    re.captures_iter(output)
        .last()
        .and_then(|cap| cap.get(1))
        .and_then(|m| m.as_str().parse::<f64>().ok())
}

/// Extracts forces on atoms.
fn parse_forces(output: &str) -> Option<Vec<[f64; 3]>> {
    // Look for the forces section
    let forces_start = output.find("Forces acting on atoms")?;
    let forces_section = &output[forces_start..];

    // Pattern: "atom    1 type  1   force =     0.00000000    0.00000000    0.00000000"
    let re = Regex::new(r"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")
        .ok()?;

    let mut forces = Vec::new();

    for cap in re.captures_iter(forces_section) {
        let fx: f64 = cap.get(1)?.as_str().parse().ok()?;
        let fy: f64 = cap.get(2)?.as_str().parse().ok()?;
        let fz: f64 = cap.get(3)?.as_str().parse().ok()?;
        forces.push([fx, fy, fz]);

        // Stop if we hit the next section
        if forces_section[cap.get(0)?.end()..].starts_with("\n\n") {
            break;
        }
    }

    if forces.is_empty() {
        None
    } else {
        Some(forces)
    }
}

/// Extracts stress tensor in kbar.
fn parse_stress(output: &str) -> Option<[[f64; 3]; 3]> {
    // Pattern for stress tensor output
    // "          total   stress  (Ry/bohr**3)                   (kbar)     P=   -0.06"
    // followed by matrix elements

    let stress_start = output.find("total   stress")?;
    let stress_section = &output[stress_start..];

    // Get the three lines of stress tensor (kbar values are on the right)
    let re = Regex::new(r"[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)").ok()?;

    let mut captures_iter = re.captures_iter(stress_section);

    let row1 = captures_iter.next()?;
    let row2 = captures_iter.next()?;
    let row3 = captures_iter.next()?;

    let parse_row = |cap: regex::Captures| -> Option<[f64; 3]> {
        Some([
            cap.get(1)?.as_str().parse().ok()?,
            cap.get(2)?.as_str().parse().ok()?,
            cap.get(3)?.as_str().parse().ok()?,
        ])
    };

    Some([parse_row(row1)?, parse_row(row2)?, parse_row(row3)?])
}

/// Counts the number of SCF iterations.
fn parse_scf_iterations(output: &str) -> Option<u32> {
    // Pattern: "iteration #  X"
    let re = Regex::new(r"iteration\s+#\s*(\d+)").ok()?;

    re.captures_iter(output)
        .last()
        .and_then(|cap| cap.get(1))
        .and_then(|m| m.as_str().parse::<u32>().ok())
}

/// Extracts total wall time in seconds.
fn parse_wall_time(output: &str) -> Option<f64> {
    // Pattern: "PWSCF        :      0.50s CPU      0.52s WALL"
    // or "PWSCF        :   1m23.45s CPU   1m25.67s WALL"

    let re = Regex::new(r"PWSCF\s+:\s+[\d.hms]+\s+CPU\s+([\d.hms]+)\s+WALL").ok()?;

    let cap = re.captures(output)?;
    let wall_str = cap.get(1)?.as_str();

    parse_time_string(wall_str)
}

/// Parses a time string like "1m23.45s" or "0.52s" or "1h 2m 3.4s" into seconds.
fn parse_time_string(s: &str) -> Option<f64> {
    let mut total_seconds = 0.0;

    // Try to match hours
    if let Some(cap) = Regex::new(r"(\d+)h").ok()?.captures(s) {
        total_seconds += cap.get(1)?.as_str().parse::<f64>().ok()? * 3600.0;
    }

    // Try to match minutes
    if let Some(cap) = Regex::new(r"(\d+)m").ok()?.captures(s) {
        total_seconds += cap.get(1)?.as_str().parse::<f64>().ok()? * 60.0;
    }

    // Try to match seconds
    if let Some(cap) = Regex::new(r"([\d.]+)s").ok()?.captures(s) {
        total_seconds += cap.get(1)?.as_str().parse::<f64>().ok()?;
    }

    if total_seconds > 0.0 {
        Some(total_seconds)
    } else {
        None
    }
}

/// Parses band energies from bands.x or nscf output.
pub fn parse_band_energies(output: &str) -> Option<Vec<Vec<f64>>> {
    // Pattern for eigenvalues block:
    // "          k = 0.0000 0.0000 0.0000 (   XXX PWs)   bands (ev):"
    // followed by eigenvalue lines

    let mut all_bands = Vec::new();
    let re_kpoint = Regex::new(r"k\s*=.*bands \(ev\):").ok()?;
    let re_energies = Regex::new(r"([-\d.]+)").ok()?;

    for mat in re_kpoint.find_iter(output) {
        let start = mat.end();
        // Find the end of this k-point's eigenvalues (next k-point or empty line)
        let remaining = &output[start..];
        let end = remaining
            .find("\n\n")
            .or_else(|| remaining.find("k ="))
            .unwrap_or(remaining.len());

        let eigenvalue_block = &remaining[..end];
        let mut bands: Vec<f64> = Vec::new();

        for cap in re_energies.captures_iter(eigenvalue_block) {
            if let Ok(val) = cap.get(1).unwrap().as_str().parse::<f64>() {
                bands.push(val);
            }
        }

        if !bands.is_empty() {
            all_bands.push(bands);
        }
    }

    if all_bands.is_empty() {
        None
    } else {
        Some(all_bands)
    }
}

/// Parses DOS data from dos.x output file.
pub fn parse_dos_file(content: &str) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut energies = Vec::new();
    let mut dos_values = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        // Skip comment lines and empty lines
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(e), Ok(d)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                energies.push(e);
                dos_values.push(d);
            }
        }
    }

    if energies.is_empty() {
        None
    } else {
        Some((energies, dos_values))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_total_energy() {
        let output = r#"
     iteration #  1     ecut=    20.00 Ry     beta= 0.70
     total energy              =     -15.12345678 Ry

     iteration #  2     ecut=    20.00 Ry     beta= 0.70
!    total energy              =     -15.85432100 Ry
"#;
        let energy = parse_total_energy(output);
        assert!((energy.unwrap() - (-15.854321)).abs() < 1e-6);
    }

    #[test]
    fn test_parse_time_string() {
        assert!((parse_time_string("1m23.45s").unwrap() - 83.45).abs() < 0.01);
        assert!((parse_time_string("0.52s").unwrap() - 0.52).abs() < 0.01);
        assert!((parse_time_string("1h2m3.4s").unwrap() - 3723.4).abs() < 0.01);
    }

    #[test]
    fn test_parse_convergence() {
        let converged = "End of self-consistent calculation";
        let not_converged = "convergence NOT achieved after 100 iterations";

        let result1 = parse_pw_output(converged);
        let result2 = parse_pw_output(not_converged);

        assert!(result1.converged);
        assert!(!result2.converged);
    }
}
