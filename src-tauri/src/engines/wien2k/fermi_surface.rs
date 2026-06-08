use serde::{Deserialize, Serialize};

use super::bands::Wien2kBandsSpinChannel;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kFermiSurfaceSettings {
    pub k_mesh: [u32; 3],
    pub spin_channel: Wien2kBandsSpinChannel,
    #[serde(default)]
    pub spin_orbit: bool,
    #[serde(default)]
    pub diagnostic_log: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kSurfaceFileData {
    pub file_name: String,
    pub size_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kFermiSurfaceResult {
    pub calculation_id: String,
    pub k_grid: [u32; 3],
    pub fermi_energy: Option<f64>,
    pub primary_file: String,
    pub bxsf_files: Vec<Wien2kSurfaceFileData>,
    pub native_output: String,
    pub diagnostics: Vec<String>,
}

pub fn validate_fermi_surface_settings(settings: &Wien2kFermiSurfaceSettings) -> Result<(), String> {
    if settings.k_mesh.iter().any(|value| *value == 0) {
        return Err("Fermi-surface k mesh values must be positive.".to_string());
    }
    if settings.k_mesh.iter().any(|value| *value > 500) {
        return Err("Fermi-surface k mesh values must be 500 or lower.".to_string());
    }
    Ok(())
}

pub fn build_regular_klist(case_name: &str, k_mesh: [u32; 3]) -> String {
    let denominator = lcm3(k_mesh[0], k_mesh[1], k_mesh[2]).max(1);
    let mut output = String::new();
    let mut index = 1_u32;
    for i in 0..k_mesh[0] {
        for j in 0..k_mesh[1] {
            for k in 0..k_mesh[2] {
                let kx = ((i as u64 * denominator as u64) / k_mesh[0] as u64) as i32;
                let ky = ((j as u64 * denominator as u64) / k_mesh[1] as u64) as i32;
                let kz = ((k as u64 * denominator as u64) / k_mesh[2] as u64) as i32;
                let label = if index == 1 {
                    case_name.chars().take(10).collect::<String>()
                } else {
                    index.to_string()
                };
                output.push_str(&format!(
                    "{:<10}{:>5}{:>5}{:>5}{:>5}{:>5.1}\n",
                    label, kx, ky, kz, denominator, 1.0_f64
                ));
                index += 1;
            }
        }
    }
    output.push_str("END\n");
    output
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    a
}

fn lcm(a: u32, b: u32) -> u32 {
    if a == 0 || b == 0 {
        0
    } else {
        a / gcd(a, b) * b
    }
}

fn lcm3(a: u32, b: u32, c: u32) -> u32 {
    lcm(lcm(a, b), c)
}
