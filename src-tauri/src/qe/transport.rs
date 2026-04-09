//! Wannier90 `postw90.x` / BoltzWann transport helpers.

use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fmt::Write;
use std::path::Path;

use super::wannier::WannierArtifact;

fn default_boltz_kmesh() -> [u32; 3] {
    [24, 24, 24]
}

fn default_relaxation_time_fs() -> f64 {
    10.0
}

fn default_false() -> bool {
    false
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransportCalculationConfig {
    pub project_id: String,
    pub source_wannier_calc_id: String,
    #[serde(default = "default_boltz_kmesh")]
    pub boltz_kmesh: [u32; 3],
    pub mu_offset_min: f64,
    pub mu_offset_max: f64,
    pub mu_offset_step: f64,
    pub temp_min: f64,
    pub temp_max: f64,
    pub temp_step: f64,
    #[serde(default = "default_relaxation_time_fs")]
    pub relaxation_time_fs: f64,
    #[serde(default = "default_false")]
    pub is_2d: bool,
    #[serde(default)]
    pub boltz_2d_dir: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransportDataset {
    pub file_name: String,
    #[serde(default)]
    pub component_labels: Vec<String>,
    #[serde(default)]
    pub values: Vec<Vec<Vec<Option<f64>>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransportTdfData {
    pub file_name: String,
    #[serde(default)]
    pub energy_values_ev: Vec<f64>,
    #[serde(default)]
    pub component_labels: Vec<String>,
    #[serde(default)]
    pub values: Vec<Vec<Option<f64>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransportResult {
    pub engine: String,
    pub seedname: String,
    pub source_wannier_calc_id: String,
    pub reference_fermi_energy_ev: f64,
    #[serde(default)]
    pub mu_values_ev: Vec<f64>,
    #[serde(default)]
    pub mu_offsets_ev: Vec<f64>,
    #[serde(default)]
    pub temperature_values_k: Vec<f64>,
    pub relaxation_time_fs: f64,
    #[serde(default)]
    pub is_2d: bool,
    #[serde(default)]
    pub boltz_2d_dir: Option<String>,
    pub conductivity: TransportDataset,
    pub sigma_s: TransportDataset,
    pub seebeck: TransportDataset,
    pub kappa: TransportDataset,
    #[serde(default)]
    pub tdf: Option<TransportTdfData>,
    #[serde(default)]
    pub notes: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub artifact_manifest: Vec<WannierArtifact>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AxisOrder {
    MuTemp,
    TempMu,
}

#[derive(Debug, Clone)]
struct ParsedGridRows {
    axis_order: AxisOrder,
    rows: Vec<(f64, f64, Vec<f64>)>,
}

pub fn validate_transport_config(config: &TransportCalculationConfig) -> Result<(), String> {
    if config.project_id.trim().is_empty() {
        return Err("Transport project_id is required.".to_string());
    }
    if config.source_wannier_calc_id.trim().is_empty() {
        return Err("Transport source_wannier_calc_id is required.".to_string());
    }
    if config.boltz_kmesh.iter().any(|value| *value == 0) {
        return Err("BoltzWann k-mesh must be positive in every direction.".to_string());
    }
    if !config.mu_offset_min.is_finite()
        || !config.mu_offset_max.is_finite()
        || !config.mu_offset_step.is_finite()
    {
        return Err("Chemical-potential scan values must be finite.".to_string());
    }
    if !config.temp_min.is_finite() || !config.temp_max.is_finite() || !config.temp_step.is_finite()
    {
        return Err("Temperature scan values must be finite.".to_string());
    }
    if config.mu_offset_step <= 0.0 {
        return Err("mu_offset_step must be positive.".to_string());
    }
    if config.temp_step <= 0.0 {
        return Err("temp_step must be positive.".to_string());
    }
    if config.mu_offset_max < config.mu_offset_min {
        return Err("mu_offset_max must be greater than or equal to mu_offset_min.".to_string());
    }
    if config.temp_max < config.temp_min {
        return Err("temp_max must be greater than or equal to temp_min.".to_string());
    }
    if !config.relaxation_time_fs.is_finite() || config.relaxation_time_fs <= 0.0 {
        return Err("relaxation_time_fs must be positive.".to_string());
    }
    if config.is_2d
        && config
            .boltz_2d_dir
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .is_none()
    {
        return Err("boltz_2d_dir is required when is_2d is enabled.".to_string());
    }
    Ok(())
}

fn normalized_keyword(line: &str) -> Option<String> {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('!') || trimmed.starts_with('#') {
        return None;
    }
    let (lhs, _) = trimmed.split_once('=')?;
    let key = lhs.trim().to_ascii_lowercase();
    if key.is_empty() {
        None
    } else {
        Some(key)
    }
}

fn format_decimal(value: f64) -> String {
    let formatted = format!("{:.12}", value);
    formatted
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

pub fn build_transport_win(
    source_win_content: &str,
    config: &TransportCalculationConfig,
    reference_fermi_energy_ev: f64,
) -> Result<String, String> {
    validate_transport_config(config)?;

    let mu_min = reference_fermi_energy_ev + config.mu_offset_min;
    let mu_max = reference_fermi_energy_ev + config.mu_offset_max;
    let mut filtered_lines: Vec<&str> = Vec::new();
    let transport_keywords: HashSet<&str> = [
        "boltzwann",
        "boltz_kmesh",
        "boltz_mu_min",
        "boltz_mu_max",
        "boltz_mu_step",
        "boltz_temp_min",
        "boltz_temp_max",
        "boltz_temp_step",
        "boltz_relax_time",
        "boltz_2d_dir",
    ]
    .into_iter()
    .collect();

    for line in source_win_content.lines() {
        if normalized_keyword(line)
            .as_deref()
            .map(|keyword| transport_keywords.contains(keyword))
            .unwrap_or(false)
        {
            continue;
        }
        filtered_lines.push(line);
    }

    let mut output = filtered_lines.join("\n");
    if !output.is_empty() && !output.ends_with('\n') {
        output.push('\n');
    }
    writeln!(output).unwrap();
    writeln!(output, "! QCortado transport block").unwrap();
    writeln!(output, "boltzwann = true").unwrap();
    writeln!(
        output,
        "boltz_kmesh = {} {} {}",
        config.boltz_kmesh[0], config.boltz_kmesh[1], config.boltz_kmesh[2]
    )
    .unwrap();
    writeln!(output, "boltz_mu_min = {}", format_decimal(mu_min)).unwrap();
    writeln!(output, "boltz_mu_max = {}", format_decimal(mu_max)).unwrap();
    writeln!(
        output,
        "boltz_mu_step = {}",
        format_decimal(config.mu_offset_step)
    )
    .unwrap();
    writeln!(
        output,
        "boltz_temp_min = {}",
        format_decimal(config.temp_min)
    )
    .unwrap();
    writeln!(
        output,
        "boltz_temp_max = {}",
        format_decimal(config.temp_max)
    )
    .unwrap();
    writeln!(
        output,
        "boltz_temp_step = {}",
        format_decimal(config.temp_step)
    )
    .unwrap();
    writeln!(
        output,
        "boltz_relax_time = {}",
        format_decimal(config.relaxation_time_fs)
    )
    .unwrap();
    if config.is_2d {
        if let Some(dir) = config
            .boltz_2d_dir
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            writeln!(output, "boltz_2d_dir = {}", dir).unwrap();
        }
    }
    Ok(output)
}

pub fn collect_transport_artifacts(work_path: &Path, seedname: &str) -> Vec<WannierArtifact> {
    let file_names = [
        format!("{}.win", seedname),
        format!("{}.wpout", seedname),
        format!("{}.werr", seedname),
        format!("{}_elcond.dat", seedname),
        format!("{}_sigmas.dat", seedname),
        format!("{}_seebeck.dat", seedname),
        format!("{}_kappa.dat", seedname),
        format!("{}_tdf.dat", seedname),
        format!("{}.chk", seedname),
        format!("{}.eig", seedname),
        format!("{}.nnkp", seedname),
        "run.sbatch".to_string(),
        "slurm.out".to_string(),
        "slurm.err".to_string(),
    ];
    let mut artifacts = Vec::new();
    for file_name in file_names {
        let path = work_path.join(&file_name);
        if let Ok(metadata) = std::fs::metadata(&path) {
            if metadata.is_file() {
                artifacts.push(WannierArtifact {
                    file_name,
                    size_bytes: metadata.len(),
                });
            }
        }
    }
    artifacts
}

fn default_component_labels(count: usize) -> Vec<String> {
    match count {
        1 => vec!["value".to_string()],
        3 => vec!["xx".to_string(), "yy".to_string(), "zz".to_string()],
        6 => vec![
            "xx".to_string(),
            "xy".to_string(),
            "yy".to_string(),
            "xz".to_string(),
            "yz".to_string(),
            "zz".to_string(),
        ],
        9 => vec![
            "xx".to_string(),
            "xy".to_string(),
            "xz".to_string(),
            "yx".to_string(),
            "yy".to_string(),
            "yz".to_string(),
            "zx".to_string(),
            "zy".to_string(),
            "zz".to_string(),
        ],
        _ => (0..count).map(|idx| format!("c{}", idx + 1)).collect(),
    }
}

fn approx_eq(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= scale * 1.0e-9
}

fn push_unique_axis(values: &mut Vec<f64>, value: f64) {
    if values.iter().any(|existing| approx_eq(*existing, value)) {
        return;
    }
    values.push(value);
}

fn axis_index(values: &[f64], value: f64) -> Option<usize> {
    values
        .iter()
        .position(|existing| approx_eq(*existing, value))
}

fn detect_axis_order(header_lines: &[String]) -> AxisOrder {
    for line in header_lines {
        let lower = line.to_ascii_lowercase();
        let mu_idx = lower.find("mu").or_else(|| lower.find("chemical"));
        let temp_idx = lower.find("temp").or_else(|| lower.find(" t("));
        if let (Some(mu_pos), Some(temp_pos)) = (mu_idx, temp_idx) {
            return if mu_pos <= temp_pos {
                AxisOrder::MuTemp
            } else {
                AxisOrder::TempMu
            };
        }
    }
    AxisOrder::MuTemp
}

fn parse_grid_rows(content: &str, file_name: &str) -> Result<ParsedGridRows, String> {
    let mut header_lines: Vec<String> = Vec::new();
    let mut rows: Vec<(f64, f64, Vec<f64>)> = Vec::new();

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        let parsed: Option<Vec<f64>> = tokens
            .iter()
            .map(|token| token.parse::<f64>().ok())
            .collect();
        if let Some(values) = parsed {
            if values.len() < 3 {
                return Err(format!(
                    "Transport table {} contains a row with fewer than 3 numeric columns.",
                    file_name
                ));
            }
            rows.push((values[0], values[1], values[2..].to_vec()));
        } else {
            header_lines.push(trimmed.to_string());
        }
    }

    if rows.is_empty() {
        return Err(format!("Transport table {} does not contain numeric rows.", file_name));
    }

    let component_count = rows[0].2.len();
    if component_count == 0 {
        return Err(format!(
            "Transport table {} does not contain any tensor columns.",
            file_name
        ));
    }
    if rows
        .iter()
        .any(|(_, _, values)| values.len() != component_count)
    {
        return Err(format!(
            "Transport table {} has inconsistent tensor column counts.",
            file_name
        ));
    }

    Ok(ParsedGridRows {
        axis_order: detect_axis_order(&header_lines),
        rows,
    })
}

fn build_dataset_from_rows(file_name: &str, parsed_rows: ParsedGridRows) -> Result<TransportDataset, String> {
    let mut mu_values_ev: Vec<f64> = Vec::new();
    let mut temperature_values_k: Vec<f64> = Vec::new();

    for (first, second, _) in &parsed_rows.rows {
        let (mu, temp) = match parsed_rows.axis_order {
            AxisOrder::MuTemp => (*first, *second),
            AxisOrder::TempMu => (*second, *first),
        };
        push_unique_axis(&mut mu_values_ev, mu);
        push_unique_axis(&mut temperature_values_k, temp);
    }

    let component_count = parsed_rows.rows[0].2.len();
    let mut values = vec![
        vec![vec![None; mu_values_ev.len()]; temperature_values_k.len()];
        component_count
    ];

    for (first, second, row_values) in parsed_rows.rows {
        let (mu, temp) = match parsed_rows.axis_order {
            AxisOrder::MuTemp => (first, second),
            AxisOrder::TempMu => (second, first),
        };
        let mu_index = axis_index(&mu_values_ev, mu)
            .ok_or_else(|| format!("Failed to locate chemical potential axis index for {}", file_name))?;
        let temp_index = axis_index(&temperature_values_k, temp)
            .ok_or_else(|| format!("Failed to locate temperature axis index for {}", file_name))?;
        for (component_index, value) in row_values.into_iter().enumerate() {
            values[component_index][temp_index][mu_index] = if value.is_finite() {
                Some(value)
            } else {
                None
            };
        }
    }

    Ok(TransportDataset {
        file_name: file_name.to_string(),
        component_labels: default_component_labels(component_count),
        values,
    })
}

fn parse_grid_dataset_file(work_path: &Path, file_name: &str) -> Result<(TransportDataset, Vec<f64>, Vec<f64>), String> {
    let path = work_path.join(file_name);
    if !path.exists() {
        return Err(format!("Required transport output not found: {}", path.display()));
    }
    let content = std::fs::read_to_string(&path)
        .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;
    let parsed = parse_grid_rows(&content, file_name)?;

    let mut mu_values_ev: Vec<f64> = Vec::new();
    let mut temperature_values_k: Vec<f64> = Vec::new();
    for (first, second, _) in &parsed.rows {
        let (mu, temp) = match parsed.axis_order {
            AxisOrder::MuTemp => (*first, *second),
            AxisOrder::TempMu => (*second, *first),
        };
        push_unique_axis(&mut mu_values_ev, mu);
        push_unique_axis(&mut temperature_values_k, temp);
    }

    let dataset = build_dataset_from_rows(file_name, parsed)?;
    Ok((dataset, mu_values_ev, temperature_values_k))
}

fn parse_optional_tdf_file(work_path: &Path, file_name: &str) -> Result<Option<TransportTdfData>, String> {
    let path = work_path.join(file_name);
    if !path.exists() {
        return Ok(None);
    }
    let content = std::fs::read_to_string(&path)
        .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        let parsed: Option<Vec<f64>> = tokens
            .iter()
            .map(|token| token.parse::<f64>().ok())
            .collect();
        if let Some(values) = parsed {
            if values.len() >= 2 {
                rows.push(values);
            }
        }
    }
    if rows.is_empty() {
        return Ok(None);
    }

    let component_count = rows[0].len() - 1;
    if rows.iter().any(|row| row.len() != component_count + 1) {
        return Err(format!(
            "Transport TDF table {} has inconsistent column counts.",
            file_name
        ));
    }

    let energy_values_ev = rows.iter().map(|row| row[0]).collect::<Vec<f64>>();
    let mut values = vec![vec![None; energy_values_ev.len()]; component_count];
    for (row_index, row) in rows.iter().enumerate() {
        for component_index in 0..component_count {
            let value = row[component_index + 1];
            values[component_index][row_index] = if value.is_finite() { Some(value) } else { None };
        }
    }

    Ok(Some(TransportTdfData {
        file_name: file_name.to_string(),
        energy_values_ev,
        component_labels: default_component_labels(component_count),
        values,
    }))
}

fn push_message(messages: &mut Vec<String>, value: &str) {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return;
    }
    if messages.iter().any(|existing| existing == trimmed) {
        return;
    }
    messages.push(trimmed.to_string());
}

fn parse_postw90_messages(work_path: &Path, seedname: &str) -> (Vec<String>, Vec<String>) {
    let mut notes = Vec::new();
    let mut warnings = Vec::new();

    for file_name in [format!("{}.wpout", seedname), format!("{}.werr", seedname)] {
        let path = work_path.join(&file_name);
        let Ok(content) = std::fs::read_to_string(&path) else {
            continue;
        };
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let lower = trimmed.to_ascii_lowercase();
            if lower.contains("warning")
                || lower.contains("error")
                || lower.contains("fatal")
                || lower.contains("failed")
            {
                push_message(&mut warnings, trimmed);
            } else if lower.contains("postw90")
                || lower.contains("boltzwann")
                || lower.contains("chemical potential")
                || lower.contains("fermi")
            {
                push_message(&mut notes, trimmed);
            }
        }
    }

    (notes, warnings)
}

pub fn parse_transport_result(
    work_path: &Path,
    config: &TransportCalculationConfig,
    seedname: &str,
    reference_fermi_energy_ev: f64,
) -> Result<TransportResult, String> {
    validate_transport_config(config)?;

    let conductivity_file = format!("{}_elcond.dat", seedname);
    let sigma_s_file = format!("{}_sigmas.dat", seedname);
    let seebeck_file = format!("{}_seebeck.dat", seedname);
    let kappa_file = format!("{}_kappa.dat", seedname);
    let tdf_file = format!("{}_tdf.dat", seedname);

    let (conductivity, mu_values_ev, temperature_values_k) =
        parse_grid_dataset_file(work_path, &conductivity_file)?;
    let (sigma_s, sigma_mu_values_ev, sigma_temperature_values_k) =
        parse_grid_dataset_file(work_path, &sigma_s_file)?;
    let (seebeck, seebeck_mu_values_ev, seebeck_temperature_values_k) =
        parse_grid_dataset_file(work_path, &seebeck_file)?;
    let (kappa, kappa_mu_values_ev, kappa_temperature_values_k) =
        parse_grid_dataset_file(work_path, &kappa_file)?;

    if sigma_mu_values_ev.len() != mu_values_ev.len()
        || seebeck_mu_values_ev.len() != mu_values_ev.len()
        || kappa_mu_values_ev.len() != mu_values_ev.len()
        || sigma_temperature_values_k.len() != temperature_values_k.len()
        || seebeck_temperature_values_k.len() != temperature_values_k.len()
        || kappa_temperature_values_k.len() != temperature_values_k.len()
    {
        return Err("Transport tables do not share a consistent mu/T grid.".to_string());
    }

    let mu_offsets_ev = mu_values_ev
        .iter()
        .map(|value| value - reference_fermi_energy_ev)
        .collect();
    let tdf = parse_optional_tdf_file(work_path, &tdf_file)?;
    let (notes, warnings) = parse_postw90_messages(work_path, seedname);

    Ok(TransportResult {
        engine: "boltzwann".to_string(),
        seedname: seedname.to_string(),
        source_wannier_calc_id: config.source_wannier_calc_id.clone(),
        reference_fermi_energy_ev,
        mu_values_ev,
        mu_offsets_ev,
        temperature_values_k,
        relaxation_time_fs: config.relaxation_time_fs,
        is_2d: config.is_2d,
        boltz_2d_dir: config.boltz_2d_dir.clone(),
        conductivity,
        sigma_s,
        seebeck,
        kappa,
        tdf,
        notes,
        warnings,
        artifact_manifest: collect_transport_artifacts(work_path, seedname),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SOURCE_WIN: &str = r#"
num_wann = 8
num_bands = 12
use_ws_distance = true
boltzwann = false
boltz_kmesh = 8 8 8
"#;

    const GRID_DATA: &str = r#"
# mu[eV] T[K] xx xy yy xz yz zz
-0.2 300 1 2 3 4 5 6
0.0 300 10 20 30 40 50 60
-0.2 600 7 8 9 10 11 12
0.0 600 70 80 90 100 110 120
"#;

    const TDF_DATA: &str = r#"
# energy xx yy zz
-1.0 1 2 3
0.0 4 5 6
"#;

    #[test]
    fn build_transport_win_rewrites_transport_keywords() {
        let config = TransportCalculationConfig {
            project_id: "p1".to_string(),
            source_wannier_calc_id: "w1".to_string(),
            boltz_kmesh: [24, 18, 12],
            mu_offset_min: -0.4,
            mu_offset_max: 0.3,
            mu_offset_step: 0.05,
            temp_min: 200.0,
            temp_max: 600.0,
            temp_step: 100.0,
            relaxation_time_fs: 12.5,
            is_2d: true,
            boltz_2d_dir: Some("z".to_string()),
        };

        let result = build_transport_win(SOURCE_WIN, &config, 5.5).unwrap();
        assert!(result.contains("boltzwann = true"));
        assert!(result.contains("boltz_kmesh = 24 18 12"));
        assert!(result.contains("boltz_mu_min = 5.1"));
        assert!(result.contains("boltz_mu_max = 5.8"));
        assert!(result.contains("boltz_relax_time = 12.5"));
        assert!(result.contains("boltz_2d_dir = z"));
        assert_eq!(result.matches("boltz_kmesh").count(), 1);
    }

    #[test]
    fn parse_transport_grid_dataset_builds_axes_and_components() {
        let temp_dir = std::env::temp_dir().join(format!("qcortado_transport_test_{}", uuid::Uuid::new_v4()));
        std::fs::create_dir_all(&temp_dir).unwrap();
        let file_name = "seed_elcond.dat";
        std::fs::write(temp_dir.join(file_name), GRID_DATA).unwrap();

        let (dataset, mu_values_ev, temp_values_k) =
            parse_grid_dataset_file(&temp_dir, file_name).unwrap();

        assert_eq!(mu_values_ev, vec![-0.2, 0.0]);
        assert_eq!(temp_values_k, vec![300.0, 600.0]);
        assert_eq!(dataset.component_labels, vec!["xx", "xy", "yy", "xz", "yz", "zz"]);
        assert_eq!(dataset.values.len(), 6);
        assert_eq!(dataset.values[0][0][0], Some(1.0));
        assert_eq!(dataset.values[5][1][1], Some(120.0));

        let _ = std::fs::remove_dir_all(&temp_dir);
    }

    #[test]
    fn parse_transport_result_collects_warnings_and_serializes() {
        let temp_dir = std::env::temp_dir().join(format!("qcortado_transport_test_{}", uuid::Uuid::new_v4()));
        std::fs::create_dir_all(&temp_dir).unwrap();
        for suffix in ["_elcond.dat", "_sigmas.dat", "_seebeck.dat", "_kappa.dat"] {
            std::fs::write(temp_dir.join(format!("seed{}", suffix)), GRID_DATA).unwrap();
        }
        std::fs::write(temp_dir.join("seed_tdf.dat"), TDF_DATA).unwrap();
        std::fs::write(
            temp_dir.join("seed.wpout"),
            "BoltzWann run started\nWarning: coarse mesh\n",
        )
        .unwrap();
        std::fs::write(temp_dir.join("seed.werr"), "error: no fatal issues\n").unwrap();

        let config = TransportCalculationConfig {
            project_id: "p1".to_string(),
            source_wannier_calc_id: "w1".to_string(),
            boltz_kmesh: [24, 24, 24],
            mu_offset_min: -0.2,
            mu_offset_max: 0.0,
            mu_offset_step: 0.2,
            temp_min: 300.0,
            temp_max: 600.0,
            temp_step: 300.0,
            relaxation_time_fs: 10.0,
            is_2d: false,
            boltz_2d_dir: None,
        };

        let result = parse_transport_result(&temp_dir, &config, "seed", 5.0).unwrap();
        assert_eq!(result.mu_offsets_ev, vec![-5.2, -5.0]);
        assert_eq!(result.mu_values_ev, vec![-0.2, 0.0]);
        assert_eq!(result.temperature_values_k, vec![300.0, 600.0]);
        assert!(result.notes.iter().any(|line| line.contains("BoltzWann")));
        assert!(result.warnings.iter().any(|line| line.contains("Warning")));
        assert!(result.tdf.is_some());
        serde_json::to_string(&result).unwrap();

        let _ = std::fs::remove_dir_all(&temp_dir);
    }
}
