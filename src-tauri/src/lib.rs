//! QCortado - A modern UI for Quantum ESPRESSO
//!
//! This is the Tauri backend providing:
//! - QE input generation and validation
//! - Process execution with streaming output
//! - Output parsing and result extraction
//! - Project management

use std::path::PathBuf;
use std::sync::Mutex;
use tauri::{AppHandle, Manager, State};

pub mod config;
pub mod projects;
pub mod qe;

use qe::{generate_pw_input, parse_pw_output, QECalculation, QEResult, QERunner};

/// Application state shared across commands.
pub struct AppState {
    /// Path to QE bin directory
    pub qe_bin_dir: Mutex<Option<PathBuf>>,
    /// Current project directory
    pub project_dir: Mutex<Option<PathBuf>>,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            qe_bin_dir: Mutex::new(None),
            project_dir: Mutex::new(None),
        }
    }
}

// ============================================================================
// Tauri Commands
// ============================================================================

/// Sets the path to the Quantum ESPRESSO bin directory.
#[tauri::command]
fn set_qe_path(app: AppHandle, path: String, state: State<AppState>) -> Result<(), String> {
    let path_buf = PathBuf::from(&path);

    // Verify pw.x exists
    if !path_buf.join("pw.x").exists() {
        return Err(format!(
            "pw.x not found in {}. Please select the QE bin directory.",
            path_buf.display()
        ));
    }

    // Update in-memory state
    *state.qe_bin_dir.lock().unwrap() = Some(path_buf);

    // Persist to disk
    config::update_qe_path(&app, Some(path))?;

    Ok(())
}

/// Gets the current QE bin directory path.
#[tauri::command]
fn get_qe_path(state: State<AppState>) -> Option<String> {
    state
        .qe_bin_dir
        .lock()
        .unwrap()
        .as_ref()
        .map(|p| p.to_string_lossy().to_string())
}

/// Checks if QE is configured and which executables are available.
#[tauri::command]
fn check_qe_executables(state: State<AppState>) -> Result<Vec<String>, String> {
    let guard = state.qe_bin_dir.lock().unwrap();
    let bin_dir = guard.as_ref().ok_or("QE path not configured")?;

    let executables = [
        "pw.x", "bands.x", "dos.x", "projwfc.x", "pp.x", "ph.x", "dynmat.x", "plotband.x",
        "neb.x", "hp.x", "turbo_lanczos.x", "xspectra.x",
    ];

    let available: Vec<String> = executables
        .iter()
        .filter(|exe| bin_dir.join(exe).exists())
        .map(|s| s.to_string())
        .collect();

    Ok(available)
}

/// Generates a pw.x input file from a calculation configuration.
#[tauri::command]
fn generate_input(calculation: QECalculation) -> Result<String, String> {
    Ok(generate_pw_input(&calculation))
}

/// Validates a calculation configuration.
#[tauri::command]
fn validate_calculation(calculation: QECalculation) -> Result<Vec<String>, String> {
    let mut warnings = Vec::new();

    // Check species count matches ntyp
    let ntyp = calculation.system.species.len();
    if ntyp == 0 {
        return Err("No atomic species defined".to_string());
    }

    // Check atom count matches nat
    let nat = calculation.system.atoms.len();
    if nat == 0 {
        return Err("No atoms defined".to_string());
    }

    // Check all atoms have valid species
    for atom in &calculation.system.atoms {
        if !calculation
            .system
            .species
            .iter()
            .any(|s| s.symbol == atom.symbol)
        {
            return Err(format!(
                "Atom '{}' has no matching species definition",
                atom.symbol
            ));
        }
    }

    // Check ecutwfc
    if calculation.system.ecutwfc <= 0.0 {
        return Err("ecutwfc must be positive".to_string());
    }
    if calculation.system.ecutwfc < 20.0 {
        warnings.push("ecutwfc < 20 Ry may give inaccurate results".to_string());
    }

    // Check ecutrho
    if let Some(ecutrho) = calculation.system.ecutrho {
        if ecutrho < calculation.system.ecutwfc {
            return Err("ecutrho cannot be less than ecutwfc".to_string());
        }
    }

    // Check convergence threshold
    if calculation.conv_thr <= 0.0 {
        return Err("conv_thr must be positive".to_string());
    }
    if calculation.conv_thr > 1e-4 {
        warnings.push("conv_thr > 1e-4 is very loose".to_string());
    }

    // Check cell parameters for ibrav=0
    if calculation.system.ibrav as i32 == 0 && calculation.system.cell_parameters.is_none() {
        return Err("CELL_PARAMETERS required when ibrav=0".to_string());
    }

    // Check celldm for ibrav != 0
    if calculation.system.ibrav as i32 != 0 && calculation.system.celldm.is_none() {
        return Err("celldm required when ibrav != 0".to_string());
    }

    Ok(warnings)
}

/// Parses QE output text and extracts results.
#[tauri::command]
fn parse_output(output: String) -> QEResult {
    parse_pw_output(&output)
}

/// Runs a pw.x calculation (blocking).
#[tauri::command]
async fn run_calculation(
    calculation: QECalculation,
    working_dir: String,
    state: State<'_, AppState>,
) -> Result<QEResult, String> {
    // Clone the path out of the guard before any await points
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };

    let runner = QERunner::new(bin_dir);
    let input = generate_pw_input(&calculation);
    let work_path = PathBuf::from(working_dir);

    // Ensure working directory exists
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    runner
        .run_pw(&input, &work_path)
        .await
        .map_err(|e| e.to_string())
}

/// Sets the current project directory.
#[tauri::command]
fn set_project_dir(path: String, state: State<AppState>) -> Result<(), String> {
    let path = PathBuf::from(&path);
    if !path.exists() {
        std::fs::create_dir_all(&path)
            .map_err(|e| format!("Failed to create project directory: {}", e))?;
    }
    *state.project_dir.lock().unwrap() = Some(path);
    Ok(())
}

/// Gets the current project directory.
#[tauri::command]
fn get_project_dir(state: State<AppState>) -> Option<String> {
    state
        .project_dir
        .lock()
        .unwrap()
        .as_ref()
        .map(|p| p.to_string_lossy().to_string())
}

/// Lists available pseudopotentials in a directory.
#[tauri::command]
fn list_pseudopotentials(pseudo_dir: String) -> Result<Vec<String>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    let mut pseudos = Vec::new();
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.ends_with(".UPF") || name.ends_with(".upf") {
            pseudos.push(name);
        }
    }
    pseudos.sort();
    Ok(pseudos)
}

/// SSSP element data from JSON
#[derive(serde::Serialize, serde::Deserialize, Clone)]
pub struct SSSPElementData {
    pub filename: String,
    pub md5: Option<String>,
    pub pseudopotential: Option<String>,
    pub cutoff_wfc: f64,
    pub cutoff_rho: f64,
}

/// Loads SSSP JSON data from the pseudo directory.
/// Looks for any file matching SSSP*.json pattern.
#[tauri::command]
fn load_sssp_data(pseudo_dir: String) -> Result<std::collections::HashMap<String, SSSPElementData>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    // Find SSSP JSON file
    let mut sssp_file: Option<PathBuf> = None;
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with("SSSP") && name.ends_with(".json") {
            sssp_file = Some(entry.path());
            break;
        }
    }

    let sssp_path = sssp_file.ok_or("No SSSP JSON file found in pseudo directory")?;
    let content = std::fs::read_to_string(&sssp_path)
        .map_err(|e| format!("Failed to read SSSP file: {}", e))?;

    let data: std::collections::HashMap<String, SSSPElementData> = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse SSSP JSON: {}", e))?;

    Ok(data)
}

// ============================================================================
// Application Entry Point
// ============================================================================

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_fs::init())
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            // Initialize projects directory on startup
            if let Err(e) = projects::ensure_projects_dir(&app.handle()) {
                eprintln!("Warning: Failed to initialize projects directory: {}", e);
            }

            // Load saved configuration
            let mut qe_bin_dir: Option<PathBuf> = None;
            match config::load_config(&app.handle()) {
                Ok(cfg) => {
                    if let Some(path) = cfg.qe_bin_dir {
                        let path_buf = PathBuf::from(&path);
                        // Only use if pw.x still exists
                        if path_buf.join("pw.x").exists() {
                            qe_bin_dir = Some(path_buf);
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Failed to load config: {}", e);
                }
            }

            // Initialize AppState with loaded config
            app.manage(AppState {
                qe_bin_dir: Mutex::new(qe_bin_dir),
                project_dir: Mutex::new(None),
            });

            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            set_qe_path,
            get_qe_path,
            check_qe_executables,
            generate_input,
            validate_calculation,
            parse_output,
            run_calculation,
            set_project_dir,
            get_project_dir,
            list_pseudopotentials,
            load_sssp_data,
            // Project management commands
            projects::list_projects,
            projects::create_project,
            projects::get_project,
            projects::add_cif_to_project,
            projects::save_calculation,
            projects::delete_project,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
