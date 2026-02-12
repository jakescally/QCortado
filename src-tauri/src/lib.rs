//! QCortado - A modern UI for Quantum ESPRESSO
//!
//! This is the Tauri backend providing:
//! - QE input generation and validation
//! - Process execution with streaming output
//! - Output parsing and result extraction
//! - Project management

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use tauri::{AppHandle, Emitter, Manager, State};

pub mod config;
pub mod process_manager;
pub mod projects;
pub mod qe;

use process_manager::ProcessManager;

use qe::{
    add_phonon_symmetry_markers, generate_matdyn_bands_input, generate_matdyn_dos_input,
    generate_ph_input, generate_q2r_input, parse_ph_output, read_phonon_dispersion_file,
    read_phonon_dos_file, MatdynCalculation, PhononPipelineConfig, PhononResult, Q2RCalculation,
    QPathPoint,
};
use qe::{
    generate_bands_x_input, generate_projwfc_input, generate_pw_input,
    parse_projwfc_projection_groups, parse_pw_output, read_bands_gnu_file, BandData, BandsXConfig,
    KPathPoint, ProjwfcConfig, QECalculation, QEResult, QERunner,
};

/// Application state shared across commands.
pub struct AppState {
    /// Path to QE bin directory
    pub qe_bin_dir: Mutex<Option<PathBuf>>,
    /// Optional command prefix prepended before all QE launches
    pub execution_prefix: Mutex<Option<String>>,
    /// Current project directory
    pub project_dir: Mutex<Option<PathBuf>>,
    /// Background process manager
    pub process_manager: ProcessManager,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            qe_bin_dir: Mutex::new(None),
            execution_prefix: Mutex::new(None),
            project_dir: Mutex::new(None),
            process_manager: ProcessManager::new(),
        }
    }
}

fn normalize_execution_prefix(prefix: Option<String>) -> Option<String> {
    prefix.and_then(|raw| {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    })
}

fn parse_execution_prefix_tokens(prefix: Option<&str>) -> Option<Vec<String>> {
    let raw = prefix?.trim();
    if raw.is_empty() {
        return None;
    }
    let tokens: Vec<String> = raw
        .split_whitespace()
        .map(|token| token.to_string())
        .collect();
    if tokens.is_empty() {
        None
    } else {
        Some(tokens)
    }
}

fn command_basename(command: &str) -> String {
    Path::new(command)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(command)
        .to_string()
}

fn is_same_command(prefix_command: &str, program: &str) -> bool {
    command_basename(prefix_command) == command_basename(program)
}

fn tokio_command_with_prefix(
    program: impl AsRef<std::ffi::OsStr>,
    execution_prefix: Option<&str>,
) -> tokio::process::Command {
    let program_os = program.as_ref();
    let program_text = program_os.to_string_lossy().to_string();

    if let Some(tokens) = parse_execution_prefix_tokens(execution_prefix) {
        let mut command = if is_same_command(&tokens[0], &program_text) {
            tokio::process::Command::new(program_os)
        } else {
            let mut cmd = tokio::process::Command::new(&tokens[0]);
            cmd.args(tokens.iter().skip(1));
            cmd.arg(program_os);
            return cmd;
        };
        command.args(tokens.iter().skip(1));
        return command;
    }

    tokio::process::Command::new(program_os)
}

fn std_command_with_prefix(
    program: impl AsRef<std::ffi::OsStr>,
    execution_prefix: Option<&str>,
) -> std::process::Command {
    let program_os = program.as_ref();
    let program_text = program_os.to_string_lossy().to_string();

    if let Some(tokens) = parse_execution_prefix_tokens(execution_prefix) {
        let mut command = if is_same_command(&tokens[0], &program_text) {
            std::process::Command::new(program_os)
        } else {
            let mut cmd = std::process::Command::new(&tokens[0]);
            cmd.args(tokens.iter().skip(1));
            cmd.arg(program_os);
            return cmd;
        };
        command.args(tokens.iter().skip(1));
        return command;
    }

    std::process::Command::new(program_os)
}

#[derive(Debug, serde::Serialize)]
pub struct TempCleanupResult {
    pub removed_paths: Vec<String>,
    pub failed_paths: Vec<String>,
    pub bytes_freed: u64,
}

fn estimate_path_size_bytes(path: &Path) -> u64 {
    let Ok(meta) = std::fs::symlink_metadata(path) else {
        return 0;
    };

    if meta.is_dir() {
        let mut total = 0u64;
        if let Ok(entries) = std::fs::read_dir(path) {
            for entry in entries.flatten() {
                total = total.saturating_add(estimate_path_size_bytes(&entry.path()));
            }
        }
        total
    } else {
        meta.len()
    }
}

fn temp_roots_for_cleanup() -> Vec<PathBuf> {
    let mut roots: Vec<PathBuf> = Vec::new();
    let mut seen: HashSet<PathBuf> = HashSet::new();
    let candidates = vec![PathBuf::from("/tmp"), std::env::temp_dir()];

    for candidate in candidates {
        if !candidate.exists() {
            continue;
        }
        let normalized = candidate.canonicalize().unwrap_or(candidate);
        if seen.insert(normalized.clone()) {
            roots.push(normalized);
        }
    }

    roots
}

/// Clears QCortado temporary working directories from system temp roots.
#[tauri::command]
fn clear_temp_storage() -> Result<TempCleanupResult, String> {
    let mut removed_paths: Vec<String> = Vec::new();
    let mut failed_paths: Vec<String> = Vec::new();
    let mut bytes_freed: u64 = 0;

    for root in temp_roots_for_cleanup() {
        let entries = match std::fs::read_dir(&root) {
            Ok(entries) => entries,
            Err(_) => continue,
        };

        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            if !name.to_ascii_lowercase().starts_with("qcortado") {
                continue;
            }

            let path = entry.path();
            let size_bytes = estimate_path_size_bytes(&path);
            let path_text = path.display().to_string();

            let remove_result = match entry.file_type() {
                Ok(file_type) if file_type.is_dir() => std::fs::remove_dir_all(&path),
                Ok(_) => std::fs::remove_file(&path),
                Err(err) => Err(err),
            };

            match remove_result {
                Ok(_) => {
                    bytes_freed = bytes_freed.saturating_add(size_bytes);
                    removed_paths.push(path_text);
                }
                Err(err) => {
                    failed_paths.push(format!("{} ({})", path_text, err));
                }
            }
        }
    }

    Ok(TempCleanupResult {
        removed_paths,
        failed_paths,
        bytes_freed,
    })
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

/// Sets a command prefix prepended to all QE launches (e.g. "mpirun").
#[tauri::command]
fn set_execution_prefix(
    app: AppHandle,
    prefix: Option<String>,
    state: State<AppState>,
) -> Result<(), String> {
    let normalized = normalize_execution_prefix(prefix);
    *state.execution_prefix.lock().unwrap() = normalized.clone();
    config::update_execution_prefix(&app, normalized)
}

/// Gets the currently configured command prefix used for QE launches.
#[tauri::command]
fn get_execution_prefix(state: State<AppState>) -> Option<String> {
    state.execution_prefix.lock().unwrap().clone()
}

/// Checks if QE is configured and which executables are available.
#[tauri::command]
fn check_qe_executables(state: State<AppState>) -> Result<Vec<String>, String> {
    let guard = state.qe_bin_dir.lock().unwrap();
    let bin_dir = guard.as_ref().ok_or("QE path not configured")?;

    let executables = [
        "pw.x",
        "bands.x",
        "dos.x",
        "projwfc.x",
        "pp.x",
        "ph.x",
        "dynmat.x",
        "plotband.x",
        "neb.x",
        "hp.x",
        "turbo_lanczos.x",
        "xspectra.x",
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

/// Gets the number of available CPU cores.
#[tauri::command]
fn get_cpu_count() -> usize {
    std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
}

/// Checks if mpirun is available on the system.
#[tauri::command]
fn check_mpi_available(state: State<AppState>) -> bool {
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();
    std_command_with_prefix("mpirun", execution_prefix.as_deref())
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Runs a pw.x calculation (blocking).
#[tauri::command]
async fn run_calculation(
    calculation: QECalculation,
    working_dir: String,
    state: State<'_, AppState>,
) -> Result<QEResult, String> {
    // Clone state needed before await points
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let runner = QERunner::new(bin_dir).with_execution_prefix(execution_prefix);
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

/// MPI configuration for parallel calculations
#[derive(Debug, Clone, serde::Deserialize)]
pub struct MpiConfig {
    /// Whether MPI is enabled
    pub enabled: bool,
    /// Number of MPI processes
    pub nprocs: u32,
}

/// Runs a pw.x calculation with streaming output via events.
#[tauri::command]
async fn run_calculation_streaming(
    app: AppHandle,
    calculation: QECalculation,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    state: State<'_, AppState>,
) -> Result<QEResult, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    // Clone state out of the guard before any await points
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let input = generate_pw_input(&calculation);
    let work_path = PathBuf::from(&working_dir);

    // Ensure working directory exists
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    let exe_path = bin_dir.join("pw.x");
    if !exe_path.exists() {
        return Err("pw.x not found".to_string());
    }

    // Build the command - with or without MPI
    let mut child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            // Use MPI
            let _ = app.emit(
                "qe-output-line",
                format!("Starting pw.x with MPI ({} processes)...", mpi.nprocs),
            );
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}. Is MPI installed?", e))?
        } else {
            // Serial mode
            tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        // No MPI config provided - serial mode
        tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    // Write input to stdin
    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write input: {}", e))?;
    }

    // Stream stdout line by line
    let stdout = child.stdout.take().ok_or("Failed to capture stdout")?;
    let mut reader = BufReader::new(stdout).lines();
    let mut full_output = String::new();

    while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
        full_output.push_str(&line);
        full_output.push('\n');
        // Emit event to frontend
        let _ = app.emit("qe-output-line", &line);
    }

    // Wait for process to complete
    let status = child.wait().await.map_err(|e| e.to_string())?;

    if !status.success() {
        // Emit the error
        let _ = app.emit(
            "qe-output-line",
            format!("\nProcess exited with code: {:?}", status.code()),
        );
        return Err(format!("pw.x failed with exit code: {:?}", status.code()));
    }

    // Parse and return the result
    Ok(parse_pw_output(&full_output))
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
fn load_sssp_data(
    pseudo_dir: String,
) -> Result<std::collections::HashMap<String, SSSPElementData>, String> {
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

    let data: std::collections::HashMap<String, SSSPElementData> =
        serde_json::from_str(&content).map_err(|e| format!("Failed to parse SSSP JSON: {}", e))?;

    Ok(data)
}

// ============================================================================
// Band Structure Commands
// ============================================================================

/// Configuration for a bands NSCF calculation
#[derive(Debug, Clone, serde::Deserialize)]
pub struct BandsCalculationConfig {
    /// Base SCF calculation to derive system settings from
    pub base_calculation: QECalculation,
    /// K-path for band structure
    pub k_path: Vec<KPathPoint>,
    /// Optional number of bands (auto if None)
    pub nbnd: Option<u32>,
    /// Project ID containing the source SCF calculation
    pub project_id: Option<String>,
    /// SCF calculation ID to get the .save directory from
    pub scf_calc_id: Option<String>,
    /// Optional projection analysis settings for fat-band rendering
    pub projections: Option<BandsProjectionOptions>,
}

#[derive(Debug, Clone, serde::Deserialize)]
pub struct BandsProjectionOptions {
    /// Run projwfc.x after bands.x and attach projection groups to the returned band data
    pub enabled: bool,
    /// projwfc.x lsym option (symmetrize projections)
    pub lsym: Option<bool>,
    /// projwfc.x diag_basis option
    pub diag_basis: Option<bool>,
    /// projwfc.x pawproj option
    pub pawproj: Option<bool>,
}

/// Runs a band structure calculation (NSCF + bands.x) with streaming output.
#[tauri::command]
async fn run_bands_calculation(
    app: AppHandle,
    config: BandsCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    state: State<'_, AppState>,
) -> Result<BandData, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    // Clone state out of the guard before any await points
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let work_path = PathBuf::from(&working_dir);

    // Ensure working directory exists
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    // Copy SCF's .save directory if project/calculation IDs are provided
    if let (Some(ref project_id), Some(ref scf_calc_id)) = (&config.project_id, &config.scf_calc_id)
    {
        let projects_dir = projects::get_projects_dir(&app)?;
        let scf_tmp_dir = projects_dir
            .join(project_id)
            .join("calculations")
            .join(scf_calc_id)
            .join("tmp");

        if scf_tmp_dir.exists() {
            let _ = app.emit(
                "qe-output-line",
                format!("SCF tmp dir: {}", scf_tmp_dir.display()),
            );

            // Check for .save directory
            let save_dir = scf_tmp_dir.join("qcortado_scf.save");
            if save_dir.exists() {
                let _ = app.emit(
                    "qe-output-line",
                    format!("Found .save directory: {}", save_dir.display()),
                );
            } else {
                let _ = app.emit("qe-output-line", "WARNING: .save directory not found!");
                // List contents of tmp dir
                if let Ok(entries) = std::fs::read_dir(&scf_tmp_dir) {
                    for entry in entries.flatten() {
                        let _ = app.emit(
                            "qe-output-line",
                            format!("  - {}", entry.file_name().to_string_lossy()),
                        );
                    }
                }
            }

            let _ = app.emit("qe-output-line", "Copying SCF data to working directory...");

            // Copy everything from the SCF tmp dir (includes .save directory)
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;

            let _ = app.emit("qe-output-line", "SCF data copied successfully.");

            // Verify copy
            let copied_save = work_path.join("qcortado_scf.save");
            if copied_save.exists() {
                let _ = app.emit(
                    "qe-output-line",
                    format!("Verified .save in working dir: {}", copied_save.display()),
                );
            } else {
                let _ = app.emit(
                    "qe-output-line",
                    "WARNING: .save not found in working dir after copy!",
                );
            }
        } else {
            return Err(format!(
                "SCF calculation tmp directory not found: {}",
                scf_tmp_dir.display()
            ));
        }
    }

    // Step 1: Create bands calculation from base SCF
    let mut bands_calc = config.base_calculation.clone();
    bands_calc.calculation = qe::CalculationType::Bands;
    bands_calc.verbosity = Some("high".to_string());

    // Convert k_path to KPoints::CrystalB
    let band_path: Vec<qe::BandPathPoint> = config
        .k_path
        .iter()
        .map(|p| qe::BandPathPoint {
            k: p.coords,
            npoints: p.npoints,
            label: Some(p.label.clone()),
        })
        .collect();
    bands_calc.kpoints = qe::KPoints::CrystalB { path: band_path };

    // Generate input
    let input = generate_pw_input(&bands_calc);

    // Save input file for reference
    std::fs::write(work_path.join("bands.in"), &input)
        .map_err(|e| format!("Failed to write input file: {}", e))?;

    // Log the K_POINTS section of the input for debugging
    let _ = app.emit("qe-output-line", "");
    let _ = app.emit("qe-output-line", "=== Generated K_POINTS section ===");
    for line in input.lines() {
        if line.contains("K_POINTS")
            || line.trim().starts_with("0.")
            || line.trim().starts_with("-0.")
            || line.trim().parse::<i32>().is_ok()
        {
            let _ = app.emit("qe-output-line", line);
        }
    }
    let _ = app.emit("qe-output-line", "=== End K_POINTS ===");
    let _ = app.emit("qe-output-line", "");

    let _ = app.emit(
        "qe-output-line",
        "=== Starting Band Structure Calculation ===",
    );

    // Log the k-path being used
    let _ = app.emit(
        "qe-output-line",
        format!("K-path has {} points:", config.k_path.len()),
    );
    for (i, point) in config.k_path.iter().enumerate() {
        let _ = app.emit(
            "qe-output-line",
            format!(
                "  {}: {} ({:.4}, {:.4}, {:.4}) -> {} points to next",
                i + 1,
                point.label,
                point.coords[0],
                point.coords[1],
                point.coords[2],
                point.npoints
            ),
        );
    }

    let projections_enabled = config
        .projections
        .as_ref()
        .map(|p| p.enabled)
        .unwrap_or(false);
    let total_steps = if projections_enabled { 3 } else { 2 };

    let _ = app.emit(
        "qe-output-line",
        format!(
            "Step 1/{}: Running NSCF calculation along k-path...",
            total_steps
        ),
    );

    // Step 2: Run pw.x for bands
    let exe_path = bin_dir.join("pw.x");
    if !exe_path.exists() {
        return Err("pw.x not found".to_string());
    }

    // Build the command - with or without MPI
    let mut child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            let _ = app.emit(
                "qe-output-line",
                format!("Using MPI with {} processes", mpi.nprocs),
            );
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    // Write input to stdin
    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write input: {}", e))?;
    }

    // Stream stdout
    let stdout = child.stdout.take().ok_or("Failed to capture stdout")?;
    let mut reader = BufReader::new(stdout).lines();
    let mut full_output = String::new();
    let mut fermi_energy: Option<f64> = None;

    while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
        full_output.push_str(&line);
        full_output.push('\n');
        let _ = app.emit("qe-output-line", &line);

        // Try to extract Fermi energy from output
        if line.contains("the Fermi energy is") {
            if let Some(idx) = line.find("the Fermi energy is") {
                let rest = &line[idx + 19..];
                if let Some(ev_idx) = rest.find("ev") {
                    if let Ok(ef) = rest[..ev_idx].trim().parse::<f64>() {
                        fermi_energy = Some(ef);
                    }
                }
            }
        }
    }

    // Wait for process to complete
    let status = child.wait().await.map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!("pw.x failed with exit code: {:?}", status.code()));
    }

    // Save output
    std::fs::write(work_path.join("bands.out"), &full_output)
        .map_err(|e| format!("Failed to write output file: {}", e))?;

    let _ = app.emit("qe-output-line", "");
    let _ = app.emit(
        "qe-output-line",
        format!("Step 2/{}: Running bands.x post-processing...", total_steps),
    );

    // Step 3: Run bands.x
    let bands_x_path = bin_dir.join("bands.x");
    if !bands_x_path.exists() {
        return Err(
            "bands.x not found. Make sure your QE installation includes bands.x".to_string(),
        );
    }

    let bands_x_config = BandsXConfig {
        prefix: bands_calc.prefix.clone(),
        outdir: bands_calc.outdir.clone(),
        filband: "bands.dat".to_string(),
        lsym: true,
    };
    let bands_x_input = generate_bands_x_input(&bands_x_config);

    // Save bands.x input for reference
    std::fs::write(work_path.join("bands_pp.in"), &bands_x_input)
        .map_err(|e| format!("Failed to write bands.x input: {}", e))?;

    let mut bands_child = tokio_command_with_prefix(&bands_x_path, execution_prefix.as_deref())
        .current_dir(&work_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start bands.x: {}", e))?;

    if let Some(mut stdin) = bands_child.stdin.take() {
        stdin
            .write_all(bands_x_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write bands.x input: {}", e))?;
    }

    // Stream bands.x output
    let bands_stdout = bands_child
        .stdout
        .take()
        .ok_or("Failed to capture bands.x stdout")?;
    let mut bands_reader = BufReader::new(bands_stdout).lines();
    let mut bands_output = String::new();

    while let Some(line) = bands_reader.next_line().await.map_err(|e| e.to_string())? {
        bands_output.push_str(&line);
        bands_output.push('\n');
        let _ = app.emit("qe-output-line", &line);
    }

    let bands_status = bands_child.wait().await.map_err(|e| e.to_string())?;
    if !bands_status.success() {
        return Err(format!(
            "bands.x failed with exit code: {:?}",
            bands_status.code()
        ));
    }

    std::fs::write(work_path.join("bands_pp.out"), &bands_output)
        .map_err(|e| format!("Failed to write bands.x output: {}", e))?;

    let _ = app.emit("qe-output-line", "");
    let _ = app.emit("qe-output-line", "Parsing band structure data...");

    // Step 4: Parse the output
    let gnu_file = work_path.join("bands.dat.gnu");
    if !gnu_file.exists() {
        return Err("bands.dat.gnu not found. bands.x may have failed.".to_string());
    }

    // Log file size for debugging
    if let Ok(metadata) = std::fs::metadata(&gnu_file) {
        let _ = app.emit(
            "qe-output-line",
            format!("bands.dat.gnu size: {} bytes", metadata.len()),
        );
    }

    let ef = fermi_energy.unwrap_or(0.0);
    let _ = app.emit(
        "qe-output-line",
        format!("Using Fermi energy: {:.4} eV", ef),
    );

    let mut band_data = read_bands_gnu_file(&gnu_file, ef)
        .map_err(|e| format!("Failed to parse band data: {}", e))?;

    // Log some stats about parsed data
    let _ = app.emit(
        "qe-output-line",
        format!(
            "Parsed: {} bands, {} k-points, energy range [{:.2}, {:.2}] eV",
            band_data.n_bands,
            band_data.n_kpoints,
            band_data.energy_range[0],
            band_data.energy_range[1]
        ),
    );

    // Add high-symmetry point markers
    qe::bands::add_symmetry_markers(&mut band_data, &config.k_path);

    if projections_enabled {
        let _ = app.emit("qe-output-line", "");
        let _ = app.emit(
            "qe-output-line",
            format!(
                "Step 3/{}: Running projwfc.x orbital projections...",
                total_steps
            ),
        );

        let projwfc_x_path = bin_dir.join("projwfc.x");
        if !projwfc_x_path.exists() {
            let _ = app.emit(
                "qe-output-line",
                "WARNING: projwfc.x not found. Skipping fat-band projection analysis.",
            );
        } else {
            let projection_options = config.projections.as_ref();
            let projwfc_config = ProjwfcConfig {
                prefix: bands_calc.prefix.clone(),
                outdir: bands_calc.outdir.clone(),
                filproj: "bands.projwfc.dat".to_string(),
                lsym: projection_options
                    .and_then(|opts| opts.lsym)
                    .unwrap_or(false),
                diag_basis: projection_options
                    .and_then(|opts| opts.diag_basis)
                    .unwrap_or(false),
                pawproj: projection_options
                    .and_then(|opts| opts.pawproj)
                    .unwrap_or(false),
            };
            let projwfc_input = generate_projwfc_input(&projwfc_config);

            std::fs::write(work_path.join("projwfc.in"), &projwfc_input)
                .map_err(|e| format!("Failed to write projwfc.x input: {}", e))?;

            let mut projwfc_child =
                tokio_command_with_prefix(&projwfc_x_path, execution_prefix.as_deref())
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start projwfc.x: {}", e))?;

            if let Some(mut stdin) = projwfc_child.stdin.take() {
                stdin
                    .write_all(projwfc_input.as_bytes())
                    .await
                    .map_err(|e| format!("Failed to write projwfc.x input: {}", e))?;
            }

            let projwfc_stdout = projwfc_child
                .stdout
                .take()
                .ok_or("Failed to capture projwfc.x stdout")?;
            let mut projwfc_reader = BufReader::new(projwfc_stdout).lines();
            let mut projwfc_output = String::new();

            while let Some(line) = projwfc_reader
                .next_line()
                .await
                .map_err(|e| e.to_string())?
            {
                projwfc_output.push_str(&line);
                projwfc_output.push('\n');
                let _ = app.emit("qe-output-line", &line);
            }

            let projwfc_status = projwfc_child.wait().await.map_err(|e| e.to_string())?;
            std::fs::write(work_path.join("projwfc.out"), &projwfc_output)
                .map_err(|e| format!("Failed to write projwfc.x output: {}", e))?;

            if !projwfc_status.success() {
                let _ = app.emit(
                    "qe-output-line",
                    format!(
                        "WARNING: projwfc.x failed with exit code {:?}. Continuing without projections.",
                        projwfc_status.code()
                    ),
                );
            } else {
                let projection_text = {
                    let filproj_path = work_path.join(&projwfc_config.filproj);
                    if filproj_path.exists() {
                        std::fs::read_to_string(&filproj_path)
                            .unwrap_or_else(|_| projwfc_output.clone())
                    } else {
                        projwfc_output.clone()
                    }
                };

                if projection_text.trim().is_empty() {
                    let _ = app.emit(
                        "qe-output-line",
                        "WARNING: projwfc output was empty. Continuing without projections.",
                    );
                } else {
                    match parse_projwfc_projection_groups(
                        &projection_text,
                        band_data.n_bands,
                        band_data.n_kpoints,
                    ) {
                        Ok(projections) => {
                            let atom_count = projections.atom_groups.len();
                            let orbital_count = projections.orbital_groups.len();
                            band_data.projections = Some(projections);
                            let _ = app.emit(
                                "qe-output-line",
                                format!(
                                    "Projection groups parsed: {} atom groups, {} orbital groups.",
                                    atom_count, orbital_count
                                ),
                            );
                        }
                        Err(parse_error) => {
                            let _ = app.emit(
                                "qe-output-line",
                                format!(
                                    "WARNING: Could not parse projwfc projections ({}). Continuing without projections.",
                                    parse_error
                                ),
                            );
                        }
                    }
                }
            }
        }
    }

    let _ = app.emit("qe-output-line", "=== Band Structure Complete ===");
    let _ = app.emit(
        "qe-output-line",
        format!(
            "  {} bands, {} k-points",
            band_data.n_bands, band_data.n_kpoints
        ),
    );
    if let Some(ref gap) = band_data.band_gap {
        let gap_type = if gap.is_direct { "direct" } else { "indirect" };
        let _ = app.emit(
            "qe-output-line",
            format!("  Band gap: {:.3} eV ({})", gap.value, gap_type),
        );
    }

    Ok(band_data)
}

// ============================================================================
// Phonon Calculation Commands
// ============================================================================

/// Runs a complete phonon calculation pipeline (ph.x → q2r.x → matdyn.x) with streaming output.
#[tauri::command]
async fn run_phonon_calculation(
    app: AppHandle,
    config: PhononPipelineConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    state: State<'_, AppState>,
) -> Result<PhononResult, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    // Get QE runtime settings
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let work_path = PathBuf::from(&working_dir);

    // Ensure working directory exists
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    let mut full_output = String::new();

    // Copy SCF's .save directory if project/calculation IDs are provided
    if let (Some(ref project_id), Some(ref scf_calc_id)) = (&config.project_id, &config.scf_calc_id)
    {
        let projects_dir = projects::get_projects_dir(&app)?;
        let scf_tmp_dir = projects_dir
            .join(project_id)
            .join("calculations")
            .join(scf_calc_id)
            .join("tmp");

        if scf_tmp_dir.exists() {
            let _ = app.emit(
                "qe-output-line",
                format!("Copying SCF data from: {}", scf_tmp_dir.display()),
            );
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;
            let _ = app.emit("qe-output-line", "SCF data copied successfully.");
        } else {
            return Err(format!(
                "SCF calculation tmp directory not found: {}",
                scf_tmp_dir.display()
            ));
        }
    }

    // ========================================================================
    // Step 1: Run ph.x
    // ========================================================================
    let _ = app.emit("qe-output-line", "");
    let _ = app.emit(
        "qe-output-line",
        "=== Step 1/4: Running ph.x Phonon Calculation ===",
    );
    let _ = app.emit(
        "qe-output-line",
        format!(
            "Q-grid: {}×{}×{}",
            config.phonon.nq[0], config.phonon.nq[1], config.phonon.nq[2]
        ),
    );

    let ph_exe = bin_dir.join("ph.x");
    if !ph_exe.exists() {
        return Err("ph.x not found. Make sure your QE installation includes ph.x".to_string());
    }

    let ph_input = generate_ph_input(&config.phonon);

    // Save ph.x input for reference
    std::fs::write(work_path.join("ph.in"), &ph_input)
        .map_err(|e| format!("Failed to write ph.x input: {}", e))?;

    // Build ph.x command with optional MPI
    let mut ph_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            let _ = app.emit(
                "qe-output-line",
                format!("Using MPI with {} processes", mpi.nprocs),
            );
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&ph_exe)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio_command_with_prefix(&ph_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start ph.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&ph_exe, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start ph.x: {}", e))?
    };

    // Write input to stdin
    if let Some(mut stdin) = ph_child.stdin.take() {
        stdin
            .write_all(ph_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write ph.x input: {}", e))?;
    }

    // Stream ph.x output
    let ph_stdout = ph_child
        .stdout
        .take()
        .ok_or("Failed to capture ph.x stdout")?;
    let mut ph_reader = BufReader::new(ph_stdout).lines();

    while let Some(line) = ph_reader.next_line().await.map_err(|e| e.to_string())? {
        full_output.push_str(&line);
        full_output.push('\n');
        let _ = app.emit("qe-output-line", &line);
    }

    let ph_status = ph_child.wait().await.map_err(|e| e.to_string())?;
    let (converged, n_qpoints) = parse_ph_output(&full_output);
    if !ph_status.success() {
        if converged {
            let _ = app.emit(
                "qe-output-line",
                format!(
                    "Warning: ph.x exited with code {:?} but output contains JOB DONE. Continuing with recovery mode.",
                    ph_status.code()
                ),
            );
        } else {
            return Err(format!(
                "ph.x failed with exit code: {:?}",
                ph_status.code()
            ));
        }
    }

    // Check convergence from stdout markers
    if !converged {
        return Err("ph.x did not converge successfully".to_string());
    }

    let _ = app.emit(
        "qe-output-line",
        format!("ph.x completed: {} q-points calculated", n_qpoints),
    );

    // ========================================================================
    // Step 2: Run q2r.x
    // ========================================================================
    let _ = app.emit("qe-output-line", "");
    let _ = app.emit(
        "qe-output-line",
        "=== Step 2/4: Running q2r.x (Force Constants) ===",
    );

    let q2r_exe = bin_dir.join("q2r.x");
    if !q2r_exe.exists() {
        return Err("q2r.x not found".to_string());
    }

    let q2r_calc = Q2RCalculation {
        fildyn: config.phonon.fildyn.clone(),
        flfrc: "force_constants".to_string(),
        zasr: config.phonon.asr.clone(),
    };
    let q2r_input = generate_q2r_input(&q2r_calc);

    std::fs::write(work_path.join("q2r.in"), &q2r_input)
        .map_err(|e| format!("Failed to write q2r.x input: {}", e))?;

    let mut q2r_child = tokio_command_with_prefix(&q2r_exe, execution_prefix.as_deref())
        .current_dir(&work_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start q2r.x: {}", e))?;

    if let Some(mut stdin) = q2r_child.stdin.take() {
        stdin
            .write_all(q2r_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write q2r.x input: {}", e))?;
    }

    let q2r_stdout = q2r_child
        .stdout
        .take()
        .ok_or("Failed to capture q2r.x stdout")?;
    let mut q2r_reader = BufReader::new(q2r_stdout).lines();

    while let Some(line) = q2r_reader.next_line().await.map_err(|e| e.to_string())? {
        let _ = app.emit("qe-output-line", &line);
    }

    let q2r_status = q2r_child.wait().await.map_err(|e| e.to_string())?;
    if !q2r_status.success() {
        return Err(format!(
            "q2r.x failed with exit code: {:?}",
            q2r_status.code()
        ));
    }

    let _ = app.emit("qe-output-line", "q2r.x completed successfully");

    // Variables to hold results
    let mut dos_data = None;
    let mut dispersion_data = None;

    // ========================================================================
    // Step 3: Run matdyn.x for DOS (if requested)
    // ========================================================================
    if config.calculate_dos {
        let _ = app.emit("qe-output-line", "");
        let _ = app.emit(
            "qe-output-line",
            "=== Step 3/4: Running matdyn.x (Phonon DOS) ===",
        );

        let matdyn_exe = bin_dir.join("matdyn.x");
        if !matdyn_exe.exists() {
            return Err("matdyn.x not found".to_string());
        }

        let dos_grid = config.dos_grid.unwrap_or([20, 20, 20]);
        let dos_delta_e = config.dos_delta_e.unwrap_or(1.0);
        let _ = app.emit(
            "qe-output-line",
            format!("DOS grid: {}×{}×{}", dos_grid[0], dos_grid[1], dos_grid[2]),
        );
        let _ = app.emit(
            "qe-output-line",
            format!("DOS deltaE: {:.4} cm^-1", dos_delta_e),
        );

        let matdyn_dos_calc = MatdynCalculation {
            flfrc: "force_constants".to_string(),
            asr: config.phonon.asr.clone(),
            dos: true,
            fldos: Some("phonon_dos".to_string()),
            nk: Some(dos_grid),
            delta_e: Some(dos_delta_e),
            q_path: None,
            flfrq: None,
        };
        let matdyn_dos_input = generate_matdyn_dos_input(&matdyn_dos_calc);

        std::fs::write(work_path.join("matdyn_dos.in"), &matdyn_dos_input)
            .map_err(|e| format!("Failed to write matdyn.x DOS input: {}", e))?;

        let mut matdyn_dos_child =
            tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start matdyn.x for DOS: {}", e))?;

        if let Some(mut stdin) = matdyn_dos_child.stdin.take() {
            stdin
                .write_all(matdyn_dos_input.as_bytes())
                .await
                .map_err(|e| format!("Failed to write matdyn.x DOS input: {}", e))?;
        }

        let matdyn_dos_stdout = matdyn_dos_child
            .stdout
            .take()
            .ok_or("Failed to capture matdyn.x stdout")?;
        let mut matdyn_dos_reader = BufReader::new(matdyn_dos_stdout).lines();

        while let Some(line) = matdyn_dos_reader
            .next_line()
            .await
            .map_err(|e| e.to_string())?
        {
            let _ = app.emit("qe-output-line", &line);
        }

        let matdyn_dos_status = matdyn_dos_child.wait().await.map_err(|e| e.to_string())?;
        if !matdyn_dos_status.success() {
            let _ = app.emit("qe-output-line", "Warning: matdyn.x DOS calculation failed");
        } else {
            // Parse DOS output
            let dos_file = work_path.join("phonon_dos");
            if dos_file.exists() {
                match read_phonon_dos_file(&dos_file) {
                    Ok(dos) => {
                        let _ = app.emit(
                            "qe-output-line",
                            format!(
                                "Phonon DOS: {} points, frequency range [{:.1}, {:.1}] cm^-1",
                                dos.frequencies.len(),
                                dos.omega_min,
                                dos.omega_max
                            ),
                        );
                        dos_data = Some(dos);
                    }
                    Err(e) => {
                        let _ = app.emit(
                            "qe-output-line",
                            format!("Warning: Failed to parse phonon DOS: {}", e),
                        );
                    }
                }
            }
        }
    } else {
        let _ = app.emit("qe-output-line", "");
        let _ = app.emit(
            "qe-output-line",
            "=== Step 3/4: Skipping DOS calculation ===",
        );
    }

    // ========================================================================
    // Step 4: Run matdyn.x for dispersion (if requested)
    // ========================================================================
    if config.calculate_dispersion {
        let _ = app.emit("qe-output-line", "");
        let _ = app.emit(
            "qe-output-line",
            "=== Step 4/4: Running matdyn.x (Phonon Dispersion) ===",
        );

        let matdyn_exe = bin_dir.join("matdyn.x");
        if !matdyn_exe.exists() {
            return Err("matdyn.x not found".to_string());
        }

        if let Some(ref q_path) = config.q_path {
            let _ = app.emit("qe-output-line", format!("Q-path: {} points", q_path.len()));

            // Create q_path with correct npoints
            let q_path_with_points: Vec<QPathPoint> = q_path
                .iter()
                .enumerate()
                .map(|(i, p)| QPathPoint {
                    label: p.label.clone(),
                    coords: p.coords,
                    npoints: if i < q_path.len() - 1 {
                        config.points_per_segment
                    } else {
                        0
                    },
                })
                .collect();

            let matdyn_bands_calc = MatdynCalculation {
                flfrc: "force_constants".to_string(),
                asr: config.phonon.asr.clone(),
                dos: false,
                fldos: None,
                nk: None,
                delta_e: None,
                q_path: Some(q_path_with_points.clone()),
                flfrq: Some("phonon_freq".to_string()),
            };
            let matdyn_bands_input = generate_matdyn_bands_input(&matdyn_bands_calc);

            std::fs::write(work_path.join("matdyn_bands.in"), &matdyn_bands_input)
                .map_err(|e| format!("Failed to write matdyn.x bands input: {}", e))?;

            let mut matdyn_bands_child =
                tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start matdyn.x for dispersion: {}", e))?;

            if let Some(mut stdin) = matdyn_bands_child.stdin.take() {
                stdin
                    .write_all(matdyn_bands_input.as_bytes())
                    .await
                    .map_err(|e| format!("Failed to write matdyn.x bands input: {}", e))?;
            }

            let matdyn_bands_stdout = matdyn_bands_child
                .stdout
                .take()
                .ok_or("Failed to capture matdyn.x stdout")?;
            let mut matdyn_bands_reader = BufReader::new(matdyn_bands_stdout).lines();

            while let Some(line) = matdyn_bands_reader
                .next_line()
                .await
                .map_err(|e| e.to_string())?
            {
                let _ = app.emit("qe-output-line", &line);
            }

            let matdyn_bands_status = matdyn_bands_child.wait().await.map_err(|e| e.to_string())?;
            if !matdyn_bands_status.success() {
                let _ = app.emit(
                    "qe-output-line",
                    "Warning: matdyn.x dispersion calculation failed",
                );
            } else {
                // Parse dispersion output.
                // matdyn writes the gnuplot-friendly data to <flfrq>.gp.
                // Keep fallback to <flfrq> for compatibility.
                let freq_gp_file = work_path.join("phonon_freq.gp");
                let freq_file = work_path.join("phonon_freq");
                let source_file = if freq_gp_file.exists() {
                    Some(freq_gp_file)
                } else if freq_file.exists() {
                    Some(freq_file)
                } else {
                    None
                };

                if let Some(source_file) = source_file {
                    match read_phonon_dispersion_file(&source_file) {
                        Ok(mut disp) => {
                            add_phonon_symmetry_markers(&mut disp, &q_path_with_points);
                            let _ = app.emit("qe-output-line", format!(
                                "Phonon dispersion: {} modes, {} q-points, frequency range [{:.1}, {:.1}] cm^-1",
                                disp.n_modes, disp.n_qpoints, disp.frequency_range[0], disp.frequency_range[1]
                            ));
                            dispersion_data = Some(disp);
                        }
                        Err(e) => {
                            let _ = app.emit(
                                "qe-output-line",
                                format!("Warning: Failed to parse phonon dispersion: {}", e),
                            );
                        }
                    }
                } else {
                    let _ = app.emit(
                        "qe-output-line",
                        "Warning: No phonon dispersion output file found",
                    );
                }
            }
        } else {
            let _ = app.emit("qe-output-line", "No Q-path specified, skipping dispersion");
        }
    } else {
        let _ = app.emit("qe-output-line", "");
        let _ = app.emit(
            "qe-output-line",
            "=== Step 4/4: Skipping dispersion calculation ===",
        );
    }

    // Calculate n_modes from dispersion or DOS data
    let n_modes = dispersion_data.as_ref().map(|d| d.n_modes).unwrap_or(0);

    let _ = app.emit("qe-output-line", "");
    let _ = app.emit("qe-output-line", "=== Phonon Calculation Complete ===");
    let _ = app.emit(
        "qe-output-line",
        format!("  {} q-points, {} modes", n_qpoints, n_modes),
    );
    if dos_data.is_some() {
        let _ = app.emit("qe-output-line", "  DOS: calculated");
    }
    if dispersion_data.is_some() {
        let _ = app.emit("qe-output-line", "  Dispersion: calculated");
    }

    Ok(PhononResult {
        converged: true,
        n_qpoints,
        n_modes,
        dos_data,
        dispersion_data,
        raw_output: full_output,
    })
}

// ============================================================================
// Background Task Commands (Process Manager)
// ============================================================================

/// Starts an SCF calculation as a background task. Returns the task_id immediately.
#[tauri::command]
async fn start_scf_calculation(
    app: AppHandle,
    calculation: QECalculation,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    // Reject if a task is already running
    if state.process_manager.has_running_tasks().await {
        return Err("A calculation is already running. Please wait for it to complete or cancel it.".to_string());
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("scf".to_string(), label).await;

    // We need to drop the Mutex guard from register before spawning
    let tid = task_id.clone();
    let app_handle = app.clone();

    tokio::spawn(async move {
        let result = run_scf_background(
            app_handle.clone(),
            &tid,
            calculation,
            working_dir,
            mpi_config,
            bin_dir,
            execution_prefix,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(qe_result) => {
                let json = serde_json::to_value(&qe_result).unwrap_or(serde_json::Value::Null);
                pm.complete(&tid, json).await;
                let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
            }
            Err(e) => {
                pm.fail(&tid, e.clone()).await;
                let _ = app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
            }
        }
    });

    Ok(task_id)
}

async fn run_scf_background(
    app: AppHandle,
    task_id: &str,
    calculation: QECalculation,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    bin_dir: PathBuf,
    execution_prefix: Option<String>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<QEResult, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let input = generate_pw_input(&calculation);
    let work_path = PathBuf::from(&working_dir);

    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    let exe_path = bin_dir.join("pw.x");
    if !exe_path.exists() {
        return Err("pw.x not found".to_string());
    }

    let mut child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            let line = format!("Starting pw.x with MPI ({} processes)...", mpi.nprocs);
            let _ = app.emit(&format!("task-output:{}", task_id), &line);
            pm.append_output(task_id, line).await;
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}. Is MPI installed?", e))?
        } else {
            tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    // Store child PID for cancellation
    if let Some(pid) = child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write input: {}", e))?;
    }

    let stdout = child.stdout.take().ok_or("Failed to capture stdout")?;
    let mut reader = BufReader::new(stdout).lines();
    let mut full_output = String::new();

    while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
        // Check cancel flag
        if cancel_flag.load(std::sync::atomic::Ordering::SeqCst) {
            return Err("Cancelled by user".to_string());
        }
        full_output.push_str(&line);
        full_output.push('\n');
        let _ = app.emit(&format!("task-output:{}", task_id), &line);
        pm.append_output(task_id, line).await;
    }

    let status = child.wait().await.map_err(|e| e.to_string())?;

    if cancel_flag.load(std::sync::atomic::Ordering::SeqCst) {
        return Err("Cancelled by user".to_string());
    }

    if !status.success() {
        let err_line = format!("\nProcess exited with code: {:?}", status.code());
        let _ = app.emit(&format!("task-output:{}", task_id), &err_line);
        pm.append_output(task_id, err_line).await;
        return Err(format!("pw.x failed with exit code: {:?}", status.code()));
    }

    Ok(parse_pw_output(&full_output))
}

/// Starts a band structure calculation as a background task.
#[tauri::command]
async fn start_bands_calculation(
    app: AppHandle,
    config: BandsCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    if state.process_manager.has_running_tasks().await {
        return Err("A calculation is already running. Please wait for it to complete or cancel it.".to_string());
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("bands".to_string(), label).await;

    let tid = task_id.clone();
    let app_handle = app.clone();

    tokio::spawn(async move {
        let result = run_bands_background(
            app_handle.clone(),
            &tid,
            config,
            working_dir,
            mpi_config,
            bin_dir,
            execution_prefix,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(band_data) => {
                let json = serde_json::to_value(&band_data).unwrap_or(serde_json::Value::Null);
                pm.complete(&tid, json).await;
                let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
            }
            Err(e) => {
                pm.fail(&tid, e.clone()).await;
                let _ = app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
            }
        }
    });

    Ok(task_id)
}

#[allow(clippy::too_many_arguments)]
async fn run_bands_background(
    app: AppHandle,
    task_id: &str,
    config: BandsCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    bin_dir: PathBuf,
    execution_prefix: Option<String>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<BandData, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let work_path = PathBuf::from(&working_dir);
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    // Helper to emit output to both the task event and the output buffer
    macro_rules! emit_line {
        ($line:expr) => {{
            let line_str: String = $line.into();
            let _ = app.emit(&format!("task-output:{}", task_id), &line_str);
            pm.append_output(task_id, line_str).await;
        }};
    }

    macro_rules! check_cancel {
        () => {
            if cancel_flag.load(std::sync::atomic::Ordering::SeqCst) {
                return Err("Cancelled by user".to_string());
            }
        };
    }

    // Copy SCF .save directory if provided
    if let (Some(ref project_id), Some(ref scf_calc_id)) = (&config.project_id, &config.scf_calc_id)
    {
        let projects_dir = projects::get_projects_dir(&app)?;
        let scf_tmp_dir = projects_dir
            .join(project_id)
            .join("calculations")
            .join(scf_calc_id)
            .join("tmp");

        if scf_tmp_dir.exists() {
            emit_line!(format!("SCF tmp dir: {}", scf_tmp_dir.display()));
            let save_dir = scf_tmp_dir.join("qcortado_scf.save");
            if save_dir.exists() {
                emit_line!(format!("Found .save directory: {}", save_dir.display()));
            } else {
                emit_line!("WARNING: .save directory not found!".to_string());
            }
            emit_line!("Copying SCF data to working directory...".to_string());
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;
            emit_line!("SCF data copied successfully.".to_string());
            let copied_save = work_path.join("qcortado_scf.save");
            if copied_save.exists() {
                emit_line!(format!("Verified .save in working dir: {}", copied_save.display()));
            } else {
                emit_line!("WARNING: .save not found in working dir after copy!".to_string());
            }
        } else {
            return Err(format!(
                "SCF calculation tmp directory not found: {}",
                scf_tmp_dir.display()
            ));
        }
    }

    check_cancel!();

    // Step 1: NSCF along k-path
    let mut bands_calc = config.base_calculation.clone();
    bands_calc.calculation = qe::CalculationType::Bands;
    bands_calc.verbosity = Some("high".to_string());

    let band_path: Vec<qe::BandPathPoint> = config
        .k_path
        .iter()
        .map(|p| qe::BandPathPoint {
            k: p.coords,
            npoints: p.npoints,
            label: Some(p.label.clone()),
        })
        .collect();
    bands_calc.kpoints = qe::KPoints::CrystalB { path: band_path };

    let input = generate_pw_input(&bands_calc);
    std::fs::write(work_path.join("bands.in"), &input)
        .map_err(|e| format!("Failed to write input file: {}", e))?;

    emit_line!("".to_string());
    emit_line!("=== Generated K_POINTS section ===".to_string());
    for line in input.lines() {
        if line.contains("K_POINTS")
            || line.trim().starts_with("0.")
            || line.trim().starts_with("-0.")
            || line.trim().parse::<i32>().is_ok()
        {
            emit_line!(line.to_string());
        }
    }
    emit_line!("=== End K_POINTS ===".to_string());
    emit_line!("".to_string());
    emit_line!("=== Starting Band Structure Calculation ===".to_string());

    emit_line!(format!("K-path has {} points:", config.k_path.len()));
    for (i, point) in config.k_path.iter().enumerate() {
        emit_line!(format!(
            "  {}: {} ({:.4}, {:.4}, {:.4}) -> {} points to next",
            i + 1, point.label, point.coords[0], point.coords[1], point.coords[2], point.npoints
        ));
    }

    let projections_enabled = config.projections.as_ref().map(|p| p.enabled).unwrap_or(false);
    let total_steps = if projections_enabled { 3 } else { 2 };
    emit_line!(format!("Step 1/{}: Running NSCF calculation along k-path...", total_steps));

    let exe_path = bin_dir.join("pw.x");
    if !exe_path.exists() {
        return Err("pw.x not found".to_string());
    }

    let mut child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!("Using MPI with {} processes", mpi.nprocs));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&exe_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    if let Some(pid) = child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(input.as_bytes()).await.map_err(|e| format!("Failed to write input: {}", e))?;
    }

    let stdout = child.stdout.take().ok_or("Failed to capture stdout")?;
    let mut reader = BufReader::new(stdout).lines();
    let mut full_output = String::new();
    let mut fermi_energy: Option<f64> = None;

    while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        full_output.push_str(&line);
        full_output.push('\n');
        emit_line!(line.clone());

        if line.contains("the Fermi energy is") {
            if let Some(idx) = line.find("the Fermi energy is") {
                let rest = &line[idx + 19..];
                if let Some(ev_idx) = rest.find("ev") {
                    if let Ok(ef) = rest[..ev_idx].trim().parse::<f64>() {
                        fermi_energy = Some(ef);
                    }
                }
            }
        }
    }

    let status = child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !status.success() {
        return Err(format!("pw.x failed with exit code: {:?}", status.code()));
    }

    std::fs::write(work_path.join("bands.out"), &full_output)
        .map_err(|e| format!("Failed to write output file: {}", e))?;

    emit_line!("".to_string());
    emit_line!(format!("Step 2/{}: Running bands.x post-processing...", total_steps));

    // Step 2: bands.x
    let bands_x_path = bin_dir.join("bands.x");
    if !bands_x_path.exists() {
        return Err("bands.x not found. Make sure your QE installation includes bands.x".to_string());
    }

    let bands_x_config = BandsXConfig {
        prefix: bands_calc.prefix.clone(),
        outdir: bands_calc.outdir.clone(),
        filband: "bands.dat".to_string(),
        lsym: true,
    };
    let bands_x_input = generate_bands_x_input(&bands_x_config);
    std::fs::write(work_path.join("bands_pp.in"), &bands_x_input)
        .map_err(|e| format!("Failed to write bands.x input: {}", e))?;

    let mut bands_child = tokio_command_with_prefix(&bands_x_path, execution_prefix.as_deref())
        .current_dir(&work_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start bands.x: {}", e))?;

    if let Some(pid) = bands_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = bands_child.stdin.take() {
        stdin.write_all(bands_x_input.as_bytes()).await.map_err(|e| format!("Failed to write bands.x input: {}", e))?;
    }

    let bands_stdout = bands_child.stdout.take().ok_or("Failed to capture bands.x stdout")?;
    let mut bands_reader = BufReader::new(bands_stdout).lines();

    while let Some(line) = bands_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        emit_line!(line);
    }

    let bands_status = bands_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !bands_status.success() {
        return Err(format!("bands.x failed with exit code: {:?}", bands_status.code()));
    }

    emit_line!("".to_string());
    emit_line!("Parsing band structure data...".to_string());

    let gnu_file = work_path.join("bands.dat.gnu");
    if !gnu_file.exists() {
        return Err("bands.dat.gnu not found. bands.x may have failed.".to_string());
    }

    if let Ok(metadata) = std::fs::metadata(&gnu_file) {
        emit_line!(format!("bands.dat.gnu size: {} bytes", metadata.len()));
    }

    let ef = fermi_energy.unwrap_or(0.0);
    emit_line!(format!("Using Fermi energy: {:.4} eV", ef));

    let mut band_data = read_bands_gnu_file(&gnu_file, ef)
        .map_err(|e| format!("Failed to parse band data: {}", e))?;

    emit_line!(format!(
        "Parsed: {} bands, {} k-points, energy range [{:.2}, {:.2}] eV",
        band_data.n_bands, band_data.n_kpoints, band_data.energy_range[0], band_data.energy_range[1]
    ));

    qe::bands::add_symmetry_markers(&mut band_data, &config.k_path);

    // Step 3: projwfc.x (optional)
    if projections_enabled {
        emit_line!("".to_string());
        emit_line!(format!("Step 3/{}: Running projwfc.x orbital projections...", total_steps));

        let projwfc_x_path = bin_dir.join("projwfc.x");
        if !projwfc_x_path.exists() {
            emit_line!("WARNING: projwfc.x not found. Skipping fat-band projection analysis.".to_string());
        } else {
            check_cancel!();
            let projection_options = config.projections.as_ref();
            let projwfc_config = ProjwfcConfig {
                prefix: bands_calc.prefix.clone(),
                outdir: bands_calc.outdir.clone(),
                filproj: "bands.projwfc.dat".to_string(),
                lsym: projection_options.and_then(|opts| opts.lsym).unwrap_or(false),
                diag_basis: projection_options.and_then(|opts| opts.diag_basis).unwrap_or(false),
                pawproj: projection_options.and_then(|opts| opts.pawproj).unwrap_or(false),
            };
            let projwfc_input = generate_projwfc_input(&projwfc_config);
            std::fs::write(work_path.join("projwfc.in"), &projwfc_input)
                .map_err(|e| format!("Failed to write projwfc.x input: {}", e))?;

            let mut projwfc_child = tokio_command_with_prefix(&projwfc_x_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start projwfc.x: {}", e))?;

            if let Some(pid) = projwfc_child.id() {
                pm.set_child_id(task_id, pid).await;
            }

            if let Some(mut stdin) = projwfc_child.stdin.take() {
                stdin.write_all(projwfc_input.as_bytes()).await.map_err(|e| format!("Failed to write projwfc.x input: {}", e))?;
            }

            let projwfc_stdout = projwfc_child.stdout.take().ok_or("Failed to capture projwfc.x stdout")?;
            let mut projwfc_reader = BufReader::new(projwfc_stdout).lines();
            let mut projwfc_output = String::new();

            while let Some(line) = projwfc_reader.next_line().await.map_err(|e| e.to_string())? {
                check_cancel!();
                projwfc_output.push_str(&line);
                projwfc_output.push('\n');
                emit_line!(line);
            }

            let projwfc_status = projwfc_child.wait().await.map_err(|e| e.to_string())?;
            std::fs::write(work_path.join("projwfc.out"), &projwfc_output)
                .map_err(|e| format!("Failed to write projwfc.x output: {}", e))?;

            if !projwfc_status.success() {
                emit_line!(format!(
                    "WARNING: projwfc.x failed with exit code {:?}. Continuing without projections.",
                    projwfc_status.code()
                ));
            } else {
                let projection_text = {
                    let filproj_path = work_path.join(&projwfc_config.filproj);
                    if filproj_path.exists() {
                        std::fs::read_to_string(&filproj_path).unwrap_or_else(|_| projwfc_output.clone())
                    } else {
                        projwfc_output.clone()
                    }
                };

                if projection_text.trim().is_empty() {
                    emit_line!("WARNING: projwfc output was empty. Continuing without projections.".to_string());
                } else {
                    match parse_projwfc_projection_groups(&projection_text, band_data.n_bands, band_data.n_kpoints) {
                        Ok(projections) => {
                            let atom_count = projections.atom_groups.len();
                            let orbital_count = projections.orbital_groups.len();
                            band_data.projections = Some(projections);
                            emit_line!(format!(
                                "Projection groups parsed: {} atom groups, {} orbital groups.",
                                atom_count, orbital_count
                            ));
                        }
                        Err(parse_error) => {
                            emit_line!(format!(
                                "WARNING: Could not parse projwfc projections ({}). Continuing without projections.",
                                parse_error
                            ));
                        }
                    }
                }
            }
        }
    }

    emit_line!("=== Band Structure Complete ===".to_string());
    emit_line!(format!("  {} bands, {} k-points", band_data.n_bands, band_data.n_kpoints));
    if let Some(ref gap) = band_data.band_gap {
        let gap_type = if gap.is_direct { "direct" } else { "indirect" };
        emit_line!(format!("  Band gap: {:.3} eV ({})", gap.value, gap_type));
    }

    Ok(band_data)
}

/// Starts a phonon calculation as a background task.
#[tauri::command]
async fn start_phonon_calculation(
    app: AppHandle,
    config: PhononPipelineConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    if state.process_manager.has_running_tasks().await {
        return Err("A calculation is already running. Please wait for it to complete or cancel it.".to_string());
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("phonon".to_string(), label).await;

    let tid = task_id.clone();
    let app_handle = app.clone();

    tokio::spawn(async move {
        let result = run_phonon_background(
            app_handle.clone(),
            &tid,
            config,
            working_dir,
            mpi_config,
            bin_dir,
            execution_prefix,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(phonon_result) => {
                let json = serde_json::to_value(&phonon_result).unwrap_or(serde_json::Value::Null);
                pm.complete(&tid, json).await;
                let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
            }
            Err(e) => {
                pm.fail(&tid, e.clone()).await;
                let _ = app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
            }
        }
    });

    Ok(task_id)
}

#[allow(clippy::too_many_arguments)]
async fn run_phonon_background(
    app: AppHandle,
    task_id: &str,
    config: PhononPipelineConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    bin_dir: PathBuf,
    execution_prefix: Option<String>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<PhononResult, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let work_path = PathBuf::from(&working_dir);
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    let mut full_output = String::new();

    macro_rules! emit_line {
        ($line:expr) => {{
            let line_str: String = $line.into();
            let _ = app.emit(&format!("task-output:{}", task_id), &line_str);
            pm.append_output(task_id, line_str).await;
        }};
    }

    macro_rules! check_cancel {
        () => {
            if cancel_flag.load(std::sync::atomic::Ordering::SeqCst) {
                return Err("Cancelled by user".to_string());
            }
        };
    }

    // Copy SCF .save directory
    if let (Some(ref project_id), Some(ref scf_calc_id)) = (&config.project_id, &config.scf_calc_id)
    {
        let projects_dir = projects::get_projects_dir(&app)?;
        let scf_tmp_dir = projects_dir
            .join(project_id)
            .join("calculations")
            .join(scf_calc_id)
            .join("tmp");

        if scf_tmp_dir.exists() {
            emit_line!(format!("Copying SCF data from: {}", scf_tmp_dir.display()));
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;
            emit_line!("SCF data copied successfully.".to_string());
        } else {
            return Err(format!(
                "SCF calculation tmp directory not found: {}",
                scf_tmp_dir.display()
            ));
        }
    }

    check_cancel!();

    // Step 1: ph.x
    emit_line!("".to_string());
    emit_line!("=== Step 1/4: Running ph.x Phonon Calculation ===".to_string());
    emit_line!(format!("Q-grid: {}×{}×{}", config.phonon.nq[0], config.phonon.nq[1], config.phonon.nq[2]));

    let ph_exe = bin_dir.join("ph.x");
    if !ph_exe.exists() {
        return Err("ph.x not found. Make sure your QE installation includes ph.x".to_string());
    }

    let ph_input = generate_ph_input(&config.phonon);
    std::fs::write(work_path.join("ph.in"), &ph_input)
        .map_err(|e| format!("Failed to write ph.x input: {}", e))?;

    let mut ph_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!("Using MPI with {} processes", mpi.nprocs));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&ph_exe)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio_command_with_prefix(&ph_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start ph.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&ph_exe, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start ph.x: {}", e))?
    };

    if let Some(pid) = ph_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = ph_child.stdin.take() {
        stdin.write_all(ph_input.as_bytes()).await.map_err(|e| format!("Failed to write ph.x input: {}", e))?;
    }

    let ph_stdout = ph_child.stdout.take().ok_or("Failed to capture ph.x stdout")?;
    let mut ph_reader = BufReader::new(ph_stdout).lines();

    while let Some(line) = ph_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        full_output.push_str(&line);
        full_output.push('\n');
        emit_line!(line);
    }

    let ph_status = ph_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();

    let (converged, n_qpoints) = parse_ph_output(&full_output);
    if !ph_status.success() {
        if converged {
            emit_line!(format!(
                "Warning: ph.x exited with code {:?} but output contains JOB DONE. Continuing with recovery mode.",
                ph_status.code()
            ));
        } else {
            return Err(format!("ph.x failed with exit code: {:?}", ph_status.code()));
        }
    }
    if !converged {
        return Err("ph.x did not converge successfully".to_string());
    }
    emit_line!(format!("ph.x completed: {} q-points calculated", n_qpoints));

    // Step 2: q2r.x
    emit_line!("".to_string());
    emit_line!("=== Step 2/4: Running q2r.x (Force Constants) ===".to_string());
    check_cancel!();

    let q2r_exe = bin_dir.join("q2r.x");
    if !q2r_exe.exists() {
        return Err("q2r.x not found".to_string());
    }

    let q2r_calc = Q2RCalculation {
        fildyn: config.phonon.fildyn.clone(),
        flfrc: "force_constants".to_string(),
        zasr: config.phonon.asr.clone(),
    };
    let q2r_input = generate_q2r_input(&q2r_calc);
    std::fs::write(work_path.join("q2r.in"), &q2r_input)
        .map_err(|e| format!("Failed to write q2r.x input: {}", e))?;

    let mut q2r_child = tokio_command_with_prefix(&q2r_exe, execution_prefix.as_deref())
        .current_dir(&work_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start q2r.x: {}", e))?;

    if let Some(pid) = q2r_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = q2r_child.stdin.take() {
        stdin.write_all(q2r_input.as_bytes()).await.map_err(|e| format!("Failed to write q2r.x input: {}", e))?;
    }

    let q2r_stdout = q2r_child.stdout.take().ok_or("Failed to capture q2r.x stdout")?;
    let mut q2r_reader = BufReader::new(q2r_stdout).lines();

    while let Some(line) = q2r_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        emit_line!(line);
    }

    let q2r_status = q2r_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !q2r_status.success() {
        return Err(format!("q2r.x failed with exit code: {:?}", q2r_status.code()));
    }
    emit_line!("q2r.x completed successfully".to_string());

    let mut dos_data = None;
    let mut dispersion_data = None;

    // Step 3: matdyn.x for DOS
    if config.calculate_dos {
        emit_line!("".to_string());
        emit_line!("=== Step 3/4: Running matdyn.x (Phonon DOS) ===".to_string());
        check_cancel!();

        let matdyn_exe = bin_dir.join("matdyn.x");
        if !matdyn_exe.exists() {
            return Err("matdyn.x not found".to_string());
        }

        let dos_grid = config.dos_grid.unwrap_or([20, 20, 20]);
        let dos_delta_e = config.dos_delta_e.unwrap_or(1.0);
        emit_line!(format!("DOS grid: {}×{}×{}", dos_grid[0], dos_grid[1], dos_grid[2]));
        emit_line!(format!("DOS deltaE: {:.4} cm^-1", dos_delta_e));

        let matdyn_dos_calc = MatdynCalculation {
            flfrc: "force_constants".to_string(),
            asr: config.phonon.asr.clone(),
            dos: true,
            fldos: Some("phonon_dos".to_string()),
            nk: Some(dos_grid),
            delta_e: Some(dos_delta_e),
            q_path: None,
            flfrq: None,
        };
        let matdyn_dos_input = generate_matdyn_dos_input(&matdyn_dos_calc);
        std::fs::write(work_path.join("matdyn_dos.in"), &matdyn_dos_input)
            .map_err(|e| format!("Failed to write matdyn.x DOS input: {}", e))?;

        let mut matdyn_dos_child = tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start matdyn.x for DOS: {}", e))?;

        if let Some(pid) = matdyn_dos_child.id() {
            pm.set_child_id(task_id, pid).await;
        }

        if let Some(mut stdin) = matdyn_dos_child.stdin.take() {
            stdin.write_all(matdyn_dos_input.as_bytes()).await.map_err(|e| format!("Failed to write matdyn.x DOS input: {}", e))?;
        }

        let matdyn_dos_stdout = matdyn_dos_child.stdout.take().ok_or("Failed to capture matdyn.x stdout")?;
        let mut matdyn_dos_reader = BufReader::new(matdyn_dos_stdout).lines();

        while let Some(line) = matdyn_dos_reader.next_line().await.map_err(|e| e.to_string())? {
            check_cancel!();
            emit_line!(line);
        }

        let matdyn_dos_status = matdyn_dos_child.wait().await.map_err(|e| e.to_string())?;
        check_cancel!();
        if !matdyn_dos_status.success() {
            emit_line!("Warning: matdyn.x DOS calculation failed".to_string());
        } else {
            let dos_file = work_path.join("phonon_dos");
            if dos_file.exists() {
                match read_phonon_dos_file(&dos_file) {
                    Ok(dos) => {
                        emit_line!(format!(
                            "Phonon DOS: {} points, frequency range [{:.1}, {:.1}] cm^-1",
                            dos.frequencies.len(), dos.omega_min, dos.omega_max
                        ));
                        dos_data = Some(dos);
                    }
                    Err(e) => {
                        emit_line!(format!("Warning: Failed to parse phonon DOS: {}", e));
                    }
                }
            }
        }
    } else {
        emit_line!("".to_string());
        emit_line!("=== Step 3/4: Skipping DOS calculation ===".to_string());
    }

    // Step 4: matdyn.x for dispersion
    if config.calculate_dispersion {
        emit_line!("".to_string());
        emit_line!("=== Step 4/4: Running matdyn.x (Phonon Dispersion) ===".to_string());
        check_cancel!();

        let matdyn_exe = bin_dir.join("matdyn.x");
        if !matdyn_exe.exists() {
            return Err("matdyn.x not found".to_string());
        }

        if let Some(ref q_path) = config.q_path {
            emit_line!(format!("Q-path: {} points", q_path.len()));

            let q_path_with_points: Vec<QPathPoint> = q_path
                .iter()
                .enumerate()
                .map(|(i, p)| QPathPoint {
                    label: p.label.clone(),
                    coords: p.coords,
                    npoints: if i < q_path.len() - 1 { config.points_per_segment } else { 0 },
                })
                .collect();

            let matdyn_bands_calc = MatdynCalculation {
                flfrc: "force_constants".to_string(),
                asr: config.phonon.asr.clone(),
                dos: false,
                fldos: None,
                nk: None,
                delta_e: None,
                q_path: Some(q_path_with_points.clone()),
                flfrq: Some("phonon_freq".to_string()),
            };
            let matdyn_bands_input = generate_matdyn_bands_input(&matdyn_bands_calc);
            std::fs::write(work_path.join("matdyn_bands.in"), &matdyn_bands_input)
                .map_err(|e| format!("Failed to write matdyn.x bands input: {}", e))?;

            let mut matdyn_bands_child = tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start matdyn.x for dispersion: {}", e))?;

            if let Some(pid) = matdyn_bands_child.id() {
                pm.set_child_id(task_id, pid).await;
            }

            if let Some(mut stdin) = matdyn_bands_child.stdin.take() {
                stdin.write_all(matdyn_bands_input.as_bytes()).await.map_err(|e| format!("Failed to write matdyn.x bands input: {}", e))?;
            }

            let matdyn_bands_stdout = matdyn_bands_child.stdout.take().ok_or("Failed to capture matdyn.x stdout")?;
            let mut matdyn_bands_reader = BufReader::new(matdyn_bands_stdout).lines();

            while let Some(line) = matdyn_bands_reader.next_line().await.map_err(|e| e.to_string())? {
                check_cancel!();
                emit_line!(line);
            }

            let matdyn_bands_status = matdyn_bands_child.wait().await.map_err(|e| e.to_string())?;
            check_cancel!();
            if !matdyn_bands_status.success() {
                emit_line!("Warning: matdyn.x dispersion calculation failed".to_string());
            } else {
                let freq_gp_file = work_path.join("phonon_freq.gp");
                let freq_file = work_path.join("phonon_freq");
                let source_file = if freq_gp_file.exists() {
                    Some(freq_gp_file)
                } else if freq_file.exists() {
                    Some(freq_file)
                } else {
                    None
                };

                if let Some(source_file) = source_file {
                    match read_phonon_dispersion_file(&source_file) {
                        Ok(mut disp) => {
                            add_phonon_symmetry_markers(&mut disp, &q_path_with_points);
                            emit_line!(format!(
                                "Phonon dispersion: {} modes, {} q-points, frequency range [{:.1}, {:.1}] cm^-1",
                                disp.n_modes, disp.n_qpoints, disp.frequency_range[0], disp.frequency_range[1]
                            ));
                            dispersion_data = Some(disp);
                        }
                        Err(e) => {
                            emit_line!(format!("Warning: Failed to parse phonon dispersion: {}", e));
                        }
                    }
                } else {
                    emit_line!("Warning: No phonon dispersion output file found".to_string());
                }
            }
        } else {
            emit_line!("No Q-path specified, skipping dispersion".to_string());
        }
    } else {
        emit_line!("".to_string());
        emit_line!("=== Step 4/4: Skipping dispersion calculation ===".to_string());
    }

    let n_modes = dispersion_data.as_ref().map(|d| d.n_modes).unwrap_or(0);

    emit_line!("".to_string());
    emit_line!("=== Phonon Calculation Complete ===".to_string());
    emit_line!(format!("  {} q-points, {} modes", n_qpoints, n_modes));
    if dos_data.is_some() {
        emit_line!("  DOS: calculated".to_string());
    }
    if dispersion_data.is_some() {
        emit_line!("  Dispersion: calculated".to_string());
    }

    Ok(PhononResult {
        converged: true,
        n_qpoints,
        n_modes,
        dos_data,
        dispersion_data,
        raw_output: full_output,
    })
}

// ============================================================================
// Task Management Commands
// ============================================================================

#[tauri::command]
async fn list_running_tasks(
    state: State<'_, AppState>,
) -> Result<Vec<process_manager::TaskSummary>, String> {
    Ok(state.process_manager.list_tasks().await)
}

#[tauri::command]
async fn get_task_info(
    task_id: String,
    state: State<'_, AppState>,
) -> Result<process_manager::TaskInfo, String> {
    state
        .process_manager
        .get_task(&task_id)
        .await
        .ok_or_else(|| format!("Task not found: {}", task_id))
}

#[tauri::command]
async fn get_task_output(
    task_id: String,
    state: State<'_, AppState>,
) -> Result<Vec<String>, String> {
    state
        .process_manager
        .get_output(&task_id)
        .await
        .ok_or_else(|| format!("Task not found: {}", task_id))
}

#[tauri::command]
async fn cancel_task(
    task_id: String,
    state: State<'_, AppState>,
) -> Result<(), String> {
    state.process_manager.cancel(&task_id).await
}

#[tauri::command]
async fn dismiss_task(
    task_id: String,
    state: State<'_, AppState>,
) -> Result<(), String> {
    state.process_manager.remove(&task_id).await;
    Ok(())
}

#[tauri::command]
async fn shutdown_and_close(
    app: AppHandle,
    state: State<'_, AppState>,
) -> Result<(), String> {
    state.process_manager.kill_all().await;
    // Give processes a moment to die
    tokio::time::sleep(std::time::Duration::from_millis(200)).await;
    app.exit(0);
    Ok(())
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
            let mut execution_prefix: Option<String> = None;
            match config::load_config(&app.handle()) {
                Ok(cfg) => {
                    if let Some(path) = cfg.qe_bin_dir {
                        let path_buf = PathBuf::from(&path);
                        // Only use if pw.x still exists
                        if path_buf.join("pw.x").exists() {
                            qe_bin_dir = Some(path_buf);
                        }
                    }
                    execution_prefix = normalize_execution_prefix(cfg.execution_prefix);
                }
                Err(e) => {
                    eprintln!("Warning: Failed to load config: {}", e);
                }
            }

            // Initialize AppState with loaded config
            app.manage(AppState {
                qe_bin_dir: Mutex::new(qe_bin_dir),
                execution_prefix: Mutex::new(execution_prefix),
                project_dir: Mutex::new(None),
                process_manager: ProcessManager::new(),
            });

            Ok(())
        })
        .on_window_event(|window, event| {
            if let tauri::WindowEvent::CloseRequested { api, .. } = event {
                let app = window.app_handle().clone();
                let pm = app.state::<AppState>().process_manager.clone();
                let window_clone = window.clone();
                // Check if tasks are running; if so, prevent close and ask user
                api.prevent_close();
                tauri::async_runtime::spawn(async move {
                    if pm.has_running_tasks().await {
                        let _ = window_clone.emit("confirm-close", ());
                    } else {
                        // No running tasks, just close
                        let _ = window_clone.destroy();
                    }
                });
            }
        })
        .invoke_handler(tauri::generate_handler![
            set_qe_path,
            get_qe_path,
            set_execution_prefix,
            get_execution_prefix,
            clear_temp_storage,
            check_qe_executables,
            generate_input,
            validate_calculation,
            parse_output,
            get_cpu_count,
            check_mpi_available,
            run_calculation,
            run_calculation_streaming,
            run_bands_calculation,
            run_phonon_calculation,
            set_project_dir,
            get_project_dir,
            list_pseudopotentials,
            load_sssp_data,
            // Background task commands
            start_scf_calculation,
            start_bands_calculation,
            start_phonon_calculation,
            list_running_tasks,
            get_task_info,
            get_task_output,
            cancel_task,
            dismiss_task,
            shutdown_and_close,
            // Project management commands
            projects::list_projects,
            projects::create_project,
            projects::get_project,
            projects::add_cif_to_project,
            projects::save_calculation,
            projects::delete_project,
            projects::delete_calculation,
            projects::set_calculation_tag,
            projects::set_last_opened_cif,
            projects::get_cif_crystal_data,
            projects::get_cif_content,
            projects::get_saved_phonon_data,
            projects::recover_phonon_calculation,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
