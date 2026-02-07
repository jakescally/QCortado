//! QCortado - A modern UI for Quantum ESPRESSO
//!
//! This is the Tauri backend providing:
//! - QE input generation and validation
//! - Process execution with streaming output
//! - Output parsing and result extraction
//! - Project management

use std::path::PathBuf;
use std::sync::Mutex;
use tauri::{AppHandle, Emitter, Manager, State};

pub mod config;
pub mod projects;
pub mod qe;

use qe::{generate_pw_input, parse_pw_output, QECalculation, QEResult, QERunner, BandData, BandsXConfig, KPathPoint, generate_bands_x_input, read_bands_gnu_file};

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

/// Gets the number of available CPU cores.
#[tauri::command]
fn get_cpu_count() -> usize {
    std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
}

/// Checks if mpirun is available on the system.
#[tauri::command]
fn check_mpi_available() -> bool {
    std::process::Command::new("mpirun")
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
    use tokio::io::{AsyncBufReadExt, BufReader, AsyncWriteExt};

    // Clone the path out of the guard before any await points
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };

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
            let _ = app.emit("qe-output-line", format!("Starting pw.x with MPI ({} processes)...", mpi.nprocs));
            tokio::process::Command::new("mpirun")
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
            tokio::process::Command::new(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        // No MPI config provided - serial mode
        tokio::process::Command::new(&exe_path)
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    // Write input to stdin
    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(input.as_bytes()).await
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
        let _ = app.emit("qe-output-line", format!("\nProcess exited with code: {:?}", status.code()));
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
    use tokio::io::{AsyncBufReadExt, BufReader, AsyncWriteExt};

    // Clone the path out of the guard before any await points
    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };

    let work_path = PathBuf::from(&working_dir);

    // Ensure working directory exists
    std::fs::create_dir_all(&work_path)
        .map_err(|e| format!("Failed to create working directory: {}", e))?;

    // Copy SCF's .save directory if project/calculation IDs are provided
    if let (Some(ref project_id), Some(ref scf_calc_id)) = (&config.project_id, &config.scf_calc_id) {
        let projects_dir = projects::get_projects_dir(&app)?;
        let scf_tmp_dir = projects_dir
            .join(project_id)
            .join("calculations")
            .join(scf_calc_id)
            .join("tmp");

        if scf_tmp_dir.exists() {
            let _ = app.emit("qe-output-line", "Copying SCF data to working directory...");

            // Copy everything from the SCF tmp dir (includes .save directory)
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;

            let _ = app.emit("qe-output-line", "SCF data copied successfully.");
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
    let band_path: Vec<qe::BandPathPoint> = config.k_path
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

    let _ = app.emit("qe-output-line", "=== Starting Band Structure Calculation ===");
    let _ = app.emit("qe-output-line", "Step 1/2: Running NSCF calculation along k-path...");

    // Step 2: Run pw.x for bands
    let exe_path = bin_dir.join("pw.x");
    if !exe_path.exists() {
        return Err("pw.x not found".to_string());
    }

    // Build the command - with or without MPI
    let mut child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            let _ = app.emit("qe-output-line", format!("Using MPI with {} processes", mpi.nprocs));
            tokio::process::Command::new("mpirun")
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio::process::Command::new(&exe_path)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        tokio::process::Command::new(&exe_path)
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    // Write input to stdin
    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(input.as_bytes()).await
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
    let _ = app.emit("qe-output-line", "Step 2/2: Running bands.x post-processing...");

    // Step 3: Run bands.x
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

    // Save bands.x input for reference
    std::fs::write(work_path.join("bands_pp.in"), &bands_x_input)
        .map_err(|e| format!("Failed to write bands.x input: {}", e))?;

    let mut bands_child = tokio::process::Command::new(&bands_x_path)
        .current_dir(&work_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start bands.x: {}", e))?;

    if let Some(mut stdin) = bands_child.stdin.take() {
        stdin.write_all(bands_x_input.as_bytes()).await
            .map_err(|e| format!("Failed to write bands.x input: {}", e))?;
    }

    // Stream bands.x output
    let bands_stdout = bands_child.stdout.take().ok_or("Failed to capture bands.x stdout")?;
    let mut bands_reader = BufReader::new(bands_stdout).lines();

    while let Some(line) = bands_reader.next_line().await.map_err(|e| e.to_string())? {
        let _ = app.emit("qe-output-line", &line);
    }

    let bands_status = bands_child.wait().await.map_err(|e| e.to_string())?;
    if !bands_status.success() {
        return Err(format!("bands.x failed with exit code: {:?}", bands_status.code()));
    }

    let _ = app.emit("qe-output-line", "");
    let _ = app.emit("qe-output-line", "Parsing band structure data...");

    // Step 4: Parse the output
    let gnu_file = work_path.join("bands.dat.gnu");
    if !gnu_file.exists() {
        return Err("bands.dat.gnu not found. bands.x may have failed.".to_string());
    }

    let ef = fermi_energy.unwrap_or(0.0);
    let mut band_data = read_bands_gnu_file(&gnu_file, ef)
        .map_err(|e| format!("Failed to parse band data: {}", e))?;

    // Add high-symmetry point markers
    qe::bands::add_symmetry_markers(&mut band_data, &config.k_path);

    let _ = app.emit("qe-output-line", format!("=== Band Structure Complete ==="));
    let _ = app.emit("qe-output-line", format!("  {} bands, {} k-points", band_data.n_bands, band_data.n_kpoints));
    if let Some(ref gap) = band_data.band_gap {
        let gap_type = if gap.is_direct { "direct" } else { "indirect" };
        let _ = app.emit("qe-output-line", format!("  Band gap: {:.3} eV ({})", gap.value, gap_type));
    }

    Ok(band_data)
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
            get_cpu_count,
            check_mpi_available,
            run_calculation,
            run_calculation_streaming,
            run_bands_calculation,
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
            projects::delete_calculation,
            projects::set_last_opened_cif,
            projects::get_cif_crystal_data,
            projects::get_cif_content,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
