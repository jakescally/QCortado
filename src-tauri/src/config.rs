//! Configuration persistence for QCortado
//!
//! Stores user settings in the app data directory.

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;
use tauri::{AppHandle, Manager};

/// Default MPI settings applied when creating new calculation runs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MpiDefaultsConfig {
    pub enabled: bool,
    pub nprocs: u32,
}

/// Application configuration stored on disk
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AppConfig {
    /// Path to QE bin directory
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qe_bin_dir: Option<String>,
    /// Optional command prefix prepended to all QE process launches
    /// Example: "mpirun" or "srun"
    #[serde(skip_serializing_if = "Option::is_none")]
    pub execution_prefix: Option<String>,
    /// Optional absolute path to the FermiSurfer executable
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fermi_surfer_path: Option<String>,
    /// Optional global defaults for MPI in calculation wizards
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mpi_defaults: Option<MpiDefaultsConfig>,
}

/// Gets the config file path
fn get_config_path(app: &AppHandle) -> Result<PathBuf, String> {
    let app_data = app
        .path()
        .app_data_dir()
        .map_err(|e| format!("Failed to get app data dir: {}", e))?;

    // Ensure directory exists
    if !app_data.exists() {
        fs::create_dir_all(&app_data)
            .map_err(|e| format!("Failed to create app data directory: {}", e))?;
    }

    Ok(app_data.join("config.json"))
}

/// Loads configuration from disk
pub fn load_config(app: &AppHandle) -> Result<AppConfig, String> {
    let config_path = get_config_path(app)?;

    if !config_path.exists() {
        return Ok(AppConfig::default());
    }

    let content =
        fs::read_to_string(&config_path).map_err(|e| format!("Failed to read config: {}", e))?;

    serde_json::from_str(&content).map_err(|e| format!("Failed to parse config: {}", e))
}

/// Saves configuration to disk
pub fn save_config(app: &AppHandle, config: &AppConfig) -> Result<(), String> {
    let config_path = get_config_path(app)?;

    let content = serde_json::to_string_pretty(config)
        .map_err(|e| format!("Failed to serialize config: {}", e))?;

    fs::write(&config_path, content).map_err(|e| format!("Failed to write config: {}", e))?;

    Ok(())
}

/// Updates a single config value and saves
pub fn update_qe_path(app: &AppHandle, path: Option<String>) -> Result<(), String> {
    let mut config = load_config(app)?;
    config.qe_bin_dir = path;
    save_config(app, &config)
}

/// Updates the global execution prefix and saves
pub fn update_execution_prefix(app: &AppHandle, prefix: Option<String>) -> Result<(), String> {
    let mut config = load_config(app)?;
    config.execution_prefix = prefix;
    save_config(app, &config)
}

/// Updates the FermiSurfer executable path and saves
pub fn update_fermi_surfer_path(app: &AppHandle, path: Option<String>) -> Result<(), String> {
    let mut config = load_config(app)?;
    config.fermi_surfer_path = path;
    save_config(app, &config)
}

/// Updates global MPI defaults and saves
pub fn update_mpi_defaults(
    app: &AppHandle,
    mpi_defaults: Option<MpiDefaultsConfig>,
) -> Result<(), String> {
    let mut config = load_config(app)?;
    config.mpi_defaults = mpi_defaults;
    save_config(app, &config)
}
