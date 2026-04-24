#![cfg_attr(feature = "viewer", allow(dead_code))]

//! QCortado - A modern UI for Quantum ESPRESSO
//!
//! This is the Tauri backend providing:
//! - QE input generation and validation
//! - Process execution with streaming output
//! - Output parsing and result extraction
//! - Project management

use base64::{engine::general_purpose::STANDARD as BASE64_STANDARD, Engine as _};
use flate2::read::GzDecoder;
use regex::Regex;
use std::collections::{HashMap, HashSet, VecDeque};
use std::ffi::{OsStr, OsString};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;
use std::time::Duration;
use tauri::{AppHandle, Emitter, Manager, State};

pub mod config;
pub mod hpc;
pub mod process_manager;
pub mod projects;
pub mod qe;
pub mod symmetry;

use process_manager::ProcessManager;

use qe::{
    add_phonon_symmetry_markers, build_epw_input, build_epw_keyword_map, build_transport_win,
    collect_epw_artifacts, epw_coarse_k_mesh, epw_coarse_q_mesh, epw_fine_k_mesh, epw_fine_q_mesh,
    export_ludwig_bundle, generate_dos_input, generate_matdyn_bands_input,
    generate_matdyn_dos_input, generate_ph_input, generate_q2r_input, parse_epw_result_v2,
    parse_ph_output, parse_transport_result, parse_wannier_hamiltonian,
    prepare_wannier_nscf_calculation, read_phonon_dispersion_file, read_phonon_dos_file,
    read_wannier_result, validate_epw_config, validate_transport_config, validate_wannier_config,
    DosCalculation, EpwArtifactManifestEntry, EpwCalculationConfig, EpwCalculationV1,
    EpwErrorRecord, EpwInputPreviewResult, EpwPrerequisiteValidation, EpwSourceRef, EpwSourcesV1,
    LudwigExportConfig, LudwigExportResult, MatdynCalculation, PhononPipelineConfig, PhononResult,
    Pw2Wannier90Config, Q2RCalculation, QPathPoint, TransportCalculationConfig, TransportResult,
    WannierCalculationConfig, WannierResult, EPW_SCHEMA_VERSION,
};
use qe::{
    generate_bands_x_input, generate_projwfc_input, generate_pw2wannier90_input, generate_pw_input,
    generate_wannier90_win, parse_dos_file, parse_projwfc_projection_groups_aligned,
    parse_pw_output, read_bands_gnu_file, BandData, BandsXConfig, KPathPoint, ProjwfcConfig,
    QECalculation, QEResult, QERunner,
};

/// Application state shared across commands.
pub struct AppState {
    /// Path to QE bin directory
    pub qe_bin_dir: Mutex<Option<PathBuf>>,
    /// Path to FermiSurfer executable
    pub fermi_surfer_path: Mutex<Option<PathBuf>>,
    /// Path to Wannier90 executable
    pub wannier90_path: Mutex<Option<PathBuf>>,
    /// Path to postw90.x executable
    pub postw90_path: Mutex<Option<PathBuf>>,
    /// Optional command prefix prepended before all QE launches
    pub execution_prefix: Mutex<Option<String>>,
    /// Optional global MPI defaults used by calculation wizards
    pub mpi_defaults: Mutex<Option<config::MpiDefaultsConfig>>,
    /// Global QE defaults used by calculation wizards/background fallbacks.
    pub qe_defaults: Mutex<config::QeDefaultsConfig>,
    /// Global project save-size mode
    pub save_size_mode: Mutex<config::SaveSizeMode>,
    /// Local or HPC execution mode
    pub execution_mode: Mutex<hpc::profile::ExecutionMode>,
    /// Saved HPC profiles
    pub hpc_profiles: Mutex<Vec<hpc::profile::HpcProfile>>,
    /// Active HPC profile ID
    pub active_hpc_profile_id: Mutex<Option<String>>,
    /// Current project directory
    pub project_dir: Mutex<Option<PathBuf>>,
    /// Background process manager
    pub process_manager: ProcessManager,
    /// Allows a programmatic exit to pass through ExitRequested interception.
    pub allow_exit: AtomicBool,
    /// Tracks whether a full shutdown workflow is currently in progress.
    pub exit_in_progress: AtomicBool,
    /// Whether automatic viewer-library publishing is enabled.
    pub viewer_auto_publish_enabled: AtomicBool,
    /// Debounce flag for queued viewer-library publishing.
    pub viewer_publish_pending: AtomicBool,
    /// Last-known viewer sync status, used by the viewer UI.
    pub viewer_sync_status: Mutex<hpc::viewer_library::ViewerSyncStatus>,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            qe_bin_dir: Mutex::new(None),
            fermi_surfer_path: Mutex::new(None),
            wannier90_path: Mutex::new(None),
            postw90_path: Mutex::new(None),
            execution_prefix: Mutex::new(None),
            mpi_defaults: Mutex::new(None),
            qe_defaults: Mutex::new(config::QeDefaultsConfig::default()),
            save_size_mode: Mutex::new(config::SaveSizeMode::Large),
            execution_mode: Mutex::new(hpc::profile::ExecutionMode::Local),
            hpc_profiles: Mutex::new(Vec::new()),
            active_hpc_profile_id: Mutex::new(None),
            project_dir: Mutex::new(None),
            process_manager: ProcessManager::new(),
            allow_exit: AtomicBool::new(false),
            exit_in_progress: AtomicBool::new(false),
            viewer_auto_publish_enabled: AtomicBool::new(true),
            viewer_publish_pending: AtomicBool::new(false),
            viewer_sync_status: Mutex::new(hpc::viewer_library::ViewerSyncStatus::default()),
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

fn normalize_mpi_defaults(
    defaults: Option<config::MpiDefaultsConfig>,
) -> Option<config::MpiDefaultsConfig> {
    defaults.map(|defaults| config::MpiDefaultsConfig {
        enabled: defaults.enabled,
        nprocs: defaults.nprocs.max(1),
    })
}

fn get_qe_smearing_default(state: &AppState) -> qe::SmearingType {
    state.qe_defaults.lock().unwrap().smearing.clone()
}

fn normalize_optional_path(path: Option<String>) -> Option<String> {
    path.and_then(|raw| {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    })
}

fn derive_postw90_path_from_wannier90_path(path: &Path) -> Option<PathBuf> {
    let parent = path.parent()?;
    Some(parent.join("postw90.x"))
}

fn resolve_local_postw90_path(state: &AppState) -> Result<PathBuf, String> {
    if let Some(wannier90_path) = state.wannier90_path.lock().unwrap().as_ref() {
        if let Some(candidate) = derive_postw90_path_from_wannier90_path(wannier90_path) {
            return Ok(candidate);
        }
    }

    // Legacy fallback for older configs that stored postw90.x separately.
    if let Some(path) = state.postw90_path.lock().unwrap().as_ref() {
        return Ok(path.clone());
    }

    Err("postw90.x path not configured. Configure the Wannier90 executable path first.".to_string())
}

fn derive_remote_postw90_path(remote_wannier90_path: Option<&str>) -> String {
    let remote_wannier90 = remote_wannier90_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or("wannier90.x");

    if remote_wannier90.contains('/') || remote_wannier90.starts_with('~') {
        let path = Path::new(remote_wannier90);
        if let Some(parent) = path.parent() {
            if parent.as_os_str().is_empty() || parent == Path::new(".") {
                "postw90.x".to_string()
            } else {
                parent.join("postw90.x").to_string_lossy().to_string()
            }
        } else {
            "postw90.x".to_string()
        }
    } else {
        "postw90.x".to_string()
    }
}

fn resolve_local_epw_path(qe_bin_dir: &Path) -> Option<PathBuf> {
    let primary = qe_bin_dir.join("epw.x");
    if primary.exists() {
        return Some(primary);
    }

    let nested_fallback = qe_bin_dir.join("EPW").join("bin").join("epw.x");
    if nested_fallback.exists() {
        return Some(nested_fallback);
    }

    let parent_fallback = qe_bin_dir
        .parent()
        .map(|parent| parent.join("EPW").join("bin").join("epw.x"));
    if let Some(path) = parent_fallback {
        if path.exists() {
            return Some(path);
        }
    }

    None
}

fn resolve_remote_epw_path_shell(qe_bin_dir: &str, remote_epw_path: Option<&str>) -> String {
    let remote_epw = remote_epw_path
        .map(str::trim)
        .filter(|value| !value.is_empty());

    if let Some(remote_epw) = remote_epw {
        if remote_epw.contains('/') || remote_epw.starts_with('~') {
            return format!(
                "tool={}; \
if [ \"$tool\" = \"~\" ]; then tool=\"$HOME\"; elif [ \"${{tool#~/}}\" != \"$tool\" ]; then tool=\"$HOME/${{tool#~/}}\"; fi; \
if [ -x \"$tool\" ]; then echo \"$tool\"; else echo missing; fi",
                shell_single_quote_local(remote_epw)
            );
        }

        return format!(
            "tool={}; if command -v \"$tool\" >/dev/null 2>&1; then command -v \"$tool\"; else echo missing; fi",
            shell_single_quote_local(remote_epw)
        );
    }

    let quoted_bin = shell_single_quote_local(qe_bin_dir);
    format!(
        "QE_BIN={bin}; \
if [ -x \"$QE_BIN/epw.x\" ]; then echo \"$QE_BIN/epw.x\"; \
elif [ -x \"$QE_BIN/EPW/bin/epw.x\" ]; then echo \"$QE_BIN/EPW/bin/epw.x\"; \
elif [ -x \"$(dirname \"$QE_BIN\")/EPW/bin/epw.x\" ]; then echo \"$(dirname \"$QE_BIN\")/EPW/bin/epw.x\"; \
else echo missing; fi",
        bin = quoted_bin
    )
}

fn normalize_hpc_text(input: &str, field: &str) -> Result<String, String> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        Err(format!("{} is required", field))
    } else {
        Ok(trimmed.to_string())
    }
}

fn now_iso() -> String {
    chrono::Utc::now().to_rfc3339()
}

fn shell_single_quote_local(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    let escaped = value.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

fn resolve_hpc_resource_type_for_resources(
    profile: &hpc::profile::HpcProfile,
    resources: Option<&hpc::profile::SlurmResourceRequest>,
) -> hpc::profile::ResourceType {
    hpc::slurm::merge_resources(profile, resources.cloned()).resource_type
}

fn build_hpc_launcher_command(
    profile: &hpc::profile::HpcProfile,
    resource_type: hpc::profile::ResourceType,
) -> String {
    let mut command = match profile.launcher {
        hpc::profile::HpcLauncher::Srun => "srun".to_string(),
        hpc::profile::HpcLauncher::Mpirun => "mpirun -np \"${SLURM_NTASKS:-1}\"".to_string(),
    };

    if let Some(extra) = profile.launcher_extra_args_for_resource(resource_type) {
        command.push(' ');
        command.push_str(extra);
    }

    command
}

fn resolve_hpc_qe_bin_dir_for_resources(
    profile: &hpc::profile::HpcProfile,
    resources: Option<&hpc::profile::SlurmResourceRequest>,
) -> String {
    let resource_type = resolve_hpc_resource_type_for_resources(profile, resources);
    profile
        .remote_qe_bin_dir_for_resource(resource_type)
        .trim_end_matches('/')
        .to_string()
}

fn command_args_include_pencil_decomposition(raw_args: &str) -> bool {
    normalize_cli_dash_text(raw_args)
        .split_whitespace()
        .any(|token| {
            let lower = token.to_ascii_lowercase();
            matches!(
                lower.as_str(),
                "-pd" | "-use_pd" | "-pencil_decomposition" | "-use_pencil_decomposition"
            ) || lower.starts_with("-pd=")
                || lower.starts_with("-use_pd=")
                || lower.starts_with("-pencil_decomposition=")
                || lower.starts_with("-use_pencil_decomposition=")
        })
}

fn qe_executable_uses_pencil_decomposition(executable: &str) -> bool {
    let executable_name = Path::new(executable)
        .file_name()
        .and_then(OsStr::to_str)
        .unwrap_or(executable)
        .to_ascii_lowercase();

    matches!(
        executable_name.as_str(),
        "pw.x"
            | "bands.x"
            | "projwfc.x"
            | "dos.x"
            | "fermi_velocity.x"
            | "ph.x"
            | "q2r.x"
            | "matdyn.x"
            | "pw2wannier90.x"
    )
}

fn build_hpc_qe_input_command(
    profile: &hpc::profile::HpcProfile,
    resource_type: hpc::profile::ResourceType,
    executable: &str,
    extra_args: Option<&str>,
    input_file: &str,
    output_file: &str,
) -> String {
    let mut command = format!(
        "{} \"$QE_BIN/{}\"",
        build_hpc_launcher_command(profile, resource_type),
        executable
    );
    let trimmed_extra_args = extra_args
        .map(|value| value.trim())
        .filter(|value| !value.is_empty());
    let has_pencil_decomposition_arg = trimmed_extra_args
        .map(command_args_include_pencil_decomposition)
        .unwrap_or(false);
    if let Some(args) = trimmed_extra_args {
        command.push(' ');
        command.push_str(args);
    }
    if qe_executable_uses_pencil_decomposition(executable) && !has_pencil_decomposition_arg {
        command.push_str(" -pd .true.");
    }
    command.push_str(&format!(" -in {} > {} 2>&1", input_file, output_file));
    command
}

fn build_hpc_logged_qe_step_command(
    profile: &hpc::profile::HpcProfile,
    resource_type: hpc::profile::ResourceType,
    step_label: &str,
    executable: &str,
    extra_args: Option<&str>,
    input_file: &str,
    output_file: &str,
) -> String {
    let qe_cmd = build_hpc_qe_input_command(
        profile,
        resource_type,
        executable,
        extra_args,
        input_file,
        output_file,
    );
    format!(
        "echo \"[QCortado] {step} started at $(date -u +%Y-%m-%dT%H:%M:%SZ)\"; \
__qcortado_step_start=$(date +%s); \
set +e; \
{cmd}; \
__qcortado_step_status=$?; \
set -e; \
__qcortado_step_end=$(date +%s); \
echo \"[QCortado] {step} finished with exit=${{__qcortado_step_status}} elapsed=$((__qcortado_step_end-__qcortado_step_start))s at $(date -u +%Y-%m-%dT%H:%M:%SZ)\"; \
[ $__qcortado_step_status -eq 0 ]",
        step = step_label,
        cmd = qe_cmd,
    )
}

fn build_hpc_logged_shell_step_command(step_label: &str, command: &str) -> String {
    format!(
        "echo \"[QCortado] {step} started at $(date -u +%Y-%m-%dT%H:%M:%SZ)\"; \
__qcortado_step_start=$(date +%s); \
set +e; \
{cmd}; \
__qcortado_step_status=$?; \
set -e; \
__qcortado_step_end=$(date +%s); \
echo \"[QCortado] {step} finished with exit=${{__qcortado_step_status}} elapsed=$((__qcortado_step_end-__qcortado_step_start))s at $(date -u +%Y-%m-%dT%H:%M:%SZ)\"; \
[ $__qcortado_step_status -eq 0 ]",
        step = step_label,
        cmd = command,
    )
}

fn normalize_cli_dash_text(input: &str) -> String {
    input
        .replace('\u{2014}', "--")
        .replace('\u{2013}', "-")
        .replace('\u{2010}', "-")
        .replace('\u{2011}', "-")
        .replace('\u{2012}', "-")
        .replace('\u{2212}', "-")
        .replace('\u{FE63}', "-")
        .replace('\u{FF0D}', "-")
}

fn sanitize_hpc_profile(
    mut profile: hpc::profile::HpcProfile,
) -> Result<hpc::profile::HpcProfile, String> {
    profile.name = normalize_hpc_text(&profile.name, "Profile name")?;
    profile.host = normalize_hpc_text(&profile.host, "SSH host")?;
    profile.username = normalize_hpc_text(&profile.username, "SSH username")?;
    profile.remote_qe_bin_dir =
        normalize_hpc_text(&profile.remote_qe_bin_dir, "Remote QE bin path")?;
    profile.remote_qe_cpu_bin_dir = sanitize_optional_hpc_field(profile.remote_qe_cpu_bin_dir);
    profile.remote_qe_gpu_bin_dir = sanitize_optional_hpc_field(profile.remote_qe_gpu_bin_dir);
    if profile.remote_qe_cpu_bin_dir.is_none() {
        profile.remote_qe_cpu_bin_dir = Some(profile.remote_qe_bin_dir.clone());
    }
    if profile.remote_qe_gpu_bin_dir.is_none() {
        profile.remote_qe_gpu_bin_dir = Some(profile.remote_qe_bin_dir.clone());
    }
    profile.remote_qe_bin_dir = profile
        .remote_qe_cpu_bin_dir
        .clone()
        .or(profile.remote_qe_gpu_bin_dir.clone())
        .unwrap_or(profile.remote_qe_bin_dir);
    profile.remote_epw_path = sanitize_optional_hpc_field(profile.remote_epw_path);
    profile.remote_wannier90_path = sanitize_optional_hpc_field(profile.remote_wannier90_path);
    profile.remote_postw90_path = None;
    profile.remote_pseudo_dir =
        normalize_hpc_text(&profile.remote_pseudo_dir, "Remote pseudo path")?;
    profile.remote_workspace_root =
        normalize_hpc_text(&profile.remote_workspace_root, "Remote workspace root")?;
    profile.remote_project_root =
        normalize_hpc_text(&profile.remote_project_root, "Remote project root")?;
    profile.launcher_extra_args = sanitize_optional_hpc_cli_field(profile.launcher_extra_args);
    profile.launcher_cpu_extra_args =
        sanitize_optional_hpc_cli_field(profile.launcher_cpu_extra_args);
    profile.launcher_gpu_extra_args =
        sanitize_optional_hpc_cli_field(profile.launcher_gpu_extra_args);
    if let Some(legacy_launcher_extra_args) = profile.launcher_extra_args.take() {
        if profile.launcher_cpu_extra_args.is_none() {
            profile.launcher_cpu_extra_args = Some(legacy_launcher_extra_args.clone());
        }
        if profile.launcher_gpu_extra_args.is_none() {
            profile.launcher_gpu_extra_args = Some(legacy_launcher_extra_args);
        }
    }
    if profile.port == 0 {
        profile.port = 22;
    }
    if profile.id.trim().is_empty() {
        profile.id = uuid::Uuid::new_v4().to_string();
        profile.created_at = now_iso();
    }
    profile.updated_at = now_iso();
    if profile.cluster.trim().is_empty() {
        profile.cluster = hpc::profile::default_cluster();
    }
    profile.default_cpu_resources = sanitize_hpc_resource_defaults(
        profile.default_cpu_resources,
        hpc::profile::ResourceType::Cpu,
    );
    profile.default_gpu_resources = sanitize_hpc_resource_defaults(
        profile.default_gpu_resources,
        hpc::profile::ResourceType::Gpu,
    );
    Ok(profile)
}

fn sanitize_optional_hpc_field(value: Option<String>) -> Option<String> {
    value.and_then(|raw| {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    })
}

fn sanitize_optional_hpc_cli_field(value: Option<String>) -> Option<String> {
    sanitize_optional_hpc_field(value).map(|value| normalize_cli_dash_text(&value))
}

fn sanitize_hpc_resource_defaults(
    mut resources: hpc::profile::SlurmResourceRequest,
    resource_type: hpc::profile::ResourceType,
) -> hpc::profile::SlurmResourceRequest {
    resources.resource_type = resource_type.clone();
    resources.partition = sanitize_optional_hpc_field(resources.partition);
    resources.walltime = sanitize_optional_hpc_field(resources.walltime);
    resources.qos = sanitize_optional_hpc_field(resources.qos);
    resources.account = sanitize_optional_hpc_field(resources.account);
    resources.constraint = sanitize_optional_hpc_field(resources.constraint);
    resources.module_preamble = sanitize_optional_hpc_cli_field(resources.module_preamble);
    resources.additional_sbatch = resources
        .additional_sbatch
        .into_iter()
        .map(|line| normalize_cli_dash_text(line.trim()))
        .filter(|line| !line.is_empty())
        .collect();
    resources.nodes = resources.nodes.map(|value| value.max(1));
    resources.ntasks = resources.ntasks.map(|value| value.max(1));
    resources.cpus_per_task = resources.cpus_per_task.map(|value| value.max(1));
    resources.memory_gb = resources.memory_gb.map(|value| value.max(1));
    resources.gpus = resources.gpus.map(|value| value.max(1));

    if resources.partition.is_none() {
        resources.partition = Some("short".to_string());
    }
    if resources.walltime.is_none() {
        resources.walltime = Some("02:00:00".to_string());
    }

    if matches!(resource_type, hpc::profile::ResourceType::Gpu) {
        if resources.gpus.is_none() {
            resources.gpus = Some(1);
        }
        if resources.nodes.is_none() {
            resources.nodes = Some(1);
        }
    } else {
        resources.gpus = Some(0);
    }

    resources
}

fn resolve_hpc_profile_from_state(
    state: &AppState,
    profile_id: Option<String>,
) -> Result<hpc::profile::HpcProfile, String> {
    let profiles = state.hpc_profiles.lock().unwrap();
    if profiles.is_empty() {
        return Err(
            "No HPC profiles are configured. Configure Andromeda in Settings first.".to_string(),
        );
    }

    let requested_id = profile_id.and_then(|raw| {
        let trimmed = raw.trim().to_string();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed)
        }
    });

    if let Some(profile_id) = requested_id {
        return profiles
            .iter()
            .find(|profile| profile.id == profile_id)
            .cloned()
            .ok_or_else(|| format!("HPC profile not found: {}", profile_id));
    }

    if let Some(active_id) = state.active_hpc_profile_id.lock().unwrap().clone() {
        if let Some(profile) = profiles.iter().find(|profile| profile.id == active_id) {
            return Ok(profile.clone());
        }
    }

    profiles
        .first()
        .cloned()
        .ok_or_else(|| "No HPC profile is available".to_string())
}

fn normalize_remote_path_for_tracking(raw: &str) -> Option<String> {
    let trimmed = raw.trim().trim_end_matches('/').to_string();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed)
    }
}

fn validate_remote_root_for_migration(path: &str, label: &str) -> Result<(), String> {
    let normalized = path.trim().trim_end_matches('/').to_string();
    if normalized.is_empty() {
        return Err(format!("{} is empty.", label));
    }
    if normalized == "/" || normalized == "." || normalized == "~" {
        return Err(format!(
            "{} '{}' is not a safe migration path.",
            label,
            path.trim()
        ));
    }
    let depth = normalized
        .split('/')
        .filter(|segment| !segment.is_empty())
        .count();
    if depth < 2 {
        return Err(format!(
            "{} '{}' is too shallow to migrate.",
            label,
            path.trim()
        ));
    }
    Ok(())
}

fn paths_overlap(path_a: &str, path_b: &str) -> bool {
    if path_a == path_b {
        return true;
    }
    path_a.starts_with(&format!("{}/", path_b)) || path_b.starts_with(&format!("{}/", path_a))
}

fn collect_saved_hpc_remote_paths(app: &AppHandle) -> Result<HashSet<String>, String> {
    let projects_dir = projects::ensure_projects_dir(app)?;
    let entries = std::fs::read_dir(&projects_dir).map_err(|e| {
        format!(
            "Failed to read projects directory {}: {}",
            projects_dir.display(),
            e
        )
    })?;
    let mut referenced_paths: HashSet<String> = HashSet::new();

    for entry in entries {
        let entry = entry.map_err(|e| format!("Failed to read projects directory entry: {}", e))?;
        let project_dir = entry.path();
        if !project_dir.is_dir() {
            continue;
        }
        let project_json = project_dir.join("project.json");
        if !project_json.exists() {
            continue;
        }

        let content = std::fs::read_to_string(&project_json).map_err(|e| {
            format!(
                "Failed to read project metadata {}: {}",
                project_json.display(),
                e
            )
        })?;
        let project: projects::Project = serde_json::from_str(&content).map_err(|e| {
            format!(
                "Failed to parse project metadata {}: {}",
                project_json.display(),
                e
            )
        })?;

        for calculation in project
            .cif_variants
            .iter()
            .flat_map(|variant| variant.calculations.iter())
        {
            let Some(parameters) = calculation.parameters.as_object() else {
                continue;
            };

            let backend = parameters
                .get("execution_backend")
                .and_then(|value| value.as_str())
                .map(|value| value.trim().to_ascii_lowercase())
                .unwrap_or_default();
            let remote_workdir = parameters
                .get("remote_workdir")
                .and_then(|value| value.as_str())
                .and_then(normalize_remote_path_for_tracking);
            let remote_project_path = parameters
                .get("remote_project_path")
                .and_then(|value| value.as_str())
                .and_then(normalize_remote_path_for_tracking);

            if backend != "hpc" && remote_workdir.is_none() && remote_project_path.is_none() {
                continue;
            }

            if let Some(path) = remote_workdir {
                referenced_paths.insert(path);
            }
            if let Some(path) = remote_project_path {
                referenced_paths.insert(path);
            }
        }
    }

    Ok(referenced_paths)
}

async fn resolve_remote_cleanup_path(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    raw_path: &str,
) -> Result<String, String> {
    resolve_remote_cleanup_path_with_timeout(profile, secret, raw_path, 120).await
}

async fn resolve_remote_cleanup_path_with_timeout(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    raw_path: &str,
    timeout_secs: u64,
) -> Result<String, String> {
    let trimmed = raw_path.trim();
    if trimmed.is_empty() {
        return Err("Remote path is empty".to_string());
    }

    if trimmed == "~" || trimmed.starts_with("~/") {
        let home = hpc::ssh::run_ssh_command_with_timeout(
            profile,
            secret,
            "printf %s \"$HOME\"",
            timeout_secs,
        )
        .await?;
        let home = home.trim();
        if home.is_empty() {
            return Err("Failed to resolve remote HOME directory".to_string());
        }
        if trimmed == "~" {
            return Ok(home.to_string());
        }
        let suffix = trimmed.trim_start_matches("~/");
        return Ok(format!("{}/{}", home.trim_end_matches('/'), suffix));
    }

    Ok(trimmed.to_string())
}

fn parse_recoverable_remote_phonon_runs(
    output: &str,
    profile_id: &str,
    location: &str,
) -> Vec<HpcRecoverableRemotePhononRun> {
    output
        .lines()
        .filter_map(|line| {
            let (raw_ts, raw_path) = line.split_once('\t')?;
            let remote_workdir = raw_path.trim().trim_end_matches('/').to_string();
            if remote_workdir.is_empty() {
                return None;
            }

            let modified_at_epoch = raw_ts.trim().parse::<i64>().unwrap_or(0);
            Some(HpcRecoverableRemotePhononRun {
                profile_id: profile_id.to_string(),
                remote_workdir,
                location: location.to_string(),
                modified_at_epoch,
            })
        })
        .collect()
}

fn build_recoverable_remote_phonon_scan_cmd(root: &str, name_filter: Option<&str>) -> String {
    let find_name_clause = name_filter
        .map(|pattern| format!(" -name {}", shell_single_quote_local(pattern)))
        .unwrap_or_default();
    format!(
        "if [ -d {root} ]; then \
            find {root} -mindepth 1 -maxdepth 1 -type d{find_name_clause} -print | \
            while IFS= read -r dir; do \
                if [ -z \"$dir\" ]; then continue; fi; \
                if ( [ -f \"$dir/phonon_freq.gp\" ] || [ -f \"$dir/phonon_freq\" ] ) && ls \"$dir\"/dynmat* >/dev/null 2>&1; then \
                    ts=$(stat -c %Y \"$dir\" 2>/dev/null || stat -f %m \"$dir\" 2>/dev/null || echo 0); \
                    printf '%s\\t%s\\n' \"$ts\" \"$dir\"; \
                fi; \
            done; \
         fi",
        root = shell_single_quote_local(root),
        find_name_clause = find_name_clause,
    )
}

fn build_remote_phonon_debug_probe_scan_cmd(root: &str, name_filter: Option<&str>) -> String {
    let find_name_clause = name_filter
        .map(|pattern| format!(" -name {}", shell_single_quote_local(pattern)))
        .unwrap_or_default();
    format!(
        "if [ -d {root} ]; then \
            find {root} -mindepth 1 -maxdepth 1 -type d{find_name_clause} -print | sort | \
            while IFS= read -r dir; do \
                if [ -z \"$dir\" ]; then continue; fi; \
                ts=$(stat -c %Y \"$dir\" 2>/dev/null || stat -f %m \"$dir\" 2>/dev/null || echo 0); \
                has_gp=0; [ -f \"$dir/phonon_freq.gp\" ] && has_gp=1; \
                has_freq=0; [ -f \"$dir/phonon_freq\" ] && has_freq=1; \
                dyn_count=$(find \"$dir\" -maxdepth 1 -type f -name 'dynmat*' | wc -l | tr -d ' '); \
                job_done=0; \
                if [ \"$dyn_count\" -gt 0 ] && grep -qs 'JOB DONE' \"$dir\"/dynmat* 2>/dev/null; then job_done=1; fi; \
                status_dynmatrix=0; \
                if [ -f \"$dir/status_run.xml\" ] && grep -qi '<current_stage>[[:space:]]*dynmatrix[[:space:]]*</current_stage>' \"$dir/status_run.xml\"; then status_dynmatrix=1; fi; \
                printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$ts\" \"$dir\" \"$has_gp\" \"$has_freq\" \"$dyn_count\" \"$job_done\" \"$status_dynmatrix\"; \
            done; \
         fi",
        root = shell_single_quote_local(root),
        find_name_clause = find_name_clause,
    )
}

fn parse_hpc_recovery_metadata(raw: &str) -> Option<HpcRecoveryJobMetadata> {
    let metadata = serde_json::from_str::<HpcRecoveryJobMetadata>(raw).ok()?;
    if metadata.schema_version != HPC_RECOVERY_METADATA_VERSION {
        return None;
    }
    if metadata.task_id.trim().is_empty() || metadata.task_kind.trim().is_empty() {
        return None;
    }
    Some(metadata)
}

fn infer_qcortado_task_kind(job_name: &str) -> Option<String> {
    let raw = job_name.trim().strip_prefix("qcortado-")?;
    let normalized = raw.trim().to_ascii_lowercase();
    if normalized.is_empty() {
        return None;
    }
    Some(normalized)
}

fn remote_path_is_under(path: &str, roots: &[String]) -> bool {
    let Some(normalized_path) = normalize_remote_path_for_tracking(path) else {
        return false;
    };
    roots.iter().any(|root| {
        normalize_remote_path_for_tracking(root)
            .map(|normalized_root| {
                normalized_path == normalized_root
                    || normalized_path
                        .starts_with(&format!("{}/", normalized_root.trim_end_matches('/')))
            })
            .unwrap_or(false)
    })
}

async fn read_remote_text_optional(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    remote_path: &str,
) -> Option<String> {
    let cmd = format!(
        "if [ -f {path} ]; then cat {path}; fi",
        path = shell_single_quote_local(remote_path)
    );
    hpc::ssh::run_ssh_command_with_timeout(profile, secret, &cmd, 8)
        .await
        .ok()
        .map(|value| value.trim_end_matches('\n').to_string())
        .filter(|value| !value.trim().is_empty())
}

async fn read_recovery_metadata_for_dir(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
) -> Option<HpcRecoveryJobMetadata> {
    let metadata_path = format!(
        "{}/{}",
        remote_workdir.trim_end_matches('/'),
        HPC_RECOVERY_METADATA_FILE
    );
    let raw = read_remote_text_optional(profile, secret, &metadata_path).await?;
    parse_hpc_recovery_metadata(&raw)
}

fn parse_squeue_headless_rows(
    output: &str,
) -> Vec<(
    String,
    String,
    String,
    Option<String>,
    Option<String>,
    Option<String>,
)> {
    output
        .lines()
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('|').collect();
            if parts.len() < 5 {
                return None;
            }
            let job_id = parts[0].trim().to_string();
            let job_name = parts[1].trim().to_string();
            let state = parts[2].trim().to_string();
            let node = parts
                .get(3)
                .map(|value| value.trim())
                .filter(|value| !value.is_empty() && *value != "(null)" && *value != "None")
                .map(str::to_string);
            let has_workdir_field = parts.len() >= 6;
            let workdir = if has_workdir_field {
                parts
                    .get(4)
                    .map(|value| value.trim())
                    .filter(|value| !value.is_empty() && *value != "N/A" && *value != "(null)")
                    .map(str::to_string)
            } else {
                None
            };
            let submitted_at_index = if has_workdir_field { 5 } else { 4 };
            let submitted_at = parts
                .get(submitted_at_index)
                .map(|value| value.trim())
                .filter(|value| !value.is_empty() && *value != "N/A")
                .map(str::to_string);
            if job_id.is_empty() || job_name.is_empty() {
                return None;
            }
            Some((job_id, job_name, state, node, workdir, submitted_at))
        })
        .collect()
}

async fn resolve_slurm_workdir(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    remote_job_id: &str,
) -> Option<String> {
    let cmd = format!(
        "scontrol show job {} 2>/dev/null | tr ' ' '\\n' | sed -n 's/^WorkDir=//p' | head -n 1",
        shell_single_quote_local(remote_job_id)
    );
    hpc::ssh::run_ssh_command_with_timeout(profile, secret, &cmd, 8)
        .await
        .ok()
        .map(|value| value.trim().trim_end_matches('/').to_string())
        .filter(|value| !value.is_empty())
}

async fn scheduler_snapshot_for_job(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    remote_job_id: &str,
) -> (String, Option<String>) {
    let squeue_cmd = format!(
        "squeue -h -j {} -o \"%T|%N\"",
        shell_single_quote_local(remote_job_id)
    );
    if let Ok(output) =
        hpc::ssh::run_ssh_command_with_timeout(profile, secret, &squeue_cmd, 8).await
    {
        if let Some(snapshot) = hpc::slurm::parse_squeue_snapshot(&output) {
            return (
                hpc::slurm::normalize_scheduler_state(&snapshot.state),
                snapshot.node,
            );
        }
    }
    let sacct_cmd = format!(
        "sacct -j {} --format=State,NodeList --parsable2 --noheader",
        shell_single_quote_local(remote_job_id)
    );
    if let Ok(output) = hpc::ssh::run_ssh_command_with_timeout(profile, secret, &sacct_cmd, 8).await
    {
        if let Some(snapshot) = hpc::slurm::parse_sacct_snapshot(&output) {
            return (
                hpc::slurm::normalize_scheduler_state(&snapshot.state),
                snapshot.node,
            );
        }
    }
    ("UNKNOWN".to_string(), None)
}

#[cfg(test)]
mod hpc_headless_recovery_tests {
    use super::*;

    #[test]
    fn parses_qcortado_squeue_rows() {
        let rows = parse_squeue_headless_rows(
            "123|qcortado-scf|RUNNING|node01|/scratch/qcortado_hpc_bundle_abc|2026-04-22T12:00:00\n\
             124|other|PENDING|(null)|N/A|N/A\n",
        );
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].0, "123");
        assert_eq!(rows[0].1, "qcortado-scf");
        assert_eq!(rows[0].3.as_deref(), Some("node01"));
        assert_eq!(
            rows[0].4.as_deref(),
            Some("/scratch/qcortado_hpc_bundle_abc")
        );
        assert_eq!(rows[1].3, None);
        assert_eq!(rows[1].4, None);
    }

    #[test]
    fn rejects_unsupported_recovery_metadata_schema() {
        let raw = r#"{
          "schema_version": 999,
          "task_id": "abc",
          "task_kind": "scf",
          "label": "SCF",
          "profile_id": "profile",
          "resource_type": "cpu",
          "slurm_job_name": "qcortado-scf",
          "submitted_at": "2026-04-22T12:00:00Z"
        }"#;
        assert!(parse_hpc_recovery_metadata(raw).is_none());
    }

    #[test]
    fn validates_remote_path_roots() {
        let roots = vec![
            "/scratch/user/qcortado".to_string(),
            "/project/user/qcortado".to_string(),
        ];
        assert!(remote_path_is_under(
            "/scratch/user/qcortado/qcortado_hpc_bundle_abc",
            &roots
        ));
        assert!(remote_path_is_under(
            "/project/user/qcortado/scf/abc",
            &roots
        ));
        assert!(!remote_path_is_under(
            "/scratch/user/qcortado-other/abc",
            &roots
        ));
        assert!(!remote_path_is_under("../qcortado_hpc_bundle_abc", &roots));
    }
}

async fn list_recoverable_remote_phonon_runs_for_profile(
    profile: &hpc::profile::HpcProfile,
    secret: Option<&str>,
    limit: usize,
) -> Result<Vec<HpcRecoverableRemotePhononRun>, String> {
    let workspace_root = resolve_remote_cleanup_path(
        profile,
        secret,
        profile.remote_workspace_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let project_root = resolve_remote_cleanup_path(
        profile,
        secret,
        profile.remote_project_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();

    let workspace_scan_cmd =
        build_recoverable_remote_phonon_scan_cmd(&workspace_root, Some("qcortado_hpc_bundle_*"));
    let workspace_output = hpc::ssh::run_ssh_command(profile, secret, &workspace_scan_cmd).await?;
    let workspace_runs =
        parse_recoverable_remote_phonon_runs(&workspace_output, &profile.id, "workspace");

    let project_phonon_root = format!("{}/phonon", project_root.trim_end_matches('/'));
    let project_scan_cmd = build_recoverable_remote_phonon_scan_cmd(&project_phonon_root, None);
    let project_output = hpc::ssh::run_ssh_command(profile, secret, &project_scan_cmd).await?;
    let project_runs =
        parse_recoverable_remote_phonon_runs(&project_output, &profile.id, "project_archive");

    let mut unique_runs: std::collections::HashMap<String, HpcRecoverableRemotePhononRun> =
        std::collections::HashMap::new();
    for run in workspace_runs.into_iter().chain(project_runs.into_iter()) {
        match unique_runs.get(&run.remote_workdir) {
            Some(existing) if existing.modified_at_epoch >= run.modified_at_epoch => {}
            _ => {
                unique_runs.insert(run.remote_workdir.clone(), run);
            }
        }
    }

    let mut runs: Vec<HpcRecoverableRemotePhononRun> = unique_runs.into_values().collect();
    runs.sort_by(|a, b| {
        b.modified_at_epoch
            .cmp(&a.modified_at_epoch)
            .then_with(|| a.remote_workdir.cmp(&b.remote_workdir))
    });
    if runs.len() > limit {
        runs.truncate(limit);
    }
    Ok(runs)
}

fn effective_execution_mode(
    state: &AppState,
    target: Option<&hpc::profile::ExecutionTarget>,
) -> hpc::profile::ExecutionMode {
    if let Some(target) = target {
        return target.mode;
    }
    *state.execution_mode.lock().unwrap()
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

fn is_explicit_command_path(command: &str) -> bool {
    Path::new(command).components().count() > 1
}

fn resolve_command_path(program: &OsStr) -> OsString {
    let program_path = Path::new(program);

    if program_path.components().count() > 1 {
        return program.to_os_string();
    }

    let Some(program_name) = program_path.file_name() else {
        return program.to_os_string();
    };

    if let Some(path_var) = std::env::var_os("PATH") {
        for path_dir in std::env::split_paths(&path_var) {
            let candidate = path_dir.join(program_name);
            if candidate.is_file() {
                return candidate.into_os_string();
            }
        }
    }

    #[cfg(target_os = "macos")]
    {
        for fallback_dir in ["/opt/homebrew/bin", "/usr/local/bin"] {
            let candidate = Path::new(fallback_dir).join(program_name);
            if candidate.is_file() {
                return candidate.into_os_string();
            }
        }
    }

    program.to_os_string()
}

fn tokio_command_with_prefix(
    program: impl AsRef<std::ffi::OsStr>,
    execution_prefix: Option<&str>,
) -> tokio::process::Command {
    let program_os = program.as_ref();
    let program_text = program_os.to_string_lossy().to_string();
    let resolved_program = resolve_command_path(program_os);

    if let Some(tokens) = parse_execution_prefix_tokens(execution_prefix) {
        let mut command = if command_basename(&tokens[0]) == command_basename(&program_text) {
            let executable = if is_explicit_command_path(&tokens[0]) {
                OsString::from(&tokens[0])
            } else {
                resolved_program.clone()
            };
            tokio::process::Command::new(executable)
        } else {
            let mut cmd =
                tokio::process::Command::new(resolve_command_path(OsStr::new(&tokens[0])));
            cmd.args(tokens.iter().skip(1));
            cmd.arg(&resolved_program);
            cmd.kill_on_drop(true);
            return cmd;
        };
        command.args(tokens.iter().skip(1));
        command.kill_on_drop(true);
        return command;
    }

    let mut command = tokio::process::Command::new(resolved_program);
    command.kill_on_drop(true);
    command
}

fn std_command_with_prefix(
    program: impl AsRef<std::ffi::OsStr>,
    execution_prefix: Option<&str>,
) -> std::process::Command {
    let program_os = program.as_ref();
    let program_text = program_os.to_string_lossy().to_string();
    let resolved_program = resolve_command_path(program_os);

    if let Some(tokens) = parse_execution_prefix_tokens(execution_prefix) {
        let mut command = if command_basename(&tokens[0]) == command_basename(&program_text) {
            let executable = if is_explicit_command_path(&tokens[0]) {
                OsString::from(&tokens[0])
            } else {
                resolved_program.clone()
            };
            std::process::Command::new(executable)
        } else {
            let mut cmd = std::process::Command::new(resolve_command_path(OsStr::new(&tokens[0])));
            cmd.args(tokens.iter().skip(1));
            cmd.arg(&resolved_program);
            return cmd;
        };
        command.args(tokens.iter().skip(1));
        return command;
    }

    std::process::Command::new(resolved_program)
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

fn can_safely_reset_working_dir(work_path: &Path) -> bool {
    work_path
        .file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.to_ascii_lowercase().starts_with("qcortado"))
        .unwrap_or(false)
}

fn prepare_working_directory(work_path: &Path, preserve_existing: bool) -> Result<(), String> {
    if !preserve_existing && work_path.exists() && can_safely_reset_working_dir(work_path) {
        std::fs::remove_dir_all(work_path).map_err(|e| {
            format!(
                "Failed to reset working directory {}: {}",
                work_path.display(),
                e
            )
        })?;
    }

    std::fs::create_dir_all(work_path).map_err(|e| {
        format!(
            "Failed to create working directory {}: {}",
            work_path.display(),
            e
        )
    })?;
    Ok(())
}

fn has_wavefunction_restart_files(save_dir: &Path) -> Result<bool, String> {
    let entries = std::fs::read_dir(save_dir).map_err(|e| {
        format!(
            "Failed to read SCF restart directory {}: {}",
            save_dir.display(),
            e
        )
    })?;

    for entry in entries {
        let entry = entry.map_err(|e| {
            format!(
                "Failed to inspect SCF restart directory {}: {}",
                save_dir.display(),
                e
            )
        })?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let file_name = entry.file_name().to_string_lossy().to_ascii_lowercase();
        if file_name.starts_with("wfc") {
            return Ok(true);
        }
    }

    Ok(false)
}

fn ensure_phonon_restart_inputs(work_path: &Path) -> Result<(), String> {
    let save_dir = work_path.join("tmp").join("qcortado_scf.save");
    if !save_dir.exists() {
        return Err(format!(
            "SCF restart directory not found at {}. \
Phonon calculations require the SCF .save data in ./tmp/qcortado_scf.save.",
            save_dir.display()
        ));
    }

    if !has_wavefunction_restart_files(&save_dir)? {
        return Err(format!(
            "SCF restart files are incomplete in {} (missing wfc* wavefunction files). \
This usually means the SCF was compacted for small storage. \
Phonon calculations require full SCF restart data. Re-run SCF with save size mode set to Large, then retry phonons.",
            save_dir.display()
        ));
    }

    Ok(())
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

/// Sets the path to the FermiSurfer executable.
#[tauri::command]
fn set_fermi_surfer_path(
    app: AppHandle,
    path: Option<String>,
    state: State<AppState>,
) -> Result<(), String> {
    let normalized = normalize_optional_path(path);

    let validated_path = if let Some(path_str) = normalized.clone() {
        let path_buf = PathBuf::from(&path_str);
        if !path_buf.exists() {
            return Err(format!(
                "FermiSurfer executable not found at {}",
                path_buf.display()
            ));
        }
        if !path_buf.is_file() {
            return Err(format!(
                "FermiSurfer path is not a file: {}",
                path_buf.display()
            ));
        }
        Some(path_buf)
    } else {
        None
    };

    *state.fermi_surfer_path.lock().unwrap() = validated_path;
    config::update_fermi_surfer_path(&app, normalized)
}

/// Gets the current FermiSurfer executable path.
#[tauri::command]
fn get_fermi_surfer_path(state: State<AppState>) -> Option<String> {
    state
        .fermi_surfer_path
        .lock()
        .unwrap()
        .as_ref()
        .map(|p| p.to_string_lossy().to_string())
}

/// Sets the path to the Wannier90 executable.
#[tauri::command]
fn set_wannier90_path(
    app: AppHandle,
    path: Option<String>,
    state: State<AppState>,
) -> Result<(), String> {
    let normalized = normalize_optional_path(path);

    let validated_path = if let Some(path_str) = normalized.clone() {
        let path_buf = PathBuf::from(&path_str);
        if !path_buf.exists() {
            return Err(format!(
                "Wannier90 executable not found at {}",
                path_buf.display()
            ));
        }
        if !path_buf.is_file() {
            return Err(format!(
                "Wannier90 path is not a file: {}",
                path_buf.display()
            ));
        }
        Some(path_buf)
    } else {
        None
    };

    *state.wannier90_path.lock().unwrap() = validated_path;
    config::update_wannier90_path(&app, normalized)
}

/// Gets the current Wannier90 executable path.
#[tauri::command]
fn get_wannier90_path(state: State<AppState>) -> Option<String> {
    state
        .wannier90_path
        .lock()
        .unwrap()
        .as_ref()
        .map(|p| p.to_string_lossy().to_string())
}

/// Sets the path to the postw90.x executable.
#[tauri::command]
fn set_postw90_path(
    app: AppHandle,
    path: Option<String>,
    state: State<AppState>,
) -> Result<(), String> {
    let normalized = normalize_optional_path(path);

    let validated_path = if let Some(path_str) = normalized.clone() {
        let path_buf = PathBuf::from(&path_str);
        if !path_buf.exists() {
            return Err(format!(
                "postw90.x executable not found at {}",
                path_buf.display()
            ));
        }
        if !path_buf.is_file() {
            return Err(format!(
                "postw90.x path is not a file: {}",
                path_buf.display()
            ));
        }
        Some(path_buf)
    } else {
        None
    };

    *state.postw90_path.lock().unwrap() = validated_path;
    config::update_postw90_path(&app, normalized)
}

/// Gets the current postw90.x executable path.
#[tauri::command]
fn get_postw90_path(state: State<AppState>) -> Option<String> {
    resolve_local_postw90_path(state.inner())
        .ok()
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

/// Sets global MPI defaults used to prefill wizard MPI controls.
#[tauri::command]
fn set_mpi_defaults(
    app: AppHandle,
    defaults: Option<config::MpiDefaultsConfig>,
    state: State<AppState>,
) -> Result<(), String> {
    let normalized = normalize_mpi_defaults(defaults);
    *state.mpi_defaults.lock().unwrap() = normalized.clone();
    config::update_mpi_defaults(&app, normalized)
}

/// Gets global MPI defaults used to prefill wizard MPI controls.
#[tauri::command]
fn get_mpi_defaults(state: State<AppState>) -> Option<config::MpiDefaultsConfig> {
    state.mpi_defaults.lock().unwrap().clone()
}

/// Sets global QE defaults used by wizard and fallback generation paths.
#[tauri::command]
fn set_qe_defaults(
    app: AppHandle,
    defaults: config::QeDefaultsConfig,
    state: State<AppState>,
) -> Result<(), String> {
    *state.qe_defaults.lock().unwrap() = defaults.clone();
    config::update_qe_defaults(&app, defaults)
}

/// Gets global QE defaults used by wizard and fallback generation paths.
#[tauri::command]
fn get_qe_defaults(state: State<AppState>) -> config::QeDefaultsConfig {
    state.qe_defaults.lock().unwrap().clone()
}

/// Sets global save-size mode used when persisting project calculation artifacts.
#[tauri::command]
fn set_save_size_mode(
    app: AppHandle,
    mode: config::SaveSizeMode,
    state: State<AppState>,
) -> Result<(), String> {
    *state.save_size_mode.lock().unwrap() = mode;
    config::update_save_size_mode(&app, mode)
}

/// Gets global save-size mode used when persisting project calculation artifacts.
#[tauri::command]
fn get_save_size_mode(state: State<AppState>) -> config::SaveSizeMode {
    *state.save_size_mode.lock().unwrap()
}

#[derive(Debug, serde::Serialize)]
struct HpcConnectionTestResult {
    success: bool,
    message: String,
}

#[derive(Debug, serde::Serialize)]
struct HpcScriptPreviewResult {
    script: String,
    sbatch_preview: String,
    validation: hpc::profile::ResourceValidation,
}

#[derive(Debug, serde::Serialize)]
struct HpcRemoteOrphanCleanupResult {
    profile_id: String,
    scanned_paths: usize,
    referenced_paths: usize,
    orphan_paths: usize,
    removed_paths: Vec<String>,
    failed_paths: Vec<String>,
}

#[derive(Debug, Clone, serde::Serialize)]
struct HpcRecoverableRemotePhononRun {
    profile_id: String,
    remote_workdir: String,
    location: String,
    modified_at_epoch: i64,
}

#[derive(Debug, serde::Serialize)]
struct HpcRemotePhononRecoveryDebugReport {
    profile_id: String,
    workspace_root: String,
    project_phonon_root: String,
    workspace_probe_output: String,
    project_probe_output: String,
    recoverable_runs: Vec<HpcRecoverableRemotePhononRun>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
struct HpcRecoveryJobMetadata {
    schema_version: u32,
    task_id: String,
    task_kind: String,
    label: String,
    profile_id: String,
    resource_type: String,
    slurm_job_name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    remote_job_id: Option<String>,
    submitted_at: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
}

const HPC_RECOVERY_METADATA_FILE: &str = ".qcortado_job.json";
const HPC_RECOVERY_METADATA_VERSION: u32 = 1;

#[derive(Debug, Clone, serde::Serialize)]
struct HpcHeadlessJobCandidate {
    profile_id: String,
    remote_job_id: String,
    job_name: String,
    task_kind: String,
    label: String,
    scheduler_state: String,
    remote_node: Option<String>,
    remote_workdir: Option<String>,
    submitted_at: Option<String>,
    metadata_status: String,
    auto_save_available: bool,
    project_id: Option<String>,
    cif_id: Option<String>,
}

#[derive(Debug, Clone, serde::Serialize)]
struct HpcAttachedJobResult {
    task_id: String,
    task_kind: String,
    remote_job_id: String,
    remote_workdir: String,
    auto_save_available: bool,
}

const HPC_PRESET_BUNDLE_KIND: &str = "qcortado_hpc_presets";
const HPC_PRESET_BUNDLE_VERSION: u32 = 1;
const IMPORTED_HPC_USERNAME_PLACEHOLDER: &str = "CHANGE_ME";
const HPC_REMOTE_PROJECT_TASK_KINDS: [&str; 7] = [
    "scf",
    "bands",
    "dos",
    "fermi_surface",
    "phonon",
    "transport",
    "epw",
];

fn default_hpc_preset_bundle_kind() -> String {
    HPC_PRESET_BUNDLE_KIND.to_string()
}

fn default_hpc_preset_bundle_port() -> u16 {
    22
}

fn normalize_hpc_profile_lookup_value(value: &str) -> String {
    value.trim().to_ascii_lowercase()
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
struct HpcPresetBundleProfile {
    #[serde(default)]
    id: String,
    name: String,
    #[serde(default = "hpc::profile::default_cluster")]
    cluster: String,
    #[serde(default = "default_hpc_preset_bundle_port")]
    port: u16,
    host: String,
    remote_qe_bin_dir: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    remote_qe_cpu_bin_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    remote_qe_gpu_bin_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    remote_epw_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    remote_wannier90_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    remote_postw90_path: Option<String>,
    remote_pseudo_dir: String,
    remote_workspace_root: String,
    remote_project_root: String,
    #[serde(default)]
    resource_mode: hpc::profile::HpcResourceMode,
    #[serde(default)]
    launcher: hpc::profile::HpcLauncher,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    launcher_cpu_extra_args: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    launcher_gpu_extra_args: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    launcher_extra_args: Option<String>,
    #[serde(default = "hpc::profile::default_cpu_resources")]
    default_cpu_resources: hpc::profile::SlurmResourceRequest,
    #[serde(default = "hpc::profile::default_gpu_resources")]
    default_gpu_resources: hpc::profile::SlurmResourceRequest,
    #[serde(default)]
    warning_acknowledged: bool,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
struct HpcPresetBundle {
    #[serde(default = "default_hpc_preset_bundle_kind")]
    kind: String,
    version: u32,
    exported_at: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    active_profile_id: Option<String>,
    #[serde(default)]
    profiles: Vec<HpcPresetBundleProfile>,
}

#[derive(Debug, serde::Serialize)]
struct HpcPresetBundleExportResult {
    bundle_path: String,
    profile_count: usize,
}

#[derive(Debug, serde::Serialize)]
struct HpcPresetBundleImportResult {
    imported_profile_count: usize,
    updated_profile_count: usize,
    created_profile_count: usize,
    profiles_requiring_username: Vec<String>,
    active_profile_id: Option<String>,
}

/// Sets the global execution mode (local or HPC).
#[tauri::command]
fn set_execution_mode(
    app: AppHandle,
    mode: hpc::profile::ExecutionMode,
    state: State<AppState>,
) -> Result<(), String> {
    *state.execution_mode.lock().unwrap() = mode;
    config::update_execution_mode(&app, mode)
}

/// Gets the current global execution mode.
#[tauri::command]
fn get_execution_mode(state: State<AppState>) -> hpc::profile::ExecutionMode {
    *state.execution_mode.lock().unwrap()
}

fn app_role() -> &'static str {
    if cfg!(feature = "viewer") {
        "viewer"
    } else {
        "research"
    }
}

#[tauri::command]
fn get_app_role() -> String {
    app_role().to_string()
}

#[tauri::command]
fn get_viewer_auto_publish_enabled(state: State<AppState>) -> bool {
    state.viewer_auto_publish_enabled.load(Ordering::SeqCst)
}

#[tauri::command]
fn set_viewer_auto_publish_enabled(
    app: AppHandle,
    state: State<AppState>,
    enabled: bool,
) -> Result<bool, String> {
    state
        .viewer_auto_publish_enabled
        .store(enabled, Ordering::SeqCst);
    config::update_viewer_auto_publish_enabled(&app, enabled)?;
    Ok(enabled)
}

async fn publish_viewer_library_with_profile(
    app: &AppHandle,
    state: &AppState,
    profile_id: Option<String>,
    project_id: Option<String>,
) -> Result<hpc::viewer_library::ViewerPublishResult, String> {
    let profile = resolve_hpc_profile_from_state(state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    match hpc::viewer_library::publish_viewer_library(
        app,
        &profile,
        secret.as_deref(),
        project_id.as_deref(),
    )
    .await
    {
        Ok(result) => {
            let _ =
                config::update_viewer_publish_status(app, Some(result.generated_at.clone()), None);
            Ok(result)
        }
        Err(err) => {
            let _ = config::update_viewer_publish_status(app, None, Some(err.clone()));
            Err(err)
        }
    }
}

pub(crate) fn queue_auto_viewer_publish(app: &AppHandle, state: &AppState) {
    if !state.viewer_auto_publish_enabled.load(Ordering::SeqCst) {
        return;
    }
    if state.viewer_publish_pending.swap(true, Ordering::SeqCst) {
        return;
    }

    let app = app.clone();
    tauri::async_runtime::spawn(async move {
        tokio::time::sleep(Duration::from_millis(1800)).await;
        let state = app.state::<AppState>();
        let _ = publish_viewer_library_with_profile(&app, state.inner(), None, None).await;
        state.viewer_publish_pending.store(false, Ordering::SeqCst);
    });
}

#[tauri::command]
async fn hpc_publish_viewer_library(
    app: AppHandle,
    profile_id: Option<String>,
    project_id: Option<String>,
    state: State<'_, AppState>,
) -> Result<hpc::viewer_library::ViewerPublishResult, String> {
    publish_viewer_library_with_profile(&app, state.inner(), profile_id, project_id).await
}

async fn sync_viewer_library_with_profile(
    app: &AppHandle,
    state: &AppState,
    profile_id: Option<String>,
) -> Result<hpc::viewer_library::ViewerSyncResult, String> {
    let profile = resolve_hpc_profile_from_state(state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    match hpc::viewer_library::sync_viewer_library(app, &profile, secret.as_deref()).await {
        Ok(result) => {
            {
                let mut status = state.viewer_sync_status.lock().unwrap();
                status.last_synced_at = Some(result.synced_at.clone());
                status.last_error = None;
                status.local_project_count = result.total_projects;
            }
            let _ = config::update_viewer_sync_status(app, Some(result.synced_at.clone()), None);
            Ok(result)
        }
        Err(err) => {
            {
                let mut status = state.viewer_sync_status.lock().unwrap();
                status.last_error = Some(err.clone());
            }
            let _ = config::update_viewer_sync_status(app, None, Some(err.clone()));
            Err(err)
        }
    }
}

#[tauri::command]
async fn viewer_sync_remote_library(
    app: AppHandle,
    profile_id: Option<String>,
    state: State<'_, AppState>,
) -> Result<hpc::viewer_library::ViewerSyncResult, String> {
    sync_viewer_library_with_profile(&app, state.inner(), profile_id).await
}

#[tauri::command]
fn viewer_get_sync_status(state: State<AppState>) -> hpc::viewer_library::ViewerSyncStatus {
    state.viewer_sync_status.lock().unwrap().clone()
}

/// Lists configured HPC profiles.
#[tauri::command]
fn hpc_list_profiles(state: State<AppState>) -> Vec<hpc::profile::HpcProfile> {
    state.hpc_profiles.lock().unwrap().clone()
}

/// Exports non-sensitive HPC presets/defaults for sharing across environments.
#[tauri::command]
fn hpc_export_preset_bundle(
    destination_path: String,
    state: State<AppState>,
) -> Result<HpcPresetBundleExportResult, String> {
    let normalized_destination_path = destination_path.trim();
    if normalized_destination_path.is_empty() {
        return Err("Destination path is required".to_string());
    }

    let profiles = state.hpc_profiles.lock().unwrap().clone();
    if profiles.is_empty() {
        return Err("No HPC profiles are configured to export.".to_string());
    }
    let active_profile_id = state.active_hpc_profile_id.lock().unwrap().clone();

    let bundle = HpcPresetBundle {
        kind: default_hpc_preset_bundle_kind(),
        version: HPC_PRESET_BUNDLE_VERSION,
        exported_at: now_iso(),
        active_profile_id,
        profiles: profiles
            .into_iter()
            .map(|profile| HpcPresetBundleProfile {
                id: profile.id,
                name: profile.name,
                cluster: profile.cluster,
                port: profile.port,
                host: profile.host,
                remote_qe_bin_dir: profile.remote_qe_bin_dir,
                remote_qe_cpu_bin_dir: profile.remote_qe_cpu_bin_dir,
                remote_qe_gpu_bin_dir: profile.remote_qe_gpu_bin_dir,
                remote_epw_path: profile.remote_epw_path,
                remote_wannier90_path: profile.remote_wannier90_path,
                remote_postw90_path: profile.remote_postw90_path,
                remote_pseudo_dir: profile.remote_pseudo_dir,
                remote_workspace_root: profile.remote_workspace_root,
                remote_project_root: profile.remote_project_root,
                resource_mode: profile.resource_mode,
                launcher: profile.launcher,
                launcher_cpu_extra_args: profile.launcher_cpu_extra_args,
                launcher_gpu_extra_args: profile.launcher_gpu_extra_args,
                launcher_extra_args: profile.launcher_extra_args,
                default_cpu_resources: profile.default_cpu_resources,
                default_gpu_resources: profile.default_gpu_resources,
                warning_acknowledged: profile.warning_acknowledged,
            })
            .collect(),
    };

    let bundle_json = serde_json::to_string_pretty(&bundle)
        .map_err(|e| format!("Failed to serialize HPC preset bundle: {}", e))?;
    let bundle_path = PathBuf::from(normalized_destination_path);

    if let Some(parent) = bundle_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                format!(
                    "Failed to create bundle directory {}: {}",
                    parent.display(),
                    e
                )
            })?;
        }
    }

    std::fs::write(&bundle_path, bundle_json).map_err(|e| {
        format!(
            "Failed to write bundle file {}: {}",
            bundle_path.display(),
            e
        )
    })?;

    Ok(HpcPresetBundleExportResult {
        bundle_path: bundle_path.display().to_string(),
        profile_count: bundle.profiles.len(),
    })
}

/// Imports non-sensitive HPC presets/defaults and merges them into local profiles.
#[tauri::command]
fn hpc_import_preset_bundle(
    app: AppHandle,
    bundle_path: String,
    state: State<AppState>,
) -> Result<HpcPresetBundleImportResult, String> {
    let normalized_bundle_path = bundle_path.trim();
    if normalized_bundle_path.is_empty() {
        return Err("Bundle path is required".to_string());
    }

    let bundle_content = std::fs::read_to_string(normalized_bundle_path).map_err(|e| {
        format!(
            "Failed to read bundle file {}: {}",
            normalized_bundle_path, e
        )
    })?;
    let bundle: HpcPresetBundle = serde_json::from_str(&bundle_content)
        .map_err(|e| format!("Failed to parse bundle file: {}", e))?;

    if bundle.kind != HPC_PRESET_BUNDLE_KIND {
        return Err(format!(
            "Unsupported HPC preset bundle kind '{}'.",
            bundle.kind
        ));
    }
    if bundle.version != HPC_PRESET_BUNDLE_VERSION {
        return Err(format!(
            "Unsupported HPC preset bundle version {}. Expected {}.",
            bundle.version, HPC_PRESET_BUNDLE_VERSION
        ));
    }
    if bundle.profiles.is_empty() {
        return Err("Bundle does not contain any HPC profiles.".to_string());
    }

    let imported_profile_count = bundle.profiles.len();
    let imported_active_profile_id = bundle
        .active_profile_id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string());
    let mut resolved_active_profile_id: Option<String> = None;
    let mut profiles_requiring_username: Vec<String> = Vec::new();
    let mut updated_profile_count = 0usize;
    let mut created_profile_count = 0usize;

    let mut profiles = state.hpc_profiles.lock().unwrap();
    let mut active_profile_id = state.active_hpc_profile_id.lock().unwrap();
    let mut used_profile_ids: HashSet<String> =
        profiles.iter().map(|profile| profile.id.clone()).collect();

    for imported_profile in bundle.profiles {
        let imported_profile_id = imported_profile.id.trim().to_string();
        let imported_cluster_lookup = if imported_profile.cluster.trim().is_empty() {
            normalize_hpc_profile_lookup_value(&hpc::profile::default_cluster())
        } else {
            normalize_hpc_profile_lookup_value(&imported_profile.cluster)
        };
        let imported_name_lookup = normalize_hpc_profile_lookup_value(&imported_profile.name);
        let imported_host_lookup = normalize_hpc_profile_lookup_value(&imported_profile.host);

        let matching_index = if !imported_profile_id.is_empty() {
            profiles
                .iter()
                .position(|profile| profile.id == imported_profile_id)
        } else {
            None
        }
        .or_else(|| {
            if imported_name_lookup.is_empty() || imported_host_lookup.is_empty() {
                return None;
            }
            profiles.iter().position(|profile| {
                normalize_hpc_profile_lookup_value(&profile.name) == imported_name_lookup
                    && normalize_hpc_profile_lookup_value(&profile.host) == imported_host_lookup
                    && normalize_hpc_profile_lookup_value(&profile.cluster)
                        == imported_cluster_lookup
            })
        });

        if let Some(index) = matching_index {
            let existing_profile = profiles[index].clone();
            let mut merged_profile = hpc::profile::HpcProfile {
                id: existing_profile.id.clone(),
                name: imported_profile.name,
                cluster: imported_profile.cluster,
                port: imported_profile.port,
                host: imported_profile.host,
                username: existing_profile.username.clone(),
                auth_method: existing_profile.auth_method.clone(),
                ssh_key_path: existing_profile.ssh_key_path.clone(),
                remote_qe_bin_dir: imported_profile.remote_qe_bin_dir,
                remote_qe_cpu_bin_dir: imported_profile.remote_qe_cpu_bin_dir,
                remote_qe_gpu_bin_dir: imported_profile.remote_qe_gpu_bin_dir,
                remote_epw_path: imported_profile.remote_epw_path,
                remote_wannier90_path: imported_profile.remote_wannier90_path,
                remote_postw90_path: imported_profile.remote_postw90_path,
                remote_pseudo_dir: imported_profile.remote_pseudo_dir,
                remote_workspace_root: imported_profile.remote_workspace_root,
                remote_project_root: imported_profile.remote_project_root,
                resource_mode: imported_profile.resource_mode,
                launcher: imported_profile.launcher,
                launcher_cpu_extra_args: imported_profile.launcher_cpu_extra_args,
                launcher_gpu_extra_args: imported_profile.launcher_gpu_extra_args,
                launcher_extra_args: imported_profile.launcher_extra_args,
                default_cpu_resources: imported_profile.default_cpu_resources,
                default_gpu_resources: imported_profile.default_gpu_resources,
                credential_persisted: existing_profile.credential_persisted,
                warning_acknowledged: imported_profile.warning_acknowledged,
                created_at: existing_profile.created_at.clone(),
                updated_at: existing_profile.updated_at.clone(),
            };
            merged_profile = sanitize_hpc_profile(merged_profile)?;
            merged_profile.credential_persisted = existing_profile.credential_persisted;
            profiles[index] = merged_profile.clone();
            updated_profile_count += 1;

            if imported_active_profile_id.as_deref() == Some(imported_profile_id.as_str()) {
                resolved_active_profile_id = Some(merged_profile.id.clone());
            }
            continue;
        }

        let mut next_profile_id = if !imported_profile_id.is_empty() {
            imported_profile_id.clone()
        } else {
            uuid::Uuid::new_v4().to_string()
        };
        while used_profile_ids.contains(&next_profile_id) {
            next_profile_id = uuid::Uuid::new_v4().to_string();
        }
        used_profile_ids.insert(next_profile_id.clone());

        let mut imported_entry = hpc::profile::HpcProfile {
            id: next_profile_id,
            name: imported_profile.name,
            cluster: imported_profile.cluster,
            port: imported_profile.port,
            host: imported_profile.host,
            username: IMPORTED_HPC_USERNAME_PLACEHOLDER.to_string(),
            auth_method: hpc::profile::HpcAuthMethod::SshKey,
            ssh_key_path: None,
            remote_qe_bin_dir: imported_profile.remote_qe_bin_dir,
            remote_qe_cpu_bin_dir: imported_profile.remote_qe_cpu_bin_dir,
            remote_qe_gpu_bin_dir: imported_profile.remote_qe_gpu_bin_dir,
            remote_epw_path: imported_profile.remote_epw_path,
            remote_wannier90_path: imported_profile.remote_wannier90_path,
            remote_postw90_path: imported_profile.remote_postw90_path,
            remote_pseudo_dir: imported_profile.remote_pseudo_dir,
            remote_workspace_root: imported_profile.remote_workspace_root,
            remote_project_root: imported_profile.remote_project_root,
            resource_mode: imported_profile.resource_mode,
            launcher: imported_profile.launcher,
            launcher_cpu_extra_args: imported_profile.launcher_cpu_extra_args,
            launcher_gpu_extra_args: imported_profile.launcher_gpu_extra_args,
            launcher_extra_args: imported_profile.launcher_extra_args,
            default_cpu_resources: imported_profile.default_cpu_resources,
            default_gpu_resources: imported_profile.default_gpu_resources,
            credential_persisted: false,
            warning_acknowledged: imported_profile.warning_acknowledged,
            created_at: now_iso(),
            updated_at: now_iso(),
        };
        imported_entry = sanitize_hpc_profile(imported_entry)?;
        imported_entry.credential_persisted = false;
        profiles_requiring_username.push(imported_entry.name.clone());
        created_profile_count += 1;

        if imported_active_profile_id.as_deref() == Some(imported_profile_id.as_str()) {
            resolved_active_profile_id = Some(imported_entry.id.clone());
        }

        profiles.push(imported_entry);
    }

    if let Some(active_id_from_bundle) = imported_active_profile_id {
        if let Some(mapped_active_profile_id) = resolved_active_profile_id {
            *active_profile_id = Some(mapped_active_profile_id);
        } else if profiles
            .iter()
            .any(|profile| profile.id == active_id_from_bundle)
        {
            *active_profile_id = Some(active_id_from_bundle);
        }
    }

    if active_profile_id.is_none() {
        *active_profile_id = profiles.first().map(|profile| profile.id.clone());
    }

    config::update_hpc_profiles(&app, profiles.clone(), active_profile_id.clone())?;

    Ok(HpcPresetBundleImportResult {
        imported_profile_count,
        updated_profile_count,
        created_profile_count,
        profiles_requiring_username,
        active_profile_id: active_profile_id.clone(),
    })
}

/// Saves or updates an HPC profile and optional credential secret.
#[tauri::command]
fn hpc_save_profile(
    app: AppHandle,
    profile: hpc::profile::HpcProfile,
    credential: Option<String>,
    persist_credential: Option<bool>,
    state: State<AppState>,
) -> Result<hpc::profile::HpcProfile, String> {
    let persist = persist_credential.unwrap_or(false);
    let mut profile = sanitize_hpc_profile(profile)?;
    profile.credential_persisted = persist;

    if let Some(secret) = credential.as_ref() {
        let trimmed = secret.trim().to_string();
        if !trimmed.is_empty() {
            hpc::credentials::set_session_secret(&profile.id, trimmed.clone());
            if persist {
                hpc::credentials::save_persisted_secret(
                    &profile.id,
                    &profile.username,
                    &profile.host,
                    &trimmed,
                )?;
            }
        }
    } else if !persist {
        hpc::credentials::clear_session_secret(&profile.id);
        let _ = hpc::credentials::delete_persisted_secret(
            &profile.id,
            &profile.username,
            &profile.host,
        );
    }

    let mut profiles = state.hpc_profiles.lock().unwrap();
    let mut active_id = state.active_hpc_profile_id.lock().unwrap();
    if let Some(existing) = profiles.iter_mut().find(|entry| entry.id == profile.id) {
        profile.created_at = existing.created_at.clone();
        *existing = profile.clone();
    } else {
        if profile.created_at.trim().is_empty() {
            profile.created_at = now_iso();
        }
        profiles.push(profile.clone());
    }

    if active_id.is_none() {
        *active_id = Some(profile.id.clone());
    }

    config::update_hpc_profiles(&app, profiles.clone(), active_id.clone())?;
    Ok(profile)
}

/// Updates only default Slurm resources for an existing HPC profile.
#[tauri::command]
fn hpc_update_profile_defaults(
    app: AppHandle,
    profile_id: String,
    resource_mode: Option<hpc::profile::HpcResourceMode>,
    launcher: Option<hpc::profile::HpcLauncher>,
    launcher_cpu_extra_args: Option<String>,
    launcher_gpu_extra_args: Option<String>,
    default_cpu_resources: hpc::profile::SlurmResourceRequest,
    default_gpu_resources: hpc::profile::SlurmResourceRequest,
    state: State<AppState>,
) -> Result<hpc::profile::HpcProfile, String> {
    let normalized_profile_id = profile_id.trim();
    if normalized_profile_id.is_empty() {
        return Err("HPC profile ID is required".to_string());
    }

    let mut profiles = state.hpc_profiles.lock().unwrap();
    let mut active_id = state.active_hpc_profile_id.lock().unwrap();

    let profile = profiles
        .iter_mut()
        .find(|entry| entry.id == normalized_profile_id)
        .ok_or_else(|| format!("HPC profile not found: {}", normalized_profile_id))?;

    if let Some(mode) = resource_mode {
        profile.resource_mode = mode;
    }
    if let Some(launcher_value) = launcher {
        profile.launcher = launcher_value;
    }
    profile.launcher_cpu_extra_args = sanitize_optional_hpc_cli_field(launcher_cpu_extra_args);
    profile.launcher_gpu_extra_args = sanitize_optional_hpc_cli_field(launcher_gpu_extra_args);
    profile.launcher_extra_args = None;
    profile.default_cpu_resources =
        sanitize_hpc_resource_defaults(default_cpu_resources, hpc::profile::ResourceType::Cpu);
    profile.default_gpu_resources =
        sanitize_hpc_resource_defaults(default_gpu_resources, hpc::profile::ResourceType::Gpu);
    profile.updated_at = now_iso();
    let updated_profile = profile.clone();

    if active_id.is_none() {
        *active_id = Some(updated_profile.id.clone());
    }

    config::update_hpc_profiles(&app, profiles.clone(), active_id.clone())?;
    Ok(updated_profile)
}

/// Migrates remote workspace/project roots for an HPC profile.
#[tauri::command]
async fn hpc_migrate_remote_roots(
    app: AppHandle,
    profile_id: String,
    new_workspace_root: String,
    new_project_root: String,
    state: State<'_, AppState>,
) -> Result<hpc::profile::HpcProfile, String> {
    let normalized_profile_id = profile_id.trim().to_string();
    if normalized_profile_id.is_empty() {
        return Err("HPC profile ID is required".to_string());
    }

    let normalized_workspace_root =
        normalize_hpc_text(&new_workspace_root, "New remote workspace root")?;
    let normalized_project_root = normalize_hpc_text(&new_project_root, "New remote project root")?;

    let current_profile = {
        let profiles = state.hpc_profiles.lock().unwrap();
        profiles
            .iter()
            .find(|entry| entry.id == normalized_profile_id)
            .cloned()
            .ok_or_else(|| format!("HPC profile not found: {}", normalized_profile_id))?
    };

    if current_profile.remote_workspace_root == normalized_workspace_root
        && current_profile.remote_project_root == normalized_project_root
    {
        return Ok(current_profile);
    }

    let secret = hpc::credentials::resolve_secret(
        &current_profile.id,
        &current_profile.username,
        &current_profile.host,
        current_profile.credential_persisted,
    )?;

    let old_workspace_resolved = resolve_remote_cleanup_path(
        &current_profile,
        secret.as_deref(),
        current_profile.remote_workspace_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let old_project_resolved = resolve_remote_cleanup_path(
        &current_profile,
        secret.as_deref(),
        current_profile.remote_project_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let new_workspace_resolved = resolve_remote_cleanup_path(
        &current_profile,
        secret.as_deref(),
        normalized_workspace_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let new_project_resolved = resolve_remote_cleanup_path(
        &current_profile,
        secret.as_deref(),
        normalized_project_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();

    validate_remote_root_for_migration(&old_workspace_resolved, "Current workspace root")?;
    validate_remote_root_for_migration(&old_project_resolved, "Current project root")?;
    validate_remote_root_for_migration(&new_workspace_resolved, "New workspace root")?;
    validate_remote_root_for_migration(&new_project_resolved, "New project root")?;

    let should_migrate_workspace = old_workspace_resolved != new_workspace_resolved;
    let should_migrate_project = old_project_resolved != new_project_resolved;

    if should_migrate_workspace && paths_overlap(&old_workspace_resolved, &new_workspace_resolved) {
        return Err("Workspace roots overlap. Choose a non-overlapping destination.".to_string());
    }
    if should_migrate_project && paths_overlap(&old_project_resolved, &new_project_resolved) {
        return Err("Project roots overlap. Choose a non-overlapping destination.".to_string());
    }
    if (should_migrate_workspace || should_migrate_project)
        && paths_overlap(&new_workspace_resolved, &new_project_resolved)
    {
        return Err(
            "New workspace root and new project root must not overlap each other.".to_string(),
        );
    }

    if should_migrate_workspace || should_migrate_project {
        let mut script_lines: Vec<String> = vec!["set -euo pipefail".to_string()];

        if should_migrate_workspace {
            script_lines.push(format!(
                "mkdir -p {}",
                shell_single_quote_local(&new_workspace_resolved)
            ));
            script_lines.push(format!(
                "if [ -d {old} ]; then cp -a {old}/. {new}/; fi",
                old = shell_single_quote_local(&old_workspace_resolved),
                new = shell_single_quote_local(&new_workspace_resolved),
            ));
        }
        if should_migrate_project {
            script_lines.push(format!(
                "mkdir -p {}",
                shell_single_quote_local(&new_project_resolved)
            ));
            script_lines.push(format!(
                "if [ -d {old} ]; then cp -a {old}/. {new}/; fi",
                old = shell_single_quote_local(&old_project_resolved),
                new = shell_single_quote_local(&new_project_resolved),
            ));
        }
        if should_migrate_workspace {
            script_lines.push(format!(
                "if [ -d {old} ]; then rm -rf -- {old}; fi",
                old = shell_single_quote_local(&old_workspace_resolved),
            ));
        }
        if should_migrate_project {
            script_lines.push(format!(
                "if [ -d {old} ]; then rm -rf -- {old}; fi",
                old = shell_single_quote_local(&old_project_resolved),
            ));
        }

        let remote_script = script_lines.join(" && ");
        hpc::ssh::run_ssh_command(&current_profile, secret.as_deref(), &remote_script)
            .await
            .map_err(|e| format!("Failed to migrate remote roots: {}", e))?;
    }

    let updated_profile = {
        let mut profiles = state.hpc_profiles.lock().unwrap();
        let mut active_id = state.active_hpc_profile_id.lock().unwrap();
        let profile = profiles
            .iter_mut()
            .find(|entry| entry.id == normalized_profile_id)
            .ok_or_else(|| format!("HPC profile not found: {}", normalized_profile_id))?;

        profile.remote_workspace_root = normalized_workspace_root;
        profile.remote_project_root = normalized_project_root;
        profile.updated_at = now_iso();
        let updated = profile.clone();

        if active_id.is_none() {
            *active_id = Some(updated.id.clone());
        }

        config::update_hpc_profiles(&app, profiles.clone(), active_id.clone())?;
        updated
    };

    Ok(updated_profile)
}

/// Deletes an HPC profile and any stored credentials.
#[tauri::command]
fn hpc_delete_profile(
    app: AppHandle,
    profile_id: String,
    state: State<AppState>,
) -> Result<(), String> {
    let mut profiles = state.hpc_profiles.lock().unwrap();
    let mut active_id = state.active_hpc_profile_id.lock().unwrap();

    let index = profiles
        .iter()
        .position(|entry| entry.id == profile_id)
        .ok_or_else(|| format!("HPC profile not found: {}", profile_id))?;
    let removed = profiles.remove(index);

    hpc::credentials::clear_session_secret(&removed.id);
    let _ =
        hpc::credentials::delete_persisted_secret(&removed.id, &removed.username, &removed.host);

    if active_id.as_ref() == Some(&profile_id) {
        *active_id = profiles.first().map(|entry| entry.id.clone());
    }

    config::update_hpc_profiles(&app, profiles.clone(), active_id.clone())?;
    Ok(())
}

/// Sets the active HPC profile ID.
#[tauri::command]
fn hpc_set_active_profile(
    app: AppHandle,
    profile_id: String,
    state: State<AppState>,
) -> Result<(), String> {
    let profiles = state.hpc_profiles.lock().unwrap();
    if !profiles.iter().any(|profile| profile.id == profile_id) {
        return Err(format!("HPC profile not found: {}", profile_id));
    }
    drop(profiles);
    *state.active_hpc_profile_id.lock().unwrap() = Some(profile_id.clone());
    config::update_hpc_profiles(
        &app,
        state.hpc_profiles.lock().unwrap().clone(),
        Some(profile_id),
    )
}

/// Gets the active HPC profile ID.
#[tauri::command]
fn hpc_get_active_profile_id(state: State<AppState>) -> Option<String> {
    state.active_hpc_profile_id.lock().unwrap().clone()
}

/// Tests SSH connection for the selected HPC profile.
#[tauri::command]
async fn hpc_test_connection(
    profile_id: Option<String>,
    state: State<'_, AppState>,
) -> Result<HpcConnectionTestResult, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;
    let output = hpc::ssh::run_ssh_command(
        &profile,
        secret.as_deref(),
        "echo QCORTADO_HPC_CONNECTION_OK",
    )
    .await;

    match output {
        Ok(result) => Ok(HpcConnectionTestResult {
            success: result.contains("QCORTADO_HPC_CONNECTION_OK"),
            message: result.trim().to_string(),
        }),
        Err(err) => Ok(HpcConnectionTestResult {
            success: false,
            message: err,
        }),
    }
}

/// Validates SSH, scheduler tools, QE binary path, and remote workspace write access.
#[tauri::command]
async fn hpc_validate_environment(
    profile_id: Option<String>,
    state: State<'_, AppState>,
) -> Result<hpc::profile::HpcEnvironmentValidation, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let mut messages: Vec<String> = Vec::new();

    let reachable = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), "echo reachable")
        .await
        .map(|value| value.contains("reachable"))
        .unwrap_or(false);
    if !reachable {
        messages.push(
            "Unable to connect to cluster via SSH. Verify VPN/campus network and credentials."
                .to_string(),
        );
    }

    let check_cmd = |tool: &str| {
        format!(
            "command -v {} >/dev/null 2>&1 && echo ok || echo missing",
            tool
        )
    };

    let sbatch_available =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &check_cmd("sbatch"))
            .await
            .map(|value| value.contains("ok"))
            .unwrap_or(false);
    if !sbatch_available {
        messages.push("sbatch not available in remote shell PATH".to_string());
    }

    let squeue_available =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &check_cmd("squeue"))
            .await
            .map(|value| value.contains("ok"))
            .unwrap_or(false);
    if !squeue_available {
        messages.push("squeue not available in remote shell PATH".to_string());
    }

    let sacct_available =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &check_cmd("sacct"))
            .await
            .map(|value| value.contains("ok"))
            .unwrap_or(false);
    if !sacct_available {
        messages.push("sacct not available in remote shell PATH".to_string());
    }

    let cpu_qe_bin_dir = profile
        .remote_qe_bin_dir_for_resource(hpc::profile::ResourceType::Cpu)
        .trim_end_matches('/')
        .to_string();
    let gpu_qe_bin_dir = profile
        .remote_qe_bin_dir_for_resource(hpc::profile::ResourceType::Gpu)
        .trim_end_matches('/')
        .to_string();
    let qe_path_checks: Vec<(String, String)> = match profile.resource_mode {
        hpc::profile::HpcResourceMode::CpuOnly => vec![("CPU".to_string(), cpu_qe_bin_dir)],
        hpc::profile::HpcResourceMode::GpuOnly => vec![("GPU".to_string(), gpu_qe_bin_dir)],
        hpc::profile::HpcResourceMode::Both => {
            if cpu_qe_bin_dir == gpu_qe_bin_dir {
                vec![("CPU/GPU".to_string(), cpu_qe_bin_dir)]
            } else {
                vec![
                    ("CPU".to_string(), cpu_qe_bin_dir),
                    ("GPU".to_string(), gpu_qe_bin_dir),
                ]
            }
        }
    };

    let mut qe_pw_available = true;
    let mut qe_epw_available = true;
    let mut qe_pw2wannier_available = true;
    let remote_epw_override = profile
        .remote_epw_path
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty());
    for (role, qe_bin_dir) in &qe_path_checks {
        let pw_available_for_role = hpc::ssh::run_ssh_command(
            &profile,
            secret.as_deref(),
            &format!("test -x {}/pw.x && echo ok || echo missing", qe_bin_dir),
        )
        .await
        .map(|value| value.contains("ok"))
        .unwrap_or(false);
        if !pw_available_for_role {
            messages.push(format!(
                "pw.x not found/executable at {} ({} QE path)",
                qe_bin_dir, role
            ));
        }
        qe_pw_available &= pw_available_for_role;

        let pw2wannier_available_for_role = hpc::ssh::run_ssh_command(
            &profile,
            secret.as_deref(),
            &format!(
                "test -x {}/pw2wannier90.x && echo ok || echo missing",
                qe_bin_dir
            ),
        )
        .await
        .map(|value| value.contains("ok"))
        .unwrap_or(false);
        if !pw2wannier_available_for_role {
            messages.push(format!(
                "pw2wannier90.x not found/executable at {} ({} QE path)",
                qe_bin_dir, role
            ));
        }
        qe_pw2wannier_available &= pw2wannier_available_for_role;
    }

    if let Some(remote_epw_override) = remote_epw_override {
        let epw_override_resolved = hpc::ssh::run_ssh_command(
            &profile,
            secret.as_deref(),
            &resolve_remote_epw_path_shell("", Some(remote_epw_override)),
        )
        .await
        .ok()
        .and_then(|value| {
            let trimmed = value.trim();
            if trimmed.is_empty() || trimmed == "missing" {
                None
            } else {
                Some(trimmed.to_string())
            }
        });
        qe_epw_available = epw_override_resolved.is_some();
        if !qe_epw_available {
            messages.push(format!(
                "epw.x not found/executable at configured path '{}'. EPW requires manual compilation.",
                remote_epw_override
            ));
        }
    } else {
        for (role, qe_bin_dir) in &qe_path_checks {
            let epw_resolved_for_role = hpc::ssh::run_ssh_command(
                &profile,
                secret.as_deref(),
                &resolve_remote_epw_path_shell(qe_bin_dir, None),
            )
            .await
            .ok()
            .and_then(|value| {
                let trimmed = value.trim();
                if trimmed.is_empty() || trimmed == "missing" {
                    None
                } else {
                    Some(trimmed.to_string())
                }
            });
            let epw_available_for_role = epw_resolved_for_role.is_some();
            if !epw_available_for_role {
                messages.push(format!(
                    "epw.x not found/executable at {} or fallback EPW/bin locations ({} QE path). EPW requires manual compilation.",
                    qe_bin_dir, role
                ));
            }
            qe_epw_available &= epw_available_for_role;
        }
    }

    let remote_wannier90 = profile
        .remote_wannier90_path
        .as_deref()
        .map(|value| value.trim())
        .filter(|value| !value.is_empty())
        .unwrap_or("wannier90.x");
    let wannier90_check = if remote_wannier90.contains('/') || remote_wannier90.starts_with('~') {
        format!(
            "tool={}; \
if [ \"$tool\" = \"~\" ]; then tool=\"$HOME\"; elif [ \"${{tool#~/}}\" != \"$tool\" ]; then tool=\"$HOME/${{tool#~/}}\"; fi; \
test -x \"$tool\" && echo ok || echo missing",
            shell_single_quote_local(remote_wannier90)
        )
    } else {
        format!(
            "command -v {} >/dev/null 2>&1 && echo ok || echo missing",
            remote_wannier90
        )
    };
    let wannier90_available =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &wannier90_check)
            .await
            .map(|value| value.contains("ok"))
            .unwrap_or(false);
    if !wannier90_available {
        messages.push(format!(
            "wannier90.x not found/executable at {}",
            remote_wannier90
        ));
    }

    let remote_postw90 = derive_remote_postw90_path(profile.remote_wannier90_path.as_deref());
    let postw90_check = if remote_postw90.contains('/') || remote_postw90.starts_with('~') {
        format!(
            "tool={}; \
if [ \"$tool\" = \"~\" ]; then tool=\"$HOME\"; elif [ \"${{tool#~/}}\" != \"$tool\" ]; then tool=\"$HOME/${{tool#~/}}\"; fi; \
test -x \"$tool\" && echo ok || echo missing",
            shell_single_quote_local(&remote_postw90)
        )
    } else {
        format!(
            "command -v {} >/dev/null 2>&1 && echo ok || echo missing",
            remote_postw90
        )
    };
    let postw90_available = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &postw90_check)
        .await
        .map(|value| value.contains("ok"))
        .unwrap_or(false);
    if !postw90_available {
        messages.push(format!(
            "postw90.x not found/executable at {}",
            remote_postw90
        ));
    }

    let probe_file = format!(
        "{}/.qcortado_probe_{}",
        profile.remote_workspace_root.trim_end_matches('/'),
        uuid::Uuid::new_v4()
    );
    let workspace_writable = hpc::ssh::run_ssh_command(
        &profile,
        secret.as_deref(),
        &format!(
            "mkdir -p {} && touch {} && rm -f {} && echo ok || echo missing",
            profile.remote_workspace_root.trim_end_matches('/'),
            probe_file,
            probe_file
        ),
    )
    .await
    .map(|value| value.contains("ok"))
    .unwrap_or(false);
    if !workspace_writable {
        messages.push(format!(
            "Workspace is not writable: {}",
            profile.remote_workspace_root
        ));
    }

    Ok(hpc::profile::HpcEnvironmentValidation {
        reachable,
        sbatch_available,
        squeue_available,
        sacct_available,
        qe_pw_available,
        qe_epw_available,
        qe_pw2wannier_available,
        wannier90_available,
        postw90_available,
        workspace_writable,
        messages,
    })
}

/// Queries cluster-wide node + queue activity snapshot for Andromeda.
#[tauri::command]
async fn hpc_get_cluster_snapshot(
    profile_id: Option<String>,
    queue_scope: Option<String>,
    include_queue: Option<bool>,
    queue_limit: Option<u32>,
    state: State<'_, AppState>,
) -> Result<hpc::cluster_snapshot::HpcClusterSnapshot, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    if !hpc::cluster_snapshot::is_andromeda_profile(&profile) {
        return Err(format!(
            "Node Activity currently supports Andromeda only (profile cluster='{}', host='{}').",
            profile.cluster, profile.host
        ));
    }

    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let resolved_scope = hpc::cluster_snapshot::normalize_queue_scope(queue_scope);
    let include_queue = include_queue.unwrap_or(true);
    let queue_limit = queue_limit.unwrap_or(1500).clamp(1, 5000) as usize;
    let remote_cmd = hpc::cluster_snapshot::build_cluster_snapshot_command(
        resolved_scope,
        include_queue,
        queue_limit,
    );
    let raw_output = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &remote_cmd).await?;
    let parsed = hpc::cluster_snapshot::parse_cluster_snapshot_output(&raw_output, include_queue)?;

    Ok(hpc::cluster_snapshot::HpcClusterSnapshot {
        captured_at: now_iso(),
        cluster: profile.cluster.trim().to_string(),
        host: profile.host.trim().to_string(),
        queue_scope: resolved_scope.as_str().to_string(),
        queue_included: include_queue,
        queue_limit,
        nodes: parsed.nodes,
        queue: parsed.queue,
        warnings: parsed.warnings,
    })
}

/// Samples structured remote utilization for the active HPC run.
#[tauri::command]
async fn hpc_sample_utilization(
    profile_id: Option<String>,
    remote_job_id: Option<String>,
    remote_node: Option<String>,
    resource_type: Option<hpc::profile::ResourceType>,
    state: State<'_, AppState>,
) -> Result<hpc::utilization::HpcUtilizationSample, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let job_id = remote_job_id.as_deref().unwrap_or("").trim().to_string();
    let node_hint = remote_node.as_deref().unwrap_or("").trim().to_string();
    let resolved_resource_type = resource_type.unwrap_or(hpc::profile::ResourceType::Cpu);
    if job_id.is_empty() {
        return Ok(hpc::utilization::HpcUtilizationSample {
            captured_at: now_iso(),
            resource_type: resolved_resource_type,
            job_id: None,
            node: if node_hint.is_empty() {
                None
            } else {
                Some(node_hint)
            },
            sources: Vec::new(),
            warnings: vec!["Waiting for remote job allocation.".to_string()],
            scheduler: None,
            cpu: None,
            memory: None,
            gpu: None,
            raw: None,
        });
    }

    let telemetry_cmd =
        hpc::utilization::build_utilization_command(&job_id, resolved_resource_type);
    let raw_output = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &telemetry_cmd).await?;
    let parsed = hpc::utilization::parse_utilization_output(&raw_output, resolved_resource_type);

    let parsed_node = parsed
        .scheduler
        .as_ref()
        .and_then(|scheduler| scheduler.node.clone());
    Ok(hpc::utilization::HpcUtilizationSample {
        captured_at: now_iso(),
        resource_type: resolved_resource_type,
        job_id: Some(job_id),
        node: parsed_node.or_else(|| (!node_hint.is_empty()).then_some(node_hint)),
        sources: parsed.sources,
        warnings: parsed.warnings,
        scheduler: parsed.scheduler,
        cpu: parsed.cpu,
        memory: parsed.memory,
        gpu: parsed.gpu,
        raw: Some(raw_output),
    })
}

/// Lists pseudopotential files (`*.UPF`/`*.upf`) from the configured remote pseudo directory.
#[tauri::command]
async fn hpc_list_remote_pseudopotentials(
    profile_id: Option<String>,
    pseudo_dir: Option<String>,
    state: State<'_, AppState>,
) -> Result<Vec<String>, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let remote_pseudo_dir = pseudo_dir
        .as_deref()
        .unwrap_or(profile.remote_pseudo_dir.as_str())
        .trim()
        .to_string();
    if remote_pseudo_dir.is_empty() {
        return Err("Remote pseudopotential directory is not configured".to_string());
    }

    let list_cmd = format!(
        "dir={dir}; \
if [ \"$dir\" = \"~\" ]; then \
  dir=\"$HOME\"; \
elif [ \"${{dir#~/}}\" != \"$dir\" ]; then \
  dir=\"$HOME/${{dir#~/}}\"; \
fi; \
if [ ! -d \"$dir\" ]; then \
  echo \"__QCORTADO_PSEUDO_DIR_MISSING__:$dir\"; \
  exit 0; \
fi; \
for file in \"$dir\"/*.UPF \"$dir\"/*.upf; do \
  if [ -f \"$file\" ]; then basename \"$file\"; fi; \
done | LC_ALL=C sort -u",
        dir = shell_single_quote_local(&remote_pseudo_dir)
    );

    let output = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &list_cmd).await?;
    let mut pseudos: Vec<String> = Vec::new();
    let mut seen = HashSet::new();
    for line in output.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(path) = trimmed.strip_prefix("__QCORTADO_PSEUDO_DIR_MISSING__:") {
            return Err(format!(
                "Remote pseudopotential directory not found: {}",
                path
            ));
        }
        if seen.insert(trimmed.to_string()) {
            pseudos.push(trimmed.to_string());
        }
    }

    pseudos.sort();
    Ok(pseudos)
}

/// Lists remote pseudopotentials and parses SOC/cutoff metadata from their headers.
#[tauri::command]
async fn hpc_list_remote_pseudopotential_metadata(
    profile_id: Option<String>,
    pseudo_dir: Option<String>,
    state: State<'_, AppState>,
) -> Result<Vec<PseudopotentialMetadata>, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let remote_pseudo_dir = pseudo_dir
        .as_deref()
        .unwrap_or(profile.remote_pseudo_dir.as_str())
        .trim()
        .to_string();
    if remote_pseudo_dir.is_empty() {
        return Err("Remote pseudopotential directory is not configured".to_string());
    }

    let list_cmd = format!(
        "dir={dir}; \
if [ \"$dir\" = \"~\" ]; then \
  dir=\"$HOME\"; \
elif [ \"${{dir#~/}}\" != \"$dir\" ]; then \
  dir=\"$HOME/${{dir#~/}}\"; \
fi; \
if [ ! -d \"$dir\" ]; then \
  echo \"__QCORTADO_PSEUDO_DIR_MISSING__:$dir\"; \
  exit 0; \
fi; \
for file in \"$dir\"/*.UPF \"$dir\"/*.upf; do \
  if [ ! -f \"$file\" ]; then \
    continue; \
  fi; \
  name=$(basename \"$file\"); \
  stem=${{name%.*}}; \
  echo \"__QCORTADO_REMOTE_METADATA_FILE__:upf_b64:$name\"; \
  head -n 200 \"$file\" | base64; \
  echo \"__QCORTADO_REMOTE_METADATA_FILE_END__\"; \
  if [ -f \"$dir/$stem.djrepo\" ]; then \
    echo \"__QCORTADO_REMOTE_METADATA_FILE__:djrepo_b64:$stem.djrepo\"; \
    base64 < \"$dir/$stem.djrepo\"; \
    echo \"__QCORTADO_REMOTE_METADATA_FILE_END__\"; \
  elif [ -f \"$dir/$stem.djrepo.gz\" ]; then \
    echo \"__QCORTADO_REMOTE_METADATA_FILE__:djrepo_gz_b64:$stem.djrepo.gz\"; \
    gzip -dc \"$dir/$stem.djrepo.gz\" | base64; \
    echo \"__QCORTADO_REMOTE_METADATA_FILE_END__\"; \
  fi; \
done",
        dir = shell_single_quote_local(&remote_pseudo_dir)
    );

    let output = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &list_cmd).await?;
    parse_remote_pseudopotential_metadata_output(&output)
}

/// Loads remote SSSP JSON data from the configured remote pseudo directory.
#[tauri::command]
async fn hpc_load_remote_sssp_data(
    profile_id: Option<String>,
    pseudo_dir: Option<String>,
    state: State<'_, AppState>,
) -> Result<std::collections::HashMap<String, SSSPElementData>, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let remote_pseudo_dir = pseudo_dir
        .as_deref()
        .unwrap_or(profile.remote_pseudo_dir.as_str())
        .trim()
        .to_string();
    if remote_pseudo_dir.is_empty() {
        return Err("Remote pseudopotential directory is not configured".to_string());
    }

    let load_cmd = format!(
        "dir={dir}; \
if [ \"$dir\" = \"~\" ]; then \
  dir=\"$HOME\"; \
elif [ \"${{dir#~/}}\" != \"$dir\" ]; then \
  dir=\"$HOME/${{dir#~/}}\"; \
fi; \
if [ ! -d \"$dir\" ]; then \
  echo \"__QCORTADO_PSEUDO_DIR_MISSING__:$dir\"; \
  exit 0; \
fi; \
sssp_file=$(find \"$dir\" -maxdepth 1 -type f \\( -name 'SSSP*.json' -o -name 'sssp*.json' \\) | LC_ALL=C sort | head -n 1); \
if [ -z \"$sssp_file\" ]; then \
  echo \"__QCORTADO_SSSP_NOT_FOUND__\"; \
  exit 0; \
fi; \
echo \"__QCORTADO_SSSP_FILE__:$sssp_file\"; \
cat \"$sssp_file\"",
        dir = shell_single_quote_local(&remote_pseudo_dir)
    );

    let output = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &load_cmd).await?;
    let mut json_lines = Vec::new();
    let mut sssp_file: Option<String> = None;

    for line in output.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            json_lines.push(line.to_string());
            continue;
        }
        if let Some(path) = trimmed.strip_prefix("__QCORTADO_PSEUDO_DIR_MISSING__:") {
            return Err(format!(
                "Remote pseudopotential directory not found: {}",
                path
            ));
        }
        if trimmed == "__QCORTADO_SSSP_NOT_FOUND__" {
            return Err("No SSSP JSON file found in remote pseudo directory".to_string());
        }
        if let Some(path) = trimmed.strip_prefix("__QCORTADO_SSSP_FILE__:") {
            sssp_file = Some(path.to_string());
            continue;
        }
        json_lines.push(line.to_string());
    }

    let json_content = json_lines.join("\n");
    if json_content.trim().is_empty() {
        return Err(match sssp_file {
            Some(path) => format!("Remote SSSP file is empty: {}", path),
            None => "No SSSP JSON data returned from remote pseudo directory".to_string(),
        });
    }

    serde_json::from_str::<std::collections::HashMap<String, SSSPElementData>>(&json_content)
        .map_err(|e| {
            if let Some(path) = sssp_file.as_deref() {
                format!("Failed to parse remote SSSP JSON ({}): {}", path, e)
            } else {
                format!("Failed to parse remote SSSP JSON: {}", e)
            }
        })
}

/// Generates a preview Slurm script for the active profile/resources.
#[tauri::command]
fn hpc_preview_slurm_script(
    profile_id: Option<String>,
    task_kind: String,
    command_lines: Vec<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    state: State<AppState>,
) -> Result<HpcScriptPreviewResult, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let commands = if command_lines.is_empty() {
        vec!["echo \"Add commands for this task in wizard settings\"".to_string()]
    } else {
        command_lines
    };
    let script = hpc::slurm::build_slurm_script(
        &profile,
        &format!("qcortado-{}", task_kind.trim().to_ascii_lowercase()),
        &commands,
        resources,
    );
    Ok(HpcScriptPreviewResult {
        script: script.script,
        sbatch_preview: script.sbatch_preview,
        validation: script.validation,
    })
}

/// Opens a read-only HPC activity popout window.
#[tauri::command]
async fn hpc_list_headless_jobs(
    profile_id: Option<String>,
    limit: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Vec<HpcHeadlessJobCandidate>, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;
    let resolved_limit = limit.unwrap_or(50).clamp(1, 200) as usize;
    let workspace_root = resolve_remote_cleanup_path_with_timeout(
        &profile,
        secret.as_deref(),
        profile.remote_workspace_root.trim_end_matches('/'),
        8,
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let mut allowed_roots = vec![workspace_root.clone()];
    if let Ok(project_root) = resolve_remote_cleanup_path_with_timeout(
        &profile,
        secret.as_deref(),
        profile.remote_project_root.trim_end_matches('/'),
        5,
    )
    .await
    {
        let trimmed = project_root.trim_end_matches('/').to_string();
        if !trimmed.is_empty() {
            allowed_roots.push(trimmed);
        }
    }
    let attached_jobs: HashSet<String> = state
        .process_manager
        .list_tasks()
        .await
        .into_iter()
        .filter_map(|task| task.remote_job_id)
        .collect();

    let mut candidates_by_key: HashMap<String, HpcHeadlessJobCandidate> = HashMap::new();
    let squeue_cmd = format!(
        "squeue -h -u {} -o \"%i|%j|%T|%N|%V\"",
        shell_single_quote_local(&profile.username)
    );
    let squeue_output =
        hpc::ssh::run_ssh_command_with_timeout(&profile, secret.as_deref(), &squeue_cmd, 10)
            .await?;
    for (job_id, job_name, state_name, node, raw_workdir, submitted_at) in
        parse_squeue_headless_rows(&squeue_output)
    {
        if attached_jobs.contains(&job_id) {
            continue;
        }
        let Some(task_kind) = infer_qcortado_task_kind(&job_name) else {
            continue;
        };
        let remote_workdir = match raw_workdir {
            Some(path) if remote_path_is_under(&path, &allowed_roots) => {
                Some(path.trim_end_matches('/').to_string())
            }
            _ => resolve_slurm_workdir(&profile, secret.as_deref(), &job_id)
                .await
                .filter(|path| remote_path_is_under(path, &allowed_roots)),
        };
        let metadata = if let Some(remote_workdir_value) = remote_workdir.as_deref() {
            read_recovery_metadata_for_dir(&profile, secret.as_deref(), remote_workdir_value).await
        } else {
            None
        };
        let label = metadata
            .as_ref()
            .map(|value| value.label.clone())
            .filter(|value| !value.trim().is_empty())
            .unwrap_or_else(|| task_kind.to_ascii_uppercase());
        let recovery_save = metadata
            .as_ref()
            .and_then(|value| value.recovery_save.as_ref());
        let candidate = HpcHeadlessJobCandidate {
            profile_id: profile.id.clone(),
            remote_job_id: job_id.clone(),
            job_name,
            task_kind: metadata
                .as_ref()
                .map(|value| value.task_kind.clone())
                .unwrap_or(task_kind),
            label,
            scheduler_state: hpc::slurm::normalize_scheduler_state(&state_name),
            remote_node: node,
            remote_workdir,
            submitted_at: metadata
                .as_ref()
                .map(|value| value.submitted_at.clone())
                .or(submitted_at),
            metadata_status: if metadata.is_some() {
                "metadata".to_string()
            } else {
                "legacy".to_string()
            },
            auto_save_available: recovery_save.is_some(),
            project_id: recovery_save.map(|value| value.project_id.clone()),
            cif_id: recovery_save.map(|value| value.cif_id.clone()),
        };
        candidates_by_key.insert(job_id, candidate);
    }

    let mut candidates: Vec<HpcHeadlessJobCandidate> = candidates_by_key.into_values().collect();
    candidates.sort_by(|a, b| {
        let rank = |state: &str| match hpc::slurm::normalize_scheduler_state(state).as_str() {
            "RUNNING" => 0,
            "PENDING" => 1,
            "COMPLETING" => 2,
            "COMPLETED" => 3,
            _ => 4,
        };
        rank(&a.scheduler_state)
            .cmp(&rank(&b.scheduler_state))
            .then_with(|| b.auto_save_available.cmp(&a.auto_save_available))
            .then_with(|| b.submitted_at.cmp(&a.submitted_at))
            .then_with(|| a.remote_job_id.cmp(&b.remote_job_id))
    });
    candidates.truncate(resolved_limit);
    Ok(candidates)
}

#[tauri::command]
async fn hpc_attach_headless_job(
    app: AppHandle,
    profile_id: Option<String>,
    remote_job_id: String,
    remote_workdir: Option<String>,
    state: State<'_, AppState>,
) -> Result<HpcAttachedJobResult, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;
    let workspace_root = resolve_remote_cleanup_path_with_timeout(
        &profile,
        secret.as_deref(),
        profile.remote_workspace_root.trim_end_matches('/'),
        8,
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let mut allowed_roots = vec![workspace_root];
    if let Ok(project_root) = resolve_remote_cleanup_path_with_timeout(
        &profile,
        secret.as_deref(),
        profile.remote_project_root.trim_end_matches('/'),
        5,
    )
    .await
    {
        allowed_roots.push(project_root.trim_end_matches('/').to_string());
    }
    let resolved_workdir = match remote_workdir {
        Some(path) if remote_path_is_under(&path, &allowed_roots) => {
            path.trim_end_matches('/').to_string()
        }
        _ => resolve_slurm_workdir(&profile, secret.as_deref(), &remote_job_id)
            .await
            .filter(|path| remote_path_is_under(path, &allowed_roots))
            .ok_or_else(|| {
                "Could not resolve a safe QCortado remote workdir for this job.".to_string()
            })?,
    };
    let metadata =
        read_recovery_metadata_for_dir(&profile, secret.as_deref(), &resolved_workdir).await;
    let mut task_kind = metadata.as_ref().map(|value| value.task_kind.clone());
    if task_kind.is_none() {
        let cmd = format!(
            "squeue -h -j {} -o \"%j\"",
            shell_single_quote_local(&remote_job_id)
        );
        task_kind = hpc::ssh::run_ssh_command_with_timeout(&profile, secret.as_deref(), &cmd, 8)
            .await
            .ok()
            .and_then(|name| infer_qcortado_task_kind(name.trim()));
    }
    let task_kind = task_kind.unwrap_or_else(|| "scf".to_string());
    let task_id = metadata
        .as_ref()
        .map(|value| value.task_id.clone())
        .filter(|value| !value.trim().is_empty())
        .unwrap_or_else(|| format!("recovered_{}", uuid::Uuid::new_v4()));
    let label = metadata
        .as_ref()
        .map(|value| value.label.clone())
        .filter(|value| !value.trim().is_empty())
        .unwrap_or_else(|| format!("Recovered {}", task_kind.to_ascii_uppercase()));
    let recovery_save_json = metadata
        .as_ref()
        .and_then(|value| value.recovery_save.as_ref())
        .and_then(|value| serde_json::to_value(value).ok());
    let resource_type = metadata
        .as_ref()
        .map(|value| value.resource_type.clone())
        .filter(|value| !value.trim().is_empty());
    let (scheduler_state, node) =
        scheduler_snapshot_for_job(&profile, secret.as_deref(), &remote_job_id).await;
    let local_sync_dir = projects::ensure_projects_dir(&app)?
        .join(".hpc_recovery")
        .join(&task_id);
    std::fs::create_dir_all(&local_sync_dir).map_err(|e| {
        format!(
            "Failed to create local recovery sync directory {}: {}",
            local_sync_dir.display(),
            e
        )
    })?;
    let run_sbatch_path = format!("{}/run.sbatch", resolved_workdir.trim_end_matches('/'));
    let slurm_script =
        read_remote_text_optional(&profile, secret.as_deref(), &run_sbatch_path).await;
    let cancel_flag = state
        .process_manager
        .attach_hpc_task(
            task_id.clone(),
            task_kind.clone(),
            label.clone(),
            metadata.as_ref().map(|value| value.submitted_at.clone()),
            Some(profile.id.clone()),
            resource_type,
            Some(remote_job_id.clone()),
            Some(scheduler_state),
            node,
            Some(resolved_workdir.clone()),
            None,
            Some(local_sync_dir.to_string_lossy().to_string()),
            recovery_save_json.clone(),
        )
        .await;

    let pm = state.process_manager.clone();
    let app_handle = app.clone();
    let attach_task_id = task_id.clone();
    let attach_kind = task_kind.clone();
    let attach_workdir = resolved_workdir.clone();
    let attach_job_id = remote_job_id.clone();
    let attach_preview = format!("reattach {}", remote_job_id);
    tokio::spawn(async move {
        let result = hpc::runner::run_attached_batch_task(
            app_handle.clone(),
            pm.clone(),
            hpc::runner::HpcAttachRequest {
                task_id: attach_task_id.clone(),
                task_kind: attach_kind.clone(),
                task_label: label,
                profile,
                secret,
                remote_job_id: attach_job_id,
                remote_workdir: attach_workdir,
                remote_project_path: None,
                slurm_script,
                sbatch_preview: Some(attach_preview),
                local_sync_dir: local_sync_dir.clone(),
                cancel_flag,
            },
        )
        .await;
        match result {
            Ok(_) => {
                let json = parse_recovered_hpc_result(&attach_kind, &local_sync_dir)
                    .unwrap_or_else(
                        |err| serde_json::json!({ "recovered": true, "parse_error": err }),
                    );
                pm.complete(&attach_task_id, json).await;
                let _ = app_handle.emit(&format!("task-complete:{}", attach_task_id), "completed");
            }
            Err(err) => {
                pm.fail(&attach_task_id, err.clone()).await;
                let _ = app_handle.emit(
                    &format!("task-status:{}", attach_task_id),
                    &format!("failed:{}", err),
                );
            }
        }
    });

    Ok(HpcAttachedJobResult {
        task_id,
        task_kind,
        remote_job_id,
        remote_workdir: resolved_workdir,
        auto_save_available: recovery_save_json.is_some(),
    })
}

fn parse_recovered_hpc_result(
    task_kind: &str,
    work_path: &Path,
) -> Result<serde_json::Value, String> {
    let kind = task_kind.trim().to_ascii_lowercase();
    match kind.as_str() {
        "scf" => {
            let output = std::fs::read_to_string(work_path.join("pw.out"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.out")))
                .map_err(|e| format!("Failed to read recovered SCF output: {}", e))?;
            serde_json::to_value(parse_pw_output(&output)).map_err(|e| e.to_string())
        }
        "phonon" => {
            let ph_output = std::fs::read_to_string(work_path.join("ph.out"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.out")))
                .unwrap_or_default();
            let (converged, n_qpoints) = parse_ph_output(&ph_output);
            let dos_data = read_phonon_dos_file(&work_path.join("phonon_dos")).ok();
            let dispersion_data = read_phonon_dispersion_file(&work_path.join("phonon_freq.gp"))
                .or_else(|_| read_phonon_dispersion_file(&work_path.join("phonon_freq")))
                .ok();
            serde_json::to_value(PhononResult {
                converged,
                n_qpoints,
                n_modes: dispersion_data
                    .as_ref()
                    .map(|data| data.n_modes)
                    .unwrap_or(0),
                dos_data,
                dispersion_data,
                raw_output: ph_output,
            })
            .map_err(|e| e.to_string())
        }
        "epw" => {
            let stdout = std::fs::read_to_string(work_path.join("epw.out"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.out")))
                .unwrap_or_default();
            let stderr = std::fs::read_to_string(work_path.join("epw.err"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.err")))
                .unwrap_or_default();
            let artifacts = collect_epw_artifacts(work_path);
            let parsed = parse_epw_result_v2(
                &format!("{}\n{}", stdout, stderr),
                work_path,
                artifacts.clone(),
                true,
            );
            serde_json::to_value(serde_json::json!({
                "schema_version": EPW_SCHEMA_VERSION,
                "result_summary": parsed.summary,
                "artifacts": artifacts,
                "recovered": true
            }))
            .map_err(|e| e.to_string())
        }
        "bands" => {
            let bands_out_text = std::fs::read_to_string(work_path.join("bands.out"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.out")))
                .unwrap_or_default();
            let fermi_energy = extract_fermi_energy_from_text(&bands_out_text).unwrap_or(0.0);
            let gnu_file = std::fs::read_dir(work_path)
                .map_err(|e| format!("Failed to scan recovered bands directory: {}", e))?
                .filter_map(Result::ok)
                .map(|entry| entry.path())
                .find(|path| {
                    path.extension()
                        .and_then(|ext| ext.to_str())
                        .map(|ext| ext.eq_ignore_ascii_case("gnu"))
                        .unwrap_or(false)
                })
                .ok_or_else(|| "No recovered bands .gnu file found.".to_string())?;
            let band_data = read_bands_gnu_file(&gnu_file, fermi_energy)
                .map_err(|e| format!("Failed to parse recovered bands data: {}", e))?;
            serde_json::to_value(band_data).map_err(|e| e.to_string())
        }
        "dos" => {
            let dos_file = work_path.join("dos.dat");
            let dos_content = std::fs::read_to_string(&dos_file).map_err(|e| {
                format!(
                    "Failed to read recovered DOS file {}: {}",
                    dos_file.display(),
                    e
                )
            })?;
            let (energies, dos_values) = parse_dos_file(&dos_content).ok_or_else(|| {
                format!(
                    "Failed to parse recovered DOS data from {}",
                    dos_file.display()
                )
            })?;
            let nscf_output = std::fs::read_to_string(work_path.join("nscf.out"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.out")))
                .unwrap_or_default();
            let energy_range = [
                energies.iter().cloned().fold(f64::INFINITY, f64::min),
                energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            ];
            let max_dos = dos_values.iter().cloned().fold(0.0_f64, f64::max);
            serde_json::to_value(ElectronicDosData {
                points: energies.len(),
                energies,
                dos: dos_values,
                fermi_energy: extract_fermi_energy_from_text(&nscf_output),
                energy_range,
                max_dos,
            })
            .map_err(|e| e.to_string())
        }
        "fermi_surface" => {
            let fermi_output = std::fs::read_to_string(work_path.join("fermi_velocity.out"))
                .or_else(|_| std::fs::read_to_string(work_path.join("slurm.out")))
                .unwrap_or_default();
            let frmsf_files = collect_frmsf_files(work_path)?;
            let primary_file = frmsf_files
                .iter()
                .find(|file| {
                    let lower = file.file_name.to_ascii_lowercase();
                    lower.ends_with("/vfermi.frmsf")
                        || lower == "vfermi.frmsf"
                        || lower.ends_with("_vfermi.frmsf")
                })
                .map(|file| file.file_name.clone())
                .unwrap_or_else(|| {
                    frmsf_files
                        .first()
                        .map(|file| file.file_name.clone())
                        .unwrap_or_default()
                });
            serde_json::to_value(FermiSurfaceData {
                k_grid: [0, 0, 0],
                fermi_energy: extract_fermi_energy_from_text(&fermi_output),
                primary_file,
                frmsf_files,
            })
            .map_err(|e| e.to_string())
        }
        _ => Ok(serde_json::json!({ "recovered": true, "task_kind": kind })),
    }
}

/// Opens a read-only HPC activity popout window.
#[tauri::command]
fn hpc_open_activity_window(app: AppHandle) -> Result<(), String> {
    let existing = app.get_webview_window("hpc-activity");
    if let Some(window) = existing {
        window
            .set_focus()
            .map_err(|e| format!("Failed to focus activity window: {}", e))?;
        return Ok(());
    }

    tauri::WebviewWindowBuilder::new(
        &app,
        "hpc-activity",
        tauri::WebviewUrl::App("index.html?hpc_activity=1".into()),
    )
    .title("QCortado Cluster Activity")
    .inner_size(760.0, 520.0)
    .resizable(true)
    .build()
    .map_err(|e| format!("Failed to open activity window: {}", e))?;
    Ok(())
}

/// Downloads artifacts for an existing HPC task into its local sync directory.
#[tauri::command]
async fn hpc_download_task_artifacts(
    app: AppHandle,
    task_id: String,
    full: Option<bool>,
    state: State<'_, AppState>,
) -> Result<hpc::runner::HpcArtifactSyncReport, String> {
    let context = state
        .process_manager
        .get_hpc_transfer_context(&task_id)
        .await
        .ok_or_else(|| format!("Task not found: {}", task_id))?;

    if context.backend.as_deref() != Some("hpc") {
        return Err("Artifact download is only available for HPC tasks.".to_string());
    }
    if full.unwrap_or(true) && matches!(context.status, process_manager::TaskStatus::Running) {
        return Err("Task is still running. Wait for completion before full download.".to_string());
    }

    let remote_workdir = if matches!(context.status, process_manager::TaskStatus::Running) {
        context.remote_workdir.or(context.remote_project_path)
    } else {
        context.remote_project_path.or(context.remote_workdir)
    }
    .ok_or_else(|| "Task does not have a remote artifact path recorded.".to_string())?;
    let local_sync_dir = context
        .local_sync_dir
        .ok_or_else(|| "Task does not have a local sync directory recorded.".to_string())?;
    let profile = resolve_hpc_profile_from_state(&state, context.hpc_profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;
    let request_full = full.unwrap_or(true);
    let mode = if request_full {
        hpc::runner::ArtifactSyncMode::Full
    } else if context.task_type.trim().eq_ignore_ascii_case("epw") {
        hpc::runner::ArtifactSyncMode::EpwReady
    } else {
        hpc::runner::ArtifactSyncMode::Minimal
    };

    let collect_line = format!("HPC_STAGE|Collecting|{} ({})", remote_workdir, mode.label());
    let _ = app.emit(&format!("task-output:{}", task_id), &collect_line);
    state
        .process_manager
        .append_output(&task_id, collect_line)
        .await;

    let report = hpc::runner::sync_remote_artifacts(
        &app,
        &state.process_manager,
        hpc::runner::HpcArtifactSyncRequest {
            task_id: context.task_id,
            task_kind: context.task_type,
            profile,
            secret,
            remote_workdir,
            local_sync_dir: PathBuf::from(local_sync_dir),
            mode,
        },
    )
    .await?;
    let remote_storage_bytes = report.downloaded_bytes.saturating_add(report.skipped_bytes);
    state
        .process_manager
        .set_remote_storage_bytes(&task_id, Some(remote_storage_bytes))
        .await;

    let saved_line = format!(
        "HPC_STAGE|Saved|{} sync complete ({} files, {:.2} MB downloaded, {} skipped, remote {:.2} MB)",
        report.mode,
        report.downloaded_files,
        report.downloaded_bytes as f64 / (1024.0 * 1024.0),
        report.skipped_files,
        remote_storage_bytes as f64 / (1024.0 * 1024.0),
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &saved_line);
    state
        .process_manager
        .append_output(&task_id, saved_line)
        .await;

    Ok(report)
}

/// Downloads artifacts for a saved HPC calculation entry.
#[tauri::command]
async fn hpc_download_calculation_artifacts(
    app: AppHandle,
    project_id: String,
    calc_id: String,
    profile_id: Option<String>,
    full: Option<bool>,
    state: State<'_, AppState>,
) -> Result<hpc::runner::HpcArtifactSyncReport, String> {
    let project = projects::get_project(app.clone(), project_id.clone())?;
    let calculation = project
        .cif_variants
        .iter()
        .flat_map(|variant| variant.calculations.iter())
        .find(|calc| calc.id == calc_id)
        .ok_or_else(|| format!("Calculation not found in project: {}", calc_id))?;

    let parameters = calculation
        .parameters
        .as_object()
        .ok_or_else(|| "Calculation parameters are unavailable.".to_string())?;

    let backend = parameters
        .get("execution_backend")
        .and_then(|value| value.as_str())
        .unwrap_or("")
        .trim()
        .to_ascii_lowercase();
    if backend != "hpc"
        && parameters
            .get("remote_workdir")
            .and_then(|value| value.as_str())
            .map(|value| value.trim().is_empty())
            .unwrap_or(true)
    {
        return Err("This calculation does not have HPC artifact metadata.".to_string());
    }

    let remote_workdir = parameters
        .get("remote_project_path")
        .or_else(|| parameters.get("remote_workdir"))
        .and_then(|value| value.as_str())
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
        .ok_or_else(|| "No remote working directory is stored for this calculation.".to_string())?;

    let saved_profile_id = parameters
        .get("hpc_profile_id")
        .and_then(|value| value.as_str())
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty());
    let resolved_profile_id = profile_id.or(saved_profile_id);
    let profile = resolve_hpc_profile_from_state(&state, resolved_profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let local_sync_dir = projects::get_projects_dir(&app)?
        .join(&project_id)
        .join("calculations")
        .join(&calc_id)
        .join("tmp");
    std::fs::create_dir_all(&local_sync_dir).map_err(|e| {
        format!(
            "Failed to create local calculation directory {}: {}",
            local_sync_dir.display(),
            e
        )
    })?;

    let request_full = full.unwrap_or(true);
    let mode = if request_full {
        hpc::runner::ArtifactSyncMode::Full
    } else if calculation.calc_type.trim().eq_ignore_ascii_case("epw") {
        hpc::runner::ArtifactSyncMode::EpwReady
    } else {
        hpc::runner::ArtifactSyncMode::Minimal
    };

    let connecting_line = format!("HPC_STAGE|Connecting|{}", profile.host);
    let _ = app.emit(&format!("task-output:{}", calc_id), &connecting_line);
    state
        .process_manager
        .append_output(&calc_id, connecting_line)
        .await;

    let collecting_line = format!("HPC_STAGE|Collecting|{} ({})", remote_workdir, mode.label());
    let _ = app.emit(&format!("task-output:{}", calc_id), &collecting_line);
    state
        .process_manager
        .append_output(&calc_id, collecting_line)
        .await;

    let report = hpc::runner::sync_remote_artifacts(
        &app,
        &state.process_manager,
        hpc::runner::HpcArtifactSyncRequest {
            task_id: calc_id.clone(),
            task_kind: calculation.calc_type.clone(),
            profile,
            secret,
            remote_workdir,
            local_sync_dir,
            mode,
        },
    )
    .await?;

    let remote_storage_bytes = report.downloaded_bytes.saturating_add(report.skipped_bytes);
    projects::refresh_calculation_artifact_metadata(
        &app,
        &project_id,
        &calc_id,
        Some(remote_storage_bytes),
        Some(report.mode.as_str()),
    )?;

    let saved_line = format!(
        "HPC_STAGE|Saved|{} sync complete ({} files, {:.2} MB downloaded, {} skipped, remote {:.2} MB)",
        report.mode,
        report.downloaded_files,
        report.downloaded_bytes as f64 / (1024.0 * 1024.0),
        report.skipped_files,
        remote_storage_bytes as f64 / (1024.0 * 1024.0),
    );
    let _ = app.emit(&format!("task-output:{}", calc_id), &saved_line);
    state
        .process_manager
        .append_output(&calc_id, saved_line)
        .await;

    Ok(report)
}

/// Lists finished remote phonon directories that look recoverable.
#[tauri::command]
async fn hpc_list_recoverable_remote_phonon_runs(
    profile_id: Option<String>,
    limit: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Vec<HpcRecoverableRemotePhononRun>, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let resolved_limit = limit.unwrap_or(30).clamp(1, 200) as usize;
    let runs = list_recoverable_remote_phonon_runs_for_profile(
        &profile,
        secret.as_deref(),
        resolved_limit,
    )
    .await?;
    eprintln!(
        "[remote-phonon-recovery] list profile={} limit={} recoverable_runs={}",
        profile.id,
        resolved_limit,
        runs.len()
    );
    for run in runs.iter().take(8) {
        eprintln!(
            "[remote-phonon-recovery] candidate location={} modified_at_epoch={} remote_workdir={}",
            run.location, run.modified_at_epoch, run.remote_workdir
        );
    }
    Ok(runs)
}

/// Returns a detailed probe report for remote phonon recovery diagnostics.
#[tauri::command]
async fn hpc_debug_remote_phonon_recovery(
    profile_id: Option<String>,
    state: State<'_, AppState>,
) -> Result<HpcRemotePhononRecoveryDebugReport, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let workspace_root = resolve_remote_cleanup_path(
        &profile,
        secret.as_deref(),
        profile.remote_workspace_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let project_root = resolve_remote_cleanup_path(
        &profile,
        secret.as_deref(),
        profile.remote_project_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let project_phonon_root = format!("{}/phonon", project_root);

    let workspace_probe_cmd =
        build_remote_phonon_debug_probe_scan_cmd(&workspace_root, Some("qcortado_hpc_bundle_*"));
    let project_probe_cmd = build_remote_phonon_debug_probe_scan_cmd(&project_phonon_root, None);
    let workspace_probe_output =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &workspace_probe_cmd).await?;
    let project_probe_output =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &project_probe_cmd).await?;

    let recoverable_runs =
        list_recoverable_remote_phonon_runs_for_profile(&profile, secret.as_deref(), 200).await?;
    let workspace_probe_rows = workspace_probe_output
        .lines()
        .filter(|line| !line.trim().is_empty())
        .count();
    let project_probe_rows = project_probe_output
        .lines()
        .filter(|line| !line.trim().is_empty())
        .count();
    eprintln!(
        "[remote-phonon-recovery] debug profile={} workspace_root={} project_phonon_root={} workspace_probe_rows={} project_probe_rows={} recoverable_runs={}",
        profile.id,
        workspace_root,
        project_phonon_root,
        workspace_probe_rows,
        project_probe_rows,
        recoverable_runs.len()
    );

    Ok(HpcRemotePhononRecoveryDebugReport {
        profile_id: profile.id.clone(),
        workspace_root,
        project_phonon_root,
        workspace_probe_output,
        project_probe_output,
        recoverable_runs,
    })
}

/// Downloads a remote phonon run and imports it into a project as a recovered calculation.
#[tauri::command]
async fn hpc_recover_remote_phonon_calculation(
    app: AppHandle,
    project_id: String,
    cif_id: String,
    remote_workdir: String,
    profile_id: Option<String>,
    source_scf_id: Option<String>,
    q_grid: Option<[u32; 3]>,
    tr2_ph: Option<f64>,
    q_path: Option<String>,
    state: State<'_, AppState>,
) -> Result<projects::CalculationRun, String> {
    let normalized_remote_workdir = remote_workdir.trim().trim_end_matches('/').to_string();
    if normalized_remote_workdir.is_empty() {
        return Err("Remote working directory is required.".to_string());
    }

    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;
    eprintln!(
        "[remote-phonon-recovery] start profile={} project_id={} cif_id={} remote_workdir={}",
        profile.id, project_id, cif_id, normalized_remote_workdir
    );

    let recovery_id = uuid::Uuid::new_v4().to_string();
    let local_sync_dir =
        std::env::temp_dir().join(format!("qcortado_remote_phonon_recovery_{}", recovery_id));
    std::fs::create_dir_all(&local_sync_dir).map_err(|e| {
        format!(
            "Failed to create temporary recovery directory {}: {}",
            local_sync_dir.display(),
            e
        )
    })?;

    let sync_result = hpc::runner::sync_remote_artifacts(
        &app,
        &state.process_manager,
        hpc::runner::HpcArtifactSyncRequest {
            task_id: format!("remote_phonon_recovery_{}", recovery_id),
            task_kind: "phonon".to_string(),
            profile: profile.clone(),
            secret,
            remote_workdir: normalized_remote_workdir.clone(),
            local_sync_dir: local_sync_dir.clone(),
            mode: hpc::runner::ArtifactSyncMode::Minimal,
        },
    )
    .await;

    let recovery_result = match sync_result {
        Ok(report) => {
            eprintln!(
                "[remote-phonon-recovery] synced profile={} remote_workdir={} downloaded_files={} downloaded_bytes={} skipped_files={} skipped_bytes={} mode={}",
                profile.id,
                normalized_remote_workdir,
                report.downloaded_files,
                report.downloaded_bytes,
                report.skipped_files,
                report.skipped_bytes,
                report.mode
            );
            let remote_storage_bytes = report.downloaded_bytes.saturating_add(report.skipped_bytes);
            let recovered = projects::recover_phonon_calculation(
                app.clone(),
                project_id.clone(),
                cif_id,
                local_sync_dir.to_string_lossy().to_string(),
                source_scf_id,
                q_grid,
                tr2_ph,
                q_path,
                Some("hpc".to_string()),
                Some(profile.id.clone()),
                Some(normalized_remote_workdir.clone()),
                None,
                Some(remote_storage_bytes),
                Some(report.mode.clone()),
            )?;
            let _ = projects::refresh_calculation_artifact_metadata(
                &app,
                &project_id,
                &recovered.id,
                Some(remote_storage_bytes),
                Some(report.mode.as_str()),
            );
            eprintln!(
                "[remote-phonon-recovery] imported calc_id={} project_id={} remote_workdir={}",
                recovered.id, project_id, normalized_remote_workdir
            );
            Ok(recovered)
        }
        Err(err) => {
            eprintln!(
                "[remote-phonon-recovery] sync failed profile={} remote_workdir={} error={}",
                profile.id, normalized_remote_workdir, err
            );
            Err(err)
        }
    };

    let _ = std::fs::remove_dir_all(&local_sync_dir);
    recovery_result
}

/// Removes remote QCortado HPC directories that are no longer referenced by local projects/tasks.
#[tauri::command]
async fn hpc_clean_remote_orphans(
    app: AppHandle,
    profile_id: Option<String>,
    state: State<'_, AppState>,
) -> Result<HpcRemoteOrphanCleanupResult, String> {
    let profile = resolve_hpc_profile_from_state(&state, profile_id)?;
    let secret = hpc::credentials::resolve_secret(
        &profile.id,
        &profile.username,
        &profile.host,
        profile.credential_persisted,
    )?;

    let workspace_root = resolve_remote_cleanup_path(
        &profile,
        secret.as_deref(),
        profile.remote_workspace_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();
    let project_root = resolve_remote_cleanup_path(
        &profile,
        secret.as_deref(),
        profile.remote_project_root.trim_end_matches('/'),
    )
    .await?
    .trim_end_matches('/')
    .to_string();

    let mut referenced_paths = collect_saved_hpc_remote_paths(&app)?;
    for task in state.process_manager.list_tasks().await {
        if task.backend.as_deref() != Some("hpc") {
            continue;
        }
        if let Some(path) = task
            .remote_workdir
            .as_deref()
            .and_then(normalize_remote_path_for_tracking)
        {
            referenced_paths.insert(path);
        }
        if let Some(path) = task
            .remote_project_path
            .as_deref()
            .and_then(normalize_remote_path_for_tracking)
        {
            referenced_paths.insert(path);
        }
    }

    let mut candidate_paths: HashSet<String> = HashSet::new();
    let workspace_scan_cmd = format!(
        "if [ -d {root} ]; then find {root} -mindepth 1 -maxdepth 1 -type d -name 'qcortado_hpc_bundle_*' -print; fi",
        root = shell_single_quote_local(&workspace_root)
    );
    let workspace_listing =
        hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &workspace_scan_cmd).await?;
    for line in workspace_listing.lines() {
        if let Some(path) = normalize_remote_path_for_tracking(line) {
            candidate_paths.insert(path);
        }
    }

    for task_kind in HPC_REMOTE_PROJECT_TASK_KINDS {
        let project_scan_cmd = format!(
            "if [ -d {root}/{kind} ]; then find {root}/{kind} -mindepth 1 -maxdepth 1 -type d -print; fi",
            root = shell_single_quote_local(&project_root),
            kind = shell_single_quote_local(task_kind),
        );
        let listing =
            hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &project_scan_cmd).await?;
        for line in listing.lines() {
            if let Some(path) = normalize_remote_path_for_tracking(line) {
                candidate_paths.insert(path);
            }
        }
    }

    let mut orphan_paths: Vec<String> = candidate_paths
        .iter()
        .filter(|path| !referenced_paths.contains(*path))
        .cloned()
        .collect();
    orphan_paths.sort();

    let mut removed_paths: Vec<String> = Vec::new();
    let mut failed_paths: Vec<String> = Vec::new();
    for path in &orphan_paths {
        let remove_cmd = format!(
            "if [ -e {target} ]; then rm -rf -- {target}; fi",
            target = shell_single_quote_local(path)
        );
        match hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &remove_cmd).await {
            Ok(_) => removed_paths.push(path.clone()),
            Err(err) => {
                eprintln!("Failed to remove orphan remote path {}: {}", path, err);
                failed_paths.push(path.clone());
            }
        }
    }

    Ok(HpcRemoteOrphanCleanupResult {
        profile_id: profile.id,
        scanned_paths: candidate_paths.len(),
        referenced_paths: referenced_paths.len(),
        orphan_paths: orphan_paths.len(),
        removed_paths,
        failed_paths,
    })
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
        "fermi_velocity.x",
        "fermi_proj.x",
        "projwfc.x",
        "pp.x",
        "ph.x",
        "epw.x",
        "pw2wannier90.x",
        "dynmat.x",
        "plotband.x",
        "neb.x",
        "hp.x",
        "turbo_lanczos.x",
        "xspectra.x",
    ];

    let available: Vec<String> = executables
        .iter()
        .filter(|exe| {
            if **exe == "epw.x" {
                resolve_local_epw_path(bin_dir.as_path()).is_some()
            } else {
                bin_dir.join(exe).exists()
            }
        })
        .map(|s| s.to_string())
        .collect();

    Ok(available)
}

/// Analyzes symmetry and returns standardized conventional/primitive cells plus
/// reciprocal-coordinate transform matrices between them.
#[tauri::command]
fn analyze_structure_symmetry(
    input: symmetry::SymmetryAnalyzeInput,
) -> Result<symmetry::SymmetryAnalyzeResult, String> {
    symmetry::analyze_structure(input)
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

    prepare_working_directory(&work_path, false)?;

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

    prepare_working_directory(&work_path, false)?;

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
                .args(["-pd", ".true."])
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

/// Metadata extracted from a pseudopotential header.
#[derive(serde::Serialize, Clone)]
pub struct PseudopotentialMetadata {
    pub filename: String,
    pub supports_soc: bool,
    pub pseudo_type: Option<String>,
    pub relativistic: Option<String>,
    pub cutoff_wfc: Option<f64>,
    pub cutoff_rho: Option<f64>,
    pub cutoff_wfc_source: Option<String>,
    pub cutoff_rho_source: Option<String>,
    #[serde(default)]
    pub available_angular_momenta: Vec<u8>,
    pub available_angular_momenta_source: Option<String>,
    pub max_angular_momentum: Option<u8>,
}

#[derive(serde::Deserialize)]
struct DjrepoCutoffHint {
    ecut: Option<f64>,
}

#[derive(serde::Deserialize)]
struct DjrepoHints {
    low: Option<DjrepoCutoffHint>,
    normal: Option<DjrepoCutoffHint>,
    high: Option<DjrepoCutoffHint>,
}

#[derive(serde::Deserialize)]
struct DjrepoMetadata {
    hints: Option<DjrepoHints>,
}

fn parse_upf_attr_map(attrs: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    let attr_re =
        Regex::new(r#"([A-Za-z_][A-Za-z0-9_-]*)\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>]+))"#)
            .expect("valid UPF attribute regex");

    for cap in attr_re.captures_iter(attrs) {
        let key = cap.get(1).map(|m| m.as_str().to_string());
        let value = cap
            .get(2)
            .or_else(|| cap.get(3))
            .or_else(|| cap.get(4))
            .map(|m| m.as_str().trim().to_string());
        if let (Some(key), Some(value)) = (key, value) {
            map.insert(key, value);
        }
    }

    map
}

fn parse_upf_bool(value: Option<&String>) -> bool {
    matches!(
        value.map(|v| v.trim().to_lowercase()),
        Some(ref v) if v == "true" || v == ".true." || v == "t" || v == "1"
    )
}

fn parse_upf_f64(value: Option<&String>) -> Option<f64> {
    value
        .and_then(|raw| raw.trim().parse::<f64>().ok())
        .filter(|value| *value > 0.0)
}

fn capture_xml_value(text: &str, tag: &str) -> Option<String> {
    let pattern = format!(
        r"(?is)<{}\b[^>]*>\s*([^<]+?)\s*</{}>",
        regex::escape(tag),
        regex::escape(tag)
    );
    capture_group(text, &pattern).map(|value| value.trim().to_string())
}

fn capture_group(text: &str, pattern: &str) -> Option<String> {
    let re = Regex::new(pattern).ok()?;
    re.captures(text)
        .and_then(|caps| caps.get(1))
        .map(|m| m.as_str().to_string())
}

fn parse_upf_cutoff(text: &str, label: &str) -> Option<f64> {
    let pattern = format!(
        r"(?im){}\s*:\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)",
        regex::escape(label),
    );
    capture_group(text, &pattern)
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| *value > 0.0)
}

fn insert_valid_angular_channel(channels: &mut std::collections::BTreeSet<u8>, raw: &str) {
    if let Ok(value) = raw.trim().parse::<i32>() {
        if (0..=6).contains(&value) {
            channels.insert(value as u8);
        }
    }
}

#[derive(Clone, Copy)]
struct AngularMomentumParseDetails {
    source: Option<&'static str>,
    max_angular_momentum: Option<u8>,
}

fn parse_upf_max_angular_momentum(content: &str) -> Option<u8> {
    capture_group(
        content,
        r#"(?is)<PP_HEADER\b[^>]*\bl_max\s*=\s*["']?([0-9]+)"#,
    )
    .or_else(|| {
        capture_group(
            content,
            r"(?im)^\s*([0-9]+)\s+Max angular momentum component\s*$",
        )
    })
    .and_then(|value| value.trim().parse::<u8>().ok())
    .filter(|value| *value <= 6)
}

fn parse_upf_angular_momentum_details(
    content: &str,
    info_text: &str,
) -> (Vec<u8>, AngularMomentumParseDetails) {
    let mut channels = std::collections::BTreeSet::new();
    let mut source: Option<&'static str> = None;
    let max_angular_momentum = parse_upf_max_angular_momentum(content);

    if let Ok(tag_re) = Regex::new(r#"(?is)<PP_(?:BETA|CHI|PSWFC|AEWFC)\b([^>]*)>"#) {
        for caps in tag_re.captures_iter(content) {
            let attrs = caps
                .get(1)
                .map(|m| parse_upf_attr_map(m.as_str()))
                .unwrap_or_default();
            if let Some(value) = attrs
                .get("angular_momentum")
                .or_else(|| attrs.get("l"))
                .or_else(|| attrs.get("lll"))
            {
                insert_valid_angular_channel(&mut channels, value);
            }
        }
        if !channels.is_empty() {
            source = Some("upf_tags");
        }
    }

    if channels.is_empty() {
        if let Ok(info_re) = Regex::new(r"(?im)\bl(?:\(\d+\))?\s*=\s*([0-9]+)") {
            for caps in info_re.captures_iter(info_text) {
                if let Some(value) = caps.get(1) {
                    insert_valid_angular_channel(&mut channels, value.as_str());
                }
            }
        }
        if !channels.is_empty() {
            source = Some("upf_info_l_lines");
        }
    }

    if channels.is_empty() {
        if let Ok(orbital_row_re) = Regex::new(
            r"(?i)^\s*(\d+[spdfghi])(?:\s+(\d+))?(?:\s+(\d+))?(?:\s+[-+0-9.eed]+){1,}.*$",
        ) {
            for line in info_text.lines().chain(content.lines().take(240)) {
                if let Some(caps) = orbital_row_re.captures(line) {
                    if let Some(third) = caps.get(3) {
                        insert_valid_angular_channel(&mut channels, third.as_str());
                        continue;
                    }
                    if let Some(second) = caps.get(2) {
                        insert_valid_angular_channel(&mut channels, second.as_str());
                    }
                }
            }
        }
        if !channels.is_empty() {
            source = Some("upf_info_orbital_rows");
        }
    }

    if channels.is_empty() {
        if let Some(max_l) = max_angular_momentum {
            for channel in 0..=max_l {
                channels.insert(channel);
            }
            if !channels.is_empty() {
                source = Some("upf_l_max_fallback");
            }
        }
    }

    (
        channels.into_iter().collect(),
        AngularMomentumParseDetails {
            source,
            max_angular_momentum,
        },
    )
}

#[cfg(test)]
fn parse_upf_available_angular_momenta(content: &str, info_text: &str) -> Vec<u8> {
    parse_upf_angular_momentum_details(content, info_text).0
}

fn parse_djrepo_wavefunction_cutoff_ry(text: &str) -> Option<f64> {
    let metadata: DjrepoMetadata = serde_json::from_str(text).ok()?;
    let hints = metadata.hints?;
    let ecut_ha = hints
        .normal
        .as_ref()
        .and_then(|hint| hint.ecut)
        .or_else(|| hints.high.as_ref().and_then(|hint| hint.ecut))
        .or_else(|| hints.low.as_ref().and_then(|hint| hint.ecut))
        .filter(|value| *value > 0.0)?;
    Some(ecut_ha * 2.0)
}

fn is_hamann_oncv_lone_rho_wavefunction_hint(
    generated: Option<&str>,
    pseudo_type: Option<&str>,
    has_wfc: bool,
    cutoff_wfc: Option<f64>,
    cutoff_rho: Option<f64>,
    info_text: &str,
) -> bool {
    let generated_lower = generated.unwrap_or_default().to_lowercase();
    let info_lower = info_text.to_lowercase();
    let pseudo_type = pseudo_type.unwrap_or_default();

    cutoff_wfc.is_none()
        && cutoff_rho.is_some()
        && !has_wfc
        && pseudo_type.eq_ignore_ascii_case("nc")
        && (generated_lower.contains("oncvpsp code by d. r. hamann")
            || info_lower.contains("oncvpsp")
            || info_lower.contains("d. r. hamann"))
}

fn parse_pseudopotential_metadata_from_content(
    filename: String,
    content: &str,
) -> PseudopotentialMetadata {
    let header_block =
        capture_group(content, r"(?is)<PP_HEADER\b[^>]*>(.*?)</PP_HEADER>").unwrap_or_default();
    let header_attrs = capture_group(content, r"(?is)<PP_HEADER\b([^>]*)>")
        .map(|attrs| parse_upf_attr_map(&attrs))
        .unwrap_or_default();
    let info_text = capture_group(content, r"(?is)<PP_INFO>(.*?)</PP_INFO>")
        .unwrap_or_else(|| content.lines().take(200).collect::<Vec<_>>().join("\n"));
    let info_lower = info_text.to_lowercase();
    let lower_content = content.to_lowercase();
    let header_has_so_xml = capture_xml_value(&header_block, "has_so");
    let header_has_wfc_xml = capture_xml_value(&header_block, "has_wfc");
    let generated = header_attrs
        .get("generated")
        .cloned()
        .or_else(|| capture_xml_value(&header_block, "generated"));

    let pseudo_type = header_attrs
        .get("pseudo_type")
        .cloned()
        .or_else(|| capture_xml_value(&header_block, "pseudo_type"))
        .or_else(|| capture_group(&info_text, r"(?im)Pseudopotential type:\s*([^\r\n<]+)"))
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty());

    let relativistic = header_attrs
        .get("relativistic")
        .cloned()
        .or_else(|| capture_xml_value(&header_block, "relativistic"))
        .or_else(|| {
            capture_group(
                &info_text,
                r"(?im)((?:non-|scalar-|fully-)?relativistic pseudopotential)",
            )
        });
    let (available_angular_momenta, angular_momentum_details) =
        parse_upf_angular_momentum_details(content, &info_text);

    let has_so = parse_upf_bool(header_attrs.get("has_so"))
        || parse_upf_bool(header_has_so_xml.as_ref())
        || lower_content.contains("<pp_spin_orb")
        || info_lower.contains("fully-relativistic")
        || info_lower.contains("spin-orbit");
    let has_wfc =
        parse_upf_bool(header_attrs.get("has_wfc")) || parse_upf_bool(header_has_wfc_xml.as_ref());

    let (raw_cutoff_wfc, raw_cutoff_wfc_source) =
        if let Some(value) = parse_upf_f64(header_attrs.get("wfc_cutoff")) {
            (Some(value), Some("upf"))
        } else if let Some(value) = capture_xml_value(&header_block, "wfc_cutoff")
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| *value > 0.0)
        {
            (Some(value), Some("upf"))
        } else if let Some(value) = parse_upf_f64(header_attrs.get("ecutwfc")) {
            (Some(value), Some("upf"))
        } else if let Some(value) =
            parse_upf_cutoff(&info_text, "Suggested minimum cutoff for wavefunctions")
        {
            (Some(value), Some("upf_info"))
        } else {
            (None, None)
        };
    let (raw_cutoff_rho, raw_cutoff_rho_source) =
        if let Some(value) = parse_upf_f64(header_attrs.get("rho_cutoff")) {
            (Some(value), Some("upf"))
        } else if let Some(value) = capture_xml_value(&header_block, "rho_cutoff")
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| *value > 0.0)
        {
            (Some(value), Some("upf"))
        } else if let Some(value) = parse_upf_f64(header_attrs.get("ecutrho")) {
            (Some(value), Some("upf"))
        } else if let Some(value) =
            parse_upf_cutoff(&info_text, "Suggested minimum cutoff for charge density")
        {
            (Some(value), Some("upf_info"))
        } else {
            (None, None)
        };
    let lone_rho_is_hartree_wfc = is_hamann_oncv_lone_rho_wavefunction_hint(
        generated.as_deref(),
        pseudo_type.as_deref(),
        has_wfc,
        raw_cutoff_wfc,
        raw_cutoff_rho,
        &info_text,
    );
    let cutoff_wfc = if lone_rho_is_hartree_wfc {
        raw_cutoff_rho.map(|value| value * 2.0)
    } else {
        raw_cutoff_wfc
    };
    let cutoff_wfc_source = if lone_rho_is_hartree_wfc {
        raw_cutoff_rho.map(|_| "upf_fallback".to_string())
    } else {
        raw_cutoff_wfc_source.map(str::to_string)
    };
    let cutoff_rho = if lone_rho_is_hartree_wfc {
        None
    } else {
        raw_cutoff_rho
    };
    let cutoff_rho_source = if lone_rho_is_hartree_wfc {
        None
    } else {
        raw_cutoff_rho_source.map(str::to_string)
    };

    PseudopotentialMetadata {
        filename,
        supports_soc: has_so
            || matches!(relativistic.as_deref(), Some(value) if value.eq_ignore_ascii_case("full")),
        pseudo_type,
        relativistic,
        cutoff_wfc,
        cutoff_rho,
        cutoff_wfc_source,
        cutoff_rho_source,
        available_angular_momenta,
        available_angular_momenta_source: angular_momentum_details.source.map(str::to_string),
        max_angular_momentum: angular_momentum_details.max_angular_momentum,
    }
}

fn parse_pseudopotential_metadata_from_sources(
    filename: String,
    content: &str,
    djrepo_text: Option<&str>,
) -> PseudopotentialMetadata {
    let mut metadata = parse_pseudopotential_metadata_from_content(filename, content);
    if let Some(djrepo_text) = djrepo_text {
        if let Some(cutoff_wfc) = parse_djrepo_wavefunction_cutoff_ry(djrepo_text) {
            metadata.cutoff_wfc = Some(cutoff_wfc);
            metadata.cutoff_rho = None;
            metadata.cutoff_wfc_source = Some("djrepo".to_string());
            metadata.cutoff_rho_source = None;
        }
    }
    metadata
}

#[cfg(test)]
mod upf_angular_momentum_tests {
    use super::parse_upf_available_angular_momenta;

    #[test]
    fn parses_old_style_upf_info_tables() {
        let content = r#"
<PP_INFO>
Generated using Vanderbilt code
nl pn  l   occ               Rcut            Rcut US             E pseu
3S  3  0  2.00      0.00000000000      1.44000000000     -9.60941084200
3P  3  1  6.00      0.00000000000      1.55000000000     -6.69550237500
3D  3  2  8.00      0.00000000000      1.50000000000     -2.08482978100
4S  4  0  0.00      0.00000000000      1.44000000000     -1.59107856400
4P  4  1  0.00      0.00000000000      1.55000000000     -1.10274995100
</PP_INFO>
<PP_HEADER>
    2                  Max angular momentum component
 Wavefunctions         nl  l   occ
                       3S  0  2.00
                       3P  1  6.00
                       3D  2  8.00
</PP_HEADER>
"#;

        let parsed = parse_upf_available_angular_momenta(content, content);
        assert_eq!(parsed, vec![0, 1, 2]);
    }

    #[test]
    fn parses_modern_upf_info_generation_tables() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    Valence configuration:
    nl pn  l   occ       Rcut    Rcut US       E pseu
    3S  1  0  2.00      1.600      1.800    -0.794728
    3P  2  1  2.00      1.600      1.800    -0.299965
    Generation configuration:
    3S  1  0  0.00      1.600      1.800     6.000000
    3P  2  1  0.00      1.600      1.800     6.000000
    3D  3  2  0.00      1.600      1.800     0.100000
  </PP_INFO>
  <PP_HEADER l_max="2" />
</UPF>
"#;

        let parsed = parse_upf_available_angular_momenta(content, content);
        assert_eq!(parsed, vec![0, 1, 2]);
    }

    #[test]
    fn falls_back_to_lmax_when_only_header_metadata_exists() {
        let content = r#"<PP_HEADER l_max="1" />"#;
        let parsed = parse_upf_available_angular_momenta(content, "");
        assert_eq!(parsed, vec![0, 1]);
    }

    #[test]
    fn metadata_tracks_confident_vs_fallback_angular_momentum_sources() {
        let si_like = r#"
<UPF version="2.0.1">
  <PP_INFO>
    3S  1  0  2.00      1.600      1.800    -0.794728
    3P  2  1  2.00      1.600      1.800    -0.299965
    3D  3  2  0.00      1.600      1.800     0.100000
  </PP_INFO>
  <PP_HEADER l_max="2" />
  <PP_BETA.1 angular_momentum="0" />
  <PP_BETA.2 angular_momentum="1" />
  <PP_BETA.3 angular_momentum="2" />
</UPF>
"#;
        let si_metadata =
            super::parse_pseudopotential_metadata_from_content("Si.UPF".to_string(), si_like);
        assert_eq!(
            si_metadata.available_angular_momenta_source.as_deref(),
            Some("upf_tags")
        );
        assert_eq!(si_metadata.max_angular_momentum, Some(2));

        let fallback_only = r#"<PP_HEADER l_max="2" />"#;
        let fallback_metadata = super::parse_pseudopotential_metadata_from_content(
            "Fallback.UPF".to_string(),
            fallback_only,
        );
        assert_eq!(
            fallback_metadata
                .available_angular_momenta_source
                .as_deref(),
            Some("upf_l_max_fallback")
        );
        assert_eq!(fallback_metadata.available_angular_momenta, vec![0, 1, 2]);
    }
}

fn decode_remote_metadata_payload(kind: &str, payload: &str) -> Result<String, String> {
    if !matches!(kind, "upf_b64" | "djrepo_b64" | "djrepo_gz_b64") {
        return Ok(payload.to_string());
    }

    let compact_payload = payload.lines().map(str::trim).collect::<String>();
    let decoded = BASE64_STANDARD
        .decode(compact_payload.as_bytes())
        .map_err(|e| format!("Failed to decode remote {} payload: {}", kind, e))?;
    Ok(String::from_utf8_lossy(&decoded).into_owned())
}

fn parse_remote_pseudopotential_metadata_output(
    output: &str,
) -> Result<Vec<PseudopotentialMetadata>, String> {
    let mut upf_contents: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();
    let mut djrepo_contents: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();
    let mut current_kind: Option<String> = None;
    let mut current_filename: Option<String> = None;
    let mut current_content: Vec<String> = Vec::new();

    let flush_current = |kind: &mut Option<String>,
                         filename: &mut Option<String>,
                         content: &mut Vec<String>,
                         upf_contents: &mut std::collections::HashMap<String, String>,
                         djrepo_contents: &mut std::collections::HashMap<String, String>|
     -> Result<(), String> {
        if let (Some(kind), Some(name)) = (kind.take(), filename.take()) {
            let joined = content.join("\n");
            let decoded = decode_remote_metadata_payload(&kind, &joined)?;
            match kind.as_str() {
                "upf" | "upf_b64" => {
                    upf_contents.insert(name, decoded);
                }
                "djrepo" | "djrepo_gz" | "djrepo_b64" | "djrepo_gz_b64" => {
                    djrepo_contents.insert(name, decoded);
                }
                _ => {}
            }
        }
        content.clear();
        Ok(())
    };

    for line in output.lines() {
        let trimmed = line.trim();
        if let Some(path) = trimmed.strip_prefix("__QCORTADO_PSEUDO_DIR_MISSING__:") {
            return Err(format!(
                "Remote pseudopotential directory not found: {}",
                path
            ));
        }
        if let Some(rest) = trimmed.strip_prefix("__QCORTADO_REMOTE_METADATA_FILE__:") {
            flush_current(
                &mut current_kind,
                &mut current_filename,
                &mut current_content,
                &mut upf_contents,
                &mut djrepo_contents,
            )?;
            if let Some((kind, name)) = rest.split_once(':') {
                current_kind = Some(kind.to_string());
                current_filename = Some(name.to_string());
            }
            continue;
        }
        if trimmed == "__QCORTADO_REMOTE_METADATA_FILE_END__" {
            flush_current(
                &mut current_kind,
                &mut current_filename,
                &mut current_content,
                &mut upf_contents,
                &mut djrepo_contents,
            )?;
            continue;
        }
        if current_kind.is_some() {
            current_content.push(line.to_string());
        }
    }

    flush_current(
        &mut current_kind,
        &mut current_filename,
        &mut current_content,
        &mut upf_contents,
        &mut djrepo_contents,
    )?;

    let mut pseudos = Vec::new();
    let mut upf_paths = upf_contents.keys().cloned().collect::<Vec<_>>();
    upf_paths.sort();
    for upf_path in upf_paths {
        let Some(content) = upf_contents.get(&upf_path) else {
            continue;
        };
        let stem = upf_path
            .rsplit_once('.')
            .map(|(prefix, _)| prefix)
            .unwrap_or(upf_path.as_str());
        let djrepo_path = format!("{}.djrepo", stem);
        let djrepo_gz_path = format!("{}.djrepo.gz", stem);
        let djrepo_text = djrepo_contents
            .get(&djrepo_path)
            .or_else(|| djrepo_contents.get(&djrepo_gz_path))
            .map(String::as_str);
        pseudos.push(parse_pseudopotential_metadata_from_sources(
            upf_path,
            content,
            djrepo_text,
        ));
    }
    pseudos.sort_by(|a, b| a.filename.cmp(&b.filename));
    Ok(pseudos)
}

fn read_optional_djrepo_companion(path: &Path) -> Result<Option<String>, String> {
    let Some(parent) = path.parent() else {
        return Ok(None);
    };
    let Some(stem) = path.file_stem().and_then(|value| value.to_str()) else {
        return Ok(None);
    };

    let plain_path = parent.join(format!("{}.djrepo", stem));
    if plain_path.exists() {
        return std::fs::read_to_string(&plain_path).map(Some).map_err(|e| {
            format!(
                "Failed to read companion metadata {}: {}",
                plain_path.display(),
                e
            )
        });
    }

    let gz_path = parent.join(format!("{}.djrepo.gz", stem));
    if gz_path.exists() {
        let file = std::fs::File::open(&gz_path).map_err(|e| {
            format!(
                "Failed to open companion metadata {}: {}",
                gz_path.display(),
                e
            )
        })?;
        let mut decoder = GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content).map_err(|e| {
            format!(
                "Failed to read companion metadata {}: {}",
                gz_path.display(),
                e
            )
        })?;
        return Ok(Some(content));
    }

    Ok(None)
}

fn read_pseudopotential_metadata(path: &Path) -> Result<PseudopotentialMetadata, String> {
    let filename = path
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| format!("Invalid pseudopotential file name: {}", path.display()))?
        .to_string();
    let file = std::fs::File::open(path)
        .map_err(|e| format!("Failed to open pseudopotential {}: {}", filename, e))?;
    let mut buffer = Vec::new();
    file.take(131_072)
        .read_to_end(&mut buffer)
        .map_err(|e| format!("Failed to read pseudopotential {}: {}", filename, e))?;
    let content = String::from_utf8_lossy(&buffer).into_owned();
    let djrepo_content = read_optional_djrepo_companion(path)?;
    Ok(parse_pseudopotential_metadata_from_sources(
        filename,
        &content,
        djrepo_content.as_deref(),
    ))
}

fn is_upf_name(name: &str) -> bool {
    name.ends_with(".UPF") || name.ends_with(".upf")
}

#[cfg(test)]
mod pseudopotential_metadata_tests {
    use base64::Engine as _;

    use super::{
        parse_djrepo_wavefunction_cutoff_ry, parse_pseudopotential_metadata_from_content,
        parse_pseudopotential_metadata_from_sources, parse_remote_pseudopotential_metadata_output,
        BASE64_STANDARD,
    };

    #[test]
    fn parses_v201_header_attributes() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    Suggested minimum cutoff for wavefunctions: 60 Ry
    Suggested minimum cutoff for charge density: 480 Ry
  </PP_INFO>
  <PP_HEADER element="Bi" pseudo_type="US" relativistic="full" has_so="T" wfc_cutoff="60.0" rho_cutoff="480.0"></PP_HEADER>
</UPF>
"#;

        let metadata =
            parse_pseudopotential_metadata_from_content("Bi.rel-pbe.UPF".to_string(), content);
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("US"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        assert_eq!(metadata.cutoff_wfc, Some(60.0));
        assert_eq!(metadata.cutoff_rho, Some(480.0));
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf"));
        assert_eq!(metadata.cutoff_rho_source.as_deref(), Some("upf"));
    }

    #[test]
    fn parses_xml_subtags_for_newer_header_style() {
        let content = r#"
<upf version="3.0.0">
  <pp_header>
    <pseudo_type>PAW</pseudo_type>
    <relativistic>full</relativistic>
    <has_so>true</has_so>
    <wfc_cutoff>45</wfc_cutoff>
    <rho_cutoff>360</rho_cutoff>
  </pp_header>
  <pp_spin_orb />
</upf>
"#;

        let metadata =
            parse_pseudopotential_metadata_from_content("Pt.rel-pbe.UPF".to_string(), content);
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("PAW"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        assert_eq!(metadata.cutoff_wfc, Some(45.0));
        assert_eq!(metadata.cutoff_rho, Some(360.0));
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf"));
        assert_eq!(metadata.cutoff_rho_source.as_deref(), Some("upf"));
    }

    #[test]
    fn interprets_lone_hamann_oncv_rho_cutoff_as_hartree_wavefunction_hint() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>
"#;

        let metadata = parse_pseudopotential_metadata_from_content("Ag.upf".to_string(), content);
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("NC"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        let cutoff_wfc = metadata
            .cutoff_wfc
            .expect("expected converted wavefunction cutoff");
        assert!((cutoff_wfc - 30.74).abs() < 1e-9);
        assert_eq!(metadata.cutoff_rho, None);
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf_fallback"));
        assert_eq!(metadata.cutoff_rho_source, None);
    }

    #[test]
    fn parses_pseudodojo_djrepo_normal_hint_as_ry_wavefunction_cutoff() {
        let djrepo = r#"
{
  "hints": {
    "high": { "ecut": 47.0 },
    "low": { "ecut": 37.0 },
    "normal": { "ecut": 41.0 }
  }
}
"#;

        let cutoff_wfc = parse_djrepo_wavefunction_cutoff_ry(djrepo)
            .expect("expected djrepo wavefunction cutoff");
        assert!((cutoff_wfc - 82.0).abs() < 1e-9);
    }

    #[test]
    fn prefers_djrepo_cutoffs_over_embedded_upf_cutoff_hints() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>
"#;
        let djrepo = r#"
{
  "hints": {
    "normal": { "ecut": 41.0 }
  }
}
"#;

        let metadata = parse_pseudopotential_metadata_from_sources(
            "Ag.upf".to_string(),
            content,
            Some(djrepo),
        );
        let cutoff_wfc = metadata
            .cutoff_wfc
            .expect("expected djrepo wavefunction cutoff");
        assert!((cutoff_wfc - 82.0).abs() < 1e-9);
        assert_eq!(metadata.cutoff_rho, None);
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("djrepo"));
        assert_eq!(metadata.cutoff_rho_source, None);
    }

    #[test]
    fn pairs_remote_djrepo_companions_even_with_banner_noise() {
        let output = r#"
Welcome to the cluster.
Authorized use only.
__QCORTADO_REMOTE_METADATA_FILE__:djrepo:Ag.djrepo
{
  "hints": {
    "normal": { "ecut": 41.0 }
  }
}
__QCORTADO_REMOTE_METADATA_FILE_END__
__QCORTADO_REMOTE_METADATA_FILE__:upf:Ag.upf
<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>
__QCORTADO_REMOTE_METADATA_FILE_END__
"#;

        let pseudos =
            parse_remote_pseudopotential_metadata_output(output).expect("expected parsed metadata");
        assert_eq!(pseudos.len(), 1);
        assert_eq!(pseudos[0].filename, "Ag.upf");
        assert!((pseudos[0].cutoff_wfc.expect("expected djrepo cutoff") - 82.0).abs() < 1e-9);
        assert_eq!(pseudos[0].cutoff_rho, None);
        assert_eq!(pseudos[0].cutoff_wfc_source.as_deref(), Some("djrepo"));
        assert_eq!(pseudos[0].cutoff_rho_source, None);
    }

    #[test]
    fn parses_base64_remote_metadata_payloads_for_djrepo_companions() {
        let djrepo = r#"{
  "hints": {
    "normal": { "ecut": 41.0 }
  }
}"#;
        let upf = r#"<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>"#;
        let output = format!(
            "Welcome to the cluster.\n__QCORTADO_REMOTE_METADATA_FILE__:djrepo_b64:Ag.djrepo\n{}\n__QCORTADO_REMOTE_METADATA_FILE_END__\n__QCORTADO_REMOTE_METADATA_FILE__:upf_b64:Ag.upf\n{}\n__QCORTADO_REMOTE_METADATA_FILE_END__\n",
            BASE64_STANDARD.encode(djrepo.as_bytes()),
            BASE64_STANDARD.encode(upf.as_bytes()),
        );

        let pseudos = parse_remote_pseudopotential_metadata_output(&output)
            .expect("expected parsed metadata");
        assert_eq!(pseudos.len(), 1);
        assert_eq!(pseudos[0].filename, "Ag.upf");
        assert!((pseudos[0].cutoff_wfc.expect("expected djrepo cutoff") - 82.0).abs() < 1e-9);
        assert_eq!(pseudos[0].cutoff_rho, None);
        assert_eq!(pseudos[0].cutoff_wfc_source.as_deref(), Some("djrepo"));
        assert_eq!(pseudos[0].cutoff_rho_source, None);
    }
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
        if !entry.path().is_file() {
            continue;
        }

        let name = entry.file_name().to_string_lossy().to_string();
        if is_upf_name(&name) {
            pseudos.push(name);
        }
    }
    pseudos.sort();
    Ok(pseudos)
}

/// Lists pseudopotentials and parses SOC/cutoff metadata from their headers.
fn list_pseudopotential_metadata_sync(
    pseudo_dir: String,
) -> Result<Vec<PseudopotentialMetadata>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    let mut pseudos = Vec::new();
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let file_path = entry.path();
        if !file_path.is_file() {
            continue;
        }

        let name = entry.file_name().to_string_lossy().to_string();
        if is_upf_name(&name) {
            pseudos.push(read_pseudopotential_metadata(&file_path)?);
        }
    }
    pseudos.sort_by(|a, b| a.filename.cmp(&b.filename));
    Ok(pseudos)
}

/// Lists pseudopotentials and parses SOC/cutoff metadata from their headers.
/// Runs in a blocking worker so the UI thread remains responsive while metadata is indexed.
#[tauri::command]
async fn list_pseudopotential_metadata(
    pseudo_dir: String,
) -> Result<Vec<PseudopotentialMetadata>, String> {
    tokio::task::spawn_blocking(move || list_pseudopotential_metadata_sync(pseudo_dir))
        .await
        .map_err(|e| format!("Failed to join pseudopotential metadata task: {}", e))?
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
fn load_sssp_data_sync(
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

/// Loads SSSP JSON data from the pseudo directory.
/// Runs in a blocking worker so the UI can show loading state immediately.
#[tauri::command]
async fn load_sssp_data(
    pseudo_dir: String,
) -> Result<std::collections::HashMap<String, SSSPElementData>, String> {
    tokio::task::spawn_blocking(move || load_sssp_data_sync(pseudo_dir))
        .await
        .map_err(|e| format!("Failed to join SSSP loading task: {}", e))?
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
    /// Optional bands.x post-processing settings
    pub bands_x: Option<BandsPostProcessingOptions>,
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
    /// projwfc.x output file name (filproj)
    pub filproj: Option<String>,
}

#[derive(Debug, Clone, serde::Deserialize)]
pub struct BandsPostProcessingOptions {
    /// bands.x output file (filband)
    pub filband: Option<String>,
    /// bands.x lsym option (print symmetry labels/order)
    pub lsym: Option<bool>,
    /// bands.x no_overlap option (overlap-based reordering toggle)
    pub no_overlap: Option<bool>,
}

/// Configuration for electronic DOS calculation (NSCF + dos.x).
#[derive(Debug, Clone, serde::Deserialize)]
pub struct ElectronicDosCalculationConfig {
    /// Base SCF calculation to derive system settings from
    pub base_calculation: QECalculation,
    /// Dense automatic k-grid for NSCF DOS run
    pub k_grid: [u32; 3],
    /// Optional DOS broadening (Ry)
    pub degauss: Option<f64>,
    /// Optional minimum energy window (eV)
    pub emin: Option<f64>,
    /// Optional maximum energy window (eV)
    pub emax: Option<f64>,
    /// Optional DOS energy step (eV)
    pub delta_e: Option<f64>,
    /// Project ID containing the source SCF calculation
    pub project_id: Option<String>,
    /// SCF calculation ID to get the .save directory from
    pub scf_calc_id: Option<String>,
}

/// Parsed electronic DOS data returned to the frontend.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ElectronicDosData {
    /// Energy samples in eV
    pub energies: Vec<f64>,
    /// DOS values for each energy sample
    pub dos: Vec<f64>,
    /// Fermi energy in eV (if detected from NSCF output)
    pub fermi_energy: Option<f64>,
    /// [min, max] sampled energy range
    pub energy_range: [f64; 2],
    /// Maximum DOS value in sampled data
    pub max_dos: f64,
    /// Number of sampled DOS points
    pub points: usize,
}

/// Configuration for Fermi-surface generation via QE's FermiSurfer pipeline.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct FermiSurfaceCalculationConfig {
    /// Base SCF calculation to derive system settings from
    pub base_calculation: QECalculation,
    /// Automatic k-grid used for the fermi_velocity.x input
    pub k_grid: [u32; 3],
    /// Optional automatic k-grid offset (0/1 for each reciprocal direction)
    pub k_offset: Option<[u32; 3]>,
    /// Optional number of bands for the generated pw-style input (auto if None)
    pub nbnd: Option<u32>,
    /// Project ID containing the source SCF calculation
    pub project_id: Option<String>,
    /// SCF calculation ID to get the .save directory from
    pub scf_calc_id: Option<String>,
}

/// Generated FermiSurfer FRMSF file metadata.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FrmsfFileData {
    /// File name in the working directory
    pub file_name: String,
    /// File size in bytes
    pub size_bytes: u64,
}

/// Parsed Fermi surface generation output returned to the frontend.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FermiSurfaceData {
    /// K-grid used for generation
    pub k_grid: [u32; 3],
    /// Fermi energy in eV (if detected from fermi_velocity.x output)
    pub fermi_energy: Option<f64>,
    /// Preferred FRMSF file to open by default
    pub primary_file: String,
    /// All generated FRMSF files
    pub frmsf_files: Vec<FrmsfFileData>,
}

fn sanitize_output_filename(raw: &str, fallback: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return fallback.to_string();
    }

    let mut sanitized: String = trimmed
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '_' || ch == '-' || ch == '.' {
                ch
            } else {
                '_'
            }
        })
        .collect();
    sanitized = sanitized.trim_matches('_').to_string();

    if sanitized.is_empty() {
        fallback.to_string()
    } else {
        sanitized
    }
}

fn missing_scf_tmp_error(path: &Path) -> String {
    format!(
        "SCF calculation tmp directory not found: {}. If this SCF ran on HPC with minimal sync, open Cluster Activity and use 'Download Full'.",
        path.display()
    )
}

#[derive(Debug, Default)]
struct HpcScfDependencyStage {
    local_bundle_copy: Option<PathBuf>,
    remote_hydration_commands: Vec<String>,
}

fn contains_qe_save_directory(root: &Path, max_depth: usize) -> bool {
    if !root.exists() || !root.is_dir() {
        return false;
    }

    let mut stack: Vec<(PathBuf, usize)> = vec![(root.to_path_buf(), 0)];
    while let Some((current, depth)) = stack.pop() {
        let entries = match std::fs::read_dir(&current) {
            Ok(entries) => entries,
            Err(_) => continue,
        };

        for entry in entries.flatten() {
            let file_type = match entry.file_type() {
                Ok(file_type) => file_type,
                Err(_) => continue,
            };
            if !file_type.is_dir() {
                continue;
            }

            let name = entry.file_name().to_string_lossy().to_string();
            if name.ends_with(".save") {
                return true;
            }

            if depth < max_depth {
                stack.push((entry.path(), depth + 1));
            }
        }
    }

    false
}

fn find_qe_save_directory(
    root: &Path,
    max_depth: usize,
    preferred_dir_name: Option<&str>,
) -> Option<PathBuf> {
    if !root.exists() || !root.is_dir() {
        return None;
    }

    let mut stack: Vec<(PathBuf, usize)> = vec![(root.to_path_buf(), 0)];
    let mut first_any: Option<PathBuf> = None;
    while let Some((current, depth)) = stack.pop() {
        let entries = match std::fs::read_dir(&current) {
            Ok(entries) => entries,
            Err(_) => continue,
        };

        for entry in entries.flatten() {
            let file_type = match entry.file_type() {
                Ok(file_type) => file_type,
                Err(_) => continue,
            };
            if !file_type.is_dir() {
                continue;
            }

            let name = entry.file_name().to_string_lossy().to_string();
            if name.ends_with(".save") {
                if preferred_dir_name
                    .map(|preferred| name == preferred)
                    .unwrap_or(false)
                {
                    return Some(entry.path());
                }
                if first_any.is_none() {
                    first_any = Some(entry.path());
                }
            }

            if depth < max_depth {
                stack.push((entry.path(), depth + 1));
            }
        }
    }

    first_any
}

fn find_preferred_wannier_save_directory(root: &Path, preferred_dir_name: &str) -> Option<PathBuf> {
    let tmp_root = root.join("tmp");
    find_qe_save_directory(&tmp_root, 2, Some(preferred_dir_name))
        .or_else(|| find_qe_save_directory(&tmp_root, 2, None))
        .or_else(|| find_qe_save_directory(root, 4, Some(preferred_dir_name)))
        .or_else(|| find_qe_save_directory(root, 4, None))
}

fn source_calc_remote_dependency_path(calc: &projects::CalculationRun) -> Option<String> {
    let params = calc.parameters.as_object()?;
    params
        .get("remote_project_path")
        .or_else(|| params.get("remote_workdir"))
        .and_then(|value| value.as_str())
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(ToString::to_string)
}

fn build_remote_scf_hydration_commands(remote_scf_dir: &str) -> Vec<String> {
    vec![
        format!("SOURCE_SCF_DIR={}", shell_single_quote_local(remote_scf_dir)),
        "echo \"HPC_STAGE|HydratingDependency|$SOURCE_SCF_DIR\"".to_string(),
        "if [ ! -d \"$SOURCE_SCF_DIR\" ]; then echo \"[QCortado] ERROR: Source SCF dependency directory not found on remote: $SOURCE_SCF_DIR\"; exit 30; fi".to_string(),
        "if [ ! -d \"$SOURCE_SCF_DIR/tmp\" ]; then echo \"[QCortado] ERROR: Source SCF tmp directory not found on remote: $SOURCE_SCF_DIR/tmp\"; exit 31; fi".to_string(),
        "mkdir -p ./tmp".to_string(),
        "cp -a \"$SOURCE_SCF_DIR/tmp/.\" ./tmp/".to_string(),
        "if [ -z \"$(find ./tmp -maxdepth 2 -type d -name '*.save' -print -quit)\" ]; then echo \"[QCortado] ERROR: No QE .save directory found after SCF dependency hydrate\"; exit 32; fi".to_string(),
    ]
}

fn resolve_hpc_scf_dependency_stage(
    app: &AppHandle,
    project_id: Option<&str>,
    scf_calc_id: Option<&str>,
) -> Result<HpcScfDependencyStage, String> {
    let project_id = project_id.map(str::trim).filter(|value| !value.is_empty());
    let scf_calc_id = scf_calc_id.map(str::trim).filter(|value| !value.is_empty());

    let (project_id, scf_calc_id) = match (project_id, scf_calc_id) {
        (Some(project_id), Some(scf_calc_id)) => (project_id, scf_calc_id),
        _ => return Ok(HpcScfDependencyStage::default()),
    };

    let projects_dir = projects::get_projects_dir(app)?;
    let scf_tmp_dir = projects_dir
        .join(project_id)
        .join("calculations")
        .join(scf_calc_id)
        .join("tmp");
    let local_has_save_payload =
        scf_tmp_dir.exists() && contains_qe_save_directory(&scf_tmp_dir, 2);

    let project = projects::get_project(app.clone(), project_id.to_string())?;
    let source_calc = project
        .cif_variants
        .iter()
        .flat_map(|variant| variant.calculations.iter())
        .find(|calc| calc.id == scf_calc_id)
        .ok_or_else(|| {
            format!(
                "Source SCF calculation {} not found in project {}.",
                scf_calc_id, project_id
            )
        })?;

    if let Some(remote_source) = source_calc_remote_dependency_path(source_calc) {
        return Ok(HpcScfDependencyStage {
            local_bundle_copy: None,
            remote_hydration_commands: build_remote_scf_hydration_commands(&remote_source),
        });
    }

    if local_has_save_payload {
        return Ok(HpcScfDependencyStage {
            local_bundle_copy: Some(scf_tmp_dir),
            remote_hydration_commands: Vec::new(),
        });
    }

    Err(missing_scf_tmp_error(&scf_tmp_dir))
}

fn insert_system_namelist_line(input: &str, line_to_insert: &str) -> Result<String, String> {
    let mut output = String::with_capacity(input.len() + line_to_insert.len() + 16);
    let mut in_system = false;
    let mut inserted = false;

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.eq_ignore_ascii_case("&SYSTEM") {
            in_system = true;
            output.push_str(line);
            output.push('\n');
            continue;
        }

        if in_system && trimmed == "/" {
            output.push_str("  ");
            output.push_str(line_to_insert);
            output.push('\n');
            inserted = true;
            in_system = false;
        }

        output.push_str(line);
        output.push('\n');
    }

    if inserted {
        Ok(output)
    } else {
        Err("Could not locate &SYSTEM namelist in pw.x input".to_string())
    }
}

fn executable_uses_pencil_decomposition(executable: &Path) -> bool {
    executable
        .file_name()
        .and_then(|value| value.to_str())
        .map(qe_executable_uses_pencil_decomposition)
        .unwrap_or(false)
}

fn stage_args_include_pencil_decomposition(args: &[&str]) -> bool {
    args.iter()
        .any(|value| command_args_include_pencil_decomposition(value))
}

#[allow(clippy::too_many_arguments)]
async fn run_local_stage_capture_stdout(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    work_path: &Path,
    executable: &Path,
    args: &[&str],
    stdin_content: Option<&str>,
    execution_prefix: Option<&str>,
    mpi_config: Option<&MpiConfig>,
    allow_mpi: bool,
) -> Result<String, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let force_pencil_decomposition = executable_uses_pencil_decomposition(executable)
        && !stage_args_include_pencil_decomposition(args);
    let mut child = if allow_mpi
        && mpi_config
            .map(|mpi| mpi.enabled && mpi.nprocs > 1)
            .unwrap_or(false)
    {
        let mpi = mpi_config.unwrap();
        let mut command = tokio_command_with_prefix("mpirun", execution_prefix);
        command.args(["-np", &mpi.nprocs.to_string()]);
        command.arg(executable);
        command.args(args);
        if force_pencil_decomposition {
            command.args(["-pd", ".true."]);
        }
        command
    } else {
        let mut command = tokio_command_with_prefix(executable, execution_prefix);
        command.args(args);
        command
    };

    child.current_dir(work_path);
    if stdin_content.is_some() {
        child.stdin(Stdio::piped());
    } else {
        child.stdin(Stdio::null());
    }
    child.stdout(Stdio::piped()).stderr(Stdio::piped());

    let mut child = child
        .spawn()
        .map_err(|e| format!("Failed to start {}: {}", executable.display(), e))?;

    if let Some(pid) = child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(content) = stdin_content {
        if let Some(mut stdin) = child.stdin.take() {
            stdin
                .write_all(content.as_bytes())
                .await
                .map_err(|e| format!("Failed to write stage input: {}", e))?;
        }
    }

    let stdout = child.stdout.take().ok_or("Failed to capture stdout")?;
    let mut reader = BufReader::new(stdout).lines();
    let mut output = String::new();
    while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
        output.push_str(&line);
        output.push('\n');
        let _ = app.emit(&format!("task-output:{}", task_id), &line);
        pm.append_output(task_id, line).await;
    }

    let status = child.wait().await.map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!(
            "{} failed with exit code {:?}",
            executable
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("process"),
            status.code()
        ));
    }

    Ok(output)
}

struct LocalStageStreams {
    stdout: String,
    stderr: String,
}

async fn run_local_stage_capture_streams(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    work_path: &Path,
    executable: &Path,
    args: &[&str],
    stdin_content: Option<&str>,
    execution_prefix: Option<&str>,
    mpi_config: Option<&MpiConfig>,
    allow_mpi: bool,
) -> Result<LocalStageStreams, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let force_pencil_decomposition = executable_uses_pencil_decomposition(executable)
        && !stage_args_include_pencil_decomposition(args);
    let mut child = if allow_mpi
        && mpi_config
            .map(|mpi| mpi.enabled && mpi.nprocs > 1)
            .unwrap_or(false)
    {
        let mpi = mpi_config.unwrap();
        let mut command = tokio_command_with_prefix("mpirun", execution_prefix);
        command.args(["-np", &mpi.nprocs.to_string()]);
        command.arg(executable);
        command.args(args);
        if force_pencil_decomposition {
            command.args(["-pd", ".true."]);
        }
        command
    } else {
        let mut command = tokio_command_with_prefix(executable, execution_prefix);
        command.args(args);
        command
    };

    child.current_dir(work_path);
    if stdin_content.is_some() {
        child.stdin(Stdio::piped());
    } else {
        child.stdin(Stdio::null());
    }
    child.stdout(Stdio::piped()).stderr(Stdio::piped());

    let mut child = child
        .spawn()
        .map_err(|e| format!("Failed to start {}: {}", executable.display(), e))?;

    if let Some(pid) = child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(content) = stdin_content {
        if let Some(mut stdin) = child.stdin.take() {
            stdin
                .write_all(content.as_bytes())
                .await
                .map_err(|e| format!("Failed to write stage input: {}", e))?;
        }
    }

    let stdout = child.stdout.take().ok_or("Failed to capture stdout")?;
    let stderr = child.stderr.take().ok_or("Failed to capture stderr")?;

    let stdout_task_id = task_id.to_string();
    let stdout_app = app.clone();
    let stdout_pm = pm.clone();
    let stdout_task = tokio::spawn(async move {
        let mut output = String::new();
        let mut reader = BufReader::new(stdout).lines();
        while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
            output.push_str(&line);
            output.push('\n');
            let _ = stdout_app.emit(&format!("task-output:{}", stdout_task_id), &line);
            stdout_pm.append_output(&stdout_task_id, line).await;
        }
        Ok::<String, String>(output)
    });

    let stderr_task_id = task_id.to_string();
    let stderr_app = app.clone();
    let stderr_pm = pm.clone();
    let stderr_task = tokio::spawn(async move {
        let mut output = String::new();
        let mut reader = BufReader::new(stderr).lines();
        while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
            output.push_str(&line);
            output.push('\n');
            let _ = stderr_app.emit(&format!("task-output:{}", stderr_task_id), &line);
            stderr_pm.append_output(&stderr_task_id, line).await;
        }
        Ok::<String, String>(output)
    });

    let status = child.wait().await.map_err(|e| e.to_string())?;
    let stdout_output = stdout_task
        .await
        .map_err(|e| format!("Failed to join stdout task: {}", e))??;
    let stderr_output = stderr_task
        .await
        .map_err(|e| format!("Failed to join stderr task: {}", e))??;

    if !status.success() {
        return Err(format!(
            "{} failed with exit code {:?}",
            executable
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("process"),
            status.code()
        ));
    }

    Ok(LocalStageStreams {
        stdout: stdout_output,
        stderr: stderr_output,
    })
}

fn collect_surface_files(
    work_path: &Path,
    extensions: &[&str],
) -> Result<Vec<FrmsfFileData>, String> {
    let mut files: Vec<FrmsfFileData> = Vec::new();
    collect_surface_files_recursive(work_path, work_path, extensions, 0, 8, &mut files)?;

    files.sort_by(|a, b| a.file_name.cmp(&b.file_name));
    Ok(files)
}

fn collect_surface_files_recursive(
    root: &Path,
    current: &Path,
    extensions: &[&str],
    depth: usize,
    max_depth: usize,
    files: &mut Vec<FrmsfFileData>,
) -> Result<(), String> {
    if depth > max_depth {
        return Ok(());
    }

    let entries = std::fs::read_dir(current).map_err(|e| {
        format!(
            "Failed to inspect working directory {}: {}",
            current.display(),
            e
        )
    })?;

    for entry in entries {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        let file_type = entry.file_type().map_err(|e| e.to_string())?;

        if file_type.is_dir() {
            let dir_name = entry.file_name().to_string_lossy().to_string();
            if dir_name.ends_with(".save") {
                continue;
            }
            collect_surface_files_recursive(root, &path, extensions, depth + 1, max_depth, files)?;
            continue;
        }

        if !file_type.is_file() {
            continue;
        }

        let matches_extension = path
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| {
                extensions
                    .iter()
                    .any(|candidate| ext.eq_ignore_ascii_case(candidate))
            })
            .unwrap_or(false);
        if !matches_extension {
            continue;
        }

        let relative_path = path.strip_prefix(root).unwrap_or(&path);
        let file_name = relative_path.to_string_lossy().to_string();
        let size_bytes = std::fs::metadata(&path).map(|meta| meta.len()).unwrap_or(0);
        files.push(FrmsfFileData {
            file_name,
            size_bytes,
        });
    }

    Ok(())
}

fn collect_frmsf_files(work_path: &Path) -> Result<Vec<FrmsfFileData>, String> {
    collect_surface_files(work_path, &["frmsf"])
}

fn collect_viewer_input_files(work_path: &Path) -> Result<Vec<FrmsfFileData>, String> {
    let mut files = collect_surface_files(work_path, &["frmsf", "bxsf"])?;
    files.sort_by(|a, b| {
        let a_ext = Path::new(&a.file_name)
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.to_ascii_lowercase())
            .unwrap_or_default();
        let b_ext = Path::new(&b.file_name)
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.to_ascii_lowercase())
            .unwrap_or_default();

        let a_rank = if a_ext == "frmsf" { 0 } else { 1 };
        let b_rank = if b_ext == "frmsf" { 0 } else { 1 };
        a_rank
            .cmp(&b_rank)
            .then_with(|| a.file_name.cmp(&b.file_name))
    });
    Ok(files)
}

fn is_safe_surface_file_name(file_name: &str) -> bool {
    let trimmed = file_name.trim();
    if trimmed.is_empty() {
        return false;
    }

    let path = Path::new(trimmed);
    if path.is_absolute() {
        return false;
    }
    if !path
        .components()
        .all(|component| matches!(component, std::path::Component::Normal(_)))
    {
        return false;
    }

    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("frmsf") || ext.eq_ignore_ascii_case("bxsf"))
        .unwrap_or(false)
}

/// Launches FermiSurfer for a saved Fermi-surface file.
#[tauri::command]
fn launch_fermi_surface_viewer(
    app: AppHandle,
    project_id: String,
    calculation_id: String,
    surface_file: Option<String>,
    state: State<AppState>,
) -> Result<(), String> {
    let fermi_surfer_path = {
        let guard = state.fermi_surfer_path.lock().unwrap();
        guard
            .as_ref()
            .ok_or("FermiSurfer path not configured. Configure it on the home page.")?
            .clone()
    };

    if !fermi_surfer_path.exists() || !fermi_surfer_path.is_file() {
        return Err(format!(
            "Configured FermiSurfer executable is invalid: {}",
            fermi_surfer_path.display()
        ));
    }

    let projects_dir = projects::get_projects_dir(&app)?;
    let calc_tmp_dir = projects_dir
        .join(&project_id)
        .join("calculations")
        .join(&calculation_id)
        .join("tmp");
    if !calc_tmp_dir.exists() {
        return Err(format!(
            "Calculation working directory not found: {}",
            calc_tmp_dir.display()
        ));
    }

    let available_files = collect_viewer_input_files(&calc_tmp_dir)?;
    if available_files.is_empty() {
        return Err(
            "No .frmsf or .bxsf files found for this Fermi-surface calculation.".to_string(),
        );
    }

    let requested_file = surface_file.and_then(|raw| {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    });

    let selected_file_name = if let Some(requested) = requested_file {
        if !is_safe_surface_file_name(&requested) {
            return Err("Invalid Fermi-surface file name.".to_string());
        }
        if !available_files
            .iter()
            .any(|file| file.file_name == requested)
        {
            return Err(format!(
                "Requested Fermi-surface file not found: {}",
                requested
            ));
        }
        requested
    } else {
        available_files[0].file_name.clone()
    };

    let surface_path = calc_tmp_dir.join(&selected_file_name);
    if !surface_path.exists() || !surface_path.is_file() {
        return Err(format!(
            "Fermi-surface file does not exist: {}",
            surface_path.display()
        ));
    }

    let launch_log_path = calc_tmp_dir.join("fermisurfer.launch.log");
    let mut launch_log = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&launch_log_path)
        .map_err(|e| {
            format!(
                "Failed to create launch log {}: {}",
                launch_log_path.display(),
                e
            )
        })?;
    use std::io::Write;
    let _ = writeln!(
        launch_log,
        "Launching FermiSurfer:\n  executable={}\n  input={}",
        fermi_surfer_path.display(),
        surface_path.display()
    );
    let stdout_log = launch_log
        .try_clone()
        .map_err(|e| format!("Failed to prepare launch log: {}", e))?;
    let stderr_log = launch_log
        .try_clone()
        .map_err(|e| format!("Failed to prepare launch log: {}", e))?;

    let mut child = std::process::Command::new(&fermi_surfer_path)
        .arg(&surface_path)
        .current_dir(&calc_tmp_dir)
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::from(stdout_log))
        .stderr(std::process::Stdio::from(stderr_log))
        .spawn()
        .map_err(|e| format!("Failed to launch FermiSurfer: {}", e))?;

    // Give GUI startup a brief window. If it exits immediately, surface logs to the user.
    std::thread::sleep(std::time::Duration::from_millis(800));
    if let Some(status) = child
        .try_wait()
        .map_err(|e| format!("Failed while checking FermiSurfer startup: {}", e))?
    {
        let launch_log_content = std::fs::read_to_string(&launch_log_path).unwrap_or_default();
        let mut excerpt = launch_log_content.trim().to_string();
        if excerpt.len() > 3000 {
            excerpt = format!("{}...", &excerpt[..3000]);
        }
        let detail = if excerpt.is_empty() {
            format!("See {}", launch_log_path.display())
        } else {
            format!("Log ({}):\n{}", launch_log_path.display(), excerpt)
        };
        return Err(format!(
            "FermiSurfer exited immediately with status {:?}. {}",
            status.code(),
            detail
        ));
    }

    Ok(())
}

#[tauri::command]
fn export_wannier_for_ludwig(
    app: AppHandle,
    mut config: LudwigExportConfig,
) -> Result<LudwigExportResult, String> {
    config.project_id = config.project_id.trim().to_string();
    config.calculation_id = config.calculation_id.trim().to_string();
    config.destination_root = config.destination_root.trim().to_string();
    if config.project_id.is_empty() {
        return Err("Project ID is required for Ludwig export.".to_string());
    }
    if config.calculation_id.is_empty() {
        return Err("Calculation ID is required for Ludwig export.".to_string());
    }
    if config.destination_root.is_empty() {
        return Err("Destination directory is required for Ludwig export.".to_string());
    }

    let calculation = projects::get_project_calculation(
        app.clone(),
        config.project_id.clone(),
        config.calculation_id.clone(),
    )?;

    if calculation.calc_type != "wannier" {
        return Err(format!(
            "Calculation {} is not a Wannier calculation.",
            config.calculation_id
        ));
    }

    let projects_dir = projects::get_projects_dir(&app)?;
    let calc_tmp_dir = projects_dir
        .join(&config.project_id)
        .join("calculations")
        .join(&config.calculation_id)
        .join("tmp");
    if !calc_tmp_dir.exists() {
        return Err(format!(
            "Saved calculation working directory not found: {}",
            calc_tmp_dir.display()
        ));
    }

    let seedname = calculation
        .result
        .as_ref()
        .and_then(|result| result.wannier_data.as_ref())
        .and_then(|payload| payload.get("seedname"))
        .and_then(|value| value.as_str())
        .or_else(|| {
            calculation
                .parameters
                .get("seedname")
                .and_then(|value| value.as_str())
        })
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or("qcortado_wannier")
        .to_string();

    let win_path = calc_tmp_dir.join(format!("{}.win", seedname));
    let hr_path = calc_tmp_dir.join(format!("{}_hr.dat", seedname));
    let wsvec_path = calc_tmp_dir.join(format!("{}_wsvec.dat", seedname));
    let wout_path = calc_tmp_dir.join(format!("{}.wout", seedname));

    let win_content = std::fs::read_to_string(&win_path)
        .map_err(|e| format!("Failed to read {}: {}", win_path.display(), e))?;
    let hr_content = std::fs::read_to_string(&hr_path)
        .map_err(|e| format!("Failed to read {}: {}", hr_path.display(), e))?;
    let wsvec_content = if wsvec_path.exists() {
        Some(
            std::fs::read_to_string(&wsvec_path)
                .map_err(|e| format!("Failed to read {}: {}", wsvec_path.display(), e))?,
        )
    } else {
        None
    };

    let hamiltonian = parse_wannier_hamiltonian(
        &seedname,
        &win_content,
        &hr_content,
        wsvec_content.as_deref(),
    )?;
    let chemical_potential_ev = config
        .chemical_potential_ev
        .or_else(|| calculation.result.as_ref().and_then(|result| result.fermi_energy))
        .ok_or_else(|| {
            "Ludwig export requires a chemical potential. No override was provided and the saved Wannier calculation has no Fermi energy.".to_string()
        })?;

    let mut provenance_files = vec![win_path, hr_path, wout_path];
    if wsvec_path.exists() {
        provenance_files.push(wsvec_path);
    }

    export_ludwig_bundle(
        &hamiltonian,
        &config,
        chemical_potential_ev,
        &provenance_files,
    )
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

    prepare_working_directory(&work_path, false)?;

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
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    // Step 1: Create bands calculation from base SCF
    let mut bands_calc = config.base_calculation.clone();
    bands_calc.calculation = qe::CalculationType::Bands;
    if bands_calc.verbosity.is_none() {
        bands_calc.verbosity = Some("high".to_string());
    }

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
    let mut input = generate_pw_input(&bands_calc);
    if let Some(nbnd) = config.nbnd {
        input = insert_system_namelist_line(&input, &format!("nbnd = {},", nbnd))?;
        let _ = app.emit("qe-output-line", format!("Requested nbnd = {}", nbnd));
    }

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
                .args(["-pd", ".true."])
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

    let bands_x_options = config.bands_x.as_ref();
    let bands_filband = bands_x_options
        .and_then(|opts| opts.filband.as_deref())
        .map(|raw| sanitize_output_filename(raw, "bands.dat"))
        .unwrap_or_else(|| "bands.dat".to_string());
    let bands_lsym = bands_x_options.and_then(|opts| opts.lsym).unwrap_or(true);
    let bands_no_overlap = bands_x_options
        .and_then(|opts| opts.no_overlap)
        .unwrap_or(true);

    let bands_x_config = BandsXConfig {
        prefix: bands_calc.prefix.clone(),
        outdir: bands_calc.outdir.clone(),
        filband: bands_filband.clone(),
        lsym: bands_lsym,
        no_overlap: bands_no_overlap,
    };
    let bands_x_input = generate_bands_x_input(&bands_x_config);

    // Save bands.x input for reference
    std::fs::write(work_path.join("bands_pp.in"), &bands_x_input)
        .map_err(|e| format!("Failed to write bands.x input: {}", e))?;

    let mut bands_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            let _ = app.emit(
                "qe-output-line",
                format!("Running bands.x with MPI ({} processes)", mpi.nprocs),
            );
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&bands_x_path)
                .args(["-pd", ".true."])
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun for bands.x: {}", e))?
        } else {
            tokio_command_with_prefix(&bands_x_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start bands.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&bands_x_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start bands.x: {}", e))?
    };

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
    let gnu_file_name = format!("{}.gnu", bands_x_config.filband);
    let gnu_file = work_path.join(&gnu_file_name);
    if !gnu_file.exists() {
        return Err(format!(
            "{} not found. bands.x may have failed.",
            gnu_file_name
        ));
    }

    // Log file size for debugging
    if let Ok(metadata) = std::fs::metadata(&gnu_file) {
        let _ = app.emit(
            "qe-output-line",
            format!("{} size: {} bytes", gnu_file_name, metadata.len()),
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
            let projection_file = projection_options
                .and_then(|opts| opts.filproj.as_deref())
                .map(|raw| sanitize_output_filename(raw, "bands.projwfc.dat"))
                .unwrap_or_else(|| "bands.projwfc.dat".to_string());
            let projwfc_config = ProjwfcConfig {
                prefix: bands_calc.prefix.clone(),
                outdir: bands_calc.outdir.clone(),
                filproj: projection_file,
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
                    match parse_projwfc_projection_groups_aligned(
                        &projection_text,
                        &band_data.energies,
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

    // Keep existing scratch only when explicit recover mode is enabled.
    prepare_working_directory(&work_path, config.phonon.recover)?;

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
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    ensure_phonon_restart_inputs(&work_path)?;

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

    let mut q2r_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            let _ = app.emit(
                "qe-output-line",
                format!("Running q2r.x with MPI ({} processes)", mpi.nprocs),
            );
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&q2r_exe)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun for q2r.x: {}", e))?
        } else {
            tokio_command_with_prefix(&q2r_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start q2r.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&q2r_exe, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start q2r.x: {}", e))?
    };

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

        let mut matdyn_dos_child = if let Some(ref mpi) = mpi_config {
            if mpi.enabled && mpi.nprocs > 1 {
                let _ = app.emit(
                    "qe-output-line",
                    format!("Running matdyn.x (DOS) with MPI ({} processes)", mpi.nprocs),
                );
                tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                    .args(["-np", &mpi.nprocs.to_string()])
                    .arg(&matdyn_exe)
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start mpirun for matdyn.x DOS: {}", e))?
            } else {
                tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start matdyn.x for DOS: {}", e))?
            }
        } else {
            tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start matdyn.x for DOS: {}", e))?
        };

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

            // Preserve per-segment interpolation from the frontend path payload.
            let q_path_with_points: Vec<QPathPoint> = q_path
                .iter()
                .enumerate()
                .map(|(i, p)| QPathPoint {
                    label: p.label.clone(),
                    coords: p.coords,
                    npoints: if i < q_path.len() - 1 { p.npoints } else { 0 },
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

            let mut matdyn_bands_child = if let Some(ref mpi) = mpi_config {
                if mpi.enabled && mpi.nprocs > 1 {
                    let _ = app.emit(
                        "qe-output-line",
                        format!(
                            "Running matdyn.x (dispersion) with MPI ({} processes)",
                            mpi.nprocs
                        ),
                    );
                    tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                        .args(["-np", &mpi.nprocs.to_string()])
                        .arg(&matdyn_exe)
                        .current_dir(&work_path)
                        .stdin(Stdio::piped())
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn()
                        .map_err(|e| {
                            format!("Failed to start mpirun for matdyn.x dispersion: {}", e)
                        })?
                } else {
                    tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                        .current_dir(&work_path)
                        .stdin(Stdio::piped())
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn()
                        .map_err(|e| format!("Failed to start matdyn.x for dispersion: {}", e))?
                }
            } else {
                tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start matdyn.x for dispersion: {}", e))?
            };

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

fn extract_fermi_energy_from_text(output: &str) -> Option<f64> {
    for line in output.lines() {
        if line.contains("the Fermi energy is") {
            if let Some(idx) = line.find("the Fermi energy is") {
                let rest = &line[idx + 19..];
                if let Some(ev_idx) = rest.find("ev") {
                    if let Ok(ef) = rest[..ev_idx].trim().parse::<f64>() {
                        return Some(ef);
                    }
                }
            }
        }
    }
    None
}

fn extract_unmapped_fermi_grid_points(output: &str) -> Option<u32> {
    for line in output.lines() {
        if !line.contains("# of elements that are not done") {
            continue;
        }
        if let Some(raw_count) = line.rsplit(':').next() {
            if let Ok(parsed) = raw_count.trim().parse::<u32>() {
                return Some(parsed);
            }
        }
    }
    None
}

async fn run_hpc_bundle_task(
    app: AppHandle,
    pm: ProcessManager,
    task_id: &str,
    task_kind: &str,
    task_label: &str,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    working_dir: &str,
    command_lines: Vec<String>,
    mut bundle_files: Vec<(String, String)>,
    bundle_copies: Vec<(PathBuf, String)>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
) -> Result<PathBuf, String> {
    // Keep a task-scoped local sync directory to prevent collisions across
    // concurrent HPC submissions from the same wizard.
    let local_sync_dir = PathBuf::from(working_dir).join("hpc_tasks").join(task_id);
    prepare_working_directory(&local_sync_dir, false)?;

    pm.set_hpc_profile_id(task_id, Some(profile.id.clone()))
        .await;
    pm.set_local_sync_dir(task_id, Some(local_sync_dir.to_string_lossy().to_string()))
        .await;

    let bundle_dir = hpc::sync::create_local_bundle_dir(task_id)?;
    let slurm_job_name = format!("qcortado-{}", task_kind);
    let script =
        hpc::slurm::build_slurm_script(&profile, &slurm_job_name, &command_lines, resources);
    let effective_resource_type = match script.effective_resources.resource_type {
        hpc::profile::ResourceType::Cpu => "cpu",
        hpc::profile::ResourceType::Gpu => "gpu",
    };
    pm.set_hpc_resource_type(task_id, Some(effective_resource_type.to_string()))
        .await;
    if let Some(save_spec) = recovery_save.as_ref() {
        let save_json = serde_json::to_value(save_spec).unwrap_or(serde_json::Value::Null);
        pm.set_recovery_save(task_id, Some(save_json)).await;
    }

    if !script.validation.warnings.is_empty() {
        for warning in &script.validation.warnings {
            let line = format!("HPC_WARNING|{}", warning);
            let _ = app.emit(&format!("task-output:{}", task_id), &line);
            pm.append_output(task_id, line).await;
        }
    }
    if !script.validation.errors.is_empty() {
        return Err(script.validation.errors.join(" "));
    }

    bundle_files.push(("run.sbatch".to_string(), script.script.clone()));
    let recovery_metadata = HpcRecoveryJobMetadata {
        schema_version: HPC_RECOVERY_METADATA_VERSION,
        task_id: task_id.to_string(),
        task_kind: task_kind.to_string(),
        label: task_label.to_string(),
        profile_id: profile.id.clone(),
        resource_type: effective_resource_type.to_string(),
        slurm_job_name,
        remote_job_id: None,
        submitted_at: chrono::Utc::now().to_rfc3339(),
        recovery_save,
    };
    let metadata_json = serde_json::to_string_pretty(&recovery_metadata)
        .map_err(|e| format!("Failed to serialize HPC recovery metadata: {}", e))?;
    bundle_files.push((HPC_RECOVERY_METADATA_FILE.to_string(), metadata_json));

    for (source_path, relative_dest) in &bundle_copies {
        hpc::sync::copy_path_into_bundle(&bundle_dir, source_path, relative_dest)?;
    }

    for (path, content) in &bundle_files {
        hpc::sync::write_bundle_text_file(&bundle_dir, path, content)?;
    }

    let run_result = hpc::runner::run_batch_task(
        app,
        pm,
        hpc::runner::HpcBatchRequest {
            task_id: task_id.to_string(),
            task_kind: task_kind.to_string(),
            task_label: task_label.to_string(),
            profile,
            secret,
            slurm_script: script.script.clone(),
            sbatch_preview: script.sbatch_preview,
            bundle_dir: bundle_dir.clone(),
            local_sync_dir: local_sync_dir.clone(),
            cancel_flag,
        },
    )
    .await;

    let _ = std::fs::remove_dir_all(&bundle_dir);
    run_result.map(|_| local_sync_dir)
}

#[derive(Debug, Clone)]
struct ValidatedTransportSource {
    calculation: projects::CalculationRun,
    source_tmp_dir: PathBuf,
    seedname: String,
    reference_fermi_energy_ev: f64,
    use_ws_distance: bool,
    source_win_content: String,
}

fn parse_wannier_win_bool(content: &str, keyword: &str) -> Option<bool> {
    let expected = keyword.trim().to_ascii_lowercase();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('!') || trimmed.starts_with('#') {
            continue;
        }
        let Some((lhs, rhs)) = line.split_once('=') else {
            continue;
        };
        let found_keyword = lhs.trim().to_ascii_lowercase();
        if found_keyword != expected {
            continue;
        }
        let normalized = rhs
            .split(['!', '#'])
            .next()
            .unwrap_or("")
            .trim()
            .trim_matches('.')
            .to_ascii_lowercase();
        return match normalized.as_str() {
            "true" | "t" => Some(true),
            "false" | "f" => Some(false),
            _ => None,
        };
    }
    None
}

fn resolve_transport_seedname(calc: &projects::CalculationRun) -> Option<String> {
    calc.parameters
        .get("seedname")
        .and_then(|value| value.as_str())
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
        .or_else(|| {
            calc.result
                .as_ref()
                .and_then(|result| result.wannier_data.as_ref())
                .and_then(|value| value.get("seedname"))
                .and_then(|value| value.as_str())
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string())
        })
}

fn validate_transport_source(
    app: &AppHandle,
    project_id: &str,
    calc_id: &str,
) -> Result<ValidatedTransportSource, String> {
    let calculation = projects::get_project_calculation(
        app.clone(),
        project_id.to_string(),
        calc_id.to_string(),
    )?;
    if calculation.calc_type.trim().to_ascii_lowercase() != "wannier" {
        return Err(format!(
            "Transport source {} is not a saved Wannier calculation.",
            calc_id
        ));
    }

    let reference_fermi_energy_ev = calculation
        .result
        .as_ref()
        .and_then(|result| result.fermi_energy)
        .ok_or_else(|| {
            "Selected Wannier calculation does not have a saved Fermi energy.".to_string()
        })?;
    let seedname = resolve_transport_seedname(&calculation)
        .ok_or_else(|| "Selected Wannier calculation is missing a saved seedname.".to_string())?;

    let source_tmp_dir = projects::get_projects_dir(app)?
        .join(project_id)
        .join("calculations")
        .join(calc_id)
        .join("tmp");
    if !source_tmp_dir.exists() {
        return Err(format!(
            "Saved Wannier artifacts are missing at {}.",
            source_tmp_dir.display()
        ));
    }

    let win_path = source_tmp_dir.join(format!("{}.win", seedname));
    let chk_path = source_tmp_dir.join(format!("{}.chk", seedname));
    let eig_path = source_tmp_dir.join(format!("{}.eig", seedname));
    for required_path in [&win_path, &chk_path, &eig_path] {
        if !required_path.exists() {
            return Err(format!(
                "Saved Wannier source is missing required artifact {}.",
                required_path.display()
            ));
        }
    }

    let source_win_content = std::fs::read_to_string(&win_path)
        .map_err(|e| format!("Failed to read {}: {}", win_path.display(), e))?;
    let use_ws_distance = parse_wannier_win_bool(&source_win_content, "use_ws_distance")
        .or_else(|| {
            calculation
                .parameters
                .get("use_ws_distance")
                .and_then(|value| value.as_bool())
        })
        .unwrap_or(false);
    if use_ws_distance {
        let wsvec_path = source_tmp_dir.join(format!("{}_wsvec.dat", seedname));
        if !wsvec_path.exists() {
            return Err(format!(
                "Saved Wannier source requires {} but it was not found.",
                wsvec_path.display()
            ));
        }
    }

    Ok(ValidatedTransportSource {
        calculation,
        source_tmp_dir,
        seedname,
        reference_fermi_energy_ev,
        use_ws_distance,
        source_win_content,
    })
}

fn copy_transport_source_file(
    source_root: &Path,
    destination_root: &Path,
    file_name: &str,
) -> Result<(), String> {
    let source_path = source_root.join(file_name);
    if !source_path.exists() {
        return Ok(());
    }
    let destination_path = destination_root.join(file_name);
    if let Some(parent) = destination_path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Failed to create transport staging directory {}: {}",
                parent.display(),
                e
            )
        })?;
    }
    std::fs::copy(&source_path, &destination_path).map_err(|e| {
        format!(
            "Failed to copy {} to {}: {}",
            source_path.display(),
            destination_path.display(),
            e
        )
    })?;
    Ok(())
}

fn stage_transport_source_files(
    source: &ValidatedTransportSource,
    destination_root: &Path,
) -> Result<(), String> {
    std::fs::create_dir_all(destination_root).map_err(|e| {
        format!(
            "Failed to create transport staging directory {}: {}",
            destination_root.display(),
            e
        )
    })?;

    let seedname = &source.seedname;
    let mut files = vec![
        format!("{}.win", seedname),
        format!("{}.chk", seedname),
        format!("{}.eig", seedname),
        format!("{}.nnkp", seedname),
        format!("{}.amn", seedname),
        format!("{}.mmn", seedname),
        format!("{}_hr.dat", seedname),
        format!("{}_wsvec.dat", seedname),
        format!("{}.wout", seedname),
    ];
    if source.use_ws_distance {
        files.push(format!("{}_wsvec.dat", seedname));
    }

    for file_name in files {
        copy_transport_source_file(&source.source_tmp_dir, destination_root, &file_name)?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct ValidatedEpwSources {
    phonon_calculation: projects::CalculationRun,
    phonon_tmp_dir: PathBuf,
    scf_calculation: Option<projects::CalculationRun>,
    scf_save_dir: Option<PathBuf>,
    wannier_calculation: Option<projects::CalculationRun>,
    wannier_tmp_dir: Option<PathBuf>,
    rebuild_wannier_nscf_save: bool,
    manifests: Vec<EpwArtifactManifestEntry>,
    warnings: Vec<String>,
}

fn normalized_rel_path(path: &Path, root: &Path) -> Option<String> {
    let rel = path.strip_prefix(root).ok()?;
    let mut components: Vec<String> = Vec::new();
    for component in rel.components() {
        let std::path::Component::Normal(segment) = component else {
            return None;
        };
        let segment = segment.to_string_lossy().to_string();
        if segment.is_empty() {
            return None;
        }
        components.push(segment);
    }
    if components.is_empty() {
        None
    } else {
        Some(components.join("/"))
    }
}

fn collect_epw_manifest_entries<F>(
    source_root: &Path,
    source_calc_id: &str,
    source_calc_type: &str,
    mut include: F,
) -> Result<Vec<EpwArtifactManifestEntry>, String>
where
    F: FnMut(&str, &str) -> bool,
{
    if !source_root.exists() || !source_root.is_dir() {
        return Ok(Vec::new());
    }

    let mut queue: VecDeque<PathBuf> = VecDeque::new();
    queue.push_back(source_root.to_path_buf());

    let mut entries: Vec<EpwArtifactManifestEntry> = Vec::new();
    while let Some(current_dir) = queue.pop_front() {
        let dir_entries = std::fs::read_dir(&current_dir).map_err(|e| {
            format!(
                "Failed to inspect EPW prerequisite directory {}: {}",
                current_dir.display(),
                e
            )
        })?;

        for entry in dir_entries {
            let entry = entry.map_err(|e| e.to_string())?;
            let path = entry.path();
            let file_type = entry.file_type().map_err(|e| e.to_string())?;
            if file_type.is_dir() {
                queue.push_back(path);
                continue;
            }
            if !file_type.is_file() {
                continue;
            }

            let Some(rel_path) = normalized_rel_path(&path, source_root) else {
                continue;
            };
            let lower_rel = rel_path.to_ascii_lowercase();
            let file_name = path
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            if !include(&lower_rel, &file_name) {
                continue;
            }
            let metadata = entry.metadata().map_err(|e| {
                format!(
                    "Failed to inspect prerequisite file {}: {}",
                    path.display(),
                    e
                )
            })?;
            entries.push(EpwArtifactManifestEntry {
                source_calc_id: source_calc_id.to_string(),
                source_calc_type: source_calc_type.to_string(),
                rel_path,
                size_bytes: metadata.len(),
            });
        }
    }

    entries.sort_by(|a, b| a.rel_path.cmp(&b.rel_path));
    Ok(entries)
}

fn resolve_epw_stage_subdir(raw: &str, fallback: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return fallback.to_string();
    }
    let normalized = trimmed
        .trim_start_matches("./")
        .trim_matches('/')
        .to_string();
    if normalized.is_empty() {
        return fallback.to_string();
    }
    let path = Path::new(&normalized);
    if path
        .components()
        .all(|component| matches!(component, std::path::Component::Normal(_)))
    {
        normalized
    } else {
        fallback.to_string()
    }
}

fn copy_epw_source_manifest_entry(
    source_root: &Path,
    destination_root: &Path,
    destination_subdir: &str,
    rel_path: &str,
) -> Result<(), String> {
    let source_path = source_root.join(rel_path);
    if !source_path.exists() || !source_path.is_file() {
        return Err(format!(
            "EPW prerequisite artifact missing: {}",
            source_path.display()
        ));
    }
    let destination_path = destination_root.join(destination_subdir).join(rel_path);
    if let Some(parent) = destination_path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Failed to create EPW staging directory {}: {}",
                parent.display(),
                e
            )
        })?;
    }
    std::fs::copy(&source_path, &destination_path).map_err(|e| {
        format!(
            "Failed to stage EPW prerequisite {} -> {}: {}",
            source_path.display(),
            destination_path.display(),
            e
        )
    })?;
    Ok(())
}

fn copy_epw_source_file_to(
    source_root: &Path,
    destination_root: &Path,
    source_rel_path: &str,
    destination_rel_path: &str,
) -> Result<(), String> {
    let source_path = source_root.join(source_rel_path);
    if !source_path.exists() || !source_path.is_file() {
        return Err(format!(
            "EPW prerequisite artifact missing: {}",
            source_path.display()
        ));
    }
    let destination_path = destination_root.join(destination_rel_path);
    if let Some(parent) = destination_path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Failed to create EPW staging directory {}: {}",
                parent.display(),
                e
            )
        })?;
    }
    std::fs::copy(&source_path, &destination_path).map_err(|e| {
        format!(
            "Failed to stage EPW prerequisite {} -> {}: {}",
            source_path.display(),
            destination_path.display(),
            e
        )
    })?;
    Ok(())
}

fn find_epw_phsave_dir(source_root: &Path, prefix: &str) -> Option<PathBuf> {
    let expected = format!("{}.phsave", prefix.trim()).to_ascii_lowercase();
    let mut fallback: Option<PathBuf> = None;
    let mut queue: VecDeque<PathBuf> = VecDeque::new();
    queue.push_back(source_root.to_path_buf());
    while let Some(current) = queue.pop_front() {
        let Ok(entries) = std::fs::read_dir(&current) else {
            continue;
        };
        for entry in entries.flatten() {
            let Ok(file_type) = entry.file_type() else {
                continue;
            };
            if !file_type.is_dir() {
                continue;
            }
            let path = entry.path();
            let name = path
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            if name == expected {
                return Some(path);
            }
            if name.ends_with(".phsave") && fallback.is_none() {
                fallback = Some(path.clone());
            }
            queue.push_back(path);
        }
    }
    fallback
}

fn stage_epw_phonon_save_directory(
    source: &ValidatedEpwSources,
    destination_root: &Path,
    config: &EpwCalculationConfig,
) -> Result<(), String> {
    let prefix = config.input.prefix.trim();
    let save_subdir = resolve_epw_stage_subdir(&config.input.dvscf_dir, "save");
    let dyn_re = Regex::new(r"(?i)(?:^|/)(?:.+\.)?dyn(\d+)(\.xml)?$").unwrap();
    let q_dir_re = Regex::new(r"(?i)(?:^|/)q_(\d+)(?:/|$)").unwrap();
    let dvscf_q_re = Regex::new(r"(?i)dvscf_q(\d+)").unwrap();

    for entry in source
        .manifests
        .iter()
        .filter(|entry| entry.source_calc_type == "phonon")
    {
        let lower = entry.rel_path.to_ascii_lowercase();
        if let Some(captures) = dyn_re.captures(&entry.rel_path) {
            let q_index = captures.get(1).map(|value| value.as_str()).unwrap_or("1");
            let extension = captures.get(2).map(|value| value.as_str()).unwrap_or("");
            copy_epw_source_file_to(
                &source.phonon_tmp_dir,
                destination_root,
                &entry.rel_path,
                &format!("{}/{}.dyn_q{}{}", save_subdir, prefix, q_index, extension),
            )?;
            continue;
        }

        if lower.contains("dvscf") {
            let q_index = q_dir_re
                .captures(&entry.rel_path)
                .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
                .or_else(|| {
                    dvscf_q_re.captures(&entry.rel_path).and_then(|captures| {
                        captures.get(1).map(|value| value.as_str().to_string())
                    })
                })
                .unwrap_or_else(|| "1".to_string());
            copy_epw_source_file_to(
                &source.phonon_tmp_dir,
                destination_root,
                &entry.rel_path,
                &format!("{}/{}.dvscf_q{}", save_subdir, prefix, q_index),
            )?;
            continue;
        }

        let file_name = Path::new(&entry.rel_path)
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        if matches!(file_name.as_str(), "ifc.q2r" | "ifc.q2r.xml") {
            copy_epw_source_file_to(
                &source.phonon_tmp_dir,
                destination_root,
                &entry.rel_path,
                &format!("{}/{}", save_subdir, file_name),
            )?;
        }
    }

    if let Some(phsave_dir) = find_epw_phsave_dir(&source.phonon_tmp_dir, prefix) {
        let destination = destination_root
            .join(&save_subdir)
            .join(format!("{}.phsave", prefix));
        projects::copy_dir_contents(&phsave_dir, &destination).map_err(|e| {
            format!(
                "Failed to stage EPW phsave directory {} -> {}: {}",
                phsave_dir.display(),
                destination.display(),
                e
            )
        })?;
    }

    Ok(())
}

fn stage_epw_scf_save_directory(
    source: &ValidatedEpwSources,
    destination_root: &Path,
    config: &EpwCalculationConfig,
) -> Result<(), String> {
    let Some(source_save_dir) = source.scf_save_dir.as_ref() else {
        return Ok(());
    };

    let outdir_subdir = resolve_epw_stage_subdir(&config.input.outdir, "tmp");
    let expected_save_dir_name = format!("{}.save", config.input.prefix.trim());
    let destination_save_dir = destination_root
        .join(&outdir_subdir)
        .join(&expected_save_dir_name);
    if let Some(parent) = destination_save_dir.parent() {
        std::fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Failed to create EPW outdir staging path {}: {}",
                parent.display(),
                e
            )
        })?;
    }

    let source_buf = source_save_dir.to_path_buf();
    let destination_buf = destination_save_dir.clone();
    projects::copy_dir_contents(&source_buf, &destination_buf).map_err(|e| {
        format!(
            "Failed to stage SCF save directory {} -> {}: {}",
            source_save_dir.display(),
            destination_save_dir.display(),
            e
        )
    })?;

    let schema_path = destination_save_dir.join("data-file-schema.xml");
    if !schema_path.exists() {
        return Err(format!(
            "Staged SCF save directory is missing {}.",
            schema_path.display()
        ));
    }

    Ok(())
}

fn stage_epw_source_files(
    source: &ValidatedEpwSources,
    destination_root: &Path,
    config: &EpwCalculationConfig,
) -> Result<(), String> {
    std::fs::create_dir_all(destination_root).map_err(|e| {
        format!(
            "Failed to create EPW staging directory {}: {}",
            destination_root.display(),
            e
        )
    })?;
    stage_epw_scf_save_directory(source, destination_root, config)?;

    let optional_wannier_subdir = config
        .input
        .wannier_dir
        .trim()
        .trim_start_matches("./")
        .trim_matches('/')
        .to_string();

    for entry in &source.manifests {
        match entry.source_calc_type.as_str() {
            "phonon" => {}
            "wannier" => {
                let Some(wannier_root) = source.wannier_tmp_dir.as_ref() else {
                    return Err(
                        "Internal EPW staging error: missing Wannier source root.".to_string()
                    );
                };
                // EPW 6.0 expects Wannier interface files in the run root (no `wannier_dir` keyword).
                copy_epw_source_manifest_entry(
                    wannier_root,
                    destination_root,
                    "",
                    &entry.rel_path,
                )?;
                // Mirror into the optional user subdir as a compatibility convenience.
                if !optional_wannier_subdir.is_empty() {
                    copy_epw_source_manifest_entry(
                        wannier_root,
                        destination_root,
                        &optional_wannier_subdir,
                        &entry.rel_path,
                    )?;
                }
            }
            _ => {}
        }
    }
    stage_epw_phonon_save_directory(source, destination_root, config)?;
    Ok(())
}

fn extract_bool_param(parameters: &serde_json::Value, key: &str) -> Option<bool> {
    parameters.get(key).and_then(|value| value.as_bool())
}

fn extract_u32_param(parameters: &serde_json::Value, key: &str) -> Option<u32> {
    parameters
        .get(key)
        .and_then(|value| value.as_u64())
        .and_then(|value| u32::try_from(value).ok())
}

fn extract_u32_triplet_param(parameters: &serde_json::Value, key: &str) -> Option<[u32; 3]> {
    let values = parameters.get(key)?.as_array()?;
    if values.len() != 3 {
        return None;
    }
    let x = u32::try_from(values[0].as_u64()?).ok()?;
    let y = u32::try_from(values[1].as_u64()?).ok()?;
    let z = u32::try_from(values[2].as_u64()?).ok()?;
    if x == 0 || y == 0 || z == 0 {
        return None;
    }
    Some([x, y, z])
}

fn extract_string_param(parameters: &serde_json::Value, key: &str) -> Option<String> {
    parameters
        .get(key)
        .and_then(|value| value.as_str())
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
}

fn validate_epw_prerequisites_detailed(
    app: &AppHandle,
    config: &EpwCalculationConfig,
) -> (EpwPrerequisiteValidation, Option<ValidatedEpwSources>) {
    let mut validation = EpwPrerequisiteValidation {
        ok: false,
        errors: Vec::new(),
        warnings: Vec::new(),
        remediation_hints: Vec::new(),
        manifests: Vec::new(),
    };

    if let Err(err) = validate_epw_config(config) {
        validation.errors.push(err);
        return (validation, None);
    }

    let project_id = config.project_id.trim();
    let phonon_calc_id = config.source_phonon_calc_id.trim();
    let wannier_calc_id = config
        .source_wannier_calc_id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string());
    let explicit_scf_calc_id = config
        .source_scf_calc_id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string());

    let phonon_calculation = match projects::get_project_calculation(
        app.clone(),
        project_id.to_string(),
        phonon_calc_id.to_string(),
    ) {
        Ok(calc) => calc,
        Err(err) => {
            validation.errors.push(format!(
                "Unable to load phonon source calculation {}: {}",
                phonon_calc_id, err
            ));
            return (validation, None);
        }
    };
    if phonon_calculation.calc_type.trim().to_ascii_lowercase() != "phonon" {
        validation.errors.push(format!(
            "EPW source {} is not a phonon calculation.",
            phonon_calc_id
        ));
        return (validation, None);
    }
    if phonon_calculation.completed_at.is_none() {
        validation.errors.push(format!(
            "Phonon source {} is incomplete. Wait for completion before running EPW.",
            phonon_calc_id
        ));
    }
    if phonon_calculation
        .result
        .as_ref()
        .map(|result| result.converged)
        .unwrap_or(false)
        == false
    {
        validation.warnings.push(format!(
            "Phonon source {} is not marked converged. EPW may fail or produce inconsistent data.",
            phonon_calc_id
        ));
    }

    let phonon_tmp_dir = match projects::get_projects_dir(app) {
        Ok(root) => root
            .join(project_id)
            .join("calculations")
            .join(phonon_calc_id)
            .join("tmp"),
        Err(err) => {
            validation
                .errors
                .push(format!("Failed to resolve projects directory: {}", err));
            return (validation, None);
        }
    };
    if !phonon_tmp_dir.exists() {
        validation.errors.push(format!(
            "Saved phonon artifacts are missing at {}.",
            phonon_tmp_dir.display()
        ));
        validation
            .remediation_hints
            .push("Re-run phonons with save policy `epw-ready` or `full-archive`.".to_string());
    }

    let phonon_manifest = collect_epw_manifest_entries(
        &phonon_tmp_dir,
        phonon_calc_id,
        "phonon",
        |lower_rel, file_name| {
            file_name.starts_with("dyn")
                || file_name.contains(".dyn")
                || file_name.contains("dvscf")
                || file_name.ends_with(".ukk")
                || matches!(
                    file_name,
                    "force_constants"
                        | "ph.in"
                        | "ph.out"
                        | "q2r.in"
                        | "q2r.out"
                        | "matdyn_dos.in"
                        | "matdyn_dos.out"
                        | "matdyn_bands.in"
                        | "matdyn_bands.out"
                        | "phonon_freq"
                        | "phonon_freq.gp"
                        | "phonon_dos"
                        | "ifc.q2r"
                        | "ifc.q2r.xml"
                )
                || (lower_rel.contains("_ph0/") && file_name.ends_with(".xml"))
        },
    )
    .unwrap_or_else(|err| {
        validation.errors.push(err);
        Vec::new()
    });
    let has_dyn = phonon_manifest.iter().any(|entry| {
        let file_name = Path::new(&entry.rel_path)
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        file_name.starts_with("dyn") || file_name.contains(".dyn")
    });
    let has_dvscf = phonon_manifest
        .iter()
        .any(|entry| entry.rel_path.to_ascii_lowercase().contains("dvscf"));
    let has_ph_metadata = phonon_manifest.iter().any(|entry| {
        let lower = entry.rel_path.to_ascii_lowercase();
        lower.ends_with("status_run.xml")
            || lower.ends_with("control_ph.xml")
            || lower.ends_with("patterns.xml")
    });
    if !has_dyn {
        validation
            .errors
            .push("Phonon source is missing required `dyn*` artifacts for EPW.".to_string());
        validation.remediation_hints.push(
            "Re-run phonons and keep `dyn*` outputs using save policy `epw-ready`.".to_string(),
        );
    }
    if !has_dvscf {
        validation.errors.push(
            "Phonon source is missing `dvscf*` perturbation artifacts required by EPW.".to_string(),
        );
        validation
            .remediation_hints
            .push("Enable EPW preparation in phonons (`fildvscf`) and rerun with save policy `epw-ready`.".to_string());
    }
    if !has_ph_metadata {
        validation.warnings.push(
            "Phonon source metadata XML files were not found (`status_run.xml`/`control_ph.xml`)."
                .to_string(),
        );
    }

    let phonon_source_scf = extract_string_param(&phonon_calculation.parameters, "source_scf_id");
    let mut resolved_scf_calc_id = explicit_scf_calc_id
        .clone()
        .or_else(|| phonon_source_scf.clone());
    let phonon_q_grid = phonon_calculation
        .parameters
        .get("q_grid")
        .and_then(|value| value.as_array())
        .and_then(|values| {
            if values.len() != 3 {
                return None;
            }
            let q1 = values[0].as_u64()?;
            let q2 = values[1].as_u64()?;
            let q3 = values[2].as_u64()?;
            Some([q1 as u32, q2 as u32, q3 as u32])
        });
    let configured_coarse_q = epw_coarse_q_mesh(&config.input);
    if let Some(saved_q) = phonon_q_grid {
        if saved_q != configured_coarse_q {
            validation.errors.push(format!(
                "Configured EPW coarse q-mesh {}x{}x{} differs from saved phonon q-grid {}x{}x{}.",
                configured_coarse_q[0],
                configured_coarse_q[1],
                configured_coarse_q[2],
                saved_q[0],
                saved_q[1],
                saved_q[2]
            ));
            validation.remediation_hints.push(
                "Select a phonon source with recorded q-grid metadata; EPW coarse q-mesh is read from the phonon run."
                    .to_string(),
            );
        }
    } else {
        validation.errors.push(
            "Selected phonon source does not record its q-grid, so EPW cannot infer the required coarse q-mesh."
                .to_string(),
        );
        validation.remediation_hints.push(
            "Re-run the phonon workflow so QCortado can save the phonon q-grid metadata."
                .to_string(),
        );
    }

    let phonon_prefix = extract_string_param(&phonon_calculation.parameters, "prefix")
        .unwrap_or_else(|| "qcortado_scf".to_string());
    if config.input.prefix.trim() != phonon_prefix {
        validation.warnings.push(format!(
            "EPW prefix '{}' differs from phonon source prefix '{}'.",
            config.input.prefix.trim(),
            phonon_prefix
        ));
    }

    let expected_save_dir_name = format!("{}.save", config.input.prefix.trim());
    let effective_wannierize = resolve_epw_wannierization_plan(config)
        .map(|plan| plan.wannierize)
        .unwrap_or(config.input.wannierize);
    let mut manifests = phonon_manifest;
    let mut scf_calculation: Option<projects::CalculationRun> = None;
    let mut scf_save_dir: Option<PathBuf> = None;
    let mut wannier_calculation: Option<projects::CalculationRun> = None;
    let mut wannier_tmp_dir: Option<PathBuf> = None;
    let mut wannier_nscf_input_path: Option<PathBuf> = None;
    let mut bundled_wannier_save_dir: Option<PathBuf> = None;
    let mut rebuild_wannier_nscf_save = false;

    if let Some(wannier_id) = wannier_calc_id.as_ref() {
        let source = match projects::get_project_calculation(
            app.clone(),
            project_id.to_string(),
            wannier_id.to_string(),
        ) {
            Ok(calc) => calc,
            Err(err) => {
                validation.errors.push(format!(
                    "Unable to load Wannier source calculation {}: {}",
                    wannier_id, err
                ));
                return (validation, None);
            }
        };
        if source.calc_type.trim().to_ascii_lowercase() != "wannier" {
            validation.errors.push(format!(
                "EPW source {} is not a Wannier calculation.",
                wannier_id
            ));
            return (validation, None);
        }
        let seedname =
            resolve_transport_seedname(&source).unwrap_or_else(|| "qcortado_wannier".to_string());
        let source_tmp_dir = projects::get_projects_dir(app)
            .map(|root| {
                root.join(project_id)
                    .join("calculations")
                    .join(wannier_id)
                    .join("tmp")
            })
            .unwrap_or_else(|_| PathBuf::from(""));
        if !source_tmp_dir.exists() {
            validation.errors.push(format!(
                "Saved Wannier artifacts are missing at {}.",
                source_tmp_dir.display()
            ));
            validation
                .remediation_hints
                .push("Re-run and save the Wannier workflow before running EPW.".to_string());
        }
        let candidate_nscf_input = source_tmp_dir.join("nscf.in");
        if candidate_nscf_input.exists() {
            wannier_nscf_input_path = Some(candidate_nscf_input);
        }
        if let Some(wannier_k_grid) = extract_u32_triplet_param(&source.parameters, "k_grid") {
            let configured_coarse_k = epw_coarse_k_mesh(&config.input);
            if wannier_k_grid != configured_coarse_k {
                let message = format!(
                    "Configured EPW coarse k-mesh {}x{}x{} differs from Wannier source NSCF k-grid {}x{}x{}.",
                    configured_coarse_k[0],
                    configured_coarse_k[1],
                    configured_coarse_k[2],
                    wannier_k_grid[0],
                    wannier_k_grid[1],
                    wannier_k_grid[2]
                );
                if effective_wannierize {
                    validation.errors.push(message);
                    validation.remediation_hints.push(
                        "Set EPW k-mesh to match the selected Wannier source k-grid (or re-run Wannier source on the desired coarse mesh)."
                            .to_string(),
                    );
                } else {
                    validation.warnings.push(message);
                }
            }
        }
        let required_wannier_files = [
            format!("{}.win", seedname),
            format!("{}.chk", seedname),
            format!("{}.eig", seedname),
        ];
        for file_name in required_wannier_files {
            if !source_tmp_dir.join(&file_name).exists() {
                validation.errors.push(format!(
                    "Wannier source is missing required artifact {}.",
                    file_name
                ));
            }
        }
        let wannier_manifest = collect_epw_manifest_entries(
            &source_tmp_dir,
            wannier_id,
            "wannier",
            |lower_rel, file_name| {
                file_name.ends_with(".win")
                    || file_name.ends_with(".chk")
                    || file_name.ends_with(".eig")
                    || file_name.ends_with(".ukk")
                    || file_name.ends_with(".nnkp")
                    || file_name.ends_with(".amn")
                    || file_name.ends_with(".mmn")
                    || file_name.ends_with("_hr.dat")
                    || file_name.ends_with("_wsvec.dat")
                    || file_name.ends_with(".wout")
                    || file_name == "nscf.in"
                    || file_name == "nscf.out"
                    || lower_rel.ends_with(".xml")
            },
        )
        .unwrap_or_else(|err| {
            validation.errors.push(err);
            Vec::new()
        });
        manifests.extend(wannier_manifest);
        if scf_save_dir.is_none() {
            let bundled_save_dir =
                find_preferred_wannier_save_directory(&source_tmp_dir, &expected_save_dir_name);
            if let Some(found_save_dir) = bundled_save_dir {
                let schema_path = found_save_dir.join("data-file-schema.xml");
                if schema_path.exists() {
                    validation.warnings.push(
                        "Using QE `.save` payload bundled with the selected Wannier source."
                            .to_string(),
                    );
                    bundled_wannier_save_dir = Some(found_save_dir.clone());
                    scf_save_dir = Some(found_save_dir);
                }
            }
        }

        let wannier_source_scf = extract_string_param(&source.parameters, "source_scf_id");
        if let (Some(phonon_scf), Some(wannier_scf)) =
            (phonon_source_scf.as_ref(), wannier_source_scf.as_ref())
        {
            if phonon_scf != wannier_scf {
                validation.errors.push(format!(
                    "Phonon source SCF ({}) and Wannier source SCF ({}) do not match.",
                    phonon_scf, wannier_scf
                ));
                validation.remediation_hints.push(
                    "Regenerate phonon/Wannier from the same SCF reference before EPW.".to_string(),
                );
            }
        }
        if resolved_scf_calc_id.is_none() {
            resolved_scf_calc_id = wannier_source_scf.clone();
        }

        // Basic spin/SOC consistency checks when metadata is available.
        let phonon_nspin = extract_u32_param(&phonon_calculation.parameters, "nspin");
        let wannier_nspin = extract_u32_param(&source.parameters, "nspin");
        if let (Some(left), Some(right)) = (phonon_nspin, wannier_nspin) {
            if left != right {
                validation.errors.push(format!(
                    "Spin-channel mismatch between phonon source (nspin={}) and Wannier source (nspin={}).",
                    left, right
                ));
            }
        }
        let phonon_soc = extract_bool_param(&phonon_calculation.parameters, "lspinorb");
        let wannier_soc = extract_bool_param(&source.parameters, "lspinorb");
        if let (Some(left), Some(right)) = (phonon_soc, wannier_soc) {
            if left != right {
                validation.errors.push(
                    "SOC mismatch between phonon and Wannier prerequisites (`lspinorb`)."
                        .to_string(),
                );
            }
        }

        wannier_calculation = Some(source);
        wannier_tmp_dir = Some(source_tmp_dir);
    } else {
        validation.warnings.push(
            "No Wannier source was provided. EPW will run with phonon-only prerequisites."
                .to_string(),
        );
    }

    if let Some(explicit_scf) = explicit_scf_calc_id.as_ref() {
        if let Some(phonon_scf) = phonon_source_scf.as_ref() {
            if phonon_scf != explicit_scf {
                validation.errors.push(format!(
                    "Configured source_scf_calc_id ({}) does not match phonon source ({})",
                    explicit_scf, phonon_scf
                ));
            }
        }
    }

    if scf_save_dir.is_none() {
        if let Some(scf_calc_id) = resolved_scf_calc_id.as_ref() {
            match projects::get_project_calculation(
                app.clone(),
                project_id.to_string(),
                scf_calc_id.to_string(),
            ) {
                Ok(source) => {
                    let source_tmp_dir = projects::get_projects_dir(app)
                        .map(|root| {
                            root.join(project_id)
                                .join("calculations")
                                .join(scf_calc_id)
                                .join("tmp")
                        })
                        .unwrap_or_else(|_| PathBuf::from(""));
                    if !source_tmp_dir.exists() {
                        validation
                            .errors
                            .push(missing_scf_tmp_error(&source_tmp_dir));
                    } else {
                        let found_save_dir = find_qe_save_directory(
                            &source_tmp_dir,
                            4,
                            Some(&expected_save_dir_name),
                        )
                        .or_else(|| find_qe_save_directory(&source_tmp_dir, 4, None));
                        if let Some(found_save_dir) = found_save_dir {
                            let schema_path = found_save_dir.join("data-file-schema.xml");
                            if !schema_path.exists() {
                                validation.errors.push(format!(
                                    "Resolved SCF save payload is missing {}.",
                                    schema_path.display()
                                ));
                                validation.remediation_hints.push(
                                "Re-run SCF/NSCF source and preserve the full `.save` directory."
                                    .to_string(),
                            );
                            } else {
                                if found_save_dir
                                    .file_name()
                                    .and_then(|value| value.to_str())
                                    .map(|value| value != expected_save_dir_name)
                                    .unwrap_or(true)
                                {
                                    validation.warnings.push(format!(
                                    "Resolved SCF save directory is '{}', while EPW prefix expects '{}'.",
                                    found_save_dir
                                        .file_name()
                                        .and_then(|value| value.to_str())
                                        .unwrap_or("<unknown>"),
                                    expected_save_dir_name
                                ));
                                }
                                scf_save_dir = Some(found_save_dir);
                                scf_calculation = Some(source);
                            }
                        } else {
                            validation.errors.push(format!(
                                "Resolved SCF source {} does not contain a QE `.save` directory.",
                                scf_calc_id
                            ));
                            validation.remediation_hints.push(
                                "Re-run SCF and keep the `.save` payload, then retry EPW."
                                    .to_string(),
                            );
                        }
                    }
                }
                Err(err) => {
                    validation.errors.push(format!(
                        "Unable to load resolved SCF source calculation {}: {}",
                        scf_calc_id, err
                    ));
                }
            }
        }
    }

    if scf_save_dir.is_none() {
        let fallback_save_dir =
            find_qe_save_directory(&phonon_tmp_dir, 4, Some(&expected_save_dir_name))
                .or_else(|| find_qe_save_directory(&phonon_tmp_dir, 4, None));
        if let Some(found_save_dir) = fallback_save_dir {
            let schema_path = found_save_dir.join("data-file-schema.xml");
            if schema_path.exists() {
                validation.warnings.push(
                    "Using QE `.save` payload discovered in phonon artifacts because resolved SCF artifacts were unavailable."
                        .to_string(),
                );
                scf_save_dir = Some(found_save_dir);
            }
        }
    }

    if effective_wannierize && wannier_calculation.is_some() && bundled_wannier_save_dir.is_none() {
        match (wannier_nscf_input_path.as_ref(), scf_save_dir.as_ref()) {
            (Some(_), Some(_)) => {
                rebuild_wannier_nscf_save = true;
                validation.warnings.push(
                    "Selected Wannier source does not retain its coarse-grid QE `.save`; QCortado will rebuild it from the saved `nscf.in` before running EPW."
                        .to_string(),
                );
            }
            (None, _) => {
                validation.errors.push(
                    "Selected Wannier source does not contain a bundled coarse-grid QE `.save` or a saved `nscf.in` for reconstruction."
                        .to_string(),
                );
                validation.remediation_hints.push(
                    "Re-run the Wannier workflow after updating QCortado so the NSCF `.save` is preserved, then retry EPW."
                        .to_string(),
                );
            }
            (Some(_), None) => {}
        }
    }

    if scf_save_dir.is_none() {
        validation.errors.push(
            "EPW requires a staged QE `.save` directory (`outdir/prefix.save`) with data-file-schema.xml."
                .to_string(),
        );
        validation.remediation_hints.push(
            "Ensure the source SCF calculation is saved locally with its `.save` directory, or run 'Download Full' for the SCF/phonon source and retry."
                .to_string(),
        );
    }

    if let Ok(wannierization_plan) = resolve_epw_wannierization_plan(config) {
        if !wannierization_plan.wannierize
            && !epw_manifest_contains_file(&manifests, &wannierization_plan.filukk)
        {
            if wannierization_plan.explicit_wannierize_override {
                validation.errors.push(format!(
                    "EPW is configured with `wannierize = .false.` but '{}' was not found among staged prerequisites.",
                    wannierization_plan.filukk
                ));
                validation.remediation_hints.push(
                    "Set `wannierize = .true.` (or remove the override) for first-pass EPW runs, or provide the required `.ukk` file."
                        .to_string(),
                );
            } else {
                validation.warnings.push(format!(
                    "No '{}' prerequisite was found. QCortado will auto-enable `wannierize = .true.` for this run.",
                    wannierization_plan.filukk
                ));
            }
        }
    }

    let mut hint_set: HashSet<String> = validation.remediation_hints.iter().cloned().collect();
    if !validation.errors.is_empty()
        && hint_set.insert(
            "If artifacts were trimmed after an HPC run, use 'Download Full' then retry EPW."
                .to_string(),
        )
    {
        validation.remediation_hints.push(
            "If artifacts were trimmed after an HPC run, use 'Download Full' then retry EPW."
                .to_string(),
        );
    }

    manifests.sort_by(|a, b| {
        a.source_calc_type
            .cmp(&b.source_calc_type)
            .then(a.source_calc_id.cmp(&b.source_calc_id))
            .then(a.rel_path.cmp(&b.rel_path))
    });
    validation.manifests = manifests.clone();
    validation.ok = validation.errors.is_empty();

    if !validation.ok {
        return (validation, None);
    }

    (
        validation.clone(),
        Some(ValidatedEpwSources {
            phonon_calculation,
            phonon_tmp_dir,
            scf_calculation,
            scf_save_dir,
            wannier_calculation,
            wannier_tmp_dir,
            rebuild_wannier_nscf_save,
            manifests,
            warnings: validation.warnings,
        }),
    )
}

fn resolve_hpc_execution(
    state: &AppState,
    execution_target: Option<&hpc::profile::ExecutionTarget>,
) -> Option<hpc::profile::HpcExecutionTarget> {
    if effective_execution_mode(state, execution_target) != hpc::profile::ExecutionMode::Hpc {
        return None;
    }
    execution_target
        .and_then(|target| target.hpc.clone())
        .or_else(|| Some(hpc::profile::HpcExecutionTarget::default()))
}

#[tauri::command]
fn validate_epw_prerequisites(
    app: AppHandle,
    config: EpwCalculationConfig,
) -> Result<EpwPrerequisiteValidation, String> {
    let (validation, _) = validate_epw_prerequisites_detailed(&app, &config);
    Ok(validation)
}

#[tauri::command]
fn build_epw_input_preview(config: EpwCalculationConfig) -> Result<EpwInputPreviewResult, String> {
    qe::build_epw_input_preview(&config)
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
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    // Local runs remain serialized. HPC submissions can run concurrently.
    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("scf".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();
        tokio::spawn(async move {
            let result = run_scf_hpc_background(
                app_handle.clone(),
                &tid,
                calculation,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
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
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

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

async fn run_scf_hpc_background(
    app: AppHandle,
    task_id: &str,
    calculation: QECalculation,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<QEResult, String> {
    let mut remote_calculation = calculation;
    remote_calculation.pseudo_dir = profile.remote_pseudo_dir.clone();
    let input = generate_pw_input(&remote_calculation);
    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let commands = vec![
        "cd \"$SLURM_SUBMIT_DIR\"".to_string(),
        format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)),
        build_hpc_qe_input_command(&profile, resource_type, "pw.x", None, "pw.in", "pw.out"),
    ];

    let work_path = run_hpc_bundle_task(
        app,
        pm.clone(),
        task_id,
        "scf",
        "SCF",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        vec![("pw.in".to_string(), input)],
        Vec::new(),
        recovery_save,
        cancel_flag,
    )
    .await?;

    let pw_out_path = work_path.join("pw.out");
    let slurm_out_path = work_path.join("slurm.out");
    let output_text = if pw_out_path.exists() {
        std::fs::read_to_string(&pw_out_path)
            .map_err(|e| format!("Failed to read {}: {}", pw_out_path.display(), e))?
    } else if slurm_out_path.exists() {
        std::fs::read_to_string(&slurm_out_path)
            .map_err(|e| format!("Failed to read {}: {}", slurm_out_path.display(), e))?
    } else {
        return Err("HPC run completed but no output file was downloaded.".to_string());
    };

    Ok(parse_pw_output(&output_text))
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

    prepare_working_directory(&work_path, false)?;

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
                .args(["-pd", ".true."])
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
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("bands".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();

        tokio::spawn(async move {
            let result = run_bands_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
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
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

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
async fn run_bands_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: BandsCalculationConfig,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<BandData, String> {
    let pipeline_start = std::time::Instant::now();
    let mut bands_calc = config.base_calculation.clone();
    bands_calc.pseudo_dir = profile.remote_pseudo_dir.clone();
    bands_calc.calculation = qe::CalculationType::Bands;
    if bands_calc.verbosity.is_none() {
        bands_calc.verbosity = Some("high".to_string());
    }
    let band_path: Vec<qe::BandPathPoint> = config
        .k_path
        .iter()
        .map(|point| qe::BandPathPoint {
            k: point.coords,
            npoints: point.npoints,
            label: Some(point.label.clone()),
        })
        .collect();
    bands_calc.kpoints = qe::KPoints::CrystalB { path: band_path };

    let mut nscf_input = generate_pw_input(&bands_calc);
    if let Some(nbnd) = config.nbnd {
        nscf_input = insert_system_namelist_line(&nscf_input, &format!("nbnd = {},", nbnd))?;
    }

    let bands_x_options = config.bands_x.as_ref();
    let bands_filband = bands_x_options
        .and_then(|opts| opts.filband.as_deref())
        .map(|raw| sanitize_output_filename(raw, "bands.dat"))
        .unwrap_or_else(|| "bands.dat".to_string());
    let bands_lsym = bands_x_options.and_then(|opts| opts.lsym).unwrap_or(true);
    let bands_no_overlap = bands_x_options
        .and_then(|opts| opts.no_overlap)
        .unwrap_or(true);
    let bands_x_config = BandsXConfig {
        prefix: bands_calc.prefix.clone(),
        outdir: bands_calc.outdir.clone(),
        filband: bands_filband.clone(),
        lsym: bands_lsym,
        no_overlap: bands_no_overlap,
    };
    let bands_x_input = generate_bands_x_input(&bands_x_config);

    let projections_enabled = config
        .projections
        .as_ref()
        .map(|entry| entry.enabled)
        .unwrap_or(false);
    let config_line = format!(
        "[QCortado] Bands HPC config: k_path_points={}, projections_enabled={}, bands_filband={}.",
        config.k_path.len(),
        projections_enabled,
        bands_filband
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &config_line);
    pm.append_output(task_id, config_line).await;

    let projection_options = config.projections.as_ref();
    let projection_file = projection_options
        .and_then(|opts| opts.filproj.as_deref())
        .map(|raw| sanitize_output_filename(raw, "bands.projwfc.dat"))
        .unwrap_or_else(|| "bands.projwfc.dat".to_string());
    let projwfc_config = ProjwfcConfig {
        prefix: bands_calc.prefix.clone(),
        outdir: bands_calc.outdir.clone(),
        filproj: projection_file,
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

    let dependency_stage = resolve_hpc_scf_dependency_stage(
        &app,
        config.project_id.as_deref(),
        config.scf_calc_id.as_deref(),
    )?;
    let has_remote_hydration = !dependency_stage.remote_hydration_commands.is_empty();
    if has_remote_hydration {
        let hydrate_line =
            "[QCortado] Remote SCF dependency hydration is enabled for this run.".to_string();
        let _ = app.emit(&format!("task-output:{}", task_id), &hydrate_line);
        pm.append_output(task_id, hydrate_line).await;
    }

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let mut commands = vec!["cd \"$SLURM_SUBMIT_DIR\"".to_string()];
    commands.extend(dependency_stage.remote_hydration_commands);
    commands.push(format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)));
    commands.push(build_hpc_logged_qe_step_command(
        &profile,
        resource_type,
        "pw.x (NSCF along k-path)",
        "pw.x",
        None,
        "bands.in",
        "bands.out",
    ));
    commands.push(build_hpc_logged_qe_step_command(
        &profile,
        resource_type,
        "bands.x post-processing",
        "bands.x",
        None,
        "bands_pp.in",
        "bands_pp.out",
    ));
    if projections_enabled {
        let projwfc_cmd = build_hpc_logged_qe_step_command(
            &profile,
            resource_type,
            "projwfc.x orbital projections",
            "projwfc.x",
            None,
            "projwfc.in",
            "projwfc.out",
        );
        commands.push(format!(
            "{} || echo \"WARNING: projwfc.x failed\"",
            projwfc_cmd
        ));
    }

    let mut bundle_files = vec![
        ("bands.in".to_string(), nscf_input),
        ("bands_pp.in".to_string(), bands_x_input),
    ];
    if projections_enabled {
        bundle_files.push(("projwfc.in".to_string(), projwfc_input));
    }

    let mut bundle_copies: Vec<(PathBuf, String)> = Vec::new();
    if let Some(local_scf_tmp_dir) = dependency_stage.local_bundle_copy {
        bundle_copies.push((local_scf_tmp_dir, ".".to_string()));
    }

    let submit_started = std::time::Instant::now();
    let submit_line = format!(
        "[QCortado] Submitting bands bundle to HPC at {}.",
        now_iso()
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &submit_line);
    pm.append_output(task_id, submit_line).await;

    let work_path = run_hpc_bundle_task(
        app.clone(),
        pm.clone(),
        task_id,
        "bands",
        "Bands",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files,
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;

    let submit_done_line = format!(
        "[QCortado] HPC run + artifact sync returned after {:.1}s. Parsing downloaded bands artifacts...",
        submit_started.elapsed().as_secs_f64()
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &submit_done_line);
    pm.append_output(task_id, submit_done_line).await;

    let bands_out_text =
        std::fs::read_to_string(work_path.join("bands.out")).unwrap_or_else(|_| {
            std::fs::read_to_string(work_path.join("slurm.out")).unwrap_or_default()
        });
    let fermi_energy = extract_fermi_energy_from_text(&bands_out_text).unwrap_or(0.0);

    let gnu_file_name = format!("{}.gnu", bands_x_config.filband);
    let gnu_file = work_path.join(&gnu_file_name);
    if !gnu_file.exists() {
        return Err(format!("{} not found after HPC bands run.", gnu_file_name));
    }

    let parse_started = std::time::Instant::now();
    let mut band_data = read_bands_gnu_file(&gnu_file, fermi_energy)
        .map_err(|e| format!("Failed to parse band data: {}", e))?;
    qe::bands::add_symmetry_markers(&mut band_data, &config.k_path);

    if projections_enabled {
        let projection_text = {
            let filproj_path = work_path.join(&projwfc_config.filproj);
            if filproj_path.exists() {
                std::fs::read_to_string(&filproj_path).unwrap_or_default()
            } else {
                std::fs::read_to_string(work_path.join("projwfc.out")).unwrap_or_default()
            }
        };
        if !projection_text.trim().is_empty() {
            if let Ok(projections) =
                parse_projwfc_projection_groups_aligned(&projection_text, &band_data.energies)
            {
                band_data.projections = Some(projections);
            }
        }
    }

    let done_line = format!(
        "[QCortado] Bands HPC pipeline complete in {:.1}s (post-download parse {:.1}s).",
        pipeline_start.elapsed().as_secs_f64(),
        parse_started.elapsed().as_secs_f64()
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &done_line);
    pm.append_output(task_id, done_line).await;

    Ok(band_data)
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
    prepare_working_directory(&work_path, false)?;
    let pipeline_start = std::time::Instant::now();

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
                emit_line!(format!(
                    "Verified .save in working dir: {}",
                    copied_save.display()
                ));
            } else {
                emit_line!("WARNING: .save not found in working dir after copy!".to_string());
            }
        } else {
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    check_cancel!();

    // Step 1: NSCF along k-path
    let mut bands_calc = config.base_calculation.clone();
    bands_calc.calculation = qe::CalculationType::Bands;
    if bands_calc.verbosity.is_none() {
        bands_calc.verbosity = Some("high".to_string());
    }

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

    let mut input = generate_pw_input(&bands_calc);
    if let Some(nbnd) = config.nbnd {
        input = insert_system_namelist_line(&input, &format!("nbnd = {},", nbnd))?;
        emit_line!(format!("Requested nbnd = {}", nbnd));
    }
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
            i + 1,
            point.label,
            point.coords[0],
            point.coords[1],
            point.coords[2],
            point.npoints
        ));
    }

    let projections_enabled = config
        .projections
        .as_ref()
        .map(|p| p.enabled)
        .unwrap_or(false);
    let total_steps = if projections_enabled { 3 } else { 2 };
    emit_line!(format!(
        "[QCortado] Bands local config: k_path_points={}, projections_enabled={}, total_steps={}, started_at={}.",
        config.k_path.len(),
        projections_enabled,
        total_steps,
        now_iso()
    ));
    emit_line!(format!(
        "Step 1/{}: Running NSCF calculation along k-path...",
        total_steps
    ));

    let exe_path = bin_dir.join("pw.x");
    if !exe_path.exists() {
        return Err("pw.x not found".to_string());
    }
    let nscf_started = std::time::Instant::now();

    let mut child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!("Using MPI with {} processes", mpi.nprocs));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&exe_path)
                .args(["-pd", ".true."])
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
        stdin
            .write_all(input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write input: {}", e))?;
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
    emit_line!(format!(
        "[QCortado] NSCF stage finished in {:.1}s.",
        nscf_started.elapsed().as_secs_f64()
    ));

    std::fs::write(work_path.join("bands.out"), &full_output)
        .map_err(|e| format!("Failed to write output file: {}", e))?;

    emit_line!("".to_string());
    emit_line!(format!(
        "Step 2/{}: Running bands.x post-processing...",
        total_steps
    ));

    // Step 2: bands.x
    let bands_x_path = bin_dir.join("bands.x");
    if !bands_x_path.exists() {
        return Err(
            "bands.x not found. Make sure your QE installation includes bands.x".to_string(),
        );
    }
    let bands_started = std::time::Instant::now();

    let bands_x_options = config.bands_x.as_ref();
    let bands_filband = bands_x_options
        .and_then(|opts| opts.filband.as_deref())
        .map(|raw| sanitize_output_filename(raw, "bands.dat"))
        .unwrap_or_else(|| "bands.dat".to_string());
    let bands_lsym = bands_x_options.and_then(|opts| opts.lsym).unwrap_or(true);
    let bands_no_overlap = bands_x_options
        .and_then(|opts| opts.no_overlap)
        .unwrap_or(true);

    let bands_x_config = BandsXConfig {
        prefix: bands_calc.prefix.clone(),
        outdir: bands_calc.outdir.clone(),
        filband: bands_filband.clone(),
        lsym: bands_lsym,
        no_overlap: bands_no_overlap,
    };
    let bands_x_input = generate_bands_x_input(&bands_x_config);
    std::fs::write(work_path.join("bands_pp.in"), &bands_x_input)
        .map_err(|e| format!("Failed to write bands.x input: {}", e))?;

    let mut bands_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!(
                "Running bands.x with MPI ({} processes)",
                mpi.nprocs
            ));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&bands_x_path)
                .args(["-pd", ".true."])
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun for bands.x: {}", e))?
        } else {
            tokio_command_with_prefix(&bands_x_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start bands.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&bands_x_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start bands.x: {}", e))?
    };

    if let Some(pid) = bands_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = bands_child.stdin.take() {
        stdin
            .write_all(bands_x_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write bands.x input: {}", e))?;
    }

    let bands_stdout = bands_child
        .stdout
        .take()
        .ok_or("Failed to capture bands.x stdout")?;
    let mut bands_reader = BufReader::new(bands_stdout).lines();

    while let Some(line) = bands_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        emit_line!(line);
    }

    let bands_status = bands_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !bands_status.success() {
        return Err(format!(
            "bands.x failed with exit code: {:?}",
            bands_status.code()
        ));
    }
    emit_line!(format!(
        "[QCortado] bands.x stage finished in {:.1}s.",
        bands_started.elapsed().as_secs_f64()
    ));

    emit_line!("".to_string());
    emit_line!("Parsing band structure data...".to_string());

    let gnu_file_name = format!("{}.gnu", bands_x_config.filband);
    let gnu_file = work_path.join(&gnu_file_name);
    if !gnu_file.exists() {
        return Err(format!(
            "{} not found. bands.x may have failed.",
            gnu_file_name
        ));
    }

    if let Ok(metadata) = std::fs::metadata(&gnu_file) {
        emit_line!(format!("{} size: {} bytes", gnu_file_name, metadata.len()));
    }

    let ef = fermi_energy.unwrap_or(0.0);
    emit_line!(format!("Using Fermi energy: {:.4} eV", ef));

    let parse_started = std::time::Instant::now();
    let mut band_data = read_bands_gnu_file(&gnu_file, ef)
        .map_err(|e| format!("Failed to parse band data: {}", e))?;

    emit_line!(format!(
        "Parsed: {} bands, {} k-points, energy range [{:.2}, {:.2}] eV",
        band_data.n_bands,
        band_data.n_kpoints,
        band_data.energy_range[0],
        band_data.energy_range[1]
    ));

    qe::bands::add_symmetry_markers(&mut band_data, &config.k_path);
    emit_line!(format!(
        "[QCortado] bands.dat.gnu parse finished in {:.1}s.",
        parse_started.elapsed().as_secs_f64()
    ));

    // Step 3: projwfc.x (optional)
    if projections_enabled {
        emit_line!("".to_string());
        emit_line!(format!(
            "Step 3/{}: Running projwfc.x orbital projections...",
            total_steps
        ));

        let projwfc_x_path = bin_dir.join("projwfc.x");
        if !projwfc_x_path.exists() {
            emit_line!(
                "WARNING: projwfc.x not found. Skipping fat-band projection analysis.".to_string()
            );
        } else {
            check_cancel!();
            let projection_options = config.projections.as_ref();
            let projection_file = projection_options
                .and_then(|opts| opts.filproj.as_deref())
                .map(|raw| sanitize_output_filename(raw, "bands.projwfc.dat"))
                .unwrap_or_else(|| "bands.projwfc.dat".to_string());
            let projwfc_config = ProjwfcConfig {
                prefix: bands_calc.prefix.clone(),
                outdir: bands_calc.outdir.clone(),
                filproj: projection_file,
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

            if let Some(pid) = projwfc_child.id() {
                pm.set_child_id(task_id, pid).await;
            }

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
                        std::fs::read_to_string(&filproj_path)
                            .unwrap_or_else(|_| projwfc_output.clone())
                    } else {
                        projwfc_output.clone()
                    }
                };

                if projection_text.trim().is_empty() {
                    emit_line!(
                        "WARNING: projwfc output was empty. Continuing without projections."
                            .to_string()
                    );
                } else {
                    match parse_projwfc_projection_groups_aligned(
                        &projection_text,
                        &band_data.energies,
                    ) {
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
    emit_line!(format!(
        "  {} bands, {} k-points",
        band_data.n_bands, band_data.n_kpoints
    ));
    if let Some(ref gap) = band_data.band_gap {
        let gap_type = if gap.is_direct { "direct" } else { "indirect" };
        emit_line!(format!("  Band gap: {:.3} eV ({})", gap.value, gap_type));
    }
    emit_line!(format!(
        "[QCortado] Bands local pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    ));

    Ok(band_data)
}

// ============================================================================
// Wannier90 Commands
// ============================================================================

/// Starts a scalar Wannier90 workflow as a background task.
#[tauri::command]
async fn start_wannier_calculation(
    app: AppHandle,
    config: WannierCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    validate_wannier_config(&config)?;
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("wannier".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();

        tokio::spawn(async move {
            let result = run_wannier_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
                cancel_flag,
                pm.clone(),
            )
            .await;

            match result {
                Ok(wannier_result) => {
                    let json =
                        serde_json::to_value(&wannier_result).unwrap_or(serde_json::Value::Null);
                    pm.complete(&tid, json).await;
                    let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
                }
                Err(e) => {
                    pm.fail(&tid, e.clone()).await;
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let qe_bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let wannier90_path = {
        let guard = state.wannier90_path.lock().unwrap();
        guard
            .as_ref()
            .ok_or("Wannier90 path not configured")?
            .clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let tid = task_id.clone();
    let app_handle = app.clone();
    tokio::spawn(async move {
        let result = run_wannier_background(
            app_handle.clone(),
            &tid,
            config,
            working_dir,
            mpi_config,
            qe_bin_dir,
            wannier90_path,
            execution_prefix,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(wannier_result) => {
                let json = serde_json::to_value(&wannier_result).unwrap_or(serde_json::Value::Null);
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
async fn run_wannier_background(
    app: AppHandle,
    task_id: &str,
    config: WannierCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    qe_bin_dir: PathBuf,
    wannier90_path: PathBuf,
    execution_prefix: Option<String>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<WannierResult, String> {
    let work_path = PathBuf::from(&working_dir);
    prepare_working_directory(&work_path, false)?;
    let pipeline_start = std::time::Instant::now();

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

    if let (Some(ref project_id), Some(ref scf_calc_id)) = (&config.project_id, &config.scf_calc_id)
    {
        let projects_dir = projects::get_projects_dir(&app)?;
        let scf_tmp_dir = projects_dir
            .join(project_id)
            .join("calculations")
            .join(scf_calc_id)
            .join("tmp");

        if scf_tmp_dir.exists() {
            emit_line!("Copying SCF data to working directory...".to_string());
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;
            emit_line!("SCF data copied successfully.".to_string());
        } else {
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    let pw2wannier_path = qe_bin_dir.join("pw2wannier90.x");
    if !pw2wannier_path.exists() {
        return Err("pw2wannier90.x not found in the configured QE bin directory.".to_string());
    }

    if !wannier90_path.exists() {
        return Err(format!(
            "Wannier90 executable not found: {}",
            wannier90_path.display()
        ));
    }

    let (nscf_calc, nscf_notes) = prepare_wannier_nscf_calculation(&config)?;

    if !matches!(nscf_calc.system.position_units, qe::PositionUnits::Crystal) {
        return Err(
            "Wannier v1 requires the base calculation to use crystal fractional atomic positions."
                .to_string(),
        );
    }

    let kpoints = match &nscf_calc.kpoints {
        qe::KPoints::Crystal { points } => points.clone(),
        _ => unreachable!(),
    };
    let win_content = generate_wannier90_win(&config, &kpoints)?;
    let nscf_input = generate_pw_input(&nscf_calc);
    let pw2_config = config
        .pw2wannier90
        .clone()
        .unwrap_or_else(Pw2Wannier90Config::default);
    let pw2wan_input = generate_pw2wannier90_input(&config, &pw2_config);

    std::fs::write(
        work_path.join(format!("{}.win", config.seedname)),
        &win_content,
    )
    .map_err(|e| format!("Failed to write Wannier .win file: {}", e))?;
    std::fs::write(work_path.join("nscf.in"), &nscf_input)
        .map_err(|e| format!("Failed to write nscf.in: {}", e))?;
    std::fs::write(work_path.join("pw2wan.in"), &pw2wan_input)
        .map_err(|e| format!("Failed to write pw2wan.in: {}", e))?;

    emit_line!(format!(
        "[QCortado] Wannier local config: seedname={}, num_wann={}, num_bands={}, mesh={}x{}x{}, started_at={}.",
        config.seedname,
        config.num_wann,
        config.num_bands,
        config.k_grid[0],
        config.k_grid[1],
        config.k_grid[2],
        now_iso()
    ));
    for note in nscf_notes {
        emit_line!(format!("[QCortado] {}", note));
    }

    emit_line!("Step 1/4: Running wannier90.x -pp preprocessing...".to_string());
    let pre_started = std::time::Instant::now();
    let pre_output = run_local_stage_capture_stdout(
        &app,
        &pm,
        task_id,
        &work_path,
        &wannier90_path,
        &["-pp", &config.seedname],
        None,
        None,
        None,
        false,
    )
    .await?;
    std::fs::write(work_path.join("wannier90_pre.out"), &pre_output)
        .map_err(|e| format!("Failed to write wannier90_pre.out: {}", e))?;
    emit_line!(format!(
        "[QCortado] wannier90.x -pp stage finished in {:.1}s.",
        pre_started.elapsed().as_secs_f64()
    ));
    check_cancel!();

    emit_line!("Step 2/4: Running pw.x NSCF on explicit full Monkhorst-Pack mesh...".to_string());
    let pw_started = std::time::Instant::now();
    let pw_output = run_local_stage_capture_stdout(
        &app,
        &pm,
        task_id,
        &work_path,
        &qe_bin_dir.join("pw.x"),
        &[],
        Some(&nscf_input),
        execution_prefix.as_deref(),
        mpi_config.as_ref(),
        true,
    )
    .await?;
    std::fs::write(work_path.join("nscf.out"), &pw_output)
        .map_err(|e| format!("Failed to write nscf.out: {}", e))?;
    let fermi_energy = extract_fermi_energy_from_text(&pw_output).unwrap_or(0.0);
    emit_line!(format!(
        "[QCortado] pw.x NSCF stage finished in {:.1}s.",
        pw_started.elapsed().as_secs_f64()
    ));
    check_cancel!();

    emit_line!("Step 3/4: Running pw2wannier90.x interface...".to_string());
    let pw2_started = std::time::Instant::now();
    let pw2_output = run_local_stage_capture_stdout(
        &app,
        &pm,
        task_id,
        &work_path,
        &pw2wannier_path,
        &[],
        Some(&pw2wan_input),
        execution_prefix.as_deref(),
        mpi_config.as_ref(),
        true,
    )
    .await?;
    std::fs::write(work_path.join("pw2wan.out"), &pw2_output)
        .map_err(|e| format!("Failed to write pw2wan.out: {}", e))?;
    emit_line!(format!(
        "[QCortado] pw2wannier90.x stage finished in {:.1}s.",
        pw2_started.elapsed().as_secs_f64()
    ));
    check_cancel!();

    emit_line!("Step 4/4: Running wannier90.x minimization and band interpolation...".to_string());
    let wan_started = std::time::Instant::now();
    let wan_output = run_local_stage_capture_stdout(
        &app,
        &pm,
        task_id,
        &work_path,
        &wannier90_path,
        &[&config.seedname],
        None,
        None,
        None,
        false,
    )
    .await?;
    std::fs::write(work_path.join("wannier90.out"), &wan_output)
        .map_err(|e| format!("Failed to write wannier90.out: {}", e))?;
    emit_line!(format!(
        "[QCortado] wannier90.x stage finished in {:.1}s.",
        wan_started.elapsed().as_secs_f64()
    ));
    check_cancel!();

    emit_line!("Parsing Wannier artifacts...".to_string());
    let parse_started = std::time::Instant::now();
    let result = read_wannier_result(&work_path, &config, fermi_energy)?;
    emit_line!(format!(
        "Parsed: {} Wannier functions, {} interpolation bands, {} k-points.",
        result.num_wann, result.band_data.n_bands, result.band_data.n_kpoints
    ));
    emit_line!(format!(
        "[QCortado] Wannier parse finished in {:.1}s.",
        parse_started.elapsed().as_secs_f64()
    ));
    emit_line!("=== Wannier Calculation Complete ===".to_string());
    emit_line!(format!(
        "[QCortado] Wannier local pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    ));

    Ok(result)
}

#[allow(clippy::too_many_arguments)]
async fn run_wannier_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: WannierCalculationConfig,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<WannierResult, String> {
    let pipeline_start = std::time::Instant::now();
    let remote_wannier90 = profile
        .remote_wannier90_path
        .as_deref()
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "wannier90.x".to_string());

    let (mut nscf_calc, nscf_notes) = prepare_wannier_nscf_calculation(&config)?;
    nscf_calc.pseudo_dir = profile.remote_pseudo_dir.clone();

    let kpoints = match &nscf_calc.kpoints {
        qe::KPoints::Crystal { points } => points.clone(),
        _ => unreachable!(),
    };
    let win_content = generate_wannier90_win(&config, &kpoints)?;
    let nscf_input = generate_pw_input(&nscf_calc);
    let pw2_config = config
        .pw2wannier90
        .clone()
        .unwrap_or_else(Pw2Wannier90Config::default);
    let pw2wan_input = generate_pw2wannier90_input(&config, &pw2_config);

    let dependency_stage = resolve_hpc_scf_dependency_stage(
        &app,
        config.project_id.as_deref(),
        config.scf_calc_id.as_deref(),
    )?;

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let pre_cmd = format!(
        "{} -pp {} > wannier90_pre.out 2>&1",
        shell_single_quote_local(&remote_wannier90),
        shell_single_quote_local(&config.seedname)
    );
    let final_cmd = format!(
        "{} {} > wannier90.out 2>&1",
        shell_single_quote_local(&remote_wannier90),
        shell_single_quote_local(&config.seedname)
    );

    let mut commands = vec!["cd \"$SLURM_SUBMIT_DIR\"".to_string()];
    commands.extend(dependency_stage.remote_hydration_commands);
    commands.push(format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)));
    commands.push(build_hpc_logged_shell_step_command(
        "wannier90.x preprocessing",
        &pre_cmd,
    ));
    commands.push(build_hpc_logged_qe_step_command(
        &profile,
        resource_type,
        "pw.x NSCF on full k-mesh",
        "pw.x",
        None,
        "nscf.in",
        "nscf.out",
    ));
    commands.push(build_hpc_logged_qe_step_command(
        &profile,
        resource_type,
        "pw2wannier90.x interface",
        "pw2wannier90.x",
        None,
        "pw2wan.in",
        "pw2wan.out",
    ));
    commands.push(build_hpc_logged_shell_step_command(
        "wannier90.x minimization",
        &final_cmd,
    ));

    let mut bundle_files = vec![
        (format!("{}.win", config.seedname), win_content),
        ("nscf.in".to_string(), nscf_input),
        ("pw2wan.in".to_string(), pw2wan_input),
    ];
    let mut bundle_copies: Vec<(PathBuf, String)> = Vec::new();
    if let Some(local_scf_tmp_dir) = dependency_stage.local_bundle_copy {
        bundle_copies.push((local_scf_tmp_dir, ".".to_string()));
    }

    let config_line = format!(
        "[QCortado] Wannier HPC config: seedname={}, num_wann={}, num_bands={}, mesh={}x{}x{}.",
        config.seedname,
        config.num_wann,
        config.num_bands,
        config.k_grid[0],
        config.k_grid[1],
        config.k_grid[2]
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &config_line);
    pm.append_output(task_id, config_line).await;
    for note in nscf_notes {
        let note_line = format!("[QCortado] {}", note);
        let _ = app.emit(&format!("task-output:{}", task_id), &note_line);
        pm.append_output(task_id, note_line).await;
    }

    let work_path = run_hpc_bundle_task(
        app.clone(),
        pm.clone(),
        task_id,
        "wannier",
        "Wannier",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files.drain(..).collect(),
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;

    let nscf_out_text = std::fs::read_to_string(work_path.join("nscf.out")).unwrap_or_else(|_| {
        std::fs::read_to_string(work_path.join("slurm.out")).unwrap_or_default()
    });
    let fermi_energy = extract_fermi_energy_from_text(&nscf_out_text).unwrap_or(0.0);
    let result = read_wannier_result(&work_path, &config, fermi_energy)?;

    let done_line = format!(
        "[QCortado] Wannier HPC pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &done_line);
    pm.append_output(task_id, done_line).await;

    Ok(result)
}

#[tauri::command]
async fn start_transport_calculation(
    app: AppHandle,
    config: TransportCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    validate_transport_config(&config)?;
    let source =
        validate_transport_source(&app, &config.project_id, &config.source_wannier_calc_id)?;
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("transport".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();
        tokio::spawn(async move {
            let result = run_transport_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                source,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
                cancel_flag,
                pm.clone(),
            )
            .await;

            match result {
                Ok(transport_result) => {
                    let json =
                        serde_json::to_value(&transport_result).unwrap_or(serde_json::Value::Null);
                    pm.complete(&tid, json).await;
                    let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
                }
                Err(e) => {
                    pm.fail(&tid, e.clone()).await;
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let postw90_path = resolve_local_postw90_path(&state)?;
    if !postw90_path.exists() {
        return Err(format!(
            "postw90.x executable not found: {}",
            postw90_path.display()
        ));
    }
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let tid = task_id.clone();
    let app_handle = app.clone();
    tokio::spawn(async move {
        let result = run_transport_background(
            app_handle.clone(),
            &tid,
            config,
            source,
            working_dir,
            mpi_config,
            postw90_path,
            execution_prefix,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(transport_result) => {
                let json =
                    serde_json::to_value(&transport_result).unwrap_or(serde_json::Value::Null);
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
async fn run_transport_background(
    app: AppHandle,
    task_id: &str,
    config: TransportCalculationConfig,
    source: ValidatedTransportSource,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    postw90_path: PathBuf,
    execution_prefix: Option<String>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<TransportResult, String> {
    let work_path = PathBuf::from(&working_dir);
    prepare_working_directory(&work_path, false)?;
    let pipeline_start = std::time::Instant::now();

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

    emit_line!(format!(
        "[QCortado] Transport local config: source_wannier={}, seedname={}, mu_offset=[{}, {}, step {}] eV, temp=[{}, {}, step {}] K, tau={} fs.",
        source.calculation.id,
        source.seedname,
        config.mu_offset_min,
        config.mu_offset_max,
        config.mu_offset_step,
        config.temp_min,
        config.temp_max,
        config.temp_step,
        config.relaxation_time_fs
    ));

    emit_line!("Preparing transport input...".to_string());
    stage_transport_source_files(&source, &work_path)?;
    let transport_win = build_transport_win(
        &source.source_win_content,
        &config,
        source.reference_fermi_energy_ev,
    )?;
    std::fs::write(
        work_path.join(format!("{}.win", source.seedname)),
        &transport_win,
    )
    .map_err(|e| format!("Failed to write transport .win file: {}", e))?;
    check_cancel!();

    emit_line!("Step 1/2: Running postw90.x / BoltzWann...".to_string());
    let run_started = std::time::Instant::now();
    let stage_output = run_local_stage_capture_streams(
        &app,
        &pm,
        task_id,
        &work_path,
        &postw90_path,
        &[&source.seedname],
        None,
        execution_prefix.as_deref(),
        mpi_config.as_ref(),
        true,
    )
    .await?;
    std::fs::write(
        work_path.join(format!("{}.wpout", source.seedname)),
        &stage_output.stdout,
    )
    .map_err(|e| format!("Failed to write {}.wpout: {}", source.seedname, e))?;
    std::fs::write(
        work_path.join(format!("{}.werr", source.seedname)),
        &stage_output.stderr,
    )
    .map_err(|e| format!("Failed to write {}.werr: {}", source.seedname, e))?;
    emit_line!(format!(
        "[QCortado] postw90.x stage finished in {:.1}s.",
        run_started.elapsed().as_secs_f64()
    ));
    check_cancel!();

    emit_line!("Step 2/2: Parsing BoltzWann output...".to_string());
    let parse_started = std::time::Instant::now();
    let result = parse_transport_result(
        &work_path,
        &config,
        &source.seedname,
        source.reference_fermi_energy_ev,
    )?;
    emit_line!(format!(
        "Parsed transport grid: {} temperatures, {} chemical potentials, {} conductivity components.",
        result.temperature_values_k.len(),
        result.mu_values_ev.len(),
        result.conductivity.component_labels.len()
    ));
    emit_line!(format!(
        "[QCortado] Transport parse finished in {:.1}s.",
        parse_started.elapsed().as_secs_f64()
    ));
    emit_line!("=== Transport Calculation Complete ===".to_string());
    emit_line!(format!(
        "[QCortado] Transport local pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    ));

    Ok(result)
}

#[allow(clippy::too_many_arguments)]
async fn run_transport_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: TransportCalculationConfig,
    source: ValidatedTransportSource,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<TransportResult, String> {
    let pipeline_start = std::time::Instant::now();
    let remote_postw90 = derive_remote_postw90_path(profile.remote_wannier90_path.as_deref());

    let stage_dir = PathBuf::from(&working_dir)
        .join("transport_source_stage")
        .join(task_id);
    prepare_working_directory(&stage_dir, false)?;
    stage_transport_source_files(&source, &stage_dir)?;
    let transport_win = build_transport_win(
        &source.source_win_content,
        &config,
        source.reference_fermi_energy_ev,
    )?;

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let launcher = build_hpc_launcher_command(&profile, resource_type);
    let transport_cmd = format!(
        "{} {} {} > {}.wpout 2> {}.werr",
        launcher,
        shell_single_quote_local(&remote_postw90),
        shell_single_quote_local(&source.seedname),
        shell_single_quote_local(&source.seedname),
        shell_single_quote_local(&source.seedname)
    );
    let commands = vec![
        "cd \"$SLURM_SUBMIT_DIR\"".to_string(),
        build_hpc_logged_shell_step_command("postw90.x / BoltzWann", &transport_cmd),
    ];
    let bundle_files = vec![(format!("{}.win", source.seedname), transport_win)];
    let bundle_copies = vec![(stage_dir.clone(), ".".to_string())];

    let config_line = format!(
        "[QCortado] Transport HPC config: source_wannier={}, seedname={}, mu_offset=[{}, {}, step {}] eV, temp=[{}, {}, step {}] K, tau={} fs.",
        source.calculation.id,
        source.seedname,
        config.mu_offset_min,
        config.mu_offset_max,
        config.mu_offset_step,
        config.temp_min,
        config.temp_max,
        config.temp_step,
        config.relaxation_time_fs
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &config_line);
    pm.append_output(task_id, config_line).await;

    let work_path = run_hpc_bundle_task(
        app.clone(),
        pm.clone(),
        task_id,
        "transport",
        "Transport",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files,
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;
    let _ = std::fs::remove_dir_all(&stage_dir);

    let result = parse_transport_result(
        &work_path,
        &config,
        &source.seedname,
        source.reference_fermi_energy_ev,
    )?;

    let done_line = format!(
        "[QCortado] Transport HPC pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &done_line);
    pm.append_output(task_id, done_line).await;

    Ok(result)
}

fn build_epw_sources_payload(
    config: &EpwCalculationConfig,
    sources: &ValidatedEpwSources,
) -> EpwSourcesV1 {
    let fallback_scf_id =
        extract_string_param(&sources.phonon_calculation.parameters, "source_scf_id");
    let resolved_scf_id = sources.scf_calculation.as_ref().map(|calc| calc.id.clone());
    let scf_ref = config
        .source_scf_calc_id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
        .or(resolved_scf_id)
        .or(fallback_scf_id)
        .map(|calc_id| EpwSourceRef {
            calc_id,
            calc_type: "scf".to_string(),
        });

    EpwSourcesV1 {
        phonon: EpwSourceRef {
            calc_id: sources.phonon_calculation.id.clone(),
            calc_type: "phonon".to_string(),
        },
        wannier: sources
            .wannier_calculation
            .as_ref()
            .map(|calculation| EpwSourceRef {
                calc_id: calculation.id.clone(),
                calc_type: "wannier".to_string(),
            }),
        scf: scf_ref,
        manifests: sources.manifests.clone(),
    }
}

fn build_epw_taxonomy_error(
    code: &str,
    message: impl Into<String>,
    hint: Option<String>,
) -> String {
    let message = message.into();
    match hint {
        Some(hint) if !hint.trim().is_empty() => {
            format!("[{}] {} Hint: {}", code, message, hint.trim())
        }
        _ => format!("[{}] {}", code, message),
    }
}

fn build_epw_validation_error(validation: &EpwPrerequisiteValidation) -> String {
    let mut lines: Vec<String> = validation
        .errors
        .iter()
        .map(|entry| format!("- {}", entry))
        .collect();
    if !validation.remediation_hints.is_empty() {
        lines.push("Remediation hints:".to_string());
        for hint in &validation.remediation_hints {
            lines.push(format!("  - {}", hint));
        }
    }
    build_epw_taxonomy_error(
        "missing-prereq",
        format!("EPW prerequisite validation failed.\n{}", lines.join("\n")),
        None,
    )
}

#[derive(Debug, Clone, Copy)]
struct EpwParallelPlan {
    npool: u32,
    nimage: u32,
    npool_explicit: bool,
}

#[derive(Debug, Clone)]
struct EpwWannierizationPlan {
    wannierize: bool,
    filukk: String,
    explicit_wannierize_override: bool,
}

fn strip_epw_keyword_value(raw_value: &str) -> String {
    raw_value
        .split_once('!')
        .map(|(value, _)| value)
        .unwrap_or(raw_value)
        .trim()
        .trim_end_matches(',')
        .trim()
        .to_string()
}

fn parse_epw_logical_keyword(raw_value: &str) -> Option<bool> {
    let normalized = strip_epw_keyword_value(raw_value)
        .trim_matches('.')
        .to_ascii_lowercase();
    match normalized.as_str() {
        "true" | "t" => Some(true),
        "false" | "f" => Some(false),
        _ => None,
    }
}

fn parse_epw_string_keyword(raw_value: &str) -> Option<String> {
    let stripped = strip_epw_keyword_value(raw_value);
    if stripped.is_empty() {
        return None;
    }
    if stripped.len() >= 2 && stripped.starts_with('\'') && stripped.ends_with('\'') {
        return Some(stripped[1..stripped.len() - 1].replace("''", "'"));
    }
    if stripped.len() >= 2 && stripped.starts_with('"') && stripped.ends_with('"') {
        return Some(stripped[1..stripped.len() - 1].to_string());
    }
    Some(stripped)
}

fn has_explicit_epw_override(config: &EpwCalculationConfig, key: &str) -> bool {
    config
        .advanced_overrides
        .iter()
        .any(|(raw_key, raw_value)| {
            raw_key.trim().eq_ignore_ascii_case(key) && !raw_value.trim().is_empty()
        })
}

fn resolve_epw_wannierization_plan(
    config: &EpwCalculationConfig,
) -> Result<EpwWannierizationPlan, String> {
    let merged = build_epw_keyword_map(config)?;
    let wannierize = merged
        .get("wannierize")
        .and_then(|raw| parse_epw_logical_keyword(raw))
        .unwrap_or(config.input.wannierize);
    let filukk = merged
        .get("filukk")
        .and_then(|raw| parse_epw_string_keyword(raw))
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| format!("{}.ukk", config.input.prefix.trim()));

    Ok(EpwWannierizationPlan {
        wannierize,
        filukk,
        explicit_wannierize_override: has_explicit_epw_override(config, "wannierize"),
    })
}

fn epw_manifest_contains_file(
    manifests: &[EpwArtifactManifestEntry],
    target_file_name: &str,
) -> bool {
    let normalized_target = target_file_name
        .trim()
        .trim_matches('/')
        .to_ascii_lowercase();
    if normalized_target.is_empty() {
        return false;
    }
    manifests.iter().any(|entry| {
        let rel = entry.rel_path.trim().trim_matches('/').to_ascii_lowercase();
        rel == normalized_target || rel.ends_with(&format!("/{}", normalized_target))
    })
}

fn parse_epw_positive_u32_keyword(raw_value: &str, key: &str) -> Result<u32, String> {
    let trimmed = strip_epw_keyword_value(raw_value);
    let parsed = trimmed.parse::<u32>().map_err(|_| {
        format!(
            "EPW keyword `{}` must be a positive integer, found '{}'.",
            key,
            raw_value.trim()
        )
    })?;
    if parsed == 0 {
        return Err(format!("EPW keyword `{}` must be greater than zero.", key));
    }
    Ok(parsed)
}

fn build_epw_parallel_plan(
    config: &EpwCalculationConfig,
    launch_ranks: u32,
) -> Result<EpwParallelPlan, String> {
    let merged = build_epw_keyword_map(config)?;
    let configured_npool = merged
        .get("npool")
        .map(|raw| parse_epw_positive_u32_keyword(raw, "npool"))
        .transpose()?;
    let configured_nimage = merged
        .get("nimage")
        .or_else(|| merged.get("images"))
        .map(|raw| parse_epw_positive_u32_keyword(raw, "nimage"))
        .transpose()?;

    let npool = configured_npool.unwrap_or(launch_ranks.max(1));
    let nimage = configured_nimage.unwrap_or(1);
    let required_ranks = npool.checked_mul(nimage).ok_or_else(|| {
        format!(
            "EPW parallel decomposition overflow: npool={} and nimage={} are too large.",
            npool, nimage
        )
    })?;
    if required_ranks != launch_ranks {
        return Err(format!(
            "EPW parallel decomposition mismatch: launch ranks={} but npool={} and nimage={} require {} ranks.",
            launch_ranks, npool, nimage, required_ranks
        ));
    }

    Ok(EpwParallelPlan {
        npool,
        nimage,
        npool_explicit: configured_npool.is_some(),
    })
}

fn resolve_local_epw_launch_ranks(mpi_config: Option<&MpiConfig>) -> u32 {
    mpi_config
        .filter(|mpi| mpi.enabled && mpi.nprocs > 1)
        .map(|mpi| mpi.nprocs)
        .unwrap_or(1)
}

fn resolve_hpc_epw_launch_ranks(resources: Option<&hpc::profile::SlurmResourceRequest>) -> u32 {
    resources.and_then(|value| value.ntasks).unwrap_or(1).max(1)
}

/// Starts an EPW (`epw.x`) workflow as a background task.
#[tauri::command]
async fn start_epw_calculation(
    app: AppHandle,
    config: EpwCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    validate_epw_config(&config)?;
    let (validation, validated_sources) = validate_epw_prerequisites_detailed(&app, &config);
    if !validation.ok {
        return Err(build_epw_validation_error(&validation));
    }
    let sources = validated_sources.ok_or_else(|| {
        build_epw_taxonomy_error(
            "missing-prereq",
            "EPW prerequisites were not resolved.",
            Some("Re-run prerequisite validation and fix reported issues.".to_string()),
        )
    })?;

    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());
    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("epw".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();
        tokio::spawn(async move {
            let result = run_epw_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                sources,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
                cancel_flag,
                pm.clone(),
            )
            .await;

            match result {
                Ok(epw_result) => {
                    let json = serde_json::to_value(&epw_result).unwrap_or(serde_json::Value::Null);
                    pm.complete(&tid, json).await;
                    let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
                }
                Err(err) => {
                    pm.fail(&tid, err.clone()).await;
                    let _ = app_handle
                        .emit(&format!("task-status:{}", tid), &format!("failed:{}", err));
                }
            }
        });
        return Ok(task_id);
    }

    let qe_bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();
    let tid = task_id.clone();
    let app_handle = app.clone();

    tokio::spawn(async move {
        let result = run_epw_background(
            app_handle.clone(),
            &tid,
            config,
            sources,
            working_dir,
            mpi_config,
            qe_bin_dir,
            execution_prefix,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(epw_result) => {
                let json = serde_json::to_value(&epw_result).unwrap_or(serde_json::Value::Null);
                pm.complete(&tid, json).await;
                let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
            }
            Err(err) => {
                pm.fail(&tid, err.clone()).await;
                let _ =
                    app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", err));
            }
        }
    });

    Ok(task_id)
}

#[tauri::command]
async fn start_epw_calculation_hpc(
    app: AppHandle,
    config: EpwCalculationConfig,
    working_dir: String,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    profile_id: Option<String>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    start_epw_calculation(
        app,
        config,
        working_dir,
        None,
        Some(hpc::profile::ExecutionTarget {
            mode: hpc::profile::ExecutionMode::Hpc,
            hpc: Some(hpc::profile::HpcExecutionTarget {
                profile_id,
                resources,
                interactive_debug: false,
                recovery_save: None,
            }),
        }),
        label,
        state,
    )
    .await
}

#[allow(clippy::too_many_arguments)]
async fn run_epw_background(
    app: AppHandle,
    task_id: &str,
    config: EpwCalculationConfig,
    sources: ValidatedEpwSources,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    qe_bin_dir: PathBuf,
    execution_prefix: Option<String>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<EpwCalculationV1, String> {
    let mut config = config;
    let work_path = PathBuf::from(&working_dir);
    prepare_working_directory(&work_path, false)?;
    let pipeline_start = std::time::Instant::now();

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

    emit_line!(format!(
        "[QCortado] EPW local config: phonon_source={}, wannier_source={}, prefix={}, coarse_k={}x{}x{}, fine_k={}x{}x{}, coarse_q={}x{}x{}, fine_q={}x{}x{}, started_at={}.",
        sources.phonon_calculation.id,
        sources
            .wannier_calculation
            .as_ref()
            .map(|calc| calc.id.as_str())
            .unwrap_or("none"),
        config.input.prefix,
        epw_coarse_k_mesh(&config.input)[0],
        epw_coarse_k_mesh(&config.input)[1],
        epw_coarse_k_mesh(&config.input)[2],
        epw_fine_k_mesh(&config.input)[0],
        epw_fine_k_mesh(&config.input)[1],
        epw_fine_k_mesh(&config.input)[2],
        epw_coarse_q_mesh(&config.input)[0],
        epw_coarse_q_mesh(&config.input)[1],
        epw_coarse_q_mesh(&config.input)[2],
        epw_fine_q_mesh(&config.input)[0],
        epw_fine_q_mesh(&config.input)[1],
        epw_fine_q_mesh(&config.input)[2],
        now_iso()
    ));
    for warning in &sources.warnings {
        emit_line!(format!("[QCortado] WARNING: {}", warning));
    }

    let wannierization_plan = resolve_epw_wannierization_plan(&config).map_err(|message| {
        build_epw_taxonomy_error(
            "incompatible-prereq",
            message,
            Some("Fix EPW keyword overrides and retry.".to_string()),
        )
    })?;
    if !wannierization_plan.wannierize
        && !epw_manifest_contains_file(&sources.manifests, &wannierization_plan.filukk)
    {
        if wannierization_plan.explicit_wannierize_override {
            return Err(build_epw_taxonomy_error(
                "missing-prereq",
                format!(
                    "EPW is configured with `wannierize = .false.` but '{}' is missing.",
                    wannierization_plan.filukk
                ),
                Some(
                    "Set `wannierize = .true.` (or remove the override), or provide the required `.ukk` file."
                        .to_string(),
                ),
            ));
        }
        config.input.wannierize = true;
        emit_line!(format!(
            "[QCortado] No '{}' prerequisite was found. Auto-enabling `wannierize = .true.` for this run.",
            wannierization_plan.filukk
        ));
    }

    let launch_ranks = resolve_local_epw_launch_ranks(mpi_config.as_ref());
    let parallel_plan = build_epw_parallel_plan(&config, launch_ranks).map_err(|message| {
        build_epw_taxonomy_error(
            "incompatible-prereq",
            message,
            Some(
                "Set EPW pools/advanced overrides so npool*nimage equals local MPI ranks (or adjust local MPI process count)."
                    .to_string(),
            ),
        )
    })?;
    if !parallel_plan.npool_explicit && launch_ranks > 1 {
        emit_line!(format!(
            "[QCortado] No EPW npool was configured; defaulting to npool={} to match local MPI ranks.",
            parallel_plan.npool
        ));
    }
    emit_line!(format!(
        "[QCortado] EPW parallel decomposition: mpi_ranks={}, npool={}, nimage={}.",
        launch_ranks, parallel_plan.npool, parallel_plan.nimage
    ));

    emit_line!("Staging EPW prerequisites...".to_string());
    stage_epw_source_files(&sources, &work_path, &config)?;
    check_cancel!();

    let total_steps = if sources.rebuild_wannier_nscf_save {
        3
    } else {
        2
    };
    if sources.rebuild_wannier_nscf_save {
        let nscf_input_path = work_path.join("nscf.in");
        if !nscf_input_path.exists() {
            return Err(build_epw_taxonomy_error(
                "missing-prereq",
                "EPW needs to rebuild the coarse-grid Wannier NSCF save, but staged `nscf.in` was not found.",
                Some(
                    "Re-run the Wannier workflow so QCortado can preserve the NSCF input and save payload."
                        .to_string(),
                ),
            ));
        }

        let pw_path = qe_bin_dir.join("pw.x");
        if !pw_path.exists() {
            return Err(build_epw_taxonomy_error(
                "executable-missing",
                format!("pw.x not found in QE bin directory: {}", pw_path.display()),
                Some("Configure QE so pw.x is available, then retry.".to_string()),
            ));
        }

        emit_line!(format!(
            "Step 1/{}: Rebuilding coarse-grid NSCF save from the selected Wannier source...",
            total_steps
        ));
        let rebuild_started = std::time::Instant::now();
        let rebuild_output = run_local_stage_capture_streams(
            &app,
            &pm,
            task_id,
            &work_path,
            &pw_path,
            &["-in", "nscf.in"],
            None,
            execution_prefix.as_deref(),
            mpi_config.as_ref(),
            true,
        )
        .await
        .map_err(|err| build_epw_taxonomy_error("run-failed", err, None))?;
        std::fs::write(work_path.join("nscf_rebuild.out"), &rebuild_output.stdout)
            .map_err(|e| format!("Failed to write nscf_rebuild.out: {}", e))?;
        std::fs::write(work_path.join("nscf_rebuild.err"), &rebuild_output.stderr)
            .map_err(|e| format!("Failed to write nscf_rebuild.err: {}", e))?;

        let rebuilt_schema_path = work_path
            .join(resolve_epw_stage_subdir(&config.input.outdir, "tmp"))
            .join(format!("{}.save", config.input.prefix.trim()))
            .join("data-file-schema.xml");
        if !rebuilt_schema_path.exists() {
            return Err(build_epw_taxonomy_error(
                "missing-prereq",
                format!(
                    "Coarse-grid NSCF rebuild finished, but {} was not created.",
                    rebuilt_schema_path.display()
                ),
                Some(
                    "Inspect nscf_rebuild.out and ensure the selected Wannier source coarse mesh matches the EPW k-mesh."
                        .to_string(),
                ),
            ));
        }
        emit_line!(format!(
            "[QCortado] coarse-grid NSCF rebuild finished in {:.1}s.",
            rebuild_started.elapsed().as_secs_f64()
        ));
        check_cancel!();
    }

    let input_text = build_epw_input(&config)?;
    std::fs::write(work_path.join("epw.in"), &input_text)
        .map_err(|e| format!("Failed to write epw.in: {}", e))?;

    let epw_path = resolve_local_epw_path(&qe_bin_dir);
    let epw_path = match epw_path {
        Some(path) => path,
        None => {
            return Err(build_epw_taxonomy_error(
                "executable-missing",
                format!(
                    "epw.x not found under QE bin directory {} (checked bin and EPW/bin fallbacks).",
                    qe_bin_dir.display()
                ),
                Some(
                    "QCortado expects epw.x in <qe>/bin or <qe>/EPW/bin. Verify your QE install layout and retry."
                        .to_string(),
                ),
            ));
        }
    };

    emit_line!(format!(
        "Step {}/{}: Running epw.x...",
        if sources.rebuild_wannier_nscf_save {
            2
        } else {
            1
        },
        total_steps
    ));
    let run_started = std::time::Instant::now();
    let stage_args = vec![
        "-npool".to_string(),
        parallel_plan.npool.to_string(),
        "-in".to_string(),
        "epw.in".to_string(),
    ];
    let stage_arg_refs: Vec<&str> = stage_args.iter().map(String::as_str).collect();
    let stage_output = run_local_stage_capture_streams(
        &app,
        &pm,
        task_id,
        &work_path,
        &epw_path,
        &stage_arg_refs,
        None,
        execution_prefix.as_deref(),
        mpi_config.as_ref(),
        true,
    )
    .await
    .map_err(|err| build_epw_taxonomy_error("run-failed", err, None))?;
    std::fs::write(work_path.join("epw.out"), &stage_output.stdout)
        .map_err(|e| format!("Failed to write epw.out: {}", e))?;
    std::fs::write(work_path.join("epw.err"), &stage_output.stderr)
        .map_err(|e| format!("Failed to write epw.err: {}", e))?;
    emit_line!(format!(
        "[QCortado] epw.x stage finished in {:.1}s.",
        run_started.elapsed().as_secs_f64()
    ));
    check_cancel!();

    emit_line!(format!(
        "Step {}/{}: Parsing EPW outputs...",
        if sources.rebuild_wannier_nscf_save {
            3
        } else {
            2
        },
        total_steps
    ));
    let parse_started = std::time::Instant::now();
    let artifacts = collect_epw_artifacts(&work_path);
    let combined_output = format!("{}\n{}", stage_output.stdout, stage_output.stderr);
    let parsed_result = parse_epw_result_v2(&combined_output, &work_path, artifacts.clone(), true);
    let mut errors: Vec<EpwErrorRecord> = Vec::new();
    if parsed_result.summary.parse_partial {
        errors.push(EpwErrorRecord {
            code: "parse-partial".to_string(),
            message: "EPW output parse is partial.".to_string(),
            hint: Some("Inspect epw.out for full context.".to_string()),
        });
    }
    emit_line!(format!(
        "[QCortado] EPW parse finished in {:.1}s.",
        parse_started.elapsed().as_secs_f64()
    ));
    emit_line!(format!(
        "[QCortado] EPW local pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    ));

    Ok(EpwCalculationV1 {
        schema_version: EPW_SCHEMA_VERSION,
        sources: build_epw_sources_payload(&config, &sources),
        input: config.input,
        goals: config.goals,
        runtime: config.runtime,
        extensions: config.extensions,
        artifacts,
        result_summary: parsed_result.summary,
        transport: parsed_result.transport,
        superconductivity: parsed_result.superconductivity,
        spectral: parsed_result.spectral,
        parsed_tables: parsed_result.parsed_tables,
        errors,
    })
}

#[allow(clippy::too_many_arguments)]
async fn run_epw_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: EpwCalculationConfig,
    sources: ValidatedEpwSources,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<EpwCalculationV1, String> {
    let mut config = config;
    let pipeline_start = std::time::Instant::now();
    let stage_dir = PathBuf::from(&working_dir)
        .join("epw_source_stage")
        .join(task_id);
    prepare_working_directory(&stage_dir, false)?;
    stage_epw_source_files(&sources, &stage_dir, &config)?;
    let wannierization_plan = resolve_epw_wannierization_plan(&config).map_err(|message| {
        build_epw_taxonomy_error(
            "incompatible-prereq",
            message,
            Some("Fix EPW keyword overrides and retry.".to_string()),
        )
    })?;
    if !wannierization_plan.wannierize
        && !epw_manifest_contains_file(&sources.manifests, &wannierization_plan.filukk)
    {
        if wannierization_plan.explicit_wannierize_override {
            return Err(build_epw_taxonomy_error(
                "missing-prereq",
                format!(
                    "EPW is configured with `wannierize = .false.` but '{}' is missing.",
                    wannierization_plan.filukk
                ),
                Some(
                    "Set `wannierize = .true.` (or remove the override), or provide the required `.ukk` file."
                        .to_string(),
                ),
            ));
        }
        config.input.wannierize = true;
        let line = format!(
            "[QCortado] No '{}' prerequisite was found. Auto-enabling `wannierize = .true.` for this run.",
            wannierization_plan.filukk
        );
        let _ = app.emit(&format!("task-output:{}", task_id), &line);
        pm.append_output(task_id, line).await;
    }

    let launch_ranks = resolve_hpc_epw_launch_ranks(resources.as_ref());
    let parallel_plan = build_epw_parallel_plan(&config, launch_ranks).map_err(|message| {
        build_epw_taxonomy_error(
            "incompatible-prereq",
            message,
            Some(
                "Set EPW pools/advanced overrides so npool*nimage equals requested HPC MPI ranks (`ntasks`)."
                    .to_string(),
            ),
        )
    })?;
    let input_text = build_epw_input(&config)?;
    let total_steps = if sources.rebuild_wannier_nscf_save {
        3
    } else {
        2
    };

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let launcher = build_hpc_launcher_command(&profile, resource_type);
    let remote_epw_override = profile
        .remote_epw_path
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty());
    let epw_cmd = format!(
        "{} \"$EPW_EXE\" -npool {} -in epw.in > epw.out 2> epw.err",
        launcher, parallel_plan.npool
    );
    let epw_resolver_cmd = if let Some(remote_epw) = remote_epw_override {
        if remote_epw.contains('/') || remote_epw.starts_with('~') {
            format!(
                "tool={}; \
if [ \"$tool\" = \"~\" ]; then tool=\"$HOME\"; elif [ \"${{tool#~/}}\" != \"$tool\" ]; then tool=\"$HOME/${{tool#~/}}\"; fi; \
if [ -x \"$tool\" ]; then EPW_EXE=\"$tool\"; \
else echo \"[QCortado] ERROR: epw.x not found/executable at configured path $tool\" >&2; exit 1; fi",
                shell_single_quote_local(remote_epw)
            )
        } else {
            format!(
                "tool={}; \
if command -v \"$tool\" >/dev/null 2>&1; then EPW_EXE=\"$(command -v \"$tool\")\"; \
else echo \"[QCortado] ERROR: epw.x command '$tool' not found in PATH\" >&2; exit 1; fi",
                shell_single_quote_local(remote_epw)
            )
        }
    } else {
        "if [ -x \"$QE_BIN/epw.x\" ]; then EPW_EXE=\"$QE_BIN/epw.x\"; \
elif [ -x \"$QE_BIN/EPW/bin/epw.x\" ]; then EPW_EXE=\"$QE_BIN/EPW/bin/epw.x\"; \
elif [ -x \"$(dirname \"$QE_BIN\")/EPW/bin/epw.x\" ]; then EPW_EXE=\"$(dirname \"$QE_BIN\")/EPW/bin/epw.x\"; \
else echo \"[QCortado] ERROR: epw.x not found in $QE_BIN or EPW/bin fallback paths\" >&2; exit 1; fi"
            .to_string()
    };
    let mut commands = vec![
        "cd \"$SLURM_SUBMIT_DIR\"".to_string(),
        format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)),
        epw_resolver_cmd,
    ];
    if sources.rebuild_wannier_nscf_save {
        if !stage_dir.join("nscf.in").exists() {
            return Err(build_epw_taxonomy_error(
                "missing-prereq",
                "EPW needs to rebuild the coarse-grid Wannier NSCF save, but staged `nscf.in` was not found.",
                Some(
                    "Re-run the Wannier workflow so QCortado can preserve the NSCF input and save payload."
                        .to_string(),
                ),
            ));
        }
        commands.push(build_hpc_logged_qe_step_command(
            &profile,
            resource_type,
            "pw.x coarse-grid NSCF rebuild",
            "pw.x",
            None,
            "nscf.in",
            "nscf_rebuild.out",
        ));
        let rebuilt_schema_path = format!(
            "{}/{}.save/data-file-schema.xml",
            resolve_epw_stage_subdir(&config.input.outdir, "tmp"),
            config.input.prefix.trim()
        );
        commands.push(format!(
            "if [ ! -f {} ]; then echo \"[QCortado] ERROR: coarse-grid NSCF rebuild did not create {}\" >&2; exit 33; fi",
            shell_single_quote_local(&rebuilt_schema_path),
            rebuilt_schema_path
        ));
    }
    commands.push(build_hpc_logged_shell_step_command("epw.x", &epw_cmd));
    let bundle_files = vec![("epw.in".to_string(), input_text)];
    let bundle_copies = vec![(stage_dir.clone(), ".".to_string())];

    let config_line = format!(
        "[QCortado] EPW HPC config: phonon_source={}, wannier_source={}, prefix={}, coarse_k={}x{}x{}, fine_k={}x{}x{}, coarse_q={}x{}x{}, fine_q={}x{}x{}, mpi_ranks={}, npool={}, nimage={}, steps={}.",
        sources.phonon_calculation.id,
        sources
            .wannier_calculation
            .as_ref()
            .map(|calc| calc.id.as_str())
            .unwrap_or("none"),
        config.input.prefix,
        epw_coarse_k_mesh(&config.input)[0],
        epw_coarse_k_mesh(&config.input)[1],
        epw_coarse_k_mesh(&config.input)[2],
        epw_fine_k_mesh(&config.input)[0],
        epw_fine_k_mesh(&config.input)[1],
        epw_fine_k_mesh(&config.input)[2],
        epw_coarse_q_mesh(&config.input)[0],
        epw_coarse_q_mesh(&config.input)[1],
        epw_coarse_q_mesh(&config.input)[2],
        epw_fine_q_mesh(&config.input)[0],
        epw_fine_q_mesh(&config.input)[1],
        epw_fine_q_mesh(&config.input)[2],
        launch_ranks,
        parallel_plan.npool,
        parallel_plan.nimage,
        total_steps,
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &config_line);
    pm.append_output(task_id, config_line).await;
    if !parallel_plan.npool_explicit && launch_ranks > 1 {
        let line = format!(
            "[QCortado] No EPW npool was configured; defaulting to npool={} to match requested MPI ranks.",
            parallel_plan.npool
        );
        let _ = app.emit(&format!("task-output:{}", task_id), &line);
        pm.append_output(task_id, line).await;
    }
    for warning in &sources.warnings {
        let warning_line = format!("[QCortado] WARNING: {}", warning);
        let _ = app.emit(&format!("task-output:{}", task_id), &warning_line);
        pm.append_output(task_id, warning_line).await;
    }

    let work_path = run_hpc_bundle_task(
        app.clone(),
        pm.clone(),
        task_id,
        "epw",
        "EPW",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files,
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;
    let _ = std::fs::remove_dir_all(&stage_dir);

    let stdout = std::fs::read_to_string(work_path.join("epw.out")).unwrap_or_else(|_| {
        std::fs::read_to_string(work_path.join("slurm.out")).unwrap_or_default()
    });
    let stderr = std::fs::read_to_string(work_path.join("epw.err")).unwrap_or_else(|_| {
        std::fs::read_to_string(work_path.join("slurm.err")).unwrap_or_default()
    });
    let artifacts = collect_epw_artifacts(&work_path);
    let parsed_result = parse_epw_result_v2(
        &format!("{}\n{}", stdout, stderr),
        &work_path,
        artifacts.clone(),
        true,
    );
    let mut errors: Vec<EpwErrorRecord> = Vec::new();
    if parsed_result.summary.parse_partial {
        errors.push(EpwErrorRecord {
            code: "parse-partial".to_string(),
            message: "EPW output parse is partial.".to_string(),
            hint: Some("Inspect epw.out/slurm.out for full context.".to_string()),
        });
    }

    let done_line = format!(
        "[QCortado] EPW HPC pipeline complete in {:.1}s.",
        pipeline_start.elapsed().as_secs_f64()
    );
    let _ = app.emit(&format!("task-output:{}", task_id), &done_line);
    pm.append_output(task_id, done_line).await;

    Ok(EpwCalculationV1 {
        schema_version: EPW_SCHEMA_VERSION,
        sources: build_epw_sources_payload(&config, &sources),
        input: config.input,
        goals: config.goals,
        runtime: config.runtime,
        extensions: config.extensions,
        artifacts,
        result_summary: parsed_result.summary,
        transport: parsed_result.transport,
        superconductivity: parsed_result.superconductivity,
        spectral: parsed_result.spectral,
        parsed_tables: parsed_result.parsed_tables,
        errors,
    })
}

/// Starts an electronic DOS calculation as a background task.
#[tauri::command]
async fn start_dos_calculation(
    app: AppHandle,
    config: ElectronicDosCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let smearing_default = get_qe_smearing_default(&state);
    let (task_id, cancel_flag) = pm.register("dos".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();

        tokio::spawn(async move {
            let result = run_dos_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
                smearing_default.clone(),
                cancel_flag,
                pm.clone(),
            )
            .await;

            match result {
                Ok(dos_data) => {
                    let json = serde_json::to_value(&dos_data).unwrap_or(serde_json::Value::Null);
                    pm.complete(&tid, json).await;
                    let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
                }
                Err(e) => {
                    pm.fail(&tid, e.clone()).await;
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let tid = task_id.clone();
    let app_handle = app.clone();

    tokio::spawn(async move {
        let result = run_dos_background(
            app_handle.clone(),
            &tid,
            config,
            working_dir,
            mpi_config,
            bin_dir,
            execution_prefix,
            smearing_default,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(dos_data) => {
                let json = serde_json::to_value(&dos_data).unwrap_or(serde_json::Value::Null);
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
async fn run_dos_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: ElectronicDosCalculationConfig,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    smearing_default: qe::SmearingType,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<ElectronicDosData, String> {
    let mut nscf_calc = config.base_calculation.clone();
    nscf_calc.pseudo_dir = profile.remote_pseudo_dir.clone();
    nscf_calc.calculation = qe::CalculationType::Nscf;
    nscf_calc.verbosity = Some("high".to_string());
    nscf_calc.kpoints = qe::KPoints::Automatic {
        grid: config.k_grid,
        offset: [0, 0, 0],
    };
    if nscf_calc.system.degauss.is_none() {
        nscf_calc.system.degauss = config.degauss;
    }
    if matches!(nscf_calc.system.occupations, qe::Occupations::Fixed) {
        nscf_calc.system.occupations = qe::Occupations::Smearing;
        nscf_calc.system.smearing = smearing_default;
        if nscf_calc.system.degauss.is_none() {
            nscf_calc.system.degauss = Some(0.02);
        }
    }

    let nscf_input = generate_pw_input(&nscf_calc);
    let dos_calc = DosCalculation {
        prefix: nscf_calc.prefix.clone(),
        outdir: nscf_calc.outdir.clone(),
        fildos: "dos.dat".to_string(),
        degauss: config.degauss.or(nscf_calc.system.degauss),
        emin: config.emin,
        emax: config.emax,
        delta_e: config.delta_e,
    };
    let dos_input = generate_dos_input(&dos_calc);

    let dependency_stage = resolve_hpc_scf_dependency_stage(
        &app,
        config.project_id.as_deref(),
        config.scf_calc_id.as_deref(),
    )?;

    let mut bundle_copies: Vec<(PathBuf, String)> = Vec::new();
    if let Some(local_scf_tmp_dir) = dependency_stage.local_bundle_copy {
        bundle_copies.push((local_scf_tmp_dir, ".".to_string()));
    }

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let mut commands = vec!["cd \"$SLURM_SUBMIT_DIR\"".to_string()];
    commands.extend(dependency_stage.remote_hydration_commands);
    commands.push(format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)));
    commands.push(build_hpc_qe_input_command(
        &profile,
        resource_type,
        "pw.x",
        None,
        "nscf.in",
        "nscf.out",
    ));
    commands.push(build_hpc_qe_input_command(
        &profile,
        resource_type,
        "dos.x",
        None,
        "dos.in",
        "dos.out",
    ));
    let bundle_files = vec![
        ("nscf.in".to_string(), nscf_input),
        ("dos.in".to_string(), dos_input),
    ];
    let work_path = run_hpc_bundle_task(
        app,
        pm,
        task_id,
        "dos",
        "Electronic DOS",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files,
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;

    let dos_file = work_path.join(&dos_calc.fildos);
    if !dos_file.exists() {
        return Err(format!(
            "DOS file not found after HPC run: {}",
            dos_file.display()
        ));
    }
    let dos_content = std::fs::read_to_string(&dos_file)
        .map_err(|e| format!("Failed to read DOS file {}: {}", dos_file.display(), e))?;
    let (energies, dos_values) = parse_dos_file(&dos_content)
        .ok_or_else(|| format!("Failed to parse DOS data from {}", dos_file.display()))?;
    if energies.is_empty() || dos_values.is_empty() {
        return Err("Parsed DOS output is empty".to_string());
    }
    let e_min = energies
        .iter()
        .fold(f64::INFINITY, |current, value| current.min(*value));
    let e_max = energies
        .iter()
        .fold(f64::NEG_INFINITY, |current, value| current.max(*value));
    let max_dos = dos_values
        .iter()
        .fold(f64::NEG_INFINITY, |current, value| current.max(*value));
    let fermi_energy = std::fs::read_to_string(work_path.join("nscf.out"))
        .ok()
        .and_then(|content| extract_fermi_energy_from_text(&content));

    Ok(ElectronicDosData {
        energies,
        dos: dos_values,
        fermi_energy,
        energy_range: [e_min, e_max],
        max_dos: if max_dos.is_finite() {
            max_dos.max(0.0)
        } else {
            0.0
        },
        points: dos_content
            .lines()
            .filter(|line| {
                let trimmed = line.trim();
                !trimmed.is_empty() && !trimmed.starts_with('#')
            })
            .count(),
    })
}

#[allow(clippy::too_many_arguments)]
async fn run_dos_background(
    app: AppHandle,
    task_id: &str,
    config: ElectronicDosCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    bin_dir: PathBuf,
    execution_prefix: Option<String>,
    smearing_default: qe::SmearingType,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<ElectronicDosData, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let work_path = PathBuf::from(&working_dir);
    prepare_working_directory(&work_path, false)?;

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
            emit_line!(format!("Copying SCF data from: {}", scf_tmp_dir.display()));
            projects::copy_dir_contents(&scf_tmp_dir, &work_path)?;
            emit_line!("SCF data copied successfully.".to_string());
        } else {
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    check_cancel!();

    let mut nscf_calc = config.base_calculation.clone();
    nscf_calc.calculation = qe::CalculationType::Nscf;
    nscf_calc.verbosity = Some("high".to_string());
    nscf_calc.kpoints = qe::KPoints::Automatic {
        grid: config.k_grid,
        offset: [0, 0, 0],
    };

    if nscf_calc.system.degauss.is_none() {
        nscf_calc.system.degauss = config.degauss;
    }
    if matches!(nscf_calc.system.occupations, qe::Occupations::Fixed) {
        nscf_calc.system.occupations = qe::Occupations::Smearing;
        nscf_calc.system.smearing = smearing_default;
        if nscf_calc.system.degauss.is_none() {
            nscf_calc.system.degauss = Some(0.02);
        }
    }

    let nscf_input = generate_pw_input(&nscf_calc);
    std::fs::write(work_path.join("nscf.in"), &nscf_input)
        .map_err(|e| format!("Failed to write NSCF input: {}", e))?;

    emit_line!("".to_string());
    emit_line!("=== Starting Electronic DOS Calculation ===".to_string());
    emit_line!(format!(
        "Dense k-grid: {}×{}×{}",
        config.k_grid[0], config.k_grid[1], config.k_grid[2]
    ));
    emit_line!("Step 1/2: Running NSCF on dense k-grid...".to_string());

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
                .args(["-pd", ".true."])
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
        stdin
            .write_all(nscf_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write NSCF input: {}", e))?;
    }

    let stdout = child.stdout.take().ok_or("Failed to capture pw.x stdout")?;
    let mut reader = BufReader::new(stdout).lines();
    let mut full_nscf_output = String::new();
    let mut fermi_energy: Option<f64> = None;

    while let Some(line) = reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        full_nscf_output.push_str(&line);
        full_nscf_output.push('\n');
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
        return Err(format!(
            "pw.x (NSCF) failed with exit code: {:?}",
            status.code()
        ));
    }

    std::fs::write(work_path.join("nscf.out"), &full_nscf_output)
        .map_err(|e| format!("Failed to write NSCF output: {}", e))?;

    emit_line!("".to_string());
    emit_line!("Step 2/2: Running dos.x post-processing...".to_string());

    let dos_path = bin_dir.join("dos.x");
    if !dos_path.exists() {
        return Err("dos.x not found. Make sure your QE installation includes dos.x".to_string());
    }

    let dos_calc = DosCalculation {
        prefix: nscf_calc.prefix.clone(),
        outdir: nscf_calc.outdir.clone(),
        fildos: "dos.dat".to_string(),
        degauss: config.degauss.or(nscf_calc.system.degauss),
        emin: config.emin,
        emax: config.emax,
        delta_e: config.delta_e,
    };

    let dos_input = generate_dos_input(&dos_calc);
    std::fs::write(work_path.join("dos.in"), &dos_input)
        .map_err(|e| format!("Failed to write dos.x input: {}", e))?;

    let mut dos_child = tokio_command_with_prefix(&dos_path, execution_prefix.as_deref())
        .current_dir(&work_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start dos.x: {}", e))?;

    if let Some(pid) = dos_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = dos_child.stdin.take() {
        stdin
            .write_all(dos_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write dos.x input: {}", e))?;
    }

    let dos_stdout = dos_child
        .stdout
        .take()
        .ok_or("Failed to capture dos.x stdout")?;
    let mut dos_reader = BufReader::new(dos_stdout).lines();
    let mut dos_output = String::new();

    while let Some(line) = dos_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        dos_output.push_str(&line);
        dos_output.push('\n');
        emit_line!(line);
    }

    let dos_status = dos_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !dos_status.success() {
        return Err(format!(
            "dos.x failed with exit code: {:?}",
            dos_status.code()
        ));
    }

    std::fs::write(work_path.join("dos.out"), &dos_output)
        .map_err(|e| format!("Failed to write dos.x output: {}", e))?;

    emit_line!("".to_string());
    emit_line!("Parsing DOS data...".to_string());

    let dos_file = work_path.join(&dos_calc.fildos);
    if !dos_file.exists() {
        return Err(format!("DOS file not found: {}", dos_file.display()));
    }

    let dos_content = std::fs::read_to_string(&dos_file)
        .map_err(|e| format!("Failed to read DOS file {}: {}", dos_file.display(), e))?;
    let (energies, dos_values) = parse_dos_file(&dos_content)
        .ok_or_else(|| format!("Failed to parse DOS data from {}", dos_file.display()))?;

    if energies.is_empty() || dos_values.is_empty() {
        return Err("Parsed DOS output is empty".to_string());
    }

    let e_min = energies
        .iter()
        .fold(f64::INFINITY, |current, value| current.min(*value));
    let e_max = energies
        .iter()
        .fold(f64::NEG_INFINITY, |current, value| current.max(*value));
    let max_dos = dos_values
        .iter()
        .fold(f64::NEG_INFINITY, |current, value| current.max(*value));
    let points = energies.len();

    let dos_data = ElectronicDosData {
        energies,
        dos: dos_values,
        fermi_energy,
        energy_range: [e_min, e_max],
        max_dos: if max_dos.is_finite() {
            max_dos.max(0.0)
        } else {
            0.0
        },
        points,
    };

    emit_line!("=== Electronic DOS Complete ===".to_string());
    emit_line!(format!(
        "  {} points, energy range [{:.3}, {:.3}] eV",
        dos_data.points, dos_data.energy_range[0], dos_data.energy_range[1]
    ));

    Ok(dos_data)
}

/// Starts a Fermi-surface generation task using QE's fermi_velocity.x workflow.
#[tauri::command]
async fn start_fermi_surface_calculation(
    app: AppHandle,
    config: FermiSurfaceCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let smearing_default = get_qe_smearing_default(&state);
    let (task_id, cancel_flag) = pm.register("fermi_surface".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();

        tokio::spawn(async move {
            let result = run_fermi_surface_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
                smearing_default.clone(),
                cancel_flag,
                pm.clone(),
            )
            .await;

            match result {
                Ok(fermi_data) => {
                    let json = serde_json::to_value(&fermi_data).unwrap_or(serde_json::Value::Null);
                    pm.complete(&tid, json).await;
                    let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
                }
                Err(e) => {
                    pm.fail(&tid, e.clone()).await;
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

    let tid = task_id.clone();
    let app_handle = app.clone();

    tokio::spawn(async move {
        let result = run_fermi_surface_background(
            app_handle.clone(),
            &tid,
            config,
            working_dir,
            mpi_config,
            bin_dir,
            execution_prefix,
            smearing_default,
            cancel_flag,
            pm.clone(),
        )
        .await;

        match result {
            Ok(fermi_data) => {
                let json = serde_json::to_value(&fermi_data).unwrap_or(serde_json::Value::Null);
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
async fn run_fermi_surface_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: FermiSurfaceCalculationConfig,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    smearing_default: qe::SmearingType,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<FermiSurfaceData, String> {
    let mut nscf_calc = config.base_calculation.clone();
    nscf_calc.pseudo_dir = profile.remote_pseudo_dir.clone();
    nscf_calc.calculation = qe::CalculationType::Nscf;
    if nscf_calc.verbosity.is_none() {
        nscf_calc.verbosity = Some("high".to_string());
    }
    let k_offset = config
        .k_offset
        .unwrap_or([0, 0, 0])
        .map(|value| u32::from(value > 0));
    nscf_calc.kpoints = qe::KPoints::Automatic {
        grid: config.k_grid,
        offset: k_offset,
    };
    if nscf_calc.system.degauss.is_none() {
        nscf_calc.system.degauss = Some(0.02);
    }
    if matches!(nscf_calc.system.occupations, qe::Occupations::Fixed) {
        nscf_calc.system.occupations = qe::Occupations::Smearing;
        nscf_calc.system.smearing = smearing_default;
        if nscf_calc.system.degauss.is_none() {
            nscf_calc.system.degauss = Some(0.02);
        }
    }

    let mut nscf_input = generate_pw_input(&nscf_calc);
    if let Some(nbnd) = config.nbnd {
        nscf_input = insert_system_namelist_line(&nscf_input, &format!("nbnd = {},", nbnd))?;
    }
    let fermi_input = nscf_input.clone();

    let dependency_stage = resolve_hpc_scf_dependency_stage(
        &app,
        config.project_id.as_deref(),
        config.scf_calc_id.as_deref(),
    )?;

    let mut bundle_copies: Vec<(PathBuf, String)> = Vec::new();
    if let Some(local_scf_tmp_dir) = dependency_stage.local_bundle_copy {
        bundle_copies.push((local_scf_tmp_dir, ".".to_string()));
    }

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let mut commands = vec!["cd \"$SLURM_SUBMIT_DIR\"".to_string()];
    commands.extend(dependency_stage.remote_hydration_commands);
    commands.push(format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)));
    commands.push(build_hpc_qe_input_command(
        &profile,
        resource_type,
        "pw.x",
        None,
        "nscf.in",
        "nscf.out",
    ));
    commands.push(build_hpc_qe_input_command(
        &profile,
        resource_type,
        "fermi_velocity.x",
        Some("-npool 1"),
        "fermi_velocity.in",
        "fermi_velocity.out",
    ));
    commands.push(
        "if [ -d ./tmp ]; then find ./tmp -type f \\( -iname '*.frmsf' -o -iname '*.bxsf' \\) -exec cp -n {} . \\; || true; rm -rf ./tmp; fi"
            .to_string(),
    );
    let bundle_files = vec![
        ("nscf.in".to_string(), nscf_input),
        ("fermi_velocity.in".to_string(), fermi_input),
    ];
    let work_path = run_hpc_bundle_task(
        app,
        pm,
        task_id,
        "fermi_surface",
        "Fermi Surface",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files,
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;

    let nscf_output = std::fs::read_to_string(work_path.join("nscf.out")).unwrap_or_default();
    let fermi_output = std::fs::read_to_string(work_path.join("fermi_velocity.out"))
        .unwrap_or_else(|_| {
            std::fs::read_to_string(work_path.join("slurm.out")).unwrap_or_default()
        });
    if let Some(missing_points) = extract_unmapped_fermi_grid_points(&fermi_output) {
        if missing_points > 0 {
            return Err(format!(
                "fermi_velocity.x reported {} unmapped k-grid points. Ensure the NSCF step and fermi_velocity input share identical K_POINTS settings.",
                missing_points
            ));
        }
    }
    let fermi_energy = extract_fermi_energy_from_text(&nscf_output)
        .or_else(|| extract_fermi_energy_from_text(&fermi_output));
    let frmsf_files = collect_frmsf_files(&work_path)?;
    if frmsf_files.is_empty() {
        return Err("No .frmsf files were generated by HPC run.".to_string());
    }
    let primary_file = frmsf_files
        .iter()
        .find(|file| {
            let lower = file.file_name.to_ascii_lowercase();
            lower.ends_with("/vfermi.frmsf")
                || lower == "vfermi.frmsf"
                || lower.ends_with("_vfermi.frmsf")
        })
        .or_else(|| {
            frmsf_files.iter().find(|file| {
                let lower = file.file_name.to_ascii_lowercase();
                lower.ends_with("/vfermi1.frmsf")
                    || lower == "vfermi1.frmsf"
                    || lower.ends_with("_vfermi1.frmsf")
            })
        })
        .map(|file| file.file_name.clone())
        .unwrap_or_else(|| frmsf_files[0].file_name.clone());

    Ok(FermiSurfaceData {
        k_grid: config.k_grid,
        fermi_energy,
        primary_file,
        frmsf_files,
    })
}

#[allow(clippy::too_many_arguments)]
async fn run_fermi_surface_background(
    app: AppHandle,
    task_id: &str,
    config: FermiSurfaceCalculationConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    bin_dir: PathBuf,
    execution_prefix: Option<String>,
    smearing_default: qe::SmearingType,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<FermiSurfaceData, String> {
    use std::process::Stdio;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    let work_path = PathBuf::from(&working_dir);
    prepare_working_directory(&work_path, false)?;

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

    // Copy SCF .save directory if provided.
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
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    check_cancel!();

    let mut nscf_calc = config.base_calculation.clone();
    nscf_calc.calculation = qe::CalculationType::Nscf;
    if nscf_calc.verbosity.is_none() {
        nscf_calc.verbosity = Some("high".to_string());
    }
    let k_offset = config
        .k_offset
        .unwrap_or([0, 0, 0])
        .map(|value| u32::from(value > 0));
    nscf_calc.kpoints = qe::KPoints::Automatic {
        grid: config.k_grid,
        offset: k_offset,
    };

    if nscf_calc.system.degauss.is_none() {
        nscf_calc.system.degauss = Some(0.02);
    }
    if matches!(nscf_calc.system.occupations, qe::Occupations::Fixed) {
        nscf_calc.system.occupations = qe::Occupations::Smearing;
        nscf_calc.system.smearing = smearing_default;
        if nscf_calc.system.degauss.is_none() {
            nscf_calc.system.degauss = Some(0.02);
        }
    }

    let mut nscf_input = generate_pw_input(&nscf_calc);
    if let Some(nbnd) = config.nbnd {
        nscf_input = insert_system_namelist_line(&nscf_input, &format!("nbnd = {},", nbnd))?;
        emit_line!(format!("Requested nbnd = {}", nbnd));
    }
    let fermi_input = nscf_input.clone();
    std::fs::write(work_path.join("nscf.in"), &nscf_input)
        .map_err(|e| format!("Failed to write NSCF input: {}", e))?;
    std::fs::write(work_path.join("fermi_velocity.in"), &fermi_input)
        .map_err(|e| format!("Failed to write fermi_velocity.x input: {}", e))?;

    emit_line!("".to_string());
    emit_line!("=== Starting Fermi Surface Generation ===".to_string());
    emit_line!(format!(
        "K-grid: {}×{}×{}",
        config.k_grid[0], config.k_grid[1], config.k_grid[2]
    ));
    emit_line!(format!(
        "K-point offset: {} {} {}",
        k_offset[0], k_offset[1], k_offset[2]
    ));
    emit_line!("Step 1/2: Running pw.x (NSCF) on Fermi grid...".to_string());

    let pw_path = bin_dir.join("pw.x");
    if !pw_path.exists() {
        return Err("pw.x not found".to_string());
    }

    let mut nscf_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!("Using MPI with {} processes", mpi.nprocs));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&pw_path)
                .args(["-pd", ".true."])
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio_command_with_prefix(&pw_path, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start pw.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&pw_path, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start pw.x: {}", e))?
    };

    if let Some(pid) = nscf_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    if let Some(mut stdin) = nscf_child.stdin.take() {
        stdin
            .write_all(nscf_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write NSCF input: {}", e))?;
    }

    let nscf_stdout = nscf_child
        .stdout
        .take()
        .ok_or("Failed to capture pw.x stdout")?;
    let mut nscf_reader = BufReader::new(nscf_stdout).lines();
    let mut nscf_output = String::new();

    while let Some(line) = nscf_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        nscf_output.push_str(&line);
        nscf_output.push('\n');
        emit_line!(line);
    }

    let nscf_status = nscf_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !nscf_status.success() {
        return Err(format!(
            "pw.x (NSCF) failed with exit code: {:?}",
            nscf_status.code()
        ));
    }

    std::fs::write(work_path.join("nscf.out"), &nscf_output)
        .map_err(|e| format!("Failed to write NSCF output: {}", e))?;

    let mut fermi_energy = extract_fermi_energy_from_text(&nscf_output);
    emit_line!("Step 2/2: Running fermi_velocity.x...".to_string());

    let fermi_velocity_path = bin_dir.join("fermi_velocity.x");
    if !fermi_velocity_path.exists() {
        return Err(
            "fermi_velocity.x not found. Build QE post-processing tools (`make pp`) and verify the QE bin path."
                .to_string(),
        );
    }

    let mut fermi_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!("Using MPI with {} processes (npool=1)", mpi.nprocs));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&fermi_velocity_path)
                .args(["-npool", "1", "-in", "fermi_velocity.in"])
                .current_dir(&work_path)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun: {}", e))?
        } else {
            tokio_command_with_prefix(&fermi_velocity_path, execution_prefix.as_deref())
                .args(["-npool", "1", "-in", "fermi_velocity.in"])
                .current_dir(&work_path)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start fermi_velocity.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&fermi_velocity_path, execution_prefix.as_deref())
            .args(["-npool", "1", "-in", "fermi_velocity.in"])
            .current_dir(&work_path)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start fermi_velocity.x: {}", e))?
    };

    if let Some(pid) = fermi_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

    let fermi_stdout = fermi_child
        .stdout
        .take()
        .ok_or("Failed to capture fermi_velocity.x stdout")?;
    let mut fermi_reader = BufReader::new(fermi_stdout).lines();
    let mut fermi_output = String::new();

    while let Some(line) = fermi_reader.next_line().await.map_err(|e| e.to_string())? {
        check_cancel!();
        fermi_output.push_str(&line);
        fermi_output.push('\n');
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

    let fermi_status = fermi_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !fermi_status.success() {
        return Err(format!(
            "fermi_velocity.x failed with exit code: {:?}",
            fermi_status.code()
        ));
    }

    std::fs::write(work_path.join("fermi_velocity.out"), &fermi_output)
        .map_err(|e| format!("Failed to write fermi_velocity.x output: {}", e))?;
    if let Some(missing_points) = extract_unmapped_fermi_grid_points(&fermi_output) {
        if missing_points > 0 {
            return Err(format!(
                "fermi_velocity.x reported {} unmapped k-grid points. Ensure the NSCF step and fermi_velocity input share identical K_POINTS settings.",
                missing_points
            ));
        }
    }

    emit_line!("".to_string());
    emit_line!("Collecting FRMSF artifacts...".to_string());

    let frmsf_files = collect_frmsf_files(&work_path)?;
    if frmsf_files.is_empty() {
        return Err(
            "No .frmsf files were generated. fermi_velocity.x completed without producing expected output."
                .to_string(),
        );
    }

    let primary_file = frmsf_files
        .iter()
        .find(|file| {
            let lower = file.file_name.to_ascii_lowercase();
            lower.ends_with("/vfermi.frmsf")
                || lower == "vfermi.frmsf"
                || lower.ends_with("_vfermi.frmsf")
        })
        .or_else(|| {
            frmsf_files.iter().find(|file| {
                let lower = file.file_name.to_ascii_lowercase();
                lower.ends_with("/vfermi1.frmsf")
                    || lower == "vfermi1.frmsf"
                    || lower.ends_with("_vfermi1.frmsf")
            })
        })
        .map(|file| file.file_name.clone())
        .unwrap_or_else(|| frmsf_files[0].file_name.clone());

    emit_line!("=== Fermi Surface Generation Complete ===".to_string());
    emit_line!(format!("  Generated {} FRMSF file(s)", frmsf_files.len()));
    emit_line!(format!("  Primary file: {}", primary_file));

    Ok(FermiSurfaceData {
        k_grid: config.k_grid,
        fermi_energy,
        primary_file,
        frmsf_files,
    })
}

/// Starts a phonon calculation as a background task.
#[tauri::command]
async fn start_phonon_calculation(
    app: AppHandle,
    config: PhononPipelineConfig,
    working_dir: String,
    mpi_config: Option<MpiConfig>,
    execution_target: Option<hpc::profile::ExecutionTarget>,
    label: String,
    state: State<'_, AppState>,
) -> Result<String, String> {
    let hpc_target = resolve_hpc_execution(&state, execution_target.as_ref());

    if hpc_target.is_none() && state.process_manager.has_running_tasks().await {
        return Err(
            "A calculation is already running. Please wait for it to complete or cancel it."
                .to_string(),
        );
    }

    let pm = state.process_manager.clone();
    let (task_id, cancel_flag) = pm.register("phonon".to_string(), label).await;

    if let Some(hpc_target) = hpc_target {
        let profile = resolve_hpc_profile_from_state(&state, hpc_target.profile_id.clone())?;
        let secret = hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        )?;
        let tid = task_id.clone();
        let app_handle = app.clone();

        tokio::spawn(async move {
            let result = run_phonon_hpc_background(
                app_handle.clone(),
                &tid,
                config,
                working_dir,
                profile,
                secret,
                hpc_target.resources,
                hpc_target.recovery_save,
                cancel_flag,
                pm.clone(),
            )
            .await;

            match result {
                Ok(phonon_result) => {
                    let json =
                        serde_json::to_value(&phonon_result).unwrap_or(serde_json::Value::Null);
                    pm.complete(&tid, json).await;
                    let _ = app_handle.emit(&format!("task-complete:{}", tid), "completed");
                }
                Err(e) => {
                    pm.fail(&tid, e.clone()).await;
                    let _ =
                        app_handle.emit(&format!("task-status:{}", tid), &format!("failed:{}", e));
                }
            }
        });
        return Ok(task_id);
    }

    let bin_dir = {
        let guard = state.qe_bin_dir.lock().unwrap();
        guard.as_ref().ok_or("QE path not configured")?.clone()
    };
    let execution_prefix = state.execution_prefix.lock().unwrap().clone();

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
async fn run_phonon_hpc_background(
    app: AppHandle,
    task_id: &str,
    config: PhononPipelineConfig,
    working_dir: String,
    profile: hpc::profile::HpcProfile,
    secret: Option<String>,
    resources: Option<hpc::profile::SlurmResourceRequest>,
    recovery_save: Option<hpc::profile::HpcRecoverySaveSpec>,
    cancel_flag: std::sync::Arc<std::sync::atomic::AtomicBool>,
    pm: ProcessManager,
) -> Result<PhononResult, String> {
    let dependency_stage = resolve_hpc_scf_dependency_stage(
        &app,
        config.project_id.as_deref(),
        config.scf_calc_id.as_deref(),
    )?;

    let mut bundle_copies: Vec<(PathBuf, String)> = Vec::new();
    if let Some(local_scf_tmp_dir) = dependency_stage.local_bundle_copy.clone() {
        bundle_copies.push((local_scf_tmp_dir, ".".to_string()));
    }

    let q2r_calc = Q2RCalculation {
        fildyn: config.phonon.fildyn.clone(),
        flfrc: "force_constants".to_string(),
        zasr: config.phonon.asr.clone(),
    };
    let q2r_input = generate_q2r_input(&q2r_calc);

    let q_path_with_points = config.q_path.as_ref().map(|q_path| {
        q_path
            .iter()
            .enumerate()
            .map(|(idx, point)| QPathPoint {
                label: point.label.clone(),
                coords: point.coords,
                npoints: if idx < q_path.len() - 1 {
                    point.npoints
                } else {
                    0
                },
            })
            .collect::<Vec<QPathPoint>>()
    });
    let matdyn_dos_input = if config.calculate_dos {
        let dos_grid = config.dos_grid.unwrap_or([20, 20, 20]);
        let dos_delta_e = config.dos_delta_e.unwrap_or(1.0);
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
        Some(generate_matdyn_dos_input(&matdyn_dos_calc))
    } else {
        None
    };
    let matdyn_bands_input = if config.calculate_dispersion {
        if let Some(q_path) = q_path_with_points.as_ref() {
            let matdyn_bands_calc = MatdynCalculation {
                flfrc: "force_constants".to_string(),
                asr: config.phonon.asr.clone(),
                dos: false,
                fldos: None,
                nk: None,
                delta_e: None,
                q_path: Some(q_path.clone()),
                flfrq: Some("phonon_freq".to_string()),
            };
            Some(generate_matdyn_bands_input(&matdyn_bands_calc))
        } else {
            None
        }
    } else {
        None
    };

    let resource_type = resolve_hpc_resource_type_for_resources(&profile, resources.as_ref());
    let qe_bin_dir = resolve_hpc_qe_bin_dir_for_resources(&profile, resources.as_ref());
    let mut commands = vec!["cd \"$SLURM_SUBMIT_DIR\"".to_string()];
    commands.extend(dependency_stage.remote_hydration_commands);
    commands.push(format!("QE_BIN={}", shell_single_quote_local(&qe_bin_dir)));
    commands.push(build_hpc_qe_input_command(
        &profile,
        resource_type,
        "ph.x",
        None,
        "ph.in",
        "ph.out",
    ));
    commands.push(build_hpc_qe_input_command(
        &profile,
        resource_type,
        "q2r.x",
        None,
        "q2r.in",
        "q2r.out",
    ));
    if config.calculate_dos {
        let matdyn_dos_cmd = build_hpc_qe_input_command(
            &profile,
            resource_type,
            "matdyn.x",
            None,
            "matdyn_dos.in",
            "matdyn_dos.out",
        );
        commands.push(format!(
            "if ! {}; then echo \"WARNING: matdyn.x DOS failed\"; fi",
            matdyn_dos_cmd
        ));
    }
    if config.calculate_dispersion && matdyn_bands_input.is_some() {
        let matdyn_bands_cmd = build_hpc_qe_input_command(
            &profile,
            resource_type,
            "matdyn.x",
            None,
            "matdyn_bands.in",
            "matdyn_bands.out",
        );
        commands.push(format!(
            "if ! {}; then echo \"WARNING: matdyn.x dispersion failed\"; fi",
            matdyn_bands_cmd
        ));
    }

    let mut bundle_files = vec![
        ("ph.in".to_string(), generate_ph_input(&config.phonon)),
        ("q2r.in".to_string(), q2r_input),
    ];
    if let Some(input) = matdyn_dos_input {
        bundle_files.push(("matdyn_dos.in".to_string(), input));
    }
    if let Some(input) = matdyn_bands_input {
        bundle_files.push(("matdyn_bands.in".to_string(), input));
    }

    let work_path = run_hpc_bundle_task(
        app,
        pm,
        task_id,
        "phonon",
        "Phonon",
        profile,
        secret,
        resources,
        &working_dir,
        commands,
        bundle_files,
        bundle_copies,
        recovery_save,
        cancel_flag,
    )
    .await?;

    // Do not re-run local `.save` validation here: minimal HPC artifact sync intentionally omits
    // heavy `tmp/*.save` payloads after the remote run has already completed.
    let ph_output = std::fs::read_to_string(work_path.join("ph.out")).unwrap_or_else(|_| {
        std::fs::read_to_string(work_path.join("slurm.out")).unwrap_or_default()
    });
    let (converged, n_qpoints) = parse_ph_output(&ph_output);
    if !converged {
        return Err("ph.x did not converge successfully".to_string());
    }

    let dos_data = {
        let dos_file = work_path.join("phonon_dos");
        if dos_file.exists() {
            read_phonon_dos_file(&dos_file).ok()
        } else {
            None
        }
    };

    let mut dispersion_data = {
        let freq_gp_file = work_path.join("phonon_freq.gp");
        let freq_file = work_path.join("phonon_freq");
        let source = if freq_gp_file.exists() {
            Some(freq_gp_file)
        } else if freq_file.exists() {
            Some(freq_file)
        } else {
            None
        };
        if let Some(source_file) = source {
            read_phonon_dispersion_file(&source_file).ok()
        } else {
            None
        }
    };

    if let (Some(ref mut disp), Some(q_path)) =
        (dispersion_data.as_mut(), q_path_with_points.as_ref())
    {
        add_phonon_symmetry_markers(disp, q_path);
    }

    let n_modes = dispersion_data
        .as_ref()
        .map(|value| value.n_modes)
        .unwrap_or(0);

    Ok(PhononResult {
        converged: true,
        n_qpoints,
        n_modes,
        dos_data,
        dispersion_data,
        raw_output: ph_output,
    })
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
    // Keep existing scratch only when explicit recover mode is enabled.
    prepare_working_directory(&work_path, config.phonon.recover)?;

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
            return Err(missing_scf_tmp_error(&scf_tmp_dir));
        }
    }

    ensure_phonon_restart_inputs(&work_path)?;

    check_cancel!();

    // Step 1: ph.x
    emit_line!("".to_string());
    emit_line!("=== Step 1/4: Running ph.x Phonon Calculation ===".to_string());
    emit_line!(format!(
        "Q-grid: {}×{}×{}",
        config.phonon.nq[0], config.phonon.nq[1], config.phonon.nq[2]
    ));

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
        stdin
            .write_all(ph_input.as_bytes())
            .await
            .map_err(|e| format!("Failed to write ph.x input: {}", e))?;
    }

    let ph_stdout = ph_child
        .stdout
        .take()
        .ok_or("Failed to capture ph.x stdout")?;
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
            return Err(format!(
                "ph.x failed with exit code: {:?}",
                ph_status.code()
            ));
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

    let mut q2r_child = if let Some(ref mpi) = mpi_config {
        if mpi.enabled && mpi.nprocs > 1 {
            emit_line!(format!("Running q2r.x with MPI ({} processes)", mpi.nprocs));
            tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                .args(["-np", &mpi.nprocs.to_string()])
                .arg(&q2r_exe)
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start mpirun for q2r.x: {}", e))?
        } else {
            tokio_command_with_prefix(&q2r_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start q2r.x: {}", e))?
        }
    } else {
        tokio_command_with_prefix(&q2r_exe, execution_prefix.as_deref())
            .current_dir(&work_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start q2r.x: {}", e))?
    };

    if let Some(pid) = q2r_child.id() {
        pm.set_child_id(task_id, pid).await;
    }

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
        check_cancel!();
        emit_line!(line);
    }

    let q2r_status = q2r_child.wait().await.map_err(|e| e.to_string())?;
    check_cancel!();
    if !q2r_status.success() {
        return Err(format!(
            "q2r.x failed with exit code: {:?}",
            q2r_status.code()
        ));
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
        emit_line!(format!(
            "DOS grid: {}×{}×{}",
            dos_grid[0], dos_grid[1], dos_grid[2]
        ));
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

        let mut matdyn_dos_child = if let Some(ref mpi) = mpi_config {
            if mpi.enabled && mpi.nprocs > 1 {
                emit_line!(format!(
                    "Running matdyn.x (DOS) with MPI ({} processes)",
                    mpi.nprocs
                ));
                tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                    .args(["-np", &mpi.nprocs.to_string()])
                    .arg(&matdyn_exe)
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start mpirun for matdyn.x DOS: {}", e))?
            } else {
                tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start matdyn.x for DOS: {}", e))?
            }
        } else {
            tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                .current_dir(&work_path)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to start matdyn.x for DOS: {}", e))?
        };

        if let Some(pid) = matdyn_dos_child.id() {
            pm.set_child_id(task_id, pid).await;
        }

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
                            dos.frequencies.len(),
                            dos.omega_min,
                            dos.omega_max
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
                    npoints: if i < q_path.len() - 1 { p.npoints } else { 0 },
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

            let mut matdyn_bands_child = if let Some(ref mpi) = mpi_config {
                if mpi.enabled && mpi.nprocs > 1 {
                    emit_line!(format!(
                        "Running matdyn.x (dispersion) with MPI ({} processes)",
                        mpi.nprocs
                    ));
                    tokio_command_with_prefix("mpirun", execution_prefix.as_deref())
                        .args(["-np", &mpi.nprocs.to_string()])
                        .arg(&matdyn_exe)
                        .current_dir(&work_path)
                        .stdin(Stdio::piped())
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn()
                        .map_err(|e| {
                            format!("Failed to start mpirun for matdyn.x dispersion: {}", e)
                        })?
                } else {
                    tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                        .current_dir(&work_path)
                        .stdin(Stdio::piped())
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn()
                        .map_err(|e| format!("Failed to start matdyn.x for dispersion: {}", e))?
                }
            } else {
                tokio_command_with_prefix(&matdyn_exe, execution_prefix.as_deref())
                    .current_dir(&work_path)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .map_err(|e| format!("Failed to start matdyn.x for dispersion: {}", e))?
            };

            if let Some(pid) = matdyn_bands_child.id() {
                pm.set_child_id(task_id, pid).await;
            }

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
                            emit_line!(format!(
                                "Warning: Failed to parse phonon dispersion: {}",
                                e
                            ));
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
async fn cancel_task(task_id: String, state: State<'_, AppState>) -> Result<(), String> {
    state.process_manager.cancel(&task_id).await
}

#[tauri::command]
async fn dismiss_task(task_id: String, state: State<'_, AppState>) -> Result<(), String> {
    state.process_manager.remove(&task_id).await;
    Ok(())
}

async fn cancel_running_hpc_jobs_for_quit(state: &AppState) -> Result<(), String> {
    let running_hpc_jobs = state.process_manager.list_running_hpc_jobs().await;
    if running_hpc_jobs.is_empty() {
        return Ok(());
    }

    let mut errors: Vec<String> = Vec::new();
    for job in running_hpc_jobs {
        let Some(remote_job_id) = job
            .remote_job_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            continue;
        };
        let Some(profile_id) = job
            .hpc_profile_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            errors.push(format!(
                "Task {} has remote job {} but no HPC profile context.",
                job.task_id, remote_job_id
            ));
            continue;
        };

        let profile = match resolve_hpc_profile_from_state(state, Some(profile_id.to_string())) {
            Ok(profile) => profile,
            Err(e) => {
                errors.push(format!(
                    "Task {} (job {}) profile lookup failed: {}",
                    job.task_id, remote_job_id, e
                ));
                continue;
            }
        };
        let secret = match hpc::credentials::resolve_secret(
            &profile.id,
            &profile.username,
            &profile.host,
            profile.credential_persisted,
        ) {
            Ok(secret) => secret,
            Err(e) => {
                errors.push(format!(
                    "Task {} (job {}) credentials failed: {}",
                    job.task_id, remote_job_id, e
                ));
                continue;
            }
        };

        let cancel_cmd = format!("scancel {}", shell_single_quote_local(remote_job_id));
        if let Err(e) = hpc::ssh::run_ssh_command(&profile, secret.as_deref(), &cancel_cmd).await {
            errors.push(format!(
                "Task {} (job {}, profile {}) cancellation failed: {}",
                job.task_id, remote_job_id, profile.id, e
            ));
        }
    }

    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "Could not safely quit because some remote jobs were not cancelled:\n{}",
            errors.join("\n")
        ))
    }
}

#[tauri::command]
async fn shutdown_and_close(app: AppHandle, state: State<'_, AppState>) -> Result<(), String> {
    if state.exit_in_progress.swap(true, Ordering::SeqCst) {
        return Ok(());
    }

    if let Err(e) = cancel_running_hpc_jobs_for_quit(state.inner()).await {
        state.exit_in_progress.store(false, Ordering::SeqCst);
        return Err(e);
    }

    state.process_manager.kill_all().await;
    // Give processes a moment to die
    tokio::time::sleep(std::time::Duration::from_millis(200)).await;
    state.allow_exit.store(true, Ordering::SeqCst);
    app.exit(0);
    Ok(())
}

fn show_main_window(app: &AppHandle) {
    if let Some(window) = app.get_webview_window("main") {
        let _ = window.show();
        let _ = window.unminimize();
        let _ = window.set_focus();
    }
}

fn request_app_quit(app: &AppHandle) {
    let state = app.state::<AppState>();
    if state.allow_exit.load(Ordering::SeqCst) || state.exit_in_progress.load(Ordering::SeqCst) {
        return;
    }

    if state.process_manager.has_running_tasks_now() {
        show_main_window(app);
        let _ = app.emit("confirm-quit", ());
    } else {
        state.allow_exit.store(true, Ordering::SeqCst);
        state.exit_in_progress.store(true, Ordering::SeqCst);
        app.exit(0);
    }
}

#[cfg(target_os = "macos")]
fn setup_macos_menu(app: &tauri::AppHandle) -> Result<(), String> {
    use tauri::menu::{AboutMetadata, Menu, MenuItem, PredefinedMenuItem, Submenu};

    let pkg_info = app.package_info();
    let config = app.config();
    let about_metadata = AboutMetadata {
        name: Some(pkg_info.name.clone()),
        version: Some(pkg_info.version.to_string()),
        copyright: config.bundle.copyright.clone(),
        authors: config.bundle.publisher.clone().map(|p| vec![p]),
        ..Default::default()
    };

    let quit_label = format!("Quit {}", pkg_info.name);
    let quit_item = MenuItem::with_id(app, "app-quit-request", quit_label, true, Some("Cmd+Q"))
        .map_err(|e| format!("Failed to create custom macOS Quit item: {}", e))?;

    let app_menu = Submenu::with_items(
        app,
        pkg_info.name.clone(),
        true,
        &[
            &PredefinedMenuItem::about(app, None, Some(about_metadata))
                .map_err(|e| format!("Failed to create About menu item: {}", e))?,
            &PredefinedMenuItem::separator(app)
                .map_err(|e| format!("Failed to create separator: {}", e))?,
            &PredefinedMenuItem::services(app, None)
                .map_err(|e| format!("Failed to create Services menu item: {}", e))?,
            &PredefinedMenuItem::separator(app)
                .map_err(|e| format!("Failed to create separator: {}", e))?,
            &PredefinedMenuItem::hide(app, None)
                .map_err(|e| format!("Failed to create Hide menu item: {}", e))?,
            &PredefinedMenuItem::hide_others(app, None)
                .map_err(|e| format!("Failed to create Hide Others menu item: {}", e))?,
            &PredefinedMenuItem::separator(app)
                .map_err(|e| format!("Failed to create separator: {}", e))?,
            &quit_item,
        ],
    )
    .map_err(|e| format!("Failed to create app menu: {}", e))?;

    let file_menu = Submenu::with_items(
        app,
        "File",
        true,
        &[&PredefinedMenuItem::close_window(app, None)
            .map_err(|e| format!("Failed to create Close Window menu item: {}", e))?],
    )
    .map_err(|e| format!("Failed to create File menu: {}", e))?;

    let edit_menu = Submenu::with_items(
        app,
        "Edit",
        true,
        &[
            &PredefinedMenuItem::undo(app, None)
                .map_err(|e| format!("Failed to create Undo menu item: {}", e))?,
            &PredefinedMenuItem::redo(app, None)
                .map_err(|e| format!("Failed to create Redo menu item: {}", e))?,
            &PredefinedMenuItem::separator(app)
                .map_err(|e| format!("Failed to create separator: {}", e))?,
            &PredefinedMenuItem::cut(app, None)
                .map_err(|e| format!("Failed to create Cut menu item: {}", e))?,
            &PredefinedMenuItem::copy(app, None)
                .map_err(|e| format!("Failed to create Copy menu item: {}", e))?,
            &PredefinedMenuItem::paste(app, None)
                .map_err(|e| format!("Failed to create Paste menu item: {}", e))?,
            &PredefinedMenuItem::select_all(app, None)
                .map_err(|e| format!("Failed to create Select All menu item: {}", e))?,
        ],
    )
    .map_err(|e| format!("Failed to create Edit menu: {}", e))?;

    let view_menu = Submenu::with_items(
        app,
        "View",
        true,
        &[&PredefinedMenuItem::fullscreen(app, None)
            .map_err(|e| format!("Failed to create Fullscreen menu item: {}", e))?],
    )
    .map_err(|e| format!("Failed to create View menu: {}", e))?;

    let window_menu = Submenu::with_items(
        app,
        "Window",
        true,
        &[
            &PredefinedMenuItem::minimize(app, None)
                .map_err(|e| format!("Failed to create Minimize menu item: {}", e))?,
            &PredefinedMenuItem::maximize(app, None)
                .map_err(|e| format!("Failed to create Maximize menu item: {}", e))?,
            &PredefinedMenuItem::separator(app)
                .map_err(|e| format!("Failed to create separator: {}", e))?,
            &PredefinedMenuItem::close_window(app, None)
                .map_err(|e| format!("Failed to create Close Window menu item: {}", e))?,
        ],
    )
    .map_err(|e| format!("Failed to create Window menu: {}", e))?;

    let menu = Menu::with_items(
        app,
        &[&app_menu, &file_menu, &edit_menu, &view_menu, &window_menu],
    )
    .map_err(|e| format!("Failed to create macOS app menu: {}", e))?;

    app.set_menu(menu)
        .map_err(|e| format!("Failed to install macOS app menu: {}", e))?;

    Ok(())
}

#[cfg(target_os = "windows")]
fn setup_windows_tray(app: &tauri::AppHandle) -> Result<(), String> {
    use tauri::menu::{Menu, MenuItem};
    use tauri::tray::TrayIconBuilder;

    let open_item = MenuItem::with_id(app, "tray-open", "Open QCortado", true, None::<&str>)
        .map_err(|e| format!("Failed to create tray Open item: {}", e))?;
    let quit_item = MenuItem::with_id(app, "tray-quit", "Quit QCortado", true, None::<&str>)
        .map_err(|e| format!("Failed to create tray Quit item: {}", e))?;
    let menu = Menu::with_items(app, &[&open_item, &quit_item])
        .map_err(|e| format!("Failed to create tray menu: {}", e))?;

    let mut tray_builder = TrayIconBuilder::with_id("qcortado-main")
        .tooltip("QCortado")
        .menu(&menu)
        .show_menu_on_left_click(true)
        .on_menu_event(|app, event| match event.id().as_ref() {
            "tray-open" => {
                show_main_window(app);
            }
            "tray-quit" => {
                request_app_quit(app);
            }
            _ => {}
        });

    if let Some(icon) = app.default_window_icon().cloned() {
        tray_builder = tray_builder.icon(icon);
    }

    tray_builder
        .build(app)
        .map_err(|e| format!("Failed to initialize Windows tray icon: {}", e))?;

    Ok(())
}

// ============================================================================
// Application Entry Point
// ============================================================================

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    let builder = tauri::Builder::default()
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
            let mut fermi_surfer_path: Option<PathBuf> = None;
            let mut wannier90_path: Option<PathBuf> = None;
            let mut postw90_path: Option<PathBuf> = None;
            let mut execution_prefix: Option<String> = None;
            let mut mpi_defaults: Option<config::MpiDefaultsConfig> = None;
            let mut qe_defaults = config::QeDefaultsConfig::default();
            let mut save_size_mode = config::SaveSizeMode::Large;
            let mut execution_mode = hpc::profile::ExecutionMode::Local;
            let mut hpc_profiles: Vec<hpc::profile::HpcProfile> = Vec::new();
            let mut active_hpc_profile_id: Option<String> = None;
            let mut viewer_auto_publish_enabled = true;
            let mut viewer_sync_status = hpc::viewer_library::ViewerSyncStatus::default();
            match config::load_config(&app.handle()) {
                Ok(cfg) => {
                    if let Some(path) = cfg.qe_bin_dir {
                        let path_buf = PathBuf::from(&path);
                        // Only use if pw.x still exists
                        if path_buf.join("pw.x").exists() {
                            qe_bin_dir = Some(path_buf);
                        }
                    }
                    if let Some(path) = cfg.fermi_surfer_path {
                        let path_buf = PathBuf::from(&path);
                        if path_buf.exists() && path_buf.is_file() {
                            fermi_surfer_path = Some(path_buf);
                        }
                    }
                    if let Some(path) = cfg.wannier90_path {
                        let path_buf = PathBuf::from(&path);
                        if path_buf.exists() && path_buf.is_file() {
                            wannier90_path = Some(path_buf);
                        }
                    }
                    if let Some(path) = cfg.postw90_path {
                        let path_buf = PathBuf::from(&path);
                        if path_buf.exists() && path_buf.is_file() {
                            postw90_path = Some(path_buf);
                        }
                    }
                    execution_prefix = normalize_execution_prefix(cfg.execution_prefix);
                    mpi_defaults = normalize_mpi_defaults(cfg.mpi_defaults);
                    qe_defaults = cfg.qe_defaults;
                    save_size_mode = cfg.save_size_mode;
                    execution_mode = cfg.execution_mode;
                    hpc_profiles = cfg.hpc_profiles;
                    active_hpc_profile_id = cfg.active_hpc_profile_id;
                    viewer_auto_publish_enabled = cfg.viewer_auto_publish_enabled;
                    viewer_sync_status.last_synced_at = cfg.viewer_last_sync_at;
                    viewer_sync_status.last_error = cfg.viewer_last_sync_error;
                }
                Err(e) => {
                    eprintln!("Warning: Failed to load config: {}", e);
                }
            }

            // Initialize AppState with loaded config
            app.manage(AppState {
                qe_bin_dir: Mutex::new(qe_bin_dir),
                fermi_surfer_path: Mutex::new(fermi_surfer_path),
                wannier90_path: Mutex::new(wannier90_path),
                postw90_path: Mutex::new(postw90_path),
                execution_prefix: Mutex::new(execution_prefix),
                mpi_defaults: Mutex::new(mpi_defaults),
                qe_defaults: Mutex::new(qe_defaults),
                save_size_mode: Mutex::new(save_size_mode),
                execution_mode: Mutex::new(execution_mode),
                hpc_profiles: Mutex::new(hpc_profiles),
                active_hpc_profile_id: Mutex::new(active_hpc_profile_id),
                project_dir: Mutex::new(None),
                process_manager: ProcessManager::new(),
                allow_exit: AtomicBool::new(false),
                exit_in_progress: AtomicBool::new(false),
                viewer_auto_publish_enabled: AtomicBool::new(viewer_auto_publish_enabled),
                viewer_publish_pending: AtomicBool::new(false),
                viewer_sync_status: Mutex::new(viewer_sync_status),
            });

            #[cfg(feature = "viewer")]
            {
                let app_handle = app.handle().clone();
                tauri::async_runtime::spawn(async move {
                    loop {
                        let state = app_handle.state::<AppState>();
                        let _ = sync_viewer_library_with_profile(&app_handle, state.inner(), None)
                            .await;
                        tokio::time::sleep(Duration::from_secs(120)).await;
                    }
                });
            }

            #[cfg(target_os = "windows")]
            {
                setup_windows_tray(&app.handle())?;
            }

            #[cfg(target_os = "macos")]
            {
                setup_macos_menu(&app.handle())?;
            }

            Ok(())
        })
        .on_window_event(|window, event| {
            // Keep the app process alive when users click X on the main window.
            if window.label() != "main" {
                return;
            }
            if let tauri::WindowEvent::CloseRequested { api, .. } = event {
                let app = window.app_handle();
                let state = app.state::<AppState>();
                if state.allow_exit.load(Ordering::SeqCst) {
                    return;
                }
                api.prevent_close();
                let _ = window.hide();
            }
        })
        .on_menu_event(|app, event| match event.id().as_ref() {
            #[cfg(target_os = "macos")]
            "app-quit-request" => {
                request_app_quit(app);
            }
            "quit" => {
                request_app_quit(app);
            }
            _ => {}
        });

    #[cfg(feature = "viewer")]
    let builder = builder.invoke_handler(tauri::generate_handler![
        get_app_role,
        hpc_list_profiles,
        hpc_import_preset_bundle,
        hpc_save_profile,
        hpc_set_active_profile,
        hpc_get_active_profile_id,
        hpc_test_connection,
        hpc_validate_environment,
        viewer_sync_remote_library,
        viewer_get_sync_status,
        list_running_tasks,
        get_task_info,
        get_task_output,
        dismiss_task,
        shutdown_and_close,
        projects::list_projects,
        projects::list_project_folders,
        projects::list_multiview_band_calculations,
        projects::get_project,
        projects::get_project_calculation,
        projects::get_project_calculation_logs,
        projects::get_project_calculation_inputs,
        projects::set_last_opened_cif,
        projects::get_cif_crystal_data,
        projects::get_cif_content,
        projects::get_saved_phonon_data,
    ]);

    #[cfg(not(feature = "viewer"))]
    let builder = builder.invoke_handler(tauri::generate_handler![
        get_app_role,
        get_viewer_auto_publish_enabled,
        set_viewer_auto_publish_enabled,
        set_qe_path,
        get_qe_path,
        set_fermi_surfer_path,
        get_fermi_surfer_path,
        set_wannier90_path,
        get_wannier90_path,
        set_postw90_path,
        get_postw90_path,
        set_execution_prefix,
        get_execution_prefix,
        set_mpi_defaults,
        get_mpi_defaults,
        set_qe_defaults,
        get_qe_defaults,
        set_save_size_mode,
        get_save_size_mode,
        set_execution_mode,
        get_execution_mode,
        hpc_list_profiles,
        hpc_export_preset_bundle,
        hpc_import_preset_bundle,
        hpc_save_profile,
        hpc_update_profile_defaults,
        hpc_migrate_remote_roots,
        hpc_delete_profile,
        hpc_set_active_profile,
        hpc_get_active_profile_id,
        hpc_test_connection,
        hpc_validate_environment,
        hpc_get_cluster_snapshot,
        hpc_sample_utilization,
        hpc_list_remote_pseudopotentials,
        hpc_list_remote_pseudopotential_metadata,
        hpc_load_remote_sssp_data,
        hpc_preview_slurm_script,
        hpc_list_headless_jobs,
        hpc_attach_headless_job,
        hpc_open_activity_window,
        hpc_download_task_artifacts,
        hpc_download_calculation_artifacts,
        hpc_list_recoverable_remote_phonon_runs,
        hpc_debug_remote_phonon_recovery,
        hpc_recover_remote_phonon_calculation,
        hpc_clean_remote_orphans,
        hpc_publish_viewer_library,
        clear_temp_storage,
        launch_fermi_surface_viewer,
        export_wannier_for_ludwig,
        check_qe_executables,
        analyze_structure_symmetry,
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
        list_pseudopotential_metadata,
        load_sssp_data,
        validate_epw_prerequisites,
        build_epw_input_preview,
        // Background task commands
        start_scf_calculation,
        start_bands_calculation,
        start_dos_calculation,
        start_fermi_surface_calculation,
        start_phonon_calculation,
        start_epw_calculation,
        start_epw_calculation_hpc,
        start_wannier_calculation,
        start_transport_calculation,
        list_running_tasks,
        get_task_info,
        get_task_output,
        cancel_task,
        dismiss_task,
        shutdown_and_close,
        // Project management commands
        projects::list_projects,
        projects::list_project_folders,
        projects::list_multiview_band_calculations,
        projects::scan_storage_inventory,
        projects::delete_storage_entries,
        projects::delete_storage_calculations,
        projects::delete_storage_selection,
        projects::lighten_storage_calculations,
        projects::create_project,
        projects::create_project_folder,
        projects::get_project,
        projects::get_project_calculation,
        projects::get_project_calculation_logs,
        projects::get_project_calculation_inputs,
        projects::update_project_metadata,
        projects::update_calculation_name,
        projects::update_calculation_band_viewer_metadata,
        projects::rename_project_folder,
        projects::delete_project_folder,
        projects::move_project_to_folder,
        projects::add_cif_to_project,
        projects::save_calculation,
        projects::export_project_archive,
        projects::cancel_project_export,
        projects::import_project_archive,
        projects::delete_project,
        projects::delete_calculation,
        projects::set_calculation_tag,
        projects::set_last_opened_cif,
        projects::get_cif_crystal_data,
        projects::get_cif_content,
        projects::get_saved_phonon_data,
        projects::recover_phonon_calculation,
        viewer_sync_remote_library,
        viewer_get_sync_status,
    ]);

    let app = builder
        .build(tauri::generate_context!())
        .expect("error while building tauri application");

    app.run(|app_handle, event| match event {
        tauri::RunEvent::ExitRequested { api, .. } => {
            let state = app_handle.state::<AppState>();
            if state.allow_exit.load(Ordering::SeqCst) {
                return;
            }

            api.prevent_exit();
            request_app_quit(app_handle);
        }
        tauri::RunEvent::Resumed => {
            let _ = app_handle.emit("app-resumed", ());
        }
        #[cfg(target_os = "macos")]
        tauri::RunEvent::Reopen {
            has_visible_windows,
            ..
        } => {
            if !has_visible_windows {
                show_main_window(app_handle);
            }
        }
        _ => {}
    });
}
