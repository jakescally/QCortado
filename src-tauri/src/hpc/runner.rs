use std::collections::{HashMap, HashSet};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::time::{Duration, Instant};

use tauri::{AppHandle, Emitter};

use crate::process_manager::ProcessManager;

use super::profile::HpcProfile;
use super::slurm::{
    is_terminal_state, normalize_scheduler_state, parse_sacct_snapshot, parse_sbatch_job_id,
    parse_squeue_snapshot, SchedulerSnapshot,
};
use super::ssh::{download_file, run_ssh_command, upload_directory};
const LIVE_LINECOUNT_MARKER: &str = "__QCORTADO_INTERNAL_LINECOUNT__=";
const LIVE_SIZE_MARKER: &str = "__QCORTADO_INTERNAL_SIZE__=";
const ENABLE_EXPERIMENTAL_REMOTE_LIVE_LOGGING: bool = false;

#[derive(Debug, Clone)]
pub struct HpcBatchRequest {
    pub task_id: String,
    pub task_kind: String,
    pub task_label: String,
    pub profile: HpcProfile,
    pub secret: Option<String>,
    pub slurm_script: String,
    pub sbatch_preview: String,
    pub bundle_dir: PathBuf,
    pub local_sync_dir: PathBuf,
    pub cancel_flag: Arc<std::sync::atomic::AtomicBool>,
}

#[derive(Debug, Clone)]
pub struct HpcAttachRequest {
    pub task_id: String,
    pub task_kind: String,
    pub task_label: String,
    pub profile: HpcProfile,
    pub secret: Option<String>,
    pub remote_job_id: String,
    pub remote_workdir: String,
    pub remote_project_path: Option<String>,
    pub slurm_script: Option<String>,
    pub sbatch_preview: Option<String>,
    pub local_sync_dir: PathBuf,
    pub cancel_flag: Arc<std::sync::atomic::AtomicBool>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct HpcBatchResult {
    pub backend: String,
    pub task_kind: String,
    pub task_label: String,
    pub remote_job_id: String,
    pub remote_workdir: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_project_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_node: Option<String>,
    pub terminal_state: String,
    pub sbatch_preview: String,
    pub slurm_script: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArtifactSyncMode {
    Minimal,
    EpwReady,
    Full,
}

impl ArtifactSyncMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Minimal => "minimal",
            Self::EpwReady => "epw-ready",
            Self::Full => "full",
        }
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, Default)]
pub struct HpcArtifactSyncReport {
    pub mode: String,
    pub downloaded_files: usize,
    pub downloaded_bytes: u64,
    pub skipped_files: usize,
    pub skipped_bytes: u64,
}

#[derive(Debug, Clone)]
pub struct HpcArtifactSyncRequest {
    pub task_id: String,
    pub task_kind: String,
    pub profile: HpcProfile,
    pub secret: Option<String>,
    pub remote_workdir: String,
    pub local_sync_dir: PathBuf,
    pub mode: ArtifactSyncMode,
}

#[derive(Debug, Clone)]
struct RemoteFileEntry {
    rel_path: String,
    size_bytes: u64,
}

#[derive(Debug, Clone, Default)]
struct LiveFileCursor {
    byte_offset: u64,
    line_count: usize,
    pending_line: String,
    warned_read_error: bool,
    warned_parse_error: bool,
}

fn shell_single_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    let escaped = value.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

fn emit_task_event_line(app: &AppHandle, task_id: &str, line: &str) {
    let _ = app.emit(&format!("task-output:{}", task_id), line);
}

async fn emit_task_line(app: &AppHandle, pm: &ProcessManager, task_id: &str, line: String) {
    emit_task_event_line(app, task_id, &line);
    pm.append_output(task_id, line).await;
}

fn push_unique_remote_path(targets: &mut Vec<String>, candidate: Option<&str>) {
    let Some(path) = candidate.map(str::trim).filter(|value| !value.is_empty()) else {
        return;
    };
    if targets.iter().any(|existing| existing == path) {
        return;
    }
    targets.push(path.to_string());
}

async fn remove_remote_paths_best_effort(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    targets: &[String],
    failure_label: &str,
) {
    for target in targets {
        let remove_cmd = format!(
            "if [ -e {target} ]; then rm -rf -- {target}; fi",
            target = shell_single_quote(target)
        );
        if let Err(err) = run_ssh_command(profile, secret, &remove_cmd).await {
            emit_task_line(
                app,
                pm,
                task_id,
                format!("HPC_WARNING|Failed to remove remote {failure_label} at {target}: {err}"),
            )
            .await;
        }
    }
}

async fn archive_remote_run_copy(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    remote_project_path: &str,
) -> Result<(), String> {
    emit_task_line(
        app,
        pm,
        task_id,
        format!("HPC_STAGE|Archiving|{}", remote_project_path),
    )
    .await;
    emit_task_line(
        app,
        pm,
        task_id,
        format!(
            "[QCortado] Remote archive copy started at {}.",
            chrono::Utc::now().to_rfc3339()
        ),
    )
    .await;
    let archive_started = Instant::now();

    let archive_cmd = format!(
        "mkdir -p {dest} && if [ -d {src} ]; then cp -a {src}/. {dest}/; fi",
        src = shell_single_quote(remote_workdir),
        dest = shell_single_quote(remote_project_path)
    );
    let archive_result = run_ssh_command(profile, secret, &archive_cmd).await;
    match archive_result {
        Ok(_) => {
            emit_task_line(
                app,
                pm,
                task_id,
                format!("HPC_STAGE|Archived|{}", remote_project_path),
            )
            .await;
            emit_task_line(
                app,
                pm,
                task_id,
                format!(
                    "[QCortado] Remote archive copy finished in {:.1}s.",
                    archive_started.elapsed().as_secs_f64()
                ),
            )
            .await;
            Ok(())
        }
        Err(err) => {
            emit_task_line(
                app,
                pm,
                task_id,
                format!(
                    "HPC_WARNING|Failed to archive run under remote project root ({}): {}",
                    remote_project_path, err
                ),
            )
            .await;
            emit_task_line(
                app,
                pm,
                task_id,
                format!(
                    "[QCortado] Remote archive copy failed after {:.1}s.",
                    archive_started.elapsed().as_secs_f64()
                ),
            )
            .await;
            Err(err)
        }
    }
}

fn normalize_rel_path(raw: &str) -> Option<String> {
    let trimmed = raw.trim().trim_start_matches("./");
    if trimmed.is_empty() {
        return None;
    }
    let path = Path::new(trimmed);
    for component in path.components() {
        match component {
            Component::Normal(_) => {}
            _ => return None,
        }
    }
    Some(trimmed.to_string())
}

fn parse_remote_manifest(output: &str) -> Vec<RemoteFileEntry> {
    output
        .lines()
        .filter_map(|line| {
            let (raw_rel, raw_size) = line.split_once('\t')?;
            let rel_path = normalize_rel_path(raw_rel)?;
            let size_bytes = raw_size.trim().parse::<u64>().unwrap_or(0);
            Some(RemoteFileEntry {
                rel_path,
                size_bytes,
            })
        })
        .collect()
}

fn is_heavy_scratch_path(rel_path: &str) -> bool {
    let lower_path = rel_path.to_ascii_lowercase();
    let path = Path::new(rel_path);
    for component in path.components() {
        let Component::Normal(segment) = component else {
            continue;
        };
        let lower = segment.to_string_lossy().to_ascii_lowercase();
        if lower == "tmp"
            || lower.starts_with("_ph")
            || lower.ends_with(".save")
            || lower == "save"
            || lower == "wfc"
        {
            return true;
        }
    }
    lower_path.ends_with(".wfc")
        || lower_path.ends_with(".wfc1")
        || lower_path.ends_with(".igk")
        || lower_path.ends_with(".hub")
        || lower_path.ends_with(".mix")
        || lower_path.ends_with(".rho")
}

fn should_download_minimal(task_kind: &str, entry: &RemoteFileEntry) -> bool {
    let lower_rel = entry.rel_path.to_ascii_lowercase();
    if lower_rel.ends_with(".frmsf") || lower_rel.ends_with(".bxsf") {
        return true;
    }

    if is_heavy_scratch_path(&entry.rel_path) {
        return false;
    }

    let file_name = Path::new(&lower_rel)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("");
    let is_top_level = !lower_rel.contains('/');

    if matches!(
        file_name,
        "run.sbatch"
            | "slurm.out"
            | "slurm.err"
            | "pw.in"
            | "pw.out"
            | "bands.in"
            | "bands.out"
            | "bands_pp.in"
            | "bands_pp.out"
            | "projwfc.in"
            | "projwfc.out"
            | "nscf.in"
            | "nscf.out"
            | "pw2wan.in"
            | "pw2wan.out"
            | "wannier90_pre.out"
            | "wannier90.out"
            | "dos.in"
            | "dos.out"
            | "fermi_velocity.in"
            | "fermi_velocity.out"
            | "ph.in"
            | "ph.out"
            | "q2r.in"
            | "q2r.out"
            | "matdyn_dos.in"
            | "matdyn_dos.out"
            | "matdyn_bands.in"
            | "matdyn_bands.out"
            | "phonon_dos"
            | "phonon_freq"
            | "phonon_freq.gp"
            | "force_constants"
    ) {
        return true;
    }

    if matches!(
        file_name,
        "run.sbatch"
            | "slurm.out"
            | "slurm.err"
            | "nscf.in"
            | "nscf.out"
            | "pw2wan.in"
            | "pw2wan.out"
            | "wannier90_pre.out"
            | "wannier90.out"
    ) && task_kind == "wannier"
    {
        return true;
    }

    if task_kind == "transport"
        && (file_name.ends_with(".win")
            || file_name.ends_with(".wpout")
            || file_name.ends_with(".werr")
            || file_name.ends_with(".chk")
            || file_name.ends_with(".eig")
            || file_name.ends_with(".nnkp")
            || file_name.ends_with("_elcond.dat")
            || file_name.ends_with("_sigmas.dat")
            || file_name.ends_with("_seebeck.dat")
            || file_name.ends_with("_kappa.dat")
            || file_name.ends_with("_tdf.dat"))
    {
        return entry.size_bytes <= 256 * 1024 * 1024;
    }

    if lower_rel.ends_with(".gnu")
        || lower_rel.ends_with(".json")
        || lower_rel.ends_with(".txt")
        || lower_rel.ends_with(".log")
    {
        return true;
    }

    if task_kind == "wannier"
        && (file_name.ends_with(".win")
            || file_name.ends_with(".nnkp")
            || file_name.ends_with(".amn")
            || file_name.ends_with(".mmn")
            || file_name.ends_with(".eig")
            || file_name.ends_with("_wsvec.dat")
            || file_name.ends_with(".wout")
            || file_name.ends_with(".chk")
            || file_name.ends_with("_hr.dat")
            || file_name.ends_with("_centres.xyz")
            || file_name.ends_with("_band.dat")
            || file_name.ends_with("_band.kpt"))
    {
        return entry.size_bytes <= 256 * 1024 * 1024;
    }

    if lower_rel.ends_with(".in") || lower_rel.ends_with(".out") || lower_rel.ends_with(".err") {
        return true;
    }

    if lower_rel.ends_with(".dat") {
        return entry.size_bytes <= 256 * 1024 * 1024;
    }

    if task_kind == "phonon" && (file_name.starts_with("dyn") || file_name.starts_with("matdyn")) {
        return entry.size_bytes <= 256 * 1024 * 1024;
    }

    if is_top_level {
        return entry.size_bytes <= 64 * 1024 * 1024;
    }

    false
}

fn should_download_epw_ready(task_kind: &str, entry: &RemoteFileEntry) -> bool {
    if should_download_minimal(task_kind, entry) {
        return true;
    }
    if task_kind != "epw" {
        return false;
    }
    if is_heavy_scratch_path(&entry.rel_path) {
        return false;
    }

    let lower_rel = entry.rel_path.to_ascii_lowercase();
    let file_name = Path::new(&lower_rel)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("");

    if matches!(
        file_name,
        "epw.in" | "epw.out" | "epw.err" | "run.sbatch" | "slurm.out" | "slurm.err"
    ) {
        return true;
    }

    if file_name.starts_with("epw")
        || file_name.ends_with(".epw")
        || file_name.ends_with(".epb")
        || file_name.ends_with(".ephmat")
        || file_name.ends_with(".a2f")
        || file_name.ends_with(".freq")
        || file_name.ends_with(".frq")
        || file_name.ends_with(".fmt")
        || file_name.ends_with(".xml")
        || file_name.ends_with(".dat")
    {
        return entry.size_bytes <= 1024 * 1024 * 1024;
    }

    false
}

async fn list_remote_files(
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
) -> Result<Vec<RemoteFileEntry>, String> {
    let manifest_cmd = format!(
        "cd {} && find . -type f -printf '%P\\t%s\\n'",
        shell_single_quote(remote_workdir)
    );
    let manifest_output = run_ssh_command(profile, secret, &manifest_cmd).await?;
    Ok(parse_remote_manifest(&manifest_output))
}

pub async fn sync_remote_artifacts(
    app: &AppHandle,
    pm: &ProcessManager,
    request: HpcArtifactSyncRequest,
) -> Result<HpcArtifactSyncReport, String> {
    std::fs::create_dir_all(&request.local_sync_dir).map_err(|e| {
        format!(
            "Failed to create local sync directory {}: {}",
            request.local_sync_dir.display(),
            e
        )
    })?;

    let mut report = HpcArtifactSyncReport {
        mode: request.mode.label().to_string(),
        ..HpcArtifactSyncReport::default()
    };
    emit_task_line(
        app,
        pm,
        &request.task_id,
        format!(
            "[QCortado] Starting {} artifact sync from {} to {}",
            request.mode.label(),
            request.remote_workdir,
            request.local_sync_dir.display()
        ),
    )
    .await;

    match request.mode {
        ArtifactSyncMode::Full => {
            emit_task_line(
                app,
                pm,
                &request.task_id,
                "[QCortado] Enumerating remote files for full sync...".to_string(),
            )
            .await;
            let scan_started = Instant::now();
            let remote_files = list_remote_files(
                &request.profile,
                request.secret.as_deref(),
                &request.remote_workdir,
            )
            .await?;
            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!(
                    "[QCortado] Full sync manifest complete: {} files discovered in {:.1}s.",
                    remote_files.len(),
                    scan_started.elapsed().as_secs_f64()
                ),
            )
            .await;
            let total_files = remote_files.len();
            let total_bytes = remote_files
                .iter()
                .fold(0u64, |sum, entry| sum.saturating_add(entry.size_bytes));

            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!("HPC_TRANSFER|start|full|{}|{}", total_files, total_bytes),
            )
            .await;

            let mut processed_files: usize = 0;
            let mut failures: Vec<String> = Vec::new();
            let mut last_emit = Instant::now();

            for entry in remote_files {
                let remote_file = format!(
                    "{}/{}",
                    request.remote_workdir.trim_end_matches('/'),
                    entry.rel_path
                );
                let local_file = request.local_sync_dir.join(&entry.rel_path);
                match download_file(
                    &request.profile,
                    request.secret.as_deref(),
                    &remote_file,
                    &local_file,
                )
                .await
                {
                    Ok(_) => {
                        report.downloaded_files += 1;
                        report.downloaded_bytes =
                            report.downloaded_bytes.saturating_add(entry.size_bytes);
                    }
                    Err(err) => {
                        report.skipped_files += 1;
                        report.skipped_bytes =
                            report.skipped_bytes.saturating_add(entry.size_bytes);
                        failures.push(format!("{} ({})", entry.rel_path, err));
                    }
                }
                processed_files += 1;

                let should_emit = processed_files == total_files
                    || processed_files % 16 == 0
                    || last_emit.elapsed() >= Duration::from_millis(900);
                if should_emit {
                    emit_task_line(
                        app,
                        pm,
                        &request.task_id,
                        format!(
                            "HPC_TRANSFER|progress|full|{}|{}|{}|{}|{}",
                            processed_files,
                            total_files,
                            report.downloaded_bytes,
                            total_bytes,
                            report.skipped_files
                        ),
                    )
                    .await;
                    last_emit = Instant::now();
                }
            }

            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!(
                    "HPC_TRANSFER|done|full|{}|{}|{}|{}|{}",
                    report.downloaded_files,
                    total_files,
                    report.downloaded_bytes,
                    total_bytes,
                    report.skipped_files
                ),
            )
            .await;

            if !failures.is_empty() {
                let preview = failures
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<String>>()
                    .join(", ");
                emit_task_line(
                    app,
                    pm,
                    &request.task_id,
                    format!(
                        "HPC_WARNING|Some files failed during full download ({}).",
                        preview
                    ),
                )
                .await;
            }
            Ok(report)
        }
        ArtifactSyncMode::Minimal | ArtifactSyncMode::EpwReady => {
            let mode_label = request.mode.label();
            let mode_title = if matches!(request.mode, ArtifactSyncMode::EpwReady) {
                "EPW-ready"
            } else {
                "Minimal"
            };
            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!(
                    "[QCortado] Enumerating remote files for {} sync...",
                    mode_label
                ),
            )
            .await;
            let scan_started = Instant::now();
            let remote_files = list_remote_files(
                &request.profile,
                request.secret.as_deref(),
                &request.remote_workdir,
            )
            .await?;
            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!(
                    "[QCortado] {} sync manifest complete: {} files discovered in {:.1}s.",
                    mode_title,
                    remote_files.len(),
                    scan_started.elapsed().as_secs_f64()
                ),
            )
            .await;
            let mut candidates: Vec<RemoteFileEntry> = Vec::new();
            for entry in remote_files {
                let include = if matches!(request.mode, ArtifactSyncMode::EpwReady) {
                    should_download_epw_ready(&request.task_kind, &entry)
                } else {
                    should_download_minimal(&request.task_kind, &entry)
                };
                if include {
                    candidates.push(entry);
                } else {
                    report.skipped_files += 1;
                    report.skipped_bytes = report.skipped_bytes.saturating_add(entry.size_bytes);
                }
            }
            let total_candidates = candidates.len();
            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!(
                    "[QCortado] {} sync selected {} files ({} skipped as heavy/scratch).",
                    mode_title, total_candidates, report.skipped_files
                ),
            )
            .await;

            let mut failures: Vec<String> = Vec::new();
            let mut processed_files: usize = 0;
            let mut last_emit = Instant::now();
            for entry in candidates {
                let remote_file = format!(
                    "{}/{}",
                    request.remote_workdir.trim_end_matches('/'),
                    entry.rel_path
                );
                let local_file = request.local_sync_dir.join(&entry.rel_path);
                match download_file(
                    &request.profile,
                    request.secret.as_deref(),
                    &remote_file,
                    &local_file,
                )
                .await
                {
                    Ok(_) => {
                        report.downloaded_files += 1;
                        report.downloaded_bytes =
                            report.downloaded_bytes.saturating_add(entry.size_bytes);
                    }
                    Err(err) => {
                        report.skipped_files += 1;
                        report.skipped_bytes =
                            report.skipped_bytes.saturating_add(entry.size_bytes);
                        failures.push(format!("{} ({})", entry.rel_path, err));
                    }
                }
                processed_files += 1;
                let should_emit = processed_files == total_candidates
                    || processed_files % 16 == 0
                    || last_emit.elapsed() >= Duration::from_millis(900);
                if should_emit {
                    emit_task_line(
                        app,
                        pm,
                        &request.task_id,
                        format!(
                            "[QCortado] {} sync progress: {}/{} files, {:.2} MB downloaded, {} skipped.",
                            mode_title,
                            processed_files,
                            total_candidates,
                            report.downloaded_bytes as f64 / (1024.0 * 1024.0),
                            report.skipped_files
                        ),
                    )
                    .await;
                    last_emit = Instant::now();
                }
            }

            if !failures.is_empty() {
                let preview = failures
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<String>>()
                    .join(", ");
                emit_task_line(
                    app,
                    pm,
                    &request.task_id,
                    format!(
                        "HPC_WARNING|Some artifacts were not downloaded during {} sync ({}).",
                        mode_label, preview
                    ),
                )
                .await;
            }

            emit_task_line(
                app,
                pm,
                &request.task_id,
                format!(
                    "[QCortado] {} sync finished: {} files downloaded ({:.2} MB), {} skipped ({:.2} MB).",
                    mode_title,
                    report.downloaded_files,
                    report.downloaded_bytes as f64 / (1024.0 * 1024.0),
                    report.skipped_files,
                    report.skipped_bytes as f64 / (1024.0 * 1024.0)
                ),
            )
            .await;

            Ok(report)
        }
    }
}

async fn update_scheduler_snapshot(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    snapshot: &SchedulerSnapshot,
) {
    pm.set_scheduler_state(task_id, Some(snapshot.state.clone()))
        .await;
    if let Some(node) = snapshot.node.as_ref() {
        pm.set_remote_node(task_id, Some(node.clone())).await;
    }
    let state_line = format!(
        "HPC_SCHED|{}|{}",
        snapshot.state,
        snapshot.node.as_deref().unwrap_or("")
    );
    emit_task_line(app, pm, task_id, state_line).await;
}

async fn emit_remote_tail_if_exists(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    file_name: &str,
    max_lines: usize,
) {
    let cmd = format!(
        "cd {} && if [ -f {} ]; then tail -n {} {}; fi",
        shell_single_quote(remote_workdir),
        shell_single_quote(file_name),
        max_lines.max(1),
        shell_single_quote(file_name)
    );

    let Ok(content) = run_ssh_command(profile, secret, &cmd).await else {
        return;
    };
    if content.trim().is_empty() {
        return;
    }

    emit_task_line(
        app,
        pm,
        task_id,
        format!("--- Remote {} (tail) ---", file_name),
    )
    .await;
    for line in content.lines() {
        emit_task_line(app, pm, task_id, line.to_string()).await;
    }
}

fn parse_redirect_target(after_redirect: &str) -> Option<String> {
    let trimmed = after_redirect.trim_start();
    if trimmed.is_empty() {
        return None;
    }

    let mut chars = trimmed.chars();
    let first = chars.next()?;
    if first == '"' || first == '\'' {
        let mut value = String::new();
        for ch in chars {
            if ch == first {
                break;
            }
            value.push(ch);
        }
        if value.is_empty() {
            None
        } else {
            Some(value)
        }
    } else {
        let token = trimmed
            .split_whitespace()
            .next()
            .unwrap_or("")
            .trim_end_matches(';')
            .to_string();
        if token.is_empty() {
            None
        } else {
            Some(token)
        }
    }
}

fn collect_redirect_targets(line: &str) -> Vec<String> {
    let mut targets = Vec::new();
    for (idx, ch) in line.char_indices() {
        if ch != '>' {
            continue;
        }
        if let Some(target) = parse_redirect_target(&line[idx + ch.len_utf8()..]) {
            targets.push(target);
        }
    }
    targets
}

fn shell_words_for_live_scan(line: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut current = String::new();
    let mut quote: Option<char> = None;
    let mut escaped = false;

    for ch in line.chars() {
        if escaped {
            current.push(ch);
            escaped = false;
            continue;
        }
        if quote != Some('\'') && ch == '\\' {
            escaped = true;
            continue;
        }
        if let Some(quote_char) = quote {
            if ch == quote_char {
                quote = None;
            } else {
                current.push(ch);
            }
            continue;
        }
        if ch == '\'' || ch == '"' {
            quote = Some(ch);
            continue;
        }
        if ch == ';' {
            if !current.is_empty() {
                words.push(std::mem::take(&mut current));
            }
            words.push(";".to_string());
            continue;
        }
        if ch.is_whitespace() {
            if !current.is_empty() {
                words.push(std::mem::take(&mut current));
            }
            continue;
        }
        current.push(ch);
    }

    if !current.is_empty() {
        words.push(current);
    }

    words
}

fn collect_option_targets(line: &str, option: &str) -> Vec<String> {
    let words = shell_words_for_live_scan(line);
    let mut targets = Vec::new();
    for pair in words.windows(2) {
        if pair[0] != option {
            continue;
        }
        let target = pair[1].trim();
        if target.is_empty() || target == ";" || target.starts_with('-') {
            continue;
        }
        targets.push(target.to_string());
    }
    targets
}

fn rank_zero_output_pattern(pattern: &str) -> Option<String> {
    let trimmed = pattern.trim();
    if trimmed.is_empty() {
        return None;
    }
    if trimmed.contains("%r") {
        return Some(trimmed.replace("%r", "0"));
    }
    if trimmed.contains("%g") {
        return Some(trimmed.replace("%g", "0"));
    }
    Some(trimmed.to_string())
}

fn is_live_output_file_name(file_name: &str) -> bool {
    file_name.ends_with(".out")
        || file_name.ends_with(".err")
        || file_name.ends_with(".wpout")
        || file_name.ends_with(".werr")
}

fn collect_live_output_files(slurm_script: &str) -> Vec<String> {
    let mut files = vec!["slurm.out".to_string(), "slurm.err".to_string()];
    let mut seen: HashSet<String> = files.iter().cloned().collect();

    for line in slurm_script.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        for target in collect_redirect_targets(trimmed) {
            if target == "&1" || target.starts_with('&') {
                continue;
            }
            if !is_live_output_file_name(&target) {
                continue;
            }
            if seen.insert(target.clone()) {
                files.push(target);
            }
        }
        if ENABLE_EXPERIMENTAL_REMOTE_LIVE_LOGGING {
            for pattern in ["-outfile-pattern", "-errfile-pattern"]
                .into_iter()
                .flat_map(|option| collect_option_targets(trimmed, option))
            {
                let Some(target) = rank_zero_output_pattern(&pattern) else {
                    continue;
                };
                if seen.insert(target.clone()) {
                    files.push(target);
                }
            }
        }
    }

    files
}

fn build_incremental_read_command(
    remote_workdir: &str,
    file_name: &str,
    byte_offset: u64,
) -> String {
    format!(
        "cd {} && if [ -f {} ]; then __qcortado_size=$(wc -c < {} | tr -d '[:space:]'); printf '{}%s\\n' \"$__qcortado_size\"; if [ \"$__qcortado_size\" -gt {} ]; then tail -c +$(({} + 1)) {}; fi; fi",
        shell_single_quote(remote_workdir),
        shell_single_quote(file_name),
        shell_single_quote(file_name),
        LIVE_SIZE_MARKER,
        byte_offset,
        byte_offset,
        shell_single_quote(file_name)
    )
}

fn build_linecount_read_command(
    remote_workdir: &str,
    file_name: &str,
    last_line_count: usize,
) -> String {
    format!(
        "cd {} && if [ -f {} ]; then awk 'NR > {} {{print}} END {{print \"{}\" NR}}' {}; fi",
        shell_single_quote(remote_workdir),
        shell_single_quote(file_name),
        last_line_count,
        LIVE_LINECOUNT_MARKER,
        shell_single_quote(file_name)
    )
}

fn build_live_probe_command(remote_workdir: &str, file_names: &[String]) -> String {
    let mut set_args = "set --".to_string();
    for file_name in file_names {
        set_args.push(' ');
        set_args.push_str(&shell_single_quote(file_name));
    }

    format!(
        "cd {} && {}; printf '[QCortado] Remote live log probe: workdir=%s\\n' \"$PWD\"; for f do if [ -e \"$f\" ]; then __qcortado_bytes=$(wc -c < \"$f\" | tr -d '[:space:]'); __qcortado_lines=$(wc -l < \"$f\" | tr -d '[:space:]'); printf '[QCortado] Remote live log probe: %s exists bytes=%s lines=%s\\n' \"$f\" \"$__qcortado_bytes\" \"$__qcortado_lines\"; else printf '[QCortado] Remote live log probe: %s missing\\n' \"$f\"; fi; done",
        shell_single_quote(remote_workdir),
        set_args
    )
}

fn parse_incremental_read_output(output: &str) -> Option<(String, u64)> {
    let (marker_line, content) = output.split_once('\n')?;
    let total_bytes = marker_line
        .strip_prefix(LIVE_SIZE_MARKER)?
        .trim()
        .parse::<u64>()
        .ok()?;
    Some((content.to_string(), total_bytes))
}

fn parse_linecount_read_output(output: &str) -> Option<(Vec<String>, usize)> {
    let mut lines: Vec<String> = output.lines().map(|line| line.to_string()).collect();
    let marker_line = lines.pop()?;
    let total_lines = marker_line
        .strip_prefix(LIVE_LINECOUNT_MARKER)?
        .trim()
        .parse::<usize>()
        .ok()?;
    Some((lines, total_lines))
}

fn take_complete_lines(buffer: &mut String) -> Vec<String> {
    let mut lines = Vec::new();
    let mut start = 0usize;
    let mut drain_to = 0usize;

    for (idx, ch) in buffer.char_indices() {
        if ch != '\n' {
            continue;
        }
        let mut line = &buffer[start..idx];
        if let Some(stripped) = line.strip_suffix('\r') {
            line = stripped;
        }
        lines.push(line.to_string());
        start = idx + ch.len_utf8();
        drain_to = start;
    }

    if drain_to > 0 {
        buffer.drain(..drain_to);
    }

    lines
}

async fn stream_remote_file_linecount(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    file_name: &str,
    last_line_count: &mut usize,
    already_announced: bool,
) -> bool {
    let read_cmd = build_linecount_read_command(remote_workdir, file_name, *last_line_count);
    let Ok(read_output) = run_ssh_command(profile, secret, &read_cmd).await else {
        return already_announced;
    };

    let Some((lines_to_emit, total_lines)) = parse_linecount_read_output(&read_output) else {
        return already_announced;
    };

    let mut announced = already_announced;
    if !lines_to_emit.is_empty() {
        if !announced {
            emit_task_line(
                app,
                pm,
                task_id,
                format!("--- Remote {} (live) ---", file_name),
            )
            .await;
            announced = true;
        }
        for line in lines_to_emit {
            emit_task_line(app, pm, task_id, line).await;
        }
    }

    *last_line_count = total_lines;
    announced
}

async fn stream_remote_file_incremental(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    file_name: &str,
    cursor: &mut LiveFileCursor,
    already_announced: bool,
    flush_partial: bool,
) -> bool {
    if !ENABLE_EXPERIMENTAL_REMOTE_LIVE_LOGGING {
        return stream_remote_file_linecount(
            app,
            pm,
            task_id,
            profile,
            secret,
            remote_workdir,
            file_name,
            &mut cursor.line_count,
            already_announced,
        )
        .await;
    }

    let read_cmd = build_incremental_read_command(remote_workdir, file_name, cursor.byte_offset);
    let read_output = match run_ssh_command(profile, secret, &read_cmd).await {
        Ok(output) => output,
        Err(err) => {
            if !cursor.warned_read_error {
                emit_task_line(
                    app,
                    pm,
                    task_id,
                    format!(
                        "HPC_WARNING|Live log read failed for {}: {}",
                        file_name, err
                    ),
                )
                .await;
                cursor.warned_read_error = true;
            }
            return already_announced;
        }
    };

    let Some((mut chunk, mut total_bytes)) = parse_incremental_read_output(&read_output) else {
        if !read_output.trim().is_empty() && !cursor.warned_parse_error {
            emit_task_line(
                app,
                pm,
                task_id,
                format!(
                    "HPC_WARNING|Live log read for {} returned an unexpected payload; skipping this poll.",
                    file_name
                ),
            )
            .await;
            cursor.warned_parse_error = true;
        }
        return already_announced;
    };

    let mut announced = already_announced;

    if total_bytes < cursor.byte_offset {
        emit_task_line(
            app,
            pm,
            task_id,
            format!(
                "HPC_WARNING|Detected log rotation/truncation for {}. Restarting live stream from beginning.",
                file_name
            ),
        )
        .await;

        let reset_cmd = build_incremental_read_command(remote_workdir, file_name, 0);
        let Ok(reset_output) = run_ssh_command(profile, secret, &reset_cmd).await else {
            return announced;
        };
        let Some((reset_chunk, reset_total)) = parse_incremental_read_output(&reset_output) else {
            return announced;
        };
        cursor.pending_line.clear();
        chunk = reset_chunk;
        total_bytes = reset_total;
    }

    if !chunk.is_empty() {
        cursor.pending_line.push_str(&chunk);
    }

    let mut lines_to_emit = take_complete_lines(&mut cursor.pending_line);
    if flush_partial && !cursor.pending_line.is_empty() {
        lines_to_emit.push(std::mem::take(&mut cursor.pending_line));
    }

    if !lines_to_emit.is_empty() {
        if file_name != "slurm.out" && !announced {
            emit_task_line(
                app,
                pm,
                task_id,
                format!("--- Remote {} (live) ---", file_name),
            )
            .await;
            announced = true;
        }
        for line in lines_to_emit {
            emit_task_line(app, pm, task_id, line).await;
        }
    }

    cursor.byte_offset = total_bytes;
    announced
}

async fn drain_remote_live_output(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    live_output_files: &[String],
    live_cursors: &mut HashMap<String, LiveFileCursor>,
    announced_live_files: &mut HashSet<String>,
) {
    for attempt in 0..3 {
        for file_name in live_output_files {
            let cursor = live_cursors.entry(file_name.clone()).or_default();
            let announced = announced_live_files.contains(file_name);
            let updated_announced = stream_remote_file_incremental(
                app,
                pm,
                task_id,
                profile,
                secret,
                remote_workdir,
                file_name,
                cursor,
                announced,
                false,
            )
            .await;
            if updated_announced {
                announced_live_files.insert(file_name.clone());
            }
        }
        if attempt < 2 {
            tokio::time::sleep(Duration::from_millis(750)).await;
        }
    }

    for file_name in live_output_files {
        let cursor = live_cursors.entry(file_name.clone()).or_default();
        let announced = announced_live_files.contains(file_name);
        let updated_announced = stream_remote_file_incremental(
            app,
            pm,
            task_id,
            profile,
            secret,
            remote_workdir,
            file_name,
            cursor,
            announced,
            true,
        )
        .await;
        if updated_announced {
            announced_live_files.insert(file_name.clone());
        }
    }
}

async fn emit_remote_live_probe(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    live_output_files: &[String],
) {
    emit_task_line(
        app,
        pm,
        task_id,
        format!(
            "[QCortado] Remote live streaming is watching: {}",
            live_output_files.join(", ")
        ),
    )
    .await;

    let probe_cmd = build_live_probe_command(remote_workdir, live_output_files);
    match run_ssh_command(profile, secret, &probe_cmd).await {
        Ok(output) => {
            for line in output.lines().filter(|line| !line.trim().is_empty()) {
                emit_task_line(app, pm, task_id, line.to_string()).await;
            }
        }
        Err(err) => {
            emit_task_line(
                app,
                pm,
                task_id,
                format!("HPC_WARNING|Remote live log probe failed: {}", err),
            )
            .await;
        }
    }
}

async fn resolve_remote_path(
    profile: &HpcProfile,
    secret: Option<&str>,
    raw_path: &str,
) -> Result<String, String> {
    let trimmed = raw_path.trim();
    if trimmed.is_empty() {
        return Err("Remote workspace path is empty".to_string());
    }

    if trimmed == "~" || trimmed.starts_with("~/") {
        let home = run_ssh_command(profile, secret, "printf %s \"$HOME\"").await?;
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

pub async fn run_batch_task(
    app: AppHandle,
    pm: ProcessManager,
    request: HpcBatchRequest,
) -> Result<HpcBatchResult, String> {
    let bundle_name = request
        .bundle_dir
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| "Invalid HPC bundle directory name".to_string())?
        .to_string();
    let workspace_root = resolve_remote_path(
        &request.profile,
        request.secret.as_deref(),
        request.profile.remote_workspace_root.trim_end_matches('/'),
    )
    .await?;
    let workspace_root = workspace_root.trim_end_matches('/').to_string();
    let remote_workdir = format!("{}/{}", workspace_root, bundle_name);
    let remote_project_root = resolve_remote_path(
        &request.profile,
        request.secret.as_deref(),
        request.profile.remote_project_root.trim_end_matches('/'),
    )
    .await
    .ok()
    .map(|value| value.trim_end_matches('/').to_string())
    .filter(|value| !value.is_empty());
    let planned_remote_project_path = remote_project_root.as_ref().map(|root| {
        format!(
            "{}/{}/{}",
            root.trim_end_matches('/'),
            request.task_kind,
            request.task_id
        )
    });
    let mut archived_remote_project_path: Option<String> = None;

    pm.set_task_backend(task_id(&request), Some("hpc".to_string()))
        .await;
    pm.set_remote_workdir(task_id(&request), Some(remote_workdir.clone()))
        .await;

    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        format!("HPC_CMD|{}", request.sbatch_preview),
    )
    .await;
    emit_task_line(&app, &pm, task_id(&request), "HPC_SCRIPT_BEGIN".to_string()).await;
    for line in request.slurm_script.lines() {
        emit_task_line(&app, &pm, task_id(&request), format!("HPC_SCRIPT|{}", line)).await;
    }
    emit_task_line(&app, &pm, task_id(&request), "HPC_SCRIPT_END".to_string()).await;

    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        format!("HPC_STAGE|Connecting|{}", request.profile.host),
    )
    .await;

    let mkdir_cmd = format!("mkdir -p {}", shell_single_quote(&workspace_root));
    run_ssh_command(&request.profile, request.secret.as_deref(), &mkdir_cmd).await?;

    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        format!("HPC_STAGE|Uploading|{}", remote_workdir),
    )
    .await;
    upload_directory(
        &request.profile,
        request.secret.as_deref(),
        &request.bundle_dir,
        &workspace_root,
    )
    .await?;

    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        format!("HPC_STAGE|Submitting|{}", request.sbatch_preview),
    )
    .await;
    let submit_cmd = format!(
        "cd {} && sbatch run.sbatch",
        shell_single_quote(&remote_workdir)
    );
    let submit_output =
        run_ssh_command(&request.profile, request.secret.as_deref(), &submit_cmd).await?;
    let job_id = parse_sbatch_job_id(&submit_output).ok_or_else(|| {
        format!(
            "Failed to parse job ID from sbatch output: {}",
            submit_output.trim()
        )
    })?;
    let job_id_record_cmd = format!(
        "cd {} && printf %s {} > .qcortado_job_id",
        shell_single_quote(&remote_workdir),
        shell_single_quote(&job_id)
    );
    let _ = run_ssh_command(
        &request.profile,
        request.secret.as_deref(),
        &job_id_record_cmd,
    )
    .await;

    pm.set_remote_job_id(task_id(&request), Some(job_id.clone()))
        .await;
    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        format!("HPC_STAGE|Submitted|{}", job_id),
    )
    .await;

    let mut last_scheduler_state: Option<String> = None;
    let live_output_files = collect_live_output_files(&request.slurm_script);
    let mut live_cursors: HashMap<String, LiveFileCursor> = HashMap::new();
    let mut announced_live_files: HashSet<String> = HashSet::new();
    let mut emitted_live_probe = false;

    let terminal_snapshot = loop {
        if request.cancel_flag.load(Ordering::SeqCst) {
            let cancel_cmd = format!("scancel {}", shell_single_quote(&job_id));
            let _ = run_ssh_command(&request.profile, request.secret.as_deref(), &cancel_cmd).await;
            return Err("Cancelled by user".to_string());
        }

        for file_name in &live_output_files {
            let cursor = live_cursors.entry(file_name.clone()).or_default();
            let announced = announced_live_files.contains(file_name);
            let updated_announced = stream_remote_file_incremental(
                &app,
                &pm,
                task_id(&request),
                &request.profile,
                request.secret.as_deref(),
                &remote_workdir,
                file_name,
                cursor,
                announced,
                false,
            )
            .await;
            if updated_announced {
                announced_live_files.insert(file_name.clone());
            }
        }

        let squeue_cmd = format!("squeue -h -j {} -o \"%T|%N\"", shell_single_quote(&job_id));
        let squeue_output =
            run_ssh_command(&request.profile, request.secret.as_deref(), &squeue_cmd).await;
        let scheduler_snapshot = match squeue_output {
            Ok(output) => parse_squeue_snapshot(&output),
            Err(_) => None,
        };

        let resolved_snapshot = if let Some(snapshot) = scheduler_snapshot {
            Some(SchedulerSnapshot {
                state: normalize_scheduler_state(&snapshot.state),
                node: snapshot.node.clone(),
                source: snapshot.source,
            })
        } else {
            let sacct_cmd = format!(
                "sacct -j {} --format=State,NodeList --parsable2 --noheader",
                shell_single_quote(&job_id)
            );
            match run_ssh_command(&request.profile, request.secret.as_deref(), &sacct_cmd).await {
                Ok(output) => parse_sacct_snapshot(&output).map(|snapshot| SchedulerSnapshot {
                    state: normalize_scheduler_state(&snapshot.state),
                    node: snapshot.node.clone(),
                    source: snapshot.source,
                }),
                Err(_) => None,
            }
        };

        if let Some(snapshot) = resolved_snapshot {
            let changed = last_scheduler_state
                .as_ref()
                .map(|state| state != &snapshot.state)
                .unwrap_or(true);
            if changed {
                update_scheduler_snapshot(&app, &pm, task_id(&request), &snapshot).await;
                last_scheduler_state = Some(snapshot.state.clone());
            }
            if ENABLE_EXPERIMENTAL_REMOTE_LIVE_LOGGING
                && snapshot.state == "RUNNING"
                && !emitted_live_probe
            {
                emit_remote_live_probe(
                    &app,
                    &pm,
                    task_id(&request),
                    &request.profile,
                    request.secret.as_deref(),
                    &remote_workdir,
                    &live_output_files,
                )
                .await;
                emitted_live_probe = true;
            }
            if is_terminal_state(&snapshot.state) {
                break snapshot;
            }
        }

        tokio::time::sleep(Duration::from_secs(4)).await;
    };
    let terminal_state = normalize_scheduler_state(&terminal_snapshot.state);
    drain_remote_live_output(
        &app,
        &pm,
        task_id(&request),
        &request.profile,
        request.secret.as_deref(),
        &remote_workdir,
        &live_output_files,
        &mut live_cursors,
        &mut announced_live_files,
    )
    .await;

    if terminal_state == "COMPLETED" {
        let completion_sync_mode = if request.task_kind == "epw" {
            ArtifactSyncMode::EpwReady
        } else {
            ArtifactSyncMode::Minimal
        };
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!(
                "HPC_STAGE|Collecting|{} ({})",
                remote_workdir,
                completion_sync_mode.label()
            ),
        )
        .await;
        let sync_started = Instant::now();
        let sync_report = sync_remote_artifacts(
            &app,
            &pm,
            HpcArtifactSyncRequest {
                task_id: request.task_id.clone(),
                task_kind: request.task_kind.clone(),
                profile: request.profile.clone(),
                secret: request.secret.clone(),
                remote_workdir: remote_workdir.clone(),
                local_sync_dir: request.local_sync_dir.clone(),
                mode: completion_sync_mode,
            },
        )
        .await?;
        let sync_title = if matches!(completion_sync_mode, ArtifactSyncMode::EpwReady) {
            "EPW-ready"
        } else {
            "Minimal"
        };
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!(
                "[QCortado] {} artifact sync finished in {:.1}s.",
                sync_title,
                sync_started.elapsed().as_secs_f64()
            ),
        )
        .await;
        let remote_storage_bytes = sync_report
            .downloaded_bytes
            .saturating_add(sync_report.skipped_bytes);
        pm.set_remote_storage_bytes(task_id(&request), Some(remote_storage_bytes))
            .await;
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!(
                "HPC_STAGE|Saved|{} sync complete ({} files, {:.2} MB downloaded, {} skipped, remote {:.2} MB)",
                sync_report.mode,
                sync_report.downloaded_files,
                sync_report.downloaded_bytes as f64 / (1024.0 * 1024.0),
                sync_report.skipped_files,
                remote_storage_bytes as f64 / (1024.0 * 1024.0),
            ),
        )
        .await;
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            if matches!(completion_sync_mode, ArtifactSyncMode::EpwReady) {
                "HPC_WARNING|Additional scratch artifacts remain remote after EPW-ready sync. Use 'Download full bundle' if needed."
                    .to_string()
            } else {
                "HPC_WARNING|Large scratch artifacts remain remote. Use 'Download full bundle' if needed."
                    .to_string()
            },
        )
        .await;
        if let Some(remote_project_path) = planned_remote_project_path.as_ref() {
            if archive_remote_run_copy(
                &app,
                &pm,
                task_id(&request),
                &request.profile,
                request.secret.as_deref(),
                &remote_workdir,
                remote_project_path,
            )
            .await
            .is_ok()
            {
                archived_remote_project_path = Some(remote_project_path.clone());
                pm.set_remote_project_path(task_id(&request), Some(remote_project_path.clone()))
                    .await;
                if remote_project_path != &remote_workdir {
                    remove_remote_paths_best_effort(
                        &app,
                        &pm,
                        task_id(&request),
                        &request.profile,
                        request.secret.as_deref(),
                        std::slice::from_ref(&remote_workdir),
                        "workspace bundle",
                    )
                    .await;
                }
                pm.set_remote_workdir(task_id(&request), Some(remote_project_path.clone()))
                    .await;
            }
        } else {
            emit_task_line(
                &app,
                &pm,
                task_id(&request),
                "HPC_WARNING|Remote project root unavailable; skipping remote archive copy"
                    .to_string(),
            )
            .await;
        }
    } else {
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!("HPC_WARNING|Remote job ended with state {}", terminal_state),
        )
        .await;

        emit_remote_tail_if_exists(
            &app,
            &pm,
            task_id(&request),
            &request.profile,
            request.secret.as_deref(),
            &remote_workdir,
            "slurm.err",
            220,
        )
        .await;
        emit_remote_tail_if_exists(
            &app,
            &pm,
            task_id(&request),
            &request.profile,
            request.secret.as_deref(),
            &remote_workdir,
            "slurm.out",
            220,
        )
        .await;
        emit_remote_tail_if_exists(
            &app,
            &pm,
            task_id(&request),
            &request.profile,
            request.secret.as_deref(),
            &remote_workdir,
            "pw.out",
            220,
        )
        .await;

        let _ = sync_remote_artifacts(
            &app,
            &pm,
            HpcArtifactSyncRequest {
                task_id: request.task_id.clone(),
                task_kind: request.task_kind.clone(),
                profile: request.profile.clone(),
                secret: request.secret.clone(),
                remote_workdir: remote_workdir.clone(),
                local_sync_dir: request.local_sync_dir.clone(),
                mode: ArtifactSyncMode::Minimal,
            },
        )
        .await;
        let mut cleanup_targets: Vec<String> = Vec::new();
        push_unique_remote_path(&mut cleanup_targets, Some(&remote_workdir));
        push_unique_remote_path(
            &mut cleanup_targets,
            archived_remote_project_path.as_deref(),
        );
        remove_remote_paths_best_effort(
            &app,
            &pm,
            task_id(&request),
            &request.profile,
            request.secret.as_deref(),
            &cleanup_targets,
            "failure artifacts",
        )
        .await;
        pm.set_remote_project_path(task_id(&request), None).await;
        pm.set_remote_workdir(task_id(&request), None).await;
        pm.set_remote_storage_bytes(task_id(&request), Some(0))
            .await;

        return Err(format!(
            "Remote job ended with state {}. Check slurm.err/slurm.out in task output.",
            terminal_state
        ));
    }

    Ok(HpcBatchResult {
        backend: "hpc".to_string(),
        task_kind: request.task_kind,
        task_label: request.task_label,
        remote_job_id: job_id,
        remote_workdir: archived_remote_project_path
            .clone()
            .unwrap_or_else(|| remote_workdir.clone()),
        remote_project_path: archived_remote_project_path,
        remote_node: terminal_snapshot.node,
        terminal_state,
        sbatch_preview: request.sbatch_preview,
        slurm_script: request.slurm_script,
    })
}

pub async fn run_attached_batch_task(
    app: AppHandle,
    pm: ProcessManager,
    request: HpcAttachRequest,
) -> Result<HpcBatchResult, String> {
    let task_id = request.task_id.as_str();
    let remote_workdir = request.remote_workdir.trim_end_matches('/').to_string();
    let job_id = request.remote_job_id.trim().to_string();
    if job_id.is_empty() {
        return Err("Remote job id is required for reattach.".to_string());
    }
    if remote_workdir.is_empty() {
        return Err("Remote working directory is required for reattach.".to_string());
    }

    emit_task_line(
        &app,
        &pm,
        task_id,
        format!("HPC_STAGE|Reattached|{} ({})", job_id, remote_workdir),
    )
    .await;
    if let Some(preview) = request
        .sbatch_preview
        .as_deref()
        .filter(|value| !value.trim().is_empty())
    {
        emit_task_line(&app, &pm, task_id, format!("HPC_CMD|{}", preview)).await;
    }
    if let Some(script) = request
        .slurm_script
        .as_deref()
        .filter(|value| !value.trim().is_empty())
    {
        emit_task_line(&app, &pm, task_id, "HPC_SCRIPT_BEGIN".to_string()).await;
        for line in script.lines() {
            emit_task_line(&app, &pm, task_id, format!("HPC_SCRIPT|{}", line)).await;
        }
        emit_task_line(&app, &pm, task_id, "HPC_SCRIPT_END".to_string()).await;
    }

    let live_output_files = request
        .slurm_script
        .as_deref()
        .map(collect_live_output_files)
        .filter(|files| !files.is_empty())
        .unwrap_or_else(|| vec!["slurm.out".to_string(), "slurm.err".to_string()]);
    let mut live_cursors: HashMap<String, LiveFileCursor> = HashMap::new();
    let mut announced_live_files: HashSet<String> = HashSet::new();
    let mut last_scheduler_state: Option<String> = None;
    let mut emitted_live_probe = false;

    let terminal_snapshot = loop {
        if request.cancel_flag.load(Ordering::SeqCst) {
            let cancel_cmd = format!("scancel {}", shell_single_quote(&job_id));
            let _ = run_ssh_command(&request.profile, request.secret.as_deref(), &cancel_cmd).await;
            return Err("Cancelled by user".to_string());
        }

        for file_name in &live_output_files {
            let cursor = live_cursors.entry(file_name.clone()).or_default();
            let announced = announced_live_files.contains(file_name);
            let updated_announced = stream_remote_file_incremental(
                &app,
                &pm,
                task_id,
                &request.profile,
                request.secret.as_deref(),
                &remote_workdir,
                file_name,
                cursor,
                announced,
                false,
            )
            .await;
            if updated_announced {
                announced_live_files.insert(file_name.clone());
            }
        }

        let squeue_cmd = format!("squeue -h -j {} -o \"%T|%N\"", shell_single_quote(&job_id));
        let scheduler_snapshot =
            run_ssh_command(&request.profile, request.secret.as_deref(), &squeue_cmd)
                .await
                .ok()
                .and_then(|output| parse_squeue_snapshot(&output));

        let resolved_snapshot = if let Some(snapshot) = scheduler_snapshot {
            Some(SchedulerSnapshot {
                state: normalize_scheduler_state(&snapshot.state),
                node: snapshot.node.clone(),
                source: snapshot.source,
            })
        } else {
            let sacct_cmd = format!(
                "sacct -j {} --format=State,NodeList --parsable2 --noheader",
                shell_single_quote(&job_id)
            );
            run_ssh_command(&request.profile, request.secret.as_deref(), &sacct_cmd)
                .await
                .ok()
                .and_then(|output| parse_sacct_snapshot(&output))
                .map(|snapshot| SchedulerSnapshot {
                    state: normalize_scheduler_state(&snapshot.state),
                    node: snapshot.node.clone(),
                    source: snapshot.source,
                })
        };

        if let Some(snapshot) = resolved_snapshot {
            let changed = last_scheduler_state
                .as_ref()
                .map(|state| state != &snapshot.state)
                .unwrap_or(true);
            if changed {
                update_scheduler_snapshot(&app, &pm, task_id, &snapshot).await;
                last_scheduler_state = Some(snapshot.state.clone());
            }
            if ENABLE_EXPERIMENTAL_REMOTE_LIVE_LOGGING
                && snapshot.state == "RUNNING"
                && !emitted_live_probe
            {
                emit_remote_live_probe(
                    &app,
                    &pm,
                    task_id,
                    &request.profile,
                    request.secret.as_deref(),
                    &remote_workdir,
                    &live_output_files,
                )
                .await;
                emitted_live_probe = true;
            }
            if is_terminal_state(&snapshot.state) {
                break snapshot;
            }
        }

        tokio::time::sleep(Duration::from_secs(4)).await;
    };

    let terminal_state = normalize_scheduler_state(&terminal_snapshot.state);
    drain_remote_live_output(
        &app,
        &pm,
        task_id,
        &request.profile,
        request.secret.as_deref(),
        &remote_workdir,
        &live_output_files,
        &mut live_cursors,
        &mut announced_live_files,
    )
    .await;
    let planned_remote_project_path = if let Some(path) = request.remote_project_path.as_ref() {
        Some(path.clone())
    } else {
        resolve_remote_path(
            &request.profile,
            request.secret.as_deref(),
            request.profile.remote_project_root.trim_end_matches('/'),
        )
        .await
        .ok()
        .map(|root| root.trim_end_matches('/').to_string())
        .filter(|value| !value.is_empty())
        .map(|root| {
            format!(
                "{}/{}/{}",
                root.trim_end_matches('/'),
                request.task_kind,
                request.task_id
            )
        })
    };
    let mut archived_remote_project_path = request.remote_project_path.clone();
    if let Some(remote_project_path) = archived_remote_project_path.as_ref() {
        pm.set_remote_project_path(task_id, Some(remote_project_path.clone()))
            .await;
    }

    if terminal_state == "COMPLETED" {
        let completion_sync_mode = if request.task_kind == "epw" {
            ArtifactSyncMode::EpwReady
        } else {
            ArtifactSyncMode::Minimal
        };
        emit_task_line(
            &app,
            &pm,
            task_id,
            format!(
                "HPC_STAGE|Collecting|{} ({})",
                remote_workdir,
                completion_sync_mode.label()
            ),
        )
        .await;
        let sync_report = sync_remote_artifacts(
            &app,
            &pm,
            HpcArtifactSyncRequest {
                task_id: request.task_id.clone(),
                task_kind: request.task_kind.clone(),
                profile: request.profile.clone(),
                secret: request.secret.clone(),
                remote_workdir: remote_workdir.clone(),
                local_sync_dir: request.local_sync_dir.clone(),
                mode: completion_sync_mode,
            },
        )
        .await?;
        let remote_storage_bytes = sync_report
            .downloaded_bytes
            .saturating_add(sync_report.skipped_bytes);
        pm.set_remote_storage_bytes(task_id, Some(remote_storage_bytes))
            .await;
        emit_task_line(
            &app,
            &pm,
            task_id,
            format!(
                "HPC_STAGE|Saved|{} sync complete ({} files, {:.2} MB downloaded, {} skipped, remote {:.2} MB)",
                sync_report.mode,
                sync_report.downloaded_files,
                sync_report.downloaded_bytes as f64 / (1024.0 * 1024.0),
                sync_report.skipped_files,
                remote_storage_bytes as f64 / (1024.0 * 1024.0),
            ),
        )
        .await;
        if let Some(remote_project_path) = planned_remote_project_path.as_ref() {
            if archive_remote_run_copy(
                &app,
                &pm,
                task_id,
                &request.profile,
                request.secret.as_deref(),
                &remote_workdir,
                remote_project_path,
            )
            .await
            .is_ok()
            {
                archived_remote_project_path = Some(remote_project_path.clone());
                pm.set_remote_project_path(task_id, Some(remote_project_path.clone()))
                    .await;
                if remote_project_path != &remote_workdir {
                    remove_remote_paths_best_effort(
                        &app,
                        &pm,
                        task_id,
                        &request.profile,
                        request.secret.as_deref(),
                        std::slice::from_ref(&remote_workdir),
                        "workspace bundle",
                    )
                    .await;
                }
                pm.set_remote_workdir(task_id, Some(remote_project_path.clone()))
                    .await;
            }
        } else {
            emit_task_line(
                &app,
                &pm,
                task_id,
                "HPC_WARNING|Remote project root unavailable; skipping remote archive copy"
                    .to_string(),
            )
            .await;
        }
    } else {
        emit_task_line(
            &app,
            &pm,
            task_id,
            format!("HPC_WARNING|Remote job ended with state {}", terminal_state),
        )
        .await;
        let _ = sync_remote_artifacts(
            &app,
            &pm,
            HpcArtifactSyncRequest {
                task_id: request.task_id.clone(),
                task_kind: request.task_kind.clone(),
                profile: request.profile.clone(),
                secret: request.secret.clone(),
                remote_workdir: remote_workdir.clone(),
                local_sync_dir: request.local_sync_dir.clone(),
                mode: ArtifactSyncMode::Minimal,
            },
        )
        .await;
        let mut cleanup_targets: Vec<String> = Vec::new();
        push_unique_remote_path(&mut cleanup_targets, Some(&remote_workdir));
        push_unique_remote_path(
            &mut cleanup_targets,
            archived_remote_project_path.as_deref(),
        );
        remove_remote_paths_best_effort(
            &app,
            &pm,
            task_id,
            &request.profile,
            request.secret.as_deref(),
            &cleanup_targets,
            "failure artifacts",
        )
        .await;
        pm.set_remote_project_path(task_id, None).await;
        pm.set_remote_workdir(task_id, None).await;
        pm.set_remote_storage_bytes(task_id, Some(0)).await;
        return Err(format!("Remote job ended with state {}.", terminal_state));
    }

    Ok(HpcBatchResult {
        backend: "hpc".to_string(),
        task_kind: request.task_kind,
        task_label: request.task_label,
        remote_job_id: job_id,
        remote_workdir: archived_remote_project_path
            .clone()
            .unwrap_or_else(|| remote_workdir.clone()),
        remote_project_path: archived_remote_project_path.take(),
        remote_node: terminal_snapshot.node,
        terminal_state,
        sbatch_preview: request.sbatch_preview.unwrap_or_default(),
        slurm_script: request.slurm_script.unwrap_or_default(),
    })
}

fn task_id(request: &HpcBatchRequest) -> &str {
    request.task_id.as_str()
}

#[cfg(test)]
mod tests {
    use super::{
        collect_live_output_files, collect_option_targets, parse_incremental_read_output,
        push_unique_remote_path, should_download_minimal, take_complete_lines, RemoteFileEntry,
        LIVE_SIZE_MARKER,
    };

    #[test]
    fn parses_incremental_payload_with_content() {
        let payload = format!("{}42\nline-1\nline-2\n", LIVE_SIZE_MARKER);
        let parsed = parse_incremental_read_output(&payload).expect("expected marker payload");
        assert_eq!(parsed.0, "line-1\nline-2\n".to_string());
        assert_eq!(parsed.1, 42);
    }

    #[test]
    fn parses_incremental_payload_for_empty_file() {
        let payload = format!("{}0\n", LIVE_SIZE_MARKER);
        let parsed = parse_incremental_read_output(&payload).expect("expected marker payload");
        assert_eq!(parsed.0, "");
        assert_eq!(parsed.1, 0);
    }

    #[test]
    fn rejects_payload_without_marker() {
        assert!(parse_incremental_read_output("line-1\nline-2\n").is_none());
    }

    #[test]
    fn keeps_partial_line_until_newline() {
        let mut buffer = "rank 0 start".to_string();
        assert!(take_complete_lines(&mut buffer).is_empty());
        assert_eq!(buffer, "rank 0 start");

        buffer.push_str(" done\nrank 1");
        assert_eq!(
            take_complete_lines(&mut buffer),
            vec!["rank 0 start done".to_string()]
        );
        assert_eq!(buffer, "rank 1");
    }

    #[test]
    fn live_output_files_include_multiple_redirects_on_one_line() {
        let files = collect_live_output_files(
            r#"
#!/bin/bash
srun "$EPW_EXE" -in epw.in > epw.out 2> epw.err
mpirun "$QE_BIN/pw.x" -in pw.in > pw.out 2>&1
postw90.x seed > seed.wpout 2> seed.werr
if mpirun -help 2>&1 | grep -q -- '-outfile-pattern'; then mpirun -outfile-pattern .qcortado-live/hp.out.rank.%r -errfile-pattern .qcortado-live/hp.out.err.rank.%r "$QE_BIN/hp.x"; fi
"#,
        );
        assert!(files.contains(&"epw.out".to_string()));
        assert!(files.contains(&"epw.err".to_string()));
        assert!(files.contains(&"pw.out".to_string()));
        assert!(files.contains(&"seed.wpout".to_string()));
        assert!(files.contains(&"seed.werr".to_string()));
        assert!(!files.contains(&".qcortado-live/hp.out.rank.0".to_string()));
        assert!(!files.contains(&".qcortado-live/hp.out.err.rank.0".to_string()));
        assert!(!files.iter().any(|file| file.contains("; then mpirun")));
    }

    #[test]
    fn shell_words_keep_quoted_option_probe_separate_from_real_option() {
        let targets = collect_option_targets(
            "if mpirun -help 2>&1 | grep -q -- '-outfile-pattern'; then mpirun -outfile-pattern '.qcortado-live/hp.out.rank.%r'; fi",
            "-outfile-pattern",
        );
        assert_eq!(targets, vec![".qcortado-live/hp.out.rank.%r".to_string()]);
    }

    #[test]
    fn live_probe_command_reports_watched_files() {
        let command = super::build_live_probe_command("/remote/work dir", &["pw.out".to_string()]);
        assert!(command.contains("cd '/remote/work dir'"));
        assert!(command.contains("set -- 'pw.out'"));
        assert!(command.contains("Remote live log probe"));
    }

    #[test]
    fn minimal_sync_keeps_surface_files_from_scratch_paths() {
        let entry = RemoteFileEntry {
            rel_path: "tmp/si.save/vfermi.frmsf".to_string(),
            size_bytes: 1_024,
        };
        assert!(should_download_minimal("fermi_surface", &entry));
    }

    #[test]
    fn minimal_sync_skips_non_surface_files_from_scratch_paths() {
        let entry = RemoteFileEntry {
            rel_path: "tmp/si.save/prefix.rho".to_string(),
            size_bytes: 1_024,
        };
        assert!(!should_download_minimal("fermi_surface", &entry));
    }

    #[test]
    fn cleanup_targets_are_unique_and_ignore_empty_paths() {
        let mut targets = Vec::new();
        push_unique_remote_path(&mut targets, Some("/remote/workdir"));
        push_unique_remote_path(&mut targets, Some("/remote/workdir"));
        push_unique_remote_path(&mut targets, Some("  "));
        push_unique_remote_path(&mut targets, None);
        push_unique_remote_path(&mut targets, Some("/remote/archive"));
        assert_eq!(
            targets,
            vec!["/remote/workdir".to_string(), "/remote/archive".to_string()]
        );
    }
}
