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
    Full,
}

impl ArtifactSyncMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Minimal => "minimal",
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
    if is_heavy_scratch_path(&entry.rel_path) {
        return false;
    }

    let lower_rel = entry.rel_path.to_ascii_lowercase();
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

    if lower_rel.ends_with(".frmsf")
        || lower_rel.ends_with(".bxsf")
        || lower_rel.ends_with(".gnu")
        || lower_rel.ends_with(".json")
        || lower_rel.ends_with(".txt")
        || lower_rel.ends_with(".log")
    {
        return true;
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

    match request.mode {
        ArtifactSyncMode::Full => {
            let remote_files = list_remote_files(
                &request.profile,
                request.secret.as_deref(),
                &request.remote_workdir,
            )
            .await?;
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
        ArtifactSyncMode::Minimal => {
            let remote_files = list_remote_files(
                &request.profile,
                request.secret.as_deref(),
                &request.remote_workdir,
            )
            .await?;
            let mut candidates: Vec<RemoteFileEntry> = Vec::new();
            for entry in remote_files {
                if should_download_minimal(&request.task_kind, &entry) {
                    candidates.push(entry);
                } else {
                    report.skipped_files += 1;
                    report.skipped_bytes = report.skipped_bytes.saturating_add(entry.size_bytes);
                }
            }

            let mut failures: Vec<String> = Vec::new();
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
                        "HPC_WARNING|Some artifacts were not downloaded during minimal sync ({}).",
                        preview
                    ),
                )
                .await;
            }

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

fn collect_live_output_files(slurm_script: &str) -> Vec<String> {
    let mut files = vec!["slurm.out".to_string(), "slurm.err".to_string()];
    let mut seen: HashSet<String> = files.iter().cloned().collect();

    for line in slurm_script.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let Some((_, redirect_rhs)) = trimmed.split_once('>') else {
            continue;
        };
        let Some(target) = parse_redirect_target(redirect_rhs) else {
            continue;
        };
        if target == "&1" || target.starts_with('&') {
            continue;
        }
        if !(target.ends_with(".out") || target.ends_with(".err")) {
            continue;
        }
        if seen.insert(target.clone()) {
            files.push(target);
        }
    }

    files
}

async fn stream_remote_file_tail_incremental(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    file_name: &str,
    max_lines: usize,
    last_lines: &mut Vec<String>,
    already_announced: bool,
) -> bool {
    let cmd = format!(
        "cd {} && if [ -f {} ]; then tail -n {} {}; fi",
        shell_single_quote(remote_workdir),
        shell_single_quote(file_name),
        max_lines.max(1),
        shell_single_quote(file_name)
    );
    let Ok(tail_output) = run_ssh_command(profile, secret, &cmd).await else {
        return already_announced;
    };

    let next_lines: Vec<String> = tail_output.lines().map(|line| line.to_string()).collect();
    let mut announced = already_announced;

    let lines_to_emit: Vec<String> = if next_lines.len() >= last_lines.len() {
        next_lines
            .iter()
            .skip(last_lines.len())
            .cloned()
            .collect::<Vec<_>>()
    } else {
        next_lines.clone()
    };

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

    *last_lines = next_lines;
    announced
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
    let mut live_tail_history: HashMap<String, Vec<String>> = HashMap::new();
    let mut announced_live_files: HashSet<String> = HashSet::new();

    let terminal_snapshot = loop {
        if request.cancel_flag.load(Ordering::SeqCst) {
            let cancel_cmd = format!("scancel {}", shell_single_quote(&job_id));
            let _ = run_ssh_command(&request.profile, request.secret.as_deref(), &cancel_cmd).await;
            return Err("Cancelled by user".to_string());
        }

        for file_name in &live_output_files {
            let previous_lines = live_tail_history
                .entry(file_name.clone())
                .or_insert_with(Vec::new);
            let announced = announced_live_files.contains(file_name);
            let updated_announced = stream_remote_file_tail_incremental(
                &app,
                &pm,
                task_id(&request),
                &request.profile,
                request.secret.as_deref(),
                &remote_workdir,
                file_name,
                500,
                previous_lines,
                announced,
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
            if is_terminal_state(&snapshot.state) {
                break snapshot;
            }
        }

        tokio::time::sleep(Duration::from_secs(4)).await;
    };
    let terminal_state = normalize_scheduler_state(&terminal_snapshot.state);

    if let Some(remote_project_path) = planned_remote_project_path.as_ref() {
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!("HPC_STAGE|Archiving|{}", remote_project_path),
        )
        .await;

        let archive_cmd = format!(
            "mkdir -p {dest} && if [ -d {src} ]; then cp -a {src}/. {dest}/; fi",
            src = shell_single_quote(&remote_workdir),
            dest = shell_single_quote(remote_project_path)
        );
        match run_ssh_command(&request.profile, request.secret.as_deref(), &archive_cmd).await {
            Ok(_) => {
                archived_remote_project_path = Some(remote_project_path.clone());
                pm.set_remote_project_path(task_id(&request), Some(remote_project_path.clone()))
                    .await;
                emit_task_line(
                    &app,
                    &pm,
                    task_id(&request),
                    format!("HPC_STAGE|Archived|{}", remote_project_path),
                )
                .await;
            }
            Err(err) => {
                emit_task_line(
                    &app,
                    &pm,
                    task_id(&request),
                    format!(
                        "HPC_WARNING|Failed to archive run under remote project root ({}): {}",
                        remote_project_path, err
                    ),
                )
                .await;
            }
        }
    } else {
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            "HPC_WARNING|Remote project root unavailable; skipping remote archive copy".to_string(),
        )
        .await;
    }

    if terminal_state == "COMPLETED" {
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!("HPC_STAGE|Collecting|{} (minimal)", remote_workdir),
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
                mode: ArtifactSyncMode::Minimal,
            },
        )
        .await?;
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
                "HPC_STAGE|Saved|Minimal sync complete ({} files, {:.2} MB downloaded, {} skipped, remote {:.2} MB)",
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
            "HPC_WARNING|Large scratch artifacts remain remote. Use 'Download full bundle' if needed."
                .to_string(),
        )
        .await;
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
        remote_workdir,
        remote_project_path: archived_remote_project_path,
        remote_node: terminal_snapshot.node,
        terminal_state,
        sbatch_preview: request.sbatch_preview,
        slurm_script: request.slurm_script,
    })
}

fn task_id(request: &HpcBatchRequest) -> &str {
    request.task_id.as_str()
}
