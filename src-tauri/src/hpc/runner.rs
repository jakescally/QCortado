use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::Ordering;
use std::time::Duration;

use tauri::{AppHandle, Emitter};

use crate::process_manager::ProcessManager;

use super::profile::HpcProfile;
use super::slurm::{
    SchedulerSnapshot, is_terminal_state, normalize_scheduler_state, parse_sacct_snapshot,
    parse_sbatch_job_id, parse_squeue_snapshot,
};
use super::ssh::{download_directory_contents, run_ssh_command, upload_directory};

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
    pub remote_node: Option<String>,
    pub terminal_state: String,
    pub sbatch_preview: String,
    pub slurm_script: String,
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

async fn update_scheduler_snapshot(
    app: &AppHandle,
    pm: &ProcessManager,
    task_id: &str,
    snapshot: &SchedulerSnapshot,
) {
    pm.set_scheduler_state(task_id, Some(snapshot.state.clone())).await;
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
    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        "HPC_SCRIPT_BEGIN".to_string(),
    )
    .await;
    for line in request.slurm_script.lines() {
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!("HPC_SCRIPT|{}", line),
        )
        .await;
    }
    emit_task_line(
        &app,
        &pm,
        task_id(&request),
        "HPC_SCRIPT_END".to_string(),
    )
    .await;

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
    let submit_output = run_ssh_command(&request.profile, request.secret.as_deref(), &submit_cmd).await?;
    let job_id = parse_sbatch_job_id(&submit_output)
        .ok_or_else(|| format!("Failed to parse job ID from sbatch output: {}", submit_output.trim()))?;

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

    if terminal_state == "COMPLETED" {
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            format!("HPC_STAGE|Collecting|{}", remote_workdir),
        )
        .await;
        download_directory_contents(
            &request.profile,
            request.secret.as_deref(),
            &remote_workdir,
            &request.local_sync_dir,
        )
        .await?;
        emit_task_line(
            &app,
            &pm,
            task_id(&request),
            "HPC_STAGE|Saved|Local sync complete".to_string(),
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

        let _ = download_directory_contents(
            &request.profile,
            request.secret.as_deref(),
            &remote_workdir,
            &request.local_sync_dir,
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
        remote_node: terminal_snapshot.node,
        terminal_state,
        sbatch_preview: request.sbatch_preview,
        slurm_script: request.slurm_script,
    })
}

fn task_id(request: &HpcBatchRequest) -> &str {
    request.task_id.as_str()
}
