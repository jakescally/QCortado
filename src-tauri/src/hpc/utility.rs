//! Lightweight scheduled operations for engine setup and validation.
//!
//! Login-node SSH is used only for staging, submission, polling and output
//! retrieval. Commands that load an engine environment run inside the Slurm
//! allocation so module-only profiles work on clusters that prohibit modules
//! on login nodes.

use std::time::{Duration, Instant};

use tauri::{AppHandle, Emitter};

use super::profile::{HpcProfile, HpcResourceMode};
use super::slurm::{
    build_slurm_script, is_terminal_state, normalize_scheduler_state, parse_sacct_snapshot,
    parse_sbatch_job_id, parse_squeue_snapshot, SlurmScript,
};
use super::ssh::{run_ssh_command_with_timeout, upload_file};

#[derive(Debug, Clone)]
pub struct ScheduledUtilityResult {
    pub job_id: String,
    pub output: String,
    pub slurm_script: String,
}

fn shell_single_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

/// Shell fragment that bootstraps a typical environment-modules or Lmod
/// installation in non-interactive batch shells.
pub fn module_environment_bootstrap_command() -> String {
    "for module_init in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do if [ -f \"$module_init\" ]; then . \"$module_init\"; break; fi; done".to_string()
}

fn emit_line(app: Option<&AppHandle>, event_name: Option<&str>, line: &str) {
    if let (Some(app), Some(event_name)) = (app, event_name) {
        let _ = app.emit(event_name, line);
    }
}

fn emit_incremental_output(
    app: Option<&AppHandle>,
    event_name: Option<&str>,
    latest: &str,
    previous: &mut String,
    prefix: &str,
) {
    let newly_visible = latest.strip_prefix(previous.as_str()).unwrap_or(latest);
    for line in newly_visible.lines() {
        if !line.trim().is_empty() {
            emit_line(app, event_name, &format!("{}{}", prefix, line));
        }
    }
    *previous = latest.to_string();
}

async fn read_remote_output(profile: &HpcProfile, secret: Option<&str>, path: &str) -> String {
    run_ssh_command_with_timeout(
        profile,
        secret,
        &format!(
            "test -f {path} && cat {path} || true",
            path = shell_single_quote(path)
        ),
        20,
    )
    .await
    .unwrap_or_default()
}

pub fn build_scheduled_utility_script(
    profile: &HpcProfile,
    job_name: &str,
    commands: &[String],
) -> Result<SlurmScript, String> {
    let resources = profile.utility_resources();
    let mut allocation_profile = profile.clone();
    allocation_profile.resource_mode = HpcResourceMode::Both;
    allocation_profile.default_cpu_resources = resources.clone();
    let script = build_slurm_script(&allocation_profile, job_name, commands, Some(resources));
    if !script.validation.errors.is_empty() {
        return Err(format!(
            "Utility-job resource settings are invalid: {}",
            script.validation.errors.join(" ")
        ));
    }
    Ok(script)
}

/// Run a short non-calculation command sequence through Slurm without
/// registering it with the task queue.
pub async fn run_scheduled_utility_operation(
    app: Option<&AppHandle>,
    event_name: Option<&str>,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_workdir: &str,
    job_name: &str,
    commands: &[String],
    timeout_secs: u64,
) -> Result<ScheduledUtilityResult, String> {
    let script = build_scheduled_utility_script(profile, job_name, commands)?;

    let token = uuid::Uuid::new_v4().to_string();
    let local_script = std::env::temp_dir().join(format!("qcortado_utility_{token}.sbatch"));
    std::fs::write(&local_script, &script.script)
        .map_err(|err| format!("Failed to write scheduled utility script: {err}"))?;
    let remote_script = format!("{}/.qcortado_utility_{}.sbatch", remote_workdir, token);
    let upload_result = upload_file(profile, secret, &local_script, &remote_script).await;
    let _ = std::fs::remove_file(&local_script);
    upload_result?;

    let stdout_path = format!("{}/slurm.out", remote_workdir);
    let stderr_path = format!("{}/slurm.err", remote_workdir);
    let clear_outputs = format!(
        "cd {dir} && rm -f slurm.out slurm.err && sbatch {script}",
        dir = shell_single_quote(remote_workdir),
        script = shell_single_quote(&remote_script),
    );
    emit_line(
        app,
        event_name,
        &format!(
            "[QCortado] Submitting lightweight scheduled operation: {}",
            script.sbatch_preview
        ),
    );
    let submit_output = run_ssh_command_with_timeout(profile, secret, &clear_outputs, 30).await?;
    let job_id = parse_sbatch_job_id(&submit_output).ok_or_else(|| {
        format!(
            "Failed to parse scheduled utility job ID from sbatch output: {}",
            submit_output.trim()
        )
    })?;
    emit_line(
        app,
        event_name,
        &format!("[QCortado] Scheduled utility job submitted: {}", job_id),
    );

    let started = Instant::now();
    let mut previous_stdout = String::new();
    let mut previous_stderr = String::new();
    let mut last_state = String::new();
    let terminal_state = loop {
        let stdout = read_remote_output(profile, secret, &stdout_path).await;
        let stderr = read_remote_output(profile, secret, &stderr_path).await;
        emit_incremental_output(app, event_name, &stdout, &mut previous_stdout, "");
        emit_incremental_output(app, event_name, &stderr, &mut previous_stderr, "[stderr] ");

        let squeue_cmd = format!("squeue -h -j {} -o \"%T|%N\"", shell_single_quote(&job_id));
        let snapshot = run_ssh_command_with_timeout(profile, secret, &squeue_cmd, 20)
            .await
            .ok()
            .and_then(|output| parse_squeue_snapshot(&output));
        let snapshot = if snapshot.is_some() {
            snapshot
        } else {
            let sacct_cmd = format!(
                "sacct -j {} --format=State,NodeList --parsable2 --noheader",
                shell_single_quote(&job_id)
            );
            run_ssh_command_with_timeout(profile, secret, &sacct_cmd, 20)
                .await
                .ok()
                .and_then(|output| parse_sacct_snapshot(&output))
        };
        if let Some(snapshot) = snapshot {
            let state = normalize_scheduler_state(&snapshot.state);
            if state != last_state {
                emit_line(
                    app,
                    event_name,
                    &format!("[QCortado] Utility job state: {}", state),
                );
                last_state = state.clone();
            }
            if is_terminal_state(&state) {
                break state;
            }
        }
        if started.elapsed() >= Duration::from_secs(timeout_secs.max(1)) {
            let cancel = format!("scancel {}", shell_single_quote(&job_id));
            let _ = run_ssh_command_with_timeout(profile, secret, &cancel, 20).await;
            return Err(
                "Scheduled utility operation timed out while waiting for Slurm.".to_string(),
            );
        }
        tokio::time::sleep(Duration::from_secs(1)).await;
    };

    let stdout = read_remote_output(profile, secret, &stdout_path).await;
    let stderr = read_remote_output(profile, secret, &stderr_path).await;
    emit_incremental_output(app, event_name, &stdout, &mut previous_stdout, "");
    emit_incremental_output(app, event_name, &stderr, &mut previous_stderr, "[stderr] ");
    let output = if stderr.trim().is_empty() {
        stdout
    } else if stdout.trim().is_empty() {
        stderr
    } else {
        format!("{}\n{}", stdout, stderr)
    };
    if terminal_state != "COMPLETED" {
        return Err(format!(
            "Scheduled utility job ended with state {}. {}",
            terminal_state,
            output.trim()
        ));
    }
    Ok(ScheduledUtilityResult {
        job_id,
        output,
        slurm_script: script.script,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hpc::profile::default_utility_resources;

    fn test_profile() -> HpcProfile {
        let mut profile: HpcProfile = serde_json::from_value(serde_json::json!({
            "id": "profile",
            "name": "Cluster",
            "host": "cluster.example.edu",
            "username": "researcher",
            "remote_qe_bin_dir": "/opt/qe/bin",
            "remote_pseudo_dir": "/opt/qe/pseudo",
            "remote_workspace_root": "/scratch/qcortado",
            "remote_project_root": "/project/qcortado"
        }))
        .expect("profile");
        let mut utility = default_utility_resources();
        utility.partition = Some("interactive".to_string());
        utility.account = Some("allocation".to_string());
        profile.utility_resources = Some(utility);
        profile
    }

    #[test]
    fn utility_scripts_allocate_compute_resources_for_engine_commands() {
        let script = build_scheduled_utility_script(
            &test_profile(),
            "qc-w2k-struct",
            &[
                "module load 'WIEN2k/24.1'".to_string(),
                "setrmt_lapw 'case' -r 0".to_string(),
            ],
        )
        .expect("utility script");

        assert!(script.script.contains("#SBATCH --partition=interactive"));
        assert!(script.script.contains("#SBATCH --time=00:10:00"));
        assert!(script.script.contains("#SBATCH --ntasks=1"));
        assert!(script.script.contains("#SBATCH --account=allocation"));
        assert!(script.script.contains("module load 'WIEN2k/24.1'"));
        assert!(script.script.contains("setrmt_lapw 'case' -r 0"));
    }
}
