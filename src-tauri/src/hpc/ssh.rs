use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::time::Duration;

use tauri::{AppHandle, Emitter};
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::process::Command;

use super::profile::{HpcAuthMethod, HpcProfile};

fn ensure_parent_dir(path: &Path) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
    }
    Ok(())
}

fn shell_single_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    let escaped = value.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

#[cfg(unix)]
fn make_executable(path: &Path) -> Result<(), String> {
    use std::os::unix::fs::PermissionsExt;
    let mut perms = fs::metadata(path)
        .map_err(|e| format!("Failed to inspect {}: {}", path.display(), e))?
        .permissions();
    perms.set_mode(0o700);
    fs::set_permissions(path, perms)
        .map_err(|e| format!("Failed to set permissions for {}: {}", path.display(), e))
}

#[cfg(not(unix))]
fn make_executable(_path: &Path) -> Result<(), String> {
    Ok(())
}

fn create_askpass_script(secret: &str) -> Result<PathBuf, String> {
    let script_path = std::env::temp_dir().join(format!(
        "qcortado_ssh_askpass_{}_{}.sh",
        std::process::id(),
        uuid::Uuid::new_v4()
    ));
    ensure_parent_dir(&script_path)?;
    let body = format!("#!/bin/sh\nprintf %s {}\n", shell_single_quote(secret));
    fs::write(&script_path, body).map_err(|e| {
        format!(
            "Failed to write askpass script {}: {}",
            script_path.display(),
            e
        )
    })?;
    make_executable(&script_path)?;
    Ok(script_path)
}

fn destination(profile: &HpcProfile) -> String {
    format!("{}@{}", profile.username.trim(), profile.host.trim())
}

fn transport_options(
    profile: &HpcProfile,
    include_batch_mode: bool,
    port_flag: &str,
) -> Vec<String> {
    let mut args = vec![
        port_flag.to_string(),
        profile.port.to_string(),
        "-o".to_string(),
        "StrictHostKeyChecking=yes".to_string(),
        "-o".to_string(),
        "ServerAliveInterval=20".to_string(),
        "-o".to_string(),
        "ServerAliveCountMax=3".to_string(),
        "-o".to_string(),
        "ConnectTimeout=15".to_string(),
    ];
    if include_batch_mode {
        args.push("-o".to_string());
        args.push("BatchMode=yes".to_string());
    }
    if let Some(path) = profile.ssh_key_path.as_ref() {
        if !path.trim().is_empty() {
            args.push("-i".to_string());
            args.push(path.trim().to_string());
        }
    }
    args
}

fn ssh_options(profile: &HpcProfile, include_batch_mode: bool) -> Vec<String> {
    transport_options(profile, include_batch_mode, "-p")
}

fn scp_options(profile: &HpcProfile, include_batch_mode: bool) -> Vec<String> {
    transport_options(profile, include_batch_mode, "-P")
}

async fn run_with_optional_askpass(
    program: &str,
    args: &[String],
    secret: Option<&str>,
    timeout_secs: u64,
) -> Result<std::process::Output, String> {
    let mut cmd = Command::new(program);
    cmd.args(args);
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::piped());
    cmd.stderr(Stdio::piped());

    let askpass_path = if let Some(value) = secret {
        let script_path = create_askpass_script(value)?;
        cmd.env("SSH_ASKPASS", &script_path);
        cmd.env("SSH_ASKPASS_REQUIRE", "force");
        cmd.env("DISPLAY", ":0");
        Some(script_path)
    } else {
        None
    };

    let output_result =
        tokio::time::timeout(Duration::from_secs(timeout_secs.max(1)), cmd.output())
            .await
            .map_err(|_| format!("{} command timed out", program))?;
    let output = output_result.map_err(|e| format!("Failed to execute {}: {}", program, e));

    if let Some(path) = askpass_path {
        let _ = fs::remove_file(path);
    }

    output
}

pub async fn run_ssh_command(
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_command: &str,
) -> Result<String, String> {
    run_ssh_command_with_timeout(profile, secret, remote_command, 120).await
}

pub async fn run_ssh_command_with_timeout(
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_command: &str,
    timeout_secs: u64,
) -> Result<String, String> {
    let use_password = matches!(profile.auth_method, HpcAuthMethod::Password);
    let mut args = ssh_options(profile, !use_password);
    if use_password {
        args.push("-o".to_string());
        args.push("PreferredAuthentications=password,keyboard-interactive".to_string());
        args.push("-o".to_string());
        args.push("PubkeyAuthentication=no".to_string());
    }
    args.push(destination(profile));
    args.push(remote_command.to_string());

    let output = run_with_optional_askpass("ssh", &args, secret, timeout_secs).await?;
    if output.status.success() {
        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        Err(format!("SSH command failed: {}", stderr.trim()))
    }
}

/// Executes a short remote command while emitting stdout/stderr lines to an
/// inline workflow event channel. This does not register a calculation task.
pub async fn run_ssh_command_streaming(
    app: &AppHandle,
    event_name: &str,
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_command: &str,
    timeout_secs: u64,
) -> Result<String, String> {
    let use_password = matches!(profile.auth_method, HpcAuthMethod::Password);
    let mut args = ssh_options(profile, !use_password);
    if use_password {
        args.push("-o".to_string());
        args.push("PreferredAuthentications=password,keyboard-interactive".to_string());
        args.push("-o".to_string());
        args.push("PubkeyAuthentication=no".to_string());
    }
    args.push(destination(profile));
    args.push(remote_command.to_string());

    let mut command = Command::new("ssh");
    command.args(&args);
    command.stdin(Stdio::null());
    command.stdout(Stdio::piped());
    command.stderr(Stdio::piped());
    let askpass_path = if let Some(value) = secret {
        let script_path = create_askpass_script(value)?;
        command.env("SSH_ASKPASS", &script_path);
        command.env("SSH_ASKPASS_REQUIRE", "force");
        command.env("DISPLAY", ":0");
        Some(script_path)
    } else {
        None
    };
    let mut child = command
        .spawn()
        .map_err(|err| format!("Failed to execute ssh: {}", err))?;
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| "Failed to capture SSH stdout.".to_string())?;
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| "Failed to capture SSH stderr.".to_string())?;
    let stdout_app = app.clone();
    let stdout_event = event_name.to_string();
    let stderr_app = app.clone();
    let stderr_event = event_name.to_string();
    let stdout_reader = tokio::spawn(async move {
        let mut lines = BufReader::new(stdout).lines();
        let mut collected = Vec::new();
        while let Ok(Some(line)) = lines.next_line().await {
            let _ = stdout_app.emit(&stdout_event, &line);
            collected.push(line);
        }
        collected
    });
    let stderr_reader = tokio::spawn(async move {
        let mut lines = BufReader::new(stderr).lines();
        let mut collected = Vec::new();
        while let Ok(Some(line)) = lines.next_line().await {
            let display = format!("[stderr] {}", line);
            let _ = stderr_app.emit(&stderr_event, &display);
            collected.push(display);
        }
        collected
    });
    let status =
        match tokio::time::timeout(Duration::from_secs(timeout_secs.max(1)), child.wait()).await {
            Ok(result) => result.map_err(|err| format!("Failed to wait for ssh: {}", err))?,
            Err(_) => {
                let _ = child.kill().await;
                if let Some(path) = askpass_path.as_ref() {
                    let _ = fs::remove_file(path);
                }
                return Err("ssh command timed out".to_string());
            }
        };
    let mut lines = stdout_reader
        .await
        .map_err(|err| format!("Failed to read SSH stdout: {}", err))?;
    lines.extend(
        stderr_reader
            .await
            .map_err(|err| format!("Failed to read SSH stderr: {}", err))?,
    );
    if let Some(path) = askpass_path {
        let _ = fs::remove_file(path);
    }
    let output = lines.join("\n");
    if status.success() {
        Ok(output)
    } else {
        Err(if output.trim().is_empty() {
            "SSH command failed.".to_string()
        } else {
            format!("SSH command failed: {}", output)
        })
    }
}

pub async fn upload_directory(
    profile: &HpcProfile,
    secret: Option<&str>,
    local_dir: &Path,
    remote_parent: &str,
) -> Result<(), String> {
    if !local_dir.exists() || !local_dir.is_dir() {
        return Err(format!(
            "Upload path does not exist or is not a directory: {}",
            local_dir.display()
        ));
    }

    let use_password = matches!(profile.auth_method, HpcAuthMethod::Password);
    let mut args = scp_options(profile, !use_password);
    args.push("-r".to_string());
    args.push(local_dir.display().to_string());
    args.push(format!("{}:{}", destination(profile), remote_parent));
    let output = run_with_optional_askpass("scp", &args, secret, 120).await?;
    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        Err(format!("Failed to upload directory: {}", stderr.trim()))
    }
}

pub async fn upload_file(
    profile: &HpcProfile,
    secret: Option<&str>,
    local_file: &Path,
    remote_file: &str,
) -> Result<(), String> {
    if !local_file.exists() || !local_file.is_file() {
        return Err(format!(
            "Upload path does not exist or is not a file: {}",
            local_file.display()
        ));
    }

    let use_password = matches!(profile.auth_method, HpcAuthMethod::Password);
    let mut args = scp_options(profile, !use_password);
    args.push(local_file.display().to_string());
    args.push(format!("{}:{}", destination(profile), remote_file));
    let output = run_with_optional_askpass("scp", &args, secret, 120).await?;
    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        Err(format!("Failed to upload file: {}", stderr.trim()))
    }
}

pub async fn download_file(
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_file: &str,
    local_file: &Path,
) -> Result<(), String> {
    ensure_parent_dir(local_file)?;

    let use_password = matches!(profile.auth_method, HpcAuthMethod::Password);
    let mut args = scp_options(profile, !use_password);
    args.push(format!("{}:{}", destination(profile), remote_file));
    args.push(local_file.display().to_string());
    let output = run_with_optional_askpass("scp", &args, secret, 120).await?;
    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        Err(format!("Failed to download file: {}", stderr.trim()))
    }
}

pub async fn download_directory_contents(
    profile: &HpcProfile,
    secret: Option<&str>,
    remote_dir: &str,
    local_dir: &Path,
) -> Result<(), String> {
    fs::create_dir_all(local_dir).map_err(|e| {
        format!(
            "Failed to create local sync directory {}: {}",
            local_dir.display(),
            e
        )
    })?;

    let use_password = matches!(profile.auth_method, HpcAuthMethod::Password);
    let mut args = scp_options(profile, !use_password);
    args.push("-r".to_string());
    args.push(format!(
        "{}:{}/.",
        destination(profile),
        remote_dir.trim_end_matches('/')
    ));
    args.push(local_dir.display().to_string());
    let output = run_with_optional_askpass("scp", &args, secret, 120).await?;
    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        Err(format!(
            "Failed to download directory {}: {}",
            remote_dir,
            stderr.trim()
        ))
    }
}
