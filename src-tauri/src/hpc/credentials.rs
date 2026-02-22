use std::collections::HashMap;
use std::process::Command;
use std::sync::{Mutex, OnceLock};

const KEYCHAIN_SERVICE: &str = "qcortado.hpc";

fn session_credentials() -> &'static Mutex<HashMap<String, String>> {
    static STORE: OnceLock<Mutex<HashMap<String, String>>> = OnceLock::new();
    STORE.get_or_init(|| Mutex::new(HashMap::new()))
}

pub fn set_session_secret(profile_id: &str, secret: String) {
    if let Ok(mut guard) = session_credentials().lock() {
        guard.insert(profile_id.to_string(), secret);
    }
}

pub fn get_session_secret(profile_id: &str) -> Option<String> {
    let Ok(guard) = session_credentials().lock() else {
        return None;
    };
    guard.get(profile_id).cloned()
}

pub fn clear_session_secret(profile_id: &str) {
    if let Ok(mut guard) = session_credentials().lock() {
        guard.remove(profile_id);
    }
}

fn keychain_account(profile_id: &str, username: &str, host: &str) -> String {
    format!("{}:{}@{}", profile_id, username, host)
}

#[cfg(target_os = "macos")]
pub fn save_persisted_secret(
    profile_id: &str,
    username: &str,
    host: &str,
    secret: &str,
) -> Result<(), String> {
    let account = keychain_account(profile_id, username, host);
    let status = Command::new("security")
        .args([
            "add-generic-password",
            "-U",
            "-a",
            &account,
            "-s",
            KEYCHAIN_SERVICE,
            "-w",
            secret,
        ])
        .status()
        .map_err(|e| format!("Failed to invoke macOS keychain command: {}", e))?;

    if status.success() {
        Ok(())
    } else {
        Err("Failed to persist secret in macOS keychain.".to_string())
    }
}

#[cfg(not(target_os = "macos"))]
pub fn save_persisted_secret(
    _profile_id: &str,
    _username: &str,
    _host: &str,
    _secret: &str,
) -> Result<(), String> {
    Err("Secure credential persistence is currently supported on macOS only.".to_string())
}

#[cfg(target_os = "macos")]
pub fn load_persisted_secret(
    profile_id: &str,
    username: &str,
    host: &str,
) -> Result<Option<String>, String> {
    let account = keychain_account(profile_id, username, host);
    let output = Command::new("security")
        .args([
            "find-generic-password",
            "-a",
            &account,
            "-s",
            KEYCHAIN_SERVICE,
            "-w",
        ])
        .output()
        .map_err(|e| format!("Failed to invoke macOS keychain command: {}", e))?;

    if output.status.success() {
        let secret = String::from_utf8_lossy(&output.stdout).trim().to_string();
        if secret.is_empty() {
            Ok(None)
        } else {
            Ok(Some(secret))
        }
    } else {
        Ok(None)
    }
}

#[cfg(not(target_os = "macos"))]
pub fn load_persisted_secret(
    _profile_id: &str,
    _username: &str,
    _host: &str,
) -> Result<Option<String>, String> {
    Ok(None)
}

#[cfg(target_os = "macos")]
pub fn delete_persisted_secret(profile_id: &str, username: &str, host: &str) -> Result<(), String> {
    let account = keychain_account(profile_id, username, host);
    let status = Command::new("security")
        .args([
            "delete-generic-password",
            "-a",
            &account,
            "-s",
            KEYCHAIN_SERVICE,
        ])
        .status()
        .map_err(|e| format!("Failed to invoke macOS keychain command: {}", e))?;

    if status.success() {
        Ok(())
    } else {
        // Treat "not found" as non-fatal because cleanup should be idempotent.
        Ok(())
    }
}

#[cfg(not(target_os = "macos"))]
pub fn delete_persisted_secret(
    _profile_id: &str,
    _username: &str,
    _host: &str,
) -> Result<(), String> {
    Ok(())
}

pub fn resolve_secret(
    profile_id: &str,
    username: &str,
    host: &str,
    persisted: bool,
) -> Result<Option<String>, String> {
    if let Some(secret) = get_session_secret(profile_id) {
        return Ok(Some(secret));
    }

    if persisted {
        return load_persisted_secret(profile_id, username, host);
    }

    Ok(None)
}
