use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use tauri::{AppHandle, Manager};

use crate::hpc::profile::HpcProfile;
use crate::hpc::ssh::{download_directory_contents, run_ssh_command, upload_directory};
use crate::projects;

pub const VIEWER_LIBRARY_SCHEMA: &str = "viewer_library_manifest.v1";
const VIEWER_LIBRARY_DIR_NAME: &str = "_qcortado_viewer_library";
const LOCAL_SYNC_MANIFEST_FILE: &str = "viewer_sync_manifest.json";
const PROJECT_FOLDERS_FILE_NAME: &str = "folders.json";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewerManifestProject {
    pub project_id: String,
    pub revision: String,
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub formula: Option<String>,
    pub calculation_count: usize,
    pub last_activity: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewerLibraryManifest {
    pub schema_version: String,
    pub generated_at: String,
    pub projects: Vec<ViewerManifestProject>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewerPublishResult {
    pub profile_id: String,
    pub remote_library_root: String,
    pub generated_at: String,
    pub published_projects: usize,
    pub requested_project_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewerSyncResult {
    pub profile_id: String,
    pub remote_library_root: String,
    pub remote_generated_at: String,
    pub synced_at: String,
    pub downloaded_projects: usize,
    pub removed_projects: usize,
    pub skipped_projects: usize,
    pub total_projects: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ViewerSyncStatus {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub last_synced_at: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub last_error: Option<String>,
    #[serde(default)]
    pub local_project_count: usize,
}

fn now_iso() -> String {
    chrono::Utc::now().to_rfc3339()
}

fn shell_single_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    let escaped = value.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

async fn resolve_remote_path(
    profile: &HpcProfile,
    secret: Option<&str>,
    raw_path: &str,
) -> Result<String, String> {
    let trimmed = raw_path.trim();
    if trimmed.is_empty() {
        return Err("Remote path is empty".to_string());
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

pub async fn resolve_remote_project_root(
    profile: &HpcProfile,
    secret: Option<&str>,
) -> Result<String, String> {
    resolve_remote_path(
        profile,
        secret,
        profile.remote_project_root.trim_end_matches('/'),
    )
    .await
}

pub fn build_remote_library_root(project_root: &str) -> String {
    format!(
        "{}/{}",
        project_root.trim_end_matches('/'),
        VIEWER_LIBRARY_DIR_NAME
    )
}

fn revision_for_summary(summary: &projects::ProjectSummary) -> String {
    format!(
        "{}:{}:{}",
        summary.last_activity.trim(),
        summary.calculation_count,
        summary.id
    )
}

fn copy_file_if_exists(src: &Path, dst: &Path) -> Result<(), String> {
    if !src.exists() || !src.is_file() {
        return Ok(());
    }
    if let Some(parent) = dst.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
    }
    fs::copy(src, dst).map_err(|e| {
        format!(
            "Failed to copy {} -> {}: {}",
            src.display(),
            dst.display(),
            e
        )
    })?;
    Ok(())
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<(), String> {
    if !src.exists() {
        return Ok(());
    }
    fs::create_dir_all(dst)
        .map_err(|e| format!("Failed to create directory {}: {}", dst.display(), e))?;

    for entry in fs::read_dir(src)
        .map_err(|e| format!("Failed to read directory {}: {}", src.display(), e))?
    {
        let entry = entry.map_err(|e| e.to_string())?;
        let src_path = entry.path();
        let dst_path = dst.join(entry.file_name());
        if src_path.is_dir() {
            copy_dir_recursive(&src_path, &dst_path)?;
        } else if src_path.is_file() {
            if let Some(parent) = dst_path.parent() {
                fs::create_dir_all(parent).map_err(|e| {
                    format!("Failed to create directory {}: {}", parent.display(), e)
                })?;
            }
            fs::copy(&src_path, &dst_path).map_err(|e| {
                format!(
                    "Failed to copy {} -> {}: {}",
                    src_path.display(),
                    dst_path.display(),
                    e
                )
            })?;
        }
    }

    Ok(())
}

fn local_sync_manifest_path(app: &AppHandle) -> Result<PathBuf, String> {
    let app_data = app
        .path()
        .app_data_dir()
        .map_err(|e| format!("Failed to get app data dir: {}", e))?;
    fs::create_dir_all(&app_data).map_err(|e| {
        format!(
            "Failed to create app data directory {}: {}",
            app_data.display(),
            e
        )
    })?;
    Ok(app_data.join(LOCAL_SYNC_MANIFEST_FILE))
}

pub fn load_local_sync_manifest(app: &AppHandle) -> Result<Option<ViewerLibraryManifest>, String> {
    let path = local_sync_manifest_path(app)?;
    if !path.exists() {
        return Ok(None);
    }
    let content = fs::read_to_string(&path).map_err(|e| {
        format!(
            "Failed to read local sync manifest {}: {}",
            path.display(),
            e
        )
    })?;
    let manifest: ViewerLibraryManifest = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse local sync manifest: {}", e))?;
    Ok(Some(manifest))
}

fn save_local_sync_manifest(
    app: &AppHandle,
    manifest: &ViewerLibraryManifest,
) -> Result<(), String> {
    let path = local_sync_manifest_path(app)?;
    let serialized = serde_json::to_string_pretty(manifest)
        .map_err(|e| format!("Failed to serialize manifest: {}", e))?;
    fs::write(&path, serialized).map_err(|e| {
        format!(
            "Failed to write local sync manifest {}: {}",
            path.display(),
            e
        )
    })
}

fn build_snapshot_for_project(
    projects_dir: &Path,
    snapshot_projects_root: &Path,
    summary: &projects::ProjectSummary,
) -> Result<ViewerManifestProject, String> {
    let src_project_dir = projects_dir.join(&summary.id);
    let src_structures_dir = src_project_dir.join("structures");
    let src_calculations_dir = src_project_dir.join("calculations");

    let dst_project_dir = snapshot_projects_root.join(&summary.id);
    let dst_structures_dir = dst_project_dir.join("structures");
    let dst_calculations_dir = dst_project_dir.join("calculations");

    fs::create_dir_all(&dst_project_dir).map_err(|e| {
        format!(
            "Failed to create snapshot directory {}: {}",
            dst_project_dir.display(),
            e
        )
    })?;

    copy_file_if_exists(
        &src_project_dir.join("project.json"),
        &dst_project_dir.join("project.json"),
    )?;

    if src_structures_dir.exists() {
        fs::create_dir_all(&dst_structures_dir).map_err(|e| {
            format!(
                "Failed to create snapshot structures directory {}: {}",
                dst_structures_dir.display(),
                e
            )
        })?;
        for entry in fs::read_dir(&src_structures_dir).map_err(|e| {
            format!(
                "Failed to read structures directory {}: {}",
                src_structures_dir.display(),
                e
            )
        })? {
            let entry = entry.map_err(|e| e.to_string())?;
            let src_path = entry.path();
            if !src_path.is_file() {
                continue;
            }
            let Some(ext) = src_path.extension().and_then(|value| value.to_str()) else {
                continue;
            };
            if ext != "cif" && ext != "json" {
                continue;
            }
            let dst_path = dst_structures_dir.join(entry.file_name());
            copy_file_if_exists(&src_path, &dst_path)?;
        }
    }

    if src_calculations_dir.exists() {
        fs::create_dir_all(&dst_calculations_dir).map_err(|e| {
            format!(
                "Failed to create snapshot calculations directory {}: {}",
                dst_calculations_dir.display(),
                e
            )
        })?;
        for entry in fs::read_dir(&src_calculations_dir).map_err(|e| {
            format!(
                "Failed to read calculations directory {}: {}",
                src_calculations_dir.display(),
                e
            )
        })? {
            let entry = entry.map_err(|e| e.to_string())?;
            let calc_dir = entry.path();
            if !calc_dir.is_dir() {
                continue;
            }
            let calc_id = entry.file_name();
            let src_calc_json = calc_dir.join("calc.json");
            let dst_calc_json = dst_calculations_dir.join(calc_id).join("calc.json");
            copy_file_if_exists(&src_calc_json, &dst_calc_json)?;
        }
    }

    Ok(ViewerManifestProject {
        project_id: summary.id.clone(),
        revision: revision_for_summary(summary),
        name: summary.name.clone(),
        description: summary.description.clone(),
        formula: summary.formula.clone(),
        calculation_count: summary.calculation_count,
        last_activity: summary.last_activity.clone(),
    })
}

fn build_local_snapshot(
    app: &AppHandle,
    stage_dir: &Path,
    requested_project_id: Option<&str>,
) -> Result<ViewerLibraryManifest, String> {
    let projects_dir = projects::ensure_projects_dir(app)?;
    let all_summaries = projects::list_projects(app.clone())?;

    let mut summaries = all_summaries;
    summaries.sort_by(|a, b| a.id.cmp(&b.id));

    let requested_project_id = requested_project_id
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string);

    if let Some(project_id) = requested_project_id.as_deref() {
        let exists = summaries.iter().any(|summary| summary.id == project_id);
        if !exists {
            return Err(format!(
                "Project not found for viewer publish: {}",
                project_id
            ));
        }
    }

    let projects_root = stage_dir.join("projects");
    fs::create_dir_all(&projects_root).map_err(|e| {
        format!(
            "Failed to create snapshot root {}: {}",
            projects_root.display(),
            e
        )
    })?;

    // Preserve project-folder organization metadata for viewer browsing.
    copy_file_if_exists(
        &projects_dir.join(PROJECT_FOLDERS_FILE_NAME),
        &projects_root.join(PROJECT_FOLDERS_FILE_NAME),
    )?;

    let mut manifest_projects: Vec<ViewerManifestProject> = Vec::new();
    for summary in &summaries {
        if let Some(project_id) = requested_project_id.as_deref() {
            if summary.id != project_id {
                continue;
            }
        }
        let entry = build_snapshot_for_project(&projects_dir, &projects_root, summary)?;
        manifest_projects.push(entry);
    }

    let manifest = ViewerLibraryManifest {
        schema_version: VIEWER_LIBRARY_SCHEMA.to_string(),
        generated_at: now_iso(),
        projects: manifest_projects,
    };

    let manifest_path = stage_dir.join("manifest.json");
    let serialized = serde_json::to_string_pretty(&manifest)
        .map_err(|e| format!("Failed to serialize viewer manifest: {}", e))?;
    fs::write(&manifest_path, serialized).map_err(|e| {
        format!(
            "Failed to write viewer manifest {}: {}",
            manifest_path.display(),
            e
        )
    })?;

    Ok(manifest)
}

pub async fn publish_viewer_library(
    app: &AppHandle,
    profile: &HpcProfile,
    secret: Option<&str>,
    requested_project_id: Option<&str>,
) -> Result<ViewerPublishResult, String> {
    let remote_project_root = resolve_remote_project_root(profile, secret).await?;
    let remote_library_root = build_remote_library_root(&remote_project_root);

    let stage_name = format!(".qcortado_viewer_stage_{}", uuid::Uuid::new_v4());
    let temp_root =
        std::env::temp_dir().join(format!("qcortado_viewer_publish_{}", uuid::Uuid::new_v4()));
    let local_stage_dir = temp_root.join(&stage_name);
    fs::create_dir_all(&local_stage_dir).map_err(|e| {
        format!(
            "Failed to create local viewer publish stage {}: {}",
            local_stage_dir.display(),
            e
        )
    })?;

    let manifest = build_local_snapshot(app, &local_stage_dir, requested_project_id)?;

    let remote_stage = format!(
        "{}/{}",
        remote_project_root.trim_end_matches('/'),
        stage_name
    );
    let remote_prev = format!("{}.prev", remote_library_root);

    let prepare_cmd = format!(
        "mkdir -p {root} && rm -rf {stage}",
        root = shell_single_quote(&remote_project_root),
        stage = shell_single_quote(&remote_stage),
    );
    run_ssh_command(profile, secret, &prepare_cmd).await?;

    upload_directory(profile, secret, &local_stage_dir, &remote_project_root).await?;

    let promote_cmd = format!(
        "set -e; rm -rf {prev}; if [ -d {final_root} ]; then mv {final_root} {prev}; fi; mv {stage} {final_root}; rm -rf {prev}",
        prev = shell_single_quote(&remote_prev),
        final_root = shell_single_quote(&remote_library_root),
        stage = shell_single_quote(&remote_stage),
    );
    run_ssh_command(profile, secret, &promote_cmd).await?;

    let _ = fs::remove_dir_all(&temp_root);

    Ok(ViewerPublishResult {
        profile_id: profile.id.clone(),
        remote_library_root,
        generated_at: manifest.generated_at,
        published_projects: manifest.projects.len(),
        requested_project_id: requested_project_id.map(|value| value.to_string()),
    })
}

fn load_remote_manifest(local_sync_root: &Path) -> Result<ViewerLibraryManifest, String> {
    let manifest_path = local_sync_root.join("manifest.json");
    if !manifest_path.exists() {
        return Err("Viewer library sync failed: manifest.json is missing.".to_string());
    }
    let content = fs::read_to_string(&manifest_path).map_err(|e| {
        format!(
            "Failed to read synced manifest {}: {}",
            manifest_path.display(),
            e
        )
    })?;
    let manifest: ViewerLibraryManifest = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse synced manifest: {}", e))?;
    if manifest.schema_version != VIEWER_LIBRARY_SCHEMA {
        return Err(format!(
            "Unsupported viewer manifest schema: {}",
            manifest.schema_version
        ));
    }
    Ok(manifest)
}

pub async fn sync_viewer_library(
    app: &AppHandle,
    profile: &HpcProfile,
    secret: Option<&str>,
) -> Result<ViewerSyncResult, String> {
    let remote_project_root = resolve_remote_project_root(profile, secret).await?;
    let remote_library_root = build_remote_library_root(&remote_project_root);

    let local_sync_root =
        std::env::temp_dir().join(format!("qcortado_viewer_sync_{}", uuid::Uuid::new_v4()));
    fs::create_dir_all(&local_sync_root).map_err(|e| {
        format!(
            "Failed to create local sync staging directory {}: {}",
            local_sync_root.display(),
            e
        )
    })?;

    download_directory_contents(profile, secret, &remote_library_root, &local_sync_root).await?;
    let remote_manifest = load_remote_manifest(&local_sync_root)?;

    let previous_manifest = load_local_sync_manifest(app)?;
    let mut previous_revisions: HashMap<String, String> = HashMap::new();
    if let Some(manifest) = previous_manifest {
        for project in manifest.projects {
            previous_revisions.insert(project.project_id, project.revision);
        }
    }

    let projects_dir = projects::ensure_projects_dir(app)?;
    let synced_projects_root = local_sync_root.join("projects");
    let mut remote_project_ids: HashSet<String> = HashSet::new();

    let mut downloaded_projects = 0usize;
    let mut skipped_projects = 0usize;

    for project in &remote_manifest.projects {
        remote_project_ids.insert(project.project_id.clone());

        let dest_project_dir = projects_dir.join(&project.project_id);
        let src_project_dir = synced_projects_root.join(&project.project_id);
        if !src_project_dir.exists() {
            return Err(format!(
                "Viewer manifest references missing snapshot for project {}",
                project.project_id
            ));
        }

        let unchanged = previous_revisions
            .get(&project.project_id)
            .map(|revision| revision == &project.revision)
            .unwrap_or(false)
            && dest_project_dir.exists();

        if unchanged {
            skipped_projects += 1;
            continue;
        }

        if dest_project_dir.exists() {
            fs::remove_dir_all(&dest_project_dir).map_err(|e| {
                format!(
                    "Failed to replace local synced project {}: {}",
                    dest_project_dir.display(),
                    e
                )
            })?;
        }

        copy_dir_recursive(&src_project_dir, &dest_project_dir)?;
        downloaded_projects += 1;
    }

    let synced_folders_metadata = synced_projects_root.join(PROJECT_FOLDERS_FILE_NAME);
    let local_folders_metadata = projects_dir.join(PROJECT_FOLDERS_FILE_NAME);
    if synced_folders_metadata.exists() {
        copy_file_if_exists(&synced_folders_metadata, &local_folders_metadata)?;
    } else if local_folders_metadata.exists() {
        fs::remove_file(&local_folders_metadata).map_err(|e| {
            format!(
                "Failed to remove stale local folders metadata {}: {}",
                local_folders_metadata.display(),
                e
            )
        })?;
    }

    let mut removed_projects = 0usize;
    for entry in fs::read_dir(&projects_dir).map_err(|e| {
        format!(
            "Failed to read local projects directory {}: {}",
            projects_dir.display(),
            e
        )
    })? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }
        let project_id = entry.file_name().to_string_lossy().to_string();
        if remote_project_ids.contains(&project_id) {
            continue;
        }
        let project_json = path.join("project.json");
        if !project_json.exists() {
            continue;
        }
        fs::remove_dir_all(&path).map_err(|e| {
            format!(
                "Failed to remove stale local project {}: {}",
                path.display(),
                e
            )
        })?;
        removed_projects += 1;
    }

    save_local_sync_manifest(app, &remote_manifest)?;

    let _ = fs::remove_dir_all(&local_sync_root);

    Ok(ViewerSyncResult {
        profile_id: profile.id.clone(),
        remote_library_root,
        remote_generated_at: remote_manifest.generated_at,
        synced_at: now_iso(),
        downloaded_projects,
        removed_projects,
        skipped_projects,
        total_projects: remote_manifest.projects.len(),
    })
}
