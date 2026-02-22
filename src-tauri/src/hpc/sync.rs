use std::fs;
use std::path::{Path, PathBuf};

fn ensure_safe_rel_path(path: &str) -> Result<PathBuf, String> {
    let rel = PathBuf::from(path);
    if rel.is_absolute() {
        return Err(format!("Refusing to write absolute bundle path: {}", path));
    }
    for component in rel.components() {
        if matches!(component, std::path::Component::ParentDir) {
            return Err(format!("Refusing to write parent-relative bundle path: {}", path));
        }
    }
    Ok(rel)
}

pub fn create_local_bundle_dir(task_id: &str) -> Result<PathBuf, String> {
    let dir = std::env::temp_dir().join(format!("qcortado_hpc_bundle_{}", task_id));
    if dir.exists() {
        fs::remove_dir_all(&dir).map_err(|e| {
            format!(
                "Failed to reset local HPC bundle directory {}: {}",
                dir.display(),
                e
            )
        })?;
    }
    fs::create_dir_all(&dir).map_err(|e| {
        format!(
            "Failed to create local HPC bundle directory {}: {}",
            dir.display(),
            e
        )
    })?;
    Ok(dir)
}

pub fn write_bundle_text_file(bundle_dir: &Path, relative_path: &str, content: &str) -> Result<(), String> {
    let rel = ensure_safe_rel_path(relative_path)?;
    let full = bundle_dir.join(rel);
    if let Some(parent) = full.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
    }
    fs::write(&full, content).map_err(|e| format!("Failed to write {}: {}", full.display(), e))
}

pub fn copy_path_into_bundle(
    bundle_dir: &Path,
    source_path: &Path,
    relative_dest: &str,
) -> Result<(), String> {
    let rel = ensure_safe_rel_path(relative_dest)?;
    let dest_path = bundle_dir.join(rel);
    if source_path.is_dir() {
        copy_dir_recursive(source_path, &dest_path)?;
        return Ok(());
    }
    if source_path.is_file() {
        if let Some(parent) = dest_path.parent() {
            fs::create_dir_all(parent)
                .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
        }
        fs::copy(source_path, &dest_path).map_err(|e| {
            format!(
                "Failed to copy {} -> {}: {}",
                source_path.display(),
                dest_path.display(),
                e
            )
        })?;
        return Ok(());
    }
    Err(format!(
        "Source path does not exist for bundle copy: {}",
        source_path.display()
    ))
}

pub fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<(), String> {
    if !dst.exists() {
        fs::create_dir_all(dst)
            .map_err(|e| format!("Failed to create directory {}: {}", dst.display(), e))?;
    }

    for entry in fs::read_dir(src).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        let dest_path = dst.join(entry.file_name());
        if path.is_dir() {
            copy_dir_recursive(&path, &dest_path)?;
        } else if path.is_file() {
            if let Some(parent) = dest_path.parent() {
                fs::create_dir_all(parent).map_err(|e| {
                    format!("Failed to create directory {}: {}", parent.display(), e)
                })?;
            }
            fs::copy(&path, &dest_path).map_err(|e| {
                format!(
                    "Failed to copy {} -> {}: {}",
                    path.display(),
                    dest_path.display(),
                    e
                )
            })?;
        }
    }

    Ok(())
}
