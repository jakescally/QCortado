//! Project management for QCortado
//!
//! Handles the hierarchy: Project -> CIF Variants -> Calculations
//! Projects are stored in the app data directory.

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;
use tauri::{AppHandle, Manager};

use crate::qe::{QEResult, read_phonon_dos_file, read_phonon_dispersion_file};

// ============================================================================
// Types
// ============================================================================

/// A project containing CIF variants and their calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Project {
    pub id: String,
    pub name: String,
    pub description: Option<String>,
    pub created_at: String,
    pub cif_variants: Vec<CifVariant>,
    /// ID of the last opened CIF variant (for restoring view state)
    #[serde(default)]
    pub last_opened_cif_id: Option<String>,
}

/// A CIF structure variant within a project
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CifVariant {
    pub id: String,
    pub filename: String,
    pub formula: String,
    pub added_at: String,
    pub calculations: Vec<CalculationRun>,
}

/// A calculation run associated with a CIF variant
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalculationRun {
    pub id: String,
    pub calc_type: String,
    pub parameters: serde_json::Value,
    pub result: Option<QEResult>,
    pub started_at: String,
    pub completed_at: Option<String>,
    /// Tags for categorizing calculations (e.g., "phonon-ready", "structure-optimized")
    #[serde(default)]
    pub tags: Vec<String>,
}

/// Summary info for project listing (lighter than full Project)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectSummary {
    pub id: String,
    pub name: String,
    pub description: Option<String>,
    pub created_at: String,
    pub formula: Option<String>,
    pub calculation_count: usize,
    pub last_activity: String,
}

/// Data needed to save a calculation to a project
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SaveCalculationData {
    pub calc_type: String,
    pub parameters: serde_json::Value,
    pub result: QEResult,
    pub started_at: String,
    pub completed_at: String,
    pub input_content: String,
    pub output_content: String,
    /// Tags for categorizing calculations (e.g., "phonon-ready", "structure-optimized")
    #[serde(default)]
    pub tags: Vec<String>,
}

/// Data about a CIF file to add to a project
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CifData {
    pub filename: String,
    pub formula: String,
    pub content: String,
    pub crystal_data: serde_json::Value,
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Gets the projects directory path
pub fn get_projects_dir(app: &AppHandle) -> Result<PathBuf, String> {
    let app_data = app
        .path()
        .app_data_dir()
        .map_err(|e| format!("Failed to get app data dir: {}", e))?;
    Ok(app_data.join("projects"))
}

/// Ensures the projects directory exists
pub fn ensure_projects_dir(app: &AppHandle) -> Result<PathBuf, String> {
    let projects_dir = get_projects_dir(app)?;
    if !projects_dir.exists() {
        fs::create_dir_all(&projects_dir)
            .map_err(|e| format!("Failed to create projects directory: {}", e))?;
    }
    Ok(projects_dir)
}

/// Generates a unique ID based on timestamp and random suffix
fn generate_id() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_millis();
    let random: u32 = rand_simple();
    format!("{}_{:08x}", timestamp, random)
}

/// Simple pseudo-random number generator (no external dependency)
fn rand_simple() -> u32 {
    use std::time::{SystemTime, UNIX_EPOCH};
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .subsec_nanos();
    nanos.wrapping_mul(1103515245).wrapping_add(12345)
}

/// Gets the current timestamp as ISO string
fn now_iso() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let duration = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap();
    let secs = duration.as_secs();
    // Format as ISO 8601 (simplified)
    let days_since_1970 = secs / 86400;
    let time_of_day = secs % 86400;
    let hours = time_of_day / 3600;
    let minutes = (time_of_day % 3600) / 60;
    let seconds = time_of_day % 60;

    // Calculate year, month, day from days since 1970
    let (year, month, day) = days_to_ymd(days_since_1970 as i64);

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, hours, minutes, seconds
    )
}

/// Convert days since 1970 to year/month/day
fn days_to_ymd(days: i64) -> (i32, u32, u32) {
    let mut remaining = days;
    let mut year = 1970i32;

    loop {
        let days_in_year = if is_leap_year(year) { 366 } else { 365 };
        if remaining < days_in_year {
            break;
        }
        remaining -= days_in_year;
        year += 1;
    }

    let days_in_months: [i64; 12] = if is_leap_year(year) {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };

    let mut month = 1u32;
    for &days_in_month in &days_in_months {
        if remaining < days_in_month {
            break;
        }
        remaining -= days_in_month;
        month += 1;
    }

    (year, month, (remaining + 1) as u32)
}

fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

// ============================================================================
// Tauri Commands
// ============================================================================

/// Lists all projects with summary information
#[tauri::command]
pub fn list_projects(app: AppHandle) -> Result<Vec<ProjectSummary>, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let mut summaries = Vec::new();

    let entries = fs::read_dir(&projects_dir)
        .map_err(|e| format!("Failed to read projects directory: {}", e))?;

    for entry in entries {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();

        if !path.is_dir() {
            continue;
        }

        let project_json = path.join("project.json");
        if !project_json.exists() {
            continue;
        }

        let content = fs::read_to_string(&project_json)
            .map_err(|e| format!("Failed to read project.json: {}", e))?;
        let project: Project = serde_json::from_str(&content)
            .map_err(|e| format!("Failed to parse project.json: {}", e))?;

        // Calculate summary info
        let calculation_count: usize = project.cif_variants.iter()
            .map(|v| v.calculations.len())
            .sum();

        let formula = project.cif_variants.first().map(|v| v.formula.clone());

        let last_activity = project.cif_variants.iter()
            .flat_map(|v| v.calculations.iter())
            .filter_map(|c| c.completed_at.as_ref())
            .max()
            .cloned()
            .unwrap_or_else(|| project.created_at.clone());

        summaries.push(ProjectSummary {
            id: project.id,
            name: project.name,
            description: project.description,
            created_at: project.created_at,
            formula,
            calculation_count,
            last_activity,
        });
    }

    // Sort by last activity (most recent first)
    summaries.sort_by(|a, b| b.last_activity.cmp(&a.last_activity));

    Ok(summaries)
}

/// Creates a new project
#[tauri::command]
pub fn create_project(
    app: AppHandle,
    name: String,
    description: Option<String>,
) -> Result<Project, String> {
    let projects_dir = ensure_projects_dir(&app)?;

    let id = generate_id();
    let project_dir = projects_dir.join(&id);

    // Create project directory structure
    fs::create_dir_all(&project_dir)
        .map_err(|e| format!("Failed to create project directory: {}", e))?;
    fs::create_dir_all(project_dir.join("structures"))
        .map_err(|e| format!("Failed to create structures directory: {}", e))?;
    fs::create_dir_all(project_dir.join("calculations"))
        .map_err(|e| format!("Failed to create calculations directory: {}", e))?;

    let project = Project {
        id: id.clone(),
        name,
        description,
        created_at: now_iso(),
        cif_variants: Vec::new(),
        last_opened_cif_id: None,
    };

    // Save project.json
    let project_json = serde_json::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))?;
    fs::write(project_dir.join("project.json"), project_json)
        .map_err(|e| format!("Failed to write project.json: {}", e))?;

    Ok(project)
}

/// Gets a project by ID
#[tauri::command]
pub fn get_project(app: AppHandle, project_id: String) -> Result<Project, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    let project_json = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;

    serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))
}

/// Adds a CIF file to a project
#[tauri::command]
pub fn add_cif_to_project(
    app: AppHandle,
    project_id: String,
    cif_data: CifData,
) -> Result<CifVariant, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    // Load existing project
    let project_json_path = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json_path)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;
    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    // Create CIF variant
    let cif_id = generate_id();
    let structures_dir = project_dir.join("structures");

    // Save CIF file
    fs::write(structures_dir.join(format!("{}.cif", cif_id)), &cif_data.content)
        .map_err(|e| format!("Failed to write CIF file: {}", e))?;

    // Save parsed crystal data JSON
    let crystal_json = serde_json::to_string_pretty(&cif_data.crystal_data)
        .map_err(|e| format!("Failed to serialize crystal data: {}", e))?;
    fs::write(structures_dir.join(format!("{}.json", cif_id)), crystal_json)
        .map_err(|e| format!("Failed to write crystal JSON: {}", e))?;

    let variant = CifVariant {
        id: cif_id,
        filename: cif_data.filename,
        formula: cif_data.formula,
        added_at: now_iso(),
        calculations: Vec::new(),
    };

    // Add to project and save
    project.cif_variants.push(variant.clone());

    let project_json = serde_json::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))?;
    fs::write(&project_json_path, project_json)
        .map_err(|e| format!("Failed to write project.json: {}", e))?;

    Ok(variant)
}

/// Saves a calculation to a project
#[tauri::command]
pub fn save_calculation(
    app: AppHandle,
    project_id: String,
    cif_id: String,
    calc_data: SaveCalculationData,
    working_dir: Option<String>,
) -> Result<CalculationRun, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    // Load existing project
    let project_json_path = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json_path)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;
    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    // Find the CIF variant
    let variant = project.cif_variants.iter_mut()
        .find(|v| v.id == cif_id)
        .ok_or_else(|| format!("CIF variant not found: {}", cif_id))?;

    // Create calculation directory
    let calc_id = generate_id();
    let calc_dir = project_dir.join("calculations").join(&calc_id);
    fs::create_dir_all(&calc_dir)
        .map_err(|e| format!("Failed to create calculation directory: {}", e))?;

    // Save input file
    fs::write(calc_dir.join("pw.in"), &calc_data.input_content)
        .map_err(|e| format!("Failed to write pw.in: {}", e))?;

    // Save output file
    fs::write(calc_dir.join("pw.out"), &calc_data.output_content)
        .map_err(|e| format!("Failed to write pw.out: {}", e))?;

    // Copy working directory contents if provided (includes .save directory, etc.)
    let mut copied_work_path: Option<PathBuf> = None;
    if let Some(work_dir) = working_dir {
        let work_path = PathBuf::from(&work_dir);
        if work_path.exists() {
            copied_work_path = Some(work_path.clone());
            let tmp_dir = calc_dir.join("tmp");
            fs::create_dir_all(&tmp_dir)
                .map_err(|e| format!("Failed to create tmp directory: {}", e))?;

            // Copy all files and directories from working dir
            copy_dir_recursive(&work_path, &tmp_dir)?;
        }
    }

    // For phonons, keep an explicit raw pipeline log and top-level output artifacts
    // so parsing fixes can be applied without depending on intermediate state.
    if calc_data.calc_type == "phonon" {
        fs::write(calc_dir.join("phonon_pipeline.out"), &calc_data.output_content)
            .map_err(|e| format!("Failed to write phonon_pipeline.out: {}", e))?;

        if let Some(work_path) = copied_work_path.as_ref() {
            for filename in ["phonon_freq", "phonon_freq.gp", "phonon_dos"] {
                let source = work_path.join(filename);
                if source.exists() {
                    fs::copy(&source, calc_dir.join(filename))
                        .map_err(|e| format!("Failed to copy {}: {}", filename, e))?;
                }
            }
        }
    }

    // Create calculation run record
    let calc_run = CalculationRun {
        id: calc_id,
        calc_type: calc_data.calc_type,
        parameters: calc_data.parameters,
        result: Some(calc_data.result),
        started_at: calc_data.started_at,
        completed_at: Some(calc_data.completed_at),
        tags: calc_data.tags,
    };

    // Save calc.json
    let calc_json = serde_json::to_string_pretty(&calc_run)
        .map_err(|e| format!("Failed to serialize calculation: {}", e))?;
    fs::write(calc_dir.join("calc.json"), calc_json)
        .map_err(|e| format!("Failed to write calc.json: {}", e))?;

    // Add to variant and save project
    variant.calculations.push(calc_run.clone());

    let project_json = serde_json::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))?;
    fs::write(&project_json_path, project_json)
        .map_err(|e| format!("Failed to write project.json: {}", e))?;

    Ok(calc_run)
}

/// Recursively copies a directory
fn copy_dir_recursive(src: &PathBuf, dst: &PathBuf) -> Result<(), String> {
    if !dst.exists() {
        fs::create_dir_all(dst)
            .map_err(|e| format!("Failed to create directory: {}", e))?;
    }

    for entry in fs::read_dir(src).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        let dest_path = dst.join(entry.file_name());

        if path.is_dir() {
            copy_dir_recursive(&path, &dest_path)?;
        } else {
            fs::copy(&path, &dest_path)
                .map_err(|e| format!("Failed to copy file: {}", e))?;
        }
    }

    Ok(())
}

/// Copies contents of a directory to another directory (public helper for bands calculation)
pub fn copy_dir_contents(src: &PathBuf, dst: &PathBuf) -> Result<(), String> {
    if !dst.exists() {
        fs::create_dir_all(dst)
            .map_err(|e| format!("Failed to create directory: {}", e))?;
    }

    for entry in fs::read_dir(src).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        let dest_path = dst.join(entry.file_name());

        if path.is_dir() {
            copy_dir_recursive(&path, &dest_path)?;
        } else {
            fs::copy(&path, &dest_path)
                .map_err(|e| format!("Failed to copy file: {}", e))?;
        }
    }

    Ok(())
}

/// Deletes a project
#[tauri::command]
pub fn delete_project(app: AppHandle, project_id: String) -> Result<(), String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    fs::remove_dir_all(&project_dir)
        .map_err(|e| format!("Failed to delete project: {}", e))?;

    Ok(())
}

/// Deletes a calculation from a project
#[tauri::command]
pub fn delete_calculation(
    app: AppHandle,
    project_id: String,
    cif_id: String,
    calc_id: String,
) -> Result<(), String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    // Load existing project
    let project_json_path = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json_path)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;
    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    // Find the CIF variant and remove the calculation
    let variant = project.cif_variants.iter_mut()
        .find(|v| v.id == cif_id)
        .ok_or_else(|| format!("CIF variant not found: {}", cif_id))?;

    let calc_index = variant.calculations.iter()
        .position(|c| c.id == calc_id)
        .ok_or_else(|| format!("Calculation not found: {}", calc_id))?;

    variant.calculations.remove(calc_index);

    // Save updated project
    let project_json = serde_json::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))?;
    fs::write(&project_json_path, project_json)
        .map_err(|e| format!("Failed to write project.json: {}", e))?;

    // Delete calculation directory
    let calc_dir = project_dir.join("calculations").join(&calc_id);
    if calc_dir.exists() {
        fs::remove_dir_all(&calc_dir)
            .map_err(|e| format!("Failed to delete calculation files: {}", e))?;
    }

    Ok(())
}

/// Sets or unsets a tag on an existing calculation
#[tauri::command]
pub fn set_calculation_tag(
    app: AppHandle,
    project_id: String,
    cif_id: String,
    calc_id: String,
    tag: String,
    enabled: bool,
) -> Result<(), String> {
    let tag = tag.trim().to_string();
    if tag.is_empty() {
        return Err("Tag cannot be empty".to_string());
    }

    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    // Load existing project
    let project_json_path = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json_path)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;
    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    // Find calculation and update tags
    let variant = project
        .cif_variants
        .iter_mut()
        .find(|v| v.id == cif_id)
        .ok_or_else(|| format!("CIF variant not found: {}", cif_id))?;
    let calculation = variant
        .calculations
        .iter_mut()
        .find(|c| c.id == calc_id)
        .ok_or_else(|| format!("Calculation not found: {}", calc_id))?;

    if enabled {
        if !calculation.tags.iter().any(|existing| existing == &tag) {
            calculation.tags.push(tag.clone());
        }
    } else {
        calculation.tags.retain(|existing| existing != &tag);
    }
    let updated_calculation = calculation.clone();

    // Save updated project
    let project_json = serde_json::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))?;
    fs::write(&project_json_path, project_json)
        .map_err(|e| format!("Failed to write project.json: {}", e))?;

    // Keep calc.json in sync when it exists
    let calc_json_path = project_dir
        .join("calculations")
        .join(&calc_id)
        .join("calc.json");
    if calc_json_path.exists() {
        let calc_json = serde_json::to_string_pretty(&updated_calculation)
            .map_err(|e| format!("Failed to serialize calculation: {}", e))?;
        fs::write(&calc_json_path, calc_json)
            .map_err(|e| format!("Failed to write calc.json: {}", e))?;
    }

    Ok(())
}

/// Sets the last opened CIF variant for a project
#[tauri::command]
pub fn set_last_opened_cif(
    app: AppHandle,
    project_id: String,
    cif_id: String,
) -> Result<(), String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    // Load existing project
    let project_json_path = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json_path)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;
    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    // Verify CIF exists
    if !project.cif_variants.iter().any(|v| v.id == cif_id) {
        return Err(format!("CIF variant not found: {}", cif_id));
    }

    // Update and save
    project.last_opened_cif_id = Some(cif_id);

    let project_json = serde_json::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))?;
    fs::write(&project_json_path, project_json)
        .map_err(|e| format!("Failed to write project.json: {}", e))?;

    Ok(())
}

/// Gets the parsed crystal data for a CIF variant
#[tauri::command]
pub fn get_cif_crystal_data(
    app: AppHandle,
    project_id: String,
    cif_id: String,
) -> Result<serde_json::Value, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    let crystal_json_path = project_dir.join("structures").join(format!("{}.json", cif_id));
    if !crystal_json_path.exists() {
        return Err(format!("Crystal data not found for CIF: {}", cif_id));
    }

    let content = fs::read_to_string(&crystal_json_path)
        .map_err(|e| format!("Failed to read crystal data: {}", e))?;

    serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse crystal data: {}", e))
}

/// Gets the raw CIF content for a CIF variant
#[tauri::command]
pub fn get_cif_content(
    app: AppHandle,
    project_id: String,
    cif_id: String,
) -> Result<String, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    let cif_path = project_dir.join("structures").join(format!("{}.cif", cif_id));
    if !cif_path.exists() {
        return Err(format!("CIF file not found: {}", cif_id));
    }

    fs::read_to_string(&cif_path)
        .map_err(|e| format!("Failed to read CIF file: {}", e))
}

/// Loads saved phonon DOS/dispersion data for a calculation, including recovery
/// from saved output artifacts when embedded data is missing.
#[tauri::command]
pub fn get_saved_phonon_data(
    app: AppHandle,
    project_id: String,
    calc_id: String,
) -> Result<serde_json::Value, String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    // Load existing project
    let project_json_path = project_dir.join("project.json");
    let content = fs::read_to_string(&project_json_path)
        .map_err(|e| format!("Failed to read project.json: {}", e))?;
    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    let mut found_calc = false;
    let mut did_recover = false;
    let mut response: Option<serde_json::Value> = None;

    'outer: for variant in project.cif_variants.iter_mut() {
        for calc in variant.calculations.iter_mut() {
            if calc.id != calc_id {
                continue;
            }
            found_calc = true;

            if calc.calc_type != "phonon" {
                return Err(format!("Calculation {} is not a phonon calculation", calc_id));
            }

            let mut dos_data = calc
                .result
                .as_ref()
                .and_then(|result| result.phonon_data.as_ref())
                .and_then(|phonon| phonon.get("dos_data"))
                .cloned();
            let mut dispersion_data = calc
                .result
                .as_ref()
                .and_then(|result| result.phonon_data.as_ref())
                .and_then(|phonon| phonon.get("dispersion_data"))
                .cloned();

            let tmp_dir = project_dir.join("calculations").join(&calc_id).join("tmp");

            let dos_missing = match dos_data.as_ref() {
                Some(v) => v.is_null(),
                None => true,
            };
            if dos_missing {
                let dos_file = tmp_dir.join("phonon_dos");
                if dos_file.exists() {
                    if let Ok(dos) = read_phonon_dos_file(&dos_file) {
                        dos_data = Some(
                            serde_json::to_value(dos)
                                .map_err(|e| format!("Failed to serialize recovered DOS data: {}", e))?,
                        );
                        did_recover = true;
                    }
                }
            }

            let dispersion_missing = match dispersion_data.as_ref() {
                Some(v) => v.is_null(),
                None => true,
            };
            if dispersion_missing {
                let dispersion_gp_file = tmp_dir.join("phonon_freq.gp");
                let dispersion_file = tmp_dir.join("phonon_freq");
                let source_file = if dispersion_gp_file.exists() {
                    Some(dispersion_gp_file)
                } else if dispersion_file.exists() {
                    Some(dispersion_file)
                } else {
                    None
                };
                if let Some(source_file) = source_file {
                    if let Ok(dispersion) = read_phonon_dispersion_file(&source_file) {
                        dispersion_data = Some(
                            serde_json::to_value(dispersion)
                                .map_err(|e| format!("Failed to serialize recovered dispersion data: {}", e))?,
                        );
                        did_recover = true;
                    }
                }
            }

            let merged = serde_json::json!({
                "dos_data": dos_data.unwrap_or(serde_json::Value::Null),
                "dispersion_data": dispersion_data.unwrap_or(serde_json::Value::Null),
            });

            if did_recover {
                if let Some(result) = calc.result.as_mut() {
                    result.phonon_data = Some(merged.clone());
                }
            }

            response = Some(merged);
            break 'outer;
        }
    }

    if !found_calc {
        return Err(format!("Calculation not found: {}", calc_id));
    }

    if did_recover {
        // Persist recovered phonon data into project and calc JSON
        let project_json = serde_json::to_string_pretty(&project)
            .map_err(|e| format!("Failed to serialize project: {}", e))?;
        fs::write(&project_json_path, project_json)
            .map_err(|e| format!("Failed to write project.json: {}", e))?;

        let calc_json_path = project_dir
            .join("calculations")
            .join(&calc_id)
            .join("calc.json");
        if calc_json_path.exists() {
            let calc = project
                .cif_variants
                .iter()
                .flat_map(|variant| variant.calculations.iter())
                .find(|calc| calc.id == calc_id)
                .ok_or_else(|| format!("Calculation not found while writing calc.json: {}", calc_id))?;
            let calc_json = serde_json::to_string_pretty(calc)
                .map_err(|e| format!("Failed to serialize calculation: {}", e))?;
            fs::write(&calc_json_path, calc_json)
                .map_err(|e| format!("Failed to write calc.json: {}", e))?;
        }
    }

    Ok(response.unwrap_or_else(|| serde_json::json!({
        "dos_data": serde_json::Value::Null,
        "dispersion_data": serde_json::Value::Null,
    })))
}
