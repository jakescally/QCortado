//! Project management for QCortado
//!
//! Handles the hierarchy: Project -> CIF Variants -> Calculations
//! Projects are stored in the app data directory.

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;
use tauri::{AppHandle, Manager};

use crate::qe::QEResult;

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
    if let Some(work_dir) = working_dir {
        let work_path = PathBuf::from(&work_dir);
        if work_path.exists() {
            let tmp_dir = calc_dir.join("tmp");
            fs::create_dir_all(&tmp_dir)
                .map_err(|e| format!("Failed to create tmp directory: {}", e))?;

            // Copy all files and directories from working dir
            copy_dir_recursive(&work_path, &tmp_dir)?;
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
