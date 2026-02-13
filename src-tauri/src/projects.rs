//! Project management for QCortado
//!
//! Handles the hierarchy: Project -> CIF Variants -> Calculations
//! Projects are stored in the app data directory.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::VecDeque;
use std::fs;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant};
use tauri::{AppHandle, Emitter, Manager};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::qe::{read_phonon_dispersion_file, read_phonon_dos_file, QEResult};

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
    /// Total on-disk bytes used by this saved calculation directory.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub storage_bytes: Option<u64>,
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

/// Metadata returned after exporting a project archive.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectArchiveExportResult {
    pub archive_path: String,
    pub project_id: String,
    pub project_name: String,
    pub archive_size_bytes: u64,
}

/// Progress payload emitted during archive export.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectArchiveExportProgress {
    pub export_id: String,
    pub project_id: String,
    pub phase: String,
    pub processed_bytes: u64,
    pub total_bytes: u64,
    pub progress_percent: f64,
}

/// Progress payload emitted during archive import.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectArchiveImportProgress {
    pub import_id: String,
    pub phase: String,
    pub processed_bytes: u64,
    pub total_bytes: u64,
    pub progress_percent: f64,
}

/// Metadata returned after importing a project archive.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectArchiveImportResult {
    pub project_id: String,
    pub project_name: String,
    pub imported_with_new_id: bool,
}

const PROJECT_ARCHIVE_MAGIC: &[u8; 7] = b"QCPROJ1";
const PROJECT_ARCHIVE_VERSION: u32 = 2;
const ARCHIVE_ENTRY_DIR: u8 = 0;
const ARCHIVE_ENTRY_FILE: u8 = 1;
const ARCHIVE_ENTRY_END: u8 = 255;
const ARCHIVE_FILE_COMPRESSION_NONE: u8 = 0;
const ARCHIVE_FILE_COMPRESSION_GZIP: u8 = 1;
const EXPORT_PROGRESS_EMIT_INTERVAL_MS: u64 = 200;
const EXPORT_PROGRESS_EMIT_BYTES: u64 = 8 * 1024 * 1024;
const IMPORT_PROGRESS_EMIT_INTERVAL_MS: u64 = 200;
const IMPORT_PROGRESS_EMIT_BYTES: u64 = 8 * 1024 * 1024;
const EXPORT_COPY_BUFFER_SIZE: usize = 256 * 1024;
const EXPORT_CANCELLED_SENTINEL: &str = "__QCORTADO_EXPORT_CANCELLED__";
const GZIP_MAGIC_PREFIX: [u8; 2] = [0x1F, 0x8B];

static EXPORT_CANCEL_FLAGS: OnceLock<Mutex<HashMap<String, Arc<AtomicBool>>>> = OnceLock::new();

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

fn export_cancel_flags() -> &'static Mutex<HashMap<String, Arc<AtomicBool>>> {
    EXPORT_CANCEL_FLAGS.get_or_init(|| Mutex::new(HashMap::new()))
}

fn register_project_export_cancel_flag(export_id: &str) -> Arc<AtomicBool> {
    let flag = Arc::new(AtomicBool::new(false));
    if let Ok(mut guard) = export_cancel_flags().lock() {
        guard.insert(export_id.to_string(), flag.clone());
    }
    flag
}

fn unregister_project_export_cancel_flag(export_id: &str) {
    if let Ok(mut guard) = export_cancel_flags().lock() {
        guard.remove(export_id);
    }
}

fn request_project_export_cancel(export_id: &str) -> bool {
    if let Ok(guard) = export_cancel_flags().lock() {
        if let Some(flag) = guard.get(export_id) {
            flag.store(true, Ordering::Relaxed);
            return true;
        }
    }
    false
}

fn export_compression_threads() -> u32 {
    let available = std::thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    let target = ((available * 8) + 9) / 10;
    target.max(1).min(available) as u32
}

fn emit_project_export_progress(
    app: &AppHandle,
    export_id: &str,
    project_id: &str,
    phase: &str,
    processed_bytes: u64,
    total_bytes: u64,
) {
    let progress_percent = if total_bytes == 0 {
        if phase == "done" {
            100.0
        } else {
            0.0
        }
    } else {
        ((processed_bytes as f64 / total_bytes as f64) * 100.0).clamp(0.0, 100.0)
    };

    let payload = ProjectArchiveExportProgress {
        export_id: export_id.to_string(),
        project_id: project_id.to_string(),
        phase: phase.to_string(),
        processed_bytes,
        total_bytes,
        progress_percent,
    };

    let _ = app.emit("project-export-progress", payload);
}

fn emit_project_import_progress(
    app: &AppHandle,
    import_id: &str,
    phase: &str,
    processed_bytes: u64,
    total_bytes: u64,
) {
    let progress_percent = if total_bytes == 0 {
        if phase == "done" {
            100.0
        } else {
            0.0
        }
    } else {
        ((processed_bytes as f64 / total_bytes as f64) * 100.0).clamp(0.0, 100.0)
    };

    let payload = ProjectArchiveImportProgress {
        import_id: import_id.to_string(),
        phase: phase.to_string(),
        processed_bytes,
        total_bytes,
        progress_percent,
    };

    let _ = app.emit("project-import-progress", payload);
}

fn write_u32_le<W: Write>(writer: &mut W, value: u32) -> Result<(), String> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|e| format!("Failed to write archive u32: {}", e))
}

fn write_u64_le<W: Write>(writer: &mut W, value: u64) -> Result<(), String> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|e| format!("Failed to write archive u64: {}", e))
}

fn read_u32_le<R: Read>(reader: &mut R) -> Result<u32, String> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|e| format!("Failed to read archive u32: {}", e))?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64, String> {
    let mut bytes = [0_u8; 8];
    reader
        .read_exact(&mut bytes)
        .map_err(|e| format!("Failed to read archive u64: {}", e))?;
    Ok(u64::from_le_bytes(bytes))
}

fn relative_path_to_archive_key(path: &Path) -> Result<String, String> {
    use std::path::Component;

    let mut parts: Vec<&str> = Vec::new();
    for component in path.components() {
        match component {
            Component::Normal(part) => {
                let text = part
                    .to_str()
                    .ok_or_else(|| format!("Path is not valid UTF-8: {}", path.display()))?;
                if text.is_empty() {
                    return Err("Archive path contains an empty segment".to_string());
                }
                parts.push(text);
            }
            Component::CurDir => {}
            _ => {
                return Err(format!(
                    "Archive path contains unsupported segment: {}",
                    path.display()
                ));
            }
        }
    }

    if parts.is_empty() {
        return Err("Archive path cannot be empty".to_string());
    }

    Ok(parts.join("/"))
}

fn archive_key_to_relative_path(raw: &str) -> Result<PathBuf, String> {
    if raw.trim().is_empty() {
        return Err("Archive entry path is empty".to_string());
    }

    let mut path = PathBuf::new();
    for part in raw.split('/') {
        if part.is_empty() || part == "." || part == ".." {
            return Err(format!("Archive entry has invalid path segment: {}", raw));
        }
        path.push(part);
    }

    if path.as_os_str().is_empty() {
        return Err("Archive entry path resolved to empty".to_string());
    }

    Ok(path)
}

fn write_archive_path<W: Write>(writer: &mut W, relative_path: &Path) -> Result<(), String> {
    let key = relative_path_to_archive_key(relative_path)?;
    let key_bytes = key.as_bytes();
    if key_bytes.len() > u32::MAX as usize {
        return Err(format!(
            "Archive entry path is too long: {}",
            relative_path.display()
        ));
    }

    write_u32_le(writer, key_bytes.len() as u32)?;
    writer
        .write_all(key_bytes)
        .map_err(|e| format!("Failed to write archive path '{}': {}", key, e))
}

fn read_archive_path<R: Read>(reader: &mut R) -> Result<PathBuf, String> {
    let path_len = read_u32_le(reader)? as usize;
    if path_len == 0 {
        return Err("Archive entry path length is zero".to_string());
    }

    let mut path_bytes = vec![0_u8; path_len];
    reader
        .read_exact(&mut path_bytes)
        .map_err(|e| format!("Failed to read archive path: {}", e))?;

    let path_key = String::from_utf8(path_bytes)
        .map_err(|e| format!("Archive path is not valid UTF-8: {}", e))?;
    archive_key_to_relative_path(&path_key)
}

fn write_archive_header<W: Write>(writer: &mut W) -> Result<(), String> {
    writer
        .write_all(PROJECT_ARCHIVE_MAGIC)
        .map_err(|e| format!("Failed to write archive header: {}", e))?;
    write_u32_le(writer, PROJECT_ARCHIVE_VERSION)
}

fn read_archive_header<R: Read>(reader: &mut R) -> Result<u32, String> {
    let mut magic = [0_u8; PROJECT_ARCHIVE_MAGIC.len()];
    reader
        .read_exact(&mut magic)
        .map_err(|e| format!("Failed to read archive header: {}", e))?;
    if &magic != PROJECT_ARCHIVE_MAGIC {
        return Err("Invalid project archive signature".to_string());
    }

    let version = read_u32_le(reader)?;
    if version != 1 && version != 2 {
        return Err(format!(
            "Unsupported archive version: {} (supported: 1, 2)",
            version
        ));
    }
    Ok(version)
}

fn write_archive_dir_entry<W: Write>(writer: &mut W, relative_path: &Path) -> Result<(), String> {
    writer
        .write_all(&[ARCHIVE_ENTRY_DIR])
        .map_err(|e| format!("Failed to write archive directory marker: {}", e))?;
    write_archive_path(writer, relative_path)
}

struct CountingReader<R: Read> {
    inner: R,
    bytes_read: Arc<AtomicU64>,
}

impl<R: Read> CountingReader<R> {
    fn new(inner: R, bytes_read: Arc<AtomicU64>) -> Self {
        Self { inner, bytes_read }
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let count = self.inner.read(buf)?;
        if count > 0 {
            self.bytes_read.fetch_add(count as u64, Ordering::Relaxed);
        }
        Ok(count)
    }
}

fn copy_reader_to_writer_with_progress<R: Read, W: Write>(
    reader: &mut R,
    writer: &mut W,
    expected_size: u64,
    on_progress: &mut dyn FnMut(u64) -> Result<(), String>,
) -> Result<u64, String> {
    let mut copied_bytes = 0_u64;
    let mut buffer = [0_u8; EXPORT_COPY_BUFFER_SIZE];

    loop {
        let read_count = reader
            .read(&mut buffer)
            .map_err(|e| format!("Failed reading archive source data: {}", e))?;
        if read_count == 0 {
            break;
        }

        writer
            .write_all(&buffer[..read_count])
            .map_err(|e| format!("Failed writing archive data: {}", e))?;
        copied_bytes = copied_bytes.saturating_add(read_count as u64);
        on_progress(read_count as u64)?;
    }

    if copied_bytes != expected_size {
        return Err(format!(
            "Copy size mismatch while writing archive data (expected {}, copied {})",
            expected_size, copied_bytes
        ));
    }

    Ok(copied_bytes)
}

#[derive(Debug, Clone)]
struct PreparedArchiveFile {
    relative_path: PathBuf,
    source_path: PathBuf,
    original_size: u64,
    compressed_temp_path: PathBuf,
    compressed_size: u64,
}

fn collect_project_entries(
    source_root: &Path,
    current_dir: &Path,
    dirs: &mut Vec<PathBuf>,
    files: &mut Vec<(PathBuf, PathBuf, u64)>,
) -> Result<(), String> {
    let mut entries: Vec<_> = fs::read_dir(current_dir)
        .map_err(|e| format!("Failed to read directory {}: {}", current_dir.display(), e))?
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    entries.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    for entry in entries {
        let path = entry.path();
        let file_type = entry
            .file_type()
            .map_err(|e| format!("Failed to inspect {}: {}", path.display(), e))?;
        let relative_path = path
            .strip_prefix(source_root)
            .map_err(|e| format!("Failed to compute archive path {}: {}", path.display(), e))?
            .to_path_buf();

        if file_type.is_dir() {
            dirs.push(relative_path);
            collect_project_entries(source_root, &path, dirs, files)?;
        } else if file_type.is_file() {
            let original_size = entry
                .metadata()
                .map_err(|e| format!("Failed to read metadata for {}: {}", path.display(), e))?
                .len();
            files.push((relative_path, path, original_size));
        } else if file_type.is_symlink() {
            return Err(format!(
                "Project archive does not support symbolic links: {}",
                path.display()
            ));
        }
    }

    Ok(())
}

fn compress_file_to_gzip(
    source_path: &Path,
    destination_path: &Path,
    cancel_flag: &AtomicBool,
    progress_tx: &mpsc::Sender<u64>,
) -> Result<u64, String> {
    if cancel_flag.load(Ordering::Relaxed) {
        return Err(EXPORT_CANCELLED_SENTINEL.to_string());
    }

    let source_file = fs::File::open(source_path)
        .map_err(|e| format!("Failed to open {}: {}", source_path.display(), e))?;
    let destination_file = fs::File::create(destination_path)
        .map_err(|e| format!("Failed to create {}: {}", destination_path.display(), e))?;

    let mut reader = BufReader::new(source_file);
    let writer = BufWriter::new(destination_file);
    let mut encoder = GzEncoder::new(writer, Compression::fast());

    let mut buffer = [0_u8; EXPORT_COPY_BUFFER_SIZE];
    let mut pending_progress = 0_u64;
    loop {
        if cancel_flag.load(Ordering::Relaxed) {
            return Err(EXPORT_CANCELLED_SENTINEL.to_string());
        }

        let read_count = reader
            .read(&mut buffer)
            .map_err(|e| format!("Failed to read {}: {}", source_path.display(), e))?;
        if read_count == 0 {
            break;
        }

        encoder
            .write_all(&buffer[..read_count])
            .map_err(|e| format!("Failed to compress {}: {}", source_path.display(), e))?;
        pending_progress = pending_progress.saturating_add(read_count as u64);
        if pending_progress >= EXPORT_PROGRESS_EMIT_BYTES / 4 {
            let _ = progress_tx.send(pending_progress);
            pending_progress = 0;
        }
    }

    if pending_progress > 0 {
        let _ = progress_tx.send(pending_progress);
    }

    let mut writer = encoder.finish().map_err(|e| {
        format!(
            "Failed to finalize compression for {}: {}",
            source_path.display(),
            e
        )
    })?;
    writer.flush().map_err(|e| {
        format!(
            "Failed to flush compressed file {}: {}",
            destination_path.display(),
            e
        )
    })?;
    drop(writer);

    let compressed_size = fs::metadata(destination_path)
        .map_err(|e| {
            format!(
                "Failed to read metadata for {}: {}",
                destination_path.display(),
                e
            )
        })?
        .len();
    Ok(compressed_size)
}

fn read_archive_entry_type<R: Read>(reader: &mut R) -> Result<u8, String> {
    let mut kind = [0_u8; 1];
    reader.read_exact(&mut kind).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            "Archive ended unexpectedly before end marker".to_string()
        } else {
            format!("Failed to read archive entry marker: {}", e)
        }
    })?;
    Ok(kind[0])
}

fn extract_archive_file_to_path<R: Read>(
    reader: &mut R,
    destination_path: &Path,
    byte_count: u64,
    on_progress: &mut dyn FnMut() -> Result<(), String>,
) -> Result<(), String> {
    if let Some(parent) = destination_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
    }

    let output_file = fs::File::create(destination_path)
        .map_err(|e| format!("Failed to create {}: {}", destination_path.display(), e))?;
    let mut writer = BufWriter::new(output_file);

    let mut limited_reader = reader.take(byte_count);
    let mut emit_progress = |_delta: u64| on_progress();
    copy_reader_to_writer_with_progress(
        &mut limited_reader,
        &mut writer,
        byte_count,
        &mut emit_progress,
    )
    .map_err(|e| {
        format!(
            "Failed to extract {} from archive: {}",
            destination_path.display(),
            e
        )
    })?;
    writer
        .flush()
        .map_err(|e| format!("Failed to flush {}: {}", destination_path.display(), e))?;
    on_progress()?;

    Ok(())
}

fn extract_gzip_archive_file_to_path<R: Read>(
    reader: &mut R,
    destination_path: &Path,
    compressed_bytes: u64,
    expected_size: u64,
    on_progress: &mut dyn FnMut() -> Result<(), String>,
) -> Result<(), String> {
    if let Some(parent) = destination_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
    }

    let output_file = fs::File::create(destination_path)
        .map_err(|e| format!("Failed to create {}: {}", destination_path.display(), e))?;
    let mut writer = BufWriter::new(output_file);

    let mut limited_reader = reader.take(compressed_bytes);
    let mut decoder = GzDecoder::new(&mut limited_reader);
    let mut copied = 0_u64;
    let mut buffer = [0_u8; EXPORT_COPY_BUFFER_SIZE];
    loop {
        let read_count = decoder.read(&mut buffer).map_err(|e| {
            format!(
                "Failed to extract gzip payload for {}: {}",
                destination_path.display(),
                e
            )
        })?;
        if read_count == 0 {
            break;
        }
        writer
            .write_all(&buffer[..read_count])
            .map_err(|e| format!("Failed to write {}: {}", destination_path.display(), e))?;
        copied = copied.saturating_add(read_count as u64);
        on_progress()?;
    }
    drop(decoder);

    // Ensure the source stream advances exactly by the declared payload size.
    if limited_reader.limit() > 0 {
        let mut drain_buffer = [0_u8; EXPORT_COPY_BUFFER_SIZE];
        while limited_reader.limit() > 0 {
            let read_count = limited_reader
                .read(&mut drain_buffer)
                .map_err(|e| format!("Failed to finalize archive payload stream: {}", e))?;
            if read_count == 0 {
                break;
            }
            on_progress()?;
        }
    }
    if limited_reader.limit() != 0 {
        return Err(format!(
            "Archive is truncated while reading compressed payload for {}",
            destination_path.display()
        ));
    }
    writer
        .flush()
        .map_err(|e| format!("Failed to flush {}: {}", destination_path.display(), e))?;

    if copied != expected_size {
        return Err(format!(
            "Archive payload size mismatch for {} (expected {}, got {})",
            destination_path.display(),
            expected_size,
            copied
        ));
    }
    on_progress()?;

    Ok(())
}

fn extract_project_archive_to_directory(
    archive_path: &Path,
    destination_dir: &Path,
    on_progress: &mut dyn FnMut(u64, u64) -> Result<(), String>,
) -> Result<(), String> {
    let mut archive_file = fs::File::open(archive_path)
        .map_err(|e| format!("Failed to open archive {}: {}", archive_path.display(), e))?;
    let archive_total_bytes = archive_file
        .metadata()
        .map_err(|e| format!("Failed to read archive metadata {}: {}", archive_path.display(), e))?
        .len();
    let mut magic = [0_u8; 4];
    let read_count = archive_file
        .read(&mut magic)
        .map_err(|e| format!("Failed to inspect archive header: {}", e))?;
    archive_file
        .seek(SeekFrom::Start(0))
        .map_err(|e| format!("Failed to rewind archive stream: {}", e))?;

    let read_counter = Arc::new(AtomicU64::new(0));
    let counting_reader = CountingReader::new(archive_file, Arc::clone(&read_counter));
    let mut decoder: Box<dyn Read> =
        if read_count >= GZIP_MAGIC_PREFIX.len() && magic[..2] == GZIP_MAGIC_PREFIX {
            Box::new(GzDecoder::new(BufReader::new(counting_reader)))
        } else {
            Box::new(BufReader::new(counting_reader))
        };
    let mut emit_stream_progress = || {
        let processed = read_counter.load(Ordering::Relaxed).min(archive_total_bytes);
        on_progress(processed, archive_total_bytes)
    };
    emit_stream_progress()?;

    let archive_version = read_archive_header(&mut decoder)?;
    emit_stream_progress()?;

    loop {
        let entry_kind = read_archive_entry_type(&mut decoder)?;
        emit_stream_progress()?;
        if entry_kind == ARCHIVE_ENTRY_END {
            break;
        }

        let relative_path = read_archive_path(&mut decoder)?;
        emit_stream_progress()?;
        let output_path = destination_dir.join(&relative_path);
        if !output_path.starts_with(destination_dir) {
            return Err(format!(
                "Archive entry escaped destination directory: {}",
                relative_path.display()
            ));
        }

        match entry_kind {
            ARCHIVE_ENTRY_DIR => {
                fs::create_dir_all(&output_path).map_err(|e| {
                    format!(
                        "Failed to create directory {}: {}",
                        output_path.display(),
                        e
                    )
                })?;
                emit_stream_progress()?;
            }
            ARCHIVE_ENTRY_FILE => {
                if archive_version == 1 {
                    let file_size = read_u64_le(&mut decoder)?;
                    extract_archive_file_to_path(
                        &mut decoder,
                        &output_path,
                        file_size,
                        &mut emit_stream_progress,
                    )?;
                } else {
                    let mut compression_method = [0_u8; 1];
                    decoder
                        .read_exact(&mut compression_method)
                        .map_err(|e| format!("Failed to read archive compression method: {}", e))?;
                    let original_size = read_u64_le(&mut decoder)?;
                    let payload_size = read_u64_le(&mut decoder)?;
                    emit_stream_progress()?;

                    match compression_method[0] {
                        ARCHIVE_FILE_COMPRESSION_NONE => {
                            extract_archive_file_to_path(
                                &mut decoder,
                                &output_path,
                                payload_size,
                                &mut emit_stream_progress,
                            )?;
                            let actual_size = fs::metadata(&output_path)
                                .map_err(|e| {
                                    format!(
                                        "Failed to read metadata for {}: {}",
                                        output_path.display(),
                                        e
                                    )
                                })?
                                .len();
                            if actual_size != original_size {
                                return Err(format!(
                                    "Archive payload size mismatch for {} (expected {}, got {})",
                                    output_path.display(),
                                    original_size,
                                    actual_size
                                ));
                            }
                        }
                        ARCHIVE_FILE_COMPRESSION_GZIP => {
                            extract_gzip_archive_file_to_path(
                                &mut decoder,
                                &output_path,
                                payload_size,
                                original_size,
                                &mut emit_stream_progress,
                            )?;
                        }
                        other => {
                            return Err(format!(
                                "Unsupported file compression method in archive: {}",
                                other
                            ));
                        }
                    }
                }
            }
            _ => {
                return Err(format!("Unsupported archive entry type: {}", entry_kind));
            }
        }
    }
    emit_stream_progress()?;

    Ok(())
}

/// Recursively computes the total size of regular files inside a directory.
fn calculate_directory_size(path: &Path) -> Result<u64, String> {
    if !path.exists() {
        return Ok(0);
    }
    if !path.is_dir() {
        return Err(format!("Path is not a directory: {}", path.display()));
    }

    let mut total = 0_u64;
    let entries = fs::read_dir(path)
        .map_err(|e| format!("Failed to read directory {}: {}", path.display(), e))?;

    for entry in entries {
        let entry = entry.map_err(|e| e.to_string())?;
        let entry_path = entry.path();
        let file_type = entry
            .file_type()
            .map_err(|e| format!("Failed to inspect {}: {}", entry_path.display(), e))?;

        if file_type.is_dir() {
            total = total.saturating_add(calculate_directory_size(&entry_path)?);
            continue;
        }
        if file_type.is_file() {
            let metadata = entry.metadata().map_err(|e| {
                format!(
                    "Failed to read metadata for {}: {}",
                    entry_path.display(),
                    e
                )
            })?;
            total = total.saturating_add(metadata.len());
        }
    }

    Ok(total)
}

/// Fills missing `storage_bytes` values for legacy calculations.
fn hydrate_missing_calculation_sizes(
    project: &mut Project,
    project_dir: &Path,
) -> Result<bool, String> {
    let mut changed = false;

    for variant in &mut project.cif_variants {
        for calculation in &mut variant.calculations {
            if calculation.storage_bytes.is_some() {
                continue;
            }

            let calc_dir = project_dir.join("calculations").join(&calculation.id);
            let size = calculate_directory_size(&calc_dir)?;
            calculation.storage_bytes = Some(size);
            changed = true;
        }
    }

    Ok(changed)
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
    let duration = SystemTime::now().duration_since(UNIX_EPOCH).unwrap();
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
        let calculation_count: usize = project
            .cif_variants
            .iter()
            .map(|v| v.calculations.len())
            .sum();

        let formula = project.cif_variants.first().map(|v| v.formula.clone());

        let last_activity = project
            .cif_variants
            .iter()
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

    let mut project: Project = serde_json::from_str(&content)
        .map_err(|e| format!("Failed to parse project.json: {}", e))?;

    if hydrate_missing_calculation_sizes(&mut project, &project_dir)? {
        let updated_project = serde_json::to_string_pretty(&project)
            .map_err(|e| format!("Failed to serialize project: {}", e))?;
        fs::write(&project_json, updated_project)
            .map_err(|e| format!("Failed to write project.json: {}", e))?;
    }

    Ok(project)
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
    fs::write(
        structures_dir.join(format!("{}.cif", cif_id)),
        &cif_data.content,
    )
    .map_err(|e| format!("Failed to write CIF file: {}", e))?;

    // Save parsed crystal data JSON
    let crystal_json = serde_json::to_string_pretty(&cif_data.crystal_data)
        .map_err(|e| format!("Failed to serialize crystal data: {}", e))?;
    fs::write(
        structures_dir.join(format!("{}.json", cif_id)),
        crystal_json,
    )
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
    let variant = project
        .cif_variants
        .iter_mut()
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
    let preserve_full_phonon_artifacts = calc_data.calc_type == "phonon"
        && calc_data
            .parameters
            .get("epw_preparation")
            .and_then(|value| value.get("preserve_full_artifacts"))
            .and_then(|value| value.as_bool())
            .or_else(|| {
                calc_data
                    .parameters
                    .get("preserve_full_artifacts")
                    .and_then(|value| value.as_bool())
            })
            .unwrap_or(false);
    if let Some(work_dir) = working_dir {
        let work_path = PathBuf::from(&work_dir);
        if work_path.exists() {
            copied_work_path = Some(work_path.clone());
            let tmp_dir = calc_dir.join("tmp");
            fs::create_dir_all(&tmp_dir)
                .map_err(|e| format!("Failed to create tmp directory: {}", e))?;

            // For phonons, default to compact artifacts unless EPW-prep explicitly asks for full data.
            if calc_data.calc_type == "phonon" && !preserve_full_phonon_artifacts {
                copy_compact_phonon_artifacts(&work_path, &tmp_dir)?;
            } else {
                // Copy all files and directories for non-phonon, or EPW-preserving phonon calculations.
                copy_dir_recursive(&work_path, &tmp_dir)?;
            }
        }
    }

    // For phonons, keep an explicit raw pipeline log and top-level output artifacts
    // so parsing fixes can be applied without depending on intermediate state.
    if calc_data.calc_type == "phonon" {
        fs::write(
            calc_dir.join("phonon_pipeline.out"),
            &calc_data.output_content,
        )
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

    let storage_bytes = Some(calculate_directory_size(&calc_dir)?);

    // Create calculation run record
    let calc_run = CalculationRun {
        id: calc_id,
        calc_type: calc_data.calc_type,
        parameters: calc_data.parameters,
        result: Some(calc_data.result),
        started_at: calc_data.started_at,
        completed_at: Some(calc_data.completed_at),
        tags: calc_data.tags,
        storage_bytes,
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
    let mut no_progress = |_delta: u64| Ok(());
    copy_dir_recursive_with_progress(src.as_path(), dst.as_path(), &mut no_progress)
}

fn copy_dir_recursive_with_progress(
    src: &Path,
    dst: &Path,
    on_progress: &mut dyn FnMut(u64) -> Result<(), String>,
) -> Result<(), String> {
    if !dst.exists() {
        fs::create_dir_all(dst)
            .map_err(|e| format!("Failed to create directory {}: {}", dst.display(), e))?;
    }

    for entry in fs::read_dir(src).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        let dest_path = dst.join(entry.file_name());

        if path.is_dir() {
            copy_dir_recursive_with_progress(&path, &dest_path, on_progress)?;
        } else if path.is_file() {
            let copied = fs::copy(&path, &dest_path)
                .map_err(|e| format!("Failed to copy {}: {}", path.display(), e))?;
            on_progress(copied)?;
        } else if path.is_symlink() {
            return Err(format!(
                "Copy operation does not support symbolic links: {}",
                path.display()
            ));
        }
    }

    Ok(())
}

/// Copies a file relative to a source root into a destination root, preserving subdirectories.
fn copy_file_relative_if_exists(
    src_root: &Path,
    dst_root: &Path,
    rel_path: &Path,
) -> Result<(), String> {
    let source = src_root.join(rel_path);
    if !source.exists() || !source.is_file() {
        return Ok(());
    }

    let destination = dst_root.join(rel_path);
    if let Some(parent) = destination.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
    }
    fs::copy(&source, &destination)
        .map_err(|e| format!("Failed to copy {}: {}", source.display(), e))?;
    Ok(())
}

fn parse_dynmat_index(name: &str) -> Option<u32> {
    let suffix = name.strip_prefix("dynmat")?;
    if suffix.is_empty() {
        return None;
    }
    suffix.parse::<u32>().ok()
}

fn find_latest_dynmat_file(tmp_dir: &Path) -> Option<PathBuf> {
    let mut latest: Option<(u32, PathBuf)> = None;

    let entries = fs::read_dir(tmp_dir).ok()?;
    for entry in entries {
        let entry = entry.ok()?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let name = entry.file_name().to_string_lossy().to_string();
        let idx = match parse_dynmat_index(&name) {
            Some(i) => i,
            None => continue,
        };
        match latest.as_ref() {
            Some((current, _)) if *current >= idx => {}
            _ => latest = Some((idx, path)),
        }
    }

    latest.map(|(_, path)| path)
}

fn count_irreducible_qpoints(tmp_dir: &Path) -> Result<u32, String> {
    let mut count = 0u32;
    let entries = fs::read_dir(tmp_dir)
        .map_err(|e| format!("Failed to read tmp directory {}: {}", tmp_dir.display(), e))?;

    for entry in entries {
        let entry = entry.map_err(|e| e.to_string())?;
        if !entry.path().is_file() {
            continue;
        }
        let name = entry.file_name().to_string_lossy().to_string();
        if let Some(idx) = parse_dynmat_index(&name) {
            // dynmat0 is metadata; dynmat1..N are irreducible q-points.
            if idx > 0 {
                count += 1;
            }
        }
    }

    Ok(count)
}

fn looks_like_completed_phonon_run(tmp_dir: &Path) -> bool {
    // QE writes this status XML when running with ldisp.
    let status_xml = tmp_dir
        .join("_ph0")
        .join("qcortado_scf.phsave")
        .join("status_run.xml");
    if let Ok(content) = fs::read_to_string(&status_xml) {
        if content.contains("STOPPED_IN") && content.contains("dynmatrix") {
            return true;
        }
        if content.contains("<RECOVER_CODE>30</RECOVER_CODE>") {
            return true;
        }
    }

    if let Some(dynmat_path) = find_latest_dynmat_file(tmp_dir) {
        if let Ok(content) = fs::read_to_string(dynmat_path) {
            if content.contains("JOB DONE") {
                return true;
            }
        }
    }

    false
}

fn copy_compact_phonon_artifacts(src_tmp_dir: &Path, staging_tmp_dir: &Path) -> Result<(), String> {
    fs::create_dir_all(staging_tmp_dir)
        .map_err(|e| format!("Failed to create staging directory: {}", e))?;

    // Always copy known top-level files if present.
    let top_level = [
        "ph.in",
        "q2r.in",
        "matdyn_dos.in",
        "matdyn_bands.in",
        "phonon_freq",
        "phonon_freq.gp",
        "phonon_dos",
        "force_constants",
        "ph_recover.stdout",
        "ph_recover.stderr",
        "q2r_recheck.stdout",
        "q2r_recheck.stderr",
        "matdyn_recheck.stdout",
        "matdyn_recheck.stderr",
    ];

    for file in top_level {
        copy_file_relative_if_exists(src_tmp_dir, staging_tmp_dir, Path::new(file))?;
    }

    // Copy all dynmatN files.
    let entries = fs::read_dir(src_tmp_dir).map_err(|e| {
        format!(
            "Failed to read tmp directory {}: {}",
            src_tmp_dir.display(),
            e
        )
    })?;
    for entry in entries {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name().to_string_lossy().to_string();
        if parse_dynmat_index(&name).is_none() {
            continue;
        }
        copy_file_relative_if_exists(src_tmp_dir, staging_tmp_dir, Path::new(&name))?;
    }

    // Copy lightweight QE recovery metadata from _ph0 (without heavy wavefunction files).
    for rel in [
        "_ph0/qcortado_scf.phsave/status_run.xml",
        "_ph0/qcortado_scf.phsave/control_ph.xml",
        "_ph0/qcortado_scf.phsave/patterns.xml",
        "_ph0/qcortado_scf.save/data-file-schema.xml",
    ] {
        copy_file_relative_if_exists(src_tmp_dir, staging_tmp_dir, Path::new(rel))?;
    }

    Ok(())
}

/// Copies contents of a directory to another directory (public helper for bands calculation)
pub fn copy_dir_contents(src: &PathBuf, dst: &PathBuf) -> Result<(), String> {
    if !dst.exists() {
        fs::create_dir_all(dst).map_err(|e| format!("Failed to create directory: {}", e))?;
    }

    for entry in fs::read_dir(src).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        let dest_path = dst.join(entry.file_name());

        if path.is_dir() {
            copy_dir_recursive(&path, &dest_path)?;
        } else {
            fs::copy(&path, &dest_path).map_err(|e| format!("Failed to copy file: {}", e))?;
        }
    }

    Ok(())
}

/// Recovers a completed phonon run from an existing temporary working directory
/// and saves it as a project calculation entry.
#[tauri::command]
pub fn recover_phonon_calculation(
    app: AppHandle,
    project_id: String,
    cif_id: String,
    working_dir: String,
    source_scf_id: Option<String>,
    q_grid: Option<[u32; 3]>,
    tr2_ph: Option<f64>,
    q_path: Option<String>,
) -> Result<CalculationRun, String> {
    let tmp_dir = PathBuf::from(&working_dir);
    if !tmp_dir.exists() || !tmp_dir.is_dir() {
        return Err(format!(
            "Recovery directory not found: {}",
            tmp_dir.display()
        ));
    }

    if !looks_like_completed_phonon_run(&tmp_dir) {
        return Err(
            "Could not confirm a completed phonon run. Expected JOB DONE in dynmat output or status_run.xml at dynmatrix stage."
                .to_string(),
        );
    }

    let dispersion_source = if tmp_dir.join("phonon_freq.gp").exists() {
        tmp_dir.join("phonon_freq.gp")
    } else if tmp_dir.join("phonon_freq").exists() {
        tmp_dir.join("phonon_freq")
    } else {
        return Err("Missing phonon_freq.gp / phonon_freq in recovery directory".to_string());
    };

    let dispersion = read_phonon_dispersion_file(&dispersion_source).map_err(|e| {
        format!(
            "Failed to parse recovered phonon dispersion from {}: {}",
            dispersion_source.display(),
            e
        )
    })?;

    let dos_data = {
        let dos_file = tmp_dir.join("phonon_dos");
        if dos_file.exists() {
            match read_phonon_dos_file(&dos_file) {
                Ok(dos) => Some(dos),
                Err(_) => None,
            }
        } else {
            None
        }
    };

    let raw_output = if let Some(dynmat_path) = find_latest_dynmat_file(&tmp_dir) {
        fs::read_to_string(&dynmat_path).unwrap_or_default()
    } else {
        fs::read_to_string(tmp_dir.join("ph_recover.stdout")).unwrap_or_default()
    };

    let n_modes = dispersion.n_modes;
    let irreducible_qpoints = count_irreducible_qpoints(&tmp_dir)?;
    let n_qpoints = if irreducible_qpoints > 0 {
        irreducible_qpoints
    } else {
        dispersion.n_qpoints as u32
    };

    let mut parameters = serde_json::json!({
        "source_scf_id": source_scf_id,
        "q_grid": q_grid,
        "tr2_ph": tr2_ph,
        "calculate_dos": dos_data.is_some(),
        "calculate_dispersion": true,
        "n_qpoints": n_qpoints,
        "n_modes": n_modes,
        "recovered_from_tmp": true,
        "recovery_source_dir": working_dir,
    });
    if let Some(path) = q_path {
        parameters["q_path"] = serde_json::Value::String(path);
    }

    let result = QEResult {
        converged: true,
        total_energy: None,
        fermi_energy: None,
        total_magnetization: None,
        forces: None,
        stress: None,
        n_scf_steps: None,
        wall_time_seconds: None,
        eigenvalues: None,
        raw_output: raw_output.clone(),
        band_data: None,
        dos_data: None,
        phonon_data: Some(serde_json::json!({
            "dos_data": dos_data,
            "dispersion_data": dispersion,
        })),
    };

    let input_content = fs::read_to_string(tmp_dir.join("ph.in")).unwrap_or_default();
    let completed_at = now_iso();
    let started_at = completed_at.clone();

    let calc_data = SaveCalculationData {
        calc_type: "phonon".to_string(),
        parameters,
        result,
        started_at,
        completed_at,
        input_content,
        output_content: raw_output,
        tags: vec![],
    };

    // Stage a compact tmp snapshot to avoid copying large wavefunction archives.
    let staging_tmp_dir =
        std::env::temp_dir().join(format!("qcortado_phonon_recovery_{}", generate_id()));
    copy_compact_phonon_artifacts(&tmp_dir, &staging_tmp_dir)?;

    let save_result = save_calculation(
        app,
        project_id,
        cif_id,
        calc_data,
        Some(staging_tmp_dir.to_string_lossy().to_string()),
    );

    // Cleanup staging directory regardless of save outcome.
    let _ = fs::remove_dir_all(&staging_tmp_dir);

    save_result
}

/// Deletes a project
#[tauri::command]
pub fn delete_project(app: AppHandle, project_id: String) -> Result<(), String> {
    let projects_dir = ensure_projects_dir(&app)?;
    let project_dir = projects_dir.join(&project_id);

    if !project_dir.exists() {
        return Err(format!("Project not found: {}", project_id));
    }

    fs::remove_dir_all(&project_dir).map_err(|e| format!("Failed to delete project: {}", e))?;

    Ok(())
}

/// Exports a project directory into a compressed `.qcproj` archive.
#[tauri::command]
pub async fn export_project_archive(
    app: AppHandle,
    project_id: String,
    destination_path: String,
    export_id: String,
) -> Result<ProjectArchiveExportResult, String> {
    let cancel_flag = register_project_export_cancel_flag(&export_id);
    let app_handle = app.clone();
    let export_id_for_cleanup = export_id.clone();

    let join_result = tauri::async_runtime::spawn_blocking(move || {
        let projects_dir = ensure_projects_dir(&app_handle)?;
        let project_dir = projects_dir.join(&project_id);
        if !project_dir.exists() || !project_dir.is_dir() {
            return Err(format!("Project not found: {}", project_id));
        }

        let destination = PathBuf::from(&destination_path);
        let destination_parent = destination
            .parent()
            .ok_or_else(|| format!("Invalid destination path: {}", destination.display()))?;
        if !destination_parent.exists() {
            return Err(format!(
                "Destination directory does not exist: {}",
                destination_parent.display()
            ));
        }

        let project_json_path = project_dir.join("project.json");
        let project_json = fs::read_to_string(&project_json_path)
            .map_err(|e| format!("Failed to read project.json: {}", e))?;
        let project: Project = serde_json::from_str(&project_json)
            .map_err(|e| format!("Failed to parse project.json: {}", e))?;

        emit_project_export_progress(&app_handle, &export_id, &project_id, "scanning", 0, 0);
        let compression_temp_root =
            std::env::temp_dir().join(format!("qcortado_export_compress_{}", generate_id()));
        fs::create_dir_all(&compression_temp_root).map_err(|e| {
            format!(
                "Failed to create temporary compression directory {}: {}",
                compression_temp_root.display(),
                e
            )
        })?;

        let export_attempt = (|| -> Result<ProjectArchiveExportResult, String> {
            if cancel_flag.load(Ordering::Relaxed) {
                return Err(EXPORT_CANCELLED_SENTINEL.to_string());
            }

            let mut dir_entries: Vec<PathBuf> = Vec::new();
            let mut source_files: Vec<(PathBuf, PathBuf, u64)> = Vec::new();
            collect_project_entries(
                &project_dir,
                &project_dir,
                &mut dir_entries,
                &mut source_files,
            )?;

            let total_bytes = source_files
                .iter()
                .fold(0_u64, |acc, (_, _, size)| acc.saturating_add(*size));
            emit_project_export_progress(
                &app_handle,
                &export_id,
                &project_id,
                "compressing",
                0,
                total_bytes,
            );

            let mut prepared_files: Vec<PreparedArchiveFile> = source_files
                .into_iter()
                .enumerate()
                .map(
                    |(idx, (relative_path, source_path, original_size))| PreparedArchiveFile {
                        relative_path,
                        source_path,
                        original_size,
                        compressed_temp_path: compression_temp_root.join(format!("{:08}.gz", idx)),
                        compressed_size: 0,
                    },
                )
                .collect();

            if !prepared_files.is_empty() {
                let thread_count = (export_compression_threads() as usize)
                    .max(1)
                    .min(prepared_files.len());
                let jobs: VecDeque<(usize, PathBuf, PathBuf)> = prepared_files
                    .iter()
                    .enumerate()
                    .map(|(idx, file)| {
                        (
                            idx,
                            file.source_path.clone(),
                            file.compressed_temp_path.clone(),
                        )
                    })
                    .collect();
                let jobs = Arc::new(Mutex::new(jobs));

                let (result_tx, result_rx) = mpsc::channel::<(usize, Result<u64, String>)>();
                let (worker_done_tx, worker_done_rx) = mpsc::channel::<()>();
                let (progress_tx, progress_rx) = mpsc::channel::<u64>();
                let mut worker_handles = Vec::with_capacity(thread_count);

                for _ in 0..thread_count {
                    let jobs = Arc::clone(&jobs);
                    let result_tx = result_tx.clone();
                    let worker_done_tx = worker_done_tx.clone();
                    let progress_tx = progress_tx.clone();
                    let cancel_flag = Arc::clone(&cancel_flag);
                    worker_handles.push(std::thread::spawn(move || {
                        loop {
                            if cancel_flag.load(Ordering::Relaxed) {
                                break;
                            }

                            let next_job = match jobs.lock() {
                                Ok(mut guard) => guard.pop_front(),
                                Err(_) => {
                                    let _ = result_tx.send((
                                        usize::MAX,
                                        Err("Export worker queue failed due to a poisoned lock."
                                            .to_string()),
                                    ));
                                    cancel_flag.store(true, Ordering::Relaxed);
                                    break;
                                }
                            };

                            let Some((index, source_path, temp_path)) = next_job else {
                                break;
                            };

                            let result = compress_file_to_gzip(
                                &source_path,
                                &temp_path,
                                cancel_flag.as_ref(),
                                &progress_tx,
                            );
                            if result.is_err() {
                                cancel_flag.store(true, Ordering::Relaxed);
                            }
                            let _ = result_tx.send((index, result));
                        }

                        let _ = worker_done_tx.send(());
                    }));
                }
                drop(result_tx);
                drop(worker_done_tx);
                drop(progress_tx);

                let mut first_error: Option<String> = None;
                let mut processed_bytes = 0_u64;
                let mut last_emitted_bytes = 0_u64;
                let mut last_emit_time = Instant::now();
                let mut completed_workers = 0_usize;

                while completed_workers < thread_count {
                    while let Ok((index, result)) = result_rx.try_recv() {
                        match result {
                            Ok(compressed_size) => {
                                if let Some(entry) = prepared_files.get_mut(index) {
                                    entry.compressed_size = compressed_size;
                                }
                            }
                            Err(err) => {
                                if first_error.is_none() {
                                    first_error = Some(err);
                                }
                                cancel_flag.store(true, Ordering::Relaxed);
                            }
                        }
                    }

                    while let Ok(delta) = progress_rx.try_recv() {
                        processed_bytes = processed_bytes.saturating_add(delta);
                    }

                    let should_emit = processed_bytes == total_bytes
                        || processed_bytes.saturating_sub(last_emitted_bytes)
                            >= EXPORT_PROGRESS_EMIT_BYTES
                        || last_emit_time.elapsed()
                            >= Duration::from_millis(EXPORT_PROGRESS_EMIT_INTERVAL_MS);
                    if should_emit {
                        last_emitted_bytes = processed_bytes;
                        last_emit_time = Instant::now();
                        emit_project_export_progress(
                            &app_handle,
                            &export_id,
                            &project_id,
                            "compressing",
                            processed_bytes,
                            total_bytes,
                        );
                    }

                    match worker_done_rx.recv_timeout(Duration::from_millis(100)) {
                        Ok(_) => completed_workers += 1,
                        Err(mpsc::RecvTimeoutError::Timeout) => {}
                        Err(mpsc::RecvTimeoutError::Disconnected) => break,
                    }
                }

                while let Ok((index, result)) = result_rx.try_recv() {
                    match result {
                        Ok(compressed_size) => {
                            if let Some(entry) = prepared_files.get_mut(index) {
                                entry.compressed_size = compressed_size;
                            }
                        }
                        Err(err) => {
                            if first_error.is_none() {
                                first_error = Some(err);
                            }
                        }
                    }
                }

                while let Ok(delta) = progress_rx.try_recv() {
                    processed_bytes = processed_bytes.saturating_add(delta);
                }

                emit_project_export_progress(
                    &app_handle,
                    &export_id,
                    &project_id,
                    "compressing",
                    processed_bytes.min(total_bytes),
                    total_bytes,
                );

                for handle in worker_handles {
                    handle
                        .join()
                        .map_err(|_| "Export worker thread panicked".to_string())?;
                }

                if let Some(err) = first_error {
                    return Err(err);
                }

                if cancel_flag.load(Ordering::Relaxed) {
                    return Err(EXPORT_CANCELLED_SENTINEL.to_string());
                }
            }

            emit_project_export_progress(
                &app_handle,
                &export_id,
                &project_id,
                "finalizing",
                total_bytes,
                total_bytes,
            );

            let archive_file = fs::File::create(&destination).map_err(|e| {
                format!("Failed to create archive {}: {}", destination.display(), e)
            })?;
            let mut writer = BufWriter::new(archive_file);

            write_archive_header(&mut writer)?;
            for relative_path in &dir_entries {
                if cancel_flag.load(Ordering::Relaxed) {
                    return Err(EXPORT_CANCELLED_SENTINEL.to_string());
                }
                write_archive_dir_entry(&mut writer, relative_path)?;
            }

            for file in &prepared_files {
                if cancel_flag.load(Ordering::Relaxed) {
                    return Err(EXPORT_CANCELLED_SENTINEL.to_string());
                }

                writer
                    .write_all(&[ARCHIVE_ENTRY_FILE])
                    .map_err(|e| format!("Failed to write archive file marker: {}", e))?;
                write_archive_path(&mut writer, &file.relative_path)?;
                writer
                    .write_all(&[ARCHIVE_FILE_COMPRESSION_GZIP])
                    .map_err(|e| format!("Failed to write archive compression method: {}", e))?;
                write_u64_le(&mut writer, file.original_size)?;
                write_u64_le(&mut writer, file.compressed_size)?;

                let compressed_file = fs::File::open(&file.compressed_temp_path).map_err(|e| {
                    format!(
                        "Failed to open temporary compressed file {}: {}",
                        file.compressed_temp_path.display(),
                        e
                    )
                })?;
                let mut compressed_reader = BufReader::new(compressed_file);
                let mut cancellation_check = |_delta: u64| -> Result<(), String> {
                    if cancel_flag.load(Ordering::Relaxed) {
                        return Err(EXPORT_CANCELLED_SENTINEL.to_string());
                    }
                    Ok(())
                };
                copy_reader_to_writer_with_progress(
                    &mut compressed_reader,
                    &mut writer,
                    file.compressed_size,
                    &mut cancellation_check,
                )?;
            }

            writer
                .write_all(&[ARCHIVE_ENTRY_END])
                .map_err(|e| format!("Failed to finalize archive stream: {}", e))?;
            writer
                .flush()
                .map_err(|e| format!("Failed to flush archive file: {}", e))?;
            drop(writer);

            let archive_size_bytes = fs::metadata(&destination)
                .map_err(|e| format!("Failed to read archive metadata: {}", e))?
                .len();
            emit_project_export_progress(
                &app_handle,
                &export_id,
                &project_id,
                "done",
                total_bytes,
                total_bytes,
            );

            Ok(ProjectArchiveExportResult {
                archive_path: destination.to_string_lossy().to_string(),
                project_id: project_id.clone(),
                project_name: project.name,
                archive_size_bytes,
            })
        })();

        let _ = fs::remove_dir_all(&compression_temp_root);

        if let Err(err) = &export_attempt {
            let _ = fs::remove_file(&destination);
            if err == EXPORT_CANCELLED_SENTINEL {
                return Err("Export canceled by user.".to_string());
            }
        }

        export_attempt
    })
    .await
    .map_err(|e| format!("Export task failed to join: {}", e));

    unregister_project_export_cancel_flag(&export_id_for_cleanup);
    join_result?
}

/// Requests cancellation of an in-flight project export.
#[tauri::command]
pub fn cancel_project_export(export_id: String) -> bool {
    request_project_export_cancel(&export_id)
}

/// Imports a compressed `.qcproj` archive into the projects directory.
#[tauri::command]
pub async fn import_project_archive(
    app: AppHandle,
    archive_path: String,
    import_id: String,
) -> Result<ProjectArchiveImportResult, String> {
    let app_handle = app.clone();
    tauri::async_runtime::spawn_blocking(move || {
        let archive = PathBuf::from(&archive_path);
        if !archive.exists() || !archive.is_file() {
            return Err(format!("Archive not found: {}", archive.display()));
        }

        let archive_total_bytes = archive
            .metadata()
            .map_err(|e| format!("Failed to read archive metadata: {}", e))?
            .len();
        emit_project_import_progress(&app_handle, &import_id, "extracting", 0, archive_total_bytes);

        let projects_dir = ensure_projects_dir(&app_handle)?;
        let extraction_dir = std::env::temp_dir().join(format!("qcortado_import_{}", generate_id()));
        fs::create_dir_all(&extraction_dir).map_err(|e| {
            format!(
                "Failed to create temporary import directory {}: {}",
                extraction_dir.display(),
                e
            )
        })?;

        let result = (|| -> Result<ProjectArchiveImportResult, String> {
            let mut last_extract_emitted = 0_u64;
            let mut last_extract_emit_time = Instant::now();
            extract_project_archive_to_directory(
                &archive,
                &extraction_dir,
                &mut |processed_bytes, total_bytes| {
                    let processed = processed_bytes.min(total_bytes);
                    let should_emit = processed == total_bytes
                        || processed.saturating_sub(last_extract_emitted) >= IMPORT_PROGRESS_EMIT_BYTES
                        || last_extract_emit_time.elapsed()
                            >= Duration::from_millis(IMPORT_PROGRESS_EMIT_INTERVAL_MS);
                    if should_emit {
                        last_extract_emitted = processed;
                        last_extract_emit_time = Instant::now();
                        emit_project_import_progress(
                            &app_handle,
                            &import_id,
                            "extracting",
                            processed,
                            total_bytes,
                        );
                    }
                    Ok(())
                },
            )?;
            emit_project_import_progress(
                &app_handle,
                &import_id,
                "extracting",
                archive_total_bytes,
                archive_total_bytes,
            );

            let project_json_path = extraction_dir.join("project.json");
            if !project_json_path.exists() {
                return Err(
                    "Archive does not contain a valid QCortado project (missing project.json)"
                        .to_string(),
                );
            }

            let project_json = fs::read_to_string(&project_json_path)
                .map_err(|e| format!("Failed to read imported project.json: {}", e))?;
            let mut project: Project = serde_json::from_str(&project_json)
                .map_err(|e| format!("Failed to parse imported project.json: {}", e))?;

            let original_id = project.id.trim().to_string();
            let mut resolved_id = if original_id.is_empty() {
                generate_id()
            } else {
                original_id.clone()
            };
            let mut imported_with_new_id = resolved_id != original_id;
            while projects_dir.join(&resolved_id).exists() {
                resolved_id = generate_id();
                imported_with_new_id = true;
            }

            if project.id != resolved_id {
                project.id = resolved_id.clone();
                let updated = serde_json::to_string_pretty(&project)
                    .map_err(|e| format!("Failed to serialize imported project: {}", e))?;
                fs::write(&project_json_path, updated)
                    .map_err(|e| format!("Failed to update imported project.json: {}", e))?;
            }

            let destination_dir = projects_dir.join(&resolved_id);
            if destination_dir.exists() {
                return Err(format!(
                    "Target project directory already exists: {}",
                    destination_dir.display()
                ));
            }

            let install_total_bytes = calculate_directory_size(&extraction_dir)?;
            emit_project_import_progress(&app_handle, &import_id, "installing", 0, install_total_bytes);

            let mut installed_bytes = 0_u64;
            let mut last_install_emitted = 0_u64;
            let mut last_install_emit_time = Instant::now();
            let import_copy_result = copy_dir_recursive_with_progress(
                extraction_dir.as_path(),
                destination_dir.as_path(),
                &mut |copied| {
                    installed_bytes = installed_bytes.saturating_add(copied);
                    let processed = installed_bytes.min(install_total_bytes);
                    let should_emit = processed == install_total_bytes
                        || processed.saturating_sub(last_install_emitted) >= IMPORT_PROGRESS_EMIT_BYTES
                        || last_install_emit_time.elapsed()
                            >= Duration::from_millis(IMPORT_PROGRESS_EMIT_INTERVAL_MS);
                    if should_emit {
                        last_install_emitted = processed;
                        last_install_emit_time = Instant::now();
                        emit_project_import_progress(
                            &app_handle,
                            &import_id,
                            "installing",
                            processed,
                            install_total_bytes,
                        );
                    }
                    Ok(())
                },
            );
            if let Err(err) = import_copy_result {
                let _ = fs::remove_dir_all(&destination_dir);
                return Err(format!("Failed to copy imported project files: {}", err));
            }

            emit_project_import_progress(
                &app_handle,
                &import_id,
                "installing",
                install_total_bytes,
                install_total_bytes,
            );

            fs::create_dir_all(destination_dir.join("structures"))
                .map_err(|e| format!("Failed to ensure structures directory exists: {}", e))?;
            fs::create_dir_all(destination_dir.join("calculations"))
                .map_err(|e| format!("Failed to ensure calculations directory exists: {}", e))?;

            emit_project_import_progress(&app_handle, &import_id, "done", 1, 1);

            Ok(ProjectArchiveImportResult {
                project_id: resolved_id,
                project_name: project.name,
                imported_with_new_id,
            })
        })();

        let _ = fs::remove_dir_all(&extraction_dir);
        result
    })
    .await
    .map_err(|e| format!("Import task failed to join: {}", e))?
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
    let variant = project
        .cif_variants
        .iter_mut()
        .find(|v| v.id == cif_id)
        .ok_or_else(|| format!("CIF variant not found: {}", cif_id))?;

    let calc_index = variant
        .calculations
        .iter()
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

    let crystal_json_path = project_dir
        .join("structures")
        .join(format!("{}.json", cif_id));
    if !crystal_json_path.exists() {
        return Err(format!("Crystal data not found for CIF: {}", cif_id));
    }

    let content = fs::read_to_string(&crystal_json_path)
        .map_err(|e| format!("Failed to read crystal data: {}", e))?;

    serde_json::from_str(&content).map_err(|e| format!("Failed to parse crystal data: {}", e))
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

    let cif_path = project_dir
        .join("structures")
        .join(format!("{}.cif", cif_id));
    if !cif_path.exists() {
        return Err(format!("CIF file not found: {}", cif_id));
    }

    fs::read_to_string(&cif_path).map_err(|e| format!("Failed to read CIF file: {}", e))
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
                return Err(format!(
                    "Calculation {} is not a phonon calculation",
                    calc_id
                ));
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
                        dos_data = Some(serde_json::to_value(dos).map_err(|e| {
                            format!("Failed to serialize recovered DOS data: {}", e)
                        })?);
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
                        dispersion_data = Some(serde_json::to_value(dispersion).map_err(|e| {
                            format!("Failed to serialize recovered dispersion data: {}", e)
                        })?);
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
                .ok_or_else(|| {
                    format!("Calculation not found while writing calc.json: {}", calc_id)
                })?;
            let calc_json = serde_json::to_string_pretty(calc)
                .map_err(|e| format!("Failed to serialize calculation: {}", e))?;
            fs::write(&calc_json_path, calc_json)
                .map_err(|e| format!("Failed to write calc.json: {}", e))?;
        }
    }

    Ok(response.unwrap_or_else(|| {
        serde_json::json!({
            "dos_data": serde_json::Value::Null,
            "dispersion_data": serde_json::Value::Null,
        })
    }))
}
