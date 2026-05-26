//! Remote engine installation records and verification helpers.
//!
//! These records describe where an engine is installed for a remote HPC
//! profile. They are runtime setup metadata, not calculation input schemas.

use serde::{Deserialize, Serialize};

use crate::engines::plugin::{EnginePlugin, EnginePluginManifest};
use crate::engines::types::{CalculationKind, EngineId, EngineImplementationStatus};
use crate::hpc::profile::HpcProfile;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineInstallation {
    pub engine_id: EngineId,
    pub hpc_profile_id: String,
    pub remote_install_root: String,
    /// Snapshot of the selected HPC profile root verified during setup.
    /// Engine run directories should derive below this as:
    /// `{root}/qcortado/{project_id}/{engine_id}/{run_id}`.
    pub remote_workspace_root: String,
    /// Snapshot of the selected HPC profile root verified during setup.
    /// Saved remote artifacts should derive below this as:
    /// `{root}/qcortado/{project_id}/{engine_id}/{calculation_id}`.
    pub remote_project_root: String,
    #[serde(default)]
    pub verified_executables: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version_hint: Option<String>,
    pub verified_at: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineInstallationRequest {
    pub engine_id: EngineId,
    pub hpc_profile_id: String,
    pub remote_install_root: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EngineInstallationVerification {
    pub success: bool,
    pub message: String,
    pub checked_executables: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version_hint: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AddEngineInstallationResult {
    pub installation: EngineInstallation,
    pub verification: EngineInstallationVerification,
}

impl EngineInstallationRequest {
    pub fn normalized_for_profile(
        self,
        profile: &HpcProfile,
    ) -> Result<ResolvedEngineInstallationRequest, String> {
        let hpc_profile_id = normalize_required_text(&self.hpc_profile_id, "HPC profile")?;
        if hpc_profile_id != profile.id {
            return Err("Resolved HPC profile did not match the requested profile.".to_string());
        }

        Ok(ResolvedEngineInstallationRequest {
            engine_id: self.engine_id,
            hpc_profile_id,
            remote_install_root: normalize_remote_path(
                &self.remote_install_root,
                "Install directory",
            )?,
            remote_workspace_root: normalize_remote_path(
                &profile.remote_workspace_root,
                "Remote workspace root",
            )?,
            remote_project_root: normalize_remote_path(
                &profile.remote_project_root,
                "Remote project root",
            )?,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedEngineInstallationRequest {
    pub engine_id: EngineId,
    pub hpc_profile_id: String,
    pub remote_install_root: String,
    pub remote_workspace_root: String,
    pub remote_project_root: String,
}

pub fn required_engine_executables(engine_id: EngineId) -> &'static [&'static str] {
    match engine_id {
        EngineId::Qe => &["pw.x", "bands.x", "dos.x", "projwfc.x"],
        EngineId::Wien2k => &["x", "init_lapw", "run_lapw", "runsp_lapw"],
    }
}

pub fn selectable_engine_manifests(
    installations: &[EngineInstallation],
) -> Vec<EnginePluginManifest> {
    let mut manifests = crate::engines::plugin::implemented_engine_manifests();
    if installations
        .iter()
        .any(|installation| installation.engine_id == EngineId::Wien2k)
    {
        let mut manifest = crate::engines::wien2k::plugin::WIEN2K_RESERVED_ENGINE_PLUGIN.manifest();
        manifest.descriptor.status = EngineImplementationStatus::Configured;
        manifest.descriptor.calculation_kinds =
            vec![CalculationKind::EngineSetup, CalculationKind::Scf];
        manifest.workflows.retain(|workflow| {
            matches!(
                workflow.kind,
                CalculationKind::EngineSetup | CalculationKind::Scf
            )
        });
        manifests.push(manifest);
    }
    manifests
}

pub fn selectable_engine_ids(installations: &[EngineInstallation]) -> Vec<EngineId> {
    selectable_engine_manifests(installations)
        .into_iter()
        .map(|manifest| manifest.descriptor.id)
        .collect()
}

pub fn is_engine_selectable(engine_id: EngineId, installations: &[EngineInstallation]) -> bool {
    selectable_engine_ids(installations).contains(&engine_id)
}

pub fn build_remote_verification_command(request: &ResolvedEngineInstallationRequest) -> String {
    match request.engine_id {
        EngineId::Qe => build_qe_verification_command(request),
        EngineId::Wien2k => build_wien2k_verification_command(request),
    }
}

pub fn parse_verification_output(
    engine_id: EngineId,
    remote_output: &str,
) -> EngineInstallationVerification {
    let success = remote_output
        .lines()
        .any(|line| line.trim() == "QCORTADO_ENGINE_VERIFY_OK");
    let message = remote_output
        .lines()
        .find_map(|line| line.strip_prefix("message="))
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .unwrap_or_else(|| {
            if success {
                "Engine installation verified."
            } else {
                "Engine installation verification failed."
            }
        })
        .to_string();
    let version_hint = remote_output
        .lines()
        .find_map(|line| line.strip_prefix("version="))
        .map(str::trim)
        .filter(|line| !line.is_empty() && *line != "unknown")
        .map(str::to_string);

    EngineInstallationVerification {
        success,
        message,
        checked_executables: required_engine_executables(engine_id)
            .iter()
            .map(|value| value.to_string())
            .collect(),
        version_hint,
    }
}

fn build_qe_verification_command(request: &ResolvedEngineInstallationRequest) -> String {
    build_engine_verification_command(
        &request.remote_install_root,
        &request.remote_workspace_root,
        &request.remote_project_root,
        required_engine_executables(EngineId::Qe),
        "if [ -x \"$root/pw.x\" ]; then bin=\"$root\"; else bin=\"$root/bin\"; fi",
        "version=$(\"$bin/pw.x\" --version 2>/dev/null | head -n 1); [ -n \"$version\" ] || version=unknown",
        "Quantum ESPRESSO",
    )
}

fn build_wien2k_verification_command(request: &ResolvedEngineInstallationRequest) -> String {
    build_engine_verification_command(
        &request.remote_install_root,
        &request.remote_workspace_root,
        &request.remote_project_root,
        required_engine_executables(EngineId::Wien2k),
        "bin=\"$root\"",
        "if [ -r \"$bin/WIEN2k_VERSION\" ]; then version=$(head -n 1 \"$bin/WIEN2k_VERSION\"); else version=unknown; fi",
        "WIEN2k",
    )
}

fn build_engine_verification_command(
    install_root: &str,
    workspace_root: &str,
    project_root: &str,
    executables: &[&str],
    bin_resolver: &str,
    version_probe: &str,
    label: &str,
) -> String {
    let executable_checks = executables
        .iter()
        .map(|executable| {
            format!(
                "if [ ! -x \"$bin/{exe}\" ]; then missing=\"$missing $bin/{exe}\"; fi",
                exe = executable
            )
        })
        .collect::<Vec<_>>()
        .join("; ");

    format!(
        "root={root}; workspace={workspace}; project={project}; ok=1; missing=''; message=''; \
if [ ! -d \"$root\" ]; then ok=0; message=\"{label} install directory not found: $root\"; fi; \
if [ \"$ok\" = 1 ]; then {bin_resolver}; {executable_checks}; if [ -n \"$missing\" ]; then ok=0; message=\"Missing executables:$missing\"; fi; fi; \
if [ \"$ok\" = 1 ]; then mkdir -p \"$workspace\" \"$project\" 2>/dev/null || ok=0; if [ \"$ok\" = 0 ]; then message=\"Failed to create remote project directories\"; fi; fi; \
if [ \"$ok\" = 1 ]; then if [ ! -w \"$workspace\" ] || [ ! -w \"$project\" ]; then ok=0; message=\"Remote workspace/project roots are not writable\"; fi; fi; \
if [ \"$ok\" = 1 ]; then {version_probe}; echo QCORTADO_ENGINE_VERIFY_OK; echo \"message={label} installation verified at $bin\"; echo \"version=$version\"; else echo QCORTADO_ENGINE_VERIFY_FAIL; echo \"message=$message\"; fi",
        root = shell_single_quote(install_root),
        workspace = shell_single_quote(workspace_root),
        project = shell_single_quote(project_root),
        label = label,
        bin_resolver = bin_resolver,
        executable_checks = executable_checks,
        version_probe = version_probe,
    )
}

fn normalize_required_text(input: &str, field: &str) -> Result<String, String> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        Err(format!("{} is required", field))
    } else {
        Ok(trimmed.to_string())
    }
}

fn normalize_remote_path(input: &str, field: &str) -> Result<String, String> {
    let trimmed = normalize_required_text(input, field)?;
    let normalized = trimmed.trim_end_matches('/').to_string();
    if normalized == "/" || normalized == "." {
        return Err(format!("{} is not a safe remote path", field));
    }
    Ok(normalized)
}

fn shell_single_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    let escaped = value.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wien2k_verifier_uses_native_top_level_scripts() {
        let request = ResolvedEngineInstallationRequest {
            engine_id: EngineId::Wien2k,
            hpc_profile_id: "andromeda".to_string(),
            remote_install_root: "/opt/WIEN2k".to_string(),
            remote_workspace_root: "/scratch/qcortado".to_string(),
            remote_project_root: "/project/qcortado".to_string(),
        };
        let command = build_remote_verification_command(&request);

        assert!(command.contains("$bin/x"));
        assert!(command.contains("$bin/init_lapw"));
        assert!(command.contains("$bin/run_lapw"));
        assert!(command.contains("$bin/runsp_lapw"));
        assert!(command.contains("WIEN2k_VERSION"));
        assert!(!command.contains("pseudo"));
    }

    #[test]
    fn verification_output_parses_success_and_version() {
        let parsed = parse_verification_output(
            EngineId::Wien2k,
            "QCORTADO_ENGINE_VERIFY_OK\nmessage=WIEN2k installation verified at /opt/WIEN2k\nversion=24.1\n",
        );

        assert!(parsed.success);
        assert_eq!(parsed.version_hint.as_deref(), Some("24.1"));
        assert!(parsed
            .checked_executables
            .contains(&"init_lapw".to_string()));
    }

    #[test]
    fn installed_wien2k_exposes_structure_setup_and_scf() {
        let installation = EngineInstallation {
            engine_id: EngineId::Wien2k,
            hpc_profile_id: "andromeda".to_string(),
            remote_install_root: "/opt/WIEN2k".to_string(),
            remote_workspace_root: "/scratch/qcortado".to_string(),
            remote_project_root: "/project/qcortado".to_string(),
            verified_executables: vec![],
            version_hint: None,
            verified_at: "2026-05-24T00:00:00Z".to_string(),
        };
        let manifest = selectable_engine_manifests(&[installation])
            .into_iter()
            .find(|entry| entry.descriptor.id == EngineId::Wien2k)
            .expect("configured WIEN2k manifest");

        assert_eq!(
            manifest.descriptor.status,
            EngineImplementationStatus::Configured
        );
        assert_eq!(
            manifest.descriptor.calculation_kinds,
            vec![CalculationKind::EngineSetup, CalculationKind::Scf]
        );
        assert_eq!(manifest.workflows.len(), 2);
        assert_eq!(manifest.workflows[0].kind, CalculationKind::EngineSetup);
        assert_eq!(manifest.workflows[1].kind, CalculationKind::Scf);
    }
}
