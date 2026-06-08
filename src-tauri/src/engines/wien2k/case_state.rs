//! Hidden WIEN2k case-state helpers.
//!
//! WIEN2k derives most file names from a case prefix and expects commands to
//! run inside that case directory. The platform should store only enough local
//! metadata to make the remote state resumable and comparable inside a QCortado
//! project.

use serde::{Deserialize, Serialize};

use super::types::{
    Wien2kCaseArtifact, Wien2kCaseArtifactRole, Wien2kCasePhase, Wien2kCaseReference,
};
use crate::engines::types::DiagnosticHint;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kCaseState {
    pub case: Wien2kCaseReference,
    pub phase: Wien2kCasePhase,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub artifacts: Vec<Wien2kCaseArtifact>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_job_id: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub diagnostics: Vec<DiagnosticHint>,
}

/// Returns a safe WIEN2k case prefix. The case name becomes the filename stem
/// for `case.struct`, `case.scf`, `case.dayfile`, and many generated files, so
/// path separators and whitespace are rejected instead of normalized silently.
pub fn normalize_case_name(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }

    let safe = trimmed
        .chars()
        .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.'));
    if safe {
        Some(trimmed.to_string())
    } else {
        None
    }
}

pub fn wien2k_case_file(case_name: &str, suffix: &str) -> String {
    format!("{case_name}.{suffix}")
}

pub fn core_case_artifacts(case_name: &str) -> Vec<Wien2kCaseArtifact> {
    vec![
        artifact(Wien2kCaseArtifactRole::Struct, case_name, "struct", true),
        artifact(
            Wien2kCaseArtifactRole::InitializationOutput,
            case_name,
            "outputnn",
            false,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationOutput,
            case_name,
            "outputst",
            false,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationOutput,
            case_name,
            "outputs",
            false,
        ),
        artifact(Wien2kCaseArtifactRole::Density, case_name, "clmsum", true),
        artifact(Wien2kCaseArtifactRole::ScfOutput, case_name, "scf", true),
        artifact(Wien2kCaseArtifactRole::Dayfile, case_name, "dayfile", false),
    ]
}

pub fn initialized_case_artifacts(case_name: &str) -> Vec<Wien2kCaseArtifact> {
    let mut artifacts = vec![
        artifact(Wien2kCaseArtifactRole::Struct, case_name, "struct", true),
        artifact(
            Wien2kCaseArtifactRole::InitializationOutput,
            case_name,
            "outputnn",
            false,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationOutput,
            case_name,
            "outputs",
            false,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationOutput,
            case_name,
            "outputst",
            false,
        ),
        artifact(Wien2kCaseArtifactRole::Density, case_name, "clmsum", true),
        artifact(Wien2kCaseArtifactRole::Dayfile, case_name, "dayfile", false),
    ];
    artifacts.extend([
        artifact(
            Wien2kCaseArtifactRole::InitializationInput,
            case_name,
            "in0",
            true,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationInput,
            case_name,
            "in1",
            true,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationInput,
            case_name,
            "in2",
            true,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationInput,
            case_name,
            "inc",
            true,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationInput,
            case_name,
            "inm",
            true,
        ),
        artifact(
            Wien2kCaseArtifactRole::InitializationInput,
            case_name,
            "klist",
            true,
        ),
    ]);
    artifacts
}

fn artifact(
    role: Wien2kCaseArtifactRole,
    case_name: &str,
    suffix: &str,
    required_for_resume: bool,
) -> Wien2kCaseArtifact {
    Wien2kCaseArtifact {
        role,
        basename: wien2k_case_file(case_name, suffix),
        required_for_resume,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn case_names_reject_paths_and_whitespace() {
        assert_eq!(
            normalize_case_name("Si_case-1.0"),
            Some("Si_case-1.0".to_string())
        );
        assert_eq!(normalize_case_name(""), None);
        assert_eq!(normalize_case_name("../case"), None);
        assert_eq!(normalize_case_name("case name"), None);
        assert_eq!(normalize_case_name("case/name"), None);
    }

    #[test]
    fn artifacts_use_wien2k_case_prefix() {
        let artifacts = initialized_case_artifacts("GaAs");
        assert!(artifacts
            .iter()
            .any(|artifact| artifact.basename == "GaAs.struct"));
        assert!(artifacts
            .iter()
            .any(|artifact| artifact.basename == "GaAs.clmsum"));
        assert!(artifacts
            .iter()
            .any(|artifact| artifact.basename == "GaAs.klist"));
    }
}
