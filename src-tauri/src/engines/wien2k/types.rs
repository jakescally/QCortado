//! Wien2k-native skeleton types.
//!
//! These types intentionally model WIEN2k's case-directory workflow instead of
//! trying to reuse QE input concepts. They are hidden and not wired to commands
//! yet; they document the future backend boundary for remote case management.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kSpinMode {
    NonSpinPolarized,
    SpinPolarized,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kCasePhase {
    Unstaged,
    StructStaged,
    InitializationRunning,
    Initialized,
    ScfRunning,
    ScfComplete,
    BandsRunning,
    BandsComplete,
    Failed,
}

impl Wien2kCasePhase {
    pub const fn is_terminal(self) -> bool {
        matches!(
            self,
            Self::ScfComplete | Self::BandsComplete | Self::Failed
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kCaseArtifactRole {
    Struct,
    InitializationInput,
    InitializationOutput,
    ScfOutput,
    Dayfile,
    Density,
    Vector,
    BandsOutput,
    Scratch,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kCommandProgram {
    InitLapw,
    RunLapw,
    RunspLapw,
    X,
}

impl Wien2kCommandProgram {
    pub const fn script_name(self) -> &'static str {
        match self {
            Self::InitLapw => "init_lapw",
            Self::RunLapw => "run_lapw",
            Self::RunspLapw => "runsp_lapw",
            Self::X => "x",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kCaseReference {
    pub case_name: String,
    pub remote_case_dir: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_scratch_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_archive_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_shadow_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub project_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cif_id: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kCaseArtifact {
    pub role: Wien2kCaseArtifactRole,
    pub basename: String,
    pub required_for_resume: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kInitializationSettings {
    /// WIEN2k `init_lapw -red` / setRMT reduction percentage.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rmt_reduction_percent: Option<f64>,
    /// WIEN2k `init_lapw -rkmax`.
    pub rkmax: f64,
    /// WIEN2k `init_lapw -gmax`.
    pub gmax: f64,
    /// WIEN2k `init_lapw -lmax`.
    pub lmax: u8,
    /// WIEN2k `init_lapw -numk 0 nx ny nz` style mesh.
    pub k_mesh: [u32; 3],
    /// WIEN2k `lstart` exchange-correlation selector passed by init_lapw.
    pub exchange_correlation: u16,
    /// WIEN2k `lstart` energy cutoff, in Ry.
    pub lstart_energy_cutoff_ry: f64,
    pub spin_mode: Wien2kSpinMode,
}

impl Default for Wien2kInitializationSettings {
    fn default() -> Self {
        Self {
            rmt_reduction_percent: None,
            rkmax: 7.0,
            gmax: 12.0,
            lmax: 10,
            k_mesh: [6, 6, 6],
            exchange_correlation: 13,
            lstart_energy_cutoff_ry: -6.0,
            spin_mode: Wien2kSpinMode::NonSpinPolarized,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kScfRunSettings {
    pub spin_mode: Wien2kSpinMode,
    pub max_iterations: u16,
    /// WIEN2k `run_lapw -ec`, in Ry.
    pub energy_convergence_ry: f64,
    /// WIEN2k `run_lapw -cc`, in electron charge units.
    pub charge_convergence: f64,
    #[serde(default)]
    pub force_convergence_mry_bohr: Option<f64>,
    #[serde(default)]
    pub parallel: bool,
}

impl Default for Wien2kScfRunSettings {
    fn default() -> Self {
        Self {
            spin_mode: Wien2kSpinMode::NonSpinPolarized,
            max_iterations: 40,
            energy_convergence_ry: 0.0001,
            charge_convergence: 0.0001,
            force_convergence_mry_bohr: None,
            parallel: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kBandsSettings {
    pub spin_mode: Wien2kSpinMode,
    #[serde(default)]
    pub spin_orbit: bool,
    #[serde(default)]
    pub parallel: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kRemoteRuntimeProfile {
    pub profile_id: String,
    pub wienroot: String,
    pub remote_workspace_root: String,
    pub remote_project_root: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_scratch_root: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kScfInput {
    pub case: Wien2kCaseReference,
    pub initialization: Wien2kInitializationSettings,
    pub run: Wien2kScfRunSettings,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kCommandPlan {
    pub program: Wien2kCommandProgram,
    pub argv: Vec<String>,
    pub working_directory: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub environment: Vec<(String, String)>,
    pub phase: Wien2kCasePhase,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub expected_artifacts: Vec<Wien2kCaseArtifact>,
}
