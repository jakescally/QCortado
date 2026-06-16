//! WIEN2k-native calculation types.
//!
//! These model WIEN2k's case-directory workflow instead of reusing QE input
//! concepts such as pseudopotentials or plane-wave cutoffs.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kSpinMode {
    NonSpinPolarized,
    SpinPolarized,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kParallelMode {
    #[default]
    Openmp,
    Kpoint,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kInitialSpinConfiguration {
    Up,
    Down,
    NonMagnetic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kDftUDoubleCounting {
    Amf,
    Sic,
    Hmf,
}

impl Default for Wien2kDftUDoubleCounting {
    fn default() -> Self {
        Self::Sic
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kFermiMethod {
    Tetra,
    Temp,
    Temps,
}

impl Default for Wien2kFermiMethod {
    fn default() -> Self {
        Self::Tetra
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Wien2kMixerMode {
    MSR1,
    MSEC3,
    MSEC4,
    MSR2,
    PRATT,
    PRAT0,
}

impl Default for Wien2kMixerMode {
    fn default() -> Self {
        Self::MSR1
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Wien2kMixerTrust {
    #[serde(rename = "default")]
    Default,
    STIFF,
    STIFFER,
    FAST,
}

impl Default for Wien2kMixerTrust {
    fn default() -> Self {
        Self::Default
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kDispersionCorrection {
    None,
    Dftd3,
    Dftd4,
}

impl Default for Wien2kDispersionCorrection {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStartingMagnetization {
    pub site_index: u32,
    pub element: String,
    pub configuration: Wien2kInitialSpinConfiguration,
    pub moment_bohr_magneton: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kHubbardTarget {
    pub site_index: u32,
    pub element: String,
    pub manifold: String,
    pub orbital_l: u8,
    pub u_ev: f64,
    pub j_ev: f64,
    #[serde(default)]
    pub recommended: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kDftUSettings {
    #[serde(default)]
    pub enabled: bool,
    #[serde(default)]
    pub double_counting: Wien2kDftUDoubleCounting,
    #[serde(default)]
    pub targets: Vec<Wien2kHubbardTarget>,
}

impl Default for Wien2kDftUSettings {
    fn default() -> Self {
        Self {
            enabled: false,
            double_counting: Wien2kDftUDoubleCounting::Sic,
            targets: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kMixerSettings {
    #[serde(default)]
    pub mode: Wien2kMixerMode,
    #[serde(default = "default_wien2k_mixer_greed")]
    pub greed: f64,
    #[serde(default = "default_wien2k_mixer_history")]
    pub history: u16,
    #[serde(default)]
    pub trust: Wien2kMixerTrust,
}

const fn default_wien2k_mixer_greed() -> f64 {
    0.2
}

const fn default_wien2k_mixer_history() -> u16 {
    8
}

impl Default for Wien2kMixerSettings {
    fn default() -> Self {
        Self {
            mode: Wien2kMixerMode::MSR1,
            greed: 0.2,
            history: 8,
            trust: Wien2kMixerTrust::Default,
        }
    }
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
        matches!(self, Self::ScfComplete | Self::BandsComplete | Self::Failed)
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
    #[serde(default)]
    pub fermi_method: Wien2kFermiMethod,
    #[serde(default)]
    pub fermi_smearing_ry: Option<f64>,
    #[serde(default)]
    pub starting_magnetization: Vec<Wien2kStartingMagnetization>,
}

impl Default for Wien2kInitializationSettings {
    fn default() -> Self {
        Self {
            rkmax: 7.0,
            gmax: 12.0,
            lmax: 10,
            k_mesh: [6, 6, 6],
            exchange_correlation: 13,
            lstart_energy_cutoff_ry: -6.0,
            spin_mode: Wien2kSpinMode::NonSpinPolarized,
            fermi_method: Wien2kFermiMethod::Tetra,
            fermi_smearing_ry: None,
            starting_magnetization: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kScfRunSettings {
    pub spin_mode: Wien2kSpinMode,
    #[serde(default)]
    pub parallel_mode: Wien2kParallelMode,
    pub max_iterations: u16,
    /// WIEN2k `run_lapw -ec`, in Ry.
    pub energy_convergence_ry: f64,
    /// WIEN2k `run_lapw -cc`, in electron charge units.
    pub charge_convergence: f64,
    #[serde(default)]
    pub force_convergence_mry_bohr: Option<f64>,
    #[serde(default)]
    pub dft_u: Wien2kDftUSettings,
    #[serde(default)]
    pub dispersion_correction: Wien2kDispersionCorrection,
    #[serde(default)]
    pub iterative_diagonalization: bool,
    #[serde(default)]
    pub force_minimization: bool,
    #[serde(default)]
    pub mixer: Wien2kMixerSettings,
}

impl Default for Wien2kScfRunSettings {
    fn default() -> Self {
        Self {
            spin_mode: Wien2kSpinMode::NonSpinPolarized,
            parallel_mode: Wien2kParallelMode::Openmp,
            max_iterations: 40,
            energy_convergence_ry: 0.0001,
            charge_convergence: 0.0001,
            force_convergence_mry_bohr: None,
            dft_u: Wien2kDftUSettings::default(),
            dispersion_correction: Wien2kDispersionCorrection::None,
            iterative_diagonalization: false,
            force_minimization: false,
            mixer: Wien2kMixerSettings::default(),
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
