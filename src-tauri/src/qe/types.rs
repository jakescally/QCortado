//! Core types for representing Quantum ESPRESSO calculations.

use serde::{Deserialize, Serialize};

/// Bravais lattice types supported by QE (ibrav parameter)
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum BravaisLattice {
    /// Free lattice (ibrav=0) - requires CELL_PARAMETERS
    Free = 0,
    /// Cubic P (sc) (ibrav=1)
    CubicP = 1,
    /// Cubic F (fcc) (ibrav=2)
    CubicF = 2,
    /// Cubic I (bcc) (ibrav=3)
    CubicI = 3,
    /// Hexagonal (ibrav=4)
    Hexagonal = 4,
    /// Trigonal R, 3-fold axis c (ibrav=5)
    TrigonalR = 5,
    /// Tetragonal P (st) (ibrav=6)
    TetragonalP = 6,
    /// Tetragonal I (bct) (ibrav=7)
    TetragonalI = 7,
    /// Orthorhombic P (ibrav=8)
    OrthorhombicP = 8,
    /// Orthorhombic base-centered (ibrav=9)
    OrthorhombicBC = 9,
    /// Orthorhombic face-centered (ibrav=10)
    OrthorhombicFC = 10,
    /// Orthorhombic body-centered (ibrav=11)
    OrthorhombicI = 11,
    /// Monoclinic P (ibrav=12)
    MonoclinicP = 12,
    /// Monoclinic base-centered (ibrav=13)
    MonoclinicBC = 13,
    /// Triclinic (ibrav=14)
    Triclinic = 14,
}

/// Calculation types supported by pw.x
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Default)]
#[serde(rename_all = "lowercase")]
pub enum CalculationType {
    /// Self-consistent field calculation
    #[default]
    Scf,
    /// Non-self-consistent calculation (bands at arbitrary k-points)
    Nscf,
    /// Band structure calculation
    Bands,
    /// Structural relaxation (fixed cell)
    Relax,
    /// Molecular dynamics
    Md,
    /// Variable-cell relaxation
    VcRelax,
    /// Variable-cell molecular dynamics
    VcMd,
}

impl CalculationType {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Scf => "scf",
            Self::Nscf => "nscf",
            Self::Bands => "bands",
            Self::Relax => "relax",
            Self::Md => "md",
            Self::VcRelax => "vc-relax",
            Self::VcMd => "vc-md",
        }
    }
}

/// Units for atomic positions
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Default)]
#[serde(rename_all = "lowercase")]
pub enum PositionUnits {
    /// Units of alat (lattice parameter a)
    Alat,
    /// Bohr atomic units
    Bohr,
    /// Angstroms
    #[default]
    Angstrom,
    /// Crystal (fractional) coordinates
    Crystal,
}

impl PositionUnits {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Alat => "alat",
            Self::Bohr => "bohr",
            Self::Angstrom => "angstrom",
            Self::Crystal => "crystal",
        }
    }
}

/// K-point specification types
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum KPoints {
    /// Gamma point only
    Gamma,
    /// Automatic Monkhorst-Pack grid
    Automatic {
        grid: [u32; 3],
        /// Offset (0 or 1 for each direction)
        offset: [u32; 3],
    },
    /// Explicit k-point list (crystal coordinates)
    Crystal { points: Vec<KPoint> },
    /// Explicit k-point list (2pi/alat units)
    Tpiba { points: Vec<KPoint> },
    /// Band structure path (crystal_b)
    CrystalB { path: Vec<BandPathPoint> },
    /// Band structure path (tpiba_b)
    TpibaB { path: Vec<BandPathPoint> },
}

impl Default for KPoints {
    fn default() -> Self {
        Self::Gamma
    }
}

/// A single k-point with weight
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct KPoint {
    pub k: [f64; 3],
    pub weight: f64,
}

/// A point in a band structure path
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BandPathPoint {
    pub k: [f64; 3],
    /// Number of points to this k-point from previous (0 = endpoint)
    pub npoints: u32,
    /// Optional label (e.g., "G", "X", "L")
    pub label: Option<String>,
}

/// An atomic species (element type)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AtomicSpecies {
    /// Element symbol (e.g., "Si", "O", "Fe")
    pub symbol: String,
    /// Atomic mass in a.m.u.
    pub mass: f64,
    /// Pseudopotential filename
    pub pseudopotential: String,
}

/// An atom in the structure
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Atom {
    /// Element symbol (must match an AtomicSpecies)
    pub symbol: String,
    /// Position [x, y, z]
    pub position: [f64; 3],
    /// Whether each coordinate is free to move (for relaxation)
    /// [true, true, true] means fully free
    #[serde(default = "default_if_pos")]
    pub if_pos: [bool; 3],
}

fn default_if_pos() -> [bool; 3] {
    [true, true, true]
}

/// 3x3 matrix for cell parameters
pub type CellMatrix = [[f64; 3]; 3];

/// Occupations scheme
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(rename_all = "lowercase")]
pub enum Occupations {
    /// Fixed occupations
    #[default]
    Fixed,
    /// Smearing (for metals)
    Smearing,
    /// From input
    FromInput,
    /// Tetrahedra
    Tetrahedra,
}

/// Smearing types for metallic systems
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(rename_all = "kebab-case")]
pub enum SmearingType {
    /// Gaussian smearing
    #[default]
    Gaussian,
    /// Methfessel-Paxton
    MethfesselPaxton,
    /// Marzari-Vanderbilt (cold smearing)
    MarzariVanderbilt,
    /// Fermi-Dirac
    FermiDirac,
}

impl SmearingType {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Gaussian => "gaussian",
            Self::MethfesselPaxton => "methfessel-paxton",
            Self::MarzariVanderbilt => "marzari-vanderbilt",
            Self::FermiDirac => "fermi-dirac",
        }
    }
}

/// Complete system definition for a QE calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QESystem {
    /// Bravais lattice type
    pub ibrav: BravaisLattice,
    /// Lattice parameters celldm(1-6) in Bohr
    /// celldm(1) = a, celldm(2) = b/a, celldm(3) = c/a, etc.
    #[serde(default)]
    pub celldm: Option<[f64; 6]>,
    /// Cell parameters (if ibrav=0)
    #[serde(default)]
    pub cell_parameters: Option<CellMatrix>,
    /// Units for cell parameters
    #[serde(default)]
    pub cell_units: Option<PositionUnits>,
    /// Atomic species definitions
    pub species: Vec<AtomicSpecies>,
    /// Atoms in the unit cell
    pub atoms: Vec<Atom>,
    /// Units for atomic positions
    #[serde(default)]
    pub position_units: PositionUnits,
    /// Wavefunction cutoff in Ry
    pub ecutwfc: f64,
    /// Charge density cutoff in Ry (default: 4 * ecutwfc)
    #[serde(default)]
    pub ecutrho: Option<f64>,
    /// Spin polarization (1 = non-magnetic, 2 = collinear magnetic)
    #[serde(default = "default_nspin")]
    pub nspin: u32,
    /// Occupations scheme
    #[serde(default)]
    pub occupations: Occupations,
    /// Smearing type (if occupations = smearing)
    #[serde(default)]
    pub smearing: SmearingType,
    /// Smearing width in Ry
    #[serde(default)]
    pub degauss: Option<f64>,
}

fn default_nspin() -> u32 {
    1
}

/// A complete QE calculation configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QECalculation {
    /// Calculation type
    pub calculation: CalculationType,
    /// Prefix for output files
    pub prefix: String,
    /// Output directory
    pub outdir: String,
    /// Pseudopotential directory
    pub pseudo_dir: String,
    /// System definition
    pub system: QESystem,
    /// K-points specification
    pub kpoints: KPoints,
    /// SCF convergence threshold
    #[serde(default = "default_conv_thr")]
    pub conv_thr: f64,
    /// Mixing beta for SCF
    #[serde(default = "default_mixing_beta")]
    pub mixing_beta: f64,
    /// Whether to calculate forces
    #[serde(default)]
    pub tprnfor: bool,
    /// Whether to calculate stress
    #[serde(default)]
    pub tstress: bool,
    /// Force convergence threshold (for relax)
    #[serde(default)]
    pub forc_conv_thr: Option<f64>,
    /// Energy convergence threshold (for relax)
    #[serde(default)]
    pub etot_conv_thr: Option<f64>,
    /// Target pressure for vc-relax (kbar)
    #[serde(default)]
    pub press: Option<f64>,
    /// Verbosity level
    #[serde(default)]
    pub verbosity: Option<String>,
}

fn default_conv_thr() -> f64 {
    1.0e-6
}

fn default_mixing_beta() -> f64 {
    0.7
}

/// Result of a completed calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QEResult {
    /// Whether the calculation converged
    pub converged: bool,
    /// Total energy in Ry
    pub total_energy: Option<f64>,
    /// Fermi energy in eV (if applicable)
    pub fermi_energy: Option<f64>,
    /// Total magnetization (if spin-polarized)
    pub total_magnetization: Option<f64>,
    /// Forces on atoms [nat][3] in Ry/Bohr
    pub forces: Option<Vec<[f64; 3]>>,
    /// Stress tensor [3][3] in kbar
    pub stress: Option<CellMatrix>,
    /// Number of SCF iterations
    pub n_scf_steps: Option<u32>,
    /// Total wall time in seconds
    pub wall_time_seconds: Option<f64>,
    /// Band energies [nk][nbnd] in eV (for bands/nscf)
    pub eigenvalues: Option<Vec<Vec<f64>>>,
    /// Raw output text (for debugging)
    pub raw_output: String,
    /// Full band structure data (for bands calculations)
    #[serde(default)]
    pub band_data: Option<serde_json::Value>,
    /// Full phonon data (for phonon calculations)
    #[serde(default)]
    pub phonon_data: Option<serde_json::Value>,
}

impl Default for QEResult {
    fn default() -> Self {
        Self {
            converged: false,
            total_energy: None,
            fermi_energy: None,
            total_magnetization: None,
            forces: None,
            stress: None,
            n_scf_steps: None,
            wall_time_seconds: None,
            eigenvalues: None,
            raw_output: String::new(),
            band_data: None,
            phonon_data: None,
        }
    }
}

// ============================================================================
// Phonon Calculation Types
// ============================================================================

/// Configuration for ph.x phonon calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononCalculation {
    /// Prefix from the SCF calculation
    pub prefix: String,
    /// Output directory containing .save from SCF
    pub outdir: String,
    /// Base name for dynamical matrix files
    pub fildyn: String,
    /// Q-point grid dimensions [nq1, nq2, nq3]
    pub nq: [u32; 3],
    /// Convergence threshold for phonons (default 1e-12)
    pub tr2_ph: f64,
    /// Whether to compute phonon dispersion
    pub ldisp: bool,
    /// Whether to recover from an interrupted calculation
    pub recover: bool,
    /// Whether to transform the representation from crystal to Cartesian
    #[serde(default = "default_trans")]
    pub trans: bool,
    /// Compute dielectric tensor and Born effective charges (for insulators)
    #[serde(default)]
    pub epsil: bool,
    /// Compute polarizability derivatives
    #[serde(default)]
    pub fpol: bool,
    /// Compute Raman coefficients (requires non-metallic setup)
    #[serde(default)]
    pub lraman: bool,
    /// File root for dvscf perturbation potentials (useful for EPW workflows)
    #[serde(default)]
    pub fildvscf: Option<String>,
    /// Electron-phonon mode (e.g. "interpolated", "lambda_tetra")
    #[serde(default)]
    pub electron_phonon: Option<String>,
    /// Smearing width for electron-phonon integration
    #[serde(default)]
    pub el_ph_sigma: Option<f64>,
    /// Number of smearing values for electron-phonon integration
    #[serde(default)]
    pub el_ph_nsigma: Option<u32>,
    /// Number of old mixes kept in phonon self-consistency
    #[serde(default)]
    pub nmix_ph: Option<u32>,
    /// Maximum linear-response SCF iterations
    #[serde(default)]
    pub niter_ph: Option<u32>,
    /// Mixing factor for linear-response potential updates
    #[serde(default)]
    pub alpha_mix: Option<f64>,
    /// Start index of q-point to compute (1-based)
    #[serde(default)]
    pub start_q: Option<u32>,
    /// Last index of q-point to compute (1-based)
    #[serde(default)]
    pub last_q: Option<u32>,
    /// Start irreducible representation index (1-based)
    #[serde(default)]
    pub start_irr: Option<u32>,
    /// Last irreducible representation index (1-based)
    #[serde(default)]
    pub last_irr: Option<u32>,
    /// Verbosity level for ph.x output
    #[serde(default)]
    pub verbosity: Option<String>,
    /// Acoustic sum rule (none, simple, crystal, one-dim, zero-dim)
    #[serde(default = "default_asr")]
    pub asr: String,
}

fn default_asr() -> String {
    "crystal".to_string()
}

fn default_trans() -> bool {
    true
}

impl Default for PhononCalculation {
    fn default() -> Self {
        Self {
            prefix: "qcortado_scf".to_string(),
            outdir: "./tmp".to_string(),
            fildyn: "dynmat".to_string(),
            nq: [4, 4, 4],
            tr2_ph: 1.0e-12,
            ldisp: true,
            recover: false,
            asr: "crystal".to_string(),
            trans: true,
            epsil: false,
            fpol: false,
            lraman: false,
            fildvscf: None,
            electron_phonon: None,
            el_ph_sigma: None,
            el_ph_nsigma: None,
            nmix_ph: None,
            niter_ph: None,
            alpha_mix: None,
            start_q: None,
            last_q: None,
            start_irr: None,
            last_irr: None,
            verbosity: None,
        }
    }
}

/// Configuration for q2r.x (interatomic force constants)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Q2RCalculation {
    /// Base name for dynamical matrix files (from ph.x)
    pub fildyn: String,
    /// Output file for force constants
    pub flfrc: String,
    /// Acoustic sum rule: "simple", "crystal", "one-dim", "zero-dim"
    pub zasr: String,
}

impl Default for Q2RCalculation {
    fn default() -> Self {
        Self {
            fildyn: "dynmat".to_string(),
            flfrc: "force_constants".to_string(),
            zasr: "crystal".to_string(),
        }
    }
}

/// Configuration for matdyn.x (phonon DOS and dispersion)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatdynCalculation {
    /// Force constants file from q2r.x
    pub flfrc: String,
    /// Acoustic sum rule
    pub asr: String,
    /// Whether to calculate DOS
    pub dos: bool,
    /// Output file for phonon DOS
    pub fldos: Option<String>,
    /// K-point grid for DOS sampling [nk1, nk2, nk3]
    pub nk: Option<[u32; 3]>,
    /// Energy step for DOS (cm^-1)
    pub delta_e: Option<f64>,
    /// Q-path for dispersion calculation
    pub q_path: Option<Vec<QPathPoint>>,
    /// Output file for frequencies along path
    pub flfrq: Option<String>,
}

impl Default for MatdynCalculation {
    fn default() -> Self {
        Self {
            flfrc: "force_constants".to_string(),
            asr: "crystal".to_string(),
            dos: true,
            fldos: Some("phonon_dos".to_string()),
            nk: Some([20, 20, 20]),
            delta_e: Some(1.0),
            q_path: None,
            flfrq: None,
        }
    }
}

/// A point in the Q-path for phonon dispersion calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QPathPoint {
    /// Label for the high-symmetry point (e.g., "Γ", "X", "L")
    pub label: String,
    /// Coordinates in reciprocal lattice units
    pub coords: [f64; 3],
    /// Number of points in the segment to the next point
    pub npoints: u32,
}

/// Phonon density of states data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononDOS {
    /// Frequencies in cm^-1
    pub frequencies: Vec<f64>,
    /// DOS values
    pub dos: Vec<f64>,
    /// Maximum frequency
    pub omega_max: f64,
    /// Minimum frequency (can be negative for instabilities)
    pub omega_min: f64,
}

/// High-symmetry point marker for phonon dispersion
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononHighSymmetryMarker {
    /// Q-path distance
    pub q_distance: f64,
    /// Label (e.g., "Γ", "X", "M")
    pub label: String,
}

/// Phonon dispersion data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononDispersion {
    /// Q-point distances along the path
    pub q_points: Vec<f64>,
    /// Frequencies for each mode at each q-point: [mode_index][q_index]
    pub frequencies: Vec<Vec<f64>>,
    /// High-symmetry point markers
    pub high_symmetry_points: Vec<PhononHighSymmetryMarker>,
    /// Number of phonon modes (3 * nat)
    pub n_modes: usize,
    /// Number of q-points
    pub n_qpoints: usize,
    /// Frequency range [min, max] in cm^-1
    pub frequency_range: [f64; 2],
}

/// Combined phonon calculation result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononResult {
    /// Whether the calculation converged
    pub converged: bool,
    /// Number of q-points calculated
    pub n_qpoints: u32,
    /// Number of phonon modes
    pub n_modes: usize,
    /// Phonon DOS data (if calculated)
    pub dos_data: Option<PhononDOS>,
    /// Phonon dispersion data (if calculated)
    pub dispersion_data: Option<PhononDispersion>,
    /// Raw output from ph.x
    pub raw_output: String,
}

impl Default for PhononResult {
    fn default() -> Self {
        Self {
            converged: false,
            n_qpoints: 0,
            n_modes: 0,
            dos_data: None,
            dispersion_data: None,
            raw_output: String::new(),
        }
    }
}

/// Configuration for a complete phonon pipeline calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononPipelineConfig {
    /// ph.x configuration
    pub phonon: PhononCalculation,
    /// Whether to calculate DOS
    pub calculate_dos: bool,
    /// DOS grid dimensions (if calculating DOS)
    pub dos_grid: Option<[u32; 3]>,
    /// DOS energy/frequency step in cm^-1 (if calculating DOS)
    #[serde(default)]
    pub dos_delta_e: Option<f64>,
    /// Whether to calculate dispersion
    pub calculate_dispersion: bool,
    /// Q-path for dispersion (if calculating dispersion)
    pub q_path: Option<Vec<QPathPoint>>,
    /// Points per segment for dispersion
    pub points_per_segment: u32,
    /// Project ID containing the source SCF
    pub project_id: Option<String>,
    /// SCF calculation ID to get the .save directory from
    pub scf_calc_id: Option<String>,
}
