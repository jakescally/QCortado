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
        }
    }
}
