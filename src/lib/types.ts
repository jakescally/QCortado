// Crystal structure types for CIF parsing and visualization

export interface LatticeParameter {
  value: number;
  uncertainty?: number;
}

export interface AtomSite {
  label: string;
  type_symbol: string;
  fract_x: number;
  fract_y: number;
  fract_z: number;
  wyckoff_symbol?: string;
  symmetry_multiplicity?: number;
  occupancy: number;
}

export interface AnisotropicParams {
  label: string;
  type_symbol: string;
  beta_11: number;
  beta_22: number;
  beta_33: number;
  beta_12: number;
  beta_13: number;
  beta_23: number;
}

export interface Citation {
  title?: string;
  journal?: string;
  year?: number;
  volume?: string;
  page_first?: string;
  page_last?: string;
  authors: string[];
}

export interface CrystalData {
  // Names
  chemical_name_common?: string;
  formula_structural?: string;
  formula_sum?: string;
  structure_type?: string;

  // Cell parameters
  cell_length_a: LatticeParameter;
  cell_length_b: LatticeParameter;
  cell_length_c: LatticeParameter;
  cell_angle_alpha: LatticeParameter;
  cell_angle_beta: LatticeParameter;
  cell_angle_gamma: LatticeParameter;
  cell_volume?: number;
  cell_formula_units_Z?: number;

  // Space group
  space_group_HM?: string;
  space_group_IT_number?: number;

  // Physical properties
  density?: number;
  measurement_temperature?: number;

  // Source/database info
  database_code?: string;
  audit_creation_date?: string;
  citation?: Citation;

  // Crystal structure
  atom_sites: AtomSite[];
  symmetry_operations: string[];
  anisotropic_params: AnisotropicParams[];
}

export type SCFPreset = "standard" | "phonon" | "relax";

export type ExecutionMode = "local" | "hpc";
export type HpcAuthMethod = "ssh_key" | "password";
export type HpcResourceType = "cpu" | "gpu";
export type HpcResourceMode = "cpu_only" | "gpu_only" | "both";
export type HpcLauncher = "srun" | "mpirun";

export interface SlurmResourceRequest {
  resource_type: HpcResourceType;
  partition?: string | null;
  walltime?: string | null;
  nodes?: number | null;
  ntasks?: number | null;
  cpus_per_task?: number | null;
  memory_gb?: number | null;
  gpus?: number | null;
  qos?: string | null;
  account?: string | null;
  constraint?: string | null;
  module_preamble?: string | null;
  additional_sbatch?: string[];
}

export interface HpcProfile {
  id: string;
  name: string;
  cluster: string;
  host: string;
  port: number;
  username: string;
  auth_method: HpcAuthMethod;
  ssh_key_path?: string | null;
  remote_qe_bin_dir: string;
  remote_pseudo_dir: string;
  remote_workspace_root: string;
  remote_project_root: string;
  resource_mode: HpcResourceMode;
  launcher: HpcLauncher;
  launcher_extra_args?: string | null;
  default_cpu_resources: SlurmResourceRequest;
  default_gpu_resources: SlurmResourceRequest;
  credential_persisted: boolean;
  warning_acknowledged: boolean;
  created_at: string;
  updated_at: string;
}

export interface HpcExecutionTarget {
  profile_id?: string | null;
  resources?: SlurmResourceRequest | null;
  interactive_debug?: boolean;
}

export interface ExecutionTarget {
  mode: ExecutionMode;
  hpc?: HpcExecutionTarget | null;
}

export interface HpcTaskMeta {
  backend?: string | null;
  remote_job_id?: string | null;
  scheduler_state?: string | null;
  remote_node?: string | null;
  remote_workdir?: string | null;
  remote_project_path?: string | null;
  remote_storage_bytes?: number | null;
}

export type QePositionUnit = "alat" | "bohr" | "angstrom" | "crystal";

export interface SavedStructureAtom {
  symbol: string;
  position: [number, number, number];
}

export interface SavedStructureData {
  position_units: QePositionUnit;
  atoms: SavedStructureAtom[];
  cell_parameters: [[number, number, number], [number, number, number], [number, number, number]] | null;
  cell_units: QePositionUnit | null;
}

export interface SavedCellSummary {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
  volume: number;
  units: string;
}

export interface OptimizedStructureOption {
  calcId: string;
  label: string;
  mode: "relax" | "vcrelax";
  completedAt: string | null;
  structure: SavedStructureData;
  cellSummary: SavedCellSummary | null;
}

// Element data for periodic table
export interface ElementInfo {
  symbol: string;
  name: string;
  atomicNumber: number;
  mass: number;
}

// Common element masses for QE calculations
export const ELEMENT_MASSES: Record<string, number> = {
  H: 1.008,
  He: 4.003,
  Li: 6.941,
  Be: 9.012,
  B: 10.81,
  C: 12.011,
  N: 14.007,
  O: 15.999,
  F: 18.998,
  Ne: 20.180,
  Na: 22.990,
  Mg: 24.305,
  Al: 26.982,
  Si: 28.086,
  P: 30.974,
  S: 32.065,
  Cl: 35.453,
  Ar: 39.948,
  K: 39.098,
  Ca: 40.078,
  Sc: 44.956,
  Ti: 47.867,
  V: 50.942,
  Cr: 51.996,
  Mn: 54.938,
  Fe: 55.845,
  Co: 58.933,
  Ni: 58.693,
  Cu: 63.546,
  Zn: 65.38,
  Ga: 69.723,
  Ge: 72.64,
  As: 74.922,
  Se: 78.96,
  Br: 79.904,
  Kr: 83.798,
  Rb: 85.468,
  Sr: 87.62,
  Y: 88.906,
  Zr: 91.224,
  Nb: 92.906,
  Mo: 95.96,
  Tc: 98.0,
  Ru: 101.07,
  Rh: 102.91,
  Pd: 106.42,
  Ag: 107.87,
  Cd: 112.41,
  In: 114.82,
  Sn: 118.71,
  Sb: 121.76,
  Te: 127.60,
  I: 126.90,
  Xe: 131.29,
  Cs: 132.91,
  Ba: 137.33,
  La: 138.91,
  Ce: 140.12,
  Pr: 140.91,
  Nd: 144.24,
  Pm: 145.0,
  Sm: 150.36,
  Eu: 151.96,
  Gd: 157.25,
  Tb: 158.93,
  Dy: 162.50,
  Ho: 164.93,
  Er: 167.26,
  Tm: 168.93,
  Yb: 173.05,
  Lu: 174.97,
  Hf: 178.49,
  Ta: 180.95,
  W: 183.84,
  Re: 186.21,
  Os: 190.23,
  Ir: 192.22,
  Pt: 195.08,
  Au: 196.97,
  Hg: 200.59,
  Tl: 204.38,
  Pb: 207.2,
  Bi: 208.98,
  Po: 209.0,
  At: 210.0,
  Rn: 222.0,
  Fr: 223.0,
  Ra: 226.0,
  Ac: 227.0,
  Th: 232.04,
  Pa: 231.04,
  U: 238.03,
};
