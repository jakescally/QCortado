//! WIEN2k `case.struct` generation and refinement metadata.
//!
//! QCortado writes an initial fixed-format structure file from the project's
//! parsed CIF structure. WIEN2k's own `setrmt_lapw`, `sgroup`, and `symmetry`
//! programs subsequently validate and refine it remotely; `cif2struct` is
//! intentionally not part of this path.

use std::collections::BTreeMap;
use std::slice;

use serde::{Deserialize, Serialize};
use spglib_sys as ffi;

use crate::symmetry::{self, Matrix3, SymmetryAnalyzeInput, SymmetryAtomInput, Vec3};

const BOHR_PER_ANGSTROM: f64 = 1.889_726_125_457_828_1;
const DEFAULT_NPT: u32 = 781;
const DEFAULT_RMT: f64 = 2.0;
const DEFAULT_SYMPREC: f64 = 1e-5;

#[derive(Debug, Clone, Deserialize)]
struct CrystalParameter {
    value: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct CrystalAtomSite {
    type_symbol: String,
    fract_x: f64,
    fract_y: f64,
    fract_z: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct CrystalRecord {
    cell_length_a: CrystalParameter,
    cell_length_b: CrystalParameter,
    cell_length_c: CrystalParameter,
    cell_angle_alpha: CrystalParameter,
    cell_angle_beta: CrystalParameter,
    cell_angle_gamma: CrystalParameter,
    atom_sites: Vec<CrystalAtomSite>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructureSiteOverride {
    pub site_index: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub npt: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r0: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rmt: Option<f64>,
}

fn default_sgroup_tolerance() -> f64 {
    1e-5
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructureControls {
    #[serde(default)]
    pub rmt_reduction_percent: f64,
    #[serde(default = "default_sgroup_tolerance")]
    pub sgroup_tolerance: f64,
    #[serde(default)]
    pub site_overrides: Vec<Wien2kStructureSiteOverride>,
}

impl Default for Wien2kStructureControls {
    fn default() -> Self {
        Self {
            rmt_reduction_percent: 0.0,
            sgroup_tolerance: default_sgroup_tolerance(),
            site_overrides: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructureSite {
    pub site_index: usize,
    pub symbol: String,
    pub atomic_number: u16,
    pub positions: Vec<Vec3>,
    pub npt: u32,
    pub r0: f64,
    pub rmt: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructureDraft {
    pub project_id: String,
    pub cif_id: String,
    pub case_name: String,
    pub lattice_type: String,
    pub spacegroup_number: i32,
    pub international_symbol: String,
    /// Lengths are in Angstrom and angles are in degrees for UI display.
    pub cell_parameters: [f64; 6],
    pub standardized_lattice: Matrix3,
    pub sites: Vec<Wien2kStructureSite>,
    pub controls: Wien2kStructureControls,
    pub struct_content: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kStructureStage {
    Rmt,
    Sgroup,
    Symmetry,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Wien2kStructureSessionPhase {
    Staged,
    RmtReady,
    SgroupReady,
    SymmetryReady,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructureSession {
    pub session_id: String,
    pub draft: Wien2kStructureDraft,
    pub remote_case_dir: String,
    pub remote_install_root: String,
    pub hpc_profile_id: String,
    pub phase: Wien2kStructureSessionPhase,
    pub current_struct: String,
    #[serde(default)]
    pub artifacts: BTreeMap<String, String>,
    #[serde(default)]
    pub transcript: Vec<String>,
    pub started_at: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructureStageResult {
    pub session_id: String,
    pub stage: Wien2kStructureStage,
    pub phase: Wien2kStructureSessionPhase,
    pub candidate_struct: String,
    pub sites: Vec<Wien2kStructureSite>,
    pub native_output: String,
    pub save_allowed: bool,
    #[serde(default)]
    pub diagnostics: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Wien2kStructurePreviewSummary {
    pub lattice_type: String,
    pub spacegroup_number: Option<i32>,
    pub spacegroup_symbol: Option<String>,
    pub cell_parameters: [f64; 6],
}

pub fn default_r0(atomic_number: u16) -> f64 {
    match atomic_number {
        1..=18 => 1e-4,
        19..=36 => 5e-5,
        37..=71 => 1e-5,
        _ => 5e-6,
    }
}

pub fn validate_controls(controls: &Wien2kStructureControls) -> Result<(), String> {
    if !controls.rmt_reduction_percent.is_finite()
        || controls.rmt_reduction_percent < 0.0
        || controls.rmt_reduction_percent >= 100.0
    {
        return Err("RMT reduction must be between 0 and 100 percent.".to_string());
    }
    if !controls.sgroup_tolerance.is_finite()
        || controls.sgroup_tolerance < 1e-7
        || controls.sgroup_tolerance > 1e-3
    {
        return Err("SGROUP tolerance must be between 1e-7 and 1e-3.".to_string());
    }
    for site in &controls.site_overrides {
        if let Some(npt) = site.npt {
            if npt == 0 || npt % 2 == 0 {
                return Err(format!(
                    "NPT for site {} must be a positive odd integer.",
                    site.site_index
                ));
            }
        }
        if let Some(r0) = site.r0 {
            if !r0.is_finite() || r0 <= 0.0 {
                return Err(format!("R0 for site {} must be positive.", site.site_index));
            }
        }
        if let Some(rmt) = site.rmt {
            if !rmt.is_finite() || rmt <= 0.0 {
                return Err(format!(
                    "RMT for site {} must be positive.",
                    site.site_index
                ));
            }
        }
        if let (Some(r0), Some(rmt)) = (site.r0, site.rmt) {
            if r0 >= rmt {
                return Err(format!(
                    "R0 for site {} must be smaller than RMT.",
                    site.site_index
                ));
            }
        }
    }
    Ok(())
}

pub fn draft_from_crystal_json(
    project_id: String,
    cif_id: String,
    case_name: String,
    crystal_json: serde_json::Value,
    controls: Wien2kStructureControls,
) -> Result<Wien2kStructureDraft, String> {
    let case_name = super::case_state::normalize_case_name(&case_name).ok_or_else(|| {
        "Case name may contain only letters, digits, '.', '_' and '-'.".to_string()
    })?;
    validate_controls(&controls)?;
    let crystal: CrystalRecord = serde_json::from_value(crystal_json)
        .map_err(|err| format!("Failed to read CIF crystal data for WIEN2k: {err}"))?;
    if crystal.atom_sites.is_empty() {
        return Err("The selected CIF has no atomic sites.".to_string());
    }
    let lattice = lattice_from_parameters(&crystal)?;
    let analyzed = symmetry::analyze_structure(SymmetryAnalyzeInput {
        lattice,
        atoms: crystal
            .atom_sites
            .iter()
            .map(|atom| SymmetryAtomInput {
                symbol: atom.type_symbol.clone(),
                position: [atom.fract_x, atom.fract_y, atom.fract_z],
            })
            .collect(),
        symprec: DEFAULT_SYMPREC,
        angle_tolerance: -1.0,
    })?;
    draft_from_standardized_result(project_id, cif_id, case_name, analyzed, controls)
}

pub fn draft_from_standardized_result(
    project_id: String,
    cif_id: String,
    case_name: String,
    analyzed: symmetry::SymmetryAnalyzeResult,
    controls: Wien2kStructureControls,
) -> Result<Wien2kStructureDraft, String> {
    validate_controls(&controls)?;
    let mut sites = grouped_standardized_sites(&analyzed)?;
    apply_site_controls(&mut sites, &controls)?;
    let cell_parameters = cell_parameters(&analyzed.standardized_conventional_lattice)?;
    let struct_content = write_struct_file(&case_name, &cell_parameters, &sites);

    Ok(Wien2kStructureDraft {
        project_id,
        cif_id,
        case_name,
        // A complete primitive-labelled draft is deliberately passed through
        // native sgroup before becoming an accepted source.
        lattice_type: "P".to_string(),
        spacegroup_number: analyzed.spacegroup_number,
        international_symbol: analyzed.international_symbol,
        cell_parameters,
        standardized_lattice: analyzed.standardized_conventional_lattice,
        sites,
        controls,
        struct_content,
    })
}

fn lattice_from_parameters(crystal: &CrystalRecord) -> Result<Matrix3, String> {
    let a = crystal.cell_length_a.value;
    let b = crystal.cell_length_b.value;
    let c = crystal.cell_length_c.value;
    let alpha = crystal.cell_angle_alpha.value.to_radians();
    let beta = crystal.cell_angle_beta.value.to_radians();
    let gamma = crystal.cell_angle_gamma.value.to_radians();
    if ![a, b, c]
        .iter()
        .all(|value| value.is_finite() && *value > 0.0)
    {
        return Err("CIF cell lengths must be positive finite values.".to_string());
    }
    if gamma.sin().abs() < 1e-12 {
        return Err("CIF cell angle gamma is invalid for WIEN2k conversion.".to_string());
    }
    let cx = c * beta.cos();
    let cy = c * (alpha.cos() - beta.cos() * gamma.cos()) / gamma.sin();
    let cz_sq = c * c - cx * cx - cy * cy;
    if cz_sq < -1e-8 {
        return Err("CIF cell angles do not define a valid cell.".to_string());
    }
    Ok([
        [a, 0.0, 0.0],
        [b * gamma.cos(), b * gamma.sin(), 0.0],
        [cx, cy, cz_sq.max(0.0).sqrt()],
    ])
}

fn grouped_standardized_sites(
    analyzed: &symmetry::SymmetryAnalyzeResult,
) -> Result<Vec<Wien2kStructureSite>, String> {
    let atoms = &analyzed.standardized_conventional_atoms;
    let mut positions: Vec<Vec3> = atoms.iter().map(|atom| atom.position).collect();
    let types: Vec<i32> = atoms.iter().map(|atom| atom.type_index).collect();
    let mut lattice = transpose(&analyzed.standardized_conventional_lattice);
    let dataset_ptr = unsafe {
        ffi::spgat_get_dataset(
            lattice.as_mut_ptr(),
            positions.as_mut_ptr(),
            types.as_ptr(),
            positions.len() as i32,
            DEFAULT_SYMPREC,
            -1.0,
        )
    };
    if dataset_ptr.is_null() {
        return Err("Failed to group standardized sites for WIEN2k.".to_string());
    }
    let dataset = unsafe { &*dataset_ptr };
    let equivalent = unsafe { slice::from_raw_parts(dataset.equivalent_atoms, positions.len()) };
    let mut by_representative: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (index, representative) in equivalent.iter().enumerate() {
        by_representative
            .entry(*representative)
            .or_default()
            .push(index);
    }
    let mut sites = Vec::new();
    for indices in by_representative.values() {
        let first = &atoms[indices[0]];
        let atomic_number = atomic_number(&first.symbol).ok_or_else(|| {
            format!(
                "Unsupported element symbol for WIEN2k structure: {}",
                first.symbol
            )
        })?;
        sites.push(Wien2kStructureSite {
            site_index: sites.len() + 1,
            symbol: first.symbol.clone(),
            atomic_number,
            positions: indices.iter().map(|index| atoms[*index].position).collect(),
            npt: DEFAULT_NPT,
            r0: default_r0(atomic_number),
            rmt: DEFAULT_RMT,
        });
    }
    unsafe { ffi::spg_free_dataset(dataset_ptr) };
    Ok(sites)
}

fn apply_site_controls(
    sites: &mut [Wien2kStructureSite],
    controls: &Wien2kStructureControls,
) -> Result<(), String> {
    for update in &controls.site_overrides {
        let site = sites
            .iter_mut()
            .find(|site| site.site_index == update.site_index)
            .ok_or_else(|| format!("No WIEN2k structure site {} exists.", update.site_index))?;
        if let Some(npt) = update.npt {
            site.npt = npt;
        }
        if let Some(r0) = update.r0 {
            site.r0 = r0;
        }
        if let Some(rmt) = update.rmt {
            site.rmt = rmt;
        }
        if site.r0 >= site.rmt {
            return Err(format!(
                "R0 for {} site {} must be smaller than RMT.",
                site.symbol, site.site_index
            ));
        }
    }
    Ok(())
}

pub fn update_site_rmts_from_struct(
    sites: &[Wien2kStructureSite],
    struct_content: &str,
) -> Vec<Wien2kStructureSite> {
    let mut updated = sites.to_vec();
    let values: Vec<f64> = struct_content
        .lines()
        .filter_map(|line| {
            let after = line.split("RMT=").nth(1)?;
            after.split_whitespace().next()?.parse::<f64>().ok()
        })
        .collect();
    for (site, rmt) in updated.iter_mut().zip(values) {
        site.rmt = rmt;
    }
    updated
}

pub fn parse_struct_summary(struct_content: &str) -> Result<Wien2kStructurePreviewSummary, String> {
    let mut lines = struct_content.lines();
    let _title = lines
        .next()
        .ok_or_else(|| "WIEN2k structure has no title line.".to_string())?;
    let header = lines
        .next()
        .ok_or_else(|| "WIEN2k structure has no lattice header.".to_string())?;
    let _mode = lines
        .next()
        .ok_or_else(|| "WIEN2k structure has no calculation mode line.".to_string())?;
    let cell_line = lines
        .next()
        .ok_or_else(|| "WIEN2k structure has no cell line.".to_string())?;
    let lattice_type = header
        .split_whitespace()
        .next()
        .ok_or_else(|| "WIEN2k lattice type is missing.".to_string())?
        .to_string();
    let fields = header.split_whitespace().collect::<Vec<_>>();
    let (spacegroup_number, spacegroup_symbol) = parse_header_spacegroup(&fields);
    let values = (0..6)
        .filter_map(|index| {
            let start = index * 10;
            cell_line
                .get(start..start + 10)
                .and_then(|value| value.trim().parse::<f64>().ok())
        })
        .collect::<Vec<_>>();
    if values.len() != 6 {
        return Err("WIEN2k structure cell parameters could not be parsed.".to_string());
    }
    Ok(Wien2kStructurePreviewSummary {
        lattice_type,
        spacegroup_number,
        spacegroup_symbol,
        cell_parameters: [
            values[0] / BOHR_PER_ANGSTROM,
            values[1] / BOHR_PER_ANGSTROM,
            values[2] / BOHR_PER_ANGSTROM,
            values[3],
            values[4],
            values[5],
        ],
    })
}

fn parse_header_spacegroup(fields: &[&str]) -> (Option<i32>, Option<String>) {
    for field in fields.iter().rev() {
        let digit_count = field
            .chars()
            .take_while(|character| character.is_ascii_digit())
            .count();
        if digit_count > 0 && digit_count < field.len() {
            return (
                field[..digit_count].parse::<i32>().ok(),
                Some(field[digit_count..].to_string()),
            );
        }
    }
    if fields.len() >= 2 {
        let symbol = fields[fields.len() - 1];
        let number = fields[fields.len() - 2].parse::<i32>().ok();
        if number.is_some() {
            return (number, Some(symbol.to_string()));
        }
    }
    (None, None)
}

fn write_struct_file(case_name: &str, cell: &[f64; 6], sites: &[Wien2kStructureSite]) -> String {
    let mut lines = Vec::new();
    lines.push(format!(
        "{:<80}",
        format!("QCortado WIEN2k structure: {case_name}")
    ));
    // The expanded standardized draft is P1 until native sgroup determines
    // the accepted WIEN2k lattice/space-group setting.
    lines.push(format!(
        "P   LATTICE,NONEQUIV.ATOMS {:3}  {:3} {:<8}",
        sites.len(),
        1,
        "P1"
    ));
    lines.push("MODE OF CALC=RELA unit=bohr".to_string());
    lines.push(format!(
        "{:10.6}{:10.6}{:10.6}{:10.6}{:10.6}{:10.6}",
        cell[0] * BOHR_PER_ANGSTROM,
        cell[1] * BOHR_PER_ANGSTROM,
        cell[2] * BOHR_PER_ANGSTROM,
        cell[3],
        cell[4],
        cell[5],
    ));
    for site in sites {
        let atom_index = -(site.site_index as i32);
        let first = site.positions[0];
        lines.push(format!(
            "ATOM{:4}: X={:10.8} Y={:10.8} Z={:10.8}",
            atom_index, first[0], first[1], first[2]
        ));
        lines.push(format!(
            "          MULT={:2}          ISPLIT=15",
            site.positions.len()
        ));
        for position in site.positions.iter().skip(1) {
            lines.push(format!(
                "    {:4}: X={:10.8} Y={:10.8} Z={:10.8}",
                atom_index, position[0], position[1], position[2]
            ));
        }
        lines.push(format!(
            "{:<10} NPT={:5}  R0={:>10} RMT={:10.5}   Z:{:10.5}",
            site.symbol,
            site.npt,
            format!("{:.9}", site.r0).trim_start_matches('0'),
            site.rmt,
            site.atomic_number as f64
        ));
        lines.push("LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000".to_string());
        lines.push("                     0.0000000 1.0000000 0.0000000".to_string());
        lines.push("                     0.0000000 0.0000000 1.0000000".to_string());
    }
    lines.push("   0      NUMBER OF SYMMETRY OPERATIONS".to_string());
    format!("{}\n", lines.join("\n"))
}

fn cell_parameters(lattice: &Matrix3) -> Result<[f64; 6], String> {
    fn norm(v: Vec3) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }
    fn angle(a: Vec3, b: Vec3) -> f64 {
        let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        (dot / (norm(a) * norm(b)))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }
    let a = norm(lattice[0]);
    let b = norm(lattice[1]);
    let c = norm(lattice[2]);
    if ![a, b, c]
        .iter()
        .all(|value| value.is_finite() && *value > 0.0)
    {
        return Err("Standardized cell is invalid.".to_string());
    }
    Ok([
        a,
        b,
        c,
        angle(lattice[1], lattice[2]),
        angle(lattice[0], lattice[2]),
        angle(lattice[0], lattice[1]),
    ])
}

fn transpose(matrix: &Matrix3) -> Matrix3 {
    [
        [matrix[0][0], matrix[1][0], matrix[2][0]],
        [matrix[0][1], matrix[1][1], matrix[2][1]],
        [matrix[0][2], matrix[1][2], matrix[2][2]],
    ]
}

fn atomic_number(symbol: &str) -> Option<u16> {
    const ELEMENTS: [&str; 118] = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
        "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
        "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
        "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ];
    let normalized = symbol.trim().trim_end_matches(|character: char| {
        character.is_ascii_digit() || character == '+' || character == '-'
    });
    ELEMENTS
        .iter()
        .position(|element| *element == normalized)
        .map(|index| index as u16 + 1)
}

pub fn stage_command(
    session: &Wien2kStructureSession,
    stage: Wien2kStructureStage,
    controls: &Wien2kStructureControls,
) -> Result<String, String> {
    validate_controls(controls)?;
    let root = shell_quote(&session.remote_install_root);
    let dir = shell_quote(&session.remote_case_dir);
    let case_name = shell_quote(&session.draft.case_name);
    let prefix =
        format!("cd {dir} && export WIENROOT={root} && export PATH=\"$WIENROOT:$PATH\" && ");
    let command = match stage {
        Wien2kStructureStage::Rmt => {
            let mut explicit = BTreeMap::<String, f64>::new();
            for override_value in &controls.site_overrides {
                if let Some(rmt) = override_value.rmt {
                    if let Some(site) = session
                        .draft
                        .sites
                        .iter()
                        .find(|site| site.site_index == override_value.site_index)
                    {
                        explicit.insert(site.symbol.clone(), rmt);
                    }
                }
            }
            let selection = if explicit.is_empty() {
                format!("-r {}", controls.rmt_reduction_percent)
            } else {
                let assignments = explicit
                    .iter()
                    .map(|(symbol, rmt)| format!("{symbol}:{rmt}"))
                    .collect::<Vec<_>>()
                    .join(",");
                format!("-a {}", shell_quote(&assignments))
            };
            format!(
                "{prefix}\"$WIENROOT/setrmt_lapw\" {case_name} {selection} && cp {case_name}.struct_setrmt {case_name}.struct && \"$WIENROOT/x\" nn"
            )
        }
        Wien2kStructureStage::Sgroup => format!(
            "{prefix}\"$WIENROOT/x\" sgroup -settol {}",
            controls.sgroup_tolerance
        ),
        Wien2kStructureStage::Symmetry => format!(
            "{prefix}test -f {case_name}.struct_sgroup && cp {case_name}.struct_sgroup {case_name}.struct && \"$WIENROOT/x\" symmetry"
        ),
    };
    Ok(command)
}

pub fn expected_candidate_suffix(stage: Wien2kStructureStage) -> &'static str {
    match stage {
        Wien2kStructureStage::Rmt => "struct",
        Wien2kStructureStage::Sgroup => "struct_sgroup",
        Wien2kStructureStage::Symmetry => "struct_st",
    }
}

pub fn validate_stage_transition(
    phase: Wien2kStructureSessionPhase,
    stage: Wien2kStructureStage,
) -> Result<Wien2kStructureSessionPhase, String> {
    match (phase, stage) {
        (Wien2kStructureSessionPhase::Staged, Wien2kStructureStage::Rmt)
        | (Wien2kStructureSessionPhase::RmtReady, Wien2kStructureStage::Rmt) => {
            Ok(Wien2kStructureSessionPhase::RmtReady)
        }
        (Wien2kStructureSessionPhase::RmtReady, Wien2kStructureStage::Sgroup)
        | (Wien2kStructureSessionPhase::SgroupReady, Wien2kStructureStage::Sgroup) => {
            Ok(Wien2kStructureSessionPhase::SgroupReady)
        }
        (Wien2kStructureSessionPhase::SgroupReady, Wien2kStructureStage::Symmetry)
        | (Wien2kStructureSessionPhase::SymmetryReady, Wien2kStructureStage::Symmetry) => {
            Ok(Wien2kStructureSessionPhase::SymmetryReady)
        }
        _ => Err(
            "Complete and accept the preceding WIEN2k structure stage before continuing."
                .to_string(),
        ),
    }
}

pub fn stage_diagnostics(stage: Wien2kStructureStage, output: &str) -> Vec<String> {
    let upper = output.to_ascii_uppercase();
    let mut diagnostics = Vec::new();
    if upper.contains("ERROR") {
        diagnostics.push(
            "WIEN2k reported an error; review the native output before retrying.".to_string(),
        );
    }
    if upper.contains("SHIFTED") || upper.contains("MOVE THE ORIGIN") {
        diagnostics.push("WIEN2k requires an origin shift; the structure cannot be saved until this is resolved.".to_string());
    }
    if stage == Wien2kStructureStage::Sgroup && upper.contains("WARNING") {
        diagnostics.push(
            "SGROUP emitted a warning; compare the proposed structure before continuing."
                .to_string(),
        );
    }
    diagnostics
}

fn shell_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn silicon_site() -> Wien2kStructureSite {
        Wien2kStructureSite {
            site_index: 1,
            symbol: "Si".to_string(),
            atomic_number: 14,
            positions: vec![[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
            npt: DEFAULT_NPT,
            r0: default_r0(14),
            rmt: DEFAULT_RMT,
        }
    }

    fn staged_session() -> Wien2kStructureSession {
        Wien2kStructureSession {
            session_id: "session".to_string(),
            draft: Wien2kStructureDraft {
                project_id: "project".to_string(),
                cif_id: "cif".to_string(),
                case_name: "Si".to_string(),
                lattice_type: "P".to_string(),
                spacegroup_number: 227,
                international_symbol: "Fd-3m".to_string(),
                cell_parameters: [5.43, 5.43, 5.43, 90.0, 90.0, 90.0],
                standardized_lattice: [[5.43, 0.0, 0.0], [0.0, 5.43, 0.0], [0.0, 0.0, 5.43]],
                sites: vec![silicon_site()],
                controls: Wien2kStructureControls::default(),
                struct_content: String::new(),
            },
            remote_case_dir: "/scratch/qcortado/project/wien2k/session".to_string(),
            remote_install_root: "/opt/WIEN2k".to_string(),
            hpc_profile_id: "remote".to_string(),
            phase: Wien2kStructureSessionPhase::Staged,
            current_struct: String::new(),
            artifacts: BTreeMap::new(),
            transcript: vec![],
            started_at: "2026-05-24T00:00:00Z".to_string(),
        }
    }

    #[test]
    fn structgen_r0_thresholds_match_wien2k_sources() {
        assert_eq!(default_r0(14), 1e-4);
        assert_eq!(default_r0(26), 5e-5);
        assert_eq!(default_r0(58), 1e-5);
        assert_eq!(default_r0(92), 5e-6);
    }

    #[test]
    fn advanced_controls_require_valid_npt_and_r0() {
        let invalid_npt = Wien2kStructureControls {
            site_overrides: vec![Wien2kStructureSiteOverride {
                site_index: 1,
                npt: Some(780),
                r0: None,
                rmt: None,
            }],
            ..Wien2kStructureControls::default()
        };
        assert!(validate_controls(&invalid_npt).is_err());
        let invalid_r0 = Wien2kStructureControls {
            site_overrides: vec![Wien2kStructureSiteOverride {
                site_index: 1,
                npt: None,
                r0: Some(2.1),
                rmt: Some(2.0),
            }],
            ..Wien2kStructureControls::default()
        };
        assert!(validate_controls(&invalid_r0).is_err());
    }

    #[test]
    fn stage_order_requires_native_review_sequence() {
        assert!(validate_stage_transition(
            Wien2kStructureSessionPhase::Staged,
            Wien2kStructureStage::Rmt
        )
        .is_ok());
        assert!(validate_stage_transition(
            Wien2kStructureSessionPhase::Staged,
            Wien2kStructureStage::Symmetry
        )
        .is_err());
        assert!(validate_stage_transition(
            Wien2kStructureSessionPhase::RmtReady,
            Wien2kStructureStage::Sgroup
        )
        .is_ok());
    }

    #[test]
    fn struct_writer_creates_fixed_format_draft_for_native_refinement() {
        let output = write_struct_file(
            "Si",
            &[5.43, 5.43, 5.43, 90.0, 90.0, 90.0],
            &[silicon_site()],
        );

        assert!(output.contains("P   LATTICE,NONEQUIV.ATOMS   1    1 P1"));
        assert!(output.contains("ATOM  -1: X="));
        assert!(output.contains("MULT= 2"));
        assert!(output.contains("ISPLIT=15"));
        assert!(output.contains("NPT=  781"));
        assert!(output.contains("R0=.000100000"));
        assert!(output.contains("RMT=   2.00000"));
        assert!(output.contains("0      NUMBER OF SYMMETRY OPERATIONS"));
    }

    #[test]
    fn native_stage_commands_are_built_from_transient_case_session() {
        let session = staged_session();
        let rmt = stage_command(
            &session,
            Wien2kStructureStage::Rmt,
            &Wien2kStructureControls::default(),
        )
        .expect("RMT command");
        let sgroup = stage_command(
            &session,
            Wien2kStructureStage::Sgroup,
            &Wien2kStructureControls::default(),
        )
        .expect("SGROUP command");
        let symmetry = stage_command(
            &session,
            Wien2kStructureStage::Symmetry,
            &Wien2kStructureControls::default(),
        )
        .expect("SYMMETRY command");

        assert!(rmt.contains("setrmt_lapw"));
        assert!(rmt.contains("-r 0"));
        assert!(rmt.contains("\"$WIENROOT/x\" nn"));
        assert!(sgroup.contains("\"$WIENROOT/x\" sgroup -settol 0.00001"));
        assert!(symmetry.contains("\"$WIENROOT/x\" symmetry"));
    }

    #[test]
    fn refined_struct_summary_reads_native_spacegroup_and_cell() {
        let refined = concat!(
            "accepted                                                                        \n",
            "F   LATTICE,NONEQUIV. ATOMS: 2  225Fm-3m   \n",
            "MODE OF CALC=RELA unit=bohr\n",
            " 10.314314 10.314314 10.314314 90.000000 90.000000 90.000000\n",
        );
        let summary = parse_struct_summary(refined).expect("refined structure summary");

        assert_eq!(summary.lattice_type, "F");
        assert_eq!(summary.spacegroup_number, Some(225));
        assert_eq!(summary.spacegroup_symbol.as_deref(), Some("Fm-3m"));
        assert!((summary.cell_parameters[0] - 5.4579).abs() < 0.001);
    }
}
