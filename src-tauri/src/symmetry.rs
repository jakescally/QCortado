//! Symmetry-driven structure normalization backed by spglib.
//!
//! This module builds standardized conventional + primitive cells and
//! explicit reciprocal-space basis transforms so frontend viewers and QE
//! k-point generation can share a single source of truth.

use serde::{Deserialize, Serialize};
use spglib_sys as ffi;
use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::{c_char, c_int};

pub type Vec3 = [f64; 3];
pub type Matrix3 = [Vec3; 3];

const DEFAULT_SYMPREC: f64 = 1e-5;
const DEFAULT_ANGLE_TOLERANCE: f64 = -1.0;

fn default_symprec() -> f64 {
    DEFAULT_SYMPREC
}

fn default_angle_tolerance() -> f64 {
    DEFAULT_ANGLE_TOLERANCE
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SymmetryAtomInput {
    pub symbol: String,
    pub position: Vec3,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SymmetryAnalyzeInput {
    /// Real-space lattice vectors in Angstrom, row-wise [a, b, c].
    pub lattice: Matrix3,
    /// Atomic positions in fractional coordinates of `lattice`.
    pub atoms: Vec<SymmetryAtomInput>,
    /// Distance tolerance for symmetry finding (Angstrom).
    #[serde(default = "default_symprec")]
    pub symprec: f64,
    /// Angle tolerance for symmetry finding (degrees). -1 lets spglib choose.
    #[serde(default = "default_angle_tolerance")]
    pub angle_tolerance: f64,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SymmetryAtomOutput {
    pub symbol: String,
    pub position: Vec3,
    pub type_index: i32,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SymmetryAnalyzeResult {
    pub spacegroup_number: i32,
    pub hall_number: i32,
    pub international_symbol: String,
    pub choice: String,
    pub input_lattice: Matrix3,
    pub standardized_conventional_lattice: Matrix3,
    pub standardized_primitive_lattice: Matrix3,
    pub standardized_conventional_atoms: Vec<SymmetryAtomOutput>,
    pub standardized_primitive_atoms: Vec<SymmetryAtomOutput>,
    /// Matrix M such that k_input = M * k_primitive (fractional reciprocal coords).
    pub primitive_to_input_reciprocal: Matrix3,
    pub input_to_primitive_reciprocal: Matrix3,
    /// Matrix M such that k_std_conv = M * k_primitive.
    pub primitive_to_standardized_conventional_reciprocal: Matrix3,
    pub standardized_conventional_to_primitive_reciprocal: Matrix3,
    /// Spglib dataset transformation matrix (input -> standardized setting basis change).
    pub transformation_matrix: Matrix3,
    pub origin_shift: Vec3,
}

#[derive(Debug)]
struct DatasetMeta {
    spacegroup_number: i32,
    hall_number: i32,
    international_symbol: String,
    choice: String,
    transformation_matrix: Matrix3,
    origin_shift: Vec3,
}

fn normalize_symbol(symbol: &str) -> String {
    symbol
        .trim()
        .trim_end_matches(|ch: char| ch.is_ascii_digit() || ch == '+' || ch == '-')
        .to_string()
}

fn wrap01(mut x: f64) -> f64 {
    x = x - x.floor();
    if x < 0.0 {
        x += 1.0;
    }
    if (x - 1.0).abs() < 1e-12 {
        0.0
    } else {
        x
    }
}

fn wrap_fractional(v: Vec3) -> Vec3 {
    [wrap01(v[0]), wrap01(v[1]), wrap01(v[2])]
}

fn finite_vec3(v: &Vec3) -> bool {
    v[0].is_finite() && v[1].is_finite() && v[2].is_finite()
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn transpose(m: &Matrix3) -> Matrix3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

fn mat_mul(a: &Matrix3, b: &Matrix3) -> Matrix3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    out
}

fn inverse3(m: &Matrix3) -> Result<Matrix3, String> {
    let a = m[0][0];
    let b = m[0][1];
    let c = m[0][2];
    let d = m[1][0];
    let e = m[1][1];
    let f = m[1][2];
    let g = m[2][0];
    let h = m[2][1];
    let i = m[2][2];

    let co00 = e * i - f * h;
    let co01 = -(d * i - f * g);
    let co02 = d * h - e * g;
    let co10 = -(b * i - c * h);
    let co11 = a * i - c * g;
    let co12 = -(a * h - b * g);
    let co20 = b * f - c * e;
    let co21 = -(a * f - c * d);
    let co22 = a * e - b * d;

    let det = a * co00 + b * co01 + c * co02;
    if det.abs() < 1e-12 {
        return Err("Singular matrix encountered while computing basis transform".to_string());
    }
    let inv_det = 1.0 / det;

    Ok([
        [co00 * inv_det, co10 * inv_det, co20 * inv_det],
        [co01 * inv_det, co11 * inv_det, co21 * inv_det],
        [co02 * inv_det, co12 * inv_det, co22 * inv_det],
    ])
}

fn reciprocal_lattice(real: &Matrix3) -> Result<Matrix3, String> {
    let a1 = real[0];
    let a2 = real[1];
    let a3 = real[2];
    let volume = dot(a1, cross(a2, a3));
    if volume.abs() < 1e-12 {
        return Err(
            "Lattice vectors are nearly singular; cannot build reciprocal basis".to_string(),
        );
    }
    let scale = (2.0 * std::f64::consts::PI) / volume;
    Ok([
        {
            let v = cross(a2, a3);
            [v[0] * scale, v[1] * scale, v[2] * scale]
        },
        {
            let v = cross(a3, a1);
            [v[0] * scale, v[1] * scale, v[2] * scale]
        },
        {
            let v = cross(a1, a2);
            [v[0] * scale, v[1] * scale, v[2] * scale]
        },
    ])
}

fn reciprocal_fractional_transform(
    source_recip: &Matrix3,
    target_recip: &Matrix3,
) -> Result<Matrix3, String> {
    // f_target = (B_target^T)^-1 * B_source^T * f_source
    let source_t = transpose(source_recip);
    let target_t = transpose(target_recip);
    let target_t_inv = inverse3(&target_t)?;
    Ok(mat_mul(&target_t_inv, &source_t))
}

fn c_char_array_to_string(arr: &[c_char]) -> String {
    let len = arr
        .iter()
        .position(|&value| value == 0)
        .unwrap_or(arr.len());
    let bytes: Vec<u8> = arr[..len].iter().map(|&value| value as u8).collect();
    String::from_utf8_lossy(&bytes).trim().to_string()
}

fn spglib_error_message() -> String {
    let err = unsafe { ffi::spg_get_error_code() };
    let message_ptr = unsafe { ffi::spg_get_error_message(err) };
    if message_ptr.is_null() {
        return format!("spglib error code {}", err);
    }
    unsafe { CStr::from_ptr(message_ptr) }
        .to_string_lossy()
        .trim()
        .to_string()
}

fn get_dataset_meta(
    lattice: &Matrix3,
    positions: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<DatasetMeta, String> {
    let mut lattice_copy = *lattice;
    let mut positions_copy = positions.to_vec();
    let types_copy = types.to_vec();

    let dataset_ptr = unsafe {
        ffi::spgat_get_dataset(
            lattice_copy.as_mut_ptr(),
            positions_copy.as_mut_ptr(),
            types_copy.as_ptr(),
            positions_copy.len() as c_int,
            symprec,
            angle_tolerance,
        )
    };

    if dataset_ptr.is_null() {
        return Err(format!(
            "spglib dataset search failed: {}",
            spglib_error_message()
        ));
    }

    let dataset = unsafe { &*dataset_ptr };
    let meta = DatasetMeta {
        spacegroup_number: dataset.spacegroup_number as i32,
        hall_number: dataset.hall_number as i32,
        international_symbol: c_char_array_to_string(&dataset.international_symbol),
        choice: c_char_array_to_string(&dataset.choice),
        transformation_matrix: dataset.transformation_matrix,
        origin_shift: dataset.origin_shift,
    };
    unsafe { ffi::spg_free_dataset(dataset_ptr) };
    Ok(meta)
}

fn standardize_cell(
    lattice: &Matrix3,
    positions: &[Vec3],
    types: &[i32],
    to_primitive: bool,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<(Matrix3, Vec<Vec3>, Vec<i32>), String> {
    let nat = positions.len();
    let capacity = if to_primitive {
        nat.max(1)
    } else {
        nat.saturating_mul(4).max(4)
    };

    let mut lattice_buffer = *lattice;
    let mut position_buffer = vec![[0.0; 3]; capacity];
    let mut type_buffer = vec![0_i32; capacity];
    for idx in 0..nat {
        position_buffer[idx] = positions[idx];
        type_buffer[idx] = types[idx];
    }

    let standardized_count = unsafe {
        ffi::spgat_standardize_cell(
            lattice_buffer.as_mut_ptr(),
            position_buffer.as_mut_ptr(),
            type_buffer.as_mut_ptr(),
            nat as c_int,
            if to_primitive { 1 } else { 0 },
            0,
            symprec,
            angle_tolerance,
        )
    };

    if standardized_count <= 0 {
        return Err(format!(
            "spglib standardization failed: {}",
            spglib_error_message()
        ));
    }

    let count = standardized_count as usize;
    position_buffer.truncate(count);
    type_buffer.truncate(count);

    for pos in &mut position_buffer {
        *pos = wrap_fractional(*pos);
    }

    Ok((lattice_buffer, position_buffer, type_buffer))
}

fn atoms_with_symbols(
    positions: &[Vec3],
    types: &[i32],
    type_to_symbol: &HashMap<i32, String>,
) -> Vec<SymmetryAtomOutput> {
    positions
        .iter()
        .zip(types.iter())
        .map(|(position, type_index)| SymmetryAtomOutput {
            symbol: type_to_symbol
                .get(type_index)
                .cloned()
                .unwrap_or_else(|| format!("T{}", type_index)),
            position: *position,
            type_index: *type_index,
        })
        .collect()
}

pub fn analyze_structure(input: SymmetryAnalyzeInput) -> Result<SymmetryAnalyzeResult, String> {
    if input.atoms.is_empty() {
        return Err("Structure has no atoms".to_string());
    }
    if input.symprec <= 0.0 || !input.symprec.is_finite() {
        return Err("symprec must be a positive finite number".to_string());
    }
    if !input.angle_tolerance.is_finite() {
        return Err("angleTolerance must be finite".to_string());
    }
    for vector in &input.lattice {
        if !finite_vec3(vector) {
            return Err("Lattice contains non-finite values".to_string());
        }
    }

    let mut symbol_to_type: HashMap<String, i32> = HashMap::new();
    let mut type_to_symbol: HashMap<i32, String> = HashMap::new();
    let mut next_type = 1_i32;
    let mut positions: Vec<Vec3> = Vec::with_capacity(input.atoms.len());
    let mut types: Vec<i32> = Vec::with_capacity(input.atoms.len());

    for atom in &input.atoms {
        let symbol = normalize_symbol(&atom.symbol);
        if symbol.is_empty() {
            return Err("Encountered atom with an empty symbol".to_string());
        }
        if !finite_vec3(&atom.position) {
            return Err(format!(
                "Atom position for '{}' contains non-finite values",
                symbol
            ));
        }
        let type_index = if let Some(existing) = symbol_to_type.get(&symbol) {
            *existing
        } else {
            let created = next_type;
            next_type += 1;
            symbol_to_type.insert(symbol.clone(), created);
            type_to_symbol.insert(created, symbol.clone());
            created
        };

        positions.push(wrap_fractional(atom.position));
        types.push(type_index);
    }

    let input_lattice = input.lattice;
    let dataset_meta = get_dataset_meta(
        &input_lattice,
        &positions,
        &types,
        input.symprec,
        input.angle_tolerance,
    )?;

    let (std_conv_lattice, std_conv_positions, std_conv_types) = standardize_cell(
        &input_lattice,
        &positions,
        &types,
        false,
        input.symprec,
        input.angle_tolerance,
    )?;
    let (std_prim_lattice, std_prim_positions, std_prim_types) = standardize_cell(
        &input_lattice,
        &positions,
        &types,
        true,
        input.symprec,
        input.angle_tolerance,
    )?;

    let input_recip = reciprocal_lattice(&input_lattice)?;
    let std_conv_recip = reciprocal_lattice(&std_conv_lattice)?;
    let std_prim_recip = reciprocal_lattice(&std_prim_lattice)?;

    let primitive_to_input = reciprocal_fractional_transform(&std_prim_recip, &input_recip)?;
    let input_to_primitive = reciprocal_fractional_transform(&input_recip, &std_prim_recip)?;
    let primitive_to_std_conv = reciprocal_fractional_transform(&std_prim_recip, &std_conv_recip)?;
    let std_conv_to_primitive = reciprocal_fractional_transform(&std_conv_recip, &std_prim_recip)?;

    Ok(SymmetryAnalyzeResult {
        spacegroup_number: dataset_meta.spacegroup_number,
        hall_number: dataset_meta.hall_number,
        international_symbol: dataset_meta.international_symbol,
        choice: dataset_meta.choice,
        input_lattice,
        standardized_conventional_lattice: std_conv_lattice,
        standardized_primitive_lattice: std_prim_lattice,
        standardized_conventional_atoms: atoms_with_symbols(
            &std_conv_positions,
            &std_conv_types,
            &type_to_symbol,
        ),
        standardized_primitive_atoms: atoms_with_symbols(
            &std_prim_positions,
            &std_prim_types,
            &type_to_symbol,
        ),
        primitive_to_input_reciprocal: primitive_to_input,
        input_to_primitive_reciprocal: input_to_primitive,
        primitive_to_standardized_conventional_reciprocal: primitive_to_std_conv,
        standardized_conventional_to_primitive_reciprocal: std_conv_to_primitive,
        transformation_matrix: dataset_meta.transformation_matrix,
        origin_shift: dataset_meta.origin_shift,
    })
}
