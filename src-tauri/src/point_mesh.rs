//! Preflight reciprocal mesh estimates, without running an SCF calculation.
use crate::symmetry::SymmetryAnalyzeInput;
use serde::{Deserialize, Serialize};
use spglib_sys as ffi;
use std::collections::HashMap;

#[derive(Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PointMeshInput {
    pub structure: SymmetryAnalyzeInput,
    pub mesh: [i32; 3],
    pub shift: [i32; 3],
    pub time_reversal: bool,
    pub use_symmetry: bool,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PointMeshResult {
    pub full_points: usize,
    pub irreducible_points: usize,
}

pub fn estimate(input: PointMeshInput) -> Result<PointMeshResult, String> {
    if input.mesh.iter().any(|n| *n < 1) || input.shift.iter().any(|s| *s != 0 && *s != 1) {
        return Err("Mesh dimensions must be positive integers and offsets must be 0 or 1.".into());
    }
    let full_points = input
        .mesh
        .iter()
        .try_fold(1usize, |n, m| n.checked_mul(*m as usize))
        .filter(|n| *n <= 2_000_000)
        .ok_or("Point-count preview supports meshes up to 2,000,000 points.")?;
    let cell = input.structure;
    let a = cell.lattice;
    let determinant = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if !a.iter().flatten().all(|v| v.is_finite())
        || determinant.abs() < 1e-12
        || cell.atoms.is_empty()
        || cell.atoms.len() > 100_000
        || !cell.symprec.is_finite()
        || cell.symprec <= 0.0
        || cell.atoms.iter().any(|atom| {
            atom.symbol.trim().is_empty() || atom.position.iter().any(|v| !v.is_finite())
        })
    {
        return Err(
            "A finite, nonsingular cell and fractional atomic positions are required.".into(),
        );
    }
    if !input.use_symmetry {
        // Only k -> -k remains. Fixed points obey 2*i + shift = 0 mod mesh.
        let fixed: usize = input
            .mesh
            .iter()
            .zip(input.shift)
            .map(|(n, shift)| {
                if n % 2 == 1 {
                    1
                } else if shift == 0 {
                    2
                } else {
                    0
                }
            })
            .product();
        return Ok(PointMeshResult {
            full_points,
            irreducible_points: if input.time_reversal {
                (full_points + fixed) / 2
            } else {
                full_points
            },
        });
    }
    let mut lattice = [
        [a[0][0], a[1][0], a[2][0]],
        [a[0][1], a[1][1], a[2][1]],
        [a[0][2], a[1][2], a[2][2]],
    ];
    let mut positions: Vec<[f64; 3]> = cell
        .atoms
        .iter()
        .map(|atom| atom.position.map(|v| v.rem_euclid(1.0)))
        .collect();
    let mut labels = HashMap::new();
    let types: Vec<i32> = cell
        .atoms
        .iter()
        .map(|atom| {
            let next = labels.len() as i32 + 1;
            *labels.entry(atom.symbol.clone()).or_insert(next)
        })
        .collect();
    let mut addresses = vec![[0i32; 3]; full_points];
    let mut mapping = vec![0i32; full_points];
    let count = unsafe {
        ffi::spg_get_ir_reciprocal_mesh(
            addresses.as_mut_ptr(),
            mapping.as_mut_ptr(),
            input.mesh.as_ptr(),
            input.shift.as_ptr(),
            i32::from(input.time_reversal),
            lattice.as_mut_ptr(),
            positions.as_mut_ptr(),
            types.as_ptr(),
            positions.len() as i32,
            cell.symprec,
        )
    };
    if count <= 0 || count as usize > full_points {
        return Err("Symmetry reduction of this mesh failed.".into());
    }
    Ok(PointMeshResult {
        full_points,
        irreducible_points: count as usize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symmetry::SymmetryAtomInput;
    fn cubic(mesh: [i32; 3], shift: [i32; 3]) -> PointMeshInput {
        PointMeshInput {
            structure: SymmetryAnalyzeInput {
                lattice: [[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]],
                atoms: vec![SymmetryAtomInput {
                    symbol: "Si".into(),
                    position: [0.; 3],
                }],
                symprec: 1e-5,
                angle_tolerance: -1.,
            },
            mesh,
            shift,
            time_reversal: true,
            use_symmetry: true,
        }
    }
    #[test]
    fn cubic_mesh_counts_and_shifts() {
        assert_eq!(
            estimate(cubic([1; 3], [0; 3])).unwrap().irreducible_points,
            1
        );
        assert_eq!(
            estimate(cubic([4; 3], [0; 3])).unwrap().irreducible_points,
            10
        );
        assert_eq!(
            estimate(cubic([4; 3], [1; 3])).unwrap().irreducible_points,
            4
        );
    }
    #[test]
    fn symmetry_and_time_reversal_are_independent() {
        let mut input = cubic([4; 3], [0; 3]);
        input.use_symmetry = false;
        assert_eq!(estimate(input.clone()).unwrap().irreducible_points, 36);
        input.time_reversal = false;
        assert_eq!(estimate(input).unwrap().irreducible_points, 64);
    }
    #[test]
    fn invalid_and_excessive_meshes_are_rejected() {
        assert!(estimate(cubic([0, 4, 4], [0; 3])).is_err());
        assert!(estimate(cubic([i32::MAX; 3], [0; 3])).is_err());
        assert!(estimate(cubic([4; 3], [2; 3])).is_err());
        let mut input = cubic([4; 3], [0; 3]);
        input.structure.lattice = [[0.; 3]; 3];
        assert!(estimate(input).is_err());
    }
}
