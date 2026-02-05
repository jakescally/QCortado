//! Input file generation for Quantum ESPRESSO.
//!
//! Generates pw.x input files from QECalculation structs.

use super::types::*;
use std::fmt::Write;

/// Generates a complete pw.x input file from a calculation configuration.
pub fn generate_pw_input(calc: &QECalculation) -> String {
    let mut output = String::new();

    // &CONTROL namelist
    write_control_namelist(&mut output, calc);

    // &SYSTEM namelist
    write_system_namelist(&mut output, calc);

    // &ELECTRONS namelist
    write_electrons_namelist(&mut output, calc);

    // &IONS namelist (if relaxation)
    if matches!(calc.calculation, CalculationType::Relax | CalculationType::VcRelax | CalculationType::Md | CalculationType::VcMd) {
        write_ions_namelist(&mut output, calc);
    }

    // &CELL namelist (if variable cell)
    if matches!(calc.calculation, CalculationType::VcRelax | CalculationType::VcMd) {
        write_cell_namelist(&mut output);
    }

    // ATOMIC_SPECIES card
    write_atomic_species(&mut output, &calc.system);

    // ATOMIC_POSITIONS card
    write_atomic_positions(&mut output, &calc.system);

    // K_POINTS card
    write_kpoints(&mut output, &calc.kpoints);

    // CELL_PARAMETERS card (if ibrav=0)
    if matches!(calc.system.ibrav, BravaisLattice::Free) {
        write_cell_parameters(&mut output, &calc.system);
    }

    output
}

fn write_control_namelist(out: &mut String, calc: &QECalculation) {
    writeln!(out, "&CONTROL").unwrap();
    writeln!(out, "  calculation = '{}',", calc.calculation.as_str()).unwrap();
    writeln!(out, "  prefix = '{}',", calc.prefix).unwrap();
    writeln!(out, "  outdir = '{}',", calc.outdir).unwrap();
    writeln!(out, "  pseudo_dir = '{}',", calc.pseudo_dir).unwrap();

    if calc.tprnfor {
        writeln!(out, "  tprnfor = .true.,").unwrap();
    }
    if calc.tstress {
        writeln!(out, "  tstress = .true.,").unwrap();
    }

    if let Some(ref verbosity) = calc.verbosity {
        writeln!(out, "  verbosity = '{}',", verbosity).unwrap();
    }

    if let Some(etot_conv) = calc.etot_conv_thr {
        writeln!(out, "  etot_conv_thr = {:.2e},", etot_conv).unwrap();
    }
    if let Some(forc_conv) = calc.forc_conv_thr {
        writeln!(out, "  forc_conv_thr = {:.2e},", forc_conv).unwrap();
    }

    writeln!(out, "/\n").unwrap();
}

fn write_system_namelist(out: &mut String, calc: &QECalculation) {
    let sys = &calc.system;

    writeln!(out, "&SYSTEM").unwrap();
    writeln!(out, "  ibrav = {},", sys.ibrav as i32).unwrap();

    // Write celldm if provided (and ibrav != 0)
    if let Some(celldm) = &sys.celldm {
        if !matches!(sys.ibrav, BravaisLattice::Free) {
            // Only write non-zero celldm values
            for (i, &val) in celldm.iter().enumerate() {
                if val != 0.0 {
                    writeln!(out, "  celldm({}) = {},", i + 1, val).unwrap();
                }
            }
        }
    }

    writeln!(out, "  nat = {},", sys.atoms.len()).unwrap();
    writeln!(out, "  ntyp = {},", sys.species.len()).unwrap();
    writeln!(out, "  ecutwfc = {},", sys.ecutwfc).unwrap();

    if let Some(ecutrho) = sys.ecutrho {
        writeln!(out, "  ecutrho = {},", ecutrho).unwrap();
    }

    if sys.nspin != 1 {
        writeln!(out, "  nspin = {},", sys.nspin).unwrap();
    }

    if !matches!(sys.occupations, Occupations::Fixed) {
        writeln!(out, "  occupations = '{}',", match sys.occupations {
            Occupations::Fixed => "fixed",
            Occupations::Smearing => "smearing",
            Occupations::FromInput => "from_input",
            Occupations::Tetrahedra => "tetrahedra",
        }).unwrap();

        if matches!(sys.occupations, Occupations::Smearing) {
            writeln!(out, "  smearing = '{}',", sys.smearing.as_str()).unwrap();
            if let Some(degauss) = sys.degauss {
                writeln!(out, "  degauss = {},", degauss).unwrap();
            }
        }
    }

    writeln!(out, "/\n").unwrap();
}

fn write_electrons_namelist(out: &mut String, calc: &QECalculation) {
    writeln!(out, "&ELECTRONS").unwrap();
    writeln!(out, "  conv_thr = {:.2e},", calc.conv_thr).unwrap();
    writeln!(out, "  mixing_beta = {},", calc.mixing_beta).unwrap();
    writeln!(out, "/\n").unwrap();
}

fn write_ions_namelist(out: &mut String, _calc: &QECalculation) {
    writeln!(out, "&IONS").unwrap();
    writeln!(out, "  ion_dynamics = 'bfgs',").unwrap();
    writeln!(out, "/\n").unwrap();
}

fn write_cell_namelist(out: &mut String) {
    writeln!(out, "&CELL").unwrap();
    writeln!(out, "  cell_dynamics = 'bfgs',").unwrap();
    writeln!(out, "/\n").unwrap();
}

fn write_atomic_species(out: &mut String, sys: &QESystem) {
    writeln!(out, "ATOMIC_SPECIES").unwrap();
    for species in &sys.species {
        writeln!(out, "  {} {} {}", species.symbol, species.mass, species.pseudopotential).unwrap();
    }
    writeln!(out).unwrap();
}

fn write_atomic_positions(out: &mut String, sys: &QESystem) {
    writeln!(out, "ATOMIC_POSITIONS {{{}}}", sys.position_units.as_str()).unwrap();
    for atom in &sys.atoms {
        let pos = atom.position;
        // Check if any position is constrained
        let has_constraints = !atom.if_pos.iter().all(|&x| x);

        if has_constraints {
            let if_pos: Vec<i32> = atom.if_pos.iter().map(|&b| if b { 1 } else { 0 }).collect();
            writeln!(out, "  {} {:12.8} {:12.8} {:12.8}  {} {} {}",
                atom.symbol, pos[0], pos[1], pos[2],
                if_pos[0], if_pos[1], if_pos[2]
            ).unwrap();
        } else {
            writeln!(out, "  {} {:12.8} {:12.8} {:12.8}",
                atom.symbol, pos[0], pos[1], pos[2]
            ).unwrap();
        }
    }
    writeln!(out).unwrap();
}

fn write_kpoints(out: &mut String, kpoints: &KPoints) {
    match kpoints {
        KPoints::Gamma => {
            writeln!(out, "K_POINTS {{gamma}}").unwrap();
        }
        KPoints::Automatic { grid, offset } => {
            writeln!(out, "K_POINTS {{automatic}}").unwrap();
            writeln!(out, "  {} {} {}  {} {} {}",
                grid[0], grid[1], grid[2],
                offset[0], offset[1], offset[2]
            ).unwrap();
        }
        KPoints::Crystal { points } => {
            writeln!(out, "K_POINTS {{crystal}}").unwrap();
            writeln!(out, "  {}", points.len()).unwrap();
            for kp in points {
                writeln!(out, "  {:12.8} {:12.8} {:12.8} {:12.8}",
                    kp.k[0], kp.k[1], kp.k[2], kp.weight
                ).unwrap();
            }
        }
        KPoints::Tpiba { points } => {
            writeln!(out, "K_POINTS {{tpiba}}").unwrap();
            writeln!(out, "  {}", points.len()).unwrap();
            for kp in points {
                writeln!(out, "  {:12.8} {:12.8} {:12.8} {:12.8}",
                    kp.k[0], kp.k[1], kp.k[2], kp.weight
                ).unwrap();
            }
        }
        KPoints::CrystalB { path } => {
            writeln!(out, "K_POINTS {{crystal_b}}").unwrap();
            writeln!(out, "  {}", path.len()).unwrap();
            for point in path {
                writeln!(out, "  {:12.8} {:12.8} {:12.8} {}",
                    point.k[0], point.k[1], point.k[2], point.npoints
                ).unwrap();
            }
        }
        KPoints::TpibaB { path } => {
            writeln!(out, "K_POINTS {{tpiba_b}}").unwrap();
            writeln!(out, "  {}", path.len()).unwrap();
            for point in path {
                writeln!(out, "  {:12.8} {:12.8} {:12.8} {}",
                    point.k[0], point.k[1], point.k[2], point.npoints
                ).unwrap();
            }
        }
    }
    writeln!(out).unwrap();
}

fn write_cell_parameters(out: &mut String, sys: &QESystem) {
    if let Some(cell) = &sys.cell_parameters {
        let units = sys.cell_units.as_ref().unwrap_or(&PositionUnits::Angstrom);
        writeln!(out, "CELL_PARAMETERS {{{}}}", units.as_str()).unwrap();
        for row in cell {
            writeln!(out, "  {:12.8} {:12.8} {:12.8}", row[0], row[1], row[2]).unwrap();
        }
        writeln!(out).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_silicon_scf_input() {
        let calc = QECalculation {
            calculation: CalculationType::Scf,
            prefix: "silicon".to_string(),
            outdir: "./tmp".to_string(),
            pseudo_dir: "./pseudo".to_string(),
            system: QESystem {
                ibrav: BravaisLattice::CubicF,
                celldm: Some([10.2, 0.0, 0.0, 0.0, 0.0, 0.0]),
                cell_parameters: None,
                cell_units: None,
                species: vec![AtomicSpecies {
                    symbol: "Si".to_string(),
                    mass: 28.086,
                    pseudopotential: "Si.pz-vbc.UPF".to_string(),
                }],
                atoms: vec![
                    Atom {
                        symbol: "Si".to_string(),
                        position: [0.0, 0.0, 0.0],
                        if_pos: [true, true, true],
                    },
                    Atom {
                        symbol: "Si".to_string(),
                        position: [0.25, 0.25, 0.25],
                        if_pos: [true, true, true],
                    },
                ],
                position_units: PositionUnits::Crystal,
                ecutwfc: 20.0,
                ecutrho: Some(80.0),
                nspin: 1,
                occupations: Occupations::Fixed,
                smearing: SmearingType::default(),
                degauss: None,
            },
            kpoints: KPoints::Automatic {
                grid: [4, 4, 4],
                offset: [1, 1, 1],
            },
            conv_thr: 1.0e-8,
            mixing_beta: 0.7,
            tprnfor: false,
            tstress: false,
            forc_conv_thr: None,
            etot_conv_thr: None,
            verbosity: None,
        };

        let input = generate_pw_input(&calc);

        // Basic sanity checks
        assert!(input.contains("calculation = 'scf'"));
        assert!(input.contains("prefix = 'silicon'"));
        assert!(input.contains("ibrav = 2"));
        assert!(input.contains("ecutwfc = 20"));
        assert!(input.contains("Si 28.086 Si.pz-vbc.UPF"));
        assert!(input.contains("K_POINTS {automatic}"));
        assert!(input.contains("4 4 4  1 1 1"));
    }
}
