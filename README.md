# QCortado

A modern, user-friendly interface for [Quantum ESPRESSO](https://www.quantum-espresso.org/) built with Tauri, Rust, and React.

## Overview

QCortado simplifies working with Quantum ESPRESSO by providing:

- **Guided workflows** for common calculations (SCF, bands, DOS, relaxation)
- **Smart input generation** with validation and sensible defaults
- **Real-time execution monitoring** with streaming output
- **Structured result parsing** - automatically extract energies, forces, and band structures
- **Project management** to organize related calculations

## Architecture

```
QCortado/
├── src/                    # React frontend
│   ├── App.tsx            # Main application component
│   └── App.css            # Styling
├── src-tauri/             # Rust backend
│   ├── src/
│   │   ├── lib.rs         # Tauri commands and app state
│   │   ├── qe/
│   │   │   ├── types.rs   # Core QE data structures
│   │   │   ├── input.rs   # Input file generation
│   │   │   ├── output.rs  # Output parsing
│   │   │   └── runner.rs  # Process execution
│   └── Cargo.toml
└── package.json
```

## Requirements

- **Quantum ESPRESSO 7.x** - compiled and installed
- **Node.js 20+**
- **Rust 1.70+**
- **Cargo**

## Development

```bash
# Install frontend dependencies
npm install

# Run in development mode
npm run tauri dev

# Build for production
npm run tauri build
```

## Tauri Commands (Rust → Frontend)

| Command | Description |
|---------|-------------|
| `set_qe_path` | Configure path to QE bin directory |
| `get_qe_path` | Get current QE path |
| `check_qe_executables` | List available QE programs |
| `generate_input` | Generate pw.x input from config |
| `validate_calculation` | Validate calculation before running |
| `run_calculation` | Execute pw.x calculation |
| `parse_output` | Parse QE output text |
| `list_pseudopotentials` | List .UPF files in a directory |

## QE Workflow Support

### Currently Implemented
- [x] pw.x input generation (SCF, NSCF, bands, relax, vc-relax)
- [x] Output parsing (energy, forces, stress, convergence)
- [x] Process execution with async support

### Planned
- [ ] Band structure workflow (SCF → NSCF → bands.x)
- [ ] DOS workflow (SCF → NSCF → dos.x)
- [ ] Phonon workflow (SCF → ph.x → dynmat.x)
- [ ] Pseudopotential management
- [ ] Structure visualization
- [ ] Result plotting

## Input File Generation

QCortado generates QE input files from strongly-typed Rust structures:

```rust
let calc = QECalculation {
    calculation: CalculationType::Scf,
    prefix: "silicon".to_string(),
    system: QESystem {
        ibrav: BravaisLattice::CubicF,
        celldm: Some([10.2, 0.0, 0.0, 0.0, 0.0, 0.0]),
        species: vec![AtomicSpecies {
            symbol: "Si".to_string(),
            mass: 28.086,
            pseudopotential: "Si.pz-vbc.UPF".to_string(),
        }],
        atoms: vec![
            Atom { symbol: "Si".to_string(), position: [0.0, 0.0, 0.0], .. },
            Atom { symbol: "Si".to_string(), position: [0.25, 0.25, 0.25], .. },
        ],
        ecutwfc: 20.0,
        ..
    },
    kpoints: KPoints::Automatic { grid: [4, 4, 4], offset: [1, 1, 1] },
    ..
};

let input_text = generate_pw_input(&calc);
```

## Output Parsing

Results are automatically parsed into structured data:

```rust
let result: QEResult = parse_pw_output(&output_text);

println!("Converged: {}", result.converged);
println!("Total energy: {} Ry", result.total_energy.unwrap());
println!("SCF iterations: {}", result.n_scf_steps.unwrap());
```

## License

MIT
