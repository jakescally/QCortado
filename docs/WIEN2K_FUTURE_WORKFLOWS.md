# Future Wien2k Workflows

## Status

Wien2k is not implemented in QCortado.

This document records future workflow boundaries so the current QE-to-engine migration does not block a later Wien2k backend. The immediate work is to isolate QE behavior behind explicit engine boundaries while preserving existing QE behavior.

## Core Difference From QE

QE uses pseudopotentials. Wien2k does not.

This has direct architecture consequences:

- Do not add pseudopotential fields to shared platform configs.
- Do not require `pseudo_dir`, UPF metadata, SSSP data, `ecutwfc`, or `ecutrho` in future Wien2k configs.
- Do not create a universal SCF config that contains both QE and Wien2k fields.
- Keep QE pseudopotential discovery and repair tools inside the QE engine.

## Execution Model

The planned Wien2k backend is remote-only.

Platform requirements:

- SSH connection.
- Slurm submission.
- Remote workspace root.
- Remote project archive root.
- Remote artifact sync.
- Live output and scheduler state polling.

Engine requirements:

- Remote Wien2k environment validation.
- Case directory staging.
- Case-specific command sequence.
- Case-specific parser and artifact policy.

The platform should be able to say that an engine supports only HPC execution:

```ts
interface EngineDescriptor {
  id: "wien2k";
  label: "Wien2k";
  executionModes: ["hpc"];
}
```

Do not expose local execution controls for a future Wien2k workflow unless local support is explicitly designed later.

## Wien2k Concepts To Keep Engine-Specific

Future Wien2k engine code should own:

- `case.struct`.
- Case directory lifecycle.
- RMT values.
- RKmax, Gmax, lmax.
- `init_lapw`.
- `x nn`.
- `x sgroup`.
- `x symmetry`.
- `lstart`.
- `kgen`.
- `dstart`.
- `run_lapw`.
- `runsp_lapw`.
- `case.scf` parsing.
- `lapw1`, `lapw2`, and `spaghetti`.
- Wien2k-specific convergence and error detection.

These should not become platform concepts unless there is a proven equivalent across engines.

## Candidate Workflow: Initialization

Future flow:

1. Start from a project CIF structure.
2. Convert structure to a Wien2k `case.struct` candidate.
3. Stage a remote case directory.
4. Run initialization sequence:
   - `x nn`
   - `x sgroup`
   - `x symmetry`
   - `lstart`
   - `kgen`
   - `dstart`
5. Parse initialization outputs.
6. Save an initialized case calculation entry.

Platform-owned pieces:

- Project and CIF selection.
- Remote job submission.
- Live output.
- Artifact sync.
- Saved calculation record.

Wien2k-owned pieces:

- Struct generation.
- Initialization command sequence.
- Input prompts or scripted choices.
- Error parsing.
- Case artifacts.

## Candidate Workflow: SCF

Future flow:

1. Select an initialized Wien2k case.
2. Choose spin mode and convergence controls.
3. Submit remote `run_lapw` or `runsp_lapw`.
4. Parse `case.scf`.
5. Save native result and normalized summary.

Possible normalized outputs:

- Total energy.
- Fermi energy.
- Convergence status.
- Iteration count.
- Magnetic moments, if available and requested.

Wien2k-specific inputs should remain separate from QE SCF inputs.

## Candidate Workflow: Bands

Future flow:

1. Select a converged Wien2k SCF case.
2. Prepare k-path input for Wien2k.
3. Run the required `lapw1`/`lapw2`/`spaghetti` sequence.
4. Parse bands output.
5. Adapt to `BandPathDataset`.

Shared viewer output:

- `BandPathDataset`.

Wien2k-owned details:

- Which files must exist in the case directory.
- How the k-path is encoded.
- Which commands are required for spin-polarized or spin-orbit cases.
- How `spaghetti` output maps to path coordinates and labels.

## Candidate Workflow: DOS

Future flow:

1. Select a converged Wien2k SCF case.
2. Run the future Wien2k DOS command sequence.
3. Parse DOS files.
4. Adapt to `DosDataset`.

Shared viewer output:

- `DosDataset`.

Wien2k-owned details:

- Command sequence.
- File names.
- Orbital/projected DOS format.
- Energy reference handling.

## Future Module Layout

Frontend:

```text
src/engines/wien2k/
  components/
    Wien2kInitWizard.tsx
    Wien2kScfWizard.tsx
    Wien2kBandsWizard.tsx
    Wien2kDosWizard.tsx
  lib/
    wien2kProgress.ts
    wien2kCaseSummary.ts
  types/
    index.ts
  adapters/
    viewerDatasets.ts
```

Backend:

```text
src-tauri/src/engines/wien2k/
  mod.rs
  types.rs
  struct_file.rs
  init.rs
  scf.rs
  bands.rs
  dos.rs
  remote.rs
  output.rs
  adapters.rs
```

This layout should be created only when Wien2k implementation work begins.

## Dependencies On Current Migration

Before Wien2k work starts, the QE migration should provide:

1. `engine_id` on saved calculations.
2. Project storage that can distinguish QE calculations from future engine calculations.
3. Shared HPC shell separated from QE toolchain fields.
4. Engine-specific remote job bundle generation.
5. Shared normalized band and DOS datasets.
6. QE pseudopotentials isolated inside the QE engine.
7. Compatibility path for old QE projects.

## Future Validation Expectations

Wien2k implementation PRs will need more than current QE validation. At minimum:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

Future backend parser modules should add Rust unit tests for representative Wien2k outputs.

Remote workflow PRs should also include manual validation notes for:

- remote case directory creation.
- Slurm submission.
- scheduler polling.
- artifact sync.
- project save and reload.

## Non-Goals For The Current Migration

- Do not implement Wien2k.
- Do not add Wien2k settings to the current QE setup UI.
- Do not add placeholder Wien2k commands that cannot run.
- Do not change QE calculation behavior.
- Do not remove QE pseudopotential features.
- Do not rewrite project storage in a way that breaks existing QE records.

