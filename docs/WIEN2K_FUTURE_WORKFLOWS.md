# Future Wien2k Workflows

## Status

QCortado implements one WIEN2k-owned workflow: remote `case.struct` source
setup (`engine_setup`). Once a remote WIEN2k installation is verified, a
project with WIEN2k selected exposes a `WIEN2k Structure` tile.

SCF initialization, SCF runs, bands, and DOS are not implemented. Their
descriptors remain reserved and are not exposed through the selectable
frontend workflow registry.

## Core Difference From QE

QE uses pseudopotentials. Wien2k does not.

This has direct architecture consequences:

- Do not add pseudopotential fields to shared platform configs.
- Do not require `pseudo_dir`, UPF metadata, SSSP data, `ecutwfc`, or `ecutrho` in future Wien2k configs.
- Do not create a universal SCF config that contains both QE and Wien2k fields.
- Keep QE pseudopotential discovery and repair tools inside the QE engine.

## Execution Model

The WIEN2k boundary is remote-only.

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

WIEN2k engine code owns:

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

## Implemented Workflow: Structure Source Setup

The structure-source miniapp deliberately precedes SCF initialization:

1. Start from the selected project's CIF data.
2. Standardize the conventional cell through shared spglib analysis inside
   the WIEN2k engine conversion boundary.
3. Write an initial fixed-format `case.struct` locally. QCortado does not call
   WIEN2k `cif2struct`.
4. Stage the draft in a transient remote directory below
   `{remote_workspace_root}/qcortado/{project_id}/wien2k/{session_id}`.
5. Run native refinement in explicit review stages:
   `setrmt_lapw` plus `x nn`, then `x sgroup -settol`, then `x symmetry`.
6. Stream native output inline and require user approval of the final
   symmetry candidate.
7. Save an `engine_setup` calculation entry only after approval.

The accepted local `<case>.struct` and its small diagnostic logs are the
canonical project source. The transient remote directory is cleaned on save
or discard when the remote profile is reachable; cleanup failure does not
invalidate an accepted local source. A future SCF workflow must restage the
accepted structure into its own remote case directory rather than depend on
the refinement session.

Platform-owned pieces:

- Project and CIF selection.
- Remote SSH transport and inline output streaming.
- Live output.
- Artifact sync.
- Saved calculation record.

WIEN2k-owned pieces:

- Struct generation.
- Structure refinement command sequence and native artifacts.
- Input prompts or scripted choices.
- Error parsing.
- Case structure artifacts.

## Future Workflow: SCF Initialization

Future SCF initialization starts from an accepted saved structure source. It
may run `lstart`, `kgen`, and `dstart` as required by the chosen WIEN2k SCF
contract. Those operations are intentionally not part of the implemented
structure miniapp.

## Candidate Workflow: SCF

Future flow:

1. Select an accepted WIEN2k structure source and initialize a new SCF case.
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
src/components/wien2k/
  Wien2kStructureWizard.tsx            # implemented engine_setup UI
src/lib/engines/wien2k/
  plugin.ts                             # structure route plus reserved contract
  structure.ts                          # implemented structure types/API
  caseState.ts                          # future case lifecycle support
  types.ts                              # future native workflow types
```

Backend:

```text
src-tauri/src/engines/wien2k/
  mod.rs
  structure.rs                          # implemented struct/session/stage ownership
  types.rs
  init.rs
  scf.rs
  bands.rs
  dos.rs
  remote.rs
  output.rs
  adapters.rs
```

Future initialization/runner/parser files should be added only as those
workflows are built.

## Dependencies On Current Migration

The structure workflow relies on the following completed platform work:

1. `engine_id` on saved calculations.
2. Project storage that can distinguish QE calculations from future engine calculations.
3. HPC profile editing that keeps shared SSH/Slurm roots alongside
   engine-specific QE and verified WIEN2k runtime paths.
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

## Non-Goals For The Structure Workflow

- Do not implement WIEN2k SCF initialization or calculation runners.
- Do not add Wien2k settings to the current QE setup UI.
- Do not add placeholder Wien2k commands that cannot run.
- Do not change QE calculation behavior.
- Do not remove QE pseudopotential features.
- Do not rewrite project storage in a way that breaks existing QE records.
