# Engine Boundaries

## Boundary Rule

The Cortado platform owns orchestration, storage, execution infrastructure, and shared viewers. Engines own the scientific input model, executable workflow, file lifecycle, and native parsers.

Engine inputs stay engine-specific. Viewer outputs may be normalized.

## What Belongs To The Platform

Platform modules can depend on stable engine identifiers and normalized datasets. They should not depend on QE namelist fields, QE pseudopotentials, Wien2k case file details, or executable-specific parser assumptions.

Platform candidates in the current repo:

- App shell:
  - `src/App.tsx`
  - `src/ViewerApp.tsx`
- Project shell:
  - `src/components/ProjectBrowser.tsx`
  - `src/components/ProjectDashboard.tsx`
  - `src-tauri/src/projects.rs`
- Tasks and live output:
  - `src/lib/TaskContext.tsx`
  - `src/components/TaskQueuePage.tsx`
  - `src/components/LiveOutputPanel.tsx`
  - `src-tauri/src/process_manager.rs`
- HPC shell:
  - `src/components/HpcSetupWizard.tsx`
  - `src/components/HpcRunSettings.tsx`
  - `src/components/HpcActivityPanel.tsx`
  - `src/components/HpcNodeActivityPage.tsx`
  - `src-tauri/src/hpc/ssh.rs`
  - `src-tauri/src/hpc/slurm.rs`
  - `src-tauri/src/hpc/runner.rs`
  - `src-tauri/src/hpc/sync.rs`
  - `src-tauri/src/hpc/utilization.rs`
  - `src-tauri/src/hpc/cluster_snapshot.rs`
- Structure and viewer foundations:
  - `src/lib/cifParser.ts`
  - `src/lib/symmetry.ts`
  - `src/lib/symmetryTransform.ts`
  - `src/lib/reciprocalLattice.ts`
  - `src/lib/brillouinZone.ts`
  - `src/lib/brillouinZoneData.ts`
  - `src/components/UnitCellViewer.tsx`
  - `src/components/BrillouinZoneViewer.tsx`
- Normalized viewers:
  - `src/components/BandPlot.tsx`
  - `src/components/ElectronicDOSPlot.tsx`
  - `src/components/TransportPlot.tsx`
  - `src/components/PhononPlot.tsx`, after its dataset contract is normalized.

## What Belongs To The QE Engine

QE code should own all QE-specific scientific behavior. These concepts should not move into generic platform types unless there is a proven cross-engine equivalent.

Current QE-owned candidates:

- Backend engine module:
  - `src-tauri/src/engines/qe/input.rs`
  - `src-tauri/src/engines/qe/output.rs`
  - `src-tauri/src/engines/qe/runner.rs`
  - `src-tauri/src/engines/qe/types.rs`
  - `src-tauri/src/engines/qe/pseudopotentials.rs`
  - `src-tauri/src/engines/qe/bands.rs`
  - `src-tauri/src/engines/qe/phonon.rs`
  - `src-tauri/src/engines/qe/hubbard.rs`
  - `src-tauri/src/engines/qe/wannier.rs`
  - `src-tauri/src/engines/qe/epw.rs`
  - `src-tauri/src/engines/qe/transport.rs`
  - `src-tauri/src/qe/mod.rs` is a compatibility shim only.
- Frontend QE workflow UI:
  - `src/components/qe/index.ts` is the explicit app import boundary.
  - `src/components/SCFWizard.tsx`
  - `src/components/BandStructureWizard.tsx`
  - `src/components/ElectronicDOSWizard.tsx`
  - `src/components/FermiSurfaceWizard.tsx`
  - `src/components/HubbardLrtWizard.tsx`
  - `src/components/PhononWizard.tsx`
  - `src/components/WannierWizard.tsx`
  - `src/components/TransportWizard.tsx`
  - `src/components/EpwWizard.tsx`
- Frontend QE helpers:
  - `src/lib/engines/qe/progress.ts`
  - `src/lib/engines/qe/bravaisInference.ts`
  - `src/lib/engines/qe/pseudopotentialCutoffs.ts`
  - `src/lib/engines/qe/pseudopotentialMetadataCache.ts`
  - `src/lib/engines/qe/hubbard.ts`
  - `src/lib/engines/qe/wannierQuality.ts`
  - `src/lib/engines/qe/epw.ts`
  - `src/lib/engines/qe/scfRunSettingsClipboard.ts`
  - `src/lib/engines/qe/scfSorting.ts`
  - `src/lib/engines/qe/phononReady.ts`
  - `src/lib/engines/qe/hpc.ts`
  - `src/lib/engines/qe/hpcProfiles.ts`
  - `src/lib/engines/qe/bandDatasetAdapter.ts`
  - Top-level files such as `src/lib/qeProgress.ts`, `src/lib/hubbard.ts`, and `src/lib/epw.ts` are compatibility shims only.

QE-specific concepts include:

- `pw.x`, `bands.x`, `projwfc.x`, `dos.x`, `ph.x`, `q2r.x`, `matdyn.x`, `hp.x`, `epw.x`.
- QE namelists and cards.
- QE `prefix`, `outdir`, `pseudo_dir`.
- QE `.save` directory lifecycle.
- Pseudopotentials, UPF parsing, SSSP metadata, and cutoff derivation.
- QE smearing names.
- QE Hubbard card syntax and `hp.x` output.
- QE-specific output recovery.

## What Belongs To Future Wien2k Engine Code

Wien2k is not implemented. These concepts should be reserved for a future engine module, not added to platform types or QE types.

- Remote-only execution.
- Case directory lifecycle.
- `case.struct`.
- RMT, RKmax, Gmax, lmax.
- `init_lapw`.
- `x nn`, `x sgroup`, `x symmetry`, `lstart`, `kgen`, `dstart`.
- `run_lapw`, `runsp_lapw`.
- `case.scf`.
- `lapw1`, `lapw2`, `spaghetti`.

Wien2k does not use pseudopotentials. Do not require pseudopotential fields, UPF metadata, SSSP data, `ecutwfc`, or `ecutrho` for Wien2k flows.

## Engine Interface Shape

The first engine boundary should be intentionally small. It should describe orchestration without pretending every engine has the same physics input.

Current frontend shape:

```ts
type EngineId = "qe" | "wien2k"; // "wien2k" is reserved, not implemented.

interface EngineDescriptor {
  id: EngineId;
  label: string;
  executionModes: Array<"local" | "hpc">;
  taskKinds: string[];
}
```

Future-compatible shape after Wien2k is introduced:

```ts
type EngineId = "qe" | "wien2k";

interface EngineDescriptor {
  id: EngineId;
  label: string;
  executionModes: Array<"local" | "hpc">;
  taskKinds: string[];
}
```

Current backend engine identity:

```rust
pub enum EngineId {
    Qe,
    // Reserved for future work. This does not implement Wien2k.
    Wien2k,
}
```

Suggested future backend orchestration shape:

```rust

pub struct EngineJobBundle {
    pub engine_id: EngineId,
    pub task_kind: String,
    pub label: String,
    pub input_files: Vec<EngineInputFile>,
    pub command_lines: Vec<String>,
    pub artifact_policy: ArtifactPolicy,
}
```

The job bundle can be shared. The task config that creates the bundle should remain engine-specific.

## Boundary Contracts

### Inputs

Inputs are not shared across engines.

- QE SCF input remains `QECalculation`.
- QE bands, DOS, phonon, Wannier, EPW, and transport configs remain QE-owned.
- Future Wien2k configs should model Wien2k case and initialization concepts directly.

Do not create:

```ts
interface UniversalScfConfig {
  ecutwfc?: number;
  pseudo_dir?: string;
  rkmax?: number;
  rmt?: Record<string, number>;
}
```

That shape mixes unrelated engines and makes both engines harder to validate.

### Outputs

Outputs may be normalized when a shared viewer benefits.

Shared candidates:

- Electronic band path dataset.
- Electronic DOS dataset.
- Phonon dispersion dataset.
- Phonon DOS dataset.
- Transport tensor dataset.
- Generic numeric tables.

Engine-native result payloads should still be stored when needed for provenance and engine-specific views.

### Projects

Project records should evolve toward:

```ts
interface CalculationRun {
  id: string;
  engine_id: EngineId;
  calc_type: string;
  parameters: unknown;
  result: unknown;
}
```

During migration, legacy records without `engine_id` should be treated as QE.

### HPC Profiles

Separate cluster access from engine toolchains.

Shared profile fields:

- host, port, username, auth method.
- remote workspace root.
- remote project root.
- Slurm defaults.
- launcher settings.

QE toolchain fields:

- remote QE bin directory.
- remote CPU/GPU QE bin directories.
- remote pseudopotential directories.
- remote EPW path.
- remote Wannier90 and postw90 paths.

Future Wien2k toolchain fields:

- remote Wien2k root or command environment.
- default scratch/case root.
- module preamble or profile-specific setup.

## Anti-Patterns To Avoid

- Moving QE files into vague generic names without a real abstraction.
- Adding Wien2k fields to QE config objects.
- Adding pseudopotential fields to platform-level structure or task configs.
- Treating `scf` as a universal task with one universal payload.
- Making viewers parse engine-native files directly.
- Renaming Tauri commands in the same PR that moves backend behavior.
- Refactoring `src-tauri/src/lib.rs`, project storage, and frontend task queue in one PR.

## PR-Sized Boundary Extraction Order

1. Add docs and a shared vocabulary.
2. Add `engine_id` metadata with QE legacy defaults.
3. Introduce normalized viewer dataset types and QE adapters.
4. Move QE-only frontend helpers behind compatibility re-exports.
5. Move backend pseudopotential code under QE ownership.
6. Extract QE task orchestration from `src-tauri/src/lib.rs` one workflow at a time.
7. Split HPC profile shell from QE toolchain fields.
8. Add a shared engine job bundle interface and adapt one QE workflow.

## Validation

Frontend-only boundary changes:

```bash
npm run build
npm run test:unit
```

Backend-only boundary changes:

```bash
cargo check --manifest-path src-tauri/Cargo.toml
```

Cross-boundary changes:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```
