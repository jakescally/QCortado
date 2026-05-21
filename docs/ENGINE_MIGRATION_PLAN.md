# Engine Migration Plan

## Purpose

QCortado currently behaves as a Quantum ESPRESSO desktop app. The migration target is an engine-based Cortado platform where shared application infrastructure is separated from engine-specific scientific workflows.

This plan is for boundary extraction only. It does not implement Wien2k, change calculation behavior, or replace QE inputs with a generic universal SCF model.

## Migration Goals

- Preserve all current QE behavior.
- Make QE an explicit engine implementation.
- Keep engine-specific inputs engine-specific.
- Normalize result datasets only where shared viewers benefit.
- Keep platform features reusable by future engines.
- Make future remote-only Wien2k support possible without mixing QE pseudopotential concepts into Wien2k code.

## Current Architecture Map

### Frontend

- `src/App.tsx`
  - Main application shell.
  - Owns routing, settings, QE executable paths, HPC profile state, recovery flows, viewer state, and wizard handoff.
- `src/ViewerApp.tsx`
  - Read-only viewer application shell.
  - Syncs remote project library and reuses project/viewer components.
- `src/components`
  - Mixes shared UI and engine-specific workflows.
  - Shared candidates: `ProjectBrowser`, `ProjectDashboard`, `TaskQueuePage`, `LiveOutputPanel`, `StorageManagerPage`, HPC pages, `UnitCellViewer`, `BrillouinZoneViewer`.
  - QE-specific workflows: `SCFWizard`, `BandStructureWizard`, `ElectronicDOSWizard`, `FermiSurfaceWizard`, `HubbardLrtWizard`, `PhononWizard`, `WannierWizard`, `TransportWizard`, `EpwWizard`.
  - Viewer candidates: `BandPlot`, `ElectronicDOSPlot`, `PhononPlot`, `TransportPlot`, `EpwViewer`.
- `src/lib`
  - Mixes platform utilities, viewer math, and QE-specific helpers.
  - Shared candidates: CIF parsing, symmetry transforms, reciprocal lattice, k-path transforms, live output, task context, theme, storage helpers.
  - QE-specific candidates: `qeProgress`, `qeBravaisInference`, `pseudopotentialCutoffs`, `pseudopotentialMetadataCache`, `hubbard`, `wannierQuality`, `epw`, QE run settings clipboard.

### Backend

- `src-tauri/src/lib.rs`
  - Main command registration and application state.
  - Current hotspot for QE settings, pseudopotential parsing, HPC commands, task orchestration, and utility helpers.
- `src-tauri/src/qe`
  - Existing QE module island.
  - Contains QE types, input generation, output parsing, runner, bands, phonons, Hubbard, Wannier, EPW, and transport logic.
- `src-tauri/src/hpc`
  - Shared SSH, Slurm, utilization, cluster snapshots, sync, and viewer library primitives.
  - Still contains QE-shaped profile fields and validation.
- `src-tauri/src/projects.rs`
  - Project storage and archive management.
  - Current saved calculation result envelope is `QEResult`.
- `src-tauri/src/config.rs`
  - Application config.
  - Current config contains QE defaults and QE executable paths.

## Target Architecture

The target shape is a platform shell plus engine modules.

```text
src/
  platform/
    app/
    projects/
    tasks/
    hpc/
    storage/
    viewer/
  results/
    bands/
    dos/
    phonons/
    transport/
    tables/
  engines/
    qe/
      components/
      lib/
      types/
      adapters/
    wien2k/
      types/
      remote/

src-tauri/src/
  platform/
    config.rs
    projects.rs
    tasks.rs
    hpc/
    storage.rs
    viewer.rs
  results/
    mod.rs
    bands.rs
    dos.rs
    phonons.rs
    transport.rs
    tables.rs
  engines/
    qe/
      mod.rs
      types.rs
      input.rs
      output.rs
      tasks.rs
      pseudopotentials.rs
      bands.rs
      phonon.rs
      hubbard.rs
      wannier.rs
      epw.rs
      transport.rs
    wien2k/
      mod.rs
      types.rs
      remote.rs
```

This is a target layout, not a single PR. Early PRs should use re-export shims so existing imports and Tauri commands keep working.

## Ownership Boundaries

### Shared platform code

Platform code may know that engines exist, but should not know QE namelists or Wien2k case files.

- Project browser and project metadata.
- Task queue and live output.
- Local app settings shell.
- HPC profile shell, SSH transport, Slurm submission, scheduler polling, and artifact sync.
- CIF parsing and crystal/unit-cell display.
- Brillouin-zone and k-path UI primitives.
- Normalized viewer datasets.

### QE engine code

QE code owns QE inputs, QE outputs, QE executable names, QE pseudopotentials, and QE artifact rules.

- `pw.x`, `bands.x`, `projwfc.x`, `dos.x`, `ph.x`, `q2r.x`, `matdyn.x`, `hp.x`, `epw.x`.
- QE namelists and cards.
- `prefix`, `outdir`, `pseudo_dir`, `.save` directories.
- Pseudopotentials, UPF parsing, SSSP metadata, cutoffs.
- QE smearing names and Hubbard card.
- QE parser behavior and QE recovery behavior.

### Future Wien2k engine code

Wien2k support is future work. It should be added as a remote-only engine after QE boundaries are explicit.

- Case directory lifecycle.
- `case.struct`.
- RMT, RKmax, Gmax, lmax.
- `init_lapw` flow.
- `x nn`, `x sgroup`, `x symmetry`, `lstart`, `kgen`, `dstart`.
- `run_lapw`, `runsp_lapw`.
- `case.scf` parsing.
- `lapw1`, `lapw2`, `spaghetti` workflows.

Wien2k must not receive QE pseudopotential fields.

## Phased Migration Steps

### Phase 1: Document and label existing behavior

PR size: docs plus minimal type annotations only.

- Add docs describing platform, engine, and dataset boundaries.
- Add an `engine_id` concept to new documentation and future type names.
- Do not change runtime behavior.

Validation:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

For docs-only changes, code validation is optional unless application files are touched.

### Phase 2: Introduce stable engine identifiers

PR size: type additions and legacy-safe defaults.

- Add frontend `EngineId = "qe"` type.
- Add backend `EngineId` or string-backed equivalent.
- Default existing saved calculations to QE when missing `engine_id`.
- Do not alter calculation payload shapes.

Validation:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

### Phase 3: Create result dataset adapters

PR size: viewer type extraction and adapters.

- Define normalized result dataset types under `src/results`.
- Add QE adapters from current `BandData`, DOS, phonon, transport, and EPW table payloads.
- Keep viewer props backward compatible during the transition.

Validation:

```bash
npm run build
npm run test:unit
```

### Phase 4: Move frontend QE helpers behind re-export shims

PR size: one helper group at a time.

- Move QE-specific helpers from `src/lib` to `src/engines/qe/lib`.
- Start with low-risk pure helpers such as `qeProgress`, `qeBravaisInference`, `hubbard`, and `wannierQuality`.
- Leave re-export files in old locations until all imports are migrated.

Validation:

```bash
npm run build
npm run test:unit
```

### Phase 5: Move backend pseudopotential code into the QE engine

PR size: backend module extraction only.

- Move UPF, SSSP, cutoff, and pseudopotential repair helpers out of `src-tauri/src/lib.rs`.
- Keep Tauri command names stable.
- Keep local and remote pseudopotential behavior unchanged.

Validation:

```bash
cargo check --manifest-path src-tauri/Cargo.toml
```

Recommended additional check:

```bash
cargo test --manifest-path src-tauri/Cargo.toml pseudopotential
```

### Phase 6: Extract QE task orchestration

PR size: one workflow family per PR.

- Move backend task bodies from `src-tauri/src/lib.rs` into QE task modules.
- Start with a workflow that has clear inputs and outputs.
- Keep command names such as `start_scf_calculation` as compatibility wrappers.

Suggested order:

1. SCF and optimization.
2. Bands and DOS.
3. Fermi surface.
4. Hubbard LRT.
5. Phonons.
6. Wannier and transport.
7. EPW.

Validation for each PR:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

### Phase 7: Split HPC profile shell from engine toolchains

PR size: type compatibility plus adapters.

- Keep cluster identity shared.
- Move QE executable and pseudopotential paths into a QE toolchain sub-structure.
- Keep legacy config loading compatible with existing profiles.
- Do not add Wien2k fields to QE structures.

Validation:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

### Phase 8: Add engine job bundle interface

PR size: interface plus one QE adapter.

- Define a shared HPC job bundle shape:
  - engine id
  - task kind
  - label
  - local bundle files
  - remote commands
  - artifact sync policy
  - result parser
- Adapt one existing QE workflow without changing behavior.
- Keep old direct code path available until the adapter is proven.

Validation:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

## Files That Should Not Be Refactored in Parallel

- `src-tauri/src/lib.rs` and `src-tauri/src/qe/*`.
- `src-tauri/src/projects.rs` and `src/lib/TaskContext.tsx`.
- `src/lib/types.ts`, `src/lib/hpcConfig.ts`, and `src-tauri/src/hpc/profile.rs`.
- `src/App.tsx` and `src/ViewerApp.tsx`.
- Workflow wizard files and their backend task payloads.

## Completion Criteria For Each Migration PR

- Existing QE workflows still run through the same UI.
- Existing saved projects still load.
- Existing command names are either preserved or wrapped.
- No universal fake SCF config is introduced.
- QE pseudopotential behavior remains inside QE-owned code.
- Shared viewers consume normalized datasets only through explicit adapters.
- Validation commands are run and failures are documented with exact logs.

