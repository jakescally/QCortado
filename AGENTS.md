# AGENTS.md

## Project summary

QCortado is currently a Tauri + React + Rust desktop app for Quantum ESPRESSO workflows. We are migrating the architecture toward an engine-based Cortado platform while preserving existing QE behavior.

The short-term goal is not to implement Wien2k. The short-term goal is to isolate QE-specific logic behind explicit engine boundaries so a Wien2k backend can later be added cleanly.

## Non-negotiable constraints

- Do not remove existing QE functionality unless explicitly instructed.
- Do not change calculation behavior during engine-boundary refactors.
- Do not mutate QE-specific files into vague generic code unless the task explicitly introduces an engine abstraction.
- Do not create a fake universal SCF config containing QE and Wien2k fields together.
- Engine-specific inputs stay engine-specific.
- Shared viewers consume normalized, engine-neutral result datasets.
- Prefer small PR-sized diffs over sweeping rewrites.
- If a task is exploratory, do not edit files.
- If a task is a refactor, preserve behavior and update imports/tests.
- If unsure whether a concept is shared or engine-specific, leave it engine-specific and document the question.

## Important architecture direction

Shared/platform concepts:
- project browser
- task queue
- live output
- HPC profile shell
- SSH/Slurm execution infrastructure
- artifact sync shell
- CIF parsing
- crystal/unit-cell viewing
- Brillouin-zone and k-path UI
- band/DOS/transport viewers
- normalized result datasets

QE-specific concepts:
- pw.x, bands.x, projwfc.x, ph.x, q2r.x, matdyn.x, hp.x, epw.x
- QE namelists/cards
- pseudopotentials and UPF parsing
- SSSP metadata
- ecutwfc/ecutrho
- QE smearing names
- QE Hubbard card
- QE .save directories
- prefix/outdir/pseudo_dir
- QE output parsing

Future Wien2k-specific concepts:
- case directory lifecycle
- case.struct
- RMT/RKmax/Gmax/lmax
- init_lapw flow
- x nn, x sgroup, x symmetry, lstart, kgen, dstart
- run_lapw/runsp_lapw
- case.scf parsing
- lapw1/lapw2/spaghetti workflows

## Validation commands

After TypeScript/frontend changes:
- npm run build
- npm run test:unit

After Rust/backend changes:
- cargo check --manifest-path src-tauri/Cargo.toml

After cross-boundary changes:
- npm run build
- npm run test:unit
- cargo check --manifest-path src-tauri/Cargo.toml

## Done means

- Existing behavior is preserved unless the task explicitly says otherwise.
- The app builds.
- Relevant tests pass or failures are explained with exact logs.
- New types/modules are documented enough for the next migration phase.
- The final response lists files changed, architectural intent, validation commands run, and remaining risks.
