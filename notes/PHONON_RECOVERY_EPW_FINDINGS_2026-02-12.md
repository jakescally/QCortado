# Phonon Recovery + EPW Findings (2026-02-12)

## Context
- User ran a long `ph.x` job (~15h wall), output showed `JOB DONE.`, then frontend reported:
  - `Error: ph.x failed with exit code: Some(1)`
- Priority was data salvage first, then understanding root cause and EPW implications.

## Data Salvage Outcome
- Original scratch run preserved at:
  - `/tmp/qcortado_phonon`
- Full backup created at:
  - `/Users/jakescally/Library/Application Support/com.qcortado.app/recovery_backups/phonon_recovery_20260211_150037`
- Recovered/saved calc in project:
  - Project ID: `1770746392230_42172df9`
  - CIF ID: `1770746397222_767ca989`
  - Recovered calc ID: `1770840797212_1a2359c5`
  - Path:
    - `/Users/jakescally/Library/Application Support/com.qcortado.app/projects/1770746392230_42172df9/calculations/1770840797212_1a2359c5`
- Verified saved dispersion payload:
  - `n_modes = 24`
  - `n_qpoints = 81`
  - frequency range `[-1354.4886, 1513.1164] cm^-1`
  - `converged = true`

## What Actually Happened in QE
- Evidence that phonon run completed despite nonzero exit:
  - `status_run.xml` showed `STOPPED_IN = dynmatrix`, `RECOVER_CODE = 30`
  - `dynmat4` contains final mode block and frequencies
  - manual `q2r.x` and `matdyn.x` reruns succeeded
- Interpretation:
  - `ph.x` exit code handling was stricter than output/convergence markers.

## Code Changes Already Made

### Backend: tolerate `JOB DONE` with nonzero `ph.x` exit
- File:
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src-tauri/src/lib.rs`
- Behavior:
  - If `ph.x` exits nonzero but parsed output is converged (`JOB DONE`), continue pipeline with warning.
  - Still error if not converged.

### Backend: compact phonon saves + explicit recovery command
- File:
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src-tauri/src/projects.rs`
- Added:
  - `recover_phonon_calculation` command
  - compact artifact copier for phonons (`copy_compact_phonon_artifacts`)
- Current compact phonon save keeps:
  - `dynmat*`, `force_constants`, `phonon_freq`, `phonon_freq.gp`, `phonon_dos` (if present)
  - `ph.in`, `q2r.in`, `matdyn_*.in`
  - small `_ph0/...phsave` and `data-file-schema.xml` metadata
- Current compact phonon save drops:
  - large wavefunction shards and most heavy scratch payloads.

### Frontend: recovery button on dashboard
- File:
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src/components/ProjectDashboard.tsx`
- Added button:
  - `Recover Phonon`
- Updated flow:
  - auto-try `/tmp/qcortado_phonon`
  - if unavailable/fails, prompt directory picker
  - show visible info/error banner

## Build / Validation Status
- `cargo check` passed for `src-tauri`.
- `npx tsc --noEmit` passed for frontend.
- `npm run build` failed due local optional Rollup package issue:
  - missing `@rollup/rollup-darwin-arm64` in this environment.

## EPW Findings (Important)

### Current recovered run is NOT EPW-ready
- `ph.in` in this run did **not** set `fildvscf` and did not request electron-phonon outputs.
- No `dvscf` files found in this scratch tree.

### How QE decides to write EPW-relevant files
- In `ph.x` input:
  - `fildvscf` enables writing potential variation for later e-ph use.
  - `ldisp=.true.` + `fildvscf` or electron-phonon mode triggers per-q subdirectories (`lqdir` behavior).
- References:
  - `qe-7.5/Doc/INPUT_PH.txt` entries for `fildvscf`, `electron_phonon`, `lqdir`, `dvscf_star`.

### EPW examples indicate a targeted subset is sufficient
- EPW example `ph.in` files use:
  - `fildvscf = 'dvscf'`, `ldisp=.true.`
- EPW example post-processing scripts copy:
  - `dyn*`, `dvscf*`, `phsave`
  - and explicitly delete `*wfc*` in q-folders.
- Implication:
  - do not blindly keep full scratch forever; keep EPW-required subset.

## Scratch Size Observations (BAs case)
- `/tmp/qcortado_phonon` total ~8.1G
- Main contributors:
  - `/tmp/qcortado_phonon/tmp/_ph0` ~5.5G
  - `/tmp/qcortado_phonon/tmp/qcortado_scf.save` ~1.4G
  - largest files are many `qcortado_scf.wfc*` shards (~390-397M each)

## Recommended Next Work (when resumed)
1. Add explicit phonon run mode for EPW:
   - set `fildvscf` in generated `ph.in`
   - optionally set relevant `electron_phonon` options based on target workflow.
2. Add save policy options for phonons:
   - `compact` (current)
   - `epw` (keep `dyn*`, `dvscf*`, required `phsave`/metadata, skip bulky unnecessary wfc shards where safe)
   - optional `full` archive for debugging.
3. Add UI labeling:
   - show whether a phonon run is `plot-only`, `EPW-ready`, or `full archive`.

## Notes
- Working tree already had unrelated edits in:
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src/App.css`
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src/App.tsx`
- Additional edits from this task are in:
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src-tauri/src/lib.rs`
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src-tauri/src/projects.rs`
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src/components/ProjectDashboard.tsx`
  - `/Users/jakescally/Applications/QuantumEspresso/QCortado/src/App.css`
