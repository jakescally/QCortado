# QE-Specific Components

This folder marks frontend ownership for Quantum ESPRESSO workflows. The existing large components remain in `src/components` for now; `index.ts` provides an explicit QE boundary for future imports without moving files.

## Current QE Owners

- `SCFWizard.tsx`
- `BandStructureWizard.tsx`
- `ElectronicDOSWizard.tsx`
- `FermiSurfaceWizard.tsx`
- `PhononWizard.tsx`
- `HubbardLrtWizard.tsx`
- `WannierWizard.tsx`
- `EpwWizard.tsx`
- `TransportWizard.tsx`
- `EpwViewer.tsx`

These components are QE-specific or QE-heavy because they work with QE executables, QE namelists/cards, pseudopotentials, QE project artifacts, or QE-derived workflow assumptions.

## Rules

- Keep QE input models engine-specific.
- Do not add Wien2k fields to QE wizard configs.
- Do not move these files until imports are straightforward and the diff remains PR-sized.
- Shared controls extracted from these files should move toward `shared/` only when they no longer depend on QE payloads.
