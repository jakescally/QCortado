# Frontend Component Boundaries

QCortado is migrating from a QE-specific application into an engine-based Cortado platform. This directory still keeps the large existing components at the top level to avoid high-risk import churn, but ownership is now documented through the `shared/` and `qe/` boundary folders.

## Boundary Folders

- `shared/` is for platform components and viewer surfaces that should be reusable across engines after result adapters normalize engine outputs.
- `qe/` is for Quantum ESPRESSO workflow components and QE-heavy UI that still depends on QE inputs, pseudopotentials, QE executables, or QE saved-result shapes.

The index files in those folders are organizational exports for new code. Existing imports may continue to reference the current top-level component files until a later PR moves files or rewires imports.

## Migration Rules

- Do not move large components just to satisfy folder structure.
- Do not turn engine-specific inputs into universal input configs.
- Prefer adapting engine outputs into normalized viewer datasets.
- If a component mixes shared shell UI with QE-specific fields, keep it in place and document the split before moving it.
- Preserve current QE labels, defaults, and behavior during organization-only changes.
