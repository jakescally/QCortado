# QCortado -> Ludwig Example

This directory contains a minimal Julia-side loader for QCortado's Ludwig export bundles.

The exported bundle includes:

- `metadata.json`
- `bands.csv`
- raw Wannier provenance files such as `.win`, `_hr.dat`, `_wsvec.dat`, and `.wout`

The loader:

1. reads the exported metadata and band grid,
2. constructs a 2D `Ludwig.Lattices.Lattice`,
3. builds periodic `ScaledInterpolation` band interpolants,
4. exposes 2D band functions suitable for `Ludwig.bz_mesh`.

QCortado does not choose a scattering model. After loading a bundle, you still need to provide either:

- your own `Weff_squared(p1, p2, p3, p4; kwargs...)` for electron-electron transport, or
- an impurity model for `electron_impurity!`.

## Julia packages

The example script assumes the following packages are available in your Julia environment:

- `Ludwig`
- `Interpolations`
- `JSON3`

Install any missing packages with:

```julia
using Pkg
Pkg.add(["Ludwig", "Interpolations", "JSON3"])
```

## Usage

```julia
julia --project=. docs/examples/ludwig/load_bundle.jl /path/to/exported_bundle
```

The script prints a short summary and builds a Ludwig Fermi-surface mesh from the exported bands. From there, insert your own scattering model where indicated in the script.
