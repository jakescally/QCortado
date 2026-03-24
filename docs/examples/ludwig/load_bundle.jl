using JSON3
using Ludwig
using Interpolations
using LinearAlgebra
using DelimitedFiles

function periodic_extend(grid::Matrix{Float64})
    nx, ny = size(grid)
    extended = Array{Float64}(undef, nx + 1, ny + 1)
    extended[1:nx, 1:ny] = grid
    extended[nx + 1, 1:ny] = grid[1, :]
    extended[1:nx, ny + 1] = grid[:, 1]
    extended[nx + 1, ny + 1] = grid[1, 1]
    return extended
end

function build_ludwig_loader(bundle_dir::AbstractString)
    metadata_path = joinpath(bundle_dir, "metadata.json")
    bands_path = joinpath(bundle_dir, "bands.csv")

    metadata = JSON3.read(read(metadata_path, String))
    rows = readdlm(bands_path, ',', Float64; skipstart = 1)

    nkx = Int(metadata.grid_shape[1])
    nky = Int(metadata.grid_shape[2])
    nbands = Int(metadata.band_count)

    a1 = collect(metadata.in_plane_lattice_angstrom[1])
    a2 = collect(metadata.in_plane_lattice_angstrom[2])
    lattice = Ludwig.Lattices.Lattice(hcat(a1, a2))

    energy_grids = [zeros(Float64, nkx, nky) for _ in 1:nbands]
    for row_index in 1:size(rows, 1)
        ix = Int(rows[row_index, 1]) + 1
        iy = Int(rows[row_index, 2]) + 1
        for band_index in 1:nbands
            energy_grids[band_index][ix, iy] = rows[row_index, 4 + band_index]
        end
    end

    x_axis = collect(range(0.0, 1.0; length = nkx + 1))
    y_axis = collect(range(0.0, 1.0; length = nky + 1))

    band_interpolants = map(energy_grids) do grid
        periodic_grid = periodic_extend(grid)
        interpolate(periodic_grid, BSpline(Cubic(Line(OnGrid())))) |>
            itp -> scale(itp, x_axis, y_axis)
    end

    rlv = Ludwig.Lattices.reciprocal_lattice_vectors(lattice)
    invrlv = inv(rlv)

    band_functions = map(band_interpolants) do itp
        k -> begin
            frac = mod.(invrlv * collect(k), 1.0)
            itp(frac[1], frac[2])
        end
    end

    return (
        metadata = metadata,
        lattice = lattice,
        band_interpolants = band_interpolants,
        band_functions = band_functions,
    )
end

function main()
    if isempty(ARGS)
        error("Usage: julia load_bundle.jl /path/to/exported_bundle")
    end

    bundle_dir = ARGS[1]
    loaded = build_ludwig_loader(bundle_dir)

    temperature_ev = 0.01
    n_levels = 24
    n_cuts = 96

    mesh = Ludwig.bz_mesh(loaded.lattice, loaded.band_functions, temperature_ev, n_levels, n_cuts)

    println("Loaded QCortado Ludwig bundle from: ", bundle_dir)
    println("Bands: ", loaded.metadata.band_count)
    println("Grid: ", loaded.metadata.grid_shape)
    println("Mode: ", loaded.metadata.mode)
    println("Chemical potential (eV): ", loaded.metadata.chemical_potential_ev)
    println("Generated Ludwig mesh with ", length(mesh.patches), " patches.")
    println()
    println("Next step:")
    println("  supply your own Weff_squared(...) or impurity model")
    println("  and use loaded.band_interpolants inside Ludwig.electron_electron(...)")
end

abspath(PROGRAM_FILE) == @__FILE__ && main()
