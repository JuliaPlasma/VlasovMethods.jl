# The reduced bases of the training sweep, their spectra, and the two files a learned model of the
# reduced field reads: the particle basis, and the reduced positions and accelerations of every
# saved state.
#
#   julia --project=scripts scripts/bump_on_tail_projections.jl [--particle-tol 1e-4] \
#       [--field-tol 1e-2] [--runs-dir dir] [--results-dir dir]
#
# Reads `<runid>.h5` from `--runs-dir`, which `bump_on_tail_training.jl` writes there, and writes
# `<runid>_projections.h5`, `<runid>_vectorfield.h5` and two spectra to `--runs-dir` and
# `--results-dir`. The tolerances are the old script's: the particle basis keeps all but `1e-4` of
# the energy of the stacked positions and velocities, the field basis all but `1e-2` of the
# accelerations'.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using ReducedBasisMethods: CotangentLiftEVD, PODSVD, h5save, spectrumplot
using VlasovMethods

const OPT = parse_arguments(ARGS; particle_tol = 1e-4, field_tol = 1e-2)

const T = Float64

fpath = runs_path("$(RUNID_BOT).h5")
SS = h5open(fpath, "r") do file
    VlasovMethods.h5load(Snapshots, file)
end

bp = particle_basis(CotangentLiftEVD(), SS; tolerance = OPT.particle_tol)
be = field_basis(PODSVD(), SS; tolerance = OPT.field_tol)
@printf("bump_on_tail_projections: kₚ = %d, kₑ = %d\n", bp.k, be.k)

# the particle basis in the file the learned model of the reduced field reads, with the same
# commit, command and seed on every file of the set
ppath = runs_path("$(RUNID_BOT)_projections.h5")
h5open(ppath, "w") do file
    h5save(file, bp)
    h5save(file, be; path = "field")
    HDF5.attributes(file)["commit"] = commit_id()
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = SEED_BOT
end
@printf("  wrote %s\n", ppath)

# the reduced positions and accelerations of every saved state, `(kₚ, nₛ, nparam)` each
X, A = SS.X, SS.A
X̃ = zeros(T, bp.k, size(X, 3), size(X, 4))
Ã = zeros(T, bp.k, size(A, 3), size(A, 4))
for p in axes(X, 4), ts in axes(X, 3)

    X̃[:, ts, p] .= bp.Φ' * vec(X[:, :, ts, p])
    Ã[:, ts, p] .= bp.Φ' * vec(A[:, :, ts, p])
end

vpath = runs_path("$(RUNID_BOT)_vectorfield.h5")
h5open(vpath, "w") do file
    file["X"] = X̃
    file["A"] = Ã
    HDF5.attributes(file)["commit"] = commit_id()
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = SEED_BOT
end
@printf("  wrote %s, X and A of size %s\n", vpath, size(X̃))

for (name, b) in (("particles", bp), ("field", be))
    figure = results_path("bump_on_tail_spectrum_$name.png")
    save(figure, spectrumplot(b.Λ; ylabel = "λᵢ", title = "$name basis, k = $(b.k)"))
    @printf("  wrote %s\n", figure)
end
