# The reduced bases of the bump-on-tail example, from the snapshots of `bump_on_tail_training.jl`:
# the cotangent-lift basis of the particles by the eigendecomposition of the Gram matrix, the POD
# basis of the accelerations, and the DEIM indices of that field basis. It also writes the reduced
# data of the training sweep, the coordinates `X` and the accelerations `A` in the particle basis,
# which the HNN of the reduced field is trained on.
#
#   julia --project=scripts scripts/bump_on_tail_projections.jl [--runid bump_on_tail] \
#       [--particle-tol 1e-8] [--field-tol 1e-4] [--runs-dir dir] [--results-dir dir]
#
# Each basis keeps the leading vectors that carry all but `tol` of the energy of its spectrum. The
# defaults are the second pair of tolerances the old script names, which give `kₚ = 120` and
# `kₑ = 345` on the default training run. The pair it ran, `1e-4` and `1e-2`, gives `kₚ = 10` and
# `kₑ = 55`; the reduced field energy of the smallest test `κ` then has three local maxima, and the
# growth-rate fit of `bump_on_tail_testing.jl` takes the second to the fifth. The run reads
# `<runid>_training.h5` from
# `--runs-dir` and writes `<runid>_projections.h5` and `<runid>_vectorfield.h5` there, and the two
# spectra to `--results-dir`.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using ReducedBasisMethods: CotangentLiftEVD, PODEVD, get_DEIM_indices, h5save, spectrumplot
using VlasovMethods

const OPT = parse_arguments(ARGS; runid = BUMP_ON_TAIL_DEFAULTS.runid, particle_tol = 1e-8,
    field_tol = 1e-4)

tpath = runs_path(OPT.runid * "_training.h5")
SS, seed = h5open(tpath, "r") do file
    Snapshots(file, "snapshots"), read(HDF5.attributes(file)["seed"])
end

# the bases
pb = particle_basis(CotangentLiftEVD(), SS; tolerance = OPT.particle_tol)
fb = field_basis(PODEVD(), SS; tolerance = OPT.field_tol)
deim = get_DEIM_indices(fb.Φ)
@printf("bump_on_tail_projections: kₚ = %d, kₑ = %d\n", pb.k, fb.k)

commit = repository_commit()
command = command_line("bump_on_tail_projections.jl", ARGS)

function provenance!(file)
    HDF5.attributes(file)["commit"] = commit
    HDF5.attributes(file)["command"] = command
    HDF5.attributes(file)["seed"] = seed
    return file
end

ppath = runs_path(OPT.runid * "_projections.h5")
h5open(ppath, "w") do file
    h5save(file, pb; path = "particle_basis")
    h5save(file, fb; path = "field_basis")
    file["deim_indices"] = deim
    provenance!(file)
end
@printf("  wrote %s\n", ppath)

# the reduced data of the training sweep, `(k, nₛ, nparam)` each: the particle basis `Φ` reduces
# the positions and the accelerations alike, as the reduced march's kick does
Φ = pb.Φ
X = zeros(pb.k, size(SS.X, 3), size(SS.X, 4))
A = zeros(pb.k, size(SS.A, 3), size(SS.A, 4))
for p in axes(SS.X, 4), ts in axes(SS.X, 3)

    mul!(view(X, :, ts, p), Φ', view(SS.X, 1, :, ts, p))
    mul!(view(A, :, ts, p), Φ', view(SS.A, 1, :, ts, p))
end

vpath = runs_path(OPT.runid * "_vectorfield.h5")
h5open(vpath, "w") do file
    file["X"] = X
    file["A"] = A
    provenance!(file)
end
@printf("  wrote %s\n", vpath)

# item 11: the two spectra
for (basis, name, what) in ((pb, "particle", "the particles, [X V]"),
    (fb, "field", "the accelerations A"))
    fig = spectrumplot(basis.Λ; title = "the spectrum of " * what, ylabel = "Λᵢ")
    figure = results_path(OPT.runid * "_spectrum_" * name * ".png")
    save(figure, fig)
    @printf("  wrote %s\n", figure)
end
