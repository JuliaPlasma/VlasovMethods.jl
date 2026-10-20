# The reduced bases of the bump-on-tail example: the particle basis and the field basis of the
# training sweep, the reduced data that the HNN of P39 reads, and the spectra of the two bases.
#
#   julia --project=scripts scripts/bump_on_tail_projections.jl \
#       [--runid BoT_Np1000_k_010_050_ns10_T25] [--particle-tol 1e-4] [--field-tol 1e-2] \
#       [--runs-dir dir] [--results-dir dir]
#
# The two files go to `--runs-dir`, the spectra to `--results-dir`. The run id names the training
# run the bases are built from, so that a run of other parameters needs its own id.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using CairoMakie
using HDF5
using LinearAlgebra
using Printf
using ReducedBasisMethods: CotangentLiftEVD, PODEVD, ReducedBasis, spectrumplot
using VlasovMethods
using VlasovMethods: field_basis, h5save, particle_basis

import ReducedBasisMethods

const OPT = parse_arguments(ARGS; runid = DEFAULT_RUNID, particle_tol = 1e-4, field_tol = 1e-2)

const RUNID = OPT.runid
const TRAINING = runs_path("$(RUNID).h5")

const SS = h5open(TRAINING, "r") do file
    Snapshots(file)
end

# the integrator parameters of the training run, which the test script marches with
const IP = h5open(TRAINING, "r") do file
    IntegratorParameters(file, "integrator")
end

# the particle basis of the positions and the field basis of the accelerations
const RB = particle_basis(CotangentLiftEVD(), SS; tolerance = OPT.particle_tol)
const FB = field_basis(PODEVD(), SS; tolerance = OPT.field_tol)

@printf("the particle basis holds %d of %d columns, the field basis %d of %d\n",
    RB.k, size(RB.Φ, 1), FB.k, size(FB.Φ, 1))

# the reduced data of P39: the positions and the accelerations of every saved state, projected onto
# the particle basis
const X = SS.X
const A = SS.A
const XR = zeros(size(RB.Φ, 2), size(X, 3), size(X, 4))
const AR = zeros(size(RB.Φ, 2), size(A, 3), size(A, 4))
for p in axes(A, 4), it in axes(A, 3)

    XR[:, it, p] .= RB.Φ' * vec(X[:, :, it, p])
    AR[:, it, p] .= RB.Φ' * vec(A[:, :, it, p])
end

const PROJECTIONS = runs_path("$(RUNID)_projections.h5")
h5open(PROJECTIONS, "w") do file
    ReducedBasisMethods.h5save(file, RB; path = "particle")
    ReducedBasisMethods.h5save(file, FB; path = "field")
    h5save(file, IP; path = "integrator")
    HDF5.attributes(file)["runid"] = RUNID
    HDF5.attributes(file)["commit"] = repository_commit()
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = 1234
end

const VECTORFIELD = runs_path("$(RUNID)_vectorfield.h5")
h5open(VECTORFIELD, "w") do file
    file["X"] = XR
    file["A"] = AR
    HDF5.attributes(file)["runid"] = RUNID
    HDF5.attributes(file)["commit"] = repository_commit()
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = 1234
end

@printf("wrote %s and %s (%d × %d × %d each)\n", PROJECTIONS, VECTORFIELD,
    size(XR, 1), size(XR, 2), size(XR, 3))

# the spectra of the two bases (item 11)
fig = spectrumplot(RB.Λ; ylabel = "λᵢ")
save(results_path("$(RUNID)_particle_spectrum.png"), fig)

fig = spectrumplot(FB.Λ; ylabel = "λᵢ")
save(results_path("$(RUNID)_field_spectrum.png"), fig)
