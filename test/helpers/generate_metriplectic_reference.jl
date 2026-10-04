# Regenerate `test/data/lenard_bernstein_metriplectic_reference.txt`, the pre-P44 result of one
# Picard step of the metriplectic Lenard–Bernstein model.
#
# The reference is the output of the `NonlinearSolve` NLsolve/Anderson solve that P44 removes, so
# this runs only against the pre-P44 source. Run it on `origin/main`, in a scratch environment
# that has `VlasovMethods` (developed at the checkout), `NLsolve` and `LineSearches`, and with
# `NLsolve` and `LineSearches` loaded so that the `NonlinearSolve` extension is active:
#
#     julia --startup-file=no --project=<scratch-env> test/helpers/generate_metriplectic_reference.jl
#
# It writes one velocity per line, `%.17e`, and fails if the source no longer carries the
# `NonlinearSolve` path, so it can never overwrite the reference with the new solver's output.
using VlasovMethods
using NLsolve
using LineSearches
using Printf

isdefined(VlasovMethods, :NonlinearSolve) || error(
    "this generator runs only against the pre-P44 source, whose solve uses `NonlinearSolve`")

const N = 64

pdist = ParticleDistribution(1, 1, N)
v0 = collect(range(-2.0, 2.0; length = N))
pdist.particles.v[1, :] .= v0
pdist.particles.w[1, :] .= 1 / N

axis = BSplineBasis(UniformMesh(16, -2.0 .. 2.0), 3, Free())
sdist = SplineDistribution(1, 1, axis,
    VlasovMethods.SimpleSplines.SplineQuadrature(axis), zeros(Float64, nbasis(axis)))
mlb = MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist))

vref = Picard_iterate_over_particles(
    zeros(N), v0, v0, zeros(N, 2), 1, 0.0, 1e-3, 3, 0.5, 3e-16 * sqrt(N), 1e-50, mlb).u

open(joinpath(@__DIR__, "..", "data", "lenard_bernstein_metriplectic_reference.txt"), "w") do io
    for x in vref
        @printf(io, "%.17e\n", x)
    end
end
