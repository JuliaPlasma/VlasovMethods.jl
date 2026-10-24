# Regenerate `test/data/lenard_bernstein_metriplectic_trajectory.txt`, ten Picard steps of the
# metriplectic Lenard–Bernstein model from a seeded state, under the solve of commit `472dd42`.
#
# The seeded state is the lattice of cell centres of [-2, 2], one per particle, each moved by a
# draw uniform in ±1/64 from `StableRNG(20261010)`, with weights 1/N. The march takes
# `ti = 1:10` at `t = (ti - 1)·Δt`, as `scripts/lenard_bernstein_metriplectic.jl` does, from
# `vn = vn_minus_one = v0` and a zero history; after each step `vn_minus_one` takes `vn`, and `vn`
# the result. It calls only `Picard_iterate_over_particles`, which `472dd42` has, so it runs
# against an exported tree of that commit in a scratch environment with `StableRNGs`:
#
#     julia --startup-file=no --project=<scratch-env> test/helpers/generate_metriplectic_trajectory.jl
#
# It writes the ten states in order, one velocity per line, `%.17e`: 640 lines.
using VlasovMethods
using StableRNGs
using Printf

const N = 64
const ΔT = 1e-3
const SEED = 20261010

pdist = ParticleDistribution(1, 1, N)
v0 = [-2 + (i - 0.5) * 4 / N for i in 1:N] .+ (rand(StableRNG(SEED), N) .- 0.5) ./ 32
pdist.particles.v[1, :] .= v0
pdist.particles.w[1, :] .= 1 / N

axis = BSplineBasis(UniformMesh(16, -2.0 .. 2.0), 3, Free())
sdist = SplineDistribution(1, 1, axis,
    VlasovMethods.SimpleSplines.SplineQuadrature(axis), zeros(Float64, nbasis(axis)))
mlb = MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist))

vn, vn_minus_one, dv_history = copy(v0), copy(v0), zeros(N, 2)
open(joinpath(@__DIR__, "..", "data", "lenard_bernstein_metriplectic_trajectory.txt"), "w") do io
    for ti in 1:10
        v = Picard_iterate_over_particles(zeros(N), vn, vn_minus_one, dv_history, ti,
            (ti - 1) * ΔT, ΔT, 3, 0.5, 3e-16 * sqrt(N), 1e-50, mlb)
        for x in v
            @printf(io, "%.17e\n", x)
        end
        global vn_minus_one, vn = vn, v
    end
end
