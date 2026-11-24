# The linear Landau damping rates of the augmented and the midpoint BUG at the rank 16; the runs
# and the root are in `test/helpers/lowrank_bug_rates.jl`. The Landau rates run in Float64 only:
# in Float32 the maxima fit of the damped log|E₁| is not monotone in τ and does not track Float64.
# The fixed-rank BUG has no Landau check: its Landau rate at the rank 16 is not monotone in τ
# either.

using Test
using VlasovMethods

include("../helpers/lowrank_bug_rates.jl")

@testset "linear Landau, $(nameof(typeof(method)))" for method in (AugmentedBUG(), MidpointBUG())
    t, y = log_e1(Float64, method, maxwellian, 0.01, 40.0)
    γ = landau_rate_fit(t, y, (0.0, 40.0)).γ
    @test γ < 0
    @test abs(γ - LANDAU_RATE) / abs(LANDAU_RATE) ≤ 0.02
end
