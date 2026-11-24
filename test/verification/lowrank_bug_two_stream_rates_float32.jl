# The two-stream growth rates of the three BUG integrators at the rank 16 in Float32; the runs and
# the root are in `test/helpers/lowrank_bug_rates.jl`. The growing mode stays far above round-off,
# so its rates run in both precisions; `lowrank_bug_two_stream_rates.jl` runs them in Float64.

using Test
using VlasovMethods

include("../helpers/lowrank_bug_rates.jl")

@testset "two-stream instability, $(nameof(typeof(method))), Float32" for method in (
    BUG(), AugmentedBUG(), MidpointBUG())
    two_stream_rate_test(Float32, method)
end
