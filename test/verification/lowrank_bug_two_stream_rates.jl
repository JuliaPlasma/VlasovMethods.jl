# The two-stream growth rates of the three BUG integrators at the rank 16 in Float64; the runs and
# the root are in `test/helpers/lowrank_bug_rates.jl`, and
# `lowrank_bug_two_stream_rates_float32.jl` runs the same cases in Float32.

using Test
using VlasovMethods

include("../helpers/lowrank_bug_rates.jl")

@testset "two-stream instability, $(nameof(typeof(method))), Float64" for method in (
    BUG(), AugmentedBUG(), MidpointBUG())
    two_stream_rate_test(Float64, method)
end
