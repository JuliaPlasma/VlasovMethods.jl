# The two-stream growth rate of the spline Galerkin model in Float64 and Float32; the runs and
# the root are in `test/helpers/spline_galerkin.jl`.

using Test
using VlasovMethods

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

@testset "two-stream rate ($T)" for T in (Float64, Float32)
    two_stream_rate_test(T)
end
