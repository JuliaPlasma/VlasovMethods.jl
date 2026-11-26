# The linear Landau damping rate of the spline Galerkin model in Float64 and Float32; the runs
# and the root are in `test/helpers/spline_galerkin.jl`.

using Test
using VlasovMethods

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

@testset "linear Landau rate ($T)" for T in (Float64, Float32)
    landau_rate_test(T)
end
