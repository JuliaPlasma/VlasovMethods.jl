# The energy and the mass of the spline Galerkin model over 10⁴ steps in Float64; the run is in
# `test/helpers/spline_galerkin.jl`, and `spline_galerkin_energy_float32.jl` runs it in Float32.

using Test
using VlasovMethods

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

@testset "energy over 10⁴ steps (Float64)" begin
    energy_drift_test(Float64)
end
