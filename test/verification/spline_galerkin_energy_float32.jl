# The energy and the mass of the spline Galerkin model over 10⁴ steps in Float32; the run is in
# `test/helpers/spline_galerkin.jl`, and `spline_galerkin_energy.jl` runs it in Float64.

using Test
using VlasovMethods

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

@testset "energy over 10⁴ steps (Float32)" begin
    energy_drift_test(Float32)
end
