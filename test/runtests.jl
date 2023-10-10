using VlasovMethods
using Test

@testset "VlasovMethods.jl" begin
    include("electric_field_tests.jl")
    include("projections_tests.jl")
    # include("2d_spline_tests.jl")
end
