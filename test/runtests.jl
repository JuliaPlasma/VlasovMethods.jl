using SafeTestsets

@safetestset SplineDistributionTests = "$(rpad("Spline Distribution",80))" begin
    include("spline_distribution_tests.jl")
end
@safetestset ParticleDistributionTests = "$(rpad("Particle Distribution",80))" begin
    include("particle_distribution_tests.jl")
end
# @safetestset ProjectionTests = "$(rpad("Projections",80))" begin
#     include("projections_tests.jl")
# end
# @safetestset ElectricFieldTests = "$(rpad("Electric Fields",80))" begin
#     include("electric_field_tests.jl")
# end
