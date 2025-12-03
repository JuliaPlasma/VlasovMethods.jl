using SafeTestsets

@safetestset SplineTests = "$(rpad("Spline",80))" begin
    include("spline_tests.jl")
end
@safetestset SplineBasisTests = "$(rpad("Spline Basis Evaluation",80))" begin
    include("spline_basis_tests.jl")
end
# @safetestset ProjectionTests = "$(rpad("Projections",80))" begin
#     include("projections_tests.jl")
# end
# @safetestset ElectricFieldTests = "$(rpad("Electric Fields",80))" begin
#     include("electric_field_tests.jl")
# end
