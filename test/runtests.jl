using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Spline Distribution" include("distributions/spline_distribution.jl")
    @safetestset "Grid Distribution" include("distributions/grid_distribution.jl")
    @safetestset "Particle Distribution" include("distributions/particle_distribution.jl")
    @safetestset "Projections" include("projections/projections.jl")
    @safetestset "Reduced Tensors" include("gridbased/reduced_tensors.jl")
end
if "broken" in GROUPS
    @safetestset "Electric Fields" include("electric_field.jl")   # issue #56
end
