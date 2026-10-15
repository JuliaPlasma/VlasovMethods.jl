using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "JET" include("quality/jet.jl")
    @safetestset "Spline Distribution" include("distributions/spline_distribution.jl")
    @safetestset "Grid Distribution" include("distributions/grid_distribution.jl")
    @safetestset "Particle Distribution" include("distributions/particle_distribution.jl")
    @safetestset "Plotting" include("plotting.jl")
    @safetestset "Projections" include("projections/projections.jl")
    @safetestset "Reduced Tensors" include("gridbased/reduced_tensors.jl")
    @safetestset "Geometric Integrator" include("methods/geometric_integrator.jl")
    @safetestset "Splitting Method" include("methods/splitting.jl")
    @safetestset "Parameter Order" include("integration/parameter_order.jl")
    @safetestset "Landau Solver" include("integration/landau_solver.jl")
    @safetestset "Metriplectic Solve" include("integration/metriplectic_solve.jl")
    @safetestset "Collision Operators" include("integration/collision_operators.jl")
    @safetestset "Vlasov–Poisson" include("integration/vlasov_poisson.jl")
    @safetestset "Vlasov with External Fields" include("models/vlasov_with_external_fields.jl")
    @safetestset "ElectromagneticFields Extension" include("integration/electromagnetic_fields_ext.jl")
    @safetestset "SODE Ensemble" include("integration/sodeensemble.jl")
    @safetestset "Names" include("integration/names.jl")
    @safetestset "Collisions" include("gridbased/collisions.jl")
    @safetestset "Plotting extension" include("integration/makie_extension.jl")
    @safetestset "Grid 1D1V Vlasov–Poisson" include("gridbased/vlasov_poisson_1d1v.jl")
    @safetestset "Snapshots" include("particles/snapshots.jl")
    @safetestset "Electric Field" include("electric_field.jl")
    @safetestset "Reduced Electric Field" include("particles/electric_field.jl")
end
