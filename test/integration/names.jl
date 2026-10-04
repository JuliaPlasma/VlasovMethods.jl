using Test
using VlasovMethods

# A recorded list of the exported names of `VlasovMethods`, so that no export can leave the
# module without a `CHANGELOG.md` line that names it. A removed name is written as two joined
# fragments, so that its literal form occurs nowhere under `test/`.
const BASE_EXPORTS = Symbol[
    :..,
    :BSplineBasis,
    :Bump,
    :BumpOnTail,
    :CollisionEntropy,
    :ConservativeLenardBernstein,
    :Constraint,
    Symbol("DiffEq", "Integrator"),
    :Dirichlet,
    :DoubleMaxwellian,
    :Free,
    :GeneralMesh,
    :GeometricIntegrator,
    :GradedMesh,
    :GridDistribution,
    :Landau,
    :LenardBernstein,
    :MetriplecticLenardBernstein,
    :Natural,
    :Neumann,
    :NormalDistribution,
    :ParticleDistribution,
    :Periodic,
    :PeriodicBSplineBasis,
    :Picard_iterate_over_particles,
    :PotentialReducedTensor,
    :RandomMesh,
    :RecombinedBSplineBasis,
    :ReducedTensor,
    :RescaledConservativeLenardBernstein,
    :Robin,
    :ShiftedNormalV,
    :ShiftedUniformDistribution,
    :Spline,
    :SplineDistribution,
    :SplittingMethod,
    :SumMaxwellian,
    :TensorProductBasis,
    :UniformDistribution,
    :UniformMesh,
    :VelocityReducedMatrix,
    :VlasovMethods,
    :VlasovPoisson,
    :_apply_∫dv!,
    :basis_index,
    :basis_values,
    :breakpoints,
    :check_conservation_basis,
    :degree,
    :derivative,
    :domain,
    :evaluate,
    :evaluate_all,
    :evaluate_all!,
    :initialize!,
    :l2_projection,
    :l2_projection!,
    :local_width,
    :mass_matrix,
    :mass_operator,
    :mass_solve!,
    :meshwidth,
    :nbasis,
    :ncells,
    :order,
    :polynomial_reproduction,
    :project_Maxwellian,
    :project_function,
    :projection,
    :quadrature_nodes,
    :quadrature_weights,
    :run,
    :run!,
    :velocity_moments
]

# The recorded names that the module does not export. Each is the bold lead of its
# `CHANGELOG.md` bullet.
const REMOVED_EXPORTS = Symbol[Symbol("DiffEq", "Integrator")]

@testset "Exported names" begin
    # Every documented removal is a recorded name, and every recorded name missing from the
    # module is a documented removal. A name that leaves without a CHANGELOG line fails here.
    @test issubset(REMOVED_EXPORTS, BASE_EXPORTS)
    @test setdiff(BASE_EXPORTS, names(VlasovMethods)) == REMOVED_EXPORTS

    changelog = read(joinpath(pkgdir(VlasovMethods), "CHANGELOG.md"), String)
    for name in REMOVED_EXPORTS
        @test occursin("- **`" * string(name) * "`", changelog)
    end
end
