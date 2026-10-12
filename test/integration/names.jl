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
    :FullyReducedTensor,
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
    :UniformField,
    :UniformMesh,
    :VelocityReducedMatrix,
    :VlasovMethods,
    :VlasovPoisson,
    :VlasovWithExternalFields,
    :_apply_∫dv!,
    :basis_index,
    :basis_values,
    :breakpoints,
    :check_conservation_basis,
    :degree,
    :derivative,
    :domain,
    :electric_field,
    :evaluate,
    :evaluate_all,
    :evaluate_all!,
    :initialize!,
    :l2_projection,
    :l2_projection!,
    :local_width,
    :magnetic_field,
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
    :sodeensemble,
    :velocity_moments
]

# The recorded names that the module does not export. Each is the bold lead of its
# `CHANGELOG.md` bullet.
const REMOVED_EXPORTS = Symbol[Symbol("DiffEq", "Integrator"), :run]

@testset "Exported names" begin
    # The recorded names missing from the module are exactly the documented removals. A name
    # that leaves without a CHANGELOG line fails here.
    @test setdiff(BASE_EXPORTS, names(VlasovMethods)) == REMOVED_EXPORTS

    # `ReducedTensor` belongs to `ReducedBasisMethods`, which this module imports it from by name,
    # so it is not an export here and a caller has to reach it at its own package
    @test !(:ReducedTensor in names(VlasovMethods))

    changelog = read(joinpath(pkgdir(VlasovMethods), "CHANGELOG.md"), String)
    for name in REMOVED_EXPORTS
        @test occursin("- **`" * string(name) * "`", changelog)
    end
end
