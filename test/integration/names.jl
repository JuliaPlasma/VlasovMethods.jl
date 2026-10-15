using Test
using VlasovMethods

# A recorded list of the exported names of `VlasovMethods`, so that no export can leave the
# module without a `CHANGELOG.md` line that names it. A removed name is written as two joined
# fragments, so that its literal form occurs nowhere under `test/`.
const BASE_EXPORTS = Symbol[
    :..,
    :AnalyticHamiltonian,
    :BSplineBasis,
    :Bump,
    :BumpOnTail,
    :CollisionEntropy,
    :ConservativeLenardBernstein,
    :Constraint,
    :DEIMElectricField,
    Symbol("DiffEq", "Integrator"),
    :Dirichlet,
    :DoubleMaxwellian,
    :ElectricField,
    :Free,
    :FullyReducedTensor,
    :GeneralMesh,
    :GeometricIntegrator,
    :GradedMesh,
    :GridDistribution,
    :IntegratorParameters,
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
    :PoissonField,
    :PotentialReducedTensor,
    :RandomMesh,
    :RecombinedBSplineBasis,
    :ReducedElectricField,
    :ReducedIntegratorCache,
    :RescaledConservativeLenardBernstein,
    :Robin,
    :ShiftedNormalV,
    :ShiftedUniformDistribution,
    :Snapshots,
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
    :bracket_tensor,
    :breakpoints,
    :check_conservation_basis,
    :degree,
    :derivative,
    :domain,
    :electric_field,
    :evaluate,
    :evaluate_all,
    :evaluate_all!,
    :field_basis,
    :initialize!,
    :l2_projection,
    :l2_projection!,
    :lift_solution,
    :local_width,
    :magnetic_field,
    :mass_matrix,
    :mass_operator,
    :mass_solve!,
    :meshwidth,
    :nbasis,
    :ncells,
    :order,
    :particle_basis,
    :plot_distribution,
    :plot_particles,
    :polynomial_reproduction,
    :project_Maxwellian,
    :project_function,
    :projection,
    :quadrature_nodes,
    :quadrature_weights,
    :reduced_integrate_vp,
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

    # The field interface exports its four types, but not the step functions: callers reach
    # `update!` and `efield!` through qualified names or explicit imports.
    @test isdisjoint((:update!, :efield!), names(VlasovMethods))

    # `ReducedTensor` belongs to `ReducedBasisMethods`, which this module imports it from by name,
    # so it is not an export here and a caller has to reach it at its own package
    @test !(:ReducedTensor in names(VlasovMethods))

    changelog = read(joinpath(pkgdir(VlasovMethods), "CHANGELOG.md"), String)
    for name in REMOVED_EXPORTS
        @test occursin("- **`" * string(name) * "`", changelog)
    end
end
