using Test
using VlasovMethods

# The export baseline covers the names `v0.2.1` exports and every name exported since that tag.
# No export can leave the module without a `CHANGELOG.md` line that names it, and
# `REMOVED_EXPORTS` carries those names. A removed name that occurs
# nowhere under `test/` is written as two joined fragments, so that its literal form does not
# appear here either.
const BASE_EXPORTS = Symbol[
    :..,
    :AnalyticHamiltonian,
    :BSplineBasis,
    :Bump,
    :BumpOnTail,
    :CLB_rhs,
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
    :LB_rhs,
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
    :Picard_iterate_over_particles!,
    :PoissonField,
    :PotentialReducedTensor,
    :RandomMesh,
    :RecombinedBSplineBasis,
    :ReducedElectricField,
    :ReducedIntegratorCache,
    :RescaledConservativeLenardBernstein,
    :Robin,
    :ScaledField,
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
    :best_and_worst,
    :bracket_tensor,
    :breakpoints,
    :check_conservation_basis,
    :degree,
    :derivative,
    :domain,
    :electric_field,
    :energy_drift,
    :evaluate,
    :evaluate_all,
    :evaluate_all!,
    :field_basis,
    :field_energy_error,
    :growth_rate,
    :initialize!,
    :integrate_vp,
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
    :plot_series,
    :polynomial_reproduction,
    :potential_error,
    :project_Maxwellian,
    :project_function,
    :projection,
    :quadrature_nodes,
    :quadrature_weights,
    :reduced_integrate_vp,
    :run,
    :run!,
    :saved_indices,
    :sodeensemble,
    :trajectory_errors,
    :velocity_moments
]

# The recorded names that the module does not export. Each is the bold lead of its
# `CHANGELOG.md` bullet.
const REMOVED_EXPORTS = Symbol[:CLB_rhs, Symbol("DiffEq", "Integrator"), :LB_rhs, :run]

# The names the unreleased section records. Each is the bold lead of its `[Unreleased]`
# bullet, so a name that ships without a line in the record fails here.
const NEW_EXPORTS = Symbol[
    :ScaledField,
    :best_and_worst,
    :energy_drift,
    :field_energy_error,
    :growth_rate,
    :integrate_vp,
    :plot_series,
    :potential_error,
    :saved_indices,
    :train_hnn,
    :trajectory_errors
]

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

@testset "The unreleased section" begin
    changelog = read(joinpath(pkgdir(VlasovMethods), "CHANGELOG.md"), String)

    # The section stands above every release heading, so a reader meets it first.
    @test occursin("## [Unreleased] — targeting 0.4.0", changelog)
    @test first(eachmatch(r"^## \[([^\]]+)\]"m, changelog)).captures[1] == "Unreleased"
    for name in NEW_EXPORTS
        @test occursin("- **`" * string(name) * "`", changelog)
        @test name in names(VlasovMethods)
    end
end

@testset "The release notes' dependency claims" begin
    changelog = read(joinpath(pkgdir(VlasovMethods), "CHANGELOG.md"), String)
    flat = replace(changelog, r"\s+" => " ")  # a sentence is matched across its line wrap

    # The compat bullet explains omitted dependency rows without denying re-exports.
    # `Spline` and `derivative` are re-exported, so a claim that nothing is re-exported is false.
    @test !occursin(r"nothing (here )?is re-exported", flat)
    @test occursin("The dependencies that left or joined `[deps]` in the same interval have no row.", flat)

    # An *Open Issues* heading exists in neither file, so neither file may cite one.
    known_issues = read(joinpath(pkgdir(VlasovMethods), "KNOWN_ISSUES.md"), String)
    @test !occursin("Open Issues", changelog) && !occursin("Open Issues", known_issues)
end
