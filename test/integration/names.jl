using Test
using VlasovMethods

# The exported names of `VlasovMethods` on `origin/main` at the base of this part, `2203010`,
# recorded here so that a name this part deletes cannot leave the module without a CHANGELOG line.
# The part removes the `DifferentialEquations` paths, so `DiffEqIntegrator` is the one name the
# branch loses; the loop below asserts that every lost name is named in `CHANGELOG.md`.
const BASE_EXPORTS = Symbol[
    :..,
    :BSplineBasis,
    :Bump,
    :BumpOnTail,
    :CollisionEntropy,
    :ConservativeLenardBernstein,
    :Constraint,
    :DiffEqIntegrator,
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

@testset "Exported names" begin
    @test :DiffEqIntegrator in BASE_EXPORTS
    @test :DiffEqIntegrator ∉ names(VlasovMethods)

    changelog = read(joinpath(pkgdir(VlasovMethods), "CHANGELOG.md"), String)
    for name in setdiff(BASE_EXPORTS, names(VlasovMethods))
        @test occursin(string(name), changelog)
    end
end
