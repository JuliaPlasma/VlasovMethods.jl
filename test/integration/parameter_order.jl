using ParticleMethods: ParticleList
using PoissonSolvers
using Random
using SimpleSplines: SplineQuadrature, TensorProductQuadrature
using Test
using VlasovMethods

# The value a type binds to the parameter called `name`, looked up by name rather than by
# position, so that the test says which parameter it means.
function typeparameter(x, name::Symbol)
    T = typeof(x)
    names = [p.name for p in Base.unwrap_unionall(T.name.wrapper).parameters]
    T.parameters[findfirst(==(name), names)]
end

# `ParticleDistribution(xdim, vdim, npart)` always fills `Float64` zeros, so the `Float32` case
# builds its `ParticleList` from a `Float32` matrix, with the same variables.
function particle_distribution(::Type{T}, xdim, vdim, npart) where {T}
    z = zeros(T, xdim + vdim + 1, npart)
    vars = (
        x = 1:xdim,
        v = (xdim + 1):(xdim + vdim),
        z = 1:(xdim + vdim),
        w = (xdim + vdim + 1):(xdim + vdim + 1)
    )
    ParticleDistribution(xdim, vdim, ParticleList(z; variables = vars))
end

# The grid constructor of `SplineDistribution` always fills `Float64` coefficients, so the
# distribution is built from a basis on a `T` mesh and `T` coefficients instead.
function spline_distribution(::Type{T}, vdim) where {T}
    axis = BSplineBasis(UniformMesh(7, -T(5) .. T(5)), 2, Free())
    if vdim == 1
        SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    else
        B = TensorProductBasis(ntuple(_ -> axis, vdim))
        SplineDistribution(1, vdim, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    end
end

@testset "Parameter order of the models, $T" for T in (Float64, Float32)
    npart = 100

    Random.seed!(1234)
    pdist1 = initialize!(particle_distribution(T, 1, 1, npart), NormalDistribution())
    sdist1 = spline_distribution(T, 1)
    ent1 = CollisionEntropy(sdist1)

    @test eltype(pdist1) == T
    @test eltype(sdist1) == T

    @testset "CollisionEntropy" begin
        @test typeparameter(ent1, :XD) === 1
        @test typeparameter(ent1, :VD) === 1
    end

    @testset "VlasovPoisson" begin
        potential = Potential(PeriodicBasisSpline((T(0), T(1)), 3, 8))
        model = VlasovPoisson(pdist1, potential)
        @test typeparameter(model, :XD) === 1
        @test typeparameter(model, :VD) === 1
    end

    @testset "$(M)" for M in (LenardBernstein, ConservativeLenardBernstein,
        RescaledConservativeLenardBernstein, MetriplecticLenardBernstein)
        model = M(pdist1, ent1)
        @test typeparameter(model, :XD) === 1
        @test typeparameter(model, :VD) === 1
    end

    @testset "Landau" begin
        pdist2 = initialize!(particle_distribution(T, 1, 2, npart), NormalDistribution())
        sdist2 = spline_distribution(T, 2)
        model = Landau(pdist2, CollisionEntropy(sdist2))
        @test typeparameter(model, :XD) === 1
        @test typeparameter(model, :VD) === 2
    end
end
