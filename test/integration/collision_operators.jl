using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis,
                     TensorProductBasis, TensorProductQuadrature
using StableRNGs
using Test
using VlasovMethods

# One vector-field evaluation of each of the five collision operators, plus one Picard step of
# the metriplectic operator (which exercises `f!`), on a StableRNGs-seeded input. The references
# in `test/data/` are written by `test/helpers/generate_collision_reference.jl`. A change that
# moves code without changing arithmetic matches them to `rtol = 64*eps(T)`; the last bits vary
# across Julia versions and platforms.
#
# The velocity input is uniform on the spline domain, so the L² projection stays positive and
# neither operator stops on the positivity check (KNOWN_ISSUES.md, K8).
const SEED = 20261010
const N = 256
const N1 = 24
const ΔT = 1e-3

reference(name) = parse.(Float64, readlines(joinpath(@__DIR__, "..", "data", name)))
rtol(T) = 64 * eps(T)

function setup_1d1v(::Type{T}) where {T}
    rng = StableRNG(SEED)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    pdist.particles.v[1, :] .= (T(4) .* rand(rng, T, N)) .- T(2)
    pdist.particles.w[1, :] .= T(1) / N
    axis = BSplineBasis(UniformMesh(17, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    return pdist, sdist
end

lb_params(model) = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)

function setup_1d2v(::Type{T}) where {T}
    rng = StableRNG(SEED)
    npart = N1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, ParticleList(zeros(T, 4, npart); variables = vars))
    pdist.particles.v .= (T(4) .* rand(rng, T, 2, npart)) .- T(2)
    pdist.particles.w[1, :] .= T(1) / npart
    axis = BSplineBasis(UniformMesh(5, T(-2) .. T(2)), 2, Free())
    B = TensorProductBasis((axis, axis))
    sdist = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    return pdist, sdist
end

@testset "Collision operators" begin
    pdist, sdist = setup_1d1v(Float64)
    v₀ = copy(pdist.particles.v[1, :])

    @testset "Lenard–Bernstein" begin
        model = LenardBernstein(pdist, CollisionEntropy(sdist))
        v = copy(v₀)
        v̇ = similar(v)
        VlasovMethods.LB_rhs!(
            v̇, v, (
                ν = model.ν, idist = model.dist, fdist = model.ent.dist), 0.0)
        @test isapprox(v̇, reference("collision_lb_reference.txt"); rtol = rtol(Float64))
    end

    @testset "Conservative Lenard–Bernstein" begin
        model = ConservativeLenardBernstein(pdist, CollisionEntropy(sdist))
        v = copy(v₀)
        v̇ = similar(v)
        VlasovMethods.CLB_rhs_GI!(v̇, 0.0, v, lb_params(model))
        @test isapprox(v̇, reference("collision_clb_reference.txt"); rtol = rtol(Float64))
    end

    @testset "Rescaled conservative Lenard–Bernstein" begin
        model = RescaledConservativeLenardBernstein(pdist, CollisionEntropy(sdist))
        v = copy(v₀)
        v̇ = similar(v)
        VlasovMethods.RCLB_rhs_GI!(v̇, 0.0, v, lb_params(model))
        @test isapprox(v̇, reference("collision_rclb_reference.txt"); rtol = rtol(Float64))
    end

    @testset "Metriplectic Lenard–Bernstein" begin
        model = MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist))
        v = copy(v₀)
        v̇ = similar(v)
        VlasovMethods.collisional_vectorfield!(v̇, v, nothing, model)
        @test isapprox(v̇, reference("collision_mlb_reference.txt"); rtol = rtol(Float64))

        # One Picard step, which exercises `f!`; `@inferred` pins its return type.
        vstep = @inferred Picard_iterate_over_particles(
            zeros(N), v, v, zeros(N, 2), 1, 0.0, ΔT, 3, 0.5, 3e-16 * sqrt(N), 1e-50, model)
        @test isapprox(vstep, reference("collision_mlb_step_reference.txt"); rtol = rtol(Float64))

        # A value pin cannot tell which `f!` method runs, so the dispatch is asserted:
        # a metriplectic call reaches the `f!` typed on `MetriplecticLenardBernstein`.
        m = which(VlasovMethods.f!,
            (Vector{Float64}, Vector{Float64}, Vector{Float64},
                NamedTuple, Float64, typeof(model)))
        @test Base.unwrap_unionall(m.sig).parameters[end] === MetriplecticLenardBernstein
    end

    @testset "Landau" begin
        pdist2, sdist2 = setup_1d2v(Float64)
        model = Landau(pdist2, CollisionEntropy(sdist2))
        v = copy(pdist2.particles.v)
        v̇ = similar(v)
        VlasovMethods.collisional_vectorfield!(v̇, v, nothing, model)
        @test isapprox(vec(v̇), reference("collision_landau_reference.txt"); rtol = rtol(Float64))
    end
end
