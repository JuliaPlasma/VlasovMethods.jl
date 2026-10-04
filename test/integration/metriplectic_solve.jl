using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: SplineQuadrature
using Test
using VlasovMethods

# One Picard step of the metriplectic Lenard–Bernstein model. The reference is the result of the
# `NonlinearSolve` path that P44 removes; `test/helpers/generate_metriplectic_reference.jl` wrote
# it, on `origin/main`, in a scratch environment with `NLsolve` and `LineSearches`.
const N = 64
const NKNOT = 17
const ΔT = 1e-3
const REFERENCE = joinpath(@__DIR__, "..", "data", "lenard_bernstein_metriplectic_reference.txt")

# The velocity domain is the support of the uniform particle cloud, so the projected
# distribution stays positive and the entropy derivative is defined everywhere.
function metriplectic_setup(::Type{T}) where {T}
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    v0 = collect(range(T(-2), T(2); length = N))
    pdist.particles.v[1, :] .= v0
    pdist.particles.w[1, :] .= T(1) / N

    axis = BSplineBasis(UniformMesh(NKNOT - 1, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))

    return MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist)), v0
end

function step_args(v0, mlb, abstol, ::Type{T}) where {T}
    (zeros(T, N), v0, v0, zeros(T, N, 2), 1, zero(T),
        T(ΔT), 3, T(0.5), abstol, 1.0e-50, mlb)
end

# The default call, so that the `maxiters` default of the solve is exercised.
function solve_step(v0, mlb, abstol, ::Type{T}) where {T}
    Picard_iterate_over_particles(step_args(v0, mlb, abstol, T)...)
end

function step_residual(v, v0, mlb, ::Type{T}) where {T}
    r = similar(v)
    VlasovMethods.f!(r, v, v0, (dist = mlb.dist, ent = mlb.entropy), T(ΔT), mlb)
    return r
end

const REFERENCE_V = parse.(Float64, readlines(REFERENCE))

@testset "Metriplectic Lenard–Bernstein solve" begin
    @testset "Float64" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = 1.0e-13
        v = solve_step(v0, mlb, abstol, Float64)
        @test maximum(abs, v .- REFERENCE_V) ≤ abstol
        @test norm(step_residual(v, v0, mlb, Float64)) ≤ abstol
    end

    @testset "Float32" begin
        mlb, v0 = metriplectic_setup(Float32)
        abstol = 1.0e-13 * eps(Float32) / eps(Float64)
        v = solve_step(v0, mlb, abstol, Float32)
        @test maximum(abs, v .- Float32.(REFERENCE_V)) ≤ abstol
        @test norm(step_residual(v, v0, mlb, Float32)) ≤ abstol
    end

    @testset "the solve writes nothing to stdout" begin
        mlb, v0 = metriplectic_setup(Float64)
        mktemp() do path, io
            redirect_stdout(io) do
                solve_step(v0, mlb, 1.0e-13, Float64)
            end
            flush(io)
            @test filesize(path) == 0
        end
    end

    @testset "a non-converging solve throws" begin
        mlb, v0 = metriplectic_setup(Float64)
        @test_throws "residual" Picard_iterate_over_particles(
            step_args(v0, mlb, 1.0e-13, Float64)...; maxiters = 1)
    end
end
