using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: SplineQuadrature
using Test
using VlasovMethods

# One Picard step of the metriplectic Lenard–Bernstein model. The reference is the result of the
# `NonlinearSolve` solve at commit `8419ed7`; `test/helpers/generate_metriplectic_reference.jl`
# writes it in a scratch environment with `NLsolve` and `LineSearches`.
const N = 64
const NKNOT = 17
const ΔT = 1e-3
const REFERENCE = joinpath(@__DIR__, "..", "data", "lenard_bernstein_metriplectic_reference.txt")

# The example tolerance `3e-16·√N`, the scale the scaling script uses for this solve. The
# `Float32` solve scales it by `eps(Float32)/eps(Float64)`.
abstol_float64() = 3.0e-16 * sqrt(N)
abstol_float32() = 3.0e-16 * sqrt(N) * eps(Float32) / eps(Float64)

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

function step_args(v0, mlb, abstol, ::Type{T};
        ti::Int = 1, dv_history = zeros(T, N, 2)) where {T}
    (zeros(T, N), v0, v0, dv_history, ti,
        zero(T), T(ΔT), 3, T(0.5), abstol, 1.0e-50, mlb)
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

# The exception the solve throws, or `nothing` when it returns.
function solve_outcome(args...; kwargs...)
    try
        Picard_iterate_over_particles(args...; kwargs...)
        return nothing
    catch e
        return e
    end
end

const REFERENCE_V = parse.(Float64, readlines(REFERENCE))

@testset "Metriplectic Lenard–Bernstein solve" begin
    @testset "Float64" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        v = solve_step(v0, mlb, abstol, Float64)
        @test maximum(abs, v .- REFERENCE_V) ≤ abstol
        @test norm(step_residual(v, v0, mlb, Float64)) ≤ abstol
    end

    @testset "Float32" begin
        mlb, v0 = metriplectic_setup(Float32)
        abstol = abstol_float32()
        v = solve_step(v0, mlb, abstol, Float32)
        @test maximum(abs, v .- Float32.(REFERENCE_V)) ≤ abstol
        @test norm(step_residual(v, v0, mlb, Float32)) ≤ abstol
    end

    @testset "the solve writes nothing to stdout" begin
        mlb, v0 = metriplectic_setup(Float64)
        mktemp() do path, io
            redirect_stdout(io) do
                solve_step(v0, mlb, abstol_float64(), Float64)
            end
            flush(io)
            @test filesize(path) == 0
        end
    end

    # The Hermite guess of the next step reads `dv_history[:, 1]`, so it must hold the field at
    # the solved iterate and not at the previous time step's velocities.
    @testset "the stored derivative is the field at the solved iterate" begin
        mlb, v0 = metriplectic_setup(Float64)
        dv_history = zeros(Float64, N, 2)
        v = Picard_iterate_over_particles(
            step_args(v0, mlb, abstol_float64(), Float64; dv_history = dv_history)...)
        expected = similar(v)
        VlasovMethods.collisional_vectorfield!(
            expected, v, (dist = mlb.dist, ent = mlb.entropy), mlb)
        @test dv_history[:, 1] == expected
    end

    @testset "a non-converging solve throws the residual and the iteration count" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        caught = solve_outcome(step_args(v0, mlb, abstol, Float64)...; maxiters = 1)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        captured = match(r"residual = ([0-9.eE+-]+), iterations = ([0-9]+)", message)
        @test captured !== nothing
        @test captured.captures[2] == "1"
        @test parse(Float64, captured.captures[1]) > abstol
    end

    # A non-finite collision frequency keeps the projected density positive while the residual
    # becomes `NaN`, so the first non-finite value is met inside the `SimpleSolvers` solve.
    @testset "a solve that meets a NaN throws the residual and the iteration count" begin
        mlb, v0 = metriplectic_setup(Float64)
        nanmlb = MetriplecticLenardBernstein(mlb.dist, mlb.entropy; ν = NaN)
        caught = solve_outcome(step_args(v0, nanmlb, abstol_float64(), Float64; ti = 4)...)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        captured = match(r"residual = ([^,]+), iterations = ([0-9]+)", message)
        @test captured !== nothing
        @test captured.captures[2] == "1"
        @test isnan(parse(Float64, captured.captures[1]))
    end

    # The catch in `Picard_iterate_over_particles` rethrows every exception that is not a
    # `NonlinearSolverException`. A large collision frequency drives the particles out of the
    # velocity support, and the projection then throws its own `DomainError` from inside the
    # solve; it must reach the caller unchanged.
    @testset "an exception from the model reaches the caller" begin
        mlb, v0 = metriplectic_setup(Float64)
        fastmlb = MetriplecticLenardBernstein(mlb.dist, mlb.entropy; ν = 1e6)
        caught = solve_outcome(step_args(v0, fastmlb, abstol_float64(), Float64; ti = 4)...)
        @test caught isa DomainError
    end
end
