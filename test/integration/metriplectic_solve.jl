using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: SplineQuadrature
using StableRNGs
using Test
using VlasovMethods

# One Picard step of the metriplectic Lenard–Bernstein model. The reference is the result of the
# `NonlinearSolve` solve at commit `8419ed7`; `test/helpers/generate_metriplectic_reference.jl`
# writes it in a scratch environment with `NLsolve` and `LineSearches`.
const N = 64
const NKNOT = 17
const ΔT = 1e-3
const REFERENCE = joinpath(@__DIR__, "..", "data", "lenard_bernstein_metriplectic_reference.txt")

# The example tolerance `3e-16·√N`, the scale `scripts/lenard_bernstein_metriplectic.jl` uses for
# this solve. The `Float32` solve scales it by `eps(Float32)/eps(Float64)`. The script fixes
# `reltol = 1e-50`, so `abstol` is the only active tolerance; the solve passes `reltol` through to
# `SimpleSolvers`.
abstol_float64() = 3.0e-16 * sqrt(N)
abstol_float32() = 3.0e-16 * sqrt(N) * eps(Float32) / eps(Float64)

# The velocity domain is the support of the uniform particle cloud, so the projected
# distribution stays positive and the entropy derivative is defined everywhere.
function metriplectic_setup(::Type{T}, n::Int = N) where {T}
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, n); variables = vars))
    v0 = collect(range(T(-2), T(2); length = n))
    pdist.particles.v[1, :] .= v0
    pdist.particles.w[1, :] .= T(1) / n

    axis = BSplineBasis(UniformMesh(NKNOT - 1, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))

    return MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist)), v0
end

function step_args(v0, mlb, abstol, ::Type{T};
        ti::Int = 1, vn_minus_one = v0, dv_history = zeros(T, N, 2)) where {T}
    (zeros(T, N), v0, vn_minus_one, dv_history, ti,
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

# The in-place step at `t = 0`, from `vn` with a zero history by default.
function step!(v, vn, mlb, abstol, ::Type{T}; ti::Int = 4, reltol = 1.0e-50,
        kwargs...) where {T}
    Picard_iterate_over_particles!(v, vn, vn, zeros(T, length(vn), 2), ti, zero(T), T(ΔT),
        abstol, reltol, mlb; kwargs...)
end

# The exception the in-place step throws, or `nothing` when it returns.
function step_outcome(args...; kwargs...)
    try
        Picard_iterate_over_particles!(args...; kwargs...)
        return nothing
    catch e
        return e
    end
end

# `@allocated` through a fixed-arity barrier whose arguments have concrete types, so that a Julia
# 1.11 closure or splat boxing is not counted. The call runs once before it is measured.
function step_bytes(f::F, v::V, vn::V, vm::V, H::HT, ti::Int, t::S, Δt::S, a::A, r::R,
        mlb::M) where {F, V, HT, S, A, R, M}
    f(v, vn, vm, H, ti, t, Δt, a, r, mlb)
    return @allocated f(v, vn, vm, H, ti, t, Δt, a, r, mlb)
end

# The wrapper behind the signature of the in-place step: the control, which allocates its result.
function wrapper_step!(v, vn, vm, H, ti, t, Δt, a, r, mlb)
    copyto!(v, Picard_iterate_over_particles(v, vn, vm, H, ti, t, Δt, 3, 0.5, a, r, mlb))
end

# The seeded state of `test/helpers/generate_metriplectic_trajectory.jl`: the cell centres of
# [-2, 2], each moved by a draw uniform in ±1/64, with weights 1/N.
function seeded_setup()
    mlb, _ = metriplectic_setup(Float64)
    v0 = [-2 + (i - 0.5) * 4 / N for i in 1:N] .+
         (rand(StableRNG(20261010), N) .- 0.5) ./ 32
    mlb.dist.particles.v[1, :] .= v0
    return mlb, v0
end

# Ten steps, `ti = 1:10` at `t = (ti - 1)·Δt`, from `vn = vn_minus_one = v0` and a zero history;
# after each step `vn_minus_one` takes `vn`, and `vn` the result of `step(vn, vm, H, ti, t)`.
function march(step, v0)
    states = Vector{Float64}[]
    vn, vm, H = copy(v0), copy(v0), zeros(N, 2)
    for ti in 1:10
        v = step(vn, vm, H, ti, (ti - 1) * ΔT)
        push!(states, copy(v))
        vm, vn = vn, v
    end
    return states
end

const TRAJECTORY = let
    stored = parse.(Float64,
        readlines(joinpath(@__DIR__, "..", "data",
            "lenard_bernstein_metriplectic_trajectory.txt")))
    [stored[((k - 1) * N + 1):(k * N)] for k in 1:(length(stored) ÷ N)]
end

@testset "Metriplectic Lenard–Bernstein solve" begin
    @testset "Float64" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        v = @inferred Picard_iterate_over_particles(step_args(v0, mlb, abstol, Float64)...)
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

    # The Hermite guess of the next step reads both columns of `dv_history`: column 1 must hold
    # the field at the solved iterate, and column 2 the column 1 that the step received.
    @testset "the stored derivative is the field at the solved iterate" begin
        mlb, v0 = metriplectic_setup(Float64)
        dv_history = zeros(Float64, N, 2)
        previous = collect(range(1.0, 2.0; length = N))
        dv_history[:, 1] .= previous
        v = Picard_iterate_over_particles(
            step_args(v0, mlb, abstol_float64(), Float64; dv_history = dv_history)...)
        expected = similar(v)
        VlasovMethods.collisional_vectorfield!(
            expected, v, (dist = mlb.dist, ent = mlb.entropy), mlb)
        @test dv_history[:, 1] == expected
        @test dv_history[:, 2] == previous
    end

    # The `ti ≥ 4` Hermite guess reads `vn_minus_one`, but the step's residual depends on `vn`
    # alone: a perturbed `vn_minus_one` changes the starting guess and not the fixed point.
    @testset "the ti ≥ 4 solve uses vn, not vn_minus_one, in the residual" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        va = Picard_iterate_over_particles(step_args(v0, mlb, abstol, Float64; ti = 4)...)
        vb = Picard_iterate_over_particles(
            step_args(v0, mlb, abstol, Float64; ti = 4, vn_minus_one = (1 - 1e-6) .* v0)...)
        @test maximum(abs, va .- REFERENCE_V) ≤ abstol
        @test maximum(abs, va .- vb) ≤ 10abstol
        @test norm(step_residual(vb, v0, mlb, Float64)) ≤ abstol
    end

    # `maxiters = 2` so the reported count is not the `1` the NaN test also reports: a message
    # that hard-codes `1` passes the NaN test but fails this one.
    @testset "a non-converging solve throws the residual and the iteration count" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        caught = solve_outcome(step_args(v0, mlb, abstol, Float64)...; maxiters = 2)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        captured = match(r"residual = ([0-9.eE+-]+), iterations = ([0-9]+)", message)
        @test captured !== nothing
        @test captured.captures[2] == "2"
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

    # The stored states are those of the solve at commit `472dd42`, which
    # `test/helpers/generate_metriplectic_trajectory.jl` writes. The pin is `==`: a change of the
    # arithmetic order in the vector field, the projection or the solve moves the last bits.
    @testset "ten wrapper steps of a seeded state give the stored trajectory" begin
        @test length(TRAJECTORY) == 10
        mlb, v0 = seeded_setup()
        abstol = abstol_float64()
        wrapper(vn, vm, H, ti, t) = Picard_iterate_over_particles(
            zeros(N), vn, vm, H, ti, t, ΔT, 3, 0.5, abstol, 1.0e-50, mlb)
        @test march(wrapper, v0) == TRAJECTORY
    end

    @testset "ten in-place steps of a seeded state give the stored trajectory" begin
        mlb, v0 = seeded_setup()
        abstol = abstol_float64()
        inplace(vn, vm, H, ti, t) = Picard_iterate_over_particles!(
            similar(vn), vn, vm, H, ti, t, ΔT, abstol, 1.0e-50, mlb)
        @test march(inplace, v0) == TRAJECTORY
    end

    # The start-up steps `ti < 4` build an `ODEProblem` for their guess and allocate; the claim
    # covers the Hermite steps `ti ≥ 4`. Each precision and length has its own model.
    @testset "the in-place step allocates nothing at ti ≥ 4" begin
        if Base.JLOptions().code_coverage == 0
            for T in (Float32, Float64), n in (64, 256)

                mlb, v0 = metriplectic_setup(T, n)
                abstol = 3.0e-16 * sqrt(n) * eps(T) / eps(Float64)
                @test step_bytes(Picard_iterate_over_particles!, similar(v0), v0, copy(v0),
                    zeros(T, n, 2), 4, zero(T), T(ΔT), abstol, 1.0e-50, mlb) == 0
                @test step_bytes(wrapper_step!, similar(v0), v0, copy(v0),
                    zeros(T, n, 2), 4, zero(T), T(ΔT), abstol, 1.0e-50, mlb) > 0
            end
        end
    end

    @testset "the in-place step and the wrapper give the same step" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        for ti in (1, 4)
            v = step!(similar(v0), v0, mlb, abstol, Float64; ti = ti)
            @test v != v0
            @test v == Picard_iterate_over_particles(
                step_args(v0, mlb, abstol, Float64; ti = ti)...)
        end
    end

    @testset "the in-place step in each precision" begin
        for (T, abstol) in ((Float64, abstol_float64()), (Float32, abstol_float32()))
            mlb, v0 = metriplectic_setup(T)
            v = step!(similar(v0), v0, mlb, abstol, T; ti = 1)
            @test eltype(v) === T
            @test maximum(abs, v .- T.(REFERENCE_V)) ≤ abstol
        end
    end

    # The solver of a model is built once and reused; a call with other options rebuilds it, and a
    # call with the first options again rebuilds it back.
    @testset "a change of maxiters rebuilds the solver of the in-place step" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        v_default = step!(similar(v0), v0, mlb, abstol, Float64)
        caught = step_outcome(
            similar(v0), v0, v0, zeros(N, 2), 4, 0.0, ΔT, abstol, 1.0e-50,
            mlb; maxiters = 2)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        @test occursin("did not converge", message)
        @test occursin("iterations = 2", message)
        @test step!(similar(v0), v0, mlb, abstol, Float64) == v_default
    end

    # A new model holds a solver of length 1, built at `abstol = reltol = 0` and `maxiters = 1`
    # for its type. A first call with those options builds a solver of its own length.
    @testset "the first call builds a solver of its own length" begin
        mlb, v0 = metriplectic_setup(Float64)
        caught = step_outcome(similar(v0), v0, v0, zeros(N, 2), 4, 0.0, ΔT, 0.0, 0.0, mlb;
            maxiters = 1)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        @test occursin("did not converge", message)
        @test occursin("iterations = 1", message)
    end

    @testset "a change of abstol or reltol rebuilds the solver of the in-place step" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        v_default = step!(similar(v0), v0, mlb, abstol, Float64)
        @test step!(similar(v0), v0, mlb, 1e-6, Float64) != v_default
        @test step!(similar(v0), v0, mlb, abstol, Float64) == v_default
        # `SimpleSolvers` applies `f_reltol` only once the iterate or the residual stops
        # changing, so a looser `reltol` alone leaves this step as it is. At `abstol = 0` the
        # solve converges through `reltol = 1e-6` and not at `reltol = 1e-50`.
        @test step_outcome(similar(v0), v0, v0, zeros(N, 2), 4, 0.0, ΔT, 0.0, 1e-6, mlb) ===
              nothing
        caught = step_outcome(
            similar(v0), v0, v0, zeros(N, 2), 4, 0.0, ΔT, 0.0, 1.0e-50, mlb)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        @test occursin("did not converge", message)
    end

    @testset "the in-place step refuses a wrong length and an aliased output" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        vn, vm, H = copy(v0), copy(v0), zeros(N, 2)
        outcome(v) = step_outcome(v, vn, vm, H, 4, 0.0, ΔT, abstol, 1.0e-50, mlb)

        caught = outcome(zeros(N - 1))
        @test caught isa DimensionMismatch
        @test caught !== nothing &&
              caught.msg == "the output `v` has 63 entries and `vn` has 64"
        caught = outcome(vn)
        @test caught isa ArgumentError
        @test caught !== nothing && caught.msg == "the output `v` shares memory with `vn`"
        caught = outcome(vm)
        @test caught isa ArgumentError
        @test caught !== nothing &&
              caught.msg == "the output `v` shares memory with `vn_minus_one`"
        caught = outcome(view(H, :, 1))
        @test caught isa ArgumentError
        @test caught !== nothing &&
              caught.msg == "the output `v` shares memory with `dv_history`"
    end

    # A model holds buffers of its particle count and does not rebuild them to another length.
    @testset "the in-place step refuses a vn of another length than the model's" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        v_first = step!(similar(v0), v0, mlb, abstol, Float64)
        w = v0[1:(N - 1)]
        caught = step_outcome(
            zeros(N - 1), w, copy(w), zeros(N - 1, 2), 4, 0.0, ΔT, abstol,
            1.0e-50, mlb)
        @test caught isa DimensionMismatch
        @test caught !== nothing &&
              caught.msg == "`vn` has 63 entries and the model has 64 particles"
        @test step!(similar(v0), v0, mlb, abstol, Float64) == v_first
    end

    @testset "the in-place step passes the failures of the solve to the caller" begin
        mlb, v0 = metriplectic_setup(Float64)
        abstol = abstol_float64()
        nanmlb = MetriplecticLenardBernstein(mlb.dist, mlb.entropy; ν = NaN)
        caught = step_outcome(
            similar(v0), v0, v0, zeros(N, 2), 4, 0.0, ΔT, abstol, 1.0e-50,
            nanmlb)
        @test caught isa ErrorException
        message = caught === nothing ? "" : sprint(showerror, caught)
        @test occursin("did not converge", message)
        @test occursin("iterations = 1", message)
        fastmlb = MetriplecticLenardBernstein(mlb.dist, mlb.entropy; ν = 1e6)
        caught = step_outcome(
            similar(v0), v0, v0, zeros(N, 2), 4, 0.0, ΔT, abstol, 1.0e-50,
            fastmlb)
        @test caught isa DomainError
    end
end
