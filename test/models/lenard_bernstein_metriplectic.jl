using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis
using Test
using VlasovMethods

const VM = VlasovMethods

# The non-uniform weights `w_γ ∝ 1 + γ/N` of `N` particles, normalised so that `n = Σ w = 1.3`.
function fixture_weights(::Type{T}, N) where {T}
    a = 1 .+ (1:N) ./ N
    return T.(1.3 .* a ./ sum(a))
end

# Twelve particles of non-uniform weight at the velocities `v_γ = −2 + 4(γ − 1)/(N − 1) + 0.1 sin γ`,
# and an arbitrary entropy gradient `S_γ = cos(2γ) w_γ`.
function fixture(::Type{T}) where {T}
    N = 12
    γ = 1:N
    w = fixture_weights(T, N)
    v = T.(-2 .+ 4 .* (γ .- 1) ./ (N - 1) .+ 0.1 .* sin.(γ))
    dS = T.(cos.(2 .* γ)) .* w
    return w, v, dS
end

# A model of `length(w)` particles with the weights `w`. The right-hand side reads the weights from
# the particle distribution and the velocities from its argument, so the particles' own velocities
# stay zero; the spline distribution is not read.
function model(::Type{T}, w; kwargs...) where {T}
    N = length(w)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    pdist.particles.w[1, :] .= w
    axis = BSplineBasis(UniformMesh(16, T(-3) .. T(3)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    return MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist); kwargs...)
end

function rhs(mlb, v, dS)
    v̇ = similar(v)
    VM.metriplectic_rhs!(v̇, v, dS, mlb)
    return v̇
end

# The moments `n`, `u` and `ε` of the weights `w` at the velocities `v`.
function moments(w, v)
    n = sum(w)
    u = sum(w .* v) / n
    ε = sum(w .* v .^ 2) / n
    return n, u, ε
end

# The matrix `𝕃` of `v̇ = 𝕃 dS` at the velocities `v`, column by column from the unit vectors.
function bracket_matrix(f, v::AbstractVector{T}) where {T}
    N = length(v)
    L = zeros(T, N, N)
    for β in 1:N
        e = zeros(T, N)
        e[β] = one(T)
        L[:, β] .= f(e)
    end
    return L
end

relerr_max(a, b) = maximum(abs, a .- b) / maximum(abs, b)

@testset "Lenard–Bernstein brackets" begin
    @test :LenardBernsteinBracket in names(VlasovMethods)
    @test :DoughertyBracket in names(VlasovMethods)
    @test :JeyakumarBracket in names(VlasovMethods)
    @test DoughertyBracket() isa LenardBernsteinBracket
    @test JeyakumarBracket() isa LenardBernsteinBracket
    @test Base.issingletontype(DoughertyBracket)
    @test Base.issingletontype(JeyakumarBracket)
    for name in (:LenardBernsteinBracket, :DoughertyBracket, :JeyakumarBracket)
        @test Base.Docs.hasdoc(VlasovMethods, name)
    end
end

@testset "Dougherty term by term ($T)" for T in (Float64, Float32)
    w, v, dS = fixture(T)
    ν = T(0.7)
    mlb = model(T, w; ν = ν, bracket = DoughertyBracket())

    # The three terms of `eq:semidiscrete_equations_of_motion`, each particle with its own weight.
    n, u, ε = moments(w, v)
    expected = ν .* (.-(ε - u^2) .* dS ./ w .+ (ε .- u .* v) .* sum(dS) ./ n .+
                (v .- u) .* sum(v .* dS) ./ n)

    @test relerr_max(rhs(mlb, v, dS), expected) ≤ 8 * eps(T)
end

@testset "Jeyakumar is the base's form ($T)" for T in (Float64, Float32)
    w, v, dS = fixture(T)
    ν = T(0.7)
    mlb = model(T, w; ν = ν, bracket = JeyakumarBracket())

    # The Jeyakumar form written out, with the moments, sums and weights computed as the model
    # computes them, so that the comparison is exact.
    pw = mlb.dist.particles.w
    n = sum(pw)
    u = dot(pw, v) / n
    eps = dot(pw, v .^ 2) / n
    dS_sum = sum(dS)
    dS_v_sum = dot(v, dS)
    w = view(pw, 1, :)
    expected = ν .* (.-n ./ w .* dS .+ (eps .- u .* v) ./ (eps - u^2) .* dS_sum .+
                (v .- u) ./ (eps - u^2) .* dS_v_sum)

    @test rhs(mlb, v, dS) == expected
end

@testset "the factor between the forms ($T)" for T in (Float64, Float32)
    w, v, dS = fixture(T)
    ν = T(0.7)
    v̇_D = rhs(model(T, w; ν = ν, bracket = DoughertyBracket()), v, dS)
    v̇_J = rhs(model(T, w; ν = ν, bracket = JeyakumarBracket()), v, dS)

    n, u, ε = moments(w, v)
    @test norm(v̇_J .- n / (ε - u^2) .* v̇_D) ≤ 64 * eps(T) * norm(v̇_J)
end

@testset "the default bracket" begin
    # Particles that fill the spline domain, so that the projection stays positive.
    T = Float64
    N = 64
    w = fill(one(T) / N, N)
    v = collect(range(T(-2.9), T(2.9); length = N))
    mlb_default = model(T, w)
    mlb_dougherty = model(T, w; bracket = DoughertyBracket())
    mlb_jeyakumar = model(T, w; bracket = JeyakumarBracket())
    @test mlb_default.bracket isa DoughertyBracket

    vf(mlb) = (v̇ = similar(v); VM.collisional_vectorfield!(v̇, v, nothing, mlb); v̇)
    @test vf(mlb_default) == vf(mlb_dougherty)
    @test vf(mlb_default) != vf(mlb_jeyakumar)
end

@testset "the bracket matrix ($T)" for T in (Float64, Float32)
    w, v, _ = fixture(T)
    n, u, ε = moments(w, v)
    scale = max(sum(w), sum(w .* abs.(v)))

    @testset "$(nameof(typeof(bracket)))" for bracket in (DoughertyBracket(), JeyakumarBracket())
        mlb = model(T, w; ν = one(T), bracket = bracket)
        L = bracket_matrix(dS -> rhs(mlb, v, dS), v)
        Lmax = maximum(abs, L)
        tol = 16 * eps(T) * Lmax

        @test maximum(abs, L .- L') ≤ tol
        @test maximum(eigvals(Symmetric((L .+ L') ./ 2))) ≤ tol
        @test maximum(abs, w' * L) ≤ tol * scale
        @test maximum(abs, (w .* v)' * L) ≤ tol * scale
        # The bracket is not zero: its symmetric part has an eigenvalue below `−max|𝕃|/10`.
        @test minimum(eigvals(Symmetric((L .+ L') ./ 2))) < -Lmax / 10
    end

    # The control: the Dougherty right-hand side without its `(ε − u v_γ) Σ S / n` term, whose
    # matrix has no `w` in its left kernel, so that the check above can fail.
    control(dS) = .-(ε - u^2) .* dS ./ w .+ (v .- u) .* sum(v .* dS) ./ n
    L_c = bracket_matrix(control, v)
    @test maximum(abs, w' * L_c) ≥ 0.05 * maximum(abs, L_c)
end

@testset "a degenerate velocity set throws ($T)" for T in (Float64, Float32)
    @testset "$(nameof(typeof(bracket)))" for bracket in (DoughertyBracket(), JeyakumarBracket())
        for N in (1, 2, 3, 12, 50, 64, 128, 256)
            for w in (fill(one(T) / N, N), fixture_weights(T, N))
                mlb = model(T, w; bracket = bracket)
                dS = ones(T, N)
                thrown = count(range(T(-3), T(3); length = 61)) do v₀
                    v = fill(v₀, N)
                    try
                        rhs(mlb, v, dS)
                        false
                    catch e
                        e isa DomainError || rethrow()
                        occursin("ε − u²", sprint(showerror, e))
                    end
                end
                @test thrown == 61
            end
        end

        # The velocities of the term-by-term fixture are spread, and do not throw.
        w, v, dS = fixture(T)
        @test all(isfinite, rhs(model(T, w; bracket = bracket), v, dS))

        # A `NaN` velocity is not refused: it propagates into the right-hand side.
        v_nan = copy(v)
        v_nan[1] = T(NaN)
        @test all(isnan, rhs(model(T, w; bracket = bracket), v_nan, dS))
    end
end

# `rhs(mlb, v, dS)`, or `nothing` where it throws the `DomainError` that names `ε − u²`.
function rhs_or_nothing(mlb, v, dS)
    try
        return rhs(mlb, v, dS)
    catch e
        e isa DomainError && occursin("ε − u²", sprint(showerror, e)) || rethrow()
        return nothing
    end
end

@testset "the threshold" begin
    # Two particles of weight 1/2 at `1 ± δ`, a dyadic `δ`, so that `u == 1` and the computed
    # `ε − u²` equals `σ² = δ²` exactly. Below `eps(T) ε` the variance clause alone throws.
    cases = (
        (Float64, 2.0^-26, true), (Float64, 2.0^-25, false), (Float32, 2.0f0^-11, false))
    @testset "$(nameof(typeof(bracket)))" for bracket in (DoughertyBracket(), JeyakumarBracket())
        for (T, δ, throws) in cases
            w = fill(one(T) / 2, 2)
            v = [one(T) - δ, one(T) + δ]
            _, u, ε = moments(w, v)
            @test u == 1
            @test ε - u^2 == δ^2
            v̇ = rhs_or_nothing(model(T, w; bracket = bracket), v, ones(T, 2))
            if throws
                @test v̇ === nothing
            else
                @test v̇ !== nothing && all(isfinite, v̇)
            end
        end
    end
end

@testset "a computed ε − u² ≤ 0 throws ($T)" for T in (Float64, Float32)
    # Twelve particles at `c ± δ`, alternating, with `σ² = δ²` a small multiple of `eps(T) ε`,
    # where the one-pass `ε − u²` of the model can be zero or negative.
    N = 12
    w = fill(one(T) / N, N)
    dS = ones(T, N)
    @testset "$(nameof(typeof(bracket)))" for bracket in (DoughertyBracket(), JeyakumarBracket())
        mlb = model(T, w; bracket = bracket)
        thrown_above = 0
        for c in range(T(0.5), T(2.9); length = 25), f in (1.1, 1.5, 2, 3, 4)

            δ = sqrt(T(f) * eps(T)) * c
            v = [c + (-1)^γ * δ for γ in 1:N]
            n, u, ε = VM.compute_moments(v, similar(v), mlb.dist, mlb)
            σ² = sum(w .* (v .- u) .^ 2) / n
            v̇ = rhs_or_nothing(mlb, v, dS)
            if v̇ === nothing
                thrown_above += σ² > eps(T) * ε
            else
                @test ε - u^2 > 0
                @test all(isfinite, v̇)
            end
        end
        @test thrown_above ≥ 1
    end
end
