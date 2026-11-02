using ForwardDiff: ForwardDiff
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis
using StableRNGs
using Test
using VlasovMethods

const VM = VlasovMethods

# `N` particles of uniform weight, drawn uniformly on the spline domain `[-2, 2]`, so that the L²
# projection stays positive (KNOWN_ISSUES.md, K8), and a cubic clamped basis on 17 nodes.
function setup(::Type{T}, N) where {T}
    rng = StableRNG(20261102)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    pdist.particles.v[1, :] .= (T(4) .* rand(rng, T, N)) .- T(2)
    pdist.particles.w[1, :] .= one(T) / max(N, 1)
    axis = BSplineBasis(UniformMesh(17, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    return pdist, sdist
end

lb_params(model) = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)

# The bytes of one call after a warm call, through barriers of fixed arity whose arguments have
# concrete types: a right-hand side `f!(v̇, t, q, params)` and a coefficient function
# `f(sdist, pdist, q, values)`.
function rhs_bytes(f!::F, v̇::A, t::B, q::C, p::D) where {F, A, B, C, D}
    f!(v̇, t, q, p)
    return @allocated f!(v̇, t, q, p)
end
function coefficient_bytes(f::F, s::A, p::B, q::C, b::D) where {F, A, B, C, D}
    f(s, p, q, b)
    return @allocated f(s, p, q, b)
end

# The right-hand side of `LenardBernstein` through the allocating `LB_rhs!`, which takes `params`
# without the model, on the barrier of the buffered right-hand sides.
function lb_allocating!(v̇, t, q, params)
    VM.LB_rhs!(v̇, q, (; params.ν, params.idist, params.fdist), t)
end

# A `Dual{Nothing, T, 12}` copy of `q`, as the Jacobian of `ImplicitMidpoint` passes.
function dual12(q::AbstractVector{T}) where {T}
    ForwardDiff.Dual{Nothing}.(q, ntuple(_ -> one(T), 12)...)
end

const MODELS = (
    (LenardBernstein, VM.LB_rhs_GI!), (ConservativeLenardBernstein, VM.CLB_rhs_GI!),
    (RescaledConservativeLenardBernstein, VM.RCLB_rhs_GI!))
const CONSERVATIVE = MODELS[2:3]
const COEFFICIENTS = ((ConservativeLenardBernstein, VM.compute_coefficients),
    (RescaledConservativeLenardBernstein, VM.compute_coefficients_rclb))

@testset "Lenard–Bernstein buffers" begin
    @testset "$(nameof(M)) allocates nothing, $T, N = $N" for (M, rhs!) in MODELS,
        T in (Float32, Float64), N in (256, 2000)
        pdist, sdist = setup(T, N)
        model = M(pdist, CollisionEntropy(sdist))
        params = lb_params(model)
        q = copy(pdist.particles.v[1, :])
        v̇ = similar(q)

        @test @inferred(rhs!(v̇, zero(T), q, params)) === v̇
        # coverage instrumentation allocates on every counted line, so under it the assertion
        # does not run rather than being weakened
        if Base.JLOptions().code_coverage == 0
            @test rhs_bytes(rhs!, v̇, zero(T), q, params) == 0
        end
        # the control: the barrier sees the allocations of the allocating `LB_rhs!`
        @test rhs_bytes(lb_allocating!, similar(q), zero(T), q, params) > 0
    end

    @testset "$(nameof(f)) allocates nothing, $T, N = $N" for (M, f) in COEFFICIENTS,
        T in (Float32, Float64), N in (256, 2000)
        pdist, sdist = setup(T, N)
        model = M(pdist, CollisionEntropy(sdist))
        cache = model.cache[T]
        q = copy(pdist.particles.v[1, :])
        projection(q, pdist, cache.sdist)

        A = @inferred f(cache.sdist, pdist, q, cache.values)
        @test A isa Tuple{T, T}
        @test A == f(cache.sdist, pdist, q)
        if Base.JLOptions().code_coverage == 0
            @test coefficient_bytes(f, cache.sdist, pdist, q, cache.values) == 0
        end
        # the control: the three-argument form allocates its own buffer
        @test coefficient_bytes((s, p, q, _) -> f(s, p, q), cache.sdist, pdist, q,
            cache.values) > 0
    end

    @testset "$(nameof(M)), Dual state, N = $N" for (M, rhs!) in CONSERVATIVE,
        N in (256, 2000)

        pdist, sdist = setup(Float64, N)
        model = M(pdist, CollisionEntropy(sdist))
        params = lb_params(model)
        q = copy(pdist.particles.v[1, :])
        qd = dual12(q)
        v̇d = similar(qd)

        # the first call of the type builds its cache
        @test @inferred(rhs!(v̇d, 0.0, qd, params)) === v̇d
        @test haskey(model.cache.caches, VM._cachehash(eltype(qd)))
        if Base.JLOptions().code_coverage == 0
            @test rhs_bytes(rhs!, v̇d, 0.0, qd, params) == 0
        end
        # the control: the same right-hand side on a copy of the state, which the barrier sees
        @test rhs_bytes((v̇, t, q, p) -> rhs!(v̇, t, copy(q), p), similar(qd), 0.0, qd,
            params) > 0

        # the value parts are the right-hand side of the values. The banded mass solve takes
        # LAPACK for a `Float64` load vector and a generic solve for a `Dual` one, so they agree
        # to the bound of the reference pins of `test/integration/collision_operators.jl`, not
        # bit for bit
        v̇ = similar(q)
        rhs!(v̇, 0.0, q, params)
        @test isapprox(ForwardDiff.value.(v̇d), v̇; rtol = 64 * eps(Float64))
    end

    @testset "the buffer of the three-argument $(nameof(f))" for (M, f) in COEFFICIENTS
        # a `Float32` state on a `Float64` spline: the buffer that the call allocates holds the
        # `Float64` evaluations, as a buffer of the caller does
        pdist, sdist = setup(Float64, 256)
        q = copy(pdist.particles.v[1, :])
        projection(q, pdist, sdist)
        q32 = Float32.(q)
        width = local_width(sdist.basis)
        @test f(sdist, pdist, q32) == f(sdist, pdist, q32, zeros(Float64, width))
        @test f(sdist, pdist, q32) != f(sdist, pdist, q32, zeros(Float32, width))
    end

    @testset "LB_rhs_GI! agrees with LB_rhs!, $T" for T in (Float32, Float64)
        pdist, sdist = setup(T, 256)
        model = LenardBernstein(pdist, CollisionEntropy(sdist))
        params = lb_params(model)
        q = copy(pdist.particles.v[1, :])
        v̇, v̇₀ = similar(q), similar(q)
        VM.LB_rhs_GI!(v̇, zero(T), q, params)
        lb_allocating!(v̇₀, zero(T), q, params)
        @test isapprox(v̇, v̇₀; rtol = 8 * eps(T))
    end

    @testset "$(nameof(M)) projects with the cache's buffers" for (M, rhs!) in MODELS
        # the projection writes the spline of the parent cache, which is the entropy's
        pdist, sdist = setup(Float64, 256)
        model = M(pdist, CollisionEntropy(sdist))
        q = copy(pdist.particles.v[1, :])
        rhs!(similar(q), 0.0, q, lb_params(model))
        cache = parent(model.cache)
        @test cache.sdist === model.ent.dist
        @test length(cache.rhs) == nbasis(sdist.basis)
        @test length(cache.values) == local_width(sdist.basis)
        ref = SplineDistribution(
            1, 1, sdist.basis, sdist.quadrature, zeros(nbasis(sdist.basis)))
        projection(q, pdist, ref)
        @test cache.sdist.coefficients == ref.coefficients
    end

    @testset "no particle and one particle, $(nameof(M))" for (M, rhs!) in MODELS
        for N in (0, 1)
            pdist, sdist = setup(Float64, N)
            model = M(pdist, CollisionEntropy(sdist))
            q = copy(pdist.particles.v[1, :])
            v̇ = similar(q)
            if M === LenardBernstein
                @test isequal(rhs!(v̇, 0.0, q, lb_params(model)),
                    lb_allocating!(similar(q), 0.0, q, lb_params(model)))
                @test length(v̇) == N
            else
                # the 2×2 system of the coefficients is singular without a velocity spread
                phrase = M === ConservativeLenardBernstein ?
                         "the coefficient system is singular" :
                         "the rescaled coefficient system is singular"
                err = try
                    rhs!(v̇, 0.0, q, lb_params(model))
                    nothing
                catch e
                    e
                end
                @test err isa ArgumentError
                @test startswith(err.msg, phrase)
            end
        end
    end
end
