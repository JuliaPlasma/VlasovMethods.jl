using ForwardDiff: ForwardDiff
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis,
                     TensorProductBasis, TensorProductQuadrature, evaluate!
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

# 24 particles of uniform weight drawn uniformly on `[-2, 2]` with the seed `0x5EED`, on the cubic
# basis of 17 knots of `metriplectic_setup` in `test/integration/metriplectic_solve.jl`. The
# projection of so few particles undershoots below zero between them.
function mlb_setup(::Type{T}) where {T}
    N = 24
    rng = StableRNG(0x5EED)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    pdist.particles.v[1, :] .= T.(4 .* rand(rng, N) .- 2)
    pdist.particles.w[1, :] .= one(T) / N
    axis = BSplineBasis(UniformMesh(16, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    return MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist)), pdist, sdist
end

# The work arrays of `compute_J!` of `MetriplecticLenardBernstein`, and the call
function mlb_J!(J, sdist, mlb)
    g = zeros(eltype(J), length(quadrature_nodes(sdist.quadrature)))
    values = zeros(eltype(J), local_width(sdist.basis))
    return VM.compute_J!(J, g, values, sdist, mlb)
end
mlb_entropy(sdist, mlb) = VM.compute_entropy(y -> sdist.spline(y), mlb)

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

# The five particle sums of a coefficient function in `BigFloat`, over the same terms in the type
# of the state: the weights, the velocities and the values of the spline at them. The sums are
# `n`, `nu`, `neps` and the two sums of `w g` and `w v g` with `g = f_s' / f_s`.
function particle_sums(sdist, pdist, q::AbstractVector{VT}, values) where {VT}
    w = pdist.particles.w
    n = nu = neps = G1 = G2 = big(0.0)
    for α in eachindex(q)
        wα, vα = big(VT(w[1, α])), big(q[α])
        g = big(evaluate!(values, sdist.spline, q[α], 1) /
                evaluate!(values, sdist.spline, q[α]))
        n += wα
        nu += wα * vα
        neps += wα * vα^2
        G1 += wα * g
        G2 += wα * vα * g
    end
    return n, nu, neps, G1, G2
end

# The coefficients of the two systems from those sums, as the functions solve them
function conservative_solve(n, nu, neps, G1, G2)
    Δ = n * neps - nu^2
    return (nu * G2 - neps * G1) / Δ, (nu * G1 - n * G2) / Δ
end
function rescaled_solve(n, nu, neps, G1, G2)
    den = nu * G1 - n * G2
    return (n * neps - nu^2) / den, (nu * G2 - neps * G1) / den
end

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

    @testset "$(nameof(f)) sums compensated, $T" for (M, f, solve) in (
            (ConservativeLenardBernstein, VM.compute_coefficients, conservative_solve),
            (RescaledConservativeLenardBernstein, VM.compute_coefficients_rclb,
            rescaled_solve)),
        T in (Float32, Float64)
        # 2000 uniform weights of 1/2000: a plain sum of a constant term loses the same rounding
        # error at every step, and the coefficients carry that error, 2e-5 in `Float32` and
        # 6e-14 in `Float64`. The compensated sums stay within a few ulp of the exact ones.
        pdist, sdist = setup(T, 2000)
        model = M(pdist, CollisionEntropy(sdist))
        cache = model.cache[T]
        q = copy(pdist.particles.v[1, :])
        projection(q, pdist, cache.sdist)

        A = @inferred f(cache.sdist, pdist, q, cache.values)
        @test A isa Tuple{T, T}
        R = solve(particle_sums(cache.sdist, pdist, q, cache.values)...)
        @test all(i -> abs(A[i] - R[i]) <= 8 * eps(T) * abs(R[i]), 1:2)
        if Base.JLOptions().code_coverage == 0
            @test coefficient_bytes(f, cache.sdist, pdist, q, cache.values) == 0
        end
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
        @test haskey(model.cache.caches, eltype(qd))
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

    @testset "$(nameof(f)), a Float32 state on a Float64 model" for (M, f) in COEFFICIENTS
        # the sums run in the type of the state, with the weights converted to it
        pdist, sdist = setup(Float64, 256)
        model = M(pdist, CollisionEntropy(sdist))
        cache = model.cache[Float32]
        q32 = Float32.(pdist.particles.v[1, :])
        projection(q32, pdist, cache.sdist)

        A = @inferred f(cache.sdist, pdist, q32, cache.values)
        @test A isa Tuple{Float32, Float32}
        if Base.JLOptions().code_coverage == 0
            @test coefficient_bytes(f, cache.sdist, pdist, q32, cache.values) == 0
        end
    end

    @testset "$(nameof(rhs!)), a Float32 state on a Float64 model" for (M, rhs!) in CONSERVATIVE
        pdist, sdist = setup(Float64, 256)
        model = M(pdist, CollisionEntropy(sdist))
        params = lb_params(model)
        q32 = Float32.(pdist.particles.v[1, :])
        v̇ = similar(q32)

        @test @inferred(rhs!(v̇, 0.0f0, q32, params)) === v̇
        @test eltype(v̇) === Float32
        if Base.JLOptions().code_coverage == 0
            @test rhs_bytes(rhs!, v̇, 0.0f0, q32, params) == 0
        end
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

    @testset "$(nameof(M)) builds its cache on a tensor-product basis" for (M, _) in MODELS
        # the buffers of a tensor-product basis are one vector per axis, not one vector
        T = Float64
        vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
        pdist = ParticleDistribution(
            1, 2, ParticleList(zeros(T, 4, 16); variables = vars))
        axis = BSplineBasis(UniformMesh(5, T(-2) .. T(2)), 2, Free())
        B = TensorProductBasis((axis, axis))
        sdist = SplineDistribution(
            1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
        model = M(pdist, CollisionEntropy(sdist))
        cache = model.cache[T]
        @test cache.values isa NTuple{2, Vector{T}}
        @test length.(cache.values) == ntuple(_ -> local_width(axis), 2)
        @test length(cache.rhs) == nbasis(B)
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

# ∂S_h/∂v_α of `compute_dS!`, from 𝕃 of `compute_J!`, is the gradient of S_h = ∫ f_s log|f_s|
# (`eq:disc_entropy`) of the projected particles, against a central difference of
# `compute_entropy`.
@testset "MP-LB gradient on a sign-changing projection" begin
    mlb, pdist, sdist = mlb_setup(Float64)
    v = copy(pdist.particles.v[1, :])
    projection(v, pdist, sdist)

    # f_s changes sign at the nodes, and no node value is close enough to zero for the
    # difference step to cross the singularity of log|f_s|
    fnodes = [sdist.spline(x) for x in quadrature_nodes(sdist.quadrature)]
    @test minimum(fnodes) < 0
    @test minimum(abs, fnodes) > 1e-8

    J = mlb_J!(zeros(nbasis(sdist.basis)), sdist, mlb)
    dS = VM.compute_dS!(
        similar(v), zeros(local_width(sdist.basis)), J, v, sdist, mlb, pdist)

    h = 1e-6
    entropy_at(u) = (projection(u, pdist, sdist); mlb_entropy(sdist, mlb))
    central = map(eachindex(v)) do α
        vp, vm = copy(v), copy(v)
        vp[α] += h
        vm[α] -= h
        (entropy_at(vp) - entropy_at(vm)) / 2h
    end
    @test maximum(abs, central - dS) ≤ 1e-6 * maximum(abs, dS)
end

@testset "compute_entropy of a sign-changing projection, $T" for T in (Float32, Float64)
    mlb, pdist, sdist = mlb_setup(T)
    projection(copy(pdist.particles.v[1, :]), pdist, sdist)
    x = quadrature_nodes(sdist.quadrature)
    w = quadrature_weights(sdist.quadrature)
    f = [sdist.spline(xr) for xr in x]
    @test minimum(f) < 0
    terms = [w[r] * f[r] * log(abs(f[r])) for r in eachindex(x)]
    S = mlb_entropy(sdist, mlb)
    @test S isa T
    @test abs(S - sum(terms)) ≤ 8 * eps(T) * sum(abs, terms)
end

# log|f_s| is -Inf at a zero of f_s and undefined at a non-finite value, so 𝕃 and the entropy
# refuse both.
@testset "a zero projection throws, $T" for T in (Float32, Float64)
    mlb, pdist, sdist = mlb_setup(T)
    sdist.coefficients .= 0
    @test_throws DomainError mlb_J!(zeros(T, nbasis(sdist.basis)), sdist, mlb)
    @test_throws DomainError mlb_entropy(sdist, mlb)
end

@testset "a non-finite projection throws, $T" for T in (Float32, Float64)
    mlb, pdist, sdist = mlb_setup(T)
    for bad in (T(NaN), T(Inf))
        sdist.coefficients .= one(T)
        sdist.coefficients[5] = bad
        @test_throws DomainError mlb_J!(zeros(T, nbasis(sdist.basis)), sdist, mlb)
        @test_throws DomainError mlb_entropy(sdist, mlb)
    end
end

# PD-LB's pointwise ratio f_s'/f_s has no meaning where f_s ≤ 0, so its coefficients refuse a
# projection that is negative at a particle.
@testset "PD-LB refuses a sign-changing projection, $T" for T in (Float32, Float64)
    pdist, sdist = setup(T, 256)
    project_function(v -> v^2 - one(T) / 4, sdist)
    q = copy(pdist.particles.v[1, :])
    @test minimum(sdist.spline(qα) for qα in q) < 0
    for (_, f) in COEFFICIENTS
        @test_throws DomainError f(sdist, pdist, q)
    end
end
