using LinearAlgebra: I, Symmetric, cholesky, dot, istriu, norm, opnorm, qr, svdvals
using PoissonSolvers: Potential
using Random: Xoshiro
using SimpleSplines: BSplineBasis, PeriodicBSplineBasis, UniformMesh, Free, Dirichlet, (..),
                     SplineQuadrature, mass_matrix, GeneralMesh, domain, leftendpoint,
                     rightendpoint
using Test
using VlasovMethods
using VlasovMethods: SplineLowRankFields1D1V, reflect_v, inner_x, inner_v, orthonormalize_x,
                     orthonormalize_v, _kfield!, _lfield!, _sfield!, _potential!,
                     _spline_1d1v_rhs!, _averaged_form, _spline_mass

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

const VM = VlasovMethods

# Every expected value comes from the test's own dense tables (`spline_tables` of the helper,
# with dense, symmetrised mass matrices and their Cholesky factors), from S1a's full-order field
# `_spline_1d1v_rhs!` applied to the coefficient matrix that the test forms from the factors, or
# from the test's own dense quadrature of the cyclically averaged form, which the test checks
# against S1a's `_averaged_form`. Seeded draws are standard normal in Float64, rounded to T.

draw(rng, T, dims...) = T.(randn(rng, dims...))
rel(a, b) = norm(a .- b) / norm(b)

# a draw of its own seed: the factors of seed s are X of s, V of s + 10, K of s + 20, L of s + 30,
# S of s + 40 and a frozen ϕ of s + 50, the draws of the measurements of the spec
sdraw(T, dims...; seed) = draw(Xoshiro(seed), T, dims...)

# the periodic cubic x-basis of [0, 4π) with nxc cells, and the clamped v-basis of degree p of
# the interval v with nvc cells
function bases(T, nxc, nvc; v = (-25, 25), p = 3)
    xb = PeriodicBSplineBasis(UniformMesh(nxc, zero(T) .. T(4π)), 3)
    vb = BSplineBasis(UniformMesh(nvc, T(v[1]) .. T(v[2])), p, Free())
    return xb, vb
end

# the space, the workspace, the dense tables and S1a's model of the cubic bases
function make_problem(T, nxc, nvc)
    xb, vb = bases(T, nxc, nvc)
    space = SplineFactorSpace(xb, vb)
    nx, nv = size(space)
    model = VlasovPoisson(TensorSplineDistribution(xb, vb, zeros(T, nx, nv)), Potential(xb))
    return (; space, w = SplineLowRankFields1D1V(space), t = spline_tables(xb, vb), model,
        nx, nv)
end
const PROBLEMS = Dict{Any, Any}()
problem(T, nxc, nvc) = get!(() -> make_problem(T, nxc, nvc), PROBLEMS, (T, nxc, nvc))

# S1a's full-order field and potential of the coefficient matrix F
fdot(P, F) = _spline_1d1v_rhs!(similar(F), P.model, F)
phi(P, F) = VM._spline_potential(P.model, F)

# B[i, j] = T̄(φᵢψⱼ, h, f) for h = ϕ ⊗ 1 + 1 ⊗ v²/2 and the f of F, the three cyclic terms of
#   T(Φ, h, f) = ∫ Φ (hₓ f_v − h_v fₓ),  T(h, f, Φ) = ∫ h (fₓ Φ_v − f_v Φₓ),
#   T(f, Φ, h) = ∫ f (Φₓ h_v − Φ_v hₓ),
# each formed on its own on the tensor quadrature
function dense_bracket(t, F, ϕ)
    f0, fx, fv = grid_values(t, F)
    h0 = t.X0' * ϕ .+ (t.vq .^ 2 ./ 2)'
    hx = t.X1' * ϕ
    hv = permutedims(t.vq)
    W = t.W
    t1 = t.X0 * (W .* (hx .* fv .- hv .* fx)) * t.V0'
    t2 = t.X0 * (W .* h0 .* fx) * t.V1' .- t.X1 * (W .* h0 .* fv) * t.V0'
    t3 = t.X1 * (W .* f0 .* hv) * t.V0' .- t.X0 * (W .* f0 .* hx) * t.V1'
    return (t1 .+ t2 .+ t3) ./ 3
end
dense_field(t, F, ϕ) = minv(t, dense_bracket(t, F, ϕ))

# the coefficient vector of g(v) in the v-basis, and of h = ϕ ⊗ 1 + 1 ⊗ v²/2
vcoef(t, g) = t.Cv \ (t.V0 * (t.wv .* g.(t.vq)))
function hcoef(t, ϕ)
    return ϕ * (t.Cv \ t.iv)' .+ (t.Cx \ t.ix) * vcoef(t, v -> v^2 / 2)'
end

# the columns of A orthonormalised in the mass matrix of the Cholesky factorisation C
orth(C, A) = C.U \ Matrix(qr(C.U * A).Q)

# the factors of the seed: X, V orthonormal, S, K, L standard normal
function factors(P, T, r, seed)
    X = orth(P.t.Cx, sdraw(T, P.nx, r; seed))
    V = orth(P.t.Cv, sdraw(T, P.nv, r; seed = seed + 10))
    return (X = X, V = V, S = sdraw(T, r, r; seed = seed + 40),
        K = sdraw(T, P.nx, r; seed = seed + 20), L = sdraw(T, P.nv, r; seed = seed + 30))
end

# the projections of the full-order field of S1a and of the test's dense T̄
kref(P, K, V, Fd) = Fd(K * V') * P.t.Mv * V
lref(P, X, L, Fd) = Fd(X * L')' * P.t.Mx * X
sref(P, X, S, V, Fd) = X' * P.t.Mx * Fd(X * S * V') * P.t.Mv * V

kfield(P, K, V) = _kfield!(similar(K), K, V, P.w)
lfield(P, X, L) = _lfield!(similar(L), X, L, P.w)
sfield(P, X, S, V) = _sfield!(similar(S), X, S, V, P.w)

maxdev(A) = maximum(abs, A - I)
morth(M, Q) = maxdev(Q' * M * Q)

# the fixture coefficient matrix Σₖ 2⁻ᵏ aₖ bₖᵀ of seed 3
function fixture(T, nx, nv)
    rng = Xoshiro(3)
    F = zeros(nx, nv)
    for k in 1:5
        F .+= 2.0^(-k) .* randn(rng, nx) .* randn(rng, nv)'
    end
    return T.(F)
end

const GRIDS = ((8, 10, 3), (16, 97, 4), (12, 98, 6), (17, 37, 1))

@testset "the guard of the space and its sizes, $T" for T in (Float64, Float32)
    for (nxc, nvc) in ((4, 1), (16, 97), (16, 98)), p in (2, 3)

        vint = (nxc, nvc) == (16, 97) ? (-6, 5.9) : (-25, 25)
        xb, vb = bases(T, nxc, nvc; v = vint, p)
        space = SplineFactorSpace(xb, vb)
        @test size(space) == (nxc, nvc + p)
        @test eltype(space) === T
        @test space.xbasis === xb && space.vbasis === vb
        Mx, Mv = mass_matrix(SplineQuadrature(xb)), mass_matrix(SplineQuadrature(vb))
        @test space.Mx == Mx && space.Mv == Mv
        @test eltype(space.Mx) === eltype(space.Ux) === eltype(space.Uv) === T
        for (M, U) in ((Mx, space.Ux), (Mv, space.Uv))
            Ud = Matrix(U)
            @test istriu(Ud)
            @test maximum(abs, Ud' * Ud - M) ≤ 64 * eps(T) * maximum(abs, M)
        end
    end

    xb, vb = bases(T, 16, 97)
    dist = TensorSplineDistribution(xb, vb, zeros(T, 16, 100))
    space = SplineFactorSpace(dist)
    @test space.xbasis === dist.xbasis && space.vbasis === dist.vbasis
    @test size(space) == (16, 100)

    mesh = UniformMesh(97, T(-25) .. T(25))
    v² = r"^ArgumentError: the v-basis of a SplineFactorSpace reproduces v²"
    @test_throws v² SplineFactorSpace(xb, BSplineBasis(mesh, 1, Free()))
    @test_throws v² SplineFactorSpace(xb, BSplineBasis(mesh, 3, Dirichlet()))
    @test_throws v² SplineFactorSpace(xb, PeriodicBSplineBasis(mesh, 3))
    @test_throws r"^ArgumentError: the x-basis of a SplineFactorSpace is periodic" SplineFactorSpace(
        BSplineBasis(UniformMesh(16, zero(T) .. T(4π)), 3, Free()), vb)
    @test SplineFactorSpace(xb, BSplineBasis(mesh, 2, Free())) isa SplineFactorSpace{T}
end

@testset "the guard of the space on two element types" begin
    xb, _ = bases(Float32, 16, 97)
    _, vb = bases(Float64, 16, 97)
    @test_throws r"^ArgumentError: the bases of a SplineFactorSpace have one element type" SplineFactorSpace(
        xb, vb)
end

@testset "inner products and orthonormal factors, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, t, nx, nv) = P
    rng = Xoshiro(11)
    for (n, M, inner, orthonormalize, other) in ((nx, t.Mx, inner_x, orthonormalize_x, nv),
        (nv, t.Mv, inner_v, orthonormalize_v, nx))
        A, B = draw(rng, T, n, 5), draw(rng, T, n, 5)
        a, b = draw(rng, T, n), draw(rng, T, n)
        G = inner(space, A, B)
        @test eltype(G) === T
        @test rel(G, A' * M * B) ≤ 8 * eps(T)
        g = inner(space, a, b)
        @test g isa T
        @test abs(g - a' * M * b) ≤ 8 * eps(T) * norm(a) * norm(M * b)
        rows = r"^DimensionMismatch: a factor of \d+ rows on the [xv] basis"
        @test_throws rows inner(space, draw(rng, T, other, 5), B)
        @test_throws rows inner(space, A, draw(rng, T, other, 5))
        @test_throws rows orthonormalize(space, draw(rng, T, other, 3))

        Q, R = orthonormalize(space, A)
        @test eltype(Q) === eltype(R) === T
        @test size(Q) == (n, 5) && size(R) == (5, 5)
        @test istriu(R)
        @test morth(M, Q) ≤ 50 * eps(T)
        @test maximum(abs, Q * R - A) ≤ 50 * eps(T) * maximum(abs, A)

        Q, R = orthonormalize(space, zeros(T, n, 3))
        @test morth(M, Q) ≤ 50 * eps(T)
        @test all(iszero, R)

        @test_throws r"^ArgumentError: cannot orthonormalise" orthonormalize(
            space, draw(rng, T, n, n + 1))
        A[2, 3] = T(NaN)
        Q, R = orthonormalize(space, A)
        @test any(isnan, R)
    end
end

@testset "the reflection and parity, $T" for T in (Float64, Float32)
    for nvc in (97, 98)
        P = problem(T, 16, nvc)
        (; space, t, nv) = P
        # the nodes of the v quadrature are symmetric, so the reversal of the nodes is v ↦ -v
        @test maximum(abs, t.vq .+ reverse(t.vq)) ≤ 256 * eps(T) * maximum(abs, t.vq)
        mirror = t.V0[end:-1:1, end:-1:1] .- t.V0
        @test maximum(abs, mirror) ≤ 256 * eps(T) * maximum(abs, t.V0)
        mirror′ = t.V1[end:-1:1, end:-1:1] .+ t.V1
        @test maximum(abs, mirror′) ≤ 256 * eps(T) * maximum(abs, t.V1)
        R = Matrix{T}(I, nv, nv)[end:-1:1, :]
        Mv = Matrix(space.Mv)
        @test maximum(abs, Mv * R - R * Mv) ≤ 256 * eps(T) * maximum(abs, Mv)

        rng = Xoshiro(5)
        B = draw(rng, T, nv, 3)
        @test reflect_v(space, B) == B[nv:-1:1, :]
        @test reflect_v(space, B[:, 1]) == B[nv:-1:1, 1]
        @test eltype(reflect_v(space, B)) === T
        rows = r"^DimensionMismatch: a factor of \d+ rows on the v basis"
        @test_throws rows reflect_v(space, draw(rng, T, nv - 1, 3))
        @test_throws rows reflect_v(space, draw(rng, T, nv + 1))
        for s in (1, -1)
            Q, _ = orthonormalize_v(space, (B .+ s .* reflect_v(space, B)) ./ 2)
            @test norm(Q .- s .* reflect_v(space, Q)) ≤ 64 * eps(T) * norm(Q)
        end
    end

    xb, vb = bases(T, 16, 97; v = (-6, 5.9))
    space = SplineFactorSpace(xb, vb)
    @test_throws r"^ArgumentError: .*symmetric" reflect_v(space, zeros(T, 100, 2))

    # a mesh of [-25, 25] that is symmetric but not uniform: its ends pass the comparison, and the
    # mesh type is refused
    half = T(25) .* (range(zero(T), one(T); length = 49)[2:end]) .^ T(1.5)
    vb = BSplineBasis(GeneralMesh(vcat(-reverse(half), zero(T), half)), 3, Free())
    space = SplineFactorSpace(xb, vb)
    @test leftendpoint(domain(vb)) == -rightendpoint(domain(vb))
    @test_throws r"^ArgumentError: .*symmetric" reflect_v(space, zeros(T, size(space)[2], 2))
end

@testset "factorisation: Eckart–Young in the mass norm, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, t, nx, nv) = P
    F = fixture(T, nx, nv)
    Ux, Uv = t.Cx.U, t.Cv.U
    σ = svdvals(Ux * F * Uv')
    scale = norm(Ux * F * Uv')
    dist = TensorSplineDistribution(space.xbasis, space.vbasis, F)
    for r in (3, 4, 6)
        d = LowRankDistribution(space, F, r)
        @test eltype(d.X) === eltype(d.S) === eltype(d.V) === T
        @test d.nfixed == 0
        @test d.space === space
        d₂ = LowRankDistribution(dist, r)
        @test d₂.X == d.X && d₂.S == d.S && d₂.V == d.V
        @test morth(t.Mx, d.X) ≤ 50 * eps(T)
        @test morth(t.Mv, d.V) ≤ 50 * eps(T)
        @test istriu(d.S)
        err = norm(Ux * (F - d.X * d.S * d.V') * Uv')
        @test abs(err - sqrt(sum(abs2, σ[(r + 1):end]))) ≤ 50 * eps(T) * scale
    end
end

@testset "factorisation: the fixed modes 1, v and v², $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, t, nx, nv) = P
    F = fixture(T, nx, nv)
    c = [vcoef(t, v -> one(v)) vcoef(t, v -> v) vcoef(t, v -> v^2)]
    for r in (3, 4, 6)
        d = LowRankDistribution(space, F, r; fixed = c)
        @test d.nfixed == 3
        for k in 1:3
            Φ, Vk = c[:, 1:k], d.V[:, 1:k]
            @test norm(Φ - Vk * (Vk' * t.Mv * Φ)) ≤ 200 * eps(T) * norm(Φ)
        end
        G = d.X * d.S * d.V'
        for k in 1:3
            m = F * t.Mv * c[:, k]
            @test maximum(abs, G * t.Mv * c[:, k] - m) ≤ 200 * eps(T) * maximum(abs, m)
        end
    end
    d = LowRankDistribution(space, F, 3)
    m = F * t.Mv * c[:, 1]
    @test maximum(abs, d.X * d.S * d.V' * t.Mv * c[:, 1] - m) >
          2000 * eps(T) * maximum(abs, m)
end

@testset "factorisation: full rank, zero and the refusals, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, t, nx, nv) = P
    F = fixture(T, nx, nv)
    c = [vcoef(t, v -> one(v)) vcoef(t, v -> v) vcoef(t, v -> v^2)]

    d = LowRankDistribution(space, F, nx)
    @test maximum(abs, d.X * d.S * d.V' - F) ≤ 100 * eps(T) * maximum(abs, F)
    @test TensorSplineDistribution(d) isa TensorSplineDistribution{T}
    @test TensorSplineDistribution(d).coefficients ≈ d.X * d.S * d.V'
    @test TensorSplineDistribution(d).xbasis === space.xbasis

    for fixed in (zeros(T, nv, 0), c)
        d = LowRankDistribution(space, zeros(T, nx, nv), 4; fixed)
        @test all(iszero, d.S)
        @test all(iszero, TensorSplineDistribution(d).coefficients)
        @test morth(t.Mx, d.X) ≤ 50 * eps(T)
        @test morth(t.Mv, d.V) ≤ 50 * eps(T)
    end

    for r in (0, min(nx, nv) + 1)
        @test_throws r"^ArgumentError: the rank" LowRankDistribution(space, F, r)
    end
    @test_throws r"^ArgumentError: 4 fixed modes exceed the rank 3" LowRankDistribution(
        space, F, 3; fixed = [c c[:, 1]])
    @test_throws r"^ArgumentError: the fixed modes are linearly dependent" LowRankDistribution(
        space, F, 4; fixed = [c[:, 1] c[:, 2] 2c[:, 2]])
    @test_throws r"^ArgumentError: the fixed modes are linearly dependent" LowRankDistribution(
        space, F, 4; fixed = zeros(T, nv, 1))
    for bad in (T(NaN), T(Inf))
        G = copy(F)
        G[3, 7] = bad
        @test_throws r"^ArgumentError: the coefficient matrix holds a non-finite" LowRankDistribution(
            space, G, 4)
        e = copy(c)
        e[5, 2] = bad
        @test_throws r"^ArgumentError: the fixed modes hold a non-finite" LowRankDistribution(
            space, F, 4; fixed = e)
    end
    @test_throws r"^DimensionMismatch: fixed modes of" LowRankDistribution(
        space, F, 4; fixed = c[1:(nv - 1), :])
    @test_throws r"^DimensionMismatch: a coefficient matrix of size" LowRankDistribution(
        space, F[:, 1:(nv - 1)], 4)
end

@testset "the dense references, $T" for T in (Float64, Float32)
    # the test's dense T̄ equals S1a's averaged form, and its field S1a's field
    P = problem(T, 16, 97)
    (; t, nx, nv, model) = P
    rng = Xoshiro(21)
    F, a = draw(rng, T, nx, nv), draw(rng, T, nx, nv)
    ϕ = phi(P, F)
    B = dense_bracket(t, F, ϕ)
    form = _averaged_form(model.distribution, a, hcoef(t, ϕ), F)
    @test abs(sum(a .* B) - form) ≤ 1024 * eps(T) * norm(a) * norm(B)
    @test rel(dense_field(t, F, ϕ), fdot(P, F)) ≤ 1024 * eps(T)
end

@testset "the factored fields agree with the full-order field, $T" for T in (Float64, Float32)
    tol = 16384 * eps(T)
    for (nxc, nvc, r) in GRIDS, seed in 1:5

        P = problem(T, nxc, nvc)
        (; t) = P
        (; X, V, S, K, L) = factors(P, T, r, seed)
        Fd1 = F -> fdot(P, F)
        Fd2 = F -> dense_field(t, F, phi(P, F))
        k̇, l̇ = kfield(P, K, V), lfield(P, X, L)
        @test eltype(k̇) === eltype(l̇) === T
        for Fd in (Fd1, Fd2)
            @test rel(k̇, kref(P, K, V, Fd)) ≤ tol
            @test rel(l̇, lref(P, X, L, Fd)) ≤ tol
        end
        if r > 1
            ṡ = sfield(P, X, S, V)
            @test eltype(ṡ) === T
            for Fd in (Fd1, Fd2)
                @test rel(ṡ, sref(P, X, S, V, Fd)) ≤ tol
            end
        end
    end
end

@testset "the smallest grid and the scale extremes, $T" for T in (Float64, Float32)
    tol = 16384 * eps(T)
    Fd(P) = F -> fdot(P, F)
    P = problem(T, 4, 1)
    @test size(P.space) == (4, 4)
    for seed in 1:5
        (; X, V, S, K, L) = factors(P, T, 2, seed)
        fields = (kfield(P, K, V), lfield(P, X, L), sfield(P, X, S, V))
        @test all(f -> all(isfinite, f), fields)
        @test rel(fields[1], kref(P, K, V, Fd(P))) ≤ tol
        @test rel(fields[2], lref(P, X, L, Fd(P))) ≤ tol
        @test rel(fields[3], sref(P, X, S, V, Fd(P))) ≤ tol
    end
    P = problem(T, 16, 97)
    for s in (T(1e-8), T(1e8)), seed in 1:5

        (; V, K) = factors(P, T, 4, seed)
        k̇ = kfield(P, s .* K, V)
        @test all(isfinite, k̇)
        @test rel(k̇, kref(P, s .* K, V, Fd(P))) ≤ tol
    end
end

@testset "the S-field vanishes at rank one, $T" for T in (Float64, Float32)
    P = problem(T, 17, 37)
    for seed in 1:5
        (; X, V, S) = factors(P, T, 1, seed)
        ṡ = sfield(P, X, S, V)
        F = X * S * V'
        B = P.t.Mx * fdot(P, F) * P.t.Mv
        @test norm(ṡ) ≤ 4096 * eps(T) * norm(X) * norm(B) * norm(V)
    end
end

@testset "the potential and a frozen potential, $T" for T in (Float64, Float32)
    tol = 4096 * eps(T)
    for (nxc, nvc, r) in GRIDS, seed in 1:5

        P = problem(T, nxc, nvc)
        (; t, w, nx) = P
        (; X, V, S, K, L) = factors(P, T, r, seed)
        ϕK = _potential!(zeros(T, nx), K, V, w)
        ϕS = _potential!(zeros(T, nx), X, S, V, w)
        @test eltype(ϕK) === eltype(ϕS) === T
        @test rel(ϕK, phi(P, K * V')) ≤ tol
        @test rel(ϕS, phi(P, X * S * V')) ≤ tol

        ϕ = sdraw(T, nx; seed = seed + 50)
        Fd = F -> dense_field(t, F, ϕ)
        @test rel(_kfield!(similar(K), K, V, w, ϕ), kref(P, K, V, Fd)) ≤ tol
        @test rel(_lfield!(similar(L), X, L, w, ϕ), lref(P, X, L, Fd)) ≤ tol
        r > 1 && @test rel(_sfield!(similar(S), X, S, V, w, ϕ), sref(P, X, S, V, Fd)) ≤ tol

        ϕL = _potential!(zeros(T, nx), X, L, w)
        @test rel(kfield(P, K, V), _kfield!(similar(K), K, V, w, ϕK)) ≤ 8 * eps(T)
        @test rel(lfield(P, X, L), _lfield!(similar(L), X, L, w, ϕL)) ≤ 8 * eps(T)
        if r > 1
            @test rel(sfield(P, X, S, V), _sfield!(similar(S), X, S, V, w, ϕS)) ≤ 8 * eps(T)
        end
    end
end

@testset "the potential and the K-field on a nonuniform x-mesh, $T" for T in (Float64, Float32)
    # the 16 cells of [0, 4π) of S1a's test, of widths between 0.52 and 1.48 times the mean, where
    # the mean shifts of the load and of ϕ change the solve by more than round-off
    L = T(4π)
    bp = collect(range(zero(T), L; length = 17))
    bp[2:16] .+= (L / 16) .* T(1.25) .* sin.(T(2π) .* (1:15) ./ 16)
    xb = PeriodicBSplineBasis(GeneralMesh(bp), 3)
    _, vb = bases(T, 16, 97)
    space = SplineFactorSpace(xb, vb)
    nx, nv = size(space)
    model = VlasovPoisson(TensorSplineDistribution(xb, vb, zeros(T, nx, nv)), Potential(xb))
    P = (; space, w = SplineLowRankFields1D1V(space),
        t = spline_tables(xb, vb), model, nx, nv)
    for seed in 1:5
        (; X, V, S, K) = factors(P, T, 4, seed)
        @test rel(_potential!(zeros(T, nx), K, V, P.w), phi(P, K * V')) ≤ 4096 * eps(T)
        @test rel(_potential!(zeros(T, nx), X, S, V, P.w), phi(P, X * S * V')) ≤
              4096 * eps(T)
        @test rel(kfield(P, K, V), kref(P, K, V, F -> fdot(P, F))) ≤ 16384 * eps(T)
    end
end

# the DLR field K̇ Vᵀ + X L̇ᵀ − X Ṡ Vᵀ at F = X S Vᵀ
function dlr_field(P, X, S, V)
    K, L = X * S, V * S'
    return kfield(P, K, V) * V' + X * lfield(P, X, L)' - X * sfield(P, X, S, V) * V'
end

@testset "the full-rank DLR field is the full-order field, $T" for T in (Float64, Float32)
    for (nxc, nvc) in ((8, 10), (8, 5), (6, 9)), seed in 1:5

        P = problem(T, nxc, nvc)
        (; X, V, S) = factors(P, T, P.nx, seed)
        F = X * S * V'
        @test rel(dlr_field(P, X, S, V), fdot(P, F)) ≤ 4096 * eps(T)
    end
end

@testset "the energy rate and the mass rate, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; t, nx, nv, model) = P
    minner(A, B) = dot(A, t.Mx * B * t.Mv)
    one_xv = (t.Cx \ t.ix) * (t.Cv \ t.iv)'
    rates = T[]
    for seed in 1:5, with in (true, false)

        X = orth(t.Cx, sdraw(T, nx, 4; seed))
        Ṽ = sdraw(T, nv, 4; seed = seed + 10)
        if with
            Ṽ[:, 1], Ṽ[:, 2] = vcoef(t, v -> one(v)), vcoef(t, v -> v^2)
        end
        V = orth(t.Cv, Ṽ)
        S = sdraw(T, 4, 4; seed = seed + 40)
        F = X * S * V'
        PF = dlr_field(P, X, S, V)
        H = hcoef(t, phi(P, F))
        rate = abs(minner(H, PF)) / sqrt(minner(H, H) * minner(PF, PF))
        if with
            @test rate ≤ 4096 * eps(T)
            Ḟ = fdot(P, F)
            @test abs(minner(one_xv, PF) - _spline_mass(model, Ḟ)) ≤
                  4096 * eps(T) * sqrt(minner(one_xv, one_xv) * minner(Ḟ, Ḟ))
        else
            push!(rates, rate)
        end
    end
    @test maximum(rates) ≥ 1e-3
end

@testset "the grading of the fields, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, t, nx, nv) = P
    for seed in 1:5
        C = sdraw(T, nx, nv; seed = seed + 1)
        for s in (1, -1)
            Ḟ = fdot(P, (C .+ s .* reflect(C)) ./ 2)
            @test norm(Ḟ .+ s .* reflect(Ḟ)) / 2 ≤ 1024 * eps(T) * norm(Ḟ)
        end
        X = orth(t.Cx, sdraw(T, nx, 3; seed = seed + 80))
        B = sdraw(T, nv, 3; seed = seed + 70)
        S = sdraw(T, 3, 3; seed = seed + 90)
        for s in (1, -1)
            V = orth(t.Cv, (B .+ s .* reflect_v(space, B)) ./ 2)
            K, L = X * S, V * S'
            ḞK = fdot(P, K * V')
            ḞS = fdot(P, X * S * V')
            @test norm(kfield(P, K, V)) ≤ 1024 * eps(T) * norm(ḞK * t.Mv) * opnorm(V)
            @test norm(sfield(P, X, S, V)) ≤
                  1024 * eps(T) * norm(X' * t.Mx * ḞS * t.Mv) * opnorm(V)
            l̇ = lfield(P, X, L)
            @test norm(l̇ .+ s .* reflect_v(space, l̇)) ≤ 4096 * eps(T) * norm(l̇)
        end
    end
end

# The bytes of one call of `f(args...)` after a warm call, the minimum over three calls, through
# a barrier with concrete argument types.
function allocations(f::F, args::A) where {F, A}
    f(args...)
    return minimum(@allocated(f(args...)) for _ in 1:3)
end

# the control: one nx × nv array, (K Vᵀ) V
control(K, V) = (K * V') * V

@testset "no field forms an nx × nv array, $T" for T in (Float64, Float32)
    # `--check-bounds=yes` changes the compiled code that the bound measures, so the measurement
    # is skipped where a run forces `--check-bounds=yes`, and under coverage instrumentation, as
    # the allocation tests of the other test files are
    if Base.JLOptions().check_bounds != 1 && Base.JLOptions().code_coverage == 0
        for (nxc, nvc, r) in ((64, 125, 4), (64, 125, 8), (512, 4093, 2))
            xb, vb = bases(T, nxc, nvc)
            space = SplineFactorSpace(xb, vb)
            w = SplineLowRankFields1D1V(space)
            nx, nv = size(space)
            rng = Xoshiro(7)
            X, V, S = draw(rng, T, nx, r), draw(rng, T, nv, r), draw(rng, T, r, r)
            K, L, ϕ = draw(rng, T, nx, r), draw(rng, T, nv, r), draw(rng, T, nx)
            bound = 160 * (nx + nv) * r * sizeof(T) + 65536

            @test allocations(_kfield!, (similar(K), K, V, w)) ≤ bound
            @test allocations(_lfield!, (similar(L), X, L, w)) ≤ bound
            @test allocations(_sfield!, (similar(S), X, S, V, w)) ≤ bound
            @test allocations(_kfield!, (similar(K), K, V, w, ϕ)) ≤ bound
            @test allocations(_lfield!, (similar(L), X, L, w, ϕ)) ≤ bound
            @test allocations(_sfield!, (similar(S), X, S, V, w, ϕ)) ≤ bound
            @test allocations(_potential!, (similar(ϕ), K, V, w)) ≤ bound
            @test allocations(_potential!, (similar(ϕ), X, S, V, w)) ≤ bound

            if nv == 4096
                @test bound < nx * nv * sizeof(T)
                @test allocations(control, (K, V)) > bound
            end
        end
    end
end

@testset "the edges of the workspace and the fields, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; w, nx, nv) = P
    r = 4
    (; X, V, S, K, L) = factors(P, T, r, 1)
    sentinel = T(7)
    out(dims...) = fill(sentinel, dims...)
    kept(A) = all(==(sentinel), A)

    dm(name) = Regex("^DimensionMismatch: " * name)
    o = out(nx + 1, r)
    @test_throws dm("_kfield!") _kfield!(o, draw(Xoshiro(1), T, nx + 1, r), V, w)
    @test kept(o)
    o = out(nv - 1, r)
    @test_throws dm("_lfield!") _lfield!(o, X, draw(Xoshiro(1), T, nv - 1, r), w)
    @test kept(o)
    o = out(r, r + 1)
    @test_throws dm("_sfield!") _sfield!(o, X, draw(Xoshiro(1), T, r, r + 1), V, w)
    @test kept(o)
    o = out(nv, r)
    @test_throws dm("_lfield!") _lfield!(o, X, L, w, draw(Xoshiro(1), T, nx + 1))
    @test kept(o)
    o = out(nx - 1)
    @test_throws dm("_potential!") _potential!(o, K, V, w)
    @test kept(o)
    @test_throws dm("_potential!") _potential!(o, X, S, V, w)
    @test kept(o)

    K₁ = copy(K)
    @test_throws r"^ArgumentError: _kfield!: the output K̇ shares memory" _kfield!(K₁, K₁, V, w)
    L₁, S₁ = copy(L), copy(S)
    @test_throws r"^ArgumentError: _lfield!: the output L̇ shares memory" _lfield!(L₁, X, L₁, w)
    @test_throws r"^ArgumentError: _sfield!: the output Ṡ shares memory" _sfield!(
        S₁, X, S₁, V, w)
    # a vector that shares the memory of a factor of one column
    K₂ = K[:, 1:1]
    @test_throws r"^ArgumentError: _potential!: the output ϕ shares memory" _potential!(
        vec(K₂), K₂, V[:, 1:1], w)
    @test_throws r"^ArgumentError: _potential!: the output ϕ shares memory" _potential!(
        vec(K₂), K₂, S[1:1, 1:1], V[:, 1:1], w)

    X0, V0, S0 = zeros(T, nx, 0), zeros(T, nv, 0), zeros(T, 0, 0)
    @test size(_kfield!(zeros(T, nx, 0), X0, V0, w)) == (nx, 0)
    @test size(_lfield!(zeros(T, nv, 0), X0, V0, w)) == (nv, 0)
    @test size(_sfield!(zeros(T, 0, 0), X0, S0, V0, w)) == (0, 0)
    @test _potential!(out(nx), X0, V0, w) == zeros(T, nx)
    @test _potential!(out(nx), X0, S0, V0, w) == zeros(T, nx)

    Z = (X = zero(X), V = zero(V), S = zero(S), K = zero(K), L = zero(L))
    @test all(==(0), kfield(P, Z.K, Z.V))
    @test all(==(0), lfield(P, Z.X, Z.L))
    @test all(==(0), sfield(P, Z.X, Z.S, Z.V))
    @test all(==(0), _potential!(out(nx), Z.K, Z.V, w))
    @test all(==(0), _potential!(out(nx), Z.X, Z.S, Z.V, w))

    for bad in (T(NaN), T(Inf))
        Kb, Lb, Sb = copy(K), copy(L), copy(S)
        Kb[2, 3] = Lb[5, 2] = Sb[1, 2] = bad
        @test any(!isfinite, kfield(P, Kb, V))
        @test any(!isfinite, lfield(P, X, Lb))
        @test any(!isfinite, sfield(P, X, Sb, V))
        @test any(!isfinite, _potential!(out(nx), Kb, V, w))
    end

    if T === Float32
        # the refusal of another element type is the absence of a method, which has no guard of
        # its own to write a message
        # ast-grep-ignore: test-throws-bare-type
        @test_throws MethodError _kfield!(
            zeros(Float64, nx, r), Float64.(K), Float64.(V), w)
        xb, vb = bases(T, 768, 97)
        @test_throws r"^ArgumentError: the matrix is singular beyond the constants" SplineLowRankFields1D1V(
            SplineFactorSpace(xb, vb))
    end
end

@testset "the workspace of a distribution, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, nx, nv) = P
    d = LowRankDistribution(space, fixture(T, nx, nv), 4)
    @test d.space === space
    @test size(d.X) == (16, 4) && size(d.V) == (100, 4)
    w₁, w₂ = SplineLowRankFields1D1V(d), SplineLowRankFields1D1V(space)
    K, L = d.X * d.S, d.V * d.S'
    @test _kfield!(similar(K), K, d.V, w₁) == _kfield!(similar(K), K, d.V, w₂)
    @test _lfield!(similar(L), d.X, L, w₁) == _lfield!(similar(L), d.X, L, w₂)
    @test _sfield!(similar(d.S), d.X, d.S, d.V, w₁) ==
          _sfield!(similar(d.S), d.X, d.S, d.V, w₂)
    d₂ = LowRankDistribution(space, d.X, d.S, d.V; nfixed = 1)
    @test d₂.space === space && d₂.nfixed == 1
end

# one step of the classical RK4 tableau of the autonomous field g ↦ rhs(g)
function rk4(rhs, y, τ)
    k1 = rhs(y)
    k2 = rhs(y .+ (τ / 2) .* k1)
    k3 = rhs(y .+ (τ / 2) .* k2)
    k4 = rhs(y .+ τ .* k3)
    return y .+ (τ / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

# the full-order steps of the two compositions at full rank
fullorder(::Val{:lie}, P, F, τ) = rk4(g -> fdot(P, g), F, τ)
function fullorder(::Val{:strang}, P, F, τ)
    ϕ0 = phi(P, F)
    ϕh = phi(P, rk4(g -> dense_field(P.t, g, ϕ0), F, τ / 2))
    return rk4(g -> dense_field(P.t, g, ϕh), F, τ)
end

# the Landau state c_x c_vᵀ, the L² projections of 1 + cos(x/2)/2 and of the Maxwellian
function landau(P)
    t = P.t
    T = eltype(t.wx)
    cx = t.Cx \ (t.X0 * (t.wx .* (1 .+ cos.(t.xq ./ 2) ./ 2)))
    cv = t.Cv \ (t.V0 * (t.wv .* exp.(-t.vq .^ 2 ./ 2) ./ sqrt(T(2π))))
    return cx * cv'
end

@testset "the projector splitting steps a spline state, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, w, nx) = P
    F₀ = landau(P)
    for c in (:lie, :strang), τ in (T(1e-4), T(-1e-4))

        Y = LowRankDistribution(space, F₀, nx)
        F = F₀
        worst = zero(T)
        for n in 0:19
            @test lowrank_step!(Y, ProjectorSplitting(c), w, n * τ, τ) === Y
            F = fullorder(Val(c), P, F, τ)
            worst = max(worst, rel(Y.X * Y.S * Y.V', F))
        end
        @test worst ≤ 1024 * eps(T)
        @test size(Y.X) == (16, 16) && size(Y.S) == (16, 16) && size(Y.V) == (100, 16)
        @test eltype(Y.X) === eltype(Y.S) === eltype(Y.V) === T
        @test Y.nfixed == 0
    end
end

@testset "a rank-4 step, and the time that the fields ignore, $T" for T in (Float64, Float32)
    P = problem(T, 16, 97)
    (; space, w, t) = P
    F₀ = landau(P)
    for c in (:lie, :strang)
        Y₀ = LowRankDistribution(space, F₀, 4)
        Y₇ = LowRankDistribution(space, copy(Y₀.X), copy(Y₀.S), copy(Y₀.V))
        lowrank_step!(Y₀, ProjectorSplitting(c), w, 0, T(1e-4))
        lowrank_step!(Y₇, ProjectorSplitting(c), w, 7, T(1e-4))
        @test all(isfinite, Y₀.X) && all(isfinite, Y₀.S) && all(isfinite, Y₀.V)
        @test morth(t.Mx, Y₀.X) ≤ 50 * eps(T)
        @test morth(t.Mv, Y₀.V) ≤ 50 * eps(T)
        @test Y₀.nfixed == 0
        @test Y₀.X == Y₇.X && Y₀.S == Y₇.S && Y₀.V == Y₇.V
    end
end
