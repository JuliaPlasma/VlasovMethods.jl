using ForwardDiff
using GeometricMachineLearning
using LinearAlgebra
using SimpleSplines
using StableRNGs
using Test
using VlasovMethods

const REAL_ELTYPES = (Float32, Float64)

maxwellian(v) = exp(-v^2 / 2) / sqrt(2 * oftype(v, π))
dmaxwellian(v) = -v * maxwellian(v)
two_stream(v) = v^2 * maxwellian(v)
dtwo_stream(v) = (2v - v^3) * maxwellian(v)

# The orthogonal parity split of `BilinearStiefelLayer`: index `k` with its mirror `n + 1 - k`,
# the even combinations on top of the odd ones.
function parity_split(::Type{T}, n) where {T}
    h = n ÷ 2
    S = zeros(T, n, n)
    for k in 1:h
        S[k, k] = S[k, n + 1 - k] = S[h + k, k] = 1 / sqrt(T(2))
        S[h + k, n + 1 - k] = -1 / sqrt(T(2))
    end
    return S
end

# The value of every basis function at every point, one row per point.
evaluation_matrix(b, xs, d = 0) = reduce(hcat, (evaluate(b, j, xs, d) for j in 1:nbasis(b)))

# The coefficients of the spline that interpolates `f` at the Greville abscissae of `b`.
interpolant(b, f) = evaluation_matrix(b, nodes(b)) \ f.(nodes(b))

# The least-squares slope of log(err) against log(1/N): the order of convergence.
function fitted_order(Ns, errs)
    lx = log.(1 ./ Ns)
    ly = log.(errs)
    return sum((lx .- sum(lx) / length(lx)) .* (ly .- sum(ly) / length(ly))) /
           sum(abs2, lx .- sum(lx) / length(lx))
end

stiefel(rng, ::Type{T}, n, m) where {T} = StiefelManifold(Matrix(qr(randn(rng, T, n, m)).Q))

# Hermite functions h₀ … h_nmax at `v`.
function hermite_functions(nmax, v)
    H = zeros(typeof(v), nmax + 2)
    H[1] = π^(-1 / 4) * exp(-v^2 / 2)
    H[2] = sqrt(2) * v * H[1]
    for n in 1:nmax
        H[n + 2] = sqrt(2 / (n + 1)) * v * H[n + 1] - sqrt(n / (n + 1)) * H[n]
    end
    return H[1:(nmax + 1)]
end
hermite(n) = v -> hermite_functions(n, v)[n + 1]

# The clamped velocity basis of degree p on N cells of [-V, V], and the periodic basis in x.
function velocity_basis(::Type{T}, N, p; V = 5) where {T}
    BSplineBasis(UniformMesh(N, (T(-V), T(V))), p)
end
function space_basis(::Type{T}, N, p) where {T}
    PeriodicBSplineBasis(UniformMesh(N, (zero(T), T(2π))), p)
end

@testset "Spline parity factors" begin
    for T in REAL_ELTYPES, p in 1:3

        N = isodd(p) ? 9 : 8
        bx, bv = space_basis(T, 12, p), velocity_basis(T, N, p)
        Nx, Nv = nbasis(bx), nbasis(bv)
        mx, mv = 3, 2
        layer = build_spline_parity_factors(bx, bv, mx, mv)
        @test layer isa BilinearStiefelLayer
        Rx, Rv = layer.metric_factors
        @test eltype(Rx) === T && eltype(Rv) === T
        Mx = Matrix(mass_matrix(SplineQuadrature(bx)))
        Mv = Matrix(mass_matrix(SplineQuadrature(bv)))
        J = reverse(Matrix{T}(I, Nv, Nv); dims = 1)
        tol = 30 * eps(T)

        # M = RᵀR for both factors, and the velocity factor commutes with the reflection
        @test maximum(abs, Rx' * Rx - Mx) ≤ tol * maximum(abs, Mx)
        @test maximum(abs, Rv' * Rv - Mv) ≤ tol * maximum(abs, Mv)
        @test maximum(abs, J * Rv * J - Rv) ≤ tol * maximum(abs, Rv)
        @test istriu(Rx)

        # The decoded modes are orthonormal in M: Φxᵀ Mx Φx = I and Φvᵀ M± Φv = I
        rng = StableRNG(1)
        ps = (Φx = stiefel(rng, T, Nx, mx), Φv = stiefel(rng, T, Nv ÷ 2, mv))
        S = parity_split(T, Nv)
        U = Rx \ Matrix(ps.Φx)
        cp = Rv \ (S[1:(Nv ÷ 2), :]' * Matrix(ps.Φv))
        cm = Rv \ (S[(Nv ÷ 2 + 1):end, :]' * Matrix(ps.Φv))
        @test opnorm(U' * Mx * U - I) ≤ tol
        @test opnorm(cp' * Mv * cp - I) ≤ tol
        @test opnorm(cm' * Mv * cm - I) ≤ tol
        @test opnorm(cp' * Mv * cm) ≤ tol

        # An even latent state decodes to an even function, an odd one to an odd function
        dec = inverse_layer(layer)
        q, p0 = randn(rng, T, mx * mv), zeros(T, mx * mv)
        Ce = dec((q = q, p = p0), ps)
        Co = dec((q = p0, p = q), ps)
        @test eltype(Ce) === T
        @test maximum(abs, Ce * J - Ce) ≤ tol * maximum(abs, Ce)
        @test maximum(abs, Co * J + Co) ≤ tol * maximum(abs, Co)

        # The layer encodes the L² inner products of the spline with the decoded modes
        C = randn(rng, T, Nx, Nv)
        z = layer(C, ps)
        @test z.q ≈ vec(U' * Mx * C * Mv * cp) rtol = 100 * eps(T)
        @test z.p ≈ vec(U' * Mx * C * Mv * cm) rtol = 100 * eps(T)
    end

    @test_throws r"^ArgumentError: the velocity domain \[-5\.0, 6\.0\] is not symmetric: v_min must equal -v_max$" build_spline_parity_factors(
        space_basis(Float64, 8, 2), BSplineBasis(UniformMesh(8, (-5.0, 6.0)), 2), 2, 2)
    @test_throws r"^ArgumentError: the velocity knot vector is not symmetric under v ↦ -v$" build_spline_parity_factors(
        space_basis(Float64, 8, 2),
        BSplineBasis(GeneralMesh([-5.0, -1.0, 0.5, 2.0, 5.0]), 2), 2, 2)
    @test_throws r"^ArgumentError: the velocity basis has 11 functions, an odd number: the parity split needs an even count, so pick N cells with N \+ p even$" build_spline_parity_factors(
        space_basis(Float64, 8, 2), velocity_basis(Float64, 9, 2), 2, 2)
    # the velocity basis must be clamped: a periodic one has no reflection j ↦ N + 1 - j
    @test_throws MethodError build_spline_parity_factors(
        space_basis(Float64, 8, 2), PeriodicBSplineBasis(UniformMesh(8, (-5.0, 5.0)), 2), 2, 2)

    # symmetric to 10 eps(T) v_max: v_max = 5 + 8 eps(T) 5 passes and 5 + 12 eps(T) 5 does
    # not; both endpoints and the knots in between are exact
    for T in REAL_ELTYPES
        vmax(f) = T(5) + f * eps(T) * T(5)
        knots(f) = GeneralMesh(T[-5, -2.5, 0, 2.5, vmax(f)])
        @test build_spline_parity_factors(
            space_basis(T, 8, 2), BSplineBasis(knots(8), 2), 2,
            2) isa BilinearStiefelLayer
        @test_throws r"^ArgumentError: the velocity domain \[-5\.0, 5\.0\d+\] is not symmetric: v_min must equal -v_max$" build_spline_parity_factors(
            space_basis(T, 8, 2), BSplineBasis(knots(12), 2), 2, 2)
    end
end

# The latent difference of f(x, v) = (1 + 0.1 cos x) M(v) given as its spline interpolant, its
# grid values and the rank-2 truncated SVD of those, after the tie. The grid has the knots in x
# and the midpoints of 2N half cells in v, so that the knots are cell faces of the grid.
function tie_errors(p, N; V = 8.0, mx = 4, mv = 2)
    bx, bv = space_basis(Float64, N, p), velocity_basis(Float64, N, p; V = V)
    Nv = nbasis(bv)
    spline = build_spline_parity_factors(bx, bv, mx, mv)
    Rx, Rv = spline.metric_factors
    qx, qv = SplineQuadrature(bx), SplineQuadrature(bv)
    # smooth modes: projections of smooth functions, orthonormal in the layer's coordinates
    X = reduce(hcat, l2_projection(qx, f) for f in (x -> 1.0, cos, sin, x -> cos(2x)))
    Ve = reduce(hcat, l2_projection(qv, hermite(n)) for n in (0, 2))
    S = parity_split(Float64, Nv)
    ps_s = (Φx = StiefelManifold(Matrix(qr(Rx * X).Q)[:, 1:mx]),
        Φv = StiefelManifold(Matrix(qr(S[1:(Nv ÷ 2), :] * (Rv * Ve)).Q)[:, 1:mv]))

    gx = range(0.0, 2π; length = N + 1)[1:N]
    gv = range(-V + V / (2N), V - V / (2N); length = 2N)
    grid = build_grid_parity_factors(gx, gv, mx, mv)
    dlra = LowRankStiefelLayer(N, 2N, mx, mv; metric_factors = grid.metric_factors)
    rng = StableRNG(2)
    ps_g = (Φx = stiefel(rng, Float64, N, mx), Φv = stiefel(rng, Float64, N, mv))
    ps_d = (Φx = stiefel(rng, Float64, N, mx), Φv = stiefel(rng, Float64, N, mv))
    tie_continuum_frontends!((spline, ps_s), (grid, ps_g), (dlra, ps_d),
        evaluation_matrix(bx, gx), evaluation_matrix(bv, gv))

    C = interpolant(bx, x -> 1 + 0.1cos(x)) * interpolant(bv, maxwellian)'
    F = [(1 + 0.1cos(x)) * maxwellian(v) for x in gx, v in gv]
    s = svd(F)
    zs = spline(C, ps_s)
    zg = grid(F, ps_g)
    zd = dlra((X = s.U[:, 1:2], S = Matrix(Diagonal(s.S[1:2])), V = s.V[:, 1:2]), ps_d)
    dist(a, b) = sqrt(sum(abs2, a.q - b.q) + sum(abs2, a.p - b.p))
    return dist(zs, zg), dist(zs, zd), sqrt(sum(abs2, zs.q))
end

@testset "Grid factors and the tie" begin
    for T in REAL_ELTYPES
        gx = range(zero(T), T(2π); length = 17)[1:16]
        gv = range(T(-4.75), T(4.75); length = 20)
        layer = build_grid_parity_factors(gx, gv, 3, 2)
        @test layer isa BilinearStiefelLayer
        Wx, Wv = layer.metric_factors
        @test eltype(Wx) === T && eltype(Wv) === T
        @test Wx ≈ Diagonal(fill(sqrt(step(gx)), 16)) rtol = 8 * eps(T)
        @test Wv ≈ Diagonal(fill(sqrt(step(gv)), 20)) rtol = 8 * eps(T)
        # each value carries the cell product √(Δx Δv)
        @test Wx[1, 1] * Wv[1, 1] ≈ sqrt(step(gx) * step(gv)) rtol = 8 * eps(T)
    end

    @test_throws r"^ArgumentError: the velocity domain \[-4\.75, 5\.0\] is not symmetric: v_min must equal -v_max$" build_grid_parity_factors(
        range(0.0, 1.0; length = 8), range(-4.75, 5.0; length = 20), 2, 2)
    @test_throws r"^ArgumentError: the velocity grid has 21 points, an odd number: the parity split needs an even count$" build_grid_parity_factors(
        range(0.0, 1.0; length = 8), range(-5.0, 5.0; length = 21), 2, 2)
    # symmetric to 10 eps(T) v_max: v_max = 5 + 8 eps(T) 5 passes and 5 + 12 eps(T) 5 does
    # not; both endpoints are exact, and `range` keeps them
    for T in REAL_ELTYPES
        vgrid(f) = range(T(-5), T(5) + f * eps(T) * T(5); length = 20)
        @test build_grid_parity_factors(range(zero(T), one(T); length = 8), vgrid(8), 2, 2) isa
              BilinearStiefelLayer
        @test_throws r"^ArgumentError: the velocity domain \[-5\.0, 5\.0\d+\] is not symmetric: v_min must equal -v_max$" build_grid_parity_factors(
            range(zero(T), one(T); length = 8), vgrid(12), 2, 2)
    end

    # The tie writes orthonormal bases into the grid and the low-rank parameters, and refuses a
    # low-rank layer of other factors, parameters of another size and evaluation matrices of
    # another size
    for T in REAL_ELTYPES
        N, p, mx, mv = 16, 2, 3, 2
        bx, bv = space_basis(T, N, p), velocity_basis(T, N, p)
        spline = build_spline_parity_factors(bx, bv, mx, mv)
        gx = range(zero(T), T(2π); length = N + 1)[1:N]
        gv = range(T(-5 + 5 / (2N)), T(5 - 5 / (2N)); length = 2N)
        grid = build_grid_parity_factors(gx, gv, mx, mv)
        dlra = LowRankStiefelLayer(N, 2N, mx, mv; metric_factors = grid.metric_factors)
        rng = StableRNG(8)
        ps_s = (
            Φx = stiefel(rng, T, nbasis(bx), mx), Φv = stiefel(rng, T, nbasis(bv) ÷ 2, mv))
        ps_g = (Φx = stiefel(rng, T, N, mx), Φv = stiefel(rng, T, N, mv))
        ps_d = (Φx = Matrix(stiefel(rng, T, N, mx)), Φv = Matrix(stiefel(rng, T, N, mv)))
        before = Matrix(ps_g.Φx)
        Ex, Ev = evaluation_matrix(bx, gx), evaluation_matrix(bv, gv)
        tie_continuum_frontends!((spline, ps_s), (grid, ps_g), (dlra, ps_d), Ex, Ev)
        @test eltype(ps_g.Φx) === T && eltype(ps_d.Φv) === T
        @test Matrix(ps_g.Φx) != before
        @test Matrix(ps_g.Φx) == ps_d.Φx && Matrix(ps_g.Φv) == ps_d.Φv
        @test opnorm(ps_g.Φx' * ps_g.Φx - I) ≤ 100 * eps(T)
        @test opnorm(ps_g.Φv' * ps_g.Φv - I) ≤ 100 * eps(T)
        wide = (Φx = Matrix(stiefel(rng, T, N, mx + 1)), Φv = ps_d.Φv)
        @test_throws r"^DimensionMismatch: ps_dlra\.Φx must have size \(16, 3\); got \(16, 4\)$" tie_continuum_frontends!(
            (spline, ps_s), (grid, ps_g), (dlra, wide), Ex, Ev)
        tall = (Φx = ps_g.Φx, Φv = stiefel(rng, T, N + 1, mv))
        @test_throws r"^DimensionMismatch: ps_grid\.Φv must have size \(16, 2\); got \(17, 2\)$" tie_continuum_frontends!(
            (spline, ps_s), (grid, tall), (dlra, ps_d), Ex, Ev)
        other = LowRankStiefelLayer(N, 2N, mx, mv)
        @test_throws r"^ArgumentError: the DLRA layer must carry the metric factors of the grid layer$" tie_continuum_frontends!(
            (spline, ps_s), (grid, ps_g), (other, ps_d), Ex, Ev)
        @test_throws r"^DimensionMismatch: E_x must have size \(16, 16\); got \(15, 16\)$" tie_continuum_frontends!(
            (spline, ps_s), (grid, ps_g), (dlra, ps_d), Ex[1:(end - 1), :], Ev)
        @test_throws r"^DimensionMismatch: E_v must have size \(32, 18\); got \(32, 17\)$" tie_continuum_frontends!(
            (spline, ps_s), (grid, ps_g), (dlra, ps_d), Ex, Ev[:, 1:(end - 1)])
    end

    # The three frontends agree on f to the order of the spline interpolant: the fitted order of
    # the latent difference is at least p + 0.5 over N = 16 … 128 cells (N + p even).
    for p in 1:3
        Ns = p == 2 ? [16, 32, 64, 128] : [17, 33, 65, 129]
        errs = [tie_errors(p, N) for N in Ns]
        eg, ed = first.(errs), getindex.(errs, 2)
        @test fitted_order(Ns, eg) ≥ p + 0.5
        @test fitted_order(Ns, ed) ≥ p + 0.5
        # the rank-2 factorisation is exact, so the grid and the low-rank frontends agree
        @test maximum(abs, ed .- eg) ≤ 1e-12
    end
end

# ∫∫ f₀ {g, h} dx dv for each column g of `G` and h of `H` (coefficients on the tensor basis,
# x fastest), by Gauß quadrature of the nodal fields, and the same sum over absolute terms.
function bracket_quadrature(bx, bv, f0, G, H)
    qx = SplineQuadrature(bx; nq = 8, dmax = 1)
    qv = SplineQuadrature(bv; nq = 24, dmax = 1)
    Bx, dBx = Matrix(basis_values(qx, 0)), Matrix(basis_values(qx, 1))
    Bv, dBv = Matrix(basis_values(qv, 0)), Matrix(basis_values(qv, 1))
    W = quadrature_weights(qx) * (quadrature_weights(qv) .* f0.(quadrature_nodes(qv)))'
    function fields(c)
        C = reshape(c, nbasis(bx), nbasis(bv))
        return dBx' * C * Bv, Bx' * C * dBv
    end
    FG, FH = fields.(eachcol(G)), fields.(eachcol(H))
    P = [sum(W .* (a[1] .* b[2] .- a[2] .* b[1])) for a in FG, b in FH]
    A = [sum(abs.(W) .* (abs.(a[1] .* b[2]) .+ abs.(a[2] .* b[1]))) for a in FG, b in FH]
    return P, A
end

@testset "Dual parity modes" begin
    for T in REAL_ELTYPES,
        (f0, df0) in ((maxwellian, dmaxwellian), (two_stream, dtwo_stream))

        bx, bv = space_basis(T, 16, 3), velocity_basis(T, 20, 3; V = 10)
        qx, qv = SplineQuadrature(bx), SplineQuadrature(bv; nq = 24)
        xm = [Spline(bx, l2_projection(qx, f))
              for f in (cos, sin, x -> cos(2x), x -> sin(2x))]
        ev = [Spline(bv, l2_projection(qv, hermite(n))) for n in (0, 2, 4)]
        od = [Spline(bv, l2_projection(qv, hermite(n))) for n in (1, 3, 5)]
        H = dual_parity_modes(xm, ev, od, df0)
        @test eltype(H) === T
        @test size(H) == (nbasis(bx) * nbasis(bv), 12)

        # the oracle: Float64 quadrature of ∫∫ f₀ {g, h} on the nodal fields
        bx64, bv64 = space_basis(Float64, 16, 3), velocity_basis(Float64, 20, 3; V = 10)
        X = Float64.(reduce(hcat, coefficients.(xm)))
        G = kron(Float64.(reduce(hcat, coefficients.(ev))), X)
        E = kron(Float64.(reduce(hcat, coefficients.(od))), X)
        B, _ = bracket_quadrature(bx64, bv64, f0, G, E)
        QP, _ = bracket_quadrature(bx64, bv64, f0, G, Float64.(H))
        PP, _ = bracket_quadrature(bx64, bv64, f0, Float64.(H), Float64.(H))
        # A_oo = ∫∫ |f₀| (|∂ₓe ∂ᵥe'| + |∂ᵥe ∂ₓe'|) on the odd modes e: the scale of the
        # round-off in {P, P}, since the bracket B_oo itself is zero by parity
        _, Aoo = bracket_quadrature(bx64, bv64, f0, E, E)
        @test opnorm(QP - I) ≤ 100 * eps(T) * cond(B)
        @test opnorm(PP) ≤ 100 * eps(T) * opnorm(inv(B))^2 * opnorm(Aoo)
        # the even modes paired with the L²-orthonormal odd modes are far from canonical
        @test opnorm(B - I) > 1
    end

    bx, bv = space_basis(Float64, 16, 3), velocity_basis(Float64, 20, 3; V = 10)
    qx, qv = SplineQuadrature(bx), SplineQuadrature(bv; nq = 24)
    xs(fs) = [Spline(bx, l2_projection(qx, f)) for f in fs]
    vs(ns) = [Spline(bv, l2_projection(qv, hermite(n))) for n in ns]

    # the keyword nq reaches the quadrature of K_v: two nodes per cell are too few, and the
    # dual modes fail the bound on {Q, P} above
    xm, ev, od = xs((cos, sin, x -> cos(2x), x -> sin(2x))), vs((0, 2, 4)), vs((1, 3, 5))
    X = reduce(hcat, coefficients.(xm))
    G = kron(reduce(hcat, coefficients.(ev)), X)
    B, _ = bracket_quadrature(
        bx, bv, maxwellian, G, kron(reduce(hcat, coefficients.(od)), X))
    QP, _ = bracket_quadrature(
        bx, bv, maxwellian, G, dual_parity_modes(xm, ev, od, dmaxwellian; nq = 2))
    @test opnorm(QP - I) > 100 * eps() * cond(B)

    @test_throws r"^ArgumentError: an x mode is constant: ∂x is not invertible on the x space$" dual_parity_modes(
        xs((x -> 1.0, cos, sin)), vs((0, 2)), vs((1, 3)), dmaxwellian)
    @test_throws r"^ArgumentError: ∂x is singular on the x space: each cos mx needs its sin mx$" dual_parity_modes(
        xs((cos, x -> cos(2x))), vs((0, 2)), vs((1, 3)), dmaxwellian)
    @test_throws r"^ArgumentError: the even and odd velocity spaces have dimensions 2 and 3: they must be equal$" dual_parity_modes(
        xs((cos, sin)), vs((0, 2)), vs((1, 3, 5)), dmaxwellian)
    @test_throws r"^ArgumentError: K_v = ∫ f₀' ψ⁺ ψ⁻ is numerically singular: cond\(K_v\) > 1/\(100 eps\)$" dual_parity_modes(
        xs((cos, sin)), vs((0, 2)), vs((1, 3)), v -> 0.0)
    bx2 = space_basis(Float64, 12, 3)
    @test_throws r"^ArgumentError: the x modes must share one basis$" dual_parity_modes(
        [xs((cos,)); Spline(bx2, l2_projection(SplineQuadrature(bx2), sin))],
        vs((0,)), vs((1,)), dmaxwellian)
    # x modes on a clamped basis
    bxc = BSplineBasis(UniformMesh(16, (0.0, 2π)), 3)
    @test_throws r"^ArgumentError: the x modes must be on a periodic basis$" dual_parity_modes(
        [Spline(bxc, l2_projection(SplineQuadrature(bxc), f)) for f in (cos, sin)],
        vs((0, 2)), vs((1, 3)), dmaxwellian)
    # odd modes on another velocity basis of the same size
    bv2 = velocity_basis(Float64, 20, 3; V = 8)
    @test_throws r"^ArgumentError: the odd velocity modes must share the basis of the even ones$" dual_parity_modes(
        xs((cos, sin)), vs((0,)),
        [Spline(bv2, l2_projection(SplineQuadrature(bv2; nq = 24), hermite(1)))], dmaxwellian)

    for T in REAL_ELTYPES
        bx, bv = space_basis(T, 16, 3), velocity_basis(T, 20, 3; V = 10)
        qx, qv = SplineQuadrature(bx), SplineQuadrature(bv; nq = 24)
        xm = [Spline(bx, l2_projection(qx, f)) for f in (cos, sin)]
        ev = [Spline(bv, l2_projection(qv, hermite(n))) for n in (0, 2)]
        c1, c3 = l2_projection(qv, hermite(1)), l2_projection(qv, hermite(3))
        # odd modes ψ₁⁻ and ψ₁⁻ + δ eps(T) ψ₃⁻: cond(K_v) · 100 eps(T) is 1.259 (Float32) and
        # 1.255 (Float64) at δ = 410, and 0.796 and 0.793 at δ = 650, so the two cases bracket
        # the threshold 1/(100 eps(T)) to within about 25 %
        od(δ) = [Spline(bv, c1), Spline(bv, c1 .+ δ * eps(T) .* c3)]
        @test_throws r"^ArgumentError: K_v = ∫ f₀' ψ⁺ ψ⁻ is numerically singular: cond\(K_v\) > 1/\(100 eps\)$" dual_parity_modes(
            xm, ev, od(410), dmaxwellian)
        @test size(dual_parity_modes(xm, ev, od(650), dmaxwellian)) ==
              (nbasis(bx) * nbasis(bv), 4)
        # a df0 that returns Float64 gives a result of the element type of the modes
        df64(v) = Float64(dmaxwellian(v))
        @test eltype(dual_parity_modes(xm, ev, od(650), df64)) === T
    end
end

# Particles of density ∝ (1 + 0.5 cos x) e^{-v²/2} on [0, 2π) × (-V, V), with positive weights
# of sum one.
function particles(rng, ::Type{T}, n; V = 5.5) where {T}
    x, v = T[], T[]
    while length(x) < n
        xa = 2π * rand(rng)
        rand(rng) < (1 + 0.5cos(xa)) / 1.5 || continue
        va = randn(rng)
        abs(va) < V || continue
        push!(x, T(xa))
        push!(v, T(va))
    end
    u = 0.5 .+ rand(rng, n)
    return x, v, T.(u ./ sum(u))
end

# Smooth observables a(x) b(v) with a closed-form bracket.
const PHI = (x -> sin(x) + 0.5cos(2x), v -> (1 + v) * exp(-v^2))
const PSI = (x -> cos(x) + 0.3sin(2x), v -> v^2 * exp(-v^2))
const CHI = (x -> cos(x) + sin(3x), v -> (0.5 + v - v^2) * exp(-v^2))

d(f) = y -> ForwardDiff.derivative(f, y)

"{φ, ψ} = ∂ₓφ ∂ᵥψ − ∂ᵥφ ∂ₓψ for separable φ = a(x) b(v), ψ = c(x) e(v), at (x, v)."
function exact_bracket(φ, ψ, x, v)
    (a, b), (c, e) = φ, ψ
    return d(a)(x) * c(x) * b(v) * d(e)(v) - a(x) * d(c)(x) * d(b)(v) * e(v)
end

# The spline machinery of the oracle, in Float64: the quadrature tables of both bases, their
# mass matrices, and the projection, evaluation and bracket of coefficient matrices.
struct Discretisation{BX, BV, QX, QV}
    bx::BX
    bv::BV
    qx::QX
    qv::QV
    Mx::Matrix{Float64}
    Mv::Matrix{Float64}
end
function Discretisation(N, p; V = 6.0)
    bx, bv = space_basis(Float64, N, p), velocity_basis(Float64, N, p; V = V)
    qx = SplineQuadrature(bx; nq = 2p + 2, dmax = 1)
    qv = SplineQuadrature(bv; nq = 2p + 2, dmax = 1)
    return Discretisation(bx, bv, qx, qv, Matrix(mass_matrix(qx)), Matrix(mass_matrix(qv)))
end
tables(q) = (Matrix(basis_values(q, 0)), Matrix(basis_values(q, 1)))
"The coefficients of the L² projection of the nodal field `F` on the quadrature grid."
function project(D, F)
    (Bx, _), (Bv, _) = tables(D.qx), tables(D.qv)
    return (D.Mx \
            (Bx * (quadrature_weights(D.qx) .* F .* quadrature_weights(D.qv)') * Bv')) /
           D.Mv
end
function project(D, φ::Tuple)
    project(D, φ[1].(quadrature_nodes(D.qx)) * φ[2].(quadrature_nodes(D.qv))')
end
"{A, B} of two splines, as a nodal field on the quadrature grid."
function bracket_field(D, A, B)
    (Bx, dBx), (Bv, dBv) = tables(D.qx), tables(D.qv)
    return (dBx' * A * Bv) .* (Bx' * B * dBv) .- (Bx' * A * dBv) .* (dBx' * B * Bv)
end
l2norm(D, C) = sqrt(dot(C, D.Mx * C * D.Mv))
"The value at each particle of the spline with coefficient matrix `C`."
function point_values(D, C, x, v)
    Ex, Ev = evaluation_matrix(D.bx, x), evaluation_matrix(D.bv, v)
    return vec(sum((Ex * C) .* Ev; dims = 2))
end

@testset "Klimontovich deposition" begin
    # (a) ⟨Π_h J, g⟩ = ⟨J, Π_h g⟩: the deposition returns the coefficients M⁻¹ c
    for T in REAL_ELTYPES, p in 1:3

        x, v, w = particles(StableRNG(3), T, 2000)
        dep = klimontovich_deposition(space_basis(T, 16, p), velocity_basis(T, 16, p; V = 6), w)
        ct = dep(x, v)
        @test eltype(ct) === T
        D = Discretisation(16, p)
        @test size(ct) == (nbasis(D.bx), nbasis(D.bv))
        Pg = project(D, PHI)
        G = D.Mx * Pg * D.Mv                    # the moments ∫ Λ g of g = φ
        lhs = dot(Float64.(ct), G)
        rhs = sum(Float64.(w) .* point_values(D, Pg, Float64.(x), Float64.(v)))
        gmax = maximum(abs, PHI[1].(range(0, 2π; length = 1001))) *
               maximum(abs, PHI[2].(range(-6, 6; length = 1001)))
        @test abs(lhs - rhs) ≤ 50 * eps(T) * sum(Float64.(w)) * gmax
        @test_throws r"^DimensionMismatch: x, v and the weights must have one entry per particle; got 1999, 2000 and 2000$" dep(
            x[1:(end - 1)], v)
        # a particle outside the velocity domain adds nothing
        out = copy(v)
        out[1:2] .= (T(-7), T(9))
        keep = [zeros(T, 2); w[3:end]]
        dep_keep = klimontovich_deposition(
            space_basis(T, 16, p), velocity_basis(T, 16, p; V = 6), keep)
        @test dep(x, out) ≈ dep_keep(x, v) rtol = 100 * eps(T)
        # the x basis must be periodic and the v basis clamped
        @test_throws MethodError klimontovich_deposition(
            velocity_basis(T, 16, p), velocity_basis(T, 16, p; V = 6), w)
        @test_throws MethodError klimontovich_deposition(
            space_basis(T, 16, p), space_basis(T, 16, p), w)
    end

    # (b) the truncation error E_h of the bracket of φ and ψ falls at order p + 1/2 or better
    x, v, w = particles(StableRNG(4), Float64, 2000)
    exact = sum(w .* exact_bracket.(Ref(PHI), Ref(PSI), x, v))
    for p in 1:3
        Ns = p == 2 ? [16, 32, 64, 128] : [17, 33, 65, 129]
        E = map(Ns) do N
            D = Discretisation(N, p)
            ct = klimontovich_deposition(D.bx, D.bv, w)(x, v)
            gh = bracket_field(D, project(D, PHI), project(D, PSI))
            G = D.Mx * project(D, gh) * D.Mv     # ∫ Λ {φ_h, ψ_h}
            abs(exact - dot(ct, G))
        end
        @test fitted_order(Ns, E) ≥ p + 0.5
    end

    # (c) the truncated bracket [φ, ψ]_h = Π_h {φ, ψ} violates the Jacobi identity at every N,
    # by more than 10³ times its round-off bound, and the defect decreases with N
    for p in 1:3
        Ns = p == 2 ? [16, 32, 64, 128] : [17, 33, 65, 129]
        rel = map(Ns) do N
            D = Discretisation(N, p)
            hb(A, B) = project(D, bracket_field(D, A, B))
            A, B, C = project(D, PHI), project(D, PSI), project(D, CHI)
            AB = hb(hb(A, B), C)
            Jh = AB .+ hb(hb(B, C), A) .+ hb(hb(C, A), B)
            bound = cond(D.Mx) * cond(D.Mv) * eps() * l2norm(D, AB)
            @test l2norm(D, Jh) > 1e3 * bound
            l2norm(D, Jh) / (l2norm(D, A) * l2norm(D, B) * l2norm(D, C))
        end
        @test all(diff(rel) .< 0)
    end
end

@testset "Klimontovich bridge loss" begin
    for T in REAL_ELTYPES
        bx, bv = space_basis(T, 12, 2), velocity_basis(T, 12, 2; V = 6)
        x, v, w = particles(StableRNG(5), T, 500)
        dep = klimontovich_deposition(bx, bv, w)
        layer = build_spline_parity_factors(bx, bv, 3, 2)
        loss = klimontovich_bridge_loss(dep, layer)
        rng = StableRNG(6)
        ps = (Φx = stiefel(rng, T, nbasis(bx), 3), Φv = stiefel(rng, T, nbasis(bv) ÷ 2, 2))
        z = (q = randn(rng, T, 6), p = randn(rng, T, 6))
        # the oracle feeds the spline frontend M_x⁻¹ c M_v⁻¹, the coefficients
        D = Discretisation(12, 2)
        c = evaluation_matrix(D.bx, Float64.(x))' * Diagonal(Float64.(w)) *
            evaluation_matrix(D.bv, Float64.(v))
        zs = layer(T.((D.Mx \ c) / D.Mv), ps)
        expected = (sum(abs2, z.q - zs.q) + sum(abs2, z.p - zs.p)) /
                   (sum(abs2, zs.q) + sum(abs2, zs.p) + eps(T))
        value = @inferred loss(z, x, v, ps)
        @test value isa T
        @test @inferred(dep(x, v)) isa Matrix{T}
        @test value ≈ expected rtol = sqrt(eps(T))
    end
end

# The linearised Vlasov–Poisson right-hand side -v ∂ₓδf + f₀' ∂ₓΦ, -∂ₓ²Φ = ρ, the free energy
# by quadrature in x, on the modes cos kx/√π, sin kx/√π and the nodes `vs` of width Δv.
function vlasov_poisson_reference(C, vs, ks, df0; nx = 64)
    Δv = vs[2] - vs[1]
    xs = 2π .* (0:(nx - 1)) ./ nx
    modes = [x -> (iseven(i) ? sin : cos)(ks[(i + 1) ÷ 2] * x) / sqrt(π)
             for i in 1:(2length(ks))]
    f(x, q) = sum(C[i, q] * modes[i](x) for i in eachindex(modes))
    ρ(x) = Δv * sum(f(x, q) for q in eachindex(vs))
    project(g) = [2π / nx * sum(g(x) * m(x) for x in xs) for m in modes]
    Φ = project(ρ) ./ [ks[(i + 1) ÷ 2]^2 for i in eachindex(modes)]
    Φx(x) = ForwardDiff.derivative(y -> sum(Φ[i] * modes[i](y) for i in eachindex(modes)), x)
    rhs = reduce(hcat,
        project(x -> -vs[q] * ForwardDiff.derivative(y -> f(y, q), x) + df0(vs[q]) * Φx(x))
        for q in eachindex(vs))
    w0 = -vs ./ df0.(vs)
    value = (Δv * sum(w0[q] * 2π / nx * sum(f(x, q)^2 for x in xs) for q in eachindex(vs)) +
             2π / nx * sum(Φx(x)^2 for x in xs)) / 2
    return rhs, value
end

@testset "Kruskal–Oberman free energy" begin
    for T in REAL_ELTYPES
        vs = range(T(-7.75), T(7.75); length = 32)
        ks = [1, 2, 3]
        H = kruskal_oberman_free_energy(dmaxwellian, vs, ks)
        C = randn(StableRNG(7), T, 6, 32) .* maxwellian.(vs)'
        r = @inferred H(C)
        @test r.value isa T && eltype(r.gradient) === T && eltype(r.rhs) === T
        rhs, value = vlasov_poisson_reference(Float64.(C), Float64.(vs), ks, dmaxwellian)
        @test r.value ≈ value rtol = 50 * eps(T)
        @test r.gradient ≈ ForwardDiff.gradient(C -> H(C).value, C) rtol = 50 * eps(T)
        # the vector field of the skew operator is the linearised Vlasov–Poisson right-hand side,
        # and the free energy is conserved along it
        @test r.rhs ≈ rhs rtol = 50 * eps(T)
        @test abs(dot(Float64.(r.gradient), rhs)) ≤
              50 * eps(T) * norm(r.gradient) * norm(rhs)
    end

    vs = range(-7.75, 7.75; length = 32)
    @test_throws r"^ArgumentError: the x space holds the constant mode k = 0, which has no conjugate: pass the wavenumbers k ≠ 0 only$" kruskal_oberman_free_energy(
        dmaxwellian, vs, [0, 1])
    @test_throws r"^ArgumentError: the wavenumbers are empty: the x space needs at least one mode$" kruskal_oberman_free_energy(
        dmaxwellian, vs, Int[])
    @test_throws r"^ArgumentError: v f₀'\(v\) ≥ 0 at a node, where the weight -v/f₀' of the free energy must be positive; for a non-monotone f₀ use dual_parity_modes$" kruskal_oberman_free_energy(
        dtwo_stream, vs, [1, 2])
    # f₀' = 0 at the nodes |v| > 7 is refused as well
    @test_throws r"^ArgumentError: v f₀'\(v\) ≥ 0 at a node, where the weight -v/f₀' of the free energy must be positive; for a non-monotone f₀ use dual_parity_modes$" kruskal_oberman_free_energy(
        v -> abs(v) > 7 ? zero(v) : dmaxwellian(v), vs, [1, 2])
    # the sign of f₀' is checked on v < 0 too: f₀' = -|v| and f₀' = 0 there are refused
    for df0 in (v -> -abs(v), v -> v < 0 ? zero(v) : dmaxwellian(v))
        @test_throws r"^ArgumentError: v f₀'\(v\) ≥ 0 at a node, where the weight -v/f₀' of the free energy must be positive; for a non-monotone f₀ use dual_parity_modes$" kruskal_oberman_free_energy(
            df0, vs, [1, 2])
    end
    # the signs are compared, not multiplied: v f₀' = 1f-46 underflows to zero in Float32
    Hs = kruskal_oberman_free_energy(v -> -v, range(-1.0f-23, 1.0f-23; length = 4), [1])
    @test Hs.w0 ≈ ones(Float32, 4) rtol = 8 * eps(Float32)
    # k and -k span the same modes, so each |k| occurs once
    for ks in ([1, 1], [1, -1], [2, 1, -2])
        @test_throws r"^ArgumentError: the wavenumbers repeat a mode: k and -k give the same modes cos kx, sin kx, so each \|k\| may occur once$" kruskal_oberman_free_energy(
            dmaxwellian, vs, ks)
    end
    # a negative wavenumber alone is a mode of its own, with the free energy of |k|
    C = randn(StableRNG(9), 4, 32) .* maxwellian.(vs)'
    @test kruskal_oberman_free_energy(dmaxwellian, vs, [-1, 2])(C).value ≈
          kruskal_oberman_free_energy(dmaxwellian, vs, [1, 2])(C).value rtol = 8 * eps()
    @test_throws r"^ArgumentError: the velocity grid has 31 points, an odd number: the parity split needs an even count$" kruskal_oberman_free_energy(
        dmaxwellian, range(-7.5, 7.5; length = 31), [1])
    @test_throws r"^ArgumentError: the velocity domain \[-7\.75, 8\.0\] is not symmetric: v_min must equal -v_max$" kruskal_oberman_free_energy(
        dmaxwellian, range(-7.75, 8.0; length = 32), [1])
    H = kruskal_oberman_free_energy(dmaxwellian, vs, [1])
    @test_throws r"^DimensionMismatch: the state must have size \(2, 32\); got \(2, 31\)$" H(zeros(2, 31))
end
