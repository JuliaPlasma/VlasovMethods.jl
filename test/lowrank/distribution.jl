using LinearAlgebra
using Random
using Test
using VlasovMethods
using VlasovMethods: AbstractFactorSpace, DistributionFunction, inner_v, inner_x,
                     orthonormalize_v, orthonormalize_x

# Every expected value comes from somewhere other than the code under test: `svdvals` of the
# weighted grid function (Eckart–Young), `velocity_moments` and `evaluate` of the
# `GridDistribution`, the identities of an orthonormal basis, or a literal. Each tolerance is a
# multiple of `eps(T)`.

maxwellian(v) = exp(-v^2 / 2) / sqrt(2π)

# a smooth grid function of numerical rank above 6 on the fixture grids: the shifted Maxwellian
# of the first term is not separable
function fA(x, v)
    k = 0.5
    maxwellian(v - 0.5sin(k * x)) * (1 + 0.3cos(k * x)) +
    0.2exp(-(v - 2)^2) * (1 + 0.4sin(2k * x)) +
    0.1v * exp(-v^2 / 2) * cos(3k * x + 0.7) + 0.05v^2 * exp(-v^2 / 2) * sin(k * x) +
    0.02exp(-(v + 1.5)^2) * cos(5k * x)
end

# a grid function of rank 2: the x factors 1 + 0.3 cos(x/2) and sin(x/2) are independent, and so
# are the v factors M(v) and v M(v)
fB(x, v) = maxwellian(v) * (1 + 0.3cos(0.5x)) + 0.1v * maxwellian(v) * sin(0.5x)

# f(x, ±vmax) is round-off, at most `eps(Float64)` of the largest value, for both fixtures
const XD = (0.0, 4π)
const VD = (-10.0, 10.0)

# a factor space that is not a grid: only its sizes, for the methods of the abstract type
struct ProbeSpace{DT} <: AbstractFactorSpace{DT}
    n::NTuple{2, Int}
end
Base.size(s::ProbeSpace) = s.n

function grid(T, f, nx, nv)
    d0 = GridDistribution(zeros(T, nx, nv), T.(XD), T.(VD))
    return GridDistribution(T[f(x, v) for x in d0.x, v in d0.v], T.(XD), T.(VD))
end

# the modes 1, v and v² on the velocity grid of `dist`, in the element type `S`
modes(dist, S = eltype(dist)) = (v = S.(collect(dist.v)); hcat(one.(v), v, v .^ 2))

# the largest entry of h QᵀQ - I, in units of eps(T)
orthogonality(Q, h) = maximum(abs, h * Q' * Q - I; init = zero(eltype(Q))) / eps(eltype(Q))

wnorm(A, dist) = norm(sqrt(step(dist.x) * step(dist.v)) .* A)

allocations(f::F, a::A) where {F, A} = (f(a); @allocated f(a))

# the exception that `f()` throws, or `nothing`
caught(f) =
    try
        f()
        nothing
    catch err
        err
    end

@testset "the fixtures are round-off at the velocity boundary" begin
    for f in (fA, fB)
        d = grid(Float64, f, 24, 40)
        @test maximum(abs, d.values[:, [1, end]]) ≤ eps() * maximum(abs, d.values)
    end
end

@testset "documented names" begin
    @test :LowRankDistribution in names(VlasovMethods)
    @test :GridFactorSpace in names(VlasovMethods)
    @test isdisjoint(
        (:AbstractFactorSpace, :inner_x, :inner_v, :orthonormalize_x,
            :orthonormalize_v),
        names(VlasovMethods))
    for name in (:LowRankDistribution, :GridFactorSpace, :AbstractFactorSpace, :inner_x,
        :inner_v, :orthonormalize_x, :orthonormalize_v)
        @test Base.Docs.hasdoc(VlasovMethods, name)
    end
end

@testset "LowRankDistribution, $T" for T in (Float64, Float32)
    e = eps(T)
    dist = grid(T, fA, 24, 40)
    hx, hv = step(dist.x), step(dist.v)
    @test hx isa T && hv isa T && !(hx ≈ hv)
    space = GridFactorSpace(dist)
    @test space.x == dist.x && space.v == dist.v
    Φ = modes(dist)

    @testset "the abstract factor space" begin
        @test space isa AbstractFactorSpace{T}
        @test size(space) == (24, 40) && eltype(space) === T
        d = LowRankDistribution(dist, 4)
        @test fieldtype(typeof(d), :space) === typeof(space)
        @test all(k -> isconcretetype(fieldtype(typeof(d), k)), fieldnames(typeof(d)))
        # the factor constructor reads the sizes through `size(space)`, on any space
        ps = ProbeSpace{T}((5, 7))
        dp = LowRankDistribution(ps, zeros(T, 5, 2), zeros(T, 2, 2), zeros(T, 7, 2))
        @test dp.space === ps && size(dp) == (5, 7) && eltype(dp) === T
        @test LowRankDistribution(ps, zeros(T, 5, 5), zeros(T, 5, 5), zeros(T, 7, 5)) isa
              LowRankDistribution{T}
        @test_throws r"^DimensionMismatch: factors of sizes" LowRankDistribution(ps,
            zeros(T, 7, 2), zeros(T, 2, 2), zeros(T, 5, 2))
        @test_throws r"^ArgumentError: the rank 6 is not in 1:5" LowRankDistribution(ps,
            zeros(T, 5, 6), zeros(T, 6, 6), zeros(T, 7, 6))
        # the grid functions have no method for a space that is not a grid
        nomethod(f) = Regex("^MethodError: no method matching (VlasovMethods\\.)?$f\\(")
        @test_throws nomethod("GridDistribution") GridDistribution(dp)
        @test_throws nomethod("velocity_moments") velocity_moments(dp)
        @test_throws nomethod("evaluate") evaluate(dp, T[0], T[0])
        @test_throws nomethod("inner_x") inner_x(ps, zeros(T, 5), zeros(T, 5))
        @test_throws nomethod("orthonormalize_v") orthonormalize_v(ps, zeros(T, 7, 2))
    end

    @testset "inner products" begin
        a, b = randn(Xoshiro(1), T, 24), randn(Xoshiro(2), T, 24)
        @test inner_x(space, a, b) ≈ hx * dot(a, b) rtol = 8e
        @test inner_x(space, a, b) isa T
        A, B = randn(Xoshiro(3), T, 40, 3), randn(Xoshiro(4), T, 40, 2)
        @test inner_v(space, A, B) ≈ hv * (A' * B) rtol = 8e
        @test eltype(inner_v(space, A, B)) === T
        C, D = randn(Xoshiro(15), T, 24, 3), randn(Xoshiro(16), T, 24, 2)
        @test inner_x(space, C, D) ≈ hx * (C' * D) rtol = 8e
        @test inner_v(space, A[:, 1], B[:, 1]) ≈ hv * dot(A[:, 1], B[:, 1]) rtol = 8e
        rows(n, axis, m) = Regex("^DimensionMismatch: a factor of $n rows on the $axis " *
                                 "grid of $m nodes")
        @test_throws rows(40, "x", 24) inner_x(space, zeros(T, 40), zeros(T, 40))
        @test_throws rows(40, "x", 24) inner_x(space, a, zeros(T, 40))
        @test_throws rows(40, "x", 24) inner_x(space, zeros(T, 40), a)
        @test_throws rows(24, "v", 40) inner_v(space, zeros(T, 24, 2), zeros(T, 24, 2))
        @test_throws rows(24, "v", 40) inner_v(space, A, zeros(T, 24, 2))
        @test_throws rows(24, "v", 40) inner_v(space, zeros(T, 24, 2), A)
        @test_throws rows(40, "x", 24) orthonormalize_x(space, zeros(T, 40, 2))
        @test_throws rows(24, "v", 40) orthonormalize_v(space, zeros(T, 24, 2))
    end

    @testset "orthonormalisation" begin
        for (orth, n, h) in ((orthonormalize_x, 24, hx), (orthonormalize_v, 40, hv))
            A = randn(Xoshiro(5), T, n, 5)
            Q, R = orth(space, A)
            @test size(Q) == (n, 5) && size(R) == (5, 5)
            @test eltype(Q) === T && eltype(R) === T
            @test orthogonality(Q, h) ≤ 50
            @test maximum(abs, Q * R - A) ≤ 50e * maximum(abs, A)
            Qz, Rz = orth(space, zeros(T, n, 3))
            @test size(Qz) == (n, 3) && size(Rz) == (3, 3)
            @test orthogonality(Qz, h) ≤ 50 && iszero(Rz)
            An = copy(A)
            An[2, 3] = T(NaN)
            _, Rn = orth(space, An)
            @test any(isnan, Rn)
            @test_throws r"^ArgumentError: cannot orthonormalise \d+ columns" orth(space,
                zeros(T, n, n + 1))
            Qs, Rs = orth(space, randn(Xoshiro(17), T, n, n))
            @test size(Qs) == (n, n) && size(Rs) == (n, n) && orthogonality(Qs, h) ≤ 50
        end
    end

    @testset "Eckart–Young without fixed modes" begin
        σ = svdvals(sqrt(hx * hv) .* dist.values)
        @test σ[7] > 1e-6 * σ[1]  # rank 6 truncates fA
        for r in (3, 4, 6)
            d = LowRankDistribution(dist, r)
            @test d isa DistributionFunction{T, 1, 1}
            @test eltype(d) === T && size(d) == (24, 40) && d.nfixed == 0
            @test size(d.X) == (24, r) && size(d.S) == (r, r) && size(d.V) == (40, r)
            @test eltype(d.X) === eltype(d.S) === eltype(d.V) === T
            @test orthogonality(d.X, hx) ≤ 50 && orthogonality(d.V, hv) ≤ 50
            G = GridDistribution(d)
            @test G.x == dist.x && G.v == dist.v && eltype(G.values) === T
            tail = sqrt(sum(abs2, σ[(r + 1):end]))
            @test abs(wnorm(dist.values - G.values, dist) - tail) ≤
                  50e * wnorm(dist.values, dist)
        end
    end

    @testset "fixed modes" begin
        # Float64 modes on a Float32 grid are converted
        for r in (3, 4, 6), fixed in (Φ, modes(dist, Float64))

            d = LowRankDistribution(dist, r; fixed)
            @test d.nfixed == 3
            @test eltype(d.X) === eltype(d.S) === eltype(d.V) === T
            @test orthogonality(d.X, hx) ≤ 50 && orthogonality(d.V, hv) ≤ 50
            for k in 1:3
                Vk, Φk = d.V[:, 1:k], Φ[:, 1:k]
                @test norm(Φk - Vk * (hv * Vk' * Φk)) ≤ 50e * norm(Φk)
            end
            # Eckart–Young for the part of f outside the fixed modes
            P = Matrix(qr(sqrt(hv) .* Φ).Q) ./ sqrt(hv)
            g = dist.values - dist.values * (hv * P * P')
            tail = sqrt(sum(abs2, svdvals(sqrt(hx * hv) .* g)[(r - 2):end]))
            F = GridDistribution(d).values
            @test abs(wnorm(dist.values - F, dist) - tail) ≤ 50e * wnorm(dist.values, dist)
        end
        # as many fixed modes as the rank
        d = LowRankDistribution(dist, 3; fixed = Φ)
        @test norm(Φ - d.V * (hv * d.V' * Φ)) ≤ 50e * norm(Φ)
    end

    @testset "the fixed modes keep the moments" begin
        m0 = velocity_moments(dist)
        for r in (3, 4, 6)
            m1 = velocity_moments(LowRankDistribution(dist, r; fixed = Φ))
            for s in 1:3
                @test eltype(m1[s]) === T && length(m1[s]) == 24
                @test maximum(abs, m1[s] - m0[s]) ≤ 1600e * maximum(abs, m0[s])
            end
        end
        # control: without the fixed modes, rank 3 changes the density far beyond the bound
        m3 = velocity_moments(LowRankDistribution(dist, 3))
        @test maximum(abs, m3.density - m0.density) > 16000e * maximum(abs, m0.density)
    end

    @testset "deterministic completion" begin
        distB = grid(T, fB, 24, 40)
        σB = svdvals(sqrt(hx * hv) .* distB.values)
        @test σB[2] > 1e-3 * σB[1] && σB[3] ≤ 100e * σB[1]
        for fixed in (zeros(T, 40, 0), modes(distB))
            Random.seed!(1)
            d1 = LowRankDistribution(distB, 6; fixed)
            after = rand()
            Random.seed!(2)
            d2 = LowRankDistribution(distB, 6; fixed)
            Random.seed!(1)
            @test rand() == after
            @test d1.X == d2.X && d1.S == d2.S && d1.V == d2.V
            @test orthogonality(d1.X, hx) ≤ 50 && orthogonality(d1.V, hv) ≤ 50
            @test maximum(abs, GridDistribution(d1).values - distB.values) ≤
                  100e * maximum(abs, distB.values)
        end
        # The completion comes from the Householder QR of the two leading columns: each column of
        # its orthogonal factor lies in the span of those columns and of the first r unit vectors,
        # where singular vectors of the round-off of fB would not.
        d = LowRankDistribution(distB, 6)
        E = hcat(d.V[:, 1:2], Matrix{T}(I, 40, 6))
        @test norm(d.V - E * (E \ d.V)) ≤ 50e * norm(d.V)
    end

    @testset "a mutable state" begin
        d = LowRankDistribution(dist, 4)
        @test ismutabletype(typeof(d))
        d6 = LowRankDistribution(dist, 6; fixed = Φ)
        d.X, d.S, d.V, d.nfixed = d6.X, d6.S, d6.V, d6.nfixed
        @test size(d.X) == (24, 6) && size(d.S) == (6, 6) && size(d.V) == (40, 6)
        @test d.nfixed == 3
        @test GridDistribution(d).values == GridDistribution(d6).values
        @test velocity_moments(d) == velocity_moments(d6)
    end

    @testset "rank min(nx, nv)" begin
        for (nx, nv) in ((24, 40), (40, 24))
            dd = grid(T, fA, nx, nv)
            d = LowRankDistribution(dd, min(nx, nv))
            @test maximum(abs, GridDistribution(d).values - dd.values) ≤
                  100e * maximum(abs, dd.values)
        end
        dd = grid(T, fA, 40, 24)
        d = LowRankDistribution(dd, 24; fixed = modes(dd))
        @test maximum(abs, GridDistribution(d).values - dd.values) ≤
              100e * maximum(abs, dd.values)
    end

    @testset "the zero matrix" begin
        dz = GridDistribution(zeros(T, 24, 40), T.(XD), T.(VD))
        for fixed in (zeros(T, 40, 0), Φ)
            d = LowRankDistribution(dz, 4; fixed)
            @test iszero(d.S) && iszero(GridDistribution(d).values)
            @test orthogonality(d.X, hx) ≤ 50 && orthogonality(d.V, hv) ≤ 50
            @test all(iszero, velocity_moments(d))
            @test d.nfixed == size(fixed, 2)
        end
        d = LowRankDistribution(dz, 4; fixed = Φ)
        @test norm(Φ - d.V * (hv * d.V' * Φ)) ≤ 50e * norm(Φ)
    end

    @testset "refused inputs" begin
        for r in (0, -1, 25)
            @test_throws Regex("^ArgumentError: the rank $r is not in 1:24") LowRankDistribution(
                dist, r)
        end
        @test LowRankDistribution(dist, 24) isa LowRankDistribution{T}
        @test_throws r"^ArgumentError: 3 fixed modes exceed the rank 2" LowRankDistribution(
            dist, 2; fixed = Φ)
        v = collect(dist.v)
        dependent = r"^ArgumentError: the fixed modes are linearly dependent"
        @test_throws dependent LowRankDistribution(dist, 4; fixed = hcat(one.(v), v, 2 .*
                                                                                     v))
        @test_throws dependent LowRankDistribution(dist, 4; fixed = hcat(one.(v), zero.(v)))
        # one fixed mode: a zero mode is refused, the constant is the first column of V
        @test_throws dependent LowRankDistribution(dist, 4; fixed = zeros(T, 40, 1))
        d1 = LowRankDistribution(dist, 4; fixed = ones(T, 40, 1))
        @test d1.nfixed == 1
        @test norm(d1.V[:, 1] * (hv * sum(d1.V[:, 1])) - ones(T, 40)) ≤ 50e * sqrt(T(40))
        @test_throws r"^DimensionMismatch: fixed modes of 39 rows" LowRankDistribution(
            dist, 4;
            fixed = Φ[1:39, :])
        for bad in (T(NaN), T(Inf))
            f = copy(dist.values)
            f[3, 4] = bad
            for fixed in (zeros(T, 40, 0), Φ)
                err = caught(() -> LowRankDistribution(
                    GridDistribution(dist.x, dist.v, f), 4;
                    fixed))
                @test err isa ArgumentError && occursin("non-finite", err.msg)
            end
            Φn = copy(Φ)
            Φn[5, 2] = bad
            err = caught(() -> LowRankDistribution(dist, 4; fixed = Φn))
            @test err isa ArgumentError && occursin("non-finite", err.msg)
        end
    end

    @testset "the factor constructor and the space refuse wrong sizes" begin
        @test_throws r"^ArgumentError: a periodic x grid needs at least one node" GridFactorSpace(
            range(zero(T), one(T); length = 0), dist.v)
        @test_throws r"^ArgumentError: a bounded velocity grid needs at least two nodes, got 1" GridFactorSpace(
            dist.x, dist.v[1:1])
        @test GridFactorSpace(dist.x[1:1], dist.v[1:2]) isa GridFactorSpace{T}
        X, S, V = zeros(T, 24, 2), zeros(T, 2, 2), zeros(T, 40, 2)
        @test LowRankDistribution(space, X, S, V).nfixed == 0
        @test LowRankDistribution(space, X, S, V; nfixed = 2).nfixed == 2
        for nfixed in (3, -1)
            @test_throws Regex("^ArgumentError: the number of fixed modes $nfixed is not in 0:2") LowRankDistribution(
                space, X, S, V; nfixed)
        end
        sizes = r"^DimensionMismatch: factors of sizes"
        @test_throws sizes LowRankDistribution(space, zeros(T, 23, 2), S, V)
        @test_throws sizes LowRankDistribution(space, X, zeros(T, 2, 3), V)
        @test_throws sizes LowRankDistribution(space, X, zeros(T, 3, 2), V)
        @test_throws sizes LowRankDistribution(space, X, S, zeros(T, 40, 3))
        @test_throws sizes LowRankDistribution(space, X, S, zeros(T, 39, 2))
        @test_throws r"^ArgumentError: the rank 0 is not in 1:24" LowRankDistribution(
            space,
            zeros(T, 24, 0), zeros(T, 0, 0), zeros(T, 40, 0))
        @test LowRankDistribution(space, zeros(T, 24, 24), zeros(T, 24, 24),
            zeros(T, 40, 24)) isa LowRankDistribution{T}
        @test_throws r"^ArgumentError: the rank 25 is not in 1:24" LowRankDistribution(
            space,
            zeros(T, 24, 25), zeros(T, 25, 25), zeros(T, 40, 25))
    end

    @testset "the reconstruction and the factored moments" begin
        X, _ = orthonormalize_x(space, randn(Xoshiro(6), T, 24, 4))
        V, _ = orthonormalize_v(space, randn(Xoshiro(7), T, 40, 4))
        S = randn(Xoshiro(8), T, 4, 4)
        d = LowRankDistribution(space, X, S, V)
        G = GridDistribution(d)
        @test G.x == space.x && G.v == space.v
        @test maximum(abs, G.values - X * S * V') ≤ 8e * maximum(abs, G.values)
        mg, mf = velocity_moments(G), velocity_moments(d)
        for s in 1:3
            @test eltype(mf[s]) === T
            @test maximum(abs, mf[s] - mg[s]) ≤ 50e * maximum(abs, mg[s])
        end
    end

    @testset "the moments do not form f" begin
        # a 512 × 512 grid holds 2 MiB of Float64 values
        big = GridDistribution(zeros(T, 512, 512), T.(XD), T.(VD))
        bs = GridFactorSpace(big)
        Xb, _ = orthonormalize_x(bs, randn(Xoshiro(9), T, 512, 8))
        Vb, _ = orthonormalize_v(bs, randn(Xoshiro(10), T, 512, 8))
        db = LowRankDistribution(bs, Xb, randn(Xoshiro(11), T, 8, 8), Vb)
        full = 512 * 512 * sizeof(T)
        @test allocations(velocity_moments, db) < full ÷ 8
        @test allocations(d -> velocity_moments(GridDistribution(d)), db) > full
    end

    @testset "point evaluation" begin
        d = LowRankDistribution(dist, 6; fixed = Φ)
        G = GridDistribution(d)
        L = T(4π)
        xs = T[0, 0.3, last(dist.x), last(dist.x) + 0.2, L - eps(L), 3L + 0.1, -2L + 1.3]
        vs = T[-10, 10, 0.1, -9.99, 2.5, first(dist.v) + 5step(dist.v)]
        # the given factors are not zero at the end nodes, so the last-cell rule at v = v_b
        # changes the value there; the factorisation of fA is round-off there
        Xg, _ = orthonormalize_x(space, randn(Xoshiro(12), T, 24, 4))
        Vg, _ = orthonormalize_v(space, randn(Xoshiro(13), T, 40, 4))
        dg = LowRankDistribution(space, Xg, randn(Xoshiro(14), T, 4, 4), Vg)
        @test minimum(abs, GridDistribution(dg).values[:, [1, end]]) > 0
        for dd in (d, dg)
            Gd = GridDistribution(dd)
            bound = 50e * maximum(abs, Gd.values)
            for x in xs, v in vs

                value = evaluate(dd, [x], [v])
                @test value isa T
                @test abs(value - evaluate(Gd, [x], [v])) ≤ bound
                @test dd(x, v) == value
                @test dd([x], [v]) == value
            end
        end
        for x in xs, v in T[-10.01, 10.01, -100, 100, NaN]

            @test evaluate(d, [x], [v]) === zero(T)
        end
        for x in T[NaN, Inf]
            # `Int` of the non-finite cell index, as in the evaluation of the grid distribution
            inexact = Regex("^InexactError: $(Int)\\($(repr(x))\\)")
            @test_throws inexact evaluate(G, [x], [T(0.1)])
            @test_throws inexact evaluate(d, [x], [T(0.1)])
        end
    end
end
