using HDF5: h5read
using LinearAlgebra: I, dot, istriu, norm, svd, svdvals
using PoissonSolvers: FFTWBasis, Potential
using Random: Xoshiro, rand, randn
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V, orthonormalize_v, orthonormalize_x

include("../helpers/dense_parity.jl")

# Every fixture lies on x ∈ [0, 4π), v ∈ [-7.5, 7.5], with `G` the Green's matrix of the model's
# FFT backend. The expected values come from the dense transcription of the algorithm in
# `test/helpers/dense_parity.jl`, from `ParitySplitting`, the full-rank parity splitting, or
# from the conservation laws of the formulation.

maxwellian(v) = exp(-v^2 / 2) / sqrt(2 * oftype(v, π))

# (1 + α cos(x/2)) M(v) (1 + β v) on nx × nv nodes, its model, its factorisation at `rank` and the
# workspace of the fields
function fixture(
        T; nx = 16, nv = 32, α = 0.2, β = 0.3, rank = (6, 6), background = nothing,
        f₀ = (x, v) -> (1 + T(α) * cos(x / 2)) * maxwellian(v) * (1 + T(β) * v))
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, nx, nv), xd, (T(-7.5), T(7.5)))
    dist.values .= [f₀(x, v) for x in dist.x, v in dist.v]
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    Y = LowRankParityDistribution(dist, rank; background)
    field = LowRankFields1D1V(Y.plus.space, model.work)
    return (; dist, model, Y, field, G = model.work)
end

rel(a, b) = norm(a .- b) / norm(b)
values_of(Y) = GridDistribution(Y).values
ranks(Y) = (size(Y.plus.S, 1), size(Y.minus.S, 1))
function factors(Y)
    (copy(Y.f0), copy(Y.plus.X), copy(Y.plus.S), copy(Y.plus.V), copy(Y.minus.X),
        copy(Y.minus.S), copy(Y.minus.V))
end
status(Y) = all(isfinite, Y.plus.S) && all(isfinite, Y.minus.S)

function orthonormality(o)
    max(norm(step(o.space.x) .* (o.X' * o.X) - I), norm(step(o.space.v) .* (o.V' * o.V) -
                                                        I))
end

even_bits(V) = V == V[end:-1:1, :]
odd_bits(V) = V == -V[end:-1:1, :]

# H = hx hv Σ f v²/2 + (hx/2) ρᵀ G ρ, and the mass hx hv Σ f
function energy(dist, G, f)
    hx, hv = step(dist.x), step(dist.v)
    ρ = hv .* vec(sum(f; dims = 2))
    return hx * hv * sum(f .* (collect(dist.v)' .^ 2 ./ 2)) + (hx / 2) * dot(ρ, G * ρ)
end

# the fixed modes span {1, v²} (f⁺) and v (f⁻)
function fixed_mode_error(Y)
    v = collect(Y.plus.space.v)
    hv = step(Y.plus.space.v)
    proj(Vf, Φ) = norm(Φ - Vf * (hv .* (Vf' * Φ))) / norm(Φ)
    return max(proj(Y.plus.V[:, 1:2], hcat(one.(v), v .^ 2)),
        proj(Y.minus.V[:, 1:1], reshape(v, :, 1)))
end

@testset "factorisation: full rank reproduces the values, low rank keeps the moments, $T" for T in (
    Float32, Float64)
    for rank in ((16, 16), (4, 3))
        (; dist, Y) = fixture(T; rank)
        f = dist.values
        @test Y isa LowRankParityDistribution
        @test Y isa VlasovMethods.DistributionFunction{T, 1, 1}
        @test eltype(Y) == T && size(Y) == size(f)
        @test eltype(Y.f0) == T
        @test ranks(Y) == rank
        @test (Y.plus.nfixed, Y.minus.nfixed) == (2, 1)
        for o in (Y.plus, Y.minus)
            @test eltype(o.X) == eltype(o.S) == eltype(o.V) == T
            @test orthonormality(o) ≤ 50 * eps(T)
            @test istriu(o.S)
        end
        @test fixed_mode_error(Y) ≤ 50 * eps(T)
        g = GridDistribution(Y)
        @test eltype(g.values) == T
        if rank == (16, 16)
            @test maximum(abs, g.values - f) ≤ 100 * eps(T) * maximum(abs, f)
        else
            for (a, b) in zip(velocity_moments(g), velocity_moments(dist))
                @test maximum(abs, a - b) ≤ 200 * eps(T) * maximum(abs, b)
            end
        end
    end
    # a state of high rank in both parts, which the ranks (4, 3) do not reproduce
    (; dist, Y) = fixture(T; rank = (4, 3),
        f₀ = (x, v) -> maxwellian(v - T(0.5) * sin(x / 2)) * (1 + T(0.1) * cos(x)))
    g = GridDistribution(Y)
    @test maximum(abs, g.values - dist.values) > 1e-6 * maximum(abs, dist.values)
    for (a, b) in zip(velocity_moments(g), velocity_moments(dist))
        @test maximum(abs, a - b) ≤ 200 * eps(T) * maximum(abs, b)
    end
    # the columns after the fixed modes are the leading right singular vectors of the part outside
    # them, so the weighted error of each part is the tail of the singular values of that part
    hx, hv = step(dist.x), step(dist.v)
    for (o, fpart, r) in zip((Y.plus, Y.minus), parity_split(dist.values, Y.f0), (4, 3))
        m = o.nfixed
        Vm = o.V[:, 1:m]
        σ = svdvals(sqrt(hx * hv) .* (fpart - (hv .* (fpart * Vm)) * Vm'))
        err = sqrt(hx * hv) * norm(fpart - o.X * o.S * o.V')
        @test norm(σ[(r - m + 1):end]) > 1e-5 * σ[1]
        @test abs(err - norm(σ[(r - m + 1):end])) ≤
              10 * eps(T) * sqrt(hx * hv) * norm(fpart)
    end
end

@testset "the parity of the bases, bit for bit, $T, nv = $nv" for T in (Float32, Float64),
    nv in (32, 33)

    for rank in ((16, 16), (4, 3))
        (; Y) = fixture(T; nv, rank)
        @test even_bits(Y.plus.V)
        @test odd_bits(Y.minus.V)
        isodd(nv) && @test all(iszero, Y.minus.V[cld(nv, 2), :])
    end
end

@testset "construction: refused inputs and the default background, $T" for T in (Float32, Float64)
    (; dist) = fixture(T)
    v = collect(dist.v)
    M = maxwellian.(v)
    nx, nv = size(dist.values)
    vals = dist.values
    @test_throws r"^ArgumentError: the low-rank parity distribution needs a velocity grid symmetric" LowRankParityDistribution(
        GridDistribution(copy(vals), (zero(T), T(4π)), (T(-5), T(6))), (4, 3))
    @test_throws r"^ArgumentError: the background f0 of the low-rank parity distribution is not even" LowRankParityDistribution(
        dist, (4, 3); background = M .* (1 .+ T(0.1) .* v))
    @test_throws r"^ArgumentError: the background f0 of the parity splitting is a function of v alone.*matrix" LowRankParityDistribution(
        dist, (4, 3); background = zeros(T, nx, nv))
    @test_throws r"^DimensionMismatch: the background of the parity splitting takes nv = \d+ values; got \d+$" LowRankParityDistribution(
        dist, (4, 3); background = M[1:(nv - 1)])
    finite = r"^ArgumentError: the values and the background of the low-rank parity distribution must be finite$"
    @test_throws finite LowRankParityDistribution(dist, (4, 3);
        background = [T(NaN); M[2:end]])
    for bad in (T(NaN), T(Inf))
        w = copy(vals)
        w[3, 5] = bad
        @test_throws finite LowRankParityDistribution(
            GridDistribution(dist.x, dist.v, w), (4, 3))
    end
    for rank in ((1, 3), (4, 0), (17, 3), (4, 17))
        @test_throws r"^ArgumentError: the rank -?\d+ of f[⁺⁻] is not in" LowRankParityDistribution(
            dist, rank)
    end

    # where cld(nv, 2) and fld(nv, 2) are below nx, they bound the ranks
    for (nv, bad) in ((24, ((13, 3), (4, 13))), (25, ((14, 3), (4, 13))))
        d = fixture(T; nv).dist
        for rank in bad
            @test_throws r"^ArgumentError: the rank 1[34] of f[⁺⁻] is not in" LowRankParityDistribution(
                d, rank)
        end
        @test ranks(LowRankParityDistribution(d, (cld(nv, 2), fld(nv, 2)))) ==
              (cld(nv, 2), fld(nv, 2))
    end

    # the default background: the x-mean of the even part, not of the values
    Y = LowRankParityDistribution(dist, (4, 3))
    meaneven = vec(sum((vals .+ vals[:, end:-1:1]) ./ 2; dims = 1)) ./ nx
    meanvals = vec(sum(vals; dims = 1)) ./ nx
    @test eltype(Y.f0) == T
    @test maximum(abs, Y.f0 - meaneven) ≤ 8 * eps(T) * maximum(abs, meaneven)
    @test maximum(abs, Y.f0 - meanvals) > 1e-3 * maximum(abs, meanvals)
    # a given background within the evenness bound is symmetrised
    b = M .+ 2 .* eps(T) .* maximum(M) .* (1:nv) ./ nv
    @test b != b[end:-1:1]
    Yb = LowRankParityDistribution(dist, (4, 3); background = b)
    @test Yb.f0 == Yb.f0[end:-1:1]
    @test maximum(abs, Yb.f0 - M) ≤ 4 * eps(T) * maximum(M)
    # the evenness bound takes the largest absolute value, so a negative background is accepted
    @test LowRankParityDistribution(dist, (4, 3); background = -M).f0 == -M
end

@testset "the method object and its guards" begin
    m = LowRankParitySplitting((6, 4))
    @test m isa LowRankParitySplitting
    @test (m.rank, m.minrank, m.tol) == ((6, 4), (6, 4), (0, 0))
    m = LowRankParitySplitting((8, 8); tol = (1e-6, 1e-5), minrank = (2, 1))
    @test (m.rank, m.minrank, m.tol) == ((8, 8), (2, 1), (1e-6, 1e-5))
    @test_throws r"^ArgumentError: LowRankParitySplitting takes a rank" LowRankParitySplitting((
        1, 1))
    @test_throws r"^ArgumentError: LowRankParitySplitting takes a rank" LowRankParitySplitting((
        2, 0))
    @test_throws r"^ArgumentError: LowRankParitySplitting takes 2 ≤ minrank" LowRankParitySplitting(
        (4, 4); minrank = (
            5, 1))
    @test_throws r"^ArgumentError: LowRankParitySplitting takes 2 ≤ minrank" LowRankParitySplitting(
        (4, 4); minrank = (
            1, 1))
    @test_throws r"^ArgumentError: LowRankParitySplitting takes 2 ≤ minrank" LowRankParitySplitting(
        (4, 4); minrank = (
            2, 0))
    @test_throws r"^ArgumentError: LowRankParitySplitting takes 2 ≤ minrank" LowRankParitySplitting(
        (4, 4); minrank = (
            2, 5))
    for bad in (-1e-3, NaN, Inf)
        @test_throws r"^ArgumentError: LowRankParitySplitting takes finite tolerances" LowRankParitySplitting(
            (4, 4); tol = (bad, 0))
        @test_throws r"^ArgumentError: LowRankParitySplitting takes finite tolerances" LowRankParitySplitting(
            (4, 4); tol = (0, bad))
    end
    for T in (Float32, Float64)
        (; Y, field) = fixture(T)
        before = factors(Y)
        for rank in ((17, 4), (4, 17))
            @test_throws r"^ArgumentError: the rank 17 of f[⁺⁻] is not in" lowrank_step!(
                Y, LowRankParitySplitting(rank; minrank = (2, 1)), field, 0, 0.05)
            @test factors(Y) == before
        end
        # on 16 × 24 the even and the odd subspace of 12 dimensions bound the ranks
        (; Y, field) = fixture(T; nv = 24)
        before = factors(Y)
        for rank in ((13, 4), (4, 13))
            @test_throws r"^ArgumentError: the rank 13 of f[⁺⁻] is not in" lowrank_step!(
                Y, LowRankParitySplitting(rank; minrank = (2, 1)), field, 0, 0.05)
            @test factors(Y) == before
        end
    end
end

@testset "one step equals the dense transcription, $T" for T in (Float32, Float64)
    bound = T == Float64 ? 1e-12 : 1e-4
    for (method, rmax, rmin, tol) in (
        (LowRankParitySplitting((6, 6)), (6, 6), (6, 6), (zero(T), zero(T))),
        (LowRankParitySplitting((8, 8); tol = (1e-6, 1e-6), minrank = (2, 1)),
        (8, 8), (2, 1), (T(1e-6), T(1e-6))),
    # a tolerance above the norm leaves each part at its own least rank
        (LowRankParitySplitting((8, 8); tol = (10, 10), minrank = (5, 3)),
        (8, 8), (5, 3), (T(10), T(10))))
        (; Y, field, G) = fixture(T)
        f₀ = values_of(Y)
        d = dense_parity_case(Y.plus.space, G, Y.f0)
        op, om = dense_objects(Y)
        op, om = dense_parity_step(d, op, om, T(0.05); tol, rmax, rmin)
        @test lowrank_step!(Y, method, field, 0, 0.05) === Y
        @test rel(Y.plus.X * Y.plus.S * Y.plus.V', dense_object(op)) ≤ bound
        @test rel(Y.minus.X * Y.minus.S * Y.minus.V', dense_object(om)) ≤ bound
        @test ranks(Y) == (size(op.S, 1), size(om.S, 1))
        @test maximum(abs, values_of(Y) - f₀) > 1000 * eps(T) * maximum(abs, f₀)
        for o in (Y.plus, Y.minus)
            @test eltype(o.X) == eltype(o.S) == eltype(o.V) == T
        end
    end
    # A state of high rank whose density and drift are sums of eight Fourier modes in x with
    # seeded coefficients, so that no span of a few x factors holds the others.
    rng = Xoshiro(7)
    a, b, φ = T.(0.05 .* randn(rng, 8)), T.(0.2 .* randn(rng, 8)), T.(2π .* rand(rng, 8))
    profile(c, x) = sum(c[k] * cos(k * x / 2 + φ[k]) for k in 1:8)
    generic(x, v) = (1 + profile(a, x)) * maxwellian(v - profile(b, x))
    # On 64 × 128 at rank (4, 4) the closed space has 17 of 64 columns, so the K-step reaches the
    # result; on 16 × 256 at rank (6, 6) the odd basis of the f⁻ substep has 90 of 128 columns,
    # so the velocity range of the increment does. The oracle state and the state of high rank
    # are stepped on both grids.
    oracle(x, v) = (1 + T(0.2) * cos(x / 2)) * maxwellian(v) * (1 + T(0.3) * v)
    for (nx, nv, rank) in ((64, 128, (4, 4)), (16, 256, (6, 6))), f₀ in (oracle, generic)

        (; Y, field, G) = fixture(T; nx, nv, rank, f₀)
        d = dense_parity_case(Y.plus.space, G, Y.f0)
        op, om = dense_objects(Y)
        op, om = dense_parity_step(d, op, om, T(0.05);
            tol = (zero(T), zero(T)), rmax = rank, rmin = rank)
        lowrank_step!(Y, LowRankParitySplitting(rank), field, 0, 0.05)
        @test rel(Y.plus.X * Y.plus.S * Y.plus.V', dense_object(op)) ≤ bound
        @test rel(Y.minus.X * Y.minus.S * Y.minus.V', dense_object(om)) ≤ bound
        @test ranks(Y) == (size(op.S, 1), size(om.S, 1)) == rank
        for o in (Y.plus, Y.minus)
            @test eltype(o.X) == eltype(o.S) == eltype(o.V) == T
        end
    end
    # The state of high rank on 16 × 64, with other ranks, tolerances and least ranks per part. A
    # truncation to a fixed rank drops the directions of a smaller trial space to second order in
    # τ, and tolerances of 1e-10 keep them; in Float32 a tolerance below eps(T) sets the rank by
    # round-off, so it takes 1e-5.
    ε = T == Float64 ? 1e-10 : 1e-5
    for (method, rmax, rmin, tol) in (
        (LowRankParitySplitting((4, 4)), (4, 4), (4, 4), (zero(T), zero(T))),
        (LowRankParitySplitting((6, 5); tol = (1e-4, 1e-3), minrank = (3, 2)),
        (6, 5), (3, 2), (T(1e-4), T(1e-3))),
        (LowRankParitySplitting((12, 14); tol = (ε, ε), minrank = (3, 2)),
        (12, 14), (3, 2), (T(ε), T(ε))))
        (; Y, field, G) = fixture(T; nv = 64, rank = (4, 4), f₀ = generic)
        d = dense_parity_case(Y.plus.space, G, Y.f0)
        op, om = dense_objects(Y)
        op, om = dense_parity_step(d, op, om, T(0.05); tol, rmax, rmin)
        lowrank_step!(Y, method, field, 0, 0.05)
        @test rel(Y.plus.X * Y.plus.S * Y.plus.V', dense_object(op)) ≤ bound
        @test rel(Y.minus.X * Y.minus.S * Y.minus.V', dense_object(om)) ≤ bound
        @test ranks(Y) == (size(op.S, 1), size(om.S, 1))
    end
end

@testset "twenty steps keep the fixed modes and the parity, $T, nv = $nv" for T in (
        Float32, Float64),
    nv in (32, 33)

    (; Y, field) = fixture(T; nv)
    method = LowRankParitySplitting((6, 6))
    for n in 1:20
        lowrank_step!(Y, method, field, 0.05 * (n - 1), 0.05)
        @test (Y.plus.nfixed, Y.minus.nfixed) == (2, 1)
        @test fixed_mode_error(Y) ≤ 50 * eps(T)
        @test even_bits(Y.plus.V)
        @test odd_bits(Y.minus.V)
        @test ranks(Y) == (6, 6)
    end
end

# an even orthonormal basis of 12 columns whose first two span {1, v²}, and the coefficients
# K = g hv V of a seeded random even g
function truncation_case(T)
    (; dist) = fixture(T)
    space = GridFactorSpace(dist)
    nx, nv = size(space)
    hv = step(space.v)
    v = collect(space.v)
    rng = Xoshiro(3)
    even(A) = (A .+ A[end:-1:1, :]) ./ 2
    g = permutedims(even(permutedims(T.(randn(rng, nx, nv)))))
    V, _ = orthonormalize_v(space, hcat(one.(v), v .^ 2, even(T.(randn(rng, nv, 10)))))
    return space, V, g * (hv .* V)
end

@testset "the truncation keeps the moments of the fixed modes, $T" for T in (Float32, Float64)
    space, V, K = truncation_case(T)
    hx, hv = step(space.x), step(space.v)
    X, S, V′ = VlasovMethods._truncate_keeping_moments(space, K, V, 2;
        tol = zero(T), minrank = 4, rank = 4)
    @test size(S) == (4, 4) && size(X, 2) == 4 && size(V′, 2) == 4
    @test eltype(X) == eltype(S) == eltype(V′) == T
    @test istriu(S)
    @test V′[:, 1:2] == V[:, 1:2]
    @test norm(hx .* (X' * X) - I) ≤ 50 * eps(T) && norm(hv .* (V′' * V′) - I) ≤ 50 * eps(T)
    em = rel(X * S * V′' * (hv .* V[:, 1:2]), K[:, 1:2])
    @test em ≤ 64 * eps(T)
    # the control: the weighted truncated SVD of the whole of K
    F = svd(sqrt(hx) .* K)
    Kc = (F.U[:, 1:4] * (F.S[1:4] .* F.Vt[1:4, :])) ./ sqrt(hx)
    @test rel(Kc * (V' * (hv .* V[:, 1:2])), K[:, 1:2]) ≥ 1e-3
end

@testset "the truncation's tolerance is relative, and its rank bounds hold, $T" for T in (
    Float32, Float64)
    space, V, K = truncation_case(T)
    hx = step(space.x)
    trunc(K; kw...) = VlasovMethods._truncate_keeping_moments(space, K, V, 2; kw...)
    σ = svdvals(sqrt(hx) .* K[:, 3:end])
    nrm = sqrt(hx) * norm(K)
    # 0.05 keeps every column, 0.3 drops some
    for tol in (T(0.05), T(0.3))
        _, S, _ = trunc(K; tol, minrank = 2, rank = 12)
        _, Ss, _ = trunc(T(1e-6) .* K; tol, minrank = 2, rank = 12)
        r = size(S, 1)
        @test size(Ss, 1) == r
        @test r > 2
        tol == T(0.3) && @test r < 12
        @test (r == 12 ? zero(T) : norm(σ[(r - 1):end])) ≤ tol * nrm
        @test norm(σ[(r - 2):end]) > tol * nrm
    end
    @test size(trunc(K; tol = T(10), minrank = 3, rank = 12)[2], 1) == 3
    @test size(trunc(K; tol = T(10), minrank = 2, rank = 12)[2], 1) == 2
    @test size(trunc(K; tol = zero(T), minrank = 2, rank = 4)[2], 1) == 4
    # with tol = 0 the tail of exact zeros is dropped: σ = (σ₁, 0, …, 0) keeps one column
    Kz = hcat(K[:, 1:3], zeros(T, size(K, 1), 9))
    @test size(trunc(Kz; tol = zero(T), minrank = 2, rank = 12)[2], 1) == 3
end

@testset "the closed space holds X̂, 1 and G X̂, $T" for T in (Float32, Float64)
    for (nx, nv, k, ncols) in ((64, 128, 8, 17), (16, 32, 12, 16))
        (; dist, G) = fixture(T; nx, nv)
        space = GridFactorSpace(dist)
        hx = step(space.x)
        X̂, _ = orthonormalize_x(space, T.(randn(Xoshiro(5), nx, k)))
        Xa = VlasovMethods._closed_space(space, X̂, G)
        @test eltype(Xa) == T
        @test size(Xa) == (nx, ncols)
        @test norm(hx .* (Xa' * Xa) - I) ≤ 50 * eps(T)
        for A in (X̂, ones(T, nx, 1), G * X̂)
            @test norm(A - Xa * (hx .* (Xa' * A))) ≤ 64 * eps(T) * norm(A)
        end
    end
end

@testset "an odd number of velocity nodes and the least ranks: energy over 400 steps, $T" for T in (
    Float32, Float64)
    bound = T == Float64 ? 1e-12 : 8 * 3 * 10^4 * eps(Float32)
    for (nv, rank) in ((33, (6, 6)), (32, (2, 1)))
        (; dist, Y, field, G) = fixture(T; nv, β = 0, rank)
        method = LowRankParitySplitting(rank)
        H₀ = energy(dist, G, values_of(Y))
        drift = zero(T)
        for n in 1:400
            lowrank_step!(Y, method, field, 0.05 * (n - 1), 0.05)
            drift = max(drift, abs(energy(dist, G, values_of(Y)) - H₀) / abs(H₀))
        end
        @test eltype(Y.plus.S) == T
        @test ranks(Y) == rank
        @test drift ≤ bound
        @test maximum(abs, values_of(Y) - dist.values) > 1e-6 * maximum(abs, dist.values)
    end
end

@testset "full rank reproduces ParitySplitting over twenty steps, $T" for T in (Float32, Float64)
    (; dist, model, Y, field) = fixture(T; β = 0, rank = (16, 16))
    nx, nv = size(dist.values)
    method = LowRankParitySplitting((16, 16))
    for n in 1:20
        lowrank_step!(Y, method, field, 0.05 * (n - 1), 0.05)
    end
    z = mktempdir() do dir
        file = joinpath(dir, "run.h5")
        run!(ParitySplitting(model, (0, 1), 0.05), file)
        h5read(file, "z")
    end
    @test size(z, 2) == 21
    f = reshape(z[:, 21], nx, nv)
    @test eltype(f) == T
    @test norm(values_of(Y) - f) / norm(f) ≤ (T == Float64 ? 1e-10 : 1e-4)
    @test norm(f - reshape(z[:, 1], nx, nv)) > 1e-4 * norm(f)
end

@testset "the decided edges: a non-finite step, $T" for T in (Float32, Float64)
    (; Y, field) = fixture(T)
    method = LowRankParitySplitting((6, 6))
    before = factors(Y)
    pattern = r"^the step τ is not finite in the element type Float(32|64) of the state; got τ = (NaN|Inf|-Inf|1\.0e39)(32)?$"
    τs = T == Float32 ? (NaN, Inf, -Inf, 1e39) : (NaN, Inf, -Inf)
    for τ in τs
        e = try
            lowrank_step!(Y, method, field, 0, τ)
            nothing
        catch err
            err
        end
        @test e isa ArgumentError
        @test e isa ArgumentError && occursin(pattern, e.msg)
        @test factors(Y) == before
    end
end

@testset "the decided edges: a field of the bounded closure, $T" for T in (Float32, Float64)
    (; Y, G) = fixture(T)
    before = factors(Y)
    bounded = LowRankFields1D1V(Y.plus.space, G; closure = BoundedVelocityClosure())
    @test_throws r"^ArgumentError: .*periodic velocity closure" lowrank_step!(
        Y, LowRankParitySplitting((6, 6)), bounded, 0, 0.05)
    @test factors(Y) == before
end

@testset "the decided edges: the float type of τ, a zero step and a homogeneous state" begin
    (; Y, field) = fixture(Float32)
    Z = deepcopy(Y)
    method = LowRankParitySplitting((6, 6))
    lowrank_step!(Y, method, field, 0, 0.05)
    lowrank_step!(Z, method, field, 0.0f0, 0.05f0)
    @test factors(Y) == factors(Z)
    @test eltype(Y.plus.S) == Float32

    for T in (Float32, Float64)
        (; Y, field) = fixture(T)
        f₀ = values_of(Y)
        lowrank_step!(Y, LowRankParitySplitting((6, 6)), field, 0, 0)
        @test maximum(abs, values_of(Y) - f₀) ≤ 8 * eps(T) * maximum(abs, f₀)

        # τ < 0 is accepted, and a step back returns close to the start
        lowrank_step!(Y, LowRankParitySplitting((6, 6)), field, 0, 0.05)
        lowrank_step!(Y, LowRankParitySplitting((6, 6)), field, 0.05, -0.05)
        @test status(Y)
        @test rel(values_of(Y), f₀) ≤ 1e-3

        h = fixture(T; f₀ = (x, v) -> maxwellian(v))
        Yh = h.Y
        M = maximum(maxwellian.(collect(h.dist.v)))
        for n in 1:10
            lowrank_step!(Yh, LowRankParitySplitting((6, 6)), h.field, 0.1 * (n - 1), 0.1)
        end
        @test eltype(Yh.plus.S) == T
        @test maximum(abs, values_of(Yh) - h.dist.values) ≤ 64 * eps(T) * M
    end
end

@testset "the decided edges: a non-finite state, $T" for T in (Float32, Float64)
    for bad in (T(NaN), T(Inf))
        (; Y, field) = fixture(T)
        Y.plus.S[1, 1] = bad
        r = ranks(Y)
        @test lowrank_step!(Y, LowRankParitySplitting((6, 6)), field, 0, 0.05) === Y
        @test !status(Y)
        @test ranks(Y) == r
    end
end

# one warm step, then the bytes of the next, through a function barrier
function step_allocations(Y, method, field, τ)
    lowrank_step!(Y, method, field, 0, τ)
    return @allocated lowrank_step!(Y, method, field, 0, τ)
end
allocations(f::F, a::A) where {F, A} = (f(a); @allocated f(a))

@testset "no step forms an nx × nv array, $T" for T in (Float32, Float64)
    # skipped where a run forces --check-bounds=yes
    if Base.JLOptions().check_bounds != 1
        method = LowRankParitySplitting((4, 4))
        bytes = map(((256, 512), (512, 1024))) do (nx, nv)
            (; Y, field) = fixture(T; nx, nv, rank = (4, 4))
            step_allocations(Y, method, field, 0.05)
        end
        @test bytes[2] / bytes[1] ≤ 3
        # Bytes of the form a nx + b nv + c have a zero mixed difference over two nx and two nv,
        # and an nx × nv array adds nx nv sizeof(T) to it, 128 × 256 × sizeof(T) here at least.
        mixed = map(((128, 256), (256, 256), (128, 512), (256, 512))) do (nx, nv)
            (; Y, field) = fixture(T; nx, nv, rank = (4, 4))
            step_allocations(Y, method, field, 0.05)
        end
        @test abs(mixed[4] - mixed[3] - mixed[2] + mixed[1]) ≤ 128 * 256 * sizeof(T) / 10
        # the control: the barrier sees an allocation
        @test allocations(collect, 1:4) > 0
    end
end
