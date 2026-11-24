using LinearAlgebra: qr
using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V, _potential!

# The Landau damping rate and the two-stream growth rate of the projector splitting at r = 16 on
# 64 × 256 nodes of [0, 4π) × [−10, 10], against the roots of the dispersion relation at k = 0.5:
# 1.415662 − 0.153359 i (Landau) and the growth rate 0.25925 (two-stream). The fits are those of
# VlasovReductionBenchmark: the least-squares line through log E₁ at its local maxima on a window
# (Landau), and through log E₁ at every step on [5, 25] (two-stream), with
# E₁ = k (2/nx) |Σᵢ ϕᵢ exp(−i k xᵢ)| of the potential ϕ = G ρ.
#
# The full-order model on this grid gives 1.26 % (Landau) and 3.28 % (two-stream) off the roots,
# so the bounds of 2 % and 5 % leave the rank and the splitting 0.74 % and 1.72 %. In Float32 the
# Landau fit ends at t = 30: after t ≈ 35 the peaks of E₁ are fewer than 30 times above its
# Float32 round-off, which then splits flat peak tops and lifts the last peaks.

const k = 0.5
const LANDAU = 0.15336
const TWO_STREAM = 0.25925

maxwellian(v) = exp(-v^2 / 2) / sqrt(oftype(v, 2π))
orth(A, h) = Matrix(qr(sqrt(h) .* A).Q)[:, 1:size(A, 2)] ./ sqrt(h)

function fourier(x, r)
    T = eltype(x)
    cols = [ones(T, length(x))]
    m = 1
    while length(cols) < r
        push!(cols, cos.(m .* x ./ 2))
        length(cols) < r && push!(cols, sin.(m .* x ./ 2))
        m += 1
    end
    return reduce(hcat, cols)
end

function problem(T; nx = 64, nv = 256, vmax = 10)
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, nx, nv), xd, (T(-vmax), T(vmax)))
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    space = GridFactorSpace(dist)
    return (space = space, w = LowRankFields1D1V(space, model.work), x = collect(dist.x),
        v = collect(dist.v), hx = step(dist.x), hv = step(dist.v))
end

# the rank-r state of f₀ = (1 + α cos(x/2)) g(v): r Fourier functions of x and vʲ M(v) of v
function state(P, α, g, r)
    T = eltype(P.x)
    f0 = (1 .+ T(α) .* cos.(P.x ./ 2)) .* g.(P.v)'
    X = orth(fourier(P.x, r), P.hx)
    V = orth(reduce(hcat, [P.v .^ j .* maxwellian.(P.v) for j in 0:(r - 1)]), P.hv)
    return LowRankDistribution(P.space, X, (P.hx * P.hv) .* (X' * f0 * V), V)
end

function logE1(Y, P)
    T = eltype(P.x)
    ϕ = _potential!(zeros(T, length(P.x)), Y.X, Y.S, Y.V, P.w)
    nx = length(P.x)
    return log(T(k) * (2 / T(nx)) * abs(sum(ϕ .* cis.(-T(k) .* P.x))))
end

function run(P, Y, method, τ, n)
    T = eltype(P.x)
    y = zeros(T, n + 1)
    y[1] = logE1(Y, P)
    t = zero(T)
    for s in 1:n
        lowrank_step!(Y, method, P.w, t, τ)
        t += τ
        y[s + 1] = logE1(Y, P)
    end
    return T.(0:n) .* τ, y
end

function slope(t, y)
    tm, ym = sum(t) / length(t), sum(y) / length(y)
    return sum((t .- tm) .* (y .- ym)) / sum(abs2, t .- tm)
end

@testset "Landau and two-stream rates of the projector splitting ($T, $c)" for T in (Float64, Float32),
    c in (:lie, :strang)

    P = problem(T)
    τ = T(0.05)
    method = ProjectorSplitting(c)

    # linear Landau damping: α = 0.01, the maxima of log E₁ on [0, 40], [0, 30] in Float32
    tend = T == Float64 ? 40 : 30
    t, y = run(P, state(P, 0.01, maxwellian, 16), method, τ, round(Int, tend / 0.05))
    @test eltype(y) === T
    n = length(y)
    peaks = [s for s in 2:(n - 1) if y[s - 1] < y[s] ≥ y[s + 1]]
    γ = slope(t[peaks], y[peaks])
    @test γ < 0
    @test abs(γ / -LANDAU - 1) ≤ 0.02

    # the two-stream instability: f₀ = (1 + 10⁻³ cos(x/2)) v² M(v), every step on [5, 25]
    t, y = run(P, state(P, 1e-3, v -> v^2 * maxwellian(v), 16), method, τ, 500)
    window = 5 .≤ t .≤ 25
    γ = slope(t[window], y[window])
    @test γ > 0
    @test abs(γ / TWO_STREAM - 1) ≤ 0.05
end
