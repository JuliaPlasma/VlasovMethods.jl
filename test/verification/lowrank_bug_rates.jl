# The dispersion rates of the three BUG integrators at the rank 16, against the roots of the
# linear dispersion relation at k = 0.5 in chapter 38 of the RCM book: the linear Landau root
# ±1.415661889 − 0.1533594669 i, so the damping rate of log|E₁| is −0.1533594669, and the
# two-stream growth rate 0.2592499317. The rates are fitted as VlasovReductionBenchmark fits them
# (`test/helpers/rate_fits.jl`).
#
# The grid is 64 × 256 nodes of x ∈ [0, 4π), v ∈ [−10, 10]. On it the full-order model itself is
# off by 1.26 % (Landau) and 3.28 % (two-stream). f at v = ±10 reaches 1e-9 to 1e-7 of its largest
# value over the runs, not round-off: the velocity resolution of this grid spreads f in v.
#
# The Landau rates run in Float64 only: in Float32 the maxima fit of the damped log|E₁| is not
# monotone in τ and does not track Float64. The growing two-stream mode stays far above
# round-off, and its rates run in both precisions. The fixed-rank BUG has no Landau check: its
# Landau rate at the rank 16 is not monotone in τ either.

using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V

include("../helpers/rate_fits.jl")

const LANDAU_RATE = -0.1533594669
const TWO_STREAM_RATE = 0.2592499317

maxwellian(v) = exp(-v^2 / 2) / sqrt(2 * oftype(v, π))
two_stream(v) = v^2 * maxwellian(v)

# |E₁| = k (2/nx) |Σᵢ ϕᵢ e^{−i k xᵢ}| of the potential ϕ = G ρ of the factors, k = 1/2
function e1(Y, G, x)
    ρ = step(Y.space.v) .* (Y.X * (Y.S * vec(sum(Y.V; dims = 1))))
    ϕ = G * ρ
    return abs(sum(ϕ[i] * exp(-im * x[i] / 2) for i in eachindex(x))) / length(x)
end

# log|E₁| at each time of t₀:τ:tend of the method from L1's factorisation at the rank 16 of
# (1 + α cos(x/2)) g(v)
function log_e1(T, method, profile, α, tend; τ = 0.05, r = 16)
    xd, vd = (zero(T), T(4π)), (T(-10), T(10))
    dist = GridDistribution(zeros(T, 64, 256), xd, vd)
    dist.values .= (1 .+ T(α) .* cos.(dist.x ./ 2)) .* profile.(dist.v)'
    G = VlasovPoisson(dist, Potential(FFTWBasis(xd, 64))).work
    Y = LowRankDistribution(dist, r)
    w = LowRankFields1D1V(Y, G)
    n = round(Int, tend / τ)
    t = [s * τ for s in 0:n]
    y = zeros(n + 1)
    y[1] = log(e1(Y, G, dist.x))
    for s in 1:n
        lowrank_step!(Y, method, w, t[s], τ)
        @test all(isfinite, Y.S)
        y[s + 1] = log(e1(Y, G, dist.x))
    end
    @test eltype(Y.S) == T
    return t, y
end

@testset "linear Landau, $(nameof(typeof(method)))" for method in (AugmentedBUG(), MidpointBUG())
    t, y = log_e1(Float64, method, maxwellian, 0.01, 40.0)
    γ = landau_rate_fit(t, y, (0.0, 40.0)).γ
    @test γ < 0
    @test abs(γ - LANDAU_RATE) / abs(LANDAU_RATE) ≤ 0.02
end

@testset "two-stream instability, $(nameof(typeof(method))), $T" for T in (Float64, Float32),
    method in (BUG(), AugmentedBUG(), MidpointBUG())

    t, y = log_e1(T, method, two_stream, 1e-3, 25.0)
    rate = growth_rate_fit(t, y, (5.0, 25.0))
    @test abs(rate - TWO_STREAM_RATE) / TWO_STREAM_RATE ≤ 0.05
end
