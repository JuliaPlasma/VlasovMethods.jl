# The runs of the dispersion-rate tests of the three BUG integrators at the rank 16, against the
# roots of the linear dispersion relation at k = 0.5 in chapter 38 of the RCM book: the linear
# Landau root ±1.415661889 − 0.1533594669 i, so the damping rate of log|E₁| is −0.1533594669, and
# the two-stream growth rate 0.2592499317. The rates are fitted as VlasovReductionBenchmark fits
# them (`rate_fits.jl`). Three test files share these runs, so that CI spreads them over its
# workers: `test/verification/lowrank_bug_landau_rates.jl` and the two-stream files
# `lowrank_bug_two_stream_rates.jl` (Float64) and `lowrank_bug_two_stream_rates_float32.jl`.
#
# The grid is 64 × 256 nodes of x ∈ [0, 4π), v ∈ [−10, 10]. On it the full-order model itself is
# off by 1.26 % (Landau) and 3.28 % (two-stream). f at v = ±10 reaches 1e-9 to 1e-7 of its largest
# value over the runs, not round-off: the velocity resolution of this grid spreads f in v.

using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V

include("rate_fits.jl")

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

# log|E₁| at each time of t₀:τ:tend of the method from the factorisation
# `LowRankDistribution(dist, r)` at the rank 16 of (1 + α cos(x/2)) g(v)
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

# the two-stream growth rate of `method` in the precision T, within 5 % of the root
function two_stream_rate_test(T, method)
    t, y = log_e1(T, method, two_stream, 1e-3, 25.0)
    rate = growth_rate_fit(t, y, (5.0, 25.0))
    @test abs(rate - TWO_STREAM_RATE) / TWO_STREAM_RATE ≤ 0.05
end
