# The linear Landau damping rate of the grid 1D1V Vlasov–Poisson scheme, and its stability under
# refinement. Run from the repository root:
#
#   julia --project=scripts scripts/landau_damping_1d1v.jl
#
# The rate is fitted to the peaks of log|E₁(t)| and checked against the textbook value
# γ ≈ −0.1533 for k = 0.5 and unit thermal velocity, to 1%. The second run halves h_v and dt, and
# its rate must agree with the first to 1.5e-4: that is the reason the 1% tolerance is defensible,
# and the run is what shows it. The stepper is RK4 — the implicit midpoint's dense Newton Jacobian
# at n_x·n_v = 64·2049 would need 137 GB.

using VlasovMethods
using PoissonSolvers
using GeometricIntegrators
using Printf

const K = 0.5
const ALPHA = 0.01
const NX = 64
const L = 4π
const VDOMAIN = (-8.6, 8.6)
const TSPAN = (0.0, 40.0)

# |E₁(t)| = k·(2/nx)·|Σᵢ ϕᵢ e^{−ikxᵢ}|, with ϕ the field of the density ρᵢ = h_v Σⱼ f_ij. The
# field energy decays at 2γ, so the fit is taken on the field's own mode and not on its energy.
function mode1(q, solver, nv, hv, x)
    ρ = hv .* vec(sum(reshape(q, NX, nv); dims = 2))
    ϕ = solve(solver, ρ)

    return K * (2 / NX) * abs(sum(ϕ .* exp.(-im * K .* x)))
end

# The local maxima of `y = log|E₁|`, each refined by the parabola through the three samples around
# it, with the refined time inside [8, 40]. A sample is a peak when `y[n-1] < y[n] ≥ y[n+1]`, so a
# plateau of two equal samples gives one peak and a minimum never passes; the vertex of a parabola
# through such a triple lies within half a step of `t[n]`, and `y[n-1] - 2y[n] + y[n+1] < 0`.
function peaks(times, y)
    ts = Float64[]
    ys = Float64[]
    for n in 2:(length(times) - 1)
        y₋, y₀, y₊ = y[n - 1], y[n], y[n + 1]
        y₋ < y₀ ≥ y₊ || continue
        denom = y₋ - 2y₀ + y₊
        δ = (y₋ - y₊) / (2 * denom)
        t = times[n] + δ * (times[n + 1] - times[n])
        8 ≤ t ≤ 40 || continue
        push!(ts, t)
        push!(ys, y₀ - (y₊ - y₋)^2 / (8 * denom))
    end

    return ts, ys
end

# The least-squares slope through the peaks; NaN for fewer than two, which the guard fails closed on.
function slope(ts, ys)
    length(ts) ≥ 2 || return NaN
    t̄ = sum(ts) / length(ts)
    ȳ = sum(ys) / length(ys)

    return sum((ts .- t̄) .* (ys .- ȳ)) / sum((ts .- t̄) .^ 2)
end

function landau_rate(nv, dt)
    basis = FFTWBasis((0.0, L), NX)
    x = basis.xgrid[1:NX]
    v = range(VDOMAIN[1], VDOMAIN[2]; length = nv)

    values = zeros(NX, nv)
    values .= (1 .+ ALPHA .* cos.(K .* x)) .* exp.(-v' .^ 2 ./ 2) ./ sqrt(2π)

    dist = GridDistribution(values, (0.0, L), VDOMAIN)
    model = VlasovPoisson(dist, Potential(basis))
    hv = step(dist.v)
    solver = model.potential.solver

    int = VlasovMethods.GeometricIntegrator(model, TSPAN, dt, RK4())
    prob = int.equation
    t₀ = GeometricIntegrators.GeometricEquations.initialtime(prob)
    Δt = GeometricIntegrators.GeometricEquations.timestep(prob)
    n = GeometricIntegrators.GeometricEquations.ntime(prob)
    solstep = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(
        int.integrator, GeometricIntegrators.GeometricEquations.initialstate(prob))

    # one solution step at a time, keeping only |E₁|: the state history of the second run is
    # 0.84 GB, and nothing else reads it
    times = Vector{Float64}(undef, n + 1)
    logE = Vector{Float64}(undef, n + 1)
    times[1] = t₀
    logE[1] = log(mode1(solstep.q, solver, nv, hv, x))
    for s in 1:n
        GeometricIntegrators.reset!(solstep, t₀ + s * Δt)
        GeometricIntegrators.integrate!(solstep, int.integrator)
        times[s + 1] = t₀ + s * Δt
        logE[s + 1] = log(mode1(solstep.q, solver, nv, hv, x))
    end

    ts, ys = peaks(times, logE)

    return slope(ts, ys), length(ts)
end

function main()
    rates = Float64[]
    for (nv, dt) in ((1025, 0.05), (2049, 0.025))
        γ, npeaks = landau_rate(nv, dt)
        @printf("nv = %d, dt = %.3f: γ = %.6f over %d peaks\n", nv, dt, γ, npeaks)
        push!(rates, γ)
    end

    γ₁, γ₂ = rates
    (abs(γ₁ + 0.1533) ≤ 1.5e-3 && abs(γ₁ - γ₂) < 1.5e-4) || error(
        "the fitted rate γ = $(γ₁) is not within 1.5e-3 of the textbook −0.1533, or halving h_v " *
        "and dt moved it by $(abs(γ₁ - γ₂)) ≥ 1.5e-4 (γ = $(γ₂))")

    return nothing
end

main()
