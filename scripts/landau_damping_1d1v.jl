# The linear Landau damping rate of the grid 1D1V Vlasov–Poisson scheme, and its stability under
# refinement. Run from the repository root:
#
#   julia --project=scripts scripts/landau_damping_1d1v.jl [--runs-dir dir] [--results-dir dir]
#
# The rate is fitted to the peaks of log|E₁(t)| and checked against the textbook value
# γ ≈ −0.1533 for k = 0.5 and unit thermal velocity, to 1%. The second run halves h_v and dt, and
# its rate must agree with the first to 1.5e-4: that is the reason the 1% tolerance is defensible,
# and the run is what shows it. The stepper is RK4 — the implicit midpoint's dense Newton Jacobian
# at n_x·n_v = 64·2049 would need 137 GB.
#
# Each run's own series goes to `--runs-dir`, the figure with both fitted rates to
# `--results-dir`, and the script fails when the two conditions above do not hold.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using GeometricIntegrators
using PoissonSolvers
using VlasovMethods

const PARSED = let (positional, options) = parse_arguments(ARGS)
    no_positional(positional, "landau_damping_1d1v.jl")
    options
end

const K = 0.5
const ALPHA = 0.01
const NX = 64
const L = 4π
const VDOMAIN = (-8.6, 8.6)
const TSPAN = (0.0, 40.0)

# The two runs of the refinement pair, and what the fit has to satisfy.
const RUNS = ((1025, 0.05), (2049, 0.025))
const TEXTBOOK = -0.1533

# |E₁(t)| = k·(2/nx)·|Σᵢ ϕᵢ e^{−ikxᵢ}|, with ϕ the field of the density ρᵢ = h_v Σⱼ f_ij. The
# field energy decays at 2γ, so the fit is taken on the field's own mode and not on its energy.
function mode1(q, solver, nv, hv, x)
    ρ = hv .* vec(sum(reshape(q, NX, nv); dims = 2))
    ϕ = solve(solver, ρ)

    return K * (2 / NX) * abs(sum(ϕᵢ * exp(-im * K * xᵢ) for (ϕᵢ, xᵢ) in zip(ϕ, x)))
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

    return sum((t - t̄) * (y - ȳ) for (t, y) in zip(ts, ys)) /
           sum(abs2(t - t̄) for t in ts)
end

# The whole run: the fitted rate, the peaks it was fitted through, and the series behind them, so
# that the caller can write both out.
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

    return slope(ts, ys), ts, ys, times, logE
end

function main()
    rates = Float64[]
    series = []

    for (nv, dt) in RUNS
        γ, ts, ys, times, logE = landau_rate(nv, dt)
        @printf("nv = %d, dt = %.3f: γ = %.6f over %d peaks\n", nv, dt, γ, length(ts))
        push!(rates, γ)
        push!(series, (
            nv = nv, dt = dt, γ = γ, ts = ts, ys = ys, times = times, logE = logE))

        # every fifth sample: the series is 801 and 1601 points, and the peaks are what the fit
        # reads, so the sampled series is written alongside them
        numbers = runs_path(@sprintf("landau_damping_1d1v_nv%d.txt", nv))
        open(numbers, "w") do io
            @printf(io, "# nv = %d, dt = %.3f, gamma = %.6f, %d peaks\n", nv, dt, γ,
                length(ts))
            println(io, "# t log|E1|")
            for (t, y) in zip(times, logE)
                @printf(io, "%.6f %.15e\n", t, y)
            end
            println(io, "# peak t, peak log|E1|")
            for (t, y) in zip(ts, ys)
                @printf(io, "%.6f %.15e\n", t, y)
            end
        end
        @printf("  wrote %s\n", numbers)
    end

    γ₁, γ₂ = rates
    ok = abs(γ₁ - TEXTBOOK) ≤ 1.5e-3 && abs(γ₁ - γ₂) < 1.5e-4

    fig = Figure(size = (1200, 500))
    ax = Axis(fig[1, 1], xlabel = "t", ylabel = "log|E₁|",
        title = @sprintf("k = %.2f, α = %.2f, nx = %d", K, ALPHA, NX))
    colours = (:steelblue, :firebrick)
    for (s, colour) in zip(series, colours)
        lines!(ax, s.times, s.logE, color = colour,
            label = @sprintf("nv = %d, dt = %.3f: γ = %.6f", s.nv, s.dt, s.γ))
        scatter!(ax, s.ts, s.ys, color = colour, markersize = 8)
    end
    axislegend(ax; position = :lb)
    # the textbook slope, laid through the first peak of the first run: the fit's own intercept is
    # what the data decides, so the comparison is the slope and not the height
    t₀, y₀ = series[1].ts[1], series[1].ys[1]
    lines!(ax, [8.0, 40.0], y₀ .+ TEXTBOOK .* ([8.0, 40.0] .- t₀), color = :black,
        linestyle = :dash, label = @sprintf("γ = %.4f", TEXTBOOK))

    figure = results_path("landau_damping_1d1v.png")
    save(figure, fig)
    @printf("  wrote %s\n", figure)

    ok || error(
        "the fitted rate γ = $(γ₁) is not within 1.5e-3 of the textbook $(TEXTBOOK), or " *
        "halving h_v and dt moved it by $(abs(γ₁ - γ₂)) ≥ 1.5e-4 (γ = $(γ₂))")

    return nothing
end

main()
