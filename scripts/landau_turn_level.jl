# What sets the level at which ln E₁ of the grid Landau run stops decaying and grows.
#
# Linear Landau damping, k = 0.5, f₀ = (1 + 0.01 cos(x/2)) M(v) on [0, 4π) × [−vmax, vmax], with the
# full-order model `VlasovPoisson(GridDistribution, Potential(FFTWBasis))` (Arakawa bracket) and
# GeometricIntegrators' RK4, in Float64. Each run records log10 E₁ at every step and log10 of the
# largest mode-1 amplitude |f̂₁| on the outermost unknown velocity nodes (j = 2 and nᵥ − 1). From
# the local maxima ("peaks") of log10 E₁ it reports:
#
#   - turn: the lowest peak after t = 5, provided a later peak lies 0.3 decades or more above it;
#     otherwise "none" up to the end time;
#   - rate: the least-squares slope of ln E₁ through the peaks in (t_turn, t_turn + 15];
#   - landau: the same slope through the peaks in [8, 30], the Landau rate (≈ −0.1534).
#
# The group `linear` builds the linearised right-hand side about the Maxwellian, restricted to the
# x-mode k (a real 2nᵥ × 2nᵥ matrix; the field is quadratic in f, so a central difference gives its
# linear part exactly). It reports the largest |Re λ|, the spectral radius ρ over all x-modes
# (RK4 is stable for τ ≤ 2.83/ρ), the spacing of the lowest frequencies over k hᵥ/√3, and the
# turn of the exact linear evolution exp(At) u₀, which has neither a time step nor a nonlinearity.
#
# The finding. Every |Re λ| is round-off, and the exact linear evolution turns at the level of the
# RK4 run, so the turn is no instability, no time-step effect and no nonlinearity. On the mode k,
# Arakawa's bracket with h = v²/2 streams f with the v-average
# (sin(k hₓ)/hₓ) (v_{j−1/2} f_{j−1} + 4 v_j f_j + v_{j+1/2} f_{j+1}) / 6. On a filament whose phase
# advances by ξ per velocity node the streaming speed is v (2 + cos ξ)/3, so ξ grows at
# k hᵥ (2 + cos ξ)/3 (the frequency spacing k hᵥ/√3 near ω = 0, and a recurrence time
# √3 · 2π/(k hᵥ)), and the filament moves outward in v at k hᵥ v sin ξ / 3: the perturbation is
# stretched in v by 3/(2 + cos ξ), up to 3. The stretching carries the perturbation of the bulk
# to ±vmax, and E₁ turns where the field that the edge of the grid returns overtakes the Landau
# decay. The turn moves later and lower with vmax at fixed hᵥ, and it comes with either velocity
# closure; nₓ and τ do not move it.
#
# Run from the repository root:
#
#   julia --startup-file=no --project=. scripts/landau_turn_level.jl [outdir] [group ...]
#
# with groups base, closure, vmax, nv, nx, tau, linear (default: all of them). All groups take
# 100 s on one thread of a laptop at load average 124: `closure` 40 s, `linear` 20 s, `vmax`
# 10 s, the rest under 10 s each, plus about 10 s of compilation per process. It
# writes, into `outdir` (default `<tempdir>/landau_turn_level`), one CSV per run
# (`t,log10E1,log10edge`) and appends to `summary.tsv`, one line per run, and prints the same lines.
using LinearAlgebra
using Printf
using GeometricIntegrators
using VlasovMethods

const PS = VlasovMethods.PoissonSolvers
const GE = GeometricIntegrators.GeometricEquations
const K = 0.5
const L = 4π
const ALPHA = 0.01

function landau_model(nx, nv, vmax)
    basis = PS.FFTWBasis((0.0, L), nx)
    x = basis.xgrid[1:nx]
    v = range(-vmax, vmax; length = nv)
    values = (1 .+ ALPHA .* cos.(K .* x)) .* exp.(-v' .^ 2 ./ 2) ./ sqrt(2π)
    dist = GridDistribution(values, (0.0, L), (-vmax, vmax))
    return VlasovPoisson(dist, PS.Potential(basis)), x
end

# the complex mode-1 amplitudes f̂₁ⱼ = (2/nₓ) Σᵢ fᵢⱼ e^{−ikxᵢ}, and |E₁| = |hᵥ Σⱼ f̂₁ⱼ| / k
mode1(f, x) = vec((2 / length(x)) .* (transpose(cis.(-K .* x)) * f))
E1(f̂, hv) = abs(hv * sum(f̂)) / K

function run_case(; nx = 64, nv = 128, vmax = 6.0, closure = :periodic, tau = 0.05,
        tend = 100.0)
    model, x = landau_model(nx, nv, vmax)
    vc = closure === :bounded ? BoundedVelocityClosure() : PeriodicVelocityClosure()
    int = VlasovMethods.GeometricIntegrator(model, (0.0, tend), tau, RK4(); closure = vc)
    prob = int.equation
    n = GE.ntime(prob)
    solstep = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(int.integrator,
        GE.initialstate(prob))
    hv = step(model.distribution.v)
    t = collect(0:n) .* tau
    y = zeros(n + 1)
    edge = zeros(n + 1)
    record(s, q) = (f̂ = mode1(reshape(q, nx, nv), x);
        y[s] = log10(E1(f̂, hv)); edge[s] = log10(max(abs(f̂[2]), abs(f̂[nv - 1]))))
    record(1, solstep.q)
    for s in 1:n
        GeometricIntegrators.reset!(solstep, s * tau)
        GeometricIntegrators.integrate!(solstep, int.integrator)
        record(s + 1, solstep.q)
    end
    return t, y, edge, hv
end

peaks(y) = [s for s in 2:(length(y) - 1) if y[s - 1] < y[s] ≥ y[s + 1]]
function slope(t, y)
    (tm = sum(t) / length(t);
        sum((t .- tm) .* (y .- sum(y) / length(y))) / sum(abs2, t .- tm))
end

# turn time, turn level (log10 E₁), ln-rate after the turn, Landau ln-rate; NaN where undefined
function analyse(t, y)
    p = filter(s -> t[s] > 5, peaks(y))
    landau = slope(t[filter(s -> 8 ≤ t[s] ≤ 30, p)], log(10) .*
                                                     y[filter(s -> 8 ≤ t[s] ≤ 30, p)])
    i = argmin(y[p])
    turned = any(y[p[(i + 1):end]] .≥ y[p[i]] + 0.3)
    turned || return (NaN, NaN, NaN, landau)
    after = filter(s -> t[p[i]] < t[s] ≤ t[p[i]] + 15, p)
    rate = length(after) ≥ 3 ? slope(t[after], log(10) .* y[after]) : NaN
    return (t[p[i]], y[p[i]], rate, landau)
end

const HEADER = "group\tnx\tnv\tvmax\thv\tclosure\ttau\tt_turn\tlog10E1_turn\tlog10edge_turn\t" *
               "rate_after\tlandau_rate\tseconds"

function case!(io, outdir, group; kw...)
    c = merge(
        (nx = 64, nv = 128, vmax = 6.0, closure = :periodic, tau = 0.05, tend = 100.0),
        kw)
    s0 = time()
    t, y, edge, hv = run_case(; c...)
    secs = time() - s0
    # a run whose E₁ stops being finite is recorded as unstable, with the first such time
    bad = findfirst(!isfinite, y)
    if bad !== nothing
        line = @sprintf("%s\t%d\t%d\t%g\t%.5f\t%s\t%g\tunstable: E1 not finite from t = %.2f",
            group, c.nx, c.nv, c.vmax, hv, c.closure, c.tau, t[bad])
        println(io, line)
        println(line)
        return nothing
    end
    tt, lt, rate, landau = analyse(t, y)
    et = isnan(tt) ? NaN : edge[findfirst(==(tt), t)]
    name = @sprintf("%s_nx%d_nv%d_vmax%g_%s_tau%g", group, c.nx, c.nv, c.vmax, c.closure,
        c.tau)
    open(joinpath(outdir, name * ".csv"), "w") do f
        println(f, "t,log10E1,log10edge")
        foreach(s -> @printf(f, "%.4f,%.10g,%.10g\n", t[s], y[s], edge[s]), eachindex(t))
    end
    line = @sprintf("%s\t%d\t%d\t%g\t%.5f\t%s\t%g\t%.2f\t%.3f\t%.3f\t%.4f\t%.4f\t%.1f",
        group,
        c.nx, c.nv, c.vmax, hv, c.closure, c.tau, tt, lt, et, rate, landau, secs)
    println(io, line)
    println(line)
    flush(io)
    flush(stdout)
end

# the linearised field about the x-independent part `fbar` of the initial state, on the x-mode
# m k: f = Σⱼ (cⱼ cos mkx + sⱼ sin mkx) eⱼ is u = [c; s]
function linear_operator(rhs!, params, fbar, x, m)
    nx, nv = size(fbar)
    C, S = cos.(m * K .* x), sin.(m * K .* x)
    A = zeros(2nv, 2nv)
    dp, dm = zeros(nx * nv), zeros(nx * nv)
    for col in 1:(2nv)
        δ = zeros(nx, nv)
        δ[:, mod1(col, nv)] .= col ≤ nv ? C : S
        rhs!(dp, 0.0, vec(fbar .+ δ), params)
        rhs!(dm, 0.0, vec(fbar .- δ), params)
        Lδ = reshape((dp .- dm) ./ 2, nx, nv)
        A[:, col] .= (2 / nx) .* [transpose(Lδ) * C; transpose(Lδ) * S]
    end
    return A
end

# The largest |Re λ| on the mode k, the spectral radius ρ over the modes m k, m = 1 … nₓ/2 − 1,
# and log10 E₁ of the exact linear evolution exp(At) u₀ on the mode k.
function linear_case(; nx = 64, nv = 128, vmax = 6.0, closure = :periodic, tend = 100.0)
    model, x = landau_model(nx, nv, vmax)
    vc = closure === :bounded ? BoundedVelocityClosure() : PeriodicVelocityClosure()
    prob = VlasovMethods.GeometricIntegrator(model, (0.0, 1.0), 0.05, RK4(); closure = vc).equation
    rhs!, params = GE.functions(prob).v, GE.parameters(prob)
    q0 = reshape(copy(GE.initialstate(prob).q), nx, nv)
    fbar = repeat(sum(q0; dims = 1) ./ nx, nx, 1)
    ρ = maximum(m -> maximum(abs, eigvals(linear_operator(rhs!, params, fbar, x, m))),
        1:(nx ÷ 2 - 1))
    F = eigen(linear_operator(rhs!, params, fbar, x, 1))
    u0 = (2 / nx) .* [transpose(q0) * cos.(K .* x); transpose(q0) * sin.(K .* x)]
    a = F.vectors \ u0
    hv = step(model.distribution.v)
    t = collect(0:0.05:tend)
    y = map(t) do τ
        u = real(F.vectors * (a .* exp.(F.values .* τ)))
        log10(abs(hv * sum(u[1:nv] .- im .* u[(nv + 1):end])) / K)
    end
    # the spacing of the two lowest distinct positive frequencies, against k hᵥ / √3, the spacing
    # of the averaged streaming v (2 + cos ξ)/3 (see the header)
    ω = unique(x -> round(x; digits = 8), sort(filter(>(1e-10), imag.(F.values))))
    return t, y, maximum(abs ∘ real, F.values), ρ, (ω[2] - ω[1]) / (K * hv / sqrt(3))
end

function main(args)
    outdir = isempty(args) || args[1] in GROUPS ? joinpath(tempdir(), "landau_turn_level") :
             args[1]
    groups = filter(in(GROUPS), args)
    isempty(groups) && (groups = GROUPS)
    mkpath(outdir)
    @printf("uptime: %s\n", readchomp(`uptime`))
    io = open(joinpath(outdir, "summary.tsv"), "a")
    println(io, "# ", readchomp(`uptime`))
    println(io, HEADER)
    println(HEADER)
    hv6 = 12 / 127
    for g in groups
        g == "base" && case!(io, outdir, g)
        # The bounded closure's end nodes see a jump of h from vmax²/2 to 0, whose x-advection
        # makes the spectral radius ρ 11-22 times the periodic one (group `linear`): RK4 needs
        # τ ≤ 2.83/ρ, 0.0093 at nᵥ = 128 on [−6, 6] and 0.0027 at nᵥ = 256 on [−10, 10].
        if g == "closure"
            case!(io, outdir, g; closure = :bounded, tau = 0.005)
            case!(io, outdir, g; nv = 256, vmax = 10.0, closure = :bounded, tau = 0.0025)
            case!(io, outdir, g; nv = 256, vmax = 10.0)
        end
        # fixed hᵥ = 12/127; the periodic ρ grows with vmax, and at vmax = 12 RK4 at τ = 0.05 is
        # unstable (recorded), so that run is repeated at τ = 0.025
        if g == "vmax"
            for vmax in (8.0, 10.0, 12.0)
                case!(io, outdir, g; vmax, nv = round(Int, 2vmax / hv6) + 1)
            end
            case!(io, outdir, g; vmax = 12.0, nv = 255, tau = 0.025)
            case!(io, outdir, g; vmax = 15.0, nv = 319, tau = 0.025, tend = 130.0)
            case!(io, outdir, g; vmax = 10.0, nv = 128)
        end
        g == "nv" && foreach(nv -> case!(io, outdir, g; nv), (192, 256, 384, 512))
        g == "nx" && foreach(nx -> case!(io, outdir, g; nx), (32, 128))
        g == "tau" && foreach(tau -> case!(io, outdir, g; tau), (0.025, 0.1))
        if g == "linear"
            for c in ((nv = 128, vmax = 6.0), (nv = 128, vmax = 6.0, closure = :bounded),
                (nv = 256, vmax = 10.0), (nv = 256, vmax = 10.0, closure = :bounded))
                s0 = time()
                t, y, reλ, ρ, dω = linear_case(; c...)
                tt, lt, rate, landau = analyse(t, y)
                line = @sprintf("linear\t64\t%d\t%g\t%.5f\t%s\texact\t%.2f\t%.3f\tmax|Reλ|=%.1e,ρ=%.1f,dω/(khv/√3)=%.3f\t%.4f\t%.4f\t%.1f",
                    c.nv, c.vmax, 2c.vmax / (c.nv - 1), get(c, :closure, :periodic), tt, lt,
                    reλ, ρ, dω, rate, landau, time() - s0)
                println(io, line)
                println(line)
            end
        end
    end
    close(io)
end

const GROUPS = ["base", "closure", "vmax", "nv", "nx", "tau", "linear"]

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
