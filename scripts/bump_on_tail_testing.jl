# The test sweep of the paper's bump-on-tail example: the ten test wave numbers of the reference's
# `bump_on_tail_3_testing.jl`, the full-order and the reduced march of each, their errors, their
# energies and their growth rates, and the figures of items 24 to 29.
#
#   julia --project=scripts scripts/bump_on_tail_testing.jl [--np 10000] [--nh 16] [--order 3] \
#       [--dt 0.05] [--final-time 25.0] [--seed 1234] [--runs-dir dir] [--results-dir dir]
#
# The default step is `0.05` and not the paper command's `0.1`: the growth-rate fit of item 28 reads
# the first four maxima of the field energy over a window of five saved states, and at `0.1` the
# ninth of the ten drawn time scales resolves so coarsely that its fit turns negative. The paper's
# command:
#
#   julia --project=scripts scripts/bump_on_tail_testing.jl --dt 0.1 --final-time 25.0
#
# The reduced model takes the particle basis that `bump_on_tail_projections.jl` wrote into the runs
# directory, so run that script first: this one finds it as the runs directory's single
# `*_projections.h5`. The growth rates go to `growth_rates.csv`, the figures to the results
# directory.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie

const OPT = parse_arguments(ARGS; np = BOT_NP, nh = BOT_NH, order = BOT_ORDER, dt = 0.05,
    final_time = BOT_T, seed = BOT_SEED)

const PARAMS = (
    np = OPT.np, nh = OPT.nh, order = OPT.order, dt = OPT.dt, T = OPT.final_time,
    κlo = KAPPA_RANGE[1], κhi = KAPPA_RANGE[2], seed = OPT.seed)

# the particle basis of the training run, and nothing else of it: the field of the reduced model is
# the scaled Poisson field, evaluated at every lifted particle, so it needs no field basis
ppath = joinpath(RUNS_DIR[],
    first(filter(f -> endswith(f, "_projections.h5"), readdir(RUNS_DIR[]))))
particle = ReducedBasis(ppath, "particle")
Ψ = particle.Φ

sweep = test_sweep()
nₛ = step_count(PARAMS.T, PARAMS.dt) + 1
nparam = length(sweep.κ)
t = collect(range(0, PARAMS.T; length = nₛ))
space(χ) = ParameterSpace(Parameter(:χ, [χ]))

# the diagnostics of every sample, one column per test wave number
x2 = zeros(nₛ, nparam);
v2 = zeros(nₛ, nparam)
xmax = zeros(nₛ, nparam);
vmax = zeros(nₛ, nparam)
ϕerr = zeros(nₛ, nparam)
Wh = zeros(nₛ, nparam);
Wr = zeros(nₛ, nparam)
Hh = zeros(nₛ, nparam);
Hr = zeros(nₛ, nparam)
dHh = zeros(nₛ, nparam);
dHr = zeros(nₛ, nparam)
γh = zeros(nparam);
γr = zeros(nparam)
xdomain = (0.0, 2π / BOT_KAPPA)
vdomain = (-10.0, 10.0)

@printf("bump_on_tail_testing: %d test samples, %d particles, k = %d, dt = %g, t ∈ [0, %g]\n",
    nparam, PARAMS.np, particle.k, PARAMS.dt, PARAMS.T)

for p in 1:nparam
    χ = sweep.χ[p]
    sample = test_sample(PARAMS, χ)
    HFM = sample.snapshots
    plt = sample.particles.particles
    ip = sample.ip

    RSS = Snapshots(Float64, 1, PARAMS.np, PARAMS.nh, nₛ, space(χ))
    rfield = ReducedElectricField(sample.field, Ψ, Ψ' * vec(plt.x))
    reduced_integrate_vp(plt, Ψ, (χ = χ,), rfield, RSS, ip,
        ReducedIntegratorCache(ip, particle.k), 1)

    err = trajectory_errors(@view(HFM.X[1, :, :, 1]), @view(HFM.V[1, :, :, 1]),
        @view(RSS.X[1, :, :, 1]), @view(RSS.V[1, :, :, 1]))
    x2[:, p] .= vec(err.x2)
    v2[:, p] .= vec(err.v2)
    xmax[:, p] .= vec(err.xmax)
    vmax[:, p] .= vec(err.vmax)

    ϕerr[:, p] .= potential_error(@view(HFM.Φ[:, :, 1]), @view(RSS.Φ[:, :, 1]))

    # a march's `W` is the energy of the scaled field, `W/χ²`, and the paper's `W` is the plain
    # one: the reference multiplies the scaled energy back by `χ²` before it plots or adds
    Wh[:, p] .= HFM.W[:, 1] .* χ^2
    Wr[:, p] .= RSS.W[:, 1] .* χ^2
    # the energy of the model is the kinetic energy of the particles plus the field energy
    Hh[:, p] .= Wh[:, p] .+ HFM.K[:, 1]
    Hr[:, p] .= Wr[:, p] .+ RSS.K[:, 1]
    dHh[:, p] .= energy_drift(@view(Hh[:, p]))
    dHr[:, p] .= energy_drift(@view(Hr[:, p]))

    γh[p] = fit_growth_rate(ip.t, @view(Wh[:, p]))
    γr[p] = fit_growth_rate(ip.t, @view(Wr[:, p]))

    @printf("  κ = %.4f, χ = %.4f: γ_HFM = %.4f, γ_ROM = %.4f\n", sweep.κ[p], χ, γh[p],
        γr[p])

    # item 29: the particles of the first sample at the initial and the final saved state
    if p == 1
        for (what, ts) in (("t0", 1), ("T", nₛ))
            figure = plot_particles(HFM, ts, 1; xdomain, vdomain)
            fpath = results_path(item_figure(29, "particles_$what"))
            save(fpath, figure)
            @printf("  wrote %s\n", fpath)
        end
    end
end

# item 28 and the growth-rate table: γ of every test sample, the two models in one axis
csv = results_path("growth_rates.csv")
open(csv, "w") do io
    println(io, "kappa,gamma_hfm,gamma_rom")
    for p in 1:nparam
        @printf(io, "%.6f,%.6f,%.6f\n", sweep.κ[p], γh[p], γr[p])
    end
end
@printf("  wrote %s\n", csv)

figure = plot_series(sweep.κ, [γh, γr]; labels = ["HFM", "ROM"], xlabel = "κ", ylabel = "γ")
fpath = results_path(item_figure(28, "growth_rate"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

# item 24: the four components of `trajectory_errors` — the root mean squares `x2` and `v2`, and
# the `L∞` maxima `xmax` and `vmax` — against time, and against κ at three saved times
figure = plot_series(t, [x2[:, p] for p in 1:nparam]; xlabel = "t", ylabel = "Δx₂")
fpath = results_path(item_figure(24, "position_error"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

figure = plot_series(t, [v2[:, p] for p in 1:nparam]; xlabel = "t", ylabel = "Δv₂")
fpath = results_path(item_figure(24, "velocity_error"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

figure = plot_series(t, [xmax[:, p] for p in 1:nparam]; xlabel = "t", ylabel = "Δmax x")
fpath = results_path(item_figure(24, "position_max"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

figure = plot_series(t, [vmax[:, p] for p in 1:nparam]; xlabel = "t", ylabel = "Δmax v")
fpath = results_path(item_figure(24, "velocity_max"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

iT, imid, ihalf = nₛ, div(3 * (nₛ - 1), 4) + 1, div(nₛ - 1, 2) + 1
figure = plot_series(sweep.κ,
    [x2[iT, :], v2[iT, :], x2[imid, :], v2[imid, :], x2[ihalf, :], v2[ihalf, :]];
    labels = ["Δx₂(T)", "Δv₂(T)", "Δx₂(0.75T)", "Δv₂(0.75T)", "Δx₂(0.5T)", "Δv₂(0.5T)"],
    xlabel = "κ", ylabel = "error")
fpath = results_path(item_figure(24, "errors_vs_kappa"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

figure = plot_series(sweep.κ,
    [xmax[iT, :], vmax[iT, :], xmax[imid, :],
        vmax[imid, :], xmax[ihalf, :], vmax[ihalf, :]];
    labels = ["Δmax x(T)", "Δmax v(T)", "Δmax x(0.75T)", "Δmax v(0.75T)",
        "Δmax x(0.5T)", "Δmax v(0.5T)"],
    xlabel = "κ", ylabel = "error")
fpath = results_path(item_figure(24, "max_errors_vs_kappa"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

# item 25: the absolute potential error against time
figure = plot_series(t, [ϕerr[:, p] for p in 1:nparam], xlabel = "t",
    ylabel = "‖Φ_ROM − Φ_HFM‖")
fpath = results_path(item_figure(25, "potential_error"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

# item 26: the field energy, full order against reduced
figure = plot_series(t, vcat([Wh[:, p] for p in 1:nparam], [Wr[:, p] for p in 1:nparam]);
    xlabel = "t", ylabel = "W")
fpath = results_path(item_figure(26, "field_energy"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

# item 27: the energy and its relative drift, full order against reduced
figure = plot_series(t, vcat([Hh[:, p] for p in 1:nparam], [Hr[:, p] for p in 1:nparam]);
    xlabel = "t", ylabel = "H")
fpath = results_path(item_figure(27, "energy"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)

figure = plot_series(t, vcat([dHh[:, p] for p in 1:nparam], [dHr[:, p] for p in 1:nparam]);
    xlabel = "t", ylabel = "ΔH/H")
fpath = results_path(item_figure(27, "energy_drift"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)
