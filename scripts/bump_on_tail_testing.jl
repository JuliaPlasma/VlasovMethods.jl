# Bump-on-tail, step 3 of 3: the full-order and the reduced model on the test sweep.
#
# Reads the bases that `bump_on_tail_projections.jl` wrote with the same options, draws the ten
# test wave numbers uniform in [0.1, 0.5] after the seed, and marches each with the full-order
# leapfrog and with two reduced leapfrogs of the particle basis: one with the DEIM field, as the
# reference marches it, and one with the full field of the lifted particles. The three marches of
# a sample start from the particles of the training sweep and take the same field. The defaults
# are the paper's parameters; the paper command is the default:
#
#     julia --project=scripts scripts/bump_on_tail_testing.jl
#
# Writes into `results_path`: the growth rates `growth_rates.csv`, with the header
# `kappa,gamma_hfm,gamma_rom,gamma_rom_full` and one row per test sample, and the figures of the
# paper's items 24 to 29, each named `item<NN>_<what>.png`; a figure of the full-field model alone
# ends in `_full`.

include("options.jl")
include("bump_on_tail_common.jl")

using CairoMakie
using ReducedBasisMethods: get_DEIM_indices

const OPT = parse_arguments(ARGS; BUMP_ON_TAIL..., samples = 10)

"The axis of a figure of `plot_series`, on a log scale."
function log_scale!(fig)
    only(c for c in fig.content if c isa Axis).yscale = log10
    return fig
end

"""
    by_model!(fig, models, n)

The figure `fig` of `plot_series` over `n` unlabelled series per model, in the order of the names
`models`, with one colour and one legend entry per model, as the reference draws all samples at
once.
"""
function by_model!(fig, models, n)
    ax = only(c for c in fig.content if c isa Axis)
    colours = [(c, 0.5) for c in (:blue, :orange, :green)][eachindex(models)]
    for (i, line) in enumerate(filter(p -> p isa Lines, ax.scene.plots))
        line.color = colours[cld(i, n)]
    end
    axislegend(ax, [LineElement(color = c) for c in colours], collect(models))
    return fig
end

function main(o)
    Φ, Ψₑ = h5open(runs_path(run_id(o) * "_projections.h5"), "r") do io
        ReducedBasis(io).Φ, ReducedBasis(io, "field").Φ
    end
    deim = get_DEIM_indices(Ψₑ)

    particles = initial_particles(o)
    poisson = potential(o)
    κ, _ = test_sweep(o.seed; n = o.samples)
    pspace = parameter_space(κ)

    nt = step_count(o.final_time, o.dt)
    nh = length(PoissonSolvers.coefficients(poisson))
    IP = IntegratorParameters(o.dt, nt, nt + 1, nh, o.particles, length(pspace))
    HFM = Snapshots(Float64, 1, o.particles, nh, nt + 1, pspace)
    # The two reduced models of the reference: the DEIM field, which its testing script marches,
    # and the full field of the lifted particles, which it keeps commented out beside it. Each
    # file of a model carries its suffix; the DEIM model's files have none.
    ROM = (DEIM = Snapshots(Float64, 1, o.particles, nh, nt + 1, pspace),
        full = Snapshots(Float64, 1, o.particles, nh, nt + 1, pspace))
    name = (DEIM = "ROM, DEIM", full = "ROM, full field")
    suffix = (DEIM = "", full = "_full")
    IC = ReducedIntegratorCache(IP, size(Φ, 2))
    X₀ = Φ' * vec(particles.x)

    for p in 1:length(pspace)
        params = pspace(p)
        println("test sample $p of $(length(pspace)), κ = $(κ[p])")
        integrate_vp(particles, params, vp_field(poisson, params), HFM, IP, p)
        reduced = (
            DEIM = DEIMElectricField(vp_field(poisson, params), Φ, Ψₑ, deim, X₀),
            full = ReducedElectricField(vp_field(poisson, params), Φ, X₀))
        for m in keys(ROM)
            reduced_integrate_vp(particles, Φ, params, reduced[m], ROM[m], IP, IC, p)
        end
    end

    t = IP.t
    n = length(κ)
    samples = ["κ = $(round(k; digits = 3))" for k in κ]
    models = ["HFM"; collect(name)]
    columns(M) = collect(eachcol(M))
    γ_hfm = fit_growth_rate(t, HFM.W)
    γ_rom = map(R -> fit_growth_rate(t, R.W), ROM)

    # The errors from the first save on, on log axes as in the reference: at t = 0 they are those
    # of the projection of the initial particles alone, and a log axis cannot draw a zero.
    later = 2:length(t)
    at = saved_indices(IP.nₛ)

    for m in keys(ROM)
        R = ROM[m]

        # item 24: the trajectory errors against time, and against κ at T, 0.75T and 0.5T
        errors = trajectory_errors(HFM.X[1, :, :, :], HFM.V[1, :, :, :],
            R.X[1, :, :, :], R.V[1, :, :, :])
        for q in keys(errors)
            save(results_path("item24_$(q)_time$(suffix[m]).png"),
                log_scale!(plot_series(
                    t[later], columns(errors[q][later, :]); labels = samples,
                    ylabel = "$(q), $(name[m])")))
            save(results_path("item24_$(q)_kappa$(suffix[m]).png"),
                log_scale!(plot_series(κ, [errors[q][i, :] for i in at];
                    labels = ["t = $(t[i])" for i in at], xlabel = "κ",
                    ylabel = "$(q), $(name[m])")))
        end

        # item 25: the absolute error of the potential against time
        save(results_path("item25_potential_error$(suffix[m]).png"),
            log_scale!(plot_series(
                t[later], columns(potential_error(HFM.Φ, R.Φ)[later, :]);
                labels = samples, ylabel = "‖Φ_ROM − Φ_HFM‖, $(name[m])")))
    end

    # item 26: the field energy of the three models, of every sample, and of the samples whose
    # growth rates of the HFM and the DEIM model agree best and worst
    save(results_path("item26_field_energy.png"),
        by_model!(
            log_scale!(plot_series(t, vcat(columns(HFM.W), (columns(R.W) for R in ROM)...);
                ylabel = "W")), models, n))
    # a sample whose energy has fewer than two maxima has no rate, and no rank
    finite = findall(p -> isfinite(γ_hfm[p]) && isfinite(γ_rom.DEIM[p]), 1:n)
    ranked = isempty(finite) ? Int[] :
             finite[[best_and_worst(γ_hfm[finite], γ_rom.DEIM[finite])...]]
    for (case, p) in zip(("best", "worst"), ranked)
        save(results_path("item26_field_energy_$(case).png"),
            log_scale!(plot_series(t, [HFM.W[:, p], (R.W[:, p] for R in ROM)...];
                labels = models .* ", " .* samples[p], ylabel = "W")))
    end

    # item 27: the energy and its relative drift of the three models, and the size of the drift
    # against κ at T, 0.75T and 0.5T
    H = hcat(HFM.K .+ HFM.W, (R.K .+ R.W for R in ROM)...)
    drift = reduce(hcat, energy_drift.(columns(H)))
    save(results_path("item27_energy.png"),
        by_model!(plot_series(t, columns(H); ylabel = "H"), models, n))
    save(results_path("item27_energy_drift.png"),
        by_model!(plot_series(t, columns(drift); ylabel = "ΔH/H"), models, n))
    save(results_path("item27_energy_drift_kappa.png"),
        log_scale!(plot_series(κ,
            [abs.(drift[i, ((j - 1) * n + 1):(j * n)]) for j in eachindex(models)
             for i in at];
            labels = ["$(model), t = $(t[i])" for model in models for i in at],
            xlabel = "κ", ylabel = "|ΔH/H|")))

    # item 28: the growth rate against κ, with the fit of the reference
    save(results_path("item28_growth_rate.png"),
        plot_series(κ, (γ_hfm, γ_rom...); labels = models, xlabel = "κ", ylabel = "γ"))
    open(results_path("growth_rates.csv"), "w") do io
        println(io, "kappa,gamma_hfm,gamma_rom,gamma_rom_full")
        for p in eachindex(κ)
            println(io, κ[p], ",", γ_hfm[p], ",", γ_rom.DEIM[p], ",", γ_rom.full[p])
        end
    end

    # item 29: the particles of the sample nearest the reference wave number, at t = 0 and at T
    p = argmin(abs.(κ .- KAPPA_REFERENCE))
    domains = (xdomain = (0.0, 2π / KAPPA_REFERENCE), vdomain = (-6.0, 8.0))
    save(results_path("item29_particles_t0.png"), plot_particles(HFM, 1, p; domains...))
    save(results_path("item29_particles_T_hfm.png"),
        plot_particles(HFM, IP.nₛ, p; domains...))
    for m in keys(ROM)
        save(results_path("item29_particles_T_rom$(suffix[m]).png"),
            plot_particles(ROM[m], IP.nₛ, p; domains...))
    end

    println("k = $(size(Φ, 2)), kₑ = $(size(Ψₑ, 2))")
    for p in eachindex(κ)
        println("κ = $(κ[p]): γ_HFM = $(γ_hfm[p]), γ_ROM = $(γ_rom.DEIM[p]), ",
            "γ_ROM full field = $(γ_rom.full[p])")
    end
    return nothing
end

main(OPT)
