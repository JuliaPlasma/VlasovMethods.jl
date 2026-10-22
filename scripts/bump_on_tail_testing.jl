# Bump-on-tail, step 3 of 3: the full-order and the reduced model on the test sweep.
#
# Reads the bases that `bump_on_tail_projections.jl` wrote with the same options, draws the ten
# test wave numbers uniform in [0.1, 0.5] after the seed, and marches each with the full-order
# leapfrog and with the reduced leapfrog of the particle basis and the DEIM field. The two marches
# of a sample start from the particles of the training sweep and take the same field. The defaults
# are the paper's parameters; the paper command is the default:
#
#     julia --project=scripts scripts/bump_on_tail_testing.jl
#
# Writes into `results_path`: the growth rates `growth_rates.csv`, with the header
# `kappa,gamma_hfm,gamma_rom` and one row per test sample, and the figures of the paper's items
# 24 to 29, each named `item<NN>_<what>.png`.

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
    ROM = Snapshots(Float64, 1, o.particles, nh, nt + 1, pspace)
    IC = ReducedIntegratorCache(IP, size(Φ, 2))
    X₀ = Φ' * vec(particles.x)

    for p in 1:length(pspace)
        params = pspace(p)
        println("test sample $p of $(length(pspace)), κ = $(κ[p])")
        integrate_vp(particles, params, vp_field(poisson, params), HFM, IP, p)
        reduced = DEIMElectricField(vp_field(poisson, params), Φ, Ψₑ, deim, X₀)
        reduced_integrate_vp(particles, Φ, params, reduced, ROM, IP, IC, p)
    end

    t = IP.t
    samples = ["κ = $(round(k; digits = 3))" for k in κ]
    hfm_rom = vcat("HFM, " .* samples, "ROM, " .* samples)
    columns(M) = collect(eachcol(M))

    # item 24: the trajectory errors against time, and against κ at T, 0.75T and 0.5T
    errors = trajectory_errors(HFM.X[1, :, :, :], HFM.V[1, :, :, :],
        ROM.X[1, :, :, :], ROM.V[1, :, :, :])
    at = saved_indices(IP.nₛ)
    for q in keys(errors)
        save(results_path("item24_$(q)_time.png"),
            plot_series(t, columns(errors[q]); labels = samples, ylabel = string(q)))
        save(results_path("item24_$(q)_kappa.png"),
            plot_series(κ, [errors[q][i, :] for i in at];
                labels = ["t = $(t[i])" for i in at], xlabel = "κ", ylabel = string(q)))
    end

    # item 25: the absolute error of the potential against time
    save(results_path("item25_potential_error.png"),
        plot_series(t, columns(potential_error(HFM.Φ, ROM.Φ)); labels = samples,
            ylabel = "‖Φ_ROM − Φ_HFM‖"))

    # item 26: the field energy, HFM against ROM
    save(results_path("item26_field_energy.png"),
        log_scale!(plot_series(t, vcat(columns(HFM.W), columns(ROM.W)); labels = hfm_rom,
            ylabel = "W")))

    # item 27: the energy and its relative drift, HFM against ROM
    H = hcat(HFM.K .+ HFM.W, ROM.K .+ ROM.W)
    save(results_path("item27_energy.png"),
        plot_series(t, columns(H); labels = hfm_rom, ylabel = "H"))
    save(results_path("item27_energy_drift.png"),
        plot_series(t, energy_drift.(columns(H)); labels = hfm_rom, ylabel = "ΔH/H"))

    # item 28: the growth rate against κ, with the fit of the reference
    γ_hfm = fit_growth_rate(t, HFM.W)
    γ_rom = fit_growth_rate(t, ROM.W)
    save(results_path("item28_growth_rate.png"),
        plot_series(κ, (γ_hfm, γ_rom); labels = ["HFM", "ROM"], xlabel = "κ", ylabel = "γ"))
    open(results_path("growth_rates.csv"), "w") do io
        println(io, "kappa,gamma_hfm,gamma_rom")
        for p in eachindex(κ)
            println(io, κ[p], ",", γ_hfm[p], ",", γ_rom[p])
        end
    end

    # item 29: the particles of the sample nearest the reference wave number, at t = 0 and at T
    p = argmin(abs.(κ .- KAPPA_REFERENCE))
    domains = (xdomain = (0.0, 2π / KAPPA_REFERENCE), vdomain = (-6.0, 8.0))
    save(results_path("item29_particles_t0.png"), plot_particles(HFM, 1, p; domains...))
    save(results_path("item29_particles_T_hfm.png"),
        plot_particles(HFM, IP.nₛ, p; domains...))
    save(results_path("item29_particles_T_rom.png"),
        plot_particles(ROM, IP.nₛ, p; domains...))

    println("k = $(size(Φ, 2)), kₑ = $(size(Ψₑ, 2))")
    for p in eachindex(κ)
        println("κ = $(κ[p]): γ_HFM = $(γ_hfm[p]), γ_ROM = $(γ_rom[p])")
    end
    return nothing
end

main(OPT)
