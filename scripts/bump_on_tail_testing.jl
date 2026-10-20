# The test sweep of the bump-on-tail example: ten wave numbers drawn uniform in `[0.1, 0.5]`, each
# marched by the full-order model and by the reduced model with the DEIM field of
# `bump_on_tail_projections.jl`, from the initial condition of `bump_on_tail_training.jl`, and the
# paper's diagnostics of the reduced model against the full one.
#
#   julia --project=scripts scripts/bump_on_tail_testing.jl [--runid bump_on_tail] \
#       [--kappa-min 0.1] [--kappa-max 0.5] [--nsamples 10] [--seed 1234] \
#       [--runs-dir dir] [--results-dir dir]
#
# The defaults are the paper's parameters; the particles, the grid and the steps are the training
# run's, read from its file. The run reads `<runid>_training.h5` and `<runid>_projections.h5` from
# `--runs-dir`, writes `<runid>_testing.h5` there, and writes to `--results-dir` the figures of the
# position and velocity errors against time and against `κ` (item 24), the potential error (25),
# the field energy (26), the energy and its drift (27), the growth rates (28) and three
# phase-space frames (29).

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using PoissonSolvers
using ReducedBasisMethods: ReducedBasis, h5save
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods

const OPT = parse_arguments(ARGS; runid = BUMP_ON_TAIL_DEFAULTS.runid,
    kappa_min = BUMP_ON_TAIL_DEFAULTS.kappa_min, kappa_max = BUMP_ON_TAIL_DEFAULTS.kappa_max,
    nsamples = BUMP_ON_TAIL_DEFAULTS.nsamples, seed = BUMP_ON_TAIL_DEFAULTS.seed)

# the training run: its steps, its initial condition and its grid
IP₀, x, v, w, ncells, degree = h5open(runs_path(OPT.runid * "_training.h5"), "r") do file
    IntegratorParameters(file, "integrator"), read(file["x"]), read(file["v"]),
    read(file["w"]), read(HDF5.attributes(file)["ncells"]),
    read(HDF5.attributes(file)["degree"])
end
P₀ = VlasovMethods.ParticleList(x, v, w)

# the bases
pb, fb, deim = h5open(runs_path(OPT.runid * "_projections.h5"), "r") do file
    ReducedBasis(file, "particle_basis"), ReducedBasis(file, "field_basis"),
    read(file["deim_indices"])
end
Φ = pb.Φ
X₀ = Φ' * vec(P₀.x)

potential = Potential(PeriodicBasisSpline(domain(KAPPA_REF), degree, ncells))

# the test sweep
κs = testing_kappas(OPT.kappa_min, OPT.kappa_max, OPT.nsamples, OPT.seed)
pspace = ParameterSpace(Parameter(:χ, [march_parameters(κ).χ for κ in κs]))
IP = IntegratorParameters(IP₀.dt, IP₀.nₜ, IP₀.nₛ, IP₀.nₕ, IP₀.nₚ, length(κs))
SSh = Snapshots(Float64, 1, IP.nₚ, IP.nₕ, IP.nₛ, pspace)
SSr = Snapshots(Float64, 1, IP.nₚ, IP.nₕ, IP.nₛ, pspace)
RIC = ReducedIntegratorCache(IP, pb.k)

for p in eachindex(κs)
    params = march_parameters(κs[p])
    @printf("bump_on_tail_testing: sample %d of %d, κ = %.4f, χ = %.4f\n",
        p, length(κs), κs[p], params.χ)
    integrate_vp(P₀, params, vp_field(potential, params), SSh, IP, p)
    rfield = DEIMElectricField(vp_field(potential, params), Φ, fb.Φ, deim, X₀)
    reduced_integrate_vp(P₀, Φ, params, rfield, SSr, IP, RIC, p)
end

fpath = runs_path(OPT.runid * "_testing.h5")
h5open(fpath, "w") do file
    h5save(file, SSh; path = "full")
    h5save(file, SSr; path = "reduced")
    file["kappa"] = κs
    HDF5.attributes(file)["commit"] = repository_commit()
    HDF5.attributes(file)["command"] = command_line("bump_on_tail_testing.jl", ARGS)
    HDF5.attributes(file)["seed"] = OPT.seed
end
@printf("  wrote %s\n", fpath)

t = IP.t
labels = [@sprintf("κ = %.3f", κ) for κ in κs]

function write_figure(fig, name; logscale = false)
    logscale && (fig.content[1].yscale = log10)
    figure = results_path(OPT.runid * "_" * name * ".png")
    save(figure, fig)
    @printf("  wrote %s\n", figure)
end

# item 24: the errors of the positions and the velocities, against time for each sample, and
# against `κ` at the final, the three-quarter and the half time
e = trajectory_errors(SSh.X[1, :, :, :], SSh.V[1, :, :, :], SSr.X[1, :, :, :], SSr.V[
    1, :, :, :])
for (q, name) in ((:x, "positions"), (:v, "velocities"))
    rms, largest = e[Symbol(q, 2)], e[Symbol(q, :max)]
    write_figure(
        plot_series(t, collect(eachcol(rms)); labels = labels,
            ylabel = "RMS error of the $(name)"),
        "errors_$(q)_time";
        logscale = true)

    s = saved_indices(IP.nₛ)
    series = [vec(E[i, :]) for E in (rms, largest) for i in s]
    names = ["$(n) at t = $(t[i])" for n in ("RMS", "max") for i in s]
    write_figure(
        plot_series(κs, series; labels = names, xlabel = "κ",
            ylabel = "error of the $(name)"),
        "errors_$(q)_parameter";
        logscale = true)
end

# item 25: the error of the potential
write_figure(
    plot_series(t, collect(eachcol(potential_error(SSh.Φ, SSr.Φ))); labels = labels,
        ylabel = "‖Φ_ROM − Φ_HFM‖"),
    "potential_error";
    logscale = true)

# item 28: the growth rates, and the samples whose rates agree best and worst
γh = growth_rate(t, SSh.W)
γr = growth_rate(t, SSr.W)
best, worst = best_and_worst(γh, γr)
write_figure(
    plot_series(κs, (γh, γr); labels = ["HFM", "ROM"], xlabel = "κ", ylabel = "γ"),
    "growth_rate")
@printf("  best κ = %.4f, worst κ = %.4f\n", κs[best], κs[worst])

# items 26 and 27: the field energy, the energy and its drift of the best and the worst sample
cases = ((best, "best"), (worst, "worst"))
case_labels = ["$(m), $(c), κ = $(@sprintf("%.3f", κs[p]))"
               for (p, c) in cases for m in ("HFM", "ROM")]
write_figure(
    plot_series(t, [S.W[:, p] for (p, _) in cases for S in (SSh, SSr)];
        labels = case_labels, ylabel = "W"),
    "field_energy";
    logscale = true)
H = [S.K[:, p] .+ S.W[:, p] for (p, _) in cases for S in (SSh, SSr)]
write_figure(plot_series(t, H; labels = case_labels, ylabel = "H"), "energy")
write_figure(plot_series(t, energy_drift.(H); labels = case_labels, ylabel = "ΔH/H"),
    "energy_drift")

# item 29: the particles of the worst sample at `t = 0`, and at the final time in both models
frames = ((SSh, 1, "particles_initial"), (SSh, IP.nₛ, "particles_final_hfm"),
    (SSr, IP.nₛ, "particles_final_rom"))
for (S, ts, name) in frames
    write_figure(
        plot_particles(S, ts, worst; xdomain = domain(KAPPA_REF),
            vdomain = (-6.0, 9.0)), name)
end
