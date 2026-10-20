# The training sweep of the self-consistent bump-on-tail example: the full-order Vlasov–Poisson
# march of one initial condition at ten wave numbers `κ ∈ [0.1, 0.5]`, each as the reference
# problem at `κ = 0.3` with the scale `χ = κ / 0.3`, and the accelerations that the field applied
# at every saved state. It is the first of three drivers; `bump_on_tail_projections.jl` builds the
# bases from its snapshots, and `bump_on_tail_testing.jl` marches the reduced model.
#
#   julia --project=scripts scripts/bump_on_tail_training.jl [--runid bump_on_tail] \
#       [--npart 10000] [--ncells 16] [--degree 3] [--tstep 0.1] [--final-time 25.0] \
#       [--kappa-min 0.1] [--kappa-max 0.5] [--nsamples 10] [--seed 1234] \
#       [--runs-dir dir] [--results-dir dir]
#
# The defaults are the paper's parameters. The run writes `<runid>_training.h5` to `--runs-dir`,
# and two figures to `--results-dir`: the field energy `W(t)` of each sample, and the growth rate
# of each sample against `κ`.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using PoissonSolvers
using ReducedBasisMethods: h5save
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods

const OPT = parse_arguments(ARGS; BUMP_ON_TAIL_DEFAULTS...)

# the initial condition, shared by every sample
dist = initial_distribution(OPT.npart, OPT.seed)
P₀ = dist.particles
w = vec(P₀.w)

# the reference problem's potential, on one wavelength at `KAPPA_REF`
potential = Potential(PeriodicBasisSpline(domain(KAPPA_REF), OPT.degree, OPT.ncells))
nh = length(VlasovMethods.coefficients(PoissonField(potential)))

# the sweep: one value of `χ` per sample, as an explicit list (K56)
κs = training_kappas(OPT.kappa_min, OPT.kappa_max, OPT.nsamples)
pspace = ParameterSpace(Parameter(:χ, [march_parameters(κ).χ for κ in κs]))

nt = step_count(OPT.final_time, OPT.tstep)
IP = IntegratorParameters(OPT.tstep, nt, nt + 1, nh, OPT.npart, length(κs))
SS = Snapshots(Float64, 1, OPT.npart, nh, nt + 1, pspace)

for p in eachindex(κs)
    params = march_parameters(κs[p])
    field = vp_field(potential, params)
    @printf("bump_on_tail_training: sample %d of %d, κ = %.4f, χ = %.4f\n",
        p, length(κs), κs[p], params.χ)
    integrate_vp(P₀, params, field, SS, IP, p)
    accelerations!(SS, field, w, IP.t, p)
end

# the snapshots, the integration parameters, the initial condition and the sweep
fpath = runs_path(OPT.runid * "_training.h5")
h5open(fpath, "w") do file
    h5save(file, SS; path = "snapshots")
    h5save(file, IP; path = "integrator")
    file["x"] = Matrix(P₀.x)
    file["v"] = Matrix(P₀.v)
    file["w"] = Matrix(P₀.w)
    file["kappa"] = κs
    HDF5.attributes(file)["commit"] = repository_commit()
    HDF5.attributes(file)["command"] = command_line("bump_on_tail_training.jl", ARGS)
    HDF5.attributes(file)["seed"] = OPT.seed
    HDF5.attributes(file)["ncells"] = OPT.ncells
    HDF5.attributes(file)["degree"] = OPT.degree
end
@printf("  wrote %s\n", fpath)

labels = [@sprintf("κ = %.3f", κ) for κ in κs]

# item 26: the field energy of each sample, which grows and saturates
fig = plot_series(IP.t, collect(eachcol(SS.W)); labels = labels, ylabel = "W")
fig.content[1].yscale = log10
figure = results_path(OPT.runid * "_training_field_energy.png")
save(figure, fig)
@printf("  wrote %s\n", figure)

# item 28: the growth rate of each sample against its wave number
γ = growth_rate(IP.t, SS.W)
fig = plot_series(κs, (γ,); labels = ["HFM"], xlabel = "κ", ylabel = "γ")
figure = results_path(OPT.runid * "_training_growth_rate.png")
save(figure, fig)
@printf("  wrote %s\n", figure)
