# The test sweep of the bump-on-tail example: the full-order march and the reduced march of the
# test samples, the diagnostics of the two against each other, and their figures.
#
#   julia --project=scripts scripts/bump_on_tail_testing.jl \
#       [--runid BoT_Np10000_k_010_050_ns10_T25] [--nsamples 10] [--seed 1234] [--order 3] \
#       [--ncells 16] [--runs-dir dir] [--results-dir dir]
#
# The snapshots go to `--runs-dir`, the figures to `--results-dir`. The order and the cell count
# rebuild the potential of the training run, and the script refuses to run when they do not give
# the same number of degrees of freedom.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using CairoMakie
using HDF5
using LinearAlgebra
using PoissonSolvers
using Printf
using ReducedBasisMethods: ReducedBasis, get_DEIM_indices
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods
using VlasovMethods: h5save, plot_particles, plot_series

const OPT = parse_arguments(
    ARGS; runid = DEFAULT_RUNID, nsamples = DEFAULT_NSAMPLES, seed = 1234,
    order = 3, ncells = 16)

const RUNID = OPT.runid
const PROJECTIONS = runs_path("$(RUNID)_projections.h5")

# the bases and the integrator parameters of the training run
const RB = ReducedBasis(PROJECTIONS, "particle")
const FB = ReducedBasis(PROJECTIONS, "field")
const IP = IntegratorParameters(PROJECTIONS, "integrator")

# the model of the test samples: the potential of the training run and the particles
const POTENTIAL = Potential(PeriodicBasisSpline(domain(), OPT.order, OPT.ncells))
const NH = potential_degrees(POTENTIAL)
NH == IP.nₕ || error(
    "the potential of the test run has $NH degrees of freedom; the training run has $(IP.nₕ)")
const DIST = initial_distribution(IP.nₚ; seed = OPT.seed)

# the test samples: the wave numbers, drawn uniformly from the range of the sweep and sorted
const KAPPA = κ_test(0.1, 0.5, OPT.nsamples; seed = OPT.seed)
const PSC = ParameterSpace(
    Parameter(:κ, 0.1, 0.5, KAPPA),
    Parameter(:ε, 0.03, 0.03, 1),
    Parameter(:α, 0.1, 0.1, 1),
    Parameter(:σ, 0.5, 0.5, 1),
    Parameter(:v₀, 4.5, 4.5, 1))

const IPTEST = IntegratorParameters(IP.dt, IP.nₜ, IP.nₛ, IP.nₕ, IP.nₚ, length(PSC))

# the bases of the two marches: the particle basis lifts the reduced coordinates, and the field
# basis is the one the DEIM interpolates
const Ψₚ = RB.Φ
const Ψₑ = FB.Φ
const INDICES = get_DEIM_indices(Ψₑ)
const X₀ = Ψₚ' * vec(DIST.particles.x)

const SSH = Snapshots(Float64, 1, IP.nₚ, NH, IP.nₛ, PSC)
const SSR = Snapshots(Float64, 1, IP.nₚ, NH, IP.nₛ, PSC)
const IC = ReducedIntegratorCache(IPTEST, size(Ψₚ, 2))

for p in eachindex(PSC)
    params = march_parameters(PSC(p).κ)

    # the full-order march, on the field of the sample
    integrate_vp(DIST.particles, params, vp_field(POTENTIAL, params), SSH, IPTEST, p)

    # the reduced march, on the DEIM field over the same scaled field
    rfield = DEIMElectricField(vp_field(POTENTIAL, params), Ψₚ, Ψₑ, INDICES, X₀)
    reduced_integrate_vp(DIST.particles, Ψₚ, params, rfield, SSR, IPTEST, IC, p)

    @printf("test sample %d of %d: κ = %g, χ = %g\n", p, length(PSC), PSC(p).κ, params.χ)
end

h5open(runs_path("$(RUNID)_test.h5"), "w") do file
    h5save(file, SSH)
    h5save(file, IPTEST; path = "integrator")
end

h5open(runs_path("$(RUNID)_reduced.h5"), "w") do file
    h5save(file, SSR)
    h5save(file, IPTEST; path = "integrator")
end

# the diagnostics of the two marches: the errors of the trajectories and of the potential (items 24
# and 25), the energies (items 26 and 27) and the growth rate (item 28)
const ERRORS = trajectory_errors(SSH.X[1, :, :, :], SSH.V[1, :, :, :],
    SSR.X[1, :, :, :], SSR.V[1, :, :, :])
const ΦERR = potential_error(SSH.Φ, SSR.Φ)
const HH = SSH.K .+ SSH.W
const HR = SSR.K .+ SSR.W
const DRIFT_H = [energy_drift(HH[:, p]) for p in axes(HH, 2)]
const DRIFT_R = [energy_drift(HR[:, p]) for p in axes(HR, 2)]

# The growth rate of the two marches: the fit runs over the local maxima of the field energy, and
# the reduced history has fewer of them than the full-order one, so the range holds the first two
# maxima, which both histories hold at the paper's parameters.
const MAXIMA = 1:2
const WINDOW = 2
const ΓH = growth_rate(IPTEST.t, SSH.W; window = WINDOW, maxima = MAXIMA)
const ΓR = growth_rate(IPTEST.t, SSR.W; window = WINDOW, maxima = MAXIMA)
const BEST, WORST = best_and_worst(ΓH, ΓR)

@printf("the growth rates agree best at the sample %d and worst at the sample %d\n", BEST,
    WORST)

# the saved times of the trends of item 24: the last one, the one nearest three quarters of the
# way and the one nearest halfway
const SAVED_TIMES = saved_indices(IPTEST.nₛ)
const TREND_LABELS = ["t = T", "t = 0.75T", "t = 0.5T"]

# the figures: the errors of the first test sample against the time (item 24), the errors of every
# test sample against the wave number at the three saved times (item 24), the potential error of
# every sample (item 25), the field energy (item 26), the total energy and its drift (item 27),
# the growth rate and the field energy of the best and the worst sample (item 28) and the
# phase-space frames (item 29)
fig = plot_series(
    IPTEST.t, (ERRORS.x2[:, 1], ERRORS.xmax[:, 1], ERRORS.v2[:, 1],
        ERRORS.vmax[:, 1]);
    labels = ["x2", "xmax", "v2", "vmax"], xlabel = "t", ylabel = "error")
save(results_path("$(RUNID)_testing_errors.png"), fig)

for (values, name) in ((ERRORS.x2, "x2"), (ERRORS.xmax, "xmax"), (ERRORS.v2, "v2"),
    (ERRORS.vmax, "vmax"))
    errfig = plot_series(KAPPA, (values[ts, :] for ts in SAVED_TIMES);
        labels = TREND_LABELS, xlabel = "κ", ylabel = "error")
    save(results_path("$(RUNID)_testing_error_$(name).png"), errfig)
end

fig = plot_series(IPTEST.t, (ΦERR[:, p] for p in axes(ΦERR, 2));
    labels = ["κ = $κ" for κ in KAPPA], xlabel = "t", ylabel = "‖Φᵣₒₘ - Φₕfₘ‖")
save(results_path("$(RUNID)_testing_potential_error.png"), fig)

fig = plot_series(IPTEST.t, (SSH.W[:, 1], SSR.W[:, 1]);
    labels = ["HFM", "ROM"], xlabel = "t", ylabel = "W")
# the field energy is positive at every saved time, so the axis takes a logarithmic scale
fig.content[1, 1].yscale = log10
save(results_path("$(RUNID)_testing_field_energy.png"), fig)

fig = plot_series(IPTEST.t, (HH[:, 1], HR[:, 1]);
    labels = ["HFM", "ROM"], xlabel = "t", ylabel = "H")
save(results_path("$(RUNID)_testing_total_energy.png"), fig)

fig = plot_series(IPTEST.t, (DRIFT_H[1], DRIFT_R[1]);
    labels = ["HFM", "ROM"], xlabel = "t", ylabel = "ΔH/H")
save(results_path("$(RUNID)_testing_energy_drift.png"), fig)

fig = plot_series(KAPPA, (ΓH, ΓR);
    labels = ["HFM", "ROM"], xlabel = "κ", ylabel = "γ")
save(results_path("$(RUNID)_testing_growth_rate.png"), fig)

for (sample, name) in ((BEST, "best"), (WORST, "worst"))
    wfig = plot_series(IPTEST.t, (SSH.W[:, sample], SSR.W[:, sample]);
        labels = ["HFM", "ROM"], xlabel = "t", ylabel = "W")
    wfig.content[1, 1].yscale = log10
    save(results_path("$(RUNID)_testing_field_energy_$(name).png"), wfig)
end

vdomain = (minimum(SSR.V), maximum(SSR.V))
for (ts, name) in ((1, "t0"), (IPTEST.nₛ, "tT"))
    frame = plot_particles(SSR, ts, 1; xdomain = domain(), vdomain = vdomain)
    save(results_path("$(RUNID)_testing_particles_$(name).png"), frame)
end
