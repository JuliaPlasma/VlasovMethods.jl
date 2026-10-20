# The test sweep of the bump-on-tail example: the full-order march and the reduced march of the
# test samples, the diagnostics of the two against each other, and their figures.
#
#   julia --project=scripts scripts/bump_on_tail_testing.jl \
#       [--runid BoT_Np1000_k_010_050_ns10_T25] [--nsamples 10] [--seed 1234] [--order 3] \
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
# and 25), the relative error of the field energy (item 26), the drift of the total energy (item
# 27) and the growth rate (item 28)
const ERRORS = trajectory_errors(SSH.X[1, :, :, :], SSH.V[1, :, :, :],
    SSR.X[1, :, :, :], SSR.V[1, :, :, :])
const ΦERR = potential_error(SSH.Φ, SSR.Φ)
const WERR = field_energy_error(SSH.W, SSR.W)
const HH = SSH.K .+ SSH.W
const HR = SSR.K .+ SSR.W
const DRIFT_H = [energy_drift(HH[:, p]) for p in axes(HH, 2)]
const DRIFT_R = [energy_drift(HR[:, p]) for p in axes(HR, 2)]
const ΓH = growth_rate(IPTEST.t, SSH.W)
const ΓR = growth_rate(IPTEST.t, SSR.W)
const BEST, WORST = best_and_worst(ΓH, ΓR)

@printf("the growth rates agree best at the sample %d and worst at the sample %d\n", BEST,
    WORST)

# the figures of the first test sample: the errors (item 24), the potential error (item 25), the
# field energy (item 26), the drift of the total energy (item 27) and the phase-space frames (item
# 29)
fig = plot_series(
    IPTEST.t, (ERRORS.x2[:, 1], ERRORS.xmax[:, 1], ERRORS.v2[:, 1],
        ERRORS.vmax[:, 1]);
    labels = ["x2", "xmax", "v2", "vmax"], xlabel = "t", ylabel = "error")
save(results_path("$(RUNID)_testing_errors.png"), fig)

fig = plot_series(IPTEST.t, (ΦERR[:, p] for p in axes(ΦERR, 2));
    labels = ["κ = $κ" for κ in KAPPA], xlabel = "t", ylabel = "‖Φᵣₒₘ - Φₕfₘ‖")
save(results_path("$(RUNID)_testing_potential_error.png"), fig)

fig = plot_series(IPTEST.t, (SSH.W[:, 1], SSR.W[:, 1]);
    labels = ["HFM", "ROM"], xlabel = "t", ylabel = "W")
# the field energy is positive at every saved time, so the axis takes a logarithmic scale
fig.content[1, 1].yscale = log10
save(results_path("$(RUNID)_testing_field_energy.png"), fig)

fig = plot_series(IPTEST.t, (DRIFT_H[1], DRIFT_R[1]);
    labels = ["HFM", "ROM"], xlabel = "t", ylabel = "ΔH/H")
save(results_path("$(RUNID)_testing_energy_drift.png"), fig)

fig = plot_series(KAPPA, (ΓH, ΓR);
    labels = ["HFM", "ROM"], xlabel = "κ", ylabel = "γ")
save(results_path("$(RUNID)_testing_growth_rate.png"), fig)

vdomain = (minimum(SSR.V), maximum(SSR.V))
for (ts, name) in ((1, "t0"), (IPTEST.nₛ, "tT"))
    frame = plot_particles(SSR, ts, 1; xdomain = domain(), vdomain = vdomain)
    save(results_path("$(RUNID)_testing_particles_$(name).png"), frame)
end
