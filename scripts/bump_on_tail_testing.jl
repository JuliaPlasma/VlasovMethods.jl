# The test sweep of the self-consistent bump-on-tail set: ten κ drawn uniformly from [0.1, 0.5]
# under seed 1234 and sorted, each marched in full and reduced order, and the paper's diagnostics
# of the two against each other. The reduced model uses the particle and field bases that
# `bump_on_tail_projections.jl` wrote, with the DEIM reduction of the field.
#
#   julia --project=scripts scripts/bump_on_tail_testing.jl [--particles 10000] [--ncells 16] \
#       [--degree 3] [--tstep 0.1] [--final-time 25.0] [--runs-dir dir] [--results-dir dir]
#
# Reads `<runid>_projections.h5` from `--runs-dir` and writes the figures of items 24–29 to
# `--results-dir`: the errors per time and against the parameter, the potential error, the field
# energies and the total energy and its drift of the best and worst sample, the growth rates, and
# one phase-space frame at t = 0 and one at T.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using PoissonSolvers: PeriodicBasisSpline, Potential
using ReducedBasisMethods: ReducedBasis, get_DEIM_indices
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods

const OPT = parse_arguments(ARGS; particles = NPARTICLES_BOT, ncells = NCELLS_BOT,
    degree = DEGREE_BOT, tstep = DT_BOT, final_time = T_BOT)

const T = Float64
const NPARTICLES = OPT.particles
const DT = OPT.tstep
const NT = round(Int, OPT.final_time / DT)
const NS = NT + 1
const CHI = test_chi()

basis = PeriodicBasisSpline((zero(T), domain_length()), OPT.degree, OPT.ncells)
potential = Potential(basis)
nh = nbasis(basis)

# the bases of the training sweep, and the DEIM rows of the field basis
ppath = runs_path("$(RUNID_BOT)_projections.h5")
bp = h5open(ppath, "r") do file
    ReducedBasis(file)
end
be = h5open(ppath, "r") do file
    ReducedBasis(file, "field")
end
Ψₚ, Ψₑ = bp.Φ, be.Φ
deim = get_DEIM_indices(Ψₑ)

dist = initial_distribution(; n = NPARTICLES)
ps = ParameterSpace(
    Parameter(:χ, first(CHI), last(CHI), CHI),
    Parameter(:ε, ε_BOT, ε_BOT, 1),
    Parameter(:a, α_BOT, α_BOT, 1),
    Parameter(:v₀, v₀_BOT, v₀_BOT, 1),
    Parameter(:σ, σ_BOT, σ_BOT, 1))

SS = Snapshots(T, 1, NPARTICLES, nh, NS, ps)
IP = IntegratorParameters(DT, NT, NS, nh, NPARTICLES, length(ps))

@printf("bump_on_tail_testing: %d particles, %d cells, degree %d, dt = %g, t ∈ [0, %g]\n",
    NPARTICLES, OPT.ncells, OPT.degree, DT, OPT.final_time)
@printf("  kₚ = %d, kₑ = %d, %d DEIM rows\n", bp.k, be.k, length(deim))

# the full-order march of every test sample, each on the sample's own scaled field, built from the
# parameter tuple the march reads
for p in eachindex(ps)
    params = march_parameters(ps(p).χ)
    integrate_vp(dist.particles, params, vp_field(potential, params), SS, IP, p)
end

# the reduced march of every test sample, on the DEIM field of that same scaled field
X0 = Ψₚ' * vec(dist.particles.x)
RSS = Snapshots(T, 1, NPARTICLES, nh, NS, ps)
RIC = ReducedIntegratorCache(IP, bp.k)
for p in eachindex(ps)
    params = march_parameters(ps(p).χ)
    efield = DEIMElectricField(vp_field(potential, params), Ψₚ, Ψₑ, deim, X0)
    reduced_integrate_vp(dist.particles, Ψₚ, params, efield, RSS, IP, RIC, p)
end

# the diagnostics of item 24–29, each on the arrays the marches wrote
t = IP.t
err = trajectory_errors(SS.X[1, :, :, :], SS.V[1, :, :, :], RSS.X[1, :, :, :],
    RSS.V[1, :, :, :])
pot = potential_error(SS.Φ, RSS.Φ)
Wh, Wr = SS.W, RSS.W
Hh = SS.K .+ Wh
Hr = RSS.K .+ Wr
γh = growth_rate(t, Wh; maxima = 1:2)
γr = growth_rate(t, Wr; maxima = 1:2)
best, worst = best_and_worst(γh, γr)

@printf("  γ_χ2 = %.4f (HFM) vs %.4f (ROM), best χ = %.4f, worst χ = %.4f\n",
    γh[2], γr[2], ps(best).χ, ps(worst).χ)

# item 24, the errors per time, at the largest over the samples
figure = results_path("bump_on_tail_errors_time.png")
save(figure,
    plot_series(t,
        (vec(maximum(err.x2; dims = 2)), vec(maximum(err.xmax; dims = 2)),
            vec(maximum(err.v2; dims = 2)), vec(maximum(err.vmax; dims = 2)));
        labels = ["RMS x", "L∞ x", "RMS v", "L∞ v"], xlabel = "t", ylabel = "error"))
@printf("  wrote %s\n", figure)

# item 24, both norms of both components against the parameter at T, 0.75T and 0.5T
times = ("T", "0.75T", "0.5T")
norms = ((:x2, "RMS x"), (:xmax, "L∞ x"), (:v2, "RMS v"), (:vmax, "L∞ v"))
ts = saved_indices(NS)
figure = results_path("bump_on_tail_errors_parameter.png")
save(figure,
    plot_series(CHI,
        [getproperty(err, name)[ts[i], :] for i in eachindex(times) for (name, _) in norms];
        labels = ["$(times[i]), $(label)" for i in eachindex(times) for (_, label) in norms],
        xlabel = "χ", ylabel = "error"))
@printf("  wrote %s\n", figure)

# item 25
figure = results_path("bump_on_tail_potential_error.png")
save(figure,
    plot_series(t, (vec(maximum(pot; dims = 2)),);
        labels = ["max over χ"], xlabel = "t", ylabel = "‖Φ_ROM − Φ_HFM‖"))
@printf("  wrote %s\n", figure)

# items 26 and 27, on the best and the worst sample of the growth rate
figure = results_path("bump_on_tail_energy.png")
save(figure,
    plot_series(t,
        (Wh[:, best], Wr[:, best], Wh[:, worst], Wr[:, worst]);
        labels = ["HFM, best", "ROM, best", "HFM, worst", "ROM, worst"],
        xlabel = "t", ylabel = "W"))
@printf("  wrote %s\n", figure)

figure = results_path("bump_on_tail_energy_total.png")
save(figure,
    plot_series(t,
        (Hh[:, best], Hr[:, best], Hh[:, worst], Hr[:, worst]);
        labels = ["HFM, best", "ROM, best", "HFM, worst", "ROM, worst"],
        xlabel = "t", ylabel = "H = K + W"))
@printf("  wrote %s\n", figure)

figure = results_path("bump_on_tail_energy_drift.png")
save(figure,
    plot_series(t,
        (energy_drift(Hh[:, best]), energy_drift(Hr[:, best]),
            energy_drift(Hh[:, worst]), energy_drift(Hr[:, worst]));
        labels = ["HFM, best", "ROM, best", "HFM, worst", "ROM, worst"],
        xlabel = "t", ylabel = "ΔH/H"))
@printf("  wrote %s\n", figure)

# item 28
figure = results_path("bump_on_tail_growth_rate.png")
save(figure, plot_series(CHI, (γh, γr);
    labels = ["HFM", "ROM"], xlabel = "χ", ylabel = "γ"))
@printf("  wrote %s\n", figure)

# item 29, the phase space of the best sample at the initial and the final saved state
xdomain = (zero(T), domain_length())
vdomain = (-10.0, 10.0)
for (name, ts) in (("t0", 1), ("T", NS))
    path = results_path("bump_on_tail_phase_$name.png")
    save(path, plot_particles(SS, ts, best; xdomain, vdomain))
    @printf("  wrote %s\n", path)
end
