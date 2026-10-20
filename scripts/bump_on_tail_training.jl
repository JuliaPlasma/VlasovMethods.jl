# The training sweep of the bump-on-tail example: the full-order march of the samples of the
# wave number, the snapshots it saves and the figures of the field energy, of the total energy and
# of the phase-space frames.
#
#   julia --project=scripts scripts/bump_on_tail_training.jl [--npart 1000] [--ncells 16] \
#       [--order 3] [--tstep 0.1] [--final-time 25.0] [--runs-dir dir] [--results-dir dir]
#
# The snapshots go to `--runs-dir`, the figures to `--results-dir`. The paper's parameters are
# `--npart 10000`; the defaults are smaller, so that the run stays inside a few minutes.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using CairoMakie
using HDF5
using PoissonSolvers
using Printf
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods
using VlasovMethods: efield!, h5save, plot_particles, plot_series, update!

const OPT = parse_arguments(
    ARGS; npart = DEFAULT_NPART, ncells = 16, order = 3, tstep = 0.1,
    final_time = DEFAULT_FINAL_TIME)

const NPART = OPT.npart
const NCELLS = OPT.ncells
const ORDER = OPT.order
const TSTEP = OPT.tstep
const FINAL_TIME = OPT.final_time

# the sweep of the wave number and the model's fixed parameters
const NSAMPLES = 10
const KAPPA_MIN = 0.1
const KAPPA_MAX = 0.5

# The step count is `round` and not the floating-point `div` of the old script, which gave 249 for
# `T = 25` and `dt = 0.1`. The integrator parameters need `nₛ - 1` to divide `nₜ`, and `nₛ = nₜ + 1`
# divides it.
const NT = step_count(FINAL_TIME, TSTEP)
const NS = NT + 1

const RUNID = bump_on_tail_runid(NPART, NSAMPLES, FINAL_TIME)

# the parameter space: the wave numbers of the sweep and the parameters the model holds fixed
const KAPPA = κ_sweep(KAPPA_MIN, KAPPA_MAX, NSAMPLES)
const PSC = ParameterSpace(
    Parameter(:κ, KAPPA_MIN, KAPPA_MAX, KAPPA),
    Parameter(:ε, 0.03, 0.03, 1),
    Parameter(:α, 0.1, 0.1, 1),
    Parameter(:σ, 0.5, 0.5, 1),
    Parameter(:v₀, 4.5, 4.5, 1))

# the model: the potential on the domain of one period of the perturbation, and the particles
const DOMAIN = domain()
const POTENTIAL = Potential(PeriodicBasisSpline(DOMAIN, ORDER, NCELLS))
const NH = potential_degrees(POTENTIAL)
const DIST = initial_distribution(NPART)

const IP = IntegratorParameters(TSTEP, NT, NS, NH, NPART, length(PSC))
const SS = Snapshots(Float64, 1, NPART, NH, NS, PSC)

# the sweep: every sample is marched at its own `χ`, and the field of the sample carries it
for p in eachindex(PSC)
    params = march_parameters(PSC(p).κ)
    efield = vp_field(POTENTIAL, params)
    integrate_vp(DIST.particles, params, efield, SS, IP, p)
    @printf("training sample %d of %d: κ = %g, χ = %g\n", p, length(PSC), PSC(p).κ,
        params.χ)
end

# The accelerations of the saved states, which the marches do not save: the field of the sample is
# updated at the positions of the state, and the scaled field divides the field by χ², so that the
# acceleration is the one the march kicked with.
for p in eachindex(PSC)
    params = march_parameters(PSC(p).κ)
    efield = vp_field(POTENTIAL, params)
    for it in 1:NS
        x = view(SS.X, 1, :, it, p)
        update!(efield, x, vec(DIST.particles.w), IP.t[it])
        efield!(efield, view(SS.A, 1, :, it, p), x)
    end
end

const H5FILE = runs_path("$(RUNID).h5")
h5open(H5FILE, "w") do file
    h5save(file, SS)
    h5save(file, IP; path = "integrator")
    HDF5.attributes(file)["runid"] = RUNID
    HDF5.attributes(file)["commit"] = repository_commit()
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = 1234
end

@printf("wrote %s (%d saved times, %d samples)\n", H5FILE, NS, length(PSC))

# the figures: the field energy of every sample (item 26), the relative drift of the total energy
# (item 27) and the phase-space frames at the first and the last saved time (item 29)
const W = SS.W
const H = SS.K .+ SS.W

fig = plot_series(IP.t, (W[:, p] for p in axes(W, 2));
    labels = ["κ = $κ" for κ in KAPPA], xlabel = "t", ylabel = "W")
# the field energy is positive at every saved time, so the axis takes a logarithmic scale. The
# block of the figure is the axis that `plot_series` drew.
fig.content[1, 1].yscale = log10
save(results_path("$(RUNID)_training_field_energy.png"), fig)

fig = plot_series(IP.t, (energy_drift(H[:, p]) for p in axes(H, 2));
    labels = ["κ = $κ" for κ in KAPPA], xlabel = "t", ylabel = "ΔH/H")
save(results_path("$(RUNID)_training_energy_drift.png"), fig)

vdomain = (minimum(SS.V), maximum(SS.V))
for (ts, name) in ((1, "t0"), (NS, "tT"))
    frame = plot_particles(SS, ts, 1; xdomain = DOMAIN, vdomain = vdomain)
    save(results_path("$(RUNID)_training_particles_$(name).png"), frame)
end

@printf("wrote the figures of %s\n", results_path("$(RUNID)_training_field_energy.png"))
