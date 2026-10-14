# The Lenard–Bernstein driver: a normal-distribution particle cloud on a cubic spline basis,
# advanced by the model's own `GeometricIntegrator`, with the right-hand side at the initial
# condition and the distribution at both ends of the run drawn as figures.
#
#   julia --project=scripts scripts/lenard_bernstein.jl [--model lb|clb|rclb] [--npart 1000] \
#       [--nknot 41] [--order 4] [--tstep 0.01] [--final-time 1.0] \
#       [--runs-dir dir] [--results-dir dir]
#
# The HDF5 trajectory goes to `--runs-dir`, the figure to `--results-dir`.
#
# `--model` picks the collision operator: `lb` is `LenardBernstein`, `clb` the conservative
# variant that solves its two coefficients from the particle momentum and energy, and `rclb` the
# rescaled conservative variant. `clb` is the default, which is what this driver built before
# the flag existed.
#
# **The defaults are dt = 0.01 over t ∈ [0, 1]**, and not the dt = 0.1 over t ∈ [0, 10] this
# driver carried: that pair fails for every one of the three models, because the particles drift
# out of the velocity domain and the projection then throws on the ones that left.
#
# The animation this driver used to record is a figure of three panels instead. CairoMakie draws
# no animated GIF without an encoder on the `PATH`, and the two ends of the run are what the
# animation was read for.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using HDF5
using VlasovMethods

const MODELS = ("lb", "clb", "rclb")

const PARSED = let (positional, options) = parse_arguments(ARGS,
        ("--model", "--npart", "--nknot", "--order", "--tstep", "--final-time"))
    no_positional(positional, "lenard_bernstein.jl")
    options
end

const MODEL = option_choice(PARSED, "--model", "clb", MODELS)
const NPART = option_int(PARSED, "--npart", 1000)
const NKNOT = option_int(PARSED, "--nknot", 41)
const ORDER = option_int(PARSED, "--order", 4)
const TSTEP = option_float(PARSED, "--tstep", 1e-2)
const TSPAN = (0.0, option_float(PARSED, "--final-time", 1.0))
const DOMAIN = (-10.0, 10.0)

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, NPART), NormalDistribution())

# create spline distribution function and entropy
sdist = SplineDistribution(1, 1, NKNOT, ORDER, DOMAIN, 0, :Dirichlet)
entropy = CollisionEntropy(sdist)

# create the collision model the flag names
model = if MODEL == "lb"
    LenardBernstein(dist, entropy)
elseif MODEL == "clb"
    ConservativeLenardBernstein(dist, entropy)
else
    RescaledConservativeLenardBernstein(dist, entropy)
end

# The right-hand side at the initial condition, saved before the run writes the final state back
# into `dist.particles.v`. This is the diagnostic the deleted `lenard_bernstein_conservative.jl`
# carried: one vector-field evaluation against `v`, which shows where the operator pushes hardest
# and whether the projection is positive. The argument order is the ODE one, which is what the
# integrator calls.
const RHS! = MODEL == "lb" ? VlasovMethods.LB_rhs_GI! :
             MODEL == "clb" ? VlasovMethods.CLB_rhs_GI! :
             VlasovMethods.RCLB_rhs_GI!

params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
v₀ = copy(dist.particles.v[1, :])
v̇₀ = zeros(NPART)
RHS!(v̇₀, 0.0, v₀, params)

# create integrator and run, writing the trajectory to the runs directory
integrator = VlasovMethods.GeometricIntegrator(model, TSPAN, TSTEP)
h5file = runs_path("lenard_bernstein_" * MODEL * ".hdf5")
run!(integrator, h5file)

@printf("lenard_bernstein: model %s, %d particles, %d knots, order %d, dt = %g, t ∈ %s\n",
    MODEL, NPART, NKNOT, ORDER, TSTEP, TSPAN)
@printf("  wrote %s (%d time slices)\n", h5file, size(h5read(h5file, "z"), 2))

# `run!(::GeometricIntegrator, h5file)` writes the velocity vector the ODE steps, one slice per
# time point: `z` is `(npart, ntime + 1)`, not the `(nd, npart, ntime + 1)` of the splitting
# method's trajectory.
z = h5read(h5file, "z")
xgrid = collect(range(DOMAIN[1], DOMAIN[2]; length = 201))

# The particle density of one time slice, by the rectangle rule of the plot's own bins. Written
# out rather than taken from `hist`, so that the two panels are the same measurement and no
# plotting-package default decides what is drawn.
function particle_density(velocities, edges)
    width = edges[2] - edges[1]
    counts = zeros(length(edges) - 1)
    for velocity in velocities
        i = searchsortedlast(edges, velocity)
        1 ≤ i ≤ length(counts) && (counts[i] += 1)
    end
    return counts ./ (length(velocities) * width)
end

edges = collect(range(DOMAIN[1], DOMAIN[2]; length = 81))
centres = (edges[1:(end - 1)] .+ edges[2:end]) ./ 2

fig = Figure(size = (1600, 450))

ax1 = Axis(fig[1, 1], xlabel = "v", ylabel = "v̇", title = "right-hand side at t = 0")
scatter!(ax1, v₀, v̇₀, markersize = 2)

for (panel, slice) in enumerate((1, size(z, 2)))
    f = projection(z[:, slice], dist, sdist)
    df = derivative(f)
    ax = Axis(fig[1, panel + 1], xlabel = "v", ylabel = "f",
        title = @sprintf("t = %.2f", TSPAN[1] + (slice - 1) * TSTEP))
    lines!(ax, centres, particle_density(z[:, slice], edges),
        color = :steelblue, label = "particles")
    lines!(ax, xgrid, f.(xgrid), color = :firebrick, label = "spline projection")
    lines!(ax, xgrid, df.(xgrid), color = :black, linestyle = :dash, label = "derivative")
    panel == 1 && axislegend(ax; position = :lt)
end

figure = results_path("lenard_bernstein_" * MODEL * ".png")
save(figure, fig)
@printf("  wrote %s\n", figure)
