# The full-order Vlasov–Poisson driver: a particle cloud on a periodic spline basis, advanced by
# the splitting method, with the phase-space scatter of the cloud drawn at four times of the run.
#
#   julia --project=scripts scripts/vlasov_poisson.jl [--npart 10000] [--ncells 16] [--order 3] \
#       [--tstep 0.1] [--final-time 20.0] [--runs-dir dir] [--results-dir dir]
#
# The HDF5 trajectory goes to `--runs-dir`, the figure to `--results-dir`. The script writes only
# under those two, so a run adds no untracked path and `git status --short` is unchanged by it. The
# two directories are ignored, so `git status --short --ignored` gains a line for each on the run
# that creates it.
#
# `BumpOnTail` is not built here: the bump on the tail belongs to the reproduction set for that
# problem, and this driver is the deposition path and the splitting method.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using HDF5
using PoissonSolvers
using VlasovMethods

const PARSED = let (positional, options) = parse_arguments(ARGS,
        ("--npart", "--ncells", "--order", "--tstep", "--final-time"))
    no_positional(positional, "vlasov_poisson.jl")
    options
end

const NPART = option_int(PARSED, "--npart", 10000)
const NCELLS = option_int(PARSED, "--ncells", 16)
const ORDER = option_int(PARSED, "--order", 3)
const TSTEP = option_float(PARSED, "--tstep", 0.1)
const TSPAN = (0.0, option_float(PARSED, "--final-time", 20.0))
const DOMAIN = (0.0, 1.0)

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, NPART), NormalDistribution())

# create electrostatic potential
potential = Potential(PeriodicBasisSpline(DOMAIN, ORDER, NCELLS))

# create Vlasov-Poisson model
model = VlasovPoisson(dist, potential)

# create integrator
integrator = SplittingMethod(model, TSPAN, TSTEP)

# integrate, writing the trajectory to the runs directory
h5file = runs_path("vlasov_poisson.hdf5")
run!(integrator, h5file)

@printf("vlasov_poisson: %d particles, %d cells, order %d, dt = %g, t ∈ %s\n",
    NPART, NCELLS, ORDER, TSTEP, TSPAN)

# `z` is `(nd, npart, ntime)`; the two rows are `x` and `v`. `run!(::SplittingMethod, h5file)`
# writes the trajectory alone, so the time of a slice is the step count times the step size.
z = h5read(h5file, "z")
@printf("  wrote %s (%d time slices)\n", h5file, size(z, 3))
slice_time(slice) = TSPAN[1] + (slice - 1) * TSTEP

# the figure shows the same cloud at four evenly spaced times, so that the filamentation the
# scheme is meant to resolve is visible rather than asserted. The two signs of `v` are drawn
# apart, split once by each particle's sign at the first slice: the filamentation of either half
# is easier to read than the sum, and a particle that crosses `v = 0` keeps the half it began in
slices = unique(round.(Int, range(1, size(z, 3); length = 4)))
positive = z[2, :, 1] .>= 0

vmax = ceil(maximum(abs.(z[2, :, begin])))

fig = Figure(size = (1500, 400))
for (panel, slice) in enumerate(slices)
    ax = Axis(fig[1, panel], xlabel = "x", ylabel = panel == 1 ? "v" : "",
        title = @sprintf("t = %.1f", slice_time(slice)), limits = (
            (0.0, 1.0), (-vmax, vmax)),
        xticks = [0.0, 0.5, 1.0])
    scatter!(ax, mod.(z[1, positive, slice], 1), z[2, positive, slice],
        color = :steelblue, markersize = 2)
    scatter!(ax, mod.(z[1, .!positive, slice], 1), z[2, .!positive, slice],
        color = :firebrick, markersize = 2)
end

figure = results_path("vlasov_poisson.png")
save(figure, fig)
@printf("  wrote %s\n", figure)
